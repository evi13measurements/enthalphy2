!> @file
!> Contains the soeo_itemroutines module
!
!> brief Routines operating on the soeoItem structure
!> author C. Nygaard
!> date 2010.07.19
module soeo_itemroutines

use soeo_typedef
use soeo_transform
use soeo_dens
use soeo_redspace
use matrix_util
use Fock_evaluator
use Matrix_Operations
use decompMod

implicit none

!Contains
!  soeoItem_init
!  soeoItem_setold
!  soeoIten_stepforth
!  soeoItem stepback
!  soeoItem_update_space
!  soeoItem_shutdown

contains

!> brief Initialize soeo-structure
!> author C. Nygaard
!> date 2010.07.07
!> param Dao Density matrix in the AO basis
!> param Fao Fock matrix in the AO basis
!> param S Overlap matrix in the AO basis
!> param decomp Structure containing information about the calculation
!> param soeo Structure containing information for SOEO calculation
!=======================================================================
subroutine soeoItem_init (Dao, Fao, S, decomp, soeo)

implicit none

type(soeoItem), intent(inout) :: soeo
type(matrix), intent(in)      :: Dao, Fao, S
type(decompItem), intent(in)  :: decomp

real(realk), allocatable      :: orbE(:), tmp(:)
integer                       :: ndim, i, j

!Space section
!-----------------------------------------------------------------------
soeo%Nbast = S%nrow
soeo%Nocc  = 0
soeo%Nact  = soeo%Nbast
!-----------------------------------------------------------------------

!Matrices section
!-----------------------------------------------------------------------
if (decomp%cfg_unres) then
  allocate (orbE(2*soeo%Nbast), tmp(2))
else
  allocate (orbE(soeo%Nbast), tmp(1))
endif

call mat_init (soeo%Dao     , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%Fao     , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%S       , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%C       , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%Dmo     , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%Fmo     , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%oldtheta, soeo%Nbast, 1         )
if (decomp%cfg_unres) then
  ndim = 2.0d0*soeo%Nact
else
  ndim = soeo%Nact
endif
allocate (soeo%nfirst(ndim),&
        & soeo%nsecond(ndim,ndim))
do i=1,ndim
  call mat_init(soeo%nfirst(i), soeo%Nbast, soeo%Nbast)
  do j=1,ndim
    call mat_init (soeo%nsecond(i,j), soeo%Nbast, soeo%Nbast)
  enddo
enddo

call mat_init (soeo%old_Dao     , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%old_Fao     , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%old_C       , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%old_Dmo     , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%old_Fmo     , soeo%Nbast, soeo%Nbast)
call mat_init (soeo%old_oldtheta, soeo%Nbast, 1    )
allocate (soeo%old_nfirst(ndim),&
        & soeo%old_nsecond(ndim,ndim))
do i=1,ndim
  call mat_init(soeo%old_nfirst(i), soeo%Nbast, soeo%Nbast)
  do j=1,ndim
    call mat_init (soeo%old_nsecond(i,j), soeo%Nbast, soeo%Nbast)
  enddo
enddo

soeo%Dao = Dao
soeo%Fao = Fao
soeo%S   = S

call mat_diag_f (Fao, S, orbE, soeo%C)

call mat_zero (soeo%Fmo)
do i=1,soeo%Nbast
  tmp(1) = orbE(i)
  if (decomp%cfg_unres) then
    tmp(2) = orbE(soeo%Nbast+i)
  endif
  call mat_create_ab_elms (i, i, tmp, soeo%Fmo)
enddo

call soeo_Dmo_init (S, soeo%C, Dao, soeo%Dmo, decomp)
call soeo_oldtheta_init (soeo%Dmo, soeo%oldtheta, decomp)

do i=1,ndim
  if (i <= soeo%Nact) then
    call soeo_get_nfirst(soeo%oldtheta, i, 'a', soeo%nfirst(i), decomp)
  else
    call soeo_get_nfirst(soeo%oldtheta, i-soeo%Nact, 'b', soeo%nfirst(i), decomp)
  endif
  do j=1,ndim
    if (i<=soeo%Nact) then
      if (j<=soeo%Nact) then
        call soeo_get_nsecond(soeo%oldtheta, i, j, 'a', 'a', soeo%nsecond(i,j), decomp)
      else
        call soeo_get_nsecond(soeo%oldtheta, i, j-soeo%Nact, 'a', 'b', soeo%nsecond(i,j), decomp)
      endif
    else
      if (j<=soeo%Nact) then
        call soeo_get_nsecond(soeo%oldtheta, i-soeo%Nact, j, 'b', 'a', soeo%nsecond(i,j), decomp)
      else
        call soeo_get_nsecond(soeo%oldtheta, i-soeo%Nact, j-soeo%Nact, 'b', 'b', soeo%nsecond(i,j), decomp)
      endif
    endif
  enddo
enddo
!-----------------------------------------------------------------------

deallocate (orbE, tmp)

!Energy section
!-----------------------------------------------------------------------
soeo%Etotal = 0.0d0
soeo%dEpred = 0.0d0
soeo%dE     = 0.0d0
!-----------------------------------------------------------------------

!Thresholds and stuff section
!-----------------------------------------------------------------------
soeo%macromaxiter = 200
soeo%macrothresh = 1.0d-5
soeo%micromaxiter = 200
soeo%microthresh = 1.0d-8
soeo%trust = 0.5d0
!-----------------------------------------------------------------------

end subroutine soeoItem_init
!=======================================================================

!> brief Saves the matrices from last iteration
!> author C. Nygaard
!> date 2010.07.14
!> param soeo The structure where the matrices are saved
!=======================================================================
subroutine soeoItem_setold (soeo)

implicit none

type(soeoItem), intent(inout) :: soeo

integer                       :: ndim, i, j

soeo%old_Dao      = soeo%Dao
soeo%old_Fao      = soeo%Fao
soeo%old_C        = soeo%C
soeo%old_Dmo      = soeo%Dmo
soeo%old_Fmo      = soeo%Fmo
soeo%old_oldtheta = soeo%oldtheta

ndim = size(soeo%nfirst)
do i=1,ndim
  soeo%old_nfirst(i) = soeo%nfirst(i)
  do j=1,ndim
    soeo%old_nsecond(i,j) = soeo%nsecond(i,j)
  enddo
enddo

end subroutine soeoItem_setold
!=======================================================================

!> brief Updates the matrices used in a SOEO calculation
!> author C. Nygaard
!> date 2010.07.08
!> param K The orbital rotation matrix
!> param deltatheta The change in occupations
!> param decomp Contains information about the calculations
!> param soeo Contains the matrices and information used in the SOEO calculation
!=======================================================================
subroutine soeoItem_stepforth (K, deltatheta, decomp, soeo)

implicit none

type(matrix), intent(in)      :: K, deltatheta
type(soeoItem), intent(inout) :: soeo
type(decompItem), intent(in)  :: decomp

integer                       :: ndim, i, j
real(realk), allocatable      :: orbE(:), tmp(:)

type(matrix)                  :: tmpmat

call soeoItem_setold (soeo)

if (decomp%cfg_unres) then
  allocate (orbE(2*soeo%Nbast), tmp(2))
else
  allocate (orbE(soeo%Nbast), tmp(1))
endif

if (K%nrow /= soeo%Nbast .or. K%ncol /= soeo%Nbast) then
  call lsquit ('Wrong dimensions of K in soeoItem_update')
endif
if (deltatheta%nrow /= soeo%Nact .or. deltatheta%ncol /= 1) then
  call lsquit ('Wrong dimensions of deltatheta in soeoItem_update')
endif

soeo%dE = soeo%Etotal

call soeo_getnew_theta (deltatheta, soeo%Nocc, soeo%Nact, soeo%oldtheta, decomp)
call soeo_getnew_Dmo (soeo%oldtheta, soeo%Dmo, decomp)
call soeo_getnew_C (K, soeo%C)

!check if C^T S C = 1
!call mat_init (tmpmat, soeo%Nbast, soeo%Nbast)
!
!write (decomp%lupri, *) 'S ='
!call mat_print (soeo%S, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
!write (decomp%lupri, *) 'C ='
!call mat_print (soeo%C, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
!call util_AO_to_MO_2 (soeo%S, soeo%C, soeo%S, tmpmat, .true.)
!write (decomp%lupri, *) 'identity?'
!call mat_print (tmpmat, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
!
!call mat_free (tmpmat)
!

call util_MO_to_AO_2 (soeo%S, soeo%C, soeo%Dmo, soeo%Dao, .false.)
!call soeo_getnew_Dao (soeo%Dmo, soeo%C, soeo%S, K, soeo%Dao)
call FCK_get_Fock (soeo%Dao, soeo%Fao, soeo%Etotal)
soeo%dE = soeo%Etotal - soeo%dE
call util_AO_to_MO_2 (soeo%S, soeo%C, soeo%Fao, soeo%Fmo, .true.)
!call mat_diag_f (soeo%Fao, soeo%S, orbE, soeo%C)
!call mat_zero (soeo%Fmo)
!do i=1,soeo%Nbast
!  tmp(1) = orbE(i)
!  if (decomp%cfg_unres) then
!    tmp(2) = orbE(soeo%Nbast+i)
!  endif
!  call mat_create_ab_elms (i, i, tmp, soeo%Fmo)
!enddo
ndim = size(soeo%nfirst)
do i=1,ndim
  if (i <= soeo%Nact) then
    call soeo_get_nfirst(soeo%oldtheta, i, 'a', soeo%nfirst(i), decomp)
  else
    call soeo_get_nfirst(soeo%oldtheta, i-soeo%Nact, 'b', soeo%nfirst(i), decomp)
  endif
  do j=1,ndim
    if (i<=soeo%Nact) then
      if (j<=soeo%Nact) then
        call soeo_get_nsecond(soeo%oldtheta, i, j, 'a', 'a', soeo%nsecond(i,j), decomp)
      else
        call soeo_get_nsecond(soeo%oldtheta, i, j-soeo%Nact, 'a', 'b', soeo%nsecond(i,j), decomp)
      endif
    else
      if (j<=soeo%Nact) then
        call soeo_get_nsecond(soeo%oldtheta, i-soeo%Nact, j, 'b', 'a', soeo%nsecond(i,j), decomp)
      else
        call soeo_get_nsecond(soeo%oldtheta, i-soeo%Nact, j-soeo%Nact, 'b', 'b', soeo%nsecond(i,j), decomp)
      endif
    endif
  enddo
enddo

deallocate (orbE, tmp)

end subroutine soeoItem_stepforth
!=======================================================================

!> brief Rejects the new matrices and replaces them by old ones
!> author C. Nygaard
!> date 2010.07.14
!> param soeo The structure containing the matrices
!=======================================================================
subroutine soeoItem_stepback (soeo)

implicit none

type(soeoItem), intent(inout) :: soeo

integer                       :: ndim, i, j

soeo%Dao      = soeo%old_Dao
soeo%Fao      = soeo%old_Fao
soeo%C        = soeo%old_C
soeo%Dmo      = soeo%old_Dmo
soeo%Fmo      = soeo%old_Fmo
soeo%oldtheta = soeo%old_oldtheta

ndim = size(soeo%nfirst)
do i=1,ndim
  soeo%nfirst(i) = soeo%old_nfirst(i)
  do j=1,ndim
    soeo%nsecond(i,j) = soeo%old_nsecond(i,j)
  enddo
enddo

soeo%Etotal = soeo%Etotal - soeo%dE

end subroutine soeoItem_stepback
!=======================================================================

!> brief Updates the active space and matrices dependent on it
!> author C. Nygaard
!> date 2010.07.08
!> param newNocc Size of the new occupied space
!> param newNact Size of the new active space
!> param soeo Contains information for the SOEO-calculation
!=======================================================================
subroutine soeoItem_update_space (newNocc, newNact, decomp, soeo)

implicit none

integer, intent(in)           :: newNocc, newNact
type(soeoItem), intent(inout) :: soeo
type(decompItem), intent(in)  :: decomp

integer                       :: ndim, i, j

ndim = size(soeo%nfirst)

do i=1,ndim
  call mat_free (soeo%nfirst(i))
  do j=1,ndim
    call mat_free (soeo%nsecond(i,j))
  enddo
enddo
deallocate (soeo%nfirst,&
          & soeo%nsecond)

do i=1,ndim
  call mat_free (soeo%old_nfirst(i))
  do j=1,ndim
    call mat_free (soeo%old_nsecond(i,j))
  enddo
enddo
deallocate (soeo%old_nfirst,&
          & soeo%old_nsecond)

soeo%Nocc = newNocc
soeo%Nact = newNact

if (decomp%cfg_unres) then
  ndim = 2.0d0*soeo%Nact
else
  ndim = soeo%Nact
endif
allocate (soeo%nfirst(ndim),&
        & soeo%nsecond(ndim,ndim))
do i=1,ndim
  call mat_init(soeo%nfirst(i), soeo%Nbast, soeo%Nbast)
  do j=1,ndim
    call mat_init (soeo%nsecond(i,j), soeo%Nbast, soeo%Nbast)
  enddo
enddo

allocate (soeo%old_nfirst(ndim),&
        & soeo%old_nsecond(ndim,ndim))
do i=1,ndim
  call mat_init(soeo%old_nfirst(i), soeo%Nbast, soeo%Nbast)
  do j=1,ndim
    call mat_init (soeo%old_nsecond(i,j), soeo%Nbast, soeo%Nbast)
  enddo
enddo

do i=1,ndim
  if (i <= soeo%Nact) then
    call soeo_get_nfirst(soeo%oldtheta, i, 'a', soeo%nfirst(i), decomp)
  else
    call soeo_get_nfirst(soeo%oldtheta, i-soeo%Nact, 'b', soeo%nfirst(i), decomp)
  endif
  do j=1,ndim
    if (i<=soeo%Nact) then
      if (j<=soeo%Nact) then
        call soeo_get_nsecond(soeo%oldtheta, i, j, 'a', 'a', soeo%nsecond(i,j), decomp)
      else
        call soeo_get_nsecond(soeo%oldtheta, i, j-soeo%Nact, 'a', 'b', soeo%nsecond(i,j), decomp)
      endif
    else
      if (j<=soeo%Nact) then
        call soeo_get_nsecond(soeo%oldtheta, i-soeo%Nact, j, 'b', 'a', soeo%nsecond(i,j), decomp)
      else
        call soeo_get_nsecond(soeo%oldtheta, i-soeo%Nact, j-soeo%Nact, 'b', 'b', soeo%nsecond(i,j), decomp)
      endif
    endif
  enddo
enddo

end subroutine soeoItem_update_space
!=======================================================================

!> brief Deallocates matrices allocated in soeoItem_init
!> author C. Nygaard
!> date 2010.07.07
!> param soeo The structure to be deallocated
!=======================================================================
subroutine soeoItem_shutdown (soeo)

implicit none

type(soeoItem), intent(inout) :: soeo

integer                       :: ndim, i, j

call mat_free (soeo%Dao     )
call mat_free (soeo%Fao     )
call mat_free (soeo%S       )
call mat_free (soeo%C       )
call mat_free (soeo%Dmo     )
call mat_free (soeo%Fmo     )
call mat_free (soeo%oldtheta)
ndim = size(soeo%nfirst)
do i=1,ndim
  call mat_free (soeo%nfirst(i))
  do j=1,ndim
    call mat_free (soeo%nsecond(i,j))
  enddo
enddo
deallocate (soeo%nfirst,&
          & soeo%nsecond)

call mat_free (soeo%old_Dao     )
call mat_free (soeo%old_Fao     )
call mat_free (soeo%old_C       )
call mat_free (soeo%old_Dmo     )
call mat_free (soeo%old_Fmo     )
call mat_free (soeo%old_oldtheta)
do i=1,ndim
  call mat_free (soeo%old_nfirst(i))
  do j=1,ndim
    call mat_free (soeo%old_nsecond(i,j))
  enddo
enddo
deallocate (soeo%old_nfirst,&
          & soeo%old_nsecond)

end subroutine soeoItem_shutdown
!=======================================================================

end module soeo_itemroutines
