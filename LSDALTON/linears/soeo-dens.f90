module soeo_dens

use Matrix_Operations
use matrix_util
use decompMod

implicit none

!Contains:
!  soeo_Dmo_init
!  soeo_oldtheta_init
!  soeo_getnew_theta
!  soeo_getnew_Dmo
!  soeo_getnew_Dao

contains

!> \brief Finds the MO density from the AO density
!> \author C. Nygaard
!> \date 2010
!> \param S The overlap matrix
!> \param C The MO orbital coefficients
!> \param Dao The density matrix in the AO basis
!> \param Dmo The density matrix in the MO basis
!=======================================================================
subroutine soeo_Dmo_init (S, C, Dao, Dmo, decomp)

implicit none

!I/O
type(decompItem), intent(in) :: decomp
type(matrix), intent(inout)  :: Dmo
type(matrix), intent(in)     :: S, C, Dao
!Other
integer                      :: Nbast, i, j
real(realk), allocatable     :: tmp(:)

if (decomp%cfg_unres) then
  allocate (tmp(2))
else
  allocate (tmp(1))
endif

Nbast = Dmo%nrow
call mat_zero (Dmo)

call util_AO_to_MO_2 (S, C, Dao, Dmo, .false.) !Dmo = C^T S Dao S C

!Cleaning up so that zero elements are really zero
do i=1,Nbast
  do j=1,Nbast
    call mat_get_ab_elms (Dmo, i, j, tmp)
    if (tmp(1) < 1.0d-5) then
      tmp(1) = 0.0d0
    endif
    if (decomp%cfg_unres) then
      if (tmp(2) < 1.0d-5) then
        tmp(2) = 0.0d0
      endif
    endif
    call mat_create_ab_elms (i, j, tmp, Dmo)
  enddo
enddo

!do i=1,Nbast
!  tmp = mat_get_elm (Dmo, i, i)
!  if (tmp < 1.0d-5 .and. tmp > -1.0d-5) then
!    call mat_create_elm (i-1, i-1, 1.00d0, Dmo)
!    call mat_create_elm (i, i, 0.90d0, Dmo)
!    exit
!  endif
!enddo

deallocate (tmp)

end subroutine soeo_Dmo_init
!=======================================================================

!> \brief Initializes the vector oldtheta from Dmo
!> \author C. Nygaard
!> \date 2010
!> \param Dmo The density matrix in the MO basis
!> \param oldtheta The starting point for the occupation angles
!
!  oldtheta(i) = acos(sqrt(Dmo(i,i)))
!
!because
!  Dmo(i,i) = cos^2(oldtheta(i))
!
!=======================================================================
subroutine soeo_oldtheta_init (Dmo, oldtheta, decomp)

implicit none

!I/O:
type(decompItem), intent(in) :: decomp
type(matrix), intent(in)     :: Dmo
type(matrix), intent(inout)  :: oldtheta
!Other:
integer                      :: i, Nbast
real(realk), allocatable     :: tmp(:)

Nbast = Dmo%nrow

if (decomp%cfg_unres) then
  allocate (tmp(2))
else
  allocate (tmp(1))
endif

do i=1,Nbast
  call mat_get_ab_elms (Dmo, i, i, tmp)

  if (tmp(1) > 1.0d0) then 
    tmp(1) = 1.0d0
  endif
  if (tmp(1) < 0.0d0) then
    tmp(1) = 0.0d0
  endif
  tmp(1) = dsqrt(tmp(1))
  tmp(1) = dacos(tmp(1))

  if (decomp%cfg_unres) then
    if (tmp(2) > 1.0d0) then 
      tmp(2) = 1.0d0
    endif
    if (tmp(2) < 0.0d0) then
      tmp(2) = 0.0d0
    endif
    tmp(2) = dsqrt(tmp(2))
    tmp(2) = dacos(tmp(2))
  endif

  call mat_create_ab_elms (i, 1, tmp, oldtheta)
enddo

deallocate (tmp)

end subroutine soeo_oldtheta_init
!=======================================================================

!> \brief Updates the occupation angles
!> \author C. Nygaard
!> \date 2010
!> \param deltatheta The change in occupation angles
!> \param Nocc The number of occupied MO orbitals
!> \param Nact The number of active orbitals
!> \param oldtheta Starting point for the occupation numbers
!=======================================================================
subroutine soeo_getnew_theta (deltatheta, Nocc, Nact, oldtheta, decomp)

implicit none

!I/O
type(decompItem), intent(in) :: decomp
type(matrix), intent(inout)  :: oldtheta
type(matrix), intent(in)     :: deltatheta
integer, intent(in)          :: Nocc, Nact
!Other
integer                      :: i
real(realk), allocatable     :: old(:), new(:)

if (decomp%cfg_unres) then
  allocate (old(2), new(2))
else
  allocate (old(1), new(1))
endif

do i=Nocc+1,Nocc+Nact
  call mat_get_ab_elms (oldtheta, i, 1, old)
  call mat_get_ab_elms (deltatheta, i-Nocc, 1, new)
  old = old + new
  call mat_create_ab_elms (i, 1, old, oldtheta)
enddo

deallocate (old, new)

end subroutine soeo_getnew_theta
!=======================================================================

!> \brief Updates the MO density matrix
!> \author C. Nygaard
!> \date 2010
!> \param theta The occupation angles
!> \param Dmo The density matrix in the MO basis
!
!  Dmo(i,i) = cos^2 (oldtheta(i))
!
!=======================================================================
subroutine soeo_getnew_Dmo (theta, Dmo, decomp)

implicit none

!I/O
type(decompItem), intent(in) :: decomp
type(matrix), intent(in)     :: theta
type(matrix), intent(inout)  :: Dmo
!Other
integer                      :: Nbast
integer                      :: i
real(realk), allocatable     :: n(:), tmp(:)

call mat_zero (Dmo)
if (decomp%cfg_unres) then
  allocate (n(2), tmp(2))
else
  allocate (n(1), tmp(1))
endif

Nbast = Dmo%nrow

do i=1,Nbast
  call mat_get_ab_elms (theta, i, 1, tmp)
  n(1) = dcos(tmp(1)) * dcos(tmp(1))
  if (decomp%cfg_unres) then
    n(2) = dcos(tmp(2)) * dcos(tmp(2))
  endif
  if (n(1) < 1.0d-5) then
     n(1) = 0.0d0
  elseif (decomp%cfg_unres) then
    if (n(2) < 1.0d-5) then
      n(2) = 0.0d0
    endif
  endif
  call mat_create_ab_elms (i, i, n, Dmo)
enddo

deallocate (n, tmp)

end subroutine soeo_getnew_Dmo
!=======================================================================

!> \brief Updates the AO density matrix
!> \author C. Nygaard
!> \date 2010
!> \param Dmo The density matrix in the MO basis
!> \param C The MO orbital coefficients
!> \param S The overlap matrix
!> \param K The MO orbital rotation matrix
!> \param Dao The density matrix in AO basis
!
!  Dao = C exp(-K) Dmo(theta) exp(K) CT
!      = C exp(-K) CT S C Dmo(theta) CT S C exp(K) CT
!      = exp(-XS) d exp(SX)
!      = d + [d,X]_S + 1/2 [[d,X]_S,X]_S + ...
!
!  X = C K CT
!  d = C Dmo(theta) CT
!
!=======================================================================
subroutine soeo_getnew_Dao (Dmo, C, S, K, Dao)

implicit none

!I/O
type(matrix), intent(in)    :: Dmo, C, S, K
type(matrix), intent(inout) :: Dao
!Other
integer                     :: Nbast, i
real(realk)                 :: fac, n, diffnorm, thresh
type(matrix)                :: X, d, tmp, tmp1, tmp2, diff

Nbast = Dmo%nrow
call mat_init (X, Nbast, Nbast)
call mat_init (d, Nbast, Nbast)
call mat_init (tmp, Nbast, Nbast)
call mat_init (tmp1, Nbast, Nbast)
call mat_init (tmp2, Nbast, Nbast)
call mat_init (diff, Nbast, Nbast)

call util_MO_to_AO_2 (S, C, K, X, .false.) ! X = C K CT
call util_MO_to_AO_2 (S, C, Dmo, d, .false.) ! d = C Dmo CT

Dao = d ; tmp = d ; tmp1 = d
fac = 1.0d0 ; i = 0 ; thresh = 1.0d-15
do
  !This terms commutator:
  !---------------------------------------------------------------------
  call ABCcommutator (Nbast, tmp1, X, S, tmp2)
  tmp1 = tmp2
  !---------------------------------------------------------------------

  !The factor in front and add the term:
  !---------------------------------------------------------------------
  if (i == 0) then
    fac = 1.0d0
  else
    fac = fac * i
  endif
  i = i + 1
  n = 1.0d0 / fac
  call mat_add (1.0d0, Dao, n, tmp2, tmp)
  !---------------------------------------------------------------------

  !Convergence check:
  !---------------------------------------------------------------------
  call mat_add (1.0d0, Dao, -1.0d0, tmp, diff)
  diffnorm = mat_sqnorm2 (diff)
  diffnorm = sqrt(diffnorm)
  if (diffnorm < thresh) exit
  Dao = tmp
  !---------------------------------------------------------------------
enddo

call mat_free (X)
call mat_free (d)
call mat_free (tmp)
call mat_free (tmp1)
call mat_free (tmp2)
call mat_free (diff)

end subroutine soeo_getnew_Dao
!=======================================================================


!> \brief Finds the new MO-coefficients
!> \author C. Nygaard
!> \date Sep. 1. 2010
!> \param C The MO orbital coefficients (inout)
!> \param Kin The MO orbital rotation matrix
!
!> C(i+1) = C(i)exp(-K)
!>        = C(i)(1 - K + 1/2 K^2 - 1/6 K^3 + ...)
!=======================================================================
subroutine soeo_getnew_C (Kin, C)

implicit none

type(matrix), intent(in)    :: Kin
type(matrix), intent(inout) :: C
type(matrix)                :: K, Cstep, tmp
integer                     :: i
real(realk)                 :: fac, n, thresh, stepnorm

call mat_init (K, Kin%nrow, Kin%ncol)
call mat_init (Cstep, C%nrow, C%ncol)
call mat_init (tmp, C%nrow, C%ncol)

K = Kin
Cstep = C
i = 1
thresh = 1.0d-15
do
  !The factor in front:
  !---------------------------------------------------------------------
  if (i == 1) then
    fac = 1.0d0
  else
    fac = fac * i
  endif
  i = i + 1
  n = 1.0d0 / fac
  !---------------------------------------------------------------------

  !This terms product:
  !---------------------------------------------------------------------
  call mat_mul (Cstep, K, 'n', 'n', -1.0d0, 0.0d0, tmp)
  Cstep = tmp
  !---------------------------------------------------------------------

  !Add the term:
  !---------------------------------------------------------------------
  call mat_daxpy (n, Cstep, C)
  !---------------------------------------------------------------------

  !Convergence check:
  !---------------------------------------------------------------------
  call mat_scal (n, tmp)
  stepnorm = mat_sqnorm2 (tmp)
  stepnorm = dsqrt(stepnorm)
  if (stepnorm < thresh) exit
  !---------------------------------------------------------------------
enddo

call mat_free (K)
call mat_free (Cstep)
call mat_free (tmp)

end subroutine soeo_getnew_C
!=======================================================================

end module soeo_dens
