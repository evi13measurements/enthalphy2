module soeo_transform

use soeo_typedef
use soeo_matop
use decompMod

implicit none

!contains:
!  soeo_get_nfirst
!  soeo_get_nsecond
!  soeo_get_gradient
!  soeo_linear_transform
!  soeo_get_GmoL
!  soeo_get_Epred

contains

!> \brief Calculates first derivative of Dmo wrt the occupation angles
!> \author C. Nygaard
!> \date 2010
!> \param oldtheta point where derivative is calculated
!> \param p Which occupation angle to derive wrt
!> \param nfirst The derivative
!
!  nfirst(p)_ij = -2 * delta_ij * delta_ip
!                    * sin(oldtheta(p)) * cos(oldtheta(p))
!
!=======================================================================
subroutine soeo_get_nfirst (oldtheta, p, part, nfirst, decomp)

implicit none

!I/O:
type(decompItem), intent(in) :: decomp
integer, intent(in)          :: p
type(matrix), intent(in)     :: oldtheta
type(matrix), intent(inout)  :: nfirst
character, intent(in)        :: part

!Other:
real(realk), allocatable     :: tmp1(:), tmp2(:)

if (decomp%cfg_unres) then
  allocate (tmp1(2), tmp2(2))
else
  allocate (tmp1(1), tmp2(1))
endif

!Initializations:
call mat_zero (nfirst)

call mat_get_ab_elms (oldtheta, p, 1, tmp1)
tmp2 = 0.0d0
if (part=='a' .or. part=='A') then !alpha-part of the matrix - elms(:)
tmp2(1) = -2.0d0 * dsin(tmp1(1)) * dcos(tmp1(1))
elseif (part=='b' .or. part=='B') then 
  if (decomp%cfg_unres) then !beta-part of the matrix - elmsb(:)
    tmp2(2) = -2.0d0 * dsin(tmp1(2)) * dcos(tmp1(2))
  else
    call lsquit ('Dimension of nfirst doesnt fit restricted calculation',-1)
  endif
else
  call lsquit ('Wrong argument in soeo_get_nfirst',-1)
endif
call mat_create_ab_elms (p, p, tmp2, nfirst)

deallocate (tmp1, tmp2)

end subroutine soeo_get_nfirst
!=======================================================================

!> \brief Calculates the second derivative of Dmo wrt the occupation angles
!> \author C. Nygaard
!> \date 2010
!> \param oldtheta Where the derivative should be calculated
!> \param p Which occupation angle to make first derivation wrt
!> \param q Which occupation angle so make second edrivation wrt
!> \param nsecons The derivative
!
!  nsecond(i,j)_pq = -2 * delta_ij * delta_ip * delta_pq
!                       * (cos^2 (oldtheta(q)) - sin^2 (oldtheta(q)))
!
!=======================================================================
subroutine soeo_get_nsecond (oldtheta, p, q, partp, partq, nsecond, decomp)

implicit none

!Input/output
type(decompItem), intent(in) :: decomp
integer, intent(in)          :: p, q
type(matrix), intent(in)     :: oldtheta
type(matrix), intent(inout)  :: nsecond
character, intent(in)        :: partp, partq

!Other:
real(realk), allocatable     :: tmp1(:), tmp2(:)

if (decomp%cfg_unres) then
  allocate (tmp1(2), tmp2(2))
else
  allocate (tmp1(1), tmp2(1))
endif

!Initializing:
call mat_zero (nsecond)

if (partp == partq) then
  if (p == q) then
    call mat_get_ab_elms (oldtheta, q, 1, tmp1)
    tmp2 = 0.0d0
    if (partq == 'a' .or. partq == 'A') then
      tmp2(1) = -2.0d0 * dcos(tmp1(1)) * dcos(tmp1(1))
      tmp2(1) = tmp2(1) + 2.0d0 * dsin(tmp1(1)) * dsin(tmp1(1))
    elseif (partq == 'b' .or. partq == 'B') then
      if (decomp%cfg_unres) then
        tmp2(2) = -2.0d0 * dcos(tmp1(2)) * dcos(tmp1(2))
        tmp2(2) = tmp2(2) + 2.0d0 * dsin(tmp1(2)) * dsin(tmp1(2))
      else
        call lsquit ('Dimension of nsecond doesnt fit restricted calculation',-1)
      endif
    else
      call lsquit ('Wrong argument in soeo_get_nsecond',-1)
    endif
    call mat_create_ab_elms (p, p, tmp2, nsecond)
  endif
endif

deallocate (tmp1, tmp2)

end subroutine soeo_get_nsecond
!=======================================================================



!> \brief Calculates left-hand-side (lhs) of the SOEO equatios
!> \author C. Nygaard
!> \date 2010
!> \param Fmo The Fock matrix in the MO basis
!> \param Dmo The density matrix in the MO basis
!> \param nfirst The first derivative of Dmo wrt the occupation angles
!> \param grad_m Matrix part of the lhs
!> \param grad_v Vector part of the lhs
!
!  grad_m    = -4*[Fmo, Dmo]
!  grad_v(p) = 2*Tr(Fmo*nfirst(p))
!
!=======================================================================
subroutine soeo_get_gradient (soeo, grad_m, grad_v, decomp)

implicit none

!Input/output:
type(decompItem), intent(in) :: decomp
type(soeoItem), intent(in)   :: soeo
type(matrix), intent(inout)  :: grad_m, grad_v
!Other:
integer                      :: i
real(realk), allocatable     :: tmp(:)

if (decomp%cfg_unres) then
  allocate (tmp(2))
else
  allocate (tmp(1))
endif

call mat_zero (grad_m)
call mat_zero (grad_v)

!The matrix-part:
!-----------------------------------------------------------------------
call commutator (-2.0d0, soeo%Fmo, soeo%Dmo, "n", "n", grad_m)
if (.not. decomp%cfg_unres) then
  call mat_scal (2.0d0, grad_m)
endif
!-----------------------------------------------------------------------

!The vector-part:
!-----------------------------------------------------------------------
do i=1,soeo%Nact
  tmp(1) = mat_TrAB (soeo%Fmo, soeo%nfirst(i))
  if (decomp%cfg_unres) then
    tmp(2) = mat_TrAB (soeo%Fmo, soeo%nfirst(i+soeo%Nact))
  else
    tmp = 2.0d0*tmp
  endif
  call mat_create_ab_elms (i, 1, tmp, grad_v)
enddo
!-----------------------------------------------------------------------

deallocate (tmp)

end subroutine soeo_get_gradient
!=======================================================================



!> \brief Calculates right-hand-side (rhs) of the SOEO equatios
!> \author C. Nygaard
!> \date 2010
!> \param Fmo The Fock matrix in the MO basis
!> \param Dmo The density matrix in the MO basis
!> \param C The MO orbital coefficients
!> \param nfirst The first derivative of Dmo wrt the occupation angles
!> \param nsecond The second derivative of Dmo wrt the occupation angles
!> \param S The overlap matrix
!> \param Dao The density matrix in AO basis
!> \param b_m Matrix part of trial vector
!> \param b_v Vector part of trial vector
!> \param sigma_m Matrix part of the rhs
!> \param sigma_v Vector part of the rhs
!
!>  sigma_m    = -2*[[b_m,Fmo],Dmo] - 2*[Fmo,[Dmo,b_m]]
!>               - 4*sum_(r)(b_v(r)*[Fmo,nfirst(r)])
!>               - 8*[gmo*L,Dmo]
!>  sigma_v(p) = 2sum_(r)(b_v(r)*Tr(Fmo*nsecond(p,r)))
!>               - 2Tr([Fmo,b_m]nfirst(p))
!>               + 4*Tr(gmo*L*nfirst(p))
!>
!>  where L = [Dmo,b_m] + sum_(r)(b_v(r)*nfirst(r))
!
!(nfirst(i) and nsecond(i,j) are matrices, not scalars)
!=======================================================================
subroutine soeo_linear_transform (soeo, b_m, b_v, sigma_m, sigma_v, decomp)

implicit none

!I/O
type(decompItem), intent(in) :: decomp
type(soeoItem), intent(in)   :: soeo
type(matrix), intent(in)     :: b_m, b_v
type(matrix), intent(inout)  :: sigma_m, sigma_v

!Other
integer                      :: ndim, i, j
type(matrix)                 :: L, GmoL, commDb, sumn
!Temporary stuff:
type(matrix)                 :: tmp_m
real(realk)                  :: val1, val2
real(realk), allocatable     :: tmp(:), tmp1(:)

!print *, 'soeo_linear_transform called and started'

if (decomp%cfg_unres) then
  allocate (tmp(2), tmp1(2))
else
  allocate (tmp(1), tmp1(1))
endif

!Initializations:
call mat_init (L, soeo%Nbast, soeo%Nbast)
call mat_init (GmoL, soeo%Nbast, soeo%Nbast)
call mat_init (commDb, soeo%Nbast, soeo%Nbast)
call mat_init (sumn, soeo%Nbast, soeo%Nbast)
call mat_init (tmp_m, soeo%Nbast, soeo%Nbast)

ndim = size(soeo%nfirst)

!Some usefull matrices
!-----------------------------------------------------------------------
call commutator (1.0d0, soeo%Dmo, b_m, "n", "n", commDb) ! [D,b_m]

call mat_zero (sumn)                                        ! sum_(r)(b_v(r) nfirst(r))
do i=1,soeo%Nact                                            !
  call mat_get_ab_elms (b_v, i, 1, tmp)                     !
  call mat_ab_daxpy (tmp, soeo%nfirst(i), sumn)             !
  if (decomp%cfg_unres) then                                ! for unrestricted
    call mat_ab_daxpy (tmp, soeo%nfirst(i+soeo%Nact), sumn) !  the sum is also
  endif                                                     !  over spin
enddo                                                       !

!Creation of L and G(L)
call mat_add (1.0d0, commDb, 1.0d0, sumn, L)
call soeo_get_GmoL (L, soeo%S, soeo%C, soeo%Dao, GmoL, decomp%lupri)
!-----------------------------------------------------------------------

!DEBUGGING

write (decomp%lupri, *)
write (decomp%lupri, *) 'Input in soeo_linear_transform:'
write (decomp%lupri, *) 'b_m ='
call mat_print (b_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'b_v ='
call mat_print (b_v, 1, soeo%Nact, 1, 1, decomp%lupri)
write (decomp%lupri, *)
write (decomp%lupri, *) 'Matrices in use'
write (decomp%lupri, *) 'Dmo ='
call mat_print (soeo%Dmo, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'sumn ='
call mat_print (sumn, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'commDb ='
call mat_print (commDb, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'L ='
call mat_print (L, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'GmoL ='
call mat_print (GmoL, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *)

!END DEBUGGING

!The matrix-part:
!-----------------------------------------------------------------------
!First term: sigma_m = 2*[[b_m,Fmo],Dmo]
! unres => *0.5
call commutator (1.0d0, b_m, soeo%Fmo, "n", "n", tmp_m)
call commutator (1.0d0, tmp_m, soeo%Dmo, "n", "n", sigma_m)
if (.not. decomp%cfg_unres) then
  call mat_scal (2.0d0, sigma_m)
endif

!Second term: sigma_m = sigma_m + 2*[Fmo,[Dmo,b_m]]
! unres => *0.5
call commutator (1.0d0, soeo%Fmo, commDb, "n", "n", tmp_m)
if (.not. decomp%cfg_unres) then
  call mat_scal (2.0d0, tmp_m)
endif
call mat_daxpy (1.0d0, tmp_m, sigma_m)

!Third term: sigma_m = sigma_m + 4*sum_(r)(b_v(r)*[Fmo,nfirst(r)])
!                    = sigma_m + 4*[Fmo,sumn]
! unres => *0.5
call commutator (1.0d0, soeo%Fmo, sumn, "n", "n", tmp_m)
if (.not. decomp%cfg_unres) then
  call mat_scal (2.0d0, tmp_m)
endif
call mat_daxpy (2.0d0, tmp_m, sigma_m)

!Fourth term: sigma_m = sigma_m + 8*[GmoL,Dmo]
! unres => *0.25
call commutator (1.0d0, GmoL, soeo%Dmo, "n", "n", tmp_m)
if (.not. decomp%cfg_unres) then
  call mat_scal (4.0d0, tmp_m)
endif
call mat_daxpy (2.0d0, tmp_m, sigma_m)

call mat_scal (-1.0d0, sigma_m)
!-----------------------------------------------------------------------

!The vector-part:
!-----------------------------------------------------------------------
call commutator (1.0d0, soeo%Fmo, b_m, "n", "n", tmp_m) ! [F,b_m]

write (decomp%lupri, *) '[Fmo, b_m] ='
call mat_print (tmp_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)

do i=1,ndim

write (decomp%lupri, *) 'i =', i

  tmp = 0.0d0
  call mat_zero (sumn)                               ! sumn = sum_(r)(b_v(r) nsecond(p,r))
  if (i <= soeo%Nact) then                           !      = b_v(p)nsecond(p,p)
    call mat_get_ab_elms (b_v, i, 1, tmp)            !
  else                                               !
    call mat_get_ab_elms (b_v, i-soeo%Nact, 1, tmp)  ! (also sum over spin)
  endif                                              !
  call mat_ab_daxpy (tmp, soeo%nsecond(i,i), sumn)   !

write (decomp%lupri, *) 'sumn = nsecond(i,i) = 0 ?'
call mat_print (sumn, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)

  !First term: sigma_v(i) = 2*sum_(r)(b_v(r)*Tr(Fmo*nsecond(p,r)))
  ! unres => *0.5
!  call mat_TrAB_ab (soeo%Fmo, sumn, tmp1)
  val1 = mat_TrAB (soeo%Fmo, sumn)
  if (.not. decomp%cfg_unres) then
    val1 = 2.0d0 * val1
  endif
  if (i <= soeo%Nact) then
    tmp(1) = val1
  else
    tmp(2) = val1
  endif
write (decomp%lupri, *) 'Tr(Fmo sumn) =', val1

  !Second term: sigma_v(i) = sigma_v(i) - 2*Tr([Fmo,b_m]*nfirst(i))
!  call mat_TrAB_ab (tmp_m, soeo%nfirst(i), tmp1)
  val1 = mat_TrAB (tmp_m, soeo%nfirst(i))
  if (.not. decomp%cfg_unres) then
    val1 = 2.0d0 * val1
  endif
  if (i <= soeo%Nact) then
    tmp(1) = tmp(1) - val1
  else
    tmp(2) = tmp(2) - val1
  endif
write (decomp%lupri, *) 'nfirst(i) ='
call mat_print (soeo%nfirst(i), 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'Tr([Fmo,b_m]*nfirst(i)) =', val1

  !Third term: sigma_v(i) = sigma_v(i) + 4*Tr(GmoL*nfirst(i))
!  call mat_TrAB_ab (GmoL, soeo%nfirst(i), tmp1)
  val1 = mat_TrAB (GmoL, soeo%nfirst(i))
  if (.not. decomp%cfg_unres) then
    val1 = 4.0d0 * val1
  endif
  if (i <= soeo%Nact) then
    tmp(1) = tmp(1) + val1
  else
    tmp(2) = tmp(2) + val1
  endif
write (decomp%lupri, *) 'Tr(GmoL nfirst(i)) =', val1

  if (i <= soeo%Nact) then
    call mat_get_ab_elms (sigma_v, i, 1, tmp1)
    tmp1(1) = tmp(1)
    call mat_create_ab_elms (i, 1, tmp1, sigma_v)
  else
    call mat_get_ab_elms (sigma_v, i-soeo%Nact, 1, tmp1)
    tmp1(2) = tmp(2)
    call mat_create_ab_elms (i-soeo%Nact, 1, tmp1, sigma_v)
  endif
enddo
!-----------------------------------------------------------------------

!Finalizations:
call mat_free (L)
call mat_free (GmoL)
call mat_free (commDb)
call mat_free (sumn)
call mat_free (tmp_m)

deallocate (tmp, tmp1)

end subroutine soeo_linear_transform
!=======================================================================

!> \brief Gets G(L) in the MO basis (see di_get_GAOL in dalton_interface.f90)
!> \author C. Nygaard
!> \date 2010
!> \param Lmo The matrix L in the MO basis
!> \param S The overlap matrix
!> \param C The MO orbital coefficients
!> \param Dao The density matrix in the AO basis
!> \param GmoL The wanted matrix G(L)
!=======================================================================
subroutine soeo_get_GmoL (Lmo, S, C, Dao, GmoL, lupri)

implicit none

!I/O
type(matrix), intent(in)    :: Lmo, S, C, Dao
integer, intent(in)         :: lupri
type(matrix), intent(inout) :: GmoL !Intent(out)
!Other
logical, external           :: do_dft
type(matrix)                :: Lao, GaoL, tmp
integer                     :: Nbast

Nbast = Lmo%nrow

call mat_init (Lao , Nbast, Nbast)
call mat_init (GaoL, Nbast, Nbast)
call mat_init (tmp , Nbast, Nbast)

write (lupri, *) 'C ='
call mat_print (C, 1, Nbast, 1, Nbast, lupri)
write (lupri, *) 'Lmo ='
call mat_print (Lmo, 1, Nbast, 1, Nbast, lupri)

!Transform L to the AO-basis:
!-----------------------------------------------------------------------
call util_MO_to_AO_2 (S, C, Lmo, Lao, .false.)
!-----------------------------------------------------------------------

write (lupri, *) 'Lao ='
call mat_print (Lao, 1, Nbast, 1, Nbast, lupri)

!Get G(L) in the AO-basis:
!-----------------------------------------------------------------------
call di_get_GAOL(Lao, GaoL)
if (do_dft()) then
!WARNING THIS SHOULD BE MODIFIED TO WORK WITH LSDALTON ONLY
  call di_get_sigma_xc_cont(Lao,Dao,S,tmp)
  call mat_daxpy (1.0d0, tmp, GaoL)
endif
!-----------------------------------------------------------------------

write (lupri, *) 'GaoL ='
call mat_print (GaoL, 1, Nbast, 1, Nbast, lupri)

!Transform G(L) back to the MO-basis:
!-----------------------------------------------------------------------
call util_AO_to_MO_2 (S, C, GaoL, GmoL, .true.)
!-----------------------------------------------------------------------

write (lupri, *) 'GmoL ='
call mat_print (GmoL, 1, Nbast, 1, Nbast, lupri)

call mat_free (Lao )
call mat_free (GaoL)
call mat_free (tmp )

end subroutine soeo_get_GmoL
!=======================================================================

!> \brief Finds the predicted energy change
!> \author C. Nygaard
!> \date 2010-07-01
!> \param Fmo The MO Fock matrix
!> \param Dmo The MO density matrix
!> \param nfirst The first derivative of Dmo wrt theta
!> \param nsecond The second derivative of Dmo wrt theta
!> \param K The orbital rotation matrix
!> \param deltatheta The change in occupation angles
!> \param S The AO overlap matrix
!> \param C The MO coefficient matrix
!> \param Dao The AO density matrix
!> \param Epred The predicted energy change
!
!> Epred = -1/2 * (K,deltatheta)^T H(K,deltatheta)
!
!=======================================================================
subroutine soeo_get_Epred (soeo, K, deltatheta, decomp)

implicit none

type(decompItem), intent(in)  :: decomp
type(soeoItem), intent(inout) :: soeo
type(matrix), intent(in)      :: K, deltatheta

type(matrix)                  :: sigma_m, sigma_v
type(matrix)                  :: grad_m, grad_v

call mat_init (sigma_m, soeo%Nbast, soeo%Nbast)
call mat_init (sigma_v, soeo%Nact, 1)
call mat_init (grad_m, soeo%Nbast, soeo%Nbast)
call mat_init (grad_v, soeo%Nact, 1)

soeo%dEpred = 0.0d0

call soeo_linear_transform (soeo, K, deltatheta, sigma_m, sigma_v, decomp) !H(v)

!soeo%dEpred = soeo_dotproduct (K, deltatheta, sigma_m, sigma_v)

!Why does it work with -1.0 and not with -0.5 as it should?
!Try making dEpred2 = v*g + 1/2 v*H(v) and see what that gives me!
!soeo%dEpred = -0.5d0 * soeo%dEpred
!soeo%dEpred = -1.0d0 * soeo%dEpred

call soeo_get_gradient (soeo, grad_m, grad_v, decomp) !g
call soeo_daxpy (0.5d0, sigma_m, sigma_v, grad_m, grad_v) !g + 1/2 H(v)
soeo%dEpred = soeo_dotproduct (K, deltatheta, grad_m, grad_v) !v*(g + 1/2 H(v))

call mat_free (sigma_m)
call mat_free (sigma_v)
call mat_free (grad_m)
call mat_free (grad_v)

end subroutine soeo_get_Epred
!=======================================================================

end module soeo_transform
