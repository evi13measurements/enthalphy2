module soeo_debug

use soeo_matop
use soeo_transform
use soeo_dens
use soeo_redspace
use Fock_evaluator
use soeo_itemroutines

implicit none

contains

!> \brief Gets full Hessian
!> \author C. Nygaard
!> \date 2010
!> \param Fmo The Fock matrix in MO basis
!> \param Dmo The density matrix in MO basis
!> \param C The MO orbital coefficients
!> \param nfirst The first derivative of Dmo wrt the occupation angles
!> \param nsecond The second derivative of Dmo wrt the occupation angles
!> \param S The overlap matrix
!> \param Dao The density matrix in AO basis
!> \param Nact Number of Active MO orbitals
!> \param oldtheta Where the Hessian is to be made
!> \param hessian The full Hessian
!=======================================================================
subroutine soeo_get_full_hessian (soeo, hessian, decomp)

implicit none

!I/O
type(soeoItem), intent(in)   :: soeo
type(decompItem), intent(in) :: decomp
real(realk), intent(inout)   :: hessian(:,:)
!Other
integer                      :: i, j, ib, Nbast, hessdim, Kdim
type(matrix)                 :: e_mat, e_mat_vec, tmpmat, hessvec
type(matrix)                 :: e_vec, tmpvec
real(realk)                  :: tmp
real(realk), allocatable     :: elm(:)
integer                      :: ndim, counter

!Initialization
!-----------------------------------------------------------------------
if (decomp%cfg_unres) then
  allocate (elm(2))
else
  allocate (elm(1))
endif

hessdim = size(hessian(:,1))
if (decomp%cfg_unres) then
  hessdim = hessdim / 2
endif
Kdim = hessdim - soeo%Nact
ndim = size(soeo%nfirst)

call mat_init (e_mat, soeo%Nbast, soeo%Nbast)
call mat_init (e_mat_vec, Kdim, 1)
call mat_init (tmpmat, soeo%Nbast, soeo%Nbast)
call mat_init (hessvec, Kdim, 1)
call mat_init (e_vec, soeo%Nact, 1)
call mat_init (tmpvec, soeo%Nact, 1)

call mat_zero (e_mat)
call mat_zero (e_mat_vec)
call mat_zero (e_vec)
hessian = 0.0d0
!-----------------------------------------------------------------------

!write (decomp%lupri, *) 'soeo_get_full_hessian started'
!write (decomp%lupri, *) 'hessdim =', hessdim, 'Kdim =', Kdim

!Setting up the matrix-part
!-----------------------------------------------------------------------
do i=1,2*Kdim
  if (i>Kdim) then
    if (.not. decomp%cfg_unres) exit
    ib = i-Kdim + hessdim
  else
    ib = i
  endif

  if (decomp%cfg_unres) then
    if (i > Kdim) then
      call mat_create_ab_elms (i-Kdim, 1, (/ 0.0d0, 1.0d0 /), e_mat_vec)
    else
      call mat_create_ab_elms (i, 1, (/ 1.0d0, 0.0d0 /), e_mat_vec)
    endif
  else
    call mat_create_elm (i, 1, 1.0d0, e_mat_vec)
  endif
  call mat_vec_to_mat ('a', e_mat_vec, e_mat)

  call soeo_linear_transform (soeo, e_mat, e_vec, tmpmat, tmpvec, decomp)

  call mat_to_vec ('a', tmpmat, hessvec)
  do j=1,Kdim
    if (decomp%cfg_unres) then
      call mat_get_ab_elms (hessvec, j, 1, elm)
      hessian(j,ib) = elm(1)
      hessian(hessdim+j,ib) = elm(2)
    else
      call mat_get_elm (hessvec, j, 1, elm(1))
      hessian(j,ib) = elm(1)
    endif
  enddo
  do j=1,soeo%Nact
    if (decomp%cfg_unres) then
      call mat_get_ab_elms (tmpvec, j, 1, elm)
      hessian(Kdim+j,ib) = elm(1)
      hessian(hessdim+Kdim+j,ib) = elm(2)
    else
      call mat_get_elm (tmpvec, j, 1, elm(1))
      hessian(Kdim+j,ib) = elm(1)
    endif
  enddo
  call mat_zero (e_mat)
  call mat_zero (e_mat_vec)
enddo
!-----------------------------------------------------------------------

!Setting up the vector-part
!-----------------------------------------------------------------------
do i=1,ndim
  if (i>soeo%Nact) then
    if (.not. decomp%cfg_unres) exit
    ib = i - soeo%Nact + hessdim + Kdim
  else
    ib = i + Kdim
  endif

  call mat_zero (e_vec)
  if (i > soeo%Nact) then
    call mat_create_ab_elms (i-soeo%Nact, 1, (/ 0.0d0, 1.0d0 /), e_vec)
  else
    call mat_create_ab_elms (i, 1, (/ 1.0d0, 0.0d0 /), e_vec)
  endif

  call soeo_linear_transform (soeo, e_mat, e_vec, tmpmat, tmpvec, decomp)

  call mat_to_vec ('a', tmpmat, hessvec)
  do j=1,Kdim
    if (decomp%cfg_unres) then
      call mat_get_ab_elms (hessvec, j, 1, elm)
      hessian(j,ib) = elm(1)
      hessian(hessdim+j,ib) = elm(2)
    else
      call mat_get_elm (hessvec, j, 1, elm(1))
      hessian(j,ib) = elm(1)
    endif
  enddo
  do j=1,soeo%Nact
    if (decomp%cfg_unres) then
      call mat_get_ab_elms (tmpvec, j, 1, elm)
      hessian(Kdim+j,ib) = elm(1)
      hessian(hessdim+Kdim+j,ib) = elm(2)
    else
      call mat_get_elm (tmpvec, j, 1, elm(1))
      hessian(Kdim+j,ib) = elm(1)
    endif
  enddo
enddo
!-----------------------------------------------------------------------

!Checking for symmetry
!-----------------------------------------------------------------------
write (decomp%lupri, *)
write (decomp%lupri, *) "Is the full hessian symmetric?"
if (decomp%cfg_unres) then
  write (decomp%lupri, *) "(alpha-part only)"
endif
write (decomp%lupri, *) "======================================"

!write (decomp%lupri, *) "Full hessian ="
!call output (hessian, 1, hessdim, 1, hessdim,&
!           & hessdim, hessdim, 1, decomp%lupri)

        
counter = 0
do i=1,hessdim
  do j=1,hessdim
    tmp = hessian(i,j) - hessian(j,i)
    if (tmp > 1.0d-10 .or. tmp < -1.0d-10) then
      counter = counter + 1
    endif
  enddo
enddo
if (counter /= 0) then
  write (decomp%lupri, *) "No, the Hessian is NOT SYMMETRIC!!!"
  write (decomp%lupri, *) "Number of elements in hessian (not diagonal):",&
                  & hessdim*hessdim-hessdim
  write (decomp%lupri, *) "Number of different elements =", counter
  write (decomp%lupri, *) "======================================"
else
  write (decomp%lupri, *) "Yes, the Hessian is SYMMETRIC!"
endif

counter = 0
do i=1,Kdim
  do j=1,Kdim
    tmp = hessian(i,j) - hessian(j,i)
    if (tmp < -1.0d-5 .or. tmp > 1.0d-5) then
      counter = counter + 1
    endif
  enddo
enddo
if (counter /= 0) then
  write (decomp%lupri, *) "KK-part of the hessian is NOT SYMMETRIC!!!"
  write (decomp%lupri, *) "Number of elements in KK-part (not diagonal):",&
                  & Kdim*Kdim-Kdim
  write (decomp%lupri, *) "Number of different elements in KK-part:",&
                  & counter
endif

counter = 0
do i=Kdim+1,hessdim
  do j=Kdim+1,hessdim
    tmp = hessian(i,j) - hessian(j,i)
    if (tmp > 1.0d-5 .or. tmp < -1.0d-5) then
      counter = counter + 1
    endif
  enddo
enddo
if (counter /= 0) then
  write (decomp%lupri, *) "tt-part of the hessian is NOT SYMMETRIC!!!"
  write (decomp%lupri, *) "Number of elements in tt-part (not diagonal):",&
                  & soeo%Nact*soeo%Nact-soeo%Nact
  write (decomp%lupri, *) "Number of different elements in tt-part:",&
                  & counter
endif

counter = 0
do i=1,Kdim
  do j=Kdim+1,hessdim
    tmp = hessian(i,j) - hessian(j,i)
    if (tmp < -1.0d-5 .or. tmp > 1.0d-5) then
      counter = counter + 1
    endif
  enddo
enddo
if (counter /= 0) then
  write (decomp%lupri, *) "Kt-part of the hessian is NOT SYMMETRIC!!!"
  write (decomp%lupri, *) "Number of elements in Kt-part:",&
                  & Kdim*soeo%Nact
  write (decomp%lupri, *) "Number of different elements in Kt-part:",&
                  & counter
endif
write (decomp%lupri, *) "======================================"
!-----------------------------------------------------------------------

!Closing
!-----------------------------------------------------------------------
call mat_free (e_mat)
call mat_free (e_mat_vec)
call mat_free (tmpmat)
call mat_free (hessvec)
call mat_free (e_vec)
call mat_free (tmpvec)

deallocate (elm)
!-----------------------------------------------------------------------

end subroutine soeo_get_full_hessian
!=======================================================================

!> \brief Gets full gradient and Hessian by Finite Difference
!> \author C. Nygaard
!> \date 2010-10-05
!=======================================================================
subroutine soeo_finite_difference (soeo, gradient, hessian, decomp)

implicit none

type(soeoItem), intent(in)   :: soeo
type(decompItem), intent(in) :: decomp
real(realk), intent(inout)   :: gradient(:), hessian(:,:)

type(matrix)                 :: theta, K
real(realk)                  :: delta
real(realk)                  :: Ep1, Em1, Ep2, Em2
real(realk)                  :: Ep1p1, Ep1m1, Em1p1, Em1m1
real(realk)                  :: Ep2p2, Ep2m2, Em2p2, Em2m2
real(realk)                  :: G1, G2
integer                      :: hessdim, Kdim, ndim
integer                      :: i, j, row, col

hessdim = size(gradient)
if (decomp%cfg_unres) then
  hessdim = hessdim / 2
endif
Kdim = hessdim - soeo%Nact
ndim = size(soeo%nfirst)

call mat_init (theta, soeo%Nbast, 1)
call mat_init (K, soeo%Nbast, soeo%Nbast)

delta = 0.05d0

!Matrix-part
!-----------------------------------------------------------------------
do i=1,2*Kdim
  if (i>Kdim) then
    if (.not. decomp%cfg_unres) exit
    row = hessdim + i - Kdim
  else
    row = i
  endif

  !Energies with change in i'th variable in K
  !---------------------------------------------------------------------
  theta = soeo%oldtheta

  call mat_zero (K)
  call soeo_fd_add_delta_to_K (K, delta, i, decomp%cfg_unres)
  call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1)

  call mat_zero (K)
  call soeo_fd_add_delta_to_K (K, -delta, i, decomp%cfg_unres)
  call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1)

  call mat_zero (K)
  call soeo_fd_add_delta_to_K (K, 2.0d0*delta, i, decomp%cfg_unres)
  call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2)

  call mat_zero (K)
  call soeo_fd_add_delta_to_K (K, -2.0d0*delta, i, decomp%cfg_unres)
  call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2)
  !---------------------------------------------------------------------

  !Matrix-part of the gradient
  ! gradient = (8Ep1 - 8Em1 - Ep2 + Em2)/(12delta)
  !---------------------------------------------------------------------
  gradient(row) = (8.0d0*Ep1 - 8.0d0*Em1 - Ep2 + Em2) / (12.0d0*delta)
  !---------------------------------------------------------------------

  !Matrix-matrix-part
  !---------------------------------------------------------------------
  do j=1,2*Kdim
    if (j>Kdim) then
      if (.not. decomp%cfg_unres) exit
      col = hessdim + j - Kdim
    else
      col = j
    endif

    !Energies with change in i'th and j'th variable in K
    !-------------------------------------------------------------------
    theta = soeo%oldtheta

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, i, decomp%cfg_unres)
    call soeo_fd_add_delta_to_K (K, delta, j, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, i, decomp%cfg_unres)
    call soeo_fd_add_delta_to_K (K, -delta, j, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, i, decomp%cfg_unres)
    call soeo_fd_add_delta_to_K (K, delta, j, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, i, decomp%cfg_unres)
    call soeo_fd_add_delta_to_K (K, -delta, j, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0d0*delta, i, decomp%cfg_unres)
    call soeo_fd_add_delta_to_K (K, 2.0d0*delta, j, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0d0*delta, i, decomp%cfg_unres)
    call soeo_fd_add_delta_to_K (K, -2.0d0*delta, j, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2m2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0d0*delta, i, decomp%cfg_unres)
    call soeo_fd_add_delta_to_K (K, 2.0d0*delta, j, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0d0*delta, i, decomp%cfg_unres)
    call soeo_fd_add_delta_to_K (K, -2.0d0*delta, j, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2m2)
    !-------------------------------------------------------------------

    !Matrix-matrix-part of the Hessian
    ! hessian = (16G1 - G2)/(48delta^2)
    !  G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    !  G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    !-------------------------------------------------------------------
    G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    hessian(row,col) = (16.0d0*G1 - G2) / (48.0d0*delta*delta)
    !-------------------------------------------------------------------
  enddo
  !---------------------------------------------------------------------

  !Matrix-vector-part
  do j=1,ndim
    if (j>soeo%Nact) then
      if (.not. decomp%cfg_unres) exit
      col = j - soeo%Nact + hessdim + Kdim
    else
      col = j + Kdim
    endif

    !Energies with change in i'th variable in K
    !                    and j'th variable in theta
    !-------------------------------------------------------------------
    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, i, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, delta, j,&
                                   & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, i, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -delta, j,&
                                   & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, i, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, delta, j,&
                                   & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, i, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -delta, j,&
                                   & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0d0*delta, i, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, 2.0d0*delta, j,&
                                   & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0d0*delta, i, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -2.0d0*delta, j,&
                                   & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2m2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0d0*delta, i, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, 2.0d0*delta, j,&
                                   & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0d0*delta, i, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -2.0d0*delta, j,&
                                   & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2m2)
    !-------------------------------------------------------------------

    !Matrix-vector-part of the Hessian
    ! hessian = (16G1 - G2)/(48delta^2)
    !  G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    !  G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    !-------------------------------------------------------------------
    G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    hessian(row,col) = (16.0d0*G1 - G2) / (48.0d0*delta*delta)
    !-------------------------------------------------------------------
  enddo
enddo
!-----------------------------------------------------------------------

!Vector-part
!-----------------------------------------------------------------------
call mat_zero (K)
do i=1,ndim
  if (i>soeo%Nact) then
    if (.not. decomp%cfg_unres) exit
    row = i - soeo%Nact + hessdim + Kdim
  else
    row = i + Kdim
  endif

  !Energies with change in i'th variable in theta
  !---------------------------------------------------------------------
  call mat_zero (K)

  theta = soeo%oldtheta
  call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
  call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1)

  theta = soeo%oldtheta
  call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
  call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1)

  theta = soeo%oldtheta
  call soeo_fd_add_delta_to_theta (theta, 2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
  call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2)

  theta = soeo%oldtheta
  call soeo_fd_add_delta_to_theta (theta, -2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
  call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2)
  !---------------------------------------------------------------------

  !Vector-part of the gradient
  ! gradient = (8Ep1 - 8Em1 - Ep2 + Em2)/(12delta)
  !---------------------------------------------------------------------
  gradient(row) = (8.0d0*Ep1 - 8.0d0*Em1 - Ep2 + Em2) / (12.0d0*delta)
  !---------------------------------------------------------------------

  !Vector-matrix-part
  !---------------------------------------------------------------------
  do j=1,2*Kdim
    if (j>Kdim) then
      if (.not. decomp%cfg_unres) exit
      col = hessdim + j - Kdim
    else
      col = j
    endif

    !Energies with change in i'th variable in theta
    !                    and j'th variable in K
    !-------------------------------------------------------------------
    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, j, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, j, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, j, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, j, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0d0*delta, j, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, 2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0d0*delta, j, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2m2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0d0*delta, j, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, 2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0d0*delta, j, decomp%cfg_unres)
    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2m2)
    !-------------------------------------------------------------------

    !Vector-matrix-part of the Hessian
    ! hessian = (16G1 - G2)/(48delta^2)
    !  G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    !  G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    !-------------------------------------------------------------------
    G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    hessian(row,col) = (16.0d0*G1 - G2) / (48.0d0*delta*delta)
    !-------------------------------------------------------------------
  enddo
  !---------------------------------------------------------------------

  !Vector-vector-part
  !---------------------------------------------------------------------
  do j=1,ndim
    if (j>soeo%Nact) then
      if (.not. decomp%cfg_unres) exit
      col = j - soeo%Nact + hessdim + Kdim
    else
      col = j + Kdim
    endif

    !Energies with change in i'th and j'th variable in theta
    !-------------------------------------------------------------------
    call mat_zero (K)

    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, delta, j,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1p1)

    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, -delta, j,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep1m1)

    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, delta, j,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1p1)

    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, -delta, j,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em1m1)

    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, 2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, 2.0d0*delta, j,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2p2)

    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, 2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, -2.0d0*delta, j,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Ep2m2)

    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, 2.0d0*delta, j,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2p2)

    theta = soeo%oldtheta
    call soeo_fd_add_delta_to_theta (theta, -2.0d0*delta, i,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, -2.0d0*delta, j,&
                                 & soeo%Nocc, soeo%Nact, decomp%cfg_unres)
    call soeo_fd_get_energy (decomp, K, theta, soeo%C, soeo%S, Em2m2)
    !-------------------------------------------------------------------

    !Vector-vector-part of the Hessian
    ! hessian = (16G1 - G2)/(48delta^2)
    !  G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    !  G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    !-------------------------------------------------------------------
    G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    hessian(row,col) = (16.0d0*G1 - G2) / (48.0d0*delta*delta)
    !-------------------------------------------------------------------
  enddo
  !---------------------------------------------------------------------  
enddo
!-----------------------------------------------------------------------

call mat_free (theta)
call mat_free (K)

end subroutine soeo_finite_difference
!=======================================================================


!> \brief Subroutine that gets the energy as function of K and theta
!> \author C. Nygaard
!> \date 2010-10-06
!=======================================================================
subroutine soeo_fd_get_energy (decomp, K, theta, Cin, S, E)

implicit none

type(decompItem), intent(in) :: decomp
type(matrix), intent(in)     :: K, theta, Cin, S
real(realk), intent(out)     :: E

integer                      :: Nbast
type(matrix)                 :: C, Dao, Dmo, Fao

Nbast = Cin%nrow

call mat_init (C, Nbast, Nbast)
call mat_init (Dao, Nbast, Nbast)
call mat_init (Dmo, Nbast, Nbast)
call mat_init (Fao, Nbast, Nbast)

C = Cin
call soeo_getnew_C (K, C)
call soeo_getnew_Dmo (theta, Dmo, decomp)
call util_MO_to_AO_2 (S, C, Dmo, Dao, .false.)
call FCK_get_Fock (Dao, Fao, E)

call mat_free (C)
call mat_free (Dao)
call mat_free (Dmo)
call mat_free (Fao)

end subroutine soeo_fd_get_energy
!=======================================================================

!> /brief Adds delta to the i'th variable in the antisymmetric matrix K
!> /author C. Nygaard
!> /date 2010-10-06
!=======================================================================
subroutine soeo_fd_add_delta_to_K (K, delta, i, unres)

implicit none

type(matrix), intent(inout) :: K
real(realk), intent(in)     :: delta
integer, intent(in)         :: i
logical, intent(in)         :: unres

integer                     :: N, Kdim
type(matrix)                :: Kvec
real(realk), allocatable    :: Kelm(:)

if (unres) then
  allocate (Kelm(2))
else
  allocate (Kelm(1))
endif
N = K%nrow
Kdim = (N*(N-1))/2
call mat_init (Kvec, Kdim, 1)

call mat_to_vec ('a', K, Kvec)
if (i <= Kdim) then
  call mat_get_ab_elms (Kvec, i, 1, Kelm)
  Kelm(1) = Kelm(1) + delta
  call mat_create_ab_elms (i, 1, Kelm, Kvec)
elseif (i > Kdim .and. i <= 2*Kdim .and. unres) then
  call mat_get_ab_elms (Kvec, i-Kdim, 1, Kelm)
  Kelm(2) = Kelm(2) + delta
  call mat_create_ab_elms (i-Kdim, 1, Kelm, Kvec)
else
  print *, 'i =', i, 'Kdim =', Kdim, 'unres =', unres
  stop "K does not contain an element in the i'th place"
endif
call mat_vec_to_mat ('a', Kvec, K)

call mat_free (Kvec)
deallocate (Kelm)

end subroutine soeo_fd_add_delta_to_K
!=======================================================================

!> \brief Adds delta to the i'th variable in the vector theta
!> \author C. Nygaard
!> \date 2010-10-06
!=======================================================================
subroutine soeo_fd_add_delta_to_theta (theta, delta, i, Nocc, Nact, unres)

implicit none

type(matrix), intent(inout) :: theta
real(realk), intent(in)     :: delta
integer, intent(in)         :: i, Nocc, Nact
logical, intent(in)         :: unres

real(realk), allocatable    :: thetaelm(:)

if (unres) then
  allocate (thetaelm(2))
else
  allocate (thetaelm(1))
endif

if (i <= Nact) then
  call mat_get_ab_elms (theta, Nocc+i, 1, thetaelm)
  thetaelm(1) = thetaelm(1) + delta
  call mat_create_ab_elms (Nocc+i, 1, thetaelm, theta)
elseif (i > Nact .and. i <= 2*Nact .and. unres) then
  call mat_get_ab_elms (theta, Nocc+i-Nact, 1, thetaelm)
  thetaelm(2) = thetaelm(2) + delta
  call mat_create_ab_elms (Nocc+i-Nact, 1, thetaelm, theta)
else
  print *, 'i =', i, 'Nact =', Nact, 'unres =', unres
  stop "theta does not contain an element in the i'th place"
endif

deallocate (thetaelm)

end subroutine soeo_fd_add_delta_to_theta
!=======================================================================

end module soeo_debug
