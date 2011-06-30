!> @file
!> Contains matrix utilities module.

!> \brief Contains common utility routines for type(matrix).
!> \author S. Host
!> \date March 2010
!>
!> This module depends only on the module matrix_operations. 
!> UNDER NO CIRCUMSTANCES ARE YOU ALLOWED TO INTRODUCE OTHER DEPENDENCIES!
!> If your routine depends on other modules than matrix_operations, it means
!> that it belongs somewhere else than here.
!> 
module matrix_util
use matrix_operations

contains 

  !> \brief Purify the density matrix D by McWeeney purification in AO basis.
  !> \date 2003
  !> \author L. Thogersen
  !> 
  !> This subroutine purifies the density: \n
  !> d_(n+1) = 3d_n*s*d_n - 2 d_n*s*d_n*s*d_n
  !>
  subroutine McWeeney_purify(S,D,failed)
    implicit none
    !> Overlap matrix
    type(Matrix), intent(in) :: S
    !> Input: AO density matrix. Output: Purified AO density matrix
    type(matrix), intent(inout) :: D
    !> True if density matrix could not be purified in 100 iterations
    logical, intent(out) :: failed
    type(Matrix) :: DSD,scr1,scr2
    integer :: i,ndim

    failed = .false.
    ndim = S%nrow
    call mat_init(DSD,ndim,ndim)
    call mat_init(scr1,ndim,ndim)
    call mat_init(scr2,ndim,ndim)

    do i = 1,100
      call mat_mul(D,S,'n','n',1.d0,0.d0,scr1)
      call mat_mul(scr1,D,'n','n',1.d0,0.d0,DSD)
      !if (debug_dd_purify) then
      !  !test for convergence
      !  call mat_add(1.d0,DSD,-1.d0,D,scr1)
      !  WRITE(LUPRI,*) 'Error in purification of it. ',i,' = ',mat_sqnorm2(scr1)
      !endif
      call mat_mul(DSD,S,'n','n',1.d0,0.d0,scr1)
      call mat_mul(scr1,D,'n','n',1.d0,0.d0,scr2)
      scr1 = D                   !scr2 = DSDSD
      call mat_add(3.d0,DSD,-2.d0,scr2,D)
!
! check for convergence
!                               scr1 is the old density 
      call mat_add(1.d0,D,-1.d0,scr1,scr2)
      if (mat_sqnorm2(scr2) < 1d-10) then
        !if (info_dens) then
        !  WRITE(LUPRI,*) 'Purification converged in ',i,' iterations'
        !endif
        exit
      endif
      if (i == 100) then !Stinne moved this inside the purify loop
         failed = .true.
         !stop ' density matrix not purified'
      end if
    enddo
    call mat_free(DSD)
    call mat_free(scr1)
    call mat_free(scr2)
 
 end subroutine McWeeney_purify

 !> \brief Normalize a matrix.
 !> \date 2005
 !> \author S. Host
 subroutine normalize(x)
 implicit none
        !> Input: Matrix to be normalized. Output: Normalized matrix.
        type(matrix), intent(inout)   :: x
        real(realk)                   :: norm

   norm = sqrt(mat_sqnorm2(x))
   call mat_scal(1.0d0/norm,x)
 end subroutine normalize

 !> \brief Calculate D(n+1) = D(X) from Taylor expansion.
 !> \date 2007
 !> \author S. Host
 !>
 !> This routine calculates D(X) in orthonormal AO (OAO) basis 
 !> from X (in OAO basis) and the current D (in OAO basis). 
 !> It is an alternative to the asymmetric BCH expansion, 
 !> and should save some matrix multiplications if ||X|| is large.
 !>
 subroutine oao_density_param(x,D,Dnew)
  implicit none
  !> Antisymmetric X in OAO basis
  type(Matrix), intent(in) :: x
  !> Old density matrix in OAO basis
  type(Matrix), intent(in) :: D
  !> New density matrix D(X) (output)
  type(matrix),intent(inout) :: Dnew 
  type(matrix) :: expx, term, scr 
  real(realk) :: fac, facinv, term_norm, test, xnorm, n
  integer :: i, ndim, j, nmatmul

  nmatmul = 0
  ndim = x%nrow
  call mat_init(expx,ndim,ndim)
  call mat_init(term,ndim,ndim)
  call mat_init(scr,ndim,ndim)

  i = 0
  n = 1.0d0
  xnorm = sqrt(mat_sqnorm2(x)) 
  do 
     i = i + 1
     test = xnorm/n
     if (ABS(test) < 0.2d0) EXIT
     n = n * 2.0d0
  enddo
  !if (info_dens) then
  !   WRITE(LUPRI,*) 'X was divided by', n, 'to obtain an xnorm of', test
  !endif

  Dnew = x   
  call mat_scal(1.0d0/n, Dnew)
  !write(lupri,*) 'x after dividing by', n
  !call mat_print(xn,1,xn%nrow,1,xn%ncol,lupri)
  fac = 1.d0
  term = Dnew
  call mat_identity(expx)
  do j = 1, 30
    fac = fac * j
    facinv = 1.d0/fac
  !write(lupri,*) 'exp(x) after iteration', j
  !call mat_print(expx,1,expx%nrow,1,expx%ncol,lupri)
    call mat_daxpy(facinv,term,expx)
    !
    ! check for convergence
    !
    term_norm = sqrt(mat_sqnorm2(term))
    !if (info_dens) then
    !  WRITE(lupri,*) 'Taylor exp.: it, term_norm',j,term_norm
    !endif
    if (term_norm < 1.0d-8) then
      !write(lupri,*) 'Taylor expansion matrix exponential converged using :', j,' terms!'  
      exit
    endif
    call mat_mul(term,Dnew,'n','n',1.d0,0.d0,scr)
    nmatmul = nmatmul + 1
    term = scr
    if (j == 30) then !Stinne moved this inside i loop
       WRITE(*,*) 'WARNING'
       WRITE(*,*) 'Taylor expansion matrix exponential not converged using 30 terms!'
       WRITE(*,*) 'Last term norm in expansion:', term_norm
       stop
    end if
  enddo

  do j = 1, i-1
     call mat_mul(expx,expx,'n','n',1.d0,0.d0,term)
     nmatmul = nmatmul + 1
     expx = term
  enddo

  !write(lupri,*) 'Total exp(x)'
  !call mat_print(expx,1,expx%nrow,1,expx%ncol,lupri)

  call mat_trans(expx, term) !term = expminusx

  !write(lupri,*) 'Total exp(-x)'
  !call mat_print(expminusx,1,expminusx%nrow,1,expminusx%ncol,lupri)

  !write(lupri,*) 'expx, density_param:'
  !call mat_print(expx,1,D%nrow,1,D%ncol,lupri)
  !write(lupri,*) 'term, density_param:'
  !call mat_print(term,1,D%nrow,1,D%ncol,lupri)
  !write(lupri,*) 'D, density_param:'
  !call mat_print(D,1,D%nrow,1,D%ncol,lupri)

  call mat_mul(term,D,'n','n',1.d0,0.d0,scr)
  call mat_mul(scr,expx,'n','n',1.d0,0.d0,Dnew)
  nmatmul = nmatmul + 2

  !if (info_dens) then
  !   WRITE(LUPRI,*) 'Number of matmuls in Taylor eval. of density:', nmatmul
  !endif 
  call mat_free(expx)
  call mat_free(term)
  call mat_free(scr)
 end subroutine oao_density_param 

 !> \brief Purify the density matrix D by McWeeney purification in OAO basis.
 !> \date 2007
 !> \author S. Host
 !> 
 !> This subroutine purifies the density: \n
 !> d_(n+1) = 3d_n*d_n - 2 d_n*d_n*d_n (in orthonormal AO basis)
 !>
 subroutine oao_purify(D,Dpur)
     implicit none
     !> Density matrix to be purified (in OAO basis)
     type(Matrix), intent(in) :: D
     !> Purified density matrix (in OAO basis) (output)
     type(matrix),intent(inout) :: Dpur 
     type(Matrix) :: DD,DDD,scr1
     integer :: i,ndim

   ndim = D%nrow
   call mat_init(DD,ndim,ndim)
   call mat_init(DDD,ndim,ndim)
   call mat_init(scr1,ndim,ndim)
! 
!  this subroutine purifies the density
!  d_(n+1) = 3d_n*d_n - 2 d_n*d_n*d_n i spin-orbitalbase
! 
   Dpur = D
   do i = 1, 100
      call mat_mul(Dpur,Dpur,'n','n',1.d0,0.d0,DD)
      call mat_mul(DD,Dpur,'n','n',1.d0,0.d0,DDD)
      !if (debug_dd_purify) then
      !  !test for convergence
      !  call mat_add(1.d0,DD,-1.d0,Dpur,scr1)
      !  WRITE(LUPRI,*) 'Error in purification of it. ',i,' = ',sqrt(mat_sqnorm2(scr1))
      !endif
      call mat_add(3.d0,DD,-2.d0,DDD,scr1)
! 
!  check for convergence
!                                 scr1 is the new density 
      call mat_add(1.d0,scr1,-1.d0,Dpur,DDD) !DDD = error matrix
      Dpur = scr1
      if (sqrt(mat_sqnorm2(DDD)) < 1.d-10) then
         !if (info_dens) then
         !   WRITE(LUPRI,*) 'Purification converged in ',i,' iterations'
         !endif
         exit
      endif

      if (i == 100) then
         WRITE(*,*) 'Density matrix not purified in 100 iterations (oao_purify)'
                stop ' Density matrix not purified (oao_purify)'
      end if
   enddo
   call mat_free(DD)
   call mat_free(DDD)
   call mat_free(scr1)
 
 end subroutine oao_purify


  !> \brief Calculates the new density D(n+1) = D(X) in AO basis from BCH expansion.
  !> \author L. Thogersen
  !> \date 2003
  !>
  !>  start D = d  , upon return D = d(x) \
  !>  d(x) = exp(-xs)*d*exp(sx) \n
  !>  = d + [d,x]_s     + 1/2 [[d,x]_s]x]_s   + 1/3*2 [[[d,x]_s]x]_sx]_s + ... \n
  !>  = d + dsx+(dsx)^T + 1/2 ( [d,x]_s sx +([d,x]_s sx)^T) + ... \n
  !>  D,S symmetric and x antisymmetric matrices
  !>
  subroutine density_param_eff(x,S,D)
    implicit none
    !> Antisymmetric matrix X (in AO basis)
    type(Matrix), intent(in) :: x
    !> Overlap matrix
    type(Matrix), intent(in) :: S
    !> Input: Old density matrix (AO basis). Output: New density matrix D(X) (AO basis).
    type(matrix), intent(inout) :: D
    type(matrix) :: Deff, scr, SX
    real(realk) :: fac, facinv, term_norm, testnorm
    integer :: ndim,i, nmatmul

    nmatmul = 0
    ndim = S%nrow
    call mat_init(Deff,ndim,ndim)
    call mat_init(scr,ndim,ndim)
    call mat_init(SX, ndim, ndim)

    call mat_mul(S,X,'n','n',1.d0,0.d0,SX)
    nmatmul = nmatmul + 1
!
! Expansion termwise (orderwise)
!
    !Deff is the latest term in the expansion
    !write(lupri,*) 'X in density param'
    !call mat_print(x,1,ndim,1,ndim,lupri)
    Deff = D
    fac = 1.d0 
    do i = 1,30
      fac = fac * i
      facinv = 1.d0/fac

      call mat_mul(Deff,SX,'n','n',1.d0,0.d0,scr)
      nmatmul = nmatmul + 1
      call mat_trans(scr,Deff)
!      call mat_mul(Deff,S,'n','n',1.d0,0.d0,scr)
!      call mat_mul(scr,x,'n','n',1.d0,0.d0,Deff)
!      call mat_trans(Deff,scr)
      call mat_daxpy(1.d0,scr,Deff)
      !
      ! Deff is now the new term
      ! this term is added to D
      !
      call mat_daxpy(facinv,Deff,D)
      !
      ! check for convergence
      !
      term_norm = sqrt(mat_sqnorm2(Deff))*facinv
      !if (info_dens) then
      !  WRITE(LUPRI,*) 'DENSITY  it., term_norm',i,term_norm
      !endif
      if (term_norm < 1.0d-8) then
        exit
      endif
      if (i == 30) then !Stinne moved this inside i loop
         WRITE(*,*) 'WARNING'
         WRITE(*,*) 'exponential expansion of density matrix not' &
         &       //' converged using 30 terms'
         WRITE(*,*) 'Last term norm in expansion:', term_norm
         stop
      end if
    enddo
    !if (info_dens) then
    !   WRITE(LUPRI,*) 'Number of matmuls in exp. param. of density:', nmatmul
    !endif
    call mat_free(Deff)
    call mat_free(scr)
    call mat_free(SX)
  end subroutine density_param_eff

   !> \brief Evaluate norm of the electronic gradient in AO basis according to the usual formula 4*(FDS - SDF).
   !> \author L. Thogersen
   !> \date 2003
   !>
   !> NOTE: F CAN BE BOTH SYMM AND ANTISYMM !
   !> 
   SUBROUTINE get_AO_gradient(F,D,S, grad)
      IMPLICIT NONE
      !> Fock/KS matrix
      TYPE(Matrix),intent(in) :: F
      !> Density matrix
      TYPE(Matrix),intent(in) :: D
      !> Overlap matrix
      TYPE(Matrix),intent(in) :: S
      !> AO gradient 4*(FDS - SDF) (output)
      TYPE(Matrix),intent(inout) :: grad
      TYPE(Matrix) :: FDS, SDF, SD
      
      CALL mat_init(FDS, F%nrow, F%ncol)
      CALL mat_init(SDF, F%nrow, F%ncol)
      call mat_init(SD , F%nrow, F%ncol)
#ifdef THE_MORE_EFFICIENT_CODE_BUT_NOT_WORKING_WITH_PURIFICATION_YET
      CALL mat_mul(S,   D,'n','n', 1D0, 0D0, SD)
      CALL mat_mul(SD, F,'n','n', 4D0, 0D0, SDF)
      call mat_mul(F,SD,'n','t',4.0d0,0.d0, FDS)
#else
      CALL mat_mul(S,   D,'n','n', 1D0, 0D0, SD)
      CALL mat_mul(SD, F,'n','n', 4D0, 0D0, SDF)
      call mat_mul(F,SD,'n','t',4.0d0,0.d0, FDS)
#endif
      CALL mat_add(1D0, SDF, -1D0, FDS, grad) ! res = SDF - tmp
      CALL mat_free(FDS)
      CALL mat_free(SDF)
      call mat_free(SD)
   END SUBROUTINE get_AO_gradient

   !> \brief Get gradient in orthonormal AO basis.
   !> \author S. Host
   !> \date 2005
   subroutine get_OAO_gradient(F, D, grad)
   implicit none
      !> Fock/KS matrix in orthonormal AO basis
      type(matrix), intent(in) :: F
      !> Density matrix in orthonormal AO basis
      type(matrix), intent(in) :: D
      !> Gradient in orthonormal AO basis (output)
      type(matrix)    :: grad 
      type(matrix)    :: wrk1
      integer         :: ndim

   ndim = F%nrow
   !Grad = -4*(D*F - F*D)

   call mat_init(wrk1,ndim,ndim)

   call MAT_MUL(D,F,'n','n',1.0d0,0.0d0,wrk1)
   grad = wrk1
   call mat_trans(grad,wrk1)
   call mat_daxpy(-1.0d0,wrk1,grad)

   call mat_scal(-4.0d0,grad)

   call mat_free(wrk1)
   end subroutine get_OAO_gradient

   !> \brief Get the "S-norm" of the matrix A, Tr(A*S*A*S)
   !> \author L. Thogersen
   !> \date 2003
   !>
   !> Output: TrASAS
   !> 
   real(realk) function util_Snorm(A,S)
     implicit none
     !> Overlap matrix
     type(Matrix), intent(in) :: S
     !> Matrix of which we want the S-norm
     type(Matrix), intent(in) :: A
     type(matrix)             :: AS, SAS ! pointer

     call mat_init(AS,A%nrow,A%ncol)
     call mat_init(SAS,A%nrow,A%ncol)
     call mat_mul(A,S,'n','n',1.d0,0.d0,AS)
     call mat_mul(S,AS,'n','n',1.d0,0.d0,SAS)
     util_Snorm = mat_dotproduct(A,SAS)
     call mat_free(AS)
     call mat_free(SAS)

   end function util_Snorm

   !> \brief Returns the symmetric part of a matrix, [M]^S = 1/2(M+MT)
   !> \author L. Thogersen
   !> \date 2003
    subroutine util_get_symm_part(A)
      implicit none
      !> Input: Matrix from which we want symmetric part. Output: Symmetric part of input matrix.
      type(Matrix),intent(inout) :: A
      type(Matrix) :: Ax,AT

      if (A%nrow /= A%ncol) then
        !WRITE(LUPRI,*) 'Trying to symmetrice matrix where nrow /= ncol'
        STOP 'Trying to symmetrice matrix where nrow /= ncol'
      endif
      call mat_init(Ax,A%nrow,A%ncol)
      call mat_init(AT,A%nrow,A%ncol)
      Ax = A
      call mat_trans(A,AT)
      call mat_add(0.5d0,Ax,0.5d0,AT,A)
      call mat_free(Ax)
      call mat_free(AT)
    end subroutine util_get_symm_part

    !> \brief Returns the antisymmetric part of a matrix, [M]^A = 1/2(M-MT)
    subroutine util_get_antisymm_part(A,B)
       implicit none
       !> Matrix from which we want antisymmetric part
       type(Matrix),intent(in)    :: A
       !> Antisymmetric part of A, B = 0.5*(A-AT)
       type(Matrix),intent(inout) :: B
       type(Matrix) :: AT
       
       if (A%nrow /= A%ncol) then
          !WRITE(LUPRI,*) 'Trying to antisymmetrize matrix where nrow /= ncol'
          STOP 'Trying to antisymmetrize matrix where nrow /= ncol'
       endif
       call mat_init(AT,A%nrow,A%ncol)
       call mat_trans(A,AT)
       call mat_add(0.5d0,A,-0.5d0,AT,B)
       call mat_free(AT)
 
    end subroutine util_get_antisymm_part

    !> \brief Convert type(matrix) from MO to AO basis
    !> \author C. Nygaard
    !> \date 2010
    !> \param S The overlap matrix
    !> \param C The MO coefficient matrix
    !> \param X_mo The input matrix in MO basis
    !> \param X_ao The output matrix in AO basis
    !> \param cov logical, .true. if X is covariant
    subroutine util_MO_to_AO_2(S,C,X_mo,X_ao,cov)
      ! if covariant
      !Transform from MO to covariant AO basis
      !  X_ao =  S C X_mo C^T S
      ! else contravariant
      !Transform from MO to contravariant AO basis
      !X_ao = C X_mo C^T
      implicit none
      logical :: cov
      type(Matrix), intent(in) :: S,X_mo,C
      type(Matrix), intent(inout) :: X_ao  !output
      type(Matrix) :: X, SC
  
      call mat_init(X,X_mo%nrow,X_mo%ncol)
      call mat_init(SC,S%nrow,S%ncol)

      if (cov) then
        call mat_mul(S,C,'n','n',1.d0,0.d0,SC)
      else
        SC = C
      endif
      call mat_mul(X_MO,SC,'n','t',1.d0,0.d0,X)
      call mat_mul(SC,X,'n','n',1.d0,0.d0,X_ao)

      call mat_free(X)
      call mat_free(SC)
      
    end subroutine util_MO_to_AO_2

    !> \brief Convert type(matrix) from AO to MO basis
    !> \author C. Nygaard
    !> \date 2010
    !> \param S The overlap matrix
    !> \param C The MO coefficient matrix
    !> \param X_ao The input matrix in AO basis
    !> \param X_mo The output matrix in MO basis
    !> \param cov logical, .true. if X is covariant 
    subroutine util_AO_to_MO_2(S,C,X_ao,X_mo,cov)
      ! if covariant
      !Transform from the AO covariant basis to the MO basis
      !X_mo = C^T X_ao C
      ! else contravariant
      !Transform from the AO contravariant basis to the MO basis
      !X_mo = C^T S X_ao S C
      implicit none
      logical :: cov
      type(Matrix), intent(in) :: S,X_ao,C
      type(Matrix), intent(inout) :: X_mo  !output
      type(Matrix) :: X,SC
      integer :: ndim

      call mat_init(X,S%nrow,S%ncol)
      call mat_init(SC,S%nrow,S%ncol)
      ndim = S%nrow
      if (cov) then
        SC = C
      else
        call mat_mul(S,C,'n','n',1.d0,0.d0,SC)
      endif
      call mat_mul(X_ao,SC,'n','n',1.d0,0.d0,X)
      call mat_mul(SC,X,'t','n',1.d0,0.d0,X_mo)

      call mat_free(X)
      call mat_free(SC)

    end subroutine util_AO_to_MO_2

   !> \brief Interface to DSYEVX for diagonalization of real symmetric matrix A*x = mu*X
   !> \author S. Host
   !> \date 2005
   subroutine util_diag(lupri,A,print_eivecs,scale_eig,low_eig_string)
   implicit none
        !> Logical unit number of output file
        integer,intent(in)       :: lupri
        !> Matrix to be diagonalized
        type(matrix), intent(in) :: A
        !> True if if eigenvectors should be printed
        logical, intent(in)      :: print_eivecs
        !> Eigenvalues are scaled by this factor before printing 
        real(realk), intent(in)  :: scale_eig 
        !> String used when printing lowest eigenvalue (for grep)
        character(*), intent(in) :: low_eig_string 
        real(realk), allocatable :: Afull(:,:), eigenval(:), eigenvec(:,:)
        real(realk), allocatable :: temp(:)
        integer, allocatable     :: itemp(:), IFAIL(:)
        real(realk)              :: VL, VU, low_eig
        integer                  :: IL, IU, neig, Ltemp, INFO, m, ndim
        logical                  :: low_eig_found 

   low_eig_found = .false.
   ndim = A%nrow
   Ltemp = 8*ndim 
   allocate(Afull(ndim,ndim),eigenval(ndim))
   allocate(eigenvec(ndim,ndim),temp(Ltemp),Itemp(5*ndim),IFAIL(ndim))

   call mat_to_full(A,1.0d0,Afull)

   !write (lupri,*) 'Afull (hessian):'
   !call OUTPUT(Afull, 1, ndim, 1, ndim, ndim, ndim, 1, lupri)

   call DSYEVX('V', 'A', 'U', ndim, Afull, ndim, VL, VU, IL, IU, &
     &  0.0d0, neig, eigenval, eigenvec, ndim, temp, Ltemp, Itemp, &
     &  IFAIL, INFO )

   if (info /= 0) STOP 'Problem in DSYEVX (subroutine util_diag)'

   write(lupri,*) "Dim of A and eigenvalues of A, # found:", ndim, neig
   write(lupri,*) "Eigenvalues are scaled by", scale_eig, "before printing!"
   do m = 1, neig
      write(lupri,'(i6,F15.7)') m, eigenval(m)*scale_eig
      if (.not. low_eig_found) then
         if (abs(eigenval(m)) > 1.0d-6) then
            low_eig = eigenval(m)*scale_eig
            low_eig_found = .true.
         endif
      endif
   enddo

   write(lupri,*) low_eig_string, low_eig

   if (print_eivecs) then
      write (lupri,*) 'Eigenvectors of A:'
      call OUTPUT(eigenvec, 1, ndim, 1, ndim, ndim, ndim, 1, lupri)
   endif

   deallocate(Afull)
   deallocate(eigenval)
   deallocate(eigenvec)
   deallocate(temp)
   deallocate(Itemp)
   deallocate(IFAIL)
   end subroutine util_diag

   !> \brief Get density matrix from canonical molecular orbitals.
   !> \author L. Thogersen
   !> \date 2003
   subroutine density_from_orbs(unres,nocc,nocca,noccb,CMO,D)
     implicit none
     !> True if unrestricted calculation
     logical,intent(in)       :: unres
     !> Number of occupied orbitals
     integer,intent(in)       :: nocc
     !> Number of occupied alpha orbitals (only referenced if unres=true)
     integer,intent(in)       :: nocca
     !> Number of occupied beta orbitals (only referenced if unres=true)
     integer,intent(in)       :: noccb
     !> Canonical molecular orbitals
     type(Matrix), intent(in) :: CMO
     !> Density matrix constructed from CMOs (output)
     type(Matrix),intent(inout) :: D 
     real(realk), allocatable :: orb(:),den(:),orb_ab(:,:),den_ab(:,:)
     type(Matrix) :: cmo_occ!,cmo_occ2
     integer :: ndim
   
     ndim=CMO%nrow
     if(unres) then
        allocate(orb_ab(ndim*ndim,2),den_ab(ndim*ndim,2))
        call mat_to_full2(CMO,1d0,orb_ab)
        CALL DGEMM('N','T',ndim,ndim,nocca,1d0,orb_ab(:,1),ndim,orb_ab(:,1),ndim, &
                  & 0d0,den_ab(:,1),ndim)
        CALL DGEMM('N','T',ndim,ndim,noccb,1d0,orb_ab(:,2),ndim,orb_ab(:,2),ndim, &
                  & 0d0,den_ab(:,2),ndim)
        call mat_set_from_full2(den_ab,1d0,D)
        deallocate(orb_ab,den_ab)
     else
#if 0
        call mat_init(cmo_occ,ndim,nocc)
        !call mat_init(cmo_occ2,ndim,nocc)
        call mat_section(CMO,1,ndim,1,nocc,cmo_occ)
        !cmo_occ2 = cmo_occ
        !call mat_mul(cmo_occ,cmo_occ2,'n','t',1.d0,0.d0,D)
        call mat_mul(cmo_occ,cmo_occ,'n','t',1.d0,0.d0,D)
        call mat_free(cmo_occ)
        !call mat_free(cmo_occ2)
#endif
        allocate(den(ndim*ndim),orb(ndim*ndim))
        call mat_to_full(CMO,1d0,orb)
        CALL DGEMM('N','T',ndim,ndim,nocc,1d0,orb,ndim,orb,ndim, &
                  & 0d0,den,ndim)
        call mat_set_from_full(den,1d0,D)
        deallocate(orb,den)
     end if
     return
   end subroutine density_from_orbs

!> \brief Calculate the commutator [A,B] = alpha * (AB - BA)
!> \author L. Thogersen
!> \date 2003
  subroutine commutator(alpha,A,B,transa,transb,Commut)
    implicit none
    !> Multiply commutator by this factor
    real(realk), intent(in) :: alpha
    !> First matrix of commutator
    type(Matrix), intent(in) :: A
    !> Second matrix of commutator
    type(Matrix), intent(in) :: B
    !> T/t if A should be transposed when constructing commutator, otherwise N/n
    character, intent(in) :: transa
    !> T/t if B should be transposed when constructing commutator, otherwise N/n
    character, intent(in) :: transb
    !> The commutator [A,B] = alpha * (AB - BA) (output)
    type(Matrix), intent(inout) :: Commut
    type(matrix) :: Prod1,Prod2
    integer :: Ndim

    Ndim = A%nrow

    call mat_init(Prod1,Ndim,Ndim)
    call mat_init(Prod2,Ndim,Ndim)

    call mat_mul(A,B,transa,transb,alpha,0.d0,Prod1)
    call mat_mul(B,A,transb,transa,alpha,0.d0,Prod2)
    call mat_add(1.d0, Prod1, -1.d0, Prod2, Commut)

    call mat_free(Prod1)
    call mat_free(Prod2)

  end subroutine commutator

  !> \brief Calculate the generalized commutator [A,B]_C = A*C*B - B*C*A
  !> \author S. Coriani
  !> \date June 2005
  subroutine ABCcommutator(ndim,A,B,C,Commutator)
    implicit none
    !> Dimension of matrices - OBSOLETE, SHOULD BE REMOVED!
    integer, intent(in) :: ndim
    !> First matrix of commutator
    type(Matrix), intent(in) :: A
    !> Second matrix of commutator
    type(Matrix), intent(in) :: B
    !> Third matrix of commutator
    type(Matrix), intent(in) :: C
    !> The commutator [A,B]_C = A*C*B - B*C*A (output)
    type(Matrix), intent(inout) :: Commutator
    type(matrix) :: Prod,Prod1
   
    call mat_init(Prod,ndim,ndim)
    call mat_init(Prod1,ndim,ndim)

    call mat_mul(A,C,'n','n',1.d0,0.d0,Prod)    
    call mat_mul(Prod,B,'n','n',1.d0,0.d0,Commutator)  
    call mat_mul(C,A,'n','n',1.d0,0.d0,Prod)    
    call mat_mul(B,Prod,'n','n',1.d0,0.d0,Prod1) 
    call mat_daxpy(-1.d0,Prod1,Commutator) 
!

    call mat_free(Prod)
    call mat_free(Prod1)
!  
  end subroutine ABCcommutator

 !> \brief Calculate the generalized commutator [A,B]_C = A*C*B - B*C*A
 !> \author Thomas Kjaergaard
 !> \date February 2008
 subroutine ABCcommutator2(ndim,alpha,A,B,C,transa,transb,transc,Commutator)
    implicit none
    !> Dimension of matrices - OBSOLETE, SHOULD BE REMOVED!
    integer, intent(in) :: ndim
    !> Multiply commutator by this factor
    real(realk), intent(in) :: alpha
    !> First matrix of commutator
    type(Matrix), intent(in) :: A
    !> Second matrix of commutator
    type(Matrix), intent(in) :: B
    !> Third matrix of commutator
    type(Matrix), intent(in) :: C
    !> T/t if A should be transposed when constructing commutator, otherwise N/n
    character, intent(in) :: transa
    !> T/t if B should be transposed when constructing commutator, otherwise N/n
    character, intent(in) :: transb
    !> T/t if C should be transposed when constructing commutator, otherwise N/n
    character, intent(in) :: transc
    !> The commutator [A,B]_C = A*C*B - B*C*A 
    type(Matrix), intent(inout) :: Commutator
    type(matrix)                :: Prod,Prod1
   
    call mat_init(Prod,ndim,ndim)
    call mat_init(Prod1,ndim,ndim)

    call mat_mul(A,C,transa,transc,1.d0,0.d0,Prod)    
    call mat_mul(Prod,B,'n',transb,alpha,0.d0,Commutator)  
    call mat_mul(C,A,transc,transa,1.d0,0.d0,Prod)    
    call mat_mul(B,Prod,transb,'n',alpha,0.d0,Prod1) 
    call mat_daxpy(-1.d0,Prod1,Commutator) 

    call mat_free(Prod)
    call mat_free(Prod1)
  
  end subroutine ABCcommutator2


  !> \brief Makes the exponential of a matrix.
  !> \author B. Jansik
  !> The algorithm scales input X by inverse
  !> frobenius norm relying on inequality nrm_f > nrm_2
  !> convergence is tested against last term of the 
  !> Taylor expansion relying on nrm_f(A+B)< nrm_f(A)+nrm_f(B)
  !> inequality. Only two temporary matrix allocations needed.
  !> Near optimal with respect to number of matrix multiply
  subroutine matrix_exponential(X,expX,thr)
  implicit none
  type(Matrix), intent(inout) :: expX
  type(Matrix), intent(in)    :: X
  real(realk),  intent(in)    :: thr

  type(Matrix) :: Xn, tmp
  real(realk) :: fac, nrmX, scal
  integer     :: i, pwr

  pwr = 3
  fac = 1d0
  scal= 1d0

  nrmX = sqrt(mat_sqnorm2(X))

  if (nrmX.gt.1d0) then
    pwr = int(ceiling(log(nrmX)/log(2d0))) + 3
  endif

  scal= 2d0**(real(-pwr))
  

  call mat_init(Xn,X%ncol,X%nrow)
  call mat_init(tmp,X%ncol,X%nrow)

  call mat_identity(Xn)
  call mat_identity(expX)

  do i=1,30
    fac = fac*(1d0/real(i))

    call mat_mul(Xn,X,'n','n',scal,0d0,tmp)
    call mat_assign(Xn,tmp)

    call mat_daxpy(fac,Xn,expX)
    
    nrmX = fac*sqrt(mat_sqnorm2(Xn))
    if (nrmX.le.thr) exit

  enddo

  call mat_free(Xn)

  do i=1,pwr
    call mat_mul(expX,expX,'n','n',1d0,0d0,tmp)
    call mat_assign(expX,tmp)
  enddo


  call mat_free(tmp)
    
  end subroutine matrix_exponential

  !> \brief Check that AO density matrix is idempotent, i.e. D=DSD
  !> \author L. Thogersen
  !> \date 2003
  subroutine check_idempotency(D,S)
    type(Matrix), intent(in)  :: D,S
    type(Matrix) :: tmp,DSD 

    call mat_init(tmp,S%nrow,S%ncol)
    call mat_init(DSD,S%nrow,S%ncol)

    call mat_mul(S,D,'n','n',1.d0,0.d0,tmp)
    call mat_mul(D,tmp,'n','n',1.d0,0.d0,DSD)
    call mat_daxpy(-1D0,D,DSD)
    norm = mat_dotproduct(DSD,DSD)
    WRITE(LUPRI,*) 'norm of DSD-D:',norm
    if (norm > 1.0d-10) then
      WRITE(LUPRI,*) 'WARNING - DSD not equal to D'
    endif
    call mat_free(tmp)
    call mat_free(DSD)
  end subroutine check_idempotency

  
  !> \brief Returns the symmetry of a matrix
  !> \author S. Reine
  !> \date 2010-04-19
  !> \param A The matrix
  !> \param thresh Elementwise threshold
  !> \param mat_get_isym The symmetry of the matrix A
  !> The mat_get_isym can return the following values:
  !>     1    A is symmetric
  !>     2    A is anti-symmetric
  !>     3    no symmetry
  !>     4    zero matrix
  function mat_get_isym(A,thresh)
  implicit none
  type(Matrix), intent(in) :: A
  real(realk),intent(in)   :: thresh
  integer                  :: mat_get_isym
! 
  Integer :: isym
  type(matrix) :: B,C
  real(realk)  :: norm
 
  norm = sqrt(mat_sqnorm2(A)/min(A%nrow,A%ncol))
! If the square-norm is below thresh it is a zero matrix
  IF (norm.LT.thresh) THEN
    isym = 4
  ELSE IF (A%nrow.EQ.A%ncol) THEN
    call mat_init(B,A%nrow,A%ncol)
    call util_get_antisymm_part(A,B)
    norm = sqrt(mat_sqnorm2(B)/A%nrow)
!   If the anti-symmetrix part of the matrix is zero (to within a rms-norm threshold) 
!   the matrix must be symmetric
    IF (norm.LT.thresh) THEN
      isym = 1
    ELSE
      call mat_init(C,A%nrow,A%ncol)
      call mat_add(1.d0,A,-1.d0,B,C)
      norm = sqrt(mat_sqnorm2(C)/A%nrow)
!     If A - anti-sym(A) = 0
      IF (norm.LT.thresh) THEN
        isym = 2
      ELSE
        isym = 3
      ENDIF
      call mat_free(C)
    ENDIF
    call mat_free(B)
  ELSE
   isym = 3
  ENDIF

  mat_get_isym = isym
    
  end function mat_get_isym

  function mat_same(A,B,thresh)
  implicit none
  type(Matrix), intent(in) :: A,B
  real(realk),intent(in)   :: thresh
  logical                  :: mat_same
!
  type(matrix) :: C
  real(realk)  :: norm
  logical      :: same

  IF ((A%nrow.EQ.B%nrow).AND.(A%ncol.EQ.B%ncol)) THEN
    call mat_init(C,A%nrow,A%ncol)
    call mat_add(1.d0,A,-1.d0,B,C)
    norm = sqrt(mat_sqnorm2(C)/min(A%nrow,A%ncol))
    IF (norm.LT.thresh) THEN
      same = .TRUE.
    ELSE
      same = .FALSE.
    ENDIF
    call mat_free(C)
  ELSE
    same = .FALSE.
  ENDIF

  mat_same = same

  end function mat_same

  !> \brief Dump F, D matrices to disk for detailed investigation of sparsity
  !> \author S. Host
  !> \date 03-05-2010
  subroutine dumpmats(iteration,D,F,optlevel)
  use files
  implicit none
     !> SCF iteration number
     integer, intent(in) :: iteration
     !> OAO density matrix from this iteration
     type(matrix), intent(in) :: D 
     !> OAO Fock/KS matrix from this iteration
     type(matrix), intent(in) :: F
     !> 2 = level 2 (valence optimization), 3 = level 3 (full optimization)
     integer, intent(in)      :: optlevel
     integer                  :: lunD, lunF, i, j
     character(len=80)        :: cit
     character(len=11)        :: filenameD, filenameF !iteration allowed to be up to 3 digits

    if (iteration < 0 .or. iteration > 999) then
       write(*,*) 'iteration =', iteration
       call lsquit('iteration out of bounds in dumpmats',-1)
    endif

    !Convert integer to character via internal file:
    write(cit,*) iteration

    !Left-justify the string:
    cit = adjustl(cit)

    !Construct filenames:
    if (optlevel == 2) then
       filenameD = 'L2_D' // trim(cit) // '.mat'
       filenameF = 'L2_F' // trim(cit) // '.mat'
    else if (optlevel == 3) then
       filenameD = 'L3_D' // trim(cit) // '.mat'
       filenameF = 'L3_F' // trim(cit) // '.mat'
    else
       call lsquit('Unsupported optlevel in dumpmats',-1)
    endif

    lunD = -1
    lunF = -1
    CALL LSOPEN(lunD,filenameD,'NEW','UNFORMATTED')
    CALL LSOPEN(lunF,filenameF,'NEW','UNFORMATTED')

    call mat_write_to_disk(lunD,D)
    call mat_write_to_disk(lunF,F)

    CALL LSCLOSE(lunD,'KEEP')
    CALL LSCLOSE(lunF,'KEEP')

  end subroutine dumpmats

end module matrix_util
