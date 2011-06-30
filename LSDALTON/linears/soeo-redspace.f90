module soeo_redspace

use files
use soeo_matop
use soeo_transform
use Matrix_Operations
use matrix_util
use dal_interface
use decompMod

implicit none

!Contains:
!  soeo_solver
!  soeo_solve_linear_system
!  soeo_solver_getX
!  soeo_solver_getres
!  soeo_solve_eigeneqs
!  soeo_find_mineval
!  soeo_binary_search_init
!  soeo_binary_search

contains

!> \brief Solves the SOEO-Newton equations in a reduced space
!> \author C. Nygaard
!> \date 2010
!> \param Fmo The Fock matrix in MO basis
!> \param Dmo The density matrix in MO basis
!> \param C The MO-orbital coefficient matrix
!> \param S The overlap matrix
!> \param Dao The density matrix in AO basis
!> \param nfirst The first derivatives of Dmo wrt the occupation angles
!> \param nsecond The second derivatives of Dmo wrt the occupation angles
!> \param orbE The MO orbital energies
!> \param maxiter Maximum number of reduced space iterations
!> \param thresh Convergence threshold for reduced space iterations
!> \param Nocc Number of occupied MO orbitals
!> \param oldtheta The starting point for the occupation angles
!> \param K The MO orbital rotation matrix
!> \param deltatheta The change in the occupation angles
!=======================================================================
subroutine soeo_solver (soeo, orbE, K, deltatheta, mu, decomp)

implicit none

!Input/output
type(decompItem), intent(in) :: decomp
type(soeoItem), intent(in)   :: soeo
real(realk), intent(in)      :: orbE(:)
type(matrix), intent(inout)  :: K, deltatheta
real(realk), intent(out)     :: mu

!Other:
integer                      :: iter, lub_m, lub_v, lusigma_m, lusigma_v
integer                      :: i, j
real(realk)                  :: val1, val2, tmp
logical                      :: err, antisymmetric
!Needed matrix-parts:
type(matrix)                 :: grad_m, b1_m, sigma1_m, bi_m,&
                               & biter_m, sigmai_m, sigmaiter_m
!Needed vector-parts:
type(matrix)                 :: grad_v, b1_v, sigma1_v, bi_v,&
                               & biter_v, sigmai_v, sigmaiter_v
!The reduced space:
real(realk)                  :: Ared(soeo%micromaxiter,soeo%micromaxiter),&
                               & gradred(soeo%micromaxiter),&
                               & xred(soeo%micromaxiter)
!Trust-region:
real(realk)                  :: mine

!Solution and residual:
type(matrix)                 :: X_m, X_v, res_m, res_v
real(realk)                  :: resnorm

!Initializations:
call mat_init (grad_m     , soeo%Nbast, soeo%Nbast)
call mat_init (grad_v     , soeo%Nact , 1    )
call mat_init (b1_m       , soeo%Nbast, soeo%Nbast)
call mat_init (b1_v       , soeo%Nact , 1    )
call mat_init (bi_m       , soeo%Nbast, soeo%Nbast)
call mat_init (bi_v       , soeo%Nact , 1    )
call mat_init (biter_m    , soeo%Nbast, soeo%Nbast)
call mat_init (biter_v    , soeo%Nact , 1    )
call mat_init (sigma1_m   , soeo%Nbast, soeo%Nbast)
call mat_init (sigma1_v   , soeo%Nact , 1    )
call mat_init (sigmai_m   , soeo%Nbast, soeo%Nbast)
call mat_init (sigmai_v   , soeo%Nact , 1    )
call mat_init (sigmaiter_m, soeo%Nbast, soeo%Nbast)
call mat_init (sigmaiter_v, soeo%Nact , 1    )
call mat_init (X_m        , soeo%Nbast, soeo%Nbast)
call mat_init (X_v        , soeo%Nact , 1    )
call mat_init (res_m      , soeo%Nbast, soeo%Nbast)
call mat_init (res_v      , soeo%Nact , 1    )

lub_m = -1 ; lub_v = -1 ; lusigma_m = -1 ; lusigma_v = -1
call lsopen (lub_m, "bmats", "unknown", "UNFORMATTED")
call lsopen (lub_v, "bvecs", "unknown", "UNFORMATTED")
call lsopen (lusigma_m, "sigmamats", "unknown", "UNFORMATTED")
call lsopen (lusigma_v, "sigmavecs", "unknown", "UNFORMATTED")

!Get the gradient (both the matrix-part and the vector-part)
!-----------------------------------------------------------------------
call soeo_get_gradient (soeo, grad_m, grad_v, decomp)
!-----------------------------------------------------------------------

!Setup the first trial-vector b1 = norm(pre(norm(grad)))
!-----------------------------------------------------------------------
b1_m = grad_m ; b1_v = grad_v

!call mat_mo_precond (CFG_NOCC, 0.0d0, orbE, b1_m)

call soeo_normalize (b1_m, b1_v)

rewind (lub_m) ; rewind (lub_v)
call mat_write_to_disk (lub_m, b1_m)
call mat_write_to_disk (lub_v, b1_v)
!-----------------------------------------------------------------------

!Building the reduced space (first entries):
!-----------------------------------------------------------------------
iter = 1
Ared = 0.0d0 ; gradred = 0.0d0
call soeo_linear_transform (soeo, b1_m, b1_v, sigma1_m, sigma1_v, decomp)

write (decomp%lupri, *) 'grad ='
call mat_print (grad_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
call mat_print (grad_v, 1, soeo%Nact, 1, 1, decomp%lupri)
write (decomp%lupri, *) 'b1 ='
call mat_print (b1_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
call mat_print (b1_v, 1, soeo%Nact, 1, 1, decomp%lupri)
write (decomp%lupri, *) 'sigma1 ='
call mat_print (sigma1_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
call mat_print (sigma1_v, 1, soeo%Nact, 1, 1, decomp%lupri)


rewind (lusigma_m) ; rewind (lusigma_v)
call mat_write_to_disk (lusigma_m, sigma1_m)
call mat_write_to_disk (lusigma_v, sigma1_v)

Ared(1,1) = soeo_dotproduct (b1_m, b1_v, sigma1_m, sigma1_v)
gradred(1) = -soeo_dotproduct (b1_m, b1_v, grad_m, grad_v)
!-----------------------------------------------------------------------

write (decomp%lupri, *) 'Ared ='
call output (Ared, 1, iter, 1, iter, 200, 200, 1, decomp%lupri)
write (decomp%lupri, *) 'gradred ='
call output (gradred, 1, iter, 1, 1, 200, 1, 1, decomp%lupri)

mu = 0.0d0
do

  !Find mu (so that step is inside trust region)
  !---------------------------------------------------------------------
!  if (iter /= 1) then
    call soeo_find_mineval (decomp%lupri, iter, Ared, mine)
    call soeo_binary_search_init (iter, soeo%Nbast, soeo%Nact, lub_m, lub_v,&
                                & Ared, gradred, mine, soeo%trust, mu, decomp)

write (decomp%lupri, *) 'iter, mineval, mu =', iter, mine, mu
write (decomp%lupri, *) 'Ared ='
call output (Ared, 1, iter, 1, iter, soeo%micromaxiter, soeo%micromaxiter, 1, decomp%lupri)

!  endif
  !---------------------------------------------------------------------

  !Solving the reduced space:
  !---------------------------------------------------------------------
  call soeo_solve_linear_system (iter, Ared, mu, gradred, xred, decomp)

  call soeo_solver_getX (iter, lub_m, lub_v, xred, X_m, X_v)
  call soeo_solver_getres (iter, lusigma_m, lusigma_v, xred, grad_m,&
                         & grad_v, res_m, res_v)
  call soeo_daxpy (-mu, X_m, X_v, res_m, res_v)
  !---------------------------------------------------------------------

  !Test for convergence:
  !---------------------------------------------------------------------
  resnorm = soeo_norm (res_m, res_v)

  if (iter == 1) then
    write (decomp%lupri, *) "Test for convergence in soeo_solver:"
  endif
  write (decomp%lupri, *) "resnorm =", resnorm

write (decomp%lupri, *) "norm(x) =", soeo_norm (X_m, X_v)

  if (resnorm < soeo%microthresh) then
    write (decomp%lupri, '("Reduced space calculation converged in ", &
                   & i4, " iterations")') iter
    exit
  elseif (iter >= soeo%micromaxiter) then
    write (decomp%lupri, '("Reduced space calculation NOT converged &
                   &after ", i4, " iterations!!!")') iter
    call lsquit ("Error in soeo reduced space",decomp%lupri)
  endif
  !---------------------------------------------------------------------

  iter = iter + 1

  !The new trial_vectors: biter = norm(pre(norm(res)))
  !---------------------------------------------------------------------
  biter_m = res_m ; biter_v = res_v
  rewind (lub_m) ; rewind (lub_v)

  call soeo_normalize (biter_m, biter_v)
!  call mat_mo_precond (CFG_NOCC, 0.0d0, orbE, biter_m)
  do i=1,iter-1
    call mat_read_from_disk (lub_m, bi_m)
    call mat_read_from_disk (lub_v, bi_v)
    call soeo_orthonormalize (bi_m, bi_v, biter_m, biter_v, err)

!print *, 'dot(bi, biter) =', soeo_dotproduct (bi_m, bi_v, biter_m, biter_v)

    if (err) then
      write (decomp%lupri, *) "Linear dependencies in soeo_solver!!!"
      write (decomp%lupri, *) "norm(biter) =", soeo_norm (biter_m, biter_v)
      write (decomp%lupri, '("biter * b", i1, " =")') i
      write (decomp%lupri, *) soeo_dotproduct (biter_m, biter_v, bi_m, bi_v)
      call lsquit ("Linear dependencies in soeo_solver!!!",decomp%lupri)
    endif
  enddo

  call mat_write_to_disk (lub_m, biter_m)
  call mat_write_to_disk (lub_v, biter_v)
  !---------------------------------------------------------------------

  !Expanding the reduced space:
  !---------------------------------------------------------------------
  call soeo_linear_transform (soeo, biter_m, biter_v, sigmaiter_m, sigmaiter_v, decomp)

  call mat_write_to_disk (lusigma_m, sigmaiter_m)
  call mat_write_to_disk (lusigma_v, sigmaiter_v)

  rewind (lub_m) ; rewind (lub_v)
  rewind (lusigma_m) ; rewind (lusigma_v)

write (decomp%lupri, *) 'iter =', iter
!write (decomp%lupri, *) 'sigmaiter_m ='
!call mat_print (sigmaiter_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
!write (decomp%lupri, *) 'sigmaiter_v ='
!call mat_print (sigmaiter_v, 1, soeo%Nact, 1, 1, decomp%lupri)
!write (decomp%lupri, *) 'biter_m ='
!call mat_print (biter_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
!write (decomp%lupri, *) 'biter_v ='
!call mat_print (biter_v, 1, soeo%Nact, 1, 1, decomp%lupri)

  do i=1,iter
    call mat_read_from_disk (lub_m, bi_m)
    call mat_read_from_disk (lub_v, bi_v)

    Ared(i,iter) = soeo_dotproduct (bi_m, bi_v, sigmaiter_m, sigmaiter_v)

    call mat_read_from_disk (lusigma_m, sigmai_m)
    call mat_read_from_disk (lusigma_v, sigmai_v)

    Ared(iter,i) = soeo_dotproduct (biter_m, biter_v, sigmai_m, sigmai_v)

write (decomp%lupri, *) 'i =', i
write (decomp%lupri, *) 'sigmai_m ='
call mat_print (sigmai_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'sigmai_v ='
call mat_print (sigmai_v, 1, soeo%Nact, 1, 1, decomp%lupri)
write (decomp%lupri, *) 'bi_m ='
call mat_print (bi_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'bi_v ='
call mat_print (bi_v, 1, soeo%Nact, 1, 1, decomp%lupri)
write (decomp%lupri, *) 'iter, i, Ared(i,iter), Ared(iter,i) :', iter, i, Ared(i,iter), Ared(iter,i)
!write (decomp%lupri, *) 'Fmo ='
!call mat_print (soeo%Fmo, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)

  enddo

  gradred(iter) = soeo_dotproduct (biter_m, biter_v, grad_m, grad_v)
  !---------------------------------------------------------------------
enddo

!Finally setting the output:
!-----------------------------------------------------------------------
K = X_m
deltatheta = X_v
!-----------------------------------------------------------------------

!Finalizing
!-----------------------------------------------------------------------
call lsclose (lub_m, "DELETE")
call lsclose (lub_v, "DELETE")
call lsclose (lusigma_m, "DELETE")
call lsclose (lusigma_v, "DELETE")

call mat_free (grad_m     )
call mat_free (grad_v     )
call mat_free (b1_m       )
call mat_free (b1_v       )
call mat_free (bi_m       )
call mat_free (bi_v       )
call mat_free (biter_m    )
call mat_free (biter_v    )
call mat_free (sigma1_m   )
call mat_free (sigma1_v   )
call mat_free (sigmai_m   )
call mat_free (sigmai_v   )
call mat_free (sigmaiter_m)
call mat_free (sigmaiter_v)
call mat_free (X_m        )
call mat_free (X_v        )
call mat_free (res_m      )
call mat_free (res_v      )
!-----------------------------------------------------------------------

end subroutine soeo_solver
!=======================================================================

!> \brief Solves linear system Ax = b
!> \author C. Nygaard
!> \date 2010
!> \param iter Size of the linear system
!> \param Ain The matrix A (in)
!> \param bin The matrix b (in)
!> \param xout The matrix x (out)
!=======================================================================
subroutine soeo_solve_linear_system (iter, Ain, mu, bin, xout, decomp)
implicit none

!I/O:
type(decompItem), intent(in) :: decomp
integer, intent(in)          :: iter
real(realk), intent(in)      :: Ain(:,:), bin(:), mu
real(realk), intent(out)     :: xout(:)
!Other:
integer                      :: INFO, i
real(realk)                  :: A(iter,iter), b(iter)
integer                      :: IPIV(iter,iter)

xout = 0.0d0
A = Ain(1:iter,1:iter)
do i=1,iter
  A(i,i) = A(i,i) - mu
enddo
b = bin(1:iter)

!See man dgesv
call dgesv (iter, 1, A, iter, IPIV, b, iter, INFO)

if (INFO < 0) then
  !INFO = -i : the i'th argument had an illegal value
  write (decomp%lupri, '(/A, i4)') "Illegal argument in dgesv, INFO =", INFO
  call lsquit(" Problem in dgesv", decomp%lupri)
elseif (INFO > 0) then
  !INFO = i : U(i,i) is exactly zero
  !           A is singular and solution cannot be computed
  write (decomp%lupri, '(/A)') "Problem in dgesv"
  write (decomp%lupri, '(/A, i4)') "A is singular, INFO =", INFO
  call lsquit(" Problem in dgesv", decomp%lupri)
else
  !INFO = 0 : successful exit of dgesv
  xout(1:iter) = b
endif

end subroutine soeo_solve_linear_system
!=======================================================================

!> \brief Finds X from the redspace solutions and trial-vectors
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Size of the reduced space
!> \param lub_m LU for file containing matrix-part of trial vectors
!> \param lub_v LU for file containing vector-part of trial vectors
!> \param xred Solution to the reduced linear system
!> \param X_m Matrix part of the solution to the linear system in the reduced space
!> \param X_v Vector part of the solution to the linear system in the reduced space
!=======================================================================
subroutine soeo_solver_getX (iter, lub_m, lub_v, xred, X_m, X_v)

implicit none

integer, intent(in)         :: iter, lub_m, lub_v
real(realk), intent(in)     :: xred(:)
type(matrix), intent(inout) :: X_m, X_v

integer                     :: i, Nbast, Nact
type(matrix)                :: bi_m, bi_v

Nbast = X_m%nrow
Nact = X_v%nrow

call mat_init (bi_m, Nbast, Nbast)
call mat_init (bi_v, Nact, 1)

call mat_zero (X_m)
call mat_zero (X_v)
rewind (lub_m) ; rewind (lub_v)

do i=1,iter
  !Finding new X's (X = sum(xred(i)*bi)):
  call mat_read_from_disk (lub_m, bi_m)
  call mat_read_from_disk (lub_v, bi_v)
  call mat_daxpy (xred(i), bi_m, X_m)
  call mat_daxpy (xred(i), bi_v, X_v)
enddo

call mat_free (bi_m)
call mat_free (bi_v)

end subroutine soeo_solver_getX
!=======================================================================

!> \brief Finds res from the redspace solutions and trial-vectors
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Size of the reduced space
!> \param lusigma_m LU for file containing matrix-part of sigma
!> \param lusigma_v LU for file containing vector-part of sigma
!> \param xred Solution to the reduced linear system
!> \param grad_m Matrix part of the gradient
!> \param grad_v Vector part of the gradient
!> \param red_m Matrix part of the residual
!> \param red_v Vector part of the residual
!=======================================================================
subroutine soeo_solver_getres (iter, lusigma_m, lusigma_v, xred, grad_m,&
                             & grad_v, res_m, res_v)

implicit none

integer, intent(in)         :: iter, lusigma_m, lusigma_v
real(realk), intent(in)     :: xred(:)
type(matrix), intent(in)    :: grad_m, grad_v
type(matrix), intent(inout) :: res_m, res_v

integer                     :: i, Nbast, Nact
type(matrix)                :: sigmai_m, sigmai_v

Nbast = res_m%nrow
Nact = res_v%nrow
call mat_init (sigmai_m, Nbast, Nbast)
call mat_init (sigmai_v, Nact, 1)

call mat_zero (res_m)
call mat_zero (res_v)
rewind (lusigma_m) ; rewind (lusigma_v)

do i=1,iter
  !and residuals (res = HX - muX + grad = sum(xred(i)*sigmai) + grad):
  call mat_read_from_disk (lusigma_m, sigmai_m)
  call mat_read_from_disk (lusigma_v, sigmai_v)
  call mat_daxpy (xred(i), sigmai_m, res_m)
  call mat_daxpy (xred(i), sigmai_v, res_v)
enddo
call mat_daxpy (1.0d0, grad_m, res_m)
call mat_daxpy (1.0d0, grad_v, res_v)

call mat_free (sigmai_m)
call mat_free (sigmai_v)

end subroutine soeo_solver_getres
!=======================================================================

!> \brief Finds eigenvalues and eigenvectors for A 
!> \author C. Nygaard
!> \date 2010
!> \param A The matrix A (in)
!> \param evecs The eigenvectors (out)
!> \param evals The eigenvalues (out)
!=======================================================================
subroutine soeo_solve_eigeneqs (A, evecs, evals, decomp)
implicit none

!I/O
type(decompItem), intent(in) :: decomp
real(realk), intent(in)      :: A(:,:)
real(realk), intent(out)     :: evecs(:,:), evals(:)
!Other
integer                      :: i, j, N, IERR
integer, allocatable         :: IV1(:)
real(realk), allocatable     :: er(:), ei(:), X(:,:), FV1(:), FV2(:)
real(realk)                  :: tmp
logical                      :: symmetric

N = size(A(:,1))

allocate (er(N), ei(N), X(N,N), FV1(N), FV2(N), IV1(N))

symmetric = .true.
do i=1,N
  do j=1,i-1
    tmp = A(i,j) - A(j,i)
    if (tmp > 1.0d-5) then
      symmetric = .false.
    endif
  enddo
enddo

ei = 0.0d0
if (symmetric) then
  call RS (N, N, A, er, 1, X, FV1, FV2, IERR)
else
  call RG (N, N, A, er, ei, 1, X, IV1, FV1, IERR)
endif

        do i=1,N
          if (ei(i) > 1.0d-5 .or. ei(i) < -1.0d-5) then
           write (decomp%lupri, '("ei(", i3, ") =")') i
           write (decomp%lupri, *) ei(i)
          endif
        enddo

do i=1,N
  if (ei(i) > 1.0d-3 .or. ei(i) < -1.0d-3) then
    write (decomp%lupri, '("Hessian has complex eigenvalues!")')
    call lsquit ("Hessian has complex eigenvalues!", decomp%lupri)
  endif
enddo
if (IERR /= 0) then
  write (decomp%lupri, '("Error in RG (soeo_solve_eigeneqs)")')
  call lsquit ("Error in RG / RS", decomp%lupri)
endif

evecs = X
evals = er

deallocate (ei, er, X, FV1, IV1)

end subroutine soeo_solve_eigeneqs
!=======================================================================

!> \brief Finds the minimum eigenvalue of a matrix Ain
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Dimension of Ain
!> \param Ain The matrix to which we want to find mineval
!> \param mineval The minimum eigenvalue of Ain
!=======================================================================
subroutine soeo_find_mineval (lupri, iter, Ain, mineval)

implicit none

!I/O
integer, intent(in)      :: iter, lupri
real(realk), intent(in)  :: Ain(:,:)
real(realk), intent(out) :: mineval
!Other
integer                  :: i, j, IERR
integer                  :: IV1(iter)
real(realk)              :: er(iter), ei(iter), X(iter,iter), FV1(iter), FV2(iter)
real(realk)              :: A(iter,iter), tmp
logical                  :: symmetric

A = Ain(1:iter,1:iter)

symmetric = .true.
do i=1,iter
  do j=1,i-1
    tmp = A(i,j) - A(j,i)
    if (tmp > 1.0d-5 .or. tmp < -1.0d-05) then
      symmetric = .false.
    endif
  enddo
enddo

ei = 0.0d0
if (symmetric) then
  call RS (iter, iter, A, er, 1, X, FV1, FV2, IERR)
else
  call RG (iter, iter, A, er, ei, 1, X, IV1, FV1, IERR)
endif

do i=1,iter
  if (ei(i) > 1.0d-3 .or. ei(i) < -1.0d-3) then
    print *, "Redspace has complex eigenvalues!"
    call lsquit ("Redspace has complex eigenvalues!", lupri)
  endif
enddo
if (IERR /= 0) then
  print *, "Error in RG / RS (soeo_find_mineval)"
  call lsquit ("Error in RG / RS")
endif

mineval = 1.0d5
do i=1,iter
  if (er(i) < mineval) then
    mineval = er(i)
  endif
enddo

end subroutine soeo_find_mineval
!=======================================================================

!> \brief Choose limits for the binary search to find mu
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Dimension of the reduced space
!> \param Nbast Size of matrix-part
!> \param Nact Size of vector-part
!> \param lub_m LU for file containing matrix-part of trial vectors
!> \param lub_v LU for file containing vector-part of trial vectors
!> \param Ared The reduced space hessian
!> \param gradred The reduced space gradient
!> \param mine Minimum eigenvalue of Ared
!> \param h The trust-radius
!> \param mu The levelshift to control steplength
!=======================================================================
subroutine soeo_binary_search_init (iter, Nbast, Nact, lub_m, lub_v,&
                                  & Ared, gradred, mine, h, mu, decomp)

implicit none

type(decompItem), intent(in) :: decomp
integer, intent(in)          :: iter, Nbast, Nact, lub_m, lub_v
real(realk), intent(in)      :: Ared(:,:), gradred(:), mine, h
real(realk), intent(out)     :: mu

if (mine >= -1.0d-3) then
  !Search from 0 and down
  write (decomp%lupri, *) 'mineval >= -1,0d-3 -- search from 0 and down'
  call soeo_binary_search (iter, Nbast, Nact, lub_m, lub_v,&
                         & Ared, gradred, 0.0d0, -16.0d0, h, mu, decomp)
  write (decomp%lupri, *) 'search done, mu =', mu
elseif (mine < -1.0d-3 .and. iter > 1) then
  !Search from mine and down
  write (decomp%lupri, *) 'mineval < -1,0d-3 and iter > 1 -- search from mineval and down'
  call soeo_binary_search (iter, Nbast, Nact, lub_m, lub_v,&
                         & Ared, gradred, mine, 10.0d0*mine, h, mu, decomp)
  write (decomp%lupri, *) 'search done, mu =', mu
else
  !Search from 2*mine and down
  write (decomp%lupri, *) 'mineval < -1.0d-3 and iter = 1 -- search from 2*mineval and down'
  call soeo_binary_search (iter, Nbast, Nact, lub_m, lub_v,&
                         & Ared, gradred, 2.0d0*mine, 10.0d0*mine, h, mu, decomp)
  write (decomp%lupri, *) 'search done, mu =', mu
endif

end subroutine soeo_binary_search_init
!=======================================================================

!> \brief Performs binary search to find mu
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Dimension of the reduced space
!> \param Nbast Size of matrix-part
!> \param Nact Size of vector-part
!> \param lub_m LU for file containing matrix-part of trial vectors
!> \param lub_v LU for file containing vector-part of trial vectors
!> \param Ared The reduced space hessian
!> \param gradred The reduced space gradient
!> \param high The high limit for the binary search
!> \param low The low limit for the binary search
!> \param h The trust-radius
!> \param mu The levelshift to control steplength
!=======================================================================
subroutine soeo_binary_search (iter, Nbast, Nact, lub_m, lub_v, Ared, gradred,&
                             & highin, lowin, h, mu, decomp)

implicit none

type(decompItem), intent(in) :: decomp
integer, intent(in)          :: iter, Nbast, Nact, lub_m, lub_v
real(realk), intent(in)      :: Ared(:,:), gradred(:), h
real(realk), intent(in)      :: highin, lowin
real(realk), intent(out)     :: mu

integer                      :: N, ct
real(realk), allocatable     :: xred(:)
real(realk)                  :: high, middle, low, normh, normm, norml
type(matrix)                 :: X_m, X_v

N = size(gradred)
allocate (xred(N))

high = highin ; low = lowin

call mat_init (X_m, Nbast, Nbast)
call mat_init (X_v, Nact, 1)

call soeo_solve_linear_system (iter, Ared, high, gradred, xred, decomp)
call soeo_solver_getX (iter, lub_m, lub_v, xred, X_m, X_v)

normh = soeo_norm (X_m, X_v)

write (decomp%lupri, *) 'h, normh =', h, normh

if (normh <= h) then
  !Step will allways fall inside trust-region
  mu = high
else

  !check if wanted mu is in the given interval
  ct = 0
  do
    ct = ct + 1
    call soeo_solve_linear_system (iter, Ared, low, gradred, xred, decomp)
    call soeo_solver_getX (iter, lub_m, lub_v, xred, X_m, X_v)
    norml = soeo_norm (X_m, X_v)
    if (norml > h) then
      high = low
      normh = norml
      low = low * 2.0d0
      write (decomp%lupri, *) 'mu is not in given interval - interval expanded'
    else
      exit
    endif
    if (ct > 100) then
      write (decomp%lupri, *) 'WARNING: mu is not in given interval!'
      exit
    endif
  enddo

  !performing the search
  ct = 0
  do
    ct = ct + 1

    middle = (high + low)*0.5d0

    call soeo_solve_linear_system (iter, Ared, middle, gradred, xred, decomp)
    call soeo_solver_getX (iter, lub_m, lub_v, xred, X_m, X_v)
    normm = soeo_norm (X_m, X_v)

write (decomp%lupri, *) 'low, normlow :', low, norml
write (decomp%lupri, *) 'middle, normmiddle :', middle, normm
write (decomp%lupri, *) 'high, normhigh :', high, normh

    if (normm - h < 1.0d-5 .and. normm - h > -1.0d-5) then
      mu = middle
      exit
    elseif (normm - h < 0.0d0) then
      low = middle
      norml = normm
    elseif (normm - h > 0.0d0) then
      high = middle
      normh = normm
    endif

    if (ct > 200) exit

  enddo
endif

call mat_free (X_m)
call mat_free (X_v)

deallocate (xred)

end subroutine soeo_binary_search
!=======================================================================

end module soeo_redspace
