module soeo_matop

use Matrix_Operations
use matrix_util
use dal_interface
use decompMod

implicit none

!Contains:
!  soeo_diag_elma
!  soeo_diag_elmb
!  soeo_dotproduct
!  soeo_norm
!  soeo_normalize
!  soeo_scal
!  soeo_orthonormalize
!  soro_orthogonalize
!  soeo_daxpy
!  soeo_get_GmoL

contains

!> \brief Gives the diagonal element (i,i) of a matrix A (alpha part)
!> \author C. Nygaard
!> \date 2010
!> \param A The matrix
!> \param i The place of the diagonal element
!=======================================================================
function soeo_diag_elma (A, i)

implicit none

real(realk) soeo_diag_elma
type(matrix), intent(in) :: A
integer, intent(in)      :: i
real(realk)              :: tmp(2)

call mat_get_ab_elms (A, i, i, tmp)
soeo_diag_elma = tmp(1)

end function soeo_diag_elma
!=======================================================================

!> \brief Gives the diagonal element (i,i) of a matrix A (beta part)
!> \author C. Nygaard
!> \date 2010
!> \param A The matrix
!> \param i The place of the diagonal element
!=======================================================================
function soeo_diag_elmb (A, i)

implicit none

real(realk) soeo_diag_elmb
type(matrix), intent(in) :: A
integer, intent(in)      :: i
real(realk)              :: tmp(2)

call mat_get_ab_elms (A, i, i, tmp)
soeo_diag_elmb = tmp(2)

end function soeo_diag_elmb
!=======================================================================

!> \brief Calculates the dot-product of a matrix-vector
!> \author C. Nygaard
!> \date 2010
!> \param a_m Matrix part of the first matrix-vector
!> \param a_v Vector part of the first matrix-vector
!> \param b_m Matrix part of the second matrix-vector
!> \param b_v Vector part of the second matrix-vector
!   a_m, a_v dot b_m, b_v
!>  soeo_dotproduct = dot(a_m, b_m) + dot(a_v, b_v)
!=======================================================================
function soeo_dotproduct (a_m, a_v, b_m, b_v)

implicit none

!I/O
real(realk) soeo_dotproduct
type(matrix), intent(in) :: a_m, b_m, a_v, b_v
!Other
real(realk)              :: dot_m, dot_v

!Works because the m-part is allways antisymmetric, so that it contains
! half as many independent variables as nonzero elements

dot_m = 0.5d0*mat_dotproduct (a_m, b_m)
dot_v = mat_dotproduct (a_v, b_v)
soeo_dotproduct = dot_m + dot_v

end function soeo_dotproduct
!=======================================================================

!> \brief Calculates the norm of a matrix-vector
!> \author C. Nygaard
!> \date 2010
!> \param x_m Matrix part of the matrix-vector
!> \param x_v Vector part of the matrix-vector
!=======================================================================
function soeo_norm (x_m, x_v)

implicit none

!I/O
real(realk) soeo_norm
type(matrix), intent(in) :: x_m, x_v
!Other
real(realk)              :: sqnorm

sqnorm = soeo_dotproduct (x_m, x_v, x_m, x_v)
soeo_norm = dsqrt(sqnorm)

end function soeo_norm
!=======================================================================

!> \brief Normalizes a matrix-vector
!> \author C. Nygaard
!> \date 2010
!> \param a_m Matrix part of the matrix-vector
!> \param a_v Vector part of the matrix-vector
!=======================================================================
subroutine soeo_normalize (a_m, a_v)

implicit none

!I/O
type(matrix), intent(inout) :: a_m, a_v
!Other
real(realk)                 :: anorm

anorm = soeo_norm (a_m, a_v)
call soeo_scal (1.0d0/anorm, a_m, a_v)

end subroutine soeo_normalize
!=======================================================================

!> \brief Multiplies a matrix-vector with a scalar alpha
!> \author C. Nygaard
!> \date 2010
!> \param alpha The scalar
!> \param a_m Matrix part of the matrix-vector
!> \param a_v Vector part of the matrix-vector
!=======================================================================
subroutine soeo_scal (alpha, a_m, a_v)

implicit none

!I/O
real(realk), intent(in)        :: alpha
type(matrix), intent(inout)    :: a_m, a_v
!Other
integer                        :: N

call mat_scal (alpha, a_m)
call mat_scal (alpha, a_v)

end subroutine soeo_scal
!=======================================================================

!> \brief Orthonormalizes matrix-vector b to matrix-vector a
!> \author C. Nygaard
!> \date 2010
!> \param a_m Matrix part of the first matrix-vector (in)
!> \param a_v Vector part of the first matrix-vector (in)
!> \param b_m Matrix part of the second matrix-vector (inout)
!> \param b_v Vector part of the second matrix-vector (inout)
!> \param err Is true if norm of b after orthogonalization is too small
!=======================================================================
subroutine soeo_orthonormalize (a_m, a_v, b_m, b_v, err)

implicit none

!I/O
type(matrix), intent(in)    :: a_m, a_v
type(matrix), intent(inout) :: b_m, b_v
logical, intent(out)        :: err
!Other
real(realk)                 :: proj, normb

err = .false.

proj = soeo_dotproduct (a_m, a_v, b_m, b_v)
call soeo_daxpy (-proj, a_m, a_v, b_m, b_v) !b = b - bproj

!print *, 'inside soeo_orthonormalize'

normb = soeo_norm (b_m, b_v)

!print *, 'normb =', normb

if (normb < 1.0d-9) then
  err = .true.
else
  call soeo_normalize (b_m, b_v)
endif

end subroutine soeo_orthonormalize
!=======================================================================

!> \brief Makes matrix-vector b orthogonal to matrix-vector a
!> \author C. Nygaard
!> \date 2010
!> \param a_m Matrix part of the first matrix-vector (in)
!> \param a_v Vector part of the first matrix-vector (in)
!> \param b_m Matrix part of the second matrix-vector (inout)
!> \param b_v Vector part of the second matrix-vector (inout)
!> \param err Is true if norm of b after orthogonalization is too small
!=======================================================================
subroutine soeo_orthogonalize (a_m, a_v, b_m, b_v, err)

implicit none

!I/O
type(matrix), intent(in)    :: a_m, a_v
type(matrix), intent(inout) :: b_m, b_v
logical, intent(out)        :: err
!Other
real(realk)                 :: num, denom, proj, normb

err = .false.

num = soeo_dotproduct (a_m, a_v, b_m, b_v)
denom = soeo_dotproduct (a_m, a_v, a_m, a_v)
proj = num / denom
call soeo_daxpy (-proj, a_m, a_v, b_m, b_v) !b = b - bproj

normb = soeo_norm (b_m, b_v)
if (normb < 1.0d-9) then
  err = .true.
endif

end subroutine soeo_orthogonalize
!=======================================================================

!> \brief Calculates y = alpha*x + y for matrix-vectors
!> \author C. Nygaard
!> \date 2010
!> \param alpha The scalar
!> \param x_m Matrix part of the first matrix-vector (in)
!> \param x_v Vector part of the first matrix-vector (in)
!> \param y_m Matrix part of the second matrix-vector (inout)
!> \param y_v Vector part of the second matrix-vector (inout)
!=======================================================================
subroutine soeo_daxpy (alpha, x_m, x_v, y_m, y_v)

implicit none

!I/O
real(realk), intent(in)     :: alpha
type(matrix), intent(in)    :: x_m, x_v
type(matrix), intent(inout) :: y_m, y_v

if (x_m%nrow == y_m%nrow .and. x_m%ncol == y_m%ncol &
                       & .and. x_v%nrow == y_v%nrow) then

  call mat_daxpy (alpha, x_m, y_m)
  call mat_daxpy (alpha, x_v, y_v)

else
  call lsquit ("Wrong dimensions in soeo_daxpy",-1)
endif

end subroutine soeo_daxpy
!=======================================================================

end module soeo_matop
