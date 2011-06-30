!> @file
!> Contains matrix type definition module.

!> \brief Matrix type definitions.
!>
!> General rules:
!> NEVER put a type(Matrix) as intent(out), this will on some platforms
!>       make the pointer 
!>       disassociated entering the routine, and memory already 
!>       allocated for the matrix will be lost. \n
!> NEVER think that e.g. A%elms = matrix will copy the matrix elements from 
!>       matrix to A%elms, it will only associate the pointer with the array
!>       matrix. \n
!> BUT type(Matrix) :: A,B; A = B SHOULD copy the matrix elements from matrix B to A
!>     see mat_assign. \n
!> ALWAYS and ONLY call mat_free on a matrix you have initialized with mat_init
!>
MODULE Matrix_module
   use precision
   !> Logical unit number of file LSDALTON.OUT used by matrix operation modules
   integer, save :: mat_lu
   !> True if various info from matrix module should be printed
   logical, save :: mat_info
   !> True if number of allocated matrices should be monitored
   logical, save :: mat_mem_monitor
   !> A simple type of sparse matrix used inside sparse1 module
   TYPE SMAT1
      !> Row index
      integer :: i
      !> Column index
      integer :: j
      !> Value of element
      real(realk) :: val
   END TYPE SMAT1
   !> The general matrix type that is used throughout LSDALTON
   TYPE Matrix
      !ajt jan09 Mark uninitialized matrices with nrow=nrol=-1
      !> number of rows
      INTEGER :: nrow = -1
      !> number of columns
      INTEGER :: ncol = -1
      !> room for various flags/information for the matrix
      integer, dimension(:),pointer :: idata
      !> bsm permutation pointer
      integer, dimension(:),pointer :: permutation
      !> pointer to double precision matrix element storage
      real(realk), pointer   :: elms(:)
      !> pointer to double precision matrix element storage for beta part
      real(realk), pointer   :: elmsb(:)
      !> pointer to complex matrix element storage
      complex(realk), pointer :: celms(:)
      !> pointer to complex matrix element storage for beta part
      complex(realk), pointer :: celmsb(:)
      !> pointer to storage scheme 1 for sparse matrices
      type(SMAT1),pointer :: selm1(:)
      !pointer to storage scheme 2 for sparse matrices
      !type(SMAT2),pointer :: selm2
      !> If only blocks of the matrix are to be stored: array of pointers to the matrix-blocks
      type(Matrix),pointer :: block(:)
      !> the indexes where the blocks start
      integer,dimension(:,:), pointer :: blockpos!(2,:)
      !> room for any integer auxiliary information
      integer, pointer     :: iaux(:)
      !> room for any real auxiliary information
      real(realk), pointer :: raux(:)
      !ajt feb09 to mark complex matrices
      !> marker for complex elements
      logical :: complex = .false.  
      real(realk),pointer :: val(:)
      integer,pointer :: col(:)
      integer,pointer :: row(:)
      integer :: nnz

#ifdef PRG_DIRAC
      !> hermiticity
      integer :: ih_sym = 1
      !> time-reversal symmetry
      integer :: tr_sym = 1
      !> irep
      integer :: irep   = 0
#endif

   END TYPE Matrix

   !> Pointer to a type(matrix). Necessary if we want arrays of derived types!
   type Matrixp
      TYPE(Matrix), pointer :: p
   end type Matrixp

END MODULE Matrix_module

