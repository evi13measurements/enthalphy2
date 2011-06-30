!> @file 
!> Contains the LSmatrix structure 
MODULE LSmatrix_type
use precision
TYPE LSMATRIX
INTEGER :: nrow
INTEGER :: ncol
real(realk), pointer   :: elms(:)
logical :: complex = .FALSE.
END TYPE LSMATRIX

type LSMatrixpointer
TYPE(LSMatrix), pointer :: p
end type LSMatrixpointer

END MODULE LSmatrix_type
