!> @file 
!> Contains matrix operations module and standalone BSM routines.

!> Contains wrapper routines that branch out to matrix routine for chosen matrix type.
!> \author L. Thogersen
!> \date 2003
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
!>     (see mat_assign). \n
!> ALWAYS and ONLY call mat_free on a matrix you have initialized with mat_init.
!>
MODULE matrix_operations
!FIXME: order routines alphabetically
!   use lstiming
   use matrix_module
!   Use matrix_operations_symm_dense
   use matrix_operations_dense
#ifndef UNITTEST
   Use matrix_operations_sparse1
#endif
   use matrix_operations_csr
!   Use matrix_operations_unres_symm_dense
   use matrix_operations_unres_dense
!   use matrix_operations_unres_sparse1
!> Matrices are symmetric and dense (not implemented)
   integer, parameter :: mtype_symm_dense = 1
!> Matrices are dense (default) 
   integer, parameter :: mtype_dense = 2
!> Matrices are block sparse (BSM)
   integer, parameter ::  mtype_sparse_block = 3
!> Matrices are compressed sparse row (CSR) sparse
   integer, parameter ::  mtype_sparse1 = 4
!> Matrices are dense and have both alpha and beta part (default for open shell)
   integer, parameter ::  mtype_unres_dense = 5
!> Matrices are CSR sparse and have both alpha and beta part (not implemented)
   integer, parameter ::  mtype_unres_sparse1 = 6
!> Matrices are compressed sparse row (CSR) 
   integer, parameter ::  mtype_csr = 7
!> Counts the total number of matrix multiplications used throughout calculation
   integer, save      :: no_of_matmuls = 0
!> Counts the number of matrices allocated, if requested
   integer, save      :: no_of_matrices = 0 
!> Tracks the maximum number of allocated matrices throughout calculation, if requested
   integer, save      :: max_no_of_matrices = 0
!> Start time for timing of matrix routines (only if requested)
   real(realk), save  :: mat_TSTR
!> End time for timing of matrix routines (only if requested)
   real(realk), save  :: mat_TEN
!*****************
!Possible matrix types - 
!(Exploiting symmetry when operating sparse matrices is probably a vaste of effort - 
!therefore these combinations are removed)
!******************
!mtype_dense, mtype_sparse1, mtype_sparse2
!mtype_cplx_dense, mtype_cplx_sparse1, mtype_cplx_sparse2
!mtype_unres_dense, mtype_unres_sparse1, mtype_unres_sparse2
!mtype_cplx_unres_dense, mtype_cplx_unres_sparse1, mtype_cplx_unres_sparse2
!mtype_symm_dense
!mtype_cplx_symm_dense
!mtype_unres_symm_dense
!mtype_cplx_unres_symm_dense
!> This is set to one of the mtype... variables to indicate chosen matrix type
   integer,save :: matrix_type = mtype_dense !default dense
!> True if timings for matrix operations are requested
   logical,save :: INFO_TIME_MAT = .false. !default no timings
   logical,save :: INFO_memory = .false. !default no memory printout
!> Overload: The '=' sign may be used to set two type(matrix) structures equal, i.e. A = B
   INTERFACE ASSIGNMENT(=)
      module procedure mat_assign
   END INTERFACE

   contains
!*** is called from config.f90
!> \brief Sets the global variable matrix_type that determines the matrix type
!> \author L. Thogersen
!> \date 2003
!> \param a Indicates the matrix type (see module documentation) 
      SUBROUTINE mat_select_type(a)
         implicit none
         INTEGER, INTENT(IN) :: a
         matrix_type = a
      END SUBROUTINE mat_select_type
!> \brief Pass info about e.g. logical unit number for LSDALTON.OUT to matrix module 
!> \author L. Thogersen
!> \date 2003
!> \param lu_info Logical unit number for LSDALTON.OUT
!> \param info_info True if various info from matrix module should be printed
!> \param mem_monitor True if number of allocated matrices should be monitored
      SUBROUTINE mat_pass_info(lu_info,info_info,mem_monitor)
         implicit none
         integer, intent(in) :: lu_info
         logical, intent(in) :: info_info, mem_monitor
         mat_lu   = lu_info
         mat_info = info_info
         mat_mem_monitor = mem_monitor
      END SUBROUTINE mat_pass_info
!> \brief If called, timings from matrix routines will be printed 
!> \author L. Thogersen
!> \date 2003
      SUBROUTINE mat_timings
         implicit none
         INFO_TIME_MAT = .true. 
      END SUBROUTINE mat_timings
!***
!> \brief Returns the number of matrix multiplications used so far 
!> \author S. Host
!> \date 2009
!> \param n Number of matrix muliplications
      SUBROUTINE mat_no_of_matmuls(n)
         implicit none
         INTEGER, INTENT(out) :: n
         n = no_of_matmuls
      END SUBROUTINE mat_no_of_matmuls
!> \brief Initialize a type(matrix)
!> \author L. Thogersen
!> \date 2003
!> \param a type(matrix) that should be initialized
!> \param nrow Number of rows for a
!> \param ncol Number of columns for a
      SUBROUTINE mat_init(a, nrow, ncol)
         implicit none
         TYPE(Matrix) :: a 
         INTEGER, INTENT(IN)   :: nrow, ncol
         if (mat_mem_monitor) then
            no_of_matrices = no_of_matrices + 1
            !write(mat_lu,*) 'Init: matrices allocated:', no_of_matrices
            if (no_of_matrices > max_no_of_matrices) max_no_of_matrices = no_of_matrices
         endif
         nullify(A%elms)
         nullify(A%elmsb)
         if (info_memory) write(mat_lu,*) 'Before mat_init: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_init(a)
         case(mtype_dense)
             call mat_dense_init(a,nrow,ncol)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_init(a,nrow,ncol)
        case(mtype_sparse_block)
#ifdef HAVE_BSM
           CALL bsm_init(a,nrow,ncol)
           a%nrow = nrow
           a%ncol = ncol
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_init(a)
         case(mtype_unres_dense)
             call mat_unres_dense_init(a,nrow,ncol)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_init(a)
         case(mtype_csr)
             call mat_csr_init(a,nrow,ncol)            
         case default
              stop "mat_init not implemented for this type of matrix"
         end select

         NULLIFY(a%iaux, a%raux)
         if (info_memory) write(mat_lu,*) 'After mat_init: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_init

!> \brief Free a type(matrix) that has been initialized with mat_init
!> \author L. Thogersen
!> \date 2003
!> \param a type(matrix) that should be freed
      SUBROUTINE mat_free(a)
         implicit none
         TYPE(Matrix) :: a 
         if (mat_mem_monitor) then
            no_of_matrices = no_of_matrices - 1
            !write(mat_lu,*) 'Free: matrices allocated:', no_of_matrices
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_free: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_free(a)
         case(mtype_dense)
             call mat_dense_free(a)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_free(a)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
             call bsm_free(a)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_free(a)
         case(mtype_unres_dense)
             call mat_unres_dense_free(a)
         case(mtype_csr)
             call mat_csr_free(a)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_free(a)
         case default
              stop "mat_free not implemented for this type of matrix"
         end select

      !free auxaliary data
      if (ASSOCIATED(a%iaux)) deallocate(a%iaux)
      if (ASSOCIATED(a%raux)) deallocate(a%raux)

      a%nrow = -1; a%ncol = -1
      if (info_memory) write(mat_lu,*) 'After mat_free: mem_allocated_global =', mem_allocated_global

      END SUBROUTINE mat_free

!> \brief Count allocated memory for type(matrix)
!> \author L. Thogersen
!> \date 2003
!> \param nsize Number of real(realk) elements that have been allocated
      SUBROUTINE stat_allocated_memory(nsize)
         implicit none
         integer, intent(in) :: nsize 
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call symm_dense_stat_allocated_memory(nsize)
         case(mtype_dense)
             call dens_stat_allocated_memory(nsize)
#ifndef UNITTEST
         case(mtype_sparse1)
             call sp1_stat_allocated_memory(nsize)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
!             call bsm_stat_allocated_memory(nsize)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_free(nsize)
         case(mtype_unres_dense)
             call unres_dens_stat_allocated_memory(nsize)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_stat_allocated_memory(nsize)
         case default
              stop "stat_allocated_memory not implemented for this type of matrix"
         end select

      END SUBROUTINE stat_allocated_memory

!> \brief Count deallocated memory for type(matrix)
!> \author L. Thogersen
!> \date 2003
!> \param nsize Number of real(realk) elements that have been deallocated
      SUBROUTINE stat_deallocated_memory(nsize)
         implicit none
         integer, intent(in) :: nsize 
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call symm_dense_stat_deallocated_memory(nsize)
         case(mtype_dense)
             call dens_stat_deallocated_memory(nsize)
#ifndef UNITTEST
         case(mtype_sparse1)
             call sp1_stat_deallocated_memory(nsize)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
!             call bsm_stat_deallocated_memory(nsize)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_free(nsize)
         case(mtype_unres_dense)
             call unres_dens_stat_deallocated_memory(nsize)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_stat_deallocated_memory(nsize)
         case default
              stop "stat_deallocated_memory not implemented for this type of matrix"
         end select

      END SUBROUTINE stat_deallocated_memory

!> \brief Convert a standard fortran matrix to a type(matrix) - USAGE DISCOURAGED!
!> \author L. Thogersen
!> \date 2003
!> \param afull Standard fortran matrix that should be converted (n x n)
!> \param alpha The output type(matrix) is multiplied by alpha
!> \param a The output type(matrix) ((2n x 2n) if unrestricted, (n x n) otherwise)
!> \param mat_label If the character label is present, sparsity will be printed if using block-sparse matrices
!> \param unres3 If present and true, mat_unres_dense_set_from_full3 will be called instead of mat_unres_dense_set_from_full
!>  
!> BE VERY CAREFUL WHEN USING mat_set_from_full AND mat_to_full!!!!!!
!> Usage of these routines should be avoided whenever possible, since
!> you have to hardcode an interface to make them work with unrestriced
!> matrices (see e.g. di_get_fock in dalton_interface.f90) This is because
!> usually, a and afull should have the same dimensions, but for unrestricted
!> a is (n x n) and afull is (2n x 2n). The exception is if
!> unres3 = true. In that case, both a and afull are (n x n), and both 
!> the alpha and beta parts of a will be set equal to the afull.
!> 
      SUBROUTINE mat_set_from_full(afull,alpha, a, mat_label,unres3)
         implicit none
         real(realk), INTENT(IN) :: afull(*)
         real(realk), intent(in) :: alpha
         TYPE(Matrix)            :: a  !output
         character(*), INTENT(IN), OPTIONAL :: mat_label
         logical,intent(in) , OPTIONAL :: unres3
         real(realk)             :: sparsity

         !write(mat_lu,*) "Usage of mat_set_from_full discouraged!!!"
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'Before mat_set_from_full: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_set_from_full(afull,alpha,a)
         case(mtype_dense)
             call mat_dense_set_from_full(afull,alpha,a)
         case(mtype_csr)
             call mat_csr_set_from_full(afull,alpha,a)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_set_from_full(afull,alpha,a)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            IF(ALPHA.NE.1D0)CALL DSCAL(a%nrow*a%ncol,ALPHA,afull,1)
            call bsm_free(a)
            CALL bsm_init_from_full(a,a%nrow,a%ncol,afull,sparsity)
            IF(ALPHA.NE.1D0)CALL DSCAL(a%nrow*a%ncol,1D0/ALPHA,afull,1)
            if(PRESENT(mat_label))&
                 &write(2,&
                 &'("BSM ",A," full->sparse, sparsity:",F6.1," % nnz:",F9.0)')&
                 & mat_label, sparsity*100D0, sparsity*a%nrow*a%ncol
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_set_from_full(afull,alpha,a)
         case(mtype_unres_dense)
            if(PRESENT(unres3))then
               if(unres3)then
                  call mat_unres_dense_set_from_full3(afull,alpha,a)
               else
                  call mat_unres_dense_set_from_full(afull,alpha,a)
               endif
            else
               call mat_unres_dense_set_from_full(afull,alpha,a)
            endif
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_set_from_full(afull,alpha,a)
         case default
              stop "mat_set_from_full not implemented for this type of matrix"
         end select
         if (info_memory) write(mat_lu,*) 'After mat_set_from_full: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('F_FULL',mat_TSTR,mat_TEN,mat_lu)
      END SUBROUTINE mat_set_from_full


!> \brief Unrestricted only! Convert a standard fortran matrix to a type(matrix) - USAGE DISCOURAGED!
!> \author L. Thogersen
!> \date 2003
!> \param afull (2n x 2n) standard fortran matrix that should be converted
!> \param alpha The output type(matrix) is multiplied by alpha
!> \param a The output (n x n) type(matrix)
!>  
      SUBROUTINE mat_set_from_full2(afull,alpha, a)
         implicit none
         real(realk), INTENT(IN) :: afull(*)
         real(realk), intent(in) :: alpha
         TYPE(Matrix)            :: a  !output
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_set_from_full2(afull,alpha,a)
!         case(mtype_dense)
!             call mat_dense_set_from_full2(afull,alpha,a)
!         case(mtype_sparse1)
!             call mat_sparse1_set_from_full2(afull,alpha,a)
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_set_from_full2(afull,alpha,a)
         case(mtype_unres_dense)
             call mat_unres_dense_set_from_full2(afull,alpha,a)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_set_from_full2(afull,alpha,a)
         case default
              stop "mat_set_from_full2 not implemented for this type of matrix"
         end select
      END SUBROUTINE mat_set_from_full2

!> \brief Convert a type(matrix) to a standard fortran matrix - USAGE DISCOURAGED!
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) that should be converted (n x n)
!> \param afull The output standard fortran matrix ((2n x 2n) if unrestricted, (n x n) otherwise)
!> \param alpha The output standard fortran matrix is multiplied by alpha
!> \param mat_label If the character label is present, sparsity will be printed if using block-sparse matrices
!>  
!> BE VERY CAREFUL WHEN USING mat_set_from_full AND mat_to_full!!!!!!
!> Usage of these routines should be avoided whenever possible, since
!> you have to hardcode an interface to make them work with unrestriced
!> matrices (see e.g. di_get_fock in dalton_interface.f90)
!>
     SUBROUTINE mat_to_full(a, alpha, afull,mat_label)

         implicit none
         TYPE(Matrix), intent(in):: a
         real(realk), intent(in) :: alpha
         real(realk), intent(out):: afull(*)  !output
         character(*), INTENT(IN), OPTIONAL :: mat_label
         real(realk)             :: sparsity
         
         !write(mat_lu,*) "Usage of mat_to_full discouraged!!!"
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         !if (SIZE(afull) < a%nrow*a%ncol) then
         !  STOP 'too small full array in mat_to_full'
         !endif
         if (info_memory) write(mat_lu,*) 'Before mat_to_full: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_to_full(a, alpha, afull)
         case(mtype_dense)
             call mat_dense_to_full(a, alpha, afull)
         case(mtype_csr)
             call mat_csr_to_full(a, alpha, afull)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_to_full(a, alpha, afull)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            call bsm_to_full(a, afull, sparsity)
            if(ALPHA.NE.1D0)CALL DSCAL(a%nrow*a%ncol, alpha, afull, 1)
            if(PRESENT(mat_label))&
                 &write(2,&
                 &'("BSM ",A," sparse->full, sparsity:",F6.1," % nnz:",F9.0)')&
                 & mat_label, sparsity*100D0, sparsity*a%nrow*a%ncol
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_to_full(a, alpha, afull)
         case(mtype_unres_dense)
             call mat_unres_dense_to_full(a, alpha, afull)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_to_full(a, alpha, afull)
         case default
              stop "mat_to_full not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('TOFULL',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_to_full: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_to_full

!> \brief Unrestricted only! Convert a type(matrix) to a standard fortran matrix - USAGE DISCOURAGED!
!> \author L. Thogersen
!> \date 2003
!> \param a The (n x n) type(matrix) that should be converted
!> \param afull The (2n x 2n) output standard fortran matrix 
!> \param alpha The output standard matrix is multiplied by alpha
!>  
      SUBROUTINE mat_to_full2(a, alpha, afull)
         implicit none
         TYPE(Matrix), intent(in):: a
         real(realk), intent(in) :: alpha
         real(realk), intent(out):: afull(*)  !output
         
         !if (SIZE(afull) < a%nrow*a%ncol) then
         !  STOP 'too small full array in mat_to_full2'
         !endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_to_full2(a, alpha, afull)
!         case(mtype_dense)
!             call mat_dense_to_full2(a, alpha, afull)
!         case(mtype_sparse1)
!             call mat_sparse1_to_full2(a, alpha, afull)
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_to_full2(a, alpha, afull)
         case(mtype_unres_dense)
             call mat_unres_dense_to_full2(a, alpha, afull)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_to_full2(a, alpha, afull)
         case default
              stop "mat_to_full2 not implemented for this type of matrix"
         end select
      END SUBROUTINE mat_to_full2

!> \brief Print a type(matrix) to file in pretty format
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) that should be printed
!> \param i_row1 Print starting from this row
!> \param i_rown Print ending at this row
!> \param j_col1 Print starting from this column
!> \param j_coln Print ending at this column
!> \param lu Print to file with this logical unit number
      SUBROUTINE mat_print(a, i_row1, i_rown, j_col1, j_coln, lu)
         implicit none
         TYPE(Matrix),intent(in) :: a
         integer, intent(in)     :: i_row1, i_rown, j_col1, j_coln, lu 
         REAL(REALK), ALLOCATABLE :: afull(:,:)
         real(realk)              :: sparsity

         if (i_row1 < 1 .or. j_col1 < 1 .or. a%nrow < i_rown .or. a%ncol < j_coln) then
           STOP 'subsection out of bounds in mat_print'
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_print: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_print(a, i_row1, i_rown, j_col1, j_coln, lu)
         case(mtype_dense)
             call mat_dense_print(a, i_row1, i_rown, j_col1, j_coln, lu)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_print(a, i_row1, i_rown, j_col1, j_coln, lu)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            ALLOCATE (afull(a%nrow,a%ncol))
            call bsm_to_full(a, afull,sparsity)
            CALL OUTPUT(afull, i_row1, i_rown, j_col1, j_coln,A%nrow,A%ncol,1, lu)
            DEALLOCATE(afull)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_print(a, i_row1, i_rown, j_col1, j_coln, lu)
         case(mtype_unres_dense)
             call mat_unres_dense_print(a, i_row1, i_rown, j_col1, j_coln, lu)
         case(mtype_csr)
             call mat_csr_print(a, lu)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_print(a, i_row1, i_rown, j_col1, j_coln, lu)
         case default
              stop "mat_print not implemented for this type of matrix"
         end select
         if (info_memory) write(mat_lu,*) 'After mat_print: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_print

!> \brief Transpose a type(matrix).
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) that should be transposed
!> \param b The transposed output type(matrix).
!>
!> Usage discouraged! If what you want is to multiply your transposed
!> matrix with something else, you should instead use mat_mul with the
!> transpose flag 'T'. This is much more efficient than transposing first 
!> and then multiplying.
!>
      SUBROUTINE mat_trans(a, b) !USAGE DISCOURAGED!!
         implicit none
         TYPE(Matrix),intent(in)     :: a
         TYPE(Matrix)                :: b !output
         REAL(REALK), ALLOCATABLE :: afull(:,:)
         
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (b%nrow /= a%ncol .or. b%ncol /= a%nrow) then
           STOP 'wrong dimensions in mat_trans'
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_trans: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_trans(a,b)
         case(mtype_dense)
             call mat_dense_trans(a,b)
         case(mtype_csr)
             call mat_csr_trans(a,b)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_trans(a,b)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
#if 0
            real(realk) :: sparsity
            ALLOCATE (afull(a%nrow,a%ncol))
            call bsm_to_full(a, afull,sparsity)
            afull = transpose(afull)
            call bsm_free(b)
            CALL bsm_init_from_full(b,a%nrow,a%ncol,afull,sparsity)
            DEALLOCATE(afull)
#else
            call bsm_transpose(a, b)
#endif
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_trans(a,b)
         case(mtype_unres_dense)
             call mat_unres_dense_trans(a,b)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_trans(a,b)
         case default
              stop "mat_trans not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('TRANS ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_trans: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_trans

!> \brief Clone a type(matrix).
!> \author B. Jansik
!> \date 2009
!> \param dest The destination type(matrix) 
!> \param src The source type(matrix) that should be cloned.
!>
!> This makes clone of matrix src to dest, similar to mat_assign, except that it is still in same
!> memory. Matrix src can then be accessed also by the dest name. All of this could be much easier, just dest=src;
!> if not for that bloody '=' operator overload with mat_assign!!
!> 
      SUBROUTINE mat_clone(dest,src)
      implicit none
      type(Matrix) :: src, dest
            dest%ncol=src%ncol; dest%nrow=src%nrow
            dest%elms=>src%elms; dest%idata=>src%idata
            dest%permutation => src%permutation
            dest%elmsb=>src%elmsb; dest%celms=>src%celms
            dest%celmsb=>src%celmsb; dest%selm1=>src%selm1
            dest%block => src%block; dest%blockpos=>src%blockpos
            dest%iaux => src%iaux; dest%raux => src%raux
            dest%complex = src%complex
            dest%val => src%val; dest%col => src%col
            dest%row => src%row; dest%nnz = src%nnz
      END SUBROUTINE mat_clone
   
!> \brief Copy a type(matrix).
!> \author L. Thogersen
!> \date 2003
!> \param a The copy output type(matrix)
!> \param b The type(matrix) that should be copied.
      SUBROUTINE mat_assign(a, b)
         implicit none
         TYPE(Matrix), INTENT(INOUT) :: a
         TYPE(Matrix), INTENT(IN)    :: b

         if (a%nrow /= b%nrow .or. a%ncol /= b%ncol) then
           STOP 'wrong dimensions in mat_assign'
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_assign: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_assign(a,b)
         case(mtype_dense)
             call mat_dense_assign(a,b)
          case(mtype_csr)
             call mat_csr_assign(a,b)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_assign(a,b)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
             call bsm_assign(a,b)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_assign(a,b)
         case(mtype_unres_dense)
             call mat_unres_dense_assign(a,b)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_assign(a,b)
         case default
              stop "mat_assign not implemented for this type of matrix"
         end select
         if (info_memory) write(mat_lu,*) 'After mat_assign: mem_allocated_global =', mem_allocated_global
       END SUBROUTINE mat_assign

#ifndef UNITTEST
!> \brief MPI broadcast a type(matrix).
!> \author T. Kjaergaard
!> \date 2010
!> \param a The type(matrix) that should be copied
!> \param slave , true if slave process 
!> \param master integer of master process
      SUBROUTINE mat_mpicopy(a, slave, master)
         implicit none
         TYPE(Matrix), INTENT(INOUT) :: a
         integer,intent(in) :: master 
         logical,intent(in) :: slave

         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_mpicopy(a,slave, master)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_mpicopy(a,slave, master)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
             call mat_bsm_mpicopy_fallback(a,slave, master)
#endif
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_mpicopy(a,slave, master)
         case default
              stop "mpicopy_typematrix not implemented for this type of matrix"
         end select
       contains
         subroutine mat_bsm_mpicopy_fallback(A,slave,master)
           use lsmpi_mod
          implicit none
          type(Matrix), intent(inout) :: A
          logical                     :: slave
          integer                     :: master
          
          real(realk), allocatable :: Afull(:,:)
          integer                  :: i, j
          
          CALL LS_MPIBCAST(A%nrow,Master)
          CALL LS_MPIBCAST(A%ncol,Master)
          allocate(Afull(A%nrow,A%ncol))
          IF(.NOT.SLAVE)call mat_to_full(A,1.0d0,Afull)
          CALL LS_MPIBCAST(Afull,A%nrow,A%ncol,Master)
          
          IF(SLAVE)THEN
             call mat_init(A,A%nrow,A%ncol)
             call mat_set_from_full(Afull,1.0d0,A) 
          ENDIF
          deallocate(Afull)
          
         end subroutine mat_bsm_mpicopy_fallback
       END SUBROUTINE mat_mpicopy
#endif

!> \brief Copy and scale a type(matrix).
!> \author L. Thogersen
!> \date 2003
!> \param alpha The scaling parameter
!> \param a The type(matrix) that should be copied
!> \param b The scaled output type(matrix).
      SUBROUTINE mat_copy(alpha,a, b) ! USAGE DISCOURAGED!
         implicit none
         REAL(REALK),  INTENT(IN)    :: alpha
         TYPE(Matrix), INTENT(IN)    :: a
         TYPE(Matrix), INTENT(INOUT) :: b

         if (b%nrow /= a%nrow .or. b%ncol /= a%ncol) then
           STOP 'wrong dimensions in mat_copy'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_copy(alpha,a,b)
         case(mtype_dense)
             call mat_dense_copy(alpha,a,b)
         case(mtype_csr)
             call mat_csr_copy(alpha,a,b)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_copy(alpha,a,b)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
             call bsm_copy(alpha,a,b)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_copy(alpha,a,b)
         case(mtype_unres_dense)
             call mat_unres_dense_copy(alpha,a,b)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_copy(alpha,a,b)
         case default
              stop "mat_copy not implemented for this type of matrix"
         end select
      END SUBROUTINE mat_copy

!> \brief Makes the trace of a square type(matrix).
!> \param a The type(matrix) we want the trace of
!> \return The trace of a
!> \author L. Thogersen
!> \date 2003
      FUNCTION mat_tr(a)
         implicit none
         TYPE(Matrix), intent(IN) :: a
         REAL(realk) :: mat_tr
#ifdef HAVE_BSM
         REAL(realk), EXTERNAL :: bsm_tr
#endif
         if (a%nrow /= a%ncol) then
           print *, 'a%nrow, a%ncol =', a%nrow, a%ncol
           STOP 'Trace is only defined for a square matrix!'
         endif
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'Before mat_tr: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             mat_Tr = mat_symm_dense_Tr(a)
         case(mtype_dense)
             mat_Tr = mat_dense_Tr(a)
         case(mtype_csr)
            mat_tr = mat_csr_Tr(a)
#ifndef UNITTEST
         case(mtype_sparse1)
             mat_Tr = mat_sparse1_Tr(a)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
             mat_Tr = bsm_tr(a)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             mat_Tr = mat_unres_symm_dense_Tr(a)
         case(mtype_unres_dense)
             mat_Tr = mat_unres_dense_Tr(a)
!         case(mtype_unres_sparse1)
!             mat_Tr = mat_unres_sparse1_Tr(a)
         case default
              stop "mat_Tr not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('TRACE ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_tr: mem_allocated_global =', mem_allocated_global
      END FUNCTION mat_tr

!> \brief Make the trace of the product of type(matrix) A and B.
!> \author L. Thogersen
!> \date 2003
!> \param a The first type(matrix) factor
!> \param b The second type(matrix) factor
!> \return Tr(a*b)
      FUNCTION mat_trAB(a,b)
         implicit none
         TYPE(Matrix), intent(IN) :: a,b
         REAL(realk) :: mat_trAB
#ifdef HAVE_BSM
         REAL(realk), EXTERNAL :: bsm_trAB
#endif
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (a%ncol /= b%nrow .or. a%nrow /= b%ncol) then
           STOP 'wrong dimensions in mat_trAB'
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_trAB: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             mat_TrAB = mat_symm_dense_TrAB(a,b)
         case(mtype_dense)
             mat_TrAB = mat_dense_TrAB(a,b)
         case(mtype_csr)
             mat_TrAB = mat_csr_TrAB(a,b)
#ifndef UNITTEST
         case(mtype_sparse1)
             mat_TrAB = mat_sparse1_TrAB(a,b)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
             mat_TrAB = bsm_trAB(a,b)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             mat_TrAB = mat_unres_symm_dense_TrAB(a,b)
         case(mtype_unres_dense)
             mat_TrAB = mat_unres_dense_TrAB(a,b)
!         case(mtype_unres_sparse1)
!             mat_TrAB = mat_unres_sparse1_TrAB(a,b)
         case default
              stop "mat_TrAB not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('TR_AB ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_trAB: mem_allocated_global =', mem_allocated_global
      END FUNCTION mat_trAB

!> \brief Traces of alpha- and beta-part of matrix-product AB
!> \author C. Nygaard
!> \date 2010-07-02
!> \param A The first matrix
!> \param B The second matrix
!> \param trace The trace of the alpha- and beta-part of the matrix-product
!=======================================================================
subroutine mat_TrAB_ab (A, B, trace)

implicit none

type(matrix), intent(in) :: A, B
real(realk), intent(out) :: trace(:)

if (size(trace) == 1) then
  trace(1) = mat_TrAB (A, B)
elseif (size(trace) == 2) then
  select case (matrix_type)
  case (mtype_dense)
    trace(1) = mat_dense_TrAB(A,B)
    trace(2) = 0.0d0
    print *, 'Warning: mat_TrAB_ab is used with mtype_dense, trace(2) = 0'
  case (mtype_unres_dense)
    call mat_unres_dense_TrAB_ab (A, B, trace)
  case default
    call lsquit ('mat_TrAB_ab not implemented for this type of matrix')
  end select
else
  call lsquit ('Wrong dimension of trace in mat_TrAB_ab')
endif


end subroutine mat_TrAB_ab
!=======================================================================

!> \brief Make c = alpha*ab + beta*c, where a,b,c are type(matrix) and alpha,beta are parameters
!> \author L. Thogersen
!> \date 2003
!> \param a The first type(matrix) factor
!> \param b The second type(matrix) factor
!> \param transa 'T'/'t' if a should be transposed, 'N'/'n' otherwise
!> \param transb 'T'/'t' if b should be transposed, 'N'/'n' otherwise
!> \param alpha The alpha parameter
!> \param beta The beta parameter
!> \param c The output type(matrix)
      SUBROUTINE mat_mul(a, b, transa, transb, alpha, beta, c)
         !c = alpha*ab + beta*c
         !transa = 'T'/'t' - transposed, 'N'/'n' - normal
         implicit none
         TYPE(Matrix), intent(IN) :: a, b
         character, intent(in)    :: transa, transb
         REAL(realk), INTENT(IN)  :: alpha, beta
         TYPE(Matrix), intent(inout):: c
         integer :: ak, bk, ci, cj
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         no_of_matmuls = no_of_matmuls + 1
         if (transa == 'n' .or. transa == 'N') then
           ak = a%ncol
           ci = a%nrow
         elseif (transa == 't' .or. transa == 'T') then
           ak = a%nrow
           ci = a%ncol
         endif
         if (transb == 'n' .or. transb == 'N') then
           bk = b%nrow
           cj = b%ncol
         elseif (transb == 't' .or. transb == 'T') then
           bk = b%ncol
           cj = b%nrow
         endif
         if (ak /= bk .or. ci /= c%nrow .or. cj /= c%ncol) then
           STOP 'wrong dimensions in mat_mul or unknown trans possibility'
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_mul: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_mul(a,b,transa, transb,alpha,beta,c)
         case(mtype_dense)
            call mat_dense_mul(a,b,transa, transb,alpha,beta,c)
#ifndef UNITTEST
         case(mtype_sparse1)
            call mat_sparse1_mul(a,b,transa, transb,alpha,beta,c)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
             call bsm_mul(a,b,transa, transb,alpha,beta,c)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_mul(a,b,transa, transb,alpha,beta,c)
         case(mtype_unres_dense)
             call mat_unres_dense_mul(a,b,transa, transb,alpha,beta,c)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_mul(a,b,transa, transb,alpha,beta,c)
         case(mtype_csr)
             call mat_csr_mul(a,b,transa, transb,alpha,beta,c)
!             call mat_csr_mul(a,b,transa,c)
         case default
              stop "mat_mul not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MATMUL ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_mul: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_mul

!> \brief Make c = alpha*a + beta*b, where a,b are type(matrix) and alpha,beta are parameters
!> \author L. Thogersen
!> \date 2003
!> \param a The first type(matrix) 
!> \param alpha The alpha parameter
!> \param b The second type(matrix) 
!> \param beta The beta parameter
!> \param c The output type(matrix)
      SUBROUTINE mat_add(alpha, a, beta, b, c)
         implicit none
         TYPE(Matrix), intent(IN) :: a, b
         REAL(realk), INTENT(IN)  :: alpha, beta
         TYPE(Matrix)             :: c

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (a%nrow /= b%nrow .or. a%ncol /= b%ncol .or. a%nrow /= c%nrow &
            &.or. a%ncol /= c%ncol) then
           print *, 'a%nrow, a%ncol, b%nrow, b%ncol, c%nrow, c%ncol', a%nrow, a%ncol, b%nrow, b%ncol, c%nrow, c%ncol
           STOP 'wrong dimensions in mat_add'
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_add: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_add(alpha,a,beta,b,c)
         case(mtype_dense)
             call mat_dense_add(alpha,a,beta,b,c)
#ifndef UNITTEST
         case(mtype_csr)
             call mat_csr_add(alpha,a,beta,b,c)
         case(mtype_sparse1)
             call mat_sparse1_add(alpha,a,beta,b,c)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
             call bsm_add(alpha,a,beta,b,c)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_add(alpha,a,beta,b,c)
         case(mtype_unres_dense)
             call mat_unres_dense_add(alpha,a,beta,b,c)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_add(alpha,a,beta,b,c)
         case default
              stop "mat_add not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('ADD   ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_add: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_add

!> \brief Make Y = alpha*X + Y where X,Y are type(matrix) and a is a parameter
!> \author L. Thogersen
!> \date 2003
!> \param alpha The alpha parameter
!> \param X The input type(matrix) 
!> \param Y The input/output type(matrix) 
      SUBROUTINE mat_daxpy(alpha, X, Y)
         implicit none
         real(realk),intent(in)       :: alpha
         TYPE(Matrix), intent(IN)     :: X
         TYPE(Matrix), intent(INOUT)  :: Y

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (x%nrow /= y%nrow .or. x%ncol /= y%ncol) then
           STOP 'wrong dimensions in mat_daxpy'
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_daxpy: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_daxpy(alpha,X,y)
         case(mtype_dense)
             call mat_dense_daxpy(alpha,x,y)
         case(mtype_csr)
             call mat_csr_daxpy(alpha,x,y)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_daxpy(alpha,x,y)
#ifdef HAVE_BSM
          case(mtype_sparse_block)
             call bsm_daxpy(alpha,x,y)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_daxpy(alpha,x,y)
         case(mtype_unres_dense)
             call mat_unres_dense_daxpy(alpha,x,y)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_daxpy(alpha,x,y)
         case default
              stop "mat_daxpy not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('DAXPY ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_daxpy: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_daxpy

!> \brief mat_daxpy if alpha is different for the two matrix parts (unres)
!> \author C. Nygaard
!> \date 2010
!> \param alpha The alpha parameters
!> \param X The input type(matrix) 
!> \param Y The input/output type(matrix)
!=======================================================================
subroutine mat_ab_daxpy (alpha, X, Y)

implicit none

real(realk), intent(in)     :: alpha(:)
type(matrix), intent(in)    :: X
type(matrix), intent(inout) :: Y

if (size(alpha) == 1) then
  call mat_daxpy (alpha(1), X, Y)
elseif (size(alpha) == 2) then
  select case (matrix_type)
  case (mtype_unres_dense)
    call mat_unres_dense_ab_daxpy (alpha, X, Y)
  case default
    stop 'mat_ab_daxpy only works for mtype_unres_dense'
  end select
else
  stop 'Wrong dimension of alpha in mat_ab_daxpy'
endif

end subroutine mat_ab_daxpy
!=======================================================================


!> \brief Make the dot product of type(matrix) a and b.
!> \author L. Thogersen
!> \date 2003
!> \param a The first type(matrix) factor
!> \param b The second type(matrix) factor
!> \return The dot product of a and b
      function mat_dotproduct(a,b)
         implicit none
         TYPE(Matrix), intent(IN) :: a,b
         REAL(realk) :: mat_dotproduct
#ifdef HAVE_BSM
         REAL(realk), EXTERNAL :: bsm_trAtransB
#endif

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (a%nrow*a%ncol /= b%nrow*b%ncol) then
           STOP 'wrong dimensions in mat_dotproduct'
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_dotproduct: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             mat_dotproduct = mat_symm_dense_dotproduct(a,b)
         case(mtype_dense)
            mat_dotproduct = mat_dense_dotproduct(a,b)
         case(mtype_csr)
            mat_dotproduct = mat_csr_dotproduct(a,b)
#ifndef UNITTEST
         case(mtype_sparse1)
             mat_dotproduct = mat_sparse1_dotproduct(a,b)
#ifdef HAVE_BSM
          case(mtype_sparse_block)
             mat_dotproduct = bsm_trAtransB(a,b)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             mat_dotproduct = mat_unres_symm_dense_dotproduct(a,b)
         case(mtype_unres_dense)
             mat_dotproduct = mat_unres_dense_dotproduct(a,b)
!         case(mtype_unres_sparse1)
!             mat_dotproduct = mat_unres_sparse1_dotproduct(a,b)
         case default
              stop "mat_dotproduct not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('DOTPRO',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_dotproduct: mem_allocated_global =', mem_allocated_global
      END FUNCTION mat_dotproduct

!> \brief Make the dot product of type(matrix) a with itself.
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) input
!> \return The dot product of a with itself
      FUNCTION mat_sqnorm2(a)
         implicit none
         TYPE(Matrix), intent(IN) :: a
         REAL(realk) :: mat_sqnorm2
#ifdef HAVE_BSM
         REAL(realk), external:: bsm_frob
#endif
         if (info_memory) write(mat_lu,*) 'Before mat_sqnorm2: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
            !         case(mtype_symm_dense)
            !             mat_sqnorm2 = mat_symm_dense_sqnorm2(a)
         case(mtype_dense)
            mat_sqnorm2 = mat_dense_sqnorm2(a)
         case(mtype_csr)
            mat_sqnorm2 = mat_csr_sqnorm2(a)
#ifndef UNITTEST
         case(mtype_sparse1)
            mat_sqnorm2 = mat_sparse1_sqnorm2(a)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            mat_sqnorm2 = bsm_frob(a)**2
#endif
#endif
!         case(mtype_unres_symm_dense)
!             mat_sqnorm2 = mat_unres_symm_dense_sqnorm2(a)
         case(mtype_unres_dense)
             mat_sqnorm2 = mat_unres_dense_sqnorm2(a)
!         case(mtype_unres_sparse1)
!             mat_sqnorm2 = mat_unres_sparse1_sqnorm2(a)
         case default
              stop "mat_sqnorm2 not implemented for this type of matrix"
         end select
         if (info_memory) write(mat_lu,*) 'After mat_sqnorm2: mem_allocated_global =', mem_allocated_global
      END FUNCTION mat_sqnorm2

!> \brief Find the absolute largest element of a type(matrix).
!> \author S. Host
!> \date 2005
!> \param a The type(matrix) input
!> \param val The absolute largest element of a
      SUBROUTINE mat_abs_max_elm(a, val) 
         implicit none
         REAL(REALK),  INTENT(OUT)   :: val
         TYPE(Matrix), INTENT(IN)    :: a
!#ifdef HAVE_BSM
!         REAL(realk), external:: bsm_max
!#endif

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'Before mat_abs_max_elm: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_abs_max_elm(a,val)
         case(mtype_dense)
             call mat_dense_abs_max_elm(a,val)
          case(mtype_csr)
             call mat_csr_abs_max_elm(a,val)
!         case(mtype_sparse1)
!            val = mat_sparse1_abs_max_elm(a)
!#ifdef HAVE_BSM
!         case(mtype_sparse_block)
!            val = bsm_abs_max(a)
!#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_abs_max_elm(a,val)
!         case(mtype_unres_dense)
!             call mat_unres_dense_abs_max_elm(a,val)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_abs_max_elm(a,val)
         case default
              stop "mat_abs_max_elm not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MAXELM',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_abs_max_elm: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_abs_max_elm

!> \brief Find the largest element of a type(matrix).
!> \author S. Host
!> \date 2005
!> \param a The type(matrix) input
!> \param val The largest element of a
      SUBROUTINE mat_max_elm(a, val) 
         implicit none
         REAL(REALK),  INTENT(OUT)   :: val
         TYPE(Matrix), INTENT(IN)    :: a
#ifdef HAVE_BSM
         REAL(realk), external:: bsm_max
#endif

         if (info_memory) write(mat_lu,*) 'Before mat_max_elm: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_max_elm(a,val)
         case(mtype_dense)
             call mat_dense_max_elm(a,val)
          case(mtype_csr)
             call mat_csr_max_elm(a,val)
#ifndef UNITTEST
         case(mtype_sparse1)
            val = mat_sparse1_max_elm(a)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            val = bsm_max(a)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_max_elm(a,val)
         case(mtype_unres_dense)
             call mat_unres_dense_max_elm(a,val)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_max_elm(a,val)
         case default
              stop "mat_max_elm not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MAXELM',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_max_elm: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_max_elm

!> \brief Find the largest element on the diagonal of a type(matrix).
!> \author S. Host
!> \date 2005
!> \param a The type(matrix) input
!> \param pos The position of the diagonal for the largest element of a
!> \param val The largest on element on the diagonal of a
      SUBROUTINE mat_max_diag_elm(a, pos, val) 
         implicit none
         TYPE(Matrix), INTENT(IN)    :: a
         INTEGER, INTENT(OUT)        :: pos
         REAL(REALK),  INTENT(OUT)   :: val

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (a%nrow /= a%ncol) then
           STOP 'matrix must be symmetric in mat_max_diag_elm'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_max_diag_elm(a,pos,val)
         case(mtype_dense)
             call mat_dense_max_diag_elm(a,pos,val)
!         case(mtype_sparse1)
!             call mat_sparse1_max_diag_elm(a,pos,val)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            CALL bsm_max_diag(a,pos,val)
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_max_diag_elm(a,pos,val)
         case(mtype_unres_dense)
             call mat_unres_dense_max_diag_elm(a,pos,val)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_max_diag_elm(a,pos,val)
         case default
              stop "mat_max_diag_elm not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MAXDIA',mat_TSTR,mat_TEN,mat_lu)
      END SUBROUTINE mat_max_diag_elm

!> \brief Squares all the off-diagonal elements in a type(matrix), and returns the sum
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) input
!> \return Sum of the squares of off-diagonal elements in a
      FUNCTION mat_outdia_sqnorm2(a)
         implicit none
         TYPE(Matrix), intent(IN) :: a
         REAL(realk) :: mat_outdia_sqnorm2
#ifdef HAVE_BSM
         REAL(realk), external:: bsm_outdia_sqnorm2
#endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             mat_outdia_sqnorm2 = mat_symm_dense_outdia_sqnorm2(a)
         case(mtype_dense)
             mat_outdia_sqnorm2 = mat_dense_outdia_sqnorm2(a)
         case(mtype_csr)
             mat_outdia_sqnorm2 = mat_csr_outdia_sqnorm2(a)
#ifndef UNITTEST
         case(mtype_sparse1)
             mat_outdia_sqnorm2 = mat_sparse1_outdia_sqnorm2(a)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
             mat_outdia_sqnorm2 = bsm_outdia_sqnorm2(a)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             mat_outdia_sqnorm2 = mat_unres_symm_dense_outdia_sqnorm2(a)
         case(mtype_unres_dense)
             mat_outdia_sqnorm2 = mat_unres_dense_outdia_sqnorm2(a)
!         case(mtype_unres_sparse1)
!             mat_outdia_sqnorm2 = mat_unres_sparse1_outdia_sqnorm2(a)
         case default
              stop "mat_outdia_sqnorm2 not implemented for this type of matrix"
         end select
!         print *, "outdia got ", mat_outdia_sqnorm2
      END FUNCTION mat_outdia_sqnorm2

!> \brief General diagonalization F*C = S*C*e
!> \author L. Thogersen
!> \date 2003
!> \param F Fock/Kohn-Sham matrix
!> \param S Overlap matrix
!> \param eival Eigenvalues
!> \param Cmo C coefficients
      SUBROUTINE mat_diag_f(F,S,eival,Cmo)
         !solves FC = SCe 
         implicit none
         TYPE(Matrix), intent(IN) :: F,S
         type(matrix)             :: Cmo  !output
         real(realk),intent(OUT)  :: eival(:)
         real(realk), allocatable :: tmp(:), eval(:), cmod(:), wrk(:)
         integer                  :: ndim
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_diag_f(F,S,eival,Cmo)
         case(mtype_dense)
             call mat_dense_diag_f(F,S,eival,Cmo)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_diag_f(F,S,eival,Cmo)
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_diag_f(F,S,eival,Cmo)
         case(mtype_unres_dense)
             call mat_unres_dense_diag_f(F,S,eival,Cmo)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_diag_f(F,S,eival,Cmo)
         case default
            print *, "FALLBACK diag_f...", S%nrow
            ndim = s%nrow
            ALLOCATE(wrk((8+2)*Ndim),tmp(Ndim*Ndim),eval(Ndim),cmod(Ndim*Ndim))
            call mat_to_full(S, 1D0, tmp)
            call mat_to_full(F, 1D0, cmod)
            call my_DSYGV(ndim,cmod,tmp,eival,wrk,size(wrk),"mat_diag_f          ")
            call mat_set_from_full(cmod, 1D0, Cmo)
            DEALLOCATE(wrk,tmp,eval,cmod)
         end select
      END SUBROUTINE mat_diag_f

!> \brief Computes dE_SCF/dmu where mu is the damping in damped roothan
!> \author L. Thogersen
!> \date 2003
!> \param Fnew New Fock/Kohn-Sham matrix
!> \param Fdamp Damped Fock/Kohn-Sham matrix
!> \param SDS S*D*S, S = overlap matrix, D = density matrix
!> \param Cmo C coefficients
!> \param nocc Number of occupied orbitals
!> \return dE_SCF/dmu
      FUNCTION mat_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
        !Find dE/dmu = sum_nu,I dE/dC_nu,I dC_nu,I/dmu
        implicit none
        type(matrix),intent(in) :: Fnew,Fdamp,SDS,Cmo
        integer, intent(in)     :: nocc
        real(realk)             :: mat_dE_dmu

         select case(matrix_type)
!         case(mtype_symm_dense)
!             mat_dE_dmu = mat_symm_dense_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
         case(mtype_dense)
             mat_dE_dmu = mat_dense_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
#ifndef UNITTEST
         case(mtype_sparse1)
             mat_dE_dmu = mat_sparse1_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
#endif
!         case(mtype_unres_symm_dense)
!             mat_dE_dmu = mat_unres_symm_dense_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
!         case(mtype_unres_dense)
!             mat_dE_dmu = mat_unres_dense_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
!         case(mtype_unres_sparse1)
!             mat_dE_dmu = mat_unres_sparse1_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
         case default
              stop "mat_dE_dmu not implemented for this type of matrix"
         end select
      END FUNCTION mat_dE_dmu
 
!> \brief Returns the sum of the elements Mat(from_row:to_row,ncol) squared
!> \author L. Thogersen
!> \date 2003
!> \param Mat Input type(matrix)
!> \param from_row Begin at this row
!> \param to_row End at this row
!> \param ncol Number of column to use
      FUNCTION mat_column_norm(Mat,ncol,from_row,to_row)
         implicit none
         type(Matrix), intent(in) :: Mat
         integer, intent(in) :: ncol, from_row, to_row
         real(realk) :: mat_column_norm

         if (to_row > Mat%nrow .or. from_row < 1 .or. ncol < 1 .or. Mat%ncol < ncol) then
           STOP 'wrong dimensions in mat_column_norm'
         endif
         if (from_row > to_row) then
           STOP 'from_row > to_row in mat_column_norm'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             mat_column_norm = mat_symm_dense_column_norm(Mat,ncol,from_row,to_row)
         case(mtype_dense)
             mat_column_norm = mat_dense_column_norm(Mat,ncol,from_row,to_row)
#ifndef UNITTEST
         case(mtype_sparse1)
             mat_column_norm = mat_sparse1_column_norm(Mat,ncol,from_row,to_row)
#endif
!         case(mtype_unres_symm_dense)
!             mat_column_norm = mat_unres_symm_dense_column_norm(Mat,ncol,from_row,to_row)
         case(mtype_unres_dense)
             mat_column_norm = mat_unres_dense_column_norm(Mat,ncol,from_row,to_row)
!         case(mtype_unres_sparse1)
!             mat_column_norm = mat_unres_sparse1_column_norm(Mat,ncol,from_row,to_row)
         case default
              stop "mat_column_norm not implemented for this type of matrix"
         end select
      END FUNCTION mat_column_norm

!> \brief Returns a section of a matrix
!> \author L. Thogersen
!> \date 2003
!> \param A Input type(matrix)
!> \param from_row Begin at this row
!> \param to_row End at this row
!> \param from_col Begin at this column
!> \param to_col End at this column
!> \param Asec The section of the type(matrix)
      subroutine mat_section(A,from_row,to_row,from_col,to_col,Asec)
         implicit none
         type(Matrix), intent(in) :: A
         integer, intent(in) :: from_row, to_row, from_col, to_col
         type(Matrix), intent(inout) :: Asec  !output

         !Check if Asec is inside A
         if (to_row > A%nrow .or. from_row < 1 .or. from_col < 1 .or. A%ncol < to_col) then
           STOP 'Asec not inside A in mat_section'
         endif
         !Check if the section size is positive
         if (from_row > to_row .or. from_col > to_col) then
           STOP 'from_row or from_col > to_row or to_col in mat_section'
         endif
         !Check if allocated space for section is the right size
         if (Asec%nrow /= to_row - from_row + 1 .or.&
            & Asec%ncol /= to_col - from_col + 1) then
            STOP 'Wrong dimensions in mat_section'
         endif

         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_section(A,from_row,to_row,from_col,to_col,Asec)
         case(mtype_dense)
             call mat_dense_section(A,from_row,to_row,from_col,to_col,Asec)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_section(A,from_row,to_row,from_col,to_col,Asec)
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_section(A,from_row,to_row,from_col,to_col,Asec)
         case(mtype_unres_dense)
             call mat_unres_dense_section(A,from_row,to_row,from_col,to_col,Asec)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_section(A,from_row,to_row,from_col,to_col,Asec)
         case default
              stop "mat_section not implemented for this type of matrix"
         end select
      END SUBROUTINE mat_section

!> \brief Puts smaller matrix A in selected position in larger matrix B
!> \author L. Thogersen
!> \date 2003
!> \param A The input type(matrix)
!> \param from_row Begin at this row in B
!> \param to_row End at this row in B
!> \param from_col Begin at this column in B
!> \param to_col End at this column in B
!> \param B The output type(matrix)
      subroutine mat_section2(A,from_row,to_row,from_col,to_col,B)
         implicit none
         type(Matrix), intent(in) :: A
         integer, intent(in) :: from_row, to_row, from_col, to_col
         type(Matrix), intent(inout) :: B  !output

         if (to_row > B%nrow .or. from_row < 1 .or. from_col < 1 .or. B%ncol < to_col) then
           STOP 'wrong dimensions in mat_section2'
         endif
         if (from_row > to_row .or. from_col > to_col) then
           STOP 'from_row or from_col > to_row or to_col in mat_section2'
         endif
         if ((to_row - from_row + 1) /= A%nrow .or. (to_col - from_col + 1) /= A%ncol) then
           STOP 'wrong dimensions in mat_section2'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_section(A,from_row,to_row,from_col,to_col,Asec)
         case(mtype_dense)
             call mat_dense_section2(A,from_row,to_row,from_col,to_col,B)
!         case(mtype_sparse1)
!             call mat_sparse1_section(A,from_row,to_row,from_col,to_col,Asec)
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_section(A,from_row,to_row,from_col,to_col,Asec)
!         case(mtype_unres_dense)
!             call mat_unres_dense_section(A,from_row,to_row,from_col,to_col,Asec)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_section(A,from_row,to_row,from_col,to_col,Asec)
         case default
              stop "mat_section2 not implemented for this type of matrix"
         end select
      END SUBROUTINE mat_section2

!> \brief Mix two matrices into a plus and a minus combination
!> \author C. Nygaard
!> \date 2010
!> \param homo One of the matrices (+ in both combinations)
!> \param lumo The other matrix (different signs in the combinations)
subroutine mat_mix_homolumo (homo, lumo)

!Made for the unrestricted part,
! removes the symmetry of the starting guess
!Mix homo=(homo^alpha , homo^beta) and lumo=(lumo^alpha , lumo^beta)
! into newhomo=(plus^alpha , minus^beta)
!  and newlumo=(minus^alpha , plus^beta)
! where plus=homo+lumo and minus=homo-lumo

implicit none

type(matrix), intent(inout) :: homo, lumo
type(matrix)                :: plus, minus
integer                     :: r, c
real(realk)                 :: Na, Nb, ya, yb

r = homo%nrow ; c = homo%ncol

if (lumo%nrow /= r .or. lumo%ncol /= c) then
  write (mat_lu, *) 'Different dimensions of homo and lumo in mat_mix_homolumo!'
  call lsquit ('Wrong dimensions in mat_mix_homolumo',mat_lu)
endif

call mat_init (plus, r, c)
call mat_init (minus, r, c)

ya = 0.1d0 ; yb = 0.05d0
Na = 1.0d0/DSQRT((1.0d0-ya)**2 + ya**2)
Nb = 1.0d0/DSQRT((1.0d0-yb)**2 + yb**2)

call mat_add (Na*(1.0d0-ya), homo, Na*ya, lumo, plus)
call mat_add (Nb*(1.0d0-yb), homo, -Nb*yb, lumo, minus)

select case (matrix_type)
  case (mtype_unres_dense)
    call mat_unres_dense_mix_homolumo (plus, minus, homo, lumo)
  case default
    stop 'mat_mix_homolumo not implemented for this type of matrix'
end select

call mat_free (plus)
call mat_free (minus)

end subroutine mat_mix_homolumo

!> \brief Preconditioning of vector x by matrix M: xprec(i) = x(i)/M(i,i)
!> \author S. Host
!> \date 2005
!> \param M The preconditioner
!> \param x Input vector to be preconditioned
!> \param xprec Preconditioned output vector
      subroutine mat_precond(M,x,xprec)
         implicit none
         type(Matrix), intent(inout) :: xprec
         type(Matrix), intent(in)    :: M, x

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_precond(M,x,xprec)
         case(mtype_dense)
             call mat_dense_precond(M,x,xprec)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            call bsm_precond(M,x,xprec)
#endif /* HAVE_BSM */
!         case(mtype_sparse1)
!             call mat_sparse1_precond(M,x,xprec)
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_precond(M,x,xprec)
         case(mtype_unres_dense)
             call mat_unres_dense_precond(M,x,xprec)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_precond(M,x,xprec)
         case default
            stop "mat_precond not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MATPRE',mat_TSTR,mat_TEN,mat_lu)
      END SUBROUTINE mat_precond

!> \brief Divide the occ-virt part of the matrix X_MO with the orbital energy difference E_a - E_i
!> \author L. Thogersen
!> \date 2005
!> \param nocc Number of occupied orbitals
!> \param omega Levelshift
!> \param Eorb_final Final orbital energies
!> \param X_MO Input/output - matrix to be preconditioned
      subroutine mat_mo_precond(nocc,omega,Eorb_final,X_MO)
         implicit none
         integer, intent(in) :: nocc
         real(realk), intent(in) :: omega
         real(realk), intent(in) :: Eorb_final(:)
         type(Matrix), intent(inout) :: X_MO
         integer :: i, j
         real(realk), ALLOCATABLE :: XMO_full(:,:)
         real(realk) :: dia

         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
         case(mtype_dense)
             call mat_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_mo_precond(nocc,omega,Eorb_final,X_MO)
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
         case(mtype_unres_dense)
             call mat_unres_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_mo_precond(nocc,omega,Eorb_final,X_MO)
         case default
            print *, "FALLBACK mo_precond"
            ! see mat_dense_mo_precond
            allocate(xmo_full(X_MO%nrow, X_MO%ncol))
            call mat_to_full(X_MO,1D0,xmo_full)
            do j = 1,nocc  !columns
               do i = nocc+1,X_MO%nrow  !rows
                  dia = 2.0d0*(Eorb_final(i) - Eorb_final(j)) + omega * 2.d0
                  xmo_full(i,j) = xmo_full(i,j) /dia
               enddo
            enddo
            !modify upper right block (ia block) 
            do j = nocc+1,X_MO%ncol  !columns
               do i = 1,nocc          !rows
                  !dia = E[2]_dia + omega*S[2]_dia
                  dia = 2.0d0*(Eorb_final(j) - Eorb_final(i)) - omega * 2.d0
                  xmo_full(i,j) = xmo_full(i,j)/dia
               enddo
            enddo
            call mat_set_from_full(xmo_full, 1D0, X_MO)
            DEALLOCATE(XMO_full)
         end select
      END SUBROUTINE mat_mo_precond

!vide the occ-virt part of the matrix X_MO with the orbital energy difference E_a - E_i
 !> \author J. Kauczor
 !> \date 2010
 !> \param nocc Number of occupied orbitals
 !> \param omega Levelshift
 !> \param Eorb_final Final orbital energies
 !> \param X_MO Input/output - matrix to be preconditioned
       subroutine mat_new_mo_precond(nocc,omega,Eorb_final,Xp_MO,xm_mo)
          implicit none
          integer, intent(in) :: nocc
          real(realk), intent(in) :: omega
          real(realk), intent(in) :: Eorb_final(:)
          type(Matrix), intent(inout) :: Xp_MO,xm_mo
          integer :: i, j,ndim
          real(realk), ALLOCATABLE :: XpMO_full(:,:),XmMO_full(:,:)
          real(realk) :: dia,aa
 
          ndim=xp_mo%nrow
          select case(matrix_type)
          case(mtype_dense)
              call mat_dense_new_mo_precond(nocc,omega,Eorb_final,Xp_MO,xm_mo)
          case default
             print *, "FALLBACK mo_precond"
             ! see mat_dense_mo_precond
             allocate(xpmo_full(Xp_MO%nrow, Xp_MO%ncol))
             allocate(xmmo_full(Xm_MO%nrow, Xm_MO%ncol))
             call mat_to_full(Xp_MO,1D0,xpmo_full)
             call mat_to_full(Xm_MO,1D0,xmmo_full)
             do j = 1,nocc  !columns
                do i = nocc+1,Xp_MO%nrow  !rows
                   aa=2d0*(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)- Eorb_final(j))-2d0*omega*omega
                 !  dia = 2.0d0*(Eorb_final(i) - Eorb_final(j)) + omega * 2.d0
                   dia=xpmo_full(i,j)
                   xpmo_full(i,j) = ((Eorb_final(i) - Eorb_final(j))*xpmo_full(i,j)-omega*xmmo_full(i,j)) /aa
                   xmmo_full(i,j) = ((Eorb_final(i) - Eorb_final(j))*xmmo_full(i,j)-omega*dia) /aa
                enddo
             enddo
             !modify upper right block (ia block) 
             do j = nocc+1,Xp_MO%ncol  !columns
                do i = 1,nocc          !rows
                   !dia = E[2]_dia + omega*S[2]_dia
                   aa=2d0*(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)-Eorb_final(j))-2d0*omega*omega
                   !dia = 2.0d0*(Eorb_final(j) - Eorb_final(i)) - omega * 2.d0
                   dia=xpmo_full(i,j)
                   xpmo_full(i,j) = ((Eorb_final(j) - Eorb_final(i))*xpmo_full(i,j)+omega*xmmo_full(i,j)) /aa
                   xmmo_full(i,j) = ((Eorb_final(j) - Eorb_final(i))*xmmo_full(i,j)+omega*dia) /aa
                  ! xmo_full(i,j) = xmo_full(i,j)/dia
                enddo
             enddo
             call mat_set_from_full(xpmo_full, 1D0, Xp_MO)
             call mat_set_from_full(xmmo_full, 1D0, Xm_MO)
             DEALLOCATE(XpMO_full)
             DEALLOCATE(XmMO_full)
          end select
       END SUBROUTINE mat_new_mo_precond
 
 
 !> \brief Divide the occ-virt part of the matrix X_MO with the orbital energy difference E_a - E_i (complex orbitals)
 !> \author J. Kauczor
 !> \date 2010
 !> \param nocc Number of occupied orbitals
 !> \param Eorb_final Final orbital energies
 !> \param X_MO Input/output - matrix to be preconditioned
       subroutine mat_new_complex_precond(nocc,omega,gammma,Eorb_final,Xp_MO,xm_mo,xpi_mo,xmi_mo)
          implicit none
          integer, intent(in) :: nocc
          real(realk), intent(in) :: omega,gammma
          real(realk), intent(in) :: Eorb_final(:)
          type(Matrix), intent(inout) :: Xp_MO,xm_mo,xpi_mo,xmi_mo
          integer :: i, j,ndim
          real(realk), ALLOCATABLE :: XpMO_full(:,:),XmMO_full(:,:)
          real(realk), ALLOCATABLE :: XpiMO_full(:,:),XmiMO_full(:,:)
          real(realk), ALLOCATABLE :: Xp(:,:),Xm(:,:)
          real(realk), ALLOCATABLE :: Xpi(:,:),Xmi(:,:)
          real(realk) :: dia,aa,ab,a,b,c,d
 
          ndim=xp_mo%nrow
          select case(matrix_type)
          case(mtype_dense)
              call mat_dense_new_complex_precond(nocc,omega,gammma,Eorb_final,Xp_MO,xm_mo,xpi_mo,xmi_mo)
          case default
             print *, "FALLBACK mo_precond"
             ! see mat_dense_mo_precond
             allocate(xpmo_full(Xp_MO%nrow, Xp_MO%ncol))
             allocate(xmmo_full(Xm_MO%nrow, Xm_MO%ncol))
             allocate(xpimo_full(Xpi_MO%nrow, Xpi_MO%ncol))
             allocate(xmimo_full(Xmi_MO%nrow, Xmi_MO%ncol))
             allocate(xp(Xp_MO%nrow, Xp_MO%ncol))
             allocate(xm(Xm_MO%nrow, Xm_MO%ncol))
             allocate(xpi(Xpi_MO%nrow, Xpi_MO%ncol))
             allocate(xmi(Xmi_MO%nrow, Xmi_MO%ncol))
             xp=0d0; xm=0d0; xpi=0d0; xmi=0d0
             call mat_to_full(Xp_MO,1D0,xpmo_full)
             call mat_to_full(Xm_MO,1D0,xmmo_full)
             call mat_to_full(Xpi_MO,1D0,xpimo_full)
             call mat_to_full(Xmi_MO,1D0,xmimo_full)
             do i=1,Xp_MO%nrow
                xp(i,i)=xpmo_full(i,i)
                xpi(i,i)=xpimo_full(i,i)
               ! xm(i,i)=xpmo_full(i,i)
               ! xmi(i,i)=xpimo_full(i,i)
             enddo
             do j = 1,nocc  !columns
                do i = nocc+1,Xp_MO%nrow  !rows
                   ab=(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)-Eorb_final(j))-(omega*omega-gammma*gammma) 
                   aa=2d0*ab*ab+8d0*(omega*omega*gammma*gammma)
                 !  dia=xpmo_full(i,j)
                   a=(Eorb_final(i) - Eorb_final(j))*ab
                   b=-omega*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
                      &-(omega*omega+gammma*gammma))
                   c=-gammma*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
                      &+(omega*omega+gammma*gammma))
                   d=2d0*omega*gammma*(Eorb_final(i) - Eorb_final(j))
                   
                   xp(i,j)  = (a*xpmo_full(i,j)+b*xmmo_full(i,j)-d*xpimo_full(i,j)-c*xmimo_full(i,j)) /aa
                   xm(i,j)  = (b*xpmo_full(i,j)+a*xmmo_full(i,j)-c*xpimo_full(i,j)-d*xmimo_full(i,j))/aa
                   xpi(i,j) = (d*xpmo_full(i,j)+c*xmmo_full(i,j)+a*xpimo_full(i,j)+b*xmimo_full(i,j))/aa
                   xmi(i,j) = (c*xpmo_full(i,j)+d*xmmo_full(i,j)+b*xpimo_full(i,j)+a*xmimo_full(i,j))/aa
                enddo
             enddo
             !modify upper right block (ia block) 
             do j = nocc+1,Xp_MO%ncol  !columns
                do i = 1,nocc          !rows
        ab=(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)-Eorb_final(j))-(omega*omega-gammma*gammma) 
         aa=2d0*ab*ab+8d0*(omega*omega*gammma*gammma)
         
         a=(Eorb_final(j) - Eorb_final(i))*ab
         b=omega*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
         &-(omega*omega+gammma*gammma))
         c=gammma*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
         &+(omega*omega+gammma*gammma))
         d=2d0*omega*gammma*(Eorb_final(j) - Eorb_final(i))
                   
                   xp(i,j)  = (a*xpmo_full(i,j)+b*xmmo_full(i,j)-d*xpimo_full(i,j)-c*xmimo_full(i,j)) /aa
                   xm(i,j)  = (b*xpmo_full(i,j)+a*xmmo_full(i,j)-c*xpimo_full(i,j)-d*xmimo_full(i,j))/aa
                   xpi(i,j) = (d*xpmo_full(i,j)+c*xmmo_full(i,j)+a*xpimo_full(i,j)+b*xmimo_full(i,j))/aa
                   xmi(i,j) = (c*xpmo_full(i,j)+d*xmmo_full(i,j)+b*xpimo_full(i,j)+a*xmimo_full(i,j))/aa
                enddo
             enddo
     
             call mat_set_from_full(xp, 1D0, Xp_MO)
             call mat_set_from_full(xm, 1D0, Xm_MO)
             call mat_set_from_full(xpi, 1D0, Xpi_MO)
             call mat_set_from_full(xmi, 1D0, Xmi_MO)
             DEALLOCATE(XpMO_full)
             DEALLOCATE(XmMO_full)
             DEALLOCATE(XpiMO_full)
             DEALLOCATE(XmiMO_full)
             DEALLOCATE(Xp)
             DEALLOCATE(Xm)
             DEALLOCATE(Xpi)
             DEALLOCATE(Xmi)
          end select
       END SUBROUTINE mat_new_complex_precond

!> \author J. Kauczor
!> \date 2009
!> \param nocc Number of occupied orbitals
!> \param Eorb_final Final orbital energies
!> \param X_MO Input/output - matrix to be preconditioned
      subroutine mat_mo_precond_complex(nocc,Eorb_final,X_MO)
         implicit none
         integer, intent(in) :: nocc
         real(realk), intent(in) :: Eorb_final(:)
         type(Matrix), intent(inout) :: X_MO
         integer :: i, j,ndim
         real(realk), ALLOCATABLE :: XMO_full(:,:)
         real(realk) :: dia

         ndim=X_MO%nrow
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_mo_precond_complex(nocc,Eorb_final,X_MO)
         case default
            print *, "FALLBACK mo_precond"
            ! see mat_dense_mo_precond
            allocate(xmo_full(X_MO%nrow, X_MO%ncol))
            call mat_to_full(X_MO,1D0,xmo_full)
            do j = 1,nocc  !columns
               do i = nocc+1,X_MO%nrow  !rows
                  dia = 2.0d0*(Eorb_final(i) - Eorb_final(j))
                  xmo_full(i,j) = xmo_full(i,j) /dia
               enddo
            enddo
            !modify upper right block (ia block) 
            do j = nocc+1,X_MO%ncol  !columns
               do i = 1,nocc          !rows
                  !dia = E[2]_dia + omega*S[2]_dia
                  dia = 2.0d0*(Eorb_final(j) - Eorb_final(i))
                  xmo_full(i,j) = xmo_full(i,j)/dia
               enddo
            enddo
            call mat_set_from_full(xmo_full, 1D0, X_MO)
            DEALLOCATE(XMO_full)
         end select
      END SUBROUTINE mat_mo_precond_complex

!> \brief Diagonal preconditioning in orthonormal AO basis 
!> \author S. Host
!> \date 2005
!> \param symm Symmetry indicator: symmetric = 1, antisymmetric = 2, nonsymmetric = 0
!> \param omega Level shift
!> \param FUP Fock/KS matrix in OAO basis, occupied part (virtual part projected out)
!> \param FUQ Fock/KS matrix in OAO basis, virtual part (occupied part projected out)
!> \param DU Density matrix in OAO basis
!> \param X_AO Matrix to be preconditioned
!> 
!> Preconditioning with orbital energy difference. This is done by using the occupied and virtual 
!> parts of the Fock/KS matrix: 
!> X_prec(i,j) = X(i,j) / [FUQ(j,j)-FUP(j,j) + FUQ(i,i)-FUP(i,i) - omega*(DU(j,j)-DU(i,i))]
!> taking care not to divide by zero and exploiting symmetry if X is symm or antisymm
!> 
      subroutine mat_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
         implicit none
         integer, intent(in) :: symm
         real(realk), intent(in) :: omega
         type(Matrix), intent(in) :: FUP, FUQ,DU
         type(Matrix), intent(inout) :: X_AO

         if (info_memory) write(mat_lu,*) 'Before mat_ao_precond: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
         case(mtype_dense)
            call mat_dense_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
         case(mtype_csr)
            call mat_csr_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
         case(mtype_unres_dense)
            call mat_unres_dense_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
#if 1
            PRINT *, "FALLBACK: mat_ao_bsm_precond half-optimal", symm
            call mat_ao_precond_bsm_fallback
#else
            call bsm_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
#endif 
#endif /* HAVE_BSM */
         case default
            print *, "FALLBACK: mat_ao_precond"
            call mat_ao_precond_fallback
         end select
         if (info_memory) write(mat_lu,*) 'After mat_ao_precond: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('AOPREC',mat_TSTR,mat_TEN,mat_lu)
       contains
         subroutine mat_ao_precond_fallback
           implicit none
           real(realk), allocatable :: FPPd(:,:), FQQd(:,:)
           real(realk), allocatable :: Dd(:,:)
           real(realk), allocatable :: X(:,:)
           integer :: n, i, j
           real(realk) :: denom

           n = FUP%nrow
           allocate(FPPd(n,n),FQQd(n,n),Dd(n,n),X(n,n))
           call mat_to_full(FUP,1d0,FPPd)
           call mat_to_full(FUQ,1d0,FQQd)
           call mat_to_full(DU,1d0,Dd)
           !call mat_to_full(SQQ,1d0,SQQd)
           call mat_to_full(X_AO,1d0,X)
              ! extract diagonals and use them instead of inefficient
              ! stride (n+1) access!
           if(symm .eq. 1 .or. symm .eq. 2) THEN
              DO i = 1, n
                 DO j = 1, n
                    denom = FQQd(j,j) + FQQd(i,i) &
                         &- FPPd(i,i) - FPPd(j,j) &
                         &- omega
                    IF(ABS(denom)>1d-9) X(i,j) = X(i,j)/(denom)  !12/10-09
                                                                 !Stinne removed factor 2 to match matop_dense
                 END DO
              END DO
           ELSE
              DO i = 1, n
                 DO j = 1, n
                    denom = FQQd(j,j) + FQQd(i,i) &
                         &- FPPd(i,i) - FPPd(j,j) &
                         &- omega*(Dd(i,i)-Dd(j,j)) 
                    IF(ABS(denom)>1d-9) X(i,j) = X(i,j)/(denom) !12/10-09
                                                                 !Stinne removed factor 2 to match matop_dense 
                 END DO
              END DO
           end if
           call mat_set_from_full(X,1d0, X_AO)
           deallocate(FPPd,FQQd,Dd,X)
         end subroutine mat_ao_precond_fallback
#ifdef HAVE_BSM
         subroutine mat_ao_precond_bsm_fallback
            implicit none
            real(realk), allocatable :: FPPdiag(:), FQQdiag(:)
            REAL(realk), ALLOCATABLE :: Ddiag(:), X(:,:)
            integer :: n, i, j
            real(realk) :: denom

            n = FUP%nrow
            allocate(FPPdiag(n),FQQdiag(n),X(n,n))
            ! extract diagonals and use them instead of inefficient
            ! stride (n+1) access!
            call bsm_extract_diag(FUP,FPPdiag)
            call bsm_extract_diag(FUQ,FQQdiag)
            call mat_to_full(X_AO,1d0,X)
            IF(symm .EQ. 1 .OR. symm .EQ. 2) THEN
               DO i = 1, n
                  DO j = 1, n
                     denom = FQQdiag(j) + FQQdiag(i) &
                          &- FPPdiag(i) - FPPdiag(j) &
                          &- omega
                     IF(ABS(denom)>1d-9) X(i,j) = X(i,j)/(2d0*denom)
                  END DO
               END DO
            ELSE
               allocate(Ddiag(n))
               CALL bsm_extract_diag(DU,Ddiag)
               DO i = 1, n
                  DO j = 1, n
                     denom = FQQdiag(j) + FQQdiag(i) &
                          &- FPPdiag(i) - FPPdiag(j) &
                          &- omega*(Ddiag(i)-Ddiag(j)) 
                     IF(ABS(denom)>1d-9) X(i,j) = X(i,j)/(2d0*denom)
                  END DO
               END DO
              DEALLOCATE(Ddiag)
            END IF
            call mat_set_from_full(X,1d0, X_AO)
            deallocate(FPPdiag,FQQdiag,X)
         end subroutine mat_ao_precond_bsm_fallback
#endif /* HAVE_BSM */
      END SUBROUTINE mat_ao_precond

!> \brief Set a type(matrix) to identity, i.e. I(i,j) = 1 for i = j, 0 otherwise
!> \author L. Thogersen
!> \date 2003
!> \param I Matrix to be set equal to identity
      subroutine mat_identity(I)
         implicit none
         type(Matrix), intent(inout) :: I
         real(realk), ALLOCATABLE    :: ifull(:,:)
         integer                     :: j

         if (info_memory) write(mat_lu,*) 'Before mat_identity: mem_allocated_global =', mem_allocated_global
         !print *, "mat_identity inefficient, use mat_add_identity instead!"
         if (I%nrow /= I%ncol) then
           STOP 'cannot make identity matrix with different ncol and nrow'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_identity(I)
         case(mtype_dense)
             call mat_dense_identity(I)
         case(mtype_csr)
            call mat_csr_identity(I)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_identity(I)
#if defined(HAVE_BSM)
         case(mtype_sparse_block)
            call bsm_identity(I)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_identity(I)
         case(mtype_unres_dense)
             call mat_unres_dense_identity(I)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_identity(I)
         case default
            print *, "FALLBACK: mat_identity"
            allocate(ifull(I%nrow, I%ncol))
            DO j = 1, I%ncol
               ifull(1:j-1,j) = 0D0
               ifull(j,j)     = 1D0
               ifull(j+1:I%nrow,j) = 0D0
            END DO
            call mat_set_from_full(ifull,1D0,I)
            deallocate(ifull)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_identity: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_identity

!> \brief Add identity to a type(matrix), i.e. C = alpha*I + beta*B  >>> NOTE: ALLOCATES A MATRIX! <<<
!> \author L. Thogersen
!> \date 2003
!> \param alpha Alpha parameter
!> \param beta Beta parameter
!> \param B Input matrix B
!> \param C Output matrix C
      SUBROUTINE mat_add_identity(alpha, beta, B, C)
         implicit none
         TYPE(Matrix), intent(IN) :: B
         REAL(realk), INTENT(IN)  :: alpha, beta
         TYPE(Matrix)             :: C
         type(matrix)             :: I

         if (info_memory) write(mat_lu,*) 'Before mat_add_identity: mem_allocated_global =', mem_allocated_global
         if (b%nrow /= c%nrow .or. b%ncol /= c%ncol) then
           STOP 'wrong dimensions in mat_add_identity'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_add(alpha,a,beta,b,c)
         case(mtype_dense)
            call mat_init(I, b%nrow, b%ncol)
            call mat_dense_identity(I)
            call mat_dense_add(alpha,I,beta,b,c)
            call mat_free(I)
         case(mtype_csr)
            call mat_csr_add_identity(alpha, beta, B, C)
#ifndef UNITTEST
         case(mtype_sparse1)
            call mat_init(I, b%nrow, b%ncol)
            call mat_sparse1_identity(I)
            call mat_sparse1_add(alpha,I,beta,b,c)
            call mat_free(I)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            c = B
            call mat_scal(beta, c)
            call bsm_add_identity(c, alpha)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_add(alpha,a,beta,b,c)
         case(mtype_unres_dense)
            call mat_init(I, b%nrow, b%ncol)
            call mat_unres_dense_identity(I)
            call mat_unres_dense_add(alpha,I,beta,b,c)
            call mat_free(I)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_add(alpha,a,beta,b,c)
         case default
              stop "mat_add_identity not implemented for this type of matrix"
         end select
         if (info_memory) write(mat_lu,*) 'After mat_add_identity: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_add_identity

!> \brief Create or overwrite element in matrix
!> \author S. Host
!> \date 2005
!> \param i Row index for element to be created
!> \param j Column index for element to be created
!> \param val Value of element to be created
!> \param A Input/output matrix
      subroutine mat_create_elm(i,j,val,A)
         implicit none
         integer, intent(in) :: i,j
         real(Realk), intent(in) :: val
         type(Matrix), intent(inout) :: A

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (i > A%nrow .or. j > A%ncol .or. i < 0 .or. j < 0) then
           WRITE(mat_lu,*) 'cannot create element, the indexes',i,j,&
                        & 'are out of the bounds nrow,ncol ',A%nrow,A%ncol
           STOP 'cannot create element, the indexes are out of bounds'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_create_elm(i,j,val,A)
         case(mtype_dense)
             call mat_dense_create_elm(i,j,val,A)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_create_elm(i,j,val,A)
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_create_elm(i,j,val,A)
         case(mtype_unres_dense)
             call mat_unres_dense_create_elm(i,j,val,A)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_create_elm(i,j,val,A)
         case default
              stop "mat_create_elm not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('CREATE',mat_TSTR,mat_TEN,mat_lu)
      END SUBROUTINE mat_create_elm

!> \brief mat_create_elm for different alpha and beta parts
!> \author C. Nygaard
!> \date 2010
!> \param r Row index for element to be created
!> \param c Column index for element to be created
!> \param elm Values of elements to be created
!> \param A Input/output matrix
!=======================================================================
  subroutine mat_create_ab_elms (r, c, elm, A)

  implicit none

  integer, intent(in)         :: r, c
  real(realk), intent(in)     :: elm(:)
  type(matrix), intent(inout) :: A

  if (size(elm) == 1) then
    call mat_create_elm (r, c, elm(1), A)
  elseif (size(elm) == 2) then
    select case (matrix_type)
    case (mtype_unres_dense)
      call mat_unres_dense_create_ab_elms (r, c, elm, A)
    case (mtype_dense)
      call mat_dense_create_elm (r, c, elm(1), A)
    case default
      stop 'mat_create_ab_elm only works for unrestricted, you want&
           & mat_create_elm'
    end select
  else
    print *, 'expected dimension = 1 or 2, actual dimension =', size(elm)
    stop 'wrong dimension of elm in mat_create_ab_elms'
  endif

  end subroutine mat_create_ab_elms
!=======================================================================

!> \brief Get element from matrix
!> \author C. Nygaard
!> \date 2010
!> \param A Input matrix
!> \param r Row index for element
!> \param c Column index for element
!> \param elm Value of element (output)
!=======================================================================
      subroutine mat_get_elm (A, r, c, elm)

      implicit none

      type(matrix), intent(in) :: A
      integer, intent(in)      :: r, c
      real(realk), intent(out) :: elm
      real(realk)              :: tmp(2)

      select case (matrix_type)
      case (mtype_dense)
        call mat_dense_get_elm (A, r, c, elm)
      case (mtype_unres_dense)
        call mat_unres_dense_get_elm (A, r, c, tmp)
        elm = tmp(1)
      case default
        stop "mat_get_elm not implemented for this matrix-type"
      end select

      end subroutine mat_get_elm
!=======================================================================

!> \brief mat_get_elm for different alpha and beta parts
!> \author C. Nygaard
!> \date 2010
!> \param A Input matrix
!> \param r Row index for element
!> \param c Column index for element
!> \param elm Value of elements (output)
!=======================================================================
subroutine mat_get_ab_elms (A, r, c, elm)

implicit none

type(matrix), intent(in) :: A
integer, intent(in)      :: r, c
real(realk), intent(out) :: elm(:)

if (size(elm) == 1) then
  call mat_get_elm (A, r, c, elm(1))
elseif (size(elm) == 2) then
  select case (matrix_type)
  case (mtype_dense)
    call mat_get_elm (A, r, c, elm(1))
    elm(2) = 0.0d0
  case (mtype_unres_dense)
    call mat_unres_dense_get_elm (A, r, c, elm)
  case default
    stop 'mat_get_ab_elm is only implemented for mtype_unres_dense'
  end select
else
  stop 'Wrong dimensions of elm in mat_get_ab_elms!'
endif

end subroutine mat_get_ab_elms
!=======================================================================

!> \brief Create or overwrite block in a type(matrix)
!> \author S. Host
!> \date 2009
!> \param A Input/output matrix where we want to create a block
!> \param fullmat Standard fortran matrix containing the block to put into A
!> \param fullrow Number of rows in fullmat
!> \param fullcol Number of columns in fullmat
!> \param insertrow Insert block in A beginning at this row
!> \param insertcol Insert block in A beginning at this col
      subroutine mat_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         implicit none
         integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
         real(Realk), intent(in) :: fullmat(fullrow,fullcol)
         type(Matrix), intent(inout) :: A

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (insertrow+fullrow > A%nrow+1 .or. &
          &  insertcol+fullcol > A%ncol+1 .or. fullrow < 1 .or. fullcol < 1 .or. &
          & insertrow < 1 .or. insertcol < 1) then
           WRITE(mat_lu,*) 'Cannot create block, the indexes', &
           & fullrow,fullcol,insertrow,insertcol, &
           & 'are out of the bounds - nrow, ncol =',A%nrow,A%ncol
           CALL lsQUIT('Cannot create block (subroutine mat_create_block)',mat_lu)
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_dense)
             call mat_dense_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            print *, "FALLBACK: mat_create_block converts to full for Block Sparse Matrices"
            call mat_create_block_bsm_fallback(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_unres_dense)
             call mat_unres_dense_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case default
              stop "mat_create_block not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('CREATE',mat_TSTR,mat_TEN,mat_lu)

       contains
         subroutine mat_create_block_bsm_fallback(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         implicit none
         integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
         real(Realk), intent(in) :: fullmat(fullrow,fullcol)
         type(Matrix), intent(inout) :: A
         real(realk), allocatable :: Afull(:,:)
         integer                  :: i, j

stop 'inside mat_create_block_bsm_fallback'
         allocate(Afull(A%nrow,A%ncol))
         call mat_to_full(A,1.0d0,Afull)

         do i = insertrow, insertrow+fullrow-1
            do j = insertcol, insertcol+fullcol-1
               Afull(i,j) = fullmat(i-insertrow+1,j-insertcol+1)
            enddo
         enddo

         call mat_set_from_full(Afull,1.0d0,A) 
         deallocate(Afull)

         end subroutine mat_create_block_bsm_fallback

      END SUBROUTINE mat_create_block

!> \brief Add block to type(matrix) - add to existing elements, don't overwrite
!> \author T. Kjaergaard
!> \date 2009
!> \param A Input/output matrix where we want to add a block
!> \param fullmat Standard fortran matrix containing the block to add to A
!> \param fullrow Number of rows in fullmat
!> \param fullcol Number of columns in fullmat
!> \param insertrow Add block to A beginning at this row
!> \param insertcol Add block to A beginning at this col
      subroutine mat_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         implicit none
         integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
         real(Realk), intent(inout) :: fullmat(fullrow,fullcol)
         type(Matrix), intent(inout) :: A

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (insertrow+fullrow > A%nrow+1 .or. &
          &  insertcol+fullcol > A%ncol+1 .or. fullrow < 1 .or. fullcol < 1 .or. &
          & insertrow < 1 .or. insertcol < 1) then
           WRITE(mat_lu,*) 'Cannot add block, the indexes', &
           & fullrow,fullcol,insertrow,insertcol, &
           & 'are out of the bounds - nrow, ncol =',A%nrow,A%ncol
           CALL lsQUIT('Cannot add block (subroutine mat_add_block)',mat_lu)
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_dense)
             call mat_dense_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            print *, "FALLBACK: mat_add_block converts to full for Block Sparse Matrices"
            call mat_add_block_bsm_fallback(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_unres_dense)
             call mat_unres_dense_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case default
              stop "mat_add_block not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('CREATE',mat_TSTR,mat_TEN,mat_lu)

       contains
         subroutine mat_add_block_bsm_fallback(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         implicit none
         integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
         real(Realk), intent(in) :: fullmat(fullrow,fullcol)
         type(Matrix), intent(inout) :: A
         real(realk), allocatable :: Afull(:,:)
         integer                  :: i, j

         allocate(Afull(A%nrow,A%ncol))
         call mat_to_full(A,1.0d0,Afull)

         do i = insertrow, insertrow+fullrow-1
            do j = insertcol, insertcol+fullcol-1
               Afull(i,j) = Afull(i,j)+fullmat(i-insertrow+1,j-insertcol+1)
            enddo
         enddo

         call mat_set_from_full(Afull,1.0d0,A) 
         deallocate(Afull)

       end subroutine mat_add_block_bsm_fallback

     END SUBROUTINE mat_add_block

!> \brief Retrieve block from type(matrix) 
!> \author T. Kjaergaard
!> \date 2009
!> \param A Input matrix from which we want to retrive a block
!> \param fullmat Return the desired block in this standard fortran matrix
!> \param fullrow Number of rows in fullmat
!> \param fullcol Number of columns in fullmat
!> \param insertrow Retrive block from A beginning at this row
!> \param insertcol Retrive block from A beginning at this col
      subroutine mat_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         implicit none
         integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
         real(Realk), intent(inout) :: fullmat(fullrow,fullcol)
         type(Matrix), intent(inout) :: A

         if (info_memory) write(mat_lu,*) 'Before mat_retrieve_block: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (insertrow+fullrow > A%nrow+1 .or. &
          &  insertcol+fullcol > A%ncol+1 .or. fullrow < 1 .or. fullcol < 1 .or. &
          & insertrow < 1 .or. insertcol < 1) then
           WRITE(mat_lu,*) 'Cannot retrieve block, the indexes', &
           & fullrow,fullcol,insertrow,insertcol, &
           & 'are out of the bounds - nrow, ncol =',A%nrow,A%ncol
           CALL lsQUIT('Cannot retrieve block (subroutine mat_retrieve_block)',mat_lu)
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_dense)
             call mat_dense_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_csr)
             call mat_csr_retrieve_block_full(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            print *, "FALLBACK: mat_retrieve_block converts to full for Block Sparse Matrices"
            call mat_retrieve_block_bsm_fallback(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_unres_dense)
             call mat_unres_dense_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case default
              stop "mat_retrieve_block not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('CREATE',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_retrieve_block: mem_allocated_global =', mem_allocated_global

       contains
         subroutine mat_retrieve_block_bsm_fallback(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         implicit none
         integer, intent(in)      :: fullrow,fullcol,insertrow,insertcol
         real(Realk), intent(out) :: fullmat(fullrow,fullcol)
         type(Matrix), intent(inout) :: A
         real(realk), allocatable :: Afull(:,:)
         integer                  :: i, j
         allocate(Afull(A%nrow,A%ncol))
         call mat_to_full(A,1.0d0,Afull)

         do i = insertrow, insertrow+fullrow-1
            do j = insertcol, insertcol+fullcol-1
                fullmat(i-insertrow+1,j-insertcol+1) = Afull(i,j)
            enddo
         enddo

         deallocate(Afull)

       end subroutine mat_retrieve_block_bsm_fallback
     END SUBROUTINE mat_retrieve_block

!> \brief Scale a type(matrix) A by a scalar alpha
!> \author L. Thogersen
!> \date 2003
!> \param alpha Scaling parameter
!> \param A Input/output matrix which we want to scale
      subroutine mat_scal(alpha,A)
         implicit none
         real(realk), intent(in) :: alpha
         type(Matrix), intent(inout) :: A

         if (info_memory) write(mat_lu,*) 'Before mat_scal: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_scal(alpha,A)
         case(mtype_dense)
             call mat_dense_scal(alpha,A)
         case(mtype_csr)
             call mat_csr_scal(alpha, A)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_scal(alpha,A)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            call bsm_scal(alpha, A)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_scal(alpha,A)
         case(mtype_unres_dense)
             call mat_unres_dense_scal(alpha,A)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_scal(alpha,A)
         case default
              stop "mat_scal not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('SCAL  ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_scal: mem_allocated_global =', mem_allocated_global
      end subroutine mat_scal

!> \brief Scale the diagonal of a type(matrix) A by a scalar alpha
!> \author L. Thogersen
!> \date 2003
!> \param alpha Scaling parameter
!> \param A Input/output matrix which we want to scale
      subroutine mat_scal_dia(alpha,A)
         implicit none
         real(realk), intent(in) :: alpha
         type(Matrix), intent(inout) :: A
         real(realk), allocatable :: afull(:,:)
         integer i

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (A%nrow /= A%ncol) then
           STOP 'cannot scale diagonal since ncol /= nrow'
         endif

         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_scal_dia(alpha,A)
         case(mtype_dense)
             call mat_dense_scal_dia(alpha,A)
#ifndef UNITTEST
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            CALL bsm_scal_dia(alpha, A)
#endif
         case(mtype_sparse1)
             call mat_sparse1_scal_dia(alpha,A)
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_scal_dia(alpha,A)
         case(mtype_unres_dense)
             call mat_unres_dense_scal_dia(alpha,A)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_scal_dia(alpha,A)
         case default
            print *, "FALLBACK scale_dia"
            allocate(afull(a%nrow, a%ncol))
            call mat_to_full(a,1D0,afull)
            do i = 1,A%nrow
               afull(i,i) = afull(i,i) * alpha
            enddo
            call mat_set_from_full(afull, 1D0, a)
            DEALLOCATE(afull)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('SCADIA',mat_TSTR,mat_TEN,mat_lu)
      end subroutine mat_scal_dia

!> \brief Set a type(matrix) A to zero
!> \author L. Thogersen
!> \date 2003
!> \param A Input/output matrix which should be set to zero
      subroutine mat_zero(A)
         implicit none
         type(Matrix), intent(inout) :: A

         if (info_memory) write(mat_lu,*) 'Before mat_zero: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_zero(A)
         case(mtype_dense)
             call mat_dense_zero(A)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_zero(A)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            call bsm_scal(0D0, A)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_zero(A)
         case(mtype_unres_dense)
             call mat_unres_dense_zero(A)
         case(mtype_csr)
             call mat_csr_zero(A)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_zero(A)
         case default
              stop "mat_zero not implemented for this type of matrix"
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('ZERO  ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_zero: mem_allocated_global =', mem_allocated_global
      end subroutine mat_zero

!> \brief Set upper or lower triangle of a type(matrix) to zero. Diagonal belongs to lower triangle.
!> \author L. Thogersen
!> \date 2003
!> \param part Indicates which part should be set to zero. 'UT'/'ut' = upper, 'LT'/'lt' = lower triangle
!> \param A Input/output matrix
      subroutine mat_zerohalf(part,A)
         implicit none
         character(len=2), intent(in) :: part
         type(Matrix), intent(inout) :: A
         real(realk), allocatable :: Afull(:,:)
         integer  i, j

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (part /= 'UT' .and. part /= 'ut' .and. part /= 'LT' .and. part /= 'lt') then
           STOP 'unknown part of matrix to zero'
         elseif (A%nrow /= A%ncol) then
           STOP 'cannot define triangles since nrow /= ncol'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_zerohalf(part,A)
         case(mtype_dense)
             call mat_dense_zerohalf(part,A)
#ifndef UNITTEST
#ifdef HAVE_BSM
          CASE(mtype_sparse_block)
             CALL bsm_zerohalf(part, A)
#endif
         case(mtype_sparse1)
             call mat_sparse1_zerohalf(part,A)
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_zerohalf(part,A)
         case(mtype_unres_dense)
             call mat_unres_dense_zerohalf(part,A)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_zerohalf(part,A)
         case default
            PRINT *, "FALLBACK zerohalf ", part
            allocate(afull(a%nrow, a%ncol))
            call mat_to_full(a,1D0,afull)
            if (part == 'ut' .or. part == 'UT') then
               !set the upper triangle to zero - diagonal is kept
               do j = 2,A%ncol
                  afull(1:j-1,j) = 0D0
               enddo
            else
               !set the lower triangle to zero - diagonal is also zeroed
               do j = 1,A%ncol
                  afull(j:A%nrow,j) = 0D0
               enddo
            endif
            call mat_set_from_full(afull, 1D0, a)
            DEALLOCATE(afull)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('ZEROHA',mat_TSTR,mat_TEN,mat_lu)
      end subroutine mat_zerohalf

!> \brief Write a type(matrix) to disk.
!> \author L. Thogersen
!> \date 2003
!> \param iunit Logical unit number of file which matrix should be written to
!> \param A Matrix which should be written on disk
      subroutine mat_write_to_disk(iunit,A)
         implicit none
         integer, intent(in) :: iunit
         type(Matrix), intent(in) :: A
         real(realk), allocatable :: afull(:,:)
#ifdef HAVE_BSM
         external mat_write_int, mat_write_real
#endif
         if (info_memory) write(mat_lu,*) 'Before mat_write_to_disk: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_write_to_disk(iunit,A)
         case(mtype_dense)
             call mat_dense_write_to_disk(iunit,A)
         case(mtype_csr)
             call mat_csr_write_to_disk(iunit,A)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_write_to_disk(iunit,A)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            call bsm_write_to_unit(iunit,A,mat_write_int,mat_write_real)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_write_to_disk(iunit,A)
         case(mtype_unres_dense)
             call mat_unres_dense_write_to_disk(iunit,A)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_write_to_disk(iunit,A)
         case default
            !print *, "FALLBACK: mat_write_to_disk"
            allocate(afull(a%nrow, a%ncol))
            call mat_to_full(a,1D0,afull)
            write(iunit) A%Nrow, A%Ncol
            write(iunit) afull
            deallocate(afull)

         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('WRITE ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_write_to_disk: mem_allocated_global =', mem_allocated_global
      end subroutine mat_write_to_disk

!> \brief Write a type(matrix) to disk in formatted form.
!> \author S. Host
!> \date 2007
!> \param iunit Logical unit number of file which matrix should be written to
!> \param A Matrix which should be written on disk
!>
!>  DEBUG ROUTINE (see debug_convert_density in debug.f90).
!>
      subroutine mat_write_to_disk2(iunit,A)
         implicit none
         integer, intent(in) :: iunit
         type(Matrix), intent(in) :: A
         real(realk), allocatable :: afull(:,:)
#ifdef HAVE_BSM
         external mat_write_int2, mat_write_real2
#endif
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_write_to_disk(iunit,A)
         case(mtype_dense)
             call mat_dense_write_to_disk2(iunit,A)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_write_to_disk(iunit,A)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            call bsm_write_to_unit(iunit,A,mat_write_int2,mat_write_real2)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_write_to_disk(iunit,A)
         case(mtype_unres_dense)
             call mat_unres_dense_write_to_disk(iunit,A)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_write_to_disk(iunit,A)
         case default
            !print *, "FALLBACK: mat_write_to_disk"
            allocate(afull(a%nrow, a%ncol))
            call mat_to_full(a,1D0,afull)
            write(iunit) A%Nrow, A%Ncol
            write(iunit) afull
            deallocate(afull)

         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('WRITE ',mat_TSTR,mat_TEN,mat_lu)
      end subroutine mat_write_to_disk2

!> \brief Add some logical auxiliary information to a type(matrix) on disk.
!> \author S. Host
!> \date June 2010
!> \param iunit Logical unit number of file containing the matrix
!> \param info Info to be written
!>
!> This is needed because the dens.restart that is dumped after calculatation
!> has ended can be either in the standard AO basis or the grand-canonical
!> basis. A dens.restart obtained from in the AO basis will not work in the
!> grand-canonical basis and vice versa. For the moment, it is only possible
!> to put true or false here, but that could easily be changed to a character
!> string saying e.g. 'AOBASIS', 'GCBASIS', or whatever other basis you could
!> think of. Should be independent of matrix type.
!> Must only be called after the matrix has been written and before 
!> file is rewinded!
!>
      subroutine mat_write_info_to_disk(iunit,info)
         implicit none
         integer, intent(in) :: iunit
         logical, intent(in) :: info

      write(iunit) info
      end subroutine mat_write_info_to_disk

!> \brief Read a type(matrix) from disk.
!> \author L. Thogersen
!> \date 2003
!> \param iunit Logical unit number of file from which matrix should be read
!> \param A Output matrix which is read from disk
      subroutine mat_read_from_disk(iunit,A)
         implicit none
         integer, intent(in) :: iunit
         type(Matrix), intent(inout) :: A  !output
         real(realk), allocatable :: afull(:,:)
         integer                  :: nrow, ncol
#ifdef HAVE_BSM
         external mat_read_int, mat_read_real
#endif
         if (info_memory) write(mat_lu,*) 'Before mat_read_from_disk: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_read_from_disk(iunit,A)
         case(mtype_dense)
             call mat_dense_read_from_disk(iunit,A)
         case(mtype_csr)
             call mat_csr_read_from_disk(iunit,A)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_read_from_disk(iunit,A)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            call bsm_read_from_unit(iunit,A,mat_read_int,mat_read_real)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_read_from_disk(iunit,A)
         case(mtype_unres_dense)
             call mat_unres_dense_read_from_disk(iunit,A)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_read_from_disk(iunit,A)
         case default
            !print *, "FALLBACK: mat_read_from_disk"
            allocate(afull(a%nrow, a%ncol))
            READ(iunit) Nrow, Ncol
            if(Nrow /= A%nrow) stop 'mat_read_from_disk: Nrow /= A%nrow'
            if(Ncol /= A%ncol) stop 'mat_read_from_disk: Ncol /= A%ncol'
            read(iunit) afull
            call mat_set_from_full(afull,1D0,a)
            deallocate(afull)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('READ  ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_read_from_disk: mem_allocated_global =', mem_allocated_global
      end subroutine mat_read_from_disk

!> \brief Read a type(matrix) from disk in formatted form.
!> \author S. Host
!> \date 2007
!> \param iunit Logical unit number of file from which matrix should be read
!> \param A Output matrix which should be read from disk
!>
!>  DEBUG ROUTINE (see debug_convert_density in debug.f90).
!>
      subroutine mat_read_from_disk2(iunit,A)
         implicit none
         integer, intent(in) :: iunit
         type(Matrix), intent(inout) :: A  !output
         real(realk), allocatable :: afull(:,:)
         integer                  :: nrow, ncol
#ifdef HAVE_BSM
         external mat_read_int2, mat_read_real2
#endif
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_read_from_disk(iunit,A)
         case(mtype_dense)
             call mat_dense_read_from_disk2(iunit,A)
#ifndef UNITTEST
         case(mtype_sparse1)
             call mat_sparse1_read_from_disk(iunit,A)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            call bsm_read_from_unit(iunit,A,mat_read_int2,mat_read_real2)
#endif
#endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_read_from_disk(iunit,A)
         case(mtype_unres_dense)
             call mat_unres_dense_read_from_disk(iunit,A)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_read_from_disk(iunit,A)
         case default
            !print *, "FALLBACK: mat_read_from_disk"
            allocate(afull(a%nrow, a%ncol))
            READ(iunit) Nrow, Ncol
            if(Nrow /= A%nrow) stop 'mat_read_from_disk: Nrow /= A%nrow'
            if(Ncol /= A%ncol) stop 'mat_read_from_disk: Ncol /= A%ncol'
            read(iunit) afull
            call mat_set_from_full(afull,1D0,a)
            deallocate(afull)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('READ  ',mat_TSTR,mat_TEN,mat_lu)
      end subroutine mat_read_from_disk2

!> \brief Read some logical auxiliary information from a type(matrix) on disk.
!> \author S. Host
!> \date June 2010
!> \param iunit Logical unit number of file containing the matrix
!> \param info Info to be read
!>
!> See description in mat_write_info_to_disk.
!>
      subroutine mat_read_info_from_disk(iunit,info)
         implicit none
         integer, intent(in)  :: iunit
         logical, intent(out) :: info

      read(iunit) info
      end subroutine mat_read_info_from_disk

!> \brief Change a vector to a matrix. The matrix must be symmetric or antisymmetric.
!> \author S. Host
!> \date 2005
!> \param symmetry Indicates symmetry of matrix, 's'/'S' = symmetric, 'a'/'A' = !antisymmetric
!> \param vec Input vector that should be transformed
!> \param mat Output matrix
!>
!> For debug purposes.
!> Symmetric: diagonal included, vecdim = matdim(matdim+1)/2
!> Antisym:   diagonal excluded, vecdim = matdim(matdim+1)/2 - matdim 
!>
    subroutine mat_vec_to_mat(symmetry, vec, mat)
    implicit none
                                                                                 
         integer :: n, m, i
         character, intent(in)      :: symmetry
         TYPE(matrix),intent(in)    :: vec
         TYPE(matrix),intent(inout) :: mat
     
    if  (symmetry /= 'a' .AND. symmetry /= 'A' .AND.       &
       & symmetry /= 's' .AND. symmetry /= 'S') then
          STOP 'Unknown symmetry possibility in mat_VEC_TO_MAT'
    endif

    if (symmetry == 's' .OR. symmetry == 'S') then     
       if (VEC%nrow /= MAT%nrow*(MAT%nrow+1)/2 .OR. &
         & VEC%nrow /= MAT%ncol*(MAT%ncol+1)/2 .OR. &
         & VEC%ncol /= 1) then
          STOP 'Wrong dimensions in mat_VEC_TO_MAT'
       endif
    endif

    if (symmetry == 'a' .OR. symmetry == 'A') then     
       if (VEC%nrow /= MAT%nrow*(MAT%nrow+1)/2 - MAT%nrow .OR. &
         & VEC%nrow /= MAT%ncol*(MAT%ncol+1)/2 - MAT%nrow .OR. &
         & VEC%ncol /= 1) then
          STOP 'Wrong dimensions in mat_VEC_TO_MAT'
       endif
    endif

    select case(matrix_type)
 !         case(mtype_symm_dense)
 !             call mat_symm_dense_vec_to_mat()
          case(mtype_dense)
              call mat_dense_vec_to_mat(symmetry, VEC, MAT)
#ifndef UNITTEST
          case(mtype_sparse1)
              call mat_sparse1_vec_to_mat(symmetry, VEC, MAT)
#endif
 !         case(mtype_unres_symm_dense)
 !             call mat_unres_symm_dense_vec_to_mat()
          case(mtype_unres_dense)
              call mat_unres_dense_vec_to_mat(symmetry, VEC, MAT)
 !         case(mtype_unres_sparse1)
 !             call mat_unres_sparse1_vec_to_mat()
          case default
               STOP "mat_VEC_TO_MAT not implemented for this type of matrix"
          end select
    end subroutine mat_VEC_TO_MAT
 
!> \brief Change a matrix to a vector. The matrix must be symmetric or antisymmetric.
!> \author S. Host
!> \date 2005
!> \param symmetry Indicates symmetry of matrix, 's'/'S' = symmetric, 'a'/'A' = !antisymmetric
!> \param mat Input matrix that should be transformed
!> \param vec Output vector
!>
!> For debug purposes.
!> Symmetric: diagonal included, vecdim = matdim(matdim+1)/2
!> Antisym:   diagonal excluded, vecdim = matdim(matdim+1)/2 - matdim 
!>
    subroutine mat_to_vec(symmetry, mat, vec)
    implicit none
                                    
         character, intent(in)      :: symmetry                                             
         TYPE(matrix),intent(inout) :: vec
         TYPE(matrix), intent(in)   :: mat

    if ( symmetry /= 'a' .AND. symmetry /= 'A' .AND.       &
       & symmetry /= 's' .AND. symmetry /= 'S') then
          STOP 'Unknown symmetry possibility in MAT_TO_VEC'
    endif

    if (symmetry == 's' .OR. symmetry == 'S') then     
       if (VEC%nrow /= MAT%nrow*(MAT%nrow+1)/2 .OR. &
         & VEC%nrow /= MAT%ncol*(MAT%ncol+1)/2 .OR. &
         & VEC%ncol /= 1) then
          STOP 'Wrong dimensions in MAT_TO_VEC'
       endif
    endif

    if (symmetry == 'a' .OR. symmetry == 'A') then     
       if (VEC%nrow /= MAT%nrow*(MAT%nrow+1)/2 - MAT%nrow .OR. &
         & VEC%nrow /= MAT%ncol*(MAT%ncol+1)/2 - MAT%nrow .OR. &
         & VEC%ncol /= 1) then
          STOP 'Wrong dimensions in MAT_TO_VEC'
       endif
    endif

    select case(matrix_type)
 !         case(mtype_symm_dense)
 !             call mat_symm_dense_mat_to_vec()
          case(mtype_dense)
              call mat_dense_mat_to_vec(symmetry, MAT, VEC)
#ifndef UNITTEST
          case(mtype_sparse1)
              call mat_sparse1_mat_to_vec(symmetry, MAT, VEC)
#endif
 !         case(mtype_unres_symm_dense)
 !             call mat_unres_symm_dense_mat_to_vec()
          case(mtype_unres_dense)
              call mat_unres_dense_mat_to_vec(symmetry, MAT, VEC)
 !         case(mtype_unres_sparse1)
 !             call mat_unres_sparse1_mat_to_vec()
          case default
               STOP "MAT_TO_VEC not implemented for this type of matrix"
          end select
    end subroutine MAT_TO_VEC

!> \brief Reports number of non-zero elements and sparsity of a type(matrix) A.
!> \author P. Salek
!> \date 2003
!> \param A Matrix for which number of non-zero elements i requested
!> \param mat_label The matrix is identified by this character string
!> \param iunit Logical unit number of file to which sparsity should be printed
!>
!> The data is printed to stream specified by iunit.
!> The matrix is identified by the mat_label string.
!>
    SUBROUTINE mat_report_sparsity(A,mat_label,nnz,iunit)
      implicit none
      TYPE(Matrix), INTENT(in) :: A
      CHARACTER(*), INTENT(IN) :: mat_label
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(out) :: nnz
      real(realk) :: sparsity
      
      if (info_memory) write(mat_lu,*) 'Before mat_report_sparsity: mem_allocated_global =', mem_allocated_global
      SELECT CASE(matrix_type)
      CASE(mtype_dense)
         CALL mat_dense_report_sparsity(A,sparsity)
      CASE(mtype_csr)
         CALL mat_csr_report_sparsity(A,sparsity)
#ifndef UNITTEST
      CASE(mtype_sparse1)
         call mat_sparse1_report_sparsity(A,sparsity)
#ifdef HAVE_BSM
      CASE(mtype_sparse_block)
         CALL bsm_get_sparsity(A, sparsity)
#endif
#endif
      CASE default
         return
      END SELECT

#ifndef UNITTEST
      WRITE(iunit,'("Matrix ",A," has nnz=",I10," sparsity: ",F10.3," %")')&
           &mat_label, INT(sparsity*A%nrow*A%ncol), sparsity*100.0
#endif
      nnz = INT(sparsity*A%nrow*A%ncol)
      if (info_memory) write(mat_lu,*) 'After mat_report_sparsity: mem_allocated_global =', mem_allocated_global
   END SUBROUTINE mat_report_sparsity

!Routines needed for purification
! - commented out because purification is not documented and no one really knows
! if purification works!

!!> \brief ?
!!> \author ?
!!> \date ?
!!> \param A ?
!!> \param B ?
!    subroutine mat_cholesky(A,B)
!      implicit none
!      type(Matrix), intent(in) :: A
!      type(Matrix), intent(inout) :: B
!
!      select case(matrix_type)
!!       case(mtype_symm_dense)
!!           call mat_symm_dense_cholesky(A,B)
!        case(mtype_dense)
!            call mat_dense_cholesky(A,B)
!#ifndef UNITTEST
!       case(mtype_sparse1)
!           call mat_sparse1_cholesky(A,B)
!#endif
!!       case(mtype_unres_symm_dense_cholseky)
!!           call mat_unres_symm_dense_cholesky(A,B)
!!       case(mtype_unres_dense_cholesky)
!!           call mat_unres_dense_cholesky(A,B)
!!       case(mtype_unres_sparse1_cholesky)
!!           call mat_unres_sparse1_cholesky(A,B)
!        case default
!            stop "mat_cholesky not implemented for this type of matrix"
!      end select
!
!      return
!    end subroutine mat_cholesky
!
!!> \brief ?
!!> \author ?
!!> \date ?
!!> \param A ?
!!> \param B ?
!    subroutine mat_inverse_triang(A,B)
!      implicit none
!      type(Matrix), intent(in) :: A
!      type(Matrix), intent(inout) :: B
!
!      select case(matrix_type)
!!       case(mtype_symm_dense)
!!           call mat_symm_dense_inverse_triang(A,B)
!        case(mtype_dense)
!            call mat_dense_inverse_triang(A,B)
!#ifndef UNITTEST
!        case(mtype_sparse1)
!            call mat_sparse1_inverse_triang(A,B)
!#endif
!!       case(mtype_unres_symm_dense_cholseky)
!!           call mat_unres_symm_dense_inverse_triang(A,B)
!!       case(mtype_unres_dense_inverse_triang)
!!           call mat_unres_dense_inverse_triang(A,B)
!!       case(mtype_unres_sparse1_inverse_triang)
!!           call mat_unres_sparse1_inverse_triang(A,B)
!        case default
!            stop "mat_inverse_triang not implemented for this type of matrix"
!      end select
!
!      return
!    end subroutine mat_inverse_triang
!
!!> \brief ?
!!> \author ?
!!> \date ?
!!> \param A ?
!!> \param B ?
!!> \param transb ?
!!> \param C ?
!    subroutine mat_simtran(A,B,transb,C)
!      implicit none
!      character, intent(in) :: transb
!      type(Matrix), intent(in) :: A,B
!      type(Matrix), intent(inout) :: C
!
!      select case(matrix_type)
!!       case(mtype_symm_dense)
!!           call mat_symm_dense_simtran(A,B,transb,C)
!        case(mtype_dense)
!            call mat_dense_simtran(A,B,transb,C)
!#ifndef UNITTEST
!        case(mtype_sparse1)
!            call mat_sparse1_simtran(A,B,transb,C)
!#endif
!!       case(mtype_unres_symm_dense_cholseky)
!!           call mat_unres_symm_dense_simtran(A,B,transb,C)
!!       case(mtype_unres_dense_simtran)
!!           call mat_unres_dense_simtran(A,B,transb,C)
!!       case(mtype_unres_sparse1_simtran)
!!           call mat_unres_sparse1_simtran(A,B,transb,C)
!        case default
!            stop "mat_simtran not implemented for this type of matrix"
!      end select
!
!      return
!    end subroutine mat_simtran
!
!!> \brief ?
!!> \author ?
!!> \date ?
!!> \param A ?
!!> \param lowcut ?
!    subroutine mat_clean(A,lowcut)
!      implicit none
!      type(Matrix), intent(inout) :: A
!      real(realk), intent(in) :: lowcut
!
!      select case(matrix_type)
!!       case(mtype_symm_dense)
!!           call mat_symm_dense_clean(A,lowcut)
!        case(mtype_dense)
!            call mat_dense_clean(A,lowcut)
!#ifndef UNITTEST
!        case(mtype_sparse1)
!            call mat_sparse1_clean(A,lowcut)
!#endif
!!       case(mtype_unres_symm_dense_cholseky)
!!           call mat_unres_symm_dense_clean(A,lowcut)
!!       case(mtype_unres_dense_clean)
!!           call mat_unres_dense_clean(A,lowcut)
!!       case(mtype_unres_sparse1_clean)
!!           call mat_unres_sparse1_clean(A,lowcut)
!        case default
!            stop "mat_clean not implemented for this type of matrix"
!      end select
!
!      return
!    end subroutine mat_clean
!
!!> \brief ?
!!> \author ?
!!> \date ?
!!> \param A ?
!!> \param B ?
!    subroutine mat_to_minus_one_half(A,B)
!      implicit none
!      type(Matrix), intent(in) :: A
!      type(Matrix), intent(inout) :: B
!
!      select case(matrix_type)
!!       case(mtype_symm_dense)
!!           call mat_symm_dense_to_minus_one_half(A,B)
!        case(mtype_dense)
!            call mat_dense_to_minus_one_half(A,B)
!!       case(mtype_sparse1)
!!           call mat_sparse1_to_minus_one_half(A,B)
!!       case(mtype_unres_symm_dense_cholseky)
!!           call mat_unres_symm_dense_to_minus_one_half(A,B)
!!       case(mtype_unres_dense_to_minus_one_half)
!!           call mat_unres_dense_to_minus_one_half(A,B)
!!       case(mtype_unres_sparse1_to_minus_one_half)
!!           call mat_unres_sparse1_to_minus_one_half(A,B)
!        case default
!            stop "mat_to_minus_one_half not implemented for this type of matrix!"
!      end select
!
!      return
!    end subroutine mat_to_minus_one_half
!
!
!!> \brief ?
!!> \author ?
!!> \date ?
!!> \param A ?
!!> \param min ?
!!> \param max ?
!    subroutine mat_gershgorin_minmax(A,min,max)
!      implicit none
!      type(Matrix), intent(in) :: A
!      real(realk), intent(out) :: min,max
!
!      select case(matrix_type)
!!       case(mtype_symm_dense)
!!           call mat_symm_dense_gershgorin_minmax(A,min,max)
!        case(mtype_dense)
!            call mat_dense_gershgorin_minmax(A,min,max)
!#ifndef UNITTEST
!        case(mtype_sparse1)
!            call mat_sparse1_gershgorin_minmax(A,min,max)
!#endif
!!       case(mtype_unres_symm_dense_cholseky)
!!           call mat_unres_symm_dense_gershgorin_minmax(A,min,max)
!!       case(mtype_unres_dense_gershgorin_minmax)
!!           call mat_unres_dense_gershgorin_minmax(A,min,max)
!!       case(mtype_unres_sparse1_gershgorin_minmax)
!!           call mat_unres_sparse1_gershgorin_minmax(A,min,max)
!        case default
!            stop "mat_gershgorin_minmax not implemented for this type of matrix!"
!      end select
!
!      return
!    end subroutine mat_gershgorin_minmax
!
!> \brief Returns sum of all elements of matrix.
!> \param A The input matrix 
!> \return The sum of all elements of matrix
    function mat_sum(A)
      implicit none
      type(Matrix), intent(in) :: A
      real(realk) :: mat_sum

      select case(matrix_type)
!       case(mtype_symm_dense)
!           mat_sum=mat_symm_dense_sum(A)
        case(mtype_dense)
            mat_sum=mat_dense_sum(A)
#ifndef UNITTEST
        case(mtype_sparse1)
            mat_sum=mat_sparse1_sum(A)
#endif
!       case(mtype_unres_symm_dense_cholseky)
!           mat_sum=mat_unres_symm_dense_sum(A)
       case(mtype_unres_dense)
           mat_sum=mat_unres_dense_sum(A)
!       case(mtype_unres_sparse1_sum)
!           mat_sum=mat_unres_sparse1_sum(A)
        case default
            stop "mat_sum not implemented for this type of matrix"
      end select

      return
    end function mat_sum

#ifndef UNITTEST
!> \brief Inquire zero cutoff - for sparse1 matrices only!! 
!> \author S. Host
!> \date 2009
!> \param cutoff The zero cutoff for sparse1 matrices
    subroutine mat_inquire_cutoff(cutoff)
      implicit none
      real(realk), intent(out) :: cutoff

      select case(matrix_type)
!       case(mtype_symm_dense)
!           call mat_symm_dense_zero_cutoff(cutoff)
!       case(mtype_dense)
!           call mat_dense_zero_cutoff(cutoff)
       case(mtype_sparse1)
           call mat_sparse1_inquire_cutoff(cutoff)
!       case(mtype_unres_symm_dense_cholseky)
!           call mat_unres_symm_dense_zero_cutoff(cutoff)
!       case(mtype_unres_dense_zero_cutoff)
!           call mat_unres_dense_zero_cutoff(cutoff)
!       case(mtype_unres_sparse1_zero_cutoff)
!           call mat_unres_sparse1_zero_cutoff(cutoff)
        case default
            stop "mat_zero_cutoff not implemented for this type of matrix"
      end select

      return
    end subroutine mat_inquire_cutoff
#endif

!> \brief Extract diagonal of A, store in dense vector vec.
!> \author B. Jansik
!> \date 2010
!> \param A The type(matrix) input
!> \param diag vector to hold the diagonal
    subroutine mat_extract_diagonal (diag,A)
      implicit none
      real(realk), intent(out) :: diag(:)
      type(Matrix), intent(in) :: A

      select case(matrix_type)
      case(mtype_dense)
           call mat_dense_extract_diagonal(diag,A)
      case default
            stop "mat_extract_diagonal not implemented for this type of matrix"
      end select

    end subroutine mat_extract_diagonal

!> \brief Make c = alpha*diag(a)b + beta*c, where a is realk(:) b,c are type(matrix) and alpha,beta are parameters
!> \author B. Jansik
!> \date 2010
!> \param a The first realk(:) diagonal 
!> \param b The second type(matrix) factor
!> \param transb 'T'/'t' if b should be transposed, 'N'/'n' otherwise
!> \param alpha The alpha parameter
!> \param beta The beta parameter
!> \param c The output type(matrix)
    subroutine mat_dmul (a, b, transb, alpha, beta, c)
         implicit none
         real(realk), intent(in)  :: a(:)
         TYPE(Matrix), intent(IN) :: b
         character, intent(in)    :: transb
         REAL(realk), INTENT(IN)  :: alpha, beta
         TYPE(Matrix), intent(inout):: c
         integer :: ak, bk, ci, cj
         select case(matrix_type)
         case(mtype_dense)
            call mat_dense_dmul(a,b,transb,alpha,beta,c)
         case default
              stop "mat_dmul not implemented for this type of matrix"
         end select

    end subroutine mat_dmul


#ifndef UNITTEST
!> \brief Set zero cutoff - for sparse1 matrices only!! 
!> \author S. Host
!> \date 2009
!> \param zero The desired zero cutoff for sparse1 matrices
    subroutine mat_zero_cutoff(zero)
      implicit none
      real(realk), intent(in) :: zero

      select case(matrix_type)
!       case(mtype_symm_dense)
!           call mat_symm_dense_zero_cutoff(zero)
!       case(mtype_dense)
!           call mat_dense_zero_cutoff(zero)
       case(mtype_sparse1)
           call mat_sparse1_zero_cutoff(zero)
!       case(mtype_unres_symm_dense_cholseky)
!           call mat_unres_symm_dense_zero_cutoff(zero)
!       case(mtype_unres_dense_zero_cutoff)
!           call mat_unres_dense_zero_cutoff(zero)
!       case(mtype_unres_sparse1_zero_cutoff)
!           call mat_unres_sparse1_zero_cutoff(zero)
        case default
            stop "mat_zero_cutoff not implemented for this type of matrix"
      end select

      return
    end subroutine mat_zero_cutoff
#endif

#ifdef VAR_MPI
!> \brief ?
!> \author ?
!> \date ?
!> \param A ?
!> \param root ?
    subroutine mat_mpixbcast(A,root)
         use infpar_module
         implicit none
         integer, intent(in) :: root
         integer :: bcnt
         type(Matrix) :: A, DAT

        !external mat_mpixbcast_int, mat_mpixbcast_real
         external mat_count_int,mat_count_real,mat_mwrite_int,mat_mwrite_real,&
                 &mat_mread_int,mat_mread_real

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)

         select case(matrix_type)

         case(mtype_dense)
           call mpixbcast(A%elms,A%nrow*A%ncol,'DOUBLE',root)

#ifdef HAVE_BSM
         case(mtype_sparse_block)
           if (root.eq.infpar%mynum) then
            !call bsm_write_to_unit (root,A,mat_mpixbcast_int,mat_mpixbcast_real)
             call bsm_write_to_unit(bcnt,A,mat_count_int,mat_count_real)
             call mat_dense_init(DAT,(bcnt/realk)+1,1)
             call mat_msetstart(DAT,DAT%elms)
             call bsm_write_to_unit(DAT,A,mat_mwrite_int,mat_mwrite_real)
             call mpixbcast(bcnt,1,'INTEGER',root)
             call mpixbcast(DAT%elms,bcnt,'CHARACTER',root)
           else  
            !call bsm_read_from_unit(root,A,mat_mpixbcast_int,mat_mpixbcast_real)
             call mpixbcast(bcnt,1,'INTEGER',root)
             call mat_dense_init(DAT,(bcnt/8)+1,1)
             call mpixbcast(DAT%elms,bcnt,'CHARACTER',root)
             call mat_msetstart(DAT,DAT%elms)
             call bsm_read_from_unit(DAT,A,mat_mread_int,mat_mread_real)
           end if

           call mat_dense_free(DAT)
#endif
         case default
           call lsquit('mat_mpixbcast() not implemented for this type of matrix',mat_lu)

         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MBCAST',mat_TSTR,mat_TEN,mat_lu)
    end subroutine mat_mpixbcast

!> \brief ?
!> \author ?
!> \date ?
!> \param A ?
!> \param source ?
!> \param tag ?
    subroutine mat_mpixrecv(A,source,tag)
         use infpar_module
         implicit none
         integer, intent(in)  :: source, tag
         integer ::  bcnt !,dest(2)
         type(Matrix) :: A, DAT

        !external mat_mpixrecv_int, mat_mpixrecv_real
         external mat_mread_int, mat_mread_real

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)

         select case(matrix_type)
         case(mtype_dense)
           call mpixrecv(A%elms,A%nrow*A%ncol,'DOUBLE',source,tag)
#ifdef HAVE_BSM
         case(mtype_sparse_block)
            !dest(1)=source; dest(2)=1000*tag
            !call bsm_read_from_unit(dest,A,mat_mpixrecv_int,mat_mpixrecv_real)
             call mpixrecv(bcnt,1,'INTEGER',source,tag)
             call mat_dense_init(DAT,(bcnt/realk)+1,1)
             call mat_msetstart(DAT,DAT%elms)
             call mpixrecv(DAT%elms,bcnt,'CHARACTER',source,tag)
             call bsm_read_from_unit(DAT,A,mat_mread_int,mat_mread_real)
             call mat_dense_free(DAT)
#endif
         case default
           call lsquit('mat_mpixrecv() not implemented for this type of matrix',mat_lu)

         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MRECV ',mat_TSTR,mat_TEN,mat_lu)
    end subroutine mat_mpixrecv

!> \brief ?
!> \author ?
!> \date ?
!> \param A ?
!> \param sendto ?
!> \param tag ?
    subroutine mat_mpixsend(A,sendto,tag)
         use infpar_module
         implicit none
         integer, intent(in) :: sendto, tag
         integer ::  bcnt !,dest(2)
         type(Matrix) :: A, DAT

         !external mat_mpixsend_int, mat_mpixsend_real
          external mat_count_int, mat_count_real, mat_mwrite_int, mat_mwrite_real

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)

         select case(matrix_type)
         case(mtype_dense)
           call mpixsend(A%elms,A%nrow*A%ncol,'DOUBLE',sendto,tag)

#ifdef HAVE_BSM
         case(mtype_sparse_block)
            !dest(1)=sendto; dest(2)=1000*tag
            !call bsm_write_to_unit(dest,A,mat_mpixsend_int,mat_mpixsend_real)
             call bsm_write_to_unit(bcnt,A,mat_count_int,mat_count_real)
             call mat_dense_init(DAT,(bcnt/realk)+1,1)
             call mat_msetstart(DAT,DAT%elms)
             call bsm_write_to_unit(DAT,A,mat_mwrite_int,mat_mwrite_real)
             call mpixsend(bcnt,1,'INTEGER',sendto,tag)
             call mpixsend(DAT%elms,bcnt,'CHARACTER',sendto,tag)
             call mat_dense_free(DAT)
#endif
         case default
           call lsquit('mat_mpixsend() not implemented for this type of matrix',mat_lu)

         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MSEND ',mat_TSTR,mat_TEN,mat_lu)
    end subroutine mat_mpixsend


#endif

END MODULE Matrix_Operations

!> \brief Standalone routine for BSM IO support
!> \author P. Salek
!> \date 2003
!> \param iunit ?
!> \param cnt ?
!> \param idata ?
subroutine mat_write_int(iunit,cnt,idata)
   use Matrix_module
   integer, intent(in) :: iunit, cnt
   integer, intent(in) :: idata(cnt)
   write(iunit) idata
end subroutine mat_write_int

!> \brief Standalone routine for BSM IO support
!> \author P. Salek
!> \date 2003
!> \param iunit ?
!> \param cnt ?
!> \param idata ?
subroutine mat_read_int(iunit,cnt,idata)
   use Matrix_module
   integer, intent(in)  :: iunit, cnt
   integer, intent(out) :: idata(cnt)
   read(iunit) idata
end subroutine mat_read_int

!> \brief Standalone routine for BSM IO support
!> \author P. Salek
!> \date 2003
!> \param iunit ?
!> \param cnt ?
!> \param data ?
subroutine mat_write_real(iunit,cnt,data)
   use Matrix_module
   integer,     intent(in) :: iunit, cnt
   real(realk), intent(in) :: data(cnt)
   write(iunit) data
end subroutine mat_write_real

!> \brief Standalone routine for BSM IO support
!> \author P. Salek
!> \date 2003
!> \param iunit ?
!> \param cnt ?
!> \param data ?
subroutine mat_read_real(iunit,cnt,data)
   use Matrix_module
   integer,     intent(in)  :: iunit, cnt
   real(realk), intent(out) :: data(cnt)
   read(iunit) data
end subroutine mat_read_real

!Hack routines - see debug_convert_density in debug.f90
subroutine mat_write_int2(iunit,cnt,idata)
   use Matrix_module
   integer, intent(in) :: iunit, cnt
   integer, intent(in) :: idata(cnt)
   write(iunit,*) idata
end subroutine mat_write_int2

subroutine mat_read_int2(iunit,cnt,idata)
   use Matrix_module
   integer, intent(in)  :: iunit, cnt
   integer, intent(out) :: idata(cnt)
   read(iunit,*) idata
end subroutine mat_read_int2

subroutine mat_write_real2(iunit,cnt,data)
   use Matrix_module
   integer,     intent(in) :: iunit, cnt
   real(realk), intent(in) :: data(cnt)
   write(iunit,*) data
end subroutine mat_write_real2

subroutine mat_read_real2(iunit,cnt,data)
   use Matrix_module
   integer,     intent(in)  :: iunit, cnt
   real(realk), intent(out) :: data(cnt)
   read(iunit,*) data
end subroutine mat_read_real2
!End hack routines

#ifdef VAR_MPI
#ifdef HAVE_BSM
#if 0

!> \brief Used for sparse matrix transfer without buffer. Elegant, but not very efficient
!> \author P. Salek
!> \date 2003
!> \param root ?
!> \param cnt ?
!> \param idata ?
    subroutine mat_mpixbcast_int(root,cnt,idata)
      use Matrix_module
      integer, intent(in) :: root, cnt
      integer, intent(in) :: idata(cnt)
!     call mpixbcast(cnt,  1,  'INTEGER',root)
      call mpixbcast(idata,cnt,'INTEGER',root)
    end subroutine mat_mpixbcast_int

!> \brief Used for sparse matrix transfer without buffer. Elegant, but not very efficient
!> \author P. Salek
!> \date 2003
!> \param root ?
!> \param cnt ?
!> \param data ?
    subroutine mat_mpixbcast_real(root,cnt,data)
       use Matrix_module
       integer, intent(in) :: root, cnt
       real(realk), intent(in) :: data(cnt)
!      call mpixbcast(cnt, 1,  'INTEGER',root)
       call mpixbcast(data,cnt,'DOUBLE', root)
    end subroutine mat_mpixbcast_real

!> \brief Used for sparse matrix transfer without buffer. Elegant, but not very efficient
!> \author P. Salek
!> \date 2003
!> \param dest ?
!> \param cnt ?
!> \param idata ?
    subroutine mat_mpixrecv_int(dest,cnt,idata)
      use Matrix_module
      integer  :: dest(2), cnt
      integer, intent(in) :: idata(cnt)
!     call mpixrecv(cnt,  1,  'INTEGER',dest(1),dest(2))
!           dest(2) = dest(2) + 1
      call mpixrecv(idata,cnt,'INTEGER',dest(1),dest(2))
            dest(2) = dest(2) + 1
    end subroutine mat_mpixrecv_int

!> \brief Used for sparse matrix transfer without buffer. Elegant, but not very efficient
!> \author P. Salek
!> \date 2003
!> \param dest ?
!> \param cnt ?
!> \param data ?
    subroutine mat_mpixrecv_real(dest,cnt,data)
       use Matrix_module
       integer :: dest(2), cnt
       real(realk), intent(in) :: data(cnt)
!      call mpixrecv(cnt, 1,  'INTEGER',dest(1),dest(2))
!           dest(2) = dest(2) + 1
       call mpixrecv(data,cnt,'DOUBLE', dest(1),dest(2))
            dest(2) = dest(2) + 1
    end subroutine mat_mpixrecv_real

!> \brief Used for sparse matrix transfer without buffer. Elegant, but not very efficient
!> \author P. Salek
!> \date 2003
!> \param dest ?
!> \param cnt ?
!> \param idata ?
    subroutine mat_mpixsend_int(dest,cnt,idata)
      use Matrix_module
      integer :: dest(2), cnt
      integer, intent(in) :: idata(cnt)
!     call mpixsend(cnt,  1,  'INTEGER',dest(1),dest(2))
!           dest(2) = dest(2) + 1
      call mpixsend(idata,cnt,'INTEGER',dest(1),dest(2))
            dest(2) = dest(2) + 1
    end subroutine mat_mpixsend_int

!> \brief Used for sparse matrix transfer without buffer. Elegant, but not very efficient
!> \author P. Salek
!> \date 2003
!> \param dest ?
!> \param cnt ?
!> \param data ?
    subroutine mat_mpixsend_real(dest,cnt,data)
       use Matrix_module
       integer  :: dest(2), cnt
       real(realk), intent(in) :: data(cnt)
!      call mpixsend(cnt, 1,  'INTEGER',dest(1),dest(2))
!           dest(2) = dest(2) + 1
       call mpixsend(data,cnt,'DOUBLE', dest(1),dest(2))
            dest(2) = dest(2) + 1
    end subroutine mat_mpixsend_real


#endif
#endif /* HAVE_BSM */
#endif /* VAR_MPI  */
