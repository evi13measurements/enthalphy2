module cholesky_integrals_class
!
   use mlcc_types
!
   type cholesky_integrals
      integer(i15) :: n_J  = 0
      integer(i15) :: n_ao = 0
      integer(i15) :: n_mo = 0
   contains
!	
      procedure :: init => init_cholesky_integrals
!
   end type cholesky_integrals
!
contains
!
   subroutine init_cholesky_integrals(chol,mo_coef,n_occ,n_vir)
!
      implicit none
!
      class(cholesky_integrals) :: chol
!
      real(dp), dimension(chol%n_mo,chol%n_mo) :: mo_coef
!
      integer(i15) :: n_occ, n_vir 
!
!     Read the Cholesky vectors from file (in the AO basis),
!     transform them to the MO basis, and save to file
!
      call get_cholesky_vectors(chol,mo_coef,n_occ,n_vir)
!
   end subroutine init_cholesky_integrals
!
end module cholesky_integrals_class
