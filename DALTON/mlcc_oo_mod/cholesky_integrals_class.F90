module cholesky_integrals_class
!
!  Cholesky integrals class
!  Written by Sarai D. Folkestad and Eirik F. Kjønstad, 21 Apr 2017
!
   use input_output
   use mlcc_types
!
   type cholesky_integrals
!
      integer(i15) :: n_J  = 0 ! Number of Cholesky vectors 
      integer(i15) :: n_ao = 0 ! Number of atomic orbitals
!
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
      integer(i15) :: n_occ,n_vir 
      real(dp), dimension(n_occ+n_vir,n_occ+n_vir) :: mo_coef
!
      write(unit_output,*) 'In init_cholesky_integrals'
      call flshfo(unit_output)
!
!     Read the Cholesky vectors from file (in the AO basis),
!     transform them to the MO basis, and save to file
!
      call read_and_transform_cholesky_vectors(chol,mo_coef,n_occ,n_vir)
!
   end subroutine init_cholesky_integrals
!
end module cholesky_integrals_class
