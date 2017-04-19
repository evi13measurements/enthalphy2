module hf_class
!
   use cholesky_integrals_class
   use printing
!
   implicit none
!
   type hartree_fock
      integer                  :: n_occ
      integer                  :: n_vir
      type(cholesky_integrals) :: cholesky
   contains
      procedure                :: init => init_hartree_fock
   end type hartree_fock
!
contains 
!
   subroutine init_hartree_fock(hf)
!
      implicit none
!
      class(hartree_fock) :: hf
!
   end subroutine init_hartree_fock
!
end module hf_class
