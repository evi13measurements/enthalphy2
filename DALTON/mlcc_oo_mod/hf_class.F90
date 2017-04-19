module hf_class
!
   use cholesky_integrals_class
   use printing
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
      hf % n_occ = 10
      hf % n_vir = 100
!
      write(unit_output,*) 'Number of occ.:', hf % n_occ
      write(unit_output,*) 'Number of vir.:', hf % n_vir 
!
   end subroutine init_hartree_fock
!
end module hf_class
