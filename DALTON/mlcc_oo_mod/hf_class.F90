module hf_class
!
   use workspace
   use cholesky_integrals_class
   use input_output
!
   implicit none
!
   type hartree_fock
      integer(i15)                          :: n_occ
      integer(i15)                          :: n_vir
      integer(i15)                          :: n_mo 
      real(dp), dimension(:,:), allocatable :: mo_coef
      type(cholesky_integrals)              :: cholesky
      character(len=10)                     :: hf_interface_program
   contains
      procedure                             :: init => init_hartree_fock
   end type hartree_fock
!
!  Private submodules and functions
!
   private :: read_sirius
!
contains 
!
   subroutine init_hartree_fock(hf)
!
      implicit none
!
      class(hartree_fock) :: hf
!
      if (hf % hf_interface_program .eq. 'DALTON    ') then
!
         call read_sirius(hf)        
!
      endif
!
!     Allocate the MO coefficients (move this to read_sirius)
!
      call allocator(hf % mo_coef, hf % n_mo, hf % n_mo) 
      hf % mo_coef = zero
!
      call hf % cholesky % init (hf % mo_coef, hf % n_occ, hf % n_vir)
!
   end subroutine init_hartree_fock
!
   subroutine read_sirius(hf)
!
!
!
      use workspace
!
      implicit none
!
      class(hartree_fock) :: hf
      integer(i15)        :: unit_identifier_sirius = -1 
      integer(i15)        :: idummy = 0
      integer(i15)        :: n_symmetries, n_basis_sym, n_orbitals_sym
      integer(i15)        :: i,j
!      
!     
!     Open Sirius Fock file
!     ---------------------
!

      open(unit=unit_identifier_sirius,file='mlcc_sirius',status='old',form='unformatted')
      rewind(unit_identifier_sirius)

     
      call generate_unit_identifier(unit_identifier_sirius)
!
   end subroutine read_sirius
end module hf_class
