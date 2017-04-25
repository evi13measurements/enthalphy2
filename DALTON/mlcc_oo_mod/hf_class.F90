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
      real(dp), dimension(:,:), allocatable :: fock_diagonal
      real(dp)                              :: nuclear_potential, scf_energy
      type(cholesky_integrals)              :: cholesky
   contains
      procedure                             :: init => init_hartree_fock
   end type hartree_fock
!
!  Private submodules and functions
!
   private :: read_hf_info
!
contains 
!
   subroutine init_hartree_fock(hf)
!
      implicit none
!
      class(hartree_fock) :: hf
!
!
      call read_hf_info(hf)        
!
      call hf % cholesky % init (hf % mo_coef, hf % n_occ, hf % n_vir)
!
   end subroutine init_hartree_fock
!
   subroutine read_hf_info(hf)
!
!
!
      use workspace
!
      implicit none
!
      class(hartree_fock) :: hf
      integer(i15)        :: unit_identifier_hf = -1 
      integer(i15)        :: n_lambda
      integer(i15)        :: i,j
!      
!     
!     Open mlcc_hf_info
!     ---------------------
!
      call generate_unit_identifier(unit_identifier_hf)
      open(unit=unit_identifier_hf,file='mlcc_hf_info',status='old',form='unformatted')
      rewind(unit_identifier_hf)
!
!     Read mlcc_hf_info
!     ---------------------
!
      read(unit_identifier_hf,*)  hf % n_mo, hf % n_occ, n_lambda, hf % nuclear_potential, hf % scf_energy
!
!     Setting n_vir
!
      hf % n_vir          = hf % n_mo - hf % n_occ
!      
!     Allocate space for Fock diagonal and coefficients.
!
      call allocator(hf % fock_diagonal,hf % n_mo,1)
      hf % fock_diagonal = zero
!      
      call allocator(hf % mo_coef,n_lambda,1)
      hf % mo_coef = zero
!
!     Read in Fock diagonal and coefficients
!
      read(unit_identifier_hf,*) (hf % fock_diagonal(i,1),i=1,hf % n_mo)
      read(unit_identifier_hf,*) (hf % mo_coef(i,1),i=1,n_lambda) 
!
!     Done with file
!     ------------------
!    
      close(unit_identifier_hf)
!
   end subroutine read_hf_info
end module hf_class
