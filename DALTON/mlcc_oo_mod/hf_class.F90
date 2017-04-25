module hf_class
!
   use workspace
   use cholesky_integrals_class
   use input_output
!
   implicit none
!
   type hartree_fock
!
      integer(i15)                          :: n_occ
      integer(i15)                          :: n_vir
      integer(i15)                          :: n_mo 
!
      real(dp), dimension(:,:), allocatable :: mo_coef
      real(dp), dimension(:,:), allocatable :: fock_diagonal
      real(dp), dimension(:,:), allocatable :: fock_matrix
!
      real(dp)                              :: nuclear_potential, scf_energy
      type(cholesky_integrals)              :: cholesky
!
   contains
!
      procedure                             :: init => init_hartree_fock
!
   end type hartree_fock
!
!  Private submodules and functions
!
   private :: read_hf_info
!
contains 
!
   subroutine init_hartree_fock(wavefn)
!
!!  Initialization of hartree-fock object
!
!   Calls read_hf_info - reads MLCC_HF_INFO
!   Initializes Cholesky vectors  
!   Allocates Fock matrix and sets it to 0. Fock matrix is constructed in derived types.
!  
      implicit none
!
      class(hartree_fock) :: wavefn
!
!     Initializing HF variables
!      
      write(unit_output,*) 'Initializing HF...'
      call read_hf_info(wavefn)          
!
!     Initializing integral and Cholesky variables
!  
      write(unit_output,*) 'Initializing Cholesky...'
      call wavefn % cholesky % init (wavefn % mo_coef, wavefn % n_occ, wavefn % n_vir)
!
!     Allocate Fock matrix and set to 0
!
      call allocator(wavefn % fock_matrix, wavefn % n_mo, wavefn % n_mo)
      wavefn % fock_matrix = zero
!
   end subroutine init_hartree_fock

   subroutine read_hf_info(wavefn)
!
!  Reads MLCC_HF_INFO and initializes HF variables n_occ, n_vir, n_mo, orbital_coef, fock_diagonal
!
!  MLCC_HF_INFO is written in the mlcc_write_sirifc subroutine 
!  and called from wr_sirifc subroutine in siropt module.
!  
      use workspace
!
      implicit none
!
      class(hartree_fock) :: wavefn
      integer(i15)        :: unit_identifier_hf = -1 
      integer(i15)        :: n_lambda
      integer(i15)        :: i,j
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
      read(unit_identifier_hf,*) wavefn % n_mo, wavefn % n_occ, n_lambda, wavefn % nuclear_potential, wavefn % scf_energy
!
!     Setting n_vir
!
      wavefn % n_vir = wavefn % n_mo - wavefn % n_occ
!      
!     Allocate space for Fock diagonal and coefficients.
!
      call allocator(wavefn % fock_diagonal, wavefn % n_mo,1)
      wavefn % fock_diagonal = zero
!      
      call allocator(wavefn % mo_coef,n_lambda,1)
      wavefn % mo_coef = zero
!
!     Read in Fock diagonal and coefficients
!
      read(unit_identifier_hf,*) (wavefn % fock_diagonal(i,1),i=1,wavefn % n_mo)
      read(unit_identifier_hf,*) (wavefn % mo_coef(i,1),i=1,n_lambda) 
!
!     Done with file
!     ------------------
!    
      close(unit_identifier_hf)
!
   end subroutine read_hf_info
!
end module hf_class
