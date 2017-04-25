module hf_class
!
!  Date:    April 2017
!  Authors: Sarai D. Folkestad and Eirik F. Kjønstad
   use workspace
   use cholesky_integrals_class
   use mlcc_types
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
      real(dp), dimension(:,:), allocatable :: fock_matrix_ij
      real(dp), dimension(:,:), allocatable :: fock_matrix_ia
      real(dp), dimension(:,:), allocatable :: fock_matrix_ai
      real(dp), dimension(:,:), allocatable :: fock_matrix_ab
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
!  Initialization of Hartree-Fock object
!
!  Calls read_hf_info - reads MLCC_HF_INFO
!  Initializes Cholesky vectors  
!  Allocates Fock matrix and sets it to 0. Fock matrix is constructed in derived types.
!  
      implicit none
!
      class(hartree_fock) :: wavefn
!
      write(unit_output,*) 'In init_hartree_fock'
      call flshfo(unit_output)
!
!     Initializing HF variables
!
      call read_hf_info(wavefn)        
!
!     Initializing integral and Cholesky variables
!  
      call wavefn % cholesky % init (wavefn % mo_coef, wavefn % n_occ, wavefn % n_vir)
!
!     Allocate Fock matrix and set to 0
!
      write(unit_output,*) 'Allocate Fock matrix blocks...'
      call allocator(wavefn % fock_matrix_ij, wavefn % n_occ, wavefn % n_occ)
      call allocator(wavefn % fock_matrix_ia, wavefn % n_occ, wavefn % n_vir)
      call allocator(wavefn % fock_matrix_ai, wavefn % n_vir, wavefn % n_occ)
      call allocator(wavefn % fock_matrix_ab, wavefn % n_vir, wavefn % n_vir)

!
      wavefn % fock_matrix_ij = zero
      wavefn % fock_matrix_ia = zero
      wavefn % fock_matrix_ai = zero
      wavefn % fock_matrix_ab = zero
!
   end subroutine init_hartree_fock

   subroutine read_hf_info(wavefn)
!
!  Date:    April 2017
!  Authors: Sarai D. Folkestad and Eirik F. Kjønstad
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
      open(unit=unit_identifier_hf,file='mlcc_hf_info',status='old',form='formatted')
      rewind(unit_identifier_hf)
!
!     Read mlcc_hf_info
!     ---------------------
!  
      read(unit_identifier_hf,*) wavefn % n_mo, wavefn % n_occ, n_lambda, wavefn % nuclear_potential, wavefn % scf_energy
!
!     Setting n_vir
!
      wavefn % n_vir = (wavefn % n_mo) - (wavefn % n_occ)

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
