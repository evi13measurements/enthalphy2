module hf_class
!
!  Date:    April 2017
!  Authors: Sarai D. Folkestad and Eirik F. Kjønstad
   use workspace
   use mlcc_types
   use input_output
!
   implicit none
!
   type hartree_fock
!
      integer(i15) :: n_occ
      integer(i15) :: n_vir
      integer(i15) :: n_mo 
      integer(i15) :: n_J
      integer(i15) :: n_ao
!
      real(dp), dimension(:,:), allocatable :: mo_coef
      real(dp), dimension(:,:), allocatable :: fock_diagonal
      real(dp), dimension(:,:), allocatable :: fock_matrix_ij
      real(dp), dimension(:,:), allocatable :: fock_matrix_ia
      real(dp), dimension(:,:), allocatable :: fock_matrix_ai
      real(dp), dimension(:,:), allocatable :: fock_matrix_ab
!
      real(dp) :: nuclear_potential
      real(dp) :: scf_energy
!
   contains
!
      procedure  :: init => init_hartree_fock
!
      procedure  :: read_cholesky_ij => read_cholesky_ij_hartree_fock
      procedure  :: read_cholesky_ia => read_cholesky_ia_hartree_fock
      procedure  :: read_cholesky_ai => read_cholesky_ai_hartree_fock
      procedure  :: read_cholesky_ab => read_cholesky_ab_hartree_fock
!
   end type hartree_fock
!
!  Private submodules and functions
!
   private :: read_hf_info
   private :: read_and_transform_to_mo_cholesky_vectors
!
contains 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!                                   !!!
!!!  Class subroutines and functions  !!!
!!!                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      call read_and_transform_to_mo_cholesky_vectors(wavefn)
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!                                    !!!
!!!  Private subroutines and functions !!!
!!!                                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
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
  subroutine read_and_transform_to_mo_cholesky_vectors(wavefn)
!
!     Read and Transform Cholesky Vectors
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 20 Apr 2017
!
!     Reads the AO Cholesky vectors from file, transforms the vectors 
!     to the MO basis, and saves the MO vectors to file
!
      use input_output
      use workspace
      use mlcc_oo_utilities
!
      implicit none
!
      class(hartree_fock) :: wavefn
!
      integer(i15) :: unit_chol_ao    = -1
      integer(i15) :: unit_chol_mo_ij = -1
      integer(i15) :: unit_chol_mo_ia = -1
      integer(i15) :: unit_chol_mo_ab = -1
!
      integer(i15) :: n_ao_sq_packed      = 0
!
      real(dp), dimension(:,:), allocatable :: chol_ao
      real(dp), dimension(:,:), allocatable :: chol_ao_sq
      real(dp), dimension(:,:), allocatable :: chol_mo_sq
!
      real(dp), dimension(:,:), allocatable :: X
!
      integer(i15) :: i,j,a,b,k 
!
!     Open Dalton file MLCC_CHOLESKY (see mlcc_write_cholesky.F)
! 
      call generate_unit_identifier(unit_chol_ao)
      open(unit=unit_chol_ao,file='mlcc_cholesky',status='old',form='formatted')
      rewind(unit_chol_ao)
!
!     Read the number of Cholesky vectors (n_J) and 
!     the number of atomic orbitals (n_ao)
!
      read(unit_chol_ao,*) wavefn%n_ao, wavefn%n_J
!
!     Open files for MO Cholesky vectors 
! 
      call generate_unit_identifier(unit_chol_mo_ij)
      call generate_unit_identifier(unit_chol_mo_ia)
      call generate_unit_identifier(unit_chol_mo_ab)
!
      open(unit_chol_mo_ij,file='cholesky_ij',status='new',form='unformatted')
      rewind(unit_chol_mo_ij)
!
      open(unit_chol_mo_ia,file='cholesky_ia',status='new',form='unformatted')
      rewind(unit_chol_mo_ia)
!
      open(unit_chol_mo_ab,file='cholesky_ab',status='new',form='unformatted')
      rewind(unit_chol_mo_ab)
!
!     Allocate packed and unpacked Cholesky AO, and 
!     unpacked Cholesky MO vectors
!
      n_ao_sq_packed = packed_size(wavefn%n_ao)
!
      call allocator(chol_ao,n_ao_sq_packed,1)
      call allocator(chol_ao_sq,wavefn%n_ao,wavefn%n_ao) 
      call allocator(chol_mo_sq,wavefn%n_mo,wavefn%n_mo)
!
      chol_ao    = zero
      chol_ao_sq = zero
      chol_mo_sq = zero
!
!     Allocate an intermediate, X
!
      call allocator(X,wavefn%n_ao,wavefn%n_mo)
!
      X = zero
!
!     Loop over the number of Cholesky vectors,
!     reading them one by one 
!
      do j = 1,wavefn%n_J
!
!        Read Cholesky AO vector
!
         read(unit_chol_ao,*) (chol_ao(i,1),i=1,n_ao_sq_packed)
!
!        Unpack/square up AO vector 
!
         call squareup(chol_ao,chol_ao_sq,wavefn%n_ao)
!
!        Transform the AO vectors to form the Cholesky MO vectors
!
         call dgemm('N','N',         &
                     wavefn%n_ao,    &
                     wavefn%n_mo,    &
                     wavefn%n_ao,    &
                     one,            &
                     chol_ao_sq,     &
                     wavefn%n_ao,    &
                     wavefn%mo_coef, &
                     wavefn%n_ao,    &
                     zero,           &
                     X,              &
                     wavefn%n_ao)
!
         call dgemm('T','N',         &
                     wavefn%n_mo,    &
                     wavefn%n_mo,    &
                     wavefn%n_ao,    &
                     one,            &
                     wavefn%mo_coef, &
                     wavefn%n_ao,    &
                     X,              &
                     wavefn%n_ao,    &
                     zero,           &
                     chol_mo_sq,     &
                     wavefn%n_mo)
!
!        Write the MO vectors to files in blocks
!
         write(unit_chol_mo_ij) ((chol_mo_sq(i,j),i=1,wavefn%n_occ),k=1,wavefn%n_occ)
         write(unit_chol_mo_ia) ((chol_mo_sq(i,a),i=1,wavefn%n_occ),a=wavefn%n_occ+1,wavefn%n_mo)
         write(unit_chol_mo_ab) ((chol_mo_sq(a,b),a=wavefn%n_occ+1,wavefn%n_mo),b=wavefn%n_occ+1,wavefn%n_mo)
!
      enddo
!
!     Close files 
!
      close(unit_chol_mo_ij)
      close(unit_chol_mo_ia)
      close(unit_chol_mo_ab)
!  
   end subroutine read_and_transform_to_mo_cholesky_vectors
end module hf_class
