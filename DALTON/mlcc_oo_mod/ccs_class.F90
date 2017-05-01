module ccs_class
!
!
!
!               Coupled cluster singles (CCS) class module                                 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017         
!                                                                           
!
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
!  General tools
!
   use types
   use utils
   use workspace
   use input_output
!
!  The ancestor class module (HF)
!
   use hf_class
!
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CCS class -::-
!  ::::::::::::::::::::::::::::::::::::: 
!
   type, extends(hartree_fock) :: cc_singles
!
!     Amplitude attributes
!
      integer(i15) :: n_t1am = 0                    ! Number of singles amplitudes
      real(dp), dimension(:,:), allocatable :: t1am ! Singles amplitude vector
!
      real(dp), dimension(:,:), allocatable :: fock_matrix_ij ! occ-occ block
      real(dp), dimension(:,:), allocatable :: fock_matrix_ia ! occ-vir block
      real(dp), dimension(:,:), allocatable :: fock_matrix_ai ! vir-occ block
      real(dp), dimension(:,:), allocatable :: fock_matrix_ab ! vir-vir block
!
!
   contains 
!
!     Initialization and driver routines
!
      procedure :: init => init_cc_singles
      procedure :: drv  => drv_cc_singles
!
!     Initialization routine for the (singles) amplitudes
!      
      procedure :: initialize_amplitudes => initialize_amplitudes_cc_singles
!
!     Initialization routine for the Fock matrix, and a Fock matrix constructor
!     (for the given T1 amplitudes)
!
      procedure                  :: initialize_fock_matrix => initialize_fock_matrix_cc_singles
      procedure, non_overridable :: fock_constructor       => fock_constructor_cc_singles
!
      procedure, non_overridable :: one_electron_t1        => one_electron_t1_cc_singles ! T1-transf. of h_pq
!
!     get Cholesky routines to calculate the occ/vir-occ/vir
!     blocks of the T1-transformed Cholesky vectors
!
      procedure, non_overridable :: get_cholesky_ij => get_cholesky_ij_cc_singles ! occ-occ
      procedure, non_overridable :: get_cholesky_ia => get_cholesky_ia_cc_singles ! occ-vir
      procedure, non_overridable :: get_cholesky_ai => get_cholesky_ai_cc_singles ! vir-occ
      procedure, non_overridable :: get_cholesky_ab => get_cholesky_ab_cc_singles ! vir-vir
!
   end type cc_singles
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCS -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface
!
!
      module subroutine get_cholesky_ij_cc_singles(wf, L_ij_J)
!
!        Get Cholesky IJ
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Calculates the T1-transformed Cholesky vector L_ij^J,
!        and places it in L_ij_J
!
         class(cc_singles) :: wf
         real(dp), dimension((wf%n_o)**2, wf%n_J) :: L_ij_J
!
      end subroutine get_cholesky_ij_cc_singles
!
!
      module subroutine get_cholesky_ia_cc_singles(wf, L_ia_J)
!
!        Get Cholesky IA
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Calculates the T1-transformed Cholesky vector L_ia^J,
!        and places it in L_ia_J
!
         class(cc_singles) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), wf%n_J) :: L_ia_J
!
      end subroutine get_cholesky_ia_cc_singles
!
!
      module subroutine get_cholesky_ai_cc_singles(wf,L_ai_J)
!
!        Get Cholesky AI
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Calculates the T1-transformed Cholesky vector L_ai^J,
!        and places it in L_ai_J 
!
         class(cc_singles) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), wf%n_J) :: L_ai_J
!
      end subroutine get_cholesky_ai_cc_singles
!
!
      module subroutine get_cholesky_ab_cc_singles(wf, L_ab_J, first, last, ab_dim, reorder)
!
!        Get Cholesky AB
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Calculates the T1-transformed Cholesky vector L_ab^J,
!        and places it in L_ab_J (with options for batching over a and b)
!
         class(cc_singles) :: wf
!
         integer(i15), intent(in) :: ab_dim
         integer(i15), intent(in) :: first
         integer(i15), intent(in) :: last
         logical, intent(in)      :: reorder
!
         real(dp), dimension(ab_dim, wf%n_J) :: L_ab_J
!
      end subroutine get_cholesky_ab_cc_singles
!
      module subroutine initialize_fock_matrix_cc_singles(wf)
!  
!        Initialize Fock matrix
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Allocates and sets Fock matrix blocks (ij, ia, ai, ab) to zero
!        before calling the Fock matrix constructor
!
         class(cc_singles) :: wf
!     
      end subroutine
!
      module subroutine fock_constructor_cc_singles(wf)
!
!        Fock Constructor
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Constructs the T1-transformed Fock matrix blocks (occ/vir-occ/vir),
!        and saves the result in the class variables fock_matrix_pq (see the
!        Hartree-Fock class for these variables)  
!
         implicit none
!
         class(cc_singles) :: wf
!
         real(dp), dimension(:,:), allocatable :: fock_ao
         real(dp), dimension(:,:), allocatable :: fock_matrix
!
         real(dp), dimension(:,:), allocatable :: h1ao ! AO basis matrix h_αβ
         real(dp), dimension(:,:), allocatable :: h1mo ! MO basis matrix h_pq  
         real(dp), dimension(:,:), allocatable :: X    ! An intermediate
!         
         integer(i15) :: unit_identifier_ao_integrals = -1 ! Unit identifier for file mlcc_aoint
!
!        Indices
!
         integer(i15) :: i = 0, j = 0, k = 0, a = 0, b = 0, ij = 0, kk = 0, ik = 0
         integer(i15) :: kj = 0, ii = 0, jj = 0, ji = 0, ai = 0, ib = 0, bi = 0, ia = 0
         integer(i15) :: aj = 0, ja = 0, ab = 0
!
!        Useful orbital information
!         
         integer(i15) :: n_ao_sq_packed = 0 ! Dimension of packed (n_ao x n_ao) matrix
!
!        Two electron integrals
!
         real(dp), dimension(:,:), allocatable :: g_ij_kl
         real(dp), dimension(:,:), allocatable :: g_ab_ij
         real(dp), dimension(:,:), allocatable :: g_ai_jb
         real(dp), dimension(:,:), allocatable :: g_ia_jk
         real(dp), dimension(:,:), allocatable :: g_ai_jk
!
!        Cholesky vectors
!
         real(dp), dimension(:,:), allocatable :: L_ij_J 
         real(dp), dimension(:,:), allocatable :: L_ia_J 
         real(dp), dimension(:,:), allocatable :: L_ai_J 
         real(dp), dimension(:,:), allocatable :: L_ab_J 
!
!        Batch settings
!
         integer(i15) :: available = 0, required = 0, max_batch_length = 0
         integer(i15) :: batch_end = 0, batch_length = 0, g_off = 0, n_batches = 0
         integer(i15) :: b_batch = 0, batch_start = 0
!
      end subroutine

      module subroutine one_electron_t1_cc_singles(wf, h1 ,h1_T1)
!
!        One-electron T1 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        T1-transforms the one-electron MO integrals h_pq
!
!              h_p_q_T1 = sum_st x_p_s * y_q_t * h_s_t
!
!              x = I - t1
!              y = I - t1^T
!
         implicit none
!
         class(cc_singles) :: wf
!
         real(dp), dimension(wf%n_mo, wf%n_mo) :: h1
         real(dp), dimension(wf%n_mo, wf%n_mo) :: h1_T1
!
         real(dp), dimension(:,:), allocatable :: x 
         real(dp), dimension(:,:), allocatable :: y 
         real(dp), dimension(:,:), allocatable :: t1
!
         real(dp), dimension(:,:), allocatable :: Z ! Intermediate for matrix multiplication
!
         integer(i15) :: p = 0, q = 0, a = 0, i = 0
!  
      end subroutine
!
   end interface 
!
!
contains
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::
!
   subroutine init_cc_singles(wf)
!
!     Initialize CCS object
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Performs the following tasks
!
!     1. Sets HF orbital and energy information by reading from file (read_hf_info)
!     2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!     3. Allocates the Fock matrix and sets it to zero (note: the matrix is constructed in 
!        the descendant classes) 
!     4. Allocates the singles amplitudes and sets them to zero, and sets associated properties 
!
      implicit none 
!
      class(cc_singles) :: wf
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
!
!     Initialize amplitudes and associated attributes
!
      call wf%initialize_amplitudes
!
!     Allocate Fock matrix and set to zero
!
      call wf%initialize_fock_matrix
!
   end subroutine init_cc_singles
!
!
   subroutine drv_cc_singles(wf)
!
!     CCS Driver
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     If called, the routine lets the user know there is no driver 
!     for CCS, then exits the program.
!
      implicit none 
!
      class(cc_singles) :: wf
!
      write(unit_output,*) 'ERROR: There is no driver for the CCS class'
      call exit
!
   end subroutine drv_cc_singles
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine initialize_amplitudes_cc_singles(wf)
!
!     Initialize Amplitudes (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Allocates the singles amplitudes, sets them to zero, and calculates
!     the number of singles amplitudes.
!
      implicit none 
!
      class(cc_singles) :: wf
!
!     Calculate the number of singles amplitudes
!
      wf%n_t1am = (wf%n_o)*(wf%n_v) 
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator(wf%t1am, wf%n_v, wf%n_o)
      wf%t1am = zero
!
   end subroutine initialize_amplitudes_cc_singles
!
!   
end module ccs_class
