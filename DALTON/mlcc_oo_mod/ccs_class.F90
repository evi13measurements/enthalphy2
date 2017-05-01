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
   type, extends(hartree_fock) :: ccs
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
      procedure :: init => init_ccs
      procedure :: drv  => drv_ccs
!
!     Initialization routine for the (singles) amplitudes
!      
      procedure :: initialize_amplitudes => initialize_amplitudes_ccs
!
!     Initialization routine for the Fock matrix, and a Fock matrix constructor
!     (for the given T1 amplitudes)
!
      procedure                  :: initialize_fock_matrix => initialize_fock_matrix_ccs
      procedure, non_overridable :: fock_constructor       => fock_constructor_ccs
!
      procedure, non_overridable :: one_electron_t1        => one_electron_t1_ccs ! T1-transf. of h_pq
!
!     get Cholesky routines to calculate the occ/vir-occ/vir
!     blocks of the T1-transformed Cholesky vectors
!
      procedure, non_overridable :: get_cholesky_ij => get_cholesky_ij_ccs ! occ-occ
      procedure, non_overridable :: get_cholesky_ia => get_cholesky_ia_ccs ! occ-vir
      procedure, non_overridable :: get_cholesky_ai => get_cholesky_ai_ccs ! vir-occ
      procedure, non_overridable :: get_cholesky_ab => get_cholesky_ab_ccs ! vir-vir
!
   end type ccs
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCS -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface
!
!
      module subroutine get_cholesky_ij_ccs(wf, L_ij_J)
!
!        Get Cholesky IJ
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Calculates the T1-transformed Cholesky vector L_ij^J,
!        and places it in L_ij_J
!
         class(ccs) :: wf
         real(dp), dimension((wf%n_o)**2, wf%n_J) :: L_ij_J
!
      end subroutine get_cholesky_ij_ccs
!
!
      module subroutine get_cholesky_ia_ccs(wf, L_ia_J)
!
!        Get Cholesky IA
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Calculates the T1-transformed Cholesky vector L_ia^J,
!        and places it in L_ia_J
!
         class(ccs) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), wf%n_J) :: L_ia_J
!
      end subroutine get_cholesky_ia_ccs
!
!
      module subroutine get_cholesky_ai_ccs(wf,L_ai_J)
!
!        Get Cholesky AI
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Calculates the T1-transformed Cholesky vector L_ai^J,
!        and places it in L_ai_J 
!
         class(ccs) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), wf%n_J) :: L_ai_J
!
      end subroutine get_cholesky_ai_ccs
!
!
      module subroutine get_cholesky_ab_ccs(wf, L_ab_J, first, last, ab_dim, reorder)
!
!        Get Cholesky AB
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Calculates the T1-transformed Cholesky vector L_ab^J,
!        and places it in L_ab_J (with options for batching over a and b)
!
         class(ccs) :: wf
!
         integer(i15), intent(in) :: ab_dim
         integer(i15), intent(in) :: first
         integer(i15), intent(in) :: last
         logical, intent(in)      :: reorder
!
         real(dp), dimension(ab_dim, wf%n_J) :: L_ab_J
!
      end subroutine get_cholesky_ab_ccs
!
      module subroutine initialize_fock_matrix_ccs(wf)
!  
!        Initialize Fock matrix
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Allocates and sets Fock matrix blocks (ij, ia, ai, ab) to zero
!        before calling the Fock matrix constructor
!
         class(ccs) :: wf
!     
      end subroutine initialize_fock_matrix_ccs
!
      module subroutine fock_constructor_ccs(wf)
!
!        Fock Constructor
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Constructs the T1-transformed Fock matrix blocks (occ/vir-occ/vir),
!        and saves the result in the class variables fock_matrix_pq (see the
!        Hartree-Fock class for these variables)  
!
         class(ccs) :: wf
!
      end subroutine fock_constructor_ccs

      module subroutine one_electron_t1_ccs(wf, h1 ,h1_T1)
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
         class(ccs) :: wf
!
         real(dp), dimension(wf%n_mo, wf%n_mo) :: h1
         real(dp), dimension(wf%n_mo, wf%n_mo) :: h1_T1
!  
      end subroutine one_electron_t1_ccs
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
   subroutine init_ccs(wf)
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
      class(ccs) :: wf
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
   end subroutine init_ccs
!
!
   subroutine drv_ccs(wf)
!
!     CCS Driver
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     If called, the routine lets the user know there is no driver 
!     for CCS, then exits the program.
!
      implicit none 
!
      class(ccs) :: wf
!
      write(unit_output,*) 'ERROR: There is no driver for the CCS class'
      call exit
!
   end subroutine drv_ccs
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine initialize_amplitudes_ccs(wf)
!
!     Initialize Amplitudes (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Allocates the singles amplitudes, sets them to zero, and calculates
!     the number of singles amplitudes.
!
      implicit none 
!
      class(ccs) :: wf
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
   end subroutine initialize_amplitudes_ccs
!
!   
end module ccs_class
