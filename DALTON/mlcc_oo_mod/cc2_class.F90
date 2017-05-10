module cc2_class
!
!
!
!            Coupled cluster perturbative doubles (CC2) class module                                 
!         Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017         
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
!  The ancestor class module (CCS)
!
   use ccs_class
!
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CC2 class -::-
!  ::::::::::::::::::::::::::::::::::::: 
!
   type, extends(ccs) :: cc2
!
   contains 
!
!     Initialization and driver routines
!
      procedure :: init => init_cc2
      procedure :: drv  => drv_cc2
!
!     Routines to construct the projection vector (omega)
!
      procedure :: construct_omega => construct_omega_cc2
!
!     Helper routines for construct_omega
!
      procedure :: omega_a1 => omega_a1_cc2 
      procedure :: omega_b1 => omega_b1_cc2 
      procedure :: omega_c1 => omega_c1_cc2 
      procedure :: omega_d1 => omega_d1_cc2      
!
!     Ground state solver helper routines
!
      procedure :: calc_energy => calc_energy_cc2
!
   end type cc2
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCS -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface
!
      module subroutine construct_omega_cc2(wf)
!
!        Construct Omega 
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!
!        Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!        for the current amplitudes of the object wf 
!
         class(cc2) :: wf 
!
      end subroutine construct_omega_cc2
!
      module subroutine omega_a1_cc2(wf, t_kc_di, c_first, c_last, c_length)
!
!        Omega A1 term
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!  
!        Calculates the A1 term, 
!  
!           sum_ckd g_adkc * u_ki^cd,
!  
!        and adds it to the singles projection vector (omeg1) of
!        the wavefunction object wfn
!
         class(cc2) :: wf
!  
         integer(i15) :: c_first, c_last, c_length
!
         real(dp), dimension(c_length*(wf%n_o),(wf%n_v)*(wf%n_o)):: t_kc_di
!
      end subroutine omega_a1_cc2
!
!
      module subroutine omega_b1_cc2(wf, t_lc_ak, c_first, c_last, c_length)
!
!        Omega B1
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!
!        Calculates the B1 term, 
!
!          - sum_ckl u_kl^ac * g_kilc,
! 
!        and adds it to the singles projection vector (omeg1) of
!        the wavefunction object wfn
!
         class(cc2) :: wf 
!
         integer(i15) :: c_first, c_last, c_length
!
         real(dp), dimension(c_length*(wf%n_o),(wf%n_v)*(wf%n_o)) :: t_lc_ak
!
      end subroutine omega_b1_cc2
!
!
      module subroutine omega_c1_cc2(wf, t_kc_ai, c_first, c_last, c_length)
!
!        C1 omega term: Omega_ai^C1 = sum_ck F_kc*u_ai_ck,
!                       u_ai_ck = 2*t_ck_ai-t_ci_ak
!        
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, March 2017
!
         class(cc2) :: wf 
!
         integer(i15) :: c_first, c_last, c_length
!
         real(dp), dimension(c_length*(wf%n_o),(wf%n_v)*(wf%n_o)) :: t_kc_ai
!
      end subroutine omega_c1_cc2
!
!
      module subroutine omega_d1_cc2(wf)
!
!        D1 omega term: Omega_ai^D1=F_ai_T1
!
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, March 2017
!
         class(cc2) :: wf
!
      end subroutine omega_d1_cc2
!
   end interface
!
contains
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::
!
   subroutine init_cc2(wf)
!
!     Initialize CC2 object
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!
!     Performs the following tasks
!
!     1. Sets HF orbital and energy information by reading from file (read_hf_info)
!     2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!     3. Allocates the Fock matrix and sets it to zero (note: the matrix is constructed in 
!        the descendant classes) 
!     4. Allocates the singles amplitudes and sets them to zero, and sets associated properties 
!     5. Allocate Omega vector
!
      implicit none
!
      class(cc2)  :: wf
!
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
!     Initialize omega vector
!
      call wf%initialize_omega
!
   end subroutine init_cc2
!
   subroutine drv_cc2(wf)
!
!     CC2 Driver
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!
      implicit none 
!
      class(cc2) :: wf
!
      call wf%ground_state_solver
!
   end subroutine drv_cc2
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine calc_energy_cc2(wf)
!
!  Calculate Energy (CC2)
!
!  Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!  Calculate the CC2 energy
!
   implicit none
!
   class(cc2) :: wf
!
!  Integrals
!
   real(dp), dimension(:,:), allocatable :: L_bj_J
   real(dp), dimension(:,:), allocatable :: L_ia_J
   real(dp), dimension(:,:), allocatable :: g_ia_bj ! = g_aibj
!
!  t2 amplitudes
!
   real(dp), dimension(:,:), allocatable :: t_ia_bj ! = g_aibj/(e_a + e_b - e_i - e_j)
!
!  Batching variables
!  
   integer(i15) :: a_batch, a_first, a_last, a_length
   integer(i15) :: required, available, n_batch, batch_dimension, max_batch_length
!
!  Indices
!  
   integer(i15) :: a = 0, b = 0
   integer(i15) :: i = 0, j = 0
!
   integer(i15) :: bj = 0, ai = 0, ia = 0
!
!  Prepare for batching over index a 
!  
   required = (2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                     !    
            + 2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                      ! Needed for g_aibj  
            + 2*((wf%n_v)**2)*(wf%n_J) + ((wf%n_o)**2)*(wf%n_J) & ! and 't2' amplitudes  
            + 2*(wf%n_v)**2*(wf%n_o)**2 &                         !
            + (wf%n_v)*(wf%n_o)*(wf%n_J) &                        ! L_ia_J
            + (wf%n_v)*(wf%n_o)**2)                               ! g_ia_jb
!     
   required = 4*required ! In words
   available = get_available()
!
   batch_dimension  = wf%n_v ! Batch over the virtual index a
   max_batch_length = 0      ! Initilization of unset variables 
   n_batch          = 0
!
   call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!  Loop over the number of a batches 
!
   do a_batch = 1, n_batch
!
!     For each batch, get the limits for the a index 
!
      call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
      a_length = a_last - a_first + 1 
!
!     Allocate L_bj_J and L_ia_J (= reordering of L_bj_J constrained to the batch)
!
      call allocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
      call allocator(L_ia_J, a_length*(wf%n_o), wf%n_J)
      call dzero(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
      call dzero(L_ia_J, a_length*(wf%n_o), wf%n_J)
!
      call wf%get_cholesky_ai(L_bj_J)
!
!     Create L_ia_J
!
      do a = 1, a_length
         do i = 1, wf%n_v
            do J = 1, wf%n_J
!
!              Calculate compound indices
!
               ia = index_two(i, a, wf%n_o)
               ai = index_two(a + a_first - 1, i, wf%n_v)
!
               L_ia_J(ia, J) = L_bj_J(ai, J)
!
            enddo
         enddo
      enddo
!
!     Allocate g_ia_bj
!
      call allocator(g_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!     Construct integral g_ia_bj (= g_aibj for the batch)
!
      call dgemm('N','T',            &
                  a_length*(wf%n_o), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  a_length*(wf%n_o), &
                  L_bj_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ia_bj,           &
                  a_length*(wf%n_o))
!
!     L_bj_J and L_ia_J
!
      call deallocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
      call deallocator(L_ia_J, a_length*(wf%n_o), wf%n_J)
!
!     Allocate t_ia_bj
!
      call allocator(t_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!     Create t2 amplitudes
!
      do i = 1, wf%n_o
         do a = 1, a_length
!
!           Calculate compound index
!
            ia = index_two(i, a, wf%n_o)
!
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculare compond index               
!
                  bj = index_two(b, j, wf%n_v)
!
                  t_ia_bj(ia, bj) = - g_ia_bj(ia, bj)/(wf%fock_diagonal(a + wf%n_o, 1) &
                                                     + wf%fock_diagonal(b + wf%n_o, 1) &
                                                     - wf%fock_diagonal(i, 1) - wf%fock_diagonal(j, 1))
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ia_bj
!
      call deallocator(g_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!     Deallocate t_ia_bj
!
      call deallocator(t_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
   enddo
!
   end subroutine calc_energy_cc2
!
end module cc2_class