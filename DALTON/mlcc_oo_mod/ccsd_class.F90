module ccsd_class
!
!
!
!           Coupled cluster singles and doubles (CCSD) class module                                 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017         
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
!  Ancestor class module (CCS)
!
   use ccs_class
!
   implicit none 
!
!  ::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CCSD class -::-
!  ::::::::::::::::::::::::::::::::::::::
!
   type, extends(cc_singles) :: cc_singles_doubles
!
!     Amplitude attributes
!
      integer(i15) :: n_t2am = 0                    ! Number of doubles amplitudes
      real(dp), dimension(:,:), allocatable :: t2am ! Doubles amplitude vector
!
!     Schrödinger equation projection attribbutes (the omega vector)
! 
!        < mu | exp(-T) H exp(T) | R >
!
      real(dp), dimension(:,:), allocatable :: omega1 ! Singles part
      real(dp), dimension(:,:), allocatable :: omega2 ! Doubles part
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_cc_singles_doubles
      procedure :: drv  => drv_cc_singles_doubles
!
!     Initialization routine for the (singles, doubles) amplitudes
!
      procedure :: initialize_amplitudes => initialize_amplitudes_cc_singles_doubles
!
!     Routine to calculate the energy from the current amplitudes
!
      procedure :: calc_energy => calc_energy_cc_singles_doubles 
!
!     Routines to construct the projection vector (omega)
!
      procedure :: construct_omega => construct_omega_cc_singles_doubles
!
!     Helper routines for construct_omega
!
      procedure :: omega_a1 => omega_a1_cc_singles_doubles
  !    procedure :: omega_b1 => omega_b1_cc_singles_doubles
!
   end type cc_singles_doubles
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCSD -::- 
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface
!
!
      module subroutine construct_omega_cc_singles_doubles(wf)
!
!        Construct Omega 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!        for the current amplitudes of the object wfn 
!
         class(cc_singles_doubles) :: wf 
!
      end subroutine construct_omega_cc_singles_doubles
!
!
      module subroutine omega_a1_cc_singles_doubles(wf)
!
!        Omega A1 term
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!  
!        Calculates the A1 term, 
!  
!           sum_ckd g_adkc * u_ki^cd,
!  
!        and adds it to the singles projection vector (omeg1) of
!        the wavefunction object wfn
!
         class(cc_singles_doubles) :: wf
!
      end subroutine omega_a1_cc_singles_doubles
!
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
   subroutine init_cc_singles_doubles(wf)
!
!     Initialize CCSD object
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Performs the following tasks:
!
!        1. Sets HF orbital and energy information by reading from file (read_hf_info)
!        2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!        3. Allocates the Fock matrix and sets it to zero
!        4. Initializes the amplitudes (sets their initial values and associated variables)
!        5. Sets the initial energy based on the initial amplitudes (in particular, the MP2
!           estimate of the doubles amplitude)
!
      implicit none 
!
      class(cc_singles_doubles) :: wf
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky 
!
!     Initialize (singles and doubles) amplitudes
!
      call wf%initialize_amplitudes
!
!     Initialize the Fock matrix (allocate and construct given the initial amplitudes)
!
      call wf%initialize_fock_matrix
!
!     Set the initial value of the energy (given the initial amplitudes) 
!
      call wf%calc_energy
!
   end subroutine init_cc_singles_doubles
!
!
   subroutine drv_cc_singles_doubles(wf)
!
      implicit none 
!
      class(cc_singles_doubles) :: wf
!
!     Print the energy    
!
      write(unit_output,*) 'The SCF energy:', wf%scf_energy
      write(unit_output,*) 'The MP2 energy:', wf%energy
!
   end subroutine drv_cc_singles_doubles
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine initialize_amplitudes_cc_singles_doubles(wf)
!
!     Initialize Amplitudes (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Allocates the amplitudes, sets them to zero, calculates
!     the number of amplitudes, and sets the doubles amplitudes
!     to the perturbative MP2 estimate.
!
      implicit none 
!
      class(cc_singles_doubles) :: wf
!
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_ia_jb
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0
      integer(i15) :: ai = 0, bj = 0, ia = 0, jb = 0, aibj = 0
!
!     Calculate the number of singles and doubles amplitudes
!
      wf%n_t1am = (wf%n_o) * (wf%n_v) 
      wf%n_t2am = (wf%n_t1am) * (wf%n_t1am + 1)/2
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator(wf%t1am, wf%n_v, wf%n_o)
      wf%t1am = zero
!
!     Allocate the doubles amplitudes and set to zero
!
      call allocator(wf%t2am, wf%n_t2am, 1)
      wf%t2am = zero
!
!
!     -::- Initialize the doubles amplitudes to the MP2 estimate -::-
!
!
!     Allocate L_ia_J and g_ia_jb
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call allocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      L_ia_J = zero
      g_ia_jb = zero
!
!     Get the Cholesky IA vector 
!
      call wf%get_cholesky_ia(L_ia_J)
!
!     Calculate g_ia_jb = g_iajb
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ia_jb,           &
                  (wf%n_o)*(wf%n_v))
!
!     Set the doubles amplitudes
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do j = 1, wf%n_o
!
!                 Get necessary indices
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
                  ia = index_two(i, a, wf%n_o)
                  jb = index_two(j, b, wf%n_o)
!
!                 Set the doubles indices
!
                  if (ai .le. bj) then ! To avoid setting the same element twice
!
                     aibj = index_packed(ai,bj)
!
                     wf%t2am(aibj, 1) = - g_ia_jb(ia,jb)/(wf%fock_diagonal(wf%n_o + a, 1) + &
                                                            wf%fock_diagonal(wf%n_o + b, 1) - &
                                                            wf%fock_diagonal(i, 1) - &
                                                            wf%fock_diagonal(j, 1))
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocations
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), (wf%n_J))
      call deallocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) 
!
   end subroutine initialize_amplitudes_cc_singles_doubles
!
!
   subroutine calc_energy_cc_singles_doubles(wf)
!
!     Calculate Energy (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the CCSD energy
!
      implicit none 
!
      class(cc_singles_doubles) :: wf 
!
      real(dp), dimension(:,:), allocatable :: L_ia_J  ! L_ia^J
      real(dp), dimension(:,:), allocatable :: g_ia_jb ! g_iajb
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0, ai = 0
      integer(i15) :: bj = 0, aibj = 0, ia = 0, jb = 0, ib = 0, ja = 0
!
!     Allocate the Cholesky vector L_ia_J = L_ia^J and set to zero 
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_ia_J = zero
!
!     Get the Cholesky vector L_ia_J 
!
      call wf%get_cholesky_ia(L_ia_J)
!
!     Allocate g_ia_jb = g_iajb and set it to zero
!
      call allocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ia_jb = zero
!
!     Calculate the integrals g_ia_jb from the Cholesky vector L_ia_J 
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ia_jb,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate the Cholesky vector L_ia_J 
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Set the initial value of the energy 
!
      wf%energy = wf%scf_energy
!
!     Add the correlation energy E = E + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a, i, wf%n_v)
                  bj   = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  ia   = index_two(i, a, wf%n_o)
                  jb   = index_two(j, b, wf%n_o)
!
                  ib   = index_two(i, b, wf%n_o)
                  ja   = index_two(j, a, wf%n_o)
!
!                 Add the correlation energy 
!
                  wf%energy = wf%energy + & 
                                 (wf%t2am(aibj,1) + (wf%t1am(a,i))*(wf%t1am(b,j)))*&
                                 (two*g_ia_jb(ia,jb) - g_ia_jb(ib,ja))
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ia_jb
!
      call deallocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine calc_energy_cc_singles_doubles
!
!
end module ccsd_class
