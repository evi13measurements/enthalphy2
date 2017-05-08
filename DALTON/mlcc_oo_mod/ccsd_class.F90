module ccsd_class
!
!
!
!           Coupled cluster singles and doubles (CCSD) class module                                 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017         
!                                                                           
!
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
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
!
!  ::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CCSD class -::-
!  ::::::::::::::::::::::::::::::::::::::
!
!
   type, extends(ccs) :: ccsd
!
!     Amplitude attributes
!
      integer(i15) :: n_t2am = 0                    ! Number of doubles amplitudes
      real(dp), dimension(:,:), allocatable :: t2am ! Doubles amplitude vector
!
!     Schrödinger equation projection vector (the omega vector)
! 
!        < mu | exp(-T) H exp(T) | R >
!
      real(dp), dimension(:,:), allocatable :: omega2 ! Doubles vector
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_ccsd
      procedure :: drv  => drv_ccsd
!
!     Initialization routine for the (singles, doubles) amplitudes
!
      procedure :: initialize_amplitudes => initialize_amplitudes_ccsd
!
!     Routine to calculate the energy from the current amplitudes
!
      procedure :: calc_energy => calc_energy_ccsd 
!
!     Routine to initialize omega (allocate and set to zero)
!
      procedure :: initialize_omega => initialize_omega_ccsd
!
!     Routines to construct the projection vector (omega)
!
      procedure :: construct_omega => construct_omega_ccsd
!
!     Helper routines for construct_omega
!
      procedure :: omega_a1 => omega_a1_ccsd 
      procedure :: omega_b1 => omega_b1_ccsd 
      procedure :: omega_c1 => omega_c1_ccsd 
      procedure :: omega_d1 => omega_d1_ccsd
!
      procedure :: omega_a2 => omega_a2_ccsd 
      procedure :: omega_b2 => omega_b2_ccsd 
      procedure :: omega_c2 => omega_c2_ccsd 
      procedure :: omega_d2 => omega_d2_ccsd 
      procedure :: omega_e2 => omega_e2_ccsd       
!
!     Ground state solver routine (helpers only, see CCS for the rest)
!
      procedure :: calc_ampeqs_norm          => calc_ampeqs_norm_ccsd
      procedure :: new_amplitudes            => new_amplitudes_ccsd
      procedure :: calc_quasi_Newton_doubles => calc_quasi_Newton_doubles_ccsd
!
!     Jacobian transformation routines 
!
      procedure :: jacobian_transformation => jacobian_transformation_ccsd
!
      procedure :: jacobian_a1 => jacobian_a1_ccsd
      procedure :: jacobian_b1 => jacobian_b1_ccsd
!
   end type ccsd
!
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCSD -::- 
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
   interface
!
!
      module subroutine initialize_omega_ccsd(wf)
!
!        Initialize Omega (CCSD)
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!        Allocates the projection vector (omega1, omega2) and sets it
!        to zero.
!
         class(ccsd) :: wf
!
      end subroutine initialize_omega_ccsd
!
!
      module subroutine construct_omega_ccsd(wf)
!
!        Construct Omega (CCSD)
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!        for the current amplitudes of the object wfn 
!
         class(ccsd) :: wf 
!
      end subroutine construct_omega_ccsd
!
!
      module subroutine omega_a1_ccsd(wf)
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
         class(ccsd) :: wf
!
      end subroutine omega_a1_ccsd
!
!
      module subroutine omega_b1_ccsd(wf)
!
!        Omega B1
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!        Calculates the B1 term, 
!
!          - sum_ckl u_kl^ac * g_kilc,
! 
!        and adds it to the singles projection vector (omeg1) of
!        the wavefunction object wfn
!
         class(ccsd) :: wf 
!
      end subroutine omega_b1_ccsd
!
!
      module subroutine omega_c1_ccsd(wf)
!
!        C1 omega term: Omega_ai^C1 = sum_ck F_kc*u_ai_ck,
!                       u_ai_ck = 2*t_ck_ai-t_ci_ak
!        
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, March 2017
!
         class(ccsd) :: wf 
!
      end subroutine omega_c1_ccsd
!
!
      module subroutine omega_d1_ccsd(wf)
!
!        D1 omega term: Omega_ai^D1=F_ai_T1
!
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, March 2017
!
            class(ccsd) :: wf
!
      end subroutine omega_d1_ccsd
!
!
      module subroutine omega_a2_ccsd(wf)
!
!        MLCC Omega A2 term: Omega A2 = g_ai_bj + sum_(cd)g_ac_bd * t_ci_dj = A2.1 + A.2.2
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, 10 Mar 2017
!
!        Structure: Batching over both a and b. If no batching is necessary L_ab_J is only read once, and g_ac_bd 
!                   is constructed and kept in memory full size. 
!                   g_ac_bd is reordered as g_ab_cd and t_ci_dj is reordered as t_cd_ij.
!                   Omega contribution for A2.2 is ordered as Omega_ab_ij, and is reordered into the packed omega2 vector.          
!
!
         class(ccsd) :: wf
!
      end subroutine omega_a2_ccsd
!
!
      module subroutine omega_b2_ccsd(wf)
!
!        MLCC Omega B2 term.  Omega B2 = sum_(kl) t_ak_bl*(g_kilj + sum_(cd) t_ci_dj * g_kc_ld)
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, 11 Mar 2017
! 
!
!        Structure: g_kilj is constructed first and reordered as g_kl_ij. 
!                   Then the contraction over cd is performed, and the results added to g_kl_ij.
!                   t_ak_bl is then reordered as t_ab_kl and the contraction over kl is performed.
!
         class(ccsd) :: wf
!
      end subroutine omega_b2_ccsd
!
!
      module subroutine omega_c2_ccsd(wf)
!
!        C2 omega term. Omega C2 = -1/2* sum_(ck)t_bk_cj*(g_ki_ac -1/2 sum_(dl)t_al_di * g_kd_lc)
!                                  - sum_(ck) t_bk_ci (g_kj_ac-sum_(dl)t_al_dj*g_kd_lc)
!
!
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!     
!
         class(ccsd) :: wf
!
      end subroutine omega_c2_ccsd
!
!
      module subroutine omega_d2_ccsd(wf)
!
         class(ccsd) :: wf
!
      end subroutine omega_d2_ccsd
!
!
      module subroutine omega_e2_ccsd(wf)
!
         class(ccsd) :: wf
!
      end subroutine omega_e2_ccsd
!
!
      module subroutine calc_ampeqs_norm_ccsd(wf, ampeqs_norm)
!
!        Calculate Amplitude Equations Norm (CCSD)
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
         class(ccsd) :: wf 
!
         real(dp) :: ampeqs_norm 
!
      end subroutine calc_ampeqs_norm_ccsd
!
!
      module subroutine new_amplitudes_ccsd(wf)
!
!        New Amplitudes (CCSD)
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!        Directs the calculation of the quasi-Newton estimate Δ t_i, 
!        and t_i + Δ t_i, and calls the DIIS routine to save & get 
!        the amplitudes for the next iteration. 
!
         class(ccsd) :: wf 
!
      end subroutine new_amplitudes_ccsd
!
!
      module subroutine calc_quasi_Newton_doubles_ccsd(wf,dt,n_variables)
!
!        Calculate quasi-Newton estimate (CCSD)
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!        Calculates the quasi-Newton estimate Δ t_i (doubbles part)
!        and places the contribution in the dt vector (of length n_variables,
!        with singles first, then doubles, etc. if inherited)
!
         class(ccsd) :: wf 
!
         integer(i15), intent(in) :: n_variables
         real(dp), dimension(n_variables, 1) :: dt
!
      end subroutine calc_quasi_Newton_doubles_ccsd
!
!
      module subroutine jacobian_transformation_ccsd(wf,c1am,c2am)
!
!        Jacobian Transformation (CCSD)
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!
!        Directs the transformation of the incoming vector c by the 
!        coupled cluster Jacobian matrix 
!
!           A_mu,nu = < mu | [e^(-T) H e^(T),tau_nu] | R >.
!
!        On exit, A*c is placed in the incoming c vector.
!
         class(ccsd) :: wf 
!
         real(dp), dimension(:,:) :: c1am 
         real(dp), dimension(:,:) :: c2am 
!
      end subroutine jacobian_transformation_ccsd
!
!
      module subroutine jacobian_a1_ccsd(wf,tr1am,c1am)
!
!        Jacobian A1 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!        Calculates the A1 term, 
!
!           sum_c F_ac c_ci - sum_k c_ak F_ki + sum_ck L_aikc c_ck
!
!        and adds it to the transformed singles vector element tr1am(a,i).
!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: tr1am 
!
         real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1am 
!
      end subroutine jacobian_a1_ccsd
!
      module subroutine jacobian_b1_ccsd(wf,tr1am,c1am)
!
!        Jacobian B1
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!        Calculates the B1 term,
!
!           sum_dl (sum_ck L_ldkc c_ck) u_li^da,
!
!        where L_ldkc = 2 * g_ldkc - g_lckd and 
!        u_li^ad = 2 * t_li^ad - t_il^ad.
!
         class(ccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: tr1am 
!
         real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1am
!
      end subroutine jacobian_b1_ccsd
!
   end interface
!
!
contains
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::-
!  ::::::::::::::::::::::::::::::::::::::::::::
!
!
   subroutine init_ccsd(wf)
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
      class(ccsd) :: wf
!
!     Set model name 
!
      wf%name = 'CCSD   '
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
!     Initialize the projection vector (omega)
!
      call wf%initialize_omega
!
   end subroutine init_ccsd
!
!
   subroutine drv_ccsd(wf)
!
      implicit none 
!
      class(ccsd) :: wf
!
      call wf%ground_state_solver
!
      call wf%jacobian_transformation(wf%omega1,wf%omega2)
!
   end subroutine drv_ccsd
!
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
!
   subroutine initialize_amplitudes_ccsd(wf)
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
      class(ccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_ia_jb
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0
      integer(i15) :: ai = 0, bj = 0, ia = 0, jb = 0, aibj = 0
!
!     Calculate the number of singles and doubles amplitudes
!
      wf%n_t1am = (wf%n_o)*(wf%n_v) 
      wf%n_t2am = (wf%n_t1am)*(wf%n_t1am + 1)/2
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
   end subroutine initialize_amplitudes_ccsd
!
!
   subroutine calc_energy_ccsd(wf)
!
!     Calculate Energy (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the CCSD energy
!
      implicit none 
!
      class(ccsd) :: wf 
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
                  wf%energy = wf%energy +                                           & 
                                 (wf%t2am(aibj,1) + (wf%t1am(a,i))*(wf%t1am(b,j)))* &
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
   end subroutine calc_energy_ccsd
!
!
end module ccsd_class
