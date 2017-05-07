submodule (cc2_class) omega
!
!
!                       -::- Omega submodule (CC2) -::-
!           Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!
!
!     Contains the following family of procedures of the CC2 class:
!
!        initialize_omega: allocates the projection vector omega1
!                          and sets it to zero.
!
!        construct_omega:  constructs the projection vector omega1
!                          for the current amplitudes t1am for the
!                          wavefunction object wf. The routine assumes that
!                          the projection vector is allocated.
!
!        omega_a1:         adds A1 term to omega1
!        omega_b1:         adds B1 term to omega1
!        omega_c1:         adds C1 term to omega1
!        omega_d1:         adds D1 term to omega1
!
   implicit none 
!
   logical :: debug = .false.
!
!
contains
!
   subroutine initialize_omega_cc2(wf)
!
!     Initialize Omega (CCSD)
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!
!     Allocates the projection vector omega1 and sets it
!     to zero.
!
      implicit none
!
      class(cc2) :: wf
!
      call allocator(wf%omega1, wf%n_v, wf%n_o)
      wf%omega1 = zero
!
   end subroutine initialize_omega_cc2
!
  subroutine construct_omega_ccsd(wf)
!
!     Construct Omega (CC2)
!     Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!
!     Directs the calculation of the omega vector
!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp) :: omega_start = zero
      real(dp) :: omega_end = zero
!
      call cpu_time(omega_start)
!
!     Set the omega vector to zero 
!
      wf%omega1 = zero
!
!     Construct singles contributions 

      call wf%omega_a1
!
      call wf%omega_b1
!
      call wf%omega_c1

      call wf%omega_d1
!
      call cpu_time(omega_end)
      write(unit_output,*)'Time in omega:', omega_end-omega_start    
!
   end subroutine construct_omega_cc2
!
!
   subroutine omega_a1_cc2(wf)
!
!     Omega A1 term
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the A1 term, 
!
!        sum_ckd g_adkc * u_ki^cd,
!
!     and adds it to the singles projection vector (omeg1) of
!     the wavefunction object wfn
!
      implicit none
!
      class(ccsd) :: wf
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_batch_length = 0, batch_dimension = 0
      integer(i15) :: n_batch = 0, a_begin = 0, a_end = 0, a_batch = 0, batch_length = 0
!
!     Indices 
!
      integer(i15) :: i = 0, j = 0, k = 0, a = 0, c = 0, d = 0
!
      integer(i15) :: ad = 0, da = 0
      integer(i15) :: ci = 0, ck = 0, di = 0, dk = 0
      integer(i15) :: kc = 0
!
      integer(i15) :: ckd = 0 
      integer(i15) :: cidk = 0, ckdi = 0
!
      integer(i15) :: ad_dim = 0   
!
      real(dp), dimension(:,:), allocatable :: L_kc_J 
      real(dp), dimension(:,:), allocatable :: L_da_J  ! L_ad^J; a is being batched over
      real(dp), dimension(:,:), allocatable :: g_da_kc ! g_adkc; a is being batched over
      real(dp), dimension(:,:), allocatable :: g_a_ckd ! reordered g_adkc; a is being batched over
      real(dp), dimension(:,:), allocatable :: u_ckd_i ! u_ki^cd
!
      logical :: reorder ! To get L_ab_J reordered, for batching over a (first index)
!
   end subroutine omega_a1_ccsd
!
!
   subroutine omega_b1_ccsd(wf)
!
!     Omega B1
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the B1 term, 
!
!        - sum_ckl u_kl^ac * g_kilc,
! 
!     and adds it to the singles projection vector (omeg1) of
!     the wavefunction object wfn
!
      implicit none
!
      class(ccsd) :: wf 
!
      integer(i15) :: a = 0, c = 0, k = 0, l = 0, ckl = 0, ki = 0
      integer(i15) :: ak = 0, akcl = 0, al = 0, alck = 0, ck = 0, ai = 0
      integer(i15) :: cl = 0, lc = 0, i = 0, j = 0
!
      real(dp), dimension(:,:), allocatable :: L_ki_J  ! L_ki^J 
      real(dp), dimension(:,:), allocatable :: L_lc_J  ! L_lc^J 
      real(dp), dimension(:,:), allocatable :: g_ki_lc ! g_kilc 
      real(dp), dimension(:,:), allocatable :: g_ckl_i ! g_kilc 
      real(dp), dimension(:,:), allocatable :: u_a_ckl ! u_kl^ac = 2 t_kl^ac - t_lk^ac

   end subroutine omega_b1_ccsd
!
!
   subroutine omega_c1_ccsd(wf)        
!
!        
!     C1 Omega 
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the C1 term, 
!
!        sum_ck F_kc*u_ai_ck,
!
!     where
!              
!        u_ai_ck = 2*t_ck_ai-t_ci_ak,
!    
!     and adds it to the projection vector (omega1) of the
!     wavefunction object wf.    
!
      implicit none
!
      class(ccsd) :: wf 
!
      real(dp), dimension(:,:), allocatable :: F_ck     
      real(dp), dimension(:,:), allocatable :: u_ai_ck  
      real(dp), dimension(:,:), allocatable :: omega1_ai
!
      integer(i15) :: i = 0, k = 0, c = 0, a = 0
      integer(i15) :: ck = 0, ai = 0, ak = 0, ci = 0, aick = 0, akci = 0
!
!     Allocation of F_ck, u_ai_ck and omega1_ai
!
      call allocator(F_ck, (wf%n_o)*(wf%n_v), 1)
      call allocator(u_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))  
      call allocator(omega1_ai, (wf%n_o)*(wf%n_v), 1)
!
   end subroutine omega_c1_ccsd
!
!
   subroutine omega_d1_ccsd(wf)
!
!     D1 omega term: Omega_ai^D1=F_ai_T1
!
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mars 2017
!
      implicit none
!
      class(ccsd) :: wf
!
!     Add F_a_i to omega
!
      call daxpy((wf%n_o)*(wf%n_v), one, wf%fock_ai, 1, wf%omega1, 1)
!
   end subroutine omega_d1_ccsd
!
!
end submodule