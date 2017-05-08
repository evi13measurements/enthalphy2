submodule (ccsd_class) omega
!
!
!                       -::- Omega submodule (CCSD) -::-
!           Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!
!     Contains the following family of procedures of the CCSD class:
!
!        initialize_omega: allocates the projection vector (omega1, omega2)
!                          and sets it to zero.
!
!        construct_omega:  constructs the projection vector (omega1, omega2) 
!                          for the current amplitudes (t1am, t2am) for the
!                          wavefunction object wf. The routine assumes that
!                          the projection vector is allocated.
!
!        omega_a1:         adds A1 term to omega1
!        omega_b1:         adds B1 term to omega1
!        omega_c1:         adds C1 term to omega1
!        omega_d1:         adds D1 term to omega1
!
!        omega_a2:         adds A2 term to omega2
!        omega_b2:         adds B2 term to omega2
!        omega_c2:         adds C2 term to omega2
!        omega_d2:         adds D2 term to omega2
!        omega_e2:         adds E2 term to omega2
!
   implicit none 
!
   logical :: debug = .false.
!
   real(dp) :: begin_timer
   real(dp) :: end_timer 
!
!
contains
!
!
   subroutine initialize_omega_ccsd(wf)
!
!     Initialize Omega (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Allocates the projection vector (omega) and sets it to zero.
!
      implicit none 
!
      class(ccsd) :: wf
!
      call allocator(wf%omega1, wf%n_v, wf%n_o)
      call dzero(wf%omega1, (wf%n_v)*(wf%n_o))
!
      call allocator(wf%omega2, wf%n_t2am, 1)
      call dzero(wf%omega2, wf%n_t2am)
!
   end subroutine initialize_omega_ccsd
!
!
   subroutine construct_omega_ccsd(wf)
!
!     Construct Omega (CCSD)
!     Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!
!     Directs the calculation of the omega vector.
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
      call dzero(wf%omega1, wf%n_t1am)
      call dzero(wf%omega2, wf%n_t2am)
!
!     Construct singles contributions 
!
      call wf%omega_a1
      call wf%omega_b1
      call wf%omega_c1
      call wf%omega_d1
!
!     Construct doubles contributions 
!
      call wf%omega_a2
      call wf%omega_b2
      call wf%omega_c2
      call wf%omega_d2
      call wf%omega_e2
!
       call cpu_time(omega_end)
       write(unit_output,*)'Time in omega:', omega_end-omega_start  
       call flshfo(unit_output)  
!
!
   end subroutine construct_omega_ccsd
!
!
   subroutine omega_a1_ccsd(wf)
!
!     Omega A1 term
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the A1 term, 
!
!        sum_ckd g_adkc * u_ki^cd,
!
!     and adds it to the singles projection vector (omeg1) of
!     the wavefunction object wf.
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
      integer(i15) :: i = 0, j = 0, a = 0, c = 0, k = 0, d = 0
      integer(i15) :: ad = 0, ad_dim = 0, ci = 0, cidk = 0, ck = 0 
      integer(i15) :: ckd = 0, ckdi = 0, di = 0, dk = 0, kc = 0, da = 0
!
      real(dp), dimension(:,:), allocatable :: L_kc_J 
      real(dp), dimension(:,:), allocatable :: L_da_J  ! L_ad^J; a is being batched over
      real(dp), dimension(:,:), allocatable :: g_da_kc ! g_adkc; a is being batched over
      real(dp), dimension(:,:), allocatable :: g_a_ckd ! reordered g_adkc; a is being batched over
      real(dp), dimension(:,:), allocatable :: u_ckd_i ! u_ki^cd
!
      logical :: reorder ! To get L_ab_J reordered, for batching over a (first index)
!
!     Allocate Cholesky vector L_kc_J
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call dzero(L_kc_J, (wf%n_o)*(wf%n_v)*(wf%n_J))
!
!     Get Cholesky vector L_kc_J
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate u_ckd_i = u_ki^cd
!
      call allocator(u_ckd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call dzero(u_ckd_i, ((wf%n_o)**2)*((wf%n_v)**2))
!
!     Calculate u_ckd_i
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do d = 1, wf%n_v
               do i = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  ckd  = index_three(c, k, d, wf%n_v, wf%n_o)
!
                  ck   = index_two(c, k, wf%n_v)
                  di   = index_two(d, i, wf%n_v)
                  ci   = index_two(c, i, wf%n_v)
                  dk   = index_two(d, k, wf%n_v)
!
                  ckdi = index_packed(ck, di)
                  cidk = index_packed(ci, dk)
!
!                 Calculate u_ckd_i
!
                  u_ckd_i(ckd, i) = two*(wf%t2am(ckdi, 1)) - wf%t2am(cidk, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Calculate the batching parameters over a = 1,2,...,n_vir,
!     for which we need to have enough room to store L_ad_J and g_ad_kc, and, 
!     later on in the same loop, g_ad_kc and g_a_ckd simultaneously
!
      available = get_available()
      required  = max(((wf%n_v)**2)*(wf%n_J) + (wf%n_v**2)*(wf%n_o)*(wf%n_v), &
                              2*((wf%n_v)**3)*(wf%n_o)) ! Eirik: redo this estimate !
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_begin, a_end, a_batch, max_batch_length, batch_dimension)
         batch_length = a_end - a_begin + 1 
!
!        Allocate the Cholesky vector L_da_J = L_ad^J
!
         ad_dim = batch_length*(wf%n_v) ! Dimension of ad for the batch over index a 
!
         call allocator(L_da_J, ad_dim, wf%n_J)
         call dzero(L_da_J, ad_dim*(wf%n_J))
!
!        Get reordered Cholesky vector L_da_J = L_ad^J 
!
         reorder = .true.
         call wf%get_cholesky_ab(L_da_J, a_begin, a_end, &
                                 ad_dim, reorder)
!
!        Allocate g_da_kc = g_adkc
!
         call allocator(g_da_kc, ad_dim, (wf%n_o)*(wf%n_v))
!
!        Calculate g_da_kc = sum_J L_da_J L_kc_J^T = sum_J L_ad^J L_kc^J = g_adkc 
!     
         call dgemm('N','T',            &
                     ad_dim,            &
                     (wf%n_o)*(wf%n_v), &
                     wf%n_J,            &
                     one,               &
                     L_da_J,            &
                     ad_dim,            &
                     L_kc_J,            &
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     g_da_kc,           &
                     ad_dim)
!
!        Deallocate the reordered Cholesky vector L_da_J
!
         call deallocator(L_da_J, ad_dim, wf%n_J)
!
!        Allocate g_a_ckd = g_adkc and set to zero
!
         call allocator(g_a_ckd, batch_length, (wf%n_o)*(wf%n_v)**2)
         call dzero(g_a_ckd, batch_length*(wf%n_o)*(wf%n_v)**2)
!
!        Reorder the integrals to g_a_ckd
!
         do c = 1, wf%n_v
            do k = 1, wf%n_o
               do d = 1, wf%n_v
                  do a = 1, batch_length
!
!                    Get the needed indices 
!
                     kc  = index_two(k, c, wf%n_o)
                     da  = index_two(d, a, wf%n_v)
                     ckd = index_three(c, k, d, wf%n_v, wf%n_o) 
!
!                    Set the value of reordered integral, g_a_ckd
!
                     g_a_ckd(a, ckd) = g_da_kc(da, kc) ! g_adkc 
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate reordered integrals g_da_kc
!
         call deallocator(g_da_kc, ad_dim, (wf%n_o)*(wf%n_v))
!
!        Calculate the A1 term (sum_ckd g_a,ckd * u_ckd,i) & add to the omega vector
! 
         call dgemm('N','N',               &
                     batch_length,         &
                     wf%n_o,               &
                     (wf%n_o)*(wf%n_v)**2, &
                     one,                  &
                     g_a_ckd,              &
                     batch_length,         &
                     u_ckd_i,              &
                     (wf%n_o)*(wf%n_v)**2, &
                     one,                  &
                     wf%omega1(a_begin,1), &
                     wf%n_v)
!
      enddo ! End of batches of the index a 
!
!     Print the omega vector 
!
      if (debug) then 
! 
         write(unit_output,*) 
         write(unit_output,*) 'Omega(a,i) after A1 term has been added:'
         write(unit_output,*)
!
         call vec_print(wf%omega1, wf%n_v, wf%n_o)
!
      endif
!
!     Deallocate vectors 
!
      call deallocator(u_ckd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call deallocator(g_a_ckd, batch_length, (wf%n_o)*(wf%n_v)**2)
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
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
!
!     Allocate Cholesky vectors L_ki,J and L_lc,J 
!
      call allocator(L_ki_J, (wf%n_o)**2, wf%n_J)
      call allocator(L_lc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call dzero(L_ki_J, (wf%n_J)*(wf%n_o)**2)
      call dzero(L_lc_J, (wf%n_o)*(wf%n_v)*(wf%n_J))
!
!     Read the Cholesky vectors L_ki_J and L_lc_J
!
      call wf%get_cholesky_ij(L_ki_J)
      call wf%get_cholesky_ia(L_lc_J)
!
!     Allocate integrals g_ki_lc = g_kilc
!
      call allocator(g_ki_lc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Calculate g_ki_lc = sum_J L_ki_J L_lc_J^T 
! 
      call dgemm('N','T',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ki_J,            &
                  (wf%n_o)**2,       &
                  L_lc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ki_lc,           &
                  (wf%n_o)**2)
!
!     Deallocate the Cholesky vectors L_ki_J and L_lc_J
!
      call deallocator(L_ki_J, (wf%n_o)**2, wf%n_J)
      call deallocator(L_lc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate reordered integrals g_ckl_i = g_kilc 
!
      call allocator(g_ckl_i, (wf%n_v)*((wf%n_o)**2), wf%n_o)
      call dzero(g_ckl_i, (wf%n_v)*((wf%n_o)**3))
!
!     Determine g_ckl_i = g_kilc 
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do l = 1, wf%n_o
               do i = 1, wf%n_o
!
!                 Calculate necessary indices
!
                  ckl = index_three(c, k, l, wf%n_v, wf%n_o) 
!
                  ki  = index_two(k, i, wf%n_o)                 
                  lc  = index_two(l, c, wf%n_o)           
!
!                 Set value of g_ckl_i
!
                  g_ckl_i(ckl, i) = g_ki_lc(ki, lc)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate unordered integrals g_ki_lc
!
      call deallocator(g_ki_lc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Allocate redordered u_a_ckl = u_kl^ac and set it to zero
!
      call allocator(u_a_ckl, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      call dzero(u_a_ckl, ((wf%n_v)**2)*(wf%n_o)**2)
!
!     Determine u_a_ckl = u_kl^ac
!
      do a = 1, wf%n_v
         do c = 1, wf%n_v
            do k = 1, wf%n_o
               do l = 1, wf%n_o
!
!                 Calculate necessary indices
!
                  ckl  = index_three(c, k, l, wf%n_v, wf%n_o) 
!
                  ak   = index_two(a, k, wf%n_v)
                  cl   = index_two(c, l, wf%n_v)
!
                  akcl = index_packed(ak, cl)
!
                  al   = index_two(a, l, wf%n_v)
                  ck   = index_two(c, k, wf%n_v)
!
                  alck = index_packed(al, ck)
!
!                 Set the value of u_a_ckl = u_kl^ac = 2*t_kl^ac - t_lk^ac = 2*t_ak,cl - t_al,ck 
!
                  u_a_ckl(a, ckl) = two*(wf%t2am(akcl, 1)) - wf%t2am(alck, 1)
!                  
               enddo
            enddo
         enddo
      enddo
!
!     Calculate the B1 term, - sum_ckl u_a_ckl g_ckl_i
!
      call dgemm('N','N',                 &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  -one,                   &
                  u_a_ckl,                &
                  wf%n_v,                 &
                  g_ckl_i,                &
                  (wf%n_v)*((wf%n_o)**2), &
                  one,                    &
                  wf%omega1,              &
                  wf%n_v) 
!
!     Print the omega vector 
!
      if (debug) then  
!
         write(unit_output,*) 
         write(unit_output,*) 'Omega(a,i) after B1 term has been added:'
         write(unit_output,*)
!
         call vec_print(wf%omega1, wf%n_v, wf%n_o)
!
      endif
!
!     Deallocate remaining vectors 
!
      call deallocator(u_a_ckl, wf%n_v, (wf%n_v)*((wf%n_o)**2))
      call deallocator(g_ckl_i, (wf%n_v)*((wf%n_o)**2), wf%n_o)
!
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
      call dzero(F_ck, (wf%n_o)*(wf%n_v))
      call dzero(u_ai_ck, ((wf%n_o)**2)*(wf%n_v)**2)
!
!     Set up u_ck_ai and virtual-occupied Fock matrix
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
!           Set up compound index
!  
            ck = index_two(c, k, wf%n_v)
!
!           Reorder MO Fock matrix
!
            F_ck(ck, 1) = wf%fock_ia(k, c)
!
            do a = 1, wf%n_v
               do i = 1, wf%n_o
!
!                 Set up compound indices
!
                  ai = index_two(a, i, wf%n_v)
                  ci = index_two(c, i, wf%n_v)
                  ak = index_two(a, k, wf%n_v)
!
                  aick = index_packed(ck, ai)
                  akci = index_packed(ci, ak)
!                    
!                 u_ck_ai
!
                  u_ai_ck(ai, ck) = two*(wf%t2am(aick, 1)) - wf%t2am(akci, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Matrix multiplication
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  F_ck,              &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  omega1_ai,         &
                  (wf%n_o)*(wf%n_v))
!
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
!           Set up compound index
!
            ai = index_two(a, i, wf%n_v)
!
!           Add to omega1
!
            wf%omega1(a,i) = wf%omega1(a, i) + omega1_ai(ai,1)
!
         enddo
      enddo
!
!     Deallocation
!
      call deallocator(F_ck, (wf%n_o)*(wf%n_v), 1)
      call deallocator(omega1_ai, (wf%n_o)*(wf%n_v), 1)
      call deallocator(u_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
   subroutine omega_a2_ccsd(wf)
!
!     MLCC Omega A2 term: Omega A2 = g_ai_bj + sum_(cd)g_ac_bd * t_ci_dj = A2.1 + A.2.2
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 10 Mar 2017
!
!     Structure: Batching over both a and b. If no batching is necessary L_ab_J is only read once, and g_ac_bd 
!                is constructed and kept in memory full size. 
!                g_ac_bd is reordered as g_ab_cd and t_ci_dj is reordered as t_cd_ij.
!                Omega contribution for A2.2 is ordered as Omega_ab_ij, and is reordered into the packed omega2 vector.          
!
      use utils
!
      implicit none
!
      class(ccsd) :: wf
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj 
      real(dp), dimension(:,:), allocatable :: g_ca_db 
      real(dp), dimension(:,:), allocatable :: g_p_ab_cd
      real(dp), dimension(:,:), allocatable :: g_m_ab_cd
      real(dp), dimension(:,:), allocatable :: L_ai_J 
      real(dp), dimension(:,:), allocatable :: L_ca_J 
      real(dp), dimension(:,:), allocatable :: L_db_J 
!
!     Reordered T2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_p_cd_ij
      real(dp), dimension(:,:), allocatable :: t_m_cd_ij
!
!     Reordered omega 2 
!  
      real(dp), dimension(:,:), allocatable :: omega2_p_ab_ij
      real(dp), dimension(:,:), allocatable :: omega2_m_ab_ij
! 
!     Indices
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: ab = 0, ca = 0, cb = 0, cd = 0, da = 0, db = 0 
      integer(i15) :: ai = 0, aj = 0, bj = 0, bi = 0, ci = 0, cj = 0, dj = 0, di = 0
      integer(i15) :: ij = 0
!
      integer(i15) :: aibj = 0, biaj = 0, cidj = 0, cjdi = 0
!
!     Batching and memory handling variables
!
      integer(i15) :: a_n_batch = 0, a_first = 0, a_last = 0, a_length = 0, a_max_length = 0, a_batch = 0
      integer(i15) :: b_n_batch = 0, b_first = 0, b_last = 0, b_length = 0, b_max_length = 0, b_batch = 0

      integer(i15) :: required = 0, available = 0
!
!     Logical for reordering in L_ab_J when batching over the last index
!
      logical :: reorder
!
!
!!!   A2.1 term   !!!
!
!
!     Create g_ai_bj
!  
      call allocator(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call allocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ai(L_ai_J)
!
!     g_ai_bj = sum_J L_ai_J*L_bj_J
!     
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               & 
                  L_ai_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ai_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ai_bj,           &
                  (wf%n_o)*(wf%n_v))
!
!
      call deallocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Add A2.1 to Omega 2
!     
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Set up compound indices
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  if(ai .ge. bj) then
!
                     aibj = index_packed(ai, bj)
!
!                    Add to omega
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + g_ai_bj(ai, bj)
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!!!   A2.2 term !!!
!
!
      required = 2*(wf%n_v)**2*(wf%n_J) &             ! Needed to get L_ca_J
               + 2*(wf%n_o)*(wf%n_v)*(wf%n_J)&        ! Needed to get L_db_J
               + (wf%n_v)**4 &                        ! Needed to get g_ac_bd
               + 2*(wf%n_v)**2*(packed_size(wf%n_v))  ! Needed to get g+-
      required = required*4                           ! Words

      available=get_available()
!
      a_max_length = 0
      call num_two_batch(required, available, a_max_length, a_n_batch, wf%n_v)
!
!     Initialize some variables for batching
!
      a_first  = 0
      a_last   = 0
      a_length = 0
!
!     Start looping over a-batches
!
      do a_batch = 1,a_n_batch
!   
         call batch_limits(a_first ,a_last ,a_batch, a_max_length, wf%n_v)
         a_length = a_last - a_first + 1     
!
!        Get cholesky vectors L_ac^J ordered as L_ca_J
!
         call allocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
         call dzero(L_ca_J, (wf%n_v)*a_length*(wf%n_J))
!
         reorder = .true.
         call wf%get_cholesky_ab(L_ca_J, a_first, a_last, (wf%n_v)*a_length, reorder)
!
!        Start looping over batches of b
!
         b_first  = 0
         b_last   = 0
         b_length = 0
!
         b_max_length = a_max_length
!
         do b_batch = 1, a_batch
!
            call batch_limits(b_first ,b_last ,b_batch, b_max_length, wf%n_v)
            b_length = b_last - b_first + 1 
!
!           Get cholesky vectors L_bd^J ordered as L_db_J
!
            call allocator(L_db_J, (wf%n_v)*b_length, wf%n_J)
            call dzero(L_db_J, (wf%n_v)*b_length*(wf%n_J))
!  
            reorder = .true.
            call wf%get_cholesky_ab(L_db_J, b_first, b_last, (wf%n_v)*b_length, reorder)
!
!           Allocate g_ca_db
!
            call allocator(g_ca_db, (wf%n_v)*a_length, (wf%n_v)*b_length)
!
!           g_ca_db = sum_J L_ca_J*L_db_J
!     
            call cpu_time(begin_timer)
            call dgemm('N','T',            &
                        (wf%n_v)*a_length, &
                        (wf%n_v)*b_length, &
                        wf%n_J,            &
                        one,               &
                        L_ca_J,            &
                        (wf%n_v)*a_length, &
                        L_db_J,            &
                        (wf%n_v)*b_length, &
                        zero,              &
                        g_ca_db,           &
                        (wf%n_v)*a_length)
            call cpu_time(end_timer)
            write(unit_output,*) 'Time to make g_acbd',end_timer-begin_timer
!
!           Deallocate L_db_J 
!
            call deallocator(L_db_J, (wf%n_v)*b_length, wf%n_J) 
!
!           Allocate for +-g, +-t
!
            call allocator(g_p_ab_cd, packed_size(a_length), packed_size(wf%n_v))
            call allocator(g_m_ab_cd, packed_size(a_length), packed_size(wf%n_v))
            call allocator(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
            call allocator(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
            call dzero(g_p_ab_cd, packed_size(a_length)*packed_size(wf%n_v))
            call dzero(g_m_ab_cd, packed_size(a_length)*packed_size(wf%n_v))
            call dzero(t_p_cd_ij, packed_size(wf%n_v)*packed_size(wf%n_o))
            call dzero(t_m_cd_ij, packed_size(wf%n_v)*packed_size(wf%n_o))
!
!           E: attempting to optimize (formaldehyde/pvtz, before: 0.73370899999999750 s)
!                                                         after:  0.53656699999999091 s)
!                 -114.33117412 new code 
!                 -114.3311740730 good...
!                 -114.33117412
!
            call cpu_time(begin_timer)
            do c = 1, wf%n_v 
               do d = 1, c
                  do b = 1, b_length
                     do a = 1, a_length
                        if ((a+a_first) .ge. (b+b_first)) then
!
!                          Calculate compound indices
!
                           ca = index_two(c, a, wf%n_v)
                           db = index_two(d, b, wf%n_v)
                           da = index_two(d, a, wf%n_v)
                           cb = index_two(c, b, wf%n_v)
!
                           ab = index_packed(a, b)
                           cd = index_packed(c, d)
!
!                          Reorder g_ca_db to g_ab_cd
!  
                           g_p_ab_cd(ab, cd) = g_ca_db(ca, db) + g_ca_db(da, cb)
                           g_m_ab_cd(ab, cd) = g_ca_db(ca, db) - g_ca_db(da, cb)
!
                          if(c .ne. d) then
                            g_p_ab_cd(ab, cd) = two*g_p_ab_cd(ab, cd)
                            g_m_ab_cd(ab, cd) = two*g_m_ab_cd(ab, cd)
                          endif
!                          
                        endif
                     enddo
                  enddo
!
                  do i = 1, wf%n_o
                     do j = 1, i
!
!                       Calculate compound indices
!     
                        ij = index_packed(i, j)
!     
                        ci = index_two(c, i, wf%n_v)
                        dj = index_two(d, j, wf%n_v)
                        cj = index_two(c, j, wf%n_v)
                        di = index_two(d, i, wf%n_v)
!  
                        cidj = index_packed(ci, dj)
                        cjdi = index_packed(cj, di)
!
!                       Reorder t_ci_dj to t_cd_ij
!  
                        t_p_cd_ij(cd, ij) = wf%t2am(cidj, 1) + wf%t2am(cjdi, 1)
                        t_m_cd_ij(cd, ij) = wf%t2am(cidj, 1) - wf%t2am(cjdi, 1)  
!
                     enddo
                  enddo
               enddo
            enddo
            call cpu_time(end_timer)
            write(unit_output,*) 'Time to reorder in A2:', end_timer-begin_timer
!
!           Dellocate g_ac_bd 
!
            call deallocator(g_ca_db, (wf%n_v)*a_length, (wf%n_v)*b_length)
!
!           Allocate omega +-
!
            call allocator(omega2_p_ab_ij, packed_size(a_length), packed_size(wf%n_o))
            call allocator(omega2_m_ab_ij, packed_size(a_length), packed_size(wf%n_o))
         !   omega2_p_ab_ij = zero (not needed)
          !  omega2_m_ab_ij = zero (not needed)
!  
!            omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
!  
            call cpu_time(begin_timer)
            call dgemm('N','N',                & 
                        packed_size(a_length), &
                        packed_size(wf%n_o),   &
                        packed_size(wf%n_v),   &
                        one/four,              &
                        g_p_ab_cd,             &
                        packed_size(a_length), &
                        t_p_cd_ij,             &
                        packed_size(wf%n_v),   &
                        zero,                  &
                        omega2_p_ab_ij,        &
                        packed_size(a_length))
!
            call dgemm('N','N',                & 
                        packed_size(a_length), &
                        packed_size(wf%n_o),   &
                        packed_size(wf%n_v),   &
                        one/four,              &
                        g_m_ab_cd,             &
                        packed_size(a_length), &
                        t_m_cd_ij,             &
                        packed_size(wf%n_v),   &
                        zero,                  &
                        omega2_m_ab_ij,        &
                        packed_size(a_length) )
            call cpu_time(end_timer)
            write(unit_output,*) 'Time for dgemm A2',end_timer-begin_timer
!
!           Deallocate +-g, +-t
!  
            call deallocator(g_p_ab_cd, packed_size(a_length), packed_size(wf%n_v))
            call deallocator(g_m_ab_cd, packed_size(a_length), packed_size(wf%n_v))
            call deallocator(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
            call deallocator(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
            do i = 1, wf%n_o
               do j = 1, i
                  do a = 1, a_length
                     do b = 1, a
!  
!                       Calculate compound indices
!     
                        Ai = index_two(a + a_first - 1, i, wf%n_v) ! A is full-space a index
                        Aj = index_two(a + a_first - 1, j, wf%n_v) ! A is full-space a index
                        Bj = index_two(b + b_first - 1, j, wf%n_v) ! B is full-space b index
                        Bi = index_two(b + b_first - 1, i, wf%n_v) ! B is full-space b index
!
!
                        ab = index_packed(a, b)
                        ij = index_packed(i, j)
!     
                        AiBj = index_packed(Ai, Bj)
                        BiAj = index_packed(Bi, Aj)
!                       
!                       Reorder into omega2_aibj
!  
                        wf%omega2(AiBj,1) = wf%omega2(AiBj, 1) + omega2_p_ab_ij(ab, ij) + omega2_m_ab_ij(ab, ij)
!
                        if (AiBj .ne. BiAj) then
                           wf%omega2(BiAj,1) = wf%omega2(BiAj, 1) + omega2_p_ab_ij(ab, ij) - omega2_m_ab_ij(ab, ij)
                        endif   
!     
                     enddo
                  enddo
               enddo
            enddo
!
!           Deallocate omega +-
!
            call deallocator(omega2_p_ab_ij, packed_size(a_length), packed_size(wf%n_o))
            call deallocator(omega2_m_ab_ij, packed_size(a_length), packed_size(wf%n_o))
!
         enddo
!
!        Deallocate L_ca_J
!
         call deallocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
      enddo
!      
!         write(unit_output,*) 
!         write(unit_output,*) 'Omega(aibj,1) after A2 term has been added:'
!         write(unit_output,*)
!         call vec_print(wf%omega2, wf%n_t2am, 1)
!
   end subroutine omega_a2_ccsd
!
!
!
   subroutine omega_b2_ccsd(wf)
!
!     MLCC Omega B2 term
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 11 Mar 2017
!
!     Omega B2 = sum_(kl) t_ak_bl*(g_kilj+sum_(cd) t_ci_dj * g_kc_ld) 
!
!     Structure: g_kilj is constructed first and reordered as g_kl_ij. 
!                Then the contraction over cd is performed, and the results added to g_kl_ij.
!                t_ak_bl is then reordered as t_ab_kl and the contraction over kl is performed.
!
      implicit none
!
      class(ccsd) :: wf 
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_kc_J     
      real(dp), dimension(:,:), allocatable :: L_ij_J  
      real(dp), dimension(:,:), allocatable :: g_kc_ld    
      real(dp), dimension(:,:), allocatable :: g_kl_cd    
      real(dp), dimension(:,:), allocatable :: g_kl_ij    
      real(dp), dimension(:,:), allocatable :: g_ki_lj 
!
!     Reordered T2 apmlitudes
!   
      real(dp), dimension(:,:), allocatable :: t_cd_ij    
      real(dp), dimension(:,:), allocatable :: t_ab_kl   
!
!     Intermediate for matrix multiplication
! 
      real(dp), dimension(:,:), allocatable :: X_kl_ij 
!
!     Reordered omega
!   
      real(dp), dimension(:,:), allocatable :: omega_ab_ij
!
!     Indices
!   
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ab = 0, cd = 0
      integer(i15) :: ai = 0, ak = 0, bj = 0, bl = 0, ci = 0, dj = 0
      integer(i15) :: kc = 0, ld = 0
      integer(i15) :: ij = 0, ki = 0, kl = 0, lj = 0
!
      integer(i15) :: aibj = 0, akbl = 0, cidj = 0 
!
!     Read cholesky vector of type L_ij_J
!
      call allocator(L_ij_J, (wf%n_o)*(wf%n_o), wf%n_J)
      L_ij_J = zero 
!
      call wf%get_cholesky_ij(L_ij_J)
!
!     Create g_ki_lj = sum_J L_li_J*L_lj_J
!
      call allocator(g_ki_lj, (wf%n_o)*(wf%n_o), (wf%n_o)*(wf%n_o)) 
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_o), &
                  (wf%n_o)*(wf%n_o), &
                  wf%n_J,            &
                  one,               &
                  L_ij_J,            &
                  (wf%n_o)*(wf%n_o), &
                  L_ij_J,            &
                  (wf%n_o)*(wf%n_o), &
                  zero,              &
                  g_ki_lj,           &
                  (wf%n_o)*(wf%n_o))
!
!
      call deallocator(L_ij_J, (wf%n_o)*(wf%n_o), wf%n_J)
!
      call allocator(g_kl_ij, (wf%n_o)*(wf%n_o),(wf%n_o)*(wf%n_o))
      g_kl_ij = zero
!
      do k = 1, wf%n_o
         do l = 1, wf%n_o
            do i = 1, wf%n_o
               do j=1, wf%n_o
!
!                 Calculate compound indices
!
                  ki = index_two(k, i, wf%n_o)
                  lj = index_two(l, j, wf%n_o)
                  kl = index_two(k, l, wf%n_o)
                  ij = index_two(i, j, wf%n_o)
!
!                  Reordering g_ki_lj to g_kl_ij
!
                  g_kl_ij(kl, ij) = g_ki_lj(ki, lj)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ki_lj, (wf%n_o)*(wf%n_o), (wf%n_o)*(wf%n_o))
!
!     Read cholesky vectors of ia-type into L_kc_J
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Create g_ck_ld = sum_(J) L_kc_J*L_ld_J
!
      call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!  
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate cholesky vectors L_ck_J
!
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Reorder g_kc_ld as g_kl_cd, also reordering t_ci_dj as t_cd_ij
!
      call allocator(t_cd_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call allocator(g_kl_cd, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
            do l = 1, wf%n_o
               do d = 1, wf%n_v
!
!                 Calculate compound indices
!
                  cd = index_two(c, d, wf%n_v)
                  kl = index_two(k, l, wf%n_o)
                  kc = index_two(k, c, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
! 
                  ci = index_two(c, k, wf%n_v)
                  dj = index_two(d, l, wf%n_v)
                  ij = index_two(k, l, wf%n_o)
!
                  cidj = index_packed(ci, dj)
!
!                 Reordering
!
                  g_kl_cd(kl, cd) = g_kc_ld(kc, ld)
                  t_cd_ij(cd, ij) = wf%t2am(cidj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_kc_ld
!
      call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_kl_cd,     &
                  (wf%n_o)**2, &
                  t_cd_ij,     &
                  (wf%n_v)**2, &
                  one,         &
                  g_kl_ij,     &
                  (wf%n_o)**2)
!
!     Deallocate t_cd_ij and g_kl_cd
!
      call deallocator(t_cd_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call deallocator(g_kl_cd, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
!     Reorder t_ak_bl to t_ab_kl
!
      call allocator(t_ab_kl, (wf%n_v)**2, (wf%n_o)**2)
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do k = 1, wf%n_o
               do l=1, wf%n_o
!
!                 Calculate compound indices
!
                  ak = index_two(a, k, wf%n_v)
                  bl = index_two(b, l, wf%n_v)
                  ab = index_two(a, b, wf%n_v)
                  kl = index_two(k, l, wf%n_o)
!
                  akbl = index_packed(ak, bl)
!
                  t_ab_kl(ab, kl) = wf%t2am(akbl, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     omega_ab_ij = sum_(kl) t_ab_kl*X_kl_ij
!
      call allocator(omega_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_ab_kl,     &
                  (wf%n_v)**2, &
                  g_kl_ij,     &
                  (wf%n_o)**2, &
                  zero,        &
                  omega_ab_ij, &
                  (wf%n_v)**2)
!
      call deallocator(t_ab_kl, (wf%n_v)**2, (wf%n_o)**2)
      call deallocator(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
!     Reorder omega
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do j = 1, wf%n_o
!
!                 Calculate compound indices
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  if (ai .ge. bj) then
!
                     ab = index_two(a, b, wf%n_v)
                     ij = index_two(i, j, wf%n_o)
!
                     aibj = index_packed(ai, bj)
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + omega_ab_ij(ab, ij)
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(omega_ab_ij,(wf%n_v)*(wf%n_v),(wf%n_o)*(wf%n_o))  
!
!     Print the omega vector, having added B2
!
      if (debug) then 
         write(unit_output,*) 
         write(unit_output,*) 'Omega(aibj,1) after B2 term has been added:'
         write(unit_output,*)
         call vec_print(wf%omega2, wf%n_t2am, 1)
      endif
!
   end subroutine omega_b2_ccsd
!
!
   subroutine omega_c2_ccsd(wf)
!
!     C2 omega term. Omega C2 = -1/2* sum_(ck)t_bk_cj*(g_ki_ac -1/2 sum_(dl)t_al_di * g_kd_lc)
!                               - sum_(ck) t_bk_ci (g_kj_ac-sum_(dl)t_al_dj*g_kd_lc)
!
!
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!     
      implicit none
!
      class(ccsd) :: wf
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_ia_J 
      real(dp), dimension(:,:), allocatable :: L_ki_J 
      real(dp), dimension(:,:), allocatable :: L_ca_J 
      real(dp), dimension(:,:), allocatable :: g_kd_lc
      real(dp), dimension(:,:), allocatable :: g_dl_ck
      real(dp), dimension(:,:), allocatable :: g_ki_ca
      real(dp), dimension(:,:), allocatable :: g_ai_ck
!
!     Reordered T2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_ai_dl
      real(dp), dimension(:,:), allocatable :: t_ck_bj
!
!     Intermediates for matrix multiplication
!
      real(dp), dimension(:,:), allocatable :: X_ai_ck
      real(dp), dimension(:,:), allocatable :: Y_ai_bj
!  
!     Indices
!     
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ca = 0
      integer(i15) :: ai = 0, aj = 0, al = 0, bi = 0, bj = 0, bk = 0, cj = 0, ck = 0, cl = 0, di = 0, dk = 0, dl = 0
      integer(i15) :: kd = 0, lc = 0
      integer(i15) :: ki = 0
!
      integer(i15) :: aldi = 0, aibj = 0, cldk = 0, bkcj = 0
!
!     Batching and memory handling
!
      integer(i15) :: required = 0, available = 0
!
      integer(i15) :: n_batch = 0, max_batch_length = 0
      integer(i15) :: a_batch = 0, a_start = 0, a_end = 0, a_length = 0 
!
!     Logical for reordering L_ab_J when batching over last index 
!
      logical :: reorder 
!
!     Allocate L_ia_J
!
      call allocator(L_ia_J,(wf%n_o)*(wf%n_v),(wf%n_J))
     ! L_ia_J=zero
!
!     Get L_ia_J
!
      call wf%get_cholesky_ia(L_ia_J)
!
!     Create g_kd_lc = sum_J L_kd_J * L_lc_J
!
      call allocator(g_kd_lc,(wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
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
                  g_kd_lc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate L_ia_J
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Reorder g_kd_lc as g_dl_ck and t_al_di as t_ai_dl
!
      call allocator(g_dl_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call allocator(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dzero(g_dl_ck, (wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v))
      call dzero(t_ai_dl, (wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v))
!
      do c = 1, wf%n_v
         do d = 1, wf%n_v
            do k = 1, wf%n_o
               do l = 1, wf%n_o
!
!                 Calculate compound indices
!
                  kd = index_two(k, d, wf%n_o)
                  lc = index_two(l, c, wf%n_o)
                  dl = index_two(d, l, wf%n_v)
                  ck = index_two(c, k, wf%n_v)
!
                  cl = index_two(c, l, wf%n_v)
                  dk = index_two(d, k, wf%n_v)
!
                  cldk = index_packed(cl, dk)
!
!                 Reorderings
!
                  g_dl_ck(dl, ck) = g_kd_lc(kd, lc)
!
                  t_ai_dl(ck, dl) = wf%t2am(cldk, 1)
!
               enddo
            enddo
         enddo
      enddo

      call deallocator(g_kd_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     -1/2*sum_(dl) t_ai_dl*g_dl_ck = X_ai_ck
!
      call allocator(X_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -half,             &
                  t_ai_dl,           &
                  (wf%n_o)*(wf%n_v), &
                  g_dl_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ai_ck,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate L_ia_J and g_dl_ck, 
!
      call deallocator(g_dl_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Constructing g_ki_ac ordered as g_ki_ca
!
!     Allocate g_ki_ca
!
      call allocator(g_ki_ca, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
      call dzero(g_ki_ca, (wf%n_o)*(wf%n_o)*(wf%n_v)*(wf%n_v))
!
!     Allocate L_ki_J
!
      call allocator(L_ki_J, (wf%n_o)*(wf%n_o), wf%n_J) 
!
!     Get cholesky vectors of ij-type
!
      call wf%get_cholesky_ij(L_ki_J)
!
!     Prepare batching over a 
!
!
!     Setup of variables needed for batching
!
      available = get_available()
      required = 2*(wf%n_v)*(wf%n_v)*(wf%n_J)*4 + 2*(wf%n_v)*(wf%n_o)*(wf%n_J)*4
      call num_batch(required, available, max_batch_length, n_batch, wf%n_v)
!
      a_start  = 1
      a_end    = 0
      a_length = 0
!
!     Start looping over batches
!
      do a_batch = 1,n_batch
!
!        Get batch limits  and  length of batch
!
         call batch_limits(a_start, a_end, a_batch, max_batch_length, wf%n_v)
         a_length = a_end - a_start + 1
!
!        Allocation for L_ac_J as L_ca_J (L_ca_J = L_acJ)
!
         call allocator(L_ca_J,(wf%n_v)*a_length, wf%n_J)
        ! L_ca_J=zero
!
!        Read Cholesky vectors
!
         reorder = .true.
         call wf%get_cholesky_ab(L_ca_J, a_start, a_end, (wf%n_v)*a_length, reorder)
!
!        g_ki_ca = sum_J L_ki_J*L_ca_J
!
         call dgemm('N','T',                                   &
                     (wf%n_o)*(wf%n_o),                        &
                     (wf%n_v)*a_length,                        &
                     wf%n_J,                                   &   
                     one,                                      &   
                     L_ki_J,                                   &
                     (wf%n_o)*(wf%n_o),                        &
                     L_ca_J,                                   &
                     (wf%n_v)*a_length,                        &
                     one,                                      &
                     g_ki_ca(1,index_two(1, a_start, wf%n_v)), &
                     (wf%n_o)*(wf%n_o))
!
!        Deallocate L_ca_J
!
         call deallocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
!
      enddo ! End of batching
!
!     Deallocate L_ki_J
!
      call deallocator(L_ki_J, (wf%n_o)*(wf%n_o), wf%n_J)
!
!     Reorder g_ki_ca to g_ai_ck
!
      call allocator(g_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      do i = 1, wf%n_o
         do k = 1, wf%n_o
            do a = 1, wf%n_v
               do c = 1, wf%n_v
!
!                 Calculate compound indices
!
                  ki = index_two(k, i, wf%n_o)
                  ca = index_two(c, a, wf%n_v)
                  ai = index_two(a, i, wf%n_v)
                  ck = index_two(c, k, wf%n_v)
!
                  g_ai_ck(ai, ck) = g_ki_ca(ki, ca)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ki_ca, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
!     X_ai_ck = X_ai_ck + g_ai_ck
!
      call daxpy((wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v), one, g_ai_ck, 1, X_ai_ck, 1)
!
!     Deallocate g_ai_kc
!
      call deallocator(g_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_bkcj_1 as t_ck_bj
!
      call allocator(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Needed indices
!
                  bk = index_two(b, k, wf%n_v)
                  cj = index_two(c, j, wf%n_v)
!
                  bkcj = index_packed(bk, cj)
!
                  bj = index_two(b, j, wf%n_v)
                  ck = index_two(c, k, wf%n_v)
!
                  t_ck_bj(ck, bj) = wf%t2am(bkcj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate intermediate Y_ai_bj
!
      call allocator(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
     ! Y_ai_bj = zero
!
!     Y_ai_bj = - sum_(ck) X_ai_ck*t_ck_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  X_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  t_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_ai_bj,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate the X intermediate
!
      call deallocator(X_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Deallocate t_ck_bj
!
      call deallocator(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Omega_aibj,1 = P_ai_bj ( 1/2*Y_ai_bj + Y_aj_bi )
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate compound indices
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  if (ai .ge. bj) then
                     aj = index_two(a, j, wf%n_v)
                     bi = index_two(b, i, wf%n_v)
!
                     aibj=index_packed(ai, bj)
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + half*Y_ai_bj(ai, bj) + Y_ai_bj(aj, bi) &
                                                               + half*Y_ai_bj(bj, ai) + Y_ai_bj(bi, aj)
!
                  endif
!  
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate intermediate Y_ai_bj
!
      call deallocator(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Print the omega vector, having added C2
!
      if (debug) then 
         write(unit_output,*) 
         write(unit_output,*) 'Omega(aibj,1) after C2 term has been added:'
         write(unit_output,*)
         call vec_print(wf%omega2, wf%n_t2am, 1)
      endif 
!
   end subroutine omega_c2_ccsd
!
!
   subroutine omega_d2_ccsd(wf)
!
!     Omega D2 
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the D2 term,
!
!        sum_ck u_jk^bc g_aikc 
!        - 1/2 * sum_ck u_jk^bc g_acki 
!        + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!
!     where 
!
!        u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!        L_ldkc  = 2 * g_ldkc  - g_lckd.
!
!     The first, second, and third terms are referred to as D2.1, D2.2, and D2.3, 
!     and comes out ordered as (ai,bj). All terms are added to the omega vector of the 
!     wavefunction object wf.
!
!     The routine adds the terms in the following order: D2.3, D2.1, D2.2
!
      implicit none 
!
      class(ccsd) :: wf 
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_batch_length = 0, batch_dimension = 0, n_batch = 0
      integer(i15) :: a_begin = 0, a_end = 0, a_batch = 0, batch_length = 0, a_full = 0, ac_dim = 0 
!
!     Indices 
!
      integer(i15) :: ai = 0, aidl = 0, al = 0, aldi = 0, a = 0, i = 0, b = 0, ca = 0
      integer(i15) :: j = 0, c = 0, d = 0, di = 0, dl = 0, k = 0, kc = 0, kd = 0, l = 0, ki = 0
      integer(i15) :: lc = 0, ld = 0, aibj = 0, bj = 0, bjck = 0, bk = 0, bkcj = 0, cj = 0, ck = 0
!
      real(dp), dimension(:,:), allocatable :: omega2_ai_bj ! For storing D2.3, D2.2 & D2.1
!
!     Vectors for D2.3 term 
!
      real(dp), dimension(:,:), allocatable :: L_kc_J  ! L_kc^J 
      real(dp), dimension(:,:), allocatable :: g_ld_kc ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: L_ld_kc ! L_ldkc = 2 * g_ldkc - g_lckd 
      real(dp), dimension(:,:), allocatable :: u_ai_ld ! u_il^ad = 2 * t_il^ad - t_li^ad 
      real(dp), dimension(:,:), allocatable :: Z_ai_kc ! An intermediate, see below
!
!     Vectors for D2.2 term 
!
      real(dp), dimension(:,:), allocatable :: g_ai_kc ! g_aikc 
      real(dp), dimension(:,:), allocatable :: u_kc_bj ! u_jk^bc
      real(dp), dimension(:,:), allocatable :: L_ai_J  ! L_ai^J 
!
!     Vectors for D2.1 term 
!
      real(dp), dimension(:,:), allocatable :: g_ai_ck ! g_acki
      real(dp), dimension(:,:), allocatable :: g_ca_ki ! g_acki; a is batched over 
      real(dp), dimension(:,:), allocatable :: L_ca_J  ! L_ac^J; a is batched over 
      real(dp), dimension(:,:), allocatable :: L_ki_J  ! L_ki^J 
      real(dp), dimension(:,:), allocatable :: u_ck_bj ! u_jk^bc
!
!     Logical for reordering L_ab_J when batching over the last index 
!
      logical :: reorder 
!
!     Allocate the Cholesky vector L_kc_J = L_kc^J and set to zero 
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
   !   L_kc_J = zero
!
!     Read the Cholesky vector L_kc_J from file
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
    !  g_ld_kc = zero 
!
!     Calculate g_ld_kc = g_ldkc = sum_J L_ld^J L_kc^J
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ld_kc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Allocate L_ld_kc = L_ldkc and set to zero    
!
      call allocator(L_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call dzero(L_ld_kc, (wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v))
!
!     Determine L_ld_kc = L_ldkc from g_ld_kc = g_ldkc 
!
      do l = 1, wf%n_o
         do d = 1, wf%n_v
            do k = 1, wf%n_o
               do c = 1, wf%n_v
!
!                 Calculate the necessary indices 
!
                  ld = index_two(l, d, wf%n_o)
                  kc = index_two(k, c, wf%n_o)
!
                  lc = index_two(l, c, wf%n_o)
                  kd = index_two(k, d, wf%n_o)
!
!                 Set the value of L_ld_kc = 2 * g_ldkc - g_lckd 
!
                  L_ld_kc(ld, kc) = two*g_ld_kc(ld, kc) - g_ld_kc(lc, kd)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ld_kc and L_kc_J
!
      call deallocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate u_ai_ld = u_il^ad and set to zero 
!
      call allocator(u_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call dzero(u_ai_ld, (wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v))
! 
!     Determine u_ai_ld = u_il^ad = 2 * t_il^ad - t_li^ad 
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do l = 1, wf%n_o
               do d = 1, wf%n_v
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a, i, wf%n_v)
                  ld   = index_two(l, d, wf%n_o)
!
                  dl   = index_two(d, l, wf%n_v)
                  aidl = index_packed(ai, dl)
!
                  al   = index_two(a, l, wf%n_v)
                  di   = index_two(d, i, wf%n_v)
!
                  aldi = index_packed(al, di)
!
!                 Set the value of u_ai_ld 
!
                  u_ai_ld(ai, ld) = two*(wf%t2am(aidl, 1)) - wf%t2am(aldi, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the intermediate Z_ai_kc = sum_dl u_ai_ld L_ld_kc and set it to zero
!
      call allocator(Z_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
    !  Z_ai_kc = zero
!
!     Form the intermediate Z_ai_kc = sum_dl u_ai_ld L_ld_kc
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_ai_ld,           &
                  (wf%n_o)*(wf%n_v), &
                  L_ld_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Z_ai_kc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate L_ld_kc
!
      call deallocator(L_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Allocate the D2.3 term omega2_ai_bj and set it to zero
!
      call allocator(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
    !  omega2_ai_bj = zero
!
!     Form the D2.3 term, 1/4 sum_kc Z_ai_kc u_kc_bj = 1/4 sum_kc Z_ai_kc(ai,kc) u_ai_ld(bj,kc)
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one/four,          &
                  Z_ai_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  u_ai_ld,           &
                  (wf%n_o)*(wf%n_v), & 
                  zero,              &
                  omega2_ai_bj,      &
                  (wf%n_o)*(wf%n_v))
!
!     Some mathematical justification for the above matrix multiplication. We have 
!
!           1/4 * sum_ck (sum_dl u_il^ad L_ldkc) u_jk^bc = 1/4 * sum_ck Z_ai,kc u_kc,bj,
!
!     where Z_ai,kc = sum_dl u_ai,ld L_ld,kc. Note that u_ai_ld(ai,ld) = u_il^ad, 
!     which means that u_ai_ld(bj,kc)^T = u_ai_ld(kc,bj) = u_kj^cb = u_jk^bc.
!
!
!     Deallocate the Z_ai_kc intermediate 
!
      call deallocator(Z_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add the D2.3 term to the omega vector 
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a, i, wf%n_v)
                  bj   = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + omega2_ai_bj(ai, bj) & 
                                                               + omega2_ai_bj(bj, ai)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the omega2_ai_bj and u_ai_ld(ai,ld) = u_il^ad vector
!
      call deallocator(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(u_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Print the omega vector, having added the D2.3 term 
!
      if (debug) then 
!
         write(unit_output,*) 
         write(unit_output,*) 'Omega(aibj,1) after D2.3 term has been added:'
         write(unit_output,*)
!
         call vec_print(wf%omega2, wf%n_t2am, 1)
!
      endif                
!
!     Allocate the L_ai_J and L_kc_J terms and set them to zero 
!
      call allocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
   !   L_ai_J = zero
   !   L_kc_J = zero
!
!     Read the Cholesky vectors from file 
!
      call wf%get_cholesky_ai(L_ai_J)
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate g_ai_kc = g_aikc and set it zero
!
      call allocator(g_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
    !  g_ai_kc = zero 
!
!     Form the g_ai_kc integrals from L_ai_J and L_kc_J
!  
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ai_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ai_kc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate the Cholesky vectors L_ai_J and L_kc_J
!
      call deallocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate u_kc_bj and set it to zero 
!
      call allocator(u_kc_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call dzero(u_kc_bj, (wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v))
!
!     Determine u_kc_bj = u_jk^bc = 2 * t_jk^bc - t_kj^bc 
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  kc   = index_two(k, c, wf%n_o)
                  bj   = index_two(b, j, wf%n_v)
                  ck   = index_two(c, k, wf%n_v)
!
                  bjck = index_packed(bj, ck)
!
                  bk   = index_two(b, k, wf%n_v)
                  cj   = index_two(c, j, wf%n_v)
!
                  bkcj = index_packed(bk, cj)
!
!                 Set the value of u_kc_bj
!     
                  u_kc_bj(kc, bj) = two*(wf%t2am(bjck, 1)) - wf%t2am(bkcj, 1)
!
               enddo
            enddo
         enddo
      enddo 
!
!     Allocate omega2_ai_bj and set it to zero 
!
      call allocator(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
    !  omega2_ai_bj = zero 
!
!     Calculate the D2.1 term sum_ck u_jk^bc g_aikc = sum_ck g_ai_kc u_kc_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_ai_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  u_kc_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  omega2_ai_bj,      &
                  (wf%n_o)*(wf%n_v))
!
!     Add the D2.1 term to the omega vector 
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a, i, wf%n_v)
                  bj   = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
!
                     wf%omega2(aibj,1) = wf%omega2(aibj,1) + omega2_ai_bj(ai,bj) & 
                                                            + omega2_ai_bj(bj,ai)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ai_kc, u_kc_bj, and the omega2_ai_bj vectors 
!
      call deallocator(g_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(u_kc_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call deallocator(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Print the omega vector, having added the D2.1 term 
!
      if (debug) then 
!
         write(unit_output,*) 
         write(unit_output,*) 'Omega(aibj,1) after D2.1 term has been added:'
         write(unit_output,*)
         call vec_print(wf%omega2, wf%n_t2am, 1)
!
      endif 
!
!     - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj 
!
!     Allocate L_ki_J and set it to zero
!
      call allocator(L_ki_J, (wf%n_o)**2, wf%n_J)
   !   L_ki_J = zero 
!
!     Read the Cholesky vector L_ki_J from file 
!
      call wf%get_cholesky_ij(L_ki_J)
!
!     Allocate the full g_ai_ck = g_acki and set it to zero 
!
      call allocator(g_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call dzero(g_ai_ck, (wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v))
!
!     Prepare for batching over the index a to calculate g_ai_ck = g_acki
!
!        To calculate this term, we need to hold L_ac^J and g_acki
!        in memory simultaneously 
!
      required = (wf%n_J)*(wf%n_v)**2 + ((wf%n_v)**2)*((wf%n_o)**2) ! Eirik: This estimate has to be updated,
                                                                      ! when we account for the T1 transformation
      available = get_available()
      batch_dimension = wf%n_v
!
!     Determine the batching variables 
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension) 
!
!     Determine g_ai_ck = g_acki successively in batches over a 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index        
!
         call batch_limits(a_begin, a_end, a_batch, max_batch_length, batch_dimension)
         batch_length = a_end - a_begin + 1 
!
!        Allocate the Cholesky vector L_ca_J = L_ac^J and set it to zero 
!
         ac_dim = batch_length*(wf%n_v) ! Dimension of ac for the batch over index a 
!
         call allocator(L_ca_J, ac_dim, wf%n_J)
      !   L_ca_J = zero
!
!        Read the Cholesky vector from file 
!
         reorder = .true.
         call wf%get_cholesky_ab(L_ca_J, a_begin, a_end, ac_dim, reorder)
!
!        Allocate the integral g_ca_ki = g_acki and set to zero 
!
         call allocator(g_ca_ki, ac_dim, (wf%n_o)**2)
         ! g_ca_ki = zero ! Not needed
!
!        Calculate g_ca_ki = g_acki from L_ca_J = L_ac^J and L_ki_J = L_ki^J
!
         call dgemm('N','T',      &
                     ac_dim,      &
                     (wf%n_o)**2, &
                     wf%n_J,      &
                     one,         &
                     L_ca_J,      &
                     ac_dim,      &
                     L_ki_J,      &
                     (wf%n_o)**2, &
                     zero,        & ! E: Changed from one to zero
                     g_ca_ki,     &
                     ac_dim)
!
!        Reorder the integrals g_ca_ki (reduced a) = g_acki = g_ai_ck (full a)
!
         do a = 1, batch_length
            do c = 1, wf%n_v
               do k = 1, wf%n_o
                  do i = 1, wf%n_o
!
!                    Calculate the necessary indices 
!
                     a_full = a - 1 + a_begin ! The full matrix index a
!
                     ai = index_two(a_full, i, wf%n_v)
                     ck = index_two(c, k, wf%n_v)
!
                     ca = index_two(c, a, wf%n_v)
                     ki = index_two(k, i, wf%n_o)
!
!                    Set the value of g_ai_ck = g_acki 
!
                     g_ai_ck(ai, ck) = g_ca_ki(ca, ki) ! g_acki
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate the g_ca_ki and L_ca_J vectors
!
         call deallocator(g_ca_ki, ac_dim, (wf%n_o)**2)
         call deallocator(L_ca_J, ac_dim, wf%n_J)
!
      enddo ! End of loop over batches of a 
!
!     Allocate the u_ck_bj = u_jk^bc vector and set it to zero 
!
      call allocator(u_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call dzero(u_ck_bj, (wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v))
!
!     Determine u_ck_bj = u_jk^bc = 2 * t_jk^bc - t_kj^bc 
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  ck = index_two(c, k, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  bjck = index_packed(bj, ck)
!
                  bk   = index_two(b, k, wf%n_v)
                  cj   = index_two(c, j, wf%n_v)
!
                  bkcj = index_packed(bk, cj)
!
!                 Set the value of u_ck_bj 
!
                  u_ck_bj(ck, bj) = two*(wf%t2am(bjck, 1)) - wf%t2am(bkcj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the D2.2 term and set it to zero 
!
      call allocator(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
 !     omega2_ai_bj = zero
!
!     Calculate the D2.2 term, - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one/two,          &
                  g_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  u_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  omega2_ai_bj,      &
                  (wf%n_o)*(wf%n_v))
!
!     Add the D2.2 term to the omega vector 
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + omega2_ai_bj(ai,bj) &
                                                               + omega2_ai_bj(bj,ai)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo      
!
!     Print the omega vector, having added the D2.2 term 
!
      if (debug) then 
!
         write(unit_output,*) 
         write(unit_output,*) 'Omega(aibj,1) after D2.2 term has been added:'
         write(unit_output,*)
!
         call vec_print(wf%omega2, wf%n_t2am, 1)
!
      endif 
!
!     Deallocate g_ai_ck, u_ck_bj, and the temporary omega2_ai_bj
!
      call deallocator(g_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(u_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call deallocator(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call deallocator(L_ki_J,(wf%n_o)**2, wf%n_J)
!
   end subroutine omega_d2_ccsd
!
! 
   subroutine omega_e2_ccsd(wf)
!
!     Omega E2
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the E2 term,
!
!        sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd) 
!        - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!
!     where
!
!        u_kl^bc = 2 * t_kl^bc - t_lk^bc.
!
!     The first term is referred to as the E2.1 term, and comes out ordered as (b,jai).
!     The second term is referred to as the E2.2 term, and comes out ordered as (aib,j).
!
!     Both are permuted added to the projection vector element omega2(ai,bj) of
!     the wavefunction object wf.
!
      implicit none 
!
      class(ccsd) :: wf 
!
!     Indices 
!
      integer(i15) :: aib = 0, aibk = 0, bk = 0, bja = 0, ibj = 0, aibj = 0, dlck = 0
      integer(i15) :: b = 0, c = 0, k = 0, d = 0, ck = 0, ckdl = 0, cl = 0, cldk = 0
      integer(i15) :: dk = 0, dl = 0, kc = 0, kdl = 0, l = 0, ld = 0, a = 0, ai = 0
      integer(i15) :: bj = 0, aicj = 0, cj = 0, i = 0, j = 0, jai = 0, dlc = 0, dkcl = 0
!
!     Vectors for E2.1 term 
!
      real(dp), dimension(:,:), allocatable :: omega2_b_jai ! For storing the E2.1 term temporarily
      real(dp), dimension(:,:), allocatable :: L_kc_J       ! L_kc^J
      real(dp), dimension(:,:), allocatable :: g_ld_kc      ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: g_kdl_c      ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: u_b_kdl      ! u_kl^bd 
      real(dp), dimension(:,:), allocatable :: X_b_c        ! An intermediate, see below for definition
      real(dp), dimension(:,:), allocatable :: t_c_jai      ! t_ij^ac 
!
!     Vectors for E2.2 term 
!
      real(dp), dimension(:,:), allocatable :: g_k_dlc      ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: u_dlc_j      ! u_lj^dc 
      real(dp), dimension(:,:), allocatable :: omega2_aib_j ! For storing the E2.2 term temporarily
      real(dp), dimension(:,:), allocatable :: Y_k_j        ! An intermediate, see below for definition 
      real(dp), dimension(:,:), allocatable :: t_aib_k      ! t_ik^ab 
!
!     Allocate the Cholesky vector L_kc_J = L_kc^J and set to zero 
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
   !   L_kc_J = zero
!
!     Read the Cholesky vector from file 
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
    !  g_ld_kc = zero
!
!     Calculate g_ld_kc = sum_J L_ld^J L_kc^J 
!
      call dgemm('N','T',            & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ld_kc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate the Cholesky vector L_kc_J
!
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate u_b_kdl = u_kl^bd and set to zero
!
      call allocator(u_b_kdl, wf%n_v, (wf%n_v)*((wf%n_o)**2))
      call dzero(u_b_kdl, (wf%n_v)*(wf%n_v)*((wf%n_o)**2))
!
!     Allocate g_kdl_c = g_ldkc and set to zero 
!
      call allocator(g_kdl_c, (wf%n_v)*((wf%n_o)**2), wf%n_v)
      call dzero(g_kdl_c, (wf%n_v)*((wf%n_o)**2)*(wf%n_v))
!
!     Determine u_b_kdl = u_kl^bd and g_kdl_c = g_ldkc
!
      do c = 1, wf%n_v ! Use as though "b" for u_b_kdl term 
         do k = 1, wf%n_o
            do d = 1, wf%n_v
               do l = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  kdl  = index_three(k, d, l, wf%n_o, wf%n_v)
!
                  ld   = index_two(l, d, wf%n_o)
                  kc   = index_two(k, c, wf%n_o)
!
                  cl   = index_two(c, l, wf%n_v)
                  ck   = index_two(c, k, wf%n_v)
                  dl   = index_two(d, l, wf%n_v)
                  dk   = index_two(d, k, wf%n_v)
!
                  ckdl = index_packed(ck, dl)
                  cldk = index_packed(cl, dk)
!
!                 Set the values of u_b_kdl and g_kdl_c
!
                  u_b_kdl(c, kdl) = two*(wf%t2am(ckdl, 1)) - wf%t2am(cldk, 1)
!
                  g_kdl_c(kdl, c) = g_ld_kc(ld, kc) ! g_ldkc 
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the unordered integrals g_ld_kc = g_ldkc
!
!        Note: It is better to reorder g_kdl_c to g_k_dlc in the 
!        calculation of the E2.2 term (less memory). For now (8 Mar 2017), 
!        I'll just keep it simple & stupid.
!
      call deallocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Allocate the intermediate X_b_c = F_bc - sum_dkl g_ldkc u_kl^bd and set to zero
!
      call allocator(X_b_c, wf%n_v, wf%n_v)
  !    X_b_c = zero 
!
!     Copy the virtual-virtual Fock matrix into the intermediate 
!
      call dcopy((wf%n_v)**2, wf%fock_ab, 1, X_b_c, 1) ! X_b_c = F_bc 
!
!     Add the second contribution, 
!     - sum_dkl g_ldkc u_kl^bd = - sum_dkl u_b_kdl * g_kdl_c, to X_b_c
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  u_b_kdl,              &
                  wf%n_v,               &
                  g_kdl_c,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  X_b_c,                &
                  wf%n_v)
!
!     Deallocate u_b_kdl and g_kdl_c
!
      call deallocator(u_b_kdl, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      call deallocator(g_kdl_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
!     Allocate t_c_jai = t_ij^ac and set to zero
!
      call allocator(t_c_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      call dzero(t_c_jai, (wf%n_v)*(wf%n_v)*(wf%n_o)**2)
!
!     Determine t_c_jai = t_ij^ac 
!
      do c = 1, wf%n_v
         do j = 1, wf%n_o
            do a = 1, wf%n_v
               do i = 1, wf%n_o
!
!                 Calculate the necessary indices
!
                  jai  = index_three(j, a, i, wf%n_o, wf%n_v)
!
                  ai   = index_two(a, i, wf%n_v)
                  cj   = index_two(c, j, wf%n_v)
!
                  aicj = index_packed(ai, cj)
!
!                 Set the value of t_c_jai 
!
                  t_c_jai(c, jai) = wf%t2am(aicj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the E2.1 term and set it to zero
!
      call allocator(omega2_b_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
    ! omega2_b_jai = zero 
!
!     Calculate the E2.1 term 
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  X_b_c,                &
                  wf%n_v,               &
                  t_c_jai,              &
                  wf%n_v,               &
                  zero,                 &
                  omega2_b_jai,         &
                  wf%n_v)
!
!     Add the E2.1 term to the omega vector 
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  jai  = index_three(j, a, i, wf%n_o, wf%n_v)
                  ibj  = index_three(i, b, j, wf%n_o, wf%n_v)
!
                  ai   = index_two(a, i, wf%n_v)
                  bj   = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai,bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
!
                     wf%omega2(aibj,1) = wf%omega2(aibj,1) + omega2_b_jai(b,jai) &
                                                            + omega2_b_jai(a,ibj)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the E2.1 term, the X intermediate, and the reordered amplitudes 
!
      call deallocator(omega2_b_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      call deallocator(X_b_c, wf%n_v, wf%n_v)
      call deallocator(t_c_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     Print the omega vector, having added the E2.1 term to it
!
      if (debug) then
!
         write(unit_output,*) 
         write(unit_output,*) 'Omega(aibj,1) after E2.1 term has been added:'
         write(unit_output,*)
!
         call vec_print(wf%omega2, wf%n_t2am, 1)
!
      endif 
!
!     Allocate the Cholesky vector L_kc_J = L_kc^J and set it to zero 
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
    !  L_kc_J = zero
!
!     Read the Cholesky vector from file 
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
    !  g_ld_kc = zero
!
!     Calculate g_ld_kc = sum_J L_ld^J L_kc^J 
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ld_kc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate the Cholesky vector L_kc_J
!
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate g_k_dlc = g_ldkc and set to zero 
!
      call allocator(g_k_dlc, wf%n_o, (wf%n_o)*(wf%n_v)**2)
      call dzero(g_k_dlc, (wf%n_o)*(wf%n_o)*(wf%n_v)**2)
!
!     Allocate u_dlc_j = u_lj^dc and set to zero
!
      call allocator(u_dlc_j, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call dzero(u_dlc_j, (wf%n_o)*(wf%n_o)*(wf%n_v)**2)
!
!     Determine g_k_dlc = g_ldkc and u_dlc_j = u_lj^dc 
!
      do k = 1, wf%n_o ! Use as though "j" for u_dlc_j term 
         do d = 1, wf%n_v
            do l = 1, wf%n_o
               do c = 1, wf%n_v
!
!                 Calculate the necessary indices 
!
                  dlc  = index_three(d, l, c, wf%n_v, wf%n_o)
!
                  ld   = index_two(l, d, wf%n_o)
                  kc   = index_two(k, c, wf%n_o)
!
                  dl   = index_two(d, l, wf%n_v)
                  ck   = index_two(c, k, wf%n_v)
!
                  dlck = index_packed(dl, ck)
!
                  dk   = index_two(d, k, wf%n_v)
                  cl   = index_two(c, l, wf%n_v)
!
                  dkcl = index_packed(dk, cl)
!
!                 Set the value of g_k_dlc and u_dlc_j 
!
                  g_k_dlc(k, dlc) = g_ld_kc(ld, kc)                           ! g_ldkc 
                  u_dlc_j(dlc, k) = two*(wf%t2am(dlck, 1)) - wf%t2am(dkcl, 1) ! u_lk^dc = 2 * t_lk^dc - t_kl^dc 
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the integrals g_ld_kc = g_ldkc 
!
      call deallocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Allocate the intermediate Y_k_j = F_kj  + sum_cdl u_lj^dc g_ldkc 
!                                     = F_k_j + sum_cdl g_k_dlc * u_dlc_j, and set it to zero
!
      call allocator(Y_k_j, wf%n_o, wf%n_o)
    !  Y_k_j = zero 
!
!     Copy the occupied-occupied Fock matrix, such that Y_k_j = F_kj 
!
      call dcopy((wf%n_o)**2, wf%fock_ij, 1, Y_k_j, 1)
!
!     Add sum_cdl g_k_dlc u_dlc_j to Y_k_j, such that 
!     Y_k_j = F_k_j + sum_cdl g_k_dlc u_dlc_j
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  g_k_dlc,              &
                  wf%n_o,               &
                  u_dlc_j,              &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  Y_k_j,                &
                  wf%n_o)
!
!     Deallocate u_dlc_j and g_k_dlc 
!
      call deallocator(u_dlc_j, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call deallocator(g_k_dlc, wf%n_o, (wf%n_o)*(wf%n_v)**2)
!
!     Allocate t_aib_k = t_ik^ab and set it to zero 
!
      call allocator(t_aib_k, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call dzero(t_aib_k, (wf%n_o)*(wf%n_o)*(wf%n_v)**2)
!
!     Determine t_aib_k = t_ik^ab 
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do k = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  aib  = index_three(a, i, b, wf%n_v, wf%n_o)
!
                  ai   = index_two(a, i, wf%n_v)
                  bk   = index_two(b, k, wf%n_v)
!
                  aibk = index_packed(ai, bk)
!
!                 Set the value of t_aib_k 
!
                  t_aib_k(aib, k) = wf%t2am(aibk, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the E2.2 term and set to zero 
!
      call allocator(omega2_aib_j, (wf%n_o)*(wf%n_v)**2, wf%n_o)
   !   omega2_aib_j = zero
!
!     Calculate the E2.2 term, 
!     - sum_k t_aib_k Y_k_j = - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc)
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  t_aib_k,              &
                  (wf%n_o)*(wf%n_v)**2, &
                  Y_k_j,                &
                  wf%n_o,               &
                  zero,                 &
                  omega2_aib_j,         &
                  (wf%n_o)*(wf%n_v)**2)
!
!     Deallocate t_aib_k and Y_k_j 
!
      call deallocator(t_aib_k, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call deallocator(Y_k_j, wf%n_o, wf%n_o)
!
!     Add the E2.2 term to the omega vector 
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a, i, wf%n_v)
                  bj   = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  aib  = index_three(a, i, b, wf%n_v, wf%n_o)
                  bja  = index_three(b, j, a, wf%n_v, wf%n_o) 
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + omega2_aib_j(aib, j) & 
                                                               + omega2_aib_j(bja, i)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the E2.2 term 
!
      call deallocator(omega2_aib_j, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     Print the omega vector, having added both the E2.1 and E2.2 terms
!
      if (debug) then 
!
         write(unit_output,*) 
         write(unit_output,*) 'Omega(aibj,1) after E2.2 term has been added:'
         write(unit_output,*)
!
         call vec_print(wf%omega2, wf%n_t2am, 1)
!
      endif 
!
   end subroutine omega_e2_ccsd
!
! 
end submodule omega