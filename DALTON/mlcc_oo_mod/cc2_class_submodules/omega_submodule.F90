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
  subroutine construct_omega_cc2(wf)
!
!     Construct Omega (CC2)
!     Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!
!     
!
      implicit none 
!
      class(cc2) :: wf
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_bj_J
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_ia_bj ! = g_aibj
!
!     t2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_ia_bj ! = g_aibj/(e_a + e_b - e_i - e_j)
!
!     Batching variables
!  
      integer(i15) :: a_batch, a_first, a_last, a_length
      integer(i15) :: required, available, n_batch, batch_dimension, max_batch_length
!
!     Indices
!  
      integer(i15) :: a = 0, b = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: bj = 0, ai = 0, ia = 0
!
!     Timing variables
!
      real(dp) :: omega_start = zero
      real(dp) :: omega_end   = zero
!
!     Start timings of omega
!
      call cpu_time(omega_start)
!
!     Set the omega vector to zero 
!
      call dzero(wf%omega1,(wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Prepare for batching over index a (Assumes A1 requires the most memory)
!  
      required = max((2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                       !    
               + 2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                            ! Needed for g_aibj  
               + 2*((wf%n_v)**2)*(wf%n_J) + ((wf%n_o)**2)*(wf%n_J) &       ! and 't2' amplitudes  
               + 2*(wf%n_v)**2*(wf%n_o)**2), &                             !
                 ((wf%n_v)**2*(wf%n_o)**2   &                              ! Needed for A1:   
               + (wf%n_v)*(wf%n_o)*(wf%n_J) &                              ! Needed for L_kc_J in A1
               + ((wf%n_v)**2)*((wf%n_o)**2) &                             ! Needed for u_ckd_i in A1
               + 2*((wf%n_v)**2)*(wf%n_J) + 2*(wf%n_v)*(wf%n_o)*(wf%n_J) & ! Needed for L_da_J in A1
               + 2*((wf%n_v)**3)*(wf%n_o)))                                ! Needed for g_da_kc and reordering in A1
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
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
         a_length = a_last - a_first + 1 
!
!        Allocate L_bj_J and L_ia_J (= reordering of L_bj_J constrained to the batch)
!
         call allocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call allocator(L_ia_J, a_length*(wf%n_o), wf%n_J)
         call dzero(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call dzero(L_ia_J, a_length*(wf%n_o), wf%n_J)
!
         call wf%get_cholesky_ai(L_bj_J)
!
!        Create L_ia_J
!
         do a = 1, a_length
            do i = 1, wf%n_v
               do J = 1, wf%n_J
!
!                 Calculate compound indices
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
!        Allocate g_ia_bj
!
         call allocator(g_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!        Construct integral g_ia_bj (= g_aibj for the batch)
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
!        L_bj_J and L_ia_J
!
         call deallocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call deallocator(L_ia_J, a_length*(wf%n_o), wf%n_J)
!
!        Allocate t_ia_bj
!
         call allocator(t_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!        Create t2 amplitudes
!
         do i = 1, wf%n_o
            do a = 1, a_length
!
!              Calculate compound index
!
               ia = index_two(i, a, wf%n_o)
!
               do b = 1, wf%n_v
                  do j = 1, wf%n_o
!
!                    Calculare compond index               
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
!        Deallocate g_ia_bj
!
         call deallocator(g_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!        Construct singles contributions to batch
!
         call wf%omega_a1(t_ia_bj, a_first, a_last, a_length)
!
         call wf%omega_b1(t_ia_bj, a_first, a_last, a_length)
!
         call wf%omega_c1(t_ia_bj, a_first, a_last, a_length)

         call wf%omega_d1()
!
!        Deallocate t_ia_bj
!
         call deallocator(t_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
      enddo
!
!     Timings
!
      call cpu_time(omega_end)
      write(unit_output,*)'Time in omega:', omega_end-omega_start    
!
   end subroutine construct_omega_cc2
!
!
   subroutine omega_a1_cc2(wf, t_kc_di, c_first, c_last, c_length)
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
      class(cc2) :: wf
!
      integer(i15) :: c_first, c_last, c_length
!
      real(dp), dimension(c_length*(wf%n_o),(wf%n_v)*(wf%n_o)) :: t_kc_di
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_batch_length = 0, batch_dimension = 0
      integer(i15) :: n_batch = 0, a_first = 0, a_last = 0, a_batch = 0, a_length = 0
!
!     Indices 
!
      integer(i15) :: i = 0, j = 0, k = 0, a = 0, c = 0, d = 0
!
      integer(i15) :: ad = 0, da = 0
      integer(i15) :: ci = 0, ck = 0, di = 0, dk = 0
      integer(i15) :: kc = 0, ic = 0 
!
      integer(i15) :: dkc = 0 
      integer(i15) :: cidk = 0, ckdi = 0
!
      integer(i15) :: ad_dim = 0   
!
      real(dp), dimension(:,:), allocatable :: L_kc_J 
      real(dp), dimension(:,:), allocatable :: L_da_J  ! L_ad^J; a is being batched over
      real(dp), dimension(:,:), allocatable :: g_da_kc ! g_adkc; a is being batched over
      real(dp), dimension(:,:), allocatable :: g_a_dkc ! reordered g_adkc; a is being batched over
      real(dp), dimension(:,:), allocatable :: u_dkc_i ! u_ki^cd
!
      logical :: reorder ! To get L_ab_J reordered, for batching over a (first index)
!
!
!     Allocate u_dkc_i = u_ki^cd
!
      call allocator(u_dkc_i, (wf%n_o)*(wf%n_v)*(c_length), wf%n_o)
      call dzero(u_dkc_i, (wf%n_o)*(wf%n_v)*(c_length), wf%n_o)
!
!     Calculate u_ckd_i
!
      do c = 1, c_length
         do k = 1, wf%n_o
            do d = 1, wf%n_v
               do i = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  dkc = index_three(d, k, c, wf%n_v, wf%n_o)
!
                  kc = index_two(k, c, wf%n_o)
                  dk = index_two(d, k, wf%n_v)
                  ic = index_two(i, c, wf%n_o)
                  di = index_two(d, i, wf%n_v) 
!
!                 Calculate u_dkc_i
!
                  u_dkc_i(dkc, i) = two*(t_kc_di(kc, di)) - t_kc_di(ic, dk)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate Cholesky vector L_kc_J for all c
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call dzero(L_kc_J, (wf%n_o)*(wf%n_v)*(wf%n_J))
!
!     Get Cholesky vector L_kc_J for all c
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Start batching over a
!
      available = get_available()
      required  = ((wf%n_v)**2)*(wf%n_J) &
                  + (wf%n_v**2)*(wf%n_o)*(wf%n_v) & 
                  + 2*((wf%n_v)**3)*(wf%n_o)
!
      required = 4*required
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
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
         a_length = a_last - a_first + 1 
!
!        Allocate the Cholesky vector L_da_J = L_ad^J
!
         ad_dim = a_length*(wf%n_v) ! Dimension of ad for the batch over index a 
!
         call allocator(L_da_J, ad_dim, wf%n_J)
         call dzero(L_da_J, ad_dim*(wf%n_J))
!
!        Get reordered Cholesky vector L_da_J = L_ad^J 
!
         reorder = .true.
         call wf%get_cholesky_ab(L_da_J, a_first, a_last, ad_dim, reorder)
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
!        Deallocate L_da_J
!
         call allocator(L_da_J, ad_dim, wf%n_J)
!
!        Allocate g_a_ckd = g_adkc and set to zero Batching both c and a
!
         call allocator(g_a_dkc, a_length, (wf%n_o)*(wf%n_v)*c_length)
         call dzero(g_a_dkc, a_length, (wf%n_o)*(wf%n_v)*c_length)
!
!        Reorder the integrals to g_a_ckd
!
         do c = 1, c_length
            do k = 1, wf%n_o
               do d = 1, wf%n_v
                  do a = 1, a_length
!
!                    Get the needed indices 
!
                     kc  = index_two(k, c + c_first - 1, wf%n_o)
                     da  = index_two(d, a, wf%n_v)
                     dkc = index_three(d, k, c, wf%n_v, wf%n_o) 
!
!                    Set the value of reordered integral, g_a_ckd
!
                     g_a_dkc(a, dkc) = g_da_kc(da, kc) ! g_adkc 
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate g_da_kc = g_adkc
!
         call deallocator(g_da_kc, ad_dim, (wf%n_o)*(wf%n_v))
!
!        Calculate the A1 term (sum_ckd g_a,ckd * u_ckd,i) & add to the omega vector
! 
         call dgemm('N','N',                     &
                     a_length,                   &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)*c_length, &
                     one,                        &
                     g_a_dkc,                    &
                     a_length,                   &
                     u_dkc_i,                    &
                     (wf%n_o)*(wf%n_v)*c_length, &
                     one,                        &
                     wf%omega1(a_first,1),       &
                     wf%n_v)
!
!        Deallocate g_a_dkc
!  
         call allocator(g_a_dkc, a_length, (wf%n_o)*(wf%n_v)*c_length)

      enddo  ! End of batches of the index a 
!
!     Deallocate u_ckd_i
!
      call deallocator(u_dkc_i, (wf%n_o)*(wf%n_v)*(c_length), wf%n_o)
!
   end subroutine omega_a1_cc2
!
!
   subroutine omega_b1_cc2(wf, t_lc_ak, c_first, c_last, c_length)
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
      class(cc2) :: wf 
!
      integer(i15) :: c_first, c_last, c_length
!
      real(dp), dimension(c_length*(wf%n_o),(wf%n_v)*(wf%n_o)) :: t_lc_ak
!
      integer(i15) :: a = 0, c = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ai = 0, ak = 0, al = 0
      integer(i15) :: kc = 0, lc = 0
      integer(i15) :: ki = 0
!
      integer(i15) :: klc = 0
      integer(i15) :: lc_batch = 0
!
      integer(i15) :: akcl = 0,  alck = 0
   
!
      real(dp), dimension(:,:), allocatable :: L_ki_J        ! L_ki^J 
      real(dp), dimension(:,:), allocatable :: L_lc_J        ! L_lc^J 
      real(dp), dimension(:,:), allocatable :: L_lc_J_batch  ! L_lc^J  with c constrained to the batch
      real(dp), dimension(:,:), allocatable :: g_ki_lc       ! g_kilc 
      real(dp), dimension(:,:), allocatable :: g_klc_i       ! g_kilc 
      real(dp), dimension(:,:), allocatable :: u_a_klc       ! u_kl^ac = 2 t_kl^ac - t_lk^ac
!
!
!     Allocate Cholesky vectors L_lc_J and L_lc_J_batch 
!
      call allocator(L_lc_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call allocator(L_lc_J_batch, (wf%n_o)*c_length, wf%n_J)
!
      call dzero(L_lc_J, (wf%n_o)*(wf%n_v)*(wf%n_J))
      call dzero(L_lc_J_batch, (wf%n_o)*c_length, wf%n_J)

!
!     Read the Cholesky vectors L_lc_J
!
      call wf%get_cholesky_ia(L_lc_J)
!
!     Constrain Cholesky vectors to c     
!
      do j = 1, wf%n_J
         do c = 1, c_length
            do l = 1, wf%n_o
!
!              Calculate compound indices
!
               lc_batch = index_two(l, c, wf%n_o)
               lc = index_two(l, c + c_first - 1, wf%n_o)
!
               L_lc_J_batch(lc_batch, J) = L_lc_J(lc, J) 
!
            enddo
         enddo
      enddo
!
!     Deallocate L_lc_J
!
      call deallocator(L_lc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate L_ki_J
!
      call allocator(L_ki_J, (wf%n_o)**2, wf%n_J)
      call dzero(L_ki_J, (wf%n_J)*(wf%n_o)**2)
!
!     Read the Cholesky vectors L_ki_J
!
      call wf%get_cholesky_ij(L_ki_J)
!
!
!     Allocate integrals g_ki_lc = g_kilc
!
      call allocator(g_ki_lc, (wf%n_o)**2, (wf%n_o)*c_length)
!
!     Calculate g_ki_lc = sum_J L_ki_J L_lc_J^T 
! 
      call dgemm('N','T',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*c_length, &
                  wf%n_J,            &
                  one,               &
                  L_ki_J,            &
                  (wf%n_o)**2,       &
                  L_lc_J,            &
                  (wf%n_o)*c_length, &
                  zero,              &
                  g_ki_lc,           &
                  (wf%n_o)**2)
!
      call deallocator(L_ki_J, (wf%n_o)**2, wf%n_J)
      call deallocator(L_lc_J_batch, (wf%n_o)*c_length, wf%n_J)
!
!
!     Allocate reordered integrals g_ckl_i = g_kilc 
!
      call allocator(g_klc_i, c_length*((wf%n_o)**2), wf%n_o)
      call dzero(g_klc_i, c_length*((wf%n_o)**3))
!
!     Determine g_ckl_i = g_kilc 
!
      do c = 1, c_length
         do k = 1, wf%n_o
            do l = 1, wf%n_o
               do i = 1, wf%n_o
!
!                 Calculate necessary indices
!
                  klc = index_three(k, l, c, wf%n_o, wf%n_o) 
!
                  ki  = index_two(k, i, wf%n_o)                 
                  lc  = index_two(l, c, wf%n_o)           
!
!                 Set value of g_ckl_i
!
                  g_klc_i(klc, i) = g_ki_lc(ki, lc)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate unordered integrals g_ki_lc
!
      call deallocator(g_ki_lc, (wf%n_o)**2, (wf%n_o)*c_length)
!
!     Allocate redordered u_a_ckl = u_kl^ac and set it to zero
!
      call allocator(u_a_klc, wf%n_v, c_length*(wf%n_o)**2)
      call dzero(u_a_klc, (wf%n_v),c_length*(wf%n_o)**2)
!
!     Determine u_a_ckl = u_kl^ac
!
      do a = 1, wf%n_v
         do c = 1, c_length
            do k = 1, wf%n_o
               do l = 1, wf%n_o
!
!                 Calculate necessary indices
!
                  klc  = index_three(k, l, c, wf%n_o, wf%n_o) 
!
                  ak   = index_two(a, k, wf%n_v)
                  lc   = index_two(l, c, wf%n_o)
!
                  al   = index_two(a, l, wf%n_v)
                  kc   = index_two(k, c, wf%n_o)
!
!
!                 Set the value of u_a_ckl = u_kl^ac = 2*t_kl^ac - t_lk^ac = 2*t_ak,cl - t_al,ck 
!
                  u_a_klc(a, klc) = two*(t_lc_ak(lc, ak)) - t_lc_ak(kc, al)
!                  
               enddo
            enddo
         enddo
      enddo
!
!     Calculate the B1 term, - sum_ckl u_a_klc g_klc_i
!
      call dgemm('N','N',                 &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  c_length*((wf%n_o)**2), &
                  -one,                   &
                  u_a_klc,                &
                  wf%n_v,                 &
                  g_klc_i,                &
                  c_length*((wf%n_o)**2), &
                  one,                    &
                  wf%omega1,              &
                  wf%n_v)
!
!     Deallocate remaining vectors 
!
      call deallocator(u_a_klc, wf%n_v, c_length*((wf%n_o)**2))
      call deallocator(g_klc_i, c_length*((wf%n_o)**2), wf%n_o) 
!
   end subroutine omega_b1_cc2
!
!
   subroutine omega_c1_cc2(wf, t_kc_ai, c_first, c_last, c_length)        
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
      class(cc2) :: wf 
!
      integer(i15) :: c_first, c_last, c_length
!
      real(dp), dimension(c_length*(wf%n_o),(wf%n_v)*(wf%n_o)) :: t_kc_ai
!
      real(dp), dimension(:,:), allocatable :: F_kc     
      real(dp), dimension(:,:), allocatable :: u_ai_kc  
      real(dp), dimension(:,:), allocatable :: omega1_ai
!
      integer(i15) :: i = 0, k = 0, c = 0, a = 0
!
      integer(i15) :: kc = 0, ai = 0, ak = 0, ic = 0
!
!     Allocation of F_ck, u_ai_ck and omega1_ai
!
      call allocator(F_kc, (wf%n_o)*c_length, 1)
      call allocator(u_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*c_length)  
      call allocator(omega1_ai, (wf%n_o)*(wf%n_v), 1)
!
      call dzero(F_kc, (wf%n_o)*c_length)
      call dzero(u_ai_kc, ((wf%n_o)**2)*(wf%n_v)*c_length)
!
!     Set up u_ck_ai and virtual-occupied Fock matrix
!
      do k = 1, wf%n_o
         do c = 1, c_length
!
!           Set up compound index
!  
            kc = index_two(k, c, wf%n_o)
!
!           Reorder MO Fock matrix
!
            F_kc(kc, 1) = wf%fock_ia(k, c)
!
            do a = 1, wf%n_v
               do i = 1, wf%n_o
!
!                 Set up compound indices
!
                  ai = index_two(a, i, wf%n_v)
                  ic = index_two(i, c, wf%n_o)
                  ak = index_two(a, k, wf%n_v)
!                    
!                 u_ck_ai
!
                  u_ai_kc(ai, kc) = two*(t_kc_ai(kc, ai)) - t_kc_ai(ic, ak)
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
                  (wf%n_o)*c_length, &
                  one,               &
                  u_ai_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  F_kc,              &
                  (wf%n_o)*c_length, &
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
!     Deallocation of F_ck, u_ai_ck and omega1_ai
!
      call deallocator(F_kc, (wf%n_o)*c_length, 1)
      call deallocator(u_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*c_length)  
      call deallocator(omega1_ai, (wf%n_o)*(wf%n_v), 1)
!
   end subroutine omega_c1_cc2
!
!
   subroutine omega_d1_cc2(wf)
!
!     D1 omega term: Omega_ai^D1=F_ai_T1
!
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mars 2017
!
      implicit none
!
      class(cc2) :: wf
!
!     Add F_a_i to omega
!
      call daxpy((wf%n_o)*(wf%n_v), one, wf%fock_ai, 1, wf%omega1, 1)
!
   end subroutine omega_d1_cc2
!
!
end submodule