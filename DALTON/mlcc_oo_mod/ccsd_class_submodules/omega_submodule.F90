submodule (ccsd_class) omega
!
!
!                       -::- Omega submodule (CCSD) -::-
!           Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!
!     Contains the following family of procedures of the CCSD class:
!
!     construct_omega(wfn): constructs the projection vector for the current amplitudes
!
contains
!
!
   subroutine construct_omega_cc_singles_doubles(wf)
!
!      Construct Omega (CCSD)
!      Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!
!      Directs the calculation of the omega vector
!
       implicit none 
!
       class(cc_singles_doubles) :: wf
!
!      Construct singles contributions 
!
       ! call wf%omega_a1
       ! call wf%omega_b1
       ! call wf%omega_c1
       ! call wf%omega_d1
!
!      Construct doubles contributions 
!
       ! call wf%omega_a2
       ! call wf%omega_b2
       ! call wf%omega_c2
       ! call wf%omega_d2
       ! call wf%omega_e2
!
   end subroutine construct_omega_cc_singles_doubles
!
!
   subroutine omega_a1_cc_singles_doubles(wf)
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
      class(cc_singles_doubles) :: wf
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
      logical :: debug = .false.
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
      L_kc_J = zero
!
!     Get Cholesky vector L_kc_J
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate u_ckd_i = u_ki^cd
!
      call allocator(u_ckd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      u_ckd_i = zero
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
!     for which we need to have enough room to store L_ad_J and g_ad_kc, and, later on in the same loop, 
!     g_ad_kc and g_a_ckd simultaneously
!
      available = get_available()
      required  = max(((wf%n_v)**2)*(wf%n_J) + (wf%n_v**2)*(wf%n_o)*(wf%n_v), &
                              2*((wf%n_v)**3)*(wf%n_o)) ! Eirik: I am not sure if this is an accurate estimate of the required memory
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
         L_da_J = zero
!
!        Get reordered Cholesky vector L_da_J = L_ad^J 
!
         reorder = .true.
         call wf%get_cholesky_ab(L_da_J,a_begin,a_end,ad_dim,reorder)
!
!        Allocate g_da_kc = g_adkc and set to zero
!
         call allocator(g_da_kc, ad_dim, (wf%n_o)*(wf%n_v))
         g_da_kc = zero 
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
         g_a_ckd = zero
!
!        Reorder the integrals to g_a_ckd
!
         do a = 1, batch_length
            do c = 1, wf%n_v
               do k = 1, wf%n_o
                  do d = 1, wf%n_v
!
!                    Calculate the necessary indices
!
                     da  = index_two(d, a, wf%n_v)
                     kc  = index_two(k, c, wf%n_o)
!
                     ckd = index_three(c, k, d, wf%n_v, wf%n_o) 
!
!                    Set the value of g_a_ckd
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
         call dgemm('N','N',&
                     batch_length,&
                     wf%n_o,&
                     (wf%n_o)*(wf%n_v)**2,&
                     one,&
                     g_a_ckd,&
                     batch_length,&
                     u_ckd_i,&
                     (wf%n_o)*(wf%n_v)**2,&
                     one,&
                     wf%omega1(a_begin,1),&
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
   end subroutine omega_a1_cc_singles_doubles
!
!
   subroutine omega_b1_cc_singles_doubles(wf)
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
      class(cc_singles_doubles) :: wf 
!
      logical :: debug = .false.
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
      L_ki_J = zero
      L_lc_J = zero
!
!     Read the Cholesky vectors L_ki_J and L_lc_J
!
      call wf%get_cholesky_ij(L_ki_J)
      call wf%get_cholesky_ia(L_lc_J)
!
!     Allocate integrals g_ki_lc = g_kilc
!
      call allocator(g_ki_lc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      g_ki_lc = zero
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
      call allocator(g_ckl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
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
      u_a_ckl = zero
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
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  u_a_ckl,              &
                  wf%n_v,               &
                  g_ckl_i,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  wf%omega1,            &
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
      call deallocator(u_a_ckl, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      call deallocator(g_ckl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
!
   end subroutine omega_b1_cc_singles_doubles
!
!
   subroutine omega_e2_cc_singles_doubles(wf)
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
      class(cc_singles_doubles) :: wf 
!
      logical :: debug = .false.
!
!     Indices 
!
      integer(i15) :: aib = 0, aibk = 0, bk = 0, bja = 0, ibj = 0, aibj = 0, dlck = 0
      integer(i15) :: b = 0, c = 0, k = 0, d = 0, ck = 0, ckdl = 0, cl = 0, cldk = 0
      integer(i15) :: dk = 0, dl = 0, kc = 0, kdl = 0, l = 0, ld = 0, a = 0, ai = 0, 
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
      L_kc_J = zero
!
!     Read the Cholesky vector from file 
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ld_kc = zero
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
      call allocator(u_b_kdl, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      u_b_kdl = zero
!
!     Allocate g_kdl_c = g_ldkc and set to zero 
!
      call allocator(g_kdl_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      g_kdl_c = zero
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
      X_b_c = zero 
!
!     Copy the virtual-virtual Fock matrix into the intermediate 
!
      call dcopy((wf%n_v)**2, F_a_b, 1, X_b_c, 1) ! X_b_c = F_bc 
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
      t_c_jai = zero 
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
      omega2_b_jai = zero 
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
         call vec_print_packed(omega2, n_t2am)
!
      endif 
!
!     Allocate the Cholesky vector L_kc_J = L_kc^J and set it to zero 
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_kc_J = zero
!
!     Read the Cholesky vector from file 
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_ld_kc,n_ov,n_ov)
      g_ld_kc = zero
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
      g_k_dlc = zero 
!
!     Allocate u_dlc_j = u_lj^dc and set to zero
!
      call allocator(u_dlc_j, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      u_dlc_j = zero 
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
      Y_k_j = zero 
!
!     Copy the occupied-occupied Fock matrix, such that Y_k_j = F_kj 
!
      call dcopy((wf%n_o)**2, F_i_j, 1, Y_k_j, 1)
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
      t_aib_k = zero
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
      omega2_aib_j = zero
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
         call vec_print_packed(wf%omega2, n_t2am)
!
      endif 
!
   end subroutine omega_e2_cc_singles_doubles
!
!
   subroutine omega_d2_cc_singles_doubles(wf)
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
!     The first term is referred to as the D2.1 term, and comes out ordered as (ai,bj) 
!     The second term is referred to as the D2.2 term, and comes out ordered as (ai,bj)
!     The third term is referred to as the D2.3 term, and comes out ordered as (ai,bj)
!
!     All terms are added to the omega vector element omega2(ai,bj) of the 
!     wavefunction object wf
!
!     The routine adds the terms in the following order: D2.3, D2.1, D2.2
!
      implicit none 
!
      class(cc_singles_doubles) :: wf 
!
      logical :: debug = .false.
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
      real(dp), dimension(:,:), allocatable :: omega2_ai_bj ! For storing the D2.3, D2.2 & D2.1 terms temporarily
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
!     Allocate the Cholesky vector L_kc_J = L_kc^J and set to zero 
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_kc_J = zero
!
!     Read the Cholesky vector L_kc_J from file
!
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ld_kc = zero 
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
      L_ld_kc = zero 
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
      u_ai_ld = zero 
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
      Z_ai_kc = zero
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
      omega2_ai_bj = zero
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
                     omega2(aibj, 1) = omega2(aibj, 1) + omega2_ai_bj(ai, bj) & 
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
         call vec_print_packed(wf%omega2, wf%n_t2am)
!
      endif                
!
!     Allocate the L_ai_J and L_kc_J terms and set them to zero 
!
      call allocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      L_ai_J = zero
      L_kc_J = zero
!
!     Read the Cholesky vectors from file 
!
      call wf%get_cholesky_ai(L_ai_J)
      call wf%get_cholesky_ia(L_kc_J)
!
!     Allocate g_ai_kc = g_aikc and set it zero
!
      call allocator(g_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ai_kc = zero 
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
      u_kc_bj = zero
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
      omega2_ai_bj = zero 
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
         call vec_print_packed(wf%omega2, wf%n_t2am)
!
      endif 
!
!     - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj 
!
!     Allocate L_ki_J and set it to zero
!
      call allocator(L_ki_J, (wf%n_o)**2, wf%n_J)
      L_ki_J = zero 
!
!     Read the Cholesky vector L_ki_J from file 
!
      call wf%get_cholesky_ij(L_ki_J)
!
!     Allocate the full g_ai_ck = g_acki and set it to zero 
!
      call allocator(g_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ai_ck = zero 
!
!     Prepare for batching over the index a to calculate g_ai_ck = g_acki
!
!        To calculate this term, we need to hold L_ac^J and g_acki
!        in memory simultaneously 
!
      required = (wf%n_J)*(wf%n_vir)**2 + ((wf%n_v)**2)*((wf%n_o)**2) ! Eirik: This estimate has to be updated,
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
         L_ca_J = zero
!
!        Read the Cholesky vector from file 
!
         call wf%get_cholesky_ab(L_ca_J, a_begin, a_end, ac_dim, .true.)
!
!        Allocate the integral g_ca_ki = g_acki and set to zero 
!
         call allocator(g_ca_ki, ac_dim, (wf%n_o)**2)
         g_ca_ki = zero
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
                     one,         &
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
      u_ck_bj = zero 
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
      omega2_ai_bj = zero
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
                  ai = index_two(a, i, n_vir)
                  bj = index_two(b, j, n_vir)
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
         call vec_print_packed(wf%omega2, wf%n_t2am)
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
   end subroutine omega_d2_cc_singles_doubles
!
!
end submodule omega