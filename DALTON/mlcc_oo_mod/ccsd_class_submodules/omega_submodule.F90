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
end submodule omega