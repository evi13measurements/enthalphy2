module mlcc_omega
!
! 	MLCC Omega 
! 	Written by Eirik F. Kjønstad and Sarai Folkestad, 28 Feb 2017
!
! 	Routines for the calculation of the omega vector < mu | exp(-T) H exp(T) | R >
!  The routine mlcc_omega_calc directs the calculation and can be called from outside the module.
!
   use mlcc_data
   use mlcc_workspace
   use mlcc_utilities
   use mlcc_cholesky
!
contains 
   subroutine mlcc_omega_calc
!
!     MLCC Omega calculation
!     Written by Eirik F. Kjønstad and Sarai Folkestad, 28 Feb 2017
!
!     Controls the calculation of the omega vector
!
      implicit none
!
!     Add the singles contributions to < mu | exp(-T) H exp(T) | R >
!
      call mlcc_omega_a1
      call mlcc_omega_b1
      call mlcc_omega_c1
      call mlcc_omega_d1     
!
!     Add the doubles contributions to < mu | exp(-T) H exp(T) | R >
!
      call mlcc_omega_e2
      call mlcc_omega_d2
      call mlcc_omega_c2
      call mlcc_omega_a2
      call mlcc_omega_b2
!
   end subroutine mlcc_omega_calc
!
   subroutine mlcc_omega_a1
!
!     MLCC Omega A1 term ( sum_ckd g_adkc * u_ki^cd )
!     Written by Eirik F. Kjønstad and Sarai Folkestad, 28 Feb 2017
!
!     NB! Needs to be rewritten with T1 transformed integrals
!     eventually (makes no difference for MP2 guess)
!
      implicit none
!
      integer :: i,j,a
!
      integer :: required,available,max_batch_length,batch_dimension,n_batch
!
      integer :: a_begin,a_end,a_batch,batch_length
!
      integer :: ad,ad_dim,c,ci,cidk,ck,ckd,ckdi,di,dk,k,kc,d,da
!
      logical :: debug = .false.
!
      real(dp), dimension(:,:), pointer :: L_kc_J  => null()
      real(dp), dimension(:,:), pointer :: L_da_J  => null()   ! L_ad^J; a is being batched over
      real(dp), dimension(:,:), pointer :: g_da_kc => null()   ! g_adkc; a is being batched over
      real(dp), dimension(:,:), pointer :: g_a_ckd => null()   ! reordered g_adkc; a is being batched over
      real(dp), dimension(:,:), pointer :: u_ckd_i => null() 
!
!     Allocate Cholesky vector L_kc_J
!
      call allocator(L_kc_J,n_ov,n_J)
!
!     Set L_kc_J to zero
!
      L_kc_J = zero
!
!     Read Cholesky vector L_kc_J
!
      call read_cholesky_ia(L_kc_J)
!
!     Allocate u_ckd_i = u_ki^cd
!
      call allocator(u_ckd_i,n_ovv,n_occ)
!
!     Calculate u_ckd_i
!
      do c=1,n_vir
         do k=1,n_occ
            do d=1,n_vir
               do i=1,n_occ
!
!                 Calculate the necessary indices 
!
                  ckd  = index_three(c,k,d,n_vir,n_occ)
                  ck   = index_two(c,k,n_vir)
                  di   = index_two(d,i,n_vir)
                  ci   = index_two(c,i,n_vir)
                  dk   = index_two(d,k,n_vir)
                  ckdi = index_packed(ck,di)
                  cidk = index_packed(ci,dk)
!
!                 Calculate u_ckd_i
!
                  u_ckd_i(ckd,i) = two*t2am(ckdi,1)-t2am(cidk,1)
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
      required        = max(n_vir*n_vir*n_J + n_vir*n_vir*n_occ*n_vir,2*n_vir*n_vir*n_occ*n_vir) ! Eirik: I am not sure if this is an accurate estimate of the required memory
      available       = get_available()
      batch_dimension = n_vir ! Batch over the virtual index a
!
      write(luprint,*) 'Required',required
      write(luprint,*) 'Available',available
!
   !   if (debug) required = 1500000
!
      max_batch_length = 0 ! Initilization of unset variables 
      n_batch = 0
!
      call n_one_batch(required,available,max_batch_length,n_batch,batch_dimension)           
!
!     Loop over the number of a batches 
!
      do a_batch=1,n_batch
!
!        For each batch, get the limits for the a index 
!
         call one_batch_limits(a_begin,a_end,a_batch,max_batch_length,batch_dimension)
         batch_length = a_end - a_begin + 1 
!
!        Allocate the Cholesky vector L_da_J = L_ad^J
!
         ad_dim = batch_length*n_vir ! Dimension of ad for the batch over index a 
         call allocator(L_da_J,ad_dim,n_J)
!
!        Read in the reordered Cholesky vector L_da_J = L_ad^J 
!
         call read_cholesky_ab_reorder(L_da_J,a_begin,a_end,ad_dim)
!
!        Allocate g_da_kc = g_adkc and set to zero
!
         call allocator(g_da_kc,ad_dim,n_ov)
         g_da_kc = zero 
!
!        Calculate g_da_kc = sum_J L_da_J L_kc_J^T = sum_J L_ad^J L_kc^J = g_adkc 
!     
         call dgemm('N','T',ad_dim,n_ov,n_J,&
                     one,L_da_J,ad_dim,L_kc_J,n_ov,&
                     zero,g_da_kc,ad_dim) 
!
!        Deallocate the reordered Cholesky vector L_da_J
!
         call deallocator(L_da_J,ad_dim,n_J)
!
!        Allocate g_a_ckd = g_adkc and set to zero
!
         call allocator(g_a_ckd,batch_length,n_ovv)
         g_a_ckd = zero
!
!        Reorder the integrals to g_a_ckd
!
         do a=1,batch_length
            do d=1,n_vir
               do k=1,n_occ
                  do c=1,n_vir
!
!                    Calculate the necessary indices
!
                     da  = index_two(d,a,n_vir)
                     kc  = index_two(k,c,n_occ)
                     ckd = index_three(c,k,d,n_vir,n_occ) 
!
!                    Set the value of g_a_ckd
!
                     g_a_ckd(a,ckd) = g_da_kc(da,kc) ! g_adkc 
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate reordered integrals g_da_kc
!
         call deallocator(g_da_kc,ad_dim,n_ov)
!
!        Calculate the A1 term (sum_ckd g_a,ckd * u_ckd,i) & add to the omega vector
!
         call dgemm('N','N',batch_length,n_occ,n_ovv,&
                     one,g_a_ckd,batch_length,u_ckd_i,n_ovv,&
                     one,omega1(a_begin,1),n_vir)
!
      enddo ! End of batches of the index a 
!
!     Print the omega vector 
!
      if (debug) then  
         write(luprint,*) 
         write(luprint,*) 'Omega(a,i) after A1 term has been added:'
         write(luprint,*)
!
         call vec_print(omega1,n_vir,n_occ)
!
      endif
!
!     Deallocate vectors 
!
      call deallocator(u_ckd_i,n_ovv,n_occ)
      call deallocator(g_a_ckd,batch_length,n_ovv)
      call deallocator(L_kc_J,n_ov,n_J)
!
   end subroutine mlcc_omega_a1
!
   subroutine mlcc_omega_b1
!
!     MLCC Omega B1 term ( - sum_ckl u_kl^ac * g_kilc )
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, 28 Feb 2017
!
!     NB! Needs to be rewritten with T1 transformed integrals
!     eventually (makes no difference for MP2 guess)
!
      implicit none
!
      integer :: i,j
!
      logical :: debug = .true.
!
      double precision factor
!
      integer :: a,c,k,l,ckl,ki,ak,akcl,al,alck,ck,ai,cl,lc
!
      real(dp), dimension(:,:), pointer :: L_ki_J        => null() 
      real(dp), dimension(:,:), pointer :: L_lc_J        => null()
      real(dp), dimension(:,:), pointer :: g_ki_lc       => null()
      real(dp), dimension(:,:), pointer :: g_ckl_i       => null() ! Reordered two-electron integrals
      real(dp), dimension(:,:), pointer :: u_a_ckl       => null() ! Reordered u_kl^ac = 2 t_kl^ac - t_lk^ac
      real(dp), dimension(:,:), pointer :: b1_a_i        => null() 
!
!     I. Calculation and reordering of g_ki,lc = sum_J L_ki^J * L_lc^J 
!
!     Allocate Cholesky vectors L_ki,J and L_lc,J 
!
      call allocator(L_ki_J,n_oo,n_J)
      call allocator(L_lc_J,n_ov,n_J)
!
      L_ki_J = zero
      L_lc_J = zero
!
!     Read L_ki,J
!
      call read_cholesky_ij(L_ki_J)
!
!     Read L_lc,J
!
      call read_cholesky_ia(L_lc_J)
!
!     Allocate integrals g_ki_lc
!
      call allocator(g_ki_lc,n_oo,n_ov)
!
      g_ki_lc = zero
!
!     Calculate g_ki_lc = sum_J L_ki,J * L_lc,J^T 
! 
      call dgemm('N','T',n_oo,n_ov,n_J,&
                  one,L_ki_J,n_oo,L_lc_J,n_ov,&
                  zero,g_ki_lc,n_oo) 
!
!     Deallocate the Cholesky vectors 
!
      call deallocator(L_ki_J,n_oo,n_J)
      call deallocator(L_lc_J,n_ov,n_J)
!
!     Allocate reordered integrals g_ckl,i 
!
      call allocator(g_ckl_i,n_oov,n_occ)
!
!     Save reordered integrals g_ckl,i
!
      do c = 1,n_vir
         do k = 1,n_occ
            do l = 1,n_occ
               do i = 1,n_occ
!
!                 Calculate necessary indices
!
                  ckl = index_three(c,k,l,n_vir,n_occ) 
                  ki  = index_two(k,i,n_occ)                 
                  lc  = index_two(l,c,n_occ)           
!
!                 Set value of g_ckl_i
!
                  g_ckl_i(ckl,i) = g_ki_lc(ki,lc)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate unordered integrals g_ki_lc
!
      call deallocator(g_ki_lc,n_oo,n_ov)
!
!     Allocate redordered u amplitudes u_a,ckl 
!
      call allocator(u_a_ckl,n_vir,n_oov)
!
!     Set u_a,ckl to zero
!
      u_a_ckl = zero
!
!     Save reordered u_a,ckl 
!
      do c = 1,n_vir
         do k = 1,n_occ
            do l = 1,n_occ
               do a = 1,n_vir
!
!                 Calculate necessary indices
!
                  ckl  = index_three(c,k,l,n_vir,n_occ) 
!
                  ak   = index_two(a,k,n_vir)
                  cl   = index_two(c,l,n_vir)
                  akcl = index_packed(ak,cl)
!
                  al   = index_two(a,l,n_vir)
                  ck   = index_two(c,k,n_vir)
                  alck = index_packed(al,ck)
!
!                 Set value of u_a_ckl = u_kl^ac = 2*t_kl^ac - t_lk^ac = 2*t_ak,cl - t_al,ck 
!
                  u_a_ckl(a,ckl) = two*t2am(akcl,1) - t2am(alck,1)
!                  
               enddo
            enddo
         enddo
      enddo
!
!     Calculate the B1 term (- sum_ckl u_a_ckl * g_ckl_i)
!
      factor = -one 
      call dgemm('N','N',n_vir,n_occ,n_oov,&
                  factor,u_a_ckl,n_vir,g_ckl_i,n_oov,&
                  one,omega1,n_vir) 
!
!     Print the omega vector 
!
      if (debug) then  
         write(luprint,*) 
         write(luprint,*) 'Omega(a,i) after B1 term has been added:'
         write(luprint,*)
         call vec_print(omega1,n_vir,n_occ)
      endif
!
!     Deallocate remaining vectors 
!
      call deallocator(u_a_ckl,n_vir,n_oov)
      call deallocator(g_ckl_i,n_oov,n_occ)
!
   end subroutine mlcc_omega_b1
!
   subroutine mlcc_omega_c1
!
!  Purpose: Calculate C1 term of Omega
!  (Omega_ai^C1=sum_ck F_kc u_ai_ck)
!
!  Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mars 2017
!
   implicit none
!
   real(dp),dimension(:,:),pointer  :: F_ck => null()
   real(dp),dimension(:,:),pointer  :: u_ck_ai => null()
   integer                          :: i,k,c,a
   integer                          :: ck,ai,ak,ci,ckai,ciak
   logical                          :: debug = .true.
!
!
!  Allocation
!
   call allocator(u_ck_ai,n_ov,n_ov)   ! 2*t_ck_ai-t_ci_ak
   call allocator(F_ck,1,n_ov) ! T1-transformed F_k_c reordered
!
! Set up u_ck_ai and MO Fock matrix virtual-occupied
!
   do k = 1,n_occ
      do c = 1,n_vir
         ck = index_two(c,k,n_vir)
!
!        MO Fock matrix
!
         F_ck(1,ck) = F_i_a(k,c)
         do i = 1,n_occ
            do a = 1,n_vir
!
!              Necessary indices
!
               ai = index_two(a,i,n_vir)
               ci = index_two(c,i,n_vir)
               ak = index_two(a,k,n_vir)
!
               ckai = index_packed(ck,ai)
               ciak = index_packed(ci,ak)
!              
!              u_ck_ai
!
               u_ck_ai(ck,ai) = two*t2am(ckai,1)-t2am(ciak,1)
            enddo
         enddo
      enddo
   enddo
!
!  Matrix multiplication
!
   call dgemm('N','N',1,n_ov,n_ov &
      ,-one,F_ck,1,u_ck_ai,n_ov &
      ,one,omega1,n_ov)
!
!  Deallocation
!
  ! call deallocator(F_ck,n_ov,1) ! wrong order; I don't think this makes any difference...
  call deallocator(F_ck,1,n_ov)
   call deallocator(u_ck_ai,n_ov,n_ov)

   end subroutine mlcc_omega_c1
!
   subroutine mlcc_omega_d1
!
!  Purpose: Calculate C1 term of Omega
!  (Omega_ai^D1=F_ai_T1
!
!  Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mars 2017
!
   implicit none
!
!  Add F_a_i to omega
!
   call daxpy(n_ov,one,F_a_i,1,omega1,1)
!
   end subroutine mlcc_omega_d1
!
   subroutine mlcc_omega_e2
!
!     MLCC Omega E2 term
!     Written by Eirik F. Kjønstad and Sarai Folkestad, 7 Mar 2017
!
!     NB! Needs to be rewritten with T1 transformed integrals
!     eventually (makes no difference for MP2 guess)
!
!     This routine calculates the E2 term,
!
!        sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd) 
!        - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!
!     where
!
!        u_kl^bc = 2 * t_kl^bc - t_lk^bc.
!
!     The first term is referred to as the E2.1 term, and comes out ordered as (b,jai) 
!     The second term is referred to as the E2.2 term, and comes out ordered as (aib,j)
!
!     Both are added to the omega vector element omega2(ai,bj)
!
      implicit none 
!
      logical :: debug = .true.
!
      integer :: b,c,k,d,ck,ckdl,cl,cldk,dk,dl,kc,kdl,l,ld,a,ai,aibj,bj,aicj,cj,i,j,jai,dlc,dkcl,dlck,aib,aibk,bk
!
      real(dp), dimension(:,:), pointer :: omega2_b_jai => null() ! For storing the E2.1 term temporarily
      real(dp), dimension(:,:), pointer :: L_kc_J       => null() ! L_kc^J
      real(dp), dimension(:,:), pointer :: g_ld_kc      => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: g_kdl_c      => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: u_b_kdl      => null() ! u_kl^bd 
      real(dp), dimension(:,:), pointer :: X_b_c        => null() ! An intermediate, see below for definition
      real(dp), dimension(:,:), pointer :: t_c_jai      => null() ! t_ij^ac 
!
      real(dp), dimension(:,:), pointer :: g_k_dlc      => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: u_dlc_j      => null() ! u_lj^dc 
      real(dp), dimension(:,:), pointer :: omega2_aib_j => null() ! For storing the E2.2 term temporarily
      real(dp), dimension(:,:), pointer :: Y_k_j        => null() ! An intermediate, see below for definition 
      real(dp), dimension(:,:), pointer :: t_aib_k      => null() ! t_ik^ab 
!
!     Allocate the Cholesky vector L_kc_J = L_kc^J and set to zero 
!
      call allocator(L_kc_J,n_ov,n_J)
      L_kc_J = zero
!
!     Read the Cholesky vector from file 
!
      call read_cholesky_ia(L_kc_J)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_ld_kc,n_ov,n_ov)
      g_ld_kc = zero
!
!     Calculate g_ld_kc = sum_J L_ld^J L_kc^J 
!
            write(luprint,*) 'Lalala 1.1'
      call flshfo(luprint)
      call dgemm('N','T',n_ov,n_ov,n_J,&
                  one,L_kc_J,n_ov,L_kc_J,n_ov,&
                  zero,g_ld_kc,n_ov)
!
!     Deallocate the Cholesky vector L_kc_J
!
      call deallocator(L_kc_J,n_ov,n_J)
!
!     Allocate u_b_kdl = u_kl^bd and set to zero
!
      call allocator(u_b_kdl,n_vir,n_oov)
      u_b_kdl = zero
!
!     Allocate g_kdl_c = g_ldkc and set to zero 
!
      call allocator(g_kdl_c,n_oov,n_vir)
      g_kdl_c = zero
!
!     Determine u_b_kdl = u_kl^bd and g_kdl_c = g_ldkc
!
      do c = 1,n_vir ! Use as though "b" for g_kdl_c term 
         do k = 1,n_occ
            do d = 1,n_vir
               do l = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  kdl  = index_three(k,d,l,n_occ,n_vir)
                  ld   = index_two(l,d,n_occ)
                  kc   = index_two(k,c,n_occ)
!
                  cl   = index_two(c,l,n_vir)
                  ck   = index_two(c,k,n_vir)
                  dl   = index_two(d,l,n_vir)
                  dk   = index_two(d,k,n_vir)
                  ckdl = index_packed(ck,dl)
                  cldk = index_packed(cl,dk)
!
!                 Set the values of u_b_kdl and g_kdl_c
!
                  u_b_kdl(c,kdl) = two*t2am(ckdl,1) - t2am(cldk,1)
!
                  g_kdl_c(kdl,c) = g_ld_kc(ld,kc) ! g_ldkc 
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the unordered integrals g_ld_kc = g_ldkc
!
      call deallocator(g_ld_kc,n_ov,n_ov) !
                                          !   It is probably better to simply reorder g_kdl_c to g_k_dlc
                                          !   in the calculation of the E2.2 term (requires less memory).
                                          !
                                          !   For now (8 Mar 2017), I'll just keep it simple & stupid.
                                          !
!
!     Allocate the intermediate X_b_c = F_bc - sum_dkl g_ldkc u_kl^bd and set to zero
!
      call allocator(X_b_c,n_vir,n_vir)
      X_b_c = zero 
!
!     Copy the virtual-virtual Fock matrix into the intermediate 
!
      call dcopy(n_vv,F_a_b,1,X_b_c,1) ! X_b_c = F_bc 
!
!     Add the second contribution, - sum_dkl g_ldkc u_kl^bd = - sum_dkl u_b_kdl * g_kdl_c, to X_b_c
!
      call dgemm('N','N',n_vir,n_vir,n_oov,&
                  -one,u_b_kdl,n_vir,g_kdl_c,n_oov,&
                  one,X_b_c,n_vir)
!
!     Deallocate u_b_kdl and g_kdl_c
!
      call deallocator(u_b_kdl,n_vir,n_oov)
      call deallocator(g_kdl_c,n_oov,n_vir)
!
!     Allocate t_c_jai = t_ij^ac and set to zero
!
      call allocator(t_c_jai,n_vir,n_oov)
      t_c_jai = zero 
!
!     Determine t_c_jai = t_ij^ac 
!
      do c = 1,n_vir
         do j = 1,n_occ
            do a = 1,n_vir
               do i = 1,n_occ
!
!                 Calculate the necessary indices
!
                  jai  = index_three(j,a,i,n_occ,n_vir)
                  ai   = index_two(a,i,n_vir)
                  cj   = index_two(c,j,n_vir)
                  aicj = index_packed(ai,cj)
!
!                 Set the value of t_c_jai 
!
                  t_c_jai(c,jai) = t2am(aicj,1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the E2.1 term and set it to zero
!
      call allocator(omega2_b_jai,n_vir,n_oov)
      omega2_b_jai = zero 
!
!     Calculate the E2.1 term 
!
      call dgemm('N','N',n_vir,n_oov,n_vir,&
                  one,X_b_c,n_vir,t_c_jai,n_vir,&
                  zero,omega2_b_jai,n_vir)
!
!     Add the E2.1 term to the omega vector 
!
      do a = 1,n_vir
         do i = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  jai  = index_three(j,a,i,n_occ,n_vir)
                  ai   = index_two(a,i,n_vir)
                  bj   = index_two(b,j,n_vir)
                  aibj = index_packed(ai,bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
                     omega2(aibj,1) = omega2(aibj,1) + omega2_b_jai(b,jai)
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the E2.1 term, the X intermediate, and the reordered amplitudes 
!
      call deallocator(omega2_b_jai,n_vir,n_oov)
      call deallocator(X_b_c,n_vir,n_vir)
      call deallocator(t_c_jai,n_vir,n_oov)
!
!     Print the omega vector, having added the E2.1 term to it
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after E2.1 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,n_ov_ov_packed)
      endif 
!
!     Allocate the Cholesky vector L_kc_J = L_kc^J and set it to zero 
!
      call allocator(L_kc_J,n_ov,n_J)
      L_kc_J = zero
!
!     Read the Cholesky vector from file 
!
      call read_cholesky_ia(L_kc_J)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_ld_kc,n_ov,n_ov)
      g_ld_kc = zero
!
!     Calculate g_ld_kc = sum_J L_ld^J L_kc^J 
!
      call dgemm('N','T',n_ov,n_ov,n_J,&
                  one,L_kc_J,n_ov,L_kc_J,n_ov,&
                  zero,g_ld_kc,n_ov)
!
!     Deallocate the Cholesky vector L_kc_J
!
      call deallocator(L_kc_J,n_ov,n_J)
!
!     Allocate g_k_dlc = g_ldkc and set to zero 
!
      call allocator(g_k_dlc,n_occ,n_ovv)
      g_k_dlc = zero 
!
!     Allocate u_dlc_j = u_lj^dc and set to zero
!
      call allocator(u_dlc_j,n_ovv,n_occ)
      u_dlc_j = zero 
!
!     Determine g_k_dlc = g_ldkc and u_dlc_j = u_lj^dc 
!
      do k = 1,n_occ ! Use as though "j" for u_dlc_j term 
         do d = 1,n_vir
            do l = 1,n_occ
               do c = 1,n_vir
!
!                 Calculate the necessary indices 
!
                  dlc  = index_three(d,l,c,n_vir,n_occ)
                  ld   = index_two(l,d,n_occ)
                  kc   = index_two(k,c,n_occ)
!
                  dl   = index_two(d,l,n_vir)
                  ck   = index_two(c,k,n_vir)
                  dlck = index_packed(dl,ck)
!
                  dk   = index_two(d,k,n_vir)
                  cl   = index_two(c,l,n_vir)
                  dkcl = index_packed(dk,cl)
!
!                 Set the value of g_k_dlc and u_dlc_j 
!
                  g_k_dlc(k,dlc) = g_ld_kc(ld,kc)                   ! g_ldkc 
                  u_dlc_j(dlc,k) = two*t2am(dlck,1)-t2am(dkcl,1)    ! u_lk^dc = 2 * t_lk^dc - t_kl^dc 
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the integrals g_ld_kc = g_ldkc 
!
      call deallocator(g_ld_kc,n_ov,n_ov)
!
!     Allocate the intermediate Y_k_j = F_kj  + sum_cdl u_lj^dc g_ldkc 
!                                     = F_k_j + sum_cdl g_k_dlc * u_dlc_j, and set it to zero
!
      call allocator(Y_k_j,n_occ,n_occ)
      Y_k_j = zero 
!
!     Copy the occupied-occupied Fock matrix, such that Y_k_j = F_kj 
!
      call dcopy(n_oo,F_i_j,1,Y_k_j,1)
!
!     Add sum_cdl g_k_dlc u_dlc_j to Y_k_j, such that Y_k_j = F_k_j + sum_cdl g_k_dlc u_dlc_j
!
      call dgemm('N','N',n_occ,n_occ,n_ovv,&
                 one,g_k_dlc,n_occ,u_dlc_j,n_ovv,&
                 one,Y_k_j,n_occ)
!
!     Deallocate u_dlc_j and g_k_dlc 
!
      call deallocator(u_dlc_j,n_ovv,n_occ)
      call deallocator(g_k_dlc,n_occ,n_ovv)
!
!     Allocate t_aib_k = t_ik^ab and set it to zero 
!
      call allocator(t_aib_k,n_ovv,n_occ)
      t_aib_k = zero
!
!     Determine t_aib_k = t_ik^ab 
!
      do a = 1,n_vir
         do i = 1,n_occ
            do b = 1,n_vir
               do k = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  aib  = index_three(a,i,b,n_vir,n_occ)
                  ai   = index_two(a,i,n_vir)
                  bk   = index_two(b,k,n_vir)
                  aibk = index_packed(ai,bk)
!
!                 Set the value of t_aib_k 
!
                  t_aib_k(aib,k) = t2am(aibk,1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the E2.2 term and set to zero 
!
      call allocator(omega2_aib_j,n_ovv,n_occ)
      omega2_aib_j = zero
!
!     Calculate the E2.2 term, - sum_k t_aib_k Y_k_j = - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc)
!
      call dgemm('N','N',n_ovv,n_occ,n_occ,&
                  -one,t_aib_k,n_ovv,Y_k_j,n_occ,&
                  zero,omega2_aib_j,n_ovv)
!
!     Deallocate t_aib_k and Y_k_j 
!
      call deallocator(t_aib_k,n_ovv,n_occ)
      call deallocator(Y_k_j,n_occ,n_occ)
!
!     Add the E2.2 term to the omega vector 
!
      do a = 1,n_vir
         do i = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a,i,n_vir)
                  bj   = index_two(b,j,n_vir)
                  aibj = index_packed(ai,bj)
!
                  aib  = index_three(a,i,b,n_vir,n_occ)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
                     omega2(aibj,1) = omega2(aibj,1) + omega2_aib_j(aib,j)
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the E2.2 term 
!
      call deallocator(omega2_aib_j,n_ovv,n_occ)
!
!     Print the omega vector, having added both the E2.1 and E2.2 terms
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after E2.2 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,n_ov_ov_packed)
      endif 
!
   end subroutine mlcc_omega_e2
!
   subroutine mlcc_omega_d2
!
!     MLCC Omega D2 term
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, 8 Mar 2017
!
!     NB! Needs to be rewritten with T1 transformed integrals
!     eventually (makes no difference for MP2 guess)
!
!     This routine calculates the D2 term,
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
!     The first term is referred to as the D2.1 term, and comes out ordered as (....) 
!     The second term is referred to as the D2.2 term, and comes out ordered as (....)
!     The third term is referred to as the D2.3 term, and comes out ordered as (ai,bj)
!
!     All terms are added to the omega vector element omega2(ai,bj)
!     The routine adds the terms in the following order: D2.3, D2.1, D2.2
!
      implicit none 
!
      logical :: debug = .true.
!
      integer :: required,available,max_batch_length,batch_dimension,n_batch
!
      integer :: a_begin,a_end,a_batch,batch_length,a_full,ac_dim 
!
      integer :: ai,aidl,al,aldi,a,i,b,j,c,d,di,dl,k,kc,kd,l,lc,ld,aibj,bj,bjck,bk,bkcj,cj,ck,ca,ki
!
      real(dp), dimension(:,:), pointer :: L_kc_J       => null() ! L_kc^J 
      real(dp), dimension(:,:), pointer :: g_ld_kc      => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: L_ld_kc      => null() ! L_ldkc = 2 * g_ldkc - g_lckd 
      real(dp), dimension(:,:), pointer :: u_ai_ld      => null() ! u_il^ad = 2 * t_il^ad - t_li^ad 
!
      real(dp), dimension(:,:), pointer :: omega2_ai_bj => null() ! For storing the D2.3 & D2.1 terms temporarily
!
      real(dp), dimension(:,:), pointer :: g_ai_kc      => null() ! g_aikc 
      real(dp), dimension(:,:), pointer :: u_kc_bj      => null() ! u_jk^bc
      real(dp), dimension(:,:), pointer :: L_ai_J       => null() ! L_ai^J 
!
      real(dp), dimension(:,:), pointer :: g_ai_ck      => null() ! g_acki
      real(dp), dimension(:,:), pointer :: g_ca_ki      => null() ! g_acki; a is batched over 
      real(dp), dimension(:,:), pointer :: L_ca_J       => null() ! L_ac^J; a is batched over 
      real(dp), dimension(:,:), pointer :: L_ki_J       => null() ! L_ki^J 
      real(dp), dimension(:,:), pointer :: u_ck_bj      => null() ! u_jk^bc
!
!     Allocate the Cholesky vector L_kc_J = L_kc^J and set to zero 
!
      call allocator(L_kc_J,n_ov,n_J)
      L_kc_J = zero
!
!     Read the Cholesky vector L_kc_J from file
!
      call read_cholesky_ia(L_kc_J)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_ld_kc,n_ov,n_ov)
      g_ld_kc = zero 
!
!     Calculate g_ld_kc = g_ldkc = sum_J L_ld^J L_kc^J
!
      call dgemm('N','T',n_ov,n_ov,n_J,&
                  one,L_kc_J,n_ov,L_kc_J,n_ov,&
                  zero,g_ld_kc,n_ov)
!
!     Allocate L_ld_kc = L_ldkc and set to zero    
!
      call allocator(L_ld_kc,n_ov,n_ov)
      L_ld_kc = zero 
!
!     Determine L_ld_kc = L_ldkc from g_ld_kc = g_ldkc 
!
      do l = 1,n_occ
         do d = 1,n_vir
            do k = 1,n_occ
               do c = 1,n_vir
!
!                 Calculate the necessary indices 
!
                  ld = index_two(l,d,n_occ)
                  kc = index_two(k,c,n_occ)
                  lc = index_two(l,c,n_occ)
                  kd = index_two(k,d,n_occ)
!
!                 Set the value of L_ld_kc = 2 * g_ldkc - g_lckd 
!
                  L_ld_kc(ld,kc) = two*g_ld_kc(ld,kc) - g_ld_kc(lc,kd)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ld_kc and L_kc_J
!
      call deallocator(g_ld_kc,n_ov,n_ov)
      call deallocator(L_kc_J,n_ov,n_J)
!
!     Allocate u_ai_ld = u_il^ad and set to zero 
!
      call allocator(u_ai_ld,n_ov,n_ov)
      u_ai_ld = zero 
! 
!     Determine u_ai_ld = u_il^ad = 2 * t_il^ad - t_li^ad 
!
      do a = 1,n_vir
         do i = 1,n_occ
            do l = 1,n_occ
               do d = 1,n_vir
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a,i,n_vir)
                  ld   = index_two(l,d,n_occ)
!
                  dl   = index_two(d,l,n_vir)
                  aidl = index_packed(ai,dl)
!
                  al   = index_two(a,l,n_vir)
                  di   = index_two(d,i,n_vir)
                  aldi = index_packed(al,di)
!
!                 Set the value of u_ai_ld 
!
                  u_ai_ld(ai,ld) = two*t2am(aidl,1) - t2am(aldi,1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the D2.3 term omega2_ai_bj and set it to zero 
!
      call allocator(omega2_ai_bj,n_ov,n_ov)
      omega2_ai_bj = zero 
!
!     Form the temporary vector sum_dl u_ai_ld L_ld_kc and place it in omega2_ai_bj
!
      call dgemm('N','N',n_ov,n_ov,n_ov,&
                  one,u_ai_ld,n_ov,L_ld_kc,n_ov,&
                  zero,omega2_ai_bj,n_ov)          ! Think of the result as an intermediate, say Z_ai_kc    
!
!     Form the D2.3 term, 1/4 sum_kc Z_ai_kc u_kc_bj = 1/4 sum_kc Z_ai_kc(ai,kc) u_ai_ld(bj,kc)
!
      call dgemm('N','T',n_ov,n_ov,n_ov,&
                  one/four,omega2_ai_bj,n_ov,u_ai_ld,n_ov,& 
                  zero,omega2_ai_bj,n_ov)
!
!     Some mathematical justification for the above matrix multiplication. We have 
!
!           1/4 * sum_ck (sum_dl u_il^ad L_ldkc) u_jk^bc = 1/4 * sum_ck Z_ai,kc u_kc,bj,
!
!     where Z_ai,kc = sum_dl u_ai,ld L_ld,kc. Note that u_ai_ld(ai,ld) = u_il^ad, 
!     which means that u_ai_ld(bj,kc)^T = u_ai_ld(kc,bj) = u_kj^cb = u_jk^bc.
!
!
!     Add the D2.3 term to the omega vector 
!
      do a = 1,n_vir
         do i = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a,i,n_vir)
                  bj   = index_two(b,j,n_vir)
                  aibj = index_packed(ai,bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
                     omega2(aibj,1) = omega2(aibj,1) + omega2_ai_bj(ai,bj)
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the omega2_ai_bj and u_ai_ld(ai,ld) = u_il^ad vector
!
      call deallocator(omega2_ai_bj,n_ov,n_ov)
      call deallocator(u_ai_ld,n_ov,n_ov) ! Eirik: note that u_ai_ld(bj,kc) = u_jk^bc, and thus this vector could be reused 
!                                                  to calculate the D2.1 term. This would require some more memory than
!                                                  I have assumed available, for now (9 Mar).
!
!     Print the omega vector, having added the D2.3 term 
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after D2.3 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,n_ov_ov_packed)
      endif                
!
!     Allocate the L_ai_J and L_kc_J terms and set them to zero 
!
      call allocator(L_ai_J,n_ov,n_J)
      call allocator(L_kc_J,n_ov,n_J)
!
      L_ai_J = zero
      L_kc_J = zero
!
!     Read the Cholesky vectors from file 
!
      call read_cholesky_ai(L_ai_J)
      call read_cholesky_ia(L_kc_J)
!
!     Allocate g_ai_kc = g_aikc and set it zero
!
      call allocator(g_ai_kc,n_ov,n_ov)
      g_ai_kc = zero 
!
!     Form the g_ai_kc integrals from L_ai_J and L_kc_J
!  
      call dgemm('N','T',n_ov,n_ov,n_J,&
                  one,L_ai_J,n_ov,L_kc_J,n_ov,&
                  zero,g_ai_kc,n_ov)
!
!     Deallocate the Cholesky vectors L_ai_J and L_kc_J
!
      call deallocator(L_ai_J,n_ov,n_J)
      call deallocator(L_kc_J,n_ov,n_J)
!
!     Allocate u_kc_bj and set it to zero 
!
      call allocator(u_kc_bj,n_ov,n_ov)
      u_kc_bj = zero
!
!     Determine u_kc_bj = u_jk^bc = 2 * t_jk^bc - t_kj^bc 
!
      do k = 1,n_occ
         do c = 1,n_vir
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  kc   = index_two(k,c,n_occ)
                  bj   = index_two(b,j,n_vir)
!  
                  ck   = index_two(c,k,n_vir)
                  bjck = index_packed(bj,ck)
!
                  bk   = index_two(b,k,n_vir)
                  cj   = index_two(c,j,n_vir)
                  bkcj = index_packed(bk,cj)
!
!                 Set the value of u_kc_bj
!     
                  u_kc_bj(kc,bj) = two*t2am(bjck,1) - t2am(bkcj,1)
!
               enddo
            enddo
         enddo
      enddo 
!
!     Allocate omega2_ai_bj and set it to zero 
!
      call allocator(omega2_ai_bj,n_ov,n_ov)
      omega2_ai_bj = zero 
!
!     Calculate the D2.1 term sum_ck u_jk^bc g_aikc = sum_ck g_ai_kc u_kc_bj
!
      call dgemm('N','N',n_ov,n_ov,n_ov,&
                  one,g_ai_kc,n_ov,u_kc_bj,n_ov,&
                  zero,omega2_ai_bj,n_ov)
!
!     Add the D2.1 term to the omega vector 
!
      do a = 1,n_vir
         do i = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a,i,n_vir)
                  bj   = index_two(b,j,n_vir)
                  aibj = index_packed(ai,bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
                     omega2(aibj,1) = omega2(aibj,1) + omega2_ai_bj(ai,bj)
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ai_kc, u_kc_bj, and the omega2_ai_bj vectors 
!
      call deallocator(g_ai_kc,n_ov,n_ov)
      call deallocator(u_kc_bj,n_ov,n_ov)
      call deallocator(omega2_ai_bj,n_ov,n_ov)
!
!     Print the omega vector, having added the D2.1 term 
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after D2.1 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,n_ov_ov_packed)
      endif 
!
! - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj 
!
!     Allocate L_ki_J and set it to zero
!
      call allocator(L_ki_J,n_oo,n_J)
      L_ki_J = zero 
!
!     Read the Cholesky vector L_ki_J from file 
!
      call read_cholesky_ij(L_ki_J)
!
!     Allocate the full g_ai_ck = g_acki and set it to zero 
!
      call allocator(g_ai_ck,n_ov,n_ov)
      g_ai_ck = zero 
!
!     Prepare for batching over the index a to calculate g_ai_ck = g_acki
!
!        To calculate this term, we need to hold L_ac^J and g_acki
!        in memory simultaneously 
!
      required = n_vir*n_vir*n_J + n_vir*n_vir*n_occ*n_occ ! Eirik: This estimate has to be updated,
                                                           !        when we account for the T1 transformation
      available = get_available()
      batch_dimension = n_vir
!
!     Determine the batching variables 
!
      call n_one_batch(required,available,max_batch_length,n_batch,batch_dimension) 
!
!     Determine g_ai_ck = g_acki successively in batches over a 
!
      do a_batch = 1,n_batch
!
!        For each batch, get the limits for the a index        
!
         call one_batch_limits(a_begin,a_end,a_batch,max_batch_length,batch_dimension)
         batch_length = a_end - a_begin + 1 
!
!        Allocate the Cholesky vector L_ca_J = L_ac^J and set it to zero 
!
         ac_dim = batch_length*n_vir         ! Dimension of ac for the batch over index a 
         call allocator(L_ca_J,ac_dim,n_J)
         L_ca_J = zero
!
!        Read the Cholesky vector from file 
!
         call read_cholesky_ab_reorder(L_ca_J,a_begin,a_end,ac_dim)
!
!        Allocate the integral g_ca_ki = g_acki and set to zero 
!
         call allocator(g_ca_ki,ac_dim,n_ov)
         g_ca_ki = zero
!
!        Calculate g_ca_ki = g_acki from L_ca_J = L_ac^J and L_ki_J = L_ki^J
!
         call dgemm('N','T',ac_dim,n_oo,n_J,&
                     one,L_ca_J,ac_dim,L_ki_J,n_oo,&
                     one,g_ca_ki,ac_dim)
!
!        Reorder the integrals g_ca_ki (reduced a) = g_acki = g_ai_ck (full a)
!
         do a = 1,batch_length
            do c = 1,n_vir
               do k = 1,n_occ
                  do i = 1,n_occ
!
!                    Calculate the necessary indices 
!
                     a_full = a - 1 + a_begin            ! The full matrix index a
                     ai     = index_two(a_full,i,n_vir)
                     ck     = index_two(c,k,n_vir)
!
                     ca     = index_two(c,a,n_vir)
                     ki     = index_two(k,i,n_occ)
!
!                    Set the value of g_ai_ck = g_acki 
!
                     g_ai_ck(ai,ck) = g_ca_ki(ca,ki) ! g_acki
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate the g_ca_ki and L_ca_J vectors
!
         call deallocator(g_ca_ki,ac_dim,n_oo)
         call deallocator(L_ca_J,ac_dim,n_J)
!
      enddo ! End of loop over batches of a 
!
! ! - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj 
!
!     Allocate the u_ck_bj = u_jk^bc vector and set it to zero 
!
      call allocator(u_ck_bj,n_ov,n_ov)
      u_ck_bj = zero 
!
!     Determine u_ck_bj = u_jk^bc = 2 * t_jk^bc - t_kj^bc 
!
      do c = 1,n_vir
         do k = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  ck   = index_two(c,k,n_vir)
                  bj   = index_two(b,j,n_vir)
!
                  bjck = index_packed(bj,ck)
!
                  bk   = index_two(b,k,n_vir)
                  cj   = index_two(c,j,n_vir)
                  bkcj = index_packed(bk,cj)
!
!                 Set the value of u_ck_bj 
!
                  u_ck_bj(ck,bj) = two*t2am(bjck,1) - t2am(bkcj,1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the D2.2 term and set it to zero 
!
      call allocator(omega2_ai_bj,n_ov,n_ov)
      omega2_ai_bj = zero
!
!     Calculate the D2.2 term, - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj
!
      call dgemm('N','N',n_ov,n_ov,n_ov,&
                  -one/two,g_ai_ck,n_ov,u_ck_bj,n_ov,&
                  zero,omega2_ai_bj,n_ov)
!
!     Add the D2.2 term to the omega vector 
!
      do a = 1,n_vir
         do i = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a,i,n_vir)
                  bj   = index_two(b,j,n_vir)
                  aibj = index_packed(ai,bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
                     omega2(aibj,1) = omega2(aibj,1) + omega2_ai_bj(ai,bj)
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
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after D2.2 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,n_ov_ov_packed)
      endif 
!
!     Deallocate g_ai_ck, u_ck_bj, and the temporary omega2_ai_bj
!
      call deallocator(g_ai_ck,n_ov,n_ov)
      call deallocator(u_ck_bj,n_ov,n_ov)
      call deallocator(omega2_ai_bj,n_ov,n_ov)
!
   end subroutine mlcc_omega_d2
!
   subroutine mlcc_omega_c2
!
!     MLCC Omega C2 term
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 8 Mar 2017
!
!     Omega C2 = -1/2* sum_(ck)t_bk_cj*(g_ki_ac -1/2 sum_(dl)t_al_di * g_kd_lc)
!                 - sum_(ck) t_bk_ci (g_kj_ac-sum_(dl)t_al_dj*g_kd_lc)
      implicit none
!
      real(dp),dimension(:,:),pointer        :: L_ia_J => null()
      real(dp),dimension(:,:),pointer        :: L_ki_J => null()
      real(dp),dimension(:,:),pointer        :: L_ca_J => null()
      real(dp),dimension(:,:),pointer        :: g_kd_lc => null()
      real(dp),dimension(:,:),pointer        :: g_dl_ck => null()
      real(dp),dimension(:,:),pointer        :: g_ki_ca => null()
      real(dp),dimension(:,:),pointer        :: g_ai_ck => null()
      real(dp),dimension(:,:),pointer        :: t_ai_dl => null()
      real(dp),dimension(:,:),pointer        :: t_ck_bj => null()
      real(dp),dimension(:,:),pointer        :: X_ai_ck => null()
      real(dp),dimension(:,:),pointer        :: Y_ai_bj => null()
      integer                                :: c,k,d,l,kd,lc,ck,dl,al,di,ai,aldi
      integer                                :: i,a,ca,ki,b,bj,bk,cj,j,bkcj,aibj,aj,bi
      integer                                :: required,available,n_batch,max_batch_length,a_start
      integer                                :: a_end,a_batch,a_length
      logical                                :: debug = .true.
!
!     Allocate L_ia_J
!
      call allocator(L_ia_J,n_ov,n_J)
      L_ia_J=zero
!
!     Get L_ia_J
!
      call get_cholesky_ia(L_ia_J)
!
!     Create g_kd_lc = sum_J L_kd_J * L_lc_J
!
      call allocator(g_kd_lc,n_ov,n_ov)
      g_kd_lc=zero
!
      call dgemm('N','T',n_ov,n_ov,n_J &
         ,one,L_ia_J,n_ov,L_ia_J,n_ov &
         ,zero,g_kd_lc,n_ov)
!
!     Reorder g_kd_lc as g_dl_ck and t_al_di as t_ai_dl
!
      call allocator(g_dl_ck,n_ov,n_ov)
      call allocator(t_ai_dl,n_ov,n_ov)
      g_dl_ck=zero
      t_ai_dl=zero
!
      do c = 1,n_vir
         do d = 1,n_vir
            do k = 1,n_occ
               do l = 1,n_occ
!
!                 Needed indices for reordering of g
!
                  kd=index_two(k,d,n_occ)
                  lc=index_two(l,c,n_occ)
                  dl=index_two(d,l,n_vir)
                  ck=index_two(c,k,n_vir)
!
!                 Needed indices for reordering of t
!
                  al=index_two(c,l,n_vir)
                  di=index_two(d,k,n_vir)
                  ai=ck
                  aldi=index_packed(al,di)
!
                  g_dl_ck(dl,ck)=g_kd_lc(kd,lc)
!
                  t_ai_dl(ai,dl)=t2am(aldi,1)
               enddo
            enddo
         enddo
      enddo
!
!     -1/2 * sum_(dl) t_ai_dl*g_dl_ck = X_ai_ck
!
      call allocator(X_ai_ck,n_ov,n_ov)
!
      call dgemm('N','N',n_ov,n_ov,n_ov &
         ,-half,t_ai_dl,n_ov,g_dl_ck,n_ov &
         ,zero,X_ai_ck,n_ov)
!
!     Deallocate L_ia_J, g_kd_lc and g_dl_ck
!
      call deallocator(L_ia_J,n_ov,n_J)
      call deallocator(g_kd_lc,n_ov,n_ov)
      call deallocator(g_dl_ck,n_ov,n_ov)
!
!     Constructing g_ki_ac ordered as g_ki_ca
!
!
!     Allocate g_ki_ca
!
      call allocator(g_ki_ca,n_oo,n_vv)
      g_ki_ca = zero
!
!     Allocate L_ki_J
!
      call allocator(L_ki_J,n_ov,n_J)
!
!     Get cholesky vectors of ij-type
!
      call get_cholesky_ij(L_ki_J)
!
!     Prepare batching over a 
!
!
!     Setup of variables needed for batching
!
      available = get_available()
      required = 2*n_vir*n_vir*n_J*4 + 2*n_vir*n_occ*n_J*4
      call n_one_batch(required,available,max_batch_length,n_batch,n_vir)
!
      a_start=1
      a_end=0
      a_length=0
!
!     Start looping over batches
!
      do a_batch = 1,n_batch
!
!        Get batch limits  and  length of batch
!
         call one_batch_limits(a_start,a_end,a_batch,max_batch_length,n_vir)
         a_length=a_end-a_start+1
!
!        Allocation for L_ac_J as L_ca_J (L_ca_J = L_acJ)
!
         call allocator(L_ca_J,n_vir*a_length,n_J)
         L_ca_J=zero
!
!        Read Cholesky vectors
!
         call get_cholesky_ab(L_ca_J,a_start,a_end,n_vir*a_length,.true.)
!
!        g_ki_ca = sum_J L_ki_J*L_ca_J
!
         call dgemm('N','T',n_oo,n_vir*a_length,n_J &
            ,one,L_ki_J,n_oo,L_ca_J,n_vir*a_length &
            ,one,g_ki_ca(1,index_two(1,a_start,n_vir)),n_oo)
!
!        Deallocate L_ca_J
!
         call deallocator(L_ca_J,n_vir*a_length,n_J)
      enddo ! End of batching
!
!     Deallocate L_ki_J
!
      call deallocator(L_ki_J,n_ov,n_J)
!
!     Reorder g_ki_ca to g_ai_ck
!
      call allocator(g_ai_ck,n_ov,n_ov)
!
      do i = 1,n_occ
         do k = 1,n_occ
            do a = 1,n_vir
               do c = 1,n_vir
!
!                 Needed indices
!
                  ki=index_two(k,i,n_occ)
                  ca=index_two(c,a,n_vir)
                  ai=index_two(a,i,n_vir)
                  ck=index_two(c,k,n_vir)
!
                  g_ai_ck(ai,ck) = g_ki_ca(ki,ca)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ki_ca,n_oo,n_vv)
!
!     X_ai_ck = X_ai_ck + g_ai_ck
!
      call daxpy(n_ov*n_ov,one,g_ai_ck,1,X_ai_ck,1)
!
!     Deallocate g_ai_kc
!
      call deallocator(g_ai_ck,n_ov,n_ov)
!
!     reorder t_bkcj_1 as t_ck_bj
!
      call allocator(t_ck_bj,n_ov,n_ov)
!
      do c = 1,n_vir
         do k = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Needed indices
!
                  bk = index_two(b,k,n_vir)
                  cj = index_two(c,j,n_vir)
                  bkcj = index_packed(bk,ck)
!
                  bj = index_two(b,j,n_vir)
                  ck = index_two(c,k,n_vir)
!
                  t_ck_bj = t1am(bkcj,1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate intermediate Y_ai_bj
!
      call allocator(Y_ai_bj,n_ov,n_ov)
      Y_ai_bj = zero
!
!     Y_ai_bj = - sum_(ck) X_ai_ck*t_ck_bj
!
      call dgemm('N','N',n_ov,n_ov,n_ov &
         ,-one,X_ai_ck,n_ov,t_ck_bj,n_ov &
         , zero,Y_ai_bj,n_ov)
!
!     Deallocate t_ck_bj
!
      call deallocator(t_ck_bj,n_ov,n_ov)
!
!     Omega_aibj,1 = 1/2*Y_ai_bj + Y_aj_bi
!
      do a = 1,n_vir
         do i = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Needed indices
!
                  ai=index_two(a,i,n_vir)
                  bj=index_two(b,j,n_vir)
!
                  if (ai .ge. bj) then
                     aj=index_two(a,j,n_vir)
                     bi=index_two(b,i,n_vir)
!
                     aibj=index_packed(ai,bj)
!
                     omega2(aibj,1)=omega2(aibj,1)+half*Y_ai_bj(ai,bj)+Y_ai_bj(aj,bi)
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
      call deallocator(Y_ai_bj,n_ov,n_ov)
!
!
!     Print the omega vector, having added D2
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after D2 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,n_ov_ov_packed)
      endif 
!
   end subroutine mlcc_omega_c2
!
   subroutine mlcc_omega_a2
   end subroutine mlcc_omega_a2
!
   subroutine mlcc_omega_b2
   end subroutine mlcc_omega_b2
!
end module mlcc_omega