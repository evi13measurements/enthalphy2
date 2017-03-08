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
!     I. The singles contribution to < mu | exp(-T) H exp(T) | R >
!
      call mlcc_omega_a1
      call mlcc_omega_b1
      call mlcc_omega_c1
      call mlcc_omega_d1 
!
!     II. The doubles contribution to < mu | exp(-T) H exp(T) | R >
!
      write(luprint,*) 'Calculating the doubles omega vector...'
      call mlcc_omega_e2
     ! call mlcc_omega_d2 ! Eirik: I am working on this - commenting to avoid errors...
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
      integer :: lucho_ia,lucho_ab
!
      integer :: i,j,a,idummy
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
!        Calculate u_ckd_i
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
!
!        Allocate the Cholesky vector L_da_J = L_ad^J
!
         batch_length = a_end - a_begin + 1 
         call allocator(L_da_J,batch_length*n_vir,n_J)
!
!        Read in the reordered Cholesky vector L_da_J = L_ad^J 
!
         ad_dim = batch_length*n_vir ! Dimension of ad for the batch over index a 
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
      call deallocator(u_ckd_i,n_ovv,n_occ)
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
   call deallocator(F_ck,n_ov,1)
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
   use mlcc_data
   use mlcc_utilities
   use mlcc_workspace
!
   implicit none
   integer :: a,i, ai
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
!     Calculates sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd) - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc)
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
      integer :: mem_left
!
      real(dp), dimension(:,:), pointer :: omega2_b_jai => null() ! To store the E2.1 term temporarily
      real(dp), dimension(:,:), pointer :: L_kc_J => null() ! L_kc^J
      real(dp), dimension(:,:), pointer :: g_ld_kc => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: g_kdl_c => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: u_b_kdl => null() ! u_kl^bd 
      real(dp), dimension(:,:), pointer :: F_b_c => null() ! F_bc, the virtual-virtual Fock matrix
      real(dp), dimension(:,:), pointer :: X_b_c => null() ! An intermediate, see below for definition
      real(dp), dimension(:,:), pointer :: t_c_jai => null() ! t_ij^ac 
!
      real(dp), dimension(:,:), pointer :: g_k_dlc => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: u_dlc_j => null() ! u_lj^dc 
      real(dp), dimension(:,:), pointer :: omega2_aib_j => null() ! To store the E2.2 term temporarily
      real(dp), dimension(:,:), pointer :: Y_k_j => null() ! An intermediate, see below for definition 
      real(dp), dimension(:,:), pointer :: t_aib_k => null() ! t_ik^ab 
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
!     Calculate u_b_kdl = u_kl^bd and g_kdl_c = g_ldkc
!
      do b = 1,n_vir ! Use as though "c" for g_kdl_c term 
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
      call deallocator(g_ld_kc,n_ov,n_ov) ! Eirik: for optimization, it may be possible to keep these integrals in memory,
                                          !        because they are also needed for the E2.2 term. 
                                          !
                                          !        It is probably even better to simply reorder g_kdl_c to g_k_dlc
                                          !        in the calculation of the E2.2 term (requires less memory).
                                          !
                                          !        For now (8 Mar 2017), I'll just keep it simple & stupid.
                                          !
!
!     Have the pointer F_b_c point to existing F_a_b (using the former for convenience)
!
      F_b_c => F_a_b
!
!     Allocate the intermediate X_b_c = F_bc - sum_dkl g_ldkc u_kl^bd and set to zero
!
      call allocator(X_b_c,n_vir,n_vir)
      X_b_c = zero 
!
!     Copy the virtual-virtual Fock matrix into the intermediate 
!
      call dcopy(n_vv,F_b_c,1,X_b_c,1) ! X_b_c = F_bc 
!
!     Add the second contribution, - sum_dkl g_ldkc u_kl^bd = - sum_dkl u_b_kdl * g_kdl_c, to X_b_c
!
      call dgemm('N','N',n_vir,n_vir,n_oov,&
                  -one,u_b_kdl,n_vir,g_kdl_c,n_oov,&
                  one,X_b_c,n_vir)
!
!     Deallocate u_b_kdl & g_kdl_c
!
      call deallocator(u_b_kdl,n_vir,n_oov)
      call deallocator(g_kdl_c,n_oov,n_vir)
!
!     Allocate t_c_jai = t_ij^ac and set to zero
!
      call allocator(t_c_jai,n_vir,n_oov)
      t_c_jai = zero 
!
!     Set the value of t_c_jai = t_ij^ac 
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
!     Allocate the E2.1 term and set to zero
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
!                 Restrict the indices to avoid adding (ai,bj) and (bj,ai), as they
!                 are identical in packed indices
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
!     Set the value of g_k_dlc = g_ldkc and u_dlc_j = u_lj^dc 
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
                  g_k_dlc(k,dlc) = g_ld_kc(ld,kc)                  ! g_ldkc 
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
!    Allocate the intermediate Y_k_j = F_kj + sum_cdl u_lj^dc g_ldkc = F_k_j + sum_cdl g_k_dlc * u_dlc_j and set to zero 
!
      call allocator(Y_k_j,n_occ,n_occ)
      Y_k_j = zero 
!
!     Copy the occupied-occupied Fock matrix, such that Y_k_j = F_kj 
!
      call dcopy(n_oo,F_i_j,1,Y_k_j,1)
!
!     Add sum_cdl g_k_dlc * u_dlc_j to Y_k_j, such that Y_k_j = F_k_j + sum_cdl g_k_dlc * u_dlc_j
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
!     Allocate t_aib_k = t_ik^ab and set to zero 
!
      call allocator(t_aib_k,n_ovv,n_occ)
      t_aib_k = zero
!
!     Set the value of t_aib_k = t_ik^ab 
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
!                 Restrict the indices to avoid adding (ai,bj) and (bj,ai), as they
!                 are identical in packed indices
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
!     Calculates sum_ck u_jk^bc g_aikc - 1/2 * sum_ck u_jk^bc g_acki + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad 
!
!     The first term is referred to as the D2.1 term, and comes out ordered as (....) 
!     The second term is referred to as the D2.2 term, and comes out ordered as (....)
!     The third term is referred to as the D2.3 term, and comes out ordered as (ai,bj)
!
!     All terms are added to the omega vector element omega2(ai,bj)
!
!     The routine adds the terms in the following order: D2.3, D2.1, D2.2
!
      implicit none 
!
      integer :: ai,aidl,al,aldi,a,i,b,j,c,d,di,dl,k,kc,kd,l,lc,ld
!
      real(dp), dimension(:,:), pointer :: L_kc_J => null() ! L_kc^J 
      real(dp), dimension(:,:), pointer :: g_ld_kc => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: L_ld_kc => null() ! L_ldkc = 2 * g_ldkc - g_lckd 
      real(dp), dimension(:,:), pointer :: u_ai_ld => null() ! u_il^ad = 2 * t_il^ad - t_li^ad 
      real(dp), dimension(:,:), pointer :: omega2_ai_bj => null() ! To store the D2.3 term temporarily
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
!     Set the value of L_ld_kc using g_ld_kc 
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
                  lc = index_two(l,d,n_occ)
                  kd = index_two(k,d,n_occ)
!
!                 Set the value of L_ld_kc
!
                  L_ld_kc(ld,kc) = two*g_ld_kc(ld,kc) - g_ld_kc(lc,kd)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ld_kc & L_kc_J ! Eirik: we should consider keeping the latter in memory if possible,
!                                   as it must be used later to form the integrals g_aikc 
!
      call deallocator(g_ld_kc,n_ov,n_ov)
      call deallocator(L_kc_J,n_ov,n_J)
!
!     Allocate u_ai_ld = u_il^ad and set to zero 
!
      call allocator(u_ai_ld,n_ov,n_ov)
      u_ai_ld = zero 
! 
!     Set the value of u_ai_ld = u_il^ad = 2 * t_il^ad - t_li^ad 
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
!     Allocate the D2.3 term omega2_ai_bj and set to zero 
!
      call allocator(omega2_ai_bj,n_ov,n_ov)
      omega2_ai_bj = zero 
!
!     Form the temporary vector sum_dl u_ai_ld L_ld_kc and place it in omega2_ai_bj ("bj = ck")
!
      call dgemm('N','N',n_ov,n_ov,n_ov,&
                  one,u_ai_ld,n_ov,L_ld_kc,n_ov,&
                  zero,omega2_ai_bj,n_ov)
!
!     Form the D2.3 term by dgemm (to do!)
!
!..... Remember: transpose u!
!
   end subroutine mlcc_omega_d2
!
   subroutine mlcc_omega_c2
   end subroutine mlcc_omega_c2
!
   subroutine mlcc_omega_a2
   end subroutine mlcc_omega_a2
!
   subroutine mlcc_omega_b2
   end subroutine mlcc_omega_b2
!
end module mlcc_omega