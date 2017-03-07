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
      if (debug) required = 1500000
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
!  (Omega_ai^C1=sum_ck F_ck_T1 u^ac_ik)
!
!  Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mars 2017
!
   use mlcc_data ! Eirik: These use statements are not necessary (included in the module already)
   use mlcc_utilities
   use mlcc_workspace
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
   call allocator(F_ck,1,n_ov) ! T1-transformed fock matrix
!
! Set up u_ck_ai and MO Fock matrix
!
   do k = 1,n_occ
      do c = 1,n_vir
         ck = index_two(c,k,n_vir)
!
!        MO Fock matrix
!
         F_ck(1,ck) = mo_fock_mat(n_occ+c,k)
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
!  T1 transformation Sarai: Should have T1-transformed Fock in mem
!
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
  ! real(dp),dimension(:,:),pointer  :: F_a_i => null() ! Eirik: I commented this 
!
!  Allocation
!
   call allocator(F_a_i,n_vir,n_occ)
!
!  MO Fock matrix
!
   ! do i = 1,n_occ ! Eirik: I commented this 
   !       do a = 1,n_vir
   !       F_a_i(a,i) = mo_fock_mat(n_occ+a,i)
   !    enddo
   ! enddo
!
!  T1 transformation
!
!
!  Add to omega
!
   call daxpy(n_ov,one,F_a_i,1,omega1,1)
!
!  Deallocation
!
   call deallocator(F_a_i,n_vir,n_occ)
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
!     Calculates sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd) - sum_k t_ik^ab (F_kj + sum_cdl g_kdlc u_lj^cd)
!
!     The first term is referred to as the E2.1 term, and comes out ordered as (b,jai) 
!     The second term is referred to as the E2.2 term, and comes out ordered as (aib,j)
!
!     Both are added to the omega vector element omega2(ai,bj)
!
      implicit none 
!
      integer :: b,c,k,d,ck,ckdl,cl,cldk,dk,dl,kc,kdl,l,ld
!
      real(dp), dimension(:,:), pointer :: omega2_b_jai => null() ! To store the E2.1 term temporarily
      real(dp), dimension(:,:), pointer :: omega2_aib_j => null() ! To store the E2.2 term temporarily
      real(dp), dimension(:,:), pointer :: L_kc_J => null() ! L_kc^J
      real(dp), dimension(:,:), pointer :: g_ld_kc => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: g_kdl_c => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: u_b_kdl => null() ! u_kl^bd 
      real(dp), dimension(:,:), pointer :: F_b_c => null() ! F_bc, the virtual-virtual Fock matrix
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
!     Allocate g_ck_dl = g_ckdl and set to zero 
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
      call deallocator(g_ld_kc,n_ov,n_ov)
!
!     Allocate F_b_c = F_bc and set to zero 
!
      call allocator(F_b_c,n_vir,n_vir)
      F_b_c = zero 
!
!     Set the virtual-virtual Fock matrix F_b_c
!
      do b = 1,n_vir
         do c = 1,n_vir
            !
         enddo
      enddo
!
   end subroutine mlcc_omega_e2
!
   subroutine mlcc_omega_d2
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