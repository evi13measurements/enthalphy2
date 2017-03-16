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
      integer :: memory_lef = 0
!
!     Add the singles contributions to < mu | exp(-T) H exp(T) | R >
!
         memory_lef = get_available()
         write(luprint,*) 'Memory 1:',memory_lef
         call flshfo(luprint)
     write(luprint,*) 'Entering omega1 A'
      call flshfo(luprint)
      call mlcc_omega_a1 ! This is the G term in the old code. -1.4259234603561466E-003 (ours) -1.4259246117271341E-003 (old)
   !                                                            Error : 1.15 * 10^(-9)
               memory_lef = get_available()
         write(luprint,*) 'Memory 2:',memory_lef
         call flshfo(luprint)
           write(luprint,*) 'Entering omega1 B'
      call flshfo(luprint)
      call mlcc_omega_b1 ! This is the H term in the old code. 2.0521348553625856E-003 (ours) 2.0521147588601551E-003 (old)
   !                                                           Error: 2.00 * 10^(-8) 
               memory_lef = get_available()
         write(luprint,*) 'Memory 3:',memory_lef
         call flshfo(luprint)
           write(luprint,*) 'Entering omega1 C'
      call flshfo(luprint)
      call mlcc_omega_c1
               memory_lef = get_available()
         write(luprint,*) 'Memory 4:',memory_lef
         call flshfo(luprint)
           write(luprint,*) 'Entering omega1 D'
      call flshfo(luprint)
      call mlcc_omega_d1 
               memory_lef = get_available()
         write(luprint,*) 'Memory 5:',memory_lef
         call flshfo(luprint)
           write(luprint,*) 'Exiting omega1'
      call flshfo(luprint)
!
      write(luprint,*) 'Omega(a,i):'
      call vec_print(omega1,n_vir,n_occ)   
!
!     Add the doubles contributions to < mu | exp(-T) H exp(T) | R >
!
      write(luprint,*) 'Entering Omega2 E'
      call flshfo(luprint)
      call mlcc_omega_e2
               memory_lef = get_available()
         write(luprint,*) 'Memory 6:',memory_lef
         call flshfo(luprint)
      write(luprint,*) 'Entering Omega2 D'
      call flshfo(luprint)
      call mlcc_omega_d2 ! Something is not deallocated here!! This should now be fixed...!
               memory_lef = get_available()
         write(luprint,*) 'Memory 7:',memory_lef
         call flshfo(luprint)
      write(luprint,*) 'Entering Omega2 C'
      call flshfo(luprint)
      call mlcc_omega_c2  ! Something is not deallocated here!! This should now be fixed...! NB this is very memory intensive, as it is now 
               memory_lef = get_available()
         write(luprint,*) 'Memory 8:',memory_lef
         call flshfo(luprint)
      write(luprint,*) 'Entering Omega2 A'
      call flshfo(luprint)
      call mlcc_omega_a2
               memory_lef = get_available()
      write(luprint,*) 'Entering Omega2 B'
      call flshfo(luprint)
      write(luprint,*) 'heiiheii5'
      call flshfo(luprint)
      call mlcc_omega_b2
               memory_lef = get_available()
         write(luprint,*) 'Memory 10:',memory_lef
         call flshfo(luprint)
            write(luprint,*) 'Exiting Omega2'
      call flshfo(luprint)
!
      write(luprint,*) 'Omega(aibj,1):'
      call vec_print_packed(omega2,n_ov_ov_packed)
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
      integer :: i=0,j=0,a=0
!
      integer :: required=0,available=0,max_batch_length=0,batch_dimension=0,n_batch=0
!
      integer :: a_begin=0,a_end=0,a_batch=0,batch_length=0
!
      integer :: ad=0,ad_dim=0,c=0,ci=0,cidk=0,ck=0,ckd=0,ckdi=0,di=0,dk=0,k=0,kc=0,d=0,da=0
!
      logical :: debug = .false.
!
      real(dp), dimension(:,:), pointer :: L_kc_J  => null()
      real(dp), dimension(:,:), pointer :: L_da_J  => null()   ! L_ad^J; a is being batched over
      real(dp), dimension(:,:), pointer :: g_da_kc => null()   ! g_adkc; a is being batched over
      real(dp), dimension(:,:), pointer :: g_a_ckd => null()   ! reordered g_adkc; a is being batched over
      real(dp), dimension(:,:), pointer :: u_ckd_i => null()   ! u_ki^cd
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
      call get_cholesky_ia(L_kc_J)
!
!     Allocate u_ckd_i = u_ki^cd
!
      call allocator(u_ckd_i,n_ovv,n_occ)
      u_ckd_i = zero
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
         call get_cholesky_ab(L_da_J,a_begin,a_end,ad_dim,.true.)
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
      logical :: debug = .false.
!
      integer :: a=0,c=0,k=0,l=0,ckl=0,ki=0,ak=0,akcl=0,al=0,alck=0,ck=0,ai=0,cl=0,lc=0,i=0,j=0
!
      real(dp), dimension(:,:), pointer :: L_ki_J  => null() ! L_ki^J 
      real(dp), dimension(:,:), pointer :: L_lc_J  => null() ! L_lc^J 
      real(dp), dimension(:,:), pointer :: g_ki_lc => null() ! g_kilc 
      real(dp), dimension(:,:), pointer :: g_ckl_i => null() ! g_kilc 
      real(dp), dimension(:,:), pointer :: u_a_ckl => null() ! u_kl^ac = 2 t_kl^ac - t_lk^ac
!
!     Allocate Cholesky vectors L_ki,J and L_lc,J 
!
      call allocator(L_ki_J,n_oo,n_J)
      call allocator(L_lc_J,n_ov,n_J)
!
      L_ki_J = zero
      L_lc_J = zero
!
!     Read the Cholesky vectors L_ki_J and L_lc_J
!
      call get_cholesky_ij(L_ki_J)
      call get_cholesky_ia(L_lc_J)
!
!     Allocate integrals g_ki_lc = g_kilc
!
      call allocator(g_ki_lc,n_oo,n_ov)
!
      g_ki_lc = zero
!
!     Calculate g_ki_lc = sum_J L_ki_J L_lc_J^T 
! 
      call dgemm('N','T',n_oo,n_ov,n_J,&
                  one,L_ki_J,n_oo,L_lc_J,n_ov,&
                  zero,g_ki_lc,n_oo) 
!
!     Deallocate the Cholesky vectors L_ki_J and L_lc_J
!
      call deallocator(L_ki_J,n_oo,n_J)
      call deallocator(L_lc_J,n_ov,n_J)
!
!     Allocate reordered integrals g_ckl_i = g_kilc 
!
      call allocator(g_ckl_i,n_oov,n_occ)
!
!     Determine g_ckl_i = g_kilc 
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
!     Allocate redordered u_a_ckl = u_kl^ac and set it to zero
!
      call allocator(u_a_ckl,n_vir,n_oov)
      u_a_ckl = zero
!
!     Determine u_a_ckl = u_kl^ac
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
!                 Set the value of u_a_ckl = u_kl^ac = 2*t_kl^ac - t_lk^ac = 2*t_ak,cl - t_al,ck 
!
                  u_a_ckl(a,ckl) = two*t2am(akcl,1) - t2am(alck,1)
!                  
               enddo
            enddo
         enddo
      enddo
!
!     Calculate the B1 term, - sum_ckl u_a_ckl g_ckl_i
!
      call dgemm('N','N',n_vir,n_occ,n_oov,&
                  -one,u_a_ckl,n_vir,g_ckl_i,n_oov,&
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
   real(dp),dimension(:,:),pointer  :: F_ck      => null()
   real(dp),dimension(:,:),pointer  :: u_ai_ck   => null()
   real(dp),dimension(:,:),pointer  :: omega1_ai => null()
   integer                          :: i=0,k=0,c=0,a=0
   integer                          :: ck=0,ai=0,ak=0,ci=0,aick=0,akci=0
   logical                          :: debug = .false.
!
!
!  Allocation
!
   call allocator(u_ai_ck,n_ov,n_ov)   ! 2*t_ck_ai-t_ci_ak
   call allocator(F_ck,n_ov,1) ! T1-transformed F_k_c reordered
   call allocator(omega1_ai,n_ov,1)
   u_ai_ck=zero
   F_ck=zero
   omega1_ai=zero
!
! Set up u_ck_ai and MO Fock matrix virtual-occupied
!
   do k = 1,n_occ
      do c = 1,n_vir
         ck = index_two(c,k,n_vir)
!
!        MO Fock matrix
!
         F_ck(ck,1) = F_i_a(k,c)
         do i = 1,n_occ
            do a = 1,n_vir
!
!              Necessary indices
!
               ai = index_two(a,i,n_vir)
               ci = index_two(c,i,n_vir)
               ak = index_two(a,k,n_vir)
!
               aick = index_packed(ck,ai)
               akci = index_packed(ci,ak)
!              
!              u_ck_ai
!
               u_ai_ck(ai,ck) = two*t2am(aick,1)-t2am(akci,1)
            enddo
         enddo
      enddo
   enddo
!
!  Matrix multiplication
!
      call dgemm('N','N',n_ov,1,n_ov &
      ,one,u_ai_ck,n_ov,F_ck,n_ov &
      ,zero,omega1_ai,n_ov)
!
!     Add to omega1
!
       do i = 1,n_occ
            do a = 1,n_vir
               ai=index_two(a,i,n_vir)
               omega1(a,i)=omega1(a,i)+omega1_ai(ai,1)
            enddo
         enddo
!
!  Deallocation
!
   call deallocator(F_ck,n_ov,1)
   call deallocator(omega1_ai,n_ov,1)
   call deallocator(u_ai_ck,n_ov,n_ov)

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
!     Both are permuted added to the omega vector element omega2(ai,bj)
!
      implicit none 
!
      logical :: debug = .false.
!
      integer :: b=0,c=0,k=0,d=0,ck=0,ckdl=0,cl=0,cldk=0,dk=0,dl=0,kc=0,kdl=0,l=0,ld=0
      integer :: a=0,ai=0,aibj=0,bj=0,aicj=0,cj=0,i=0,j=0,jai=0,dlc=0,dkcl=0,dlck=0,aib=0,aibk=0,bk=0,bja=0,ibj=0
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
      call get_cholesky_ia(L_kc_J)
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
!     Determine u_b_kdl = u_kl^bd and g_kdl_c = g_ldkc
!
            write(luprint,*) 'Lalala 1.2'
      call flshfo(luprint)
      do c = 1,n_vir ! Use as though "b" for u_b_kdl term 
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
                  write(luprint,*) 'Lalala 1.3'
      call flshfo(luprint)
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
            write(luprint,*) 'Lalala 1.4'
      call flshfo(luprint)
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
                  write(luprint,*) 'Lalala 1.5'
      call flshfo(luprint)
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
            write(luprint,*) 'Lalala 1.6'
      call flshfo(luprint)
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
                  write(luprint,*) 'Lalala 1.7'
      call flshfo(luprint)
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
            write(luprint,*) 'Lalala 1.8'
      call flshfo(luprint)
      do a = 1,n_vir
         do i = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  jai  = index_three(j,a,i,n_occ,n_vir)
                  ibj  = index_three(i,b,j,n_occ,n_vir)
!
                  ai   = index_two(a,i,n_vir)
                  bj   = index_two(b,j,n_vir)
                  aibj = index_packed(ai,bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
                     omega2(aibj,1) = omega2(aibj,1) + omega2_b_jai(b,jai) + omega2_b_jai(a,ibj)
                  endif
!
               enddo
            enddo
         enddo
      enddo
                  write(luprint,*) 'Lalala 1.9'
      call flshfo(luprint)
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
      call get_cholesky_ia(L_kc_J)
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
            write(luprint,*) 'Lalala 2.0'
      call flshfo(luprint)
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
                  write(luprint,*) 'Lalala 2.1'
      call flshfo(luprint)
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
            write(luprint,*) 'Lalala 2.2'
      call flshfo(luprint)
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
                  write(luprint,*) 'Lalala 2.3'
      call flshfo(luprint)
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
      write(luprint,*)'Dgemm in E2 done'
      call flshfo(luprint)
!     Deallocate t_aib_k and Y_k_j 
!
      call deallocator(t_aib_k,n_ovv,n_occ)
      call deallocator(Y_k_j,n_occ,n_occ)
!
!     Add the E2.2 term to the omega vector 
!
            write(luprint,*) 'Lalala 2.4'
      call flshfo(luprint)
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
                  bja  = index_three(b,j,a,n_vir,n_occ) 
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
                     omega2(aibj,1) = omega2(aibj,1) + omega2_aib_j(aib,j) + omega2_aib_j(bja,i)
                  endif
!
               enddo
            enddo
         enddo
      enddo
                  write(luprint,*) 'Lalala 2.5'
      call flshfo(luprint)
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
!     The first term is referred to as the D2.1 term, and comes out ordered as (ai,bj) 
!     The second term is referred to as the D2.2 term, and comes out ordered as (ai,bj)
!     The third term is referred to as the D2.3 term, and comes out ordered as (ai,bj)
!
!     All terms are added to the omega vector element omega2(ai,bj)
!     The routine adds the terms in the following order: D2.3, D2.1, D2.2
!
      implicit none 
!
      logical :: debug = .false.
!
      integer :: required=0,available=0,max_batch_length=0,batch_dimension=0,n_batch=0
!
      integer :: a_begin=0,a_end=0,a_batch=0,batch_length=0,a_full=0,ac_dim=0 
!
      integer :: ai=0,aidl=0,al=0,aldi=0,a=0,i=0,b=0,j=0,c=0,d=0,di=0,dl=0,k=0,kc=0,kd=0,l=0
      integer :: lc=0,ld=0,aibj=0,bj=0,bjck=0,bk=0,bkcj=0,cj=0,ck=0,ca=0,ki=0
!
      real(dp), dimension(:,:), pointer :: L_kc_J       => null() ! L_kc^J 
      real(dp), dimension(:,:), pointer :: g_ld_kc      => null() ! g_ldkc 
      real(dp), dimension(:,:), pointer :: L_ld_kc      => null() ! L_ldkc = 2 * g_ldkc - g_lckd 
      real(dp), dimension(:,:), pointer :: u_ai_ld      => null() ! u_il^ad = 2 * t_il^ad - t_li^ad 
      real(dp), dimension(:,:), pointer :: Z_ai_kc      => null() ! An intermediate, see below
!
      real(dp), dimension(:,:), pointer :: omega2_ai_bj => null() ! For storing the D2.3, D2.2 & D2.1 terms temporarily
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
      call get_cholesky_ia(L_kc_J)
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
!
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
!     Allocate the D2.3 term omega2_ai_bj and set it to zero (remove, I think)
!
     ! call allocator(omega2_ai_bj,n_ov,n_ov)
     ! omega2_ai_bj = zero  ! Debugging, Eirik, 16 Mar 2017
!
!     Allocate the intermediate Z_ai_kc = sum_dl u_ai_ld L_ld_kc and set it to zero
!
      call allocator(Z_ai_kc,n_ov,n_ov)
      Z_ai_kc = zero
!
!     Form the intermediate Z_ai_kc = sum_dl u_ai_ld L_ld_kc
!
      call dgemm('N','N',n_ov,n_ov,n_ov,&
                  one,u_ai_ld,n_ov,L_ld_kc,n_ov,&
                  zero,Z_ai_kc,n_ov)
!
!     Deallocate L_ld_kc
!
      call deallocator(L_ld_kc,n_ov,n_ov) ! Remove the one below...!
!
!     Allocate the D2.3 term omega2_ai_bj and set it to zero
!
      call allocator(omega2_ai_bj,n_ov,n_ov)
      omega2_ai_bj = zero
!
!     Form the D2.3 term, 1/4 sum_kc Z_ai_kc u_kc_bj = 1/4 sum_kc Z_ai_kc(ai,kc) u_ai_ld(bj,kc)
!
      call dgemm('N','T',n_ov,n_ov,n_ov,&
                  one/four,Z_ai_kc,n_ov,u_ai_ld,n_ov,& 
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
!     Deallocate the Z_ai_kc intermediate 
!
      call deallocator(Z_ai_kc,n_ov,n_ov)
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
!
                  aibj = index_packed(ai,bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
                     omega2(aibj,1) = omega2(aibj,1) + omega2_ai_bj(ai,bj) + omega2_ai_bj(bj,ai)
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the omega2_ai_bj and u_ai_ld(ai,ld) = u_il^ad vector
!
     ! call deallocator(L_ld_kc,n_ov,n_ov) ! Eirik: debugging (15 Mar, 2017)
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
      call get_cholesky_ai(L_ai_J)
      call get_cholesky_ia(L_kc_J)
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
                     omega2(aibj,1) = omega2(aibj,1) + omega2_ai_bj(ai,bj) + omega2_ai_bj(bj,ai)
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
      call get_cholesky_ij(L_ki_J)
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
         call get_cholesky_ab(L_ca_J,a_begin,a_end,ac_dim,.true.)
!
!        Allocate the integral g_ca_ki = g_acki and set to zero 
!
         call allocator(g_ca_ki,ac_dim,n_oo) ! Eirik: debugging: n_oo was n_ov!! ( 15 Mar, 2017)
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
                     omega2(aibj,1) = omega2(aibj,1) + omega2_ai_bj(ai,bj) + omega2_ai_bj(bj,ai)
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
      call deallocator(L_ki_J,n_oo,n_J) ! Eirik: debugging (15 Mar, 2017)
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
      integer                                :: c=0,k=0,d=0,l=0,kd=0,lc=0,ck=0,dl=0,al=0,di=0,ai=0,aldi=0,cl=0,dk=0,cldk=0
      integer                                :: i=0,a=0,ca=0,ki=0,b=0,bj=0,bk=0,cj=0,j=0,bkcj=0,aibj=0,aj=0,bi=0
      integer                                :: required=0,available=0,n_batch=0,max_batch_length=0,a_start=0
      integer                                :: a_end=0,a_batch=0,a_length=0
      logical                                :: debug = .false.
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
!     Deallocate
!
      call deallocator(L_ia_J,n_ov,n_J)
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
                  cl=index_two(c,l,n_vir)
                  dk=index_two(d,k,n_vir)
                  cldk=index_packed(cl,dk)
!
                  g_dl_ck(dl,ck)=g_kd_lc(kd,lc)
!
                  t_ai_dl(ck,dl)=t2am(cldk,1)
               enddo
            enddo
         enddo
      enddo

      call deallocator(g_kd_lc,n_ov,n_ov)
!
!     -1/2 * sum_(dl) t_ai_dl*g_dl_ck = X_ai_ck
!
      call allocator(X_ai_ck,n_ov,n_ov)
!
      call dgemm('N','N',n_ov,n_ov,n_ov &
         ,-half,t_ai_dl,n_ov,g_dl_ck,n_ov &
         ,zero,X_ai_ck,n_ov)
!
!     Deallocate L_ia_J and g_dl_ck, 
!
      call deallocator(g_dl_ck,n_ov,n_ov)
      call deallocator(t_ai_dl,n_ov,n_ov)
!
!     Constructing g_ki_ac ordered as g_ki_ca
!
!     Allocate g_ki_ca
!
      call allocator(g_ki_ca,n_oo,n_vv)
      g_ki_ca = zero
!
!     Allocate L_ki_J
!
      call allocator(L_ki_J,n_oo,n_J) 
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
      call deallocator(L_ki_J,n_oo,n_J)
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
                  bkcj = index_packed(bk,cj)
!
                  bj = index_two(b,j,n_vir)
                  ck = index_two(c,k,n_vir)
!
                  t_ck_bj(ck,bj) = t2am(bkcj,1)
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
!     Deallocate the X intermediate
!
      call deallocator(X_ai_ck,n_ov,n_ov)
!
!     Deallocate t_ck_bj
!
      call deallocator(t_ck_bj,n_ov,n_ov)
!
!     Omega_aibj,1 =P_ai_bj( 1/2*Y_ai_bj + Y_aj_bi)
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
                     omega2(aibj,1)=omega2(aibj,1)+half*Y_ai_bj(ai,bj)+Y_ai_bj(aj,bi)+half*Y_ai_bj(bj,ai)+Y_ai_bj(bi,aj)
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
!     Print the omega vector, having added C2
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after C2 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,n_ov_ov_packed)
      endif 
!
   end subroutine mlcc_omega_c2
!
   subroutine mlcc_omega_a2
!
!     MLCC Omega A2 term
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 10 Mar 2017
!
!     Omega A2 = g_ai_bj + sum_(cd)g_ac_bd * t_ci_dj = A2.1 + A.2.2
!
!     Structure: Batching over both a and b. If no batching is necessary L_ab_J is only read once, and g_ac_bd 
!                is constructed and kept in memory full size. 
!                g_ac_bd is reordered as g_ab_cd and t_ci_dj is reordered as t_cd_ij.
!                Omega contribution for A2.2 is ordered as Omega_ab_ij, and is reordered into the packed omega2 vector.          
!
      implicit none
!
      real(dp),dimension(:,:),pointer     :: g_ai_bj => null()
      real(dp),dimension(:,:),pointer     :: g_ca_db => null()
      real(dp),dimension(:,:),pointer     :: g_ab_cd => null()
      real(dp),dimension(:,:),pointer     :: t_cd_ij => null()
      real(dp),dimension(:,:),pointer     :: omega2_ab_ij => null()
      real(dp),dimension(:,:),pointer     :: L_ai_J => null()
      real(dp),dimension(:,:),pointer     :: L_ca_J => null()
      real(dp),dimension(:,:),pointer     :: L_db_J => null()
      integer                             :: a=0,i=0,b=0,j=0,ai=0,bj=0,aibj=0,c=0,d=0,ca=0,db=0,ab=0,cd=0,ci=0,cidj=0,ij=0,dj=0
      integer                             :: a_n_batch=0,b_n_batch=0,a_start=0,a_end=0,b_start=0,b_end=0,a_length=0,b_length=0
      integer                             :: required=0,available=0,a_max_length=0,b_max_length=0,a_batch=0,b_batch=0
      logical                             :: debug = .false.
!
!!!   A2.1 term   !!!
!
!
!     Create g_ai_bj 
!  
      call allocator(g_ai_bj,n_ov,n_ov)
      call allocator(L_ai_J,n_ov,n_J)
!
      call get_cholesky_ai(L_ai_J)
!
!     g_ai_bj=sum_J L_ai_J*L_bj_J
!
      call dgemm('N','T',n_ov,n_ov,n_J &
         ,one, L_ai_J,n_ov,L_ai_J,n_ov &
         ,zero,g_ai_bj,n_ov)
!
      call deallocator(L_ai_J,n_ov,n_J)
!
!     Add A2.1 to Omega 2
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
                  if(ai .ge. bj) then
!
                     aibj=index_packed(ai,bj)
!
                     omega2(aibj,1)=omega2(aibj,1)+g_ai_bj(ai,bj)
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ai_bj,n_ov,n_ov)
!
!
!!!   A2.2 term !!!
!
!
!     Create g_ac_bd by batching over both a and b
!
!
!     Prepare for batching over a
!
!
!     How many a-batches?
!
      required = 2*n_vv*n_J*4 + 2*n_ov*n_J*4 + 2*n_vv*n_vv*4 ! Needed to get cholesky of ab-type and for g_ab_bd
      available=get_available()
!
      a_max_length=0
!
      call n_one_batch(required,available,a_max_length,a_n_batch,n_vir)
!
!     Initialize some variables for batching
!
      a_start=0
      a_end=0
      a_length=0
!
!     Start looping over a-batches
!
      do a_batch = 1,a_n_batch
!
         call one_batch_limits(a_start,a_end,a_batch,a_max_length,n_vir)
         a_length=a_end-a_start+1
!
!        Get cholesky vectors L_ac^J ordered as L_ca_J
!
         call allocator(L_ca_J,n_vir*a_length,n_J)
!
         call get_cholesky_ab(L_ca_J,a_start,a_end,n_vir*a_length,.true.) 
!        
!        Prepare for batching over b
!
!
!        How many batches ?
!
         required = 2*n_vv*n_J*4 + 2*n_ov*n_J*4 + 2*n_vv*a_length*n_vir*4 ! Needed to get cholesky of ab-type 
                                                                          ! and to generate g_ac_bd.
         available=get_available()
!  
         b_max_length=0
!  
         call n_one_batch(required,available,b_max_length,b_n_batch,n_vir)
!  

         if (a_n_batch .eq. 1 .and. b_n_batch .eq. 1 ) then ! No need to read ab-type cholesky twice. 
                                                            ! Note that this means that a_length=b_length=n_vir
!
!           Allocate full g_ca_db because we have room for it!
!
            call allocator(g_ca_db,n_vv,n_vv)
!
            call dgemm('N','T',n_vv,n_vv,n_J &
               ,one,L_ca_J,n_vv,L_ca_J,n_vv &
               ,zero,g_ca_db,n_vv)
!
!           Allocate for reordered g_ca_db(g_ab_cd) and t_ci_dj (t_cd_ij)
!
            call allocator(g_ab_cd,n_vv,n_vv)
            call allocator(t_cd_ij,n_vv,n_oo)
            do c = 1,n_vir
               do d = 1,n_vir
!
!                 Reorder g_ca_db to g_ab_cd
!
                  do a = 1,n_vir
                     do b = 1,n_vir
!
!                       Needed indices
!  
                        ca=index_two(c,a,n_vir)
                        db=index_two(d,b,n_vir)
                        ab=index_two(a,b,n_vir)
                        cd=index_two(c,d,n_vir)
!
                        g_ab_cd(ab,cd)=g_ca_db(ca,db)
                     enddo
                  enddo
!
!                 Reorder t_ci_dj to t_cd_ij
!
                  do i = 1,n_occ
                     do j = 1,n_occ
!
!                       Needed indices
!  
                        cd=index_two(c,d,n_vir)
                        ij=index_two(i,j,n_occ)
                        ci=index_two(c,i,n_vir)
                        dj=index_two(d,j,n_vir)
!
                        cidj=index_packed(ci,dj)
!
                        t_cd_ij(cd,ij)=t2am(cidj,1)
!
                     enddo
                  enddo
               enddo
            enddo
            call deallocator(g_ca_db,n_vv,n_vv)
!
            call allocator(omega2_ab_ij,n_vv,n_oo)
!
!           omega2_ab_ij= sum_(cd) g_ab_cd*t_cd_ij
!
            call dgemm('N','N',n_vv,n_oo,n_vv &
               ,one,g_ab_cd,n_vv,t_cd_ij,n_vv &
               ,zero,omega2_ab_ij,n_vv)
!
!           Deallocate t_cd_ij and g_ab_cd
!
            call deallocator(g_ab_cd,n_vv,n_vv)
            call deallocator(t_cd_ij,n_vv,n_oo)
!
!           Reorder into omega2
!
               do i = 1,n_occ
                  do j = 1,n_occ
                     do a = 1,n_vir
                        do b = 1,n_vir
!
!                          Needed indices
!  
                           ai=index_two(a,i,n_vir)
                           bj=index_two(b,j,n_vir)
                           if (ai .ge. bj) then
                              ab=index_two(a,b,n_vir)
                              ij=index_two(i,j,n_occ)
!  
                              aibj=index_packed(ai,bj)
                              omega2(aibj,1)=omega2(aibj,1)+omega2_ab_ij(ab,ij)
!  
                           endif
                        enddo
                     enddo
                  enddo
               enddo
!
            call deallocator(omega2_ab_ij,n_vv,n_oo)
!
         else
!
!           Initialize some variables for batching over b
!
            b_start=0
            b_end=0
            b_length=0  
!
!           Start looping over batches of b
!
            do b_batch = 1,b_n_batch
!
               call one_batch_limits(b_start,b_end,b_batch,b_max_length,n_vir)
               b_length=b_end-b_start+1
!
!              Get cholesky vectors L_bd^J ordered as L_db_J
!
               call allocator(L_db_J,n_vir*b_length,n_J)
!  
               call get_cholesky_ab(L_db_J,b_start,b_end,n_vir*b_length,.true.) 
!
!              Allocate g_ac_bd for batches of a and b, ordered as g_ca_db
!
               call allocator(g_ca_db,n_vir*a_length,n_vir*b_length)
!
!              g_ca_db = sum_J L_ca_J*L_db_J
!
               call dgemm('N','T',n_vir*a_length,n_vir*b_length,n_J &
                  ,one,L_ca_J,n_vir*a_length,L_db_J,n_vir*b_length &
                  ,zero,g_ca_db,n_vir*a_length)
!
               call deallocator(L_db_J,n_vir*b_length,n_J)
!
!              Reorder g_ca_db into g_ab_cd and t_ci_dj into t_cd_ij
!
               call allocator(g_ab_cd,a_length*b_length,n_vir*n_vir)
               call allocator(t_cd_ij,n_vv,n_oo)
!
               do c = 1,n_vir 
                  do d = 1,n_vir
!
!                    Reorder g_ca_db to g_ab_cd
!
                     do a = 1,a_length
                        do b = 1,b_length
!
!                          Needed indices
!
                           ca=index_two(c,a,n_vir)
                           db=index_two(d,b,n_vir)
                           ab=index_two(a,b,n_vir)
                           cd=index_two(c,d,n_vir)
!  
                           g_ab_cd(ab,cd)=g_ca_db(ca,db)
                        enddo
                     enddo
!
!                 Reorder t_ci_dj to t_cd_ij
!
                  do i = 1,n_occ
                     do j = 1,n_occ
!
!                          Needed indices
!     
                           cd=index_two(c,d,n_vir)
                           ij=index_two(i,j,n_occ)
                           ci=index_two(c,i,n_vir)
                           dj=index_two(d,j,n_vir)
!  
                           cidj=index_packed(ci,dj)
!  
                           t_cd_ij(cd,ij)=t2am(cidj,1)
!   
                        enddo
                     enddo
                  enddo
               enddo
!
               call deallocator(g_ca_db,n_vir*a_length,n_vir*b_length)
!
!
               call allocator(omega2_ab_ij,a_length*b_length,n_oo)
!  
!              omega2_ab_ij= sum_(cd) g_ab_cd*t_cd_ij
!  
               call dgemm('N','N',a_length*b_length,n_oo,n_vv &
                  ,one,g_ab_cd,a_length*b_length,t_cd_ij,n_vv &
                  ,zero,omega2_ab_ij,a_length*b_length)
!  
!              Deallocate t_cd_ij and g_ab_cd
!
               call deallocator(t_cd_ij,n_vv,n_oo)
               call deallocator(g_ab_cd,a_length*b_length,n_vir*n_vir)
!
!              Reorder into omega2_aibj
!  
               do i = 1,n_occ
                  do j = 1,n_occ
                     do a = 1,a_length
                        do b = 1,b_length
!  
!                          Needed indices
!     
                           Ai=index_two(a+a_start-1,i,n_vir) ! A is full-space a index
                           Bj=index_two(b+b_start-1,j,n_vir) ! B is full-space b index
                           if (Ai .ge. Bj) then
                              ab=index_two(a,b,a_length)
                              ij=index_two(i,j,n_occ)
!     
                              AiBj=index_packed(Ai,Bj)
                              omega2(AiBj,1)=omega2(AiBj,1)+omega2_ab_ij(ab,ij)
!     
                           endif
                        enddo
                     enddo
                  enddo
               enddo
!
               call deallocator(omega2_ab_ij,n_vv,n_oo)
!
            enddo
!
         endif
!
         call deallocator(L_ca_J,n_vv,n_J)
!
      enddo
!
!
!     Print the omega vector, having added A2
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after A2 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,n_ov_ov_packed)
      endif
!
   end subroutine mlcc_omega_a2
!
   subroutine mlcc_omega_b2
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
      real(dp),dimension(:,:),pointer     :: g_kc_ld     => null()
      real(dp),dimension(:,:),pointer     :: g_kl_cd     => null()
      real(dp),dimension(:,:),pointer     :: g_kl_ij     => null()
      real(dp),dimension(:,:),pointer     :: g_ki_lj     => null()
      real(dp),dimension(:,:),pointer     :: t_cd_ij     => null()
      real(dp),dimension(:,:),pointer     :: t_ab_kl     => null()
      real(dp),dimension(:,:),pointer     :: X_kl_ij     => null()
      real(dp),dimension(:,:),pointer     :: omega_ab_ij => null()
      real(dp),dimension(:,:),pointer     :: L_kc_J      => null()
      real(dp),dimension(:,:),pointer     :: L_ij_J      => null()
      integer                             :: c=0,d=0,k=0,l=0,i=0,j=0,kl=0,ij=0,ci=0,dj=0,kc=0,ld=0,cidj=0,cd=0,ki=0
      integer                             :: lj=0,ak=0,bl=0,akbl=0,ab=0,b=0,a=0,ai=0,bj=0,aibj=0
      logical                             :: debug = .false.
!
!     Read cholesky vector of type L_ij_J
!
      call allocator(L_ij_J,n_oo,n_J)
      L_ij_J = zero ! Eirik: debug
!
      call get_cholesky_ij(L_ij_J)
!
!     Create g_ki_lj = sum_J L_li_J*L_lj_J
!
      call allocator(g_ki_lj,n_oo,n_oo)
      g_ki_lj = zero ! Eirik: debug 
!
      call dgemm('N','T',n_oo,n_oo,n_J &
         ,one,L_ij_J,n_oo,L_ij_J,n_oo &
         ,zero,g_ki_lj,n_oo)
!
      call deallocator(L_ij_J,n_oo,n_J)
!
!     Reordering g_ki_lj to g_kl_ij
!
      call allocator(g_kl_ij,n_oo,n_oo)
      g_kl_ij = zero
!
      do k = 1,n_occ
         do l = 1,n_occ
            do i = 1,n_occ
               do j=1,n_occ
!
!                 Needed indices
!
                  ki=index_two(k,i,n_occ)
                  lj=index_two(l,j,n_occ)
                  kl=index_two(k,l,n_occ)
                  ij=index_two(i,j,n_occ)
!
                  g_kl_ij(kl,ij)=g_ki_lj(ki,lj)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ki_lj,n_oo,n_oo)
!
!     Read cholesky vectors of ia-type into L_kc_J
!
      call allocator(L_kc_J,n_ov,n_J)
!
      call get_cholesky_ia(L_kc_J)
!
!     Create g_ck_ld = sum_(J) L_kc_J*L_ld_J
!
      call allocator(g_kc_ld,n_ov,n_ov)
!
      call dgemm('N','T',n_ov,n_ov,n_J &
         ,one,L_kc_J,n_ov,L_kc_J,n_ov &
         ,zero,g_kc_ld,n_ov)
!
!     Deallocate cholesky vectors L_ck_J
!
      call deallocator(L_kc_J,n_ov,n_J)
!
!     Reorder g_kc_ld as g_kl_cd, also reordering t_ci_dj as t_cd_ij
!
      call allocator(t_cd_ij,n_vv,n_oo)
      call allocator(g_kl_cd,n_oo,n_vv)
!
      do k = 1,n_occ
         do c = 1,n_vir
            do l = 1,n_occ
               do d = 1,n_vir
!
!                 Needed indices
!

                  cd=index_two(c,d,n_vir)
                  kl=index_two(k,l,n_occ)
                  kc=index_two(k,c,n_occ)
                  ld=index_two(l,d,n_occ)
!
                  ci=index_two(c,k,n_vir)
                  dj=index_two(d,l,n_vir)
                  ij=index_two(k,l,n_occ)
                  cidj=index_packed(ci,dj)
!
                  g_kl_cd(kl,cd)=g_kc_ld(kc,ld)
                  t_cd_ij(cd,ij)=t2am(cidj,1)
               enddo
            enddo
         enddo
      enddo
      call deallocator(g_kc_ld,n_ov,n_ov)
!
      call dgemm('N','N',n_oo,n_oo,n_vv &
         ,one,g_kl_cd,n_oo,t_cd_ij,n_vv &
         ,one,g_kl_ij,n_oo)
!
      call deallocator(t_cd_ij,n_vv,n_oo)
      call deallocator(g_kl_cd,n_oo,n_vv)
!
!     Reorder t_ak_bl to t_ab_kl
!
      call allocator(t_ab_kl,n_vv,n_oo)
!
      do a = 1,n_vir
         do b = 1,n_vir
            do k = 1,n_occ
               do l=1,n_occ
!
!                 Needed indices
!
                  ak=index_two(a,k,n_vir)
                  bl=index_two(b,l,n_vir)
                  ab=index_two(a,b,n_vir)
                  kl=index_two(k,l,n_occ)
!
                  akbl=index_packed(ak,bl)
!
                  t_ab_kl(ab,kl)=t2am(akbl,1)
!
               enddo
            enddo
         enddo
      enddo
!
!     omega_ab_ij= sum_(kl) t_ab_kl*X_kl_ij
!
      call allocator(omega_ab_ij,n_vv,n_oo)
!
      call dgemm('N','N',n_vv,n_oo,n_oo &
         ,one,t_ab_kl,n_vv,g_kl_ij,n_oo &
         ,zero,omega_ab_ij,n_vv)
!
      call deallocator(t_ab_kl,n_vv,n_oo)
      call deallocator(g_kl_ij,n_oo,n_oo)
!
!     Reorder omega
!
      do a = 1,n_vir
         do b = 1,n_vir
            do i = 1,n_occ
               do j = 1,n_occ
!
!                 Needed indices
!
                  ai=index_two(a,i,n_vir)
                  bj=index_two(b,j,n_vir)
                  if (ai .ge. bj) then
                     ab=index_two(a,b,n_vir)
                     ij=index_two(i,j,n_occ)
                     aibj=index_packed(ai,bj)
                     omega2(aibj,1)=omega2(aibj,1)+omega_ab_ij(ab,ij)
                  endif
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(omega_ab_ij,n_vv,n_oo)  
!   
!
!
!     Print the omega vector, having added B2
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after B2 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,n_ov_ov_packed)
      endif
!
   end subroutine mlcc_omega_b2
!
end module mlcc_omega