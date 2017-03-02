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
      integer :: i,j,idummy
!
      integer required,available,max_batch_length,batch_dimension,n_batch
!
      integer a_begin,a_end
!
      real(dp), dimension(:,:), pointer :: L_kc_J  => null()
      real(dp), dimension(:,:), pointer :: L_ad_J  => null()   ! Here, a is being batched over
      real(dp), dimension(:,:), pointer :: g_ad_kc => null()   ! Here, a is being batched over
      real(dp), dimension(:,:), pointer :: g_a_ckd => null()   ! Here, a is being batched over 
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
!     Calculate the batching parameters over a = 1,2,...,n_vir,
!     for which we need to have enough room to store L_ad_J and g_ad_kc
!
      required = n_vir*n_vir*n_J + n_vir*n_vir*n_occ*n_vir
      available = get_available()
      batch_dimension = n_vir ! Batch over the virtual index a
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
!        For each batch, loop over the a indices in the given batch 
!
         a_begin = 1 + (a_batch-1)*max_batch_length
         a_end   = min(max_batch_length+(a_batch-1)*max_batch_length,n_vir)
!
      enddo
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
      integer :: lucho_ij,lucho_ia,idummy,j,i
!
      logical :: debug = .true.
!
      integer :: a,c,k,l,ckl,ki,cl,ak,akcl,al,alck,ck,ai
!
      real(dp), dimension(:,:), pointer :: L_ki_J_packed => null() ! (ki) index is packed
      real(dp), dimension(:,:), pointer :: L_lc_J        => null()
      real(dp), dimension(:,:), pointer :: g_ki_lc       => null()
      real(dp), dimension(:,:), pointer :: g_ckl_i       => null() ! Reordered two-electron integrals
      real(dp), dimension(:,:), pointer :: u_a_ckl       => null() ! Reordered u_kl^ac = 2 t_kl^ac - t_lk^ac
      real(dp), dimension(:,:), pointer :: b1_a_i        => null() 
!
      write(ml_lupri,*) 'In mlcc_omega_b1'
!
!     I. Calculation and reordering of g_ki,lc = sum_J L_ki^J * L_lc^J 
!
!     Allocate Cholesky vectors L_ki,J and L_lc,J 
!
      call allocator(L_ki_J_packed,n_occ*(n_occ+1)/2,n_J)
      call allocator(L_lc_J,n_occ*n_vir,n_J)
!
      L_ki_J_packed = zero
      L_lc_J = zero
!
!     Read L_ki,J
!
      call read_cholesky_ij(L_ki_J_packed)
!
!     Read L_lc,J
!
      call read_cholesky_ia(L_lc_J)
!
!     Allocate integrals g_ki_lc
!
      call allocator(g_ki_lc,n_oo_packed,n_ov) ! (ki) index is packed 
!
      g_ki_lc = zero
!
!     Calculate g_ki_lc = sum_J L_ki,J * L_lc,J^T 
! 
      call dgemm('N','T',n_oo_packed,n_ov,n_J,one,L_ki_J_packed,n_oo_packed,L_lc_J,n_ov,zero,g_ki_lc,n_oo_packed) 
!
!     Deallocate the Cholesky vectors 
!
      call deallocator(L_ki_J_packed,n_oo,n_J)
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
                  ki  = index_packed(k,i)                    ! Packed 
                  cl  = index_two(c,l,n_vir)             ! Not packed
!
!                 Set value of g_ckl_i
!
                  g_ckl_i(ckl,i) = g_ki_lc(ki,cl)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate unordered integrals g_ki_lc
!
      call deallocator(g_ki_lc,n_oo_packed,n_ov)
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
!                 Set value of u_a_ckl
!
                  u_a_ckl(a,ckl) = two*t2am(akcl,1) - t2am(alck,1)
!                  
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the B1 term
!
      call allocator(b1_a_i,n_vir,n_occ)
!
!     Calculate the B1 term (b1_a_i = - sum_ckl u_a_ckl * g_ckl_i)
!
      call dgemm('N','N',n_vir,n_occ,n_oov,-one,u_a_ckl,n_vir,g_ckl_i,n_oov,zero,omega1,n_vir) 
!
!     Print the omega vector 
!
      if (debug) call vec_print(omega1,n_vir,n_occ)
!
   end subroutine mlcc_omega_b1
!
   subroutine mlcc_omega_c1
!
!  Purpose: Calculate C1 term of Omega
!  (Omega_ai^C1=sum_ck F_kc u^ac_ik)
!
!  Written by Sarai D. Folkestad and Eirik F. Kjønstad, 28. Feb 2017
!
   use mlcc_data
   use mlcc_utilities
   use mlcc_workspace
   implicit none
!
   real(dp),dimension(:,:),pointer  :: F_kc
   real(dp),dimension(:,:),pointer  :: u_ck_ai
   integer                          :: i,k,c,a
   integer                          :: nck,nai,nak,nci,nckai,nciak
!
!
!  Allocation
!
   call allocator(u_ck_ai,n_ov,n_ov)   ! 2*t_ck_ai-t_ci_ak
!
! Set up u_ck_ai
!
   do i=1,n_occ
      do a=1,n_vir
         do k=1,n_occ
            do c=1,n_vir
               nck=index_two(c,k,n_vir)
               nai=index_two(a,i,n_vir)
               if (nck .le. nai) then
                  nci=index_two(c,i,n_vir)
                  nak=index_two(a,k,n_vir)
                  nckai=index_packed(nck,nai)
                  nciak=index_packed(nci,nak)
                  u_ck_ai(nck,nai)=two*t2am(nckai,1)-t2am(nciak,1)
               endif
            enddo
         enddo
      enddo
   enddo
!
! Allocation 
!
   call allocator(F_kc,n_ov,1) ! T1-transformed fock matrix
!
!  IO - Read fock matrix
!
!
!  T1 transformation
!
!
!  Deallocation
!
   call deallocator(F_kc,n_ov,1)
   call deallocator(u_ck_ai,n_ov,n_ov)

   end subroutine mlcc_omega_c1
!
   subroutine mlcc_omega_d1
   end subroutine mlcc_omega_d1
!
   subroutine mlcc_omega_e2
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