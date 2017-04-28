module ccs_class
!
   use input_output
   use workspace
   use mlcc_types
   use hf_class
!
   type, extends(hartree_fock) :: cc_singles
!
      integer(i15)                          :: n_t1am = 0 ! Number of singles amplitudes
      real(dp), dimension(:,:), allocatable :: t1am       ! Singles amplitudes
!
   contains 
!
      procedure, public  :: init               => init_cc_singles
      procedure, public  :: drv                => drv_cc_singles
!
      procedure :: fock_constructor   => fock_constructor_cc_singles
      procedure :: initialize_singles => initialize_singles_cc_singles
!
      procedure :: get_cholesky_ij => get_cholesky_ij_cc_singles
      procedure :: get_cholesky_ia => get_cholesky_ia_cc_singles
      procedure :: get_cholesky_ai => get_cholesky_ai_cc_singles
      procedure :: get_cholesky_ab => get_cholesky_ab_cc_singles
!
      procedure :: h_one_electron_mo_t1 => ham_one_electron_mo_t1_cc_singles
!
   end type cc_singles
!
!  private routines and variables for ccs_class
!
 !  private :: h1mo_T1
!
contains
!
   subroutine init_cc_singles(wfn)
!
      implicit none 
!
      class(cc_singles) :: wfn
!
      write(unit_output,*) 'In init_cc_singles'
!
!     Read Hartree-Fock info from SIRIUS
!
      call wfn % read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wfn % read_transform_cholesky
!
!     Initialize singles amplitudes and associated attributes
!
      call wfn % initialize_singles
!
   end subroutine init_cc_singles
!
   subroutine drv_cc_singles(wfn)
!
      implicit none 
!
      class(cc_singles) :: wfn
!
      write(unit_output,*) 'In drv_cc_singles'
!
   end subroutine drv_cc_singles
!
   subroutine initialize_singles_cc_singles(wfn)
!
      implicit none 
!
      class(cc_singles) :: wfn
!
!     Calculate the number of singles amplitudes
!
      wfn % n_t1am = (wfn % n_occ)*(wfn % n_vir) 
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator ( wfn % t1am, wfn % n_t1am, 1)
      wfn % t1am = zero
!
   end subroutine initialize_singles_cc_singles
!
   subroutine fock_constructor_cc_singles(wfn)
!
      use mlcc_oo_utilities
      use workspace
!
      implicit none
!
      class(cc_singles)                          :: wfn
!
      real(dp), dimension(:,:), allocatable      :: fock_ao
      real(dp), dimension(:,:), allocatable      :: fock_matrix
      real(dp), dimension(:,:), allocatable      :: ao_int 
      real(dp), dimension(:,:), allocatable      :: h1mo 
      integer                                    :: n_ao_sq_packed  
      real(dp), dimension(:,:), allocatable      :: X      
      integer                                    :: unit_identifier_ao_integrals = -1
      integer                                    :: i=0,j=0,k=0,a=0,b=0,ij=0,kk=0,ik=0,kj=0,ii=0,jj=0
      integer                                    :: ji=0,ai=0,ib=0,bi=0,ia=0,aj=0,ja=0,ab=0
      integer                                    :: n_oo, n_ov, n_vv
      real(dp), dimension(:,:), allocatable      :: g_ij_kl
      real(dp), dimension(:,:), allocatable      :: g_ab_ij
      real(dp), dimension(:,:), allocatable      :: g_ai_jb
      real(dp), dimension(:,:), allocatable      :: g_ia_jk
      real(dp), dimension(:,:), allocatable      :: g_ai_jk
      real(dp), dimension(:,:), allocatable      :: L_ij_J 
      real(dp), dimension(:,:), allocatable      :: L_ia_J 
      real(dp), dimension(:,:), allocatable      :: L_ai_J 
      real(dp), dimension(:,:), allocatable      :: L_ab_J 
      integer                                    :: available=0,required=0,max_batch_length=0,n_batches=0,batch_start=0
      integer                                    :: batch_end=0,batch_length=0,g_off=0
      integer                                    :: b_batch = 0
!
!     Allocate one-electron integrals in MO basis
!
      call allocator(h1mo, wfn % n_mo, wfn % n_mo)
      call allocator(fock_matrix, wfn % n_mo, wfn % n_mo)
      h1mo = zero
!
!
!
!!!   ONE-ELECTRON CONTRIBUTION !!!
!
!
!     Allocate for one-electron ao integrals
!
      n_ao_sq_packed = packed_size(wfn % n_ao)
      call allocator(ao_int,n_ao_sq_packed,1)
!
!     Open mlcc_aoint
!
      open(unit=unit_identifier_ao_integrals,file='mlcc_aoint',status='old',form='formatted')
      rewind(unit_identifier_ao_integrals)
!
!     Read one-electron ao integrals
!
      read(unit_identifier_ao_integrals,*)(ao_int(i,1),i=1,n_ao_sq_packed)
!
!     Close mlcc_aoint
!
      close(unit_identifier_ao_integrals)
!
!
!     Allocate ao fock matrix, add one-electron contributions
!
      call allocator(fock_ao, wfn % n_ao, wfn % n_ao)
      call squareup(ao_int,fock_ao,wfn % n_ao)   
!
!     Deallocation of ao integrals
!   
      call deallocator(ao_int,n_ao_sq_packed,1) 
!
!     Transform to one-electron part to mo and add to fock_matrix 
!
      call allocator(X, wfn % n_ao, wfn % n_mo)
!
      call dgemm('N','N',         &
                  wfn%n_ao,    &
                  wfn%n_mo,    &
                  wfn%n_ao,    &
                  one,            &
                  fock_ao,        &
                  wfn%n_ao,    &
                  wfn%mo_coef, &
                  wfn%n_ao,    &
                  zero,           &
                  X,              &
                  wfn%n_ao)
!
      call dgemm('T','N',         &        
                  wfn%n_mo,    &
                  wfn%n_mo,    &
                  wfn%n_ao,    &
                  one,            &
                  wfn%mo_coef, &
                  wfn%n_ao,    &
                  X,              &
                  wfn%n_ao,    &
                  zero,           &
                  h1mo,           &
                  wfn%n_mo)
! 
!     T1-transformation of one-electron integrals in mo basis
!
      call wfn % ham_one_electron_mo_t1
 !     call h1mo_T1(h1mo,fock_matrix, wfn % n_occ, wfn % n_vir, wfn % n_mo, wfn % t1am, wfn % n_t1am)
      call deallocator(h1mo,wfn % n_mo,wfn % n_mo)
!
!     Deallocate intermediate X and fock_ao
!
      call deallocator(X,wfn % n_ao, wfn % n_mo)
      call deallocator(fock_ao, wfn % n_ao , wfn % n_ao)
!
!
!
!!! TWO-ELECTRON CONTRIBUTION !!!
!
!
!
!!  Occupied-Occupied block: F_ij = h_ij + sum_k (2*g_ijkk - g_ikkj) !!
!  
!
!     Allocation for L_ij_J_pack
! 
      n_oo = (wfn % n_occ) * (wfn % n_occ)
!
      call allocator(L_ij_J,n_oo, wfn % n_J)
      call allocator(g_ij_kl,n_oo,n_oo)
      L_ij_J=zero
      g_ij_kl=zero
!
!     Read cholesky vectors
!
!      call get_cholesky_ij(L_ij_J)
!
      call dgemm('N','T',     &
                  n_oo,       &
                  n_oo,       &
                  wfn%n_J, &
                  one,        &
                  L_ij_J,     &
                  n_oo,       &
                  L_ij_J,     &
                  n_oo,       &
                  zero,       &
                  g_ij_kl,    &
                  n_oo)
!
!     Add two-electron contributions to occupied-occupied block
!
      do i=1, wfn % n_occ
         do j=1,wfn % n_occ
            ij=index_two(i,j,wfn % n_occ)
            do k = 1, wfn % n_occ
               kk=index_two(k,k, wfn % n_occ)
               ik=index_two(i,k, wfn % n_occ)
               kj=index_two(k,j, wfn % n_occ)
               fock_matrix(i,j)=fock_matrix(i,j)+two*g_ij_kl(ij,kk)-g_ij_kl(ik,kj)
            enddo
         enddo
      enddo
!
!     Deallocate g_ij_kl
! 
      call deallocator(g_ij_kl,n_oo,n_oo)
!
!!    Occupied-vacant blocks
!
!
!     Allocation for g_ia_jk 
!
      n_ov = (wfn % n_occ) * (wfn % n_vir)
!
      call allocator(L_ia_J,n_ov, wfn % n_J)
      call allocator(g_ia_jk,n_ov,n_oo)
!
!     Reading Cholesky vector L_ia_J
!
!      call get_cholesky_ia(L_ia_J)
!
!     g_ia_jk
!
      call dgemm('N','T',     &
                  n_ov,       &
                  n_oo,       &
                  wfn%n_J, &
                  one,        &
                  L_ia_J,     &
                  n_ov,       &
                  L_ij_J,     &
                  n_oo,       &
                  zero,       &
                  g_ia_jk,    &
                  n_ov)
!
!     Dealllocate L_ia_J
!
      call deallocator(L_ia_J,n_ov, wfn % n_J)
!
!
!     Allocation for g_ai_jk 
!
      call allocator(L_ai_J,n_ov, wfn  % n_J)
      call allocator(g_ai_jk,n_ov,n_oo)
!
!     Reading Cholesky vector L_ai_J
!
!      call get_cholesky_ai(L_ai_J)
!
!     g_ai_jk
!
      call dgemm('N','T',     &
                  n_ov,       &
                  n_oo,       &
                  wfn%n_J, &
                  one,        &  
                  L_ai_J,     &
                  n_ov,       &
                  L_ij_J,     &
                  n_oo,       &
                  zero,       &
                  g_ai_jk,    &
                  n_ov)
!
!     Deallocate L_ai_J
!
      call deallocator(L_ai_J,n_ov, wfn % n_J)
!
!     Adding terms to Fock matrix
!
      do i=1, wfn % n_occ
         do a=1, wfn % n_vir
            do j=1, wfn % n_occ
!
!              Needed indices
!
               ia=index_two(i,a, wfn % n_occ)
               ja=index_two(j,a, wfn % n_occ)
!
               ai=index_two(a,i, wfn % n_vir)
               aj=index_two(a,j, wfn % n_vir)
!
               jj=index_two(j,j, wfn % n_occ)
               ji=index_two(j,i, wfn % n_occ)
               ij=index_two(i,j, wfn % n_occ)
! 
               fock_matrix(i,a + wfn%n_occ)=fock_matrix(i,a + wfn%n_occ)+two*g_ia_jk(ia,jj)-g_ia_jk(ja,ij) ! g_ia_jk(ja,ij) = g_jaij = g_ijja
               fock_matrix(a + wfn%n_occ,i)=fock_matrix(a + wfn%n_occ,i)+two*g_ai_jk(ai,jj)-g_ai_jk(aj,ji)
            enddo
         enddo
      enddo
      call deallocator(g_ia_jk,n_ov,n_oo)
      call deallocator(g_ai_jk,n_ov,n_oo)
!
!!    Vacant-vacant block F_ab = h_ab + sum_k (2*g_abkk - g_akkb) !!
!
      n_vv = (wfn % n_vir) * (wfn % n_vir)
      call allocator(g_ab_ij,n_vv,n_oo)
      g_ab_ij=zero
!
!     Batching over a
!
!
!     Setup of variables needed for batching
!
      available = get_available()
      required = 2*(wfn%n_vir)*(wfn%n_vir)*(wfn%n_J)*4 + 2*(wfn%n_vir)*(wfn%n_occ)*(wfn%n_J)*4
      call num_batch(required,available,max_batch_length,n_batches,wfn%n_vir)
!
      batch_start=1
      batch_end=0
      batch_length=0
!
!     Start batchig loop
!
      do b_batch = 1,n_batches
!
!        Get batch limits  and  length of batch
!
         call batch_limits(batch_start,batch_end,b_batch,max_batch_length,wfn%n_vir)
         batch_length=batch_end-batch_start+1
!
!        Allocation of L_ab_J
!
         call allocator(L_ab_J,wfn%n_vir*batch_length,wfn%n_J)
         L_ab_J=zero
!
!        Read Cholesky vectors
!
!         call get_cholesky_ab(L_ab_J,batch_start,batch_end,(wfn%n_vir)*batch_length,.false.)
!
!        g_ab_ij=sum_J L_ab_J* L_ij_J
!
         g_off = index_two(1,batch_start,wfn%n_vir)
!
         call dgemm('N','T',                      &
                     (wfn%n_vir)*batch_length, &
                     n_oo,                        &
                     wfn%n_J,                  &
                     one,                         &
                     L_ab_J,                      &
                     (wfn%n_vir)*batch_length, &
                     L_ij_J,                      &
                     n_oo,                        &
                     one,                         &
                     g_ab_ij(g_off,1),            &
                     n_vv)
!
!
!        Deallocation of L_ab_J
!
         call deallocator(L_ab_J,batch_length*(wfn%n_vir),wfn%n_J)
!
      enddo ! batching done
!
!     Deallocation of L_ij_J
!
      call deallocator(L_ij_J,n_oo,wfn%n_J)
!
!     Allocate for g_ai_jb
!
      call allocator(g_ai_jb,n_ov,n_ov)
      call allocator(L_ai_J,n_ov,wfn%n_J)
      call allocator(L_ia_J,n_ov,wfn%n_J)
!
!
!     Reading Cholesky vector L_ia_J and L_ai_J
!
!      call get_cholesky_ai(L_ai_J)
!
      call dgemm('N','T',     &
                  n_ov,       &
                  n_ov,       &
                  wfn%n_J, &
                  one,        &
                  L_ai_J,     &
                  n_ov,       &
                  L_ia_J,     &
                  n_ov,       &
                  zero,       &
                  g_ai_jb,    &
                  n_ov)
!
!     Deallocate L_ia_J
!
     call deallocator(L_ia_J,n_ov,wfn%n_J)
     call deallocator(L_ai_J,n_ov,wfn%n_J)
!
!     Calculation of two-electron terms for virtual-virtual blocks
!
      do a = 1,wfn%n_vir
         do b = 1,wfn%n_vir
            ab=index_two(a,b, wfn%n_vir)
            do i = 1,wfn%n_occ
               ii=index_two(i,i, wfn%n_occ)
               ai=index_two(a,i, wfn%n_vir)
               bi=index_two(b,i, wfn%n_vir)
               ia=index_two(i,a, wfn%n_occ)
               ib=index_two(i,b, wfn%n_occ)
               fock_matrix(wfn%n_occ+a,wfn%n_occ+b)=fock_matrix(wfn%n_occ+a,wfn%n_occ+b) &
                           +two*g_ab_ij(ab,ii)-g_ai_jb(ai,ib)
            enddo
         enddo 
      enddo
     call deallocator(g_ab_ij,n_vv,n_oo)
     call deallocator(g_ai_jb,n_ov,n_ov)
!
!     Save the blocks of the Fock matrix in memory (ij,ia,ai,ab)
!
!
      do i = 1,wfn%n_occ
         do j = 1,wfn%n_occ
            wfn%fock_matrix_ij(i,j) = fock_matrix(i,j)
         enddo
      enddo
!
      do i = 1,wfn%n_occ
         do a = 1,wfn%n_vir
            wfn%fock_matrix_ia(i,a) = fock_matrix(i,wfn%n_occ+a)
            wfn%fock_matrix_ai(a,i) = fock_matrix(wfn%n_occ+a,i)
         enddo
      enddo
!
      do a = 1,wfn%n_vir
         do b = 1,wfn%n_vir
            wfn%fock_matrix_ab(a,b) = fock_matrix(wfn%n_occ+a,wfn%n_occ+b)
         enddo
      enddo
!
   call deallocator(fock_matrix,wfn%n_mo,wfn%n_mo)
!
   end subroutine fock_constructor_cc_singles
!
   subroutine ham_one_electron_mo_t1_cc_singles(wfn,h1,h1_T1)
!
!  Purpose: T1-transform of one-electron mo integrals (h1mo)
!
!           h_p_q_T1= sum_(st) x_p_s * y_q_t * h_s_t
!           x = I-t1
!           y = I-t1^T
!
!
   use workspace
   use mlcc_types
!
   implicit none
!
   class(cc_singles) :: wfn
!
   double precision h1(wfn % n_mo, wfn % n_mo)
   double precision h1_T1(wfn % n_mo, wfn % n_mo)
!
   real(dp), dimension(:,:), allocatable     :: x 
   real(dp), dimension(:,:), allocatable     :: y 
   real(dp), dimension(:,:), allocatable     :: t1
   real(dp), dimension(:,:), allocatable     :: Z ! Intermediate for matrix multiplication
   integer(i15)                              :: p=0,q=0,a=0,i=0
!
!  Create t1, x, and y
!
   call allocator(t1, wfn % n_mo, wfn % n_mo)
   t1 = zero
!
   call allocator(y, wfn % n_mo, wfn % n_mo)
   call allocator(x, wfn % n_mo, wfn % n_mo)
!
   x = zero
   y = zero
!
!  t1_p_q = t1am_p_q for p virtual and q occupied, 0 otherwize
!
   do a = 1,n_vir
      do i = 1,n_occ
         t1(n_occ+a,i)=t1am(a,i)
      enddo
   enddo
!
   do p = 1,n_orbitals
      do q = 1,n_orbitals
         if (p .eq. q) then
            x(p,q) = 1
            y(p,q) = 1
         else
            x(p,q)=x(p,q)-t1(p,q)
            y(p,q)=y(p,q)+t1(q,p)
         endif
      enddo
   enddo
!
!  Deallocate t1
!
   call deallocator(t1, wfn % n_mo, wfn % n_mo)
!
!  Allocate Intermediate
!   
   call allocator(Z, wfn % n_mo, wfn % n_mo)
!
!  h1_T1 = x*h1*y^T = x*Z
!
   call dgemm('N','T',        &
               wfn % n_mo, &
               wfn % n_mo, &
               wfn % n_mo, &
               one,           &
               h1,            &
               wfn % n_mo, &
               y,             &
               wfn % n_mo, &
               zero,          &
               Z,             &
               wfn % n_mo)
!
   call dgemm('N','N',        &
               wfn % n_mo, &
               wfn % n_mo, &
               wfn % n_mo, &
               one,           &
               x,             &
               wfn % n_mo, &
               Z,             &
               wfn % n_mo, &
               zero,          &
               h1_T1,         &
               wfn % n_mo)
!
!  Deallocations
!
   call deallocator(Z, wfn % n_mo, wfn % n_mo)
   call deallocator(y, wfn % n_mo, wfn % n_mo)
   call deallocator(x, wfn % n_mo, wfn % n_mo)
!
   end subroutine ham_one_electron_mo_t1_cc_singles
!
   subroutine get_cholesky_ij_cc_singles(wfn,L_ij_J)
!
      implicit none 
!
      class(cc_singles) :: wfn
!
      real(dp), dimension((wfn%n_occ)**2,wfn%n_J) :: L_ij_J
!
      real(dp), dimension(:,:), allocatable :: L_ia_J 
      real(dp), dimension(:,:), allocatable :: L_iJ_a ! L_ia^J
      real(dp), dimension(:,:), allocatable :: L_iJ_k ! L_ik^J 
!
      integer(i15) :: i=0,J=0,a=0,ij=0,ia=0,ik=0,k=0
!
      integer(i15) :: n_ov 
!
!     Calculate n_occ * n_vir = n_ov
!
      n_ov = (wfn % n_occ)*(wfn % n_vir)
!
!     Allocate
!
      call allocator(L_ia_J,n_ov,wfn%n_J)
      call allocator(L_iJ_a,(wfn%n_occ)*(wfn%n_J),wfn%n_vir))

   end subroutine get_cholesky_ij_cc_singles
!
   subroutine get_cholesky_ia_cc_singles(wfn,L_ia_J)
!
      implicit none 
!
      class(cc_singles) :: wfn
!
      real(dp), dimension((wfn % n_occ)*(wfn % n_vir), wfn % n_J) :: L_ia_J
!
   end subroutine get_cholesky_ia_cc_singles
!
   subroutine get_cholesky_ai_cc_singles(wfn)
!
      implicit none 
!
      class(cc_singles) :: wfn
!
   end subroutine get_cholesky_ai_cc_singles
!
   subroutine get_cholesky_ab_cc_singles(wfn)
!
      implicit none 
!
      class(cc_singles) :: wfn
!
   end subroutine get_cholesky_ab_cc_singles
!
end module ccs_class
