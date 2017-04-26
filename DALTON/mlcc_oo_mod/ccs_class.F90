module ccs_class
!
   use input_output
   use workspace
   use mlcc_types
   use hf_class
!
   type, extends(hartree_fock) :: cc_singles
!
      integer(i15) :: n_t1am = 0
      real(dp), dimension(:,:), allocatable :: t1am ! Singles amplitudes
!
   contains 
!
      procedure :: init             => init_cc_singles
      procedure :: drv              => drv_cc_singles
      procedure :: fock_constructor => fock_constructor_cc_singles
!
   end type cc_singles
!
!  private routines and variables for ccs_class
!
   private :: h1mo_T1
!
contains
!
   subroutine init_cc_singles(wavefn)
!
      implicit none 
!
      class(cc_singles) :: wavefn
!
      write(unit_output,*) 'In init_cc_singles'
!
!     Initialize the Hartree-Fock-specific quantities
!
      call init_hartree_fock(wavefn)
!
!     Calculate the number of singles amplitudes
!
      wavefn % n_t1am = (wavefn % n_occ)*(wavefn % n_vir) 
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator ( wavefn % t1am, wavefn % n_t1am, 1)
      wavefn % t1am = zero
!
   end subroutine init_cc_singles
!
   subroutine drv_cc_singles(wavefn)
!
      implicit none 
!
      class(cc_singles) :: wavefn
!
      write(unit_output,*) 'In drv_cc_singles'
!
   end subroutine drv_cc_singles
!
   subroutine fock_constructor_cc_singles(wavefn)
!
      use mlcc_oo_utilities
      use workspace
!
      implicit none
!
      class(cc_singles)                          :: wavefn
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
      call allocator(h1mo, wavefn % n_mo, wavefn % n_mo)
      call allocator(fock_matrix, wavefn % n_mo, wavefn % n_mo)
      h1mo = zero
!
!
!
!!!   ONE-ELECTRON CONTRIBUTION !!!
!
!
!     Allocate for one-electron ao integrals
!
      n_ao_sq_packed = packed_size(wavefn % n_ao)
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
      call allocator(fock_ao, wavefn % n_ao, wavefn % n_ao)
      call squareup(ao_int,fock_ao,wavefn % n_ao)   
!
!     Deallocation of ao integrals
!   
      call deallocator(ao_int,n_ao_sq_packed,1) 
!
!     Transform to one-electron part to mo and add to fock_matrix 
!
      call allocator(X, wavefn % n_ao, wavefn % n_mo)
!
      call dgemm('N','N',         &
                  wavefn%n_ao,    &
                  wavefn%n_mo,    &
                  wavefn%n_ao,    &
                  one,            &
                  fock_ao,        &
                  wavefn%n_ao,    &
                  wavefn%mo_coef, &
                  wavefn%n_ao,    &
                  zero,           &
                  X,              &
                  wavefn%n_ao)
!
      call dgemm('T','N',         &        
                  wavefn%n_mo,    &
                  wavefn%n_mo,    &
                  wavefn%n_ao,    &
                  one,            &
                  wavefn%mo_coef, &
                  wavefn%n_ao,    &
                  X,              &
                  wavefn%n_ao,    &
                  zero,           &
                  h1mo,           &
                  wavefn%n_mo)
! 
!     T1-transformation of one-electron integrals in mo basis
!
      call h1mo_T1(h1mo,fock_matrix, wavefn % n_occ, wavefn % n_vir, wavefn % n_mo, wavefn % t1am, wavefn % n_t1am)
      call deallocator(h1mo,wavefn % n_mo,wavefn % n_mo)
!
!     Deallocate intermediate X and fock_ao
!
      call deallocator(X,wavefn % n_ao, wavefn % n_mo)
      call deallocator(fock_ao, wavefn % n_ao , wavefn % n_ao)
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
      n_oo = (wavefn % n_occ) * (wavefn % n_occ)
!
      call allocator(L_ij_J,n_oo, wavefn % n_J)
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
                  wavefn%n_J, &
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
      do i=1, wavefn % n_occ
         do j=1,wavefn % n_occ
            ij=index_two(i,j,wavefn % n_occ)
            do k = 1, wavefn % n_occ
               kk=index_two(k,k, wavefn % n_occ)
               ik=index_two(i,k, wavefn % n_occ)
               kj=index_two(k,j, wavefn % n_occ)
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
      n_ov = (wavefn % n_occ) * (wavefn % n_vir)
!
      call allocator(L_ia_J,n_ov, wavefn % n_J)
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
                  wavefn%n_J, &
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
      call deallocator(L_ia_J,n_ov, wavefn % n_J)
!
!
!     Allocation for g_ai_jk 
!
      call allocator(L_ai_J,n_ov, wavefn  % n_J)
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
                  wavefn%n_J, &
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
      call deallocator(L_ai_J,n_ov, wavefn % n_J)
!
!     Adding terms to Fock matrix
!
      do i=1, wavefn % n_occ
         do a=1, wavefn % n_vir
            do j=1, wavefn % n_occ
!
!              Needed indices
!
               ia=index_two(i,a, wavefn % n_occ)
               ja=index_two(j,a, wavefn % n_occ)
!
               ai=index_two(a,i, wavefn % n_vir)
               aj=index_two(a,j, wavefn % n_vir)
!
               jj=index_two(j,j, wavefn % n_occ)
               ji=index_two(j,i, wavefn % n_occ)
               ij=index_two(i,j, wavefn % n_occ)
! 
               fock_matrix(i,a + wavefn%n_occ)=fock_matrix(i,a + wavefn%n_occ)+two*g_ia_jk(ia,jj)-g_ia_jk(ja,ij) ! g_ia_jk(ja,ij) = g_jaij = g_ijja
               fock_matrix(a + wavefn%n_occ,i)=fock_matrix(a + wavefn%n_occ,i)+two*g_ai_jk(ai,jj)-g_ai_jk(aj,ji)
            enddo
         enddo
      enddo
      call deallocator(g_ia_jk,n_ov,n_oo)
      call deallocator(g_ai_jk,n_ov,n_oo)
!
!!    Vacant-vacant block F_ab = h_ab + sum_k (2*g_abkk - g_akkb) !!
!
      n_vv = (wavefn % n_vir) * (wavefn % n_vir)
      call allocator(g_ab_ij,n_vv,n_oo)
      g_ab_ij=zero
!
!     Batching over a
!
!
!     Setup of variables needed for batching
!
      available = get_available()
      required = 2*(wavefn%n_vir)*(wavefn%n_vir)*(wavefn%n_J)*4 + 2*(wavefn%n_vir)*(wavefn%n_occ)*(wavefn%n_J)*4
      call num_batch(required,available,max_batch_length,n_batches,wavefn%n_vir)
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
         call batch_limits(batch_start,batch_end,b_batch,max_batch_length,wavefn%n_vir)
         batch_length=batch_end-batch_start+1
!
!        Allocation of L_ab_J
!
         call allocator(L_ab_J,wavefn%n_vir*batch_length,wavefn%n_J)
         L_ab_J=zero
!
!        Read Cholesky vectors
!
!         call get_cholesky_ab(L_ab_J,batch_start,batch_end,(wavefn%n_vir)*batch_length,.false.)
!
!        g_ab_ij=sum_J L_ab_J* L_ij_J
!
         g_off = index_two(1,batch_start,wavefn%n_vir)
!
         call dgemm('N','T',                      &
                     (wavefn%n_vir)*batch_length, &
                     n_oo,                        &
                     wavefn%n_J,                  &
                     one,                         &
                     L_ab_J,                      &
                     (wavefn%n_vir)*batch_length, &
                     L_ij_J,                      &
                     n_oo,                        &
                     one,                         &
                     g_ab_ij(g_off,1),            &
                     n_vv)
!
!
!        Deallocation of L_ab_J
!
         call deallocator(L_ab_J,batch_length*(wavefn%n_vir),wavefn%n_J)
!
      enddo ! batching done
!
!     Deallocation of L_ij_J
!
      call deallocator(L_ij_J,n_oo,wavefn%n_J)
!
!     Allocate for g_ai_jb
!
      call allocator(g_ai_jb,n_ov,n_ov)
      call allocator(L_ai_J,n_ov,wavefn%n_J)
      call allocator(L_ia_J,n_ov,wavefn%n_J)
!
!
!     Reading Cholesky vector L_ia_J and L_ai_J
!
!      call get_cholesky_ai(L_ai_J)
!
      call dgemm('N','T',     &
                  n_ov,       &
                  n_ov,       &
                  wavefn%n_J, &
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
     call deallocator(L_ia_J,n_ov,wavefn%n_J)
     call deallocator(L_ai_J,n_ov,wavefn%n_J)
!
!     Calculation of two-electron terms for virtual-virtual blocks
!
      do a = 1,wavefn%n_vir
         do b = 1,wavefn%n_vir
            ab=index_two(a,b, wavefn%n_vir)
            do i = 1,wavefn%n_occ
               ii=index_two(i,i, wavefn%n_occ)
               ai=index_two(a,i, wavefn%n_vir)
               bi=index_two(b,i, wavefn%n_vir)
               ia=index_two(i,a, wavefn%n_occ)
               ib=index_two(i,b, wavefn%n_occ)
               fock_matrix(wavefn%n_occ+a,wavefn%n_occ+b)=fock_matrix(wavefn%n_occ+a,wavefn%n_occ+b) &
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
      do i = 1,wavefn%n_occ
         do j = 1,wavefn%n_occ
            wavefn%fock_matrix_ij(i,j) = fock_matrix(i,j)
         enddo
      enddo
!
      do i = 1,wavefn%n_occ
         do a = 1,wavefn%n_vir
            wavefn%fock_matrix_ia(i,a) = fock_matrix(i,wavefn%n_occ+a)
            wavefn%fock_matrix_ai(a,i) = fock_matrix(wavefn%n_occ+a,i)
         enddo
      enddo
!
      do a = 1,wavefn%n_vir
         do b = 1,wavefn%n_vir
            wavefn%fock_matrix_ab(a,b) = fock_matrix(wavefn%n_occ+a,wavefn%n_occ+b)
         enddo
      enddo
!
   call deallocator(fock_matrix,wavefn%n_mo,wavefn%n_mo)
!
   end subroutine fock_constructor_cc_singles
!
   subroutine h1mo_T1(h1,h1_T1,n_vir,n_occ,n_orbitals,t1am,n_t1am)
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
   double precision h1(n_orbitals,n_orbitals)
   double precision h1_T1(n_orbitals,n_orbitals)
   double precision t1am(n_t1am,1)
!
   real(dp),dimension(:,:),allocatable     :: x 
   real(dp),dimension(:,:),allocatable     :: y 
   real(dp),dimension(:,:),allocatable     :: t1
   real(dp),dimension(:,:),allocatable     :: Z ! Intermediate for matrix multiplication
   integer(i15)                            :: p=0,q=0,a=0,i=0
   integer(i15)                            :: n_vir,n_occ,n_t1am,n_orbitals
!
!  Create t1, x, and y
!
   call allocator(t1,n_orbitals,n_orbitals)
   call allocator(y,n_orbitals,n_orbitals)
   call allocator(x,n_orbitals,n_orbitals)
   t1 = zero
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
   call deallocator(t1,n_orbitals,n_orbitals)
!
!  Allocate Intermediate
!   
   call allocator(Z,n_orbitals,n_orbitals)
!
!  h1_T1 = x*h1*y^T = x*Z
!
   call dgemm('N','T',n_orbitals,n_orbitals,n_orbitals &
      ,one,h1,n_orbitals,y,n_orbitals &
      ,zero,Z,n_orbitals)
!
   call dgemm('N','N',n_orbitals,n_orbitals,n_orbitals &
      ,one,x,n_orbitals,Z,n_orbitals &
      ,zero,h1_T1,n_orbitals)
!
!  Deallocations
!
   call deallocator(Z,n_orbitals,n_orbitals)
   call deallocator(y,n_orbitals,n_orbitals)
   call deallocator(x,n_orbitals,n_orbitals)
!
   end subroutine h1mo_T1
!
end module ccs_class
