module mlcc_fock
!
use mlcc_types
use mlcc_data
use mlcc_workspace
use mlcc_utilities
!
contains
!
   subroutine mlcc_get_fock()
!  Purpose: Construct Fock matrix in MO basis
!
!  One-electron part obtained by reading one-electron integrals from MLCC_AOINT.
!  MLCC_AOINT is written to file from the mlcc_wr_aoint routine in sirius/mlcc_wr_aoint.F
!
!  Two-electron contributions calculated using Cholesky vectors
      use mlcc_cholesky
!
      implicit none
      real(dp), dimension(:,:), allocatable      :: fock_ao
      real(dp), dimension(:,:), allocatable      :: ao_int 
      real(dp), dimension(:,:), allocatable      :: h1mo   
      real(dp), dimension(:,:), allocatable      :: X      
      integer                                    :: luaoin = -1
      integer                                    :: idummy=0,i=0,j=0,k=0,a=0,b=0,ij=0,kk=0,ik=0,kj=0,ii=0,jj=0
      integer                                    :: ji=0,ai=0,ib=0,bi=0,ia=0,aj=0,ja=0,ab=0
      real(dp), dimension(:,:), allocatable      :: g_ij_kl
      real(dp), dimension(:,:), allocatable      :: g_ab_ij
      real(dp), dimension(:,:), allocatable      :: g_ai_jb
      real(dp), dimension(:,:), allocatable      :: g_ia_jk
      real(dp), dimension(:,:), allocatable      :: g_ai_jk
      real(dp), dimension(:,:), allocatable      :: L_ij_J 
      real(dp), dimension(:,:), allocatable      :: L_ia_J 
      real(dp), dimension(:,:), allocatable      :: L_ai_J 
      real(dp), dimension(:,:), allocatable      :: L_ab_J 
      integer                                    :: available=0,required=0,max_batch_length=0,n_batch=0,batch_start=0
      integer                                    :: batch_end=0,batch_length=0,g_off=0
      integer                                    :: b_batch = 0
      logical                                    :: debug = .true.
!
   call allocator_n(mo_fock_mat,n_orbitals,n_orbitals)
   call allocator_n(h1mo,n_orbitals,n_orbitals)
   mo_fock_mat = zero
   h1mo = zero
!
!
!
!!! ONE-ELECTRON CONTRIBUTION !!!
!
!
!  Allocate for one-electron ao integrals
!
      call allocator_n(ao_int,n_basis_2_pack,1)
!
!  Read one-electron ao integrals
!
      call gpopen(luaoin,'MLCC_AOINT','UNKNOWN','SEQUENTIAL','FORMATTED',idummy,.false.)
      rewind(luaoin)
!
      read(luaoin,*)(ao_int(i,1),i=1,n_basis_2_pack)
!
      call gpclose(luaoin,'KEEP')
!
!  Allocate ao fock matrix, add one-electron contributions
!
      call allocator_n(fock_ao,n_basis,n_basis)
      call squareup(ao_int,fock_ao,n_basis)   
!
!  Deallocation of ao integrals
!   
      call deallocator_n(ao_int,n_basis_2_pack,1) 
!
!  Transform to one-electron part to mo and add to mo_fock_mat 
!
      call allocator_n(X,n_basis,n_orbitals)
!
      call dgemm('N','N',n_basis,n_orbitals,n_basis,one,fock_ao,n_basis,orb_coefficients,n_basis,zero,X,n_basis)
      call dgemm('T','N',n_orbitals,n_orbitals,n_basis,one,orb_coefficients,n_basis,X,n_basis,zero,h1mo,n_orbitals)
! 
!     T1-transformation of one-electron integrals in mo basis
!
      call h1mo_T1(h1mo,mo_fock_mat)
      call deallocator_n(h1mo,n_orbitals,n_orbitals)
!
      call deallocator_n(X,n_basis,n_orbitals)
      call deallocator_n(fock_ao,n_basis,n_basis)
!
!
!!! TWO-ELECTRON CONTRIBUTION !!!
!
!
!
!!  occupied-occupied block: F_ij = h_ij + sum_k (2*g_ijkk - g_ikkj) !!
!  
!
!  Allocation for L_ij_J_pack
! 
   call allocator_n(L_ij_J,n_oo,n_J)
   call allocator_n(g_ij_kl,n_oo,n_oo)
   L_ij_J=zero
   g_ij_kl=zero
!
!  Read cholesky vectors
!
   call get_cholesky_ij(L_ij_J)
!
!  g_ij_kl = sum_J(L_ij_J*L_J_kl)
!
   call dgemm('N','T',n_oo,n_oo,n_J &
      ,one,L_ij_J,n_oo,L_ij_J,n_oo &
      ,zero,g_ij_kl,n_oo)
!
!  Add two-electron contributions to occupied-occupied block
!
   do i=1,n_occ
      do j=1,n_occ
         ij=index_two(i,j,n_occ)
         do k = 1,n_occ
            kk=index_two(k,k,n_occ)
            ik=index_two(i,k,n_occ)
            kj=index_two(k,j,n_occ)
            mo_fock_mat(i,j)=mo_fock_mat(i,j)+two*g_ij_kl(ij,kk)-g_ij_kl(ik,kj)
         enddo
      enddo
   enddo
!
!     Deallocate g_ij_kl
! 
      call deallocator_n(g_ij_kl,n_oo,n_oo)
!
!!    Occupied-vacant blocks F_ai=F_ia=0 because this Fock matrix satisfies the HF equations !! OBS T1!
!
!
!     Allocation for g_ia_jk 
!
     call allocator_n(L_ia_J,n_ov,n_J)
     call allocator_n(g_ia_jk,n_ov,n_oo)
!
!     Reading Cholesky vector L_ia_J
!
      call get_cholesky_ia(L_ia_J)
!
!     g_ia_jk
!
      call dgemm('N','T',n_ov,n_oo,n_J,one,L_ia_J,n_ov,L_ij_J,n_oo,zero,g_ia_jk,n_ov)
      call deallocator_n(L_ia_J,n_ov,n_J)
!
!     Allocation for g_ai_jk 
!
      call allocator_n(L_ai_J,n_ov,n_J)
      call allocator_n(g_ai_jk,n_ov,n_oo)
!
!     Reading Cholesky vector L_ai_J
!
      call get_cholesky_ai(L_ai_J)
!
!     g_ai_jk
!
      call dgemm('N','T',n_ov,n_oo,n_J,one,L_ai_J,n_ov,L_ij_J,n_oo,zero,g_ai_jk,n_ov)
      call deallocator_n(L_ai_J,n_ov,n_J)
!
!     Adding terms to Fock matrix
!
      do i=1,n_occ
         do a=1,n_vir
            do j=1,n_occ
!
!              Needed indices
!
               ia=index_two(i,a,n_occ)
               ja=index_two(j,a,n_occ)
!
               ai=index_two(a,i,n_vir)
               aj=index_two(a,j,n_vir)
!
               jj=index_two(j,j,n_occ)
               ji=index_two(j,i,n_occ)
               ij=index_two(i,j,n_occ)
!
               mo_fock_mat(i,a+n_occ)=mo_fock_mat(i,a+n_occ)+two*g_ia_jk(ia,jj)-g_ia_jk(ja,ij)
               mo_fock_mat(a+n_occ,i)=mo_fock_mat(a+n_occ,i)+two*g_ai_jk(ai,jj)-g_ai_jk(aj,ji)
            enddo
         enddo
      enddo
      call deallocator_n(g_ia_jk,n_ov,n_oo)
      call deallocator_n(g_ai_jk,n_ov,n_oo)
!
!!    Vacant-vacant block F_ab = h_ab + sum_k (2*g_abkk - g_akkb) !!
!
      call allocator_n(g_ab_ij,n_vv,n_oo)
      g_ab_ij=zero
!
!     Batching over a
!
!
!     Setup of variables needed for batching
!
      available = get_available()
      required = 2*n_vir*n_vir*n_J*4 + 2*n_vir*n_occ*n_J*4
      call n_one_batch(required,available,max_batch_length,n_batch,n_vir)
!
      batch_start=1
      batch_end=0
      batch_length=0
!
!     Start batchig loop
!
      do b_batch = 1,n_batch
!
!        Get batch limits  and  length of batch
!
         call one_batch_limits(batch_start,batch_end,b_batch,max_batch_length,n_vir)
         batch_length=batch_end-batch_start+1
!
!        Allocation of L_ab_J
!
         call allocator_n(L_ab_J,n_vir*batch_length,n_J)
         L_ab_J=zero
!
!        Read Cholesky vectors
!
         call get_cholesky_ab(L_ab_J,batch_start,batch_end,n_vir*batch_length,.false.)
!
!        g_ab_ij=sum_J L_ab_J* L_ij_J
!
         g_off = index_two(1,batch_start,n_vir)
!
         call dgemm('N','T',n_vir*batch_length,n_oo,n_J &
            ,one,L_ab_J,n_vir*batch_length,L_ij_J,n_oo &
            ,one,g_ab_ij(g_off,1),n_vv)
!
!        Deallocation of L_ab_J
!
         call deallocator_n(L_ab_J,batch_length*n_vir,n_J)
!
      enddo ! batching done
!
!     Deallocation of L_ij_J
!
      call deallocator_n(L_ij_J,n_oo,n_J)
!
!     Allocate for g_ai_jb
!
      call allocator_n(g_ai_jb,n_ov,n_ov)
      call allocator_n(L_ai_J,n_ov,n_J)
      call allocator_n(L_ia_J,n_ov,n_J)
!
!
!     Reading Cholesky vector L_ia_J and L_ai_J
!
      call get_cholesky_ia(L_ia_J)
      call get_cholesky_ai(L_ai_J)
      call dgemm('N','T',n_ov,n_ov,n_J,one,L_ai_J,n_ov,L_ia_J,n_ov,zero,g_ai_jb,n_ov)
!
!     Deallocate L_ia_J
!
     call deallocator_n(L_ia_J,n_ov,n_J)
     call deallocator_n(L_ai_J,n_ov,n_J)
!
!     Calculation of two-electron terms for virtual-virtual blocks
!
      do a = 1,n_vir
         do b = 1,n_vir
            ab=index_two(a,b,n_vir)
            do i = 1,n_occ
               ii=index_two(i,i,n_occ)
               ai=index_two(a,i,n_vir)
               bi=index_two(b,i,n_vir)
               ia=index_two(i,a,n_occ)
               ib=index_two(i,b,n_occ)
               mo_fock_mat(n_occ+a,n_occ+b)=mo_fock_mat(n_occ+a,n_occ+b)+two*g_ab_ij(ab,ii)-g_ai_jb(ai,ib)
            enddo
         enddo 
      enddo
     call deallocator_n(g_ab_ij,n_vv,n_oo)
     call deallocator_n(g_ai_jb,n_ov,n_ov)
!
!    Clean-up of Fock matrix
!
     call mlcc_cleanup(mo_fock_mat,n_orbitals,n_orbitals)
!
!     Prints
!
      if (debug) then
         do i=1,n_occ
            write(luprint,*)(mo_fock_mat(i,j),j=1,n_occ)
         enddo
         do i=1,n_vir
            write(luprint,*)(mo_fock_mat(j,i+n_occ),j=1,n_occ)
         enddo
         do i=1,n_vir
            write(luprint,*)(mo_fock_mat(i+n_occ,j),j=1,n_occ)
         enddo
         do i=1,5
            write(luprint,*)(mo_fock_mat(i+n_occ,j+n_occ),j=1,5) 
         enddo
      endif
!
!     Save the blocks of the Fock matrix in memory (ij,ia,ai,ab)
!
      F_i_j = zero
      F_i_a = zero
      F_a_i = zero
      F_a_b = zero
!
      do i = 1,n_occ
         do j = 1,n_occ
            F_i_j(i,j) = mo_fock_mat(i,j)
         enddo
      enddo
!
      do i = 1,n_occ
         do a = 1,n_vir
            F_i_a(i,a) = mo_fock_mat(i,n_occ+a)
            F_a_i(a,i) = mo_fock_mat(n_occ+a,i)
         enddo
      enddo
!
      do a = 1,n_vir
         do b = 1,n_vir
            F_a_b(a,b) = mo_fock_mat(n_occ+a,n_occ+b)
         enddo
      enddo
!

   call deallocator_n(mo_fock_mat,n_orbitals,n_orbitals)
!
   end subroutine mlcc_get_fock
!
   subroutine h1mo_T1(h1,h1_T1)
!
!  Purpose: T1-transform of one-electron mo integrals (h1mo)
!
!           h_p_q_T1= sum_(st) x_p_s * y_q_t * h_s_t
!           x = I-t1
!           y = I-t1^T
!
!
   implicit none
!
   double precision h1(n_orbitals,n_orbitals)
   double precision h1_T1(n_orbitals,n_orbitals)
!
   real(dp),dimension(:,:),allocatable     :: x 
   real(dp),dimension(:,:),allocatable     :: y 
   real(dp),dimension(:,:),allocatable     :: t1
   real(dp),dimension(:,:),allocatable     :: Z ! Intermediate for matrix multiplication
   integer                                 :: p=0,q=0,a=0,i=0
!
!  Create t1, x, and y
!
   call allocator_n(t1,n_orbitals,n_orbitals)
   call allocator_n(y,n_orbitals,n_orbitals)
   call allocator_n(x,n_orbitals,n_orbitals)
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
   call deallocator_n(t1,n_orbitals,n_orbitals)
!
!  Allocate Intermediate
!   
   call allocator_n(Z,n_orbitals,n_orbitals)
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
   call deallocator_n(Z,n_orbitals,n_orbitals)
   call deallocator_n(y,n_orbitals,n_orbitals)
   call deallocator_n(x,n_orbitals,n_orbitals)
!
   end subroutine h1mo_T1
end module mlcc_fock