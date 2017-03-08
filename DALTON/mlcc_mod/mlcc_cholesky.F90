module mlcc_cholesky
!
    use mlcc_data
    use mlcc_utilities
    use mlcc_types
    use mlcc_workspace
!
contains
   subroutine read_cholesky_ia(L_ia_J)
!
!     Purpose: Read Cholesky vectors L_ia^J from file and place them 
!              in the incoming vector  
!
      implicit none
!
      double precision L_ia_J(n_ov,n_J)
!
      integer :: lucho_ia 
      integer :: i,j,idummy
!
      lucho_ia = -1
      call gpopen(lucho_ia,'CHOLESKY_IA','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
      rewind(lucho_ia)
!
      do j=1,n_J
         read(lucho_ia) (L_ia_J(i,j), i=1,n_ov)
      enddo
!
      call gpclose(lucho_ia,'KEEP')      
!
   end subroutine read_cholesky_ia
!
   subroutine read_cholesky_ai(L_ai_J)
!
!     Purpose: Read Cholesky vectors L_ia^J from file and place them 
!              in the incoming vector  
!
      use mlcc_workspace
      implicit none
!
      double precision L_ai_J(n_ov,n_J)
!
      integer :: lucho_ia 
      integer :: i,a,ia,ai,j
      real(dp),dimension(:,:),pointer :: L_ia_J
!
!     Allocation
!
      call allocator(L_ia_J,n_ov,n_J)
      L_ia_J=zero
!     
!     IO
!
      call read_cholesky_ia(L_ia_J)
!
!     Reorder
!      
      do i = 1,n_occ
         do a = 1,n_vir
!
!           Needed indices
!
            ai=index_two(a,i,n_vir)
            ia=index_two(i,a,n_occ)
!
            do j=1,n_J
               L_ai_J(ai,j)=L_ia_J(ia,j)
            enddo
         enddo
      enddo
!
!     Deallocation
!
      call deallocator(L_ia_J,n_ov,n_J)

   end subroutine read_cholesky_ai
!
   subroutine read_cholesky_ij(L_ij_J)
!
!     Purpose: Read Cholesky vectors L_ij^J from file and place them 
!              in the incoming vector  
!
      implicit none
!
      double precision L_ij_J(n_oo,n_J)
!
      integer :: lucho_ij
      integer :: i,j,idummy
!
      lucho_ij = -1
      call gpopen(lucho_ij,'CHOLESKY_IJ','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
      rewind(lucho_ij)
!
      do j = 1,n_J
         read(lucho_ij) (L_ij_J(i,j), i=1,n_oo)
      enddo
!
      call gpclose(lucho_ij,'KEEP')    
!      
   end subroutine read_cholesky_ij
!
   subroutine read_cholesky_ab(L_ab_J,b_start,b_end,ab_dim)
!
!  Purpose: Read Cholesky vectors L_ab^J from file and place them 
!           in the incoming vector. If b_start .ne. 1 and b_end .eq. n_vir
!           we batch
! 
!           b_start - first element to be read
!           b_end   - last element to be read
!
!           ab_dim  - dimension over batching variables 
   implicit none
!
   integer :: lucho_ab,ab_dim
   integer :: a,b,j,idummy,i
   integer :: b_start,b_end
   real(dp),dimension(ab_dim,n_J) :: L_ab_J
   real(dp) :: dummy
   integer  :: batch_length 
!
   batch_length = b_end-b_start+1
!
   lucho_ab = -1
   call gpopen(lucho_ab,'CHOLESKY_AB','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
   rewind(lucho_ab)
!
   if (b_start .ne. 1) then
!
!     Calculate index of last element to throw away
!
      idummy=index_two(n_vir,b_start-1,n_vir)
!
!     Read from a_start
!
      do j = 1,n_J
         read(lucho_ab)(dummy,i=1,idummy),((L_ab_J(index_two(a,b,n_vir),j),b=1,batch_length),a=1,n_vir)
 !       read(lucho_ab)(dummy,i=1,idummy),(L_ab_J(a,j),a=1,ab_dim)
      enddo
   else
!
!     Read from start
!
      do j = 1,n_J
         read(lucho_ab)((L_ab_J(index_two(a,b,n_vir),j),b=1,batch_length),a=1,n_vir)
 !       read(lucho_ab)(L_ab_J(a,j),a=1,ab_dim)
      enddo
   endif
   !
   call gpclose(lucho_ab,'KEEP')    
!   
   end subroutine read_cholesky_ab
!
   subroutine read_cholesky_ab_reorder(L_ba_J,a_start,a_end,ab_dim)
!
!  Purpose: Read Cholesky vectors L_ab^J from file and place them 
!           in the incoming vector L_ba^J (note the different order,
!           used so that batching over a provides a single block matrix). 
!           If a_start .ne. 1 and a_end .eq. n_vir, we batch
! 
!           a_start - first element to be read
!           a_end   - last element to be read
!
!           ab_dim  - dimension over batching variables 
!
   implicit none
!
   integer :: lucho_ab,ab_dim
   integer :: a,b,j,idummy,i,k
   integer :: a_start,a_end
   real(dp),dimension(ab_dim,n_J) :: L_ba_J
   real(dp) :: dummy
   integer  :: batch_length 
!
   batch_length = a_end-a_start+1
!
   lucho_ab = -1
   call gpopen(lucho_ab,'CHOLESKY_AB','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
   rewind(lucho_ab)
!
!  Calculate number of elements to throw away
!
   idummy = index_two(n_vir,a_start-1,n_vir)
!
!  Loop over all Cholesky vectors
!
   do j = 1,n_J
!
!     Read in L_ba_J (which is equal to L_ab_J by symmetry)
!
      if (a_start .eq. 1) then 
         read(lucho_ab)((L_ba_J(index_two(b,a,n_vir),j),b=1,n_vir),a=1,batch_length)
      else
         read(lucho_ab)(dummy,i=1,idummy),((L_ba_J(index_two(b,a,n_vir),j),b=1,n_vir),a=1,batch_length)
      endif
!
   enddo
!
   call gpclose(lucho_ab,'KEEP')    
!   
   end subroutine read_cholesky_ab_reorder
!
   subroutine get_cholesky_ia(L_ia_J)
!
!  Purpose: Read and T1-transform ia cholesky vectors
!           L_ia_J_T1 = L_ia_J => only reading
!
      implicit none
!
      double precision L_ia_J(n_ov,n_J)
!
      call read_cholesky_ia(L_ia_J)
!
  end subroutine get_cholesky_ia
!
   subroutine get_cholesky_ai(L_ai_J)
!
!  Purpose: Read and T1-transform ia cholesky vectors
!           L_ai_J_T1 = L_ia_J - sum_j t_aj*L_ji_J + sum_b t_bi*L_ab_J
!                       - sum_(bj)t_aj*t_bi*L_jb_J
!
      implicit none
!
      double precision L_ai_J(n_ov,n_J)
!
      integer                          :: required,available,max_batch_length,n_batch
      integer                          :: a_batch,batch_start,batch_end,batch_length
      integer                          :: a,b,J,i,ai,Ja,ba,k,ik,iJ,kb,kJ
      real(dp),dimension(:,:),pointer  :: L_ba_J => null()
      real(dp),dimension(:,:),pointer  :: L_Ja_b => null()
      real(dp),dimension(:,:),pointer  :: L_Ja_i => null()
      real(dp),dimension(:,:),pointer  :: L_ik_J => null()
      real(dp),dimension(:,:),pointer  :: L_k_iJ => null()
      real(dp),dimension(:,:),pointer  :: L_a_iJ => null()
      real(dp),dimension(:,:),pointer  :: L_kJ_b => null()
      real(dp),dimension(:,:),pointer  :: L_kJ_i => null()
      real(dp),dimension(:,:),pointer  :: L_kb_J => null()
!
      call read_cholesky_ai(L_ai_J)
!
!
!!!     L_ab_J contributions    !!!
!
!
!
!     Allocate L_Ja_i
!
      call allocator(L_Ja_i,n_J*n_vir,n_occ)
      L_Ja_i=zero
!
!     Batching setup
!
      required = 2*n_vv*n_J*4
      available = get_available()
      max_batch_length = 0
      n_batch = 0
      a_batch = 0
      batch_length = 0
      batch_start = 0
      batch_end = 0 
!
      call n_one_batch(required,available,max_batch_length,n_batch,n_vir)
!
      do a_batch=1,n_batch
!
!        Allocate L_ab_J and L_Ja_b
!
         call allocator(L_ba_J,n_vir*batch_length,n_J)
         call allocator(L_Ja_b,n_vir*n_J,batch_length)
         L_ba_J=zero
         L_Ja_b=zero
!
!        Get start, end and length of batch
!
         call one_batch_limits(batch_start,batch_end,a_batch,max_batch_length,n_vir)
         batch_length=batch_end-batch_start+1
!
!        Read ab cholesky vectors, batching over a
!  
         call read_cholesky_ab_reorder(L_ba_J,batch_start,batch_end,n_vir*batch_length)
!
!        Reorder
!
         do a = 1,batch_length
            do b = 1,n_vir
               do J = 1,n_J
!
!                 Needed indices
!
                  ba=index_two(b,a,n_vir)
                  Ja=index_two(J,a,n_J)
!
                  L_Ja_b(Ja,b)=L_ba_J(ba,J)
!
               enddo
            enddo
         enddo
!
!        Add to L_Ja_i
!
         call dgemm('N','N',batch_length*n_J,n_occ,n_vir &
            ,one,L_Ja_b,batch_length*n_J,t1am,n_vir &
            , one, L_Ja_i(index_two(1,batch_start,n_J),1),n_vir*n_J)
!
!        Deallocate L_ab_J and L_Ja_b
!
         call deallocator(L_ba_J,n_vir*batch_length,n_J)
         call deallocator(L_Ja_b,n_vir*n_J,batch_length)
!
      enddo
!
!     Add terms to T1-transformed ai cholesky vector 
!
      do i = 1,n_occ
         do a = 1,n_vir
            do J=1,n_J
!
!              Needed indices
!
               Ja=index_two(J,a,n_J)
               ai=index_two(a,i,n_J)
!
               L_ai_J(ai,J)=L_ai_J(ai,J)+L_Ja_i(Ja,i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_Ja_i
!
      call deallocator(L_Ja_i,n_J*n_vir,n_occ)
!
!
!!!     L_ij_J contributions    !!!
!
!  
!
!     Allocate L_a_iJ, L_ik_J, L_k_iJ
!
      call allocator(L_a_iJ,n_vir,n_J*n_occ)
      call allocator(L_ik_J,n_oo,n_J)
      call allocator(L_k_iJ,n_occ,n_occ*n_J)
      L_a_iJ = zero   
      L_ik_J = zero 
      L_k_iJ = zero
!  
!     Read ij cholesky vectors
!
      call read_cholesky_ij(L_ik_J)
!
!     Reorder ij cholesky vectors
!
      do i = 1,n_occ
         do k=1,n_occ
            do J=1,n_J
!
!              Needed indices
!
               ik=index_two(i,k,n_occ)
               iJ=index_two(i,J,n_occ)
!
               L_k_iJ(k,iJ)=L_ik_J(ik,J)
!
            enddo
         enddo
      enddo
      call dgemm('N','N',n_vir,n_occ*n_J,n_occ &
         ,-one,t1am,n_vir,L_k_iJ,n_occ &
         ,zero,L_a_iJ,n_vir)
!
!     Add terms to T1-transformation of L_ai_J
!
      do i = 1,n_occ
         do a = 1,n_vir
            do J = 1,n_J
!
!              Needed indices
!
               ai=index_two(a,i,n_vir)
               iJ=index_two(i,J,n_occ)
!
               L_ai_J(ai,J)=L_ai_J(ai,J)+L_a_iJ(a,iJ)
            enddo
         enddo
      enddo
!
!     Deallocate L_a_iJ, L_ik_J, L_k_iJ
!
      call deallocator(L_a_iJ,n_vir,n_J*n_occ)      
      call deallocator(L_ik_J,n_oo,n_J)
      call deallocator(L_k_iJ,n_occ,n_occ*n_J)
!
!
!!!    L_jb_J contributions    !!!     
!
!
      call allocator(L_kJ_b,n_occ*n_J,n_vir)
      call allocator(L_kb_J,n_occ*n_vir,n_J)
!
!     Read L_kb_J
!
      call read_cholesky_ia(L_kb_J)
!
!     Reorder L_kb_J L_kJ_b
!
      do k = 1,n_occ
         do b=1,n_vir
            do J=1,n_J
!
!              Needed indices
!
               kb=index_two(k,b,n_occ)
               kJ=index_two(k,J,n_occ)
!
               L_kJ_b(kJ,b)=L_kb_J(kb,J)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_kb_J
!
      call deallocator(L_kb_J,n_occ*n_vir,n_J)
!
!     Allocate L_kJ_i for dgemm
!
      call allocator(L_kJ_i,n_occ*n_J,n_occ)
!
!     sum_b L_kJ_b*t_b_i = L_kJ_i
!
      call dgemm('N','N',n_occ*n_J,n_occ,n_vir &
         ,one,L_kJ_b,n_occ*n_J,t1am,n_vir &
         ,zero,L_kJ_i,n_occ*n_J)
!
!     Deallocate L_kJ_b
!
      call deallocator(L_kJ_b,n_occ*n_J,n_vir)
!
!     Allocate L_k_iJ
!  
      call allocator(L_k_iJ,n_occ,n_occ*n_J)
!
!     Reorder L_kJ_i to L_k_iJ    
!
      do i = 1,n_occ
         do k = 1,n_occ
            do J = 1,n_J
!
!              Needed indices
!
               kJ=index_two(k,J,n_occ)
               iJ=index_two(i,J,n_occ)
!
               L_k_iJ(k,iJ)=L_kJ_i(kJ,i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_kJ_i
!
      call deallocator(L_kJ_i,n_occ*n_J,n_occ)
!
!     Allocate L_a_iJ for dgemm
!
      call allocator(L_a_iJ,n_vir,n_occ*n_J)
!      
!     sum_k t_a_k*L_k_iJ = L_a_iJ
!
      call dgemm('N','N',n_vir,n_occ*n_J,n_occ &
         ,-one,t1am,n_vir,L_k_iJ,n_occ &
         ,zero,L_a_iJ,n_vir)
!
!     Add contribution to L ai_J
!
      do a = 1,n_vir
         do i = 1,n_occ
            do J = 1,n_J
!
!              Needed indices
!
               iJ=index_two(i,J,n_occ)
               ai=index_two(a,i,n_vir)
!
               L_ai_J(ai,J)=L_ai_J(ai,J)+L_a_iJ(a,iJ)
            enddo
         enddo
      enddo

      call deallocator(L_a_iJ,n_vir,n_occ*n_J)
      call deallocator(L_k_iJ,n_occ,n_occ*n_J)
!
  end subroutine get_cholesky_ai
!
   subroutine get_cholesky_ij(L_ij_J)
!
!  Purpose: Read and T1-transform ia cholesky vectors
!           L_ij_J_T1 = L_ij_J + sum_a t_aj*L_ia_J
!
!
!  Needed memory for routine: (n_occ*n_vir*n_J*2)
!
      implicit none
!
      double precision L_ij_J(n_oo,n_J)
!
      real(dp),dimension(:,:),pointer     :: L_ia_J => null()
      real(dp),dimension(:,:),pointer     :: L_iJ_a => null()
      real(dp),dimension(:,:),pointer     :: L_iJ_k => null()
      integer                             :: i,J,a,ij,ia,ik,k
!
!     Allocation
!    
      call allocator(L_ia_J,n_ov,n_J)
      call allocator(L_iJ_a,n_occ*n_J,n_vir)
!
!     Reading needed cholesky vectors
!
      call read_cholesky_ij(L_ij_J)
      call read_cholesky_ia(L_ia_J)
!
!     Reorder L_ia_J to L_iJ_a
!
      do i=1,n_occ
         do J=1,n_J
            do a=1,n_vir
!              
!              Needed indices
!
               iJ=index_two(i,J,n_occ)
               ia=index_two(i,a,n_occ)
!
               L_iJ_a(iJ,a)=L_ia_J(ia,J)
            enddo
         enddo
      enddo
!
!     Deallocate L_ia_J
!
      call deallocator(L_ia_J,n_ov,n_J)
!
!     Allocate L_iJ_k
!
      call allocator(L_iJ_k,n_occ*n_J,n_occ)
!
!     T1-transformation
!
      call dgemm('N','N',n_occ*n_J,n_occ,n_vir &
         , one,L_iJ_a,n_occ*n_J,t1am,n_vir &
         , one,L_iJ_k,n_occ*n_J)
!
!     Place terms from L_iJ_k into L_ij_J
!
      do i=1,n_occ
         do k=1,n_occ
            do J=1,n_J
!              
!              Needed indices
!
               iJ=index_two(i,J,n_occ)
               ik=index_two(i,k,n_occ)
!
               L_ij_J(ik,J)=L_ij_J(ik,J)+L_iJ_k(iJ,k)
            enddo
         enddo
      enddo
!
!     Deallocate L_iJ_k and L_iJ_a
!
      call deallocator(L_iJ_k,n_occ*n_J,n_occ)
      call deallocator(L_iJ_a,n_occ*n_J,n_vir)
!     
  end subroutine get_cholesky_ij
!
   subroutine get_cholesky_ab(L_ab_J,start,end,ab_dim,reorder)
!
!  Purpose: Read and T1-transform ia cholesky vectors
!           L_ab_J_T1 = L_ab_J - sum_i t_ai*L_ib_J
!
!   reorder = .true. => batch over first index
!   Required memory: n_J*batch_length*n_vir + n_vir*n_occ*n_J*2
      implicit none
!
      integer :: lucho_ab,ab_dim
      integer :: a,b,J,i,ia,aJ,ib,Jb,ab
      integer :: start,end
      logical,intent(in) :: reorder
!
      real(dp),dimension(:,:),pointer     :: L_ib_J => null()
      real(dp),dimension(:,:),pointer     :: L_Jb_i => null()
      real(dp),dimension(:,:),pointer     :: L_Jb_a => null()
      real(dp),dimension(:,:),pointer     :: L_a_Jb => null()
      real(dp),dimension(:,:),pointer     :: L_i_Jb => null()
      integer                             :: batch_length
      double precision L_ab_J(ab_dim,n_J)
      batch_length = end-start+1
!
!     Testing which index is batched
!
      if (reorder) then !! Batching over a !!
!
!        Allocate L_ib_J
!     
         call allocator(L_ib_J,n_ov,n_J)
         L_ib_J=zero
!
!        Read L_ia_J
!  
         call read_cholesky_ia(L_ib_J) ! OBS: Using L_ia_J insted of L_ai_J to avoid two reorderings!
                                    ! Possible because of symmetry L_ai_J(ai,J)==L_ia_J(ia,J)
!
!        Read L_ab_J for batch of a
!
         call read_cholesky_ab_reorder(L_ab_J,start,end,ab_dim)
!
!        Allocate L_i,Jb
!
         call allocator(L_i_Jb,n_occ,n_J*n_vir)
         L_i_Jb=zero
!
!        Reorder L_ib_J to L_i_Jb
!
         do i=1,n_occ
            do b=1,n_vir
               do J=1,n_J
!
!                 Needed indices
!
                  ib=index_two(i,b,n_occ)
                  Jb=index_two(J,b,n_J)
!
                  L_i_Jb(i,Jb)=L_ib_J(ib,J)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator(L_ib_J,n_ov,n_J)
!
!        Allocate L_a_Jb for batch of a
!  
         call allocator(L_a_Jb,batch_length,n_vir*n_J)
!
!        T1-transformation
!
         call dgemm('N','N',batch_length,n_vir*n_J,n_occ &
            ,-one,t1am(start,1),n_vir,L_i_Jb &
            ,zero, L_a_Jb,batch_length)
!
!        Add terms of L_a_Jb to L_ab_J
!
         do a=1,batch_length
            do b=1,n_vir
               do J=1,n_J
!
!                 Needed indices
!
                  Jb=index_two(J,b,n_J)
                  ab=index_two(a,b,n_vir)
!
                  L_ab_J(ab,J)=L_ab_J(ab,J)+L_a_Jb(a,Jb)
               enddo
            enddo
         enddo
!
!        Dellocate L_i_Jb and L_a_Jb for batch of a
!
         call deallocator(L_a_Jb,batch_length,n_J*n_vir)
         call deallocator(L_i_Jb,n_J*n_vir,n_occ)
!
      else  !! Batching over b !!
!
!        Allocate L_ib_J
!     
         call allocator(L_ib_J,n_ov,n_J)
         L_ib_J=zero
!
!        Read L_ia_J
!  
         call read_cholesky_ia(L_ib_J) ! OBS: Using L_ia_J insted of L_ai_J to avoid two reorderings!
                                    ! Possible because of symmetry L_ai_J(ai,J)==L_ia_J(ia,J)
!
!        Read L_ab_J for batch of b
!
         call read_cholesky_ab(L_ab_J,start,end,ab_dim)
!
!        Allocate L_Jb,i for batch of b
!
         call allocator(L_Jb_i,n_J*batch_length,n_occ)
         L_Jb_i=zero
!
!        Reorder L_ib_J to L_Jb_i
!
         do i=1,n_occ
            do b=1,batch_length
               do J=1,n_J
!
!                 Needed indices
!
                  ib=index_two(i,b+start-1,n_occ) ! OBS: in L_ib_J we have all b, not the case for L_Jb_i
                  Jb=index_two(J,b,n_J)
!
                  L_Jb_i(Jb,i)=L_ib_J(ib,J)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator(L_ib_J,n_ov,n_J)
!
!        Allocate L_Jb_a for batch of b
!  
         call allocator(L_Jb_a,n_J*batch_length,n_vir)
!
!        T1-transformation
!
         call dgemm('N','T',n_J*batch_length,n_vir,n_occ &
            ,-one,L_Jb_i,n_J*batch_length,t1am,n_vir &
            ,zero,L_Jb_a,batch_length*n_J)
!
!        Add terms of L_Jb_a to L_ab_J
!
         do a=1,n_vir
            do b=1,batch_length
               do J=1,n_J
!
!                 Needed indices
!
                  Jb=index_two(J,b,n_J)
                  ab=index_two(a,b,n_vir)
!
                  L_ab_J(ab,J)=L_ab_J(ab,J)+L_Jb_a(Jb,a)
               enddo
            enddo
         enddo
!
!        Dellocate L_Jb,i and L_Jb_a for batch of b
!
         call deallocator(L_Jb_a,n_J*batch_length,n_vir)
         call deallocator(L_Jb_i,n_J*batch_length,n_occ)
!
      endif
!
  end subroutine get_cholesky_ab
end module mlcc_cholesky