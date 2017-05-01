submodule (ccs_class) cholesky
!
!
!                       -::- Cholesky submodule (CCS) -::-
!           Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!
!     Contains the following family of procedures of the CCS class:
!
!        get_cholesky_ij(L_ij_J):         returns the T1-transformed Cholesky vector L_ij^J 
!        get_cholesky_ia(L_ia_J):         returns the T1-transformed Cholesky vector L_ia^J 
!        get_cholesky_ai(L_ai_J):         returns the T1-transformed Cholesky vector L_ai^J
!        get_cholesky_ab(L_ab_J, ...):    returns the T1-transformed Cholesky vector L_ab^J,
!                                         but has options (...) for batching over the two 
!                                         indices a and b
!
   implicit none 
!
!
contains
!
!
   module subroutine get_cholesky_ij_cc_singles(wf, L_ij_J)
!
!     Get Cholesky IJ
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Reads and T1-transforms the IA Cholesky vectors:
!     
!        L_ij_J_T1 = L_ij_J + sum_a t_aj * L_ia_J
!
!     Needed memory by routine: (n_o*n_v*n_J*2)
!
      implicit none 
!
      class(cc_singles) :: wf
!
      real(dp), dimension((wf%n_o)**2, wf%n_J) :: L_ij_J
!
      real(dp), dimension(:,:), allocatable :: L_ia_J ! L_ia^J
      real(dp), dimension(:,:), allocatable :: L_iJ_a ! L_ia^J reordered
      real(dp), dimension(:,:), allocatable :: L_iJ_k ! L_ik^J reordered
!
      integer(i15) :: i = 0, J = 0, a = 0, ij = 0, ia = 0, ik = 0, k = 0
!
!     Allocate
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call allocator(L_iJ_a, (wf%n_o)*(wf%n_J), wf%n_v)
!
      L_ia_J = zero
      L_iJ_a = zero
!
!     Read the untransformed Cholesky vectors 
!
      call wf%read_cholesky_ij(L_ij_J)
      call wf%read_cholesky_ia(L_ia_J)
!
!     Reorder L_ia_J to L_iJ_a 
!
      do i = 1, wf%n_o
         do J = 1, wf%n_J
            do a = 1, wf%n_v
!              
!              Needed indices
!
               iJ = index_two(i, J, wf%n_o)
               ia = index_two(i, a, wf%n_o)
!
               L_iJ_a(iJ, a) = L_ia_J(ia, J)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_ia_J
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate L_iJ_k
!
      call allocator(L_iJ_k, (wf%n_o)*(wf%n_J), wf%n_o)
      L_iJ_k = zero
!
!     T1-transformation
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_J), &
                  wf%n_o,            &
                  wf%n_v,            &
                  one,               &
                  L_iJ_a,            &
                  (wf%n_o)*(wf%n_J), &
                  wf%t1am,           &
                  wf%n_v,            &
                  zero,              &
                  L_iJ_k,            &
                  (wf%n_o)*(wf%n_J))
!
!     Place terms from L_iJ_k into L_ij_J
!
      do i = 1, wf%n_o
         do k = 1, wf%n_o
            do J = 1, wf%n_J
!              
!              Needed indices
!
               iJ = index_two(i, J, wf%n_o)
               ik = index_two(i, k, wf%n_o)
!
               L_ij_J(ik, J) = L_ij_J(ik, J) + L_iJ_k(iJ, k)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_iJ_k and L_iJ_a
!
      call deallocator(L_iJ_k, (wf%n_o)*(wf%n_J), wf%n_o)
      call deallocator(L_iJ_a, (wf%n_o)*(wf%n_J), wf%n_v)
!    
   end subroutine get_cholesky_ij_cc_singles
!
!
   module subroutine get_cholesky_ia_cc_singles(wf, L_ia_J)
!
!     Get Cholesky IA
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Reads and T1-transforms IA Cholesky vectors
!
!        L_ia_J_T1 = L_ia_J (only reading necessary)
!
      implicit none 
!
      class(cc_singles) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), wf%n_J) :: L_ia_J
!
      call wf%read_cholesky_ia(L_ia_J)
!
   end subroutine get_cholesky_ia_cc_singles
!
!
   module subroutine get_cholesky_ai_cc_singles(wf, L_ai_J)
!
!     Get Cholesky AI
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Read and T1-transform Cholesky AI vectors:
!     
!        L_ai_J_T1 = L_ia_J - sum_j  t_aj*L_ji_J 
!                           + sum_b  t_bi*L_ab_J
!                           - sum_bj t_aj*t_bi*L_jb_J
!
!     Required memory by routine: 2*n_v*n_J*batch_length 
!                                 + n_v*n_o*n_J
!
      implicit none 
!
      class(cc_singles) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), wf%n_J) :: L_ai_J 
!
      logical :: reorder ! Reorder or not, when reading Cholesky AB 
!
!     Batch variables
!
      integer(i15) :: required = 0, available = 0, max_batch_length = 0, n_batch = 0, L_off = 0
      integer(i15) :: a_batch = 0, batch_start = 0, batch_end = 0, batch_length = 0
!
!     Indices
!
      integer(i15) :: a = 0, b = 0, J = 0, i = 0, ai = 0, Ja = 0
      integer(i15) :: ba = 0, k = 0, ik = 0, iJ = 0, kb = 0, kJ = 0
!
!     Cholesky vectors (in many different orderings)
!
      real(dp), dimension(:,:), allocatable :: L_ba_J
      real(dp), dimension(:,:), allocatable :: L_Ja_b
      real(dp), dimension(:,:), allocatable :: L_Ja_i
      real(dp), dimension(:,:), allocatable :: L_ik_J
      real(dp), dimension(:,:), allocatable :: L_k_iJ
      real(dp), dimension(:,:), allocatable :: L_a_iJ
      real(dp), dimension(:,:), allocatable :: L_kJ_b
      real(dp), dimension(:,:), allocatable :: L_kJ_i
      real(dp), dimension(:,:), allocatable :: L_kb_J
!
!     Read L_ai^J from file 
!
      call wf%read_cholesky_ai(L_ai_J)
!
!                       
!     -::- L_ab_J contributions -::-
!
!
!     Allocate L_Ja_i
!
      call allocator(L_Ja_i, (wf%n_J)*(wf%n_v), wf%n_o)
      L_Ja_i = zero
!
!     Set batching variables 
!
      required = 2*(wf%n_v)**2*(wf%n_J)*4
      available = get_available()
      max_batch_length = 0
!
      n_batch = 0
      a_batch = 0
!
      batch_length = 0
      batch_start = 0
      batch_end = 0
!
!     Calculate the number of batches 
!
      call num_batch(required, available, max_batch_length, n_batch, wf%n_v)
!
      do a_batch = 1, n_batch
!
!        Get start, end, and length of batch
!
         call batch_limits(batch_start, batch_end, a_batch, max_batch_length, wf%n_v)
         batch_length = batch_end - batch_start + 1
!
!        Allocate L_ab_J and L_Ja_b
!
         call allocator(L_ba_J, (wf%n_v)*batch_length, wf%n_J) ! L_ab^J = L_ba_J(ba,J)
         call allocator(L_Ja_b, batch_length*(wf%n_J), wf%n_v)
!
         L_ba_J = zero
         L_Ja_b = zero 
!
!        Read Cholesky AB vectors, batching over a
! 
         reorder = .true.
         call wf%read_cholesky_ab(L_ba_J, batch_start, batch_end, (wf%n_v)*batch_length, reorder)
!
!        Reorder the Cholesky array L_ba_J
!
         do a = 1, batch_length
            do b = 1, wf%n_v
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  ba = index_two(b, a, wf%n_v)
                  Ja = index_two(J, a, wf%n_J)
!
                  L_Ja_b(Ja, b) = L_ba_J(ba, J) ! L_ab^J
!
               enddo
            enddo
         enddo        
!
!        Calculate sum_b L_Ja_b*t_b_i = L_Ja_i 
!        
         L_off = index_two(1, batch_start, wf%n_J)
!
         call dgemm('N','N',                &
                     batch_length*(wf%n_J), &
                     wf%n_o,                &
                     wf%n_v,                &
                     one,                   &
                     L_Ja_b,                &
                     batch_length*(wf%n_J), &
                     wf%t1am,               &
                     wf%n_v,                &
                     one,                   &
                     L_Ja_i(L_off, 1),      &
                     (wf%n_v)*(wf%n_J))
!
!        Deallocate L_ab_J and L_Ja_b
!
         call deallocator(L_ba_J, (wf%n_v)*batch_length, wf%n_J)
         call deallocator(L_Ja_b, batch_length*(wf%n_J), wf%n_v)
!
      enddo ! batching over a 
!
!     Add terms to T1-transformed Cholesky AI vector 
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do J = 1, wf%n_J
!
!              Needed indices
!
               Ja = index_two(J, a, wf%n_J)
               ai = index_two(a, i, wf%n_v)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_Ja_i(Ja, i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_Ja_i
!
      call deallocator(L_Ja_i, (wf%n_J)*(wf%n_v), wf%n_o)
!
!
!     -::- L_ij_J contributions -::-
!
!
!     Allocate L_a_iJ, L_ik_J, L_k_iJ
!
      call allocator(L_a_iJ, wf%n_v, (wf%n_J)*(wf%n_o))
      call allocator(L_k_iJ, wf%n_o, (wf%n_o)*(wf%n_J))
!
      call allocator(L_ik_J, (wf%n_o)**2, wf%n_J)
!
      L_a_iJ = zero   
      L_k_iJ = zero
!
      L_ik_J = zero 
!  
!     Read Cholesky IJ vectors
!
      call wf%read_cholesky_ij(L_ik_J) ! L_ik_J(ik,J) = L_ik^J 
!
!     Reorder IJ Cholesky vectors
!
      do i = 1, wf%n_o
         do k = 1, wf%n_o
            do J = 1, wf%n_J
!
!              Needed indices
!
               ik = index_two(i, k, wf%n_o)
               iJ = index_two(i, J, wf%n_o)
!
               L_k_iJ(k, iJ) = L_ik_J(ik, J) ! L_k_iJ(k,iJ) = L_ik^J
!
            enddo
         enddo
      enddo
!
!     Calculate -sum_k t_a_k*L_k_iJ = L_a_iJ  ! Here we assume L_ik^J = L_ki^J 
!
      call dgemm('N','N',            &
                  wf%n_v,            &
                  (wf%n_o)*(wf%n_J), &
                  wf%n_o,            &
                  -one,              &
                  wf%t1am,           &
                  wf%n_v,            &
                  L_k_iJ,            &
                  wf%n_o,            &
                  zero,              &
                  L_a_iJ,            &
                  wf%n_v)
!
!     Add terms to T1-transformation of L_ai_J
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do J = 1, wf%n_J
!
!              Needed indices
!
               ai = index_two(a, i, wf%n_v)
               iJ = index_two(i, J, wf%n_o)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_a_iJ, L_ik_J, L_k_iJ
!
      call deallocator(L_a_iJ, wf%n_v, (wf%n_J)*(wf%n_o))      
      call deallocator(L_k_iJ, wf%n_o, (wf%n_o)*(wf%n_J))
!
      call deallocator(L_ik_J, (wf%n_o)**2, wf%n_J)
!
!
!     -::- L_jb_J contributions -::-    
!
!
      call allocator(L_kJ_b, (wf%n_o)*(wf%n_J), wf%n_v)
      call allocator(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Read the Cholesky vector L_kb_J
!
      call wf%read_cholesky_ia(L_kb_J)
!
!     Reorder L_kb_J to L_kJ_b
!
      do k = 1, wf%n_o
         do b = 1, wf%n_v
            do J = 1, wf%n_J
!
!              Needed indices
!
               kb = index_two(k, b, wf%n_o)
               kJ = index_two(k, J, wf%n_o)
!
               L_kJ_b(kJ, b) = L_kb_J(kb, J)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_kb_J
!
      call deallocator(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate L_kJ_i 
!
      call allocator(L_kJ_i, (wf%n_o)*(wf%n_J), wf%n_o)
      L_kJ_i = zero
!
!     Calculate sum_b L_kJ_b*t_b_i = L_kJ_i
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_J), &
                  wf%n_o,            &
                  wf%n_v,            &
                  one,               &
                  L_kJ_b,            &
                  (wf%n_o)*(wf%n_J), &
                  wf%t1am,           &
                  wf%n_v,            &
                  zero,              &
                  L_kJ_i,            &
                  (wf%n_o)*(wf%n_J))
!
!     Deallocate L_kJ_b
!
      call deallocator(L_kJ_b, (wf%n_o)*(wf%n_J), wf%n_v)
!
!     Allocate L_k_iJ
!  
      call allocator(L_k_iJ, (wf%n_o), (wf%n_o)*(wf%n_J))
      L_k_iJ = zero
!
!     Reorder L_kJ_i to L_k_iJ    
!
      do i = 1, wf%n_o
         do k = 1, wf%n_o
            do J = 1, wf%n_J
!
!              Needed indices
!
               kJ = index_two(k, J, wf%n_o)
               iJ = index_two(i, J, wf%n_o)
!
               L_k_iJ(k, iJ) = L_kJ_i(kJ, i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_kJ_i
!
      call deallocator(L_kJ_i, (wf%n_o)*(wf%n_J), wf%n_o)
!
!     Allocate L_a_iJ
!
      call allocator(L_a_iJ, wf%n_v, (wf%n_o)*(wf%n_J))
      L_a_iJ = zero
!      
!     Calculate sum_k t_a_k*L_k_iJ = L_a_iJ
!
      call dgemm('N','N',            &
                  wf%n_v,            &
                  (wf%n_o)*(wf%n_J), &
                  wf%n_o,            &
                  -one,              &
                  wf%t1am,           &
                  wf%n_v,            &
                  L_k_iJ,            &
                  wf%n_o,            &
                  zero,              &
                  L_a_iJ,            &
                  wf%n_v)
!
!     Add contribution to L ai_J
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do J = 1, wf%n_J
!
!              Needed indices
!
               iJ = index_two(i, J, wf%n_o)
               ai = index_two(a, i, wf%n_v)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
            enddo
         enddo
      enddo
!
!     Deallocations
!
      call deallocator(L_a_iJ, wf%n_v, (wf%n_o)*(wf%n_J))
      call deallocator(L_k_iJ, wf%n_o, (wf%n_o)*(wf%n_J))
!
   end subroutine get_cholesky_ai_cc_singles
!
!
   module subroutine get_cholesky_ab_cc_singles(wf, L_ab_J, first, last, ab_dim, reorder)
!
!     Get Cholesky AB
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Reads and T1-transforms the IA Cholesky vectors:
!
!           L_ab_J_T1 = L_ab_J - sum_i t_ai*L_ib_J
!
!     If reorder = .true.,  L_ba_J is returned with batching over a
!     If reorder = .false., L_ab_J is returned with batching over b
!
!     Required memory: n_J*batch_length*n_v + n_v*n_o*n_J*2
!
      implicit none
!
      class(cc_singles) :: wf
!
      integer(i15), intent(in) :: ab_dim
      integer(i15), intent(in) :: first
      integer(i15), intent(in) :: last
!
      logical, intent(in) :: reorder
!
      real(dp), dimension(ab_dim, wf%n_J) :: L_ab_J
!
      integer(i15) :: memory_lef = 0
!
      integer :: unit_chol_ab = -1 ! Unit identifier for cholesky_ab file
!
      integer :: a = 0, b = 0, J = 0, i = 0, ia = 0, aJ = 0, ib = 0, Jb = 0, ab = 0, ba = 0
!
      real(dp), dimension(:,:), allocatable :: L_ib_J
      real(dp), dimension(:,:), allocatable :: L_Jb_i
      real(dp), dimension(:,:), allocatable :: L_Jb_a
      real(dp), dimension(:,:), allocatable :: L_a_Jb
      real(dp), dimension(:,:), allocatable :: L_i_Jb

      integer(i15) :: batch_length = 0
!
      batch_length = last - first + 1
!
!     Testing which index is batched
!
      if (reorder) then !! Batching over a !!
!
!        Allocate L_ib_J
!     
         call allocator(L_ib_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ib_J = zero
!
!        Read L_ia_J
!  
         call wf%read_cholesky_ia(L_ib_J) ! Note: using L_ia_J instead of L_ai_J, here, to avoid two reorderings.
                                          ! This is possible because of the symmetry L_ai_J(ai,J) == L_ia_J(ia,J).
!
!        Read L_ab_J for batch of a
!
         call wf%read_cholesky_ab(L_ab_J, first, last, ab_dim, reorder) ! L_ab_J(ba) = L_ab^J 
!
!        Allocate L_i,Jb
!
         call allocator(L_i_Jb, wf%n_o, (wf%n_J)*(wf%n_v))
         L_i_Jb = zero
!
!        Reorder L_ib_J to L_i_Jb
!
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  ib = index_two(i, b, wf%n_o)
                  Jb = index_two(J, b, wf%n_J)
!
                  L_i_Jb(i, Jb) = L_ib_J(ib, J)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator(L_ib_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Allocate L_a_Jb for batch of a
!  
         call allocator(L_a_Jb, batch_length, (wf%n_v)*(wf%n_J))
!
!        Calculate  -t1_a_i * L_i_Jb = L_a_Jb
!
         call dgemm('N','N',            &
                     batch_length,      &
                     (wf%n_v)*(wf%n_J), &
                     wf%n_o,            &
                     -one,              &
                     wf%t1am(first,1),  &
                     wf%n_v,            &
                     L_i_Jb,            &
                     wf%n_o,            &
                     zero,              &
                     L_a_Jb,            &
                     batch_length)
!
!        Add terms of L_a_Jb to L_ab_J
!
         do a = 1, batch_length
            do b = 1, wf%n_v
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  Jb = index_two(J, b, wf%n_J)
                  ba = index_two(b, a, wf%n_v)
!
                  L_ab_J(ba, J) = L_ab_J(ba, J) + L_a_Jb(a, Jb) 
!
               enddo
            enddo
         enddo
!
!        Dellocate L_i_Jb and L_a_Jb for batch of a
!
         call deallocator(L_a_Jb, batch_length, (wf%n_J)*(wf%n_v))
         call deallocator(L_i_Jb, (wf%n_J)*(wf%n_v), (wf%n_o))
!
      else  !! Batching over b !!
!
!        Allocate L_ib_J
!     
         call allocator(L_ib_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ib_J = zero
!
!        Read L_ia_J
!  
         call wf%read_cholesky_ia(L_ib_J) ! Note: using L_ia_J instead of L_ai_J, here, to avoid two reorderings.
                                             ! This is possible because of the symmetry L_ai_J(ai,J) == L_ia_J(ia,J)
!
!        Read L_ab_J for batch of b
!
         call wf%read_cholesky_ab(L_ab_J, first, last, ab_dim, reorder)
!
!        Allocate L_Jb,i for batch of b
!
         call allocator(L_Jb_i, (wf%n_J)*batch_length, wf%n_o)
         L_Jb_i = zero
!
!        Reorder L_ib_J to L_Jb_i
!
         do i = 1, wf%n_o
            do b = 1, batch_length
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  ib = index_two(i, b + first - 1, wf%n_o) ! Note: in L_ib_J we have all b, not the case for L_Jb_i
                  Jb = index_two(J, b, wf%n_J)
!
                  L_Jb_i(Jb, i) = L_ib_J(ib, J)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator(L_ib_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Allocate L_Jb_a for batch of b
!  
         call allocator(L_Jb_a, (wf%n_J)*batch_length, wf%n_v)
!
!        T1-transformation
!
         call dgemm('N','T',                &
                     (wf%n_J)*batch_length, &
                     wf%n_v,                &
                     wf%n_o,                &
                     -one,                  &
                     L_Jb_i,                &
                     (wf%n_J)*batch_length, &
                     wf%t1am,               &
                     wf%n_v,                &
                     zero,                  &
                     L_Jb_a,                &
                     batch_length*(wf%n_J))
!
!        Add terms of L_Jb_a to L_ab_J
!
         do a = 1, wf%n_v
            do b = 1, batch_length
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  Jb = index_two(J, b, wf%n_J)
                  ab = index_two(a, b, wf%n_v)
!
                  L_ab_J(ab, J) = L_ab_J(ab, J) + L_Jb_a(Jb, a)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_Jb,i and L_Jb_a for batch of b
!
         call deallocator(L_Jb_a, (wf%n_J)*batch_length, wf%n_v)
         call deallocator(L_Jb_i, (wf%n_J)*batch_length, wf%n_o)
!
      endif
!
   end subroutine get_cholesky_ab_cc_singles
!
!
end submodule cholesky
