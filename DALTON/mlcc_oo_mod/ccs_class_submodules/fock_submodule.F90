submodule (ccs_class) fock
!
!
!                       -::- Fock submodule (CCS) -::-
!            Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!
!  Contains the following family of procedures of the CCS class:
!  
!  initialize_fock_matrix_ccs(wf):		Allocates and sets Fock matrix to 0        
!  construct_fock_ccs(wf):			      Constructs T1_transformed mo Fock matrix 
!  one_electron_t1_ccs(wf, h1 ,h1_T1):	T1-transformation of one-electron mo integrals
!
   implicit none 
!
!
contains
!
!
    subroutine initialize_fock_matrix_ccs(wf)
!
!     Initialize Fock Matrix (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Allocates and constructs the Fock matrix 
!
      implicit none
!  
      class(ccs) :: wf   
!
      call allocator(wf%fock_ij, wf%n_o, wf%n_o)
      call allocator(wf%fock_ia, wf%n_o, wf%n_v)
      call allocator(wf%fock_ai, wf%n_v, wf%n_o)
      call allocator(wf%fock_ab, wf%n_v, wf%n_v)

!
      wf%fock_ij = zero
      wf%fock_ia = zero
      wf%fock_ai = zero
      wf%fock_ab = zero
!
      call wf%construct_fock
!
   end subroutine initialize_fock_matrix_ccs
!
!
   subroutine construct_fock_ccs(wf)
!
!     Construct Fock
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Constructs the T1-transformed Fock matrix blocks (occ/vir-occ/vir),
!     and saves the result in the class variables fock_matrix_pq (see the
!     Hartree-Fock class for these variables)  
!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: fock_ao
      real(dp), dimension(:,:), allocatable :: fock_matrix
!
      real(dp), dimension(:,:), allocatable :: h1ao ! AO basis matrix h_αβ
      real(dp), dimension(:,:), allocatable :: h1mo ! MO basis matrix h_pq  
      real(dp), dimension(:,:), allocatable :: X    ! An intermediate
!      
      integer(i15) :: unit_identifier_ao_integrals = -1 ! Unit identifier for file mlcc_aoint
!
!     Indices
!
      integer(i15) :: i = 0, j = 0, k = 0, a = 0, b = 0
      integer(i15) :: kj = 0, ii = 0, ij = 0, kk = 0, ik = 0, jj = 0, ji = 0, ai = 0, ib = 0, bi = 0, ia = 0
      integer(i15) :: aj = 0, ja = 0, ab = 0
!
!     Useful orbital information
!      
      integer(i15) :: n_ao_sq_packed = 0 ! Dimension of packed (n_ao x n_ao) matrix
!
!     Two electron integrals
!
      real(dp), dimension(:,:), allocatable :: g_ij_kl
      real(dp), dimension(:,:), allocatable :: g_ab_ij
      real(dp), dimension(:,:), allocatable :: g_ai_jb
      real(dp), dimension(:,:), allocatable :: g_ia_jk
      real(dp), dimension(:,:), allocatable :: g_ai_jk
!
!     Cholesky vectors
!
      real(dp), dimension(:,:), allocatable :: L_ij_J 
      real(dp), dimension(:,:), allocatable :: L_ia_J 
      real(dp), dimension(:,:), allocatable :: L_ai_J 
      real(dp), dimension(:,:), allocatable :: L_ab_J 
!
!     Batch settings
!
      integer(i15) :: available = 0, required = 0, max_batch_length = 0
      integer(i15) :: batch_end = 0, batch_length = 0, g_off = 0, n_batches = 0
      integer(i15) :: b_batch = 0, batch_start = 0
!
!     Allocate one-electron MO integrals
!
      call allocator(h1mo, wf%n_mo, wf%n_mo)
      h1mo = zero
!
      call allocator(fock_matrix, wf%n_mo, wf%n_mo)
      fock_matrix = zero
!
!
!     -::- One-electron contribution -::-
!
!
!     Allocate for one-electron ao integrals
!
      n_ao_sq_packed = packed_size(wf%n_ao)
!
      call allocator(h1ao, n_ao_sq_packed, 1)
      h1ao = zero
!
!     Open mlcc_aoint file
!
      call generate_unit_identifier(unit_identifier_ao_integrals)
      open(unit=unit_identifier_ao_integrals,file='mlcc_aoint',status='old',form='formatted')
      rewind(unit_identifier_ao_integrals)
!
!     Read in one-electron AO integrals
!
      read(unit_identifier_ao_integrals,*) (h1ao(i,1), i = 1, n_ao_sq_packed)
!
!     Close mlcc_aoint
!
      close(unit_identifier_ao_integrals)
!
!     Allocate the AO Fock matrix and add the one-electron contributions
!
      call allocator(fock_ao, wf%n_ao, wf%n_ao)
      fock_ao = zero
!
      call squareup(h1ao, fock_ao, wf%n_ao)  
!
!     Deallocation of one-electron AO integrals
!   
      call deallocator(h1ao, n_ao_sq_packed, 1) 
!
!     Transform to one-electron part to MO basis and save it
!     in the fock_matrix 
!
      call allocator(X, wf%n_ao, wf%n_mo)
!
      call dgemm('N','N',     &
                  wf%n_ao,    &
                  wf%n_mo,    &
                  wf%n_ao,    &
                  one,        &
                  fock_ao,    &
                  wf%n_ao,    &
                  wf%mo_coef, &
                  wf%n_ao,    &
                  zero,       &
                  X,          &
                  wf%n_ao)
!
      call dgemm('T','N',     &        
                  wf%n_mo,    &
                  wf%n_mo,    &
                  wf%n_ao,    &
                  one,        &
                  wf%mo_coef, &
                  wf%n_ao,    &
                  X,          &
                  wf%n_ao,    &
                  zero,       &
                  h1mo,       &
                  wf%n_mo)
! 
!     T1-transformation of one-electron integrals in MO basis
!
      call wf%one_electron_t1(h1mo, fock_matrix)
      call deallocator(h1mo, wf%n_mo, wf%n_mo)
!
!     Deallocate intermediate X and fock_ao
!
      call deallocator(X,wf%n_ao, wf%n_mo)
      call deallocator(fock_ao, wf%n_ao, wf%n_ao)
!
!
!     -::- Two-electron contribution -::-
!
!
!     :: Occupied-occupied block: F_ij = h_ij + sum_k (2*g_ijkk - g_ikkj) ::
!  
!
!     Allocation for L_ij_J
! 
      call allocator(L_ij_J, (wf%n_o)**2, wf%n_J)
      call allocator(g_ij_kl, (wf%n_o)**2, (wf%n_o)**2)
!
      L_ij_J  = zero
      g_ij_kl = zero
!
!     Read Cholesky IJ vector
!
      call wf%get_cholesky_ij(L_ij_J)
!
!     Calculate g_ij_kl
! 
      call dgemm('N','T',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  wf%n_J,      &
                  one,         &
                  L_ij_J,      &
                  (wf%n_o)**2, &
                  L_ij_J,      &
                  (wf%n_o)**2, &
                  zero,        &
                  g_ij_kl,     &
                  (wf%n_o)**2)
!
!     Add two-electron contributions to occupied-occupied block
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            ij = index_two(i, j, wf%n_o)
!
            do k = 1, wf%n_o
!
               kk = index_two(k, k, wf%n_o)
               ik = index_two(i, k, wf%n_o)
               kj = index_two(k, j, wf%n_o)
!
               fock_matrix(i,j) = fock_matrix(i,j) + & 
                                    two*g_ij_kl(ij,kk) - &
                                    g_ij_kl(ik,kj)
!
            enddo
!
         enddo
!
      enddo
!
!     Deallocate g_ij_kl
! 
      call deallocator(g_ij_kl, (wf%n_o)**2, (wf%n_o)**2)
!
!
!    :: Occupied-virtual blocks ::
!
!
!     Allocation for g_ia_jk 
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call allocator(g_ia_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      L_ia_J  = zero
      g_ia_jk = zero
!
!     Read Cholesky vector L_ia_J
!
      call wf%get_cholesky_ia(L_ia_J)
!
!     Calculate g_ia_jk
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)**2,       &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ij_J,            &
                  (wf%n_o)**2,       & 
                  zero,              &
                  g_ia_jk,           &
                  (wf%n_o)*(wf%n_v))
!
!     Dealllocate L_ia_J
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocation for g_ai_jk 
!
      call allocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call allocator(g_ai_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      L_ai_J  = zero
      g_ai_jk = zero
!
!     Get Cholesky AI vector
!
      call wf%get_cholesky_ai(L_ai_J)
!
!     Calculate g_ai_jk
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)**2,       &
                  wf%n_J,            &
                  one,               &  
                  L_ai_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ij_J,            &
                  (wf%n_o)**2,       &
                  zero,              &
                  g_ai_jk,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate L_ai_J
!
      call deallocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!     Add terms to Fock matrix
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
!
!              Needed indices
!
               ia = index_two(i, a, wf%n_o)
               ja = index_two(j, a, wf%n_o)
!
               ai = index_two(a, i, wf%n_v)
               aj = index_two(a, j, wf%n_v)
!
               jj = index_two(j, j, wf%n_o)
               ji = index_two(j, i, wf%n_o)
               ij = index_two(i, j, wf%n_o)
! 
!              Set the blocks of the Fock matrix
!
               fock_matrix(i, a + wf%n_o) = fock_matrix(i, a + wf%n_o) + & 
                                                 two*g_ia_jk(ia,jj) - g_ia_jk(ja,ij) ! g_ia_jk(ja,ij) = g_jaij = g_ijja
!
               fock_matrix(a + wf%n_o, i) = fock_matrix(a + wf%n_o, i) + & 
                                                 two*g_ai_jk(ai,jj) - g_ai_jk(aj,ji)
!
            enddo
         enddo
      enddo
!
      call deallocator(g_ia_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      call deallocator(g_ai_jk, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
!
!     :: Virtual-virtual block F_ab = h_ab + sum_k (2*g_abkk - g_akkb) ::
!
!
      call allocator(g_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
      g_ab_ij = zero
!
!     Batch over index b
!
      available = get_available()
!
      required  = 2*(wf%n_v)*(wf%n_v)*(wf%n_J)*4 + &
                  2*(wf%n_v)*(wf%n_o)*(wf%n_J)*4
!
      call num_batch(required, available, max_batch_length, n_batches, wf%n_v)
!
      batch_start  = 1
      batch_end    = 0
      batch_length = 0
!
!     Loop over the batches
!
      do b_batch = 1, n_batches 
!
!        Get batch limits and length of batch
!
         call batch_limits(batch_start, batch_end, b_batch, max_batch_length, wf%n_v)
         batch_length = batch_end - batch_start + 1
!
!        Allocate L_ab_J
!
         call allocator(L_ab_J, (wf%n_v)*batch_length, wf%n_J)
         L_ab_J = zero
!
!        Read Cholesky vectors
!
         call wf%get_cholesky_ab(L_ab_J, batch_start, batch_end, (wf%n_v)*batch_length, .false.)
!
!        Calculate g_ab_ij = sum_J L_ab_J*L_ij_J
!
         g_off = index_two(1, batch_start, wf%n_v)
!
         call dgemm('N','T',                &
                     (wf%n_v)*batch_length, &
                     (wf%n_o)**2,           &
                     wf%n_J,                &
                     one,                   &
                     L_ab_J,                &
                     (wf%n_v)*batch_length, &
                     L_ij_J,                &
                     (wf%n_o)**2,           &
                     one,                   &
                     g_ab_ij(g_off,1),      &
                     (wf%n_v)**2)
!
!        Deallocate L_ab_J
!
         call deallocator(L_ab_J, batch_length*(wf%n_v), wf%n_J)
!
      enddo ! batching done
!
!     Deallocate L_ij_J
!
      call deallocator(L_ij_J, (wf%n_o)**2, wf%n_J)
!
!     Allocate for g_ai_jb
!
      call allocator(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ai_jb = 0
!
      call allocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_ai_J = 0
      L_ia_J = 0
!
!     Read Cholesky vectors L_ia_J and L_ai_J
!
      call wf%get_cholesky_ia(L_ia_J)
      call wf%get_cholesky_ai(L_ai_J)
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ai_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ai_jb,           &
                  (wf%n_o)*(wf%n_v))
!
!     Calculate two-electron terms for virtual-virtual blocks
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ai_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ai_jb,           &
                  (wf%n_o)*(wf%n_v))
!
!    Deallocate L_ia_J
!
     call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
     call deallocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            ab = index_two(a, b, wf%n_v)
!
            do i = 1, wf%n_o
!
               ii = index_two(i, i, wf%n_o)
               ai = index_two(a, i, wf%n_v)
               bi = index_two(b, i, wf%n_v)
               ia = index_two(i, a, wf%n_o)
               ib = index_two(i, b, wf%n_o)
!
               fock_matrix(wf%n_o + a, wf%n_o + b) = fock_matrix(wf%n_o + a, wf%n_o + b) &
                                                       + two*g_ab_ij(ab,ii) &
                                                       - g_ai_jb(ai,ib)
!
            enddo
!
         enddo 
      enddo
!
     call deallocator(g_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
     call deallocator(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Save the blocks of the Fock matrix in memory (ij,ia,ai,ab)
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            wf%fock_ij(i,j) = fock_matrix(i,j)
!
         enddo
      enddo
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%fock_ia(i,a) = fock_matrix(i, wf%n_o + a)
            wf%fock_ai(a,i) = fock_matrix(wf%n_o + a, i)
!
         enddo
      enddo
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            wf%fock_ab(a,b) = fock_matrix(wf%n_o + a, wf%n_o + b)
!
         enddo
      enddo
!
      call deallocator(fock_matrix, wf%n_mo, wf%n_mo)
!
   end subroutine construct_fock_ccs
!
!
   subroutine one_electron_t1_ccs(wf, h1 ,h1_T1)
!
!     One-electron T1 
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     T1-transforms the one-electron MO integrals h_pq
!
!           h_p_q_T1 = sum_st x_p_s * y_q_t * h_s_t
!
!           x = I - t1
!           y = I - t1^T
!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo) :: h1
      real(dp), dimension(wf%n_mo, wf%n_mo) :: h1_T1
!
      real(dp), dimension(:,:), allocatable :: x 
      real(dp), dimension(:,:), allocatable :: y 
      real(dp), dimension(:,:), allocatable :: t1
!
      real(dp), dimension(:,:), allocatable :: Z ! Intermediate for matrix multiplication
!
      integer(i15) :: p = 0, q = 0, a = 0, i = 0
!
!     Allocate the arrays t1, x, and y
!
      call allocator(t1, wf%n_mo, wf%n_mo)
      t1 = zero
!
      call allocator(y, wf%n_mo, wf%n_mo)
      call allocator(x, wf%n_mo, wf%n_mo)
!
      y = zero
      x = zero
!
!     Set t1_p_q = t1am_p_q for p virtual and q occupied, 0 otherwise
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            t1(wf%n_o + a, i) = wf%t1am(a, i)
!
         enddo
      enddo
!
!     Form the x and y arrays
!
      do p = 1, wf%n_mo
         do q = 1, wf%n_mo
!
            if (p .eq. q) then
!
               x(p,q) = 1
               y(p,q) = 1
!
            else
!
               x(p,q) = x(p,q) - t1(p,q)
               y(p,q) = y(p,q) + t1(q,p)
!
            endif
!
         enddo
      enddo
!
!     Deallocate t1 (only x and y are needed below)
!
      call deallocator(t1, wf%n_mo, wf%n_mo)
!
!     Allocate Z intermediate
!   
      call allocator(Z, wf%n_mo, wf%n_mo)
!
!     Calculate h1_T1 = x*h1*y^T = x*Z
!
      call dgemm('N','T',  &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_mo, &
                  one,     &
                  h1,      &
                  wf%n_mo, &
                  y,       &
                  wf%n_mo, &
                  zero,    &
                  Z,       &
                  wf%n_mo)
!
      call dgemm('N','N',  &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_mo, &
                  one,     &
                  x,       &
                  wf%n_mo, &
                  Z,       &
                  wf%n_mo, &
                  zero,    &
                  h1_T1,   &
                  wf%n_mo)
!
!     Deallocate x and y, and the intermediate Z
!
      call deallocator(Z, wf%n_mo, wf%n_mo)
      call deallocator(y, wf%n_mo, wf%n_mo)
      call deallocator(x, wf%n_mo, wf%n_mo)
!
   end subroutine one_electron_t1_ccs
!
!
end submodule
