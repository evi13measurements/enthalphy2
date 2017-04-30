module ccs_class
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!                                                                 
!           Coupled cluster singles (CCS) class module                                 
!  Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017  
!                                                                 
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!             -::- Modules used by the class -::-
!
!  General tools
!
   use types
   use utils
   use workspace
   use input_output
!
!  The ancestor class module (HF)
!
   use hf_class
!
   implicit none 
!
!            -::- Definition of the CCS class -::- 
!
   type, extends(hartree_fock) :: cc_singles
!
!     Amplitude attributes
!
      integer(i15) :: n_t1am = 0                    ! Number of singles amplitudes
      real(dp), dimension(:,:), allocatable :: t1am ! Singles amplitude vector
!
   contains 
!
!     Initialization and driver routines
!
      procedure :: init => init_cc_singles
      procedure :: drv  => drv_cc_singles
!
!     Initialization routine for the (singles) amplitudes
!      
      procedure :: initialize_amplitudes => initialize_amplitudes_cc_singles
!
!     A Fock constructor to calculate the T1-transformed Fock  matrix for the current singles 
!     amplitudes (see the HF class for the Fock matrix itself)
!
      procedure, non_overridable :: fock_constructor => fock_constructor_cc_singles
!
!     get Cholesky routines to calculate the occ/vir-occ/vir
!     blocks of the T1-transformed Cholesky vectors
!
      procedure, non_overridable :: get_cholesky_ij => get_cholesky_ij_cc_singles ! occ-occ
      procedure, non_overridable :: get_cholesky_ia => get_cholesky_ia_cc_singles ! occ-vir
      procedure, non_overridable :: get_cholesky_ai => get_cholesky_ai_cc_singles ! vir-occ
      procedure, non_overridable :: get_cholesky_ab => get_cholesky_ab_cc_singles ! vir-vir
!
!     T1-transformator of the one-electron integrals h_pq
!
      procedure, non_overridable :: one_electron_t1 => one_electron_t1_cc_singles
!
   end type cc_singles
!
contains
!
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!              Initialization and driver routines
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! 
   subroutine init_cc_singles(wf)
!
!     Initialize CCS object
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Performs the following tasks
!
!     1. Sets HF orbital and energy information by reading from file (read_hf_info)
!     2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!     3. Allocates the Fock matrix and sets it to zero (note: the matrix is constructed in 
!        the descendant classes) 
!     4. Allocates the singles amplitudes and sets them to zero, and sets associated properties 
!
      implicit none 
!
      class(cc_singles) :: wf
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf % read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf % read_transform_cholesky
!
!     Allocate Fock matrix and set to zero
!
      call wf % initialize_fock_matrix
!
!     Initialize amplitudes and associated attributes
!
      call wf % initialize_amplitudes
!
   end subroutine init_cc_singles
!
!
   subroutine drv_cc_singles(wf)
!
!     CCS Driver
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Lets the user know there is no driver for CCS and exits
!     the program if called.
!
      implicit none 
!
      class(cc_singles) :: wf
!
      write(unit_output,*) 'ERROR: There is no driver for the CCS class.'
      call exit
!
   end subroutine drv_cc_singles
!
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!               Class subroutines and functions 
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
   subroutine initialize_amplitudes_cc_singles(wf)
!
!     Initialize Amplitudes (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Allocates the singles amplitudes, sets them to zero, and calculates
!     the number of singles amplitudes.
!
      implicit none 
!
      class(cc_singles) :: wf
!
!     Calculate the number of singles amplitudes
!
      wf % n_t1am = (wf % n_o)*(wf % n_v) 
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator(wf % t1am, wf % n_t1am, 1)
      wf % t1am = zero
!
   end subroutine initialize_amplitudes_cc_singles
!
!
   subroutine fock_constructor_cc_singles(wf)
!
!     Fock Constructor
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Constructs the T1-transformed Fock matrix blocks (occ/vir-occ/vir),
!     and saves the result in the class variables fock_matrix_pq (see the
!     Hartree-Fock class for these variables)  
!
      implicit none
!
      class(cc_singles) :: wf
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
      integer(i15) :: i = 0, j = 0, k = 0, a = 0, b = 0, ij = 0, kk = 0, ik = 0
      integer(i15) :: kj = 0, ii = 0, jj = 0, ji = 0, ai = 0, ib = 0, bi = 0, ia = 0
      integer(i15) :: aj = 0, ja = 0, ab = 0
!
!     Useful orbital information
!      
      integer(i15) :: n_oo = 0, n_ov = 0, n_vv = 0 ! n_o*n_o, n_o*n_v, n_v*n_v
      integer(i15) :: n_ao_sq_packed = 0           ! Dimension of packed (n_ao x n_ao) matrix
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
      call allocator(h1mo, wf % n_mo, wf % n_mo)
!
      h1mo = zero
!
      call allocator(fock_matrix, wf % n_mo, wf % n_mo)
!
!
!     -::- One-electron contribution -::-
!
!
!     Allocate for one-electron ao integrals
!
      n_ao_sq_packed = packed_size(wf % n_ao)
!
      call allocator(h1ao, n_ao_sq_packed, 1)
      h1ao = zero
!
!     Open mlcc_aoint file
!
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
      call allocator(fock_ao, wf % n_ao, wf % n_ao)
      call squareup(h1ao, fock_ao, wf % n_ao)   
!
!     Deallocation of one-electron AO integrals
!   
      call deallocator(h1ao, n_ao_sq_packed, 1) 
!
!     Transform to one-electron part to MO basis and save it
!     in the fock_matrix 
!
      call allocator(X, wf % n_ao, wf % n_mo)
!
      call dgemm('N','N',        &
                  wf % n_ao,    &
                  wf % n_mo,    &
                  wf % n_ao,    &
                  one,           &
                  fock_ao,       &
                  wf % n_ao,    &
                  wf % mo_coef, &
                  wf % n_ao,    &
                  zero,          &
                  X,             &
                  wf % n_ao)
!
      call dgemm('T','N',        &        
                  wf % n_mo,    &
                  wf % n_mo,    &
                  wf % n_ao,    &
                  one,           &
                  wf % mo_coef, &
                  wf % n_ao,    &
                  X,             &
                  wf % n_ao,    &
                  zero,          &
                  h1mo,          &
                  wf % n_mo)
! 
!     T1-transformation of one-electron integrals in MO basis
!
      call wf % one_electron_t1(h1mo)
      call deallocator(h1mo, wf % n_mo, wf % n_mo)
!
!     Deallocate intermediate X and fock_ao
!
      call deallocator(X,wf % n_ao, wf % n_mo)
      call deallocator(fock_ao, wf % n_ao , wf % n_ao)
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
      n_oo = (wf % n_o) * (wf % n_o)
!
      call allocator(L_ij_J, n_oo, wf % n_J)
      call allocator(g_ij_kl, n_oo, n_oo)
!
      L_ij_J  = zero
      g_ij_kl = zero
!
!     Read Cholesky IJ vector
!
      call wf % get_cholesky_ij(L_ij_J)
!
!     Calculate g_ij_kl
! 
      call dgemm('N','T',    &
                  n_oo,      &
                  n_oo,      &
                  wf % n_J, &
                  one,       &
                  L_ij_J,    &
                  n_oo,      &
                  L_ij_J,    &
                  n_oo,      &
                  zero,      &
                  g_ij_kl,   &
                  n_oo)
!
!     Add two-electron contributions to occupied-occupied block
!
      do i = 1, wf % n_o
         do j = 1, wf % n_o
!
            ij = index_two(i, j, wf % n_o)
!
            do k = 1, wf % n_o
!
               kk = index_two(k,k, wf % n_o)
               ik = index_two(i,k, wf % n_o)
               kj = index_two(k,j, wf % n_o)
!
               fock_matrix(i,j) = fock_matrix(i,j) + two*g_ij_kl(ij,kk) - g_ij_kl(ik,kj)
!
            enddo
!
         enddo
!
      enddo
!
!     Deallocate g_ij_kl
! 
      call deallocator(g_ij_kl, n_oo, n_oo)
!
!
!    :: Occupied-virtual blocks ::
!
!
!     Allocation for g_ia_jk 
!
      n_ov = (wf % n_o) * (wf % n_v)
!
      call allocator(L_ia_J, n_ov, wf % n_J)
      call allocator(g_ia_jk, n_ov, n_oo)
!
!     Read Cholesky vector L_ia_J
!
      call wf % get_cholesky_ia(L_ia_J)
!
!     Calculate g_ia_jk
!
      call dgemm('N','T',    &
                  n_ov,      &
                  n_oo,      &
                  wf % n_J, &
                  one,       &
                  L_ia_J,    &
                  n_ov,      &
                  L_ij_J,    &
                  n_oo,      &
                  zero,      &
                  g_ia_jk,   &
                  n_ov)
!
!     Dealllocate L_ia_J
!
      call deallocator(L_ia_J, n_ov, wf % n_J)
!
!
!     Allocation for g_ai_jk 
!
      call allocator(L_ai_J, n_ov, wf % n_J)
      call allocator(g_ai_jk, n_ov, n_oo)
!
!     Get Cholesky AI vector
!
      call wf % get_cholesky_ai(L_ai_J)
!
!     Calculate g_ai_jk
!
      call dgemm('N','T',    &
                  n_ov,      &
                  n_oo,      &
                  wf % n_J, &
                  one,       &  
                  L_ai_J,    &
                  n_ov,      &
                  L_ij_J,    &
                  n_oo,      &
                  zero,      &
                  g_ai_jk,   &
                  n_ov)
!
!     Deallocate L_ai_J
!
      call deallocator(L_ai_J, n_ov, wf % n_J)
!
!     Add terms to Fock matrix
!
      do i = 1, wf % n_o
         do a = 1, wf % n_v
            do j = 1, wf % n_o
!
!              Needed indices
!
               ia = index_two(i, a, wf % n_o)
               ja = index_two(j, a, wf % n_o)
!
               ai = index_two(a, i, wf % n_v)
               aj = index_two(a, j, wf % n_v)
!
               jj = index_two(j, j, wf % n_o)
               ji = index_two(j, i, wf % n_o)
               ij = index_two(i, j, wf % n_o)
! 
!              Set the blocks of the Fock matrix
!
               fock_matrix(i, a + wf % n_o) = fock_matrix(i,a + wf%n_o) + & 
                                                 two*g_ia_jk(ia,jj) - g_ia_jk(ja,ij) ! g_ia_jk(ja,ij) = g_jaij = g_ijja
!
               fock_matrix(a + wf % n_o, i) = fock_matrix(a + wf % n_o, i) + & 
                                                 two*g_ai_jk(ai,jj) - g_ai_jk(aj,ji)
!
            enddo
         enddo
      enddo
!
      call deallocator(g_ia_jk,n_ov,n_oo)
      call deallocator(g_ai_jk,n_ov,n_oo)
!
!
!     :: Virtual-virtual block F_ab = h_ab + sum_k (2*g_abkk - g_akkb) ::
!
!
      n_vv = (wf % n_v)*(wf % n_v)
!
      call allocator(g_ab_ij, n_vv, n_oo)
      g_ab_ij = zero
!
!     Batch over index a
!
      available = get_available()
!
      required  = 2*(wf % n_v)*(wf % n_v)*(wf % n_J)*4 + &
                  2*(wf % n_v)*(wf % n_o)*(wf % n_J)*4
!
      call num_batch(required, available, max_batch_length, n_batches, wf % n_v)
!
      batch_start  = 1
      batch_end    = 0
      batch_length = 0
!
!     Loop over the batches
!
      do b_batch = 1, n_batches ! E: is this a batch over a or b?
!
!        Get batch limits and length of batch
!
         call batch_limits(batch_start, batch_end, b_batch, max_batch_length, wf % n_v)
         batch_length = batch_end - batch_start + 1
!
!        Allocate L_ab_J
!
         call allocator(L_ab_J, (wf % n_v)*batch_length, wf % n_J)
         L_ab_J = zero
!
!        Read Cholesky vectors
!
         call wf % get_cholesky_ab(L_ab_J, batch_start, batch_end, (wf % n_v)*batch_length, .false.)
!
!        Calculate g_ab_ij = sum_J L_ab_J*L_ij_J
!
         g_off = index_two(1, batch_start, wf % n_v)
!

         call dgemm('N','T',                     &
                     (wf % n_v)*batch_length, &
                     n_oo,                       &
                     wf % n_J,                  &
                     one,                        &
                     L_ab_J,                     &
                     (wf % n_v)*batch_length, &
                     L_ij_J,                     &
                     n_oo,                       &
                     one,                        &
                     g_ab_ij(g_off,1),           &
                     n_vv)
!
!        Deallocate L_ab_J
!
         call deallocator(L_ab_J, batch_length*(wf % n_v), wf % n_J)
!
      enddo ! batching done
!
!     Deallocate L_ij_J
!
      call deallocator(L_ij_J, n_oo, wf % n_J)
!
!     Allocate for g_ai_jb
!
      call allocator(g_ai_jb, n_ov, n_ov)
      call allocator(L_ai_J, n_ov, wf % n_J)
      call allocator(L_ia_J, n_ov, wf % n_J)
!
!     Read Cholesky vectors L_ia_J and L_ai_J
!
      call wf % get_cholesky_ai(L_ai_J)
!
      call dgemm('N','T',    &
                  n_ov,      &
                  n_ov,      &
                  wf % n_J, &
                  one,       &
                  L_ai_J,    &
                  n_ov,      &
                  L_ia_J,    &
                  n_ov,      &
                  zero,      &
                  g_ai_jb,   &
                  n_ov)
!
!     Deallocate L_ia_J
!
     call deallocator(L_ia_J, n_ov, wf % n_J)
     call deallocator(L_ai_J, n_ov, wf % n_J)
!
!     Calculate two-electron terms for virtual-virtual blocks
!
      call dgemm('N','T',  &
                  n_ov,    &
                  n_ov,    &
                  wf%n_J, &
                  one,     &
                  L_ai_J,  &
                  n_ov,    &
                  L_ia_J,  &
                  n_ov,    &
                  zero,    &
                  g_ai_jb, &
                  n_ov)
!
      do a = 1, wf % n_v
         do b = 1, wf % n_v
!
            ab = index_two(a, b, wf%n_v)
!
            do i = 1, wf % n_o
!
               ii = index_two(i, i, wf % n_o)
               ai = index_two(a, i, wf % n_v)
               bi = index_two(b, i, wf % n_v)
               ia = index_two(i, a, wf % n_o)
               ib = index_two(i, b, wf % n_o)
!
               fock_matrix(wf % n_o + a, wf % n_o + b) = fock_matrix(wf % n_o + a, wf % n_o + b) &
                                                             + two*g_ab_ij(ab,ii) - g_ai_jb(ai,ib)
!
            enddo
!
         enddo 
      enddo
!
     call deallocator(g_ab_ij, n_vv, n_oo)
     call deallocator(g_ai_jb, n_ov, n_ov)
!
!     Save the blocks of the Fock matrix in memory (ij,ia,ai,ab)
!
      do i = 1, wf % n_o
         do j = 1, wf % n_o
!
            wf % fock_matrix_ij(i,j) = fock_matrix(i,j)
!
         enddo
      enddo
!
      do i = 1, wf % n_o
         do a = 1, wf % n_v
!
            wf % fock_matrix_ia(i,a) = fock_matrix(i, wf % n_o + a)
            wf % fock_matrix_ai(a,i) = fock_matrix(wf % n_o + a, i)
!
         enddo
      enddo
!
      do a = 1, wf % n_v
         do b = 1, wf % n_v
!
            wf % fock_matrix_ab(a,b) = fock_matrix(wf % n_o + a, wf % n_o + b)
!
         enddo
      enddo
!
      call deallocator(fock_matrix, wf % n_mo, wf % n_mo)
!
   end subroutine fock_constructor_cc_singles
!
!
   subroutine one_electron_t1_cc_singles(wf,h1)
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
      class(cc_singles) :: wf
!
      real(dp), dimension(wf % n_mo, wf % n_mo) :: h1
      real(dp), dimension(:,:), allocatable       :: h1_T1
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
      call allocator(t1, wf % n_mo, wf % n_mo)
      t1 = zero
!
      call allocator(y, wf % n_mo, wf % n_mo)
      call allocator(x, wf % n_mo, wf % n_mo)
!
      y = zero
      x = zero
!
!     Set t1_p_q = t1am_p_q for p virtual and q occupied, 0 otherwise
!
      do a = 1, wf % n_v
         do i = 1, wf % n_o
!
            t1(wf % n_o + a, i) = wf % t1am(a, i)
!
         enddo
      enddo
!
!     Form the x and y arrays
!
      do p = 1, wf % n_mo
         do q = 1, wf % n_mo
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
      call deallocator(t1, wf % n_mo, wf % n_mo)
!
!     Allocate Z intermediate
!   
      call allocator(Z, wf % n_mo, wf % n_mo)
!
!     Calculate h1_T1 = x*h1*y^T = x*Z
!
      call dgemm('N','T',     &
                  wf % n_mo, &
                  wf % n_mo, &
                  wf % n_mo, &
                  one,        &
                  h1,         &
                  wf % n_mo, &
                  y,          &
                  wf % n_mo, &
                  zero,       &
                  Z,          &
                  wf % n_mo)
!
      call dgemm('N','N',     &
                  wf % n_mo, &
                  wf % n_mo, &
                  wf % n_mo, &
                  one,        &
                  x,          &
                  wf % n_mo, &
                  Z,          &
                  wf % n_mo, &
                  zero,       &
                  h1_T1,      &
                  wf % n_mo)
!
!     Deallocate x and y, and the intermediate Z
!
      call deallocator(Z, wf % n_mo, wf % n_mo)
      call deallocator(y, wf % n_mo, wf % n_mo)
      call deallocator(x, wf % n_mo, wf % n_mo)
!
   end subroutine one_electron_t1_cc_singles
!
!
   subroutine get_cholesky_ij_cc_singles(wf,L_ij_J)
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
      real(dp), dimension((wf % n_o)**2, wf % n_J) :: L_ij_J
!
      real(dp), dimension(:,:), allocatable :: L_ia_J ! L_ia^J
      real(dp), dimension(:,:), allocatable :: L_iJ_a ! L_ia^J reordered
      real(dp), dimension(:,:), allocatable :: L_iJ_k ! L_ik^J reordered
!
      integer(i15) :: i = 0, J = 0, a = 0, ij = 0, ia = 0, ik = 0, k = 0
!
      integer(i15) :: n_ov = 0
!
!     Calculate n_o * n_v = n_ov
!
      n_ov = (wf % n_o)*(wf % n_v)
!
!     Allocate
!
      call allocator(L_ia_J, n_ov, wf % n_J)
      call allocator(L_iJ_a, (wf % n_o)*(wf % n_J), wf % n_v)
!
      L_ia_J = zero
      L_iJ_a = zero
!
!     Read the untransformed Cholesky vectors 
!
      call wf % read_cholesky_ij(L_ij_J)
      call wf % read_cholesky_ia(L_ia_J)
!
!     Reorder L_ia_J to L_iJ_a 
!
      do i = 1, wf % n_o
         do J = 1, wf % n_J
            do a = 1, wf % n_v
!              
!              Needed indices
!
               iJ = index_two(i, J, wf % n_o)
               ia = index_two(i, a, wf % n_o)
!
               L_iJ_a(iJ, a) = L_ia_J(ia, J)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_ia_J
!
      call deallocator(L_ia_J, n_ov, wf % n_J)
!
!     Allocate L_iJ_k
!
      call allocator(L_iJ_k, (wf % n_o)*(wf % n_J), wf % n_o)
      L_iJ_k = zero
!
!     T1-transformation
!
      call dgemm('N','N',                    &
                  (wf % n_o)*(wf % n_J), &
                  wf % n_o,               &
                  wf % n_v,               &
                  one,                       &
                  L_iJ_a,                    &
                  (wf % n_o)*(wf % n_J), &
                  wf % t1am,                &
                  wf % n_v,               &
                  zero,                      &
                  L_iJ_k,                    &
                  (wf % n_o)*(wf % n_J))
!
!     Place terms from L_iJ_k into L_ij_J
!
      do i = 1, wf % n_o
         do k = 1, wf % n_o
            do J = 1, wf % n_J
!              
!              Needed indices
!
               iJ = index_two(i, J, wf % n_o)
               ik = index_two(i, k, wf % n_o)
!
               L_ij_J(ik, J) = L_ij_J(ik, J) + L_iJ_k(iJ, k)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_iJ_k and L_iJ_a
!
      call deallocator(L_iJ_k, (wf % n_o)*(wf % n_J), wf % n_o)
      call deallocator(L_iJ_a, (wf % n_o)*(wf % n_J), wf % n_v)
!    
   end subroutine get_cholesky_ij_cc_singles
!
   subroutine get_cholesky_ia_cc_singles(wf,L_ia_J)
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
      real(dp), dimension((wf % n_o)*(wf % n_v), wf % n_J) :: L_ia_J
!
      call wf % read_cholesky_ia(L_ia_J)
!
   end subroutine get_cholesky_ia_cc_singles
!
!
   subroutine get_cholesky_ai_cc_singles(wf,L_ai_J)
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
      real(dp), dimension((wf % n_v)*(wf % n_o), wf % n_J) :: L_ai_J 
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
      call wf % read_cholesky_ai(L_ai_J)
!
!                       
!     -::- L_ab_J contributions -::-
!
!
!     Allocate L_Ja_i
!
      call allocator(L_Ja_i, (wf % n_J)*(wf % n_v), wf % n_o)
      L_Ja_i = zero
!
!     Set batching variables 
!
      required = 2*(wf % n_v)**2*(wf % n_J)*4
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
      call num_batch(required, available, max_batch_length, n_batch, wf % n_v)
!
      do a_batch = 1, n_batch
!
!        Get start, end, and length of batch
!
         call batch_limits(batch_start, batch_end, a_batch, max_batch_length, wf % n_v)
         batch_length = batch_end - batch_start + 1
!
!        Allocate L_ab_J and L_Ja_b
!
         call allocator(L_ba_J, (wf % n_v)*batch_length, wf % n_J) ! L_ab^J = L_ba_J(ba,J)
         call allocator(L_Ja_b, batch_length*(wf % n_J), wf % n_v)
!
         L_ba_J = zero
         L_Ja_b = zero 
!
!        Read Cholesky AB vectors, batching over a
! 
         reorder = .true.
         call wf % read_cholesky_ab(L_ba_J, batch_start, batch_end,&
                                     (wf % n_v)*batch_length, reorder)
!
!        Reorder the Cholesky array L_ba_J
!
         do a = 1, batch_length
            do b = 1, wf % n_v
               do J = 1, wf % n_J
!
!                 Needed indices
!
                  ba = index_two(b, a, wf % n_v)
                  Ja = index_two(J, a, wf % n_J)
!
                  L_Ja_b(Ja, b) = L_ba_J(ba, J) ! L_ab^J
!
               enddo
            enddo
         enddo         
!
!        Calculate sum_b L_Ja_b*t_b_i = L_Ja_i 
!        
         L_off = index_two(1, batch_start, wf % n_J)
!
         call dgemm('N','N',                  &
                     batch_length*(wf % n_J), &
                     wf % n_o,                &
                     wf % n_v,                & 
                     one,                     &
                     L_Ja_b,                  &
                     batch_length*(wf % n_J), &
                     wf % t1am,               &
                     wf % n_v,                &
                     one,                     &
                     L_Ja_i(L_off, 1),        &
                     (wf % n_v)*(wf % n_J))
!
!        Deallocate L_ab_J and L_Ja_b
!
         call deallocator(L_ba_J, (wf % n_v)*batch_length, wf % n_J)
         call deallocator(L_Ja_b, batch_length*(wf % n_J), wf % n_v)
!
      enddo ! batching over a 
!
!     Add terms to T1-transformed Cholesky AI vector 
!
      do i = 1, wf % n_o
         do a = 1, wf % n_v
            do J=1, wf % n_J
!
!              Needed indices
!
               Ja = index_two(J, a, wf % n_J)
               ai = index_two(a, i, wf % n_v)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_Ja_i(Ja, i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_Ja_i
!
      call deallocator(L_Ja_i, (wf % n_J)*(wf % n_v), wf % n_o)
!
!
!     -::- L_ij_J contributions -::-
!
!
!     Allocate L_a_iJ, L_ik_J, L_k_iJ
!
      call allocator(L_a_iJ, wf % n_v, (wf % n_J)*(wf % n_o))
      call allocator(L_ik_J, (wf % n_o)**2, wf % n_J)
      call allocator(L_k_iJ, wf % n_o, (wf % n_o)*(wf % n_J))
!
      L_a_iJ = zero   
      L_ik_J = zero 
      L_k_iJ = zero
!  
!     Read Cholesky IJ vectors
!
      call wf % read_cholesky_ij(L_ik_J) ! L_ik_J(ik,J) = L_ik^J 
!
!     Reorder IJ Cholesky vectors
!
      do i = 1, wf % n_o
         do k = 1, wf % n_o
            do J = 1, wf % n_J
!
!              Needed indices
!
               ik = index_two(i, k, wf % n_o)
               iJ = index_two(i, J, wf % n_o)
!
               L_k_iJ(k, iJ) = L_ik_J(ik, J) ! L_k_iJ(k,iJ) = L_ik^J
!
            enddo
         enddo
      enddo
!
!     Calculate -sum_k t_a_k*L_k_iJ = L_a_iJ  ! Here we assume L_ik^J = L_ki^J - this is safe, right?
!
      call dgemm('N','N',                &
                  wf % n_v,              &
                  (wf % n_o)*(wf % n_J), &
                  wf % n_o               &
                  -one,                  &
                  wf % t1am,             &
                  wf % n_v,              &
                  L_k_iJ,                &
                  wf % n_o,              &
                  zero,                  &
                  L_a_iJ,                &
                  wf % n_v)
!
!     Add terms to T1-transformation of L_ai_J
!
      do i = 1, wf % n_o
         do a = 1, wf % n_v
            do J = 1, wf % n_J
!
!              Needed indices
!
               ai = index_two(a, i, wf % n_v)
               iJ = index_two(i, J, wf % n_o)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_a_iJ, L_ik_J, L_k_iJ
!
      call deallocator(L_a_iJ, wf % n_v, (wf % n_J)*(wf % n_o))      
      call deallocator(L_ik_J, (wf % n_o)**2, wf % n_J)
      call deallocator(L_k_iJ, wf % n_o, (wf % n_o)*(wf % n_J))
!
!
!     -::- L_jb_J contributions -::-    
!
!
      call allocator(L_kJ_b, (wf % n_o)*(wf % n_J), wf % n_v)
      call allocator(L_kb_J, (wf % n_o)*(wf % n_v), wf % n_J)
!
!     Read the Cholesky vector L_kb_J
!
      call wf % read_cholesky_ia(L_kb_J)
!
!     Reorder L_kb_J to L_kJ_b
!
      do k = 1, wf % n_o
         do b = 1, wf % n_v
            do J = 1, wf % n_J
!
!              Needed indices
!
               kb = index_two(k, b, wf % n_o)
               kJ = index_two(k, J, wf % n_o)
!
               L_kJ_b(kJ, b) = L_kb_J(kb, J)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_kb_J
!
      call deallocator(L_kb_J, (wf % n_o)*(wf % n_v), wf % n_J)
!
!     Allocate L_kJ_i 
!
      call allocator(L_kJ_i, (wf % n_o)*(wf % n_J), wf % n_o)
      L_kJ_i = zero
!
!     Calculate sum_b L_kJ_b*t_b_i = L_kJ_i
!
      call dgemm('N','N',                &
                  (wf % n_o)*(wf % n_J), &
                  wf % n_o,              &
                  wf % n_v,              &
                  one,                   &
                  L_kJ_b,                &
                  (wf % n_o)*(wf % n_J), &
                  wf % t1am,             &
                  wf % n_v,              &
                  zero,                  &
                  L_kJ_i,                &
                  (wf % n_o)*(wf % n_J))
!
!     Deallocate L_kJ_b
!
      call deallocator(L_kJ_b, (wf % n_o)*(wf % n_J), wf % n_v)
!
!     Allocate L_k_iJ
!  
      call allocator(L_k_iJ, (wf % n_o), (wf % n_o)*(wf % n_J))
      L_k_iJ = zero
!
!     Reorder L_kJ_i to L_k_iJ    
!
      do i = 1, wf % n_o
         do k = 1, wf % n_o
            do J = 1, wf % n_J
!
!              Needed indices
!
               kJ = index_two(k, J, wf % n_o)
               iJ = index_two(i, J, wf % n_o)
!
               L_k_iJ(k, iJ) = L_kJ_i(kJ, i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_kJ_i
!
      call deallocator(L_kJ_i, (wf % n_o)*(wf % n_J), wf % n_o)
!
!     Allocate L_a_iJ
!
      call allocator(L_a_iJ, wf % n_v, (wf % n_o)*(wf % n_J))
      L_a_iJ = zero
!      
!     Calculate sum_k t_a_k*L_k_iJ = L_a_iJ
!
      call dgemm('N','N',                &
                  wf % n_v,              &
                  (wf % n_o)*(wf % n_J), &
                  wf % n_o,              &
                  -one,                  &
                  wf % t1am,             &
                  wf % n_v,              &
                  L_k_iJ,                &
                  wf % n_o,              &
                  zero,                  &
                  L_a_iJ,                &
                  wf % n_v)
!
!     Add contribution to L ai_J
!
      do a = 1, wf % n_v
         do i = 1, wf % n_o
            do J = 1, wf % n_J
!
!              Needed indices
!
               iJ = index_two(i, J, wf % n_o)
               ai = index_two(a, i, wf % n_v)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
            enddo
         enddo
      enddo
!
!     Deallocations
!
      call deallocator(L_a_iJ, wf % n_v, (wf % n_o)*(wf % n_J))
      call deallocator(L_k_iJ, wf % n_o, (wf % n_o)*(wf % n_J))
!
   end subroutine get_cholesky_ai_cc_singles
!
!
   subroutine get_cholesky_ab_cc_singles(wf,L_ab_J,first,last,ab_dim,reorder)
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
      real(dp), dimension(ab_dim, wf % n_J) :: L_ab_J
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
         call allocator(L_ib_J, (wf % n_o)*(wf % n_v), wf % n_J)
         L_ib_J = zero
!
!        Read L_ia_J
!  
         call wf % read_cholesky_ia(L_ib_J) ! Note: using L_ia_J instead of L_ai_J, here, to avoid two reorderings.
                                          ! This is possible because of the symmetry L_ai_J(ai,J) == L_ia_J(ia,J).
!
!        Read L_ab_J for batch of a
!
         call wf % read_cholesky_ab(L_ab_J, first, last, ab_dim, reorder) ! Eirik: L_ab_J(ba) = L_ab^J 
!
!        Allocate L_i,Jb
!
         call allocator(L_i_Jb, wf % n_o, (wf % n_J)*(wf % n_v))
         L_i_Jb = zero
!
!        Reorder L_ib_J to L_i_Jb
!
         do i = 1, wf % n_o
            do b = 1, wf % n_v
               do J = 1, wf % n_J
!
!                 Needed indices
!
                  ib = index_two(i, b, wf % n_o)
                  Jb = index_two(J, b, wf % n_J)
!
                  L_i_Jb(i, Jb)=L_ib_J(ib, J)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator(L_ib_J, (wf % n_o)*(wf % n_v), wf % n_J)
!
!        Allocate L_a_Jb for batch of a
!  
         call allocator(L_a_Jb, batch_length, (wf % n_v)*(wf % n_J))
!
!        Calculate  -t1_a_i * L_i_Jb = L_a_Jb
!
         call dgemm('N','N',                &
                     batch_length,          &
                     (wf % n_v)*(wf % n_J), &
                     wf % n_o,              &
                     -one,                  &
                     wf % t1am(first,1),    &
                     wf % n_v,              &
                     L_i_Jb,                &
                     wf % n_o,              &
                     zero,                  &
                     L_a_Jb,                &
                     batch_length)
!
!        Add terms of L_a_Jb to L_ab_J
!
         do a = 1, batch_length
            do b = 1, wf % n_v
               do J = 1, wf % n_J
!
!                 Needed indices
!
                  Jb = index_two(J, b, wf % n_J)
                  ba = index_two(b, a, wf % n_v)
!
                  L_ab_J(ba, J) = L_ab_J(ba, J) + L_a_Jb(a, Jb) 
!
               enddo
            enddo
         enddo
!
!        Dellocate L_i_Jb and L_a_Jb for batch of a
!
         call deallocator(L_a_Jb, batch_length, (wf % n_J)*(wf % n_v))
         call deallocator(L_i_Jb, (wf % n_J)*(wf % n_v),(wf % n_o))
!
      else  !! Batching over b !!
!
!        Allocate L_ib_J
!     
         call allocator(L_ib_J, (wf % n_o)*(wf % n_v), wf % n_J)
         L_ib_J = zero
!
!        Read L_ia_J
!  
         call wf % read_cholesky_ia(L_ib_J) ! Note: using L_ia_J instead of L_ai_J, here, to avoid two reorderings.
                                             ! This is possible because of the symmetry L_ai_J(ai,J) == L_ia_J(ia,J)
!
!        Read L_ab_J for batch of b
!
         call wf % read_cholesky_ab(L_ab_J, first, last, ab_dim, reorder)
!
!        Allocate L_Jb,i for batch of b
!
         call allocator(L_Jb_i, (wf % n_J)*batch_length, wf % n_o)
         L_Jb_i = zero
!
!        Reorder L_ib_J to L_Jb_i
!
         do i = 1, wf % n_o
            do b = 1, batch_length
               do J = 1, wf % n_J
!
!                 Needed indices
!
                  ib = index_two(i, b + first - 1, wf % n_o) ! Note: in L_ib_J we have all b, not the case for L_Jb_i
                  Jb = index_two(J, b, wf % n_J)
!
                  L_Jb_i(Jb, i) = L_ib_J(ib, J)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator(L_ib_J, (wf % n_o)*(wf % n_v), wf % n_J)
!
!        Allocate L_Jb_a for batch of b
!  
         call allocator(L_Jb_a, (wf % n_J)*batch_length, wf % n_v)
!
!        T1-transformation
!
         call dgemm('N','T',&
                     (wf % n_J)*batch_length, &
                     wf % n_v,                &
                     wf % n_o,                &
                     -one,                    &
                     L_Jb_i,                  &
                     (wf % n_J)*batch_length, &
                     wf % t1am,               &
                     wf % n_v,                &
                     zero,                    &
                     L_Jb_a,                  &
                     batch_length*(wf % n_J))
!
!        Add terms of L_Jb_a to L_ab_J
!
         do a = 1, wf % n_v
            do b = 1, batch_length
               do J = 1, wf % n_J
!
!                 Needed indices
!
                  Jb = index_two(J, b, wf % n_J)
                  ab = index_two(a, b, wf % n_v)
!
                  L_ab_J(ab, J) = L_ab_J(ab, J) + L_Jb_a(Jb, a)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_Jb,i and L_Jb_a for batch of b
!
         call deallocator(L_Jb_a, (wf % n_J)*batch_length, wf % n_v)
         call deallocator(L_Jb_i, (wf % n_J)*batch_length, wf % n_o)
!
      endif
!
      integer(i15) :: memory_lef
!
      real(dp) :: L_ab_J(ab_dim,wfn%n_J)
!
      integer(i15) :: ba=0
      integer(i15) :: lucho_ab,ab_dim
      integer(i15) :: a=0,b=0,J=0,i=0,ia=0,aJ=0,ib=0,Jb=0,ab=0
      integer(i15) :: n_ov
      integer(i15) :: start,end
      logical,intent(in) :: reorder
!
      real(dp),dimension(:,:),allocatable     :: L_ib_J
      real(dp),dimension(:,:),allocatable     :: L_Jb_i
      real(dp),dimension(:,:),allocatable     :: L_Jb_a
      real(dp),dimension(:,:),allocatable     :: L_a_Jb
      real(dp),dimension(:,:),allocatable     :: L_i_Jb
!
      integer(i15) :: batch_length=0

!
!     calculating length of batch
!
      batch_length = end-start+1
!
!     Usefull variables
!
      n_ov = (wfn%n_occ)*(wfn%n_vir)
!
!
!     Testing which index is batched
!
!
      if (reorder) then !! Batching over a !!
!
!        Allocate L_ib_J
!     
         call allocator(L_ib_J,n_ov,wfn%n_J)
         L_ib_J=zero
!
!        Read L_ia_J
!  
         call read_cholesky_ia(L_ib_J) ! OBS: Using L_ia_J insted of L_ai_J to avoid two reorderings!
                                       ! Possible because of symmetry L_ai_J(ai,J)==L_ia_J(ia,J)
!
!        Read L_ab_J for batch of a
!
         call wfn%read_cholesky_ab(L_ab_J,start,end,ab_dim,.true.) ! L_ab_J(ba) = L_ab^J 
!
!        Allocate L_i,Jb
!
         call allocator(L_i_Jb,wfn%n_occ,(wfn%n_J)*(wfn%n_vir))
         L_i_Jb=zero
!
!        Reorder L_ib_J to L_i_Jb
!
         do i=1,wfn%n_occ
            do b=1,wfn%n_vir
               do J=1,wfn%n_J
!
!                 Needed indices
!
                  ib=index_two(i,b,wfn%n_occ)
                  Jb=index_two(J,b,wfn%n_J)
!
                  L_i_Jb(i,Jb)=L_ib_J(ib,J)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator_n(L_ib_J,n_ov,wfn%n_J)
!
!        Allocate L_a_Jb for batch of a
!  
         call allocator_n(L_a_Jb,batch_length,(wfn%n_vir)*(wfn%n_J))
!
!        - t1_a_i * L_i_Jb = L_a_Jb
!
         call dgemm('N','N',       &
                     batch_length,          &
                     (wfn%n_vir)*(wfn%n_J), &
                     wfn%n_occ,             &
                     -one,                  &
                     wfn%t1am(start,1),     &
                     wfn%n_vir,             &
                     L_i_Jb,                &
                     wfn%n_occ,             &
                     zero,                  &
                     L_a_Jb,                &
                     batch_length)
!
!        Add terms of L_a_Jb to L_ab_J
!
         do a=1,batch_length
            do b=1,wfn%n_vir
               do J=1,wfn%n_J
!
!                 Needed indices
!
                  Jb=index_two(J,b,wfn%n_J)
!
                  ba = index_two(b,a,wfn%n_vir)
                  L_ab_J(ba,J) = L_ab_J(ba,J) + L_a_Jb(a,Jb) 
!
               enddo
            enddo
         enddo
!
!        Dellocate L_i_Jb and L_a_Jb for batch of a
!
         call deallocator_n(L_a_Jb,batch_length,(wfn%n_J)*(wfn%n_vir))
         call deallocator_n(L_i_Jb,(wfn%n_J)*(wfn%n_vir),(wfn%n_occ))
!
      else  !! Batching over b !!
!
!        Allocate L_ib_J
!     
         call allocator_n(L_ib_J,n_ov,(wfn%n_J))
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
         call allocator_n(L_Jb_i,(wfn%n_J)*batch_length,(wfn%n_occ))
         L_Jb_i=zero
!
!        Reorder L_ib_J to L_Jb_i
!
         do i=1,wfn%n_occ
            do b=1,batch_length
               do J=1,wfn%n_J
!
!                 Needed indices
!
                  ib=index_two(i,b+start-1,(wfn%n_occ)) ! OBS: in L_ib_J we have all b, not the case for L_Jb_i
                  Jb=index_two(J,b,wfn%n_J)
!
                  L_Jb_i(Jb,i)=L_ib_J(ib,J) ! The quantity on the right is the entire L_ib^J matrix; the left the batched. OK
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator_n(L_ib_J,n_ov,wfn%n_J)
!
!        Allocate L_Jb_a for batch of b
!  
         call allocator_n(L_Jb_a,wfn%n_J*batch_length,wfn%n_vir)
!
!        T1-transformation
!
         call dgemm('N','T',&
                     wfn%n_J*batch_length, &
                     wfn%n_vir,            &
                     wfn%n_occ,            &
                     -one,                 &
                     L_Jb_i,               &
                     wfn%n_J*batch_length, &
                     wfn%t1am,             &
                     wfn%n_vir,            &
                     zero,                 &
                     L_Jb_a,               &
                     batch_length*wfn%n_J)
!
!        Add terms of L_Jb_a to L_ab_J
!
         do a=1,wfn%n_vir
            do b=1,batch_length
               do J=1,wfn%n_J
!
!                 Needed indices
!
                  Jb=index_two(J,b,wfn%n_J)
                  ab=index_two(a,b,wfn%n_vir)
!
                  L_ab_J(ab,J)=L_ab_J(ab,J)+L_Jb_a(Jb,a)
               enddo
            enddo
         enddo
!
!        Dellocate L_Jb,i and L_Jb_a for batch of b
!
         call deallocator_n(L_Jb_a,wfn%n_J*batch_length,wfn%n_vir)
         call deallocator_n(L_Jb_i,wfn%n_J*batch_length,wfn%n_occ)
!
!
      endif
!
   end subroutine get_cholesky_ab_cc_singles
!
!
end module ccs_class
