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
      procedure :: fock_constructor => fock_constructor_cc_singles
!
!     get Cholesky routines to calculate the occ/vir-occ/vir
!     blocks of the T1-transformed Cholesky vectors
!
      procedure :: get_cholesky_ij => get_cholesky_ij_cc_singles ! occ-occ
      procedure :: get_cholesky_ia => get_cholesky_ia_cc_singles ! occ-vir
      procedure :: get_cholesky_ai => get_cholesky_ai_cc_singles ! vir-occ
      procedure :: get_cholesky_ab => get_cholesky_ab_cc_singles ! vir-vir
!
!     T1-transformator of the one-electron integrals h_pq
!
      procedure :: one_electron_t1 => one_electron_t1_cc_singles
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
   subroutine init_cc_singles(wfn)
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
      class(cc_singles) :: wfn
!
!     Read Hartree-Fock info from SIRIUS
!
      call wfn % read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wfn % read_transform_cholesky
!
!     Allocate Fock matrix and set to zero
!
      call wfn % allocate_fock_matrix
!
!     Initialize amplitudes and associated attributes
!
      call wfn % initialize_amplitudes
!
   end subroutine init_cc_singles
!
!
   subroutine drv_cc_singles(wfn)
!
!     CCS Driver
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Lets the user know there is no driver for CCS and exits
!     the program if called.
!
      implicit none 
!
      class(cc_singles) :: wfn
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
   subroutine initialize_amplitudes_cc_singles(wfn)
!
!     Initialize Amplitudes (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Allocates the singles amplitudes, sets them to zero, and calculates
!     the number of singles amplitudes.
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
      call allocator(wfn % t1am, wfn % n_t1am, 1)
      wfn % t1am = zero
!
   end subroutine initialize_amplitudes_cc_singles
!
!
   subroutine fock_constructor_cc_singles(wfn)
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
      class(cc_singles) :: wfn
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
      integer(i15) :: n_oo = 0, n_ov = 0, n_vv = 0 ! n_occ*n_occ, n_occ*n_vir, n_vir*n_vir
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
      call allocator(h1mo, wfn % n_mo, wfn % n_mo)
      h1mo = zero
!
      call allocator(fock_matrix, wfn % n_mo, wfn % n_mo)
!
!
!     -::- One-electron contribution -::-
!
!
!     Allocate for one-electron ao integrals
!
      n_ao_sq_packed = packed_size(wfn % n_ao)
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
      call allocator(fock_ao, wfn % n_ao, wfn % n_ao)
      call squareup(h1ao, fock_ao, wfn % n_ao)   
!
!     Deallocation of one-electron AO integrals
!   
      call deallocator(h1ao, n_ao_sq_packed, 1) 
!
!     Transform to one-electron part to MO basis and save it
!     in the fock_matrix 
!
      call allocator(X, wfn % n_ao, wfn % n_mo)
!
      call dgemm('N','N',        &
                  wfn % n_ao,    &
                  wfn % n_mo,    &
                  wfn % n_ao,    &
                  one,           &
                  fock_ao,       &
                  wfn % n_ao,    &
                  wfn % mo_coef, &
                  wfn % n_ao,    &
                  zero,          &
                  X,             &
                  wfn % n_ao)
!
      call dgemm('T','N',        &        
                  wfn % n_mo,    &
                  wfn % n_mo,    &
                  wfn % n_ao,    &
                  one,           &
                  wfn % mo_coef, &
                  wfn % n_ao,    &
                  X,             &
                  wfn % n_ao,    &
                  zero,          &
                  h1mo,          &
                  wfn % n_mo)
! 
!     T1-transformation of one-electron integrals in MO basis
!
      call wfn % one_electron_t1(h1mo)
      call deallocator(h1mo, wfn % n_mo, wfn % n_mo)
!
!     Deallocate intermediate X and fock_ao
!
      call deallocator(X,wfn % n_ao, wfn % n_mo)
      call deallocator(fock_ao, wfn % n_ao , wfn % n_ao)
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
      n_oo = (wfn % n_occ) * (wfn % n_occ)
!
      call allocator(L_ij_J, n_oo, wfn % n_J)
      call allocator(g_ij_kl, n_oo, n_oo)
!
      L_ij_J  = zero
      g_ij_kl = zero
!
!     Read Cholesky IJ vector
!
      call wfn % get_cholesky_ij(L_ij_J)
!
!     Calculate g_ij_kl
! 
      call dgemm('N','T',    &
                  n_oo,      &
                  n_oo,      &
                  wfn % n_J, &
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
      do i = 1, wfn % n_occ
         do j = 1, wfn % n_occ
!
            ij = index_two(i, j, wfn % n_occ)
!
            do k = 1, wfn % n_occ
!
               kk = index_two(k,k, wfn % n_occ)
               ik = index_two(i,k, wfn % n_occ)
               kj = index_two(k,j, wfn % n_occ)
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
      n_ov = (wfn % n_occ) * (wfn % n_vir)
!
      call allocator(L_ia_J, n_ov, wfn % n_J)
      call allocator(g_ia_jk, n_ov, n_oo)
!
!     Read Cholesky vector L_ia_J
!
      call wfn % get_cholesky_ia(L_ia_J)
!
!     Calculate g_ia_jk
!
      call dgemm('N','T',    &
                  n_ov,      &
                  n_oo,      &
                  wfn % n_J, &
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
      call deallocator(L_ia_J, n_ov, wfn % n_J)
!
!
!     Allocation for g_ai_jk 
!
      call allocator(L_ai_J, n_ov, wfn % n_J)
      call allocator(g_ai_jk, n_ov, n_oo)
!
!     Get Cholesky AI vector
!
      call wfn % get_cholesky_ai(L_ai_J)
!
!     Calculate g_ai_jk
!
      call dgemm('N','T',    &
                  n_ov,      &
                  n_oo,      &
                  wfn % n_J, &
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
      call deallocator(L_ai_J, n_ov, wfn % n_J)
!
!     Add terms to Fock matrix
!
      do i = 1, wfn % n_occ
         do a = 1, wfn % n_vir
            do j = 1, wfn % n_occ
!
!              Needed indices
!
               ia = index_two(i, a, wfn % n_occ)
               ja = index_two(j, a, wfn % n_occ)
!
               ai = index_two(a, i, wfn % n_vir)
               aj = index_two(a, j, wfn % n_vir)
!
               jj = index_two(j, j, wfn % n_occ)
               ji = index_two(j, i, wfn % n_occ)
               ij = index_two(i, j, wfn % n_occ)
! 
!              Set the blocks of the Fock matrix
!
               fock_matrix(i, a + wfn % n_occ) = fock_matrix(i,a + wfn%n_occ) + & 
                                                 two*g_ia_jk(ia,jj) - g_ia_jk(ja,ij) ! g_ia_jk(ja,ij) = g_jaij = g_ijja
!
               fock_matrix(a + wfn % n_occ, i) = fock_matrix(a + wfn % n_occ, i) + & 
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
      n_vv = (wfn % n_vir)*(wfn % n_vir)
!
      call allocator(g_ab_ij, n_vv, n_oo)
      g_ab_ij = zero
!
!     Batch over index a
!
      available = get_available()
!
      required  = 2*(wfn % n_vir)*(wfn % n_vir)*(wfn % n_J)*4 + &
                  2*(wfn % n_vir)*(wfn % n_occ)*(wfn % n_J)*4
!
      call num_batch(required, available, max_batch_length, n_batches, wfn % n_vir)
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
         call batch_limits(batch_start, batch_end, b_batch, max_batch_length, wfn % n_vir)
         batch_length = batch_end - batch_start + 1
!
!        Allocate L_ab_J
!
         call allocator(L_ab_J, (wfn % n_vir)*batch_length, wfn % n_J)
         L_ab_J = zero
!
!        Read Cholesky vectors
!
         call wfn % get_cholesky_ab(L_ab_J, batch_start, batch_end, (wfn % n_vir)*batch_length, .false.)
!
!        Calculate g_ab_ij = sum_J L_ab_J*L_ij_J
!
         g_off = index_two(1, batch_start, wfn % n_vir)
!
         call dgemm('N','T',                     &
                     (wfn % n_vir)*batch_length, &
                     n_oo,                       &
                     wfn % n_J,                  &
                     one,                        &
                     L_ab_J,                     &
                     (wfn % n_vir)*batch_length, &
                     L_ij_J,                     &
                     n_oo,                       &
                     one,                        &
                     g_ab_ij(g_off,1),           &
                     n_vv)
!
!        Deallocate L_ab_J
!
         call deallocator(L_ab_J, batch_length*(wfn % n_vir), wfn % n_J)
!
      enddo ! batching done
!
!     Deallocate L_ij_J
!
      call deallocator(L_ij_J, n_oo, wfn % n_J)
!
!     Allocate for g_ai_jb
!
      call allocator(g_ai_jb, n_ov, n_ov)
      call allocator(L_ai_J, n_ov, wfn % n_J)
      call allocator(L_ia_J, n_ov, wfn % n_J)
!
!     Read Cholesky vectors L_ia_J and L_ai_J
!
      call wfn % get_cholesky_ai(L_ai_J)
!
      call dgemm('N','T',    &
                  n_ov,      &
                  n_ov,      &
                  wfn % n_J, &
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
     call deallocator(L_ia_J, n_ov, wfn % n_J)
     call deallocator(L_ai_J, n_ov, wfn % n_J)
!
!     Calculate two-electron terms for virtual-virtual blocks
!
      do a = 1, wfn % n_vir
         do b = 1, wfn % n_vir
!
            ab = index_two(a, b, wfn%n_vir)
!
            do i = 1, wfn % n_occ
!
               ii = index_two(i, i, wfn % n_occ)
               ai = index_two(a, i, wfn % n_vir)
               bi = index_two(b, i, wfn % n_vir)
               ia = index_two(i, a, wfn % n_occ)
               ib = index_two(i, b, wfn % n_occ)
!
               fock_matrix(wfn % n_occ + a, wfn % n_occ + b) = fock_matrix(wfn % n_occ + a, wfn % n_occ + b) &
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
      do i = 1, wfn % n_occ
         do j = 1, wfn % n_occ
!
            wfn % fock_matrix_ij(i,j) = fock_matrix(i,j)
!
         enddo
      enddo
!
      do i = 1, wfn % n_occ
         do a = 1, wfn % n_vir
!
            wfn % fock_matrix_ia(i,a) = fock_matrix(i, wfn % n_occ + a)
            wfn % fock_matrix_ai(a,i) = fock_matrix(wfn % n_occ + a, i)
!
         enddo
      enddo
!
      do a = 1, wfn % n_vir
         do b = 1, wfn % n_vir
!
            wfn % fock_matrix_ab(a,b) = fock_matrix(wfn % n_occ + a, wfn % n_occ + b)
!
         enddo
      enddo
!
      call deallocator(fock_matrix, wfn % n_mo, wfn % n_mo)
!
   end subroutine fock_constructor_cc_singles
!
!
   subroutine one_electron_t1_cc_singles(wfn,h1)
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
      class(cc_singles) :: wfn
!
      real(dp), dimension(wfn % n_mo, wfn % n_mo) :: h1
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
      call allocator(t1, wfn % n_mo, wfn % n_mo)
      t1 = zero
!
      call allocator(y, wfn % n_mo, wfn % n_mo)
      call allocator(x, wfn % n_mo, wfn % n_mo)
!
      y = zero
      x = zero
!
!     Set t1_p_q = t1am_p_q for p virtual and q occupied, 0 otherwise
!
      do a = 1, wfn % n_vir
         do i = 1, wfn % n_occ
!
            t1(wfn % n_occ + a, i) = wfn % t1am(a, i)
!
         enddo
      enddo
!
!     Form the x and y arrays
!
      do p = 1, wfn % n_mo
         do q = 1, wfn % n_mo
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
      call deallocator(t1, wfn % n_mo, wfn % n_mo)
!
!     Allocate Z intermediate
!   
      call allocator(Z, wfn % n_mo, wfn % n_mo)
!
!     Calculate h1_T1 = x*h1*y^T = x*Z
!
      call dgemm('N','T',     &
                  wfn % n_mo, &
                  wfn % n_mo, &
                  wfn % n_mo, &
                  one,        &
                  h1,         &
                  wfn % n_mo, &
                  y,          &
                  wfn % n_mo, &
                  zero,       &
                  Z,          &
                  wfn % n_mo)
!
      call dgemm('N','N',     &
                  wfn % n_mo, &
                  wfn % n_mo, &
                  wfn % n_mo, &
                  one,        &
                  x,          &
                  wfn % n_mo, &
                  Z,          &
                  wfn % n_mo, &
                  zero,       &
                  h1_T1,      &
                  wfn % n_mo)
!
!     Deallocate x and y, and the intermediate Z
!
      call deallocator(Z, wfn % n_mo, wfn % n_mo)
      call deallocator(y, wfn % n_mo, wfn % n_mo)
      call deallocator(x, wfn % n_mo, wfn % n_mo)
!
   end subroutine one_electron_t1_cc_singles
!
!
   subroutine get_cholesky_ij_cc_singles(wfn,L_ij_J)
!
!     Get Cholesky IJ
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Reads and T1-transforms the IA Cholesky vectors:
!     
!        L_ij_J_T1 = L_ij_J + sum_a t_aj * L_ia_J
!
!     Needed memory by routine: (n_occ*n_vir*n_J*2)
!
      implicit none 
!
      class(cc_singles) :: wfn
!
      real(dp), dimension((wfn % n_occ)**2, wfn % n_J) :: L_ij_J
!
      real(dp), dimension(:,:), allocatable :: L_ia_J ! L_ia^J
      real(dp), dimension(:,:), allocatable :: L_iJ_a ! L_ia^J reordered
      real(dp), dimension(:,:), allocatable :: L_iJ_k ! L_ik^J reordered
!
      integer(i15) :: i = 0, J = 0, a = 0, ij = 0, ia = 0, ik = 0, k = 0
!
      integer(i15) :: n_ov = 0
!
!     Calculate n_occ * n_vir = n_ov
!
      n_ov = (wfn % n_occ)*(wfn % n_vir)
!
!     Allocate
!
      call allocator(L_ia_J, n_ov, wfn % n_J)
      call allocator(L_iJ_a, (wfn % n_occ)*(wfn % n_J), wfn % n_vir)
!
      L_ia_J = zero
      L_iJ_a = zero
!
!     Read the untransformed Cholesky vectors 
!
      call wfn % read_cholesky_ij(L_ij_J)
      call wfn % read_cholesky_ia(L_ia_J)
!
!     Reorder L_ia_J to L_iJ_a 
!
      do i = 1, wfn % n_occ
         do J = 1, wfn % n_J
            do a = 1, wfn % n_vir
!              
!              Needed indices
!
               iJ = index_two(i, J, wfn % n_occ)
               ia = index_two(i, a, wfn % n_occ)
!
               L_iJ_a(iJ, a) = L_ia_J(ia, J)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_ia_J
!
      call deallocator(L_ia_J, n_ov, wfn % n_J)
!
!     Allocate L_iJ_k
!
      call allocator(L_iJ_k, (wfn % n_occ)*(wfn % n_J), wfn % n_occ)
      L_iJ_k = zero
!
!     T1-transformation
!
      call dgemm('N','N',                    &
                  (wfn % n_occ)*(wfn % n_J), &
                  wfn % n_occ,               &
                  wfn % n_vir,               &
                  one,                       &
                  L_iJ_a,                    &
                  (wfn % n_occ)*(wfn % n_J), &
                  wfn % t1am,                &
                  wfn % n_vir,               &
                  zero,                      &
                  L_iJ_k,                    &
                  (wfn % n_occ)*(wfn % n_J))
!
!     Place terms from L_iJ_k into L_ij_J
!
      do i = 1, wfn % n_occ
         do k = 1, wfn % n_occ
            do J = 1, wfn % n_J
!              
!              Needed indices
!
               iJ = index_two(i, J, wfn % n_occ)
               ik = index_two(i, k, wfn % n_occ)
!
               L_ij_J(ik, J) = L_ij_J(ik, J) + L_iJ_k(iJ, k)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_iJ_k and L_iJ_a
!
      call deallocator(L_iJ_k, (wfn % n_occ)*(wfn % n_J), wfn % n_occ)
      call deallocator(L_iJ_a, (wfn % n_occ)*(wfn % n_J), wfn % n_vir)
!    
   end subroutine get_cholesky_ij_cc_singles
!
   subroutine get_cholesky_ia_cc_singles(wfn,L_ia_J)
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
      class(cc_singles) :: wfn
!
      real(dp), dimension((wfn % n_occ)*(wfn % n_vir), wfn % n_J) :: L_ia_J
!
      call wfn % read_cholesky_ia(L_ia_J)
!
   end subroutine get_cholesky_ia_cc_singles
!
!
   subroutine get_cholesky_ai_cc_singles(wfn,L_ai_J)
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
!     Required memory by routine: 2*n_vir*n_J*batch_length 
!                                 + n_vir*n_occ*n_J
!
      implicit none 
!
      class(cc_singles) :: wfn
!
      real(dp), dimension((wfn % n_vir)*(wfn % n_occ), wfn % n_J) :: L_ai_J 
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
      call wfn % read_cholesky_ai(L_ai_J)
!
!                       
!     -::- L_ab_J contributions -::-
!
!
!     Allocate L_Ja_i
!
      call allocator(L_Ja_i, (wfn % n_J)*(wfn % n_vir), wfn % n_occ)
      L_Ja_i = zero
!
!     Set batching variables 
!
      required = 2*(wfn % n_vir)**2*(wfn % n_J)*4
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
      call num_batch(required, available, max_batch_length, n_batch, wfn % n_vir)
!
      do a_batch = 1, n_batch
!
!        Get start, end, and length of batch
!
         call batch_limits(batch_start, batch_end, a_batch, max_batch_length, wfn % n_vir)
         batch_length = batch_end - batch_start + 1
!
!        Allocate L_ab_J and L_Ja_b
!
         call allocator(L_ba_J, (wfn % n_vir)*batch_length, wfn % n_J) ! L_ab^J = L_ba_J(ba,J)
         call allocator(L_Ja_b, batch_length*(wfn % n_J), wfn % n_vir)
!
         L_ba_J = zero
         L_Ja_b = zero 
!
!        Read Cholesky AB vectors, batching over a
! 
         reorder = .true.
         call wfn % read_cholesky_ab(L_ba_J, batch_start, batch_end,&
                                     (wfn % n_vir)*batch_length, reorder)
!
!        Reorder the Cholesky array L_ba_J
!
         do a = 1, batch_length
            do b = 1, wfn % n_vir
               do J = 1, wfn % n_J
!
!                 Needed indices
!
                  ba = index_two(b, a, wfn % n_vir)
                  Ja = index_two(J, a, wfn % n_J)
!
                  L_Ja_b(Ja, b) = L_ba_J(ba, J) ! L_ab^J
!
               enddo
            enddo
         enddo         
!
!        Calculate sum_b L_Ja_b*t_b_i = L_Ja_i 
!        
         L_off = index_two(1, batch_start, wfn % n_J)
!
         call dgemm('N','N',                   &
                     batch_length*(wfn % n_J), &
                     wfn % n_occ,              &
                     wfn % n_vir,              &
                     one,                      &
                     L_Ja_b,                   &
                     batch_length*(wfn % n_J), &
                     wfn % t1am,               &
                     wfn % n_vir,              &
                     one,                      &
                     L_Ja_i(L_off, 1),         &
                     (wfn % n_vir)*(wfn % n_J))
!
!        Deallocate L_ab_J and L_Ja_b
!
         call deallocator(L_ba_J, (wfn % n_vir)*batch_length, wfn % n_J)
         call deallocator(L_Ja_b, batch_length*(wfn % n_J), wfn % n_vir)
!
      enddo ! batching over a 
!
!     Add terms to T1-transformed Cholesky AI vector 
!
      do i = 1, wfn % n_occ
         do a = 1, wfn % n_vir
            do J=1, wfn % n_J
!
!              Needed indices
!
               Ja = index_two(J, a, wfn % n_J)
               ai = index_two(a, i, wfn % n_vir)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_Ja_i(Ja, i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_Ja_i
!
      call deallocator(L_Ja_i, (wfn % n_J)*(wfn % n_vir), wfn % n_occ)
!
!
!     -::- L_ij_J contributions -::-
!
!
!     Allocate L_a_iJ, L_ik_J, L_k_iJ
!
      call allocator(L_a_iJ, wfn % n_vir, (wfn % n_J)*(wfn % n_occ))
      call allocator(L_ik_J, (wfn % n_occ)**2, wfn % n_J)
      call allocator(L_k_iJ, wfn % n_occ, (wfn % n_occ)*(wfn % n_J))
!
      L_a_iJ = zero   
      L_ik_J = zero 
      L_k_iJ = zero
!  
!     Read Cholesky IJ vectors
!
      call wfn % read_cholesky_ij(L_ik_J) ! L_ik_J(ik,J) = L_ik^J 
!
!     Reorder IJ Cholesky vectors
!
      do i = 1, wfn % n_occ
         do k = 1, wfn % n_occ
            do J = 1, wfn % n_J
!
!              Needed indices
!
               ik = index_two(i, k, wfn % n_occ)
               iJ = index_two(i, J, wfn % n_occ)
!
               L_k_iJ(k, iJ) = L_ik_J(ik, J) ! L_k_iJ(k,iJ) = L_ik^J
!
            enddo
         enddo
      enddo
!
!     Calculate -sum_k t_a_k*L_k_iJ = L_a_iJ  ! Here we assume L_ik^J = L_ki^J - this is safe, right?
!
      call dgemm('N','N',                    &
                  wfn % n_vir,               &
                  (wfn % n_occ)*(wfn % n_J), &
                  wfn % n_occ                &
                  -one,                      &
                  wfn % t1am,                &
                  wfn % n_vir,               &
                  L_k_iJ,                    &
                  wfn % n_occ,               &
                  zero,                      &
                  L_a_iJ,                    &
                  wfn % n_vir)
!
!     Add terms to T1-transformation of L_ai_J
!
      do i = 1, wfn % n_occ
         do a = 1, wfn % n_vir
            do J = 1, wfn % n_J
!
!              Needed indices
!
               ai = index_two(a, i, wfn % n_vir)
               iJ = index_two(i, J, wfn % n_occ)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_a_iJ, L_ik_J, L_k_iJ
!
      call deallocator(L_a_iJ, wfn % n_vir, (wfn % n_J)*(wfn % n_occ))      
      call deallocator(L_ik_J, (wfn % n_occ)**2, wfn % n_J)
      call deallocator(L_k_iJ, wfn % n_occ, (wfn % n_occ)*(wfn % n_J))
!
!
!     -::- L_jb_J contributions -::-    
!
!
      call allocator(L_kJ_b, (wfn % n_occ)*(wfn % n_J), wfn % n_vir)
      call allocator(L_kb_J, (wfn % n_occ)*(wfn % n_vir), wfn % n_J)
!
!     Read the Cholesky vector L_kb_J
!
      call wfn % read_cholesky_ia(L_kb_J)
!
!     Reorder L_kb_J to L_kJ_b
!
      do k = 1, wfn % n_occ
         do b = 1, wfn % n_vir
            do J = 1, wfn % n_J
!
!              Needed indices
!
               kb = index_two(k, b, wfn % n_occ)
               kJ = index_two(k, J, wfn % n_occ)
!
               L_kJ_b(kJ, b) = L_kb_J(kb, J)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_kb_J
!
      call deallocator(L_kb_J, (wfn % n_occ)*(wfn % n_vir), wfn % n_J)
!
!     Allocate L_kJ_i 
!
      call allocator(L_kJ_i, (wfn % n_occ)*(wfn % n_J), wfn % n_occ)
      L_kJ_i = zero
!
!     Calculate sum_b L_kJ_b*t_b_i = L_kJ_i
!
      call dgemm('N','N',                    &
                  (wfn % n_occ)*(wfn % n_J), &
                  wfn % n_occ,               &
                  wfn % n_vir,               &
                  one,                       &
                  L_kJ_b,                    &
                  (wfn % n_occ)*(wfn % n_J), &
                  wfn % t1am,                &
                  wfn % n_vir,               &
                  zero,                      &
                  L_kJ_i,                    &
                  (wfn % n_occ)*(wfn % n_J))
!
!     Deallocate L_kJ_b
!
      call deallocator(L_kJ_b, (wfn % n_occ)*(wfn % n_J), wfn % n_vir)
!
!     Allocate L_k_iJ
!  
      call allocator(L_k_iJ, (wfn % n_occ), (wfn % n_occ)*(wfn % n_J))
      L_k_iJ = zero
!
!     Reorder L_kJ_i to L_k_iJ    
!
      do i = 1, wfn % n_occ
         do k = 1, wfn % n_occ
            do J = 1, wfn % n_J
!
!              Needed indices
!
               kJ = index_two(k, J, wfn % n_occ)
               iJ = index_two(i, J, wfn % n_occ)
!
               L_k_iJ(k, iJ) = L_kJ_i(kJ, i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_kJ_i
!
      call deallocator(L_kJ_i, (wfn % n_occ)*(wfn % n_J), wfn % n_occ)
!
!     Allocate L_a_iJ
!
      call allocator(L_a_iJ, wfn % n_vir, (wfn % n_occ)*(wfn % n_J))
      L_a_iJ = zero
!      
!     Calculate sum_k t_a_k*L_k_iJ = L_a_iJ
!
      call dgemm('N','N',                    &
                  wfn % n_vir,               &
                  (wfn % n_occ)*(wfn % n_J), &
                  wfn % n_occ,               &
                  -one,                      &
                  wfn % t1am,                &
                  wfn % n_vir,               &
                  L_k_iJ,                    &
                  wfn % n_occ,               &
                  zero,                      &
                  L_a_iJ,                    &
                  wfn % n_vir)
!
!     Add contribution to L ai_J
!
      do a = 1, wfn % n_vir
         do i = 1, wfn % n_occ
            do J = 1, wfn % n_J
!
!              Needed indices
!
               iJ = index_two(i, J, wfn % n_occ)
               ai = index_two(a, i, wfn % n_vir)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
            enddo
         enddo
      enddo
!
!     Deallocations
!
      call deallocator(L_a_iJ, wfn % n_vir, (wfn % n_occ)*(wfn % n_J))
      call deallocator(L_k_iJ, wfn % n_occ, (wfn % n_occ)*(wfn % n_J))
!
   end subroutine get_cholesky_ai_cc_singles
!
!
   subroutine get_cholesky_ab_cc_singles(wfn,L_ab_J,first,last,ab_dim,reorder)
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
!     Required memory: n_J*batch_length*n_vir + n_vir*n_occ*n_J*2
!
      implicit none
!
      class(cc_singles) :: wfn
!
      integer(i15), intent(in) :: ab_dim
      integer(i15), intent(in) :: first
      integer(i15), intent(in) :: last
!
      logical, intent(in) :: reorder
!
      real(dp), dimension(ab_dim, wfn % n_J) :: L_ab_J
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
         call allocator(L_ib_J, (wfn % n_occ)*(wfn % n_vir), wfn % n_J)
         L_ib_J = zero
!
!        Read L_ia_J
!  
         call wfn % read_cholesky_ia(L_ib_J) ! Note: using L_ia_J instead of L_ai_J, here, to avoid two reorderings.
                                             ! This is possible because of the symmetry L_ai_J(ai,J) == L_ia_J(ia,J).
!
!        Read L_ab_J for batch of a
!
         call wfn % read_cholesky_ab(L_ab_J, first, last, ab_dim, reorder) ! Eirik: L_ab_J(ba) = L_ab^J 
!
!        Allocate L_i,Jb
!
         call allocator(L_i_Jb, wfn % n_occ, (wfn % n_J)*(wfn % n_vir))
         L_i_Jb = zero
!
!        Reorder L_ib_J to L_i_Jb
!
         do i = 1, wfn % n_occ
            do b = 1, wfn % n_vir
               do J = 1, wfn % n_J
!
!                 Needed indices
!
                  ib = index_two(i, b, wfn % n_occ)
                  Jb = index_two(J, b, wfn % n_J)
!
                  L_i_Jb(i, Jb)=L_ib_J(ib, J)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator(L_ib_J, (wfn % n_occ)*(wfn % n_vir), wfn % n_J)
!
!        Allocate L_a_Jb for batch of a
!  
         call allocator(L_a_Jb, batch_length, (wfn % n_vir)*(wfn % n_J))
!
!        Calculate  -t1_a_i * L_i_Jb = L_a_Jb
!
         call dgemm('N','N',                    &
                     batch_length,              &
                     (wfn % n_vir)*(wfn % n_J), &
                     wfn % n_occ,               &
                     -one,                      &
                     wfn % t1am(first,1),       &
                     wfn % n_vir,               &
                     L_i_Jb,                    &
                     wfn % n_occ,               &
                     zero,                      &
                     L_a_Jb,                    &
                     batch_length)
!
!        Add terms of L_a_Jb to L_ab_J
!
         do a = 1, batch_length
            do b = 1, wfn % n_vir
               do J = 1, wfn % n_J
!
!                 Needed indices
!
                  Jb = index_two(J, b, wfn % n_J)
                  ba = index_two(b, a, wfn % n_vir)
!
                  L_ab_J(ba, J) = L_ab_J(ba, J) + L_a_Jb(a, Jb) 
!
               enddo
            enddo
         enddo
!
!        Dellocate L_i_Jb and L_a_Jb for batch of a
!
         call deallocator(L_a_Jb, batch_length, (wfn % n_J)*(wfn % n_vir))
         call deallocator(L_i_Jb, (wfn % n_J)*(wfn % n_vir),(wfn % n_occ))
!
      else  !! Batching over b !!
!
!        Allocate L_ib_J
!     
         call allocator(L_ib_J, (wfn % n_occ)*(wfn % n_vir), wfn % n_J)
         L_ib_J = zero
!
!        Read L_ia_J
!  
         call wfn % read_cholesky_ia(L_ib_J) ! Note: using L_ia_J instead of L_ai_J, here, to avoid two reorderings.
                                             ! This is possible because of the symmetry L_ai_J(ai,J) == L_ia_J(ia,J)
!
!        Read L_ab_J for batch of b
!
         call wfn % read_cholesky_ab(L_ab_J, first, last, ab_dim, reorder)
!
!        Allocate L_Jb,i for batch of b
!
         call allocator(L_Jb_i, (wfn % n_J)*batch_length, wfn % n_occ)
         L_Jb_i = zero
!
!        Reorder L_ib_J to L_Jb_i
!
         do i = 1, wfn % n_occ
            do b = 1, batch_length
               do J = 1, wfn % n_J
!
!                 Needed indices
!
                  ib = index_two(i, b + first - 1, wfn % n_occ) ! Note: in L_ib_J we have all b, not the case for L_Jb_i
                  Jb = index_two(J, b, wfn % n_J)
!
                  L_Jb_i(Jb, i) = L_ib_J(ib, J)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_ib_J
!  
         call deallocator(L_ib_J, (wfn % n_occ)*(wfn % n_vir), wfn % n_J)
!
!        Allocate L_Jb_a for batch of b
!  
         call allocator(L_Jb_a, (wfn % n_J)*batch_length, wfn % n_vir)
!
!        T1-transformation
!
         call dgemm('N','T',&
                     (wfn % n_J)*batch_length, &
                     wfn % n_vir,              &
                     wfn % n_occ,              &
                     -one,                     &
                     L_Jb_i,                   &
                     (wfn % n_J)*batch_length, &
                     wfn % t1am,               &
                     wfn % n_vir,              &
                     zero,                     &
                     L_Jb_a,                   &
                     batch_length*(wfn % n_J))
!
!        Add terms of L_Jb_a to L_ab_J
!
         do a = 1, wfn % n_vir
            do b = 1, batch_length
               do J = 1, wfn % n_J
!
!                 Needed indices
!
                  Jb = index_two(J, b, wfn % n_J)
                  ab = index_two(a, b, wfn % n_vir)
!
                  L_ab_J(ab, J) = L_ab_J(ab, J) + L_Jb_a(Jb, a)
!
               enddo
            enddo
         enddo
!
!        Dellocate L_Jb,i and L_Jb_a for batch of b
!
         call deallocator(L_Jb_a, (wfn % n_J)*batch_length, wfn % n_vir)
         call deallocator(L_Jb_i, (wfn % n_J)*batch_length, wfn % n_occ)
!
      endif
!
   end subroutine get_cholesky_ab_cc_singles
!
!
end module ccs_class
