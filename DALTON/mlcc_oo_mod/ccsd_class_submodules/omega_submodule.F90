submodule (ccsd_class) omega
!
!
!                       -::- Cholesky submodule (CCS) -::-
!           Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!
!     Contains the following family of procedures of the CCSD class:
!
!     construct_omega(wfn):         constructs the projection vector for the current amplitudes
!     omega_c1_cc_singles_doubles:  C1 contribution to omega  
!
contains
!
!
      subroutine construct_omega_cc_singles_doubles(wf)
!
         implicit none 
!
         class(cc_singles_doubles) :: wf 
!
         ! call omega_a1
         ! call omega_b1 ...
!
      end subroutine construct_omega_cc_singles_doubles
!
!
      subroutine omega_c1_cc_singles_doubles(wf)        
!
!        C1 omega term: Omega_ai^C1 = sum_ck F_kc*u_ai_ck,
!                       u_ai_ck = 2*t_ck_ai-t_ci_ak
!        
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mars 2017
!
!      
         implicit none
!
         class(cc_singles_doubles) :: wf 
!
         real(dp),dimension(:,:),allocatable  :: F_ck      => null()
         real(dp),dimension(:,:),allocatable  :: u_ai_ck   => null()
         real(dp),dimension(:,:),allocatable  :: omega1_ai => null()
!
         integer(i15) :: i=0,k=0,c=0,a=0
         integer(i15) :: ck=0,ai=0,ak=0,ci=0,aick=0,akci=0
!
         logical :: debug = .false.
!
!        Allocation of F_ck, u_ai_ck and omega1_ai
!
         call allocator(F_ck,n_ov,1)
         call allocator(u_ai_ck,n_ov,n_ov)  
         call allocator(omega1_ai,n_ov,1)
!
         F_ck      = zero
         u_ai_ck   = zero
         omega1_ai = zero
!
!        Set up u_ck_ai and virtual-occupied Fock matrix
!
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
!              Set up compound index
!  
               ck = index_two(c, k, wf%n_v)
!
!              Reorder MO Fock matrix
!
               F_ck(ck, 1) = F_i_a(k, c)
!
               do a = 1, wf%n_v
                  do i = 1, wf%n_o
!
!                    Set up compound indices
!
                     ai = index_two(a, i, wf%n_v)
                     ci = index_two(c, i, wf%n_v)
                     ak = index_two(a, k, wf%n_v)
!
                     aick = index_packed(ck, ai)
                     akci = index_packed(ci, ak)
!                    
!                    u_ck_ai
!
                     u_ai_ck(ai, ck) = two*t2am(aick, 1) - t2am(akci, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Matrix multiplication
!
         call dgemm('N','N',            &
                     (wf%n_o)*(wf%n_v), &
                     1,                 &
                     (wf%n_o)*(wf%n_v), &
                     one,               &
                     u_ai_ck,           &
                     (wf%n_o)*(wf%n_v), &
                     F_ck,              &
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     omega1_ai,         &
                     (wf%n_o)*(wf%n_v))
!
!
         do a = 1, wf%n_v
            do i = 1, wf%n_o
!
!              Set up compound index
!
               ai = index_two(a, i, wf%n_v)
!
!              Add to omega1
!
               omega1(a,i)=omega1(a,i)+omega1_ai(ai,1)
!
            enddo
         enddo
!
!        Deallocation
!
         call deallocator(F_ck,n_ov,1)
         call deallocator(omega1_ai,n_ov,1)
         call deallocator(u_ai_ck,n_ov,n_ov)
!
   end subroutine omega_c1_cc_singles_doubles
!
!
   subroutine mlcc_omega_d1(wf)
!
!     D1 omega term: Omega_ai^D1=F_ai_T1
!
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mars 2017
!
      implicit none
!
      class(cc_singles_doubles) :: wf
!
!     Add F_a_i to omega
!
      call daxpy((wf%n_o)*(wf%n_v),one,F_a_i,1,omega1,1)
!
   end subroutine mlcc_omega_d1
!
!
   subroutine omega_a2_cc_singles_doubles(wf)
!
!     MLCC Omega A2 term: Omega A2 = g_ai_bj + sum_(cd)g_ac_bd * t_ci_dj = A2.1 + A.2.2
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 10 Mar 2017
!
!     Structure: Batching over both a and b. If no batching is necessary L_ab_J is only read once, and g_ac_bd 
!                is constructed and kept in memory full size. 
!                g_ac_bd is reordered as g_ab_cd and t_ci_dj is reordered as t_cd_ij.
!                Omega contribution for A2.2 is ordered as Omega_ab_ij, and is reordered into the packed omega2 vector.          
!
      implicit none
!
      class(cc_singles_doubles) :: wf
!
!     Integrals
!
      real(dp),dimension(:,:),allocatable :: g_ai_bj 
      real(dp),dimension(:,:),allocatable :: g_ca_db 
      real(dp),dimension(:,:),allocatable :: g_ab_cd
      real(dp),dimension(:,:),allocatable :: L_ai_J 
      real(dp),dimension(:,:),allocatable :: L_ca_J 
      real(dp),dimension(:,:),allocatable :: L_db_J 
!
!     Reordered T2 amplitudes
!
      real(dp),dimension(:,:),allocatable :: t_cd_ij
!
!     Reordered omega 2 
!  
      real(dp),dimension(:,:),allocatable :: omega2_ab_ij
! 
!     Indices
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: ab = 0, ca = 0, cd = 0, db = 0 
      integer(i15) :: ai = 0, bj = 0, ci = 0, dj = 0
      integer(i15) :: ij = 0
!
      integer(i15) :: aibj = 0, cidj = 0
!
!     Batching and memory handling variables
!
      integer(i15) :: a_n_batch = 0, a_start = 0, a_end = 0, a_length = 0, a_max_length = 0, a_batch = 0
      integer(i15) :: b_n_batch = 0, b_start = 0, b_end = 0, b_length = 0, b_max_length = 0, b_batch = 0

      integer(i15) :: required = 0, available = 0
!
!     Debug 
!
      logical :: debug = .false.
!
!!!   A2.1 term   !!!
!
!
!     Create g_ai_bj
!  
      call allocator(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call allocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call get_cholesky_ai(L_ai_J)
!
!     g_ai_bj = sum_J L_ai_J*L_bj_J
!     
      call dgemm('N','T',   &
         (wf%n_o)*(wf%n_v), &
         (wf%n_o)*(wf%n_v), &
         wf%n_J,            &
         one,               & 
         L_ai_J,            &
         (wf%n_o)*(wf%n_v), &
         L_ai_J,            &
         (wf%n_o)*(wf%n_v), &
         zero,              &
         g_ai_bj,           &
         (wf%n_o)*(wf%n_v))
!
!
      call deallocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Add A2.1 to Omega 2
!     
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Set up compound indices
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  if(ai .ge. bj) then
!
                     aibj = index_packed(ai, bj)
!
!                    Add to omega
!
                     omega2(aibj, 1) = omega2(aibj, 1) + g_ai_bj(ai, bj)
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
!     Determine number of a batches
!
      required = 2*(wf%n_v)*(wf%n_v)*(wf%n_J)*4 + 2*(wf%n_o)*(wf%n_v)*(wf%n_J)*4 + 2*(wf%n_v)*(wf%n_v)*(wf%n_v)*(wf%n_v)*4 ! Needed to get cholesky of ab-type and for g_ab_bd
      available=get_available()
!
      a_max_length=0
!
      call num_batch(required, available, a_max_length, a_n_batch, wf%n_v)
!
!     Initialize some variables for batching
!
      a_start  =0
      a_end    =0
      a_length =0
!
!     Start looping over a-batches
!
      do a_batch = 1, a_n_batch
!
         call one_batch_limits(a_start, a_end, a_batch, a_max_length, wf%n_v)
         a_length = a_end - a_start + 1
!
!        Get cholesky vectors L_ac^J ordered as L_ca_J
!
         call allocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
!
         call get_cholesky_ab(L_ca_J, a_start, a_end, (wf%n_v)*a_length, .true.) 
!        
!        Prepare for batching over b
!
!        Determine number of a batches
!
         required = 2*(wf%n_v)*(wf%n_v)*(wf%n_J)*4 + 2*(wf%n_o)*(wf%n_v)*(wf%n_J)*4 &  ! Needed to get cholesky of ab-type 
                  + 2*(wf%n_v)*(wf%n_v)*a_length*(wf%n_v)*4                            ! and to generate g_ac_bd.
                                                                          
         available = get_available()
!  
         b_max_length = 0
!  
         call num_batch(required, available, b_max_length, b_n_batch, wf%n_v)
!  

         if (a_n_batch .eq. 1 .and. b_n_batch .eq. 1 ) then ! No need to read ab-type cholesky twice. 
                                                            ! Note that this means that a_length=b_length=n_vir
!
!           Allocate full g_ca_db because we have room for it!
!
            call allocator(g_ca_db, (wf%n_v)*(wf%n_v),(wf%n_v)*(wf%n_v))
!
            call dgemm('N','T',   &
               (wf%n_v)*(wf%n_v), &
               (wf%n_v)*(wf%n_v), &
               wf%n_J,            &
               one,               &
               L_ca_J,            &
               (wf%n_v)*(wf%n_v), &
               L_ca_J,            &
               (wf%n_v)*(wf%n_v), &
               zero,              &
               g_ca_db,           &
               (wf%n_v)*(wf%n_v))
!
!           Allocate for reordered g_ca_db(g_ab_cd) and t_ci_dj (t_cd_ij)
!
            call allocator(g_ab_cd, (wf%n_v)*(wf%n_v), (wf%n_v)*(wf%n_v))
            call allocator(t_cd_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
            do c = 1, wf%n_v
               do d = 1, wf%n_v
                  do a = 1, wf%n_v
                     do b = 1, wf%n_v
!
!                       Calculate compound indices
!  
                        ca = index_two(c, a, wf%n_v)
                        db = index_two(d, b, wf%n_v)
                        ab = index_two(a, b, wf%n_v)
                        cd = index_two(c, d, wf%n_v)
!
!                 Reorder g_ca_db to g_ab_cd
!
                        g_ab_cd(ab,cd) = g_ca_db(ca,db)
                     enddo
                  enddo
!
                  do i = 1,n_occ
                     do j = 1,n_occ
!
!                       Calculate compound indices
!  
                        cd = index_two(c, d, wf%n_v)
                        ij = index_two(i, j, wf%n_o)
                        ci = index_two(c, i, wf%n_v)
                        dj = index_two(d, j, wf%n_v)
!
                        cidj = index_packed(ci, dj)
!
!                       Reorder t_ci_dj to t_cd_ij
!
                        t_cd_ij(cd, ij) = t2am(cidj, 1)
!
                     enddo
                  enddo
               enddo
            enddo
!
!           Deallocate g_ca_db
!
            call deallocator(g_ca_db, (wf%n_v)*(wf%n_v), (wf%n_v)*(wf%n_v))
!
!           Allocate omega2_ab_ij
!
            call allocator(omega2_ab_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
!           omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
!
            call dgemm('N','N',            &
                        (wf%n_v)*(wf%n_v), &
                        (wf%n_o)*(wf%n_o), &
                        (wf%n_v)*(wf%n_v), &
                        one,               &
                        g_ab_cd,           &
                        (wf%n_v)*(wf%n_v), &
                        t_cd_ij,           &
                        (wf%n_v)*(wf%n_v), &
                        zero,              &
                        omega2_ab_ij,      &
                        (wf%n_v)*(wf%n_v))
!
!           Deallocate t_cd_ij and g_ab_cd
!
            call deallocator(g_ab_cd, (wf%n_v)*(wf%n_v), (wf%n_v)*(wf%n_v))
            call deallocator(t_cd_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
               do i = 1, wf%n_o
                  do j = 1, wf%n_o
                     do a = 1, wf%n_v
                        do b = 1, wf%n_v
!
!                          Calculate compound indices
!  
                           ai = index_two(a, i, wf%n_v)
                           bj = index_two(b, j, wf%n_v)
!
                           if (ai .ge. bj) then
!
                              ab = index_two(a, b, wf%n_v)
                              ij = index_two(i, j, wf%n_o)
!  
                              aibj = index_packed(ai, bj)
!
!                             Reorder into omega2
!
                              omega2(aibj, 1) = omega2(aibj, 1) + omega2_ab_ij(ab, ij)
!  
                           endif
                        enddo
                     enddo
                  enddo
               enddo
!
!           Deallocate reordered omega 2
!
            call deallocator(omega2_ab_ij,n_vv,n_oo)
!
         else
!
!           Initialize some variables for batching over b
!
            b_start  =0
            b_end    =0
            b_length =0  
!
!           Start looping over batches of b
!
            do b_batch = 1, b_n_batch
!
               call one_batch_limits(b_start, b_end, b_batch, b_max_length, wf%n_v)
               b_length=b_end-b_start+1
!
!              Get cholesky vectors L_bd^J ordered as L_db_J
!
               call allocator(L_db_J, (wf%n_v)*b_length, wf%n_J)
!  
               call get_cholesky_ab(L_db_J, b_start, b_end, (wf%n_v)*b_length, .true.) 
!
!              Allocate g_ac_bd for batches of a and b, ordered as g_ca_db
!
               call allocator(g_ca_db, (wf%n_v)*a_length, (wf%n_v)*b_length)
!
!              g_ca_db = sum_J L_ca_J*L_db_J
!  
               call dgemm('N','T',            &
                           (wf%n_v)*a_length, &
                           (wf%n_v)*b_length, &
                           wf%n_J,            &
                           one,               &
                           L_ca_J,            &
                           (wf%n_v)*a_length, &
                           L_db_J,            &
                           (wf%n_v)*b_length, &
                           zero,              &
                           g_ca_db,           &
                           (wf%n_v)*a_length)
!
!              Deallocate L_db_J       
!
               call deallocator(L_db_J, (wf%n_v)*b_length, wf%n_J)
!
!              Reorder g_ca_db into g_ab_cd and t_ci_dj into t_cd_ij
!
               call allocator(g_ab_cd, a_length*b_length, (wf%n_v)*(wf%n_v))
               call allocator(t_cd_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
               do c = 1, wf%n_v 
                  do d = 1, wf%n_v
                     do a = 1, a_length
                        do b = 1, b_length
!
!                          Calculate compound indices
!
                           ca = index_two(c, a, wf%n_v)
                           db = index_two(d, b, wf%n_v)
                           ab = index_two(a, b, wf%n_v)
                           cd = index_two(c, d, wf%n_v)
!
!                          Reorder g_ca_db to g_ab_cd
!  
                           g_ab_cd(ab,cd) = g_ca_db(ca,db)
!
                        enddo
                     enddo
!
                  do i = 1, wf%n_o
                     do j = 1, wf%n_o
!
!                          Calculate compound indices
!     
                           cd = index_two(c, d, wf%n_v)
                           ij = index_two(i, j, wf%n_o)
                           ci = index_two(c, i, wf%n_v)
                           dj = index_two(d, j, wf%n_v)
!  
                           cidj = index_packed(ci, dj)
!
!                          Reorder t_ci_dj to t_cd_ij
!  
                           t_cd_ij(cd, ij) = t2am(cidj, 1)
!   
                        enddo
                     enddo
                  enddo
               enddo
!
               call deallocator(g_ca_db, (wf%n_v)*a_length, (wf%n_v)*b_length)
!
!
               call allocator(omega2_ab_ij, a_length*b_length, (wf%n_o)*(wf%n_o))
!  
!              omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
!  
               call dgemm('N','N',            & 
                           a_length*b_length, &
                           (wf%n_o)*(wf%n_o), &
                           (wf%n_v)*(wf%n_v), &
                           one,               &
                           g_ab_cd,           &
                           a_length*b_length, &
                           t_cd_ij,           &
                           (wf%n_v)*(wf%n_v), &
                           zero,              &
                           omega2_ab_ij,      &
                           a_length*b_length)
!  
!              Deallocate t_cd_ij and g_ab_cd
!
               call deallocator(t_cd_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
               call deallocator(g_ab_cd, a_length*b_length, (wf%n_v)*(wf%n_v))
!
               do i = 1, wf%n_o
                  do j = 1, wf%n_o
                     do a = 1, a_length
                        do b = 1, b_length
!  
!                          Calculate compound indices
!     
                           Ai = index_two(a + a_start - 1, i, wf%n_v) ! A is full-space a index
                           Bj = index_two(b + b_start - 1, j, wf%n_v) ! B is full-space b index
!
                           if (Ai .ge. Bj) then
!
                              ab = index_two(a, b, a_length)
                              ij = index_two(i, j, wf%n_o)
!     
                              AiBj = index_packed(Ai, Bj)
!
!                             Reorder into omega2_aibj
!  
                              omega2(AiBj,1) = omega2(AiBj, 1) + omega2_ab_ij(ab, ij)
!     
                           endif
!
                        enddo
                     enddo
                  enddo
               enddo
!
!              Deallocate reorderer omega
!
               call deallocator(omega2_ab_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
            enddo
!
         endif
!
!        Deallocate L_ca_J
!
         call deallocator(L_ca_J, (wf%n_v)*(wf%n_v), wf%n_J)
!
      enddo
!
!     Print the omega vector, having added A2
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after A2 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,wf%n_t2am)
      endif
!
   end subroutine omega_a2_cc_singles_doubles
!
!
!
   subroutine omega_b2_cc_singles_doubles
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
!
!     Integrals
!
      real(dp),dimension(:,:),allocatable :: L_kc_J     
      real(dp),dimension(:,:),allocatable :: L_ij_J  
      real(dp),dimension(:,:),allocatable :: g_kc_ld    
      real(dp),dimension(:,:),allocatable :: g_kl_cd    
      real(dp),dimension(:,:),allocatable :: g_kl_ij    
      real(dp),dimension(:,:),allocatable :: g_ki_lj 
!
!     Reordered T2 apmlitudes
!   
      real(dp),dimension(:,:),allocatable :: t_cd_ij    
      real(dp),dimension(:,:),allocatable :: t_ab_kl   
!
!     Intermediate for matrix multiplication
! 
      real(dp),dimension(:,:),allocatable :: X_kl_ij 
!
!     Reordered omega
!   
      real(dp),dimension(:,:),allocatable :: omega_ab_ij
!
!     Indices
!   
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ab = 0, cd = 0
      integer(i15) :: ai = 0, ak = 0, bj = 0, bl = 0, ci = 0, dj = 0
      integer(i15) :: kc = 0, ld = 0
      integer(i15) :: ij = 0, ki = 0, kl = 0, lj = 0
!
      integer(i15) :: aibj = 0, akbl = 0, cidj = 0 
!
      logical :: debug = .false.
!
!     Read cholesky vector of type L_ij_J
!
      call allocator(L_ij_J, (wf%n_o)*(wf%n_o), wf%n_J)
      L_ij_J = zero 
!
      call get_cholesky_ij(L_ij_J)
!
!     Create g_ki_lj = sum_J L_li_J*L_lj_J
!
      call allocator(g_ki_lj, (wf%n_o)*(wf%n_o), (wf%n_o)*(wf%n_o))
      g_ki_lj = zero 
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_o), &
                  (wf%n_o)*(wf%n_o), &
                  wf%n_J,            &
                  one,               &
                  L_ij_J,            &
                  (wf%n_o)*(wf%n_o), &
                  L_ij_J,            &
                  (wf%n_o)*(wf%n_o), &
                  zero,              &
                  g_ki_lj,           &
                  (wf%n_o)*(wf%n_o))
!
!
      call deallocator(L_ij_J, (wf%n_o)*(wf%n_o), wf%n_J)
!
      call allocator(g_kl_ij, (wf%n_o)*(wf%n_o),(wf%n_o)*(wf%n_o))
      g_kl_ij = zero
!
      do k = 1, wf%n_o
         do l = 1, wf%n_o
            do i = 1, wf%n_o
               do j=1, wf%n_o
!
!                 Calculate compound indices
!
                  ki = index_two(k, i, wf%n_o)
                  lj = index_two(l, j, wf%n_o)
                  kl = index_two(k, l, wf%n_o)
                  ij = index_two(i, j, wf%n_o)
!
!                  Reordering g_ki_lj to g_kl_ij
!
                  g_kl_ij(kl, ij) = g_ki_lj(ki, lj)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ki_lj, (wf%n_o)*(wf%n_o), (wf%n_o)*(wf%n_o))
!
!     Read cholesky vectors of ia-type into L_kc_J
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call get_cholesky_ia(L_kc_J)
!
!     Create g_ck_ld = sum_(J) L_kc_J*L_ld_J
!
      call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!  
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate cholesky vectors L_ck_J
!
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Reorder g_kc_ld as g_kl_cd, also reordering t_ci_dj as t_cd_ij
!
      call allocator(t_cd_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call allocator(g_kl_cd, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
            do l = 1, wf%n_o
               do d = 1, wf%n_v
!
!                 Calculate compound indices
!

                  cd = index_two(c, d, wf%n_v)
                  kl = index_two(k, l, wf%n_o)
                  kc = index_two(k, c, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
! 
                  ci = index_two(c, k, wf%n_v)
                  dj = index_two(d, l, wf%n_v)
                  ij = index_two(k, l, wf%n_o)
!
                  cidj = index_packed(ci, dj)
!
!                 Reordering
!
                  g_kl_cd(kl, cd) = g_kc_ld(kc, ld)
                  t_cd_ij(cd, ij) = t2am(cidj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_kc_ld
!
      call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',
                  (wf%n_o)*(wf%n_o), &
                  (wf%n_o)*(wf%n_o), &
                  (wf%n_v)*(wf%n_v), &
                  one,               &
                  g_kl_cd,           &
                  n_oo,              &
                  t_cd_ij,           &
                  (wf%n_v)*(wf%n_v), &
                  one,               &
                  g_kl_ij,           &
                  (wf%n_o)*(wf%n_o))
!
!     Deallocate t_cd_ij and g_kl_cd
!
      call deallocator(t_cd_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call deallocator(g_kl_cd, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
!     Reorder t_ak_bl to t_ab_kl
!
      call allocator(t_ab_kl, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do k = 1, wf%n_o
               do l=1, wf%n_o
!
!                 Calculate compound indices
!
                  ak = index_two(a, k, wf%n_v)
                  bl = index_two(b, l, wf%n_v)
                  ab = index_two(a, b, wf%n_v)
                  kl = index_two(k, l, wf%n_o)
!
                  akbl = index_packed(ak, bl)
!
                  t_ab_kl(ab, kl) = t2am(akbl, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     omega_ab_ij = sum_(kl) t_ab_kl*X_kl_ij
!
      call allocator(omega_ab_ij, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_v), &
                  (wf%n_o)*(wf%n_o), &
                  (wf%n_o)*(wf%n_o), &
                  one,               &
                  t_ab_kl,           &
                  (wf%n_v)*(wf%n_v), &
                  g_kl_ij,           &
                  (wf%n_o)*(wf%n_o), &
                  zero,              &
                  omega_ab_ij,       &
                  (wf%n_v)*(wf%n_v))
!
      call deallocator(t_ab_kl, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call deallocator(g_kl_ij, (wf%n_o)*(wf%n_o), (wf%n_o)*(wf%n_o))
!
!     Reorder omega
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do j = 1, wf%n_o
!
!                 Calculate compound indices
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  if (ai .ge. bj) then
!
                     ab = index_two(a, b, wf%n_v)
                     ij = index_two(i, j, wf%n_o)
!
                     aibj = index_packed(ai, bj)
!
                     omega2(aibj, 1) = omega2(aibj, 1) + omega_ab_ij(ab, ij)
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(omega_ab_ij,(wf%n_v)*(wf%n_v),(wf%n_o)*(wf%n_o))  
!
!     Print the omega vector, having added B2
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after B2 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,wf%n_t2am)
      endif
!
   end subroutine omega_b2_cc_singles_doubles
!
!
   subroutine omega_c2_cc_singles_doubles(wf)
!
!     C2 omega term. Omega C2 = -1/2* sum_(ck)t_bk_cj*(g_ki_ac -1/2 sum_(dl)t_al_di * g_kd_lc)
!                               - sum_(ck) t_bk_ci (g_kj_ac-sum_(dl)t_al_dj*g_kd_lc)
!
!
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!     
      implicit none
!
      class(cc_singles_doubles) :: wf
!
!     Integrals
!
      real(dp),dimension(:,:),allocatable :: L_ia_J 
      real(dp),dimension(:,:),allocatable :: L_ki_J 
      real(dp),dimension(:,:),allocatable :: L_ca_J 
      real(dp),dimension(:,:),allocatable :: g_kd_lc
      real(dp),dimension(:,:),allocatable :: g_dl_ck
      real(dp),dimension(:,:),allocatable :: g_ki_ca
      real(dp),dimension(:,:),allocatable :: g_ai_ck
!
!     Reordered T2 amplitudes
!
      real(dp),dimension(:,:),allocatable :: t_ai_dl
      real(dp),dimension(:,:),allocatable :: t_ck_bj
!
!     Intermediates for matrix multiplication
!
      real(dp),dimension(:,:),allocatable :: X_ai_ck
      real(dp),dimension(:,:),allocatable :: Y_ai_bj
!  
!     Indices
!     
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ca = 0
      integer(i15) :: ai = 0, aj = 0, al = 0, bi = 0, bj = 0, bk = 0, cj = 0, ck = 0, cl = 0, di = 0, dk = 0, dl = 0
      integer(i15) :: kd = 0, lc = 0
      integer(i15) :: ki = 0
!
      integer(i15) :: aldi = 0, aibj = 0, cldk = 0, bkcj = 0
!
!     Batching and memory handling
!
      integer(i15) :: required = 0, available = 0
!
      integer(i15) :: n_batch = 0, max_batch_length = 0
      integer(i15) :: a_batch = 0, a_start = 0, a_end = 0, a_length = 0 
!
!     Debug
!
      logical :: debug = .false.
!
!
!     Allocate L_ia_J
!
      call allocator(L_ia_J,(wf%n_o)*(wf%n_v),(wf%n_J))
      L_ia_J=zero
!
!     Get L_ia_J
!
      call get_cholesky_ia(L_ia_J)
!
!     Create g_kd_lc = sum_J L_kd_J * L_lc_J
!
      call allocator(g_kd_lc,(wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
      g_kd_lc = zero
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &   
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kd_lc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate L_ia_J
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Reorder g_kd_lc as g_dl_ck and t_al_di as t_ai_dl
!
      call allocator(g_dl_ck,n_ov,n_ov)
      call allocator(t_ai_dl,n_ov,n_ov)
      g_dl_ck = zero
      t_ai_dl = zero
!
      do c = 1, wf%n_v
         do d = 1, wf%n_v
            do k = 1, wf%n_o
               do l = 1, wf%n_o
!
!                 Calculate compound indices
!
                  kd = index_two(k, d, wf%n_o)
                  lc = index_two(l, c, wf%n_o)
                  dl = index_two(d, l, wf%n_v)
                  ck = index_two(c, k, wf%n_v)
!
                  cl = index_two(c, l, wf%n_v)
                  dk = index_two(d, k, wf%n_v)
!
                  cldk = index_packed(cl, dk)
!
!                 Reorderings
!
                  g_dl_ck(dl, ck) = g_kd_lc(kd, lc)
!
                  t_ai_dl(ck, dl) = t2am(cldk, 1)
!
               enddo
            enddo
         enddo
      enddo

      call deallocator(g_kd_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     -1/2*sum_(dl) t_ai_dl*g_dl_ck = X_ai_ck
!
      call allocator(X_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -half,             &
                  t_ai_dl,           &
                  (wf%n_o)*(wf%n_v), &
                  g_dl_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ai_ck,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate L_ia_J and g_dl_ck, 
!
      call deallocator(g_dl_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Constructing g_ki_ac ordered as g_ki_ca
!
!     Allocate g_ki_ca
!
      call allocator(g_ki_ca, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
      g_ki_ca = zero
!
!     Allocate L_ki_J
!
      call allocator(L_ki_J, (wf%n_o)*(wf%n_o), wf%n_J) 
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
      required = 2*(wf%n_v)*(wf%n_v)*(wf%n_J)*4 + 2*(wf%n_v)*(wf%n_o)*(wf%n_J)*4
      call n_one_batch(required, available, max_batch_length, n_batch, wf%n_v)
!
      a_start  = 1
      a_end    = 0
      a_length = 0
!
!     Start looping over batches
!
      do a_batch = 1,n_batch
!
!        Get batch limits  and  length of batch
!
         call one_batch_limits(a_start, a_end, a_batch, max_batch_length, wf%n_v)
         a_length = a_end - a_start + 1
!
!        Allocation for L_ac_J as L_ca_J (L_ca_J = L_acJ)
!
         call allocator(L_ca_J,(wf%n_v)*a_length, wf%n_J)
         L_ca_J=zero
!
!        Read Cholesky vectors
!
         call get_cholesky_ab(L_ca_J, a_start, a_end, (wf%n_v)*a_length, .true.)
!
!        g_ki_ca = sum_J L_ki_J*L_ca_J
!
         call dgemm('N','T',                          &
            (wf%n_o)*(wf%n_o),                        &
            (wf%n_v)*a_length,                        &
            wf%n_J,                                   &   
            one,                                      &   
            L_ki_J,                                   &
            (wf%n_o)*(wf%n_o),                        &
            L_ca_J,                                   &
            (wf%n_v)*a_length,                        &
            one,                                      &
            g_ki_ca(1,index_two(1, a_start, wf%n_v)), &
            (wf%n_o)*(wf%n_o))
!
!        Deallocate L_ca_J
!
         call deallocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
!
      enddo ! End of batching
!
!     Deallocate L_ki_J
!
      call deallocator(L_ki_J, (wf%n_o)*(wf%n_o), wf%n_J)
!
!     Reorder g_ki_ca to g_ai_ck
!
      call allocator(g_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      do i = 1, wf%n_o
         do k = 1, wf%n_o
            do a = 1, wf%n_v
               do c = 1, wf%n_v
!
!                 Calculate compound indices
!
                  ki = index_two(k, i, wf%n_o)
                  ca = index_two(c, a, wf%n_v)
                  ai = index_two(a, i, wf%n_v)
                  ck = index_two(c, k, wf%n_v)
!
                  g_ai_ck(ai, ck) = g_ki_ca(ki, ca)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ki_ca, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
!     X_ai_ck = X_ai_ck + g_ai_ck
!
      call daxpy((wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v), one, g_ai_ck, 1, X_ai_ck, 1)
!
!     Deallocate g_ai_kc
!
      call deallocator(g_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     reorder t_bkcj_1 as t_ck_bj
!
      call allocator(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Needed indices
!
                  bk = index_two(b, k, wf%n_v)
                  cj = index_two(c, j, wf%n_v)
!
                  bkcj = index_packed(bk, cj)
!
                  bj = index_two(b, j, wf%n_v)
                  ck = index_two(c, k, wf%n_v)
!
                  t_ck_bj(ck, bj) = t2am(bkcj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate intermediate Y_ai_bj
!
      call allocator(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      Y_ai_bj = zero
!
!     Y_ai_bj = - sum_(ck) X_ai_ck*t_ck_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  X_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  t_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_ai_bj,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate the X intermediate
!
      call deallocator(X_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Deallocate t_ck_bj
!
      call deallocator(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Omega_aibj,1 =P_ai_bj( 1/2*Y_ai_bj + Y_aj_bi)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate compound indices
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  if (ai .ge. bj) then
                     aj = index_two(a, j, wf%n_v)
                     bi = index_two(b, i, wf%n_v)
!
                     aibj=index_packed(ai, bj)
!
                     omega2(aibj, 1)=omega2(aibj, 1) + half*Y_ai_bj(ai, bj) + Y_ai_bj(aj, bi) &
                                                     + half*Y_ai_bj(bj, ai) + Y_ai_bj(bi, aj)
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
      call deallocator(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Print the omega vector, having added C2
!
      if (debug) then 
         write(luprint,*) 
         write(luprint,*) 'Omega(aibj,1) after C2 term has been added:'
         write(luprint,*)
         call vec_print_packed(omega2,wf%n_t2am)
      endif 
!
   end subroutine omega_c2_cc_singles_doubles
!
!  
end submodule omega