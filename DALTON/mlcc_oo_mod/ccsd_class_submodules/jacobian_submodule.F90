submodule (ccsd_class) jacobian
!
!
!                     -::- Jacobian submodule (CCSD) -::-                             
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017         
!                                                                           
!
contains
!
!
   module subroutine jacobian_transformation_ccsd(wf,c1am,c2am)
!
!     Jacobian Transformation (CCSD)
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!
!     Directs the transformation of the incoming vector c by the 
!     coupled cluster Jacobian matrix, 
!
!        A_mu,nu = < mu | [e^(-T) H e^(T),tau_nu] | R >.
!
!     On exit, A*c is placed in the incoming c vector.
!
!     Reads doubles amplitudes from file, but assumes the singles are held in 
!     memory from a ground state calculation. The routine assumes that the omega 
!     vector and the doubles amplitudes are deallocated, in terms of memory 
!     presumed available (~ 5 * v^2 * o^2).
!
      class(ccsd) :: wf 
!
      real(dp), dimension(:,:) :: c1am 
      real(dp), dimension(:,:) :: c2am 
!
      real(dp), dimension(:,:), allocatable :: tr1am ! Transformed vector, singles 
      real(dp), dimension(:,:), allocatable :: tr2am ! Transformed vector, doubles 
!
!     Allocate the transformed vector and set it to zero
!
      call allocator(tr1am, wf%n_v, wf%n_o)
      call allocator(tr2am, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dzero(tr1am, (wf%n_v)*(wf%n_o))
      call dzero(tr2am, (wf%n_v)**2*(wf%n_o)**2)
!
!     Calculate the singles contributions 
!
      call wf%jacobian_a1(tr1am,c1am)
!
      write(unit_output,*) 'tr1am(a,i)'
      call vec_print(tr1am, wf%n_v, wf%n_o)
!
      call wf%jacobian_b1(tr1am,c1am)
!
   end subroutine jacobian_transformation_ccsd
!
!
   module subroutine jacobian_a1_ccsd(wf,tr1am,c1am)
!
!     Jacobian A1 
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Calculates the A1 term, 
!
!        sum_c F_ac c_ci - sum_k c_ak F_ki + sum_ck (2 g_aikc - g_acki) c_ck
!
!     and adds it to the transformed singles vector element tr1am(a,i).
!
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: tr1am 
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1am 
!
      real(dp), dimension(:,:), allocatable :: L_ai_J ! L_ai^J 
      real(dp), dimension(:,:), allocatable :: L_kc_J ! L_kc^J
!
      real(dp), dimension(:,:), allocatable :: g_ai_kc ! g_aikc 
      real(dp), dimension(:,:), allocatable :: g_ai_ck ! g_aikc reordered
                                                
!
      real(dp), dimension(:,:), allocatable :: L_ki_J ! L_ki^J
      real(dp), dimension(:,:), allocatable :: L_ca_J ! L_ac^J
!
      real(dp), dimension(:,:), allocatable :: g_ca_ki ! g_acki 
      real(dp), dimension(:,:), allocatable :: g_ia_ck ! g_acki 
!
      real(dp), dimension(:,:), allocatable :: tr1am_tmp ! Temporary for a batching
!
!     Indices
!
      integer(i15) :: i = 0, a = 0, k = 0, c = 0
      integer(i15) :: ai = 0, kc = 0, ck = 0, ca = 0
      integer(i15) :: ki = 0, ia = 0
!
!     Batch variables
!
      integer(i15) :: n_batch = 0, max_batch_length = 0
      integer(i15) :: batch_dimension = 0, available = 0
      integer(i15) :: required = 0, a_batch = 0, a_begin = 0
      integer(i15) :: a_end = 0, batch_length = 0, ac_dim = 0
!
      logical :: reorder
!
!     Construct the Fock matrix to make sure it is up to date
!
      call wf%construct_fock
!
!     :: Calculate sum_c F_ac c_ci ::
!
      call dgemm('N','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ab, &
                  wf%n_v,     &
                  c1am,       &
                  wf%n_v,     &
                  one,        &
                  tr1am,      &
                  wf%n_v)
!
!     :: Calculate - sum_k c_ak F_ki ::
!
      call dgemm('N','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  one,        &
                  c1am,       &
                  wf%n_v,     &
                  wf%fock_ij, &
                  wf%n_o,     &
                  -one,       &
                  tr1am,      &
                  wf%n_v)
!
!     :: Calculate 2 * sum_ck g_aikc c_ck :: 
!
!     Construct g_ai_kc = g_aikc
!
      call allocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ai(L_ai_J)
      call wf%get_cholesky_ia(L_kc_J)
!
      call allocator(g_ai_kc, (wf%n_v)*(wf%n_o), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ai_J,            &
                  (wf%n_v)*(wf%n_o), &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ai_kc,           &
                  (wf%n_v)*(wf%n_o))
!
!     Deallocate Cholesky vectors 
!
      call deallocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate and calculate the reordered integral g_ai_ck = g_aikc 
!
      call allocator(g_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!
            do k = 1, wf%n_o
               do c = 1, wf%n_v
!
                  kc = index_two(k, c, wf%n_o)
                  ck = index_two(c, k, wf%n_v)
!
                  g_ai_ck(ai, ck) = g_ai_kc(ai, kc)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate unordered integral
!
      call deallocator(g_ai_kc, (wf%n_v)*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!     Calculate and add 2 * sum_ck g_aikc c_ck = 2 * sum_ck g_ai_ck c_ck
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  g_ai_ck,           &
                  (wf%n_v)*(wf%n_o), &
                  c1am,              &
                  (wf%n_v)*(wf%n_o), &
                  two,               &
                  tr1am,             &
                  (wf%n_v)*(wf%n_o))
!
!     Deallocate the reordered integral
!
      call deallocator(g_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Calculate - sum_ck g_acki c_ck ::
!
!     Allocate and get Cholesky vector L_ki_J = L_ki^J
!
      call allocator(L_ki_J, (wf%n_o)**2, wf%n_J)
!
      call wf%get_cholesky_ij(L_ki_J)
!
!     Set up batching variables
!
!     Note: in forming g_acki, we need to hold L_ac_J and two g_acki's 
!     in memory, requiring a total memory of v^2 * n_J + 2 v^2 o^2
!
      required = (wf%n_v)**2*(wf%n_J) + 2*(wf%n_v)**2*(wf%n_o)**2
      available = get_available()
!
      batch_dimension = wf%n_v ! a in L_ac^J and g_acki 
!
      call num_batch(required, available, max_batch_length, &
                     n_batch, batch_dimension)
!
!     Allocate vector to temporarily hold result
!     tr1am_tmp(i,a) is added to tr1am(a,i) after the
!     calculation over batches
!
      call allocator(tr1am_tmp, wf%n_o, wf%n_v)
      call dzero(tr1am_tmp, (wf%n_o)*(wf%n_v))
!
!     Begin batching over a 
!
      do a_batch = 1, n_batch 
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_begin, a_end, a_batch, &
                           max_batch_length, batch_dimension)
!
         batch_length = a_end - a_begin + 1
         ac_dim = batch_length*(wf%n_v)
!
!        Form the Cholesky vector L_ca_J = L_ac^J
!
         call allocator(L_ca_J, batch_length*(wf%n_v), wf%n_J)
!
         reorder = .true.
         call wf%get_cholesky_ab(L_ca_J, a_begin, a_end, &
                                 ac_dim, reorder)
!
!        Form the integral g_ca_ki = g_acki
!
         call allocator(g_ca_ki, ac_dim, (wf%n_o)**2)
!
         call dgemm('N','T',      &
                     ac_dim,      &
                     (wf%n_o)**2, &
                     wf%n_J,      &
                     one,         &
                     L_ca_J,      &
                     ac_dim,      &
                     L_ki_J,      &
                     (wf%n_o)**2, &
                     zero,        &
                     g_ca_ki,     &
                     ac_dim)
!
!        Deallocate Cholesky vector L_ca_J 
!
         call deallocator(L_ca_J, ac_dim, wf%n_J)
!
!        Form the reordered integral g_ia_ck(ia,ck) = g_ca_ki(ca,ki) = g_acki
!
         call allocator(g_ia_ck, batch_length*(wf%n_o), (wf%n_v)*(wf%n_o))
!
         do a = 1, batch_length
            do i = 1, wf%n_o
!
               ia = index_two(i, a, wf%n_o)
!
               do k = 1, wf%n_o
!
                  ki = index_two(k, i, wf%n_o)
!
                  do c = 1, wf%n_v
!
                     ck = index_two(c, k, wf%n_v)
                     ca = index_two(c, a, wf%n_v)
!
                     g_ia_ck(ia,ck) = g_ca_ki(ca,ki)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate the g_ca_ki ordered integral
!
         call deallocator(g_ca_ki, (wf%n_v)*batch_length, (wf%n_o)**2)
!
!        Calculate - sum_ck g_acki c_ck = - sum_ck g_ia_ck c_ck 
!
         call dgemm('N','N',                &
                     batch_length*(wf%n_o), &
                     1,                     &
                     (wf%n_v)*(wf%n_o),     &
                     one,                   &
                     g_ia_ck,               &
                     (wf%n_o)*batch_length, &
                     c1am,                  &
                     (wf%n_v)*(wf%n_o),     &
                     -one,                  &
                     tr1am_tmp(1,a_begin),  &
                     batch_length*(wf%n_o))
!
!        Deallocations
!
         call deallocator(g_ia_ck, batch_length*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      enddo ! end of a batches
!
!     Add the temporary transposed transformed vector
!     into the actual vector 
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            tr1am(a,i) = tr1am(a,i) + tr1am_tmp(i,a)
!
         enddo
      enddo
!
!     Deallocations
!
      call deallocator(tr1am_tmp, wf%n_o, wf%n_v)
      call deallocator(L_ki_J, (wf%n_o)**2, wf%n_J)
!
   end subroutine jacobian_a1_ccsd
!
!
   module subroutine jacobian_b1_ccsd(wf,tr1am,c1am)
!
!     Jacobian B1
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Calculates the B1 term,
!
!       sum_dl (sum_ck L_ldkc c_ck) u_li^da,
!
!     where L_ldkc = 2 * g_ldkc - g_lckd and 
!     u_li^ad = 2 * t_li^ad - t_il^ad.
!
      class(ccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: tr1am 
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1am
!
      real(dp), dimension(:,:), allocatable :: L_ld_J ! L_ld^J
!
      real(dp), dimension(:,:), allocatable :: g_ld_kc ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: L_dl_ck ! L_ldkc reordered
!
      real(dp), dimension(:,:), allocatable :: X ! An intermediate; see below
!
!     Indices
!
      integer(i15) :: l = 0, d = 0, k = 0, c = 0
      integer(i15) :: dl = 0, ld = 0, kd = 0, kc = 0, ck = 0, lc = 0
!
!     :: Calculate the intermediate X_dl = sum_ck L_ldkc c_ck ::
!
!     Form g_ld_kc = g_ldkc = sum_J L_ld^J L_kc^J
!
      call allocator(L_ld_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_ld_J)
!
      call allocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  L_ld_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ld_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ld_kc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate Cholesky vector L_ld_J
!
      call deallocator(L_ld_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Form L_dl_ck = L_ldkc = 2 * g_ld_kc - g_lc_kd
!
      call allocator(L_dl_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      do l = 1, wf%n_o
         do d = 1, wf%n_v
!
            dl = index_two(d, l, wf%n_v)
            ld = index_two(l, d, wf%n_o)
!
            do k = 1, wf%n_o
!
               kd = index_two(k, d, wf%n_o)
!
               do c = 1, wf%n_v
!
                  kc = index_two(k, c, wf%n_o)
                  ck = index_two(c, k, wf%n_v)
                  lc = index_two(l, c, wf%n_o)
!
                  L_dl_ck(dl,ck) = two*g_ld_kc(ld,kc) - g_ld_kc(lc,kd)
!
               enddo
!
            enddo
!
         enddo
      enddo
!
!     Deallocate the unordered integral g_ld_kc
!
      call deallocator(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_dl = sum_ck L_dl_ck c_ck 
!
      call allocator(X, (wf%n_v)*(wf%n_o), 1)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  L_dl_ck,           &
                  (wf%n_v)*(wf%n_o), &
                  c1am,              &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X,                 &
                  (wf%n_v)*(wf%n_o))
!
!     Deallocate integral L_dl_ck 
!
      call deallocator(L_dl_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Calculate sum_dl u_il^ad X_dl = sum_dl u_ai_dl X_dl ::
!
! .... Get the doubles amplitudes from file, & form u.
!
   end subroutine jacobian_b1_ccsd
!
!
end submodule jacobian 