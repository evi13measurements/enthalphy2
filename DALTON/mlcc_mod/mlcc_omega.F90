module mlcc_omega
!
! 	MLCC Omega 
! 	Written by Eirik F. Kjønstad and Sarai Folkestad, February 2017
!
! 	Routines for the calculation of the < mu | exp(-T) H exp(T) | R >, the omega vector
!
   use mlcc_data
!
contains 
   subroutine mlcc_omega_calc
!
!     MLCC Omega calculation
!     Written by Eirik F. Kjønstad and Sarai Folkestad, February 2017
!
!     Controls the calculation of the omega vector
!
      implicit none
!
!     I. The singles contribution to < mu | exp(-T) H exp(T) | R >
!
      call mlcc_omega_a1
      call mlcc_omega_b1
      call mlcc_omega_c1
      call mlcc_omega_d1
!
!     II. The doubles contribution to < mu | exp(-T) H exp(T) | R >
!
      call mlcc_omega_e2
      call mlcc_omega_d2
      call mlcc_omega_c2
      call mlcc_omega_a2
      call mlcc_omega_b2
!
   end subroutine mlcc_omega_calc
!
   subroutine mlcc_omega_a1
!
!     MLCC Omega A1 term ( sum_ckd g_adkc * u_ki^cd )
!     Written by Eirik F. Kjønstad and Sarai Folkestad, February 2017
!
!     NB! Needs to be rewritten with T1 transformed integrals
!     eventually (makes no difference for MP2 guess)
!
      implicit none
!
      write(ml_lupri,*) 'In mlcc_omega_a1'
!
   end subroutine mlcc_omega_a1
!
   subroutine mlcc_omega_b1
!
!     MLCC Omega B1 term ( - sum_ckl u_kl^ac * g_kilc )
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, February 2017
!
!     NB! Needs to be rewritten with T1 transformed integrals
!     eventually (makes no difference for MP2 guess)
!
      use mlcc_workspace
      use mlcc_utilities
!
      implicit none
!
      integer :: n_ki, n_lc ! Number of elements in the two Cholesky vectors (for a particular J)
      integer :: n_ckl, n_i ! Number of rows and columns in g_ckl_i (see below)
!
      integer :: lucho_ij,lucho_ia,idummy,j,i
!
      integer :: c,k,l,ckl,ki,cl
!
      real(dp), dimension(:,:), pointer :: L_ki_J_packed => null() ! (ki) index is packed
      real(dp), dimension(:,:), pointer :: L_lc_J => null()
      real(dp), dimension(:,:), pointer :: g_ki_lc => null()
      real(dp), dimension(:,:), pointer :: g_ckl_i => null() ! Reordered two-electron integrals
!
      write(ml_lupri,*) 'In mlcc_omega_b1'
!
!     I. Calculation and reordering of g_ki,lc = sum_J L_ki^J * L_lc^J 
!
!     Allocate Cholesky vectors L_ki,J and L_lc,J 
!
      call allocator(L_ki_J_packed,n_occ*(n_occ+1)/2,n_J)
      call allocator(L_lc_J,n_occ*n_vir,n_J)
!
      L_ki_J_packed = zero
      L_lc_J = zero
!
!     Read L_ki,J
!
      lucho_ij = -1
      call gpopen(lucho_ij,'CHOLESKY_IJ','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
      rewind(lucho_ij)
!
      do j = 1,n_J
         read(lucho_ij) (L_ki_J_packed(i,j), i=1,n_oo_packed)
      enddo
!
      call gpclose(lucho_ij,'KEEP')
!
!     Read L_lc,J
!
      lucho_ia = -1
      call gpopen(lucho_ia,'CHOLESKY_IA','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
      rewind(lucho_ia)
!
      do j=1,n_J
         read(lucho_ia) (L_lc_J(i,j), i=1,n_ov)
      enddo
!
      call gpclose(lucho_ia,'KEEP')
!
!     Allocate integrals g_ki_lc
!
      call allocator(g_ki_lc,n_oo_packed,n_ov) ! (ki) index is packed 
!
      g_ki_lc = zero
!
!     Calculate g_ki_lc = sum_J L_ki,J * L_lc,J^T 
! 
      call dgemm('N','T',n_oo_packed,n_ov,n_J,one,L_ki_J_packed,n_oo_packed,L_lc_J,n_ov,zero,g_ki_lc,n_oo_packed) 
!
!     Deallocate the Cholesky vectors 
!
      call deallocator(L_ki_J_packed,n_ki,n_J)
      call deallocator(L_lc_J,n_lc,n_J)
!
!     Allocate reordered integrals g_ckl,i 
!
      call allocator(g_ckl_i,n_oov,n_occ)
!
!     Save reordered integrals g_ckl,i
!
      do c = 1,n_vir
         do k = 1,n_occ
            do l = 1,n_occ
               do i = 1,n_occ
!
!                 Calculate necessary indices
!
                  ckl = index_three(c,k,l,n_vir,n_occ) 
                  ki  = index_packed(k,i)                    ! Packed 
                  cl  = index_two(c,l,n_vir)             ! Not packed
!
!                 Set value of g_ckl_i
!
                  g_ckl_i(ckl,i) = g_ki_lc(ki,cl)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate unordered integrals g_ki_lc
!
      call deallocator(g_ki_lc,n_oo_packed,n_ov)
!
   end subroutine mlcc_omega_b1
!
   subroutine mlcc_omega_c1
!
!  Purpose: Calculate C1 term of Omega
!  (Omega_ai^C1=sum_ck F_kc u^ac_ik)
!
!  Written by Sarai D. Folkestad and Eirik F. Kjønstad, 28. Feb 2017
!
   use mlcc_data
   use mlcc_utilities
   use mlcc_workspace
   implicit none
!
   real(dp),dimension(:,:),pointer  :: F_kc
   real(dp),dimension(:,:),pointer  :: u_ck_ai
   integer                          :: i,k,c,a
   integer                          :: nck,nai,nak,nci,nckai,nciak
!
!
!  Allocation
!
   call allocator(u_ck_ai,n_ov,n_ov)   ! 2*t_ck_ai-t_ci_ak
!
! Set up u_ck_ai
!
   do i=1,n_occ
      do a=1,n_vir
         do k=1,n_occ
            do c=1,n_vir
               nck=index_two(c,k,n_vir)
               nai=index_two(a,i,n_vir)
               if (nck .le. nai) then
                  nci=index_two(c,i,n_vir)
                  nak=index_two(a,k,n_vir)
                  nckai=index_packed(nck,nai)
                  nciak=index_packed(nci,nak)
                  u_ck_ai(nck,nai)=two*t2am(nckai,1)-t2am(nciak,1)
               endif
            enddo
         enddo
      enddo
   enddo
!
! Allocation 
!
   call allocator(F_kc,n_ov,1) ! T1-transformed fock matrix
!
!  IO - Read fock matrix
!
!
!  T1 transformation
!
!
!  Deallocation
!
   call deallocator(F_kc,n_ov,1)
   call deallocator(u_ck_ai,n_ov,n_ov)

   end subroutine mlcc_omega_c1
!
   subroutine mlcc_omega_d1
   end subroutine mlcc_omega_d1
!
   subroutine mlcc_omega_e2
   end subroutine mlcc_omega_e2
!
   subroutine mlcc_omega_d2
   end subroutine mlcc_omega_d2
!
   subroutine mlcc_omega_c2
   end subroutine mlcc_omega_c2
!
   subroutine mlcc_omega_a2
   end subroutine mlcc_omega_a2
!
   subroutine mlcc_omega_b2
   end subroutine mlcc_omega_b2
!
end module mlcc_omega