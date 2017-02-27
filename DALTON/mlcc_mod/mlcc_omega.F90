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
!     Written by Eirik F. Kjønstad, February 2017
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
!     Written by Eirik F. Kjønstad, February 2017
!
!     NB! Needs to be rewritten with T1 transformed integrals
!     eventually (makes no difference for MP2 guess)
!
      use mlcc_workspace
!
      implicit none
!
      integer lucho_ij,lucho_ia,idummy,J,i 
!
      real(dp), dimension(:,:), pointer :: L_ki_J
      real(dp), dimension(:,:), pointer :: L_lc_J
      real(dp), dimension(:,:), pointer :: g_ki_lc
!
!     I. Calculation of g_ki,lc = sum_J L_ki^J * L_lc^J 
!
!     Allocate Cholesky vectors L_ki,J and L_lc,J 
!
      call allocator(L_ki_J,n_occ*n_occ,n_J)
      call allocator(L_lc_J,n_occ*n_vir,n_J)
!
!     Allocate integrals g_ki_lc
!
      call allocator(g_ki_lc,n_occ*n_occ,n_occ*n_vir)
!
!     Read L_ki,J
!
      lucho_ij = -1
      call gpopen(lucho_ij,'CHOLESKY_IJ','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
      rewind(lucho_ij)
!
      do J=1,n_J
         read(lucho_ij) (L_ki_J(i,j),i=1,n_occ*n_occ)
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
      do J=1,n_J
         read(lucho_ia) (L_lc_J(i,j),i=1,n_occ*n_vir)
      enddo
!
      call gpclose(lucho_ia,'KEEP')
!
!     Calculate ... TODO
! 

!
   end subroutine mlcc_omega_b1
!
   subroutine mlcc_omega_c1
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