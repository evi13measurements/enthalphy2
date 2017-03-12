module mlcc_energy
!
!  MLCC CCSD Ground State Module
!  Written by Eirik F. Kjønstad and Sarai D. Folkestad, 10 Mar 2017
!
   use mlcc_workspace
   use mlcc_data
   use mlcc_omega 
!
contains
   subroutine mlcc_energy_drv
!
!     CCSD Energy Driver
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, 10 Mar 2017
!
!     Controls the iterative loop for the calculation of the coupled 
!     cluster singles and doubles ground state energy 
!
      implicit none 
!
      logical :: debug = .true.
!
      logical :: converged           = .false.
      logical :: converged_energy    = .false.
      logical :: converged_solution  = .false.
!
      integer :: max_iterations = 1
      integer :: iteration = 1
!
      real(dp) :: energy_threshold = 1.0D-8
      real(dp) :: solution_threshold = 1.0D-8
      real(dp) :: energy = zero 
      real(dp) :: prev_energy = zero
      real(dp) :: omega_norm
!
!     Enter the iterative loop 
!
		do while ((.not. converged) .and. (iteration .le. max_iterations))
!
!        Calculate the current coupled cluster energy 
! 
         prev_energy = energy 
         call mlcc_ccsd_energy(energy)
!
!        Reset the omega vector and re-calculate it
!  
         omega1 = zero
         omega2 = zero
!         
         call mlcc_omega_calc
!
!        Test for convergence of the omega vector and the energy 
!
         call mlcc_norm(omega_norm,omega1,omega2)
!
         converged_energy   = abs(energy-prev_energy) .lt. energy_threshold
         converged_solution = omega_norm              .lt. solution_threshold
!
         if (converged_energy .and. converged_solution) converged = .true.
         if (converged_energy) write(luprint,*) 'Energy has converged!'
         if (converged_solution) write(luprint,*) 'Solution has converged!'
!
!        Find the next set of amplitudes, by a DIIS step,
!        if the equations have not yet converged
!
         if (.not. converged) then 
            call mlcc_ccsd_update_amplitudes
            iteration = iteration + 1
         endif
!
!        Print some information necessary for debug purposes
!
         if (debug) then 
!
            write(luprint,*)
            write(luprint,*) 'Iteration nr.:',iteration-1
            write(luprint,*)
            write(luprint,*) 'Current CCSD energy:',energy 
            write(luprint,*) 'Energy difference:',abs(energy-prev_energy)
            write(luprint,*) 'Norm of solution:',omega_norm
            write(luprint,*) 
            write(luprint,*) 'Energy convergence:',converged_energy
            write(luprint,*) 'Solution convergence:',converged_solution
            write(luprint,*) 'Convergence:',converged
            write(luprint,*)
!
         endif
!
		enddo
!
	end subroutine mlcc_energy_drv
!
   subroutine mlcc_ccsd_energy(energy)
!
!     CCSD Energy 
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, 10 Mar 2017
!
!     Calculates the coupled cluster energy 
!
!        E = E_HF + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb,
!
!     where L_iajb = 2 * g_iajb - g_ibja
!
      implicit none 
!
      real(dp) :: energy 
!
      integer :: a,i,b,j,ai,bj,aibj,ia,jb,ib,ja
!
      real(dp), dimension(:,:), pointer :: L_ia_J  => null() ! L_ia_J
      real(dp), dimension(:,:), pointer :: g_ia_jb => null() ! g_iajb
!
!     Allocate the Cholesky vector L_ia_J = L_ia^J and set to zero 
!
      call allocator(L_ia_J,n_ov,n_J)
      L_ia_J = zero
!
!     Read the Cholesky vector L_ia_J from file 
!
      call read_cholesky_ia(L_ia_J)
!
!     Allocate g_ia_jb = g_iajb and set it to zero
!
      call allocator(g_ia_jb,n_ov,n_ov)
      g_ia_jb = zero
!
!     Calculate the integrals g_ia_jb from the Cholesky vector L_ia_J 
!
      call dgemm('N','T',n_ov,n_ov,n_J,&
                  one,L_ia_J,n_ov,L_ia_J,n_ov,&
                  zero,g_ia_jb,n_ov)
!
!     Deallocate the Cholesky vector L_ia_J 
!
      call deallocator(L_ia_J,n_ov,n_J)
!
!     Set the initial value of the energy 
!
      energy = scf_energy
!
!     Add the correlation energy E = E + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!
      do i = 1,n_occ
         do a = 1,n_vir
            do j = 1,n_occ
               do b = 1,n_vir
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a,i,n_vir)
                  bj   = index_two(b,j,n_vir)
                  aibj = index_packed(ai,bj)
!
                  ia   = index_two(i,a,n_occ)
                  jb   = index_two(j,b,n_occ)
!
                  ib   = index_two(i,b,n_occ)
                  ja   = index_two(j,a,n_occ)
!
!                 Add the correlation energy 
!
                  energy = energy + (t2am(aibj,1)+t1am(a,i)*t1am(b,j))*(two*g_ia_jb(ia,jb) - g_ia_jb(ib,ja))
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine mlcc_ccsd_energy
!
   subroutine mlcc_norm(norm,vec1,vec2)
!
!     MLCC Norm
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, 10 Mar 2017
!
!     Calculates the norm of the singles-packed-doubles vector vec 
!
      implicit none 
!
      integer :: a,i,b,j,ai,bj,aibj
!
      real(dp) :: norm,norm1,norm2,ddot
!
      real(dp), dimension(n_vir,n_occ)      :: vec1
      real(dp), dimension(n_ov_ov_packed,1) :: vec2
!
!     Calculate the norm 
!
      norm1 = ddot(n_ov,vec1,1,vec1,1)             ! vec1^2
      norm2 = ddot(n_ov_ov_packed,vec2,1,vec2,1)   ! vec2^2
      norm  = sqrt(norm1 + norm2)                  ! norm = sqrt(vec1^2 + vec2^2)
!
   end subroutine mlcc_norm
!
   subroutine mlcc_ccsd_update_amplitudes
!
!     CCSD Update amplitudes
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, 10 Mar 2017
!
!     Performs a DIIS update of the amplitudes for the next iteration in the solver
! !
!       implicit none
! !
! !
! !
!    ! FROM CCSD_NXTAM
!    ! SINGLES BIT
!             DO 100 ISYMI = 1,NSYM
!             ISYMA = MULD2H(ISYMT,ISYMI)
!             DO 110 I = 1,NRHF(ISYMI)
!                KOFFI = IRHF(ISYMI) + I
!                DO 120 A = 1,NVIR(ISYMA)
! C
!                   KOFFA = IVIR(ISYMA) + A
!                   NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
! C
!                   OMEGA1(NAI) = T1AM(NAI) + OMEGA1(NAI)/
!      *                       (FREQ + FCDIAG(KOFFI) - FCDIAG(KOFFA))
! C
!   120          CONTINUE
!   110       CONTINUE
!   100    CONTINUE
!   ! DOUBLES BIT
!        DO 200 ISYMBJ = 1,NSYM
!          ISYMAI = MULD2H(ISYMBJ,ISYMT)
!          IF (ISYMAI .LE. ISYMBJ) THEN
!          DO 210 ISYMJ = 1,NSYM
!             ISYMB = MULD2H(ISYMJ,ISYMBJ)
!             DO 220 ISYMI = 1,NSYM
!                ISYMA = MULD2H(ISYMI,ISYMAI)
!                DO 230 J = 1,NRHF(ISYMJ)
!                   KOFFJ = IRHF(ISYMJ) + J
!                   DO 240 B = 1,NVIR(ISYMB)
!                      KOFFB = IVIR(ISYMB) + B
!                      NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
!                      DO 250 I = 1,NRHF(ISYMI)
!                         KOFFI = IRHF(ISYMI) + I
!                         DO 260 A = 1,NVIR(ISYMA)
!                            KOFFA = IVIR(ISYMA) + A
!                            NAI = IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+A
! C
!                            IF (ISYMAI.EQ.ISYMBJ .AND. NAI.LE.NBJ) THEN
!                               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
!      *                            + INDEX(NAI,NBJ)
!                               OMEGA2(NAIBJ) = T2AM(NAIBJ)+OMEGA2(NAIBJ)/
!      *                                   (FREQ + 
!      *                                    FCDIAG(KOFFI) + FCDIAG(KOFFJ)
!      *                                  - FCDIAG(KOFFA) - FCDIAG(KOFFB))
!                            ELSE IF (ISYMAI.LT.ISYMBJ) THEN
!                               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
!      *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
!                               OMEGA2(NAIBJ) = T2AM(NAIBJ)+OMEGA2(NAIBJ)/
!      *                                   (FREQ + 
!      *                                    FCDIAG(KOFFI) + FCDIAG(KOFFJ)
!      *                                  - FCDIAG(KOFFA) - FCDIAG(KOFFB))
!                            ENDIF
! C
! C
!   260                   CONTINUE
!   250                CONTINUE
!   240             CONTINUE
!   230          CONTINUE
!   220       CONTINUE
!   210    CONTINUE
!          ENDIF
!   200 CONTINUE
!
! FROM CCSD ENERGY. 
! CALL CCSD_DIIS(WORK(KOMEG1),WORK(KTAMP1),NTAMP,ITER)
!
! FROM CCSD_DIIS
!
 !     SUBROUTINE CCSD_DIIS(S,T,NDIM,NITER) ! S = OMEG1 ... OMEG2 ... T = T1AM ... T2AM ... NDIM = Nt1am + nt2am, ITER = iteration 
!
! This is an unreadable routine!! But I will read it, of course.
!
   end subroutine mlcc_ccsd_update_amplitudes
!
end module mlcc_energy
!
!
!
!
!
