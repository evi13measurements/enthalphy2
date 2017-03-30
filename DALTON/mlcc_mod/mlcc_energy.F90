module mlcc_energy
!
!  MLCC CCSD Ground State Module
!  Written by Eirik F. Kjønstad and Sarai D. Folkestad, 10 Mar 2017
!
   use mlcc_workspace
   use mlcc_data
   use mlcc_omega 
   use mlcc_fock
   use mlcc_utilities
!
!  Some DIIS specific variables 
!
   integer :: maxdiis = 8
   integer :: ludiis_dt = -1, ludiis_t_dt = -1
   real(dp), dimension(:,:), pointer :: G           => null() ! The DIIS matrix, G * w = H, G = G(maxdiis,maxdiis)
   real(dp), dimension(:,:), pointer :: copy_of_G   => null() ! Copy for purposes...
   real(dp), dimension(:,:), pointer :: H           => null() ! The DIIS vector, G * w = H, H = H(maxdiis,1)
   integer,  dimension(:,:), pointer :: lu_integers => null() ! An integer array from LU factorization of G by dgetrf routine
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
      integer :: idum=0
!
      logical :: converged           = .false.
      logical :: converged_energy    = .false.
      logical :: converged_solution  = .false.
!
      integer :: max_iterations = 50
      integer :: iteration = 1
!
      integer :: memory_lef
!
      real(dp) :: energy_threshold = 1.0D-12
      real(dp) :: solution_threshold = 1.0D-12
      real(dp) :: energy = zero 
      real(dp) :: prev_energy = zero
      real(dp) :: omega_norm
!
!      Allocate G matrix
!
      call allocator(G,maxdiis+1,maxdiis+1)
!
!     Enter the iterative loop 
!
      do while ((.not. converged) .and. (iteration .le. max_iterations))
!
!        Calculate the current coupled cluster energy 
! 
         prev_energy = energy 
!
         call mlcc_ccsd_energy(energy) !  Fixed memory leak (Eirik,15 Mar 2017)
         write(luprint,*) 'THE ENERGY:::',energy
!
!        Reset the omega vector and re-calculate it
!  
         omega1 = zero
         omega2 = zero
!         
         write(luprint,*)'Entering mlcc_omega_calc'
         call mlcc_omega_calc
         write(luprint,*)'Exiting mlcc_omega_calc'
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
         if (iteration .eq. 1) then 
!
!        Open the files that stores dt and t + dt 
!  
            ludiis_dt   = -1  ! dt
            call gpopen(ludiis_dt,'DIIS_DT','UNKNOWN','SEQUENTIAL','FORMATTED',idum,.false.)
            rewind(ludiis_dt)
!  
            ludiis_t_dt = -1  ! t + dt
            call gpopen(ludiis_t_dt,'DIIS_T_DT','UNKNOWN','SEQUENTIAL','FORMATTED',idum,.false.)
            rewind(ludiis_t_dt)
!
         endif 
!
!
!        Find the next set of amplitudes, by a DIIS step,
!        if the equations have not yet converged
!
         if (.not. converged) then 
            call mlcc_ccsd_update_amplitudes(iteration)
            iteration = iteration + 1
         endif
!
!        Recalculate the (T1-transformed) Fock matrix
!
         call mlcc_get_fock
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
!     Deallocate DIIS matrix 
!
      call deallocator(G,maxdiis+1,maxdiis+1)
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
      integer :: a=0,i=0,b=0,j=0,ai=0,bj=0,aibj=0,ia=0,jb=0,ib=0,ja=0
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
      call get_cholesky_ia(L_ia_J)
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
!     Deallocate g_ia_jb (Eirik: debugging, 15 Mar 2017)
!
      call deallocator(g_ia_jb,n_ov,n_ov)
!
!     Print the t1 amplitudes 
!
!       write(luprint,*) 't1am(a,i):'
!       call vec_print(t1am,n_vir,n_occ)
! !
!       write(luprint,*) 't2am(aibj,1):'
!       call vec_print_packed(t2am,n_ov_ov_packed)
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
      real(dp) :: norm,norm_1,norm_2,ddot
!
      real(dp), dimension(n_vir,n_occ)      :: vec1
      real(dp), dimension(n_ov_ov_packed,1) :: vec2
!
!     Calculate the norm 
!
      norm_1 = ddot(n_ov,vec1,1,vec1,1)            ! vec1^2
      norm_2 = ddot(n_ov_ov_packed,vec2,1,vec2,1)  ! vec2^2
      norm  = sqrt(norm_1 + norm_2)                ! norm = sqrt(vec1^2 + vec2^2)
!
   end subroutine mlcc_norm
!
   subroutine mlcc_ccsd_update_amplitudes(iteration)
!
!     CCSD Update amplitudes
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, 10 Mar 2017
!
!     Performs a DIIS update of the amplitudes for the next iteration in the solver
!
!     The next amplitudes are 
!
!        t_n+1 = sum_k w_k (t_k + dt_k), 
! 
!     where the weights w_k in front of 
!
!        dt_k = - omega_mu/eps_mu,  eps_ai   = eps_a - eps_i, 
!                                   eps_aibj = eps_a + eps_b - eps_i - eps_j 
!
!     are determined so as to minimize 
!
!        f(w_k) = sum_k w_k dt_k, with the constraint that g(w_k) = sum_k w_k - 1 = 0. (***)
!
!     The routine proceeds as follows:
!
!        1. Determine the quasi-Newton amplitude correction (dt_k)
!        2. Write to file the quasi-Newton amplitude correction (dt_k) & the shifted amplitudes (t_k + dt_k)
!        3. Set up the least-squares eigenvalue problem A w = B associated with (***) by reading the amplitude corrections dt_k from file 
!        4. Use the solution of 3., w, to update the amplitudes --- i.e., t_k+1 =  sum_k w_k (t_k + d_k) ---
!           by reading the shifted amplitudes from file 
!        
      implicit none
!
      integer :: idum=0,iteration,dim_G
!
      integer :: lu_error=0
!
      integer :: a=0,i=0,b=0,j=0,ai=0,bj=0,aibj=0,p=0,current_index=0
!
      real(dp) :: ddot,norm_of_solution
!
      real(dp), dimension(:,:), pointer :: dt1_k => null() ! kth perturbation-corrected quasi-Newton correction, singles 
      real(dp), dimension(:,:), pointer :: dt2_k => null() ! kth perturbation-corrected quasi-Newton correction, doubles
!
      real(dp), dimension(:,:), pointer :: dt1_j => null() ! jth perturbation-correction quasi-Newton correction, singles (j = 1,2,3,... by reading)
      real(dp), dimension(:,:), pointer :: dt2_j => null() ! jth perturbation-correction quasi-Newton correction, doubles (j = 1,2,3,... by reading)
!
      real(dp), dimension(:,:), pointer :: tdt1_j => null() ! jth t + dt vector on file, singles (j = 1,2,3,... by reading)
      real(dp), dimension(:,:), pointer :: tdt2_j => null() ! jth t + dt vector on file, doubles (j = 1,2,3,... by reading)
!
!     Allocate the current (k) & varying (j = 1,2,3,...) quasi-Newton amplitudes correction dt_k and dt_j and set them to zero
!
      call allocator(dt1_k,n_vir,n_occ)
      call allocator(dt2_k,n_ov_ov_packed,1) ! Time could be saved by overwriting the omega vector instead,
                                             ! but this makes things difficult to read & understand (see old DIIS routine)
      call allocator(dt1_j,n_vir,n_occ)
      call allocator(dt2_j,n_ov_ov_packed,1) ! Time could be saved by overwriting the omega vector instead,
                                             ! but this makes things difficult to read & understand (see old DIIS routine)
!
      dt1_j = zero
      dt2_j = zero
!
      dt1_k = zero
      dt2_k = zero 
!
!     Calculate the current d_k vector
!
      do a = 1,n_vir
         do i = 1,n_occ
            do b = 1,n_vir
               do j = 1,n_occ
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a,i,n_vir)
                  bj   = index_two(b,j,n_vir)
                  aibj = index_packed(ai,bj) 
!
!                 Set the value of d_k = - omega_mu/eps_mu
!
                  dt1_k(a,i)  = -omega1(a,i)/(fock_diagonal(n_occ+a,1)-&
                                              fock_diagonal(i,1))
!
                  dt2_k(aibj,1) = -omega2(aibj,1)/(fock_diagonal(n_occ+a,1)+&
                                                   fock_diagonal(n_occ+b,1)-&
                                                   fock_diagonal(i,1)-&
                                                   fock_diagonal(j,1))
!
               enddo
            enddo
         enddo
      enddo
!
!     Calculate the current index for overwriting the G matrix (G(:,current_index) and G(current_index,:) is overwritten)
!
      current_index = modulo(iteration,maxdiis)
      write(luprint,*) 'Current_index,iteration,maxdiis',current_index,iteration,maxdiis
      if (current_index .eq. 0) current_index = maxdiis
!
      if (current_index .eq. 1) then
         write(luprint,*) 'Rewinding both files'
         rewind(ludiis_dt)
         rewind(ludiis_t_dt)
      endif
!
!     Write dt_k to file (singles, then doubles)
!
      write(luprint,*) 'Writing current dt_k to file'
      call flshfo(luprint)
      write(ludiis_dt,*) ((dt1_k(a,i),a=1,n_vir),i=1,n_occ),(dt2_k(p,1),p=1,n_ov_ov_packed)
!
!     Add the quasi-Newton amplitude correction to the amplitudes (t_k <- t_k + dt_k)
!
      call daxpy(n_ov,one,dt1_k,1,t1am,1)
      call daxpy(n_ov_ov_packed,one,dt2_k,1,t2am,1)
!
!     Write t_k + dt_k to file 
!
      write(luprint,*) 'Writing current t_k + dt_k to file'
      call flshfo(luprint)
      write(ludiis_t_dt,*) ((t1am(a,i),a=1,n_vir),i=1,n_occ),(t2am(p,1),p=1,n_ov_ov_packed)
!
!     If the first iteration, then allocate the matrices 
!
      if (iteration .eq. 1) then 
!
!        Allocate the G matrix 
!
!         call allocator(G,maxdiis+1,maxdiis+1) ! Sarai: I think we should always allocate and deallocate in same routine.
!                                               !        Moved this to driver, outsidel loop.
!
!        Set the G matrix and the LU integers array 
!
         G = zero ! Is altered later
!
      endif
!
!     Allocate temporary matrices 
!
      call allocator(copy_of_G,maxdiis+1,maxdiis+1)
      call allocator(H,maxdiis+1,1)
      call allocator_int(lu_integers,maxdiis+1,1)
!
      lu_integers = 0 ! Is altered later 
      copy_of_G = zero
      H = zero ! Fixed throughout the calculation, but is overwritten by dgetrs & must be reset in every iteration 
!
      if (current_index .eq. 1) G = zero
!
!     Calculate the effective dimensionality of G & set its values    
!
      dim_G = current_index
      rewind(ludiis_dt)
!
      do j = 1,dim_G
!
!        Read the jth entry of the file containing the dt's
!
         write(luprint,*) 'Reading dt_j from file, j = ', j 
         call flshfo(luprint)
         read(ludiis_dt,*) ((dt1_j(a,i),a=1,n_vir),i=1,n_occ),(dt2_j(p,1),p=1,n_ov_ov_packed) ! Reads the jth entry of the file 
!
!        Calculate G(current_index,j) = dt_current_index * dt_j + 1
!
         G(current_index,j) = ddot(n_ov,dt1_k,1,dt1_j,1)
         G(current_index,j) = G(current_index,j) + ddot(n_ov_ov_packed,dt2_k,1,dt2_j,1)
         G(j,current_index) = G(current_index,j)
         G(current_index+1,j) = -one
         G(j,current_index+1) = -one 
!
      enddo
      H(dim_G+1,1) = -one
!
      write(luprint,*) 'The G matrix'
      write(luprint,*) G 
      write(luprint,*) 'The H vector'
      write(luprint,*) H
!
!     Solve the eigenvalue problem G * w = H 
!
      copy_of_G = G
      lu_error = -1
      lu_integers = 0
      call dgetrf(maxdiis+1,maxdiis+1,copy_of_G,maxdiis+1,lu_integers,lu_error)
!
      if (lu_error .eq. 0) write(luprint,*) 'Successful LU factorization'
!
      lu_error = -1
      call dgetrs('N',maxdiis+1,1,copy_of_G,maxdiis+1,lu_integers,H,maxdiis+1,lu_error) ! Solution is placed in H
!
      if (lu_error .eq. 0) write(luprint,*) 'Successful solution of G * omega = H'
!
      write(luprint,*) 'The integers',lu_integers
  !    if (current_index .eq. 1) H = 1 ! This is the only solution for this case (hack)
!
      write(luprint,*) 'The H vector (ie, the solution)'
      write(luprint,*) H
!
!     Deallocate the temporary dt vector
!
      call deallocator(dt1_j,n_vir,n_occ)
      call deallocator(dt2_j,n_ov_ov_packed,1)
!
!     Allocate the temporary t + dt vector
!
      call allocator(tdt1_j,n_vir,n_occ)
      call allocator(tdt2_j,n_ov_ov_packed,1)
!
      tdt1_j = zero
      tdt2_j = zero 
!
!     Update the amplitudes 
!
      call dzero(t1am,n_ov)
      call dzero(t2am,n_ov_ov_packed)
      rewind(ludiis_t_dt)
!
      do j = 1,dim_G
!
!        Read the jth t + dt contribution on file 
!
         write(luprint,*) 'Reading t_j + dt_j from file, j = ', j 
         call flshfo(luprint)
         read(ludiis_t_dt,*) ((tdt1_j(a,i),a=1,n_vir),i=1,n_occ),(tdt2_j(p,1),p=1,n_ov_ov_packed) ! Reads the jth entry of the file
!
!        Add the contributions w_j * (t_j + dt_j) to the amplitudes 
!
         write(luprint,*) 'Adding prefactor',H(j,1)
!
         call daxpy(n_ov,H(j,1),tdt1_j,1,t1am,1)
         call daxpy(n_ov_ov_packed,H(j,1),tdt2_j,1,t2am,1)
!
      enddo
!
!     Deallocate vectors 
!
      call deallocator(tdt1_j,n_vir,n_occ)
      call deallocator(tdt2_j,n_ov_ov_packed,1)
      call deallocator(dt1_k,n_vir,n_occ)
      call deallocator(dt2_k,n_ov_ov_packed,1)
      call deallocator(copy_of_G,maxdiis+1,maxdiis+1)
      call deallocator(H,maxdiis+1,1)
      call deallocator_int(lu_integers,maxdiis+1,1)
!
   end subroutine mlcc_ccsd_update_amplitudes
!
end module mlcc_energy
!
!
!
!
!

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
