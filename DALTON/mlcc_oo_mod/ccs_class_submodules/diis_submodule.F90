submodule (ccs_class) diis
!
!
!                          -::- DIIS submodule (CCS) -::-
!             Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!
!   Consists of the following subroutines of the CCS module:
!
!   ground_state_solver:        controls the iterative loop, calling in turn
!                               the calculation of the energy, the amplitude equations 
!                               (and its norm), and the new_amplitudes routine
!   new_amplitudes:             calculates the quasi-Newton estimate and passes the 
!                               information needed by the DIIS routine.
!   diis_ccs:                   This routine saves the quasi-Newton estimate Δ t and
!                               t + Δ t to file. It uses the previous estimates to
!                               select the amplitudes t for the next iteration.
!    
!
!   calc_ampeqs:                Updates the amplitude equations for the current amplitudes.
!   calc_ampeqs_norm:           Calculates the norm of the amplitude equations.
!   calc_quasi_Newton_singles:  Calculates the singles part of the quasi-Newton estimate.
!
!  Can be inherited by models of the same level (e.g. CC2) without modification.
!
!  When inherited by higher level models (e.g. CCSD), the new_amplitudes and calc_ampeqs_norm
!  routines should be overridden to account for the doubles quasi-Newton estimate, amplitudes, 
!  and projection vector.
!
!
   implicit none
!
!  Some variables available to all routines of the module
!
   integer(i15) :: iteration = -1 
!
   integer(i15) :: unit_dt          = -1 ! Unit identifier for Δ t_i file 
   integer(i15) :: unit_t_dt        = -1 ! Unit identifier for t_i + Δ t_i file
   integer(i15) :: unit_diis_matrix = -1 ! Unit identifier for DIIS matrix file
!
   integer(i15), parameter :: diis_dim = 9 ! The maximum dimension of the DIIS matrix plus 1 
!
!
contains
!
!
   module subroutine ground_state_solver_ccs(wf)
!
!     Ground State Solver 
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Directs the solution of the ground state amplitude equations
!     using a DIIS algorithm. The problem the routine solves is 
!
!        X_mu(t) = 0, where t = { t_mu }_mu 
!
!     For standard coupled cluster theories, the vector X is the
!     projection vector (omega).
!
      implicit none
!
      class(ccs) :: wf 
!
      real(dp) :: prev_energy
      real(dp) :: ampeqs_norm
!
      real(dp) :: last_time, first_time
!
      logical :: converged_energy = .false.
      logical :: converged_ampeqs = .false.
!
      logical :: converged = .false. ! True iff both the energy and the equations have converged 
!
!     Let the user know the ground state solver is running
!
      write(unit_output,'(/T3,A)')  ':: Ground state solver (DIIS)'
      write(unit_output,'(T3,A/)')  ':: S. D. Folkestad, E. F. Kjønstad, May 2017'
      write(unit_output,'(T3,A,1X,A/)') &
                                    'Requested the ground state for:', wf%name
!
      write(unit_output,'(T3,A)')   'Iter.    Norm of amplitude eq.'
      write(unit_output,'(T3,A)')   '------------------------------'    
!
!     Make sure the initial energy is up to date 
!
      call wf%calc_energy
!
!     Open DIIS files 
!
      call generate_unit_identifier(unit_dt)
      open(unit=unit_dt,file='diis_dt',status='unknown',form='unformatted')
!
      call generate_unit_identifier(unit_t_dt)
      open(unit=unit_t_dt,file='diis_t_dt',status='unknown',form='unformatted')
!
      call generate_unit_identifier(unit_diis_matrix)
      open(unit=unit_diis_matrix,file='diis_matrix',status='unknown',form='unformatted')
!
!     Enter iterative loop
!
      iteration = 1
!
      call cpu_time(first_time)
!
      do while ((.not. converged) .and. (iteration .le. wf%settings%ampeqs_max_iterations))
!
!        Save the previous energy 
!
         prev_energy = wf%energy 
!
!        Update the energy 
!
         call wf%calc_energy
!
!        Update the Fock matrix 
!
         call wf%construct_fock
!
!        Construct the current amplitude equations vector,
!        and calculate the norm of the amplitude equations
!
         call wf%calc_ampeqs
         call wf%calc_ampeqs_norm(ampeqs_norm)
!
!        Check for convergence of the energy and the amplitude equations
!
         converged_energy = abs(wf%energy-prev_energy) .lt. wf%settings%energy_threshold
         converged_ampeqs = ampeqs_norm                .lt. wf%settings%ampeqs_threshold
!
!        Print information to output 
!
         write(unit_output,'(T3,I2,7X,E10.4)') iteration, ampeqs_norm 
!
!        Perform DIIS update if convergence hasn't been reached
!
         if (converged_energy .and. converged_ampeqs) then
!
            converged = .true.
!
            write(unit_output,'(/T3,A,I2,A/)') 'Converged in ', iteration, ' iterations!'
            write(unit_output,'(T3,A,F14.8/)') 'Total energy:', wf%energy
!
         else
!
!           Update the amplitudes for the next iteration 
!
            call wf%new_amplitudes
            iteration = iteration + 1
!
         endif
!
      enddo
!
!     Close the DIIS files
!
      close(unit_dt)
      close(unit_t_dt)
      close(unit_diis_matrix)
!
      call cpu_time(last_time)
!
      write(unit_output,*) 'Total time:',last_time-first_time
!
   end subroutine ground_state_solver_ccs
!
!
   module subroutine calc_ampeqs_ccs(wf)
!
!     Calculate Amplitude Equations (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Constructs the amplitude equations vector (the projection vector 
!     in CCS) for the amplitudes of the current iteration of the ground 
!     state solver. It also calculates the norm of the amplitude equations, 
!     which is zero when the equations are exactly solved.
!
      implicit none 
!
      class(ccs) :: wf 
!
!     Update the projection vector 
!
      call wf%construct_omega
!
   end subroutine calc_ampeqs_ccs
!
!
   module subroutine calc_ampeqs_norm_ccs(wf, ampeqs_norm)
!
!     Calculate Amplitude Equations Norm (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
      implicit none 
!
      class(ccs) :: wf 
!
      real(dp) :: ampeqs_norm 
!
      real(dp) :: ddot ! For dot product
!
      ampeqs_norm = zero
      ampeqs_norm = ddot(wf%n_t1am, wf%omega1, 1, wf%omega1, 1)
      ampeqs_norm = sqrt(ampeqs_norm)
!
   end subroutine calc_ampeqs_norm_ccs
!
!
   module subroutine new_amplitudes_ccs(wf)
!
!     New Amplitudes (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Directs the calculation of the quasi-Newton estimate Δ t_i, 
!     and t_i + Δ t_i, and calls the DIIS routine to save & get 
!     the amplitudes for the next iteration.
!
      implicit none 
!
      class(ccs) :: wf 
!
      real(dp), dimension(:,:), allocatable :: dt   ! Δ t_i
      real(dp), dimension(:,:), allocatable :: t_dt ! t_i + Δ t_i
!
      integer(i15) :: n_variables = 0
!
      n_variables = wf%n_t1am 
!
!     Allocate Δ t_i and t_i + Δ t_i vectors 
! 
      call allocator(dt, n_variables, 1)
      call allocator(t_dt, n_variables, 1)
!
      dt   = zero 
      t_dt = zero 
!
!     Calculate Δ t_i
!
      call wf%calc_quasi_Newton_singles(dt, n_variables)
!
!     Set t_i + Δ t_i 
!
      call dcopy(wf%n_t1am, dt, 1, t_dt, 1)           ! t_dt = Δ t_i 
      call daxpy(wf%n_t1am, one, wf%t1am, 1, t_dt, 1) ! t_dt = t_i + Δ t_i
!
!     Save estimates to file and get the next amplitudes
!     (they are placed in dt on exit from diis) 
!
      call wf%diis(dt, t_dt, n_variables)
!
!     Set the new amplitudes 
!
      call dcopy(wf%n_t1am, dt, 1, wf%t1am, 1)
!
!     Deallocate vectors 
!
      call deallocator(dt, n_variables, 1)
      call deallocator(t_dt, n_variables, 1)
!
   end subroutine new_amplitudes_ccs
!
!
   module subroutine calc_quasi_Newton_singles_ccs(wf,dt,n_variables)
!
!     Calculate quasi-Newton estimate (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Calculates the quasi-Newton estimate Δ t_i (singles part)
!     and places the contribution in the dt vector (of length n_variables,
!     with singles first, then doubles, etc. if inherited)
!
      implicit none 
!
      class(ccs) :: wf 
!
      integer(i15), intent(in) :: n_variables
      real(dp), dimension(n_variables, 1) :: dt
!
      integer(i15) :: a = 0, i = 0, ai = 0
!
!     Calculate the singles Δ t_i contribbution
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = index_two(a, i, wf%n_v)
!
            dt(ai, 1) = - wf%omega1(a, i)/(wf%fock_diagonal(wf%n_o + a, 1) - &
                                             wf%fock_diagonal(i, 1))
!
         enddo
      enddo
!
   end subroutine calc_quasi_Newton_singles_ccs
!
!
   module subroutine diis_ccs(wf, dt, t_dt, n_variables)
!
!     DIIS routine (CCS)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     The next amplitudes are 
!
!        t_n+1 = sum_k w_k (t_k + dt_k), 
! 
!     where the weights w_k in front of the quasi-Newton estimate dt_k
!     are determined so as to minimize 
!
!        f(w_k) = sum_k w_k dt_k, 
!
!     with the constraint that g(w_k) = sum_k w_k - 1 = 0.
!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      integer(i15), intent(in) :: n_variables
!
      real(dp), dimension(n_variables, 1) :: dt 
      real(dp), dimension(n_variables, 1) :: t_dt 
!
      real(dp), dimension(:,:), allocatable :: dt_i ! To hold previous Δ t_i temporarily
!
      real(dp) :: ddot
!
      integer(i15) :: i = 0, j = 0
!
      integer      :: info = -1         ! Error integer for dgesv routine (LU factorization)
      integer(i15) :: current_index = 0 ! Progressing as follows: 1,2,...,7,8,1,2,...
!
      real(dp), dimension(:,:), allocatable :: diis_vector
      real(dp), dimension(:,:), allocatable :: diis_matrix
!
      integer(i15), dimension(diis_dim) :: ipiv = 0 ! Pivot integers (see dgesv routine)
!
!     Set the current index 
!
      current_index = iteration - (diis_dim-1)*((iteration-1)/(diis_dim-1)) ! 1,2,...,7,8,1,2,...
!
!     :: Save (Δ t_i) and (t_i + Δ t_i) to files ::
!
      if (current_index .eq. 1) then  
         rewind(unit_dt)
         rewind(unit_t_dt)
      endif
!
      write(unit_dt)   (dt(i,1), i = 1, n_variables)
      write(unit_t_dt) (t_dt(i,1), i = 1, n_variables)
!
!     :: Solve the least squares problem, G * w = H ::
!
!        G : DIIS matrix, G_ij = Δ t_i Δ t_j,
!        H : DIIS vector,  H_i = 0,
!
!     where i, j = 1, 2, ..., current_index. To enforce normality
!     of the solution, G is extended with a row & column of -1's 
!     and H with a -1 at the end.
!
!     First set the DIIS vector to one 
!
      call allocator(diis_vector,current_index+1,1)
      diis_vector = zero 
!
!     Allocate the DIIS matrix and read in previous matrix elements
!
      call allocator(diis_matrix, current_index+1, current_index+1)
      diis_matrix = zero 
!
      if (current_index .gt. 1) then 
!
         rewind(unit_diis_matrix)
!
         do j = 1, current_index - 1
            do i = 1, current_index - 1
!
               read(unit_diis_matrix) diis_matrix(i,j)
!
            enddo
         enddo
!
      endif
!
!     Get the parts of the DIIS matrix G not constructed in 
!     the previous iterations 
!
      call allocator(dt_i, n_variables, 1) ! Allocate temporary holder of quasi-Newton estimates
      dt_i = zero 
!
      rewind(unit_dt)
!
      do i = 1, current_index
!
         read(unit_dt) (dt_i(j,1), j = 1, n_variables) 
!
         diis_matrix(current_index,i) = ddot(n_variables, dt, 1, dt_i, 1) 
         diis_matrix(i,current_index) = diis_matrix(current_index,i)
!
         diis_matrix(current_index+1,i) = -one
         diis_matrix(i,current_index+1) = -one 
!
      enddo
!
      diis_vector(current_index+1,1) = -one
!
!     Write the current DIIS matrix to file 
!
      rewind(unit_diis_matrix)
!
      do j = 1, current_index
         do i = 1, current_index
!
            write(unit_diis_matrix) diis_matrix(i,j)
!
         enddo
      enddo
!
!     Solve the DIIS equation 
!
!     Note: on exit, the solution is in the diis_vector,
!     provided info = 0 (see LAPACK documentation for more)
!
      call dgesv(current_index+1,  &
                  1,               &
                  diis_matrix,     &
                  current_index+1, &
                  ipiv,            &
                  diis_vector,     &
                  current_index+1, &
                  info)
!
!     :: Update the amplitudes (placed in dt on exit) ::
!
      call dzero(dt, n_variables)
!
      rewind(unit_t_dt)
!
      do i = 1, current_index
!
!        Read the t_i + Δ t_i vector 
!
         call dzero(t_dt, n_variables)
         read(unit_t_dt) (t_dt(j, 1), j = 1, n_variables)
!
!        Add w_i (t_i + Δ t_i) to the amplitudes 
!
         call daxpy(n_variables, diis_vector(i, 1), t_dt, 1, dt, 1)
!
      enddo
!
!     Deallocations 
!
      call deallocator(dt_i, n_variables, 1)
      call deallocator(diis_vector, current_index + 1, 1)
      call deallocator(diis_matrix, current_index + 1, current_index+1)
!
   end subroutine diis_ccs 
!
!
end submodule diis