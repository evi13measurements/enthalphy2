submodule (ccsd_class) diis
!
!
!                             -::- DIIS submodule (CCSD) -::-
!                 Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!
!     Consists of the following subroutines of the CCSD module:
! 
!     new_amplitudes:             Calculates the quasi-Newton estimate and passes the 
!                                 information needed by the DIIS routine.
!     calc_ampeqs_norm:           Calculates the norm of the amplitude equations.
!     calc_quasi_Newton_doubles:  Calculates the doubles part of the quasi-Newton estimate.
!
!     Can be inherited by models of the same level (e.g. CC3) without modification.
!
!     When inherited by higher level models (e.g. CCSDT), the new_amplitudes and calc_ampeqs_norm
!     routines should be overridden to account for the triples quasi-Newton estimate, amplitudes, 
!     and projection vector.
!
!
   implicit none 
!
!
contains
!
!
   module subroutine calc_ampeqs_norm_ccsd(wf, ampeqs_norm)
!
!     Calculate Amplitude Equations Norm (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp) :: ampeqs_norm 
!
      real(dp) :: ddot ! For dot product
!
      ampeqs_norm = zero
      ampeqs_norm = ddot(wf%n_t1am, wf%omega1, 1, wf%omega1, 1)
      ampeqs_norm = ddot(wf%n_t2am, wf%omega2, 1, wf%omega2, 1) + ampeqs_norm
      ampeqs_norm = sqrt(ampeqs_norm)
!
   end subroutine calc_ampeqs_norm_ccsd
!
!
   module subroutine new_amplitudes_ccsd(wf)
!
!     New Amplitudes (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Directs the calculation of the quasi-Newton estimate Δ t_i, 
!     and t_i + Δ t_i, and calls the DIIS routine to save & get 
!     the amplitudes for the next iteration.
!
      implicit none 
!
      class(ccsd) :: wf 
!
      integer(i15) :: i = 0 ! for debug purposes
!
      real(dp), dimension(:,:), allocatable :: dt   ! Δ t_i
      real(dp), dimension(:,:), allocatable :: t_dt ! t_i + Δ t_i
!
      integer(i15) :: n_variables = 0
!
      n_variables = wf%n_t1am + wf%n_t2am
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
      call wf%calc_quasi_Newton_doubles(dt, n_variables)
!
!     Set t_i + Δ t_i 
!
      call dcopy(n_variables, dt, 1, t_dt, 1) ! t_dt = Δ t_i 
!
      call daxpy(wf%n_t1am, one, wf%t1am, 1, t_dt, 1)                   ! t_dt = t_i + Δ t_i singles 
      call daxpy(wf%n_t2am, one, wf%t2am, 1, t_dt(wf%n_t1am + 1, 1), 1) ! t_dt = t_i + Δ t_i doubles     
!
!     Save estimates to file and get the next amplitudes
!     (they are placed in dt on exit from diis) 
!
      call wf%diis(dt, t_dt, n_variables)
!
!     Set the new amplitudes 
!
      call dcopy(wf%n_t1am, dt, 1, wf%t1am, 1)
      call dcopy(wf%n_t2am, dt(wf%n_t1am + 1, 1), 1, wf%t2am, 1)
!
!     Deallocate vectors 
!
      call deallocator(dt, n_variables, 1)
      call deallocator(t_dt, n_variables, 1)
!
   end subroutine new_amplitudes_ccsd
!
!
   module subroutine calc_quasi_Newton_doubles_ccsd(wf,dt,n_variables)
!
!     Calculate quasi-Newtoni doubles estimate (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Calculates the quasi-Newton estimate Δ t_i (doubbles part)
!     and places the contribution in the dt vector (of length n_variables,
!     with singles first, then doubles, etc. if inherited)
!
      implicit none 
!
      class(ccsd) :: wf 
!
      integer(i15), intent(in) :: n_variables
      real(dp), dimension(n_variables, 1) :: dt
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0
      integer(i15) :: ai = 0, bj = 0, aibj = 0, offset = 0
!
!     Calculate the doubles Δ t_i contribbution
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate the necessary indices 
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj) 
!
                  offset = wf%n_t1am + aibj ! dt has singles first, then doubles 
!
                  dt(offset,1) = - wf%omega2(aibj, 1)/(wf%fock_diagonal(wf%n_o + a, 1) + &
                                                       wf%fock_diagonal(wf%n_o + b, 1) - &
                                                       wf%fock_diagonal(i, 1) -          &
                                                       wf%fock_diagonal(j, 1))
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine calc_quasi_Newton_doubles_ccsd
!
!
end submodule diis 
