submodule (cc2_class) omega
!
!
!                       -::- Omega submodule (CC2) -::-
!           Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!
!
!     Contains the following family of procedures of the CC2 class:
!
!        initialize_omega: allocates the projection vector omega1
!                          and sets it to zero.
!
!        construct_omega:  constructs the projection vector omega1
!                          for the current amplitudes t1am for the
!                          wavefunction object wf. The routine assumes that
!                          the projection vector is allocated.
!
!        omega_a1:         adds A1 term to omega1
!        omega_b1:         adds B1 term to omega1
!        omega_c1:         adds C1 term to omega1
!        omega_d1:         adds D1 term to omega1
!
   implicit none 
!
   logical :: debug = .false.
!
!
contains
!
end submodule