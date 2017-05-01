submodule (ccsd_class) omega
!
!
!                       -::- Cholesky submodule (CCS) -::-
!           Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!
!     Contains the following family of procedures of the CCSD class:
!
!     construct_omega(wfn): constructs the projection vector for the current amplitudes
!
contains
!
!
      subroutine construct_omega_cc_singles_doubles(wfn)
!
         implicit none 
!
         class(cc_singles_doubles) :: wfn 
!
         ! call omega_a1
         ! call omega_b1 ...
!
      end subroutine construct_omega_cc_singles_doubles
!
!
end submodule omega