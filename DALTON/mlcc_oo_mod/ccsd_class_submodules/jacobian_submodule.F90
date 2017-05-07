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
!     coupled cluster Jacobian matrix 
!
!        A_mu,nu = < mu | [e^(-T) H e^(T),tau_nu] | R >.
!
!     On exit, A*c is placed in the incoming c vector.
!
      class(ccsd) :: wf 
!
      real(dp), dimension(:,:) :: c1am 
      real(dp), dimension(:,:) :: c2am 
!
      real(dp), dimension(:,:), allocatable :: tr1am ! Transformed vector, singles 
      real(dp), dimension(:,:), allocatable :: tr2am ! Transformed vector, doubles 
!
!     Calculate the singles contributions 
!
     ! call wf%jacobian_a1(c1am,c2am)
!
   end subroutine jacobian_transformation_ccsd
!
!
end submodule jacobian 