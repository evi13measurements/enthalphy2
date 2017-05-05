module cc2_class
!
!
!
!               Coupled cluster singles (CC2) class module                                 
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017         
!                                                                           
!
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
!  General tools
!
   use types
   use utils
   use workspace
   use input_output
!
!  The ancestor class module (ccs)
!
   use ccs_class
!
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CC2 class -::-
!  ::::::::::::::::::::::::::::::::::::: 
!
   type, extends(ccs) :: cc2
!
      real(dp), dimension(:,:), allocatable :: omega1
   contains 
!
!     Initialization and driver routines
!
      procedure :: init => init_cc2
      procedure :: drv  => drv_cc2
      procedure :: initialize_omega => initialize_omega_cc2
!      
   end type cc2
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCS -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface
   !
   !
   end interface
!
contains
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::
!
   subroutine init_cc2(wf)
!
!     Initialize CC2 object
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!
!     Performs the following tasks
!
!     1. Sets HF orbital and energy information by reading from file (read_hf_info)
!     2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!     3. Allocates the Fock matrix and sets it to zero (note: the matrix is constructed in 
!        the descendant classes) 
!     4. Allocates the singles amplitudes and sets them to zero, and sets associated properties 
!     5. Allocate Omega vector
!
      implicit none
!
      class(cc2)  :: wf
!
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
!
!     Initialize amplitudes and associated attributes
!
      call wf%initialize_amplitudes
!
!     Allocate Fock matrix and set to zero
!
      call wf%initialize_fock_matrix
!
!     Initialize omega vector
!
      call wf%initialize_omega
!
   end subroutine init_cc2
!
   subroutine drv_cc2(wf)
!
!     CC2 Driver
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!
      implicit none 
!
      class(cc2) :: wf
!
   end subroutine drv_cc2
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
!
end module cc2_class