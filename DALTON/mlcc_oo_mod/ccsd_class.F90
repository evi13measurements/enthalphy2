module ccsd_class
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!                                                                 
!      Coupled cluster singles and doubles (CCSD) class module                                 
!   Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017  
!                                                                 
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Use general tools 
!
   use types
   use utils
   use workspace
   use input_output
!
!  Use ancestor class module (CCS)
!
   use ccs_class
!
   implicit none 
!
!  Definition of the CCSD class 
!
   type, extends(cc_singles) :: cc_singles_doubles
!
!     Amplitude attributes
!
      integer(i15) :: n_t2am = 0                    ! Number of doubles amplitudes
      real(dp), dimension(:,:), allocatable :: t2am ! Doubles amplitude vector
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_cc_singles_doubles
      procedure :: drv  => drv_cc_singles_doubles
!
!     Initialization routine for the (singles, doubles) amplitudes
!
      procedure :: initialize_amplitudes => initialize_amplitudes_cc_singles_doubles
!
   end type cc_singles_doubles
!
contains
!
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!              Initialization and driver routines
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
   subroutine init_cc_singles_doubles(wfn)
!
!     Initialize CCSD object
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Performs the following tasks
!
!     1. Sets HF orbital and energy information by reading from file (read_hf_info)
!     2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!     3. Allocates the Fock matrix and sets it to zero (note: the matrix is constructed in 
!        the descendant classes)
!     4. Allocates the singles and doubles amplitudes and sets them to zero, and sets associated
!        properties
!
      implicit none 
!
      class(cc_singles_doubles) :: wfn

!     Read Hartree-Fock info from SIRIUS
!
      call wfn % read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wfn % read_transform_cholesky 
!
!     Initialize (singles and doubles) amplitudes
!
      call wfn % initialize_amplitudes
!
   end subroutine init_cc_singles_doubles
!
!
   subroutine drv_cc_singles_doubles(wfn)
!
      implicit none 
!
      class(cc_singles_doubles) :: wfn
!
      write(unit_output,*) 'In drv_cc_singles_doubles'
!
!     Call the CCS driver
!
!        This driver handles calculations that are not specific
!        to CCSD, such as the amplitude equations
!
      call drv_cc_singles(wfn)
!
   end subroutine drv_cc_singles_doubles
!
!
   subroutine initialize_amplitudes_cc_singles_doubles(wfn)
!
      implicit none 
!
      class(cc_singles_doubles) :: wfn
!
!     Calculate the number of singles and doubles amplitudes
!
      wfn % n_t1am = (wfn % n_occ)*(wfn % n_vir) 
      wfn % n_t2am = (wfn % n_t1am)*(wfn % n_t1am + 1)/2
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator(wfn % t1am, wfn % n_t1am, 1)
      wfn % t1am = zero
!
!     Allocate the doubles amplitudes and set to zero
!
      call allocator (wfn % t2am, wfn % n_t2am, 1)
      wfn % t2am = zero
!
   end subroutine initialize_amplitudes_cc_singles_doubles
!
!
end module ccsd_class
