module ccsd_class
!
!  Coupled cluster singles and doubles (CCSD) 
!
   use ccs_class
!
   implicit none 
!
   type, extends(cc_singles) :: cc_singles_doubles
!
      integer(i15), private                          :: n_t2am = 0
      real(dp), dimension(:,:), allocatable, private :: t2am
!
   contains
!
      procedure, public :: init => init_cc_singles_doubles
      procedure, public :: drv  => drv_cc_singles_doubles
!
      procedure, private :: initialize_doubles => initialize_doubles_cc_singles_doubles
!
   end type cc_singles_doubles
!
contains
!
   subroutine init_cc_singles_doubles(wavefn)
!
      implicit none 
!
      class(cc_singles_doubles) :: wavefn
!
      write(unit_output,*) 'In init_cc_singles_doubles'
!
!     Read Hartree-Fock info from SIRIUS
!
      wavefn % read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      wavefn % read_transform_cholesky 
!
!     Initialize singles amplitudes
!
      wavefn % initialize_singles 
!
!     Initialize doubles amplitudes 
!
      wavefn % initialize_doubles
!
   end subroutine init_cc_singles_doubles
!
   subroutine drv_cc_singles_doubles(wavefn)
!
      implicit none 
!
      class(cc_singles_doubles) :: wavefn
!
      write(unit_output,*) 'In drv_cc_singles_doubles'
!
!     Call the CCS driver
!
!        This driver handles calculations that are not specific
!        to CCSD, such as the amplitude equations
!
      call drv_cc_singles(wavefn)
!
   end subroutine drv_cc_singles_doubles
!
   subroutine initialize_doubles_cc_singles_doubles(wavefn)
!
      implicit none 
!
      class(cc_singles_doubles) :: wavefn
!
!     Calculate the number of singles amplitudes
!
      wavefn % n_t2am = (wavefn % n_t1am)*(wavefn % n_t1am + 1)/2
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator ( wavefn % t2am, wavefn % n_t2am, 1)
      wavefn % t2am = zero
!
   end subroutine initialize_doubles_cc_singles_doubles
!
end module ccsd_class
