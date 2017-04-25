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
      integer(i15) :: n_t2am = 0
      real(dp), dimension(:,:), allocatable :: t2am
!
   contains
!
      procedure :: init => init_cc_singles_doubles
      procedure :: drv  => drv_cc_singles_doubles
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
!     Initializing CCS-specific quantities
!
      call init_cc_singles(wavefn)
!
! Stuff...  
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
end module ccsd_class
