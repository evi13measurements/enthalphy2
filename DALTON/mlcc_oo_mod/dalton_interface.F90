module dalton_interface
!
!   use cholesky_integrals_class
contains
!
   subroutine dalton_interface_driver()
!
!  Dalton interface driver
!
!
      use input_output
!
      implicit none
!
!      class(cholesky_integrals) :: cholesky
      integer                   :: unit_cholesky_ao
!
!     Setting n_J      
!
      call generate_unit_identifier(unit_cholesky_ao)
!
!      write(unit_output,*)'Cholesky unit identifier', unit_cholesky_ao  
   end subroutine dalton_interface_driver
!
end module dalton_interface