module printing
!
   integer :: unit_output
!
contains
   subroutine init_output_file
!
      implicit none
!
      integer :: unit_output = 1
      open(unit=unit_output,file='mlcc.out',status='new',form='formatted')
      close(unit=unit_output)
!
   end subroutine init_output_file
!
end module printing