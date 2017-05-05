module input_output
! Write doc!!
   use types
!
   implicit none
!
   integer(i15) :: unit_output = 0 
   integer, private :: n_files = 200 ! To not overwrite the DALTON.OUT identifier (for debug against old code)
!
contains
!
   subroutine generate_unit_identifier(unit_identifier)
!
      implicit none
!
      integer(i15) :: unit_identifier
!
      n_files = n_files + 1
      unit_identifier = n_files
!
   end subroutine generate_unit_identifier
!
   subroutine init_output_file
!
      implicit none
!
      call  generate_unit_identifier(unit_output)
!
      open(unit=unit_output,file='mlcc.out',status='unknown',form='formatted')
      close(unit=unit_output)
!
   end subroutine init_output_file
!
end module input_output
