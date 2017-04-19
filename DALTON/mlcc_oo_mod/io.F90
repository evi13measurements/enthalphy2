module io
!
   implicit none
!
   private
   integer  :: n_files = 0
!
contains
!
   subroutine generate_unit_identifier(unit_identifier)
!
      implicit none
!
      integer :: unit_identifier
!
      n_files = n_files + 1
      unit_identifier = n_files
!
   end subroutine generate_unit_identifier
!
end module io