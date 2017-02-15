module mlcc_utilities

use mlcc_data
contains
integer function index_t2(nai,nbj)
!
! Purpose: Calculates t2 indices.
!
!   
   implicit none
!
   integer, intent(in)     :: nai, nbj
!
   index_t2 = max(nai,nbj)*(max(nai,nbj)-3)/2+nai+nbj
!
end function index_t2
!
!
integer function index_t1(i,a)
!
! Purpose: Calculates t1 indices.
!
!   
   implicit none
!
   integer, intent(in)     :: i, a
!
   index_t1 = n_vir*(i-1) + a
!
end function index_t1
!
integer function packed_size()
! Purpose: Returns size of packed vectors alpha >= beta
!
    use mlcc_data
!   
    implicit none
!
    packed_size = n_orbitals*(n_orbitals+1)/2
end function packed_size
subroutine squareup(packed,unpacked)
!
! Purpose: Squareup of packed vectors.
! 
    use mlcc_data
    use mlcc_types
!
    implicit none
!
!
    real(dp),intent(in)     ::    packed
    real(dp)                ::    unpacked

end subroutine
end module mlcc_utilities