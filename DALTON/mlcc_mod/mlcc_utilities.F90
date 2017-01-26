module mlcc_utilities

use mlcc_data
contains
integer function index_t2(nai,nbj)
!   
   implicit none
!
   integer, intent(in)     :: nai, nbj
!
   index_t2 = max(nai,nbj)*(max(nai,nbj)-3)/2+nai+nbj
!
end function index_t2
!
integer function index_t1(i,a)
!   
   implicit none
!
   integer, intent(in)     :: i, a
!
   index_t1 = n_virt*(i-1) + a
!
end function index_t1
!
end module mlcc_utilities