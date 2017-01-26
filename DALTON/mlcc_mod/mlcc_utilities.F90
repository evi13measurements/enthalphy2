module mlcc_utilities
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


end module mlcc_utilities