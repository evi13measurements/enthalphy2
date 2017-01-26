module mlcc_utilities
function index_t2(nai,nbj)
!   
   implicit none
!
   integer, intent(in)     :: nai, nbj
   integer, intent(out)    :: index_t2

   index_t2 = max(nai,nbj)*(max(nai,nbj)-3)/2+nai+nbj

end function


end module mlcc_utilities