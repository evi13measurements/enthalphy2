module mlcc_mod_work
!
!
!  mlcc work definitions
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: define pointers needed for memory management
!
!  use other modules?
!   
   implicit none
!
   integer, private                       :: work_length = 0
   integer, private                       :: work_remains = 0
   integer, private                       :: work_used = 0
!
!
contains
!
subroutine work_init(mem)
!
   implicit none
!
   integer, intent(in)                    :: mem 
!
   work_length = mem
   work_remains = mem
   work_used = 0
!
!   
end subroutine work_init
!
!
subroutine allocator(elm,M,N)
   implicit none
   real*8, pointer                        :: elm(:,:)
   integer, intent(in)                    :: M,N
   integer                                :: size
   integer                                :: stat
   integer                                :: error
!
   size = M*N
!
   allocate(elm(M,N),stat=error)
!
   if (stat.ne.0) then
      print*,"error: couldn't allocate memory for array, size=",size
      stop
   endif
!    
!
   work_remains = work_remains-size
   work_used = work_used+size
   if (work_remains .lt. 0) then
      print*,"error: User specified memory too small"
      stop
   endif
end subroutine allocator
end module mlcc_mod_work