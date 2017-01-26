module mlcc_workspace
!
!
!  mlcc mod work definitions
!  Authors Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad
!  January 2017
!
!  Purpose: define pointers needed for memory management
!
   use mlcc_types
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
subroutine work_init(mem,lupri)
!
!  Authors Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad
!  January 2017
!
!  Purpose: set up memory management
!
   implicit none
!
   integer, intent(in)                    :: mem
   integer, intent(in)                    :: lupri
!
   write(lupri,*)
   write(lupri,*) 'In work_int'
   write(lupri,*)
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
!
!
!  Authors Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad
!  January 2017
!
!  Purpose: allocation and update of memory info
!
   implicit none
!  
   real(dp), dimension(:,:), pointer        :: elm
   integer, intent(in)                      :: M,N
   integer                                  :: size
   integer                                  :: stat, error
!
   size = M*N
!
   allocate(elm(M,N),stat = error)
!
   if (stat .ne. 0) then
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
!
!
subroutine deallocator(elm,M,N)
!
!
!  Authors Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad
!  January 2017
!
!  Purpose: dallocation and update of memory info
!
   implicit none
!
   real(dp), dimension(:,:), pointer       :: elm
   integer                                 :: stat, error
   integer, intent(in)                     :: M, N
   integer                                 :: size
!
   size = M*N
!
   deallocate(elm,stat = error)  
   if (stat .ne. 0) then
      print*,"error: couldn't deallocate array"
      stop
   endif
! 
   work_remains = work_remains+size
   work_used = work_used-size
end subroutine deallocator
!
!
end module mlcc_work