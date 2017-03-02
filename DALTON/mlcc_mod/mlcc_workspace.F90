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
   use mlcc_data
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
!  Authors Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad
!  January 2017
!
!  Purpose: set up memory management
!
   implicit none
!
   integer, intent(in)                    :: mem
!
   write(luprint,*)
   write(luprint,*) 'In work_int'
   write(luprint,*)
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
   integer, intent(in)                      :: M,N
   real(dp), dimension(:,:), pointer        :: elm
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
   work_remains = work_remains-4*size
   work_used = work_used+4*size
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
   work_remains = work_remains+4*size
   work_used = work_used-4*size
!
end subroutine deallocator
!
!
subroutine allocator_int(elm,M,N)
!
!
!  Authors Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad
!  January 2017
!
!  Purpose: allocation and update of memory info
!
   implicit none
!  
   integer, intent(in)                      :: M,N
   integer, dimension(:,:), pointer         :: elm
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
   work_remains = work_remains-2*size
   work_used = work_used+2*size
   if (work_remains .lt. 0) then
      print*,"error: User specified memory too small"
      stop
   endif
end subroutine allocator_int
!
!
subroutine deallocator_int(elm,M,N)
!
!
!  Authors Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad
!  January 2017
!
!  Purpose: dallocation and update of memory info
!
   implicit none
!
   integer, dimension(:,:), pointer        :: elm
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
   work_remains = work_remains+2*size
   work_used = work_used-2*size
end subroutine deallocator_int
!
subroutine vec_allocator(elm,N)
!
!  (Eirik: These vector allocators will not be necessary... To remove...)
!  Authors Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad
!  January 2017
!
!  Purpose: allocation of vector and update of memory info
!
   implicit none
!  
   integer, intent(in)                      :: N
   real(dp), dimension(:), pointer          :: elm
   integer                                  :: size
   integer                                  :: stat, error
!
   size = N
!
   allocate(elm(N),stat = error)
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
!
end subroutine vec_allocator
!
subroutine vec_deallocator(elm,N)
!
!
!  Authors Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad
!  January 2017
!
!  Purpose: dallocation of vector and update of memory info
!
   implicit none
!
   real(dp), dimension(:), pointer         :: elm
   integer                                 :: stat, error
   integer, intent(in)                     :: N
   integer                                 :: size
!
   size = N
!
   deallocate(elm,stat = error)  
   if (stat .ne. 0) then
      print*,"error: couldn't deallocate array"
      stop
   endif
! 
   work_remains = work_remains+size
   work_used = work_used-size
!
end subroutine vec_deallocator
!
integer function get_available()
   get_available=work_remains
end function get_available

end module mlcc_workspace