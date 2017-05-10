module workspace
!
!  Workspace 
!  Written by Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad, Jan 2017
!
!  A collection of subroutines for memory management and monitoring
!
   use types
!   
   implicit none
!
   integer, private :: work_length  = 0
   integer, private :: work_remains = 0
   integer, private :: work_used    = 0
   integer(i15)     :: mem          = 10000000000 ! Sarai: This should be set from user input
!
contains
!
   subroutine work_init
!
!     Work initilization 
!     Written by Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad, Jan 2017
!  
!     Sets up memory management
!
      implicit none
!  
      work_length  = mem
      work_remains = mem
      work_used    = 0  
!   
   end subroutine work_init
!
   subroutine allocator(elm,M,N)
!
!     Allocator 
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!  
!     Allocates array and updates memory variables
!  
      implicit none
!     
      integer, intent(in)                    :: M,N
      real(dp), dimension(:,:), allocatable  :: elm
      integer                                :: size
      integer                                :: stat, error
!  
      size = M*N
!  
      allocate(elm(M,N),stat = error)
!  
      if (stat .ne. 0) then
         print*,"Error: couldn't allocate memory for array of size",size
         stop
      endif
!       
      work_remains = work_remains - 4*size
      work_used    = work_used + 4*size
!
      if (work_remains .lt. 0) then
         print*,"Error: user-specified memory too small"
         stop
      endif
!
   end subroutine allocator
!
   subroutine deallocator(elm,M,N)
!
!     Deallocator 
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!  
!     Deallocation and update of memory information
!  
      implicit none
!  
      real(dp), dimension(:,:), allocatable   :: elm
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
      work_remains = work_remains + 4*size
      work_used    = work_used - 4*size
!
   end subroutine deallocator
!
   subroutine allocator_int(elm,M,N)
!
!     Allocator for integer arrays  
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!  
!     Allocation and update of memory information
!  
      implicit none
!     
      integer, intent(in)                  :: M,N
      integer, dimension(:,:), allocatable :: elm
      integer                              :: size
      integer                              :: stat, error
!  
      size = M*N
!  
      allocate(elm(M,N),stat = error)
!  
      if (stat .ne. 0) then
         print*,"Error: couldn't allocate memory for array of size",size
         stop
      endif 
!  
      work_remains = work_remains-2*size
      work_used = work_used+2*size
      if (work_remains .lt. 0) then
         print*,"Error: user-specified memory too small"
         stop
      endif
!
   end subroutine allocator_int
!
   subroutine deallocator_int(elm,M,N)
!
!     Deallocator for integer arrays 
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!  
!     Deallocation and update of memory information
!  
      implicit none
!  
      integer, dimension(:,:), allocatable :: elm
      integer                              :: stat, error
      integer, intent(in)                  :: M, N
      integer                              :: size
!  
      size = M*N
!  
      deallocate(elm,stat = error)  
      if (stat .ne. 0) then
         print*,"error: couldn't deallocate array"
         stop
      endif
!  
      work_remains = work_remains + 2*size
      work_used    = work_used - 2*size
!
   end subroutine deallocator_int
!
   integer function get_available()
!
!     Get available memory 
!     Written by Eirik F. Kjønstad and Sarai F. Folkestad, Mar 2017
!  
!     Returns the available memory 
!  
      get_available = work_remains
!
   end function get_available
!
end module workspace
