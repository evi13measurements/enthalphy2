module mlcc_utilities
!
!   MLCC Utilities 
!   Written by Sarai D. Folkstad and Eirik F. Kjønstad, 28 Feb 2017
!
    use mlcc_data
!
contains
   integer function index_packed(i,j)
!
!     Purpose: Returns packed index (ij) for symmetric arrays
! 
      implicit none
!
      integer, intent(in) :: i,j
!
      index_packed = (max(i,j)*(max(i,j)-3)/2) + i + j
!
   end function index_packed
!
   integer function packed_size(N)
!
!     Purpose: Returns size of packed vectors alpha >= beta
!
      use mlcc_data
!   
      implicit none
!
      integer :: N
!
      packed_size = N*(N+1)/2
!
   end function packed_size
!
   subroutine squareup(packed,unpacked,N)
!
!     Purpose: Squares up to full dimension (N x N) of packed vectors.
! 
      use mlcc_data
      use mlcc_types
!
      implicit none
!
      real(dp),dimension(:,:),intent(in)     ::    packed
      real(dp),dimension(:,:)                ::    unpacked
      integer                                ::    i,j,N
!
      do i=1,N
         do j=1,N
            unpacked(i,j) = packed(index_packed(i,j),1)
         enddo
      enddo
!
      ! do i = 1,N
      !    do j = i+1,N
      !       unpacked(i,j)=packed(index_packed(i,j),1) ! Eirik: I am rewriting this.
      !       unpacked(j,i)=unpacked(i,j)
      !    enddo
      !    unpacked(i,i)=packed(index_packed(i,i),1)
      ! enddo
!
   end subroutine
!
   subroutine packin(packed,unpacked,N)
!
!     Purpose: Pack down full square matrix of dimension (N x N).
! 
      use mlcc_data
      use mlcc_types
!
      implicit none
!
      real(dp),dimension(:,:)              ::    packed
      real(dp),dimension(:,:),intent(in)   ::    unpacked
      integer                              ::    i,j,N
!
      do i = 1,N
         do j = i,N
            packed(index_packed(i,j),1)=unpacked(i,j)
         enddo
      enddo
!
   end subroutine
!
   integer function index_three(p,q,r,dim_p,dim_q)
!
!     Purpose: Returns the compound index (pqr)
!
      implicit none
!
      integer :: p,q,r,dim_p,dim_q
!
      index_three = dim_p*(dim_q*(r-1)+q-1)+p
!
   end function index_three
!
   integer function index_two(p,q,dim_p)
!
!     Purpose: Returns the compound index (pq)
!
      implicit none
!
      integer :: p,q,dim_p
!
      index_two = dim_p*(q-1)+p
!
   end function index_two
!
end module mlcc_utilities