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
   index_t2 = (max(nai,nbj)*(max(nai,nbj)-3)/2)+nai+nbj
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
integer function packed_size(N)
! Purpose: Returns size of packed vectors alpha >= beta
!
    use mlcc_data
!   
    implicit none
    integer         ::   N
!
    packed_size = N*(N+1)/2
!
end function packed_size
!
!
!
subroutine squareup(packed,unpacked,N)
!
! Purpose: Squareup to full dimension (NxN) of packed vectors.
! 
    use mlcc_data
    use mlcc_types
!
    implicit none
!
!
    real(dp),dimension(:,:),intent(in)     ::    packed
    real(dp),dimension(:,:)                ::    unpacked
    integer                                ::    i,j,N
!
    do i = 1,N
        do j = i+1,N
            unpacked(i,j)=packed(index_t2(i,j),1)
            unpacked(j,i)=unpacked(i,j)
        enddo
        unpacked(i,i)=packed(index_t2(i,i),1)
    enddo
!
end subroutine
!
!
subroutine packin(packed,unpacked,N)
!
! Purpose: Pack down of full square matrix of dimension (NxN).
! 

    use mlcc_data
    use mlcc_types
!
    implicit none
!
!
    real(dp),dimension(:,:)              ::    packed
    real(dp),dimension(:,:),intent(in)   ::    unpacked
    integer                              ::    i,j,N
!
    do i = 1,N
        do j = i,N
          packed(index_t2(i,j),1)=unpacked(i,j)
        enddo
    enddo
!
end subroutine
!
integer function three_i(p,q,r,dim_p,dim_q)
!
!   Three index integer function
!
!   Calculates the compound index (pqr)
!
    implicit none
!
    integer :: p,q,r,dim_p,dim_q
!
    three_i = dim_p*(dim_q*(r-1)+q-1)+p
!
end function three_i
!
integer function two_i(p,q,dim_p)
!
!   Two index integer function
!
!   Calculates the compound 
!
    implicit none
!
    integer :: p,q,dim_p
!
    two_i = dim_p*(q-1)+p
!
end function two_i
!
end module mlcc_utilities