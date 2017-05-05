module utils
!
!   Utilities module 
!   Written by Sarai D. Folkstad and Eirik F. Kjønstad, 28 Feb 2017
!
   use input_output
   use types
!
contains
!
!
   integer(i15) function index_packed(i,j)
!
!     Purpose: Returns packed index (ij) for symmetric arrays
! 
      implicit none
!
      integer(i15), intent(in) :: i,j
!
      index_packed = (max(i,j)*(max(i,j)-3)/2) + i + j
!
   end function index_packed
!
!
   integer(i15) function packed_size(N)
!
!     Purpose: Returns size of packed symmetric matrices
!              of dimension N x N (triangular elements) 
!   
      implicit none
!
      integer(i15), intent(in) :: N
!
      packed_size = N*(N+1)/2
!
   end function packed_size
!
!
   subroutine squareup(packed,unpacked,N)
!
!     Purpose: Squares up to full dimension (N x N) 
!              of packed matrices.
!
      implicit none
!
      integer(i15), intent(in) :: N
!
      real(dp), dimension(:,:), intent(in) :: packed
      real(dp), dimension(:,:)             :: unpacked
!
      integer(i15) :: i = 0, j = 0
!
      do i = 1, N
         do j = 1, N
            unpacked(i, j) = packed(index_packed(i,j), 1)
         enddo
      enddo
!
   end subroutine
!
!
   subroutine packin(packed,unpacked,N)
!
!     Purpose: Pack down full square matrix of dimension N x N.
!
      implicit none
!
      integer(i15), intent(in) :: N
!
      real(dp), dimension(:,:) :: packed
      real(dp), dimension(:,:),intent(in) :: unpacked
!
      integer(i15) :: i = 0, j = 0
!
      do i = 1, N
         do j = 1, N
            packed(index_packed(i, j), 1) = unpacked(i, j)
         enddo
      enddo
!
   end subroutine
!
!
   integer(i15) function index_three(p,q,r,dim_p,dim_q)
!
!     Purpose: Returns the compound index (pqr)
!
      implicit none
!
      integer(i15), intent(in) :: p, q, r, dim_p, dim_q
!
      index_three = dim_p*(dim_q*(r-1)+q-1)+p
!
   end function index_three
!
!
   integer(i15) function index_two(p,q,dim_p)
!
!     Purpose: Returns the compound index (pq)
!
      implicit none
!
      integer(i15), intent(in) :: p,q,dim_p
!
      index_two = dim_p*(q-1)+p
!
   end function index_two
!
!
   subroutine num_batch(required,available,max_batch_length,n_batch,batch_dimension)
!  Purpose: Calculate number of batches
!
!  required          =  required memory (words)
!  available         =  available memory (words)
!  max_batch_length  =  length of batch 
!  n_batch           =  number of batches
!  batch_dimension   =  original size of dimension that we batch over
!
!  Batching structure will be:
!  With rest:     (n_batch-1)*(max_batch_length) + rest = required
!  Without rest:  (n_batch)*(max_batch_length) = required
!
      implicit none
!
!     
      integer(i15), intent(in)           :: required, available, batch_dimension
      integer(i15)                       :: max_batch_length,n_batch
!
      if (required .lt. available) then
         n_batch = 1
         max_batch_length = batch_dimension
         return
      endif
!
!  Max batch size
!
      max_batch_length = available/(required/batch_dimension)
!
!  Number of full batches
!
      n_batch=batch_dimension/max_batch_length
!
!  Test for rest
!
      if (n_batch*max_batch_length .lt. required) then
         n_batch = n_batch+1
      endif
!
   end subroutine num_batch
!
   subroutine num_two_batch(required,available,max_batch_length,n_batch,batch_dimension)
!  Purpose: Calculate number of batches when batching over two variables of same dimension.
!
!  required          =  required memory (words)
!  available         =  available memory (words)
!  max_batch_length  =  length of batch 
!  n_batch           =  number of batches
!  batch_dimension   =  original size of dimension that we batch over
!
!  Batching structure will be:
!  With rest:     (n_batch-1)*(max_batch_length) + rest = required
!  Without rest:  (n_batch)*(max_batch_length) = required
!
      implicit none
!
!     
      integer(i15), intent(in)           :: required, available, batch_dimension
      integer(i15)                       :: max_batch_length,n_batch,i
!
   if (required .lt. available) then
         n_batch = 1
         max_batch_length = batch_dimension
         return
   endif
!  
   do i = 1, batch_dimension
      if (available .gt. required/i**2) then
         n_batch = i
         max_batch_length = batch_dimension/n_batch
         return
      endif
   enddo
!
   end subroutine num_two_batch
!
   subroutine batch_limits(first,last,batch_number,max_batch_length,batch_dimension)
!
!     Purpose: Find batch limits (first and last) 
!
!        batch_number: the current batch (1,2,...,n_batch)
!        max_batch_length: the length of each batch (except the last, which may be a rest, see n_one_batch routine)
!        batch_dimension: the dimensionality of the batching variable (e.g., n_vir for a virtual index)
!
      implicit none 
!
      integer(i15) :: first,last
      integer(i15), intent(in) :: batch_number,max_batch_length,batch_dimension
!
      first = 1 + (batch_number-1)*max_batch_length
      last  = min(max_batch_length+(batch_number-1)*max_batch_length,batch_dimension)
!
   end subroutine batch_limits
!
!
   subroutine vec_print(vec,dim_1,dim_2)
!
!     Purpose: prints a vector with a compound index (p q) of dimension (dim_1 x dim_2)
!        (for debugging, remove or replace later on)
!
      implicit none
!
      integer(i15) :: p = 0, q = 0, pq = 0
!
      integer(i15), intent(in) :: dim_1,dim_2
      real(dp), dimension(dim_1, dim_2), intent(in) :: vec
!
      do q = 1, dim_2
         do p = 1, dim_1
!
            write(unit_output,*) p, q, vec(p,q)
!
         enddo
      enddo
!
   end subroutine vec_print
!
!
end module utils
