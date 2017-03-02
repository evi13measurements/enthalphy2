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
!
   subroutine n_one_batch(required,available,max_batch_length,n_batch,batch_dimension)
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
      integer, intent(in)           :: required, available, batch_dimension
      integer                       :: max_batch_length,n_batch
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
      n_batch=required/max_batch_length
!
!  Test for rest
!
      if (n_batch*max_batch_length .lt. required) then
         n_batch = n_batch+1
      endif
!
   end subroutine n_one_batch
!
   subroutine one_batch_limits(first,last,batch_number,max_batch_length,batch_dimension)
!
!     Purpose: Find batch limits (first and last) 
!
!        batch_number: the current batch (1,2,...,n_batch)
!        max_batch_length: the length of each batch (except the last, which may be a rest, see n_one_batch routine)
!        batch_dimension: the dimensionality of the batching variable (e.g., n_vir for a virtual index)
!
      implicit none 
!
      integer :: first,last
      integer, intent(in) :: batch_number,max_batch_length,batch_dimension
!
      first = 1 + (batch_number-1)*max_batch_length
      last  = min(max_batch_length+(batch_number-1)*max_batch_length,batch_dimension)
!
   end subroutine one_batch_limits
!
   subroutine read_cholesky_ia(L_ia_J)
!
!     Purpose: Read Cholesky vectors L_ia^J from file and place them 
!              in the incoming vector  
!
      implicit none
!
      double precision L_ia_J(n_ov,n_J)
!
      integer :: lucho_ia 
      integer :: i,j,idummy
!
      lucho_ia = -1
      call gpopen(lucho_ia,'CHOLESKY_IA','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
      rewind(lucho_ia)
!
      do j=1,n_J
         read(lucho_ia) (L_ia_J(i,j), i=1,n_ov)
      enddo
!
      call gpclose(lucho_ia,'KEEP')      
!
   end subroutine read_cholesky_ia
!
   subroutine read_cholesky_ij(L_ij_J)
!
!     Purpose: Read Cholesky vectors L_ij^J from file and place them 
!              in the incoming vector  
!
      implicit none
!
      double precision L_ij_J(n_oo_packed,n_J)
!
      integer :: lucho_ij
      integer :: i,j,idummy
!
      lucho_ij = -1
      call gpopen(lucho_ij,'CHOLESKY_IJ','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
      rewind(lucho_ij)
!
      do j = 1,n_J
         read(lucho_ij) (L_ij_J(i,j), i=1,n_oo_packed)
      enddo
!
      call gpclose(lucho_ij,'KEEP')    
!      
   end subroutine read_cholesky_ij
   subroutine read_cholesky_ab(L_ab_J)
!
!  Purpose: Read Cholesky vectors L_ab^J from file and place them 
!           in the incoming vector  
!
   implicit none
!
   double precision L_ab_J(n_vv,n_J)
!
   integer :: lucho_ab
   integer :: a,j,idummy
!
   lucho_ab = -1
   call gpopen(lucho_ab,'CHOLESKY_AB','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
   rewind(lucho_ab)
!
   do j = 1,n_J
      read(lucho_ab) (L_ab_J(a,j), a=1,n_vv)
   enddo
!
   call gpclose(lucho_ab,'KEEP')    
!   
   end subroutine read_cholesky_ab
!
   subroutine vec_print(vec,dim_1,dim_2)
!
!     Purpose: prints a vector with a compound index (p q) of dimension (dim_1 x dim_2)
!        (for debugging, remove or replace later on)
!
      implicit none
!
      integer p,q,pq
      integer dim_1,dim_2
      double precision vec(dim_1,dim_2)
!
      do q = 1,dim_2
         do p = 1,dim_1
            write(luprint,*) p,q,vec(p,q)
         enddo
      enddo
!
   end subroutine vec_print
end module mlcc_utilities