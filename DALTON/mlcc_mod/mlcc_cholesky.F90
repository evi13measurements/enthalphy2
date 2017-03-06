module mlcc_cholesky
!
    use mlcc_data
    use mlcc_utilities
    use mlcc_types
!
contains
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
   subroutine read_cholesky_ai(L_ai_J)
!
!     Purpose: Read Cholesky vectors L_ia^J from file and place them 
!              in the incoming vector  
!
      use mlcc_workspace
      implicit none
!
      double precision L_ai_J(n_ov,n_J)
!
      integer :: lucho_ia 
      integer :: i,a,ia,ai,j
      real(dp),dimension(:,:),pointer :: L_ia_J
!
!     Allocation
!
      call allocator(L_ia_J,n_ov,n_J)
      L_ia_J=zero
!     
!     IO
!
      call read_cholesky_ia(L_ia_J)
!
!     Reorder
!      
      do i = 1,n_occ
         do a = 1,n_vir
!
!           Needed indices
!
            ai=index_two(a,i,n_vir)
            ia=index_two(i,a,n_occ)
!
            do j=1,n_J
               L_ai_J(ai,j)=L_ia_J(ia,j)
            enddo
         enddo
      enddo
!
!     Deallocation
!
      call deallocator(L_ia_J,n_ov,n_J)

   end subroutine read_cholesky_ai
!
   subroutine read_cholesky_ij(L_ij_J)
!
!     Purpose: Read Cholesky vectors L_ij^J from file and place them 
!              in the incoming vector  
!
      implicit none
!
      double precision L_ij_J(n_oo,n_J)
!
      integer :: lucho_ij
      integer :: i,j,idummy
!
      lucho_ij = -1
      call gpopen(lucho_ij,'CHOLESKY_IJ','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
      rewind(lucho_ij)
!
      do j = 1,n_J
         read(lucho_ij) (L_ij_J(i,j), i=1,n_oo)
      enddo
!
      call gpclose(lucho_ij,'KEEP')    
!      
   end subroutine read_cholesky_ij
!
   subroutine read_cholesky_ab(L_ab_J,a_start,a_end,ab_dim)
!
!  Purpose: Read Cholesky vectors L_ab^J from file and place them 
!           in the incoming vector. If b_start .ne. 1 and b_end .eq. n_vir
!           we batch
! 
!           b_start - first element to be read
!           b_end   - last element to be read
!
!           ab_dim  - dimension over batching variables 
   implicit none
!
   integer :: lucho_ab,ab_dim
   integer :: a,b,j,idummy,i
   integer :: a_start,a_end
   real(dp),dimension(ab_dim,n_J) :: L_ab_J
   real(dp) :: dummy
   integer  :: batch_length 
!
   batch_length = a_end-a_start+1
!
   lucho_ab = -1
   call gpopen(lucho_ab,'CHOLESKY_AB','UNKNOWN','SEQUENTIAL','UNFORMATTED',idummy,.false.)
   rewind(lucho_ab)
!
   if (a_start .ne. 1) then
!
!     Calculate index of last element to throw away
!
      idummy=index_two(n_vir,a_start-1,n_vir)
!
!     Read from a_start
!
      do j = 1,n_J
        read(lucho_ab)(dummy,i=1,idummy),(L_ab_J(a,j),a=1,ab_dim)
      enddo
   else
!
!     Read from start
!
      do j = 1,n_J
        read(lucho_ab)(L_ab_J(a,j),a=1,ab_dim)
      enddo
   endif
   !
   call gpclose(lucho_ab,'KEEP')    
!   
   end subroutine read_cholesky_ab
!
end module mlcc_cholesky