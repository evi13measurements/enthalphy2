!Interface routines to new MLCC code
!Requires a call to MLCC_input, mlcc_driver, mlcc_HF_reader and mlcc_ao_reader
!Dalton input is a mess, where is the best way to call the routine SIR_NEWINP?
!
!
subroutine mlcc_interface_input(lunit,word,lupri)
!
!  Read in input from  DALTON.INP, starting with *MLCC until next * and store in input string
!  Pass on input to mlcc_input in mlcc_module
!
!  Rolf H. Myhre 
!  December 2016
!
   use mlcc_drive, only: mlcc_input
!
!
   implicit none
!
   integer, intent(in)              :: lunit !input unit
   integer, intent(in)              :: lupri !general output unit
!
   character(len=7), intent(inout)  :: word
!
   character(len=254)               :: buffer
   character(len=254)               :: in_string
!
   integer                          :: r_off, length, last, actual
   integer                          :: i_count,i, r_len
!
   r_off=1
!
   if (word(1:7) .eq. '*MLCC  ') then
!
      in_string(1:7) = word
      i_count = i_count + 1
      r_off = r_off+7
!
      do
!
!        Read next keyword      
         read (lunit,'(a)') buffer
!
!        Skip comments
         do while (buffer(1:1) .eq. '!' .or. buffer(1:1) .eq. '#' )
            read (lunit,'(a)') buffer
         end do
!            
!        Stop when hitting the next *
         if (buffer(1:1) .eq.'*') then
            word = buffer(1:7)
            exit
         end if
!
!
         r_len = len(trim(buffer))
!
         write(in_string(r_off:),'(a)') trim(buffer)
         r_off = r_off + r_len + 1
         i_count = i_count + 1
!
      end do
!
!     Fortran adds whitespace to integers.
!     Strip off excess whitespaces and count number of records
!
      i_count = 1
      actual = 1
      last = 1
      length = len(trim(in_string))
!
      do while (actual .lt. length)
!
         if (in_string(last:last+1) .eq. "  ") then
            in_string(last+1:) = in_string(last+2:)
            actual = actual+1
         else if (in_string(last:last) .eq. " ") then
            i_count = i_count + 1
            last=last+1
            actual=actual+1
         else
            last=last+1
            actual=actual+1
         end if
!
      end do
!
      write(lupri,*) 'Input string: ', trim(in_string)
      call flush(lupri)
!
   else 
      write(lupri,*) 'mlcc_input called without *MLCC'
      call flshfo(lupri)
      call quit('mlcc_input, something is wrong')
   end if
!
   call mlcc_input(in_string,i_count,lupri)
!
   write(lupri,*)
   write(lupri,*) 'left mlcc_input'
   write(lupri,*)
!
end subroutine mlcc_interface_input
!
!
subroutine mlcc_interface_drv(work,lwork,lupri)
!
!  Interface routine for the MLCC driver. 
!
!  Rolf H. Myhre 
!  December 2016
!
   use mlcc_drive, only: mlcc_drv
!
   implicit none
!
   integer, intent(in)                    :: lupri !general output unit
   integer, intent(in)                    :: lwork !free space in work
   real*8, intent(in), dimension(lwork)   :: work !work statis array
!
   call mlcc_drv(work,lwork,lupri)
!
end subroutine mlcc_interface_drv
!
subroutine mlcc_iajb(vec)
!
!  MLCC iajb MO integral calculator
!  Authors Sarai Folkestad, Eirik Kjønstad, January 2017
!
!  Purpose: read in atomic orbitals and calculate iajb MO integrals
!
   use mlcc_data
   use mlcc_types
!
   implicit none
!
   real(dp), dimension(n_t2am_pack,1) :: vec
!
   vec = zero
!
end subroutine mlcc_iajb

subroutine mlcc_get_cholesky()
!
   use mlcc_data
   use mlcc_workspace
!
   implicit none
!
   integer       :: ludiag,idum,i,lucho,lucho_ij,lusec,dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
   integer       :: n_cho_diag,n_cho
   real(dp), dimension(:,:), pointer      :: cho_diag
   real(dp), dimension(:,:), pointer      :: cho_ao
!
   call get_n_cho_diag(n_cho_diag,n_orbitals)
!
!  Allocation for diagonal elements
   call allocator(cho_diag,n_cho_diag,1)
   cho_diag=zero
!
!  Read diagonal
!
   ludiag = -1
!
   call gpopen(ludiag,'CHODIAG','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
   rewind(ludiag)
   read(ludiag) (cho_diag(i,1),i=1,n_cho_diag)
   call gpclose(ludiag,'KEEP')
!
   do i = 1,n_cho_diag
      write(ml_lupri,*)i, cho_diag(i,1) 
   enddo
!  
! Skip reduction. Need n_reduced for reading? 
!
! Read number of cholesky vectors
!  lusec = -1
!  call gpopen(lusec,'CHOLESKY.SEC','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
!  read(lusec) dummy1,dummy2,dummy3,dummy4,n_cho,dummy5,dummy6
!  call gpclose(lusec,'KEEP')
!
!  write(ml_lupri,*)
!  write(ml_lupri,*)n_cho THIS PRINTS OUT 55! SO IT IS = n_cho_diag NOT J. dummy5 looks like reduced dim (=54)
!  write(ml_lupri,*)
! Read Cholesky AO
!
   lucho = -2
   call gpopen(lucho,'CHOLES_1_0','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
   call gpclose(lucho,'KEEP')
! ij file io
   lucho_ij = -3
   call gpopen(lucho_ij,'CHOLESKY_IJ','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
   call gpclose(lucho_ij,'KEEP')
   
end subroutine mlcc_get_cholesky

subroutine get_n_cho_diag(n_cho_diag,n_orbitals)

   implicit none
!
   integer, intent(in)    ::  n_orbitals
   integer                ::  n_cho_diag
   integer                ::  i,j,counter
!
   counter=0
   do i=1,n_orbitals
      do j=1,n_orbitals
         if (i .ge. j) then
            counter = counter+1
         endif
      enddo
   enddo
!
   n_cho_diag=counter  
!
end subroutine get_n_cho_diag
