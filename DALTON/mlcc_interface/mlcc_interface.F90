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
!  MLCC Cholesky vector reader, and transformator
!  Authors:  Eirik Kjønstad and Sarai Folkestad, January 2017
!
!  Purpose: Read Cholesky vectors in AO basis, transform to MO basis and write to files CHOLESKY_IJ,
!           CHOLESKY_IA, CHOLESKY_AI and CHOLESKY_AB according to PQ-type.
!
   use mlcc_data
   use mlcc_workspace
!
   implicit none
!
   integer       :: ludiag,idum,i,lucho,lucho_ij,lusec,dummy
   integer       :: n_cho,n_J
   real(dp), dimension(:,:), pointer      :: cho_diag
   real(dp), dimension(:,:), pointer      :: cho_ao
   integer       :: lencho
!
! Read number (n_J) and length (n_cho) of Cholesky vectors
!
  lusec = -1
  call gpopen(lusec,'CHOLESKY.RST','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
  rewind(lusec)
  read(lusec) dummy,dummy,dummy,dummy,n_cho,n_J,dummy
  call gpclose(lusec,'KEEP')
!
!
  write(ml_lupri,*)
  write(ml_lupri,*)n_cho, n_J
  write(ml_lupri,*)
!
!  Allocation for diagonal elements
   call allocator(cho_diag,n_cho,1)
   cho_diag=zero
!
!  Read diagonal
!
   ludiag = -1
!
   call gpopen(ludiag,'CHODIAG','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
   rewind(ludiag)
   read(ludiag) (cho_diag(i,1),i=1,n_cho)
   call gpclose(ludiag,'KEEP')
!
   do i = 1,n_cho
      write(ml_lupri,*)i, cho_diag(i,1) 
   enddo
!  
! Skip reduction. Need n_reduced for reading? 
!
! Read Cholesky AO
!
   lucho = -2
   call gpopen(lucho,'CHOLES_1_0','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
   rewind(lucho)

   call gpclose(lucho,'KEEP')
! ij file io
   lucho_ij = -3
   call gpopen(lucho_ij,'CHOLESKY_IJ','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
   rewind(lucho_ij)
   call gpclose(lucho_ij,'KEEP')
   
end subroutine mlcc_get_cholesky

subroutine hf_reader
   use mlcc_data
   use mlcc_workspace
!
!  Hartree-Fock reader routine
!  Authors Henrik Koch, Rolf H. Myhre, Sarai Folkestad, Eirik Kjønstad
!  January 2017
!
!  Purpose: read in data from LUSIFC file
!
      implicit none
!
!
      integer  :: lusifc = -1 
      integer  :: idummy = 0
      integer  :: n_symmetries, n_basis_sym, n_orbitals_sym
      integer  :: i
!      
!
!     
!     Open Sirius Fock file
!     ---------------------
!
      call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',idummy,'.FALSE.')
      rewind(lusifc)
!
!
!     Read in various stuff from Sirius Fock file. Things depending on symmetry is mostly
!     discarded at the end of the subroutine as we do not use symmetry. Information in 
!     file should be in Cholesky orbital format. If Cholesky orbitals has been generated,
!     SIRIFC will contain the data in the Cholesky basis
!
      call mollab('TRCCINT ',lusifc,ml_lupri)
!      
      read(lusifc) n_symmetries, n_orbitals, n_basis, n_lambda, n_occ, &
      &            n_orbitals_sym, n_basis_sym, nuclear_potential, scf_energy
!
!      
      if (n_symmetries /= 1) call quit('error in mlcc_mod_init: not implemented with symmetry')
!      
!      
!     Calculate number of virtuals and amplitudes
!
      n_vir          = n_orbitals - n_occ
      n_t1am         = n_vir*n_occ
      n_t2am         = n_t1am*n_t1am
      n_t2am_pack    = n_t1am*(n_t1am+1)/2
      write(ml_lupri,*) n_t1am,n_t2am,n_t2am_pack
!      
!     Allocate space for Fock diagonal and coefficients.
!
      call allocator(fock_diagonal,n_orbitals,1)
!      
      call allocator(orb_coefficients,n_lambda,1)
!
!     Read in Fock diagonal and coefficients
!
      read(lusifc) (fock_diagonal(i,1),i=1,n_orbitals)
      read(lusifc) (orb_coefficients(i,1),i=1,n_lambda)
!
!     Done with file
!
      call gpclose(lusifc,'KEEP')
!
   end subroutine hf_reader

