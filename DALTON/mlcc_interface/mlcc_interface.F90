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
   use mlcc_input_data
   use mlcc_workspace
!
!
   implicit none
!
   integer, intent(in)              :: lunit !input unit
   integer, intent(in)              :: lupri !general output unit
!
   character(len=7), intent(inout)  :: word
   character(len=24)                :: buffer
!
   integer                          :: r_off, length, last, actual
   integer                          :: i_count,i, r_len
!
   r_off=1
!
   if (word(1:7) .eq. '*MLCC  ') then
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
         select case (trim(buffer))
!        
            case('.MLACTIVE')
!        
               mlcc_active=.true.
!        
            case('.PRINT')
!        
               read (lunit,*) print_mlcc
!        
            case default
               write(lupri,*) 'Keyword ', trim(buffer), ' not recognized in mlcc_input'
               stop
!        
         end select
!
!
      end do
!
   else 
      write(lupri,*) 'mlcc_input called without *MLCC'
      call flshfo(lupri)
      call quit('mlcc_input, something is wrong')
   end if
!
   if(print_mlcc .ge. 3) then
!
      write(lupri,*)
      write(lupri,*) 'Output from mlcc3_input'
      write(lupri,*) 'print_mlcc', print_mlcc
      write(lupri,*) 'mlcc_active', mlcc_active
      write(lupri,*)
!
   end if
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
!
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
   use mlcc_utilities
!
   implicit none
!
   integer       :: ludiag,lucho,lucho_ij,lucho_ia,lucho_ab,lumlch
   integer       :: i,j,k,idum
   integer       :: n_twoel_diag
   integer, dimension(:,:), pointer       :: index_reduced
   real(dp), dimension(:,:), pointer      :: cho_diag
   real(dp), dimension(:,:), pointer      :: cho_ao
   integer       :: lencho
   real(dp), dimension(:,:), pointer      :: cholesky_ao,cholesky_ao_sq
   real(dp), dimension(:,:), pointer      :: cholesky_mo_sq,cholesky_mo,X
!
   n_twoel_diag = packed_size(n_basis)
!
! Open MLCC_CHOLESKY - See mlcc_write_cholesky.F
!  
   lumlch = -1
   call gpopen(lumlch,'MLCC_CHOLESKY','OLD','SEQUENTIAL','FORMATTED',idum,.false.)
   rewind(lumlch)
!
! Read number (n_J) and length (n_reduced) of Cholesky vectors
!
   read(lumlch,*)n_reduced,n_J
!
! Read reduced index array
!
   call allocator_int(index_reduced,n_reduced,1)
   index_reduced=0
!
   read(lumlch,*)(index_reduced(i,1),i=1,n_reduced)
!
!========================================
!
! Preparation of files for cholesky in MO basis:
! ij-type:
  lucho_ij = -1
  call gpopen(lucho_ij,'CHOLESKY_IJ','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
  rewind(lucho_ij)
! ai-type:
  lucho_ia = -1
  call gpopen(lucho_ia,'CHOLESKY_IA','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
  rewind(lucho_ia)
! ab-type:
  lucho_ab = -1
  call gpopen(lucho_ab,'CHOLESKY_AB','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
  rewind(lucho_ab)
!
!
!
! Allocation
! Cholesky AO array, intermediate X and cholesky MO array
!
   call allocator(cholesky_ao,n_twoel_diag,1)
   call allocator(cholesky_ao_sq,n_basis,n_basis)
   call allocator(X,n_basis,n_orbitals)
   call allocator(cholesky_mo_sq,n_orbitals,n_orbitals)
!
   cholesky_ao = zero
   cholesky_ao_sq = zero
   cholesky_mo_sq = zero
   X = zero
   do j=1,n_J
!
! Read Cholesky AO:
!
      read(lumlch,*) (cholesky_ao(i,1),i=1,n_twoel_diag)
! Square up of cholesky_ao:
!
      call squareup(cholesky_ao,cholesky_ao_sq,n_basis)
!
! Transform to MO
!  
      call dgemm('N','N',n_basis,n_orbitals,n_basis,one,cholesky_ao_sq,n_basis,orb_coefficients,n_basis,zero,X,n_basis)
      call dgemm('T','N',n_orbitals,n_orbitals,n_basis,one,orb_coefficients,n_basis,X,n_basis,zero,cholesky_mo_sq,n_orbitals)
!
! Packing MO Cholesky and save to file
!
      write(lucho_ij)((cholesky_mo_sq(i,k),k=1,i),i=1,n_occ)
      write(lucho_ia)((cholesky_mo_sq(i,k),k=n_occ+1,n_orbitals),i=1,n_occ)
      write(lucho_ab)((cholesky_mo_sq(i,k),k=n_occ+1,i),i=n_occ+1,n_orbitals)
!
  enddo
!
!  Close all files
!
   call gpclose(lumlch,'KEEP')
   call gpclose(lucho_ij,'KEEP')
   call gpclose(lucho_ia,'KEEP')
   call gpclose(lucho_ab,'KEEP')
!
!  Deallocation
!
   call deallocator(cholesky_ao,n_twoel_diag,1)
   call deallocator(cholesky_ao,n_basis,n_basis)
   call deallocator_int(index_reduced,n_reduced,1)
   call deallocator(cholesky_mo_sq,n_orbitals,n_orbitals)
   call deallocator(X,n_basis,n_orbitals)
!
!================================
end subroutine mlcc_get_cholesky
!
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
      use mlcc_data
      use mlcc_workspace
!
      implicit none
!
!
      integer  :: lusifc = -1 
      integer  :: idummy = 0
      integer  :: n_symmetries, n_basis_sym, n_orbitals_sym
      integer  :: i,j
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
!      
!     Allocate space for Fock diagonal and coefficients.
!
      call allocator(fock_diagonal,n_orbitals,1)
      fock_diagonal = zero
!      
      call allocator(orb_coefficients,n_lambda,1)
      orb_coefficients = zero

!
!     Read in Fock diagonal and coefficients
!
      read(lusifc) (fock_diagonal(i,1),i=1,n_orbitals)
      write(*,*) (fock_diagonal(i,1),i=1,n_orbitals)
      read(lusifc) (orb_coefficients(i,1),i=1,n_lambda) 
!
!     Done with file
!
      call gpclose(lusifc,'KEEP')
!

   end subroutine hf_reader

