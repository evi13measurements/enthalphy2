module cholesky_integrals_class
!
   type cholesky_integrals
      integer           :: n_J  = 0
      integer           :: n_ao = 0
      integer           :: n_mo = 0
      character(len=10) :: integral_program
   contains
!	
      procedure :: init => init_cholesky_integrals
!
   end type cholesky_integrals
!
contains
!
   subroutine init_cholesky_integrals(cholesky)
!
      use dalton_interface
!
      implicit none
!
      class(cholesky_integrals)   :: cholesky
!
      cholesky % integral_program = 'DALTON    '
!
      if (cholesky % integral_program .eq. 'DALTON    ') then
!
         call dalton_interface_driver(cholesky)
!
      endif
!
   end subroutine init_cholesky_integrals
!
   subroutine dalton_interface_driver(cholesky)
!
!  Dalton interface driver
!
      use input_output
!
      implicit none
!
      class(cholesky_integrals) :: cholesky
      integer                   :: unit_cholesky_ao = -1
!
!     Getting identifier for AO integrals file    
!
      call generate_unit_identifier(unit_cholesky_ao)
!
!   n_twoel_diag = packed_size(n_basis)
!
!     Open MLCC_CHOLESKY - See mlcc_write_cholesky.F
!  
!      call gpopen(unit_cholesky_ao,'MLCC_CHOLESKY','OLD','SEQUENTIAL','FORMATTED',idum,.false.)
      open(unit=unit_cholesky_ao,file='MLCC_CHOLESKY',status='old',form='formatted')
      rewind(unit_cholesky_ao)
!
!     Read number (n_J) and length (n_reduced) of Cholesky vectors
!
      read(unit_cholesky_ao,*) n_reduced, cholesky % n_J
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
   end subroutine dalton_interface_driver
!
end module cholesky_integrals_class
