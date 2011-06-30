!> @file
!> Contains module that reads the dalton input and initializes lsitem
MODULE DALTONINFO
use files
use precision
use typedef
use integraldriver
use READMOLEFILE
use BuildBasisSet
use lstiming
use fragment
use molecule_module

CONTAINS 
!> \brief initiate the lsitem structure which contain all info about integralevaluation schemes, molecule and basisset.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ls_init(ls,lupri,luerr,nbast,integral,dodft,doprint)
implicit none
!> contains everything needed for integralevaluation (incl. molecule,basis)
TYPE(LSITEM) :: ls
!> the logical unit number for the output file
Integer,intent(in)  :: lupri
!> the logical unit number for the error file
Integer,intent(in)  :: luerr
!> the number of basisfunctions
integer,intent(out) :: nbast
!> a collection of logicals read from the DALTON.INP file
type(integralconfig) :: integral 
!> Should we print stuff to output file
Logical      :: doprint
!> is this a dft calculation or a HF calculation
Logical      :: dodft

ls%lupri = lupri
ls%luerr = luerr
CALL dalton_init(ls%input,lupri,luerr,nbast,integral,ls%fopen,dodft,doprint)
CALL II_init_setting(ls%setting)
CALL II_set_default_setting(ls%setting,ls%input)

END SUBROUTINE ls_init

!> \brief frees the lsitem structure.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ls_free(ls)
implicit none
!> contains everything needed for integralevaluation (incl. molecule,basis)
TYPE(LSITEM) :: ls
CALL dalton_finalize(ls%input,ls%lupri,ls%luerr,ls%fopen)
CALL II_free_setting(ls%setting)
END SUBROUTINE ls_free

!> \brief build the daltoninput structure associated with the lsitem 
!> \author T. Kjaergaard
!> \date 2010
!>
!> build the daltoninput structure associated with the lsitem, which means 
!> initialize the IOitem (filenames) read MOLECULE.INP file and build molecule
!> and finally build the basis. These informations should not be changed and
!> the name ls%input indicate that it is input information, while
!> ls%setting is the current info about molecule,basis, etc.
!>
SUBROUTINE dalton_init(intinp,LUPRI,LUERR,nbast,integral,fopen,dodft,doprint)
use infpar_module
implicit none
!> contains info about input from DALTON.INP(only integral info) and MOLECULE.INP
TYPE(DALTONINPUT)    :: intinp
!> the logical unit number for the output file
integer,intent(in)   :: LUPRI
!> the logical unit number for the error file
integer,intent(in)   :: LUERR
!> number of basisfunctions
integer, intent(out) :: nbast
!> information about the integral evaluation info read from DALTON.INP
type(integralconfig)      :: integral
!> nolonger used
logical              :: fopen
!> is it a DFT run
logical              :: dodft
!> should we print to output file
LOGICAL              :: doprint
!
TYPE(BASISSETLIBRARYITEM) :: LIBRARY
integer              :: LUCME,IDUMMY,I
CHARACTER(len=80)    :: BASISSETNAME
CHARACTER(len=9)     :: BASISLABEL
real(realk)          :: Tim1,Tim2,TIME_BUILD_BASIS,TIME_READ_DALTON
real(realk)          :: TIME_AUXBUILD_BASIS,TIME_READ_MOL,TSTART
Integer              :: numAtoms, numNodes
logical              :: lopen
LUCME=-1
!CAREFUL


!IF(LINSCARUN)THEN
!  INQUIRE(FILE='LSDALTON.OUT',OPENED=lopen)
!  fopen = (.NOT. lopen).AND.(LUPRI.EQ.-1)
!  IF (fopen.AND.(LUERR.NE.-1)) CALL LSQUIT('dalton_init LUERR and LUPRI not both -1')
!  IF (fopen) THEN
!     CALL LSOPEN(LUPRI,'LSDALTON.OUT','NEW','FORMATTED')
!     CALL LSOPEN(LUERR,'LSDALTON.ERR','UNKNOWN','FORMATTED')
!  ENDIF
!ELSE
!  INQUIRE(FILE='newLSDALTON.OUT',OPENED=lopen)
!  fopen = (.NOT. lopen).AND.(LUPRI.EQ.-1)
!  IF (fopen.AND.(LUERR.NE.-1)) CALL LSQUIT('dalton_init LUERR and LUPRI not both -1')
!  IF (fopen) THEN
!     CALL LSOPEN(LUPRI,'newLSDALTON.OUT','NEW','FORMATTED')
!     CALL LSOPEN(LUERR,'newLSDALTON.ERR','UNKNOWN','FORMATTED')
!  ENDIF
!ENDIF
intinp%nfock = 0 !number of fock matrices calculated

NULLIFY(intinp%MOLECULE)
NULLIFY(intinp%BASIS)
ALLOCATE(intinp%MOLECULE)
ALLOCATE(intinp%BASIS)

IF (doprint) THEN
   !CALL PRINT_INTRO(LUPRI)
   CALL LS_TSTAMP(' ',LUPRI)
   CALL PRINT_MOL_FILE(LUPRI)
   CALL PRINT_DALTON_FILE(LUPRI)
ENDIF

!Initialize IO
CALL io_init(intinp%IO)

!*************************************************
!*
!*  READ THE LINSCA INTEGRAL SECTION OF THE DALTON
!*  INPUT-FILE 
!*
!*************************************************
CALL LSTIMER('START',TIM1,TIM2,LUPRI)

intinp%DALTON = integral

IF(intinp%DALTON%TRILEVEL) intinp%DALTON%NOSEGMENT = .true.
IF(intinp%DALTON%TRILEVEL) intinp%DALTON%NOFAMILY = .true.

IF(intinp%DALTON%TIMINGS) CALL LSTIMER('READ DALTONFILE',TIM1,TIM2,LUPRI)
!*************************************************
!*
!*  READ THE MOLECULE.INP FILE AND BUILD THE MOLECULE 
!*  STRUCTURE
!*
!*************************************************

CALL READ_MOL_FILE_AND_BUILD_MOLECULE(LUPRI,intinp%MOLECULE,LIBRARY,doprint,intinp%dalton%molprint,&
     &intinp%dalton%DoSpherical,intinp%dalton%IntegralThreshold,intinp%dalton%auxbasis)

integral%nelectrons = intinp%MOLECULE%nelectrons 
integral%molcharge = INT(intinp%MOLECULE%charge)
numAtoms = intinp%MOLECULE%nAtoms
#ifdef VAR_LSMPI
  intinp%numNodes     = infpar%nodtot
  intinp%numFragments = min(numAtoms,intinp%numNodes)
#else
  intinp%numNodes     = 1
  intinp%numFragments = 1
#endif
IF (doprint) THEN
  IF (intinp%numFragments.GT.1) WRITE(LUPRI,*) 'Integrals calculated using ',intinp%numFragments, ' fragments'
ENDIF

IF(intinp%DALTON%TIMINGS) CALL LSTIMER('READ MOLFILE',TIM1,TIM2,LUPRI)

!*************************************************
!*
!*  CALL LINSCF WITH DEFAULT BASIS AS DEFINED IN INPUT
!*
!*************************************************

IF (doprint) THEN
  WRITE(LUPRI,'(2X,A)')' '
  WRITE(LUPRI,'(2X,A)')'CALLING BUILD BASIS WITH DEFAULT REGULAR BASIS'
  WRITE(LUPRI,'(2X,A)')' '
ENDIF

BASISLABEL='REGULAR  '
CALL Build_basis(LUPRI,intinp%DALTON%BASPRINT,&
     &intinp%MOLECULE,intinp%BASIS%REGULAR,LIBRARY,&
     &BASISLABEL,intinp%DALTON%UNCONT,intinp%DALTON%NOSEGMENT,doprint)

nbast = getNbasis('Regular','Contracted',intinp%MOLECULE,LUPRI)

IF(intinp%DALTON%TIMINGS) CALL LSTIMER('BUILD BASIS',TIM1,TIM2,LUPRI)

IF (doprint) THEN
  CALL PRINT_BASISSETLIBRARY(LUPRI,LIBRARY)
ENDIF

IF(intinp%DALTON%AUXBASIS)THEN

   IF (doprint) THEN
      WRITE(LUPRI,'(2X,A)')' '
      WRITE(LUPRI,'(2X,A)')'CALLING BUILD BASIS WITH DEFAULT AUXILIARY BASIS'
      WRITE(LUPRI,'(2X,A)')' '
   ENDIF
   BASISLABEL='AUXILIARY'
   CALL Build_basis(LUPRI,intinp%DALTON%BASPRINT,&
     &intinp%MOLECULE,intinp%BASIS%AUXILIARY,LIBRARY,&
     &BASISLABEL,intinp%DALTON%UNCONT,intinp%DALTON%NOSEGMENT,doprint)
   IF(intinp%DALTON%TIMINGS) CALL LSTIMER('BUILD AUXBASIS',TIM1,TIM2,LUPRI)
   IF (doprint) THEN
      CALL PRINT_MOLECULE_AND_BASIS(LUPRI,intinp%MOLECULE,intinp%BASIS%AUXILIARY)
   ENDIF
   
ELSE
intinp%BASIS%AUXILIARY%natomtypes = 0
ENDIF

intinp%BASIS%HUCKEL%natomtypes = 0

IF (doprint) THEN
   CALL PRINT_MOLECULEINFO(LUPRI,intinp%MOLECULE,intinp%BASIS)
   CALL PRINT_MOLECULE_AND_BASIS(LUPRI,intinp%MOLECULE,intinp%BASIS%REGULAR)
ENDIF

intinp%DO_DFT = dodft
!CALL TSTAMP(' ',LUPRI)
!CALL LSCLOSE(LUCME,'KEEP')
!CALL LSCLOSE(LUMOLDEN,'KEEP')

END SUBROUTINE dalton_init

!> \brief finalize the daltoninput structure
!> \author T. Kjaergaard
!> \date 2010
!>
!> free molecule info
!> free basis
!> free io
!>
SUBROUTINE dalton_finalize(intinp,LUPRI,LUERR,FOPEN)
implicit none
!> contains info about input from DALTON.INP(only integral info) and MOLECULE.INP
TYPE(DALTONINPUT)    :: intinp
!> the logical unit number for the output file
integer   :: LUPRI
!> the logical unit number for the error file
integer   :: LUERR
!> should the files LUPRI and LUERR be closed
LOGICAL :: FOPEN

IF (FOPEN) THEN
  CALL LSCLOSE(LUERR,'KEEP')
  CALL LSCLOSE(LUPRI,'KEEP')
ENDIF

call free_Moleculeinfo(intinp%MOLECULE)
IF(intinp%BASIS%REGULAR%natomtypes .NE.0)THEN
   call free_basissetinfo(intinp%BASIS%REGULAR)
ENDIF
IF(intinp%BASIS%AUXILIARY%natomtypes .NE.0)THEN
   call free_basissetinfo(intinp%BASIS%AUXILIARY)
ENDIF
IF(intinp%BASIS%HUCKEL%natomtypes .NE.0)THEN
   call free_basissetinfo(intinp%BASIS%HUCKEL)
ENDIF
DEALLOCATE(intinp%MOLECULE)
DEALLOCATE(intinp%BASIS)
call io_free(intinp%io)

END SUBROUTINE dalton_finalize

!> \brief free the daltoninput structure
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE dalton_free(intinp)
implicit none
!> contains info about input from DALTON.INP(only integral info) and MOLECULE.INP
TYPE(DALTONINPUT) :: intinp

call free_Moleculeinfo(intinp%MOLECULE)
call free_basissetinfo(intinp%BASIS%REGULAR)
IF(intinp%DALTON%AUXBASIS)THEN
   call free_basissetinfo(intinp%BASIS%AUXILIARY)
ENDIF
DEALLOCATE(intinp%MOLECULE)
DEALLOCATE(intinp%BASIS)
NULLIFY(intinp%MOLECULE)
NULLIFY(intinp%BASIS)
call io_free(intinp%io)

END SUBROUTINE dalton_free
 
!> \brief print dalton file
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE PRINT_DALTON_FILE(LUPRI)
!> the logical unit number for the output file
implicit none
INTEGER            :: LUCMD !Logical unit number for the daltoninput
INTEGER            :: IDUMMY,LUPRI
character(len=40)  :: WORD
character(len=2)   :: PROMPT
LOGICAL            :: DONE_LINSCA,file_exist

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,*) '-----------------------------------------'
WRITE(LUPRI,*) '         PRINTING THE DALTON.INP FILE '
WRITE(LUPRI,*) '-----------------------------------------'
WRITE(LUPRI,*) '                     '

INQUIRE(file='DALTON.INP',EXIST=file_exist) 
IF(file_exist)THEN
  LUCMD=-1
  CALL LSOPEN(LUCMD,'DALTON.INP','OLD','FORMATTED')
ELSE
  CALL LSQUIT('DALTON.INP does not exist',lupri)
ENDIF
rewind(LUCMD)
DO
   READ (LUCMD, '(A40)') WORD
   PROMPT = WORD(1:2)
   IF(PROMPT(1:1) .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
   SELECT CASE(WORD) 
   CASE ('*END OF INPUT')
      WRITE(LUPRI,'(2X,A40)') WORD   
      EXIT
   CASE DEFAULT
      WRITE(LUPRI,'(2X,A40)') WORD
   END SELECT
ENDDO
WRITE(LUPRI,*)' '

CALL LSCLOSE(LUCMD,'KEEP')

END SUBROUTINE PRINT_DALTON_FILE

!> \brief print molecule file
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE PRINT_MOL_FILE(LUPRI)
implicit none
!> the logical unit number for the output file
integer            :: lupri
!
integer            :: LUINFO,IDUMMY,I,J,IPOS,IPOS2
logical            :: file_exist,Angstrom,Symmetry,dopbc
logical            :: ATOMBASIS,BASIS,ATOMDF,AUXBASIS
CHARACTER(len=80)  :: WORD
integer            :: Atomtypes,natoms,Molecularcharge,ios
integer,ALLOCATABLE:: BasisSetCharge(:) !Charge of each type 
real(realk),PARAMETER :: AngstromConvert = 0.5291772083D0
real(realk)        :: IntegralThreshold,Q
logical            :: DoOwn,DoCartesian
CHARACTER(len=1) :: KASYM(3,3),CRT
CHARACTER(len=80),ALLOCATABLE  :: ATOMBASISSET(:)
CHARACTER(len=2) :: SYMTXT,ID3

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,*) '-----------------------------------------'
WRITE(LUPRI,*) '         PRINTING THE MOLECULE.INP FILE '
WRITE(LUPRI,*) '-----------------------------------------'
WRITE(LUPRI,*) '                     '

INQUIRE(file='MOLECULE.INP',EXIST=file_exist) 
IF(file_exist)THEN
  LUINFO=-1
  CALL LSOPEN(LUINFO,'MOLECULE.INP','OLD','FORMATTED')
ELSE
  CALL LSQUIT('MOLECULE.INP does not exist',lupri)
ENDIF

rewind(LUINFO)
DO
   READ (LUINFO, '(A80)') WORD
   IF (WORD(1:1) == '!' .or. WORD(1:1) == '#') CYCLE
   SELECT CASE(WORD)
   CASE('BASIS') 
      WRITE(LUPRI,'(2X,A40)') WORD 
      DO
         READ (LUINFO, '(A80)') WORD
         IF (WORD(1:1) == '!' .or. WORD(1:1) == '#') CYCLE
         WRITE(LUPRI,'(2X,A40)') WORD 
         EXIT
      ENDDO
   CASE DEFAULT
      WRITE(LUPRI,'(2X,A40)') WORD     
   END SELECT
   READ (LUINFO, '(A80)') WORD
   WRITE(LUPRI,'(2X,A40)') WORD     
   READ (LUINFO, '(A80)') WORD
   WRITE(LUPRI,'(2X,A40)') WORD     
   EXIT
ENDDO

READ (LUINFO, '(a80)') WORD
! Number of different atom types
IPOS = INDEX(WORD,'Ato')
IF (IPOS .NE. 0) THEN
   IPOS2 = INDEX(WORD(IPOS:80),'=')
   IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 10)) THEN
      WRITE (LUPRI,'(2X,A50)') 'Incorrect input for # atomtypes'
      WRITE (LUPRI,'(2X,A50)') 'Format is "Atomtypes=?"'
      CALL LSQUIT('Incorrect input for # atomtypes',lupri)
   ELSE
      READ (WORD((IPOS+IPOS2):80),*) Atomtypes
   ENDIF

   IF (Atomtypes.EQ.0) THEN
      WRITE (LUPRI,'(A)')&
           &' You have specified a molecule with zero atoms,&
           & thus all answers to all your input are zero! &
           & (or you made an input error in the .mol file)'
      CALL LSQUIT('No atoms according to .mol input!',lupri)
   ELSE IF (Atomtypes.LT.0) THEN
      WRITE (LUPRI,'(/A,I6)')&
           &' >>> READI1 error, no. of atomic types negative:',Atomtypes
      CALL LSQUIT('Negative number of atoms according to .mol input',lupri)
   ENDIF
ELSE
   !*******************************************************************
   !*  OLD INDPUT FORMAT
   !******************************************************************
   READ (WORD,'(BN,A1,I4,I3,A2,10A1,D10.2,6I5)',IOSTAT=ios) CRT,&
        & Atomtypes,MolecularCharge,SYMTXT,((KASYM(I,J),I=1,3),J=1,3),&
        & ID3,IntegralThreshold
   IF(IOS .NE. 0) THEN
      WRITE (LUPRI,'(A)') ' Error in the determination of the number &
           & of atomic types'
      WRITE (LUPRI,*) "Correct input structure is: Atomtypes=???"
      CALL LSQUIT('Error in determining the number of different &
           & atom types',lupri)
   ENDIF
ENDIF

   WRITE(LUPRI,'(2X,A60)') WORD
   
   DO I=1,Atomtypes
      
      READ (LUINFO, '(a80)') WORD
      IPOS = INDEX(WORD,'Ato')
      IF (IPOS .NE. 0) THEN
         IPOS2 = INDEX(WORD(IPOS:),'=')
         IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
            WRITE (LUPRI,'(2X,A50)') 'Incorrect input for # of atoms'
            WRITE (LUPRI,'(2X,A50)') 'Format is "Atoms=?"'
            CALL LSQUIT('Incorrect input for # of atoms',lupri)
         ELSE
            READ (WORD((IPOS+IPOS2):),*) nAtoms
         ENDIF
      ELSE ! OLD INPUT STYLE
         READ (WORD,'(BN,F10.0,I5)',IOSTAT=ios) Q,nAtoms
         IF(IOS .NE. 0) THEN
            WRITE (LUPRI,'(2X,A50,I4)') ' Error in the determination of the number &
                 & of atoms  for type ',I
            WRITE (LUPRI,'(2X,A50)') "Correct input structure is: Atomtypes=???"
            CALL LSQUIT('Error in determining the number of atoms',lupri)
         ENDIF
      ENDIF
      WRITE(LUPRI,'(2X,A60)')  WORD
      
      IF(nAtoms.GT.10)THEN
         DO J=1,10
            READ (LUINFO, '(a80)') WORD
            WRITE(LUPRI,'(2X,A60)') WORD
         ENDDO
         WRITE(LUPRI,'(2X,A60)') 'WARNING: SINCE YOU HAVE MORE THAN 10 ATOMS&
              & OF THIS'
         WRITE(LUPRI,'(2X,A60)') 'TYPE I WILL NOT PRINT THEM ALL TO LIMIT OUTPUT'
         DO J=11,nAtoms
            READ (LUINFO, '(a80)') WORD
         ENDDO
      ELSE
         DO J=1,nAtoms
            READ (LUINFO, '(a80)') WORD
            WRITE(LUPRI,'(2X,A60)') WORD
         ENDDO
      ENDIF
   ENDDO
   
   CALL LSCLOSE(LUINFO,'KEEP')

END SUBROUTINE PRINT_MOL_FILE

!SUBROUTINE ccfragment(DALTON,ATOMS,NATOMS,FRAGMENT,LUPRI,IPRINT)
!implicit none
!TYPE(DALTONINPUT) :: dalton
!TYPE(DALTONINPUT) :: fragment
!integer           :: LUPRI,NATOMS,IPRINT
!integer           :: ATOMS(NATOMS)
!CHARACTER(len=9)  :: BASISLABEL
!
!CALL io_init(FRAGMENT%IO)
!CALL READ_DALTON_FILE_INTEGRALS_LINSCA(LUPRI,FRAGMENT%DALTON,.FALSE.)
!
!NULLIFY(FRAGMENT%BASIS)
!ALLOCATE(FRAGMENT%BASIS)
!CALL Copy_basissetinfo(dalton%BASIS%REGULAR,fragment%BASIS%REGULAR)
!
!!The values of daltonitem which are set in readmolfile
!FRAGMENT%DALTON%AUXBASIS = DALTON%DALTON%AUXBASIS
!FRAGMENT%DALTON%DoSpherical  = DALTON%DALTON%DoSpherical 
!
!IF(FRAGMENT%DALTON%AUXBASIS)THEN
!   CALL Copy_basissetinfo(dalton%BASIS%AUXILIARY,fragment%BASIS%AUXILIARY)
!ELSE
!   fragment%BASIS%AUXILIARY%natomtypes=0
!ENDIF
!fragment%BASIS%HUCKEL%natomtypes=0
!
!NULLIFY(FRAGMENT%MOLECULE)
!ALLOCATE(FRAGMENT%MOLECULE)
!CALL BUILD_FRAGMENT(DALTON%MOLECULE,FRAGMENT%MOLECULE,fragment%BASIS,&
!     &FRAGMENT%DALTON%AUXBASIS,ATOMS,nATOMS,LUPRI)
!
!IF(IPRINT .GT. 0)THEN
!   CALL PRINT_MOLECULEINFO(LUPRI,FRAGMENT%MOLECULE,FRAGMENT%BASIS)
!   CALL PRINT_MOLECULE_AND_BASIS(LUPRI,FRAGMENT%MOLECULE,FRAGMENT%BASIS%REGULAR)
!   IF(FRAGMENT%DALTON%AUXBASIS)THEN
!      CALL PRINT_MOLECULE_AND_BASIS(LUPRI,FRAGMENT%MOLECULE,FRAGMENT%BASIS%AUXILIARY)
!   ENDIF
!ENDIF
!
!FRAGMENT%DALTON%exchangeFactor = 1.0d0
!
!END SUBROUTINE ccfragment
!
!SUBROUTINE build_ccfragment_set_setting(LS,ATOMS,NATOMS,LUPRI,IPRINT)
!implicit none
!TYPE(LSITEM)       :: LS
!TYPE(MOLECULEINFO),target :: fragment
!TYPE(BASISINFO),target    :: fragmentbasis
!integer            :: LUPRI,NATOMS,IPRINT,iAO
!integer            :: ATOMS(NATOMS)
!CHARACTER(len=9)   :: BASISLABEL
!
!CALL Copy_basissetinfo(ls%input%BASIS%REGULAR,fragmentBASIS%REGULAR)
!
!CALL BUILD_FRAGMENT(ls%input%MOLECULE,FRAGMENT,fragmentBASIS,.FALSE.,ATOMS,nATOMS,LUPRI)
!
!CALL II_set_default_setting(LS%SETTING,LS%INPUT)
!
!DO iAO=1,LS%SETTING%nAO
!   LS%SETTING%MOLECULE(iAO)%p => fragment
!   LS%SETTING%BASIS(iAO)%p => fragmentBASIS
!   LS%SETTING%FRAGMENT(iAO)%p => fragment
!ENDDO
!LS%SETTING%sameMOL  = .TRUE.
!LS%SETTING%sameBAS  = .TRUE.
!LS%SETTING%sameFRAG = .TRUE.
!
!END SUBROUTINE build_ccfragment_set_setting
!
!> \brief build a lsitem as a fragment of a full lsitem, used in DEC
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE build_ccfragmentlsitem(LSFULL,FRAGMENT,ATOMS,NATOMS,LUPRI,IPRINT)
implicit none
!> full lsitem structure 
TYPE(LSITEM)       :: LSFULL
!> fragment lsitem structure to be build 
TYPE(LSITEM)       :: FRAGMENT
!> list of atoms in the fragment
integer            :: ATOMS(NATOMS)
!> number of atoms in the fragment
integer            :: NATOMS
!> the logical unit number for the output file
integer            :: LUPRI
!> the printlevel integer, determining how much output should be generated
integer            :: IPRINT
!
integer            :: IAO
CHARACTER(len=9)   :: BASISLABEL

CALL io_init(FRAGMENT%INPUT%IO)
FRAGMENT%INPUT%DALTON = lsfull%INPUT%DALTON

NULLIFY(FRAGMENT%INPUT%BASIS)
ALLOCATE(FRAGMENT%INPUT%BASIS)
CALL Copy_basissetinfo(lsfull%INPUT%BASIS%REGULAR,fragment%INPUT%BASIS%REGULAR)

!The values of daltonitem which are set in readmolfile
fragment%input%DALTON%AUXBASIS = lsfull%input%DALTON%AUXBASIS
FRAGMENT%input%dalton%DoSpherical  = lsfull%input%DALTON%DoSpherical 

IF(FRAGMENT%INPUT%DALTON%AUXBASIS)THEN
   CALL Copy_basissetinfo(lsfull%input%BASIS%AUXILIARY,fragment%input%BASIS%AUXILIARY)
ELSE
   fragment%input%BASIS%AUXILIARY%natomtypes=0
ENDIF
fragment%input%BASIS%HUCKEL%natomtypes=0

NULLIFY(FRAGMENT%INPUT%MOLECULE)
ALLOCATE(FRAGMENT%INPUT%MOLECULE)
CALL BUILD_FRAGMENT(lsfull%input%MOLECULE,FRAGMENT%input%MOLECULE,fragment%input%BASIS,&
     &FRAGMENT%INPUT%DALTON%AUXBASIS,ATOMS,nATOMS,LUPRI)

IF(IPRINT .GT. 0)THEN
   CALL PRINT_MOLECULEINFO(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS)
   CALL PRINT_MOLECULE_AND_BASIS(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS%REGULAR)
   IF(FRAGMENT%INPUT%DALTON%AUXBASIS)THEN
      CALL PRINT_MOLECULE_AND_BASIS(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS%AUXILIARY)
   ENDIF
ENDIF

FRAGMENT%input%DALTON%exchangeFactor = 1.0d0

CALL II_init_setting(FRAGMENT%setting)
CALL II_set_default_setting(FRAGMENT%SETTING,FRAGMENT%INPUT)
fragment%fopen = .false.

END SUBROUTINE build_ccfragmentlsitem

END MODULE DALTONINFO
