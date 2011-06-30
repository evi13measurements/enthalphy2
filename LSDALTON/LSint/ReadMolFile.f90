!> @file
!> Read molecule input
MODULE READMOLEFILE
  use precision
  use TYPEDEF
  use molecule_module
  use fundamental
contains
!> \brief read the molecule file and build the molecule structure 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param molecule the molecule structure to be built
!> \param BASISSETLIBRARY the info about basisset read from mol file
!> \param doprint if the information should be printet
!> \param iprint the printlevel integer, determining how much output should be generated
!> \param DoSpherical if the basis should be a Spherical or cartesian basis
!> \param integralthreshold Obsolete not really used 
!> \param auxbasis if an auxilliary basis is specified in mol file
!>
!> THIS IS THE MAIN DRIVER ROUTINE WHICH READS THE MOLECULE INPUT FILE
!> AND BUILDS THE MOLECULE OBJECT 
!>
!> THE RECIPE:
!> 1. OPEN MOLECULE.INP FILE -> LUINFO
!> 2. READ THE FIRST LINE CONTAINING KEYWORD BASIS,ATOMBASIS,...
!> 3. IN CASE OF BASIS READ THE SECOND LINE CONTAINING BASISSET
!> 4. READ THE 2 LINES OF COMMENTS
!> 5. READ THE FOURTH LINE CONTAINING NUMBER OF ATOMS, MOLECULARCHARGE,..    
!> 6. CALL 'READ_GEOMETRY' WHICH READS THE FIFTH LINE A NUMBER OF TIME 
!>    AND READS THE X,Y,Z COORDINATES OF THE INDIVIDUAL ATOMS
!>
SUBROUTINE READ_MOL_FILE_AND_BUILD_MOLECULE(LUPRI,MOLECULE,&
     &BASISSETLIBRARY,doprint,iprint,DoSpherical,IntegralThreshold,auxbasis)
implicit none
INTEGER,intent(in)               :: LUPRI,iprint
TYPE(MOLECULEINFO),intent(inout) :: MOLECULE
LOGICAL,intent(inout) :: AUXBASIS,DoSpherical
LOGICAL,intent(in) :: doprint
REAL(REALK),intent(inout)        :: IntegralThreshold
TYPE(BASISSETLIBRARYITEM),intent(inout) :: BASISSETLIBRARY
!
integer            :: LUINFO
logical            :: PRINTATOMCOORD,file_exist,Angstrom,Symmetry,dopbc
logical            :: ATOMBASIS,BASIS
CHARACTER(len=80)  :: BASISSET,AUXBASISSET
integer            :: MolecularCharge,Atomtypes,Totalnatoms
real(realk),PARAMETER :: AngstromConvert = 0.5291772083D0
!CHARACTER(len=80),ALLOCATABLE  :: ATOMBASISSET(:)

IF (doprint) WRITE(LUPRI,'(2X,A18,2X,I6)')'MOLPRINT IS SET TO',IPRINT
BASIS=.FALSE.
ATOMBASIS=.FALSE.
AUXBASIS=.FALSE.

INQUIRE(file='MOLECULE.INP',EXIST=file_exist) 
IF(file_exist)THEN
  LUINFO=-1
  CALL LSOPEN(LUINFO,'MOLECULE.INP','OLD','FORMATTED')
ELSE
  CALL LSQUIT('MOLECULE.INP does not exist',lupri)
ENDIF

CALL Obtain_Totalnatoms(LUPRI,LUINFO,Totalnatoms)

CALL init_MoleculeInfo(Molecule,Totalnatoms,'INPUT-Molecule________')

PRINTATOMCOORD = .FALSE.
IF(Totalnatoms .LE. 11)PRINTATOMCOORD = .TRUE.
CALL READ_LINE1(LUPRI,LUINFO,BASIS,ATOMBASIS)

IF(BASIS)THEN
  CALL READ_LINE2(LUPRI,LUINFO,IPRINT,BASISSET,AUXBASISSET,AUXBASIS)
ENDIF

IF(BASIS .AND. AUXBASIS) THEN
   BASISSETLIBRARY%BASISSETNAME(1)=BASISSET
   BASISSETLIBRARY%BASISSETNAME(2)=AUXBASISSET
   BASISSETLIBRARY%nbasissets=2
ELSEIF(BASIS)THEN
   BASISSETLIBRARY%BASISSETNAME(1)=BASISSET
   BASISSETLIBRARY%nbasissets=1
ELSE !ATOMBASIS
   !We dont know yet 
ENDIF

CALL READ_COMMENTS(LUPRI,LUINFO,.FALSE.)

CALL READ_LINE4(LUPRI,LUINFO,Atomtypes,DoSpherical,MolecularCharge&
                           &,Angstrom,Symmetry,IntegralThreshold,doprint)

Molecule%charge = MolecularCharge

IF (Angstrom.AND.doprint)THEN
  WRITE (LUPRI,'(2X,A)')'Coordinates are entered in Angstroms&
                         & and converted to atomic units.'
  WRITE (LUPRI,'(2X,A,F11.8,A2)')'Conversion factor : 1 bohr ='&
                               &,AngstromConvert,' A'
ENDIF

!************************************************
!*****       Read geometry input data       *****
!************************************************
DOPBC=.FALSE.

CALL READ_GEOMETRY(LUPRI,LUINFO,IPRINT,BASISSETLIBRARY,Atomtypes,dopbc,&
     &ATOMBASIS,BASIS,AUXBASIS,Angstrom,MOLECULE,PRINTATOMCOORD,doprint)

CALL DETERMINE_NELECTRONS(Molecule,Molecule%charge,Molecule%nelectrons)

!CALL PRINT_MOLECULEINFO(LUPRI,MOLECULE,BASISSETLIBRARY)

!CALL PRINT_BASISSETLIBRARY(LUPRI,BASISSETLIBRARY)

CALL LSCLOSE(LUINFO,'KEEP')

END SUBROUTINE READ_MOL_FILE_AND_BUILD_MOLECULE

!> \brief determine the total number of atoms
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param totalnatoms the number of total atoms to be determined
SUBROUTINE Obtain_Totalnatoms(LUPRI,LUINFO,Totalnatoms)
implicit none
integer          :: Totalnatoms,ipos,ipos2,natoms,I,J
integer          :: LUINFO,Atomtypes,numbertypes
CHARACTER(len=80):: LINE
LOGICAL          :: FOUND,BASIS
INTEGER          :: LUPRI
INTEGER          :: MolecularCharge,ios
CHARACTER(len=2) :: SYMTXT,ID3
CHARACTER(len=1) :: KASYM(3,3),CRT
real(realk)      :: IntegralThreshold,AtomicCharge

rewind(LUINFO)
BASIS=.FALSE.

DO
  READ (LUINFO,'(a80)') LINE
  IF (LINE(1:1) == '!' .or. LINE(1:1) == '#') CYCLE
  EXIT
ENDDO

IPOS = INDEX(LINE,'ATOM')
IF (IPOS .NE. 0) THEN !ATOMBASIS
  IPOS2 = INDEX(LINE(IPOS:80),'BASIS')
  IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. (IPOS+10) )) THEN
     WRITE (LUPRI,'(2X,a10,a70)') ' Keyword ',LINE
     WRITE (LUPRI,'(2X,a80)') ' not recognized in READ_MOLECULE_FILE.'
     WRITE (LUPRI,'(2X,a80)') ' maybe it is not implemented yet.'
     WRITE (LUPRI,'(2X,a80)') ' format is ATOMBASIS'
     CALL LSQUIT('Only ATOMBASIS and BASIS is implemented',lupri)
  ENDIF
ELSE
  IPOS2 = INDEX(LINE,'BASIS')
  IF (IPOS2 .EQ. 0) THEN
     WRITE (LUPRI,'(2X,a10,a70)') ' Keyword ',LINE
     WRITE (LUPRI,'(2X,a80)') ' not recognized in READ_MOLECULE_FILE.'
     WRITE (LUPRI,'(2X,a80)') ' maybe it is not implemented yet.'
     WRITE (LUPRI,'(2X,a80)') ' format is BASIS'
     CALL LSQUIT('Only ATOMBASIS and BASIS is implemented',lupri)
  ELSE
    BASIS=.TRUE.
  ENDIF
ENDIF

IF(BASIS) THEN
   DO
      READ (LUINFO,'(a80)') LINE  !BASISSETLINE
      IF (LINE(1:1) == '!' .or. LINE(1:1) == '#') CYCLE
      EXIT
   ENDDO
ENDIF

CALL READ_COMMENTS(LUPRI,LUINFO,.FALSE.)

Totalnatoms=0
FOUND=.FALSE.
DO   
  IF (FOUND) EXIT
  READ (LUINFO,'(a80)') LINE
  IPOS = INDEX(LINE,'Atomtypes')
  IF (IPOS .EQ. 0) IPOS = INDEX(LINE,'AtomTypes')
  IF (IPOS .NE. 0) THEN
    IPOS2 = INDEX(LINE(IPOS:80),'=')
    IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 10)) THEN
      WRITE (LUPRI,'(2X,A40)') 'Incorrect input for # atomtypes'
      WRITE (LUPRI,'(2X,A40)') 'Format is "Atomtypes=?"'
      CALL LSQUIT('Incorrect input for # atomtypes',lupri)
    ELSE
      READ (LINE((IPOS+IPOS2):80),*) Atomtypes
      FOUND=.TRUE.
    ENDIF
  ELSE !ASSUME OLD INPUT
     READ (LINE,'(BN,A1,I4,I3,A2,10A1,D10.2,6I5)',IOSTAT=ios) CRT,&
          & Atomtypes,MolecularCharge,SYMTXT,((KASYM(I,J),I=1,3),J=1,3),&
          & ID3,IntegralThreshold
     IF(IOS .NE. 0) THEN
        WRITE (LUPRI,'(2X,A40)') ' Error in the determination of the number &
             & of atomic types'
        WRITE (LUPRI,'(2X,A60)') "Correct input structure is: Atomtypes=???"
        CALL LSQUIT('Error in determining the number of different &
             & atom types',lupri)
     ELSE
     ENDIF
     EXIT
  ENDIF
ENDDO

numbertypes=1
DO 
  IF(numbertypes>Atomtypes) EXIT
  READ (LUINFO,'(a80)') LINE
  IPOS = INDEX(LINE,'Atoms')
  IF (IPOS .NE. 0) THEN
    IPOS2 = INDEX(LINE(IPOS:),'=')
    IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
      WRITE (LUPRI,'(2X,A40)') 'Incorrect input for # of atoms'
      WRITE (LUPRI,'(2X,A40)') 'Format is "Atoms=?"'
      CALL LSQUIT('Incorrect input for # of atoms',lupri)
    ELSE
      READ (LINE((IPOS+IPOS2):),*) nAtoms
      numbertypes=numbertypes+1
    ENDIF
    Totalnatoms=Totalnatoms+nAtoms
    DO J=1,nAtoms
       READ (LUINFO, '(a80)') LINE
    ENDDO
  ELSE !OLD INPUT
     READ (LINE,'(BN,F10.0,I5)',IOSTAT=ios) AtomicCharge, nAtoms
     IF(IOS .NE. 0) THEN
        WRITE (LUPRI,'(2X,A40)') ' Error in the determination of the number of atoms'
        WRITE (LUPRI,'(2X,A40)') "Correct input structure is: Atomtypes=???"
        CALL LSQUIT('Error in determining the number of atoms',lupri)
     ENDIF
     numbertypes=numbertypes+1
     Totalnatoms=Totalnatoms+nAtoms
     DO J=1,nAtoms
        READ (LUINFO, '(a80)') LINE
     ENDDO
  ENDIF
ENDDO
!CALL LSCLOSE(LUINFO,'KEEP')

END SUBROUTINE Obtain_Totalnatoms

!> \brief read the first line of molecule file (BASIS or ATOMBASIS)
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param basis if the BASIS is used
!> \param atombasis if the ATOMBASIS is used
SUBROUTINE READ_LINE1(LUPRI,LUINFO,BASIS,ATOMBASIS)
implicit none
CHARACTER(len=80)  :: WORD
LOGICAL            :: BASIS,ATOMBASIS
INTEGER            :: LUINFO,IPOS,IPOS2,LUPRI

rewind(LUINFO)
DO
  READ (LUINFO,'(a80)') WORD
  IF(WORD(1:1) == '!' .OR. WORD(1:1) == '#') CYCLE
  EXIT
ENDDO

BASIS=.FALSE.
ATOMBASIS=.FALSE.

IPOS = INDEX(WORD,'ATOM')
IF (IPOS .NE. 0) THEN !ATOMBASIS
  IPOS2 = INDEX(WORD(IPOS:80),'BASIS')
  IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. (IPOS+10) )) THEN
     WRITE (LUPRI,'(2X,a10,a70)') ' Keyword ',WORD
     WRITE (LUPRI,'(2X,a80)') ' not recognized in READ_MOLECULE_FILE.'
     WRITE (LUPRI,'(2X,a80)') ' maybe it is not implemented yet.'
     WRITE (LUPRI,'(2X,a80)') ' format is ATOMBASIS'
     CALL LSQUIT('Only ATOMBASIS and BASIS is implemented',lupri)
  ELSE
     ATOMBASIS=.TRUE.
  ENDIF
ELSE
  IPOS2 = INDEX(WORD,'BASIS')
  IF (IPOS2 .NE. 0) THEN
     BASIS=.TRUE.
  ELSE
     WRITE (LUPRI,'(2X,a10,a70)') ' Keyword ',WORD
     WRITE (LUPRI,'(2X,a80)') ' not recognized in READ_MOLECULE_FILE.'
     WRITE (LUPRI,'(2X,a80)') ' maybe it is not implemented yet.'
     WRITE (LUPRI,'(2X,a80)') ' format is BASIS'
     CALL LSQUIT('Only ATOMBASIS and BASIS is implemented',lupri)
  ENDIF
ENDIF

END SUBROUTINE READ_LINE1

!> \brief read the second line of molecule file (the name of basisset) if BASIS is used 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param iprint the printlevel integer, determining how much output should be generated
!> \param basisset the name of the basisset 
!> \param auxbasisset the name of the auxilliary basisset 
!> \param auxbasis if an auxilliary basisset is given
SUBROUTINE READ_LINE2(LUPRI,LUINFO,IPRINT,BASISSET,AUXBASISSET,AUXBASIS)
implicit none
CHARACTER(len=80)  :: WORD,BASISSET,AUXBASISSET
INTEGER            :: N,I,K,Itmp,M,L
LOGICAL            :: NEXT,AUXBASIS
INTEGER            :: LUINFO
INTEGER            :: LUPRI,IPRINT
AUXBASIS=.FALSE.
N=0
DO I=1,80
  BASISSET(I:I) = ' '
  AUXBASISSET(I:I) = ' '
ENDDO
READ (LUINFO, '(a80)') WORD
NEXT = .FALSE.
DO I=1,80
  IF (WORD(I:I) .NE. ' ' .AND. .NOT. NEXT) THEN
    NEXT=.TRUE.  !INSIDE WORD
    k=i          !Start of word
  ELSEIF(WORD(I:I) .EQ. ' ' .AND. NEXT) THEN
    L=I-1        !end of WORD
    NEXT=.FALSE. !OUT OF WORD
    N=N+1
    IF(N==1)THEN
      BASISSET(1:I-K)=WORD(K:L) !FOUND BASISSET
    ENDIF
    IF(N==2) THEN !FOUND SECOND BASISSET
      IF(WORD(K:K) .EQ. '!')EXIT !not a basisset, a commet
      IF(WORD(K:K) .EQ. '#')EXIT !not a basisset, a commet
      IF(WORD(K:K+2) .EQ. 'Aux')THEN
        AUXBASIS=.TRUE.
      ENDIF
!     Read a portion of the form "Aux = AUXBASNAM", with arbitrarily many
!     blanks. This is assumed to end the line, so we RETURN afterwards.
!     1. find '='
      Itmp = K + 3
      DO WHILE (WORD(Itmp:Itmp) .ne. '=' .and. Itmp .lt. 80)
         Itmp = itmp + 1
      ENDDO
      Itmp = Itmp + 1
      IF (Itmp .ge. 80) CALL LSquit('"Aux" given but no basis provided!',lupri)
!     2. find non-blank
      DO WHILE (WORD(Itmp:Itmp) .eq. ' ' .and. itmp .lt. 80)
        Itmp = Itmp + 1
      ENDDO
      If (itmp .ge. 80) CALL LSquit('"Aux" given but no basis provided!',lupri)
!     3. copy to next blank
      M = 1
      DO WHILE(WORD(itmp:itmp) .ne. ' ' .and. itmp .lt. 80)
        AUXBASISSET(M:M) = WORD(itmp:itmp)
        M=M+1
        Itmp = Itmp + 1
      ENDDO
    ENDIF
    IF(N>2)THEN
        WRITE (LUPRI,'(/A)') ' Found string"'//WORD(K:L)&
        &//'" but not correct format in MOLECULE.INP.'
    ENDIF
  ENDIF
ENDDO
IF (IPRINT .GT. 1) THEN
  WRITE (LUPRI,'(/A)') ' Basis set "'//BASISSET(1:I-K)&
        &//'" from the basis set library will be used.'
  IF(AUXBASIS)THEN
     WRITE (LUPRI,'(/A)') ' AuxBasis set "'//AUXBASISSET(1:M)&
          &//'" from the basis set library will be used.'
  ENDIF
END IF
END SUBROUTINE READ_LINE2

!> \brief find N strings (seperated by empty ' ' characters) in the word
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param N the number of strings in the word
!> \param word the character string
SUBROUTINE FIND_N_STRING(LUPRI,N,WORD)
CHARACTER(len=80)  :: RWORD,WORD
INTEGER            :: I,K,L,N,M
LOGICAL            :: NEXT
INTEGER            :: LUPRI

DO I=1,80
  RWORD(I:I) = WORD(I:I)
  WORD(I:I) = ' '
ENDDO
M=0
NEXT=.FALSE.
DO I=1,80
  IF (RWORD(I:I) .NE. ' ' .AND. .NOT. NEXT) THEN
    NEXT=.TRUE.  !INSIDE WORD
    K=I          !Start of word
  ELSEIF(RWORD(I:I) .EQ. ' ' .AND. NEXT) THEN
    L=I-1        !end of WORD
    NEXT=.FALSE. !OUT OF WORD
    M=M+1   
    IF(N==M)THEN
      WORD(1:I-K)=RWORD(K:L)
    ENDIF
    RETURN
  ENDIF
END DO
IF (N<M) THEN
 WRITE(LUPRI,*)'TRYING TO READ MORE LABELS THEN IN STRING'
 CALL LSQUIT('TO FEW LABELS IN FIND WORD.',lupri)
ENDIF

END SUBROUTINE FIND_N_STRING

!> \brief read the 2 comment lines in the molecule input file
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param printing should the comment lines be printed to output file
SUBROUTINE READ_COMMENTS(LUPRI,LUINFO,PRINTING)
implicit none
INTEGER  :: LUINFO,LUPRI
LOGICAL  :: PRINTING
CHARACTER(len=60):: LINE 

READ (LUINFO, '(2X,a80)') LINE
IF (PRINTING) WRITE(LUPRI,*)LINE
READ (LUINFO, '(2X,a80)') LINE
IF (PRINTING) WRITE(LUPRI,*)LINE

END SUBROUTINE READ_COMMENTS

!> \brief read lin 4 in the molecule input file (Atomtypes= Charge= and so on)
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param Atomtypes the number of different atomtypes
!> \param DoSpherical should the basis be spherical or cartesian
!> \param molecularcharge the molecular charge
!> \param angstrom is the coordinates given in angstrom or bohr (au)
!> \param symmetry if we should use symmetry - obsolete 
!> \param IntegralThreshold integral threshold - obsolete
!> \param doprint if we should print this to output files
SUBROUTINE READ_LINE4(LUPRI,LUINFO,Atomtypes,DoSpherical,MolecularCharge&
                                    &,Angstrom,Symmetry,IntegralThreshold,doprint)
! Read in the fourth molecule.inp line using the new (and old) input scheme
implicit none
INTEGER          :: IPOS,IPOS2,Atomtypes,MolecularCharge,IPSO
LOGICAL          :: DoSpherical,Angstrom,Symmetry,doprint
CHARACTER(len=2) :: SYMTXT,ID3
CHARACTER(len=80):: LINE 
INTEGER          :: LUINFO,IOS,i,j
CHARACTER(len=1) :: KASYM(3,3),CRT
real(realk)      :: IntegralThreshold
INTEGER          :: LUPRI

READ (LUINFO, '(a80)') LINE
  ! Number of different atom types
  IPOS = INDEX(LINE,'Ato')
  IF (IPOS .NE. 0) THEN
    IPOS2 = INDEX(LINE(IPOS:80),'=')
    IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 10)) THEN
      WRITE (LUPRI,'(2X,A40)') 'Incorrect input for # atomtypes'
      WRITE (LUPRI,'(2X,A40)') 'Format is "Atomtypes=?"'
      CALL LSQUIT('Incorrect input for # atomtypes',lupri)
    ELSE
      READ (LINE((IPOS+IPOS2):80),*) Atomtypes
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

!
!     Kind of basis set
!     

    DoSpherical = .TRUE.
    IPOS = INDEX(LINE,'Car')
    IF (IPOS .NE. 0) DoSpherical = .FALSE.

!    DOSPHERICAL = .TRUE.
!    IPOS = INDEX(LINE,'NoSph')
!    IF (IPOS .EQ. 0) THEN
!      DOSPHERICAL = .FALSE.
!    ENDIF


!
!     Charge of molecule
!     
    IPOS = INDEX(LINE,'Cha')
    IF (IPOS .NE. 0) THEN
      IPOS2 = INDEX(LINE(IPOS:80),'=')
      IF ((IPOS2 .EQ. 0) .OR. (IPOS2 .GT. 7)) THEN
        WRITE (LUPRI,'(2X,A40)') 'Incorrect input for molecular charge'
        WRITE (LUPRI,'(2X,A40)') 'Format is "Charge=?"'
        CALL LSQUIT('Incorrect input for molecular charge',lupri)
      ELSE
        READ (LINE(IPOS+IPOS2:80),*) MolecularCharge
      ENDIF
    ELSE
        MolecularCharge = 0
    ENDIF

!     
!     Angstrom?
!
    IPOS = INDEX(LINE,'Ang') 
    IF (IPOS .NE. 0) THEN
      Angstrom = .TRUE.
    ELSE
      Angstrom = .FALSE.
    ENDIF

!     
!     Symmetry generators is disregarded at the moment but the old code to read
!     is in abacus/herrdn.F in subroutine LINE4
!   Symmetry turned off always
    SYMMETRY = .FALSE.
!
!     Change of integral threshold?
!
    IPOS = INDEX(LINE,'Int')
    IF (IPOS .NE. 0) THEN
      IPOS2 = INDEX(LINE(IPOS:80),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 11)) THEN
        WRITE (LUPRI,'(2X,A40)') 'Incorrect input for integral threshold'
        WRITE (LUPRI,'(2X,A40)') 'Format is "Integrals=?"'
        CALL LSQUIT('Incorrect input for integral threshold',lupri)
      ELSE
        !SUBROUTINE BY KENNETH THAT SIMULATES A FREEFORMAT READ
!        CALL FREFRM(LINE,IPOS+IPOS2,IDUMMY,IntegralThreshold,'REA',IERR)
        READ (LINE((IPOS+IPOS2):),*) IntegralThreshold
        IF(IntegralThreshold>1.00D+0)THEN
        WRITE (LUPRI,'(2X,A40)') 'Incorrect input for integral threshold'
        WRITE (LUPRI,'(2X,A40)') 'Format is "Integrals=1.00D-15"'
        WRITE (LUPRI,'(2X,A40)') 'Your choice was: ',IntegralThreshold
        CALL LSQUIT('Incorrect input for integral threshold',lupri)
        ENDIF       
      ENDIF
    ELSE
      IF (doprint) WRITE(LUPRI,'(2X,A49)')'Integral threshold is set to the default 1.00D-15'
      IntegralThreshold=1.00D-15
    ENDIF
  ELSE

    
!*******************************************************************
!*
!*  OLD INDPUT FORMAT
!*
!******************************************************************

    READ (LINE,'(BN,A1,I4,I3,A2,10A1,D10.2,6I5)',IOSTAT=ios) CRT,&
    & Atomtypes,MolecularCharge,SYMTXT,((KASYM(I,J),I=1,3),J=1,3),&
    & ID3,IntegralThreshold
    IF(IOS .NE. 0) THEN
      WRITE (LUPRI,'(2X,A)') ' Error in the determination of the number &
     & of atomic types'
      WRITE (LUPRI,*) "Correct input structure is: Atomtypes=???"
      CALL LSQUIT('Error in determining the number of different &
     & atom types',lupri)
    ENDIF
    IF (LINE(21:30) .EQ. '          ') THEN
      IF (doprint) WRITE(LUPRI,'(2X,A55)') &
     &'Integral threshold is set to the default 1.00D-15'
      IntegralThreshold=1.00D-15
    ENDIF
    IF (SYMTXT(1:1) .EQ. '0')THEN
       SYMMETRY = .FALSE.
    ELSE IF(SYMTXT(2:2) .EQ. '0') THEN
       SYMMETRY = .FALSE.
    ELSE
       SYMMETRY = .FALSE.       
    ENDIF
    IF(CRT .EQ. 'C' .OR. CRT .EQ. 'c') THEN
      DoSpherical = .FALSE.
    ELSE
      DoSpherical = .true.
    ENDIF
    IF((ID3 .EQ. 'A' .OR. ID3 .EQ. 'X').OR. ID3 .EQ. '*' ) THEN
      Angstrom = .TRUE.
    ELSE
      Angstrom = .FALSE.
    ENDIF
  ENDIF
  IF(.NOT.DoSpherical)call lsquit('at the moment only spherical harmonic basis functionas have been implemented.',lupri)

END SUBROUTINE READ_LINE4

!> \brief read line 5 in the molecule input file CONTAINING ATOMICCHARGE,NUMBER OF ATOMS,ATOMICBASISSET and then x,y,z coordinates.
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param iprint the printlevel, determining how much output should be generated
!> \param BASISSETLIBRARY the info about the basisset required
!> \param Atomtypes the number of different atomtypes
!> \param if the periodic boundry conditions is to be used
!> \param atombasis is the atombasis keyword used in 1. line
!> \param basis is the basis keyword used in 1. line
!> \param auxbasis is the auxilliary basis given in 1. line
!> \param angstrom is the coordinates given in angstrom or bohr (au)
!> \param molecule the molecule to be built
!> \param PRINTATOMCOORD should the coordinates be printed
!> \param doprint if we should print this to output files
!>
!> THIS IS THE MAIN ROUTINE WHICH READS THE ATOMIC INFORMATION 
!>
!> THE RECIPE:
!> 1. READ LINE5 CONTAINING ATOMICCHARGE,NUMBER OF ATOMS,ATOMICBASISSET,..
!> 2. READS THE X,Y,Z COORDINATES OF THE INDIVIDUAL ATOMS
!> 3. IF THE FILE IS WRITTEN IN ANGSTROM WE CONVERT TO ATOMIC UNIT 
!>
SUBROUTINE READ_GEOMETRY(LUPRI,LUINFO,IPRINT,BASISSETLIBRARY,Atomtypes,dopbc,&
     &ATOMBASIS,BASIS,AUXBASIS,Angstrom,MOLECULE,PRINTATOMCOORD,doprint)
  use ls_util
implicit none
TYPE(MOLECULEINFO) :: MOLECULE
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
INTEGER            :: Atomtypes,LUINFO,IPOS,atomnumber
INTEGER            :: I,J,natoms
CHARACTER(len=4)   :: StringFormat
LOGICAL            :: ATOMBASIS,BASIS,Angstrom,dopbc,AUXBASIS,doprint
real(realk)        :: AtomicCharge
real(realk),PARAMETER :: AngstromConvert = 0.5291772083D0
CHARACTER(len=80)  :: LINE,Atomicbasisset,Auxbasisset 
CHARACTER(len=1)   :: CHRXYZ(3)=(/'x','y','z'/)
INTEGER            :: LUPRI,basissetnumber,basissetnumber1,basissetnumber2
INTEGER            :: unique1,unique2,uniqueCharge,IPRINT
LOGICAL            :: PRINTATOMCOORD
atomnumber=0
basissetnumber=0
DO I=1,Atomtypes
 CALL READ_LINE5(LUPRI,LUINFO,AtomicCharge,nAtoms,AtomicBasisset,ATOMBASIS,BASIS,AUXBASIS,AUXBASISSET)

 IF(ATOMBASIS)THEN
    CALL DETERMINE_UNIQUE_BASIS(BASISSETLIBRARY,ATOMICBASISSET,&
                                             &unique1,basissetnumber)
    IF(unique1 == 0)THEN !found new basisset
       basissetnumber=basissetnumber+1
       IF(basissetnumber .GT. maxBasisSetInLIB) THEN
          WRITE(LUPRI,*)'You use many different basisset (which is okay),&
               & but you need to increase maxBasisSetInLIB in TYPE-DEF.f90&
               & to a number equal to or greater than the number of different&
               & basissets (including auxiliary basissets). At the moment it&
               & is set to ', maxBasisSetInLIB,'Which is clearly not enough for you.&
               & Thomas Kjaergaard'
          CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
       ENDIF
       BASISSETLIBRARY%BASISSETNAME(basissetnumber)=ATOMICBASISSET
       BASISSETLIBRARY%Charges(basissetnumber,1)=INT(AtomicCharge)
       BASISSETLIBRARY%nCharges(basissetnumber)=1 
       basissetnumber1=basissetnumber
       BASISSETLIBRARY%nbasissets=basissetnumber
    ELSE !old basisset
       CALL DETERMINE_UNIQUE_CHARGE(BASISSETLIBRARY,unique1,AtomicCharge,uniqueCharge)
       IF(uniqueCharge /= 0)THEN !found new charge

          IF(uniqueCharge .GT. maxNumberOfChargesinLIB) THEN
             WRITE(LUPRI,*)'You use many different charges (which is okay),&
                  & but you need to increase maxNumberOfChargesinLIB in TYPE-DEF.f90&
                  & to a number equal to or greater than the number of different&
                  & charges in the same basisset. So if you have 18 different atoms&
                  & with 6-31G basis and 15 different atoms with STO-3G you need&
                  & to increase maxNumberOfChargesinLIB to 18 or higher. &
                  & maxNumberOfChargesinLIB is set to ', maxNumberOfChargesinLIB,&
                  &'Which is clearly not enough for you.&
                  & Thomas Kjaergaard'
             CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
          ENDIF

          BASISSETLIBRARY%Charges(unique1,uniqueCharge)=INT(AtomicCharge)
          BASISSETLIBRARY%nCharges(unique1)=uniqueCharge
       ELSE
       ENDIF !old basisset old charge
    ENDIF
 
    IF(AUXBASIS)THEN
       CALL DETERMINE_UNIQUE_BASIS(BASISSETLIBRARY,AUXBASISSET,&
                                             &unique2,basissetnumber)
       IF(unique2 == 0)THEN !found new basisset
          basissetnumber=basissetnumber+1
          IF(basissetnumber .GT. maxBasisSetInLIB) THEN
             WRITE(LUPRI,*)'You use many different basisset (which is okay),&
                  & but you need to increase maxBasisSetInLIB in TYPE-DEF.f90&
                  & to a number equal to or greater than the number of different&
                  & basissets (including auxiliary basissets). At the moment it&
                  & is set to ', maxBasisSetInLIB,'Which is clearly not enough for you.&
                  & Thomas Kjaergaard'
             CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
          ENDIF
          BASISSETLIBRARY%BASISSETNAME(basissetnumber)=AUXBASISSET
          BASISSETLIBRARY%Charges(basissetnumber,1)=INT(AtomicCharge)
          BASISSETLIBRARY%nCharges(basissetnumber)=1 
          basissetnumber2=basissetnumber
          BASISSETLIBRARY%nbasissets=basissetnumber
       ELSE !old basisset
          CALL DETERMINE_UNIQUE_CHARGE(BASISSETLIBRARY,unique2,AtomicCharge,uniqueCharge)

          IF(uniqueCharge /= 0)THEN !found new charge
             IF(uniqueCharge .GT. maxNumberOfChargesinLIB) THEN
                WRITE(LUPRI,*)'You use many different charges (which is okay),&
                     & but you need to increase maxNumberOfChargesinLIB in TYPE-DEF.f90&
                     & to a number equal to or greater than the number of different&
                     & charges in the same basisset. So if you have 18 different atoms&
                     & with 6-31G basis and 15 different atoms with STO-3G you need&
                     & to increase maxNumberOfChargesinLIB to 18 or higher. &
                     & maxNumberOfChargesinLIB is set to ', maxNumberOfChargesinLIB,&
                     &'Which is clearly not enough for you.&
                     & Thomas Kjaergaard'
                CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
             ENDIF
             BASISSETLIBRARY%Charges(unique2,uniqueCharge)=INT(AtomicCharge)
             BASISSETLIBRARY%nCharges(unique2)=uniqueCharge
          ELSE
          ENDIF !old basisset old charge
       ENDIF
    ENDIF
 ENDIF


 IF(BASIS)THEN
    unique1=1
    basissetnumber1=1
    unique2=2
    basissetnumber2=2
    IF(I==1)THEN
       BASISSETLIBRARY%Charges(unique1,1)=INT(AtomicCharge)
       BASISSETLIBRARY%nCharges(unique1)=1
       IF(AUXBASIS)THEN
          BASISSETLIBRARY%Charges(unique2,1)=INT(AtomicCharge)
          BASISSETLIBRARY%nCharges(unique2)=1
       ENDIF
    ELSE
       CALL DETERMINE_UNIQUE_CHARGE(BASISSETLIBRARY,unique1,AtomicCharge,uniqueCharge)
       IF(uniqueCharge /= 0)THEN !found new charge
          IF(uniqueCharge .GT. maxNumberOfChargesinLIB) THEN
             WRITE(LUPRI,*)'You use many different atoms/charges (which is okay),&
                  & but you need to increase maxNumberOfChargesinLIB in TYPE-DEF.f90&
                  & to a number equal to or greater than the number of different&
                  & charges. So if you have 18 different atoms&
                  & with 6-31G basis, you need&
                  & to increase maxNumberOfChargesinLIB to 18 or higher. &
                  & maxNumberOfChargesinLIB is set to ', maxNumberOfChargesinLIB,&
                  &'Which is clearly not enough for you.&
                  & Thomas Kjaergaard'
             CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
          ENDIF
          BASISSETLIBRARY%Charges(unique1,uniqueCharge)=INT(AtomicCharge)
          BASISSETLIBRARY%nCharges(unique1)=uniqueCharge
          IF(AUXBASIS)THEN
             BASISSETLIBRARY%Charges(unique2,uniqueCharge)=INT(AtomicCharge)
             BASISSETLIBRARY%nCharges(unique2)=uniqueCharge
          ENDIF
       ENDIF
    ENDIF

 ENDIF
 DO J=1,nAtoms
    atomnumber=atomnumber+1
    MOLECULE%ATOM(atomnumber)%GHOST = .FALSE.
    MOLECULE%ATOM(atomnumber)%Charge=AtomicCharge
    CALL ATTACH_BASISINDEX_AND_BASISLABEL(LUPRI,MOLECULE,&
         &atomnumber,basissetnumber1,basissetnumber2,unique1,unique2,&
         &BASIS,ATOMBASIS,AUXBASIS)
    ! Set atomic numbers, isotopes and masses
    MOLECULE%ATOM(atomnumber)%Isotope = 1
    MOLECULE%ATOM(atomnumber)%Atomic_number = NINT(AtomicCharge)
    MOLECULE%ATOM(atomnumber)%Mass = &
    Isotopes(MOLECULE%ATOM(Atomnumber)%Atomic_number, &
    MOLECULE%ATOM(AtomNumber)%Isotope,'MASS',LUPRI)
    !READ_ATOMCOORD 

    READ (LUINFO, '(a80)') LINE

    IPOS = INDEX(LINE,' ')
    IF (IPOS .LT. 2) THEN
      WRITE (LUPRI,*) 'Atom name must start in first column'
      CALL LSQUIT('Error in placement of atom name. See output',lupri)
    ELSEIF (IPOS .EQ. 2) THEN
      StringFormat = '(A1)'
    ELSEIF (IPOS .EQ. 3) THEN
      StringFormat = '(A2)'
    ELSEIF (IPOS .EQ. 4) THEN
      StringFormat = '(A3)'
    ELSEIF (IPOS .EQ. 5) THEN
      StringFormat = '(A4)'
    ELSE
      WRITE (LUPRI,*) 'Atom name must be less then 4 letters'
      CALL LSQUIT('Error in atom name. See output',lupri)
    END IF

    READ (LINE,StringFormat) MOLECULE%ATOM(atomnumber)%Name

!Read x,y and z coordinate 
    READ (LINE(IPOS:80),*) MOLECULE%ATOM(atomnumber)%CENTER(1), &
    & MOLECULE%ATOM(atomnumber)%CENTER(2), MOLECULE%ATOM(atomnumber)%CENTER(3)

!#if defined(SYS_CRAY) && defined (VAR_NOFREE) && defined (SYS_T3D) && defined (SYS_T90)
!         READ (LINE(IPOS:80),*,ERR=101) &
!                & MOLECULE%ATOM(atomnumber)%CENTER(1)&  !X
!                & MOLECULE%ATOM(atomnumber)%CENTER(2)&  !Y
!                & MOLECULE%ATOM(atomnumber)%CENTER(3)   !Z
!         GOTO 104
!  101    READ (LINE(IPOS:80),'(BN,3F20.0)',ERR=102) & 
!                & MOLECULE%ATOM(atomnumber)%CENTER(1)&  !X
!                & MOLECULE%ATOM(atomnumber)%CENTER(2)&  !Y
!                & MOLECULE%ATOM(atomnumber)%CENTER(3)   !Z
!         GO TO 104
!  102    READ (LINE(IPOS:80),'(BN,3F10.0)',ERR=103)& 
!                & MOLECULE%ATOM(atomnumber)%CENTER(1)&  !X
!                & MOLECULE%ATOM(atomnumber)%CENTER(2)&  !Y
!                & MOLECULE%ATOM(atomnumber)%CENTER(3)   !Z
!         GO TO 104
!  103    CONTINUE
!            WRITE(LUPRI,'(/A,I5/A,I5,A)')&
!     &      ' ERROR: Unable to read Cartesian coordinates of atom no.',&
!     &      atomnumber,' from the following line in the MOLECULE input file:'
!            WRITE(LUPRI,'(A)') LINE
!         CALL QUIT('ERROR reading atomic coordinates in MOLECULE input')
!  104    CONTINUE
!#else
!    CALL FREFRM(LINE,IPOS,IDUMMY,MOLECULE%ATOM(atomnumber)%CENTER(1),'REA',IERR)
!    CALL FREFRM(LINE,IPOS,IDUMMY,MOLECULE%ATOM(atomnumber)%CENTER(2),'REA',IERR)
!    CALL FREFRM(LINE,IPOS,IDUMMY,MOLECULE%ATOM(atomnumber)%CENTER(3),'REA',IERR)
!#endif

    IF (dopbc) THEN
      write(LUPRI,*) 'debug: read atom position'
      write(LUPRI,*) 'PBC not implemented'
!     Here we need to do two things: First, check that
!     the given atom position is within the cell. Second,
!     if the position was given in lattice coordinates we
!     need to transform to standard coordinates, since
!     subsequent calculations assume standard coordinates.
!     ErikT

!     call pbc_transf_atom_coord(NUCIND,CORD(1,NUCIND),atom_coord_lat)
    ENDIF

    IPOS = INDEX(LINE,'Isotope=')
    IF (IPOS .NE. 0) THEN
      IPOS = IPOS + 8
      READ (LINE(IPOS:80),'(I3)') MOLECULE%ATOM(atomnumber)%Isotope
    ELSE
      MOLECULE%ATOM(atomnumber)%Isotope = 1
    END IF

    IF(Angstrom) THEN
      MOLECULE%ATOM(atomnumber)%CENTER(1) = MOLECULE%ATOM(atomnumber)%CENTER(1)/AngstromConvert 
      MOLECULE%ATOM(atomnumber)%CENTER(2) = MOLECULE%ATOM(atomnumber)%CENTER(2)/AngstromConvert 
      MOLECULE%ATOM(atomnumber)%CENTER(3) = MOLECULE%ATOM(atomnumber)%CENTER(3)/AngstromConvert 
    ENDIF

! IF(ABS(CORDinates).GT.CORMAX) THEN
!        WRITE (LUPRI,'(A,E12.6,A/A/A,E12.6)')
!     &    ' Atomic coordinate ',CORD(J,NUCIND),
!     &    ' too large in CNTINP.',
!     &    ' Note: Program is unstable for large coordinates.',
!     &    ' Maximum coordinate value:',CORMAX
!        CALL QUIT('*** ERROR: Atomic coordinate too large in CNTINP')
!ENDIF

!IF(IPRINT) THEN
!  WRITE(LUPRI,'(6X,A4,3F20.15)') NAMN(NUCIND),(CORD(J,NUCIND), J = 1,3)
!ENDIF

!!$C
!!$C     To avoid problems in later stages, we now rewrite the coordinates
!!$C     of the nuclei in traditional Dalton format, with 4 characters for the
!!$C     name of the atom
!!$C
!!$         IF (IPOS .EQ. 0) THEN
!!$            WRITE (MLINE(NMLINE),'(A4,3F20.10,2X,A14)') NAMN(NUCIND),
!!$     &           (CORD(J,NUCIND),J = 1, 3), '                '
!!$         ELSE
!!$            IF (MASSNM .LT. 10) THEN
!!$               WRITE (MLINE(NMLINE),'(A4,3F20.10,A10,I1,A5)') 
!!$     &              NAMN(NUCIND),(CORD(J,NUCIND),J = 1, 3), 
!!$     &              '  Isotope=', MASSNM, '     '
!!$            ELSE IF (MASSNM .LT. 100) THEN
!!$               WRITE (MLINE(NMLINE),'(A4,3F20.10,A10,I2,A4)') 
!!$     &              NAMN(NUCIND),(CORD(J,NUCIND),J = 1, 3), 
!!$     &              '  Isotope=', MASSNM, '    '
!!$            ELSE
!!$               WRITE (MLINE(NMLINE),'(A4,3F20.10,A10,I3,A3)') 
!!$     &              NAMN(NUCIND),(CORD(J,NUCIND),J = 1, 3), 
!!$     &              '  Isotope=', MASSNM, '   '
!!$            END IF
!!$         END IF


!C*****************************************************************************
!C       MULBSI  - basis-set identifier (WK/UniKA/04-11-2002).
!C       CHARGE  - charge of center
!C       NOORBT  - TRUE: no orbitals on this center
!C       GNUEXP  - exponent of Gaussian nuclear charge distribution
!C*****************************************************************************

!         MULBSI(NUCIND) = MAX(MBSI,1)
!         IF (MBSI .GT. 1) THEN
!            CHARGE(NUCIND) = D0
!            IZATOM(NUCIND) = -1
!C           ... code -1 to tell MBSI (value of 0 is floating orbitals)
!            GNUEXP(NUCIND) = D0
!         ELSE
!            CHARGE(NUCIND) = Q
!            IZATOM(NUCIND) = NQ
!C           ... if point charge this is reset outside to -2
!C               to tell this is a point charge /hjaaj Nov 2003
!            GNUEXP(NUCIND) = GEXP
!         END IF
  ENDDO
ENDDO

IF ((IPRINT .GT. 0).AND.doprint) THEN
   WRITE (LUPRI,'(2X,A,I3)')' Total number of atoms:',MOLECULE%natoms
ENDIF

IF(IPRINT .GT. -1 .AND. PRINTATOMCOORD.AND.DOPRINT) THEN
   CALL LSHEADER(LUPRI,'Cartesian Coordinates Linsca (au)')
   CALL PRINT_GEOMETRY(MOLECULE,LUPRI)
ENDIF
END SUBROUTINE READ_GEOMETRY

!> \brief read line 5 in the molecule input file: Read the atom-specific information in new input style
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param AtomicCharge the charge of the atom
!> \param natoms the number of atoms of this type
!> \param Atomicbasisset the basisset for this atom
!> \param atombasis if ATOMBASIS keyword is given in line 1
!> \param basis if BASIS keyword is given in line 1
!> \param auxbasis is the auxilliary basis given in 1. line
!> \param auxbasisset the auxilliary basisset for this atom
SUBROUTINE READ_LINE5(LUPRI,LUINFO,AtomicCharge,nAtoms,AtomicBasisset,ATOMBASIS,BASIS,AuxBASIS,Auxbasisset)
implicit none
real(realk)        :: AtomicCharge
INTEGER            :: IPOS,IPOS2,IPOS3,nAtoms,LUINFO
CHARACTER(len=80)  :: TEMPLINE,Atomicbasisset,Auxbasisset
LOGICAL            :: ATOMBASIS,BASIS,AUXBASIS
CHARACTER(len=5)   :: StringFormat
INTEGER            :: LUPRI,ios

READ (LUINFO, '(a80)') TEMPLINE
IPOS = INDEX(TEMPLINE,'Cha')
IF (IPOS .NE. 0 ) THEN
   IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
   IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 7)) THEN
      WRITE (LUPRI,'(2X,A40)') 'Incorrect input for atomic charge'
      WRITE (LUPRI,'(2X,A40)') 'Format is "Charge=?"'
      CALL LSQUIT('Incorrect input for atomic charge',lupri)
   ELSE
      READ (TEMPLINE((IPOS+IPOS2):),*) AtomicCharge
      IPOS = INDEX(TEMPLINE,'Ato')
      IF (IPOS .NE. 0) THEN
         IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
         IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
            WRITE (LUPRI,'(2X,A40)') 'Incorrect input for # of atoms'
            WRITE (LUPRI,'(2X,A40)') 'Format is "Atoms=?"'
            CALL LSQUIT('Incorrect input for # of atoms',lupri)
         ELSE
            READ (TEMPLINE((IPOS+IPOS2):),*) nAtoms
         ENDIF
      ENDIF
   ENDIF
ELSE !OLD INPUT STYLE
   READ (TEMPLINE,'(BN,F10.0,I5)',IOSTAT=ios) AtomicCharge, nAtoms
   IF(IOS .NE. 0) THEN
      WRITE (LUPRI,'(2X,A50)') ' Error in the determination of the number of atoms'
      WRITE (LUPRI,'(2X,A50)') "Correct input structure is: Atomtypes=???"
      CALL LSQUIT('Error in determining the number of atoms',lupri)
   ENDIF
ENDIF

!     Multiple basis sets used?
!         
!  IF (LMULBS) THEN
!    IPOS = INDEX(TEMPLINE,'Set')
!    IF (IPOS .NE. 0) THEN
!      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
!      IF (IPOS2 .EQ. 0 .OR. (IPOS2. GT. 5)) THEN
!        WRITE (2,*) 'Incorrect input for # of basis sets'
!        WRITE (2,*) 'Format is "Sets=?"'
!        CALL QUIT('Incorrect input for # of basis sets')
!      ELSE
!        READ (TEMPLINE((IPOS+IPOS2):),*) MBSI
!      END IF
!    END IF
!  END IF
! Read in basis set information
! Integral ignored for now
IF (.NOT. (BASIS .OR. ATOMBASIS)) THEN
   IPOS = INDEX(TEMPLINE,'Blo')
   IF (IPOS .NE. 0) THEN
      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 7)) THEN
         WRITE (LUPRI,'(2X,A40)') 'Incorrect input for # of integral blocks'
         WRITE (LUPRI,'(2X,A40)') 'Format is "Blocks=? ? ?"'
         CALL LSQUIT('Incorrect input for # of integral blocks',lupri)
      ELSE
         !        ISTART = IPOS + IPOS2
         !        CALL FREFRM(TEMPLINE,ISTART,IQM,DUMMY,'INT',IERR)
         !        DO IQMLOP = 1, IQM
         !          CALL FREFRM(TEMPLINE,ISTART,JCO(IQMLOP),DUMMY,'INT',IERR)
         !        END DO
         WRITE(LUPRI,*) 'integral blocks not implemented'
         CALL LSQUIT('# of integral blocks not implemented',lupri)       
      END IF
   END IF
ENDIF

IF (ATOMBASIS) THEN
   IPOS = INDEX(TEMPLINE,'Bas')
   IF (IPOS .NE. 0) THEN
      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
         WRITE (LUPRI,'(2X,A40)') 'Incorrect input for choice of atomic basis set'
         WRITE (LUPRI,'(2X,A40)') 'Format is "Basis=? ? ?"'
         CALL LSQUIT('Incorrect input for choice of atomic basis set',lupri)
      ELSE
         IPOS3 = INDEX(TEMPLINE((IPOS+IPOS2):),' ')
         IF (IPOS3 .LT. 10) THEN
            WRITE (StringFormat,'(A2,I1,A1,1X)') '(A',IPOS3 - 1,')'
         ELSE
            WRITE (StringFormat,'(A2,I2,A1)') '(A',(IPOS3 - 1),')'
         ENDIF
         READ (TEMPLINE((IPOS + IPOS2):),StringFormat) AtomicBasisset
      ENDIF
   ELSE
      WRITE (LUPRI,*) 'ATOMBASIS selected, but no atomic basis&
           & set specified for one atom type'
      CALL LSQUIT( 'ATOMBASIS selected, but no atomic basis set &
           &specified for one atom type',lupri)
   ENDIF

   !Auxiliary basisset
   IPOS = INDEX(TEMPLINE,'Aux')
   IF (IPOS .NE. 0) THEN
      AUXBASIS=.TRUE.
      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 4)) THEN
         WRITE (LUPRI,*) 'Incorrect input for choice of auxiliary atomic basis set'
         WRITE (LUPRI,*) ' Format is "Aux=? ? ?"'
         CALL LSQUIT('Incorrect input for choice of auxiliary atomic basis set',lupri)
      ELSE
         IPOS3 = INDEX(TEMPLINE((IPOS+IPOS2):),' ')
         IF (IPOS3 .LT. 10) THEN
            WRITE (StringFormat,'(A2,I1,A1,1X)') '(A',IPOS3 - 1,')'
         ELSE
            WRITE (StringFormat,'(A2,I2,A1)') '(A',(IPOS3 - 1),')'
         ENDIF
         READ (TEMPLINE((IPOS + IPOS2):),StringFormat) AuxBasisset
      ENDIF
   ENDIF
ENDIF
END SUBROUTINE READ_LINE5

!> \brief determine if the basis is unique
!> \author T. Kjaergaard
!> \date 2010
!> \param BASISSETLIBRARY the info about the basisset required
!> \param BASISSET the basis that may or may not be unique
!> \param unique output determining if the basis is unique
!> \param nsets the number of different basissets in library
!>
!> if the basiset is not contained in basissetlibrary unique is set to
!> zero otherwise it is set to the index of the basisset in the basisset 
!> library 
!>
SUBROUTINE DETERMINE_UNIQUE_BASIS(BASISSETLIBRARY,BASISSET,unique,nsets)
implicit none
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
CHARACTER(len=80)         :: basisset
INTEGER                   :: unique,nsets,I

unique=0

IF(nsets/=0)THEN
   DO I=1,nsets 
      IF(BASISSETLIBRARY%BASISSETNAME(I) .EQ. basisset)THEN
         unique=I
      ENDIF
   ENDDO
ENDIF

END SUBROUTINE DETERMINE_UNIQUE_BASIS

!> \brief determine if the charge is unique
!> \author T. Kjaergaard
!> \date 2010
!> \param BASISSETLIBRARY the info about the basisset required
!> \param basisnumber the number in the basislibrary
!> \param atomicCharge that may or may not be unique
!> \param uniqueCharge 0 if unique
!>
!> if the charge is unique in the basissetlibrary
!> uniqueCharge is set to zero otherwise it is set to the index
!>
SUBROUTINE DETERMINE_UNIQUE_CHARGE(BASISSETLIBRARY,basisnumber,&
     &atomicCharge,uniqueCharge)
implicit none
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
INTEGER               :: basisnumber,uniqueCharge,I
real(realk)           :: atomicCharge
LOGICAL               :: FOUND
uniqueCharge=0
FOUND=.FALSE.
DO I=1,BASISSETLIBRARY%nCharges(basisnumber)
  IF(BASISSETLIBRARY%Charges(basisnumber,I) .EQ. INT(atomicCharge) )THEN
     FOUND=.TRUE.
  ENDIF
ENDDO
IF(.NOT. FOUND) uniqueCharge=BASISSETLIBRARY%nCharges(basisnumber)+1

END SUBROUTINE DETERMINE_UNIQUE_CHARGE

!> \brief ATTACH_BASISINDEX AND BASISLABEL to the atom
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param moleclue the molecule structure
!> \param I the atom index in the molecule structure
!> \param basissetnumber1 the basisindex
!> \param basissetnumber1 the auxilliary basisindex
!> \param unique1 basis index 
!> \param unique2 auxilliary basis index
!> \param basis if BASIS keyword is given in line 1
!> \param atombasis if ATOMBASIS keyword is given in line 1
!> \param auxbasis is the auxilliary basis given in 1. line
!>
!> IF YOU ADD MORE BASISSET TO MOLECULE FILE YOU NEED TO 
!> INCREASE THE SIZE IN THE ATOM DERIVED TYPE IN TYPE-DEF.f90. 
!>
SUBROUTINE ATTACH_BASISINDEX_AND_BASISLABEL(LUPRI,MOLECULE,&
     &I,basissetnumber1,basissetnumber2,unique1,unique2,&
     &BASIS,ATOMBASIS,AUXBASIS)
implicit none
TYPE(MOLECULEINFO) :: MOLECULE 
LOGICAL            :: BASIS,ATOMBASIS,AUXBASIS
INTEGER            :: unique1,unique2,I,LUPRI !I=atomnumber
INTEGER            :: basissetnumber1,basissetnumber2


IF(BASIS)THEN
   MOLECULE%ATOM(I)%basisindex(1)=1
   MOLECULE%ATOM(I)%basislabel(1)='REGULAR  '
   MOLECULE%ATOM(I)%nbasis=1
   IF(AUXBASIS)THEN
      MOLECULE%ATOM(I)%basisindex(2)=2
      MOLECULE%ATOM(I)%basislabel(2)='AUXILIARY'
      MOLECULE%ATOM(I)%nbasis=2
   ENDIF
ELSEIF(ATOMBASIS)THEN
   IF(unique1 /= 0)THEN
      MOLECULE%ATOM(I)%basisindex(1)=unique1
      MOLECULE%ATOM(I)%basislabel(1)='REGULAR  '
   ELSE
      MOLECULE%ATOM(I)%basisindex(1)=basissetnumber1
      MOLECULE%ATOM(I)%basislabel(1)='REGULAR  '
   ENDIF
   MOLECULE%ATOM(I)%nbasis=1
   IF(AUXBASIS)THEN
      IF(unique2 /= 0)THEN
         MOLECULE%ATOM(I)%basisindex(2)=unique2
         MOLECULE%ATOM(I)%basislabel(2)='AUXILIARY'
      ELSE
         MOLECULE%ATOM(I)%basisindex(2)=basissetnumber2
         MOLECULE%ATOM(I)%basislabel(2)='AUXILIARY'
      ENDIF
   MOLECULE%ATOM(I)%nbasis=2
   ENDIF
ELSE
   CALL LSQUIT('Something wrong in determining basisetsets in READ_GEOMETRY',-1)
ENDIF

END SUBROUTINE ATTACH_BASISINDEX_AND_BASISLABEL

!> \brief build distinct atoms 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param ndatoms the number of distinct atoms
!> \param natoms the total number of atoms
!> \param natomlist for a given distinct atom a atomindex i full molecule
!> \param uatomtype for a given distinct atom the atomtype in full molecule
!> \param ls contains molecule and basis structure
!> \param iprint the printlevel, determining how much output should be generate
SUBROUTINE BUILD_DISTINCT_ATOMS(LUPRI,NDATOMS,NATOMS,NATOMLIST,UATOMTYPE,LS,IPRINT)
use BUILDBASISSET
use molecule_module
IMPLICIT NONE

INTEGER  :: LUPRI,NDATOMS,NATOMS,IPRINT,I
TYPE(LSITEM)        :: LS
INTEGER             :: CHARGE(NDATOMS),NATOMLIST(NDATOMS),UATOMTYPE(NDATOMS)
INTEGER             :: NDATOMS2
Character(len=80)   :: NAME(NDATOMS)
INTEGER             :: itype,TYPEN(NDATOMS),icharge
INTEGER             :: R
INTEGER,pointer :: NEWATOM(:)

R = LS%INPUT%BASIS%REGULAR%Labelindex
ALLOCATE(NEWATOM(LS%INPUT%BASIS%REGULAR%natomtypes))
NEWATOM=0
NDATOMS2 = 0
IF(R.EQ.0)THEN
   DO I=1,LS%INPUT%MOLECULE%nATOMs
      icharge = INT(LS%INPUT%MOLECULE%ATOM(I)%Charge)
      itype = LS%INPUT%BASIS%REGULAR%Chargeindex(icharge)
      IF(NEWATOM(itype).EQ.0)THEN
         NEWATOM(itype)=1

         NDATOMS2 = NDATOMS2+1
         TYPEN(NDATOMS2) = itype
         CHARGE(NDATOMS2) = icharge
         NAME(NDATOMS2)= LS%INPUT%BASIS%REGULAR%ATOMTYPE(itype)%NAME
         NATOMLIST(NDATOMS2) = I
         UATOMTYPE(NDATOMS2) = itype
!         UATOMLIST(I)=NDATOMS2
      ENDIF
   ENDDO
ELSE
   DO I=1,LS%INPUT%MOLECULE%nATOMs
      icharge = INT(LS%INPUT%MOLECULE%ATOM(I)%Charge)
      itype = LS%INPUT%MOLECULE%ATOM(I)%IDtype(R)
      IF(NEWATOM(itype).EQ.0)THEN
         NEWATOM(itype)=1

         NDATOMS2 = NDATOMS2+1
         TYPEN(NDATOMS2) = itype
         CHARGE(NDATOMS2) = icharge
         NAME(NDATOMS2)= LS%INPUT%BASIS%REGULAR%ATOMTYPE(itype)%NAME
         NATOMLIST(NDATOMS2) = I
         UATOMTYPE(NDATOMS2) = itype
!         UATOMLIST(I)=NDATOMS2        
      ENDIF
   ENDDO
ENDIF

IF(NDATOMS .NE. NDATOMS2)CALL LSQUIT('DISTINCT_ATOMS and BUILD_DISTINCT_ATOMS not consistent',lupri)

!DO I=1,NDATOMS
!   NATOMLIST(I)=ATOMINDEX(I)
!ENDDO
DEALLOCATE(NEWATOM)

END SUBROUTINE BUILD_DISTINCT_ATOMS
!!$FUNCTION ISOMAS(CHARGE,MASSNumber)
!!$!  Function to switch from mass number to isotope number sorted
!!$!  according to abundance
!!$implicit none
!!$INTEGER    ::  IORD,I,QMASS,MASSNumber
!!$IORD = 0
!!$DO I = 1, 5
!!$  QMASS = DISOTP(ICHARG,I,'MASS')
!!$  IF (ANINT(QMASS) .EQ. MASSNM) IORD = I
!!$END DO
!!$IF (IORD .EQ. 0) THEN
!!$  WRITE (LUPRI,'(/A,I4,A,I4)') 'ERROR: unknown mass',MASSNM,
!!$  &' for atom with charge ',ICHARG
!!$  CALL QUIT('Unknown mass for chosen atomic charge')
!!$ELSE
!!$  ISOMAS = IORD
!!$END IF
!!$END FUNCTION ISOMAS
!!$
!!$FUNCTION DISOTP(IATOM,ISOTOP,TYPE)
!!$! NOTE: Isotopes are sorted according to abundance,
!!$! i.e. DISOTP(IATOM,1,TYPE) will return the most abundant
!!$! isotope etc.
!!$! Proton mass and electron charge: 1986 CODATA Recommended Values
!!$!
!!$! Nuclear masses:
!!$!   A. H. Wapstra and K. Bos,
!!$!   Atomic Data and Nuclear Tables 19 (1977) 177
!!$!
!!$! Abundancies:
!!$!   Handbook of Chemistry and Physics, 73rd Edition
!!$!        
!!$! Nuclear moments and spins:
!!$!   P. Raghavan,
!!$!   Atomic Data and Nuclear Data Tables 42 (1989) 189
!!$!
!!$! Quadrupole moments:
!!$!   P.Pykkoe and J.Li
!!$!   Report HUKI 1-92
!!$!   Updated  to P.Pyykkoe, Mol.Phys. (2001) by K.Ruud, Aug.2001
!!$!
!!$!   Nuclear masses, Abundancies, nuclear moments, spins 
!!$!   and quadrupole moments for Z= 55 to Z = 86:
!!$!
!!$!   I. Mills, T. Cvitas, K. Homann, N. Kallay, and K. Kuchitsu
!!$!   Quantities, Units and Symbols in Physical Chemistry
!!$!   (IUPAC, Blackwell Scientific, Oxford, 1988)
!!$!
!!$CHARACTER*(*) TYPE
!!$INTEGER           :: MAXISO,MAXCHR
!!$real(realk)       :: DATNUC(5,6,86),PMASS,EMASS,XFAMU,THRESH,DMP
!!$real(realk)       :: FOUR,ZERO
!!$PMASS = 1.007276470D0
!!$EMASS = 9.10938188D-31
!!$XFAMU = 1822.88848D0
!!$THRESH = 1.0D-10
!!$DMP = PMASS*XFAMU*EMASS
!!$MAXISO = 6
!!$MAXCHR = 86
!!$FOUR = 4.0D0
!!$ZERO = 0.0D0

!!$! H - Ne
!!$! ======
!!$     DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=0,10) /
!!$! Dummy:
!!$     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
!!$     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
!!$     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
!!$     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
!!$     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
!!$     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
!!$! H:
!!$     &  1.007825D0, 99.985000D0,   .500000D0,  2.792847D0,   .000000D0,
!!$     &  2.014102D0,   .015000D0,  1.000000D0,   .857438D0,   .002860D0,
!!$     &  3.016049D0,   .000000D0,   .500000D0,  2.978962D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! He:
!!$     &  4.002603D0, 99.999870D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &  3.016029D0,   .000130D0,   .500000D0, -2.127625D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! Li:
!!$     &  7.016005D0, 92.500000D0,  1.500000D0,  3.256427D0,  -.040100D0,
!!$     &  6.015123D0,  7.500000D0,  1.000000D0,   .822047D0,  -.000808D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! Be:
!!$     &  9.012183D0,100.000000D0,  1.500000D0, -1.177800D0,   .052880D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! B:
!!$     & 11.009305D0, 80.100000D0,  1.500000D0,  2.688649D0,   .040590D0,
!!$     & 10.012938D0, 19.900000D0,  3.000000D0,  1.800645D0,   .084590D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$!  C:
!!$     & 12.000000D0, 98.900000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 13.003355D0,  1.100000D0,   .500000D0,   .702412D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! N:
!!$     & 14.003074D0, 99.630000D0,  1.000000D0,   .403761D0,   .020440D0,
!!$     & 15.000109D0,   .370000D0,   .500000D0,  -.283189D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 12.000000D0, 98.900000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! O:
!!$     & 15.994915D0, 99.760000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 17.999159D0,   .200000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 16.999131D0,   .040000D0,  2.500000D0, -1.893790D0,  -.025580D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 14.003074D0, 99.630000D0,  1.000000D0,   .403761D0,   .020200D0,
!!$! F:
!!$     & 18.998403D0,100.000000D0,   .500000D0,  2.628868D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 15.994915D0, 99.760000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! Ne:
!!$     & 19.992439D0, 90.480000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 21.991384D0,  9.250000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 20.993845D0,   .270000D0,  1.500000D0,  -.661797D0,   .101550D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0/
!!$! Na - Ar
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=11,18) /
!!$! Na:
!!$     & 22.989770D0,100.000000D0,  1.500000D0,  2.217656D0,   .104000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! Mg:
!!$     & 23.985045D0, 78.990000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 25.982595D0, 11.010000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 24.985839D0, 10.000000D0,  2.500000D0,  -.855450D0,   .199400D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! Al:
!!$     & 26.981541D0,100.000000D0,  2.500000D0,  3.641507D0,   .146600D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! Si:
!!$     & 27.976928D0, 92.230000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 28.976496D0,  4.670000D0,   .500000D0,  -.555290D0,   .000000D0,
!!$     & 29.973772D0,  3.100000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! P:
!!$     & 30.973763D0,100.000000D0,   .500000D0,  1.131600D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! S:
!!$     & 31.972072D0, 95.020000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 33.967868D0,  4.210000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 32.971459D0,   .750000D0,  1.500000D0,   .643821D0,  -.067800D0,
!!$     & 35.967079D0,   .020000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! Cl:
!!$     & 34.968853D0, 75.770000D0,  1.500000D0,   .821874D0,  -.081650D0,
!!$     & 36.965903D0, 24.230000D0,  1.500000D0,   .684124D0,  -.064350D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! Ar:
!!$     & 39.962383D0, 99.600000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 35.967546D0,   .337000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 37.962732D0,   .063000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0/
!!$! K - Ca
!!$! ======
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=19,20) /
!!$! K:
!!$     & 38.963708D0, 93.258100D0,  1.500000D0,   .391507D0,   .058500D0,
!!$     & 40.961825D0,  6.730200D0,  1.500000D0,   .214893D0,   .071100D0,
!!$     & 39.963999D0,   .011700D0,  4.000000D0, -1.298100D0,  -.073000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$! Ca:
!!$     & 39.962591D0, 96.941000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 43.955485D0,  2.086000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 41.958622D0,   .647000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 47.952532D0,   .187000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 42.958770D0,   .135000D0,  3.500000D0, -1.317643D0,  -.040800D0,
!!$     & 45.953689D0,   .004000D0,   .000000D0,   .000000D0,   .000000D0/
!!$! Sc - Zn
!!$! =======
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=21,30) /
!!$! Sc:
!!$     & 44.955914D0,100.000000D0,  3.500000D0,  4.756487D0,  -.220000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Ti:
!!$C
!!$     & 47.947947D0, 73.800000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 45.952633D0,  8.000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 46.951765D0,  7.300000D0,  2.500000D0,  -.788480D0,   .302000D0,
!!$     & 48.947871D0,  5.500000D0,  3.500000D0, -1.104170D0,   .247000D0,
!!$     & 49.944786D0,  5.400000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     V:
!!$C
!!$     & 50.943963D0, 99.750000D0,  3.500000D0,  5.148706D0,  -.052000D0,
!!$     & 49.947161D0,   .250000D0,  6.000000D0,  3.345689D0,   .210000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Cr:
!!$C
!!$     & 51.940510D0, 83.790000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 52.940651D0,  9.500000D0,  1.500000D0,  -.474540D0,  -.150000D0,
!!$     & 49.946046D0,  4.345000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 53.938882D0,  2.365000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Mn:
!!$C
!!$     & 54.938046D0,100.000000D0,  2.500000D0,  3.468719D0,   .330000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Fe:
!!$C
!!$     & 55.934939D0, 91.720000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 53.939612D0,  5.900000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 56.935396D0,  2.100000D0,   .500000D0,   .090623D0,   .000000D0,
!!$     & 57.933278D0,   .280000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Co:
!!$C
!!$     & 58.933198D0,100.000000D0,  3.500000D0,  4.627000D0,   .420000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Ni:
!!$C
!!$     & 57.935347D0, 68.077000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 59.930789D0, 26.223000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 61.928346D0,  3.634000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 60.931059D0,  1.140000D0,  1.500000D0,  -.750020D0,   .162000D0,
!!$     & 63.927968D0,  0.926000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Cu:
!!$C
!!$     & 62.929599D0, 69.170000D0,  1.500000D0,  2.227206D0,  -.220000D0,
!!$     & 64.927792D0, 30.830000D0,  1.500000D0,  2.381610D0,  -.204000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Zn:
!!$C
!!$     & 63.929145D0, 48.600000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 65.926035D0, 27.900000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 67.924846D0, 18.800000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 66.927129D0,  4.100000D0,  2.500000D0,   .875479D0,   .150000D0,
!!$     & 69.925325D0,   .600000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0/
!!$C
!!$C     Ga - Kr
!!$C     =======
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=31,36) /
!!$C
!!$C     Ga:
!!$C
!!$     & 68.925581D0, 60.108000D0,  1.500000D0,  2.016589D0,   .171000D0,
!!$     & 70.924701D0, 39.892000D0,  1.500000D0,  2.562266D0,   .107000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Ge:
!!$C
!!$     & 73.921179D0, 35.940000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 71.922080D0, 27.660000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 69.924250D0, 21.240000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 72.923464D0,  7.720000D0,  4.500000D0,  -.879468D0,  -.196000D0,
!!$     & 75.921403D0,  7.440000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     As:
!!$C
!!$     & 74.921596D0,100.000000D0,  1.500000D0,  1.439475D0,   .314000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Se:
!!$C
!!$     & 79.916521D0, 49.610000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 77.917304D0, 23.770000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 75.919207D0,  9.360000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 81.916709D0,  8.740000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 76.919908D0,  7.630000D0,   .500000D0,   .535042D0,   .000000D0,
!!$     & 73.922477D0,  0.890000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Br:
!!$C
!!$     & 78.918336D0, 50.690000D0,  1.500000D0,  2.106400D0,   .313000D0,
!!$     & 80.916290D0, 49.310000D0,  1.500000D0,  2.270562D0,   .261500D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Kr:
!!$C
!!$     & 83.911506D0, 57.000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 85.910614D0, 17.300000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 81.913483D0, 11.600000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 82.914134D0, 11.500000D0,  4.500000D0,  -.970669D0,   .259000D0,
!!$     & 79.916375D0,  2.250000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 77.920397D0,  0.350000D0,   .000000D0,   .000000D0,   .000000D0/
!!$C
!!$C     Rb:
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=37,45) /
!!$C
!!$     & 84.911800D0, 72.170000D0,  2.500000D0,  1.353352D0,   .276000D0,
!!$     & 86.909184D0, 27.830000D0,  1.500000D0,  2.751818D0,   .133500D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Sr:
!!$C
!!$     & 87.905625D0, 82.580000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 85.909273D0,  9.860000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 86.908890D0,  7.000000D0,  4.500000D0, -1.093603D0,   .335000D0,
!!$     & 83.913428D0,   .560000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Y:
!!$C
!!$     & 88.905856D0,100.000000D0,   .500000D0,  -.137415D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Zr:
!!$C
!!$     & 89.904708D0, 51.450000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 93.906319D0, 17.380000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 91.905039D0, 17.150000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 90.905644D0, 11.220000D0,  2.500000D0, -1.303620D0,  -.176000D0,
!!$     & 95.908272D0,  2.800000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Nb:
!!$C
!!$     & 92.906378D0,100.000000D0,  4.500000D0,  6.170500D0,  -.320000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Mo:
!!$C
!!$     & 97.905405D0, 24.130000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 95.904676D0, 16.680000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 94.905838D0, 15.920000D0,  2.500000D0,  -.914200D0,  -.022000D0,
!!$     & 93.905086D0, 14.840000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 99.907473D0,  9.630000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 96.906018D0,  9.550000D0,  2.500000D0,  -.933500D0,  0.255000D0,
!!$C
!!$C     Tc:
!!$C
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Ru:
!!$C
!!$     &101.904348D0, 31.600000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &103.905422D0, 18.600000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &100.905581D0, 17.100000D0,  2.500000D0,  -.718800D0,   .457000D0,
!!$     & 98.905937D0, 12.700000D0,  2.500000D0,  -.641300D0,   .079000D0,
!!$     & 99.904218D0, 12.600000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     & 95.907596D0,  5.540000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Rh:
!!$C
!!$     &102.905503D0,100.000000D0,   .500000D0,  -.088400D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0/
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=46,54) /
!!$C
!!$C     Pd:
!!$C
!!$     &105.903475D0, 27.330000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &107.903894D0, 26.460000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &104.905075D0, 22.330000D0,  2.500000D0,  -.642000D0,   .660000D0,
!!$     &109.905169D0, 11.720000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &103.904026D0, 11.140000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &101.905609D0,  1.020000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Ag:
!!$C
!!$     &106.905095D0, 51.839000D0,   .500000D0,  -.113570D0,   .000000D0,
!!$     &108.904754D0, 48.161000D0,   .500000D0,   .130563D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Cd:
!!$C
!!$     &113.903361D0, 28.730000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &111.902761D0, 24.130000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &110.904182D0, 12.800000D0,   .500000D0,  -.594886D0,   .000000D0,
!!$     &109.903007D0, 12.490000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &112.904401D0, 12.220000D0,   .500000D0,  -.622301D0,   .000000D0,
!!$     &115.904758D0,  7.490000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     In:
!!$C
!!$     &114.903875D0, 95.700000D0,  4.500000D0,  5.540800D0,   .810000D0,
!!$     &112.904056D0,  4.300000D0,  4.500000D0,  5.528900D0,   .799000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Sn:
!!$C
!!$     &119.902199D0, 32.590000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &117.901607D0, 24.220000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &115.901744D0, 14.530000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &118.903310D0,  8.580000D0,   .500000D0, -1.047280D0,   .000000D0,
!!$     &116.902954D0,  7.680000D0,   .500000D0, -1.001040D0,   .000000D0,
!!$     &123.905271D0,  5.790000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Sb:
!!$C
!!$     &120.903824D0, 57.360000D0,  2.500000D0,  3.363400D0,  -.360000D0,
!!$     &122.904222D0, 42.640000D0,  3.500000D0,  2.549800D0,  -.490000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Te:
!!$C
!!$     &129.906229D0, 33.870000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &127.904464D0, 31.700000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &125.903310D0, 18.930000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &124.904435D0,  7.120000D0,   .500000D0,  -.888505D0,   .000000D0,
!!$     &123.902825D0,  4.790000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &121.903055D0,  2.590000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     I:
!!$C
!!$     &126.904477D0,100.000000D0,  2.500000D0,  2.813273D0,  -.710000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
!!$C
!!$C     Xe:
!!$C
!!$     &131.904148D0, 26.900000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &128.904780D0, 26.400000D0,   .500000D0,  -.777976D0,   .000000D0,
!!$     &130.905076D0, 21.200000D0,  1.500000D0,   .691862D0,  -.114000D0,
!!$     &133.905395D0, 10.400000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &135.907219D0,  8.900000D0,   .000000D0,   .000000D0,   .000000D0,
!!$     &129.903510D0,  4.100000D0,   .000000D0,   .000000D0,   .000000D0/
!!$C
!!$C
!!$C
!!$C
!!$C    Cs - Rn*
!!$C    =======
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=55,64) /
!!$C
!!$C     Cs:
!!$C
!!$     &132.905429D0,100.000000D0,  3.50000D0,   2.582025D0,  -0.00343D0,
!!$     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
!!$     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
!!$     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
!!$     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
!!$     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
!!$C
!!$C     Ba:
!!$C   
!!$     &137.905232D0, 71.70000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &136.905812D0, 11.23000D0,   1.50000D0,  0.937365D0,   0.245000D0,
!!$     &135.904553D0,  7.85400D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &134.905665D0,  6.59200D0,   1.50000D0,  0.837943D0,   0.160000D0,
!!$     &133.904486D0,  2.41700D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &131.905042D0,  0.10100D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$C
!!$C     La:
!!$C
!!$     &138.906347D0, 99.90980D0,   3.50000D0,  2.7830455D0,  0.200000D0,
!!$     &137.907105D0,  0.09020D0,   5.00000D0,  3.7136460D0,  0.450000D0,
!!$     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$C
!!$C     Ce
!!$C
!!$     &139.905433D0, 88.48000D0,   0.00000D0,  0.000000D0,    0.00000D0,
!!$     &141.909241D0, 11.08000D0,   0.00000D0,  0.000000D0,    0.00000D0,
!!$     &137.905985D0,  0.25000D0,   0.00000D0,  0.000000D0,    0.00000D0,
!!$     &135.907140D0,  0.19000D0,   0.00000D0,  0.000000D0,    0.00000D0,
!!$     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$C
!!$C     Pr:
!!$C
!!$     &140.907647D0,100.0000D0,    2.50000D0,  4.275400D0,  -0.058900D0,
!!$     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$C
!!$C     Nd:
!!$C
!!$     &141.907719D0, 27.130000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &143.910083D0, 23.800000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &145.913113D0, 17.190000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$     &142.909810D0, 12.180000D0,   3.50000D0, -0.065000D0,  -0.630000D0,
!!$     &144.912570D0,  8.300000D0,   3.50000D0, -1.065000D0,  -0.330000D0,
!!$     &147.916889D0,  5.760000D0,   0.00000D0,  0.000000D0,   0.000000D0,
!!$C
!!$C     Pm:
!!$C
!!$     &144.912743D0,100.000000D0,   2.50000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$C
!!$C     Sm:
!!$C
!!$     &151.919728D0, 26.700000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &153.922205D0, 22.700000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &146.914894D0, 15.000000D0,   3.500000D0,-0.814800D0,  -0.259000D0,
!!$     &148.917180D0, 13.800000D0,   3.500000D0,-0.671700D0,   0.075000D0,
!!$     &147.914819D0, 11.300000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &149.917273D0,  7.400000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Eu:
!!$C
!!$     &152.921225D0, 52.200000D0,   2.500000D0, 1.533000D0,   2.412000D0,
!!$     &150.919702D0, 47.800000D0,   2.500000D0, 3.471700D0,   0.903000D0,
!!$     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
!!$C
!!$C     Gd:
!!$C
!!$     &157.924019D0, 24.840000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &159.927049D0, 21.860000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &155.922118D0, 20.470000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &156.923956D0, 15.650000D0,   1.500000D0,-0.337260D0,   1.350000D0,
!!$     &154.922618D0, 14.800000D0,   1.500000D0,-0.257230D0,   1.270000D0,
!!$     &153.920861D0,  2.180000D0,   0.000000D0, 0.000000D0,   0.000000D0/
!!$C
!!$C     Tb:
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=65,74) /
!!$C
!!$     &158.925342D0,100.000000D0,   1.500000D0, 2.014000D0,   1.432000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Dy:
!!$C
!!$     &163.929171D0, 28.200000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &161.926795D0, 25.500000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &162.928728D0, 24.900000D0,   2.500000D0, 0.672600D0,   2.648000D0,
!!$     &160.926930D0, 18.900000D0,   2.500000D0,-0.480300D0,   2.507000D0,
!!$     &159.925193D0,  2.340000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &157.924277D0,  0.100000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Ho:
!!$C
!!$     &164.930319D0,100.000000D0,   3.500000D0, 4.173000D0,   3.580000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Er:
!!$C
!!$     &165.930290D0, 33.600000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &167.932368D0, 26.800000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &166.932368D0, 22.950000D0,   3.500000D0,-0.563850D0,   3.565000D0,
!!$     &169.935461D0, 14.900000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &163.929198D0,  1.610000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &161.928775D0,  0.140000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Tm:
!!$C
!!$     &168.934212D0,100.000000D0,   0.500000D0,-0.231600D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Yb:
!!$C
!!$     &173.938859D0, 31.800000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &171.936378D0, 21.900000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &172.938208D0, 16.120000D0,   2.500000D0,-0.679890D0,   2.800000D0,
!!$     &170.936323D0, 14.300000D0,   0.500000D0, 0.493670D0,   0.000000D0,
!!$     &175.942564D0, 12.700000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &169.934759D0,  3.050000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Lu:
!!$C
!!$     &174.940770D0, 97.410000D0,   3.500000D0, 2.232700D0,   3.490000D0,
!!$     &175.942679D0,  2.590000D0,   7.000000D0, 3.169200D0,   4.970000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Hf:
!!$C
!!$     &179.9465457D0,35.100000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &177.943696D0, 27.297000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &176.943217D0, 18.606000D0,   3.500000D0, 0.793500D0,   3.365000D0,
!!$     &178.9458122D0,13.629000D0,   4.500000D0,-0.640900D0,   3.793000D0,
!!$     &175.941406D0,  5.206000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &173.940044D0,  0.162000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Ta:
!!$C
!!$     &180.947992D0, 99.988000D0,   3.500000D0, 2.370500D0,   3.170000D0,
!!$     &179.947462D0,  0.012000D0,   8.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     W:
!!$C
!!$     &183.950928D0, 30.670000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &185.954357D0, 28.600000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &181.948202D0, 26.300000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &182.950928D0, 14.300000D0,   0.500000D0, 0.11778476,   0.000000D0,
!!$     &179.947462D0,  0.162000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0/
!!$C
!!$C     Re:
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=75,80) /
!!$C
!!$     &186.955744D0, 62.600000D0,   2.500000D0, 3.219700D0,   2.070000D0,
!!$     &184.952951D0, 37.400000D0,   2.500000D0, 3.187100D0,   2.180000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Os:
!!$C
!!$     &191.961467D0, 41.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &189.958436D0, 26.400000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &188.958436D0, 16.100000D0,   1.500000D0, 0.659933D0,   0.856000D0,
!!$     &187.955830D0, 13.300000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &186.955741D0,  1.600000D0,   0.500000D0, 0.06465189,   0.000000D0,
!!$     &185.953830D0,  1.580000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Ir:
!!$C
!!$     &192.962917D0, 62.600000D0,   1.500000D0, 0.163700D0,   0.751000D0,
!!$     &190.960584D0, 37.400000D0,   1.500000D0, 0.150700D0,   0.816000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Pt:
!!$C
!!$     &194.964766D0, 33.800000D0,   0.500000D0, 0.609520D0,   0.000000D0,
!!$     &193.962655D0, 32.900000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &195.964926D0, 25.300000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &197.967869D0,  7.200000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &191.961019D0,  0.790000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &189.959917D0,  0.010000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Au:
!!$C   
!!$     &196.966543D0,100.000000D0,   1.500000D0, 0.148158D0,   0.547000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Hg:
!!$C
!!$     &201.970617D0, 29.860000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &199.968300D0, 23.100000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &198.968254D0, 16.870000D0,   0.500000D0, 0.50588549,   0.000000D0,
!!$     &200.970277D0, 13.180000D0,   1.500000D0,-0.5602257 ,   0.386000D0,
!!$     &197.966743D0,  9.970000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &203.973467D0,  6.870000D0,   0.000000D0, 0.000000D0,   0.000000D0/
!!$C
!!$C     Tl:
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=81,86) /
!!$C
!!$     &204.974401D0, 70.476000D0,   0.500000D0, 1.63831461D0, 0.000000D0,
!!$     &202.972320D0, 29.524000D0,   0.500000D0, 1.62225787D0, 0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Pb:
!!$C   
!!$     &207.976627D0, 52.400000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &205.975872D0, 24.100000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &206.975872D0, 22.100000D0,   0.500000D0, 0.582583D0,   0.000000D0,
!!$     &203.973020D0,  1.400000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Bi:
!!$C
!!$     &208.980374D0,100.000000D0,   4.500000D0, 4.110600D0,  -0.516000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Po:
!!$C
!!$     &208.982404D0,  0.000000D0,   0.500000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     At:
!!$C
!!$     &209.987126D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$C
!!$C     Rn:
!!$C
!!$     &222.017571D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
!!$     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0/

!!$IF (ISOTOP .GT. 6) THEN
!!$  WRITE (LUPRI,'(//,A,2(/,A,I5),A)')'ISOTOP too large in DISOTP. ',&
!!$  &'Input value: ',ISOTOP,' Maximum value:',6,' Program cannot continue.'
!!$  CALL QUIT('MAXISO exceeded in DISOTP')
!!$END IF
!!$IF (IATOM .GT. 86) THEN
!!$  WRITE (LUPRI,'(//,A,2(/,A,I5),A)')' IATOM too large in DISOTP. ',&
!!$  &'Input value: ',IATOM,' Maximum value:',86,' Program cannot continue.'
!!$  CALL QUIT('MAXCHR exceeded in DISOTP')
!!$END IF
!!$
!!$IF (IATOM .LE. 0) THEN
!!$!  This is a floating orbital, a point charge,
!!$!   or an auxiliary basis set /Mar 2004 hjaaj
!!$  DISOTP = D0
!!$ELSE IF (TYPE .EQ. 'MASS') THEN
!!$  DISOTP = DATNUC(1,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'A') THEN
!!$  DISOTP = NINT(DATNUC(1,ISOTOP,IATOM))
!!$ELSE IF (TYPE .EQ. 'ABUNDANCE') THEN
!!$  DISOTP = DATNUC(LUPRI,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'SPIN') THEN
!!$  DISOTP = DATNUC(3,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'MMOM') THEN
!!$  DISOTP = DATNUC(4,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'GVAL') THEN
!!$  SPIN = DATNUC(3,ISOTOP,IATOM)
!!$  IF (SPIN .GT. THRESH) THEN
!!$    DISOTP = DATNUC(4,ISOTOP,IATOM)/SPIN
!!$  ELSE
!!$    DISOTP = 0.0D0
!!$  END IF
!!$ELSE IF (TYPE .EQ. 'LARMOR') THEN
!!$  SPIN = DATNUC(3,ISOTOP,IATOM)
!!$  IF (SPIN .GT. THRESH) THEN
!!$    DISOTP = ABS(ECHARGE*DATNUC(4,ISOTOP,IATOM)/(D4*PI*SPIN*DMP))
!!$  ELSE
!!$    DISOTP = 0.0D0
!!$  END IF
!!$ELSE IF (TYPE .EQ. 'QMOM') THEN
!!$  DISOTP = DATNUC(5,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'NEUTRONS') THEN
!!$  DISOTP = FLOAT(NINT(DATNUC(1,ISOTOP,IATOM)-IATOM))
!!$ELSE
!!$  WRITE (LUPRI,'(//,3A,/,A)')'Keyword ',TYPE,' unknown in DISOTP. ',&
!!$  &'Program cannot continue.'
!!$  CALL QUIT('Illegal keyword in DISOTP')
!!$END IF
!!$
!!$END FUNCTION DISOTP

END MODULE
