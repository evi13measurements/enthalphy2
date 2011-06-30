!> @file
!> Module containing subroutines related to the molecule

!> Standard molecule module. Contains also the orbital information.
!> \author S. Reine
!> \date 2010-02-21
MODULE molecule_module
use precision
use typedef
use memory_handling
use molecule_type
CONTAINS


!*****************************************
!*
!*  BASISSETINFO INITIATION ROUTINES
!*
!*****************************************

!> \brief Returns number of atoms and regular and auxiliary orbitals of given molecule
!> \author S. Reine
!> \date 2010-02-21
!> \param molecule Contains the information about the molecule
!> \param nAtoms The number of atoms
!> \param nBastReg The number of regular orbitals
!> \param nBastAux The number of auxiliary orbitals
SUBROUTINE getMolecularDimensions(MOLECULE,nAtoms,nBastReg,nBastAux)
implicit none
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
Integer,intent(out)           :: nAtoms,nBastReg,nBastAux

nAtoms   = MOLECULE%nAtoms
nBastReg = MOLECULE%nBastREG
nBastAux = MOLECULE%nBastAUX
!
END SUBROUTINE getMolecularDimensions

!> \brief Sets up the orbital information for a given molecule
!> \author S. Reine
!> \date 2010-02-21
!> \param Molecule Contains the information about the molecule
!> \param orbitalInfo Contains orbital indeces for the different atoms
SUBROUTINE setMolecularOrbitalInfo(MOLECULE,orbitalInfo)
implicit none
TYPE(MOLECULEINFO),intent(in)            :: MOLECULE
TYPE(MOLECULARORBITALINFO),intent(inout) :: orbitalInfo
!
integer :: iAtom,iReg,iAux,nOrbReg,nOrbAux


CALL getMolecularDimensions(MOLECULE,orbitalInfo%nAtoms,orbitalInfo%nBastReg,&
   &                        orbitalInfo%nBastAux)

CALL initMolecularOrbitalInfo(orbitalInfo,orbitalInfo%nAtoms)

iReg = 0
iAux = 0
DO iAtom=1,orbitalInfo%nAtoms
   nOrbReg = MOLECULE%ATOM(iAtom)%nContOrbREG
   nOrbAux = MOLECULE%ATOM(iAtom)%nContOrbAUX
   orbitalInfo%numAtomicOrbitalsReg(iAtom)=nOrbReg
   orbitalInfo%numAtomicOrbitalsAux(iAtom)=nOrbAux
   orbitalInfo%startAtomicOrbitalsReg(iAtom) = iReg+1
   orbitalInfo%startAtomicOrbitalsAux(iAtom) = iAux+1
   iReg = iReg + nOrbReg
   iAux = iAux + nOrbAux
   orbitalInfo%endAtomicOrbitalsReg(iAtom) = iReg
   orbitalInfo%endAtomicOrbitalsAux(iAtom) = iAux
ENDDO
END SUBROUTINE setMolecularOrbitalInfo

!> \brief Returns the orbital information of a given atom
!> \author S. Reine
!> \date 2010-02-19
!> \param orbitalInfo Contains the orbital-information of a given molecule
!> \param iAtom Specifies the atomic number in question
!> \param nReg The number of regular basis functions on given atom
!> \param startReg The starting orbital index of given atom
!> \param endReg The last orbital index of given atom
!> \param nAux The number of auxiliary basis functions on given atom
!> \param startAux The starting auxiliary orbital index of given atom
!> \param endAux The last auxiliary orbital index of given atom
SUBROUTINE getAtomicOrbitalInfo(orbitalInfo,iAtom,nReg,startReg,endReg,nAux,startAux,endAux)
use typedef
implicit none
TYPE(MOLECULARORBITALINFO),intent(IN) :: orbitalInfo
Integer,intent(IN)  :: iAtom
Integer,intent(OUT) :: nReg,startReg,endReg,nAux,startAux,endAux
!
nReg     = orbitalInfo%numAtomicOrbitalsReg(iAtom)
startReg = orbitalInfo%startAtomicOrbitalsReg(iAtom)
endReg   = orbitalInfo%endAtomicOrbitalsReg(iAtom)
nAux     = orbitalInfo%numAtomicOrbitalsAux(iAtom)
startAux = orbitalInfo%startAtomicOrbitalsAux(iAtom)
endAux   = orbitalInfo%endAtomicOrbitalsAux(iAtom)
!
END SUBROUTINE getAtomicOrbitalInfo

!> \brief Initialize orbitalInfo type
!> \author S. Reine
!> \date 2010-02-21
!> \param orbitalInfo Contains orbital indeces for the different atoms
!> \param nAtoms The number of atoms
SUBROUTINE initMolecularOrbitalInfo(orbitalInfo,nAtoms)
implicit none
TYPE(MOLECULARORBITALINFO),intent(out) :: orbitalInfo
integer,intent(in)                     :: nAtoms

CALL mem_alloc(orbitalInfo%numAtomicOrbitalsReg,nAtoms)
CALL mem_alloc(orbitalInfo%startAtomicOrbitalsReg,nAtoms)
CALL mem_alloc(orbitalInfo%endAtomicOrbitalsReg,nAtoms)
CALL mem_alloc(orbitalInfo%numAtomicOrbitalsAux,nAtoms)
CALL mem_alloc(orbitalInfo%startAtomicOrbitalsAux,nAtoms)
CALL mem_alloc(orbitalInfo%endAtomicOrbitalsAux,nAtoms)

END SUBROUTINE initMolecularOrbitalInfo

 
!> \brief Frees orbitalInfo type
!> \author S. Reine
!> \date 2010-02-21
!> \param orbitalInfo Contains orbital indeces for the different atoms
SUBROUTINE freeMolecularOrbitalInfo(orbitalInfo)
implicit none
TYPE(MOLECULARORBITALINFO),INTENT(INOUT) :: orbitalInfo
!
CALL mem_dealloc(orbitalInfo%numAtomicOrbitalsReg)
CALL mem_dealloc(orbitalInfo%startAtomicOrbitalsReg)
CALL mem_dealloc(orbitalInfo%endAtomicOrbitalsReg)
CALL mem_dealloc(orbitalInfo%numAtomicOrbitalsAux)
CALL mem_dealloc(orbitalInfo%startAtomicOrbitalsAux)
CALL mem_dealloc(orbitalInfo%endAtomicOrbitalsAux)
END SUBROUTINE freeMolecularOrbitalInfo

!> \brief Determined the number of electrons for a given molecule
!> \author T. Kjaergaard
!> \date 2010-02-21
!> \param Molecule Contains the information about the molecule
!> \param Moleculecharge The charge of the molecule
!> \param nelectrons The number of electrons
SUBROUTINE DETERMINE_NELECTRONS(Molecule,Moleculecharge,nelectrons)
implicit none
TYPE(MOLECULEINFO),intent(IN) :: MOLECULE
real(realk),intent(IN)        :: Moleculecharge
integer,intent(OUT)           :: Nelectrons
!
integer :: I,NCHARGE

NCHARGE = 0
DO I = 1,MOLECULE%NATOMS
   NCHARGE = INT(MOLECULE%ATOM(I)%CHARGE)+NCHARGE
ENDDO

NELECTRONS = NCHARGE - INT(Moleculecharge)

END SUBROUTINE DETERMINE_NELECTRONS

!> \brief Divide a molecule into molecular fragments (by setting up indices)
!> \author S. Reine
!> \date 2010-02-05
!> \param Molecule The molecule to be fragmented
!> \param fragmentIndex Indices specifying for each atom in Molecule which fragment it belongs to
!> \param numFragments The number of fragments the molecule should be divied into
!> \param lupri Default output unit
SUBROUTINE fragmentMolecule(Molecule,fragmentIndex,numFragments,lupri)
implicit none
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
Integer,intent(in)            :: numFragments,lupri
Integer,intent(out)           :: fragmentIndex(MOLECULE%nAtoms)
!
Integer :: numOrbitals,numFragOrbitals,iFragment,I,totOrb
logical :: Increased


IF (numFragments.GT.MOLECULE%nAtoms) THEN
  CALL LSQUIT('ERROR: fragmentMolecule entered with numFragments > nAtoms',lupri)
ELSEIF (numFragments.EQ.MOLECULE%nAtoms) THEN
  TOTorb=0
  iFragment = 0
  DO I=1,MOLECULE%nAtoms
    numOrbitals = MOLECULE%ATOM(I)%nContOrbREG
    TOTorb = TOTorb + MOLECULE%ATOM(I)%nContOrbREG
    iFragment = iFragment + 1
    fragmentIndex(I) = iFragment
  ENDDO
ELSE
! Divide the molecule into fragments with approximately similar number of
! orbitals

! First get the average number of orbitals per fragment
  numFragOrbitals = MOLECULE%nbastREG/numFragments

! Then partition the molecule into fragments of about this number of orbitails
  TOTorb=0
  numOrbitals = 0
  iFragment   = 1
  Increased = .FALSE.
  DO I=1,MOLECULE%nAtoms
    numOrbitals = numOrbitals + MOLECULE%ATOM(I)%nContOrbREG
    TOTorb = TOTorb + MOLECULE%ATOM(I)%nContOrbREG
    fragmentIndex(I) = iFragment
    Increased = .TRUE.
    IF((TOTorb .GE. iFragment*numFragOrbitals .AND. .NOT. (ifragment.EQ.numFragments) ))THEN
      iFragment = iFragment + 1
      numOrbitals = 0
      Increased = .FALSE.
    ENDIF
  ENDDO
  IF(.NOT.Increased) ifragment=ifragment-1
  IF(iFragment .NE. numFragments) THEN
     TOTorb=0
     numOrbitals = 0
     iFragment   = numFragments
     Increased = .FALSE.
     DO I=MOLECULE%nAtoms,1,-1
        numOrbitals = numOrbitals + MOLECULE%ATOM(I)%nContOrbREG
        TOTorb = TOTorb + MOLECULE%ATOM(I)%nContOrbREG
        fragmentIndex(I) = iFragment
        Increased = .TRUE.
        IF (TOTorb .GE. (numFragments-iFragment+1)*numFragOrbitals .AND. .NOT. ((numFragments-ifragment+1).EQ.numFragments)) THEN
           iFragment = iFragment - 1
           numOrbitals = 0
           Increased = .FALSE.
        ENDIF
     ENDDO
     IF(.NOT.Increased)ifragment=ifragment+1
     ifragment = numFragments-iFragment+1
  ENDIF
  IF(iFragment .NE. numFragments) THEN
     WRITE(LUPRI,*)'FRAGEMENT ',iFragment
     WRITE(LUPRI,*)'NODES     ',numFragments
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     WRITE(LUPRI,*)'WARNING WARNING WANING '
     CALL LSQUIT('ifrag not equal to number of nodes',lupri)
  ENDIF
ENDIF

END SUBROUTINE fragmentMolecule

!> \brief Builds a molecular fragment from a subset of the atoms in the original molecule
!> \author S. Reine and T. Kjaergaard
!> \date 2010-02-21
!> \param DALMOL The original molecule
!> \param FRAGMOL The fragment molecule
!> \param FRAGBASIS the basisinfo 
!> \param AUXBASIS logical = true if AUXBASIS is given
!> \param ATOMS List of atoms to be included in the fragment
!> \param nATOMS The number of atoms to be included
!> \param lupri Default output unit
SUBROUTINE BUILD_FRAGMENT(DALMOL,FRAGMOL,FRAGBASIS,AUXBASIS,ATOMS,nATOMS,lupri)
implicit none
INTEGER,intent(IN)                    :: NATOMS,lupri
INTEGER,intent(IN)                    :: ATOMS(NATOMS)
TYPE(MOLECULEINFO),intent(IN)         :: DALMOL
TYPE(MOLECULEINFO),intent(INOUT)      :: FRAGMOL
TYPE(BASISINFO),intent(INOUT)         :: FRAGBASIS
LOGICAL,intent(in)                    :: AUXBASIS 
!
INTEGER            :: I
Character(len=22)  :: FRAGMENTNAME

IF ((NATOMS.GT.999999).OR.(ATOMS(1).GT.999999).OR.(ATOMS(1).GT.999999)) &
     & CALL LSQUIT('Error in BUILD_FRAGMENT -> FRAGMENTNAME',-1)
write(FRAGMENTNAME,'(A4,3I6)') 'FRAG',NATOMS,ATOMS(1),ATOMS(NATOMS)
CALL init_MoleculeInfo(FRAGMOL,natoms,FRAGMENTNAME)
FRAGMOL%charge = DALMOL%charge
DO I = 1,nAtoms
   CALL COPY_ATOM(DALMOL,ATOMS(I),FRAGMOL,I,lupri)
ENDDO
FRAGMOL%nbastREG=0
FRAGMOL%nprimbastREG=0
FRAGMOL%nbastAUX=0
FRAGMOL%nprimbastAUX=0
FRAGMOL%nbastHUC=0
FRAGMOL%nprimbastHUC=0

CALL DETERMINE_NBAST(FRAGMOL,FRAGBASIS%REGULAR)
IF(AUXBASIS)THEN
   CALL DETERMINE_NBAST(FRAGMOL,FRAGBASIS%AUXILIARY)
ENDIF
END SUBROUTINE BUILD_FRAGMENT

!> \brief 
!> \author
!> \date
!> \param 
SUBROUTINE buildFragmentFromFragmentIndex(FRAGMENT,MOLECULE,FragmentIndex,iFrag,lupri)
implicit none
TYPE(MOLECULEINFO) :: FRAGMENT
TYPE(MOLECULEINFO),intent(IN)  :: MOLECULE
Integer,intent(IN)             :: FragmentIndex(MOLECULE%nAtoms)
Integer,intent(IN)             :: iFrag
Integer,intent(IN)             :: lupri
!
Integer :: iAtom
Integer :: nAtoms
!
nAtoms=0
Do iAtom=1,MOLECULE%nATOMS
  IF(FragmentIndex(iAtom) .EQ. iFrag)THEN
    nAtoms=nAtoms+1
    CALL COPY_ATOM(MOLECULE,iAtom,FRAGMENT,nAtoms,lupri)
  ENDIF
ENDDO
FRAGMENT%nAtoms = nAtoms
!
END SUBROUTINE buildFragmentFromFragmentIndex

!> \brief 
!> \author
!> \date
!> \param 
SUBROUTINE DETERMINE_NBAST(MOLECULE,BASINFO,spherical,UNCONTRACTED)
implicit none
TYPE(BASISSETINFO)  :: BASINFO
TYPE(MOLECULEINFO)  :: MOLECULE
INTEGER             :: I,TOTcont,TOTprim,R,K,type,lupri,TEMP1,TEMP2,icharge
LOGICAL,OPTIONAL    :: spherical,UNCONTRACTED
!
Logical :: spher, uncont,REG,AUX,HUC
!
! Defaults
spher  = .true.
uncont = .false.
! Optional settings
IF (present(spherical)) spher = spherical
IF (present(UNCONTRACTED)) uncont = UNCONTRACTED
REG = .FALSE.
AUX = .FALSE.
HUC = .FALSE.
IF(BASINFO%label(1:9) .EQ. 'REGULAR  ') REG = .TRUE.
IF(BASINFO%label(1:9) .EQ. 'AUXILIARY') AUX = .TRUE.
TOTcont=0
TOTprim=0
R = BASINFO%Labelindex
DO I=1,MOLECULE%nAtoms
   IF(R.EQ.0)THEN
      icharge = INT(MOLECULE%ATOM(I)%charge)
      type = BASINFO%chargeindex(icharge) 
   ELSE
      type=MOLECULE%ATOM(I)%IDtype(R)
   ENDIF
    IF(uncont)THEN
       IF(spher)THEN
          TOTcont=TOTcont+BASINFO%ATOMTYPE(type)%Totnprim
          IF(REG) MOLECULE%ATOM(I)%nprimOrbREG=BASINFO%ATOMTYPE(type)%Totnprim
          IF(AUX) MOLECULE%ATOM(I)%nprimOrbAUX=BASINFO%ATOMTYPE(type)%Totnprim
          IF(HUC) MOLECULE%ATOM(I)%nprimOrbHUC=BASINFO%ATOMTYPE(type)%Totnprim
          IF(REG) MOLECULE%ATOM(I)%ncontOrbREG=BASINFO%ATOMTYPE(type)%Totnprim
          IF(AUX) MOLECULE%ATOM(I)%ncontOrbAUX=BASINFO%ATOMTYPE(type)%Totnprim
          IF(HUC) MOLECULE%ATOM(I)%ncontOrbHUC=BASINFO%ATOMTYPE(type)%Totnprim
       ELSE
          TEMP1=0
          DO K=1,BASINFO%ATOMTYPE(type)%nAngmom
             TEMP1=TEMP1+BASINFO%ATOMTYPE(type)%SHELL(K)%nprim*K*(k+1)/2
          ENDDO
          TOTcont=TOTcont+TEMP1
          IF(REG) MOLECULE%ATOM(I)%nprimOrbREG=TEMP1
          IF(REG) MOLECULE%ATOM(I)%ncontOrbREG=TEMP1
          IF(AUX) MOLECULE%ATOM(I)%ncontOrbAUX=TEMP1
          IF(AUX) MOLECULE%ATOM(I)%nprimOrbAUX=TEMP1
          IF(HUC) MOLECULE%ATOM(I)%ncontOrbHUC=TEMP1
          IF(HUC) MOLECULE%ATOM(I)%nprimOrbHUC=TEMP1
       ENDIF
       TOTprim=TOTcont
    ELSE !DEFAULT
      IF(spher)THEN !DEFAULT
         TOTcont=TOTcont+BASINFO%ATOMTYPE(type)%Totnorb      
         TOTprim=TOTprim+BASINFO%ATOMTYPE(type)%Totnprim      
         IF(REG) MOLECULE%ATOM(I)%nprimOrbREG=BASINFO%ATOMTYPE(type)%Totnprim      
         IF(AUX) MOLECULE%ATOM(I)%nprimOrbAUX=BASINFO%ATOMTYPE(type)%Totnprim      
         IF(HUC) MOLECULE%ATOM(I)%nprimOrbHUC=BASINFO%ATOMTYPE(type)%Totnprim      
         IF(REG) MOLECULE%ATOM(I)%ncontOrbREG=BASINFO%ATOMTYPE(type)%Totnorb      
         IF(AUX) MOLECULE%ATOM(I)%ncontOrbAUX=BASINFO%ATOMTYPE(type)%Totnorb      
         IF(HUC) MOLECULE%ATOM(I)%ncontOrbHUC=BASINFO%ATOMTYPE(type)%Totnorb      
      ELSE
         TEMP1=0
         TEMP2=0
         DO K=1,BASINFO%ATOMTYPE(type)%nAngmom
            TEMP1=TEMP1+BASINFO%ATOMTYPE(type)%SHELL(K)%nprim*K*(k+1)/2
            TEMP2=TEMP2+BASINFO%ATOMTYPE(type)%SHELL(K)%norb*K*(k+1)/2
         ENDDO
            TOTprim=TOTprim+TEMP1
            TOTcont=TOTprim+TEMP2
            IF(REG) MOLECULE%ATOM(I)%nprimOrbREG=TEMP1
            IF(AUX) MOLECULE%ATOM(I)%nprimOrbAUX=TEMP1
            IF(HUC) MOLECULE%ATOM(I)%nprimOrbHUC=TEMP1
            IF(REG) MOLECULE%ATOM(I)%ncontOrbREG=TEMP2
            IF(AUX) MOLECULE%ATOM(I)%ncontOrbAUX=TEMP2
            IF(HUC) MOLECULE%ATOM(I)%ncontOrbHUC=TEMP2
      ENDIF
    ENDIF
ENDDO

BASINFO%nbast=TOTcont
BASINFO%nprimbast=TOTprim
IF(REG)MOLECULE%nbastREG=TOTcont
IF(REG)MOLECULE%nprimbastREG=TOTprim

IF(AUX)MOLECULE%nbastAUX=TOTcont
IF(AUX)MOLECULE%nprimbastAUX=TOTprim

IF(HUC)MOLECULE%nbastHUC=TOTcont
IF(HUC)MOLECULE%nprimbastHUC=TOTprim


END SUBROUTINE DETERMINE_NBAST

!> \brief 
!> \author
!> \date
!> \param 
INTEGER FUNCTION getNbasis(AOtype,intType,MOLECULE,LUPRI)
implicit none
Character*(*)     :: AOtype,intType
Integer           :: LUPRI
!Type(DaltonInput) :: DALTON
TYPE(MOLECULEINFO):: MOLECULE
!
integer :: np,nc

IF (AOtype.EQ.'Regular') THEN
  nc = MOLECULE%nbastREG
  np = MOLECULE%nprimbastREG
ELSEIF (AOtype.EQ.'DF-Aux') THEN
  nc = MOLECULE%nbastAUX
  np = MOLECULE%nprimbastAUX
ELSEIF (AOtype.EQ.'Huckel') THEN
  nc = MOLECULE%nbastHUC
  np = MOLECULE%nprimbastHUC
ELSEIF (AOtype.EQ.'Empty') THEN
  nc = 1
  np = 1
ELSEIF (AOtype.EQ.'Nuclear') THEN
  nc = MOLECULE%nAtoms
  np = MOLECULE%nAtoms
ELSEIF (AOtype.EQ.'Nuclei') THEN
  nc = MOLECULE%nAtoms
  np = MOLECULE%nAtoms
ELSE
  WRITE(LUPRI,'(1X,2A)') 'Error in getNbasis. Not valid AOtype =', AOtype
  CALL LSQUIT('AOtype not valid in getNbasis',lupri)
ENDIF
IF (intType.EQ.'Contracted') THEN
  getNbasis = nc
ELSEIF (intType.EQ.'Primitive') THEN
  getNbasis = np
ELSE
  WRITE(LUPRI,'(1X,2A)') 'Error in getNbasis. Not valid intType =',intType
  CALL LSQUIT('Error: intType not valid in getNbasis',lupri)
ENDIF
END FUNCTION getNbasis

SUBROUTINE GET_GEOMETRY(LUPRI,IPRINT,MOLECULE,natoms,X,Y,Z)
IMPLICIT NONE
INTEGER            :: LUPRI,IPRINT,natoms
TYPE(MOLECULEINFO) :: MOLECULE
REAL(REALK)        :: X(natoms),Y(natoms),Z(natoms)
!
integer :: I

DO I=1,nAtoms
     X(I) = MOLECULE%ATOM(I)%CENTER(1)
     Y(I) = MOLECULE%ATOM(I)%CENTER(2)
     Z(I) = MOLECULE%ATOM(I)%CENTER(3)
ENDDO

END SUBROUTINE GET_GEOMETRY

SUBROUTINE PRINT_GEOMETRY(MOLECULE,LUPRI)
IMPLICIT NONE
INTEGER :: LUPRI
INTEGER :: I
CHARACTER(len=1)   :: CHRXYZ(3)=(/'x','y','z'/)
TYPE(MOLECULEINFO) :: MOLECULE
   WRITE (LUPRI,'(2X,A,I3)')' Total number of coordinates:',3*MOLECULE%natoms
   DO I=1,MOLECULE%nAtoms
      WRITE (LUPRI,'(/I4,3X,A,5X,A,3X,F15.10)')&
           &  (3*I-2), MOLECULE%ATOM(I)%Name,CHRXYZ(1),&
           & MOLECULE%ATOM(I)%CENTER(1)
      WRITE (LUPRI,'(I4,12X,A,3X,F15.10)')&
           &  3*I-1, CHRXYZ(2), MOLECULE%ATOM(I)%CENTER(2),&
           &  3*I, CHRXYZ(3), MOLECULE%ATOM(I)%CENTER(3)
   ENDDO
END SUBROUTINE PRINT_GEOMETRY 

END MODULE molecule_module
