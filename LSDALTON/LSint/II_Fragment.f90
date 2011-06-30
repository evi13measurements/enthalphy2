!> @file
!> Module containing subroutines for dividing molecule into fragments and for building framgents.
MODULE Fragment
  use precision
  use TYPEDEF
  use BUILDBASISSET

CONTAINS
!> \brief initialise the FRAGMENTINFO structure
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE initFragmentInfo(fragment,num_fragments,nAtoms)
implicit none
TYPE(FRAGMENTINFO),intent(inout) :: fragment
Integer,intent(in)               :: num_fragments, nAtoms

fragment%numFragments     = min(num_fragments,nAtoms)
fragment%atomsInMolecule  = nAtoms
fragment%numberOrbialsSet = .FALSE.
NULLIFY(fragment%natoms)
NULLIFY(fragment%nContOrb)
NULLIFY(fragment%nPrimOrb)
NULLIFY(fragment%nStartContOrb)
NULLIFY(fragment%nStartPrimOrb)
NULLIFY(fragment%fragmentIndex)
ALLOCATE(fragment%nAtoms(num_fragments))
ALLOCATE(fragment%nContOrb(num_fragments,3))
ALLOCATE(fragment%nPrimOrb(num_fragments,3))
ALLOCATE(fragment%nStartContOrb(num_fragments,3))
ALLOCATE(fragment%nStartPrimOrb(num_fragments,3))
ALLOCATE(fragment%fragmentIndex(nAtoms))
END SUBROUTINE initFragmentInfo

!> \brief free the FRAGMENTINFO structure
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE freeFragmentInfo(fragment)
implicit none
TYPE(FRAGMENTINFO) :: fragment
DEALLOCATE(fragment%natoms)
DEALLOCATE(fragment%nContOrb)
DEALLOCATE(fragment%nPrimOrb)
DEALLOCATE(fragment%nStartContOrb)
DEALLOCATE(fragment%nStartPrimOrb)
DEALLOCATE(fragment%fragmentIndex)
NULLIFY(fragment%natoms)
NULLIFY(fragment%nContOrb)
NULLIFY(fragment%nPrimOrb)
NULLIFY(fragment%nStartContOrb)
NULLIFY(fragment%nStartPrimOrb)
NULLIFY(fragment%fragmentIndex)

END SUBROUTINE freeFragmentInfo

!> \brief Builds the four fragments from the different molecules and fragment-info
!> \author S. Reine
!> \date 2010-02-05
!> \param Setting Contains information about the integral setting
!> \param I Select the LHS-fragments according to block number I
!> \param J Select the RHS-fragments according to block number J
!> \param sameFragmentLHS Specifies if the two LHS fragments are identical
!> \param sameFragmentRHS Specifies if the two RHS fragments are identical
!> \param lupri Default output unit
SUBROUTINE SetDaltonFragments(SETTING,I,J,sameFragmentLHS,sameFragmentRHS,lupri)
implicit none
Type(LSSETTING),intent(inout) :: SETTING
Integer,intent(in)            :: I,J,lupri
Logical,intent(out)           :: sameFragmentLHS,sameFragmentRHS
!
SETTING%FRAGMENTS%iLHSblock = I
SETTING%FRAGMENTS%iRHSblock = J

sameFragmentLHS = SETTING%FRAGMENTS%LHSblock%blocks(I)%sameFragments
sameFragmentRHS = SETTING%FRAGMENTS%RHSblock%Blocks(J)%sameFragments

CALL buildFragments(SETTING%FRAGMENT,I,SETTING%FRAGMENTS%LHSBlock,J,SETTING%FRAGMENTS%RHSBlock,&
    &               SETTING,SETTING%MOLECULE,setting%fragments,lupri)

END SUBROUTINE SetDaltonFragments

!> \brief Free the dalton-fragments
!> \autor S. Reine
!> \param Setting Contains the integral settings
SUBROUTINE freeDaltonFragments(SETTING)
implicit none
Type(LSSETTING),intent(inout) :: SETTING
!
Integer :: indAO

CALL freeFragments(SETTING%FRAGMENT,SETTING%fragBuild,setting%nAO)

! Restore to defaul settings
DO indAO=1,setting%nAO
  setting%fragment(indAO)%p => setting%molecule(indAO)%p
ENDDO

!Set up fragments
END SUBROUTINE FreeDaltonFragments

!> \brief This subroutine sets up the fragmentinfo for each molecule in setting
!> \author S. Reine
!> \date 2010-02-05
!> \param Setting Contains information about integral setting
!> \param AO1 Determines basis-set used for AO 1
!> \param AO2 Determines basis-set used for AO 2
!> \param AO3 Determines basis-set used for AO 3
!> \param AO4 Determines basis-set used for AO 4
!> \param lupri Default output unit
!> \param luerr Default unit for error output
!>
!> The fragmentinfo is used for (some) ls-integral subroutines if the .FRAGMENT 
!> keyword is set. This info is used to construct the molecular fragments
!> in subroutine SetDaltonFragments - again used in the setInputAO to set up
!> the AO-batches used by the Thermite routines.
SUBROUTINE buildFragmentInfoAndBlocks(SETTING,AO1,AO2,AO3,AO4,LUPRI,LUERR)
use molecule_module
implicit none
Type(LSSETTING),intent(inout) :: SETTING
Character*(*),intent(in)      :: AO1,AO2,AO3,AO4
Integer            :: LUPRI,LUERR
!
integer       :: iAO,jAO
logical       :: isset(setting%nAO)
Character(80) :: AOstring(4)

AOstring(1) = AO1
AOstring(2) = AO2
AOstring(3) = AO3
AOstring(4) = AO4

IF (setting%nAO.NE.4) CALL LSQUIT('Error in buildFragmentInfoAndBlocks. nAO .NE. 4',lupri)

Setting%Fragments%identical = .FALSE.
DO iAO=1,setting%nAO
  Setting%Fragments%identical(iAO,iAO) = .TRUE.
ENDDO

isset = .FALSE.
DO iAO=1,setting%nAO
  IF (.not.isset(iAO)) THEN
    NULLIFY(Setting%Fragments%Info(iAO)%p)
    ALLOCATE(Setting%Fragments%Info(iAO)%p)
!Simen Add possibility to have different number of fragments for the different AO's
    CALL initFragmentInfo(SETTING%FRAGMENTS%INFO(iAO)%p,SETTING%numFragments,&
     &                    setting%MOLECULE(iAO)%p%nAtoms)
    Setting%Fragments%infoAllocated(iAO) = .TRUE.
    CALL FragmentMolecule(Setting%MOLECULE(iAO)%p,setting%fragments%Info(iAO)%p%fragmentIndex,&
     &                    setting%numFragments,lupri)
    CALL BuildFragmentInfo(setting%Fragments%Info(iAO)%p,&
     &                    setting%fragments%Info(iAO)%p%fragmentIndex,setting%molecule(iAO)%p)
  ENDIF
  DO jAO=iAO+1,setting%nAO
    IF (setting%sameMol(iAO,jAO).AND.setting%sameBas(iAO,jAO)&
   &    .AND.(AOstring(iAO).EQ.AOstring(jAO))) THEN ! AND if .EQ. # of fragments
      Setting%Fragments%Info(jAO)%p => Setting%Fragments%Info(iAO)%p
      isset(jAO) = .TRUE.
      Setting%Fragments%identical(iAO,jAO) = .TRUE.
      Setting%Fragments%identical(jAO,iAO) = .TRUE.
    ENDIF
  ENDDO
ENDDO
CALL BuildDaltonBlocks(SETTING,AO1,AO2,AO3,AO4,LUPRI,LUERR)

!
END SUBROUTINE buildFragmentInfoAndBlocks

!> \brief free the FRAGMENTINFO structure and blocks
!> \author S. Reine
!> \date 2010
SUBROUTINE freeFragmentInfoAndBlocks(SETTING)
implicit none
Type(LSSETTING)   :: SETTING
!
integer :: iAO
!
DO iAO=1,4
  IF (SETTING%FRAGMENTS%infoAllocated(iAO)) THEN
    CALL freeFragmentInfo(SETTING%FRAGMENTS%INFO(iAO)%p)
    SETTING%FRAGMENTS%infoAllocated(iAO) = .FALSE.
    DEALLOCATE(SETTING%FRAGMENTS%INFO(iAO)%p)
  ENDIF
ENDDO
! Reset to default when freeing fragments
setting%sameFrag = setting%sameMol
DEALLOCATE(SETTING%FRAGMENTS%LHSblock%blocks)
DEALLOCATE(SETTING%FRAGMENTS%RHSblock%blocks)
NULLIFY(SETTING%FRAGMENTS%LHSblock%blocks)
NULLIFY(SETTING%FRAGMENTS%RHSblock%blocks)

END SUBROUTINE freeFragmentInfoAndBlocks

!> \brief build the FRAGMENTINFO structure
!> \author S. Reine
!> \date 2010
SUBROUTINE BuildFragmentinfo(fragment,fragmentIndex,MOLECULE)
use molecule_type
implicit none
Type(MOLECULEINFO)   :: MOLECULE
TYPE(FRAGMENTINFO)   :: fragment
Integer              :: fragmentIndex(MOLECULE%nAtoms)
Integer              :: LUPRI,LUERR,J,I,ifragment

J=1 ! REGULAR
fragment%numberOrbialsSet = .TRUE.

ifragment=1
fragment%nAtoms(ifragment)=0
fragment%nContOrb(ifragment,J)=0
fragment%nPrimOrb(ifragment,J)=0
fragment%nStartContOrb(ifragment,J)=1
fragment%nStartPrimOrb(ifragment,J)=1

DO I=1,MOLECULE%nAtoms
   IF(fragmentIndex(I) .NE. ifragment)THEN
      ifragment = ifragment+1
      fragment%nAtoms(ifragment)=0
      fragment%nContOrb(ifragment,J)=0
      fragment%nPrimOrb(ifragment,J)=0
      fragment%nStartContOrb(ifragment,J)=&
           &fragment%nStartContOrb(ifragment-1,J)+fragment%nContOrb(ifragment-1,J)
      fragment%nStartPrimOrb(ifragment,J)=&
           &fragment%nStartPrimOrb(ifragment-1,J)+fragment%nPrimOrb(ifragment-1,J)
   ENDIF
   fragment%nContOrb(ifragment,J)=fragment%nContOrb(ifragment,J)&
        &+MOLECULE%ATOM(I)%nContOrbREG
   fragment%nPrimOrb(ifragment,J)=fragment%nPrimOrb(ifragment,J)&
        &+MOLECULE%ATOM(I)%nPrimOrbREG
   fragment%nAtoms(ifragment)=fragment%nAtoms(ifragment)+1
ENDDO
J=2!AUX

ifragment=1
fragment%nContOrb(ifragment,J)=0
fragment%nPrimOrb(ifragment,J)=0
fragment%nStartContOrb(ifragment,J)=1
fragment%nStartPrimOrb(ifragment,J)=1

DO I=1,MOLECULE%nAtoms
   IF(FragmentIndex(I) .NE. ifragment)THEN
      ifragment = ifragment+1
      fragment%nContOrb(ifragment,J)=0
      fragment%nPrimOrb(ifragment,J)=0
      fragment%nStartContOrb(ifragment,J)=&
           &fragment%nStartContOrb(ifragment-1,J)+fragment%nContOrb(ifragment-1,J)
      fragment%nStartPrimOrb(ifragment,J)=&
           &fragment%nStartPrimOrb(ifragment-1,J)+fragment%nPrimOrb(ifragment-1,J)
   ENDIF
   fragment%nContOrb(ifragment,J)=fragment%nContOrb(ifragment,J)+MOLECULE%ATOM(I)%nContOrbAUX
   fragment%nPrimOrb(ifragment,J)=fragment%nPrimOrb(ifragment,J)+MOLECULE%ATOM(I)%nPrimOrbAUX
ENDDO

J=3!Huckel

ifragment=1
fragment%nContOrb(ifragment,J)=0
fragment%nPrimOrb(ifragment,J)=0
fragment%nStartContOrb(ifragment,J)=1
fragment%nStartPrimOrb(ifragment,J)=1

DO I=1,MOLECULE%nAtoms
   IF(FragmentIndex(I) .NE. ifragment)THEN
      ifragment = ifragment+1
      fragment%nContOrb(ifragment,J)=0
      fragment%nPrimOrb(ifragment,J)=0
      fragment%nStartContOrb(ifragment,J)=&
           &fragment%nStartContOrb(ifragment-1,J)+fragment%nContOrb(ifragment-1,J)
      fragment%nStartPrimOrb(ifragment,J)=&
           &fragment%nStartPrimOrb(ifragment-1,J)+fragment%nPrimOrb(ifragment-1,J)
   ENDIF
      fragment%nContOrb(ifragment,J)=fragment%nContOrb(ifragment,J)+MOLECULE%ATOM(I)%nContOrbHUC
      fragment%nPrimOrb(ifragment,J)=fragment%nPrimOrb(ifragment,J)+MOLECULE%ATOM(I)%nPrimOrbHUC
ENDDO


END SUBROUTINE BUILDFRAGMENTINFO

!> \brief build the LHS and RHS BLOCKINFO structure
!> \author S. Reine
!> \date 2010
SUBROUTINE BuildDaltonBlocks(SETTING,AO1,AO2,AO3,AO4,LUPRI,LUERR)
implicit none
Character*(*)        :: AO1,AO2,AO3,AO4
TYPE(LSSETTING)   :: SETTING
Integer              :: LUPRI,LUERR
!
CALL BuildBlock(SETTING%FRAGMENTS%LHSblock,SETTING,&
     &          SETTING%FRAGMENTS%INFO(1)%p,SETTING%FRAGMENTS%INFO(2)%p,&
     &          Setting%fragments%identical(1,2),AO1,AO2,LUPRI,LUERR)
CALL BuildBlock(SETTING%FRAGMENTS%RHSblock,SETTING,&
     &          SETTING%FRAGMENTS%INFO(3)%p,SETTING%FRAGMENTS%INFO(4)%p,&
     &          Setting%fragments%identical(3,4),AO3,AO4,LUPRI,LUERR)
END SUBROUTINE BuildDaltonBlocks

!> \brief build the BLOCKINFO structure
!> \author S. Reine
!> \date 2010
SUBROUTINE BuildBlock(Block,SETTING,fragment1,fragment2,identical,AO1,AO2,LUPRI,LUERR)
implicit none
Character*(*)        :: AO1,AO2
TYPE(FRAGMENTINFO)   :: fragment1,fragment2
TYPE(LSSETTING)   :: SETTING
TYPE(BLOCKINFO)      :: Block
Integer              :: LUPRI,LUERR
Logical, intent(in)  :: identical
!
Integer       :: iAO,iFrag1,iFrag2,iBlock,startFrag
Integer       :: basisSetIdentifier(2)
Integer       :: numFragments(2)
Character(80) :: AOstring(2)
logical       :: sameAOs

AOstring(1) = AO1    
AOstring(2) = AO2    

sameAOs = (AOstring(1).EQ.AOstring(2)).AND.identical
Block%sameAOs = sameAOs

numFragments(1) = fragment1%numFragments
numFragments(2) = fragment2%numFragments
DO iAO=1,2
  IF (AOstring(iAO).EQ.'Regular') THEN
    basisSetIdentifier(iAO) = 1
  ELSE IF (AOstring(iAO).EQ.'DF-Aux') THEN
    basisSetIdentifier(iAO) = 2
  ELSE IF (AOstring(iAO).EQ.'Huckel') THEN
    basisSetIdentifier(iAO) = 3
  ELSE IF (AOstring(iAO).EQ.'Empty') THEN
    basisSetIdentifier(iAO) = 0
    numFragments(iAO)       = 1
  ELSE IF (AOstring(iAO).EQ.'Nuclear') THEN
    basisSetIdentifier(iAO) = 0
  ELSE
    WRITE(LUPRI,'(1X,A,I1,2A)') 'Error in BuildBlock AOstring(',iAO,',)=',AOstring(iAO)
    CALL LSQUIT('Programming error. Wrong AOstring in BuildBlock',lupri)
  ENDIF
ENDDO

Block%numBlocks = numFragments(1)*numFragments(2)
IF (sameAOs) Block%numBlocks = numFragments(1)*(numFragments(1)+1)/2

NULLIFY(Block%blocks)
ALLOCATE(Block%blocks(Block%numBlocks))

iBlock = 0
DO iFrag1=1,numFragments(1)
  startFrag = 1
  IF (sameAOs) startFrag = iFrag1
  DO iFrag2=startFrag,numFragments(2)
    iBlock = iBlock + 1
    Block%blocks(iBlock)%node = MOD(iBlock-1,SETTING%numnodes)
    Block%blocks(iBlock)%fragment1     = iFrag1
    Block%blocks(iBlock)%fragment2     = iFrag2
    Block%blocks(iBlock)%sameFragments = sameAOs .AND. (iFrag1.EQ.iFrag2)
    IF(basisSetIdentifier(1).EQ. 0)THEN
      Block%blocks(iBlock)%nbast1    = 1
      Block%blocks(iBlock)%startOrb1 = 1
      Block%blocks(iBlock)%nprimbast1    = 1
      Block%blocks(iBlock)%startprimOrb1 = 1
      Block%blocks(iBlock)%nAtoms1 = fragment1%nAtoms(iFrag1)            
    ELSE
      Block%blocks(iBlock)%nbast1    = fragment1%nContOrb(iFrag1,basisSetIdentifier(1))
      Block%blocks(iBlock)%startOrb1 = fragment1%nStartContOrb(iFrag1,basisSetIdentifier(1))
      Block%blocks(iBlock)%nprimbast1    = fragment1%nprimOrb(iFrag1,basisSetIdentifier(1))
      Block%blocks(iBlock)%startprimOrb1 = fragment1%nStartprimOrb(iFrag1,basisSetIdentifier(1))
      Block%blocks(iBlock)%nAtoms1    = fragment1%nAtoms(iFrag1)
    ENDIF
    IF(basisSetIdentifier(2).EQ. 0)THEN
      Block%blocks(iBlock)%nbast2    = 1
      Block%blocks(iBlock)%startOrb2 = 1
      Block%blocks(iBlock)%nprimbast2    = 1
      Block%blocks(iBlock)%startprimOrb2 = 1
      Block%blocks(iBlock)%nAtoms2    = fragment2%nAtoms(iFrag2)
    ELSE
      Block%blocks(iBlock)%nbast2    = fragment2%nContOrb(iFrag2,basisSetIdentifier(2))
      Block%blocks(iBlock)%startOrb2 = fragment2%nStartContOrb(iFrag2,basisSetIdentifier(2))
      Block%blocks(iBlock)%nprimbast2    = fragment2%nprimOrb(iFrag2,basisSetIdentifier(2))
      Block%blocks(iBlock)%startprimOrb2 = fragment2%nStartprimOrb(iFrag2,basisSetIdentifier(2))
      Block%blocks(iBlock)%nAtoms2    = fragment2%nAtoms(iFrag2)
    ENDIF
  ENDDO
ENDDO
    
END SUBROUTINE BuildBlock

!> \brief Get Block Orbital Info
!> \author S. Reine
!> \date 2010
SUBROUTINE GetBlockOrbitalInfo(B,n1,n2,s1,s2,e1,e2)
implicit none
TYPE(BLOCK) :: B
Integer     :: n1,n2,s1,s2,e1,e2

n1=B%nbast1
n2=B%nbast2
s1 = B%startOrb1
s2 = B%startOrb2
e1 = s1 + n1 - 1
e2 = s2 + n2 - 1

END SUBROUTINE GetBlockOrbitalInfo

!> \brief Get number of basis function on fragment
!> \author S. Reine
!> \date 2010
FUNCTION getNbasisFragment(fragment,ifragment,basisSet,basisType)
implicit none
Real(realk)        :: getNbasisFragment
Character*(*)      :: basisSet,basisType
TYPE(FRAGMENTINFO) :: fragment
Integer            :: ifragment
IF (basisSet.EQ.'Regular') THEN
  IF (basisType.EQ.'Contracted') THEN
    getNbasisFragment = fragment%nContOrb(ifragment,1)
  ELSE IF (basisType.EQ.'Primitive') THEN
    getNbasisFragment = fragment%nPrimOrb(ifragment,1)
  ELSE 
    CALL LSQUIT('Programming error. Wrong basisType in getNbasisFragment!',-1)
  ENDIF
ELSE IF (basisSet.EQ.'DF-Aux') THEN
  IF (basisType.EQ.'Contracted') THEN
    getNbasisFragment = fragment%nContOrb(ifragment,2)
  ELSE IF (basisType.EQ.'Primitive') THEN
    getNbasisFragment = fragment%nPrimOrb(ifragment,2)
  ELSE 
    CALL LSQUIT('Programming error. Wrong basisType in getNbasisFragment!',-1)
  ENDIF
ELSE IF (basisSet.EQ.'Huckel') THEN
  IF (basisType.EQ.'Contracted') THEN
    getNbasisFragment = fragment%nContOrb(ifragment,3)
  ELSE IF (basisType.EQ.'Primitive') THEN
    getNbasisFragment = fragment%nPrimOrb(ifragment,3)
  ELSE 
    CALL LSQUIT('Programming error. Wrong basisType in getNbasisFragment!',-1)
  ENDIF
ELSE IF (basisSet.EQ.'Empty') THEN
  getNbasisFragment = 1 
ELSE
  CALL LSQUIT('Programming error. Wrong basisSet in getNbasisFragment!',-1)
ENDIF
END FUNCTION getNbasisFragment

!> \brief Get dalton orbital info
!> \author S. Reine
!> \date 2010
SUBROUTINE GetDaltonOrbitalInfo(SETTING,iLHS,iRHS,nbast1,nbast2,nbast3,nbast4,&
     &                          start1,start2,start3,start4,end1,end2,end3,end4)
implicit none
Type(LSSETTING)   :: SETTING
Integer           :: iLHS,iRHS
Integer           :: nbast1,nbast2,nbast3,nbast4,start1,start2,start3,start4
Integer           :: end1,end2,end3,end4
!
CALL GetBlockOrbitalInfo(SETTING%FRAGMENTS%LHSblock%blocks(iLHS),nbast1,nbast2,&
     &                   start1,start2,end1,end2)
CALL GetBlockOrbitalInfo(SETTING%FRAGMENTS%RHSblock%blocks(iRHS),nbast3,nbast4,&
     &                   start3,start4,end3,end4)
END SUBROUTINE GetDaltonOrbitalInfo

!> \brief Build molecular fragments (from framentinfo) to be used in setInputAO
!> \author S. Reine
!> \date 2010-02-05
!> \param Fragment The four molecular fragments to be used in setInputAO
!> \param Molecule The four molecules from which the fragments are constructed
!> \param iLHS LHS block number - indexing the two LHS fragments to be used
!> \param iRHS RHS block number - indexing the two RHS fragments to be used
!> \param LHSBlocks The LHS-block information
!> \param RHSBlocks The RHS-block information
!> \param Setting Specifies the integral settings
!> \param fragments Holds the fragment-information for the four molecules
!> \param lupri Default output unit
SUBROUTINE buildFragments(FRAGMENT,iLHS,LHSBlocks,iRHS,RHSblocks,SETTING,MOLECULE,&
     &                    fragments,lupri)
use molecule_module
use molecule_type
implicit none
TYPE(MOLECULE_PT),pointer      :: FRAGMENT(:),MOLECULE(:)
TYPE(FRAGMENTITEM),intent(in)  :: fragments
TYPE(BLOCKINFO),intent(in)     :: LHSBLOCKS,RHSBLOCKS
TYPE(LSSETTING),intent(inout)  :: SETTING
INTEGER,intent(in)             :: IRHS,ILHS 
INTEGER :: iAO,jAO,iFrag,jFrag,nAtoms
INTEGER :: LUPRI
LOGICAL :: isset(setting%nAO)   
Character(len=22) :: FRAGMENTNAME

isset = .FALSE.
IF (setting%nAO.NE.4) CALL LSQUIT('Error in buildFragments. nAO .NE. 4',lupri)

!Initialize all fragments to be different
setting%sameFrag = .FALSE.
DO iAO=1,setting%nAO
  setting%sameFrag(iAO,iAO) = .TRUE.
ENDDO

DO iAO=1,setting%nAO
! Set up the number of atoms in the fragment and get the fragment number
  IF (iAO.EQ.1) THEN
    nAtoms = LHSBlocks%blocks(iLHS)%nAtoms1
    iFrag  = LHSBlocks%blocks(iLHS)%fragment1
  ELSEIF (iAO.EQ.2) THEN
    nAtoms = LHSBlocks%blocks(iLHS)%nAtoms2
    iFrag  = LHSBlocks%blocks(iLHS)%fragment2
  ELSEIF (iAO.EQ.3) THEN
    nAtoms = RHSBlocks%blocks(iRHS)%nAtoms1
    iFrag  = RHSBlocks%blocks(iRHS)%fragment1
  ELSEIF (iAO.EQ.4) THEN
    nAtoms = RHSBlocks%blocks(iRHS)%nAtoms2
    iFrag  = RHSBlocks%blocks(iRHS)%fragment2
  ELSE 
    CALL LSQUIT('Error in buildFragments. Wring case for iAO',lupri)
  ENDIF
  IF (.NOT. isset(iAO)) THEN
    IF (MOLECULE(iAO)%p%nATOMS.NE.size(fragments%info(iAO)%p%fragmentIndex)) THEN
      write(lupri,'(1X,A,I1,A,I3,A,I3)') &
         &  'Input inconsistency in buildFragments. Wrong dimension of fragmentIndex(',iAO,&
         &  '). nAtoms = ',MOLECULE(iAO)%p%nATOMS, ' and dimension of fragmentIndex =',&
         &  size(fragments%info(iAO)%p%fragmentIndex)
      CALL LSQUIT('Input inconsistency in buildFragments. Wrong dimension of fragmentIndex.',lupri)
    ENDIF
    NULLIFY(FRAGMENT(iAO)%p)
    ALLOCATE(FRAGMENT(iAO)%p) 
!   If the number of atoms in the fragment is equal to the number of atoms in the
!   molecule they are identical and we do not build the fragment
    IF (MOLECULE(iAO)%p%nATOMS.EQ.nAtoms) THEN
      FRAGMENT(iAO)%p => MOLECULE(iAO)%p
!   Otherwise build the fragment
    ELSE
      IF ((iFrag.GT.999999999).OR.(nAtoms.GT.999999999)) CALL LSQUIT('Error in buildFragments -> FRAGMENTNAME',-1)
      write(FRAGMENTNAME,'(A4,2I9)') 'FRAG',iFrag,nAtoms
      CALL init_MoleculeInfo(FRAGMENT(iAO)%p,nAtoms,FRAGMENTNAME)
      IF (nAtoms.GE.0) THEN
        CALL buildFragmentFromFragmentIndex(FRAGMENT(iAO)%p,MOLECULE(iAO)%p,&
     &                                      fragments%info(iAO)%p%FragmentIndex,iFrag,lupri)
        setting%fragBuild(iAO) = .TRUE.
      ELSE
        setting%fragBuild(iAO) = .FALSE.
      ENDIF
    ENDIF
  ENDIF
  DO jAO=iAO+1,setting%nAO
    IF (jAO.EQ.2) THEN
      jFrag  = LHSBlocks%blocks(iLHS)%fragment2
    ELSEIF (jAO.EQ.3) THEN
      jFrag  = RHSBlocks%blocks(iRHS)%fragment1
    ELSEIF (jAO.EQ.4) THEN
      jFrag  = RHSBlocks%blocks(iRHS)%fragment2
    ELSE
      CALL LSQUIT('Error in buildFragments. Wring case for jAO',lupri)
    ENDIF
    IF ((iFrag.EQ.jFrag).AND.fragments%identical(iAO,jAO)) THEN
       FRAGMENT(jAO)%p => FRAGMENT(iAO)%p
       isset(jAO) = .TRUE.
       setting%sameFrag(iAO,jAO) = .TRUE.
       setting%sameFrag(jAO,iAO) = .TRUE.
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE buildFragments

!> \brief free the molecule fragments 
!> \author S. Reine
!> \date 2010
SUBROUTINE freeFragments(FRAGMENT,fragBuild,nAO)
use molecule_module
use molecule_type
implicit none
INTEGER, intent(in)             :: nAO
TYPE(MOLECULE_PT),intent(inout) :: FRAGMENT(nAO)
LOGICAL,intent(inout)           :: fragBuild(nAO)
!
integer :: indAO

DO indAO = 1,nAO
  IF (fragBuild(indAO)) THEN
    CALL free_MoleculeInfo(FRAGMENT(indAO)%p)
    DEALLOCATE(FRAGMENT(indAO)%p)
    NULLIFY(FRAGMENT(indAO)%p)
   fragBuild(indAO) = .FALSE.
  ENDIF
ENDDO
END SUBROUTINE freeFragments

END MODULE Fragment
