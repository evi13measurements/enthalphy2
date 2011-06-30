!> @file
!> Contains the Atomic Orbital batch structure and associated subroutines

MODULE AO_Type
use precision
use LSmatrix_type
use lsmatrix_operations_dense
!*****************************************
!*
!* 
!*
!*****************************************
!> maximum AO angular moments
INTEGER, PARAMETER      :: maxAOangmom=10
real(realk),parameter   :: ExpThr = 10E-9
TYPE AOBATCHPOINTER 
TYPE(AOBATCH),pointer  :: p
END TYPE AOBATCHPOINTER

!> OBJECT CONTAINING INFORMATION ABOUT THE AOBATCHES
TYPE AOBATCH
!The type specifies the kind of primitive basis-functions
Logical                 :: TYPE_Empty
Logical                 :: TYPE_Hermite
Logical                 :: TYPE_Cartesian
Logical                 :: TYPE_Nucleus
Logical                 :: spherical
real(realk)             :: CENTER(3)
INTEGER                 :: nPrimitives
INTEGER                 :: atom
INTEGER                 :: batch
INTEGER                 :: maxContracted
integer                 :: maxAngmom
TYPE(LSMatrix),pointer    :: pExponents
!nAngmom greater than one is for family basis-sets
integer                 :: nAngmom     
real(realk)             :: extent
INTEGER                 :: ANGMOM(maxAOangmom)
INTEGER                 :: nContracted(maxAOangmom)
INTEGER                 :: startOrbital(maxAOangmom)
INTEGER                 :: startprimOrbital(maxAOangmom)
INTEGER                 :: nOrbComp(maxAOangmom)
INTEGER                 :: nPrimOrbComp(maxAOangmom)
INTEGER                 :: nOrbitals(maxAOangmom)
TYPE(LSMatrixpointer)     :: pCC(maxAOangmom)    
INTEGER                 :: CCindex(maxAOangmom) 
END TYPE AOBATCH

TYPE AOITEMPOINTER 
TYPE(AOITEM),pointer  :: p
END TYPE AOITEMPOINTER

!> OBJECT CONTAINING collection of AObatches 
TYPE AOITEM
LOGICAL                      :: EMPTY
TYPE(AOBATCH),pointer        :: BATCH(:)
INTEGER                      :: nbatches
TYPE(LSMatrix),pointer         :: CC(:) !contractioncoefficients 
INTEGER                      :: nCC
TYPE(LSMatrix),pointer         :: Exponents(:) !contractioncoefficients 
INTEGER                      :: nExp
INTEGER                      :: natoms
INTEGER                      :: nbast
INTEGER,pointer              :: ATOMICnORB(:)!size = natoms
INTEGER,pointer              :: ATOMICnBATCH(:)!size = natoms
END TYPE AOITEM

TYPE BATCHORBITALINFO
Integer :: nBatches
Integer :: nBast
Integer,pointer :: orbToBatch(:)
END TYPE BATCHORBITALINFO

Contains
!> \brief Print the AOITEM structure
!> \author T. Kjaergaard
!> \date 2010
!> \param LUPRI the logical unit number for the output file
!> \param AO the Atomic Orbital item to be printet
SUBROUTINE PRINT_AO(LUPRI,AO)
implicit none
TYPE(AOITEM)        :: AO
INTEGER             :: I,J,K,L,nAngmom
INTEGER             :: LUPRI,nPrimitives,ncol,JK
!

WRITE(LUPRI,*) '                     '
   WRITE(LUPRI,'(A)')'PRINTING AO BATCH '
   WRITE(LUPRI,'(3X,A18,I7)')'# ATOMS           ',AO%natoms
   WRITE(LUPRI,'(3X,A18,I7)')'# BASISFUNCTIONS  ',AO%nbast
   WRITE(LUPRI,'(3X,A6,2X,A10,2X,A11)')'ATOM  ','# ORBITALS','# AOBATCHES'
DO I=1,AO%natoms
   WRITE(LUPRI,'(3X,I6,6X,I6,7X,I6)')I,AO%ATOMICnORB(I),AO%ATOMICnBATCH(I)
ENDDO
DO I=1,AO%nbatches
   WRITE(LUPRI,'(5X,A,I4)')'AOBATCH NUMBER',I
   CALL PRINT_AOBATCH(AO%BATCH(I),LUPRI)
ENDDO

WRITE(LUPRI,*)' '

END SUBROUTINE PRINT_AO

!> \brief Print the AOBATCH structure
!> \author T. Kjaergaard
!> \date 2010
!> \param IUNIT the logical unit number for the output file
!> \param AOB the Atomic Orbital batch to be printet
SUBROUTINE PRINT_AOBATCH(AOB,IUNIT)
 IMPLICIT NONE
 TYPE(AOBATCH)  :: AOB
 Integer        :: IUNIT
!
 Integer        :: i,k,nangmom,nprimitives,ncol
!
 WRITE(IUNIT,'(5X,A,L)') 'Type of basis-function: Empty    :',AOB%type_empty
 WRITE(IUNIT,'(5X,A,L)') 'Type of basis-function: Hermite  :',AOB%type_hermite
 WRITE(IUNIT,'(5X,A,L)') 'Type of basis-function: Cartesian:',AOB%type_cartesian
 WRITE(IUNIT,'(5X,A,L)') 'Type of basis-function: Nucleus  :',AOB%type_nucleus
 IF (AOB%spherical) THEN
   WRITE(IUNIT,'(5X,A)') 'The orbitals are spherical'
 ELSE
   WRITE(IUNIT,'(5X,A)') 'The orbitals are non-spherical'
 ENDIF
 WRITE(IUNIT,'(5X,A,2X,3F12.8)')'Cartesian center (A):',&
       &AOB%CENTER(1),AOB%CENTER(2),AOB%CENTER(3)
 nPrimitives=AOB%nPrimitives
 WRITE(IUNIT,'(5X,A,I3)')    'Number of primitives,       nPrimitives   = ', AOB%nPrimitives
 WRITE(IUNIT,'(5X,A,I3)')    'Max. # of contracted,       maxContracted = ', AOB%maxContracted
 WRITE(IUNIT,'(5X,A,F12.8)') 'The max. extent of the AO-batch,   extent = ', AOB%extent
 WRITE(IUNIT,'(5X,A,I3)')    'Max. angular momentum,      maxAngmom     = ', AOB%maxAngmom
 WRITE(IUNIT,'(5X,A,I3)')    '# of ang. mom. blocks,      nAngmom       = ', AOB%nAngmom
 WRITE(IUNIT,'(3X,A)')       '------------------ Orbital block information -------------------------------'
 WRITE(IUNIT,'(3X,A)')       '-    Block#  Ang.mom.  #cont.   1.orb.   #orb.   #orb.c.  #pri.o.c  1.porb -'
 DO I=1,AOB%nAngmom
   WRITE(IUNIT,'(3X,8I9)') I,AOB%angmom(I),AOB%nContracted(I),AOB%startOrbital(I),&
      & AOB%nOrbitals(I),AOB%nOrbComp(I),AOB%nPrimOrbComp(I),AOB%startprimOrbital(I)
 ENDDO
 WRITE(IUNIT,'(3X,A)')       '----------------------------------------------------------------------------'
 WRITE(IUNIT,'(3X,A)')       '----------- Exponents -----------'
 call LSMAT_dense_PRINT(AOB%pExponents, 1, nPrimitives, 1, 1,IUNIT)
 WRITE(IUNIT,'(3X,A)')       '---------------------------------'
 DO K=1,AOB%nAngmom
    nAngmom=AOB%Angmom(K)
    WRITE(IUNIT,'(3X,A,I3,A,I3,A)') '--- Contraction coefficient for ang.mom: ',nAngmom,' - CCindex:',AOB%CCindex(K),' ---'
    ncol=AOB%pCC(K)%p%ncol
    CALL LSMAT_DENSE_PRINT(AOB%pCC(K)%p,1,nPrimitives,1,ncol,IUNIT)
 WRITE(IUNIT,'(3X,A)')         '-----------------------------------------------------------------------'
 ENDDO

END SUBROUTINE PRINT_AOBATCH

!> \brief Initialize a batchorbitalinfo type
!> \author S. Reine
!> \date 18-03-2010
!> \param BO The batchorbitalinfo
!> \param nbast The number of orbitals/basis functions
SUBROUTINE initBatchOrbitalInfo(BO,nBast)
use memory_handling
implicit none
TYPE(BATCHORBITALINFO),intent(INOUT) :: BO
!this is really intent out but this will diassociate the pointer
Integer,intent(IN)                 :: nBast

BO%nBatches = 0
BO%nBast = nBast
call mem_alloc(BO%orbToBatch,nBast)

END SUBROUTINE initBatchOrbitalInfo

!> \brief Set up the batchorbitalinfo
!> \author S. Reine
!> \date 18-03-2010
!> \param BO The batchorbitalinfo
!> \param AO The AO-batch
!> \param lupri Deafult output unit
SUBROUTINE setBatchOrbitalInfo(BO,AO,lupri)
implicit none
TYPE(BATCHORBITALINFO),intent(INOUT) :: BO 
!this is really intent out but this will diassociate the pointer
TYPE(AOITEM),intent(IN)            :: AO
Integer,intent(IN)                 :: lupri
!
Integer :: iBatch,iAngmom,startOrb,numOrb,iOrb

IF (BO%nBast.NE.AO%nBast) THEN
  WRITE(LUPRI,'(1X,A,I8,A,I8)') 'Error in setBatchOrbitalInfo. AO%nBast = ',AO%nBast,&
    & ' and BO%nBast = ',BO%nBast
  CALL LSQUIT('Error in setBatchOrbitalInfo. nBast mismatch!',-1)
ELSE
  BO%nBatches = AO%nBatches
  DO iBatch=1,AO%nBatches
    DO iAngmom=1,AO%Batch(iBatch)%nAngmom
      startOrb = AO%Batch(iBatch)%startOrbital(iAngmom)
      numOrb   = AO%Batch(iBatch)%nOrbitals(iAngmom)
      DO iOrb=startOrb,startOrb+numOrb-1
        BO%orbToBatch(iOrb) = iBatch
      ENDDO
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE setBatchOrbitalInfo

!> \brief Free a batchorbitalinfo type
!> \author S. Reine
!> \date 18-03-2010
!> \param BO The batchorbitalinfo
SUBROUTINE freeBatchOrbitalInfo(BO)
use memory_handling
implicit none
TYPE(BATCHORBITALINFO),intent(INOUT) :: BO

IF (.NOT.ASSOCIATED(BO%orbToBatch)) THEN
  CALL LSQUIT('Error in freeBatchOrbitalInfo. orbToBatch not associated!',-1)
ELSE
  BO%nBatches = 0
  BO%nBast = 0
  call mem_dealloc(BO%orbToBatch)
ENDIF
END SUBROUTINE freeBatchOrbitalInfo

END MODULE AO_Type

