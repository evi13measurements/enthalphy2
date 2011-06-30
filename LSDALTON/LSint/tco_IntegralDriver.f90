!> @file
!> Contains module for three-center overlap integrals. Need to be updated for performance!
MODULE tco_Integral_Driver
use precision
use typedef
use sphcart_matrices
use lstiming
use ls_util

!********** INTEGRAINPUT and INTEGRALOUTPUT
TYPE TCO_INTEGRALINPUT
TYPE(AOITEMPOINTER) :: AO(3)
INTEGER             :: AOdim(3)
LOGICAL             :: same12AOs,same13AOs,same23AOs
REAL(REALK),pointer :: VECTOR(:)
REAL(REALK),pointer :: MATRIX(:,:)
LOGICAL             :: CONT_VECTOR
LOGICAL             :: CONT_MATRIX
END TYPE TCO_INTEGRALINPUT

TYPE TCO_INTEGRALOUTPUT
REAL(REALK),pointer :: ResultMat(:,:,:)
INTEGER             :: ndim(3)
END TYPE TCO_INTEGRALOUTPUT

!********** ODBATCHES
TYPE TCO_ODITEM
LOGICAL                   :: same12AOs,same13AOs,same23AOs
INTEGER                   :: nbatches
TYPE(TCO_ODBATCH),pointer :: batch(:)
END TYPE TCO_ODITEM

TYPE TCO_ODBATCH
TYPE(AOBATCHPOINTER)      :: AO(3)
INTEGER                   :: nPrimitives
INTEGER                   :: maxContracted
INTEGER                   :: nAngmom
LOGICAL                   :: same12AO,same13AO,same23AO
END TYPE TCO_ODBATCH

!******** INTEGRAND
TYPE TCO_ORBITAL
CHARACTER(80)       :: type
LOGICAL             :: spherical
REAL(REALK)         :: center(3)
INTEGER             :: nPrimitives
INTEGER             :: maxContracted
INTEGER             :: maxAngmom
INTEGER             :: nAngmom
INTEGER,pointer     :: angmom(:)
INTEGER,pointer     :: nContracted(:)
INTEGER,pointer     :: startOrbital(:)
INTEGER,pointer     :: startprimOrbital(:)
INTEGER,pointer     :: nOrbComp(:)
INTEGER,pointer     :: nPrimOrbComp(:)
INTEGER,pointer     :: nOrbitals(:)
REAL(REALK),pointer :: exponents(:)
REAL(REALK),pointer :: CC(:,:,:)
TYPE(SPHMATPOINTER),pointer :: SPH_MAT(:)
END TYPE TCO_ORBITAL

TYPE TCO_INTEGRAND
CHARACTER(80)       :: type
TYPE(TCO_ORBITAL)   :: orbital1,orbital2,orbital3
LOGICAL             :: same12AO,same13AO,same23AO
INTEGER             :: nPrimitives
INTEGER             :: maxContracted
INTEGER             :: nAngmom
REAL(REALK),pointer :: integralPrefactor(:)
REAL(REALK),pointer :: preExpFactor(:)
REAL(REALK),pointer :: halfPinv(:)
REAL(REALK),pointer :: vectorP1(:,:)
REAL(REALK),pointer :: vectorP2(:,:)
REAL(REALK),pointer :: vectorP3(:,:)
INTEGER,pointer     :: iPrim1(:)
INTEGER,pointer     :: iPrim2(:)
INTEGER,pointer     :: iPrim3(:)
INTEGER,pointer     :: nContracted(:)
INTEGER,pointer     :: nOrbComp(:)
INTEGER,pointer     :: nOrbitals(:)
INTEGER,pointer     :: iAng1(:)
INTEGER,pointer     :: iAng2(:)
INTEGER,pointer     :: iAng3(:)
END TYPE TCO_INTEGRAND

TYPE TCO_ALLOCITEM
INTEGER             :: max_orb_nPrimitives
INTEGER             :: max_orb_maxContracted
INTEGER             :: max_orb_nAngmom
INTEGER             :: max_orb_maxAngmom
INTEGER             :: max_nPrimitives
INTEGER             :: max_maxContracted
INTEGER             :: max_nAngmom
INTEGER             :: max_E0IJK
INTEGER             :: max_INOUT
INTEGER             :: max_WORK
END TYPE TCO_ALLOCITEM

TYPE TCO_INTEGRALITEM
INTEGER             :: size1
INTEGER             :: size2
REAL(REALK),pointer :: E0IJK(:)  ! nPrim, 0:max_l1, 0:max_l2, 0:max_l3, 4
REAL(REALK),pointer :: IN(:)     ! size1, size2
REAL(REALK),pointer :: OUT(:)    ! size1, size2
REAL(REALK),pointer :: WORK(:)
END TYPE TCO_INTEGRALITEM

CONTAINS

!*************************************************************************
!*
!*                Initialization of TCO_INTEGRALINPUT
!*
!*************************************************************************
SUBROUTINE init_TCO_integralinput(INPUT,SETTING)
implicit none
TYPE(TCO_INTEGRALINPUT) :: INPUT
TYPE(LSSETTING)         :: SETTING

NULLIFY(INPUT%AO(1)%p)
NULLIFY(INPUT%AO(2)%p)
NULLIFY(INPUT%AO(3)%p)
INPUT%AOdim(1)=0
INPUT%AOdim(2)=0
INPUT%AOdim(3)=0
INPUT%same12aos=.FALSE.
INPUT%same13aos=.FALSE.
INPUT%same23aos=.FALSE.
NULLIFY(INPUT%VECTOR)
NULLIFY(INPUT%MATRIX)
INPUT%CONT_VECTOR=.FALSE.
INPUT%CONT_MATRIX=.FALSE.

END SUBROUTINE init_TCO_integralinput

SUBROUTINE attach_TCO_matrixToInput(INPUT,Matr,n1,n2)
implicit none
TYPE(TCO_INTEGRALINPUT) :: INPUT
REAL(REALK)             :: Matr(n1,n2)
INTEGER                 :: n1,n2

IF(INPUT%CONT_VECTOR) THEN
  NULLIFY(INPUT%VECTOR)
  ALLOCATE(INPUT%VECTOR(n1))
  CALL DCOPY(n1,Matr,1,INPUT%VECTOR,1)
ELSEIF(INPUT%CONT_MATRIX) THEN
  NULLIFY(INPUT%MATRIX)
  ALLOCATE(INPUT%MATRIX(n1,n2))
  CALL DCOPY(n1*n2,Matr,1,INPUT%MATRIX,1)
ELSE
  CALL LSQUIT('Programming error, wrong contraction type in TCO_CONTRACT',-1)
ENDIF

END SUBROUTINE attach_TCO_matrixToInput

SUBROUTINE deattach_TCO_matrixFromInput(INPUT)
implicit none
TYPE(TCO_INTEGRALINPUT) :: INPUT

IF(INPUT%CONT_VECTOR) THEN
  DEALLOCATE(INPUT%VECTOR)
  NULLIFY(INPUT%VECTOR)
ELSEIF(INPUT%CONT_MATRIX) THEN
  DEALLOCATE(INPUT%MATRIX)
  NULLIFY(INPUT%MATRIX)
ELSE
  CALL LSQUIT('Programming error, wrong contraction type in TCO_CONTRACT',-1)
ENDIF

END SUBROUTINE deattach_TCO_matrixFromInput

!*************************************************************************
!*
!*                      Creation of TCO_ODBATCHES
!*
!*************************************************************************
FUNCTION get_TCO_number(n1,n2,n3,same12,same13,same23)
implicit none
INTEGER                 :: get_TCO_number
INTEGER                 :: n1,n2,n3
LOGICAL                 :: same12,same13,same23

IF(same12.AND.same13.AND.same23) THEN              ! T T T
  get_TCO_number=n1*(n1+1)*(n1+2)/6
ELSEIF(same12.AND.(same13.EQV.same23)) THEN        ! T F F
  get_TCO_number=n1*(n1+1)/2*n3
ELSEIF(same13.AND.(same12.EQV.same23)) THEN        ! F T F
  get_TCO_number=n1*(n1+1)/2*n2
ELSEIF(same23.AND.(same12.EQV.same13)) THEN        ! F F T
  get_TCO_number=n2*(n2+1)/2*n1
ELSEIF(.NOT.(same12.OR.same13.OR.same23)) THEN     ! F F F
  get_TCO_number=n1*n2*n3
ELSE                                               ! T T F,  T F T,  F T T
  CALL LSQUIT('Programming error, problem with determining same AObatches pairs',-1)
ENDIF

END FUNCTION get_TCO_number

SUBROUTINE create_TCO_ODbatches(OD,INPUT,LUPRI,IPRINT)
implicit none
TYPE(TCO_ODITEM)        :: OD
TYPE(TCO_INTEGRALINPUT) :: INPUT
INTEGER                 :: LUPRI,IPRINT
!
INTEGER                 :: nbatches,np,mc,na,IOD,I1,I2,I3,I2start,I3start

IF(IPRINT.GT.5) WRITE(LUPRI,*) 'CREATE TCO_OD CALLED'

OD%same12AOs = INPUT%same12AOs
OD%same13AOs = INPUT%same13AOs
OD%same23AOs = INPUT%same23AOs

nbatches = get_TCO_number(INPUT%AO(1)%p%nbatches,&
                          INPUT%AO(2)%p%nbatches,&
                          INPUT%AO(3)%p%nbatches,&
                          OD%same12AOs,OD%same13AOs,OD%same23AOs)
OD%nbatches = nbatches
NULLIFY(OD%batch)
ALLOCATE(OD%batch(nbatches))

IOD = 0
DO I1=1,INPUT%AO(1)%p%nbatches
  I2start = 1
  IF(OD%same12AOs) I2start = I1
  DO I2=I2start,INPUT%AO(2)%p%nbatches
    I3start = 1
    IF(OD%same13AOs) I3start = I1
    IF(OD%same23AOs) I3start = I2
    DO I3=I3start,INPUT%AO(3)%p%nbatches
      IOD = IOD + 1
      OD%batch(IOD)%AO(1)%p => INPUT%AO(1)%p%batch(I1)
      OD%batch(IOD)%AO(2)%p => INPUT%AO(2)%p%batch(I2)
      OD%batch(IOD)%AO(3)%p => INPUT%AO(3)%p%batch(I3)
      
      OD%batch(IOD)%same12AO = (OD%same12AOs .AND. I1.EQ.I2)
      OD%batch(IOD)%same13AO = (OD%same13AOs .AND. I1.EQ.I3)
      OD%batch(IOD)%same23AO = (OD%same23AOs .AND. I2.EQ.I3)
      
      np = get_TCO_number(OD%batch(IOD)%AO(1)%p%nPrimitives,&
                          OD%batch(IOD)%AO(2)%p%nPrimitives,&
			  OD%batch(IOD)%AO(3)%p%nPrimitives,&
			  OD%batch(IOD)%same12AO,&
			  OD%batch(IOD)%same13AO,&
			  OD%batch(IOD)%same23AO)
      mc = get_TCO_number(OD%batch(IOD)%AO(1)%p%maxContracted,&
                          OD%batch(IOD)%AO(2)%p%maxContracted,&
			  OD%batch(IOD)%AO(3)%p%maxContracted,&
			  .FALSE.,.FALSE.,.FALSE.)
      na = get_TCO_number(OD%batch(IOD)%AO(1)%p%nAngmom,&
                          OD%batch(IOD)%AO(2)%p%nAngmom,&
			  OD%batch(IOD)%AO(3)%p%nAngmom,&
			  OD%batch(IOD)%same12AO,&
			  OD%batch(IOD)%same13AO,&
			  OD%batch(IOD)%same23AO)

      OD%batch(IOD)%nPrimitives   =  np
      OD%batch(IOD)%maxContracted =  mc
      OD%batch(IOD)%nAngmom       =  na
    ENDDO
  ENDDO
ENDDO
IF(IOD.NE.OD%nbatches) THEN
  CALL LSQUIT('Programming error, incorrect number of TCO_ODbatches',lupri)
ENDIF

CALL print_TCO_OD(OD,LUPRI,IPRINT)

END SUBROUTINE create_TCO_ODbatches

SUBROUTINE print_TCO_OD(OD,LUPRI,IPRINT)
implicit none
TYPE(TCO_ODITEM)        :: OD
INTEGER                 :: LUPRI,IPRINT
!
INTEGER                 :: IOD

IF(IPRINT.GT.2) THEN
  WRITE(LUPRI,'(1X,A,I8)') 'Number of batches:',OD%nbatches
  IF(OD%same12AOs.AND.OD%same13AOs.AND.OD%same23AOs) THEN              ! T T T
    WRITE(LUPRI,'(1X,A)') 'Batches build from three identical AO-batches'
  ELSEIF(OD%same12AOs.AND.(OD%same13AOs.EQV.OD%same23AOs)) THEN        ! T F F
    WRITE(LUPRI,'(1X,A)') 'Batches build from identical AO-batches 1 and 2, and different AO-batch 3'
  ELSEIF(OD%same13AOs.AND.(OD%same12AOs.EQV.OD%same23AOs)) THEN        ! F T F
    WRITE(LUPRI,'(1X,A)') 'Batches build from identical AO-batches 1 and 3, and different AO-batch 2'
  ELSEIF(OD%same23AOs.AND.(OD%same12AOs.EQV.OD%same13AOs)) THEN        ! F F T
    WRITE(LUPRI,'(1X,A)') 'Batches build from identical AO-batches 2 and 3, and different AO-batch 1'
  ELSEIF(.NOT.(OD%same12AOs.OR.OD%same13AOs.OR.OD%same23AOs)) THEN     ! F F F
    WRITE(LUPRI,'(1X,A)') 'Batches build from three different AO-batches'
  ELSE                                               ! T T F,  T F T,  F T T
    CALL LSQUIT('Programming error, problem with determining same AObatches pairs',lupri)
  ENDIF
  IF(IPRINT.GT.3) THEN
    DO IOD=1,OD%nbatches
      WRITE(LUPRI,'(1X,A,I4,A)') '*** BATCH number ',IOD,' ***'
      CALL print_TCO_ODbatch(OD%batch(IOD),LUPRI,IPRINT)
      WRITE(LUPRI,'(1X,A)') '*************************'
    ENDDO
  ENDIF
ENDIF

END SUBROUTINE print_TCO_OD

SUBROUTINE print_TCO_ODbatch(ODB,LUPRI,IPRINT)
implicit none
TYPE(TCO_ODBATCH)       :: ODB
INTEGER                 :: LUPRI,IPRINT

WRITE(LUPRI,'(5X,A,I3)') 'Number of primitives,      nPrimitives   = ', ODB%nPrimitives
WRITE(LUPRI,'(5X,A,I3)') 'Max. number of contracted, maxContracted = ', ODB%maxContracted
WRITE(LUPRI,'(5X,A,I3)') 'Number of angular blocks,  nAngmom       = ', ODB%nAngmom
WRITE(LUPRI,'(5X,A,L)')  'Same center and primitives, batch 1 and 2= ', ODB%same12AO
WRITE(LUPRI,'(5X,A,L)')  'Same center and primitives, batch 1 and 3= ', ODB%same13AO
WRITE(LUPRI,'(5X,A,L)')  'Same center and primitives, batch 2 and 3= ', ODB%same23AO

IF(IPRINT.GT.5) THEN
  WRITE(LUPRI,'(3X,A)') '* First AO-batch  *'
  CALL PRINT_AOBATCH(ODB%AO(1)%p,LUPRI)
  WRITE(LUPRI,'(3X,A)') '* Second AO-batch *'
  CALL PRINT_AOBATCH(ODB%AO(2)%p,LUPRI)
  WRITE(LUPRI,'(3X,A)') '* Third AO-batch  *'
  CALL PRINT_AOBATCH(ODB%AO(3)%p,LUPRI)
ENDIF

END SUBROUTINE print_TCO_ODbatch

!*************************************************************************
!*
!*                 Three Center Overlap Integral Driver
!*
!*************************************************************************
SUBROUTINE TCO_INTEGRALDRIVER(LUPRI,IPRINT,INPUT,OUTPUT)
implicit none
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
INTEGER                  :: LUPRI,IPRINT
!
TYPE(TCO_ODITEM)         :: OD
REAL(REALK)              :: TSTART,TEND

CALL LSTIMER('START ',TSTART,TEND,LUPRI)
CALL LSHEADER(LUPRI,'THE THREE_CENTER_OVERLAP DRIVER')

CALL create_TCO_ODbatches(OD,INPUT,LUPRI,IPRINT)

CALL TCO_DirectCalculation(OD,INPUT,OUTPUT,LUPRI,IPRINT)

DEALLOCATE(OD%batch)
NULLIFY(OD%batch)

CALL LSTIMER('TCO  ',TSTART,TEND,LUPRI)

END SUBROUTINE TCO_INTEGRALDRIVER

SUBROUTINE TCO_DirectCalculation(OD,INPUT,OUTPUT,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
TYPE(TCO_ODITEM)         :: OD
INTEGER                  :: LUPRI,IPRINT
!
TYPE(TCO_ALLOCITEM)      :: ALLOC
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
TYPE(SPHMAT),pointer     :: SPH_MAT(:)
INTEGER                  :: IOD,iAngmom,nSPHMAT

IF (IPRINT.GT.5) THEN
  CALL LSHEADER(LUPRI,'TCO Direct Calculation')
  WRITE(LUPRI,'(1X,A)') ''
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') '***                   TCO Direct Calculation                   ***'
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') ''
ENDIF

CALL set_TCO_alloc(ALLOC,OD,LUPRI,IPRINT)

nSPHMAT = ALLOC%max_orb_maxAngmom+1
CALL setPrecalculatedSphmat(SPH_MAT,nSPHMAT,LUPRI,IPRINT)
CALL init_TCO_integrand(P,ALLOC,LUPRI,IPRINT)
CALL init_TCO_integral(INTEGRAL,ALLOC,LUPRI,IPRINT)

DO IOD=1,OD%nbatches
  CALL set_TCO_integrand(P,OD%batch(IOD),SPH_MAT,LUPRI,IPRINT)
  CALL build_TCO_E0IJK(P,INTEGRAL,LUPRI,IPRINT)
  DO iAngmom=1,P%nAngmom
    CALL build_TCO_E0ABC(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
    CALL build_TCO_integral(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
    CALL TCO_contractBasis(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
    CALL TCO_transposition(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
    CALL TCO_sphericalTransform(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
    CALL TCO_distributeIntegrals(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
  ENDDO
ENDDO

CALL free_TCO_integral(INTEGRAL)
CALL free_TCO_integrand(P)
CALL freePrecalculatedSphmat(SPH_MAT,nSPHMAT)

END SUBROUTINE TCO_DirectCalculation

SUBROUTINE set_TCO_alloc(ALLOC,OD,LUPRI,IPRINT)
implicit none
TYPE(TCO_ALLOCITEM)      :: ALLOC
TYPE(TCO_ODITEM)         :: OD
INTEGER                  :: LUPRI,IPRINT
!
INTEGER                  :: I,I1,I2,I3,I2start,I3start
INTEGER                  :: n1,n2,n3,np,mc,na,ml
INTEGER                  :: l1,l2,l3,l,ijk1,ijk2,ijk3,ijk,lm1,lm2,lm3,lm
INTEGER                  :: e0ijk,inout,work1,work2,work3,work

ALLOC%max_orb_nPrimitives   = 0
ALLOC%max_orb_maxContracted = 0
ALLOC%max_orb_nAngmom       = 0
ALLOC%max_orb_maxAngmom     = 0
ALLOC%max_nPrimitives       = 0
ALLOC%max_maxContracted     = 0
ALLOC%max_nAngmom           = 0
ALLOC%max_E0IJK             = 0
ALLOC%max_INOUT             = 0
ALLOC%max_WORK              = 0

DO I=1,OD%nbatches
  n1 = OD%batch(I)%AO(1)%p%nPrimitives
  n2 = OD%batch(I)%AO(2)%p%nPrimitives
  n3 = OD%batch(I)%AO(3)%p%nPrimitives
  np = MAX(n1,n2,n3)
  ALLOC%max_orb_nPrimitives = MAX(ALLOC%max_orb_nPrimitives,np)
  n1 = OD%batch(I)%AO(1)%p%maxContracted
  n2 = OD%batch(I)%AO(2)%p%maxContracted
  n3 = OD%batch(I)%AO(3)%p%maxContracted
  mc = MAX(n1,n2,n3)
  ALLOC%max_orb_maxContracted = MAX(ALLOC%max_orb_maxContracted,mc)
  n1 = OD%batch(I)%AO(1)%p%nAngmom
  n2 = OD%batch(I)%AO(2)%p%nAngmom
  n3 = OD%batch(I)%AO(3)%p%nAngmom
  na = MAX(n1,n2,n3)
  ALLOC%max_orb_nAngmom = MAX(ALLOC%max_orb_nAngmom,na)
  n1 = OD%batch(I)%AO(1)%p%maxAngmom
  n2 = OD%batch(I)%AO(2)%p%maxAngmom
  n3 = OD%batch(I)%AO(3)%p%maxAngmom
  ml = MAX(n1,n2,n3)
  ALLOC%max_orb_maxAngmom = MAX(ALLOC%max_orb_maxAngmom,ml)
ENDDO

DO I=1,OD%nbatches
  np = OD%batch(I)%nPrimitives
  ALLOC%max_nPrimitives = MAX(ALLOC%max_nPrimitives,np)
  mc = OD%batch(I)%maxContracted
  ALLOC%max_maxContracted = MAX(ALLOC%max_maxContracted,mc)
  na = OD%batch(I)%nAngmom
  ALLOC%max_nAngmom = MAX(ALLOC%max_nAngmom,na)
  DO I1=1,OD%batch(I)%AO(1)%p%nAngmom
    I2start=1
    IF(OD%batch(I)%same12AO) I2start=I1
    DO I2=I2start,OD%batch(I)%AO(2)%p%nAngmom
      I3start=1
      IF(OD%batch(I)%same13AO) I3start=I1
      IF(OD%batch(I)%same23AO) I3start=I2
      DO I3=I3start,OD%batch(I)%AO(3)%p%nAngmom
        l1 = OD%batch(I)%AO(1)%p%angmom(I1)
        l2 = OD%batch(I)%AO(2)%p%angmom(I2)
        l3 = OD%batch(I)%AO(3)%p%angmom(I3)

        ijk1 = (l1+1)*(l1+2)/2
        ijk2 = (l2+1)*(l2+2)/2
        ijk3 = (l3+1)*(l3+2)/2
        ijk  = ijk1*ijk2*ijk3
        lm1  = 2*l1+1
        lm2  = 2*l2+1
        lm3  = 2*l3+1
        lm   = lm1*lm2*lm3
  
        e0ijk = np*(l1+1)*(l2+1)*(l3+1)*3
        inout = MAX(np*ijk,mc*ijk,mc*lm)
        work1 = np*(l1+1)*(l2+1)*(l3+1)
        work2 = np*mc
        work3 = ijk*lm
        work  = MAX(work1,work2,work3)
  
        ALLOC%max_E0IJK = MAX(ALLOC%max_E0IJK,e0ijk)
        ALLOC%max_INOUT = MAX(ALLOC%max_INOUT,inout)
        ALLOC%max_WORK  = MAX(ALLOC%max_WORK ,work )
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE set_TCO_alloc

SUBROUTINE init_TCO_integrand(P,ALLOC,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_ALLOCITEM)      :: ALLOC
INTEGER                  :: LUPRI,IPRINT
!
INTEGER                  :: np,na

CALL init_TCO_orbital(P%orbital1,ALLOC,LUPRI)
CALL init_TCO_orbital(P%orbital2,ALLOC,LUPRI)
CALL init_TCO_orbital(P%orbital3,ALLOC,LUPRI)

np = ALLOC%max_nPrimitives
na = ALLOC%max_nAngmom

NULLIFY(P%integralPrefactor)
NULLIFY(P%preExpFactor)
NULLIFY(P%halfPinv)
NULLIFY(P%vectorP1)
NULLIFY(P%vectorP2)
NULLIFY(P%vectorP3)
NULLIFY(P%iPrim1)
NULLIFY(P%iPrim2)
NULLIFY(P%iPrim3)

ALLOCATE(P%integralPrefactor(np))
ALLOCATE(P%preExpFactor(np))
ALLOCATE(P%halfPinv(np))
ALLOCATE(P%vectorP1(np,3))
ALLOCATE(P%vectorP2(np,3))
ALLOCATE(P%vectorP3(np,3))
ALLOCATE(P%iPrim1(np))
ALLOCATE(P%iPrim2(np))
ALLOCATE(P%iPrim3(np))

NULLIFY(P%nContracted)
NULLIFY(P%nOrbComp)
NULLIFY(P%nOrbitals)
NULLIFY(P%iAng1)
NULLIFY(P%iAng2)
NULLIFY(P%iAng3)

ALLOCATE(P%nContracted(na))
ALLOCATE(P%nOrbComp(na))
ALLOCATE(P%nOrbitals(na))
ALLOCATE(P%iAng1(na))
ALLOCATE(P%iAng2(na))
ALLOCATE(P%iAng3(na))

END SUBROUTINE init_TCO_integrand

SUBROUTINE free_TCO_integrand(P)
implicit none
TYPE(TCO_INTEGRAND)      :: P

CALL free_TCO_orbital(P%orbital1)
CALL free_TCO_orbital(P%orbital2)
CALL free_TCO_orbital(P%orbital3)

DEALLOCATE(P%integralPrefactor)
DEALLOCATE(P%preExpFactor)
DEALLOCATE(P%halfPinv)
DEALLOCATE(P%vectorP1)
DEALLOCATE(P%vectorP2)
DEALLOCATE(P%vectorP3)
DEALLOCATE(P%iPrim1)
DEALLOCATE(P%iPrim2)
DEALLOCATE(P%iPrim3)

NULLIFY(P%integralPrefactor)
NULLIFY(P%preExpFactor)
NULLIFY(P%halfPinv)
NULLIFY(P%vectorP1)
NULLIFY(P%vectorP2)
NULLIFY(P%vectorP3)
NULLIFY(P%iPrim1)
NULLIFY(P%iPrim2)
NULLIFY(P%iPrim3)

DEALLOCATE(P%nContracted)
DEALLOCATE(P%nOrbComp)
DEALLOCATE(P%nOrbitals)
DEALLOCATE(P%iAng1)
DEALLOCATE(P%iAng2)
DEALLOCATE(P%iAng3)

NULLIFY(P%nContracted)
NULLIFY(P%nOrbComp)
NULLIFY(P%nOrbitals)
NULLIFY(P%iAng1)
NULLIFY(P%iAng2)
NULLIFY(P%iAng3)

END SUBROUTINE free_TCO_integrand

SUBROUTINE init_TCO_orbital(ORB,ALLOC,LUPRI)
implicit none
TYPE(TCO_ORBITAL)        :: ORB
TYPE(TCO_ALLOCITEM)      :: ALLOC
INTEGER                  :: LUPRI
!
INTEGER                  :: np,mc,na

np = ALLOC%max_orb_nPrimitives
mc = ALLOC%max_orb_maxContracted
na = ALLOC%max_orb_nAngmom

NULLIFY(ORB%angmom)
NULLIFY(ORB%nContracted)
NULLIFY(ORB%startOrbital)
NULLIFY(ORB%startPrimOrbital)
NULLIFY(ORB%nOrbComp)
NULLIFY(ORB%nPrimOrbComp)
NULLIFY(ORB%nOrbitals)
NULLIFY(ORB%exponents)
NULLIFY(ORB%CC)
NULLIFY(ORB%SPH_MAT)

ALLOCATE(ORB%angmom(na))
ALLOCATE(ORB%nContracted(na))
ALLOCATE(ORB%startOrbital(na))
ALLOCATE(ORB%startPrimOrbital(na))
ALLOCATE(ORB%nOrbComp(na))
ALLOCATE(ORB%nPrimOrbComp(na))
ALLOCATE(ORB%nOrbitals(na))
ALLOCATE(ORB%exponents(np))
ALLOCATE(ORB%CC(np,mc,na))
ALLOCATE(ORB%SPH_MAT(na))

END SUBROUTINE init_TCO_orbital

SUBROUTINE free_TCO_orbital(ORB)
implicit none
TYPE(TCO_ORBITAL)        :: ORB

DEALLOCATE(ORB%angmom)
DEALLOCATE(ORB%nContracted)
DEALLOCATE(ORB%startOrbital)
DEALLOCATE(ORB%startPrimOrbital)
DEALLOCATE(ORB%nOrbComp)
DEALLOCATE(ORB%nPrimOrbComp)
DEALLOCATE(ORB%nOrbitals)
DEALLOCATE(ORB%exponents)
DEALLOCATE(ORB%CC)
DEALLOCATE(ORB%SPH_MAT)

NULLIFY(ORB%angmom)
NULLIFY(ORB%nContracted)
NULLIFY(ORB%startOrbital)
NULLIFY(ORB%startPrimOrbital)
NULLIFY(ORB%nOrbComp)
NULLIFY(ORB%nPrimOrbComp)
NULLIFY(ORB%nOrbitals)
NULLIFY(ORB%exponents)
NULLIFY(ORB%CC)
NULLIFY(ORB%SPH_MAT)

END SUBROUTINE free_TCO_orbital

SUBROUTINE init_TCO_integral(INTEGRAL,ALLOC,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
TYPE(TCO_ALLOCITEM)      :: ALLOC
INTEGER                  :: LUPRI,IPRINT

NULLIFY(INTEGRAL%E0IJK)
NULLIFY(INTEGRAL%IN)
NULLIFY(INTEGRAL%OUT)
NULLIFY(INTEGRAL%WORK)

ALLOCATE(INTEGRAL%E0IJK(ALLOC%max_E0IJK))
ALLOCATE(INTEGRAL%IN(ALLOC%max_INOUT))
ALLOCATE(INTEGRAL%OUT(ALLOC%max_INOUT))
ALLOCATE(INTEGRAL%WORK(ALLOC%max_WORK))

END SUBROUTINE init_TCO_integral

SUBROUTINE free_TCO_integral(INTEGRAL)
implicit none
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL

DEALLOCATE(INTEGRAL%E0IJK)
DEALLOCATE(INTEGRAL%IN)
DEALLOCATE(INTEGRAL%OUT)
DEALLOCATE(INTEGRAL%WORK)

NULLIFY(INTEGRAL%E0IJK)
NULLIFY(INTEGRAL%IN)
NULLIFY(INTEGRAL%OUT)
NULLIFY(INTEGRAL%WORK)

END SUBROUTINE free_TCO_integral

SUBROUTINE set_TCO_integrand(P,ODB,SPH_MAT,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_ODBATCH)        :: ODB
TYPE(SPHMAT),pointer     :: SPH_MAT(:)
INTEGER                  :: LUPRI,IPRINT
!
REAL(REALK),parameter    :: PI=3.14159265358979323846D0
REAL(REALK),parameter    :: Half=0.5d0, OneHalf=1.5d0
INTEGER                  :: idim,i1,i2,i3,i123,i2start,i3start
REAL(REALK)              :: vector12(3),vector13(3),vector23(3)
REAL(REALK)              :: dist12_squared,dist13_squared,dist23_squared
REAL(REALK)              :: e1,e2,e3,e,Kexp

CALL set_TCO_orbital(P%orbital1,ODB%AO(1)%p,SPH_MAT,LUPRI)
CALL set_TCO_orbital(P%orbital2,ODB%AO(2)%p,SPH_MAT,LUPRI)
CALL set_TCO_orbital(P%orbital3,ODB%AO(3)%p,SPH_MAT,LUPRI)

!IF(P%orbital1%type_Empty .AND.&
!   P%orbital2%type_Empty .AND.&
!   P%orbital3%type_Empty) THEN
!  P%type_empty = .TRUE.
!  P%type_hermite = .FALSE.
!  P%type_cartesian = .FALSE.
!  P%type_nucleus = .FALSE.
!ELSE
!!  P%type = 'NonEmpty'
!  P%type_empty = .FALSE.
!  P%type_hermite = .FALSE.
!  P%type_cartesian = .FALSE.
!  P%type_nucleus = .FALSE.
!ENDIF

P%same12AO = ODB%same12AO
P%same13AO = ODB%same13AO
P%same23AO = ODB%same23AO

P%nPrimitives   = ODB%nPrimitives
P%maxContracted = ODB%maxContracted
P%nAngmom       = ODB%nAngmom

dist12_squared = 0.d0
dist13_squared = 0.d0
dist23_squared = 0.d0
DO idim=1,3
  vector12(idim) = P%orbital1%center(idim) - P%orbital2%center(idim)
  vector13(idim) = P%orbital1%center(idim) - P%orbital3%center(idim)
  vector23(idim) = P%orbital2%center(idim) - P%orbital3%center(idim)
  dist12_squared = dist12_squared + vector12(idim)**2
  dist13_squared = dist13_squared + vector13(idim)**2
  dist23_squared = dist23_squared + vector23(idim)**2
ENDDO

i123=0
DO i1=1,P%orbital1%nPrimitives
  i2start=1
  IF(P%same12AO) i2start=i1
  DO i2=i2start,P%orbital2%nPrimitives
    i3start=1
    IF(P%same13AO) i3start=i1
    IF(P%same23AO) i3start=i2
    DO i3=i3start,P%orbital3%nPrimitives
      i123 = i123+1
      e1 = P%orbital1%exponents(i1)
      e2 = P%orbital2%exponents(i2)
      e3 = P%orbital3%exponents(i3)
      e  = e1 + e2 + e3
      Kexp = (e1*e2*dist12_squared + e1*e3*dist13_squared + e2*e3*dist23_squared)/e
      
      P%integralPrefactor(i123) = (PI/e)**OneHalf
      P%preExpFactor(i123)      = exp(-Kexp)
      P%halfPinv(i123)          = Half/e
      DO idim=1,3
        P%vectorP1(i123,idim) = (-e2*vector12(idim)-e3*vector13(idim))/e
	P%vectorP2(i123,idim) = ( e1*vector12(idim)-e3*vector23(idim))/e
	P%vectorP3(i123,idim) = ( e1*vector13(idim)+e2*vector23(idim))/e
      ENDDO
      P%iPrim1(i123) = i1
      P%iPrim2(i123) = i2
      P%iPrim3(i123) = i3
    ENDDO
  ENDDO
ENDDO
IF(i123.NE.P%nPrimitives) THEN
  CALL LSQUIT('Programming error, wrong nPrimitives in Integrand',lupri)
ENDIF

i123=0
DO i1=1,P%orbital1%nAngmom
  i2start=1
  IF(P%same12AO) i2start=i1
  DO i2=i2start,P%orbital2%nAngmom
    i3start=1
    IF(P%same13AO) i3start=i1
    IF(P%same23AO) i3start=i2
    DO i3=i3start,P%orbital3%nAngmom
      i123=i123+1
      P%nContracted(i123) = P%orbital1%nContracted(i1)*&
                            P%orbital2%nContracted(i2)*&
        	            P%orbital3%nContracted(i3)
      P%nOrbComp(i123)    = P%orbital1%nOrbComp(i1)*&
                            P%orbital2%nOrbComp(i2)*&
			    P%orbital3%nOrbComp(i3)
      P%nOrbitals(i123)   = P%orbital1%nOrbitals(i1)*&
                            P%orbital2%nOrbitals(i2)*&
			    P%orbital3%nOrbitals(i3)
      P%iAng1(i123)       = i1
      P%iAng2(i123)       = i2
      P%iAng3(i123)       = i3
    ENDDO
  ENDDO
ENDDO
IF(i123.NE.P%nAngmom) THEN
  CALL LSQUIT('Programming error, wrong nAngmom in Integrand',lupri)
ENDIF

IF(IPRINT.GT.15) THEN
  CALL print_TCO_integrand(P,LUPRI,IPRINT)
ENDIF

END SUBROUTINE set_TCO_integrand

SUBROUTINE set_TCO_orbital(ORB,AOB,SPH_MAT,LUPRI)
implicit none
TYPE(TCO_ORBITAL)        :: ORB
TYPE(AOBATCH)            :: AOB
TYPE(SPHMAT),pointer     :: SPH_MAT(:)
INTEGER                  :: LUPRI
!
INTEGER                  :: np,nc,mc,na,ia

!ORB%type          = AOB%type
ORB%spherical     = AOB%spherical
ORB%center        = AOB%center
ORB%nPrimitives   = AOB%nPrimitives
ORB%maxContracted = AOB%maxContracted
ORB%maxAngmom     = AOB%maxAngmom
ORB%nAngmom       = AOB%nAngmom

np = ORB%nPrimitives
mc = ORB%maxContracted
na = ORB%nAngmom

ORB%angmom(1:na)           = AOB%angmom(1:na)
ORB%nContracted(1:na)      = AOB%nContracted(1:na)
ORB%startOrbital(1:na)     = AOB%startOrbital(1:na)
ORB%startPrimOrbital(1:na) = AOB%startPrimOrbital(1:na)
ORB%nOrbComp(1:na)         = AOB%nOrbComp(1:na)
ORB%nPrimOrbComp(1:na)     = AOB%nPrimOrbComp(1:na)
ORB%nOrbitals(1:na)        = AOB%nOrbitals(1:na)
CALL DCOPY(np,AOB%pExponents%elms,1,ORB%exponents,1)
DO ia=1,na
  nc = ORB%nContracted(ia)
  CALL DCOPY(np*nc,AOB%pCC(ia)%p%elms,1,ORB%CC(1:np,1:nc,ia),1)
ENDDO
DO ia=1,na
  ORB%SPH_MAT(ia)%p => SPH_MAT(ORB%angmom(ia)+1)
ENDDO

END SUBROUTINE set_TCO_orbital

SUBROUTINE print_TCO_integrand(P,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
INTEGER                  :: LUPRI,IPRINT
!
INTEGER                  :: I,J

CALL LSHEADER(LUPRI,'TCO_Integrand output')
WRITE(LUPRI,'(1X,A)') ''
WRITE(LUPRI,'(1X,A)') '***************************'
WRITE(LUPRI,'(1X,A)') '***    TCO_Integrand    ***'
WRITE(LUPRI,'(1X,A)') '***************************'
WRITE(LUPRI,'(1X,A)') ''

WRITE(LUPRI,'(3X,A)') '-------------- orbital 1 --------------'
CALL print_TCO_orbital(P%orbital1,LUPRI,IPRINT)
WRITE(LUPRI,'(3X,A)') '-------------- orbital 2 --------------'
CALL print_TCO_orbital(P%orbital2,LUPRI,IPRINT)
WRITE(LUPRI,'(3X,A)') '-------------- orbital 3 --------------'
CALL print_TCO_orbital(P%orbital3,LUPRI,IPRINT)

WRITE(LUPRI,'(3X,A)') '--- three center overlap  specifics ---'
!WRITE(LUPRI,'(5X,A,A55)')'Type             = ', P%type
WRITE(LUPRI,'(5X,A,L3)') 'same12AO         = ', P%same12AO
WRITE(LUPRI,'(5X,A,L3)') 'same13AO         = ', P%same13AO
WRITE(LUPRI,'(5X,A,L3)') 'same23AO         = ', P%same23AO
WRITE(LUPRI,'(5X,A,I3)') 'nPrimitives      = ', P%nPrimitives
WRITE(LUPRI,'(5X,A,I3)') 'maxContracted    = ', P%maxContracted
WRITE(LUPRI,'(5X,A,I3)') 'nAngmom          = ', P%nAngmom
WRITE(LUPRI,'(5X,A,8I3)')'No.                ', (I,I=1,P%nAngmom)
WRITE(LUPRI,'(5X,A,8I3)')'iAng1            = ', (P%iAng1(I),I=1,P%nAngmom)
WRITE(LUPRI,'(5X,A,8I3)')'iAng2            = ', (P%iAng2(I),I=1,P%nAngmom)
WRITE(LUPRI,'(5X,A,8I3)')'iAng3            = ', (P%iAng3(I),I=1,P%nAngmom)
WRITE(LUPRI,'(5X,A,8I3)')'nContracted      = ', (P%nContracted(I),I=1,P%nAngmom)
WRITE(LUPRI,'(5X,A,8I3)')'nOrbComp         = ', (P%nOrbComp(I),I=1,P%nAngmom)
WRITE(LUPRI,'(5X,A,8I3)')'nOrbitals        = ', (P%nOrbitals(I),I=1,P%nAngmom)
WRITE(LUPRI,'(3X,A)')    '----------------------------------------------------------'
WRITE(LUPRI,'(3X,A)')    '- iPrim  Int     Kfac    1/2p     iprim1  iprim2  iprim3 -'
DO I=1,P%nPrimitives
  WRITE(LUPRI,'(5X,I4,3F8.4,3I8)') I,P%integralPrefactor(I),P%preExpFactor(I),P%halfPinv(I),&
                                      P%iprim1(I),P%iprim2(I),P%iprim3(I)
ENDDO
WRITE(LUPRI,'(3X,A)')    '--------------------------------------------------------------------------------'
WRITE(LUPRI,'(3X,A)')    '- iPrim  P1x     P1y     P1z     P2x     P2y     P2z     P3x     P3y     P3z   -'
DO I=1,P%nPrimitives
  WRITE(LUPRI,'(5X,I4,9F8.4)') I,(P%vectorP1(I,J),J=1,3),(P%vectorP2(I,J),J=1,3),(P%vectorP3(I,J),J=1,3)
ENDDO
WRITE(LUPRI,'(3X,A)')    '--------------------------------------------------------------------------------'

END SUBROUTINE print_TCO_integrand

SUBROUTINE print_TCO_orbital(ORB,LUPRI,IPRINT)
implicit none
TYPE(TCO_ORBITAL)        :: ORB
INTEGER                  :: LUPRI,IPRINT
!
INTEGER                  :: I,J,K,L

!WRITE(LUPRI,'(5X,A,A55)') 'Type          = ', ORB%type
WRITE(LUPRI,'(5X,A,1L2)') 'Spherical     = ', ORB%spherical
WRITE(LUPRI,'(5X,A,3F10.6)') 'center (A)    = ', (ORB%center(I),I=1,3)
WRITE(LUPRI,'(5X,A,I3)')  'nPrimitives   = ', ORB%nPrimitives
WRITE(LUPRI,'(5X,A,I3)')  'maxContracted = ', ORB%maxContracted
WRITE(LUPRI,'(5X,A,I3)')  'maxAngmom     = ', ORB%maxAngmom
WRITE(LUPRI,'(5X,A,I3)')  'nAngmom       = ', ORB%nAngmom
WRITE(LUPRI,'(3X,A)')     '-------------------------------------------------------------------'
WRITE(LUPRI,'(3X,A)')     '-  Block  angmom  #cont.  1.orb.  1.porb    #oc  #pri.oc   #orb.  -'
DO I=1,ORB%nAngmom
  WRITE(LUPRI,'(1X,8I8)') I,ORB%angmom(I),ORB%nContracted(I),&
                            ORB%startOrbital(I),ORB%startPrimOrbital(I),&
			    ORB%nOrbComp(I),ORB%nPrimOrbComp(I),&
			    ORB%nOrbitals(I)
ENDDO
WRITE(LUPRI,'(3X,A)')     '-------------------------------------------------------------------'
WRITE(LUPRI,'(3X,A)')     '----------------- Exponents -----------------'
WRITE(LUPRI,'(5X,5F12.6)') (ORB%exponents(I),I=1,ORB%nPrimitives)
WRITE(LUPRI,'(3X,A)')     '----------------- Contraction Coefficients -----------------'
DO K=1,ORB%nAngmom
  L = ORB%angmom(K)
  WRITE(LUPRI,'(5X,A,I3,A,I3)') '* Angular block number',K,' with angular momentum', L
  DO J=1,ORB%nContracted(K)
    WRITE(LUPRI,'(2X,I3,5F12.6,/(5X,5F12.6))') J,(ORB%CC(I,J,K),I=1,ORB%nPrimitives)
  ENDDO
  IF(L .GE. 2)THEN
    WRITE(LUPRI,'(7X,A,I3)')'Spherical transformation matrix for angular momentum', L
    CALL OUTPUT(ORB%SPH_MAT(K)%p%elms, 1,2*L+1, 1,(L+1)*(L+2)/2,&
		                       2*L+1,(L+1)*(L+2)/2, 1,LUPRI)
  ENDIF
  WRITE(LUPRI,*)
ENDDO

END SUBROUTINE print_TCO_orbital

SUBROUTINE build_TCO_E0IJK(P,INTEGRAL,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
INTEGER                  :: LUPRI,IPRINT
!
INTEGER                  :: nPrim,ml1,ml2,ml3,len,idim

nPrim = P%nPrimitives
ml1   = P%orbital1%maxAngmom
ml2   = P%orbital2%maxAngmom
ml3   = P%orbital3%maxAngmom
len   = nPrim*(ml1+1)*(ml2+1)*(ml3+1)

DO idim=1,3
  CALL cart_TCO_E0IJKcoeffs(INTEGRAL%E0IJK((idim-1)*len+1:idim*len),&
                            INTEGRAL%WORK,nPrim,ml1,ml2,ml3,&
			    P%halfPinv,P%vectorP1(:,idim),&
			    P%vectorP2(:,idim),P%vectorP3(:,idim),LUPRI)
  IF(IPRINT.GT.20) THEN
    CALL print_TCO_E0IJKcoeffs(INTEGRAL%E0IJK,nPrim,ml1,ml2,ml3,idim,LUPRI)
  ENDIF
ENDDO

END SUBROUTINE build_TCO_E0IJK

SUBROUTINE cart_TCO_E0IJKcoeffs(E0IJK,E1IJK,nPrim,MAXI,MAXJ,MAXK,&
                                HPINV,PA,PB,PC,LUPRI)
implicit none
REAL(REALK)              :: E0IJK(nPrim,0:MAXI,0:MAXJ,0:MAXK)
REAL(REALK)              :: E1IJK(nPrim,0:MAXI,0:MAXJ,0:MAXK)
INTEGER                  :: nPrim,MAXI,MAXJ,MAXK
REAL(REALK)              :: HPINV(nPrim),PA(nPrim),PB(nPrim),PC(nPrim)
INTEGER                  :: LUPRI
!
REAL(REALK),parameter    :: ONE=1.d0 ,ZERO=0.d0
INTEGER                  :: iPrim,I,J,K,IM1,JM1,KM1

DO iPrim=1,nPrim
  E0IJK(iPrim,0,0,0) = ONE
ENDDO
DO iPrim=1,nPrim
  E1IJK(iPrim,0,0,0) = ZERO
ENDDO


DO I=1,MAXI
  IM1=I-1
  DO iPrim=1,nPrim
    E0IJK(iPrim,I,0,0)=PA(iPrim)*E0IJK(iPrim,IM1,0,0)+E1IJK(iPrim,IM1,0,0)
  ENDDO
  DO iPrim=1,nPrim
    E1IJK(iPrim,I,0,0)=HPINV(iPrim)*I*E0IJK(iPrim,IM1,0,0)
  ENDDO
ENDDO

DO J=1,MAXJ
  JM1=J-1
  DO iPrim=1,nPrim
    E0IJK(iPrim,0,J,0)=PB(iPrim)*E0IJK(iPrim,0,JM1,0)+E1IJK(iPrim,0,JM1,0)
  ENDDO
  DO iPrim=1,nPrim
    E1IJK(iPrim,0,J,0)=HPINV(iPrim)*J*E0IJK(iPrim,0,JM1,0)
  ENDDO
ENDDO

DO K=1,MAXK
  KM1=K-1
  DO iPrim=1,nPrim
    E0IJK(iPrim,0,0,K)=PC(iPrim)*E0IJK(iPrim,0,0,KM1)+E1IJK(iPrim,0,0,KM1)
  ENDDO
  DO iPrim=1,nPrim
    E1IJK(iPrim,0,0,K)=HPINV(iPrim)*K*E0IJK(iPrim,0,0,KM1)
  ENDDO
ENDDO


IF(MAXI.GT.0) THEN
  DO J=1,MAXJ
    JM1=J-1
    DO I=1,MAXI
      IM1=I-1
      DO iPrim=1,nPrim
        E0IJK(iPrim,I,J,0)=PA(iPrim)*E0IJK(iPrim,IM1,J,0)+E1IJK(iPrim,IM1,J,0)
      ENDDO
      DO iPrim=1,nPrim
        E1IJK(iPrim,I,J,0)=HPINV(iPrim)*(I*E0IJK(iPrim,IM1,J,0)+J*E0IJK(iPrim,I,JM1,0))
      ENDDO
    ENDDO
  ENDDO

  DO K=1,MAXK
    KM1=K-1
    DO I=1,MAXI
      IM1=I-1
      DO iPrim=1,nPrim
        E0IJK(iPrim,I,0,K)=PA(iPrim)*E0IJK(iPrim,IM1,0,K)+E1IJK(iPrim,IM1,0,K)
      ENDDO
      DO iPrim=1,nPrim
        E1IJK(iPrim,I,0,K)=HPINV(iPrim)*(I*E0IJK(iPrim,IM1,0,K)+K*E0IJK(iPrim,I,0,KM1))
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF(MAXJ.GT.0) THEN
  DO K=1,MAXK
    KM1=K-1
    DO J=1,MAXJ
      JM1=J-1
      DO iPrim=1,nPrim
        E0IJK(iPrim,0,J,K)=PB(iPrim)*E0IJK(iPrim,0,JM1,K)+E1IJK(iPrim,0,JM1,K)
      ENDDO
      DO iPrim=1,nPrim
        E1IJK(iPrim,0,J,K)=HPINV(iPrim)*(J*E0IJK(iPrim,0,JM1,K)+K*E0IJK(iPrim,0,J,KM1))
      ENDDO
    ENDDO
  ENDDO
ENDIF


IF(MAXI.GT.0.AND.MAXJ.GT.0) THEN
  DO K=1,MAXK
    KM1=K-1
    DO J=1,MAXJ
      JM1=J-1
      DO I=1,MAXI
        IM1=I-1
	DO iPrim=1,nPrim
          E0IJK(iPrim,I,J,K)=PA(iPrim)*E0IJK(iPrim,IM1,J,K)+E1IJK(iPrim,IM1,J,K)
	ENDDO
	DO iPrim=1,nPrim
          E1IJK(iPrim,I,J,K)=HPINV(iPrim)*(I*E0IJK(iPrim,IM1,J,K)+J*E0IJK(iPrim,I,JM1,K)+K*E0IJK(iPrim,I,J,KM1))
	ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE cart_TCO_E0IJKcoeffs

SUBROUTINE print_TCO_E0IJKcoeffs(E0IJK,nPrim,MAXI,MAXJ,MAXK,idim,LUPRI)
implicit none
REAL(REALK)              :: E0IJK(nPrim,0:MAXI,0:MAXJ,0:MAXK,3)
INTEGER                  :: nPrim,MAXI,MAXJ,MAXK,idim,LUPRI
!
CHARACTER(7)             :: WORD
INTEGER                  :: I,J,K,iPrim

WRITE(LUPRI,*)'Output from TCO E0_IJK coeffs'
WRITE(LUPRI,*)'-----------------------------'
WRITE (LUPRI,'(2X,A,3I5)') 'MAXI,MAXJ,MAXK   ',MAXI,MAXJ,MAXK

IF(idim.EQ.1) WORD = 'E0IJK_X'
IF(idim.EQ.2) WORD = 'E0IJK_Y'
IF(idim.EQ.3) WORD = 'E0IJK_Z'

DO I=0,MAXI 
  DO J=0,MAXJ 
    DO K=0,MAXK
      WRITE(LUPRI,'(/,2X,A7,A1,I1,A1,I1,A1,I1,A1)') WORD,'(',I,',',J, ',',K,')'
      WRITE(LUPRI,'(1X,6F12.8)') (E0IJK(iPrim,I,J,K,idim),iPrim=1,nPrim)
    ENDDO
  ENDDO
ENDDO
WRITE(LUPRI,*)

END SUBROUTINE print_TCO_E0IJKcoeffs

SUBROUTINE print_TCO_matrix(A,n1,n2,WORD,order,LUPRI)
implicit none
REAL(REALK)              :: A(n1,n2)
CHARACTER(5)             :: WORD
INTEGER                  :: n1,n2,order,LUPRI
!
INTEGER                  :: i1,i2

IF(order.eq.1) THEN
  DO i2=1,n2
    WRITE(LUPRI,'(/,2X,A5,A1,i4,A1)') WORD,'(',i2,')'
    WRITE(LUPRI,'(1X,6F12.8)') (A(i1,i2),i1=1,n1)
  ENDDO
  WRITE(LUPRI,*)
ELSEIF(order.eq.2) THEN
  DO i1=1,n1
    WRITE(LUPRI,'(/,2X,A5,A1,i4,A1)') WORD,'(',i1,')'
    WRITE(LUPRI,'(1X,6F12.8)') (A(i1,i2),i2=1,n2)
  ENDDO
  WRITE(LUPRI,*)
ELSE
  WRITE(LUPRI,'(A)') 'WARNING!! Wrong option in print_TCO_matrix'
ENDIF

END SUBROUTINE print_TCO_matrix

SUBROUTINE build_TCO_E0ABC(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: nPrim,ml1,ml2,ml3,l1,l2,l3,ijk1,ijk2,ijk3,ijk

nPrim = P%nPrimitives

ml1 = P%orbital1%maxAngmom
ml2 = P%orbital2%maxAngmom
ml3 = P%orbital3%maxAngmom
l1  = P%orbital1%angmom(P%iAng1(iAngmom))
l2  = P%orbital2%angmom(P%iAng2(iAngmom))
l3  = P%orbital3%angmom(P%iAng3(iAngmom))
ijk1 = (l1+1)*(l1+2)/2
ijk2 = (l2+1)*(l2+2)/2
ijk3 = (l3+1)*(l3+2)/2
ijk  = ijk1*ijk2*ijk3

CALL build_TCO_E0ABC1(INTEGRAL%E0IJK,INTEGRAL%IN,P%preExpFactor,&
                      nPrim,ijk,l1,l2,l3,ml1,ml2,ml3,LUPRI)
INTEGRAL%size1 = nPrim
INTEGRAL%size2 = ijk

IF(IPRINT.GT.20) THEN
  WRITE(LUPRI,*)'Output from TCO E0_ABC coeffs'
  WRITE(LUPRI,*)'-----------------------------'
  WRITE(LUPRI,'(2X,A,3I5)') 'l1,l2,l3   ',l1,l2,l3
  CALL print_TCO_matrix(INTEGRAL%IN,nPrim,ijk,'E0ABC',1,LUPRI)
ENDIF

END SUBROUTINE build_TCO_E0ABC

SUBROUTINE build_TCO_E0ABC1(E0IJK,E0ABC,preExpFactor,nPrim,nijk,&
                            l1,l2,l3,MAXI,MAXJ,MAXK,LUPRI)
implicit none
REAL(REALK)              :: E0IJK(nPrim,0:MAXI,0:MAXJ,0:MAXK,3)
REAL(REALK)              :: E0ABC(nPrim,nijk)
REAL(REALK)              :: preExpFactor(nPrim)
INTEGER                  :: nPrim,nijk,l1,l2,l3,MAXI,MAXJ,MAXK,LUPRI
!
INTEGER                  :: ijk,i1,j1,k1,i2,j2,k2,i3,j3,k3,iPrim

ijk=0
DO i3=l3,0,-1
DO j3=l3-i3,0,-1
k3=l3-i3-j3
  DO i2=l2,0,-1
  DO j2=l2-i2,0,-1
  k2=l2-i2-j2
    DO i1=l1,0,-1
    DO j1=l1-i1,0,-1
    k1=l1-i1-j1
      ijk=ijk+1
      DO iPrim=1,nPrim
        E0ABC(iPrim,ijk) = preExpFactor(iPrim)*E0IJK(iPrim,i1,i2,i3,1)
      ENDDO
      DO iPrim=1,nPrim
        E0ABC(iPrim,ijk) = E0ABC(iPrim,ijk)*E0IJK(iPrim,j1,j2,j3,2)
      ENDDO
      DO iPrim=1,nPrim
        E0ABC(iPrim,ijk) = E0ABC(iPrim,ijk)*E0IJK(iPrim,k1,k2,k3,3)
      ENDDO
    ENDDO
    ENDDO
  ENDDO
  ENDDO
ENDDO
ENDDO
IF(ijk.NE.nijk) THEN
  CALL LSQUIT('Programming error, wrong number of terms in E0IJK->E0ABC',lupri)
ENDIF

END SUBROUTINE build_TCO_E0ABC1

SUBROUTINE build_TCO_integral(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: nPrim,l1,l2,l3,ijk1,ijk2,ijk3,ijk

nPrim = P%nPrimitives

l1  = P%orbital1%angmom(P%iAng1(iAngmom))
l2  = P%orbital2%angmom(P%iAng2(iAngmom))
l3  = P%orbital3%angmom(P%iAng3(iAngmom))
ijk1 = (l1+1)*(l1+2)/2
ijk2 = (l2+1)*(l2+2)/2
ijk3 = (l3+1)*(l3+2)/2
ijk  = ijk1*ijk2*ijk3

IF(nPrim.NE.INTEGRAL%size1 .OR. ijk.NE.INTEGRAL%size2) THEN
  CALL LSQUIT('Programming error, wrong dimensions in build integral',lupri)
ENDIF

CALL build_TCO_integral1(INTEGRAL%IN,INTEGRAL%OUT,P%integralPrefactor,&
                         nPrim,ijk,LUPRI)
CALL swap_realPointers(INTEGRAL%IN,INTEGRAL%OUT)

IF(IPRINT.GT.20) THEN
  WRITE(LUPRI,*)'Output from TCO uncontracted integrals'
  WRITE(LUPRI,*)'--------------------------------------'
  WRITE(LUPRI,'(2X,A,3I5)') 'l1,l2,l3   ',l1,l2,l3
  CALL print_TCO_matrix(INTEGRAL%IN,nPrim,ijk,'INT_U',1,LUPRI)
ENDIF

END SUBROUTINE build_TCO_integral

SUBROUTINE build_TCO_integral1(E0ABC,INT,integralPrefactor,nPrim,nijk,LUPRI)
implicit none
REAL(REALK)              :: E0ABC(nPrim,nijk)
REAL(REALK)              :: INT(nPrim,nijk)
REAL(REALK)              :: integralPrefactor(nPrim)
INTEGER                  :: nPrim,nijk,LUPRI
!
INTEGER                  :: iPrim,ijk

DO ijk=1,nijk
  DO iPrim=1,nPrim
    INT(iPrim,ijk)=E0ABC(iPrim,ijk)*integralPrefactor(iPrim)
  ENDDO
ENDDO

END SUBROUTINE build_TCO_integral1

SUBROUTINE TCO_contractBasis(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: nPrim,nCont,ijk,iA1,iA2,iA3,nC1,nC2,nC3

nPrim = P%nPrimitives
iA1 = P%iAng1(iAngmom)
iA2 = P%iAng2(iAngmom)
iA3 = P%iAng3(iAngmom)
nC1 = P%orbital1%nContracted(iA1)
nC2 = P%orbital2%nContracted(iA2)
nC3 = P%orbital3%nContracted(iA3)
nCont = nC1*nC2*nC3
ijk   = INTEGRAL%size2

CALL TCO_contractBasis1(INTEGRAL%IN,INTEGRAL%OUT,INTEGRAL%WORK,P,&
                        nPrim,nCont,ijk,iA1,iA2,iA3,nC1,nC2,nC3,LUPRI,IPRINT)
CALL swap_realPointers(INTEGRAL%IN,INTEGRAL%OUT)
INTEGRAL%size1 = nCont

IF(IPRINT.GT.20) THEN
  WRITE(LUPRI,*)'Output from TCO contracted integrals'
  WRITE(LUPRI,*)'------------------------------------'
  WRITE(LUPRI,'(2X,A,3I5)') 'l1,l2,l3   ',&
        P%orbital1%angmom(iA1),P%orbital2%angmom(iA2),P%orbital3%angmom(iA3)
  CALL print_TCO_matrix(INTEGRAL%IN,nCont,ijk,'INT_C',1,LUPRI)
ENDIF

END SUBROUTINE TCO_contractBasis

SUBROUTINE TCO_contractBasis1(primINT,contINT,CC,P,nPrim,nCont,ijk,&
                              iA1,iA2,iA3,nC1,nC2,nC3,LUPRI,IPRINT)
implicit none
REAL(REALK)              :: primINT(nPrim,ijk),contINT(nCont,ijk),CC(nPrim,nCont)
TYPE(TCO_INTEGRAND)      :: P
INTEGER                  :: nPrim,nCont,ijk,iA1,iA2,iA3,nC1,nC2,nC3,LUPRI,IPRINT

CALL TCO_constructContraction(CC,P,nPrim,iA1,iA2,iA3,nC1,nC2,nC3,LUPRI,IPRINT)
CALL DGEMM('T','N',nCont,ijk,nPrim,&
           1.d0,CC,nPrim,primINT,nPrim,0.d0,contINT,nCont)

END SUBROUTINE TCO_contractBasis1

SUBROUTINE TCO_constructContraction(CC,P,nPrim,iA1,iA2,iA3,nC1,nC2,nC3,LUPRI,IPRINT)
implicit none
REAL(REALK)              :: CC(nPrim,nC1,nC2,nC3)
TYPE(TCO_INTEGRAND)      :: P
INTEGER                  :: nPrim,iA1,iA2,iA3,nC1,nC2,nC3,LUPRI,IPRINT
!
INTEGER                  :: iC1,iC2,iC3,iPrim,iP1,iP2,iP3

IF(P%same12AO.AND.P%same13AO.AND.P%same23AO) THEN              ! T T T
  DO iC3=1,nC3
    DO iC2=1,nC2
      DO iC1=1,nC1
        DO iPrim=1,nPrim
          iP1 = P%iPrim1(iPrim)
          iP2 = P%iPrim2(iPrim)
          iP3 = P%iPrim3(iPrim)
	  IF(iP1.NE.iP2) THEN
	    IF(iP2.NE.iP3) THEN
	      CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                              P%orbital2%CC(iP2,iC2,iA2)*&
		                      P%orbital3%CC(iP3,iC3,iA3)+&
				      P%orbital1%CC(iP3,iC1,iA1)*&
	                              P%orbital2%CC(iP1,iC2,iA2)*&
		                      P%orbital3%CC(iP2,iC3,iA3)+&
				      P%orbital1%CC(iP2,iC1,iA1)*&
	                              P%orbital2%CC(iP3,iC2,iA2)*&
		                      P%orbital3%CC(iP1,iC3,iA3)+&
				      P%orbital1%CC(iP1,iC1,iA1)*&
	                              P%orbital2%CC(iP3,iC2,iA2)*&
		                      P%orbital3%CC(iP2,iC3,iA3)+&
				      P%orbital1%CC(iP2,iC1,iA1)*&
	                              P%orbital2%CC(iP1,iC2,iA2)*&
		                      P%orbital3%CC(iP3,iC3,iA3)+&
				      P%orbital1%CC(iP3,iC1,iA1)*&
	                              P%orbital2%CC(iP2,iC2,iA2)*&
		                      P%orbital3%CC(iP1,iC3,iA3)
	    ELSE
	      CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                              P%orbital2%CC(iP2,iC2,iA2)*&
		                      P%orbital3%CC(iP3,iC3,iA3)+&
				      P%orbital1%CC(iP2,iC1,iA1)*&
	                              P%orbital2%CC(iP1,iC2,iA2)*&
		                      P%orbital3%CC(iP3,iC3,iA3)+&
				      P%orbital1%CC(iP2,iC1,iA1)*&
	                              P%orbital2%CC(iP3,iC2,iA2)*&
		                      P%orbital3%CC(iP1,iC3,iA3)
	    ENDIF
	  ELSE
	    IF(iP2.NE.iP3) THEN
	      CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                              P%orbital2%CC(iP2,iC2,iA2)*&
		                      P%orbital3%CC(iP3,iC3,iA3)+&
				      P%orbital1%CC(iP1,iC1,iA1)*&
	                              P%orbital2%CC(iP3,iC2,iA2)*&
		                      P%orbital3%CC(iP2,iC3,iA3)+&
				      P%orbital1%CC(iP3,iC1,iA1)*&
	                              P%orbital2%CC(iP1,iC2,iA2)*&
		                      P%orbital3%CC(iP2,iC3,iA3)
	    ELSE
	      CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                              P%orbital2%CC(iP2,iC2,iA2)*&
		                      P%orbital3%CC(iP3,iC3,iA3)
            ENDIF
	  ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(P%same12AO.AND.(P%same13AO.EQV.P%same23AO)) THEN        ! T F F
  DO iC3=1,nC3
    DO iC2=1,nC2
      DO iC1=1,nC1
        DO iPrim=1,nPrim
          iP1 = P%iPrim1(iPrim)
          iP2 = P%iPrim2(iPrim)
          iP3 = P%iPrim3(iPrim)
	  IF(iP1.NE.iP2) THEN
	    CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                            P%orbital2%CC(iP2,iC2,iA2)*&
		                    P%orbital3%CC(iP3,iC3,iA3)+&
				    P%orbital1%CC(iP2,iC1,iA1)*&
	                            P%orbital2%CC(iP1,iC2,iA2)*&
		                    P%orbital3%CC(iP3,iC3,iA3)
	  ELSE
	    CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                            P%orbital2%CC(iP2,iC2,iA2)*&
		                    P%orbital3%CC(iP3,iC3,iA3)
	  ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(P%same13AO.AND.(P%same12AO.EQV.P%same23AO)) THEN        ! F T F
  DO iC3=1,nC3
    DO iC2=1,nC2
      DO iC1=1,nC1
        DO iPrim=1,nPrim
          iP1 = P%iPrim1(iPrim)
          iP2 = P%iPrim2(iPrim)
          iP3 = P%iPrim3(iPrim)
	  IF(iP1.NE.iP3) THEN
	    CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                            P%orbital2%CC(iP2,iC2,iA2)*&
		                    P%orbital3%CC(iP3,iC3,iA3)+&
				    P%orbital1%CC(iP3,iC1,iA1)*&
	                            P%orbital2%CC(iP2,iC2,iA2)*&
		                    P%orbital3%CC(iP1,iC3,iA3)
	  ELSE
	    CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                            P%orbital2%CC(iP2,iC2,iA2)*&
		                    P%orbital3%CC(iP3,iC3,iA3)
	  ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(P%same23AO.AND.(P%same12AO.EQV.P%same13AO)) THEN        ! F F T
  DO iC3=1,nC3
    DO iC2=1,nC2
      DO iC1=1,nC1
        DO iPrim=1,nPrim
          iP1 = P%iPrim1(iPrim)
          iP2 = P%iPrim2(iPrim)
          iP3 = P%iPrim3(iPrim)
	  IF(iP2.NE.iP3) THEN
	    CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                            P%orbital2%CC(iP2,iC2,iA2)*&
		                    P%orbital3%CC(iP3,iC3,iA3)+&
				    P%orbital1%CC(iP1,iC1,iA1)*&
	                            P%orbital2%CC(iP3,iC2,iA2)*&
		                    P%orbital3%CC(iP2,iC3,iA3)
	  ELSE
	    CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                            P%orbital2%CC(iP2,iC2,iA2)*&
		                    P%orbital3%CC(iP3,iC3,iA3)
	  ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(.NOT.(P%same12AO.OR.P%same13AO.OR.P%same23AO)) THEN     ! F F F
  DO iC3=1,nC3
    DO iC2=1,nC2
      DO iC1=1,nC1
        DO iPrim=1,nPrim
          iP1 = P%iPrim1(iPrim)
          iP2 = P%iPrim2(iPrim)
          iP3 = P%iPrim3(iPrim)
	  CC(iPrim,iC1,iC2,iC3) = P%orbital1%CC(iP1,iC1,iA1)*&
	                          P%orbital2%CC(iP2,iC2,iA2)*&
		                  P%orbital3%CC(iP3,iC3,iA3)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE                                               ! T T F,  T F T,  F T T
  CALL LSQUIT('Programming error, problem with determining same AObatches pairs',lupri)
ENDIF

IF(IPRINT.GT.30) THEN
  WRITE(LUPRI,'(5X,A)') '***************************************************'
  WRITE(LUPRI,'(5X,A)') '***                 Overlap  CC                 ***'
  WRITE(LUPRI,'(5X,A)') '***************************************************'
  WRITE(LUPRI,'(1X,A)') '  C1  C2  C3'
  DO iC1=1,nC1
    DO iC2=1,nC2
      DO iC3=1,nC3
        WRITE(LUPRI,'(1X,3I4,6F10.4/,(13X,6F10.4))')&
	            iC1,iC2,iC3,(CC(iPrim,iC1,iC2,iC3),iPrim=1,nPrim)
      ENDDO
    ENDDO
  ENDDO
  WRITE(LUPRI,*)
ENDIF

END SUBROUTINE TCO_constructContraction

SUBROUTINE TCO_transposition(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: nCont,ijk

nCont = INTEGRAL%size1
ijk   = INTEGRAL%size2

CALL transpose_matrix(INTEGRAL%IN,INTEGRAL%OUT,nCont,ijk)
CALL swap_realPointers(INTEGRAL%IN,INTEGRAL%OUT)
INTEGRAL%size1 = ijk
INTEGRAL%size2 = nCont

IF(IPRINT.GT.200) THEN
  WRITE(LUPRI,*)'Output from TCO contracted integrals (transposed)'
  WRITE(LUPRI,*)'-------------------------------------------------'
  WRITE(LUPRI,'(2X,A,3I5)') 'l1,l2,l3   ',P%orbital1%angmom(P%iAng1(iAngmom)),&
        P%orbital2%angmom(P%iAng2(iAngmom)),P%orbital3%angmom(P%iAng3(iAngmom))
  CALL print_TCO_matrix(INTEGRAL%IN,ijk,nCont,'INT_C',2,LUPRI)
ENDIF

END SUBROUTINE TCO_transposition

SUBROUTINE TCO_sphericalTransform(P,INTEGRAL,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: nCont,iA1,iA2,iA3,l1,l2,l3
LOGICAL                  :: sph1,sph2,sph3
INTEGER                  :: ijk1,ijk2,ijk3,ijk,lm1,lm2,lm3,lm

iA1 = P%iAng1(iAngmom)
iA2 = P%iAng2(iAngmom)
iA3 = P%iAng3(iAngmom)
l1  = P%orbital1%angmom(iA1)
l2  = P%orbital2%angmom(iA2)
l3  = P%orbital3%angmom(iA3)
sph1 = P%orbital1%spherical.AND.(l1.GT.1)
sph2 = P%orbital2%spherical.AND.(l2.GT.1)
sph3 = P%orbital3%spherical.AND.(l3.GT.1)

IF(sph1.OR.sph2.OR.sph3) THEN
  ijk1 = (l1+1)*(l1+2)/2
  ijk2 = (l2+1)*(l2+2)/2
  ijk3 = (l3+1)*(l3+2)/2
  ijk  = ijk1*ijk2*ijk3
  lm1  = 2*l1+1
  lm2  = 2*l2+1
  lm3  = 2*l3+1
  lm   = lm1*lm2*lm3

  nCont = INTEGRAL%size2

  IF(ijk.NE.INTEGRAL%size1) THEN
    CALL LSQUIT('Programming error, wrong dimensions in build integral',lupri)
  ENDIF

  CALL TCO_sphericalTransform1(INTEGRAL%IN,INTEGRAL%OUT,INTEGRAL%WORK,P,&
                               ijk,lm,nCont,iA1,iA2,iA3,lm1,lm2,lm3,&
			       ijk1,ijk2,ijk3,LUPRI,IPRINT)
  CALL swap_realPointers(INTEGRAL%IN,INTEGRAL%OUT)
  INTEGRAL%size1 = lm
ENDIF

IF(IPRINT.GT.20) THEN
  WRITE(LUPRI,*)'Output from TCO contracted spherical integrals'
  WRITE(LUPRI,*)'----------------------------------------------'
  WRITE(LUPRI,'(2X,A,3I5)') 'l1,l2,l3   ',l1,l2,l3
! size1 = lm,  size2 = nCont
  CALL print_TCO_matrix(INTEGRAL%IN,INTEGRAL%size1,INTEGRAL%size2,&
                        'INTCS',2,LUPRI)
ENDIF

END SUBROUTINE TCO_sphericalTransform

SUBROUTINE TCO_sphericalTransform1(cartINT,sphINT,spherical,P,ijk,lm,nCont,&
                                   iA1,iA2,iA3,lm1,lm2,lm3,ijk1,ijk2,ijk3,&
				   LUPRI,IPRINT)
implicit none
REAL(REALK)              :: cartINT(ijk,nCont),sphINT(lm,nCont)
REAL(REALK)              :: spherical(ijk,lm)
TYPE(TCO_INTEGRAND)      :: P
INTEGER                  :: ijk,lm,nCont
INTEGER                  :: iA1,iA2,iA3,lm1,lm2,lm3,ijk1,ijk2,ijk3
INTEGER                  :: LUPRI,IPRINT

CALL TCO_constructSphericalTransformation(spherical,&
                           P%orbital1%SPH_MAT(iA1)%p%elms,lm1,ijk1,&
	                   P%orbital2%SPH_MAT(iA2)%p%elms,lm2,ijk2,&
		           P%orbital3%SPH_MAT(iA3)%p%elms,lm3,ijk3,LUPRI,IPRINT)
CALL DGEMM('T','N',lm,nCont,ijk,&
           1.d0,spherical,ijk,cartINT,ijk,0.d0,sphINT,lm)

END SUBROUTINE TCO_sphericalTransform1

SUBROUTINE TCO_constructSphericalTransformation(spherical,&
               spher1,lm1,ijk1,spher2,lm2,ijk2,spher3,lm3,ijk3,LUPRI,IPRINT)
implicit none
REAL(REALK)              :: spherical(ijk1,ijk2,ijk3,lm1,lm2,lm3)
REAL(REALK)              :: spher1(lm1,ijk1),spher2(lm2,ijk2),spher3(lm3,ijk3)
INTEGER                  :: lm1,lm2,lm3,ijk1,ijk2,ijk3
INTEGER                  :: LUPRI,IPRINT
!
INTEGER                  :: i_ijk1,i_ijk2,i_ijk3,i_lm1,i_lm2,i_lm3
REAL(REALK)              :: temp1,temp2,temp3

DO i_lm3=1,lm3
  DO i_lm2=1,lm2
    DO i_lm1=1,lm1
      DO i_ijk3=1,ijk3
        temp3 = spher3(i_lm3,i_ijk3)
        DO i_ijk2=1,ijk2
	  temp2 = spher2(i_lm2,i_ijk2)*temp3
	  DO i_ijk1=1,ijk1
	    temp1 = spher1(i_lm1,i_ijk1)*temp2
	    spherical(i_ijk1,i_ijk2,i_ijk3,i_lm1,i_lm2,i_lm3) = temp1
	  ENDDO
	ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO

IF(IPRINT.GT.50) THEN
  WRITE(LUPRI,'(5X,A)') '***************************************************'
  WRITE(LUPRI,'(5X,A)') '***          Spherical  Transformation          ***'
  WRITE(LUPRI,'(5X,A)') '***************************************************'
  DO i_lm3=1,lm3
    DO i_lm2=1,lm2
      DO i_lm1=1,lm1
        WRITE(LUPRI,'(5X,3(A,I3))') 'lm1 =',i_lm1,' lm2 =',i_lm2,' lm3 =',i_lm3
	DO i_ijk1=1,ijk1
	  DO i_ijk2=1,ijk2
            WRITE(LUPRI,'(1X,2I4,6F10.4/,(9X,6F10.4))') i_ijk1,i_ijk2,&
	    (spherical(i_ijk1,i_ijk2,i_ijk3,i_lm1,i_lm2,i_lm3),i_ijk3=1,ijk3)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  WRITE(LUPRI,*)
ENDIF

END SUBROUTINE TCO_constructSphericalTransformation

SUBROUTINE TCO_distributeIntegrals(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
INTEGER                  :: iAngmom,LUPRI,IPRINT

IF(INPUT%same12AOs.AND.INPUT%same13AOs.AND.INPUT%same23AOs) THEN         ! T T T
  CALL DistributeNormal1(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
ELSEIF(INPUT%same12AOs.AND.(INPUT%same13AOs.EQV.INPUT%same23AOs)) THEN   ! T F F
  CALL DistributeNormal2(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
ELSEIF(INPUT%same13AOs.AND.(INPUT%same12AOs.EQV.INPUT%same23AOs)) THEN   ! F T F
  CALL DistributeNormal3(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
ELSEIF(INPUT%same23AOs.AND.(INPUT%same12AOs.EQV.INPUT%same13AOs)) THEN   ! F F T
  CALL DistributeNormal4(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
ELSEIF(.NOT.(INPUT%same12AOs.OR.INPUT%same13AOs.OR.INPUT%same23AOs)) THEN! F F F
  CALL DistributeNormal5(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
ELSE                                               ! T T F,  T F T,  F T T
  CALL LSQUIT('Programming error, impossible combination of same_xy_AOs',lupri)
ENDIF

END SUBROUTINE TCO_distributeIntegrals

SUBROUTINE DistributeNormal1(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: iA1,iA2,iA3
INTEGER                  :: nAng1,nAng2,nAng3,nCont1,nCont2,nCont3
INTEGER                  :: iAng1,iAng2,iAng3,iCont1,iCont2,iCont3
INTEGER                  :: start1,start2,start3,n1,n2,n3,i1,i2,i3
LOGICAL                  :: same12,same23
INTEGER                  :: int123
REAL(REALK)              :: v2,v3,m23,value,value1,value2,value3,value12,value13

iA1 = P%iAng1(iAngmom)
iA2 = P%iAng2(iAngmom)
iA3 = P%iAng3(iAngmom)
nAng1 = P%orbital1%nOrbComp(iA1)
nAng2 = P%orbital2%nOrbComp(iA2)
nAng3 = P%orbital3%nOrbComp(iA3)
nCont1 = P%orbital1%nContracted(iA1)
nCont2 = P%orbital2%nContracted(iA2)
nCont3 = P%orbital3%nContracted(iA3)
start1 = P%orbital1%startOrbital(iA1)
start2 = P%orbital2%startOrbital(iA2)
start3 = P%orbital3%startOrbital(iA3)
same12 = (start1.EQ.start2)
same23 = (start2.EQ.start3)

IF(INPUT%CONT_VECTOR) THEN

IF(same12.AND.same23) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+&
	                                INTEGRAL%IN(int123)*v3
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(same12) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    v2=INPUT%VECTOR(i2)
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      value2=value*v2
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+value*v3
	      OUTPUT%ResultMat(i1,i3,1)=OUTPUT%ResultMat(i1,i3,1)+value2
	      OUTPUT%ResultMat(i3,i1,1)=OUTPUT%ResultMat(i3,i1,1)+value2
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(same23) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      value3=value*v3
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+value3
	      OUTPUT%ResultMat(i2,i1,1)=OUTPUT%ResultMat(i2,i1,1)+value3
	      OUTPUT%ResultMat(i2,i3,1)=OUTPUT%ResultMat(i2,i3,1)+&
	                                value*INPUT%VECTOR(i1)
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    v2=INPUT%VECTOR(i2)
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      value1=value*INPUT%VECTOR(i1)
	      value2=value*v2
	      value3=value*v3
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+value3
	      OUTPUT%ResultMat(i3,i1,1)=OUTPUT%ResultMat(i3,i1,1)+value2
	      OUTPUT%ResultMat(i2,i3,1)=OUTPUT%ResultMat(i2,i3,1)+value1
	      OUTPUT%ResultMat(i1,i3,1)=OUTPUT%ResultMat(i1,i3,1)+value2
	      OUTPUT%ResultMat(i2,i1,1)=OUTPUT%ResultMat(i2,i1,1)+value3
	      OUTPUT%ResultMat(i3,i2,1)=OUTPUT%ResultMat(i3,i2,1)+value1
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ELSEIF(INPUT%CONT_MATRIX) THEN

IF(same12.AND.same23) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+&
	                               INTEGRAL%IN(int123)*m23
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(same12) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)+INPUT%MATRIX(i3,i2)
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+value*m23
	      OUTPUT%ResultMat(i3,1,1)=OUTPUT%ResultMat(i3,1,1)+&
	                               value*INPUT%MATRIX(i1,i2)
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(same23) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      value13=value*(INPUT%MATRIX(i1,i3)+INPUT%MATRIX(i3,i1))
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+value*m23
	      OUTPUT%ResultMat(i2,1,1)=OUTPUT%ResultMat(i2,1,1)+value13
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)+INPUT%MATRIX(i3,i2)
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      value12=value*(INPUT%MATRIX(i1,i2)+INPUT%MATRIX(i2,i1))
	      value13=value*(INPUT%MATRIX(i1,i3)+INPUT%MATRIX(i3,i1))
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+value*m23
	      OUTPUT%ResultMat(i3,1,1)=OUTPUT%ResultMat(i3,1,1)+value12
	      OUTPUT%ResultMat(i2,1,1)=OUTPUT%ResultMat(i2,1,1)+value13
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ELSE

IF(same12.AND.same23) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,i2,i3)=value
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(same12) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,i2,i3)=value
	      OUTPUT%ResultMat(i1,i3,i2)=value
	      OUTPUT%ResultMat(i3,i1,i2)=value
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(same23) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,i2,i3)=value
	      OUTPUT%ResultMat(i2,i1,i3)=value
	      OUTPUT%ResultMat(i2,i3,i1)=value
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,i2,i3)=value
	      OUTPUT%ResultMat(i3,i1,i2)=value
	      OUTPUT%ResultMat(i2,i3,i1)=value
	      OUTPUT%ResultMat(i1,i3,i2)=value
	      OUTPUT%ResultMat(i2,i1,i3)=value
	      OUTPUT%ResultMat(i3,i2,i1)=value
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ENDIF

END SUBROUTINE DistributeNormal1

SUBROUTINE DistributeNormal2(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: iA1,iA2,iA3
INTEGER                  :: nAng1,nAng2,nAng3,nCont1,nCont2,nCont3
INTEGER                  :: iAng1,iAng2,iAng3,iCont1,iCont2,iCont3
INTEGER                  :: start1,start2,start3,n1,n2,n3,i1,i2,i3
LOGICAL                  :: same12
INTEGER                  :: int123
REAL(REALK)              :: v3,m23,value,value3

iA1 = P%iAng1(iAngmom)
iA2 = P%iAng2(iAngmom)
iA3 = P%iAng3(iAngmom)
nAng1 = P%orbital1%nOrbComp(iA1)
nAng2 = P%orbital2%nOrbComp(iA2)
nAng3 = P%orbital3%nOrbComp(iA3)
nCont1 = P%orbital1%nContracted(iA1)
nCont2 = P%orbital2%nContracted(iA2)
nCont3 = P%orbital3%nContracted(iA3)
start1 = P%orbital1%startOrbital(iA1)
start2 = P%orbital2%startOrbital(iA2)
start3 = P%orbital3%startOrbital(iA3)
same12 = (start1.EQ.start2)

IF(INPUT%CONT_VECTOR) THEN

IF(same12) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+&
	                                INTEGRAL%IN(int123)*v3
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value3=INTEGRAL%IN(int123)*v3
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+value3
	      OUTPUT%ResultMat(i2,i1,1)=OUTPUT%ResultMat(i2,i1,1)+value3
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ELSEIF(INPUT%CONT_MATRIX) THEN

IF(same12) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+&
	                               INTEGRAL%IN(int123)*m23
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+value*m23
	      OUTPUT%ResultMat(i2,1,1)=OUTPUT%ResultMat(i2,1,1)+&
	                               value*INPUT%MATRIX(i1,i3)
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ELSE

IF(same12) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,i2,i3)=INTEGRAL%IN(int123)
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,i2,i3)=value
	      OUTPUT%ResultMat(i2,i1,i3)=value
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ENDIF

END SUBROUTINE DistributeNormal2

SUBROUTINE DistributeNormal3(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: iA1,iA2,iA3
INTEGER                  :: nAng1,nAng2,nAng3,nCont1,nCont2,nCont3
INTEGER                  :: iAng1,iAng2,iAng3,iCont1,iCont2,iCont3
INTEGER                  :: start1,start2,start3,n1,n2,n3,i1,i2,i3
LOGICAL                  :: same13
INTEGER                  :: int123
REAL(REALK)              :: v3,m23,value

iA1 = P%iAng1(iAngmom)
iA2 = P%iAng2(iAngmom)
iA3 = P%iAng3(iAngmom)
nAng1 = P%orbital1%nOrbComp(iA1)
nAng2 = P%orbital2%nOrbComp(iA2)
nAng3 = P%orbital3%nOrbComp(iA3)
nCont1 = P%orbital1%nContracted(iA1)
nCont2 = P%orbital2%nContracted(iA2)
nCont3 = P%orbital3%nContracted(iA3)
start1 = P%orbital1%startOrbital(iA1)
start2 = P%orbital2%startOrbital(iA2)
start3 = P%orbital3%startOrbital(iA3)
same13 = (start1.EQ.start3)

IF(INPUT%CONT_VECTOR) THEN

IF(same13) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+&
	                                INTEGRAL%IN(int123)*v3
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+value*v3
	      OUTPUT%ResultMat(i3,i2,1)=OUTPUT%ResultMat(i3,i2,1)+&
	                                value*INPUT%VECTOR(i1)
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ELSEIF(INPUT%CONT_MATRIX) THEN

IF(same13) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+&
	                               INTEGRAL%IN(int123)*m23
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+value*m23
	      OUTPUT%ResultMat(i3,1,1)=OUTPUT%ResultMat(i3,1,1)+&
	                               value*INPUT%MATRIX(i2,i1)
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ELSE

IF(same13) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,i2,i3)=INTEGRAL%IN(int123)
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,i2,i3)=value
	      OUTPUT%ResultMat(i3,i2,i1)=value
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ENDIF

END SUBROUTINE DistributeNormal3

SUBROUTINE DistributeNormal4(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: iA1,iA2,iA3
INTEGER                  :: nAng1,nAng2,nAng3,nCont1,nCont2,nCont3
INTEGER                  :: iAng1,iAng2,iAng3,iCont1,iCont2,iCont3
INTEGER                  :: start1,start2,start3,n1,n2,n3,i1,i2,i3
LOGICAL                  :: same23
INTEGER                  :: int123
REAL(REALK)              :: v2,v3,m23,value

iA1 = P%iAng1(iAngmom)
iA2 = P%iAng2(iAngmom)
iA3 = P%iAng3(iAngmom)
nAng1 = P%orbital1%nOrbComp(iA1)
nAng2 = P%orbital2%nOrbComp(iA2)
nAng3 = P%orbital3%nOrbComp(iA3)
nCont1 = P%orbital1%nContracted(iA1)
nCont2 = P%orbital2%nContracted(iA2)
nCont3 = P%orbital3%nContracted(iA3)
start1 = P%orbital1%startOrbital(iA1)
start2 = P%orbital2%startOrbital(iA2)
start3 = P%orbital3%startOrbital(iA3)
same23 = (start2.EQ.start3)

IF(INPUT%CONT_VECTOR) THEN

IF(same23) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+&
	                                INTEGRAL%IN(int123)*v3
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
	  v3=INPUT%VECTOR(i3)
          i2=n2
          DO iAng2=1,nAng2
	    v2=INPUT%VECTOR(i2)
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+value*v3
	      OUTPUT%ResultMat(i1,i3,1)=OUTPUT%ResultMat(i1,i3,1)+value*v2
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ELSEIF(INPUT%CONT_MATRIX) THEN

IF(same23) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+&
	                               INTEGRAL%IN(int123)*m23
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    m23=INPUT%MATRIX(i2,i3)+INPUT%MATRIX(i3,i2)
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+&
	                               INTEGRAL%IN(int123)*m23
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ELSE

IF(same23) THEN
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      OUTPUT%ResultMat(i1,i2,i3)=INTEGRAL%IN(int123)
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  int123=1
  DO iCont3=1,nCont3
    n3=start3+(iCont3-1)*nAng3
    DO iCont2=1,nCont2
      n2=start2+(iCont2-1)*nAng2
      DO iCont1=1,nCont1
        n1=start1+(iCont1-1)*nAng1
        i3=n3
        DO iAng3=1,nAng3
          i2=n2
          DO iAng2=1,nAng2
	    i1=n1
	    DO iAng1=1,nAng1
	      value=INTEGRAL%IN(int123)
	      OUTPUT%ResultMat(i1,i2,i3)=value
	      OUTPUT%ResultMat(i1,i3,i2)=value
	      int123=int123+1
	      i1=i1+1
	    ENDDO
	    i2=i2+1
	  ENDDO
	  i3=i3+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

ENDIF

END SUBROUTINE DistributeNormal4

SUBROUTINE DistributeNormal5(P,INTEGRAL,INPUT,OUTPUT,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TCO_INTEGRAND)      :: P
TYPE(TCO_INTEGRALITEM)   :: INTEGRAL
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
INTEGER                  :: iAngmom,LUPRI,IPRINT
!
INTEGER                  :: iA1,iA2,iA3
INTEGER                  :: nAng1,nAng2,nAng3,nCont1,nCont2,nCont3
INTEGER                  :: iAng1,iAng2,iAng3,iCont1,iCont2,iCont3
INTEGER                  :: start1,start2,start3,n1,n2,n3,i1,i2,i3
INTEGER                  :: int123
REAL(REALK)              :: v3,m23

iA1 = P%iAng1(iAngmom)
iA2 = P%iAng2(iAngmom)
iA3 = P%iAng3(iAngmom)
nAng1 = P%orbital1%nOrbComp(iA1)
nAng2 = P%orbital2%nOrbComp(iA2)
nAng3 = P%orbital3%nOrbComp(iA3)
nCont1 = P%orbital1%nContracted(iA1)
nCont2 = P%orbital2%nContracted(iA2)
nCont3 = P%orbital3%nContracted(iA3)
start1 = P%orbital1%startOrbital(iA1)
start2 = P%orbital2%startOrbital(iA2)
start3 = P%orbital3%startOrbital(iA3)

IF(INPUT%CONT_VECTOR) THEN

int123=1
DO iCont3=1,nCont3
  n3=start3+(iCont3-1)*nAng3
  DO iCont2=1,nCont2
    n2=start2+(iCont2-1)*nAng2
    DO iCont1=1,nCont1
      n1=start1+(iCont1-1)*nAng1
      i3=n3
      DO iAng3=1,nAng3
        v3=INPUT%VECTOR(i3)
        i2=n2
        DO iAng2=1,nAng2
	  i1=n1
	  DO iAng1=1,nAng1
	    OUTPUT%ResultMat(i1,i2,1)=OUTPUT%ResultMat(i1,i2,1)+&
	                              INTEGRAL%IN(int123)*v3
	    int123=int123+1
	    i1=i1+1
	  ENDDO
	  i2=i2+1
	ENDDO
	i3=i3+1
      ENDDO
    ENDDO
  ENDDO
ENDDO

ELSEIF(INPUT%CONT_MATRIX) THEN

int123=1
DO iCont3=1,nCont3
  n3=start3+(iCont3-1)*nAng3
  DO iCont2=1,nCont2
    n2=start2+(iCont2-1)*nAng2
    DO iCont1=1,nCont1
      n1=start1+(iCont1-1)*nAng1
      i3=n3
      DO iAng3=1,nAng3
        i2=n2
        DO iAng2=1,nAng2
	  m23=INPUT%MATRIX(i2,i3)
	  i1=n1
	  DO iAng1=1,nAng1
	    OUTPUT%ResultMat(i1,1,1)=OUTPUT%ResultMat(i1,1,1)+&
	                             INTEGRAL%IN(int123)*m23
	    int123=int123+1
	    i1=i1+1
	  ENDDO
	  i2=i2+1
	ENDDO
	i3=i3+1
      ENDDO
    ENDDO
  ENDDO
ENDDO

ELSE

int123=1
DO iCont3=1,nCont3
  n3=start3+(iCont3-1)*nAng3
  DO iCont2=1,nCont2
    n2=start2+(iCont2-1)*nAng2
    DO iCont1=1,nCont1
      n1=start1+(iCont1-1)*nAng1
      i3=n3
      DO iAng3=1,nAng3
        i2=n2
        DO iAng2=1,nAng2
	  i1=n1
	  DO iAng1=1,nAng1
	    OUTPUT%ResultMat(i1,i2,i3)=INTEGRAL%IN(int123)
	    int123=int123+1
	    i1=i1+1
	  ENDDO
	  i2=i2+1
	ENDDO
	i3=i3+1
      ENDDO
    ENDDO
  ENDDO
ENDDO

ENDIF

END SUBROUTINE DistributeNormal5

SUBROUTINE swap_realPointers(a,b)
implicit none
REAL(REALK),pointer      :: a(:),b(:),c(:)

c => a
a => b
b => c

END SUBROUTINE swap_realPointers

SUBROUTINE transpose_matrix(IN,OUT,n1,n2)
implicit none
INTEGER                  :: n1,n2
REAL(REALK)              :: IN(n1,n2)
REAL(REALK)              :: OUT(n2,n1)
!
INTEGER                  :: i1,i2

DO i1=1,n1
  DO i2=1,n2
    OUT(i2,i1) = IN(i1,i2)
  ENDDO
ENDDO

END SUBROUTINE transpose_matrix

END MODULE tco_Integral_Driver
