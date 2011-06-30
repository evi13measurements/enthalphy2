MODULE integraldriver
  use TYPEDEF
  use READMOLEFILE
  use BuildBasisSet
  use READDALTONFILE
  use Matrix_module
  use Matrix_Operations
  use ODbatches
  use sphcart_matrices
  use precision
  use lstiming
!  integer,parameter :: maxPasses = 40
!****** ORBITAL
SAVE
real(realk) :: CPUTIME_Build_Integrand
real(realk) :: WALLTIME_Build_Integrand
real(realk) :: CPUTIME_Build_HermiteTUV
real(realk) :: WALLTIME_Build_HermiteTUV
real(realk) :: CPUTIME_Contract_Q
real(realk) :: WALLTIME_Contract_Q
real(realk) :: CPUTIME_Contract_P
real(realk) :: WALLTIME_Contract_P

real(realk) :: CPUTIME_DistributeHermiteQ
real(realk) :: CPUTIME_ContractEcoeffQ
real(realk) :: CPUTIME_ContractBasisQ
real(realk) :: CPUTIME_SphericalTransformQ
real(realk) :: CPUTIME_AddToTUVQ
real(realk) :: WALLTIME_DistributeHermiteQ
real(realk) :: WALLTIME_ContractEcoeffQ
real(realk) :: WALLTIME_ContractBasisQ
real(realk) :: WALLTIME_SphericalTransformQ
real(realk) :: WALLTIME_AddToTUVQ

real(realk) :: CPUTIME_DistributeHermiteP
real(realk) :: CPUTIME_ContractEcoeffP
real(realk) :: CPUTIME_ContractBasisP
real(realk) :: CPUTIME_SphericalTransformP
real(realk) :: CPUTIME_AddToPQ
real(realk) :: WALLTIME_DistributeHermiteP
real(realk) :: WALLTIME_ContractEcoeffP
real(realk) :: WALLTIME_ContractBasisP
real(realk) :: WALLTIME_SphericalTransformP
real(realk) :: WALLTIME_AddToPQ

TYPE Orbital
Character(len=80)       :: type
Logical                 :: spherical
Integer                 :: maxAngmom
Integer                 :: nAngmom
Integer                 :: nPrimitives
Integer                 :: nPasses
Integer                 :: maxContracted
Integer                 :: CCidentifier
Real(realk)             :: center(3)
!Dimensions according to nAngmom
integer,pointer     :: angmom(:) !1:nangmom
integer,pointer     :: nContracted(:)
integer,pointer     :: startOrbital(:,:)
integer,pointer     :: startprimOrbital(:,:)
integer,pointer     :: nOrbComp(:)
integer,pointer     :: nPrimOrbComp(:)
integer,pointer     :: nOrbitals(:)
!Dimensions according to nPrimitives
real(realk),pointer :: exponents(:)
!Dimensions according to nPrimitives,maxContracted,nAngmom
!real(realk),pointer :: CC(:,:,:)
type(matrixpointer),pointer :: CC(:)
TYPE(SPHMATPOINTER),pointer   :: SPH_MAT(:)
END TYPE Orbital


!****** OVERLAPPOINTER
TYPE Overlappointer
TYPE(Overlap),pointer :: p
END TYPE Overlappointer
!****** OVERLAP
TYPE Overlap
Character(len=80)   :: type
Logical             :: hermiteSingle
TYPE(Orbital)       :: orbital1,orbital2
logical             :: sameAO
Integer             :: nAngmom
Integer             :: nPrimitives
Integer             :: totOrbitals
Integer             :: nTUV
Integer             :: maxContracted
Real(realk)         :: maxGab
Real(realk)         :: maxPrimGab
!minAngmom and maxAngmom represents the actual angular momenta of the orbitals
Integer             :: minAngmom
Integer             :: maxAngmom
!startAngmom and endAngmom represents the angular momenta of the inner 
!Hermite integrals (will typically differ from min and max when dealing
!with differentiated integrals, etc.)
Integer             :: startAngmom
Integer             :: endAngmom
!ToDo: for ETUV's not attached, remember to increase dimension here, and
!      update the building of Ecoefficients
Real(realk),pointer :: distance12(:,:)
Real(realk)         :: squaredDistance
!Dimensions according to nPrimitives*nPasses
Real(realk),pointer :: center(:,:)
Real(realk),pointer :: exponents(:)
Real(realk),pointer :: reducedExponents(:)
Real(realk),pointer :: preExpFac(:)
Integer,pointer     :: iprim1(:)    !poniter to primitive of orbital 1
Integer,pointer     :: iprim2(:)    !poniter to primitive of orbital 2
!Dimensions according to nAngmom
Integer,pointer     :: angmom(:)
Integer,pointer     :: indexAng1(:)
Integer,pointer     :: indexAng2(:)
Integer,pointer     :: nOrbitals(:)
Integer,pointer     :: nOrbComp(:)
Integer,pointer     :: nContracted(:)
Integer,pointer     :: ETUVindex(:)
Integer             :: batchindex1
Integer             :: batchindex2


!McMurcie-Davidson E-coefficients (and/or F-coefficients in case of J-engine) 
Logical             :: ETUVisSet
Logical             :: sphericalETUV
Integer             :: lenETUV
Real(realk),pointer :: ETUV(:) ! ntuv, ijk, nprim, iAngmom
Real(realk),pointer :: FTUV(:,:,:)

!Pass specific variables (passes are used to increase perforance by 
!collecting OD-batches of similar structure - this decrease the time
!spent in the different contraction/matrix multiplies steps, as 
!well as reducing the overhead related to memory/cache loadings and
!subroutine calls).
Integer             :: nPasses
END TYPE Overlap

!****** INTEGRAND
TYPE Integrand
Character(len=80)    :: Operator
TYPE(OverlapPointer) :: P,Q
Integer              :: nAngmom
Integer              :: nPrimitives
Integer              :: maxContracted
Integer              :: minAngmom
Integer              :: maxAngmom
Integer              :: startAngmom
Integer              :: endAngmom
!Dimensions according to nPrimitives
Integer,pointer      :: iprimP(:)
Integer,pointer      :: iprimQ(:)
Real(realk),pointer  :: distance(:,:)
Real(realk),pointer  :: squaredDistance(:)
Real(realk),pointer  :: exponents(:)
Real(realk),pointer  :: reducedExponents(:)
Real(realk),pointer  :: integralPrefactor(:)
Real(realk)          :: ORIGO(3)
LOGICAL              :: samePQ
LOGICAL              :: kinetic
END TYPE Integrand

!****** INTEGRAL
TYPE Allocitem
Integer                :: maxPrimLHS    !used for integrals+overlap
Integer                :: maxPrimTUVLHS !used for integrals+overlap
Integer                :: maxPrimLHSpass    !used for integrals+overlap
Integer                :: maxPrimTUVLHSpass !used for integrals+overlap
Integer                :: maxnAngLHS    !used for overlap
Integer                :: maxContLHS       !used for orbital
Integer                :: maxAngmomLHS     !used for orbital
Integer                :: maxTUVLHS        !used for F-type overlap
Integer                :: maxPrimTUVijkLHS !used for overlap (Ecoeff)
Integer                :: maxTotOrbLHS     !used for temporary array (PQ)
Integer                :: maxijkLHS     !used for temporary array (PQ)

Integer                :: maxPrimRHS
Integer                :: maxPrimTUVRHS
Integer                :: maxPrimRHSpass
Integer                :: maxPrimTUVRHSpass
Integer                :: maxnAngRHS
Integer                :: maxContRHS
Integer                :: maxAngmomRHS
Integer                :: maxTUVRHS
Integer                :: maxPrimTUVijkRHS
Integer                :: maxTotOrbRHS
Integer                :: maxijkRHS
END TYPE Allocitem

TYPE Integralitem
Real(realk),pointer    :: IN(:)
Real(realk),pointer    :: OUT(:)
Real(realk),pointer    :: RTUV(:)
Real(realk),pointer    :: integrals(:)
!
Real(realk),pointer    :: Wtuv(:,:)       !nPrimitives:ntuv
Real(realk),pointer    :: tuvTUV(:,:,:) !ntuvQ(i):ntuvP:nPrimPQ
!Real(realk),pointer    :: TUVQ(:,:,:)     !ntuvP:nPrimP:nOrbQ
Real(realk),pointer    :: PQ(:,:)         !nOrbQ:nOrbP
!Generic integral intermadiates used of contraction with Ecoefficients, 
!Spherical transformation, and contraction to contracted basis
Real(realk),pointer    :: IntegralIN(:,:,:)
Real(realk),pointer    :: IntegralOUT(:,:,:)
INTEGER                :: nTUV
INTEGER                :: nEFG !multipole orders
INTEGER                :: nOrb
INTEGER                :: nAng
INTEGER                :: nPrim
LOGICAL                :: Jengine
TYPE(TUVitem),pointer  :: TUV
END TYPE Integralitem

TYPE TUVitem
INTEGER                :: nTABFJW1
INTEGER                :: nTABFJW2
Real(realk),pointer    :: TABFJW(:,:)
INTEGER,pointer        :: Tindex(:)
INTEGER,pointer        :: Uindex(:)
INTEGER,pointer        :: Vindex(:)
INTEGER,pointer        :: TUVindex(:,:,:)
TYPE(SPHMAT),pointer   :: SPH_MAT(:)
INTEGER                :: nSPHMAT
END TYPE TUVitem

TYPE LINKshell
integer,pointer :: elms(:)
INTEGER         :: DIM
END TYPE LINKshell

CONTAINS

!*************************************************************************
!*
!*                 Main integral driver
!*
!*************************************************************************
SUBROUTINE MAIN_INTEGRAL_DRIVER(LUPRI,IPRINT,INPUT,OUTPUT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
TYPE(ODmatrixitem)   :: ODmat
TYPE(ODmatrixitem),allocatable   :: ODRESmat(:)
Character(len=80)    :: IDENTIFIER
!
Real(realk)          :: TSTART,TEND
INTEGER              :: idmat
!
#if MOD_FMM
!Simen temporary to interface FMM
Integer,parameter :: LWRK = 1
Real(realk) :: WRK(LWRK)
!Simen
#endif

!
IDENTIFIER='THERMIT'

IF(IPRINT .GT. 5) WRITE(LUPRI,*)'CALLING CREATE OD'
CALL Create_ODbatches(OD_LHS,IDENTIFIER,INPUT,'LHS',LUPRI)
CALL Create_ODbatches(OD_RHS,IDENTIFIER,INPUT,'RHS',LUPRI)

!if calc Fock or exchange
IF (Input%DO_FOCK .AND. INPUT%DO_LINK) THEN
   CALL Create_ODbatchmatrices(ODmat,INPUT,'RHS',INPUT%DMAT_RHS(:,:,1),INPUT%AOdim(3),INPUT%AOdim(4),LUPRI)
   ALLOCATE(ODRESmat(Input%NDMAT_RHS))
   DO idmat = 1,Input%NDMAT_RHS
      CALL Create_ODbatchmatrices(ODRESmat(idmat),INPUT,'LHS',Output%ResultMat(:,:,1,1,idmat),INPUT%AOdim(1),INPUT%AOdim(2),LUPRI)
   ENDDO
!   DEALLOCATE(Int_Output%ResultMat)
!   Nullify(Int_Output%ResultMat)
ENDIF

CALL PRINT_OD(OD_LHS,LUPRI,IPRINT)
CALL PRINT_OD(OD_RHS,LUPRI,IPRINT)

#if MOD_FMM
  IF (Input%DO_FMM) &
    & CALL QUIT('FMM for **INTEGRALS/.LINSCA under development. Not yet working!')
! Initialize spherical multipole moments, and put to file
!IF (INPUT%DO_FMM) CALL SETMM(INPUT%DMAT_LHS,INPUT%DMAT_RHS,1,.FALSE.,.FALSE.,&
!     &                       IPRINT,f77_work,len_f77_work)
#endif

IF (INPUT%DO_JENGINE) THEN
  CALL Jengine(OD_LHS,OD_RHS,INPUT,OUTPUT,LUPRI,IPRINT)
ELSEIF(INPUT%DO_LINK) THEN
   IF(INPUT%DO_PASSES)THEN
      CALL LINKwPASS(OD_LHS,OD_RHS,ODmat,ODRESmat,INPUT,OUTPUT,LUPRI,IPRINT)
   ELSE
      CALL LINK(OD_LHS,OD_RHS,ODmat,ODRESmat,INPUT,OUTPUT,LUPRI,IPRINT)
   ENDIF
ELSE
   !ToDo Move INPUT%DO_PASSES assignment to proper places
!   INPUT%DO_PASSES = .TRUE.
  CALL McMurchieDavidson(OD_LHS,OD_RHS,INPUT,OUTPUT,LUPRI,IPRINT)
ENDIF

#if MOD_FMM
IF (INPUT%DO_FMM) CALL FMMFCK(Output%ResultMat,INPUT%DMAT_RHS,WRK,LWRK)
#endif


 IF (Input%DO_FOCK .AND. INPUT%DO_LINK ) THEN
    DO idmat = 1,Input%NDMAT_RHS
       CALL from_ODbatchmatrices_to_full(ODRESmat(idmat),INPUT,'LHS',&
            &Output%ResultMat(:,:,1,1,idmat),INPUT%AOdim(1),INPUT%AOdim(2),LUPRI)
    ENDDO
 ENDIF

 DeAllocate(OD_LHS%BATCH)
 DeAllocate(OD_RHS%BATCH)
 Nullify(OD_LHS%BATCH)
 Nullify(OD_RHS%BATCH)
 IF (Input%DO_FOCK .AND. INPUT%DO_LINK ) THEN
    CALL free_ODbatchmatrices(ODmat)
    DO idmat = 1,Input%NDMAT_RHS
       CALL free_ODbatchmatrices(ODRESmat(idmat))       
    ENDDO
 ENDIF

END SUBROUTINE MAIN_INTEGRAL_DRIVER

!*************************************************************************
!*
!*                 McMurchie-Davidson based integral driver
!*
!*************************************************************************
SUBROUTINE McMurchieDavidson(OD_LHS,OD_RHS,INPUT,OUTPUT,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
!
Integer               :: ILHS,IRHS,Start_RHS,End_RHS
TYPE(TUVitem),target  :: SharedTUV
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)       :: Alloc
TYPE(Overlap)         :: P,Q2
TYPE(Overlap),pointer :: Q(:)
TYPE(Overlap)         :: PassQ
TYPE(Integrand)       :: PQ
Integer               :: nPrimP,nPrimQ2
Integer,pointer       :: nPrimQ(:)
Integer,pointer       :: ODpassesIndex(:)
logical               :: dopasses,SET_ORBITAL
integer               :: maxpasses,nPassTypes,iPassType,numPasses,TPRINT
real(realk),allocatable :: WORK(:)
integer :: i,iend,LWORK,WORKLENGTH
Integer,pointer :: Maxpassfortype(:)

IF (IPRINT.GT.5) THEN
  CALL LSHEADER(LUPRI,'McMurchieDavidson')
  WRITE(LUPRI,'(1X,A)') ''
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') '***                     McMurchieDavidson                      ***'
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') ''
ENDIF
CALL initTUVitem(sharedTUV,Integral,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)

doPasses  = INPUT%DO_PASSES
maxPasses = INPUT%maxPasses

IF(doPasses) WRITE(LUPRI,*)'Integrals are calculated using passes with maxpasses=',maxPasses

IF(INPUT%setOVERLAPoutside)THEN !THIS IS THE MOST EFFICIENT BUT NOT LINEAR IN MEMMORY
   CALL initAlloc(Alloc,LUPRI,IPRINT,'RHS')
   NULLIFY(Q)
   ALLOCATE(Q(OD_RHS%nbatches))
   NULLIFY(nPrimQ)
   ALLOCATE(nPrimQ(OD_RHS%nbatches))
   DO IRHS=1,OD_RHS%nbatches
      CALL SET_Overlap(Q(IRHS),nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS')
   ENDDO

   NULLIFY(ODpassesIndex)
   ALLOCATE(ODpassesIndex(OD_RHS%nbatches))
   CALL SelectODPassTypes(ODpassesIndex,Q,OD_RHS%nbatches,nPassTypes,IPRINT,LUPRI)

   CALL initAlloc(Alloc,LUPRI,IPRINT,'LHS')
   CALL SET_ALLOC(Alloc,Input,OD_LHS,'LHS',IPRINT,LUPRI)
   
   !$OMP PARALLEL PRIVATE(integral,P,PQ,ILHS,IRHS,iPassType,numpasses,Start_RHS,End_RHS,nPrimP,PassQ,WORK,WORKLENGTH,LWORK)    
   integral%TUV => sharedTUV
   CALL INIT_OVERLAP(P,Alloc,Input,OD_LHS,'LHS',1,IPRINT,LUPRI)
   CALL allocateIntegrals(PQ,Integral,Alloc,maxPasses)
   CALL INIT_PASS(PassQ,Alloc,Input,OD_RHS,'RHS',2,maxPasses,IPRINT,LUPRI)
   WORKLENGTH = MAX(Alloc%maxprimRHS*Alloc%maxContRHS*Alloc%maxContRHS,&
        & Alloc%maxprimLHS*Alloc%maxContLHS*Alloc%maxContLHS)
   WORKLENGTH = MAX(WORKLENGTH,Alloc%maxijkLHS*Alloc%maxijkLHS,Alloc%maxijkRHS*Alloc%maxijkRHS)
   ALLOCATE(WORK(WORKLENGTH))
   LWORK = 1
   !$OMP DO SCHEDULE(DYNAMIC,1)
   DO ILHS = 1,OD_LHS%nbatches
      CALL SET_INITIALIZED_Overlap(P,nPrimP,Input,SharedTUV,Integral,Alloc,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS',.FALSE.,.TRUE.)
      IF (nPrimP .GT. 0) THEN
         CALL Determine_RHS_loop(INPUT,OD_RHS,ILHS,Start_RHS,End_RHS)
         DO iPassType=1,nPassTypes
            numPasses = 0
            DO IRHS = Start_RHS,End_RHS
               IF (ODpassesIndex(IRHS).EQ.iPassType) THEN
                  IF(OD_LHS%BATCH(ILHS)%maxGAB*OD_RHS%BATCH(IRHS)%maxGAB &
                       & .GT. INPUT%CS_THRESHOLD  .OR. .NOT. INPUT%CS_SCREEN)THEN
                     IF (IPRINT.GT.3) WRITE(LUPRI,'(1X,A,I3,A,I3)') 'Overlap ditributions P',ILHS,' and Q',IRHS
                     IF(nPrimQ(IRHS) .GT. 0)THEN
                        IF (doPasses .AND..NOT. ((ILHS.EQ.IRHS).AND.INPUT%sameODs)) THEN
                           numPasses = numPasses + 1
                           CALL AddOverlapToPass(PassQ,Q(IRHS),numPasses,maxPasses,LUPRI,IPRINT)
                           IF (numPasses.EQ.maxPasses) THEN
                              CALL ExplicitIntegrals(Integral,PQ,P,PassQ,INPUT,OUTPUT,ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                              CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
                              numPasses = 0
                           ENDIF
                        ELSE
                           CALL ExplicitIntegrals(Integral,PQ,P,Q(IRHS),INPUT,OUTPUT,ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                           CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
                        ENDIF
                     ENDIF! PS-screening
                  ENDIF! CS-screening
               ENDIF ! Correct passtype
            ENDDO !RHS
            IF (numPasses.GT.0) THEN
               CALL ExplicitIntegrals(Integral,PQ,P,PassQ,INPUT,OUTPUT,ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
               CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
            ENDIF
         ENDDO !Passtypes
      ENDIF
   ENDDO
   !$OMP END DO 
   
   CALL deallocateIntegrals(PQ,Integral)
   CALL FreePass(PassQ)
   CALL FREE_OVERLAP(P)
   DEALLOCATE(WORK)
   
   !$OMP END PARALLEL
   
   DO IRHS=1,OD_RHS%nbatches
      IF (nPrimQ(IRHS) .GT. 0) CALL FREE_OVERLAP(Q(IRHS))
   ENDDO
   
   DEALLOCATE(nPrimQ)
   DEALLOCATE(Q)
   DEALLOCATE(ODpassesIndex)
ELSE  !THIS IS LINEAR IN MEMMORY BUT YOU ARE SETTING EACH OVERLAP OD_LHS%nbatches times
   NULLIFY(ODpassesIndex)
   ALLOCATE(ODpassesIndex(OD_RHS%nbatches))
   NULLIFY(nPrimQ)
   ALLOCATE(nPrimQ(OD_RHS%nbatches)) 
   CALL SelectODPassTypesFromODbatch(ODpassesIndex,nPrimQ,OD_RHS,OD_RHS%nbatches,nPassTypes,&
        &maxpasses,maxpassfortype,INPUT,IPRINT,LUPRI,'RHS')

   CALL initAlloc(Alloc,LUPRI,IPRINT,'Both')
   CALL SET_ALLOC2(Alloc,Input,OD_RHS,'RHS',IPRINT,LUPRI,OD_RHS%nbatches,ODpassesIndex,nPassTypes,maxPassfortype,maxpasses)
   CALL SET_ALLOC(Alloc,Input,OD_LHS,'LHS',IPRINT,LUPRI)
   DEALLOCATE(maxPassfortype)
   NULLIFY(maxPassfortype)

   !$OMP PARALLEL PRIVATE(integral,P,Q2,PQ,ILHS,IRHS,iPassType,numpasses,Start_RHS,End_RHS,nPrimP,PassQ,WORK,WORKLENGTH,LWORK) 
   
   integral%TUV => sharedTUV
   CALL INIT_OVERLAP(P,Alloc,Input,OD_LHS,'LHS',1,IPRINT,LUPRI)
   CALL INIT_OVERLAP(Q2,Alloc,Input,OD_RHS,'RHS',2,IPRINT,LUPRI)
   CALL allocateIntegrals(PQ,Integral,Alloc,maxPasses)
   CALL INIT_PASS(PassQ,Alloc,Input,OD_RHS,'RHS',2,maxPasses,IPRINT,LUPRI)
   WORKLENGTH = MAX(Alloc%maxprimRHS*Alloc%maxContRHS*Alloc%maxContRHS,&
        & Alloc%maxprimLHS*Alloc%maxContLHS*Alloc%maxContLHS)
   WORKLENGTH = MAX(WORKLENGTH,Alloc%maxijkLHS*Alloc%maxijkLHS,Alloc%maxijkRHS*Alloc%maxijkRHS)
   ALLOCATE(WORK(WORKLENGTH))
   LWORK = 1
   !$OMP DO SCHEDULE(DYNAMIC,1)
   DO ILHS = 1,OD_LHS%nbatches
      ! CALL initAlloc(Alloc,LUPRI,IPRINT,'LHS')
      ! CALL SET_Overlap(P,nPrimP,Input,SharedTUV,Integral,Alloc,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS')
      CALL SET_INITIALIZED_Overlap(P,nPrimP,Input,SharedTUV,Integral,Alloc,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS',.FALSE.,.TRUE.)
      IF (nPrimP .GT. 0) THEN
         CALL Determine_RHS_loop(INPUT,OD_RHS,ILHS,Start_RHS,End_RHS)
         DO iPassType=1,nPassTypes
            numPasses = 0
            SET_ORBITAL = .TRUE.
            DO IRHS = Start_RHS,End_RHS
               IF (ODpassesIndex(IRHS).EQ.iPassType) THEN
                  IF(OD_LHS%BATCH(ILHS)%maxGAB*OD_RHS%BATCH(IRHS)%maxGAB &
                       & .GT. INPUT%CS_THRESHOLD  .OR. .NOT. INPUT%CS_SCREEN)THEN
                     IF (IPRINT.GT.3) WRITE(LUPRI,'(1X,A,I3,A,I3)') 'Overlap ditributions P',ILHS,' and Q',IRHS
                     IF(nPrimQ(IRHS) .GT. 0)THEN
                        IF (doPasses .AND..NOT. ((ILHS.EQ.IRHS).AND.INPUT%sameODs)) THEN
                           CALL SET_INITIALIZED_Overlap(Q2,nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,OD_RHS%BATCH(IRHS),2,&
                                & LUPRI,IPRINT,'RHS',.TRUE.,SET_ORBITAL)
                           SET_ORBITAL=.FALSE.!all overlaps of this type have the same orbitals -> set once
                           numPasses = numPasses + 1
                           CALL AddOverlapToPass(PassQ,Q2,numPasses,maxPasses,LUPRI,IPRINT)
                           IF (numPasses.EQ.maxPasses) THEN
                              CALL ExplicitIntegrals(Integral,PQ,P,PassQ,INPUT,OUTPUT,ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                              CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)

                              numPasses = 0
                           ENDIF
                        ELSE
                           CALL SET_INITIALIZED_Overlap(Q2,nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,OD_RHS%BATCH(IRHS),2,&
                                & LUPRI,IPRINT,'RHS',.TRUE.,.TRUE.)
                           CALL ExplicitIntegrals(Integral,PQ,P,Q2,INPUT,OUTPUT,ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                           CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)

                        ENDIF
                     ENDIF! PS-screening
                  ENDIF! CS-screening
               ENDIF ! Correct passtype
            ENDDO !RHS
            IF (numPasses.GT.0) THEN
               CALL ExplicitIntegrals(Integral,PQ,P,PassQ,INPUT,OUTPUT,ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
               CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
            ENDIF
         ENDDO !Passtypes
      ENDIF
   ENDDO
   !$OMP END DO 
      CALL deallocateIntegrals(PQ,Integral)
   CALL FREE_OVERLAP(P)
   CALL FREE_OVERLAP(Q2)
   CALL FreePass(PassQ)
   DEALLOCATE(WORK)
   
   !$OMP END PARALLEL
   
   DEALLOCATE(ODpassesIndex)
   DEALLOCATE(nPrimQ)
ENDIF

   CALL freeTUVitem(sharedTUV,Input)

END SUBROUTINE McMurchieDavidson


SUBROUTINE ExplicitIntegrals(Integral,PQ,P,Q,INPUT,OUTPUT,ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
implicit none
Integer              :: ILHS,IRHS,LUPRI,IPRINT,LWORK,WORKLENGTH
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
TYPE(Integralitem)   :: Integral
TYPE(Overlap)        :: P,Q
TYPE(Integrand)      :: PQ
REAL(REALK)          :: WORK(WORKLENGTH)
REAL(REALK)          :: CPU1,CPU2,WALL1,WALL2

CALL GETTIM(CPU1,WALL1)
CALL Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,LUPRI,IPRINT)
CALL GETTIM(CPU2,WALL2)
CPUTIME_Build_Integrand = CPUTIME_Build_Integrand+CPU2-CPU1
WALLTIME_Build_Integrand = WALLTIME_Build_Integrand+ WALL2-WALL1
!   Hermite 2-electron integral over general operator w
CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
CALL GETTIM(CPU1,WALL1)
CPUTIME_Build_HermiteTUV = CPUTIME_Build_HermiteTUV+CPU1-CPU2
WALLTIME_Build_HermiteTUV = WALLTIME_Build_HermiteTUV + WALL1-WALL2
CALL Contract_Q(INTEGRAL,PQ,Input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
CALL GETTIM(CPU2,WALL2)
CPUTIME_Contract_Q = CPUTIME_Contract_Q + CPU2-CPU1
WALLTIME_Contract_Q = WALLTIME_Contract_Q + WALL2-WALL1
CALL Contract_P(INTEGRAL,PQ,input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
CALL GETTIM(CPU1,WALL1)
CPUTIME_Contract_P = CPUTIME_Contract_P + CPU1-CPU2
WALLTIME_Contract_P = WALLTIME_Contract_P + WALL1-WALL2

END SUBROUTINE ExplicitIntegrals
!*************************************************************************
!*
!*                 End of McMurchie-Davidson
!*
!*************************************************************************

!*************************************************************************
!*
!*                 Jengine based integral driver
!*
!*************************************************************************
SUBROUTINE Jengine(OD_LHS,OD_RHS,INPUT,OUTPUT,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
!
Integer               :: ILHS,IRHS,Start_LHS,End_LHS,NFTUVbatches
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)       :: Alloc
TYPE(TUVitem),target  :: SharedTUV
TYPE(Overlap)         :: P
TYPE(Overlap),pointer :: Q(:)
TYPE(Overlap),pointer :: F(:)
TYPE(Integrand)       :: PQ
Integer               :: nPrimP, nPrimQ,LWORK,WORKLENGTH,maxpasses
real(realk),allocatable :: WORK(:)

IF (IPRINT.GT.5) THEN
  CALL LSHEADER(LUPRI,'Jengine')
  WRITE(LUPRI,'(1X,A)') ''
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') '***                          Jengine                           ***'
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') ''
ENDIF

maxpasses = 1
CALL initTUVitem(sharedTUV,Integral,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
CALL initAlloc(Alloc,LUPRI,IPRINT,'RHS')

NULLIFY(Q)
ALLOCATE(Q(OD_RHS%nbatches))
CALL SET_FTUVbatches(F,NFTUVbatches,OD_RHS,Q,Input,SharedTUV,Integral,Alloc,LUPRI,IPRINT)
DEALLOCATE(Q)
Q=>F
NULLIFY(F)

CALL initAlloc(Alloc,LUPRI,IPRINT,'LHS')
CALL SET_ALLOC(Alloc,Input,OD_LHS,'LHS',IPRINT,LUPRI)

!$OMP PARALLEL PRIVATE(integral,P,PQ,ILHS,IRHS,Start_LHS,End_LHS,nPrimP,nPrimQ,WORK,LWORK,WORKLENGTH) 
integral%TUV => sharedTUV
CALL INIT_OVERLAP(P,Alloc,Input,OD_LHS,'LHS',1,IPRINT,LUPRI)
WORKLENGTH = MAX(Alloc%maxprimRHS*Alloc%maxContRHS*Alloc%maxContRHS,&
     & Alloc%maxprimLHS*Alloc%maxContLHS*Alloc%maxContLHS)
WORKLENGTH = MAX(WORKLENGTH,Alloc%maxijkLHS*Alloc%maxijkLHS,Alloc%maxijkRHS*Alloc%maxijkRHS)
ALLOCATE(WORK(WORKLENGTH))
LWORK = 1
CALL allocateIntegrals(PQ,Integral,Alloc,1)

!$OMP DO SCHEDULE(DYNAMIC,1)
DO ILHS = 1,OD_LHS%nbatches
   CALL SET_INITIALIZED_Overlap(P,nPrimP,Input,SharedTUV,Integral,Alloc,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS',.FALSE.,.TRUE.)
 IF (nPrimP.GT.0) THEN
!  NULLIFY(Integral%TUVQ)
!  ALLOCATE(Integral%TUVQ(P%nPrimitives,P%nTUV,1))
  CALL DZERO(Integral%integrals,P%nPrimitives*P%nTUV*1)
  DO IRHS = 1,NFTUVbatches
!Simen  The Overlap object should inherit maxGab, and should be passed on to FTUV-batches
!     IF(OD_LHS%BATCH(ILHS)%maxGAB*OD_RHS%BATCH(IRHS)%maxGAB &
!          & .GT. INPUT%CS_THRESHOLD  .OR. .NOT. INPUT%CS_SCREEN)THEN
     IF (.TRUE.) THEN
        IF (IPRINT.GT.3) WRITE(LUPRI,'(1X,A,I3,A,I3)') 'Overlap ditributions P',ILHS,' and Q',IRHS
        CALL Build_Integrand(PQ,P,Q(IRHS),INPUT,ILHS,IRHS,LUPRI,IPRINT)
        !   Hermite 2-electron integral over general operator w 
        CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
        CALL ContractFTUV(INTEGRAL,P,Q(IRHS),LUPRI,IPRINT)
        !    CALL CLEAR_integral(INTEGRAL)
     ENDIF
  ENDDO
  CALL Contract_P(INTEGRAL,PQ,input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
  CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
!  DEALLOCATE(Integral%TUVQ)
 ENDIF
ENDDO
!$OMP END DO

CALL deallocateIntegrals(PQ,Integral)
CALL FREE_OVERLAP(P)
DEALLOCATE(WORK)

!$OMP END PARALLEL

DO IRHS = 1,NFTUVbatches
  IF (Q(IRHS)%nPrimitives.GT.0) CALL FREE_FTUVOVERLAP(Q(IRHS))
ENDDO
DEALLOCATE(Q)
NULLIFY(Q)

CALL freeTUVitem(sharedTUV,Input)

END SUBROUTINE Jengine
!*************************************************************************
!*
!*                 End of Jengine
!*
!*************************************************************************

!*************************************************************************
!*
!*                 Link based integral driver    JCP 109,1663 
!*
!*************************************************************************
SUBROUTINE LINK(OD_LHS,OD_RHS,ODmat,ODRESmat,INPUT,OUTPUT,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
TYPE(ODmatrixitem)   :: ODmat,ODRESmat(:)
!
Integer               :: ILHS,IRHS,Start_RHS,End_RHS,i,j
TYPE(TUVitem),target  :: SharedTUV
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)       :: Alloc
TYPE(Overlap),pointer :: P(:)
TYPE(Overlap),pointer :: Q(:)
TYPE(Integrand)       :: PQ
Integer,pointer       :: nPrimP(:)
Integer,pointer       :: nPrimQ(:)

Real(realk),pointer   :: RED_GAB_LHS(:,:)
Real(realk),pointer   :: RED_GAB_RHS(:,:)
Real(realk),pointer   :: RED_DMAT_LHS(:,:,:)
Real(realk),pointer   :: RED_DMAT_RHS(:,:,:)
Integer               :: dim1,dim2,dim3,dim4,A,B,C,D,idmat,endC
Integer               :: LISTSIZE!,LISTSIZE1,LISTSIZE2
Integer               :: nB,nC,nD,OLDI,RHSINDEX,startB,startD,I2,JK
Integer,allocatable   :: LIST(:)!,LIST2(:),LIST1(:)
LOGICAL,allocatable   :: DoINT(:)
real(realk),allocatable :: maxLHSGAB(:),maxRHSGAB(:)
Integer                 :: LWORK,WORKLENGTH,STARTC
real(realk),allocatable :: WORK(:)
!real(realk),allocatable :: SORTING(:)
!integer,allocatable     :: bshell(:)
real(realk)           :: maxLHSELM,maxRHSELM,DMATELM1,DMATELM2,MAXDMAT,CPUTIME1,CPUTIME2
!LOGICAL               :: A_LOOP_DONE,B_LOOP_DONE,C_LOOP_DONE,D_LOOP_DONE,UNIQUE
TYPE(LINKshell),pointer  :: brashell(:),ketshell(:),ML(:,:)
real(realk)           :: CPUTimeLINK,WALLTIME1,WALLTIME2,WALLTIMELINK
real(realk)           :: CPUTimeEXCHANGE,WALLTimeEXCHANGE
real(realk)           :: CPUTimeEXPLICIT,WALLTimeEXPLICIT
real(realk)           :: CPUTIMESTART,WALLTIMESTART,CPUTIMEEND,WALLTIMEEND
real(realk) :: CPUTimeMERGE,WALLTimeMERGE
logical        :: NOELEMENTSADDED

CALL GETTIM(CPUTIMESTART,WALLTIMESTART)

CPUTIMELINK = 0.D0
WALLTIMELINK = 0.D0
CPUTIMEEXCHANGE = 0.D0
WALLTIMEEXCHANGE = 0.D0
CPUTIMEEXPLICIT = 0.D0
WALLTIMEEXPLICIT = 0.D0
CPUTIME_Build_Integrand = 0.D0
WALLTIME_Build_Integrand = 0.D0
CPUTIME_Build_HermiteTUV = 0.D0
WALLTIME_Build_HermiteTUV = 0.D0
CPUTIME_Contract_Q = 0.D0
WALLTIME_Contract_Q = 0.D0
CPUTIME_Contract_P = 0.D0
WALLTIME_Contract_P = 0.D0
CPUTimeMERGE = 0.D0
WALLTimeMERGE = 0.D0

CPUTIME_DistributeHermiteQ = 0.D0
CPUTIME_ContractEcoeffQ = 0.D0
CPUTIME_ContractBasisQ = 0.D0
CPUTIME_SphericalTransformQ = 0.D0
CPUTIME_AddToTUVQ = 0.D0
WALLTIME_DistributeHermiteQ = 0.D0
WALLTIME_ContractEcoeffQ = 0.D0
WALLTIME_ContractBasisQ = 0.D0
WALLTIME_SphericalTransformQ = 0.D0
WALLTIME_AddToTUVQ = 0.D0

CPUTIME_DistributeHermiteP = 0.D0
CPUTIME_ContractEcoeffP = 0.D0
CPUTIME_ContractBasisP = 0.D0
CPUTIME_SphericalTransformP = 0.D0
CPUTIME_AddToPQ = 0.D0
WALLTIME_DistributeHermiteP = 0.D0
WALLTIME_ContractEcoeffP = 0.D0
WALLTIME_ContractBasisP = 0.D0
WALLTIME_SphericalTransformP = 0.D0
WALLTIME_AddToPQ = 0.D0

IF(.NOT. INPUT%CS_SCREEN)THEN
   WRITE(LUPRI,*)'Link requires screening, but you have turned off screening'
   CALL QUIT('Link requires screening, but you have turned off screening')
ENDIF

IF (IPRINT.GT.5) THEN
  IF(INPUT%DO_DALINK)THEN
     CALL LSHEADER(LUPRI,'DaLinK')
     WRITE(LUPRI,'(1X,A)') ''
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') '***                         DaLINK                             ***'
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') ''
  ELSE
     CALL LSHEADER(LUPRI,'LinK')
     WRITE(LUPRI,'(1X,A)') ''
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') '***                           LINK                             ***'
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') ''
  ENDIF
ENDIF

!REDUCE THE GAB MATRIX TO AOBATCH INDEXES

dim1=INPUT%AO(1)%p%nbatches
dim2=INPUT%AO(2)%p%nbatches
dim3=INPUT%AO(3)%p%nbatches
dim4=INPUT%AO(4)%p%nbatches
NULLIFY(RED_GAB_LHS)
NULLIFY(RED_GAB_RHS)
ALLOCATE(RED_GAB_LHS(dim1,dim2))
ALLOCATE(RED_GAB_RHS(dim3,dim4))
CALL DZERO(RED_GAB_LHS,dim1*dim2)
CALL DZERO(RED_GAB_RHS,dim3*dim4)
CALL REDUCE_MAT(LUPRI,INPUT%GAB_LHS,INPUT%AOdim(1),INPUT%AOdim(2),&
     &RED_GAB_LHS,dim1,dim2,INPUT,'LHS',.TRUE.)
CALL REDUCE_MAT(LUPRI,INPUT%GAB_RHS,INPUT%AOdim(3),INPUT%AOdim(4),&
     &RED_GAB_RHS,dim3,dim4,INPUT,'RHS',.TRUE.)

!CALCULATE MAXGAB VECTOR AND MAX ELEMENT IN GAB
ALLOCATE(maxLHSGAB(dim1))
CALL MAXGABELM(RED_GAB_LHS,dim1,dim2,maxLHSGAB,maxLHSELM)
ALLOCATE(maxRHSGAB(dim3))
CALL MAXGABELM(RED_GAB_RHS,dim3,dim4,maxRHSGAB,maxRHSELM)

! determine significant bra shell pairs

NULLIFY(brashell)
ALLOCATE(brashell(dim1))
CALL DETERMINE_SHELL_PAIRS(dim1,dim2,RED_GAB_LHS,INPUT%CS_THRESHOLD,&
        &maxRHSELM,brashell,INPUT%sameLHSaos)

NULLIFY(ketshell)
ALLOCATE(ketshell(dim3))
CALL DETERMINE_SHELL_PAIRS(dim3,dim4,RED_GAB_RHS,INPUT%CS_THRESHOLD,&
        &maxLHSELM,ketshell,.FALSE.)

NULLIFY(RED_DMAT_RHS)
ALLOCATE(RED_DMAT_RHS(dim3,dim4,Input%NDMAT_RHS))
CALL DZERO(RED_DMAT_RHS,dim3*dim4*Input%NDMAT_RHS)
DO idmat=1,Input%NDMAT_RHS
   CALL REDUCE_MAT(LUPRI,INPUT%DMAT_RHS(:,:,idmat),INPUT%AOdim(3),&
        &INPUT%AOdim(4),RED_DMAT_RHS(:,:,idmat),dim3,dim4,INPUT,'RHS',.FALSE.)
ENDDO

IF(INPUT%DO_DALINK)THEN
   NULLIFY(RED_DMAT_LHS)
   ALLOCATE(RED_DMAT_LHS(dim1,dim2,Input%NDMAT_LHS))
   CALL DZERO(RED_DMAT_LHS,dim1*dim2*Input%NDMAT_LHS)
   DO idmat=1,Input%NDMAT_LHS
      CALL REDUCE_MAT(LUPRI,INPUT%DMAT_LHS(:,:,idmat),INPUT%AOdim(3),&
           &INPUT%AOdim(4),RED_DMAT_LHS(:,:,idmat),dim1,dim2,INPUT,'LHS',.FALSE.)
   ENDDO
ENDIF

NULLIFY(ML)
ALLOCATE(ML(dim1,Input%NDMAT_RHS))
CALL DETERMINE_BRAKET_PAIRS(dim1,dim2,dim3,Input%NDMAT_RHS,RED_DMAT_RHS,&
     &maxLHSGAB,maxRHSGAB,INPUT%CS_THRESHOLD,ML)

CALL initTUVitem(sharedTUV,Integral,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
CALL initAlloc(Alloc,LUPRI,IPRINT,'Both')

NULLIFY(Q)
ALLOCATE(Q(OD_RHS%nbatches))
NULLIFY(nPrimQ)
ALLOCATE(nPrimQ(OD_RHS%nbatches))
DO IRHS=1,OD_RHS%nbatches
  CALL SET_Overlap(Q(IRHS),nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS')
ENDDO

NULLIFY(P)
ALLOCATE(P(OD_LHS%nbatches))
NULLIFY(nPrimP)
ALLOCATE(nPrimP(OD_LHS%nbatches))
DO ILHS=1,OD_LHS%nbatches
  CALL SET_Overlap(P(ILHS),nPrimP(ILHS),Input,SharedTUV,Integral,Alloc,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS')
ENDDO

IF(INPUT%DO_DALINK)THEN
   IF(INPUT%sameLHSaos) THEN
      DO C=1,dim3
         DO nD=1,brashell(C)%DIM
            D=ketshell(C)%elms(nD)
            IRHS=(C-1)*dim1-(C*(C-1)/2)+D 
            Q(IRHS)%batchindex1=C
            Q(IRHS)%batchindex2=D
         ENDDO
      ENDDO
   ELSE
      DO C=1,dim3
         DO nD=1,brashell(C)%DIM
            D=ketshell(C)%elms(nD)
            IRHS=(C-1)*dim1+D
            Q(IRHS)%batchindex1=C
            Q(IRHS)%batchindex2=D
         ENDDO
      ENDDO
   ENDIF
ENDIF

idmat=1!,Input%NDMAT_RHS !this needs to be addressed

!$OMP PARALLEL PRIVATE(integral,PQ,IRHS,RHSINDEX,A,nB,B,ILHS,I,nC,C,nD,D,&
!$OMP LIST,DoINT,LISTSIZE,WORK,LWORK,WORKLENGTH,STARTC,DMATELM1,DMATELM2,MAXDMAT,&
!$OMP NOELEMENTSADDED)

CALL allocateIntegrals(PQ,Integral,Alloc,1)
ALLOCATE(LIST(dim3*dim4))
ALLOCATE(DoINT(dim3*dim4))
WORKLENGTH = MAX(Alloc%maxprimRHS*Alloc%maxContRHS*Alloc%maxContRHS,&
     & Alloc%maxprimLHS*Alloc%maxContLHS*Alloc%maxContLHS)
WORKLENGTH = MAX(WORKLENGTH,Alloc%maxijkLHS*Alloc%maxijkLHS,Alloc%maxijkRHS*Alloc%maxijkRHS)
ALLOCATE(WORK(WORKLENGTH))
LWORK = 1
!$OMP DO SCHEDULE(DYNAMIC,1)
DO A=1,dim1
   integral%TUV => sharedTUV
   DO nB=1,brashell(A)%DIM
      B=brashell(A)%elms(nB)
      ILHS=(A-1)*dim1+B
      IF(INPUT%sameLHSaos) ILHS=(A-1)*dim1-(A*(A-1)/2)+B 
      IF(nPrimP(ILHS) .GT. 0)THEN
         !CREATING DoINT LIST from the AC element of the density matrix
         DO IRHS=1,ILHS
            DoINT(IRHS) = .FALSE.
         ENDDO
         DO nC=1,ML(A,idmat)%DIM
            C = ML(A,idmat)%elms(nC) 
            NOELEMENTSADDED = .TRUE.
            STARTC = (C-1)*dim1-(C*(C-1)/2)
            IF(INPUT%sameRHSaos)THEN
               DO nD=1,ketshell(C)%DIM
                  D=ketshell(C)%elms(nD)
                  IF(RED_DMAT_RHS(A,C,idmat)*RED_GAB_RHS(A,B)&
                       &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD )EXIT
                  NOELEMENTSADDED = .FALSE.
                  IF(D .GE. C)THEN
                     IRHS = STARTC+D
                     DoINT(IRHS) = .TRUE.
                  ELSE
                     IRHS=(D-1)*dim1-(D*(D-1)/2)+C
                     DoINT(IRHS) = .TRUE.
                  ENDIF
               ENDDO
            ELSE
               DO nD=1,ketshell(C)%DIM
                  D=ketshell(C)%elms(nD)
                  IF(RED_DMAT_RHS(A,C,idmat)*RED_GAB_RHS(A,B)&
                       &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD )EXIT
                  NOELEMENTSADDED = .FALSE.
                  DoINT((C-1)*dim1+D)=.TRUE.
               ENDDO
            ENDIF
            IF(NOELEMENTSADDED)EXIT !NO ELEMENTS ADDED IN D LOOP
         ENDDO
         !Add to DoINT LIST from the BC element of the density matrix
         DO nC=1,ML(A,idmat)%DIM
            C = ML(A,idmat)%elms(nC)
            STARTC=(C-1)*dim1-(C*(C-1)/2)
            NOELEMENTSADDED = .TRUE.
            IF(INPUT%sameRHSaos)THEN
               DO nD=1,ketshell(C)%DIM
                  D=ketshell(C)%elms(nD)
                  IF(RED_DMAT_RHS(B,C,idmat)*RED_GAB_RHS(A,B)&
                       &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD ) EXIT
                  NOELEMENTSADDED = .FALSE.
                  IF(D .GE. C)THEN
                     IRHS = STARTC+D
                     DoINT(IRHS)=.TRUE.
                  ELSE
                     IRHS = (D-1)*dim1-(D*(D-1)/2)+C
                     DoINT(IRHS)=.TRUE.
                  ENDIF
               ENDDO
            ELSE
               DO nD=1,ketshell(C)%DIM
                  D=ketshell(C)%elms(nD)
                  IF(RED_DMAT_RHS(B,C,idmat)*RED_GAB_RHS(A,B)&
                       &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD ) EXIT
                  NOELEMENTSADDED = .FALSE.
                  Doint((C-1)*dim1+D)=.TRUE.
               ENDDO
            ENDIF
            IF(NOELEMENTSADDED)EXIT !NO ELEMENTS ADDED IN D LOOP
         ENDDO
         CALL GETTIM(CPUTIME1,WALLTIME1)
         CALL BUILD_LIST(LUPRI,LIST,DoINT,LISTSIZE,dim3*dim4,ILHS) !assume sameOD
         CALL GETTIM(CPUTIME2,WALLTIME2)
         CPUTIMEMERGE = CPUTIMEMERGE + CPUTIME2-CPUTIME1
         WALLTIMEMERGE = WALLTIMEMERGE + WALLTIME2-WALLTIME1
!         CALL allocateIntegrals(PQ,Integral,Alloc,1)
         DO RHSINDEX=1,LISTSIZE
            IRHS=LIST(RHSINDEX)
            IF(nPrimQ(IRHS) .GT. 0)THEN
               IF(INPUT%DO_DALINK )THEN
                  !EXTRA SCREENING
                  C=Q(IRHS)%batchindex2
                  D=Q(IRHS)%batchindex2
                  DMATELM1=RED_DMAT_RHS(A,C,idmat)*RED_DMAT_LHS(B,D,idmat)
                  DMATELM2=RED_DMAT_RHS(B,C,idmat)*RED_DMAT_LHS(A,D,idmat)
                  MAXDMAT=MAX(DMATELM1,DMATELM2)
                  IF(MAXDMAT*RED_GAB_LHS(A,B)*RED_GAB_RHS(C,D) .GE. INPUT%CS_THRESHOLD )THEN 
                     CALL GETTIM(CPUTIME1,WALLTIME1)
                     CALL ExplicitIntegrals(Integral,PQ,P(ILHS),Q(IRHS),INPUT,OUTPUT,&
                          & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                     CALL GETTIM(CPUTIME2,WALLTIME2)
                     CPUTIMEEXPLICIT = CPUTIMEEXPLICIT + CPUTIME2-CPUTIME1
                     WALLTIMEEXPLICIT = WALLTIMEEXPLICIT + WALLTIME2-WALLTIME1
                     !$OMP CRITICAL
                     IF(INPUT%DRHS_SYM)THEN
                        CALL Distribute_Exchange3SYM(INTEGRAL,PQ,Integral%PQ,&
                             & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                             & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                     ELSE
                        CALL Distribute_Exchange3(INTEGRAL,PQ,Integral%PQ,&
                             & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                             & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                     ENDIF
                     !$OMP END CRITICAL                     
                     CALL GETTIM(CPUTIME1,WALLTIME1)
                     CPUTIMEEXCHANGE = CPUTIMEEXCHANGE + CPUTIME1-CPUTIME2
                     WALLTIMEEXCHANGE = WALLTIMEEXCHANGE + WALLTIME1-WALLTIME2
                  ENDIF
               ELSE
                  CALL GETTIM(CPUTIME1,WALLTIME1)
                  CALL ExplicitIntegrals(Integral,PQ,P(ILHS),Q(IRHS),INPUT,OUTPUT,&
                       & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                  CALL GETTIM(CPUTIME2,WALLTIME2)
                  CPUTIMEEXPLICIT = CPUTIMEEXPLICIT + CPUTIME2-CPUTIME1
                  WALLTIMEEXPLICIT = WALLTIMEEXPLICIT + WALLTIME2-WALLTIME1
                  !$OMP CRITICAL
                  IF(INPUT%DRHS_SYM)THEN
                     CALL Distribute_Exchange3SYM(INTEGRAL,PQ,Integral%PQ,&
                          & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                          & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT)
                  ELSE
                     CALL Distribute_Exchange3(INTEGRAL,PQ,Integral%PQ,&
                          & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                          & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT)
                  ENDIF
                  !$OMP END CRITICAL
                  CALL GETTIM(CPUTIME1,WALLTIME1)
                  CPUTIMEEXCHANGE = CPUTIMEEXCHANGE + CPUTIME1-CPUTIME2
                  WALLTIMEEXCHANGE = WALLTIMEEXCHANGE + WALLTIME1-WALLTIME2
               ENDIF
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDDO
!$OMP END DO 

DEALLOCATE(LIST)
DEALLOCATE(DoINT)
CALL deallocateIntegrals(PQ,Integral)
DEALLOCATE(WORK)

!$OMP END PARALLEL

DO IRHS=1,OD_RHS%nbatches
  IF (nPrimQ(IRHS) .GT. 0) CALL FREE_OVERLAP(Q(IRHS))
ENDDO
DO ILHS=1,OD_LHS%nbatches
  IF (nPrimP(ILHS) .GT. 0) CALL FREE_OVERLAP(P(ILHS))
ENDDO

CALL freeTUVitem(sharedTUV,Input)
DEALLOCATE(nPrimQ)
DEALLOCATE(Q)
DEALLOCATE(nPrimP)
DEALLOCATE(P)

DEALLOCATE(RED_GAB_LHS)
DEALLOCATE(RED_GAB_RHS)
NULLIFY(RED_GAB_LHS)
NULLIFY(RED_GAB_RHS)
DEALLOCATE(maxLHSGAB)
DEALLOCATE(maxRHSGAB)
DEALLOCATE(brashell)
NULLIFY(brashell)
DEALLOCATE(ketshell)
NULLIFY(ketshell)
DEALLOCATE(ML)
NULLIFY(ML)
DEALLOCATE(RED_DMAT_RHS)
NULLIFY(RED_DMAT_RHS)
IF(INPUT%DO_DALINK)THEN
   DEALLOCATE(RED_DMAT_LHS)
   NULLIFY(RED_DMAT_LHS)
ENDIF
CALL GETTIM(CPUTIMEEND,WALLTIMEEND)
CPUTIMELINK = CPUTIMELINK + CPUTIMEEND-CPUTIMESTART
WALLTIMELINK = WALLTIMELINK + WALLTIMEEND-WALLTIMESTART

WRITE(lupri,*)'OVERALL WALL TIMINGS '
CALL TIMTXT('>>>  WALL Time used in Link is             ',WALLTIMELINK,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  WALL Time used in merge is            ',WALLTIMEMERGE,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Exchange3 is        ',WALLTIMEEXCHANGE,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Explicit int is     ',WALLTIMEEXPLICIT,LUPRI)
WRITE(lupri,*)'----------------------------------------------------------------------   '
WALLTIMEEND = WALLTIMEEXPLICIT+WALLTIMEEXCHANGE+WALLTIMEMERGE
CALL TIMTXT('>>>  WALL Time used in Total               ',WALLTIMEEND,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'EVEN MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  WALL Time used in Build_Integrand is  ',WALLTIME_Build_Integrand,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Build_HermiteTUV is ',WALLTIME_Build_HermiteTUV,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Contract_Q is       ',WALLTIME_Contract_Q,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Contract_P is       ',WALLTIME_Contract_P,LUPRI)
WALLTIMEEND = WALLTIME_Contract_P+WALLTIME_Contract_Q+WALLTIME_Build_HermiteTUV+WALLTIME_Build_Integrand
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  WALL Time used in Total in Explicit   ',WALLTIMEEND,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'EVEN MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  WALL Time used in DistributeHermiteQ  ',WALLTIME_DistributeHermiteQ,LUPRI)
CALL TIMTXT('>>>  WALL Time used in ContractEcoeffQ     ',WALLTIME_ContractEcoeffQ,LUPRI)
CALL TIMTXT('>>>  WALL Time used in ContractBasisQ      ',WALLTIME_ContractBasisQ,LUPRI)
CALL TIMTXT('>>>  WALL Time used in SphericalTransformQ ',WALLTIME_SphericalTransformQ,LUPRI)
CALL TIMTXT('>>>  WALL Time used in AddToTUVQ           ',WALLTIME_AddToTUVQ,LUPRI)
WALLTIMEEND = WALLTIME_DistributeHermiteQ+ WALLTIME_ContractEcoeffQ+WALLTIME_ContractBasisQ+WALLTIME_SphericalTransformQ+WALLTIME_AddToTUVQ
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  WALL Time used in Total in Contract_Q ',WALLTIMEEND,LUPRI)
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  WALL Time used in DistributeHermiteP  ',WALLTIME_DistributeHermiteP,LUPRI)
CALL TIMTXT('>>>  WALL Time used in ContractEcoeffP     ',WALLTIME_ContractEcoeffP,LUPRI)
CALL TIMTXT('>>>  WALL Time used in ContractBasisP      ',WALLTIME_ContractBasisP,LUPRI)
CALL TIMTXT('>>>  WALL Time used in SphericalTransformP ',WALLTIME_SphericalTransformP,LUPRI)
CALL TIMTXT('>>>  WALL Time used in AddToPQ             ',WALLTIME_AddToPQ,LUPRI)
WALLTIMEEND = WALLTIME_DistributeHermiteP+ WALLTIME_ContractEcoeffP+WALLTIME_ContractBasisP+WALLTIME_SphericalTransformP+WALLTIME_AddToPQ
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  WALL Time used in Total in Contract_P ',WALLTIMEEND,LUPRI)
WRITE(lupri,*)'   '

WRITE(lupri,*)'OVERALL CPU TIMINGS '
CALL TIMTXT('>>>  CPU Time used in Link is             ',CPUTIMELINK,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  CPU Time used in merge is            ',CPUTIMEMERGE,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Exchange3 is        ',CPUTIMEEXCHANGE,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Explicit int is     ',CPUTIMEEXPLICIT,LUPRI)
WRITE(lupri,*)'----------------------------------------------------------------------   '
CPUTIMEEND = CPUTIMEEXPLICIT+CPUTIMEEXCHANGE+CPUTIMEMERGE
CALL TIMTXT('>>>  CPU Time used in Total               ',CPUTIMEEND,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'EVEN MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  CPU Time used in Build_Integrand is  ',CPUTIME_Build_Integrand,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Build_HermiteTUV is ',CPUTIME_Build_HermiteTUV,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Contract_Q is       ',CPUTIME_Contract_Q,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Contract_P is       ',CPUTIME_Contract_P,LUPRI)
CPUTIMEEND = CPUTIME_Contract_P+CPUTIME_Contract_Q+CPUTIME_Build_HermiteTUV+CPUTIME_Build_Integrand
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  CPU Time used in Total in Explicit   ',CPUTIMEEND,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'EVEN MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  CPU Time used in DistributeHermiteQ  ',CPUTIME_DistributeHermiteQ,LUPRI)
CALL TIMTXT('>>>  CPU Time used in ContractEcoeffQ     ',CPUTIME_ContractEcoeffQ,LUPRI)
CALL TIMTXT('>>>  CPU Time used in ContractBasisQ      ',CPUTIME_ContractBasisQ,LUPRI)
CALL TIMTXT('>>>  CPU Time used in SphericalTransformQ ',CPUTIME_SphericalTransformQ,LUPRI)
CALL TIMTXT('>>>  CPU Time used in AddToTUVQ           ',CPUTIME_AddToTUVQ,LUPRI)
CPUTIMEEND = CPUTIME_DistributeHermiteQ+ CPUTIME_ContractEcoeffQ+CPUTIME_ContractBasisQ+CPUTIME_SphericalTransformQ+CPUTIME_AddToTUVQ
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  CPU Time used in Total in Contract_Q ',CPUTIMEEND,LUPRI)
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  CPU Time used in DistributeHermiteP  ',CPUTIME_DistributeHermiteP,LUPRI)
CALL TIMTXT('>>>  CPU Time used in ContractEcoeffP     ',CPUTIME_ContractEcoeffP,LUPRI)
CALL TIMTXT('>>>  CPU Time used in ContractBasisP      ',CPUTIME_ContractBasisP,LUPRI)
CALL TIMTXT('>>>  CPU Time used in SphericalTransformP ',CPUTIME_SphericalTransformP,LUPRI)
CALL TIMTXT('>>>  CPU Time used in AddToPQ             ',CPUTIME_AddToPQ,LUPRI)
CPUTIMEEND = CPUTIME_DistributeHermiteP+ CPUTIME_ContractEcoeffP+CPUTIME_ContractBasisP+CPUTIME_SphericalTransformP+CPUTIME_AddToPQ
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  CPU Time used in Total in Contract_P ',CPUTIMEEND,LUPRI)
WRITE(lupri,*)'   '
END SUBROUTINE LINK

!*************************************************************************
!*
!*                 Link with Pass based integral driver    JCP 109,1663 
!*
!*************************************************************************
SUBROUTINE LINKwPASS(OD_LHS,OD_RHS,ODmat,ODRESmat,INPUT,OUTPUT,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
TYPE(ODmatrixitem)   :: ODmat,ODRESmat(:)
!
Integer               :: ILHS,IRHS,Start_RHS,End_RHS,i,j
TYPE(TUVitem),target  :: SharedTUV
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)       :: Alloc
TYPE(Overlap)         :: P
TYPE(Overlap)         :: Q,PassQ
TYPE(Integrand)       :: PQ
Integer,pointer       :: nPrimP(:)
Integer,pointer       :: nPrimQ(:)
Integer,pointer       :: ODpassesIndex(:)
Integer               :: nPassTypes
Integer,pointer :: Maxpassfortype(:)
Real(realk),pointer   :: RED_GAB_LHS(:,:)
Real(realk),pointer   :: RED_GAB_RHS(:,:)
Real(realk),pointer   :: RED_DMAT_LHS(:,:,:)
Real(realk),pointer   :: RED_DMAT_RHS(:,:,:)
Integer               :: dim1,dim2,dim3,dim4,A,B,C,D,idmat,endC
Integer               :: LISTSIZE!,LISTSIZE1,LISTSIZE2
Integer               :: nB,nC,nD,OLDI,RHSINDEX,startB,startD,I2,JK
Integer,allocatable   :: LIST(:)!,LIST2(:),LIST1(:)
LOGICAL,allocatable   :: DoINT(:)
real(realk),allocatable :: maxLHSGAB(:),maxRHSGAB(:)
Integer                 :: LWORK,WORKLENGTH
real(realk),allocatable :: WORK(:)
!real(realk),allocatable :: SORTING(:)
!integer,allocatable     :: bshell(:)
real(realk)           :: maxLHSELM,maxRHSELM,DMATELM1,DMATELM2,MAXDMAT,CPUTIME1,CPUTIME2
!LOGICAL               :: A_LOOP_DONE,B_LOOP_DONE,C_LOOP_DONE,D_LOOP_DONE,UNIQUE
TYPE(LINKshell),pointer  :: brashell(:),ketshell(:),ML(:,:)
real(realk)           :: CPUTimeLINK,WALLTIME1,WALLTIME2,WALLTIMELINK
real(realk)           :: CPUTimeEXCHANGE,WALLTimeEXCHANGE
real(realk)           :: CPUTimeEXPLICIT,WALLTimeEXPLICIT
real(realk)           :: CPUTIMESTART,WALLTIMESTART,CPUTIMEEND,WALLTIMEEND
logical               :: dopasses,SET_ORBITAL
integer               :: maxpasses,iPassType,numpasses,startc
integer,allocatable   :: batchindex1(:),batchindex2(:)
real(realk) :: CPUTimeMERGE,WALLTimeMERGE
logical        :: NOELEMENTSADDED

CALL GETTIM(CPUTIMESTART,WALLTIMESTART)

CPUTIMELINK = 0.D0
WALLTIMELINK = 0.D0
CPUTIMEEXCHANGE = 0.D0
WALLTIMEEXCHANGE = 0.D0
CPUTIMEEXPLICIT = 0.D0
WALLTIMEEXPLICIT = 0.D0
CPUTIME_Build_Integrand = 0.D0
WALLTIME_Build_Integrand = 0.D0
CPUTIME_Build_HermiteTUV = 0.D0
WALLTIME_Build_HermiteTUV = 0.D0
CPUTIME_Contract_Q = 0.D0
WALLTIME_Contract_Q = 0.D0
CPUTIME_Contract_P = 0.D0
WALLTIME_Contract_P = 0.D0
CPUTIMEMERGE = 0.D0
WALLTIMEMERGE = 0.D0
CPUTIME_DistributeHermiteQ = 0.D0
CPUTIME_ContractEcoeffQ = 0.D0
CPUTIME_ContractBasisQ = 0.D0
CPUTIME_SphericalTransformQ = 0.D0
CPUTIME_AddToTUVQ = 0.D0
WALLTIME_DistributeHermiteQ = 0.D0
WALLTIME_ContractEcoeffQ = 0.D0
WALLTIME_ContractBasisQ = 0.D0
WALLTIME_SphericalTransformQ = 0.D0
WALLTIME_AddToTUVQ = 0.D0

CPUTIME_DistributeHermiteP = 0.D0
CPUTIME_ContractEcoeffP = 0.D0
CPUTIME_ContractBasisP = 0.D0
CPUTIME_SphericalTransformP = 0.D0
CPUTIME_AddToPQ = 0.D0
WALLTIME_DistributeHermiteP = 0.D0
WALLTIME_ContractEcoeffP = 0.D0
WALLTIME_ContractBasisP = 0.D0
WALLTIME_SphericalTransformP = 0.D0
WALLTIME_AddToPQ = 0.D0
IF(.NOT. INPUT%CS_SCREEN)THEN
   WRITE(LUPRI,*)'Link requires screening, but you have turned off screening'
   CALL QUIT('Link requires screening, but you have turned off screening')
ENDIF

IF (IPRINT.GT.5) THEN
  IF(INPUT%DO_DALINK)THEN
     CALL LSHEADER(LUPRI,'DaLinK')
     WRITE(LUPRI,'(1X,A)') ''
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') '***                         DaLINK                             ***'
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') ''
  ELSE
     CALL LSHEADER(LUPRI,'LinK')
     WRITE(LUPRI,'(1X,A)') ''
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') '***                           LINK                             ***'
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') ''
  ENDIF
ENDIF

!REDUCE THE GAB MATRIX TO AOBATCH INDEXES

dim1=INPUT%AO(1)%p%nbatches
dim2=INPUT%AO(2)%p%nbatches
dim3=INPUT%AO(3)%p%nbatches
dim4=INPUT%AO(4)%p%nbatches
NULLIFY(RED_GAB_LHS)
NULLIFY(RED_GAB_RHS)
ALLOCATE(RED_GAB_LHS(dim1,dim2))
ALLOCATE(RED_GAB_RHS(dim3,dim4))
CALL DZERO(RED_GAB_LHS,dim1*dim2)
CALL DZERO(RED_GAB_RHS,dim3*dim4)
CALL REDUCE_MAT(LUPRI,INPUT%GAB_LHS,INPUT%AOdim(1),INPUT%AOdim(2),&
     &RED_GAB_LHS,dim1,dim2,INPUT,'LHS',.TRUE.)
CALL REDUCE_MAT(LUPRI,INPUT%GAB_RHS,INPUT%AOdim(3),INPUT%AOdim(4),&
     &RED_GAB_RHS,dim3,dim4,INPUT,'RHS',.TRUE.)

!CALCULATE MAXGAB VECTOR AND MAX ELEMENT IN GAB
ALLOCATE(maxLHSGAB(dim1))
CALL MAXGABELM(RED_GAB_LHS,dim1,dim2,maxLHSGAB,maxLHSELM)
ALLOCATE(maxRHSGAB(dim3))
CALL MAXGABELM(RED_GAB_RHS,dim3,dim4,maxRHSGAB,maxRHSELM)

! determine significant bra shell pairs

NULLIFY(brashell)
ALLOCATE(brashell(dim1))
CALL DETERMINE_SHELL_PAIRS(dim1,dim2,RED_GAB_LHS,INPUT%CS_THRESHOLD,&
        &maxRHSELM,brashell,INPUT%sameLHSaos)

NULLIFY(ketshell)
ALLOCATE(ketshell(dim3))
CALL DETERMINE_SHELL_PAIRS(dim3,dim4,RED_GAB_RHS,INPUT%CS_THRESHOLD,&
        &maxLHSELM,ketshell,.FALSE.)

NULLIFY(RED_DMAT_RHS)
ALLOCATE(RED_DMAT_RHS(dim3,dim4,Input%NDMAT_RHS))
CALL DZERO(RED_DMAT_RHS,dim3*dim4*Input%NDMAT_RHS)
DO idmat=1,Input%NDMAT_RHS
   CALL REDUCE_MAT(LUPRI,INPUT%DMAT_RHS(:,:,idmat),INPUT%AOdim(3),&
        &INPUT%AOdim(4),RED_DMAT_RHS(:,:,idmat),dim3,dim4,INPUT,'RHS',.FALSE.)
ENDDO

IF(INPUT%DO_DALINK)THEN
   NULLIFY(RED_DMAT_LHS)
   ALLOCATE(RED_DMAT_LHS(dim1,dim2,Input%NDMAT_LHS))
   CALL DZERO(RED_DMAT_LHS,dim1*dim2*Input%NDMAT_LHS)
   DO idmat=1,Input%NDMAT_LHS
      CALL REDUCE_MAT(LUPRI,INPUT%DMAT_LHS(:,:,idmat),INPUT%AOdim(3),&
           &INPUT%AOdim(4),RED_DMAT_LHS(:,:,idmat),dim1,dim2,INPUT,'LHS',.FALSE.)
   ENDDO
ENDIF

NULLIFY(ML)
ALLOCATE(ML(dim1,Input%NDMAT_RHS))
CALL DETERMINE_BRAKET_PAIRS(dim1,dim2,dim3,Input%NDMAT_RHS,RED_DMAT_RHS,&
     &maxLHSGAB,maxRHSGAB,INPUT%CS_THRESHOLD,ML)

IF(INPUT%DO_DALINK)THEN
   ALLOCATE(batchindex1(OD_RHS%nbatches))
   ALLOCATE(batchindex2(OD_RHS%nbatches))
   IF(INPUT%sameLHSaos) THEN
      DO C=1,dim3
         DO nD=1,brashell(C)%DIM
            D=ketshell(C)%elms(nD)
            IRHS=(C-1)*dim1-(C*(C-1)/2)+D 
            batchindex1(IRHS)=C
            batchindex2(IRHS)=D
         ENDDO
      ENDDO
   ELSE
      DO C=1,dim3
         DO nD=1,brashell(C)%DIM
            D=ketshell(C)%elms(nD)
            IRHS=(C-1)*dim1+D
            batchindex1(IRHS)=C
            batchindex2(IRHS)=D
         ENDDO
      ENDDO
   ENDIF
ENDIF

CALL initTUVitem(sharedTUV,Integral,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
doPasses  = INPUT%DO_PASSES
maxPasses = INPUT%maxPasses
CALL initAlloc(Alloc,LUPRI,IPRINT,'Both')
IF(doPasses) WRITE(LUPRI,*)'Integrals are calculated using passes with maxpasses=',maxPasses

NULLIFY(ODpassesIndex)
ALLOCATE(ODpassesIndex(OD_RHS%nbatches))
NULLIFY(nPrimQ)
ALLOCATE(nPrimQ(OD_RHS%nbatches)) 
CALL SelectODPassTypesFromODbatch(ODpassesIndex,nPrimQ,OD_RHS,OD_RHS%nbatches,nPassTypes,&
     &maxpasses,maxpassfortype,INPUT,IPRINT,LUPRI,'RHS')
CALL SET_ALLOC2(Alloc,Input,OD_RHS,'RHS',IPRINT,LUPRI,OD_RHS%nbatches,ODpassesIndex,nPassTypes,maxPassfortype,maxpasses)
DEALLOCATE(maxPassfortype)
NULLIFY(maxPassfortype)

NULLIFY(nPrimP)
ALLOCATE(nPrimP(OD_LHS%nbatches)) 
DO ILHS = 1,OD_LHS%nbatches
   CALL GET_NPRIM_FROM_ODBATCH(nPrimP(ILHS),OD_LHS%BATCH(ILHS),INPUT,'LHS')  
ENDDO

CALL SET_ALLOC(Alloc,Input,OD_LHS,'LHS',IPRINT,LUPRI)
idmat = 1

!$OMP PARALLEL PRIVATE(integral,PQ,IRHS,RHSINDEX,A,nB,B,ILHS,I,nC,C,OLDI,nD,D,&
!$OMP LIST,LISTSIZE,WORK,LWORK,WORKLENGTH,STARTC,DoINT,NOELEMENTSADDED,DMATELM1,DMATELM2,MAXDMAT)

integral%TUV => sharedTUV
CALL INIT_OVERLAP(P,Alloc,Input,OD_LHS,'LHS',1,IPRINT,LUPRI)
CALL INIT_OVERLAP(Q,Alloc,Input,OD_RHS,'RHS',2,IPRINT,LUPRI)
CALL allocateIntegrals(PQ,Integral,Alloc,maxpasses)
CALL INIT_PASS(PassQ,Alloc,Input,OD_RHS,'RHS',2,maxPasses,IPRINT,LUPRI)
WORKLENGTH = MAX(Alloc%maxprimRHS*Alloc%maxContRHS*Alloc%maxContRHS,&
     & Alloc%maxprimLHS*Alloc%maxContLHS*Alloc%maxContLHS)
WORKLENGTH = MAX(WORKLENGTH,Alloc%maxijkLHS*Alloc%maxijkLHS,Alloc%maxijkRHS*Alloc%maxijkRHS)
LWORK = 1
ALLOCATE(WORK(WORKLENGTH))
ALLOCATE(LIST(dim3*dim4))
ALLOCATE(DoINT(dim3*dim4))

!$OMP DO SCHEDULE(DYNAMIC,1)
DO A=1,dim1
   DO nB=1,brashell(A)%DIM
      B=brashell(A)%elms(nB)
      ILHS=(A-1)*dim1+B
      IF(INPUT%sameLHSaos) ILHS=(A-1)*dim1-(A*(A-1)/2)+B !this refer to OVERLAP which has not yet been set up
      IF(nPrimP(ILHS) .GT. 0)THEN
         !CREATING DoINT LIST from the AC element of the density matrix
         DO IRHS=1,ILHS
            DoINT(IRHS) = .FALSE.
         ENDDO
         DO nC=1,ML(A,idmat)%DIM
            C = ML(A,idmat)%elms(nC) 
            NOELEMENTSADDED = .TRUE.
            STARTC = (C-1)*dim1-(C*(C-1)/2)
            IF(INPUT%sameRHSaos)THEN
               DO nD=1,ketshell(C)%DIM
                  D=ketshell(C)%elms(nD)
                  IF(RED_DMAT_RHS(A,C,idmat)*RED_GAB_RHS(A,B)&
                       &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD )EXIT
                  NOELEMENTSADDED = .FALSE.
                  IF(D .GE. C)THEN
                     IRHS = STARTC+D
                     DoINT(IRHS) = .TRUE.
                  ELSE
                     IRHS=(D-1)*dim1-(D*(D-1)/2)+C
                     DoINT(IRHS) = .TRUE.
                  ENDIF
               ENDDO
            ELSE
               DO nD=1,ketshell(C)%DIM
                  D=ketshell(C)%elms(nD)
                  IF(RED_DMAT_RHS(A,C,idmat)*RED_GAB_RHS(A,B)&
                       &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD )EXIT
                  NOELEMENTSADDED = .FALSE.
                  DoINT((C-1)*dim1+D)=.TRUE.
               ENDDO
            ENDIF
            IF(NOELEMENTSADDED)EXIT !NO ELEMENTS ADDED IN D LOOP
         ENDDO
        !Add to DoINT from the BC element of the density matrix
         DO nC=1,ML(A,idmat)%DIM
            C = ML(A,idmat)%elms(nC)
            NOELEMENTSADDED = .TRUE.
            STARTC=(C-1)*dim1-(C*(C-1)/2)
            IF(INPUT%sameRHSaos)THEN
               DO nD=1,ketshell(C)%DIM
                  D=ketshell(C)%elms(nD)
                  IF(RED_DMAT_RHS(B,C,idmat)*RED_GAB_RHS(A,B)&
                       &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD ) EXIT
                  NOELEMENTSADDED = .FALSE.
                  IF(D .GE. C)THEN
                     IRHS = (C-1)*dim1-(C*(C-1)/2)+D
                     DoINT(IRHS)=.TRUE.
                  ELSE
                     IRHS = (D-1)*dim1-(D*(D-1)/2)+C
                     DoINT(IRHS)=.TRUE.
                  ENDIF
               ENDDO
            ELSE
               DO nD=1,ketshell(C)%DIM
                  D=ketshell(C)%elms(nD)
                  IF(RED_DMAT_RHS(B,C,idmat)*RED_GAB_RHS(A,B)&
                       &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD ) EXIT
                  NOELEMENTSADDED = .FALSE.
                  Doint((C-1)*dim1+D)=.TRUE.
               ENDDO
            ENDIF
            IF(NOELEMENTSADDED)EXIT !NO ELEMENTS ADDED IN D LOOP
         ENDDO
         CALL GETTIM(CPUTIME1,WALLTIME1)
         CALL BUILD_LIST(LUPRI,LIST,DoINT,LISTSIZE,dim3*dim4,ILHS)
         CALL GETTIM(CPUTIME2,WALLTIME2)
         CPUTIMEMERGE = CPUTIMEMERGE + CPUTIME2-CPUTIME1
         WALLTIMEMERGE = WALLTIMEMERGE + WALLTIME2-WALLTIME1
         IF(LISTSIZE .GT. 0)THEN
            CALL SET_INITIALIZED_Overlap(P,nPrimP(ILHS),Input,SharedTUV,Integral,Alloc,&
                 &OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS',.TRUE.,.TRUE.)
            DO iPassType=1,nPassTypes
               numPasses = 0
               SET_ORBITAL = .TRUE.
               DO RHSINDEX=1,LISTSIZE
                  IRHS=LIST(RHSINDEX)
                  IF (ODpassesIndex(IRHS).EQ.iPassType) THEN
                     IF(nPrimQ(IRHS) .GT. 0)THEN
                        IF(INPUT%DO_DALINK )THEN
                           !EXTRA DENSITY ACCELERATED SCREENING
                           C=batchindex1(IRHS)
                           D=batchindex2(IRHS)
                           DMATELM1=RED_DMAT_RHS(A,C,idmat)*RED_DMAT_LHS(B,D,idmat)
                           DMATELM2=RED_DMAT_RHS(B,C,idmat)*RED_DMAT_LHS(A,D,idmat)
                           MAXDMAT=MAX(DMATELM1,DMATELM2)
                           IF(MAXDMAT*RED_GAB_LHS(A,B)*RED_GAB_RHS(C,D) .GE. INPUT%CS_THRESHOLD )THEN 
                              !DO NOT NEED TO SET ORBITALS IF FIRST
                              IF (doPasses .AND..NOT. ((ILHS.EQ.IRHS).AND.INPUT%sameODs)) THEN
                                 CALL SET_INITIALIZED_Overlap(Q,nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,&
                                      &OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS',.TRUE.,SET_ORBITAL)
                                 SET_ORBITAL = .FALSE.
                                 numPasses = numPasses + 1
                                 CALL AddOverlapToPass(PassQ,Q,numPasses,maxPasses,LUPRI,IPRINT)
                                 IF (numPasses.EQ.maxPasses) THEN
                                    CALL GETTIM(CPUTIME1,WALLTIME1)
                                    CALL ExplicitIntegrals(Integral,PQ,P,PassQ,INPUT,OUTPUT,&
                                         & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                                    CALL GETTIM(CPUTIME2,WALLTIME2)
                                    CPUTIMEEXPLICIT = CPUTIMEEXPLICIT + CPUTIME2-CPUTIME1
                                    WALLTIMEEXPLICIT = WALLTIMEEXPLICIT + WALLTIME2-WALLTIME1
                                    !$OMP CRITICAL
                                    IF(INPUT%DRHS_SYM)THEN
                                       CALL Distribute_Exchange3SYM(INTEGRAL,PQ,Integral%PQ,&
                                            & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                                            & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                                    ELSE
                                       CALL Distribute_Exchange3(INTEGRAL,PQ,Integral%PQ,&
                                            & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                                            & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                                    ENDIF
                                    !$OMP END CRITICAL
                                    CALL GETTIM(CPUTIME1,WALLTIME1)
                                    CPUTIMEEXCHANGE = CPUTIMEEXCHANGE + CPUTIME1-CPUTIME2
                                    WALLTIMEEXCHANGE = WALLTIMEEXCHANGE + WALLTIME1-WALLTIME2
                                    numPasses = 0
                                 ENDIF
                              ELSE
                                 CALL SET_INITIALIZED_Overlap(Q,nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,&
                                      &OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS',.TRUE.,.TRUE.)
                                 CALL GETTIM(CPUTIME1,WALLTIME1)
                                 CALL ExplicitIntegrals(Integral,PQ,P,Q,INPUT,OUTPUT,&
                                      & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                                 CALL GETTIM(CPUTIME2,WALLTIME2)
                                 CPUTIMEEXPLICIT = CPUTIMEEXPLICIT + CPUTIME2-CPUTIME1
                                 WALLTIMEEXPLICIT = WALLTIMEEXPLICIT + WALLTIME2-WALLTIME1
                                 !$OMP CRITICAL
                                 IF(INPUT%DRHS_SYM)THEN
                                    CALL Distribute_Exchange3SYM(INTEGRAL,PQ,Integral%PQ,&
                                         & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                                         & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                                 ELSE
                                    CALL Distribute_Exchange3(INTEGRAL,PQ,Integral%PQ,&
                                         & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                                         & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                                 ENDIF
                                    !$OMP END CRITICAL
                                 CALL GETTIM(CPUTIME1,WALLTIME1)
                                 CPUTIMEEXCHANGE = CPUTIMEEXCHANGE + CPUTIME1-CPUTIME2
                                 WALLTIMEEXCHANGE = WALLTIMEEXCHANGE + WALLTIME1-WALLTIME2
                              ENDIF
                           ENDIF
                        ELSE
                           IF (doPasses .AND..NOT. ((ILHS.EQ.IRHS).AND.INPUT%sameODs)) THEN
                              CALL SET_INITIALIZED_Overlap(Q,nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,&
                                   &OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS',.TRUE.,SET_ORBITAL)
                              SET_ORBITAL = .FALSE.
                              numPasses = numPasses + 1
                              CALL AddOverlapToPass(PassQ,Q,numPasses,maxPasses,LUPRI,IPRINT)
                              IF (numPasses.EQ.maxPasses) THEN
                                 CALL GETTIM(CPUTIME1,WALLTIME1)
                                 CALL ExplicitIntegrals(Integral,PQ,P,PassQ,INPUT,OUTPUT,&
                                      & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                                 CALL GETTIM(CPUTIME2,WALLTIME2)
                                 CPUTIMEEXPLICIT = CPUTIMEEXPLICIT + CPUTIME2-CPUTIME1
                                 WALLTIMEEXPLICIT = WALLTIMEEXPLICIT + WALLTIME2-WALLTIME1
                                 !$OMP CRITICAL
                                 IF(INPUT%DRHS_SYM)THEN
                                    CALL Distribute_Exchange3SYM(INTEGRAL,PQ,Integral%PQ,&
                                         & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                                         & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                                 ELSE
                                    CALL Distribute_Exchange3(INTEGRAL,PQ,Integral%PQ,&
                                         & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                                         & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                                 ENDIF
                                 !$OMP END CRITICAL
                                 CALL GETTIM(CPUTIME1,WALLTIME1)
                                 CPUTIMEEXCHANGE = CPUTIMEEXCHANGE + CPUTIME1-CPUTIME2
                                 WALLTIMEEXCHANGE = WALLTIMEEXCHANGE + WALLTIME1-WALLTIME2
                                 numPasses = 0
                              ENDIF
                           ELSE
                              CALL SET_INITIALIZED_Overlap(Q,nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,&
                                   &OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS',.TRUE.,.TRUE.)
                              CALL GETTIM(CPUTIME1,WALLTIME1)
                              CALL ExplicitIntegrals(Integral,PQ,P,Q,INPUT,OUTPUT,&
                                   & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                              CALL GETTIM(CPUTIME2,WALLTIME2)
                              CPUTIMEEXPLICIT = CPUTIMEEXPLICIT + CPUTIME2-CPUTIME1
                              WALLTIMEEXPLICIT = WALLTIMEEXPLICIT + WALLTIME2-WALLTIME1
                              !$OMP CRITICAL
                              IF(INPUT%DRHS_SYM)THEN
                                 CALL Distribute_Exchange3SYM(INTEGRAL,PQ,Integral%PQ,&
                                      & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                                      & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                              ELSE
                                 CALL Distribute_Exchange3(INTEGRAL,PQ,Integral%PQ,&
                                      & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                                      & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                              ENDIF
                              !$OMP END CRITICAL
                              CALL GETTIM(CPUTIME1,WALLTIME1)
                              CPUTIMEEXCHANGE = CPUTIMEEXCHANGE + CPUTIME1-CPUTIME2
                              WALLTIMEEXCHANGE = WALLTIMEEXCHANGE + WALLTIME1-WALLTIME2
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
               IF(numPasses.GT.0) THEN
                  CALL GETTIM(CPUTIME1,WALLTIME1)
                  CALL ExplicitIntegrals(Integral,PQ,P,PassQ,INPUT,OUTPUT,&
                       & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                  CALL GETTIM(CPUTIME2,WALLTIME2)
                  CPUTIMEEXPLICIT = CPUTIMEEXPLICIT + CPUTIME2-CPUTIME1
                  WALLTIMEEXPLICIT = WALLTIMEEXPLICIT + WALLTIME2-WALLTIME1
                  !$OMP CRITICAL
                  IF(INPUT%DRHS_SYM)THEN
                     CALL Distribute_Exchange3SYM(INTEGRAL,PQ,Integral%PQ,&
                          & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                          & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                  ELSE
                     CALL Distribute_Exchange3(INTEGRAL,PQ,Integral%PQ,&
                          & PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals,&
                          & INPUT,ODmat,ODRESmat,OUTPUT,LUPRI,IPRINT) 
                  ENDIF
                  !$OMP END CRITICAL
                  CALL GETTIM(CPUTIME1,WALLTIME1)
                  CPUTIMEEXCHANGE = CPUTIMEEXCHANGE + CPUTIME1-CPUTIME2
                  WALLTIMEEXCHANGE = WALLTIMEEXCHANGE + WALLTIME1-WALLTIME2
               ENDIF
            ENDDO
         ENDIF
      ENDIF ! PS-screening
   ENDDO
ENDDO
!$OMP END DO 

DEALLOCATE(LIST)
DEALLOCATE(DoINT)
CALL deallocateIntegrals(PQ,Integral)
CALL FREE_OVERLAP(P)
CALL FREE_OVERLAP(Q)
CALL FreePass(PassQ)
DEALLOCATE(WORK)

!$OMP END PARALLEL

DEALLOCATE(ODpassesIndex)
DEALLOCATE(nPrimQ)
DEALLOCATE(nPrimP)
CALL freeTUVitem(sharedTUV,Input)

DEALLOCATE(RED_GAB_LHS)
DEALLOCATE(RED_GAB_RHS)
NULLIFY(RED_GAB_LHS)
NULLIFY(RED_GAB_RHS)
DEALLOCATE(maxLHSGAB)
DEALLOCATE(maxRHSGAB)
DEALLOCATE(brashell)
NULLIFY(brashell)
DEALLOCATE(ketshell)
NULLIFY(ketshell)
DEALLOCATE(ML)
NULLIFY(ML)
DEALLOCATE(RED_DMAT_RHS)
NULLIFY(RED_DMAT_RHS)
IF(INPUT%DO_DALINK)THEN
   DEALLOCATE(RED_DMAT_LHS)
   NULLIFY(RED_DMAT_LHS)
ENDIF

CALL GETTIM(CPUTIMEEND,WALLTIMEEND)
CPUTIMELINK = CPUTIMELINK + CPUTIMEEND-CPUTIMESTART
WALLTIMELINK = WALLTIMELINK + WALLTIMEEND-WALLTIMESTART


WRITE(lupri,*)'OVERALL WALL TIMINGS '
CALL TIMTXT('>>>  WALL Time used in LinkwPASS is        ',WALLTIMELINK,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  WALL Time used in merge is            ',WALLTIMEMERGE,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Exchange3 is        ',WALLTIMEEXCHANGE,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Explicit int is     ',WALLTIMEEXPLICIT,LUPRI)
WRITE(lupri,*)'----------------------------------------------------------------------   '
WALLTIMEEND = WALLTIMEEXPLICIT+WALLTIMEEXCHANGE+WALLTIMEMERGE
CALL TIMTXT('>>>  WALL Time used in Total               ',WALLTIMEEND,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'EVEN MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  WALL Time used in Build_Integrand is  ',WALLTIME_Build_Integrand,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Build_HermiteTUV is ',WALLTIME_Build_HermiteTUV,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Contract_Q is       ',WALLTIME_Contract_Q,LUPRI)
CALL TIMTXT('>>>  WALL Time used in Contract_P is       ',WALLTIME_Contract_P,LUPRI)
WALLTIMEEND = WALLTIME_Contract_P+WALLTIME_Contract_Q+WALLTIME_Build_HermiteTUV+WALLTIME_Build_Integrand
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  WALL Time used in Total in Explicit   ',WALLTIMEEND,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'EVEN MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  WALL Time used in DistributeHermiteQ  ',WALLTIME_DistributeHermiteQ,LUPRI)
CALL TIMTXT('>>>  WALL Time used in ContractEcoeffQ     ',WALLTIME_ContractEcoeffQ,LUPRI)
CALL TIMTXT('>>>  WALL Time used in ContractBasisQ      ',WALLTIME_ContractBasisQ,LUPRI)
CALL TIMTXT('>>>  WALL Time used in SphericalTransformQ ',WALLTIME_SphericalTransformQ,LUPRI)
CALL TIMTXT('>>>  WALL Time used in AddToTUVQ           ',WALLTIME_AddToTUVQ,LUPRI)
WALLTIMEEND = WALLTIME_DistributeHermiteQ+ WALLTIME_ContractEcoeffQ+WALLTIME_ContractBasisQ+WALLTIME_SphericalTransformQ+WALLTIME_AddToTUVQ
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  WALL Time used in Total in Contract_Q ',WALLTIMEEND,LUPRI)
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  WALL Time used in DistributeHermiteP  ',WALLTIME_DistributeHermiteP,LUPRI)
CALL TIMTXT('>>>  WALL Time used in ContractEcoeffP     ',WALLTIME_ContractEcoeffP,LUPRI)
CALL TIMTXT('>>>  WALL Time used in ContractBasisP      ',WALLTIME_ContractBasisP,LUPRI)
CALL TIMTXT('>>>  WALL Time used in SphericalTransformP ',WALLTIME_SphericalTransformP,LUPRI)
CALL TIMTXT('>>>  WALL Time used in AddToPQ             ',WALLTIME_AddToPQ,LUPRI)
WALLTIMEEND = WALLTIME_DistributeHermiteP+ WALLTIME_ContractEcoeffP+WALLTIME_ContractBasisP+WALLTIME_SphericalTransformP+WALLTIME_AddToPQ
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  WALL Time used in Total in Contract_P ',WALLTIMEEND,LUPRI)
WRITE(lupri,*)'   '

WRITE(lupri,*)'OVERALL CPU TIMINGS '
CALL TIMTXT('>>>  CPU Time used in Link is             ',CPUTIMELINK,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  CPU Time used in merge is            ',CPUTIMEMERGE,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Exchange3 is        ',CPUTIMEEXCHANGE,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Explicit int is     ',CPUTIMEEXPLICIT,LUPRI)
WRITE(lupri,*)'----------------------------------------------------------------------   '
CPUTIMEEND = CPUTIMEEXPLICIT+CPUTIMEEXCHANGE+CPUTIMEMERGE
CALL TIMTXT('>>>  CPU Time used in Total               ',CPUTIMEEND,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'EVEN MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  CPU Time used in Build_Integrand is  ',CPUTIME_Build_Integrand,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Build_HermiteTUV is ',CPUTIME_Build_HermiteTUV,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Contract_Q is       ',CPUTIME_Contract_Q,LUPRI)
CALL TIMTXT('>>>  CPU Time used in Contract_P is       ',CPUTIME_Contract_P,LUPRI)
CPUTIMEEND = CPUTIME_Contract_P+CPUTIME_Contract_Q+CPUTIME_Build_HermiteTUV+CPUTIME_Build_Integrand
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  CPU Time used in Total in Explicit   ',CPUTIMEEND,LUPRI)
WRITE(lupri,*)'   '
WRITE(lupri,*)'EVEN MORE EXPLICIT TIMINGS   '
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  CPU Time used in DistributeHermiteQ  ',CPUTIME_DistributeHermiteQ,LUPRI)
CALL TIMTXT('>>>  CPU Time used in ContractEcoeffQ     ',CPUTIME_ContractEcoeffQ,LUPRI)
CALL TIMTXT('>>>  CPU Time used in ContractBasisQ      ',CPUTIME_ContractBasisQ,LUPRI)
CALL TIMTXT('>>>  CPU Time used in SphericalTransformQ ',CPUTIME_SphericalTransformQ,LUPRI)
CALL TIMTXT('>>>  CPU Time used in AddToTUVQ           ',CPUTIME_AddToTUVQ,LUPRI)
CPUTIMEEND = CPUTIME_DistributeHermiteQ+ CPUTIME_ContractEcoeffQ+CPUTIME_ContractBasisQ+CPUTIME_SphericalTransformQ+CPUTIME_AddToTUVQ
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  CPU Time used in Total in Contract_Q ',CPUTIMEEND,LUPRI)
WRITE(lupri,*)'   '
CALL TIMTXT('>>>  CPU Time used in DistributeHermiteP  ',CPUTIME_DistributeHermiteP,LUPRI)
CALL TIMTXT('>>>  CPU Time used in ContractEcoeffP     ',CPUTIME_ContractEcoeffP,LUPRI)
CALL TIMTXT('>>>  CPU Time used in ContractBasisP      ',CPUTIME_ContractBasisP,LUPRI)
CALL TIMTXT('>>>  CPU Time used in SphericalTransformP ',CPUTIME_SphericalTransformP,LUPRI)
CALL TIMTXT('>>>  CPU Time used in AddToPQ             ',CPUTIME_AddToPQ,LUPRI)
CPUTIMEEND = CPUTIME_DistributeHermiteP+ CPUTIME_ContractEcoeffP+CPUTIME_ContractBasisP+CPUTIME_SphericalTransformP+CPUTIME_AddToPQ
WRITE(lupri,*)'----------------------------------------------------------------------   '
CALL TIMTXT('>>>  CPU Time used in Total in Contract_P ',CPUTIMEEND,LUPRI)
WRITE(lupri,*)'   '
END SUBROUTINE LINKWPASS

SUBROUTINE REDUCE_MAT(IUNIT,GAB,nbast1,nbast2,RED_GAB,dim1,dim2,INPUT,SIDE,POS)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
INTEGER              :: nbast1,nbast2,dim1,dim2,IUNIT
REAL(REALK)          :: GAB(nbast1,nbast2),RED_GAB(dim1,dim2)
Character*(*)        :: SIDE
!
Integer              :: IA,Ja,KA,IB,JB,KB,sA,sB,CA,CB,AOA,AOB,nKA,nKB
logical              :: POS !if all elements are positive or zero
SELECT CASE(SIDE)
 CASE('LHS')
!    sameAOs=INPUT%sameLHSaos
    AOA = 1
    AOB = 2
 CASE('RHS')
!    sameAOs=INPUT%sameRHSaos
    AOA = 3
    AOB = 4
 CASE DEFAULT
    WRITE(IUNIT,'(1X,2A)') 'Wrong case in Reduce mat =',SIDE
    CALL QUIT('Wrong case in reduce_mat')
 END SELECT

IF(POS)THEN
  DO IA=1,dim1
   DO IB=1,dim2 
    DO JA=1,INPUT%AO(AOA)%p%BATCH(IA)%nAngmom
     sA=INPUT%AO(AOA)%p%BATCH(IA)%startOrbital(JA)
     DO JB=1,INPUT%AO(AOB)%p%BATCH(IB)%nAngmom
      sB=INPUT%AO(AOB)%p%BATCH(IB)%startOrbital(JB)
      nKA=INPUT%AO(AOA)%p%BATCH(IA)%nOrbComp(JA)
      nKB=INPUT%AO(AOB)%p%BATCH(IB)%nOrbComp(JB)
      DO CA=1,INPUT%AO(AOA)%p%BATCH(IA)%nContracted(JA)
       DO CB=1,INPUT%AO(AOB)%p%BATCH(IB)%nContracted(JB)
        DO KA=1,nKA
         DO KB=1,nKB
          RED_GAB(IA,IB)=MAX(RED_GAB(IA,IB),GAB(sA+KA+(CA-1)*nKA-1,sB+KB+(CB-1)*nKB-1))
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
ELSE !WE NEED TO TAKE THE ABS() OF THE ELEMENT
  DO IA=1,dim1
   DO IB=1,dim2 
    DO JA=1,INPUT%AO(AOA)%p%BATCH(IA)%nAngmom
     sA=INPUT%AO(AOA)%p%BATCH(IA)%startOrbital(JA)
     DO JB=1,INPUT%AO(AOB)%p%BATCH(IB)%nAngmom
      sB=INPUT%AO(AOB)%p%BATCH(IB)%startOrbital(JB)
      nKA=INPUT%AO(AOA)%p%BATCH(IA)%nOrbComp(JA)
      nKB=INPUT%AO(AOB)%p%BATCH(IB)%nOrbComp(JB)
      DO CA=1,INPUT%AO(AOA)%p%BATCH(IA)%nContracted(JA)
       DO CB=1,INPUT%AO(AOB)%p%BATCH(IB)%nContracted(JB)
        DO KA=1,nKA
         DO KB=1,nKB
          RED_GAB(IA,IB)=MAX(RED_GAB(IA,IB),ABS(GAB(sA+KA+(CA-1)*nKA-1,sB+KB+(CB-1)*nKB-1)))
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
ENDIF

END SUBROUTINE REDUCE_MAT

SUBROUTINE BUILD_LIST(LUPRI,LIST,DoINT,LISTSIZE,SIZE,ILHS)
IMPLICIT NONE
Integer  :: LISTSIZE,SIZE,LUPRI,ILHS
Integer  :: LIST(SIZE)
Logical  :: DoINT(SIZE)
!
Integer  :: K,J

K=0
DO J=1,ILHS
   IF(DoINT(J))THEN
      K=K+1
      LIST(K)=J
   ENDIF
ENDDO
LISTSIZE=K

END SUBROUTINE BUILD_LIST
!!$SUBROUTINE MERGE_LIST(LUPRI,LIST,LIST1,LIST2,LISTSIZE,LISTSIZE1,LISTSIZE2,SIZE,ILHS,sameODs)
!!$IMPLICIT NONE
!!$Integer  :: SIZE,LISTSIZE1,LISTSIZE2,LISTSIZE
!!$Integer  :: LIST(SIZE),LIST1(LISTSIZE1),LIST2(LISTSIZE2),LUPRI
!!$!
!!$Integer  :: I,J,K,ILHS,TMPJ
!!$logical  :: UNIQUE,sameODs
!!$
!!$K=0
!!$DO J=1,LISTSIZE1
!!$   IF(LIST1(J) .LE. ILHS)THEN
!!$      TMPJ = LIST1(J)
!!$      UNIQUE=.TRUE.
!!$      DO I=1,K
!!$         IF(LIST(I).EQ.TMPJ) THEN
!!$            UNIQUE=.FALSE.
!!$            EXIT
!!$         ENDIF
!!$      ENDDO
!!$      IF(UNIQUE) THEN
!!$         K=K+1
!!$         LIST(K)=TMPJ
!!$      ENDIF
!!$   ENDIF
!!$ENDDO
!!$
!!$!SECOND LIST
!!$DO J=1,LISTSIZE2
!!$   IF(LIST2(J) .LE. ILHS)THEN
!!$      TMPJ = LIST2(J) 
!!$      UNIQUE=.TRUE.
!!$      DO I=1,K
!!$         IF(LIST(I).EQ.TMPJ)THEN
!!$            UNIQUE=.FALSE.
!!$            EXIT
!!$         ENDIF
!!$      ENDDO
!!$      IF(UNIQUE) THEN
!!$         K=K+1
!!$         LIST(K)=TMPJ
!!$      ENDIF
!!$   ENDIF
!!$ENDDO
!!$LISTSIZE=K
!!$
!!$END SUBROUTINE MERGE_LIST
SUBROUTINE MAXGABELM(RED_GAB,dim1,dim2,maxGAB,maxELM)
IMPLICIT NONE
INTEGER      :: dim1,dim2 
REAL(REALK)  :: maxELM,RED_GAB(dim1,dim2),maxGAB(dim1)
!
INTEGER      :: A,B

maxELM=0
DO A=1,dim1
   maxGAB(A)= 0
   DO B=1,dim2
      maxGAB(A)=MAX(maxGAB(A),RED_GAB(A,B)) 
   ENDDO
   maxELM=MAX(MAXELM,maxGAB(A))
ENDDO

END SUBROUTINE MAXGABELM

SUBROUTINE DETERMINE_SHELL_PAIRS(dim1,dim2,RED_GAB_LHS,CS_THRESHOLD,maxRHSELM,brashell,sameAOs)
IMPLICIT NONE
integer                 :: dim1,dim2
real(realk)             :: RED_GAB_LHS(dim1,dim2),CS_THRESHOLD,maxRHSELM
TYPE(LINKshell)         :: brashell(dim1)
logical                 :: sameAOs
!
real(realk),allocatable :: SORTING(:)
integer,allocatable     :: bshell(:)
integer                 :: A,startB,I,B,C

ALLOCATE(SORTING(dim2))
ALLOCATE(bshell(dim2))
DO A=1,dim1
   I=1
   StartB=1
   IF(sameAOs)StartB=A
   DO B=StartB,dim2
      IF(RED_GAB_LHS(A,B) .GT. CS_THRESHOLD/maxRHSELM)THEN
         bshell(I)=B
         SORTING(I)=RED_GAB_LHS(A,B)
         I=I+1
      ENDIF
   ENDDO
   IF(I.EQ.1)THEN
      brashell(A)%DIM=0
      NULLIFY(brashell(A)%elms)
      ALLOCATE(brashell(A)%elms(I))
   ELSE
      brashell(A)%DIM=I-1
      CALL Qsort(SORTING(1:I-1),bshell(1:I-1))
      NULLIFY(brashell(A)%elms)
      ALLOCATE(brashell(A)%elms(I-1))
      DO C=1,I-1
         brashell(A)%elms(C)=bSHELL(C)
      ENDDO
   ENDIF
ENDDO
DEALLOCATE(bshell)
DEALLOCATE(SORTING)
END SUBROUTINE DETERMINE_SHELL_PAIRS

SUBROUTINE DETERMINE_BRAKET_PAIRS(dim1,dim2,dim3,NDMAT,RED_DMAT_RHS,maxLHSGAB,maxRHSGAB,CS_THRESHOLD,ML)
IMPLICIT NONE
integer                 :: dim1,dim2,dim3,NDMAT
real(realk)             :: RED_DMAT_RHS(dim1,dim2,NDMAT),CS_THRESHOLD,maxRHSELM
real(realk)             :: maxLHSGAB(dim1),maxRHSGAB(dim3)
TYPE(LINKshell)         :: ML(dim1,NDMAT)
logical                 :: sameAOs
!
real(realk),allocatable :: SORTING(:)
integer,allocatable     :: bshell(:)
integer                 :: A,B,C,I,idmat
ALLOCATE(SORTING(dim3))
ALLOCATE(bshell(dim3))
DO idmat=1,NDMAT
   DO A=1,dim1
      I=0
!      StartC=1
!      IF(sameODs)StartC=A
      DO C=1,A!dim3
         IF(RED_DMAT_RHS(A,C,idmat)*maxLHSGAB(A)*maxRHSGAB(C).GT. CS_THRESHOLD)THEN
            I=I+1
            bshell(I)=C
            SORTING(I)=RED_DMAT_RHS(A,C,idmat)*maxRHSGAB(C)
         ENDIF
      ENDDO
      IF(I .EQ. 0)THEN
         ML(A,idmat)%DIM=0
         NULLIFY(ML(A,idmat)%elms)
         ALLOCATE(ML(A,idmat)%elms(1))
      ELSE
         ML(A,idmat)%DIM=I
         CALL Qsort(SORTING(1:I),bshell(1:I))
         NULLIFY(ML(A,idmat)%elms)
         ALLOCATE(ML(A,idmat)%elms(I))
         DO C=1,I
            ML(A,idmat)%elms(C)=bSHELL(C)
         ENDDO
      ENDIF
   ENDDO
ENDDO
DEALLOCATE(SORTING)
DEALLOCATE(bshell)
END SUBROUTINE DETERMINE_BRAKET_PAIRS

!*************************************************************************
!*
!*                 End of Link
!*
!*************************************************************************

SUBROUTINE DistributeIntegrals(Integral,PQ,Input,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT,dimQ,dimP

dimQ=PQ%Q%p%totOrbitals
dimP=PQ%P%p%totOrbitals

!$OMP CRITICAL
IF (Input%CS_int) THEN
   IF (INPUT%sameLHSaos)THEN
      CALL DistributePQ_CS1(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT)
   ELSE
      CALL DistributePQ_CS2(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT)
   ENDIF
ELSE
   IF (Input%DO_FOCK) THEN
      IF (Input%DO_Coulomb) THEN
         IF (Input%DO_JENGINE) THEN
            IF(INPUT%sameLHSaos)THEN
               CALL Distribute_Jengine1(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT) 
            ELSE
               CALL Distribute_Jengine2(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT) 
            ENDIF
         ELSE
            IF (INPUT%sameLHSaos.AND.INPUT%sameRHSaos.AND.INPUT%sameODs)THEN
               CALL Distribute_Coulomb1(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT) 
            ELSE IF (INPUT%sameLHSaos .AND. INPUT%sameRHSaos)THEN
               CALL Distribute_Coulomb2(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT) 
            ELSE IF (INPUT%sameLHSaos)THEN
               CALL Distribute_Coulomb3(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT) 
            ELSE IF (INPUT%sameRHSaos)THEN
               CALL Distribute_Coulomb4(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT) 
            ELSE IF (INPUT%sameODs .AND. (INPUT%sameLHSaos.NEQV.INPUT%sameRHSaos))THEN
               CALL QUIT('Progamming error, sameODs but sameLHSaos.NE.sameRHSaos!')
            ELSE IF (INPUT%sameODs)THEN
               CALL Distribute_Coulomb5(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT) 
            ELSE                                
               CALL Distribute_Coulomb6(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT) 
            ENDIF
         ENDIF
      ENDIF
      IF(Input%DO_Exchange)THEN
         IF((INPUT%sameLHSaos.AND.INPUT%sameRHSaos).AND.INPUT%sameODs )THEN
            CALL Distribute_Exchange(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT) 
         ELSE
            CALL QUIT('Progamming error Distribute_Exchange, only works if all four AOs are identical!')
         ENDIF
      ENDIF
   ELSE
      ! No pre-contraction with densities (either full integrals, or contraction
      ! with full integrals)
      CALL DistributePQint(INTEGRAL,PQ,Integral%PQ,dimQ,dimP,INPUT,OUTPUT,LUPRI,IPRINT)
   ENDIF
ENDIF
!$OMP END CRITICAL

END SUBROUTINE DistributeIntegrals

SUBROUTINE Distribute_Coulomb1(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
!sameLHS = TRUE sameRHS = TRUE sameOD = TRUE
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: STARTA1,STARTB1,STARTC1,STARTD1
Integer :: iC,iD,iOrbQ,idmat
real(realk) :: JINT,JFAC,dtemp_lhs,dtemp_rhs
integer :: PASSPOFFSET,PASSQOFFSET
real(realk),allocatable :: TEMP(:,:)
integer :: iPassP,iPassQ

IF (IPRINT.GT.50) THEN
  ALLOCATE(TEMP(Output%ndim(1),Output%ndim(2)))
  CALL DCOPY(Output%ndim(1)*Output%ndim(2),Output%resultMat,1,TEMP,1)
ENDIF

JFAC=INPUT%CoulombFactor
nAngmomP = PQ%P%p%nAngmom
!nderivP=Input%nDerivP
!totOrbP = PQ%P%p%totOrbitals
!
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
   startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
   startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)
   IF(startA.NE.startB)THEN
      iOrbP=iOrbitalP-1+PASSPOFFSET
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=iAngB+STARTB1
               DO iAngA=1,nAngA
                  iA1=iAngA+STARTA1
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                     DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=STARTD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=STARTC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS
                                          dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat)+Input%DMAT_RHS(iC1,iD1,idmat)
                                          dtemp_lhs = Input%DMAT_RHS(iB1,iA1,idmat)+Input%DMAT_RHS(iA1,iB1,idmat) 
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                               & +  dtemp_rhs*Jint
                                          Output%ResultMat(iB1,iA1,1,1,idmat) = Output%ResultMat(iB1,iA1,1,1,idmat)&
                                               & +  dtemp_rhs*Jint
                                          Output%ResultMat(iC1,iD1,1,1,idmat) = Output%ResultMat(iC1,iD1,1,1,idmat)&
                                               & +  dtemp_lhs*Jint
                                          Output%ResultMat(iD1,iC1,1,1,idmat) = Output%ResultMat(iD1,iC1,1,1,idmat)&
                                               & +  dtemp_lhs*Jint
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
!                           iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE  !startC=startD
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          dtemp_lhs = Input%DMAT_RHS(iB1,iA1,idmat) + Input%DMAT_RHS(iA1,iB1,idmat)     
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                               & +  Input%DMAT_RHS(iD1,iC1,idmat)*Jint
                                          Output%ResultMat(iB1,iA1,1,1,idmat) = Output%ResultMat(iB1,iA1,1,1,idmat)&
                                               & +  Input%DMAT_RHS(iD1,iC1,idmat)*Jint
                                          Output%ResultMat(iC1,iD1,1,1,idmat) = Output%ResultMat(iC1,iD1,1,1,idmat)&
                                               & +  dtemp_lhs*Jint
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
!                           iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS
                                          dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat) + Input%DMAT_RHS(iC1,iD1,idmat)
!                                          dtemp_lhs = Input%DMAT_RHS(iB1,iA1,idmat) + Input%DMAT_RHS(iA1,iB1,idmat)  
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                               & +  dtemp_rhs*Jint
                                          Output%ResultMat(iB1,iA1,1,1,idmat) = Output%ResultMat(iB1,iA1,1,1,idmat)&
                                               & +  dtemp_rhs*Jint
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
!                           iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE  !startC = startD
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat)
!                                          dtemp_lhs = Input%DMAT_RHS(iB1,iA1,idmat) + Input%DMAT_RHS(iA1,iB1,idmat)           
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                               & +  dtemp_rhs*Jint
                                          Output%ResultMat(iB1,iA1,1,1,idmat) = Output%ResultMat(iB1,iA1,1,1,idmat)&
                                               & +  dtemp_rhs*Jint
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
!                           iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ENDIF
!                     iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
                     iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!      iOrbitalP = iOrbitalP + nOrbP
   ELSE !NOT (startA NE startB)
      iOrbP=iOrbitalP-1+PASSPOFFSET
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=startB1+iAngB
               DO iAngA=1,nAngA
                  iA1=startA1+iAngA
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                     DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET

                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat) + Input%DMAT_RHS(iC1,iD1,idmat)
                                          dtemp_lhs = Input%DMAT_RHS(iB1,iA1,idmat)
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                               & +  dtemp_rhs*Jint
                                          Output%ResultMat(iC1,iD1,1,1,idmat) = Output%ResultMat(iC1,iD1,1,1,idmat)&
                                               & +  dtemp_lhs*Jint
                                          Output%ResultMat(iD1,iC1,1,1,idmat) = Output%ResultMat(iD1,iC1,1,1,idmat)&
                                               & +  dtemp_lhs*Jint
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                        !   iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE !startC.NE.startD
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET

                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                               & +  Input%DMAT_RHS(iD1,iC1,idmat)*Jint
                                          Output%ResultMat(iC1,iD1,1,1,idmat) = Output%ResultMat(iC1,iD1,1,1,idmat)&
                                               & +  Input%DMAT_RHS(iB1,iA1,idmat)*Jint
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                          ! iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ELSE   !not ((startA.NE.startC).OR.(startB.NE.startD)) 
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET

                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat) + Input%DMAT_RHS(iC1,iD1,idmat)
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                               & +  dtemp_rhs*Jint
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                         !  iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET

                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat)
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                               & +  dtemp_rhs*Jint
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                         !  iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ENDIF
!                     iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
                  iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!      iOrbitalP = iOrbitalP + nOrbP
   ENDIF
ENDDO
   iOrbitalP = iOrbitalP + nOrbP
ENDDO

IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(1X,A,5I5)') 'ADDED Mat in Distribute_Coulomb1',Output%ndim(1),Output%ndim(2),&
&Output%ndim(3),Output%ndim(4),Output%ndim(5)
  DO i5=1,Output%ndim(5)
    DO i4=1,Output%ndim(4)
      DO i3=1,Output%ndim(3)
        DO i2=1,Output%ndim(2)
          WRITE(LUPRI,'(3X,4I5)') i2,i3,i4,i5
          WRITE(LUPRI,'(5X,5F12.8)') (Output%resultMat(i1,i2,i3,i4,i5)-TEMP(i1,i2),i1=1,Output%ndim(1))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE(TEMP)
ENDIF

END SUBROUTINE Distribute_Coulomb1

SUBROUTINE Distribute_Coulomb2(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
!sameLHS = TRUE sameRHS = TRUE sameOD = FALSE
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: STARTA1,STARTB1,STARTC1,STARTD1
Integer :: iC,iD,iOrbQ,idmat
real(realk) :: JINT,JFAC,dtemp_lhs,dtemp_rhs
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ

JFAC=INPUT%CoulombFactor
nAngmomP = PQ%P%p%nAngmom

iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
   startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
   startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
   IF(startA.NE.startB)THEN
      iOrbP=iOrbitalP-1+PASSPOFFSET
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=startB1+iAngB
               DO iAngA=1,nAngA
                  iA1=startA1+iAngA
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                  DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     IF(startC.NE.startD)THEN
                        iOrbQ=iOrbitalQ-1+PASSQOFFSET
                        DO iContD=1,nContD
                           STARTD1=startD-1+(iContD-1)*nAngD
                           DO iContC=1,nContC
                              STARTC1=startC-1+(iContC-1)*nAngC
                              DO iAngD=1,nAngD
                                 iD1=startD1+iAngD
                                 DO iAngC=1,nAngC
                                    iC1=startC1+iAngC
                                    iOrbQ=iOrbQ+1
                                    Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                    DO idmat=1,Input%NDMAT_RHS
                                       dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat) + Input%DMAT_RHS(iC1,iD1,idmat)
                                       Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                            &+ dtemp_rhs*Jint
                                       Output%ResultMat(iB1,iA1,1,1,idmat) = Output%ResultMat(iB1,iA1,1,1,idmat)&
                                            & + dtemp_rhs*Jint
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                        !
                        !iOrbitalQ = iOrbitalQ + nOrbQ
                     ELSE ! startC = startD
                        iOrbQ=iOrbitalQ-1+PASSQOFFSET
                        DO iContD=1,nContD
                           STARTD1=startD-1+(iContD-1)*nAngD
                           DO iContC=1,nContC
                              STARTC1=startC-1+(iContC-1)*nAngC
                              DO iAngD=1,nAngD
                                 iD1=startD1+iAngD
                                 DO iAngC=1,nAngC
                                    iC1=startC1+iAngC
                                    iOrbQ=iOrbQ+1
                                    Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                    DO idmat=1,Input%NDMAT_RHS
                                       dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat)
                                       Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                            & + dtemp_rhs*Jint
                                       Output%ResultMat(iB1,iA1,1,1,idmat) = Output%ResultMat(iB1,iA1,1,1,idmat)&
                                            & + dtemp_rhs*Jint
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                        !
                        !iOrbitalQ = iOrbitalQ + nOrbQ
                     ENDIF
                  ENDDO
                     iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !iOrbitalP = iOrbitalP + nOrbP
   ELSE
      iOrbP=iOrbitalP-1+PASSPOFFSET
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=startB1+iAngB
               DO iAngA=1,nAngA
                  iA1=startA1+iAngA
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                  DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     IF(startC.NE.startD)THEN
                        iOrbQ=iOrbitalQ-1+PASSQOFFSET
                        DO iContD=1,nContD
                           STARTD1=startD-1+(iContD-1)*nAngD
                           DO iContC=1,nContC
                              STARTC1=startC-1+(iContC-1)*nAngC
                              DO iAngD=1,nAngD
                                 iD1=startD1+iAngD
                                 DO iAngC=1,nAngC
                                    iC1=startC1+iAngC
                                    iOrbQ=iOrbQ+1
                                    Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                    DO idmat=1,Input%NDMAT_RHS
                                       dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat)+Input%DMAT_RHS(iC1,iD1,idmat)
                                       Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                            & + dtemp_rhs*Jint
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                        !iOrbitalQ = iOrbitalQ + nOrbQ
                     ELSE !startC=startD
                        iOrbQ=iOrbitalQ-1+PASSQOFFSET
                        DO iContD=1,nContD
                           STARTD1=startD-1+(iContD-1)*nAngD
                           DO iContC=1,nContC
                              STARTC1=startC-1+(iContC-1)*nAngC
                              DO iAngD=1,nAngD
                                 iD1=startD1+iAngD
                                 DO iAngC=1,nAngC
                                    iC1=startC1+iAngC
                                    iOrbQ=iOrbQ+1
                                    Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                    DO idmat=1,Input%NDMAT_RHS
                                       Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                            & + Input%DMAT_RHS(iD1,iC1,idmat)*Jint
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                        !iOrbitalQ = iOrbitalQ + nOrbQ
                     ENDIF
                  ENDDO
                  iOrbitalQ = iOrbitalQ + nOrbQ
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !iOrbitalP = iOrbitalP + nOrbP
   ENDIF
ENDDO
   iOrbitalP = iOrbitalP + nOrbP
ENDDO

END SUBROUTINE Distribute_Coulomb2

SUBROUTINE Distribute_Coulomb3(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
!sameLHS = TRUE sameRHS = FALSE sameOD = FALSE
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: STARTA1,STARTB1,STARTC1,STARTD1
Integer :: iC,iD,iOrbQ,idmat
real(realk) :: JINT,JFAC,dtemp_lhs,dtemp_rhs
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ

JFAC=INPUT%CoulombFactor

nAngmomP = PQ%P%p%nAngmom
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
   startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
   startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
   IF(startA.NE.startB)THEN
      iOrbP=iOrbitalP-1+PASSPOFFSET
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=startB1+iAngB
               DO iAngA=1,nAngA
                  iA1=startA1+iAngA
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                  DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO iContD=1,nContD
                        STARTD1=startD-1+(iContD-1)*nAngD
                        DO iContC=1,nContC
                           STARTC1=startC-1+(iContC-1)*nAngC
                           DO iAngD=1,nAngD
                              iD1=startD1+iAngD
                              DO iAngC=1,nAngC
                                 iC1=startC1+iAngC
                                 iOrbQ=iOrbQ+1
                                 Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                 DO idmat=1,Input%NDMAT_RHS
                                    dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat)
                                    Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                         & + dtemp_rhs*Jint
                                    Output%ResultMat(iB1,iA1,1,1,idmat) = Output%ResultMat(iB1,iA1,1,1,idmat)&
                                         & + dtemp_rhs*Jint
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                     !
                  ENDDO

                     iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !iOrbitalP = iOrbitalP + nOrbP
   ELSE
      iOrbP=iOrbitalP-1+PASSPOFFSET
      iContP=0
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=startB1+iAngB
               DO iAngA=1,nAngA
                  iA1=startA1+iAngA
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                  DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO iContD=1,nContD
                        STARTD1=startD-1+(iContD-1)*nAngD
                        DO iContC=1,nContC
                           STARTC1=startC-1+(iContC-1)*nAngC
                           DO iAngD=1,nAngD
                              iD1=startD1+iAngD
                              DO iAngC=1,nAngC
                                 iC1=startC1+iAngC
                                 iOrbQ=iOrbQ+1
                                 Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                 DO idmat=1,Input%NDMAT_RHS
                                    Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                         & + Input%DMAT_RHS(iD1,iC1,idmat)*Jint
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                     iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !iOrbitalP = iOrbitalP + nOrbP
   ENDIF
ENDDO
iOrbitalP = iOrbitalP + nOrbP
ENDDO

END SUBROUTINE Distribute_Coulomb3

SUBROUTINE Distribute_Coulomb4(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
!sameLHS = FALSE sameRHS = TRUE sameOD = FALSE
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: STARTA1,STARTB1,STARTC1,STARTD1
Integer :: iC,iD,iOrbQ,idmat
real(realk) :: JINT,JFAC,dtemp_lhs,dtemp_rhs
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ

JFAC=INPUT%CoulombFactor

nAngmomP = PQ%P%p%nAngmom
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
   startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
   startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
   iOrbP=iOrbitalP-1+PASSPOFFSET
   DO iContB=1,nContB
      STARTB1=startB-1+(iContB-1)*nAngB
      DO iContA=1,nContA
         STARTA1=startA-1+(iContA-1)*nAngA
         DO iAngB=1,nAngB
            iB1=startB1+iAngB
            DO iAngA=1,nAngA
               iA1=startA1+iAngA
               iOrbP=iOrbP+1
               iOrbitalQ = 1
               endAngmomQ = PQ%Q%p%nAngmom
               IF (PQ%samePQ) endAngmomQ = iAngmomP
               DO iAngmomQ=1,endAngmomQ
               DO iPassQ = 1, PQ%Q%p%nPasses
                  nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                  iC = PQ%Q%p%indexAng1(iAngmomQ)
                  iD = PQ%Q%p%indexAng2(iAngmomQ)
                  nContC = PQ%Q%p%orbital1%nContracted(iC)
                  nContD = PQ%Q%p%orbital2%nContracted(iD)
                  nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                  nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                  PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                  startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                  startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                  IF(startC.NE.startD)THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO iContD=1,nContD
                        STARTD1=startD-1+(iContD-1)*nAngD
                        DO iContC=1,nContC
                           STARTC1=startC-1+(iContC-1)*nAngC
                           DO iAngD=1,nAngD
                              iD1=startD1+iAngD
                              DO iAngC=1,nAngC
                                 iC1=startC1+iAngC
                                 iOrbQ=iOrbQ+1
                                 Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                 DO idmat=1,Input%NDMAT_RHS
                                    dtemp_rhs = Input%DMAT_RHS(iD1,iC1,idmat) + Input%DMAT_RHS(iC1,iD1,idmat)
                                    Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                         & + dtemp_rhs*Jint
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                     !
                     !iOrbitalQ = iOrbitalQ + nOrbQ
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO iContD=1,nContD
                        STARTD1=startD-1+(iContD-1)*nAngD
                        DO iContC=1,nContC
                           STARTC1=startC-1+(iContC-1)*nAngC
                           DO iAngD=1,nAngD
                              iD1=startD1+iAngD
                              DO iAngC=1,nAngC
                                 iC1=startC1+iAngC
                                 iOrbQ=iOrbQ+1
                                 Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                 DO idmat=1,Input%NDMAT_RHS
                                    Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                         & + Input%DMAT_RHS(iD1,iC1,idmat)*Jint
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                     !
                     !iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDIF
               ENDDO
                  iOrbitalQ = iOrbitalQ + nOrbQ
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
   iOrbitalP = iOrbitalP + nOrbP
ENDDO

END SUBROUTINE Distribute_Coulomb4

SUBROUTINE Distribute_Coulomb5(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
!sameLHS = FALSE sameRHS = FALSE sameOD = TRUE
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: STARTA1,STARTB1,STARTC1,STARTD1
Integer :: iC,iD,iOrbQ,idmat
real(realk) :: JINT,JFAC,dtemp_lhs,dtemp_rhs
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ

JFAC=INPUT%CoulombFactor

nAngmomP = PQ%P%p%nAngmom
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
   startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
   startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
   iOrbP=iOrbitalP-1+PASSPOFFSET
   DO iContB=1,nContB
      STARTB1=startB-1+(iContB-1)*nAngB
      DO iContA=1,nContA
         STARTA1=startA-1+(iContA-1)*nAngA
         DO iAngB=1,nAngB
            iB1=startB1+iAngB
            DO iAngA=1,nAngA
               iA1=startA1+iAngA
               iOrbP=iOrbP+1
               iOrbitalQ = 1
               endAngmomQ = PQ%Q%p%nAngmom
               IF (PQ%samePQ) endAngmomQ = iAngmomP
               DO iAngmomQ=1,endAngmomQ
               DO iPassQ = 1, PQ%Q%p%nPasses
                  nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                  iC = PQ%Q%p%indexAng1(iAngmomQ)
                  iD = PQ%Q%p%indexAng2(iAngmomQ)
                  nContC = PQ%Q%p%orbital1%nContracted(iC)
                  nContD = PQ%Q%p%orbital2%nContracted(iD)
                  nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                  nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                  PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                  startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                  startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                  IF ((startA.NE.startC).OR.(startB.NE.startD))THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO iContD=1,nContD
                        STARTD1=startD-1+(iContD-1)*nAngD
                        DO iContC=1,nContC
                           STARTC1=startC-1+(iContC-1)*nAngC
                           DO iAngD=1,nAngD
                              iD1=startD1+iAngD
                              DO iAngC=1,nAngC
                                 iC1=startC1+iAngC
                                 iOrbQ=iOrbQ+1
                                 Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                 DO idmat=1,Input%NDMAT_RHS
                                    Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                         & + Input%DMAT_RHS(iD1,iC1,idmat)*Jint
                                    Output%ResultMat(iC1,iD1,1,1,idmat) = Output%ResultMat(iC1,iD1,1,1,idmat)&
                                         & + Input%DMAT_RHS(iB1,iA1,idmat)*Jint
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                     !
                     !iOrbitalQ = iOrbitalQ + nOrbQ
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO iContD=1,nContD
                        STARTD1=startD-1+(iContD-1)*nAngD
                        DO iContC=1,nContC
                           STARTC1=startC-1+(iContC-1)*nAngC
                           DO iAngD=1,nAngD
                              iD1=startD1+iAngD
                              DO iAngC=1,nAngC
                                 iC1=startC1+iAngC
                                 iOrbQ=iOrbQ+1
                                 Jint = JFAC*QPmat(iOrbQ,iOrbP)
                                 DO idmat=1,Input%NDMAT_RHS
                                    Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                         & + Input%DMAT_RHS(iD1,iC1,idmat)*Jint
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                     !
                     !iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDIF
               ENDDO
                  iOrbitalQ = iOrbitalQ + nOrbQ
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
   iOrbitalP = iOrbitalP + nOrbP
ENDDO

END SUBROUTINE Distribute_Coulomb5

SUBROUTINE Distribute_Coulomb6(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
!sameLHS = FALSE sameRHS = FALSE sameOD = FALSE
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: STARTA1,STARTB1,STARTC1,STARTD1
Integer :: iC,iD,iOrbQ,idmat
real(realk) :: JINT,JFAC,dtemp_lhs,dtemp_rhs
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ

JFAC=INPUT%CoulombFactor

nAngmomP = PQ%P%p%nAngmom
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
   startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
   startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
   iOrbP=iOrbitalP-1+PASSPOFFSET
   DO iContB=1,nContB
      STARTB1=startB-1+(iContB-1)*nAngB
      DO iContA=1,nContA
         STARTA1=startA-1+(iContA-1)*nAngA
         DO iAngB=1,nAngB
            iB1=startB1+iAngB
            DO iAngA=1,nAngA
               iA1=startA1+iAngA
               iOrbP=iOrbP+1
               iOrbitalQ = 1
               endAngmomQ = PQ%Q%p%nAngmom
               IF (PQ%samePQ) endAngmomQ = iAngmomP
               DO iAngmomQ=1,endAngmomQ
               DO iPassQ = 1, PQ%Q%p%nPasses
                  nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                  iC = PQ%Q%p%indexAng1(iAngmomQ)
                  iD = PQ%Q%p%indexAng2(iAngmomQ)
                  nContC = PQ%Q%p%orbital1%nContracted(iC)
                  nContD = PQ%Q%p%orbital2%nContracted(iD)
                  nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                  nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                  PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                  startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                  startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                  iOrbQ=iOrbitalQ-1+PASSQOFFSET
                  DO iContD=1,nContD
                     STARTD1=startD-1+(iContD-1)*nAngD
                     DO iContC=1,nContC
                        STARTC1=startC-1+(iContC-1)*nAngC
                        DO iAngD=1,nAngD
                           iD1=startD1+iAngD
                           DO iAngC=1,nAngC
                              iC1=startC1+iAngC
                              iOrbQ=iOrbQ+1
                              Jint = JFAC*QPmat(iOrbQ,iOrbP)
                              DO idmat=1,Input%NDMAT_RHS
                                 Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                      & + Input%DMAT_RHS(iD1,iC1,idmat)*Jint
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  !
               ENDDO
                  iOrbitalQ = iOrbitalQ + nOrbQ
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
   iOrbitalP = iOrbitalP + nOrbP
ENDDO

END SUBROUTINE Distribute_Coulomb6

SUBROUTINE Distribute_Exchange(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: iC,iD,iOrbQ,idmat
Integer :: startA1,startB1,startC1,startD1
real(realk) :: KINT,KFAC
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ

!
KFAC=INPUT%exchangeFactor
nAngmomP = PQ%P%p%nAngmom
!nderivP=Input%nDerivP
!totOrbP = PQ%P%p%totOrbitals
!
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
   startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
   startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
   IF(startA.NE.startB)THEN
      iOrbP=iOrbitalP-1+PASSPOFFSET
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=iAngB+STARTB1
               DO iAngA=1,nAngA
                  iA1=iAngA+STARTA1
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                  DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=STARTD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=STARTC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS
                                          Output%ResultMat(iA1,iC1,1,1,idmat) = Output%ResultMat(iA1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iD1,idmat)
                                          Output%ResultMat(iB1,iC1,1,1,idmat) = Output%ResultMat(iB1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iA1,iD1,idmat)
                                          Output%ResultMat(iA1,iD1,1,1,idmat) = Output%ResultMat(iA1,iD1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iC1,idmat)
                                          Output%ResultMat(iB1,iD1,1,1,idmat) = Output%ResultMat(iB1,iD1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iA1,iC1,idmat)
                                          Output%ResultMat(iC1,iA1,1,1,idmat) = Output%ResultMat(iC1,iA1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iD1,iB1,idmat)
                                          Output%ResultMat(iD1,iA1,1,1,idmat) = Output%ResultMat(iD1,iA1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iC1,iB1,idmat)
                                          Output%ResultMat(iD1,iB1,1,1,idmat) = Output%ResultMat(iD1,iB1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iC1,iA1,idmat)
                                          Output%ResultMat(iC1,iB1,1,1,idmat) = Output%ResultMat(iC1,iB1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iD1,iA1,idmat)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          Output%ResultMat(iA1,iC1,1,1,idmat) = Output%ResultMat(iA1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iD1,idmat)
                                          Output%ResultMat(iB1,iC1,1,1,idmat) = Output%ResultMat(iB1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iA1,iD1,idmat)
                                          Output%ResultMat(iC1,iA1,1,1,idmat) = Output%ResultMat(iC1,iA1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iD1,iB1,idmat)
                                          Output%ResultMat(iC1,iB1,1,1,idmat) = Output%ResultMat(iC1,iB1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iD1,iA1,idmat)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS
                                          Output%ResultMat(iA1,iC1,1,1,idmat) = Output%ResultMat(iA1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iD1,idmat)
                                          Output%ResultMat(iB1,iC1,1,1,idmat) = Output%ResultMat(iB1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iA1,iD1,idmat)
                                          Output%ResultMat(iA1,iD1,1,1,idmat) = Output%ResultMat(iA1,iD1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iC1,idmat)
                                          Output%ResultMat(iB1,iD1,1,1,idmat) = Output%ResultMat(iB1,iD1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iA1,iC1,idmat)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          Output%ResultMat(iA1,iC1,1,1,idmat) = Output%ResultMat(iA1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iD1,idmat)
                                          Output%ResultMat(iB1,iC1,1,1,idmat) = Output%ResultMat(iB1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iA1,iD1,idmat)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ENDIF
                  ENDDO
                  iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!      iOrbitalP = iOrbitalP + nOrbP
   ELSE !NOT (startA NE startB)
      iOrbP=iOrbitalP-1+PASSPOFFSET
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=startB1+iAngB
               DO iAngA=1,nAngA
                  iA1=startA1+iAngA
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                  DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          Output%ResultMat(iA1,iC1,1,1,idmat) = Output%ResultMat(iA1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iD1,idmat)
                                          Output%ResultMat(iA1,iD1,1,1,idmat) = Output%ResultMat(iA1,iD1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iC1,idmat)
                                          Output%ResultMat(iC1,iA1,1,1,idmat) = Output%ResultMat(iC1,iA1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iD1,iB1,idmat)
                                          Output%ResultMat(iD1,iA1,1,1,idmat) = Output%ResultMat(iD1,iA1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iC1,iB1,idmat)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          Output%ResultMat(iA1,iC1,1,1,idmat) = Output%ResultMat(iA1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iD1,idmat)
                                          Output%ResultMat(iC1,iA1,1,1,idmat) = Output%ResultMat(iC1,iA1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iD1,iB1,idmat)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ELSE   !not ((startA.NE.startC).OR.(startB.NE.startD)) 
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          Output%ResultMat(iA1,iC1,1,1,idmat) = Output%ResultMat(iA1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iD1,idmat)
                                          Output%ResultMat(iA1,iD1,1,1,idmat) = Output%ResultMat(iA1,iD1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iC1,idmat)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          Output%ResultMat(iA1,iC1,1,1,idmat) = Output%ResultMat(iA1,iC1,1,1,idmat)&
                                               & -  kint*Input%DMAT_RHS(iB1,iD1,idmat)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ENDIF
                  ENDDO
                     iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!      iOrbitalP = iOrbitalP + nOrbP
   ENDIF
   iOrbitalP = iOrbitalP + nOrbP
ENDDO
ENDDO
IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(1X,A,5I5)') 'ResultMat in Distribute_Exchange',Output%ndim(1),Output%ndim(2),&
&Output%ndim(3),Output%ndim(4),Output%ndim(5)
  DO i5=1,Output%ndim(5)
    DO i4=1,Output%ndim(4)
      DO i3=1,Output%ndim(3)
        DO i2=1,Output%ndim(2)
          WRITE(LUPRI,'(3X,4I5)') i2,i3,i4,i5
          WRITE(LUPRI,'(5X,5F12.8)') (Output%resultMat(i1,i2,i3,i4,i5),i1=1,Output%ndim(1))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE Distribute_Exchange

SUBROUTINE Distribute_Exchange2(Integral,PQ,QPmat,dimQ,dimP,Input,ODmat,ODRESmat,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
TYPE(ODmatrixitem),intent(in)   :: ODmat
TYPE(ODmatrixitem)   :: ODRESmat(:)
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: iC,iD,iOrbQ,idmat
Integer :: startA1,startB1,startC1,startD1
real(realk) :: KINT,KFAC
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ
Integer :: BD_IOD, AD_IOD, BC_IOD,AC_IOD, DB_IOD, CB_IOD, CA_IOD, DA_IOD
Integer :: BD_FIOD(Input%NDMAT_RHS), AD_FIOD(Input%NDMAT_RHS), BC_FIOD(Input%NDMAT_RHS)
Integer :: AC_FIOD(Input%NDMAT_RHS), DB_FIOD(Input%NDMAT_RHS), CB_FIOD(Input%NDMAT_RHS)
Integer :: CA_FIOD(Input%NDMAT_RHS), DA_FIOD(Input%NDMAT_RHS)
Integer :: iAC,iAD,iBC,iBD,iCA,iCB,iDA,iDB

!
KFAC=INPUT%exchangeFactor
nAngmomP = PQ%P%p%nAngmom
!
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
   startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
   startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
   IF(startA.NE.startB)THEN
      iOrbP=iOrbitalP-1+PASSPOFFSET
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=iAngB+STARTB1
               DO iAngA=1,nAngA
                  iA1=iAngA+STARTA1
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                  DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     BD_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startD))
                     AD_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startD))
                     BC_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startC))
                     AC_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startC))
                     DB_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startB))
                     CB_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startB))
                     CA_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startA))
                     DA_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startA))
                     DO idmat=1,Input%NDMAT_RHS
                        BD_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startD))
                        AD_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startD))
                        BC_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startC))
                        AC_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startC))
                        DB_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startB))
                        CB_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startB))
                        CA_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startA))
                        DA_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startA))
                     ENDDO
                     IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=STARTD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=STARTC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
                                       iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
                                       iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA
                                       iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB
                                       iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC
                                       iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD
                                       iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
                                       iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC
                                       DO idmat=1,Input%NDMAT_RHS
ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)= ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)-  kint*&!D(iB1,iD1)
     &ODmat%BATCH(BD_IOD)%matAB%elms(iBD)
ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)= ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)-  kint*&!D(iA1,iD1)
     &ODmat%BATCH(AD_IOD)%matAB%elms(iAD)
ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)= ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)-  kint*&!D(iB1,iC1)
     &ODmat%BATCH(BC_IOD)%matAB%elms(iBC)
ODRESmat(idmat)%BATCH(BD_FIOD(idmat))%matAB%elms(iBD)= ODRESmat(idmat)%BATCH(BD_FIOD(idmat))%matAB%elms(iBD)-  kint*&!D(iA1,iC1)
     &ODmat%BATCH(AC_IOD)%matAB%elms(iAC)
ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)= ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)-  kint*&!D(iD1,iB1)
     &ODmat%BATCH(DB_IOD)%matAB%elms(iDB)
ODRESmat(idmat)%BATCH(DA_FIOD(idmat))%matAB%elms(iDA)= ODRESmat(idmat)%BATCH(DA_FIOD(idmat))%matAB%elms(iDA)-  kint*&!D(iC1,iB1)
     &ODmat%BATCH(CB_IOD)%matAB%elms(iCB)
ODRESmat(idmat)%BATCH(DB_FIOD(idmat))%matAB%elms(iDB)= ODRESmat(idmat)%BATCH(DB_FIOD(idmat))%matAB%elms(iDB)-  kint*&!D(iC1,iA1)
     &ODmat%BATCH(CA_IOD)%matAB%elms(iCA)
ODRESmat(idmat)%BATCH(CB_FIOD(idmat))%matAB%elms(iCB)= ODRESmat(idmat)%BATCH(CB_FIOD(idmat))%matAB%elms(iCB)-  kint*&!D(iD1,iA1)
     &ODmat%BATCH(DA_IOD)%matAB%elms(iDA)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
                                       iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
                                       iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA
                                       iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB
                                       iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC
                                       iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD
                                       iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
                                       iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC
                                       DO idmat=1,Input%NDMAT_RHS 
ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)= ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)-  kint*&!D(iB1,iD1)
     &ODmat%BATCH(BD_IOD)%matAB%elms(iBD)
ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)= ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)-  kint*&!D(iA1,iD1)
     &ODmat%BATCH(AD_IOD)%matAB%elms(iAD)
ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)= ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)-  kint*&!D(iD1,iB1)
     &ODmat%BATCH(DB_IOD)%matAB%elms(iDB)
ODRESmat(idmat)%BATCH(CB_FIOD(idmat))%matAB%elms(iCB)= ODRESmat(idmat)%BATCH(CB_FIOD(idmat))%matAB%elms(iCB)-  kint*&!D(iD1,iA1)
     &ODmat%BATCH(DA_IOD)%matAB%elms(iDA)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
                                       iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
                                       iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA
                                       iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB
                                       iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC
                                       iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD
                                       iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
                                       iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC
                                       DO idmat=1,Input%NDMAT_RHS
ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)= ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)-  kint*&!D(iB1,iD1)
     &ODmat%BATCH(BD_IOD)%matAB%elms(iBD)
ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)= ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)-  kint*&!D(iA1,iD1)
     &ODmat%BATCH(AD_IOD)%matAB%elms(iAD)
ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)= ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)-  kint*&!D(iB1,iC1)
     &ODmat%BATCH(BC_IOD)%matAB%elms(iBC)
ODRESmat(idmat)%BATCH(BD_FIOD(idmat))%matAB%elms(iBD)= ODRESmat(idmat)%BATCH(BD_FIOD(idmat))%matAB%elms(iBD)-  kint*&!D(iA1,iC1)
     &ODmat%BATCH(AC_IOD)%matAB%elms(iAC)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
                                       iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
                                       iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA
                                       iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB
                                       iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC
                                       iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD
                                       iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
                                       iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC
                                       DO idmat=1,Input%NDMAT_RHS 
ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)= ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)-  kint*&!D(iB1,iD1)
     &ODmat%BATCH(BD_IOD)%matAB%elms(iBD)
ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)= ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)-  kint*&!D(iA1,iD1)
     &ODmat%BATCH(AD_IOD)%matAB%elms(iAD)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ENDIF
                  ENDDO
                  iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!      iOrbitalP = iOrbitalP + nOrbP
   ELSE !NOT (startA NE startB)
      iOrbP=iOrbitalP-1+PASSPOFFSET
      DO iContB=1,nContB
         STARTB1=startB-1+(iContB-1)*nAngB
         DO iContA=1,nContA
            STARTA1=startA-1+(iContA-1)*nAngA
            DO iAngB=1,nAngB
               iB1=startB1+iAngB
               DO iAngA=1,nAngA
                  iA1=startA1+iAngA
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  DO iAngmomQ=1,endAngmomQ
                  DO iPassQ = 1, PQ%Q%p%nPasses
                     nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                     iC = PQ%Q%p%indexAng1(iAngmomQ)
                     iD = PQ%Q%p%indexAng2(iAngmomQ)
                     nContC = PQ%Q%p%orbital1%nContracted(iC)
                     nContD = PQ%Q%p%orbital2%nContracted(iD)
                     nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                     nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                     PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                     startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                     startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)


                     BD_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startD))
!                     AD_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startD))
                     BC_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startC))
!                     AC_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startC))
                     DB_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startB))
                     CB_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startB))
!                     CA_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startA))
!                     DA_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startA))
                     DO idmat=1,Input%NDMAT_RHS
!                        BD_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startD))
                        AD_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startD))
!                        BC_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startC))
                        AC_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startC))
!                        DB_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startB))
!                        CB_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startB))
                        CA_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startA))
                        DA_FIOD(idmat) = ODRESmat(idmat)%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startA))
                     ENDDO
                     IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
                                       iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
                                       iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA
                                       iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB
                                       iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC
                                       iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD
                                       iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
                                       iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC
                                       DO idmat=1,Input%NDMAT_RHS 
ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)= ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)-  kint*&!D(iB1,iD1)
     &ODmat%BATCH(BD_IOD)%matAB%elms(iBD)
ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)= ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)-  kint*&!D(iB1,iC1)
     &ODmat%BATCH(BC_IOD)%matAB%elms(iBC)
ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)= ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)-  kint*&!D(iD1,iB1)
     &ODmat%BATCH(DB_IOD)%matAB%elms(iDB)
ODRESmat(idmat)%BATCH(DA_FIOD(idmat))%matAB%elms(iDA)= ODRESmat(idmat)%BATCH(DA_FIOD(idmat))%matAB%elms(iDA)-  kint*&!D(iC1,iB1)
     &ODmat%BATCH(CB_IOD)%matAB%elms(iCB)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
                                       iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
                                       iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA
                                       iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB
                                       iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC
                                       iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD
                                       iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
                                       iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC
                                       DO idmat=1,Input%NDMAT_RHS 
ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)= ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)-  kint*&!D(iB1,iD1)
     &ODmat%BATCH(BD_IOD)%matAB%elms(iBD)
ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)= ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)-  kint*&!D(iD1,iB1)
     &ODmat%BATCH(DB_IOD)%matAB%elms(iDB)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ELSE   !not ((startA.NE.startC).OR.(startB.NE.startD)) 
                        IF (startC.NE.startD)THEN
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
                                       iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
                                       iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA
                                       iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB
                                       iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC
                                       iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD
                                       iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
                                       iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC
                                       DO idmat=1,Input%NDMAT_RHS 
ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)= ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)-  kint*&!D(iB1,iD1)
     &ODmat%BATCH(BD_IOD)%matAB%elms(iBD)
ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)= ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)-  kint*&!D(iB1,iC1)
     &ODmat%BATCH(BC_IOD)%matAB%elms(iBC)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ELSE
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              STARTD1=startD-1+(iContD-1)*nAngD
                              DO iContC=1,nContC
                                 STARTC1=startC-1+(iContC-1)*nAngC
                                 DO iAngD=1,nAngD
                                    iD1=startD1+iAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC1+iAngC
                                       iOrbQ=iOrbQ+1
                                       kint = KFAC*QPmat(iOrbQ,iOrbP)
                                       iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
                                       iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
                                       iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA
                                       iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB
                                       iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC
                                       iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD
                                       iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
                                       iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC
                                       DO idmat=1,Input%NDMAT_RHS 
ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)= ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms(iAC)-  kint*&!D(iB1,iD1)
     &ODmat%BATCH(BD_IOD)%matAB%elms(iBD)
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDIF
                     ENDIF
                  ENDDO
                     iOrbitalQ = iOrbitalQ + nOrbQ
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!      iOrbitalP = iOrbitalP + nOrbP
   ENDIF
   iOrbitalP = iOrbitalP + nOrbP
ENDDO
ENDDO
IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(1X,A,5I5)') 'ResultMat in Distribute_Exchange2',Output%ndim(1),Output%ndim(2),&
&Output%ndim(3),Output%ndim(4),Output%ndim(5)
  DO i5=1,Output%ndim(5)
    DO i4=1,Output%ndim(4)
      DO i3=1,Output%ndim(3)
        DO i2=1,Output%ndim(2)
          WRITE(LUPRI,'(3X,4I5)') i2,i3,i4,i5
          WRITE(LUPRI,'(5X,5F12.8)') (Output%resultMat(i1,i2,i3,i4,i5),i1=1,Output%ndim(1))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE Distribute_Exchange2

SUBROUTINE Distribute_Exchange3(Integral,PQ,QPmat,dimQ,dimP,Input,ODmat,ODRESmat,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
TYPE(ODmatrixitem),intent(in)   :: ODmat
TYPE(ODmatrixitem)   :: ODRESmat(:)
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: iC,iD,iOrbQ,idmat
Integer :: startA1,startB1,startC1,startD1
real(realk) :: KINT,KFAC
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ
Integer :: BD_IOD, AD_IOD, BC_IOD,AC_IOD, DB_IOD, CB_IOD, CA_IOD, DA_IOD
Integer :: BD_FIOD(Input%NDMAT_RHS), AD_FIOD(Input%NDMAT_RHS), BC_FIOD(Input%NDMAT_RHS)
Integer :: AC_FIOD(Input%NDMAT_RHS), DB_FIOD(Input%NDMAT_RHS), CB_FIOD(Input%NDMAT_RHS)
Integer :: CA_FIOD(Input%NDMAT_RHS), DA_FIOD(Input%NDMAT_RHS)
Integer :: iAC,iAD,iBC,iBD,iCA,iCB,iDA,iDB

!
KFAC=INPUT%exchangeFactor
nAngmomP = PQ%P%p%nAngmom
!
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
      nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
      iA = PQ%P%p%indexAng1(iAngmomP)
      iB = PQ%P%p%indexAng2(iAngmomP)
      nContA = PQ%P%p%orbital1%nContracted(iA)
      nContB = PQ%P%p%orbital2%nContracted(iB)
      nAngA = PQ%P%p%orbital1%nOrbComp(iA)
      nAngB = PQ%P%p%orbital2%nOrbComp(iB)
      PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
      startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
      startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
      IF(startA.NE.startB)THEN
         iOrbP=iOrbitalP-1+PASSPOFFSET
         iOrbitalQ = 1
         endAngmomQ = PQ%Q%p%nAngmom
         IF (PQ%samePQ) endAngmomQ = iAngmomP
         DO iAngmomQ=1,endAngmomQ
            DO iPassQ = 1, PQ%Q%p%nPasses
               nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
               iC = PQ%Q%p%indexAng1(iAngmomQ)
               iD = PQ%Q%p%indexAng2(iAngmomQ)
               nContC = PQ%Q%p%orbital1%nContracted(iC)
               nContD = PQ%Q%p%orbital2%nContracted(iD)
               nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
               nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
               PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
               startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
               startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
               BD_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startD))
               AD_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startD))
               BC_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startC))
               AC_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startC))
               DB_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startB))
               CB_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startB))
               CA_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startA))
               DA_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startA))
               IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                  IF (startC.NE.startD)THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_FIOD(idmat))%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,KFAC)
!                        CALL EXCHANGE_3BCAD(ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,KFAC)
!                        CALL EXCHANGE_3ADBC(ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,KFAC)
!                        CALL EXCHANGE_3BDAC(ODRESmat(idmat)%BATCH(BD_FIOD(idmat))%matAB%elms,ODmat%BATCH(AC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,KFAC)
!                        CALL EXCHANGE_3CADB(ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,KFAC)
!                        CALL EXCHANGE_3DACB(ODRESmat(idmat)%BATCH(DA_FIOD(idmat))%matAB%elms,ODmat%BATCH(CB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,KFAC)
!                        CALL EXCHANGE_3DBCA(ODRESmat(idmat)%BATCH(DB_IOD)%matAB%elms,ODmat%BATCH(CA_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,KFAC)
!                        CALL EXCHANGE_3CBDA(ODRESmat(idmat)%BATCH(CB_IOD)%matAB%elms,ODmat%BATCH(DA_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,KFAC)
                        CALL EXCHANGE_CASE1(&
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,& !ACBD
                             ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,& !BCAD
                             ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,& !ADBC
                             ODRESmat(idmat)%BATCH(BD_IOD)%matAB%elms,ODmat%BATCH(AC_IOD)%matAB%elms,& !BDAC
                             ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,& !CADB
                             ODRESmat(idmat)%BATCH(DA_IOD)%matAB%elms,ODmat%BATCH(CB_IOD)%matAB%elms,& !DACB
                             ODRESmat(idmat)%BATCH(DB_IOD)%matAB%elms,ODmat%BATCH(CA_IOD)%matAB%elms,& !DBCA
                             ODRESmat(idmat)%BATCH(CB_IOD)%matAB%elms,ODmat%BATCH(DA_IOD)%matAB%elms,& !CBDA
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                     ENDDO
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3BCAD(ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3CADB(ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3CBDA(ODRESmat(idmat)%BATCH(CB_IOD)%matAB%elms,ODmat%BATCH(DA_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE2(&   ! 1 2 5 8
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(CB_IOD)%matAB%elms,ODmat%BATCH(DA_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)

                     ENDDO
                  ENDIF
               ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                  IF (startC.NE.startD)THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3BCAD(ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3ADBC(ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3BDAC(ODRESmat(idmat)%BATCH(BD_IOD)%matAB%elms,ODmat%BATCH(AC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE3(&   ! 1 2 3 4
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(BD_IOD)%matAB%elms,ODmat%BATCH(AC_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,1.D0)

                     ENDDO
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3BCAD(ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE4(&   ! 1 2
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,1.D0)
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            iOrbitalQ = iOrbitalQ + nOrbQ
         ENDDO
      ELSE !NOT (startA NE startB)
         iOrbP=iOrbitalP-1+PASSPOFFSET
         iOrbitalQ = 1
         endAngmomQ = PQ%Q%p%nAngmom
         IF (PQ%samePQ) endAngmomQ = iAngmomP
         DO iAngmomQ=1,endAngmomQ
            DO iPassQ = 1, PQ%Q%p%nPasses
               nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
               iC = PQ%Q%p%indexAng1(iAngmomQ)
               iD = PQ%Q%p%indexAng2(iAngmomQ)
               nContC = PQ%Q%p%orbital1%nContracted(iC)
               nContD = PQ%Q%p%orbital2%nContracted(iD)
               nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
               nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
               PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
               startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
               startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                     BD_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startD))
                     BC_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startC))
                     DB_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startB))
                     CB_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startB))
                     AD_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startD))
                     AC_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startC))
                     CA_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startA))
                     DA_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startA))
               IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                  IF (startC.NE.startD)THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3ADBC(ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3CADB(ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3DACB(ODRESmat(idmat)%BATCH(DA_IOD)%matAB%elms,ODmat%BATCH(CB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE5(&   ! 1 3 5 6
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(DA_IOD)%matAB%elms,ODmat%BATCH(CB_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,1.D0)
                     ENDDO
                     
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3CADB(ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE6(&   ! 1 5
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,1.D0)

                     ENDDO
                  ENDIF
               ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                  IF (startC.NE.startD)THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3ADBC(ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE7(&   ! 1 3
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,1.D0)

                     ENDDO
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE8(&   ! 1
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,1.D0)

                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            iOrbitalQ = iOrbitalQ + nOrbQ
         ENDDO
      ENDIF
   ENDDO
      iOrbitalP = iOrbitalP + nOrbP
ENDDO

IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(1X,A,5I5)') 'ResultMat in Distribute_Exchange2',Output%ndim(1),Output%ndim(2),&
&Output%ndim(3),Output%ndim(4),Output%ndim(5)
  DO i5=1,Output%ndim(5)
    DO i4=1,Output%ndim(4)
      DO i3=1,Output%ndim(3)
        DO i2=1,Output%ndim(2)
          WRITE(LUPRI,'(3X,4I5)') i2,i3,i4,i5
          WRITE(LUPRI,'(5X,5F12.8)') (Output%resultMat(i1,i2,i3,i4,i5),i1=1,Output%ndim(1))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE Distribute_Exchange3

SUBROUTINE Distribute_Exchange3SYM(Integral,PQ,QPmat,dimQ,dimP,Input,ODmat,ODRESmat,Output,LUPRI,IPRINT)
! We assume that the density matrix and thereby the Fock matrix is symmetrix
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
TYPE(ODmatrixitem),intent(in)   :: ODmat
TYPE(ODmatrixitem)   :: ODRESmat(:)
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: iC,iD,iOrbQ,idmat
Integer :: startA1,startB1,startC1,startD1
real(realk) :: KINT,KFAC
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ
Integer :: BD_IOD, AD_IOD, BC_IOD,AC_IOD, DB_IOD, CB_IOD, CA_IOD, DA_IOD
Integer :: BD_FIOD(Input%NDMAT_RHS), AD_FIOD(Input%NDMAT_RHS), BC_FIOD(Input%NDMAT_RHS)
Integer :: AC_FIOD(Input%NDMAT_RHS), DB_FIOD(Input%NDMAT_RHS), CB_FIOD(Input%NDMAT_RHS)
Integer :: CA_FIOD(Input%NDMAT_RHS), DA_FIOD(Input%NDMAT_RHS)
Integer :: iAC,iAD,iBC,iBD,iCA,iCB,iDA,iDB
LOGICAL :: DoDistibute(8)
!
DoDistibute = .FALSE.
KFAC=INPUT%exchangeFactor
nAngmomP = PQ%P%p%nAngmom
!
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
      nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
      iA = PQ%P%p%indexAng1(iAngmomP)
      iB = PQ%P%p%indexAng2(iAngmomP)
      nContA = PQ%P%p%orbital1%nContracted(iA)
      nContB = PQ%P%p%orbital2%nContracted(iB)
      nAngA = PQ%P%p%orbital1%nOrbComp(iA)
      nAngB = PQ%P%p%orbital2%nOrbComp(iB)
      PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
      startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
      startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
      IF(startA.NE.startB)THEN
         iOrbP=iOrbitalP-1+PASSPOFFSET
         iOrbitalQ = 1
         endAngmomQ = PQ%Q%p%nAngmom
         IF (PQ%samePQ) endAngmomQ = iAngmomP
         DO iAngmomQ=1,endAngmomQ
            DO iPassQ = 1, PQ%Q%p%nPasses
               nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
               iC = PQ%Q%p%indexAng1(iAngmomQ)
               iD = PQ%Q%p%indexAng2(iAngmomQ)
               nContC = PQ%Q%p%orbital1%nContracted(iC)
               nContD = PQ%Q%p%orbital2%nContracted(iD)
               nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
               nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
               PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
               startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
               startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
               BD_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startD))
               AD_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startD))
               BC_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startC))
               AC_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startC))
               DB_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startB))
               CB_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startB))
               CA_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startA))
               DA_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startA))
               IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                  IF (startC.NE.startD)THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                        DO idmat=1,Input%NDMAT_RHS
                           CALL EXCHANGE_CASE1SYM(&
                                ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,& !ACBD      A<C
                                ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,& !BCAD      B<C 
                                ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,& !ADBC      A<D 
                                ODRESmat(idmat)%BATCH(BD_IOD)%matAB%elms,ODmat%BATCH(AC_IOD)%matAB%elms,& !BDAC      B<D 
                                !ODRESmat(idmat)%BATCH(CA_FIOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,& !CADB      C<A 
                                !ODRESmat(idmat)%BATCH(DA_FIOD)%matAB%elms,ODmat%BATCH(CB_IOD)%matAB%elms,& !DACB      D<A
                                !ODRESmat(idmat)%BATCH(DB_FIOD)%matAB%elms,ODmat%BATCH(CA_IOD)%matAB%elms,& !DBCA      D<B
                                !ODRESmat(idmat)%BATCH(CB_FIOD)%matAB%elms,ODmat%BATCH(DA_IOD)%matAB%elms,& !CBDA      C<B
                                iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,-KFAC)
                        ENDDO
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3BCAD(ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3CADB(ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3CBDA(ODRESmat(idmat)%BATCH(CB_IOD)%matAB%elms,ODmat%BATCH(DA_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE2SYM(&   ! 1 2 5 8
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
                             !ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
                             !ODRESmat(idmat)%BATCH(CB_IOD)%matAB%elms,ODmat%BATCH(DA_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,-KFAC)

                     ENDDO
                  ENDIF
               ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                  IF (startC.NE.startD)THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3BCAD(ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3ADBC(ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3BDAC(ODRESmat(idmat)%BATCH(BD_IOD)%matAB%elms,ODmat%BATCH(AC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE3(&   ! 1 2 3 4
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(BD_IOD)%matAB%elms,ODmat%BATCH(AC_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,-KFAC)

                     ENDDO
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3BCAD(ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE4(&   ! 1 2
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(BC_IOD)%matAB%elms,ODmat%BATCH(AD_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,-KFAC)
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            iOrbitalQ = iOrbitalQ + nOrbQ
         ENDDO
      ELSE ! (startA = startB)
         iOrbP=iOrbitalP-1+PASSPOFFSET
         iOrbitalQ = 1
         endAngmomQ = PQ%Q%p%nAngmom
         IF (PQ%samePQ) endAngmomQ = iAngmomP
         DO iAngmomQ=1,endAngmomQ
            DO iPassQ = 1, PQ%Q%p%nPasses
               nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
               iC = PQ%Q%p%indexAng1(iAngmomQ)
               iD = PQ%Q%p%indexAng2(iAngmomQ)
               nContC = PQ%Q%p%orbital1%nContracted(iC)
               nContD = PQ%Q%p%orbital2%nContracted(iD)
               nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
               nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
               PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
               startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
               startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
               BD_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startD))
               BC_IOD = ODmat%batchindex(ODmat%Aindex(startB),ODmat%Bindex(startC))
               DB_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startB))
               CB_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startB))
               AD_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startD))
               AC_IOD = ODmat%batchindex(ODmat%Aindex(startA),ODmat%Bindex(startC))
               CA_IOD = ODmat%batchindex(ODmat%Aindex(startC),ODmat%Bindex(startA))
               DA_IOD = ODmat%batchindex(ODmat%Aindex(startD),ODmat%Bindex(startA))
               IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                  IF (startC.NE.startD)THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3ADBC(ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3CADB(ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3DACB(ODRESmat(idmat)%BATCH(DA_IOD)%matAB%elms,ODmat%BATCH(CB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE5SYM(&   ! 1 3 5 6
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
                             !ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
                             !ODRESmat(idmat)%BATCH(DA_IOD)%matAB%elms,ODmat%BATCH(CB_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,-KFAC)

                     ENDDO
                     
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3CADB(ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE6(&   ! 1 5
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(CA_IOD)%matAB%elms,ODmat%BATCH(DB_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,-KFAC)

                     ENDDO
                  ENDIF
               ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                  IF (startC.NE.startD)THEN
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!                        CALL EXCHANGE_3ADBC(ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
                        CALL EXCHANGE_CASE7(&   ! 1 3
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
                             ODRESmat(idmat)%BATCH(AD_IOD)%matAB%elms,ODmat%BATCH(BC_IOD)%matAB%elms,&
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,-KFAC)

                     ENDDO
                  ELSE
                     iOrbQ=iOrbitalQ-1+PASSQOFFSET
                     DO idmat=1,Input%NDMAT_RHS
!                        CALL EXCHANGE_3ACBD(ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&
!                             & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)

                        !Fill in the entire Fmat diagonal block - this will be handled in the 
                        CALL EXCHANGE_CASE8(&   ! 1
                             ODRESmat(idmat)%BATCH(AC_IOD)%matAB%elms,ODmat%BATCH(BD_IOD)%matAB%elms,&      
                             iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat,-KFAC)

                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            iOrbitalQ = iOrbitalQ + nOrbQ
         ENDDO
      ENDIF
   ENDDO
      iOrbitalP = iOrbitalP + nOrbP
ENDDO

IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(1X,A,5I5)') 'ResultMat in Distribute_Exchange2',Output%ndim(1),Output%ndim(2),&
&Output%ndim(3),Output%ndim(4),Output%ndim(5)
  DO i5=1,Output%ndim(5)
    DO i4=1,Output%ndim(4)
      DO i3=1,Output%ndim(3)
        DO i2=1,Output%ndim(2)
          WRITE(LUPRI,'(3X,4I5)') i2,i3,i4,i5
          WRITE(LUPRI,'(5X,5F12.8)') (Output%resultMat(i1,i2,i3,i4,i5),i1=1,Output%ndim(1))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE Distribute_Exchange3SYM

SUBROUTINE EXCHANGE_CASE1(F1,D1,F2,D2,F3,D3,F4,D4,F5,D5,F6,D6,F7,D7,F8,D8,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint
REAL(REALK)    :: F1(:),D1(:),F2(:),D2(:),F3(:),D3(:),F4(:),D4(:)
REAL(REALK)    :: F5(:),D5(:),F6(:),D6(:),F7(:),D7(:),F8(:),D8(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
   TMP_B5=(iContB-1)*nAngB*nAngD*nContD
   TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
      TMP_A6=(iContA-1)*nAngD*nAngA*nContD
      TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
                  iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
                  iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
                     iDB5 = iDB5+1
                     iDA6 = iDA6+1
                     Delm1 = D1(iBD1)  
                     Delm2 = D2(iAD2)  
                     Delm5 = D5(iDB5)  
                     Delm8 = D8(iDA6)
                     SUM3=0
                     SUM7=0
                     SUM6=0
                     SUM4=0
                     iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
                     iCB6=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B6
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
                        iCA5 = iCA5+1
                        iCB6 = iCB6+1
                        SUM3 = SUM3+kint*D3(iBC2)  
                        SUM4 = SUM4+kint*D4(iAC1)  
                        SUM6 = SUM6+kint*D6(iCB6)  
                        SUM7 = SUM7+kint*D7(iCA5)  
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F2(iBC2) = F2(iBC2)+kint*Delm2
                        F5(iCA5) = F5(iCA5)+kint*Delm5
                        F8(iCB6) = F8(iCB6)+kint*Delm8
                     ENDDO
                     F3(iAD2) = F3(iAD2)+SUM3
                     F4(iBD1) = F4(iBD1)+SUM4  
                     F6(iDA6) = F6(iDA6)+SUM6
                     F7(iDB5) = F7(iDB5)+SUM7
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE1

SUBROUTINE EXCHANGE_CASE2(F1,D1,F2,D2,F5,D5,F8,D8,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint
REAL(REALK)    :: F1(:),D1(:),F2(:),D2(:)
REAL(REALK)    :: F5(:),D5(:),F8(:),D8(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
   TMP_B5=(iContB-1)*nAngB*nAngD*nContD
   TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
      TMP_A6=(iContA-1)*nAngD*nAngA*nContD
      TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
                  iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
                  iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
                     iDB5 = iDB5+1
                     iDA6 = iDA6+1
                     Delm1 = D1(iBD1)  
                     Delm2 = D2(iAD2)  
                     Delm5 = D5(iDB5)  
                     Delm8 = D8(iDA6)
                     iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
                     iCB6=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B6
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
                        iCA5 = iCA5+1
                        iCB6 = iCB6+1
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F2(iBC2) = F2(iBC2)+kint*Delm2
                        F5(iCA5) = F5(iCA5)+kint*Delm5
                        F8(iCB6) = F8(iCB6)+kint*Delm8
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE2

SUBROUTINE EXCHANGE_CASE3(F1,D1,F2,D2,F3,D3,F4,D4,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F2(:),D2(:),F3(:),D3(:),F4(:),D4(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
   TMP_B5=(iContB-1)*nAngB*nAngD*nContD
   TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
      TMP_A6=(iContA-1)*nAngD*nAngA*nContD
      TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
                  iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
                  iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
                     iDB5 = iDB5+1
                     iDA6 = iDA6+1
                     Delm1 = D1(iBD1)  
                     Delm2 = D2(iAD2)  
                     SUM3=0
                     SUM4=0
                     iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
                     iCB6=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B6
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
                        iCA5 = iCA5+1
                        iCB6 = iCB6+1
                        SUM3 = SUM3+kint*D3(iBC2)  
                        SUM4 = SUM4+kint*D4(iAC1)  
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F2(iBC2) = F2(iBC2)+kint*Delm2
                     ENDDO
                     F3(iAD2) = F3(iAD2)+SUM3
                     F4(iBD1) = F4(iBD1)+SUM4  
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE3

SUBROUTINE EXCHANGE_CASE4(F1,D1,F2,D2,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F2(:),D2(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
   TMP_B5=(iContB-1)*nAngB*nAngD*nContD
   TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
      TMP_A6=(iContA-1)*nAngD*nAngA*nContD
      TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
                  iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
                  iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
                     Delm1 = D1(iBD1)  
                     Delm2 = D2(iAD2)  
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F2(iBC2) = F2(iBC2)+kint*Delm2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE4

SUBROUTINE EXCHANGE_CASE5(F1,D1,F3,D3,F5,D5,F6,D6,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F3(:),D3(:)
REAL(REALK)    :: F5(:),D5(:),F6(:),D6(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
   TMP_B5=(iContB-1)*nAngB*nAngD*nContD
   TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
      TMP_A6=(iContA-1)*nAngD*nAngA*nContD
      TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
                  iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
                  iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
                     iDB5 = iDB5+1
                     iDA6 = iDA6+1
                     Delm1 = D1(iBD1)  
                     Delm5 = D5(iDB5)  
                     SUM3=0
                     SUM6=0
                     iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
                     iCB6=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B6
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
                        iCA5 = iCA5+1
                        iCB6 = iCB6+1
                        SUM3 = SUM3+kint*D3(iBC2)  
                        SUM6 = SUM6+kint*D6(iCB6)  
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F5(iCA5) = F5(iCA5)+kint*Delm5
                     ENDDO
                     F3(iAD2) = F3(iAD2)+SUM3
                     F6(iDA6) = F6(iDA6)+SUM6
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE5

SUBROUTINE EXCHANGE_CASE6(F1,D1,F5,D5,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:)
REAL(REALK)    :: F5(:),D5(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
   TMP_B5=(iContB-1)*nAngB*nAngD*nContD
   TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
      TMP_A6=(iContA-1)*nAngD*nAngA*nContD
      TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
                  iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
                  iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iDB5 = iDB5+1
                     Delm1 = D1(iBD1)  
                     Delm5 = D5(iDB5)  
                     iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iCA5 = iCA5+1
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F5(iCA5) = F5(iCA5)+kint*Delm5
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE6

SUBROUTINE EXCHANGE_CASE7(F1,D1,F3,D3,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F3(:),D3(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
 !  TMP_B5=(iContB-1)*nAngB*nAngD*nContD
 !  TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
 !     TMP_A6=(iContA-1)*nAngD*nAngA*nContD
 !     TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
!                  iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
!                  iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     Delm1 = D1(iBD1)  
                     SUM3=0
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
                        SUM3 = SUM3+kint*D3(iBC2)  
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                     ENDDO
                     F3(iAD2) = F3(iAD2)+SUM3
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE7

SUBROUTINE EXCHANGE_CASE8(F1,D1,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B1=(iContB-1)*nAngB*nAngD
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC 
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     Delm1 = D1(iBD1)  
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE8

SUBROUTINE EXCHANGE_CASE1SYM(F1,D1,F2,D2,F3,D3,F4,D4,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F2(:),D2(:),F3(:),D3(:),F4(:),D4(:)
!REAL(REALK)    :: F5(:),D5(:),F6(:),D6(:),F7(:),D7(:),F8(:),D8(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
!   TMP_B5=(iContB-1)*nAngB*nAngD*nContD
!   TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
!      TMP_A6=(iContA-1)*nAngD*nAngA*nContD
!      TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
 !                 iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
 !                 iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
 !                    iDB5 = iDB5+1
 !                    iDA6 = iDA6+1
                     Delm1 = D1(iBD1)  
                     Delm2 = D2(iAD2)  
 !                    Delm5 = D5(iDB5)  
 !                    Delm8 = D8(iDA6)
                     SUM3=0.D0
 !                    SUM7=0
 !                    SUM6=0
                     SUM4=0.D0
 !                    iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
 !                    iCB6=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B6
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = KFAC*2*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
 !                       iCA5 = iCA5+1
 !                       iCB6 = iCB6+1
                        SUM3 = SUM3+kint*D3(iBC2)  
                        SUM4 = SUM4+kint*D4(iAC1)  
!                        SUM6 = SUM6+kint*D6(iCB6)  
!                        SUM7 = SUM7+kint*D7(iCA5)  
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F2(iBC2) = F2(iBC2)+kint*Delm2
!                        F5(iCA5) = F5(iCA5)+kint*Delm5
!                        F8(iCB6) = F8(iCB6)+kint*Delm8
                     ENDDO
                     F3(iAD2) = F3(iAD2)+SUM3
                     F4(iBD1) = F4(iBD1)+SUM4  
!                     F6(iDA6) = F6(iDA6)+SUM6
!                     F7(iDB5) = F7(iDB5)+SUM7
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE1SYM

SUBROUTINE EXCHANGE_CASE2SYM(F1,D1,F2,D2,&!F5,D5,F8,D8,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F2(:),D2(:)
!REAL(REALK)    :: F5(:),D5(:),F8(:),D8(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
 !  TMP_B5=(iContB-1)*nAngB*nAngD*nContD
 !  TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
 !     TMP_A6=(iContA-1)*nAngD*nAngA*nContD
 !     TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
 !                 iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
 !                 iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
!                     iDB5 = iDB5+1
!                     iDA6 = iDA6+1
                     Delm1 = D1(iBD1)  
                     Delm2 = D2(iAD2)  
!                     Delm5 = D5(iDB5)  
!                     Delm8 = D8(iDA6)
!                     iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
!                     iCB6=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B6
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = 2*KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
!                        iCA5 = iCA5+1
!                        iCB6 = iCB6+1
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F2(iBC2) = F2(iBC2)+kint*Delm2
!                        F5(iCA5) = F5(iCA5)+kint*Delm5
!                        F8(iCB6) = F8(iCB6)+kint*Delm8
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE2SYM

SUBROUTINE EXCHANGE_CASE3SYM(F1,D1,F2,D2,F3,D3,F4,D4,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F2(:),D2(:),F3(:),D3(:),F4(:),D4(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
   !TMP_B5=(iContB-1)*nAngB*nAngD*nContD
   !TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
      !TMP_A6=(iContA-1)*nAngD*nAngA*nContD
      !TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
                  !iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
                  !iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
                     !iDB5 = iDB5+1
                     !iDA6 = iDA6+1
                     Delm1 = D1(iBD1)  
                     Delm2 = D2(iAD2)  
                     SUM3=0
                     SUM4=0
                     !iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
                     !iCB6=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B6
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = 2*KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
                        !iCA5 = iCA5+1
                        !iCB6 = iCB6+1
                        SUM3 = SUM3+kint*D3(iBC2)  
                        SUM4 = SUM4+kint*D4(iAC1)  
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F2(iBC2) = F2(iBC2)+kint*Delm2
                     ENDDO
                     F3(iAD2) = F3(iAD2)+SUM3
                     F4(iBD1) = F4(iBD1)+SUM4  
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE3SYM

SUBROUTINE EXCHANGE_CASE4SYM(F1,D1,F2,D2,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F2(:),D2(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
 !  TMP_B5=(iContB-1)*nAngB*nAngD*nContD
 !  TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
 !     TMP_A6=(iContA-1)*nAngD*nAngA*nContD
 !     TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
!                  iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
!                  iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
                     Delm1 = D1(iBD1)  
                     Delm2 = D2(iAD2)  
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = 2*KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                        F2(iBC2) = F2(iBC2)+kint*Delm2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE4SYM

SUBROUTINE EXCHANGE_CASE5SYM(F1,D1,F3,D3,&!F5,D5,F6,D6,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F3(:),D3(:)
!REAL(REALK)    :: F5(:),D5(:),F6(:),D6(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
 !  TMP_B5=(iContB-1)*nAngB*nAngD*nContD
 !  TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
  !    TMP_A6=(iContA-1)*nAngD*nAngA*nContD
  !    TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
  !                iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
  !                iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     iAD2 = (iAngD-1)*nAngA+TMP_D2
  !                   iDB5 = iDB5+1
  !                   iDA6 = iDA6+1
                     Delm1 = D1(iBD1)  
  !                   Delm5 = D5(iDB5)  
                     SUM3=0
                     SUM6=0
 !                    iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
 !                    iCB6=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B6
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = 2*KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
 !                       iCA5 = iCA5+1
 !                       iCB6 = iCB6+1
                        SUM3 = SUM3+kint*D3(iBC2)  
!                        SUM6 = SUM6+kint*D6(iCB6)  
                        F1(iAC1) = F1(iAC1)+kint*Delm1
!                        F5(iCA5) = F5(iCA5)+kint*Delm5
                     ENDDO
                     F3(iAD2) = F3(iAD2)+SUM3
!                     F6(iDA6) = F6(iDA6)+SUM6
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE5SYM

SUBROUTINE EXCHANGE_CASE6SYM(F1,D1,&!F5,D5,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:)
!REAL(REALK)    :: F5(:),D5(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
!   TMP_B5=(iContB-1)*nAngB*nAngD*nContD
!   TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
!      TMP_A6=(iContA-1)*nAngD*nAngA*nContD
!      TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
  !                iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
  !                iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
 !                    iDB5 = iDB5+1
                     Delm1 = D1(iBD1)  
 !                    Delm5 = D5(iDB5)  
 !                    iCA5=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A5
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint =2*KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
!                        iCA5 = iCA5+1
                        F1(iAC1) = F1(iAC1)+kint*Delm1
!                        F5(iCA5) = F5(iCA5)+kint*Delm5
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE6SYM

SUBROUTINE EXCHANGE_CASE7SYM(F1,D1,F3,D3,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:),F3(:),D3(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B2=(iContB-1)*nAngB*nAngC    !tmp_b3=tmp_b2
   TMP_B1=(iContB-1)*nAngB*nAngD
!   TMP_B5=(iContB-1)*nAngB*nAngD*nContD
!   TMP_B6=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      TMP_A2=(iContA-1)*nAngA*nAngD
!      TMP_A6=(iContA-1)*nAngD*nAngA*nContD
!      TMP_A5=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D2 = (iContD-1)*nAngA*nAngD*nContA+TMP_A2+iAngA
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  TMP_C2=(iContC-1)*nAngB*nAngC*nContB+TMP_B2+iAngB
!                  iDB5=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B5
!                  iDA6=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A6
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     Delm1 = D1(iBD1)  
                     SUM3=0
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = 2*KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        iBC2 = (iAngC-1)*nAngB+TMP_C2   
                        SUM3 = SUM3+kint*D3(iBC2)  
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                     ENDDO
                     F3(iAD2) = F3(iAD2)+SUM3
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE7SYM

SUBROUTINE EXCHANGE_CASE8SYM(F1,D1,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat,KFAC)
IMPLICIT NONE
REAL(REALK)    :: QPmat(:,:),kint,KFAC
REAL(REALK)    :: F1(:),D1(:)
INTEGER        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
REAL(REALK)    :: Delm1,Delm2,Delm5,Delm8,SUM3,SUM4,SUM6,SUM7
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: TMP_B1,TMP_A1,TMP_D1,TMP_C1,iOrbP,iOrbQ,iBD5,iDA6,iBD1
INTEGER        :: TMP_B2,TMP_A2,TMP_D2,TMP_C2,iAD2,iDB5
INTEGER        :: TMP_B5,TMP_A5,iCA5,iCB6
INTEGER        :: TMP_B6,TMP_A6,iAC1,iBC2

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB


!Fill in the entire Fmat this w
iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B1=(iContB-1)*nAngB*nAngD
   DO iContA=1,nContA
      TMP_A1=(iContA-1)*nAngA*nAngC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D1 = (iContD-1)*nAngB*nAngD*nContB+TMP_B1+iAngB
               DO iContC=1,nContC
                  TMP_C1=(iContC-1)*nAngA*nAngC*nContA+TMP_A1+iAngA
                  DO iAngD=1,nAngD
                     iBD1 = (iAngD-1)*nAngB+TMP_D1
                     Delm1 = D1(iBD1)  
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = 2*KFAC*QPmat(iOrbQ,iOrbP)
                        iAC1 = (iAngC-1)*nAngA+TMP_C1   
                        F1(iAC1) = F1(iAC1)+kint*Delm1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_CASE8SYM
!-----------------------------
SUBROUTINE EXCHANGE_3ACBD(F,D,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
IMPLICIT NONE
REAL(REALK)        :: F(:),D(:),QPmat(:,:),KFAC,Delm,kint
INTEGER            :: START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
INTEGER            :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER            :: TMP_B,TMP_A,iOrbP,iOrbQ,TMP_D,TMP_C,iBD,iAC
!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B=(iContB-1)*nAngB*nAngD
   DO iContA=1,nContA
      TMP_A=(iContA-1)*nAngA*nAngC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D = (iContD-1)*nAngB*nAngD*nContB+TMP_B+iAngB
               DO iContC=1,nContC
                  TMP_C=(iContC-1)*nAngA*nAngC*nContA+TMP_A+iAngA
                  DO iAngD=1,nAngD
                     iBD = (iAngD-1)*nAngB+TMP_D
                     Delm = D(iBD)  
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)                              
                        iAC = (iAngC-1)*nAngA+TMP_C   
                        F(iAC) = F(iAC)+kint*Delm
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_3ACBD

SUBROUTINE EXCHANGE_3BCAD(F,D,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!!$ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)= ODRESmat(idmat)%BATCH(BC_FIOD(idmat))%matAB%elms(iBC)-  kint*&!D(iA1,iD1)
!!$     &ODmat%BATCH(AD_IOD)%matAB%elms(iAD)
IMPLICIT NONE
REAL(REALK)        :: F(:),D(:),QPmat(:,:),KFAC,Delm,kint
INTEGER            :: START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
INTEGER            :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER            :: TMP_B,TMP_A,iOrbP,iOrbQ,TMP_D,TMP_C,iAD,iBC

!!$  iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
!!$  iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B=(iContB-1)*nAngB*nAngC
   DO iContA=1,nContA
      TMP_A=(iContA-1)*nAngA*nAngD
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D = (iContD-1)*nAngA*nAngD*nContA+TMP_A+iAngA
               DO iContC=1,nContC
                  TMP_C=(iContC-1)*nAngB*nAngC*nContB+TMP_B+iAngB
                  DO iAngD=1,nAngD
                     iAD = (iAngD-1)*nAngA+TMP_D
                     Delm = D(iAD)  
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)                              
                        iBC = (iAngC-1)*nAngB+TMP_C   
                        F(iBC) = F(iBC)+kint*Delm
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_3BCAD

SUBROUTINE EXCHANGE_3ADBC(F,D,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!!$ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)= ODRESmat(idmat)%BATCH(AD_FIOD(idmat))%matAB%elms(iAD)-  kint*&!D(iB1,iC1)
!!$     &ODmat%BATCH(BC_IOD)%matAB%elms(iBC)
IMPLICIT NONE
REAL(REALK)        :: F(:),D(:),QPmat(:,:),KFAC,SUM,kint
INTEGER            :: START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
INTEGER            :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER            :: TMP_B,TMP_A,iOrbP,iOrbQ,TMP_D,TMP_C,iAD,iBC


!!$   iAD = iAngA+(iAngD-1)*nAngA+(iContA-1)*nAngA*nAngD+(iContD-1)*nAngA*nAngD*nContA
!!$   iBC = iAngB+(iAngC-1)*nAngB+(iContB-1)*nAngB*nAngC+(iContC-1)*nAngB*nAngC*nContB
iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B=(iContB-1)*nAngB*nAngC
   DO iContA=1,nContA
      TMP_A=(iContA-1)*nAngA*nAngD
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D = (iContD-1)*nAngA*nAngD*nContA+TMP_A+iAngA
               DO iContC=1,nContC
                  TMP_C=(iContC-1)*nAngB*nAngC*nContB+TMP_B+iAngB
                  DO iAngD=1,nAngD
                     iAD = (iAngD-1)*nAngA+TMP_D
                     SUM=0
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)                              
                        iBC = (iAngC-1)*nAngB+TMP_C   
                        SUM = SUM+kint*D(iBC)  
                     ENDDO
                     F(iAD) = F(iAD)+SUM
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_3ADBC

SUBROUTINE EXCHANGE_3BDAC(F,D,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
IMPLICIT NONE
REAL(REALK)        :: F(:),D(:),QPmat(:,:),KFAC,SUM,kint
INTEGER            :: START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
INTEGER            :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER            :: TMP_B,TMP_A,iOrbP,iOrbQ,TMP_D,TMP_C,iBD,iAC

!!$  iAC = iAngA+(iAngC-1)*nAngA+(iContA-1)*nAngA*nAngC+(iContC-1)*nAngA*nAngC*nContA
!!$  iBD = iAngB+(iAngD-1)*nAngB+(iContB-1)*nAngB*nAngD+(iContD-1)*nAngB*nAngD*nContB
iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B=(iContB-1)*nAngB*nAngD
   DO iContA=1,nContA
      TMP_A=(iContA-1)*nAngA*nAngC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               TMP_D = (iContD-1)*nAngB*nAngD*nContB+TMP_B+iAngB
               DO iContC=1,nContC
                  DO iAngD=1,nAngD
                     iBD = (iAngD-1)*nAngB+TMP_D
                     TMP_C = (iContC-1)*nAngA*nAngC*nContA+TMP_A+iAngA    
                     SUM = 0
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)                              
                        iAC = (iAngC-1)*nAngA+TMP_C
                        SUM=SUM+kint*D(iAC)  
                     ENDDO
                        F(iBD) = F(iBD)+SUM  
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_3BDAC

SUBROUTINE EXCHANGE_3CADB(F,D,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!!$ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)= ODRESmat(idmat)%BATCH(CA_FIOD(idmat))%matAB%elms(iCA)-  kint*&!D(iD1,iB1)
!!$     &ODmat%BATCH(DB_IOD)%matAB%elms(iDB)
IMPLICIT NONE
REAL(REALK)        :: F(:),D(:),QPmat(:,:),KFAC,Delm,kint
INTEGER            :: START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
INTEGER            :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER            :: TMP_B,TMP_A,iOrbP,iOrbQ,TMP_D,TMP_C,iCA,iDB

!!$ iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC
!!$ iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B=(iContB-1)*nAngB*nAngD*nContD
   DO iContA=1,nContA
      TMP_A=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               DO iContC=1,nContC
                  iDB=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B
                  DO iAngD=1,nAngD
                     iDB = iDB+1
                     iCA=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A
                     Delm = D(iDB)  
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)                              
                        iCA = iCA+1
                        F(iCA) = F(iCA)+kint*Delm
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_3CADB

SUBROUTINE EXCHANGE_3DACB(F,D,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!!$ODRESmat(idmat)%BATCH(DA_FIOD(idmat))%matAB%elms(iDA)= ODRESmat(idmat)%BATCH(DA_FIOD(idmat))%matAB%elms(iDA)-  kint*&!D(iC1,iB1)
!!$     &ODmat%BATCH(CB_IOD)%matAB%elms(iCB)
IMPLICIT NONE
REAL(REALK)        :: F(:),D(:),QPmat(:,:),KFAC,SUM,kint
INTEGER            :: START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
INTEGER            :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER            :: TMP_B,TMP_A,iOrbP,iOrbQ,TMP_D,TMP_C,iDA,iCB

!!$  iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD
!!$  iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A=(iContA-1)*nAngD*nAngA*nContD
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               DO iContC=1,nContC
                  iDA=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A
                  DO iAngD=1,nAngD
                     iDA = iDA+1
                     iCB=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B
                     SUM=0
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)      
                        iCB = iCB+1
                        SUM=SUM+kint*D(iCB)  
                     ENDDO
                     F(iDA) = F(iDA)+SUM
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_3DACB

SUBROUTINE EXCHANGE_3DBCA(F,D,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!!$ODRESmat(idmat)%BATCH(DB_FIOD(idmat))%matAB%elms(iDB)= ODRESmat(idmat)%BATCH(DB_FIOD(idmat))%matAB%elms(iDB)-  kint*&!D(iC1,iA1)
!!$     &ODmat%BATCH(CA_IOD)%matAB%elms(iCA)
IMPLICIT NONE
REAL(REALK)        :: F(:),D(:),QPmat(:,:),KFAC,kint,SUM
INTEGER            :: START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
INTEGER            :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER            :: TMP_B,TMP_A,iOrbP,iOrbQ,TMP_D,TMP_C,iDB,iCA

!!$   iDB = iAngD+(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+(iContB-1)*nAngD*nAngB*nContD
!!$   iCA = iAngC+(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+(iContA-1)*nAngC*nAngA*nContC

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B=(iContB-1)*nAngB*nAngD*nContD
   DO iContA=1,nContA
      TMP_A=(iContA-1)*nAngA*nAngC*nContC
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               DO iContC=1,nContC
                  iDB=(iAngB-1)*nAngD+(iContD-1)*nAngD*nAngB+TMP_B
                  DO iAngD=1,nAngD
                     iDB = iDB+1
                     iCA=(iAngA-1)*nAngC+(iContC-1)*nAngC*nAngA+TMP_A
                     SUM=0
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)                              
                        iCA = iCA+1
                        SUM = SUM +kint*D(iCA)  
                     ENDDO
                     F(iDB) = F(iDB)+SUM
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_3DBCA

SUBROUTINE EXCHANGE_3CBDA(F,D,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat)
!!$ODRESmat(idmat)%BATCH(CB_FIOD(idmat))%matAB%elms(iCB)= ODRESmat(idmat)%BATCH(CB_FIOD(idmat))%matAB%elms(iCB)-  kint*&!D(iD1,iA1)
!!$     &ODmat%BATCH(DA_IOD)%matAB%elms(iDA)
IMPLICIT NONE
REAL(REALK)        :: F(:),D(:),QPmat(:,:),KFAC,SUM,kint,Delm
INTEGER            :: START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC
!
INTEGER            :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER            :: TMP_B,TMP_A,iOrbP,iOrbQ,TMP_D,TMP_C,iDA,iCB

!!$  iCB = iAngC+(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+(iContB-1)*nAngC*nAngB*nContC
!!$  iDA = iAngD+(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+(iContA-1)*nAngD*nAngA*nContD

iOrbP = START_iOrbP
DO iContB=1,nContB
   TMP_B=(iContB-1)*nAngC*nAngB*nContC
   DO iContA=1,nContA
      TMP_A=(iContA-1)*nAngD*nAngA*nContD
      DO iAngB=1,nAngB
         DO iAngA=1,nAngA
            iOrbP=iOrbP+1
            iOrbQ = START_iOrbQ
            DO iContD=1,nContD
               DO iContC=1,nContC
                  iDA=(iAngA-1)*nAngD+(iContD-1)*nAngD*nAngA+TMP_A
                  DO iAngD=1,nAngD
                     iDA = iDA+1
                     iCB=(iAngB-1)*nAngC+(iContC-1)*nAngC*nAngB+TMP_B
                     Delm = D(iDA)
                     DO iAngC=1,nAngC
                        iOrbQ=iOrbQ+1
                        kint = QPmat(iOrbQ,iOrbP)                              
                        iCB = iCB+1
                        F(iCB) = F(iCB) + kint*Delm
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_3CBDA

SUBROUTINE Distribute_Jengine1(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: iC,iD,iOrbQ,idmat
real(realk) :: Value
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ

!
nAngmomP = PQ%P%p%nAngmom

nderivP=Input%nDerivP
totOrbP = PQ%P%p%totOrbitals
!
iOrbitalP = 1
DO iDerivP=1,NDerivP
   DO iAngmomP=1,nAngmomP
      DO iPassP = 1, PQ%P%p%nPasses
      nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
      totOrbQ = PQ%Q%p%totOrbitals
      endOrbP = iOrbitalP+nOrbP-1
!
      iA = PQ%P%p%indexAng1(iAngmomP)
      iB = PQ%P%p%indexAng2(iAngmomP)
      nContA = PQ%P%p%orbital1%nContracted(iA)
      nContB = PQ%P%p%orbital2%nContracted(iB)
      nAngA = PQ%P%p%orbital1%nOrbComp(iA)
      nAngB = PQ%P%p%orbital2%nOrbComp(iB)
      PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
      startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
      startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
      IF(startA.NE.startB)THEN
         iOrbP=iOrbitalP-1+PASSPOFFSET
         iContP=0
         DO iContB=1,nContB
            DO iContA=1,nContA
               iContP=iContP+1
               DO iAngB=1,nAngB
                  iB1=startB-1+iAngB+(iContB-1)*nAngB
                  DO iAngA=1,nAngA
                     iA1=startA-1+iAngA+(iContA-1)*nAngA
                     iOrbP=iOrbP+1
                     iOrbitalQ = 1
                     endAngmomQ = PQ%Q%p%nAngmom
                     IF (PQ%samePQ) endAngmomQ = iAngmomP
                     nderivQ=input%nderivQ
                     DO iderivQ =1,nderivQ
                        DO iAngmomQ=1,endAngmomQ
                        DO iPassQ = 1, PQ%Q%p%nPasses
                           nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                           endOrbQ = iOrbitalQ+nOrbQ-1
                           !  ideriv=iderivP+(iderivQ-1)*NderivP
                           iC = PQ%Q%p%indexAng1(iAngmomQ)
                           iD = PQ%Q%p%indexAng2(iAngmomQ)
                           nContC = PQ%Q%p%orbital1%nContracted(iC)
                           nContD = PQ%Q%p%orbital2%nContracted(iD)
                           nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                           nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                           PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                           startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                           startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                           !
                           ideriv=iderivP+(iderivQ-1)*NderivP
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              DO iContC=1,nContC
                                 iContQ = (iContD-1)*nContC + iContC
                                 DO iAngD=1,nAngD
                                    iD1=startD-1+iAngD+(iContD-1)*nAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC-1+iAngC+(iContC-1)*nAngC
                                       iOrbQ=iOrbQ+1
                                       Value = QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat)&
                                               & + Value
                                          Output%ResultMat(iB1,iA1,1,1,idmat) = Output%ResultMat(iB1,iA1,1,1,idmat)&
                                               & + Value
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !
                        ENDDO
                           iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
        ! iOrbitalP = iOrbitalP + nOrbP
      ELSE
         iOrbP=iOrbitalP-1+PASSPOFFSET
         iContP=0
         DO iContB=1,nContB
            DO iContA=1,nContA
               iContP=iContP+1
               DO iAngB=1,nAngB
                  iB1=startB-1+iAngB+(iContB-1)*nAngB
                  DO iAngA=1,nAngA
                     iA1=startA-1+iAngA+(iContA-1)*nAngA
                     iOrbP=iOrbP+1
                     iOrbitalQ = 1
                     endAngmomQ = PQ%Q%p%nAngmom
                     IF (PQ%samePQ) endAngmomQ = iAngmomP
                     nderivQ=input%nderivQ
                     DO iderivQ =1,nderivQ
                        DO iAngmomQ=1,endAngmomQ
                        DO iPassQ = 1, PQ%Q%p%nPasses
                           nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                           endOrbQ = iOrbitalQ+nOrbQ-1
                           !  ideriv=iderivP+(iderivQ-1)*NderivP
                           iC = PQ%Q%p%indexAng1(iAngmomQ)
                           iD = PQ%Q%p%indexAng2(iAngmomQ)
                           nContC = PQ%Q%p%orbital1%nContracted(iC)
                           nContD = PQ%Q%p%orbital2%nContracted(iD)
                           nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                           nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                           PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                           startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                           startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                           !
                           ideriv=iderivP+(iderivQ-1)*NderivP
                           iOrbQ=iOrbitalQ-1+PASSQOFFSET
                           DO iContD=1,nContD
                              DO iContC=1,nContC
                                 iContQ = (iContD-1)*nContC + iContC
                                 DO iAngD=1,nAngD
                                    iD1=startD-1+iAngD+(iContD-1)*nAngD
                                    DO iAngC=1,nAngC
                                       iC1=startC-1+iAngC+(iContC-1)*nAngC
                                       iOrbQ=iOrbQ+1
                                       Value = QPmat(iOrbQ,iOrbP)
                                       DO idmat=1,Input%NDMAT_RHS 
                                          Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat) + Value
                                       ENDDO
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                           !
                        ENDDO
                           iOrbitalQ = iOrbitalQ + nOrbQ
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !iOrbitalP = iOrbitalP + nOrbP
      ENDIF
   ENDDO
      iOrbitalP = iOrbitalP + nOrbP
   ENDDO
ENDDO
IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(1X,A,5I5)') 'ResultMat in DistributePQint',Output%ndim(1),Output%ndim(2),&
&Output%ndim(3),Output%ndim(4),Output%ndim(5)
  DO i5=1,Output%ndim(5)
    DO i4=1,Output%ndim(4)
      DO i3=1,Output%ndim(3)
        DO i2=1,Output%ndim(2)
          WRITE(LUPRI,'(3X,4I5)') i2,i3,i4,i5
          WRITE(LUPRI,'(5X,5F12.8)') (Output%resultMat(i1,i2,i3,i4,i5),i1=1,Output%ndim(1))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE Distribute_Jengine1

SUBROUTINE Distribute_Jengine2(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
Integer :: iA1,iB1,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ,iC1,iD1
Integer :: iContQ,iContC,iContD,iAngC,iAngD,nContC,nContD,nAngC,nAngD,startC,startD
Integer :: iC,iD,iOrbQ,idmat
real(realk) :: Value
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ

!
nAngmomP = PQ%P%p%nAngmom

nderivP=Input%nDerivP
totOrbP = PQ%P%p%totOrbitals
!
iOrbitalP = 1
DO iDerivP=1,NDerivP
   DO iAngmomP=1,nAngmomP
      DO iPassP = 1, PQ%P%p%nPasses
      nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
      totOrbQ = PQ%Q%p%totOrbitals
      endOrbP = iOrbitalP+nOrbP-1
      iA = PQ%P%p%indexAng1(iAngmomP)
      iB = PQ%P%p%indexAng2(iAngmomP)
      nContA = PQ%P%p%orbital1%nContracted(iA)
      nContB = PQ%P%p%orbital2%nContracted(iB)
      nAngA = PQ%P%p%orbital1%nOrbComp(iA)
      nAngB = PQ%P%p%orbital2%nOrbComp(iB)
      PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
      startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
      startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
      iOrbP=iOrbitalP-1+PASSPOFFSET
      iContP=0
      DO iContB=1,nContB
         DO iContA=1,nContA
            iContP=iContP+1
            DO iAngB=1,nAngB
               iB1=startB-1+iAngB+(iContB-1)*nAngB
               DO iAngA=1,nAngA
                  iA1=startA-1+iAngA+(iContA-1)*nAngA
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = PQ%Q%p%nAngmom
                  IF (PQ%samePQ) endAngmomQ = iAngmomP
                  nderivQ=input%nderivQ
                  DO iderivQ =1,nderivQ
                     DO iAngmomQ=1,endAngmomQ
                     DO iPassQ = 1, PQ%Q%p%nPasses
                        nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
                        endOrbQ = iOrbitalQ+nOrbQ-1
                        !  ideriv=iderivP+(iderivQ-1)*NderivP
                        iC = PQ%Q%p%indexAng1(iAngmomQ)
                        iD = PQ%Q%p%indexAng2(iAngmomQ)
                        nContC = PQ%Q%p%orbital1%nContracted(iC)
                        nContD = PQ%Q%p%orbital2%nContracted(iD)
                        nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
                        nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
                        PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                        startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
                        startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
                        ideriv=iderivP+(iderivQ-1)*NderivP
                        iOrbQ=iOrbitalQ-1+PASSQOFFSET
                        DO iContD=1,nContD
                           DO iContC=1,nContC
                              iContQ = (iContD-1)*nContC + iContC
                              DO iAngD=1,nAngD
!                                 iD=startD-1+iAngD+(iContD-1)*nAngD
                                 DO iAngC=1,nAngC
!                                    iC=startC-1+iAngC+(iContC-1)*nAngC
                                    iOrbQ=iOrbQ+1
                                    Value = QPmat(iOrbQ,iOrbP)
                                    DO idmat=1,Input%NDMAT_RHS 
                                       Output%ResultMat(iA1,iB1,1,1,idmat) = Output%ResultMat(iA1,iB1,1,1,idmat) + Value
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                        iOrbitalQ = iOrbitalQ + nOrbQ
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
      iOrbitalP = iOrbitalP + nOrbP
   ENDDO
ENDDO
IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(1X,A,5I5)') 'ResultMat in DistributePQint',Output%ndim(1),Output%ndim(2),&
&Output%ndim(3),Output%ndim(4),Output%ndim(5)
  DO i5=1,Output%ndim(5)
    DO i4=1,Output%ndim(4)
      DO i3=1,Output%ndim(3)
        DO i2=1,Output%ndim(2)
          WRITE(LUPRI,'(3X,4I5)') i2,i3,i4,i5
          WRITE(LUPRI,'(5X,5F12.8)') (Output%resultMat(i1,i2,i3,i4,i5),i1=1,Output%ndim(1))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE Distribute_Jengine2

SUBROUTINE DistributePQ_CS1(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
real(realk) :: Value
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ

nAngmomP = PQ%P%p%nAngmom

nderivP=Input%nDerivP
totOrbP = PQ%P%p%totOrbitals
!
iOrbitalP = 1
DO iDerivP=1,NDerivP
   DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
      nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
      totOrbQ = PQ%Q%p%totOrbitals
      endOrbP = iOrbitalP+nOrbP-1
      iA = PQ%P%p%indexAng1(iAngmomP)
      iB = PQ%P%p%indexAng2(iAngmomP)
      nContA = PQ%P%p%orbital1%nContracted(iA)
      nContB = PQ%P%p%orbital2%nContracted(iB)
      nAngA = PQ%P%p%orbital1%nOrbComp(iA)
      nAngB = PQ%P%p%orbital2%nOrbComp(iB)
      PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
      startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
      startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
      iOrbP=iOrbitalP-1+PASSPOFFSET
      iContP=0
      DO iContB=1,nContB
         DO iContA=1,nContA
            iContP=iContP+1
            DO iAngB=1,nAngB
               iB=startB-1+iAngB+(iContB-1)*nAngB
               DO iAngA=1,nAngA
                  iA=startA-1+iAngA+(iContA-1)*nAngA
                  iOrbP=iOrbP+1                 
                  nderivQ=input%nderivQ
                  DO iderivQ =1,nderivQ
                     ideriv=iderivP+(iderivQ-1)*NderivP
                     Value = sqrt(QPmat(iOrbP,iorbP))
                     Output%ResultMat(iA,iB,1,1,ideriv) = Value
                     Output%ResultMat(iB,iA,1,1,ideriv) = Value
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
      iOrbitalP = iOrbitalP + nOrbP
   ENDDO
ENDDO

END SUBROUTINE DistributePQ_CS1

SUBROUTINE DistributePQ_CS2(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iContB,iContA,iAngB,iAngA,nderivQ,iderivQ,ideriv
integer :: iPassP,iPassQ
integer :: PASSPOFFSET,PASSQOFFSET
!
nAngmomP = PQ%P%p%nAngmom

nderivP=Input%nDerivP
totOrbP = PQ%P%p%totOrbitals
!
iOrbitalP = 1
DO iDerivP=1,NDerivP
   DO iAngmomP=1,nAngmomP
   DO iPassP = 1, PQ%P%p%nPasses
      nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
      totOrbQ = PQ%Q%p%totOrbitals
      endOrbP = iOrbitalP+nOrbP-1
      iA = PQ%P%p%indexAng1(iAngmomP)
      iB = PQ%P%p%indexAng2(iAngmomP)
      nContA = PQ%P%p%orbital1%nContracted(iA)
      nContB = PQ%P%p%orbital2%nContracted(iB)
      nAngA = PQ%P%p%orbital1%nOrbComp(iA)
      nAngB = PQ%P%p%orbital2%nOrbComp(iB)
      PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
      startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
      startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
      iOrbP=iOrbitalP-1+PASSPOFFSET
      iContP=0
      DO iContB=1,nContB
         DO iContA=1,nContA
            iContP=iContP+1
            DO iAngB=1,nAngB
               iB=startB-1+iAngB+(iContB-1)*nAngB
               DO iAngA=1,nAngA
                  iA=startA-1+iAngA+(iContA-1)*nAngA
                  iOrbP=iOrbP+1                 
                  nderivQ=input%nderivQ
                  DO iderivQ =1,nderivQ
                     ideriv=iderivP+(iderivQ-1)*NderivP
                     Output%ResultMat(iA,iB,1,1,ideriv) = &
                          &sqrt(QPmat(iOrbP,iorbP))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
      iOrbitalP = iOrbitalP + nOrbP
   ENDDO
ENDDO

END SUBROUTINE DistributePQ_CS2

SUBROUTINE DistributePQint(Integral,PQ,QPmat,dimQ,dimP,Input,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,i1,i2,i3,i4,i5
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
integer :: iPassP,iPassQ
integer :: PASSPOFFSET,PASSQOFFSET
!
nAngmomP = PQ%P%p%nAngmom

nderivP=Input%nDerivP
totOrbP = PQ%P%p%totOrbitals
!
iOrbitalP = 1
DO iDerivP=1,NDerivP
   DO iAngmomP=1,nAngmomP
      DO iPassP = 1, PQ%P%p%nPasses
      nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
      totOrbQ = PQ%Q%p%totOrbitals
      endOrbP = iOrbitalP+nOrbP-1
!
      iA = PQ%P%p%indexAng1(iAngmomP)
      iB = PQ%P%p%indexAng2(iAngmomP)
      nContA = PQ%P%p%orbital1%nContracted(iA)
      nContB = PQ%P%p%orbital2%nContracted(iB)
      nAngA = PQ%P%p%orbital1%nOrbComp(iA)
      nAngB = PQ%P%p%orbital2%nOrbComp(iB)
      PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
      startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
      startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)     
      CALL DistributePQin1(QPmat,dimQ,dimP,iOrbitalP,&
           & totOrbQ,totOrbP,nOrbP,nContA,nContB,nAngA,nAngB,startA,startB, &
           & PQ%Q%p,Input,Output,iAngmomP,iDerivP,NderivP,PASSPOFFSET,&
           & PQ%samePQ,LUPRI,IPRINT)
   ENDDO
      iOrbitalP = iOrbitalP + nOrbP
   ENDDO
ENDDO
IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(1X,A,5I5)') 'ResultMat in DistributePQint',Output%ndim(1),Output%ndim(2),&
&Output%ndim(3),Output%ndim(4),Output%ndim(5)
  DO i5=1,Output%ndim(5)
    DO i4=1,Output%ndim(4)
      DO i3=1,Output%ndim(3)
        DO i2=1,Output%ndim(2)
          WRITE(LUPRI,'(3X,4I5)') i2,i3,i4,i5
          WRITE(LUPRI,'(5X,5F12.8)') (Output%resultMat(i1,i2,i3,i4,i5),i1=1,Output%ndim(1))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE DistributePQint

SUBROUTINE DistributePQin1(QPmat,dimQ,dimP,iOrbitalP,totOrbQ,totOrbP,nOrbP,nContA,nContB,nAngA,nAngB,startA,startB,&
           & Q,Input,Int_Output,iAngmomP,iDerivP,NderivP,PASSPOFFSET,samePQ,LUPRI,IPRINT) 
implicit none
integer              :: totOrbQ,nOrbP,nContA,nContB,nAngA,nAngB,startA,startB
integer              :: iOrbitalP,totOrbP,iAngmomP
integer              :: iderivP,NderivP,lupri,iprint
Real(realk)          :: QP(totOrbQ,totOrbP)!QP(totOrbQ,nOrbP)! = QP(totOrbQ,nAngA,nAngB,nContA*nContB)
Type(Overlap)        :: Q
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Int_Output
Logical              :: samePQ
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer                 :: iA,iContA,iAngA,iAngB,iB,iContB,iContP,iOrbP
Integer                 :: nOrbQ,nderivQ,endOrbQ
Integer                 :: iAngmomQ,iOrbitalQ,iderivQ
Integer                 :: endAngmomQ
Integer                 :: iC,iD,nContC,nContD,nAngC,nAngD,startC,startD
integer :: iPassP,iPassQ
integer :: PASSPOFFSET,PASSQOFFSET

iOrbP=iOrbitalP-1+PASSPOFFSET
iContP=0
DO iContB=1,nContB
   DO iContA=1,nContA
      iContP=iContP+1
      DO iAngB=1,nAngB
         iB=startB-1+iAngB+(iContB-1)*nAngB
         DO iAngA=1,nAngA
            iA=startA-1+iAngA+(iContA-1)*nAngA
            iOrbP=iOrbP+1

            iOrbitalQ = 1
            endAngmomQ = Q%nAngmom
            IF (samePQ) endAngmomQ = iAngmomP
            nderivQ=input%nderivQ
            DO iderivQ =1,nderivQ
               DO iAngmomQ=1,endAngmomQ
               DO iPassQ = 1, Q%nPasses
                  nOrbQ = Q%nOrbitals(iAngmomQ)
                  endOrbQ = iOrbitalQ+nOrbQ-1
                  
                  iC = Q%indexAng1(iAngmomQ)
                  iD = Q%indexAng2(iAngmomQ)
                  nContC = Q%orbital1%nContracted(iC)
                  nContD = Q%orbital2%nContracted(iD)
                  nAngC = Q%orbital1%nOrbComp(iC)
                  nAngD = Q%orbital2%nOrbComp(iD)
                  PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngC*nAngD
                  startC  = Q%orbital1%startOrbital(iC,iPassQ)
                  startD  = Q%orbital2%startOrbital(iD,iPassQ)
                  
                  CALL DistributePQin2(QPmat,dimQ,dimP&
                       &      ,iOrbitalQ,iOrbitalP,totOrbQ,totOrbP,iOrbP,&
                       &      nOrbQ,nContC,nContD,nAngC,nAngD,startC,startD,iderivP,&
                       &      iderivQ,NderivP,NderivQ,PASSQOFFSET,int_Output,INPUT,&
                       &      startA,startB,iAngA,iAngB,iContA,iContB,iContP,iA,iB,nAngA,nAngB)
               ENDDO
               iOrbitalQ = iOrbitalQ + nOrbQ
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE DISTRIBUTEPQIN1

SUBROUTINE DistributePQin2(QPmat,dimQ,dimP,iOrbitalQ,iOrbitalP,totOrbQ,totOrbP,iOrbP,nOrbQ,&
     &nContC,nContD,nAngC,nAngD,startC,startD,&
     &iderivP,iderivQ,NderivP,NderivQ,PASSQOFFSET,Output,INPUT,startA,startB,iAngA,&
     &iAngB,iContA,iContB,iContP,iA,iB,nAngA,nAngB)
implicit none
Integer              :: nContC,nContD,nAngC,nAngD,startC,startD,iderivP,iderivQ
Integer              :: nderivP,nderivQ,nOrbQ,nangA,nangB,iOrbitalQ,iOrbitalP
integer              :: totOrbP,iOrbP,totOrbQ
Integer              :: iContA,iContB,iAngA,iAngB,startA,startB,iContP,iA,iB
Real(realk)          :: QP(totOrbQ,totOrbP)
Type(IntegralOutput) :: Output
Type(IntegralInput)  :: Input
Integer              :: dimQ,dimP
REAL(REALK)          :: QPMAT(dimQ,dimP)
!
Integer     :: iC,iD,iContC,iContD,iAngC,iAngD,iContQ,iOrbQ
Real(realk) :: integral,dtemp_lhs,dtemp_rhs
Integer     :: idmat,ideriv
Integer     :: n1,n2
integer :: iPassP,iPassQ,PASSQOFFSET
!
ideriv=iderivP+(iderivQ-1)*NderivP
!     Regular case
iOrbQ=iOrbitalQ-1+PASSQOFFSET
DO iContD=1,nContD
   DO iContC=1,nContC
      iContQ = (iContD-1)*nContC + iContC
      DO iAngD=1,nAngD
         iD=startD-1+iAngD+(iContD-1)*nAngD
         DO iAngC=1,nAngC
            iC=startC-1+iAngC+(iContC-1)*nAngC
            iOrbQ=iOrbQ+1
            integral = QPmat(iOrbQ,iOrbP)
            !        Contract with density
!            IF (Input%DO_FOCK) THEN
!               DO idmat=1,Input%NDMAT_RHS
!                  !          idmat = 1
!                  IF (Input%DO_Coulomb) THEN
!                     !simen change Dmat dimensions
!                     n1 = Output%ndim(1)
!                     n2 = Output%ndim(2)
!                     CALL AddToCoulombMat(Output%ResultMat(:,:,1,1,idmat),n1,n2,&
!                          &                            Input%DO_JENGINE,INPUT%sameLHSaos,INPUT%sameRHSaos,&
!                          &                            INPUT%sameODs,INPUT%CoulombFactor,integral,iA,iB,iC,iD,&
!                          &                            startA,startB,startC,startD,&
!                          &                            Input%DMAT_RHS(:,:,idmat),n1,n2)
!                  ENDIF
!                  IF (Input%DO_Exchange) THEN
!                     !simen change Dmat dimensions
!                     n1 = Output%ndim(1)
!                     n2 = Output%ndim(2)
!                     CALL AddToExchangeMat(Output%ResultMat(:,:,1,1,idmat),n1,n2,&
!                          &                            INPUT%sameLHSaos,INPUT%sameRHSaos,INPUT%sameODs,&
!                          &                            INPUT%exchangeFactor,integral,iA,iB,iC,iD,startA,&
!                          &                            startB,startC,startD,Input%DMAT_RHS(:,:,idmat),n1,n2)
!                  ENDIF
!               ENDDO
!               !        No contraction with density
!            ELSE
               IF(INPUT%AddToIntegral)THEN  
                  IF (INPUT%sameLHSaos.AND.INPUT%sameRHSaos.AND.INPUT%sameODs)THEN
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = &
                          &Output%ResultMat(iA,iB,iC,iD,ideriv)+integral
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = &
                          &Output%ResultMat(iA,iB,iC,iD,ideriv)+integral
                     Output%ResultMat(iB,iA,iC,iD,ideriv) = &
                          &Output%ResultMat(iB,iA,iC,iD,ideriv)+integral
                     Output%ResultMat(iA,iB,iD,iC,ideriv) = &
                          &Output%ResultMat(iA,iB,iD,iC,ideriv)+integral
                     Output%ResultMat(iB,iA,iD,iC,ideriv) = &
                          &Output%ResultMat(iB,iA,iD,iC,ideriv)+integral
                     Output%ResultMat(iC,iD,iA,iB,ideriv) = &
                          &Output%ResultMat(iC,iD,iA,iB,ideriv)+integral
                     Output%ResultMat(iD,iC,iA,iB,ideriv) = &
                          &Output%ResultMat(iD,iC,iA,iB,ideriv)+integral
                     Output%ResultMat(iC,iD,iB,iA,ideriv) = &
                          &Output%ResultMat(iC,iD,iB,iA,ideriv)+integral
                     Output%ResultMat(iD,iC,iB,iA,ideriv) = &
                          &Output%ResultMat(iD,iC,iB,iA,ideriv)+integral
                  ELSE IF (INPUT%sameLHSaos.AND.INPUT%sameRHSaos)THEN
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = &
                          &Output%ResultMat(iA,iB,iC,iD,ideriv)+integral
                     Output%ResultMat(iB,iA,iC,iD,ideriv) = &
                          &Output%ResultMat(iB,iA,iC,iD,ideriv)+integral
                     Output%ResultMat(iA,iB,iD,iC,ideriv) = &
                          &Output%ResultMat(iA,iB,iD,iC,ideriv)+integral
                     Output%ResultMat(iB,iA,iD,iC,ideriv) = &
                          &Output%ResultMat(iB,iA,iD,iC,ideriv)+integral
                  ELSE IF (INPUT%sameLHSaos)THEN
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = &
                          &Output%ResultMat(iA,iB,iC,iD,ideriv)+integral
                     Output%ResultMat(iB,iA,iC,iD,ideriv) = &
                          &Output%ResultMat(iB,iA,iC,iD,ideriv)+integral
                  ELSE IF (INPUT%sameRHSaos)THEN
                     Output%ResultMat(iA,iB,iD,iC,ideriv) = &
                          &Output%ResultMat(iA,iB,iD,iC,ideriv)+integral
                     Output%ResultMat(iB,iA,iD,iC,ideriv) = &
                          &Output%ResultMat(iB,iA,iD,iC,ideriv)+integral
                  ELSE IF (INPUT%sameODs .AND. (INPUT%sameLHSaos.NEQV.INPUT%sameRHSaos))THEN
                     CALL QUIT('Progamming error, sameODs but sameLHSaos.NE.sameRHSaos!')
                  ELSE IF (INPUT%sameODs)THEN
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = &
                          &Output%ResultMat(iA,iB,iC,iD,ideriv)+integral
                     Output%ResultMat(iC,iD,iA,iB,ideriv) = &
                          &Output%ResultMat(iC,iD,iA,iB,ideriv)+integral
                  ELSE
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = &
                          &Output%ResultMat(iA,iB,iC,iD,ideriv)+integral
                  ENDIF
               ELSE
                  IF (INPUT%sameLHSaos.AND.INPUT%sameRHSaos.AND.INPUT%sameODs)THEN
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = integral
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = integral
                     Output%ResultMat(iB,iA,iC,iD,ideriv) = integral
                     Output%ResultMat(iA,iB,iD,iC,ideriv) = integral
                     Output%ResultMat(iB,iA,iD,iC,ideriv) = integral
                     Output%ResultMat(iC,iD,iA,iB,ideriv) = integral
                     Output%ResultMat(iD,iC,iA,iB,ideriv) = integral
                     Output%ResultMat(iC,iD,iB,iA,ideriv) = integral
                     Output%ResultMat(iD,iC,iB,iA,ideriv) = integral
                  ELSE IF (INPUT%sameLHSaos.AND.INPUT%sameRHSaos)THEN
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = integral
                     Output%ResultMat(iB,iA,iC,iD,ideriv) = integral
                     Output%ResultMat(iA,iB,iD,iC,ideriv) = integral
                     Output%ResultMat(iB,iA,iD,iC,ideriv) = integral
                  ELSE IF (INPUT%sameLHSaos)THEN
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = integral
                     Output%ResultMat(iB,iA,iC,iD,ideriv) = integral
                  ELSE IF (INPUT%sameRHSaos)THEN
                     Output%ResultMat(iA,iB,iD,iC,ideriv) = integral
                     Output%ResultMat(iB,iA,iD,iC,ideriv) = integral
                  ELSE IF (INPUT%sameODs .AND. (INPUT%sameLHSaos.NEQV.INPUT%sameRHSaos))THEN
                     CALL QUIT('Progamming error, sameODs but sameLHSaos.NE.sameRHSaos!')
                  ELSE IF (INPUT%sameODs)THEN
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = integral
                     Output%ResultMat(iC,iD,iA,iB,ideriv) = integral
                  ELSE
                     Output%ResultMat(iA,iB,iC,iD,ideriv) = integral
                  ENDIF
               ENDIF
            !ENDIF
            !               iD = iD+1
         ENDDO
      ENDDO
      !         iC = iC+1
   ENDDO
ENDDO


END SUBROUTINE DistributePQin2

SUBROUTINE SET_OVERLAP(P,np,Input,sharedTUV,Integral,Alloc,ODB,IELECTRON,LUPRI,IPRINT,SIDE)
Implicit none
TYPE(IntegralInput) :: Input
TYPE(IntegralItem)  :: Integral
TYPE(AllocItem)     :: Alloc
TYPE(Overlap)       :: P
TYPE(ODBATCH)       :: ODB
Character*(*)       :: side
Integer             :: IELECTRON,LUPRI,IPRINT,nderiv,np
TYPE(TUVitem)       :: SharedTUV
!
Integer             :: idir,i1,i2,i12,start1,start2,end1,end2,a,b,a1,a2,orb1,orb2
integer             :: l1,l2,ijk1,ijk2,ijk,maxijk
Real(realk)         :: e1,e2,d2,maxGab,maxPrimGab,maxPrimGabElm,signP
Character(len=80)   :: t1,t2
Real(realk),pointer :: GAB(:,:),primGAB(:,:) 
LOGICAL             :: LHS 
!
SELECT CASE(SIDE)
CASE('LHS')
   IF(INPUT%PS_SCREEN)THEN
      primGAB => INPUT%pGAB_LHS
      maxPrimGabElm=INPUT%PS_MAXELM_RHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
      GAB => INPUT%GAB_LHS
   ENDIF
   LHS=.TRUE.
   NDERIV=INPUT%NDERIVP
CASE('RHS')
   IF(INPUT%PS_SCREEN)THEN
      primGAB => INPUT%pGAB_RHS
      maxPrimGabElm=INPUT%PS_MAXELM_LHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
      GAB => INPUT%GAB_RHS
   ENDIF
   LHS=.FALSE.
   NDERIV=INPUT%NDERIVQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_OVERLAP =',SIDE
   CALL QUIT('Wrong case in SET_OVERLAP')
END SELECT

 CALL SET_ORBITAL(P%orbital1,ODB%AO(1)%p,integral,LUPRI)
 CALL SET_ORBITAL(P%orbital2,ODB%AO(2)%p,integral,LUPRI)

!Determine overlap type !I want to remove type
 t1 = P%orbital1%type
 t2 = P%orbital2%type
 P%hermiteSingle = .FALSE.
 IF ( t1.EQ.'Empty'.AND.t2.EQ.'Empty') THEN
   P%type = 'Empty'
 ELSE IF ( (t1.EQ.'Hermite'.AND.t2.EQ.'Empty').OR.(t1.EQ.'Empty'  .AND.t2.EQ.'Hermite') ) THEN
   P%type = 'Hermite-single'
   P%hermiteSingle = .TRUE.
 ELSE IF (t1.EQ.'Hermite'.AND.t2.EQ.'Hermite') THEN
   P%type = 'Hermite'
 ELSE IF ( (t1.EQ.'Cartesian'.AND.t2.EQ.'Empty').OR.(t1.EQ.'Empty'  .AND.t2.EQ.'Cartesian') ) THEN
   P%type = 'Cartesian-singel'
 ELSE IF (t1.EQ.'Cartesian'.AND.t2.EQ.'Cartesian') THEN
   P%type = 'Cartesian'
 ELSE
   WRITE(LUPRI,'(1X,A)') 'Not a proper combination of orbital types in SET_OVERLAP:'
   WRITE(LUPRI,'(5X,4A)') 'orbital1%type =', t1, 'and orbital2%type =', t2
   CALL QUIT('Not a proper combination of orbital types in SET_OVERLAP')
 ENDIF

 P%sameAO = ODB%sameAO
 P%ETUVisSet = .FALSE.
 P%sphericalETUV = .TRUE.
 CALL GET_NPRIMITIVES(np,P,Input,Side)!can be done outside
 P%nPrimitives   = np
! Default is to build OD-batches first, and then later collect into passes
 P%nPasses       = 1
 IF (np.EQ.0) THEN
   CALL FREE_ORBITAL(P%orbital1)
   CALL FREE_ORBITAL(P%orbital2)
   RETURN
 ENDIF

 P%nAngmom       = ODB%nAngmom
 P%maxContracted = ODB%maxContracted
!
 Nullify(P%distance12)
 Allocate(P%distance12(3,1))
 d2 = 0.0d0
 DO idir=1,3
   P%distance12(idir,1) = P%orbital1%center(idir)-P%orbital2%center(idir)
   d2 = d2 + P%distance12(idir,1) * P%distance12(idir,1)
 ENDDO
 P%squaredDistance = d2
!
 Nullify(P%center)
 Nullify(P%exponents)
 Nullify(P%reducedExponents)
 Nullify(P%preExpFac)
 Nullify(P%iprim1)
 Nullify(P%iprim2)
 Allocate(P%center(3,np))
 Allocate(P%exponents(np))
 Allocate(P%reducedExponents(np))
 Allocate(P%preExpFac(np))
 Allocate(P%iprim1(np))
 Allocate(P%iprim2(np))


 i12 = 0
 DO i1=1,P%orbital1%nPrimitives
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nPrimitives
     maxPrimGab=0.0d0
     IF(INPUT%PS_SCREEN)THEN
        DO a1=1,P%orbital1%nAngmom 
           DO a2=1,P%orbital2%nAngmom
              a = P%orbital1%startprimOrbital(a1,1)+(i1-1)*P%orbital1%nOrbComp(a1)
              b = P%orbital2%startprimOrbital(a2,1)+(i2-1)*P%orbital2%nOrbComp(a2)
              maxPrimGab=MAX(maxPrimGab,primGAB(a,b))
           ENDDO
        ENDDO
     ENDIF
     IF(INPUT%PS_SCREEN)THEN
       IF(maxPrimGab .GT. INPUT%PS_THRESHOLD/maxPrimGabElm)THEN
        i12 = i12 + 1
        e1  = P%orbital1%exponents(i1)
        e2  = P%orbital2%exponents(i2)
        P%exponents(i12)        = e1 + e2
        IF (P%exponents(i12) .NE. 0.0d0) THEN
          P%reducedExponents(i12) = e1*e2/(e1+e2)
          P%center(1,i12)         = (e1*P%orbital1%center(1) + e2*P%orbital2%center(1))/(e1+e2)
          P%center(2,i12)         = (e1*P%orbital1%center(2) + e2*P%orbital2%center(2))/(e1+e2)
          P%center(3,i12)         = (e1*P%orbital1%center(3) + e2*P%orbital2%center(3))/(e1+e2)
          P%preExpFac(i12)        = exp(-P%reducedExponents(i12)*d2)!fast way, maybe afterwards - P%reducedExponents(i12) sequential in mem
        ELSEIF(Input%operator(1:6) .EQ. 'Nucrep')THEN
          P%reducedExponents(i12) = e1
          P%center(1,i12)         = P%orbital1%center(1)
          P%center(2,i12)         = P%orbital1%center(2)
          P%center(3,i12)         = P%orbital1%center(3)
          P%preExpFac(i12)        = exp(-e1*d2)
        ELSE
          P%reducedExponents(i12) = 0.0d0
          P%center(1,i12)         = 0.0d0
          P%center(2,i12)         = 0.0d0
          P%center(3,i12)         = 0.0d0
          P%preExpFac(i12)        = 1.0d0
        ENDIF
        IF ( t1.EQ.'Empty'.OR.t2.EQ.'Empty') THEN 
          P%preExpFac(i12)        = 1.0d0
        ELSE
          IF (P%exponents(i12) .NE. 0.0d0) THEN
             P%preExpFac(i12)        = exp(-P%reducedExponents(i12)*d2)
          ELSEIF(Input%operator(1:6) .EQ. 'Nucrep')THEN
             P%preExpFac(i12)        = exp(-e1*d2)
          ELSE
             P%preExpFac(i12)        = 1.0d0
          ENDIF
        ENDIF
        P%iprim1(i12)           = i1
        P%iprim2(i12)           = i2
       ENDIF
     ELSE
        i12 = i12 + 1
        e1  = P%orbital1%exponents(i1)
        e2  = P%orbital2%exponents(i2)
        P%exponents(i12)        = e1 + e2
        IF (P%exponents(i12) .NE. 0.0d0) THEN
          P%reducedExponents(i12) = e1*e2/(e1+e2)
          P%center(1,i12)         = (e1*P%orbital1%center(1) + e2*P%orbital2%center(1))/(e1+e2)
          P%center(2,i12)         = (e1*P%orbital1%center(2) + e2*P%orbital2%center(2))/(e1+e2)
          P%center(3,i12)         = (e1*P%orbital1%center(3) + e2*P%orbital2%center(3))/(e1+e2)
          P%preExpFac(i12)        = exp(-P%reducedExponents(i12)*d2)!fast way, maybe afterwards - P%reducedExponents(i12) sequential in mem
        ELSEIF(Input%operator(1:6) .EQ. 'Nucrep')THEN
          P%reducedExponents(i12) = e1
          P%center(1,i12)         = P%orbital1%center(1)
          P%center(2,i12)         = P%orbital1%center(2)
          P%center(3,i12)         = P%orbital1%center(3)
          P%preExpFac(i12)        = exp(-e1*d2)
        ELSE
          P%reducedExponents(i12) = 0.0d0
          P%center(1,i12)         = 0.0d0
          P%center(2,i12)         = 0.0d0
          P%center(3,i12)         = 0.0d0
          P%preExpFac(i12)        = 1.0d0
        ENDIF
        IF ( t1.EQ.'Empty'.OR.t2.EQ.'Empty') THEN 
          P%preExpFac(i12)        = 1.0d0
        ELSE
          IF (P%exponents(i12) .NE. 0.0d0) THEN
             P%preExpFac(i12)        = exp(-P%reducedExponents(i12)*d2) !are you doing this twice?
          ELSEIF(Input%operator(1:6) .EQ. 'Nucrep')THEN
             P%preExpFac(i12)        = exp(-e1*d2)
          ELSE
             P%preExpFac(i12)        = 1.0d0
          ENDIF
        ENDIF
        P%iprim1(i12)           = i1
        P%iprim2(i12)           = i2
     ENDIF   
   ENDDO
 ENDDO
 IF (np .NE. i12) THEN
   WRITE(LUPRI,'(1X,A,2I5)') 'Error in set_overlap. Mismatch in number of primitives, (np,i12)=',&
     &  np,i12
   CALL QUIT('Error: Mismatch between get_nPrimitives and i12 in set_overlap')
 ENDIF

IF(INPUT%PS_SCREEN)THEN
   P%maxPrimGab = maxPrimGab
   NULLIFY(primGAB)
ENDIF

 Nullify(P%angmom)
 Allocate(P%angmom(P%nAngmom))
 Nullify(P%indexAng1)
 Allocate(P%indexAng1(P%nAngmom))
 Nullify(P%indexAng2)
 Allocate(P%indexAng2(P%nAngmom))
 Nullify(P%nOrbitals)
 Allocate(P%nOrbitals(P%nAngmom))
 Nullify(P%nOrbComp)
 Allocate(P%nOrbComp(P%nAngmom))
 Nullify(P%nContracted)
 Allocate(P%nContracted(P%nAngmom))

!CAREFUL
!Initialize to twice some maximal value + 1
 P%minAngmom   = 99
!CAREFUL
 P%maxAngmom = 0
 P%totOrbitals = 0
 i12 = 0
 maxGab = 0.0d0
 maxijk = 0
 DO i1=1,P%orbital1%nAngmom
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nAngmom
     i12  = i12 + 1
     l1   = P%orbital1%angmom(i1)
     l2   = P%orbital2%angmom(i2)
     ijk1 = (l1+1)*(l1+2)/2
     ijk2 = (l2+1)*(l2+2)/2
     ijk  = ijk1*ijk2
     maxijk = max(maxijk,ijk)
     P%angmom(i12)      = P%orbital1%angmom(i1) + P%orbital2%angmom(i2)
     P%indexAng1(i12)   = i1
     P%indexAng2(i12)   = i2
     P%nOrbitals(i12)   = P%orbital1%nOrbitals(i1)*P%orbital2%nOrbitals(i2)
     P%totOrbitals      = P%totOrbitals + P%nOrbitals(i12)*NDERIV
     P%nOrbComp(i12)    = P%orbital1%nOrbComp(i1)*P%orbital2%nOrbComp(i2)
     P%nContracted(i12) = P%orbital1%nContracted(i1)*P%orbital2%nContracted(i2)
     P%minAngmom        = min(P%minAngmom,P%angmom(i12))
     P%maxAngmom        = max(P%maxAngmom,P%angmom(i12))
     start1 = P%orbital1%startOrbital(i1,1)
     end1   = start1 + P%orbital1%nOrbitals(i1) - 1
     start2 = P%orbital2%startOrbital(i2,1)
     end2   = start2 + P%orbital2%nOrbitals(i2) - 1
     IF (Input%CS_SCREEN) THEN
       DO orb1=start1,end1
         DO orb2=start2,end2
           maxGab = MAX(maxGab,GAB(orb1,orb2))
         ENDDO
       ENDDO
     ENDIF
   ENDDO
 ENDDO
 IF (Input%CS_SCREEN) THEN
   P%maxGab = maxGab
   NULLIFY(GAB)
 ENDIF

 P%startAngmom = 0
 IF (P%hermiteSingle) P%startAngmom = P%minAngmom

 IF( INPUT%operator(1:7) .EQ. 'Kinetic' .AND. LHS) THEN
    P%startAngmom = P%startAngmom + 2
    P%endAngmom   = P%maxAngmom + 2
 ELSE
    P%endAngmom   = P%maxAngmom
 ENDIF  

 P%nTUV = (P%endAngmom+1)*(P%endAngmom+2)*(P%endAngmom+3)/6 - &
    &     P%startAngmom*(P%startAngmom+1)*(P%startAngmom+2)/6

 IF (IELECTRON.EQ.1) THEN
   signP=1.0d0
   Alloc%maxPrimLHS    = max(np,Alloc%maxPrimLHS)
   Alloc%maxPrimTUVLHS = max(np*P%nTUV,P%totOrbitals,np*maxijk,Alloc%maxPrimTUVLHS)
   Alloc%maxnAngLHS     = max(P%nAngmom,Alloc%maxnAngLHS)
   Alloc%maxContLHS    = max(P%Orbital1%maxContracted,P%Orbital2%maxContracted,Alloc%maxContLHS)
   Alloc%maxangmomLHS = max(P%Orbital1%maxAngmom,P%Orbital2%maxAngmom,Alloc%maxangmomLHS)
   Alloc%maxTUVLHS = max(P%nTUV,Alloc%maxTUVLHS)
   Alloc%maxTotOrbLHS = max(P%totOrbitals,Alloc%maxTotOrbLHS)
   Alloc%maxijkLHS = max(maxijk,Alloc%maxijkLHS)
 ELSE
   signP=-1.0d0
   Alloc%maxPrimRHS    = max(np,Alloc%maxPrimRHS)
   Alloc%maxPrimTUVRHS = max(np*P%nTUV,P%totOrbitals,np*maxijk,Alloc%maxPrimTUVRHS)
   Alloc%maxnAngRHS     = max(P%nAngmom,Alloc%maxnAngRHS)
   Alloc%maxContRHS    = max(P%Orbital1%maxContracted,P%Orbital2%maxContracted,Alloc%maxContRHS)
   Alloc%maxangmomRHS = max(P%Orbital1%maxAngmom,P%Orbital2%maxAngmom,Alloc%maxangmomRHS)
   Alloc%maxTUVRHS = max(P%nTUV,Alloc%maxTUVRHS)
   Alloc%maxTotOrbRHS = max(P%totOrbitals,Alloc%maxTotOrbRHS)
   Alloc%maxijkRHS = max(maxijk,Alloc%maxijkRHS)
 ENDIF

IF (INPUT%DO_JENGINE .AND. IELECTRON.EQ.2) THEN
  CALL SET_FTUV(P,INPUT,Integral,LUPRI,IPRINT)
! Change Overlap P to confirm with FTUV
  P%TYPE = 'FTUV'
  P%nAngmom = 1
  P%orbital1%nAngmom = 1
  P%orbital2%nAngmom = 1
  P%angmom(1) = P%maxAngmom
  P%orbital1%angmom(1) = 0
  P%orbital2%angmom(1) = 0
  P%orbital1%nContracted(1)  = 1
  P%orbital2%nContracted(1)  = 1
  P%orbital1%nOrbComp(1)     = 1
  P%orbital2%nOrbComp(1)     = 1
  P%orbital1%startOrbital(1,1) = 1
  P%orbital2%startOrbital(1,1) = 1
  P%indexAng1(1)   = 1
  P%indexAng2(1)   = 1
  P%nOrbitals(1)   = 1
  P%totOrbitals    = 1
  P%nOrbComp(1)    = 1
  P%nContracted(1) = 1
ELSE
!GET E COEFFICIENTS
 IF (.NOT.P%hermiteSingle) THEN
   IF (INPUT%setETUVoutside) THEN
    CALL AttachETUVtoOverlap(P,signP,SharedTUV,LUPRI,IPRINT)
    IF(IELECTRON.EQ.2)THEN
       Alloc%maxPrimTUVijkRHS = max(P%lenETUV,Alloc%maxPrimTUVijkRHS) 
    ELSE
       Alloc%maxPrimTUVijkLHS = max(P%lenETUV,Alloc%maxPrimTUVijkLHS) 
    ENDIF
   ENDIF
 ENDIF
ENDIF

IF(INPUT%CS_SCREEN)THEN
   NULLIFY(GAB)
ENDIF
END SUBROUTINE SET_OVERLAP

SUBROUTINE AddOverlapToPass(PassP,P,numPass,maxPasses,LUPRI,IPRINT)
Implicit none
TYPE(Overlap)       :: PassP,P
Integer             :: numPass,maxPasses,LUPRI,IPRINT
!
Integer :: np,na,mc,ma,mp,a1,a2,ia
Integer :: iAngmom,l1,l2,l,ijk1,ijk2,ijk,nTUV,length,indexETUV,passIndexETUV
!
integer :: i

!IF (.not. P%ETUVisSet) CALL QUIT('ETUV not set in AddOverlapToPass')
IF (P%TYPE.EQ.'FTUV') CALL QUIT('Cannot call AddOverlapToPass for FTUV-type overlap')

! Initial allocations and setting
IF (numPass.EQ.1) THEN

   ! Orbital 1
   np = P%orbital1%nPrimitives
   na = P%orbital1%nAngmom
   mc = P%orbital1%maxContracted
   ma = P%orbital1%maxAngmom

!   NULLIFY(PassP%Orbital1%exponents)
!   ALLOCATE(PassP%Orbital1%exponents(np))
!!$   NULLIFY(PassP%Orbital1%angmom)
!!$   ALLOCATE(PassP%Orbital1%angmom(na))
!!$   NULLIFY(PassP%Orbital1%nContracted)
!!$   ALLOCATE(PassP%Orbital1%nContracted(na))
!!$   NULLIFY(PassP%Orbital1%startOrbital)
!!$   ALLOCATE(PassP%Orbital1%startOrbital(na,maxpasses))
!!$   NULLIFY(PassP%Orbital1%startprimOrbital)
!!$   ALLOCATE(PassP%Orbital1%startprimOrbital(na,maxpasses))
!!$   NULLIFY(PassP%Orbital1%nOrbComp)
!!$   ALLOCATE(PassP%Orbital1%nOrbComp(na))
!!$   NULLIFY(PassP%Orbital1%nPrimOrbComp)
!!$   ALLOCATE(PassP%Orbital1%nPrimOrbComp(na))
!!$   NULLIFY(PassP%Orbital1%nOrbitals)
!!$   ALLOCATE(PassP%Orbital1%nOrbitals(na))
!   NULLIFY(PassP%Orbital1%CC)
!   ALLOCATE(PassP%Orbital1%CC(np,mc,na))
!   NULLIFY(PassP%Orbital1%SPH_MAT)
!   ALLOCATE(PassP%Orbital1%SPH_MAT(ma+1))

!   NULLIFY(PassP%orbital1%startOrbital)
!   NULLIFY(PassP%orbital1%startprimOrbital)
!   ALLOCATE(PassP%orbital1%startOrbital(na,maxPasses))
!   ALLOCATE(PassP%orbital1%startprimOrbital(na,maxPasses))
   
   PassP%orbital1%type          = P%orbital1%type          
   PassP%orbital1%spherical     = P%orbital1%spherical     
   PassP%orbital1%center        = P%orbital1%center        
   PassP%orbital1%nPrimitives   = np
   PassP%orbital1%maxContracted = mc
   PassP%orbital1%maxAngmom     = ma
   PassP%orbital1%nAngmom       = na
   PassP%orbital1%CCidentifier  = P%orbital1%CCidentifier  
   !replace with dcopy
   PassP%orbital1%exponents(1:np) = P%orbital1%exponents(1:np)
   PassP%orbital1%angmom(1:na)       = P%orbital1%angmom(1:na)
   PassP%orbital1%nContracted(1:na)  = P%orbital1%nContracted(1:na)
   PassP%orbital1%nOrbComp(1:na)     = P%orbital1%nOrbComp(1:na)
   PassP%orbital1%nPrimOrbComp(1:na) = P%orbital1%nPrimOrbComp(1:na)
   PassP%orbital1%nOrbitals(1:na)    = P%orbital1%nOrbitals(1:na)
   !maybe replace with do loop
!   PassP%orbital1%CC(1:np,1:mc,1:na) = P%orbital1%CC(1:np,1:mc,1:na)
!   NULLIFY(PassP%orbital1%CC)
!   ALLOCATE(PassP%orbital1%CC(na))
   DO ia = 1,na
      PassP%orbital1%CC(ia)%p => P%orbital1%CC(ia)%p
   ENDDO
   PassP%orbital1%SPH_MAT(1:ma+1) = P%orbital1%SPH_MAT(1:ma+1)
   
   ! Orbital 2
   np = P%orbital2%nPrimitives
   na = P%orbital2%nAngmom
   mc = P%orbital2%maxContracted
   ma = P%orbital2%maxAngmom
   
!   NULLIFY(PassP%Orbital2%exponents)
!   ALLOCATE(PassP%Orbital2%exponents(np))
!!$   NULLIFY(PassP%Orbital2%angmom)
!!$   ALLOCATE(PassP%Orbital2%angmom(na))
!!$   NULLIFY(PassP%Orbital2%nContracted)
!!$   ALLOCATE(PassP%Orbital2%nContracted(na))
!!$   NULLIFY(PassP%Orbital2%startOrbital)
!!$   ALLOCATE(PassP%Orbital2%startOrbital(na,maxpasses))
!!$   NULLIFY(PassP%Orbital2%startprimOrbital)
!!$   ALLOCATE(PassP%Orbital2%startprimOrbital(na,maxpasses))
!!$   NULLIFY(PassP%Orbital2%nOrbComp)
!!$   ALLOCATE(PassP%Orbital2%nOrbComp(na))
!!$   NULLIFY(PassP%Orbital2%nPrimOrbComp)
!!$   ALLOCATE(PassP%Orbital2%nPrimOrbComp(na))
!!$   NULLIFY(PassP%Orbital2%nOrbitals)
!!$   ALLOCATE(PassP%Orbital2%nOrbitals(na))
!   NULLIFY(PassP%Orbital2%CC)
!   ALLOCATE(PassP%Orbital2%CC(np,mc,na))
!   NULLIFY(PassP%Orbital2%SPH_MAT)
!   ALLOCATE(PassP%Orbital2%SPH_MAT(ma+1))

!   NULLIFY(PassP%orbital2%startOrbital)
!   NULLIFY(PassP%orbital2%startprimOrbital)
!   ALLOCATE(PassP%orbital2%startOrbital(na,maxPasses))
!   ALLOCATE(PassP%orbital2%startprimOrbital(na,maxPasses))

   PassP%orbital2%type          = P%orbital2%type          
   PassP%orbital2%spherical     = P%orbital2%spherical     
   PassP%orbital2%center        = P%orbital2%center        
   PassP%orbital2%nPrimitives   = np
   PassP%orbital2%maxContracted = mc
   PassP%orbital2%maxAngmom     = ma
   PassP%orbital2%nAngmom       = na
   PassP%orbital2%CCidentifier  = P%orbital2%CCidentifier  
   PassP%orbital2%exponents(1:np) = P%orbital2%exponents(1:np)
   PassP%orbital2%angmom(1:na)       = P%orbital2%angmom(1:na)
   PassP%orbital2%nContracted(1:na)  = P%orbital2%nContracted(1:na)
   PassP%orbital2%nOrbComp(1:na)     = P%orbital2%nOrbComp(1:na)
   PassP%orbital2%nPrimOrbComp(1:na) = P%orbital2%nPrimOrbComp(1:na)
   PassP%orbital2%nOrbitals(1:na)    = P%orbital2%nOrbitals(1:na)
!   PassP%orbital2%CC(1:np,1:mc,1:na) = P%orbital2%CC(1:np,1:mc,1:na)
!   NULLIFY(PassP%orbital2%CC)
!   ALLOCATE(PassP%orbital2%CC(na))
   DO ia = 1,na
      PassP%orbital2%CC(ia)%p => P%orbital2%CC(ia)%p
   ENDDO

   PassP%orbital2%SPH_MAT(1:ma+1) = P%orbital2%SPH_MAT(1:ma+1)

   ! Overlap
   mp = P%nPrimitives*maxPasses
   np = P%nPrimitives
   na = P%nAngmom
!   Nullify(PassP%center)
!   Nullify(PassP%exponents)
!   Nullify(PassP%reducedExponents)
!   Nullify(PassP%preExpFac)
!   Nullify(PassP%iprim1)
!   Nullify(PassP%iprim2)
!   Nullify(PassP%nOrbitals)
!   Nullify(PassP%ETUVindex)
!   Allocate(PassP%center(3,mp))
!   Allocate(PassP%exponents(mp))
!   Allocate(PassP%reducedExponents(mp))
!   Allocate(PassP%preExpFac(mp))
!   Allocate(PassP%iprim1(mp))
!   Allocate(PassP%iprim2(mp))
!!$   Allocate(PassP%Angmom(na))
!!$   Allocate(PassP%indexAng1(na))
!!$   Allocate(PassP%indexAng2(na))
!!$   Allocate(PassP%nOrbComp(na))
!!$   Allocate(PassP%nContracted(na))
!!$   Allocate(PassP%nOrbitals(na))
!!$   Allocate(PassP%ETUVindex(na))
   !ToDo sort according to iprim1 and iprim2
   ! Allocate(PassP%iprim1(np))
   ! Allocate(PassP%iprim2(np))

   PassP%hermiteSingle   = P%hermiteSingle
   PassP%type            = P%type
   PassP%sameAO          = P%sameAO
   PassP%ETUVisSet       = P%ETUVisSet
   PassP%sphericalETUV   = P%sphericalETUV
   PassP%nPrimitives     = P%nPrimitives
   PassP%maxContracted   = P%maxContracted
   PassP%squaredDistance = P%squaredDistance
   PassP%maxPrimGab      = P%maxPrimGab 
   PassP%maxGab          = P%maxGab 
   PassP%nAngmom         = P%nAngmom      
   PassP%angmom(1:na)    = P%angmom(1:na)
   PassP%indexAng1(1:na) = P%indexAng1(1:na)
   PassP%indexAng2(1:na) = P%indexAng2(1:na)
   PassP%nOrbComp(1:na)  = P%nOrbComp(1:na)
   PassP%nContracted(1:na) = P%nContracted(1:na)
   PassP%minAngmom       = P%minAngmom
   PassP%maxAngmom       = P%maxAngmom
   PassP%startAngmom     = P%startAngmom
   PassP%endAngmom       = P%endAngmom
   PassP%nTUV            = P%nTUV
!   IF (PassP%ETUVisSet) THEN
!      NULLIFY(PassP%ETUV)
!      ALLOCATE(PassP%ETUV(P%lenETUV*maxPasses))
!   ENDIF
   PassP%distance12(:,1)      = P%distance12(:,1)
ENDIF

! Settings for all passes

! Orbital indexes
na = P%orbital1%nAngmom
PassP%orbital1%startPrimOrbital(1:na,numPass) = P%orbital1%startPrimOrbital(1:na,1)
PassP%orbital1%startOrbital(1:na,numPass)     = P%orbital1%startOrbital(1:na,1)
na = P%orbital2%nAngmom
PassP%orbital2%startPrimOrbital(1:na,numPass) = P%orbital2%startPrimOrbital(1:na,1)
PassP%orbital2%startOrbital(1:na,numPass)     = P%orbital2%startOrbital(1:na,1)
PassP%orbital1%nPasses                     = numPass
PassP%orbital2%nPasses                     = numPass

! Overlap specifics
IF (numPass.EQ.1) THEN
   PassP%maxPrimGab  = P%maxPrimGab
   PassP%maxGab      = P%maxGab
ELSE
   PassP%maxPrimGab  = MAX (PassP%maxPrimGab,P%maxPrimGab)
   PassP%maxGab      = MAX (PassP%maxGab,P%maxGab)
ENDIF

PassP%nPasses     = numPass
PassP%totOrbitals = P%totOrbitals*numPass
na = P%nAngmom
PassP%nOrbitals(1:na)   = P%nOrbitals(1:na)*numPass
np=P%nPrimitives
CALL AddODtoPass(PassP%center,PassP%exponents,PassP%reducedExponents,&
     &           PassP%preExpFac,PassP%iprim1,PassP%iprim2,&
     &           P%center(:,1:np),P%exponents(1:np),P%reducedExponents(1:np),P%preExpFac(1:np),&
     &           P%iprim1(1:np),P%iprim2(1:np),numPass,maxPasses,P%nPrimitives)

PassP%distance12(:,numpass)      = P%distance12(:,1)

IF (PassP%ETUVisSet) THEN
   DO iAngmom=1,P%nAngmom
      ! Set up indeces and dimensions
      indexETUV = P%ETUVindex(iAngmom)
      passIndexETUV = (indexETUV-1)*maxPasses+1
      PassP%ETUVindex(iAngmom) = passIndexETUV
      a1     = P%indexAng1(iAngmom)
      a2     = P%indexAng2(iAngmom)
      l1     = P%orbital1%angmom(a1)
      l2     = P%orbital2%angmom(a2)
      l      = l1+l2
      ijk1   = (l1+1)*(l1+2)/2
      ijk2   = (l2+1)*(l2+2)/2
      ijk    = ijk1*ijk2
      nTUV   = (l+1)*(l+2)*(l+3)/6
      length = nTUV*ijk*P%nPrimitives
      CALL AddETUVtoPass(PassP%ETUV,passIndexETUV,P%ETUV,indexETUV,&
           &                length,numPass,maxPasses)
   ENDDO
ELSE
   PassP%distance12(:,numpass)      = P%distance12(:,1)
ENDIF
!IF(IPRINT.GT.15)THEN
!   CALL PRINT_OVERLAP(PassP,LUPRI,IPRINT)
!   CALL PRINT_OVERLAP(PassP,LUPRI,300)
!ENDIF

END SUBROUTINE AddOverlapToPass

SUBROUTINE AddETUVtoPass(PassETUV,passIndexETUV,ETUV,indexETUV,&
     &                   length,numPass,maxPasses)
implicit none
Integer             :: passIndexETUV,indexETUV,length,numPass,maxPasses
Real(realk),pointer :: PassETUV(:),ETUV(:)
!
Integer :: startPassIndex,lenm1

startPassIndex = passIndexETUV + (numPass-1)*length
lenm1 = length - 1
PassETUV(startPassIndex:startPassIndex + lenm1) = ETUV(indexETUV:indexETUV+lenm1)
 
END SUBROUTINE AddETUVtoPass

SUBROUTINE AddODtoPass(PassPcenter,PassPexponents,PassPreducedExponents,&
     &            PassPpreExpFac,PassPiprim1,PassPiprim2,&
     &            Pcenter,Pexponents,PreducedExponents,PpreExpFac,&
     &            Piprim1,Piprim2,numPass,maxPasses,nPrim)
implicit none
Integer             :: numPass,maxPasses,nPrim
Integer             :: PassPiprim1(:),PassPiprim2(:)
Integer             :: Piprim1(nPrim),Piprim2(nPrim)
Real(realk)         :: PassPcenter(:,:),PassPexponents(:)
Real(realk)         :: PassPreducedExponents(:),PassPpreExpFac(:)
Real(realk)         :: Pcenter(3,nPrim),Pexponents(nPrim)
Real(realk)         :: PreducedExponents(nPrim),PpreExpFac(nPrim)
!
Integer :: startPrim,endOrbital
!
startPrim  = 1+(numPass-1)*nPrim
endOrbital = numPass*nPrim
!
PassPiPrim1(startPrim:endOrbital)           = Piprim1
PassPiPrim2(startPrim:endOrbital)           = Piprim2
PassPcenter(1:3,startPrim:endOrbital)       = Pcenter
PassPexponents(startPrim:endOrbital)        = Pexponents
PassPreducedExponents(startPrim:endOrbital) = PreducedExponents
PassPpreExpFac(startPrim:endOrbital)        = PpreExpFac

END SUBROUTINE AddODtoPass

SUBROUTINE AttachETUVtoOverlap(P,signP,TUV,LUPRI,IPRINT)
Implicit none
TYPE(Overlap) :: P
TYPE(TUVitem) :: TUV
Integer       :: LUPRI,IPRINT
Real(realk)   :: signP
!
Integer :: iA1,iA2,l1,l2,l,ijk1,ijk2,ijk,nTUV,indexETUV,nETUV,iAngmom
integer :: LENGTH

NULLIFY(P%ETUVindex)
ALLOCATE(P%ETUVindex(P%nAngmom))
indexETUV = 1
DO iAngmom=1,P%nAngmom
  P%ETUVindex(iAngmom) = indexETUV
  iA1 = P%indexAng1(iangmom)
  iA2 = P%indexAng2(iangmom)
  l1 = P%orbital1%angmom(iA1)
  l2 = P%orbital2%angmom(iA2)
  l  = l1+l2
  ijk1 = (l1+1)*(l1+2)/2
  ijk2 = (l2+1)*(l2+2)/2
  ijk  = ijk1*ijk2
  nTUV = (l+1)*(l+2)*(l+3)/6
  indexETUV = indexETUV + nTUV*ijk*P%nPrimitives
ENDDO
nETUV = indexETUV - 1
P%lenETUV = nETUV

NULLIFY(P%ETUV)
ALLOCATE(P%ETUV(nETUV))

DO iAngmom=1,P%nAngmom
  indexETUV = P%ETUVindex(iAngmom)
  iA1 = P%indexAng1(iangmom)
  iA2 = P%indexAng2(iangmom)
  l1 = P%orbital1%angmom(iA1)
  l2 = P%orbital2%angmom(iA2)
  l  = l1+l2
  ijk1 = (l1+1)*(l1+2)/2
  ijk2 = (l2+1)*(l2+2)/2
  ijk  = ijk1*ijk2
  nTUV = (l+1)*(l+2)*(l+3)/6
  LENGTH=ijk*nTUV*P%nPrimitives
  CALL BuildEcoeffTensor2(TUV,P,signP,P%ETUV(indexETUV:indexETUV+LENGTH-1),ijk,nTUV,P%nPrimitives,iAngmom,1,LUPRI,IPRINT)
ENDDO

P%ETUVisSet = .TRUE.

END SUBROUTINE AttachETUVtoOverlap

SUBROUTINE SET_FTUV(P,Input,Integral,LUPRI,IPRINT)
Implicit none
TYPE(Overlap)       :: P
TYPE(IntegralInput) :: Input
TYPE(IntegralItem)  :: Integral
Integer             :: LUPRI,IPRINT
!
Logical             :: Sph1,Sph2
Integer             :: iAngmom,l1,l2,l,ijk1,ijk2,ijk,ioff,nTUV,lm1,lm2,lm
Integer             :: iPrim,iA1,iA2,nCont,nC1,nC2,dim1,iAng1,iAng2
Integer             :: idmat,start1,start2,i1,i2,iC1,iC2,indlm,iTUV
Real(realk),pointer :: Ecoeffs(:,:,:),Spherical(:,:),SpherEcoeff(:,:,:)
Real(realk),pointer :: CC(:,:,:),ContractedDmat(:,:,:)
Real(realk)         :: DM1=-1.0d0,D1=1.0d0,D0=0.0d0,dtemp
Integer             :: eS1,eS2

NULLIFY(P%FTUV)
ALLOCATE(P%FTUV(P%nPrimitives,P%nTUV,Input%NDMAT_RHS))
CALL DZERO(P%FTUV,P%nTUV*P%nPrimitives*Input%NDMAT_RHS)

DO iAngmom = 1,P%nAngmom
  IF (P%hermiteSingle) THEN
    CALL QUIT('ToDo Hermite-single SET_FTUV')
  ELSE
    iA1 = P%indexAng1(iangmom)
    iA2 = P%indexAng2(iangmom)
    l1 = P%orbital1%angmom(iA1)
    l2 = P%orbital2%angmom(iA2)
    l  = l1+l2
    ijk1 = (l1+1)*(l1+2)/2
    ijk2 = (l2+1)*(l2+2)/2
    ijk  = ijk1*ijk2
    nTUV = (l+1)*(l+2)*(l+3)/6
    NULLIFY(Ecoeffs)
    ALLOCATE(Ecoeffs(nTUV,ijk,P%nPrimitives))
    CALL BuildEcoeffTensor1(Integral%TUV,P,DM1,Ecoeffs,nTUV,ijk,P%nPrimitives,iAngmom,1,LUPRI,IPRINT)

!   Spherical transformation of E-coefficients
    Sph1 = P%orbital1%spherical.AND.(l1.GT.1)
    Sph2 = P%orbital2%spherical.AND.(l2.GT.1)
    lm1 = ijk1
    IF (Sph1) lm1 = 2*l1+1
    lm2 = ijk2
    IF (Sph2) lm2 = 2*l2+1
    lm  = lm1*lm2
    IF (Sph1.OR.Sph2) THEN
      NULLIFY(Spherical)
      ALLOCATE(Spherical(ijk,lm))
      CALL ContructSphericalTransformation(Spherical,P%orbital1%SPH_MAT(l1+1)%p%elms,&
     &                                     P%orbital2%SPH_MAT(l2+1)%p%elms,&
     &                                     lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)
      NULLIFY(SpherEcoeff)
      ALLOCATE(SpherEcoeff(nTUV,lm,P%nPrimitives))
      DO iPrim=1,P%nPrimitives


        CALL DGEMM('N','N',nTUV,lm,ijk,1.0d0,Ecoeffs(1,1,iPrim),nTUV,Spherical,ijk,&
     &             0.0d0,SpherEcoeff(1,1,iPrim),nTUV)
      ENDDO
      DEALLOCATE(Ecoeffs)
      DEALLOCATE(Spherical)
    ELSE
      NULLIFY(SpherEcoeff)
      SpherEcoeff => Ecoeffs

    ENDIF
!   Set up contraction matrix
    iA1 = P%indexAng1(iangmom)
    iA2 = P%indexAng2(iangmom)
    nC1 = P%orbital1%nContracted(iA1)
    nC2 = P%orbital2%nContracted(iA2)
    nCont = nC1*nC2
    NULLIFY(CC)
    ALLOCATE(CC(P%nPrimitives,nC1,nC2))
    CALL ConstructContraction(CC,P,1,P%nPrimitives,nC1,nC2,iA1,iA2,LUPRI,IPRINT)

!   Make contraction of density-matrix and contraction coefficients
    NULLIFY(ContractedDmat)
    ALLOCATE(ContractedDmat(P%nPrimitives,lm,Input%NDMAT_RHS))
    CALL DZERO(ContractedDmat,lm*P%nPrimitives*Input%NDMAT_RHS)
    start1 = P%orbital1%startOrbital(iA1,1)
    start2 = P%orbital2%startOrbital(iA2,1)
    DO idmat=1,Input%nDMAT_RHS
      i1 = start1
      DO iC1=1,nC1
        DO iAng1=1,lm1
          i2 = start2
          DO iC2=1,nC2
            DO iAng2=1,lm2
              indlm = (iAng2-1)*lm1 + iAng1
              dtemp = Input%DMAT_RHS(i1,i2,idmat)
              IF ((start1.NE.start2).AND.Input%sameRHSaos) &
     &            dtemp = dtemp + Input%DMAT_RHS(i2,i1,idmat)
              dtemp = dtemp*Input%CoulombFactor
              DO iPrim=1,P%nPrimitives
                ContractedDmat(iPrim,indlm,idmat) = ContractedDmat(iPrim,indlm,idmat) &
     &                                              + dtemp * CC(iPrim,iC1,iC2)
              ENDDO
              i2=i2+1
            ENDDO
          ENDDO
          i1=i1+1
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(CC)

!   Set up FTUV's
    DO idmat=1,Input%NDMAT_RHS
      DO iTUV=1,nTUV
        DO indlm=1,lm
          DO iPrim=1,P%nPrimitives
            P%FTUV(iPrim,iTUV,idmat) = P%FTUV(iPrim,iTUV,idmat) &
     &               + ContractedDmat(iPrim,indlm,idmat)*SpherEcoeff(iTUV,indlm,iPrim)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(ContractedDmat)
    IF (Sph1.OR.Sph2) THEN
      DEALLOCATE(SpherEcoeff)
    ELSE
      NULLIFY(SpherEcoeff)
      DEALLOCATE(Ecoeffs)
    ENDIF
  ENDIF
ENDDO
END SUBROUTINE SET_FTUV

!Finds the number of significant product primitives
SUBROUTINE GET_NPRIMITIVES(nPrimitives,P,Input,Side)
Implicit none
TYPE(IntegralInput) :: Input
TYPE(Overlap)       :: P
Integer             :: nPrimitives
Character*(*)       :: Side
!
Integer              :: i1,i2,a1,a2,a,b,start2
Real(realk)          :: maxgab,maxelm
Real(realk), pointer :: GAB(:,:)

IF (Side.EQ.'LHS') THEN
  GAB => INPUT%pGAB_LHS
  MAXELM=INPUT%PS_MAXELM_RHS
ELSEIF (Side.EQ.'RHS') THEN
  GAB => INPUT%pGAB_RHS
  MAXELM=INPUT%PS_MAXELM_LHS
ELSE
  WRITE(*,'(1X,2A)') 'Error in GET_NPRIMITIVES. Side =',Side
  CALL QUIT('Programming error. Wrong Side in GET_NPRIMITIVES')
ENDIF
nPrimitives = 0
IF (P%type(1:4) .NE. 'Empt')THEN
   DO i1=1,P%orbital1%nPrimitives
      start2 = 1
      IF (P%sameAO) start2 = i1
      DO i2=start2,P%orbital2%nPrimitives
         !   Find maxgab for given primitive
         maxgab = 0.0d0
         IF(INPUT%PS_SCREEN)THEN
            DO a1=1,P%orbital1%nAngmom 
               DO a2=1,P%orbital2%nAngmom
                  a = P%orbital1%startprimOrbital(a1,1)+(i1-1)*P%orbital1%nOrbComp(a1)
                  b = P%orbital2%startprimOrbital(a2,1)+(i2-1)*P%orbital2%nOrbComp(a2)
                  maxgab=MAX(maxgab,GAB(a,b))
               ENDDO
            ENDDO
         ENDIF
         !   Add if maxgab greater than threshold
         IF(INPUT%PS_SCREEN)THEN
            IF( MAXGAB .GT. INPUT%PS_THRESHOLD/MAXELM)THEN
               nPrimitives = nPrimitives + 1
            ENDIF
         ELSE
            nPrimitives = nPrimitives + 1
         ENDIF
      ENDDO
   ENDDO
ELSE
   nPrimitives = 1
ENDIF

END SUBROUTINE GET_NPRIMITIVES

SUBROUTINE SET_ORBITAL(Orb,AOB,integral,LUPRI)
IMPLICIT NONE
TYPE(Orbital) :: Orb
TYPE(AOBATCH) :: AOB
TYPE(IntegralItem)  :: Integral
integer             :: lupri
!
Integer :: ia,na,np,maxc,dim1,dim2,i,L,nc
!
 np   = AOB%nPrimitives
 na   = AOB%nAngmom

 Orb%type          = AOB%type
 Orb%spherical     = AOB%spherical
 Orb%center        = AOB%center
 Orb%nPrimitives   = np 
 Orb%nPasses       = 1 
 Orb%maxContracted = AOB%maxContracted
 Orb%maxAngmom     = AOB%maxAngmom
 Orb%nAngmom       = na
 Orb%CCidentifier  = AOB%CCindex(1)

 NULLIFY(Orb%exponents)
 ALLOCATE(Orb%exponents(np))
 Orb%exponents     = AOB%pExponents%elms(1:np)

 NULLIFY(Orb%angmom)
 ALLOCATE(Orb%angmom(na))
 Orb%angmom(:)     = AOB%angmom(1:na)

 NULLIFY(Orb%nContracted)
 ALLOCATE(Orb%nContracted(na))
 Orb%nContracted   = AOB%nContracted(1:na)

 NULLIFY(Orb%startOrbital)
 ALLOCATE(Orb%startOrbital(na,1))
 Orb%startOrbital(:,1)  = AOB%startOrbital(1:na)

 NULLIFY(Orb%startprimOrbital)
 ALLOCATE(Orb%startprimOrbital(na,1))
 Orb%startprimOrbital(:,1)  = AOB%startprimOrbital(1:na)

 NULLIFY(Orb%nOrbComp)
 ALLOCATE(Orb%nOrbComp(na))
 Orb%nOrbComp  = AOB%nOrbComp(1:na)

 NULLIFY(Orb%nPrimOrbComp)
 ALLOCATE(Orb%nPrimOrbComp(na))
 Orb%nPrimOrbComp  = AOB%nPrimOrbComp(1:na)

 NULLIFY(Orb%nOrbitals)
 ALLOCATE(Orb%nOrbitals(na))
 Orb%nOrbitals  = AOB%nOrbitals(1:na)

! NULLIFY(Orb%CC)
! ALLOCATE(Orb%CC(np,AOB%maxContracted,na))
! DO ia=1,na
!   nc = Orb%nContracted(ia)
!   CALL DCOPY(np*nc,AOB%pCC(ia)%p%elms,1,Orb%CC(1;np,1:nc,ia),1)
! ENDDO
 NULLIFY(Orb%CC)
 ALLOCATE(Orb%CC(na))
 DO ia=1,na
    Orb%CC(ia)%p => AOB%pCC(ia)%p
 ENDDO

 NULLIFY(Orb%SPH_MAT)
 ALLOCATE(Orb%SPH_MAT(Orb%maxAngmom+1))
 DO I=1,Orb%maxAngmom+1
    Orb%SPH_MAT(I)%p => integral%TUV%SPH_MAT(I)
 ENDDO

END SUBROUTINE SET_ORBITAL

SUBROUTINE PRINT_OVERLAP(P,IUNIT,IPRINT,SIDE)
IMPLICIT NONE
TYPE(Overlap) :: P
Integer       :: IUNIT,IPRINT
Character*(*)      :: Side
!
integer :: i,ipass,ip
CALL LSHEADER(IUNIT,'OVERLAP')
IF (SIDE.EQ.'RHS') THEN
   WRITE(IUNIT,'(1X,A)') ''
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') '***   OVERLAP RHS     ***'
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') ''
ELSEIF (SIDE.EQ.'LHS') THEN
   WRITE(IUNIT,'(1X,A)') ''
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') '***   OVERLAP LHS     ***'
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') ''
ENDIF
IF(P%type .NE. 'FTUV')THEN
   WRITE(IUNIT,'(3X,A)') '-------------- orbital 1 --------------'
   CALL PRINT_ORBITAL(P%orbital1,IUNIT,IPRINT)
   WRITE(IUNIT,'(3X,A)') '-------------- orbital 2 --------------'
   CALL PRINT_ORBITAL(P%orbital2,IUNIT,IPRINT)
ENDIF
 WRITE(IUNIT,'(3X,A)') '--------- overlap specifics -----------'
 WRITE(IUNIT,'(5X,A,A)')  'type             =', P%type
 WRITE(IUNIT,'(5X,A,I3)') 'nAngmom          =', P%nAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'nPrimitives      =', P%nPrimitives
 WRITE(IUNIT,'(5X,A,I3)') 'nPasses          =', P%nPasses
 WRITE(IUNIT,'(5X,A,I3)') 'totOrbitals      =', P%totOrbitals
 WRITE(IUNIT,'(5X,A,I3)') 'nTUV             =', P%nTUV
 WRITE(IUNIT,'(5X,A,I3)') 'maxContracted    =', P%maxContracted
 WRITE(IUNIT,'(5X,A,I3)') 'minAngmom        =', P%minAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'maxAngmom        =', P%maxAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'startAngmom      =', P%startAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'endAngmom        =', P%endAngmom
IF(P%type .NE. 'FTUV')THEN
   DO ipass=1,P%nPasses
      WRITE(IUNIT,'(5X,A,3F8.4,A,I3)')'distance12(A) =', (P%distance12(i,ipass), i=1,3),' for pass ',ipass
   ENDDO
      WRITE(IUNIT,'(5X,A,F8.4)') 'squaredDistance  =', P%squaredDistance
ENDIF
 WRITE(IUNIT,'(5X,A,8I3)')'angmom           =', (P%angmom(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'indexAng1        =', (P%indexAng1(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'indexAng2        =', (P%indexAng2(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'nOrbitals        =', (P%nOrbitals(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'nOrbComp         =', (P%nOrbComp(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'nContracted      =', (P%nContracted(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,1L2)')  'sameAO           =  ', P%sameAO
 WRITE(IUNIT,'(3X,A)')    '------------------------------------------------------------------------'
IF(P%type .NE. 'FTUV')THEN
   DO ipass=1,P%nPasses
      WRITE(IUNIT,'(3X,A,I5)') '*** Pass number ',ipass
      WRITE(IUNIT,'(3X,A)')    '------------------------------------------------------------------------'
      WRITE(IUNIT,'(3X,A)')    '- iPrim   Px      Py      Pz      p       mu      pre    iprim1 iprim2 -'
      ip = (ipass-1)*P%nPrimitives
      DO i=1+ip,P%nPrimitives+ip
         WRITE(IUNIT,'(5X,I4,6F8.4,2I7)') i, P%center(1,i), P%center(2,i), P%center(3,i), &
              &  P%exponents(i), P%reducedExponents(i), P%preExpFac(i), P%iprim1(i), P%iprim2(i)
      ENDDO
   ENDDO
ELSE
   DO ipass=1,P%nPasses
      WRITE(IUNIT,'(3X,A,I5)') '*** Pass number ',ipass
      WRITE(IUNIT,'(3X,A)')    '------------------------------------------------------------------------'
      WRITE(IUNIT,'(3X,A)')    '- iPrim   Px      Py      Pz      p       mu      pre'
      ip = (ipass-1)*P%nPrimitives
      DO i=1+ip,P%nPrimitives+ip
         WRITE(IUNIT,'(5X,I4,6F8.4)') i, P%center(1,i), P%center(2,i), P%center(3,i), &
              &  P%exponents(i), P%reducedExponents(i), P%preExpFac(i)
      ENDDO
   ENDDO
ENDIF
WRITE(IUNIT,'(3X,A)')    '------------------------------------------------------------------------'
END SUBROUTINE PRINT_OVERLAP

SUBROUTINE PRINT_ORBITAL(Orb,IUNIT,IPRINT)
TYPE(Orbital) :: Orb
Integer       :: IUNIT,IPRINT
!
Integer :: I,J,K
 WRITE(IUNIT,'(5X,A,A55)')'Type          = ', Orb%type
 WRITE(IUNIT,'(5X,A,1L2)')  'Spherical     = ', Orb%spherical
 WRITE(IUNIT,'(5X,A,I3)') 'maxAngmom     =', Orb%maxAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'nAngmom       =', Orb%nAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'nPrimitives   =', Orb%nPrimitives
 WRITE(IUNIT,'(5X,A,I3)') 'maxContracted =', Orb%maxContracted
 WRITE(IUNIT,'(5X,A,3F10.6)') 'center (A)    =', (Orb%center(i),i=1,3)
 WRITE(IUNIT,'(3X,A)')    '-----------------------------------------------------------'
 WRITE(IUNIT,'(3X,A)')    '-  Block   angmom  #cont.  1.orb.   #orb.   #oc   #pri.oc  1.porb  -'
 DO I=1,Orb%nAngmom
   WRITE(IUNIT,'(1X,8I8)')  I, Orb%angmom(I), Orb%nContracted(I), Orb%startOrbital(I,1), &
     & Orb%nOrbitals(I), Orb%nOrbComp(I), Orb%nPrimOrbComp(I),Orb%startprimOrbital(I,1)
 ENDDO
 WRITE(IUNIT,'(3X,A)')    '-----------------------------------------------------------'
 WRITE(IUNIT,'(3X,A)')    '----------------- Exponents -----------------'
 WRITE(IUNIT,'(5X,5F12.6)') (Orb%exponents(I),I=1,Orb%nPrimitives)
 WRITE(IUNIT,'(3X,A)')    '----------------- Contraction Coefficients -----------------'
 DO K=1,Orb%nAngmom
   WRITE(IUNIT,'(5X,A,I3,A,I3)') '* Angular block number',K,' with angular momentum', Orb%angmom(K)
   DO J=1,Orb%nContracted(K)
!     WRITE(IUNIT,'(5X,I3,5F12.6,/(8X,5F12.6))') J,(Orb%CC(I,J,K),I=1,Orb%nPrimitives)
     WRITE(IUNIT,'(5X,I3,5F12.6,/(8X,5F12.6))') J,(Orb%CC(K)%p%elms(I+(J-1)*Orb%nPrimitives),I=1,Orb%nPrimitives)
   ENDDO
   IF(Orb%angmom(K) .GE. 2)THEN
      WRITE(IUNIT,'(7X,A,I3)')'Spherical transformation matrix for angular momentum', Orb%angmom(K)
      CALL OUTPUT(orb%SPH_MAT(Orb%angmom(K)+1)%p%elms,1,2*Orb%angmom(K)+1,&
         &1,(Orb%angmom(K)+1)*(Orb%angmom(K)+2)/2,2*Orb%angmom(K)+1,&
         &(Orb%angmom(K)+1)*(Orb%angmom(K)+2)/2,1,IUNIT)
   ENDIF
      WRITE(IUNIT,*)' '
 ENDDO

END SUBROUTINE PRINT_ORBITAL

SUBROUTINE Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,IUNIT,IPRINT)
TYPE(Integrand)       :: PQ
TYPE(Overlap),target  :: P,Q
TYPE(IntegralInput)   :: INPUT
Integer               :: IUNIT,IPRINT,ILHS,IRHS
!
Integer     :: ipq,ip,iq,idir
Real(realk) :: d2, pref
Real(realk), parameter :: PI=3.14159265358979323846D0, Nill = 0.0d0, OneHalf=1.5d0
Real(realk), parameter :: Two = 2.0d0, TwoHalf=2.5d0 
logical :: coulomb,nucrep
!
integer :: i
 PQ%Operator = INPUT%operator
 PQ%kinetic  = PQ%Operator(1:7).EQ.'Kinetic'
 PQ%ORIGO(1) = INPUT%ORIGO(1)
 PQ%ORIGO(2) = INPUT%ORIGO(2)
 PQ%ORIGO(3) = INPUT%ORIGO(3)
 PQ%P%p            => P
 PQ%Q%p            => Q
!simen Refine!
 PQ%nPrimitives    = P%nPrimitives*P%nPasses   * Q%nPrimitives*Q%nPasses
 PQ%maxContracted  = P%maxContracted * Q%maxContracted
 PQ%minAngmom      = P%minAngmom     + Q%minAngmom
 PQ%maxAngmom      = P%maxAngmom     + Q%maxAngmom
 PQ%startAngmom    = P%startAngmom   + Q%startAngmom 
 PQ%endAngmom      = P%endAngmom     + Q%endAngmom 
 PQ%samePQ         = INPUT%sameODs .AND. (ILHS .EQ. IRHS)
 PQ%nAngmom        = P%nAngmom       * Q%nAngmom 
 IF (PQ%samePQ) THEN
  PQ%nAngmom = P%nAngmom*(P%nAngmom+1)/2
!  PQ%nPrimitives = P%nPrimitives*(P%nPrimitives+1)/2
 ENDIF
!Dimensions according to nPrimitives
!Nullify(PQ%distance)
!Nullify(PQ%squaredDistance)
!Nullify(PQ%exponents)
!Nullify(PQ%reducedExponents)
!Nullify(PQ%integralPrefactor)
!Nullify(PQ%iprimP)
!Nullify(PQ%iprimQ)
!Allocate(PQ%distance(PQ%nPrimitives,3))
!Allocate(PQ%squaredDistance(PQ%nPrimitives))
!Allocate(PQ%exponents(PQ%nPrimitives))
!Allocate(PQ%reducedExponents(PQ%nPrimitives))
!Allocate(PQ%integralPrefactor(PQ%nPrimitives))
!Allocate(PQ%iprimP(PQ%nPrimitives))
!Allocate(PQ%iprimQ(PQ%nPrimitives))
coulomb = (PQ%Operator.EQ.'Coulomb').OR.(PQ%Operator(1:3).EQ.'Erf')
nucrep = PQ%Operator(1:6).EQ.'Nucrep'

CALL SetIntegrand(PQ%iprimP,PQ%iprimQ,PQ%distance,PQ%squaredDistance,&
     &PQ%exponents,PQ%reducedExponents,coulomb,nucrep,PQ%integralPrefactor,&
     &PQ%nPrimitives,P%exponents,Q%exponents,P%center,Q%center,&
     &P%nPrimitives*P%nPasses,Q%nPrimitives*Q%nPasses,IUNIT)

IF(INPUT%ATTFACTOR) CALL DSCAL(PQ%nPrimitives,1.d0+INPUT%ATTOMEGA/PI,PQ%integralPrefactor,1)

 IF (IPRINT.GT.15) THEN
    CALL PRINT_OVERLAP(P,IUNIT,IPRINT,'LHS')
    CALL PRINT_OVERLAP(Q,IUNIT,IPRINT,'RHS')
    CALL PRINT_Integrand(PQ,IUNIT)
 ENDIF
END SUBROUTINE Build_Integrand

SUBROUTINE SetIntegrand(iprimP,iprimQ,rPQ,squaredDistance,exponents,&
     &                  reducedExponents,coulomb,nucrep,integralPrefactor,npq,&
     &                  pexps,qexps,pcent,qcent,np,nq,IUNIT)
implicit none
integer     :: np,nq,npq,IUNIT
integer     :: iprimP(npq),iprimQ(npq)
real(realk) :: rPQ(3,npq),squaredDistance(npq)
real(realk) :: exponents(npq),reducedExponents(npq),integralPrefactor(npq)
real(realk) :: pexps(np),pcent(3,np)
real(realk) :: qexps(nq),qcent(3,nq)
real(realk) :: qx,qy,qz,pqx,pqy,pqz
logical     :: nucrep
!
integer :: ipq,ip,iq,idir,startq
real(realk) :: pqdist(3),pqdist2,d2,pref,p,q,p_q
Real(realk), parameter :: PI=3.14159265358979323846D0, Nill = 0.0d0, OneHalf=1.5d0
Real(realk), parameter :: Two = 2.0d0, TwoHalf=2.5d0 
Real(realk), parameter :: PIFAC = 34.986836655249725D0 !Two*PI**TwoHalf
Real(realk), parameter :: TWOPI = 6.2831853071795862D0 
!
logical :: coulomb

ipq = 0
DO iq=1, nq
   qx = qcent(1,iq)
   qy = qcent(2,iq)
   qz = qcent(3,iq)
   q  = qexps(iq)
   DO ip=1, np
      !    TEST here
      ipq = ipq+1
      iprimP(ipq) = ip
      iprimQ(ipq) = iq
      pqx = pcent(1,ip) - qx
      pqy = pcent(2,ip) - qy
      pqz = pcent(3,ip) - qz

      rPQ(1,ipq) = pqx
      rPQ(2,ipq) = pqy
      rPQ(3,ipq) = pqz
       
      squaredDistance(ipq) = pqx*pqx+pqy*pqy+pqz*pqz
      p  = pexps(ip)
      p_q = p + q
      exponents(ipq)          = p + q
      reducedExponents(ipq)   = p*q/p_q

      IF (coulomb) THEN
        pref = PIFAC/(p*q*SQRT(p_q))
      ELSEIF(nucrep)THEN
        pref = TWOPI/p
      ELSE
        pref = (PI/p_q)**(OneHalf)
      ENDIF
      integralPrefactor(ipq)  = pref
   ENDDO
ENDDO
END SUBROUTINE SetIntegrand

SUBROUTINE PRINT_Integrand(PQ,IUNIT)
TYPE(Integrand)       :: PQ
Integer               :: IUNIT
!
Integer :: i
!
CALL LSHEADER(IUNIT,'Integrand output')
WRITE(IUNIT,'(1X,A)') ''
WRITE(IUNIT,'(1X,A)') '*************************'
WRITE(IUNIT,'(1X,A)') '***     Integrand     ***'
WRITE(IUNIT,'(1X,A)') '*************************'
WRITE(IUNIT,'(1X,A)') ''
WRITE(IUNIT,'(5X,2A)')   'Operator      = ', PQ%Operator
WRITE(IUNIT,'(5X,A,I4)') 'nAngmom       =', PQ%nAngmom
WRITE(IUNIT,'(5X,A,I4)') 'nPrimitives   =', PQ%nPrimitives
WRITE(IUNIT,'(5X,A,I4)') 'maxContracted =', PQ%maxContracted
WRITE(IUNIT,'(5X,A,I4)') 'minAngmom     =', PQ%minAngmom
WRITE(IUNIT,'(5X,A,I4)') 'maxAngmom     =', PQ%maxAngmom
WRITE(IUNIT,'(5X,A,I4)') 'startAngmom   =', PQ%startAngmom
WRITE(IUNIT,'(5X,A,I4)') 'endAngmom     =', PQ%endAngmom
WRITE(IUNIT,'(3X,A)')    '--------------------------------------------------------------------------'
WRITE(IUNIT,'(3X,A)')    '- prim  iP  iQ   X_PQ    Y_PQ    Z_PQ   R_PQ^2    exp   redExp  intPre   -'
DO i=1,PQ%nPrimitives
WRITE(IUNIT,'(5X,3I4,7F8.4)') i,PQ%iprimP(i),PQ%iprimQ(i),&
    & PQ%distance(1,i),PQ%distance(2,i),PQ%distance(3,i),&
    & PQ%squaredDistance(i),PQ%exponents(i),PQ%reducedExponents(i),PQ%integralPrefactor(i)
ENDDO
WRITE(IUNIT,'(3X,A)')    '--------------------------------------------------------------------------'

END SUBROUTINE PRINT_Integrand

SUBROUTINE Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
implicit none     
TYPE(Integrand),target  :: PQ
TYPE(Integralitem)      :: integral
TYPE(Integralinput)     :: input
INTEGER                 :: LUPRI,IPRINT

SELECT CASE(PQ%Operator)
CASE ('Overlap')
  CALL GET_WTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT,'Overlap')
CASE ('Coulomb')
  CALL GET_WTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT,'Coulomb')
CASE ('Erfc   ')
  CALL GET_WTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT,'Erfc   ')
CASE ('Erf    ')
  CALL GET_WTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT,'Erf    ')
CASE ('Nucrep')
  CALL GET_WTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT,'Nucrep ')
CASE ('Kinetic')
  CALL GET_WTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT,'Overlap')
CASE ('Carmom')
  !Cartesian multipole order
  CALL GET_MOMTUV(INTEGRAL,PQ,LUPRI,IPRINT,INPUT%DERIVORDER)
CASE DEFAULT
  WRITE(LUPRI,'(1X,A,A16)') 'Programming error! Not a case in Build_HermiteTUV:',PQ%Operator 
END SELECT
END SUBROUTINE Build_HermiteTUV

SUBROUTINE GET_WTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT,OPERATORLABEL)
! The final S(T,U,V) integrals are arranged as follows:
! S(000)
! S(100) S(010) S(001)
! S(200) S(110) S(101) S(020) S(011) S(002)
! S(300) S(210) S(201) S(120) S(111) S(102) S(030) S(021) S(012) S(003)
implicit none
Character(len=7)        :: Operatorlabel
TYPE(Integrand)         :: PQ
TYPE(Integralitem)      :: integral
TYPE(Integralinput)     :: input
real(realk),ALLOCATABLE :: INTTEMP(:,:)!nPrim,TUV
real(realk),ALLOCATABLE :: INTTEMP2(:,:)!nPrim,TUV
!real(realk),ALLOCATABLE :: SJ000(:,:)
real(realk),ALLOCATABLE :: SJ0002(:,:)
INTEGER                 :: LUPRI,SUM,J,K,T,U,V,TUV,IOFF
INTEGER                 :: nPrim,IPRINT,ntuv,L,I
INTEGER,ALLOCATABLE     :: TUVindex(:,:,:)
INTEGER                 :: zeroX,zeroY,zeroZ,Jmax,Jstart
real(realk)             :: X0,Y0,Z0
!
real(realk),allocatable :: wtuvnew(:,:),Rpq(:,:),inttemp000(:,:)
real(realk) :: diff
real(4) :: tarr(2),etime,tnew,told,tstart

NPrim=PQ%nPrimitives
JMAX=PQ%endAngmom
Jstart=PQ%startAngmom

SELECT CASE(OPERATORLABEL)
CASE ('Overlap')
   CALL BUILD_SJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT)
CASE ('Coulomb')
   CALL BUILD_RJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral)
CASE ('Erfc   ')
   ALLOCATE(SJ0002(0:JMAX,nPrim))
   CALL BUILD_RJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral)
   CALL BUILD_ERF_RJ000(SJ0002,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,INPUT%ATTomega)
   CALL DAXPY(NPrim*(JMAX+1),INPUT%ATTbeta,SJ0002,1,Integral%IN,1)
   DEALLOCATE(SJ0002)
CASE ('Erf    ')
   CALL BUILD_ERF_RJ000(Integral%IN,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,INPUT%ATTomega)
CASE ('Nucrep')
   CALL BUILD_NUCLEAR_RJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral)
CASE DEFAULT
  WRITE(LUPRI,'(1X,A,A16)') 'Programming error! Not a case in GET_WTUV:',OPERATORLABEL
END SELECT

IF (IPRINT .GE. 10) THEN
   CALL LSHEADER(LUPRI,'Output from W000')
   DO I=1,NPrim
      DO J=0,Jmax
         WRITE(LUPRI,'(2X,A6,I4,A1,I2,A2,F16.9)')'SJ000(',J,',',I,')=',Integral%IN(1+J+(I-1)*(Jmax+1))
      ENDDO
   ENDDO
END IF

ntuv=(JMAX+1)*(JMAX+2)*(JMAX+3)/6-Jstart*(Jstart+1)*(Jstart+2)/6
INTEGRAL%nTUV=ntuv
Integral%nEFG=1 !no cartesian multipole moments

CALL WTUVrecurrence(Integral%Rtuv,Integral%IN,Integral%TUV,PQ%distance,&
     &                jStart,jMax,nPrim,ntuv,lupri,iprint)

IOFF = 0

IF (IPRINT .GE. 10) THEN
   CALL LSHEADER(LUPRI,'Output from WTUVrecurrence')
   WRITE (LUPRI,'(2X,A,I10)') 'JSTART', JSTART
   WRITE (LUPRI,'(2X,A,I10)') 'JMAX  ', JMAX
   WRITE (LUPRI,'(2X,A,I10)') 'NPrim  ', NPrim
   WRITE (LUPRI,'(2X,A,I10)') 'NTUV  ', nTUV
   IF (IPRINT .GE. 20) THEN
      CALL LSHEADER(LUPRI,'Hermite integrals S(t,u,v)')
      DO J = Jstart, JMAX
         DO T = J,0,-1
            DO U = J-T,0,-1
               V=J-T-U
               TUV=Integral%TUV%tuvIndex(T,U,V)
               WRITE (LUPRI,'(2X,A2,I1,A1,I1,A1,I1,A1,2X,5F12.8/,(12X,5F12.8))')&
                    & 'W(',T,',',U,',',V,')', &
                    &(INTEGRAL%Rtuv(TUV+(i-1)*nTUV),I=1,NPrim)
               WRITE (LUPRI,*) ' '
            ENDDO
         ENDDO
      ENDDO
   END IF
END IF

Integral%nPrim=nPrim

END SUBROUTINE GET_WTUV

SUBROUTINE WTUVrecurrence(WTUV,WJ000,sharedTUV,Rpq,jmin,jmax,&
     &                      nPrim,ntuv,lupri,iprint)
implicit none
Integer           :: jmin,jmax,nPrim,ntuv,lupri,iprint
INTEGER,PARAMETER :: MAXJ = 50
TYPE(TUVitem)     :: SharedTUV
Real(realk)       :: WTUV(ntuv,nPrim),WJ000(0:jmax,nPrim),Rpq(3,nPrim)
!
Real(realk) :: Xtemp,Ytemp,Ztemp
Integer     :: i,n,j,t,u,v,ntuvfull,ituv,ituvfull,ioff,j2,nx,ijk,jstart,jend,dir

integer::k
Real(realk),pointer :: CUR(:,:),OLD(:,:),TEMP(:,:)
Real(realk) :: W100,W010,W001,W200,W110,W101,W020,W011,W002,WJ,WJ1,WJ2,WJ3
integer :: tm1,um1,vm1,tm2,um2,vm2,m1,ituvm,ituvm2,ioffp


IF (jmax.EQ.0) THEN
  CALL DCOPY(nPrim,WJ000,1,WTUV,1)
ELSE IF (jmax.EQ.1) THEN
  IF (jmin.EQ.1) THEN
    DO n=1,nPrim
      WJ1 = WJ000(1,n)
      WTUV(1,n) = Rpq(1,n)*WJ1
      WTUV(2,n) = Rpq(2,n)*WJ1
      WTUV(3,n) = Rpq(3,n)*WJ1
    ENDDO
  ELSE
    DO n=1,nPrim
      WTUV(1,n) = WJ000(0,n)
      WJ1 = WJ000(1,n)
      WTUV(2,n) = Rpq(1,n)*WJ1
      WTUV(3,n) = Rpq(2,n)*WJ1
      WTUV(4,n) = Rpq(3,n)*WJ1
    ENDDO
  ENDIF
ELSE IF (jmax.EQ.2) THEN
   IF (jmin.EQ.2) THEN
    DO n=1,nPrim
      Xtemp = Rpq(1,n)
      Ytemp = Rpq(2,n)
      Ztemp = Rpq(3,n)
      WJ1 = WJ000(1,n)
      WJ2 = WJ000(2,n)
      WTUV(1,n) = WJ1 + Xtemp*Xtemp*WJ2
      WTUV(2,n) =       Xtemp*Ytemp*WJ2
      WTUV(3,n) =       Xtemp*Ztemp*WJ2
      WTUV(4,n) = WJ1 + Ytemp*Ytemp*WJ2
      WTUV(5,n) =       Ytemp*Ztemp*WJ2
      WTUV(6,n) = WJ1 + Ztemp*Ztemp*WJ2
    ENDDO
   ELSE IF (jmin.EQ.1) THEN
    DO n=1,nPrim
      Xtemp = Rpq(1,n)
      Ytemp = Rpq(2,n)
      Ztemp = Rpq(3,n)
      WJ1 = WJ000(1,n)
      WJ2 = WJ000(2,n)
      WTUV(1,n) = Xtemp*WJ1
      WTUV(2,n) = Ytemp*WJ1
      WTUV(3,n) = Ztemp*WJ1
      WTUV(4,n) = WJ1 + Xtemp*Xtemp*WJ2
      WTUV(5,n) =       Xtemp*Ytemp*WJ2
      WTUV(6,n) =       Xtemp*Ztemp*WJ2
      WTUV(7,n) = WJ1 + Ytemp*Ytemp*WJ2
      WTUV(8,n) =       Ytemp*Ztemp*WJ2
      WTUV(9,n) = WJ1 + Ztemp*Ztemp*WJ2
    ENDDO
   ELSE
    DO n=1,nPrim
      Xtemp = Rpq(1,n)
      Ytemp = Rpq(2,n)
      Ztemp = Rpq(3,n)
      WJ  = WJ000(0,n)
      WJ1 = WJ000(1,n)
      WJ2 = WJ000(2,n)
      WTUV(1,n)  = WJ
      WTUV(2,n)  = Xtemp*WJ1
      WTUV(3,n)  = Ytemp*WJ1
      WTUV(4,n)  = Ztemp*WJ1
      WTUV(5,n)  = WJ1 + Xtemp*Xtemp*WJ2
      WTUV(6,n)  =       Xtemp*Ytemp*WJ2
      WTUV(7,n)  =       Xtemp*Ztemp*WJ2
      WTUV(8,n)  = WJ1 + Ytemp*Ytemp*WJ2
      WTUV(9,n)  =       Ytemp*Ztemp*WJ2
      WTUV(10,n) = WJ1 + Ztemp*Ztemp*WJ2
    ENDDO
   ENDIF
ELSE  ! J > 2
  ntuvfull = (jmax+1)*(jmax+2)*(jmax+3)/6
  ioffp    = jmin*(jmin+1)*(jmin+2)/6+1
  NULLIFY(CUR)
  NULLIFY(OLD)
  ALLOCATE(CUR(ntuvfull,nPrim))
  ALLOCATE(OLD(ntuvfull,nPrim))

  DO j=jmax-3,0,-1
    TEMP => CUR
    CUR  => OLD
    OLD  => TEMP
    DO n=1,nPrim
      Xtemp = Rpq(1,n)
      Ytemp = Rpq(2,n)
      Ztemp = Rpq(3,n)
      WJ   = WJ000(j,n)
      WJ1  = WJ000(j+1,n)
      WJ2  = WJ000(j+2,n)
      WJ3  = WJ000(j+3,n)
      W100 = Xtemp*WJ2
      W010 = Ytemp*WJ2
      W001 = Ztemp*WJ2
      W200 = WJ2 + Xtemp*Xtemp*WJ3
      W110 =       Xtemp*Ytemp*WJ3
      W101 =       Xtemp*Ztemp*WJ3
      W020 = WJ2 + Ytemp*Ytemp*WJ3
      W011 =       Ytemp*Ztemp*WJ3
      W002 = WJ2 + Ztemp*Ztemp*WJ3
      CUR(1,n)  = WJ                           !000
      CUR(2,n)  = Xtemp*WJ1                    !100
      CUR(3,n)  = Ytemp*WJ1                    !010
      CUR(4,n)  = Ztemp*WJ1                    !001
      CUR(5,n)  = WJ1 + Xtemp*W100             !200
      CUR(6,n)  =       Ytemp*W100             !110
      CUR(7,n)  =       Ztemp*W100             !101
      CUR(8,n)  = WJ1 + Ytemp*W010             !020
      CUR(9,n)  =       Ztemp*W010             !011
      CUR(10,n) = WJ1 + Ztemp*W001             !002
      CUR(11,n) = 2.0d0*W100 + Xtemp*W200      !300
      CUR(12,n) =       W010 + Xtemp*W110      !210
      CUR(13,n) =       W001 + Xtemp*W101      !201
      CUR(14,n) =       W100 + Ytemp*W110      !120
      CUR(15,n) =              Xtemp*W011      !111
      CUR(16,n) =       W100 + Ztemp*W101      !102
      CUR(17,n) = 2.0d0*W010 + Ytemp*W020      !030
      CUR(18,n) =       W001 + Ytemp*W011      !021
      CUR(19,n) =       W010 + Ztemp*W011      !012
      CUR(20,n) = 2.0d0*W001 + Ztemp*W002      !003
    ENDDO
    ituv = 20
    DO k=4,jmax-j
      DO t=k,0,-1
        DO u=k-t,0,-1
          v=k-t-u
          IF ((t.ge.u).AND.(t.ge.v)) THEN
            tm1 = t-1
            um1 = u
            vm1 = v
            tm2 = t-2
            um2 = u
            vm2 = v
            m1  = tm1
            dir = 1
          ELSEIF (u.ge.v) THEN
            tm1 = t
            um1 = u-1
            vm1 = v
            tm2 = t
            um2 = u-2
            vm2 = v
            m1  = um1
            dir = 2
          ELSE
            tm1 = t
            um1 = u
            vm1 = v-1
            tm2 = t
            um2 = u
            vm2 = v-2
            m1  = vm1
            dir = 3
          ENDIF
          ituv   = ituv + 1
          ituvm  = sharedTUV%tuvIndex(tm1,um1,vm1)
          ituvm2 = sharedTUV%tuvIndex(tm2,um2,vm2)
          DO n=1,nPrim
            CUR(ituv,n) = m1*OLD(ituvm2,n) + Rpq(dir,n)*OLD(ituvm,n)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO n=1,nPrim
    CALL DCOPY(ntuv,CUR(ioffp,n),1,WTUV(1,n),1)
  ENDDO
  DEALLOCATE(CUR)
  DEALLOCATE(OLD)
ENDIF

END SUBROUTINE WTUVrecurrence 

SUBROUTINE RECURRENCE(CUR,OLD,JVAL,RJ000,Xpq,Ypq,Zpq,JMAX,&
     &                  NPrim,NTUV,TUVindex,zeroX,zeroY,zeroZ)
IMPLICIT NONE
real(realk)             :: CUR(nPrim,NTUV)
real(realk)             :: OLD(nPrim,NTUV)
real(realk)             :: RJ000(nPrim,0:Jmax)
INTEGER                 :: J,K,T,U,V,TUV,JVAL,JMAX,NTUV
INTEGER                 :: nPrim
INTEGER,PARAMETER       :: MAXJ = 50
INTEGER                 :: TUVindex(0:MAXJ,0:MAXJ,0:MAXJ)
real(realk)             :: Xpq(nPrim),Ypq(nPrim),Zpq(nPrim)
INTEGER                 :: MAXT,MAXU,MAXV,M1T,M2T,M1U,M2U
INTEGER                 :: M1V,M2V,IUMAX,I
INTEGER                 :: zeroX,zeroY,zeroZ,M,MP1
LOGICAL                 :: NOTSAMEX,NOTSAMEY,NOTSAMEZ
Real(realk)             :: TMIN1,UMIN1,VMIN1

NOTSAMEX=zeroX .EQ. 1
NOTSAMEY=zeroY .EQ. 1
NOTSAMEZ=zeroZ .EQ. 1

IF (JVAL .EQ. 1) THEN
   !     JVAL = 1
   !     ========
   CALL DCOPY(NPrim,RJ000(1,JMAX-1),1,CUR(1,1),1)
   !LOOP CAN BE PREVENTED IF SAME CENTER
   IF(NOTSAMEX) THEN
!      DO I = 1, Nprim
!         CUR(I,2) = Xpq(I)*RJ000(I,JMAX)
!      ENDDO

      M = MOD(NPrim,4)
      IF(M .EQ. 0) THEN
         DO I = 1,NPrim,4
            CUR(I,2)   = Xpq(I)*RJ000(I,JMAX)
            CUR(I+1,2) = Xpq(I+1)*RJ000(I+1,JMAX)
            CUR(I+2,2) = Xpq(I+2)*RJ000(I+2,JMAX)
            CUR(I+3,2) = Xpq(I+3)*RJ000(I+3,JMAX) 
         ENDDO
      ELSE
         DO I = 1,M
            CUR(I,2) = Xpq(I)*RJ000(I,JMAX) 
         ENDDO
         IF(NPrim .GT. 4)THEN
            MP1=M+1
            DO I = MP1,NPrim,4
               CUR(I,2)   = Xpq(I)*RJ000(I,JMAX) 
               CUR(I+1,2) = Xpq(I+1)*RJ000(I+1,JMAX) 
               CUR(I+2,2) = Xpq(I+2)*RJ000(I+2,JMAX) 
               CUR(I+3,2) = Xpq(I+3)*RJ000(I+3,JMAX) 
            ENDDO
         ENDIF
      ENDIF
   ENDIF
   IF(NOTSAMEY) THEN
!      DO I = 1, Nprim
!         CUR(I,3) = Ypq(I)*RJ000(I,JMAX)
!      ENDDO
      M = MOD(NPrim,4)
      IF(M .EQ. 0) THEN
         DO I = 1,NPrim,4
            CUR(I,3)   = Ypq(I)*RJ000(I,JMAX)
            CUR(I+1,3) = Ypq(I+1)*RJ000(I+1,JMAX)
            CUR(I+2,3) = Ypq(I+2)*RJ000(I+2,JMAX)
            CUR(I+3,3) = Ypq(I+3)*RJ000(I+3,JMAX) 
         ENDDO
      ELSE
         DO I = 1,M
            CUR(I,3) = Ypq(I)*RJ000(I,JMAX) 
         ENDDO
         IF(NPrim .GT. 4)THEN
            MP1=M+1
            DO I = MP1,NPrim,4
               CUR(I,3)   = Ypq(I)*RJ000(I,JMAX) 
               CUR(I+1,3) = Ypq(I+1)*RJ000(I+1,JMAX) 
               CUR(I+2,3) = Ypq(I+2)*RJ000(I+2,JMAX) 
               CUR(I+3,3) = Ypq(I+3)*RJ000(I+3,JMAX) 
            ENDDO
         ENDIF
      ENDIF
   ENDIF
   IF(NOTSAMEZ) THEN
!      DO I = 1, Nprim
!         CUR(I,4) = Zpq(I)*RJ000(I,JMAX)
!      ENDDO
      M = MOD(NPrim,4)
      IF(M .EQ. 0) THEN
         DO I = 1,NPrim,4
            CUR(I,4)   = Zpq(I)*RJ000(I,JMAX)
            CUR(I+1,4) = Zpq(I+1)*RJ000(I+1,JMAX)
            CUR(I+2,4) = Zpq(I+2)*RJ000(I+2,JMAX)
            CUR(I+3,4) = Zpq(I+3)*RJ000(I+3,JMAX) 
         ENDDO
      ELSE
         DO I = 1,M
            CUR(I,4) = Zpq(I)*RJ000(I,JMAX) 
         ENDDO
         IF(NPrim .GT. 4)THEN
            MP1=M+1
            DO I = MP1,NPrim,4
               CUR(I,4)   = Zpq(I)*RJ000(I,JMAX) 
               CUR(I+1,4) = Zpq(I+1)*RJ000(I+1,JMAX) 
               CUR(I+2,4) = Zpq(I+2)*RJ000(I+2,JMAX) 
               CUR(I+3,4) = Zpq(I+3)*RJ000(I+3,JMAX) 
            ENDDO
         ENDIF
      ENDIF
   ENDIF
ELSE
   
   !     JVAL > 1
   !     ========
   
   MAXT   = JMAX
   MAXU   = JMAX
   MAXV   = JMAX
   
   !        R(0,0,0)
   
   CALL DCOPY(Nprim,RJ000(1,JMAX-JVAL),1,CUR,1)
   
   !        R(T,0,0)
   
   IF (NOTSAMEX) THEN !NOT SAME CENTER
      DO I = 1, NPrim
         CUR(I,2) = Xpq(I)*OLD(I,1)
      ENDDO
      DO T = 2, MIN(MAXT,JVAL)
!         TMIN1 = DFLOAT(T - 1)
         TMIN1 = T - 1
         TUV   = TUVINDEX(T  ,0,0)
         M1T   = TUVINDEX(T-1,0,0)
         M2T   = TUVINDEX(T-2,0,0)
         DO I = 1, Nprim
            CUR(I,TUV) = Xpq(I)*OLD(I,M1T) + TMIN1*OLD(I,M2T)
         ENDDO
      ENDDO
   ELSE
      DO T = 2, MIN(MAXT,JVAL), 2
!        TMIN1 = DFLOAT(T - 1)
         TMIN1 = T - 1
         TUV   = TUVINDEX(T,0,0)
         M2T   = TUVINDEX(T-2,0,0)

!         CALL DCOPY(Nprim,OLD(1,M2T),1,CUR(1,TUV),1)
!         CALL DSCAL(Nprim,TMIN1,CUR(1,TUV),1) 
!         DO I = 1, NPrim
!            CUR(I,TUV) = TMIN1*OLD(I,M2T)
!         ENDDO

         M = MOD(NPrim,4)
         IF(M .EQ. 0) THEN
            DO I = 1,NPrim,4
               CUR(I,TUV)   = TMIN1*OLD(I,M2T)   
               CUR(I+1,TUV) = TMIN1*OLD(I+1,M2T)   
               CUR(I+2,TUV) = TMIN1*OLD(I+2,M2T)   
               CUR(I+3,TUV) = TMIN1*OLD(I+3,M2T)   
            ENDDO
         ELSE
            DO I = 1,M
               CUR(I,TUV) = TMIN1*OLD(I,M2T)
            ENDDO
            IF(NPrim .GT. 4)THEN
               MP1=M+1
               DO I = MP1,NPrim,4
                  CUR(I,TUV)   = TMIN1*OLD(I,M2T)   
                  CUR(I+1,TUV) = TMIN1*OLD(I+1,M2T)   
                  CUR(I+2,TUV) = TMIN1*OLD(I+2,M2T)   
                  CUR(I+3,TUV) = TMIN1*OLD(I+3,M2T)   
               ENDDO
            ENDIF
         ENDIF

      ENDDO
   END IF
   
   !        R(T,U,0)
   
   IF (NOTSAMEY) THEN !NOT SAME CENTER
      DO T = 0, MIN(MAXT,JVAL - 1), zeroX
         TUV = TUVINDEX(T,1,0)
         M1U = TUVINDEX(T,0,0)
         DO I = 1, NPrim
            CUR(I,TUV) = Ypq(I)*OLD(I,M1U)
         ENDDO
      ENDDO
      DO U = 2, MIN(MAXU,JVAL)
!        UMIN1  = DFLOAT(U - 1)
         UMIN1  = U - 1
         DO T = 0, MIN(MAXT,JVAL - U), zeroX
            TUV = TUVINDEX(T,U  ,0)
            M1U = TUVINDEX(T,U-1,0)
            M2U = TUVINDEX(T,U-2,0)
            DO I = 1, Nprim
               CUR(I,TUV) = Ypq(I)*OLD(I,M1U) + UMIN1*OLD(I,M2U)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO U = 2, MIN(MAXU,JVAL), 2
!        UMIN1  = DFLOAT(U - 1)
         UMIN1  = U - 1
         DO T = 0, MIN(MAXT,JVAL - U), zeroX
            TUV = TUVINDEX(T,U  ,0)
            M2U = TUVINDEX(T,U-2,0)
!            CALL DCOPY(Nprim,OLD(1,M2U),1,CUR(1,TUV),1)
!            CALL DSCAL(nPrim,UMIN1,CUR(1,TUV),1) 
!            DO I = 1, NPrim
!               CUR(I,TUV) = UMIN1*OLD(I,M2U)
!            ENDDO

            M = MOD(NPrim,4)
            IF(M .EQ. 0) THEN
               DO I = 1,NPrim,4
                  CUR(I,TUV)   = UMIN1*OLD(I,M2U)   
                  CUR(I+1,TUV) = UMIN1*OLD(I+1,M2U)   
                  CUR(I+2,TUV) = UMIN1*OLD(I+2,M2U)   
                  CUR(I+3,TUV) = UMIN1*OLD(I+3,M2U)   
               ENDDO
            ELSE
               DO I = 1,M
                  CUR(I,TUV) = UMIN1*OLD(I,M2U)
               ENDDO
               IF(NPrim .GT. 4)THEN
                  MP1=M+1
                  DO I = MP1,NPrim,4
                     CUR(I,TUV)   = UMIN1*OLD(I,M2U)   
                     CUR(I+1,TUV) = UMIN1*OLD(I+1,M2U)   
                     CUR(I+2,TUV) = UMIN1*OLD(I+2,M2U)   
                     CUR(I+3,TUV) = UMIN1*OLD(I+3,M2U)   
                  ENDDO
               ENDIF
            ENDIF
            
         ENDDO
      ENDDO
   END IF
   
   !        R(T,U,V)
   
   IF (NOTSAMEZ) THEN !IF NOT SAME CENTER
      IUMAX  = JVAL - 1
      DO U = 0, MIN(MAXU,IUMAX), zeroY
         DO T = 0, MIN(MAXT,IUMAX - U), zeroX
            TUV = TUVINDEX(T,U,1)
            M1V = TUVINDEX(T,U,0)
            DO I = 1, NPrim
               CUR(I,TUV) = Zpq(I)*OLD(I,M1V)
            ENDDO
         ENDDO
      ENDDO
      DO V = 2, MIN(MAXV,JVAL)
!        VMIN1  = DFLOAT(V - 1)
         VMIN1  = V - 1
         IUMAX  = JVAL - V
         DO U = 0, MIN(MAXU,IUMAX), zeroY
            DO T = 0, MIN(MAXT,IUMAX - U), zeroX
               TUV = TUVINDEX(T,U,V  )
               M1V = TUVINDEX(T,U,V-1)
               M2V = TUVINDEX(T,U,V-2)
               DO I = 1, NPrim
                  CUR(I,TUV) = Zpq(I)*OLD(I,M1V)+VMIN1*OLD(I,M2V)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO V = 2, MIN(MAXV,JVAL), 2
!        VMIN1  = DFLOAT(V - 1)
         VMIN1  = V - 1
         IUMAX  = JVAL - V
         DO U = 0, MIN(MAXU,IUMAX), zeroY
            DO T = 0, MIN(MAXT,IUMAX - U), zeroX
               TUV = TUVINDEX(T,U,V  )
               M2V = TUVINDEX(T,U,V-2)
!               CALL DCOPY(Nprim,OLD(1,M2V),1,CUR(1,TUV),1)
!               CALL DSCAL(nPrim,VMIN1,CUR(1,TUV),1) 
!               DO I = 1, Nprim
!                  CUR(I,TUV) = VMIN1*OLD(I,M2V)
!               ENDDO


            M = MOD(NPrim,4)
            IF(M .EQ. 0) THEN
               DO I = 1,NPrim,4
                  CUR(I,TUV)   = VMIN1*OLD(I,M2V)   
                  CUR(I+1,TUV) = VMIN1*OLD(I+1,M2V)   
                  CUR(I+2,TUV) = VMIN1*OLD(I+2,M2V)   
                  CUR(I+3,TUV) = VMIN1*OLD(I+3,M2V)   
               ENDDO
            ELSE
               DO I = 1,M
                  CUR(I,TUV) = VMIN1*OLD(I,M2V)
               ENDDO
               IF(NPrim .GT. 4)THEN
                  MP1=M+1
                  DO I = MP1,NPrim,4
                     CUR(I,TUV)   = VMIN1*OLD(I,M2V)   
                     CUR(I+1,TUV) = VMIN1*OLD(I+1,M2V)   
                     CUR(I+2,TUV) = VMIN1*OLD(I+2,M2V)   
                     CUR(I+3,TUV) = VMIN1*OLD(I+3,M2V)   
                  ENDDO
               ENDIF
            ENDIF

            ENDDO
         ENDDO
      ENDDO
   END IF
END IF

END SUBROUTINE RECURRENCE

SUBROUTINE BUILD_SJ000(SJ000,PQ,nPrim,Jmax,LUPRI,IPRINT)
IMPLICIT NONE
TYPE(Integrand)        :: PQ
!real(realk)            :: SJ000(nPrim,PQ%startangmom:PQ%maxAngmom)
real(realk)            :: SJ000(0:Jmax,nPrim)
INTEGER                :: nPrim,J,I,LUPRI,IPRINT,Jmax

CALL BUILD_SJ000_OP(SJ000,Nprim,Jmax,PQ%reducedExponents&
     &,PQ%exponents,PQ%squaredDistance)

IF(IPRINT .GT. 25) THEN
   DO I=1,NPrim
      DO J=0,Jmax
         WRITE(LUPRI,'(2X,A5,I4,A1,I2,A2,F16.9)')'SJ000(',I,',',J,')=',SJ000(J,I)
      ENDDO
   ENDDO
ENDIF
   
END SUBROUTINE BUILD_SJ000

SUBROUTINE BUILD_SJ000_OP(SJ000,nPrim,Jmax,Alpha,P,R2)
IMPLICIT NONE
!real(realk)            :: SJ000(nPrim,0:Jmax)
real(realk)            :: SJ000(0:Jmax,nPrim)
Real(realk), parameter :: PI=3.14159265358979323846D0
Real(realk), parameter :: OneHalf=1.5d0, Two = 2.0d0
INTEGER                :: nPrim,J,I,Jmax
REAL(REALK)            :: ALPHA(nPrim),P(nPrim),R2(nprim)
!
REAL(REALK)            :: TEMP(nprim),TEMP2(nprim)

DO I=1,NPrim
   TEMP(I)=Alpha(I)*R2(I)
ENDDO
DO I=1,NPrim
   TEMP2(I)=EXP(-TEMP(I))
ENDDO
DO I=1,NPrim
   DO J=0,Jmax
      SJ000(J,I)=((-Two*Alpha(I))**J)&
           &*((PI/P(I))**OneHalf)*TEMP2(I)
   ENDDO
ENDDO
END SUBROUTINE BUILD_SJ000_OP

SUBROUTINE Contract_P(INTEGRAL,PQ,input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: INTEGRAL
TYPE(IntegralInput):: INPUT
TYPE(Integrand)    :: PQ
Integer            :: LUPRI,IPRINT,WORKLENGTH,LWORK 
Real(REALK)        :: WORK(WORKLENGTH)
!
Integer  :: iAngmom,startMomQ,iOrbital,ideriv,nderiv
!
REAL(REALK)          :: CPU1,CPU2,WALL1,WALL2

nderiv=INPUT%NderivP

CALL swapRealPointers(Integral%RTUV,Integral%integrals)

iOrbital = 1
!
!Loop over angular contributions sharing the set of primitive functions
!DO iDerivP=1,nderivP NOT YET IMPLEMENTED
DO iAngmom=1,PQ%P%p%nAngmom
   ! Distribute Hermite 2-el integrals WTUV to <tuv1|w|tuv2>
   CALL GETTIM(CPU1,WALL1)
   CALL DistributeHermiteP(Integral,PQ,iAngmom,LUPRI,IPRINT)
   CALL GETTIM(CPU2,WALL2)
   CPUTIME_DistributeHermiteP = CPUTIME_DistributeHermiteP+CPU2-CPU1
   WALLTIME_DistributeHermiteP = WALLTIME_DistributeHermiteP+ WALL2-WALL1
   CALL ContractEcoeff(Integral,PQ%P%p,1.d0,iAngmom,LUPRI,IPRINT)
   CALL GETTIM(CPU1,WALL1)
   CPUTIME_ContractEcoeffP = CPUTIME_ContractEcoeffP+CPU1-CPU2
   WALLTIME_ContractEcoeffP = WALLTIME_ContractEcoeffP+ WALL1-WALL2
   CALL ContractBasis(Integral,PQ%P%p,iAngmom,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
   CALL GETTIM(CPU2,WALL2)
   CPUTIME_ContractBasisP = CPUTIME_ContractBasisP+CPU2-CPU1
   WALLTIME_ContractBasisP = WALLTIME_ContractBasisP+ WALL2-WALL1
   CALL SphericalTransform(Integral,PQ%P%p,iAngmom,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
   CALL GETTIM(CPU1,WALL1)
   CPUTIME_SphericalTransformP = CPUTIME_SphericalTransformP+CPU1-CPU2
   WALLTIME_SphericalTransformP = WALLTIME_SphericalTransformP+ WALL1-WALL2 
   CALL AddToPQ(Integral,PQ,iAngmom,iOrbital,LUPRI,IPRINT)
   CALL GETTIM(CPU2,WALL2)
   CPUTIME_AddToPQ = CPUTIME_AddToPQ+CPU2-CPU1
   WALLTIME_AddToPQ = WALLTIME_AddToPQ+ WALL2-WALL1
   iOrbital = iOrbital + PQ%P%p%nOrbitals(iAngmom)
ENDDO
!ENDDO

END SUBROUTINE Contract_P

SUBROUTINE AddToPQ(Integral,PQ,iAngmom,iOrbital,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: iAngmom,iOrbital,LUPRI,IPRINT
!
Integer :: nQ,nAng,nCont,endOrb,nP

nQ     = PQ%Q%p%totOrbitals
nP     = PQ%P%p%totOrbitals
nAng   = PQ%P%p%nOrbComp(iAngmom)
nCont  = PQ%P%p%nContracted(iAngmom)*PQ%P%p%nPasses
CALL AddToPQ1(Integral%PQ,Integral%IN,iOrbital,nP,nQ,nAng,nCont,LUPRI,IPRINT)
END SUBROUTINE AddToPQ

SUBROUTINE AddToPQ1(PQ,AddPQ,iOrbital,nP,nQ,nAng,nCont,LUPRI,IPRINT)
implicit none
Real(realk) :: PQ(nQ,nP),AddPQ(nAng,nQ,nCont)
Integer     :: iOrbital,nP,nQ,nAng,nCont,LUPRI,IPRINT
!
Integer     :: iQ,iAng,iCont,iP
!
iP = iOrbital
DO iCont=1,nCont
  DO iAng=1,nAng
    DO iQ=1,nQ
      PQ(iQ,iP) = AddPQ(iAng,iQ,iCont)
    ENDDO
    iP=iP+1
  ENDDO
ENDDO
CALL PrintPQ(PQ,nQ,nP,IOrbital,nAng,nCont,LUPRI,IPRINT)
END SUBROUTINE AddToPQ1

SUBROUTINE PrintPQ(PQ,nQ,nP,iOrbital,nAng,nCont,LUPRI,IPRINT)
implicit none
Integer     :: nQ,nAng,nCont,LUPRI,IPRINT,iorbital,nP
Real(realk) :: PQ(nQ,nP)
!
Integer     :: iQ,iAng,iCont
!
IF (IPRINT.GT.5) THEN
  CALL LSHEADER(LUPRI,'PQ contribution')
  DO iCont=1,nCont
    DO iAng=1,nAng
      WRITE(LUPRI,'(5X,A,I3,A,I3)') 'iCont =',iCont,' iAng =',iAng
      WRITE(LUPRI,'(5X,F16.9)') (PQ(iQ,iOrbital+iAng+(iCont-1)*nAng-1),iQ=1,nQ)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintPQ

SUBROUTINE DistributeHermiteP(Integral,PQ,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: LUPRI,IPRINT,iAngmom
!
Integer :: ntuvP,nP,nOrbQ,startOrbQ,endOrbQ
Integer :: startP,fullSP,endP,ioffP,fullOP,jP,tP,uP,vP,ituvP,ifullP,iOrbQ,iPrimP
Integer :: ifullp1,ifullp2,ifullp3,end2P

!
startP = 0
IF (PQ%P%p%hermiteSingle) startP = PQ%P%p%angmom(iAngmom)

fullSP   = PQ%P%p%startAngmom
endP     = PQ%P%p%angmom(iAngmom)
nP       = PQ%P%p%nPrimitives
ioffP    = startP*(startP+1)*(startP+2)/6
fullOP   = fullSP*(fullSP+1)*(fullSP+2)/6
ntuvP    = (endP+1)*(endP+2)*(endP+3)/6 - ioffP
nOrbQ    = PQ%Q%p%totOrbitals

Integral%nAng  = ntuvP
Integral%nPrim = np
Integral%nOrb  = nOrbQ
IF(PQ%kinetic)THEN
   CALL DistributeHermiteP_kinetic(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,fullSP,&
     &                             startP,endP,nP,ioffP,fullOP,PQ%P%p%nTUV,ntuvP,nOrbQ,LUPRI,IPRINT)
ELSE
   CALL DistributeHermiteP_regular(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,fullSP,&
     &                             startP,endP,nP,ioffP,fullOP,PQ%P%p%nTUV,ntuvP,nOrbQ,LUPRI,IPRINT)
ENDIF

IF (IPRINT.GT.20) THEN
   CALL PrintTensor(Integral%IN,'Primitive           ',&
        &ntuvP,nOrbQ,nP,Lupri,'iTUV  ','iOrb  ','Prim  ',3)
ENDIF

END SUBROUTINE DistributeHermiteP

!!$SUBROUTINE DistributeHermiteP_regular(IntegralIN,TUVQ,TUVindex,fullSP,startP,endP,nP,ioffP,&
!!$     &                                fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT)
!!$implicit none
!!$Real(realk)     :: IntegralIN(ntuvP,nOrbQ,nP)
!!$Real(realk)     :: TUVQ(ntuvFull,np,nOrbQ)
!!$Integer,pointer :: TUVindex(:,:,:)
!!$INTEGER         :: TUVvector(ntuvP)
!!$INTEGER         :: TUVvector1(ntuvP),TUVvector2(ntuvFULL)
!!$Integer         :: fullSP,startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
!!$!
!!$Integer :: jP,tP,uP,vP,ituvP,ifullP,iOrbQ,iPrimP,iTUVvector,nTUVvector,TMPTUV
!!$   
!!$IF((ioffP.EQ. 0) .AND. (fullOP .EQ. 0))THEN
!!$   iOrbQ=1
!!$   iPrimP=1
!!$   iTUVvector=1
!!$   DO jP = startP,endP
!!$      DO tP=jP,0,-1
!!$         DO uP=jP-tP,0,-1
!!$            vP=jP-tP-uP
!!$            ituvP  = TUVindex(tP,uP,vP)
!!$            TUVvector(iTUVvector) = ituvP
!!$            iTUVvector = iTUVvector+1
!!$            IntegralIN(ituvP,iOrbQ,iPrimP) = TUVQ(ituvP,iPrimP,iOrbQ)
!!$         ENDDO
!!$      ENDDO
!!$   ENDDO
!!$   nTUVvector=iTUVvector-1
!!$   DO iPrimP=2,nP
!!$      DO iTUVvector=1,nTUVvector
!!$         ituvP = TUVvector(iTUVvector)
!!$         IntegralIN(ituvP,iOrbQ,iPrimP) = TUVQ(ituvP,iPrimP,iOrbQ)
!!$      ENDDO
!!$   ENDDO
!!$   DO iOrbQ=2,nOrbQ
!!$      DO iPrimP=1,nP
!!$         DO iTUVvector=1,nTUVvector
!!$            ituvP = TUVvector(iTUVvector)
!!$            IntegralIN(ituvP,iOrbQ,iPrimP) = TUVQ(ituvP,iPrimP,iOrbQ)
!!$         ENDDO
!!$      ENDDO
!!$   ENDDO
!!$ELSE
!!$   iTUVvector = 1
!!$   iOrbQ=1
!!$   iPrimP=1
!!$   DO jP = startP,endP
!!$      DO tP=jP,0,-1
!!$         DO uP=jP-tP,0,-1
!!$            vP=jP-tP-uP
!!$            TMPTUV = TUVindex(tP,uP,vP)
!!$            ituvP  = TMPTUV-ioffP
!!$            TUVvector1(iTUVvector) = ituvP
!!$            ifullP = TMPTUV-fullOP
!!$            TUVvector2(iTUVvector) = ifullP
!!$            iTUVvector = iTUVvector +1
!!$            IntegralIN(ituvP,iOrbQ,iPrimP) = TUVQ(ifullP,iPrimP,iOrbQ)
!!$         ENDDO
!!$      ENDDO
!!$   ENDDO
!!$   nTUVvector=iTUVvector-1   
!!$   DO iPrimP=2,nP
!!$      DO iTUVvector=1,nTUVvector
!!$         ituvP = TUVvector1(iTUVvector) 
!!$         ifullP = TUVvector2(iTUVvector)
!!$         IntegralIN(ituvP,iOrbQ,iPrimP) = TUVQ(ifullP,iPrimP,iOrbQ)
!!$      ENDDO
!!$   ENDDO
!!$   DO iOrbQ=2,nOrbQ
!!$      DO iPrimP=1,nP
!!$         DO iTUVvector=1,nTUVvector
!!$            ituvP = TUVvector1(iTUVvector) 
!!$            ifullP = TUVvector2(iTUVvector)
!!$            IntegralIN(ituvP,iOrbQ,iPrimP) = TUVQ(ifullP,iPrimP,iOrbQ)
!!$         ENDDO
!!$      ENDDO
!!$   ENDDO
!!$ENDIF
!!$
!!$END SUBROUTINE DistributeHermiteP_regular
!!$
!!$SUBROUTINE DistributeHermiteP_kinetic(IntegralIN,TUVQ,TUVindex,fullSP,startP,endP,nP,ioffP,&
!!$     &                                fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT)
!!$implicit none
!!$Real(realk)     :: IntegralIN(ntuvP,nOrbQ,nP)
!!$Real(realk)     :: TUVQ(ntuvFull,np,nOrbQ)
!!$Integer,pointer :: TUVindex(:,:,:)
!!$Integer         :: fullSP,startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
!!$!
!!$Integer :: jP,tP,uP,vP,ituvP,iOrbQ,iPrimP
!!$Integer :: ifullp1,ifullp2,ifullp3
!!$Real(realk), parameter :: Half=0.5d0
!!$
!!$DO iOrbQ=1,nOrbQ
!!$   DO jP = startP,endP
!!$      DO tP=jP,0,-1
!!$         DO uP=jP-tP,0,-1
!!$            vP=jP-tP-uP
!!$            ituvP  = TUVindex(tP,uP,vP)-ioffP
!!$            ifullP1 = TUVindex(tP+2,uP,vP)-fullOP
!!$            ifullP2 = TUVindex(tP,uP+2,vP)-fullOP
!!$            ifullP3 = TUVindex(tP,uP,vP+2)-fullOP             
!!$            DO iPrimP=1,nP
!!$               IntegralIN(ituvP,iOrbQ,iPrimP) = &
!!$                    &- HALF*TUVQ(ifullP1,iPrimP,iOrbQ) &
!!$                    &- HALF*TUVQ(ifullP2,iPrimP,iOrbQ) &
!!$                    &- HALF*TUVQ(ifullP3,iPrimP,iOrbQ)
!!$            ENDDO
!!$         ENDDO
!!$      ENDDO
!!$   ENDDO
!!$ENDDO
!!$END SUBROUTINE DistributeHermiteP_kinetic

SUBROUTINE DistributeHermiteP_regular(IntegralIN,TUVQ,TUVindex,fullSP,startP,endP,nP,ioffP,&
     &                                fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT)
implicit none
Real(realk)     :: IntegralIN(ntuvP,nOrbQ,nP)
Real(realk)     :: TUVQ(np,ntuvFull,nOrbQ)
Integer,pointer :: TUVindex(:,:,:)
Integer         :: fullSP,startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
!
Integer :: jP,tP,uP,vP,ituvP,ifullP,iOrbQ,iPrimP
   
DO iOrbQ=1,nOrbQ
   DO jP = startP,endP
      DO tP=jP,0,-1
         DO uP=jP-tP,0,-1
            vP=jP-tP-uP
            ituvP  = TUVindex(tP,uP,vP)-ioffP
            ifullP = TUVindex(tP,uP,vP)-fullOP
            DO iPrimP=1,nP
               IntegralIN(ituvP,iOrbQ,iPrimP) = TUVQ(iPrimP,ifullP,iOrbQ)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE DistributeHermiteP_regular

SUBROUTINE DistributeHermiteP_kinetic(IntegralIN,TUVQ,TUVindex,fullSP,startP,endP,nP,ioffP,&
     &                                fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT)
implicit none
Real(realk)     :: IntegralIN(ntuvP,nOrbQ,nP)
Real(realk)     :: TUVQ(np,ntuvFull,nOrbQ)
Integer,pointer :: TUVindex(:,:,:)
Integer         :: fullSP,startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
!
Integer :: jP,tP,uP,vP,ituvP,iOrbQ,iPrimP
Integer :: ifullp1,ifullp2,ifullp3
Real(realk), parameter :: Half=0.5d0

DO iOrbQ=1,nOrbQ
   DO jP = startP,endP
      DO tP=jP,0,-1
         DO uP=jP-tP,0,-1
            vP=jP-tP-uP
            ituvP  = TUVindex(tP,uP,vP)-ioffP
            ifullP1 = TUVindex(tP+2,uP,vP)-fullOP
            ifullP2 = TUVindex(tP,uP+2,vP)-fullOP
            ifullP3 = TUVindex(tP,uP,vP+2)-fullOP             
            DO iPrimP=1,nP
               IntegralIN(ituvP,iOrbQ,iPrimP) = &
                    &- HALF*TUVQ(iPrimP,ifullP1,iOrbQ) &
                    &- HALF*TUVQ(iPrimP,ifullP2,iOrbQ) &
                    &- HALF*TUVQ(iPrimP,ifullP3,iOrbQ)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE DistributeHermiteP_kinetic

SUBROUTINE Contract_Q(INTEGRAL,PQ,Input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)  :: INTEGRAL
TYPE(Integrand)     :: PQ
TYPE(IntegralInput) :: Input
Integer             :: LUPRI,IPRINT,LWORK,WORKLENGTH
REAL(REALK)         :: WORK(WORKLENGTH)
!
Integer  :: iAngmom,startMomQ,iOrbital,ideriv,nderiv
!
REAL(REALK)          :: CPU1,CPU2,WALL1,WALL2

nderiv=Input%NderivQ
!NULLIFY(Integral%TUVQ)
!ALLOCATE(Integral%TUVQ(PQ%P%p%nPrimitives,PQ%P%p%nTUV,PQ%Q%p%totOrbitals))
!
iOrbital = 1
!
!Loop over angular contributions sharing the set of primitive functions
DO iDeriv=1,nderiv
   DO iAngmom=1,PQ%Q%p%nAngmom
      ! Distribute Hermite 2-el integrals WTUV to <tuv1|w|tuv2>

      CALL GETTIM(CPU1,WALL1)
      CALL DistributeHermiteQ(Integral,PQ,iAngmom,iDeriv,LUPRI,IPRINT)
      CALL GETTIM(CPU2,WALL2)
      CPUTIME_DistributeHermiteQ = CPUTIME_DistributeHermiteQ+CPU2-CPU1
      WALLTIME_DistributeHermiteQ = WALLTIME_DistributeHermiteQ+ WALL2-WALL1
      CALL ContractEcoeff(Integral,PQ%Q%p,-1.d0,iAngmom,LUPRI,IPRINT)
      CALL GETTIM(CPU1,WALL1)
      CPUTIME_ContractEcoeffQ = CPUTIME_ContractEcoeffQ+CPU1-CPU2
      WALLTIME_ContractEcoeffQ = WALLTIME_ContractEcoeffQ+ WALL1-WALL2
      CALL ContractBasis(Integral,PQ%Q%p,iAngmom,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
      CALL GETTIM(CPU2,WALL2)
      CPUTIME_ContractBasisQ = CPUTIME_ContractBasisQ+CPU2-CPU1
      WALLTIME_ContractBasisQ = WALLTIME_ContractBasisQ+ WALL2-WALL1
      CALL SphericalTransform(Integral,PQ%Q%p,iAngmom,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
      CALL GETTIM(CPU1,WALL1)
      CPUTIME_SphericalTransformQ = CPUTIME_SphericalTransformQ+CPU1-CPU2
      WALLTIME_SphericalTransformQ = WALLTIME_SphericalTransformQ+ WALL1-WALL2
      CALL AddToTUVQ(Integral,PQ,iAngmom,iOrbital,LUPRI,IPRINT)
      CALL GETTIM(CPU2,WALL2)
      CPUTIME_AddToTUVQ = CPUTIME_AddToTUVQ+CPU2-CPU1
      WALLTIME_AddToTUVQ = WALLTIME_AddToTUVQ+ WALL2-WALL1
      iOrbital = iOrbital + PQ%Q%p%nOrbitals(iAngmom)
   ENDDO
ENDDO
    
END SUBROUTINE Contract_Q

SUBROUTINE AddToTUVQ(Integral,PQ,iAngmom,iOrb,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: iAngmom,iOrb,LUPRI,IPRINT 
!
Integer :: nPrimP,nTUVP,nCompQ,nContQ,endOrb,totOrbQ

nPrimP  = PQ%P%p%nPrimitives*PQ%P%p%nPasses
nTUVP   = PQ%P%p%nTUV
nCompQ  = PQ%Q%p%nOrbComp(iAngmom)
nContQ  = PQ%Q%p%nContracted(iAngmom)*PQ%Q%p%nPasses
endOrb  = iOrb + nCompQ*nContQ - 1
totOrbQ = PQ%Q%p%totOrbitals

CALL AddToTUVQ0(Integral%integrals,Integral%IN,totOrbQ,iOrb,endOrb,nTUVP,nPrimP,&
     &          nCompQ,nContQ,LUPRI,IPRINT)

END SUBROUTINE AddToTUVQ

!!$SUBROUTINE AddToTUVQ0(TUVQ,OldTUVQ,totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
!!$implicit none
!!$Integer     :: totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
!!$Real(realk) :: TUVQ(nTUVP,nPrimP,totOrbQ)
!!$Real(realk) :: OldTUVQ(nCompQ,nTUVP,nPrimP,nContQ)
!!$
!!$
!!$CALL AddToTUVQ1(TUVQ(:,:,startOrb:endOrb),OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
!!$
!!$END SUBROUTINE AddToTUVQ0
!!$
!!$SUBROUTINE AddToTUVQ1(TUVQ,OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
!!$implicit none
!!$Real(realk) :: TUVQ(nTUVP,nPrimP,nCompQ,nContQ)
!!$Real(realk) :: OldTUVQ(nCompQ,nTUVP,nPrimP,nContQ)
!!$Integer     :: nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
!!$!
!!$Integer :: iTUVP,iPrimP,iCompQ,iContQ
!!$!
!!$DO iContQ=1,nContQ
!!$   DO iPrimP=1,nPrimP
!!$      DO iCompQ=1,nCompQ
!!$         DO iTUVP=1,nTUVP
!!$            TUVQ(iTUVP,iPrimP,iCompQ,iContQ) = OldTUVQ(iCompQ,iTUVP,iPrimP,iContQ)
!!$         ENDDO
!!$      ENDDO
!!$   ENDDO
!!$!  CALL transposition(TUVQ(:,:,:,iContQ),OldTUVQ(:,:,:,iContQ),nPrimP*nTUVP,nCompQ)
!!$ENDDO
!!$
!!$CALL PrintTUVQ(TUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
!!$END SUBROUTINE AddToTUVQ1
SUBROUTINE AddToTUVQ0(TUVQ,OldTUVQ,totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
implicit none
Integer     :: totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
Real(realk) :: TUVQ(nPrimP,nTUVP,totOrbQ)
Real(realk) :: OldTUVQ(nCompQ,nTUVP,nPrimP,nContQ)


CALL AddToTUVQ1(TUVQ(:,:,startOrb:endOrb),OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)

END SUBROUTINE AddToTUVQ0

SUBROUTINE AddToTUVQ1(TUVQ,OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
implicit none
Real(realk) :: TUVQ(nPrimP,nTUVP,nCompQ,nContQ)
Real(realk) :: OldTUVQ(nCompQ,nTUVP,nPrimP,nContQ)
Integer     :: nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
!
Integer :: iTUVP,iPrimP,iCompQ,iContQ
!
DO iContQ=1,nContQ
  DO iCompQ=1,nCompQ
    DO iTUVP=1,nTUVP
      DO iPrimP=1,nPrimP
        TUVQ(iPrimP,iTUVP,iCompQ,iContQ) = OldTUVQ(iCompQ,iTUVP,iPrimP,iContQ)
      ENDDO
    ENDDO
  ENDDO
!  CALL transposition(TUVQ(:,:,:,iContQ),OldTUVQ(:,:,:,iContQ),nPrimP*nTUVP,nCompQ)
ENDDO

CALL PrintTUVQ(TUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
END SUBROUTINE AddToTUVQ1

!!$SUBROUTINE PrintTUVQ(TUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
!!$implicit none
!!$Real(realk) :: TUVQ(nTUVP,nPrimP,nCompQ,nContQ)
!!$Integer     :: nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
!!$!
!!$Integer :: iTUVP,iPrimP,iCompQ,iContQ
!!$!
!!$IF (IPRINT.GT.10) THEN
!!$  CALL LSHEADER(LUPRI,'TUVQ')
!!$  DO iContQ=1,nContQ
!!$    DO iCompQ=1,nCompQ
!!$      DO iPrimP=1,nPrimP
!!$        WRITE(LUPRI,'(5X,A,I3,A,I3,A,I3)') 'iPrimP =',iPrimP,' iCompQ =',iCompQ,' iContQ =',iContQ
!!$        WRITE(LUPRI,'(5X,5F10.4)')  (TUVQ(iTUVP,iPrimP,iCompQ,iContQ), iTUVP=1,nTUVP)
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDIF
!!$END SUBROUTINE PrintTUVQ

SUBROUTINE PrintTUVQ(TUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
implicit none
Real(realk) :: TUVQ(nPrimP,nTUVP,nCompQ,nContQ)
Integer     :: nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
!
Integer :: iTUVP,iPrimP,iCompQ,iContQ
!
IF (IPRINT.GT.10) THEN
  CALL LSHEADER(LUPRI,'TUVQ')
  DO iContQ=1,nContQ
    DO iCompQ=1,nCompQ
      DO iTUVP=1,nTUVP
        WRITE(LUPRI,'(5X,A,I3,A,I3,A,I3)') 'iTUVP =',iTUVP,' iCompQ =',iCompQ,' iContQ =',iContQ
        WRITE(LUPRI,'(5X,5F10.4)')  (TUVQ(iPrimP,iTUVP,iCompQ,iContQ), iPrimP=1,nPrimP)
      ENDDO
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintTUVQ

SUBROUTINE SphericalTransform(Integral,P,iAngmom,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,iAngmom,LWORK,WORKLENGTH
!
Real(realk)         :: WORK(WORKLENGTH)
!Real(realk),pointer :: Spherical(:,:)
Integer             :: lm1,lm2,ijk1,ijk2,iA1,iA2,ang1,ang2,lm,ijk
Integer             :: dim1
Logical             :: Sph1,Sph2
Real(realk)         :: D1=1.0d0,D0=0.0d0
!
iA1  = P%indexAng1(iangmom)
iA2  = P%indexAng2(iangmom)
ang1 = P%orbital1%angmom(iA1)
ang2 = P%orbital2%angmom(iA2)
Sph1 = P%orbital1%spherical.AND.(ang1.GT.1)
Sph2 = P%orbital2%spherical.AND.(ang2.GT.1)
IF (Sph1.OR.Sph2) THEN
  ijk1 = (ang1+1)*(ang1+2)/2 
  ijk2 = (ang2+1)*(ang2+2)/2
  ijk  = ijk1 * ijk2
  lm1 = ijk1
  IF (Sph1) lm1 = 2*ang1+1 
  lm2 = ijk2
  IF (Sph2) lm2 = 2*ang2+1
  lm  = lm1 * lm2
  
!  NULLIFY(Spherical)
!  ALLOCATE(Spherical(ijk1*ijk2,lm1*lm2))
   IF(LWORK+ijk1*ijk2*lm1*lm2-1 .GT. WORKLENGTH)THEN
      print*,'LWORK                ',LWORK
      print*,'ijk1*ijk2*lm1*lm2    ',ijk1*ijk2*lm1*lm2
      CALL QUIT('MEM in SphericalTransform')
   ENDIF

  CALL ContructSphericalTransformation1(WORK(LWORK:LWORK+ijk1*ijk2*lm1*lm2-1),&
       & ijk1,ijk2,lm1,lm2,ang1,ang2,P,lm,integral,lupri,iprint)

!  DEALLOCATE(Spherical)
!  NULLIFY(Spherical)
  Integral%nAng = lm
  call swapRealPointers(Integral%IN,Integral%OUT)
ENDIF
END SUBROUTINE SphericalTransform

SUBROUTINE ContructSphericalTransformation1(Spherical,ijk1,ijk2,lm1,lm2,ang1,ang2,P,lm,Integral,lupri,iprint)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,iAngmom
!
Real(realk)         :: Spherical(ijk1*ijk2,lm1*lm2)
Integer             :: lm1,lm2,ijk1,ijk2,iA1,iA2,ang1,ang2,lm,ijk
Integer             :: dim1
Logical             :: Sph1,Sph2
Real(realk)         :: D1=1.0d0,D0=0.0d0

! Sets up spherical transformation matrix
  CALL ContructSphericalTransformation(Spherical,P%orbital1%SPH_MAT(ang1+1)%p%elms,&
     &              P%orbital2%SPH_MAT(ang2+1)%p%elms,lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)

  
! Spherical transformation (with both centers beloning to one electron simultaneously)
  dim1 = Integral%nOrb*Integral%nPrim
  CALL DGEMM('T','N',lm,dim1,ijk1*ijk2,D1,Spherical,ijk1*ijk2,&
     &     Integral%IN,ijk1*ijk2,D0,Integral%OUT,lm)

END SUBROUTINE CONTRUCTSPHERICALTRANSFORMATION1


SUBROUTINE ContructSphericalTransformation(Spherical,Spher1,Spher2,lm1,lm2,&
     &                                     ijk1,ijk2,LUPRI,IPRINT)
Real(realk) :: Spherical(ijk1,ijk2,lm1,lm2)
Real(realk) :: Spher1(lm1,ijk1)
Real(realk) :: Spher2(lm2,ijk2)
Integer     :: ijk1,ijk2,lm1,lm2
!
Integer     :: indijk1,indijk2,indlm1,indlm2
!

DO indijk2=1,ijk2
   DO indijk1=1,ijk1
      DO indlm2=1,lm2
         DO indlm1=1,lm1
            Spherical(indijk1,indijk2,indlm1,indlm2) = Spher1(indlm1,indijk1)*Spher2(indlm2,indijk2)

         ENDDO
      ENDDO
   ENDDO
ENDDO
CALL PrintSphericalTransformation(Spherical,lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)
END SUBROUTINE ContructSphericalTransformation

SUBROUTINE PrintSphericalTransformation(Spherical,lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)
Real(realk) :: Spherical(ijk1,ijk2,lm1,lm2)
Integer     :: ijk1,ijk2,lm1,lm2
!
Integer     :: indijk1,indijk2,indlm1,indlm2
!
IF (IPRINT.GT.50) THEN
  CALL LSHEADER(LUPRI,'SphericalTransformation')
  DO indlm2=1,lm2
    DO indlm1=1,lm1
      WRITE(LUPRI,'(5X,A,I3,A,I3)') 'lm1 =',indlm1,' lm2 =',indlm2
      DO indijk1=1,ijk1
        WRITE(LUPRI,'(5X,5F10.4)') &
       &       (Spherical(indijk1,indijk2,indlm1,indlm2),indijk2=1,ijk2)
      ENDDO
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintSphericalTransformation

SUBROUTINE ContractBasis(Integral,P,iAngmom,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,iAngmom,LWORK,WORKLENGTH
!
Integer             :: iA1,iA2,nCont,nC1,nC2,dim1
Real(realk)         :: WORK(WORKLENGTH)
Real(realk),pointer :: CC(:,:,:)
Real(realk)         :: D1=1.0d0,D0=0.0d0
iA1 = P%indexAng1(iangmom)
iA2 = P%indexAng2(iangmom)
nC1 = P%orbital1%nContracted(iA1)
nC2 = P%orbital2%nContracted(iA2)
nCont = nC1*nC2
dim1 = Integral%nAng*Integral%nOrb

IF(LWORK+P%nPrimitives*nC1*nC2-1 .GT. WORKLENGTH)THEN
   print*,'LWORK                ',LWORK
   print*,'P%nPrimitives        ',P%nPrimitives
   print*,'nC1*nC2              ',nC1*nC2
   print*,'P%nPrimitives*nC1*nC2',P%nPrimitives*nC1*nC2
   CALL QUIT('MEM in ContractBasis')
ENDIF
!NULLIFY(CC)
!ALLOCATE(CC(P%nPrimitives,nC1,nC2))

CALL ContractBasis1(Integral%IN,Integral%OUT,&
     &      WORK(LWORK:LWORK+P%nPrimitives*nC1*nC2-1),P,P%nPrimitives,&
     &      nCont,P%nPasses,dim1,nC1,nC2,iA1,iA2,LUPRI,IPRINT)

!DEALLOCATE(CC)

IF (IPRINT.GT.50) THEN
   CALL PrintTensor(Integral%OUT,'Contracted          ',&
        &Integral%nAng,Integral%nOrb,nCont*P%nPasses,Lupri,'iAngP ','iOrb  ','iCont ',3)
ENDIF

CALL swapRealPointers(Integral%IN,Integral%OUT)
Integral%nPrim = nCont*P%nPasses

END SUBROUTINE ContractBasis

SUBROUTINE ContractBasis1(PrimInt,ContInt,CC,P,nPrim,nCont,nPasses,nDim,&
    &                     nC1,nC2,iA1,iA2,LUPRI,IPRINT)
implicit none
Integer       :: nPrim,nCont,nPasses,nDim,nC1,nC2,iA1,iA2,LUPRI,IPRINT
Real(realk)   :: PrimInt(nDim,nPrim,nPasses),ContInt(nDim,nCont,nPasses)
Real(realk)   :: CC(nPrim,nC1,nC2)
TYPE(Overlap) :: P
!
Integer     :: iPass
Real(realk) :: D1=1.0d0,D0=0.0d0

iPass = 1
  CALL ConstructContraction(CC,P,iPass,nPrim,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
DO iPass=1,nPasses
  CALL DGEMM('N','N',nDim,nCont,nPrim,D1,PrimInt(1,1,iPass),nDim,&
       &     CC,nPrim,D0,ContInt(1,1,iPass),nDim)
ENDDO

END SUBROUTINE ContractBasis1

SUBROUTINE PrintContracted(Contracted,nAng,nOrb,nPrim,LUPRI,IPRINT)
Real(realk) :: Contracted(nAng,nOrb,nPrim)
Integer     :: nAng,nOrb,nPrim
!
Integer :: iAng,iOrb,iPrim
IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(5X,A)') '**********************************************************'
  WRITE(LUPRI,'(5X,A)') '***                   Contracted'
  WRITE(LUPRI,'(5X,A)') '**********************************************************'
  DO iAng=1,nAng
    DO iOrb=1,nOrb
      WRITE(LUPRI,'(5X,A,I3,A,I3)') 'iAngP =',iAng,' iOrb =',iOrb
      WRITE(LUPRI,'(17X,5F10.4)') (Contracted(iAng,iOrb,iPrim),iPrim=1,nPrim)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintContracted

SUBROUTINE ConstructContraction(CC,P,iPass,nP,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
Implicit none
Type(Overlap) :: P
Integer       :: iPass,iA1,iA2,nP,nC1,nC2,LUPRI,IPRINT
Real(realk)   :: CC(nP,nC1,nC2)
!
Integer       :: iP,i1,i2,iC1,iC2,nP1,nP2

nP1 = P%orbital1%nPrimitives
nP2 = P%orbital2%nPrimitives
IF(P%sameAO)THEN
   DO iC2=1,nC2
      DO iC1=1,nC1
         DO iP=1,nP
            i1 = P%iprim1(iP+(iPass-1)*nP)
            i2 = P%iprim2(iP+(iPass-1)*nP)
            IF(i1 .NE. i2)THEN
!            CC(iP,iC1,iC2)=P%orbital1%CC(i1,iC1,iA1)*P%orbital2%CC(i2,iC2,iA2)&
!                 &        +P%orbital1%CC(i2,iC1,iA1)*P%orbital2%CC(i1,iC2,iA2)
            CC(iP,iC1,iC2)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2)&
                 &        +P%orbital1%CC(iA1)%p%elms(i2+(iC1-1)*nP2)*P%orbital2%CC(iA2)%p%elms(i1+(iC2-1)*nP1)


            ELSE
!            CC(iP,iC1,iC2)=P%orbital1%CC(i1,iC1,iA1)*P%orbital2%CC(i2,iC2,iA2)
            CC(iP,iC1,iC2)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ELSE
   DO iC2=1,nC2
      DO iC1=1,nC1
         DO iP=1,nP
            i1 = P%iprim1(iP+(iPass-1)*nP)
            i2 = P%iprim2(iP+(iPass-1)*nP)
!            CC(iP,iC1,iC2)=P%orbital1%CC(i1,iC1,iA1)*P%orbital2%CC(i2,iC2,iA2) 
            CC(iP,iC1,iC2)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2) 
         ENDDO
      ENDDO
   ENDDO
ENDIF

IF (IPRINT.GT.30) THEN
   CALL PrintTensor(CC,'Overlap CC          ',nP,nC1,nC2,Lupri,&
        & 'prim  ','C1    ','C2    ',1)
ENDIF

END SUBROUTINE ConstructContraction

SUBROUTINE PrintOverlapContraction(CC,nP,nC1,nC2,LUPRI,IPRINT)
Implicit none
Integer       :: nP,nC1,nC2,LUPRI,IPRINT
Real(realk)   :: CC(nP,nC1,nC2)
!
Integer       :: iP,iC1,iC2
IF (IPRINT.GT.30) THEN
  WRITE(LUPRI,'(5X,A)') '***************************************************'
  WRITE(LUPRI,'(5X,A)') '***                  Overlap CC'
  WRITE(LUPRI,'(5X,A)') '***************************************************'
  WRITE(LUPRI,'(7X,A)') 'C1  C2 CC(prim12)'
  DO iC2=1,nC2
    DO iC1=1,nC1
      WRITE(LUPRI,'(5X,2I4,5F10.4/,(13X,5F10.4))') iC1,iC2,(CC(iP,iC1,iC2),iP=1,nP)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintOverlapContraction

SUBROUTINE ContractEcoeff(Integral,P,signP,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)  :: Integral
TYPE(Overlap)       :: P
Integer             :: LUPRI,IPRINT,iAngmom
Real(realk)         :: signP
!
Real(realk),pointer :: Ecoeffs(:)
Integer             :: l1,l2,ijk1,ijk2,ijk,iPrimP
Real(realk)         :: D1=1.0d0,D0=0.0d0
!
IF (P%hermiteSingle) THEN
  CALL SingleHermiteEcoeff(Integral%In,P,signP,iAngmom,Integral%nAng,Integral%nOrb,LUPRI,IPRINT)
ELSE
!  CALL QUIT('Only Hermite implemented in ContractEcoeff')
  l1 = P%orbital1%angmom(P%indexAng1(iangmom))
  l2 = P%orbital2%angmom(P%indexAng2(iangmom))
  ijk1 = (l1+1)*(l1+2)/2
  ijk2 = (l2+1)*(l2+2)/2
  ijk  = ijk1*ijk2
  IF (P%ETUVisSet) THEN
    CALL ContractEcoeff1(Integral%IN,Integral%OUT,P%ETUV,P%ETUVindex(iAngmom), &
       &                 Integral%nOrb,Integral%nAng,ijk,Integral%nPrim)
  ELSE
    NULLIFY(Ecoeffs)
    ALLOCATE(Ecoeffs(ijk*Integral%nAng*Integral%nPrim))
    CALL BuildEcoeffTensor2(integral%TUV,P,signP,Ecoeffs,ijk,Integral%nAng,&
         &Integral%nPrim,iAngmom,P%nPasses,LUPRI,IPRINT)
    CALL ContractEcoeff1(Integral%IN,Integral%OUT,Ecoeffs,1, &
         &                 Integral%nOrb,Integral%nAng,ijk,Integral%nPrim)
    DEALLOCATE(Ecoeffs)
  ENDIF
  Integral%nAng = ijk
  CALL swapRealPointers(Integral%IN,Integral%OUT)
ENDIF

IF (IPRINT.GT.50) THEN
   CALL PrintTensor(Integral%IN,'ContractEcoeff      ',&
        &Integral%nAng,Integral%nOrb,Integral%nPrim,Lupri,&
        &'ijk   ','iOrb  ','Prim  ',3)
ENDIF

END SUBROUTINE ContractEcoeff

SUBROUTINE ContractEcoeff1(IntegralIN,IntegralOUT,Ecoeffs,startE,nOrb,nAng,ijk,nPrim)
implicit none
Integer             :: nOrb,nAng,ijk,nPrim,startE
Real(realk)         :: IntegralIN(nAng,nOrb,nPrim)
Real(realk)         :: IntegralOUT(ijk,nOrb,nPrim)
Real(realk),pointer :: Ecoeffs(:)
!
Integer :: iPrimP,iE
Real(realk)         :: D1=1.0d0,D0=0.0d0

IF (ijk.EQ.1) THEN
   CALL ContractEcoeff1ss(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nPrim)
ELSE
   iE = startE
   DO iPrimP = 1,nPrim
      CALL DGEMM('N','N',ijk,nOrb,nAng,D1,Ecoeffs(iE),ijk,&
           &          IntegralIN(1,1,iPrimP),nAng,D0,IntegralOUT(1,1,iPrimP),ijk)
      iE = iE+nAng*ijk
   ENDDO
ENDIF

END SUBROUTINE ContractEcoeff1

SUBROUTINE ContractEcoeff1ss(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nPrim)
implicit none
Integer     :: nOrb,nPrim
Real(realk) :: IntegralIN(nOrb,nPrim)
Real(realk) :: IntegralOUT(nOrb,nPrim)
Real(realk) :: Ecoeffs(nPrim)
!
Integer :: iPrim,iOrb
Real(realk)         :: D1=1.0d0,D0=0.0d0,E

DO iPrim = 1,nPrim
  E = Ecoeffs(iPrim)
  DO iOrb=1,nOrb
    IntegralOUT(iOrb,iPrim) = E*IntegralIN(iOrb,iPrim)
  ENDDO
ENDDO

END SUBROUTINE ContractEcoeff1ss

SUBROUTINE ContractFTUV(Integral,P,Q,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P,Q
Integer            :: LUPRI,IPRINT
IF (P%type.EQ.'Hermite-single') THEN
  CALL QUIT('Jengine not implemented with hermite-single yet')
!  CALL SingleHermiteEcoeff(Integral%IntegralIn,P,iAngmom,Integral%nAng,Integral%nOrb,LUPRI,IPRINT)
ELSE
  CALL ContractFTUV1(Integral%integrals,Integral%RTUV,Q%FTUV,Integral%TUV%TUVindex,&
     &           Integral%TUV%Tindex,Integral%TUV%Uindex,Integral%TUV%Vindex,&
     &           Integral%ntuv,Integral%nPrim,P%nTUV,Q%nTUV,P%nPrimitives,Q%nPrimitives,1,lupri)
ENDIF

Integral%nAng  = 1
Integral%nPrim = 1
Integral%nOrb  = P%nPrimitives*P%nTUV
CALL PrintTUVQ(Integral%integrals,P%nTUV,P%nPrimitives,1,1,LUPRI,IPRINT)

END SUBROUTINE ContractFTUV

SUBROUTINE ContractFTUV1(TUV,WTUV,FTUV,TUVindex,Tindex,Uindex,Vindex,&
     &                   ntuvPQ,nPrimPQ,ntuvP,ntuvQ,nPrimP,nPrimQ,nDmat,lupri)
implicit none
Integer             :: ntuvP,ntuvQ,nPrimP,nPrimQ,ntuvPQ,nPrimPQ,nDmat
Real(realk)         :: TUV(nPrimP,ntuvP,nDmat)!nprimP,ntuvP,ndmat
Real(realk),pointer :: FTUV(:,:,:)
Real(realk)         :: WTUV(nTUVPQ,nPrimPQ)
Integer,pointer     :: TUVindex(:,:,:),Tindex(:),Uindex(:),Vindex(:)
!
Integer     :: idmat,tuvP,tP,uP,vP,tuvQ,tQ,uQ,vQ,jQ,tuvPQ
Integer     :: iPrimP,iPrimQ,iPrimPQ,lupri
Real(realk) :: D1=1.0d0,DM1=-1.0d0,signQ,fsum
!Simen
Real(realk) :: newfsum(nPrimP)

!CALL DZERO(TUV,nPrimP*ntuvP*nDmat)
! replace by loop
idmat = 1
DO tuvP=1,ntuvP
  tP = Tindex(tuvP)
  uP = Uindex(tuvP)
  vP = Vindex(tuvP)
  DO tuvQ=1,ntuvQ
    tQ = Tindex(tuvQ)
    uQ = Uindex(tuvQ)
    vQ = Vindex(tuvQ)
    jQ = tQ+uQ+vQ
    signQ = D1
    IF (mod(jQ,2).EQ.1) signQ = DM1
    tuvPQ = TUVindex(tP+tQ,uP+uQ,vP+vQ)
!simen
    IF (.TRUE.) THEN
       DO iPrimP=1,nPrimP
          newfsum(iPrimP) = WTUV(tuvPQ,iPrimP)*FTUV(1,tuvQ,idmat)
       ENDDO
       iPrimPQ = nPrimP
       DO iPrimQ=2,nPrimQ
          DO iPrimP=1,nPrimP
             iPrimPQ = iPrimPQ + 1
             newfsum(iPrimP) = newfsum(iPrimP) + WTUV(tuvPQ,iPrimPQ)*FTUV(iPrimQ,tuvQ,idmat)
          ENDDO
       ENDDO
       DO iPrimP=1,nPrimP
          TUV(iPrimP,tuvP,idmat) = TUV(iPrimP,tuvP,idmat) + newfsum(iPrimP)!*signQ         
       ENDDO
!simen orig
    ELSE
    iPrimPQ = 1
    DO iPrimP=1,nPrimP
      fsum = WTUV(tuvPQ,iPrimPQ)*FTUV(1,tuvQ,idmat)
      iPrimPQ = iPrimPQ+1
      DO iPrimQ=2,nPrimQ
        fsum = fsum + WTUV(tuvPQ,iPrimPQ)*FTUV(iPrimQ,tuvQ,idmat)
        iPrimPQ = iPrimPQ+1
      ENDDO
      TUV(iPrimP,tuvP,idmat) = TUV(iPrimP,tuvP,idmat) + fsum!*signQ
    ENDDO
    ENDIF
!simen
  ENDDO
ENDDO

END SUBROUTINE ContractFTUV1

SUBROUTINE BuildEcoeffTensor1(TUV,Q,signQ,Ecoeffs,x,y,z,iAngmom,nPasses,LUPRI,IPRINT)
!Ordering of Ecoeffs(nAng,ijk,nprim)
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: Q
Integer            :: iAngmom,LUPRI,IPRINT,nPasses
Integer            :: ijk
Real(realk)        :: signQ
!Allocatable for ECOEFFS - Dimensions are nPrimitives, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
real(realk),allocatable  :: ETIJ(:,:,:,:,:) 
integer             :: startQ,endQ,nQ,ioffQ,ntuvQ,ijk1,ijk2,jQ,Q1,Q2,tQ,uQ,vQ
integer             :: ituvQ,iQ1,jQ1,kQ1,iQ2,jQ2,kQ2,iPrimQ
integer             :: iorb1,iorb2,iang,iprim,x,y,z,l1,l2
Real(realk)         :: Ecoeffs(x,y,z) !nAng,ijk1*ijk2,nprim
Real(realk)         :: signijk

l1 = Q%orbital1%angmom(Q%indexAng1(iangmom))
l2 = Q%orbital2%angmom(Q%indexAng2(iangmom))

!ALLOCATE(ETIJ(Q%nPrimitives,0:l1+l2,0:l1,0:l2,3))
ALLOCATE(ETIJ(z,0:l1+l2,0:l1,0:l2,3))

CALL GET_ECOEFF(ETIJ,z,l1+l2,l1,l2,Q%nPrimitives,nPasses,Q,LUPRI,IPRINT)

CALL DZERO(Ecoeffs,x*y*z)

ijk=0
DO Q2 = 0,l2
   DO iQ2=Q2,0,-1
      DO jQ2=Q2-iQ2,0,-1
         kQ2=Q2-iQ2-jQ2
         DO Q1 = 0,l1
            DO iQ1=Q1,0,-1
               DO jQ1=Q1-iQ1,0,-1
                  kQ1=Q1-iQ1-jQ1
                  IF(iQ1+jQ1+kQ1+iQ2+jQ2+kQ2 .EQ. l1+l2) then
                  ijk=ijk+1
                     DO tQ=0,iQ1+iQ2
                        DO uQ=0,jQ1+jQ2
                           DO vQ=0,kQ1+kQ2
                              signijk = signQ**(tQ+uQ+vQ)
                              ituvQ=TUV%TUVindex(tQ,uQ,vQ)
                              DO iPrimQ=1,z!Q%nPrimitives
                                 Ecoeffs(ituvQ,ijk,iPrimQ) = &
                                      &ETIJ(iPrimQ,tQ,iQ1,iQ2,1)&
                                      &*ETIJ(iPrimQ,uQ,jQ1,jQ2,2)&
                                      &*ETIJ(iPrimQ,vQ,kQ1,kQ2,3)&
                                      &*Q%preExpFac(iPrimQ)*signijk
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDo
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

DEALLOCATE(ETIJ)

END SUBROUTINE BuildEcoeffTensor1

SUBROUTINE BuildEcoeffTensor2(TUV,Q,signQ,Ecoeffs,x,y,z,iAngmom,nPasses,LUPRI,IPRINT)
!different ordering of Ecoeffs(ijk,nAng,nprim)
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: Q
Integer            :: iAngmom,LUPRI,IPRINT,nPasses
Integer            :: ijk
Real(realk)        :: signQ
!Allocatable for ECOEFFS - Dimensions are nPrimitives, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
real(realk),allocatable  :: ETIJ(:,:,:,:,:) 
integer             :: startQ,endQ,nQ,ioffQ,ntuvQ,ijk1,ijk2,jQ,Q1,Q2,tQ,uQ,vQ
integer             :: ituvQ,iQ1,jQ1,kQ1,iQ2,jQ2,kQ2,iPrimQ
integer             :: iorb1,iorb2,iang,iprim,x,y,z,l1,l2
Real(realk)         :: Ecoeffs(x,y,z) !ijk1*ijk2,nAng,nprim
Real(realk)         :: signijk

l1 = Q%orbital1%angmom(Q%indexAng1(iangmom))
l2 = Q%orbital2%angmom(Q%indexAng2(iangmom))

!ALLOCATE(ETIJ(Q%nPrimitives,0:l1+l2,0:l1,0:l2,3))
ALLOCATE(ETIJ(z,0:l1+l2,0:l1,0:l2,3))

CALL GET_ECOEFF(ETIJ,z,l1+l2,l1,l2,Q%nPrimitives,nPasses,Q,LUPRI,IPRINT)

CALL DZERO(Ecoeffs,x*y*z)

ijk=0
DO Q2 = 0,l2
   DO iQ2=Q2,0,-1
      DO jQ2=Q2-iQ2,0,-1
         kQ2=Q2-iQ2-jQ2
         DO Q1 = 0,l1
            DO iQ1=Q1,0,-1
               DO jQ1=Q1-iQ1,0,-1
                  kQ1=Q1-iQ1-jQ1
                  IF(iQ1+jQ1+kQ1+iQ2+jQ2+kQ2 .EQ. l1+l2) then
                  ijk=ijk+1
                     DO tQ=0,iQ1+iQ2
                        DO uQ=0,jQ1+jQ2
                           DO vQ=0,kQ1+kQ2
                              signijk = signQ**(tQ+uQ+vQ)
                              ituvQ=TUV%TUVindex(tQ,uQ,vQ)
                              DO iPrimQ=1,z!Q%nPrimitives
                                 Ecoeffs(ijk,ituvQ,iPrimQ) = &
                                      &ETIJ(iPrimQ,tQ,iQ1,iQ2,1)&
                                      &*ETIJ(iPrimQ,uQ,jQ1,jQ2,2)&
                                      &*ETIJ(iPrimQ,vQ,kQ1,kQ2,3)&
                                      &*Q%preExpFac(iPrimQ)*signijk
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDo
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

DEALLOCATE(ETIJ)

END SUBROUTINE BuildEcoeffTensor2

SUBROUTINE SingleHermiteEcoeff(HerInt,P,signP,iAngmom,nAng,nOrb,LUPRI,IPRINT)
implicit none
Real(realk)   :: HerInt(:)
TYPE(Overlap) :: P
Integer       :: iAngmom,nAng,nOrb,LUPRI,IPRINT
Real(realk)   :: signP
!
Integer     :: nPrim1,nPrim2,ang1,ang2
Integer     :: j1,j2,i1,i2
Real(realk)   :: sign
Real(realk), parameter :: D1 = 1.0d0
!
i1 = P%indexAng1(iAngmom)
i2 = P%indexAng2(iAngmom)
j1 = P%orbital1%angmom(i1)
j2 = P%orbital2%angmom(i2)
IF(j1.EQ.0 .AND. j2 .EQ. 0)THEN
   sign = D1
ELSE
   sign = signP**(j1+j2)
ENDIF
nPrim1 = P%orbital1%nPrimitives
nPrim2 = P%orbital2%nPrimitives
CALL SingleHermiteEcoeff1(HerInt,P,sign,j1,j2,nPrim1,nPrim2,nAng,nOrb,LUPRI,IPRINT)
END SUBROUTINE SingleHermiteEcoeff

SUBROUTINE SingleHermiteEcoeff1(HerInt,P,signP,j1,j2,nPrim1,nPrim2,nAng,nOrb,LUPRI,IPRINT)
implicit none
Real(realk)   :: HerInt(nAng*nOrb,nPrim1,nPrim2)
TYPE(Overlap) :: P
Integer       :: nAng,nOrb,LUPRI,IPRINT
Integer       :: nPrim1,nPrim2,j1,j2
Real(realk)   :: signP
!
Real(realk)   :: pref2,pref12,D1=1.0d0,D2=2.0d0
Integer       :: iPrim1,iPrim2,iAngOrb
DO iPrim2=1,nPrim2
  pref2=signP/(D2*P%orbital2%exponents(iPrim2))**j2
  DO iPrim1=1,nPrim1
    pref12=D1/(D2*P%orbital1%exponents(iPrim1))**j1*pref2
    DO iAngOrb=1,nAng*nOrb
      HerInt(iAngOrb,iPrim1,iPrim2) = HerInt(iAngOrb,iPrim1,iPrim2) * pref12
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE SingleHermiteEcoeff1

SUBROUTINE PrintContractEcoeff(EcoeffCont,nAng,nOrb,nPrim,LUPRI,IPRINT)
implicit none
Real(realk) :: EcoeffCont(nAng,nOrb,nPrim)
Integer     :: nAng,nOrb,nPrim,LUPRI,IPRINT
!
Integer :: iAng,iOrb,iPrim
IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(3X,A)') '***************************************************************'
  WRITE(LUPRI,'(3X,A)') '***                       ContractEcoeff'
  WRITE(LUPRI,'(3X,A)') '***************************************************************'
  DO iPrim=1,nPrim
    DO iOrb=1,nOrb
        WRITE(LUPRI,'(5X,A,I3,A,I3)') 'iPrim =',iPrim,', iOrb =',iOrb
        WRITE(LUPRI,'(7X,5F10.4/,(7X,5F10.4))') (EcoeffCont(iAng,iOrb,iPrim),iAng=1,nAng)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintContractEcoeff

SUBROUTINE DistributeHermiteQ(Integral,PQ,iAngmom,ideriv,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: LUPRI,IPRINT,iAngmom,ideriv
!
Integer                 :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ
Integer                 :: nP,nQ,ioffP,ioffQ,sPQ,ioffPQ
Integer                 :: iPrimPQ,iP,iQ,iTUV,maxQ,minQ,nEFG
Integer                 :: ntuvP,ntuvQ,startP,endP,startQ,endQ,endPQ,ntuvPQ
Real(realk)             :: signQ,DM1 = -1.0d0,sign
Real(realk),pointer     :: intermediate1(:,:,:),intermediate2(:,:)

real(4) :: tarr(2),etime,tnew,told,tstart

startP = PQ%P%p%startAngmom
endP   = PQ%P%p%endAngmom
startQ = 0
IF (PQ%Q%p%hermiteSingle) startQ = PQ%Q%p%angmom(iAngmom)
maxQ   = PQ%Q%p%endAngmom
minQ   = PQ%Q%p%startAngmom
endQ   = PQ%Q%p%angmom(iAngmom)
nP     = PQ%P%p%nPrimitives*PQ%P%p%nPasses
nQ     = PQ%Q%p%nPrimitives*PQ%Q%p%nPasses
ioffP  = startP*(startP+1)*(startP+2)/6
ioffQ  = startQ*(startQ+1)*(startQ+2)/6
sPQ    = startP+PQ%Q%p%startAngmom
ioffPQ = sPQ*(sPQ+1)*(sPQ+2)/6
endPQ  = maxQ + PQ%P%p%endAngmom

ntuvP  = PQ%P%p%nTUV
ntuvQ  = (endQ+1)*(endQ+2)*(endQ+3)/6-ioffQ
ntuvPQ = (endPQ+1)*(endPQ+2)*(endPQ+3)/6-ioffPQ
Integral%nOrb  = ntuvP*nP
Integral%nAng  = ntuvQ
Integral%nPrim = nQ
nEFG=Integral%nEFG
IF (((maxQ.EQ.0).OR.((endP.EQ.0).AND.((maxQ.EQ.endQ).AND.(minQ.EQ.startQ)))).AND.(nEFG.EQ.1)) THEN
   call swapRealPointers(Integral%IN,Integral%RTUV)
ELSE
  CALL DistributeHermiteQ1(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
     &                     startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQ*nEFG,ntuvPQ,ntuvP,ntuvQ,ideriv,lupri)
ENDIF
!ENDDO
!tnew = etime(tarr)-tstart
!WRITE(LUPRI,'(A,F14.6,7I5)') 'Distributetimings',tnew,startP,endP,startQ,endQ,minQ,maxQ,nP*nQ

!NULLIFY(Integral%tuvTUV)
!ALLOCATE(Integral%tuvTUV(ntuvQ,ntuvP,nP*nQ))
!CALL DZERO(Integral%tuvTUV,ntuvQ*ntuvP*nP*nQ)
!DO iPrimPQ=1,nP*nQ
!  DO jQ = startQ,endQ
!     signQ = DM1**jQ
!     DO tQ=jQ,0,-1
!        DO uQ=jQ-tQ,0,-1
!           vQ=jQ-tQ-uQ
!           ituvQ=Integral%TUV%TUVindex(tQ,uQ,vQ)-ioffQ
!           DO jP = startP,endP
!              sign = DM1**(jP-jQ)
!              DO tP=jP,0,-1
!                 DO uP=jP-tP,0,-1
!                    vP=jP-tP-uP
!                    ituvP=Integral%TUV%TUVindex(tP,uP,vP)-ioffP
!                    ituvPQ=Integral%TUV%TUVindex(tP+tQ,uP+uQ,vP+vQ)-ioffPQ+(ideriv-1)*ntuvPQ
!                    Integral%tuvTUV(ituvQ,ituvP,iPrimPQ) = &
!       &                signQ*Integral%WTUV(ituvPQ,iPrimPQ)
!!       &                signQ*Integral%WTUV(iPrimPQ,ituvPQ)
!                 ENDDO
!              ENDDO
!           ENDDO
!        ENDDO
!     ENDDO
!  ENDDO
!ENDDO

IF (IPRINT.GT.50) THEN
   write(lupri,*)'  Number of primitives',nP*nQ
   CALL PrintHermitePQ(Integral%IN,Integral%TUV%TUVindex,nP*nQ,&
  & startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvP,ntuvQ,ideriv,lupri)
ENDIF
!
!NULLIFY(Integral%IN)
!ALLOCATE(Integral%IN(ntuvQ,ntuvP*nP,nQ))
!!CALL transposition(Integral%IN,Integral%tuvTUV,ntuvP*ntuvQ,nq*np)
!CALL DCOPY(ntuvP*nP*nQ*ntuvQ,Integral%tuvTUV,1,Integral%IN,1)
!DEALLOCATE(Integral%tuvTUV)
!NULLIFY(Integral%tuvTUV)

END SUBROUTINE DistributeHermiteQ

SUBROUTINE DistributeHermiteQ1(tuvTUV,WTUV,TUVindex,nPrim,startP,endP,startQ,endQ,&
     &                         ioffP,ioffQ,ioffPQ,nTUVEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,lupri)
implicit none
Real(realk)     :: tuvTUV(ntuvQ,ntuvP,nPrim),WTUV(ntuvEFGPQ,nPrim)
Integer         :: nPrim,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ
Integer         :: ntuvEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv
Integer,pointer :: TUVindex(:,:,:)
!
Integer     :: TUVPQindex(nTUVP*nTUVQ)
Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ
Integer     :: iPrimPQ,iOFF,lupri,iituvP,iituvQ,inTUVPQ

iPrimPQ=1
iOFF=(ideriv-1)*ntuvPQ-ioffPQ
ituvP = 0
inTUVPQ = 0
DO jP = startP,endP
   DO tP=jP,0,-1
      DO uP=jP-tP,0,-1
         vP=jP-tP-uP
         ituvP = ituvP+1
         ituvQ = 0
         DO jQ = startQ,endQ
            DO tQ=jQ,0,-1
               DO uQ=jQ-tQ,0,-1
                  vQ=jQ-tQ-uQ
                  ituvQ=ituvQ+1
                  ituvPQ=TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
                  inTUVPQ = inTUVPQ+1
                  TUVPQindex(inTUVPQ)=ituvPQ 
                  tuvTUV(ituvQ,ituvP,iPrimPQ) = WTUV(ituvPQ,iPrimPQ)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
DO iPrimPQ=2,nPrim
   inTUVPQ = 0
   DO iituvP = 1,ituvP
      DO iituvQ = 1,ituvQ
         inTUVPQ = inTUVPQ+1
         ituvPQ = TUVPQindex(inTUVPQ)
         tuvTUV(iituvQ,iituvP,iPrimPQ) = WTUV(ituvPQ,iPrimPQ)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE DistributeHermiteQ1

SUBROUTINE PrintHermitePQ(tuvTUV,TUVindex,nPrim,startP,endP,startQ,endQ,&
    &                   ioffP,ioffQ,ioffPQ,ntuvP,ntuvQ,ideriv,lupri)
implicit none
Integer         :: nPrim,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ
Real(realk)     :: tuvTUV(ntuvQ,ntuvP,nPrim)
Integer         :: ntuvP,ntuvQ,ideriv,lupri
Integer,pointer :: TUVindex(:,:,:)
!
Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ
Integer     :: iPrimPQ

ituvP = 0
WRITE(LUPRI,'(3X,A)') '***************************************************************'
WRITE(LUPRI,'(3X,A)') '***                         HermitePQ'
WRITE(LUPRI,'(3X,A)') '***************************************************************'
DO jP = startP,endP
   DO tP=jP,0,-1
      DO uP=jP-tP,0,-1
         vP=jP-tP-uP
         ituvP = ituvP+1
         ituvQ = 0
         DO jQ = startQ,endQ
            DO tQ=jQ,0,-1
               DO uQ=jQ-tQ,0,-1
                  vQ=jQ-tQ-uQ
                  ituvQ=ituvQ+1
                  WRITE(LUPRI,'(5X,A,I1,A,I1,A,I1,A,I1,A,I1,A,I1,A)') &
                       &            'W(',tP,',',uP,',',vP,'|',tQ,',',uQ,',',vQ,') ='
                  
                  WRITE(LUPRI,'(5X,6F10.4/,(5X,6F10.4))') &
                       &           (tuvTUV(ituvQ,ituvP,iPrimPQ),iPrimPQ=1,nPrim)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE PrintHermitePQ

!!$SUBROUTINE PrintHermitePQ(tuvtuv,np,nq,startP,endP,startQ,endQ,LUPRI,IPRINT)
!!$implicit none
!!$TYPE(Integralitem) :: Integral
!!$Integer            :: np,nq,startP,endP,startQ,endQ
!!$Integer            :: LUPRI,IPRINT
!!$!
!!$Integer            :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ioffP,ioffQ,iPrimP,iPrimQ
!!$!
!!$IF (IPRINT.GT.50) THEN
!!$  WRITE(LUPRI,'(3X,A)') '***************************************************************'
!!$  WRITE(LUPRI,'(3X,A)') '***                         HermitePQ'
!!$  WRITE(LUPRI,'(3X,A)') '***************************************************************'
!!$  ioffP  = startP*(startP+1)*(startP+2)/6
!!$  ioffQ  = startQ*(startQ+1)*(startQ+2)/6
!!$  DO jQ = startQ,endQ
!!$    DO tQ=jQ,0,-1
!!$      DO uQ=jQ-tQ,0,-1
!!$        vQ=jQ-tQ-uQ
!!$        ituvQ=Integral%TUV%TUVindex(tQ,uQ,vQ)-ioffQ
!!$!       WRITE(LUPRI,'(3X,A,3I2)') 'tQ,uQ,vQ',tQ,uQ,vQ
!!$        DO jP = startP,endP
!!$          DO tP=jP,0,-1
!!$            DO uP=jP-tP,0,-1
!!$            vP=jP-tP-uP
!!$            ituvP=Integral%TUV%TUVindex(tP,uP,vP)-ioffP
!!$            WRITE(LUPRI,'(5X,A,I1,A,I1,A,I1,A,I1,A,I1,A,I1,A)') &
!!$      &            'W(',tP,',',uP,',',vP,'|',tQ,',',uQ,',',vQ,') ='
!!$            DO iPrimQ=1,nq
!!$               WRITE(LUPRI,'(17X,I3,5F10.4/,(20X,5F10.4))') iPrimQ,&
!!$      &           (tuvTUV(ituvQ,ituvP,(iPrimP-1)*nq+iPrimQ),iPrimP=1,np)
!!$!      &           (Integral%tuvTUV(ituvQ,ituvP,iPrimP,iPrimQ),iPrimQ=1,nq)
!!$            ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDIF
!!$END SUBROUTINE PrintHermitePQ

!!$SUBROUTINE SPHERICALTRANSMAT(LUPRI,IPRINT,Orb,maxAngmom)
!!$IMPLICIT NONE
!!$TYPE(Orbital) :: Orb
!!$INTEGER  :: LUPRI,IPRINT,JMAX,maxAngmom
!!$Real(realk), parameter :: DM1 = -1.0d0, DO = 0.0d0, D1 = 1.0d0, D2 = 2.0d0
!!$INTEGER  :: M1,L,MADR,MABS,V0, NDER, IOFF
!!$INTEGER  :: EE,FF,GG,BX,BY,BZ,II,JJ,KK,PX,PY,PZ
!!$REAL(realk)  :: FACNRM,FAC1,FAC2, FAC3,FAC4, FACTOR
!!$REAL(realk)  :: SMATRIX(1),PMATRIX(9),DMATRIX(30)  
!!$INTEGER  :: T,U,V,A,B,C,P,Q,R,X,Y,Z,TOT,IADR,AX0,AY0,AZ0,NSIZE
!!$INTEGER  :: M,Ncol,Nrow,INDEX,NDIM,jkh,Ladr,Location(0:maxAngmom)
!!$
!!$
!!$real(realk),pointer :: SPHMAT(:)
!!$integer,pointer     :: SPHMATindex(:)  !1:nangmom
!!$
!!$NSIZE=0
!!$DO L = 0, maxAngmom
!!$   Ladr=0
!!$   DO jkh=1,Orb%nangmom
!!$      IF(L .EQ. Orb%angmom(jkh))Ladr=jkh
!!$   ENDDO
!!$   IF(Ladr .NE. 0)THEN
!!$      NRow = 2*L+1
!!$      NCol = (L+1)*(L+2)/2
!!$      Orb%SPHMATindex(Ladr)=NSIZE+1
!!$      NSIZE=NSIZE+Nrow*Ncol
!!$   ENDIF
!!$ENDDO
!!$
!!$NULLIFY(Orb%SPHMAT)
!!$ALLOCATE(Orb%SPHMAT(NSIZE))
!!$CALL DZERO(Orb%SPHMAT,NSIZE)
!!$DO L = 0, maxAngmom
!!$   Ladr=0
!!$   DO jkh=1,Orb%nangmom
!!$      IF(L .EQ. Orb%angmom(jkh))Ladr=jkh
!!$   ENDDO
!!$   IF(Ladr .NE. 0)THEN
!!$      NRow = 2*L+1
!!$      NCol = (L+1)*(L+2)/2
!!$      DO M1 = 0, 2*L 
!!$         M = M1 - L
!!$         IF (L.EQ.1) THEN
!!$            IF (M .EQ. -1) MADR =  0  
!!$            IF (M .EQ.  0) MADR =  1 
!!$            IF (M .EQ.  1) MADR = -1 
!!$         ELSE
!!$            MADR = M
!!$         END IF
!!$         MABS = ABS(M)
!!$         V0 = 0
!!$         IF (M .LT. 0) V0 = 1 
!!$         FACNRM = D1
!!$         IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(LUPRI,L+MABS)*FACULT(LUPRI,L-MABS))&
!!$           &/(FACULT(LUPRI,L)*(D2**MABS))
!!$         FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
!!$         FACNRM = FACNRM/SQRT(FACUL2(LUPRI,2*L-1))
!!$         DO T = 0, L - MABS, 2
!!$         DO U = 0, T, 2
!!$         DO V = V0, MABS, 2
!!$            !        almost 6.4.48 in the book
!!$            FAC3 = FACNRM*BINOM(LUPRI,L,T/2)*BINOM(LUPRI,L-T/2,MABS+T/2)&
!!$                 &                    *BINOM(LUPRI,T/2,U/2)*BINOM(LUPRI,MABS,V)
!!$            DO A = 0, MIN(0,T+MABS-U-V) 
!!$            DO B = 0, MIN(0,U+V)
!!$            DO C = 0, MIN(0,L-T-MABS)
!!$               !           6.4.47 in the book
!!$               DO P = 0, - A, 2
!!$               DO Q = 0, - B, 2
!!$               DO R = 0, - C, 2
!!$                  FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
!!$                       &   D2**(-A-B-C-P-Q-R-T)*FAC3
!!$                  X = T+MABS-U-V-2*A-P
!!$                  Y = U+V-2*B-Q
!!$                  Z = L-T-MABS-2*C-R
!!$                  TOT = X + Y + Z
!!$                  IADR = 1 + (2*L+1)*(NCRT(X,Y,Z)-1) + L + MADR
!!$                  Orb%SPHMAT(IADR+Orb%SPHMATindex(Ladr)-1) = &
!!$                       &Orb%SPHMAT(IADR+Orb%SPHMATindex(Ladr)-1) + FACTOR 
!!$               ENDDO
!!$               ENDDO
!!$               ENDDO
!!$            ENDDO
!!$            ENDDO
!!$            ENDDO
!!$         ENDDO
!!$         ENDDO
!!$         ENDDO
!!$      ENDDO
!!$   ENDIF
!!$ENDDO
!!$
!!$END SUBROUTINE SPHERICALTRANSMAT


!************************************************************************************************
! ANDY : ECOEFFS - TO BE MOVED LATER  TO A SEPARATE FILE ? TYPEDEFS ????!
!************************************************************************************************
SUBROUTINE GET_ECOEFF(ETIJ,nprimpass,maxang12,maxang1,maxang2,nprim,npass,P,LUPRI,IPRINT)
implicit none
REAL(REALK),PARAMETER    :: D1 = 1.0D0, D2 = 2.0D0, DHALF = 0.5D0
integer                  :: IPRINT, LUPRI, I, J, K, CARTDIR
integer                  :: nprim,maxang12,maxang1,maxang2,l
TYPE(Overlap)            :: P
integer                  :: nprimpass,npass
real(realk) :: PINV(nPrimPass), APINV(nPrimPass), BPINV(nPrimPass), HPINV(nPrimPass), PA(nPrimPass,3), PB(nPrimPass,3)
real(realk) :: TWOA(nPrimPass), TWOB(nPrimPass)
!Allocatable for ECOEFFS - Dimensions are No. Pairs, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
!real(realk),allocatable  :: ETIJ(:,:,:,:,:) 
real(realk)  :: ETIJ(nprimpass,0:maxang12,0:maxang1,0:maxang2,3),Xdist,Ydist,Zdist 

!IF (IPRINT.GT.50) CALL PRN_ECINFO(P,LUPRI)

!CAMT WHERE DO WE DEALLOCATE THIS ?
!ALLOCATE(ETIJ(P%nPrimitives,0:P%orbital1%maxAngmom+P%orbital2%maxAngmom,&
!&             0:P%orbital1%maxAngmom,0:P%orbital2%maxAngmom,3))

!CAMT SETUP REQUIRED EXPONENT RELATED QUANTITIES
IF (P%type(1:4) .NE. 'Empt')THEN
   DO K=1,nPrimPass
      I=P%iprim1(K)
      J=P%iprim2(K)
      PINV(K)     =  D1 / (P%orbital1%exponents(I)& 
           &                 + P%orbital2%exponents(J))    ! 1/p
      APINV(K)    =  P%orbital1%exponents(I) * PINV(K)      ! a/p
      BPINV(K)    =  P%orbital2%exponents(J) * PINV(K)      ! b/p
      HPINV(K)    =  DHALF*PINV(K)                     ! 1/2p
      TWOA(K)     =  D2*P%orbital1%exponents(I)             ! 2a
      TWOB(K)     =  D2*P%orbital2%exponents(J)             ! 2b
   ENDDO
   DO I = 1,nPass
      Xdist=P%distance12(1,I)
      Ydist=P%distance12(2,I)
      Zdist=P%distance12(3,I)
      DO K=1,nPrim
         J = K+(I-1)*nPrim
         !CAMT SETUP REQUIRED DISTANCES
         PA(J,1) = -BPINV(J)*Xdist
         PA(J,2) = -BPINV(J)*Ydist
         PA(J,3) = -BPINV(J)*Zdist
         PB(J,1) =  APINV(J)*Xdist
         PB(J,2) =  APINV(J)*Ydist
         PB(J,3) =  APINV(J)*Zdist
      END DO
   ENDDO
ENDIF

!AMT THEN CALL ROUTINES FOR CARTESIAN OR HERMITE ECOEFFS 
!AMT LOOP OVER CARTESIAN DIRECTONS (1,2,3 -> X,Y,Z)
DO CARTDIR = 1,3
      IF (P%type(1:4) .EQ. 'Herm') THEN
           CALL HERM_ECOEFFS(ETIJ(1,0,0,0,CARTDIR),nPrimPass,maxAng1,&
          &             maxAng2,PA(:,CARTDIR),PB(:,CARTDIR),PINV,APINV,BPINV,HPINV,       &
          &             CARTDIR,TWOA,TWOB,LUPRI)
      ELSE IF (P%type(1:4) .EQ. 'Cart') THEN
           CALL CART_ECOEFFS(ETIJ(1,0,0,0,CARTDIR),nPrimPass,maxAng1,&
          &             maxAng2,PA(:,CARTDIR),PB(:,CARTDIR),PINV,APINV,BPINV,HPINV,       &
          &             CARTDIR,LUPRI)
      ELSE IF (P%type(1:4) .EQ. 'Empt') THEN
           DO i=1,nPrimPass
              ETIJ(i,0,0,0,CARTDIR) = D1
           ENDDO
      ELSE
           WRITE(LUPRI,*)'ERROR : UNRECOGNISED TYPE !!!', P%type
      ENDIF
      IF (IPRINT.GT.20) CALL PRN_ECOEFFS(ETIJ,nPrimPass,maxAng1,&
           &            maxAng2,CARTDIR,LUPRI)
ENDDO
! DEALLOCATE MEMORY
END SUBROUTINE GET_ECOEFF

SUBROUTINE CART_ECOEFFS(ETIJ,nPrimPairs,MAXI,MAXJ,PA,PB,&
                      & PINV,APINV,BPINV,HPINV,X,LUPRI)
!CAMT SUBROUTINE TO CALCULATE E's USING THE USUAL CARTESIAN 
!    RECURANNCE RELATIONS - SEE P354 OF BOOK
!CAMT SEE ERIECF ROUTINE for F77 EQUIVALENT
implicit none
REAL(REALK),PARAMETER  ::  D1 = 1.0D0, D2 = 2.0D0
INTEGER                ::  I,J,K
INTEGER                ::  T, X, nPrimPairs,MAXI,MAXJ,IJ,LUPRI
REAL(REALK)            ::  T1
REAL(REALK)            ::  PA(nPrimPairs), PB(nPrimPairs) 
REAL(REALK)            ::  ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
REAL(REALK)            ::  PINV(nPrimPairs),APINV(nPrimPairs),BPINV(nPrimPairs),HPINV(nPrimPairs)

! VARIABLES
! =========      
! nPrimPairs   NUMBER OF PRIMITIVE PAIRS IN OVERLAP DISTRIBUTION
! MAXI,MAXJ    MAXIMUM I,J VALUES IN ETIJ's (MAX ANGMOM OF G_i and G_j RESPECTIVELY)
! I,J,K,IJ     COUNTERS     
! PA           DISTANCE BETWEEN CENTRE OF CHARGE P AND A (THE CENTRE OF PRIMITIVE GAUSSIAN G_i)           
! PB           DISTANCE BETWEEN CENTRE OF CHARGE P AND B (THE CENTRE OF PRIMITIVE GAUSSIAN G_j)
! HPINV        1/2p (p = a + b I.E. SUM OF EXPONENTS OF PRIM G_i and G_j)           
! PINV         1/p
! APINV, BPINV a/p b/p
! ETIJ         E COEFFICIENTS 
! X            INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
! T            COUNTER INDEX USED TO COUNT UPTO SUM OF ANGULAR MOMENTA OF G_i AND G_j 

!C
!C Run over I (J = 0)                                     
!C ==================
!C
   DO I = 0, MAXI                                          
     IF (I .LE. 2) THEN
      CALL ECOEF_ICASES(nPrimPairs,MAXI,MAXJ,ETIJ,PA,HPINV,I)
     ELSE                                                 
!C                                                               
!C            E(I,0)                                            
!C                                                               
      DO K = 1, nPrimPairs                                   
        ETIJ(K,  0,I,0) = PA(K)*ETIJ(K,  0,I-1,0)+ ETIJ(K,  1,I-1,0)            
        ETIJ(K,I-1,I,0) = HPINV(K)*ETIJ(K,I-2,I-1,0)+ PA(K)*ETIJ(K,I-1,I-1,0)        
        ETIJ(K,  I,I,0) = HPINV(K)*ETIJ(K,I-1,I-1,0)           
      END DO                                             
      DO T = 1, I - 2                                    
!       T1 = DFLOAT(T + 1)                           
        T1 = T + 1
        DO K = 1, nPrimPairs                                
          ETIJ(K,T,I,0) = HPINV(K)*ETIJ(K,T-1,I-1,0)+ PA(K)*ETIJ(K,  T,I-1,0)+ T1*ETIJ(K,T+1,I-1,0)    
        END DO                                         
      END DO                                            
     END IF                                               
!C                                                               
!C   Run over J                                          
!C   ==========                                          
!C                                                               
     DO J = 1, MAXJ                                      
       IJ = I + J                                        
       IF (IJ .LE. 2) THEN
         CALL ECOEF_IJCASES(nPrimPairs,MAXI,MAXJ,ETIJ,PA,PB,HPINV,I,IJ)
       ELSE                                              
!C                                                                
!C             E(I,J)                                            
!C                                                               
         DO K = 1, nPrimPairs                                
           ETIJ(K,   0,I,J) = PB(K)*ETIJ(K,   0,I,J-1)+ ETIJ(K,   1,I,J-1)   
           ETIJ(K,IJ-1,I,J) = HPINV(K)*ETIJ(K,IJ-2,I,J-1)+ PB(K)*ETIJ(K,IJ-1,I,J-1)   
           ETIJ(K,  IJ,I,J) = HPINV(K)*ETIJ(K,IJ-1,I,J-1)  
         END DO                                         
         DO T = 1, IJ - 2                               
!          T1 = DFLOAT(T + 1)                       
           T1 = T + 1
           DO K = 1, nPrimPairs                            
             ETIJ(K,T,I,J)=HPINV(K)*ETIJ(K,T-1,I,J-1)+ PB(K)*ETIJ(K,  T,I,J-1)+ T1*ETIJ(K,T+1,I,J-1)   
           END DO                                      
         END DO                                        
       END IF                                            
     END DO                                               
   END DO                                                  
RETURN
END SUBROUTINE CART_ECOEFFS

SUBROUTINE ECOEF_ICASES(nPrimPairs,MAXI,MAXJ,ETIJ,PA,HPINV,I)
implicit none
REAL(REALK),PARAMETER::  D1 = 1.0D0, D2 = 2.0D0
INTEGER     :: I,nPrimPairs,MAXI,MAXJ,K
REAL(REALK) :: PA(nPrimPairs),HPINV(nPrimPairs)
REAL(REALK) :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
!C                                                                               
!C           E(0,0)                                              
!C
  IF (I .EQ. 0) THEN                               
    DO K = 1, nPrimPairs                               
      ETIJ(K,0,0,0) = D1                          
    END DO
  ELSE IF (I .EQ. 1) THEN                          
!C                                                             
!C           E(1,0)                                            
!C
    DO K = 1, nPrimPairs                               
      ETIJ(K,0,1,0) = PA(K)                       
      ETIJ(K,1,1,0) = HPINV(K)                      
    END DO                                         
  ELSE IF (I .EQ. 2) THEN                          
!C                                                
!C           E(2,0)                                              
!C                                                              
    DO K = 1, nPrimPairs                                
      ETIJ(K,0,2,0) = PA(K)*PA(K) + HPINV(K)     
      ETIJ(K,1,2,0) = D2*PA(K)*HPINV(K)             
      ETIJ(K,2,2,0) = HPINV(K)*HPINV(K)            
    END DO                                         
  ENDIF
END SUBROUTINE ECOEF_ICASES 

SUBROUTINE ECOEF_IJCASES(nPrimPairs,MAXI,MAXJ,ETIJ,PA,PB,HPINV,I,IJ)
implicit none
REAL(REALK),PARAMETER::  D1 = 1.0D0, D2 = 2.0D0
INTEGER     :: I,IJ,nPrimPairs,MAXI,MAXJ,K
REAL(REALK) :: PA(nPrimPairs),PB(nPrimPairs),HPINV(nPrimPairs)
REAL(REALK) :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
  IF (IJ .EQ. 1) THEN                               
!C                                                                 
!C              E(0,1)                                           
!C                                                               
    DO K = 1, nPrimPairs                             
      ETIJ(K,0,0,1) = PB(K)                    
      ETIJ(K,1,0,1) = HPINV(K)                   
    END DO                                      
  ELSE IF (IJ .EQ. 2) THEN                         
!C                                                
!C                 E(0,2)                        
!C                                              
    IF (I .EQ. 0) THEN                          
      DO K = 1, nPrimPairs                          
        ETIJ(K,0,0,2) = PB(K)*PB(K) + HPINV(K)    
        ETIJ(K,1,0,2) = D2*PB(K)*HPINV(K)            
        ETIJ(K,2,0,2) = HPINV(K)*HPINV(K)               
      END DO                                   
    ELSE                                        
!C                                                               
!C                    E(1,1)                                    
!C                                                             
      DO K = 1, nPrimPairs                          
        ETIJ(K,0,1,1) = PA(K)*PB(K) + HPINV(K)    
        ETIJ(K,1,1,1) = (PA(K) + PB(K))*HPINV(K)      
        ETIJ(K,2,1,1) = HPINV(K)*HPINV(K)               
      END DO                                  
    END IF                                     
  ENDIF 
END SUBROUTINE ECOEF_IJCASES

SUBROUTINE PRN_ECOEFFS(ETIJ,nPrimPairs,MAXI,MAXJ,X,LUPRI)
implicit none
INTEGER          ::  MAXI,MAXJ,X,I,J,T,K,nPrimPairs,LUPRI
CHARACTER(len=4) ::  WORD
REAL(REALK)      :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
!C
!C     *********************************
!C     ***** PRINT E COEFFICIENTS  *****
!C     *********************************
!C
         WRITE(LUPRI,*)'Output from LSint ECFs'
         WRITE(LUPRI,*)'----------------------'
         WRITE (LUPRI,'(2X,A,2I5)') 'MAXI,MAXJ   ', MAXI, MAXJ 
            IF (X .EQ. 1) WORD = 'EX00'
            IF (X .EQ. 2) WORD = 'EY00'
            IF (X .EQ. 3) WORD = 'EZ00'
            DO I = 0, MAXI 
            DO J = 0, MAXJ 
            DO T = 0, I + J
               WRITE (LUPRI,'(/,2X,A4,A1,I1,A1,I1,A1,I1,A1,/)')&
     &              WORD, '(', I, ',', J, ', ',  T, ')'
               WRITE (LUPRI,'(1X,6F12.8)') (ETIJ(K,T,I,J),K=1,nPrimPairs)
            END DO
            END DO
            END DO
END SUBROUTINE PRN_ECOEFFS

SUBROUTINE PRN_ECINFO(P,LUPRI)
implicit none
integer                  :: LUPRI, I, J 
TYPE(Overlap)            :: P

WRITE(LUPRI,*)'HERE IS THE RELEVANT INFORMATION FROM THE P'
WRITE(LUPRI,*)'-------------------------------------------'
WRITE(LUPRI,*)'nAngmom',P%nAngmom
WRITE(LUPRI,*)'nPrimitives',P%nPrimitives
WRITE(LUPRI,*)'maxContracted',P%maxContracted
WRITE(LUPRI,*)'Distance Between A and B:'
DO I=1,P%npasses
   WRITE(LUPRI,*)'Pass nr:',I
   WRITE(LUPRI,*)'X',P%distance12(1,I)
   WRITE(LUPRI,*)'Y',P%distance12(2,I)
   WRITE(LUPRI,*)'Z',P%distance12(3,I)
enddo
WRITE(LUPRI,*)'Squared Distance',P%squaredDistance

WRITE(LUPRI,*)'HERE IS THE INFORMATION RELEVANT FOR ORBITAL 1'
WRITE(LUPRI,*)'----------------------------------------------'
WRITE(LUPRI,*)'TYPE',P%orbital1%type
WRITE(LUPRI,*)'SPHERICAL ?',P%orbital1%Spherical
WRITE(LUPRI,*)'Maximum Angular Momentum',P%orbital1%maxAngmom
WRITE(LUPRI,*)'Angular Momentum',P%orbital1%nAngmom
WRITE(LUPRI,*)'Number of Primitives',P%orbital1%nPrimitives
WRITE(LUPRI,*)'Max Contracted',P%orbital1%maxContracted
WRITE(LUPRI,*)'Centre of orbital 1:'
WRITE(LUPRI,*)'X',P%orbital1%center(1)
WRITE(LUPRI,*)'Y',P%orbital1%center(2)
WRITE(LUPRI,*)'Z',P%orbital1%center(3)
WRITE(LUPRI,*)'Primitive Exponents in Orbital1'
WRITE(LUPRI,*)'-------------------------------'
DO I=1,P%orbital1%nPrimitives
   WRITE(LUPRI,*)'Primitive',I,P%orbital1%exponents(I)
ENDDO


WRITE(LUPRI,*)'HERE IS THE INFORMATION RELEVANT FOR ORBITAL 2'
WRITE(LUPRI,*)'----------------------------------------------'
WRITE(LUPRI,*)'TYPE',P%orbital2%type
WRITE(LUPRI,*)'SPHERICAL ?',P%orbital2%Spherical
WRITE(LUPRI,*)'Maximum Angular Momentum',P%orbital2%maxAngmom
WRITE(LUPRI,*)'Angular Momentum',P%orbital2%nAngmom
WRITE(LUPRI,*)'Number of Primitives',P%orbital2%nPrimitives
WRITE(LUPRI,*)'Max Contracted',P%orbital2%maxContracted
WRITE(LUPRI,*)'Centre of orbital 1:'
WRITE(LUPRI,*)'X',P%orbital2%center(1)
WRITE(LUPRI,*)'Y',P%orbital2%center(2)
WRITE(LUPRI,*)'Z',P%orbital2%center(3)
WRITE(LUPRI,*)'Primitive Exponents in Orbital2'
WRITE(LUPRI,*)'-------------------------------'
DO I=1,P%orbital2%nPrimitives
   WRITE(LUPRI,*)'Primitive',I,P%orbital2%exponents(I)
ENDDO
END SUBROUTINE PRN_ECINFO


!*** HERMITE E COEFFICIENTS

SUBROUTINE HERM_ECOEFFS(ETIJ,nPrimPairs,MAXI,MAXJ,PA,PB,&
                      & PINV,APINV,BPINV,HPINV,X,TWOA,TWOB,LUPRI)
!CAMT SUBROUTINE TO CALCULATE E's USING THE USUAL CARTESIAN 
!    RECURANNCE RELATIONS - SEE P354 OF BOOK
!CAMT SEE ERIECF ROUTINE for F77 EQUIVALENT
implicit none
REAL(REALK),PARAMETER  ::  D1 = 1.0D0, D2 = 2.0D0
INTEGER                ::  I,J,K
INTEGER                ::  T, X, nPrimPairs,MAXI,MAXJ,IJ,LUPRI
REAL(REALK)            ::  T1, TIM, TJM
REAL(REALK)            ::  PA(nPrimPairs), PB(nPrimPairs)
REAL(REALK)            ::  ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
REAL(REALK)            ::  PINV(nPrimPairs),APINV(nPrimPairs),BPINV(nPrimPairs),HPINV(nPrimPairs)
REAL(REALK)            ::  TWOA(nPrimPairs),TWOB(nPrimPairs)

! VARIABLES
! =========      
! nPrimPairs   NUMBER OF PRIMITIVE PAIRS IN OVERLAP DISTRIBUTION
! MAXI,MAXJ    MAXIMUM I,J VALUES IN ETIJ's (MAX ANGMOM OF G_i and G_j RESPECTIVELY)
! I,J,K,IJ     COUNTERS     
! PA           DISTANCE BETWEEN CENTRE OF CHARGE P AND A (THE CENTRE OF PRIMITIVE GAUSSIAN G_i)           
! PB           DISTANCE BETWEEN CENTRE OF CHARGE P AND B (THE CENTRE OF PRIMITIVE GAUSSIAN G_j)
! HPINV        1/2p (p = a + b I.E. SUM OF EXPONENTS OF PRIM G_i and G_j)           
! PINV         1/p
! APINV, BPINV a/p b/p
! TWOA, TWOB   2a, 2b
! ETIJ         E COEFFICIENTS 
! X            INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
! T            COUNTER INDEX USED TO COUNT UPTO SUM OF ANGULAR MOMENTA OF G_i AND G_j 

!C
!C Run over I (J = 0)                                     
!C ==================
!C
   DO I = 0, MAXI
!    TIM = DFLOAT(I-1)
     TIM = I-1
     IF (I .LE. 2) THEN
      CALL ECOEF_ICASES_HERM(nPrimPairs,MAXI,MAXJ,ETIJ,PA,HPINV,TWOA,I)
     ELSE
!C                                                               
!C            E(I,0)                                            
!C                                                               
      DO K = 1, nPrimPairs
        ETIJ(K,  0,I,0) = PA(K)*ETIJ(K,  0,I-1,0)+ ETIJ(K,  1,I-1,0)&
                        & - TIM*ETIJ(K, 0,I-2,0)/TWOA(K)
        ETIJ(K,I-1,I,0) = HPINV(K)*ETIJ(K,I-2,I-1,0)+ PA(K)*ETIJ(K,I-1,I-1,0)
        ETIJ(K,  I,I,0) = HPINV(K)*ETIJ(K,I-1,I-1,0)
      END DO                                   
      DO T = 1, I - 2                                    
!       T1 = DFLOAT(T + 1)
        T1 = T + 1
        DO K = 1, nPrimPairs                                
          ETIJ(K,T,I,0) = HPINV(K)*ETIJ(K,T-1,I-1,0)+ PA(K)*ETIJ(K,  T,I-1,0)&
                        & + T1*ETIJ(K,T+1,I-1,0) - TIM*ETIJ(K, T,I-2,0)/TWOA(K)
        END DO
      END DO         
     END IF             
!C                                                               
!C   Run over J                                          
!C   ==========                                          
!C                                                               
     DO J = 1, MAXJ
       IJ = I + J
!      TJM = DFLOAT(J-1)
       TJM = J-1

       IF (IJ .LE. 2) THEN
         CALL ECOEF_IJCASES_HERM(nPrimPairs,MAXI,MAXJ,ETIJ,PA,PB,HPINV,TWOB,I,IJ)
       ELSE
!C                                                                
!C             E(I,J)                                            
!C                                                               
          IF(J.LT.2)THEN
             DO K = 1, nPrimPairs
                ETIJ(K,   0,I,J) = PB(K)*ETIJ(K,   0,I,J-1)+ ETIJ(K,   1,I,J-1)
                ETIJ(K,IJ-1,I,J) = HPINV(K)*ETIJ(K,IJ-2,I,J-1)+ PB(K)*ETIJ(K,IJ-1,I,J-1)
                ETIJ(K,  IJ,I,J) = HPINV(K)*ETIJ(K,IJ-1,I,J-1)
             END DO
             DO T = 1, IJ - 2
                T1 = T + 1
                DO K = 1, nPrimPairs
                   ETIJ(K,T,I,J)= HPINV(K)*ETIJ(K,T-1,I,J-1)+ PB(K)*ETIJ(K,  T,I,J-1)&
                        & + T1*ETIJ(K,T+1,I,J-1)! - TJM*ETIJ(K,  T,I,J-2)/TWOB(K)
                END DO
             END DO
          ELSE
             DO K = 1, nPrimPairs
                ETIJ(K,   0,I,J) = PB(K)*ETIJ(K,   0,I,J-1)+ ETIJ(K,   1,I,J-1)&
                     & -TJM*ETIJ(K,  0,I,J-2)/TWOB(K)
                ETIJ(K,IJ-1,I,J) = HPINV(K)*ETIJ(K,IJ-2,I,J-1)+ PB(K)*ETIJ(K,IJ-1,I,J-1)
                ETIJ(K,  IJ,I,J) = HPINV(K)*ETIJ(K,IJ-1,I,J-1)
             END DO
             DO T = 1, IJ - 2
                T1 = T + 1
                DO K = 1, nPrimPairs
                   ETIJ(K,T,I,J)= HPINV(K)*ETIJ(K,T-1,I,J-1)+ PB(K)*ETIJ(K,  T,I,J-1)&
                        & + T1*ETIJ(K,T+1,I,J-1) - TJM*ETIJ(K,  T,I,J-2)/TWOB(K)
                END DO
             END DO
          ENDIF
       END IF
     END DO
   END DO
RETURN
END SUBROUTINE HERM_ECOEFFS

SUBROUTINE ECOEF_ICASES_HERM(nPrimPairs,MAXI,MAXJ,ETIJ,PA,HPINV,TWOA,I)
implicit none
REAL(REALK),PARAMETER::  D1 = 1.0D0, D2 = 2.0D0
INTEGER     :: I,nPrimPairs,MAXI,MAXJ,K 
REAL(REALK) :: PA(nPrimPairs),HPINV(nPrimPairs)
REAL(REALK) :: TWOA(nPrimPairs) 
REAL(REALK) :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
!C                                                                               
!C           E(0,0)                                              
!C
  IF (I .EQ. 0) THEN
    DO K = 1, nPrimPairs
      ETIJ(K,0,0,0) = D1
    END DO
  ELSE IF (I .EQ. 1) THEN
!C                                                             
!C           E(1,0)                                            
!C
    DO K = 1, nPrimPairs
      ETIJ(K,0,1,0) = PA(K)
      ETIJ(K,1,1,0) = HPINV(K)
    END DO
  ELSE IF (I .EQ. 2) THEN
!C                                                
!C           E(2,0)                                              
!C                                                              
    DO K = 1, nPrimPairs
      ETIJ(K,0,2,0) = PA(K)*PA(K) + HPINV(K) - D1/TWOA(K)
      ETIJ(K,1,2,0) = D2*PA(K)*HPINV(K)
      ETIJ(K,2,2,0) = HPINV(K)*HPINV(K)
    END DO
  ENDIF
END SUBROUTINE ECOEF_ICASES_HERM


SUBROUTINE ECOEF_IJCASES_HERM(nPrimPairs,MAXI,MAXJ,ETIJ,PA,PB,HPINV,TWOB,I,IJ)
implicit none
REAL(REALK),PARAMETER::  D1 = 1.0D0, D2 = 2.0D0
INTEGER     :: I,IJ,nPrimPairs,MAXI,MAXJ,K
REAL(REALK) :: PA(nPrimPairs),PB(nPrimPairs),HPINV(nPrimPairs)
REAL(REALK) :: TWOB(nPrimPairs)
REAL(REALK) :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
  IF (IJ .EQ. 1) THEN
!C                                                                 
!C              E(0,1)                                           
!C                                                               
    DO K = 1, nPrimPairs
      ETIJ(K,0,0,1) = PB(K)    
      ETIJ(K,1,0,1) = HPINV(K) 
    END DO
  ELSE IF (IJ .EQ. 2) THEN     
!C                                                
!C                 E(0,2)                        
!C                                              
    IF (I .EQ. 0) THEN
      DO K = 1, nPrimPairs     
        ETIJ(K,0,0,2) = PB(K)*PB(K) + HPINV(K) - D1/TWOB(K)
        ETIJ(K,1,0,2) = D2*PB(K)*HPINV(K)
        ETIJ(K,2,0,2) = HPINV(K)*HPINV(K)
      END DO
    ELSE
!C                                                               
!C                    E(1,1)                                    
!C                                                             
      DO K = 1, nPrimPairs     
        ETIJ(K,0,1,1) = PA(K)*PB(K) + HPINV(K)
        ETIJ(K,1,1,1) = (PA(K) + PB(K))*HPINV(K)
        ETIJ(K,2,1,1) = HPINV(K)*HPINV(K)
      END DO
    END IF
  ENDIF
END SUBROUTINE ECOEF_IJCASES_HERM

!*******************************************************************************************
! ANDY END OF ECOEFFS ?
!*******************************************************************************************

SUBROUTINE BUILD_RJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral)
IMPLICIT NONE
INTEGER         :: nPrim,LUPRI,IPRINT,Jmax
!REAL(REALK)     :: RJ000(nPrim,0:Jmax)
REAL(REALK)     :: RJ000(0:Jmax,nPrim)
TYPE(Integrand) :: PQ
TYPE(integralitem) :: integral

CALL BUILD_RJ000_OP(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
     &PQ%reducedExponents,PQ%squaredDistance,&
     &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
   
END SUBROUTINE BUILD_RJ000

SUBROUTINE BUILD_NUCLEAR_RJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral)
IMPLICIT NONE
INTEGER         :: nPrim,LUPRI,IPRINT,Jmax
!REAL(REALK)     :: RJ000(nPrim,0:Jmax)
REAL(REALK)     :: RJ000(0:Jmax,nPrim)
TYPE(Integrand) :: PQ
TYPE(integralitem) :: integral

CALL BUILD_RJ000_OP(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
     &PQ%Exponents,PQ%squaredDistance,&
     &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)

END SUBROUTINE BUILD_NUCLEAR_RJ000

SUBROUTINE BUILD_RJ000_OP(RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
IMPLICIT NONE
INTEGER         :: nPrim,Jmax,Lupri,Iprint,nTABFJW1,nTABFJW2
REAL(REALK)     :: RJ000(0:Jmax,nPrim)
REAL(REALK)     :: alpha(nPrim),R2(nprim),Prefactor(nprim)
REAL(REALK)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)
INTEGER         :: nprim1,nprim2,nprim3,INDADS(nPrim,3),I
REAL(REALK)     :: D2JP36,WVALU(Nprim),WVAL,WVALS(Nprim,3)
REAL(REALK),PARAMETER :: HALF =0.5d0,D1=1.d0,D2 = 2.D0, D4 = 4.D0, D10=10.d0
Real(realk),parameter :: D12 = 12.D0, TENTH = 0.1D0
REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6.D0, COEF4 = D1/24.D0
REAL(REALK), PARAMETER :: COEF5 = - D1/120.D0, COEF6 = D1/720.D0
Integer :: IPNT,J
Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA
REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092D0
REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686D0
REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909D0
REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346D0
Real(realk), parameter :: PI=3.14159265358979323846D0
REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730D00
REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
Real(realk) :: W2,W3,W4,W5,W6,R
LOGICAL :: maxJgt0
REAL(REALK), PARAMETER :: SMALL = 0.000001D0

Nprim1 = 0
Nprim2 = 0
Nprim3 = 0
D2JP36 = 2*JMAX + 36

DO I = 1, Nprim
   WVAL = alpha(I)*R2(I)
!  0 < WVAL < 0.000001
   IF (WVAL .LT. SMALL) THEN         
      RJ000(0,I) = D1
      DO J=1,JMAX
         RJ000(J,I)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
      ENDDO
!  0 < WVAL < 12 
   ELSE IF (WVAL .LT. D12) THEN
     IPNT = NINT(D10*WVAL)
     WDIFF = WVAL - TENTH*IPNT
     W2    = WDIFF*WDIFF
     W3    = W2*WDIFF
     W4    = W3*WDIFF
     W5    = W4*WDIFF
     W6    = W5*WDIFF
     R = TABFJW(JMAX,IPNT)
     R = R -TABFJW(JMAX+1,IPNT)*WDIFF
     R = R + COEF2*TABFJW(JMAX+2,IPNT)*W2
     R = R + COEF3*TABFJW(JMAX+3,IPNT)*W3
     R = R + COEF4*TABFJW(JMAX+4,IPNT)*W4
     R = R + COEF5*TABFJW(JMAX+5,IPNT)*W5
     R = R + COEF6*TABFJW(JMAX+6,IPNT)*W6
       
     RJ000(JMAX,I) = R
     IF (JMAX.GT.0) THEN
       REXPW = HALF*EXP(-WVAL)
       DO J=JMAX-1,0,-1
         RJ000(J,I) = (WVAL*RJ000(J+1,I) + REXPW)*D2/(2*J + 1)
       ENDDO
     ENDIF
!  12 < WVAL <= (2J+36) 
   ELSE IF (WVAL.LE.D2JP36) THEN
     REXPW = HALF*EXP(-WVAL)
     RWVAL = D1/WVAL
     GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
     RJ000(0,I) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
     DO J=1,JMAX
       RJ000(J,I) = RWVAL*((J - HALF)*RJ000(J-1,I)-REXPW)
     ENDDO
!  (2J+36) < WVAL 
   ELSE
     RWVAL = PID4/WVAL
     RJ000(0,I) = SQRT(RWVAL)
     RWVAL = RWVAL*PID4I
     DO J = 1, JMAX
       RJ000(J,I) = RWVAL*(J - HALF)*RJ000(J-1,I)
     ENDDO
   END IF
ENDDO

! Scaling
DO I=1,nPrim
  PREF = Prefactor(I)
  RJ000(0,I) = PREF*RJ000(0,I)
  IF (jmax.GT.0) THEN
    D2MALPHA = -2*alpha(I)
    DO j=1,jmax
      PREF = PREF*D2MALPHA
      RJ000(J,I) = PREF*RJ000(J,I)
    ENDDO
  ENDIF
ENDDO

END SUBROUTINE BUILD_RJ000_OP

!!$SUBROUTINE BUILD_RJ000_OP_OLD(RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW)
!!$IMPLICIT NONE
!!$INTEGER         :: nPrim,Jmax,Lupri,Iprint
!!$REAL(REALK)     :: RJ000(nPrim,0:Jmax)
!!$REAL(REALK)     :: alpha(nPrim),R2(nprim),Prefactor(nprim)
!!$REAL(REALK),POINTER :: TABFJW(:,:)
!!$INTEGER         :: nprim1,nprim2,nprim3,INDADS(nPrim,3),I
!!$REAL(REALK)     :: D2JP36,WVALU(Nprim),WVAL,WVALS(Nprim,3)
!!$REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0
!!$
!!$Nprim1 = 0
!!$Nprim2 = 0
!!$Nprim3 = 0
!!$!D2JP36 = DFLOAT(2*JMAX + 36)
!!$D2JP36 = 2*JMAX + 36
!!$
!!$!DO I = 1, Nprim
!!$!   WVALU(I)= alpha(I)*R2(I)
!!$!ENDDO
!!$
!!$DO I = 1, Nprim
!!$   WVAL = alpha(I)*R2(I)
!!$   IF (WVAL .LT. D12) THEN          ! 0 < WVAL < 12 
!!$      Nprim1 = Nprim1 + 1
!!$      INDADS(Nprim1,1) = I
!!$      WVALS(Nprim1,1)  = WVAL
!!$   ELSE IF (WVAL .LE. D2JP36) THEN  ! 12 < WVAL < (2J+36) 
!!$      Nprim2 = Nprim2 + 1
!!$      INDADS(Nprim2,2) = I
!!$      WVALS(Nprim2,2)  = WVAL
!!$   ELSE                             ! (2J+36) < WVAL 
!!$      Nprim3 = Nprim3 + 1
!!$      INDADS(Nprim3,3) = I
!!$      WVALS(Nprim3,3)  = WVAL
!!$   END IF
!!$ENDDO
!!$IF (Nprim1 .GT. 0) THEN
!!$!  WVAL < 12
!!$   CALL GETGAM1(RJ000,nPrim1,nPrim,INDADS(:,1),WVALS(:,1),Jmax,TABFJW)
!!$END IF
!!$IF (Nprim2 .GT. 0) THEN
!!$!  Near asymptotic region   12 < WVAL < (2J+36) 
!!$   CALL GETGAM2(RJ000,nPrim2,nPrim,INDADS(:,2),WVALS(:,2),Jmax)
!!$END IF
!!$IF (Nprim3 .GT. 0) THEN
!!$!  Asymptotic region      (2J+36) < WVAL 
!!$   CALL GETGAM3(RJ000,nPrim3,nPrim,INDADS(:,3),WVALS(:,3),Jmax)
!!$ENDIF
!!$
!!$CALL SCALE_RJ000(RJ000,prefactor,alpha,nPrim,JMAX)
!!$
!!$END SUBROUTINE BUILD_RJ000_OP_OLD

!!$SUBROUTINE GETGAM1(RJ000,nPrim1,nPrim,INDADS,WVALS,Jmax,TABFJW)
!!$IMPLICIT NONE
!!$INTEGER         :: nPrim,Jmax,nPrim1,INDADS(nPrim1)
!!$REAL(REALK)     :: RJ000(nPrim,0:Jmax),TABFJW(:),WVALS(Nprim1)
!!$REAL(REALK)     :: WVAL,WDIF,INTTEMP(nPrim1,0:Jmax),REXPW(Nprim1),FCT
!!$INTEGER         :: ISTRT0,I,IPNT,ISTART,J,M,MP1
!!$REAL(REALK), PARAMETER :: HALF = 0.5D0
!!$REAL(REALK), PARAMETER :: D1 = 1.D0, D10 = 10.D0
!!$REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0
!!$REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6.D0, COEF4 = D1/24.D0
!!$REAL(REALK), PARAMETER :: COEF5 = - D1/120.D0, COEF6 = D1/720.D0
!!$
!!$ISTRT0 = 1 + 121*JMAX
!!$DO I = 1, Nprim1
!!$   WVAL = WVALS(I)
!!$   IPNT = NINT(D10*WVAL)
!!$!  WDIF = WVAL - TENTH*DFLOAT(IPNT)
!!$   WDIF = WVAL - TENTH*IPNT
!!$   ISTART = ISTRT0 + IPNT
!!$   INTTEMP(I,JMAX) = (((((COEF6*TABFJW(ISTART + 726)*WDIF&
!!$        &                         + COEF5*TABFJW(ISTART + 605))*WDIF&
!!$        &                         + COEF4*TABFJW(ISTART + 484))*WDIF&
!!$        &                         + COEF3*TABFJW(ISTART + 363))*WDIF&
!!$        &                         + COEF2*TABFJW(ISTART + 242))*WDIF&
!!$        &                         - TABFJW(ISTART + 121))*WDIF + TABFJW(ISTART)
!!$ENDDO
!!$IF (JMAX .GT. 0) THEN
!!$!   DO I = 1, Nprim1
!!$!      REXPW(I) = HALF*EXP(-WVALS(I))
!!$!   ENDDO
!!$
!!$   M = MOD(NPrim1,4)
!!$   IF(M .EQ. 0) THEN
!!$      DO I = 1,NPrim1,4
!!$         REXPW(I)   = HALF*EXP(-WVALS(I))
!!$         REXPW(I+1) = HALF*EXP(-WVALS(I+1))
!!$         REXPW(I+2) = HALF*EXP(-WVALS(I+2))
!!$         REXPW(I+3) = HALF*EXP(-WVALS(I+3))
!!$      ENDDO
!!$   ELSE
!!$      DO I = 1,M
!!$         REXPW(I) = HALF*EXP(-WVALS(I))
!!$      ENDDO
!!$      IF(NPrim1 .GT. 4)THEN
!!$         MP1=M+1
!!$         DO I = MP1,NPrim1,4
!!$            REXPW(I)   = HALF*EXP(-WVALS(I))
!!$            REXPW(I+1) = HALF*EXP(-WVALS(I+1))
!!$            REXPW(I+2) = HALF*EXP(-WVALS(I+2))
!!$            REXPW(I+3) = HALF*EXP(-WVALS(I+3))
!!$         ENDDO
!!$      ENDIF
!!$   ENDIF
!!$
!!$
!!$   DO J = JMAX - 1, 0, -1
!!$!     FCT = D2/DFLOAT(2*J + 1)
!!$      FCT = D2/(2*J + 1)
!!$      DO I = 1, Nprim1
!!$         INTTEMP(I,J) = FCT*(WVALS(I)*INTTEMP(I,J+1) + REXPW(I))
!!$      ENDDO
!!$   ENDDO
!!$   DO J = 1, JMAX
!!$!      DO I = 1, Nprim1
!!$!         RJ000(INDADS(I),J) = INTTEMP(I,J)
!!$!      ENDDO
!!$
!!$      M = MOD(NPrim1,4)
!!$      IF(M .EQ. 0) THEN
!!$         DO I = 1,NPrim1,4
!!$            RJ000(INDADS(I),J) = INTTEMP(I,J)
!!$            RJ000(INDADS(I+1),J) = INTTEMP(I+1,J)
!!$            RJ000(INDADS(I+2),J) = INTTEMP(I+2,J)
!!$            RJ000(INDADS(I+3),J) = INTTEMP(I+3,J)
!!$         ENDDO
!!$      ELSE
!!$         DO I = 1,M
!!$            RJ000(INDADS(I),J) = INTTEMP(I,J)
!!$         ENDDO
!!$         IF(NPrim1 .GT. 4)THEN
!!$            MP1=M+1
!!$            DO I = MP1,NPrim1,4
!!$               RJ000(INDADS(I),J) = INTTEMP(I,J)
!!$               RJ000(INDADS(I+1),J) = INTTEMP(I+1,J)
!!$               RJ000(INDADS(I+2),J) = INTTEMP(I+2,J)
!!$               RJ000(INDADS(I+3),J) = INTTEMP(I+3,J)
!!$            ENDDO
!!$         ENDIF
!!$      ENDIF
!!$      
!!$   ENDDO
!!$END IF
!!$!DO I = 1, Nprim1
!!$!   RJ000(INDADS(I),0) =  INTTEMP(I,0)
!!$!ENDDO
!!$M = MOD(NPrim1,4)
!!$IF(M .EQ. 0) THEN
!!$   DO I = 1,NPrim1,4
!!$      RJ000(INDADS(I),0) = INTTEMP(I,0)
!!$      RJ000(INDADS(I+1),0) = INTTEMP(I+1,0)
!!$      RJ000(INDADS(I+2),0) = INTTEMP(I+2,0)
!!$      RJ000(INDADS(I+3),0) = INTTEMP(I+3,0)
!!$   ENDDO
!!$ELSE
!!$   DO I = 1,M
!!$      RJ000(INDADS(I),0) = INTTEMP(I,0)
!!$   ENDDO
!!$   IF(NPrim1 .GT. 4)THEN
!!$      MP1=M+1
!!$      DO I = MP1,NPrim1,4
!!$         RJ000(INDADS(I),0) = INTTEMP(I,0)
!!$         RJ000(INDADS(I+1),0) = INTTEMP(I+1,0)
!!$         RJ000(INDADS(I+2),0) = INTTEMP(I+2,0)
!!$         RJ000(INDADS(I+3),0) = INTTEMP(I+3,0)
!!$      ENDDO
!!$   ENDIF
!!$ENDIF
!!$
!!$END SUBROUTINE GETGAM1
!!$
!!$SUBROUTINE GETGAM2(RJ000,nPrim2,nPrim,INDADS,WVALS,Jmax)
!!$IMPLICIT NONE
!!$!     Near asymptotic region   12 < WVAL < (2J+36) 
!!$INTEGER         :: nPrim,Jmax,nPrim2,INDADS(nPrim2)
!!$REAL(REALK)     :: RJ000(nPrim,0:Jmax),WVALS(Nprim2)
!!$INTEGER         :: I,J
!!$REAL(REALK)     :: REXPW(Nprim2),GVAL,RWVAL
!!$REAL(REALK)     :: INTTEMP(nPrim2,0:Jmax),FCT
!!$REAL(REALK), PARAMETER :: D1 = 1.D0, D10 = 10.D0
!!$REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0
!!$REAL(REALK), PARAMETER :: HALF = 0.5D0
!!$REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092D0 
!!$REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686D0   
!!$REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909D0
!!$REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346D0
!!$REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730D00
!!$REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
!!$
!!$  DO I = 1, Nprim2
!!$     REXPW(I)   = HALF*EXP(-WVALS(I))
!!$  ENDDO
!!$  DO I = 1, Nprim2
!!$     WVALS(I) = D1/WVALS(I)
!!$  ENDDO
!!$  DO I = 1, Nprim2
!!$     RWVAL = WVALS(I)
!!$     GVAL = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
!!$     INTTEMP(I,0) = SQRPIH*SQRT(RWVAL) - REXPW(I)*GVAL*RWVAL
!!$  ENDDO
!!$  DO I = 1, Nprim2
!!$     RJ000(INDADS(I),0) = INTTEMP(I,0)
!!$  ENDDO
!!$  DO J = 1, JMAX
!!$!    FCT = DFLOAT(J) - HALF
!!$     FCT = J - HALF
!!$     DO I = 1, Nprim2
!!$        INTTEMP(I,J) = (FCT*INTTEMP(I,J-1) - REXPW(I))*WVALS(I)
!!$     ENDDO
!!$     DO I = 1, Nprim2
!!$        RJ000(INDADS(I),J) = INTTEMP(I,J)
!!$     ENDDO
!!$  ENDDO
!!$END SUBROUTINE GETGAM2
!!$
!!$SUBROUTINE GETGAM3(RJ000,nPrim3,nPrim,INDADS,WVALS,Jmax)
!!$!     Asymptotic region      (2J+36) < WVAL 
!!$IMPLICIT NONE
!!$INTEGER         :: nPrim,Jmax,nPrim3,INDADS(nPrim3)
!!$REAL(REALK)     :: RJ000(nPrim,0:Jmax),WVALS(Nprim3)
!!$INTEGER         :: I,J
!!$REAL(REALK)     :: REXPW(Nprim3),GVAL,RWVAL,FACTOR
!!$REAL(REALK)     :: INTTEMP(nPrim3,0:Jmax),FCT
!!$Real(realk), parameter :: PI=3.14159265358979323846D0
!!$REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0
!!$REAL(REALK), PARAMETER :: HALF = 0.5D0
!!$REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
!!$
!!$  DO I = 1, Nprim3
!!$     WVALS(I) = PID4/WVALS(I)
!!$  ENDDO
!!$  DO I = 1, Nprim3
!!$     INTTEMP(I,0) = SQRT(WVALS(I))
!!$  ENDDO
!!$  DO I = 1, Nprim3
!!$     RJ000(INDADS(I),0) = INTTEMP(I,0)
!!$  ENDDO
!!$  DO J = 1, JMAX
!!$!    FACTOR = PID4I*(DFLOAT(J) - HALF)
!!$     FACTOR = PID4I*(J - HALF)
!!$     DO I = 1, Nprim3
!!$        INTTEMP(I,J) = FACTOR*INTTEMP(I,J-1)*WVALS(I)
!!$     ENDDO
!!$     DO I = 1, Nprim3
!!$        RJ000(INDADS(I),J) = INTTEMP(I,J)
!!$     ENDDO
!!$  ENDDO
!!$
!!$END SUBROUTINE GETGAM3
!!$
!!$SUBROUTINE SCALE_RJ000(RJ000,Prefactor,alpha,nprim,Jmax)
!!$IMPLICIT NONE
!!$INTEGER         :: I,nprim,Jmax,J
!!$REAL(REALK)     :: RJ000(nPrim,0:Jmax),prefac(Nprim),prefac2(Nprim)
!!$REAL(REALK)     :: Prefactor(nprim)
!!$REAL(REALK)     :: alpha(nprim)
!!$REAL(REALK), PARAMETER :: D2 = 2.D0
!!$
!!$IF (JMAX .GT. 0) THEN
!!$   CALL DCOPY(nprim,prefactor,1,prefac2,1)
!!$   DO I = 1, Nprim
!!$      prefac(I)= -D2*alpha(I)
!!$   ENDDO
!!$
!!$   DO J = 0, JMAX-1
!!$      DO I = 1, Nprim
!!$         RJ000(I,J) = prefac2(I)*RJ000(I,J)
!!$      ENDDO
!!$      DO I = 1, Nprim
!!$         prefac2(I) = prefac(I)*prefac2(I)
!!$      ENDDO
!!$   ENDDO
!!$   DO I = 1, Nprim
!!$      RJ000(I,JMAX) = prefac2(I)*RJ000(I,JMAX)
!!$   ENDDO
!!$ELSE
!!$   DO I = 1, Nprim
!!$      RJ000(I,0) = prefactor(I)*RJ000(I,0)
!!$   ENDDO
!!$ENDIF
!!$
!!$END SUBROUTINE SCALE_RJ000

SUBROUTINE GAMMA_TABULATION(LUPRI,JMX,sharedTUV)
use TYPEDEF
!
!     ***** Tabulation of incomplete gamma function *****
!
!     For J = JMX a power series expansion is used, see for
!     example Eq.(39) given by V. Saunders in "Computational
!     Techniques in Quantum Chemistry and Molecular Physics",
!     Reidel 1975.  For J < JMX the values are calculated
!     using downward recursion in J.
!
IMPLICIT NONE
INTEGER           :: JMX,MAXJ0,IADR,IPOINT,IORDER,JADR,JMAX,I,J,LUPRI
INTEGER,PARAMETER :: MXQN=13
INTEGER,PARAMETER :: MAXJ = 4*(MXQN - 1) + 2
TYPE(TUVitem):: sharedTUV
REAL(REALK)       :: DENOM,D2MAX1,R2MAX1,TERM,SUM,REXPW,WVAL,D2WAL
REAL(REALK), PARAMETER :: HALF = 0.5D0,  TEN6 = 1.0D6
REAL(REALK), PARAMETER :: D1 = 1.D0, D10 = 10.D0
REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0
REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6.D0, COEF4 = D1/24.D0
REAL(REALK), PARAMETER :: COEF5 = - D1/120.D0, COEF6 = D1/720.D0

REAL(REALK), PARAMETER :: PI    = 3.14159265358979323846D00
REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730D00
REAL(REALK), PARAMETER :: R2PI52 = 5.91496717279561287782D00
REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI

REAL(REALK), PARAMETER :: GFAC30 =  .4999489092D0 
REAL(REALK), PARAMETER :: GFAC31 = -.2473631686D0
REAL(REALK), PARAMETER :: GFAC32 =  .321180909D0
REAL(REALK), PARAMETER :: GFAC33 = -.3811559346D0
REAL(REALK), PARAMETER :: GFAC20 = .4998436875D0
REAL(REALK), PARAMETER :: GFAC21 = -.24249438D0
REAL(REALK), PARAMETER :: GFAC22 =  .24642845D0
REAL(REALK), PARAMETER :: GFAC10 =  .499093162D0
REAL(REALK), PARAMETER :: GFAC11 = -.2152832D0
REAL(REALK), PARAMETER :: GFAC00 =  .490D0

IF (JMX .GT. MAXJ) THEN
   WRITE (LUPRI,'(//A,I5,A,I3)')&
   &      ' GAMTAB ERROR: JMX =',JMX,', which is greater than',MAXJ
   CALL QUIT('GAMTAB ERROR: JMX greater than limit.')
END IF
JMAX = JMX + 6
MAXJ0 = JMAX
!
!     WVAL = 0.0
!
IADR = 1
DENOM = D1
DO J = 0,JMAX
   SHAREDTUV%TABFJW(J,0) = D1/DENOM
   IADR = IADR + 121
   DENOM = DENOM + D2
ENDDO
!
!     WVAL = 0.1, 0.2, 0.3,... 12.0
!
IADR = IADR - 121
!D2MAX1 = DFLOAT(2*JMAX + 1)
D2MAX1 = 2*JMAX + 1
R2MAX1 = D1/D2MAX1
DO IPOINT = 1,120
!  WVAL = TENTH*DFLOAT(IPOINT)
   WVAL = TENTH*IPOINT
   D2WAL = WVAL + WVAL
   IADR = IADR + 1
   TERM = R2MAX1
   SUM = TERM
   DENOM = D2MAX1
   DO IORDER = 2,200
      DENOM = DENOM + D2
      TERM = TERM*D2WAL/DENOM
      SUM = SUM + TERM
      IF (TERM .LE. 1.0D-15) EXIT
   ENDDO
   REXPW = EXP(-WVAL)
   SHAREDTUV%TABFJW(JMAX,IPOINT) = REXPW*SUM
   DENOM = D2MAX1
   JADR = IADR
   DO J = 1,JMAX
      DENOM = DENOM - D2
      SHAREDTUV%TABFJW(JMAX-J,IPOINT) = (SHAREDTUV%TABFJW(JMAX-J+1,IPOINT)*D2WAL + REXPW)/DENOM
      JADR = JADR - 121
   ENDDO
ENDDO
END SUBROUTINE GAMMA_TABULATION

SUBROUTINE PrintTensor(Tensor,label,dim1,dim2,dim3,Lupri,Label1,Label2,Label3,option)
integer          :: dim1,dim2,dim3,lupri,option
real(realk)      :: Tensor(dim1,dim2,dim3)
character(len=20):: Label
character(len=6) :: Label1,Label2,label3

  WRITE(LUPRI,'(5X,A)') '***************************************************'
  WRITE(LUPRI,'(5X,A,A20)') '***               ',LABEL
  WRITE(LUPRI,'(5X,A)') '***************************************************'
IF(option .EQ. 1)THEN
   WRITE(LUPRI,'(1X,A6,1X,A6,12X,A6 )') Label2,label3,label1
   DO j=1,dim2
      DO k=1,dim3
         WRITE(LUPRI,'(1X,I4,2X,I4,5F9.4/,(11X,5F9.4))') j,k,(Tensor(i,j,k),i=1,dim1)
      ENDDO
   ENDDO
ELSEIF(option .EQ. 2)THEN
   WRITE(LUPRI,'(1X,A6,1X,A6,12X,A6 )') Label1,label3,label2
   DO i=1,dim1
      DO k=1,dim3
         WRITE(LUPRI,'(1X,I4,2X,I4,5F9.4/,(11X,5F9.4))') i,k,(Tensor(i,j,k),j=1,dim2)
      ENDDO
   ENDDO
ELSEIF(option .EQ. 3)THEN
   WRITE(LUPRI,'(1X,A6,1X,A6,12X,A6 )') Label1,label2,label3
   DO i=1,dim1
      DO j=1,dim2
         WRITE(LUPRI,'(1X,I4,2X,I4,5F9.4/,(11X,5F9.4))') i,j,(Tensor(i,j,k),k=1,dim3)
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE PrintTensor

SUBROUTINE Determine_RHS_loop(INPUT,OD_RHS,ILHS,Start_RHS,End_RHS)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(ODITEM)         :: OD_RHS
Integer              :: ILHS,Start_RHS,End_RHS
! Screening-integral loop
  IF (INPUT%CS_int) THEN
    Start_RHS = ILHS
    End_RHS   = ILHS
! Triangular loop
  ELSE IF (INPUT%sameODs) THEN
    Start_RHS = 1
    End_RHS   = ILHS
! Regular loop
  ELSE
    Start_RHS = 1
    End_RHS   = OD_RHS%nbatches
  ENDIF
END SUBROUTINE Determine_RHS_loop

SUBROUTINE FREE_OVERLAP(P)
Implicit none
TYPE(Overlap) :: P

CALL FREE_ORBITAL(P%orbital1)
CALL FREE_ORBITAL(P%orbital2)
 
DEAllocate(P%center)
DEAllocate(P%exponents)
DEAllocate(P%reducedExponents)
DEAllocate(P%preExpFac)
DEAllocate(P%iprim1)
DEAllocate(P%iprim2)
Nullify(P%center)
Nullify(P%exponents)
Nullify(P%reducedExponents)
Nullify(P%preExpFac)
Nullify(P%iprim1)
Nullify(P%iprim2)
  
DEAllocate(P%distance12)
nullify(P%distance12)

DEAllocate(P%angmom)
Nullify(P%angmom)
DEAllocate(P%indexAng1)
Nullify(P%indexAng1)
DEAllocate(P%indexAng2)
Nullify(P%indexAng2)

DEAllocate(P%nOrbitals)
Nullify(P%nOrbitals)
DEAllocate(P%nOrbComp)
Nullify(P%nOrbComp)
DEAllocate(P%nContracted)
Nullify(P%nContracted)

IF (P%ETUVisSet) THEN
  DEALLOCATE(P%ETUV)
  DEALLOCATE(P%ETUVindex)
ENDIF

END SUBROUTINE FREE_OVERLAP

SUBROUTINE FREE_ORBITAL(Orb)
IMPLICIT NONE
TYPE(Orbital) :: Orb
 
DEALLOCATE(Orb%exponents)
NULLIFY(Orb%exponents)
 
DEALLOCATE(Orb%angmom)
NULLIFY(Orb%angmom)
 
DEALLOCATE(Orb%nContracted)
NULLIFY(Orb%nContracted)
 
DEALLOCATE(Orb%startOrbital)
NULLIFY(Orb%startOrbital)

DEALLOCATE(Orb%startprimOrbital)
NULLIFY(Orb%startprimOrbital)
 
DEALLOCATE(Orb%nOrbComp)
NULLIFY(Orb%nOrbComp)
 
DEALLOCATE(Orb%nPrimOrbComp)
NULLIFY(Orb%nPrimOrbComp)
 
DEALLOCATE(Orb%nOrbitals)
NULLIFY(Orb%nOrbitals)

DEALLOCATE(Orb%CC)
NULLIFY(Orb%CC)

IF(ASSOCIATED(Orb%SPH_MAT)) DEALLOCATE(Orb%SPH_MAT)
NULLIFY(Orb%SPH_MAT)

END SUBROUTINE FREE_ORBITAL

SUBROUTINE CLEAR_integral(INTEGRAL)
TYPE(Integralitem)      :: integral
 
!DEALLOCATE(INTEGRAL%Wtuv)
NULLIFY(INTEGRAL%Wtuv)
 
END SUBROUTINE

SUBROUTINE AddToCoulombMat(JMat,Jdim1,Jdim2,Jengine,sameLHSaos,sameRHSaos,sameODs,Jfac,integral,&
     &                     iA,iB,iC,iD,startA,startB,startC,startD,Dmat,Ddim1,Ddim2)
implicit none
Real(realk)          :: integral
Integer              :: Jdim1,Jdim2,Ddim1,Ddim2
Real(realk)          :: Jmat(Jdim1,Jdim2),Dmat(Ddim1,Ddim2)
Real(realk)          :: Jfac
Logical              :: Jengine,sameLHSaos,sameRHSaos,sameODs
Integer              :: iA,iB,iC,iD,startA,startB,startC,startD
!
Real(realk) :: dtemp_lhs,dtemp_rhs

IF (Jengine) THEN
  JMat(iA,iB) = JMat(iA,iB) + integral
  IF (sameLHSaos) THEN
    IF (startA.NE.startB) JMAT(iB,iA) = JMAT(iB,iA) + integral
  ENDIF
ELSE
  dtemp_rhs = Dmat(iD,iC)
  if (startC.NE.startD.AND.sameRHSaos) dtemp_rhs=dtemp_rhs+Dmat(iC,iD)
  dtemp_rhs = dtemp_rhs*Jfac
  
  IF (sameODs) THEN
    dtemp_lhs = Dmat(iB,iA)
    if (startA.NE.startB.AND.sameLHSaos) dtemp_lhs=dtemp_lhs+Dmat(iA,iB)
    dtemp_lhs = dtemp_lhs*Jfac
  ENDIF
  
  
  IF (sameLHSaos.AND.sameRHSaos.AND.sameODs)THEN
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
    IF (startA.NE.startB) JMAT(iB,iA) = JMAT(iB,iA) + dtemp_rhs*integral
    IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
      JMAT(iC,iD) = JMAT(iC,iD) + dtemp_lhs*integral
      IF (startC.NE.startD) JMAT(iD,iC) = JMAT(iD,iC) + dtemp_lhs*integral
    ENDIF
  ELSE IF (sameLHSaos)THEN
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
    IF (startA.NE.startB) JMAT(iB,iA) = JMAT(iB,iA) + dtemp_rhs*integral
  ELSE IF (sameRHSaos)THEN
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
  ELSE IF (sameODs .AND. (sameLHSaos.NEQV.sameRHSaos))THEN
    CALL QUIT('Progamming error, sameODs but sameLHSaos.NE.sameRHSaos!')
  ELSE IF (sameODs)THEN
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
    IF ((startA.NE.startC).OR.(startB.NE.startD)) JMAT(iC,iD) = JMAT(iC,iD) + dtemp_lhs*integral
  ELSE                                
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
  ENDIF
ENDIF
END SUBROUTINE AddToCoulombMat

SUBROUTINE AddToExchangeMat(KMat,Kdim1,Kdim2,sameLHSaos,sameRHSaos,sameODs,Kfac,integral,&
     &                      iA,iB,iC,iD,startA,startB,startC,startD,Dmat,Ddim1,Ddim2)
implicit none
Real(realk) :: integral
Integer     :: Kdim1,Kdim2,Ddim1,Ddim2
Real(realk) :: KMat(Kdim1,Kdim2),Dmat(Ddim1,Ddim2)
Real(realk) :: Kfac
Logical     :: sameLHSaos,sameRHSaos,sameODs
Integer              :: iA,iB,iC,iD,startA,startB,startC,startD
!
Real(realk) :: dtemp_lhs,dtemp_rhs,kint

kint = integral*Kfac

IF (sameLHSaos.AND.sameRHSaos.AND.sameODs)THEN
  KMat(iA,iC) = KMat(iA,iC) - kint*Dmat(iB,iD)
  IF (startA.NE.startB) KMat(iB,iC) = KMat(iB,iC) - kint*Dmat(iA,iD)
  IF (startC.NE.startD) KMat(iA,iD) = KMat(iA,iD) - kint*Dmat(iB,iC)
  IF ((startA.NE.startB).AND.(startC.NE.startD)) KMat(iB,iD) = KMat(iB,iD) - kint*Dmat(iA,iC)
  IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
    KMat(iC,iA) = KMat(iC,iA) - kint*DMat(iD,iB)
    IF (startA.NE.startB) KMat(iC,iB) = KMat(iC,iB) - kint*DMat(iD,iA)
    IF (startC.NE.startD) KMat(iD,iA) = KMat(iD,iA) - kint*DMat(iC,iB)
    IF ((startA.NE.startB).AND.(startC.NE.startD)) KMat(iD,iB) = KMat(iD,iB) - kint*DMat(iC,iA)
  ENDIF
!ELSE IF (sameRHSaos) THEN
!  KMat(iA,iC) = KMat(iA,iC) - kint*Dmat(iB,iD)
!  IF (startC.NE.startD) KMat(iA,iD) = KMat(iA,iD) - kint*Dmat(iB,iC)
!ELSE IF (sameLHSaos) THEN
!  KMat(iA,iC) = KMat(iA,iC) - kint*Dmat(iB,iD)
!  IF (startA.NE.startB) KMat(iB,iC) = KMat(iB,iC) - kint*DMat(iA,iD)
ELSE
!  KMat(iA,iC) = KMat(iA,iC) - kint*Dmat(iB,iD)
  CALL QUIT('Progamming error AddToExchangeMat, only works if all four AOs are identical!')
ENDIF
END SUBROUTINE AddToExchangeMat

SUBROUTINE initTUVitem(sharedTUV,Integral,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
implicit none
TYPE(TUVitem),target :: sharedTUV
TYPE(Integralitem)   :: Integral
TYPE(IntegralInput)  :: Input
TYPE(ODITEM)         :: OD_LHS,OD_RHS
Integer              :: LUPRI,IPRINT
!
LOGICAL              :: TABULATE_BOYS
Integer,parameter :: MAXAOJ=12,MAXDER=2
Integer,parameter :: MAXJ=4*MAXAOJ+MAXDER
!Integer,parameter :: MAXAOJ=12,MAXDER=2
!Integer,parameter :: MAXJ=12
Integer,parameter :: NODES=121,ORDER=7
Integer,parameter :: NTABFJ000=NODES*(MAXJ+ORDER)
Integer,parameter :: NTUVMAX = (MAXJ+1)*(MAXJ+2)*(MAXJ+3)/6
!
TABULATE_BOYS = .FALSE.
IF(INPUT%Operator(1:7) .EQ. 'Coulomb') TABULATE_BOYS = .TRUE.
IF(INPUT%Operator(1:6).EQ.'Nucrep') TABULATE_BOYS = .TRUE.
IF(INPUT%Operator(1:4) .EQ. 'Erfc') TABULATE_BOYS = .TRUE.
IF(INPUT%Operator(1:4) .EQ. 'Erf ') TABULATE_BOYS = .TRUE.

IF(TABULATE_BOYS)THEN
   !Pretabulate Gammafunction/Boysfunction 
   NULLIFY(sharedTUV%TABFJW)  !MXQN=13     MAXJ = 4*(MXQN - 1) + 2 
   ALLOCATE(SharedTUV%TABFJW(0:MAXJ+6,0:120))     !121*(MAXJ + 7) = 6897
   SharedTUV%nTABFJW1=MAXJ+6
   SharedTUV%nTABFJW2=120
   CALL GAMMA_TABULATION(lupri,MAXJ,sharedTUV)    !Jmax is set to 10 is that ok?
ENDIF
NULLIFY(SharedTUV%TUVindex)
NULLIFY(SHAREDTUV%Tindex)
NULLIFY(SHAREDTUV%Uindex)
NULLIFY(SHAREDTUV%Vindex)
ALLOCATE(SharedTUV%TUVindex(0:MAXJ,0:MAXJ,0:MAXJ))
ALLOCATE(SHAREDTUV%Tindex(NTUVMAX))
ALLOCATE(SHAREDTUV%Uindex(NTUVMAX))
ALLOCATE(SHAREDTUV%Vindex(NTUVMAX))

CALL integralTUVindex(SharedTUV,MAXJ,LUPRI,IPRINT)

CALL Build_PRECALCULATED_SPHMAT(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)

integral%TUV => sharedTUV

END SUBROUTINE initTUVitem

SUBROUTINE Build_PRECALCULATED_SPHMAT(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)
IMPLICIT NONE
TYPE(ODITEM)     :: OD_LHS,OD_RHS
TYPE(TUVitem)    :: SharedTUV
INTEGER          :: LUPRI,IPRINT
!
!Real(realk), parameter :: DM1 = -1.0d0, DO = 0.0d0, D1 = 1.0d0, D2 = 2.0d0
!INTEGER  :: M1,L,MADR,MABS,V0, NDER, IOFF,I,MAXANGMOM
!INTEGER  :: EE,FF,GG,BX,BY,BZ,II,JJ,KK,PX,PY,PZ
!REAL(realk)  :: FACNRM,FAC1,FAC2, FAC3,FAC4, FACTOR
!INTEGER  :: T,U,V,A,B,C,P,Q,R,X,Y,Z,TOT,IADR,AX0,AY0,AZ0,NSIZE
!INTEGER  :: M,Ncol,Nrow,INDEX,NDIM,nangmom
INTEGER  :: I,MAXANGMOM

MAXANGMOM = 0
DO I=1,OD_LHS%nbatches
   MAXANGMOM = MAX(MAXANGMOM,OD_LHS%BATCH(I)%AO(1)%p%maxangmom)
   MAXANGMOM = MAX(MAXANGMOM,OD_LHS%BATCH(I)%AO(2)%p%maxangmom)
ENDDO
DO I=1,OD_RHS%nbatches
   MAXANGMOM = MAX(MAXANGMOM,OD_RHS%BATCH(I)%AO(1)%p%maxangmom)
   MAXANGMOM = MAX(MAXANGMOM,OD_RHS%BATCH(I)%AO(2)%p%maxangmom)
ENDDO

SharedTUV%nSPHMAT = MAXANGMOM + 1
CALL set_Precalculated_SPHMAT(SharedTUV%SPH_MAT,SharedTUV%nSPHMAT,LUPRI,IPRINT)

!nANGMOM=MAXANGMOM+1 
!
!NULLIFY(SharedTUV%SPH_MAT)
!ALLOCATE(SharedTUV%SPH_MAT(nANGMOM))
!SharedTUV%nSPHMAT=nANGMOM
!DO I=1,nANGMOM
!   IF(I.EQ. 1)THEN ! S
!      NULLIFY(SharedTUV%SPH_MAT(I)%elms)
!      ALLOCATE(SharedTUV%SPH_MAT(I)%elms(1))
!      SharedTUV%SPH_MAT(I)%elms(1)=1.d0
!   ELSEIF(I.EQ. 2)THEN ! P
!      NULLIFY(SharedTUV%SPH_MAT(I)%elms)
!      ALLOCATE(SharedTUV%SPH_MAT(I)%elms(9))
!      CALL DZERO(SharedTUV%SPH_MAT(I)%elms,9)
!      SharedTUV%SPH_MAT(I)%elms(1)=1.d0
!      SharedTUV%SPH_MAT(I)%elms(5)=1.d0
!      SharedTUV%SPH_MAT(I)%elms(9)=1.d0
!   ELSEIF(I .GT. 2)THEN
!      L = I-1 !angmom
!      NRow = 2*L+1
!      NCol = (L+1)*(L+2)/2
!      NSIZE= Nrow*Ncol !(2*L+1)*(L+1)*(L+2)/2
!      NULLIFY(SharedTUV%SPH_MAT(I)%elms)
!      ALLOCATE(SharedTUV%SPH_MAT(I)%elms(NSIZE))
!      CALL DZERO(SharedTUV%SPH_MAT(I)%elms,NSIZE)
!      DO M1 = 0, 2*L 
!         M = M1 - L
!         IF (L.EQ.1) THEN
!            IF (M .EQ. -1) MADR =  0  
!            IF (M .EQ.  0) MADR =  1 
!            IF (M .EQ.  1) MADR = -1 
!         ELSE
!            MADR = M
!         END IF
!         MABS = ABS(M)
!         V0 = 0
!         IF (M .LT. 0) V0 = 1 
!         FACNRM = D1
!         IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(LUPRI,L+MABS)*&
!              &FACULT(LUPRI,L-MABS))/(FACULT(LUPRI,L)*(D2**MABS))
!         FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
!         FACNRM = FACNRM/SQRT(FACUL2(LUPRI,2*L-1))
!         DO T = 0, L - MABS, 2
!         DO U = 0, T, 2
!         DO V = V0, MABS, 2
!            !        almost 6.4.48 in the book
!            FAC3 = FACNRM*BINOM(LUPRI,L,T/2)*BINOM(LUPRI,L-T/2,MABS+T/2)&
!                 &                    *BINOM(LUPRI,T/2,U/2)*BINOM(LUPRI,MABS,V)
!            DO A = 0, MIN(0,T+MABS-U-V) 
!            DO B = 0, MIN(0,U+V)
!            DO C = 0, MIN(0,L-T-MABS)
!               !           6.4.47 in the book
!               DO P = 0, - A, 2
!               DO Q = 0, - B, 2
!               DO R = 0, - C, 2
!                  FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
!                       &   D2**(-A-B-C-P-Q-R-T)*FAC3
!                  X = T+MABS-U-V-2*A-P
!                  Y = U+V-2*B-Q
!                  Z = L-T-MABS-2*C-R
!                  TOT = X + Y + Z
!                  IADR = 1 + (2*L+1)*(NCRT(X,Y,Z)-1) + L + MADR
!                  SharedTUV%SPH_MAT(I)%elms(IADR) = &
!                       & SharedTUV%SPH_MAT(I)%elms(IADR) + FACTOR 
!               ENDDO
!               ENDDO
!               ENDDO
!            ENDDO
!            ENDDO
!            ENDDO
!         ENDDO
!         ENDDO
!         ENDDO
!      ENDDO
!   ENDIF
!ENDDO

END SUBROUTINE BUILD_PRECALCULATED_SPHMAT

SUBROUTINE integralTUVindex(sharedTUV,MAXJ,LUPRI,IPRINT)
implicit none
TYPE(TUVitem)  :: sharedTUV
Integer        :: LUPRI,IPRINT,MAXJ
!
Integer :: TUV,T,U,V,J
TUV=0
DO J = 0, MAXJ
   DO T = J,0,-1       
      DO U = J-T,0,-1
         TUV=TUV+1
         V=J-T-U
         sharedTUV%TUVindex(T,U,V)=TUV
         sharedTUV%Tindex(TUV)=T
         sharedTUV%Uindex(TUV)=U
         sharedTUV%Vindex(TUV)=V
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE integralTUVindex

SUBROUTINE freeTUVitem(sharedTUV,Input)
implicit none
TYPE(TUVItem)  :: sharedTUV
TYPE(IntegralInput) :: Input
INTEGER             :: I
LOGICAL              :: TABULATE_BOYS
!
TABULATE_BOYS = .FALSE.
IF(INPUT%Operator(1:7) .EQ. 'Coulomb') TABULATE_BOYS = .TRUE.
IF(INPUT%Operator(1:6).EQ.'Nucrep') TABULATE_BOYS = .TRUE.
IF(INPUT%Operator(1:4) .EQ. 'Erfc') TABULATE_BOYS = .TRUE.
IF(INPUT%Operator(1:4) .EQ. 'Erf ') TABULATE_BOYS = .TRUE.

IF(TABULATE_BOYS) DEALLOCATE(sharedTUV%TABFJW)
DEALLOCATE(SharedTUV%TUVindex)
DEALLOCATE(SHAREDTUV%Tindex)
DEALLOCATE(SHAREDTUV%Uindex)
DEALLOCATE(SHAREDTUV%Vindex)

!DO I=1,SharedTUV%nSPHMAT
!   DEALLOCATE(SharedTUV%SPH_MAT(I)%elms)
!   NULLIFY(SharedTUV%SPH_MAT(I)%elms)
!ENDDO
!
!DEALLOCATE(SharedTUV%SPH_MAT)
!NULLIFY(SharedTUV%SPH_MAT)

CALL free_Precalculated_SPHMAT(SharedTUV%SPH_MAT,SharedTUV%nSPHMAT)

END SUBROUTINE freeTUVitem

! Subroutine that returns A, the transpose of B
SUBROUTINE transposition(A,B,n1,n2)
implicit none
integer     :: n1,n2
Real(realk) :: A(n1,n2)
Real(realk) :: B(n2,n1)
A = transpose(B)
END SUBROUTINE transposition

SUBROUTINE SET_FTUVbatches(F,NFTUVbatches,OD,Q,Input,SharedTUV,Integral,Alloc,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)   :: INPUT
TYPE(Overlap),pointer :: F(:)
TYPE(Overlap),pointer :: Q(:)
Integer               :: NFTUVbatches,LUPRI,IPRINT
TYPE(ODITEM)          :: OD 
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)    :: Alloc
TYPE(TUVitem)         :: SharedTUV
!
Integer :: I,J,nUnique,ftuvindex,np,ndmat
Integer :: ODtoFTUVindex(OD%nbatches),nPrimOD(OD%nbatches),Identifier(OD%nbatches)
Integer :: UniqeIdentfiers(OD%nbatches),FTUVprim(OD%nbatches)
Integer :: FTUVminAng(OD%nbatches),FTUVmaxAng(OD%nbatches),FTUVntuv(OD%nbatches)
Integer :: offPrim(OD%nbatches)
Logical :: unique
Integer :: nPrimQ,nTUV
!Simen  maxPrim is an empirical parameter fixed to reduce CPU-timings. 
!Simen  Two opposing effects occur: 
!Simen     1. timings in ContractFTUV reduce with incresing number of primitives
!Simen        per FTUV-batch.
!Simen     2. timings in GET_WTUV increase with increasing number of primitives.
!Simen  Parameter was fixed using ifort/SUN X2200 Opteron/2.2 Ghz (titan.uio.no)
Integer, parameter :: maxPrim = 64

ndmat = 1

nUnique = 0
DO I = 1,OD%nbatches
 CALL SET_Overlap(Q(I),nPrimQ,Input,SharedTUV,Integral,Alloc,OD%BATCH(I),2,LUPRI,IPRINT,'RHS')
 IF (nPrimQ.GT.0) THEN
  np = nPrimQ
  nPrimOD(I) = np
!Simen  Add center information to this identifier - maybe some integer value of 
!Simen  a weighted OD-center divided by BOXSIZE. Note the larger the boxsize, the 
!Simen  less effective CS-screening becomes, whereas the FTUV-contraction step becomes
!Simen  faster. I other words there is a tradeoff between two opposing effects here.
!Simen  Note also that Maxgab should be inherited from OD-batches to FTUV-batches
  IF (Q(I)%maxAngmom.GT.98) CALL QUIT('Non-unique angmom-identifier in SET_FTUVbatches')
  Identifier(I) = (Q(I)%maxAngmom+1)*100 + (Q(I)%minAngmom+1)
  unique = .true.
  DO J=nUnique,1,-1
    IF (Identifier(I).EQ.UniqeIdentfiers(J)) THEN
      IF ((FTUVprim(J)+np).GT.maxPrim) EXIT
      ftuvIndex = J
      unique = .false.
      exit
     ENDIF
  ENDDO
  IF (unique) THEN
    nUnique = nUnique + 1
    ftuvIndex = nUnique
    UniqeIdentfiers(nUnique) = Identifier(I)
    FTUVprim(ftuvIndex) = Q(I)%nPrimitives
    FTUVminAng(ftuvIndex) = Q(I)%startAngmom
    FTUVmaxAng(ftuvIndex) = Q(I)%endAngmom
    FTUVntuv(ftuvIndex)   = (Q(I)%endAngmom+1)*(Q(I)%endAngmom+2)*(Q(I)%endAngmom+3)/6 &
     &                     - Q(I)%startAngmom*(Q(I)%startAngmom+1)*(Q(I)%startAngmom+2)/6
  ELSE
    FTUVprim(ftuvIndex) = FTUVprim(ftuvIndex) + Q(I)%nPrimitives
  ENDIF
  ODtoFTUVindex(I) = ftuvIndex
 ELSE
  ODtoFTUVindex(I) = 0
 ENDIF
ENDDO
NFTUVbatches = nUnique

Alloc%maxPrimRHS    = 0
Alloc%maxPrimTUVRHS = 0
DO I=1,NFTUVbatches
  np =  FTUVprim(I)
  nTUV =  FTUVnTUV(I)
  Alloc%maxPrimRHS    = max(np,Alloc%maxPrimRHS)
  Alloc%maxPrimTUVRHS = max(np*nTUV,Alloc%maxPrimTUVRHS)
ENDDO

!Initialize FTUV-batches
ALLOCATE(F(NFTUVbatches))
DO J=1,NFTUVbatches
  CALL InitFTUVbatches(F(J),FTUVprim(J),FTUVminAng(J),FTUVmaxAng(J),FTUVntuv(J),ndmat)
  offPrim(J) = 0
ENDDO

!Add OD-batches to corresponding FTUV-batches
DO I=1,OD%nbatches
 J = ODtoFTUVindex(I)
 IF (J.GT.0) THEN
  np = nPrimOD(I)
  CALL AddODtoFTUV(F(J),np,offPrim(J),Q(I),ndmat,LUPRI,IPRINT)
  CALL FREE_OVERLAP(Q(I))
  DEALLOCATE(Q(I)%FTUV)
  NULLIFY(Q(I)%FTUV)
  offPrim(J) = offPrim(J) + np
 ENDIF
ENDDO

END SUBROUTINE SET_FTUVbatches

SUBROUTINE AddODtoFTUV(F,np,offp,Q,ndmat,LUPRI,IPRINT)
implicit none
TYPE(Overlap) :: F,Q
Integer       :: np,offp,ndmat,LUPRI,IPRINT
!
Integer :: ip,idir,idmat,tuv
IF (F%nTUV.NE.Q%nTUV) THEN
  CALL QUIT('Programming error: nTUV maismatch in AddODtoFTUV')
ENDIF
DO ip=1,np
  F%preExpFac(ip+offp) = Q%preExpFac(ip)
  DO idir=1,3
    F%center(idir,ip+offp) = Q%center(idir,ip)
  ENDDO
  F%reducedExponents(ip+offp) = Q%reducedExponents(ip)
  F%exponents(ip+offp) = Q%exponents(ip)
ENDDO
DO idmat = 1,ndmat
  DO tuv=1,F%nTUV
    DO ip=1,np
      F%FTUV(ip+offp,tuv,idmat) = Q%FTUV(ip,tuv,idmat)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE AddODtoFTUV

SUBROUTINE InitFTUVbatches(F,nPrim,minAng,maxAng,nTUV,ndmat)
implicit none
TYPE(Overlap) :: F
Integer       :: nPrim,minAng,maxAng,nTUV,ndmat
!
NULLIFY(F%angmom)
NULLIFY(F%orbital1%angmom)
NULLIFY(F%orbital2%angmom)
NULLIFY(F%orbital1%nContracted)
NULLIFY(F%orbital2%nContracted)
NULLIFY(F%orbital1%nOrbComp)
NULLIFY(F%orbital2%nOrbComp)
NULLIFY(F%orbital1%startOrbital)
NULLIFY(F%orbital2%startOrbital)
NULLIFY(F%indexAng1)
NULLIFY(F%indexAng2)
NULLIFY(F%nOrbitals)
NULLIFY(F%nOrbComp)
NULLIFY(F%nContracted)
ALLOCATE(F%angmom(1))
ALLOCATE(F%orbital1%angmom(1))
ALLOCATE(F%orbital2%angmom(1))
ALLOCATE(F%orbital1%nContracted(1))
ALLOCATE(F%orbital2%nContracted(1))
ALLOCATE(F%orbital1%nOrbComp(1))
ALLOCATE(F%orbital2%nOrbComp(1))
ALLOCATE(F%orbital1%startOrbital(1,1))
ALLOCATE(F%orbital2%startOrbital(1,1))
ALLOCATE(F%indexAng1(1))
ALLOCATE(F%indexAng2(1))
ALLOCATE(F%nOrbitals(1))
ALLOCATE(F%nOrbComp(1))
ALLOCATE(F%nContracted(1))

F%TYPE                     = 'FTUV'
F%nAngmom                  = 1
F%orbital1%nAngmom         = 1
F%orbital2%nAngmom         = 1
F%nPrimitives              = nPrim
F%nPasses                  = 1
F%minAngmom                = minAng
F%maxAngmom                = maxAng
F%startAngmom              = minAng
F%endAngmom                = maxAng
F%nTUV                     = nTUV
F%totOrbitals              = 1
F%angmom(1)                = maxAng
F%orbital1%angmom(1)       = 0
F%orbital2%angmom(1)       = 0
F%orbital1%nContracted(1)  = 1
F%orbital2%nContracted(1)  = 1
F%orbital1%nOrbComp(1)     = 1
F%orbital2%nOrbComp(1)     = 1
F%orbital1%startOrbital(1,1) = 1
F%orbital2%startOrbital(1,1) = 1
F%indexAng1(1)             = 1
F%indexAng2(1)             = 1
F%nOrbitals(1)             = 1
F%nOrbComp(1)              = 1
F%nContracted(1)           = 1

NULLIFY(F%preExpFac)
NULLIFY(F%center)
NULLIFY(F%reducedExponents)
NULLIFY(F%exponents)
NULLIFY(F%FTUV)
ALLOCATE(F%preExpFac(nPrim))
ALLOCATE(F%center(3,nPrim))
ALLOCATE(F%reducedExponents(nPrim))
ALLOCATE(F%exponents(nPrim))
ALLOCATE(F%FTUV(nPrim,nTUV,ndmat))

NULLIFY(F%distance12)
ALLOCATE(F%distance12(3,1))

END SUBROUTINE InitFTUVbatches

SUBROUTINE FREE_FTUVOverlap(F)
implicit none
TYPE(Overlap) :: F
!
DEALLOCATE(F%angmom)
DEALLOCATE(F%orbital1%angmom)
DEALLOCATE(F%orbital2%angmom)
DEALLOCATE(F%orbital1%nContracted)
DEALLOCATE(F%orbital2%nContracted)
DEALLOCATE(F%orbital1%nOrbComp)
DEALLOCATE(F%orbital2%nOrbComp)
DEALLOCATE(F%orbital1%startOrbital)
DEALLOCATE(F%orbital2%startOrbital)
DEALLOCATE(F%indexAng1)
DEALLOCATE(F%indexAng2)
DEALLOCATE(F%nOrbitals)
DEALLOCATE(F%nOrbComp)
DEALLOCATE(F%nContracted)
DEALLOCATE(F%preExpFac)
DEALLOCATE(F%center)
DEALLOCATE(F%reducedExponents)
DEALLOCATE(F%exponents)
DEALLOCATE(F%FTUV)
DEALLOCATE(F%distance12)
NULLIFY(F%distance12)

END SUBROUTINE FREE_FTUVOverlap

SUBROUTINE MAXDISTANCE(MAXDIST,distance,nPrim)
IMPLICIT NONE
INTEGER      :: nPrim,I
REAL(REALK)  :: distance(nPrim),MAXDIST

MAXDIST = sum ( abs ( distance(1:1+(nPrim-1)*1:1) ) )

END SUBROUTINE MAXDISTANCE

SUBROUTINE GET_MOMTUV(INTEGRAL,PQ,LUPRI,IPRINT,CARTORDER)

! The final MOM(T,U,V) integrals are arranged as follows:
! S(000)
! S(100) S(010) S(001)
! S(200) S(110) S(101) S(020) S(011) S(002)
! S(300) S(210) S(201) S(120) S(111) S(102) S(030) S(021) S(012) S(003)
implicit none
TYPE(Integrand)         :: PQ
TYPE(Integralitem)      :: integral
real(realk),ALLOCATABLE :: INTTEMP(:,:)!TUV*EFG,nPrim
real(realk),ALLOCATABLE :: INTTEMP2(:,:)!TUV,nPrim
real(realk),ALLOCATABLE :: M000(:)
INTEGER                 :: LUPRI,SUM,J,K,T,U,V,TUV,IOFF
INTEGER                 :: nPrim,IPRINT,ntuv,L,I
INTEGER,ALLOCATABLE     :: TUVindex(:,:,:)
real(realk),ALLOCATABLE :: DISTANCE(:,:)
INTEGER                 :: zeroX,zeroY,zeroZ,Jmax,Jstart,e,f,g,efg,nefg
INTEGER                 :: Xorder,Yorder,Zorder,CARTORDER,C
real(realk)             :: X0,Y0,Z0

!CARTORDER=Xorder+Yorder+Zorder
Xorder=CARTORDER
Yorder=CARTORDER
Zorder=CARTORDER

NPrim=PQ%nPrimitives
JMAX=PQ%endAngmom!+CARTORDER
Jstart=PQ%startAngmom !always zero
ALLOCATE(M000(nPrim))
CALL BUILD_M000(M000,PQ,nPrim,LUPRI,IPRINT)

ntuv=(JMAX+1)*(JMAX+2)*(JMAX+3)/6
nEFG=(CARTORDER+1)*(CARTORDER+2)*(CARTORDER+3)/6
INTEGRAL%nTUV=ntuv
INTEGRAL%nEFG=nefg

IF (JMAX+CARTORDER .EQ. 0) THEN
   !Special case JMAX = 0 A SIMPEL OVERLAP MATRIX WITH ONLY S ORBITALS
   !WRITE(LUPRI,*)'SPECIAL CASE OF CARTISIAN MULTIPOLE MOMENTS - OVERLAP'
   NULLIFY(INTEGRAL%Rtuv)
   ALLOCATE(INTEGRAL%Rtuv(Nprim))
   CALL DCOPY(nPrim,M000,1,Integral%Rtuv,1)
ELSE

   IF (IPRINT .GE. 10) THEN
      CALL LSHEADER(LUPRI,'Output from Recurrence M(0,0,0)')
      J = 0; T=0; U=0; V=0
      WRITE (LUPRI,'(2X,3(A,I1),A,2X,5F12.8/,(12X,5F12.8))')&
           & 'M000(',T,',',U,',',V,')', &
           &(M000(I),I=1,NPrim)
      WRITE (LUPRI,*)' '
   END IF
   
   IF(CARTORDER .EQ. 0) THEN
!!$      !SIMPEL OVERLAP
!!$      NULLIFY(INTEGRAL%Wtuv)
!!$      ALLOCATE(INTEGRAL%Wtuv(Nprim,nTUV))
!!$      DO J = Jstart, PQ%endAngmom
!!$         DO T = J,0,-1
!!$            DO U = J-T,0,-1
!!$               V=J-T-U
!!$               DO I=1,nPrim
!!$                  INTEGRAL%Wtuv(I,INTEGRAL%TUVINDEX(T,U,V))=INTTEMP2(I,INTEGRAL%TUVindex(T,U,V))
!!$               ENDDO
!!$            ENDDO
!!$         ENDDO
!!$      ENDDO
   ELSE
 !     WRITE(LUPRI,'(2X,A,3F16.9)')'CALCULATE: DISTANCE FROM ',PQ%ORIGO(1),PQ%ORIGO(2),PQ%ORIGO(3)
      
      ALLOCATE(DISTANCE(nprim,3))

      DO I = 1,Nprim
         DISTANCE(I,1) = PQ%P%p%center(1,I)-PQ%ORIGO(1)
         DISTANCE(I,2) = PQ%P%p%center(2,I)-PQ%ORIGO(2)
         DISTANCE(I,3) = PQ%P%p%center(3,I)-PQ%ORIGO(3)
 !        WRITE(LUPRI,'(2X,A,3F16.9)')'DISTANCE =',DISTANCE(I,1),DISTANCE(I,2),DISTANCE(I,3)
      ENDDO

!      NULLIFY(INTEGRAL%Rtuv)
!      ALLOCATE(INTEGRAL%Rtuv(nTUV*Nprim*nEFG))
      CALL DZERO(INTEGRAL%Rtuv,nTUV*nPrim*nEFG)

 !     WRITE(LUPRI,*)'QQ ALLOCATION IS nprim=',nprim
 !     WRITE(LUPRI,*)'QQ ALLOCATION IS nTUV =',nTUV
 !     WRITE(LUPRI,*)'QQ ALLOCATION IS nEFG =',nEFG

 !     WRITE(LUPRI,*)'CARTORDER',CARTORDER,'CALL Multipole recurrence'
      CALL Multipole_Recurrence(INTEGRAL%Rtuv,M000,Jstart,JMAX,CARTORDER,&
           & distance(:,1),distance(:,2),distance(:,3),&
           & Xorder,Yorder,Zorder,NPrim,nTUV,nEFG,Integral%TUV%TUVindex,&
           & zeroX, zeroY, zeroZ,PQ%P%p%Exponents,lupri)
      DEALLOCATE(DISTANCE)
   ENDIF

   !    Print section
   !    =============
   
   IF (IPRINT .GE. 10) THEN
      CALL LSHEADER(LUPRI,'Output from ERITUV')
      WRITE (LUPRI,'(2X,A13,I10)') 'MAX angmom ', JMAX
      WRITE (LUPRI,'(2X,A13,I10)') 'NPrim  ', NPrim
      WRITE (LUPRI,'(2X,A13,I10)') 'NTUV  ', nTUV
      WRITE (LUPRI,'(2X,A13,I10)') 'NEFG  ', nEFG
      IF (IPRINT .GE. 20) THEN
         CALL LSHEADER(LUPRI,'Hermite integrals M(t,u,v)')
         efg=0
         DO C=0,CARTORDER
            DO e = C,0,-1
               DO f = C-e,0,-1
                  g=C-e-f
                  efg=efg+1
                  WRITE(LUPRI,*)'MULTIPOLE MOMENT LEVEL = (',e,',',f,',',g,')   EFG=',efg
                  DO J = 0, JMAX
                     DO T = J,0,-1
                        DO U = J-T,0,-1
                           V=J-T-U
                           WRITE(LUPRI,*)'MULTIPOLE MOMENT LEVEL = (',T,',',U,',',V,')   TUV=',INTEGRAL%TUV%TUVINDEX(T,U,V)


                           WRITE (LUPRI,'(2X,A7,I1,A1,I1,A1,I1,A1,2X,5F12.8/,(12X,5F12.8))')&
                                & 'CARMOM(',T,',',U,',',V,')', &
                                &(INTEGRAL%Rtuv(INTEGRAL%TUV%TUVINDEX(T,U,V)+(efg-1)*nTUV+(I-1)*nTUV*nEFG),I=1,NPrim)
                           WRITE (LUPRI,*) ' '
                           CALL FLUSH(LUPRI)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      END IF
   END IF
ENDIF      
DEALLOCATE(M000)
Integral%nPrim=nPrim

END SUBROUTINE GET_MOMTUV

SUBROUTINE Multipole_Recurrence(OUTPUT,M000,Jstart,JMAX,CARTORDER,&
        & Xpq,Ypq,Zpq,Xorder,Yorder,Zorder,NPrim,nTUV,nEFG,TUVindex,&
        & zeroX, zeroY, zeroZ,EXPONENTS,lupri)
IMPLICIT NONE
INTEGER                 :: J,K,T,U,V,TUV,JVAL,JMAX,NTUV,nEFG
INTEGER                 :: nPrim
real(realk)             :: OUTPUT(NTUV*nPrim*nEFG),M000(nPrim)
INTEGER,PARAMETER       :: MAXJ = 50
INTEGER                 :: TUVindex(0:MAXJ,0:MAXJ,0:MAXJ)
real(realk)             :: Xpq(nPrim),Ypq(nPrim),Zpq(nPrim),EXPONENTS(nPrim)
real(realk)             :: P(nPrim)
INTEGER                 :: MAXT,MAXU,MAXV,TMIN1,M1T,M2T,M1U,M2U
INTEGER                 :: M1V,M2V,VMIN1,UMIN1,IUMAX,I,e,f,g,efg
INTEGER                 :: zeroX,zeroY,zeroZ,Xorder,Yorder,Zorder
INTEGER                 :: Tm1,Tp1,jstart,CARTORDER,lupri,C
real(realk),allocatable :: MX(:,:,:),MXY(:,:,:,:,:),MXYZ(:,:,:,:,:,:,:)

!CALL LSHEADER(LUPRI,'Multipole recurrence')

DO I = 1, Nprim
   P(I)=1/(2.d0*EXPONENTS(I))   
ENDDO
! CARTORDER = 0 is a special case already taken care of

ALLOCATE(MX(nPrim,0:Xorder,0:Xorder))

CALL DZERO(MX,nPrim*(Xorder+1)*(Xorder+1))
IF(Xorder .GE. 1)THEN
   IF(Xorder .GT. 1)THEN
 !     WRITE(LUPRI,*)'XORDER GREATER THAN 1 SO STANDARD LOOP'
      CALL STANDARD_LOOPX(M000,MX,Xpq,Xorder,nPrim,P,JMAX,lupri)
   ELSEIF(Xorder .EQ. 1)THEN
      CALL ORDER_EQ_ONE_LOOP(M000,MX,Xpq,Xorder,nPrim,JMAX,lupri)
   ENDIF
   !M(nprim,e,T) done
ELSE
   CALL DCOPY(nPrim,M000,1,MX(1,0,0),1)
   !M(nprim,e,T) done
ENDIF


ALLOCATE(MXY(nPrim,0:Xorder,0:Xorder,0:Yorder,0:Yorder))
CALL DZERO(MXY,nPrim*(Xorder+1)*(Xorder+1)*(Yorder+1)*(Yorder+1))
IF(Yorder .GE. 1)THEN
   DO e = 0,Xorder
      DO T = 0,MIN(e,JMAX)
         !      DO U = J-T,0,-1
         IF(Yorder .GT. 1)THEN
            CALL STANDARD_LOOPX(MX(:,e,T),MXY(:,e,T,:,:),Ypq,Yorder,nPrim,P,JMAX-T,lupri)
         ELSEIF(Yorder .EQ. 1)THEN
            CALL ORDER_EQ_ONE_LOOP(MX(:,e,T),MXY(:,e,T,:,:),Ypq,Yorder,nPrim,JMAX-T,lupri)
         ENDIF
      ENDDO
   ENDDO
ELSE
   DO e = 1,Xorder
      DO T = 0,MIN(e,JMAX)
         DO I=1,Nprim
            MXY(I,e,T,0,0)=MX(I,e,T)
         enddo
      ENDDO
   ENDDO
ENDIF

ALLOCATE(MXYZ(nPrim,0:Xorder,0:Xorder,0:Yorder,0:Yorder,0:Zorder,0:Zorder))
CALL DZERO(MXYZ,nPrim*(Xorder+1)*(Xorder+1)*(Yorder+1)*(Yorder+1)*(Zorder+1)*(Zorder+1))
IF(Zorder .GE. 1)THEN
   DO e = 0,Xorder
      DO T = 0,MIN(e,JMAX)
         DO f = 0,MIN(Yorder,CARTORDER-e)
            DO U = 0,MIN(f,JMAX-T)
               !      DO U = J-T,0,-1
               IF(Zorder .GT. 1)THEN
                  CALL STANDARD_LOOPX(MXY(:,e,T,f,U),MXYZ(:,e,T,f,U,:,:),Zpq,Zorder,nPrim,P,JMAX-T-U,lupri)
               ELSEIF(Zorder .EQ. 1)THEN
                  CALL ORDER_EQ_ONE_LOOP(MXY(:,e,T,f,U),MXYZ(:,e,T,f,U,:,:),Zpq,Zorder,nPrim,JMAX-T-U,lupri)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ELSE
   DO e = 0,Xorder
      DO T = 0,MIN(e,JMAX)
         DO f = 0,Yorder
            DO U = 0,MIN(f,JMAX-T)
               DO I=1,Nprim
                  MXYZ(I,e,T,f,U,0,0)=MXY(I,e,T,f,U)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF
efg=0
!Nefg=(CARTORDER+1)*(CARTORDER+2)*(CARTORDER+3)/6
DO C=0,CARTORDER
   DO e = C,0,-1
      DO f = C-e,0,-1
         g=C-e-f
         efg=efg+1
         TUV=0
         DO J = 0, JMAX
            DO T = J,0,-1
               DO U = J-T,0,-1
                  V=J-T-U
                  TUV=TUV+1
                  IF(T .LE. e)THEN
                     IF(U .LE. f)THEN
                        IF(V .LE. g)THEN
                           DO I=1,nPrim
                              OUTPUT(TUV+(efg-1)*nTUV+(I-1)*nTUV*nEFG)=MXYZ(I,e,T,f,U,g,V)
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
CALL FLUSH(LUPRI)
DEALLOCATE(MX,MXY,MXYZ)

END SUBROUTINE Multipole_Recurrence

SUBROUTINE BUILD_M000(M000,PQ,nPrim,LUPRI,IPRINT)
IMPLICIT NONE
TYPE(Integrand)        :: PQ
real(realk)            :: M000(nPrim)
INTEGER                :: nPrim,J,I,LUPRI,IPRINT

CALL BUILD_M000_OP(M000,Nprim,PQ%exponents,PQ%P%p%reducedExponents&
     &,PQ%squaredDistance)

IF(IPRINT .GT. 25) THEN
   DO I=1,NPrim
      WRITE(LUPRI,'(2X,A5,I4,A1,I2,A2,F16.9)')'M000(',I,',',0,')=',M000(I)
   ENDDO
ENDIF

END SUBROUTINE BUILD_M000

SUBROUTINE BUILD_M000_OP(M000,nPrim,P,alpha,R2)
IMPLICIT NONE
real(realk)            :: M000(nPrim)
Real(realk), parameter :: PI=3.14159265358979323846D0
Real(realk), parameter :: OneHalf=1.5d0
INTEGER                :: nPrim,I
REAL(REALK)            :: P(nPrim)!exponents
REAL(REALK)            :: ALPHA(nPrim),R2(nprim)

DO I=1,NPrim
   M000(I)=((PI/P(I))**OneHalf)!*EXP(-Alpha(I)*R2(I))
ENDDO

END SUBROUTINE BUILD_M000_OP

SUBROUTINE STANDARD_LOOPX(INPUT,M,Xpq,Xorder,nPrim,P,LIMIT,lupri)
!written for X loop (therefore Xpq and Xorder) but general for any direction
IMPLICIT NONE
INTEGER                 :: nPrim,LIMIT,C,Xorder,T,lupri,I
real(realk)             :: Xpq(nPrim),P(nPrim)
real(realk)             :: M(nPrim,0:Xorder,0:Xorder),INPUT(nPrim)

CALL DCOPY(nPrim,INPUT,1,M(1,0,0),1)

!M(prim,e,t)
DO C=1,Xorder-1
   IF(C .EQ. 1 )THEN
      DO I = 1, Nprim
         M(I,1,0)=Xpq(I)*INPUT(I)  
      ENDDO
      DO I = 1, Nprim
         M(I,1,1)=INPUT(I)                        
      ENDDO
   ELSEIF(C .EQ. 2)THEN
      DO I = 1, Nprim
         M(I,2,0)=Xpq(I)*M(I,1,0)+P(I)*M(I,1,1)     
      ENDDO
      DO I = 1, Nprim
         M(I,2,1)=M(I,1,0)+Xpq(I)*M(I,1,1)          
      ENDDO
      DO I = 1, Nprim
         M(I,2,2)=2*M(I,1,1)                      
      ENDDO
   ELSE
      DO I = 1, Nprim
         M(I,C,0)= Xpq(I)*M(I,C-1,0)+P(I)*M(I,C-1,1)     
      ENDDO
      DO T = 1,C-2
         DO I = 1, Nprim
            M(I,C,T)=T*M(I,C-1,T-1)+Xpq(I)*M(I,C-1,T)+P(I)*M(I,C-1,T+1)
         ENDDO
      ENDDO
      !        T=C-1
      DO I = 1, Nprim
         M(I,C,C-1)=(C-1)*M(I,C-1,C-2)+Xpq(I)*M(I,C-1,C-1)
      ENDDO
      !        T=C
      DO I = 1, Nprim
         M(I,C,C)=C*M(I,C-1,C-1)
      ENDDO
   ENDIF
ENDDO
C=Xorder
IF(C .EQ. 2)THEN
   DO I = 1, Nprim
      M(I,2,0)=Xpq(I)*M(I,1,0)+P(I)*M(I,1,1)     
   ENDDO
   IF(LIMIT .GT. 0)THEN
      DO I = 1, Nprim
         M(I,2,1)=M(I,1,0)+Xpq(I)*M(I,1,1)          
      ENDDO
      IF(LIMIT .GT. 1)THEN
         DO I = 1, Nprim
            M(I,2,2)=2*M(I,1,1)                      
         ENDDO
      ENDIF
   ENDIF
ELSE
   DO I = 1, Nprim
      M(I,C,0)= Xpq(I)*M(I,C-1,0)+P(I)*M(I,C-1,1)     
   ENDDO
   DO T = 1,MIN(LIMIT,C-2)
      DO I = 1, Nprim
         M(I,C,T)=T*M(I,C-1,T-1)+Xpq(I)*M(I,C-1,T)+P(I)*M(I,C-1,T+1)
      ENDDO
   ENDDO
   IF(LIMIT .GT. C-2)THEN
      !        T=C-1
      DO I = 1, Nprim
         M(I,C,C-1)=(C-1)*M(I,C-1,C-2)+Xpq(I)*M(I,C-1,C-1)
      ENDDO
      IF(LIMIT .GT. C-1)THEN
         !        T=C
         DO I = 1, Nprim
            M(I,C,C)=C*M(I,C-1,C-1)
         ENDDO
      ENDIF
   ENDIF
ENDIF

END SUBROUTINE STANDARD_LOOPX

SUBROUTINE ORDER_EQ_ONE_LOOP(INPUT,M,Xpq,Xorder,nPrim,LIMIT,lupri)
!written for X loop (therefore Xpq and Xorder) but general for any direction
IMPLICIT NONE
INTEGER                 :: nPrim,LIMIT,C,Xorder,T,lupri,I
real(realk)             :: Xpq(nPrim),M(nPrim,0:Xorder,0:Xorder),INPUT(nPrim)
  CALL DCOPY(nPrim,INPUT,1,M(1,0,0),1)
  DO I = 1, Nprim
     M(I,1,0)=Xpq(I)*INPUT(I)  
  ENDDO
  IF(LIMIT .GT. 0)THEN
     DO I = 1, Nprim
        M(I,1,1)=INPUT(I)                        
     ENDDO
  ENDIF
END SUBROUTINE ORDER_EQ_ONE_LOOP

RECURSIVE SUBROUTINE Qsort(a,INDEXES)
IMPLICIT NONE
REAL(REALK), INTENT(INOUT) :: a(:)
INTEGER :: INDEXES(:)
INTEGER :: split

  IF(size(a) > 7) THEN 
     CALL Partition(a,INDEXES,split)
     CALL Qsort(a(:split-1),INDEXES(:split-1))
     CALL Qsort(a(split:),INDEXES(split:))
  ELSEIF(size(a) > 1)THEN
     CALL INSERTION(a,INDEXES)
  END IF
 
END SUBROUTINE Qsort

SUBROUTINE INSERTION(a,INDEXES)
IMPLICIT NONE
REAL(REALK), INTENT(INOUT) :: a(:)
INTEGER                    :: I,J, INDEXES(:),tempI
real(realk) :: temp

Do I = 2, Size(a)
   temp = a(I)
   tempI = INDEXES(I)
   DO J = I-1, 1, -1
      IF(temp .GT. a(J)) Then
         a(J+1) = a(J)
         INDEXES(J+1) = INDEXES(J)
      ELSE
         Exit
      ENDIF
   ENDDO
   a(J+1) = temp
   INDEXES(J+1) = tempI
End Do

END SUBROUTINE INSERTION

SUBROUTINE Partition(a, INDEXES, marker)
IMPLICIT NONE
REAL(REALK), INTENT(INOUT) :: a(:)
INTEGER, INTENT(OUT)       :: marker
INTEGER                    :: left, right, INDEXES(:),tempI,I
real(realk) :: pivot, temp
 
pivot = (a(1) + a(size(a))) / 2 ! Average of first and last elements 
                                !to prevent quadratic behavior with
                                ! sorted or reverse sorted data
IF(a(1) .EQ. pivot)THEN
   DO I =size(a)-1,2,-1
      IF(a(I) .NE. pivot)THEN
         pivot = (a(1) + a(I)) / 2
         EXIT
      ENDIF
   ENDDO
ENDIF

left = 0                         ! behavior with sorted or reverse sorted data
right = size(a) + 1

DO WHILE (left < right)
   right = right - 1
   DO WHILE (a(right) < pivot)
      right = right-1
   END DO
   left = left + 1
   DO WHILE (a(left) > pivot)
      left = left + 1
   END DO
   IF (left < right) THEN 
      temp = a(right)
      a(right) = a(left)
      a(left) = temp
      tempI = INDEXES(right)
      INDEXES(right) = INDEXES(left)
      INDEXES(left) = tempI
   END IF
END DO

IF (left == right) THEN
   marker = left + 1
ELSE
   marker = left
END IF

END SUBROUTINE Partition

!!$ WILL SORT THE SMALLEST FIRST
!!$SUBROUTINE Partition(a, INDEXES, marker)
!!$IMPLICIT NONE
!!$REAL(REALK), INTENT(INOUT) :: a(:)
!!$INTEGER, INTENT(OUT)       :: marker
!!$INTEGER                    :: left, right, INDEXES(:),tempI
!!$real(realk) :: pivot, temp
!!$ 
!!$  pivot = (a(1) + a(size(a))) / 2  ! Average of first and last elements to prevent quadratic 
!!$  left = 0                         ! behavior with sorted or reverse sorted data
!!$  right = size(a) + 1
!!$
!!$  DO WHILE (left < right)
!!$     right = right - 1
!!$     DO WHILE (a(right) > pivot)
!!$        right = right-1
!!$     END DO
!!$     left = left + 1
!!$     DO WHILE (a(left) < pivot)
!!$        left = left + 1
!!$     END DO
!!$     IF (left < right) THEN 
!!$        temp = a(left)
!!$        a(left) = a(right)
!!$        a(right) = temp
!!$        tempI = INDEXES(left)
!!$        INDEXES(left) = INDEXES(right)
!!$        INDEXES(right) = tempI
!!$     END IF
!!$  END DO
!!$
!!$  IF (left == right) THEN
!!$     marker = left + 1
!!$  ELSE
!!$     marker = left
!!$  END IF
!!$
!!$END SUBROUTINE Partition

SUBROUTINE initAlloc(Alloc,LUPRI,IPRINT,SIDE)
implicit none
TYPE(Allocitem)    :: Alloc
Integer            :: LUPRI,IPRINT
Character*(*)      :: Side
IF (SIDE.EQ.'Both') THEN
   Alloc%maxPrimLHS    = 0
   Alloc%maxPrimTUVLHS = 0
   Alloc%maxPrimLHSpass    = 0
   Alloc%maxPrimTUVLHSpass = 0
   Alloc%maxnAngLHS     = 0
   Alloc%maxContLHS    = 0
   Alloc%maxAngmomLHS  = 0
   Alloc%maxAngmomLHS = 0
   Alloc%maxTUVLHS = 0
   Alloc%maxPrimTUVijkLHS = 0
   Alloc%maxTotOrbLHS = 0
   Alloc%maxijkLHS = 0
   Alloc%maxPrimRHS    = 0
   Alloc%maxPrimTUVRHS = 0
   Alloc%maxPrimRHSpass    = 0
   Alloc%maxPrimTUVRHSpass = 0
   Alloc%maxnAngRHS     = 0
   Alloc%maxContRHS    = 0
   Alloc%maxAngmomRHS  = 0
   Alloc%maxAngmomRHS = 0
   Alloc%maxTUVRHS = 0
   Alloc%maxPrimTUVijkRHS = 0
   Alloc%maxTotOrbRHS = 0
   Alloc%maxijkRHS = 0
ELSE IF (SIDE.EQ.'LHS') THEN
   Alloc%maxPrimLHS    = 0
   Alloc%maxPrimTUVLHS = 0
   Alloc%maxPrimLHSpass    = 0
   Alloc%maxPrimTUVLHSpass = 0
   Alloc%maxnAngLHS     = 0
   Alloc%maxContLHS    = 0
   Alloc%maxAngmomLHS  = 0
   Alloc%maxAngmomLHS = 0
   Alloc%maxTUVLHS = 0
   Alloc%maxPrimTUVijkLHS = 0
   Alloc%maxTotOrbLHS = 0
   Alloc%maxijkLHS = 0
ELSE IF (SIDE.EQ.'RHS') THEN
   Alloc%maxPrimRHS    = 0 !used for integrals
   Alloc%maxPrimTUVRHS = 0 !used for integrals
   Alloc%maxPrimRHSpass    = 0 !used for integrals
   Alloc%maxPrimTUVRHSpass = 0 !used for integrals
   Alloc%maxnAngRHS     = 0
   Alloc%maxContRHS    = 0
   Alloc%maxAngmomRHS  = 0
   Alloc%maxAngmomRHS = 0
   Alloc%maxTUVRHS = 0
   Alloc%maxPrimTUVijkRHS = 0
   Alloc%maxTotOrbRHS = 0
   Alloc%maxijkRHS = 0
ELSE
  WRITE(LUPRI,'(1X,2A)') 'Programming error. Side =',SIDE
  CALL QUIT('Programming error. Side not an option in initIntegral')
ENDIF
END SUBROUTINE initAlloc

SUBROUTINE allocateIntegrals(PQ,Integral,Alloc,maxpasses)
implicit none
TYPE(Integrand)    :: PQ
TYPE(Integralitem) :: Integral
TYPE(Allocitem)    :: Alloc
Integer            :: maxPasses
!
Integer :: maxTUVdim,maxPrim

Nullify(PQ%distance) 
Nullify(PQ%squaredDistance)
Nullify(PQ%exponents)
Nullify(PQ%reducedExponents)
Nullify(PQ%integralPrefactor)
Nullify(PQ%iprimP)
Nullify(PQ%iprimQ)

IF(Alloc%maxPrimRHSpass .NE. 0)THEN
   maxPrim   = Alloc%maxPrimLHS*Alloc%maxPrimRHSpass
ELSE
   maxPrim   = Alloc%maxPrimLHS*Alloc%maxPrimRHS*maxPasses
ENDIF
Allocate(PQ%distance(3,maxPrim))
Allocate(PQ%squaredDistance(maxPrim))
Allocate(PQ%exponents(maxPrim))
Allocate(PQ%reducedExponents(maxPrim))
Allocate(PQ%integralPrefactor(maxPrim))
Allocate(PQ%iprimP(maxPrim))
Allocate(PQ%iprimQ(maxPrim))

NULLIFY(Integral%IN)
NULLIFY(Integral%OUT)
NULLIFY(Integral%RTUV)
NULLIFY(Integral%integrals)
NULLIFY(Integral%PQ)

IF(Alloc%maxPrimTUVRHSpass .NE. 0) THEN
    maxTUVdim = Alloc%maxPrimTUVRHSpass*Alloc%maxPrimTUVLHS
ELSE
    maxTUVdim = Alloc%maxPrimTUVRHS*Alloc%maxPrimTUVLHS*maxPasses
ENDIF
ALLOCATE(Integral%IN(maxTUVdim))
ALLOCATE(Integral%OUT(maxTUVdim))
ALLOCATE(Integral%RTUV(maxTUVdim))
ALLOCATE(Integral%integrals(maxTUVdim))
ALLOCATE(Integral%PQ(Alloc%maxTotOrbRHS*maxpasses,Alloc%maxTotOrbLHS))

END SUBROUTINE allocateIntegrals

SUBROUTINE deallocateIntegrals(PQ,Integral)
IMPLICIT NONE
TYPE(Integrand)    :: PQ
TYPE(Integralitem) :: Integral
 
deallocate(PQ%distance)
deallocate(PQ%squaredDistance)
deallocate(PQ%exponents)
deallocate(PQ%reducedExponents)
deallocate(PQ%integralPrefactor)
deallocate(PQ%iprimP)
deallocate(PQ%iprimQ)
deallocate(Integral%IN)
deallocate(Integral%OUT)
deallocate(Integral%RTUV)
deallocate(Integral%integrals)
DEALLOCATE(Integral%PQ)
Nullify(PQ%distance)
Nullify(PQ%squaredDistance)
Nullify(PQ%exponents)
Nullify(PQ%reducedExponents)
Nullify(PQ%integralPrefactor)
Nullify(PQ%iprimP)
Nullify(PQ%iprimQ)
NULLIFY(Integral%IN)
NULLIFY(Integral%OUT)
NULLIFY(Integral%RTUV)
NULLIFY(Integral%integrals)
NULLIFY(Integral%PQ)

END SUBROUTINE deallocateIntegrals

SUBROUTINE swapRealPointers(new,old)
implicit none
Real(realk),pointer :: new(:),old(:),temp(:)
temp => old
old  => new
new  => temp
END SUBROUTINE swapRealPointers

SUBROUTINE SelectODPassTypes(ODpassesIndex,P,nBatches,nPassTypes,IPRINT,LUPRI)
implicit none
Integer       :: nBatches,nPassTypes,IPRINT,LUPRI
Integer       :: ODpassesIndex(nBatches)
Type(Overlap) :: P(nBatches)

Integer       :: I,J,nPrim,nCont,numAng1,minAng1,maxAng1,numAng2,minAng2,maxAng2
Integer       :: nPrim1,maxCont1,nPrim2,maxCont2
!Integer(long) :: UniqeIdentfiers(nBatches)
Character(len=57)  :: UniqeIdentfiers(nBatches)
Integer       :: angmomIdentifier1,angmomIdentifier2,CCidentifier1,CCidentifier2
Logical       :: spher1,spher2,hermiteSingle,unique,sameAO
Integer       :: EXP

nPassTypes = 0
DO I=1,nBatches
  nPrim   = P(I)%nPrimitives
  IF(nPrim.GT.0)THEN
     !ToDo select according to the same iprim1 and iprim2 set-ups 
     !(used for primitive screening)
     numAng1 = P(I)%orbital1%nAngmom
     numAng2 = P(I)%orbital2%nAngmom
     angmomIdentifier1 = 0
     angmomIdentifier2 = 0
     DO J=1,numAng1
        angmomIdentifier1 = angmomIdentifier1 + 2**P(I)%orbital1%angmom(J)
     ENDDO
     DO J=1,numAng2
        angmomIdentifier2 = angmomIdentifier2 + 2**P(I)%orbital2%angmom(J)
     ENDDO
  
     CCidentifier1 = P(I)%orbital1%CCidentifier
     CCidentifier2 = P(I)%orbital2%CCidentifier
     spher1        = P(I)%orbital1%spherical
     spher2        = P(I)%orbital2%spherical
     hermiteSingle = P(I)%hermiteSingle
     sameAO        = P(I)%sameAO
     
     IF (nPrim.GT.999999999) THEN
        CALL QUIT('nPrim>999999999 in SelectODPassTypes')
     ELSE IF (angmomIdentifier1.GT.999999999) THEN
        CALL QUIT('angmomIdentifier1>999999999 in SelectODPassTypes')
     ELSE IF (angmomIdentifier2.GT.999999999) THEN
        CALL QUIT('angmomIdentifier2>999999999 in SelectODPassTypes')
     ELSE IF (CCidentifier1.GT.999999999) THEN
        CALL QUIT('CCidentifier1>999999999 in SelectODPassTypes')
     ELSE IF (CCidentifier2.GT.999999999) THEN
        CALL QUIT('CCidentifier2>999999999 in SelectODPassTypes')
     ENDIF
     WRITE(UniqeIdentfiers(I),'(5I9,4L3)') nPrim,angmomIdentifier2,angmomIdentifier1,&
          & CCidentifier1,CCidentifier2,spher1,spher2,hermiteSingle,sameAO
     unique = .TRUE.
     DO J=I-1,1,-1
        IF (UniqeIdentfiers(I).EQ.UniqeIdentfiers(J)) THEN
           ODpassesIndex(I)=ODpassesIndex(J)
           unique = .FALSE.
           EXIT
        ENDIF
     ENDDO
     IF (unique) THEN
        nPassTypes = nPassTypes + 1
        ODpassesIndex(I) = nPassTypes
     ENDIF
  ENDIF
ENDDO

IF (IPRINT.GT.0) THEN
  CALL LSHEADER(LUPRI,'Output from SelectODPassTypes')
  WRITE(LUPRI,'(1X,A,I5)') 'Number of OD-batches', nBatches
  WRITE(LUPRI,'(1X,A,I5)') 'Number of pass-types', nPassTypes
  IF (IPRINT.GT.5) THEN
    WRITE(LUPRI,'(3X,A)') 'Batch Pass      nPrim  angInd1  angInd2   CCind1    CCind2 s1 s2 hS AO'
    DO I=1,nBatches 
      WRITE(LUPRI,'(3X,2I5,2X,1A57)') I,ODpassesIndex(I),UniqeIdentfiers(I)
    ENDDO
  ENDIF
ENDIF

END SUBROUTINE SelectODPassTypes

SUBROUTINE SelectODPassTypesFromODbatch(ODpassesIndex,nPrim,ODB,nBatches,&
     &nPassTypes,maxpasses,maxpassfortype,INPUT,IPRINT,LUPRI,SIDE)
implicit none
Integer       :: nBatches,nPassTypes,IPRINT,LUPRI
Integer       :: ODpassesIndex(nBatches),maxpasses,nPrim(nBatches)
!Type(Overlap) :: P(nBatches)
TYPE(ODITEM)        :: ODB
TYPE(IntegralInput) :: Input
Character*(*)       :: side
Integer,pointer :: Maxpassfortype(:)

Integer       :: I,J,nCont,numAng1,minAng1,maxAng1,numAng2,minAng2,maxAng2
Integer       :: nPrim1,maxCont1,nPrim2,maxCont2
!Integer(long) :: UniqeIdentfiers(nBatches)
Character(len=57)  :: UniqeIdentfiers(nBatches)
Integer       :: angmomIdentifier1,angmomIdentifier2,CCidentifier1,CCidentifier2
Logical       :: spher1,spher2,hermiteSingle,unique,sameAO
Integer       :: EXP,i1,i2,start2,a1,a2,a,b
real(realk)   :: MAXGAB,MAXELM
Real(realk),pointer :: GAB(:,:),primGAB(:,:) 
Character(len=80)   :: t1,t2

IF (Side.EQ.'LHS') THEN
  GAB => INPUT%pGAB_LHS
  MAXELM=INPUT%PS_MAXELM_RHS
ELSEIF (Side.EQ.'RHS') THEN
  GAB => INPUT%pGAB_RHS
  MAXELM=INPUT%PS_MAXELM_LHS
ELSE
  WRITE(*,'(1X,2A)') 'Error in GET_NPRIMITIVES. Side =',Side
  CALL QUIT('Programming error. Wrong Side in GET_NPRIMITIVES')
ENDIF

nPassTypes = 0
DO I=1,nBatches
   !DETERMINE number of Primitives
   CALL GET_NPRIM_FROM_ODBATCH(nPrim(I),ODB%BATCH(I),INPUT,SIDE)  
  
  IF(nPrim(I).GT.0)THEN
     !ToDo select according to the same iprim1 and iprim2 set-ups 
     !(used for primitive screening)
     numAng1 = ODB%BATCH(I)%AO(1)%p%nAngmom
     numAng2 = ODB%BATCH(I)%AO(2)%p%nAngmom
     angmomIdentifier1 = 0
     angmomIdentifier2 = 0
     DO J=1,numAng1
        angmomIdentifier1 = angmomIdentifier1 + 2**ODB%BATCH(I)%AO(1)%p%angmom(J)
     ENDDO
     DO J=1,numAng2
        angmomIdentifier2 = angmomIdentifier2 + 2**ODB%BATCH(I)%AO(2)%p%angmom(J)
     ENDDO
  
     CCidentifier1 = ODB%BATCH(I)%AO(1)%p%CCindex(1)
     CCidentifier2 = ODB%BATCH(I)%AO(2)%p%CCindex(1)
     spher1        = ODB%BATCH(I)%AO(1)%p%spherical
     spher2        = ODB%BATCH(I)%AO(2)%p%spherical

     hermiteSingle = .FALSE. !difficult
     t1 = ODB%BATCH(I)%AO(1)%p%type
     t2 = ODB%BATCH(I)%AO(2)%p%type
     IF( (t1.EQ.'Hermite'.AND.t2.EQ.'Empty').OR.(t1.EQ.'Empty'  .AND.t2.EQ.'Hermite') ) THEN
        hermiteSingle = .TRUE.
     ENDIF

     sameAO        = ODB%BATCH(I)%sameAO
     
     IF (nPrim(I).GT.999999999) THEN
        CALL QUIT('nPrim>999999999 in SelectODPassTypes')
     ELSE IF (angmomIdentifier1.GT.999999999) THEN
        CALL QUIT('angmomIdentifier1>999999999 in SelectODPassTypes')
     ELSE IF (angmomIdentifier2.GT.999999999) THEN
        CALL QUIT('angmomIdentifier2>999999999 in SelectODPassTypes')
     ELSE IF (CCidentifier1.GT.999999999) THEN
        CALL QUIT('CCidentifier1>999999999 in SelectODPassTypes')
     ELSE IF (CCidentifier2.GT.999999999) THEN
        CALL QUIT('CCidentifier2>999999999 in SelectODPassTypes')
     ENDIF
     WRITE(UniqeIdentfiers(I),'(5I9,4L3)') nPrim(I),angmomIdentifier2,angmomIdentifier1,&
          & CCidentifier1,CCidentifier2,spher1,spher2,hermiteSingle,sameAO
     unique = .TRUE.
     DO J=I-1,1,-1
        IF (UniqeIdentfiers(I).EQ.UniqeIdentfiers(J)) THEN
           ODpassesIndex(I)=ODpassesIndex(J)
           unique = .FALSE.
           EXIT
        ENDIF
     ENDDO
     IF (unique) THEN
        nPassTypes = nPassTypes + 1
        ODpassesIndex(I) = nPassTypes
     ENDIF
  ELSE
     ODpassesIndex(I) = 0
  ENDIF
ENDDO

NULLIFY(Maxpassfortype)
ALLOCATE(Maxpassfortype(nPassTypes))
DO I=1,nPassTypes
   Maxpassfortype(I)=1
ENDDO
DO I=1,nBatches
   unique = .TRUE.
   DO J=I-1,1,-1
      IF(ODpassesIndex(I).EQ. 0)EXIT
      IF (UniqeIdentfiers(I).EQ.UniqeIdentfiers(J)) THEN
         Maxpassfortype(ODpassesIndex(I))=Maxpassfortype(ODpassesIndex(I))+1
         EXIT
      ENDIF
   ENDDO
ENDDO
DO I=1,nPassTypes
   Maxpassfortype(I)=MIN(maxpasses,Maxpassfortype(I))
ENDDO

IF (IPRINT.GT.0) THEN
  CALL LSHEADER(LUPRI,'Output from SelectODPassTypes')
  WRITE(LUPRI,'(1X,A,I5)') 'Number of OD-batches', nBatches
  WRITE(LUPRI,'(1X,A,I5)') 'Number of pass-types', nPassTypes
  IF (IPRINT.GT.5) THEN
    WRITE(LUPRI,'(3X,A)') 'Batch Pass      nPrim  angInd1  angInd2   CCind1    CCind2 s1 s2 hS AO'
    DO I=1,nBatches 
      WRITE(LUPRI,'(3X,2I5,2X,1A57)') I,ODpassesIndex(I),UniqeIdentfiers(I)
    ENDDO
  ENDIF
ENDIF

END SUBROUTINE SelectODPassTypesFromODbatch

SUBROUTINE GET_NPRIM_FROM_ODBATCH(nPrim,ODB,input,side)
IMPLICIT NONE
INTEGER :: nPrim
Character*(*)       :: side
TYPE(IntegralInput) :: Input
TYPE(ODBATCH)       :: ODB
!
integer             :: i1,start2,i2,a1,a2,a,b
Character(len=80)   :: t1,t2
real(realk)         :: maxgab,MAXELM
Real(realk),pointer :: GAB(:,:)

IF (Side.EQ.'LHS') THEN
  GAB => INPUT%pGAB_LHS
  MAXELM=INPUT%PS_MAXELM_RHS
ELSEIF (Side.EQ.'RHS') THEN
  GAB => INPUT%pGAB_RHS
  MAXELM=INPUT%PS_MAXELM_LHS
ELSE
  WRITE(*,'(1X,2A)') 'Error in GET_NPRIM_FROM_ODBATCH. Side =',Side
  CALL QUIT('Programming error. Wrong Side in GET_NPRIM_FROM_ODBATCH')
ENDIF

nPrim = 0
t1 = ODB%AO(1)%p%type
t2 = ODB%AO(2)%p%type
IF (.NOT. (t1.EQ.'Empty'.AND.t2.EQ.'Empty') ) THEN
   DO i1=1,ODB%AO(1)%p%nPrimitives
      start2 = 1
      IF (ODB%sameAO) start2 = i1 
      DO i2=start2,ODB%AO(2)%p%nPrimitives
         !   Find maxgab for given primitive
         maxgab = 0.0d0
         IF(INPUT%PS_SCREEN)THEN
            DO a1=1,ODB%AO(1)%p%nAngmom
               DO a2=1,ODB%AO(2)%p%nAngmom
                  a = ODB%AO(1)%p%startprimOrbital(a1)+(i1-1)*ODB%AO(1)%p%nOrbComp(a1)
                  b = ODB%AO(2)%p%startprimOrbital(a2)+(i2-1)*ODB%AO(2)%p%nOrbComp(a2)
                  maxgab=MAX(maxgab,GAB(a,b))
               ENDDO
            ENDDO
         ENDIF
         !   Add if maxgab greater than threshold
         IF(INPUT%PS_SCREEN)THEN
            IF( MAXGAB .GT. INPUT%PS_THRESHOLD/MAXELM)THEN
               nPrim = nPrim + 1
            ENDIF
         ELSE
            nPrim = nPrim + 1
         ENDIF
      ENDDO
   ENDDO
ELSE
   nPrim = 1
ENDIF
END SUBROUTINE GET_NPRIM_FROM_ODBATCH

SUBROUTINE BUILD_ERF_RJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,OMEGA)
IMPLICIT NONE
INTEGER         :: nPrim,LUPRI,IPRINT,Jmax,ZeroX,ZeroY,ZeroZ
REAL(REALK)     :: RJ000(0:Jmax,nPrim),OMEGA,X0,Y0,Z0
REAL(REALK)     :: Prefactor(nPrim),alpha(nPrim)
TYPE(Integrand) :: PQ
TYPE(integralitem) :: integral

CALL DCOPY(nPrim,PQ%integralPrefactor,1,Prefactor,1)
CALL DCOPY(nPrim,PQ%reducedExponents,1,alpha,1)
!COPY BECAUSE BUILD_ERF2_RJ000 changes the values of prefactor and ALPHA
CALL BUILD_ERF2_RJ000(OMEGA,RJ000,nPrim,Jmax,LUPRI,IPRINT,&
     &alpha,PQ%squaredDistance,&
     &Prefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)

END SUBROUTINE BUILD_ERF_RJ000

SUBROUTINE BUILD_ERF2_RJ000(OMEGA,RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
!  Two-electron integrals over the operator: erf(omega r12)/r12 
!  with beta = omega^2/(alpha+omega^2) 
IMPLICIT NONE
INTEGER         :: nPrim,Jmax,Lupri,Iprint,nTABFJW1,nTABFJW2
REAL(REALK)     :: RJ000(0:Jmax,nPrim),Prefactor(nPrim)
!REAL(REALK)     :: RJ000(nPrim,0:Jmax),Prefactor(nPrim)
REAL(REALK)     :: alpha(nPrim),R2(nprim)
REAL(REALK)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)
INTEGER         :: nprim1,nprim2,nprim3,INDADS(nPrim,3),I,NODS

REAL(REALK)     :: D2JP36,WVALU(Nprim),WVAL,WVALS(Nprim,3)
REAL(REALK),PARAMETER :: HALF =0.5d0,D1=1.d0,D2 = 2.D0, D4 = 4.D0, D10=10.d0
Real(realk),parameter :: D12 = 12.D0, TENTH = 0.1D0
REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6.D0, COEF4 = D1/24.D0
REAL(REALK), PARAMETER :: COEF5 = - D1/120.D0, COEF6 = D1/720.D0
REAL(REALK), PARAMETER :: SMALL = 0.000001D0
Integer :: IPNT,J
Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA
REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092D0
REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686D0
REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909D0
REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346D0
Real(realk), parameter :: PI=3.14159265358979323846D0
REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730D00
REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
REAL(REALK)     :: FACINT(nPRim),prefac(Nprim),prefac2(Nprim)
REAL(REALK)     :: OMEGA,OMEGA2,BETAAT,SQBETA,FAC,THRSH

Real(realk) :: W2,W3,W4,W5,W6,R
LOGICAL :: maxJgt0

OMEGA2=OMEGA*OMEGA
NODS  = 0
DO I = 1, NPrim
   !      IF (ABS(Prefactor(I)) .GT. THRSH) THEN
   NODS = NODS + 1
   BETAAT = OMEGA2 / (ALPHA(I) + OMEGA2)
   SQBETA = SQRT(BETAAT)
   Prefactor(I) = Prefactor(I)*SQBETA
   ALPHA(I) = ALPHA(I)*BETAAT         
   WVALU(NODS) = ALPHA(I)*R2(I)
   !         INDADR(NODS) = I
   !      ELSE
   !         WRITE(LUPRI,*)'no done iPrim',I
   !         Prefactor(I) = D0
   !         DO J = 0, JMAX
   !            RJ000(I,J) = D0
   !         END DO
   !      END IF
ENDDO

!  Calculate gamma function
!  ========================
!
Nprim1 = 0
Nprim2 = 0
Nprim3 = 0
D2JP36 = 2*JMAX + 36

DO I = 1, NODS
   WVAL = WVALU(I)
   !  0 < WVAL < 0.000001
   IF (WVAL .LT. SMALL) THEN         
      RJ000(0,I) = D1
      DO J=1,JMAX
         RJ000(J,I)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
      ENDDO
      !  0 < WVAL < 12 
   ELSE IF (WVAL .LT. D12) THEN
      IPNT = NINT(D10*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W4    = W3*WDIFF
      W5    = W4*WDIFF
      W6    = W5*WDIFF
      R = TABFJW(JMAX,IPNT)
      R = R -TABFJW(JMAX+1,IPNT)*WDIFF
      R = R + COEF2*TABFJW(JMAX+2,IPNT)*W2
      R = R + COEF3*TABFJW(JMAX+3,IPNT)*W3
      R = R + COEF4*TABFJW(JMAX+4,IPNT)*W4
      R = R + COEF5*TABFJW(JMAX+5,IPNT)*W5
      R = R + COEF6*TABFJW(JMAX+6,IPNT)*W6
      
      RJ000(JMAX,I) = R
      IF (JMAX.GT.0) THEN
         REXPW = HALF*EXP(-WVAL)
         DO J=JMAX-1,0,-1
            RJ000(J,I) = (WVAL*RJ000(J+1,I) + REXPW)*D2/(2*J + 1)
         ENDDO
      ENDIF
      !  12 < WVAL <= (2J+36) 
   ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = HALF*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0,I) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      DO J=1,JMAX
         RJ000(J,I) = RWVAL*((J - HALF)*RJ000(J-1,I)-REXPW)
      ENDDO
      !  (2J+36) < WVAL 
   ELSE
      RWVAL = PID4/WVAL
      RJ000(0,I) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      DO J = 1, JMAX
         RJ000(J,I) = RWVAL*(J - HALF)*RJ000(J-1,I)
      ENDDO
   END IF
ENDDO

! Scaling
DO I=1,nPrim
  PREF = Prefactor(I)
  RJ000(0,I) = PREF*RJ000(0,I)
  IF (jmax.GT.0) THEN
    D2MALPHA = -2*alpha(I)
    DO j=1,jmax
      PREF = PREF*D2MALPHA
      RJ000(J,I) = PREF*RJ000(J,I)
    ENDDO
  ENDIF
ENDDO
!
!     *************************
!     ***** PRINT SECTION *****
!     *************************
!
IF (IPRINT .GT. 10) THEN
   CALL LSHEADER(lupri,'Output from Erf_RJ000 ')
   WRITE (LUPRI,'(2X,A,F16.8)') 'Omega2:  ', OMEGA2
   WRITE (LUPRI,'(2X,A,I10)') 'JMAX:  ', JMAX
   WRITE (LUPRI,'(2X,A,I10)') 'NPrim:  ', NPrim
   IF (IPRINT .GT. 20) THEN
      CALL LSHEADER(lupri,'Scaled incomplete gamma function in ERF_RJ000')
      CALL OUTPUT(RJ000,1,JMAX+1,1,NPrim,JMAX+1,NPrim,1,LUPRI)
   END IF
END IF

END SUBROUTINE BUILD_ERF2_RJ000

SUBROUTINE SET_ALLOC(Alloc,Input,ODbat,SIDE,IPRINT,LUPRI)
Implicit none
INTEGER             :: nbatches,lupri,IPRINT,maxpasses
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
TYPE(ODITEM)        :: ODbat
Character*(*)       :: side
!Integer,optional    :: maxpassfortype(nPasstype)
!TYPE(ODBATCH)       :: ODbat(nbatches)
!
Integer             :: i1,i2,l1,l2,ijk1,ijk2,ijk,maxijk,I,nprim
Integer             :: maxangmom,totOrbitals,maxtotorbitals,NDERIV
Integer             :: maxprim,endangmom,maxtuv,maxcont,maxnangmom
Integer             :: MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK,MAXnTUV
Integer             :: MAXPRIMTUVIJK,maxA,maxB,start,maxijk2
LOGICAL             :: LHS,hermiteSingle 

SELECT CASE(SIDE)
CASE('LHS')
   LHS=.TRUE.
   NDERIV=INPUT%NDERIVP
CASE('RHS')
   LHS=.FALSE.
   NDERIV=INPUT%NDERIVQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_OVERLAP =',SIDE
   CALL QUIT('Wrong case in SET_OVERLAP')
END SELECT

maxPrim = 0
maxAngmom = 0
maxnAngmom = 0
maxCont = 0
maxijk = 0
MAXPRIMTUV = 0
MAXPRIMIJK = 0
MAXnTUV = 0
MAXPRIMTUVIJK = 0
maxA = 0
maxB = 0
maxTotorb=0
maxijk2=0

DO I=1,ODbat%nbatches
   CALL GET_NPRIM_FROM_ODBATCH(nPrim,ODBat%BATCH(I),INPUT,SIDE)
   maxPrim = MAX(maxPrim,nPrim)
   maxnAngmom = MAX(maxnAngmom,ODbat%BATCH(I)%nAngmom)
   maxCont = MAX(maxCont,ODbat%BATCH(I)%maxContracted)
   maxA = MAX(maxA,ODbat%BATCH(I)%AO(1)%p%maxAngmom)
   maxB = MAX(maxB,ODbat%BATCH(I)%AO(2)%p%maxAngmom)
   totOrbitals = 0
   maxijk = 0
   maxangmom = 0
   DO I1=1,ODbat%BATCH(I)%AO(1)%p%nAngmom
      start=1
      IF(ODbat%BATCH(I)%sameAO)start=I1
      DO I2=start,ODbat%BATCH(I)%AO(2)%p%nAngmom
         l1 = ODbat%BATCH(I)%AO(1)%p%angmom(I1)
         l2 = ODbat%BATCH(I)%AO(2)%p%angmom(I2)
         ijk1 = (l1+1)*(l1+2)/2
         ijk2 = (l2+1)*(l2+2)/2
         ijk  = ijk1*ijk2
         maxijk = max(maxijk,ijk)
         maxAngmom = max(maxAngmom,l1+l2)
         totOrbitals = totOrbitals + ODbat%BATCH(I)%AO(1)%p%nOrbitals(i1)*ODbat%BATCH(I)%AO(2)%p%nOrbitals(i2)*NDERIV
      ENDDO
   ENDDO
   MAXPRIMIJK = MAX(MAXPRIMIJK,maxijk*nprim)
   maxTotorb = MAX(maxTotorb,totOrbitals)
   endAngmom  = maxAngmom
   IF(INPUT%operator(1:7).EQ.'Kinetic'.AND.LHS)endAngmom  = maxAngmom + 2
   maxijk2 = MAX(maxijk2,maxijk)
   maxTUV = (endAngmom+1)*(endAngmom+2)*(endAngmom+3)/6
   MAXPRIMTUV = MAX(MAXPRIMTUV,maxTUV*nPrim)
   MAXnTUV=MAX(MAXnTUV,maxTUV)
   MAXPRIMTUVIJK = MAX(MAXPRIMTUVIJK,maxijk*nPrim*maxTUV*ODbat%BATCH(I)%nAngmom) 
ENDDO

IF(LHS)THEN
   Alloc%maxPrimLHS    = maxprim
   Alloc%maxPrimTUVLHS = max(MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK)
   Alloc%maxnAngLHS     = maxnAngmom
   Alloc%maxContLHS    = maxCont
   Alloc%maxAngmomLHS  = max(maxA,maxB,Alloc%maxAngmomLHS)
   Alloc%maxTUVLHS        = max(maxTUV,Alloc%maxTUVLHS)
   Alloc%maxPrimTUVijkLHS = max(MAXPRIMTUVIJK,Alloc%maxPrimTUVijkLHS) 
   Alloc%maxTotOrbLHS = max(maxTotorb,Alloc%maxTotOrbLHS)
   Alloc%maxijkLHS = max(maxijk2,Alloc%maxijkLHS)
ELSE
   Alloc%maxPrimRHS    = maxprim
   Alloc%maxPrimTUVRHS = max(MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK)
   Alloc%maxnAngRHS     = maxnAngmom
   Alloc%maxContRHS    = maxCont
   Alloc%maxAngmomRHS  = max(maxA,maxB,Alloc%maxAngmomRHS)
   Alloc%maxTUVRHS        = max(maxTUV,Alloc%maxTUVRHS)
   Alloc%maxPrimTUVijkRHS = max(MAXPRIMTUVIJK,Alloc%maxPrimTUVijkRHS) 
   Alloc%maxTotOrbRHS = max(maxTotorb,Alloc%maxTotOrbRHS)
   Alloc%maxijkRHS = max(maxijk2,Alloc%maxijkRHS)
ENDIF

END SUBROUTINE SET_ALLOC

SUBROUTINE SET_ALLOC2(Alloc,Input,ODbat,SIDE,IPRINT,LUPRI,nbatches,ODpassesIndex,nPasstype,maxpassfortype,maxpasses)
Implicit none
INTEGER             :: nbatches,lupri,IPRINT,maxpasses,nPasstype
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
TYPE(ODITEM)        :: ODbat
Character*(*)       :: side
Integer             :: ODpassesIndex(nbatches)
Integer             :: maxpassfortype(nPasstype)
!TYPE(ODBATCH)       :: ODbat(nbatches)
!
Integer             :: i1,i2,l1,l2,ijk1,ijk2,ijk,maxijk,I,nprim
Integer             :: maxangmom,totOrbitals,maxtotorbitals,NDERIV
Integer             :: maxprim,endangmom,maxtuv,maxcont,maxnangmom
Integer             :: MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK,MAXnTUV
Integer             :: MAXPRIMTUVIJK,maxA,maxB,start,maxijk2,maxprimpass
LOGICAL             :: LHS,hermiteSingle,Dopass 
Integer             :: MAXPRIMIJKpass, MAXTOTORBpass, MAXPRIMTUVpass

!same as SET_ALLOC but this one also incoorporate the passes.
if(maxpasses .GT. 1)dopass = .TRUE.

SELECT CASE(SIDE)
CASE('LHS')
   LHS=.TRUE.
   NDERIV=INPUT%NDERIVP
CASE('RHS')
   LHS=.FALSE.
   NDERIV=INPUT%NDERIVQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_OVERLAP =',SIDE
   CALL QUIT('Wrong case in SET_OVERLAP')
END SELECT

maxPrim = 0
maxPrimpass = 0
maxAngmom = 0
maxnAngmom = 0
maxCont = 0
maxijk = 0
MAXPRIMTUV = 0
MAXPRIMIJK = 0
MAXnTUV = 0
MAXPRIMTUVIJK = 0
maxA = 0
maxB = 0
maxTotorb=0
maxijk2=0
MAXPRIMIJKpass = 0
MAXTOTORBpass = 0
MAXPRIMTUVpass = 0 

DO I=1,ODbat%nbatches
   CALL GET_NPRIM_FROM_ODBATCH(nPrim,ODBat%BATCH(I),INPUT,SIDE)
   maxPrim = MAX(maxPrim,nPrim)
   IF(ODpassesIndex(I).EQ.0)THEN
      maxPrimpass = MAX(maxPrimpass,nPrim)
   ELSE
      maxPrimpass = MAX(maxPrimpass,nPrim*maxpassfortype(ODpassesIndex(I)))
   ENDIF
   maxnAngmom = MAX(maxnAngmom,ODbat%BATCH(I)%nAngmom)
   maxCont = MAX(maxCont,ODbat%BATCH(I)%maxContracted)
   maxA = MAX(maxA,ODbat%BATCH(I)%AO(1)%p%maxAngmom)
   maxB = MAX(maxB,ODbat%BATCH(I)%AO(2)%p%maxAngmom)
   totOrbitals = 0
   maxijk = 0
   maxangmom = 0
   DO I1=1,ODbat%BATCH(I)%AO(1)%p%nAngmom
      start=1
      IF(ODbat%BATCH(I)%sameAO)start=I1
      DO I2=start,ODbat%BATCH(I)%AO(2)%p%nAngmom
         l1 = ODbat%BATCH(I)%AO(1)%p%angmom(I1)
         l2 = ODbat%BATCH(I)%AO(2)%p%angmom(I2)
         ijk1 = (l1+1)*(l1+2)/2
         ijk2 = (l2+1)*(l2+2)/2
         ijk  = ijk1*ijk2
         maxijk = max(maxijk,ijk)
         maxAngmom = max(maxAngmom,l1+l2)
         totOrbitals = totOrbitals + ODbat%BATCH(I)%AO(1)%p%nOrbitals(i1)*ODbat%BATCH(I)%AO(2)%p%nOrbitals(i2)*NDERIV
      ENDDO
   ENDDO
   MAXPRIMIJK = MAX(MAXPRIMIJK,maxijk*nprim)
   MAXTOTORB = MAX(MAXTOTORB,totOrbitals)
   endAngmom  = maxAngmom
   IF(INPUT%operator(1:7).EQ.'Kinetic'.AND.LHS)endAngmom  = maxAngmom + 2
   maxijk2 = MAX(maxijk2,maxijk)
   maxTUV = (endAngmom+1)*(endAngmom+2)*(endAngmom+3)/6
   MAXPRIMTUV = MAX(MAXPRIMTUV,maxTUV*nPrim)
   MAXnTUV=MAX(MAXnTUV,maxTUV)
   MAXPRIMTUVIJK = MAX(MAXPRIMTUVIJK,maxijk*nPrim*maxTUV*ODbat%BATCH(I)%nAngmom) 
   IF(DoPass)THEN
      IF(ODpassesIndex(I).EQ.0)THEN
      MAXPRIMIJKpass = MAX(MAXPRIMIJKpass,maxijk*nprim)
      MAXTOTORBpass = MAX(MAXTOTORBpass,totOrbitals)
      MAXPRIMTUVpass = MAX(MAXPRIMTUVpass,maxTUV*nPrim)
      ELSE
      MAXPRIMIJKpass = MAX(MAXPRIMIJKpass,maxijk*nprim*maxpassfortype(ODpassesIndex(I)))
      MAXTOTORBpass = MAX(MAXTOTORBpass,totOrbitals*maxpassfortype(ODpassesIndex(I)))
      MAXPRIMTUVpass = MAX(MAXPRIMTUVpass,maxTUV*nPrim*maxpassfortype(ODpassesIndex(I)))
      ENDIF
   ENDIF
ENDDO

IF(LHS)THEN
   Alloc%maxPrimLHS    = maxprim
   Alloc%maxPrimTUVLHS = max(MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK)
   Alloc%maxnAngLHS     = maxnAngmom
   Alloc%maxContLHS    = maxCont
   Alloc%maxAngmomLHS  = max(maxA,maxB,Alloc%maxAngmomLHS)
   Alloc%maxTUVLHS        = max(maxTUV,Alloc%maxTUVLHS)
   Alloc%maxPrimTUVijkLHS = max(MAXPRIMTUVIJK,Alloc%maxPrimTUVijkLHS) 
   Alloc%maxTotOrbLHS = max(maxTotorb,Alloc%maxTotOrbLHS)
   Alloc%maxijkLHS = max(maxijk2,Alloc%maxijkLHS)
   IF(DoPass)THEN
      Alloc%maxPrimTUVLHSpass = max(MAXPRIMTUVpass,MAXTOTORBpass,MAXPRIMIJKpass)
      Alloc%maxPrimLHSpass = maxprimpass
   ENDIF
ELSE
   Alloc%maxPrimRHS    = maxprim
   Alloc%maxPrimTUVRHS = max(MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK)
   Alloc%maxnAngRHS     = maxnAngmom
   Alloc%maxContRHS    = maxCont
   Alloc%maxAngmomRHS  = max(maxA,maxB,Alloc%maxAngmomRHS)
   Alloc%maxTUVRHS        = max(maxTUV,Alloc%maxTUVRHS)
   Alloc%maxPrimTUVijkRHS = max(MAXPRIMTUVIJK,Alloc%maxPrimTUVijkRHS) 
   Alloc%maxTotOrbRHS = max(maxTotorb,Alloc%maxTotOrbRHS)
   Alloc%maxijkRHS = max(maxijk2,Alloc%maxijkRHS)
   IF(DoPass)THEN
      Alloc%maxPrimTUVRHSpass = max(MAXPRIMTUVpass,MAXTOTORBpass,MAXPRIMIJKpass)
      Alloc%maxPrimRHSpass = maxprimpass
   ENDIF
ENDIF

END SUBROUTINE SET_ALLOC2

SUBROUTINE INIT_OVERLAP(P,Alloc,Input,ODbat,SIDE,IELECTRON,IPRINT,LUPRI)
Implicit none
INTEGER             :: nbatches,lupri,IPRINT,IELECTRON
TYPE(Overlap)       :: P
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
TYPE(ODITEM)        :: ODbat
!TYPE(Integralitem)  :: Integral
Character*(*)       :: side
!TYPE(ODBATCH)       :: ODbat(nbatches)
!
Character(len=80)   :: t1,t2
Integer             :: i1,i2,l1,l2,ijk1,ijk2,ijk,maxijk,I,nprim
Integer             :: maxangmom,totOrbitals,maxtotorbitals,NDERIV
Integer             :: maxprim,endangmom,maxtuv,maxcont,maxnangmom
Integer             :: MAXPRIMTUV,MAXPRIMTOTORB,MAXPRIMIJK,MAXnTUV
Integer             :: MAXPRIMTUVIJK,maxA,maxB
LOGICAL             :: LHS,hermiteSingle 

SELECT CASE(SIDE)
CASE('LHS')
   LHS=.TRUE.
   NDERIV=INPUT%NDERIVP
CASE('RHS')
   LHS=.FALSE.
   NDERIV=INPUT%NDERIVQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in INIT_OVERLAP =',SIDE
   CALL QUIT('Wrong case in INIT_OVERLAP')
END SELECT

IF(LHS)THEN
   maxprim = Alloc%maxPrimLHS
!   MAX(MAXPRIMTUV,MAXPRIMTOTORB,MAXPRIMIJK) = Alloc%maxPrimTUVLHS
   maxnAngmom = Alloc%maxnAngLHS
   maxCont = Alloc%maxContLHS
   maxA = Alloc%maxAngmomLHS
   maxTUV = Alloc%maxTUVLHS
   MAXPRIMTUVIJK = Alloc%maxPrimTUVijkLHS
   IF(maxnAngmom .EQ. 0)maxnAngmom=1
   IF(MAXPRIMTUVIJK .EQ. 0)MAXPRIMTUVIJK=1
ELSE
   maxprim = Alloc%maxPrimRHS
!   max(MAXPRIMTUV,MAXPRIMTOTORB,MAXPRIMIJK) = Alloc%maxPrimTUVRHS
   maxnAngmom = Alloc%maxnAngRHS
   maxCont = Alloc%maxContRHS
   maxA = Alloc%maxAngmomRHS
   maxTUV = Alloc%maxTUVRHS
   MAXPRIMTUVIJK = Alloc%maxPrimTUVijkRHS
   IF(maxnAngmom .EQ. 0)maxnAngmom=1
   IF(MAXPRIMTUVIJK .EQ. 0)MAXPRIMTUVIJK=1
ENDIF

CALL INIT_ORBITAL(P%orbital1,maxprim,maxnangmom,maxA)
CALL INIT_ORBITAL(P%orbital2,maxprim,maxnangmom,maxA)

 t1 = ODbat%BATCH(1)%AO(1)%p%type
 t2 = ODbat%BATCH(1)%AO(2)%p%type
 hermiteSingle = .FALSE.
 IF ( (t1.EQ.'Hermite'.AND.t2.EQ.'Empty').OR.(t1.EQ.'Empty'  .AND.t2.EQ.'Hermite') ) THEN
   hermiteSingle = .TRUE.
 ENDIF

Nullify(P%center)
Nullify(P%exponents)
Nullify(P%reducedExponents)
Nullify(P%preExpFac)
Nullify(P%iprim1)
Nullify(P%iprim2)
Allocate(P%center(3,maxprim))
Allocate(P%exponents(maxprim))
Allocate(P%reducedExponents(maxprim))
Allocate(P%preExpFac(maxprim))
Allocate(P%iprim1(maxprim))
Allocate(P%iprim2(maxprim))

NULLIFY(P%distance12)
ALLOCATE(P%distance12(3,1))

Nullify(P%angmom)
Allocate(P%angmom(maxnAngmom))
Nullify(P%indexAng1)
Allocate(P%indexAng1(maxnAngmom))
Nullify(P%indexAng2)
Allocate(P%indexAng2(maxnAngmom))
Nullify(P%nOrbitals)
Allocate(P%nOrbitals(maxnAngmom))
Nullify(P%nOrbComp)
Allocate(P%nOrbComp(maxnAngmom))
Nullify(P%nContracted)
Allocate(P%nContracted(maxnAngmom))

IF (INPUT%DO_JENGINE .AND. IELECTRON.EQ.2 ) THEN
   NULLIFY(P%FTUV)
   ALLOCATE(P%FTUV(maxPrim,maxTUV,Input%NDMAT_RHS))
   CALL QUIT('NOT DONE FOR JENGINE YET')
ELSE
   IF(.NOT. hermitesingle)THEN
      IF (INPUT%setETUVoutside) THEN
         NULLIFY(P%ETUVindex)
         ALLOCATE(P%ETUVindex(maxnAngmom))
         NULLIFY(P%ETUV)
         ALLOCATE(P%ETUV(MAXPRIMTUVIJK)) !dimensions
         P%ETUVisSet = .TRUE.
      ENDIF
   ENDIF
ENDIF
END SUBROUTINE INIT_OVERLAP

SUBROUTINE INIT_ORBITAL(Orb,maxprim,maxangmom,maxA)
IMPLICIT NONE
TYPE(Orbital) :: Orb
TYPE(IntegralItem)  :: Integral
Integer :: maxprim,maxangmom,maxA

 NULLIFY(Orb%exponents)
 ALLOCATE(Orb%exponents(maxprim))
 NULLIFY(Orb%angmom)
 ALLOCATE(Orb%angmom(maxangmom))
 NULLIFY(Orb%nContracted)
 ALLOCATE(Orb%nContracted(maxangmom))
 NULLIFY(Orb%startOrbital)
 ALLOCATE(Orb%startOrbital(maxangmom,1))
 NULLIFY(Orb%startprimOrbital)
 ALLOCATE(Orb%startprimOrbital(maxangmom,1))
 NULLIFY(Orb%nOrbComp)
 ALLOCATE(Orb%nOrbComp(maxangmom))
 NULLIFY(Orb%nPrimOrbComp)
 ALLOCATE(Orb%nPrimOrbComp(maxangmom))
 NULLIFY(Orb%nOrbitals)
 ALLOCATE(Orb%nOrbitals(maxangmom))
 NULLIFY(Orb%CC)
 ALLOCATE(Orb%CC(maxangmom))
 NULLIFY(Orb%SPH_MAT)
 ALLOCATE(Orb%SPH_MAT(maxA+1))

END SUBROUTINE INIT_ORBITAL

SUBROUTINE SET_INITIALIZED_OVERLAP(P,np,Input,sharedTUV,Integral,Alloc,ODB,IELECTRON,LUPRI,IPRINT,SIDE,nPrimAsInput,SETORBITAL)
Implicit none
TYPE(IntegralInput) :: Input
TYPE(IntegralItem)  :: Integral
TYPE(AllocItem)     :: Alloc
TYPE(Overlap)       :: P
TYPE(ODBATCH)       :: ODB
Character*(*)       :: side
Integer             :: IELECTRON,LUPRI,IPRINT,nderiv,np
TYPE(TUVitem)       :: SharedTUV
LOGICAL             :: nPrimAsInput,SETORBITAL
!
Integer             :: idir,i1,i2,i12,start1,start2,end1,end2,a,b,a1,a2,orb1,orb2
integer             :: l1,l2,ijk1,ijk2,ijk,maxijk,na
Real(realk)         :: e1,e2,d2,maxGab,maxPrimGab,maxPrimGabElm,signP
Character(len=80)   :: t1,t2
Real(realk),pointer :: GAB(:,:),primGAB(:,:) 
LOGICAL             :: LHS 
!
SELECT CASE(SIDE)
CASE('LHS')
   IF(INPUT%PS_SCREEN)THEN
      primGAB => INPUT%pGAB_LHS
      maxPrimGabElm=INPUT%PS_MAXELM_RHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
      GAB => INPUT%GAB_LHS
   ENDIF
   LHS=.TRUE.
   NDERIV=INPUT%NDERIVP
CASE('RHS')
   IF(INPUT%PS_SCREEN)THEN
      primGAB => INPUT%pGAB_RHS
      maxPrimGabElm=INPUT%PS_MAXELM_LHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
      GAB => INPUT%GAB_RHS
   ENDIF
   LHS=.FALSE.
   NDERIV=INPUT%NDERIVQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_INITIALIZED_OVERLAP =',SIDE
   CALL QUIT('Wrong case in SET_INITIALIZED_OVERLAP')
END SELECT

IF(SETORBITAL)THEN
 CALL SET_INITIALIZED_ORBITAL(P%orbital1,ODB%AO(1)%p,integral,LUPRI)
 CALL SET_INITIALIZED_ORBITAL(P%orbital2,ODB%AO(2)%p,integral,LUPRI)
ELSE
 !ORBITALS ALREADY SET BUT THE ORBITALS HAVE A NEW CENTER AND OTHER indices
 na   = ODB%AO(1)%p%nAngmom
 P%orbital1%center                    = ODB%AO(1)%p%center
 P%orbital1%startOrbital(1:na,1)      = ODB%AO(1)%p%startOrbital(1:na)
 P%orbital1%startprimOrbital(1:na,1)  = ODB%AO(1)%p%startprimOrbital(1:na)
 na   = ODB%AO(2)%p%nAngmom
 P%orbital2%center                    = ODB%AO(2)%p%center
 P%orbital2%startOrbital(1:na,1)      = ODB%AO(2)%p%startOrbital(1:na)
 P%orbital2%startprimOrbital(1:na,1)  = ODB%AO(2)%p%startprimOrbital(1:na)
ENDIF

!Determine overlap type
 t1 = P%orbital1%type
 t2 = P%orbital2%type
 P%hermiteSingle = .FALSE.
 IF ( t1.EQ.'Empty'.AND.t2.EQ.'Empty') THEN
   P%type = 'Empty'
 ELSE IF ( (t1.EQ.'Hermite'.AND.t2.EQ.'Empty').OR.(t1.EQ.'Empty'  .AND.t2.EQ.'Hermite') ) THEN
   P%type = 'Hermite-single'
   P%hermiteSingle = .TRUE.
 ELSE IF (t1.EQ.'Hermite'.AND.t2.EQ.'Hermite') THEN
   P%type = 'Hermite'
 ELSE IF ( (t1.EQ.'Cartesian'.AND.t2.EQ.'Empty').OR.(t1.EQ.'Empty'  .AND.t2.EQ.'Cartesian') ) THEN
   P%type = 'Cartesian-singel'
 ELSE IF (t1.EQ.'Cartesian'.AND.t2.EQ.'Cartesian') THEN
   P%type = 'Cartesian'
 ELSE
   WRITE(LUPRI,'(1X,A)') 'Not a proper combination of orbital types in SET_INITIALIZED_OVERLAP:'
   WRITE(LUPRI,'(5X,4A)') 'orbital1%type =', t1, 'and orbital2%type =', t2
   CALL QUIT('Not a proper combination of orbital types in SET_INITIALIZED_OVERLAP')
 ENDIF

 P%sameAO = ODB%sameAO
 P%ETUVisSet = .FALSE.
 P%sphericalETUV = .TRUE.
IF(.NOT.nPrimAsInput)THEN
 CALL GET_NPRIMITIVES(np,P,Input,Side)
ENDIF
 P%nPrimitives   = np
! Default is to build OD-batches first, and then later collect into passes
 P%nPasses       = 1
 IF (np.EQ.0) THEN
   RETURN
 ENDIF

 P%nAngmom       = ODB%nAngmom
 P%maxContracted = ODB%maxContracted

!
 d2 = 0.0d0
 DO idir=1,3
   P%distance12(idir,1) = P%orbital1%center(idir)-P%orbital2%center(idir)
   d2 = d2 + P%distance12(idir,1) * P%distance12(idir,1)
 ENDDO
 P%squaredDistance = d2
!
 !Nullify(P%center)
 !Nullify(P%exponents)
 !Nullify(P%reducedExponents)
 !Nullify(P%preExpFac)
 !Nullify(P%iprim1)
 !Nullify(P%iprim2)
 !Allocate(P%center(3,np))
 !Allocate(P%exponents(np))
 !Allocate(P%reducedExponents(np))
 !Allocate(P%preExpFac(np))
 !Allocate(P%iprim1(np))
 !Allocate(P%iprim2(np))
  
 i12 = 0
 DO i1=1,P%orbital1%nPrimitives
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nPrimitives
     maxPrimGab=0.0d0
     IF(INPUT%PS_SCREEN)THEN
        DO a1=1,P%orbital1%nAngmom 
           DO a2=1,P%orbital2%nAngmom
              a = P%orbital1%startprimOrbital(a1,1)+(i1-1)*P%orbital1%nOrbComp(a1)
              b = P%orbital2%startprimOrbital(a2,1)+(i2-1)*P%orbital2%nOrbComp(a2)
              maxPrimGab=MAX(maxPrimGab,primGAB(a,b))
           ENDDO
        ENDDO
     ENDIF
     IF(INPUT%PS_SCREEN)THEN
        IF(maxPrimGab .GT. INPUT%PS_THRESHOLD/maxPrimGabElm)THEN
         i12 = i12 + 1
         e1  = P%orbital1%exponents(i1)
         e2  = P%orbital2%exponents(i2)
         P%exponents(i12)        = e1 + e2
         IF (P%exponents(i12) .NE. 0.0d0) THEN
            P%reducedExponents(i12) = e1*e2/(e1+e2)
            P%center(1,i12)         = (e1*P%orbital1%center(1) + e2*P%orbital2%center(1))/(e1+e2)
            P%center(2,i12)         = (e1*P%orbital1%center(2) + e2*P%orbital2%center(2))/(e1+e2)
            P%center(3,i12)         = (e1*P%orbital1%center(3) + e2*P%orbital2%center(3))/(e1+e2)
            P%preExpFac(i12)        = exp(-P%reducedExponents(i12)*d2)
         ELSEIF(Input%operator(1:6) .EQ. 'Nucrep')THEN
            P%reducedExponents(i12) = e1
            P%center(1,i12)         = P%orbital1%center(1)
            P%center(2,i12)         = P%orbital1%center(2)
            P%center(3,i12)         = P%orbital1%center(3)
            P%preExpFac(i12)        = exp(-e1*d2)
         ELSE
           P%reducedExponents(i12) = 0.0d0
           P%center(1,i12)         = 0.0d0
           P%center(2,i12)         = 0.0d0
           P%center(3,i12)         = 0.0d0
           P%preExpFac(i12)        = 1.0d0
         ENDIF
         IF ( t1.EQ.'Empty'.OR.t2.EQ.'Empty') THEN 
            P%preExpFac(i12)        = 1.0d0
         ELSE
            IF (P%exponents(i12) .NE. 0.0d0) THEN
               P%preExpFac(i12)        = exp(-P%reducedExponents(i12)*d2)
            ELSEIF(Input%operator(1:6) .EQ. 'Nucrep')THEN
               P%preExpFac(i12)        = exp(-e1*d2)
            ELSE
               P%preExpFac(i12)        = 1.0d0
            ENDIF
         ENDIF
         P%iprim1(i12)           = i1
         P%iprim2(i12)           = i2
        ENDIF 
     ELSE
         i12 = i12 + 1
         e1  = P%orbital1%exponents(i1)
         e2  = P%orbital2%exponents(i2)
         P%exponents(i12)        = e1 + e2
         IF (P%exponents(i12) .NE. 0.0d0) THEN
            P%reducedExponents(i12) = e1*e2/(e1+e2)
            P%center(1,i12)         = (e1*P%orbital1%center(1) + e2*P%orbital2%center(1))/(e1+e2)
            P%center(2,i12)         = (e1*P%orbital1%center(2) + e2*P%orbital2%center(2))/(e1+e2)
            P%center(3,i12)         = (e1*P%orbital1%center(3) + e2*P%orbital2%center(3))/(e1+e2)
            P%preExpFac(i12)        = exp(-P%reducedExponents(i12)*d2)
         ELSEIF(Input%operator(1:6) .EQ. 'Nucrep')THEN
            P%reducedExponents(i12) = e1
            P%center(1,i12)         = P%orbital1%center(1)
            P%center(2,i12)         = P%orbital1%center(2)
            P%center(3,i12)         = P%orbital1%center(3)
            P%preExpFac(i12)        = exp(-e1*d2)
         ELSE
           P%reducedExponents(i12) = 0.0d0
           P%center(1,i12)         = 0.0d0
           P%center(2,i12)         = 0.0d0
           P%center(3,i12)         = 0.0d0
           P%preExpFac(i12)        = 1.0d0
         ENDIF
         IF ( t1.EQ.'Empty'.OR.t2.EQ.'Empty') THEN 
            P%preExpFac(i12)        = 1.0d0
         ELSE
            IF (P%exponents(i12) .NE. 0.0d0) THEN
               P%preExpFac(i12)        = exp(-P%reducedExponents(i12)*d2)
            ELSEIF(Input%operator(1:6) .EQ. 'Nucrep')THEN
               P%preExpFac(i12)        = exp(-e1*d2)
            ELSE
               P%preExpFac(i12)        = 1.0d0
            ENDIF
         ENDIF
         P%iprim1(i12)           = i1
         P%iprim2(i12)           = i2
     ENDIF
   ENDDO
 ENDDO
 IF (np .NE. i12) THEN
   WRITE(LUPRI,'(1X,A,2I5)') 'Error in set_overlap. Mismatch in number of primitives, (np,i12)=',&
     &  np,i12
   CALL QUIT('Error: Mismatch between get_nPrimitives and i12 in set_INITIALIZED_overlap')
 ENDIF

IF(INPUT%PS_SCREEN)THEN
   P%maxPrimGab = maxPrimGab
   NULLIFY(primGAB)
ENDIF

 !Nullify(P%angmom)
 !Allocate(P%angmom(P%nAngmom))
 !Nullify(P%indexAng1)
 !Allocate(P%indexAng1(P%nAngmom))
 !Nullify(P%indexAng2)
 !Allocate(P%indexAng2(P%nAngmom))
 !Nullify(P%nOrbitals)
 !Allocate(P%nOrbitals(P%nAngmom))
 !Nullify(P%nOrbComp)
 !Allocate(P%nOrbComp(P%nAngmom))
 !Nullify(P%nContracted)
 !Allocate(P%nContracted(P%nAngmom))
!CAREFUL
!Initialize to twice some maximal value + 1
 P%minAngmom   = 99
!CAREFUL
 P%maxAngmom = 0
 P%totOrbitals = 0
 i12 = 0
 maxGab = 0.0d0
 maxijk = 0
 DO i1=1,P%orbital1%nAngmom
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nAngmom
     i12  = i12 + 1
     l1   = P%orbital1%angmom(i1)
     l2   = P%orbital2%angmom(i2)
     ijk1 = (l1+1)*(l1+2)/2
     ijk2 = (l2+1)*(l2+2)/2
     ijk  = ijk1*ijk2
     maxijk = max(maxijk,ijk)
     P%angmom(i12)      = P%orbital1%angmom(i1) + P%orbital2%angmom(i2)
     P%indexAng1(i12)   = i1
     P%indexAng2(i12)   = i2
     P%nOrbitals(i12)   = P%orbital1%nOrbitals(i1)*P%orbital2%nOrbitals(i2)
     P%totOrbitals      = P%totOrbitals + P%nOrbitals(i12)*NDERIV
     P%nOrbComp(i12)    = P%orbital1%nOrbComp(i1)*P%orbital2%nOrbComp(i2)
     P%nContracted(i12) = P%orbital1%nContracted(i1)*P%orbital2%nContracted(i2)
     P%minAngmom        = min(P%minAngmom,P%angmom(i12))
     P%maxAngmom        = max(P%maxAngmom,P%angmom(i12))
     start1 = P%orbital1%startOrbital(i1,1)
     end1   = start1 + P%orbital1%nOrbitals(i1) - 1
     start2 = P%orbital2%startOrbital(i2,1)
     end2   = start2 + P%orbital2%nOrbitals(i2) - 1
     IF (Input%CS_SCREEN) THEN
       DO orb1=start1,end1
         DO orb2=start2,end2
           maxGab = MAX(maxGab,GAB(orb1,orb2))
         ENDDO
       ENDDO
     ENDIF
   ENDDO
 ENDDO
 IF (Input%CS_SCREEN) THEN
   P%maxGab = maxGab
   NULLIFY(GAB)
 ENDIF

 P%startAngmom = 0
 IF (P%hermiteSingle) P%startAngmom = P%minAngmom

 IF( INPUT%operator(1:7) .EQ. 'Kinetic' .AND. LHS) THEN
    P%startAngmom = P%startAngmom + 2
    P%endAngmom   = P%maxAngmom + 2
 ELSE
    P%endAngmom   = P%maxAngmom
 ENDIF  

 P%nTUV = (P%endAngmom+1)*(P%endAngmom+2)*(P%endAngmom+3)/6 - &
    &     P%startAngmom*(P%startAngmom+1)*(P%startAngmom+2)/6

 IF (IELECTRON.EQ.1) THEN
   signP=1.0d0
!   Alloc%maxPrimLHS    = max(np,Alloc%maxPrimLHS)
!   Alloc%maxPrimTUVLHS = max(np*P%nTUV,np*P%totOrbitals,np*maxijk,Alloc%maxPrimTUVLHS)
!   Alloc%maxAngLHS     = max(P%maxAngmom,Alloc%maxAngLHS)
 ELSE
   signP=-1.0d0
!   Alloc%maxPrimRHS    = max(np,Alloc%maxPrimRHS)
!   Alloc%maxPrimTUVRHS = max(np*P%nTUV,np*P%totOrbitals,np*maxijk,Alloc%maxPrimTUVRHS)
!   Alloc%maxAngRHS     = max(P%maxAngmom,Alloc%maxAngRHS)
 ENDIF

IF (INPUT%DO_JENGINE .AND. IELECTRON.EQ.2) THEN
  CALL SET_INITIALIZED_FTUV(P,INPUT,Integral,LUPRI,IPRINT)
! Change Overlap P to confirm with FTUV
  P%TYPE = 'FTUV'
  P%nAngmom = 1
  P%orbital1%nAngmom = 1
  P%orbital2%nAngmom = 1
  P%angmom(1) = P%maxAngmom
  P%orbital1%angmom(1) = 0
  P%orbital2%angmom(1) = 0
  P%orbital1%nContracted(1)  = 1
  P%orbital2%nContracted(1)  = 1
  P%orbital1%nOrbComp(1)     = 1
  P%orbital2%nOrbComp(1)     = 1
  P%orbital1%startOrbital(1,1) = 1
  P%orbital2%startOrbital(1,1) = 1
  P%indexAng1(1)   = 1
  P%indexAng2(1)   = 1
  P%nOrbitals(1)   = 1
  P%totOrbitals    = 1
  P%nOrbComp(1)    = 1
  P%nContracted(1) = 1
ELSE
!GET E COEFFICIENTS
 IF (.NOT.P%hermiteSingle) THEN
   IF (INPUT%setETUVoutside) THEN
    CALL AttachETUVtoInitializedOverlap(P,signP,SharedTUV,LUPRI,IPRINT)
   ENDIF
 ENDIF
ENDIF

IF(INPUT%CS_SCREEN)THEN
   NULLIFY(GAB)
ENDIF
END SUBROUTINE SET_INITIALIZED_OVERLAP

SUBROUTINE SET_INITIALIZED_ORBITAL(Orb,AOB,integral,LUPRI)
IMPLICIT NONE
TYPE(Orbital) :: Orb
TYPE(AOBATCH) :: AOB
TYPE(IntegralItem)  :: Integral
integer             :: lupri
!
Integer :: ia,na,np,maxc,dim1,dim2,i,L,nc
!
 np   = AOB%nPrimitives
 na   = AOB%nAngmom

 Orb%type          = AOB%type
 Orb%spherical     = AOB%spherical
 Orb%nPrimitives   = np 
 Orb%nPasses       = 1 
 Orb%maxContracted = AOB%maxContracted
 Orb%maxAngmom     = AOB%maxAngmom
 Orb%nAngmom       = na
 Orb%CCidentifier  = AOB%CCindex(1)
 Orb%center        = AOB%center
 Orb%exponents(1:np)     = AOB%pExponents%elms(1:np)
 Orb%angmom(1:na)     = AOB%angmom(1:na)
 Orb%nContracted(1:na)   = AOB%nContracted(1:na)
 Orb%startOrbital(1:na,1)  = AOB%startOrbital(1:na)
 Orb%startprimOrbital(1:na,1)  = AOB%startprimOrbital(1:na)
 Orb%nOrbComp(1:na)  = AOB%nOrbComp(1:na)
 Orb%nPrimOrbComp(1:na)  = AOB%nPrimOrbComp(1:na)
 Orb%nOrbitals(1:na)  = AOB%nOrbitals(1:na)
 DO ia=1,na
    Orb%CC(ia)%p => AOB%pCC(ia)%p
 ENDDO
 DO I=1,Orb%maxAngmom+1
    Orb%SPH_MAT(I)%p => integral%TUV%SPH_MAT(I)
 ENDDO

END SUBROUTINE SET_INITIALIZED_ORBITAL

SUBROUTINE AttachETUVtoInitializedOverlap(P,signP,TUV,LUPRI,IPRINT)
Implicit none
TYPE(Overlap) :: P
TYPE(TUVitem) :: TUV
Integer       :: LUPRI,IPRINT
Real(realk)   :: signP
!
Integer :: iA1,iA2,l1,l2,l,ijk1,ijk2,ijk,nTUV,indexETUV,nETUV,iAngmom
integer :: LENGTH

!NULLIFY(P%ETUVindex)
!ALLOCATE(P%ETUVindex(P%nAngmom))
indexETUV = 1
DO iAngmom=1,P%nAngmom
  P%ETUVindex(iAngmom) = indexETUV
  iA1 = P%indexAng1(iangmom)
  iA2 = P%indexAng2(iangmom)
  l1 = P%orbital1%angmom(iA1)
  l2 = P%orbital2%angmom(iA2)
  l  = l1+l2
  ijk1 = (l1+1)*(l1+2)/2
  ijk2 = (l2+1)*(l2+2)/2
  ijk  = ijk1*ijk2
  nTUV = (l+1)*(l+2)*(l+3)/6
  indexETUV = indexETUV + nTUV*ijk*P%nPrimitives
ENDDO
nETUV = indexETUV - 1
P%lenETUV = nETUV

!NULLIFY(P%ETUV)
!ALLOCATE(P%ETUV(nETUV))
DO iAngmom=1,P%nAngmom
  indexETUV = P%ETUVindex(iAngmom)
  iA1 = P%indexAng1(iangmom)
  iA2 = P%indexAng2(iangmom)
  l1 = P%orbital1%angmom(iA1)
  l2 = P%orbital2%angmom(iA2)
  l  = l1+l2
  ijk1 = (l1+1)*(l1+2)/2
  ijk2 = (l2+1)*(l2+2)/2
  ijk  = ijk1*ijk2
  nTUV = (l+1)*(l+2)*(l+3)/6
  LENGTH=ijk*nTUV*P%nPrimitives
  CALL BuildEcoeffTensor2(TUV,P,signP,P%ETUV(indexETUV:indexETUV+LENGTH-1),ijk,nTUV,P%nPrimitives,iAngmom,1,LUPRI,IPRINT)
ENDDO

P%ETUVisSet = .TRUE.

END SUBROUTINE AttachETUVtoInitializedOverlap

SUBROUTINE SET_INITIALIZED_FTUV(P,Input,Integral,LUPRI,IPRINT)
Implicit none
TYPE(Overlap)       :: P
TYPE(IntegralInput) :: Input
TYPE(IntegralItem)  :: Integral
Integer             :: LUPRI,IPRINT
!
Logical             :: Sph1,Sph2
Integer             :: iAngmom,l1,l2,l,ijk1,ijk2,ijk,ioff,nTUV,lm1,lm2,lm
Integer             :: iS1,iS2,iPrim,iA1,iA2,nCont,nC1,nC2,dim1,iAng1,iAng2
Integer             :: idmat,start1,start2,i1,i2,iC1,iC2,indlm,iTUV
Real(realk),pointer :: Ecoeffs(:,:,:),Spherical(:,:),SpherEcoeff(:,:,:)
Real(realk),pointer :: CC(:,:,:),ContractedDmat(:,:,:)
Real(realk)         :: DM1=-1.0d0,D1=1.0d0,D0=0.0d0,dtemp
Integer             :: eS1,eS2

!NULLIFY(P%FTUV)
!ALLOCATE(P%FTUV(P%nPrimitives,P%nTUV,Input%NDMAT_RHS))
CALL DZERO(P%FTUV,P%nTUV*P%nPrimitives*Input%NDMAT_RHS)

DO iAngmom = 1,P%nAngmom
  IF (P%hermiteSingle) THEN
    CALL QUIT('ToDo Hermite-single SET_FTUV')
  ELSE
    iA1 = P%indexAng1(iangmom)
    iA2 = P%indexAng2(iangmom)
    l1 = P%orbital1%angmom(iA1)
    l2 = P%orbital2%angmom(iA2)
    l  = l1+l2
    ijk1 = (l1+1)*(l1+2)/2
    ijk2 = (l2+1)*(l2+2)/2
    ijk  = ijk1*ijk2
    nTUV = (l+1)*(l+2)*(l+3)/6
    NULLIFY(Ecoeffs)
    ALLOCATE(Ecoeffs(nTUV,ijk,P%nPrimitives))

    CALL BuildEcoeffTensor1(Integral%TUV,P,DM1,Ecoeffs,nTUV,ijk,P%nPrimitives,iAngmom,1,LUPRI,IPRINT)

!   Spherical transformation of E-coefficients
    Sph1 = P%orbital1%spherical.AND.(l1.GT.1)
    Sph2 = P%orbital2%spherical.AND.(l2.GT.1)
    lm1 = ijk1
    IF (Sph1) lm1 = 2*l1+1
    lm2 = ijk2
    IF (Sph2) lm2 = 2*l2+1
    lm  = lm1*lm2
    IF (Sph1.OR.Sph2) THEN
      NULLIFY(Spherical)
      ALLOCATE(Spherical(ijk,lm))

!     CALL ContructSphericalTransformation(Spherical,Reshape(P%orbital1%SPHMAT(iS1:eS1),(/lm1,ijk1/)),&
!    &                                     Reshape(P%orbital2%SPHMAT(iS2:eS2),(/lm2,ijk2/)),&
!    &                                     lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)
      CALL ContructSphericalTransformation(Spherical,P%orbital1%SPH_MAT(l1+1)%p%elms,&
     &                                     P%orbital2%SPH_MAT(l2+1)%p%elms,&
     &                                     lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)
      NULLIFY(SpherEcoeff)
      ALLOCATE(SpherEcoeff(nTUV,lm,P%nPrimitives+ijk*lm+nTUV*ijk*P%nPrimitives))
!      DO iPrim=1,P%nPrimitives
      DO iPrim=1,P%nPrimitives


        CALL DGEMM('N','N',nTUV,lm,ijk,1.0d0,Ecoeffs(1,1,iPrim),nTUV,Spherical,ijk,&
     &             0.0d0,SpherEcoeff(1,1,iPrim),nTUV)
      ENDDO
      DEALLOCATE(Ecoeffs)
      DEALLOCATE(Spherical)
    ELSE
      NULLIFY(SpherEcoeff)
      SpherEcoeff => Ecoeffs

    ENDIF
!   Set up contraction matrix
    iA1 = P%indexAng1(iangmom)
    iA2 = P%indexAng2(iangmom)
    nC1 = P%orbital1%nContracted(iA1)
    nC2 = P%orbital2%nContracted(iA2)
    nCont = nC1*nC2
    NULLIFY(CC)
    ALLOCATE(CC(P%nPrimitives,nC1,nC2))
    CALL ConstructContraction(CC,P,1,P%nPrimitives,nC1,nC2,iA1,iA2,LUPRI,IPRINT)

!   Make contraction of density-matrix and contraction coefficients
    NULLIFY(ContractedDmat)
    ALLOCATE(ContractedDmat(P%nPrimitives,lm,Input%NDMAT_RHS))
    CALL DZERO(ContractedDmat,lm*P%nPrimitives*Input%NDMAT_RHS)
    start1 = P%orbital1%startOrbital(iA1,1)
    start2 = P%orbital2%startOrbital(iA2,1)
    DO idmat=1,Input%nDMAT_RHS
      i1 = start1
      DO iC1=1,nC1
        DO iAng1=1,lm1
          i2 = start2
          DO iC2=1,nC2
            DO iAng2=1,lm2
              indlm = (iAng2-1)*lm1 + iAng1
              dtemp = Input%DMAT_RHS(i1,i2,idmat)
              IF ((start1.NE.start2).AND.Input%sameRHSaos) &
     &            dtemp = dtemp + Input%DMAT_RHS(i2,i1,idmat)
              dtemp = dtemp*Input%CoulombFactor
              DO iPrim=1,P%nPrimitives
                ContractedDmat(iPrim,indlm,idmat) = ContractedDmat(iPrim,indlm,idmat) &
     &                                              + dtemp * CC(iPrim,iC1,iC2)
              ENDDO
              i2=i2+1
            ENDDO
          ENDDO
          i1=i1+1
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(CC)

!   Set up FTUV's
    DO idmat=1,Input%NDMAT_RHS
      DO iTUV=1,nTUV
        DO indlm=1,lm
          DO iPrim=1,P%nPrimitives
            P%FTUV(iPrim,iTUV,idmat) = P%FTUV(iPrim,iTUV,idmat) &
     &               + ContractedDmat(iPrim,indlm,idmat)*SpherEcoeff(iTUV,indlm,iPrim)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(ContractedDmat)
    IF (Sph1.OR.Sph2) THEN
      DEALLOCATE(SpherEcoeff)
    ELSE
      NULLIFY(SpherEcoeff)
      DEALLOCATE(Ecoeffs)
    ENDIF
  ENDIF
ENDDO
END SUBROUTINE SET_INITIALIZED_FTUV

SUBROUTINE INIT_PASS(PassP,Alloc,Input,ODbat,SIDE,IELECTRON,maxpasses,IPRINT,LUPRI)
Implicit none
INTEGER             :: nbatches,lupri,IPRINT,IELECTRON
TYPE(Overlap)       :: PassP
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
TYPE(ODITEM)        :: ODbat
!TYPE(Integralitem)  :: Integral
Character*(*)       :: side
!TYPE(ODBATCH)       :: ODbat(nbatches)
!
Character(len=80)   :: t1,t2
Integer             :: i1,i2,l1,l2,ijk1,ijk2,ijk,maxijk,I,nprim
Integer             :: maxangmom,totOrbitals,maxtotorbitals,NDERIV
Integer             :: maxprim,endangmom,maxtuv,maxcont,maxnangmom
Integer             :: MAXPRIMTUV,MAXPRIMTOTORB,MAXPRIMIJK,MAXnTUV
Integer             :: MAXPRIMTUVIJK,maxA,maxB,maxpasses
Integer             :: lenETUV,ma,mc,mp,na,np
LOGICAL             :: LHS,hermiteSingle 

SELECT CASE(SIDE)
CASE('LHS')
   LHS=.TRUE.
!ORBITAL
   np = Alloc%maxPrimLHS 
   na = Alloc%maxnAngLHS 
   mc = Alloc%maxContLHS
   ma = Alloc%maxAngmomLHS
CASE('RHS')
   LHS=.FALSE.
   np = Alloc%maxPrimRHS 
   na = Alloc%maxnAngRHS 
   mc = Alloc%maxContRHS
   ma = Alloc%maxAngmomRHS
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_INITIALIZED_OVERLAP =',SIDE
   CALL QUIT('Wrong case in SET_INITIALIZED_OVERLAP')
END SELECT

   NULLIFY(PassP%Orbital1%exponents)
   ALLOCATE(PassP%Orbital1%exponents(np))
   NULLIFY(PassP%Orbital1%angmom)
   ALLOCATE(PassP%Orbital1%angmom(na))
   NULLIFY(PassP%Orbital1%nContracted)
   ALLOCATE(PassP%Orbital1%nContracted(na))
   NULLIFY(PassP%Orbital1%startOrbital)
   ALLOCATE(PassP%Orbital1%startOrbital(na,maxpasses))
   NULLIFY(PassP%Orbital1%startprimOrbital)
   ALLOCATE(PassP%Orbital1%startprimOrbital(na,maxpasses))
   NULLIFY(PassP%Orbital1%nOrbComp)
   ALLOCATE(PassP%Orbital1%nOrbComp(na))
   NULLIFY(PassP%Orbital1%nPrimOrbComp)
   ALLOCATE(PassP%Orbital1%nPrimOrbComp(na))
   NULLIFY(PassP%Orbital1%nOrbitals)
   ALLOCATE(PassP%Orbital1%nOrbitals(na))
!   NULLIFY(PassP%Orbital1%CC)
!   ALLOCATE(PassP%Orbital1%CC(np,mc,na))
 
  NULLIFY(PassP%orbital1%CC)
  ALLOCATE(PassP%orbital1%CC(na))
 
   NULLIFY(PassP%Orbital1%SPH_MAT)
   ALLOCATE(PassP%Orbital1%SPH_MAT(ma+1))

!IF(LHS)THEN !ORBITAL 2 
!   np = Alloc%maxPrimLHS 
!   na = Alloc%maxnAngLHS 
!   mc = Alloc%maxContLHS
!   ma = Alloc%maxAngmomLHS
!ELSE
!   np = Alloc%maxPrimRHS 
!   na = Alloc%maxnAngRHS 
!   mc = Alloc%maxContRHS
!   ma = Alloc%maxAngmomRHS
!ENDIF
  
   NULLIFY(PassP%Orbital2%exponents)
   ALLOCATE(PassP%Orbital2%exponents(np))
   NULLIFY(PassP%Orbital2%angmom)
   ALLOCATE(PassP%Orbital2%angmom(na))
   NULLIFY(PassP%Orbital2%nContracted)
   ALLOCATE(PassP%Orbital2%nContracted(na))
   NULLIFY(PassP%Orbital2%startOrbital)
   ALLOCATE(PassP%Orbital2%startOrbital(na,maxpasses))
   NULLIFY(PassP%Orbital2%startprimOrbital)
   ALLOCATE(PassP%Orbital2%startprimOrbital(na,maxpasses))
   NULLIFY(PassP%Orbital2%nOrbComp)
   ALLOCATE(PassP%Orbital2%nOrbComp(na))
   NULLIFY(PassP%Orbital2%nPrimOrbComp)
   ALLOCATE(PassP%Orbital2%nPrimOrbComp(na))
   NULLIFY(PassP%Orbital2%nOrbitals)
   ALLOCATE(PassP%Orbital2%nOrbitals(na))
!   NULLIFY(PassP%Orbital2%CC)
!   ALLOCATE(PassP%Orbital2%CC(np,mc,na))
  NULLIFY(PassP%orbital2%CC)
  ALLOCATE(PassP%orbital2%CC(na))

   NULLIFY(PassP%Orbital2%SPH_MAT)
   ALLOCATE(PassP%Orbital2%SPH_MAT(ma+1))

IF(LHS)THEN !OVERLAP
   mp = Alloc%maxPrimLHS*maxPasses
   np = Alloc%maxPrimLHS 
   na = Alloc%maxnAngLHS 
   lenETUV = Alloc%maxPrimTUVijkLHS
ELSE
   mp = Alloc%maxPrimRHS*maxPasses
   np = Alloc%maxPrimRHS 
   na = Alloc%maxnAngRHS 
   lenETUV = Alloc%maxPrimTUVijkRHS
ENDIF

   Nullify(PassP%center)
   Nullify(PassP%exponents)
   Nullify(PassP%reducedExponents)
   Nullify(PassP%preExpFac)
   Nullify(PassP%iprim1)
   Nullify(PassP%iprim2)

   Allocate(PassP%center(3,mp))
   Allocate(PassP%exponents(mp))
   Allocate(PassP%reducedExponents(mp))
   Allocate(PassP%preExpFac(mp))
   Allocate(PassP%iprim1(mp))
   Allocate(PassP%iprim2(mp))

   Nullify(PassP%Angmom)
   Nullify(PassP%indexAng1)
   Nullify(PassP%indexAng2)
   Nullify(PassP%nOrbComp)
   Nullify(PassP%nContracted)
   Nullify(PassP%nOrbitals)
   Nullify(PassP%ETUVindex)

   Allocate(PassP%Angmom(na))
   Allocate(PassP%indexAng1(na))
   Allocate(PassP%indexAng2(na))
   Allocate(PassP%nOrbComp(na))
   Allocate(PassP%nContracted(na))
   Allocate(PassP%nOrbitals(na))
   Allocate(PassP%ETUVindex(na))

   t1 = ODbat%BATCH(1)%AO(1)%p%type
   t2 = ODbat%BATCH(1)%AO(2)%p%type
   hermiteSingle = .FALSE.
   IF ( (t1.EQ.'Hermite'.AND.t2.EQ.'Empty').OR.(t1.EQ.'Empty'  .AND.t2.EQ.'Hermite') ) THEN
      hermiteSingle = .TRUE.
   ENDIF
   NULLIFY(PassP%distance12)
   ALLOCATE(PassP%distance12(3,maxpasses))
   IF(.NOT. hermitesingle)THEN
      IF (INPUT%setETUVoutside) THEN
         NULLIFY(PassP%ETUV)
         ALLOCATE(PassP%ETUV(lenETUV*maxPasses))
      ELSE
!         NULLIFY(PassP%distance12)
!         ALLOCATE(PassP%distance12(3,maxpasses))
      ENDIF
   ENDIF

 END SUBROUTINE INIT_PASS

SUBROUTINE FreePass(PassP)
implicit none
TYPE(Overlap) :: PassP

   DEALLOCATE(PassP%Orbital1%exponents)
   DEALLOCATE(PassP%Orbital1%angmom)
   DEALLOCATE(PassP%Orbital1%nContracted)
   DEALLOCATE(PassP%Orbital1%startOrbital)
   DEALLOCATE(PassP%Orbital1%startprimOrbital)
   DEALLOCATE(PassP%Orbital1%nOrbComp)
   DEALLOCATE(PassP%Orbital1%nPrimOrbComp)
   DEALLOCATE(PassP%Orbital1%nOrbitals)
   DEALLOCATE(PassP%Orbital1%CC)
   DEALLOCATE(PassP%Orbital1%SPH_MAT)

   NULLIFY(PassP%Orbital1%exponents)
   NULLIFY(PassP%Orbital1%angmom)
   NULLIFY(PassP%Orbital1%nContracted)
   NULLIFY(PassP%Orbital1%startOrbital)
   NULLIFY(PassP%Orbital1%startprimOrbital)
   NULLIFY(PassP%Orbital1%nOrbComp)
   NULLIFY(PassP%Orbital1%nPrimOrbComp)
   NULLIFY(PassP%Orbital1%nOrbitals)
   NULLIFY(PassP%Orbital1%CC)
   NULLIFY(PassP%Orbital1%SPH_MAT)

   DEALLOCATE(PassP%Orbital2%exponents)
   DEALLOCATE(PassP%Orbital2%angmom)
   DEALLOCATE(PassP%Orbital2%nContracted)
   DEALLOCATE(PassP%Orbital2%startOrbital)
   DEALLOCATE(PassP%Orbital2%startprimOrbital)
   DEALLOCATE(PassP%Orbital2%nOrbComp)
   DEALLOCATE(PassP%Orbital2%nPrimOrbComp)
   DEALLOCATE(PassP%Orbital2%nOrbitals)
   DEALLOCATE(PassP%Orbital2%CC)
   DEALLOCATE(PassP%Orbital2%SPH_MAT)

   NULLIFY(PassP%Orbital2%exponents)
   NULLIFY(PassP%Orbital2%angmom)
   NULLIFY(PassP%Orbital2%nContracted)
   NULLIFY(PassP%Orbital2%startOrbital)
   NULLIFY(PassP%Orbital2%startprimOrbital)
   NULLIFY(PassP%Orbital2%nOrbComp)
   NULLIFY(PassP%Orbital2%nPrimOrbComp)
   NULLIFY(PassP%Orbital2%nOrbitals)
   NULLIFY(PassP%Orbital2%CC)
   NULLIFY(PassP%Orbital2%SPH_MAT)

   DeAllocate(PassP%center)
   DeAllocate(PassP%exponents)
   DeAllocate(PassP%reducedExponents)
   DeAllocate(PassP%preExpFac)
   DeAllocate(PassP%iprim1)
   DeAllocate(PassP%iprim2)
   DeAllocate(PassP%Angmom)
   DeAllocate(PassP%indexAng1)
   DeAllocate(PassP%indexAng2)
   DeAllocate(PassP%nOrbComp)
   DeAllocate(PassP%nContracted)
   DeAllocate(PassP%nOrbitals)
   DeAllocate(PassP%ETUVindex)

   Nullify(PassP%center)
   Nullify(PassP%exponents)
   Nullify(PassP%reducedExponents)
   Nullify(PassP%preExpFac)
   Nullify(PassP%iprim1)
   Nullify(PassP%iprim2)
   Nullify(PassP%nOrbitals)
   Nullify(PassP%Angmom)
   Nullify(PassP%indexAng1)
   Nullify(PassP%indexAng2)
   Nullify(PassP%nOrbComp)
   Nullify(PassP%nContracted)
   Nullify(PassP%nOrbitals)
   Nullify(PassP%ETUVindex)

   DEALLOCATE(PassP%distance12)
   NULLIFY(PassP%distance12)

   IF (PassP%ETUVisSet) THEN
      DEALLOCATE(PassP%ETUV)
      NULLIFY(PassP%ETUV)
   ENDIF

END SUBROUTINE FreePass

END MODULE integraldriver
