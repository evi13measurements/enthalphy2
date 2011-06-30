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
!****** ORBITAL
TYPE Orbital
Character(len=80)       :: type
Logical                 :: spherical
Integer                 :: maxAngmom
Integer                 :: nAngmom
Integer                 :: nPrimitives
Integer                 :: maxContracted
Real(realk)             :: center(3)
!Dimensions according to nAngmom
integer,pointer     :: angmom(:) !1:nangmom
integer,pointer     :: nContracted(:)
integer,pointer     :: startOrbital(:)
integer,pointer     :: startprimOrbital(:)
integer,pointer     :: nOrbComp(:)
integer,pointer     :: nPrimOrbComp(:)
integer,pointer     :: nOrbitals(:)
!Dimensions according to nPrimitives
real(realk),pointer :: exponents(:)
!Dimensions according to nPrimitives,maxContracted,nAngmom
real(realk),pointer :: CC(:,:,:)
TYPE(SPHMATPOINTER),pointer   :: SPH_MAT(:)
END TYPE Orbital


!****** OVERLAPPOINTER
TYPE Overlappointer
TYPE(Overlap),pointer :: p
END TYPE Overlappointer
!****** OVERLAP
TYPE Overlap
Character(len=80)   :: type
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
Real(realk)         :: distance12(3)
Real(realk)         :: squaredDistance
!Dimensions according to nPrimitives
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

!McMurcie-Davidson E-coefficients (and/or F-coefficients in case of J-engine) 
Logical             :: ETUVisSet
Logical             :: sphericalETUV
Real(realk),pointer :: ETUV(:) ! ntuv, ijk, nprim, iAngmom
Real(realk),pointer :: FTUV(:,:,:)
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
LOGICAL              :: sameOD
END TYPE Integrand

!****** INTEGRAL
TYPE Integralitem
Real(realk),pointer    :: Wtuv(:,:)       !nPrimitives:ntuv
Real(realk),pointer    :: tuvTUV(:,:,:,:) !ntuvQ(i):ntuvP:nPrimP:nPrimQ
Real(realk),pointer    :: TUVQ(:,:,:)     !ntuvP:nPrimP:nOrbQ
Real(realk),pointer    :: PQ(:,:)         !nOrbQ:nOrbP
!Generic integral intermadiates used of contraction with Ecoefficients, 
!Spherical transformation, and contraction to contracted basis
Real(realk),pointer    :: IntegralIN(:,:,:)
Real(realk),pointer    :: IntegralOUT(:,:,:)
INTEGER                :: nTUV
INTEGER                :: nOrb
INTEGER                :: nAng
INTEGER                :: nPrim
LOGICAL                :: Jengine
TYPE(TUVitem),pointer  :: TUV
END TYPE Integralitem

TYPE TUVitem
INTEGER                :: nTABFJW
Real(realk),pointer    :: TABFJW(:)
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
Character(len=80)    :: IDENTIFIER
!
Real(realk)          :: TSTART,TEND
!
#if MOD_FMM
!Simen temporary to interface FMM
Integer,parameter :: LWRK = 1
Real(realk) :: WRK(LWRK)
!Simen
#endif

CALL LSTIMER('START ',TSTART,TEND,LUPRI)
!
IDENTIFIER='ANDY'

!IPRINT = 51

CALL LSHEADER(LUPRI,'THE MAIN_INTEGRAL DRIVER')

IF(IPRINT .GT. 5) WRITE(LUPRI,*)'CALLING CREATE OD'
CALL Create_ODbatches(OD_LHS,IDENTIFIER,INPUT,'LHS',LUPRI)
CALL Create_ODbatches(OD_RHS,IDENTIFIER,INPUT,'RHS',LUPRI)

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
   CALL LINK(OD_LHS,OD_RHS,INPUT,OUTPUT,LUPRI,IPRINT)
ELSE
  CALL McMurchieDavidson(OD_LHS,OD_RHS,INPUT,OUTPUT,LUPRI,IPRINT)
ENDIF

#if MOD_FMM
IF (INPUT%DO_FMM) CALL FMMFCK(Output%ResultMat,INPUT%DMAT_RHS,WRK,LWRK)
#endif

 DeAllocate(OD_LHS%BATCH)
 DeAllocate(OD_RHS%BATCH)
 Nullify(OD_LHS%BATCH)
 Nullify(OD_RHS%BATCH)

CALL LSTIMER('LSINT ',TSTART,TEND,LUPRI)

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
TYPE(Overlap)         :: P
TYPE(Overlap),pointer :: Q(:)
TYPE(Integrand)       :: PQ
Integer               :: nPrimP
Integer,pointer       :: nPrimQ(:)

IF (IPRINT.GT.5) THEN
  CALL LSHEADER(LUPRI,'McMurchieDavidson')
  WRITE(LUPRI,'(1X,A)') ''
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') '***                     McMurchieDavidson                      ***'
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') ''
ENDIF

CALL initTUVitem(sharedTUV,Input,LUPRI,IPRINT)
CALL Build_PRECALCULATED_SPHMAT(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)
integral%TUV => sharedTUV

NULLIFY(Q)
ALLOCATE(Q(OD_RHS%nbatches))
NULLIFY(nPrimQ)
ALLOCATE(nPrimQ(OD_RHS%nbatches))
DO IRHS=1,OD_RHS%nbatches
  CALL SET_Overlap(Q(IRHS),nPrimQ(IRHS),Input,Integral,OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS')
  IF ((nPrimQ(IRHS) .GT. 0).AND.INPUT%setETUVoutside) &
  &  CALL AttachETUVtoOverlap(Q(IRHS),SharedTUV,LUPRI,IPRINT)
ENDDO


!$OMP PARALLEL DO PRIVATE(integral,P,PQ,ILHS,IRHS,Start_RHS,End_RHS,nPrimP) SCHEDULE(DYNAMIC,1)
DO ILHS = 1,OD_LHS%nbatches
 integral%TUV => sharedTUV
 CALL SET_Overlap(P,nPrimP,Input,Integral,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS')
 IF (nPrimP .GT. 0) THEN
  IF (INPUT%setETUVoutside) &
    &  CALL AttachETUVtoOverlap(P,SharedTUV,LUPRI,IPRINT)
  CALL Determine_RHS_loop(INPUT,OD_RHS,ILHS,Start_RHS,End_RHS)
  DO IRHS = Start_RHS,End_RHS
     IF(OD_LHS%BATCH(ILHS)%maxGAB*OD_RHS%BATCH(IRHS)%maxGAB &
          & .GT. INPUT%CS_THRESHOLD  .OR. .NOT. INPUT%CS_SCREEN)THEN
        IF (IPRINT.GT.3) WRITE(LUPRI,'(1X,A,I3,A,I3)') 'Overlap ditributions P',ILHS,' and Q',IRHS
        IF(nPrimQ(IRHS) .GT. 0)THEN
          CALL Build_Integrand(PQ,P,Q(IRHS),INPUT,ILHS,IRHS,LUPRI,IPRINT)
          !   Hermite 2-electron integral over general operator w 
          CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
          CALL Contract_Q(INTEGRAL,PQ,Input,LUPRI,IPRINT)
          CALL Contract_P(INTEGRAL,PQ,input,LUPRI,IPRINT)
!$OMP CRITICAL            
          CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
!$OMP END CRITICAL            
          CALL CLEAR_integral(INTEGRAL)
          DEALLOCATE(INTEGRAL%TUVQ)
          CALL FREE_integrand(PQ)
        ENDIF
     ENDIF
  ENDDO
  CALL FREE_OVERLAP(P)
 ENDIF
ENDDO
!$OMP END PARALLEL DO

DO IRHS=1,OD_RHS%nbatches
  IF (nPrimQ(IRHS) .GT. 0) CALL FREE_OVERLAP(Q(IRHS))
ENDDO

CALL freeTUVitem(sharedTUV,Input)
DEALLOCATE(nPrimQ)
DEALLOCATE(Q)

END SUBROUTINE McMurchieDavidson
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
TYPE(TUVitem),target  :: SharedTUV
TYPE(Overlap)         :: P
TYPE(Overlap),pointer :: Q(:)
TYPE(Overlap),pointer :: F(:)
TYPE(Integrand)       :: PQ
Integer               :: nPrimP, nPrimQ

IF (IPRINT.GT.5) THEN
  CALL LSHEADER(LUPRI,'Jengine')
  WRITE(LUPRI,'(1X,A)') ''
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') '***                          Jengine                           ***'
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') ''
ENDIF

CALL initTUVitem(sharedTUV,Input,LUPRI,IPRINT)
CALL Build_PRECALCULATED_SPHMAT(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)
integral%TUV => sharedTUV

NULLIFY(Q)
ALLOCATE(Q(OD_RHS%nbatches))
CALL SET_FTUVbatches(F,NFTUVbatches,OD_RHS,Q,Input,Integral,LUPRI,IPRINT)
DEALLOCATE(Q)
Q=>F
NULLIFY(F)

!$OMP PARALLEL DO PRIVATE(integral,P,PQ,ILHS,IRHS,Start_LHS,End_LHS,nPrimP,nPrimQ) SCHEDULE(DYNAMIC,1)
DO ILHS = 1,OD_LHS%nbatches
 integral%TUV => sharedTUV
 CALL SET_Overlap(P,nPrimP,Input,Integral,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS')
 IF (nPrimP.GT.0) THEN
  NULLIFY(Integral%TUVQ)
  ALLOCATE(Integral%TUVQ(P%nPrimitives,P%nTUV,1))
  CALL DZERO(Integral%TUVQ,P%nPrimitives*P%nTUV*1)
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
        CALL CLEAR_integral(INTEGRAL)
        CALL FREE_integrand(PQ)
     ENDIF
  ENDDO
  CALL Contract_P(INTEGRAL,PQ,input,LUPRI,IPRINT)
!$OMP CRITICAL            
  CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
!$OMP END CRITICAL            
  DEALLOCATE(Integral%TUVQ)
  CALL FREE_OVERLAP(P)
 ENDIF
ENDDO
!$OMP END PARALLEL DO

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
SUBROUTINE LINK(OD_LHS,OD_RHS,INPUT,OUTPUT,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
!
Integer               :: ILHS,IRHS,Start_RHS,End_RHS,i,j
TYPE(TUVitem),target  :: SharedTUV
TYPE(Integralitem)    :: Integral
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
Integer               :: LISTSIZE,LISTSIZE1,LISTSIZE2
Integer               :: nB,nC,nD,OLDI,RHSINDEX,startB,startD,I2,JK
Integer,allocatable   :: LIST(:),LIST2(:),LIST1(:)
real(realk),allocatable :: maxLHSGAB(:),maxRHSGAB(:)
!real(realk),allocatable :: SORTING(:)
!integer,allocatable     :: bshell(:)
real(realk)           :: maxLHSELM,maxRHSELM,DMATELM1,DMATELM2,MAXDMAT
!LOGICAL               :: A_LOOP_DONE,B_LOOP_DONE,C_LOOP_DONE,D_LOOP_DONE,UNIQUE
TYPE(LINKshell),pointer  :: brashell(:),ketshell(:),ML(:,:)
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
IF(.NOT.INPUT%DO_DALINK)THEN ! DO LINK
   NULLIFY(brashell)
   ALLOCATE(brashell(dim1))
   CALL DETERMINE_SHELL_PAIRS(dim1,dim2,RED_GAB_LHS,INPUT%CS_THRESHOLD,&
        &maxRHSELM,brashell,INPUT%sameLHSaos)
!   WRITE(LUPRI,*)'BRASHELL'
!   DO A=1,dim1
!      DO B=1,brashell(A)%DIM
!         WRITE(LUPRI,*)'(A,B)=(',A,',',B,')=',brashell(A)%elms(B)
!      ENDDO
!   ENDDO
   NULLIFY(ketshell)
   ALLOCATE(ketshell(dim3))
   CALL DETERMINE_SHELL_PAIRS(dim3,dim4,RED_GAB_RHS,INPUT%CS_THRESHOLD,&
        &maxLHSELM,ketshell,.FALSE.)
!   WRITE(LUPRI,*)'KETSHELL'
!   DO A=1,dim3
!      DO B=1,ketshell(A)%DIM
!         WRITE(LUPRI,*)'(C,D)=(',A,',',B,')=',ketshell(A)%elms(B)
!      ENDDO
!   ENDDO
ENDIF

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

IF(.NOT.INPUT%DO_DALINK)THEN ! DO LINK
   NULLIFY(ML)
   ALLOCATE(ML(dim1,Input%NDMAT_RHS))
   CALL DETERMINE_BRAKET_PAIRS(dim1,dim2,dim3,Input%NDMAT_RHS,RED_DMAT_RHS,&
        &maxLHSGAB,maxRHSGAB,INPUT%CS_THRESHOLD,ML)
ENDIF

CALL initTUVitem(sharedTUV,Input,LUPRI,IPRINT)
CALL Build_PRECALCULATED_SPHMAT(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)
integral%TUV => sharedTUV

NULLIFY(Q)
ALLOCATE(Q(OD_RHS%nbatches))
NULLIFY(nPrimQ)
ALLOCATE(nPrimQ(OD_RHS%nbatches))
DO IRHS=1,OD_RHS%nbatches
   CALL SET_Overlap(Q(IRHS),nPrimQ(IRHS),Input,Integral,OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS')
   IF ((nPrimQ(IRHS) .GT. 0).AND.INPUT%setETUVoutside) &
        &  CALL AttachETUVtoOverlap(Q(IRHS),SharedTUV,LUPRI,IPRINT)
ENDDO

NULLIFY(P)
ALLOCATE(P(OD_LHS%nbatches))
NULLIFY(nPrimP)
ALLOCATE(nPrimP(OD_LHS%nbatches))
DO ILHS=1,OD_LHS%nbatches
  CALL SET_Overlap(P(ILHS),nPrimP(ILHS),Input,Integral,OD_LHS%BATCH(ILHS),2,LUPRI,IPRINT,'LHS')
  IF ((nPrimP(ILHS) .GT. 0).AND.INPUT%setETUVoutside) &
  &  CALL AttachETUVtoOverlap(P(ILHS),SharedTUV,LUPRI,IPRINT)
ENDDO

IF(.NOT.INPUT%DO_DALINK)THEN
idmat=1!,Input%NDMAT_RHS !this needs to be addressed
!$OMP PARALLEL DO PRIVATE(integral,PQ,IRHS,RHSINDEX,A,nB,B,ILHS,I,nC,C,OLDI,nD,D,&
!$OMP LIST,LIST1,LIST2,LISTSIZE1,LISTSIZE2,LISTSIZE) SCHEDULE(DYNAMIC,1)
   DO A=1,dim1
      ALLOCATE(LIST(dim3*dim4))
      ALLOCATE(LIST1(dim3*dim4))
      ALLOCATE(LIST2(dim3*dim4))
      DO nB=1,brashell(A)%DIM
         B=brashell(A)%elms(nB)
         ILHS=(A-1)*dim1+B
         IF(INPUT%sameLHSaos) ILHS=(A-1)*dim1-(A*(A-1)/2)+B 
         IF(nPrimP(ILHS) .GT. 0)THEN
            !CREATING FIRST LIST from the AC element of the density matrix
            I=0
            DO nC=1,ML(A,idmat)%DIM
               C = ML(A,idmat)%elms(nC) 
               OLDI=I
               IF(INPUT%sameRHSaos)THEN
                  DO nD=1,ketshell(C)%DIM
                     D=ketshell(C)%elms(nD)
                     IF(RED_DMAT_RHS(A,C,idmat)*RED_GAB_RHS(A,B)&
                          &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD )EXIT
                     IF(D .GE. C)THEN
                        I=I+1
                        LIST1(I)=(C-1)*dim1-(C*(C-1)/2)+D !=IRHS
                     ELSE
                        I=I+1
                        LIST1(I)=(D-1)*dim1-(D*(D-1)/2)+C !=IRHS     
                     ENDIF
                  ENDDO
               ELSE
                  DO nD=1,ketshell(C)%DIM
                     D=ketshell(C)%elms(nD)
                     IF(RED_DMAT_RHS(A,C,idmat)*RED_GAB_RHS(A,B)&
                          &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD )EXIT
                     I=I+1
                     LIST1(I)=(C-1)*dim1+D
                  ENDDO
               ENDIF
               IF(I .EQ. OLDI)EXIT !NO ELEMENTS ADDED IN D LOOP
            ENDDO
            !THIS LIST CAN CONTAIN THE SAME ELEMENT 2 TIMES IF sameRHSaos is true
            !THIS WIL BE HANDLED IN MERGELIST
            LISTSIZE1=I
            !CREATING SECOND LIST from the BC element of the density matrix
            I=0
            DO nC=1,ML(A,idmat)%DIM
               C = ML(A,idmat)%elms(nC)
               OLDI=I
               IF(INPUT%sameRHSaos)THEN
                  DO nD=1,ketshell(C)%DIM
                     D=ketshell(C)%elms(nD)
                     IF(RED_DMAT_RHS(B,C,idmat)*RED_GAB_RHS(A,B)&
                          &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD ) EXIT
                     IF(D .GE. C)THEN
                        I=I+1
                        LIST2(I)=(C-1)*dim1-(C*(C-1)/2)+D
                     ELSE
                        I=I+1
                        LIST2(I)=(D-1)*dim1-(D*(D-1)/2)+C
                     ENDIF
                  ENDDO
               ELSE
                  DO nD=1,ketshell(C)%DIM
                     D=ketshell(C)%elms(nD)
                     IF(RED_DMAT_RHS(B,C,idmat)*RED_GAB_RHS(A,B)&
                          &*RED_GAB_RHS(C,D) .LE. INPUT%CS_THRESHOLD ) EXIT
                     I=I+1
                     LIST2(I)=(C-1)*dim1+D
                  ENDDO
               ENDIF
               IF(I .EQ. OLDI)EXIT !NO ELEMENTS ADDED IN D LOOP
            ENDDO
            LISTSIZE2=I
            CALL MERGE_LIST(LUPRI,LIST,LIST1(1:LISTSIZE1),LIST2(1:LISTSIZE2),&
                 &LISTSIZE,LISTSIZE1,LISTSIZE2,dim3*dim4,ILHS,input%sameODs)
            DO RHSINDEX=1,LISTSIZE
               integral%TUV => sharedTUV
               IRHS=LIST(RHSINDEX)
               IF(nPrimQ(IRHS) .GT. 0)THEN
                  CALL Build_Integrand(PQ,P(ILHS),Q(IRHS),INPUT,ILHS,IRHS,LUPRI,IPRINT)
                  !   Hermite 2-electron integral over general operator w 
                  CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
                  CALL Contract_Q(INTEGRAL,PQ,Input,LUPRI,IPRINT)
                  CALL Contract_P(INTEGRAL,PQ,input,LUPRI,IPRINT)
                  !$OMP CRITICAL            
                  CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
                  !$OMP END CRITICAL            
                  CALL CLEAR_integral(INTEGRAL)
                  DEALLOCATE(INTEGRAL%TUVQ)
                  CALL FREE_integrand(PQ)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      DEALLOCATE(LIST)
      DEALLOCATE(LIST1)
      DEALLOCATE(LIST2)
   ENDDO
!$OMP END PARALLEL DO
!ENDDO
ELSE
   idmat=1
!$OMP PARALLEL DO PRIVATE(integral,PQ,IRHS,A,STARTB,B,ILHS,C,startD,D,DMATELM1,&
!$OMP DMATELM2,MAXDMAT) SCHEDULE(DYNAMIC,1)
   DO A=1,dim1
     STARTB=1
     IF(INPUT%sameLHSaos) STARTB=A 
     DO B=STARTB,dim2
       ILHS=(A-1)*dim1+B
       IF(INPUT%sameLHSaos) ILHS=(A-1)*dim1-(A*(A-1)/2)+B 
       IF(nPrimP(ILHS) .GT. 0)THEN
         DO C=1,dim3
           STARTD=1
           IF(INPUT%sameRHSaos) STARTD=C
           DO D=STARTD,dim4
             IRHS=(C-1)*dim1+D
             IF(INPUT%sameRHSaos) IRHS=(C-1)*dim1-(C*(C-1)/2)+D 
             IF(IRHS .LE. ILHS)THEN
               DMATELM1=RED_DMAT_RHS(A,C,idmat)*RED_DMAT_LHS(B,D,idmat)
               DMATELM2=RED_DMAT_RHS(B,C,idmat)*RED_DMAT_LHS(A,D,idmat)
               MAXDMAT=MAX(DMATELM1,DMATELM2)
               IF(MAXDMAT*RED_GAB_LHS(A,B)*RED_GAB_RHS(C,D) .GE. INPUT%CS_THRESHOLD )THEN
                 IF(nPrimQ(IRHS) .GT. 0)THEN
                   CALL Build_Integrand(PQ,P(ILHS),Q(IRHS),INPUT,ILHS,IRHS,LUPRI,IPRINT)
                   !   Hermite 2-electron integral over general operator w 
                   CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
                   CALL Contract_Q(INTEGRAL,PQ,Input,LUPRI,IPRINT)
                   CALL Contract_P(INTEGRAL,PQ,input,LUPRI,IPRINT)
                   !$OMP CRITICAL            
                   CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
                   !$OMP END CRITICAL            
                   CALL CLEAR_integral(INTEGRAL)
                   DEALLOCATE(INTEGRAL%TUVQ)
                   CALL FREE_integrand(PQ)
                 ENDIF
               ENDIF
             ENDIF
           ENDDO
         ENDDO
       ENDIF
     ENDDO
   ENDDO
!$OMP END PARALLEL DO
ENDIF

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
IF(.NOT.INPUT%DO_DALINK)THEN ! DO LINK
   DEALLOCATE(brashell)
   NULLIFY(brashell)
   DEALLOCATE(ketshell)
   NULLIFY(ketshell)
   DEALLOCATE(ML)
   NULLIFY(ML)
ENDIF
DEALLOCATE(RED_DMAT_RHS)
NULLIFY(RED_DMAT_RHS)
IF(INPUT%DO_DALINK)THEN
   DEALLOCATE(RED_DMAT_LHS)
   NULLIFY(RED_DMAT_LHS)
ENDIF

END SUBROUTINE LINK

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
          RED_GAB(IA,IB)=MAX(RED_GAB(IA,IB),ABS(GAB(sA+KA+(CA-1)*KA-1,sB+KB+(CB-1)*KB-1)))
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

SUBROUTINE MERGE_LIST(LUPRI,LIST,LIST1,LIST2,LISTSIZE,LISTSIZE1,LISTSIZE2,SIZE,ILHS,sameODs)
IMPLICIT NONE
Integer  :: SIZE,LISTSIZE1,LISTSIZE2,LISTSIZE
Integer  :: LIST(SIZE),LIST1(LISTSIZE1),LIST2(LISTSIZE2),LUPRI
!
Integer  :: I,J,K,ILHS
logical  :: UNIQUE,sameODs

K=0
DO J=1,LISTSIZE1
   IF(LIST1(J) .LE. ILHS)THEN
      UNIQUE=.TRUE.
      DO I=1,K
         IF(LIST(I).EQ.LIST1(J)) UNIQUE=.FALSE.
      ENDDO
      IF(UNIQUE) THEN
         K=K+1
         LIST(K)=LIST1(J)
      ENDIF
   ENDIF
ENDDO

!SECOND LIST
DO J=1,LISTSIZE2
   IF(LIST2(J) .LE. ILHS)THEN
      UNIQUE=.TRUE.
      DO I=1,K
         IF(LIST(I).EQ.LIST2(J)) UNIQUE=.FALSE.
      ENDDO
      IF(UNIQUE) THEN
         K=K+1
         LIST(K)=LIST2(J)
      ENDIF
   ENDIF
ENDDO
LISTSIZE=K

END SUBROUTINE MERGE_LIST

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
   I=0
   StartB=1
   IF(sameAOs)StartB=A
   DO B=StartB,dim2
      IF(RED_GAB_LHS(A,B) .GT. CS_THRESHOLD/maxRHSELM)THEN
         I=I+1
         bshell(I)=B
         SORTING(I)=RED_GAB_LHS(A,B)
      ENDIF
   ENDDO
   brashell(A)%DIM=I
   CALL Qsort(SORTING(1:I),bshell(1:I))
   NULLIFY(brashell(A)%elms)
   ALLOCATE(brashell(A)%elms(I))
   DO C=1,I
      brashell(A)%elms(C)=bSHELL(C)
   ENDDO
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
      DO C=1,dim3
         IF(RED_DMAT_RHS(A,C,idmat)*maxLHSGAB(A)*maxRHSGAB(C).GT. CS_THRESHOLD)THEN
            I=I+1
            bshell(I)=C
            SORTING(I)=RED_DMAT_RHS(A,C,idmat)*maxRHSGAB(C)
         ENDIF
      ENDDO
      ML(A,idmat)%DIM=I
      CALL Qsort(SORTING(1:I),bshell(1:I))
      NULLIFY(ML(A,idmat)%elms)
      ALLOCATE(ML(A,idmat)%elms(dim3))
      DO C=1,I
         ML(A,idmat)%elms(C)=bSHELL(C)
      ENDDO
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
Integer              :: LUPRI,IPRINT
!simen Select different cases here: contraction with densities etc.

  ! No pre-contraction with densities (either full integrals, or contraction
  ! with full integrals)
  CALL DistributePQint(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
  DEALLOCATE(INTEGRAL%PQ)
  NULLIFY(INTEGRAL%PQ)

END SUBROUTINE DistributeIntegrals

SUBROUTINE DistributePQint(Integral,PQ,Input,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: i1,i2,i3,i4,i5
!
nAngmomP = PQ%P%p%nAngmom

nderivP=Input%nDerivP
!

iOrbitalP = 1
DO iDerivP=1,NDerivP
   DO iAngmomP=1,nAngmomP
      nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
      totOrbQ = PQ%Q%p%totOrbitals
      endOrbP = iOrbitalP+nOrbP-1
      CALL DistributePQin1(Integral%PQ(:,iOrbitalP:endOrbP),nOrbP,totOrbQ, &
           &               PQ%P%p,PQ%Q%p,Input,Output,iAngmomP,iDerivP,NderivP,&
           &               PQ%sameOD,LUPRI,IPRINT)
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

SUBROUTINE DistributePQin1(QP,nOrbP,totOrbQ,P,Q,Input,Int_Output,iAngmomP,iDerivP,NderivP,samePQ,LUPRI,IPRINT)  
implicit none
Integer              :: nOrbP,totOrbQ,iAngmomP,LUPRI,IPRINT,iderivP,nderivP
Real(realk)          :: QP(totOrbQ,nOrbP)
Type(Overlap)        :: P,Q
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Int_Output
Logical              :: samePQ
!
Integer                 :: nAngmomQ,nOrbQ,nderivQ
Integer                 :: iAngmomQ,iOrbitalQ,iderivQ
Real(realk)             :: PQ(nOrbP,totOrbQ)
Integer                 :: iOrbP,iOrbQ,endAngmomQ
!
!NULLIFY(PQ)
!ALLOCATE(PQ(nOrbP,totOrbQ))
! Transposition
DO iOrbP=1,nOrbP
  DO iOrbQ=1,totOrbQ
    PQ(iOrbP,iOrbQ) = QP(iOrbQ,iOrbP)
  ENDDO
ENDDO
!
iOrbitalQ = 1
endAngmomQ = Q%nAngmom
IF (samePQ) endAngmomQ = iAngmomP
nderivQ=input%nderivQ
DO iderivQ =1,nderivQ
   DO iAngmomQ=1,endAngmomQ
      nOrbQ = Q%nOrbitals(iAngmomQ)
      CALL DistributePQin2(PQ(1,iOrbitalQ),P,Q,nOrbP,nOrbQ,iAngmomP,iAngmomQ,iDerivP,NderivP,iDerivQ,Input,Int_Output,LUPRI,IPRINT)
      iOrbitalQ = iOrbitalQ + nOrbQ
   ENDDO
ENDDO
!DEALLOCATE(PQ)
END SUBROUTINE DistributePQin1

!SUBROUTINE DistributePQin2(PQ(1,iOrbitalQ),P,Q,iAngmomP,iAngmomQ,Input,Output,LUPRI,IPRINT)

SUBROUTINE DistributePQin2(PQ,P,Q,nOrbP,nOrbQ,iAngmomP,iAngmomQ,iDerivP,NderivP,iDerivQ,Input,Output,LUPRI,IPRINT)
implicit none
Real(realk)          :: PQ(nOrbP,nOrbQ)
Type(Overlap)        :: P,Q
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: nOrbP,nOrbQ,iAngmomP,iAngmomQ,LUPRI,IPRINT,iderivP,iderivQ
!
Integer :: iA,iB,iC,iD,NderivP
Integer :: nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD
Integer :: startA,startB,startC,startD
!
iA = P%indexAng1(iAngmomP)
iB = P%indexAng2(iAngmomP)
iC = Q%indexAng1(iAngmomQ)
iD = Q%indexAng2(iAngmomQ)
nContA = P%orbital1%nContracted(iA)
nContB = P%orbital2%nContracted(iB)
nContC = Q%orbital1%nContracted(iC)
nContD = Q%orbital2%nContracted(iD)
nAngA = P%orbital1%nOrbComp(iA)
nAngB = P%orbital2%nOrbComp(iB)
nAngC = Q%orbital1%nOrbComp(iC)
nAngD = Q%orbital2%nOrbComp(iD)
startA  = P%orbital1%startOrbital(iA)
startB  = P%orbital2%startOrbital(iB)
startC  = Q%orbital1%startOrbital(iC)
startD  = Q%orbital2%startOrbital(iD)
CALL DistributePQin3(PQ,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,&
      &      startA,startB,startC,startD,iderivP,NderivP,iderivQ,Output,INPUT)
END SUBROUTINE DistributePQin2

SUBROUTINE DistributePQin3(PQ,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,&
      &                    startA,startB,startC,startD,iderivP,NderivP,iderivQ,Output,INPUT)
implicit none
Real(realk)          :: PQ(nAngA,nAngB,nContA*nContB,nAngC,nAngD,nContC*nContD)
!Real(realk)          :: abcd(nAngA*nContA,nAngB*nContB,nAngC*nContC,nAngD*nContD)
Type(IntegralOutput) :: Output
Type(IntegralInput)  :: Input
Integer              :: nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD
Integer              :: startA,startB,startC,startD,iderivP,iderivQ
!
Integer     :: iA,iB,iC,iD,NderivP
Integer     :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
Integer     :: iContP,iContQ
Real(realk) :: integral,dtemp_lhs,dtemp_rhs
Integer     :: idmat,ideriv
Integer     :: n1,n2
!
!IF(INPUT%sameLHSaos .AND. INPUT%sameRHSaos) THEN
ideriv=iderivP+(iderivQ-1)*NderivP

iA=startA
DO iContA=1,nContA
 DO iAngA=1,nAngA
  iB=startB
  DO iContB=1,nContB
   iContP = (iContB-1)*nContA + iContA
   DO iAngB=1,nAngB
!   Calculation of screening integrals
    IF (Input%CS_int) THEN
      IF ((startA.EQ.startC).AND.(startB.EQ.startD)) THEN
        integral = sqrt(PQ(iAngA,iAngB,iContP,iAngA,iAngB,iContP))
        Output%ResultMat(iA,iB,1,1,ideriv) = integral
        IF (INPUT%sameLHSaos) Output%ResultMat(iB,iA,1,1,ideriv) = integral
      ENDIF
!   Regular case
    ELSE
     iC=startC
     DO iContC=1,nContC
      DO iAngC=1,nAngC
       iD=startD
       DO iContD=1,nContD
        iContQ = (iContD-1)*nContC + iContC
        DO iAngD=1,nAngD
         integral = PQ(iAngA,iAngB,iContP,iAngC,iAngD,iContQ)
!        Contract with density
         IF (Input%DO_FOCK) THEN
          DO idmat=1,Input%NDMAT_RHS
!          idmat = 1
           IF (Input%DO_Coulomb) THEN
!simen change Dmat dimensions
             n1 = Output%ndim(1)
             n2 = Output%ndim(2)
             CALL AddToCoulombMat(Output%ResultMat(:,:,1,1,idmat),n1,n2,&
     &                            INPUT,integral,iA,iB,iC,iD,startA,startB,startC,startD,&
     &                            Input%DMAT_RHS(:,:,idmat),n1,n2)
           ENDIF
           IF (Input%DO_Exchange) THEN
!simen change Dmat dimensions
             n1 = Output%ndim(1)
             n2 = Output%ndim(2)
             CALL AddToExchangeMat(Output%ResultMat(:,:,1,1,idmat),n1,n2,&
     &            INPUT,integral,iA,iB,iC,iD,startA,startB,startC,startD,&
     &            Input%DMAT_RHS(:,:,idmat),n1,n2)
           ENDIF
          ENDDO
!        No contraction with density
         ELSE
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
         ENDIF
         iD = iD+1
        ENDDO
       ENDDO
       iC = iC+1
      ENDDO
     ENDDO
    ENDIF  !Screening or regular integrals
    iB = iB+1
   ENDDO
  ENDDO
  iA = iA+1
 ENDDO
ENDDO

#if 0
ELSEIF(INPUT%sameRHSaos)THEN
   iA=startA
   DO iContA=1,nContA
      DO iAngA=1,nAngA
         iB=startB
         DO iContB=1,nContB
            iContP = (iContB-1)*nContA + iContA
            DO iAngB=1,nAngB
               iC=startC
               DO iContC=1,nContC
                  DO iAngC=1,nAngC
                     iD=startD
                     DO iContD=1,nContD
                        iContQ = (iContD-1)*nContC + iContC
                        DO iAngD=1,nAngD
                           Output%ResultMat(iA,iB,iC,iD,ideriv) = &
                                PQ(iAngA,iAngB,iContP,iAngC,iAngD,iContQ)
                           IF(iC .NE. iD)THEN
                              Output%ResultMat(iA,iB,iD,iC,ideriv) = &
                                   PQ(iAngA,iAngB,iContP,iAngC,iAngD,iContQ)
                           ENDIF
                           iD = iD+1
                        ENDDO
                     ENDDO
                     iC = iC+1
                  ENDDO
               ENDDO
               iB = iB+1
            ENDDO
         ENDDO
         iA = iA+1
      ENDDO
   ENDDO
elseif(INPUT%sameLHSaos)THEN
   iA=startA
   DO iContA=1,nContA
      DO iAngA=1,nAngA
         iB=startB
         DO iContB=1,nContB
            iContP = (iContB-1)*nContA + iContA
            DO iAngB=1,nAngB
               iC=startC
               DO iContC=1,nContC
                  DO iAngC=1,nAngC
                     iD=startD
                     DO iContD=1,nContD
                        iContQ = (iContD-1)*nContC + iContC
                        DO iAngD=1,nAngD
                           Output%ResultMat(iA,iB,iC,iD,ideriv) = &
                                PQ(iAngA,iAngB,iContP,iAngC,iAngD,iContQ)
                           IF(iA .NE. iB)THEN
                              Output%ResultMat(iB,iA,iC,iD,ideriv) = &
                                   PQ(iAngA,iAngB,iContP,iAngC,iAngD,iContQ)
                           ENDIF
                           iD = iD+1
                        ENDDO
                     ENDDO
                     iC = iC+1
                  ENDDO
               ENDDO
               iB = iB+1
            ENDDO
         ENDDO
         iA = iA+1
      ENDDO
   ENDDO
else
   iA=startA
   DO iContA=1,nContA
      DO iAngA=1,nAngA
         iB=startB
         DO iContB=1,nContB
            iContP = (iContB-1)*nContA + iContA
            DO iAngB=1,nAngB
               iC=startC
               DO iContC=1,nContC
                  DO iAngC=1,nAngC
                     iD=startD
                     DO iContD=1,nContD
                        iContQ = (iContD-1)*nContC + iContC
                        DO iAngD=1,nAngD
                           Output%ResultMat(iA,iB,iC,iD,ideriv) = &
                                PQ(iAngA,iAngB,iContP,iAngC,iAngD,iContQ)
                           iD = iD+1
                        ENDDO
                     ENDDO
                     iC = iC+1
                  ENDDO
               ENDDO
               iB = iB+1
            ENDDO
         ENDDO
         iA = iA+1
      ENDDO
   ENDDO
ENDIF
#endif

END SUBROUTINE DistributePQin3

SUBROUTINE SET_OVERLAP(P,np,Input,Integral,ODB,IELECTRON,IUNIT,IPRINT,SIDE)
Implicit none
TYPE(IntegralInput) :: Input
TYPE(IntegralItem)  :: Integral
TYPE(Overlap)       :: P
TYPE(ODBATCH)       :: ODB
Character*(*)       :: side
Integer             :: IELECTRON,IUNIT,IPRINT,nderiv,np
!
Integer             :: idir,i1,i2,i12,start1,start2,end1,end2,a,b,a1,a2,orb1,orb2
Real(realk)         :: e1,e2,d2,maxGab,maxPrimGab,maxPrimGabElm
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
   WRITE(IUNIT,'(1X,2A)') 'Wrong case in SET_OVERLAP =',SIDE
   CALL QUIT('Wrong case in SET_OVERLAP')
END SELECT

 CALL SET_ORBITAL(P%orbital1,ODB%AO(1)%p,integral,IUNIT)
 CALL SET_ORBITAL(P%orbital2,ODB%AO(2)%p,integral,IUNIT)

!Determine overlap type
 t1 = P%orbital1%type
 t2 = P%orbital2%type
 IF ( t1.EQ.'Empty'.AND.t2.EQ.'Empty') THEN
   P%type = 'Empty'
 ELSE IF ( (t1.EQ.'Hermite'.AND.t2.EQ.'Empty').OR.(t1.EQ.'Empty'  .AND.t2.EQ.'Hermite') ) THEN
   P%type = 'Hermite-single'
 ELSE IF (t1.EQ.'Hermite'.AND.t2.EQ.'Hermite') THEN
   P%type = 'Hermite'
 ELSE IF ( (t1.EQ.'Cartesian'.AND.t2.EQ.'Empty').OR.(t1.EQ.'Empty'  .AND.t2.EQ.'Cartesian') ) THEN
   P%type = 'Cartesian-singel'
 ELSE IF (t1.EQ.'Cartesian'.AND.t2.EQ.'Cartesian') THEN
   P%type = 'Cartesian'
 ELSE
   WRITE(IUNIT,'(1X,A)') 'Not a proper combination of orbital types in SET_OVERLAP:'
   WRITE(IUNIT,'(5X,4A)') 'orbital1%type =', t1, 'and orbital2%type =', t2
   CALL QUIT('Not a proper combination of orbital types in SET_OVERLAP')
 ENDIF

 P%sameAO = ODB%sameAO
 P%ETUVisSet = .FALSE.
 P%sphericalETUV = .TRUE.
 CALL GET_NPRIMITIVES(np,P,Input,Side)
 P%nPrimitives   = np
 IF (np.EQ.0) THEN
   CALL FREE_ORBITAL(P%orbital1)
   CALL FREE_ORBITAL(P%orbital2)
   RETURN
 ENDIF

 P%nAngmom       = ODB%nAngmom
 P%maxContracted = ODB%maxContracted
!
 d2 = 0.0d0
 DO idir=1,3
   P%distance12(idir) = P%orbital1%center(idir)-P%orbital2%center(idir)
   d2 = d2 + P%distance12(idir) * P%distance12(idir)
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
              a = P%orbital1%startprimOrbital(a1)+(i1-1)*P%orbital1%nOrbComp(a1)
              b = P%orbital2%startprimOrbital(a2)+(i2-1)*P%orbital2%nOrbComp(a2)
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
   WRITE(IUNIT,'(1X,A,2I5)') 'Error in set_overlap. Mismatch in number of primitives, (np,i12)=',&
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
 DO i1=1,P%orbital1%nAngmom
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nAngmom
     i12 = i12 + 1
     P%angmom(i12)      = P%orbital1%angmom(i1) + P%orbital2%angmom(i2)
     P%indexAng1(i12)   = i1
     P%indexAng2(i12)   = i2
     P%nOrbitals(i12)   = P%orbital1%nOrbitals(i1)*P%orbital2%nOrbitals(i2)
     P%totOrbitals      = P%totOrbitals + P%nOrbitals(i12)*NDERIV
     P%nOrbComp(i12)    = P%orbital1%nOrbComp(i1)*P%orbital2%nOrbComp(i2)
     P%nContracted(i12) = P%orbital1%nContracted(i1)*P%orbital2%nContracted(i2)
     P%minAngmom        = min(P%minAngmom,P%angmom(i12))
     P%maxAngmom        = max(P%maxAngmom,P%angmom(i12))
     start1 = P%orbital1%startOrbital(i1)
     end1   = start1 + P%orbital1%nOrbitals(i1) - 1
     start2 = P%orbital2%startOrbital(i2)
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
 IF (P%type.EQ.'Hermite-single') P%startAngmom = P%minAngmom

 IF( INPUT%operator(1:7) .EQ. 'Kinetic' .AND. LHS) THEN
    P%startAngmom = P%startAngmom + 2
    P%endAngmom   = P%maxAngmom + 2
 ELSE
    P%endAngmom   = P%maxAngmom
 ENDIF  

 P%nTUV = (P%endAngmom+1)*(P%endAngmom+2)*(P%endAngmom+3)/6 - &
    &     P%startAngmom*(P%startAngmom+1)*(P%startAngmom+2)/6

IF (INPUT%DO_JENGINE .AND. IELECTRON.EQ.2) THEN
  CALL SET_FTUV(P,INPUT,Integral,IUNIT,IPRINT)
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
  P%orbital1%startOrbital(1) = 1
  P%orbital2%startOrbital(1) = 1
  P%indexAng1(1)   = 1
  P%indexAng2(1)   = 1
  P%nOrbitals(1)   = 1
  P%totOrbitals    = 1
  P%nOrbComp(1)    = 1
  P%nContracted(1) = 1
ELSE
!GET E COEFFICIENTS
 IF (P%type.NE.'Hermite-single') THEN
!  CALL GET_ECOEFF(P,IUNIT,IPRINT)
 ENDIF
ENDIF

IF(INPUT%CS_SCREEN)THEN
   NULLIFY(GAB)
ENDIF
END SUBROUTINE SET_OVERLAP

SUBROUTINE AttachETUVtoOverlap(P,TUV,LUPRI,IPRINT)
Implicit none
TYPE(Overlap) :: P
TYPE(TUVitem) :: TUV
Integer       :: LUPRI,IPRINT
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
  LENGTH=nTUV*ijk*P%nPrimitives
  CALL BuildEcoeffTensor(TUV,P,P%ETUV(indexETUV:indexETUV+LENGTH-1),nTUV,ijk,P%nPrimitives,iAngmom,LUPRI,IPRINT)
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
Real(realk)         :: D1=1.0d0,D0=0.0d0,dtemp
Integer             :: eS1,eS2

NULLIFY(P%FTUV)
ALLOCATE(P%FTUV(P%nPrimitives,P%nTUV,Input%NDMAT_RHS))
CALL DZERO(P%FTUV,P%nTUV*P%nPrimitives*Input%NDMAT_RHS)

DO iAngmom = 1,P%nAngmom
  IF (P%type.EQ.'Hermite-single') THEN
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
    CALL BuildEcoeffTensor(Integral%TUV,P,Ecoeffs,nTUV,ijk,P%nPrimitives,iAngmom,LUPRI,IPRINT)

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
    CALL ConstructContraction(CC,P,P%nPrimitives,nC1,nC2,iA1,iA2,LUPRI,IPRINT)

!   Make contraction of density-matrix and contraction coefficients
    NULLIFY(ContractedDmat)
    ALLOCATE(ContractedDmat(lm,P%nPrimitives,Input%NDMAT_RHS))
    CALL DZERO(ContractedDmat,lm*P%nPrimitives*Input%NDMAT_RHS)
    start1 = P%orbital1%startOrbital(iA1)
    start2 = P%orbital2%startOrbital(iA2)
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
                ContractedDmat(indlm,iPrim,idmat) = ContractedDmat(indlm,iPrim,idmat) &
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
      DO iPrim=1,P%nPrimitives
        DO indlm=1,lm
          DO iTUV=1,nTUV
            P%FTUV(iPrim,iTUV,idmat) = P%FTUV(iPrim,iTUV,idmat) &
     &               + ContractedDmat(indlm,iPrim,idmat)*SpherEcoeff(iTUV,indlm,iPrim)
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
DO i1=1,P%orbital1%nPrimitives
  start2 = 1
  IF (P%sameAO) start2 = i1
  DO i2=start2,P%orbital2%nPrimitives
!   Find maxgab for given primitive
    maxgab = 0.0d0
    IF(INPUT%PS_SCREEN)THEN
       DO a1=1,P%orbital1%nAngmom 
          DO a2=1,P%orbital2%nAngmom
             a = P%orbital1%startprimOrbital(a1)+(i1-1)*P%orbital1%nOrbComp(a1)
             b = P%orbital2%startprimOrbital(a2)+(i2-1)*P%orbital2%nOrbComp(a2)
             maxgab=MAX(maxgab,GAB(a,b))
          ENDDO
       ENDDO
    ENDIF
!   Add if maxgab greater than threshold
    IF(.NOT.INPUT%PS_SCREEN .OR. MAXGAB .GT. INPUT%PS_THRESHOLD/MAXELM)THEN
      nPrimitives = nPrimitives + 1
    ENDIF  
  ENDDO
ENDDO
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
 Orb%maxContracted = AOB%maxContracted
 Orb%maxAngmom     = AOB%maxAngmom
 Orb%nAngmom       = na

 NULLIFY(Orb%exponents)
 ALLOCATE(Orb%exponents(np))
 Orb%exponents     = AOB%pExponents%elms(1:np)

 NULLIFY(Orb%angmom)
 ALLOCATE(Orb%angmom(na))
 Orb%angmom        = AOB%angmom(1:na)

 NULLIFY(Orb%nContracted)
 ALLOCATE(Orb%nContracted(na))
 Orb%nContracted   = AOB%nContracted(1:na)

 NULLIFY(Orb%startOrbital)
 ALLOCATE(Orb%startOrbital(na))
 Orb%startOrbital  = AOB%startOrbital(1:na)

 NULLIFY(Orb%startprimOrbital)
 ALLOCATE(Orb%startprimOrbital(na))
 Orb%startprimOrbital  = AOB%startprimOrbital(1:na)

 NULLIFY(Orb%nOrbComp)
 ALLOCATE(Orb%nOrbComp(na))
 Orb%nOrbComp  = AOB%nOrbComp(1:na)

 NULLIFY(Orb%nPrimOrbComp)
 ALLOCATE(Orb%nPrimOrbComp(na))
 Orb%nPrimOrbComp  = AOB%nPrimOrbComp(1:na)

 NULLIFY(Orb%nOrbitals)
 ALLOCATE(Orb%nOrbitals(na))
 Orb%nOrbitals  = AOB%nOrbitals(1:na)

 NULLIFY(Orb%CC)
 ALLOCATE(Orb%CC(np,AOB%maxContracted,na))
 DO ia=1,na
   nc = Orb%nContracted(ia)
   CALL DCOPY(np*nc,AOB%pCC(ia)%p%elms,1,Orb%CC(1,1,ia),1)
 ENDDO

 NULLIFY(Orb%SPH_MAT)
 ALLOCATE(Orb%SPH_MAT(Orb%maxAngmom+1))
 DO I=1,Orb%maxAngmom+1
    Orb%SPH_MAT(I)%p => integral%TUV%SPH_MAT(I)
 ENDDO

END SUBROUTINE SET_ORBITAL

SUBROUTINE PRINT_OVERLAP(P,IUNIT,IPRINT)
TYPE(Overlap) :: P
Integer       :: IUNIT,IPRINT
!
integer :: i
CALL LSHEADER(IUNIT,'OVERLAP')
WRITE(IUNIT,'(1X,A)') ''
WRITE(IUNIT,'(1X,A)') '*************************'
WRITE(IUNIT,'(1X,A)') '***      OVERLAP      ***'
WRITE(IUNIT,'(1X,A)') '*************************'
WRITE(IUNIT,'(1X,A)') ''
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
 WRITE(IUNIT,'(5X,A,I3)') 'totOrbitals      =', P%totOrbitals
 WRITE(IUNIT,'(5X,A,I3)') 'nTUV             =', P%nTUV
 WRITE(IUNIT,'(5X,A,I3)') 'maxContracted    =', P%maxContracted
 WRITE(IUNIT,'(5X,A,I3)') 'minAngmom        =', P%minAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'maxAngmom        =', P%maxAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'startAngmom      =', P%startAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'endAngmom        =', P%endAngmom
IF(P%type .NE. 'FTUV')THEN
   WRITE(IUNIT,'(5X,A,3F8.4)')'distance12 (A)   =', (P%distance12(i), i=1,3)
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
   WRITE(IUNIT,'(3X,A)')    '- iPrim   Px      Py      Pz      p       mu      pre    iprim1 iprim2 -'
   DO i=1,P%nPrimitives
      WRITE(IUNIT,'(5X,I4,6F8.4,2I7)') i, P%center(1,i), P%center(2,i), P%center(3,i), &
           &  P%exponents(i), P%reducedExponents(i), P%preExpFac(i), P%iprim1(i), P%iprim2(i)
   ENDDO
ELSE
   WRITE(IUNIT,'(3X,A)')    '- iPrim   Px      Py      Pz      p       mu      pre'
   DO i=1,P%nPrimitives
      WRITE(IUNIT,'(5X,I4,6F8.4)') i, P%center(1,i), P%center(2,i), P%center(3,i), &
           &  P%exponents(i), P%reducedExponents(i), P%preExpFac(i)
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
   WRITE(IUNIT,'(1X,8I8)')  I, Orb%angmom(I), Orb%nContracted(I), Orb%startOrbital(I), &
     & Orb%nOrbitals(I), Orb%nOrbComp(I), Orb%nPrimOrbComp(I),Orb%startprimOrbital(I)
 ENDDO
 WRITE(IUNIT,'(3X,A)')    '-----------------------------------------------------------'
 WRITE(IUNIT,'(3X,A)')    '----------------- Exponents -----------------'
 WRITE(IUNIT,'(5X,5F12.6)') (Orb%exponents(I),I=1,Orb%nPrimitives)
 WRITE(IUNIT,'(3X,A)')    '----------------- Contraction Coefficients -----------------'
 DO K=1,Orb%nAngmom
   WRITE(IUNIT,'(5X,A,I3,A,I3)') '* Angular block number',K,' with angular momentum', Orb%angmom(K)
   DO J=1,Orb%nContracted(K)
     WRITE(IUNIT,'(5X,I3,5F12.6,/(8X,5F12.6))') J,(Orb%CC(I,J,K),I=1,Orb%nPrimitives)
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
 PQ%ORIGO(1) = INPUT%ORIGO(1)
 PQ%ORIGO(2) = INPUT%ORIGO(2)
 PQ%ORIGO(3) = INPUT%ORIGO(3)
 PQ%P%p            => P
 PQ%Q%p            => Q
!simen Refine!
 PQ%nPrimitives    = P%nPrimitives   * Q%nPrimitives
 PQ%maxContracted  = P%maxContracted * Q%maxContracted
 PQ%minAngmom      = P%minAngmom     + Q%minAngmom
 PQ%maxAngmom      = P%maxAngmom     + Q%maxAngmom
 PQ%startAngmom    = P%startAngmom   + Q%startAngmom 
 PQ%endAngmom      = P%endAngmom     + Q%endAngmom 
 PQ%sameOD         = INPUT%sameODs .AND. (ILHS .EQ. IRHS)
 PQ%nAngmom        = P%nAngmom       * Q%nAngmom 
 IF (PQ%sameOD) THEN
  PQ%nAngmom = P%nAngmom*(P%nAngmom+1)/2
  PQ%nPrimitives = P%nPrimitives*(P%nPrimitives+1)/2
ENDIF
!Dimensions according to nPrimitives
Nullify(PQ%distance)
Nullify(PQ%squaredDistance)
Nullify(PQ%exponents)
Nullify(PQ%reducedExponents)
Nullify(PQ%integralPrefactor)
Nullify(PQ%iprimP)
Nullify(PQ%iprimQ)
Allocate(PQ%distance(PQ%nPrimitives,3))
Allocate(PQ%squaredDistance(PQ%nPrimitives))
Allocate(PQ%exponents(PQ%nPrimitives))
Allocate(PQ%reducedExponents(PQ%nPrimitives))
Allocate(PQ%integralPrefactor(PQ%nPrimitives))
Allocate(PQ%iprimP(PQ%nPrimitives))
Allocate(PQ%iprimQ(PQ%nPrimitives))
coulomb = (PQ%Operator.EQ.'Coulomb').OR.(PQ%Operator(1:3).EQ.'Erf')
nucrep = PQ%Operator(1:6).EQ.'Nucrep'

CALL SetIntegrand(PQ%iprimP,PQ%iprimQ,PQ%distance,PQ%squaredDistance,PQ%exponents,&
     &  PQ%reducedExponents,coulomb,nucrep,PQ%integralPrefactor,PQ%nPrimitives,PQ%sameOD,&
!     &  PQ%reducedExponents,PQ%operator,PQ%integralPrefactor,PQ%nPrimitives,PQ%sameOD,&
     &  P%distance12,P%squaredDistance,Q%distance12,Q%squaredDistance,&
     &  P%exponents,Q%exponents,P%center,Q%center,P%nPrimitives,Q%nPrimitives,&
     &  P%type .EQ. 'Empty',Q%type .EQ. 'Empty',IUNIT)

IF(INPUT%ATTFACTOR) CALL DSCAL(PQ%nPrimitives,1.d0+INPUT%ATTOMEGA/PI,PQ%integralPrefactor,1)

 IF (IPRINT.GT.15) THEN
    CALL PRINT_OVERLAP(P,IUNIT,IPRINT)
    CALL PRINT_OVERLAP(Q,IUNIT,IPRINT)
    CALL PRINT_Integrand(PQ,IUNIT)
 ENDIF
END SUBROUTINE Build_Integrand

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
    & PQ%distance(i,1),PQ%distance(i,2),PQ%distance(i,3),&
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
real(realk),ALLOCATABLE :: SJ000(:,:),SJ0002(:,:)
INTEGER                 :: LUPRI,SUM,J,K,T,U,V,TUV,IOFF
INTEGER                 :: nPrim,IPRINT,ntuv,L,I
INTEGER,ALLOCATABLE     :: TUVindex(:,:,:)
INTEGER                 :: zeroX,zeroY,zeroZ,Jmax,Jstart,IPQXYZ
real(realk)             :: X0,Y0,Z0

NPrim=PQ%nPrimitives
JMAX=PQ%endAngmom
Jstart=PQ%startAngmom
ALLOCATE(SJ000(nPrim,0:JMAX))

CALL MAXDISTANCE(X0,PQ%distance(:,1),nPrim)
CALL MAXDISTANCE(Y0,PQ%distance(:,2),nPrim)
CALL MAXDISTANCE(Z0,PQ%distance(:,3),nPrim)
zeroX=1
zeroY=1
zeroZ=1
IF(X0 .LE. 1.00D-15) zeroX=2
IF(Y0 .LE. 1.00D-15) zeroY=2
IF(Z0 .LE. 1.00D-15) zeroZ=2
IPQXYZ=0
IF(zeroX+zeroY+zeroZ .EQ. 6) IPQXYZ=7
SELECT CASE(OPERATORLABEL)
CASE ('Overlap')
   CALL BUILD_SJ000(SJ000,PQ,nPrim,JMAX,LUPRI,IPRINT,IPQXYZ)
CASE ('Coulomb')
   CALL BUILD_RJ000(SJ000,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,IPQXYZ)
CASE ('Erfc   ')
   ALLOCATE(SJ0002(nPrim,0:JMAX))
   CALL BUILD_RJ000(SJ000,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,IPQXYZ)
   CALL BUILD_ERF_RJ000(SJ0002,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,INPUT%ATTomega,IPQXYZ)
!   CALL DSCAL(NPrim*(JMAX+1),INPUT%ATTalpha,SJ000,1)
   CALL DAXPY(NPrim*(JMAX+1),INPUT%ATTbeta,SJ0002,1,SJ000,1)
   DEALLOCATE(SJ0002)
CASE ('Erf    ')
   CALL BUILD_ERF_RJ000(SJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,INPUT%ATTomega,IPQXYZ)
CASE ('Nucrep')
   CALL BUILD_NUCLEAR_RJ000(SJ000,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,IPQXYZ)
CASE DEFAULT
  WRITE(LUPRI,'(1X,A,A16)') 'Programming error! Not a case in GET_WTUV:',OPERATORLABEL
END SELECT

ntuv=(JMAX+1)*(JMAX+2)*(JMAX+3)/6
INTEGRAL%nTUV=ntuv

IF (JMAX .EQ. 0) THEN
   !Special case JMAX = 0
   NULLIFY(INTEGRAL%Wtuv)
   ALLOCATE(INTEGRAL%Wtuv(Nprim,nTUV))
   CALL DCOPY(nPrim,SJ000,1,Integral%Wtuv,1)
ELSE
   IF(Jstart .EQ. 0)THEN
      ALLOCATE(INTTEMP(nPrim,nTUV))
      NULLIFY(INTEGRAL%Wtuv)
      ALLOCATE(INTEGRAL%Wtuv(Nprim,nTUV))
      IF( zeroX+zeroY+zeroZ .NE. 3) CALL DZERO(INTEGRAL%Wtuv,nPrim*nTUV)
      DO J=1,JMAX
         IF (MOD(JMAX-J,2).EQ.0) THEN
            CALL Recurrence(INTEGRAL%Wtuv,INTTEMP,J,SJ000,&
                 & PQ%distance(:,1),PQ%distance(:,2),&
                 & PQ%distance(:,3),JMAX,NPrim,nTUV,Integral%TUV%TUVindex,&
                 & zeroX, zeroY, zeroZ)
         ELSE
            CALL Recurrence(INTTEMP,INTEGRAL%Wtuv,J,SJ000,&
                 & PQ%distance(:,1),PQ%distance(:,2),&
                 & PQ%distance(:,3),JMAX,NPrim,nTUV,Integral%TUV%TUVindex,&
                 & zeroX, zeroY, zeroZ)
         ENDIF
      ENDDO
   ELSE
      ALLOCATE(INTTEMP(nPrim,nTUV))
      ALLOCATE(INTTEMP2(nPrim,nTUV))
      IF( zeroX+zeroY+zeroZ .NE. 3)CALL DZERO(INTTEMP2,nPrim*nTUV)
      DO J=1,JMAX
         IF (MOD(JMAX-J,2).EQ.0) THEN
            CALL Recurrence(INTTEMP2,INTTEMP,J,SJ000,&
                 & PQ%distance(:,1),PQ%distance(:,2),&
                 & PQ%distance(:,3),JMAX,NPrim,nTUV,Integral%TUV%TUVindex,&
                 & zeroX, zeroY, zeroZ)
         ELSE
            CALL Recurrence(INTTEMP,INTTEMP2,J,SJ000,&
                 & PQ%distance(:,1),PQ%distance(:,2),&
                 & PQ%distance(:,3),JMAX,NPrim,nTUV,Integral%TUV%TUVindex,&
                 & zeroX, zeroY, zeroZ)
         END IF
      ENDDO
   ENDIF
ENDIF


IF (Jstart .NE. 0 .AND. JMAX .NE. 0) THEN
   IF (IPRINT .GE. 10) THEN
   CALL LSHEADER(LUPRI,'Output from Recurrence')
   WRITE (LUPRI,'(2X,A,I10)') 'NPrim  ', NPrim
   WRITE (LUPRI,'(2X,A,I10)') 'NTUV  ', nTUV
      IF (IPRINT .GE. 20) THEN
      CALL LSHEADER(LUPRI,'Hermite integrals S(t,u,v)')
         DO J = 0, JMAX
            DO T = J,0,-1
               DO U = J-T,0,-1
                  V=J-T-U
                  WRITE (LUPRI,'(2X,3(A,I1),A,2X,5F12.8/,(12X,5F12.8))')&
                       & 'W(',T,',',U,',',V,')', &
                       &(INTTEMP2(I,Integral%TUV%TUVINDEX(T,U,V)),I=1,NPrim)
                  WRITE (LUPRI,*)' '
               ENDDO
            ENDDO
         ENDDO
      END IF
   END IF

   IOFF= Jstart*(Jstart+1)*(Jstart+2)/6
   NULLIFY(INTEGRAL%Wtuv)
   ALLOCATE(INTEGRAL%Wtuv(Nprim,nTUV-IOFF))
   DO J = Jstart, JMAX
      DO T = J,0,-1
         DO U = J-T,0,-1
            V=J-T-U
            DO I=1,nPrim
               INTEGRAL%Wtuv(I,INTEGRAL%TUV%TUVINDEX(T,U,V)-IOFF)=INTTEMP2(I,INTEGRAL%TUV%TUVindex(T,U,V))
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(INTTEMP2) 
ENDIF


!    Print section
!    =============

IF (IPRINT .GE. 10) THEN
   IOFF= Jstart*(Jstart+1)*(Jstart+2)/6
   CALL LSHEADER(LUPRI,'Output from ERITUV')
   WRITE (LUPRI,'(2X,A,I10)') 'MAX angmom ', JMAX
   WRITE (LUPRI,'(2X,A,I10)') 'NPrim  ', NPrim
   WRITE (LUPRI,'(2X,A,I10)') 'NTUV  ', nTUV
   IF (IPRINT .GE. 20) THEN
      CALL LSHEADER(LUPRI,'Hermite integrals S(t,u,v)')
      DO J = Jstart, JMAX
         DO T = J,0,-1
            DO U = J-T,0,-1
               V=J-T-U
               WRITE (LUPRI,'(2X,A2,I1,A1,I1,A1,I1,A1,2X,5F12.8/,(12X,5F12.8))')&
                    & 'W(',T,',',U,',',V,')', &
                    &(INTEGRAL%Wtuv(I,INTEGRAL%TUV%TUVINDEX(T,U,V)-IOFF),I=1,NPrim)
               WRITE (LUPRI,*) ' '
            ENDDO
         ENDDO
      ENDDO
   END IF
END IF

DEALLOCATE(SJ000)
IF (JMAX .NE. 0) THEN
   DEALLOCATE(INTTEMP)
ENDIF

END SUBROUTINE GET_WTUV

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

SUBROUTINE BUILD_SJ000(SJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,IPQXYZ)
IMPLICIT NONE
TYPE(Integrand)        :: PQ
!real(realk)            :: SJ000(nPrim,PQ%startangmom:PQ%maxAngmom)
real(realk)            :: SJ000(nPrim,0:Jmax)
INTEGER                :: nPrim,J,I,LUPRI,IPRINT,Jmax,IPQXYZ
IF(IPQXYZ .EQ. 7)THEN !Same center -> squaredDistance=0
   CALL BUILD_SJ000_ZERO(SJ000,Nprim,Jmax,PQ%reducedExponents&
        &,PQ%exponents)
ELSE
   CALL BUILD_SJ000_OP(SJ000,Nprim,Jmax,PQ%reducedExponents&
        &,PQ%exponents,PQ%squaredDistance)
ENDIF

IF(IPRINT .GT. 25) THEN
   DO J=0,Jmax
      DO I=1,NPrim
         WRITE(LUPRI,'(2X,A5,I4,A1,I2,A2,F16.9)')'SJ000(',I,',',J,')=',SJ000(I,J)
      ENDDO
   ENDDO
ENDIF
   
END SUBROUTINE BUILD_SJ000

SUBROUTINE BUILD_SJ000_OP(SJ000,nPrim,Jmax,Alpha,P,R2)
IMPLICIT NONE
real(realk)            :: SJ000(nPrim,0:Jmax)
Real(realk), parameter :: PI=3.14159265358979323846D0
Real(realk), parameter :: OneHalf=1.5d0, Two = 2.0d0
INTEGER                :: nPrim,J,I,Jmax
REAL(REALK)            :: ALPHA(nPrim),P(nPrim),R2(nprim)
DO J=0,Jmax
   DO I=1,NPrim
      SJ000(I,J)=((-Two*Alpha(I))**J)&
           &*((PI/P(I))**OneHalf)*EXP(-Alpha(I)*R2(I))
   ENDDO
ENDDO
END SUBROUTINE BUILD_SJ000_OP

SUBROUTINE BUILD_SJ000_ZERO(SJ000,nPrim,Jmax,Alpha,P)
IMPLICIT NONE
real(realk)            :: SJ000(nPrim,0:Jmax)
Real(realk), parameter :: PI=3.14159265358979323846D0
Real(realk), parameter :: OneHalf=1.5d0, Two = 2.0d0
INTEGER                :: nPrim,J,I,Jmax
REAL(REALK)            :: ALPHA(nPrim),P(nPrim)
DO J=0,Jmax
   DO I=1,NPrim
      SJ000(I,J)=((-Two*Alpha(I))**J)*((PI/P(I))**OneHalf)
   ENDDO
ENDDO
 END SUBROUTINE BUILD_SJ000_ZERO

SUBROUTINE Contract_P(INTEGRAL,PQ,input,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: INTEGRAL
TYPE(IntegralInput):: INPUT
TYPE(Integrand)    :: PQ
Integer            :: LUPRI,IPRINT 
!
Integer  :: iAngmom,startMomQ,iOrbital,ideriv,nderiv
!
nderiv=INPUT%NderivP
NULLIFY(Integral%PQ)
ALLOCATE(Integral%PQ(PQ%Q%p%totOrbitals,PQ%P%p%totOrbitals))

iOrbital = 1
!
!Loop over angular contributions sharing the set of primitive functions
!DO iDerivP=1,nderivP NOT YET IMPLEMENTED
   DO iAngmom=1,PQ%P%p%nAngmom
      ! Distribute Hermite 2-el integrals WTUV to <tuv1|w|tuv2>
      CALL DistributeHermiteP(Integral,PQ,iAngmom,LUPRI,IPRINT)
      CALL ContractEcoeff(Integral,PQ%P%p,iAngmom,LUPRI,IPRINT)
      CALL ContractBasis(Integral,PQ%P%p,iAngmom,LUPRI,IPRINT)
      CALL SphericalTransform(Integral,PQ%P%p,iAngmom,LUPRI,IPRINT)
      CALL AddToPQ(Integral,PQ,iAngmom,iOrbital,LUPRI,IPRINT)
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
Integer :: nQ,nAng,nCont,endOrb

nQ     = PQ%Q%p%totOrbitals
nAng   = PQ%P%p%nOrbComp(iAngmom)
nCont  = PQ%P%p%nContracted(iAngmom)
endOrb = iOrbital + nQ*nCont*nAng-1

CALL AddToPQ1(Integral%PQ(:,iOrbital:endOrb),&
     &        Integral%IntegralIN,nQ,nAng,nCont,LUPRI,IPRINT)
DEALLOCATE(Integral%IntegralIN)
NULLIFY(Integral%IntegralIN)
END SUBROUTINE AddToPQ

SUBROUTINE AddToPQ1(PQ,AddPQ,nQ,nAng,nCont,LUPRI,IPRINT)
implicit none
Real(realk) :: PQ(nQ,nAng,nCont),AddPQ(nAng,nQ,nCont)
Integer     :: nQ,nAng,nCont,LUPRI,IPRINT
!
Integer     :: iQ,iAng,iCont
!
DO iCont=1,nCont
  call transposition(PQ(1,1,iCont),AddPQ(1,1,iCont),nQ,nAng)
ENDDO
CALL PrintPQ(PQ,nQ,nAng,nCont,LUPRI,IPRINT)
END SUBROUTINE AddToPQ1

SUBROUTINE PrintPQ(PQ,nQ,nAng,nCont,LUPRI,IPRINT)
implicit none
Real(realk) :: PQ(nQ,nAng,nCont)
Integer     :: nQ,nAng,nCont,LUPRI,IPRINT
!
Integer     :: iQ,iAng,iCont
!
IF (IPRINT.GT.5) THEN
  CALL LSHEADER(LUPRI,'PQ contribution')
  DO iCont=1,nCont
    DO iAng=1,nAng
      WRITE(LUPRI,'(5X,A,I3,A,I3)') 'iCont =',iCont,' iAng =',iAng
      WRITE(LUPRI,'(5X,F16.9)') (PQ(iQ,iAng,iCont),iQ=1,nQ)
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
Real(realk), parameter :: Half=0.5d0

!
startP = 0
IF (PQ%P%p%type.EQ.'Hermite-single') startP = PQ%P%p%angmom(iAngmom)

IF(PQ%Operator(1:7) .EQ. 'Kinetic')THEN
   fullSP = PQ%P%p%startAngmom
   endP   = PQ%P%p%angmom(iAngmom)
   nP     = PQ%P%p%nPrimitives
   ioffP  = startP*(startP+1)*(startP+2)/6
   fullOP = fullSP*(fullSP+1)*(fullSP+2)/6
   ntuvP  = (endP+1)*(endP+2)*(endP+3)/6 - ioffP
   nOrbQ  = PQ%Q%p%totOrbitals
  
   Integral%nAng  = ntuvP
   Integral%nPrim = np
   Integral%nOrb  = nOrbQ

   NULLIFY(Integral%IntegralIN)
   ALLOCATE(Integral%IntegralIN(ntuvP,nOrbQ,nP))
   
   DO iOrbQ=1,nOrbQ
      DO jP = startP,endP
         DO tP=jP,0,-1
            DO uP=jP-tP,0,-1
               vP=jP-tP-uP
               ituvP  = Integral%TUV%TUVindex(tP,uP,vP)-ioffP
               ifullP1 = Integral%TUV%TUVindex(tP+2,uP,vP)-fullOP
               ifullP2 = Integral%TUV%TUVindex(tP,uP+2,vP)-fullOP
               ifullP3 = Integral%TUV%TUVindex(tP,uP,vP+2)-fullOP             
               DO iPrimP=1,nP
                  Integral%IntegralIN(ituvP,iOrbQ,iPrimP) = &
                       &- HALF*Integral%TUVQ(iPrimP,ifullP1,iOrbQ) &
                       &- HALF*Integral%TUVQ(iPrimP,ifullP2,iOrbQ) &
                       &- HALF*Integral%TUVQ(iPrimP,ifullP3,iOrbQ)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ELSE
   fullSP = PQ%P%p%startAngmom
   endP   = PQ%P%p%angmom(iAngmom)
   nP     = PQ%P%p%nPrimitives
   ioffP  = startP*(startP+1)*(startP+2)/6
   fullOP = fullSP*(fullSP+1)*(fullSP+2)/6
   ntuvP  = (endP+1)*(endP+2)*(endP+3)/6 - ioffP
   nOrbQ  = PQ%Q%p%totOrbitals
   
   Integral%nAng  = ntuvP
   Integral%nPrim = np
   Integral%nOrb  = nOrbQ
   
   NULLIFY(Integral%IntegralIN)
   ALLOCATE(Integral%IntegralIN(ntuvP,nOrbQ,nP))
   
   DO iOrbQ=1,nOrbQ
      DO jP = startP,endP
         DO tP=jP,0,-1
            DO uP=jP-tP,0,-1
               vP=jP-tP-uP
               ituvP  = Integral%TUV%TUVindex(tP,uP,vP)-ioffP
               ifullP = Integral%TUV%TUVindex(tP,uP,vP)-fullOP
               DO iPrimP=1,nP
                  Integral%IntegralIN(ituvP,iOrbQ,iPrimP) = Integral%TUVQ(iPrimP,ifullP,iOrbQ)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF

IF (IPRINT.GT.20) THEN
   CALL PrintTensor(Integral%IntegralIN,'Primitive           ',&
        &ntuvP,nOrbQ,nP,Lupri,'iTUV  ','iOrb  ','Prim  ',3)
ENDIF

END SUBROUTINE DistributeHermiteP

SUBROUTINE Contract_Q(INTEGRAL,PQ,Input,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)  :: INTEGRAL
TYPE(Integrand)     :: PQ
TYPE(IntegralInput) :: Input
Integer             :: LUPRI,IPRINT
!
Integer  :: iAngmom,startMomQ,iOrbital,ideriv,nderiv
!
nderiv=Input%NderivQ
NULLIFY(Integral%TUVQ)
ALLOCATE(Integral%TUVQ(PQ%P%p%nPrimitives,PQ%P%p%nTUV,PQ%Q%p%totOrbitals))
!
iOrbital = 1
!
!Loop over angular contributions sharing the set of primitive functions
DO iDeriv=1,nderiv
   DO iAngmom=1,PQ%Q%p%nAngmom
      ! Distribute Hermite 2-el integrals WTUV to <tuv1|w|tuv2>
      CALL DistributeHermiteQ(Integral,PQ,iAngmom,iDeriv,LUPRI,IPRINT)
      CALL ContractEcoeff(Integral,PQ%Q%p,iAngmom,LUPRI,IPRINT)
      CALL ContractBasis(Integral,PQ%Q%p,iAngmom,LUPRI,IPRINT)
      CALL SphericalTransform(Integral,PQ%Q%p,iAngmom,LUPRI,IPRINT)
      CALL AddToTUVQ(Integral,PQ,iAngmom,iOrbital,LUPRI,IPRINT)
      iOrbital = iOrbital + PQ%Q%p%nOrbitals(iAngmom)
   ENDDO
ENDDO
!IF(INPUT%OPERATOR .EQ. 'mom 100' )CALL QUIT('TEST FINE - done contract_q ')
    
END SUBROUTINE Contract_Q

SUBROUTINE AddToTUVQ(Integral,PQ,iAngmom,iOrbital,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: iAngmom,iOrbital,LUPRI,IPRINT 
!
Integer :: nPrimP,nTUVP,nCompQ,nContQ,endOrb
!
nPrimP = PQ%P%p%nPrimitives
nTUVP  = PQ%P%p%nTUV
nCompQ = PQ%Q%p%nOrbComp(iAngmom)
nContQ = PQ%Q%p%nContracted(iAngmom)
endOrb = iOrbital + nCompQ*nContQ - 1

!CALL AddToTUVQ1(Reshape(Integral%TUVQ(:,:,iOrbital:endOrb),(/nPrimP,nTUVP,nCompQ,nContQ/)),&
CALL AddToTUVQ1(Integral%TUVQ(:,:,iOrbital:endOrb),&
     &          Integral%IntegralIN,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
DEALLOCATE(Integral%IntegralIN)
NULLIFY(Integral%IntegralIN)
END SUBROUTINE AddToTUVQ

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
ENDDO

CALL PrintTUVQ(TUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
END SUBROUTINE AddToTUVQ1

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

SUBROUTINE SphericalTransform(Integral,P,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,iAngmom
!
Real(realk),pointer :: Spherical(:,:)
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
  
  NULLIFY(Spherical)
  ALLOCATE(Spherical(ijk1*ijk2,lm1*lm2))
  NULLIFY(Integral%IntegralOUT)
  ALLOCATE(Integral%IntegralOUT(lm,Integral%nOrb,Integral%nPrim))

  CALL ContructSphericalTransformation(Spherical,P%orbital1%SPH_MAT(ang1+1)%p%elms,&
     &              P%orbital2%SPH_MAT(ang2+1)%p%elms,lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)

! Spherical transformation (with both centers beloning to one electron simultaneously)
  dim1 = Integral%nOrb*Integral%nPrim
  CALL DGEMM('T','N',lm,dim1,ijk,D1,Spherical,ijk,&
     &     Integral%IntegralIN,ijk,D0,Integral%IntegralOUT,lm)

  DEALLOCATE(Spherical)
  NULLIFY(Spherical)
  DEALLOCATE(Integral%IntegralIN)
  Integral%nAng = lm
  Integral%IntegralIN => Integral%IntegralOUT
  NULLIFY(Integral%IntegralOUT)
ENDIF
END SUBROUTINE SphericalTransform

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

SUBROUTINE ContractBasis(Integral,P,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,iAngmom
!
Integer             :: iA1,iA2,nCont,nC1,nC2,dim1
Real(realk),pointer :: CC(:,:,:)
Real(realk)         :: D1=1.0d0,D0=0.0d0
iA1 = P%indexAng1(iangmom)
iA2 = P%indexAng2(iangmom)
nC1 = P%orbital1%nContracted(iA1)
nC2 = P%orbital2%nContracted(iA2)
nCont = nC1*nC2

NULLIFY(Integral%IntegralOUT)
ALLOCATE(Integral%IntegralOUT(Integral%nAng,Integral%nOrb,nCont))
NULLIFY(CC)
ALLOCATE(CC(P%nPrimitives,nC1,nC2))

CALL ConstructContraction(CC,P,Integral%nPrim,nC1,nC2,iA1,iA2,LUPRI,IPRINT)

!Actual contraction (with both centers on Q simultaneously)
dim1 = Integral%nAng*Integral%nOrb

CALL DGEMM('N','N',dim1,nCont,Integral%nPrim,D1,Integral%IntegralIN,dim1,&
     &     CC,Integral%nPrim,D0,Integral%IntegralOUT,dim1)


IF (IPRINT.GT.50) THEN
   CALL PrintTensor(Integral%IntegralOUT,'Contracted          ',&
        &Integral%nAng,Integral%nOrb,nCont,Lupri,'iAngP ','iOrb  ','iCont ',3)
ENDIF

Integral%nPrim = nCont
DEALLOCATE(Integral%IntegralIN)
Integral%IntegralIN => Integral%IntegralOUT
NULLIFY(Integral%IntegralOUT)
DEALLOCATE(CC)
END SUBROUTINE ContractBasis

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

SUBROUTINE ConstructContraction(CC,P,nP,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
Implicit none
Type(Overlap) :: P
Integer       :: iA1,iA2,nP,nC1,nC2,LUPRI,IPRINT
Real(realk)   :: CC(nP,nC1,nC2)
!
Integer       :: iP,i1,i2,iC1,iC2

IF(P%sameAO)THEN
   DO iC2=1,nC2
      DO iC1=1,nC1
         DO iP=1,nP
            i1 = P%iprim1(iP)
            i2 = P%iprim2(iP)
            IF(i1 .NE. i2)THEN
            CC(iP,iC1,iC2)=P%orbital1%CC(i1,iC1,iA1)*P%orbital2%CC(i2,iC2,iA2)&
                 &        +P%orbital1%CC(i2,iC1,iA1)*P%orbital2%CC(i1,iC2,iA2)
            ELSE
            CC(iP,iC1,iC2)=P%orbital1%CC(i1,iC1,iA1)*P%orbital2%CC(i2,iC2,iA2)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ELSE
   DO iC2=1,nC2
      DO iC1=1,nC1
         DO iP=1,nP
            i1 = P%iprim1(iP)
            i2 = P%iprim2(iP)
            CC(iP,iC1,iC2)=P%orbital1%CC(i1,iC1,iA1)*P%orbital2%CC(i2,iC2,iA2) 
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

SUBROUTINE ContractEcoeff(Integral,P,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)  :: Integral
TYPE(Overlap)       :: P
Integer             :: LUPRI,IPRINT,iAngmom
Real(realk),pointer :: Ecoeffs(:,:,:)
Integer :: l1,l2,ijk1,ijk2,ijk,iPrimP
Real(realk)         :: D1=1.0d0,D0=0.0d0
Real(realk),external:: ddot
IF (P%type.EQ.'Hermite-single') THEN
  CALL SingleHermiteEcoeff(Integral%IntegralIn,P,iAngmom,Integral%nAng,Integral%nOrb,LUPRI,IPRINT)
ELSE
!  CALL QUIT('Only Hermite implemented in ContractEcoeff')
  l1 = P%orbital1%angmom(P%indexAng1(iangmom))
  l2 = P%orbital2%angmom(P%indexAng2(iangmom))
  ijk1 = (l1+1)*(l1+2)/2
  ijk2 = (l2+1)*(l2+2)/2
  ijk  = ijk1*ijk2
  NULLIFY(Integral%IntegralOUT)
  ALLOCATE(Integral%IntegralOUT(ijk,Integral%nOrb,Integral%nPrim))

IF (P%ETUVisSet) THEN
  CALL ContractEcoeff1(Integral%IntegralIN,Integral%IntegralOUT,&
       &P%ETUV(P%ETUVindex(iAngmom):P%ETUVindex(iAngmom)+Integral%nPrim*Integral%nAng*ijk-1), &
     &                 Integral%nOrb,Integral%nAng,ijk,Integral%nPrim)
ELSE
  NULLIFY(Ecoeffs)
  ALLOCATE(Ecoeffs(Integral%nAng,ijk,Integral%nPrim))
  CALL BuildEcoeffTensor(integral%TUV,P,Ecoeffs,Integral%nAng,ijk,Integral%nPrim,&
            &            iAngmom,LUPRI,IPRINT)
  CALL ContractEcoeff1(Integral%IntegralIN,Integral%IntegralOUT,Ecoeffs, &
     &                 Integral%nOrb,Integral%nAng,ijk,Integral%nPrim)
  DEALLOCATE(Ecoeffs)
ENDIF

  Integral%nAng = ijk
  DEALLOCATE(Integral%IntegralIN)
  Integral%IntegralIN => Integral%IntegralOUT
  NULLIFY(Integral%IntegralOUT)
ENDIF

IF (IPRINT.GT.50) THEN
   CALL PrintTensor(Integral%IntegralIN,'ContractEcoeff      ',&
        &Integral%nAng,Integral%nOrb,Integral%nPrim,Lupri,&
        &'ijk   ','iOrb  ','Prim  ',3)
ENDIF

END SUBROUTINE ContractEcoeff

SUBROUTINE ContractEcoeff1(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nAng,ijk,nPrim)
implicit none
Integer     :: nOrb,nAng,ijk,nPrim
Real(realk) :: IntegralIN(nAng,nOrb,nPrim)
Real(realk) :: IntegralOUT(ijk,nOrb,nPrim)
Real(realk) :: Ecoeffs(nAng,ijk,nPrim)
!
Integer :: iPrimP
Real(realk)         :: D1=1.0d0,D0=0.0d0

DO iPrimP = 1,nPrim
   CALL DGEMM('T','N',ijk,nOrb,nAng,D1,Ecoeffs(1,1,iPrimP),nAng,&
   &          IntegralIN(1,1,iPrimP),nAng,D0,IntegralOUT(1,1,iPrimP),ijk)
ENDDO

END SUBROUTINE ContractEcoeff1

SUBROUTINE ContractFTUV(Integral,P,Q,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P,Q
Integer            :: LUPRI,IPRINT
IF (P%type.EQ.'Hermite-single') THEN
  CALL QUIT('Jengine not implemented with hermite-single yet')
!  CALL SingleHermiteEcoeff(Integral%IntegralIn,P,iAngmom,Integral%nAng,Integral%nOrb,LUPRI,IPRINT)
ELSE
  CALL ContractFTUV1(Integral%tuvQ,Integral%WTUV,Q%FTUV,Integral%TUV%TUVindex,&
     &           Integral%TUV%Tindex,Integral%TUV%Uindex,Integral%TUV%Vindex,&
     &               P%nTUV,Q%nTUV,P%nPrimitives,Q%nPrimitives)
ENDIF

Integral%nAng  = 1
Integral%nPrim = 1
Integral%nOrb  = P%nPrimitives*P%nTUV
CALL PrintTUVQ(Integral%tuvQ,P%nTUV,P%nPrimitives,1,1,LUPRI,IPRINT)

END SUBROUTINE ContractFTUV

SUBROUTINE ContractFTUV1(TUV,WTUV,FTUV,TUVindex,Tindex,Uindex,Vindex,&
     &                   ntuvP,ntuvQ,nPrimP,nPrimQ)
implicit none
Real(realk),pointer :: TUV(:,:,:),WTUV(:,:),FTUV(:,:,:)
Integer,pointer     :: TUVindex(:,:,:),Tindex(:),Uindex(:),Vindex(:)
Integer             :: ntuvP,ntuvQ,nPrimP,nPrimQ
!
Integer     :: idmat,tuvP,tP,uP,vP,tuvQ,tQ,uQ,vQ,jQ,tuvPQ
Integer     :: iPrimP,iPrimQ,iPrimPQ
Real(realk) :: D1=1.0d0,DM1=-1.0d0,signQ,fsum

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
    iPrimPQ = 1
    DO iPrimP=1,nPrimP
      fsum = WTUV(iPrimPQ,tuvPQ)*FTUV(1,tuvQ,idmat)
      iPrimPQ = iPrimPQ+1
      DO iPrimQ=2,nPrimQ
        fsum = fsum + WTUV(iPrimPQ,tuvPQ)*FTUV(iPrimQ,tuvQ,idmat)
        iPrimPQ = iPrimPQ+1
      ENDDO
      TUV(iPrimP,tuvP,idmat) = TUV(iPrimP,tuvP,idmat) + fsum*signQ
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE ContractFTUV1

SUBROUTINE BuildEcoeffTensor(TUV,Q,Ecoeffs,x,y,z,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: Q
Integer            :: iAngmom,LUPRI,IPRINT
Integer            :: ijk
!Allocatable for ECOEFFS - Dimensions are nPrimitives, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
real(realk),allocatable  :: ETIJ(:,:,:,:,:) 
integer             :: startQ,endQ,nQ,ioffQ,ntuvQ,ijk1,ijk2,jQ,Q1,Q2,tQ,uQ,vQ
integer             :: ituvQ,iQ1,jQ1,kQ1,iQ2,jQ2,kQ2,iPrimQ
integer             :: iorb1,iorb2,iang,iprim,x,y,z,l1,l2
Real(realk)         :: Ecoeffs(x,y,z) !nAng,ijk1*ijk2,nprim

l1 = Q%orbital1%angmom(Q%indexAng1(iangmom))
l2 = Q%orbital2%angmom(Q%indexAng2(iangmom))

ALLOCATE(ETIJ(Q%nPrimitives,0:l1+l2,0:l1,0:l2,3))

CALL GET_ECOEFF(ETIJ,Q%nPrimitives,l1+l2,l1,l2,Q,LUPRI,IPRINT)

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
                              ituvQ=TUV%TUVindex(tQ,uQ,vQ)
                              DO iPrimQ=1,Q%nPrimitives
                                 Ecoeffs(ituvQ,ijk,iPrimQ) = &
                                      &ETIJ(iPrimQ,tQ,iQ1,iQ2,1)&
                                      &*ETIJ(iPrimQ,uQ,jQ1,jQ2,2)&
                                      &*ETIJ(iPrimQ,vQ,kQ1,kQ2,3)&
                                      &*Q%preExpFac(iPrimQ)
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

END SUBROUTINE BuildEcoeffTensor

SUBROUTINE SingleHermiteEcoeff(HerInt,P,iAngmom,nAng,nOrb,LUPRI,IPRINT)
implicit none
Real(realk)   :: HerInt(:,:,:)
TYPE(Overlap) :: P
Integer       :: iAngmom,nAng,nOrb,LUPRI,IPRINT
!
Integer     :: nPrim1,nPrim2,ang1,ang2
Integer     :: j1,j2,i1,i2
!
i1 = P%indexAng1(iAngmom)
i2 = P%indexAng2(iAngmom)
j1 = P%orbital1%angmom(i1)
j2 = P%orbital2%angmom(i2)
nPrim1 = P%orbital1%nPrimitives
nPrim2 = P%orbital2%nPrimitives
CALL SingleHermiteEcoeff1(HerInt,P,j1,j2,nPrim1,nPrim2,nAng,nOrb,LUPRI,IPRINT)
END SUBROUTINE SingleHermiteEcoeff

SUBROUTINE SingleHermiteEcoeff1(HerInt,P,j1,j2,nPrim1,nPrim2,nAng,nOrb,LUPRI,IPRINT)
implicit none
Real(realk)   :: HerInt(nAng*nOrb,nPrim1,nPrim2)
TYPE(Overlap) :: P
Integer       :: nAng,nOrb,LUPRI,IPRINT
Integer       :: nPrim1,nPrim2,j1,j2
!
Real(realk)   :: pref2,pref12,D1=1.0d0,D2=2.0d0
Integer       :: iPrim1,iPrim2,iAngOrb
DO iPrim2=1,nPrim2
  pref2=D1/(D2*P%orbital2%exponents(iPrim2))**j2
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
Integer                 :: iPrimPQ,iP,iQ,iTUV
Integer                 :: ntuvP,ntuvQ,startP,endP,startQ,endQ,endPQ,ntuvPQ
Real(realk)             :: signQ,DM1 = -1.0d0,sign
Real(realk),pointer     :: intermediate1(:,:,:),intermediate2(:,:)
startP = PQ%P%p%startAngmom
endP   = PQ%P%p%endAngmom
startQ = 0
IF (PQ%Q%p%type.EQ.'Hermite-single') startQ = PQ%Q%p%angmom(iAngmom)
endQ   = PQ%Q%p%angmom(iAngmom)
nP     = PQ%P%p%nPrimitives
nQ     = PQ%Q%p%nPrimitives
ioffP  = startP*(startP+1)*(startP+2)/6
ioffQ  = startQ*(startQ+1)*(startQ+2)/6
sPQ    = startP+PQ%Q%p%startAngmom
ioffPQ = sPQ*(sPQ+1)*(sPQ+2)/6
endPQ  = endQ + PQ%P%p%endAngmom

ntuvP  = PQ%P%p%nTUV
ntuvQ  = (endQ+1)*(endQ+2)*(endQ+3)/6-ioffQ
ntuvPQ = (endPQ+1)*(endPQ+2)*(endPQ+3)/6-ioffPQ
Integral%nOrb  = ntuvP*nP
Integral%nAng  = ntuvQ
Integral%nPrim = nQ
NULLIFY(Integral%tuvTUV)
ALLOCATE(Integral%tuvTUV(ntuvQ,ntuvP,nP,nQ))
!CALL DZERO(Integral%tuvTUV,ntuvQ*ntuvP*nP*nQ)
DO jQ = startQ,endQ
   signQ = DM1**jQ
   DO tQ=jQ,0,-1
      DO uQ=jQ-tQ,0,-1
         vQ=jQ-tQ-uQ
         ituvQ=Integral%TUV%TUVindex(tQ,uQ,vQ)-ioffQ
         DO jP = startP,endP
            sign = DM1**(jP-jQ)
            DO tP=jP,0,-1
               DO uP=jP-tP,0,-1
                  vP=jP-tP-uP
                  ituvP=Integral%TUV%TUVindex(tP,uP,vP)-ioffP
                  ituvPQ=Integral%TUV%TUVindex(tP+tQ,uP+uQ,vP+vQ)-ioffPQ+(ideriv-1)*ntuvPQ
                  DO iPrimPQ=1,PQ%nPrimitives
                     Integral%tuvTUV(ituvQ,ituvP,PQ%iPrimP(iPrimPQ),PQ%iPrimQ(iPrimPQ)) &
                          &  = signQ*Integral%WTUV(iPrimPQ,ituvPQ)
                     IF(PQ%sameOD) &
                          & Integral%tuvTUV(ituvQ,ituvP,PQ%iPrimQ(iPrimPQ),PQ%iPrimP(iPrimPQ))&
                          &  = sign*signQ*Integral%WTUV(iPrimPQ,ituvPQ)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
CALL PrintHermitePQ(Integral,np,nq,startP,endP,startQ,endQ,LUPRI,IPRINT)

NULLIFY(Integral%IntegralIN)
ALLOCATE(Integral%IntegralIN(ntuvQ,ntuvP*nP,nQ))
CALL DCOPY(ntuvP*nP*nQ*ntuvQ,Integral%tuvTUV,1,Integral%IntegralIN,1)
DEALLOCATE(Integral%tuvTUV)
NULLIFY(Integral%tuvTUV)

END SUBROUTINE DistributeHermiteQ

SUBROUTINE PrintHermitePQ(Integral,np,nq,startP,endP,startQ,endQ,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
Integer            :: np,nq,startP,endP,startQ,endQ
Integer            :: LUPRI,IPRINT
!
Integer            :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ioffP,ioffQ,iPrimP,iPrimQ
!
IF (IPRINT.GT.50) THEN
  WRITE(LUPRI,'(3X,A)') '***************************************************************'
  WRITE(LUPRI,'(3X,A)') '***                         HermitePQ'
  WRITE(LUPRI,'(3X,A)') '***************************************************************'
  ioffP  = startP*(startP+1)*(startP+2)/6
  ioffQ  = startQ*(startQ+1)*(startQ+2)/6
  DO jQ = startQ,endQ
    DO tQ=jQ,0,-1
      DO uQ=jQ-tQ,0,-1
        vQ=jQ-tQ-uQ
        ituvQ=Integral%TUV%TUVindex(tQ,uQ,vQ)-ioffQ
!       WRITE(LUPRI,'(3X,A,3I2)') 'tQ,uQ,vQ',tQ,uQ,vQ
        DO jP = startP,endP
          DO tP=jP,0,-1
            DO uP=jP-tP,0,-1
            vP=jP-tP-uP
            ituvP=Integral%TUV%TUVindex(tP,uP,vP)-ioffP
            WRITE(LUPRI,'(5X,A,I1,A,I1,A,I1,A,I1,A,I1,A,I1,A)') &
      &            'W(',tP,',',uP,',',vP,'|',tQ,',',uQ,',',vQ,') ='
            DO iPrimP=1,np
               WRITE(LUPRI,'(17X,I3,5F15.4/,(20X,5F15.4))') iPrimP,&
      &           (Integral%tuvTUV(ituvQ,ituvP,iPrimP,iPrimQ),iPrimQ=1,nq)
            ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintHermitePQ

SUBROUTINE Build_PRECALCULATED_SPHMAT(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)
IMPLICIT NONE
TYPE(ODITEM)     :: OD_LHS,OD_RHS
TYPE(TUVitem)    :: SharedTUV
INTEGER          :: LUPRI,IPRINT
!
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

END SUBROUTINE BUILD_PRECALCULATED_SPHMAT

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
SUBROUTINE GET_ECOEFF(ETIJ,nprim,maxang12,maxang1,maxang2,P,LUPRI,IPRINT)
implicit none
REAL(REALK),PARAMETER    :: D1 = 1.0D0, D2 = 2.0D0, DHALF = 0.5D0
integer                  :: IPRINT, LUPRI, I, J, K, CARTDIR
integer                  :: nprim,maxang12,maxang1,maxang2,l
TYPE(Overlap)            :: P
real(realk),allocatable  :: PINV(:), APINV(:), BPINV(:), HPINV(:), PA(:,:), PB(:,:)
real(realk),allocatable  :: TWOA(:), TWOB(:)
!Allocatable for ECOEFFS - Dimensions are No. Pairs, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
!real(realk),allocatable  :: ETIJ(:,:,:,:,:) 
real(realk)  :: ETIJ(nprim,0:maxang12,0:maxang1,0:maxang2,3) 

IF (IPRINT.GT.50) CALL PRN_ECINFO(P,LUPRI)

ALLOCATE(PINV(P%nPrimitives))
ALLOCATE(APINV(P%nPrimitives))
ALLOCATE(BPINV(P%nPrimitives))
ALLOCATE(HPINV(P%nPrimitives))
ALLOCATE(PA(P%nPrimitives,3))
ALLOCATE(PB(P%nPrimitives,3))
ALLOCATE(TWOA(P%nPrimitives))
ALLOCATE(TWOB(P%nPrimitives))
!CAMT WHERE DO WE DEALLOCATE THIS ?
!ALLOCATE(ETIJ(P%nPrimitives,0:P%orbital1%maxAngmom+P%orbital2%maxAngmom,&
!&             0:P%orbital1%maxAngmom,0:P%orbital2%maxAngmom,3))

!CAMT SETUP REQUIRED EXPONENT RELATED QUANTITIES
DO K=1,P%nPrimitives
   I=P%iprim1(K)
   J=P%iprim2(K)
   PINV(K)     =  D1 / (P%orbital1%exponents(I)& 
        &                           + P%orbital2%exponents(J))    ! 1/p
   APINV(K)    =  P%orbital1%exponents(I) * PINV(K)      ! a/p
   BPINV(K)    =  P%orbital2%exponents(J) * PINV(K)      ! b/p
   HPINV(K)    =  DHALF*PINV(K)                     ! 1/2p
   TWOA(K)     =  D2*P%orbital1%exponents(I)             ! 2a
   TWOB(K)     =  D2*P%orbital2%exponents(J)             ! 2b

!CAMT SETUP REQUIRED DISTANCES
   PA(K,1) = -BPINV(K)*P%distance12(1)
   PA(K,2) = -BPINV(K)*P%distance12(2)
   PA(K,3) = -BPINV(K)*P%distance12(3)
   PB(K,1) =  APINV(K)*P%distance12(1)
   PB(K,2) =  APINV(K)*P%distance12(2)
   PB(K,3) =  APINV(K)*P%distance12(3)
END DO

!AMT THEN CALL ROUTINES FOR CARTESIAN OR HERMITE ECOEFFS 
!AMT LOOP OVER CARTESIAN DIRECTONS (1,2,3 -> X,Y,Z)
DO CARTDIR = 1,3
      IF (P%type(1:4) .EQ. 'Herm') THEN
           CALL HERM_ECOEFFS(ETIJ(1,0,0,0,CARTDIR),P%nPrimitives,maxAng1,&
          &             maxAng2,PA(:,CARTDIR),PB(:,CARTDIR),PINV,APINV,BPINV,HPINV,       &
          &             CARTDIR,TWOA,TWOB,LUPRI)
      ELSE IF (P%type(1:4) .EQ. 'Cart') THEN
           CALL CART_ECOEFFS(ETIJ(1,0,0,0,CARTDIR),P%nPrimitives,maxAng1,&
          &             maxAng2,PA(:,CARTDIR),PB(:,CARTDIR),PINV,APINV,BPINV,HPINV,       &
          &             CARTDIR,LUPRI)
      ELSE IF (P%type(1:4) .EQ. 'Empt') THEN
           DO i=1,P%nPrimitives
              ETIJ(i,0,0,0,CARTDIR) = D1
           ENDDO
      ELSE
           WRITE(LUPRI,*)'ERROR : UNRECOGNISED TYPE !!!', P%type
      ENDIF
      IF (IPRINT.GT.20) CALL PRN_ECOEFFS(ETIJ,P%nPrimitives,maxAng1,&
                     &            maxAng2,CARTDIR,LUPRI)
ENDDO
! DEALLOCATE MEMORY
DEALLOCATE(PINV)
DEALLOCATE(APINV)
DEALLOCATE(BPINV)
DEALLOCATE(HPINV)
DEALLOCATE(PA)
DEALLOCATE(PB)
DEALLOCATE(TWOA)
DEALLOCATE(TWOB)
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
WRITE(LUPRI,*)'X',P%distance12(1)
WRITE(LUPRI,*)'Y',P%distance12(2)
WRITE(LUPRI,*)'Z',P%distance12(3)
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
                        & + T1*ETIJ(K,T+1,I,J-1)!-TJM*ETIJ(K,  T,I,J-2)/TWOB(K)
                END DO
             END DO
          ELSE ! J.GE.2
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

SUBROUTINE BUILD_RJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,IPQXYZ)
IMPLICIT NONE
INTEGER         :: nPrim,LUPRI,IPRINT,Jmax,IPQXYZ
REAL(REALK)     :: RJ000(nPrim,0:Jmax)
TYPE(Integrand) :: PQ
TYPE(integralitem) :: integral

IF(IPQXYZ .EQ. 7)THEN !Same center -> squaredDistance=0
   CALL BUILD_RJ000_ZERO(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
        &PQ%reducedExponents,PQ%integralPrefactor)
ELSE
   CALL BUILD_RJ000_OP(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
        &PQ%reducedExponents,PQ%squaredDistance,&
        &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW)
ENDIF
   
END SUBROUTINE BUILD_RJ000

SUBROUTINE BUILD_NUCLEAR_RJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,IPQXYZ)
IMPLICIT NONE
INTEGER         :: nPrim,LUPRI,IPRINT,Jmax,IPQXYZ
REAL(REALK)     :: RJ000(nPrim,0:Jmax)
TYPE(Integrand) :: PQ
TYPE(integralitem) :: integral

IF(IPQXYZ .EQ. 7)THEN !Same center -> squaredDistance=0
   CALL BUILD_RJ000_ZERO(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
        &PQ%Exponents,PQ%integralPrefactor)
ELSE
   CALL BUILD_RJ000_OP(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
        &PQ%Exponents,PQ%squaredDistance,&
        &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW)
ENDIF

END SUBROUTINE BUILD_NUCLEAR_RJ000

SUBROUTINE BUILD_RJ000_OP(RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW)
IMPLICIT NONE
INTEGER         :: nPrim,Jmax,Lupri,Iprint,nTABFJW
REAL(REALK)     :: RJ000(nPrim,0:Jmax)
REAL(REALK)     :: alpha(nPrim),R2(nprim),Prefactor(nprim),TABFJW(nTABFJW)
INTEGER         :: nprim1,nprim2,nprim3,I,J
INTEGER         :: INDADS1(nPrim),INDADS2(nPrim),INDADS3(nPrim)
REAL(REALK)     :: D2JP36,WVALU(Nprim),WVAL,WVALS1(Nprim),WVALS2(Nprim),WVALS3(Nprim)
REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0
REAL(REALK), PARAMETER :: D1 = 1.D0, SMALL = 0.000001D0

Nprim1 = 0
Nprim2 = 0
Nprim3 = 0
!D2JP36 = DFLOAT(2*JMAX + 36)
D2JP36 = 2*JMAX + 36

!DO I = 1, Nprim
!   WVALU(I)= alpha(I)*R2(I)
!ENDDO

DO I = 1, Nprim
   WVAL = alpha(I)*R2(I)
   IF (WVAL .LT. SMALL) THEN          ! 0 < WVAL < 0.000001
      RJ000(I,0) = D1
      DO J=1,JMAX
         RJ000(I,J)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
      ENDDO
   ELSEIF (WVAL .LT. D12) THEN          ! 0.000001 < WVAL < 12 
      Nprim1 = Nprim1 + 1
      INDADS1(Nprim1) = I
      WVALS1(Nprim1)  = WVAL
   ELSE IF (WVAL .LE. D2JP36) THEN  ! 12 < WVAL < (2J+36) 
      Nprim2 = Nprim2 + 1
      INDADS2(Nprim2) = I
      WVALS2(Nprim2)  = WVAL
   ELSE                             ! (2J+36) < WVAL 
      Nprim3 = Nprim3 + 1
      INDADS3(Nprim3) = I
      WVALS3(Nprim3)  = WVAL
   END IF
ENDDO

IF (Nprim1 .GT. 0) THEN
!  WVAL < 12
   CALL GETGAM1(RJ000,nPrim1,nPrim,INDADS1,WVALS1,Jmax,TABFJW,nTABFJW)
END IF
IF (Nprim2 .GT. 0) THEN
!  Near asymptotic region   12 < WVAL < (2J+36) 
   CALL GETGAM2(RJ000,nPrim2,nPrim,INDADS2,WVALS2,Jmax)
END IF
IF (Nprim3 .GT. 0) THEN
!  Asymptotic region      (2J+36) < WVAL 
   CALL GETGAM3(RJ000,nPrim3,nPrim,INDADS3,WVALS3,Jmax)
ENDIF

CALL SCALE_RJ000(RJ000,prefactor,alpha,nPrim,JMAX)

END SUBROUTINE BUILD_RJ000_OP

SUBROUTINE BUILD_RJ000_ZERO(RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,Prefactor)
IMPLICIT NONE
INTEGER         :: I,nprim,Jmax,J,Lupri,Iprint
REAL(REALK)     :: RJ000(nPrim,0:Jmax),prefac(Nprim),prefac2(Nprim)
REAL(REALK)     :: Prefactor(nprim)
REAL(REALK)     :: alpha(nprim),FAC
REAL(REALK), PARAMETER :: D2 = 2.D0,D1 = 1.D0

CALL DCOPY(nprim,prefactor,1,prefac2,1)
CALL DCOPY(NPrim,Prefac2,1,RJ000(1:nPrim,0),1)
IF (JMAX .GT. 0) THEN      
   DO I = 1, Nprim
      prefac(I)= -D2*alpha(I)
   ENDDO
   DO J = 1, JMAX-1
      FAC = D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
      DO I = 1, Nprim
         RJ000(I,J) = FAC*prefac(I)*prefac2(I)
      ENDDO
      DO I = 1, Nprim
         prefac2(I) = prefac(I)*prefac2(I)
      ENDDO
   ENDDO
   FAC = D1/(2*JMAX + 1)      
   DO I = 1, Nprim
      RJ000(I,JMAX) = FAC*prefac(I)*prefac2(I)
   ENDDO
END IF

END SUBROUTINE BUILD_RJ000_ZERO

SUBROUTINE GETGAM1(RJ000,nPrim1,nPrim,INDADS,WVALS,Jmax,TABFJW,nTABFJW)
IMPLICIT NONE
INTEGER         :: nPrim,Jmax,nPrim1,INDADS(nPrim1),nTABFJW
REAL(REALK)     :: RJ000(nPrim,0:Jmax),TABFJW(nTABFJW),WVALS(Nprim1)
REAL(REALK)     :: WVAL,WDIF,INTTEMP(nPrim1,0:Jmax),REXPW(Nprim1),FCT
INTEGER         :: ISTRT0,I,IPNT,ISTART,J,M,MP1
REAL(REALK), PARAMETER :: HALF = 0.5D0
REAL(REALK), PARAMETER :: D1 = 1.D0, D10 = 10.D0
REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0
REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6.D0, COEF4 = D1/24.D0
REAL(REALK), PARAMETER :: COEF5 = - D1/120.D0, COEF6 = D1/720.D0

ISTRT0 = 1 + 121*JMAX
DO I = 1, Nprim1
   WVAL = WVALS(I)
   IPNT = NINT(D10*WVAL)
!  WDIF = WVAL - TENTH*DFLOAT(IPNT)
   WDIF = WVAL - TENTH*IPNT
   ISTART = ISTRT0 + IPNT
   INTTEMP(I,JMAX) = (((((COEF6*TABFJW(ISTART + 726)*WDIF&
        &                         + COEF5*TABFJW(ISTART + 605))*WDIF&
        &                         + COEF4*TABFJW(ISTART + 484))*WDIF&
        &                         + COEF3*TABFJW(ISTART + 363))*WDIF&
        &                         + COEF2*TABFJW(ISTART + 242))*WDIF&
        &                         - TABFJW(ISTART + 121))*WDIF + TABFJW(ISTART)
ENDDO
IF (JMAX .GT. 0) THEN
!   DO I = 1, Nprim1
!      REXPW(I) = HALF*EXP(-WVALS(I))
!   ENDDO

   M = MOD(NPrim1,4)
   IF(M .EQ. 0) THEN
      DO I = 1,NPrim1,4
         REXPW(I)   = HALF*EXP(-WVALS(I))
         REXPW(I+1) = HALF*EXP(-WVALS(I+1))
         REXPW(I+2) = HALF*EXP(-WVALS(I+2))
         REXPW(I+3) = HALF*EXP(-WVALS(I+3))
      ENDDO
   ELSE
      DO I = 1,M
         REXPW(I) = HALF*EXP(-WVALS(I))
      ENDDO
      IF(NPrim1 .GT. 4)THEN
         MP1=M+1
         DO I = MP1,NPrim1,4
            REXPW(I)   = HALF*EXP(-WVALS(I))
            REXPW(I+1) = HALF*EXP(-WVALS(I+1))
            REXPW(I+2) = HALF*EXP(-WVALS(I+2))
            REXPW(I+3) = HALF*EXP(-WVALS(I+3))
         ENDDO
      ENDIF
   ENDIF


   DO J = JMAX - 1, 0, -1
!     FCT = D2/DFLOAT(2*J + 1)
      FCT = D2/(2*J + 1)
      DO I = 1, Nprim1
         INTTEMP(I,J) = FCT*(WVALS(I)*INTTEMP(I,J+1) + REXPW(I))
      ENDDO
   ENDDO
   DO J = 1, JMAX
!      DO I = 1, Nprim1
!         RJ000(INDADS(I),J) = INTTEMP(I,J)
!      ENDDO

      M = MOD(NPrim1,4)
      IF(M .EQ. 0) THEN
         DO I = 1,NPrim1,4
            RJ000(INDADS(I),J) = INTTEMP(I,J)
            RJ000(INDADS(I+1),J) = INTTEMP(I+1,J)
            RJ000(INDADS(I+2),J) = INTTEMP(I+2,J)
            RJ000(INDADS(I+3),J) = INTTEMP(I+3,J)
         ENDDO
      ELSE
         DO I = 1,M
            RJ000(INDADS(I),J) = INTTEMP(I,J)
         ENDDO
         IF(NPrim1 .GT. 4)THEN
            MP1=M+1
            DO I = MP1,NPrim1,4
               RJ000(INDADS(I),J) = INTTEMP(I,J)
               RJ000(INDADS(I+1),J) = INTTEMP(I+1,J)
               RJ000(INDADS(I+2),J) = INTTEMP(I+2,J)
               RJ000(INDADS(I+3),J) = INTTEMP(I+3,J)
            ENDDO
         ENDIF
      ENDIF
      
   ENDDO
END IF
!DO I = 1, Nprim1
!   RJ000(INDADS(I),0) =  INTTEMP(I,0)
!ENDDO
M = MOD(NPrim1,4)
IF(M .EQ. 0) THEN
   DO I = 1,NPrim1,4
      RJ000(INDADS(I),0) = INTTEMP(I,0)
      RJ000(INDADS(I+1),0) = INTTEMP(I+1,0)
      RJ000(INDADS(I+2),0) = INTTEMP(I+2,0)
      RJ000(INDADS(I+3),0) = INTTEMP(I+3,0)
   ENDDO
ELSE
   DO I = 1,M
      RJ000(INDADS(I),0) = INTTEMP(I,0)
   ENDDO
   IF(NPrim1 .GT. 4)THEN
      MP1=M+1
      DO I = MP1,NPrim1,4
         RJ000(INDADS(I),0) = INTTEMP(I,0)
         RJ000(INDADS(I+1),0) = INTTEMP(I+1,0)
         RJ000(INDADS(I+2),0) = INTTEMP(I+2,0)
         RJ000(INDADS(I+3),0) = INTTEMP(I+3,0)
      ENDDO
   ENDIF
ENDIF

END SUBROUTINE GETGAM1

SUBROUTINE GETGAM2(RJ000,nPrim2,nPrim,INDADS,WVALS,Jmax)
IMPLICIT NONE
!     Near asymptotic region   12 < WVAL < (2J+36) 
INTEGER         :: nPrim,Jmax,nPrim2,INDADS(nPrim2)
REAL(REALK)     :: RJ000(nPrim,0:Jmax),WVALS(Nprim2)
INTEGER         :: I,J
REAL(REALK)     :: REXPW(Nprim2),GVAL,RWVAL
REAL(REALK)     :: INTTEMP(nPrim2,0:Jmax),FCT
REAL(REALK), PARAMETER :: D1 = 1.D0, D10 = 10.D0
REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0
REAL(REALK), PARAMETER :: HALF = 0.5D0
REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092D0 
REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686D0   
REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909D0
REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346D0
REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730D00
REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2

  DO I = 1, Nprim2
     REXPW(I)   = HALF*EXP(-WVALS(I))
  ENDDO
  DO I = 1, Nprim2
     WVALS(I) = D1/WVALS(I)
  ENDDO
  DO I = 1, Nprim2
     RWVAL = WVALS(I)
     GVAL = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
     INTTEMP(I,0) = SQRPIH*SQRT(RWVAL) - REXPW(I)*GVAL*RWVAL
  ENDDO
  DO I = 1, Nprim2
     RJ000(INDADS(I),0) = INTTEMP(I,0)
  ENDDO
  DO J = 1, JMAX
!    FCT = DFLOAT(J) - HALF
     FCT = J - HALF
     DO I = 1, Nprim2
        INTTEMP(I,J) = (FCT*INTTEMP(I,J-1) - REXPW(I))*WVALS(I)
     ENDDO
     DO I = 1, Nprim2
        RJ000(INDADS(I),J) = INTTEMP(I,J)
     ENDDO
  ENDDO
END SUBROUTINE GETGAM2

SUBROUTINE GETGAM3(RJ000,nPrim3,nPrim,INDADS,WVALS,Jmax)
!     Asymptotic region      (2J+36) < WVAL 
IMPLICIT NONE
INTEGER         :: nPrim,Jmax,nPrim3,INDADS(nPrim3)
REAL(REALK)     :: RJ000(nPrim,0:Jmax),WVALS(Nprim3)
INTEGER         :: I,J
REAL(REALK)     :: REXPW(Nprim3),GVAL,RWVAL,FACTOR
REAL(REALK)     :: INTTEMP(nPrim3,0:Jmax),FCT
Real(realk), parameter :: PI=3.14159265358979323846D0
REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0
REAL(REALK), PARAMETER :: HALF = 0.5D0
REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI

  DO I = 1, Nprim3
     WVALS(I) = PID4/WVALS(I)
  ENDDO
  DO I = 1, Nprim3
     INTTEMP(I,0) = SQRT(WVALS(I))
  ENDDO
  DO I = 1, Nprim3
     RJ000(INDADS(I),0) = INTTEMP(I,0)
  ENDDO
  DO J = 1, JMAX
!    FACTOR = PID4I*(DFLOAT(J) - HALF)
     FACTOR = PID4I*(J - HALF)
     DO I = 1, Nprim3
        INTTEMP(I,J) = FACTOR*INTTEMP(I,J-1)*WVALS(I)
     ENDDO
     DO I = 1, Nprim3
        RJ000(INDADS(I),J) = INTTEMP(I,J)
     ENDDO
  ENDDO

END SUBROUTINE GETGAM3

SUBROUTINE SCALE_RJ000(RJ000,Prefactor,alpha,nprim,Jmax)
IMPLICIT NONE
INTEGER         :: I,nprim,Jmax,J
REAL(REALK)     :: RJ000(nPrim,0:Jmax),prefac(Nprim),prefac2(Nprim)
REAL(REALK)     :: Prefactor(nprim)
REAL(REALK)     :: alpha(nprim)
REAL(REALK), PARAMETER :: D2 = 2.D0

IF (JMAX .GT. 0) THEN
   CALL DCOPY(nprim,prefactor,1,prefac2,1)
   DO I = 1, Nprim
      prefac(I)= -D2*alpha(I)
   ENDDO

   DO J = 0, JMAX-1
      DO I = 1, Nprim
         RJ000(I,J) = prefac2(I)*RJ000(I,J)
      ENDDO
      DO I = 1, Nprim
         prefac2(I) = prefac(I)*prefac2(I)
      ENDDO
   ENDDO
   DO I = 1, Nprim
      RJ000(I,JMAX) = prefac2(I)*RJ000(I,JMAX)
   ENDDO
ELSE
   DO I = 1, Nprim
      RJ000(I,0) = prefactor(I)*RJ000(I,0)
   ENDDO
ENDIF

END SUBROUTINE SCALE_RJ000

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
   SHAREDTUV%TABFJW(IADR) = D1/DENOM
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
   SHAREDTUV%TABFJW(IADR) = REXPW*SUM
   DENOM = D2MAX1
   JADR = IADR
   DO J = 1,JMAX
      DENOM = DENOM - D2
      SHAREDTUV%TABFJW(JADR - 121) = (SHAREDTUV%TABFJW(JADR)*D2WAL + REXPW)/DENOM
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
 
DEALLOCATE(Orb%SPH_MAT)
NULLIFY(Orb%SPH_MAT)

END SUBROUTINE FREE_ORBITAL


SUBROUTINE FREE_Integrand(PQ)
IMPLICIT NONE
TYPE(Integrand)       :: PQ
 
DeAllocate(PQ%distance)
DeAllocate(PQ%squaredDistance)
DeAllocate(PQ%exponents)
DeAllocate(PQ%reducedExponents)
DeAllocate(PQ%integralPrefactor)
DeAllocate(PQ%iprimP)
DeAllocate(PQ%iprimQ)
Nullify(PQ%distance)
Nullify(PQ%squaredDistance)
Nullify(PQ%exponents)
Nullify(PQ%reducedExponents)
Nullify(PQ%integralPrefactor)
Nullify(PQ%iprimP)
Nullify(PQ%iprimQ)

END SUBROUTINE free_Integrand

SUBROUTINE CLEAR_integral(INTEGRAL)
TYPE(Integralitem)      :: integral
 
DEALLOCATE(INTEGRAL%Wtuv)
NULLIFY(INTEGRAL%Wtuv)
 
END SUBROUTINE


SUBROUTINE AddToCoulombMat(JMat,Jdim1,Jdim2,INPUT,integral,iA,iB,iC,iD,startA,startB,startC,startD,Dmat,Ddim1,Ddim2)
implicit none
Type(IntegralInput)  :: Input
Real(realk)          :: integral
Integer              :: Jdim1,Jdim2,Ddim1,Ddim2
Real(realk)          :: Jmat(Jdim1,Jdim2),Dmat(Ddim1,Ddim2)
Integer              :: iA,iB,iC,iD,startA,startB,startC,startD
!
Real(realk) :: dtemp_lhs,dtemp_rhs

IF (Input%DO_JENGINE) THEN
  JMat(iA,iB) = JMat(iA,iB) + integral
  IF (INPUT%sameLHSaos) THEN
    IF (startA.NE.startB) JMAT(iB,iA) = JMAT(iB,iA) + integral
  ENDIF
ELSE
  dtemp_rhs = Dmat(iD,iC)
  if (startC.NE.startD.AND.INPUT%sameRHSaos) dtemp_rhs=dtemp_rhs+Dmat(iC,iD)
  dtemp_rhs = dtemp_rhs*Input%CoulombFactor
  
  IF (INPUT%sameODs) THEN
    dtemp_lhs = Dmat(iB,iA)
    if (startA.NE.startB.AND.INPUT%sameLHSaos) dtemp_lhs=dtemp_lhs+Dmat(iA,iB)
    dtemp_lhs = dtemp_lhs*Input%CoulombFactor
  ENDIF
  
  
  IF (INPUT%sameLHSaos.AND.INPUT%sameRHSaos.AND.INPUT%sameODs)THEN
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
    IF (startA.NE.startB) JMAT(iB,iA) = JMAT(iB,iA) + dtemp_rhs*integral
    IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
      JMAT(iC,iD) = JMAT(iC,iD) + dtemp_lhs*integral
      IF (startC.NE.startD) JMAT(iD,iC) = JMAT(iD,iC) + dtemp_lhs*integral
    ENDIF
  ELSE IF (INPUT%sameLHSaos)THEN
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
    IF (startA.NE.startB) JMAT(iB,iA) = JMAT(iB,iA) + dtemp_rhs*integral
  ELSE IF (INPUT%sameRHSaos)THEN
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
  ELSE IF (INPUT%sameODs .AND. (INPUT%sameLHSaos.NEQV.INPUT%sameRHSaos))THEN
    CALL QUIT('Progamming error, sameODs but sameLHSaos.NE.sameRHSaos!')
  ELSE IF (INPUT%sameODs)THEN
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
    IF ((startA.NE.startC).OR.(startB.NE.startD)) JMAT(iC,iD) = JMAT(iC,iD) + dtemp_lhs*integral
  ELSE                                
    JMAT(iA,iB) = JMAT(iA,iB) + dtemp_rhs*integral
  ENDIF
ENDIF
END SUBROUTINE AddToCoulombMat

SUBROUTINE AddToExchangeMat(KMat,Kdim1,Kdim2,INPUT,integral,iA,iB,iC,iD,startA,startB,startC,startD,Dmat,Ddim1,Ddim2)
implicit none
Type(IntegralInput)  :: Input
Real(realk)          :: integral
Integer              :: Kdim1,Kdim2,Ddim1,Ddim2
Real(realk)          :: KMat(Kdim1,Kdim2),Dmat(Ddim1,Ddim2)
Integer              :: iA,iB,iC,iD,startA,startB,startC,startD
!
Real(realk) :: dtemp_lhs,dtemp_rhs

integral = integral*INPUT%exchangeFactor

IF (INPUT%sameLHSaos.AND.INPUT%sameRHSaos.AND.INPUT%sameODs)THEN
  KMat(iA,iC) = KMat(iA,iC) - integral*Dmat(iB,iD)
  IF (startA.NE.startB) KMat(iB,iC) = KMat(iB,iC) - integral*Dmat(iA,iD)
  IF (startC.NE.startD) KMat(iA,iD) = KMat(iA,iD) - integral*Dmat(iB,iC)
  IF ((startA.NE.startB).AND.(startC.NE.startD)) KMat(iB,iD) = KMat(iB,iD) - integral*Dmat(iA,iC)
  IF ((startA.NE.startC).OR.(startB.NE.startD)) THEN
    KMat(iC,iA) = KMat(iC,iA) - integral*DMat(iD,iB)
    IF (startA.NE.startB) KMat(iC,iB) = KMat(iC,iB) - integral*DMat(iD,iA)
    IF (startC.NE.startD) KMat(iD,iA) = KMat(iD,iA) - integral*DMat(iC,iB)
    IF ((startA.NE.startB).AND.(startC.NE.startD)) KMat(iD,iB) = KMat(iD,iB) - integral*DMat(iC,iA)
  ENDIF
!ELSE IF (INPUT%sameRHSaos) THEN
!  KMat(iA,iC) = KMat(iA,iC) - integral*Dmat(iB,iD)
!  IF (startC.NE.startD) KMat(iA,iD) = KMat(iA,iD) - integral*Dmat(iB,iC)
!ELSE IF (INPUT%sameLHSaos) THEN
!  KMat(iA,iC) = KMat(iA,iC) - integral*Dmat(iB,iD)
!  IF (startA.NE.startB) KMat(iB,iC) = KMat(iB,iC) - integral*DMat(iA,iD)
ELSE
!  KMat(iA,iC) = KMat(iA,iC) - integral*Dmat(iB,iD)
  CALL QUIT('Progamming error AddToExchangeMat, only works if all four AOs are identical!')
ENDIF
END SUBROUTINE AddToExchangeMat

SUBROUTINE initTUVitem(sharedTUV,Input,LUPRI,IPRINT)
implicit none
TYPE(TUVitem)       :: sharedTUV
TYPE(IntegralInput) :: Input
Integer             :: LUPRI,IPRINT
!
Integer,parameter :: MAXAOJ=12,MAXDER=2
Integer,parameter :: MAXJ=4*MAXAOJ+MAXDER
Integer,parameter :: NODES=121,ORDER=7
Integer,parameter :: NTABFJ000=NODES*(MAXJ+ORDER)
Integer,parameter :: NTUVMAX = (MAXJ+1)*(MAXJ+2)*(MAXJ+3)/6
LOGICAL           :: TABULATE_BOYS
!
TABULATE_BOYS = .FALSE.
IF(INPUT%Operator(1:7) .EQ. 'Coulomb') TABULATE_BOYS = .TRUE.
IF(INPUT%Operator(1:6).EQ.'Nucrep') TABULATE_BOYS = .TRUE.
IF(INPUT%Operator(1:4) .EQ. 'Erfc') TABULATE_BOYS = .TRUE.
IF(INPUT%Operator(1:4) .EQ. 'Erf ') TABULATE_BOYS = .TRUE.

IF(TABULATE_BOYS)THEN
   !Pretabulate Gammafunction/Boysfunction 
   NULLIFY(sharedTUV%TABFJW)  !MXQN=13     MAXJ = 4*(MXQN - 1) + 2 
   SharedTUV%nTABFJW = NTABFJ000
   ALLOCATE(SharedTUV%TABFJW(NTABFJ000))       !121*(MAXJ + 7) = 6897
   CALL GAMMA_TABULATION(lupri,MAXJ,sharedTUV) !Jmax is set to 10 is that ok?
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
END SUBROUTINE initTUVitem

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

SUBROUTINE SetIntegrand(iprimP,iprimQ,distance,squaredDistance,exponents,&
     &                  reducedExponents,coulomb,nucrep,integralPrefactor,npq,sameOD,&
     &                  pdist,p2dist,qdist,q2dist,pexps,qexps,pcent,qcent,&
     &                  np,nq,pempty,qempty,IUNIT)
implicit none
integer     :: np,nq,npq,IUNIT
integer     :: iprimP(npq),iprimQ(npq)
real(realk) :: distance(npq,3),squaredDistance(npq)
real(realk) :: exponents(npq),reducedExponents(npq),integralPrefactor(npq)
real(realk) :: pdist(3),p2dist,pexps(np),pcent(3,np)
real(realk) :: qdist(3),q2dist,qexps(nq),qcent(3,nq),pq(3)
logical     :: sameOD,pempty,qempty,nucrep
!
integer :: ipq,ip,iq,idir,startq
real(realk) :: pqdist(3),pqdist2,d2,pref,p,q,p_q
Real(realk), parameter :: PI=3.14159265358979323846D0, Nill = 0.0d0, OneHalf=1.5d0
Real(realk), parameter :: Two = 2.0d0, TwoHalf=2.5d0 
Real(realk), parameter :: PIFAC = 34.986836655249725D0 !Two*PI**TwoHalf
Real(realk), parameter :: TWOPI = 6.2831853071795862D0 
!
logical :: coulomb

IF(pempty.AND.qempty) THEN
  DO idir=1,3
    pqdist(idir) = 0.0d0
  ENDDO
  pqdist2 = 0.0d0
ELSE IF (pempty) THEN
  DO idir=1,3
    pqdist(idir) = qdist(idir)
  ENDDO
  pqdist2 = q2dist
ELSE IF (qempty) THEN
  DO idir=1,3
    pqdist(idir) = pdist(idir)
  ENDDO
  pqdist2 = p2dist
ENDIF

ipq = 0
DO ip=1, np
   startq = 1
   IF (sameOD) startq = ip
   DO iq=startq, nq
      !    TEST here
      ipq = ipq+1
      iprimP(ipq) = ip
      iprimQ(ipq) = iq
      IF (pempty.OR.qempty) THEN
        distance(ipq,1) = pqdist(1)
        distance(ipq,2) = pqdist(2)
        distance(ipq,3) = pqdist(3)
        squaredDistance(ipq) = pqdist2
      ELSE
        pq(1) = pcent(1,ip) - qcent(1,iq)
        pq(2) = pcent(2,ip) - qcent(2,iq)
        pq(3) = pcent(3,ip) - qcent(3,iq)
        
        distance(ipq,1) = pq(1)
        distance(ipq,2) = pq(2)
        distance(ipq,3) = pq(3)
        squaredDistance(ipq) = pq(1)*pq(1)+pq(2)*pq(2)+pq(3)*pq(3)
      ENDIF
      p  = pexps(ip)
      q  = qexps(iq)
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

SUBROUTINE SET_FTUVbatches(F,NFTUVbatches,OD,Q,Input,Integral,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)   :: INPUT
TYPE(Overlap),pointer :: F(:)
TYPE(Overlap),pointer :: Q(:)
Integer               :: NFTUVbatches,LUPRI,IPRINT
TYPE(ODITEM)          :: OD 
TYPE(Integralitem)    :: Integral
!
Integer :: I,J,nUnique,ftuvindex,np,ndmat
Integer :: ODtoFTUVindex(OD%nbatches),nPrimOD(OD%nbatches),Identifier(OD%nbatches)
Integer :: UniqeIdentfiers(OD%nbatches),FTUVprim(OD%nbatches)
Integer :: FTUVminAng(OD%nbatches),FTUVmaxAng(OD%nbatches),FTUVntuv(OD%nbatches)
Integer :: offPrim(OD%nbatches)
Logical :: unique
Integer :: nPrimQ
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
 CALL SET_Overlap(Q(I),nPrimQ,Input,Integral,OD%BATCH(I),2,LUPRI,IPRINT,'RHS')
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
ALLOCATE(F%orbital1%startOrbital(1))
ALLOCATE(F%orbital2%startOrbital(1))
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
F%orbital1%startOrbital(1) = 1
F%orbital2%startOrbital(1) = 1
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
real(realk),ALLOCATABLE :: INTTEMP(:,:)!nPrim,TUV
real(realk),ALLOCATABLE :: INTTEMP2(:,:)!nPrim,TUV
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

IF (JMAX+CARTORDER .EQ. 0) THEN
   !Special case JMAX = 0 A SIMPEL OVERLAP MATRIX WITH ONLY S ORBITALS
   !WRITE(LUPRI,*)'SPECIAL CASE OF CARTISIAN MULTIPOLE MOMENTS - OVERLAP'
   NULLIFY(INTEGRAL%Wtuv)
   ALLOCATE(INTEGRAL%Wtuv(Nprim,nTUV))
   CALL DCOPY(nPrim,M000,1,Integral%Wtuv,1)
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

      NULLIFY(INTEGRAL%Wtuv)
      ALLOCATE(INTEGRAL%Wtuv(Nprim,nTUV*nEFG))
      CALL DZERO(INTEGRAL%Wtuv,nPrim*nTUV*nEFG)

 !     WRITE(LUPRI,*)'QQ ALLOCATION IS nprim=',nprim
 !     WRITE(LUPRI,*)'QQ ALLOCATION IS nTUV =',nTUV
 !     WRITE(LUPRI,*)'QQ ALLOCATION IS nEFG =',nEFG

 !     WRITE(LUPRI,*)'CARTORDER',CARTORDER,'CALL Multipole recurrence'
      CALL Multipole_Recurrence(INTEGRAL%Wtuv,M000,Jstart,JMAX,CARTORDER,&
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
                                &(INTEGRAL%Wtuv(I,INTEGRAL%TUV%TUVINDEX(T,U,V)+(efg-1)*nTUV),I=1,NPrim)
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

END SUBROUTINE GET_MOMTUV

SUBROUTINE Multipole_Recurrence(OUTPUT,M000,Jstart,JMAX,CARTORDER,&
        & Xpq,Ypq,Zpq,Xorder,Yorder,Zorder,NPrim,nTUV,nEFG,TUVindex,&
        & zeroX, zeroY, zeroZ,EXPONENTS,lupri)
IMPLICIT NONE
INTEGER                 :: nEFG,NTUV
real(realk)             :: OUTPUT(nPrim,NTUV*nEFG),M000(nPrim)
INTEGER                 :: J,K,T,U,V,TUV,JVAL,JMAX
INTEGER                 :: nPrim
INTEGER,PARAMETER       :: MAXJ = 50
INTEGER                 :: TUVindex(0:MAXJ,0:MAXJ,0:MAXJ)
real(realk)             :: Xpq(nPrim),Ypq(nPrim),Zpq(nPrim),EXPONENTS(nPrim)
real(realk)             :: P(nPrim)
INTEGER                 :: MAXT,MAXU,MAXV,TMIN1,M1T,M2T,M1U,M2U
INTEGER                 :: M1V,M2V,VMIN1,UMIN1,IUMAX,I,e,f,g,efg
INTEGER                 :: zeroX,zeroY,zeroZ,Xorder,Yorder,Zorder
INTEGER                 :: Tm1,Tp1,jstart,CARTORDER,lupri,C,nefg2
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
Nefg2=(CARTORDER+1)*(CARTORDER+2)*(CARTORDER+3)/6
IF(Nefg.ne.nefg2)CALL QUIT('something wrong in Multipole recurrence')
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
                              OUTPUT(I,TUV+(efg-1)*nTUV)=MXYZ(I,e,T,f,U,g,V)
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

SUBROUTINE BUILD_ERF_RJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,OMEGA,IPQXYZ)
IMPLICIT NONE
INTEGER         :: nPrim,LUPRI,IPRINT,Jmax,IPQXYZ,ZeroX,ZeroY,ZeroZ
REAL(REALK)     :: RJ000(nPrim,0:Jmax),OMEGA,X0,Y0,Z0
REAL(REALK)     :: Prefactor(nPrim),alpha(nPrim)
TYPE(Integrand) :: PQ
TYPE(integralitem) :: integral

CALL DCOPY(nPrim,PQ%integralPrefactor,1,Prefactor,1)
CALL DCOPY(nPrim,PQ%reducedExponents,1,alpha,1)
!COPY BECAUSE BUILD_ERF2_RJ000 changes the values of prefactor and ALPHA
CALL BUILD_ERF2_RJ000(OMEGA,RJ000,nPrim,Jmax,LUPRI,IPRINT,&
     &alpha,PQ%squaredDistance,&
     &Prefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW,IPQXYZ)

END SUBROUTINE BUILD_ERF_RJ000

SUBROUTINE BUILD_ERF2_RJ000(OMEGA,RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW,IPQXYZ)
!  Two-electron integrals over the operator: erf(omega r12)/r12 
!  with beta = omega^2/(alpha+omega^2) 
IMPLICIT NONE
INTEGER         :: nPrim,Jmax,Lupri,Iprint,nTABFJW
REAL(REALK)     :: RJ000(nPrim,0:Jmax),Prefactor(nPrim)
REAL(REALK)     :: alpha(nPrim),R2(nprim),TABFJW(nTABFJW)!,Prefactor(nprim)
!INTEGER         :: nprim1,nprim2,nprim3,I
REAL(REALK)     :: D2JP36,WVAL!,WVALU(Nprim)
INTEGER         :: INDADS1(nPrim),INDADS2(nPrim),INDADS3(nPrim),NODS
INTEGER         :: INDADR(nPrim),IPQXYZ,I,J,nprim1,nprim2,nprim3
REAL(REALK)     :: FACINT(nPRim),WVALS1(Nprim),WVALS2(Nprim),WVALS3(Nprim),prefac(Nprim)
REAL(REALK)     :: OMEGA,OMEGA2,BETAAT,SQBETA,FAC,THRSH,WVALU(Nprim)
REAL(REALK), PARAMETER :: DP5 = -0.5D0, D1 = 1.D0, D0 = 0.D0,D12 = 12.D0
REAL(REALK), PARAMETER :: D2 = 2.D0
REAL(REALK), PARAMETER :: SMALL = 0.000001D0

! One-center Integrals
! ====================
! Note: There should be no testing for small integrals since
! this may in the case of one-center integrals introduce
! numerical instabilities for large exponents.

OMEGA2=OMEGA*OMEGA
IF (IPQXYZ.EQ.7) THEN
   DO I = 1, NPrim
      BETAAT = OMEGA2 / (ALPHA(I) + OMEGA2)
      SQBETA = SQRT(BETAAT)
!      WRITE(LUPRI,*)'Rj000(',I,',',0,')=',SQBETA,'*',prefactor(I)
      Prefactor(I) = Prefactor(I)*SQBETA
      ALPHA(I) = ALPHA(I)*BETAAT
   ENDDO

   CALL BUILD_RJ000_ZERO(RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,Prefactor)

!     Multicenter Integrals
!     =====================
!

ELSE

   NODS  = 0
   DO I = 1, NPrim
!      IF (ABS(Prefactor(I)) .GT. THRSH) THEN
         NODS = NODS + 1
         BETAAT = OMEGA2 / (ALPHA(I) + OMEGA2)
         SQBETA = SQRT(BETAAT)
         Prefactor(I) = Prefactor(I)*SQBETA
         ALPHA(I) = ALPHA(I)*BETAAT         
         WVALU(NODS) = ALPHA(I)*R2(I)
         INDADR(NODS) = I
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
      IF (WVAL .LT. SMALL) THEN          ! 0 < WVAL < 0.000001
         RJ000(I,0) = D1
         DO J=1,JMAX
            RJ000(I,J)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
         ENDDO
      ELSEIF (WVAL .LT. D12) THEN          ! 0.000001 < WVAL < 12 
         Nprim1 = Nprim1 + 1
         INDADS1(Nprim1) = INDADR(I)
         WVALS1(Nprim1)  = WVAL
      ELSE IF (WVAL .LE. D2JP36) THEN  ! 12 < WVAL < (2J+36) 
         Nprim2 = Nprim2 + 1
         INDADS2(Nprim2) = INDADR(I)
         WVALS2(Nprim2)  = WVAL
      ELSE                             ! (2J+36) < WVAL 
         Nprim3 = Nprim3 + 1
         INDADS3(Nprim3) = INDADR(I)
         WVALS3(Nprim3)  = WVAL
      END IF
   ENDDO
   IF (Nprim1 .GT. 0) THEN
      !  WVAL < 12
      CALL GETGAM1(RJ000,nPrim1,nPrim,INDADS1,WVALS1,Jmax,TABFJW,nTABFJW)
   END IF
   IF (Nprim2 .GT. 0) THEN
      !  Near asymptotic region   12 < WVAL < (2J+36) 
      CALL GETGAM2(RJ000,nPrim2,nPrim,INDADS2,WVALS2,Jmax)
   END IF
   IF (Nprim3 .GT. 0) THEN
      !  Asymptotic region      (2J+36) < WVAL 
      CALL GETGAM3(RJ000,nPrim3,nPrim,INDADS3,WVALS3,Jmax)
   ENDIF

   CALL SCALE_RJ000(RJ000,prefactor,alpha,nPrim,JMAX)
END IF
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
      CALL OUTPUT(RJ000,1,NPrim,1,JMAX+1,NPrim,JMAX+1,1,LUPRI)
   END IF
END IF

END SUBROUTINE BUILD_ERF2_RJ000

END MODULE integraldriver

