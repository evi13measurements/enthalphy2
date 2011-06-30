!> @file
!> Contains the main Thermite integral drivers

!> \brief Main Thermite drivers for the calculation of integrals based on the McMurchie-Davidson scheme, using (Turbo-)hermite functions. (PCCP 2007 9, 4771) 
!> \author T. Kjaergaard and S. Reine
!> \date 2008 
MODULE integraldriver
  use TYPEDEF
  use READMOLEFILE
  use BuildBasisSet
  use ODbatches
  use OD_type
  use sphcart_matrices
  use precision
  use lstiming
  use thermite_integrals
  use thermite_OD
  use ls_util
SAVE
TYPE LINKshell
integer,pointer :: Belms(:)
integer,pointer :: IODelms(:)
INTEGER         :: DIM
END TYPE LINKshell
CONTAINS
!> \brief Main integral driver
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!>  This is the Main integral driver, which sets op the overlap distributions 
!>  (call ODbatches) from the Atomic orbitals given in input INPUT%AO.
!>  The routine then calls either the main McMurchie-Davidson driver or 
!>  specialised routines. For the calculation of exchange integrals 
!>  a specialised LinK routine is used and for the calculation of coulomb 
!>  integrals, using the jengine algorithm, another specialised routine 
!>  have been implemented. 
!> 
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param INTOUT the integral output specifications, determines how the output should be given
SUBROUTINE MAIN_INTEGRAL_DRIVER(LUPRI,IPRINT,INPUT,INTOUT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: INTOUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM),target  :: OD_LHSt,OD_RHSt
TYPE(ODITEM),pointer :: OD_LHS,OD_RHS
!Character(len=80)    :: IDENTIFIER
!
Real(realk)         :: TS,TE
INTEGER             :: idmat
!
!IDENTIFIER='THERMIT'
!Quick Exit if possible
IF (INPUT%sameODs .AND. (INPUT%sameLHSaos.NEQV.INPUT%sameRHSaos))THEN
   CALL LSQUIT('Progamming error, sameODs but sameLHSaos.NE.sameRHSos!',lupri)
ENDIF

IF(IPRINT .GT. 5) WRITE(LUPRI,*)'CALLING CREATE OD'
CALL Create_ODbatches(OD_LHSt,INPUT,'LHS',LUPRI)
OD_LHS => OD_LHSt
IF(OD_LHS%nbatches .NE. 0)THEN
   IF(INPUT%sameODs)THEN
      OD_RHS => OD_LHSt
   ELSE
      CALL Create_ODbatches(OD_RHSt,INPUT,'RHS',LUPRI)
      OD_RHS => OD_RHSt
   ENDIF
   IF(OD_RHS%nbatches .NE. 0)THEN      
      CALL PRINT_OD(OD_LHS,LUPRI,IPRINT)
      CALL PRINT_OD(OD_RHS,LUPRI,IPRINT)
      
      IF (INPUT%DO_JENGINE) THEN
         CALL LSTIMER('START ',TS,TE,LUPRI)
         CALL Jengine(OD_LHS,OD_RHS,INPUT,INTOUT,LUPRI,IPRINT)
         IF (INPUT%DO_FMM .AND. (.NOT.INPUT%SKIP_FMM )) THEN
            CALL LSTIMER('J-non  ',TS,TE,LUPRI)
            CALL JmatClassical(INTOUT,INTOUT%ndim(1),INTOUT%ndim(2),INPUT)
            CALL LSTIMER('J-cls  ',TS,TE,LUPRI)
         ELSE
            CALL LSTIMER('Jengine',TS,TE,LUPRI)
         ENDIF
      ELSEIF(INPUT%DO_LINK) THEN
        CALL LINK(OD_LHS,OD_RHS,INPUT,INTOUT,LUPRI,IPRINT)
      ELSE
         CALL LSTIMER('START ',TS,TE,LUPRI)
         CALL McMurchieDavidson(OD_LHS,OD_RHS,INPUT,INTOUT,LUPRI,IPRINT)
         IF (INPUT%DO_FMM) THEN
            CALL LSTIMER('ea-non',TS,TE,LUPRI)
            CALL electronNuclearClassic(INTOUT,INTOUT%ndim(1),INTOUT%ndim(2),INPUT)
            CALL LSTIMER('ea-cls  ',TS,TE,LUPRI)
         ENDIF
      ENDIF
   ENDIF
   IF(.NOT.Input%sameODs)THEN
      CALL FREE_ODitem(OD_RHS)
   ENDIF
   CALL FREE_ODitem(OD_LHS)
ENDIF

IF (Input%DO_GRADIENT) IPRINT = 0

END SUBROUTINE MAIN_INTEGRAL_DRIVER

!> \brief McMurchie-Davidson based integral driver
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!>  This is the Main McMurchie-Davidson based integral driver, which 
!>  sets up overlaps P and Q and loops over the left hand side (P) overlaps
!>  and 
!>  The routine then calls either the main McMurchie-Davidson driver or 
!>  specialised routines. For the calculation of exchange integrals 
!>  a specialised LinK routine is used and for the calculation of coulomb 
!>  integrals, using the jengine algorithm, another specialised routine 
!>  have been implemented. 
!> 
!>  \param OD_LHS the ODbatches belonging to the Left hand side
!>  \param OD_RHS the ODbatches belonging to the Reft hand side
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param OUTPUT the integral output specifications, determines how the output should be given
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
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
Integer               :: nPrimP
Integer,pointer       :: nPrimQ(:)
Integer,pointer       :: ODpassesIndex(:)
logical               :: dopasses,setOrbital,screen
integer               :: maxpasses,nPassTypes,iPassType,numpasses
real(realk),pointer   :: WORK(:)
integer :: LWORK,WORKLENGTH
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
IF(doPasses)THEN
   maxPasses = INPUT%maxPasses
ELSE
   maxPasses = 1
ENDIF

call mem_alloc(nPrimQ,OD_RHS%nbatches)
call mem_alloc(ODpassesIndex,OD_RHS%nbatches)
CALL initAlloc(Alloc,LUPRI,IPRINT,'Both')

NULLIFY(Q)
ALLOCATE(Q(OD_RHS%nbatches))
DO IRHS=1,OD_RHS%nbatches
   CALL SET_Overlap(Q(IRHS),nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS')
ENDDO
CALL SelectODPassTypes(ODpassesIndex,Q,OD_RHS%nbatches,nPassTypes,IPRINT,LUPRI)
CALL initAlloc(Alloc,LUPRI,IPRINT,'LHS')
CALL SET_ALLOC(Alloc,Input,OD_LHS,'LHS',IPRINT,LUPRI)

!$OMP PARALLEL PRIVATE(integral,P,PQ,ILHS,IRHS,iPassType,numpasses,Start_RHS,End_RHS,nPrimP,&
!$OMP                  setOrbital,PassQ,WORK,WORKLENGTH,LWORK,screen)    
integral%TUV => sharedTUV
CALL INIT_OVERLAP(P,Alloc,Input,OD_LHS,'LHS',1,IPRINT,LUPRI)
CALL allocateIntegrals(PQ,Integral,Input,Alloc,maxPasses,1)
CALL INIT_PASS(PassQ,Alloc,Input,OD_RHS,'RHS',maxPasses,IPRINT,LUPRI)

WORKLENGTH = MAX(Alloc%maxprimRHS*Alloc%maxContRHS,Alloc%maxprimLHS*Alloc%maxContLHS) !contractbasis needs
WORKLENGTH = MAX(WORKLENGTH,Alloc%maxijkLHS*Alloc%maxijkLHS,Alloc%maxijkRHS*Alloc%maxijkRHS) !sphericaltransform needs
IF(.NOT.Input%orderAngPrim)THEN
   WORKLENGTH = MAX(WORKLENGTH,Alloc%maxETUVlenRHS*maxpasses) !needed by DirectcontractEQ
ENDIF
call mem_alloc(WORK,WORKLENGTH)
LWORK = 1


!* LHS loop startes here
!$OMP DO SCHEDULE(DYNAMIC,1)
DO ILHS = 1,OD_LHS%nbatches
   CALL SET_INITIALIZED_Overlap(P,nPrimP,Input,SharedTUV,Integral,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS',.FALSE.,.TRUE.)
   IF (nPrimP .GT. 0) THEN
      CALL Determine_RHS_loop(INPUT,OD_RHS,ILHS,Start_RHS,End_RHS)
      DO iPassType=1,nPassTypes
         numPasses = 0
         setOrbital = .TRUE.
         !********** RHS loop startes here
         DO IRHS = Start_RHS,End_RHS
            IF(nPrimQ(IRHS) .GT. 0)THEN
               IF (ODpassesIndex(IRHS).EQ.iPassType) THEN
                  screen = getScreening(P,Q(IRHS),INPUT,LUPRI,IPRINT)
                  IF (.NOT.screen) THEN
                     IF (IPRINT.GT.3) WRITE(LUPRI,'(1X,A,I3,A,I3)') 'Overlap distributions P',ILHS,' and Q',IRHS
                     IF (doPasses .AND..NOT. ((ILHS.EQ.IRHS).AND.INPUT%sameODs)) THEN
                        numPasses = numPasses + 1
                        IF(numPasses.EQ.1)CALL SetPassOrbitals(PassQ,Q(IRHS),maxPasses,LUPRI,IPRINT)
                        CALL AddOverlapToPass(PassQ,Q(IRHS),numPasses,maxPasses,LUPRI,IPRINT)
                        CALL FinalizePass(PassQ,Q(IRHS),numPasses,LUPRI,IPRINT)
                        IF (numPasses.EQ.maxPasses) THEN
                           CALL ExplicitIntegrals(Integral,PQ,P,PassQ,INPUT,OUTPUT,ILHS,IRHS,WORK,&
                                &LWORK,WORKLENGTH,LUPRI,IPRINT)
                           CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
                           numPasses = 0
                        ENDIF
                     ELSE
                        CALL ExplicitIntegrals(Integral,PQ,P,Q(IRHS),INPUT,OUTPUT,ILHS,IRHS,&
                             &                 WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
                        CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
                     ENDIF
                  ENDIF! CS-screening
               ENDIF ! Correct passtype
            ENDIF! PS-screening
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
CALL Free_Overlap(PassQ)
call mem_dealloc(WORK)

!$OMP END PARALLEL

DO IRHS=1,OD_RHS%nbatches
   IF (nPrimQ(IRHS) .GT. 0) CALL FREE_OVERLAP(Q(IRHS))
ENDDO
DEALLOCATE(Q)
call mem_dealloc(ODpassesIndex)
call mem_dealloc(nPrimQ)
!DEALLOCATE(ODpassesIndex)
!DEALLOCATE(nPrimQ)

CALL freeTUVitem(sharedTUV,Input)
END SUBROUTINE McMurchieDavidson

!> \brief driver routine to calculate the explicit integral for a given P and Q overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!>  This is a routine which calculated the explicit 4 center integral from the 
!>  overlap P and Q given in input. It first build the integrand structure PQ
!>  which contains stuff like the reduced exponent, the integral prefactor, etc. 
!>  It then calculates the integral over zero'th order primitive hermite functions 
!>  Eq. 9.9.14 (and for instance 9.5.41) in the book. and performs 
!>  recurrence relations (see for instance Eq. 9.9.18-9.9.20) in order to obtain the full     
!>  integral over primitive hermite functions. 
!>  Finally it does the 
!>     Ecoefficient contraction  (see Eq. 9.5.1)
!>     Primitive Contraction (from primitve AO integrals to contracted AO integrals)    
!>     Spherical transformation (from cartesian to spherical harmonic AOs)
!>  both for electron 2 (contract_Q) and electron 1 (contract_P)
!>
!>  \param Integral contains arrays to store intermidiates and final integrals
!>  \param PQ contains integrand info like reduced exponents, integral prefactor, etc.
!>  \param P contains overlap distribution for the left hand side, electron 1
!>  \param Q contains overlap distribution for the right hand side, electron 2
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param OUTPUT the integral output specifications, determines how the output should be given
!>  \param ILHS the index for which LHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param IRHS the index for which RHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param WORK a temporary array to store keep temporary information 
!>  \param LWORK last unused index in work array
!>  \param WORKLENGTH the size of WORK
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
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

CALL Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,LUPRI,IPRINT)
!   Hermite 2-electron integral over general operator w
CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
CALL Contract_Q(INTEGRAL,PQ,Input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
CALL Contract_P(INTEGRAL,PQ,input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)

END SUBROUTINE ExplicitIntegrals

!> \brief Jengine based integral driver (CPL 323, 425 and JCP 114, 6572)
!> \author S. Reine
!> \date 2010
!>
!> \param OD_LHS the ODbatches belonging to the Left hand side
!> \param OD_RHS the ODbatches belonging to the Reft hand side
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param OUTPUT the integral output specifications, determines how the output should be given
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
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
TYPE(Overlap)         :: PassF
TYPE(Overlap),pointer :: Q(:)
TYPE(Overlap),pointer :: F(:)
TYPE(Integrand)       :: PQ
Integer               :: nPrimP, nPrimQ,LWORK,WORKLENGTH,maxPrimPass,nPrimPass,PassType,doLHS
real(realk),pointer   :: WORK(:)
Logical               :: screen,dopasses

IF (IPRINT.GT.5) THEN
  CALL LSHEADER(LUPRI,'Jengine')
  WRITE(LUPRI,'(1X,A)') ''
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') '***                          Jengine                           ***'
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') ''
ENDIF

CALL initTUVitem(sharedTUV,Integral,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
CALL initAlloc(Alloc,LUPRI,IPRINT,'RHS')

doPasses  = INPUT%DO_PASSES

IF(doPasses) THEN
  maxPrimPass   = 64
  WRITE(LUPRI,*)'Jengine integrals are calculated using passes with maxPrimPass=',maxPrimPass
ELSE
  maxPrimPass = 1
ENDIF

NULLIFY(Q)
ALLOCATE(Q(OD_RHS%nbatches))
CALL SET_FTUVbatches(F,NFTUVbatches,OD_RHS,Q,Input,SharedTUV,Integral,Alloc,&
     &               INPUT%NDMAT_RHS,maxPrimPass,LUPRI,IPRINT)
DEALLOCATE(Q)
Q=>F
NULLIFY(F)

CALL initAlloc(Alloc,LUPRI,IPRINT,'LHS')
CALL SET_ALLOC(Alloc,Input,OD_LHS,'LHS',IPRINT,LUPRI)

!$OMP PARALLEL PRIVATE(integral,P,PQ,ILHS,IRHS,Start_LHS,End_LHS,nPrimP,nPrimQ,WORK,LWORK,WORKLENGTH,&
!$OMP screen,PassF,nPrimPass,passType,doLHS) 
integral%TUV => sharedTUV
CALL INIT_OVERLAP(P,Alloc,Input,OD_LHS,'LHS',1,IPRINT,LUPRI)
WORKLENGTH = MAX(Alloc%maxprimRHS*Alloc%maxContRHS,&
     & Alloc%maxprimLHS*Alloc%maxContLHS) !needed in contractbasis
WORKLENGTH = MAX(WORKLENGTH,Alloc%maxijkLHS*Alloc%maxijkLHS,Alloc%maxijkRHS*Alloc%maxijkRHS) !needed in sphericaltransform
call mem_alloc(WORK,WORKLENGTH)
!ALLOCATE(WORK(WORKLENGTH))
LWORK = 1
CALL allocateIntegrals(PQ,Integral,Input,Alloc,1,INPUT%NDMAT_RHS)

!$OMP DO SCHEDULE(DYNAMIC,1)
DO ILHS = 1,OD_LHS%nbatches
   CALL SET_INITIALIZED_Overlap(P,nPrimP,Input,SharedTUV,Integral,OD_LHS%BATCH(ILHS),&
     &                          1,LUPRI,IPRINT,'LHS',.FALSE.,.TRUE.)
 IF (nPrimP.GT.0) THEN
  CALL LS_DZERO(Integral%integrals,P%nPrimitives*P%nTUV*input%ndmat_rhs)
  nPrimPass = 0
  passType  = 0
  doLHS = 0 ! this is an integer that keeps track of the RHS loop if all the RHS loop is screened away
            ! then doLHS remains 0 and the contract_P should not be called.
  DO IRHS = 1,NFTUVbatches
      screen = getScreening(P,Q(IRHS),INPUT,LUPRI,IPRINT)
      IF (.NOT.screen) THEN
        IF (IPRINT.GT.3) WRITE(LUPRI,'(1X,A,I3,A,I3)') 'Overlap ditributions P',ILHS,' and Q',IRHS
!       ******************************* Passes *****************************
!       Integrals are calculated using passes (collecting similar FTUV-batches 
!       before performing integration to reduce computation overhead)
        IF (doPasses) THEN
!         If new pass-type calculate contribution from previous pass type
          IF ((passType.NE.Q(IRHS)%passType).AND.(nPrimPass.GT.0)) THEN
              CALL JengineInnerCont(PQ,P,PassF,INTEGRAL,INPUT,ILHS,IRHS, &
     &                              INPUT%NDMAT_RHS,Input%orderAngPrim,LUPRI,IPRINT)
              nPrimPass = 0
              CALL free_Overlap(PassF)
          ENDIF
!         If the OD-batch has more primivites that the maximum do not add to pass
!         but calculate directly
          IF (Q(IRHS)%nPrimitives.GE.maxPrimPass) THEN
              CALL JengineInnerCont(PQ,P,Q(IRHS),INTEGRAL,INPUT,ILHS,IRHS, &
     &                              INPUT%NDMAT_RHS,Input%orderAngPrim,LUPRI,IPRINT)
          ELSE
            IF (nPrimPass.EQ.0) THEN
              passType = Q(IRHS)%passType
              CALL InitFTUVbatches(PassF,passType,maxPrimPass,Q(IRHS)%startAngmom,&
                                   Q(IRHS)%endAngmom,Q(IRHS)%nTUV,INPUT%NDMAT_RHS)
            ENDIF

!           Calcluate old pass first if number of primitives would exceed
!           the maximal number
            IF ((Q(IRHS)%nPrimitives+nPrimPass).GT.maxPrimPass) THEN
              CALL JengineInnerCont(PQ,P,PassF,INTEGRAL,INPUT,ILHS,IRHS, &
     &                              INPUT%NDMAT_RHS,Input%orderAngPrim,LUPRI,IPRINT)
              nPrimPass = 0
              PassF%nPrimitives = 0
            ENDIF
            CALL AddODtoFTUV(PassF,Q(IRHS)%nPrimitives,nPrimPass,Q(IRHS),INPUT%NDMAT_RHS,LUPRI,IPRINT)
            nPrimPass = nPrimPass + Q(IRHS)%nPrimitives
          ENDIF
!       ****************************** NO Passes ***************************
!       Integrals are calculated batchwise
        ELSE
            CALL JengineInnerCont(PQ,P,Q(IRHS),INTEGRAL,INPUT,ILHS,IRHS, &
     &                            INPUT%NDMAT_RHS,Input%orderAngPrim,LUPRI,IPRINT)
        ENDIF
        doLHS = 1
     ENDIF
  ENDDO
  IF(doLHS .EQ. 1)THEN
     IF (doPasses.AND.(nPrimPass.GT.0)) THEN
        CALL JengineInnerCont(PQ,P,PassF,INTEGRAL,INPUT,ILHS,IRHS, &
             &                     INPUT%NDMAT_RHS,Input%orderAngPrim,LUPRI,IPRINT)
     ENDIF
     CALL Contract_P(INTEGRAL,PQ,input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
     CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
     IF (doPasses.AND.(nPrimPass.GT.0)) THEN
        nPrimPass = 0
        CALL free_Overlap(PassF)
     ENDIF
  ENDIF
 ENDIF
ENDDO
!$OMP END DO

CALL deallocateIntegrals(PQ,Integral)
CALL FREE_OVERLAP(P)
call mem_dealloc(WORK)
!DEALLOCATE(WORK)

!$OMP END PARALLEL

DO IRHS = 1,NFTUVbatches
  IF (Q(IRHS)%nPrimitives.GT.0) CALL FREE_OVERLAP(Q(IRHS))
ENDDO
DEALLOCATE(Q)
NULLIFY(Q)

CALL freeTUVitem(sharedTUV,Input)

END SUBROUTINE Jengine

!> \brief The inner contraction required in a Jengine based integration
!> \author S. Reine
!> \date 2010
!>
!>  \param PQ contains integrand info like reduced exponents, integral prefactor, etc.
!>  \param P contains overlap distribution for the left hand side, electron 1
!>  \param Q contains overlap distribution for the right hand side, electron 2
!>  \param Integral contains arrays to store intermidiates and final integrals
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param ILHS the index for which LHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param IRHS the index for which RHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param NDMAT the number of density matrices
!>  \param orderAngPrim the ordering of the Wtuv(ntuvPQ,nPrimPQ) integrals (DEFAULT(hardcoded) .false.)
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE JengineInnerCont(PQ,P,Q,INTEGRAL,INPUT,ILHS,IRHS,NDMAT,orderAngPrim,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT) :: INPUT
integer             :: LUPRI,IPRINT,ILHS,IRHS,NDMAT
TYPE(Integralitem)  :: Integral
TYPE(Overlap)       :: P
TYPE(Overlap)       :: Q
TYPE(Integrand)     :: PQ
Logical             :: orderAngPrim
!
CALL Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,LUPRI,IPRINT)
!   Hermite 2-electron integral over general operator w 
CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
CALL ContractFTUV(INTEGRAL,P,Q,NDMAT,orderAngPrim,LUPRI,IPRINT)
END SUBROUTINE JengineInnerCont

!> \brief LinK based integral driver (JCP 109,1663) 
!> \author T. Kjaergaard
!> \date 2010
!>
!>  \param OD_LHS the ODbatches belonging to the Left hand side
!>  \param OD_RHS the ODbatches belonging to the Reft hand side
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param LSOUTPUT the integral output specifications, determines how the output should be given
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE LINK(OD_LHS,OD_RHS,INPUT,LSOUTPUT,LUPRI,IPRINT)
  use thermite_distribute_exchange
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: LSOUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
!
Integer,target        :: ILHS,STATICSIZEOFDOINT
Integer               :: IRHS,Start_RHS,End_RHS,i,j
TYPE(TUVitem),target  :: SharedTUV
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)       :: Alloc
TYPE(Overlap),pointer :: P(:),Q(:)
TYPE(Overlap),pointer :: PassQ(:)
TYPE(Integrand)       :: PQ
Integer,pointer       :: nPrimP(:)
Integer,pointer       :: nPrimQ(:)
Integer,pointer       :: ODpassesIndex(:),ODpassesIndex2(:)
Integer               :: nPassTypes,nRHSoverlaps,nLHSoverlaps
Integer,pointer       :: Maxpassfortype(:),numpasses(:)
Real(realk),pointer   :: RED_GAB_LHS(:,:)
Real(realk),pointer   :: RED_GAB_RHS(:,:)
Real(realk),pointer   :: RED_DMAT_LHS(:,:)
Real(realk),pointer   :: RED_DMAT_RHS(:,:)
Integer               :: dim1,dim2,dim3,dim4,A,B,C,D,idmat,endC
Integer               :: LISTSIZE!,LISTSIZE1,LISTSIZE2
Integer               :: nB,nC,nD,OLDI,RHSINDEX,startB,startD,I2,JK
Integer,pointer   :: LIST(:)!,LIST2(:),LIST1(:)
LOGICAL,pointer   :: DoINT(:)
real(realk),pointer :: maxLHSGAB(:),maxRHSGAB(:)
Integer                 :: LWORK,WORKLENGTH
real(realk),pointer :: WORK(:)
!real(realk),allocatable :: SORTING(:)
!integer,allocatable     :: bshell(:)
real(realk)           :: maxLHSELM,maxRHSELM,DMATELM1,DMATELM2,MAXDMAT,CPUTIME1,CPUTIME2
!LOGICAL               :: A_LOOP_DONE,B_LOOP_DONE,C_LOOP_DONE,D_LOOP_DONE,UNIQUE
TYPE(LINKshell),pointer  :: brashell(:),ketshell(:),ML(:)
real(realk)           :: CPUTimeLINK,WALLTIMELINK
real(realk)           :: CPUTIMESTART,WALLTIMESTART,CPUTIMEEND,WALLTIMEEND
logical               :: dopasses,set_orbital1
integer               :: maxpasses,iPassType,startc
integer,pointer   :: batchindex1(:),batchindex2(:)
integer,pointer   :: SIZEOFDOINT
logical        :: NOELEMENTSADDED,LHS,screen,dalink,doneBuild,DRHS_SYM,DLHS_SYM,sameODs
real(realk)    :: CS_THRESHOLD
type(lstensor) :: redtensor
CALL LS_GETTIM(CPUTIMESTART,WALLTIMESTART)
IF(.NOT. INPUT%CS_SCREEN)THEN
   WRITE(LUPRI,*)'Link requires screening, but you have turned off screening'
   CALL LSQUIT('Link requires screening, but you have turned off screening',lupri)
ENDIF
dalink = INPUT%DO_DALINK
DRHS_SYM = INPUT%DRHS_SYM
DLHS_SYM = INPUT%DLHS_SYM
CS_THRESHOLD = INPUT%CS_THRESHOLD 
sameODs = INPUT%sameODs
IF (IPRINT.GT.5) THEN
  IF(DALINK)THEN
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
call mem_alloc(RED_GAB_LHS,dim1,dim2)
call mem_alloc(RED_GAB_RHS,dim3,dim4)
CALL LS_DZERO(RED_GAB_LHS,dim1*dim2)
CALL LS_DZERO(RED_GAB_RHS,dim3*dim4)
call Build_reduce_lstensor(redtensor,input%LST_GAB_LHS,dim1,dim2,1,1)
call Build_full_2dim_from_lstensor(redtensor,RED_GAB_LHS,dim1,dim2)
call lstensor_free(redtensor)
call Build_reduce_lstensor(redtensor,input%LST_GAB_RHS,dim3,dim4,1,1)
call Build_full_2dim_from_lstensor(redtensor,RED_GAB_RHS,dim3,dim4)
call lstensor_free(redtensor)

!CALCULATE MAXGAB VECTOR AND MAX ELEMENT IN GAB
call mem_alloc(maxLHSGAB,dim1)
CALL MAXGABELM(RED_GAB_LHS,dim1,dim2,maxLHSGAB,maxLHSELM)
call mem_alloc(maxRHSGAB,dim3)
CALL MAXGABELM(RED_GAB_RHS,dim3,dim4,maxRHSGAB,maxRHSELM)

doPasses  = INPUT%DO_PASSES
IF(doPasses)THEN
   maxPasses = INPUT%maxPasses
ELSE
   maxPasses = 1
ENDIF
!IF(doPasses) WRITE(LUPRI,*)'Integrals are calculated using passes with maxpasses=',maxPasses

CALL initTUVitem(sharedTUV,Integral,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
CALL initAlloc(Alloc,LUPRI,IPRINT,'Both')

call mem_alloc(nPrimQ,OD_RHS%nbatches)
NULLIFY(Q)
ALLOCATE(Q(OD_RHS%nbatches))
nRHSoverlaps = 0
DO IRHS=1,OD_RHS%nbatches
   nRHSoverlaps = nRHSoverlaps + 1
   CALL SET_Overlap(Q(nRHSoverlaps),nPrimQ(IRHS),Input,SharedTUV,Integral,Alloc,OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,'RHS')
   IF(nPrimQ(IRHS) .EQ. 0)THEN
      nRHSoverlaps = nRHSoverlaps - 1 !IRHS NOT INCLUDED 
   ENDIF
ENDDO
call mem_alloc(ODpassesIndex,nRHSoverlaps)
CALL SelectODPassTypes(ODpassesIndex,Q,nRHSoverlaps,nPassTypes,IPRINT,LUPRI)

call mem_alloc(nPrimP,OD_LHS%nbatches)
NULLIFY(P)
ALLOCATE(P(OD_LHS%nbatches))
nLHSoverlaps = 0
DO ILHS=1,OD_LHS%nbatches
   nLHSoverlaps = nLHSoverlaps + 1
   CALL SET_Overlap(P(nLHSoverlaps),nPrimP(ILHS),Input,SharedTUV,Integral,Alloc,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,'LHS')
   IF(nPrimP(ILHS) .EQ. 0)THEN
      nLHSoverlaps = nLHSoverlaps - 1 !IRHS NOT INCLUDED 
   ENDIF
ENDDO

! determine significant bra shell pairs

LHS=.TRUE.
NULLIFY(brashell)
ALLOCATE(brashell(dim1))
CALL DETERMINE_SHELL_PAIRS(dim1,dim2,RED_GAB_LHS,CS_THRESHOLD,&
        &maxRHSELM,brashell,OD_LHS,LHS,INPUT%sameRHSaos,nPrimP)

LHS=.FALSE.
NULLIFY(ketshell)
ALLOCATE(ketshell(dim3))
CALL DETERMINE_SHELL_PAIRS(dim3,dim4,RED_GAB_RHS,CS_THRESHOLD,&
        &maxLHSELM,ketshell,OD_RHS,LHS,INPUT%sameRHSaos,nPrimQ)

call mem_alloc(RED_DMAT_RHS,dim1,dim3)
call Build_reduce_lstensor(redtensor,input%LST_DRHS,dim1,dim3,1,1)
call Build_full_2dim_from_lstensor(redtensor,RED_DMAT_RHS,dim1,dim3)
call lstensor_free(redtensor)

IF(DALINK)THEN
   IF(SameODs)THEN
      RED_DMAT_LHS => RED_DMAT_RHS 
   ELSE
      call mem_alloc(RED_DMAT_LHS,dim2,dim4)
      call Build_reduce_lstensor(redtensor,input%LST_DLHS,dim2,dim4,1,1)
      call Build_full_2dim_from_lstensor(redtensor,RED_DMAT_LHS,dim2,dim4)
      call lstensor_free(redtensor)
   ENDIF
   call mem_alloc(batchindex1,nRHSoverlaps)
   call mem_alloc(batchindex2,nRHSoverlaps)
   DO C=1,dim3
      DO nD=1,ketshell(C)%DIM
         D=ketshell(C)%belms(nD)
         IRHS=ketshell(C)%IODelms(nD)
         batchindex1(IRHS)=C
         batchindex2(IRHS)=D
      ENDDO
   ENDDO
ENDIF

NULLIFY(ML)
ALLOCATE(ML(dim1))
CALL DETERMINE_BRAKET_PAIRS(dim1,dim2,dim3,RED_DMAT_RHS,&
     &maxLHSGAB,maxRHSGAB,CS_THRESHOLD,ML,sameODs)

!$OMP PARALLEL PRIVATE(integral,PQ,IRHS,RHSINDEX,A,nB,B,ILHS,nC,C,nD,D,&
!$OMP LIST,LISTSIZE,WORK,LWORK,WORKLENGTH,DoINT,NOELEMENTSADDED,DMATELM1,&
!$OMP DMATELM2,MAXDMAT,iPasstype,numpasses,PassQ,doneBuild,ODpassesIndex2,&
!$OMP SIZEOFDOINT,STATICSIZEOFDOINT)  

call mem_alloc(ODpassesIndex2,nPassTypes)
NULLIFY(PassQ)
ALLOCATE(PassQ(nPassTypes))
DO iPassType=1,nPassTypes
   doneBUILD = .FALSE.
   DO IRHS=1,nRHSoverlaps
      IF(ODpassesIndex(IRHS).EQ.iPassType) THEN
         CALL INIT_PASS(PassQ(iPassType),Alloc,Input,OD_RHS,'RHS',maxPasses,IPRINT,LUPRI)
         CALL SetPassOrbitals(PassQ(iPassType),Q(IRHS),maxPasses,LUPRI,IPRINT)
         ODpassesIndex2(iPassType) = IRHS
         doneBUILD = .TRUE.
      ENDIF
      IF(doneBUILD)EXIT
   ENDDO
ENDDO
call mem_alloc(numpasses,nPassTypes)
numpasses = 0
integral%TUV => sharedTUV
CALL allocateIntegrals(PQ,Integral,Input,Alloc,maxpasses,INPUT%NDMAT_RHS)
!CALL INIT_PASS(PassQ,Alloc,Input,OD_RHS,'RHS',maxPasses,IPRINT,LUPRI)
WORKLENGTH = MAX(Alloc%maxprimRHS*Alloc%maxContRHS,&
     & Alloc%maxprimLHS*Alloc%maxContLHS) !needed for contractbasis
WORKLENGTH = MAX(WORKLENGTH,Alloc%maxijkLHS*Alloc%maxijkLHS,Alloc%maxijkRHS*Alloc%maxijkRHS) !needed in sphericaltransform
IF(.NOT.Input%orderAngPrim)THEN
   WORKLENGTH = MAX(WORKLENGTH,Alloc%maxETUVlenRHS*maxpasses) !needed by DirectcontractEQ
ENDIF
LWORK = 1
call mem_alloc(WORK,WORKLENGTH)
call mem_alloc(LIST,dim3*dim4)
call mem_alloc(DoINT,dim3*dim4)
NULLIFY(SIZEOFDOINT)
IF(sameODs)THEN
   SIZEOFDOINT => ILHS
ELSE
   STATICSIZEOFDOINT = dim3*dim4
   SIZEOFDOINT => STATICSIZEOFDOINT
ENDIF
!$OMP DO SCHEDULE(DYNAMIC,1)
DO A = 1,dim1
   DO nB=1,brashell(A)%DIM
      B=brashell(A)%Belms(nB)
      ILHS=brashell(A)%IODelms(nB)
      !CREATING DoINT LIST from the AC element of the density matrix
      DO IRHS=1,SIZEOFDOINT!ILHS
         DoINT(IRHS) = .FALSE.
      ENDDO
      DO nC=1,ML(A)%DIM
         C = ML(A)%Belms(nC) 
         NOELEMENTSADDED = .TRUE.
         DO nD=1,ketshell(C)%DIM
            D=ketshell(C)%Belms(nD)
            IF(RED_DMAT_RHS(A,C)*RED_GAB_RHS(A,B)&
                 &*RED_GAB_RHS(C,D) .LE. CS_THRESHOLD )EXIT
            NOELEMENTSADDED = .FALSE.
            IRHS = ketshell(C)%IODelms(nD)
            DoINT(IRHS) = .TRUE.
         ENDDO
         IF(NOELEMENTSADDED)EXIT !NO ELEMENTS ADDED IN D LOOP
      ENDDO
      IF(INPUT%sameLHSaos)then
         !Add to DoINT from the BC element of the density matrix
         DO nC=1,ML(B)%DIM
            C = ML(B)%Belms(nC)
            NOELEMENTSADDED = .TRUE.
            DO nD=1,ketshell(C)%DIM
               D=ketshell(C)%Belms(nD)
               IF(RED_DMAT_RHS(B,C)*RED_GAB_RHS(A,B)&
                    &*RED_GAB_RHS(C,D) .LE. CS_THRESHOLD ) EXIT
               NOELEMENTSADDED = .FALSE.
               IRHS = ketshell(C)%IODelms(nD)
               DoINT(IRHS)=.TRUE.
            ENDDO
            IF(NOELEMENTSADDED)EXIT !NO ELEMENTS ADDED IN D LOOP
         ENDDO
      ENDIF
      CALL BUILD_LIST(LUPRI,LIST,DoINT,LISTSIZE,dim3*dim4,SIZEOFDOINT)
      DO RHSINDEX=1,LISTSIZE
         IRHS=LIST(RHSINDEX)
         iPassType = ODpassesIndex(IRHS)
         IF(DALINK)THEN
            !EXTRA DENSITY ACCELERATED SCREENING
            C=batchindex1(IRHS)
            D=batchindex2(IRHS)
            DMATELM1=RED_DMAT_RHS(A,C)*RED_DMAT_LHS(B,D)
            DMATELM2=RED_DMAT_RHS(B,C)*RED_DMAT_LHS(A,D)
            MAXDMAT=MAX(DMATELM1,DMATELM2)
            IF(MAXDMAT*RED_GAB_LHS(A,B)*RED_GAB_RHS(C,D) .LT. CS_THRESHOLD ) CYCLE
         ENDIF
         IF (doPasses .AND..NOT. ((ILHS.EQ.IRHS).AND.sameODs)) THEN
            numPasses(iPassType) = numPasses(iPassType) + 1
            CALL AddOverlapToPass(PassQ(iPassType),Q(IRHS),numPasses(iPassType),maxPasses,LUPRI,IPRINT)
            IF (numPasses(iPassType).EQ.maxPasses) THEN
               CALL FinalizePass(PassQ(iPassType),Q(IRHS),maxPasses,LUPRI,IPRINT)
               CALL ExplicitIntegrals(Integral,PQ,P(ILHS),PassQ(iPassType),INPUT,LSOUTPUT,&
                    & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
               CALL DistributeExchangeLinK(PQ,Integral%PQ,DRHS_SYM,PQ%Q%p%totOrbitals,&
                       & PQ%P%p%totOrbitals,INPUT,LSOUTPUT,LUPRI,IPRINT)
               numPasses(iPassType) = 0
            ENDIF
         ELSE
            CALL ExplicitIntegrals(Integral,PQ,P(ILHS),Q(IRHS),INPUT,LSOUTPUT,&
                 & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
            CALL DistributeExchangeLinK(PQ,Integral%PQ,DRHS_SYM,PQ%Q%p%totOrbitals,&
                 & PQ%P%p%totOrbitals,INPUT,LSOUTPUT,LUPRI,IPRINT) 
         ENDIF
      ENDDO
      DO iPassType=1,nPassTypes
         IF(numPasses(iPassType).GT.0) THEN
            IRHS = ODpassesIndex2(iPassType)
            CALL FinalizePass(PassQ(iPassType),Q(IRHS),numPasses(iPassType),LUPRI,IPRINT)
            IRHS=ILHS+1
            !IRHS has been used to build Q(IRHS) which was then added to 
            !PassQ and for all Q(IRHS), IRHS was different from ILHS.
            !which was treated as a special case.  
            !IRHS do no longer have any meaning, but is not allowed be accidentally be equal ILHS
            !because then a triangular loop will be used which PassQ was not built for
            CALL ExplicitIntegrals(Integral,PQ,P(ILHS),PassQ(iPassType),INPUT,LSOUTPUT,&
                 & ILHS,IRHS,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
            CALL DistributeExchangeLinK(PQ,Integral%PQ,DRHS_SYM,PQ%Q%p%totOrbitals,&
                    & PQ%P%p%totOrbitals,INPUT,LSOUTPUT,LUPRI,IPRINT)
            numPasses(iPassType) = 0
         ENDIF
      ENDDO
      !numpass is now 0 for all ipasstype, and ready for next loop
   ENDDO
ENDDO
!$OMP END DO 
call mem_dealloc(LIST)
call mem_dealloc(DoINT)
call mem_dealloc(WORK)
CALL deallocateIntegrals(PQ,Integral)
DO iPassType=1,nPassTypes
   CALL Free_overlap(PassQ(iPassType))
ENDDO
DEALLOCATE(PassQ)
call mem_dealloc(numpasses)
call mem_dealloc(ODpassesIndex2)
!$OMP END PARALLEL

DO IRHS = 1,nRHSoverlaps
  CALL FREE_OVERLAP(Q(IRHS))
ENDDO
DEALLOCATE(Q)
DO ILHS = 1,nLHSoverlaps
  CALL FREE_OVERLAP(P(ILHS))
ENDDO
DEALLOCATE(P)
call mem_dealloc(ODpassesIndex)
call mem_dealloc(nPrimQ)
call mem_dealloc(nPrimP)
CALL freeTUVitem(sharedTUV,Input)

call mem_dealloc(RED_GAB_LHS)
call mem_dealloc(RED_GAB_RHS)
call mem_dealloc(maxLHSGAB)
call mem_dealloc(maxRHSGAB)
DO A=1,dim1
   IF(brashell(A)%DIM.NE.0)THEN
      call linkshell_free(brashell(A))
   ENDIF
ENDDO
DEALLOCATE(brashell)
NULLIFY(brashell)
DO C=1,dim3
   IF(ketshell(C)%DIM.NE.0)THEN
      call linkshell_free(ketshell(C))
   ENDIF
ENDDO
DEALLOCATE(ketshell)
NULLIFY(ketshell)

DO A=1,dim1
   call linkshell_free(ML(A))
ENDDO
DEALLOCATE(ML)
NULLIFY(ML)

IF(DALINK)THEN
   IF(.NOT.SameODs)call mem_dealloc(RED_DMAT_LHS)
   call mem_dealloc(batchindex1)
   call mem_dealloc(batchindex2)
ENDIF
call mem_dealloc(RED_DMAT_RHS)

CALL LS_GETTIM(CPUTIMEEND,WALLTIMEEND)
CPUTIMELINK = CPUTIMEEND-CPUTIMESTART
WALLTIMELINK = WALLTIMEEND-WALLTIMESTART

WRITE(lupri,*)'OVERALL WALL TIMINGS '
CALL ls_TIMTXT('>>>  WALL Time used in Link is             ',WALLTIMELINK,LUPRI)
WRITE(lupri,*)'OVERALL CPU TIMINGS '
CALL ls_TIMTXT('>>>  CPU Time used in Link is             ',CPUTIMELINK,LUPRI)
END SUBROUTINE LINK

!> \brief reduce the input matrix given with dimension of basisfunctions to ODbatch dimension 
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param IUNIT the logical unit number for the output file
!> \param GAB the input matrix with dimension nbast1,nbast2
!> \param nbast1 first dimension of GAB
!> \param nbast2 second dimension of GAB
!> \param RED_GAB the output matrix with reduced dimension dim1,dim2
!> \param dim1 first dimension of RED_GAB
!> \param dim2 second dimension of RED_GAB
!> \param SIDE Label RHS or LHS
!> \param POS logical to specify that GAB is positive and no absolute value is required. 
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
    CALL LSQUIT('Wrong case in reduce_mat',-1)
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

!!$SUBROUTINE REDUCE_DMAT(IUNIT,DMAT,nbast1,nbast2,RED_DMAT,dim1,dim2,INPUT,ndmat,SIDE,SYMMETRIC)
!!$implicit none
!!$TYPE(INTEGRALINPUT)  :: INPUT
!!$INTEGER              :: nbast1,nbast2,dim1,dim2,IUNIT,ndmat
!!$REAL(REALK)          :: DMAT(nbast1,nbast2,ndmat),RED_DMAT(dim1,dim2),TMP
!!$Character*(*)        :: SIDE
!!$logical              :: SYMMETRIC
!!$!
!!$Integer              :: IA,Ja,KA,IB,JB,KB,sA,sB,CA,CB,AOA,AOB,nKA,nKB,idmat
!!$ print*,'DMAT',DMAT(nbast1,nbast2,ndmat),'ndmat',ndmat
!!$SELECT CASE(SIDE)
!!$ CASE('LHS')
!!$!    sameAOs=INPUT%sameLHSaos
!!$    AOA = 1
!!$    AOB = 2
!!$ CASE('RHS')
!!$!    sameAOs=INPUT%sameRHSaos
!!$    AOA = 3
!!$    AOB = 4
!!$ CASE DEFAULT
!!$    WRITE(IUNIT,'(1X,2A)') 'Wrong case in Reduce mat =',SIDE
!!$    CALL LSQUIT('Wrong case in reduce_mat')
!!$ END SELECT
!!$ DO IDMAT = 1,ndmat
!!$    DO IA=1,dim1
!!$       DO IB=1,dim2 
!!$          DO JA=1,INPUT%AO(AOA)%p%BATCH(IA)%nAngmom
!!$             sA=INPUT%AO(AOA)%p%BATCH(IA)%startOrbital(JA)
!!$             DO JB=1,INPUT%AO(AOB)%p%BATCH(IB)%nAngmom
!!$                sB=INPUT%AO(AOB)%p%BATCH(IB)%startOrbital(JB)
!!$                nKA=INPUT%AO(AOA)%p%BATCH(IA)%nOrbComp(JA)
!!$                nKB=INPUT%AO(AOB)%p%BATCH(IB)%nOrbComp(JB)
!!$                DO CA=1,INPUT%AO(AOA)%p%BATCH(IA)%nContracted(JA)
!!$                   DO CB=1,INPUT%AO(AOB)%p%BATCH(IB)%nContracted(JB)
!!$                      DO KA=1,nKA
!!$                         DO KB=1,nKB
!!$                            RED_DMAT(IA,IB)=MAX(RED_DMAT(IA,IB),ABS(DMAT(sA+KA+(CA-1)*nKA-1,sB+KB+(CB-1)*nKB-1,idmat)))
!!$                         ENDDO
!!$                      ENDDO
!!$                   ENDDO
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$ ENDDO
!!$!THE REDUCED MATRIX MUST BE SYMMETRIC FOR LINK TO WORK
!!$IF(.NOT.SYMMETRIC .AND. (dim1.EQ.dim2))THEN
!!$    DO IA=1,dim1
!!$       DO IB=1,IA-1
!!$          TMP = MAX(RED_DMAT(IA,IB),RED_DMAT(IB,IA))
!!$          RED_DMAT(IA,IB) = TMP
!!$          RED_DMAT(IB,IA) = TMP
!!$       ENDDO
!!$    ENDDO
!!$ENDIF
!!$END SUBROUTINE REDUCE_DMAT

!> \brief build list of RHS ODbatches to calculate, based on logical array
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param LUPRI the logical unit number for the output file
!> \param LIST the list to be built
!> \param DOINT the logical array from which the LIST should be built
!> \param LISTSIZE the size of the list to be built
!> \param SIZE the allocated size of LIST and Doint
!> \param ILHS the ODbatch index for the LHS OD batch
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

!> \brief determines the maximum gab matrix element, and array of max gab elements for given bacthindex
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param RED_GAB the gab matrix reduced to odbatches
!> \param dim1 the size of the first dimension
!> \param dim2 the size of the second dimension
!> \param maxGab the array of maximum gab matrix elements
!> \param maxElm the maximum gab matrix element
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

!> \brief determines the significant shell pairs (used in LinK) see JCP 109,1663
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param dim1 the size of the first dimension
!> \param dim2 the size of the second dimension
!> \param RED_GAB_LHS the gab matrix reduced to odbatches
!> \param CS_THRESHOLD the Cauchy-Schwarz screening threshold
!> \param maxRHSElm the maximum rhs gab matrix element
!> \param brashell the significant shell pairs 
!> \param OD the ODbatch
!> \param LHS flag to determine if this is a left hand side (LHS) shell or RHS
!> \param sameRHS flag to determine if the RHS atomic orbitals are the same
!> \param nPrimOD number of primitives for each ODbatch
SUBROUTINE DETERMINE_SHELL_PAIRS(dim1,dim2,RED_GAB_LHS,CS_THRESHOLD,maxRHSELM,brashell,OD,LHS,sameRHS,nPrimOD)
IMPLICIT NONE
integer                 :: dim1,dim2
real(realk)             :: RED_GAB_LHS(dim1,dim2),CS_THRESHOLD,maxRHSELM
TYPE(LINKshell)         :: brashell(:)
TYPE(ODITEM)            :: OD
INTEGER                 :: nPrimOD(:)
logical                 :: LHS,sameRHS
!
real(realk)             :: SORTING(dim2,dim1)
integer                 :: bshell(dim2,dim1),IODindex(dim1,dim2),nIODs
integer                 :: A,startB,I,B,C,number_bshell(dim1),IOD,nB,IRHS
IF(LHS)THEN
   number_bshell = 0
   nIODs = 0
   DO IOD=1,OD%nbatches
      IF(nPrimOD(IOD).GT.0)THEN
         nIODs = nIODs + 1
         A = OD%BATCH(IOD)%IA
         B = OD%BATCH(IOD)%IB
         IF(RED_GAB_LHS(A,B) .GT. CS_THRESHOLD/maxRHSELM)THEN
            number_bshell(A) = number_bshell(A)+1
            bshell(number_bshell(A),A) = B
            SORTING(number_bshell(A),A) = RED_GAB_LHS(A,B)
            IODindex(A,B) = nIODs
         ENDIF
      ENDIF
   ENDDO
   DO A=1,dim1
      IF(number_bshell(A).GT.0)THEN
         I = number_bshell(A)
         call linkshell_init(brashell(A),I)
         CALL Qsort(SORTING(1:I,A),bshell(1:I,A))
         DO C = 1,I
            brashell(A)%Belms(C)=bSHELL(C,A)
            brashell(A)%IODelms(C)=IODindex(A,bSHELL(C,A))
         ENDDO
      ELSE
         brashell(A)%DIM = 0
      ENDIF
   ENDDO
ELSE
   number_bshell = 0
   nIODs = 0
   DO IOD=1,OD%nbatches
      IF(nPrimOD(IOD).GT.0)THEN
         nIODs = nIODs + 1
         A = OD%BATCH(IOD)%IA
         B = OD%BATCH(IOD)%IB
         IF(B.EQ.A)THEN
            IF(RED_GAB_LHS(A,B) .GT. CS_THRESHOLD/maxRHSELM)THEN
               number_bshell(A) = number_bshell(A)+1
               bshell(number_bshell(A),A) = B
               SORTING(number_bshell(A),A) = RED_GAB_LHS(A,B)
               IODindex(A,B) = nIODs
            ENDIF
         ELSE
            IF(RED_GAB_LHS(A,B) .GT. CS_THRESHOLD/maxRHSELM)THEN
               number_bshell(A) = number_bshell(A)+1
               bshell(number_bshell(A),A) = B
               SORTING(number_bshell(A),A) = RED_GAB_LHS(A,B)
               IODindex(A,B) = nIODs
               IF(sameRHS)then
                  number_bshell(B) = number_bshell(B)+1
                  bshell(number_bshell(B),B) = A
                  SORTING(number_bshell(B),B) = RED_GAB_LHS(A,B)
                  IODindex(B,A) = nIODs
               endif
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   DO A=1,dim1
      IF(number_bshell(A).GT.0)THEN
         I = number_bshell(A)
         call linkshell_init(brashell(A),I)
         CALL Qsort(SORTING(1:I,A),bshell(1:I,A))
         DO C = 1,I
            brashell(A)%Belms(C)=bSHELL(C,A)
            brashell(A)%IODelms(C)=IODindex(A,bSHELL(C,A))
         ENDDO
      ELSE
         brashell(A)%DIM = 0
      ENDIF
   ENDDO
ENDIF

END SUBROUTINE DETERMINE_SHELL_PAIRS

!> \brief determines the significant braket pairs (used in LinK)see JCP 109,1663
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param dim1 the number of Atomic orbital shells on center 1
!> \param dim2 the number of Atomic orbital shells on center 2
!> \param dim3 the number of Atomic orbital shells on center 3
!> \param RED_DMAT_RHS the density matrix in shell dimensions
!> \param maxLHSGAB the maximum lhs gab matrix elements for given shell index
!> \param maxRHSGAB the maximum rhs gab matrix elements for given shell index
!> \param CS_THRESHOLD the Cauchy-Schwarz screening threshold
!> \param ML the significant braket pairs 
!> \param sameODs flag to determine if the LHS and RHS atomic orbitals are the same
SUBROUTINE DETERMINE_BRAKET_PAIRS(dim1,dim2,dim3,RED_DMAT_RHS,maxLHSGAB,maxRHSGAB,CS_THRESHOLD,ML,sameODs)
IMPLICIT NONE
integer                 :: dim1,dim2,dim3,NDMAT
real(realk)             :: RED_DMAT_RHS(dim1,dim3),CS_THRESHOLD,maxRHSELM
real(realk)             :: maxLHSGAB(dim1),maxRHSGAB(dim3)
TYPE(LINKshell)         :: ML(dim1)
logical                 :: sameODs
!
real(realk)             :: SORTING(dim3)
integer                 :: bshell(dim3)
integer                 :: A,B,C,I,endC
DO A=1,dim1
   I=0
   endC=dim3
!   IF(sameODs)endC=A
   DO C=1,endC
      IF(RED_DMAT_RHS(A,C)*maxLHSGAB(A)*maxRHSGAB(C).GT. CS_THRESHOLD)THEN
         I=I+1
         bshell(I)=C
         SORTING(I)=RED_DMAT_RHS(A,C)*maxRHSGAB(C)
      ENDIF
   ENDDO
   IF(I .EQ. 0)THEN
      call linkshell_init(ML(A),I)
   ELSE
      ML(A)%DIM=I
      CALL Qsort(SORTING(1:I),bshell(1:I))
      call linkshell_init(ML(A),I)
      DO C=1,I
         ML(A)%Belms(C)=bSHELL(C)
      ENDDO
   ENDIF
ENDDO
END SUBROUTINE DETERMINE_BRAKET_PAIRS

!*************************************************************************
!*
!*                 End of Link
!*
!*************************************************************************

!SUBROUTINE th_init_IntCalcCounter()
!  use thermite_distribute_exchange
!implicit none
!
!call init_IntCalcCounter()
!
!END SUBROUTINE th_init_IntCalcCounter

!SUBROUTINE th_get_IntCalcCounter(nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO)
!  use thermite_distribute_exchange
!implicit none
!integer :: nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO
!
!call get_IntCalcCounter(nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO)
!
!END SUBROUTINE th_get_IntCalcCounter

!> \brief wrapper routine which determines which integral distribution routine to call 
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains the integral which should be distributed
!> \param PQ contains info about the integrand 
!> \param Input contains info about the requested integral evaluation
!> \param Output contains output lstensor 
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DistributeIntegrals(Integral,PQ,Input,Output,LUPRI,IPRINT)
  use thermite_distribute
  use thermite_distribute_exchange
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT,dimQ,dimP,nMAT,nSPHMAT,nDERIV

dimQ=PQ%Q%p%totOrbitals
dimP=PQ%P%p%totOrbitals

IF (Input%CS_int) THEN
   IF(Input%PS_int) THEN
      !Cauchy-Schwarz screening integrals
!$OMP CRITICAL (distributePS)
      call distributelstensorPS(OUTPUT%resultmat2,PQ,Integral%PQ,dimQ,dimP,&
           &Input,Output,LUPRI,IPRINT)
!$OMP END CRITICAL (distributePS)
   ELSE
      !Cauchy-Schwarz screening integrals
!$OMP CRITICAL (distributeCS)
      call distributelstensorCS(OUTPUT%resultmat2,PQ,Integral%PQ,dimQ,dimP,&
           &Input,Output,LUPRI,IPRINT)      
!$OMP END CRITICAL (distributeCS)
   ENDIF
ELSE
   IF (Input%DO_FOCK) THEN
      IF (Input%DO_Coulomb) THEN
         IF (Input%DO_JENGINE) THEN
!$OMP CRITICAL (distributeJengine)
            CALL distributeJengineLstensor(OUTPUT%resultmat2,PQ,Integral%PQ,dimQ,dimP,&
                 &Input,output,LUPRI,IPRINT)
!$OMP END CRITICAL (distributeJengine)
         ELSE
            if(output%docontrule)then
!$OMP CRITICAL (distributePQ)
               CALL GeneraldistributePQ(OUTPUT%contrule(1),OUTPUT%resultmat2,PQ,&
                    &Integral%PQ,dimQ,dimP,Input,output,LUPRI,IPRINT)
!$OMP END CRITICAL (distributePQ)
            endif
         ENDIF
      ENDIF
      IF(Input%DO_Exchange)THEN
         IF(Input%DO_GRADIENT)THEN
            call lstensorExchangeGrad(PQ,Integral%PQ,dimQ,dimP,Input,input%LST_DRHS,&
                 &OUTPUT%resultmat2,LUPRI,IPRINT)
         ELSE
            if(output%docontrule)then
!$OMP CRITICAL (distributePQ)
               CALL GeneraldistributePQ(OUTPUT%contrule(OUTPUT%ncontrule),OUTPUT%resultmat2,&
                    &PQ,Integral%PQ,dimQ,dimP,&
                    &Input,output,LUPRI,IPRINT)
!$OMP END CRITICAL (distributePQ)
            endif
         ENDIF
      ENDIF
   ELSE
      ! No pre-contraction with densities (either full integrals, or contraction
      ! with full integrals)         
      IF(INPUT%DO_MULMOM)THEN
         nMAT =input%nderivQ
         nderiv = INPUT%DERIVORDER
         nSPHMAT=(nderiv+1)*(nderiv+1)
!$OMP CRITICAL (distributeMM)
         CALL PrintMMtoFile(PQ,Integral%PQ,dimQ,dimP,nMAT,nSPHMAT,INPUT,OUTPUT,LUPRI,IPRINT)
!$OMP END CRITICAL (distributeMM)
      ELSE
         if(OUTPUT%dograd)then
            call distributeGrad(OUTPUT%resultmat2,PQ,integral%PQ,dimQ,dimP,Input,LUPRI,IPRINT)
         else
            if(output%docontrule)then
!$OMP CRITICAL (distributePQ)
               CALL GeneraldistributePQ(OUTPUT%contrule(1),OUTPUT%resultmat2,PQ,Integral%PQ,dimQ,dimP,Input,output,LUPRI,IPRINT)
!$OMP END CRITICAL (distributePQ)
            endif
         endif
      ENDIF
   ENDIF
ENDIF

END SUBROUTINE DistributeIntegrals

!> \brief build the integrand structure PQ from P and Q overlap distributions
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param PQ contains info about the integrand, to be built. 
!> \param P the left hand side overlap distribution
!> \param Q the right hand side overlap distribution
!> \param Input contains info about the requested integral evaluation
!> \param ILHS the left hand side ODbatch index 
!> \param IRHS the right hand side ODbatch index 
!> \param IUNIT the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
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
logical :: coulomb,nucrep,mulmom
!
integer :: i
 PQ%Operator = INPUT%operator
 PQ%kinetic  = PQ%Operator(1:7).EQ.'Kinetic'
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
coulomb = (PQ%Operator.EQ.'Coulomb').OR.(PQ%Operator(1:3).EQ.'Erf')
nucrep = PQ%Operator(1:6).EQ.'Nucrep'
mulmom = PQ%Operator(1:6).EQ.'Mulmom'

!CALL SetIntegrand(PQ%iprimP,PQ%iprimQ,PQ%distance,PQ%squaredDistance,&
CALL SetIntegrand(PQ%distance,PQ%squaredDistance,&
     &PQ%exponents,PQ%reducedExponents,coulomb,nucrep,PQ%integralPrefactor,&
     &PQ%nPrimitives,P%exponents,Q%exponents,P%center,Q%center,&
     &P%nPrimitives*P%nPasses,Q%nPrimitives*Q%nPasses,&
     &Q%nPrimitives,Q%nPasses,Input%orderPQ,IUNIT)

IF(INPUT%ATTFACTOR) CALL DSCAL(PQ%nPrimitives,1.d0+INPUT%ATTOMEGA/PI,PQ%integralPrefactor,1)

 IF (IPRINT.GT.15) THEN
    CALL PRINT_OVERLAP(P,IUNIT,IPRINT,'LHS')
    CALL PRINT_OVERLAP(Q,IUNIT,IPRINT,'RHS')
    CALL PRINT_Integrand(PQ,IUNIT,INPUT)
 ENDIF
END SUBROUTINE Build_Integrand

!> \brief set integrand wrapper routine which branch out depending on the requested order of primitive AOs 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param rPQ x,y and z distance between overlap P and Q, (3,nPrimP*nPrimQ)
!> \param squaredDistance the squared distance between overlap P and Q
!> \param exponents the exponent = p + q (p = exponent of LHS overlap P,etc) 
!> \param reducedExponents the reduced exponent p*q/p_q
!> \param coulomb flag which tells if it is a coulomb repulsion integral
!> \param nucrep flag which tells if it is a electron-nuclear attraction integral
!> \param integralprefactor the integralprefactor which depend on the type of integral 
!> \param npq the number of primitives AOs (nprimitives for P * nprimitives for Q)
!> \param pexp the exponents of overlap distribution P
!> \param qexp the exponents of overlap distribution P
!> \param pcent the centers of overlap distribution P
!> \param qcent the centers of overlap distribution P
!> \param np nprimitives*nPasses for P
!> \param nq nprimitives*nPasses for Q
!> \param nPrimQ nprimitives for Q
!> \param nPassQ number of Passes for Q
!> \param orderPQ the requested order of primitives (nPrimP*nPrimQ or nPrimQ*nPrimP)
!> \param IUNIT the logical unit number for the output file
SUBROUTINE SetIntegrand(rPQ,squaredDistance,exponents,&
     &                  reducedExponents,coulomb,nucrep,integralPrefactor,npq,&
     &                  pexps,qexps,pcent,qcent,np,nq,nprimQ,nPassQ,orderPQ,IUNIT)
implicit none
integer     :: np,nq,npq,IUNIT,nprimQ,nPassQ
!integer     :: iprimP(:),iprimQ(:)
real(realk) :: rPQ(:,:),squaredDistance(:)
real(realk) :: exponents(:),reducedExponents(:),integralPrefactor(:)
real(realk) :: pexps(:),pcent(:,:)
real(realk) :: qexps(:),qcent(:,:)
real(realk) :: px,py,pz,qx,qy,qz,pqx,pqy,pqz
logical     :: coulomb,nucrep,orderPQ
IF (orderPQ) THEN
  CALL SetIntegrandPQ(rPQ,squaredDistance,exponents,&
     &                reducedExponents,coulomb,nucrep,integralPrefactor,npq,&
     &                pexps,qexps,pcent,qcent,np,nq,nprimQ,nPassQ,IUNIT)
ELSE
  CALL SetIntegrandQP(rPQ,squaredDistance,exponents,&
     &                reducedExponents,coulomb,nucrep,integralPrefactor,npq,&
     &                pexps,qexps,pcent,qcent,np,nq,nprimQ,nPassQ,IUNIT)
ENDIF
END SUBROUTINE SetIntegrand

!> \brief set the integrand structure PQ from P and Q overlap distributions, using nPrimP*nPrimQ ordering
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param rPQ x,y and z distance between overlap P and Q, (3,nPrimP*nPrimQ)
!> \param squaredDistance the squared distance between overlap P and Q
!> \param exponents the exponent = p + q (p = exponent of LHS overlap P,etc) 
!> \param reducedExponents the reduced exponent p*q/p_q
!> \param coulomb flag which tells if it is a coulomb repulsion integral
!> \param nucrep flag which tells if it is a electron-nuclear attraction integral
!> \param integralprefactor the integralprefactor which depend on the type of integral 
!> \param npq the number of primitives AOs (nprimitives for P * nprimitives for Q)
!> \param pexp the exponents of overlap distribution P
!> \param qexp the exponents of overlap distribution P
!> \param pcent the centers of overlap distribution P
!> \param qcent the centers of overlap distribution P
!> \param np nprimitives*nPasses for P
!> \param nq nprimitives*nPasses for Q
!> \param nPrimQ nprimitives for Q
!> \param nPassQ number of Passes for Q
!> \param IUNIT the logical unit number for the output file
SUBROUTINE SetIntegrandPQ(rPQ,squaredDistance,exponents,&
     &                    reducedExponents,coulomb,nucrep,integralPrefactor,npq,&
     &                    pexps,qexps,pcent,qcent,np,nq,nprimQ,nPassQ,IUNIT)
implicit none
integer     :: np,nq,npq,IUNIT,nprimQ,nPassQ
real(realk) :: rPQ(3,npq),squaredDistance(npq)
real(realk) :: exponents(npq),reducedExponents(npq),integralPrefactor(npq)
real(realk) :: pexps(np),pcent(3,np)
real(realk) :: qexps(nq),qcent(3,nq)
real(realk) :: px,py,pz,qx,qy,qz,pqx,pqy,pqz
logical     :: nucrep,mulmom
!
integer :: ipq,ip,iq,idir,startq,ipassQ,tmpipq
real(realk) :: pqdist(3),pqdist2,d2,p,q,p_q
Real(realk), parameter :: PI=3.14159265358979323846D0, Nill = 0.0d0, OneHalf=1.5d0
Real(realk), parameter :: Two = 2.0d0, TwoHalf=2.5d0 
Real(realk), parameter :: PIFAC = 34.986836655249725D0 !Two*PI**TwoHalf
Real(realk), parameter :: TWOPI = 6.2831853071795862D0 
!
logical :: coulomb

!Primitives ordered according to (P,Q)
!WE BUILD THE PART WHICH IS THE SAME FOR ALL PASSES
ipq = 0
DO iq=1, nPrimQ
   q  = qexps(iq)
   DO ip=1, np
      ipq = ipq+1
      p  = pexps(ip)
      p_q = p + q
      exponents(ipq)          = p_q
      reducedExponents(ipq)   = p*q/p_q
      IF (coulomb) THEN
         integralPrefactor(ipq) = PIFAC/(p*q*SQRT(p_q))
      ELSEIF(nucrep)THEN
         integralPrefactor(ipq) = TWOPI/p
      ELSE
         integralPrefactor(ipq) = (PI/p_q)**(OneHalf)
      ENDIF
   ENDDO
ENDDO
!COPY TO ALL PASSES
DO ipassq=2,nPassQ
   tmpipq = 0
   DO iq=1,nPrimQ
      q  = qexps(iq)
      DO ip=1, np
         tmpipq = tmpipq+1
         ipq = ipq+1
         exponents(ipq)         = exponents(tmpipq)
         reducedExponents(ipq)  = reducedExponents(tmpipq)
         integralPrefactor(ipq) = integralPrefactor(tmpipq)
      ENDDO
   ENDDO
ENDDO
!BUILD THE PART WHICH IS UNIQUE FOR EACH OVERLAP
ipq = 0
DO iq=1, nq
   qx = qcent(1,iq)
   qy = qcent(2,iq)
   qz = qcent(3,iq)
   DO ip=1, np
      ipq = ipq+1
!      iprimP(ipq) = ip
!      iprimQ(ipq) = iq
      pqx = pcent(1,ip) - qx
      pqy = pcent(2,ip) - qy
      pqz = pcent(3,ip) - qz
      rPQ(1,ipq) = pqx
      rPQ(2,ipq) = pqy
      rPQ(3,ipq) = pqz
      squaredDistance(ipq) = pqx*pqx+pqy*pqy+pqz*pqz      
   ENDDO
ENDDO

END SUBROUTINE SetIntegrandPQ

!> \brief set the integrand structure PQ from P and Q overlap distributions, using nPrimQ*nPrimP ordering
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param rPQ x,y and z distance between overlap P and Q, (3,nPrimP*nPrimQ)
!> \param squaredDistance the squared distance between overlap P and Q
!> \param exponents the exponent = p + q (p = exponent of LHS overlap P,etc) 
!> \param reducedExponents the reduced exponent p*q/p_q
!> \param coulomb flag which tells if it is a coulomb repulsion integral
!> \param nucrep flag which tells if it is a electron-nuclear attraction integral
!> \param integralprefactor the integralprefactor which depend on the type of integral 
!> \param npq the number of primitives AOs (nprimitives for P * nprimitives for Q)
!> \param pexp the exponents of overlap distribution P
!> \param qexp the exponents of overlap distribution P
!> \param pcent the centers of overlap distribution P
!> \param qcent the centers of overlap distribution P
!> \param np nprimitives*nPasses for P
!> \param nq nprimitives*nPasses for Q
!> \param nPrimQ nprimitives for Q
!> \param nPassQ number of Passes for Q
!> \param IUNIT the logical unit number for the output file
SUBROUTINE SetIntegrandQP(rPQ,squaredDistance,exponents,&
     &                    reducedExponents,coulomb,nucrep,integralPrefactor,npq,&
     &                    pexps,qexps,pcent,qcent,np,nq,nprimQ,nPassQ,IUNIT)
implicit none
integer     :: np,nq,npq,IUNIT,nprimQ,nPassQ
!integer     :: iprimP(npq),iprimQ(npq)
real(realk) :: rPQ(npq,3),squaredDistance(npq)
real(realk) :: exponents(npq),reducedExponents(npq),integralPrefactor(npq)
real(realk) :: pexps(np),pcent(3,np)
real(realk) :: qexps(nq),qcent(3,nq)
real(realk) :: px,py,pz,qx,qy,qz,pqx,pqy,pqz
logical     :: nucrep
!
integer :: ipq,ip,iq,idir,startq,ipassQ,tmpipq
real(realk) :: pqdist(3),pqdist2,d2,p,q,p_q,pq_inv
Real(realk), parameter :: PI=3.14159265358979323846D0, Nill = 0.0d0, OneHalf=1.5d0
Real(realk), parameter :: Two = 2.0d0, TwoHalf=2.5d0 
Real(realk), parameter :: PIFAC = 34.986836655249725D0 !Two*PI**TwoHalf
Real(realk), parameter :: TWOPI = 6.2831853071795862D0 
!
logical :: coulomb

!Primitives ordered according to (Q,P)
ipq = 0
DO ip=1, np
   p  = pexps(ip)
   !WE BUILD THE PART WHICH IS THE SAME FOR ALL PASSES
   DO iq=1, nPrimQ
      q  = qexps(iq)
      p_q = p + q
      exponents(ipq+iq) = p_q  
      reducedExponents(ipq+iq) = p*q/p_q
      IF (coulomb) THEN
         integralPrefactor(ipq+iq) = PIFAC/(p*q*SQRT(p_q))
      ELSEIF(nucrep)THEN
        integralPrefactor(ipq+iq) = TWOPI/p
      ELSE
        integralPrefactor(ipq+iq) = (PI/p_q)**(OneHalf)
      ENDIF
   ENDDO   
   !COPY TO ALL PASSES
   tmpipq = ipq
   ipq = ipq+nPrimQ
   do ipassQ=2,npassQ
      DO iq=1, nPrimQ
         exponents(ipq+iq) = exponents(tmpipq+iq)
      ENDDO
      ipq = ipq+nPrimQ
   ENDDO
   ipq=tmpipq
   ipq = ipq+nPrimQ
   do ipassQ=2,npassQ
      DO iq=1, nPrimQ
         reducedExponents(ipq+iq) = reducedExponents(tmpipq+iq) 
      ENDDO
      ipq = ipq+nPrimQ
   ENDDO
   ipq=tmpipq
   ipq = ipq+nPrimQ
   do ipassQ=2,npassQ
      DO iq=1, nPrimQ
         integralPrefactor(ipq+iq) = integralPrefactor(tmpipq+iq)
      ENDDO
      ipq = ipq+nPrimQ
   ENDDO
ENDDO
!BUILD THE PART WHICH IS UNIQUE FOR EACH OVERLAP
ipq = 0
DO ip=1, np
   px = pcent(1,ip)
   py = pcent(2,ip)
   pz = pcent(3,ip)
   DO iq=1, nq
      !    TEST here
      ipq = ipq+1
!      iprimP(ipq) = ip
!      iprimQ(ipq) = iq
      pqx = px - qcent(1,iq)
      pqy = py - qcent(2,iq)
      pqz = pz - qcent(3,iq)

      rPQ(ipq,1) = pqx
      rPQ(ipq,2) = pqy
      rPQ(ipq,3) = pqz
       
      squaredDistance(ipq) = pqx*pqx+pqy*pqy+pqz*pqz
   ENDDO
ENDDO

END SUBROUTINE SetIntegrandQP

!> \brief Print the integrand structure PQ
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param PQ contains info about the integrand, to be built. 
!> \param IUNIT the logical unit number for the output file
!> \param INPUT contains info about the requested integral evaluation
SUBROUTINE PRINT_Integrand(PQ,IUNIT,INPUT)
TYPE(Integrand)       :: PQ
Integer               :: IUNIT
TYPE(IntegralInput)   :: INPUT

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
WRITE(IUNIT,'(3X,A)')    '----------------------------------------------------------------------------------'
!WRITE(IUNIT,'(3X,A)')    '  prim  iP  iQ   X_PQ    Y_PQ    Z_PQ   R_PQ^2     exp     redExp     intPre     '
WRITE(IUNIT,'(3X,A)')    '  prim   X_PQ    Y_PQ    Z_PQ   R_PQ^2     exp     redExp     intPre     '
IF(Input%orderPQ)THEN
   DO i=1,PQ%nPrimitives
!      WRITE(IUNIT,'(5X,3I4,4F8.4,2F10.3,F12.7)') i,PQ%iprimP(i),PQ%iprimQ(i),&
      WRITE(IUNIT,'(5X,1I4,4F8.4,2F10.3,F12.7)') i,&
           & PQ%distance(1,i),PQ%distance(2,i),PQ%distance(3,i),&
           & PQ%squaredDistance(i),PQ%exponents(i),PQ%reducedExponents(i),PQ%integralPrefactor(i)
   ENDDO
ELSE
   DO i=1,PQ%nPrimitives
!      WRITE(IUNIT,'(5X,3I4,7F8.4)') i,PQ%iprimP(i),PQ%iprimQ(i),&
      WRITE(IUNIT,'(5X,I4,7F8.4)') i,&
           & PQ%distance(i,1),PQ%distance(i,2),PQ%distance(i,3),&
           & PQ%squaredDistance(i),PQ%exponents(i),PQ%reducedExponents(i),PQ%integralPrefactor(i)
   ENDDO
ENDIF
WRITE(IUNIT,'(3X,A)')    '----------------------------------------------------------------------------------'

END SUBROUTINE PRINT_Integrand

!> \brief contract the LHS overlap P (primitive to contracted functions, Ecoeff contraction,..)
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built. 
!> \param INPUT contains info about the requested integral evaluation
!> \param WORK array to store temporary things
!> \param LWORK last used index of WORK array
!> \param WORKLENGTH the allocated size of the WORK array
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Contract_P(INTEGRAL,PQ,input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: INTEGRAL
TYPE(IntegralInput):: INPUT
TYPE(Integrand)    :: PQ
Integer            :: LUPRI,IPRINT,WORKLENGTH,LWORK 
Real(REALK)        :: WORK(WORKLENGTH)
!
Integer  :: iAngmom,startMomQ,iOrbital,ideriv,nderiv
logical  :: orderAP
!
REAL(REALK)          :: CPU1,CPU2,WALL1,WALL2
Real(realk),dimension(:),pointer :: ptemp

!nderiv=INPUT%NderivP

orderAP = .TRUE.

!Swap real pointers
ptemp => Integral%integrals
Integral%integrals => Integral%RTUV
Integral%RTUV      => ptemp
!call swapRealPointers(Integral%integrals,Integral%RTUV)

iOrbital = 1
!
!Loop over angular contributions sharing the set of primitive functions
!DO iDerivP=1,nderivP NOT YET IMPLEMENTED
DO iAngmom=1,PQ%P%p%nAngmom
   ! Distribute Hermite 2-el integrals WTUV to <tuv1|w|tuv2>
   CALL DistributeHermiteP(Integral,PQ,iAngmom,LUPRI,IPRINT)
   CALL ContractEcoeff(Integral,PQ%P%p,1.d0,iAngmom,LUPRI,IPRINT)
   CALL ContractBasis(Integral,PQ%P%p,iAngmom,WORK,LWORK,WORKLENGTH,orderAP,LUPRI,IPRINT)
   IF (PQ%P%p%derivOrder.GT.0) CALL extractDifferentiated(Integral,PQ%P%p,iAngmom,orderAP,LUPRI,IPRINT)
   CALL SphericalTransform(Integral,PQ%P%p,iAngmom,WORK,LWORK,WORKLENGTH,orderAP,LUPRI,IPRINT)
   CALL AddToPQ(Integral,PQ,iAngmom,iOrbital,LUPRI,IPRINT)
   iOrbital = iOrbital + PQ%P%p%nOrbitals(iAngmom)*PQ%P%p%nDerivComp
ENDDO
!ENDDO

END SUBROUTINE Contract_P

!> \brief Extracts the directional derivative components from the Hermite primitive integrals components
!> \author S. Reine
!> \date 2010-03-07
!> \param integral Contains the information about the integrals
!> \param P Contains the information about the product overlap 
!> \param iAngmom The angular component of the overlap
!> \param orderAP True if the angular components come before the primitives
!> \param LUPRI Default output pint unit
!> \param IPRINT Print level (0 no output - high lots of output)
SUBROUTINE extractDifferentiated(Integral,P,iAngmom,orderAP,LUPRI,IPRINT)
implicit none
TYPE(Integralitem),intent(INOUT) :: INTEGRAL
TYPE(Overlap),intent(IN)         :: P
Logical,intent(IN)               :: orderAP
Integer,intent(IN)               :: iAngmom,LUPRI,IPRINT
!
Integer :: iA1,iA2,l1,l2,ijk,lm,ijkdiff,ijkcart,dim1,ijk1,ijk2
Real(realk),dimension(:),pointer :: ptemp
!
iA1 = P%indexAng1(iangmom)
iA2 = P%indexAng2(iangmom)
l1 = P%orbital1%angmom(iA1)
l2 = P%orbital2%angmom(iA2)
CALL GET_IJK(l1,l2,ijk1,ijk2,lm,ijkdiff,P%sphericalEcoeff,P%derivOrder,P%single)
ijk = ijk1*ijk2
dim1 = Integral%nOrb*Integral%nPrim
IF (orderAP) THEN
  CALL extractDifferentiated_AP(Integral%IN,INTEGRAL%OUT,l1,l2,ijkdiff,ijk,dim1,&
     & P%nDerivComp,P%derivOrder,P%single,P%orbital1%TYPE_empty,P%orbital2%TYPE_empty,&
     & LUPRI,IPRINT)
ELSE
  CALL LSQUIT('Error in extractDifferentiated. orderAP = .FALSE.',lupri)
ENDIF
!
Integral%nAng = ijk
Integral%nDeriv = P%nDerivComp

!Swap pointers IN and OUT
ptemp => Integral%IN
Integral%IN  => Integral%OUT
Integral%OUT => ptemp
!
END SUBROUTINE extractDifferentiated

!> \brief Extracts the directional derivative components from the Hermite primitive integrals components orderAngPrim
!> \author S. Reine
!> \date 2010-03-07
!> \param HermiteDiff Differentiated integrals in Hermite-primitive basis (including expoentital perfactors)
!> \param DirectionalDiff Extracted Hermite components in the different derivative directions
!> \param l1 Angular momentum of orbital 1
!> \param l2 Angular momentum of orbital 2
!> \param ijkdiff Number of Hermite differentiated ijk-components
!> \param ijk Number of undifferentiated ijk-components
!> \param dim1 The number of contracted functions multiplied by the number of orbitals (of the other electron)
!> \param derivComp The number of directional derivative components
!> \param derivOrder The order of differentiation (of current electron)
!> \param single True if only one AO
!> \param empty1 True if orbital 1 is type_empty
!> \param empty2 True if orbital 2 is type_empty
!> \param LUPRI Default output pint unit
!> \param IPRINT Print level (0 no output - high lots of output)
SUBROUTINE extractDifferentiated_AP(HermiteDiff,DirectionalDiff,l1,l2,ijkdiff,ijk,dim1,&
     &                              derivComp,derivOrder,single,empty1,empty2,LUPRI,IPRINT)
implicit none
Integer,intent(IN)      :: l1,l2,ijkdiff,dim1,ijk,derivComp,derivOrder
Real(realk),intent(IN)  :: HermiteDiff(ijkdiff,dim1)
Real(realk),intent(OUT) :: DirectionalDiff(ijk,derivComp,dim1)
Logical,intent(IN)      :: single,empty1,empty2
Integer,intent(IN)      :: LUPRI,IPRINT
!
Integer :: increment(2),icent,ncent,iComp
Integer :: iX,iY,iZ,jX,jY,jZ
Integer :: ijkIndex(0:l1+derivOrder,0:l1+derivOrder,0:l1+derivOrder,0:l2+derivOrder,0:l2+derivOrder,0:l2+derivOrder)
Integer :: i,iDiff,iDim,P1,iP1,jP1,kP1,P2,iP2,jP2,kP2

ncent = 1 + derivOrder
IF (single) ncent = 1

iComp = 0
iDiff=0
DO icent=1,ncent
  increment(1) = derivOrder+1-icent
  increment(2) = icent - 1
  IF (single) THEN
    IF (empty1) THEN
      increment(1) = 0
      increment(2) = derivOrder
    ELSEIF (empty2) THEN
      increment(1) = derivOrder
      increment(2) = 0
    ELSE
      CALL LSQUIT('Error in extractDifferentiated_AP. single and not empty1 or empty2',-1)
    ENDIF
  ENDIF

  P2 = l2+increment(2)
  P1 = l1+increment(1)
! Set up differentiated ijk index
  DO iP2=P2,0,-1
    DO jP2=P2-iP2,0,-1
      kP2 = P2-iP2-jP2
      DO iP1=P1,0,-1
        DO jP1=P1-iP1,0,-1
          kP1=P1-iP1-jP1
          iDiff = iDiff+1
          ijkIndex(iP1,jP1,kP1,iP2,jP2,kP2) = iDiff
        ENDDO
      ENDDO
    ENDDO
  ENDDO

! Loop over the three Cartesian directions for both orbitals
  DO jX=increment(2),0,-1
    DO jY=increment(2)-jX,0,-1
      jZ=increment(2)-jX-jY
      DO iX=increment(1),0,-1
        DO iY=increment(1)-iX,0,-1
          iZ = increment(1)-iX-iY
          iComp = iComp+1
!         Regular ijk-loop
          i = 0
          DO iP2=l2,0,-1
            DO jP2=l2-iP2,0,-1
              kP2 = l2-iP2-jP2
              DO iP1=l1,0,-1
                DO jP1=l1-iP1,0,-1
                  kP1=l1-iP1-jP1
                  i=i+1
                  iDiff = ijkIndex(iP1+iX,jP1+iY,kP1+iZ,iP2+jX,jP2+jY,kP2+jZ)
                  DO iDim=1,dim1
                    DirectionalDiff(i,iComp,iDim) = HermiteDiff(iDiff,iDim)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
!         Regular ijk-loop ends here
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
IF (IPRINT.GT.20) THEN
  CALL PrintTensor(DirectionalDiff,'Directional deriv.  ',ijk,derivComp,dim1,&
     &lupri,'ijk   ','deriv ','dim1  ',3)
ENDIF
!
END SUBROUTINE extractDifferentiated_AP

!> \brief wrapper routine to add the calculated integral into integral%PQ for storage until distribution
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built. 
!> \param iAngmom the index of angular momentum 
!> \param iOrbital the current index in the integral%PQ array  
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddToPQ(Integral,PQ,iAngmom,iOrbital,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: iAngmom,iOrbital,LUPRI,IPRINT
!
Integer :: nQ,nAng,nCont,nDerivP,endOrb,nP

nQ      = PQ%Q%p%totOrbitals
nP      = PQ%P%p%totOrbitals
nAng    = PQ%P%p%nOrbComp(iAngmom)
nCont   = PQ%P%p%nContracted(iAngmom)*PQ%P%p%nPasses
nDerivP = PQ%P%p%nDerivComp
CALL AddToPQ1(Integral%PQ,Integral%IN,iOrbital,nP,nQ,nAng,nCont,nDerivP,LUPRI,IPRINT)
END SUBROUTINE AddToPQ

!> \brief add the calculated integral into integral%PQ for storage until distribution
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param PQ the array to store the integrals 
!> \param addPQ the array to add to PQ
!> \param iOrbital the current index in the PQ array  
!> \param nP the number of AO functions on the LHS overlap P
!> \param nQ the number of AO functions on the RHS overlap Q
!> \param nAng the number of angular momentum on P
!> \param nCont the number of contracted functions on P
!> \param nDerivP the number of derivatives for LHS overlap P
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddToPQ1(PQ,AddPQ,iOrbital,nP,nQ,nAng,nCont,nDerivP,LUPRI,IPRINT)
implicit none
Real(realk) :: PQ(nQ,nP),AddPQ(nAng,nDerivP,nQ,nCont)
Integer     :: iOrbital,nP,nQ,nAng,nCont,nDerivP,LUPRI,IPRINT
!
Integer     :: iQ,iAng,iCont,iP,iDerivP,nOrbP
!
nOrbP = nCont*nAng
!
DO iDerivP=1,nDerivP
  iP = iOrbital
  DO iCont=1,nCont
    DO iAng=1,nAng
      DO iQ=1,nQ
        PQ(iQ,iP+(iDerivP-1)*nOrbP) = AddPQ(iAng,iDerivP,iQ,iCont)
      ENDDO
      iP=iP+1
    ENDDO
  ENDDO
ENDDO
CALL PrintPQ(PQ,nQ,nP,IOrbital,nAng,nCont,nDerivP,LUPRI,IPRINT)
END SUBROUTINE AddToPQ1

!> \brief Print the PQ array in integral%PQ
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param PQ the array to store the integrals, until distribution 
!> \param nQ the total number of AO functions on the RHS overlap Q
!> \param nP the total number of AO functions on the LHS overlap P
!> \param iOrbital the current index in the PQ array  
!> \param nAng the number of angular momentum on P
!> \param nCont the number of contracted functions on P
!> \param nDerivP the number of derivatives for LHS overlap P
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE PrintPQ(PQ,nQ,nP,iOrbital,nAng,nCont,nDerivP,LUPRI,IPRINT)
implicit none
Integer     :: nQ,nAng,nCont,LUPRI,IPRINT,iorbital,nP,nDerivP
Real(realk) :: PQ(nQ,nP)
!
Integer     :: iP,iQ,iAng,iCont,iDerivP,nOrbP
!
IF (IPRINT.GT.5) THEN
  CALL LSHEADER(LUPRI,'PQ contribution')
  nOrbP = nCont*nAng
  DO iDerivP=1,nDerivP
    IF (nDerivP.GT.1) WRITE(LUPRI,'(3X,A,I3)') 'Derivative component iDerivP =',iDerivP
    iP = iOrbital
    DO iCont=1,nCont
      DO iAng=1,nAng
        WRITE(LUPRI,'(5X,A,I3,A,I3)') 'iCont =',iCont,' iAng =',iAng
        WRITE(LUPRI,'(5X,F16.9)') (PQ(iQ,iP+(iDerivP-1)*nOrbP),iQ=1,nQ)
        iP = iP + 1
      ENDDO
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintPQ

!> \brief wrapper to distribute primitive hermite integrals to integral%IN 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> distribute primitive hermite integrals to integral%IN 
!> after the contraction of the RHS overlap Q have been performed 
!>
!> \param integral the storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built
!> \param iAngmom the index of angular momentum
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DistributeHermiteP(Integral,PQ,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: LUPRI,IPRINT,iAngmom
!
Integer :: ntuvP,nP,nOrbQ,startOrbQ,endOrbQ
Integer :: startP,fullSP,endP,ioffP,fullOP,jP,tP,uP,vP,ituvP,ifullP,iOrbQ,iPrimP
Integer :: ifullp1,ifullp2,ifullp3,end2P,der
!
der = PQ%P%p%derivOrder

startP = 0
IF (PQ%P%p%type_hermite_single) startP = PQ%P%p%angmom(iAngmom) + der

fullSP   = PQ%P%p%startAngmom
endP     = PQ%P%p%angmom(iAngmom) + der
nP       = PQ%P%p%nPrimitives
ioffP    = startP*(startP+1)*(startP+2)/6
fullOP   = fullSP*(fullSP+1)*(fullSP+2)/6
ntuvP    = (endP+1)*(endP+2)*(endP+3)/6 - ioffP
nOrbQ    = PQ%Q%p%totOrbitals
Integral%nAng  = ntuvP
Integral%nPrim = np
Integral%nOrb  = nOrbQ
Integral%nDeriv  = 1

IF(PQ%kinetic)THEN
   CALL DistributeHermiteP_kinetic(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,&
     &                             startP,endP,nP,ioffP,fullOP,PQ%P%p%nTUV,ntuvP,nOrbQ,der.EQ.1,LUPRI,IPRINT)
ELSE
   CALL DistributeHermiteP_regular(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,&
     &                             startP,endP,nP,ioffP,fullOP,PQ%P%p%nTUV,ntuvP,nOrbQ,LUPRI,IPRINT)
ENDIF

IF (IPRINT.GT.20) THEN
   CALL PrintTensor(Integral%IN,'DistributeHermiteP  ',&
        &ntuvP,nOrbQ,nP,Lupri,'iTUV  ','iOrb  ','Prim  ',3)
ENDIF

END SUBROUTINE DistributeHermiteP

!> \brief distribute primitive hermite integrals to integral%IN, the default worker routine 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> distribute primitive hermite integrals to integral%IN 
!> after the contraction of the RHS overlap Q have been performed 
!>
!> \param integralIN the storage of the primitive hermite integrals
!> \param TUVQ the intermidiate integrals (at this point the RHS contraction have been performed)
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param startp start angular momentum  on overlap P
!> \param endp end angular momentum on overlap P 
!> \param nP number of primitve(*nPasses) for the P overlap
!> \param ioffP offset for the overlap P tuv index 
!> \param fullOP offset for the full tuv index 
!> \param ntuvfull number of tuv indexes for the full P and Q overlap
!> \param ntuvP  number of tuv indexes for the P overlap 
!> \param nOrbQ the total number of orbitals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DistributeHermiteP_regular(IntegralIN,TUVQ,TUVindex,startP,endP,nP,ioffP,&
     &                                fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT)
implicit none
Real(realk)     :: IntegralIN(ntuvP,nOrbQ,nP)
Real(realk)     :: TUVQ(np,ntuvFull,nOrbQ)
Integer,pointer :: TUVindex(:,:,:)
Integer         :: startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
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

!> \brief distribute primitive hermite integrals to integral%IN, special case
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> distribute primitive hermite integrals to integral%IN 
!> after the contraction of the RHS overlap Q have been performed 
!> this is a special case for the kinetic energy integral, where we apply 
!> nabla to the right hand side overlap Q, and we therefor need to sum the x,y and z componets
!>
!> \param integralIN the storage of the primitive hermite integrals
!> \param TUVQ the intermidiate integrals (at this point the RHS contraction have been performed)
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param startp start angular momentum  on overlap P
!> \param endp end angular momentum on overlap P
!> \param nP number of primitve(*nPasses) for the P overlap
!> \param ioffP offset for the overlap P tuv index 
!> \param fullOP offset for the full tuv index 
!> \param ntuvfull number of tuv indexes for the full P and Q overlap
!> \param ntuvP  number of tuv indexes for the P overlap 
!> \param nOrbQ the total number of orbitals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DistributeHermiteP_kinetic(IntegralIN,TUVQ,TUVindex,startP,endP,nP,ioffP,&
     &                                fullOP,ntuvFull,ntuvP,nOrbQ,deriv,LUPRI,IPRINT)
implicit none
Real(realk)     :: IntegralIN(ntuvP,nOrbQ,nP)
Real(realk)     :: TUVQ(np,ntuvFull,nOrbQ)
Integer,pointer :: TUVindex(:,:,:)
logical         :: deriv
Integer         :: startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
!
Integer :: jP,tP,uP,vP,ituvP,iOrbQ,iPrimP
Integer :: ifullp1,ifullp2,ifullp3
Real(realk), parameter :: Half=0.5d0,One=1.d0
Real(realk) :: factor

IF (deriv) THEN
  factor = One
ELSE
  factor = Half
ENDIF

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
                    &- factor*TUVQ(iPrimP,ifullP1,iOrbQ) &
                    &- factor*TUVQ(iPrimP,ifullP2,iOrbQ) &
                    &- factor*TUVQ(iPrimP,ifullP3,iOrbQ)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE DistributeHermiteP_kinetic

!> \brief contract the RHS overlap Q (primitive to contracted functions, Ecoeff contraction,..)
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built. 
!> \param INPUT contains info about the requested integral evaluation
!> \param WORK array to store temporary things
!> \param LWORK last used index of WORK array
!> \param WORKLENGTH the allocated size of the WORK array
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Contract_Q(INTEGRAL,PQ,Input,WORK,LWORK,WORKLENGTH,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)  :: INTEGRAL
TYPE(Integrand)     :: PQ
TYPE(IntegralInput) :: Input
Integer             :: LUPRI,IPRINT,LWORK,WORKLENGTH
REAL(REALK)         :: WORK(WORKLENGTH)
!
Integer  :: iAngmom,startMomQ,iOrbital,ideriv,nderiv
logical  :: orderAP
!
REAL(REALK)          :: CPU1,CPU2,WALL1,WALL2

orderAP = Input%orderAngPrim

nderiv=Input%NderivQ
!NULLIFY(Integral%TUVQ)
!ALLOCATE(Integral%TUVQ(PQ%P%p%nPrimitives,PQ%P%p%nTUV,PQ%Q%p%totOrbitals))
!
iOrbital = 1
!
!Loop over angular contributions sharing the set of primitive functions
DO iDeriv=1,nderiv
   DO iAngmom=1,PQ%Q%p%nAngmom
      CALL ContractEcoeffQ(Integral,PQ,iAngmom,iDeriv,WORK,WORKLENGTH,orderAP,LUPRI,IPRINT)
      CALL ContractBasis(Integral,PQ%Q%p,iAngmom,WORK,LWORK,WORKLENGTH,orderAP,LUPRI,IPRINT)
      CALL SphericalTransform(Integral,PQ%Q%p,iAngmom,WORK,LWORK,WORKLENGTH,orderAP,LUPRI,IPRINT)
      CALL AddToTUVQ(Integral,PQ,iAngmom,iOrbital,orderAP,LUPRI,IPRINT)
      iOrbital = iOrbital + PQ%Q%p%nOrbitals(iAngmom)*PQ%Q%p%nDerivComp
   ENDDO
ENDDO
    
END SUBROUTINE Contract_Q

!> \brief wrapper contract Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built.
!> \param iAngmom the index of angular momentum
!> \param ideriv the derivative index
!> \param WORK array to store temporary things
!> \param WORKLENGTH the allocated size of the WORK array
!> \param orderAP True if the angular components come before the primitives
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ContractEcoeffQ(Integral,PQ,iAngmom,iDeriv,WORK,WORKLENGTH,orderAP,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)  :: INTEGRAL
TYPE(Integrand)     :: PQ
Integer             :: LUPRI,IPRINT,WORKLENGTH
REAL(REALK)         :: WORK(WORKLENGTH)
Integer             :: iAngmom,iDeriv
logical             :: orderAP
!
IF(orderAP)THEN
!  Distribute Hermite 2-el integrals WTUV to <tuv1|w|tuv2>
   CALL DistributeHermiteQ(Integral,PQ,iAngmom,iDeriv,LUPRI,IPRINT)
   CALL ContractEcoeff(Integral,PQ%Q%p,-1.d0,iAngmom,LUPRI,IPRINT)
ELSE
   CALL DirectcontractEQ(Integral,PQ,PQ%Q%p,-1.d0,iAngmom,ideriv,WORK,WORKLENGTH,LUPRI,IPRINT)
ENDIF
!
END SUBROUTINE ContractEcoeffQ

!> \brief contract Ecoefficients directly for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built.
!> \param Q contains overlap distribution for the right hand side, electron 2
!> \param signQ the sign (-1)
!> \param iAngmom the index of angular momentum
!> \param ideriv the derivative index
!> \param WORK array to store temporary things
!> \param WORKLENGTH the allocated size of the WORK array
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DirectcontractEQ(Integral,PQ,Q,signQ,iAngmom,ideriv,WORK,WORKLENGTH,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
TYPE(Overlap)      :: Q
Integer            :: LUPRI,IPRINT,iAngmom,ideriv,WORKLENGTH
Real(realk)        :: WORK(:)
!
Integer             :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ,nPrimQ,nPassQ
Integer             :: nP,nQ,ioffP,ioffQ,sPQ,ioffPQ,i1,i2,j1,j2,nprim1,nprim2
Integer             :: iPrimPQ,iP,iQ,iTUV,maxQ,minQ,nEFG,l1,l2,ijk1,ijk2,ijk,ijkcart
Integer             :: ntuvP,ntuvQ,startP,endP,startQ,endQ,endPQ,ntuvPQ,startE
Integer             :: start1,start2,nAng,ipassQ,iangQ,iprimq
Real(realk)         :: signQ,DM1 = -1.0d0,sign
!Real(realk),pointer :: intermediate1(:,:,:),intermediate2(:,:)
!Real(realk),pointer :: Ecoeffs(:)

real(4) :: tarr(2),etime,tnew,told,tstart
Real(realk),dimension(:),pointer :: ptemp
Real(realk),parameter :: D2=2.0d0,D1=1.0d0

startP = PQ%P%p%startAngmom
endP   = PQ%P%p%endAngmom
startQ = 0
IF (PQ%Q%p%type_hermite_single) startQ = PQ%Q%p%angmom(iAngmom)
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
Integral%nDeriv = 1
nEFG=Integral%nEFG
IF ((((maxQ.EQ.0).OR.((endP.EQ.0).AND.((maxQ.EQ.endQ).AND.(minQ.EQ.startQ)))).AND.(nEFG.EQ.1)).AND.Q%type_hermite_Single) THEN
   !special case
   ptemp => Integral%IN
   Integral%IN => Integral%RTUV
   Integral%RTUV => ptemp
   l1 = Q%orbital1%angmom(Q%indexAng1(iangmom))
   l2 = Q%orbital2%angmom(Q%indexAng2(iangmom))
   IF(l1.EQ.0 .AND. l2 .EQ. 0)THEN
      sign = D1
   ELSE
      sign = signQ**(l1+l2)
   ENDIF
   nPrim1 = Q%orbital1%nPrimitives
   nPrim2 = Q%orbital2%nPrimitives
   CALL SingleHermiteEcoeff1_PA(Integral%IN,Q,sign,l1,l2,nPrim1,nPrim2,PQ%Q%p%nPasses,ntuvQ,ntuvP*nP,LUPRI,IPRINT)
ELSE
   IF (Q%type_hermite_single) THEN
      l1 = Q%orbital1%angmom(Q%indexAng1(iangmom))
      l2 = Q%orbital2%angmom(Q%indexAng2(iangmom))
      IF(l1.EQ.0 .AND. l2 .EQ. 0)THEN
         sign = D1
      ELSE
         sign = signQ**(l1+l2)
      ENDIF
      nPrim1 = Q%orbital1%nPrimitives
      nPrim2 = Q%orbital2%nPrimitives
      CALL DirectSingleHermiteEcoeff1(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
           &                  nP,nQ,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQ*nEFG,&
           &                  ntuvPQ,ntuvP,ntuvQ,ideriv,Q,sign,l1,l2,nPrim1,nPrim2,PQ%Q%p%nPasses,lupri)
   ELSE
      l1 = Q%orbital1%angmom(Q%indexAng1(iangmom))
      l2 = Q%orbital2%angmom(Q%indexAng2(iangmom))
      CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,Q%sphericalEcoeff,Q%derivOrder,Q%single)
      IF (Q%ETUVisSet) THEN
!         call mem_alloc(Ecoeffs,ijk*Integral%nAng*Integral%nPrim)
         startE = Q%ETUVindex(iAngmom)-1
         nPrimQ = PQ%Q%p%nPrimitives
         nPassQ = PQ%Q%p%nPasses
         nAng = ntuvQ*ijk
         IF(nPrimQ*nPassQ*ntuvQ*ijk .GT. WORKLENGTH)CALL LSQUIT('SOMETHING WRONG NOT ENOUGH MEM i WORK',-1)
         DO ipassQ = 1,nPassQ
            do iAngQ = 1,nAng
               START1 = (ipassQ-1)*nPrimQ+(iAngQ-1)*nPrimQ*nPassQ
               START2 = startE+(iAngQ-1)*nPrimQ+(ipassQ-1)*nPrimQ*nAng
               DO IprimQ = 1,nPrimQ
                  WORK(iprimQ+START1) = Q%ETUV(iprimQ+START2)
               ENDDO
            ENDDO
         ENDDO
         CALL DirectcontractEQ12(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
              &                  nP,nQ,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQ*nEFG,&
              &                  ntuvPQ,ntuvP,ntuvQ,ideriv,WORK,1,ijk,l1,l2,lupri)
!         call mem_dealloc(Ecoeffs)
      ELSE
         IF(ijk*Integral%nAng*Integral%nPrim .GT. WORKLENGTH)CALL LSQUIT('SOMETHING WRONG NOT ENOUGH MEM i WORK',-1)
!         call mem_alloc(Ecoeffs,ijk*Integral%nAng*Integral%nPrim)
         CALL BuildEcoeffTensor_PA(integral%TUV,Q,signQ,WORK,ijk,Integral%nAng,&
              &Integral%nPrim,iAngmom,Q%nPasses,LUPRI,IPRINT)
         CALL DirectcontractEQ12(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
              &                  nP,nQ,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQ*nEFG,&
              &                  ntuvPQ,ntuvP,ntuvQ,ideriv,WORK,1,ijk,l1,l2,lupri)
!         call mem_dealloc(Ecoeffs)
      ENDIF
      Integral%nAng = ijk
   ENDIF
ENDIF
END SUBROUTINE DirectcontractEQ

!> \brief contraction of single hermite Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Herint storage of hermite integrals
!> \param P contains overlap distribution for the left or right hand side
!> \param signP the sign (+1 or -1)
!> \param j1 angular momentum of orbital 1
!> \param j2 angular momentum of orbital 2
!> \param nPrim1 number of primitive for orbital 1
!> \param nPrim2 number of primitive for orbital 2
!> \param nPasses number of passes 
!> \param nAng number of angular components 
!> \param nOrb number of orbitals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SingleHermiteEcoeff1_PA(Herint,P,signP,j1,j2,nPrim1,nPrim2,nPasses,nAng,nOrb,LUPRI,IPRINT)
implicit none
Real(realk)   :: HerInt(nPrim1,nPrim2,nPasses,nAng*nOrb)
TYPE(Overlap) :: P
Integer       :: nAng,nOrb,LUPRI,IPRINT
Integer       :: nPrim1,nPrim2,j1,j2,nPasses
Real(realk)   :: signP
!
Real(realk)   :: pref2,pref12
Real(realk),parameter   :: D1=1.0d0,D2=2.0d0
Integer       :: iPrim1,iPrim2,iAngOrb,iPasses

!outside or inside?  
DO iAngOrb=1,nAng*nOrb
   DO iPasses=1,nPasses
      DO iPrim2=1,nPrim2
         pref2=signP/(D2*P%orbital2%exponents(iPrim2))**j2
         DO iPrim1=1,nPrim1
            pref12=D1/(D2*P%orbital1%exponents(iPrim1))**j1*pref2
            HerInt(iPrim1,iPrim2,iPasses,iAngOrb) = HerInt(iPrim1,iPrim2,iPasses,iAngOrb) * pref12
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SINGLEHERMITEECOEFF1_PA

!> \brief direct contraction of single hermite Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param OUT the output intermidiates
!> \param WTUV2 the primitive hermite integrals 
!> \param nPrim total number of primitive
!> \param nPrimP number of primitive for overlap P
!> \param nPrimQ number of primitive for overlap Q
!> \param startP starting angular momentum for overlap P
!> \param endP end angular momentum for overlap P
!> \param startQ starting angular momentum for overlap Q
!> \param endQ end angular momentum for overlap Q
!> \param ioffP offset for tuvindexing for overlap P
!> \param ioffQ offset for tuvindexing for overlap Q
!> \param ioffPQ offset for tuvindexing for overlap PQ
!> \param nTUVEFGPQ total number of angular components
!> \param ntuvPQ number of angular components for PQ
!> \param ntuvP number of angular components for P
!> \param ntuvQ number of angular components for Q
!> \param ideriv derivative index
!> \param Q the overlap distribution
!> \param signQ the sign for the Q overlap
!> \param l1 angular momentum for orbital 1
!> \param l2 angular momentum for orbital 2
!> \param nPrim1 number of primitives for orbital 1
!> \param nPrim2 number of primitives for orbital 2
!> \param nPasses number of passes
!> \param LUPRI the logical unit number for the output file
SUBROUTINE DirectSingleHermiteEcoeff1(OUT,WTUV2,TUVindex,nPrim,nPrimP,nPrimQ,startP,endP,startQ,endQ,&
     &           ioffP,ioffQ,ioffPQ,nTUVEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,Q,signQ,l1,l2,nPrim1,nPrim2,nPasses,lupri)
implicit none
TYPE(Overlap)      :: Q
Integer         :: nPrim,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ
Integer         :: ntuvEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,nPrimP,nPrimQ,STARTE
Integer         :: j1,j2,nprim1,nprim2,npasses,ipasses,iprim1,iprim2
!Real(realk)     :: WTUV(ntuvEFGPQ,nPrim),tuvTUV(ntuvQ,ntuvP,nPrim)
Real(realk)     :: WTUV2(nPrim,ntuvEFGPQ)!,tuvTUV(ntuvQ,ntuvP,nPrim)
!Real(realk)     :: OUT(ntuvQ,ntuvP,nPrim),signQ,pref2,pref12
Real(realk)     :: OUT(nPrimQ,nPrimP,ntuvP,ntuvQ),signQ,pref2,pref12
Integer,pointer :: TUVindex(:,:,:)
!Real(realk)     :: Ec(nPrimQ)
Integer             :: l1,l2,ijk
!
Integer     :: TUVQPindex(nTUVP*nTUVQ),ituvQP,ntuvp2,ntuvq2
Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ,iOrbP,ijkQ,M
Integer     :: iPrimPQ,iOFF,lupri,iituvP,iituvQ,inTUVPQ,iE,iPrimQ,iprimP,iprimQP
Real(realk),parameter :: D2=2.0d0,D1=1.0d0
Real(realk),pointer  :: TMP(:)
CALL MEM_ALLOC(TMP,nPrim)
iPrimQ=0
DO iPasses=1,nPasses
   DO iPrim2=1,nPrim2
      pref2=signQ/(D2*Q%orbital2%exponents(iPrim2))**l2
      DO iPrim1=1,nPrim1
         iPrimQ = iPrimQ + 1
         pref12=D1/(D2*Q%orbital1%exponents(iPrim1))**l1*pref2
         DO iPrimP=1,nPrimP
            iPrimQP = iPrimQ+(iPrimP-1)*nPrimQ
            TMP(iPrimQP)=pref12
         ENDDO
      ENDDO
   ENDDO
ENDDO

!FIRST WE FIND THE TUVQP INDEX BECAUSE FINDING THINGS IN TUVindex IS EXPENSIVE
iOFF = (ideriv-1)*ntuvPQ-ioffPQ
ituvQP = 0
ituvP = 0
DO jP = startP,endP
   DO tP=jP,0,-1
      DO uP=jP-tP,0,-1
         vP=jP-tP-uP
         ituvP=ituvP+1
         ituvQ=0
         DO jQ = startQ,endQ
            DO tQ=jQ,0,-1
               DO uQ=jQ-tQ,0,-1
                  vQ=jQ-tQ-uQ
                  ituvQ=ituvQ+1
                  iTUVQP = iTUVQP+1 
                  TUVQPindex(iTUVQP) = TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
ntuvp2=ituvp
ntuvq2=ituvq

DO ituvP=1,ntuvP2
   DO ituvQ=1,ntuvQ2
      ituvQP=TUVQPindex(ituvQ+(ituvP-1)*ntuvQ2) 
      DO iPrimP=1,nPrimP
         iPrimQP = (iPrimP-1)*nPrimQ
         DO iPrimQ=1,nPrimQ
            OUT(iPrimQ,iPrimP,ituvP,ituvQ) = WTUV2(iPrimQP+iPrimQ,ituvQP)*TMP(iPrimQP+iPrimQ)
         ENDDO
      ENDDO
   ENDDO
ENDDO
CALL MEM_DEALLOC(TMP)

END SUBROUTINE DirectSingleHermiteEcoeff1

!> \brief direct contraction of hermite Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param OUT the output intermidiates
!> \param WTUV3 the primitive hermite integrals 
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param nPrim total number of primitive
!> \param nPrimP number of primitive for overlap P
!> \param nPrimQ number of primitive for overlap Q
!> \param startP starting angular momentum for overlap P
!> \param endP end angular momentum for overlap P
!> \param startQ starting angular momentum for overlap Q
!> \param endQ end angular momentum for overlap Q
!> \param ioffP offset for tuvindexing for overlap P
!> \param ioffQ offset for tuvindexing for overlap Q
!> \param ioffPQ offset for tuvindexing for overlap PQ
!> \param nTUVEFGPQ total number of angular components
!> \param ntuvPQ number of angular components for PQ
!> \param ntuvP number of angular components for P
!> \param ntuvQ number of angular components for Q
!> \param ideriv derivative index
!> \param Ecoeffs the ecoefficients 
!> \param startE the starting index for Ecoeffs
!> \param ijk cartesian anugular components
!> \param l1 angular momentum for orbital 1
!> \param l2 angular momentum for orbital 2
!> \param LUPRI the logical unit number for the output file
SUBROUTINE DirectcontractEQ12(OUT,WTUV3,TUVindex,nPrim,nPrimP,nPrimQ,startP,endP,startQ,endQ,&
     &                         ioffP,ioffQ,ioffPQ,nTUVEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,Ecoeffs,startE,ijk,l1,l2,lupri)
implicit none
Integer         :: nPrim,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,nPrimQ
Integer         :: ntuvEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,nPrimP,STARTE
Real(realk)     :: WTUV3(nPrim,ntuvEFGPQ)
Real(realk)     :: OUT(nPrimQ,nPrimP,ntuvP,ijk)
Integer,pointer :: TUVindex(:,:,:)
Real(realk)     :: Ecoeffs(:)
Real(realk)     :: Ec(nPrimQ)
Integer             :: l1,l2,ijk
!
Real(realk) :: maxEcont,THRESHOLD,Etmp
Integer     :: startec,TUVQPindex(ntuvP*ntuvQ),ntuvp2,ntuvq2
Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,iOrbP,ijkQ,M,ituvQP
Integer     :: iPrimPQ,iOFF,lupri,iituvP,iituvQ,inTUVPQ,iE,iPrimQ,iprimP,iprimQP
logical     :: addcontribution

!FIRST WE FIND THE TUVQP INDEX BECAUSE FINDING THINGS IN TUVindex IS EXPENSIVE
iOFF = (ideriv-1)*ntuvPQ-ioffPQ
ituvQP = 0
ituvQ = 0
DO jQ = startQ,endQ
   DO tQ=jQ,0,-1
      DO uQ=jQ-tQ,0,-1
         vQ=jQ-tQ-uQ
         ituvQ = ituvQ+1
         ituvP = 0
         DO jP = startP,endP
            DO tP=jP,0,-1
               DO uP=jP-tP,0,-1
                  vP=jP-tP-uP
                  ituvP = ituvP+1
                  iTUVQP = iTUVQP+1 
                  TUVQPindex(iTUVQP) = TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
ntuvp2=ituvp
ntuvq2=ituvq
   
THRESHOLD = 1.0D-15
DO ijkQ = 1, ijk
   !ituvQ = 1
   STARTEc= startE-1+(ijkQ-1)*nPrimQ
   maxEcont = 0.d0
   DO iPrimQ=1,nPrimQ
      Etmp = Ecoeffs(startEc+iPrimQ)
      Ec(iPrimQ) = Etmp
      maxEcont = MAX(maxEcont,ABS(Etmp))
   ENDDO
   IF(maxEcont .GT. THRESHOLD)THEN
      DO ituvP =1,ntuvP2
         ituvQP=TUVQPindex(ituvP)
         iPrimQP = 0
         DO iPrimP=1,nPrimP
            DO iPrimQ=1,nPrimQ
               OUT(iPrimQ,iPrimP,ituvP,ijkQ) = WTUV3(iPrimQP+iPrimQ,ituvQP)*Ec(iPrimQ)
            ENDDO
            iPrimQP = iPrimQP+nPrimQ
         ENDDO
      ENDDO
   ELSE
      CALL LS_DZERO(OUT(:,:,:,ijkQ),nPrimQ*nPrimP*ntuvP2)
   ENDIF

   DO ituvQ = 2,ntuvQ2
      STARTEc= startE-1+(ijkQ-1)*nPrimQ+(ituvQ-1)*ijk*nPrimQ
      maxEcont = 0.d0
      DO iPrimQ=1,nPrimQ
         Etmp = Ecoeffs(startEc+iPrimQ)
         Ec(iPrimQ) = Etmp
         maxEcont = MAX(maxEcont,ABS(Etmp))
      ENDDO
      IF(maxEcont .LT. THRESHOLD)CYCLE
      IOFF = (ituvQ-1)*ntuvP2
      DO ituvP = 1,ntuvP2
         ituvQP=TUVQPindex(ituvP+IOFF)!TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
         iPrimQP=0
         DO iPrimP=1,nPrimP 
            DO iPrimQ=1,nPrimQ
               OUT(iPrimQ,iPrimP,ituvP,ijkQ) = OUT(iPrimQ,iPrimP,ituvP,ijkQ)+WTUV3(iPrimQ+iPrimQP,ituvQP)*Ec(iPrimQ)
            ENDDO
            iPrimQP = iPrimQP+nPrimQ
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE DirectcontractEQ12

!> \brief add intermidiate to Integral%integrals  
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param PQ contains integrand info like reduced exponents, integral prefactor, etc.
!> \param iAngmom the index of angular momentum
!> \param iOrb the current index in the PQ array
!> \param orderAP True if the angular components come before the primitives
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddToTUVQ(Integral,PQ,iAngmom,iOrb,orderAP,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: iAngmom,iOrb,LUPRI,IPRINT 
logical            :: orderAP
!
Integer :: nPrimP,nTUVP,nCompQ,nContQ,endOrb,totOrbQ
!
Integer :: i

nPrimP  = PQ%P%p%nPrimitives*PQ%P%p%nPasses
nTUVP   = PQ%P%p%nTUV
nCompQ  = PQ%Q%p%nOrbComp(iAngmom)
nContQ  = PQ%Q%p%nContracted(iAngmom)*PQ%Q%p%nPasses
endOrb  = iOrb + nCompQ*nContQ - 1
totOrbQ = PQ%Q%p%totOrbitals

IF(orderAP)THEN
   CALL AddToTUVQ_AP(Integral%integrals,Integral%IN,totOrbQ,iOrb,endOrb,nTUVP,nPrimP,&
        &          nCompQ,nContQ,LUPRI,IPRINT)
ELSE
   CALL AddToTUVQ_PA(Integral%integrals,Integral%IN,totOrbQ,iOrb,endOrb,nTUVP,nPrimP,&
        &          nCompQ,nContQ,LUPRI,IPRINT)
ENDIF
END SUBROUTINE AddToTUVQ

!> \brief Perform Spherical Transformation of the intermidiate integrals
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param P the overlap distribution P or Q
!> \param iAngmom the index of angular momentum
!> \param WORK array to store temporary things
!> \param LWORK last used index of WORK array
!> \param WORKLENGTH the allocated size of the WORK array
!> \param orderAP True if the angular components come before the primitives
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SphericalTransform(Integral,P,iAngmom,WORK,LWORK,WORKLENGTH,orderAP,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,iAngmom,LWORK,WORKLENGTH
logical            :: orderAP
!
Real(realk)         :: WORK(WORKLENGTH)
!Real(realk),pointer :: Spherical(:,:)
Integer             :: lm1,lm2,ijk1,ijk2,iA1,iA2,ang1,ang2,lm,ijk
Integer             :: dim1
Logical             :: Sph1,Sph2,spherical
Real(realk),parameter :: D1=1.0d0,D0=0.0d0
Real(realk),dimension(:),pointer :: ptemp
!
! Do not perform a spherical transformation if the E-coefficients are spherical
IF (P%sphericalEcoeff) RETURN

!Consistency testing for derivative case
IF ((P%derivOrder.GT.0).AND.(.NOT.P%orbital1%spherical.OR..NOT.P%orbital2%spherical)) &
     & CALL LSQUIT('Error in SphericalTransform. derivOrder>0 and not spherical!',-1)

iA1  = P%indexAng1(iangmom)
iA2  = P%indexAng2(iangmom)
ang1 = P%orbital1%angmom(iA1)
ang2 = P%orbital2%angmom(iA2)
Sph1 = P%orbital1%spherical.AND.(ang1.GT.1)
Sph2 = P%orbital2%spherical.AND.(ang2.GT.1)
spherical = Sph1.OR.Sph2
IF (spherical) THEN
  ijk1 = (ang1+1)*(ang1+2)/2 
  ijk2 = (ang2+1)*(ang2+2)/2
  CALL GET_IJK(ang1,ang2,lm1,lm2,lm,ijk,.TRUE.,0,P%single)
  
  IF ((LWORK-1+ijk*lm).GT. WORKLENGTH) THEN
     print*,'LWORK     ',LWORK
     print*,'ijk*lm    ',ijk*lm
     CALL LSQUIT('MEM in SphericalTransform',-1)
  ENDIF
  CALL SphericalTransformation(WORK(LWORK:LWORK-1+ijk*lm),&
       & ijk1,ijk2,lm1,lm2,ang1,ang2,P,Integral%IN,Integral%OUT,Integral%nPrim*Integral%nOrb*Integral%nDeriv,&
       & orderAP,lupri,iprint)
  Integral%nAng = lm
  
! Swap pointers IN and OUT
  ptemp => Integral%IN
  Integral%IN  => Integral%OUT
  Integral%OUT => ptemp
ENDIF
END SUBROUTINE SphericalTransform

!> \brief Perform contraction with contraction coefficients to obtain contracted integrals
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param P the overlap distribution P or Q
!> \param iAngmom the index of angular momentum
!> \param WORK array to store temporary things
!> \param LWORK last used index of WORK array
!> \param WORKLENGTH the allocated size of the WORK array
!> \param orderAP True if the angular components come before the primitives
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ContractBasis(Integral,P,iAngmom,WORK,LWORK,WORKLENGTH,orderAP,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,iAngmom,LWORK,WORKLENGTH,IWORK
Real(realk),dimension(WORKLENGTH) :: WORK
logical            :: orderAP
!
Integer             :: iA1,iA2,nCont,nC1,nC2,dim1
Real(realk),parameter  :: D1=1.0d0,D0=0.0d0
Real(realk),dimension(:),pointer :: ptemp
iA1 = P%indexAng1(iangmom)
iA2 = P%indexAng2(iangmom)
nC1 = P%orbital1%nContracted(iA1)
nC2 = P%orbital2%nContracted(iA2)
nCont = nC1*nC2
dim1 = Integral%nAng*Integral%nOrb

IWORK = LWORK+P%nPrimitives*nC1*nC2-1
IF(IWORK .GT. WORKLENGTH)THEN
   print*,'nC1                  ',nC1
   print*,'nC2                  ',nC2
   print*,'WORKLENGTH           ',WORKLENGTH
   print*,'LWORK                ',LWORK
   print*,'IWORK                ',IWORK
   print*,'P%nPrimitives        ',P%nPrimitives
   print*,'nC1*nC2              ',nC1*nC2
   print*,'P%nPrimitives*nC1*nC2',P%nPrimitives*nC1*nC2
   CALL LSQUIT('MEM in ContractBasis',-1)
ENDIF
IF(orderAP)THEN
   CALL ContractBasis_AP(Integral%IN,Integral%OUT,&
        &      WORK(LWORK:IWORK),P,P%nPrimitives,&
        &      nCont,P%nPasses,dim1,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
   IF (IPRINT.GT.50) THEN
      CALL PrintTensor(Integral%OUT,'Contracted          ',&
           &Integral%nAng,Integral%nOrb,nCont*P%nPasses,Lupri,'iAngP ','iOrb  ','iCont ',3)
   ENDIF
ELSE
   CALL ContractBasis_PA(Integral%IN,Integral%OUT,&
        &      WORK(LWORK:IWORK),P,P%nPrimitives,&
        &      nCont,P%nPasses,dim1,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
ENDIF

ptemp => Integral%IN
Integral%IN  => Integral%OUT
Integral%OUT => ptemp
!call swapRealPointers(Integral%IN,Integral%OUT)
Integral%nPrim = nCont*P%nPasses

END SUBROUTINE ContractBasis

!> \brief Perform contraction with Ecoefficients 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param P the overlap distribution P or Q
!> \param signP if Q overlap signP=-1 else signP=1
!> \param iAngmom the index of angular momentum
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ContractEcoeff(Integral,P,signP,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)  :: Integral
TYPE(Overlap)       :: P
Integer             :: LUPRI,IPRINT,iAngmom
Real(realk)         :: signP
!
Real(realk),pointer :: Ecoeffs(:)
Integer             :: l1,l2,ijk1,ijk2,ijk,ijkcart,iPrimP
Real(realk),parameter  :: D1=1.0d0,D0=0.0d0
Real(realk),dimension(:),pointer :: ptemp
!
IF (P%type_hermite_single) THEN
  CALL SingleHermiteEcoeff(Integral%In,P,signP,iAngmom,P%nPasses,Integral%nAng,Integral%nOrb,LUPRI,IPRINT)
ELSE
!  CALL LSQUIT('Only Hermite implemented in ContractEcoeff')
  l1 = P%orbital1%angmom(P%indexAng1(iangmom))
  l2 = P%orbital2%angmom(P%indexAng2(iangmom))
  CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,P%sphericalEcoeff,P%derivOrder,P%single)
  IF (P%ETUVisSet) THEN
!    IF (IPRINT.GT.50) THEN
!       CALL PrintTensor(P%ETUV(P%ETUVindex(iAngmom)),'E-coefficients      ',&
!            &ijk,Integral%nAng,Integral%nPrim,Lupri,&
!            &'ijk   ','tuv   ','prim  ',3)
!    ENDIF
    CALL ContractEcoeff1(Integral%IN,Integral%OUT,P%ETUV,P%ETUVindex(iAngmom), &
       &                 Integral%nOrb,Integral%nAng,ijk,Integral%nPrim)
  ELSE
    call mem_alloc(Ecoeffs,ijk*Integral%nAng*Integral%nPrim)
    CALL BuildEcoeffTensor(integral%TUV,P,signP,Ecoeffs,ijk,ijkcart,Integral%nAng,&
         &Integral%nPrim,.TRUE.,iAngmom,P%nPasses,LUPRI,IPRINT)
    CALL ContractEcoeff1(Integral%IN,Integral%OUT,Ecoeffs,1, &
         &                 Integral%nOrb,Integral%nAng,ijk,Integral%nPrim)
    call mem_dealloc(Ecoeffs)
  ENDIF
  Integral%nAng = ijk
  ptemp => Integral%IN
  Integral%IN  => Integral%OUT
  Integral%OUT => ptemp
!  call swapRealPointers(Integral%IN,Integral%OUT)
ENDIF

IF (IPRINT.GT.50) THEN
   CALL PrintTensor(Integral%IN,'ContractEcoeff      ',&
        &Integral%nAng,Integral%nOrb,Integral%nPrim,Lupri,&
        &'ijk   ','iOrb  ','Prim  ',3)
ENDIF

END SUBROUTINE ContractEcoeff

!> \brief wrapper routine to perform hermite integrals with FTUVs
!> \author S. Reine
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param P the left hand side overlap distribution 
!> \param Q the right hand side overlap distribution 
!> \param ndmat number of density matrices
!> \param orderAngPrim True if the angular components come before the primitives
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ContractFTUV(Integral,P,Q,ndmat,orderAngPrim,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P,Q
Integer            :: LUPRI,IPRINT,ndmat,tuvOffP,startAngP,tuvOffQ,startAngQ,startAngPQ,tuvOffPQ
Logical            :: orderAngPrim
!
startAngP  = P%startAngmom
tuvOffP    = startAngP*(startAngP+1)*(startAngP+2)/6
startAngQ  = Q%startAngmom
tuvOffQ    = startAngQ*(startAngQ+1)*(startAngQ+2)/6
startAngPQ = startAngP + startAngQ
tuvOffPQ   = startAngPQ*(startAngPQ+1)*(startAngPQ+2)/6
IF (orderAngPrim) THEN
  CALL ContractFTUVAP(Integral%integrals,Integral%RTUV,Q%FTUV,Integral%TUV%TUVindex,&
     &           Integral%TUV%Tindex,Integral%TUV%Uindex,Integral%TUV%Vindex,&
     &           Integral%ntuv,Integral%nPrim,P%nTUV,Q%nTUV,tuvOffP,tuvOffQ,tuvOffPQ,&
     &           P%nPrimitives,Q%nPrimitives,ndmat,lupri)
ELSE
  CALL ContractFTUVPA(Integral%integrals,Integral%RTUV,Q%FTUV,Integral%TUV%TUVindex,&
     &           Integral%TUV%Tindex,Integral%TUV%Uindex,Integral%TUV%Vindex,&
     &           Integral%ntuv,Integral%nPrim,P%nTUV,Q%nTUV,tuvOffP,tuvOffQ,tuvOffPQ,&
     &           P%nPrimitives,Q%nPrimitives,ndmat,lupri)
ENDIF

Integral%nAng  = 1
Integral%nPrim = ndmat
Integral%nOrb  = P%nPrimitives*P%nTUV
CALL PrintTUVQ(Integral%integrals,P%nTUV,P%nPrimitives,1,ndmat,LUPRI,IPRINT)

END SUBROUTINE ContractFTUV

!> \brief Perform hermite integrals with FTUVs ordering (nTUVPQ,nPrimPQ)
!> \author S. Reine
!> \date 2010
!>
!> \param TUV the output intermediate integral
!> \param WTUV the primitive hermite integrals 
!> \param FTUV the FTUV batches 
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param Tindex for a given TUV index it gives the T index
!> \param Uindex for a given TUV index it gives the U index
!> \param Vindex for a given TUV index it gives the V index
!> \param ntuvPQ number of angular components for PQ
!> \param nPrimPQ total number of primitive
!> \param ntuvP number of angular components for P
!> \param ntuvQ number of angular components for Q
!> \param tuvoffP offset for tuvindexing for overlap P
!> \param tuvoffQ offset for tuvindexing for overlap Q
!> \param tuvoffPQ offset for tuvindexing for overlap PQ
!> \param nPrimP number of primitive for overlap P
!> \param nPrimQ number of primitive for overlap Q
!> \param nDmat the number of density matrices
!> \param LUPRI the logical unit number for the output file
SUBROUTINE ContractFTUVAP(TUV,WTUV,FTUV,TUVindex,Tindex,Uindex,Vindex,&
     &                    ntuvPQ,nPrimPQ,ntuvP,ntuvQ,tuvOffP,tuvOffQ,tuvOffPQ,&
     &                    nPrimP,nPrimQ,nDmat,lupri)
implicit none
Integer             :: ntuvP,ntuvQ,nPrimP,nPrimQ,ntuvPQ,nPrimPQ,nDmat
Integer             :: tuvOffP,tuvOffQ,tuvOffPQ
Real(realk)         :: TUV(nPrimP,ntuvP,nDmat)
Real(realk),pointer :: FTUV(:,:,:)
Real(realk)         :: WTUV(nTUVPQ,nPrimPQ)
Integer,pointer     :: TUVindex(:,:,:),Tindex(:),Uindex(:),Vindex(:)
!
Integer     :: idmat,tuvP,tP,uP,vP,tuvQ,tQ,uQ,vQ,jQ,tuvPQ
Integer     :: iPrimP,iPrimQ,iPrimPQ,lupri
Real(realk),parameter :: D1=1.0d0,DM1=-1.0d0
Real(realk) :: signQ,fsum
!Simen
Real(realk) :: newfsum(nPrimP)

! replace by loop
DO idmat = 1,ndmat
DO tuvP=1,ntuvP
  tP = Tindex(tuvP+tuvOffP)
  uP = Uindex(tuvP+tuvOffP)
  vP = Vindex(tuvP+tuvOffP)
  DO tuvQ=1,ntuvQ
    tQ = Tindex(tuvQ+tuvOffQ)
    uQ = Uindex(tuvQ+tuvOffQ)
    vQ = Vindex(tuvQ+tuvOffQ)
    jQ = tQ+uQ+vQ
!   signQ = D1
!   IF (mod(jQ,2).EQ.1) signQ = DM1
    tuvPQ = TUVindex(tP+tQ,uP+uQ,vP+vQ)-tuvOffPQ
!simen
    IF (.FALSE.) THEN
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
          TUV(iPrimP,tuvP,idmat) = TUV(iPrimP,tuvP,idmat) + newfsum(iPrimP)
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
      TUV(iPrimP,tuvP,idmat) = TUV(iPrimP,tuvP,idmat) + fsum
    ENDDO
    ENDIF
!simen
  ENDDO
ENDDO
ENDDO
END SUBROUTINE ContractFTUVAP

!> \brief Perform hermite integrals with FTUVs ordering (nPrimPQ,nTUVPQ)
!> \author S. Reine
!> \date 2010
!>
!> \param TUV the output intermediate integral
!> \param WTUV the primitive hermite integrals 
!> \param FTUV the FTUV batches 
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param Tindex for a given TUV index it gives the T index
!> \param Uindex for a given TUV index it gives the U index
!> \param Vindex for a given TUV index it gives the V index
!> \param ntuvPQ number of angular components for PQ
!> \param nPrimPQ total number of primitive
!> \param ntuvP number of angular components for P
!> \param ntuvQ number of angular components for Q
!> \param tuvoffP offset for tuvindexing for overlap P
!> \param tuvoffQ offset for tuvindexing for overlap Q
!> \param tuvoffPQ offset for tuvindexing for overlap PQ
!> \param nPrimP number of primitive for overlap P
!> \param nPrimQ number of primitive for overlap Q
!> \param nDmat the number of density matrices
!> \param LUPRI the logical unit number for the output file
SUBROUTINE ContractFTUVPA(TUV,WTUV,FTUV,TUVindex,Tindex,Uindex,Vindex,&
     &                    ntuvPQ,nPrimPQ,ntuvP,ntuvQ,tuvOffP,tuvOffQ,tuvOffPQ,&
     &                    nPrimP,nPrimQ,nDmat,lupri)
implicit none
Integer             :: ntuvP,ntuvQ,nPrimP,nPrimQ,ntuvPQ,nPrimPQ,nDmat
Integer             :: tuvOffP,tuvOffQ,tuvOffPQ
Real(realk)         :: TUV(nPrimP,ntuvP,nDmat)
Real(realk),pointer :: FTUV(:,:,:)
Real(realk)         :: WTUV(nPrimPQ,nTUVPQ)
Integer,pointer     :: TUVindex(:,:,:),Tindex(:),Uindex(:),Vindex(:)
!
Integer     :: idmat,tuvP,tP,uP,vP,tuvQ,tQ,uQ,vQ,jQ,tuvPQ
Integer     :: iPrimP,iPrimQ,iPrimPQ,lupri
Real(realk) :: fsum

DO idmat = 1,ndmat
DO tuvP=1,ntuvP
  tP = Tindex(tuvP+tuvOffP)
  uP = Uindex(tuvP+tuvOffP)
  vP = Vindex(tuvP+tuvOffP)
  DO tuvQ=1,ntuvQ
    tQ = Tindex(tuvQ+tuvOffQ)
    uQ = Uindex(tuvQ+tuvOffQ)
    vQ = Vindex(tuvQ+tuvOffQ)
    jQ = tQ+uQ+vQ
    tuvPQ = TUVindex(tP+tQ,uP+uQ,vP+vQ)-tuvOffPQ
    iPrimPQ = 1
    DO iPrimP=1,nPrimP
      fsum = WTUV(iPrimPQ,tuvPQ)*FTUV(1,tuvQ,idmat)
      iPrimPQ = iPrimPQ+1
      DO iPrimQ=2,nPrimQ
        fsum = fsum + WTUV(iPrimPQ,tuvPQ)*FTUV(iPrimQ,tuvQ,idmat)
        iPrimPQ = iPrimPQ+1
      ENDDO
      TUV(iPrimP,tuvP,idmat) = TUV(iPrimP,tuvP,idmat) + fsum
    ENDDO
  ENDDO
ENDDO
ENDDO
END SUBROUTINE ContractFTUVPA

!> \brief wrapper routine to Distribute the hermite integrals in order to later use DGEMM for the Ecoefficient contraction. 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built. 
!> \param iAngmom the index of angular momentum
!> \param ideriv derivative index
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
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
Real(realk),dimension(:),pointer :: ptemp

real(4) :: tarr(2),etime,tnew,told,tstart

startP = PQ%P%p%startAngmom
endP   = PQ%P%p%endAngmom
startQ = 0
IF (PQ%Q%p%type_hermite_single) startQ = PQ%Q%p%angmom(iAngmom)
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
Integral%nDeriv  = 1
nEFG=Integral%nEFG
IF (((maxQ.EQ.0).OR.((endP.EQ.0).AND.((maxQ.EQ.endQ).AND.(minQ.EQ.startQ)))).AND.(nEFG.EQ.1)) THEN
  ptemp => Integral%IN
  Integral%IN   => Integral%RTUV
  Integral%RTUV => ptemp
!  call swapRealPointers(Integral%IN,Integral%RTUV)
ELSE
!  CALL DistributeHermiteQ1(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
!     &                     startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQ*nEFG,ntuvPQ,ntuvP,ntuvQ,ideriv,lupri)
  CALL DistributeHermiteQ1(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
     &                     startP,endP,startQ,endQ,ioffPQ,ntuvPQ*nEFG,ntuvPQ,ntuvP,ntuvQ,ideriv,lupri)
ENDIF

IF (IPRINT.GT.50) THEN
   write(lupri,*)'  Number of primitives',nP*nQ
   CALL PrintHermitePQ(Integral%IN,nP*nQ,startP,endP,startQ,endQ,ntuvP,ntuvQ,lupri)
ENDIF

END SUBROUTINE DistributeHermiteQ

!> \brief determine the start and end index of the RHS loop 
!> \author S. Reine
!> \date 2010
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param OD_RHS the ODbatches belonging to the Reft hand side
!> \param ILHS the left hand side 
!> \param start_RHS the start index of the RHS loop
!> \param end_RHS the end index of the RHS loop
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

!> \brief clear the integrals
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Integral contains arrays to store intermidiates and final integrals
SUBROUTINE CLEAR_integral(INTEGRAL)
TYPE(Integralitem)      :: integral
 
!DEALLOCATE(INTEGRAL%Wtuv)
NULLIFY(INTEGRAL%Wtuv)
 
END SUBROUTINE

!> \brief Subroutine that returns A, the transpose of B
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param A the output matrix, to be transpose of B
!> \param B the input matrix
!> \param n1 dimension 1 of the output matrix
!> \param n2 dimension 2 of the output matrix
SUBROUTINE transposition(A,B,n1,n2)
implicit none
integer     :: n1,n2
Real(realk) :: A(n1,n2)
Real(realk) :: B(n2,n1)
A = transpose(B)
END SUBROUTINE transposition

!> \brief determine the maximum distance
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param MAXDIST the maximum distance
!> \param distance the input distance
!> \param nprim the size of distance
SUBROUTINE MAXDISTANCE(MAXDIST,distance,nPrim)
IMPLICIT NONE
INTEGER      :: nPrim,I
REAL(REALK)  :: distance(nPrim),MAXDIST

MAXDIST = sum ( abs ( distance(1:1+(nPrim-1)*1:1) ) )

END SUBROUTINE MAXDISTANCE

!> \brief recursive quick sort routine, sorts largest element first 
!> \author T. Kjaergaard
!> \date 2010
!> \param A the array to be sorted 
!> \param INDEXES an array of indexes which should be reordered in the same was as A
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

!> \brief insertion sort routine, sorts largest element first 
!> \author T. Kjaergaard
!> \date 2010
!> \param A the array to be sorted 
!> \param INDEXES an array of indexes which should be reordered in the same was as A
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

!> \brief the actual sorting routine, used by qsort routine, sorts largest element first 
!> \author T. Kjaergaard
!> \date 2010
!> \param A the array to be sorted 
!> \param INDEXES an array of indexes which should be reordered in the same was as A
!> \param marker the index in the array
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

!> \brief initialise the allocitem 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Alloc the allocitem, which contain info about maximum values required for allocation
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param SIDE either both if both side should be initialised or LHS,RHS if only one side
SUBROUTINE initAlloc(Alloc,LUPRI,IPRINT,SIDE)
implicit none
TYPE(Allocitem)    :: Alloc
Integer            :: LUPRI,IPRINT
Character*(*)      :: Side
IF (SIDE.EQ.'Both') THEN
!  LHS
   Alloc%maxJLHS = 0
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
   Alloc%maxETUVlenLHS = 0
!  RHS
   Alloc%maxJRHS = 0
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
   Alloc%maxETUVlenRHS = 0
ELSE IF (SIDE.EQ.'LHS') THEN
   Alloc%maxJLHS = 0
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
   Alloc%maxETUVlenLHS = 0
ELSE IF (SIDE.EQ.'RHS') THEN
   Alloc%maxJRHS = 0
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
   Alloc%maxETUVlenRHS = 0
ELSE
  WRITE(LUPRI,'(1X,2A)') 'Programming error. Side =',SIDE
  CALL LSQUIT('Programming error. Side not an option in initIntegral',lupri)
ENDIF
END SUBROUTINE initAlloc

!> \brief allocate the integral item
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param PQ contains integrand info, to be allocated, using values in alloc
!> \param Integral contains arrays to store intermidiates and final integrals, to be allocated
!> \param input the integral input specifications, contains all info about what to do
!> \param Alloc the allocitem, which contain info about maximum values required for allocation
!> \param maxpasses the maximum number of passes
!> \param ndmat number of density matrix
SUBROUTINE allocateIntegrals(PQ,Integral,Input,Alloc,maxpasses,ndmat)
  use memory_handling
implicit none
TYPE(Integrand)    :: PQ
TYPE(Integralitem) :: Integral
TYPE(Integralinput):: Input
TYPE(Allocitem)    :: Alloc
Integer            :: maxPasses,ndmat
!
Integer :: maxTUVdim,maxPrim

IF(Alloc%maxPrimRHSpass .NE. 0)THEN
   maxPrim   = Alloc%maxPrimLHS*Alloc%maxPrimRHSpass
ELSE
   maxPrim   = Alloc%maxPrimLHS*Alloc%maxPrimRHS*maxPasses
ENDIF
IF (Input%orderAngPrim) THEN
  call mem_alloc(PQ%distance,3,maxPrim)
ELSE
  call mem_alloc(PQ%distance,maxPrim,3)
ENDIF
call mem_alloc(PQ%squaredDistance,maxPrim)
Call Mem_alloc(PQ%exponents,maxPrim)
Call Mem_alloc(PQ%reducedExponents,maxPrim)
Call Mem_alloc(PQ%integralPrefactor,maxPrim)
!Call Mem_alloc(PQ%iprimP,maxPrim)
!Call Mem_alloc(PQ%iprimQ,maxPrim)

IF(Alloc%maxPrimTUVRHSpass .NE. 0) THEN
   maxTUVdim = Alloc%maxPrimTUVRHSpass*Alloc%maxPrimTUVLHS
ELSE
   maxTUVdim = Alloc%maxPrimTUVRHS*Alloc%maxPrimTUVLHS*maxPasses
ENDIF
CALL MEM_ALLOC(Integral%IN,maxTUVdim+10000)
CALL MEM_ALLOC(Integral%OUT,maxTUVdim+10000)
CALL MEM_ALLOC(Integral%RTUV,maxTUVdim+10000)
CALL MEM_ALLOC(Integral%integrals,maxTUVdim+10000)
CALL MEM_ALLOC(Integral%PQ,Alloc%maxTotOrbRHS*maxpasses*ndmat,Alloc%maxTotOrbLHS+10000)
!$OMP CRITICAL (memory)
mem_allocated_integrand = mem_allocated_integrand + &
     & mem_realsize*size(PQ%distance)  + mem_realsize*size(PQ%squaredDistance) + &
     & mem_realsize*size(PQ%exponents) + mem_realsize*size(PQ%reducedExponents)+ &
     & mem_realsize*size(PQ%integralPrefactor) !+ mem_intsize*size(PQ%iprimP) + &
 !    & mem_intsize*size(PQ%iprimQ)

max_mem_used_integrand = MAX(max_mem_used_integrand, mem_allocated_integrand)

mem_allocated_integralitem = mem_allocated_integralitem + &
     & mem_realsize*size(Integral%IN)  + mem_realsize*size(Integral%OUT) + &
     & mem_realsize*size(Integral%RTUV) + mem_realsize*size(Integral%integrals)+ &
     & mem_realsize*size(Integral%PQ)

max_mem_used_integralitem = MAX(max_mem_used_integralitem, mem_allocated_integralitem)
!$OMP END CRITICAL (memory)

END SUBROUTINE allocateIntegrals

!> \brief deallocate the integral item
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param PQ contains integrand info to be deallocated
!> \param Integral contains arrays to store intermidiates and final integrals, to be deallocated
SUBROUTINE deallocateIntegrals(PQ,Integral)
  use memory_handling
IMPLICIT NONE
TYPE(Integrand)    :: PQ
TYPE(Integralitem) :: Integral
 
!$OMP CRITICAL (memory)
mem_allocated_integrand = mem_allocated_integrand - &
     & mem_realsize*size(PQ%distance)  - mem_realsize*size(PQ%squaredDistance) - &
     & mem_realsize*size(PQ%exponents) - mem_realsize*size(PQ%reducedExponents)- &
     & mem_realsize*size(PQ%integralPrefactor) !- mem_intsize*size(PQ%iprimP) - &
     !& mem_intsize*size(PQ%iprimQ)
if (mem_allocated_integrand < 0) then
   call LSQUIT('Error in deallocateintegrand using mem_allocated_integrand - probably integer overflow!',-1)
endif
mem_allocated_integralitem = mem_allocated_integralitem - &
     & mem_realsize*size(Integral%IN)  - mem_realsize*size(Integral%OUT) - &
     & mem_realsize*size(Integral%RTUV) - mem_realsize*size(Integral%integrals)- &
     & mem_realsize*size(Integral%PQ)
if (mem_allocated_integralitem < 0) then
   call LSQUIT('Error in deallocateintegrals using mem_allocated_integralitem - probably integer overflow!',-1)
endif
!$OMP END CRITICAL (memory)

call mem_dealloc(PQ%distance)
call mem_dealloc(PQ%squaredDistance)
call mem_dealloc(PQ%exponents)
call mem_dealloc(PQ%reducedExponents)
call mem_dealloc(PQ%integralPrefactor)
!call mem_dealloc(PQ%iprimP)
!call mem_dealloc(PQ%iprimQ)
call mem_dealloc(Integral%IN)
call mem_dealloc(Integral%OUT)
call mem_dealloc(Integral%RTUV)
call mem_dealloc(Integral%integrals)
CALL MEM_DEALLOC(Integral%PQ)

END SUBROUTINE deallocateIntegrals

!> \brief calculate maximum values and store them in the alloc item
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Alloc the allocitem, which contain info about maximum values required for allocation
!> \param input the integral input specifications, contains all info about what to do
!> \param ODbat the ODbatch
!> \param side LHS or RHS character 
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param LUPRI the logical unit number for the output file
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
Integer             :: MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK,MAXnTUV,l,nTUV
Integer             :: MAXPRIMTUVIJK,maxA,maxB,start,maxijk2,maxlenETUV
Integer             :: nDer,iDer,ijk1der,ijk2der,ijkcart
LOGICAL             :: LHS,hermiteSingle,single

SELECT CASE(SIDE)
CASE('LHS')
   LHS=.TRUE.
   NDERIV = INPUT%NDERIVP
CASE('RHS')
   LHS=.FALSE.
   NDERIV = INPUT%NDERIVQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_OVERLAP =',SIDE
   CALL LSQUIT('Wrong case in SET_OVERLAP',lupri)
END SELECT

nDer = 0
IF (INPUT%DO_GRADIENT) nDer = INPUT%derivOrder

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
maxlenETUV=0
DO I=1,ODbat%nbatches
   CALL GET_NPRIM_FROM_ODBATCH(nPrim,ODBat%BATCH(I),INPUT,SIDE,lupri)
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
         single = (ODbat%BATCH(I)%AO(1)%p%TYPE_Empty.OR.ODbat%BATCH(I)%AO(2)%p%TYPE_Empty).AND..NOT.&
     &            (ODbat%BATCH(I)%AO(1)%p%TYPE_Empty.AND.ODbat%BATCH(I)%AO(2)%p%TYPE_Empty)
         CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,.FALSE.,nDer,single)
         l  = l1+l2+nDer
         nTUV = (l+1)*(l+2)*(l+3)/6
         maxlenETUV = MAX(maxlenETUV,nTUV*ijk*nPrim)
         maxijk = max(maxijk,ijk)
         maxAngmom = max(maxAngmom,l1+l2)
         totOrbitals = totOrbitals + ODbat%BATCH(I)%AO(1)%p%nOrbitals(i1)*ODbat%BATCH(I)%AO(2)%p%nOrbitals(i2)*NDERIV
      ENDDO
   ENDDO
   MAXPRIMIJK = MAX(MAXPRIMIJK,maxijk*nprim)
   maxTotorb = MAX(maxTotorb,totOrbitals)
   maxijk2 = MAX(maxijk2,maxijk)
   endAngmom  = maxAngmom
   IF(INPUT%operator(1:7).EQ.'Kinetic'.AND.LHS) endAngmom  = maxAngmom + 2
   endAngmom = endAngmom + nDer ! maybe to maxAngmom?
   maxTUV = (endAngmom+1)*(endAngmom+2)*(endAngmom+3)/6
   MAXPRIMTUV = MAX(MAXPRIMTUV,maxTUV*nPrim)
   MAXnTUV=MAX(MAXnTUV,maxTUV)
   MAXPRIMTUVIJK = MAX(MAXPRIMTUVIJK,maxijk*nPrim*maxTUV*ODbat%BATCH(I)%nAngmom) 
ENDDO

IF(LHS)THEN
   Alloc%maxJLHS = maxAngmom 
   Alloc%maxPrimLHS    = maxprim
   Alloc%maxPrimTUVLHS = max(MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK)
   Alloc%maxnAngLHS     = maxnAngmom
   Alloc%maxContLHS    = maxCont
   Alloc%maxAngmomLHS  = max(maxA,maxB,Alloc%maxAngmomLHS)
   Alloc%maxTUVLHS        = max(maxTUV,Alloc%maxTUVLHS)
   Alloc%maxPrimTUVijkLHS = max(MAXPRIMTUVIJK,Alloc%maxPrimTUVijkLHS) 
   Alloc%maxTotOrbLHS = max(maxTotorb,Alloc%maxTotOrbLHS)
   Alloc%maxijkLHS = max(maxijk2,Alloc%maxijkLHS)
   Alloc%maxETUVlenLHS = maxlenETUV
ELSE
   Alloc%maxJRHS = maxAngmom 
   Alloc%maxPrimRHS    = maxprim
   Alloc%maxPrimTUVRHS = max(MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK)
   Alloc%maxnAngRHS     = maxnAngmom
   Alloc%maxContRHS    = maxCont
   Alloc%maxAngmomRHS  = max(maxA,maxB,Alloc%maxAngmomRHS)
   Alloc%maxTUVRHS        = max(maxTUV,Alloc%maxTUVRHS)
   Alloc%maxPrimTUVijkRHS = max(MAXPRIMTUVIJK,Alloc%maxPrimTUVijkRHS) 
   Alloc%maxTotOrbRHS = max(maxTotorb,Alloc%maxTotOrbRHS)
   Alloc%maxijkRHS = max(maxijk2,Alloc%maxijkRHS)
   Alloc%maxETUVlenRHS = maxlenETUV
ENDIF

END SUBROUTINE SET_ALLOC

!> \brief wrapper routine which branch out and build the primitive hermite integral 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param integral contains arrays to store intermidiates and final integrals
!> \param input the integral input specifications, contains all info about what to do
!> \param PQ contains integrand info like reduced exponents, integral prefactor, etc
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
implicit none     
TYPE(Integrand)         :: PQ
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
  !Cartesian multipole order. Center at input%origo
  !Default is 0,0,0 center of coordinat system.
  CALL GETMOMTUV(INTEGRAL,INPUT,PQ,PQ%nPrimitives,LUPRI,IPRINT,INPUT%DERIVORDER,.TRUE.)
CASE ('Mulmom')
  !multipolemoments. Same as Carmom but center of the multipole 
  !moment integral is the center of the contracted ODbatch.
  CALL GETMOMTUV(INTEGRAL,INPUT,PQ,PQ%nPrimitives,LUPRI,IPRINT,INPUT%DERIVORDER,.FALSE.)
CASE DEFAULT
  WRITE(LUPRI,'(1X,A,A16)') 'Programming error! Not a case in Build_HermiteTUV:',PQ%Operator 
END SELECT
END SUBROUTINE Build_HermiteTUV

!> \brief builds the primitive hermite integral  
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param integral contains arrays to store intermidiates and final integrals
!> \param input the integral input specifications, contains all info about what to do
!> \param PQ contains integrand info like reduced exponents, integral prefactor, etc
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param OPERATORLABEL the label describing the operator 
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
!real(realk),pointer :: SJ000(:,:)
real(realk),pointer :: SJ0002(:,:)
INTEGER                 :: LUPRI,SUM,J,K,T,U,V,TUV,IOFF
INTEGER                 :: nPrim,IPRINT,ntuv,L,I
INTEGER                 :: zeroX,zeroY,zeroZ,Jmax,Jstart
real(realk)             :: X0,Y0,Z0
!
Logical :: orderAP

orderAP = Input%orderAngPrim

NPrim=PQ%nPrimitives
JMAX=PQ%endAngmom
Jstart=PQ%startAngmom

SELECT CASE(OPERATORLABEL)
CASE ('Overlap')
   CALL buildSJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT)
CASE ('Coulomb')
   CALL buildRJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,INPUT%HIGH_RJ000_ACCURACY)
CASE ('Erfc   ')
   !FIXME: mem_alloc cannot be called with 0:JMAX!!!!!
   !call mem_alloc(SJ0002,0:JMAX,nPrim)
   ALLOCATE(SJ0002(0:JMAX,nPrim))
   CALL buildRJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,INPUT%HIGH_RJ000_ACCURACY)
   CALL buildErfRJ000(SJ0002,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,INPUT%ATTomega,INPUT%HIGH_RJ000_ACCURACY)
   CALL DAXPY(NPrim*(JMAX+1),INPUT%ATTbeta,SJ0002,1,Integral%IN,1)
   !call mem_dealloc(SJ0002)
   DEALLOCATE(SJ0002)
CASE ('Erf    ')
   CALL buildErfRJ000(Integral%IN,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,INPUT%ATTomega,INPUT%HIGH_RJ000_ACCURACY)
CASE ('Nucrep')
   CALL buildNuclearRJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,INPUT%HIGH_RJ000_ACCURACY)
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
     &              jStart,jMax,nPrim,ntuv,orderAP,lupri,iprint)

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
               TUV=Integral%TUV%tuvIndex(T,U,V)-Jstart*(Jstart+1)*(Jstart+2)/6
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

!> \brief build the Classical Coulomb-matrix contribution through FMM
!> \author S. Reine
!> \date 2010
!> \param INTOUT the integral output specifications, determines how the output should be given
!> \param n1 the size of the first dimension
!> \param n2 the size of the second dimension
!> \param input the integral input specifications, contains all info about what to do
SUBROUTINE JmatClassical(INTOUT,n1,n2,Input)
use ks_settings
implicit none
Integer     :: n1,n2
TYPE(INTEGRALOUTPUT) :: INTOUT
!Real(realk) :: Jmat(n1,n2)
TYPE(INTEGRALINPUT)  :: INPUT

! Dummy argument for mm_get_J_matrix - should be removed
Integer,parameter :: LWRK = 1
Real(realk)       :: WRK(LWRK)
Real(realk),pointer :: Jmat(:,:)
!
INTEGER     :: NLEVEL,LMAX,TLMAX,LEXTRA,ALGORITHM,T_CONTRACTOR
REAL(REALK) :: GRAIN,RPQMIN,SCREEN
LOGICAL     :: USEUMAT,BRFREE,GRSRCH,INCNN
LOGICAL     :: DYNLMAX,ALLSQR,TRSRCH,NOBOXP,CONTRACTED
LOGICAL :: full


!
    call mem_alloc(Jmat,n1,n2)
    Jmat = 0.d0
!   Set default values for the mm driver
    CALL mm_init_defs(NLEVEL,GRAIN,TLMAX,ALGORITHM,                            &
     &                USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                    &
     &                ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,SCREEN)

!Simen find appropriate variable to determine if the basis set is contracted/uncontracted
contracted = .TRUE.
    
!   Change default settings based on LSint-settings and input
   CALL mm_init_scheme(NLEVEL,GRAIN,Input%MM_LMAX,Input%MM_TLMAX,ALGORITHM,                     &
                       USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                  &
                       ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,        &
                       Input%MM_SCREENTHR,contracted)
! Calculate moments here?

    full = .TRUE.
    IF (do_increment) full = .FALSE.
    IF (Input%MM_NOONE) full = .FALSE.
    full = .FALSE. !Problem for the energy for closed shell systems

    IF (full) THEN
!     Both two-electron repulsion and nuclear-electron attraction contributions to the classical Coulomb matrix
!Simen Hack to cope with the (strange) order for density-fitting FMM
      CALL mm_get_J_matrix('FULL_J',Jmat,n2,n1,WRK,LWRK)
!     CALL mm_get_J_matrix('FULL_J',Jmat,n1,n2,WRK,LWRK)
    ELSE
!     Only two-electron contribution
!Simen Hack to cope with the (strange) order for density-fitting FMM
      CALL mm_get_J_matrix('TWO_EL',Jmat,n2,n1,WRK,LWRK)
!     CALL mm_get_J_matrix('TWO_EL',Jmat,n1,n2,WRK,LWRK)
    ENDIF
    call add_full_2dim_to_lstensor(intout%resultmat2,Jmat,n1,n2,1)
    call mem_dealloc(Jmat)

END SUBROUTINE JmatClassical

!> \brief build the Classical electron-nuclear attraction contribution through FMM
!> \author S. Reine
!> \date 2010
!> \param INTOUT the integral output specifications, determines how the output should be given
!> \param n1 the size of the first dimension
!> \param n2 the size of the second dimension
!> \param input the integral input specifications, contains all info about what to do
SUBROUTINE electronNuclearClassic(INTOUT,n1,n2,Input)
implicit none
Integer     :: n1,n2
TYPE(INTEGRALOUTPUT) :: INTOUT
TYPE(INTEGRALINPUT)  :: INPUT

! Dummy argument for mm_get_J_matrix - should be removed
Integer,parameter :: LWRK = 1
Real(realk)       :: WRK(LWRK)
Real(realk),pointer :: H1(:,:)
!
INTEGER     :: NLEVEL,LMAX,TLMAX,LEXTRA,ALGORITHM,T_CONTRACTOR
REAL(REALK) :: GRAIN,RPQMIN,SCREEN
LOGICAL     :: USEUMAT,BRFREE,GRSRCH,INCNN
LOGICAL     :: DYNLMAX,ALLSQR,TRSRCH,NOBOXP,CONTRACTED

!
    call mem_alloc(H1,n1,n2)
    H1 = 0.d0
!   Set default values for the mm driver
    CALL mm_init_defs(NLEVEL,GRAIN,TLMAX,ALGORITHM,                            &
     &                USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                    &
     &                ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,SCREEN)

!Simen find appropriate variable to determine if the basis set is contracted/uncontracted
contracted = .TRUE.

!   Change default settings based on LSint-settings and input
   CALL mm_init_scheme(NLEVEL,GRAIN,Input%MM_LMAX,Input%MM_TLMAX,ALGORITHM,                     &
                       USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                  &
                       ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,        &
                       Input%MM_SCREENTHR,contracted)
! Calculate moments here?

!   Calssical nuclear-electron attraction contribution
    CALL mm_get_J_matrix('ONE_EL',H1,n1,n2,WRK,LWRK)
    call add_full_2dim_to_lstensor(intout%resultmat2,H1,n1,n2,1)
    call mem_dealloc(H1)

END SUBROUTINE electronNuclearClassic

!> \brief determine if the integral should be calculated or screened away
!> \author S. Reine
!> \date 2010
!> \param P contains overlap distribution for the left hand side, electron 1
!> \param Q contains overlap distribution for the right hand side, electron 2
!> \param input the integral input specifications, contains all info about what to do
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
LOGICAL FUNCTION getScreening(P,Q,INPUT,LUPRI,IPRINT)
implicit none
TYPE(Overlap)        :: P
TYPE(Overlap)        :: Q
TYPE(INTEGRALINPUT)  :: INPUT
integer              :: LUPRI,IPRINT
!
logical     :: screen
Real(realk) :: distance,distance2,extentSum
Integer     :: dir

screen = .FALSE.
getScreening = .FALSE.
!Cauchy-Schwarz screening
IF (INPUT%CS_SCREEN) THEN
  screen = (P%maxGab*Q%maxGab.LE.INPUT%CS_THRESHOLD)
  IF (screen) THEN
    getScreening = .TRUE.
    RETURN
  ENDIF
ENDIF
!Classical or extent based screening
IF (INPUT%NonClassical_SCREEN.OR.INPUT%OE_SCREEN) THEN
  distance2 = 0.0d0
  DO dir=1,3
    distance = P%ODcenter(dir) - Q%ODcenter(dir)
    distance2 = distance2 + distance*distance
  ENDDO
  extentSum = P%ODextent + Q%ODextent
  screen = distance2.GT.(extentSum*extentSum)
  IF (screen) THEN
    getScreening = .TRUE.
    RETURN
  ENDIF
ENDIF

END FUNCTION getScreening

!> \brief determine if the integral should be calculated or screened away
!> \author S. Reine
!> \date 2010
!> \param P contains overlap distribution for the left hand side, electron 1
!> \param Q contains overlap distribution for the right hand side, electron 2
!> \param input the integral input specifications, contains all info about what to do
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
LOGICAL FUNCTION getScreeningOD(P,Q,INPUT,LUPRI,IPRINT)
implicit none
TYPE(ODBATCH)        :: P
TYPE(ODBATCH)        :: Q
TYPE(INTEGRALINPUT)  :: INPUT
integer              :: LUPRI,IPRINT
!
logical     :: screen
Real(realk) :: distance,distance2,extentSum
Integer     :: dir

screen = .FALSE.
getScreeningOD = .FALSE.
!Cauchy-Schwarz screening
IF (INPUT%CS_SCREEN) THEN
  screen = (P%maxGab*Q%maxGab.LE.INPUT%CS_THRESHOLD)
  IF (screen) THEN
    getScreeningOD = .TRUE.
    RETURN
  ENDIF
ENDIF
!Classical screening
IF (INPUT%NonClassical_SCREEN.OR.INPUT%OE_SCREEN) THEN
  distance2 = 0.0d0
  DO dir=1,3
    distance = P%ODcenter(dir) - Q%ODcenter(dir)
    distance2 = distance2 + distance*distance
  ENDDO
  extentSum = P%ODextent + Q%ODextent
  screen = distance2.GT.(extentSum*extentSum)
  IF (screen) THEN
    getScreeningOD = .TRUE.
    RETURN
  ENDIF
ENDIF

END FUNCTION getScreeningOD

!> \brief initialise the Linkshell item used in LinK
!> \author T. Kjaergaard
!> \date 2010
!> \param A the linkshell
!> \param I the size of the linkshell
subroutine LINKshell_init(A,I)
  use memory_handling
implicit none
type(LINKshell), intent(inout) :: A
integer, intent(in)            :: I
   !NULLIFY(A%Belms)   done in mem_alloc
   !NULLIFY(A%IODelms)
   A%DIM = I
   if (I == 0) then
      call mem_alloc(A%Belms,1)
      call mem_alloc(A%IODelms,1)
   else
      call mem_alloc(A%Belms,I)
      call mem_alloc(A%IODelms,I)
   endif
   mem_allocated_linkshell = mem_allocated_linkshell + &
                           & size(A%Belms)*mem_intsize + size(A%IODelms)*mem_intsize + &
                           & mem_intsize
   max_mem_used_linkshell = MAX(max_mem_used_linkshell, mem_allocated_linkshell)
   !ALLOCATE(A%Belms(I))
   !ALLOCATE(A%IODelms(I))
end subroutine LINKshell_init

!> \brief free the Linkshell item used in LinK
!> \author T. Kjaergaard
!> \date 2010
!> \param A the linkshell
subroutine LINKshell_free(A)
  use memory_handling
implicit none
type(LINKshell), intent(inout) :: A
   mem_allocated_linkshell = mem_allocated_linkshell - &
                           & size(A%Belms)*mem_intsize - size(A%IODelms)*mem_intsize - &
                           & mem_intsize
   call mem_dealloc(A%Belms)
   call mem_dealloc(A%IODelms)
end subroutine LINKshell_free

END MODULE integraldriver
