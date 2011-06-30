!> @file
!> Module containing soubroutines for settting up input for the Thermite driver
MODULE Integralinfo
use precision
use integral_type
use typedef
CONTAINS
!> \brief initiate the INTEGRALINPUT structure 
!> \author T. Kjaergaard and S. Reine
!> \date 2009 
!>
!> This routine initiate the INTEGRALINPUT structure and should always be 
!> call before the MAIN_INTEGRAL_DRIVER is called. In most cases the 
!> values are set equal to values contained in the setting%scheme
!> but other values are set to the default values and can then be changed
!> after this call
!>
SUBROUTINE init_integral_INPUT(INTINPUT,SETTING)
implicit none
!> the INTEGRALINPUT to initiate
TYPE(INTEGRALINPUT) :: INTINPUT
!> contains info about thresholds and integral evaluation schemes 
TYPE(LSSETTING)   :: SETTING

Nullify(INTINPUT%DMAT_LHS)
Nullify(INTINPUT%DMAT_RHS)
Nullify(INTINPUT%AO(1)%p)
Nullify(INTINPUT%AO(2)%p)
Nullify(INTINPUT%AO(3)%p)
Nullify(INTINPUT%AO(4)%p)
INTINPUT%exchangeFactor = SETTING%SCHEME%exchangeFactor
INTINPUT%HIGH_RJ000_ACCURACY = SETTING%SCHEME%HIGH_RJ000_ACCURACY
INTINPUT%currentFragment=1
INTINPUT%NDMAT_LHS=0
INTINPUT%NDMAT_RHS=0
INTINPUT%DRHS_SYM=.FALSE.
INTINPUT%DLHS_SYM=.FALSE.
INTINPUT%sameODs=.FALSE.
INTINPUT%sameRHSaos=.FALSE.
INTINPUT%sameLHSaos=.FALSE.
INTINPUT%RHS_DMAT=.FALSE.
INTINPUT%LHS_DMAT=.FALSE.
INTINPUT%CENTERS=0
INTINPUT%ndim=0
INTINPUT%DO_PASSES=SETTING%SCHEME%DOPASS
INTINPUT%sphericalEcoeff=.NOT.SETTING%SCHEME%nonSphericalETUV
INTINPUT%maxpasses=SETTING%SCHEME%maxpasses
INTINPUT%DO_FOCK=.FALSE.
INTINPUT%DO_COULOMB=.FALSE.
INTINPUT%DO_EXCHANGE=.FALSE.
INTINPUT%DO_ENERGY=.FALSE.
INTINPUT%DO_JENGINE=.FALSE.
INTINPUT%DO_LINK=.FALSE.
INTINPUT%SKIP_FMM=.FALSE.
INTINPUT%DO_FMM=.FALSE.
INTINPUT%PUREFMM=SETTING%SCHEME%PUREFMM
INTINPUT%DO_MULMOM=.FALSE.
INTINPUT%LU_MMDATA=-1
INTINPUT%LU_MMDATR=-1
INTINPUT%MMUnique_iD1 = 0
INTINPUT%MMUnique_iD2 = 0
INTINPUT%MMstartA = 0
INTINPUT%MMstartB = 0
INTINPUT%MMindex = 0
! Screen integrals based on non-classical extent
INTINPUT%NonClassical_SCREEN=.FALSE.
!!!!!!!!INTINPUT%FMM_lmaxLocal=0
!!!!!!!!INTINPUT%FMM_lmaxTranslationFarField=0
!!!!!!!!INTINPUT%FMM_screenThreshold=0.0d0
INTINPUT%MM_SCREENTHR = SETTING%SCHEME%MM_SCREEN*SETTING%SCHEME%THRESHOLD
INTINPUT%MM_NOONE = SETTING%SCHEME%MM_NO_ONE
INTINPUT%MM_TLMAX = SETTING%SCHEME%MM_TLMAX
INTINPUT%MM_LMAX  = SETTING%SCHEME%MM_LMAX
!Specifies calculation of CS-integrals
INTINPUT%CS_int=.FALSE.
!Specifies calculation of PS-integrals
INTINPUT%PS_int=.FALSE.
!Default is to set ETUV-tensor outside main integral loops. 
!Can be turned off by input-keyword .ETUVIL
INTINPUT%setETUVoutside=.NOT.SETTING%SCHEME%ETUVinsideLOOP
!Default is to set OVERLAP type outside main integral loops. 
!Can be turned off by input-keyword .OVERLAPIL
INTINPUT%orderAngPrim = SETTING%SCHEME%orderAngprim
INTINPUT%orderPQ      = .TRUE.
IF(.NOT.SETTING%SCHEME%orderAngprim)INTINPUT%orderPQ      = .FALSE.
!Specifies usage of CS-screening
INTINPUT%CS_SCREEN=.FALSE.
INTINPUT%CS_THRESHOLD=SETTING%SCHEME%CS_THRESHOLD*SETTING%SCHEME%THRESHOLD
! PRIMITIVE SCREENING
INTINPUT%PS_SCREEN=.FALSE.
INTINPUT%PS_THRESHOLD=SETTING%SCHEME%PS_THRESHOLD*SETTING%SCHEME%THRESHOLD
! Screen OD-batches based on AO-batch extent
INTINPUT%OD_SCREEN=.FALSE.
INTINPUT%OD_THRESHOLD=SETTING%SCHEME%OD_THRESHOLD*SETTING%SCHEME%THRESHOLD
!Specifies overlap integral screening based on OD-batch extents
INTINPUT%OE_SCREEN=.FALSE.
INTINPUT%OE_THRESHOLD=SETTING%SCHEME%OE_THRESHOLD*SETTING%SCHEME%THRESHOLD
!INTINPUT%ORIGO(1) = 0.d0
!INTINPUT%ORIGO(1) = 0.d0
!INTINPUT%ORIGO(1) = 0.d0
INTINPUT%NDERIVP = 1
INTINPUT%NDERIVQ = 1
INTINPUT%derOrderP = 0
INTINPUT%derOrderQ = 0
INTINPUT%DERIVORDER = 0
INTINPUT%AddToIntegral = .FALSE.
INTINPUT%DO_GRADIENT=.FALSE.
!Default factor for closed shells
INTINPUT%CoulombFactor = 2.0d0
INTINPUT%ATTomega = SETTING%SCHEME%CAMmu
INTINPUT%ATTalpha = SETTING%SCHEME%CAMalpha 
INTINPUT%ATTbeta = SETTING%SCHEME%CAMbeta/SETTING%SCHEME%CAMalpha 
INTINPUT%ATTFACTOR = .FALSE.
INTINPUT%FTUVmaxprim = SETTING%SCHEME%FTUVmaxprim
INTINPUT%decpacked = .FALSE.
INTINPUT%uselst_DRHS = .FALSE.
INTINPUT%uselst_DLHS = .FALSE.
NULLIFY(INTINPUT%LST_GAB_LHS)
NULLIFY(INTINPUT%LST_GAB_RHS)
NULLIFY(INTINPUT%LST_pGAB_LHS)
NULLIFY(INTINPUT%LST_pGAB_RHS)

END SUBROUTINE init_integral_INPUT

!!$SUBROUTINE setIntegralInputDensities(intinput,setting)
!!$implicit none
!!$TYPE(INTEGRALINPUT) :: INTINPUT
!!$TYPE(LSSETTING)   :: SETTING
!!$!
!!$INTEGER :: LHS1,LHS2,RHS1,RHS2,nbastLHS1,nbastLHS2,nbastRHS1,nbastRHS2
!!$INTEGER :: nDmatLHS,nDmatRHS,lupri,idmat
!!$LOGICAL :: useAO1,useAO2
!!$
!!$LHS1 = setting%LHSdmatAOindex1
!!$LHS2 = setting%LHSdmatAOindex2
!!$RHS1 = setting%RHSdmatAOindex1
!!$RHS2 = setting%RHSdmatAOindex2
!!$nbastLHS1 = INTINPUT%AOdim(LHS1)
!!$nbastLHS2 = INTINPUT%AOdim(LHS2)
!!$nbastRHS1 = INTINPUT%AOdim(RHS1)
!!$nbastRHS2 = INTINPUT%AOdim(RHS2)
!!$nDmatLHS = setting%nDmatLHS
!!$nDmatRHS = setting%nDmatRHS
!!$lupri = 6
!!$
!!$IF ((INTINPUT%operator.EQ.'Kinetic') .AND.setting%RHSdfull)&
!!$     &CALL LSQUIT('Error in setIntegralInputDensities. Kinetic and RHS density')
!!$
!!$IF (setting%LHSdfull) THEN
!!$   CALL attachDmatToInput(INTINPUT,Setting%DfullLHS,nbastLHS1,nbastLHS2,nDmatLHS,'LHS')
!!$   useAO1 = .NOT.INTINPUT%AO(LHS1)%p%empty
!!$   useAO2 = .NOT.INTINPUT%AO(LHS2)%p%empty
!!$   CALL Build_lstensor_from_full_3dim(INTINPUT%lst_DLHS,Setting%DfullLHS,INTINPUT%AO(LHS1)%p,&
!!$        &INTINPUT%AO(LHS2)%p,nbastLHS1,nbastLHS2,nDmatLHS,useAO1,useAO2,lupri)
!!$ELSEIF(setting%LHSdmat)THEN
!!$   IF((INTINPUT%DO_EXCHANGE .AND. INTINPUT%DO_DALINK).AND. matrix_type .EQ. mtype_unres_dense)THEN
!!$      !unrestricted is special for exchange
!!$      call mem_alloc(setting%DfullLHS,nbastLHS1,nbastLHS2,2*nDmatRHS)
!!$      DO idmat=1,2*nDmatRHS,2
!!$         CALL DCOPY(nbastLHS1*nbastLHS2,setting%DmatLHS(idmat)%p%elms, 1,setting%DfullLHS(:,:,idmat),  1)
!!$         CALL DCOPY(nbastLHS1*nbastLHS2,setting%DmatLHS(idmat)%p%elmsb, 1,setting%DfullLHS(:,:,idmat+1),  1)
!!$      ENDDO
!!$      CALL attachDmatToInput(INTINPUT,Setting%DfullLHS,nbastLHS1,nbastLHS2,2*nDmatLHS,'LHS')
!!$      CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dLHS,Setting%DfullLHS,INTINPUT%AO(LHS1)%p,&
!!$           &INTINPUT%AO(LHS2)%p,nbastLHS1,nbastLHS2,2*nDmatLHS,useAO1,useAO2,lupri)
!!$   ELSE
!!$      IF(setting%DmatLHS(1)%p%nrow .NE. nbastLHS1)CALL LSQUIT('dimension mismatch in setIntegralInputDensities 1')
!!$      IF(setting%DmatLHS(1)%p%ncol .NE. nbastLHS2)CALL LSQUIT('dimension mismatch in setIntegralInputDensities 2')
!!$      call mem_alloc(setting%DfullLHS,nbastLHS1,nbastLHS2,setting%nDmatLHS)
!!$      CALL ls_mat_retrive_block(setting%DmatLHS,setting%DfullLHS,ndmatLHS,nbastLHS1,nbastLHS2,1,1,.FALSE.)
!!$      CALL attachDmatToInput(INTINPUT,Setting%DfullLHS,nbastLHS1,nbastLHS2,nDmatLHS,'LHS')
!!$      CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dLHS,Setting%DfullLHS,INTINPUT%AO(LHS1)%p,&
!!$           &INTINPUT%AO(LHS2)%p,nbastLHS1,nbastLHS2,nDmatLHS,useAO1,useAO2,lupri)
!!$   ENDIF
!!$ENDIF
!!$   
!!$IF (setting%RHSdfull)THEN
!!$   CALL attachDmatToInput(INTINPUT,Setting%DfullRHS,nbastRHS1,nbastRHS2,setting%nDmatRHS,'RHS')
!!$   useAO1 = .NOT.INTINPUT%AO(RHS1)%p%empty
!!$   useAO2 = .NOT.INTINPUT%AO(RHS2)%p%empty
!!$   CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dRHS,Setting%DfullRHS,INTINPUT%AO(RHS1)%p,&
!!$        &INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,lupri)
!!$ELSEIF(setting%RHSdmat)THEN
!!$   IF(INTINPUT%DO_EXCHANGE.AND. (matrix_type .EQ. mtype_unres_dense))THEN
!!$      !unrestricted is special for exchange
!!$      call mem_alloc(setting%DfullLHS,nbastLHS1,nbastLHS2,2*nDmatRHS)
!!$      DO idmat=1,2*nDmatRHS,2
!!$         CALL DCOPY(nbastLHS1*nbastLHS2,setting%DmatLHS(idmat)%p%elms, 1,setting%DfullLHS(:,:,idmat),  1)
!!$         CALL DCOPY(nbastLHS1*nbastLHS2,setting%DmatLHS(idmat)%p%elmsb, 1,setting%DfullLHS(:,:,idmat+1),  1)
!!$      ENDDO
!!$      CALL attachDmatToInput(INTINPUT,Setting%DfullRHS,nbastRHS1,nbastRHS2,2*nDmatRHS,'RHS')
!!$      CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dRHS,Setting%DfullRHS,INTINPUT%AO(RHS1)%p,&
!!$           &INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,2*nDmatRHS,useAO1,useAO2,lupri)
!!$   ELSE
!!$      IF(setting%DmatRHS(1)%p%nrow .NE. nbastRHS1)CALL LSQUIT('dimension mismatch in setIntegralInputDensities 3')
!!$      IF(setting%DmatRHS(1)%p%ncol .NE. nbastRHS2)CALL LSQUIT('dimension mismatch in setIntegralInputDensities 4')
!!$      call mem_alloc(setting%DfullRHS,nbastRHS1,nbastRHS2,setting%nDmatRHS)
!!$      CALL ls_mat_retrive_block(setting%DmatRHS,setting%DfullRHS,setting%ndmatRHS,&
!!$           &                      nbastRHS1,nbastRHS2,1,1,.FALSE.)
!!$      CALL attachDmatToInput(INTINPUT,Setting%DfullRHS,nbastRHS1,nbastRHS2,nDmatRHS,'RHS')
!!$      CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dRHS,Setting%DfullRHS,INTINPUT%AO(RHS1)%p,&
!!$           &INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,lupri)
!!$   ENDIF
!!$ENDIF
!!$
!!$END SUBROUTINE setIntegralInputDensities

END MODULE INTEGRALINFO
