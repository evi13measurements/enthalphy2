!> @file
!> Contains soubroutines that bridges integral interface routines to the main Thermite driver
MODULE ls_Integral_Interface
  use TYPEDEF
  use Integralinfo
  use integraldriver
  use BUILDAOBATCH
  use lstiming
  use BUILDBASISSET
  use io
  use precision
  use Matrix_module
  use Matrix_Operations
  use Fragment
  use molecule_module

  INTERFACE ls_attachDmatToSetting
     MODULE PROCEDURE ls_attachDmatToSetting_matsingle,&
                    & ls_attachDmatToSetting_matsinglep,&
                    & ls_attachDmatToSetting_matarray,&
                    & ls_attachDmatToSetting_matarrayp,&
                    & ls_attachDmatToSetting_full
  END INTERFACE

CONTAINS
!> \brief Does Sherical Transformation 
!> \author T. Kjaergaard
!> \date 2010
!> \param SubBlockInt the integrals to be transformed
!> \param nbast1 number of basis functions for the 1. center
!> \param nbast2 number of basis functions for the 2. center
!> \param ndMAT1 number of matrices before the transformation
!> \param ndmat2 number of matrices after the transformation
!> \param nderiv level of derivative 
!> \param iprint level of output printing  
!> \param lupri Default print unit
SUBROUTINE ShericalMomentTransformation(SubBlockInt,nbast1,nbast2,ndMAT1,ndmat2,nderiv,&
     &                                  iprint,lupri)
implicit none
Type(LSSETTING)     :: SETTING
Integer             :: nbast1,nbast2,ndMAT1,ndmat2,nderiv,iprint,lupri
Real(realk),pointer :: SubBlockInt(:,:,:,:,:)
!
Real(realk),pointer :: Integrals2(:,:,:,:,:)
!
call mem_alloc(Integrals2,nbast1,nbast2,1,1,ndmat2)
CALL LS_DZERO(Integrals2,nbast1*nbast2*ndmat2)
CALL SPHERICAL_TRANSFORMATION(SubBlockInt,Integrals2,nbast1,nbast2,ndMAT1,ndmat2,nderiv,&
     &                        iprint,lupri)
call mem_dealloc(SubBlockInt)
SubBlockInt => Integrals2
END SUBROUTINE ShericalMomentTransformation

!> \brief Adds a sub block to a full tensor
!> \author S. Reine
!> \date 2010
!> \param integrals the full tensir
!> \param SubBlockInt the subblock to be added
!> \param sameAOsLHS
!> \param sameFragmentLHS
!> \param sameAOsRHS
!> \param sameFragmentRHS
!> \param nbast1 number of basis functions for the 1. center
!> \param nbast2 number of basis functions for the 2. center
!> \param nbast3 number of basis functions for the 3. center
!> \param nbast4 number of basis functions for the 4. center
!> \param start1 the start index for the 1. dimension 
!> \param start2 the start index for the 2. dimension 
!> \param start3 the start index for the 3. dimension 
!> \param start4 the start index for the 4. dimension 
!> \param end1 the end index for the 1. dimension 
!> \param end2 the end index for the 2. dimension 
!> \param end3 the end index for the 3. dimension 
!> \param end4 the end index for the 4. dimension 
!> \param ndmat the number of matrices (the fifth dimension)
SUBROUTINE addSubBlocks(integrals,SubBlockInt,sameAOsLHS,sameFragmentLHS,&
     &                  sameAOsRHS,sameFragmentRHS,nbast1,nbast2,nbast3,nbast4,&
     &                  start1,start2,start3,start4,end1,end2,end3,end4,ndmat)
implicit none
Real(realk), pointer :: Integrals(:,:,:,:,:)
Real(realk), pointer :: SubBlockInt(:,:,:,:,:)
Integer              :: ndmat
type(matrix)         :: tmp
Integer              :: nbast1,nbast2,nbast3,nbast4
Integer              :: start1,start2,start3,start4
Integer              :: end1,end2,end3,end4
Logical              :: sameAOsLHS,sameAOsRHS
Logical              :: sameFragmentLHS,sameFragmentRHS

! Using full integral matrices - insert sub-blocks into final integrals
  CALL addSubBlockFull(integrals,SubBlockInt,sameAOsLHS,sameFragmentLHS,&
     &                 sameAOsRHS,sameFragmentRHS,nbast1,nbast2,nbast3,nbast4,&
     &                 start1,start2,start3,start4,end1,end2,end3,end4,ndmat)

END SUBROUTINE addSubBlocks

!> \brief Generalized wrapper routine to get explicit integrals for given operator Oper (in case of MPI) 
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param ndmat the size of the fifth dimension - number of matrices
!> \param nderiv the order of derivative 
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> Generalized routine to get explicit integrals for given operator Oper 
!>
!>     'Coulomb'       1/r_12
!>     'Overlap'       dirac_delta(r_1-r_2)
!>     'Kinetic'       -1/2 nabla
!>
!> and with given AOs according to AO1,AO2,AO3 and AO4:
!>
!>     'Regular'       Regular basis
!>     'Huckel'        Hückel basis
!>     'DF-Aux'        Auxiliary basis for density-fitting
!>     'Empty'         Empty, used for two and three-center integrals
!>
!> AO1 and AO2 belong to electron 1, and AO3 and AO4 to electron 2
!>
!> First four dimensions of Integrals should correspond to the dimensions of 
!> AO1-4, the fifth dimension is an added index allowing for example for 
!> different derivative integrals.
!>
SUBROUTINE ls_getIntegrals(AO1,AO2,AO3,AO4,ndmat,nderiv,&
     &                     Oper,Spec,intType,SETTING,LUPRI,LUERR)
#ifdef VAR_LSMPI
  use infpar_module
  use lsmpi_mod
#endif
implicit none
Integer              :: LUPRI,LUERR,ndmat,nderiv
Real(realk), pointer :: Integrals(:,:,:,:,:)
Real(realk), pointer :: Integrals2(:,:,:,:,:)
Real(realk), pointer :: SubBlockInt(:,:,:,:,:)
type(matrix)         :: tmp
Character*(*)        :: AO1,AO2,AO3,AO4,Oper,intType,Spec
Character(len=12)    :: Oper1
integer              :: lenstring
Type(LSSETTING)      :: SETTING
Integer              :: I,J,DIMVEC(2)
Integer              :: I1,I2,I3,I4,iNode
Integer              :: nbast1,nbast2,nbast3,nbast4,n1,n2,s1,s2,s3,s4
Integer              :: start1,start2,start3,start4,ndim
Integer              :: end1,end2,end3,end4,imat,ndmat1,ndmat2
real(realk),parameter :: ONE = 1.D0
logical              :: sphmom,unres
Logical              :: sameAOsLHS,sameAOsRHS,sameODs
Logical              :: sameFragmentLHS,sameFragmentRHS
Real(realk),pointer  :: SubBlock(:,:),SubBlockTrans(:,:)
Real(realk),pointer  :: DblockLHS(:,:,:),DblockRHS(:,:,:)
Logical             :: FRAGMENT_SAVE,permuteLHS,permuteRHS
Real(realk)          :: tstart,tend
integer :: ndim2(5),ierr,itask
ndim2 = setting%output%ndim
CALL LSTIMER('START ',tstart,tend,lupri)

CALL ls_setDefaultFragments(setting)

ndmat1 = ndmat
ndmat2 = ndmat
unres = matrix_type .EQ. mtype_unres_dense !global variable
sphmom = Oper .EQ. 'Sphmom'
lenstring = len(Oper)
IF (sphmom) THEN
   ndmat1 = (nderiv+1)*(nderiv+2)*(nderiv+3)/6 !# of cartesian matrices
   Oper1(1:6) = 'Carmom'
ELSE
   Oper1(1:lenstring) = Oper
ENDIF

#ifdef VAR_LSMPI
! ***************************************************************************
! *                                MPI Specific                             *
! ***************************************************************************
  IF (infpar%mynum.EQ.infpar%master) THEN
!   Brano: Spawn here!
    call ls_mpibcast('LSGETINT',8,infpar%master)
    call lsmpi_getIntegrals_masterToSlave(AO1,AO2,AO3,AO4,ndmat,nderiv,&
     &                                     Oper,Spec,intType,SETTING,LUPRI,LUERR)
  ENDIF
  allocate(integrals(ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5)))
  call LS_dzero(integrals,ndim2(1)*ndim2(2)*ndim2(3)*ndim2(4)*ndim2(5))
  SETTING%SCHEME%PUREFMM=.TRUE.
  CALL buildFragmentInfoAndBlocks(SETTING,AO1,AO2,AO3,AO4,LUPRI,LUERR)
  CALL GetSymmetries(sameAOsLHS,sameAOsRHS,sameODs,AO1,AO2,AO3,AO4)
  itask = 1
  DO I=1,SETTING%FRAGMENTS%LHSblock%numBlocks
  permuteLHS = SETTING%FRAGMENTS%LHSblock%sameAOs.AND.&
     &     (.NOT.SETTING%FRAGMENTS%LHSblock%Blocks(I)%sameFragments)
  DO J=1,SETTING%FRAGMENTS%RHSblock%numBlocks
  itask = itask + 1
  write(*,*) 'debug:IJ',infpar%mynum,I,J,itask,mod(itask,infpar%nodtot),infpar%nodtot
  IF (mod(itask,infpar%nodtot).NE.infpar%mynum) CYCLE
  permuteRHS = SETTING%FRAGMENTS%RHSblock%sameAOs.AND.&
     &     (.NOT.SETTING%FRAGMENTS%RHSblock%Blocks(J)%sameFragments)
  CALL SetDaltonFragments(SETTING,I,J,sameFragmentLHS,sameFragmentRHS,lupri)
  CALL GetDaltonOrbitalInfo(SETTING,I,J,nbast1,nbast2,nbast3,nbast4,&
     &           start1,start2,start3,start4,end1,end2,end3,end4)
                 
  ! Copy density-matrix sub-block
  CALL ls_setDblock(setting,nbast1,nbast2,nbast3,nbast4,&
     &            start1,start2,start3,start4,permuteLHS,permuteRHS)
  call initIntegralOutputDims(setting%Output,nbast1,nbast2,nbast3,nbast4,ndmat1)
#else
! ***************************************************************************
! *                                  Serial                                 *
! ***************************************************************************
  IF(sphmom)call initIntegralOutputDims(setting%output,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndmat1)
!  DO iAO=1,setting%nAO
!    setting%fragment(iAO)%p => setting%molecule(iAO)%p
!  ENDDO
#endif

! ***************************************************************************
! *                     Actual calculation of integrals                     *
! ***************************************************************************
  CALL ls_getIntegrals1(AO1,AO2,AO3,AO4,ndmat1,nderiv,Oper1(1:lenstring),&
     &                  Spec,intType,SETTING,LUPRI,LUERR)


#ifndef VAR_LSMPI
! ***************************************************************************
! *                                  Serial                                 *
! ***************************************************************************
  Call ls_freeDmatFromSetting(setting)
  !Free full size density matrices
  !   CALL ls_freeDfull(setting)
  IF(sphmom)THEN
     call mem_alloc(SubBlockInt,ndim2(1),ndim2(2),1,1,ndmat1)
     CALL retrieve_Output(lupri,setting,SubBlockInt)
     call initIntegralOutputDims(setting%output,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndmat2)
     call mem_alloc(Integrals2,ndim2(1),ndim2(2),1,1,ndmat2)
     CALL LS_DZERO(Integrals2,ndim2(1)*ndim2(2)*ndmat2)
     CALL SPHERICAL_TRANSFORMATION(SubBlockInt,Integrals2,&
          &ndim2(1),ndim2(2),ndMAT1,ndmat2,nderiv,&
          &SETTING%SCHEME%INTPRINT,lupri)
     call mem_dealloc(SubBlockInt)
     call mem_alloc(SubBlockInt,ndim2(1),ndim2(1),1,1,ndmat2)
     CALL DCOPY(ndim2(1)*ndim2(2)*ndmat2,Integrals2,1,SubBlockInt,1)
     call mem_dealloc(Integrals2)
     call ls_Build_lstensor_from_full_5dim(setting%output%resultmat2,setting,AO1,AO2,AO3,AO4,ndmat2,intType,&
          &LUPRI,LUERR,SubBlockInt)
     call mem_dealloc(SubBlockInt)
  ENDIF
#else
! ***************************************************************************
! *                                MPI Specific                             *
! ***************************************************************************
  call mem_alloc(SubBlockInt,nbast1,nbast2,nbast3,nbast4,ndmat1)
  Call LS_DZERO(SubBlockInt,nbast1*nbast2*nbast3*nbast4*ndmat1)
  CALL retrieve_Output(lupri,setting,SubBlockInt)
  CALL FreeDaltonFragments(SETTING)
  CALL ls_freeDblock(setting)
  setting%RHSdfull = .FALSE.
  setting%LHSdfull = .FALSE.
               
! Special for spherical multipole-moment integrals               
  IF (sphmom) CALL ShericalMomentTransformation(SubBlockInt,nbast1,nbast2,ndMAT1,ndmat2,&
     &                           nderiv,SETTING%SCHEME%INTPRINT,lupri)
               
! Place sub-block integrals into the correct matrix, utilizing permutational symmetry
  CALL addSubBlocks(integrals,SubBlockInt,sameAOsLHS,sameFragmentLHS,&
     &              sameAOsRHS,sameFragmentRHS,nbast1,nbast2,nbast3,nbast4,&
     &              start1,start2,start3,start4,end1,end2,end3,end4,ndmat2)               
  call mem_dealloc(SubBlockInt)
  ENDDO !J
  ENDDO !I
write(6,*) 'debug:node',infpar%mynum
call lsheader(6,'Overlap matrix')
call output(integrals,1,ndim2(1),1,ndim2(1),ndim2(1),ndim2(1),1,6)
!ToDo Retrieve output from all the nodes and collect them to the master

setting%output%ndim = ndim2
call initIntegralOutputDims(setting%Output,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5))
call ls_Build_lstensor_from_full_5dim(setting%output%resultmat2,setting,AO1,AO2,AO3,AO4,ndmat,intType,LUPRI,LUERR,integrals)

CALL freeFragmentInfoAndBlocks(SETTING)
IF(Oper .EQ. 'Nucrep' .AND. SETTING%SCHEME%FMM)THEN
  FRAGMENT_SAVE = SETTING%SCHEME%FRAGMENT
  SETTING%SCHEME%FRAGMENT = .FALSE.
  Ndim=SIZE(Integrals)/2
  CALL ls_electronNuclearClassic(AO1,AO2,AO3,AO4,ndmat,nderiv,Oper1(1:lenstring),intType,SETTING,LUPRI,LUERR)
  SETTING%SCHEME%FRAGMENT = FRAGMENT_SAVE
ENDIF
SETTING%SCHEME%PUREFMM=.FALSE.
deallocate(integrals)
#endif

!If density-matrices have been assigned to settings free them
CALL ls_freeDmatFromSetting(setting)
!CALL LSTIMER('ls_getIntegrals',tstart,tend,lupri)

END SUBROUTINE ls_getIntegrals

!> \brief Generalized routine to get explicit integrals for given operator Oper 
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param ndmat the size of the fifth dimension - number of matrices
!> \param nderiv the order of derivative 
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!>
!> Generalized routine to get explicit integrals for given operator Oper 
!>
!>     'Coulomb'       1/r_12
!>     'Overlap'       dirac_delta(r_1-r_2)
!>     'Kinetic'       -1/2 nabla
!>
!> and with given AOs according to AO1,AO2,AO3 and AO4:
!>
!>     'Regular'       Regular basis
!>     'Huckel'        Hückel basis
!>     'DF-Aux'        Auxiliary basis for density-fitting
!>     'Empty'         Empty, used for two and three-center integrals
!>
!> AO1 and AO2 belong to electron 1, and AO3 and AO4 to electron 2
!>
!> First four dimensions of Integrals should correspond to the dimensions of 
!> AO1-4, the fifth dimension is an added index allowing for example for 
!> different derivative integrals.
!>
SUBROUTINE ls_getIntegrals1(AO1,AO2,AO3,AO4,ndmat,nderiv,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
Character*(*)        :: AO1,AO2,AO3,AO4,Oper,Spec,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR,ndmat,nderiv
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,nMAT,nAtoms,iMol
Logical              :: Classical,DO4CENTERCOULOMB, buildGab,doscreen,primScreen,batchOnly
Real(realk),pointer  :: derivativeIntegrals(:,:,:,:,:)
integer :: ndim2(5)
ndim2 = setting%output%ndim 
buildGab = .FALSE.
CALL init_integral_input(INT_INPUT,SETTING)

INT_INPUT%operator = Oper
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%NonClassical_SCREEN = .FALSE.

! Specifies if we calculate four-center integrals
DO4CENTERCOULOMB = (Oper .EQ. 'Coulomb').AND.(AO1 .EQ. 'Regular').AND.(AO2 .NE. 'Regular') &
     & .AND.(AO3 .NE. 'Regular').AND.(AO4 .NE. 'Regular')

! Specifies if integrals are limited to one ore more AO-batches
batchOnly = (setting%batchindex(1).NE.0).OR.(setting%batchindex(2).NE.0) &
     &      .OR.(setting%batchindex(3).NE.0).OR.(setting%batchindex(4).NE.0)

! Screen only when calculationg three- and four-center integrals
! Currently turned off when calculating integrals only for a specific batch
doscreen = (((AO1.NE.'Empty').AND.(AO2.NE.'Empty').AND.((AO3.NE.'Empty').OR.(AO4.NE.'Empty'))) &
     & .OR.((AO1.NE.'Empty').OR.(AO2.NE.'Empty')).AND.((AO3.NE.'Empty').AND.(AO4.NE.'Empty'))) & 
     & .AND. (.NOT. batchOnly)

primScreen = (.NOT.SETTING%SCHEME%DECPACKED).AND.(.NOT.DO4CENTERCOULOMB)
IF (doscreen) THEN
  Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
  ! Primitive screening turned off for DECPACKED
  IF (primScreen) Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
  CALL get_screening_matrices(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
  setting%output%ndim = ndim2
ENDIF
IF(Oper .EQ. 'Nucrep') THEN

  INT_INPUT%AddToIntegral = .TRUE.
  
  Classical = SETTING%SCHEME%FMM .AND. (.NOT. SETTING%SCHEME%MM_NO_ONE)
  IF (Classical) THEN
    Int_input%NonClassical_SCREEN = SETTING%SCHEME%FMM
    Int_input%DO_FMM              = SETTING%SCHEME%FMM
    INT_INPUT%MM_TLMAX            = SETTING%SCHEME%MM_TLMAX
    INT_INPUT%MM_LMAX             = SETTING%SCHEME%MM_LMAX
    INT_INPUT%MM_SCREENTHR        = SETTING%SCHEME%MM_SCREEN*SETTING%SCHEME%THRESHOLD
  ENDIF
  IF(SETTING%SCHEME%PUREFMM) Int_input%DO_FMM =.FALSE.
ENDIF

IF (DO4CENTERCOULOMB) setting%output%ndim = ndim2

CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

IF(Oper .EQ. 'Carmom') THEN
   INT_INPUT%ORIGO(1) = 0.d0
   INT_INPUT%ORIGO(2) = 0.d0
   INT_INPUT%ORIGO(3) = 0.d0
   INT_INPUT%NDERIVQ = ndmat
   INT_INPUT%DERIVORDER = nderiv
ENDIF

SELECT CASE(Spec)
CASE('Regular')
   !Default undifferentiated case      
   setting%Output%dograd = .false.
   IF(SETTING%SCHEME%DECPACKED)THEN
      INT_INPUT%DECPACKED = .TRUE.
      CALL initIntegralOutputDims(setting%Output,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5))
      call init_lstensor_5dim(setting%output%resultmat2,Int_Input%AO(1)%p,Int_Input%AO(4)%p,&
      & Int_Input%AO(2)%p,Int_Input%AO(3)%p,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5),&
      &.TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,lupri)
!      &.TRUE.,.TRUE.,.TRUE.,.TRUE.,INT_INPUT%OD_SCREEN,lupri)
      call add_contrule(setting%output,1,4,2,3,5,0,0,0,0,0,0,0,0,0,0,1.d0,.false.)
   ELSE
      CALL initIntegralOutputDims(setting%Output,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5))
      call init_lstensor_5dim(setting%output%resultmat2,Int_Input%AO(1)%p,Int_Input%AO(2)%p,Int_Input%AO(3)%p,&
      & Int_Input%AO(4)%p,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5),&
      & .TRUE.,.TRUE.,.TRUE.,.TRUE.,INT_INPUT%OD_SCREEN,lupri)
!      & .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,lupri)
      call add_contrule(setting%output,1,2,3,4,5,0,0,0,0,0,0,0,0,0,0,1.d0,.false.)
   ENDIF
CASE('Gradient')
   setting%Output%dograd = .true.
   nAtoms =  setting%molecule(1)%p%nAtoms
   CALL initIntegralOutputDims(setting%Output,3,nAtoms,1,1,ndmat)
   setting%Output%resultmat2%gradienttensor = .true.   
   call init_gradientlstensor(setting%output%resultmat2,nAtoms,setting%Output%ndim(5),lupri)
!   call init_gradientlstensor(setting%output%resultmat2,Int_Input%AO(1)%p,ndmat,lupri)
   DO iMol=2,4
      IF (setting%molecule(iMol)%p%nAtoms.NE.nAtoms) CALL LSQUIT('Error in ls_getIntegrals1. nAtoms inconsistency!',-1)
   ENDDO

   !Molecular gradient
   INT_INPUT%DO_GRADIENT = .TRUE.
   INT_INPUT%sphericalEcoeff = .FALSE.
   INT_INPUT%derOrderP  = 1
   INT_INPUT%derOrderQ  = 0
   INT_INPUT%sameODs  = .FALSE.
   IF ((AO1.EQ.'Empty').AND.(AO2.EQ.'Empty')) THEN
      CALL LSQUIT('Error in ls_getIntegrals1. DO_GRADIENT and AO1=AO2="Empty"',-1)
   ELSEIF ((AO1.EQ.'Empty').OR.(AO2.EQ.'Empty')) THEN
      INT_INPUT%NDERIVP = 3
   ELSE
      INT_INPUT%NDERIVP = 6
   ENDIF
   INT_INPUT%NDERIVQ = 1
   INT_INPUT%DERIVORDER = max(INT_INPUT%derOrderP,INT_INPUT%derOrderQ)
   IF (Oper.EQ.'Kinetic') INT_INPUT%sameODs = .FALSE.
END SELECT

CALL ls_setIntegralInputDensities(INT_INPUT,SETTING) !takes the D in setting and plug into input

call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
call ls_freeIntegralInputDensities(INT_INPUT)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
IF(doscreen) CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)
IF(Spec(1:7).EQ.'Regular') call free_contrule(setting%output)
setting%output%ndim = ndim2

END SUBROUTINE ls_getIntegrals1

!> \brief Adds a sub block to a full tensor
!> \author S. Reine
!> \date 2010
!> \param integrals the full tensir
!> \param SubBlock the subblock to be added
!> \param sameAOsLHS if both LHS AOs are the same
!> \param sameFragmentLHS if both LHS fragments are the same
!> \param sameAOsRHS if both RHS AOs are the same
!> \param sameFragmentRHS if both RHS fragments are the same
!> \param nbast1 number of basis functions for the 1. center
!> \param nbast2 number of basis functions for the 2. center
!> \param nbast3 number of basis functions for the 3. center
!> \param nbast4 number of basis functions for the 4. center
!> \param start1 the start index for the 1. dimension 
!> \param start2 the start index for the 2. dimension 
!> \param start3 the start index for the 3. dimension 
!> \param start4 the start index for the 4. dimension 
!> \param end1 the end index for the 1. dimension 
!> \param end2 the end index for the 2. dimension 
!> \param end3 the end index for the 3. dimension 
!> \param end4 the end index for the 4. dimension 
!> \param ndmat the number of matrices (the fifth dimension)
SUBROUTINE addSubBlockFull(integrals,SubBlock,sameAOsLHS,sameFragmentLHS,&
     &                  sameAOsRHS,sameFragmentRHS,nbast1,nbast2,nbast3,nbast4,&
     &                  start1,start2,start3,start4,end1,end2,end3,end4,nMat)
implicit none
integer             :: nbast1,nbast2,nbast3,nbast4,nMat
integer             :: start1,start2,start3,start4
integer             :: end1,end2,end3,end4
Real(realk),pointer :: integrals(:,:,:,:,:)
Real(realk),pointer :: SubBlock(:,:,:,:,:)
Real(realk),pointer :: SubBlockPerm(:,:,:,:,:)
Logical             :: sameAOsLHS,sameAOsRHS,sameFragmentLHS,sameFragmentRHS
!
Integer :: I1,I2,I3,I4,s2,s4,orb1,orb2,orb3,orb4,iMat
Logical :: permuteLHS,permuteRHS

permuteLHS = sameAOsLHS.AND..NOT.sameFragmentLHS
permuteRHS = sameAOsRHS.AND..NOT.sameFragmentRHS

! Regular block
integrals(start1:end1,start2:end2,start3:end3,start4:end4,1:nMat) = &
     &  integrals(start1:end1,start2:end2,start3:end3,start4:end4,1:nMat) &
     &+ SubBlock(:,:,:,:,1:nMat)

! Integral blocks using permutational symmetry
IF (permuteLHS) THEN
  call mem_alloc(SubBlockPerm,nbast2,nbast1,nbast3,nbast4,nMat)
  DO iMat=1,nMat
    DO I3=1,nbast3
      orb3 = start3 + I3 -1
      s4 = 1
      IF (sameFragmentRHS) s4=I3
      DO I4=s4,nbast4
        orb4 = start4 + I4 - 1                     
        DO I1=1,nbast1
          orb1 = start1 + I1 - 1
          s2 = 1
          IF (sameFragmentLHS) s2=I1
          DO I2=s2,nbast2
            orb2 = start2 + I2 - 1
            integrals(orb2,orb1,orb3,orb4,iMat) = integrals(orb2,orb1,orb3,orb4,iMat) + SubBlock(I1,I2,I3,I4,iMat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  call mem_dealloc(SubBlockPerm)
ENDIF
IF (permuteRHS) THEN
  call mem_alloc(SubBlockPerm,nbast2,nbast1,nbast3,nbast4,nMat)
  DO iMat=1,nMat
    DO I3=1,nbast3
      orb3 = start3 + I3 - 1
      s4 = 1
      IF (sameFragmentRHS) s4=I3
      DO I4=s4,nbast4
        orb4 = start4 + I4 - 1
        DO I1=1,nbast1
          orb1 = start1 + I1 - 1
          s2 = 1
          IF (sameFragmentLHS) s2=I1
          DO I2=s2,nbast2
            orb2 = start2 + I2 - 1
            integrals(orb1,orb2,orb4,orb3,iMat) = integrals(orb1,orb2,orb4,orb3,iMat) + SubBlock(I1,I2,I3,I4,iMat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  call mem_dealloc(SubBlockPerm)
ENDIF
IF (permuteLHS.AND.permuteRHS) THEN
  call mem_alloc(SubBlockPerm,nbast2,nbast1,nbast3,nbast4,nMat)
  DO iMat=1,nMat
    DO I3=1,nbast3
      orb3 = start3 + I3 - 1
      s4 = 1
      IF (sameFragmentRHS) s4=I3
      DO I4=s4,nbast4
        orb4 = start4 + I4 - 1 
        DO I1=1,nbast1
          orb1 = start1 + I1 - 1 
          s2 = 1
          IF (sameFragmentLHS) s2=I1
          DO I2=s2,nbast2
            orb2 = start2 + I2 - 1
            integrals(orb2,orb1,orb4,orb3,iMat) = integrals(orb2,orb1,orb4,orb3,iMat) + SubBlock(I1,I2,I3,I4,iMat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  call mem_dealloc(SubBlockPerm)
ENDIF

END SUBROUTINE addSubBlockFull

!> \brief Calculate the classical part (in an FMM sense) of the electron nuclear attraction contribution to the fock matrix
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param ndmat the size of the fifth dimension - number of matrices
!> \param nderiv the order of derivative 
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_electronNuclearClassic(AO1,AO2,AO3,AO4,ndmat,nderiv,Oper,intType,SETTING,LUPRI,LUERR)
implicit none
Character*(*)        :: AO1,AO2,AO3,AO4,Oper,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR,ndmat,nderiv
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,nMAT
Logical              :: Classical

CALL init_integral_input(INT_INPUT,SETTING)

INT_INPUT%operator = Oper
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%NonClassical_SCREEN = .FALSE.
  
Classical = SETTING%SCHEME%FMM .AND. (.NOT. SETTING%SCHEME%MM_NO_ONE)
IF (Classical) THEN
   Int_input%NonClassical_SCREEN = SETTING%SCHEME%FMM
   Int_input%DO_FMM              = SETTING%SCHEME%FMM
   INT_INPUT%MM_TLMAX            = SETTING%SCHEME%MM_TLMAX
   INT_INPUT%MM_LMAX             = SETTING%SCHEME%MM_LMAX
   INT_INPUT%MM_SCREENTHR        = SETTING%SCHEME%MM_SCREEN*SETTING%SCHEME%THRESHOLD
ENDIF

CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

CALL ls_setIntegralInputDensities(INT_INPUT,SETTING)

CALL initIntegralOutputDims(setting%Output,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,ndmat)
IF (INT_INPUT%DO_FMM) THEN
   call electronNuclearClassic(setting%Output,setting%Output%ndim(1),setting%Output%ndim(2),INT_INPUT)
ENDIF
call ls_freeIntegralInputDensities(INT_INPUT)
Call ls_freeDmatFromSetting(setting)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)

END SUBROUTINE ls_electronNuclearClassic

!> \brief Calculate the exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param Dsym true if the density matrix is symmetric
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_get_exchange_mat(Dsym,AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
Character*(*)        :: AO1,AO2,AO3,AO4,Oper,Spec,intType
logical              :: Dsym !symmetric Dmat
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,I,J,IDMAT,nAtoms,iMol
REAL(REALK)          :: AVERAG
logical              :: fragment
Real(realk),pointer :: TMP(:,:)
Real(realk),pointer :: TMP2(:,:,:,:,:)
Logical             :: FRAGMENT_SAVE
integer :: ndim2(5)
CALL ls_setDefaultFragments(setting)
ndim2 = setting%output%ndim 
do i=1,4
   setting%fragment(i)%p => setting%MOLECULE(i)%p
enddo
SETTING%sameFrag = SETTING%sameMOL

FRAGMENT_SAVE = SETTING%SCHEME%FRAGMENT
SETTING%SCHEME%FRAGMENT = .FALSE.
CALL init_integral_input(INT_INPUT,SETTING)
INT_INPUT%DO_FOCK = .TRUE.
INT_INPUT%exchangeFactor = SETTING%SCHEME%exchangeFactor
INT_INPUT%DO_EXCHANGE = (INT_INPUT%exchangeFactor.GT.0.0d0)&
     & .OR. (SETTING%SCHEME%CAMbeta .GT. 0.0d0)

INT_INPUT%DO_LINK = SETTING%SCHEME%LINK
INT_INPUT%DO_DALINK = SETTING%SCHEME%DALINK

IF(INT_INPUT%DO_LINK.AND. Dsym)THEN
   INT_INPUT%DRHS_SYM=.TRUE.
   INT_INPUT%DLHS_SYM=.TRUE.
ENDIF

IF (INT_INPUT%DO_EXCHANGE) THEN
  Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN

  Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
  Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
  CALL get_screening_matrices(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)           
  setting%output%ndim = ndim2
  IF(Oper .EQ. 'Nucrep') INT_INPUT%AddToIntegral = .TRUE.
  INT_INPUT%operator = Oper

  CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)
  CALL ls_setintegralinputdensities(INT_INPUT,SETTING)

  SELECT CASE(Spec)
  CASE ('Regular')
    CALL initIntegralOutputDims(setting%Output,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5))
    call init_lstensor_5dim(setting%output%resultmat2,Int_Input%AO(2)%p,Int_Input%AO(4)%p,&
         & Int_Input%AO(1)%p,Int_Input%AO(1)%p,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5),&
         & .TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,lupri)
    call add_contrule(setting%output,2,4,0,0,5,0,0,0,0,0,1,3,0,0,0,1.d0,.true.)
  CASE ('Gradient')
    nAtoms =  setting%molecule(1)%p%nAtoms
    DO iMol=2,4
      IF (setting%molecule(iMol)%p%nAtoms.NE.nAtoms) CALL LSQUIT('Error in ls_get_exchange_mat. nAtoms inconsistency!',lupri)
    ENDDO
    call init_gradientlstensor(setting%output%resultmat2,nAtoms,setting%Output%ndim(5),lupri)
    !Molecular gradient
    INT_INPUT%DO_GRADIENT = .TRUE.
    INT_INPUT%sphericalEcoeff = .FALSE.
    INT_INPUT%derOrderP  = 1
    INT_INPUT%derOrderQ  = 0
    IF ((AO1.EQ.'Empty').AND.(AO2.EQ.'Empty')) THEN
      CALL LSQUIT('Error in ls_getIntegrals1. DO_GRADIENT and AO1=AO2="Empty"',lupri)
    ELSEIF ((AO1.EQ.'Empty').OR.(AO2.EQ.'Empty')) THEN
      INT_INPUT%NDERIVP = 3
    ELSE
      INT_INPUT%NDERIVP = 6
    ENDIF
    INT_INPUT%NDERIVQ = 1
    INT_INPUT%DERIVORDER = max(INT_INPUT%derOrderP,INT_INPUT%derOrderQ)
    INT_INPUT%sameODs = .FALSE.
  CASE DEFAULT
    WRITE(LUPRI,'(1X,2A)') 'Error: Wrong case in ls_get_exchange_mat =',Spec
    CALL LSQUIT('Wrong case in ls_get_exchange_mat',lupri)
  END SELECT
  INT_INPUT%orderAngPrim = .FALSE.
  INT_INPUT%orderPQ      = .FALSE.

  call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)

  IF(INT_INPUT%DRHS_SYM .AND. Spec(1:7).EQ.'Regular')THEN
    ! SYMMETRIZE LSTENSOR
     call mem_alloc(TMP2,INT_INPUT%AOdim(2),INT_INPUT%AOdim(4),1,1,SETTING%nDmatRHS)
     call Build_full_5dim_from_lstensor(setting%Output%resultmat2,TMP2,&
          &INT_INPUT%AOdim(2),INT_INPUT%AOdim(4),1,1,SETTING%nDmatRHS)         
     call lstensor_free(setting%Output%resultmat2)
     call mem_alloc(TMP,INT_INPUT%AOdim(2),INT_INPUT%AOdim(4))
     DO IDMAT=1,SETTING%nDmatRHS
        DO J =1,INT_INPUT%AOdim(2)
           DO I =1,INT_INPUT%AOdim(4)
              TMP(I,J)= 0.5D0*TMP2(J,I,1,1,idmat)+0.5D0*TMP2(I,J,1,1,idmat)
           ENDDO
        ENDDO
        CALL DCOPY(setting%Output%ndim(1)*setting%Output%ndim(2),TMP,1,TMP2(:,:,1,1,idmat),1)
     ENDDO
     call mem_dealloc(TMP)
     call Build_lstensor_from_full_5dim(setting%Output%resultmat2,TMP2,Int_Input%AO(2)%p,Int_Input%AO(4)%p,&
          & Int_Input%AO(2)%p,Int_Input%AO(4)%p,INT_INPUT%AOdim(2),INT_INPUT%AOdim(4),1,1,SETTING%ndmatRHS,&
          & .TRUE.,.TRUE.,.FALSE.,.FALSE.,lupri)
     call mem_dealloc(TMP2)
  ENDIF
  IF(Spec(1:7).EQ.'Regular') call free_contrule(setting%output)
  CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
  CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)
  call ls_freeIntegralInputDensities(INT_INPUT)
  Call ls_freeDmatFromSetting(setting)
ENDIF
SETTING%SCHEME%FRAGMENT = FRAGMENT_SAVE
setting%output%ndim = ndim2

END SUBROUTINE ls_get_exchange_mat

!> \brief Calculate the coulomb matrix
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_get_coulomb_mat(AO1,AO2,AO3,AO4,Oper,intType,SETTING,LUPRI,LUERR)
implicit none
Character*(*)        :: AO1,AO2,AO3,AO4,Oper,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds
Logical             :: FRAGMENT_SAVE
integer :: ndim2(5)
CALL ls_setDefaultFragments(setting)
ndim2 = setting%output%ndim 

CALL ls_setDfull(setting)

FRAGMENT_SAVE = SETTING%SCHEME%FRAGMENT
SETTING%SCHEME%FRAGMENT = .FALSE.

CALL init_integral_input(INT_INPUT,SETTING)

INT_INPUT%DO_FOCK = .TRUE.
INT_INPUT%DO_Coulomb  = .TRUE.
INT_INPUT%exchangeFactor = SETTING%SCHEME%exchangeFactor

Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
CALL get_screening_matrices(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)           
setting%output%ndim = ndim2
INT_INPUT%operator = Oper

CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)


call init_lstensor_5dim(setting%output%resultmat2,Int_Input%AO(1)%p,Int_Input%AO(2)%p,&
     & Int_Input%AO(3)%p,Int_Input%AO(4)%p,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5),&
     & .TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,lupri)
!     & .TRUE.,.TRUE.,.FALSE.,.FALSE.,Int_input%OD_SCREEN,lupri)

CALL ls_setintegralinputdensities(INT_INPUT,SETTING)
call add_contrule(setting%output,1,2,0,0,5,0,0,0,0,0,3,4,0,0,0,1.d0,.true.)
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
call ls_freeIntegralInputDensities(INT_INPUT)
Call ls_freeDmatFromSetting(setting)
CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)

SETTING%SCHEME%FRAGMENT = FRAGMENT_SAVE
call free_contrule(setting%output)
CALL ls_freeDfull(setting)
setting%output%ndim = ndim2

END SUBROUTINE ls_get_coulomb_mat

!> \brief Calculate the coulomb and exchange matrix
!> \author T.Kjaergaard and S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_get_coulomb_and_exchange_mat(AO1,AO2,AO3,AO4,Oper,intType,SETTING,LUPRI,LUERR)
implicit none
Character*(*)        :: AO1,AO2,AO3,AO4,Oper,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: ndmat,nAObuilds,I
Logical             :: FRAGMENT_SAVE
integer :: ndim2(5)
CALL ls_setDefaultFragments(setting)
ndim2 = setting%output%ndim 

FRAGMENT_SAVE = SETTING%SCHEME%FRAGMENT
SETTING%SCHEME%FRAGMENT = .FALSE.
CALL init_integral_input(INT_INPUT,SETTING)
INT_INPUT%DO_FOCK = .TRUE.
INT_INPUT%DO_Coulomb  = .TRUE.
INT_INPUT%exchangeFactor = SETTING%SCHEME%exchangeFactor
INT_INPUT%DO_EXCHANGE = (INT_INPUT%exchangeFactor.GT.0.0d0)&
     & .OR. (SETTING%SCHEME%CAMbeta .GT. 0.0d0)

Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
CALL get_screening_matrices(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)           
setting%output%ndim = ndim2
INT_INPUT%operator = Oper
CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

call init_lstensor_5dim(setting%output%resultmat2,Int_Input%AO(1)%p,Int_Input%AO(2)%p,&
     & Int_Input%AO(3)%p,Int_Input%AO(4)%p,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5),&
     & .TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,lupri)

CALL ls_setintegralinputdensities(INT_INPUT,SETTING)
CALL initIntegralOutputDims(setting%Output,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,SETTING%ndmatRHS)
call add_contrule(setting%output,1,2,0,0,5,0,0,0,0,0,3,4,0,0,0,int_input%CoulombFactor,.true.)
call add_contrule(setting%output,1,3,0,0,5,0,0,0,0,0,2,4,0,0,0,-SETTING%SCHEME%exchangeFactor,.true.)
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
call ls_freeIntegralInputDensities(INT_INPUT)
Call ls_freeDmatFromSetting(setting)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)
call free_contrule(setting%output)

SETTING%SCHEME%FRAGMENT = FRAGMENT_SAVE
setting%output%ndim = ndim2

END SUBROUTINE ls_get_coulomb_and_exchange_mat

!> \brief Calculate the coulomb matrix using the jengine method
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_jengine(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
Integer              :: LUPRI,LUERR
Type(LSSETTING)      :: SETTING
Character*(*)        :: AO1,AO2,AO3,AO4,Oper,Spec,intType
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds

Real(realk),pointer  :: MatBlock(:,:,:,:,:)
Real(realk),pointer  :: MatBlock2(:,:)
Integer              :: nbast1,nbast2,nbast3,nbast4
Integer              :: start1,start2,start3,start4,ndim
Integer              :: I1,I2,I3,I4,I,J,nbast,iNode,jNode
Integer              :: end1,end2,end3,end4,idmatRHS,numFragments
real(realk),parameter :: ONE = 1.D0
Logical             :: permuteLHS,permuteRHS,FRAGMENT_SAVE
Logical             :: sameFragmentLHS,sameFragmentRHS
logical             :: l1,l2
integer             :: n1,n2,n3,n4,n5,nAtoms,iMol
!
!Integer              :: ndmatRHS
type(matrix),pointer :: JMat(:)
integer :: ndim2(5)
ndim2 = setting%output%ndim 

CALL ls_setDefaultFragments(setting)

IF(SETTING%SCHEME%FRAGMENT )THEN
allocate(Jmat(setting%nDmatRHS))
do I=1,setting%nDmatRHS
   call mat_init(Jmat(I),ndim2(1),ndim2(2))
   call mat_zero(Jmat(I))
enddo
SETTING%SCHEME%PUREFMM=.TRUE.
IF (SETTING%numNodes.NE.1) CALL LSQUIT('numNodes not equal to one in ls_jengine',lupri)
DO jNode=1,SETTING%numNodes
  SETTING%NODE = jNode
   CALL buildFragmentInfoAndBlocks(SETTING,AO1,AO2,AO3,AO4,LUPRI,LUERR)
   DO I=1,SETTING%FRAGMENTS%LHSblock%numBlocks
     permuteLHS = SETTING%FRAGMENTS%LHSblock%sameAOs.AND.&
     &            (.NOT.SETTING%FRAGMENTS%LHSblock%Blocks(I)%sameFragments)
     iNode = SETTING%FRAGMENTS%LHSblock%blocks(I)%node
     DO J=1,SETTING%FRAGMENTS%RHSblock%numBlocks
       permuteRHS = SETTING%FRAGMENTS%RHSblock%sameAOs.AND.&
     &              (.NOT.SETTING%FRAGMENTS%RHSblock%Blocks(J)%sameFragments)
       !Only J-block belonging to current node
       IF (SETTING%FRAGMENTS%RHSblock%blocks(J)%node.EQ.SETTING%NODE) THEN  

         CALL SetDaltonFragments(SETTING,I,J,sameFragmentLHS,sameFragmentRHS,lupri)
         CALL GetDaltonOrbitalInfo(SETTING,I,J,nbast1,nbast2,nbast3,nbast4,&
     &                             start1,start2,start3,start4,end1,end2,end3,end4)

!        Copy density-matrix sub-block
         IF (setting%LHSdmat) THEN
            call mem_alloc(setting%DfullLHS,nbast1,nbast2,setting%nDmatLHS)
            setting%LHSdalloc = .TRUE.
            CALL ls_mat_retrive_block(setting%DmatLHS,setting%DfullLHS,setting%ndmatLHS,&
     &                                nbast1,nbast2,start1,start2,permuteLHS)
            setting%LHSdfull = .TRUE.
         ENDIF                    
         IF (setting%RHSdmat) THEN
            call mem_alloc(setting%DfullRHS,nbast3,nbast4,setting%nDmatRHS)
            setting%RHSdalloc = .TRUE.
            CALL ls_mat_retrive_block(setting%DmatRHS,setting%DfullRHS,setting%ndmatRHS,&
     &                                nbast3,nbast4,start3,start4,permuteRHS)
            setting%RHSdfull = .TRUE.
         ENDIF

         call mem_alloc(MatBlock,nbast1,nbast2,1,1,setting%nDmatRHS)
         CALL LS_DZERO(MatBlock,nbast1*nbast2*setting%nDmatRHS)

!        *** CALCULATE INTEGRALS: matrix sub-block ***
         call initIntegralOutputDims(setting%Output,nbast1,nbast2,1,1,setting%nDmatRHS)
         CALL ls_jengine1(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
         CALL retrieve_Output(lupri,setting,MatBlock)
         setting%RHSdfull = .FALSE.
         setting%LHSdfull = .FALSE.
!        Copy matrix sub-block into full matrix
         DO idmatRHS = 1, setting%nDmatRHS
            MatBlock2 => MatBlock(:,:,1,1,iDmatRHS)
            call mat_add_block(Jmat(idmatRHS),MatBlock2,nbast1,nbast2,start1,start2)
         ENDDO
         IF (permuteLHS) THEN
            call mem_alloc(MatBlock2,nbast2,nbast1)
            DO idmatRHS = 1, setting%nDmatRHS
               DO I1=1,nbast1
                  DO I2=1,nbast2
                     MatBlock2(I2,I1) = MatBlock(I1,I2,1,1,idmatRHS)
                  ENDDO
               ENDDO
               call mat_add_block(Jmat(idmatRHS),MatBlock2,nbast2,nbast1,start2,start1)
            ENDDO
            call mem_dealloc(MatBlock2)
         ENDIF
         call mem_dealloc(MatBlock)

         CALL FreeDaltonFragments(SETTING)
         CALL ls_freeDblock(setting)
       ENDIF
     ENDDO
   ENDDO
   CALL freeFragmentInfoAndBlocks(SETTING)
!CAREFUL
ENDDO
Call ls_freeDmatFromSetting(setting)

setting%output%ndim = ndim2
call initIntegralOutputDims(setting%Output,ndim2(1),ndim2(2),1,1,ndim2(5))
call mem_alloc(MatBlock,ndim2(1),ndim2(2),1,1,ndim2(5))
if(setting%nDmatRHS.NE.ndim2(5))CALL LSQUIT('error ndim2(5) .NE.setting%nDmatRHS as they should be ',lupri)
do I=1,setting%nDmatRHS
   call mat_to_full(Jmat(I),1D0,MatBlock(:,:,1,1,I))
   call mat_free(Jmat(I))
enddo
deallocate(Jmat)
call add_contrule(setting%output,1,2,0,0,5,0,0,0,0,0,3,4,0,0,0,1.d0,.false.)
call ls_Build_lstensor_from_full_5dim(setting%output%resultmat2,setting,AO1,AO2,AO3,AO4,ndim2(5),intType,LUPRI,LUERR,MatBlock)
call free_contrule(setting%output)
call mem_dealloc(MatBlock)

IF(SETTING%SCHEME%FMM)THEN
   FRAGMENT_SAVE = SETTING%SCHEME%FRAGMENT
   SETTING%SCHEME%FRAGMENT = .FALSE.
   CALL ls_setDfull(setting)
   CALL ls_jengineClassical(AO1,AO2,AO3,AO4,Oper,intType,SETTING,LUPRI,LUERR)
   CALL ls_freeDfull(setting)
   SETTING%SCHEME%FRAGMENT = FRAGMENT_SAVE
ENDIF
SETTING%SCHEME%PUREFMM=.FALSE.
!CAREFUL
ELSE   !NO FRAGMENT
   CALL ls_setDfull(setting)
   CALL ls_jengine1(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
   CALL ls_freeDfull(setting)
ENDIF

END SUBROUTINE ls_jengine

!> \brief retrive a full block from a type matrix
!> \author S. Reine
!> \date 2010
!> \param Dmat the matrix type
!> \param Dblock the output block 
!> \param ndmat the number of matrices
!> \param n1 the number of functions in the small block dim 1
!> \param n2 the number of functions in the small block dim 2
!> \param s1 the start index of the small block dim 1
!> \param s2 the start index of the small block dim 2
!> \param permute if the permutational symmetry is used 
SUBROUTINE ls_mat_retrive_block(Dmat,Dblock,ndmat,n1,n2,s1,s2,permute)
implicit none
Integer       :: ndmat,n1,n2,s1,s2
type(matrixp) :: Dmat(ndmat)
Real(realk)   :: Dblock(n1,n2,ndmat)
Logical       :: permute
!
Integer :: idmat
DO idmat = 1, ndmat
  call mat_retrieve_block(Dmat(idmat)%p,Dblock(:,:,iDmat),n1,n2,s1,s2)
  IF (permute) THEN
    call ls_add_sym_block(Dmat(idmat)%p,Dblock(:,:,iDmat),n1,n2,s1,s2)
  ENDIF
ENDDO
END SUBROUTINE ls_mat_retrive_block

!> \brief add symmetry block from full block 
!> \author S. Reine
!> \date 2010
!> \param Dmat the matrix type
!> \param Dblock the output block 
!> \param n1 the number of functions in the small block dim 1
!> \param n2 the number of functions in the small block dim 2
!> \param s1 the start index of the small block dim 1
!> \param s2 the start index of the small block dim 2
SUBROUTINE ls_add_sym_block(Dmat,Dblock,n1,n2,s1,s2)
implicit none
TYPE(Matrix) :: Dmat
Integer      :: n1,n2,s1,s2
Real(realk)  :: Dblock(n1,n2)
!
Real(realk),pointer :: Dtrans(:,:)
Integer                 :: i1,i2

call mem_alloc(Dtrans,n2,n1)
call mat_retrieve_block(Dmat,Dtrans,n2,n1,s2,s1)
DO i2=1,n2
  DO i1=1,n1
     Dblock(i1,i2) = Dblock(i1,i2) + Dtrans(i2,i1)
  ENDDO
ENDDO
call mem_dealloc(Dtrans)
END SUBROUTINE ls_add_sym_block

!> \brief calculate the coulomb matrix using the jengine method
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_jengine1(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
!Real(realk)          :: Mat(:,:,:,:,:)
Character*(*)        :: AO1,AO2,AO3,AO4,Oper,Spec,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,idmat,nAtoms,iMol,nmat
integer :: ndim2(5)
ndim2 = setting%output%ndim 

nmat = setting%output%ndim(5)
CALL init_integral_input(INT_INPUT,SETTING)

Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
CALL get_screening_matrices(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
setting%output%ndim = ndim2
Int_input%DO_FMM = SETTING%SCHEME%FMM .AND. (Oper.EQ.'Coulomb')
Int_input%NonClassical_SCREEN = Int_input%DO_FMM
Int_input%OE_SCREEN = SETTING%SCHEME%OE_SCREEN .AND. (Oper.EQ.'Overlap')

INT_INPUT%operator = Oper

CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

INT_INPUT%sameODs      = .FALSE.
INT_INPUT%DO_FOCK      = .TRUE.
INT_INPUT%DO_COULOMB   = .TRUE.
INT_INPUT%DO_JENGINE   = .TRUE.
INT_INPUT%orderAngPrim = .FALSE.
INT_INPUT%orderPQ      = .FALSE.

CALL ls_setintegralinputdensities(INT_INPUT,SETTING)

SELECT CASE(Spec)
CASE('Regular')
   call init_lstensor_5dim(setting%output%resultmat2,Int_Input%AO(1)%p,Int_Input%AO(2)%p,Int_Input%AO(3)%p,&
!        & Int_Input%AO(4)%p,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5),.TRUE.,.TRUE.,.FALSE.,.FALSE.,Int_input%OD_SCREEN,lupri)
        & Int_Input%AO(4)%p,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5),.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,lupri)
   call add_contrule(setting%output,1,2,0,0,5,0,0,0,0,0,3,4,0,0,0,1.d0,.false.)
  !Default undifferentiated case
CASE('Gradient')
  setting%Output%dograd = .true.
  nAtoms =  setting%molecule(1)%p%nAtoms
  call init_gradientlstensor(setting%output%resultmat2,nAtoms,nmat,lupri)
  DO iMol=2,4
    IF (setting%molecule(iMol)%p%nAtoms.NE.nAtoms) CALL LSQUIT('Error in ls_jengine1. nAtoms inconsistency!',lupri)
  ENDDO
  !Molecular gradient
  INT_INPUT%DO_GRADIENT = .TRUE.
  INT_INPUT%sphericalEcoeff = .FALSE.
  INT_INPUT%derOrderP  = 1
  INT_INPUT%derOrderQ  = 0
  IF ((AO1.EQ.'Empty').AND.(AO2.EQ.'Empty')) THEN
    CALL LSQUIT('Error in ls_getIntegrals1. DO_GRADIENT and AO1=AO2="Empty"',lupri)
  ELSEIF ((AO1.EQ.'Empty').OR.(AO2.EQ.'Empty')) THEN
    INT_INPUT%NDERIVP = 3
  ELSE
    INT_INPUT%NDERIVP = 6
  ENDIF
  INT_INPUT%NDERIVQ = 1
  INT_INPUT%DERIVORDER = max(INT_INPUT%derOrderP,INT_INPUT%derOrderQ)
CASE DEFAULT
  CALL LSQUIT('Error in ls_jengine1. Wrong Spec case',lupri)
END SELECT

IF(SETTING%SCHEME%PUREFMM)Int_input%SKIP_FMM = .TRUE. 
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
call ls_freeIntegralInputDensities(INT_INPUT)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)
IF(Spec(1:7).EQ.'Regular') call free_contrule(setting%output)
setting%output%ndim = ndim2

END SUBROUTINE ls_jengine1

!> \brief calculate the classical part of the coulomb matrix using the jengine method (in a FMM sense)
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_jengineClassical(AO1,AO2,AO3,AO4,Oper,intType,SETTING,LUPRI,LUERR)
implicit none
Character*(*)        :: AO1,AO2,AO3,AO4,Oper,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR,ndmatRHS
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,idmat

CALL init_integral_input(INT_INPUT,SETTING)

Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN

Int_input%DO_FMM = SETTING%SCHEME%FMM .AND. (Oper.EQ.'Coulomb')
IF(Int_input%DO_FMM)THEN
   Int_input%NonClassical_SCREEN = Int_input%DO_FMM
   Int_input%OE_SCREEN = SETTING%SCHEME%OE_SCREEN .AND. (Oper.EQ.'Overlap')
   INT_INPUT%operator = Oper
   CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)
   INT_INPUT%sameODs      = .FALSE.
   INT_INPUT%DO_FOCK      = .TRUE.
   INT_INPUT%DO_COULOMB   = .TRUE.
   INT_INPUT%DO_JENGINE   = .TRUE.
   INT_INPUT%orderAngPrim = .FALSE.
   INT_INPUT%orderPQ      = .FALSE.
   
   CALL ls_setintegralinputdensities(INT_INPUT,SETTING)
   CALL initIntegralOutputDims(setting%Output,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,ndmatRHS)
   CALL JmatClassical(setting%Output,setting%Output%ndim(1),setting%Output%ndim(2),INT_INPUT)
   CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
   call ls_freeIntegralInputDensities(INT_INPUT)
   Call ls_freeDmatFromSetting(setting)
ENDIF

END SUBROUTINE ls_jengineClassical

!> \brief wrapper to calculate the screening matrices required for quadratic or better scaling
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param Input the integral input specification
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE get_screening_matrices(Input,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
implicit none
TYPE(INTEGRALINPUT)  :: Input
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!logical              :: family
Character*(*)        :: AO1,AO2,AO3,AO4,Oper
!
! Contracted screening matrix/matrices
IF(Input%CS_SCREEN) THEN
  IF(SETTING%SCHEME%INTPRINT .GT. 10) WRITE(LUPRI,'(2X,A)')&
          &'call GET_GAB_MATRIX - the  screening matrix (THE BOOK 9.12.26)'
  IF (Oper.EQ.'Nucrep') THEN
    CALL get_screen_matrix('LHS','Contracted',Input,AO1,AO2,'Coulomb',SETTING,LUPRI,LUERR)
    CALL get_screen_matrix('RHS','Contracted',Input,AO3,AO4,'Nuclei',SETTING,LUPRI,LUERR)
  ELSE
    CALL get_screen_matrix('LHS','Contracted',Input,AO1,AO2,Oper,SETTING,LUPRI,LUERR)
    CALL get_screen_matrix('RHS','Contracted',Input,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
  ENDIF
ENDIF
!  
! Primitive screening matrix/matrices
IF(Input%PS_SCREEN)THEN
  IF (Oper.EQ.'Nucrep') THEN
    CALL get_screen_matrix('LHS','Primitive',Input,AO1,AO2,'Coulomb',SETTING,LUPRI,LUERR)
    CALL get_screen_matrix('RHS','Primitive',Input,AO3,AO4,'Nuclei',SETTING,LUPRI,LUERR)
  ELSE
    CALL get_screen_matrix('LHS','Primitive',Input,AO1,AO2,Oper,SETTING,LUPRI,LUERR)
    CALL get_screen_matrix('RHS','Primitive',Input,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
  ENDIF
ENDIF
!
END SUBROUTINE get_screening_matrices

!> \brief free the calculated screening matrices
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param Input the integral input specification
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE free_screening_matrices(Input,SETTING,LUPRI,LUERR)
implicit none
TYPE(INTEGRALINPUT)  :: Input
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
IF(Input%CS_SCREEN) THEN
  IF(.NOT.Input%CS_RHSGABusePointer)THEN
     call lstensor_free(input%LST_GAB_RHS)
     deallocate(input%LST_GAB_RHS)
  ENDIF
  nullify(input%LST_GAB_RHS)
  call lstensor_free(input%LST_GAB_LHS)
  deallocate(input%LST_GAB_LHS)
  nullify(input%LST_GAB_LHS)
ENDIF
IF(Input%PS_SCREEN)THEN
  IF(.NOT.Input%PS_RHSGABusePointer)THEN
     call lstensor_free(input%LST_pGAB_RHS)
     deallocate(input%LST_pGAB_RHS)
  ENDIF
  nullify(input%LST_pGAB_RHS)
  call lstensor_free(input%LST_pGAB_LHS)
  deallocate(input%LST_pGAB_LHS)
  nullify(input%LST_pGAB_LHS)
ENDIF
END SUBROUTINE free_screening_matrices

!> \brief calculate the screening matrix
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param IntegralSide The side of the screening matrix (LHS or RHS) 
!> \param intType the label for primitive or contracted calc
!> \param Input the integral input specification
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1/3
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2/4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE get_screen_matrix(integralSide,intType,Input,AO1,AO2,Oper,SETTING,LUPRI,LUERR)
implicit none
TYPE(INTEGRALINPUT),target  :: Input
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
Character*(*)        :: AO1,AO2,Oper,integralSide,intType
!
type(lstensor)       :: TMP_lstensor 
type(lstensor),pointer :: GAB
Integer             :: iprint,nc1,nc2,np1,np2,n1,n2,i1,i2,s1,s2
Logical             :: uncont,intnorm,nuclei,LHS,RHS
Real(realk)         :: maxGabElm,maxgabelm2
Character(80)       :: Filename
Character(49)       :: identifier
Real(realk),pointer :: screenMat(:,:)
TYPE(MOLECULE_PT),pointer :: FRAGMENTS(:)
integer              :: f1,f2,iAO,AOA,AOB
logical             :: RHSGABusePointer
!
NULLIFY(FRAGMENTS)
ALLOCATE(FRAGMENTS(4))
iprint = SETTING%SCHEME%intPrint

LHS = integralSide.EQ.'LHS'
AOA=1
AOB=2
RHS = integralSide.EQ.'RHS'
IF(RHS)AOA=3
IF(RHS)AOB=4

IF (.NOT.(LHS.OR.RHS)) THEN
  WRITE(LUPRI,'(1X,2A)') 'Input error in get_screen_matrix: integralSide = ',integralSide
  CALL LSQUIT('Input error in get_screen_matrix',lupri)
ENDIF
IF(SETTING%SCHEME%FRAGMENT)THEN
    DO iAO=1,setting%nAO
      FRAGMENTS(iAO)%p => SETTING%FRAGMENT(iAO)%p
    ENDDO
    IF (LHS) THEN
       FRAGMENTS(3)%p => FRAGMENTS(1)%p
       FRAGMENTS(4)%p => FRAGMENTS(2)%p
       f1 = 1
       f2 = 2
       s1 = SETTING%FRAGMENTS%LHSblock%blocks(SETTING%FRAGMENTS%iLHSblock)%startOrb1-1
       s2 = SETTING%FRAGMENTS%LHSblock%blocks(SETTING%FRAGMENTS%iLHSblock)%startOrb2-1
    ELSE
       FRAGMENTS(1)%p => FRAGMENTS(3)%p
       FRAGMENTS(2)%p => FRAGMENTS(4)%p
       f1 = 3
       f2 = 4
       s1 = SETTING%FRAGMENTS%RHSblock%blocks(SETTING%FRAGMENTS%iRHSblock)%startOrb1-1
       s2 = SETTING%FRAGMENTS%RHSblock%blocks(SETTING%FRAGMENTS%iRHSblock)%startOrb2-1
    ENDIF
    n1 = getNbasis(AO1,intType,FRAGMENTS(f1)%p,LUPRI)
    n2 = getNbasis(AO2,intType,FRAGMENTS(f2)%p,LUPRI)
    IF (AO1.EQ.'Nuclear') THEN
      s1 = n1
    ENDIF
    IF (AO2.EQ.'Nuclear') THEN
      s2 = n2
    ENDIF
ELSE
   n1 = getNbasis(AO1,intType,SETTING%MOLECULE(AOA)%p,LUPRI)
   n2 = getNbasis(AO2,intType,SETTING%MOLECULE(AOB)%p,LUPRI)
   IF((SETTING%batchindex(AOA).NE.0).AND.SETTING%molID(AOA).NE.0)THEN
      WRITE(LUPRI,'(1X,A)') 'Input error in get_screen_matrix: only made for either batchindex use og molID use '
      CALL LSQUIT('Input error in get_screen_matrix',lupri)      
   ENDIF
   s1 = SETTING%molID(AOA) !not 0, when the molecule pointers are different molecules
   s2 = SETTING%molID(AOB)
   IF(SETTING%batchindex(AOA).NE.0)THEN
      s1=SETTING%batchindex(AOA) !not 0, only when a single AObatch is computed
      IF (intType.EQ.'Primitive') THEN
         CALL LSQUIT('Primitive screening using batchindex not implemented',lupri)
      ELSE
         n1=SETTING%batchdim(AOA) 
      ENDIF
   ENDIF
   IF(SETTING%batchindex(AOB).NE.0)THEN
      s2=SETTING%batchindex(AOB)
      IF (intType.EQ.'Primitive') THEN
         CALL LSQUIT('Primitive screening using batchindex not implemented',lupri)
      ELSE
         n2=SETTING%batchdim(AOB) 
      ENDIF
   ENDIF
ENDIF
!
IF (iprint.gt.5) THEN
  WRITE(LUPRI,'(1X,A)') 'Output from get_screen_matrix'
  WRITE(LUPRI,'(3X,2A)')   'AO1:          ',AO1
  WRITE(LUPRI,'(3X,2A)')   'AO2:          ',AO2
  WRITE(LUPRI,'(3X,2A)')   'intType:      ',intType
  WRITE(LUPRI,'(3X,2A)')   'integralSide: ',integralSide
  WRITE(LUPRI,'(3X,A,I3)') 'n1:           ',n1
  WRITE(LUPRI,'(3X,A,I3)') 'n2:           ',n2
  WRITE(LUPRI,'(3X,A,I3)') 's1:           ',s1
  WRITE(LUPRI,'(3X,A,I3)') 's2:           ',s2
ENDIF
!
!determine if we can use a pointer to the LHS GAB for the RHS GAB 
RHSGABusePointer = .TRUE.
!we only do this for the right hand side. 
!So we always build for LHS and maybe set pointer for RHS
IF(LHS)RHSGABusePointer = .FALSE.
IF(.NOT.Input%sameODs)RHSGABusePointer = .FALSE.
IF (intType.EQ.'Primitive') THEN
   IF(.NOT.ASSOCIATED(input%LST_pGAB_LHS))RHSGABusePointer = .FALSE.
ELSE
   IF(.NOT.ASSOCIATED(input%LST_GAB_LHS))RHSGABusePointer = .FALSE.
ENDIF
IF(Oper.EQ.'Nuclei')RHSGABusePointer = .FALSE.

IF (intType.EQ.'Primitive') THEN
   IF(LHS)THEN
      allocate(input%LST_pGAB_LHS)
      GAB => input%LST_pGAB_LHS
   ELSE
      allocate(input%LST_pGAB_RHS)
      GAB => input%LST_pGAB_RHS
   ENDIF
ELSE
   IF(LHS)THEN
      allocate(input%LST_GAB_LHS)
      GAB => input%LST_GAB_LHS
   ELSE
      allocate(input%LST_GAB_RHS)
      GAB => input%LST_GAB_RHS
   ENDIF
ENDIF

IF(RHSGABusePointer)THEN
   IF (intType.EQ.'Primitive') THEN
      input%LST_pGAB_RHS => input%LST_pGAB_LHS
      Input%PS_MAXELM_RHS = Input%PS_MAXELM_LHS 
      Input%PS_RHSGABusePointer = .TRUE.
   ELSE
      input%LST_GAB_RHS => input%LST_GAB_LHS
      Input%CS_MAXELM_RHS = Input%CS_MAXELM_LHS
      Input%CS_RHSGABusePointer = .TRUE.
   ENDIF
ELSE
   IF (intType.EQ.'Primitive') THEN
      Input%PS_RHSGABusePointer = .FALSE.
   ELSE
      Input%CS_RHSGABusePointer = .FALSE.
   ENDIF
   write(identifier,'(A3,A22,A1,A22,A1)') 'CS_',Setting%FRAGMENT(AOA)%p%label,'_',Setting%FRAGMENT(AOB)%p%label,'_'
   CALL io_get_filename(Filename,identifier,AO1,AO2,AO1,AO2,s1,s2,Oper,intType,SETTING%SCHEME%FRAGMENT,LUPRI,LUERR)
   IF (io_file_exist(Filename,SETTING%IO))THEN
      CALL io_read_lstensor(GAB,Filename,SETTING%IO,LUPRI,LUERR)
   ELSE
      IF ((AO1.EQ.'Empty').AND.(AO2.EQ.'Empty')) THEN
         call mem_alloc(screenMat,n1,n2)
         screenMat(1,1) = 1.0d0
         call build_Nuclearlstensor(ScreenMat,GAB,1)
         call mem_dealloc(screenMat)
      ELSEIF (Oper.EQ.'Nuclei') THEN
         call mem_alloc(screenMat,n1,n2)
         CALL ls_getNucScreenIntegrals(ScreenMat,Filename,AO1,AO2,Oper,intType,SETTING,INPUT,LUPRI,LUERR,LHS)
         call mem_dealloc(screenMat)
      ELSE
         call initIntegralOutputDims(setting%Output,n1,n2,1,1,1)
         CALL ls_getScreenIntegrals(Filename,AO1,AO2,Oper,intType,SETTING,LUPRI,LUERR,LHS)
         IF(intType.EQ.'Primitive') THEN
            call free_primassistarrays(setting%Output,setting%Output%resultmat2)
         ENDIF
         CALL retrieve_output(lupri,setting,GAB) 
      ENDIF
      CALL io_add_filename(SETTING%IO,Filename,LUPRI)
      CALL io_write_lstensor(GAB,Filename,SETTING%IO,LUPRI,LUERR)
   ENDIF
   call determine_maxelm(maxGABElm,GAB)
   IF (intType.EQ.'Primitive') THEN
      IF(LHS)THEN
         Input%PS_MAXELM_LHS = maxGabElm
      ELSE
         Input%PS_MAXELM_RHS = maxGabElm
      ENDIF
   ELSE
      IF(LHS)THEN
         Input%CS_MAXELM_LHS = maxGabElm
      ELSE
         Input%CS_MAXELM_RHS = maxGabElm
      ENDIF
   ENDIF
ENDIF

DEALLOCATE(FRAGMENTS)
END SUBROUTINE get_screen_matrix

!> \brief calculate the screening integrals
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param Filename the filename used to store the screening matrix
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1/3
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2/4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param LHSGAB if this is a LHS screening matrix
SUBROUTINE ls_getScreenIntegrals(Filename,AO1,AO2,Oper,intType,SETTING,LUPRI,LUERR,LHSGAB)
implicit none
!Real(realk)          :: Integrals(:,:)
Character*(*)        :: Filename
Character*(*)        :: AO1,AO2,Oper,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
LOGICAL              :: LHSGAB
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds
Real(realk)          :: tstart,tend
Integer              :: iAtom,iFragment,iAO
TYPE(MOLECULE_PT),pointer   :: FRAGMENTS(:)
integer :: ndim2(5)
!--------------------------------------
TYPE(INTEGRALINPUT)  :: CONT_INPUT
TYPE(AOITEM),target  :: CONTAObuild(4)
Integer              :: CONTnAObuilds
!--------------------------------------

ndim2 = setting%output%ndim 
NULLIFY(FRAGMENTS)
ALLOCATE(FRAGMENTS(4))
CALL LSTIMER('START ',tstart,tend,lupri)
CALL ls_setDefaultFragments(setting)
CALL init_integral_input(INT_INPUT,SETTING)
INT_INPUT%operator = Oper
INT_INPUT%CS_int=.TRUE.
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
INT_INPUT%CS_SCREEN = .FALSE.
CALL SetInputAO(INT_INPUT,AO1,AO2,AO1,AO2,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
CALL initIntegralOutputDims(setting%Output,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,1)

IF (Oper.EQ.'Nuclei') THEN
   CALL LSQUIT('error in getScreenIntegrals',lupri)
ELSEIF ((AO1.EQ.'Empty').AND.(AO2.EQ.'Empty')) THEN
   CALL LSQUIT('result is screenMat(1,1) = 1.0d0 but this should have been done a different place',lupri)
!    Integrals(1,1) = 1.0d0
ELSE
   IF(intType.EQ.'Primitive')THEN
!      call init_lstensor_5dim(setting%output%resultmat2,Int_Input%AO(1)%p,Int_Input%AO(2)%p,Int_Input%AO(1)%p,&
!           & Int_Input%AO(1)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,1,.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,lupri)
      !make contracted AO, which is used because we want to store the primlstensor in a special format 
      CALL init_integral_input(CONT_INPUT,SETTING)
      CONT_INPUT%operator = Oper
      CONT_INPUT%CS_int=.TRUE.
      CONT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
      CONT_INPUT%CS_SCREEN = .FALSE.
      CALL SetInputAO(CONT_INPUT,AO1,AO2,AO1,AO2,'Contracted',CONTAObuild,CONTnAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
      Int_Input%PS_int = .TRUE.
      call init_primGablstensor(setting%output,setting%output%resultmat2,Int_Input%AO(1)%p,&
           &Int_Input%AO(1)%p,CONT_Input%AO(1)%p,CONT_Input%AO(2)%p,lupri)
      CALL FreeInputAO(CONTAObuild,CONTnAObuilds,LUPRI)
   ELSE
      call init_lstensor_5dim(setting%output%resultmat2,Int_Input%AO(1)%p,Int_Input%AO(2)%p,Int_Input%AO(1)%p,&
           & Int_Input%AO(1)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,1,.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,lupri)
   ENDIF
  call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
ENDIF

CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
!DEALLOCATE(setting%Output%ResultMat)

CALL LSTIMER(Filename(1:15),tstart,tend,lupri)
DEALLOCATE(FRAGMENTS)
NULLIFY(FRAGMENTS)
setting%output%ndim = ndim2
END SUBROUTINE ls_getScreenIntegrals

!> \brief calculate the nuclear screening integrals
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param Integrals the output matrix of nuclear screening integrals
!> \param Filename the filename used to store the screening matrix
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1/3
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2/4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param OUTPUT_INPUT the integral input 
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param LHSGAB if this is a LHS screening matrix
SUBROUTINE ls_getNucScreenIntegrals(Integrals,Filename,AO1,AO2,Oper,intType,SETTING,OUTPUT_INPUT,LUPRI,LUERR,LHSGAB)
implicit none
Real(realk)          :: Integrals(:,:)
Character*(*)        :: Filename
Character*(*)        :: AO1,AO2,Oper,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
LOGICAL              :: LHSGAB
TYPE(INTEGRALINPUT)  :: OUTPUT_INPUT
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,AOA,AOB
Real(realk)          :: tstart,tend
Integer              :: iAtom,iFragment,iAO,natoms
TYPE(MOLECULE_PT),pointer   :: FRAGMENTS(:)
CALL ls_setDefaultFragments(setting)
NULLIFY(FRAGMENTS)
ALLOCATE(FRAGMENTS(4))
IF(LHSGAB)THEN
   AOA=1
   AOB=2
ELSE
   AOA=3
   AOB=4
ENDIF
CALL LSTIMER('START ',tstart,tend,lupri)
CALL init_integral_input(INT_INPUT,SETTING)
INT_INPUT%operator = Oper
INT_INPUT%CS_int=.TRUE.
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
INT_INPUT%CS_SCREEN = .FALSE.
CALL SetInputAO(INT_INPUT,AO1,AO2,AO1,AO2,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
!CALL initIntegralOutput(Int_Output,INT_INPUT%AOdim(AOA),INT_INPUT%AOdim(AOB),1,1,1)

IF(SETTING%SCHEME%FRAGMENT)THEN
   DO iAO=1,setting%nAO
      FRAGMENTS(iAO)%p => SETTING%FRAGMENT(iAO)%p
   ENDDO
   IF (LHSGAB) THEN
      FRAGMENTS(3)%p => FRAGMENTS(1)%p
      FRAGMENTS(4)%p => FRAGMENTS(2)%p
      IF (AO1.EQ.'Empty') THEN
         iFragment = 2
      ELSEIF (AO2.EQ.'Empty') THEN
         iFragment = 1
      ELSE
         CALL LSQUIT('Error in ls_getScreenIntegrals, FRAGMENT,LHSGAB',lupri)
      ENDIF
   ELSE
      FRAGMENTS(1)%p => FRAGMENTS(3)%p
      FRAGMENTS(2)%p => FRAGMENTS(4)%p
      IF (AO1.EQ.'Empty') THEN
         iFragment = 4
      ELSEIF (AO2.EQ.'Empty') THEN
         iFragment = 3
      ELSE
         CALL LSQUIT('Error in ls_getScreenIntegrals, FRAGMENT,RHSGAB',lupri)
      ENDIF
   ENDIF
   DO iAtom=1,FRAGMENTS(iFragment)%p%nAtoms
      Integrals(iAtom,1) = FRAGMENTS(iFragment)%p%ATOM(iAtom)%Charge
   ENDDO
   nAtoms = FRAGMENTS(iFragment)%p%nAtoms
   IF (intType.EQ.'Primitive') THEN
      IF(LHSGAB)THEN
         call build_Nuclearlstensor(Integrals,OUTPUT_INPUT%LST_pGAB_LHS,nAtoms)
      ELSE
         call build_Nuclearlstensor(Integrals,OUTPUT_INPUT%LST_pGAB_RHS,nAtoms)
      ENDIF
   ELSE
      IF(LHSGAB)THEN
         call build_Nuclearlstensor(Integrals,OUTPUT_INPUT%LST_GAB_LHS,nAtoms)
      ELSE
         call build_Nuclearlstensor(Integrals,OUTPUT_INPUT%LST_GAB_RHS,nAtoms)
      ENDIF         
   ENDIF
ELSE
   DO iAtom=1,SETTING%MOLECULE(1)%p%nAtoms
      Integrals(iAtom,1) = SETTING%MOLECULE(1)%p%ATOM(iAtom)%Charge
   ENDDO
   nAtoms = SETTING%MOLECULE(1)%p%nAtoms
   IF (intType.EQ.'Primitive') THEN
      IF(LHSGAB)THEN
         call build_Nuclearlstensor(Integrals,OUTPUT_INPUT%LST_pGAB_LHS,nAtoms)
      ELSE
         call build_Nuclearlstensor(Integrals,OUTPUT_INPUT%LST_pGAB_RHS,nAtoms)
      ENDIF
   ELSE
      IF(LHSGAB)THEN
         call build_Nuclearlstensor(Integrals,OUTPUT_INPUT%LST_GAB_LHS,nAtoms)
      ELSE
         call build_Nuclearlstensor(Integrals,OUTPUT_INPUT%LST_GAB_RHS,nAtoms)
      ENDIF         
   ENDIF
ENDIF

CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
!DEALLOCATE(Int_Output%ResultMat)

CALL LSTIMER(Filename(1:15),tstart,tend,lupri)
DEALLOCATE(FRAGMENTS)

END SUBROUTINE ls_getNucScreenIntegrals

!> \brief Sets up the input AOs
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param INT_INPUT the integral input 
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param intType the label for primitive or contracted calc
!> \param AObuild the list of AOITEMS to be build
!> \param nAObuilds the number of AOITEMS build
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param LHSGAB if this is a LHS screening matrix
SUBROUTINE SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
implicit none
Character*(*)        :: AO1,AO2,AO3,AO4,intType
Type(LSSETTING)      :: SETTING
TYPE(INTEGRALINPUT)  :: INT_INPUT
Integer              :: LUPRI,LUERR,nAObuilds
TYPE(AOITEM),target  :: AObuild(4)
TYPE(MOLECULE_PT),pointer  :: FRAGMENTS(:)
LOGICAL              :: LHSGAB !THIS ONLY HAS AN EFFECT WHEN USING FRAGMENTS 
!
Integer                    :: IAO,JAO,indAO
Character(80)              :: AOstring(4),AOtype
TYPE(BASISSETINFO),pointer :: AObasis
Logical                    :: uniqueAO,emptyAO,uncont,intnrm,sameFrag(4,4)
Integer                    :: ndim(4),indexUnique(4),AObatchdim,batchindex(4)

IF (setting%nAO.NE.4) CALL LSQUIT('Error in SetInputAO. nAO .NE. 4',lupri)
NULLIFY(FRAGMENTS)
ALLOCATE(FRAGMENTS(4))

!Used for testing if the AOs are identical (depends also in the AO-strings)
sameFrag = setting%sameFrag
batchindex = SETTING%batchindex
AOtype = 'Default'
IF (INT_INPUT%operator .EQ. 'Mulmom') AOtype = 'Cartesian'

AOstring(1) = AO1
AOstring(2) = AO2
AOstring(3) = AO3
AOstring(4) = AO4

DO iAO=1,setting%nAO
  FRAGMENTS(iAO)%p => SETTING%FRAGMENT(iAO)%p
ENDDO

INT_INPUT%sameLHSaos = (AO1.EQ.AO2) .AND. (.NOT. AO1.EQ.'Empty').AND.samefrag(1,2)
INT_INPUT%sameRHSaos = (AO3.EQ.AO4) .AND. (.NOT. AO3.EQ.'Empty').AND.samefrag(3,4)
INT_INPUT%sameODs    = (AO1.EQ.AO3) .AND. (AO2.EQ.AO4).AND.samefrag(1,3).AND.samefrag(2,4)

!Specical settings for Cauchy-Schwarz screening integrals
IF(INT_INPUT%CS_int)THEN
   IF(LHSGAB)THEN
      FRAGMENTS(3)%p => FRAGMENTS(1)%p
      FRAGMENTS(4)%p => FRAGMENTS(2)%p
      AOstring(3) = AO1
      AOstring(4) = AO2
      sameFrag(3,4) = sameFrag(1,2)
      sameFrag(4,3) = sameFrag(2,1)
      INT_INPUT%sameRHSaos = INT_INPUT%sameLHSaos
      batchindex(3) = SETTING%batchindex(1) 
      batchindex(4) = SETTING%batchindex(2) 
   ELSE
      FRAGMENTS(1)%p => FRAGMENTS(3)%p
      FRAGMENTS(2)%p => FRAGMENTS(4)%p
      AOstring(1) = AO3
      AOstring(2) = AO4
      sameFrag(1,2) = sameFrag(3,4)
      sameFrag(2,1) = sameFrag(4,3)
      INT_INPUT%sameLHSaos = INT_INPUT%sameRHSaos
      batchindex(1) = SETTING%batchindex(3) 
      batchindex(2) = SETTING%batchindex(4) 
   ENDIF
   sameFrag(1,3) = .TRUE.
   sameFrag(2,4) = .TRUE.
   sameFrag(3,1) = .TRUE.
   sameFrag(4,2) = .TRUE.
   IF(sameFrag(1,2))then
      sameFrag(1,4) = .TRUE.
      sameFrag(2,3) = .TRUE.
      sameFrag(3,2) = .TRUE.
      sameFrag(4,1) = .TRUE.
   ELSE
      sameFrag(1,4) = .FALSE.
      sameFrag(2,3) = .FALSE.
      sameFrag(3,2) = .FALSE.
      sameFrag(4,1) = .FALSE.
   ENDIF
   INT_INPUT%sameODs = .TRUE.
ENDIF
!ENDIF
   
IF (intType.EQ.'Primitive') THEN
  uncont = .TRUE.
  intnrm = .TRUE.
ELSEIF (intType.EQ.'Contracted') THEN
  uncont = .FALSE.
  intnrm = .FALSE.
ELSE
  WRITE(LUPRI,'(1X,2A)') 'Wrong case in SetInputAO, intType =',intType
  CALL LSQUIT('Error - wrong intType in SetInputAO',lupri)
ENDIF

!Simen Should be set another place!
IF (AOstring(3).NE.'Regular'.OR.AOstring(4).NE.'Regular') INT_INPUT%CoulombFactor = 1.0d0

nAObuilds = 0
DO IAO=1,4
   uniqueAO = .true.
   DO JAO=1,IAO-1
      IF ((AOstring(IAO).EQ.AOstring(JAO)).AND.sameFrag(iAO,jAO)) THEN
         uniqueAO = .false.
         indAO=indexUnique(JAO)
        EXIT
     ENDIF
  ENDDO
  IF (SETTING%SCHEME%FRAGMENT) uniqueAO = .true.
  IF (uniqueAO) THEN
    nAObuilds = nAObuilds+1
    indAO   = nAObuilds
    indexUnique(IAO) = indAO
    CALL SetAObatch(AObuild(indAO),batchindex(iAO),ndim(indAO),AOstring(iAO),intType,AOtype,SETTING%Scheme,&
     &              FRAGMENTS(iAO)%p,Setting%BASIS(iAO)%p,LUPRI,LUERR)
    IF ((AOstring(IAO).EQ.'Nuclear').OR.(AOstring(IAO).EQ.'Nuclei')) THEN
      INT_INPUT%sameLHSaos = .FALSE.
      INT_INPUT%sameRHSaos = .FALSE.
    ENDIF
  ENDIF
  Int_Input%AO(IAO)%p => AObuild(indAO)
  Int_Input%AOdim(IAO) = ndim(indAO)
ENDDO

DEALLOCATE(FRAGMENTS)

END SUBROUTINE SetInputAO

!> \brief Sets up the AO-batch
!> \author S. Reine
!> \date 18-03-2010
!> \param AObatch The AO-batch
!> \param AO Specifying what basis set to use: 'Regular', 'DF-Aux', ...
!> \param intType Specifying contracted or primitive basis
!> \param AOtype Specifying basis-function type: 'Hermite', 'Cartesian', 'Deafult'
!> \param Scheme Specifies integral scheme
!> \param Molecule The molecule
!> \param Basis The basis set
!> \param LUPRI Default output unit
!> \param LUERR Deafult error unit
SUBROUTINE SetAObatch(AObatch,batchindex,nDim,AO,intType,AOtype,Scheme,Molecule,Basis,LUPRI,LUERR)
implicit none
TYPE(AOITEM),intent(OUT)          :: AObatch
Integer,intent(IN)                :: batchindex
Integer,intent(OUT)               :: nDim
Character*(*),intent(IN)          :: AO,intType,AOtype
Type(LSINTSCHEME),intent(IN)      :: Scheme
Type(MOLECULEINFO),intent(IN)     :: Molecule
Type(BASISINFO),intent(IN),target :: Basis
Integer,intent(IN)                :: LUPRI,LUERR
!
TYPE(BASISSETINFO),pointer :: AObasis
Logical :: uncont,intnrm,emptyAO
integer :: AObatchdim

IF (intType.EQ.'Primitive') THEN
  uncont = .TRUE.
  intnrm = .TRUE.
ELSEIF (intType.EQ.'Contracted') THEN
  uncont = .FALSE.
  intnrm = .FALSE.
ELSE
  WRITE(LUPRI,'(1X,2A)') 'Wrong case in SetAObatch, intType =',intType
  CALL LSQUIT('Error - wrong intType in SetAObatch',lupri)
ENDIF

emptyAO = .false.
SELECT CASE(AO)
CASE ('Regular')
  AObasis => Basis%REGULAR
CASE ('Huckel')
  AObasis => Basis%HUCKEL
CASE ('DF-Aux')
  AObasis => Basis%AUXILIARY
CASE ('Empty')
  emptyAO = .true.
  CALL BUILD_EMPTY_AO(AObatch,LUPRI)
  nDim = 1
CASE ('Nuclear')
  emptyAO = .true.
  CALL BUILD_EMPTY_NUCLEAR_AO(AObatch,Molecule,LUPRI)
  nDim = 1
CASE ('Nuclei')
  emptyAO = .true.
  CALL BUILD_EMPTY_NUCLEAR_AO(AObatch,Molecule,LUPRI)
  nDim = Molecule%nAtoms
CASE DEFAULT
  print*,'case: ',AO
  WRITE(lupri,*) 'case: ',AO
  WRITE(luerr,*) 'case: ',AO
  CALL LSQuit('Programming error: Not a case in SetAObatch!',lupri)
END SELECT
IF (.not.emptyAO) THEN
   IF(batchindex.EQ.0)THEN
      CALL BUILD_AO(LUPRI,SCHEME,SCHEME%AOPRINT,&
           &              Molecule,AObasis,AObatch,AOtype,&
           &              uncont,intnrm)
      nDim = getNbasis(AO,intType,Molecule,LUPRI)
   ELSE
      CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SCHEME,&
           & SCHEME%AOPRINT,molecule,AObasis,AObatch,&
           & AOtype,uncont,intnrm,batchindex,AObatchdim)
      nDim = AObatchdim
   ENDIF
ENDIF

END SUBROUTINE SetAObatch

!> \brief frees the input AOs
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param AObuild the list of AOITEMS to be build
!> \param nAObuilds the number of AOITEMS build
!> \param LUPRI logical unit number of the Default output file
SUBROUTINE FreeInputAO(AObuild,nAObuilds,LUPRI)
implicit none
Integer              :: LUPRI,nAObuilds
TYPE(AOITEM),target  :: AObuild(4)
!
integer :: iAO

DO iAO=1,nAObuilds
  CALL free_aoitem(lupri,AObuild(iAO))
ENDDO
 
END SUBROUTINE FreeInputAO

!> \brief calculates the multipole moments required for FMM
!> \author T. Kjaergaard
!> \date 2010
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param SETTING Integral evalualtion settings
!> \param nbast the number of basis functions 
!> \param nbastaux the number of auxillary basis functions 
!> \param ndim1 the size of the 1. dimension 
!> \param ndim2 the size of the 2. dimension 
!> \param ndim3 the size of the 3. dimension 
!> \param ndim4 the size of the 4. dimension 
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param COULOMB if this is a coulomb calc (or overlap)
SUBROUTINE ls_multipolemoment(LUPRI,LUERR,SETTING,nbast,nbastaux,ndim1,ndim2,ndim3,ndim4,AO1,AO2,AO3,AO4,Oper,intType,COULOMB)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,ndim1,ndim2,ndim3,ndim4
TYPE(LSSETTING)       :: SETTING
Character*(*)         :: AO1,AO2,AO3,AO4,Oper,intType
!TYPE(MATRIX)          :: Dmat,DmatRHS
Real(realk),pointer   :: Dfull(:,:),RHSDENS(:,:)
Real(realk),pointer   :: integrals(:,:,:,:,:)
INTEGER               :: nderiv,nmat,nsphmat,nbast,I,nbastaux
LOGICAL               :: LHSDENSFIT,RHSDENSFIT,FOURcenter
TYPE(AOITEM),target   :: AObuild(4)
Integer               :: nAObuilds,LU_DENS,LUINTM2
Integer               :: READA,READB,READC,READD,READE,READF
LOGICAL               :: L1,L2
Integer               :: MMunique_ID1,MMunique_ID2,IDUMMY,N,T,J,IDER,start1,start2
Integer               :: nrowLHS,ncolLHS,nrowRHS,ncolRHS
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
real(realk),pointer :: WORK(:)
real(realk),parameter :: TWO=2.D0,HALF=0.5D0
LOGICAL               :: COULOMB,fragment2
REAL(REALK)           :: TS,TE
IDUMMY=1
LHSDENSFIT = .FALSE.
RHSDENSFIT = .FALSE.
FOURcenter = .FALSE.
IF(AO1(1:6) .EQ.'DF-Aux') LHSDENSFIT = .TRUE.
IF(AO3(1:6) .EQ.'DF-Aux') RHSDENSFIT = .TRUE.
IF(.NOT.LHSDENSFIT .AND. .NOT.RHSDENSFIT)THEN
   IF(AO1(1:7) .EQ. 'Regular' .AND. AO2(1:7) .EQ. 'Regular')THEN
      IF(AO3(1:7) .EQ. 'Regular' .AND. AO4(1:7) .EQ. 'Regular') FOURcenter = .TRUE.
   ENDIF
ENDIF
!print*,'LHSDENSFIT',LHSDENSFIT,'RHSDENSFIT',RHSDENSFIT,'FOUR',FOURcenter

IF(SETTING%SCHEME%NO_MMFILES)THEN
   WRITE(LUPRI,*)'You have chosen not to print the MM_DATA and the'
   WRITE(LUPRI,*)'MM_CNTS files which means that you take these files from the'
   WRITE(LUPRI,*)'traditionel dalton and not the new Integral-interface'
ELSE
   !BUILD MM_CNT0
   LUINTM2 = -1
   CALL OPENMMFILE(LUINTM2,'MM_CNT0',LUPRI)
   WRITE (LUINTM2) LHSDENSFIT, RHSDENSFIT
   call lsclose(LUINTM2,'KEEP')

   IF(.NOT.SETTING%SCHEME%CREATED_MMFILES)THEN
      CALL LSTIMER('START ',TS,TE,LUPRI)
      ! get info if buffer should be used
      INT_OUTPUT%USEBUFMM = SETTING%SCHEME%USEBUFMM
      IF(INT_OUTPUT%USEBUFMM) THEN
         ! IDER SHOULD BE SET TO 2 OR HIGHER FOR GRADIENT AND HIGHER DERIAVATIVES
         ! AS THESE ARE NOT YET IMPLEMENTED WE SET IDER = 1
         IDER = 1
         ! initialize buffer, allocate memory
         call LS_INITMMBUF(INT_OUTPUT,IDER)
      END IF
      !
      SETTING%SCHEME%LU_LUINTM = -1
      CALL OPENMMFILE(SETTING%SCHEME%LU_LUINTM,'MM_DATA',LUPRI)
      IF (INT_OUTPUT%USEBUFMM)THEN
         SETTING%SCHEME%LU_LUINTR = -1
         CALL OPENMMFILE(SETTING%SCHEME%LU_LUINTR,'MM_DATR',LUPRI)
      ENDIF
      MMunique_ID1 = 0
      MMunique_ID2 = 0
      start1 = 0
      start2 = 0
      IF(FOURCENTER)THEN
         CALL MM_calculation(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
         IF(INT_OUTPUT%USEBUFMM) THEN
            ! we empty the buffer and write the stop signal
            IF (INT_OUTPUT%IBUFI .GT. 1) THEN
               CALL LS_EMPTYIBUF(INT_OUTPUT,INT_OUTPUT%IBUF,SETTING%SCHEME%LU_LUINTM)
               CALL LS_EMPTYRBUF(INT_OUTPUT,INT_OUTPUT%RBUF,SETTING%SCHEME%LU_LUINTR)
            END IF
            WRITE(SETTING%SCHEME%LU_LUINTM) -1
            WRITE(SETTING%SCHEME%LU_LUINTR) -1
            ! we buffer also the saving of the nuclear information
            INT_OUTPUT%IBUFN = 1
            DO I = 1, SETTING%MOLECULE(1)%p%nATOMS
               IF (INT_OUTPUT%IBUFN .GE. INT_OUTPUT%MMBUFLEN-1) THEN
                  CALL LS_EMPTYNUCBUF(INT_OUTPUT,INT_OUTPUT%NBUF,SETTING%SCHEME%LU_LUINTR)
               END IF
               CALL LS_FILLNUCBUF(INT_OUTPUT,SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE,  &
                    &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1),SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2),& 
                    &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3))
            END DO
            CALL LS_EMPTYNUCBUF(INT_OUTPUT,INT_OUTPUT%NBUF,SETTING%SCHEME%LU_LUINTR)
            WRITE(SETTING%SCHEME%LU_LUINTR) -2
         ELSE
            WRITE(SETTING%SCHEME%LU_LUINTM) -1,0,0,0,0,0,0,0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
            DO I = 1,SETTING%MOLECULE(1)%p%nATOMS
               WRITE(SETTING%SCHEME%LU_LUINTM) SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE,&
                    &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(:)
            ENDDO
         END IF
         CALL LSTIMER('MULMOM',TS,TE,LUPRI)
      ELSEIF(RHSDENSFIT)THEN
         CALL MM_calculation(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
         CALL MM_calculation(AO3,AO4,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
         IF(INT_OUTPUT%USEBUFMM) THEN
            ! we empty the buffer
            IF (INT_OUTPUT%IBUFI .GT. 1) THEN
               CALL LS_EMPTYIBUF(INT_OUTPUT,INT_OUTPUT%IBUF,SETTING%SCHEME%LU_LUINTM)
               CALL LS_EMPTYRBUF(INT_OUTPUT,INT_OUTPUT%RBUF,SETTING%SCHEME%LU_LUINTR)
            END IF
            WRITE(SETTING%SCHEME%LU_LUINTM) -1
            WRITE(SETTING%SCHEME%LU_LUINTR) -1
            ! we buffer also the saving of the nuclear information
            INT_OUTPUT%IBUFN = 1
            DO I = 1, SETTING%MOLECULE(1)%p%nATOMS
               IF (INT_OUTPUT%IBUFN .GE. INT_OUTPUT%MMBUFLEN-1) THEN
                  CALL LS_EMPTYNUCBUF(INT_OUTPUT,INT_OUTPUT%NBUF,SETTING%SCHEME%LU_LUINTR)
               END IF
               CALL LS_FILLNUCBUF(INT_OUTPUT,SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE,  &
                         &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1),SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2),&
                         &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3))
            END DO
            CALL LS_EMPTYNUCBUF(INT_OUTPUT,INT_OUTPUT%NBUF,SETTING%SCHEME%LU_LUINTR)
            WRITE(SETTING%SCHEME%LU_LUINTR) -2
         ELSE
            WRITE(SETTING%SCHEME%LU_LUINTM) -1,0,0,0,0,0,0,0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
            DO I = 1,SETTING%MOLECULE(1)%p%nATOMS
               WRITE(SETTING%SCHEME%LU_LUINTM) SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE,&
                    &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(:)
            ENDDO
         END IF
         CALL LSTIMER('MULMOM-RHSDENS',TS,TE,LUPRI)
      ELSEIF(LHSDENSFIT)THEN
         CALL MM_calculation(AO3,AO4,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
         CALL MM_calculation(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
         IF(INT_OUTPUT%USEBUFMM) THEN
            ! we empty the buffer
            IF (INT_OUTPUT%IBUFI .GT. 1) THEN
               CALL LS_EMPTYIBUF(INT_OUTPUT,INT_OUTPUT%IBUF,SETTING%SCHEME%LU_LUINTM)
               CALL LS_EMPTYRBUF(INT_OUTPUT,INT_OUTPUT%RBUF,SETTING%SCHEME%LU_LUINTR)
            END IF
            WRITE(SETTING%SCHEME%LU_LUINTM) -1
            WRITE(SETTING%SCHEME%LU_LUINTR) -1
            ! we buffer also the saving of the nuclear information
            INT_OUTPUT%IBUFN = 1
            DO I = 1, SETTING%MOLECULE(1)%p%nATOMS
               IF (INT_OUTPUT%IBUFN .GE. INT_OUTPUT%MMBUFLEN-1) THEN
                  CALL LS_EMPTYNUCBUF(INT_OUTPUT,INT_OUTPUT%NBUF,SETTING%SCHEME%LU_LUINTR)
               END IF
               CALL LS_FILLNUCBUF(INT_OUTPUT,SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE,  &
                         &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1),SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2),& 
                         &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3))
            END DO
            CALL LS_EMPTYNUCBUF(INT_OUTPUT,INT_OUTPUT%NBUF,SETTING%SCHEME%LU_LUINTR)
            WRITE(SETTING%SCHEME%LU_LUINTR) -2
         ELSE
            WRITE(SETTING%SCHEME%LU_LUINTM) -1,0,0,0,0,0,0,0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
            DO I = 1,SETTING%MOLECULE(1)%p%nATOMS
               WRITE(SETTING%SCHEME%LU_LUINTM) SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE,&
                    &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(:)
            ENDDO
         END IF
         CALL LSTIMER('MULMOM-LHSDENS',TS,TE,LUPRI)
      ENDIF
      call lsclose(SETTING%SCHEME%LU_LUINTM,'KEEP')
      IF(INT_OUTPUT%USEBUFMM)  THEN
         call lsclose(SETTING%SCHEME%LU_LUINTR,'KEEP')
         call LS_FREEMMBUF(INT_OUTPUT)
      ENDIF

      !BUILD MM_CNTS
      LUINTM2 = -1
      CALL OPENMMFILE(LUINTM2,'MM_CNTS',LUPRI)
      WRITE (LUINTM2) SETTING%SCHEME%MM_LMAX, NBAST, SETTING%SCHEME%MMunique_ID1,&
           &SETTING%MOLECULE(1)%p%NATOMS, IDUMMY,nbastaux,LHSDENSFIT,RHSDENSFIT
      call lsclose(LUINTM2,'KEEP')
      IF(SETTING%SCHEME%CFG_LSDALTON .AND. COULOMB)THEN
         SETTING%SCHEME%CREATED_MMFILES = .TRUE.
      ELSE
         !WE HAVE TO OVERWRITE EVERY TIME BECAUSE "FCK3" would overwrite the file
      ENDIF
   ELSE
      !BUILD MM_CNTS
      LUINTM2 = -1
      CALL OPENMMFILE(LUINTM2,'MM_CNTS',LUPRI)
      WRITE (LUINTM2) SETTING%SCHEME%MM_LMAX, NBAST, SETTING%SCHEME%MMunique_ID1,&
           &SETTING%MOLECULE(1)%p%NATOMS, IDUMMY,nbastaux,LHSDENSFIT,RHSDENSFIT
      call lsclose(LUINTM2,'KEEP')
   ENDIF
   
   IF(SETTING%LHSdmat)THEN
    nrowLHS = SETTING%DmatLHS(1)%p%nrow
    ncolLHS = SETTING%DmatLHS(1)%p%ncol
    call mem_alloc(Dfull,nrowLHS,nrowLHS)
    IF(matrix_type .EQ. mtype_unres_dense)THEN
     CALL DCOPY(nrowLHS*ncolLHS,SETTING%DmatLHS(1)%p%elms,1,Dfull,1)
     CALL DAXPY(nrowLHS*ncolLHS,1.d0,SETTING%DmatLHS(1)%p%elmsb,1,Dfull,1)
     CALL DSCAL(nrowLHS*ncolLHS,0.5d0,Dfull,1)
    ELSE
     call mat_to_full(SETTING%DmatLHS(1)%p,1D0,Dfull(:,:))
    ENDIF
   ELSEIF(SETTING%LHSdfull)THEN
    nrowLHS = ndim1
    ncolLHS = ndim2
    call mem_alloc(Dfull,nrowLHS,nrowLHS)
    IF(matrix_type .EQ. mtype_unres_dense)THEN
     call lsquit('this use of ls_multipole have not been tested')
     CALL DCOPY(nrowLHS*ncolLHS,SETTING%DfullLHS,1,Dfull,1)
!     CALL DAXPY(nrowLHS*ncolLHS,1.d0,SETTING%DfullLHS,1,Dfull,1)
     CALL DSCAL(nrowLHS*ncolLHS,0.5d0,Dfull,1)
    ELSE
       CALL DCOPY(nrowLHS*ncolLHS,SETTING%DfullLHS,1,Dfull,1)
    ENDIF
   ELSE
      call lsquit('LHS density matrix not attached to setting in ls_multipolemoment')
   ENDIF

   IF(SETTING%RHSdmat)THEN
    nrowRHS = SETTING%DmatRHS(1)%p%nrow
    ncolRHS = SETTING%DmatRHS(1)%p%ncol
    call mem_alloc(RHSDENS,nrowRHS,ncolRHS)
    IF(matrix_type .EQ. mtype_unres_dense)THEN
     CALL DCOPY(nrowRHS*ncolRHS,SETTING%DmatRHS(1)%p%elms,1,RHSDENS,1)
     CALL DAXPY(nrowRHS*ncolRHS,1.d0,SETTING%DmatRHS(1)%p%elmsb,1,RHSDENS,1)
     CALL DSCAL(nrowRHS*ncolRHS,0.5d0,RHSDENS,1)
    ELSE
     call mat_to_full(SETTING%DmatRHS(1)%p,1D0,RHSDENS(:,:))
    ENDIF
   ELSEIF(SETTING%RHSdfull)THEN
    nrowRHS = ndim3
    ncolRHS = ndim4
    call mem_alloc(RHSDENS,nrowRHS,nrowRHS)
    IF(matrix_type .EQ. mtype_unres_dense)THEN
     call lsquit('this use of ls_multipole have not been tested')
     CALL DCOPY(nrowRHS*ncolRHS,SETTING%DfullRHS,1,RHSDENS,1)
!     CALL DAXPY(nrowRHS*ncolRHS,1.d0,SETTING%DfullRHS,1,Dfull,1)
     CALL DSCAL(nrowRHS*ncolRHS,0.5d0,RHSDENS,1)
    ELSE
     call dcopy(nrowRHS*ncolRHS,SETTING%DfullRHS,1,RHSDENS,1)
    ENDIF
   ELSE
      call lsquit('RHS density matrix not attached to setting in ls_multipolemoment')
   ENDIF
   
   LU_DENS = -1
   IF(FOURCENTER)THEN
      CALL OPENMMFILE(LU_DENS,'MM_DENS',LUPRI)
      N = (NBAST*(NBAST-1))/2+NBAST
      call mem_alloc(WORK,N)
      !ALLOCATE(WORK(N))
      T = 0
      DO I = 1, NBAST
         DO J = 1, I
            WORK(T+J) = TWO*DFULL(I,J)+TWO*DFULL(J,I)
         END DO
         T = T + I
         WORK(I*(I-1)/2+I)=TWO*DFULL(I,I)
      END DO
      REWIND(LU_DENS)
      WRITE(LU_DENS) NBAST
      CALL ls_write(LU_DENS,WORK,N)
      call mem_dealloc(WORK)
      !DEALLOCATE(WORK)
      call lsclose(LU_DENS,'KEEP')
   ELSEIF(RHSDENSFIT)THEN
      CALL OPENMMFILE(LU_DENS,'MM_DENR',LUPRI)
      N = (NBAST*(NBAST-1))/2+NBAST
      call mem_alloc(WORK,N)
      !ALLOCATE(WORK(N))
      T = 0
      DO I = 1, NBAST
         DO J = 1, I
            WORK(T+J) = TWO*DFULL(I,J)+TWO*DFULL(J,I)
         END DO
         T = T + I
         WORK(I*(I-1)/2+I)=TWO*DFULL(I,I)
      END DO
      REWIND(LU_DENS)
      WRITE(LU_DENS) NBAST
      CALL ls_write(LU_DENS,WORK,N)
      call mem_dealloc(WORK)
      !DEALLOCATE(WORK)
      call lsclose(LU_DENS,'KEEP')
      CALL OPENMMFILE(LU_DENS,'MM_FRDN',LUPRI)
      REWIND(LU_DENS)
      WRITE(LU_DENS) NBASTAUX
      CALL ls_write(LU_DENS,RHSDENS,NBASTAUX)
      call lsclose(LU_DENS,'KEEP')
   ELSEIF(LHSDENSFIT)THEN
      CALL OPENMMFILE(LU_DENS,'MM_DENL',LUPRI)
      N = (NBAST*(NBAST-1))/2+NBAST
      call mem_alloc(WORK,N)
      !ALLOCATE(WORK(N))
      T = 0
      DO I = 1, NBAST
         DO J = 1, I
            WORK(T+J) = TWO*DFULL(I,J)+TWO*DFULL(J,I)
         END DO
         T = T + I
         WORK(I*(I-1)/2+I)=TWO*DFULL(I,I)
      END DO
      REWIND(LU_DENS)
      WRITE(LU_DENS) NBAST
      CALL ls_write(LU_DENS,WORK,N)
      call mem_dealloc(WORK)
      !DEALLOCATE(WORK)
      call lsclose(LU_DENS,'KEEP')
      CALL OPENMMFILE(LU_DENS,'MM_FLDN',LUPRI)
      call mem_alloc(WORK,NBASTAUX)
      !ALLOCATE(WORK(NBASTAUX))
      DO J = 1, NBASTAUX
         WORK(J) = 0.D0
      ENDDO
      REWIND(LU_DENS)
      WRITE(LU_DENS) NBASTAUX
      CALL ls_write(LU_DENS,RHSDENS,NBASTAUX)
      call mem_dealloc(WORK)
      !DEALLOCATE(WORK)
      call lsclose(LU_DENS,'KEEP')
   ENDIF
   call mem_dealloc(Dfull)
   call mem_dealloc(RHSDENS)
ENDIF

END SUBROUTINE ls_multipolemoment
      

!> \brief opens the file to write the multipoles, which is then read by the FMM driver.
!> \author T. Kjaergaard
!> \date 2010
!> \param LUINTM logical unit number of the MM file
!> \param FILENAME the name of the MM file
!> \param LUPRI logical unit number of the Default output file
SUBROUTINE OPENMMFILE(LUINTM,FILENAME,LUPRI)
Character*(*)   :: FILENAME
INTEGER         :: LUINTM!logic unit number
INTEGER         :: LUPRI,IDUMMY
LOGICAL         :: fileexist  

INQUIRE(file=FILENAME,EXIST=fileexist) 
IF(fileexist)THEN
   CALL lsOPEN(LUINTM,FILENAME,'UNKNOWN','UNFORMATTED')
   call lsclose(LUINTM,'DELETE')
   INQUIRE(file=FILENAME,EXIST=fileexist) 
   IF(fileexist)THEN
      WRITE(LUPRI,*)'ERROR',FILENAME ,' not deleted in ls_getmultipolemom'
      WRITE(*,*)'ERROR',FILENAME ,' not deleted in ls_getmultipolemom'
      CALL LSQUIT('ERROR',FILENAME ,' not deleted in ls_getmultipolemom',lupri)
   ELSE
!     WRITE(LUPRI,*)FILENAME ,' exist so we overwrite'
      CALL lsOPEN(LUINTM,FILENAME,'UNKNOWN','UNFORMATTED')
   ENDIF
ELSE
   WRITE(LUPRI,*)FILENAME ,' do not exist so we make a one'
   CALL lsOPEN(LUINTM,FILENAME,'UNKNOWN','UNFORMATTED')
ENDIF

END SUBROUTINE OPENMMFILE

!> \brief wrapper for the calculation of the multipoles moments used by FMM
!> \author T. Kjaergaard
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param start1 the start index for the 1 dimension
!> \param start2 the start index for the 2 dimension
!> \param MMunique_ID1 a unique identifier used in the file,required by FMM code
!> \param MMunique_ID2 a unique identifier used in the file,required by FMM code
!> \param integral_output the output specifications of the integral eval
SUBROUTINE MM_calculation(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR
TYPE(LSSETTING)       :: SETTING
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
Character*(*)         :: AO1,AO2,intType
Integer               :: start1,start2,MMunique_ID1,MMunique_ID2
Integer               :: s1,s2,I,J,inode
logical               :: samefragment,Ldummy

IF(SETTING%SCHEME%FRAGMENT) THEN
   CALL buildFragmentInfoAndBlocks(SETTING,AO1,AO2,'Empty','Empty',LUPRI,LUERR)
   SETTING%FRAGMENTS%iRHSblock = 1
   DO I=1,SETTING%FRAGMENTS%LHSblock%numBlocks
      iNode = SETTING%FRAGMENTS%LHSblock%blocks(I)%node
      CALL SetDaltonFragments(SETTING,I,1,sameFragment,Ldummy,lupri)
      SETTING%FRAGMENTS%iLHSblock = I
      s1 = SETTING%FRAGMENTS%LHSblock%blocks(I)%startOrb1-1
      s2 = SETTING%FRAGMENTS%LHSblock%blocks(I)%startOrb2-1
      if((s2 .GE. s1).OR.(AO1.EQ.'Empty').OR.(AO2.EQ.'Empty'))THEN
         Call MM_kernel(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1+s1,start2+s2,&
     &                  MMunique_ID1,MMunique_ID2,INT_OUTPUT)
      endif
      CALL FreeDaltonFragments(SETTING)
   ENDDO
   CALL freeFragmentInfoAndBlocks(SETTING)
ELSE
   Call MM_kernel(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,&
     &            MMunique_ID1,MMunique_ID2,INT_OUTPUT)
ENDIF
SETTING%SCHEME%MMunique_ID1 = MMunique_ID1

end subroutine MM_calculation

!> \brief calculate the multipoles and write them to a file, which is then read by the FMM driver.
!> \author T. Kjaergaard
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param sA the start index for the 1 dimension
!> \param sB the start index for the 2 dimension
!> \param MMunique_ID1 a unique identifier used in the file,required by FMM code
!> \param MMunique_ID2 a unique identifier used in the file,required by FMM code
!> \param integral_output the output specifications of the integral eval
subroutine MM_kernel(AO1,AO2,intType,SETTING,LUPRI,LUERR,sA,sB,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR
TYPE(LSSETTING)       :: SETTING
Character*(*)         :: AO1,AO2,intType
Real(realk),pointer   :: integrals(:,:,:,:,:)
INTEGER               :: nderiv,nmat
TYPE(AOITEM),target   :: AObuild(4)
Integer               :: nAObuilds
Integer               :: MMunique_ID1,MMunique_ID2,sA,sB
TYPE(INTEGRALINPUT)   :: INT_INPUT
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT

CALL init_integral_input(INT_INPUT,SETTING)
INT_INPUT%LU_MMDATA = SETTING%SCHEME%LU_LUINTM
IF (INT_OUTPUT%USEBUFMM) INT_INPUT%LU_MMDATR = SETTING%SCHEME%LU_LUINTR
INT_INPUT%operator = 'Mulmom'
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN        
CALL SetInputAO(INT_INPUT,AO1,AO2,'Empty','Empty',intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)
! We screen the multipole moments because if the 
! integral (ab|cd)D_{cd} is not calculated then there is no need to 
! build the multipole moments.
! very generally then we do not know if we have the expansion
! M_{ab}*M_{cd}*D_{cd}
! or 
! M_{ab}*M_{alpha}*C_{alpha}
! so when screening we set GAB_RHS = 1 and build M_{ab} if
! G_{ab} > modifiedThreshold = Threshold*1.D-2
! so we set Gab_{RHS} = 1 inside SET_OVERLAP inside MAININTDRIVER
IF(sA .EQ. sB)THEN
   INT_INPUT%sameLHSaos = (AO1.EQ.AO2)
ELSE
   INT_INPUT%sameLHSaos = .FALSE.
ENDIF

NULLIFY(integrals)
ALLOCATE(integrals(1,1,1,1,1))
Int_Output%ndim = 1
nderiv = SETTING%SCHEME%MM_LMAX
nMAT = (nderiv+1)*(nderiv+2)*(nderiv+3)/6
INT_INPUT%NDERIVQ = nmat
INT_INPUT%DERIVORDER = nderiv
!INT_INPUT%sameLHSaos = .TRUE.
!INT_INPUT%sameRHSaos = .TRUE.
INT_INPUT%do_mulmom = .TRUE.
Int_Output%ResultMat => Integrals
INT_INPUT%MMunique_ID1 = MMunique_ID1
INT_INPUT%MMunique_ID2 = MMunique_ID2
INT_INPUT%MMstartA = sA
INT_INPUT%MMstartB = sB
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,INT_OUTPUT)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
MMunique_ID1 = INT_INPUT%MMunique_ID1
MMunique_ID2 = INT_INPUT%MMunique_ID2
DEALLOCATE(integrals)
NULLIFY(integrals)
end subroutine MM_kernel

!> \brief determine the symmetries: Is both the LHS AOs the same? RHS? are both Overlap distributions the same
!> \author S. Reine
!> \date 2010
!> \param sameAOsLHS if the LHS AOS are the same
!> \param sameAOsRHS if the RHS AOS are the same
!> \param sameODs if both ODs are the same
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
SUBROUTINE GetSymmetries(sameAOsLHS,sameAOsRHS,sameODs,AO1,AO2,AO3,AO4)
implicit none
Logical            :: sameAOsLHS,sameAOsRHS,sameODs
Character*(*)      :: AO1,AO2,AO3,AO4

sameAOsLHS = (AO1.EQ.AO2).AND..NOT.(AO1.EQ.'Empty')
sameAOsRHS = (AO3.EQ.AO4).AND..NOT.(AO3.EQ.'Empty')
sameODs    = (AO1.EQ.AO3) .AND. (AO2.EQ.AO4)

END SUBROUTINE GetSymmetries

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_matarrayp(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
use matrix_util, only : mat_get_isym
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrixp),intent(IN),target :: Dmat(ndmat)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
integer                :: idmat
real(realk), parameter :: thresh = 1.d-15

call ls_attachDmatToSetting2matp(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)

END SUBROUTINE ls_attachDmatToSetting_matarrayp

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_matarray(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
use matrix_util, only : mat_get_isym
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrix),intent(IN),target  :: Dmat(ndmat)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
TYPE(matrixp)  :: Dmatp(ndmat)
integer                :: idmat
call ls_attachDmatToSetting2mat(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)

END SUBROUTINE ls_attachDmatToSetting_matarray

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_matsinglep(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
use matrix_util, only : mat_get_isym
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrixp),intent(IN),target :: Dmat
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
!TYPE(matrixp)  :: Dmatp(1)
integer                :: idmat
real(realk), parameter :: thresh = 1.d-15

call ls_attachDmatToSetting2matp((/Dmat/),ndmat,setting,side,AOindex1,AOindex2,lupri)

END SUBROUTINE ls_attachDmatToSetting_matsinglep

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_matsingle(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
use matrix_util, only : mat_get_isym
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrix),intent(IN),target  :: Dmat
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
!TYPE(matrixp)          :: Dmatp(1)
integer                :: idmat

call ls_attachDmatToSetting2matsingle(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)

END SUBROUTINE ls_attachDmatToSetting_matsingle

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_full(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
use matrix_util, only : mat_get_isym
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
real(realk),intent(IN)          :: Dmat(:,:,:)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
call ls_attachDmatToSetting2full(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)

END SUBROUTINE ls_attachDmatToSetting_full

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting2matp(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
use matrix_util, only : mat_get_isym
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrixp),intent(IN),target :: Dmat(ndmat)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
integer                :: idmat
real(realk), parameter :: thresh = 1.d-15

SELECT CASE(side)
CASE('LHS')
  IF (setting%LHSdmat) CALL LSQUIT('Error in ls_attachDmatToSetting. LHS',lupri)
  setting%nDmatLHS = ndmat
  setting%LHSdmat  = .TRUE.
  nullify(setting%DmatLHS)
  allocate(setting%DmatLHS(ndmat))
  do idmat=1,ndmat
     setting%DmatLHS(idmat)%p => Dmat(idmat)%p
  enddo
  call mem_alloc(setting%DsymLHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymLHS(idmat) = mat_get_isym(Dmat(idmat)%p,thresh)
  ENDDO
  setting%LHSdmatAOindex1 = AOindex1
  setting%LHSdmatAOindex2 = AOindex2
CASE('RHS')
  IF (setting%RHSdmat) CALL LSQUIT('Error in ls_attachDmatToSetting. RHS',lupri)
  setting%nDmatRHS = ndmat
  setting%RHSdmat  = .TRUE.
  nullify(setting%DmatRHS)
  allocate(setting%DmatRHS(ndmat))
  do idmat=1,ndmat
     setting%DmatRHS(idmat)%p => Dmat(idmat)%p
  enddo
  call mem_alloc(setting%DsymRHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymRHS(idmat) = mat_get_isym(Dmat(idmat)%p,thresh)
  ENDDO
  setting%RHSdmatAOindex1 = AOindex1
  setting%RHSdmatAOindex2 = AOindex2
CASE DEFAULT
  WRITE(LUPRI,'(1X,2A)') 'Error in ls_attachDmatToSetting. Side =',side
  CALL LSQUIT('Error in ls_attachDmatToSetting. Wrong side',lupri)
END SELECT

END SUBROUTINE ls_attachDmatToSetting2matp

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting2mat(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
use matrix_util, only : mat_get_isym
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrix),intent(IN),target  :: Dmat(ndmat)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
integer                :: idmat
real(realk), parameter :: thresh = 1.d-15

SELECT CASE(side)
CASE('LHS')
  IF (setting%LHSdmat) CALL LSQUIT('Error in ls_attachDmatToSetting. LHS',lupri)
  setting%nDmatLHS = ndmat
  setting%LHSdmat  = .TRUE.
  nullify(setting%DmatLHS)
  allocate(setting%DmatLHS(ndmat))
  do idmat=1,ndmat
     setting%DmatLHS(idmat)%p => Dmat(idmat)
  enddo
  call mem_alloc(setting%DsymLHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymLHS(idmat) = mat_get_isym(Dmat(idmat),thresh)
  ENDDO
  setting%LHSdmatAOindex1 = AOindex1
  setting%LHSdmatAOindex2 = AOindex2
CASE('RHS')
  IF (setting%RHSdmat) CALL LSQUIT('Error in ls_attachDmatToSetting. RHS',lupri)
  setting%nDmatRHS = ndmat
  setting%RHSdmat  = .TRUE.
  nullify(setting%DmatRHS)
  allocate(setting%DmatRHS(ndmat))
  do idmat=1,ndmat
     setting%DmatRHS(idmat)%p => Dmat(idmat)
  enddo
  call mem_alloc(setting%DsymRHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymRHS(idmat) = mat_get_isym(Dmat(idmat),thresh)
  ENDDO
  setting%RHSdmatAOindex1 = AOindex1
  setting%RHSdmatAOindex2 = AOindex2
CASE DEFAULT
  WRITE(LUPRI,'(1X,2A)') 'Error in ls_attachDmatToSetting. Side =',side
  CALL LSQUIT('Error in ls_attachDmatToSetting. Wrong side',lupri)
END SELECT

END SUBROUTINE ls_attachDmatToSetting2mat

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting2matsingle(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
use matrix_util, only : mat_get_isym
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrix),intent(IN),target  :: Dmat
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
integer                :: idmat
real(realk), parameter :: thresh = 1.d-15

SELECT CASE(side)
CASE('LHS')
  IF (setting%LHSdmat) CALL LSQUIT('Error in ls_attachDmatToSetting. LHS',lupri)
  setting%nDmatLHS = ndmat
  setting%LHSdmat  = .TRUE.
  nullify(setting%DmatLHS)
  allocate(setting%DmatLHS(ndmat))
  do idmat=1,ndmat
     setting%DmatLHS(idmat)%p => Dmat
  enddo
  call mem_alloc(setting%DsymLHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymLHS(idmat) = mat_get_isym(Dmat,thresh)
  ENDDO
  setting%LHSdmatAOindex1 = AOindex1
  setting%LHSdmatAOindex2 = AOindex2
CASE('RHS')
  IF (setting%RHSdmat) CALL LSQUIT('Error in ls_attachDmatToSetting. RHS',lupri)
  setting%nDmatRHS = ndmat
  setting%RHSdmat  = .TRUE.
  nullify(setting%DmatRHS)
  allocate(setting%DmatRHS(ndmat))
  do idmat=1,ndmat
     setting%DmatRHS(idmat)%p => Dmat
  enddo
  call mem_alloc(setting%DsymRHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymRHS(idmat) = mat_get_isym(Dmat,thresh)
  ENDDO
  setting%RHSdmatAOindex1 = AOindex1
  setting%RHSdmatAOindex2 = AOindex2
CASE DEFAULT
  WRITE(LUPRI,'(1X,2A)') 'Error in ls_attachDmatToSetting. Side =',side
  CALL LSQUIT('Error in ls_attachDmatToSetting. Wrong side',lupri)
END SELECT

END SUBROUTINE ls_attachDmatToSetting2matsingle

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting2full(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
use matrix_util, only : mat_get_isym
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
real(realk),intent(IN),target   :: Dmat(:,:,:)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
integer                :: idmat
real(realk), parameter :: thresh = 1.d-15

SELECT CASE(side)
CASE('LHS')
  IF (setting%LHSdfull) CALL LSQUIT('Error in ls_attachDmatToSetting. LHS',lupri)
  setting%nDmatLHS = ndmat
  setting%LHSdfull  = .TRUE.
  setting%DfullLHS => Dmat
  setting%LHSdmatAOindex1 = AOindex1
  setting%LHSdmatAOindex2 = AOindex2
  call mem_alloc(setting%DsymLHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymLHS(idmat) = 1 !sym !3 !no symmetry   VERY TEMP
  ENDDO
CASE('RHS')
  IF (setting%RHSdfull) CALL LSQUIT('Error in ls_attachDmatToSetting. RHS',lupri)
  setting%nDmatRHS = ndmat
  setting%RHSdfull  = .TRUE.
  setting%DfullRHS => Dmat
  setting%RHSdmatAOindex1 = AOindex1
  setting%RHSdmatAOindex2 = AOindex2
  call mem_alloc(setting%DsymRHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymRHS(idmat) = 1! sym 3 !no symmetry tmp 
  ENDDO
CASE DEFAULT
  WRITE(LUPRI,'(1X,2A)') 'Error in ls_attachDmatToSetting. Side =',side
  CALL LSQUIT('Error in ls_attachDmatToSetting. Wrong side',lupri)
END SELECT

END SUBROUTINE ls_attachDmatToSetting2full

!> \brief Deassosiate density matrices from setting (if they have been associated)
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param setting Integral settings
SUBROUTINE ls_freeDmatFromSetting(setting)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
!
IF (setting%LHSdmat) THEN
  setting%LHSdmat = .FALSE.
  setting%nDmatLHS = 1
  IF (.NOT.ASSOCIATED(setting%DsymLHS)) THEN
     CALL LSQUIT('Error in call to ls_freeDmatFromSetting. DsymLHS',-1)
  ELSE
     call mem_dealloc(setting%DsymLHS)
  ENDIF
  if(associated(setting%DmatLHS))then
     deallocate(setting%DmatLHS)
     nullify(setting%DmatLHS)
  endif
ENDIF
!
IF (setting%RHSdmat) THEN
  setting%RHSdmat = .FALSE.
  setting%nDmatRHS = 1
  IF (.NOT.ASSOCIATED(setting%DsymRHS)) THEN
     CALL LSQUIT('Error in call to ls_freeDmatFromSetting. DsymRHS',-1)
  ELSE
     call mem_dealloc(setting%DsymRHS)
  ENDIF
  if(associated(setting%DmatRHS))then
     deallocate(setting%DmatRHS)
     nullify(setting%DmatRHS)
  endif
ENDIF
!

if(setting%LHSdalloc)then
   call mem_dealloc(Setting%DfullLHS)
   setting%LHSdalloc=.FALSE.
endif
IF (setting%LHSdfull) THEN
   setting%LHSdfull = .FALSE.
   setting%nDmatLHS = 1
   nullify(setting%DfullLHS)
ENDIF
IF (ASSOCIATED(setting%DsymLHS)) call mem_dealloc(setting%DsymLHS)
!
if(setting%RHSdalloc)then
   call mem_dealloc(Setting%DfullRHS)
   setting%RHSdalloc=.FALSE.
endif
IF (setting%RHSdfull) THEN
   setting%RHSdfull = .FALSE.
   setting%nDmatRHS = 1
   nullify(setting%DfullRHS)
ENDIF
IF (ASSOCIATED(setting%DsymRHS)) call mem_dealloc(setting%DsymRHS)

END SUBROUTINE ls_freeDmatFromSetting

!> \brief Determines whether the matrices are identical (one-by-one)
!> \author S. Reine
!> \date 2010-04-19
!> \param D First set of matrices
!> \param P Second set of matrices
!> \param nd The number of D matrices
!> \param np The number of P matrices
!> \param ls_same_mats True if D and P matrices are identical
FUNCTION ls_same_mats(D,P,nd,np)
use matrix_util, only : mat_same
implicit none
Integer,intent(IN)       :: nd,np
TYPE(matrixP),intent(IN) :: D(nd),P(np)
Logical                  :: ls_same_mats
!
logical :: same
integer :: id
real(realk),parameter :: thresh = 1.d-15

IF (nd.EQ.np) THEN
  same = .TRUE.
  DO id = 1,nd
    IF (.NOT.mat_same(D(id)%p,P(id)%p,thresh)) same = .FALSE.
  ENDDO
ELSE
  same = .FALSE.
ENDIF

ls_same_mats = same

END FUNCTION ls_same_mats

!> \brief Determines the symmetry of a matrix
!> \author S. Reine
!> \date 2010-07-01
!> \param D Matrix
!> \param ls_mat_sym 1 if symmetric, 2 if anti-symmetric and 3 if no symmetry, 4 zero matrix
FUNCTION ls_mat_sym(D)
use matrix_util, only : mat_get_isym
implicit none
TYPE(matrix),intent(IN) :: D
Integer                 :: ls_mat_sym
!
real(realk),parameter :: thresh = 1.d-15

ls_mat_sym = mat_get_isym(D,thresh)

END FUNCTION ls_mat_sym

!> \brief Allocate and extract a (full) sub-block of the density matrices (when applicable)
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param setting Integral settings
!> \param nbast1 The number of (consecutive) elements to be extracted for dimension 1
!> \param nbast2 The number of (consecutive) elements to be extracted for dimension 2
!> \param nbast3 The number of (consecutive) elements to be extracted for dimension 3
!> \param nbast4 The number of (consecutive) elements to be extracted for dimension 4
!> \param start1 The index of the first element to be extracted for dimension 1
!> \param start2 The index of the first element to be extracted for dimension 2
!> \param start3 The index of the first element to be extracted for dimension 3
!> \param start4 The index of the first element to be extracted for dimension 4
!> \param permuteLHS True if one should add the transpose of the LHS density-matrix
!> \param permuteRHS True if one should add the transpose of the RHS density-matrix
SUBROUTINE ls_setDblock(setting,nbast1,nbast2,nbast3,nbast4,start1,start2,start3,start4,&
     &                 permuteLHS,permuteRHS)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
Integer,intent(IN)            :: nbast1,nbast2,nbast3,nbast4
Integer,intent(IN)            :: start1,start2,start3,start4
Logical,intent(IN)            :: permuteLHS,permuteRHS

IF (setting%LHSdmat) THEN
   call mem_alloc(setting%DfullLHS,nbast1,nbast2,setting%nDmatLHS)
   CALL ls_mat_retrive_block(setting%DmatLHS,setting%DfullLHS,setting%ndmatLHS,&
      &                      nbast1,nbast2,start1,start2,permuteLHS)
   setting%LHSdalloc = .TRUE.
   setting%LHSdfull = .TRUE.
ENDIF
!
IF (setting%RHSdmat) THEN
   call mem_alloc(setting%DfullRHS,nbast3,nbast4,setting%nDmatRHS)
   CALL ls_mat_retrive_block(setting%DmatRHS,setting%DfullRHS,setting%ndmatRHS,&
      &                      nbast3,nbast4,start3,start4,permuteRHS)
   setting%RHSdalloc = .TRUE.
   setting%RHSdfull = .TRUE.
ENDIF
!
END SUBROUTINE ls_setDblock

!> \brief Allocate the full density matrices and convert from matrix type (when applicable)
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param setting Integral settings
!> Added unrestricted case 2010-04-06 S.Reine
SUBROUTINE ls_setDfull(setting)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
!
Integer :: n1,n2,n3,n4,I

n1 = 1
n2 = 1
n3 = 1
n4 = 1
IF (setting%LHSdmat) THEN
  n1 = setting%DmatLHS(1)%p%nrow
  n2 = setting%DmatLHS(1)%p%ncol
ENDIF
IF (setting%RHSdmat) THEN
  n3 = setting%DmatRHS(1)%p%nrow
  n4 = setting%DmatRHS(1)%p%ncol
ENDIF
!
!
IF(matrix_type .NE. mtype_unres_dense)THEN
  CALL ls_setDblock(setting,n1,n2,n3,n4,1,1,1,1,.FALSE.,.FALSE.)
ELSE ! Unrestricted
  IF (setting%LHSdmat) THEN
    CALL LSQUIT('Error in ls_setDfull. Unrestricted LHSdmat not yet implemented',-1)
  ENDIF
  IF (setting%RHSdmat) THEN
    call mem_alloc(setting%DfullRHS,n3,n4,setting%nDmatRHS)
    setting%RHSdalloc = .TRUE.
    setting%RHSdfull = .TRUE.
    DO I =1,setting%ndmatRHS
      CALL DCOPY(n3*n4,setting%DmatRHS(I)%p%elms,1,setting%DfullRHS(:,:,I),1)
      CALL DAXPY(n3*n4,1.d0,setting%DmatRHS(I)%p%elmsb,1,setting%DfullRHS(:,:,I),1)
      CALL DSCAL(n3*n4,0.5d0,setting%DfullRHS(:,:,I),1)
    ENDDO
  ENDIF
ENDIF
!
END SUBROUTINE ls_setDfull

!> \brief Frees the density matrix subblock
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param setting Integral settings
SUBROUTINE ls_freeDblock(setting)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
!
IF (setting%LHSdalloc) THEN
   CALL mem_dealloc(setting%DfullLHS)
   setting%LHSdalloc=.FALSE.
ENDIF

IF (setting%RHSdalloc) THEN
   CALL mem_dealloc(setting%DfullRHS)
   setting%RHSdalloc=.FALSE.   
ENDIF
!
END SUBROUTINE ls_freeDblock

!> \brief Frees the full density matrix
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param setting Integral settings
SUBROUTINE ls_freeDfull(setting)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
CALL ls_freeDblock(setting)
END SUBROUTINE ls_freeDfull

!> \brief build an lstensor structure form a full 5 dim array
!> \author T. Kjaeergaard
!> \date 2010
!> \param tensor the output lstensor
!> \param SETTING Integral evalualtion settings
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param ndmat the size of the fifth dimension - number of matrices
!> \param intType the label for primitive or contracted calc
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param fullmat the full 5 dimensional array
SUBROUTINE ls_Build_lstensor_from_full_5dim(tensor,setting,AO1,AO2,AO3,AO4,ndmat,intType,LUPRI,LUERR,fullmat)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
type(lstensor)       :: tensor
Character*(*)        :: AO1,AO2,AO3,AO4,intType
real(realk)          :: fullmat(:,:,:,:,:)
Integer              :: LUPRI,LUERR,ndmat,dim
!
Character(80)        :: AOstring(4)
TYPE(AOITEM)         :: AO(4),lAO(4)
Integer              :: indAO,ndim(4),ndim2(4),IAO
logical              :: useAO(4)

AOstring(1)=AO1;AOstring(2)=AO2;AOstring(3)=AO3;AOstring(4)=AO4;
indAO = 0
DO IAO=1,4
   CALL SetAObatch(AO(iAO),setting%batchindex(iAO),ndim(iAO),AOstring(iAO),&
        & intType,'Default',SETTING%Scheme,setting%FRAGMENT(iAO)%p,&
        & Setting%BASIS(iAO)%p,LUPRI,LUERR)
ENDDO

useAO=.TRUE.
do IAO=1,4
   IF(setting%output%docontrule)THEN
      IF(setting%output%contrule(1)%OutIntCont(IAO).NE.0)THEN
         ndim2(IAO) = ndim(setting%output%contrule(1)%OutIntCont(IAO))
         lAO(IAO) = AO(setting%output%contrule(1)%OutIntCont(IAO))
      ELSE
         ndim2(IAO) = 1
         useAO(IAO)=.FALSE.
      ENDIF      
   ELSE
      ndim2(IAO) = ndim(IAO)      
      lAO(IAO) = AO(IAO)      
   ENDIF
enddo
call Build_lstensor_from_full_5dim(tensor,fullmat,lAO(1),lAO(2),lAO(3),lAO(4),ndim2(1),ndim2(2),&
     &ndim2(3),ndim2(4),ndmat,useAO(1),useAO(2),useAO(3),useAO(4),lupri)
DO IAO=1,4
   CALL free_aoitem(lupri,AO(iAO))
ENDDO

END SUBROUTINE ls_Build_lstensor_from_full_5dim

!> \brief build an lstensor structure form a full 2 dim array
!> \author T. Kjaeergaard
!> \date 2010
!> \param tensor the output lstensor
!> \param SETTING Integral evalualtion settings
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param intType the label for primitive or contracted calc
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param fullmat the full 2 dimensional array
SUBROUTINE ls_Build_lstensor_from_full_2dim(tensor,setting,AO1,AO2,intType,LUPRI,LUERR,fullmat)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
type(lstensor)       :: tensor
Character*(*)        :: AO1,AO2,intType
real(realk)          :: fullmat(:,:)
Integer              :: LUPRI,LUERR,dim
!
Character(80)        :: AOstring(2)
TYPE(AOITEM)         :: AO(2),lAO(2)
Integer              :: indAO,ndim(2),ndim2(2),IAO
logical              :: useAO(2)

AOstring(1)=AO1;AOstring(2)=AO2
indAO = 0
DO IAO=1,2
   CALL SetAObatch(AO(iAO),setting%batchindex(iAO),ndim(iAO),AOstring(iAO),&
        & intType,'Default',SETTING%Scheme,setting%FRAGMENT(iAO)%p,&
        & Setting%BASIS(iAO)%p,LUPRI,LUERR)
ENDDO
useAO=.TRUE.
do IAO=1,2
   ndim2(IAO) = ndim(IAO)      
   lAO(IAO) = AO(IAO)      
enddo
call Build_lstensor_from_full_2dim(tensor,fullmat,lAO(1),lAO(2),ndim2(1),ndim2(2),&
     & useAO(1),useAO(2),lupri)
DO IAO=1,2
   CALL free_aoitem(lupri,AO(iAO))
ENDDO

END SUBROUTINE ls_Build_lstensor_from_full_2dim

!> \brief init a lstensor which describes a 5 dimensional array
!> \author T. Kjaeergaard
!> \date 2010
!> \param SETTING Integral evalualtion settings
!> \param AO1 the type of the basis ('Regular','DF-Aux','Empty') on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux','Empty') on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux','Empty') on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux','Empty') on center 4
!> \param ndmat the size of the fifth dimension - number of matrices
!> \param intType the label for primitive or contracted calc
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_init_lstensor_5dim(setting,AO1,AO2,AO3,AO4,ndmat,intType,LUPRI,LUERR)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
Character*(*)        :: AO1,AO2,AO3,AO4,intType
Integer              :: LUPRI,LUERR,ndmat
!
Character(80)        :: AOstring(4)
TYPE(AOITEM)         :: AO(4),lAO(4)
Integer              :: indAO,ndim(4),ndim2(4),IAO
logical              :: useAO(4)

AOstring(1)=AO1;AOstring(2)=AO2;AOstring(3)=AO3;AOstring(4)=AO4;
indAO = 0
DO IAO=1,4
   CALL SetAObatch(AO(iAO),setting%batchindex(iAO),ndim(iAO),AOstring(iAO),&
        & intType,'Default',SETTING%Scheme,setting%FRAGMENT(iAO)%p,&
        & Setting%BASIS(iAO)%p,LUPRI,LUERR)
ENDDO

useAO=.TRUE.
do IAO=1,4
   IF(setting%output%docontrule)THEN
      IF(setting%output%contrule(1)%OutIntCont(IAO).NE.0)THEN
         ndim2(IAO) = ndim(setting%output%contrule(1)%OutIntCont(IAO))
         lAO(IAO) = AO(setting%output%contrule(1)%OutIntCont(IAO))
      ELSE
         ndim2(IAO) = 1
         useAO(IAO)=.FALSE.
      ENDIF      
   ELSE
      ndim2(IAO) = ndim(IAO)      
      lAO(IAO) = AO(IAO)      
   ENDIF
enddo

call init_lstensor_5dim(setting%Output%ResultMat2,lAO(1),lAO(2),lAO(3),lAO(4),ndim2(1),ndim2(2),&
     &ndim2(3),ndim2(4),ndmat,useAO(1),useAO(2),useAO(3),useAO(4),.FALSE.,lupri)

DO IAO=1,4
   CALL free_aoitem(lupri,AO(iAO))
ENDDO

END SUBROUTINE ls_init_lstensor_5dim

!> \brief set the densities in the integralinput structure correctly
!> \author S. Reine and T. Kjaeergaard
!> \date 2010
!> \param intinput the integral input to be set 
!> \param SETTING Integral evalualtion settings containing the densities
SUBROUTINE ls_setIntegralInputDensities(intinput,setting)
implicit none
TYPE(INTEGRALINPUT) :: INTINPUT
TYPE(LSSETTING)   :: SETTING
!
INTEGER :: LHS1,LHS2,RHS1,RHS2,nbastLHS1,nbastLHS2,nbastRHS1,nbastRHS2
INTEGER :: nDmatLHS,nDmatRHS,lupri,idmat
LOGICAL :: useAO1,useAO2

lupri = 6

IF ((INTINPUT%operator.EQ.'Kinetic') .AND.setting%RHSdfull)&
     &CALL LSQUIT('Error in ls_setintegralinputdensities. Kinetic and RHS density',-1)

IF (setting%LHSdfull) THEN
   LHS1 = setting%LHSdmatAOindex1
   LHS2 = setting%LHSdmatAOindex2
   nbastLHS1 = INTINPUT%AOdim(LHS1)
   nbastLHS2 = INTINPUT%AOdim(LHS2)
   nDmatLHS = setting%nDmatLHS
   CALL attachDmatToInput(INTINPUT,Setting%DfullLHS,nbastLHS1,nbastLHS2,nDmatLHS,'LHS')
   useAO1 = .NOT.INTINPUT%AO(LHS1)%p%empty
   useAO2 = .NOT.INTINPUT%AO(LHS2)%p%empty
   INTINPUT%uselst_DLHS = .TRUE.
   CALL Build_lstensor_from_full_3dim(INTINPUT%lst_DLHS,Setting%DfullLHS,INTINPUT%AO(LHS1)%p,&
        &INTINPUT%AO(LHS2)%p,nbastLHS1,nbastLHS2,nDmatLHS,useAO1,useAO2,lupri)
ELSEIF(setting%LHSdmat)THEN
   LHS1 = setting%LHSdmatAOindex1
   LHS2 = setting%LHSdmatAOindex2
   nbastLHS1 = INTINPUT%AOdim(LHS1)
   nbastLHS2 = INTINPUT%AOdim(LHS2)
   nDmatLHS = setting%nDmatLHS
   useAO1 = .NOT.INTINPUT%AO(LHS1)%p%empty
   useAO2 = .NOT.INTINPUT%AO(LHS2)%p%empty
   IF((INTINPUT%DO_EXCHANGE .AND. INTINPUT%DO_DALINK).AND. matrix_type .EQ. mtype_unres_dense)THEN
      call LSquit(' excahnge error unres in ls_atta',-1)
   ELSE
      call mem_alloc(setting%DfullLHS,nbastLHS1,nbastLHS2,setting%nDmatLHS)
      setting%LHSdalloc = .TRUE.
      CALL ls_mat_retrive_block(setting%DmatLHS,setting%DfullLHS,ndmatLHS,nbastLHS1,nbastLHS2,1,1,.FALSE.)
      CALL attachDmatToInput(INTINPUT,Setting%DfullLHS,nbastLHS1,nbastLHS2,nDmatLHS,'LHS')
      INTINPUT%uselst_DLHS = .TRUE.
      CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dLHS,Setting%DfullLHS,INTINPUT%AO(LHS1)%p,&
           &INTINPUT%AO(LHS2)%p,nbastLHS1,nbastLHS2,nDmatLHS,useAO1,useAO2,lupri)
   ENDIF
ENDIF
   
IF (setting%RHSdfull)THEN
   RHS1 = setting%RHSdmatAOindex1
   RHS2 = setting%RHSdmatAOindex2
   nbastRHS1 = INTINPUT%AOdim(RHS1)
   nbastRHS2 = INTINPUT%AOdim(RHS2)
   nDmatRHS = setting%nDmatRHS
   CALL attachDmatToInput(INTINPUT,Setting%DfullRHS,nbastRHS1,nbastRHS2,setting%nDmatRHS,'RHS')
   useAO1 = .NOT.INTINPUT%AO(RHS1)%p%empty
   useAO2 = .NOT.INTINPUT%AO(RHS2)%p%empty
   INTINPUT%uselst_DRHS = .TRUE.
   CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dRHS,Setting%DfullRHS,INTINPUT%AO(RHS1)%p,&
        &INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,lupri)
ELSEIF(setting%RHSdmat)THEN
   RHS1 = setting%RHSdmatAOindex1
   RHS2 = setting%RHSdmatAOindex2
   nbastRHS1 = INTINPUT%AOdim(RHS1)
   nbastRHS2 = INTINPUT%AOdim(RHS2)
   nDmatRHS = setting%nDmatRHS
   useAO1 = .NOT.INTINPUT%AO(RHS1)%p%empty
   useAO2 = .NOT.INTINPUT%AO(RHS2)%p%empty
   IF(INTINPUT%DO_EXCHANGE.AND. (matrix_type .EQ. mtype_unres_dense))THEN
      call LSquit(' excahnge error unres in ls_atta',-1)
   ELSE
      call mem_alloc(setting%DfullRHS,nbastRHS1,nbastRHS2,setting%nDmatRHS)
      setting%RHSdalloc = .TRUE.
      CALL ls_mat_retrive_block(setting%DmatRHS,setting%DfullRHS,setting%ndmatRHS,&
           &                      nbastRHS1,nbastRHS2,1,1,.FALSE.)
      CALL attachDmatToInput(INTINPUT,Setting%DfullRHS,nbastRHS1,nbastRHS2,nDmatRHS,'RHS')
      INTINPUT%uselst_DRHS = .TRUE.
      CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dRHS,Setting%DfullRHS,INTINPUT%AO(RHS1)%p,&
           &INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,lupri)
   ENDIF
ENDIF

END SUBROUTINE ls_setIntegralInputDensities

!> \brief free the densities in the integralinput structure
!> \author S. Reine and T. Kjaeergaard
!> \date 2010
!> \param intinput the integral input
SUBROUTINE ls_freeIntegralInputDensities(intinput)
implicit none
TYPE(INTEGRALINPUT) :: INTINPUT
!
IF (INTINPUT%uselst_DLHS) THEN
   CALL lstensor_free(INTINPUT%lst_DLHS)
   INTINPUT%uselst_DLHS = .FALSE.
ENDIF
IF (INTINPUT%uselst_DRHS)THEN
   CALL lstensor_free(INTINPUT%lst_dRHS)
   INTINPUT%uselst_DRHS = .FALSE.
ENDIF

END SUBROUTINE ls_freeIntegralInputDensities

!SUBROUTINE ls_init_IntCalcCounter()
!  use thermite_distribute_exchange
!implicit none
!
!call th_init_IntCalcCounter()
!
!END SUBROUTINE ls_init_IntCalcCounter

!SUBROUTINE ls_get_IntCalcCounter(nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO)
!  use thermite_distribute_exchange
!implicit none
!integer :: nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO
!
!call th_get_IntCalcCounter(nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO)
!
!END SUBROUTINE ls_get_IntCalcCounter

!> \brief Sets the four fragments for point to the four molecules
!> \author S. Reine
!> \date 2010-08-23
!> \param setting The ls-settings
SUBROUTINE ls_setDefaultFragments(setting)
implicit none
TYPE(LSSETTING),intent(inout) :: setting
Integer :: iao

DO iao=1,4
  setting%fragment(iao)%p => setting%molecule(iao)%p
ENDDO
setting%sameFrag = setting%sameMol
END SUBROUTINE ls_setDefaultFragments

END MODULE ls_Integral_Interface

#ifdef VAR_LSMPI
SUBROUTINE lsmpi_getIntegrals_masterToSlave(AO1,AO2,AO3,AO4,ndmat,nderiv,Oper,Spec,intType,SETTING,LUPRI,LUERR)
use lsmpi_mod
use infpar_module
use typedef
implicit none
Integer,intent(in)            :: LUPRI,LUERR,ndmat,nderiv
Character*(*)                 :: AO1,AO2,AO3,AO4,Oper,Spec,intType
Type(LSSETTING),intent(inout) :: SETTING
!
Integer :: dim

dim = len(AO1)
CALL ls_mpibcast(dim,infpar%master)
CALL ls_mpibcast(AO1,dim,infpar%master)
dim = len(AO2)
CALL ls_mpibcast(dim,infpar%master)
CALL ls_mpibcast(AO2,dim,infpar%master)
dim = len(AO3)
CALL ls_mpibcast(dim,infpar%master)
CALL ls_mpibcast(AO3,dim,infpar%master)
dim = len(AO4)
CALL ls_mpibcast(dim,infpar%master)
CALL ls_mpibcast(AO4,dim,infpar%master)
dim = len(Oper)
CALL ls_mpibcast(dim,infpar%master)
CALL ls_mpibcast(Oper,dim,infpar%master)
dim = len(Spec)
CALL ls_mpibcast(dim,infpar%master)
CALL ls_mpibcast(Spec,dim,infpar%master)
dim = len(intType)
CALL ls_mpibcast(dim,infpar%master)
CALL ls_mpibcast(intType,dim,infpar%master)

CALL ls_mpibcast(ndmat,infpar%master)
CALL ls_mpibcast(nderiv,infpar%master)
CALL ls_mpibcast(LUPRI,infpar%master)
CALL ls_mpibcast(LUERR,infpar%master)
CALL mpicopy_setting(setting)

END SUBROUTINE lsmpi_getIntegrals_masterToSlave

SUBROUTINE lsmpi_getIntegrals_Slave
use infpar_module
use lsmpi_mod
use typedef
use ls_Integral_Interface
implicit none
Integer         :: LUPRI,LUERR,ndmat,nderiv
Character*(80)  :: AO1,AO2,AO3,AO4,Oper,Spec,intType
Type(LSSETTING) :: SETTING
!
Integer :: dimAO1,dimAO2,dimAO3,dimAO4,dimOper,dimSpec,dimIntType

write(*,*) 'debug:entered lsmpi_getIntegrals_Slave',infpar%mynum

CALL ls_mpibcast(dimAO1,infpar%master)
CALL ls_mpibcast(AO1,dimAO1,infpar%master)
CALL ls_mpibcast(dimAO2,infpar%master)
CALL ls_mpibcast(AO2,dimAO2,infpar%master)
CALL ls_mpibcast(dimAO3,infpar%master)
CALL ls_mpibcast(AO3,dimAO3,infpar%master)
CALL ls_mpibcast(dimAO4,infpar%master)
CALL ls_mpibcast(AO4,dimAO4,infpar%master)
CALL ls_mpibcast(dimOper,infpar%master)
CALL ls_mpibcast(Oper,dimOper,infpar%master)
CALL ls_mpibcast(dimSpec,infpar%master)
CALL ls_mpibcast(Spec,dimSpec,infpar%master)
CALL ls_mpibcast(dimIntType,infpar%master)
CALL ls_mpibcast(intType,dimIntType,infpar%master)

CALL ls_mpibcast(ndmat,infpar%master)
CALL ls_mpibcast(nderiv,infpar%master)
CALL ls_mpibcast(LUPRI,infpar%master)
CALL ls_mpibcast(LUERR,infpar%master)
write(*,*) 'debug:calling mpicopy_setting',infpar%mynum
CALL mpicopy_setting(setting)
write(*,*) 'debug:calling ls_getIntegrals',infpar%mynum
call ls_getIntegrals(AO1(1:dimAO1),AO2(1:dimAO2),AO3(1:dimAO3),AO4(1:dimAO4),ndmat,nderiv,&
     &               Oper(1:dimOper),Spec(1:dimSpec),intType(1:dimIntType),SETTING,LUPRI,LUERR)
write(*,*) 'debug:ending ls_getIntegrals',infpar%mynum

END SUBROUTINE lsmpi_getIntegrals_Slave
#endif

