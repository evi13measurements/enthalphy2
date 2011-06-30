!> @file 
!> Contains integralinput structure which contain all info required to do an integral evaluation.  
MODULE integral_type
use precision
use AO_type
use lstensor_operationsmod
!*****************************************
!*
!* OBJECT CONTAINING INFORMATION ABOUT
!* THE MAIN-INTEGRAL SETTINGS OBJECT
!*   - used in ThermiteDriver.f90
!*
!*****************************************
TYPE INTEGRALINPUT
REAL(REALK) :: Integralthreshold
TYPE(AOITEMPOINTER)   :: AO(4)
Integer               :: AOdim(4)
Integer               :: currentFragment(4)
Real(realk),pointer   :: DMAT_LHS(:,:,:)
Real(realk),pointer   :: DMAT_RHS(:,:,:)
type(lstensor)        :: LST_DLHS
type(lstensor)        :: LST_DRHS
type(lstensor),pointer :: LST_GAB_LHS
type(lstensor),pointer :: LST_GAB_RHS
type(lstensor),pointer :: LST_pGAB_LHS
type(lstensor),pointer :: LST_pGAB_RHS
logical               :: PS_RHSGABusePointer
logical               :: CS_RHSGABusePointer
logical               :: uselst_DRHS
logical               :: uselst_DLHS
Real(realk),pointer   :: GAB_LHS(:,:)
Real(realk),pointer   :: GAB_RHS(:,:)
Real(realk),pointer   :: pGAB_LHS(:,:)
Real(realk),pointer   :: pGAB_RHS(:,:)
!TYPE(LSMATRIX),pointer  :: DMAT_LHS
!TYPE(LSMATRIX),pointer  :: DMAT_RHS
LOGICAL       :: LHS_DMAT,RHS_DMAT
INTEGER       :: FTUVmaxprim
INTEGER       :: LU_MMDATA,MMINDEX
INTEGER       :: LU_MMDATR
INTEGER       :: NDMAT_LHS,NDMAT_RHS
LOGICAL       :: DRHS_SYM,DLHS_SYM
INTEGER       :: NDIM_LHS(2),NDIM_RHS(2)
LOGICAL       :: sameLHSaos,sameRHSaos,sameODs
INTEGER       :: CENTERS,ndim
LOGICAL       :: DO_PASSES
integer       :: maxPasses
LOGICAL       :: HIGH_RJ000_ACCURACY
integer       :: MMunique_ID1
integer       :: MMunique_ID2
integer       :: MMstartA
integer       :: MMstartB
LOGICAL       :: PUREFMM,SKIP_FMM
LOGICAL       :: DO_COULOMB,DO_EXCHANGE,DO_FMM,DO_MULMOM
LOGICAL       :: DO_FOCK,DO_ENERGY,DO_JENGINE,DO_LINK,DO_DALINK
LOGICAL       :: CS_int !Cauchy-Schwarz integrals
LOGICAL       :: PS_int !Primitive Cauchy-Schwarz integrals
LOGICAL       :: setETUVoutside !Sets up ETUV-tensor before integral-loop
LOGICAL       :: sphericalEcoeff !ETUV-tensor in solid-harmonical form
LOGICAL       :: orderAngPrim !Performs integration with angular components first and primitive components last (Ang,Prim)
LOGICAL       :: orderPQ      !Performs integration with primitives P before primtive Q (P,Q)
!Cauchy-Schwarz screening
LOGICAL       :: CS_SCREEN 
Real(realk)   :: CS_THRESHOLD
!Cauchy-Schwarz screening
LOGICAL       :: PS_SCREEN 
Real(realk)   :: PS_THRESHOLD
!Overlap-extent screening (for overlap integrals)
LOGICAL       :: OE_SCREEN 
Real(realk)   :: OE_THRESHOLD
!Screen OD-batches by AO-batch extent
LOGICAL       :: OD_SCREEN 
Real(realk)   :: OD_THRESHOLD
!Screen integrals by their non-classical extent
LOGICAL       :: NonClassical_SCREEN 
Real(realk)   :: MM_SCREENTHR = 1.d-10
Integer       :: MM_TLMAX     = 8
Integer       :: MM_LMAX      = 20
Logical       :: MM_NOONE
!Primitive-Screening maximum element of prim GAB matrix
Real(realk)   :: CS_MAXELM_LHS
Real(realk)   :: CS_MAXELM_RHS
Real(realk)   :: PS_MAXELM_LHS
Real(realk)   :: PS_MAXELM_RHS
Real(realk)   :: CoulombFactor
Real(realk)   :: exchangeFactor
Real(realk)   :: ORIGO(3)
!Derivative info
INTEGER       :: nderivP
INTEGER       :: nderivQ
INTEGER       :: derOrderP
INTEGER       :: derOrderQ
INTEGER       :: derivOrder
Logical       :: DO_GRADIENT

LOGICAL       :: AddToIntegral
!Attenuation parameters 
Real(realk)   :: ATTomega 
Real(realk)   :: ATTalpha
Real(realk)   :: ATTbeta
LOGICAL       :: ATTFACTOR
Character(len=80)     :: operator
LOGICAL       :: DECPACKED

END TYPE INTEGRALINPUT

END MODULE integral_type
