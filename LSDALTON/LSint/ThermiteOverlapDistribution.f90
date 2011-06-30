!> @file
!> Contains the orbital and overlap types, and subroutines to set these up
MODULE Thermite_OD
use typedef
use SphCart_matrices
!use ODbatches
use OD_type
use ls_util

TYPE Allocitem
Integer                :: maxJLHS       !used for WORK
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
Integer                :: maxETUVlenLHS

Integer                :: maxJRHS    
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
Integer                :: maxETUVlenRHS
END TYPE Allocitem

!Simen Should remove all explicit dependence on IntegralItem from this file
!      and the move IntegralItem to ThermiteIntegrals.f90
!****** INTEGRAL
TYPE Integralitem
Real(realk),dimension(:),pointer    :: IN
Real(realk),dimension(:),pointer    :: OUT
Real(realk),dimension(:),pointer    :: RTUV
Real(realk),dimension(:),pointer    :: integrals
!
Real(realk),dimension(:,:),pointer     :: Wtuv       !nPrimitives:ntuv
Real(realk),dimension(:,:,:),pointer   :: tuvTUV     !ntuvQ(i):ntuvP:nPrimPQ
!Real(realk),dimension(:,:,:),pointer    :: TUVQ(:,:,:)     !ntuvP:nPrimP:nOrbQ
Real(realk),dimension(:,:),pointer     :: PQ         !nOrbQ:nOrbP
!Generic integral intermadiates used of contraction with Ecoefficients, 
!Spherical transformation, and contraction to contracted basis
Real(realk),dimension(:,:,:),pointer    :: IntegralIN
Real(realk),dimension(:,:,:),pointer    :: IntegralOUT
INTEGER                :: nTUV
INTEGER                :: nEFG !multipole orders
INTEGER                :: nOrb
INTEGER                :: nAng
INTEGER                :: nPrim
INTEGER                :: nDeriv
LOGICAL                :: Jengine
TYPE(TUVitem),pointer  :: TUV
END TYPE Integralitem

TYPE Orbital
!Only one of these types can be true
Logical                       :: TYPE_Empty
Logical                       :: TYPE_Hermite
Logical                       :: TYPE_Cartesian
Logical                       :: TYPE_Nucleus
!-----------------------------------
Logical                       :: spherical
Logical                       :: FTUVorb
Integer                       :: maxAngmom
Integer                       :: nAngmom
Integer                       :: nPrimitives
Integer                       :: nPasses
Integer                       :: totOrbitals
Integer                       :: maxContracted
Integer                       :: CCidentifier
Real(realk)                   :: center(3)       !Cartesian center
!Dimension according to nPasses
integer,pointer               :: atom(:)
!Dimensions according to nAngmom
integer,pointer               :: angmom(:) !1:nangmom
integer,pointer               :: nContracted(:)
integer,pointer               :: startOrbital(:,:)
integer,pointer               :: startLocOrb(:)
integer,pointer               :: batch(:)
integer,pointer               :: startprimOrbital(:,:)
integer,pointer               :: nOrbComp(:)
integer,pointer               :: nPrimOrbComp(:)
integer,pointer               :: nOrbitals(:)
!Dimensions according to nPrimitives
real(realk),pointer           :: exponents(:)
!Dimensions according to nPrimitives,maxContracted,nAngmom
type(lsmatrixpointer),pointer :: CC(:)
TYPE(SPHMATPOINTER),pointer   :: SPH_MAT(:)
END TYPE Orbital

!****** OVERLAPPOINTER
TYPE Overlappointer
TYPE(Overlap),pointer :: p
END TYPE Overlappointer


!****** OVERLAP
TYPE Overlap
!only one of these types can be true
Logical             :: TYPE_Empty
Logical             :: TYPE_Hermite
Logical             :: TYPE_Hermite_single
Logical             :: TYPE_Cartesian
Logical             :: TYPE_Cartesian_single
Logical             :: TYPE_Nucleus
Logical             :: TYPE_FTUV
!
Logical             :: single
real(realk)         :: ODcenter(3)
real(realk)         :: ODextent
Logical             :: sphericalEcoeff
TYPE(Orbital)       :: orbital1,orbital2
logical             :: sameAO
Integer             :: nAngmom
Integer             :: nPrimitives
Integer             :: totOrbitals
Integer             :: nTUV
Integer             :: maxContracted
Real(realk)         :: maxGab
!Real(realk)         :: maxPrimGab
!minAngmom and maxAngmom represents the actual angular momenta of the orbitals
Integer             :: minAngmom
Integer             :: maxAngmom
!The order of differentiation (0 for energy, 1 for gradients, 2 for hessians etc.)
Integer             :: derivOrder
!The number of differentiated components (for gradients: 3 for a single AO gradient 
!and 6 for a product of two AOs)
Integer             :: nDerivComp 
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
Logical             :: ETUVisAlloc
Logical             :: FTUVisAlloc
Integer,pointer     :: lenETUV(:)
Real(realk),pointer :: ETUV(:) ! ntuv, ijk, nprim, iAngmom
Real(realk),pointer :: FTUV(:,:,:)

!Pass specific variables (passes are used to increase perforance by 
!collecting OD-batches of similar structure - this decrease the time
!spent in the different contraction/matrix multiplies steps, as 
!well as reducing the overhead related to memory/cache loadings and
!subroutine calls).
Integer             :: nPasses
Integer             :: passType
END TYPE Overlap

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

CONTAINS

!> \brief One of two routines to set up the overlap (and the two corresponding orbitals) from the OD-batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009-??-??
!> \param P The overlap
!> \param np The number of (significant) primitives
!> \param Input The integral input
!> \param sharedTUV The TUV indeces
!> \param Integral The integral specifications
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param ODB The OD-batch information
!> \param iElectron The electron (1 = 'LHS', 2 = 'RHS')
!> \param lupri The default print unit
!> \param iprint The print level (the higher the more information)
!> \param side The side ('LHS' or 'RHS')
!> Simen: either remove side or iElectron!!!
SUBROUTINE SET_OVERLAP(P,np,Input,sharedTUV,Integral,Alloc,ODB,IELECTRON,LUPRI,IPRINT,SIDE)
use memory_handling
Implicit none
TYPE(IntegralInput),target,intent(in) :: Input
TYPE(IntegralItem),intent(inout)      :: Integral
TYPE(AllocItem),intent(inout)         :: Alloc
TYPE(Overlap),intent(out)             :: P
TYPE(ODBATCH),intent(in)              :: ODB
Character*(*)                         :: side
Integer,intent(in)                    :: IELECTRON,LUPRI,IPRINT
Integer,intent(inout)                 :: np
TYPE(TUVitem),intent(in)              :: SharedTUV
!
Integer             :: nderiv,idir,i1,i2,i12,start1,start2,end1,end2,orb1,orb2,start
integer             :: l1,l2,ijk1,ijk2,ijk,ijkcart,maxijk,iangmom
Real(realk)         :: e1,e2,d2,maxGab,maxPrimGab,maxPrimGabElm,signP,DMATmax
!Character(len=80)   :: t1,t2
Type(lstensor),pointer :: DMAT2,GAB2,primGAB2
LOGICAL             :: LHS, DMATscreen,useFTUV,doscreen,screen
Integer             :: ndmat,idmat,indexETUV,l,nETUV,nTUV,maxangmom,minangmom,maxnp
Integer,pointer     :: lenEtuv(:)
Integer :: batchA,batchB,atomA,atomB,sAA,sA,sB,iSA,iSB,Dindex,elms,dimA,dimB,Gindex
real(realk),pointer :: primgabmat(:,:)
!
DMATscreen = .FALSE.
NULLIFY(DMAT2)
NULLIFY(GAB2)
NULLIFY(primGAB2)
!
SELECT CASE(SIDE)
CASE('LHS')
   ndmat = INPUT%NDMAT_LHS
   IF(INPUT%PS_SCREEN)THEN
      primGAB2 => INPUT%LST_pGAB_LHS
      maxPrimGabElm=INPUT%PS_MAXELM_RHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
      GAB2 => INPUT%LST_GAB_LHS      
   ENDIF
   LHS=.TRUE.
   NDERIV=INPUT%NDERIVP
   P%derivOrder = INPUT%derOrderP
CASE('RHS')
   ndmat = INPUT%NDMAT_RHS
   IF(INPUT%PS_SCREEN)THEN
      primGAB2 => INPUT%LST_pGAB_RHS
      maxPrimGabElm=INPUT%PS_MAXELM_LHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
      GAB2 => INPUT%LST_GAB_RHS
      IF (Input%DO_JENGINE) THEN
        DMATscreen = .TRUE.
        DMAT2 => INPUT%LST_DRHS
      ENDIF
   ENDIF
   LHS=.FALSE.
   NDERIV=INPUT%NDERIVQ
   P%derivOrder = INPUT%derOrderQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_OVERLAP =',SIDE
   CALL LSQUIT('Wrong case in SET_OVERLAP',lupri)
END SELECT

 CALL ALLOCATE_ORBITAL(P%orbital1,ODB%AO(1)%p%nprimitives,1,ODB%AO(1)%p%nAngmom,ODB%AO(1)%p%maxAngmom,.FALSE.)
 CALL ALLOCATE_ORBITAL(P%orbital2,ODB%AO(2)%p%nprimitives,1,ODB%AO(2)%p%nAngmom,ODB%AO(2)%p%maxAngmom,.FALSE.)
 CALL SET_ORBITAL(P%orbital1,ODB%AO(1)%p,integral,LUPRI)
 CALL SET_ORBITAL(P%orbital2,ODB%AO(2)%p,integral,LUPRI)

!Determine overlap type
 P%type_Empty = P%orbital1%type_Empty .AND. P%orbital2%type_Empty 
 IF(P%type_Empty)THEN
    IF(Input%operator(1:6) .EQ. 'Mulmom' .AND. (.NOT.LHS) )THEN
       ! this is done in order to use screening on multipole moments used for FMM
       P%maxGab = 1.D0
    ENDIF
 ENDIF
 P%type_Hermite_single = (P%orbital1%type_Hermite .AND. P%orbital2%type_Empty).OR.&
      &(P%orbital1%type_Empty .AND. P%orbital2%type_Hermite)
 P%type_Hermite = P%orbital1%type_Hermite .AND. P%orbital2%type_Hermite
 P%type_Cartesian_single = (P%orbital1%type_Cartesian .AND. P%orbital2%type_Empty).OR.&
      &(P%orbital1%type_Empty .AND. P%orbital2%type_Cartesian)
 P%type_Cartesian = P%orbital1%type_Cartesian .AND. P%orbital2%type_Cartesian
 P%type_Nucleus = (P%orbital1%type_Nucleus .AND. P%orbital2%type_Empty)&
      &.OR.(P%orbital1%type_Empty .AND. P%orbital2%type_Nucleus)
 P%type_FTUV = .FALSE.

 P%single = P%type_Hermite_single.OR.P%type_Cartesian_single.OR.P%type_Nucleus
 CALL getDerivComp(P%nDerivComp,P%derivOrder,P%type_empty,P%single)
   
 P%sphericalEcoeff = Input%sphericalEcoeff.AND.(.NOT.P%type_hermite_single)
 useFTUV = INPUT%DO_JENGINE .AND. IELECTRON.EQ.2
 P%sphericalEcoeff = P%sphericalEcoeff .AND.(.NOT.useFTUV).AND.&
     &               (.NOT.INPUT%operator(1:7) .EQ. 'Kinetic' .AND. LHS) 
 P%ODcenter = ODB%ODcenter
 P%ODextent = ODB%ODextent
 P%sameAO = ODB%sameAO
 P%ETUVisSet = .FALSE.
 CALL GET_NPRIMITIVES(np,P,Input,Side,lupri)!can be done outside
 maxnp = max(np,P%orbital1%nPrimitives,P%orbital2%nPrimitives)
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
 indexETUV = 1
 minAngmom   = 99
 maxangmom = 0
 i12 = 0
 CALL MEM_ALLOC(lenETUV,P%orbital1%nAngmom*P%orbital2%nAngmom)
 DO i1=1,P%orbital1%nAngmom
    start2 = 1
    IF (P%sameAO) start2 = i1
    DO i2=start2,P%orbital2%nAngmom
       i12  = i12 + 1
       l1   = P%orbital1%angmom(i1)
       l2   = P%orbital2%angmom(i2)
       l  = l1+l2+P%derivOrder
       CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,P%sphericalEcoeff,P%derivOrder,P%single)
       nTUV = (l+1)*(l+2)*(l+3)/6
       lenETUV(i12) =  nTUV*ijk*P%nPrimitives
       indexETUV = indexETUV + nTUV*ijk*P%nPrimitives
       maxangmom = MAX(maxangmom,l)
       minangmom = MIN(minangmom,l)
    ENDDO
 ENDDO

 nETUV = indexETUV - 1
 P%maxAngmom = maxAngmom
 P%minAngmom = minAngmom 
 P%startAngmom = 0
 IF (P%type_hermite_single) P%startAngmom = minAngmom + P%derivOrder
 IF( INPUT%operator(1:7) .EQ. 'Kinetic' .AND. LHS) THEN
    P%startAngmom = P%startAngmom + 2
    P%endAngmom   = maxAngmom + 2 + P%derivOrder
 ELSE
    P%endAngmom   = maxAngmom + P%derivOrder
 ENDIF  
 P%nTUV = (P%endAngmom+1)*(P%endAngmom+2)*(P%endAngmom+3)/6 - &
    &     P%startAngmom*(P%startAngmom+1)*(P%startAngmom+2)/6

 CALL ALLOCATE_OVERLAP(P,np,np,1,P%nAngmom,useFTUV,P%nTUV,Input%NDMAT_RHS,P%type_hermite_single,&
      & nETUV,INPUT%setETUVoutside)

 DO iAngmom = 1,P%nAngmom
    P%lenETUV(iAngmom) = lenETUV(iAngmom)  
 ENDDO
 CALL MEM_DEALLOC(lenETUV)
 d2 = 0.0d0
 DO idir=1,3
   P%distance12(idir,1) = P%orbital1%center(idir)-P%orbital2%center(idir)
   d2 = d2 + P%distance12(idir,1) * P%distance12(idir,1)
 ENDDO
 P%squaredDistance = d2
!
! IF((INPUT%PS_SCREEN).AND.(.NOT.P%type_Empty)) THEN
!   AtomA = P%orbital1%atom(1)
!   AtomB = P%orbital2%atom(1)
!   BatchA = P%orbital1%batch(1)
!   BatchB = P%orbital2%batch(1)
!   Gindex = primGAB2%INDEX(atomA,atomB,1,1)
! ENDIF

 doscreen = ((INPUT%PS_SCREEN).AND.(.NOT.P%type_Empty))
! maxprimgab = 0.d0
 IF(doscreen)then
  call mem_alloc(primgabmat,P%orbital1%nPrimitives,P%orbital2%nPrimitives)
  AtomA = P%orbital1%atom(1)
  AtomB = P%orbital2%atom(1)
  BatchA = P%orbital1%batch(1)
  BatchB = P%orbital2%batch(1)
  Gindex = primGAB2%INDEX(P%orbital1%atom(1),P%orbital2%atom(1),1,1)
  DO i1=1,P%orbital1%nPrimitives
   DO i2=1,P%orbital2%nPrimitives
    primgabmat(i1,i2) = &
         &  primGAB2%LSAO(Gindex)%BATCH(BatchA,BatchB,1,1)%elms(i1 + (i2-1)*P%orbital1%nPrimitives)
   ENDDO
  ENDDO
 ENDIF
 i12 = 0
 DO i1=1,P%orbital1%nPrimitives
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nPrimitives
     IF (doscreen) THEN
!       maxprimgab = MAX(maxprimgab,primgabmat(i1,i2))
       screen = primgabmat(i1,i2) .LT. INPUT%PS_THRESHOLD/maxPrimGabElm
     ELSE
       screen = .FALSE.
     ENDIF
     IF (.NOT.screen) THEN
       i12 = i12 + 1
       e1  = P%orbital1%exponents(i1)
       e2  = P%orbital2%exponents(i2)
       P%exponents(i12)        = e1 + e2
       IF(Input%operator(1:6) .EQ. 'Nucrep')THEN
         P%reducedExponents(i12) = e1
         P%center(1,i12)         = P%orbital1%center(1)
         P%center(2,i12)         = P%orbital1%center(2)
         P%center(3,i12)         = P%orbital1%center(3)
       ELSEIF (P%exponents(i12) .NE. 0.0d0) THEN
         P%reducedExponents(i12) = e1*e2/(e1+e2)
         P%center(1,i12)         = (e1*P%orbital1%center(1) + e2*P%orbital2%center(1))/(e1+e2)
         P%center(2,i12)         = (e1*P%orbital1%center(2) + e2*P%orbital2%center(2))/(e1+e2)
         P%center(3,i12)         = (e1*P%orbital1%center(3) + e2*P%orbital2%center(3))/(e1+e2)
       ELSE
         P%reducedExponents(i12) = 0.0d0
         P%center(1,i12)         = 0.0d0
         P%center(2,i12)         = 0.0d0
         P%center(3,i12)         = 0.0d0
       ENDIF 
       IF ( P%single) THEN 
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
   CALL LSQUIT('Error: Mismatch between get_nPrimitives and i12 in set_overlap',lupri)
 ENDIF

IF(doscreen) call mem_dealloc(primgabmat)
IF(INPUT%PS_SCREEN)THEN
!   P%maxPrimGab = maxPrimGab
   NULLIFY(primGAB2)
ENDIF


!CAREFUL
!Initialize to twice some maximal value + 1
 P%minAngmom   = 99
!CAREFUL
 P%maxAngmom = 0
 P%totOrbitals = 0
 i12 = 0
 maxijk = 0
 DO i1=1,P%orbital1%nAngmom
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nAngmom
     i12  = i12 + 1
     l1   = P%orbital1%angmom(i1)
     l2   = P%orbital2%angmom(i2)
     CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,P%sphericalEcoeff,P%derivOrder,P%single)
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
  ENDDO
ENDDO


IF (Input%CS_SCREEN .AND. (.NOT.P%type_Empty)) THEN
   AtomA = P%orbital1%atom(1)
   BatchA = P%orbital1%batch(1)
   AtomB = P%orbital2%atom(1)
   BatchB = P%orbital2%batch(1)
   Gindex = GAB2%INDEX(atomA,atomB,1,1)
   IF (DMATscreen) Dindex = DMAT2%INDEX(atomA,atomB,1,1)
   dimA = P%orbital1%nOrbComp(1)*P%orbital1%nContracted(1)
   DO SAA=2,P%orbital1%nAngmom
      dimA = dimA+P%orbital1%nOrbComp(SAA)*P%orbital1%nContracted(SAA)
   ENDDO
   dimB = P%orbital2%nOrbComp(1)*P%orbital2%nContracted(1)
   DO SAA=2,P%orbital2%nAngmom
      dimB = dimB+P%orbital2%nOrbComp(SAA)*P%orbital2%nContracted(SAA)
   ENDDO
   maxGab = 0.0d0
   DO i1=1,P%orbital1%nAngmom
      start = 1
      IF (P%sameAO) start = i1
      DO i2=start,P%orbital2%nAngmom
         start1 = P%orbital1%startOrbital(i1,1)
         end1   = start1 + P%orbital1%nOrbitals(i1) - 1
         start2 = P%orbital2%startOrbital(i2,1)
         end2   = start2 + P%orbital2%nOrbitals(i2) - 1
         DO orb1=start1,end1
            DO orb2=start2,end2
               sA = 0
               DO SAA=2,i1
                  sA = sA+P%orbital1%nOrbComp(SAA-1)*P%orbital1%nContracted(SAA-1)
               ENDDO
               sB = 0
               DO SAA=2,i2
                  sB = sB+P%orbital2%nOrbComp(SAA-1)*P%orbital2%nContracted(SAA-1)
               ENDDO
               iSA = (orb1-start1)+1+sA
               iSB = (orb2-start2)+1+sB
               IF (DMATscreen) THEN
                  DMATmax = 0.0d0
                  DO idmat=1,ndmat
                     elms = iSA + (iSB-1)*dimA + (idmat-1)*dimA*dimB
                     DMATmax = max(DMATmax,abs(DMAT2%LSAO(Dindex)%BATCH(batchA,batchB,1,1)%elms(elms)))
                  ENDDO
                  elms = iSA + (iSB-1)*dimA
                  maxGab = MAX(maxGab,GAB2%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%elms(elms)*DMATmax)
               ELSE
                  elms = iSA + (iSB-1)*dimA
                  maxGab = MAX(maxGab,GAB2%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%elms(elms))
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   P%maxGab = maxGab
!   NULLIFY(GAB)
   NULLIFY(GAB2)
ENDIF

 P%startAngmom = 0
 IF (P%type_hermite_single) P%startAngmom = P%minAngmom + P%derivOrder

 IF( INPUT%operator(1:7) .EQ. 'Kinetic' .AND. LHS) THEN
    P%startAngmom = P%startAngmom + 2
    P%endAngmom   = P%maxAngmom + 2 + P%derivOrder
 ELSE
    P%endAngmom   = P%maxAngmom + P%derivOrder
 ENDIF  

 P%nTUV = (P%endAngmom+1)*(P%endAngmom+2)*(P%endAngmom+3)/6 - &
    &     P%startAngmom*(P%startAngmom+1)*(P%startAngmom+2)/6

 IF (IELECTRON.EQ.1) THEN
   signP=1.0d0
   Alloc%maxJLHS    = max(P%maxAngmom,Alloc%maxJLHS)
   Alloc%maxPrimLHS    = max(maxnp,Alloc%maxPrimLHS)
   Alloc%maxPrimTUVLHS = max(maxnp*P%nTUV,P%totOrbitals,maxnp*maxijk,Alloc%maxPrimTUVLHS)
   Alloc%maxnAngLHS     = max(P%nAngmom,Alloc%maxnAngLHS)
   Alloc%maxContLHS    = max(P%Orbital1%maxContracted*P%Orbital2%maxContracted,Alloc%maxContLHS)
   Alloc%maxangmomLHS = max(P%Orbital1%maxAngmom,P%Orbital2%maxAngmom,Alloc%maxangmomLHS)
   Alloc%maxTUVLHS = max(P%nTUV,Alloc%maxTUVLHS)
   Alloc%maxTotOrbLHS = max(P%totOrbitals,Alloc%maxTotOrbLHS)
   Alloc%maxijkLHS = max(maxijk,Alloc%maxijkLHS)
   Alloc%maxETUVlenLHS = max(P%lenETUV(P%nAngmom),Alloc%maxETUVlenLHS)
 ELSE
   signP=-1.0d0
   Alloc%maxJRHS    = max(P%maxAngmom,Alloc%maxJRHS)
   Alloc%maxPrimRHS    = max(maxnp,Alloc%maxPrimRHS)
   Alloc%maxPrimTUVRHS = max(maxnp*P%nTUV,P%totOrbitals,maxnp*maxijk,Alloc%maxPrimTUVRHS)
   Alloc%maxnAngRHS     = max(P%nAngmom,Alloc%maxnAngRHS)
   Alloc%maxContRHS    = max(P%Orbital1%maxContracted*P%Orbital2%maxContracted,Alloc%maxContRHS)
   Alloc%maxangmomRHS = max(P%Orbital1%maxAngmom,P%Orbital2%maxAngmom,Alloc%maxangmomRHS)
   Alloc%maxTUVRHS = max(P%nTUV,Alloc%maxTUVRHS)
   Alloc%maxTotOrbRHS = max(P%totOrbitals,Alloc%maxTotOrbRHS)
   Alloc%maxijkRHS = max(maxijk,Alloc%maxijkRHS)
   Alloc%maxETUVlenRHS = max(P%lenETUV(P%nAngmom),Alloc%maxETUVlenRHS)
 ENDIF

IF (INPUT%DO_JENGINE .AND. IELECTRON.EQ.2) THEN
  CALL SET_FTUV(P,INPUT,Integral,LUPRI,IPRINT)
! Change Overlap P to confirm with FTUV
  P%TYPE_FTUV = .TRUE.
  P%TYPE_Empty = .FALSE.
  P%TYPE_Hermite = .FALSE.
  P%TYPE_Hermite_single = .FALSE.
  P%TYPE_Cartesian = .FALSE.
  P%TYPE_Cartesian_single = .FALSE.
  P%TYPE_Nucleus = .FALSE.
  P%sphericalEcoeff = .FALSE.
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
  P%totOrbitals    = ndmat
  P%nOrbComp(1)    = 1
  P%nContracted(1) = 1
ELSE
!GET E COEFFICIENTS
 IF (.NOT.P%type_hermite_single) THEN
   IF (INPUT%setETUVoutside) THEN
    IF(LHS)THEN
       CALL AttachETUVtoOverlap(P,signP,SharedTUV,.TRUE.,LUPRI,IPRINT)
    ELSE
       CALL AttachETUVtoOverlap(P,signP,SharedTUV,Input%orderAngPrim,LUPRI,IPRINT)
    ENDIF
    IF(IELECTRON.EQ.2)THEN
       Alloc%maxPrimTUVijkRHS = max(nETUV,Alloc%maxPrimTUVijkRHS) 
    ELSE
       Alloc%maxPrimTUVijkLHS = max(nETUV,Alloc%maxPrimTUVijkLHS) 
    ENDIF
   ENDIF
 ENDIF
ENDIF

IF(INPUT%CS_SCREEN)THEN
   NULLIFY(GAB2)
ENDIF
END SUBROUTINE SET_OVERLAP

!> \brief Returns the number of two-center derivative components for a given derivative order
!> \author S. Reine
!> \date 2010-03-05
!> \param nDerivComp The number of derivative components
!> \param derOrder The derivative order
!> \param empty True if empty overlap
!> \param single True if only a single AO
SUBROUTINE getDerivComp(nDerivComp,derOrder,empty,single)
implicit none
Integer,intent(OUT) :: nDerivComp
Integer,intent(IN)  :: derOrder
Logical,intent(IN)  :: empty,single
!
Integer :: iDer,jDer
!
!Consistency testing
IF (empty.AND.single) CALL LSQUIT('Error in getDerivComp. Single and empty!',-1)
IF (derOrder.LT.0) CALL LSQUIT('Error in getDerivComp. derOrder<0',-1)
!
IF (empty) THEN
  nDerivComp = 1
ELSE IF (single) THEN
  nDerivComp = (derOrder+1)*(derOrder+2)/2
ELSE 
  nDerivComp = 0
  DO iDer=0,derOrder
    jDer = derOrder-iDer
    nDerivComp = nDerivComp + (iDer+1)*(iDer+2)*(jDer+1)*(jDer+2)/4
  ENDDO
ENDIF
!
END SUBROUTINE getDerivComp

!> \brief initialise the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P The overlap
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param Input contains info about the integral evaluation requested
!> \param ODbat The OD-batch information
!> \param side The side ('LHS' or 'RHS') 
!> \param iElectron The electron (1 = 'LHS', 2 = 'RHS')
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
SUBROUTINE INIT_OVERLAP(P,Alloc,Input,ODbat,SIDE,IELECTRON,IPRINT,LUPRI)
Implicit none
INTEGER             :: lupri,IPRINT,IELECTRON
TYPE(Overlap)       :: P
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
TYPE(ODITEM)        :: ODbat
Character*(*)       :: side
!TYPE(ODBATCH)       :: ODbat(nbatches)
!
!Character(len=80)   :: t1,t2
Integer             :: NDERIV
Integer             :: maxprim,maxtuv,maxcont,maxnangmom
Integer             :: MAXPRIMTUVIJK,maxA
LOGICAL           :: LHS,hermiteSingle,useFTUV 

SELECT CASE(SIDE)
CASE('LHS')
   LHS=.TRUE.
   NDERIV=INPUT%NDERIVP
CASE('RHS')
   LHS=.FALSE.
   NDERIV=INPUT%NDERIVQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in INIT_OVERLAP =',SIDE
   CALL LSQUIT('Wrong case in INIT_OVERLAP',lupri)
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

CALL ALLOCATE_ORBITAL(P%orbital1,maxprim,1,maxnangmom,maxA,.FALSE.)
CALL ALLOCATE_ORBITAL(P%orbital2,maxprim,1,maxnangmom,maxA,.FALSE.)
 hermiteSingle = (ODbat%BATCH(1)%AO(1)%p%type_Empty .AND. ODbat%BATCH(1)%AO(2)%p%type_Hermite)&
      &.OR.(ODbat%BATCH(1)%AO(1)%p%type_Hermite .AND. ODbat%BATCH(1)%AO(2)%p%type_Empty)
 useFTUV = INPUT%DO_JENGINE .AND. IELECTRON.EQ.2
 CALL ALLOCATE_OVERLAP(P,maxprim,maxprim,1,maxnAngmom,useFTUV,maxTUV,Input%NDMAT_RHS,hermitesingle,&
      & MAXPRIMTUVIJK,INPUT%setETUVoutside)

END SUBROUTINE INIT_OVERLAP

!> \brief allocate the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P The overlap
!> \param nprim1 the number of primitive functions of orbital 1
!> \param nprim2 the number of primitive functions of orbital 2
!> \param npass the number of passes 
!> \param nAng the number of angular moments
!> \param allocFTUV flag to allocate the FTUV
!> \param nTUV the number of angular components 
!> \param ndmat the number of density matrices
!> \param hermitesingle flag to determine if this is a single orbital or 2 orbitals
!> \param MAXPRIMTUVIJK the maximum size of Ecoeffs
!> \param setETUVoutside flag to determine if the Ecoeffs should be done outside the loop
SUBROUTINE ALLOCATE_OVERLAP(P,nprim1,nprim2,npass,nAng,allocFTUV,nTUV,ndmat,hermitesingle,MAXPRIMTUVIJK,setETUVoutside)
use memory_handling
IMPLICIT NONE
TYPE(Overlap)    :: P
INTEGER          :: nprim1,nprim2,npass,nAng,nTUV,ndmat,MAXPRIMTUVIJK
LOGICAL          :: allocFTUV,hermitesingle,setETUVoutside

Call Mem_alloc(P%center,3,nprim1)
Call Mem_alloc(P%exponents,nprim1)        !call mem_alloc(A%Belms,I)
Call Mem_alloc(P%reducedExponents,nprim1)
Call Mem_alloc(P%preExpFac,nprim1)
Call Mem_alloc(P%iprim1,nprim2)
Call Mem_alloc(P%iprim2,nprim2)
CALL MEM_ALLOC(P%distance12,3,npass)
Call Mem_alloc(P%angmom,nAng)
Call Mem_alloc(P%lenETUV,nAng)
Call Mem_alloc(P%indexAng1,nAng)
Call Mem_alloc(P%indexAng2,nAng)
Call Mem_alloc(P%nOrbitals,nAng)
Call Mem_alloc(P%nOrbComp,nAng)
Call Mem_alloc(P%nContracted,nAng)

!$OMP CRITICAL (memory)
mem_allocated_overlap = mem_allocated_overlap + mem_realsize*size(P%center)&
     &+ mem_realsize*size(P%exponents) + mem_realsize*size(P%reducedExponents) &
     &+ mem_realsize*size(P%preExpFac) + mem_intsize*size(P%iprim1) + mem_intsize*size(P%iprim2)&
     &+ mem_intsize*size(P%lenETUV)
mem_allocated_overlap = mem_allocated_overlap + mem_realsize*size(P%distance12)
mem_allocated_overlap = mem_allocated_overlap + mem_intsize*size(P%angmom)&
     & + mem_intsize*size(P%indexAng1) + mem_intsize*size(P%indexAng2) &
     & + mem_intsize*size(P%nOrbitals) + mem_intsize*size(P%nOrbComp)&
     & + mem_intsize*size(P%nContracted)
!$OMP END CRITICAL (memory)

P%ETUVisAlloc = .FALSE.
P%FTUVisAlloc = .FALSE.
IF (allocFTUV) THEN
   CALL MEM_ALLOC(P%FTUV,nPrim1,nTUV,ndmat)
   P%FTUVisAlloc = .TRUE.
!$OMP CRITICAL (memory)
   mem_allocated_overlap = mem_allocated_overlap + mem_realsize*size(P%FTUV)
!$OMP END CRITICAL (memory)
ELSE
   IF(.NOT. hermitesingle)THEN
      IF (setETUVoutside) THEN
         CALL MEM_ALLOC(P%ETUVindex,nAng)
         CALL MEM_ALLOC(P%ETUV,MAXPRIMTUVIJK) !dimensions
         P%ETUVisAlloc = .TRUE.
!$OMP CRITICAL (memory)
         mem_allocated_overlap = mem_allocated_overlap + mem_realsize*size(P%ETUV)&
              & + mem_intsize*size(P%ETUVindex)
!$OMP END CRITICAL (memory)
      ENDIF
   ENDIF
ENDIF

!$OMP CRITICAL (memory)
max_mem_used_overlap = MAX(max_mem_used_overlap, mem_allocated_overlap)
!$OMP END CRITICAL (memory)

END SUBROUTINE ALLOCATE_OVERLAP

!> \brief allocate the orbital
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Orb The orbital, to allocate
!> \param maxprim the maximum number of primitives
!> \param npass the number of passes 
!> \param maxangmom the maximum number of angular momentum
!> \param maxA the maximum number of spherical transformation matrices
!> \param FTUVorb flag to determine if this is a FTUV orbital
SUBROUTINE ALLOCATE_ORBITAL(Orb,maxprim,npass,maxangmom,maxA,FTUVorb)
use memory_handling
IMPLICIT NONE
TYPE(Orbital) :: Orb
Integer :: maxprim,maxangmom,maxA,npass
logical :: FTUVorb

CALL MEM_ALLOC(Orb%angmom,maxangmom)
CALL MEM_ALLOC(Orb%nContracted,maxangmom)
CALL MEM_ALLOC(Orb%nOrbComp,maxangmom)
CALL MEM_ALLOC(Orb%startOrbital,maxangmom,npass)
CALL MEM_ALLOC(Orb%startLocOrb,maxangmom)
CALL MEM_ALLOC(Orb%atom,npass)
CALL MEM_ALLOC(Orb%batch,npass)

!$OMP CRITICAL (memory)
mem_allocated_overlap = mem_allocated_overlap + mem_intsize*size(Orb%angmom)&
     & + mem_intsize*size(Orb%nContracted) +  mem_intsize*size(Orb%nOrbComp)&
     & + mem_intsize*size(Orb%startOrbital)+ mem_intsize*size(Orb%atom)& 
     & + mem_intsize*size(Orb%startLocOrb) + mem_intsize*size(Orb%batch)
!$OMP END CRITICAL (memory)

IF(.NOT.FTUVorb)THEN
 CALL MEM_ALLOC(Orb%exponents,maxprim)
 CALL MEM_ALLOC(Orb%startprimOrbital,maxangmom,npass)
 CALL MEM_ALLOC(Orb%nPrimOrbComp,maxangmom)
 CALL MEM_ALLOC(Orb%nOrbitals,maxangmom)
!$OMP CRITICAL (memory)
 mem_allocated_overlap = mem_allocated_overlap + mem_realsize*size(Orb%exponents)
 mem_allocated_overlap = mem_allocated_overlap + mem_intsize*size(Orb%startprimOrbital)&
      & + mem_intsize*size(Orb%nPrimOrbComp) + mem_intsize*size(Orb%nOrbitals)
!$OMP END CRITICAL (memory)
 NULLIFY(Orb%CC)
 ALLOCATE(Orb%CC(maxangmom))
 NULLIFY(Orb%SPH_MAT)
 ALLOCATE(Orb%SPH_MAT(maxA+1))
 Orb%FTUVorb = .FALSE.
ELSE
 Orb%FTUVorb = .TRUE.
ENDIF
!$OMP CRITICAL (memory)
max_mem_used_overlap = MAX(max_mem_used_overlap, mem_allocated_overlap)
!$OMP END CRITICAL (memory)

END SUBROUTINE ALLOCATE_ORBITAL

!> \brief initialise the already allocated overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P The overlap, to be initialised 
!> \param np the number of primitves functions
!> \param Input contains info about the integral evaluation requested
!> \param sharedTUV contains TUVindexing 
!> \param Integral storage of integrals and intermediates
!> \param ODB The OD-batch information
!> \param iElectron The electron (1 = 'LHS', 2 = 'RHS')
!> \param lupri The default print unit
!> \param iprint The print level (the higher the more information)
!> \param side The side ('LHS' or 'RHS') 
!> \param nPrimAsInput the number of primitive is either given as input or calc
!> \param SETORBITAL if the ORBITALS have ALREADY been SET
SUBROUTINE SET_INITIALIZED_OVERLAP(P,np,Input,sharedTUV,Integral,ODB,IELECTRON,LUPRI,IPRINT,SIDE,nPrimAsInput,SETORBITAL)
Implicit none
TYPE(IntegralInput),target :: Input
TYPE(IntegralItem)  :: Integral
TYPE(Overlap)       :: P
TYPE(ODBATCH)       :: ODB
Character*(*)       :: side
Integer             :: IELECTRON,LUPRI,IPRINT,nderiv,np
TYPE(TUVitem)       :: SharedTUV
LOGICAL             :: nPrimAsInput,SETORBITAL
!
Integer             :: idir,i1,i2,i12,start1,start2,end1,end2,orb1,orb2
integer             :: l1,l2,ijk1,ijk2,ijk,ijkcart,maxijk,na
Real(realk)         :: e1,e2,d2,maxGab,maxPrimGab,maxPrimGabElm,signP,DMATmax
!Character(len=80)   :: t1,t2
!Real(realk),pointer :: GAB(:,:)
!Real(realk),pointer :: primGAB(:,:)
type(lstensor),pointer :: DMAT2,GAB2,primGAB2
!Real(realk),pointer :: DMAT(:,:,:)
LOGICAL             :: LHS, DMATscreen,useFTUV,doscreen,screen
Integer             :: ndmat,idmat,atomA,atomB,batchA,batchB,sA,dimA,sAA,sB,dimB
Integer             :: Dindex,iSA,iSB,elms,Gindex
real(realk),pointer :: primgabmat(:,:)
!
DMATscreen = .FALSE.
NULLIFY(DMAT2)
NULLIFY(GAB2)
NULLIFY(primGAB2)
SELECT CASE(SIDE)
CASE('LHS')
   IF(INPUT%PS_SCREEN)THEN
      primGAB2 => INPUT%LST_pGAB_LHS
      maxPrimGabElm=INPUT%PS_MAXELM_RHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
      GAB2 => INPUT%LST_GAB_LHS
   ENDIF
   LHS=.TRUE.
   ndmat=1
   NDERIV=INPUT%NDERIVP
   P%derivOrder = INPUT%derOrderP
CASE('RHS')
   IF(INPUT%PS_SCREEN)THEN
      primGAB2 => INPUT%LST_pGAB_RHS
      maxPrimGabElm=INPUT%PS_MAXELM_LHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
      GAB2 => INPUT%LST_GAB_RHS
      ndmat = INPUT%NDMAT_RHS
      IF (Input%DO_JENGINE) THEN
        DMATscreen = .TRUE.
        DMAT2 => INPUT%LST_DRHS
      ENDIF
   ENDIF
   LHS=.FALSE.
   NDERIV=INPUT%NDERIVQ
   P%derivOrder = INPUT%derOrderQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_INITIALIZED_OVERLAP =',SIDE
   CALL LSQUIT('Wrong case in SET_INITIALIZED_OVERLAP',lupri)
END SELECT

IF(SETORBITAL)THEN
 CALL SET_ORBITAL(P%orbital1,ODB%AO(1)%p,integral,LUPRI)
 CALL SET_ORBITAL(P%orbital2,ODB%AO(2)%p,integral,LUPRI)
ELSE
 !ORBITALS ALREADY SET BUT THE ORBITALS HAVE A NEW CENTER AND OTHER indices
 na   = ODB%AO(1)%p%nAngmom
 P%orbital1%center                    = ODB%AO(1)%p%center
 P%orbital1%startOrbital(1:na,1)      = ODB%AO(1)%p%startOrbital(1:na)
 P%orbital1%atom(1)                   = ODB%AO(1)%p%atom
 P%orbital1%batch(1)                  = ODB%AO(1)%p%batch
 P%orbital1%startprimOrbital(1:na,1)  = ODB%AO(1)%p%startprimOrbital(1:na)
 na   = ODB%AO(2)%p%nAngmom
 P%orbital2%center                    = ODB%AO(2)%p%center
 P%orbital2%startOrbital(1:na,1)      = ODB%AO(2)%p%startOrbital(1:na)
 P%orbital2%atom(1)                   = ODB%AO(2)%p%atom
 P%orbital2%batch(1)                  = ODB%AO(2)%p%batch
 P%orbital2%startprimOrbital(1:na,1)  = ODB%AO(2)%p%startprimOrbital(1:na)
ENDIF

 P%type_Empty = P%orbital1%type_Empty .AND. P%orbital2%type_Empty 
 IF(P%type_Empty)THEN
    IF(Input%operator(1:6) .EQ. 'Mulmom' .AND. (.NOT.LHS) )THEN
       ! this is done in order to use screening on multipole moments used for FMM
       P%maxGab = 1.D0
    ENDIF
 ENDIF
 !ONLY ONE OF THESE TYPES MUST BE TRUE
 P%type_Hermite_single = (P%orbital1%type_Hermite .AND. P%orbital2%type_Empty)&
      &.OR.(P%orbital1%type_Empty .AND. P%orbital2%type_Hermite)
 P%type_Hermite = P%orbital1%type_Hermite .AND. P%orbital2%type_Hermite
 P%type_Cartesian_single = (P%orbital1%type_Cartesian .AND. P%orbital2%type_Empty)&
      &.OR.(P%orbital1%type_Empty .AND. P%orbital2%type_Cartesian)
 P%type_Cartesian = P%orbital1%type_Cartesian .AND. P%orbital2%type_Cartesian
 P%type_Nucleus = (P%orbital1%type_Nucleus .AND. P%orbital2%type_Empty)&
      &.OR.(P%orbital1%type_Empty .AND. P%orbital2%type_Nucleus)
 P%type_FTUV = .FALSE.
 P%sphericalEcoeff = Input%sphericalEcoeff.AND.(.NOT.P%type_hermite_single)
 useFTUV = INPUT%DO_JENGINE .AND. IELECTRON.EQ.2
 P%sphericalEcoeff = P%sphericalEcoeff .AND.(.NOT.useFTUV).AND.&
     &               (.NOT.INPUT%operator(1:7) .EQ. 'Kinetic' .AND. LHS) 

 P%single = P%type_Hermite_single.OR.P%type_Cartesian_single.OR.P%type_Nucleus
 CALL getDerivComp(P%nDerivComp,P%derivOrder,P%type_empty,P%single)

 P%ODcenter = ODB%ODcenter
 P%ODextent = ODB%ODextent
 P%sameAO = ODB%sameAO
 P%ETUVisSet = .FALSE.
IF(.NOT.nPrimAsInput)THEN
 CALL GET_NPRIMITIVES(np,P,Input,Side,lupri)
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
 doscreen = ((INPUT%PS_SCREEN).AND.(.NOT.P%type_Empty))
 IF(doscreen) THEN
  call mem_alloc(primgabmat,P%orbital1%nPrimitives,P%orbital2%nPrimitives)
  AtomA = P%orbital1%atom(1)
  AtomB = P%orbital2%atom(1)
  BatchA = P%orbital1%batch(1)
  BatchB = P%orbital2%batch(1)
  Gindex = primGAB2%INDEX(P%orbital1%atom(1),P%orbital2%atom(1),1,1)
  DO i1=1,P%orbital1%nPrimitives
   DO i2=1,P%orbital2%nPrimitives
    primgabmat(i1,i2) = &
         &  primGAB2%LSAO(Gindex)%BATCH(BatchA,BatchB,1,1)%elms(i1 + (i2-1)*P%orbital1%nPrimitives)
   ENDDO
  ENDDO
 ENDIF
 i12 = 0
! maxPrimGab = 0.d0
 DO i1=1,P%orbital1%nPrimitives
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nPrimitives
     IF (doscreen) THEN
!       maxPrimGab = MAX(maxPrimGab,primgabmat(i1,i2))
       screen = primgabmat(i1,i2) .LT. INPUT%PS_THRESHOLD/maxPrimGabElm
     ELSE
       screen = .FALSE.
     ENDIF
!    Primitive screening
     IF (.NOT.screen) THEN
       i12 = i12 + 1
       e1  = P%orbital1%exponents(i1)
       e2  = P%orbital2%exponents(i2)
       P%exponents(i12)        = e1 + e2
       IF (P%exponents(i12) .NE. 0.0d0) THEN
          P%reducedExponents(i12) = e1*e2/(e1+e2)
          P%center(1,i12)         = (e1*P%orbital1%center(1) + e2*P%orbital2%center(1))/(e1+e2)
          P%center(2,i12)         = (e1*P%orbital1%center(2) + e2*P%orbital2%center(2))/(e1+e2)
          P%center(3,i12)         = (e1*P%orbital1%center(3) + e2*P%orbital2%center(3))/(e1+e2)
       ELSEIF(Input%operator(1:6) .EQ. 'Nucrep')THEN
          P%reducedExponents(i12) = e1
          P%center(1,i12)         = P%orbital1%center(1)
          P%center(2,i12)         = P%orbital1%center(2)
          P%center(3,i12)         = P%orbital1%center(3)
       ELSE
         P%reducedExponents(i12) = 0.0d0
         P%center(1,i12)         = 0.0d0
         P%center(2,i12)         = 0.0d0
         P%center(3,i12)         = 0.0d0
       ENDIF
       IF ( P%single) THEN 
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
   WRITE(LUPRI,'(1X,A,2I5,L2)') 'Error in set_INITIALIZED_overlap. Mismatch in number of primitives, (np,i12)=',&
     &  np,i12,nPrimAsInput
   CALL LSQUIT('Error: Mismatch between get_nPrimitives and i12 in set_INITIALIZED_overlap',lupri)
 ENDIF

IF(INPUT%PS_SCREEN)THEN
!   P%maxPrimGab = maxPrimGab
   NULLIFY(primGAB2)
ENDIF
IF(doscreen)call mem_dealloc(primgabmat)

!CAREFUL
!Initialize to twice some maximal value + 1
 P%minAngmom   = 99
!CAREFUL
 P%maxAngmom = 0
 P%totOrbitals = 0
 i12 = 0

 IF (Input%CS_SCREEN) THEN
    AtomA = P%orbital1%atom(1)
    BatchA = P%orbital1%batch(1)
    AtomB = P%orbital2%atom(1)
    BatchB = P%orbital2%batch(1)
    Gindex = GAB2%INDEX(atomA,atomB,1,1)
    IF (DMATscreen) Dindex = DMAT2%INDEX(atomA,atomB,1,1)
    dimA = P%orbital1%nOrbComp(1)*P%orbital1%nContracted(1)
    DO SAA=2,P%orbital1%nAngmom
       dimA = dimA+P%orbital1%nOrbComp(SAA)*P%orbital1%nContracted(SAA)
    ENDDO
    dimB = P%orbital2%nOrbComp(1)*P%orbital2%nContracted(1)
    DO SAA=2,P%orbital2%nAngmom
       dimB = dimB+P%orbital2%nOrbComp(SAA)*P%orbital2%nContracted(SAA)
    ENDDO
ENDIF

 maxGab = 0.0d0
 maxijk = 0
 DO i1=1,P%orbital1%nAngmom
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nAngmom
     i12  = i12 + 1
     l1   = P%orbital1%angmom(i1)
     l2   = P%orbital2%angmom(i2)
     CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,P%sphericalEcoeff,P%derivOrder,P%single)
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
          sA = 0
          DO SAA=2,i1
             sA = sA+P%orbital1%nOrbComp(SAA-1)*P%orbital1%nContracted(SAA-1)
          ENDDO
          sB = 0
          DO SAA=2,i2
             sB = sB+P%orbital2%nOrbComp(SAA-1)*P%orbital2%nContracted(SAA-1)
          ENDDO
          iSA = (orb1-start1)+1+sA
          iSB = (orb2-start2)+1+sB            
          IF (DMATscreen) THEN
              DMATmax = 0.0d0
              DO idmat=1,ndmat
                 elms = iSA + (iSB-1)*dimA + (idmat-1)*dimA*dimB
                 DMATmax = max(DMATmax,abs(DMAT2%LSAO(Dindex)%BATCH(batchA,batchB,1,1)%elms(elms)))
              ENDDO
              elms = iSA + (iSB-1)*dimA
              maxGab = MAX(maxGab,GAB2%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%elms(elms)*DMATmax)
          ELSE
             elms = iSA + (iSB-1)*dimA
             maxGab = MAX(maxGab,GAB2%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%elms(elms))
          ENDIF
         ENDDO
       ENDDO
     ENDIF
   ENDDO
 ENDDO
 IF (Input%CS_SCREEN) THEN
   P%maxGab = maxGab
   NULLIFY(GAB2)
 ENDIF

 P%startAngmom = 0
 IF (P%type_hermite_single) P%startAngmom = P%minAngmom + P%derivOrder

 IF( INPUT%operator(1:7) .EQ. 'Kinetic' .AND. LHS) THEN
    P%startAngmom = P%startAngmom + 2
    P%endAngmom   = P%maxAngmom + 2 + P%derivOrder
 ELSE
    P%endAngmom   = P%maxAngmom + P%derivOrder
 ENDIF  

 P%nTUV = (P%endAngmom+1)*(P%endAngmom+2)*(P%endAngmom+3)/6 - &
    &     P%startAngmom*(P%startAngmom+1)*(P%startAngmom+2)/6

 IF (IELECTRON.EQ.1) THEN
   signP=1.0d0
 ELSE
   signP=-1.0d0
 ENDIF

IF (INPUT%DO_JENGINE .AND. IELECTRON.EQ.2) THEN
  CALL SET_FTUV(P,INPUT,Integral,LUPRI,IPRINT)
! Change Overlap P to confirm with FTUV
  P%type_FTUV = .TRUE.
  P%PassType = 1
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
  P%totOrbitals    = ndmat
  P%nOrbComp(1)    = 1
  P%nContracted(1) = 1
ELSE
!GET E COEFFICIENTS
 IF (.NOT.P%type_hermite_single) THEN
   IF (INPUT%setETUVoutside) THEN
      IF(LHS)THEN
         CALL AttachETUVtoOverlap(P,signP,SharedTUV,.TRUE.,LUPRI,IPRINT)
      ELSE
         CALL AttachETUVtoOverlap(P,signP,SharedTUV,Input%orderAngPrim,LUPRI,IPRINT)
      ENDIF
   ENDIF
 ENDIF
ENDIF

END SUBROUTINE SET_INITIALIZED_OVERLAP

!> \brief set the orbitals when passes are used 
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap, to be built
!> \param P the overlap
!> \param maxpasses the maximum number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SetPassOrbitals(PassP,P,maxPasses,LUPRI,IPRINT)
Implicit none
TYPE(Overlap)       :: PassP,P
Integer             :: maxPasses,LUPRI,IPRINT
!
Integer :: np,na,mc,ma,mp,ia
!
integer :: i,startPrim,endOrbital,ip

! Orbital 1
np = P%orbital1%nPrimitives
na = P%orbital1%nAngmom
mc = P%orbital1%maxContracted
ma = P%orbital1%maxAngmom

PassP%orbital1%TYPE_Empty = P%orbital1%type_Empty
PassP%orbital1%TYPE_Hermite = P%orbital1%type_Hermite
PassP%orbital1%TYPE_Cartesian = P%orbital1%type_Cartesian
PassP%orbital1%TYPE_Nucleus = P%orbital1%type_Nucleus

PassP%orbital1%spherical     = P%orbital1%spherical     
PassP%orbital1%center        = P%orbital1%center        
PassP%orbital1%nPrimitives   = np
PassP%orbital1%maxContracted = mc
PassP%orbital1%maxAngmom     = ma
PassP%orbital1%nAngmom       = na
PassP%orbital1%totOrbitals   = P%orbital1%totOrbitals
PassP%orbital1%CCidentifier  = P%orbital1%CCidentifier  
!replace with dcopy
PassP%orbital1%exponents(1:np) = P%orbital1%exponents(1:np)
PassP%orbital1%angmom(1:na)       = P%orbital1%angmom(1:na)
PassP%orbital1%nContracted(1:na)  = P%orbital1%nContracted(1:na)
PassP%orbital1%nOrbComp(1:na)     = P%orbital1%nOrbComp(1:na)
PassP%orbital1%nPrimOrbComp(1:na) = P%orbital1%nPrimOrbComp(1:na)
PassP%orbital1%nOrbitals(1:na)    = P%orbital1%nOrbitals(1:na)
PassP%orbital1%startLocOrb(1:na)  = P%orbital1%startLocOrb(1:na)
!maybe replace with do loop
!   PassP%orbital1%CC(1:np,1:mc,1:na) = P%orbital1%CC(1:np,1:mc,1:na)
DO ia = 1,na
   PassP%orbital1%CC(ia)%p => P%orbital1%CC(ia)%p
ENDDO
PassP%orbital1%SPH_MAT(1:ma+1) = P%orbital1%SPH_MAT(1:ma+1)

! Orbital 2
np = P%orbital2%nPrimitives
na = P%orbital2%nAngmom
mc = P%orbital2%maxContracted
ma = P%orbital2%maxAngmom

PassP%orbital2%TYPE_Empty = P%orbital2%type_Empty
PassP%orbital2%TYPE_Hermite = P%orbital2%type_Hermite
PassP%orbital2%TYPE_Cartesian = P%orbital2%type_Cartesian
PassP%orbital2%TYPE_Nucleus = P%orbital2%type_Nucleus

PassP%orbital2%spherical     = P%orbital2%spherical     
PassP%orbital2%center        = P%orbital2%center        
PassP%orbital2%nPrimitives   = np
PassP%orbital2%maxContracted = mc
PassP%orbital2%maxAngmom     = ma
PassP%orbital2%nAngmom       = na
PassP%orbital2%totOrbitals   = P%orbital2%totOrbitals
PassP%orbital2%CCidentifier  = P%orbital2%CCidentifier  
PassP%orbital2%exponents(1:np) = P%orbital2%exponents(1:np)
PassP%orbital2%angmom(1:na)       = P%orbital2%angmom(1:na)
PassP%orbital2%nContracted(1:na)  = P%orbital2%nContracted(1:na)
PassP%orbital2%nOrbComp(1:na)     = P%orbital2%nOrbComp(1:na)
PassP%orbital2%nPrimOrbComp(1:na) = P%orbital2%nPrimOrbComp(1:na)
PassP%orbital2%nOrbitals(1:na)    = P%orbital2%nOrbitals(1:na)
PassP%orbital2%startLocOrb(1:na)  = P%orbital2%startLocOrb(1:na)
DO ia = 1,na
   PassP%orbital2%CC(ia)%p => P%orbital2%CC(ia)%p
ENDDO

PassP%orbital2%SPH_MAT(1:ma+1) = P%orbital2%SPH_MAT(1:ma+1)

!Common Overlap
mp = P%nPrimitives*maxPasses
np = P%nPrimitives
na = P%nAngmom

PassP%type_Empty            = P%type_Empty
PassP%type_Hermite          = P%type_Hermite
PassP%type_Hermite_single   = P%type_Hermite_single
PassP%type_Cartesian        = P%type_Cartesian
PassP%type_Cartesian_single = P%type_Cartesian_single
PassP%type_Nucleus          = P%type_Nucleus
PassP%type_FTUV             = P%type_FTUV

PassP%single                = P%single
PassP%nDerivComp            = P%nDerivComp

PassP%sphericalEcoeff   = P%sphericalEcoeff
PassP%sameAO            = P%sameAO
PassP%ETUVisSet         = P%ETUVisSet
PassP%nPrimitives       = P%nPrimitives
PassP%maxContracted     = P%maxContracted
PassP%squaredDistance   = P%squaredDistance
!PassP%maxPrimGab        = P%maxPrimGab 
PassP%maxGab            = P%maxGab 
PassP%nAngmom           = P%nAngmom      
PassP%angmom(1:na)      = P%angmom(1:na)
PassP%indexAng1(1:na)   = P%indexAng1(1:na)
PassP%indexAng2(1:na)   = P%indexAng2(1:na)
PassP%nOrbComp(1:na)    = P%nOrbComp(1:na)
PassP%nContracted(1:na) = P%nContracted(1:na)
PassP%minAngmom         = P%minAngmom
PassP%maxAngmom         = P%maxAngmom
PassP%startAngmom       = P%startAngmom
PassP%derivOrder        = P%derivOrder
PassP%endAngmom         = P%endAngmom
PassP%nTUV              = P%nTUV

startPrim = 0
endOrbital = nP-1
DO I = 1,maxpasses
   DO IP = 1,np
      PassP%exponents(startPrim+IP) = P%exponents(IP)
   ENDDO
   startPrim = startPrim+np
   endOrbital = endOrbital+np
ENDDO

END SUBROUTINE SetPassOrbitals

!> \brief add the overlap P to an overlap, which we call PassP
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap, to be built
!> \param P the overlap
!> \param numPass the current number of passes in PassP 
!> \param maxpasses the maximum number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddOverlapToPass(PassP,P,numPass,maxPasses,LUPRI,IPRINT)
Implicit none
TYPE(Overlap)       :: PassP,P
Integer             :: numPass,maxPasses,LUPRI,IPRINT
!
Integer :: np,na
Integer :: iAngmom,length,indexETUV,passIndexETUV
!
integer :: i,startPassIndex,startPrim,endOrbital

! Orbital indexes
na = P%orbital1%nAngmom
!PassP%orbital1%startPrimOrbital(1:na,numPass) = P%orbital1%startPrimOrbital(1:na,1)
!PassP%orbital2%startPrimOrbital(1:na,numPass) = P%orbital2%startPrimOrbital(1:na,1)
PassP%orbital1%startOrbital(1:na,numPass)     = P%orbital1%startOrbital(1:na,1)
na = P%orbital2%nAngmom
PassP%orbital2%startOrbital(1:na,numPass)     = P%orbital2%startOrbital(1:na,1)
PassP%orbital1%atom(numPass)     = P%orbital1%atom(1)
PassP%orbital2%atom(numPass)     = P%orbital2%atom(1)
PassP%orbital1%batch(numPass)     = P%orbital1%batch(1)
PassP%orbital2%batch(numPass)     = P%orbital2%batch(1)

np=P%nPrimitives
startPrim  = 1+(numPass-1)*np
endOrbital = numPass*nP
PassP%center(1:3,startPrim:endOrbital)       = P%center(:,1:np)
!PassP%exponents(startPrim:endOrbital)        = P%exponents(1:np)
!PassP%reducedExponents(startPrim:endOrbital) = P%reducedExponents(1:np)

IF (PassP%ETUVisSet) THEN
   DO iAngmom=1,P%nAngmom
      ! Set up indeces and dimensions
      indexETUV = P%ETUVindex(iAngmom)
      passIndexETUV = (indexETUV-1)*maxPasses+1
      PassP%ETUVindex(iAngmom) = passIndexETUV
      length = P%lenETUV(iAngmom)
      startPassIndex = passIndexETUV + (numPass-1)*length
      !The Ordering of P%ETUV is (tuvP,ijkP,nPrimP,iAngmom)
      !The Ordering of PassP%ETUV is (tuvP,ijkP,nPrimP,nPassP,iAngmom)
      DO I=0,length-1
         PassP%ETUV(startPassIndex+I) = P%ETUV(indexETUV+I)
      ENDDO
!      PassP%ETUV(startPassIndex:startPassIndex + length - 1) = P%ETUV(indexETUV:indexETUV+length - 1)
   ENDDO
ELSE
   PassP%iPrim1(startPrim:endOrbital)           = P%iprim1(1:np)
   PassP%iPrim2(startPrim:endOrbital)           = P%iprim2(1:np)
   PassP%preExpFac(startPrim:endOrbital)        = P%preExpFac(1:np)
   PassP%distance12(:,numpass)      = P%distance12(:,1)
ENDIF

END SUBROUTINE AddOverlapToPass

!> \brief Set the Final Passp overlap values
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap, to be built
!> \param P the overlap
!> \param numPass the current number of passes in PassP 
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE FinalizePass(PassP,P,numPass,LUPRI,IPRINT)
Implicit none
TYPE(Overlap)       :: PassP,P
Integer             :: numPass,LUPRI,IPRINT
!
integer :: na,np

IF (PassP%ETUVisSet) THEN
   np=P%nPrimitives
   PassP%iPrim1(1:np)           = P%iprim1(1:np)
   PassP%iPrim2(1:np)           = P%iprim2(1:np)
ENDIF

PassP%orbital1%nPasses  = numPass 
PassP%orbital2%nPasses  = numPass 
PassP%nPasses = numPass                                                             
PassP%totOrbitals = P%totOrbitals*numPass                                           
na = P%nAngmom
PassP%nOrbitals(1:na)   = P%nOrbitals(1:na)*numPass   

END SUBROUTINE FinalizePass

!> \brief attach Ecoeff to overlap 
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P the overlap
!> \param signP the sign (1 for P, -1 for Q)
!> \param TUV The TUV indeces
!> \param orderAP True if the angular components come before the primitives
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AttachETUVtoOverlap(P,signP,TUV,orderAp,LUPRI,IPRINT)
Implicit none
TYPE(Overlap) :: P
TYPE(TUVitem) :: TUV
Integer       :: LUPRI,IPRINT
Real(realk)   :: signP
logical       :: orderAP
!
Integer :: iA1,iA2,l1,l2,l,ijk1,ijk2,ijk,lm,nTUV,indexETUV,iAngmom
integer :: LENGTH

indexETUV = 1
DO iAngmom=1,P%nAngmom
  P%ETUVindex(iAngmom) = indexETUV
  iA1 = P%indexAng1(iangmom)
  iA2 = P%indexAng2(iangmom)
  l1 = P%orbital1%angmom(iA1)
  l2 = P%orbital2%angmom(iA2)
  l  = l1+l2+P%derivOrder
  CALL GET_IJK(l1,l2,ijk1,ijk2,lm,ijk,P%sphericalEcoeff,P%derivOrder,P%single)
  nTUV = (l+1)*(l+2)*(l+3)/6
  LENGTH=lm*nTUV*P%nPrimitives
!  LENGTH = P%lenETUV(iAngmom) 
  CALL BuildEcoeffTensor(TUV,P,signP,P%ETUV(indexETUV:indexETUV+LENGTH-1),lm,ijk,nTUV,P%nPrimitives,&
     &                   orderAP,iAngmom,1,LUPRI,IPRINT)
  indexETUV = indexETUV + LENGTH
ENDDO
P%ETUVisSet = .TRUE.

END SUBROUTINE AttachETUVtoOverlap

!> \brief set up the FTUV batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P the overlap
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SET_FTUV(P,Input,Integral,LUPRI,IPRINT)
use memory_handling
Implicit none
TYPE(Overlap)       :: P
TYPE(IntegralInput) :: Input
TYPE(IntegralItem)  :: Integral
Integer             :: LUPRI,IPRINT
!
Logical             :: Sph1,Sph2
Integer             :: iAngmom,l1,l2,l,ijk1,ijk2,ijk,nTUV,lm1,lm2,lm
Integer             :: iPrim,iA1,iA2,nCont,nC1,nC2,iAng1,iAng2
Integer             :: idmat,start1,start2,i1,i2,iC1,iC2,indlm,iTUV
Integer             :: startAng, tuvOff
Real(realk),pointer :: Ecoeffs(:,:,:),Spherical(:,:),SpherEcoeff(:,:,:)
Real(realk),pointer :: CC(:,:,:),ContractedDmat(:,:,:)
Real(realk),parameter :: DM1=-1.0d0
real(realk)          :: dtemp
Integer             :: dimC,dimD,SAA,atomC,atomD,batchC,batchD,Dcdindex
Integer             :: Ddcindex,sC,sD,iSC,iSD,elms


CALL LS_DZERO(P%FTUV,P%nTUV*P%nPrimitives*Input%NDMAT_RHS)

dimC = P%orbital1%nOrbComp(1)*P%orbital1%nContracted(1)
DO SAA=2,P%orbital1%nAngmom
   dimC = dimC+P%orbital1%nOrbComp(SAA)*P%orbital1%nContracted(SAA)
ENDDO
dimD = P%orbital2%nOrbComp(1)*P%orbital2%nContracted(1)
DO SAA=2,P%orbital2%nAngmom
   dimD = dimD+P%orbital2%nOrbComp(SAA)*P%orbital2%nContracted(SAA)
ENDDO

DO iAngmom = 1,P%nAngmom
  iA1 = P%indexAng1(iangmom)
  iA2 = P%indexAng2(iangmom)
  l1 = P%orbital1%angmom(iA1)
  l2 = P%orbital2%angmom(iA2)
  l  = l1+l2
  ijk1 = (l1+1)*(l1+2)/2
  ijk2 = (l2+1)*(l2+2)/2
  ijk  = ijk1*ijk2
  startAng = P%startAngmom
  tuvOff   = startAng*(startAng+1)*(startAng+2)/6
  nTUV = (l+1)*(l+2)*(l+3)/6 - tuvOff
  call mem_alloc(Ecoeffs,nTUV,ijk,P%nPrimitives)
  !NULLIFY(Ecoeffs)
  !ALLOCATE(Ecoeffs(nTUV,ijk,P%nPrimitives))
  CALL BuildEcoeffTensor1_AP(Integral%TUV,P,DM1,Ecoeffs,nTUV,ijk,P%nPrimitives,iAngmom,1,LUPRI,IPRINT)

! Spherical transformation of E-coefficients
  Sph1 = P%orbital1%spherical.AND.(l1.GT.1)
  Sph2 = P%orbital2%spherical.AND.(l2.GT.1)
  lm1 = ijk1
  IF (Sph1) lm1 = 2*l1+1
  lm2 = ijk2
  IF (Sph2) lm2 = 2*l2+1
  lm  = lm1*lm2
  IF (Sph1.OR.Sph2) THEN
    call mem_alloc(Spherical,ijk,lm)
    !NULLIFY(Spherical)
    !ALLOCATE(Spherical(ijk,lm))
    CALL Buildsphericaltransformation(Spherical,P%orbital1%SPH_MAT(l1+1)%p%elms,&
   &                                      P%orbital2%SPH_MAT(l2+1)%p%elms,&
   &                                      lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)
    call mem_alloc(SpherEcoeff,nTUV,lm,P%nPrimitives)
    !NULLIFY(SpherEcoeff)
    !ALLOCATE(SpherEcoeff(nTUV,lm,P%nPrimitives))
    DO iPrim=1,P%nPrimitives


      CALL DGEMM('N','N',nTUV,lm,ijk,1.0d0,Ecoeffs(1,1,iPrim),nTUV,Spherical,ijk,&
   &             0.0d0,SpherEcoeff(1,1,iPrim),nTUV)
    ENDDO
    call mem_dealloc(Ecoeffs)
    call mem_dealloc(Spherical)
    !DEALLOCATE(Ecoeffs)
    !DEALLOCATE(Spherical)
  ELSE
    NULLIFY(SpherEcoeff)
    SpherEcoeff => Ecoeffs

  ENDIF
! Set up contraction matrix
  iA1 = P%indexAng1(iangmom)
  iA2 = P%indexAng2(iangmom)
  nC1 = P%orbital1%nContracted(iA1)
  nC2 = P%orbital2%nContracted(iA2)
  nCont = nC1*nC2
  call mem_alloc(CC,P%nPrimitives,nC1,nC2)
  !NULLIFY(CC)
  !ALLOCATE(CC(P%nPrimitives,nC1,nC2))
  CALL ConstructContraction_AP(CC,P,P%nPrimitives,nC1,nC2,iA1,iA2,LUPRI,IPRINT)

! Make contraction of density-matrix and contraction coefficients
  call mem_alloc(ContractedDmat,P%nPrimitives,lm,Input%NDMAT_RHS)
  !NULLIFY(ContractedDmat)
  !ALLOCATE(ContractedDmat(P%nPrimitives,lm,Input%NDMAT_RHS))
  CALL LS_DZERO(ContractedDmat,lm*P%nPrimitives*Input%NDMAT_RHS)
  start1 = P%orbital1%startOrbital(iA1,1)
  start2 = P%orbital2%startOrbital(iA2,1)

  AtomC = P%orbital1%atom(1)
  BatchC = P%orbital1%batch(1)
  AtomD = P%orbital2%atom(1)
  BatchD = P%orbital2%batch(1)
  Dcdindex = INPUT%LST_DRHS%INDEX(atomC,atomD,1,1)
  IF(Input%sameRHSaos) Ddcindex = INPUT%LST_DRHS%INDEX(atomD,atomC,1,1)
  
  sC = 0
  DO SAA=2,iA1
     sC = sC+P%orbital1%nOrbComp(SAA-1)*P%orbital1%nContracted(SAA-1)
  ENDDO
  sD = 0
  DO SAA=2,iA2
     sD = sD+P%orbital2%nOrbComp(SAA-1)*P%orbital2%nContracted(SAA-1)
  ENDDO

  DO idmat=1,Input%nDMAT_RHS
    i1 = start1
    DO iC1=1,nC1
      DO iAng1=1,lm1
        i2 = start2
        DO iC2=1,nC2
          DO iAng2=1,lm2
            indlm = (iAng2-1)*lm1 + iAng1

            iSC = (i1-start1)+1+sC
            iSD = (i2-start2)+1+sD
            elms = iSC + (iSD-1)*dimC + (idmat-1)*dimC*dimD
            dtemp = INPUT%LST_DRHS%LSAO(Dcdindex)%BATCH(batchC,batchD,1,1)%elms(elms)
!            dtemp = Input%DMAT_RHS(i1,i2,idmat)
!            IF ((start1.NE.start2).AND.Input%sameRHSaos) &
!   &            dtemp = dtemp + Input%DMAT_RHS(i2,i1,idmat)
            IF ((start1.NE.start2).AND.Input%sameRHSaos) THEN
               elms = iSD + (iSC-1)*dimD + (idmat-1)*dimD*dimC
               dtemp = dtemp + INPUT%LST_DRHS%LSAO(Ddcindex)%BATCH(batchD,batchC,1,1)%elms(elms)
!               dtemp = dtemp + Input%DMAT_RHS(i2,i1,idmat)
            ENDIF
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
  call mem_dealloc(CC)
  !DEALLOCATE(CC)

! Set up FTUV's
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
  call mem_dealloc(ContractedDmat)
  !DEALLOCATE(ContractedDmat)
  IF (Sph1.OR.Sph2) THEN
    call mem_dealloc(SpherEcoeff)
    !DEALLOCATE(SpherEcoeff)
  ELSE
    NULLIFY(SpherEcoeff)
    call mem_dealloc(Ecoeffs)
    !DEALLOCATE(Ecoeffs)
  ENDIF
ENDDO

IF (IPRINT.GT.20) THEN
  CALL LSHEADER(LUPRI,'SET_FTUV')
  WRITE(LUPRI,'(1X,A,I5)') 'Number of primitives', P%nPrimitives
  WRITE(LUPRI,'(1X,A,I5)') "Number of TUV's     ", nTUV
  WRITE(LUPRI,'(1X,A,I5)') "Number of D-mat     ", Input%NDMAT_RHS
  DO idmat=1,Input%NDMAT_RHS
    DO iTUV=1,nTUV
      WRITE(LUPRI,'(1X,A,I5,A,I5)') 'iTUV =', iTUV, ' iDmat=',idmat
      WRITE(LUPRI,'(3X,5F12.8)') (P%FTUV(iPrim,iTUV,idmat),iPrim=1,P%nPrimitives)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE SET_FTUV

!> \brief Finds the number of significant product primitives
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param nPrimitives the number of nPrimitives for this overlap 
!> \param P the overlap
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param side LHS or RHS
!> \param LUPRI the logical unit number for the output file
SUBROUTINE GET_NPRIMITIVES(nPrimitives,P,Input,Side,lupri)
Implicit none
TYPE(IntegralInput),target :: Input
TYPE(Overlap)       :: P
Integer             :: nPrimitives,lupri
Character*(*)       :: Side
!
Integer              :: i1,i2,start2
Real(realk)          :: maxgab,maxelm
!Real(realk), pointer :: GAB(:,:)
type(lstensor),pointer :: primGAB2
integer :: atomA,atomB,batchA,batchB,Gindex

IF (Side.EQ.'LHS') THEN
  primGAB2 => INPUT%LST_pGAB_LHS
!  GAB => INPUT%pGAB_LHS
  MAXELM=INPUT%PS_MAXELM_RHS
ELSEIF (Side.EQ.'RHS') THEN
  primGAB2 => INPUT%LST_pGAB_RHS
!  GAB => INPUT%pGAB_RHS
  MAXELM=INPUT%PS_MAXELM_LHS
ELSE
  WRITE(*,'(1X,2A)') 'Error in GET_NPRIMITIVES. Side =',Side
  CALL LSQUIT('Programming error. Wrong Side in GET_NPRIMITIVES',-1)
ENDIF
nPrimitives = 0
IF(INPUT%PS_SCREEN.AND.(.NOT.P%type_Empty))THEN
  AtomA = P%orbital1%atom(1)
  AtomB = P%orbital2%atom(1)
  BatchA = P%orbital1%batch(1)
  BatchB = P%orbital2%batch(1)
  Gindex = primGAB2%INDEX(atomA,atomB,1,1)
  DO i1=1,P%orbital1%nPrimitives
     start2 = 1
     IF (P%sameAO) start2 = i1
     DO i2=start2,P%orbital2%nPrimitives
!        maxGab = getMaxPrimGab(primGAB2%LSAO(Gindex)%BATCH(BatchA,BatchB,1,1)%elms,i1,i2,P,lupri)
        maxGab = primGAB2%LSAO(Gindex)%BATCH(BatchA,BatchB,1,1)%elms(i1 + (i2-1)*P%orbital1%nPrimitives)
        IF( MAXGAB .GT. INPUT%PS_THRESHOLD/MAXELM)THEN
            nPrimitives = nPrimitives + 1
        ENDIF
      ENDDO
   ENDDO
ELSE
   i1 = P%orbital1%nPrimitives
   i2 = P%orbital2%nPrimitives
   nPrimitives = i1*i2
   IF (P%sameAO) nPrimitives = i1*(i1+1)/2
ENDIF

END SUBROUTINE GET_NPRIMITIVES

!> \brief set the values of the orbital
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Orb the orbital to be set
!> \param AOB the Atomic batch
!> \param INTEGRAL storage of integrals and intermediates
!> \param LUPRI the logical unit number for the output file
SUBROUTINE SET_ORBITAL(Orb,AOB,integral,LUPRI)
IMPLICIT NONE
TYPE(Orbital) :: Orb
TYPE(AOBATCH) :: AOB
TYPE(IntegralItem)  :: Integral
integer             :: lupri
!
Integer :: ia,na,np,i
!
 np   = AOB%nPrimitives
 na   = AOB%nAngmom

 Orb%TYPE_Empty = AOB%type_Empty
 Orb%TYPE_Hermite = AOB%type_Hermite
 Orb%TYPE_Cartesian = AOB%type_Cartesian
 Orb%TYPE_Nucleus = AOB%type_Nucleus

 Orb%spherical     = AOB%spherical
 Orb%center        = AOB%center
 Orb%nPrimitives   = np 
 Orb%nPasses       = 1 
 Orb%maxContracted = AOB%maxContracted
 Orb%maxAngmom     = AOB%maxAngmom
 Orb%nAngmom       = na
 Orb%CCidentifier  = AOB%CCindex(1)
 Orb%exponents(1:np) = AOB%pExponents%elms(1:np)
 Orb%angmom(1:na) = AOB%angmom(1:na)
 Orb%nContracted(1:na)  = AOB%nContracted(1:na)
 Orb%startOrbital(1:na,1) = AOB%startOrbital(1:na)
 Orb%atom(1) = AOB%atom
 Orb%batch(1) = AOB%batch
 Orb%startprimOrbital(1:na,1) = AOB%startprimOrbital(1:na)
 Orb%nOrbComp(1:na) = AOB%nOrbComp(1:na)
 Orb%nPrimOrbComp(1:na) = AOB%nPrimOrbComp(1:na)
 Orb%nOrbitals(1:na) = AOB%nOrbitals(1:na)
 Orb%totOrbitals   = 0 
 DO ia=1,na
    Orb%startLocOrb(ia) = Orb%totOrbitals + 1
    Orb%totOrbitals = Orb%totOrbitals + Orb%nOrbitals(ia)
 ENDDO
 DO ia=1,na
    Orb%CC(ia)%p => AOB%pCC(ia)%p
 ENDDO
 DO I=1,Orb%maxAngmom+1
    Orb%SPH_MAT(I)%p => integral%TUV%SPH_MAT(I)
 ENDDO

END SUBROUTINE SET_ORBITAL

!> \brief print the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P the overlap distribution
!> \param IUNIT the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param side LHS or RHS
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
ELSE
   WRITE(IUNIT,'(1X,A)') ''
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') '***   OVERLAP         ***'
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') ''
ENDIF
IF(.NOT.P%type_FTUV)THEN
   WRITE(IUNIT,'(3X,A)') '-------------- orbital 1 --------------'
   CALL PRINT_ORBITAL(P%orbital1,P%nPasses,IUNIT,IPRINT)
   WRITE(IUNIT,'(3X,A)') '-------------- orbital 2 --------------'
   CALL PRINT_ORBITAL(P%orbital2,P%nPasses,IUNIT,IPRINT)
ENDIF
 WRITE(IUNIT,'(3X,A)') '--------- overlap specifics -----------'
 WRITE(IUNIT,'(5X,A,L1)')  'type empty            =', P%type_Empty
 WRITE(IUNIT,'(5X,A,L1)')  'type cartesian-signle =', P%type_Cartesian_single
 WRITE(IUNIT,'(5X,A,L1)')  'type cartesian        =', P%type_Cartesian
 WRITE(IUNIT,'(5X,A,L1)')  'type hermite          =', P%type_hermite
 WRITE(IUNIT,'(5X,A,L1)')  'type hermite_single   =', P%type_hermite_single
 WRITE(IUNIT,'(5X,A,L1)')  'type ftuv             =', P%type_FTUV
 WRITE(IUNIT,'(5X,A,L1)')  'type Nucleus          =', P%type_Nucleus
 WRITE(IUNIT,'(5X,A,L1)')  'Single                =', P%single
 WRITE(IUNIT,'(5X,A,L1)')  'Spherical-Ecoeff      =', P%sphericalEcoeff

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
 WRITE(IUNIT,'(5X,A,I3)') 'derivOrder       =', P%derivOrder
 WRITE(IUNIT,'(5X,A,I3)') 'nDerivComp       =', P%nDerivComp
IF(.NOT. P%type_FTUV)THEN
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
IF(.NOT.P%type_FTUV)THEN
   DO ipass=1,P%nPasses
      WRITE(IUNIT,'(3X,A,I5)') '*** Pass number ',ipass
      WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------------'
      WRITE(IUNIT,'(3X,A)')    '  iPrim   Px      Py      Pz        p         mu    pre   iprim1 iprim2  '
      ip = (ipass-1)*P%nPrimitives
      DO i=1+ip,P%nPrimitives+ip
         WRITE(IUNIT,'(5X,I4,3F8.4,2F10.3,1F5.1,2I7)') i, P%center(1,i), P%center(2,i), P%center(3,i), &
              &  P%exponents(i), P%reducedExponents(i), P%preExpFac(i), P%iprim1(i), P%iprim2(i)
      ENDDO
   ENDDO
ELSE
   DO ipass=1,P%nPasses
      WRITE(IUNIT,'(3X,A,I5)') '*** Pass number ',ipass
      WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------------'
      WRITE(IUNIT,'(3X,A)')    '  iPrim   Px      Py      Pz        p         mu    pre'
      ip = (ipass-1)*P%nPrimitives
      DO i=1+ip,P%nPrimitives+ip
         WRITE(IUNIT,'(5X,I4,3F8.4,2F10.3,1F5.1)') i, P%center(1,i), P%center(2,i), P%center(3,i), &
              &  P%exponents(i), P%reducedExponents(i), P%preExpFac(i)
      ENDDO
   ENDDO
ENDIF
WRITE(IUNIT,'(3X,A)')    '------------------------------------------------------------------------'
END SUBROUTINE PRINT_OVERLAP

!> \brief print the orbital
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Orb the orbital to be printet
!> \param nPasses the number of passes
!> \param IUNIT the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE PRINT_ORBITAL(Orb,nPasses,IUNIT,IPRINT)
TYPE(Orbital),intent(IN) :: Orb
Integer,intent(IN)       :: IUNIT,IPRINT,nPasses
!
Integer :: I,J,K,firstAtom,iCent
Logical :: printCenter

 WRITE(IUNIT,'(5X,A17,1L2)')'Type Empty     = ', Orb%type_Empty
 WRITE(IUNIT,'(5X,A17,1L2)')'Type Cartesian = ', Orb%type_Cartesian
 WRITE(IUNIT,'(5X,A17,1L2)')'Type hermite   = ', Orb%type_hermite
 WRITE(IUNIT,'(5X,A17,1L2)')'Type nucleus   = ', Orb%type_Nucleus

 WRITE(IUNIT,'(5X,A17,1L2)')'Spherical      = ', Orb%spherical
 WRITE(IUNIT,'(5X,A17,I3)') 'maxAngmom      = ', Orb%maxAngmom
 WRITE(IUNIT,'(5X,A17,I3)') 'nAngmom        = ', Orb%nAngmom
 WRITE(IUNIT,'(5X,A17,I3)') 'nPrimitives    = ', Orb%nPrimitives
 WRITE(IUNIT,'(5X,A17,I3)') 'totOrbitals    = ', Orb%totOrbitals
 WRITE(IUNIT,'(5X,A17,I3)') 'maxContracted  = ', Orb%maxContracted
!Print center and atomic information (special case for passes)
 IF (.NOT.Orb%type_Empty) THEN
   IF (nPasses.EQ.1) THEN
     WRITE(IUNIT,'(5X,A17,3F10.6)') 'center (A)     = ', (Orb%center(i),i=1,3)
   ELSE
     printCenter = .TRUE.
     firstAtom = Orb%Atom(1)
     DO iCent=2,nPasses
       IF (firstAtom.NE.Orb%Atom(iCent)) THEN
         printCenter = .FALSE.
         EXIT
       ENDIF
     ENDDO
     IF (printCenter) THEN
       WRITE(IUNIT,'(5X,A17,3F10.6)') 'Pass-center (A) =', (Orb%center(i),i=1,3)
     ELSE
       WRITE(IUNIT,'(5X,A)')        'Pass contains orbitals that do not share a common center'
     ENDIF
   ENDIF
   WRITE(IUNIT,'(5X,A15,10I5,/(20X,10I5))') 'Atomic index = ',(Orb%Atom(iCent),iCent=1,nPasses)
 ENDIF
 WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------'
 WRITE(IUNIT,'(3X,A)')    '-  Block   angmom  #cont.  1.orb.  1.local #orb.   #oc   #pri.oc  1.porb  -'
 DO I=1,Orb%nAngmom
   WRITE(IUNIT,'(1X,9I8)')  I, Orb%angmom(I), Orb%nContracted(I), Orb%startOrbital(I,1), &
     & Orb%startLocOrb(I), Orb%nOrbitals(I), Orb%nOrbComp(I), Orb%nPrimOrbComp(I),Orb%startprimOrbital(I,1)
 ENDDO
 WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------'
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

!> \brief build Ecoefficient tensor (ordering nAng,ijk,nprim) call from set_ftuv
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param TUV The TUV indeces
!> \param Q the overlap distribution
!> \param signQ the sign (1 for P and -1 for Q)
!> \param Ecoeffs the ecoefficients to be built
!> \param ntuv the number of hermite angular components
!> \param Nijk the number of cartesian angular components
!> \param nPrim the number of primitives
!> \param iAngmom the angular moment index
!> \param nPasses the number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE BuildEcoeffTensor1_AP(TUV,Q,signQ,Ecoeffs,ntuv,Nijk,nPrim,iAngmom,nPasses,LUPRI,IPRINT)
!Ordering of Ecoeffs(nAng,ijk,nprim)
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: Q
Integer            :: iAngmom,LUPRI,IPRINT,nPasses
Integer            :: ntuv,Nijk,nPrim
Real(realk)        :: signQ
!Allocatable for ECOEFFS - Dimensions are nPrimitives, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
real(realk),pointer  :: ETIJ(:,:,:,:,:) 
integer             :: Q1,Q2,tQ,uQ,vQ
integer             :: ituvQ,iQ1,jQ1,kQ1,iQ2,jQ2,kQ2,iPrimQ
integer             :: l1,l2
integer             :: startAng, endAng, tuvOff,ijk
Real(realk)         :: Ecoeffs(ntuv,Nijk,nPrim) !nAng,ijk1*ijk2,nprim
Real(realk)         :: signijk

l1 = Q%orbital1%angmom(Q%indexAng1(iangmom))
l2 = Q%orbital2%angmom(Q%indexAng2(iangmom))

startAng = Q%startAngmom
tuvOff   = startAng*(startAng+1)*(startAng+2)/6
endAng   = Q%endAngmom
IF (endAng.LT.l1+l2) THEN
  CALL LSQUIT('Error in BuildEcoeffTensor1',lupri)
ENDIF

!ALLOCATE(ETIJ(nprim,0:l1+l2,0:l1,0:l2,3))
call mem_alloc(ETIJ,nprim,l1+l2,l1,l2,3,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)

CALL GET_ECOEFF(ETIJ,nprim,l1+l2,l1,l2,Q%nPrimitives,nPasses,Q,LUPRI,IPRINT)

CALL LS_DZERO(Ecoeffs,ntuv*Nijk*nPrim)

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
                              IF (tQ+uQ+vQ.GE.startAng) THEN
                                 signijk = signQ**(tQ+uQ+vQ)
                                 ituvQ=TUV%TUVindex(tQ,uQ,vQ)-tuvOff
                                 DO iPrimQ=1,nprim!Q%nPrimitives
                                    Ecoeffs(ituvQ,ijk,iPrimQ) = &
                                         &ETIJ(iPrimQ,tQ,iQ1,iQ2,1)&
                                         &*ETIJ(iPrimQ,uQ,jQ1,jQ2,2)&
                                         &*ETIJ(iPrimQ,vQ,kQ1,kQ2,3)&
                                         &*Q%preExpFac(iPrimQ)*signijk
                                 ENDDO
                              ENDIF
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

call mem_dealloc(ETIJ)

END SUBROUTINE BuildEcoeffTensor1_AP

!> \brief wrapper routine that brach out and build Ecoefficient tensor
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param TUV The TUV indeces
!> \param P the overlap distribution
!> \param signP the sign (1 for P and -1 for Q)
!> \param Ecoeffs the ecoefficients to be built
!> \param lm the number of spherical angular components
!> \param Nijk the number of cartesian angular components
!> \param ntuvP the number of hermite angular components
!> \param nP the number of primitives
!> \param orderAP the order
!> \param iAngmom the angular moment index
!> \param nPasses the number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE BuildEcoeffTensor(TUV,P,signP,Ecoeffs,lm,Nijk,tuvP,nP,orderAP,iAngmom,nPasses,LUPRI,IPRINT)
use memory_handling
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: P
Integer            :: iAngmom,nPasses,LUPRI,IPRINT
Integer            :: lm,Nijk,tuvP,nP
Real(realk)        :: signP
Real(realk),target :: Ecoeffs(lm,tuvP,nP) !ijk1*ijk2,nAng,nprim
Logical            :: orderAP
!
Real(realk),pointer     :: EcoeffN(:,:,:)
Real(realk),pointer :: SpherMat(:,:)
Integer                 :: ijk1,ijk2,lm1,lm2,ang1,ang2,iA1,iA2
Logical                 :: Sph1,Sph2,Spherical

iA1  = P%indexAng1(iangmom)
iA2  = P%indexAng2(iangmom)
ang1 = P%orbital1%angmom(iA1)
ang2 = P%orbital2%angmom(iA2)
Sph1 = P%orbital1%spherical.AND.(ang1.GT.1)
Sph2 = P%orbital2%spherical.AND.(ang2.GT.1)
Spherical = P%sphericalEcoeff.AND.(Sph1.OR.Sph2)
!Test consistency for derivatice case
IF (Spherical.AND.P%derivORder.GT.0) CALL LSQUIT('Error in BuildEcoeffTensor. Spherical and deriv>0!',lupri)
IF (Spherical) THEN
  NULLIFY(EcoeffN)
  CALL mem_alloc(EcoeffN,Nijk,tuvP,nP)
  CALL mem_alloc(SpherMat,Nijk,lm)
ELSE
  EcoeffN => Ecoeffs
ENDIF
IF(orderAP)THEN
   CALL BuildEcoeffTensor_AP(TUV,P,signP,EcoeffN,Nijk,tuvP,nP,iAngmom,nPasses,LUPRI,IPRINT)
ELSE
   CALL BuildEcoeffTensor_PA(TUV,P,signP,EcoeffN,Nijk,tuvP,nP,iAngmom,nPasses,LUPRI,IPRINT)
ENDIF
IF (Spherical) THEN
  ijk1 = (ang1+1)*(ang1+2)/2
  ijk2 = (ang2+1)*(ang2+2)/2
  lm1 = 2*ang1+1
  lm2 = 2*ang2+1
  CALL SphericalTransformation(SpherMat,ijk1,ijk2,lm1,lm2,ang1,ang2,P,&
    &         EcoeffN,Ecoeffs,tuvP*nP,orderAP,LUPRI,IPRINT)
  CALL mem_dealloc(EcoeffN)
  CALL mem_dealloc(SpherMat)
ENDIF

END SUBROUTINE BuildEcoeffTensor

!> \brief build Ecoefficient tensor (ordering ijk,nAng,nprim)
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param TUV The TUV indeces
!> \param P the overlap distribution P or Q
!> \param signP the sign (1 for P and -1 for Q)
!> \param Ecoeffs the ecoefficients to be built
!> \param Nijk the number of cartesian angular components
!> \param ntuv the number of hermite angular components
!> \param nPrim the number of primitives
!> \param iAngmom the angular moment index
!> \param nPasses the number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE BuildEcoeffTensor_AP(TUV,P,signP,Ecoeffs,Nijk,ntuv,nprim,iAngmom,nPasses,LUPRI,IPRINT)
!different ordering of Ecoeffs(ijk,nAng,nprim)
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: P
Integer            :: iAngmom,LUPRI,IPRINT,nPasses
Integer            :: Nijk,ntuv,nprim
Real(realk)        :: signP
!Allocatable for ECOEFFS - Dimensions are nPrimitives, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
real(realk),pointer  :: ETIJ(:,:,:,:,:) 
integer             :: P1,P2,tP,uP,vP
integer             :: ituvP,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
integer             :: l1,l2
Real(realk)         :: Ecoeffs(Nijk,ntuv,nprim) !ijk1*ijk2,nAng,nprim
Real(realk)         :: signijk
Integer             :: der,increment(2),icent,ncent,iPrim1,iPrim2,ijk
Real(realk),pointer :: pref(:)

l1 = P%orbital1%angmom(P%indexAng1(iangmom))
l2 = P%orbital2%angmom(P%indexAng2(iangmom))
der= P%derivOrder

IF (der.GT.0) THEN
  IF (.NOT.P%type_Hermite) CALL LSQUIT('Error in BuildEcoeffTensor_AP. derivOrder>0 and not type hermite',lupri)
ENDIF

!For derivative case the ETIJ can be split into different components if time-critical, i.e.
!ETIJ(l1+l2+1,l1+1,l2) and ETIJ(l1+l2+1,l1,l2+1) for gradients etc.
!ALLOCATE(ETIJ(nprim,0:l1+l2+2*der,0:l1+der,0:l2+der,3))
call mem_alloc(ETIJ,nprim,l1+l2+2*der,l1+der,l2+der,3,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)

CALL GET_ECOEFF(ETIJ,nprim,l1+l2+2*der,l1+der,l2+der,P%nPrimitives,nPasses,P,LUPRI,IPRINT)

CALL LS_DZERO(Ecoeffs,Nijk*ntuv*nprim)

!Derivative settings
ncent = 1 + der
IF (P%single) ncent = 1

pref => P%preExpFac
IF (der.GT.0) THEN
  CALL mem_alloc(pref,P%nPrimitives)
ELSE
  pref => P%preExpFac
ENDIF

ijk=0
!Derivative loop
DO icent=1,ncent
! Nuclear derivative settings (using Hermite gaussians accoring PCCP,2007,9,4771-4779)
! Set up increments
  increment(1) = der+1-icent
  increment(2) = icent - 1
  IF (P%single) THEN
    IF (P%orbital1%TYPE_empty) THEN
      increment(1) = 0
      increment(2) = der
    ELSEIF (P%orbital2%TYPE_empty) THEN
      increment(1) = der
      increment(2) = 0
    ELSE
      CALL LSQUIT('Error in BuildEcoeffTensor_AP. overlap%single and not orbital%empty',lupri)
    ENDIF
  ENDIF

! Set up prefactors that include the (2a)^derA * (2b)^derB
  IF (der.GT.0) THEN
    DO iPrimP=1,P%nPrimitives
      iPrim1 = P%iprim1(iPrimP)
      iPrim2 = P%iprim2(iPrimP)
      pref(iPrimP) = P%preExpFac(iPrimP)*(2.d0*P%orbital1%exponents(iPrim1))**increment(1)*&
     &               (2.d0*P%orbital2%exponents(iPrim2))**increment(2)
    ENDDO
  ENDIF
 
! Regular loop starts here
  DO P2 = 0,l2+increment(2)
     DO iP2=P2,0,-1
        DO jP2=P2-iP2,0,-1
           kP2=P2-iP2-jP2
           DO P1 = 0,l1+increment(1)
              DO iP1=P1,0,-1
                 DO jP1=P1-iP1,0,-1
                    kP1=P1-iP1-jP1
                    IF(iP1+jP1+kP1+iP2+jP2+kP2 .EQ. l1+l2+der) then
                    ijk=ijk+1
                       DO tP=0,iP1+iP2
                          DO uP=0,jP1+jP2
                             DO vP=0,kP1+kP2
                                signijk = signP**(tP+uP+vP)
                                ituvP=TUV%TUVindex(tP,uP,vP)
                                DO iPrimP=1,nprim
                                   Ecoeffs(ijk,ituvP,iPrimP) = &
                                        &ETIJ(iPrimP,tP,iP1,iP2,1)&
                                        &*ETIJ(iPrimP,uP,jP1,jP2,2)&
                                        &*ETIJ(iPrimP,vP,kP1,kP2,3)&
                                        &*pref(iPrimP)*signijk
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
 
ENDDO

IF (der.GT.0) THEN
! CALL PrintTensor(Ecoeffs,'E-coefficients build',&
!      &Nijk,ntuv,nPrim,2,'ijk   ','tuv   ','prim  ',3)
  CALL mem_dealloc(pref)
ENDIF

call mem_dealloc(ETIJ)

END SUBROUTINE BuildEcoeffTensor_AP

!> \brief build Ecoefficient tensor (ordering nPrim,Nijk,ntuv)
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param TUV The TUV indeces
!> \param P the overlap distribution P or Q
!> \param signP the sign (1 for P and -1 for Q)
!> \param Ecoeffs the ecoefficients to be built
!> \param Nijk the number of cartesian angular components
!> \param ntuv the number of hermite angular components
!> \param nPrim the number of primitives
!> \param iAngmom the angular moment index
!> \param nPasses the number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE BuildEcoeffTensor_PA(TUV,P,signP,Ecoeffs,Nijk,ntuv,nprim,iAngmom,nPasses,LUPRI,IPRINT)
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: P
Integer            :: iAngmom,LUPRI,IPRINT,nPasses
Integer            :: Nijk,ntuv,nprim
Real(realk)        :: signP
!Allocatable for ECOEFFS - Dimensions are nPrimitives, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
real(realk),pointer  :: ETIJ(:,:,:,:,:) 
integer             :: P1,P2,tP,uP,vP
integer             :: ituvP,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
integer             :: l1,l2,ijk
Real(realk)         :: Ecoeffs(nPrim,Nijk,ntuv)
Real(realk)         :: signijk

IF (P%derivOrder.GT.0) CALL LSQUIT('BuildEcoeffTensor_PA for derivOrder>0',lupri)

l1 = P%orbital1%angmom(P%indexAng1(iangmom))
l2 = P%orbital2%angmom(P%indexAng2(iangmom))

!ALLOCATE(ETIJ(nprim,0:l1+l2,0:l1,0:l2,3))
call mem_alloc(ETIJ,nprim,l1+l2,l1,l2,3,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)

CALL GET_ECOEFF(ETIJ,nprim,l1+l2,l1,l2,P%nPrimitives,nPasses,P,LUPRI,IPRINT)

CALL LS_DZERO(Ecoeffs,nPrim*Nijk*ntuv)

ijk=0
DO P2 = 0,l2
   DO iP2=P2,0,-1
      DO jP2=P2-iP2,0,-1
         kP2=P2-iP2-jP2
         DO P1 = 0,l1
            DO iP1=P1,0,-1
               DO jP1=P1-iP1,0,-1
                  kP1=P1-iP1-jP1
                  IF(iP1+jP1+kP1+iP2+jP2+kP2 .EQ. l1+l2) then
                  ijk=ijk+1
                     DO tP=0,iP1+iP2
                        DO uP=0,jP1+jP2
                           DO vP=0,kP1+kP2
                              signijk = signP**(tP+uP+vP)
                              ituvP=TUV%TUVindex(tP,uP,vP)
                              DO iPrimP=1,nPrim!P%nPrimitives
                                 Ecoeffs(iPrimP,ijk,ituvP) = &
                                      &ETIJ(iPrimP,tP,iP1,iP2,1)&
                                      &*ETIJ(iPrimP,uP,jP1,jP2,2)&
                                      &*ETIJ(iPrimP,vP,kP1,kP2,3)&
                                      &*P%preExpFac(iPrimP)*signijk
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

call mem_dealloc(ETIJ)

END SUBROUTINE BuildEcoeffTensor_PA

!> \brief wrapper routine for contract with Ecoefficients, special case for single hermite ODbatch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Herint the hermite integrals
!> \param P the overlap distribution P or Q
!> \param signP the sign (1 for P and -1 for Q)
!> \param iAngmom the angular moment index
!> \param nPasses the number of passes
!> \param nAng the number of angular components
!> \param nOrb the number of orbital components
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SingleHermiteEcoeff(HerInt,P,signP,iAngmom,nPasses,nAng,nOrb,LUPRI,IPRINT)
implicit none
Real(realk)   :: HerInt(:)
TYPE(Overlap) :: P
Integer       :: iAngmom,nPasses,nAng,nOrb,LUPRI,IPRINT
Real(realk)   :: signP
!
Integer     :: nPrim1,nPrim2
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
CALL SingleHermiteEcoeff1(HerInt,P,sign,j1,j2,nPrim1,nPrim2,nPasses,nAng,nOrb,LUPRI,IPRINT)
END SUBROUTINE SingleHermiteEcoeff

!> \brief contract with Ecoefficients, special case for single hermite ODbatch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Herint the hermite integrals
!> \param P the overlap distribution P or Q
!> \param signP the sign (1 for P and -1 for Q)
!> \param j1 angular moment for orbital 1
!> \param j2 angular moment for orbital 2
!> \param nPrim1 number of primitives for orbital 1
!> \param nPrim2 number of primitives for orbital 2
!> \param nPasses the number of passes
!> \param nAng the number of angular components
!> \param nOrb the number of orbital components
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SingleHermiteEcoeff1(HerInt,P,signP,j1,j2,nPrim1,nPrim2,nPasses,nAng,nOrb,LUPRI,IPRINT)
implicit none
Real(realk)   :: HerInt(nAng*nOrb,nPrim1,nPrim2,nPasses)
TYPE(Overlap) :: P
Integer       :: nAng,nOrb,LUPRI,IPRINT
Integer       :: nPrim1,nPrim2,j1,j2,nPasses
Real(realk)   :: signP
!
Real(realk)   :: pref2,pref12,D1=1.0d0,D2=2.0d0
Integer       :: iPrim1,iPrim2,iAngOrb,iPasses
DO iPasses=1,nPasses
  DO iPrim2=1,nPrim2
    pref2=signP/(D2*P%orbital2%exponents(iPrim2))**j2
    DO iPrim1=1,nPrim1
      pref12=D1/(D2*P%orbital1%exponents(iPrim1))**j1*pref2
      DO iAngOrb=1,nAng*nOrb
        HerInt(iAngOrb,iPrim1,iPrim2,iPasses) = HerInt(iAngOrb,iPrim1,iPrim2,iPasses) * pref12
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE SingleHermiteEcoeff1

!> \brief print Ecoefficients
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param EcoeffCont the Ecoefficient to be printed
!> \param nAng the number of angular components
!> \param nOrb the number of orbital components
!> \param nPrim the number of primitives
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
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

!> \brief wrapper routine for spherical transformation of integrals
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spherical the Spherical transformation matrix 
!> \param ijk1 the number of cartesian angular components of orbital 1
!> \param ijk2 the number of cartesian angular components of orbital 2
!> \param lm1 the number of spherical angular components of orbital 1
!> \param lm2 the number of spherical angular components of orbital 2
!> \param ang1 the angular momentum for orbital 1
!> \param ang2 the angular momentum for orbital 2
!> \param P the overlap distribution either P or Q
!> \param integralNonSp the input array
!> \param integralSpher the output array
!> \param dim1 the dimension of the array to be transformed
!> \param orderAngPrim flag which determines the order of the intermidiate 
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SphericalTransformation(Spherical,ijk1,ijk2,lm1,lm2,ang1,ang2,P,&
    &         integralNonSp,integralSpher,dim1,orderAngPrim,lupri,iprint)
implicit none
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,dim1
Logical            :: orderAngPrim
!
Real(realk)         :: Spherical(ijk1*ijk2,lm1*lm2)
Real(realk)         :: integralNonSp(*)
Real(realk)         :: integralSpher(*)
Integer             :: lm1,lm2,ijk1,ijk2,ang1,ang2

! Sets up spherical transformation matrix
CALL Buildsphericaltransformation(Spherical,P%orbital1%SPH_MAT(ang1+1)%p%elms,&
     &              P%orbital2%SPH_MAT(ang2+1)%p%elms,lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)

! Spherical transformation either on first or second dimension
IF (orderAngPrim) THEN
  ! First dimension
  CALL SphericalTransformation_AP(Spherical,ijk1*ijk2,lm1*lm2,integralNonSp,integralSpher,&
     &                            dim1,lupri,iprint)
ELSE
  ! Second dimension
  CALL SphericalTransformation_PA(Spherical,ijk1*ijk2,lm1*lm2,integralNonSp,integralSpher,&
     &                            dim1,lupri,iprint)
ENDIF
!  
END SUBROUTINE SphericalTransformation

!> \brief spherical transformation of integrals order AP
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spherical the Spherical transformation matrix 
!> \param ijk the number of cartesian angular components
!> \param lm the number of spherical angular components
!> \param integralNonSp the input array
!> \param integralSpher the output array
!> \param dim1 the dimension of the array to be transformed
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SphericalTransformation_AP(Spherical,ijk,lm,integralNonSp,integralSpher,&
     &                                dim1,lupri,iprint)
implicit none
Integer            :: LUPRI,IPRINT
Integer            :: lm,ijk,dim1
Real(realk)        :: Spherical(ijk,lm)
Real(realk)        :: integralNonSp(ijk,dim1)
Real(realk)        :: integralSpher(lm,dim1)
!
Real(realk)         :: D1=1.0d0,D0=0.0d0

! Spherical transformation (with both centers beloning to one electron simultaneously)
  CALL DGEMM('T','N',lm,dim1,ijk,D1,Spherical,ijk,&
     &     integralNonSp,ijk,D0,integralSpher,lm)

END SUBROUTINE SphericalTransformation_AP

!> \brief spherical transformation of integrals order PA
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spherical the Spherical transformation matrix 
!> \param ijk the number of cartesian angular components
!> \param lm the number of spherical angular components
!> \param integralNonSp the input array
!> \param integralSpher the output array
!> \param dim1 the dimension of the array to be transformed
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SphericalTransformation_PA(Spherical,ijk,lm,integralNonSp,integralSpher,&
     &                                dim1,lupri,iprint)
implicit none
Integer            :: LUPRI,IPRINT
Integer            :: lm,ijk, dim1
Real(realk)        :: Spherical(ijk,lm)
Real(realk)        :: integralNonSp(dim1,ijk)
Real(realk)        :: integralSpher(dim1,lm)
!
Real(realk)        :: D1=1.0d0,D0=0.0d0

!! Spherical transformation (with both centers beloning to one electron simultaneously)
  CALL DGEMM('N','N',dim1,lm,ijk,D1,integralNonSp,dim1,&
     &     Spherical,ijk,D0,integralSpher,dim1)

END SUBROUTINE SphericalTransformation_PA

!> \brief construct the spherical transformation matrix
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spherical the Spherical transformation matrix 
!> \param Spher1 the Spherical transformation matrix for orbital 1 
!> \param Spher2 the Spherical transformation matrix for orbital 2
!> \param lm1 the number of spherical angular components for orbital 1 
!> \param lm2 the number of spherical angular components for orbital 2 
!> \param ijk1 the number of cartesian angular components for orbital 1 
!> \param ijk2 the number of cartesian angular components for orbital 2 
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Buildsphericaltransformation(Spherical,Spher1,Spher2,lm1,lm2,&
     &                                      ijk1,ijk2,LUPRI,IPRINT)
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
END SUBROUTINE Buildsphericaltransformation

!> \brief print the spherical transformation matrix
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spherical the Spherical transformation matrix 
!> \param lm1 the number of spherical angular components for orbital 1 
!> \param lm2 the number of spherical angular components for orbital 2 
!> \param ijk1 the number of cartesian angular components for orbital 1 
!> \param ijk2 the number of cartesian angular components for orbital 2 
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
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

!> \brief construct the contraction matrix, order AP
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param CC the contraction matrix
!> \param P the overlap distribution P or Q
!> \param nP the number of primitives
!> \param nC1 the number of contracted functions on orbital 1
!> \param nC2 the number of contracted functions on orbital 2
!> \param iA1 the angular moment index on orbital 1
!> \param iA2 the angular moment index on orbital 2
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ConstructContraction_AP(CC,P,nP,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
Implicit none
Type(Overlap) :: P
Integer       :: iA1,iA2,nP,nC1,nC2,LUPRI,IPRINT
Real(realk)   :: CC(nP,nC1,nC2)
!
Integer       :: iP,i1,i2,iC1,iC2,nP1,nP2

nP1 = P%orbital1%nPrimitives
nP2 = P%orbital2%nPrimitives
IF(P%sameAO)THEN
   DO iC2=1,nC2
      DO iC1=1,nC1
         DO iP=1,nP
            i1 = P%iprim1(iP)
            i2 = P%iprim2(iP)
            IF(i1 .NE. i2)THEN
            CC(iP,iC1,iC2)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2)&
                 &        +P%orbital1%CC(iA1)%p%elms(i2+(iC1-1)*nP2)*P%orbital2%CC(iA2)%p%elms(i1+(iC2-1)*nP1)


            ELSE
            CC(iP,iC1,iC2)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2)
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
            CC(iP,iC1,iC2)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2) 
         ENDDO
      ENDDO
   ENDDO
ENDIF

IF (IPRINT.GT.30) THEN
   CALL PrintTensor(CC,'Overlap CC          ',nP,nC1,nC2,Lupri,&
        & 'prim  ','C1    ','C2    ',1)
ENDIF

END SUBROUTINE ConstructContraction_AP

!> \brief construct the contraction matrix order PA
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param CC the contraction matrix
!> \param P the overlap distribution P or Q
!> \param nP the number of primitives
!> \param nC1 the number of contracted functions on orbital 1
!> \param nC2 the number of contracted functions on orbital 2
!> \param iA1 the angular moment index on orbital 1
!> \param iA2 the angular moment index on orbital 2
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ConstructContraction_PA(CC,P,nP,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
Implicit none
Type(Overlap) :: P
Integer       :: iA1,iA2,nP,nC1,nC2,LUPRI,IPRINT
Real(realk)   :: CC(nC1,nC2,nP)
!
Integer       :: iP,i1,i2,iC1,iC2,nP1,nP2

nP1 = P%orbital1%nPrimitives
nP2 = P%orbital2%nPrimitives
IF(P%sameAO)THEN
   DO iP=1,nP
      i1 = P%iprim1(iP)
      i2 = P%iprim2(iP)
      IF(i1 .NE. i2)THEN
         DO iC2=1,nC2
            DO iC1=1,nC1
               CC(iC1,iC2,iP)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2)&
                    &        +P%orbital1%CC(iA1)%p%elms(i2+(iC1-1)*nP2)*P%orbital2%CC(iA2)%p%elms(i1+(iC2-1)*nP1)
            ENDDO
         ENDDO
      ELSE
         DO iC2=1,nC2
            DO iC1=1,nC1
               CC(iC1,iC2,iP)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2)
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ELSE
   DO iP=1,nP
      i1 = P%iprim1(iP)
      i2 = P%iprim2(iP)
      DO iC2=1,nC2
         DO iC1=1,nC1
            CC(iC1,iC2,iP)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2) 
         ENDDO
      ENDDO
   ENDDO
ENDIF
!IF (IPRINT.GT.30) THEN
!   CALL PrintTensor(CC,'Overlap CC          ',nP,nC1,nC2,Lupri,&
!        & 'prim  ','C1    ','C2    ',1)
!ENDIF

END SUBROUTINE ConstructContraction_PA

!> \brief Print the contraction matrix
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param CC the contraction matrix
!> \param nP the number of primitives
!> \param nC1 the number of contracted functions on orbital 1
!> \param nC2 the number of contracted functions on orbital 2
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
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

!> \brief build the single cartesian Ecoeff (ETIJ)
!> \author A. Teale
!> \date 2009
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimpass the number of primitives * number of npasses
!> \param maxang12 the combined maximum angular momentum
!> \param maxang1 the maximum angular momentum for orbital 1
!> \param maxang2 the maximum angular momentum for orbital 2
!> \param nprim the number of primitives
!> \param npass the number of passes
!> \param P the overlap distribution P or Q
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE GET_ECOEFF(ETIJ,nprimpass,maxang12,maxang1,maxang2,nprim,npass,P,LUPRI,IPRINT)
implicit none
REAL(REALK),PARAMETER    :: D1 = 1.0D0, D2 = 2.0D0, DHALF = 0.5D0
integer                  :: IPRINT, LUPRI, I, J, K, CARTDIR
integer                  :: nprim,maxang12,maxang1,maxang2
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
IF (.NOT.P%type_Empty)THEN
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
      IF (P%type_hermite .OR. P%type_hermite_single) THEN
           CALL HERM_ECOEFFS(ETIJ(1,0,0,0,CARTDIR),nPrimPass,maxAng1,&
          &             maxAng2,PA(:,CARTDIR),PB(:,CARTDIR),HPINV,TWOA,TWOB,LUPRI)
      ELSE IF (P%type_cartesian .OR. P%type_cartesian_single) THEN
           CALL CART_ECOEFFS(ETIJ(1,0,0,0,CARTDIR),nPrimPass,maxAng1,&
          &             maxAng2,PA(:,CARTDIR),PB(:,CARTDIR),HPINV,LUPRI)
      ELSE IF (P%type_Empty) THEN
           DO i=1,nPrimPass
              ETIJ(i,0,0,0,CARTDIR) = D1
           ENDDO
      ELSE IF (P%type_Nucleus) THEN
           DO i=1,nPrimPass
              ETIJ(i,0,0,0,CARTDIR) = D1
           ENDDO
      ELSE
           WRITE(LUPRI,*)'ERROR : UNRECOGNISED TYPE OF OVERLAP!!!'
      ENDIF
      IF (IPRINT.GT.20) CALL PRN_ECOEFFS(ETIJ,nPrimPass,maxAng1,&
           &            maxAng2,CARTDIR,LUPRI)
ENDDO
! DEALLOCATE MEMORY
END SUBROUTINE GET_ECOEFF

!> \brief SUBROUTINE TO CALCULATE E's USING THE USUAL CARTESIAN RECURANNCE RELATIONS - SEE P354 OF BOOK
!> \author A. Teale
!> \date 2009
!>
!> SEE ERIECF ROUTINE for F77 EQUIVALENT
!>
!> VARIABLES
!> =========      
!> nPrimPairs NUMBER OF PRIMITIVE PAIRS IN OVERLAP DISTRIBUTION
!> MAXI,MAXJ  MAXIMUM I,J VALUES IN ETIJ's (MAX ANGMOM OF G_i and G_j RESPECTIVELY)
!> I,J,K,IJ   COUNTERS     
!> PA         DISTANCE BETWEEN CENTRE OF CHARGE P AND A (THE CENTRE OF PRIMITIVE GAUSSIAN G_i)
!> PB         DISTANCE BETWEEN CENTRE OF CHARGE P AND B (THE CENTRE OF PRIMITIVE GAUSSIAN G_j)
!> HPINV      1/2p (p = a + b I.E. SUM OF EXPONENTS OF PRIM G_i and G_j)           
!> PINV       1/p
!> APINV      a/p
!> BPINV      b/p
!> ETIJ       E COEFFICIENTS 
!> X          INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> T          COUNTER INDEX USED TO COUNT UPTO SUM OF ANGULAR MOMENTA OF G_i AND G_j 
!>
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param LUPRI the logical unit number for the output file
SUBROUTINE CART_ECOEFFS(ETIJ,nPrimPairs,MAXI,MAXJ,PA,PB,HPINV,LUPRI)
implicit none
REAL(REALK),PARAMETER  ::  D1 = 1.0D0, D2 = 2.0D0
INTEGER                ::  I,J,K
INTEGER                ::  T,nPrimPairs,MAXI,MAXJ,IJ,LUPRI
REAL(REALK)            ::  T1
REAL(REALK)            ::  PA(nPrimPairs), PB(nPrimPairs) 
REAL(REALK)            ::  ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
REAL(REALK)            ::  HPINV(nPrimPairs)


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

!> \brief helper routine for building the ETIJ, J=0  and I .LE. 2
!> \author A. Teale
!> \date 2009
!>
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param HPINV the argument 1/2p
!> \param LUPRI the logical unit number for the output file
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

!> \brief helper routine for building the ETIJ
!> \author A. Teale
!> \date 2009
!>
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param I counter
!> \param IJ counter
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

!> \brief print the ecoeff tensor
!> \author A. Teale
!> \date 2009
!>
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param X an INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> \param LUPRI the logical unit number for the output file
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

!> \brief print the ecoeff info
!> \author A. Teale
!> \date 2009
!>
!> \param P the overlap distribution P or Q
!> \param LUPRI the logical unit number for the output file
SUBROUTINE PRN_ECINFO(P,LUPRI)
implicit none
integer                  :: LUPRI, I
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
WRITE(LUPRI,*)'TYPE hermite  ',P%orbital1%TYPE_Hermite
WRITE(LUPRI,*)'TYPE empty    ',P%orbital1%TYPE_Empty
WRITE(LUPRI,*)'TYPE cartesian',P%orbital1%TYPE_Cartesian
WRITE(LUPRI,*)'TYPE nucleus  ',P%orbital1%TYPE_Nucleus

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
WRITE(LUPRI,*)'TYPE hermite  ',P%orbital2%TYPE_Hermite
WRITE(LUPRI,*)'TYPE empty    ',P%orbital2%TYPE_Empty
WRITE(LUPRI,*)'TYPE cartesian',P%orbital2%TYPE_Cartesian
WRITE(LUPRI,*)'TYPE nucleus  ',P%orbital2%TYPE_Nucleus

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

!> \brief SUBROUTINE TO CALCULATE E's USING THE HERMITE RECURANNCE RELATIONS
!> \author A. Teale
!> \date 2009
!>
!> VARIABLES
!> =========      
!> nPrimPairs   NUMBER OF PRIMITIVE PAIRS IN OVERLAP DISTRIBUTION
!> MAXI,MAXJ    MAXIMUM I,J VALUES IN ETIJ's (MAX ANGMOM OF G_i and G_j RESPECTIVELY)
!> I,J,K,IJ     COUNTERS     
!> PA           DISTANCE BETWEEN CENTRE OF CHARGE P AND A (THE CENTRE OF PRIMITIVE GAUSSIAN G_i)           
!> PB           DISTANCE BETWEEN CENTRE OF CHARGE P AND B (THE CENTRE OF PRIMITIVE GAUSSIAN G_j)
!> HPINV        1/2p (p = a + b I.E. SUM OF EXPONENTS OF PRIM G_i and G_j)           
!> PINV         1/p
!> APINV, BPINV a/p b/p
!> TWOA, TWOB   2a, 2b
!> ETIJ         E COEFFICIENTS 
!> X            INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> T            COUNTER INDEX USED TO COUNT UPTO SUM OF ANGULAR MOMENTA OF G_i AND G_j 
!>
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param TWOA 2a (a is exponent on orbital 1)
!> \param TWOB 2b (b is exponent on orbital 2)
!> \param LUPRI the logical unit number for the output file
SUBROUTINE HERM_ECOEFFS(ETIJ,nPrimPairs,MAXI,MAXJ,PA,PB,&
                      & HPINV,TWOA,TWOB,LUPRI)
implicit none
REAL(REALK),PARAMETER  ::  D1 = 1.0D0, D2 = 2.0D0
INTEGER                ::  I,J,K
INTEGER                ::  T,nPrimPairs,MAXI,MAXJ,IJ,LUPRI
REAL(REALK)            ::  T1, TIM, TJM
REAL(REALK)            ::  PA(nPrimPairs), PB(nPrimPairs)
REAL(REALK)            ::  ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
REAL(REALK)            ::  HPINV(nPrimPairs)
REAL(REALK)            ::  TWOA(nPrimPairs),TWOB(nPrimPairs)
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

!> \brief helper routine for building the hermite ETIJ
!> \author A. Teale
!> \date 2009
!>
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param HPINV the argument 1/2p
!> \param TWOA 2a (a is exponent on orbital 1)
!> \param I counter
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

!> \brief second helper routine for building the hermite ETIJ
!> \author A. Teale
!> \date 2009
!>
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param TWOB 2b (b is exponent on orbital 2)
!> \param I counter
!> \param IJ counter
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

!> \brief free the overlap structure
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P the overlap structure
SUBROUTINE FREE_OVERLAP(P)
use memory_handling
Implicit none
TYPE(Overlap) :: P

CALL FREE_ORBITAL(P%orbital1)
CALL FREE_ORBITAL(P%orbital2)
 
!$OMP CRITICAL (memory) 
mem_allocated_overlap = mem_allocated_overlap - mem_realsize*size(P%center)&
     &- mem_realsize*size(P%exponents) - mem_realsize*size(P%reducedExponents) &
     &- mem_realsize*size(P%preExpFac) - mem_intsize*size(P%iprim1) - mem_intsize*size(P%iprim2)&
     &- mem_intsize*size(P%lenETUV)
mem_allocated_overlap = mem_allocated_overlap - mem_realsize*size(P%distance12)
mem_allocated_overlap = mem_allocated_overlap - mem_intsize*size(P%angmom)&
     & - mem_intsize*size(P%indexAng1) - mem_intsize*size(P%indexAng2) &
     & - mem_intsize*size(P%nOrbitals) - mem_intsize*size(P%nOrbComp)&
     & - mem_intsize*size(P%nContracted)
if (mem_allocated_overlap < 0) then
   call LSQUIT('Error using mem_allocated_overlap - probably integer overflow!',-1)
endif
!$OMP END CRITICAL (memory) 
Call Mem_dealloc(P%center)
Call Mem_dealloc(P%exponents)  
Call Mem_dealloc(P%reducedExponents)
Call Mem_dealloc(P%preExpFac)
Call Mem_dealloc(P%lenETUV)
Call Mem_dealloc(P%iprim1)
Call Mem_dealloc(P%iprim2)
Call Mem_dealloc(P%distance12)
Call Mem_dealloc(P%angmom)
Call Mem_dealloc(P%indexAng1)
Call Mem_dealloc(P%indexAng2)
Call Mem_dealloc(P%nOrbitals)
Call Mem_dealloc(P%nOrbComp)
Call Mem_dealloc(P%nContracted)

IF (P%ETUVisALLOC) THEN
!$OMP CRITICAL (memory) 
  mem_allocated_overlap = mem_allocated_overlap - mem_realsize*size(P%ETUV)&
       & - mem_intsize*size(P%ETUVindex)
if (mem_allocated_overlap < 0) then
   call LSQUIT('Error using mem_allocated_overlap - probably integer overflow!',-1)
endif
!$OMP END CRITICAL (memory) 
  CALL MEM_DEALLOC(P%ETUV)
  CALL MEM_DEALLOC(P%ETUVindex)
  P%ETUVisALLOC = .FALSE.
ENDIF
IF (P%FTUVisALLOC) THEN
!$OMP CRITICAL (memory) 
   mem_allocated_overlap = mem_allocated_overlap - mem_realsize*size(P%FTUV)
if (mem_allocated_overlap < 0) then
   call LSQUIT('Error using mem_allocated_overlap - probably integer overflow!',-1)
endif
!$OMP END CRITICAL (memory) 
  CALL MEM_DEALLOC(P%FTUV)
  P%FTUVisALLOC = .FALSE.
ENDIF

END SUBROUTINE FREE_OVERLAP

!> \brief free the orbital structure
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Orb the orbital structure
SUBROUTINE FREE_ORBITAL(Orb)
use memory_handling
IMPLICIT NONE
TYPE(Orbital) :: Orb
 
!$OMP CRITICAL (memory) 
mem_allocated_overlap = mem_allocated_overlap - mem_intsize*size(Orb%angmom)&
     & - mem_intsize*size(Orb%nContracted) -  mem_intsize*size(Orb%nOrbComp)&
     & - mem_intsize*size(Orb%startOrbital) - mem_intsize*size(Orb%atom)&
     & - mem_intsize*size(Orb%startLocOrb) - mem_intsize*size(Orb%batch)
if (mem_allocated_overlap < 0) then
   call LSQUIT('Error using mem_allocated_overlap - probably integer overflow!',-1)
endif
!$OMP END CRITICAL (memory) 

IF (.NOT.ASSOCIATED(Orb%angmom)) CALL LSQUIT('Error in FREE_ORBITAL - angmom',-1)
IF (.NOT.ASSOCIATED(Orb%nContracted)) CALL LSQUIT('Error in FREE_ORBITAL - nContracted',-1)
IF (.NOT.ASSOCIATED(Orb%startOrbital)) CALL LSQUIT('Error in FREE_ORBITAL - startOrbital',-1)
IF (.NOT.ASSOCIATED(Orb%startLocOrb)) CALL LSQUIT('Error in FREE_ORBITAL - startLocOrb',-1)
IF (.NOT.ASSOCIATED(Orb%nOrbComp)) CALL LSQUIT('Error in FREE_ORBITAL - nOrbComp',-1)
IF (.NOT.ASSOCIATED(Orb%atom)) CALL LSQUIT('Error in FREE_ORBITAL - atom',-1)
IF (.NOT.ASSOCIATED(Orb%batch)) CALL LSQUIT('Error in FREE_ORBITAL - batch',-1)
CALL MEM_DEALLOC(Orb%angmom)
CALL MEM_DEALLOC(Orb%nContracted)
CALL MEM_DEALLOC(Orb%startOrbital)
CALL MEM_DEALLOC(Orb%startLocOrb)
CALL MEM_DEALLOC(Orb%nOrbComp)
CALL MEM_DEALLOC(Orb%atom)
CALL MEM_DEALLOC(Orb%batch)
 
IF(.NOT. Orb%FTUVorb)THEN
   !$OMP CRITICAL (memory) 
   mem_allocated_overlap = mem_allocated_overlap - mem_realsize*size(Orb%exponents)
   mem_allocated_overlap = mem_allocated_overlap - mem_intsize*size(Orb%startprimOrbital)&
        & - mem_intsize*size(Orb%nPrimOrbComp) - mem_intsize*size(Orb%nOrbitals)
   if (mem_allocated_overlap < 0) then
      call LSQUIT('Error using mem_allocated_overlap - probably integer overflow!',-1)
   endif
   !$OMP END CRITICAL (memory) 
   IF (.NOT.ASSOCIATED(Orb%exponents)) CALL LSQUIT('Error in FREE_ORBITAL - exponents',-1)
   IF (.NOT.ASSOCIATED(Orb%startprimOrbital)) CALL LSQUIT('Error in FREE_ORBITAL - startprimOrbital',-1)
   IF (.NOT.ASSOCIATED(Orb%nPrimOrbComp)) CALL LSQUIT('Error in FREE_ORBITAL - nPrimOrbComp',-1)
   IF (.NOT.ASSOCIATED(Orb%nOrbitals)) CALL LSQUIT('Error in FREE_ORBITAL - nOrbitals',-1)
   CALL MEM_DEALLOC(Orb%exponents)
   CALL MEM_DEALLOC(Orb%startprimOrbital)
   CALL MEM_DEALLOC(Orb%nPrimOrbComp)
   CALL MEM_DEALLOC(Orb%nOrbitals)
   IF (.NOT.ASSOCIATED(Orb%CC)) CALL LSQUIT('Error in FREE_ORBITAL - CC',-1)
   DEALLOCATE(Orb%CC)
   NULLIFY(Orb%CC)   
   IF(ASSOCIATED(Orb%SPH_MAT)) THEN
      DEALLOCATE(Orb%SPH_MAT)
      NULLIFY(Orb%SPH_MAT)
   ELSE
      NULLIFY(Orb%SPH_MAT)
   ENDIF
ENDIF

END SUBROUTINE FREE_ORBITAL

!> \brief set the FTUV batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param F The FTUV overlap
!> \param NFTUVbatches the number of FTUV batches
!> \param OD The OD-batch
!> \param Q The right or left hand side overlap distribution
!> \param Input The integral input
!> \param sharedTUV The TUV indeces
!> \param Integral The integral specifications
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param ndmat The number of density matrix
!> \param maxPrimPass maximum number of primitive*npasses
!> \param lupri The default print unit
!> \param iprint The print level (the higher the more information)
SUBROUTINE SET_FTUVbatches(F,NFTUVbatches,OD,Q,Input,SharedTUV,Integral,Alloc,ndmat,maxPrimPass,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)   :: INPUT
TYPE(Overlap),pointer :: F(:)
TYPE(Overlap),pointer :: Q(:)
Integer               :: NFTUVbatches,maxPrimPass,LUPRI,IPRINT
TYPE(ODITEM)          :: OD 
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)    :: Alloc
TYPE(TUVitem)         :: SharedTUV
!
Integer :: I,J,K,nUnique,ftuvindex,np,ndmat
Integer  ,allocatable :: ODtoFTUVindex(:),nPrimOD(:),FTUVtoPassIndex(:)
Integer(KIND=8),allocatable :: Identifier(:), UniqeIdentfiers(:)
Integer  ,allocatable :: FTUVprim(:),FTUVPassType(:)
Integer  ,allocatable :: FTUVminAng(:),FTUVmaxAng(:),FTUVntuv(:)
Integer  ,allocatable :: offPrim(:)
Logical :: unique,uniquePassType
Integer :: nTUV,nPassTypes,iPassType
Integer :: CSscreenLOG
!Simen  maxPrim is an empirical parameter fixed to reduce CPU-timings. 
!Simen  Two opposing effects occur: 
!Simen     1. timings in ContractFTUV reduce with incresing number of primitives
!Simen        per FTUV-batch.
!Simen     2. timings in GET_WTUV increase with increasing number of primitives.
!Simen  Parameter was fixed using ifort/SUN X2200 Opteron/2.2 Ghz (titan.uio.no)
Integer :: maxPrim = 64

ALLOCATE(ODtoFTUVindex(OD%nbatches))
ALLOCATE(nPrimOD(OD%nbatches))
ALLOCATE(FTUVtoPassIndex(OD%nbatches))
ALLOCATE(Identifier(OD%nbatches))
ALLOCATE(UniqeIdentfiers(OD%nbatches))
ALLOCATE(FTUVprim(OD%nbatches))
ALLOCATE(FTUVPassType(OD%nbatches))
ALLOCATE(FTUVminAng(OD%nbatches))
ALLOCATE(FTUVmaxAng(OD%nbatches))
ALLOCATE(FTUVntuv(OD%nbatches))
ALLOCATE(offPrim(OD%nbatches))

IF(INPUT%DO_FMM.OR.INPUT%OE_SCREEN)THEN
   maxprim = 1
ELSE
   maxprim = Input%FTUVmaxprim
ENDIF

nUnique = 0
nPassTypes=0
DO I = 1,OD%nbatches
 CALL SET_Overlap(Q(I),nP,Input,SharedTUV,Integral,Alloc,OD%BATCH(I),2,LUPRI,IPRINT,'RHS')
 IF (nP.GT.0) THEN
  nPrimOD(I) = np
!Simen  Add center information to this identifier - maybe some integer value of 
!Simen  a weighted OD-center divided by BOXSIZE. Note the larger the boxsize, the 
!Simen  less effective CS-screening becomes, whereas the FTUV-contraction step becomes
!Simen  faster. I other words there is a tradeoff between two opposing effects here.
!Simen  Note also that Maxgab should be inherited from OD-batches to FTUV-batches
  IF (Q(I)%maxAngmom.GT.98) CALL LSQUIT('Non-unique angmom-identifier in SET_FTUVbatches',lupri)
  IF (INPUT%CS_SCREEN) THEN
    CSscreenLOG = max(0,10-int(log10(Q(i)%maxGab)))
  ELSE
    CSscreenLOG = 0
  ENDIF
  IF (CSscreenLOG.GT.998) CALL LSQUIT('Non-unique CSscreenLOG-identifier in SET_FTUVbatches',lupri)
  Identifier(I) = (CSscreenLOG+1)*10000+(Q(I)%maxAngmom+1)*100 + (Q(I)%minAngmom+1)
  unique = .true.
  uniquePassType = .true.
  DO J=nUnique,1,-1
    IF (Identifier(I).EQ.UniqeIdentfiers(J)) THEN
      iPassType = FTUVpassType(J)
      uniquePassType = .false.
      IF ((FTUVprim(J)+np).GT.maxPrim) EXIT
      ftuvIndex = J
      unique = .false.
      exit
     ENDIF
  ENDDO
  IF (uniquePassType) THEN
    nPassTypes = nPassTypes + 1
    iPassType  = nPassTypes
  ENDIF
  IF (unique) THEN
    nUnique = nUnique + 1
    ftuvIndex = nUnique
    UniqeIdentfiers(nUnique) = Identifier(I)
    FTUVpassType(ftuvIndex) = iPassType
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
  Alloc%maxPrimRHS    = max(np,Alloc%maxPrimRHS,maxPrimPass)
  Alloc%maxPrimTUVRHS = max(np*nTUV,maxPrimPass*nTUV,Alloc%maxPrimTUVRHS)
ENDDO

!Initialize FTUV-batches
ALLOCATE(F(NFTUVbatches))
I = 0
DO iPassType=1,nPassTypes
  DO J=1,NFTUVbatches
    IF (FTUVpassType(J).EQ.iPassType) THEN
      I = I + 1
      FTUVtoPassIndex(J) = I
      CALL InitFTUVbatches(F(I),iPassType,FTUVprim(J),FTUVminAng(J),FTUVmaxAng(J),FTUVntuv(J),ndmat)
      offPrim(I) = 0
    ENDIF
  ENDDO
ENDDO

!Add OD-batches to corresponding FTUV-batches
DO I=1,OD%nbatches
 J = ODtoFTUVindex(I)
 IF (J.GT.0) THEN
  K = FTUVtoPassIndex(J)
  np = nPrimOD(I)
  CALL AddODtoFTUV(F(K),np,offPrim(K),Q(I),ndmat,LUPRI,IPRINT)
  CALL FREE_OVERLAP(Q(I))
  offPrim(K) = offPrim(K) + np
 ENDIF
ENDDO

DEALLOCATE(ODtoFTUVindex,nPrimOD,FTUVtoPassIndex,Identifier,UniqeIdentfiers,FTUVprim)
DEALLOCATE(FTUVPassType,FTUVminAng,FTUVmaxAng,FTUVntuv,offPrim)

END SUBROUTINE SET_FTUVbatches

!> \brief add Overlap distribution to FTUV
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param F The FTUV overlap
!> \param np number of primitive functions 
!> \param offp an offset 
!> \param Q The right or left hand side overlap distribution
!> \param ndmat The number of density matrix
!> \param lupri The default print unit
!> \param iprint The print level (the higher the more information)
SUBROUTINE AddODtoFTUV(F,np,offp,Q,ndmat,LUPRI,IPRINT)
implicit none
TYPE(Overlap) :: F,Q
Integer       :: np,offp,ndmat,LUPRI,IPRINT
!
Integer :: ip,idir,idmat,tuv
IF (F%nTUV.NE.Q%nTUV) THEN
  CALL LSQUIT('Programming error: nTUV maismatch in AddODtoFTUV',lupri)
ENDIF
F%maxGab = max(F%maxGab,Q%maxGab)
F%ODextent = Q%ODextent
F%ODcenter = Q%ODcenter
F%nPrimitives = F%nPrimitives + np

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

!> \brief init the FTUVbatches
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param F The FTUV overlap
!> \param iPasstype the type of this overlap
!> \param nPrim number of primitives
!> \param minAng the minimum angular momentum
!> \param maxAng the maximum angular momentum
!> \param ntuv number of angular components
!> \param ndmat The number of density matrix
SUBROUTINE InitFTUVbatches(F,iPassType,nPrim,minAng,maxAng,nTUV,ndmat)
implicit none
TYPE(Overlap) :: F
Integer       :: nPrim,minAng,maxAng,nTUV,ndmat,iPassType
!
CALL ALLOCATE_ORBITAL(F%orbital1,1,1,1,1,.TRUE.)
CALL ALLOCATE_ORBITAL(F%orbital2,1,1,1,1,.TRUE.)
CALL ALLOCATE_OVERLAP(F,nprim,1,1,1,.TRUE.,nTUV,ndmat,.FALSE.,nTUV,.FALSE.)

F%TYPE_Empty = .FALSE.
F%TYPE_Hermite = .FALSE.
F%TYPE_Hermite_single = .FALSE.
F%TYPE_Cartesian = .FALSE.
F%TYPE_Cartesian_single = .FALSE.
F%TYPE_Nucleus = .FALSE.
F%TYPE_FTUV = .TRUE.

F%single = .FALSE.

F%passType                 = iPassType
F%nAngmom                  = 1
F%orbital1%nAngmom         = 1
F%orbital2%nAngmom         = 1
F%nPrimitives              = 0
F%nPasses                  = 1
F%minAngmom                = minAng
F%maxAngmom                = maxAng
F%startAngmom              = minAng
F%derivOrder               = 0
F%endAngmom                = maxAng
F%nTUV                     = nTUV
F%totOrbitals              = ndmat
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
F%maxGab                   = 0.0d0

END SUBROUTINE InitFTUVbatches

!> \brief select overlap distribution pass types from overlap structure
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param ODpassesIndex for a given batch this gives the passtype 
!> \param P the overlap 
!> \param nBatches the number of batches
!> \param nPassTypes the number of passtypes 
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
SUBROUTINE SelectODPassTypes(ODpassesIndex,P,nBatches,nPassTypes,IPRINT,LUPRI)
implicit none
Integer       :: nBatches,nPassTypes,IPRINT,LUPRI
Integer       :: ODpassesIndex(nBatches)
Type(Overlap) :: P(nBatches)

Integer       :: I,J,nPrim,numAng1,numAng2
Character(len=60),allocatable  :: UniqeIdentfiers(:)
Integer       :: angmomIdentifier1,angmomIdentifier2,CCidentifier1,CCidentifier2
Logical       :: spher1,spher2,type_hermite_single,unique,sameAO,spherE

ALLOCATE(UniqeIdentfiers(nBatches))

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
     type_hermite_single = P(I)%type_hermite_single
     spherE        = P(I)%sphericalEcoeff
     sameAO        = P(I)%sameAO
     
     IF (nPrim.GT.999999999) THEN
        CALL LSQUIT('nPrim>999999999 in SelectODPassTypes',lupri)
     ELSE IF (angmomIdentifier1.GT.999999999) THEN
        CALL LSQUIT('angmomIdentifier1>999999999 in SelectODPassTypes',lupri)
     ELSE IF (angmomIdentifier2.GT.999999999) THEN
        CALL LSQUIT('angmomIdentifier2>999999999 in SelectODPassTypes',lupri)
     ELSE IF (CCidentifier1.GT.999999999) THEN
        CALL LSQUIT('CCidentifier1>999999999 in SelectODPassTypes',lupri)
     ELSE IF (CCidentifier2.GT.999999999) THEN
        CALL LSQUIT('CCidentifier2>999999999 in SelectODPassTypes',lupri)
     ENDIF
     WRITE(UniqeIdentfiers(I),'(5I9,5L3)') nPrim,angmomIdentifier2,angmomIdentifier1,&
          & CCidentifier1,CCidentifier2,spher1,spher2,type_hermite_single,spherE,sameAO
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

DEALLOCATE(UniqeIdentfiers)

END SUBROUTINE SelectODPassTypes

!> \brief select overlap distribution pass types from OD Batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param ODpassesIndex for a given batch this gives the passtype 
!> \param nPrim the number of primitives
!> \param ODB the overlap batch
!> \param nbatches the number of batches in ODB
!> \param nPassTypes the number of passtypes 
!> \param maxpasses the number of maximum passes 
!> \param the number of passes for each type
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
!> \param side LHS or RHS specification
SUBROUTINE SelectODPassTypesFromODbatch(ODpassesIndex,nPrim,ODB,nBatches,&
     &nPassTypes,maxpasses,maxpassfortype,INPUT,IPRINT,LUPRI,SIDE)
use memory_handling
implicit none
Integer       :: nBatches,nPassTypes,IPRINT,LUPRI
Integer       :: ODpassesIndex(nBatches),maxpasses,nPrim(nBatches)
!Type(Overlap) :: P(nBatches)
TYPE(ODITEM)        :: ODB
TYPE(IntegralInput) :: Input
Character*(*)       :: side
Integer,pointer :: Maxpassfortype(:)

Integer       :: I,J,numAng1,numAng2
!Integer(long) :: UniqeIdentfiers(nBatches)
Character(len=60)  :: UniqeIdentfiers(nBatches)
Integer       :: angmomIdentifier1,angmomIdentifier2,CCidentifier1,CCidentifier2
Logical       :: spher1,spher2,type_hermite_single,spherE,unique,sameAO
!real(realk)   :: MAXGAB
!Real(realk),pointer :: GAB(:,:),primGAB(:,:) 
!Character(len=80)   :: t1,t2

nPassTypes = 0
DO I=1,nBatches
   !DETERMINE number of Primitives
   CALL GET_NPRIM_FROM_ODBATCH(nPrim(I),ODB%BATCH(I),INPUT,SIDE,lupri)  
  
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

     type_hermite_single = (ODB%BATCH(I)%AO(1)%p%type_hermite.AND.ODB%BATCH(I)%AO(2)%p%type_empty)&
          & .OR. (ODB%BATCH(I)%AO(1)%p%type_empty.AND.ODB%BATCH(I)%AO(2)%p%type_hermite)

     spherE = Input%sphericalEcoeff

     sameAO        = ODB%BATCH(I)%sameAO
     
     IF (nPrim(I).GT.999999999) THEN
        CALL LSQUIT('nPrim>999999999 in SelectODPassTypes',lupri)
     ELSE IF (angmomIdentifier1.GT.999999999) THEN
        CALL LSQUIT('angmomIdentifier1>999999999 in SelectODPassTypes',lupri)
     ELSE IF (angmomIdentifier2.GT.999999999) THEN
        CALL LSQUIT('angmomIdentifier2>999999999 in SelectODPassTypes',lupri)
     ELSE IF (CCidentifier1.GT.999999999) THEN
        CALL LSQUIT('CCidentifier1>999999999 in SelectODPassTypes',lupri)
     ELSE IF (CCidentifier2.GT.999999999) THEN
        CALL LSQUIT('CCidentifier2>999999999 in SelectODPassTypes',lupri)
     ENDIF
     WRITE(UniqeIdentfiers(I),'(5I9,5L3)') nPrim(I),angmomIdentifier2,angmomIdentifier1,&
          & CCidentifier1,CCidentifier2,spher1,spher2,type_hermite_single,spherE,sameAO
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

call mem_alloc(Maxpassfortype,nPassTypes)
!NULLIFY(Maxpassfortype)
!ALLOCATE(Maxpassfortype(nPassTypes))
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

!> \brief determine the number of primitives from OD Batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param nPrim the number of primitives
!> \param ODB the overlap batch
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param side LHS or RHS specification
!> \param lupri The default print unit
SUBROUTINE GET_NPRIM_FROM_ODBATCH(nPrim,ODB,input,side,lupri)
IMPLICIT NONE
INTEGER :: nPrim,lupri
Character*(*)       :: side
TYPE(IntegralInput),target :: Input
TYPE(ODBATCH)       :: ODB
!
integer             :: i1,start2,i2
!Character(len=80)   :: t1,t2
real(realk)         :: maxgab,MAXELM
!Real(realk),pointer :: GAB(:,:)
type(lstensor),pointer :: primGAB2
integer             :: AtomA,atomB,batchA,batchB,Gindex

IF (Side.EQ.'LHS') THEN
!  GAB => INPUT%pGAB_LHS
  primGAB2 => INPUT%LST_pGAB_LHS
  MAXELM=INPUT%PS_MAXELM_RHS
ELSEIF (Side.EQ.'RHS') THEN
!  GAB => INPUT%pGAB_RHS
  primGAB2 => INPUT%LST_pGAB_RHS
  MAXELM=INPUT%PS_MAXELM_LHS
ELSE
  WRITE(*,'(1X,2A)') 'Error in GET_NPRIM_FROM_ODBATCH. Side =',Side
  CALL LSQUIT('Programming error. Wrong Side in GET_NPRIM_FROM_ODBATCH',-1)
ENDIF

nPrim = 0
IF (.NOT.(ODB%AO(1)%p%type_empty .AND.ODB%AO(2)%p%type_empty)) THEN
   IF(INPUT%PS_SCREEN)THEN   
      AtomA = ODB%AO(1)%p%atom
      AtomB = ODB%AO(2)%p%atom
      BatchA = ODB%AO(1)%p%batch
      BatchB = ODB%AO(2)%p%batch
      Gindex = primGAB2%INDEX(atomA,atomB,1,1)
      DO i1=1,ODB%AO(1)%p%nPrimitives
         start2 = 1
         IF (ODB%sameAO) start2 = i1 
         DO i2=start2,ODB%AO(2)%p%nPrimitives
 !           maxGab = getMaxPrimGabOD(primGAB2%LSAO(Gindex)%BATCH(BatchA,BatchB,1,1)%elms,i1,i2,ODB,lupri)
            maxGab = primGAB2%LSAO(Gindex)%BATCH(BatchA,BatchB,1,1)%elms(i1 &
                 &+ (i2-1)*ODB%AO(1)%p%nPrimitives)
            IF( MAXGAB .GT. INPUT%PS_THRESHOLD/MAXELM)THEN
               !Add if maxgab greater than threshold
               nPrim = nPrim + 1
            ENDIF
         ENDDO
      ENDDO
   ELSE
      DO i1=1,ODB%AO(1)%p%nPrimitives
         start2 = 1
         IF (ODB%sameAO) start2 = i1 
         DO i2=start2,ODB%AO(2)%p%nPrimitives
            nPrim = nPrim + 1
         ENDDO
      ENDDO
   ENDIF
ELSE
   nPrim = 1
ENDIF
END SUBROUTINE GET_NPRIM_FROM_ODBATCH

!> \brief initialise the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param Input contains info about the integral evaluation requested
!> \param ODbat The OD-batch information
!> \param side The side ('LHS' or 'RHS') 
!> \param maxpasses the maximum number of passes
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
SUBROUTINE INIT_PASS(PassP,Alloc,Input,ODbat,SIDE,maxpasses,IPRINT,LUPRI)
Implicit none
INTEGER             :: lupri,IPRINT,IELECTRON
TYPE(Overlap)       :: PassP
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
TYPE(ODITEM)        :: ODbat
!TYPE(Integralitem)  :: Integral
Character*(*)       :: side
!TYPE(ODBATCH)       :: ODbat(nbatches)
!
!Character(len=80)   :: t1,t2
!Integer             :: NDERIV
Integer             :: endangmom
Integer             :: maxpasses
Integer             :: lenETUV,ma,mc,mp,na,np
LOGICAL             :: LHS,type_hermite_single 

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
   CALL LSQUIT('Wrong case in SET_INITIALIZED_OVERLAP',lupri)
END SELECT

CALL ALLOCATE_ORBITAL(PassP%Orbital1,np,maxpasses,na,ma,.FALSE.)
CALL ALLOCATE_ORBITAL(PassP%Orbital2,np,maxpasses,na,ma,.FALSE.)
  
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

   type_hermite_single = (ODbat%BATCH(1)%AO(1)%p%type_hermite .AND. ODbat%BATCH(1)%AO(2)%p%type_empty) &
        &.OR.(ODbat%BATCH(1)%AO(1)%p%type_empty .AND. ODbat%BATCH(1)%AO(2)%p%type_hermite)
   CALL ALLOCATE_OVERLAP(PassP,mp,mp,maxpasses,na,.FALSE.,np,1,type_hermite_single,lenETUV*maxPasses,INPUT%setETUVoutside)

 END SUBROUTINE INIT_PASS

!> \brief print a general 3 dim tensor
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param tensor the tensor to be printet
!> \param label printing label
!> \param dim1 the size of dimension 1
!> \param dim2 the size of dimension 2
!> \param dim3 the size of dimension 3
!> \param lupri The default print unit
!> \param label1 character string for dimension 1
!> \param label2 character string for dimension 2
!> \param label3 character string for dimension 3
!> \param option determines which dimension it should print in the innermost loop
SUBROUTINE PrintTensor(Tensor,label,dim1,dim2,dim3,Lupri,Label1,Label2,Label3,option)
integer          :: dim1,dim2,dim3,lupri,option
real(realk),dimension(dim1,dim2,dim3) :: Tensor
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

!> \brief determine the number of cartesian or spherical angular components
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param l1 the angular momentum for orbital 1
!> \param l2 the angular momentum for orbital 2
!> \param ijk1 the number of angular components for orbital 1
!> \param ijk2 the number of angular components for orbital 2
!> \param ijk the total number of angular components
!> \param ijkcart the total number of cartesian angular components
!> \param spherical describes if the spherical components should be determined
!> \param derivorder the order of geometrical derivative
!> \param single if only orbital 1 or 2 is used  
SUBROUTINE GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,spherical,derivOrder,single)
implicit none
Integer :: l1,l2,ijk1,ijk2,ijk,ijkcart,derivOrder
Logical :: spherical,single
!
Integer :: iDer,ijk1der,ijk2der,der1,der2
!
!Consistency check
IF ((derivOrder.GT.0).AND.spherical) CALL LSQUIT('Error in get_ijk. derivOrder>0 and spherical',-1)

! Regular case
IF (derivOrder.EQ.0) THEN
  ijk1 = (l1 + 1)*(l1 + 2)/2
  ijk2 = (l2 + 1)*(l2 + 2)/2
  ijkcart = ijk1*ijk2
  IF (spherical) THEN
   ijk1 = 2*l1 + 1
   ijk2 = 2*l2 + 1
  ENDIF
  ijk = ijk1*ijk2
! Derivative case. Note here that we assume Hermite Guassians following PCCP,2007,9,4771-4779
ELSEIF (derivOrder.GT.0) THEN
  ijk1 = (l1 + 1)*(l1 + 2)/2
  ijk2 = (l2 + 1)*(l2 + 2)/2
  IF (single) THEN
    IF (l1.GE.l2) THEN
      ijk1der = (l1 + 1 + derivOrder)*(l1 + 2 + derivOrder)/2
      ijkcart = ijk1der*ijk2
    ELSE
      ijk2der = (l2 + 1 + derivOrder)*(l2 + 2 + derivOrder)/2
      ijkcart = ijk1*ijk2der
    ENDIF
  ELSE
    ijk = 0
    ijkcart = 0
    DO iDer=0,derivOrder
      der1 = iDer
      der2 = derivOrder - iDer
      ijk1der = (l1 + 1 + der1)*(l1 + 2 + der1)/2
      ijk2der = (l2 + 1 + der2)*(l2 + 2 + der2)/2
      ijkcart = ijkcart + ijk1der*ijk2der
    ENDDO
  ENDIF
  ijk = ijkcart
ELSE
  CALL LSQUIT('Error in GET_IJK. derivOrder<0',-1)
ENDIF
END SUBROUTINE GET_IJK
 
END MODULE Thermite_OD
