!> @file 
!> Contains integral distribution routines that places calculated integrals in the proper output
!> Thermite integral distribution module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE thermite_distribute
  use TYPEDEF
  use READMOLEFILE
  use BuildBasisSet
!  use ODbatches
  use OD_type
  use sphcart_matrices
  use precision
  use lstiming
  use thermite_integrals
  use thermite_OD
  use LSTENSOR_OPERATIONSMOD

TYPE derivativeInfo
Integer :: Atom(4)
!For each nDerivP this returns the directional component (grad x,y,z=1,2,3, hessian xx,xy,xz,yy,yz,zz=1,2,3,4,5,6, etc.)
Integer :: nDeriv
Integer :: derivComp
Integer, pointer :: dirComp(:) 
!For each nDerivP this returns the AO the contribution is added to (1-4 mean A,B,C,D, 1,4 means AD, etc.)
Integer, pointer :: AO(:,:)
Integer :: translate      !0 if no translational symmetry, n to translate other contibutions to n
END TYPE derivativeInfo

CONTAINS
!> \brief new distributePQ to lstensor
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param RES contains the result lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributeJengineLstensor(RES,PQ,QPmat2,dimQ,dimP,Input,Lsoutput,LUPRI,IPRINT)
  implicit none 
  Type(integrand),intent(in)      :: PQ
  Type(IntegralInput),intent(in)  :: Input
  Type(IntegralOutput),intent(inout) :: Lsoutput
  Integer,intent(in)              :: LUPRI,IPRINT
  Integer,intent(in)              :: dimQ,dimP
  REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
  TYPE(lstensor),intent(inout)    :: RES
  !
  Integer :: nAngmomP,nOrbP,nderivP,nMAT
  Integer :: iAngmomP,iOrbitalP,iderivP
  Integer :: iA,iB,sA,sB,sC,sD,SAA,SBB
  Integer :: nContA,nContB,nAngA,nAngB,startA,startB
  Integer :: iOrbP,idmat,iPassP,n1,n2,nA,nB,nC,nD
  logical :: SameLHSaos,dograd,dopermutation
  integer :: AtomA,atomB,derivOrder,AB,BA,dimA,dimB
  real(realk) :: dummy(Input%NDMAT_RHS)
!
  Type(derivativeInfo) :: derivInfo
  Integer     :: iAtom,iDer,passPoffset,batchA,batchB,localA,localB,Dab,Dba
  Real(realk) :: derCont

  
  DUMMY=0.d0
  derivOrder = Input%derivOrder
  dograd = INPUT%DO_GRADIENT
  nMAT=Input%NDMAT_RHS
  nAngmomP = PQ%P%p%nAngmom
  SameLHSaos = INPUT%SameLHSaos
  nderivP=Input%nDerivP
  iA = PQ%P%p%indexAng1(nAngmomP)
  dimA = 0
  DO SAA=1,iA
     dimA = dimA+PQ%P%p%orbital1%nOrbComp(SAA)*PQ%P%p%orbital1%nContracted(SAA)
  ENDDO
  iB = PQ%P%p%indexAng2(nAngmomP)
  dimB = 0
  DO SBB=1,iB
     dimB = dimB+PQ%P%p%orbital2%nOrbComp(SBB)*PQ%P%p%orbital2%nContracted(SBB)
  ENDDO
  IF (derivOrder.GT.0) THEN
     CALL initDerivativeOverlapInfo(derivInfo,PQ,Input)
     IF (IPRINT.GT.50) THEN
        CALL printDerivativeOverlapInfo(derivInfo,LUPRI)
     ENDIF
  ENDIF
  iOrbitalP = 1
  DO iAngmomP=1,nAngmomP
     sA = 0
     DO SAA=2,PQ%P%p%indexAng1(iAngmomP)
        sA = sA+PQ%P%p%orbital1%nOrbComp(SAA-1)*PQ%P%p%orbital1%nContracted(SAA-1)
     ENDDO
     sB = 0
     DO SBB=2,PQ%P%p%indexAng2(iAngmomP)
        sB = sB+PQ%P%p%orbital2%nOrbComp(SBB-1)*PQ%P%p%orbital2%nContracted(SBB-1)
     ENDDO
     nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
     DO iDerivP=1,NDerivP
        DO iPassP = 1, PQ%P%p%nPasses
           CALL getOverlapInfo(nContA,nContB,nAngA,nAngB,startA,startB,localA,localB,atomA,atomB,PQ%P%p,iAngmomP,iPassP)
           dopermutation = SameLHSaos.AND.(startA.NE.startB)
           PassPoffset = (iPassP-1)*nContA*nContB*nAngA*nAngB
           iOrbP=iOrbitalP-1+PASSPOFFSET
           IF (dograd) THEN
              derivInfo%Atom(1)=atomA
              derivInfo%Atom(2)=atomB
              iAtom = derivInfo%Atom(derivInfo%AO(1,iDerivP))
              iDer  = derivInfo%dirComp(iDerivP)
              AB = iAtom
              
              Dab = Input%LST_DLHS%INDEX(atomA,atomB,1,1)
              batchA = PQ%P%p%orbital1%batch(iPassP)
              batchB = PQ%P%p%orbital2%batch(iPassP)
              IF(dopermutation)THEN
                 Dba = Input%LST_DLHS%INDEX(atomB,atomA,1,1)              
                 call distributeJengineLstensorgrad(&
                      & Lsoutput%ResultMat2%LSAO(AB)%BATCH(1,1,1,1)%elms,&
                      & Input%LST_DLHS%LSAO(Dab)%BATCH(batchA,batchB,1,1)%elms,&
                      & Input%LST_DLHS%LSAO(Dba)%BATCH(batchB,batchA,1,1)%elms,&
                      & dimA,dimB,sA,sB,iOrbP,startA,startB,nAngA,nAngB,&
                      & nContA,nContB,QPmat2,dimQ,dimP,nmat,&
                      & dopermutation,lupri,iprint,iDer,iAtom,derivOrder,input)
              ELSE
                 call distributeJengineLstensorgrad(&
                      & Lsoutput%ResultMat2%LSAO(AB)%BATCH(1,1,1,1)%elms,&
                      & Input%LST_DLHS%LSAO(Dab)%BATCH(batchA,batchB,1,1)%elms,&
                      & Input%LST_DLHS%LSAO(Dab)%BATCH(batchA,batchB,1,1)%elms,&
                      & dimA,dimB,sA,sB,iOrbP,startA,startB,nAngA,nAngB,&
                      & nContA,nContB,QPmat2,dimQ,dimP,nmat,&
                      & dopermutation,lupri,iprint,iDer,iAtom,derivOrder,input)
              ENDIF
           ELSE
              batchA = PQ%P%p%orbital1%batch(iPassP)
              batchB = PQ%P%p%orbital2%batch(iPassP)           
              AB = RES%INDEX(atoma,atomb,1,1)
              if(dopermutation)THEN
                 BA = RES%INDEX(atomb,atoma,1,1)
                 call distributeJengineLstensor2(&
                      & Lsoutput%ResultMat2%LSAO(AB)%BATCH(batchA,batchB,1,1)%elms,&
                      & Lsoutput%ResultMat2%LSAO(BA)%BATCH(batchB,batchA,1,1)%elms,&
                      & dimA, dimB,dimA,dimB,sA,sB,iOrbP,startA,startB,nAngA,nAngB,&
                      & nContA,nContB,QPmat2,dimQ,dimP,nmat,&
                      & dopermutation,lupri,iprint,input)
              else
                 call distributeJengineLstensor2(&
                      & Lsoutput%ResultMat2%LSAO(AB)%BATCH(batchA,batchB,1,1)%elms,&
                      & DUMMY,dimA, dimB,1,1,sA,sB,iOrbP,startA,startB,nAngA,nAngB,&
                      & nContA,nContB,QPmat2,dimQ,dimP,nmat,&
                      & dopermutation,lupri,iprint,input)              
              endif
           ENDIF
        ENDDO
        iOrbitalP = iOrbitalP + nOrbP
     ENDDO
  ENDDO
  IF (derivOrder.GT.0) CALL freeDerivativeOverlapInfo(derivInfo)

end SUBROUTINE distributeJengineLstensor

!> \brief new distributePQ to lstensor
!> \author T. Kjaergaard
!> \date 2009-02-05
!> \param Jab the coulomb matrix \latexonly J_{ab} \endlatexonly
!> \param Jba the coulomb matrix \latexonly J_{ba} \endlatexonly
!> \param dimA the dimension of the A batch for Jab
!> \param dimB the dimension of the B batch for Jab
!> \param dimA2 the dimension of the A batch for Jba
!> \param dimB2 the dimension of the B batch for Jba
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param iOrbP the P overlap distribution orbital index
!> \param startA the starting orbital index of the full A center
!> \param startB the starting orbital index of the full B center
!> \param nAngA number of angular momentums for the A batch
!> \param nAngB number of angular momentums for the B batch
!> \param nContA number of contracted functions on the A batch
!> \param nContB number of contracted functions on the B batch 
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param nmat the number of density matrices
!> \param dopermutation flag if the the permutational symmetry should be used
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param input integral input contains integral specification 
subroutine distributeJengineLstensor2(Jab,Jba,dimA,dimB,dimA2,dimB2,sA,sB,iOrbP,startA,startB,&
     &nAngA,nAngB,nContA,nContB,QPmat,dimQ,dimP,nmat,dopermutation,lupri,iprint,input)
  implicit none
  Type(IntegralInput),intent(in)  :: Input
  integer,intent(inout)     :: iOrbP
  integer,intent(in)        :: sA,sB,startA,startB,nAngA,nAngB,nContA,nContB
  integer,intent(in)        :: nmat,lupri,iprint,dimQ,dimP,dimA,dimB,dimA2,dimB2
  real(realk),intent(inout) :: Jab(dimA,dimB,nmat), Jba(dimB2,dimA2,nmat)
  real(realk),intent(in)    :: QPmat(dimQ,dimP)
  logical,intent(in)        :: dopermutation
!
  integer                   :: tmpB,tmpA,iContB,iContA,iB1,iA1,idmat
  real(realk)               :: Value

  DO iContB=1,nContB
   tmpB = sB+(iContB-1)*nAngB
   DO iContA=1,nContA
    tmpA = sA+(iContA-1)*nAngA
    DO iB1=1+tmpB,nAngB+tmpB
     DO iA1=1+tmpA,nAngA+tmpA
      iOrbP=iOrbP+1
      DO idmat=1,nmat
       Value = QPmat(idmat,iOrbP)
       !Regular contribution
       Jab(iA1,iB1,idmat) = Jab(iA1,iB1,idmat) + Value
       IF(dopermutation)THEN
          Jba(iB1,iA1,idmat) = Jba(iB1,iA1,idmat) + Value
       ENDIF
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO

end subroutine distributeJengineLstensor2

!> \brief new distributePQ to lstensor
!> \author S. Reine
!> \date 2009-02-05
!> \param grad the output coulomb gradient 
!> \param Dab the density matrix \latexonly D_{ab} \endlatexonly
!> \param Dba the density matrix \latexonly D_{ba} \endlatexonly
!> \param dimA the dimension of the A batch for Jab
!> \param dimB the dimension of the B batch for Jab
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param iOrbP the P overlap distribution orbital index
!> \param startA the starting orbital index of the full A center
!> \param startB the starting orbital index of the full B center
!> \param nAngA number of angular momentums for the A batch
!> \param nAngB number of angular momentums for the B batch
!> \param nContA number of contracted functions on the A batch
!> \param nContB number of contracted functions on the B batch 
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param nmat the number of density matrices
!> \param dopermutation flag if the the permutational symmetry should be used
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param ider derivative index
!> \param iatom atom index
!> \param derivorder the order of the derivate contribution requested
!> \param input integral input contains integral specification 
subroutine distributeJengineLstensorgrad(grad,Dab,Dba,dimA,dimB,sA,sB,iOrbP,startA,startB,&
     &nAngA,nAngB,nContA,nContB,QPmat,dimQ,dimP,nmat,dopermutation,lupri,iprint,ider,iatom,derivorder,input)
  implicit none
  Type(IntegralInput),intent(in)  :: Input
  integer,intent(inout)     :: iOrbP,dimA,dimB,sA,sB
  integer,intent(in)        :: startA,startB,nAngA,nAngB,nContA,nContB,derivOrder
  integer,intent(in)        :: nmat,lupri,iprint,dimQ,dimP,iDer,iAtom
  real(realk),intent(inout) :: grad(3,nmat)
  real(realk),intent(in)    :: QPmat(dimQ,dimP),Dab(dimA,dimB,nmat),Dba(dimB,dimA,nmat)
  logical,intent(in)        :: dopermutation
!
  integer                   :: tmpB,tmpA,iContB,iContA,iB1,iA1,idmat
  real(realk)               :: Value,Dabfac,Dercont(nmat)

  Dercont = 0.d0
  DO iContB=1,nContB
   tmpB = sB+(iContB-1)*nAngB
   DO iContA=1,nContA
    tmpA = sA+(iContA-1)*nAngA
    DO iB1=1+tmpB,nAngB+tmpB
     DO iA1=1+tmpA,nAngA+tmpA
      iOrbP=iOrbP+1
      DO idmat=1,nmat
       Value = QPmat(idmat,iOrbP)
       !Gradient contribution
       Dabfac = Dab(iA1,iB1,idmat)
       IF(dopermutation) Dabfac = Dabfac + Dba(iB1,iA1,idmat)!Input%Dmat_LHS(iB1,iA1,idmat)
       derCont(idmat) = derCont(idmat) + Value * Dabfac * 0.5d0
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  IF (DerivOrder.GT.0) THEN
      DO idmat=1,nmat
!        grad(iDer,idmat) = grad(iDer,idmat) + derCont(idmat)
         grad(iDer,1) = grad(iDer,1) + derCont(idmat)
      enddo
  ENDIF

end subroutine distributeJengineLstensorgrad

!> \brief new distributePQ to lstensor
!> \author T. Kjaergaard
!> \date 2009-02-05
!> \param RES the output lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat2 matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param input the integral input contains integral specification 
!> \param output the integral output specification 
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributelstensorCS(RES,PQ,QPmat2,dimQ,dimP,Input,Output,LUPRI,IPRINT)
implicit none !QQQQQQQQQQ
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: Output
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
TYPE(lstensor),intent(inout)    :: RES
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP,nMAT
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,sA,sB,sC,sD,SAA,SBB,SCC,SDD
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,nderivQ,iderivQ,ideriv
Integer :: iOrbQ,idmat
Integer :: startA1,startB1
real(realk) :: THRESHOLD,TMP
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP
Integer :: n1,n2,dimA,dimB,nA,nB
Integer :: AtomA,atomB!,atomC,atomD
!Integer :: ABCD,BACD,ABDC,BADC,CDAB,CDBA,DCAB,DCBA,nderiv
Integer :: nderiv
logical :: samePQ,SameLHSaos,dopermutation
integer :: batchA,batchB,AB,BA
real(realk),pointer :: dummy(:)

SameLHSaos = INPUT%SameLHSaos
nAngmomP = PQ%P%p%nAngmom

iA = PQ%P%p%indexAng1(nAngmomP)
dimA = 0
DO SAA=1,iA
   dimA = dimA+PQ%P%p%orbital1%nOrbComp(SAA)*PQ%P%p%orbital1%nContracted(SAA)
ENDDO
iB = PQ%P%p%indexAng2(nAngmomP)
dimB = 0
DO SBB=1,iB
   dimB = dimB+PQ%P%p%orbital2%nOrbComp(SBB)*PQ%P%p%orbital2%nContracted(SBB)
ENDDO
if(.NOT.sameLHSaos)call mem_alloc(dummy,dimA*dimB)
iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   sA = 0
   DO SAA=2,iA
      sA = sA+PQ%P%p%orbital1%nOrbComp(SAA-1)*PQ%P%p%orbital1%nContracted(SAA-1)
   ENDDO
   sB = 0
   DO SBB=2,iB
      sB = sB+PQ%P%p%orbital2%nOrbComp(SBB-1)*PQ%P%p%orbital2%nContracted(SBB-1)
   ENDDO
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   DO iPassP = 1, PQ%P%p%nPasses
      PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
      startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
      startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)
      atomA = PQ%P%p%orbital1%atom(iPassP)
      atomB = PQ%P%p%orbital2%atom(iPassP)
!      if(PQ%P%p%orbital1%type_nucleus)atomA=1
!      if(PQ%P%p%orbital2%type_nucleus)atomB=1
      batchA = PQ%P%p%orbital1%batch(iPassP)
      batchB = PQ%P%p%orbital2%batch(iPassP)
      iOrbP=iOrbitalP-1+PASSPOFFSET
      AB = RES%INDEX(ATOMA,ATOMB,1,1)
!      print*,'atomAB,batchABmstartAB,dimAB',atomA,atomB,batchA,batchB,startA,startB,dimA,dimB
!      print*,'B RES%LSAO(AB)%BATCH(batchA,batchB,1,1)%elms',RES%LSAO(AB)%BATCH(batchA,batchB,1,1)%elms
      IF (sameLHSaos)THEN
         BA = RES%INDEX(ATOMB,ATOMA,1,1)
         CALL distributelstensorCS2(& 
              & RES%LSAO(AB)%BATCH(batchA,batchB,1,1)%elms,& 
              & RES%LSAO(BA)%BATCH(batchB,batchA,1,1)%elms,& 
              & iOrbP,dimA,dimB,nAngA,nContA,nAngB,nContB,sA,sB,&
              &QPmat2,dimQ,dimP,sameLHSaos,lupri)
      ELSE
         CALL distributelstensorCS2(& 
              & RES%LSAO(AB)%BATCH(batchA,batchB,1,1)%elms(:),dummy,& 
              & iOrbP,dimA,dimB,nAngA,nContA,nAngB,nContB,sA,sB,&
              & QPmat2,dimQ,dimP,sameLHSaos,lupri)
      ENDIF
   ENDDO
   iOrbitalP = iOrbitalP + nOrbP
ENDDO
if(.NOT.sameLHSaos)call mem_dealloc(dummy)

END SUBROUTINE distributelstensorCS

!> \brief new distributePQ to lstensor
!> \author T. Kjaergaard
!> \date 2009-02-05
!> \param RES the output lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat2 matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param input the integral input contains integral specification 
!> \param output the integral output specification 
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributelstensorPS(RES3,PQ,QPmat2,dimQ,dimP,Input,Output,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: Output
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
TYPE(lstensor),intent(inout)    :: RES3
!
Integer :: nAngmomP,nOrbP,nderivP,totOrbQ,endOrbP,totOrbP,nMAT
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,sA,sB,sC,sD,SAA,SBB,SCC,SDD
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,nderivQ,iderivQ,ideriv
Integer :: iOrbQ,idmat
Integer :: startA1,startB1
real(realk) :: THRESHOLD,TMP
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP
Integer :: n1,n2,dimA,dimB,nA,nB
Integer :: AtomA,atomB!,atomC,atomD
!Integer :: ABCD,BACD,ABDC,BADC,CDAB,CDBA,DCAB,DCBA,nderiv
Integer :: nderiv
logical :: samePQ,SameLHSaos,dopermutation
integer :: batchA,batchB,AB,BA,cBatchA,cBatchB
real(realk),pointer :: dummy(:)

SameLHSaos = INPUT%SameLHSaos
nAngmomP = PQ%P%p%nAngmom

iOrbitalP = 1
DO iAngmomP=1,nAngmomP
   nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   nContA = PQ%P%p%orbital1%nContracted(iA)
   nContB = PQ%P%p%orbital2%nContracted(iB)
   nAngA = PQ%P%p%orbital1%nOrbComp(iA)
   nAngB = PQ%P%p%orbital2%nOrbComp(iB)
   DO iPassP = 1, PQ%P%p%nPasses
      PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
      atomA = PQ%P%p%orbital1%atom(iPassP)
      atomB = PQ%P%p%orbital2%atom(iPassP)
!      print*,'atomA,atomB',atomA,atomB
      batchA = PQ%P%p%orbital1%batch(iPassP)
      batchB = PQ%P%p%orbital2%batch(iPassP)
      iOrbP=iOrbitalP-1+PASSPOFFSET
      AB = RES3%INDEX(ATOMA,ATOMB,1,1)
!      print*,'RES3%INDEX(ATOMA,ATOMB,1,1)',AB
      cBatchA = Output%ContractedBatchA(ATOMA)%BATCH(batchA)
      cBatchB = Output%ContractedBatchB(ATOMB)%BATCH(batchB)
      dimA = Output%CnPrimitivesA(ATOMA)%BATCH(cBatchA)
      dimB = Output%CnPrimitivesB(ATOMB)%BATCH(cBatchB)
      sA = Output%startBatchA(ATOMA)%BATCH(BatchA)
      sB = Output%startBatchB(ATOMB)%BATCH(BatchB)
!      print*,'batchA,batchB',batchA,batchB
!      print*,'cbatchA,cbatchB',cbatchA,cbatchB
!      print*,'dimA,dimB',dimA,dimB
!      print*,'sA,sB',sA,sB
      IF (sameLHSaos)THEN
         BA = RES3%INDEX(ATOMB,ATOMA,1,1)
!         print*,'RES3%INDEX(ATOMB,ATOMA,1,1)',BA
         CALL distributelstensorPS2(& 
              & RES3%LSAO(AB)%BATCH(cbatchA,cbatchB,1,1)%elms,& 
              & RES3%LSAO(BA)%BATCH(cbatchB,cbatchA,1,1)%elms,& 
              & iOrbP,dimA,dimB,nAngA,nContA,nAngB,nContB,sA,sB,&
              &QPmat2,dimQ,dimP,sameLHSaos,lupri)
      ELSE
         call mem_alloc(dummy,dimA*dimB)
         CALL distributelstensorPS2(& 
              & RES3%LSAO(AB)%BATCH(cbatchA,cbatchB,1,1)%elms(:),dummy,& 
              & iOrbP,dimA,dimB,nAngA,nContA,nAngB,nContB,sA,sB,&
              & QPmat2,dimQ,dimP,sameLHSaos,lupri)
         call mem_dealloc(dummy)
      ENDIF
   ENDDO
   iOrbitalP = iOrbitalP + nOrbP
ENDDO

END SUBROUTINE distributelstensorPS

!> \brief distribute the Cauchy-Schwarz (the diagonal 4center eri) integrals to output lstensor
!> \author T. Kjaergaard
!> \date 2009-02-05
!> \param CS_AB the output \latexonly CS_{ab} \endlatexonly in lstensor format
!> \param CS_BA the output \latexonly CS_{ba} \endlatexonly in lstensor format
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param dimA the dimension of the A batch for Jab
!> \param dimB the dimension of the B batch for Jab
!> \param nAngA number of angular momentums for the A batch
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B batch
!> \param nContB number of contracted functions on the B batch 
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param QPmat4 matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param sameLHSaos flag to specify if the Left hand side AOs are the same
!> \param LUPRI logical unit number for output printing
SUBROUTINE distributelstensorPS2(CS_AB,CS_BA,START_iOrbP,dimA,dimB,nAngA,nContA,nAngB,&
     &nContB,sA,sB,QPmat4,dimQ,dimP,sameLHSaos,lupri)
  IMPLICIT NONE
  INTEGER,intent(in)        :: dimA,dimB,lupri,dimQ,dimP
  REAL(REALK),intent(in)    :: QPmat4(dimQ,dimP) 
  REAL(REALK),intent(inout) :: CS_AB(dimA,dimB)
  REAL(REALK),intent(inout) :: CS_BA(dimB,dimA)
  INTEGER,intent(in)        :: START_iOrbP,sA,nContA,nAngA
  INTEGER,intent(in)        :: sB,nContB,nAngB
  logical                   :: sameLHSaos
  !
  REAL(REALK)    :: Int
  INTEGER        :: iAngA,iAngB,iContA,iContB
  INTEGER        :: iOrbP,iA,iB

!  print*,'dimA',dimA,'dimB',dimB
  iOrbP = START_iOrbP
  DO iContB=1,nContB
   DO iContA=1,nContA
    DO iAngB=1,nAngB
     DO iAngA=1,nAngA
      iOrbP=iOrbP+1
      int = sqrt(QPmat4(iOrbP,iOrbP))
 !     print*,'iContA+sA',iContA,sA,'iContB+sB',iContB,sB
      CS_AB(iContA+sA,iContB+sB) = MAX(CS_AB(iContA+sA,iContB+sB),int)  
      IF (sameLHSaos)then
         CS_BA(iContB+sB,iContA+sA) = MAX(CS_BA(iContB+sB,iContA+sA),int)  
      endif
     ENDDO
    ENDDO
   ENDDO
  ENDDO

END SUBROUTINE distributelstensorPS2

!> \brief distribute the Cauchy-Schwarz (the diagonal 4center eri) integrals to output lstensor
!> \author T. Kjaergaard
!> \date 2009-02-05
!> \param CS_AB the output \latexonly CS_{ab} \endlatexonly in lstensor format
!> \param CS_BA the output \latexonly CS_{ba} \endlatexonly in lstensor format
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param dimA the dimension of the A batch for Jab
!> \param dimB the dimension of the B batch for Jab
!> \param nAngA number of angular momentums for the A batch
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B batch
!> \param nContB number of contracted functions on the B batch 
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param QPmat4 matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param sameLHSaos flag to specify if the Left hand side AOs are the same
!> \param LUPRI logical unit number for output printing
SUBROUTINE distributelstensorCS2(CS_AB,CS_BA,START_iOrbP,dimA,dimB,nAngA,nContA,nAngB,&
     &nContB,sA,sB,QPmat4,dimQ,dimP,sameLHSaos,lupri)
  IMPLICIT NONE
  INTEGER,intent(in)        :: dimA,dimB,lupri,dimQ,dimP
  REAL(REALK),intent(in)    :: QPmat4(dimQ,dimP) 
  REAL(REALK),intent(inout) :: CS_AB(dimA,dimB)
  REAL(REALK),intent(inout) :: CS_BA(dimB,dimA)
  INTEGER,intent(in)        :: START_iOrbP,sA,nContA,nAngA
  INTEGER,intent(in)        :: sB,nContB,nAngB
  logical                   :: sameLHSaos
  !
  REAL(REALK)    :: Int
  INTEGER        :: iAngA,iAngB,iContA,iContB
  INTEGER        :: iOrbP,iA,iB

  iOrbP = START_iOrbP
  DO iContB=1,nContB
   DO iContA=1,nContA
    DO iAngB=1,nAngB
     iB=iAngB+(iContB-1)*nAngB
     DO iAngA=1,nAngA
      iA=iAngA+(iContA-1)*nAngA
      iOrbP=iOrbP+1
      int = sqrt(QPmat4(iOrbP,iOrbP))
      CS_AB(iA+sA,iB+sB) = int  
!      print*,'CS_AB                    (',iA+sA,',',iB+sB,')=',CS_AB(iA+sA,iB+sB)
!      print*,'should correspond to FULL(',iA+startA-1,',',iB+startB-1,')'
      IF (sameLHSaos)then
         CS_BA(iB+sB,iA+sA) = int  
!         print*,'CS_BA(',iB+sB,',',iA+sA,')=',CS_BA(iB+sB,iA+sA)
!         print*,'should correspond to FULL(',iB+startB-1,',',iA+startA-1,')'
      endif
     ENDDO
    ENDDO
   ENDDO
  ENDDO

END SUBROUTINE distributelstensorCS2

!> \brief print calculated multipole moments to file
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param PQ contain info about the overlap distributions
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param nMAT = dimQ = number of cartesian MM components
!> \param nSPHMAT number of spherical MM components
!> \param Input contain info about the requested integral 
!> \param Output2 contain info about whether or not to use the buffer to store integral before dumping to disk
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE printMMtoFile(PQ,QPmat,dimQ,dimP,nMAT,nSPHMAT,INPUT,OUTPUT2,LUPRI,IPRINT)
!use mm_interface_mod
  implicit none
  Type(integrand),intent(in)      :: PQ
  Type(IntegralInput),intent(inout)  :: Input
  Type(IntegralOutput),intent(inout) :: Output2
  Integer,intent(in)              :: LUPRI,IPRINT
  Integer,intent(in)              :: dimQ,dimP,nMAT,nSPHMAT
  REAL(REALK),intent(in)          :: QPMAT(dimQ,dimP)
  !
  Integer :: nAngmomP,nOrbP,nderivP,endOrbP,nderiv
  Integer :: iAngmomP,iOrbitalP,iderivP
  Integer :: iA,iB,I,J
  Integer :: nContA,nContB,nAngA,nAngB,startA,startB
  integer :: iPassP,iOrbP,iContP,iContB,iContA,iAngB,iAngA,imat,iderivQ,m
  integer :: PASSPOFFSET,iAng1A,iAng1B,ml,II_START,nderivQ,IBATCH
  INTEGER,allocatable :: INDDSTR(:,:,:,:)
  REAL(REALK)          :: THRESHOLD
  REAL(REALK) :: TEMP(dimP,nMAT),TEMP2(dimP,nSPHMAT)
  logical             :: AddBatchCont(dimP), AddBatch
  integer             :: IBATCHCont(dimP), NBATCH,INPUTSTARTA,INPUTSTARTB

  INPUTSTARTA = INPUT%MMstartA
  INPUTSTARTB = INPUT%MMstartB

  AddBatch = .FALSE.
  THRESHOLD = INPUT%MM_SCREENTHR!1.0D-15     
  THRESHOLD = INPUT%CS_THRESHOLD*1.d-1
  CALL mm_set_screen_threshold(THRESHOLD)
  II_START = INPUT%MMunique_ID1 !unique identifier ()start with 0 
!  IBATCH = INPUT%MMunique_ID2+1 !BATCH INDEX starts with 0 
  NBATCH = INPUT%MMunique_ID2 !BATCH INDEX starts with 0 
  !nMAT =input%nderivQ
  nderiv = INPUT%DERIVORDER
  !nSPHMAT=(nderiv+1)*(nderiv+1)
  nAngmomP = PQ%P%p%nAngmom
  nderivP=Input%nDerivP
  nderivQ=Input%nDerivQ
!  WRITE(lupri,*)'THE QP mat dimQ',dimQ,'dimP',dimP
!  CALL OUTPUT(QPmat,1,dimQ,1,dimP,dimQ,dimP,1,lupri)

  DO I = 1,dimP
     DO J = 1,nMAT
        TEMP(I,J) = QPmat(J,I)
     ENDDO
  ENDDO
  CALL LS_DZERO(TEMP2,dimP*nSPHMAT)
  CALL SPHERICAL_TRANSFORMATION2(TEMP,TEMP2,dimP,nMAT,nSPHMAT,nderiv,iprint,lupri)
  !
  iOrbitalP = 1
  DO iDerivP=1,NDerivP
     DO iAngmomP=1,nAngmomP
        nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
        iAng1A = PQ%P%p%indexAng1(iAngmomP)
        iAng1B = PQ%P%p%indexAng2(iAngmomP)
        nContA = PQ%P%p%orbital1%nContracted(iAng1A)
        nContB = PQ%P%p%orbital2%nContracted(iAng1B)
        nAngA = PQ%P%p%orbital1%nOrbComp(iAng1A)
        nAngB = PQ%P%p%orbital2%nOrbComp(iAng1B)
        !call mem_alloc_nonp(INDDSTR,nContA,nContB,nAngA,nAngB)
        ALLOCATE(INDDSTR(nContA,nContB,nAngA,nAngB))
        INDDSTR = 0
        DO iPassP = 1, PQ%P%p%nPasses
           AddBatchCont = .FALSE.
           endOrbP = iOrbitalP+nOrbP-1
           !
           PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
           startA  = INPUTSTARTA+PQ%P%p%orbital1%startOrbital(iAng1A,iPassP)
           startB  = INPUTSTARTB+PQ%P%p%orbital2%startOrbital(iAng1B,iPassP)     
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
                       imat=0
                       DO iderivQ =0,nderiv!-1
                          DO m = 0,iderivQ
                             ! empty buffer if necessary
                             if (OUTPUT2%USEBUFMM .and. (OUTPUT2%IBUFI.GE. OUTPUT2%MMBUFLEN-1)) then
                                CALL LS_EMPTYIBUF(OUTPUT2,OUTPUT2%IBUF,INPUT%lu_mmdata)
                                CALL LS_EMPTYRBUF(OUTPUT2,OUTPUT2%RBUF,INPUT%lu_mmdatr)
                             endif
                             IF(m.eq.0) THEN
                                ml = 0
                                imat=imat+1
                                IF(ABS(TEMP2(iOrbP,imat)).GT.THRESHOLD)THEN
                                   IF(iB .LE. iA .OR. startA.LT.startB)THEN
                                      AddBatch = .TRUE.
                                      IF (.NOT.AddBatchCont(iContP)) THEN
                                         AddBatchCont(iContP) = .TRUE.
                                         NBATCH = NBATCH + 1
                                         IBATCHCont(iContP) = NBATCH
                                      ENDIF
                                      IBATCH = IBATCHCont(iContP)
                                      IF ( INDDSTR( ICONTA, ICONTB, iAngA,iAngB) .EQ. 0 )THEN
                                         II_START = II_START + 1
                                         INDDSTR( ICONTA, ICONTB, iAngA,iAngB) = II_START
                                      ENDIF
                                      ! fill buffer
                                      if (OUTPUT2%USEBUFMM)then
                                         CALL LS_FILLBUFFER(OUTPUT2,iderivQ,ml,iB,iA,&
                                              &PQ%P%p%orbital2%ANGMOM(iAng1B),PQ%P%p%orbital1%ANGMOM(iAng1A),&
                                              &IBATCH,INDDSTR( ICONTA, ICONTB, iAngA,iAngB),PQ%P%p%ODextent,&
                                              &PQ%P%p%ODcenter(1),PQ%P%p%ODcenter(2),PQ%P%p%ODcenter(3),&
                                              &-TEMP2(iOrbP,imat),0.0D0)
                                      else
                                         WRITE(INPUT%lu_mmdata)iderivQ,ml,iB,iA,&
                                              &PQ%P%p%orbital2%ANGMOM(iAng1B),&
                                              &PQ%P%p%orbital1%ANGMOM(iAng1A),IBATCH,INDDSTR( ICONTA, ICONTB, iAngA,iAngB),&
                                              &PQ%P%p%ODextent,PQ%P%p%ODcenter(1),&
                                              &PQ%P%p%ODcenter(2),PQ%P%p%ODcenter(3),&
                                              &-TEMP2(iOrbP,imat),0.D0!4*INPUT%DMAT_LHS(iA,iB,1)
                                      endif
                                   ENDIF
                                ENDIF
                             ELSE
                                ml = m
                                imat=imat+1
                                IF(ABS(TEMP2(iOrbP,imat)).GT.THRESHOLD)THEN
                                   IF(iB .LE. iA .OR. startA.LT.startB)THEN
                                      AddBatch = .TRUE.
                                      IF (.NOT.AddBatchCont(iContP)) THEN
                                         AddBatchCont(iContP) = .TRUE.
                                         NBATCH = NBATCH + 1
                                         IBATCHCont(iContP) = NBATCH
                                      ENDIF
                                      IBATCH = IBATCHCont(iContP)
                                      IF ( INDDSTR( ICONTA, ICONTB, iAngA,iAngB) .EQ. 0 )THEN
                                         II_START = II_START + 1
                                         INDDSTR( ICONTA, ICONTB, iAngA,iAngB) = II_START
                                      ENDIF
                                      ! fill buffer
                                      if (OUTPUT2%USEBUFMM)then
                                         CALL LS_FILLBUFFER(OUTPUT2,iderivQ,ml,iB,iA,&
                                              &PQ%P%p%orbital2%ANGMOM(iAng1B),PQ%P%p%orbital1%ANGMOM(iAng1A),&
                                              &IBATCH,INDDSTR( ICONTA, ICONTB, iAngA,iAngB),PQ%P%p%ODextent,& 
                                              &PQ%P%p%ODcenter(1),PQ%P%p%ODcenter(2),PQ%P%p%ODcenter(3),& 
                                              &-TEMP2(iOrbP,imat),0.0D0)
                                      else
                                         WRITE(INPUT%lu_mmdata)iderivQ,ml,iB,iA,&
                                              &PQ%P%p%orbital2%ANGMOM(iAng1B),&
                                              &PQ%P%p%orbital1%ANGMOM(iAng1A),IBATCH,INDDSTR( ICONTA, ICONTB, iAngA,iAngB),&
                                              &PQ%P%p%ODextent,PQ%P%p%ODcenter(1),&
                                              &PQ%P%p%ODcenter(2),PQ%P%p%ODcenter(3),&
                                              &-TEMP2(iOrbP,imat),0.D0!4*INPUT%DMAT_LHS(iA,iB,1)
                                      endif
                                   ENDIF
                                ENDIF
                                ! empty buffer if necessary
                                if (OUTPUT2%USEBUFMM .and. (OUTPUT2%IBUFI.GE. OUTPUT2%MMBUFLEN-1)) then
                                   CALL LS_EMPTYIBUF(OUTPUT2,OUTPUT2%IBUF,INPUT%lu_mmdata)
                                   CALL LS_EMPTYRBUF(OUTPUT2,OUTPUT2%RBUF,INPUT%lu_mmdatr)
                                endif
                                ml = -m
                                imat=imat+1
                                IF(ABS(TEMP2(iOrbP,imat)).GT.THRESHOLD)THEN
                                   IF(iB .LE. iA .OR. startA.LT.startB)THEN
                                      AddBatch = .TRUE.
                                      IF (.NOT.AddBatchCont(iContP)) THEN
                                         AddBatchCont(iContP) = .TRUE.
                                         NBATCH = NBATCH + 1
                                         IBATCHCont(iContP) = NBATCH
                                      ENDIF
                                      IBATCH = IBATCHCont(iContP)
                                      IF ( INDDSTR( ICONTA, ICONTB, iAngA,iAngB) .EQ. 0 )THEN
                                         II_START = II_START + 1
                                         INDDSTR( ICONTA, ICONTB, iAngA,iAngB) = II_START
                                      ENDIF
                                      ! fill buffer
                                      if (OUTPUT2%USEBUFMM)then
                                         CALL LS_FILLBUFFER(OUTPUT2,iderivQ,ml,iB,iA,&
                                              &PQ%P%p%orbital2%ANGMOM(iAng1B),PQ%P%p%orbital1%ANGMOM(iAng1A),&
                                              &IBATCH,INDDSTR( ICONTA, ICONTB, iAngA,iAngB),PQ%P%p%ODextent,& 
                                              &PQ%P%p%ODcenter(1),PQ%P%p%ODcenter(2),PQ%P%p%ODcenter(3),& 
                                              &-TEMP2(iOrbP,imat),0.0D0)
                                      else
                                         WRITE(INPUT%lu_mmdata)iderivQ,ml,iB,iA,&
                                              &PQ%P%p%orbital2%ANGMOM(iAng1B),&
                                              &PQ%P%p%orbital1%ANGMOM(iAng1A),IBATCH,INDDSTR( ICONTA, ICONTB, iAngA,iAngB),&
                                              &PQ%P%p%ODextent,PQ%P%p%ODcenter(1),&
                                              &PQ%P%p%ODcenter(2),PQ%P%p%ODcenter(3),&
                                              &-TEMP2(iOrbP,imat),0.D0!4*INPUT%DMAT_LHS(iA,iB,1)
                                      endif
                                   ENDIF
                                ENDIF
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        !call mem_dealloc_nonp(INDDSTR)
        DEALLOCATE(INDDSTR)
        iOrbitalP = iOrbitalP + nOrbP
     ENDDO
  ENDDO
  INPUT%MMunique_ID1 = II_START
  IF (AddBatch) INPUT%MMunique_ID2 = NBATCH

END SUBROUTINE printMMtoFile

!> \brief Get some key information about overlap type for given angular index and pass number
!> \author S. Reine
!> \date 2010-03-23
!> \param nContA The number of contracted functions for orbital 1
!> \param nContB The number of contracted functions for orbital 2
!> \param nAngA The number of angular components for orbital 1
!> \param nAngB The number of angular components for orbital 2
!> \param startA The starting orbital number for orbital 1
!> \param startB The starting orbital number for orbital 2
!> \param atomA The atom-center index of orbital 1
!> \param atomB The atom-center index of orbital 2
!> \param P The overlap type
!> \param iAngmom The angular index (allways 1 for .NOFAMILY)
!> \param iPass The pass number
SUBROUTINE getOverlapInfo(nContA,nContB,nAngA,nAngB,startA,startB,localA,localB,atomA,atomB,P,iAngmom,iPass)
implicit none
TYPE(overlap),intent(IN) :: P
Integer,intent(OUT)      :: nContA,nContB,nAngA,nAngB,startA,startB,localA,localB,atomA,atomB
Integer,intent(IN)       :: iAngmom,iPass
!
CALL getOrbitalInfo(nContA,nAngA,startA,localA,atomA,P%orbital1,iPass,P%indexAng1(iAngmom))
CALL getOrbitalInfo(nContB,nAngB,startB,localB,atomB,P%orbital2,iPass,P%indexAng2(iAngmom))

END SUBROUTINE getOverlapInfo

!> \brief Get some key information about an orbital type for given pass number
!> \author S. Reine
!> \date 2010-03-23
!> \param nCont The number of contracted functions
!> \param nAng The number of angular components
!> \param start The starting orbital number
!> \param local The local starting orbital number
!> \param atom The atom-center index
!> \param orb The orbital type
!> \param iPass The pass number
!> \param iAng The angular index (allways 1 for .NOFAMILY)
SUBROUTINE getOrbitalInfo(nCont,nAng,start,local,atom,orb,iPass,iAng)
implicit none
integer,intent(OUT)      :: nCont,nAng,start,local,atom
TYPE(orbital),intent(IN) :: orb
Integer,intent(IN)       :: iPass,iAng
!
nCont = orb%nContracted(iAng)
nAng  = orb%nOrbComp(iAng)
start = orb%startOrbital(iAng,iPass)
local = orb%startLocOrb(iAng)
atom  =orb%atom(iPass)
!
END SUBROUTINE getOrbitalInfo

!> \brief Initialize the derivative info
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param PQ The integrand
!> \param Input The integral input
SUBROUTINE initDerivativeOverlapInfo(derivInfo,PQ,Input)
implicit none
Type(integrand),intent(in)          :: PQ
Type(IntegralInput),intent(in)      :: Input
Type(derivativeInfo), intent(INOUT) :: derivInfo
!
Logical :: emptyA,emptyB,emptyC,emptyD
Integer :: nP,nQ,n,derP,derQ,nDer,iDer,iDeriv
!
derivInfo%atom = 0
!
!Set up translation (if no translation translate=0)
IF (PQ%P%p%type_nucleus.OR.PQ%Q%p%type_nucleus) THEN
!  If one side is equal to nuclei we exploit translational symmetry
  IF (PQ%P%p%type_nucleus) THEN
    IF (PQ%P%p%orbital1%type_nucleus) THEN
      derivInfo%translate = 1
    ELSE
      derivInfo%translate = 2
    ENDIF
  ELSE
    IF (PQ%Q%p%orbital1%type_nucleus) THEN
      derivInfo%translate = 3
    ELSE
      derivInfo%translate = 4
    ENDIF
  ENDIF
ELSE
  derivInfo%translate = 0
ENDIF
!
emptyA=PQ%P%p%orbital1%type_empty
emptyB=PQ%P%p%orbital2%type_empty
emptyC=PQ%Q%p%orbital1%type_empty
emptyD=PQ%Q%p%orbital2%type_empty

nP = Input%nDerivP
nQ = Input%nDerivQ
n  = nP*nQ
derP = Input%derOrderP
derQ = Input%derOrderQ
nDer = derP + derQ
derivInfo%nDeriv = nDer
derivInfo%derivComp = n
call mem_alloc(derivInfo%dirComp,n)
call mem_alloc(derivInfo%AO,nDer,n)

iDeriv = 0
IF (emptyA.AND.emptyB) THEN
  CALL setDerivativeComponents(derivInfo,nDer,iDeriv,0,0,derQ,emptyC,emptyD)
ELSEIF (emptyB) THEN
  CALL setDerivativeComponents(derivInfo,nDer,iDeriv,derP,0,derQ,emptyC,emptyD)
ELSEIF (emptyA) THEN
  CALL setDerivativeComponents(derivInfo,nDer,iDeriv,0,derP,derQ,emptyC,emptyD)
ELSE
  DO iDer=0,derP
    CALL setDerivativeComponents(derivInfo,nDer,iDeriv,derP-iDer,iDer,derQ,emptyC,emptyD)
  ENDDO
ENDIF
END SUBROUTINE initDerivativeOverlapInfo

!> \brief Print the derivative info
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param LUPRI Default output unit
SUBROUTINE printDerivativeOverlapInfo(derivInfo,LUPRI)
implicit none
TYPE(derivativeInfo), intent(IN) :: derivInfo
Integer,intent(IN)               :: LUPRI
!
Integer :: iDer,nDeriv,derComp

CALL LSHEADER(LUPRI,'Derivative Info')
IF (derivInfo%atom(1).EQ.0) THEN
  write(LUPRI,'(3X,A)') 'Atomic numbers not set'
ELSE
  write(LUPRI,'(3X,A,4I5)') 'Atomic numbers: ',derivInfo%atom
ENDIF

IF (derivInfo%translate.EQ.0) THEN
  write(LUPRI,'(3X,A)') 'No translational symmetry employed'
ELSE
  write(LUPRI,'(3X,A,I2)') 'Translation symmetry employed. Contributions subtracted to AO index number',derivInfo%translate
ENDIF

nDeriv = derivInfo%nDeriv
derComp = derivInfo%derivComp
WRITE(LUPRI,'(3X,A,I1,A,I4)') 'Derivative order ',nDeriv,', number of derivative components', derComp
WRITE(LUPRI,'(5X,A)') 'Derivative index   directional component   AO derivatives'
DO iDer=1,derComp
  WRITE(LUPRI,'(5X,I10,I20,15X,9I2)') iDer,derivInfo%dirComp(iDer),derivInfo%AO(:,iDer)
ENDDO

END SUBROUTINE printDerivativeOverlapInfo

!> \brief Sets up the derivative component information
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param nDer The derivative order
!> \param iDeriv The derivative component
!> \param derA The derivative order of orbital A
!> \param derB The derivative order of orbital B
!> \param derQ The derivative order of RHS overlap
!> \param emptyC Specifies if orbital C is empty
!> \param emptyD Specifies if orbital D is empty
SUBROUTINE setDerivativeComponents(derivInfo,nDer,iDeriv,derA,derB,derQ,emptyC,emptyD)
implicit none
TYPE(derivativeInfo), intent(INOUT) :: derivInfo
Integer,intent(INOUT)               :: iDeriv
Integer,intent(IN)                  :: nDer,derA,derB,derQ
Logical,intent(IN)                  :: emptyC,emptyD
!
Integer :: xA,yA,zA,xB,yB,zB
DO xA=derA,0,-1
  DO yA=derA-xA,0,-1
    zA=derA-xA-yA
    DO xB=derB,0,-1
      DO yB=derB-xB,0,-1
        zB=derB-xB-yB
        CALL setDerivativeComp1(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,derQ,emptyC,emptyD)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE setDerivativeComponents

!> \brief Sets up the derivative component information - continued 1
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param nDer The derivative order
!> \param iDeriv The derivative component
!> \param derA The derivative order of orbital A
!> \param xA The x-directional cartesian component derivative order of orbital A
!> \param yA The y-directional cartesian component derivative order of orbital A
!> \param zA The z-directional cartesian component derivative order of orbital A
!> \param derB The derivative order of orbital B
!> \param xB The x-directional cartesian component derivative order of orbital B
!> \param yB The y-directional cartesian component derivative order of orbital B
!> \param zB The z-directional cartesian component derivative order of orbital B
!> \param derQ The derivative order of RHS overlap
!> \param emptyC Specifies if orbital C is empty
!> \param emptyD Specifies if orbital D is empty
SUBROUTINE setDerivativeComp1(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,derQ,emptyC,emptyD)
implicit none
TYPE(derivativeInfo), intent(INOUT) :: derivInfo
Integer,intent(INOUT)               :: iDeriv
Integer,intent(IN)                  :: nDer,derA,xA,yA,zA,derB,xB,yB,zB,derQ
Logical,intent(IN)                  :: emptyC,emptyD
!
Integer :: iDer
!
IF (emptyC.AND.emptyD) THEN
  CALL setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,0,0)
ELSEIF (emptyD) THEN
  CALL setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,derQ,0)
ELSEIF (emptyC) THEN
  CALL setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,0,derQ)
ELSE
  DO iDer=0,derQ
    CALL setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,iDer,derQ-iDer)
  ENDDO
ENDIF
END SUBROUTINE setDerivativeComp1

!> \brief Sets up the derivative component information - continued 2
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param nDer The derivative order
!> \param iDeriv The derivative component
!> \param derA The derivative order of orbital A
!> \param xA The x-directional cartesian component derivative order of orbital A
!> \param yA The y-directional cartesian component derivative order of orbital A
!> \param zA The z-directional cartesian component derivative order of orbital A
!> \param derB The derivative order of orbital B
!> \param xB The x-directional cartesian component derivative order of orbital B
!> \param yB The y-directional cartesian component derivative order of orbital B
!> \param zB The z-directional cartesian component derivative order of orbital B
!> \param derC The derivative order of orbital C
!> \param derD The derivative order of orbital D
SUBROUTINE setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,derC,derD)
implicit none
TYPE(derivativeInfo), intent(INOUT) :: derivInfo
Integer,intent(INOUT)               :: iDeriv
Integer,intent(IN)                  :: nDer,derA,xA,yA,zA,derB,xB,yB,zB,derC,derD
!
Integer :: xC,yC,zC,xD,yD,zD,direction(0:nDer,0:nDer,0:nDer),dir,iX,iY,iZ,nX,nY,nZ,iAO
!
dir = 0
DO iX=nDer,0,-1
  DO iY=nDer-iX,0,-1
    iZ=nDer-iX-iY
    dir = dir + 1
    direction(iX,iY,iZ) = dir
  ENDDO
ENDDO
!
DO xC=derC,0,-1
  DO yC=derC-xC,0,-1
    zC=derC-xC-yC
    DO xD=derD,0,-1
      DO yD=derD-xD,0,-1
        zD=derD-xD-yD
        iDeriv = iDeriv + 1
        nX = xA+xB+xC+xD
        nY = yA+yB+yC+yD
        nZ = zA+zB+zC+zD
        derivInfo%dirComp(iDeriv) = direction(nX,nY,nZ)
        iAO = 0
        CALL setDerivAO(derivInfo%AO,iAO,xA,iDeriv,1)
        CALL setDerivAO(derivInfo%AO,iAO,xB,iDeriv,2)
        CALL setDerivAO(derivInfo%AO,iAO,xC,iDeriv,3)
        CALL setDerivAO(derivInfo%AO,iAO,xD,iDeriv,4)
        CALL setDerivAO(derivInfo%AO,iAO,yA,iDeriv,1)
        CALL setDerivAO(derivInfo%AO,iAO,yB,iDeriv,2)
        CALL setDerivAO(derivInfo%AO,iAO,yC,iDeriv,3)
        CALL setDerivAO(derivInfo%AO,iAO,yD,iDeriv,4)
        CALL setDerivAO(derivInfo%AO,iAO,zA,iDeriv,1)
        CALL setDerivAO(derivInfo%AO,iAO,zB,iDeriv,2)
        CALL setDerivAO(derivInfo%AO,iAO,zC,iDeriv,3)
        CALL setDerivAO(derivInfo%AO,iAO,zD,iDeriv,4)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE setDerivComp2

!> \brief Sets up information about the directional cartesian components for given deriviative index
!> \author S. Reine
!> \date 2010-03-23
!> \param AO The AO index (1-4) of each directional derivative cartesian component and deriviative index
!> \param iAO The directional derivative component (1 for gradients, 2 for Hessians, etc.)
!> \param order The cartesian component derivative order
!> \param iDeriv The (full) derivative index
!> \param AOindex The AO index
SUBROUTINE setDerivAO(AO,iAO,order,iDeriv,AOindex)
implicit none
Integer,pointer       :: AO(:,:)
Integer,intent(INOUT) :: iAO
Integer,intent(IN)    :: order,iDeriv,AOindex
!
Integer :: i
DO i=1,order
  iAO=iAO+1
  AO(iAO,iDeriv) = AOindex
ENDDO
END SUBROUTINE setDerivAO

!> \brief Frees the derivative info
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
SUBROUTINE freeDerivativeOverlapInfo(derivInfo)
implicit none
TYPE(derivativeInfo), intent(INOUT) :: derivInfo
derivInfo%translate = 0
derivInfo%nDeriv = 0
derivInfo%derivComp = 0
call mem_dealloc(derivInfo%dirComp)
call mem_dealloc(derivInfo%AO)
END SUBROUTINE freeDerivativeOverlapInfo

!> \brief new distributePQ to lstensor
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param RES contains the result lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param Output contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributeGrad(RES,PQ,QPmat2,dimQ,dimP,Input,LUPRI,IPRINT)
implicit none !QQQQQQQQQQ
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
TYPE(lstensor),intent(inout)    :: RES
!
Integer :: nAngmomP,nOrbP,nderivP,nMAT
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,sA,sB,sC,sD,SAA,SBB,SCC,SDD
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,nderivQ,iderivQ,ideriv
Integer :: iOrbitalQ,endAngmomQ,iAngmomQ,nOrbQ
Integer :: nContC,nContD,nAngC,nAngD,startC,startD
Integer :: iC,iD,iOrbQ
!Integer :: startA1,startB1,startC1,startD1
real(realk) :: THRESHOLD!,TMP
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP,iPassQ
Integer :: dimA,dimB,dimC,dimD
Integer :: AtomA,atomB,atomC,atomD
!Integer :: ABCD,BACD,ABDC,BADC,CDAB,CDBA,DCAB,DCBA
Integer :: nderiv
logical :: samePQ,SameRHSaos,SameLHSaos,SameODs,dograd,translate
integer :: batchA,batchB,batchC,batchD,derivorder,iatom,ider
integer :: localA,localB,localC,localD,D_AB,D_BA,D_CA,D_AC,itrans
!real(realk),pointer :: DUMMY(:)
Type(derivativeInfo) :: derivInfo

derivOrder = Input%derivOrder
IF (derivOrder.GT.0) THEN
  CALL initDerivativeOverlapInfo(derivInfo,PQ,Input)
  IF (IPRINT.GT.50) THEN
    CALL printDerivativeOverlapInfo(derivInfo,LUPRI)
  ENDIF
ENDIF
dograd = INPUT%DO_GRADIENT
nMAT=Input%NDMAT_LHS
nAngmomP = PQ%P%p%nAngmom
THRESHOLD = INPUT%CS_THRESHOLD
SamePQ = PQ%samePQ
SameLHSaos = INPUT%SameLHSaos
SameRHSaos = INPUT%SameRHSaos
SameODs = INPUT%SameODs
nderivQ=input%nderivQ
nderivP=Input%nDerivP
nderiv = nderivP*nderivQ

iA = PQ%P%p%indexAng1(nAngmomP)
dimA = 0
DO SAA=1,iA
   dimA = dimA+PQ%P%p%orbital1%nOrbComp(SAA)*PQ%P%p%orbital1%nContracted(SAA)
ENDDO
iB = PQ%P%p%indexAng2(nAngmomP)
dimB = 0
DO SBB=1,iB
   dimB = dimB+PQ%P%p%orbital2%nOrbComp(SBB)*PQ%P%p%orbital2%nContracted(SBB)
ENDDO
iC = PQ%Q%p%indexAng1(PQ%Q%p%nAngmom)
dimC = 0
DO SCC=1,iC
   dimC = dimC+PQ%Q%p%orbital1%nOrbComp(SCC)*PQ%Q%p%orbital1%nContracted(SCC)
ENDDO
iD = PQ%Q%p%indexAng2(PQ%Q%p%nAngmom)
dimD = 0
DO SDD=1,iD
   dimD = dimD+PQ%Q%p%orbital2%nOrbComp(SDD)*PQ%Q%p%orbital2%nContracted(SDD)
ENDDO

iOrbitalP = 1
DO iAngmomP=1,nAngmomP
 sA = 0
 DO SAA=2,PQ%P%p%indexAng1(iAngmomP)
    sA = sA+PQ%P%p%orbital1%nOrbComp(SAA-1)*PQ%P%p%orbital1%nContracted(SAA-1)
 ENDDO
 
 sB = 0
 DO SBB=2,PQ%P%p%indexAng2(iAngmomP)
    sB = sB+PQ%P%p%orbital2%nOrbComp(SBB-1)*PQ%P%p%orbital2%nContracted(SBB-1)
 ENDDO 
 nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
 DO iDerivP=1,NDerivP
  DO iPassP = 1, PQ%P%p%nPasses
   CALL getOverlapInfo(nContA,nContB,nAngA,nAngB,startA,startB,localA,localB,atomA,atomB,PQ%P%p,iAngmomP,iPassP)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
!   if(PQ%P%p%orbital1%type_nucleus)atomA=1
!   if(PQ%P%p%orbital2%type_nucleus)atomB=1
   batchA = PQ%P%p%orbital1%batch(iPassP)
   batchB = PQ%P%p%orbital2%batch(iPassP)
   iOrbP=iOrbitalP-1+PASSPOFFSET
   derivInfo%Atom(1)=atomA
   derivInfo%Atom(2)=atomB
   D_AB = Input%LST_DLHS%INDEX(ATOMA,ATOMB,1,1)
   D_BA = Input%LST_DLHS%INDEX(ATOMB,ATOMA,1,1)
   
   iOrbitalQ = 1
   endAngmomQ = PQ%Q%p%nAngmom
   IF (samePQ) endAngmomQ = iAngmomP
   DO iAngmomQ=1,endAngmomQ
    nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
    sC = 0
    DO SCC=2,PQ%Q%p%indexAng1(iAngmomQ)
       sC = sC+PQ%Q%p%orbital1%nOrbComp(SCC-1)*PQ%Q%p%orbital1%nContracted(SCC-1)
    ENDDO
    sD = 0
    DO SDD=2,PQ%Q%p%indexAng2(iAngmomQ)
       sD = sD+PQ%Q%p%orbital2%nOrbComp(SDD-1)*PQ%Q%p%orbital2%nContracted(SDD-1)
    ENDDO
    DO iderivQ =1,nderivQ
     DO iPassQ = 1, PQ%Q%p%nPasses
      CALL getOverlapInfo(nContC,nContD,nAngC,nAngD,startC,startD,localC,localD,atomC,atomD,PQ%Q%p,iAngmomQ,iPassQ)
      PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
      iOrbQ=iOrbitalQ-1+PASSQOFFSET
!      if(PQ%Q%p%orbital1%type_nucleus)atomC=1
!      if(PQ%Q%p%orbital2%type_nucleus)atomD=1
      batchC = PQ%Q%p%orbital1%batch(iPassQ)
      batchD = PQ%Q%p%orbital2%batch(iPassQ)
      derivInfo%atom(3) = atomC
      derivInfo%atom(4) = atomD
      D_AC = Input%LST_DLHS%INDEX(ATOMA,ATOMC,1,1)
      D_CA = Input%LST_DLHS%INDEX(ATOMC,ATOMA,1,1)

      ideriv=iderivP+(iderivQ-1)*NderivP
      iAtom = derivInfo%Atom(derivInfo%AO(1,iDeriv))
      iDer  = derivInfo%dirComp(iDeriv)
      translate = derivInfo%translate.GT.0
      iTrans = 1
      IF (translate) iTrans = derivInfo%Atom(derivInfo%translate)
      call distributeGrad2(&
           & RES%LSAO(iatom)%BATCH(1,1,1,1)%elms,&
           & RES%LSAO(itrans)%BATCH(1,1,1,1)%elms,&
           & Input%LST_DLHS%LSAO(D_AC)%BATCH(batchA,batchC,1,1)%elms,&
           & Input%LST_DLHS%LSAO(D_CA)%BATCH(batchC,batchA,1,1)%elms,&
           & Input%LST_DLHS%LSAO(D_AB)%BATCH(batchA,batchB,1,1)%elms,&
           & Input%LST_DLHS%LSAO(D_BA)%BATCH(batchB,batchA,1,1)%elms,&
           & iOrbP,iOrbQ,startA,startB,startC,startD,nAngA,nAngB,nAngC,nAngD,&
           & nContA,nContB,nContC,nContD,QPmat2,dimQ,dimP,nmat,dimA,&
           & dimB,dimC,dimD,sA,sB,sC,sD,lupri,iprint,iDer,iAtom,&
           & INPUT%sameLHSaos,INPUT%SameRHSaos,INPUT%SameODs,&
           & INPUT%AddToIntegral,derivOrder,derivinfo,ideriv,translate)
     ENDDO
     iOrbitalQ = iOrbitalQ + nOrbQ
    ENDDO
   ENDDO
  ENDDO
  iOrbitalP = iOrbitalP + nOrbP
 ENDDO
ENDDO
!write(lupri,*)'TK grad'
!call print_lstensor(RES,lupri,.TRUE.)
IF (derivOrder.GT.0) CALL freeDerivativeOverlapInfo(derivInfo)

END SUBROUTINE DISTRIBUTEGRAD

!> \brief distribute the gradient integrals to the lstensor format
!> \author T. Kjaergaard
!> \date 2009-02-05
!> \param grad the gradient 
!> \param grad2 the translational gradient 
!> \param D_AC the density \latexonly D_{ac} \endlatexonly 
!> \param D_CA the density \latexonly D_{ca} \endlatexonly 
!> \param D_AB the density \latexonly D_{ab} \endlatexonly 
!> \param D_BA the density \latexonly D_{ba} \endlatexonly 
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param startA the starting orbital index of the full A AO
!> \param startB the starting orbital index of the full B AO
!> \param startC the starting orbital index of the full C AO
!> \param startD the starting orbital index of the full D AO
!> \param nAngA number of angular momentums for the A batch
!> \param nAngB number of angular momentums for the B batch
!> \param nAngC number of angular momentums for the C batch
!> \param nAngD number of angular momentums for the D batch
!> \param nContA number of contracted functions on the A batch
!> \param nContB number of contracted functions on the B batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nContD number of contracted functions on the D batch 
!> \param QPmat2 matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat2
!> \param dimP dimension 2 of QPmat2
!> \param nmat the number of density matrices
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param ider the derivative index 
!> \param iatom the atomic index
!> \param sameLHSaos flag to indicate that the two left hand side AOs are the same
!> \param sameRHSaos flag to indicate that the two right hand side AOs are the same
!> \param sameODs flag to indicate that the two overlap distributions are the same
!> \param addtotintegral if the integral should be placed in (false) or added to the output grad (true)
!> \param derivorder the order of the derivative order requested
!> \param derinfo contains info about the derivation
!> \param ideriv the derivative index
!> \param translate flag to indicate that we use translational symmetry
subroutine distributeGrad2(grad,grad2,D_AC,D_CA,D_AB,D_BA,start_iOrbP,start_iOrbQ,&
     & startA,startB,startC,startD,nAngA,nAngB,nAngC,nAngD,nContA,nContB,&
     & nContC,nContD,QPmat2,dimQ,dimP,nmat,dimA,dimB,dimC,dimD,sA,sB,sC,sD,&
     & lupri,iprint,iDer,iAtom,sameLHSaos,SameRHSaos,SameODs,AddtoIntegral,derivOrder,derinfo,ideriv,translate)
implicit none
integer,intent(in)        :: nmat,lupri,iprint,dimQ,dimP,iDer,iAtom,dimA,dimB,dimC,dimD
real(realk),intent(inout) :: grad(3,nmat),grad2(3,nmat)
real(realk),intent(in)    :: D_AB(dimA,dimB,nmat),D_BA(dimB,dimA,nmat)
real(realk),intent(in)    :: D_AC(dimA,dimC,nmat),D_CA(dimC,dimA,nmat)
integer,intent(inout)     :: start_iOrbP,start_iOrbQ,ideriv
integer,intent(in)        :: startA,startB,startC,startD,nAngA,nAngB,nAngC,nAngD
integer,intent(in)        :: nContA,nContB,nContC,nContD,derivOrder,sA,sB,sC,sD
real(realk),intent(in)    :: QPmat2(dimQ,dimP)
Type(derivativeInfo),intent(in) :: derInfo
logical,intent(in) :: sameLHSaos,SameRHSaos,SameODs,AddtoIntegral,translate
!
integer                   :: iContB,iContA,iContC,iContD,iOrbP,iOrbQ
integer                   :: iAngB,iAngA,iAngC,iAngD,iB,iA,iC,iD,idmat,imat
real(realk)               :: integral,Dab,Dac,Dercont

do imat=1,nmat
 derCont = 0.d0
 iOrbP = START_iOrbP
 DO iContB=1,nContB
  DO iContA=1,nContA
   DO iAngB=1,nAngB
    iB=sB+iAngB+(iContB-1)*nAngB
    DO iAngA=1,nAngA
     iA=sA+iAngA+(iContA-1)*nAngA
     iOrbP=iOrbP+1            
     Dab = D_ab(iA,iB,imat)
     IF (sameLHSaos.AND.iA-sA+startA-1.NE.iB-sB+startB-1) Dab = Dab + D_ba(iB,iA,imat)
     iOrbQ = START_iOrbQ
     DO iContD=1,nContD
      DO iContC=1,nContC
       DO iAngD=1,nAngD
        iD=sD+iAngD+(iContD-1)*nAngD
        DO iAngC=1,nAngC
         iC=sC+iAngC+(iContC-1)*nAngC
         iOrbQ=iOrbQ+1
         integral  = QPmat2(iOrbQ,iOrbP)
         Dac = D_AC(iA,iC,imat)
         IF (sameODs.AND.iA-sA+startA-1.NE.iC-sC+startC-1) Dac = Dac + D_CA(iC,iA,imat)
         IF(AddtoIntegral)THEN
            IF (sameLHSaos.AND.sameRHSaos.AND.sameODs)THEN
               CALL LSQUIT('Implement!',lupri)
            ELSE IF (sameLHSaos.AND.sameRHSaos)THEN
               CALL LSQUIT('Implement!',lupri)
            ELSE IF (sameLHSaos)THEN
               CALL LSQUIT('Implement!',lupri)
            ELSE IF (sameRHSaos)THEN
               CALL LSQUIT('Implement!',lupri)
            ELSE IF (sameODs .AND. (sameLHSaos.NEQV.sameRHSaos))THEN
               CALL LSQUIT('Progamming error, sameODs but sameLHSaos.NE.sameRHSaos!',lupri)
            ELSE IF (sameODs)THEN
               CALL LSQUIT('Implement!',lupri)
            ELSE
               derCont = derCont + integral*Dab
            ENDIF
         ELSE
            IF (sameLHSaos.AND.sameRHSaos.AND.sameODs)THEN
               CALL LSQUIT('Implement!',lupri)
            ELSE IF (sameLHSaos.AND.sameRHSaos)THEN
               CALL LSQUIT('Implement!',lupri)
            ELSE IF (sameLHSaos)THEN
               derCont = derCont + integral*dab
            ELSE IF (sameRHSaos)THEN
               CALL LSQUIT('Implement!',lupri)
            ELSE IF (sameODs .AND. (sameLHSaos.NEQV.sameRHSaos))THEN
               CALL LSQUIT('Progamming error, sameODs but sameLHSaos.NE.sameRHSaos!',lupri)
            ELSE IF (sameODs)THEN
               derCont = derCont + integral*Dac
            ELSE
               derCont = derCont + integral*Dac
            ENDIF
         ENDIF
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
!$OMP CRITICAL(distributeGrd) 
 grad(iDer,imat) = grad(iDer,imat) + derCont
 IF (translate) THEN
    grad2(iDer,imat) = grad2(iDer,imat) - derCont
 ENDIF
!$OMP END CRITICAL(distributeGrd) 
ENDDO

end subroutine distributeGrad2

!> \brief new distributePQ to lstensor
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param RES contains the result lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE GeneraldistributePQ(contrule,RES,PQ,QPmat2,dimQ,dimP,Input,Lsoutput,LUPRI,IPRINT)
  implicit none !QQQQQQQQQQ
  Type(integrand),intent(in)      :: PQ
  Type(IntegralInput),intent(in)  :: Input
  Type(IntegralOutput),intent(inout) :: Lsoutput
  Integer,intent(in)              :: LUPRI,IPRINT
  Integer,intent(in)              :: dimQ,dimP
  REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
  TYPE(lstensor),intent(inout)    :: RES
  type(contractionrule)  :: contrule
  !
  Integer :: nAngmomP,nOrbP,nderivP,nMAT
  Integer :: iAngmomP,iOrbitalP,iderivP
  Integer :: iA,iB,sA,sB,sC,sD,SAA,SBB,SCC,SDD
  Integer :: nContA,nContB,nAngA,nAngB,startA,startB
  Integer :: iOrbP,nderivQ,iderivQ,ideriv
  Integer :: iOrbitalQ,endAngmomQ,iAngmomQ,nOrbQ
  Integer :: nContC,nContD,nAngC,nAngD,startC,startD
  Integer :: iC,iD,iOrbQ
!  Integer :: startA1,startB1,startC1,startD1
  real(realk) :: Factor
!  real(realk),pointer :: TMP5(:,:,:,:,:)
  integer :: PASSPOFFSET,PASSQOFFSET
  integer :: iPassP,iPassQ
  Integer :: dimABCD(5)!dimA,dimB,dimC,dimD,nderiv
!  Integer :: ABCD,BACD,ABDC,BADC,CDAB,CDBA,DCAB,DCBA
  Integer :: nderiv
!  Integer :: RHS_ABCD,RHS_BACD,RHS_ABDC,RHS_BADC
!  Integer :: RHS_CDAB,RHS_CDBA,RHS_DCAB,RHS_DCBA
  logical :: samePQ,SameRHSaos,SameLHSaos,SameODs,RHScontraction
  integer,target :: batch(5)
  integer,target :: atom(0:4),AtomA,atomB,atomC,atomD
  integer :: nAng(4),nCont(4)
  integer :: dimOABCD(5),dimOBACD(5),dimOABDC(5),dimOBADC(5)
  integer :: dimOCDAB(5),dimOCDBA(5),dimODCAB(5),dimODCBA(5)
  integer :: dimRABCD(5),dimRBACD(5),dimRABDC(5),dimRBADC(5)
  integer :: dimRCDAB(5),dimRCDBA(5),dimRDCAB(5),dimRDCBA(5)
  integer :: i1,startintegral(4),ODsym(5),ABsym(5),CDsym(5)
  !integer :: BADCsym(4),CDBAsym(4),DCABsym(4),DCABsym(4)
! type(integerp) :: AtomO(4),AtomR(4)
  integer :: batchOABCD(4),batchOABDC(4),batchOBACD(4),batchOBADC(4)
  integer :: batchOCDAB(4),batchODCAB(4),batchOCDBA(4),batchODCBA(4)
  !RHS
  integer :: batchRABCD(4),batchRABDC(4),batchRBACD(4),batchRBADC(4)
  integer :: batchRCDAB(4),batchRDCAB(4),batchRCDBA(4),batchRDCBA(4)
  !
  real(realk),pointer :: DUMMY(:)

  if(.NOT.lsoutput%docontrule)call LSquit('ended up in GeneraldistributePQ but lsoutput%docontrule is false',lupri)
  ODsym(1)=3
  ODsym(2)=4
  ODsym(3)=1
  ODsym(4)=2
  ODsym(5)=5

  ABsym(1)=2
  ABsym(2)=1
  ABsym(3)=3
  ABsym(4)=4
  ABsym(5)=5

  CDsym(1)=1
  CDsym(2)=2
  CDsym(3)=4
  CDsym(4)=3
  CDsym(5)=5

  nAngmomP = PQ%P%p%nAngmom
!  THRESHOLD = INPUT%CS_THRESHOLD
  SamePQ = PQ%samePQ
  SameLHSaos = INPUT%SameLHSaos
  SameRHSaos = INPUT%SameRHSaos
  SameODs = INPUT%SameODs
  nderivQ=input%nderivQ
  nderivP=Input%nDerivP
  nderiv = nderivP*nderivQ
  Factor =  contrule%factor
  RHScontraction = contrule%RHScontraction
  nmat=1
  if(RHScontraction)nMAT=Input%NDMAT_RHS
  iA = PQ%P%p%indexAng1(nAngmomP)
  dimABCD(1) = 0
  DO SAA=1,iA
     dimABCD(1) = dimABCD(1)+PQ%P%p%orbital1%nOrbComp(SAA)*PQ%P%p%orbital1%nContracted(SAA)
  ENDDO
  iB = PQ%P%p%indexAng2(nAngmomP)
  dimABCD(2) = 0
  DO SBB=1,iB
     dimABCD(2) = dimABCD(2)+PQ%P%p%orbital2%nOrbComp(SBB)*PQ%P%p%orbital2%nContracted(SBB)
  ENDDO
  iC = PQ%Q%p%indexAng1(PQ%Q%p%nAngmom)
  dimABCD(3) = 0
  DO SCC=1,iC
     dimABCD(3) = dimABCD(3)+PQ%Q%p%orbital1%nOrbComp(SCC)*PQ%Q%p%orbital1%nContracted(SCC)
  ENDDO
  iD = PQ%Q%p%indexAng2(PQ%Q%p%nAngmom)
  dimABCD(4) = 0
  DO SDD=1,iD
     dimABCD(4) = dimABCD(4)+PQ%Q%p%orbital2%nOrbComp(SDD)*PQ%Q%p%orbital2%nContracted(SDD)
  ENDDO
  IF(.NOT.RHScontraction.OR.(.NOT.(sameLHSaos.AND.sameRHSaos.AND.sameODs)))THEN
     NULLIFY(DUMMY)
     ALLOCATE(DUMMY(dimABCD(1)*dimABCD(2)*dimABCD(3)*dimABCD(4)*nderiv))
  ENDIF
  dimABCD(5) = nderiv
  atom(0)=1
  ! set up Output contraction loops
  do i1=1,5
     if(contrule%OutIntCont(i1).NE.0) then
        dimOABCD(i1) = dimABCD(contrule%OutIntCont(i1))
        dimOBACD(i1) = dimABCD(ABsym(contrule%OutIntCont(i1)))
        dimOABDC(i1) = dimABCD(CDsym(contrule%OutIntCont(i1)))
        dimOBADC(i1) = dimABCD(CDsym(ABsym(contrule%OutIntCont(i1))))

        dimOCDAB(i1) = dimABCD(ODsym(contrule%OutIntCont(i1)))
        dimOCDBA(i1) = dimABCD(ABsym(ODsym(contrule%OutIntCont(i1))))
        dimODCAB(i1) = dimABCD(CDsym(ODsym(contrule%OutIntCont(i1))))
        dimODCBA(i1) = dimABCD(CDsym(ABsym(ODsym(contrule%OutIntCont(i1)))))
     else
        dimOABCD(i1) = 1; dimOBACD(i1) = 1; dimOABDC(i1) = 1; dimOBADC(i1) = 1
        dimOCDAB(i1) = 1; dimOCDBA(i1) = 1; dimODCAB(i1) = 1; dimODCBA(i1) = 1
     endif
  enddo

  IF(RHScontraction)THEN
   ! set up RHS contraction loops
   do i1=1,5
     if(contrule%RhsIntCont(i1).NE.0)then
        dimRABCD(i1) = dimABCD(contrule%RhsIntCont(i1))
        dimRBACD(i1) = dimABCD(ABsym(contrule%RhsIntCont(i1)))
        dimRABDC(i1) = dimABCD(CDsym(contrule%RhsIntCont(i1)))
        dimRBADC(i1) = dimABCD(CDsym(ABsym(contrule%RhsIntCont(i1))))
        
        dimRCDAB(i1) = dimABCD(ODsym(contrule%RhsIntCont(i1)))
        dimRCDBA(i1) = dimABCD(ABsym(ODsym(contrule%RhsIntCont(i1))))
        dimRDCAB(i1) = dimABCD(CDsym(ODsym(contrule%RhsIntCont(i1))))
        dimRDCBA(i1) = dimABCD(ABsym(CDsym(ODsym(contrule%RhsIntCont(i1)))))
     else
        dimRABCD(i1) = 1; dimRBACD(i1) = 1; dimRABDC(i1) = 1; dimRBADC(i1) = 1
        dimRCDAB(i1) = 1; dimRCDBA(i1) = 1; dimRDCAB(i1) = 1; dimRDCBA(i1) = 1
     endif
   enddo
  ENDIF

  iOrbitalP = 1
  DO iDerivP=1,NDerivP
   DO iAngmomP=1,nAngmomP
!    nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
!    DO iPassP = 1, PQ%P%p%nPasses
!     CALL getOverlapInfo(nContA,nContB,nAngA,nAngB,startA,startB,localA,localB,atomA,atomB,PQ%P%p,iAngmomP,iPassP)
    nOrbP  = PQ%P%p%nOrbitals(iAngmomP)
    iA = PQ%P%p%indexAng1(iAngmomP)
    iB = PQ%P%p%indexAng2(iAngmomP)
    nContA = PQ%P%p%orbital1%nContracted(iA)
    nContB = PQ%P%p%orbital2%nContracted(iB)
    nAngA = PQ%P%p%orbital1%nOrbComp(iA)
    nAngB = PQ%P%p%orbital2%nOrbComp(iB)
    DO iPassP = 1, PQ%P%p%nPasses
     PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
     startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
     startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)
     atomA = PQ%P%p%orbital1%atom(iPassP)
     atomB = PQ%P%p%orbital2%atom(iPassP)
     if(PQ%P%p%orbital1%type_nucleus)atomA=1
     if(PQ%P%p%orbital2%type_nucleus)atomB=1
     atom(1) = AtomA
     atom(2) = AtomB
     batch(1) = PQ%P%p%orbital1%batch(iPassP)
     batch(2) = PQ%P%p%orbital2%batch(iPassP)
     iOrbP=iOrbitalP-1+PASSPOFFSET
     iOrbitalQ = 1
     endAngmomQ = PQ%Q%p%nAngmom
     IF (samePQ) endAngmomQ = iAngmomP
     DO iDerivQ=1,NDerivQ
      ideriv=iderivP+(iderivQ-1)*NderivP
      DO iAngmomQ=1,endAngmomQ
!       nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
!       DO iPassQ = 1, PQ%Q%p%nPasses
!        CALL getOverlapInfo(nContC,nContD,nAngC,nAngD,startC,startD,localC,localD,atomC,atomD,PQ%Q%p,iAngmomQ,iPassQ)
       nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
       iC = PQ%Q%p%indexAng1(iAngmomQ)
       iD = PQ%Q%p%indexAng2(iAngmomQ)
       nContC = PQ%Q%p%orbital1%nContracted(iC)
       nContD = PQ%Q%p%orbital2%nContracted(iD)
       nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
       nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
       DO iPassQ = 1, PQ%Q%p%nPasses
        PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
        iOrbQ=iOrbitalQ-1+PASSQOFFSET
        startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
        startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
        atomC = PQ%Q%p%orbital1%atom(iPassQ)
        atomD = PQ%Q%p%orbital2%atom(iPassQ)
        if(PQ%Q%p%orbital1%type_nucleus)atomC=1
        if(PQ%Q%p%orbital2%type_nucleus)atomD=1
        atom(3)=atomC
        atom(4)=atomD
        batch(3) = PQ%Q%p%orbital1%batch(iPassQ)
        batch(4) = PQ%Q%p%orbital2%batch(iPassQ)
        call determineBatch(contrule%OutIntCont,batchOABCD,batchOBACD,batchOABDC,batchOBADC,&
             & batchOCDAB,batchOCDBA,batchODCAB,batchODCBA,batch,ABsym,CDsym,ODsym)        
        IF(RHScontraction)THEN
           call determineBatch(contrule%RhsIntCont,batchRABCD,batchRBACD,batchRABDC,&
                & batchRBADC,batchRCDAB,batchRCDBA,batchRDCAB,batchRDCBA,batch,ABsym,CDsym,ODsym)
        ENDIF
        
        sA = 0
        DO SAA=2,iA
           sA = sA+PQ%P%p%orbital1%nOrbComp(SAA-1)*PQ%P%p%orbital1%nContracted(SAA-1)
        ENDDO
        
        sB = 0
        DO SBB=2,iB
           sB = sB+PQ%P%p%orbital2%nOrbComp(SBB-1)*PQ%P%p%orbital2%nContracted(SBB-1)
        ENDDO
        
        sC = 0
        DO SCC=2,iC
           sC = sC+PQ%Q%p%orbital1%nOrbComp(SCC-1)*PQ%Q%p%orbital1%nContracted(SCC-1)
        ENDDO
        
        sD = 0
        DO SDD=2,iD
           sD = sD+PQ%Q%p%orbital2%nOrbComp(SDD-1)*PQ%Q%p%orbital2%nContracted(SDD-1)
        ENDDO
        startintegral(1) = sA
        startintegral(2) = sB
        startintegral(3) = sC
        startintegral(4) = sD
        nAng(1) = nAngA
        nAng(2) = nAngB
        nAng(3) = nAngC
        nAng(4) = nAngD
        nCont(1) = nContA
        nCont(2) = nContB
        nCont(3) = nContC
        nCont(4) = nContD
        CALL GeneraldistributePQ2(RES,input%LST_DRHS,contrule,batchOABCD,batchOBACD,batchOABDC,batchOBADC,&
             & batchOCDAB,batchOCDBA,batchODCAB,batchODCBA,batchRABCD,batchRBACD,batchRABDC,batchRBADC,batchRCDAB,batchRCDBA,&
             & batchRDCAB,batchRDCBA,dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
             & dimRABCD,dimRBACD,dimRABDC,dimRBADC,dimRCDAB,dimRDCAB,dimRCDBA,dimRDCBA,RHScontraction,&
             & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
             & INPUT%AddToIntegral,SameLHSaos,sameRHSaos,sameODs,Factor,DUMMY,atom,atomA,atomB,atomC,atomD,lupri,iprint)
       ENDDO
       iOrbitalQ = iOrbitalQ + nOrbQ
      ENDDO
     ENDDO
    ENDDO
    iOrbitalP = iOrbitalP + nOrbP
   ENDDO
  ENDDO
  IF(.NOT.RHScontraction.OR.(.NOT.(sameLHSaos.AND.sameRHSaos.AND.sameODs)))THEN
     DEALLOCATE(DUMMY)
     NULLIFY(DUMMY)
  ENDIF
END SUBROUTINE GENERALDISTRIBUTEPQ

subroutine determineBatch(OutIntCont,batchOABCD,batchOBACD,batchOABDC,batchOBADC,&
     &batchOCDAB,batchOCDBA,batchODCAB,batchODCBA,batch,ABsym,CDsym,ODsym)
  implicit none 
  integer,intent(in) :: OutIntCont(5),batch(4),ABsym(5),CDsym(5),ODsym(5)
  integer,intent(inout) :: batchOABCD(4),batchOBACD(4),batchOABDC(4),batchOBADC(4)
  integer,intent(inout) :: batchOCDAB(4),batchOCDBA(4),batchODCAB(4),batchODCBA(4)
  integer :: i1
  
  do i1=1,4
     if(OutIntCont(i1).NE.0) then
        batchOABCD(i1) = batch(OutIntCont(i1))
        batchOBACD(i1) = batch(ABsym(OutIntCont(i1)))
        batchOABDC(i1) = batch(CDsym(OutIntCont(i1)))
        batchOBADC(i1) = batch(CDsym(ABsym(OutIntCont(i1))))
        
        batchOCDAB(i1) = batch(ODsym(OutIntCont(i1)))
        batchOCDBA(i1) = batch(ABsym(ODsym(OutIntCont(i1))))
        batchODCAB(i1) = batch(CDsym(ODsym(OutIntCont(i1))))
        batchODCBA(i1) = batch(ABsym(CDsym(ODsym(OutIntCont(i1)))))
     else
        batchOABCD(i1) = 1; batchOBACD(i1) = 1; batchOABDC(i1) = 1
        batchOBADC(i1) = 1; batchOCDAB(i1) = 1; batchOCDBA(i1) = 1
        batchODCAB(i1) = 1; batchODCBA(i1) = 1
     endif
  enddo
end subroutine determineBatch

subroutine GeneraldistributePQ2(RES,DRHS,contrule,batchOABCD,batchOBACD,batchOABDC,batchOBADC,&
             & batchOCDAB,batchOCDBA,batchODCAB,batchODCBA,batchRABCD,batchRBACD,batchRABDC,batchRBADC,batchRCDAB,batchRCDBA,&
             & batchRDCAB,batchRDCBA,dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
             & dimRABCD,dimRBACD,dimRABDC,dimRBADC,dimRCDAB,dimRDCAB,dimRCDBA,dimRDCBA,RHScontraction,&
             & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
             & AddToIntegral,SameLHSaos,sameRHSaos,sameODs,Factor,DUMMY,atom,atomA,atomB,atomC,atomD,lupri,iprint)
  implicit none
  type(lstensor),intent(inout) :: RES
  type(lstensor),intent(in)    :: DRHS
  type(contractionrule),intent(in)    :: contrule
  integer,intent(in) :: batchOABCD(4),batchOABDC(4),batchOBACD(4),batchOBADC(4)
  integer,intent(in) :: batchOCDAB(4),batchODCAB(4),batchOCDBA(4),batchODCBA(4)
  integer,intent(in) :: batchRABCD(4),batchRABDC(4),batchRBACD(4),batchRBADC(4)
  integer,intent(in) :: batchRCDAB(4),batchRDCAB(4),batchRCDBA(4),batchRDCBA(4)
  integer,intent(in) :: dimOABCD(5),dimOBACD(5),dimOABDC(5),dimOBADC(5)
  integer,intent(in) :: dimOCDAB(5),dimOCDBA(5),dimODCAB(5),dimODCBA(5)
  integer,intent(in) :: dimRABCD(5),dimRBACD(5),dimRABDC(5),dimRBADC(5)
  integer,intent(in) :: dimRCDAB(5),dimRCDBA(5),dimRDCAB(5),dimRDCBA(5),lupri,iprint
  integer,intent(in) :: atomA,atomB,atomC,atomD
  logical,intent(in) :: RHScontraction,AddToIntegral,SameLHSaos,sameRHSaos,sameODs
  integer,intent(inout) :: iOrbP,iOrbQ
  integer,intent(in) :: nderiv,startA,startB,startC,startD,nAng(4),nCont(4),nmat,ideriv,startintegral(4)
  integer,target :: atom(0:4)
  real(realk),intent(in) :: QPmat2(:,:),Factor
  real(realk),intent(inout) :: DUMMY(:)
  !
  type(integerp) :: atomO(4),atomR(4)
  Integer :: ABCD,BACD,ABDC,BADC,CDAB,CDBA,DCAB,DCBA
  Integer :: RHS_ABCD,RHS_BACD,RHS_ABDC,RHS_BADC
  Integer :: RHS_CDAB,RHS_CDBA,RHS_DCAB,RHS_DCBA

  atom(0)=1
  atom(1)=atomA
  atom(2)=atomB
  atom(3)=atomC
  atom(4)=atomD

  atomO(1)%p => atom(contrule%OutIntCont(1))
  atomO(2)%p => atom(contrule%OutIntCont(2))
  atomO(3)%p => atom(contrule%OutIntCont(3))
  atomO(4)%p => atom(contrule%OutIntCont(4))
  IF(RHScontraction)then
     atomR(1)%p => atom(contrule%RhsIntCont(1))
     atomR(2)%p => atom(contrule%RhsIntCont(2))
     atomR(3)%p => atom(contrule%RhsIntCont(3))
     atomR(4)%p => atom(contrule%RhsIntCont(4))
  ENDIF
  
  ABCD = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
  IF(RHScontraction)then
   RHS_ABCD = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
   IF (sameLHSaos.AND.sameRHSaos.AND.sameODs)THEN
    atom(1) = atomB; atom(2)=atomA;
    BACD = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_BACD = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomA; atom(2)=atomB
    atom(3) = atomD; atom(4)=atomC;
    ABDC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_ABDC = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomB; atom(2)=atomA;
    BADC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_BADC = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomC; atom(2)=atomD; atom(3)=atomA; atom(4)=atomB
    CDAB = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_CDAB = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomD; atom(2)=atomC;
    DCAB = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_DCAB = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(3)=atomB; atom(4)=atomA
    DCBA = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_DCBA = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomC; atom(2)=atomD;
    CDBA = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_CDBA = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomA; atom(2)=atomB;atom(3) = atomC; atom(4)=atomD; !back to normal
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & RES%LSAO(BACD)%BATCH(batchOBACD(1),batchOBACD(2),batchOBACD(3),batchOBACD(4))%elms,& 
         & RES%LSAO(ABDC)%BATCH(batchOABDC(1),batchOABDC(2),batchOABDC(3),batchOABDC(4))%elms,& 
         & RES%LSAO(BADC)%BATCH(batchOBADC(1),batchOBADC(2),batchOBADC(3),batchOBADC(4))%elms,& 
         & RES%LSAO(CDAB)%BATCH(batchOCDAB(1),batchOCDAB(2),batchOCDAB(3),batchOCDAB(4))%elms,& 
         & RES%LSAO(DCAB)%BATCH(batchODCAB(1),batchODCAB(2),batchODCAB(3),batchODCAB(4))%elms,& 
         & RES%LSAO(CDBA)%BATCH(batchOCDBA(1),batchOCDBA(2),batchOCDBA(3),batchOCDBA(4))%elms,& 
         & RES%LSAO(DCBA)%BATCH(batchODCBA(1),batchODCBA(2),batchODCBA(3),batchODCBA(4))%elms,& 
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DRHS%LSAO(RHS_ABCD)%BATCH(batchRABCD(1),batchRABCD(2),batchRABCD(3),batchRABCD(4))%elms,&
         & DRHS%LSAO(RHS_BACD)%BATCH(batchRBACD(1),batchRBACD(2),batchRBACD(3),batchRBACD(4))%elms,&
         & DRHS%LSAO(RHS_ABDC)%BATCH(batchRABDC(1),batchRABDC(2),batchRABDC(3),batchRABDC(4))%elms,&
         & DRHS%LSAO(RHS_BADC)%BATCH(batchRBADC(1),batchRBADC(2),batchRBADC(3),batchRBADC(4))%elms,&
         & DRHS%LSAO(RHS_CDAB)%BATCH(batchRCDAB(1),batchRCDAB(2),batchRCDAB(3),batchRCDAB(4))%elms,&
         & DRHS%LSAO(RHS_DCAB)%BATCH(batchRDCAB(1),batchRDCAB(2),batchRDCAB(3),batchRDCAB(4))%elms,&
         & DRHS%LSAO(RHS_CDBA)%BATCH(batchRCDBA(1),batchRCDBA(2),batchRCDBA(3),batchRCDBA(4))%elms,&
         & DRHS%LSAO(RHS_DCBA)%BATCH(batchRDCBA(1),batchRDCBA(2),batchRDCBA(3),batchRDCBA(4))%elms,&
         & dimRABCD,dimRBACD,dimRABDC,dimRBADC,dimRCDAB,dimRDCAB,dimRCDBA,dimRDCBA,&
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE IF (sameLHSaos.AND.sameRHSaos)THEN
    atom(1) = atomB; atom(2)=atomA;! atom(3)=atomC; atom(4)=atomD
    BACD = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_BACD = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomA; atom(2)=atomB; atom(3) = atomD; atom(4)=atomC;
    ABDC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_ABDC = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomB; atom(2)=atomA;
    BADC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_BADC = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomA; atom(2)=atomB; atom(3)=atomC; atom(4)=atomD !back to normal  
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & RES%LSAO(BACD)%BATCH(batchOBACD(1),batchOBACD(2),batchOBACD(3),batchOBACD(4))%elms,& 
         & RES%LSAO(ABDC)%BATCH(batchOABDC(1),batchOABDC(2),batchOABDC(3),batchOABDC(4))%elms,& 
         & RES%LSAO(BADC)%BATCH(batchOBADC(1),batchOBADC(2),batchOBADC(3),batchOBADC(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DRHS%LSAO(RHS_ABCD)%BATCH(batchRABCD(1),batchRABCD(2),batchRABCD(3),batchRABCD(4))%elms,&
         & DRHS%LSAO(RHS_BACD)%BATCH(batchRBACD(1),batchRBACD(2),batchRBACD(3),batchRBACD(4))%elms,&
         & DRHS%LSAO(RHS_ABDC)%BATCH(batchRABDC(1),batchRABDC(2),batchRABDC(3),batchRABDC(4))%elms,&
         & DRHS%LSAO(RHS_BADC)%BATCH(batchRBADC(1),batchRBADC(2),batchRBADC(3),batchRBADC(4))%elms,&
         & DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimRABCD,dimRBACD,dimRABDC,dimRBADC,dimRCDAB,dimRDCAB,dimRCDBA,dimRDCBA,&
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE IF (sameLHSaos)THEN   
    atom(1) = atomB; atom(2)=atomA;! atom(3)=atomC; atom(4)=atomD
    BACD = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_BACD = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomA; atom(2)=atomB; !back to normal
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & RES%LSAO(BACD)%BATCH(batchOBACD(1),batchOBACD(2),batchOBACD(3),batchOBACD(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DRHS%LSAO(RHS_ABCD)%BATCH(batchRABCD(1),batchRABCD(2),batchRABCD(3),batchRABCD(4))%elms,&
         & DRHS%LSAO(RHS_BACD)%BATCH(batchRBACD(1),batchRBACD(2),batchRBACD(3),batchRBACD(4))%elms,&
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimRABCD,dimRBACD,dimRABDC,dimRBADC,dimRCDAB,dimRDCAB,dimRCDBA,dimRDCBA,&
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE IF (sameRHSaos)THEN   
    atom(3) = atomD; atom(4)=atomC;
    ABDC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_ABDC = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(3)=atomC; atom(4)=atomD  !back to normal
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & DUMMY,&
         & RES%LSAO(ABDC)%BATCH(batchOABDC(1),batchOABDC(2),batchOABDC(3),batchOABDC(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DRHS%LSAO(RHS_ABCD)%BATCH(batchRABCD(1),batchRABCD(2),batchRABCD(3),batchRABCD(4))%elms,&
         & DUMMY,&
         & DRHS%LSAO(RHS_ABDC)%BATCH(batchRABDC(1),batchRABDC(2),batchRABDC(3),batchRABDC(4))%elms,&
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimRABCD,dimRBACD,dimRABDC,dimRBADC,dimRCDAB,dimRDCAB,dimRCDBA,dimRDCBA,&
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE IF (sameODs)THEN
    atom(1) = atomC; atom(2)=atomD; atom(3)=atomA; atom(4)=atomB
    CDAB = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    RHS_CDAB = RES%INDEX(ATOMR(1)%p,ATOMR(2)%p,ATOMR(3)%p,ATOMR(4)%p)
    atom(1) = atomA; atom(2)=atomB; atom(3)=atomC; atom(4)=atomD !back to normal
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,&
         & RES%LSAO(CDAB)%BATCH(batchOCDAB(1),batchOCDAB(2),batchOCDAB(3),batchOCDAB(4))%elms,&
         & DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DRHS%LSAO(RHS_ABCD)%BATCH(batchRABCD(1),batchRABCD(2),batchRABCD(3),batchRABCD(4))%elms,&
         & DUMMY,DUMMY,DUMMY,&
         & DRHS%LSAO(RHS_CDAB)%BATCH(batchRCDAB(1),batchRCDAB(2),batchRCDAB(3),batchRCDAB(4))%elms,&
         & DUMMY,DUMMY,DUMMY,&
         & dimRABCD,dimRBACD,dimRABDC,dimRBADC,dimRCDAB,dimRDCAB,dimRCDBA,dimRDCBA,&
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE 
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DRHS%LSAO(RHS_ABCD)%BATCH(batchRABCD(1),batchRABCD(2),batchRABCD(3),batchRABCD(4))%elms,&
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimRABCD,dimRBACD,dimRABDC,dimRBADC,dimRCDAB,dimRDCAB,dimRCDBA,dimRDCBA,&
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ENDIF
  ELSE! no RHScontraction
   IF (sameLHSaos.AND.sameRHSaos.AND.sameODs)THEN
    atom(1) = atomB; atom(2)=atomA;
    BACD = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomA; atom(2)=atomB
    atom(3) = atomD; atom(4)=atomC;
    ABDC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomB; atom(2)=atomA;
    BADC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomC; atom(2)=atomD
    atom(3)=atomA; atom(4)=atomB
    CDAB = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomD; atom(2)=atomC;
    DCAB = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(3)=atomB; atom(4)=atomA
    DCBA = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomC; atom(2)=atomD;
    CDBA = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomA; atom(2)=atomB;atom(3) = atomC; atom(4)=atomD; !back to normal
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & RES%LSAO(BACD)%BATCH(batchOBACD(1),batchOBACD(2),batchOBACD(3),batchOBACD(4))%elms,& 
         & RES%LSAO(ABDC)%BATCH(batchOABDC(1),batchOABDC(2),batchOABDC(3),batchOABDC(4))%elms,& 
         & RES%LSAO(BADC)%BATCH(batchOBADC(1),batchOBADC(2),batchOBADC(3),batchOBADC(4))%elms,& 
         & RES%LSAO(CDAB)%BATCH(batchOCDAB(1),batchOCDAB(2),batchOCDAB(3),batchOCDAB(4))%elms,& 
         & RES%LSAO(DCAB)%BATCH(batchODCAB(1),batchODCAB(2),batchODCAB(3),batchODCAB(4))%elms,& 
         & RES%LSAO(CDBA)%BATCH(batchOCDBA(1),batchOCDBA(2),batchOCDBA(3),batchOCDBA(4))%elms,& 
         & RES%LSAO(DCBA)%BATCH(batchODCBA(1),batchODCBA(2),batchODCBA(3),batchODCBA(4))%elms,& 
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&      
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,& !dummy input
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE IF (sameLHSaos.AND.sameRHSaos)THEN
    atom(1) = atomB; atom(2)=atomA;! atom(3)=atomC; atom(4)=atomD
    BACD = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomA; atom(2)=atomB; atom(3) = atomD; atom(4)=atomC;
    ABDC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomB; atom(2)=atomA;
    BADC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomA; atom(2)=atomB; atom(3)=atomC; atom(4)=atomD !back to normal  
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & RES%LSAO(BACD)%BATCH(batchOBACD(1),batchOBACD(2),batchOBACD(3),batchOBACD(4))%elms,& 
         & RES%LSAO(ABDC)%BATCH(batchOABDC(1),batchOABDC(2),batchOABDC(3),batchOABDC(4))%elms,& 
         & RES%LSAO(BADC)%BATCH(batchOBADC(1),batchOBADC(2),batchOBADC(3),batchOBADC(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,& !dummy input
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE IF (sameLHSaos)THEN   
    atom(1) = atomB; atom(2)=atomA;! atom(3)=atomC; atom(4)=atomD
    BACD = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomA; atom(2)=atomB; !back to normal
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & RES%LSAO(BACD)%BATCH(batchOBACD(1),batchOBACD(2),batchOBACD(3),batchOBACD(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,& !dummy input
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE IF (sameRHSaos)THEN   
    atom(3) = atomD; atom(4)=atomC;
    ABDC = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(3)=atomC; atom(4)=atomD  !back to normal
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & DUMMY,&
         & RES%LSAO(ABDC)%BATCH(batchOABDC(1),batchOABDC(2),batchOABDC(3),batchOABDC(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,& !dummy input
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE IF (sameODs)THEN
    atom(1) = atomC; atom(2)=atomD; atom(3)=atomA; atom(4)=atomB
    CDAB = RES%INDEX(ATOMO(1)%p,ATOMO(2)%p,ATOMO(3)%p,ATOMO(4)%p)
    atom(1) = atomA; atom(2)=atomB; atom(3)=atomC; atom(4)=atomD !back to normal
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,&
         & RES%LSAO(CDAB)%BATCH(batchOCDAB(1),batchOCDAB(2),batchOCDAB(3),batchOCDAB(4))%elms,&
         & DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,& !dummy input
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ELSE
    CALL generaldistributePQ3(& 
         & RES%LSAO(ABCD)%BATCH(batchOABCD(1),batchOABCD(2),batchOABCD(3),batchOABCD(4))%elms,& 
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
         & DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,&
         & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,& !dummy input
         & contrule%OutIntCont,contrule%RhsIntCont,RHScontraction,&
         & iOrbP,iOrbQ,nderiv,startA,startB,startC,startD,&
         & nAng,nCont,startintegral,QPmat2,nmat,ideriv,&
         & Addtointegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
   ENDIF
  ENDIF
END subroutine GENERALDISTRIBUTEPQ2

!!$!> \brief ge
!!$!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!!$!> \date 2009-02-05
!!$!> \param OUT_ABCD the (A,B,C,D,ideriv) tensor block
!!$!> \param OUT_BACD the (B,A,C,D,ideriv) tensor block
!!$!> \param OUT_ABDC the (A,B,D,C,ideriv) tensor block
!!$!> \param OUT_BADC the (B,A,D,C,ideriv) tensor block
!!$!> \param OUT_CDAB the (C,D,A,B,ideriv) tensor block
!!$!> \param OUT_DCAB the (D,C,A,B,ideriv) tensor block
!!$!> \param OUT_CDBA the (C,D,B,A,ideriv) tensor block
!!$!> \param OUT_DCBA the (D,C,B,A,ideriv) tensor block
!!$!> \param START_iOrbP start index for orbital P (index 2 in QPmat)
!!$!> \param START_iOrbQ start index for orbital Q (index 1 in QPmat)
!!$!> \param dimA dim of lstensors A dimension 
!!$!> \param dimB dim of lstensors B dimension 
!!$!> \param dimC dim of lstensors C dimension 
!!$!> \param dimD dim of lstensors D dimension 
!!$!> \param nderiv dim of lstensors 5. dimension 
!!$!> \param nAngA number of orbital components for atomic orbital on center A
!!$!> \param nContA number of contracted orbitals for atomic orbital on center A
!!$!> \param nAngB number of orbital components for atomic orbital on center B
!!$!> \param nContB number of contracted orbitals for atomic orbital on center B
!!$!> \param nAngC number of orbital components for atomic orbital on center C
!!$!> \param nContC number of contracted orbitals for atomic orbital on center C
!!$!> \param nAngD number of orbital components for atomic orbital on center D
!!$!> \param nContD number of contracted orbitals for atomic orbital on center D
!!$!> \param startA start Orbital index for A index 
!!$!> \param startB start Orbital index for B index 
!!$!> \param startC start Orbital index for C index 
!!$!> \param startD start Orbital index for D index 
!!$!> \param sA start lstensorelm index for A index 
!!$!> \param sB start lstensorelm index for B index 
!!$!> \param sC start lstensorelm index for C index 
!!$!> \param sD start lstensorelm index for D index 
!!$!> \param QPmat4 matrix containing calculated integrals
!!$!> \param nmat the number of density matrices 
SUBROUTINE generaldistributePQ3(OUT_ABCD,OUT_BACD,OUT_ABDC,OUT_BADC,&
     & OUT_CDAB,OUT_DCAB,OUT_CDBA,OUT_DCBA,&
     & dimOABCD,dimOBACD,dimOABDC,dimOBADC,dimOCDAB,dimODCAB,dimOCDBA,dimODCBA,&
     & RHS_ABCD,RHS_BACD,RHS_ABDC,RHS_BADC,&
     & RHS_CDAB,RHS_DCAB,RHS_CDBA,RHS_DCBA,&
     & dimRABCD,dimRBACD,dimRABDC,dimRBADC,dimRCDAB,dimRDCAB,dimRCDBA,dimRDCBA,&
     & OutIntCont,RhsIntCont,RHScontraction,&
     & START_iOrbP,START_iOrbQ,nderiv,startA,startB,startC,startD,&
     & nAng,nCont,startintegral,QPmat4,nmat,ideriv,&
     & AddToIntegral,SameLHSaos,sameRHSaos,sameODs,Factor,lupri,iprint)
IMPLICIT NONE
REAL(REALK),intent(in)    :: QPmat4(:,:) !
!  INTEGER,intent(in)        :: dimA,dimB,dimC,dimD
INTEGER,intent(in)        :: lupri,iprint
LOGICAL,intent(in)        :: RHScontraction
INTEGER,intent(in)        :: dimOABCD(5),dimOBACD(5),dimOABDC(5),dimOBADC(5)
INTEGER,intent(in)        :: dimOCDAB(5),dimODCAB(5),dimOCDBA(5),dimODCBA(5)
INTEGER,intent(in)        :: dimRABCD(5),dimRBACD(5),dimRABDC(5),dimRBADC(5)
INTEGER,intent(in)        :: dimRCDAB(5),dimRDCAB(5),dimRCDBA(5),dimRDCBA(5)
REAL(REALK),intent(inout) :: OUT_ABCD(dimOABCD(1),dimOABCD(2),dimOABCD(3),dimOABCD(4),dimOABCD(5)*nmat)
REAL(REALK),intent(inout) :: OUT_BACD(dimOBACD(1),dimOBACD(2),dimOBACD(3),dimOBACD(4),dimOBACD(5)*nmat)
REAL(REALK),intent(inout) :: OUT_ABDC(dimOABDC(1),dimOABDC(2),dimOABDC(3),dimOABDC(4),dimOABDC(5)*nmat)
REAL(REALK),intent(inout) :: OUT_BADC(dimOBADC(1),dimOBADC(2),dimOBADC(3),dimOBADC(4),dimOBADC(5)*nmat)
REAL(REALK),intent(inout) :: OUT_CDAB(dimOCDAB(1),dimOCDAB(2),dimOCDAB(3),dimOCDAB(4),dimOCDAB(5)*nmat)
REAL(REALK),intent(inout) :: OUT_DCAB(dimODCAB(1),dimODCAB(2),dimODCAB(3),dimODCAB(4),dimODCAB(5)*nmat)
REAL(REALK),intent(inout) :: OUT_CDBA(dimOCDBA(1),dimOCDBA(2),dimOCDBA(3),dimOCDBA(4),dimOCDBA(5)*nmat)
REAL(REALK),intent(inout) :: OUT_DCBA(dimODCBA(1),dimODCBA(2),dimODCBA(3),dimODCBA(4),dimODCBA(5)*nmat)
!
REAL(REALK),intent(in) :: RHS_ABCD(dimRABCD(1),dimRABCD(2),dimRABCD(3),dimRABCD(4),dimRABCD(5)*nmat)
REAL(REALK),intent(in) :: RHS_BACD(dimRBACD(1),dimRBACD(2),dimRBACD(3),dimRBACD(4),dimRBACD(5)*nmat) 
REAL(REALK),intent(in) :: RHS_ABDC(dimRABDC(1),dimRABDC(2),dimRABDC(3),dimRABDC(4),dimRABDC(5)*nmat) 
REAL(REALK),intent(in) :: RHS_BADC(dimRBADC(1),dimRBADC(2),dimRBADC(3),dimRBADC(4),dimRBADC(5)*nmat) 
REAL(REALK),intent(in) :: RHS_CDAB(dimRCDAB(1),dimRCDAB(2),dimRCDAB(3),dimRCDAB(4),dimRCDAB(5)*nmat) 
REAL(REALK),intent(in) :: RHS_DCAB(dimRDCAB(1),dimRDCAB(2),dimRDCAB(3),dimRDCAB(4),dimRDCAB(5)*nmat) 
REAL(REALK),intent(in) :: RHS_CDBA(dimRCDBA(1),dimRCDBA(2),dimRCDBA(3),dimRCDBA(4),dimRCDBA(5)*nmat) 
REAL(REALK),intent(in) :: RHS_DCBA(dimRDCBA(1),dimRDCBA(2),dimRDCBA(3),dimRDCBA(4),dimRDCBA(5)*nmat) 
REAL(REALK),intent(in) :: Factor
!
INTEGER,intent(in)        :: startA,startB,startC,startD
INTEGER,intent(in)        :: OutIntCont(5),RhsIntCont(5)
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ
INTEGER,intent(in)        :: startintegral(4)
INTEGER,intent(in)        :: nCont(4),nAng(4)
INTEGER,intent(in)        :: nmat,ideriv,nderiv
LOGICAL                   :: AddToIntegral,SameLHSaos,sameRHSaos,sameODs
!
REAL(REALK)    :: TMPintegral(nAng(1)*nCont(1),nAng(2)*nCont(2),nAng(3)*nCont(3),nAng(4)*nCont(4))
REAL(REALK)    :: Int
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iA,iB,iC,iD,iOrbQ
INTEGER        :: nA,nB,nC,nD

INTEGER  :: loopO1,loopO2,loopO3,loopO4,loopO5
INTEGER  :: loopR1,loopR2,loopR3,loopR4,loopR5
INTEGER,pointer  :: A,B,C,D
!  INTEGER,pointer  :: O1,O2,O3,O4
!  INTEGER,pointer  :: R1,R2,R3,R4
type(INTEGERP)   :: O(5),R(5)!,OCDAB(5),RCDAB(5)
INTEGER,target :: OUT(5),INTEGRAL(5),RHSI(5),startI(5)
type(INTEGERP)  :: sO(5),sR(5)
!  INTEGER  :: sOCDAB(5),sRCDAB(5)
!  INTEGER  :: sO1,sO2,sO3,sO4
!  INTEGER  :: sR1,sR2,sR3,sR4
INTEGER  :: i1,i2,i3,i4,i9,i10,i11,i12,i13,i14,imat

loopO1 = 1; loopO2 = 1; loopO3 = 1; loopO4 = 1; loopO5 = 1
loopR1 = 1; loopR2 = 1; loopR3 = 1; loopR4 = 1; loopR5 = 1
startI(1) = startintegral(1); startI(2) = startintegral(2); 
startI(3) = startintegral(3); startI(4) = startintegral(4); 
startI(5) = 1
nA = nAng(1)*nCont(1)
nB = nAng(2)*nCont(2)
nC = nAng(3)*nCont(3)
nD = nAng(4)*nCont(4)
A => INTEGRAL(1)
B => INTEGRAL(2)
C => INTEGRAL(3)
D => INTEGRAL(4)
INTEGRAL(5) = ideriv

iOrbP = START_iOrbP
DO iContB=1,nCont(2)
 DO iContA=1,nCont(1)
  DO iAngB=1,nAng(2)
   iB=iAngB+(iContB-1)*nAng(2)
   DO iAngA=1,nAng(1)
    iA=iAngA+(iContA-1)*nAng(1)
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nCont(4)
     DO iContC=1,nCont(3)
      DO iAngD=1,nAng(4)
       iD=iAngD+(iContD-1)*nAng(4)
       DO iAngC=1,nAng(3)
        iC=iAngC+(iContC-1)*nAng(3)
        iOrbQ=iOrbQ+1
        TMPintegral(iA,iB,iC,iD) = Factor*QPmat4(iOrbQ,iOrbP)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

IF(RHScontraction)THEN
 ! set up Output contraction loops
 do i1=1,4
    if(OutIntCont(i1).NE.0) then
       O(i1)%p => INTEGRAL(OutIntCont(i1)); sO(i1)%p => startI(OutIntCont(i1))
    else
       nullify(O(i1)%p); allocate(O(i1)%p); O(i1)%p = 1; 
       nullify(sO(i1)%p); allocate(sO(i1)%p); sO(i1)%p = 0
    endif
 enddo
 nullify(O(5)%p); allocate(O(5)%p); O(5)%p = 1

 ! set up RHS contraction loops
 do i1=1,4
    if(RhsIntCont(i1).NE.0)then
       R(i1)%p => INTEGRAL(RhsIntCont(i1)); sR(i1)%p => startI(RhsIntCont(i1)) 
    else
       nullify(R(i1)%p); allocate(R(i1)%p); R(i1)%p = 1
       nullify(sR(i1)%p); allocate(sR(i1)%p); sR(i1)%p = 0
    endif
 enddo

 nullify(R(5)%p); allocate(R(5)%p); R(5)%p = 1

 !OUT(O1,O2,O3,O4) = OUT(O1,O2,O3,O4)+IN(A,B,C,D)*RHS(R1,R2,R3,R4)
 ! ASSUME NO PERMUTATIONAL SYMMETRY
 !  write(lupri,*)'nA-D',nA,nB,nC,nD
 DO i1=1,loopO1
  OUT(1) = i1
  DO i2=1,loopO2
   OUT(2) = i2
   DO i3=1,loopO3
    OUT(3) = i3
    DO i4=1,loopO4
     OUT(4)=i4
     DO i13=1,loopO5
      OUT(5)=i13
      !       write(lupri,*)'inside O loops'
      DO i9=1,loopR1
       RHSI(1)=i9
       DO i10=1,loopR2
        RHSI(2)=10
        DO i11=1,loopR3
         RHSI(3)=i11
         DO i12=1,loopR4
          RHSI(4)=i12
          DO i14=1,loopR5
           RHSI(5)=i14
           !            write(lupri,*)'inside O and R loops'
           DO iA=1,nA
            INTEGRAL(1) = iA
            DO iB=1,nB
             INTEGRAL(2)= iB
             DO iC=1,nC
              INTEGRAL(3)= iC
              DO iD=1,nD
               !                  write(lupri,*)'inside ABCD loops'
               INTEGRAL(4)= iD
               int = TMPintegral(A,B,C,D)
               !                  write(lupri,*)'startA-1+R(1)%p+sR(1)%p',startA-1+iA+startintegral(1)
               !                  write(lupri,*)'startB-1+R(2)%p+sR(2)%p',startB-1+iB+startintegral(2)
               !                  IF((startA-1+iA+startintegral(1) .EQ.1).AND.(startB-1+iB+startintegral(2) .EQ.1))&
               !&                  write(lupri,*)'TMPintegral2(',A,B,C,D,')=',int
               do imat=0,(nmat-1)*nderiv,nderiv
               IF(sameLHSaos.AND.(startA.NE.startB))THEN
                IF (sameODs.AND.((startA.NE.startC).OR.(startB.NE.startD))) THEN
                 IF (sameRHSaos.AND.(startC.NE.startD))THEN
                    OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_ABCD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BACD
                    startI(1) = startintegral(2); startI(2) = startintegral(1);                            
                    OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_BACD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                    startI(1) = startintegral(1); startI(2) = startintegral(2);
                    startI(3) = startintegral(4); startI(4) = startintegral(3);
                    OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_ABDC(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BADC
                    startI(1) = startintegral(2); startI(2) = startintegral(1); 
                    OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_BADC(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                    startI(1) = startintegral(3); startI(2) = startintegral(4); 
                    startI(3) = startintegral(1); startI(4) = startintegral(2); 
                    OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_CDAB(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iD; INTEGRAL(2) = iC;                                     !DCAB
                    startI(1) = startintegral(4); startI(2) = startintegral(3); 
                    OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_DCAB(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iB; INTEGRAL(4) = iA  !CDBA
                    startI(1) = startintegral(3); startI(2) = startintegral(4); 
                    startI(3) = startintegral(2); startI(4) = startintegral(1); 
                    OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_CDBA(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iD; INTEGRAL(2) = iC;                                     !DCBA
                    startI(1) = startintegral(4); startI(2) = startintegral(3); 
                    OUT_DCBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_DCBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_DCBA(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD; !Back to default ABCD
                    startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                    startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                 ELSE  !startC=startD no CD permutation
                    OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_ABCD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat) 
                    INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BACD
                    startI(1) = startintegral(2); startI(2) = startintegral(1);                            
                    OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p +imat)+&
                         &int*RHS_BACD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                    startI(1) = startintegral(3); startI(2) = startintegral(4);                            
                    startI(3) = startintegral(1); startI(4) = startintegral(2); 
                    OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_CDAB(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(3) = iB; INTEGRAL(4) = iA                                      !CDBA
                    startI(3) = startintegral(2); startI(4) = startintegral(1); 
                    OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_CDBA(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                    startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                    startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                 ENDIF
                ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN !no OD permutation
                 IF (sameRHSaos.AND.(startC.NE.startD))THEN
                    OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_ABCD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BACD
                    startI(1) = startintegral(2); startI(2) = startintegral(1); 
                    OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_BACD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                    startI(1) = startintegral(1); startI(2) = startintegral(2); 
                    startI(3) = startintegral(4); startI(4) = startintegral(3); 
                    OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_ABDC(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BADC
                    startI(1) = startintegral(2); startI(2) = startintegral(1); 
                    OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_BADC(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                    startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                    startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                 ELSE  !startC = startD   no OD or CD permutation
                    OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                         &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+&
                         &int*RHS_ABCD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iB; INTEGRAL(2) = iA; !INTEGRAL(3) = iC; INTEGRAL(4) = iD  !BACD
                    startI(1) = startintegral(2); startI(2) = startintegral(1); 
                    !                           startI(3) = startintegral(1); startI(4) = startintegral(2); 
                    OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                         &OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+&
                         &int*RHS_BACD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iA; INTEGRAL(2) = iB; !Back to default ABCD
                    startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                 ENDIF
                ENDIF
               ELSE !NOT (startA NE startB)     !no AB symmetry
                IF(sameODs.AND.((startA.NE.startC).OR.(startB.NE.startD))) THEN
                 IF (sameRHSaos.AND.(startC.NE.startD))THEN
                    OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_ABCD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                    startI(3) = startintegral(4); startI(4) = startintegral(3); 
                    OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_ABDC(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                    startI(1) = startintegral(3); startI(2) = startintegral(4); 
                    startI(3) = startintegral(1); startI(4) = startintegral(2); 
                    OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+&
                         &int*RHS_CDAB(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iD; INTEGRAL(2) = iC; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !DCAB
                    startI(1) = startintegral(4); startI(2) = startintegral(3); 
                    OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                         &OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+&
                         &int*RHS_DCAB(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                    startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                    startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                 ELSE !startC.NE.startD     !no CD symmetry 
                    OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                         &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+&
                         &int*RHS_ABCD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                    startI(1) = startintegral(3); startI(2) = startintegral(4); 
                    startI(3) = startintegral(1); startI(4) = startintegral(2); 
                    OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                         &OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+&
                         &int*RHS_CDAB(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                    startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                    startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                 ENDIF
                ELSE   !not ((startA.NE.startC).OR.(startB.NE.startD)) ! no od
                 IF (sameRHSaos.AND.(startC.NE.startD))THEN
                    OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                         &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+&
                         &int*RHS_ABCD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                    startI(3) = startintegral(4); startI(4) = startintegral(3); 
                    OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                         &OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+&
                         &int*RHS_ABDC(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                    INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                    startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                 ELSE !no permutational symmetry         ! only ABCD
                    OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                         &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+&
                         &int*RHS_ABCD(R(1)%p+sR(1)%p,R(2)%p+sR(2)%p,R(3)%p+sR(3)%p,R(4)%p+sR(4)%p,R(5)%p+imat)  
                 ENDIF
                ENDIF
               ENDIF
               enddo
              ENDDO !iD
             ENDDO !iC
            ENDDO !iB
           ENDDO !iA
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO

 do i1=1,4
    if(OutIntCont(i1).EQ.0) then
       deallocate(O(i1)%p); nullify(O(i1)%p); 
       deallocate(sO(i1)%p); nullify(sO(i1)%p); 
    endif
 enddo
 deallocate(O(5)%p); nullify(O(5)%p); 
 do i1=1,4
    if(RhsIntCont(i1).EQ.0)then
       deallocate(R(i1)%p); nullify(R(i1)%p); 
       deallocate(sR(i1)%p); nullify(sR(i1)%p); 
    endif
 enddo
 deallocate(R(5)%p); nullify(R(5)%p); 

ELSE !not RHS contraction
 ! set up Output contraction loops
 do i1=1,4
    if(OutIntCont(i1).NE.0) then
       O(i1)%p => INTEGRAL(OutIntCont(i1)); sO(i1)%p => startI(OutIntCont(i1))
    else
       nullify(O(i1)%p); allocate(O(i1)%p); O(i1)%p = 1; 
       nullify(sO(i1)%p); allocate(sO(i1)%p); sO(i1)%p = 0
    endif
 enddo
 nullify(O(5)%p); allocate(O(5)%p); O(5)%p = ideriv; 
 IF(AddToIntegral)THEN
  DO iA=1,nA
     INTEGRAL(1) = iA
     DO iB=1,nB
        INTEGRAL(2)= iB
        DO iC=1,nC
           INTEGRAL(3)= iC
           DO iD=1,nD
              INTEGRAL(4)= iD
              int = TMPintegral(A,B,C,D)
              do imat=0,(nmat-1)*nderiv,nderiv
               IF(sameLHSaos.AND.(startA.NE.startB))THEN
                 IF (sameODs.AND.((startA.NE.startC).OR.(startB.NE.startD))) THEN
                    IF (sameRHSaos.AND.(startC.NE.startD))THEN
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int  
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BACD
                       startI(1) = startintegral(2); startI(2) = startintegral(1);                            
                       OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                       startI(1) = startintegral(1); startI(2) = startintegral(2); 
                       startI(3) = startintegral(4); startI(4) = startintegral(3);
                       OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BADC
                       startI(1) = startintegral(2); startI(2) = startintegral(1); 
                       OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                       startI(1) = startintegral(3); startI(2) = startintegral(4); 
                       startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iD; INTEGRAL(2) = iC;                                     !DCAB
                       startI(1) = startintegral(4); startI(2) = startintegral(3); 
                       OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iB; INTEGRAL(4) = iA  !CDBA
                       startI(1) = startintegral(3); startI(2) = startintegral(4); 
                       startI(3) = startintegral(2); startI(4) = startintegral(1); 
                       OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iD; INTEGRAL(2) = iC;                                     !DCBA
                       startI(1) = startintegral(4); startI(2) = startintegral(3); 
                       OUT_DCBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_DCBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD; !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ELSE  !startC=startD no CD permutation
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BACD
                       startI(1) = startintegral(2); startI(2) = startintegral(1);                            
                       OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                       startI(1) = startintegral(3); startI(2) = startintegral(4);                            
                       startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(3) = iB; INTEGRAL(4) = iA                                      !CDBA
                       startI(3) = startintegral(2); startI(4) = startintegral(1); 
                       OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ENDIF
                 ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN !no OD permutation
                    IF (sameRHSaos.AND.(startC.NE.startD))THEN
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BACD
                       startI(1) = startintegral(2); startI(2) = startintegral(1); 
                       OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                       startI(1) = startintegral(1); startI(2) = startintegral(2); 
                       startI(3) = startintegral(4); startI(4) = startintegral(3); 
                       OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BADC
                       startI(1) = startintegral(2); startI(2) = startintegral(1); 
                       OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ELSE  !startC = startD   no OD or CD permutation
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                            &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA; !INTEGRAL(3) = iC; INTEGRAL(4) = iD  !BACD
                       startI(1) = startintegral(2); startI(2) = startintegral(1); 
                       !                           startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                            &OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                    ENDIF
                 ENDIF
               ELSE !NOT (startA NE startB)     !no AB symmetry
                 IF (sameODs.AND.((startA.NE.startC).OR.(startB.NE.startD))) THEN
                    IF (sameRHSaos.AND.(startC.NE.startD))THEN
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                       startI(3) = startintegral(4); startI(4) = startintegral(3); 
                       OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                       startI(1) = startintegral(3); startI(2) = startintegral(4); 
                       startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat )+int
                       INTEGRAL(1) = iD; INTEGRAL(2) = iC; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !DCAB
                       startI(1) = startintegral(4); startI(2) = startintegral(3); 
                       OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  &
                            &OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ELSE !startC.NE.startD     !no CD symmetry 
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                            &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                       startI(1) = startintegral(3); startI(2) = startintegral(4); 
                       startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                            &OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ENDIF
                 ELSE   !not ((startA.NE.startC).OR.(startB.NE.startD)) ! no od
                    IF (sameRHSaos.AND.(startC.NE.startD))THEN
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                            &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+int
                       INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                       startI(3) = startintegral(4); startI(4) = startintegral(3); 
                       OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                            &OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+int
                       INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ELSE !no permutational symmetry         ! only ABCD
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  &
                            &OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat)+int
                    ENDIF
                 ENDIF
                ENDIF
              enddo
           ENDDO !iD
        ENDDO !iC
     ENDDO !iB
  ENDDO !iA
 ELSE !no addintegrals
  DO iA=1,nA
     INTEGRAL(1) = iA
     DO iB=1,nB
        INTEGRAL(2)= iB
        DO iC=1,nC
           INTEGRAL(3)= iC
           DO iD=1,nD
              INTEGRAL(4)= iD
              int = TMPintegral(A,B,C,D)
              do imat=0,(nmat-1)*nderiv,nderiv
              IF(sameLHSaos.AND.(startA.NE.startB))THEN
                 IF (sameODs.AND.((startA.NE.startC).OR.(startB.NE.startD))) THEN
                    IF (sameRHSaos.AND.(startC.NE.startD))THEN
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int  
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BACD
                       startI(1) = startintegral(2); startI(2) = startintegral(1);                            
                       OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                       startI(1) = startintegral(1); startI(2) = startintegral(2); 
                       startI(3) = startintegral(4); startI(4) = startintegral(3);
                       OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BADC
                       startI(1) = startintegral(2); startI(2) = startintegral(1); 
                       OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                       startI(1) = startintegral(3); startI(2) = startintegral(4); 
                       startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iD; INTEGRAL(2) = iC;                                     !DCAB
                       startI(1) = startintegral(4); startI(2) = startintegral(3); 
                       OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iB; INTEGRAL(4) = iA  !CDBA
                       startI(1) = startintegral(3); startI(2) = startintegral(4); 
                       startI(3) = startintegral(2); startI(4) = startintegral(1); 
                       OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iD; INTEGRAL(2) = iC;                                     !DCBA
                       startI(1) = startintegral(4); startI(2) = startintegral(3); 
                       OUT_DCBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD; !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ELSE  !startC=startD no CD permutation
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BACD
                       startI(1) = startintegral(2); startI(2) = startintegral(1);                            
                       OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                       startI(1) = startintegral(3); startI(2) = startintegral(4);                            
                       startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(3) = iB; INTEGRAL(4) = iA                                      !CDBA
                       startI(3) = startintegral(2); startI(4) = startintegral(1); 
                       OUT_CDBA(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ENDIF
                 ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN !no OD permutation
                    IF (sameRHSaos.AND.(startC.NE.startD))THEN
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BACD
                       startI(1) = startintegral(2); startI(2) = startintegral(1); 
                       OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                       startI(1) = startintegral(1); startI(2) = startintegral(2); 
                       startI(3) = startintegral(4); startI(4) = startintegral(3); 
                       OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA;                                     !BADC
                       startI(1) = startintegral(2); startI(2) = startintegral(1); 
                       OUT_BADC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ELSE  !startC = startD   no OD or CD permutation
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  int
                       INTEGRAL(1) = iB; INTEGRAL(2) = iA; !INTEGRAL(3) = iC; INTEGRAL(4) = iD  !BACD
                       startI(1) = startintegral(2); startI(2) = startintegral(1); 
                       !                           startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_BACD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                    ENDIF
                 ENDIF
              ELSE !NOT (startA NE startB)     !no AB symmetry
                 IF (sameODs.AND.((startA.NE.startC).OR.(startB.NE.startD))) THEN
                    IF (sameRHSaos.AND.(startC.NE.startD))THEN
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                       startI(3) = startintegral(4); startI(4) = startintegral(3); 
                       OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                       startI(1) = startintegral(3); startI(2) = startintegral(4); 
                       startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iD; INTEGRAL(2) = iC; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !DCAB
                       startI(1) = startintegral(4); startI(2) = startintegral(3); 
                       OUT_DCAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat ) =  int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ELSE !startC.NE.startD     !no CD symmetry 
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  int
                       INTEGRAL(1) = iC; INTEGRAL(2) = iD; INTEGRAL(3) = iA; INTEGRAL(4) = iB  !CDAB
                       startI(1) = startintegral(3); startI(2) = startintegral(4); 
                       startI(3) = startintegral(1); startI(4) = startintegral(2); 
                       OUT_CDAB(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  int
                       INTEGRAL(1) = iA; INTEGRAL(2) = iB; INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(1) = startintegral(1); startI(2) = startintegral(2); !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ENDIF
                 ELSE   !not ((startA.NE.startC).OR.(startB.NE.startD)) ! no od
                    IF (sameRHSaos.AND.(startC.NE.startD))THEN
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  int
                       INTEGRAL(3) = iD; INTEGRAL(4) = iC  !ABDC
                       startI(3) = startintegral(4); startI(4) = startintegral(3); 
                       OUT_ABDC(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  int
                       INTEGRAL(3) = iC; INTEGRAL(4) = iD  !Back to default ABCD
                       startI(3) = startintegral(3); startI(4) = startintegral(4); !Back to default ABCD
                    ELSE !no permutational symmetry         ! only ABCD
                       OUT_ABCD(O(1)%p+sO(1)%p,O(2)%p+sO(2)%p,O(3)%p+sO(3)%p,O(4)%p+sO(4)%p,O(5)%p+imat) =  int
                    ENDIF
                 ENDIF
              ENDIF
              enddo
           ENDDO !iD
        ENDDO !iC
     ENDDO !iB
  ENDDO !iA
 ENDIF
 do i1=1,4
    if(OutIntCont(i1).EQ.0) then
       deallocate(O(i1)%p); nullify(O(i1)%p)
       deallocate(sO(i1)%p); nullify(sO(i1)%p)
    endif
 enddo
 deallocate(O(5)%p); nullify(O(5)%p)
ENDIF

END SUBROUTINE generaldistributePQ3

SUBROUTINE BuildTheta(theta,P,Q,Input)
implicit none
TYPE(Overlap),intent(in)       :: P,Q
Type(IntegralInput),intent(in) :: Input
Real(realk),intent(inout)      :: theta(P%orbital1%totOrbitals,P%orbital2%totOrbitals,&
     &                                  Q%orbital1%totOrbitals,Q%orbital2%totOrbitals,&
     &                                  Q%nPasses)
!
Integer :: idmat,iPassQ,iAngA,iAngB,iAngC,iAngD,iAngP,iAngQ,iA,iB,iC,iD
Real(realk),pointer :: LAC(:,:,:,:),LAD(:,:,:,:),LBC(:,:,:,:),LBD(:,:,:,:)
Real(realk),pointer :: RAC(:,:,:,:),RAD(:,:,:,:),RBC(:,:,:,:),RBD(:,:,:,:)
Integer :: nA,nB,nC,nD,ndmat,nPass,localA,localB,localC,localD
Integer :: startA,startB,startC,startD,endA,endB,endC,endD

IF (Input%NDMAT_RHS.NE.Input%NDMAT_LHS) CALL LSQUIT('Error in BuildTheta.Dmat inconsitency!',-1)
IF (P%nPasses.GT.1) CALL LSQUIT('Error in BuildTheta.PassP inconsitency!',-1)

!Dimensions
ndmat = Input%NDMAT_RHS
nA    = P%orbital1%totOrbitals
nB    = P%orbital2%totOrbitals
nC    = Q%orbital1%totOrbitals
nD    = Q%orbital2%totOrbitals
nPass = Q%nPasses

!Allocations
call mem_alloc(LAC,nA,nC,nPass,ndmat)
call mem_alloc(LAD,nA,nD,nPass,ndmat)
call mem_alloc(LBC,nB,nC,nPass,ndmat)
call mem_alloc(LBD,nB,nD,nPass,ndmat)
call mem_alloc(RAC,nA,nC,nPass,ndmat)
call mem_alloc(RAD,nA,nD,nPass,ndmat)
call mem_alloc(RBC,nB,nC,nPass,ndmat)
call mem_alloc(RBD,nB,nD,nPass,ndmat)

!Set up local/reduced density-matrices
DO idmat=1,ndmat
 DO iPassQ=1,Q%nPasses
!  C contributions
   localC = 0
   DO iAngC=1,Q%orbital1%nAngmom
     startC = Q%orbital1%startOrbital(iAngC,iPassQ)
     endC   = startC + Q%orbital1%nOrbitals(iAngC) - 1
     DO iC=startC,endC
       localC = localC+1
!      AC-density
       localA = 0
       DO iAngA=1,P%orbital1%nAngmom
         startA = P%orbital1%startOrbital(iAngA,1)
         endA   = startA + P%orbital1%nOrbitals(iAngA) - 1
         DO iA=startA,endA
           localA = localA+1
           LAC(localA,localC,iPassQ,idmat) = Input%DMAT_LHS(iA,iC,idmat)
           RAC(localA,localC,iPassQ,idmat) = Input%DMAT_RHS(iA,iC,idmat)
         ENDDO
       ENDDO
!      BC-density
       localB = 0
       DO iAngB=1,P%orbital2%nAngmom
         startB = P%orbital2%startOrbital(iAngB,1)
         endB   = startB + P%orbital2%nOrbitals(iAngB) - 1
         DO iB=startB,endB
           localB = localB+1
           LBC(localB,localC,iPassQ,idmat) = Input%DMAT_LHS(iB,iC,idmat)
           RBC(localB,localC,iPassQ,idmat) = Input%DMAT_RHS(iB,iC,idmat)
         ENDDO
       ENDDO
     ENDDO
   ENDDO

!  D contributions
   localD = 0
   DO iAngD=1,Q%orbital2%nAngmom
     startD = Q%orbital2%startOrbital(iAngD,iPassQ)
     endD   = startD + Q%orbital2%nOrbitals(iAngD) - 1
     DO iD=startD,endD
       localD = localD+1
!      AD-density
       localA = 0
       DO iAngA=1,P%orbital1%nAngmom
         startA = P%orbital1%startOrbital(iAngA,1)
         endA   = startA + P%orbital1%nOrbitals(iAngA) - 1
         DO iA=startA,endA
           localA = localA+1
           LAD(localA,localD,iPassQ,idmat) = Input%DMAT_LHS(iA,iD,idmat)
           RAD(localA,localD,iPassQ,idmat) = Input%DMAT_RHS(iA,iD,idmat)
         ENDDO
       ENDDO
!      BD-density
       localB = 0
       DO iAngB=1,P%orbital2%nAngmom
         startB = P%orbital2%startOrbital(iAngB,1)
         endB   = startB + P%orbital2%nOrbitals(iAngB) - 1
         DO iB=startB,endB
           localB = localB+1
           LBD(localB,localD,iPassQ,idmat) = Input%DMAT_LHS(iB,iD,idmat)
           RBD(localB,localD,iPassQ,idmat) = Input%DMAT_RHS(iB,iD,idmat)
         ENDDO
       ENDDO
     ENDDO
   ENDDO
 ENDDO
ENDDO

! Build theta for given pass
idmat=1
  DO iPassQ=1,nPass
    DO iD=1,nD
      DO iC=1,nC
        DO iB=1,nB
          DO iA=1,nA
            theta(iA,iB,iC,iD,iPassQ) = LAC(iA,iC,iPassQ,idmat)*RBD(iB,iD,iPassQ,idmat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

DO idmat=2,ndmat
  DO iPassQ=1,nPass
    DO iD=1,nD
      DO iC=1,nC
        DO iB=1,nB
          DO iA=1,nA
            theta(iA,iB,iC,iD,iPassQ) = theta(iA,iB,iC,iD,iPassQ) + LAC(iA,iC,iPassQ,idmat)*RBD(iB,iD,iPassQ,idmat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO

IF (Input%SameLHSaos.AND.Input%SameRHSaos.AND.Input%SameODs) THEN
  CALL LSQUIT('Error in BuildTheta. Not implemented case TTT',-1)
ELSE IF (Input%SameLHSaos.AND.Input%SameRHSaos) THEN
! Integral-type loop
  DO idmat=1,ndmat
    DO iPassQ=1,Q%nPasses
      DO iAngQ=1,Q%nAngmom
        iAngC = Q%indexAng1(iAngQ)
        iAngD = Q%indexAng2(iAngQ)
        startD = Q%orbital2%startOrbital(iAngD,iPassQ)
        localD = Q%orbital2%startLocOrb(iAngD)
        endD   = localD + Q%orbital2%nOrbitals(iAngD) - 1
        startC = Q%orbital1%startOrbital(iAngC,iPassQ)
        localC = Q%orbital1%startLocOrb(iAngC)
        endC   = localC + Q%orbital1%nOrbitals(iAngC) - 1
        IF (startC.NE.startD) THEN
          DO iD=localD,endD
            DO iC=localC,endC
              DO iAngP=1,P%nAngmom
                iAngA = P%indexAng1(iAngP)
                iAngB = P%indexAng2(iAngP)
                startB = P%orbital2%startOrbital(iAngB,1)
                localB = P%orbital2%startLocOrb(iAngB)
                endB   = localB + P%orbital2%nOrbitals(iAngB) - 1
                startA = P%orbital1%startOrbital(iAngA,1)
                localA = P%orbital1%startLocOrb(iAngA)
                endA   = localA + P%orbital1%nOrbitals(iAngA) - 1
                IF (startA.NE.startB) THEN
                  DO iB=localB,endB
                    DO iA=localA,endA
                      theta(iA,iB,iC,iD,iPassQ) = theta(iA,iB,iC,iD,iPassQ) + LAD(iA,iD,iPassQ,idmat)*RBC(iB,iC,iPassQ,idmat)
                      theta(iA,iB,iC,iD,iPassQ) = theta(iA,iB,iC,iD,iPassQ) + LBC(iB,iC,iPassQ,idmat)*RAD(iA,iD,iPassQ,idmat)
                      theta(iA,iB,iC,iD,iPassQ) = theta(iA,iB,iC,iD,iPassQ) + LBD(iB,iD,iPassQ,idmat)*RAC(iA,iC,iPassQ,idmat)
                    ENDDO !A
                  ENDDO !B
                ELSE
                  DO iB=localB,endB
                    DO iA=localA,endA
                      theta(iA,iB,iC,iD,iPassQ) = theta(iA,iB,iC,iD,iPassQ) + LAD(iA,iD,iPassQ,idmat)*RBC(iB,iC,iPassQ,idmat)
                    ENDDO !A
                  ENDDO !B
                ENDIF
              ENDDO !AngP
            ENDDO !C
          ENDDO !D
        ELSE ! startC .EQ. startD
          DO iD=localD,endD
            DO iC=localC,endC
              DO iAngP=1,P%nAngmom
                iAngA = P%indexAng1(iAngP)
                iAngB = P%indexAng2(iAngP)
                startB = P%orbital2%startOrbital(iAngB,1)
                localB = P%orbital2%startLocOrb(iAngB)
                endB   = localB + P%orbital2%nOrbitals(iAngB) - 1
                startA = P%orbital1%startOrbital(iAngA,1) 
                localA = P%orbital1%startLocOrb(iAngA) 
                endA   = localA + P%orbital1%nOrbitals(iAngA) - 1
                IF (startA.NE.startB) THEN
                  DO iB=localB,endB
                    DO iA=localA,endA 
                      theta(iA,iB,iC,iD,iPassQ) = theta(iA,iB,iC,iD,iPassQ) + LBC(iB,iC,iPassQ,idmat)*RAD(iA,iD,iPassQ,idmat)
                    ENDDO !A
                  ENDDO !B
                ENDIF
              ENDDO !AngP
            ENDDO !C
          ENDDO !D
        ENDIF !startC startD
      ENDDO !AngQ
    ENDDO !PassQ
  ENDDO !ndmat
ELSE IF (Input%SameLHSaos) THEN
  CALL LSQUIT('Error in BuildTheta. Not implemented case TFF',-1)
ELSE IF (Input%SameRHSaos) THEN
  CALL LSQUIT('Error in BuildTheta. Not implemented case FTF',-1)
ELSEIF (Input%SameODs) THEN
  CALL LSQUIT('Error in BuildTheta. Not implemented case FFT',-1)
ENDIF
  
!Integral-type loop
!DO idmat=1,Input%NDMAT_RHS
!  DO iPassQ=1,Q%nPasses
!    DO iAngQ=1,Q%nAngmom
!      iAngC = Q%indexAng1(iAngQ)
!      iAngD = Q%indexAng2(iAngQ)
!      startD = Q%orbital1%startOrbital(iAngD,iPassQ)
!      endD   = startD + Q%orbital1%nOrbitals(iAngD) - 1
!      startC = Q%orbital1%startOrbital(iAngC,iPassQ)
!      endC   = startC + Q%orbital1%nOrbitals(iAngC) - 1
!      DO iD=startD,endD
!        DO iC=startC,endC
!          DO iAngP=1,Q%nAngmom
!            iAngA = P%indexAng1(iAngP)
!            iAngB = P%indexAng2(iAngP)
!            startB = P%orbital1%startOrbital(iAngB,1)
!            endB   = startB + P%orbital1%nOrbitals(iAngB) - 1
!            startA = P%orbital1%startOrbital(iAngA,1)
!            endA   = startA + P%orbital1%nOrbitals(iAngA) - 1
!            DO iB=startB,endB
!              DO iA=startA,endA
!
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
!    ENDDO
!  ENDDO
!ENDDO

call mem_dealloc(LAC)
call mem_dealloc(LAD)
call mem_dealloc(LBC)
call mem_dealloc(LBD)
call mem_dealloc(RAC)
call mem_dealloc(RAD)
call mem_dealloc(RBC)
call mem_dealloc(RBD)

END SUBROUTINE BuildTheta

SUBROUTINE BuildDerivIntegral(integral,QPmat,P,Q,samePQ,Input)
implicit none
TYPE(Overlap),intent(in)       :: P,Q
Logical,intent(in)             :: samePQ
Type(IntegralInput),intent(in) :: Input
Real(realk),intent(inout)      :: integral(P%orbital1%totOrbitals,P%orbital2%totOrbitals,&
     &                                     Q%orbital1%totOrbitals,Q%orbital2%totOrbitals,&
     &                                     Q%nPasses,Input%nDerivP)
Real(realk),intent(in)         :: QPmat(Q%totOrbitals,P%totOrbitals)
!
Integer :: iPassP,iPassQ,iAngA,iAngB,iAngC,iAngD,iA,iB,iC,iD
Integer :: nA,nB,nC,nD,ndmat,nPass,localA,localB,localC,localD
Integer :: startA,startB,startC,startD,atomA,atomB,atomC,atomD
Integer :: iOrbitalP,iOrbitalQ,iOrbP,iOrbQ,endOrbP,endOrbQ,endAngmomQ,iAngmomP,iAngmomQ
Integer :: iContA,iContB,iContC,iContD,nAngA,nAngB,nAngC,nAngD,iDerivP,iDerivQ
Integer :: nContA,nContB,nContC,nContD,nDeriv,iContP,iContQ,nDerivP,nDerivQ,nOrbP,nOrbQ

IF (P%nPasses.GT.1) CALL LSQUIT('Error in BuildDerivIntegral. PassP>1 not implemented!',-1)
IF (Input%nDerivQ.GT.1) CALL LSQUIT('Error in BuildDerivIntegral. nDerivQ>1 not implemented!',-1)

!Dimensions
ndmat   = Input%NDMAT_RHS
nA      = P%orbital1%totOrbitals
nB      = P%orbital2%totOrbitals
nC      = Q%orbital1%totOrbitals
nD      = Q%orbital2%totOrbitals
nDerivP = Input%nDerivP
nDerivQ = Input%nDerivQ
nDeriv  = nDerivP*nDerivQ
nPass   = Q%nPasses

! Integral-type loop
iOrbitalP = 1
DO iAngmomP=1,P%nAngmom
  DO iDerivP=1,NDerivP
    DO iPassP = 1, P%nPasses
      nOrbP  = P%nOrbitals(iAngmomP)
      endOrbP = iOrbitalP+nOrbP-1

      CALL getOverlapInfo(nContA,nContB,nAngA,nAngB,startA,startB,localA,localB,atomA,atomB,P,iAngmomP,iPassP)
      iOrbP=iOrbitalP-1+(iPassP-1)*nContA*nContB*nAngA*nAngB
      iContP=0
      DO iContB=1,nContB
         DO iContA=1,nContA
            iContP=iContP+1
            DO iAngB=1,nAngB
               iB=localB-1+iAngB+(iContB-1)*nAngB
               DO iAngA=1,nAngA
                  iA=localA-1+iAngA+(iContA-1)*nAngA
                  iOrbP=iOrbP+1
                  iOrbitalQ = 1
                  endAngmomQ = Q%nAngmom
                  IF (samePQ) endAngmomQ = iAngmomP
                  nderivQ=input%nderivQ
                  DO iAngmomQ=1,endAngmomQ
                    DO iderivQ =1,nderivQ
                     DO iPassQ = 1, Q%nPasses
                        nOrbQ = Q%nOrbitals(iAngmomQ)
                        endOrbQ = iOrbitalQ+nOrbQ-1
      
                        CALL getOverlapInfo(nContC,nContD,nAngC,nAngD,startC,startD,localC,localD,atomC,atomD,Q,iAngmomQ,iPassQ)
                        iOrbQ=iOrbitalQ-1+(iPassQ-1)*nContC*nContD*nAngC*nAngD
                        DO iContD=1,nContD
                           DO iContC=1,nContC
                              iContQ = (iContD-1)*nContC + iContC
                              DO iAngD=1,nAngD
                                 iD=localD-1+iAngD+(iContD-1)*nAngD
                                 DO iAngC=1,nAngC
                                    iC=localC-1+iAngC+(iContC-1)*nAngC
                                    iOrbQ=iOrbQ+1
                                    integral(iA,iB,iC,iD,iPassQ,iDerivP) = QPmat(iOrbQ,iOrbP)
                                 ENDDO !iAngC
                              ENDDO !iAngD
                            ENDDO !iContC
                          ENDDO !iContD
                     ENDDO !iPassQ
                     iOrbitalQ = iOrbitalQ + nOrbQ
                    ENDDO !iderivQ
                  ENDDO !iAngmomQ
               ENDDO !iAngA
            ENDDO !iAngB
         ENDDO !iContA
      ENDDO !iContB
    ENDDO !iPassP
    iOrbitalP = iOrbitalP + nOrbP
  ENDDO !iDerivP
ENDDO !iAngmomP

END SUBROUTINE BuildDerivIntegral

END MODULE thermite_distribute
