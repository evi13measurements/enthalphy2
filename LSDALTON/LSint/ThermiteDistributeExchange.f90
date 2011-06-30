!> @file 
!> Contains integral distribution routines that places calculated integrals in the proper output
!> Thermite integral distribution module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE thermite_distribute_Exchange
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
  use thermite_distribute_ExchangeW
  use thermite_distribute, only: derivativeInfo, initderivativeoverlapinfo, &
       &printderivativeoverlapinfo, getoverlapinfo, freederivativeoverlapinfo

SAVE
CONTAINS
!> \brief distribute Integral to exchange contribution for LinK or exchange gradient
!> \author T. Kjaergaard and S. Reine
!> \date 2010-04-23
!> \param PQ contain info about the overlap distributions
!> \param QPmat3 matrix containing calculated integrals
!> \param Dsym Flag indicating if the density matrix is symmetric
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param LSOutput output structure contains exchange matrix in lstensor format
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributeExchangeLinK(PQ,QPmat3,Dsym,dimQ,dimP,Input,LSOutput,LUPRI,IPRINT)
! We assume that the density matrix and thereby the Fock matrix is symmetrix
implicit none
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: LSOutput
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
REAL(REALK),intent(in)          :: QPMAT3(dimQ,dimP)
Logical,intent(in)              :: Dsym
real(realk),pointer :: tmpD1(:,:),tmpD2(:,:,:,:,:)
!
IF(input%do_gradient)THEN
   call lstensorExchangeGrad(PQ,QPmat3,dimQ,dimP,Input,input%LST_DRHS,LSOUTPUT%resultmat2,LUPRI,IPRINT)
ELSE
!$OMP CRITICAL (distributeLinK) 
   CALL lstensorExchange(PQ,QPmat3,dimQ,dimP,Input,input%LST_DRHS,LSOutput%ResultMat2,Dsym,LUPRI,IPRINT)
!$OMP END CRITICAL (distributeLinK) 
ENDIF
END SUBROUTINE distributeExchangeLinK

!> \brief distribute Integral to exchange contribution for LinK
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param PQ contain info about the overlap distributions
!> \param QPmat2 matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param Dmat the density matrix in a lstensor format
!> \param Kmat the exchange matrix in a lstensor format
!> \param Dsym Flag indicating if the density matrix is symmetric
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE lstensorExchange(PQ,QPmat2,dimQ,dimP,Input,Dmat,Kmat,Dsym,LUPRI,IPRINT)
implicit none
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
TYPE(lstensor),intent(in)       :: Dmat
TYPE(lstensor),intent(inout)    :: Kmat
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
Logical                         :: Dsym
!
Integer :: nmat,nAngmomP,nAngmomQ,nOrbP,iAngmomP,iOrbitalP
Integer :: dimA,dimB,dimC,dimD,SAA,SBB,SCC,SDD,atomA,atomB,atomC,atomD
Integer :: iA,iB,iC,iD,nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ
Integer :: iContQ,nContC,nContD,nAngC,nAngD,startC,startD,batchA
Integer :: iOrbQ,idmat,startA1,startB1,startC1,startD1,batchB,batchC,batchD
real(realk) :: THRESHOLD,TMP
integer :: PASSPOFFSET,PASSQOFFSET,iPassP,iPassQ,sA,sB,sC,sD,n1,n2,nPassP,nPassQ
logical :: samePQ,SameRHSaos,SameLHSaos,SameODs
Integer :: Kac,Kbc,Kad,Kbd,Kca,Kcb,Kda,Kdb,Dbd,Dad,Dac,Dbc,Dca,Dcb,Dda,Ddb,iopt
integer,pointer :: option(:)

nMAT=Kmat%nmat
nAngmomP = PQ%P%p%nAngmom
nAngmomQ = PQ%Q%p%nAngmom
nPassP = PQ%P%p%nPasses
nPassQ = PQ%Q%p%nPasses
THRESHOLD = INPUT%CS_THRESHOLD
!intthreshold = INPUT%CS_THRESHOLD
SamePQ = PQ%samePQ
SameLHSaos = INPUT%SameLHSaos
SameRHSaos = INPUT%SameRHSaos
SameODs = INPUT%SameODs
allocate(option(nAngmomP*nPassP*nAngmomQ*nPassQ))
!------------
iopt=0
DO iAngmomP=1,nAngmomP
   iA = PQ%P%p%indexAng1(iAngmomP)
   iB = PQ%P%p%indexAng2(iAngmomP)
   DO iPassP = 1, nPassP
      startA  = PQ%P%p%orbital1%startOrbital(iA,iPassP)
      startB  = PQ%P%p%orbital2%startOrbital(iB,iPassP)
      endAngmomQ = nAngmomQ
      IF (samePQ) endAngmomQ = iAngmomP
      DO iAngmomQ=1,endAngmomQ
         iC = PQ%Q%p%indexAng1(iAngmomQ)
         iD = PQ%Q%p%indexAng2(iAngmomQ)
         DO iPassQ = 1, nPassQ
            iopt = iopt+1
            startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
            startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
            IF((startA.NE.startB).AND.SameLHSaos)THEN
               IF (((startA.NE.startC).OR.(startB.NE.startD)).AND.SameODs ) THEN
                  IF ((startC.NE.startD).AND.SameRHSaos)THEN
                     IF(Dsym)THEN
                        option(iopt) = 1
                     ELSE
                        option(iopt) = 2
                     ENDIF
                  ELSE ! NOT CD SYMMETRY
                     IF(Dsym)THEN
                        option(iopt) = 3
                     ELSE
                        option(iopt) = 4
                     ENDIF
                  ENDIF
               ELSE !NOT ((startA.NE.startC).AND.(startB.NE.startD)) THEN
                  IF ((startC.NE.startD).AND.SameRHSaos)THEN
                        option(iopt) = 5
                  ELSE
                        option(iopt) = 6
                  ENDIF
               ENDIF
            ELSE ! (startA = startB)
               IF (((startA.NE.startC).OR.(startB.NE.startD)).AND.SameODs) THEN
                  IF ((startC.NE.startD).AND.SameRHSaos)THEN
                     IF(Dsym)THEN
                        option(iopt) = 7
                     ELSE
                        option(iopt) = 8
                     ENDIF
                  ELSE
                        option(iopt) = 9
                  ENDIF
               ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN
                  IF ((startC.NE.startD).AND.SameRHSaos)THEN
                        option(iopt) = 10
                  ELSE
                        option(iopt) = 11
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!-------------

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
iC = PQ%Q%p%indexAng1(nAngmomQ)
dimC = 0
DO SCC=1,iC
   dimC = dimC+PQ%Q%p%orbital1%nOrbComp(SCC)*PQ%Q%p%orbital1%nContracted(SCC)
ENDDO
iD = PQ%Q%p%indexAng2(nAngmomQ)
dimD = 0
DO SDD=1,iD
   dimD = dimD+PQ%Q%p%orbital2%nOrbComp(SDD)*PQ%Q%p%orbital2%nContracted(SDD)
ENDDO
!
iOrbitalP = 1
iopt=0
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
  batchA = PQ%P%p%orbital1%batch(iPassP)
  batchB = PQ%P%p%orbital2%batch(iPassP)
  iOrbP=iOrbitalP-1+PASSPOFFSET
  iOrbitalQ = 1
  endAngmomQ = nAngmomQ
  IF (samePQ) endAngmomQ = iAngmomP
  DO iAngmomQ=1,endAngmomQ
   sC = 0
   DO SCC=2,PQ%Q%p%indexAng1(iAngmomQ)
    sC = sC+PQ%Q%p%orbital1%nOrbComp(SCC-1)*PQ%Q%p%orbital1%nContracted(SCC-1)
   ENDDO
   sD = 0
   DO SDD=2,PQ%Q%p%indexAng2(iAngmomQ)
      sD = sD+PQ%Q%p%orbital2%nOrbComp(SDD-1)*PQ%Q%p%orbital2%nContracted(SDD-1)
   ENDDO
   nOrbQ = PQ%Q%p%nOrbitals(iAngmomQ)
   iC = PQ%Q%p%indexAng1(iAngmomQ)
   iD = PQ%Q%p%indexAng2(iAngmomQ)
   nContC = PQ%Q%p%orbital1%nContracted(iC)
   nContD = PQ%Q%p%orbital2%nContracted(iD)
   nAngC = PQ%Q%p%orbital1%nOrbComp(iC)
   nAngD = PQ%Q%p%orbital2%nOrbComp(iD)
   DO iPassQ = 1, PQ%Q%p%nPasses
    iopt = iopt+1
    PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
    iOrbQ=iOrbitalQ-1+PASSQOFFSET
    TMP=0.D0
    DO n2=1,nOrbP
       DO n1=1,nAngC*nAngD*nContC*nContD
          TMP = TMP+ABS(QPmat2(iOrbQ+n1,iOrbP+n2))
       ENDDO
    ENDDO
    IF(ABS(TMP).LT.THRESHOLD)CYCLE
    startC  = PQ%Q%p%orbital1%startOrbital(iC,iPassQ)
    startD  = PQ%Q%p%orbital2%startOrbital(iD,iPassQ)
    atomC = PQ%Q%p%orbital1%atom(iPassQ)
    atomD = PQ%Q%p%orbital2%atom(iPassQ)
    batchC = PQ%Q%p%orbital1%batch(iPassQ)
    batchD = PQ%Q%p%orbital2%batch(iPassQ)
    KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
    DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
    SELECT CASE(option(iopt))
    CASE(1) 
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KAD = KMAT%INDEX(ATOMA,ATOMD,1,1)
       DBC = Dmat%INDEX(ATOMB,ATOMC,1,1)
       KBC = KMAT%INDEX(ATOMB,ATOMC,1,1)
       DAD = Dmat%INDEX(ATOMA,ATOMD,1,1)
       KAC = KMAT%INDEX(ATOMA,ATOMC,1,1)
       DBD = Dmat%INDEX(ATOMB,ATOMD,1,1)
       !FULL PERMUTATIONEL SYMMETRY
       !(ABCD)=(BACD)=(ABDC)=(BADC)=(CDAB)=(DCAB)=(DCBA)=(CDBA)
       CALL exchangePermutTTTSYM(&
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC} 
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,& 
            & Kmat%LSAO(KAD)%BATCH(batchA,batchD,1,1)%elms,&!K_AD=(BACD)D_{BC}
            & Dmat%LSAO(DBC)%BATCH(batchB,batchC,1,1)%elms,& 
            & Kmat%LSAO(KBC)%BATCH(batchB,batchC,1,1)%elms,&!K_BC=(ABDC)D_{AD} 
            & Dmat%LSAO(DAD)%BATCH(batchA,batchD,1,1)%elms,& 
            & Kmat%LSAO(KAC)%BATCH(batchA,batchC,1,1)%elms,&!K_AC=(BADC)D_{BD}
            & Dmat%LSAO(DBD)%BATCH(batchB,batchD,1,1)%elms,&
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,&
            & nContA,nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT)
     CASE(2)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KAD = KMAT%INDEX(ATOMA,ATOMD,1,1)
       DBC = Dmat%INDEX(ATOMB,ATOMC,1,1)
       KBC = KMAT%INDEX(ATOMB,ATOMC,1,1)
       DAD = Dmat%INDEX(ATOMA,ATOMD,1,1)
       KAC = KMAT%INDEX(ATOMA,ATOMC,1,1)
       DBD = Dmat%INDEX(ATOMB,ATOMD,1,1)
       KDB = KMAT%INDEX(ATOMD,ATOMB,1,1)
       DCA = Dmat%INDEX(ATOMC,ATOMA,1,1)
       KCB = KMAT%INDEX(ATOMC,ATOMB,1,1)
       DDA = Dmat%INDEX(ATOMD,ATOMA,1,1)
       KCA = KMAT%INDEX(ATOMC,ATOMA,1,1)
       DDB = Dmat%INDEX(ATOMD,ATOMB,1,1)
       KDA = KMAT%INDEX(ATOMD,ATOMA,1,1)
       DCB = Dmat%INDEX(ATOMC,ATOMB,1,1)
       CALL exchangePermutTTT(&
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC}
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,&
            & Kmat%LSAO(KAD)%BATCH(batchA,batchD,1,1)%elms,&!K_AD=(BACD)D_{BC}
            & Dmat%LSAO(DBC)%BATCH(batchB,batchC,1,1)%elms,& 
            & Kmat%LSAO(KBC)%BATCH(batchB,batchC,1,1)%elms,&!K_BC=(ABDC)D_{AD} 
            & Dmat%LSAO(DAD)%BATCH(batchA,batchD,1,1)%elms,& 
            & Kmat%LSAO(KAC)%BATCH(batchA,batchC,1,1)%elms,&!K_AC=(BADC)D_{BD}
            & Dmat%LSAO(DBD)%BATCH(batchB,batchD,1,1)%elms,& 
            & Kmat%LSAO(KDB)%BATCH(batchD,batchB,1,1)%elms,&!K_DB=(CDAB)D_{CA} 
            & Dmat%LSAO(DCA)%BATCH(batchC,batchA,1,1)%elms,& 
            & Kmat%LSAO(KCB)%BATCH(batchC,batchB,1,1)%elms,&!K_CB=(DCAB)D_{DA} 
            & Dmat%LSAO(DDA)%BATCH(batchD,batchA,1,1)%elms,& 
            & Kmat%LSAO(KCA)%BATCH(batchC,batchA,1,1)%elms,&!K_CA=(DCBA)D_{DB}
            & Dmat%LSAO(DDB)%BATCH(batchD,batchB,1,1)%elms,& 
            & Kmat%LSAO(KDA)%BATCH(batchD,batchA,1,1)%elms,&!K_DA=(CDBA)D_{CB}
            & Dmat%LSAO(DCB)%BATCH(batchC,batchB,1,1)%elms,& 
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,&
            & iOrbP,iOrbQ,nContB,nContA,nAngB,nAngA,nContD,&
            & nContC,nAngD,nAngC,QPmat2,nMAT,startA.EQ.startC,startB.EQ.startD)      
    CASE(3)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KAD = KMAT%INDEX(ATOMA,ATOMD,1,1)
       DBC = Dmat%INDEX(ATOMB,ATOMC,1,1)
       CALL exchangePermutTFTsym(&   ! 1 2 5 8
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC} 
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,& 
            & Kmat%LSAO(KAD)%BATCH(batchA,batchD,1,1)%elms,&!K_AD=(BACD)D_{BC}
            & Dmat%LSAO(DBC)%BATCH(batchB,batchC,1,1)%elms,& 
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,nContA,&
            & nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT)
    CASE(4)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KAD = KMAT%INDEX(ATOMA,ATOMD,1,1)
       DBC = Dmat%INDEX(ATOMB,ATOMC,1,1)
       KDB = KMAT%INDEX(ATOMD,ATOMB,1,1)
       DCA = Dmat%INDEX(ATOMC,ATOMA,1,1)
       KDA = KMAT%INDEX(ATOMD,ATOMA,1,1)
       DCB = Dmat%INDEX(ATOMC,ATOMB,1,1)
       CALL exchangePermutTFT(&   ! 1 2 5 8
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC}
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,&
            & Kmat%LSAO(KAD)%BATCH(batchA,batchD,1,1)%elms,&!K_AD=(BACD)D_{BC}
            & Dmat%LSAO(DBC)%BATCH(batchB,batchC,1,1)%elms,&
            & Kmat%LSAO(KDB)%BATCH(batchD,batchB,1,1)%elms,&!K_DB=(CDAB)D_{CA} 
            & Dmat%LSAO(DCA)%BATCH(batchC,batchA,1,1)%elms,&
            & Kmat%LSAO(KDA)%BATCH(batchD,batchA,1,1)%elms,&!K_DA=(CDBA)D_{CB}
            & Dmat%LSAO(DCB)%BATCH(batchC,batchB,1,1)%elms,&
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,nContA,&
            & nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT)
    CASE(5)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KAD = KMAT%INDEX(ATOMA,ATOMD,1,1)
       DBC = Dmat%INDEX(ATOMB,ATOMC,1,1)
       KBC = KMAT%INDEX(ATOMB,ATOMC,1,1)
       DAD = Dmat%INDEX(ATOMA,ATOMD,1,1)
       KAC = KMAT%INDEX(ATOMA,ATOMC,1,1)
       DBD = Dmat%INDEX(ATOMB,ATOMD,1,1)
       CALL exchangePermutTTF(&   ! 1 2 3 4
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC} 
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,&
            & Kmat%LSAO(KAD)%BATCH(batchA,batchD,1,1)%elms,&!K_AD=(BACD)D_{BC}
            & Dmat%LSAO(DBC)%BATCH(batchB,batchC,1,1)%elms,&
            & Kmat%LSAO(KBC)%BATCH(batchB,batchC,1,1)%elms,&!K_BC=(ABDC)D_{AD} 
            & Dmat%LSAO(DAD)%BATCH(batchA,batchD,1,1)%elms,&
            & Kmat%LSAO(KAC)%BATCH(batchA,batchC,1,1)%elms,&!K_AC=(BADC)D_{BD}
            & Dmat%LSAO(DBD)%BATCH(batchB,batchD,1,1)%elms,&
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,nContA,&
            & nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT)
    CASE(6)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KAD = KMAT%INDEX(ATOMA,ATOMD,1,1)
       DBC = Dmat%INDEX(ATOMB,ATOMC,1,1)
       CALL exchangePermutTFF(&   ! 1 2
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC} 
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,&
            & Kmat%LSAO(KAD)%BATCH(batchA,batchD,1,1)%elms,&!K_AD=(BACD)D_{BC}
            & Dmat%LSAO(DBC)%BATCH(batchB,batchC,1,1)%elms,& 
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,nContA,&
            & nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT)
    CASE(7)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KBC = KMAT%INDEX(ATOMB,ATOMC,1,1)
       DAD = Dmat%INDEX(ATOMA,ATOMD,1,1)
       CALL exchangePermutFTTsym(&   ! 1 3 5 6
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC} 
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,& 
            & Kmat%LSAO(KBC)%BATCH(batchB,batchC,1,1)%elms,&!K_BC=(ABDC)D_{AD} 
            & Dmat%LSAO(DAD)%BATCH(batchA,batchD,1,1)%elms,&
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,nContA,&
            & nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT)
    CASE(8)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KBC = KMAT%INDEX(ATOMB,ATOMC,1,1)
       DAD = Dmat%INDEX(ATOMA,ATOMD,1,1)
       KDB = KMAT%INDEX(ATOMD,ATOMB,1,1)
       DCA = Dmat%INDEX(ATOMC,ATOMA,1,1)
       KCB = KMAT%INDEX(ATOMC,ATOMB,1,1)
       DDA = Dmat%INDEX(ATOMD,ATOMA,1,1)
       CALL exchangePermutFTT(&   ! 1 3 5 6
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC} 
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,&
            & Kmat%LSAO(KBC)%BATCH(batchB,batchC,1,1)%elms,&!K_BC=(ABDC)D_{AD} 
            & Dmat%LSAO(DAD)%BATCH(batchA,batchD,1,1)%elms,&
            & Kmat%LSAO(KDB)%BATCH(batchD,batchB,1,1)%elms,&!K_DB=(CDAB)D_{CA} 
            & Dmat%LSAO(DCA)%BATCH(batchC,batchA,1,1)%elms,&
            & Kmat%LSAO(KCB)%BATCH(batchC,batchB,1,1)%elms,&!K_CB=(DCAB)D_{DA} 
            & Dmat%LSAO(DDA)%BATCH(batchD,batchA,1,1)%elms,&
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,nContA,&
            & nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT,&
            & startB.EQ.startC,startB.EQ.startD)
    CASE(9)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KDB = KMAT%INDEX(ATOMD,ATOMB,1,1)
       DCA = Dmat%INDEX(ATOMC,ATOMA,1,1)
       CALL exchangePermutFFT(&   ! 1 5
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC}
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,&
            & Kmat%LSAO(KDB)%BATCH(batchD,batchB,1,1)%elms,&!K_DB=(CDAB)D_{CA} 
            & Dmat%LSAO(DCA)%BATCH(batchC,batchA,1,1)%elms,&
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,nContA,&
            & nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT)
    CASE(10)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       KBC = KMAT%INDEX(ATOMB,ATOMC,1,1)
       DAD = Dmat%INDEX(ATOMA,ATOMD,1,1)
       CALL exchangePermutFTF(&   ! 1 3
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,&!K_BD=(ABCD)D_{AC} 
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,& 
            & Kmat%LSAO(KBC)%BATCH(batchB,batchC,1,1)%elms,&!K_BC=(ABDC)D_{AD} 
            & Dmat%LSAO(DAD)%BATCH(batchA,batchD,1,1)%elms,&
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,nContA,&
            & nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT)
       CASE(11)
       KBD = KMAT%INDEX(ATOMB,ATOMD,1,1)
       DAC = Dmat%INDEX(ATOMA,ATOMC,1,1)
       CALL exchangePermutFFF(&   ! 1
            & Kmat%LSAO(KBD)%BATCH(batchB,batchD,1,1)%elms,& !K_BD=(ABCD)D_{AC}
            & Dmat%LSAO(DAC)%BATCH(batchA,batchC,1,1)%elms,&
            & dimA,dimB,dimC,dimD,sA,sB,sC,sD,iOrbP,iOrbQ,nContB,nContA,&
            & nAngB,nAngA,nContD,nContC,nAngD,nAngC,QPmat2,nMAT)
     CASE DEFAULT
        WRITE (LUPRI,'(A,I5,A)') ' Option ',option(iopt),&
             & ' not recognized in lstensorExchange'
     END SELECT
   ENDDO
   iOrbitalQ = iOrbitalQ + nOrbQ
  ENDDO
 ENDDO
 iOrbitalP = iOrbitalP + nOrbP
ENDDO
deallocate(option)

END SUBROUTINE lstensorExchange

!> \brief distribute Integral to exchange gradient
!> \author S. Reine
!> \date 2010-04-23
!> \param PQ contain info about the overlap distributions
!> \param QPmat2 matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param Dmat the density matrix in a lstensor format
!> \param Kmat the exchange gradient in a lstensor format
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE lstensorExchangeGrad(PQ,QPmat2,dimQ,dimP,Input,Dmat,Kmat,LUPRI,IPRINT)
implicit none
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
TYPE(lstensor),intent(in)       :: Dmat
TYPE(lstensor),intent(inout)    :: Kmat
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
!
Integer :: nmat,nAngmomP,nAngmomQ,nOrbP,iAngmomP,iOrbitalP,nDerivP,nDerivQ
Integer :: dimA,dimB,dimC,dimD,SAA,SBB,SCC,SDD,atomA,atomB,atomC,atomD
Integer :: iA,iB,iC,iD,nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iOrbitalQ,endAngmomQ,endAngmom,iAngmomQ,nOrbQ,endOrbQ
Integer :: iContQ,nContC,nContD,nAngC,nAngD,startC,startD,batchA
Integer :: iOrbQ,idmat,startA1,startB1,startC1,startD1,batchB,batchC,batchD
integer :: PASSPOFFSET,PASSQOFFSET,iPassP,iPassQ,sA,sB,sC,sD,n1,n2
integer :: iderivP,iderivQ,localA,localB,localC,localD,iContB,iContA,iAngA,iAngB
integer :: iB1,iTrans,iA1,ideriv,iContD,iContC,iAngC,iAngD,iD1,iC1
logical :: samePQ,SameRHSaos,SameLHSaos,SameODs,translate
!Integer :: Kac,Kbc,Kad,Kbd,Kca,Kcb,Kda,Kdb,Dbd,Dad,Dac,Dbc,Dca,Dcb,Dda,Ddb
Type(derivativeInfo) :: derivInfo
Integer     :: iAtom,iDer,DL_AC,DL_BC,DR_BD,DR_AD,DL_AD,DL_BD,DR_BC,DR_AC
Real(realk) :: derCont,Dab,kint,Dfac
Real(realk) :: DUMMY(1)

IF (Input%derivOrder.GT.0) THEN
  CALL initDerivativeOverlapInfo(derivInfo,PQ,Input)
  IF (IPRINT.GT.50) THEN
    CALL printDerivativeOverlapInfo(derivInfo,LUPRI)
  ENDIF
  IF (Input%NDMAT_LHS.NE.Input%NDMAT_RHS) CALL LSQUIT('Error in DistributeExchangeGrad. #LHS .NE. #RHS DMAT',lupri)
ENDIF

nMAT=Input%NDMAT_RHS
nAngmomP = PQ%P%p%nAngmom
nAngmomQ = PQ%Q%p%nAngmom
SamePQ = PQ%samePQ
SameLHSaos = INPUT%SameLHSaos
SameRHSaos = INPUT%SameRHSaos
SameODs = INPUT%SameODs
nDerivP = Input%nDerivP
nDerivQ = Input%nDerivQ

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
!
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
 DO iDerivP=1,nDerivP
  DO iPassP = 1, PQ%P%p%nPasses
   CALL getOverlapInfo(nContA,nContB,nAngA,nAngB,startA,startB,localA,localB,atomA,atomB,PQ%P%p,iAngmomP,iPassP)
   PASSPOFFSET = (iPassP-1)*nContA*nContB*nAngA*nAngB
   batchA = PQ%P%p%orbital1%batch(iPassP)
   batchB = PQ%P%p%orbital2%batch(iPassP)
   iOrbP=iOrbitalP-1+PASSPOFFSET
   derivInfo%Atom(1)=atomA
   derivInfo%Atom(2)=atomB
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
    DO iDerivQ = 1, nDerivQ
     DO iPassQ = 1, PQ%Q%p%nPasses
      CALL getOverlapInfo(nContC,nContD,nAngC,nAngD,startC,startD,localC,localD,atomC,atomD,PQ%Q%p,iAngmomQ,iPassQ)
      PASSQOFFSET = (iPassQ-1)*nContC*nContD*nAngD*nAngC
      iOrbQ=iOrbitalQ-1+PASSQOFFSET
      batchC = PQ%Q%p%orbital1%batch(iPassQ)
      batchD = PQ%Q%p%orbital2%batch(iPassQ)
      derivInfo%atom(3) = atomC
      derivInfo%atom(4) = atomD

      iDeriv=iderivP+(iderivQ-1)*NderivP
      iAtom = derivInfo%Atom(derivInfo%AO(1,iDeriv))
      iDer  = derivInfo%dirComp(iDeriv)
      translate = derivInfo%translate.GT.0
      iTrans = 0
      IF (translate) iTrans = derivInfo%Atom(derivInfo%translate)

      IF((startA.NE.startB).AND.SameLHSaos)THEN
       IF (((startA.NE.startC).OR.(startB.NE.startD)).AND.SameODs)THEN
        IF (INPUT%DO_GRADIENT) CALL LSQUIT('Implement grad-exchange 1!',lupri)
       ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN      
        DL_AC = Input%LST_DLHS%INDEX(ATOMA,ATOMC,1,1)
        DL_BC = Input%LST_DLHS%INDEX(ATOMB,ATOMC,1,1)
        DR_BD = Input%LST_DRHS%INDEX(ATOMB,ATOMD,1,1)
        DR_AD = Input%LST_DRHS%INDEX(ATOMA,ATOMD,1,1)
        IF ((startC.NE.startD).AND.SameRHSaos)THEN
           DL_AD = Input%LST_DLHS%INDEX(ATOMA,ATOMD,1,1)
           DL_BD = Input%LST_DLHS%INDEX(ATOMB,ATOMD,1,1)
           DR_BC = Input%LST_DRHS%INDEX(ATOMB,ATOMC,1,1)
           DR_AC = Input%LST_DRHS%INDEX(ATOMA,ATOMC,1,1)
           call lstensorExchangeGrad2(&
           & Input%LST_DLHS%LSAO(DL_AC)%BATCH(batchA,batchC,1,1)%elms,&
           & Input%LST_DLHS%LSAO(DL_AD)%BATCH(batchA,batchD,1,1)%elms,&
           & Input%LST_DLHS%LSAO(DL_BC)%BATCH(batchB,batchC,1,1)%elms,&
           & Input%LST_DLHS%LSAO(DL_BD)%BATCH(batchB,batchD,1,1)%elms,&
           & Input%LST_DRHS%LSAO(DR_BD)%BATCH(batchB,batchD,1,1)%elms,&
           & Input%LST_DRHS%LSAO(DR_BC)%BATCH(batchB,batchC,1,1)%elms,&
           & Input%LST_DRHS%LSAO(DR_AD)%BATCH(batchA,batchD,1,1)%elms,&
           & Input%LST_DRHS%LSAO(DR_AC)%BATCH(batchA,batchC,1,1)%elms,&
           & dercont,iOrbP,iOrbQ,nAngA,nAngB,nAngC,nAngD,nContA,nContB,&
           & nContC,nContD,QPmat2,dimQ,dimP,nmat,dimA,dimB,dimC,dimD,&
           & sA,sB,sC,sD,lupri,iprint,1)
        ELSE
           call lstensorExchangeGrad2(&
           & Input%LST_DLHS%LSAO(DL_AC)%BATCH(batchA,batchC,1,1)%elms,DUMMY,&
           & Input%LST_DLHS%LSAO(DL_BC)%BATCH(batchB,batchC,1,1)%elms,DUMMY,&
           & Input%LST_DRHS%LSAO(DR_BD)%BATCH(batchB,batchD,1,1)%elms,DUMMY,&
           & Input%LST_DRHS%LSAO(DR_AD)%BATCH(batchA,batchD,1,1)%elms,DUMMY,&
           & dercont,iOrbP,iOrbQ,nAngA,nAngB,nAngC,nAngD,nContA,nContB,&
           & nContC,nContD,QPmat2,dimQ,dimP,nmat,dimA,dimB,dimC,dimD,&
           & sA,sB,sC,sD,lupri,iprint,2)
        ENDIF
       ENDIF
      ELSE
       IF (((startA.NE.startC).OR.(startB.NE.startD)).AND.SameODs)THEN
        IF (INPUT%DO_GRADIENT) CALL LSQUIT('Implement grad-exchange 1!',lupri)
       ELSE !NOT ((startA.NE.startC).OR.(startB.NE.startD)) THEN
        DL_AC = Input%LST_DLHS%INDEX(ATOMA,ATOMC,1,1)
        DR_BD = Input%LST_DRHS%INDEX(ATOMB,ATOMD,1,1)
        IF ((startC.NE.startD).AND.SameRHSaos)THEN
           DL_AD = Input%LST_DLHS%INDEX(ATOMA,ATOMD,1,1)
           DR_BC = Input%LST_DRHS%INDEX(ATOMB,ATOMC,1,1)
           call lstensorExchangeGrad2(&
           & Input%LST_DLHS%LSAO(DL_AC)%BATCH(batchA,batchC,1,1)%elms,&
           & Input%LST_DLHS%LSAO(DL_AD)%BATCH(batchA,batchD,1,1)%elms,DUMMY,DUMMY,&
           & Input%LST_DRHS%LSAO(DR_BD)%BATCH(batchB,batchD,1,1)%elms,&
           & Input%LST_DRHS%LSAO(DR_BC)%BATCH(batchB,batchC,1,1)%elms,DUMMY,DUMMY,&
           & dercont,iOrbP,iOrbQ,nAngA,nAngB,nAngC,nAngD,nContA,nContB,&
           & nContC,nContD,QPmat2,dimQ,dimP,nmat,dimA,dimB,dimC,dimD,&
           & sA,sB,sC,sD,lupri,iprint,3)
        ELSE
           call lstensorExchangeGrad2(&
           & Input%LST_DLHS%LSAO(DL_AC)%BATCH(batchA,batchC,1,1)%elms,DUMMY,DUMMY,DUMMY,&
           & Input%LST_DRHS%LSAO(DR_BD)%BATCH(batchB,batchD,1,1)%elms,DUMMY,DUMMY,DUMMY,&
           & dercont,iOrbP,iOrbQ,nAngA,nAngB,nAngC,nAngD,nContA,nContB,&
           & nContC,nContD,QPmat2,dimQ,dimP,nmat,dimA,dimB,dimC,dimD,&
           & sA,sB,sC,sD,lupri,iprint,4)
        ENDIF
       ENDIF
      ENDIF
!$OMP CRITICAL (distributeEXCHANGEGRAD) 
      Kmat%LSAO(iatom)%BATCH(1,1,1,1)%elms(ider) = Kmat%LSAO(iatom)%BATCH(1,1,1,1)%elms(ider) - 0.5d0 * derCont
      IF (translate) THEN    
         Kmat%LSAO(iTrans)%BATCH(1,1,1,1)%elms(ider) = Kmat%LSAO(iTrans)%BATCH(1,1,1,1)%elms(ider) + 0.5d0 * derCont
      ENDIF
!$OMP END CRITICAL (distributeEXCHANGEGRAD) 
     ENDDO
     iOrbitalQ = iOrbitalQ + nOrbQ
    ENDDO
   ENDDO
  ENDDO
  iOrbitalP = iOrbitalP + nOrbP
 ENDDO
ENDDO
!write(lupri,*)'TK grad'
!call print_lstensor(RES,lupri)
IF (Input%derivOrder.GT.0) call freeDerivativeOverlapInfo(derivInfo)

END SUBROUTINE LSTENSOREXCHANGEGRAD

!> \brief distribute Integral to exchange gradient
!> \author S. Reine
!> \date 2010-04-23
!> \param DLHSac left hand side Density matrix \latexonly D_{ac} \endlatexonly
!> \param DLHSad left hand side Density matrix \latexonly D_{ad} \endlatexonly
!> \param DLHSbc left hand side Density matrix \latexonly D_{bc} \endlatexonly
!> \param DLHSbd left hand side Density matrix \latexonly D_{bd} \endlatexonly
!> \param DRHSbd right hand side Density matrix \latexonly D_{bd} \endlatexonly
!> \param DRHSbc right hand side Density matrix \latexonly D_{bc} \endlatexonly
!> \param DRHSad right hand side Density matrix \latexonly D_{ad} \endlatexonly
!> \param DRHSac right hand side Density matrix \latexonly D_{ac} \endlatexonly
!> \param Dercont derivate contribution to gradient
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nContA number of contracted functions on the A batch
!> \param nContB number of contracted functions on the B batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nContD number of contracted functions on the D batch 
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
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
!> \param option flag indicating the permutational symmetry 1 for full permutational sym, 2 not same AOs on RHS, 3 not same AOs on LHS, 4 no permut sym.  
subroutine lstensorExchangeGrad2(DLHSac,DLHSad,DLHSbc,DLHSbd,DRHSbd,DRHSbc,DRHSad,DRHSac,Dercont,start_iOrbP,start_iOrbQ,&
     & nAngA,nAngB,nAngC,nAngD,nContA,nContB,nContC,nContD,QPmat,dimQ,dimP,nmat,dimA,dimB,dimC,dimD,sA,sB,sC,sD,&
     & lupri,iprint,option)
implicit none
integer,intent(in)        :: nmat,lupri,iprint,dimQ,dimP,dimA,dimB,dimC,dimD
real(realk) :: DLHSac(dimA,dimC,nmat),DLHSad(dimA,dimD,nmat),DLHSbc(dimB,dimC,nmat),DLHSbd(dimB,dimD,nmat)
real(realk) :: DRHSbd(dimB,dimD,nmat),DRHSbc(dimB,dimC,nmat),DRHSad(dimA,dimD,nmat),DRHSac(dimA,dimC,nmat)
integer,intent(inout)     :: start_iOrbP,start_iOrbQ
integer,intent(in)        :: nAngA,nAngB,nAngC,nAngD
integer,intent(in)        :: nContA,nContB,nContC,nContD,sA,sB,sC,sD,option
real(realk),intent(in)    :: QPmat(dimQ,dimP)
real(realk),intent(inout) :: Dercont
!
integer                   :: iContB,iContA,iContC,iContD,iOrbP,iOrbQ
integer                   :: iAngB,iAngA,iAngC,iAngD,iB,iA,iC,iD,idmat
real(realk)               :: kint,Dfac

Dercont=0
iorbP = start_iorbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB=sB+iAngB+(iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA=sA+iAngA+(iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD=sD+iAngD+(iContD-1)*nAngD
       DO iAngC=1,nAngC
        iC=sC+iAngC+(iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint  = QPmat(iOrbQ,iOrbP)
        Dfac = 0.d0
        IF(option.EQ.1)THEN
           DO idmat=1,nmat
            Dfac = Dfac + DLHSac(iA,iC,idmat) * DRHSbd(iB,iD,idmat)
            Dfac = Dfac + DLHSad(iA,iD,idmat) * DRHSbc(iB,iC,idmat)
            Dfac = Dfac + DLHSbc(iB,iC,idmat) * DRHSad(iA,iD,idmat)
            Dfac = Dfac + DLHSbd(iB,iD,idmat) * DRHSac(iA,iC,idmat)
           ENDDO
        ELSEIF(option.EQ.2)THEN
           DO idmat=1,nmat
            Dfac = Dfac + DLHSac(iA,iC,idmat) * DRHSbd(iB,iD,idmat)
            Dfac = Dfac + DLHSbc(iB,iC,idmat) * DRHSad(iA,iD,idmat)
           ENDDO
        ELSEIF(option.EQ.3)THEN
           DO idmat=1,nmat
            Dfac = Dfac + DLHSac(iA,iC,idmat) * DRHSbd(iB,iD,idmat)
            Dfac = Dfac + DLHSad(iA,iD,idmat) * DRHSbc(iB,iC,idmat)
           ENDDO
        ELSE
           DO idmat=1,nmat
            Dfac = Dfac + DLHSac(iA,iC,idmat) * DRHSbd(iB,iD,idmat)
           ENDDO
        ENDIF
        derCont = derCont + kint*Dfac
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE lstensorExchangeGrad2

END MODULE thermite_distribute_Exchange
