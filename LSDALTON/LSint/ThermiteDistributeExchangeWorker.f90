!> @file 
!> Contains integral distribution routines that places calculated integrals in the proper output
!> Thermite integral distribution module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE thermite_distribute_ExchangeW
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
  use thermite_distribute, only: derivativeInfo, initderivativeoverlapinfo, &
       &printderivativeoverlapinfo, getoverlapinfo, freederivativeoverlapinfo

SAVE
!Integer     :: nCalcInt
!Integer     :: nCalcIntZero
!Integer     :: nCalcIntZeroContrib
!real(realk) :: intthreshold
CONTAINS
!> \brief distribute exchange contribution all permutational symmetry
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat1 matrix containing calculated integrals
!> \param nmat number of density matrices
!> \param AeqC is startA.EQ.startC
!> \param AeqC is startB.EQ.startD
SUBROUTINE exchangePermutTTT(Kbd,Dac,Kad,Dbc,Kbc,Dad,Kac,Dbd,Kdb,Dca,Kcb,Dda,Kca,Ddb,Kda,Dcb,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat1,nMat,AeqC,BeqD)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat1(:,:)
REAL(REALK),intent(in)    :: Dbd(dimB,dimD,nmat),Dad(dimA,dimD,nmat),Dbc(dimB,dimC,nmat),Dac(dimA,dimC,nmat)
REAL(REALK),intent(in)    :: Ddb(dimD,dimB,nmat),Dda(dimD,dimA,nmat),Dcb(dimC,dimB,nmat),Dca(dimC,dimA,nmat)
REAL(REALK),intent(inout) :: Kac(dimA,dimC,nmat),Kbc(dimB,dimC,nmat),Kad(dimA,dimD,nmat),Kbd(dimB,dimD,nmat)
REAL(REALK),intent(inout) :: Kca(dimC,dimA,nmat),Kcb(dimC,dimB,nmat),Kda(dimD,dimA,nmat),Kdb(dimD,dimB,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
LOGICAL        :: zeroInt,zeroCont,AeqC,BeqD
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       do imat = 1,nmat
        tDbd(imat) = Dbd(iB,iD,imat)
        tDad(imat) = Dad(iA,iD,imat)
        tDdb(imat) = Ddb(iD,iB,imat)
        tDda(imat) = Dda(iD,iA,imat)
       enddo
       sum_ad = 0
       sum_bd = 0
       sum_da = 0
       sum_db = 0
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = QPmat1(iOrbQ,iOrbP)
        !nCalcInt = nCalcInt+1
        !IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          Kac(iA,iC,imat) = Kac(iA,iC,imat)+kint*tDbd(imat)
          Kbc(iB,iC,imat) = Kbc(iB,iC,imat)+kint*tDad(imat)
          IF(AeqC)THEN
             Kac(iC,iA,imat) = Kac(iC,iA,imat)+kint*tDdb(imat)
          ELSE
             Kca(iC,iA,imat) = Kca(iC,iA,imat)+kint*tDdb(imat)
          ENDIF
          Kcb(iC,iB,imat) = Kcb(iC,iB,imat)+kint*tDda(imat)           
          sum_ad(imat) = sum_ad(imat) + kint*Dbc(iB,iC,imat)
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          sum_da(imat) = sum_da(imat) + kint*Dcb(iC,iB,imat)
          sum_db(imat) = sum_db(imat) + kint*Dca(iC,iA,imat)
          !--------------------------------------------------------
         ! IF(ABS(kint*tDbd(imat)).GT.intthreshold)zeroCont = .FALSE.
         ! IF(ABS(kint*tDad(imat)).GT.intthreshold)zeroCont = .FALSE.
         ! IF(ABS(kint*tDdb(imat)).GT.intthreshold)zeroCont = .FALSE.
         ! IF(ABS(kint*tDda(imat)).GT.intthreshold)zeroCont = .FALSE.
         ! IF(ABS(kint*Dbc(iB,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
         ! IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
         ! IF(ABS(kint*Dcb(iC,iB,imat)).GT.intthreshold)zeroCont = .FALSE.
         ! IF(ABS(kint*Dca(iC,iA,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kad(iA,iD,imat) = Kad(iA,iD,imat)+sum_ad(imat)
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
          Kda(iD,iA,imat) = Kda(iD,iA,imat)+sum_da(imat)
          IF(BeqD)THEN
             Kbd(iD,iB,imat) = Kbd(iD,iB,imat)+sum_db(imat)
          ELSE
             Kdb(iD,iB,imat) = Kdb(iD,iB,imat)+sum_db(imat)
          ENDIF
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutTTT

!> \brief distribute exchange contribution all permutational symmetry, assuming symmetric density matrix
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat2 matrix containing calculated integrals
!> \param nmat number of density matrices
SUBROUTINE exchangePermutTTTsym(Kbd,Dac,Kad,Dbc,Kbc,Dad,Kac,Dbd,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat2,nMat)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat2(:,:)
REAL(REALK),intent(in)    :: Dbd(dimB,dimD,nmat),Dad(dimA,dimD,nmat),Dbc(dimB,dimC,nmat),Dac(dimA,dimC,nmat)
REAL(REALK),intent(inout) :: Kac(dimA,dimC,nmat),Kbc(dimB,dimC,nmat),Kad(dimA,dimD,nmat),Kbd(dimB,dimD,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
real(realk),parameter :: D2=2.d0
!LOGICAL        :: zeroInt,zeroCont
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       do imat = 1,nmat
        tDbd(imat) = Dbd(iB,iD,imat)
        tDad(imat) = Dad(iA,iD,imat)
       enddo
       sum_ad = 0
       sum_bd = 0
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = D2*QPmat2(iOrbQ,iOrbP)
!        nCalcInt = nCalcInt+1
!        IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          Kac(iA,iC,imat) = Kac(iA,iC,imat)+kint*tDbd(imat)
          Kbc(iB,iC,imat) = Kbc(iB,iC,imat)+kint*tDad(imat)
          sum_ad(imat) = sum_ad(imat) + kint*Dbc(iB,iC,imat)
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*tDbd(imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*tDad(imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dbc(iB,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kad(iA,iD,imat) = Kad(iA,iD,imat)+sum_ad(imat)
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutTTTsym

!> \brief distribute exchange contribution TFT (AB sym, CD sym, ODsym) no CD sym
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat3 matrix containing calculated integrals
!> \param nmat number of density matrices
SUBROUTINE exchangePermutTFT(Kbd,Dac,Kad,Dbc,Kdb,Dca,Kda,Dcb,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat3,nMat)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat3(:,:)
REAL(REALK),intent(in)    :: Dbc(dimB,dimC,nmat),Dac(dimA,dimC,nmat)
REAL(REALK),intent(in)    :: Dcb(dimC,dimB,nmat),Dca(dimC,dimA,nmat)
REAL(REALK),intent(inout) :: Kad(dimA,dimD,nmat),Kbd(dimB,dimD,nmat)
REAL(REALK),intent(inout) :: Kda(dimD,dimA,nmat),Kdb(dimD,dimB,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
!LOGICAL        :: zeroInt,zeroCont
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = QPmat3(iOrbQ,iOrbP)
 !       nCalcInt = nCalcInt+1
 !       IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          sum_ad(imat) = sum_ad(imat) + kint*Dbc(iB,iC,imat)
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          sum_da(imat) = sum_da(imat) + kint*Dcb(iC,iB,imat)
          sum_db(imat) = sum_db(imat) + kint*Dca(iC,iA,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*Dbc(iB,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dcb(iC,iB,imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dca(iC,iA,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kad(iA,iD,imat) = Kad(iA,iD,imat)+sum_ad(imat)
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
          Kda(iD,iA,imat) = Kda(iD,iA,imat)+sum_da(imat)
          Kdb(iD,iB,imat) = Kdb(iD,iB,imat)+sum_db(imat)
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutTFT

!> \brief distribute exchange contribution TFT (AB sym, CD sym, ODsym) assuming symmetric density matrix (no CD sym)
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat4 matrix containing calculated integrals
!> \param nmat number of density matrices
SUBROUTINE exchangePermutTFTsym(Kbd,Dac,Kad,Dbc,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat4,nMat)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat4(:,:)
REAL(REALK),intent(in)    :: Dbc(dimB,dimC,nmat),Dac(dimA,dimC,nmat)
REAL(REALK),intent(inout) :: Kad(dimA,dimD,nmat),Kbd(dimB,dimD,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
real(realk),parameter :: D2=2.d0
!LOGICAL        :: zeroInt,zeroCont
!zeroInt = .TRUE.
!zeroCont = .TRUE.
iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       sum_ad = 0
       sum_bd = 0
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = D2*QPmat4(iOrbQ,iOrbP)
!        nCalcInt = nCalcInt+1
!        IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          sum_ad(imat) = sum_ad(imat) + kint*Dbc(iB,iC,imat)
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*Dbc(iB,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kad(iA,iD,imat) = Kad(iA,iD,imat)+sum_ad(imat)
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutTFTsym

!> \brief distribute exchange contribution TTF (AB sym, CD sym, ODsym) no ODsym
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat5 matrix containing calculated integrals
!> \param nmat number of density matrices
SUBROUTINE exchangePermutTTF(Kbd,Dac,Kad,Dbc,Kbc,Dad,Kac,Dbd,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat5,nMat)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat5(:,:)
REAL(REALK),intent(in)    :: Dbd(dimB,dimD,nmat),Dad(dimA,dimD,nmat),Dbc(dimB,dimC,nmat),Dac(dimA,dimC,nmat)
REAL(REALK),intent(inout) :: Kac(dimA,dimC,nmat),Kbc(dimB,dimC,nmat),Kad(dimA,dimD,nmat),Kbd(dimB,dimD,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
!LOGICAL        :: zeroInt,zeroCont
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       do imat = 1,nmat
        tDbd(imat) = Dbd(iB,iD,imat)
        tDad(imat) = Dad(iA,iD,imat)
       enddo
       sum_ad = 0
       sum_bd = 0
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = QPmat5(iOrbQ,iOrbP)
!        nCalcInt = nCalcInt+1
!        IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          Kac(iA,iC,imat) = Kac(iA,iC,imat)+kint*tDbd(imat)
          Kbc(iB,iC,imat) = Kbc(iB,iC,imat)+kint*tDad(imat)
          sum_ad(imat) = sum_ad(imat) + kint*Dbc(iB,iC,imat)
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*tDbd(imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*tDad(imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dbc(iB,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kad(iA,iD,imat) = Kad(iA,iD,imat)+sum_ad(imat)
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutTTF

!> \brief distribute exchange contribution TFF (AB sym, CD sym, ODsym) only AB sym 
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat6 matrix containing calculated integrals
!> \param nmat number of density matrices
SUBROUTINE exchangePermutTFF(Kbd,Dac,Kad,Dbc,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat6,nMat)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat6(:,:)
REAL(REALK),intent(in)    :: Dac(dimA,dimC,nmat),Dbc(dimB,dimC,nmat)
REAL(REALK),intent(inout) :: Kbd(dimB,dimD,nmat),Kad(dimA,dimD,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
!LOGICAL        :: zeroInt,zeroCont
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       sum_ad = 0
       sum_bd = 0
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = QPmat6(iOrbQ,iOrbP)
!        nCalcInt = nCalcInt+1
!        IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          sum_ad(imat) = sum_ad(imat) + kint*Dbc(iB,iC,imat)
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*Dbc(iB,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kad(iA,iD,imat) = Kad(iA,iD,imat)+sum_ad(imat)
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutTFF

!> \brief distribute exchange contribution FTT (AB sym, CD sym, ODsym) no AB sym 
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat7 matrix containing calculated integrals
!> \param nmat number of density matrices
!> \param BeqC Bstart equal to Cstart?
!> \param BeqD Bstart equal to Dstart?
SUBROUTINE exchangePermutFTT(Kbd,Dac,Kbc,Dad,Kdb,Dca,Kcb,Dda,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat7,nMat,BeqC,BeqD)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat7(:,:)
REAL(REALK),intent(in)    :: Dac(dimA,dimC,nmat),Dad(dimA,dimD,nmat),Dca(dimC,dimA,nmat),Dda(dimD,dimA,nmat)
REAL(REALK),intent(inout) :: Kbd(dimB,dimD,nmat),Kbc(dimB,dimC,nmat),Kdb(dimD,dimB,nmat),Kcb(dimC,dimB,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
LOGICAL        :: zeroInt,zeroCont,BeqC,BeqD
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       do imat = 1,nmat
        tDad(imat) = Dad(iA,iD,imat)
        tDda(imat) = Dda(iD,iA,imat)
       enddo
       sum_bd = 0
       sum_db = 0
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = QPmat7(iOrbQ,iOrbP)
!        nCalcInt = nCalcInt+1
!        IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          Kbc(iB,iC,imat) = Kbc(iB,iC,imat)+kint*tDad(imat)
          IF(BeqC)THEN
             Kbc(iC,iB,imat) = Kbc(iC,iB,imat)+kint*tDda(imat)
          ELSE
             Kcb(iC,iB,imat) = Kcb(iC,iB,imat)+kint*tDda(imat)
          ENDIF
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          sum_db(imat) = sum_db(imat) + kint*Dca(iC,iA,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*tDad(imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*tDda(imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dca(iC,iA,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
          IF(BeqD)THEN
             Kbd(iB,iD,imat) = Kdb(iB,iD,imat)+sum_db(imat)
          ELSE
             Kdb(iD,iB,imat) = Kdb(iD,iB,imat)+sum_db(imat)
          ENDIF
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutFTT

!> \brief distribute exchange contribution FTT (AB sym, CD sym, ODsym) no AB sym, assume Sym D mat 
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat8 matrix containing calculated integrals
!> \param nmat number of density matrices
SUBROUTINE exchangePermutFTTsym(Kbd,Dac,Kbc,Dad,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat8,nMat)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat8(:,:)
REAL(REALK),intent(in)    :: Dac(dimA,dimC,nmat),Dad(dimA,dimD,nmat)
REAL(REALK),intent(inout) :: Kbd(dimB,dimD,nmat),Kbc(dimB,dimC,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
real(realk),parameter :: D2=2.d0
!LOGICAL        :: zeroInt,zeroCont
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       do imat = 1,nmat
        tDad(imat) = Dad(iA,iD,imat)
       enddo
       sum_ad = 0
       sum_bd = 0
       sum_da = 0
       sum_db = 0
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = D2*QPmat8(iOrbQ,iOrbP)
!        nCalcInt = nCalcInt+1
!        IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          Kbc(iB,iC,imat) = Kbc(iB,iC,imat)+kint*tDad(imat)
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*tDad(imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutFTTsym

!> \brief distribute exchange contribution FFT (AB sym, CD sym, ODsym) only OD sym
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat9 matrix containing calculated integrals
!> \param nmat number of density matrices
SUBROUTINE exchangePermutFFT(Kbd,Dac,Kdb,Dca,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat9,nMat)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat9(:,:)
REAL(REALK),intent(inout) :: Kbd(dimB,dimD,nmat),Kdb(dimD,dimB,nmat)
REAL(REALK),intent(in)    :: Dac(dimA,dimC,nmat),Dca(dimC,dimA,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
!LOGICAL        :: zeroInt,zeroCont
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       sum_bd = 0
       sum_db = 0
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = QPmat9(iOrbQ,iOrbP)
!        nCalcInt = nCalcInt+1
!        IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          sum_db(imat) = sum_db(imat) + kint*Dca(iC,iA,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dca(iC,iA,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
          Kdb(iD,iB,imat) = Kdb(iD,iB,imat)+sum_db(imat)
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutFFT

!> \brief distribute exchange contribution FTF (AB sym, CD sym, ODsym) only CD sym
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat10 matrix containing calculated integrals
!> \param nmat number of density matrices
SUBROUTINE exchangePermutFTF(Kbd,Dac,Kbc,Dad,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,&
     & nContD,nContC,nAngD,nAngC,QPmat10,nMat)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat10(:,:)
REAL(REALK),intent(in)    :: Dac(dimA,dimC,nmat),Dad(dimA,dimD,nmat)
REAL(REALK),intent(inout) :: Kbd(dimB,dimD,nmat),KBC(dimB,dimC,nmat)
!
REAL(REALK)    :: kint,tDbd(nmat),tDad(nmat),tDdb(nmat),tDda(nmat)
REAL(REALK)    :: sum_ad(nmat),sum_bd(nmat),sum_da(nmat),sum_db(nmat)
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
!LOGICAL        :: zeroInt,zeroCont
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       do imat = 1,nmat
        tDad(imat) = Dad(iA,iD,imat)
       enddo
       sum_bd = 0
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = QPmat10(iOrbQ,iOrbP)
!        nCalcInt = nCalcInt+1
!        IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
          Kbc(iB,iC,imat) = Kbc(iB,iC,imat)+kint*tDad(imat)
          sum_bd(imat) = sum_bd(imat) + kint*Dac(iA,iC,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*tDad(imat)).GT.intthreshold)zeroCont = .FALSE.
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       enddo
       do imat = 1,nmat
          Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+sum_bd(imat)
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutFTF

!> \brief distribute exchange contribution no permutational sym
!> \author T. Kjaergaard
!> \date 2010-04-23
!> \param Kbd the exchange matrix \latexonly K_{bd} = (ab|cd)D_{ac}\endlatexonly
!> \param Dac the exchange matrix \latexonly D_{ac}\endlatexonly
!> \param Kad the exchange matrix \latexonly K_{ad} = (ba|cd)D_{bc}\endlatexonly
!> \param Dbc the exchange matrix \latexonly D_{bc}\endlatexonly
!> \param Kbc the exchange matrix \latexonly K_{bc} = (ab|dc)D_{ad}\endlatexonly
!> \param Dad the exchange matrix \latexonly D_{ad}\endlatexonly
!> \param Kac the exchange matrix \latexonly K_{ac} = (ba|dc)D_{bd}\endlatexonly
!> \param Dbd the exchange matrix \latexonly D_{bd}\endlatexonly
!> \param Kbd the exchange matrix \latexonly K_{db} = (ad|cb)D_{ca}\endlatexonly
!> \param Dca the exchange matrix \latexonly D_{ca}\endlatexonly
!> \param Kcb the exchange matrix \latexonly K_{cb} = (ad|cb)D_{da}\endlatexonly
!> \param Dda the exchange matrix \latexonly D_{da}\endlatexonly
!> \param Kca the exchange matrix \latexonly K_{ca} = (ba|dc)D_{db}\endlatexonly
!> \param Ddb the exchange matrix \latexonly D_{db}\endlatexonly
!> \param Kda the exchange matrix \latexonly K_{da} = (ba|dc)D_{cb}\endlatexonly
!> \param Dcb the exchange matrix \latexonly D_{cb}\endlatexonly
!> \param dimA the dimension of the A batch
!> \param dimB the dimension of the B batch
!> \param dimC the dimension of the C batch
!> \param dimD the dimension of the D batch
!> \param sA the starting orbital index of the A batch  
!> \param sB the starting orbital index of the B batch  
!> \param sC the starting orbital index of the C batch  
!> \param sD the starting orbital index of the D batch  
!> \param start_iOrbP starting orbital index for the P overlapdistribution
!> \param start_iOrbQ starting orbital index for the Q overlapdistribution
!> \param nContB number of contracted functions on the B batch 
!> \param nContA number of contracted functions on the A batch
!> \param nAngB number of angular momentums for the B center/B batch
!> \param nAngA number of angular momentums for the A center/A batch (1. center)
!> \param nContD number of contracted functions on the D batch 
!> \param nContC number of contracted functions on the C batch 
!> \param nAngD number of angular momentums for the D center/D batch
!> \param nAngC number of angular momentums for the C center/C batch
!> \param QPmat11 matrix containing calculated integrals
!> \param nmat number of density matrices
SUBROUTINE exchangePermutFFF(Kbd,Dac,&
     & dimA,dimB,dimC,dimD,sA,sB,sC,sD,&
     &START_iOrbP,START_iOrbQ,nContB,nContA,nAngB,nAngA,nContD,nContC,&
     &nAngD,nAngC,QPmat11,nmat)
IMPLICIT NONE
INTEGER,intent(in)        :: dimA,dimB,dimC,dimD,sA,sB,sC,sD,nmat
INTEGER,intent(in)        :: START_iOrbP,START_iOrbQ,nContB,nContA
INTEGER,intent(in)        :: nAngB,nAngA,nContD,nContC,nAngD,nAngC
REAL(REALK),intent(in)    :: QPmat11(:,:)
REAL(REALK),intent(in)    :: Dac(dimA,dimC,nmat)
REAL(REALK),intent(inout) :: Kbd(dimB,dimD,nmat)
!
REAL(REALK)    :: kint,tDbd,tDad,tDdb,tDda,sum_ad,sum_bd,sum_da,sum_db
INTEGER        :: iAngA,iAngB,iAngC,iAngD,iContA,iContB,iContC,iContD
INTEGER        :: iOrbP,iOrbQ,iA,iB,iC,iD,imat
!LOGICAL        :: zeroInt,zeroCont
!zeroInt = .TRUE.
!zeroCont = .TRUE.

iOrbP = START_iOrbP
DO iContB=1,nContB
 DO iContA=1,nContA
  DO iAngB=1,nAngB
   iB = sB + iAngB + (iContB-1)*nAngB
   DO iAngA=1,nAngA
    iA = sA + iAngA + (iContA-1)*nAngA
    iOrbP=iOrbP+1
    iOrbQ = START_iOrbQ
    DO iContD=1,nContD
     DO iContC=1,nContC
      DO iAngD=1,nAngD
       iD = sD + iAngD + (iContD-1)*nAngD
       DO iAngC=1,nAngC
        iC = sC + iAngC + (iContC-1)*nAngC
        iOrbQ=iOrbQ+1
        kint = QPmat11(iOrbQ,iOrbP)
!        nCalcInt = nCalcInt+1
!        IF(ABS(kint) .GT.intthreshold) zeroInt = .FALSE.
        do imat = 1,nmat
           Kbd(iB,iD,imat) = Kbd(iB,iD,imat)+kint*Dac(iA,iC,imat)
          !--------------------------------------------------------
!          IF(ABS(kint*Dac(iA,iC,imat)).GT.intthreshold)zeroCont = .FALSE.
          !--------------------------------------------------------
        enddo
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
!IF(zeroInt) nCalcIntZero = nCalcIntZero+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC
!IF(zeroCont) nCalcIntZeroContrib = nCalcIntZeroContrib+nContB*nAngB*nContA*nAngA*nContD*nAngD*nContC*nAngC

END SUBROUTINE exchangePermutFFF
!SUBROUTINE init_IntCalcCounter()
!implicit none
!
!nCalcInt = 0
!nCalcIntZero = 0
!nCalcIntZeroContrib = 0
!
!END SUBROUTINE init_IntCalcCounter

!SUBROUTINE get_IntCalcCounter(nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO)
!implicit none
!integer :: nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO
!
!nCalcIntO = nCalcInt 
!nCalcIntZeroO = nCalcIntZero
!nCalcIntZeroContribO = nCalcIntZeroContrib
!
!END SUBROUTINE get_IntCalcCounter

END MODULE thermite_distribute_ExchangeW
