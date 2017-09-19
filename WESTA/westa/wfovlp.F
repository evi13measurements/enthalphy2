#ifdef UNDEF
/* Comdeck todo */
(1) OCEINT ???
(2) jan 2005: check if zero vectors in SMO in BIORAS, seen to
    occur with bad active space dividing deg. orbitals, where
    sometimes pi_x, sometimes pi_y was included. /hjaaj
#endif
#  /* Deck wfovlp */
      FUNCTION WFOVLP(CVECA,CMOA,CVECB,CMOB,SAO,INDXCI,WRK,LWRK)
C
C     WESTA.WFOLP 17-May-1986 hjaaj
C                 29-Apr-1989 ha. RAS version
C
C     (calculate overlap between two RAS wave functions)
C
C Assumptions:
C   NCONF  in /INFVAR/ is correct
C
#include "implicit.h"
      DIMENSION CVECA(*),CMOA(*),CVECB(*),CMOB(*),
     *          SAO(*),WRK(*), INDXCI(*)
C
#include "iratdef.h"
C
      PARAMETER (D1  =  1.0D0)
C
C  PRIUNIT : LUPRI
C  INFORB : NORBT, NBAST
C  INFVAR : NCONF
C  WESINF : IPRDEF
C
#include "priunit.h"
#include "maxorb.h"
#include "inforb.h"
#include "infvar.h"
#include "wesinf.h"
C
C     CORE ALLOCATION
C     ..BIORAS..
C
C...BIORAS
C
      KVCA1  = 1
      KVCA2  = KVCA1  + NORBT*NBAST
      KTRA1A = KVCA2  + NORBT*NBAST
      KTRA2A = KTRA1A + NASHT * NASHT
      KTIN1  = KTRA2A + NASHT * NASHT
      KTIN2  = KTIN1  + NISHT
      KW1    = KTIN2  + NISHT
      LW1    = LWRK   - KW1
C
      IF (KW1 .GT. LWRK) CALL ERRWRK('WFOVLP.BIORAS',KW1,LWRK)
C
C
CDBG  iprbio = 99
      IF (IPRBIO.GE.5) THEN
         CALL HEADER('Biorthogonal orbitals from SYMORT_WFOVLP',1)
      END IF
C
C  We take only active part of TRA-matrices out of BIORAS +
C  the inactive diagonal parts for scaling.
C
      CALL DZERO(WRK(KTRA1A),2*NASHT*NASHT)
      CALL BIORAS(CMOA,CMOB,WRK(KVCA1),WRK(KVCA2),
     *            SAO,WRK(KTIN1),WRK(KTIN2),WRK(KTRA1A),
     *            WRK(KTRA2A),WRK(KW1),LW1)
C     CALL BIORAS(CMO1,CMO2,VC1,VC2,SAO12,TIN1,TIN2,TRA1A,
C    *            TRA2A,WRK,LWRK)
C
CDBG  write (lupri,*) 'TRA1A matrix'
CDBG  call output(wrk(ktra1a),1,nasht,1,nasht,nasht,nasht,1,lupri)
CDBG  write (lupri,*) 'TRA2A matrix'
CDBG  call output(wrk(ktra2a),1,nasht,1,nasht,nasht,nasht,1,lupri)
CDBG  oa1b1 = ddot(nconf,cveca,1,cvecb,1)
CDBG  write (lupri,*) 'CVECA, CVECB before WESTRACI', oa1b1
CDBG  sum = 0.0d0
CDBG  do i = 1,nconf
CDBG     sum = sum + cveca(i)*cvecb(i)
CDBG     if (abs(cveca(i)) + abs(cvecb(i)) .gt. 1.0D-2) then
CDBG        WRITE (lupri,'(I8,3F20.14)') I, CVECA(I), CVECB(I), SUM
CDBG     end if
CDBG  end do

C
C     CALL WESTRACI(NCVEC,CVEC,LCVEC,UMO,XNDXCI,WRK,LFREE)
C
      IF ( NCONF .GT. 1 ) THEN
         CALL WESTRACI(1,CVECA,NCONF,WRK(KTRA1A),INDXCI,WRK(KW1),LW1)
         CALL WESTRACI(1,CVECB,NCONF,WRK(KTRA2A),INDXCI,WRK(KW1),LW1)
CDBG     cnrma = dnrm2(nconf,cveca,1)
CDBG     cnrmb = dnrm2(nconf,cvecb,1)
CDBG     oa1b1 = ddot(nconf,cveca,1,cvecb,1)
CDBG     write (lupri,*) 'CVECA, CVECB after WESTRACI',cnrma,cnrmb,oa1b1
CDBG  sum = 0.0d0
CDBG  do i = 1,nconf
CDBG     sum = sum + cveca(i)*cvecb(i)
CDBG     if (abs(cveca(i)) + abs(cvecb(i)) .gt. 1.0D-2) then
CDBG        WRITE (lupri,'(I8,3F20.14)') I, CVECA(I), CVECB(I), sum
CDBG     end if
CDBG  end do
      END IF
C
C Scale CI vectors by (Tin1)**2 resp. (Tin2)**2
C
      FAC1 = D1
      FAC2 = D1
      DO I = 1, NISHT
        FAC1 = FAC1 * WRK(KTIN1+I-1) * WRK(KTIN1+I-1)
        FAC2 = FAC2 * WRK(KTIN2+I-1) * WRK(KTIN2+I-1)
      END DO

      IF (NCONF .EQ. 1 .AND. NASHT .GT. 0) THEN
      ! treat HSRHF cases (high spin RHF)
         DO I = 1, NASHT
            FAC1 = FAC1 * WRK(KTRA1A + (NASHT+1)*(I-1))
            FAC2 = FAC2 * WRK(KTRA2A + (NASHT+1)*(I-1))
         END DO
      END IF
C
      CALL DSCAL(NCONF,FAC1,CVECA,1)
      CALL DSCAL(NCONF,FAC2,CVECB,1)
C
C
      WFOVLP = DDOT(NCONF,CVECA,1,CVECB,1)
 
      IF (IPRBIO.GE.10) THEN
         WRITE (LUPRI,'(/A)') 'WFOVLP: inactive factors vector 1 and 2'
         CALL OUTPUT (WRK(KTIN1),1,NISHT,1,2,NISHT,2,1,LUPRI)
      END IF
      IF (IPRBIO.GE.5) WRITE(LUPRI,'(/A,3F20.8)')
     &   'WFOVLP: TIN1, TIN2, WFOVLP:', FAC1, FAC2, WFOVLP

C
      IF (IPRBIO .GE. 50) THEN
         WRITE (LUPRI,'(//A,I8)')
     *      'CI vector no. 1 and 2 at end SYMORT-WFOVLP, length =',NCONF
         sum = 0.0d0
         do i = 1,nconf
            sum = sum + cveca(i)*cvecb(i)
            if (abs(cveca(i)) + abs(cvecb(i)) .gt. 1.0D-2) then
               WRITE (lupri,'(I8,3F20.14)') I, CVECA(I), CVECB(I), sum
            end if
         end do
      END IF
C
      RETURN
C     END OF WESTA.WFOVLP
      END
#  /* Deck wesswf */
      SUBROUTINE WESSWF(SAO,SWF,TSWF,WRK,LWRK)
C
C 20-May-1986
C 30-Apr-1989 RAS-version (ha)
C
C Purpose:
C      To compute a matrix (SWF) containing all overlap elements
C      and the corrresponding transformation matrix
C      (= (SWF)-1/2) for a spectrum containing NSTAT2 wave functions.
C
C Input:
C      NSTAT2: number of states (wave functions) (common)
C      LUSPEC: unit numbers for wave functions in spectrum (common)
C      SAO   : atomic overlap matrix
C
C
C      Common: NCONF (NCONF = NCONFX),  LUSYMV(2)
C
C Output:
C      SWF   : overlap matrix between all states
C      TSWF  : transformation matrix corresponding to SWF
C
C Scratch:
C      WRK   : CI and MO files for 2 states
C
#include "implicit.h"
      DIMENSION SAO(*),SWF(*),TSWF(NSTAT2,*), WRK(LWRK)
#include "iratdef.h"
C
      PARAMETER ( D1 = 1.0D0 , THRESH = 1.D-4 )
C
C  PRIUNIT: LUPRI
C  WESINF: NSTAT2,NCONF2, IPRDEF
C  WESTAP: LUSPCW,LUSYMV(2)
C  INFORB: N2ASHX,NCMOT
C  INFVAR: NCONF
C  INFDIM: LCINDX
C  INFDIM: LSYM
C
#include "priunit.h"
#include "maxorb.h"
#include "wesinf.h"
#include "westap.h"
#include "inforb.h"
#include "infvar.h"
#include "infdim.h"
#include "infinp.h"
C
C
      IPRECV = IPRDEF-10
C
C Obtain NCONFX and LCINDX information
C we have identical symmetries..
C
      LSYM  = LSYM2
      NOSYM = 0
      CALL SETCI(NCONFX,NCDETS,LSYM,WRK,LWRK,NOSYM)
      NCONF = NCONFX
C
C ***** core allocation *****
C
      NNSTA2 = NSTAT2 * ( NSTAT2 + 1 ) / 2
      N2STA2 = NSTAT2 * NSTAT2
C
C CI vectors, orbitals..
      KCI1EX = 1
      KCI2EX = KCI1EX + NCONFX
      KCMO1  = KCI2EX + NCONFX
      KCMO2  = KCMO1  + NCMOT
      KCVEC1 = KCMO2  + NCMOT
      KCVEC2 = KCVEC1 + NCONF2
      KCINDX = KCVEC2 + NCONF2
      KW1    = KCINDX + LCINDX
      LW1    = LWRK   - KW1
      IF (KW1 .GT. LWRK) CALL ERRWRK('WESSWF.WESECV',KW1,LWRK)
C
C ....S-1/2 section
      KSWF1 = 1
      KU1   = KSWF1 + N2STA2
      KBIG  = KU1   + N2STA2
      KJBIG = KBIG  + NSTAT2
      KW2   = KJBIG + NSTAT2
      IF (KW2 .GT. LWRK) CALL ERRWRK('WESSWF.JACO',KW2,LWRK)
C
      IF (NASHT .LT. 1) THEN
         WRITE (LUPRI,'(//A)') ' STOP IN WESSWF: NASHT .LT. 1'
         CALL QUIT( ' WESSWF: NASHT .LT. 1')
      ENDIF
C
C Obtain string information for expanded graph
C
      CALL GETCIX(WRK(KCINDX),LSYM,LSYM,WRK(KW1),LW1,NOSYM)
C
C Start double loop over final states
C
      REWIND LUSPCW
      READ (LUSPCW)
      READ (LUSPCW)
      JI = 1
      SWF(JI) = D1
      DO 200 I = 2,NSTAT2
         CALL READT(LUSPCW,NCONF2,WRK(KCVEC1))
         CALL READT(LUSPCW,NCMOT,WRK(KCMO1))
C
         CALL EXPCID(WRK(KCVEC1),WRK(KCI1EX),LSYM2,1,LUSYMV(2),
     &               WRK(KCINDX),WRK(KW1),LW1,IPRECV)
         IF (IPRDEF .GE. 50) THEN
            WRITE (LUPRI,'(//A,I8)')
     &         ' CI vector no. 2 in WESSWF, length =',NCONFX
            WRITE (LUPRI,'(2(1X,5F15.12/))')
     &         (WRK(KCI2EX-1+III),III=1,NCONFX)
         END IF
C
         REWIND LUSPCW
         DO 100 J = 1,(I-1)
            CALL READT(LUSPCW,NCONF2,WRK(KCVEC2))
            CALL READT(LUSPCW,NCMOT,WRK(KCMO2))
C
            CALL EXPCID(WRK(KCVEC2),WRK(KCI2EX),LSYM2,1,LUSYMV(2),
     &                  WRK(KCINDX),WRK(KW1),LW1,IPRECV)
C
            JI = JI + 1
            SWF(JI) = WFOVLP(WRK(KCI1EX),WRK(KCMO1),WRK(KCI2EX),
     &                       WRK(KCMO2),SAO,WRK(KCINDX),WRK(KW1),LW1)
  100    CONTINUE
         JI = JI + 1
         SWF(JI) = D1
         READ (LUSPCW)
         READ (LUSPCW)
  200 CONTINUE
C
C     Print out SWF
C
      IF (IPRDEF .GE. 10) THEN
         WRITE (LUPRI,'(//A)') ' WESSWF-TEST: Spectral overlap matrix:'
         CALL OUTPAK (SWF,NSTAT2,1,LUPRI)
      END IF
C
C Calculation of TSWF = SWF**(-1/2)
C
C Do Ut (SWF) U = D and find TSWF = U (D**(-1/2)) Ut
C
C 1. GET D**(-1/2)   (DMAT)
C
      CALL DCOPY(NNSTA2,SWF,1,WRK(KSWF1),1)
      CALL DUNIT(WRK(KU1),NSTAT2)
      CALL JACO(WRK(KSWF1),WRK(KU1),NSTAT2,NSTAT2,NSTAT2,
     *          WRK(KBIG),WRK(KJBIG))
C
      CALL DZERO(TSWF,N2STA2)
      II = 0
      DO 600 I = 1,NSTAT2
         II = II + I
         DRT2 = SQRT( WRK((KSWF1-1) + II) )
         IF (ABS(DRT2) .LE. THRESH) GO TO 900
         TSWF(I,I) = D1 / DRT2
  600 CONTINUE
C
C GET S**(-1/2)
C
      CALL DGEMM('N','T',NSTAT2,NSTAT2,NSTAT2,1.D0,
     &           TSWF,      NSTAT2,
     &           WRK(KU1),  NSTAT2,0.D0,
     &           WRK(KSWF1),NSTAT2)
      CALL DGEMM('N','N',NSTAT2,NSTAT2,NSTAT2,1.D0,
     &           WRK(KU1),  NSTAT2,
     &           WRK(KSWF1),NSTAT2,0.D0,
     &           TSWF,      NSTAT2)
C
C     print out TSWF
C
      IF (IPRTRM .GE. 10) THEN
         WRITE (LUPRI,'(//A)') ' WESSWF-TEST: S-1/2 tranformation:'
         CALL OUTPUT (TSWF,1,NSTAT2,1,NSTAT2,NSTAT2,NSTAT2,1,LUPRI)
      END IF
C
      RETURN
C
  900 CONTINUE
         WRITE (LUPRI,'(///2A,1PD10.2/)')
     1   ' >>> WESSWF: SWF MATRIX TOO CLOSE TO SINGULARITY,',
     2   ' ORTHOGONALIZATION GIVES EIGENVALUE =',DRT2
         CALL QUIT( 'WESSHK: ERROR IN WESSWF')
C
      END
#  /* Deck bioras */
      SUBROUTINE BIORAS(CMO1,CMO2,VC1,VC2,SAO12,TIN1,TIN2,TRA1A,
     *                  TRA2A,WRK,LWRK)
C
C   20-OCT-1988 HA.   *** RAS VERSION ***
C
C Purpose:
C      TO COMPUTE PARTITION MATRICES AND BIORTHOGONAL ORBITALS
C      USING PER-AKE MALMQUISTS SCHEME FOR RAS WAVE FUNCTIONS.
C
C Input:  Molecular orbitals for two states (symmetry ordered)
C         (DISPLCED) OVERLAP MATRIX
C
C OUTPUT: PARTTITION MATRICES FOR WESTRACI, PSEUDOCORRESPONDING ORBITALS
C         FOR WESDM
C
#include "implicit.h"
      DIMENSION CMO1(*),CMO2(*),VC1(*),VC2(*),TIN1(*),TIN2(*),
     *          SAO12(*),TRA1A(*),TRA2A(*),WRK(*)
#include "iratdef.h"
C
C  PRIUNIT: LUPRI
C  WESINF: IPRBIO,OCEINT,...
C  INFORB: N2ASHX...
C
#include "priunit.h"
#include "wesinf.h"
#include "inforb.h"
C
C ***** core allocation *****
C
      N2OCCX = NOCCT * NOCCT
C
C...TRORAS KSMO is occupied-occupied block of total SMO
      KSMO  = 1
      KWTT  = KSMO  + N2OCCX
      LWT   = LWRK  - KWTT
      IF (KWTT.GE.LWRK) CALL ERRWRK('BIORAS-1',KWTT,LWRK)
C...FULL PARTITION MATRICES AND PSEUDOCORRESPONDING ORBITALS
      KCXA  = KWTT
      KCYB  = KCXA  + N2OCCX
      KVSS  = KCYB  + N2OCCX
      KTRA1 = KVSS  + NOCCT*NBAST
      KTRA2 = KTRA1 + N2OCCX
      KWP   = KTRA2 + N2OCCX
      LWP   = LWRK  - KWP
      IF (KWP.GE.LWRK) CALL ERRWRK('BIORAS-2',KWP,LWRK)
C
C
C *** SECTION 1 ***
C
C--First switch from CMO order to VC order
C
C  If OCE -calculation construct MO-verlap matrix (SMO) outside TRORAS
C
      IF (OCEINT) THEN
         WRITE(LUPRI,'(/A)') ' OCEINT NOT IMPLEMENTED (BIORAS)'
                 CALL QUIT ( ' OCEINT NOT IMPLEMENTED, BIORAS ')
c        CALL OCEOVL(CMO1,CMO2,SAO12,WRK(KSMO),WRK(KWTT),LWT)
C        CALL OCEOVL(CMO1,CMO2,SAO,SMO,WRK,LWRK)
      ENDIF
C
      CALL SYMMO(1,CMO1,VC1)
      CALL SYMMO(1,CMO2,VC2)
C
C
C   test of MO:s
C
      IF (IPRBIO .GE. 80) THEN
         WRITE(LUPRI,'(A/)') ' test write out of CMO1 before TRORAS '
            CALL PRORB(CMO1,.FALSE.,LUPRI)
         WRITE(LUPRI,'(A/)') ' test write out of VC1 before TRORAS '
           CALL OUTPUT(VC1,1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A/)') ' test write out of CMO2 before TRORAS '
            CALL PRORB(CMO2,.FALSE.,LUPRI)
         WRITE(LUPRI,'(A/)') ' test write out of VC2 before TRORAS '
           CALL OUTPUT(VC2,1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
      ENDIF
C
C
      CALL TRORAS(WRK(KSMO),SAO12,VC1,VC2,WRK(KWTT),LWT)
C
C     CALL TRORAS(SMO,SAO,VC1,VC2,WRK,LWRK)
C
C
C *** SECTION 2 ***
C
C PARTITION MATRICES
C
CDBG  write (LUPRI,*) 'nsym',nsym
CDBG  write (LUPRI,*) 'nish',nish
CDBG  write (LUPRI,*) 'nocc',nocc
CDBG  write (LUPRI,*) 'nas1',nas1
CDBG  write (LUPRI,*) 'nas2',nas2
CDBG  write (LUPRI,*) 'nas3',nas3
      CALL PART(WRK(KSMO),WRK(KTRA1),WRK(KTRA2),NSYM,NISH,NOCC,
     *          NAS1,NAS2,NAS3,LWP,WRK(KWP))
C
C     SUBROUTINE PART(SXY,TRA1,TRA2,NSYM,NISH,NOSH,NRAS1,NRAS2,
C    *                NRAS3,LWRK,WRK)
C
      IF (IPRBIO .GE. 10) THEN
         WRITE (LUPRI,'(//A)') ' BIORAS: TRA1 matrix '
         CALL OUTPUT(WRK(KTRA1),1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,LUPRI)
         WRITE (LUPRI,'(//A)') ' BIORAS: TRA2 matrix '
         CALL OUTPUT(WRK(KTRA2),1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,LUPRI)
      ENDIF
C
C  Take out active part of TRA matrices and the inactive diagonals
C
        DO 210 II = 1, NISHT
           TIN1(II) = WRK(KTRA1 + (II-1)*NOCCT + II -1)
           TIN2(II) = WRK(KTRA2 + (II-1)*NOCCT + II -1)
  210   CONTINUE
      IF (NASHT .GT. 0) THEN
        IOFF  = (NISHT * NOCCT) + NISHT
        IOFFT = 1
        DO 220 II = 1, NASHT
           CALL DCOPY(NASHT,WRK(KTRA1 + IOFF),1,TRA1A(IOFFT),1)
           CALL DCOPY(NASHT,WRK(KTRA2 + IOFF),1,TRA2A(IOFFT),1)
           IOFF  = IOFF  + NOCCT
           IOFFT = IOFFT + NASHT
  220   CONTINUE
      ENDIF
C
C *** SECTION 3 ***
C
C TRANSFORMATION MATRICES C AND PSEUDOCORRESPONDING ORBITALS
C
      CALL MKCXA(NSYM,NOCC,N2OCCX,WRK(KTRA1),WRK(KCXA))
      CALL MKCXA(NSYM,NOCC,N2OCCX,WRK(KTRA2),WRK(KCYB))
C
      IF (IPRBIO .GE. 10) THEN
         WRITE (LUPRI,'(//A)') ' BIORAS: CXA matrix '
         CALL OUTPUT(WRK(KCXA),1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,LUPRI)
         WRITE (LUPRI,'(//A)') ' BIORAS: CYB matrix '
         CALL OUTPUT(WRK(KCYB),1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,LUPRI)
      ENDIF
C
         CALL DGEMM('N','N',NBAST,NOCCT,NOCCT,1.D0,
     &              VC1,      NBAST,
     &              WRK(KCXA),NOCCT,0.D0,
     &              WRK(KVSS),NBAST)
C
      CALL DCOPY((NOCCT*NBAST),WRK(KVSS),1,VC1,1)
C
         CALL DGEMM('N','N',NBAST,NOCCT,NOCCT,1.D0,
     &              VC2,      NBAST,
     &              WRK(KCYB),NOCCT,0.D0,
     &              WRK(KVSS),NBAST)
C
       CALL DCOPY((NOCCT*NBAST),WRK(KVSS),1,VC2,1)
C
C--Reorder VC order to CMO order
C
      CALL SYMMO(2,CMO1,VC1)
      CALL SYMMO(2,CMO2,VC2)
C
C
      RETURN
      END
#  /* Deck troras */
      SUBROUTINE TRORAS(SMO,SAO,VC1,VC2,WRK,LWRK)
C
C 20-OCT 1988 HA: - Guga TROVLP MODIFIED FOR RAS
C Revised 25-Apr-1989 hjaaj
C
C Purpose:
C       Transform geometry displaced overlap matrix from AO to MO
C       basis.
C       RETURNS NOW ONLY THE NON-SECONDARY PART OF THE MO OVERLAP MATRIX
C
#include "implicit.h"
      DIMENSION VC1(NBAST,*),VC2(NBAST,*),SAO(*),WRK(*),SMO(NOCCT,*)
C
C  PRIUNIT: LUPRI
C  INFORB: NBAST,NORBT,...
C  INFDIM: NISHDI,NASHDI,...
C  WESINF: OCEINT,IPRBIO
C
#include "priunit.h"
#include "inforb.h"
#include "infdim.h"
#include "wesinf.h"
C
C
C TRANSFORM THE A-B SQUARE BLOCK OF SAO TO SMO
C (SAO IS TRIANGULAR WITH SIDE = 2 * NBAST? NO IT IS SQUARE BLOCKED
C
      IF (N2BASX.GE.LWRK) CALL ERRWRK('TROVLP',N2BASX,LWRK)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C  IF OCEINT THEN ORBITALS S ASSUMED TO BE OBTAINED FROM OCEOVL
C  (AS INPUT IN SMO) OCEINT NOT IMPLEMENTED YET (21/10 1988)
C
      IF ( .NOT. OCEINT) THEN
         CALL DGEMM('N','N',NBAST,NOCCT,NBAST,1.D0,
     &              SAO, NBAST,
     &              VC2, NBAST,0.D0,
     &              WRK, NBAST)
         CALL DGEMM('T','N',NOCCT,NOCCT,NBAST,1.D0,
     &              VC1, NBAST,
     &              WRK, NBAST,0.D0,
     &              SMO, NOCCT)
     
C
      ENDIF
      IF (IPRBIO .GE. 10) THEN
         WRITE(LUPRI,'(/A)')' TROVLP-TEST: SMO12 matrix (occ-occ block)'
         CALL OUTPUT(SMO,1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,LUPRI)
         call flshfo(lupri)
      END IF
      IF (IPRBIO .GE. 50) THEN
         WRITE(LUPRI,'(/A)')' TROVLP-TEST: VC1 matrix (bas-occ block)'
         CALL OUTPUT(VC1,1,NBAST,1,NOCCT,NBAST,NOCCT,1,LUPRI)
         WRITE(LUPRI,'(/A)')' TROVLP-TEST: VC2 matrix (bas-occ block)'
         CALL OUTPUT(VC2,1,NBAST,1,NOCCT,NBAST,NOCCT,1,LUPRI)
         WRITE(LUPRI,'(/A)')
     &   ' TROVLP-TEST: SAO*VC2  matrix (bas-occ block)'
         CALL OUTPUT(WRK,1,NBAST,1,NOCCT,NBAST,NOCCT,1,LUPRI)
      END IF
C
C
      RETURN
      END
#  /* Deck part */
      SUBROUTINE PART(SXY,TRA1,TRA2,NSYM,NISH,NOSH,NRAS1,NRAS2,
     *                NRAS3,LWRK,WRK)
C
C Uppsala 24/11 1988 (HA)
C
C MALMQUISTS ROUTINE, SLIGHTLY MODIFIED TO TAKE RAS
C In this version we assume no symmetry
C
C PURPOSE:
C   SXY CONTAINS THE NONSECONDARY PART OF THE MO OVERLAP MATRIX.
C   UPON RETURN, TRA1 AND TRA2 WILL CONTAIN THE COEFFICIENTS
C   FOR SEQUENTIAL SINGLE-ORBITAL TRANSFORMATIONS (VI.2, MY IJQC ARTICLE)
C   TO BIORTHONORMAL ORBITALS. SXY, TRA1 AND TRA2 ARE SYMMETRY-BLOCKED.
C                                     (MALMQUIST 84-04-04)
C                                     (UPDATED   84-05-08)
C                                     (IBM VERS  87-02-21)
C
#include "implicit.h"
      DIMENSION SXY(*),TRA1(*),TRA2(*), WRK(LWRK)
      DIMENSION NISH(8),NOSH(8),NRAS1(8),NRAS2(8),NRAS3(8)
#include "priunit.h"
C
C sum orbital data over symmetries as if only one IR
C
      NR1 = 0
      NR2 = 0
      NR3 = 0
      NO  = 0
      DO 10 ISYM  = 1, NSYM
        NO  = NO  + NOSH(ISYM)
        NR1 = NR1 + NRAS1(ISYM)
        NR2 = NR2 + NRAS2(ISYM)
        NR3 = NR3 + NRAS3(ISYM)
  10  CONTINUE
      IF(NO .EQ. 0) THEN
          WRITE (LUPRI, *) ' STOP in PART: No occupied orbitals '
          CALL QUIT(' No occupied orbitals in PART')
      ENDIF
C
C CORE ALLOCATE
C
      KSCR = 1
      KLW  = KSCR + NO**2 + 3*NO
      IF (KLW .GT. LWRK)  CALL ERRWRK('BIORAS.PART',KLW,LWRK)
C
C
      CALL PART1(NO,NR1,NR2,NR3,SXY(1),TRA1(1),TRA2(1),WRK(KSCR))
C
      RETURN
      END
#  /* Deck part1 */
      SUBROUTINE PART1(NDIM,NR1,NR2,NR3,SXY,B,A,SCR)
C
C 20/10 1988 (HA) SLIGHTLY MODIFIED FOR RAS
C
#include "implicit.h"
      DIMENSION SXY(NDIM,NDIM),A(NDIM,NDIM),B(NDIM,NDIM)
      DIMENSION SCR(NDIM,NDIM+3)
      DIMENSION ISTA(5),NSIZE(4)
#include "priunit.h"
C
C     NOTE: NSIZE IS A TEMPORARY MEASURE TO GIVE COMPATIBILITY
C     WITH LU2.
C
C
C     PURPOSE: SEE SUBROUTINE PART.
C     SUBDIVISION INTO TWO LEVELS OF ROUTINE CALLS IS MERELY TO
C     FACILITATE HANDLING OF SYMMETRY AND INDEXING.
C                                        (MALMQUIST 84-04-04)
C -- COPY THIS SYMMETRY BLOCK OF OVERLAP MATRIX INTO SCRATCH
C -- PUT UNIT MATRIX INTO A AND B:
      DO 10 J=1,NDIM
        DO  5 I=1,NDIM
          SCR(I,J)=SXY(I,J)
          A(I,J)=0.0D00
   5      B(I,J)=0.0D00
        A(J,J)=1.0D00
  10    B(J,J)=1.0D00
C -- SOLVE A:=SXY(INVERSE)*A
      CALL DOOL (NDIM,NDIM,NDIM,NDIM,SCR,A,DET,SCR(1,NDIM+1),
     &           SCR(1,NDIM+2),SCR(1,NDIM+3))
        write (lupri,*)'Sxy inverse'
        call output(a,1,ndim,1,ndim,ndim,ndim,1,lupri)
        call flshfo(lupri)
COLD  NA=NDIM-NI
      NI = NDIM - NR1 - NR2 - NR3
C -- DEFINE BEGINNING OF THE INACTIVE PARTITIONING:
      ISTA(1)=1
      NBLOCK=0
      IF(NI.GT.0) THEN
        NBLOCK=1
        NSIZE(1)=NI
        ISTA(2)=NI+1
      END IF
C -- DEFINE BEGINNING OF THE RAS PARTITIONINGS:
      IF(NR1.GT.0) THEN
        NBLOCK=NBLOCK+1
        NSIZE(NBLOCK)=NR1
        ISTA(NBLOCK+1)=ISTA(NBLOCK)+NR1
      END IF
      IF(NR2.GT.0) THEN
        NBLOCK=NBLOCK+1
        NSIZE(NBLOCK)=NR2
        ISTA(NBLOCK+1)=ISTA(NBLOCK)+NR2
      END IF
      IF(NR3.GT.0) THEN
        NBLOCK=NBLOCK+1
        NSIZE(NBLOCK)=NR3
        ISTA(NBLOCK+1)=ISTA(NBLOCK)+NR3
      END IF
C -- PERFORM A BLOCK LU PARTITIONING OF SYX IN MATRICES A AND B.
C -- FIRST COPY SYX INTO A, UNIT MATRIX IS ALREADY IN B:
C -- LOOP BACKWARDS OVER RAS PARTITIONS:
      DO 70 IBLK=NBLOCK,2,-1
C -- COPY DIAG BLOCK OF A INTO SCRATCH SPACE:
        DO 30 I=ISTA(IBLK),ISTA(IBLK+1)-1
          DO 30 J=ISTA(IBLK),ISTA(IBLK+1)-1
            SCR(J,I)=A(J,I)
  30    CONTINUE
C -- SOLVE (BLOCK ROW OF A):=(DIAG BLOCK OF A, INVERSE)*(BL-ROW OF A)
        NSZ=NSIZE(IBLK)
        IST=ISTA(IBLK)
        NL=IST-1
        CALL DOOL(NDIM,NDIM,NSZ,NL,SCR(IST,IST),A(IST,1),DET,
     &          SCR(1,NDIM+1),SCR(1,NDIM+2),SCR(1,NDIM+3))
        write (lupri,*)'Syx(inverse) Syx; iblk,ist,nsz',iblk,ist,nsz
        call output(a,ist,ist-1+nsz,1,nl,ndim,ndim,1,lupri)
        call flshfo(lupri)
C -- UPDATE EARLIER BLOCKS OF A:
        DO 50 J=1,ISTA(IBLK)-1
          DO 50 I=1,ISTA(IBLK)-1
            SUM=0.0D00
            DO 40 K=ISTA(IBLK),ISTA(IBLK+1)-1
              SUM=SUM+A(I,K)*A(K,J)
  40        CONTINUE
            A(I,J)=A(I,J)-SUM
  50    CONTINUE
C -- FINALLY, MOVE A ROW BLOCK OF A INTO B_transpose, REPLACING WITH ZEROES:
        DO 60 J=1,ISTA(IBLK)-1
          DO 60 I=ISTA(IBLK),ISTA(IBLK+1)-1
            B(J,I)=A(I,J)
            A(I,J)=0.0D00
  60    CONTINUE
  70  CONTINUE
C -- WE HAVE A BLOCK LU PARTITIONING NOW, OF THE FORM A*B=SYX.
C -- LU2 REQUIRES TENTATIVE CXA AND CYB WHERE CYB*(CXA TRANSPOSE)=SYX.
C -- SO, B was transposed above, BEFORE LU2 CALL:
      write (lupri,*) 'before LU2, A-matrix=CYB',ndim,nblock
      call output(a,1,ndim,1,ndim,ndim,ndim,1,lupri)
      write (lupri,*) 'before LU2, B-matrix=CXA',ndim,nblock
      call output(b,1,ndim,1,ndim,ndim,ndim,1,lupri)
      call flshfo(lupri)
      CALL LU2(NDIM,NBLOCK,NSIZE,B,A,SCR)
C
C     LU PARTITIONING OF THE MATRICES B AND A WAS DONE IN-PLACE.
C     NOW CHANGE SIGN OF LOWER-TRIANGULAR PARTS AND
C     INVERT UPPER-TRIANGULAR PARTS, AS INDICATED IN (VI.6):
C
      DO 200 I=2,NDIM
        DO 200 J=1,I-1
          A(I,J)=-A(I,J)
          B(I,J)=-B(I,J)
 200  CONTINUE
      DO 240 L=NDIM,1,-1
        A(L,L)=1.0D0/A(L,L)
        B(L,L)=1.0D0/B(L,L)
        DO 210 M=L+1,NDIM
          A(L,M)=A(L,L)*A(L,M)
 210      B(L,M)=B(L,L)*B(L,M)
        DO 230 K=1,L-1
          DO 220 M=L+1,NDIM
            A(K,M)=A(K,M)-A(K,L)*A(L,M)
 220        B(K,M)=B(K,M)-B(K,L)*B(L,M)
          A(K,L)=-A(L,L)*A(K,L)
          B(K,L)=-B(L,L)*B(K,L)
 230      CONTINUE
 240    CONTINUE
      RETURN
      END
#  /* Deck dool */
      SUBROUTINE DOOL(NDIM,MDIM,N,M,A,B,DET,IPIV,JPIV,BUF)
C
C     SOLVES THE MATRIX EQUATION AX=B BY DOOLITTLE´S METHOD
C     ACTUAL DIMENSIONS ARE N*N AND N*M
C     ALLOCATED DIMENSIONS ARE NDIM*NDIM AND NDIM*MDIM
C     A AND B ARE DESTROYED, AND X IS RETURNED AS MATRIX B
C                                       (MALMQUIST 82-11-12)
C                                        (UPDATE 83-09-28)
C
#include "implicit.h"
      DIMENSION A(NDIM,NDIM),B(NDIM,MDIM),IPIV(NDIM),JPIV(NDIM)
      DIMENSION BUF(NDIM)
C
C  EQUATION IS SOLVED BY FACTORIZING A=L*R IN SAME SPACE AS A.
C  PIVOTING IS ACHIEVED BY INDIRECT INDEXING.
C  FIRST PIVOTING INDICES ARE ASSIGNED START VALUES.
C
      DO I=1,N
        IPIV(I)=I
        JPIV(I)=I
      END DO
      DET=1.0D0
      DO 6 I=1,N
C     ... NOW FIND BETTER PIVOT ELEMENT:
        AMAX=-1.0D0
        DO 2 L=I,N
        DO 2 K=I,N
          AM=ABS(A(IPIV(K),JPIV(L)))
          IF(AM.GT.AMAX) THEN
            AMAX=AM
            IP=K
            JP=L
          END IF
 2      CONTINUE
        IF(IP.NE.I) THEN
          DET=-DET
          IDUM=IPIV(I)
          IPIV(I)=IPIV(IP)
          IPIV(IP)=IDUM
         END IF
        IF(JP.NE.I) THEN
          DET=-DET
          IDUM=JPIV(I)
          JPIV(I)=JPIV(JP)
          JPIV(JP)=IDUM
        END IF
        IP=IPIV(I)
        JP=JPIV(I)
        DIAG=A(IP,JP)
        BUF(I)=DIAG
        DET=DET*DIAG
        DO 5 K=I+1,N
          KP=IPIV(K)
          C=A(KP,JP)/DIAG
          A(KP,JP)=C
          DO 5 L=I+1,N
            LP=JPIV(L)
            A(KP,LP)=A(KP,LP)-C*A(IP,LP)
 5      CONTINUE
 6    CONTINUE
C
C  FIRST RESUBSTITUTION STEP:
C
      DO 7 J=1,M
      DO 7 I=2,N
        IP=IPIV(I)
        SUM=B(IP,J)
        DO K=1,I-1
          SUM=SUM-A(IP,JPIV(K))*B(IPIV(K),J)
        END DO
        B(IP,J)=SUM
 7    CONTINUE
C
C  SECOND RESUBSTITUTION STEP:
C
      DO 9 J=1,M
      DO 9 I=N,1,-1
        IP=IPIV(I)
        SUM=B(IP,J)
        DO 8 K=I+1,N
 8        SUM=SUM-A(IP,JPIV(K))*B(IPIV(K),J)
 9      B(IP,J)=SUM/BUF(I)
C
C  REORGANIZATION PART:
C
      DO 12 J=1,M
        DO 10 I=1,N
10      BUF(I)=B(IPIV(I),J)
        DO 11 I=1,N
11      B(JPIV(I),J)=BUF(I)
12    CONTINUE
      RETURN
      END
#  /* Deck lu2 */
      SUBROUTINE LU2 (NDIM,NBLOCK,NSIZE,CXA,CYB,SCR)
C
C  GIVES A SIMULTANEOUS LU-PARTITIONING OF MATRICES CXA,CYB IN THE
C  SENSE THAT CXA*X = L1*U1 AND CYB*X = L2*U2, WHERE X IS A BLOCK
C  UNITARY MATRIX. X IS NEVER FORMED, BUT AT EACH STEP SUCH A
C  TRANSFORMATION IS APPLIED THAT THE LU PARTITIONING CAN BE
C  APPLIED IN AN OPTIMAL WAY -- THIS IS REFERRED TO AS A
C  UNITARY PSEUDO-PIVOTATION. MATRICES CXA AND CYB ARE DESTROYED,
C  AND WILL CONTAIN THE NONTRIVIAL ELEMENTS OF THE TRIANGULAR
C  MATRICES.
C                                         ( MALMQUIST 84-01-16 )
C
#include "implicit.h"
      DIMENSION CXA(NDIM,NDIM),CYB(NDIM,NDIM),SCR(NDIM)
      DIMENSION NSIZE(NBLOCK)
C threshold put to 10-10  (18/9 -89)
      PARAMETER (THRDIAG=1.0D-10)
#include "priunit.h"
C
      NTOT=0
      DO 5 IBLOCK=1,NBLOCK
  5     NTOT=NTOT+NSIZE(IBLOCK)
      IF(NTOT.NE.NDIM) CALL QUIT('LU2: ERROR; NTOT.NE.NDIM.')

      IEND=0
      DO 210 IBLOCK=1,NBLOCK
        ISTA=IEND+1
        IEND=IEND+NSIZE(IBLOCK)
        DO 200 II=ISTA,IEND
          S1=0.0D0
          S2=0.0D0
          S3=0.0D0
          DO 10 J=II,IEND
            S1=S1+CXA(II,J)**2
            S2=S2+CYB(II,J)**2
            S3=S3+CXA(II,J)*CYB(II,J)
  10      CONTINUE
          IF((S1.LT.THRDIAG).OR.(S2.LT.THRDIAG)) THEN
            WRITE(LUPRI,*)' LU PARTITIONING IS TROUBLESOME.'
            WRITE(LUPRI,*)
     &      ' DIAGONAL ELEMENT no.',II,' block',IBLOCK,' IS TOO SMALL:'
            WRITE(LUPRI,*)' IN MATRIX CXA IT IS', SQRT(S1)
            WRITE(LUPRI,*)' IN MATRIX CYB IT IS', SQRT(S2)
            WRITE(LUPRI,*)' EVEN AFTER OPTIMAL PIVOT-TRANSFORMATION.'
            IF(S1.LT.THRDIAG) THEN
              WRITE(LUPRI,*)' MATRIX CXA IS'
              call output(cxa,1,ndim,1,ndim,ndim,ndim,1,LUPRI)
            END IF
            IF(S2.LT.THRDIAG) THEN
              WRITE(LUPRI,*)' MATRIX CYB IS'
              call output(cyb,1,ndim,1,ndim,ndim,ndim,1,LUPRI)
            END IF
            CALL QUIT('LU2: LU partitioning failed')
          END IF
C
          X1=1.0D0/SQRT(S1)
          X2=1.0D0/SIGN(SQRT(S2),S3)
          DO 20 I=II,IEND
  20        SCR(I)=X1*CXA(II,I)+X2*CYB(II,I)
          S=2*(1.0D0+X1*X2*S3)
          X=1.0D0/SIGN(SQRT(S),SCR(II))
          DO 30 I=II,IEND
  30        SCR(I)=X*SCR(I)
          X=1.0D0/(1.0D0+SCR(II))
          DO 70 I=1,IEND
            S2=0.0D0
            DO 50 J=II+1,IEND
  50          S2=S2+CXA(I,J)*SCR(J)
            S =S2+CXA(I,II)*SCR(II)
            S2=CXA(I,II)+X*S2
            CXA(I,II)=S
            DO 60 J=II+1,IEND
  60          CXA(I,J)=CXA(I,J)-S2*SCR(J)
  70      CONTINUE
          DO 110 I=1,IEND
            S2=0.0D0
            DO 90 J=II+1,IEND
  90          S2=S2+CYB(I,J)*SCR(J)
            S =S2+CYB(I,II)*SCR(II)
            S2=CYB(I,II)+X*S2
            CYB(I,II)=S
            DO 100 J=II+1,IEND
 100          CYB(I,J)=CYB(I,J)-S2*SCR(J)
 110      CONTINUE
          X=1.0D0/CXA(II,II)
          DO 120 I=II+1,IEND
            CXA(I,II)=X*CXA(I,II)
            DO 120 J=II+1,NTOT
  120         CXA(I,J)=CXA(I,J)-CXA(I,II)*CXA(II,J)
          X=1.0D0/CYB(II,II)
          DO 130 I=II+1,IEND
            CYB(I,II)=X*CYB(I,II)
            DO 130 J=II+1,NTOT
  130         CYB(I,J)=CYB(I,J)-CYB(I,II)*CYB(II,J)

  200   CONTINUE
  210 CONTINUE
      RETURN
      END
#  /* Deck mkcxa */
      SUBROUTINE MKCXA(NSYM,NOSH,NCXA,TRA,CXA)
C
C  changed for call to MKCXAL with one symmetry (ha 14/12 1988)
C
#include "implicit.h"
      DIMENSION TRA(NCXA),CXA(NCXA)
      DIMENSION NOSH(NSYM)
#include "priunit.h"
C
      NDIM = 0
      DO 10 ISYM = 1, NSYM
         NDIM = NDIM + NOSH(ISYM)
   10 CONTINUE
      IF (NDIM.EQ.0) THEN
         WRITE(LUPRI,*) ' STOP in MKCXA, no occupied orbitals'
         CALL QUIT(' No occ. orbitals in MKCXA ')
      ENDIF
C
      CALL MKCXAL(NDIM,TRA,CXA)
C
      RETURN
      END
#  /* Deck mkcxal */
      SUBROUTINE MKCXAL(NDIM,TRAL,CXAL)
C -- NOTE: THE FOLLOWING ROUTINE MAY NOT BE VECTORIZED.
C -- THERE IS A COMPILER BUG ON FORTRAN VERSION 2.2.0 (JUNE 1987).
C -- WRONG RESULTS PRODUCED EVEN ON VECTORIZED LEVEL 1.
#include "implicit.h"
      DIMENSION TRAL(NDIM,NDIM),CXAL(NDIM,NDIM)
      DO 20 J=1,NDIM
        DO 10 I=1,J-1
  10      CXAL(I,J)=0.0D00
  20    CXAL(J,J)=1.0D00

      DO 60 K=1,NDIM
        DO 40 I=1,K-1
          SUM=0.0D00
          DO 30 J=1,K-1
  30        SUM=SUM+CXAL(I,J)*TRAL(J,K)
  40      CXAL(I,K)=-SUM/TRAL(K,K)
        DO 50 I=K,NDIM
          SUM=TRAL(I,K)
          IF (I.EQ.K) SUM=-1.0D00
          DO 45 J=1,K-1
  45        SUM=SUM+CXAL(I,J)*TRAL(J,K)
  50      CXAL(I,K)=-SUM/TRAL(K,K)
  60  CONTINUE
      RETURN
      END
