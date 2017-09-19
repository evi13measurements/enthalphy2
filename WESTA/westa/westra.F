!
!  File: westa/westra.F
!
#  /* Deck westraci */
      SUBROUTINE WESTRACI(NCVEC,CVEC,LCVEC,UMO,XNDXCI,WRK,LFREE)
C
C New WESTRACI 040923 hjaaj, based on TRACI.
C
C Driver for TRACI2 using determinant routines
C Version with CSF business
C
C WESTA version: UMO on input is single orbital transformation matrix.
C                (in TRACI SOT matrix is made inside routine)½
C
C  JULY 15 1988 Jeppe Olsen
C
C
C MOTECC-90: The algorithms used in this module, TRACI, are
C            described in Chapter 8 Section D.5 of MOTECC-90
C            "Counter Rotations of CI coefficients"
C
#include "implicit.h"
      DIMENSION CVEC(LCVEC,NCVEC), UMO(*), XNDXCI(*), WRK(LFREE)
C
C Used from common blocks:
C  PRIUNIT : LUPRI
C  WESINF : IPRDEF
C  INFINP : MCTYPE,LSYM
C  INFORB : MULD2H(8,8), NASHT, N2ASHX
C  INFVAR : NCONF
C  DETBAS : KIASTR,KIBSTR, ...
C  STRNUM : NAEL,NBEL, ...
C  CIINFO : NDTASM
C  CSFBAS : Pointers for CSF information
C  CBREOR : SLREOR
C
#include "wesinf.h"
#include "priunit.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "detbas.h"
#include "mxpdim.h"
#include "strnum.h"
#include "csfbas.h"
#include "ciinfo.h"
#include "cbreor.h"
C..
      LOGICAL CSFEXP
C..
C
      IPRTCI = IPRDEF
C
C
      CSFEXP = .NOT.FLAG(27)
      IF (CSFEXP) THEN
        NCDET = NDTASM(LSYM)
        IF (NCONF .NE. NCSASM(LSYM)) THEN
           WRITE (LUPRI,*) 'WESTRACI ERROR, NCONF .ne. NCSASM(LSYM):',
     *        NCONF, NCSASM(LSYM)
           CALL QUIT('WESTRACI ERROR, NCONF .ne. NCSASM(LSYM)')
        END IF
      ELSE
        NCDET = NCONF
        IF (NCONF .NE. NDTASM(LSYM)) THEN
           WRITE (LUPRI,*) 'WESTRACI ERROR, NCONF .ne. NDTASM(LSYM):',
     *        NCONF, NDTASM(LSYM)
           CALL QUIT('WESTRACI ERROR, NCONF .ne. NDTASM(LSYM)')
        END IF
      END IF
c Local memory
      KLSOT  = 1
      KLC2   = KLSOT + N2ASHX
      KLFREE = KLC2  + MAX(N2ASHX, NCDET )
      IF(CSFEXP) THEN
        KLDET1 = KLFREE
        KLDET2 = KLDET1 + NCDET
        KLFREE = KLDET2 + NCDET
      END IF
C..   note: C2 only needs max of nia*nib (see TRACI1)
      IF ( KLFREE-1 .GT. LFREE ) CALL ERRWRK('WESTRACI',KLFREE-1,LFREE)
c orbital tranformation matrix
      IF (.NOT. SLREOR) THEN
         CALL DCOPY(N2ASHX,UMO,1,WRK(KLSOT),1)
      ELSE
C        reorder from Sirius order to Lunar order
         CALL REORMT(UMO,WRK(KLSOT),NASHT,NASHT,XNDXCI(KLTSOB),
     &               XNDXCI(KLTSOB) )
      END IF
C..
c     IF (IPRTCI .GT. 10 ) THEN
         WRITE(LUPRI,'(/A)')
     &   ' WESTRACI: SINGLE ORBITAL TRANSFORMATION MATRIX '
         CALL WRTMAT(WRK(KLSOT),NASHT,NASHT,NASHT,NASHT,0)
c     END IF
C
C.. TRANSPOSE SOT MATRIX IN ACCORDANCE WITH LUNAR DESCRIPTION
      CALL DGETRN(WRK(KLSOT),NASHT,NASHT)
C
C.. Store CI offsets in C arrays
      IPRXXX = MAX(0,IPRTCI - 10)
      CALL CIOFF(LSYM,1,XNDXCI,IPRXXX)
C     CALL CIOFF(IREFSM,ICORHC,XNDXCI,NTEST)
C
      DO 450 IVEC = 1, NCVEC
        IF(CSFEXP) THEN
c
c CSF to DET transformation, rotate in DET basis, DET to CSF transformation
c
          CALL COPVEC(CVEC(1,IVEC),WRK(KLDET1),NCONF)
          CALL CSDTVC(WRK(KLDET1),WRK(KLDET2),1,XNDXCI(KDTOC),
     &                XNDXCI(KICTS(1)),LSYM,0,IPRXXX)
C         CSDTVC(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,IREFSM,ICOPY,NTEST)
          CALL TRACI2(WRK(KLSOT),NASHT,LSYM,NAEL,XNDXCI(KIASTR),
     &                XNDXCI(KTAIJ),XNDXCI(KTATO),XNDXCI(KTASYM),
     &                XNDXCI(KSTASA),XNDXCI(KSTBAA),NBEL,
     &                XNDXCI(KIBSTR),XNDXCI(KTBIJ),
     &                XNDXCI(KTBTO),XNDXCI(KTBSYM),XNDXCI(KSTASB),
     &                XNDXCI(KSTBAB),NASTR,NBSTR,XNDXCI(KIANNI),
     &                WRK(KLDET2),XNDXCI(KCOFF),WRK(KLC2),
     &                NOCTPA,XNDXCI(KNSSOA),XNDXCI(KISSOA),
     &                XNDXCI(KTPFSA),NOCTPB,XNDXCI(KNSSOB),
     &                XNDXCI(KISSOB),
     &                XNDXCI(KTPFSB),XNDXCI(KIOCOC),XNDXCI(KICOOS),
     &                XNDXCI(KCDTAS),MULD2H,MAXSYM,IPRTCI)
          CALL CSDTVC(WRK(KLDET1),WRK(KLDET2),2,XNDXCI(KDTOC),
     &                XNDXCI(KICTS(1)),LSYM,0,IPRXXX)
C         CSDTVC(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,IREFSM,ICOPY,NTEST)
          CALL COPVEC(WRK(KLDET1),CVEC(1,IVEC),NCONF)
        ELSE
C
C.. Rotate in DET expansion
C
          CALL TRACI2(WRK(KLSOT),NASHT,LSYM,NAEL,XNDXCI(KIASTR),
     &                XNDXCI(KTAIJ),XNDXCI(KTATO),XNDXCI(KTASYM),
     &                XNDXCI(KSTASA),XNDXCI(KSTBAA),NBEL,
     &                XNDXCI(KIBSTR),XNDXCI(KTBIJ),
     &                XNDXCI(KTBTO),XNDXCI(KTBSYM),XNDXCI(KSTASB),
     &                XNDXCI(KSTBAB),NASTR,NBSTR,XNDXCI(KIANNI),
     &                CVEC(1,IVEC),XNDXCI(KCOFF),WRK(KLC2),
     &                NOCTPA,XNDXCI(KNSSOA),XNDXCI(KISSOA),
     &                XNDXCI(KTPFSA),NOCTPB,XNDXCI(KNSSOB),
     &                XNDXCI(KISSOB),
     &                XNDXCI(KTPFSB),XNDXCI(KIOCOC),XNDXCI(KICOOS),
     &                XNDXCI(KCDTAS),MULD2H,MAXSYM,IPRTCI)
        END IF
  450 CONTINUE
C
      RETURN
      END
