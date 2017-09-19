C
C
C VERSION : $Revision: 1.1 $
C DATE    : $Date: 2012/08/14 13:01:22 $
C FILE    : dirana.F
C
C
      SUBROUTINE MPAB(A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT A TIMES B = C.
C  Written by George D. Purvis 1983
C  Revised 6-Nov-1984 Hans Joergen Aa. Jensen
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NCOLA.NE.NROWB) THEN
         WRITE (*,9000) NCOLA,NROWB
         CALL QUIT('ERROR, inconsistent matrix dimensions in MPAB')
      ENDIF
 9000 FORMAT(/' MPAB error: NCOLA .ne. NROWB, values =',2I10)
C
      IF (NROWB .EQ. 0) RETURN
      DO 40 J = 1,NCOLB
        DO 10 I = 1,NROWA
   10     C(I,J) = A(I,1)*B(1,J)
        DO 30 K = 2,NROWB
          IF (B(K,J).EQ.ZERO) GO TO 30
          BKJ = B(K,J)
          DO 20 I = 1,NROWA
   20       C(I,J) = BKJ*A(I,K) + C(I,J)
   30   CONTINUE
   40 CONTINUE
C
      RETURN
      END
C  /* Deck mpabt */
      SUBROUTINE MPABT (A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                  B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                  C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT A TIMES B(transposed) = C.
C  Written by George D. Purvis 1983
C  Revised 6-Nov-1984 Hans Joergen Aa. Jensen
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NCOLA.NE.NCOLB) THEN
        WRITE (*,9000) NCOLA,NCOLB
        CALL QUIT('ERROR, inconsistent matrix dimensions in MPABT')
      ENDIF
 9000 FORMAT(/' MPABT error: NCOLA .ne. NCOLB, values =',2I10)
C
      IF (NCOLB .EQ. 0) RETURN
      DO 10 J = 1,NROWB
        DO 10 I = 1,NROWA
   10     C(I,J) = A(I,1)*B(J,1)
C
      DO 40 K = 2,NCOLB
        DO 30 J = 1,NROWB
          IF (B(J,K).EQ.ZERO) GO TO 30
          DO 20 I = 1,NROWA
   20       C(I,J) = A(I,K)*B(J,K)+C(I,J)
   30   CONTINUE
   40 CONTINUE
      RETURN
      END
C  /* Deck mpatb */
      SUBROUTINE MPATB(A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                 B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                 C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT A(transpose) TIMES B = C.
C  Written by George D. Purvis 1983
C  Revised 6-Nov-1984 Hans Joergen Aa. Jensen
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NROWA.NE.NROWB) THEN
        WRITE (*,9000) NROWA,NROWB
        CALL QUIT('ERROR, inconsistent matrix dimensions in MPATB')
      ENDIF
 9000 FORMAT(/' MPATB error: NROWA .ne. NROWB, values =',2I10)
C
      DO 30 J = 1,NCOLB
        DO 20 I = 1,NCOLA
          CIJ = ZERO
          DO 10 K=1,NROWA
   10       CIJ = A(K,I)*B(K,J) + CIJ
   20     C(I,J) = CIJ
   30 CONTINUE
C
      RETURN
      END
C  /* Deck mpatbt */
      SUBROUTINE MPATBT(A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                  B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                  C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT A(transpose) TIMES B(transpose) = C.
C  6-April-1987 Hans Joergen Aa. Jensen
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NROWA .NE. NCOLB) THEN
         WRITE (*,9000) NROWA,NCOLB
         CALL QUIT('ERROR, inconsistent matrix dimensions in MPATBT')
      ENDIF
 9000 FORMAT(/' MPATBT error: NROWA .ne. NCOLB, values =',2I10)
C
      DO 30 J = 1,NROWB
         DO 20 I = 1,NCOLA
            CIJ = ZERO
            DO 10 K = 1,NROWA
   10          CIJ = A(K,I)*B(J,K) + CIJ
            C(I,J) = CIJ
   20    CONTINUE
   30 CONTINUE
C
      RETURN
      END
C  /* Deck mpapb */
      SUBROUTINE MPAPB(NDIMA,NCOLB,AP,B,C)
C
C  18-May-1987 hjaaj
C
C  C = AP * B
C
C  where C and B are (NDIMA,NCOLB) matrices
C  and AP is a packed, symmetric NDIMA by NDIMA matrix.
C
#include <implicit.h>
      DIMENSION AP(*),B(NDIMA,NCOLB),C(NDIMA,NCOLB)
      PARAMETER ( D0 = 0.0D0 )
C
      DO 500 K = 1,NCOLB
         IROW = 0
         DO 300 I = 1,NDIMA
            BIK = B(I,K)
            SUM = D0
            DO 100 J = 1,I-1
               SUM    = SUM    + AP(IROW+J) * B(J,K)
               C(J,K) = C(J,K) + AP(IROW+J) * BIK
  100       CONTINUE
            IROW   = IROW + I
            C(I,K) = SUM  + AP(IROW) * BIK
  300    CONTINUE
  500 CONTINUE
C
      RETURN
      END
C  /* Deck mpapv */
      SUBROUTINE MPAPV(N,AP,VIN,VOUT)
C
C  28-Apr-1985 hjaaj
C
C  VOUT = AP * VIN
C
C  where VOUT and VIN are N vectors
C  and AP is a packed, symmetric N by N matrix.
C
#include <implicit.h>
      DIMENSION AP(*),VIN(*),VOUT(*)
      PARAMETER (D0 = 0.0D0)
C
      KI = 0
      DO 300 I = 1,N
         VINI = VIN(I)
         SUM = D0
         DO 100 J = 1,I-1
            SUM = SUM + AP(KI+J) * VIN(J)
            VOUT(J) = VOUT(J) + AP(KI+J) * VINI
  100    CONTINUE
         KI = KI + I
         VOUT(I) = SUM + AP(KI) * VINI
  300 CONTINUE
C
      RETURN
      END
C  /* Deck ampab */
      SUBROUTINE AMPAB(A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                 B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                 C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT
C
C C = C + A TIMES B
C
C 28-Mar-1987 Hans Joergen Aa. Jensen, based on MPAB
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NCOLA.NE.NROWB) THEN
        WRITE (*,9000) NCOLA,NROWB
        CALL QUIT('ERROR, inconsistent matrix dimensions in AMPAB')
      ENDIF
 9000 FORMAT(/' AMPAB error: NCOLA .ne. NROWB, values =',2I10)
C
      DO 40 J = 1,NCOLB
         DO 30 K = 1,NROWB
            IF (B(K,J) .NE. ZERO) THEN
               BKJ = B(K,J)
               DO 20 I = 1,NROWA
   20             C(I,J) = BKJ*A(I,K) + C(I,J)
            END IF
   30    CONTINUE
   40 CONTINUE
C
      RETURN
      END
C  /* Deck ampabt */
      SUBROUTINE AMPABT (A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                   B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                   C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT
C
C C = C + A TIMES B(transpose)
C
C 28-Mar-1987 Hans Joergen Aa. Jensen, based on MPABT
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NCOLA.NE.NCOLB) THEN
         WRITE (*,9000) NCOLA,NCOLB
         CALL QUIT('ERROR, inconsistent matrix dimensions in AMPABT')
      ENDIF
 9000 FORMAT(/' AMPABT error: NCOLA .ne. NCOLB, values =',2I10)
C
      DO 40 K = 1,NCOLB
         DO 30 J = 1,NROWB
            IF (B(J,K) .NE. ZERO) THEN
               DO 20 I = 1,NROWA
   20             C(I,J) = A(I,K)*B(J,K) + C(I,J)
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
C  /* Deck ampatb */
      SUBROUTINE AMPATB(A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                  B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                  C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT
C
C C = C + A(transpose) TIMES B
C
C 28-Mar-1987 Hans Joergen Aa. Jensen, based on MPATB
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NROWA.NE.NROWB) THEN
         WRITE (*,9000) NROWA,NROWB
         CALL QUIT('ERROR, inconsistent matrix dimensions in AMPATB')
      ENDIF
 9000 FORMAT(/' AMPATB error: NROWA .ne. NROWB, values =',2I10)
C
      DO 30 J = 1,NCOLB
         DO 20 I = 1,NCOLA
            CIJ = ZERO
            DO 10 K = 1,NROWA
   10          CIJ = A(K,I)*B(K,J) + CIJ
   20       C(I,J) = CIJ + C(I,J)
   30 CONTINUE
C
      RETURN
      END
C  /* Deck ampapv */
      SUBROUTINE AMPAPV(N,AP,VIN,VOUT)
C
C 28-Mar-1987 hjaaj
C
C  VOUT = VOUT + AP * VIN
C
C  where VOUT and VIN are N vectors
C  and AP is a packed, symmetric N by N matrix
C
#include <implicit.h>
      DIMENSION AP(*),VIN(*),VOUT(*)
      PARAMETER (D0 = 0.0D0)
C
      KI = 0
      DO 300 I = 1,N
         VINI = VIN(I)
         SUM = D0
         DO 100 J = 1,I-1
            SUM = SUM + AP(KI+J) * VIN(J)
            VOUT(J) = VOUT(J) + AP(KI+J) * VINI
  100    CONTINUE
         KI = KI + I
         VOUT(I) = SUM + AP(KI) * VINI + VOUT(I)
  300 CONTINUE
C
      RETURN
      END
C  /* Deck smpab */
      SUBROUTINE SMPAB(A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                 B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                 C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT
C
C C = C - A TIMES B
C
C 28-Mar-1987 Hans Joergen Aa. Jensen, based on MPAB
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NCOLA.NE.NROWB) THEN
        WRITE (*,9000) NCOLA,NROWB
        CALL QUIT('ERROR, inconsistent matrix dimensions in SMPAB')
      ENDIF
 9000 FORMAT(/' SMPAB error: NCOLA .ne. NROWB, values =',2I10)
C
      DO 40 J = 1,NCOLB
         DO 30 K = 1,NROWB
            IF (B(K,J) .NE. ZERO) THEN
               BKJ = B(K,J)
               DO 20 I = 1,NROWA
   20             C(I,J) = C(I,J) - BKJ*A(I,K)
            END IF
   30    CONTINUE
   40 CONTINUE
C
      RETURN
      END
C  /* Deck smpabt */
      SUBROUTINE SMPABT (A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                   B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                   C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT
C
C C = C - A TIMES B(transpose)
C
C 28-Mar-1987 Hans Joergen Aa. Jensen, based on MPABT
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NCOLA.NE.NCOLB) THEN
         WRITE (*,9000) NCOLA,NCOLB
         CALL QUIT('ERROR, inconsistent matrix dimensions in SMPABT')
      ENDIF
 9000 FORMAT(/' SMPABT error: NCOLA .ne. NCOLB, values =',2I10)
C
      DO 40 K = 1,NCOLB
         DO 30 J = 1,NROWB
            IF (B(J,K) .NE. ZERO) THEN
               DO 20 I = 1,NROWA
   20             C(I,J) = C(I,J) - A(I,K)*B(J,K)
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
C  /* Deck smpatb */
      SUBROUTINE SMPATB(A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                  B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                  C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT
C
C C = C - A(transpose) TIMES B
C
C 28-Mar-1987 Hans Joergen Aa. Jensen, based on MPATB
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NROWA.NE.NROWB) THEN
         WRITE (*,9000) NROWA,NROWB
         CALL QUIT('ERROR, inconsistent matrix dimensions in SMPATB')
      ENDIF
 9000 FORMAT(/' SMPATB error: NROWA .ne. NROWB, values =',2I10)
C
      DO 30 J = 1,NCOLB
         DO 20 I = 1,NCOLA
            CIJ = ZERO
            DO 10 K = 1,NROWA
   10          CIJ = A(K,I)*B(K,J) + CIJ
   20       C(I,J) = C(I,J) - CIJ
   30 CONTINUE
C
      RETURN
      END
C  /* Deck smpapv */
      SUBROUTINE SMPAPV(N,AP,VIN,VOUT)
C
C 28-Mar-1987 hjaaj
C
C  VOUT = VOUT - AP * VIN
C
C  where VOUT and VIN are N vectors
C  and AP is a packed, symmetric N by N matrix
C
#include <implicit.h>
      DIMENSION AP(*),VIN(*),VOUT(*)
      PARAMETER (D0 = 0.0D0)
C
      KI = 0
      DO 300 I = 1,N
         VINI = VIN(I)
         SUM  = D0
         DO 100 J = 1,I-1
            SUM     = SUM     + AP(KI+J) * VIN(J)
            VOUT(J) = VOUT(J) - AP(KI+J) * VINI
  100    CONTINUE
         KI = KI + I
         VOUT(I) = VOUT(I) - (SUM + AP(KI) * VINI)
  300 CONTINUE
C
      RETURN
      END
C  /* Deck xmpab */
      SUBROUTINE XMPAB(X,
     1                 A,NROWA,NCOLA,NRDIMA,NCDIMA,
     2                 B,NROWB,NCOLB,NRDIMB,NCDIMB,
     3                 C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT
C
C C = C + SCALAR * A TIMES B
C
C 28-Mar-1987 Hans Joergen Aa. Jensen, based on MPAB
C 16-Mar-1990 Olav Vahtras, based on MPAB,AMPAB,SMPAB
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NCOLA.NE.NROWB) THEN
        WRITE (*,9000) NCOLA,NROWB
        CALL QUIT('ERROR, inconsistent matrix dimensions in XMPAB')
      ENDIF
 9000 FORMAT(/' XMPAB error: NCOLA .ne. NROWB, values =',2I10)
C
      DO 40 J = 1,NCOLB
         DO 30 K = 1,NROWB
            IF (B(K,J) .NE. ZERO) THEN
               XBKJ = X*B(K,J)
               DO 20 I = 1,NROWA
   20             C(I,J) = XBKJ*A(I,K) + C(I,J)
            END IF
   30    CONTINUE
   40 CONTINUE
C
      RETURN
      END
C  /* Deck xmpabt */
      SUBROUTINE XMPABT (X,
     1                   A,NROWA,NCOLA,NRDIMA,NCDIMA,
     2                   B,NROWB,NCOLB,NRDIMB,NCDIMB,
     3                   C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT
C
C C = C + SCALAR * A TIMES B(transpose)
C
C 28-Mar-1987 Hans Joergen Aa. Jensen, based on MPABT
C 16-Mar-1990 Olav Vahtras, based on MPABT,AMPABT,SMPABT
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NCOLA.NE.NCOLB) THEN
         WRITE (*,9000) NCOLA,NCOLB
         CALL QUIT('ERROR, inconsistent matrix dimensions in XMPABT')
      ENDIF
 9000 FORMAT(/' XMPABT error: NCOLA .ne. NCOLB, values =',2I10)
C
      DO 40 K = 1,NCOLB
         DO 30 J = 1,NROWB
            IF (B(J,K) .NE. ZERO) THEN
               XBJK = X*B(J,K)
               DO 20 I = 1,NROWA
   20             C(I,J) = XBJK*A(I,K) + C(I,J)
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
C  /* Deck xmpatb */
      SUBROUTINE XMPATB(X,
     1                  A,NROWA,NCOLA,NRDIMA,NCDIMA,
     2                  B,NROWB,NCOLB,NRDIMB,NCDIMB,
     3                  C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT
C
C C = C + SCALAR * A(transpose) TIMES B
C
C 28-Mar-1987 Hans Joergen Aa. Jensen, based on MPATB
C 16-Mar-1990 Olav Vahtras, based on MPATB,AMPATB,SMPATB
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NROWA.NE.NROWB) THEN
         WRITE (*,9000) NROWA,NROWB
         CALL QUIT('ERROR, inconsistent matrix dimensions in XMPATB')
      ENDIF
 9000 FORMAT(/' XMPATB error: NROWA .ne. NROWB, values =',2I10)
C
      DO 30 J = 1,NCOLB
         DO 20 I = 1,NCOLA
            CIJ = ZERO
            DO 10 K = 1,NROWA
   10          CIJ = A(K,I)*B(K,J) + CIJ
   20       C(I,J) = X*CIJ + C(I,J)
   30 CONTINUE
C
      RETURN
      END
C  /* Deck mcopy */
      SUBROUTINE MCOPY(NROWA,NCOLA,A,NRDIMA,B,NRDIMB)
C
C  1-Nov-1989 HJAaJ
C
C Copy A to B
C
#include <implicit.h>
      DIMENSION A(NRDIMA,NCOLA), B(NRDIMB,NCOLA)
C
      DO 200 J = 1,NCOLA
         DO 100 I = 1,NROWA
            B(I,J) = A(I,J)
  100    CONTINUE
  200 CONTINUE
      RETURN
      END
C  /* Deck mtrsp */
      SUBROUTINE MTRSP(NROWA,NCOLA,A,NRDIMA,B,NRDIMB)
C
C 30-Jul-1986 HJAaJ
C New version 1-Nov-1989 hjaaj
C 900108-hjaaj: block with NBLK for reduced paging
C when virtual memory
C
C Calculate B = A(transposed)
C
#include <implicit.h>
      DIMENSION A(NRDIMA,NCOLA), B(NRDIMB,NROWA)
      PARAMETER (NBLK = 128)
C
      DO 400 JBLK = 1,NCOLA,NBLK
         JEND = MIN(NCOLA,JBLK-1+NBLK)
         DO 300 IBLK = 1,NROWA,NBLK
            IEND = MIN(NROWA,IBLK-1+NBLK)
            DO 200 J = JBLK,JEND
               DO 100 I = IBLK,IEND
                  B(J,I) = A(I,J)
  100          CONTINUE
  200       CONTINUE
  300    CONTINUE
  400 CONTINUE
      RETURN
      END
C  /* Deck uthu */
      SUBROUTINE UTHU(H,HT,U,WRK,NAO,NMO)
C
C /VER 2/ 28-Apr-1985 Hans Jorgen Aa. Jensen
C
C CALLS: MPAPV, DDOT
C
C This subroutine transforms the symmetric, packed matrix H of
C dimension NAO to the symmetric, packed HT of dimension NMO.
C U(NAO,NMO) is the transformation matrix.
C
C                                      T
C  HT( ij ) =  SUM (k,l = 1,...,NAO)  U (i,k) H(kl) U(l,j)
C
C This routine is optimized with respect to vector operations.
C
#include <implicit.h>
      DIMENSION H(*),HT(*),U(NAO,*),WRK(*)
C
      J1 = 1
      DO 400 J = 1,NMO
         CALL MPAPV(NAO,H,U(1,J),WRK)
         CALL MPATB(U     ,NAO,J,NAO,NMO,
     &              WRK   ,NAO,1,NAO,1,
     &              HT(J1),J   ,1)
C        DO 200 K = 1,J
C           HT(J1-1+K) = DDOT(NAO,WRK(1),1,U(1,K),1)
C 200    CONTINUE
         J1 = J1 + J
  400 CONTINUE
C
      RETURN
      END
C  /* Deck authu */
      SUBROUTINE AUTHU(H,HT,U,WRK,NAO,NMO)
C
C /VER 2/ 28-Apr-1985 Hans Jorgen Aa. Jensen
C
C CALLS: MPAPV, DDOT
C
C "Add UT H U" to HT.
C This subroutine transforms the symmetric, packed matrix H of
C dimension NAO and adds the result to the symmetric, packed HT
C of dimension NMO.
C U(NAO,NMO) is the transformation matrix.
C
C                                      T
C  HT( ij ) =  SUM (k,l = 1,...,NAO)  U (i,k) H(kl) U(l,j)
C
C This routine is optimized with respect to vector operations.
C
#include <implicit.h>
      DIMENSION H(*),HT(*),U(NAO,*),WRK(*)
C
      J1 = 1
      DO 400 J = 1,NMO
         CALL MPAPV (NAO,H,U(1,J),WRK)
         CALL AMPATB(U     ,NAO,J,NAO,NMO,
     &               WRK   ,NAO,1,NAO,1,
     &               HT(J1),J   ,1)
C        DO 200 K = 1,J
C           HT(J1-1+K) = HT(J1-1+K) + DDOT(NAO,WRK(1),1,U(1,K),1)
C 200    CONTINUE
         J1 = J1 + J
  400 CONTINUE
C
      RETURN
      END
C  /* Deck utau */
      SUBROUTINE UTAU(AP,APT,U,WRK,NAO,NMO)
C
C 18-Feb-1987 Hans Jorgen Aa. Jensen
C
C CALLS: DDOT
C
C This subroutine transforms the antisymmetric, packed matrix AP of
C dimension NAO to the antisymmetric, packed APT of dimension NMO.
C U(NAO,NMO) is the transformation matrix.
C
C                                       T
C  APT( ij ) =  SUM (k,l = 1,...,NAO)  U (i,k) AP(kl) U(l,j)
C
C This routine is optimized with respect to vector operations.
C
#include <implicit.h>
      DIMENSION AP(*),APT(*),U(NAO,*),WRK(*)
      PARAMETER ( D0 = 0.0D0 )
C
      J1 = 1
      DO 400 JMO = 1,NMO
C        ... First UT * AP
         KI = 0
         DO 200 I = 1,NAO
            UIJMO = U(I,JMO)
            SUM   = D0
            DO 100 J = 1,I-1
               SUM    = SUM    + AP(KI+J) * U(J,JMO)
               WRK(J) = WRK(J) + AP(KI+J) * UIJMO
  100       CONTINUE
            WRK(I) = -SUM
C
            KI = KI + I
            IF (AP(KI) .NE. D0) THEN
               CALL QUIT('UTAU: AP not antisymmetric '//
     &                   '(diagonal not zero).')
            END IF
  200    CONTINUE
C        ... then (UT AP) * U
         CALL MPATB(U      ,NAO,JMO,NAO,NMO,
     &              WRK    ,NAO,1,NAO,1,
     &              APT(J1),JMO ,1)
C        DO 300 K = 1,JMO
C           APT(J1-1+K) = DDOT(NAO,WRK(1),1,U(1,K),1)
C 300    CONTINUE
         J1 = J1 + JMO
  400 CONTINUE
C
      RETURN
      END
C  /* Deck uhut */
      SUBROUTINE UHUT(H,HT,U,NAO,NMO)
C
C /VER 1/ 28-Mar-1987 Hans Jorgen Aa. Jensen
C 900720-hjaaj: NOSYNC compiler option.  Corrected fatal error
C (this routine can never have been used! HT would have been wrong).
C
C This subroutine backtransforms the symmetric, packed matrix H of
C dimension NMO to the symmetric, packed HT of dimension NAO
C (contravariant basis).
C U(NAO,NMO) is the transformation matrix.
C The transformation is here performed as one NAO**2*NMO**2/4 process
C instead of two successive N**3 transformations.
C
C                                                   T
C  HT( ij ) =  SUM (k,l = 1,...,NMO)  U(i,k) H(kl) U (l,j)
C
#include <implicit.h>
      DIMENSION H(*),HT(*),U(NAO,*)
      PARAMETER (D0 = 0.0D0)
C
      NNAO = NAO*(NAO+1)/2
      DO 100 IJ = 1,NNAO
         HT(IJ) = D0
  100 CONTINUE
C
      KL = 0
      DO 420 K = 1,NMO
         DO 410 L = 1,K
            KL   = KL + 1
            HKL  = H(KL)
#include <nosync.h>
            DO 240 I = 1,NAO
               IROW = (I*I-I)/2
               IF (K .NE. L) THEN
                  UHIK = U(I,L) * HKL
                  UHIL = U(I,K) * HKL
                  DO 210 J = 1,I
                     HT(IROW + J) = HT(IROW + J) + UHIL * U(J,L)
     *                                           + UHIK * U(J,K)
  210             CONTINUE
               ELSE
                  UHIK = U(I,K) * HKL
                  DO 220 J = 1,I
                     HT(IROW + J) = HT(IROW + J) + UHIK * U(J,K)
  220             CONTINUE
               END IF
  240       CONTINUE
  410    CONTINUE
  420 CONTINUE
C
      RETURN
      END
C  /* Deck auhut */
      SUBROUTINE AUHUT(H,HT,U,NAO,NMO)
C
C /VER 1/ 28-Mar-1987 Hans Jorgen Aa. Jensen
C 900720-hjaaj: NOSYNC compiler option.  Corrected fatal error
C (this routine can never have been used! HT would have been wrong).
C
C This subroutine backtransforms the symmetric, packed matrix H of
C dimension NMO to the symmetric, packed HT of dimension NAO
C (contravariant basis).
C U(NAO,NMO) is the transformation matrix.
C The transformation is here performed as one NAO**2*NMO**2/4 process
C instead of two successive N**3 transformations.
C
C                                                              T
C  HT( ij ) =  HT( ij ) + SUM (k,l = 1,...,NMO)  U(i,k) H(kl) U (l,j)
C
#include <implicit.h>
      DIMENSION H(*),HT(*),U(NAO,*)
C
      KL = 0
      DO 420 K = 1,NMO
         DO 410 L = 1,K
            KL   = KL + 1
            HKL  = H(KL)
#include <nosync.h>
            DO 240 I = 1,NAO
               IROW = (I*I-I)/2
               IF (K .NE. L) THEN
                  UHIK = U(I,L) * HKL
                  UHIL = U(I,K) * HKL
                  DO 210 J = 1,I
                     HT(IROW + J) = HT(IROW + J) + UHIL * U(J,L)
     *                                           + UHIK * U(J,K)
  210             CONTINUE
               ELSE
                  UHIK = U(I,K) * HKL
                  DO 220 J = 1,I
                     HT(IROW + J) = HT(IROW + J) + UHIK * U(J,K)
  220             CONTINUE
               END IF
  240       CONTINUE
  410    CONTINUE
  420 CONTINUE
C
      RETURN
      END
C  /* Deck uthub */
      SUBROUTINE UTHUB(HAO,HMO,U,WRK,NSYM,NAO,NMO)
C
C 20-May-1993 Hans Jorgen Aa. Jensen (based on UTHU)
C
C CALLS: MPAPV, DDOT
C
C This subroutine transforms the NSYM symmetric, packed matrix HAO
C of dimension NAO(isym) to the symmetric, packed HMO of dimension
C NMO(isym).
C U(NAO(isym),NMO(isym)) is the transformation matrix.
C
C                                      T
C  HMO( ij ) =  SUM (k,l = 1,...,NAO)  U (i,k) HAO(kl) U(l,j)
C
C This routine is optimized with respect to vector operations.
C
#include <implicit.h>
      DIMENSION HAO(*),HMO(*),U(*),WRK(*), NAO(NSYM), NMO(NSYM)
C
      IUOFF  = 1
      IHAOFF = 1
      IHMOFF = 1
      DO 600 ISYM = 1,NSYM
         NMOI = NMO(ISYM)
         NAOI = NAO(ISYM)
C
         J1 = IHMOFF
         DO 400 J = 1,NMOI
            CALL MPAPV(NAOI,HAO(IHAOFF),U(IUOFF+(J-1)*NAOI),WRK)
            CALL MPATB(U(IUOFF),NAOI,J,NAOI,NMOI,
     &                 WRK     ,NAOI,1,NAOI,1,
     &                 HMO(J1) ,J   ,1)
C           DO 200 K = 1,J
C              HMO(J1-1+K) = DDOT(NAOI,WRK(1),1,U(IUOFF+(K-1)*NAOI),1)
C 200       CONTINUE
            J1 = J1 + J
  400    CONTINUE
C
         IUOFF  = IUOFF + NAOI*NMOI
         IHAOFF = IHAOFF + (NAOI*(NAOI+1))/2
         IHMOFF = IHMOFF + (NMOI*(NMOI+1))/2
  600 CONTINUE
C
      RETURN
      END
C  /* Deck dv3dot */
      FUNCTION DV3DOT(N,V1,V2,V3)
C
C     7-Aug-1986 hjaaj
C
#include <implicit.h>
      DIMENSION V1(*), V2(*), V3(*)
      PARAMETER ( D0 = 0.0D0 )
      T = D0
      DO 100 K = 1,N
         T = T + V1(K) * V2(K) * V3(K)
  100 CONTINUE
      DV3DOT = T
      RETURN
      END
C  /* Deck ddotm */
      FUNCTION DDOTM(NU,NV,NRDIMA,U,A,V)
C
C 28-Feb-1989 Hans Joergen Aa. Jensen
C
C Calculate dot product with metric: DDOTM = ut A v
C
#include <implicit.h>
      DIMENSION U(NU), A(NRDIMA,NV), V(NV)
C
      SUM = 0.0D0
      DO 300 L = 1,NV
         SUM1 = 0.0D0
         DO 200 K = 1,NU
            SUM1 = SUM1 + U(K) * A(K,L)
  200    CONTINUE
         SUM = SUM + V(L)*SUM1
  300 CONTINUE
      DDOTM = SUM
C
      RETURN
      END
