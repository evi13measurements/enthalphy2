C     FILE : wesnad.F
C
      SUBROUTINE WESNAD(WRK,LWRK)
C
C  (c) Hans Joergen Aa. Jensen

C  (original version May 1986)
C  (revised Sep 2004 hjaaj)
C
C  WESTA control module for
C  calculation of non-adiabatic coupling elements (NSTAT2.eq.1) or
C  calculation of adiabatic correction elements (NSTAT2.eq.0).
C
C
C  ********************* OUTLAY OF MODULE *****************************
C
C
C  The module performes a numerical calculation of the matrix elements
C
C  < WF(i) | d/dQ WF(j) >  and  < d/dQ WF(i) | d/dQ WF(j) >.
C
C  It requires the results of 4 MCSCF calculations;
C  for two states at two geometries.
C
C  Each overlap element is obtained though a 4-step calculation
C
C    Part 1: Atomic orbital overlap over displaced geometries
C
C    Part 2: Transformation to biorthogonal orbitals
C
C    Part 3: Configuration interaction in biorthogonal basis
C
C    Part 4: MO scaling
C
C  The coupling elements are obtained form the overlap elements and
C  the displacements. Compensation factors are calculated if symmetric
C  orthogonalization is required ("SYMORT")
C
C  ********************************************************************
C
C Declarations
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION WRK(LWRK)
#include "codata.h"
      PARAMETER (D0 = 0.D0, D1 = 1.0D0, D2 = 2.0D0,
     &           DP5 = 0.5D0, DM1 = -1.0D0)
C
C WESINF : DISPL,NGEOM,SYMORT, ?
C INFORB : N2BASX...
C INFVAR : NCONF..
C INFDIM : LCINDX..
C INFINP : NACTEL
C
#include "maxorb.h"
#include "wesinf.h"
#include "inforb.h"
#include "infvar.h"
#include "infdim.h"
#include "infinp.h"
C
      IROW(I) = ( I * (I-1) ) / 2
C
C
C
      CALL HEADER('>>> OUTPUT FROM WESNAD <<<',-1)
      IF (NGEOM .NE. 2) THEN
         WRITE (LUPRI, '(///A/A//)')
     *      ' *** ERROR (WESNAD), non-adiabatic coupling elements',
     *      '                     require two geometries. Nothing done.'
         RETURN
      END IF
C
C     Core allocation
C
      NBSG   = NBAST * NGEOM
      NNBSG  = NBSG * (NBSG + 1) /2
C...NADOVL...READMW...BIORAS
      KSAOAA = 1
      KSAOAB = KSAOAA + N2BASX
      KSAOBB = KSAOAB + N2BASX
C   NADOVL
      KW0    = KSAOBB + N2BASX
      LW0    = LWRK   - KW0
      IF (KW0 .GT. LWRK) CALL ERRWRK('WESNAD.NADOVL',KW0,LWRK)
C   WESGMO...WFOVLP
      KCMOA1 = KW0
      KCMOA2 = KCMOA1 + NCMOT
      KCMOB1 = KCMOA2 + NSTAT2*NCMOT
      KCMOB2 = KCMOB1 + NCMOT
      KW1    = KCMOB2 + NSTAT2*NCMOT
      LW1    = LWRK   - KW1
      IF (KW1 .GT. LWRK) CALL ERRWRK('WESNAD.BIORAS',-KW1,LWRK)
C... WFOVLP
      KCI1   = KW1
      KCI2   = KCI1   + NCONF
      KCINDX = KCI2   + NCONF
      KW2    = KCINDX + LCINDX
      LW2    = LWRK   - KW2
      IF (KW2 .GT. LWRK) CALL ERRWRK('WESNAD.WFOVLP',-KW2,LWRK)
C
CHJ-860825: Possibility for loop over different DELQ's here
C           e.g.: move left atom, move right atom and in between,
C                 "center-of-mass"
C
C *** CALCULATE GEOMETRY DISPLACED ATOMIC OVERLAP INTEGRALS ***
C
      CALL NADOVL(WRK(KSAOAA),WRK(KSAOAB),WRK(KSAOBB),WRK(KW1),LW1)
C
C
C *** TEST OF INPUT:
C     TEST IF CMOA1,CMOA2,CMOB1,CMOB2 ARE ORTHONORMAL MO'S
C
C
      CALL WESTST(WRK(KSAOAA),WRK(KSAOBB),WRK(KCMOA1),WRK(KW1),LW1)
C
C
C *** CONTROL OF NONADIABATIC COUPLING CALCULATION ***
C
C     Get CI index information
C
      IF (NCONF .GT. 1) THEN
         NOSYM = 0
         CALL GETCIX(WRK(KCINDX),LSYM1,LSYM1,WRK(KW2),LW2,NOSYM)
      END IF
C
      IF (NSTAT2.EQ.0) GOTO 1000
C     if (NSTAT2.eq.0) only adiabatic correction
C
C Obtain first OA1B2 and OA2B1 elements needed in any case ***
C
      CALL HEADER('Cross overlaps for non-adiabatic coupling.',1)
      CALL WESGMO(WRK(KCMOA1),WRK(KCMOA2),WRK(KCMOB1),WRK(KCMOB2),
     *            WRK(KW1),LW1)
C
C     CALL WESGCI(IGEOM1,ISTAT1,NCVEC1,CVEC1,
C    *            IGEOM2,ISTAT2,NCVEC2,CVEC2)
C     OA1B2 = WFOVLP(CVECA,CMOA,CVECB,CMOB,SAO,
C    *               INDXCI,WRK,LWRK()
C
      CALL WESGCI(1,1,NCONF,WRK(KCI1),2,2,NCONF,WRK(KCI2))
      OA1B2 = WFOVLP(WRK(KCI1),WRK(KCMOA1),WRK(KCI2),WRK(KCMOB2),
     *               WRK(KSAOAB),WRK(KCINDX),WRK(KW2),LW2)
C
      CALL WESGCI(1,2,NCONF,WRK(KCI1),2,1,NCONF,WRK(KCI2))
      OA2B1 = WFOVLP(WRK(KCI1),WRK(KCMOA2),WRK(KCI2),WRK(KCMOB1),
     *               WRK(KSAOAB),WRK(KCINDX),WRK(KW2),LW2)
C
C
C
C **********************************************************************
C **********************************************************************
C **********************************************************************
C
      IF (SYMORT .AND. IPRNAD .GE. 1) THEN
         CALL HEADER('RESULTS WITHOUT SYMMETRIC ORTHOGONALIZATION',-1)
         WRITE(LUPRI,'(//A,F12.8)') ' Displacement A to B: ',DISPL
         WRITE(LUPRI,'(/,2(/A,F20.12)/)')
     *      ' 1: <a1 / b2> : ',OA1B2,
     *      ' 2: <a2 / b1> : ',OA2B1
         IF (DISPL .NE. D0) THEN
            COUPL1 = DP5 * (OA1B2 - OA2B1) / DISPL
            COUPL2 = - (OA1B2 + OA2B1) / (DISPL*DISPL)
            WRITE (LUPRI,'(/A,F20.8)')
     *         '  <    WF(1)   /  d/dQ WF(2) > =',COUPL1,
     *         '  < d/dQ WF(1) /  d/dQ WF(2) > =',COUPL2
         END IF
      END IF
C
C   We need four extra overlap elements: OA1B1, OA2B2, OA1A2, OB1B2
C
 1000 CONTINUE
C
C   First <a1 / b1> :
C
      CALL WESGMO(WRK(KCMOA1),WRK(KCMOA2),WRK(KCMOB1),WRK(KCMOB2),
     *            WRK(KW1),LW1)
C
      CALL WESGCI(1,1,NCONF,WRK(KCI1),2,1,NCONF,WRK(KCI2))
      OA1B1 = WFOVLP(WRK(KCI1),WRK(KCMOA1),WRK(KCI2),WRK(KCMOB1),
     *               WRK(KSAOAB),WRK(KCINDX),WRK(KW2),LW2)
C
C   Then <a2 / b2>, <a1 / a2> and <b1 / b2> (if NSTAT2.eq.1):
C
      IF (SYMORT .AND. NSTAT2 .EQ. 1) THEN
         CALL WESGCI(1,2,NCONF,WRK(KCI1),2,2,NCONF,WRK(KCI2))
         OA2B2 = WFOVLP(WRK(KCI1),WRK(KCMOA2),WRK(KCI2),WRK(KCMOB2),
     *                  WRK(KSAOAB),WRK(KCINDX),WRK(KW2),LW2)
C
         CALL WESGMO(WRK(KCMOA1),WRK(KCMOA2),WRK(KCMOB1),WRK(KCMOB2),
     *               WRK(KW1),LW1)
C
         CALL WESGCI(1,1,NCONF,WRK(KCI1),1,2,NCONF,WRK(KCI2))
         OA1A2 = WFOVLP(WRK(KCI1),WRK(KCMOA1),WRK(KCI2),WRK(KCMOA2),
     *                  WRK(KSAOAA),WRK(KCINDX),WRK(KW2),LW2)
C
         CALL WESGCI(2,1,NCONF,WRK(KCI1),2,2,NCONF,WRK(KCI2))
         OB1B2 = WFOVLP(WRK(KCI1),WRK(KCMOB1),WRK(KCI2),WRK(KCMOB2),
     *                  WRK(KSAOBB),WRK(KCINDX),WRK(KW2),LW2)
C
C        Write out overlap elements with many decimals
C
         IF (IPRNAD .GE. 5) WRITE(LUPRI,'(///A//,(A,F20.14))')
     *         ' Calculated wave function overlaps:',
     *         ' OA1B2 = ', OA1B2, ' OA2B1 = ', OA2B1,
     *         ' OA1B1 = ', OA1B1, ' OA2B2 = ', OA2B2,
     *         ' OA1A2 = ', OA1A2, ' OB1B2 = ', OB1B2
      ELSE
         OA1B2 = 1.0D0
         OA2B1 = 1.0D0
         OA2B2 = 1.0D0
         OA1A2 = 1.0D0
         OB1B2 = 1.0D0
         IF (abs(OA1B1) .GT. 1.0D0) IPRNAD = 100
         IF (IPRNAD .GE. 0) WRITE(LUPRI,'(///A//,(A,F20.14))')
     *         ' Calculated wave function overlaps:',
     *         ' OA1B1 = ', OA1B1
         IF (abs(OA1B1) .GT. 1.0D0) THEN
            CALL QUIT('FATAL ERROR, wave function overlap > 1.0 !')
         END IF
      END IF
C
C        Fix phase of B1 and B2 such that OA1B1 .gt. 0 and OA2B2 .gt. 0
C
         IPHFIX = 0
         IF (OA1B1 .LT. D0) THEN
            IPHFIX = 1
            OA1B1 = -OA1B1
            OA2B1 = -OA2B1
            OB1B2 = -OB1B2
         END IF
         IF (OA2B2 .LT. D0) THEN
            IPHFIX = 1
            OA2B2 = -OA2B2
            OA1B2 = -OA1B2
            OB1B2 = -OB1B2
         END IF
         IF (IPRNAD .GE. 0 .AND. IPHFIX .GT. 0) THEN
          IF (SYMORT .AND. NSTAT2.EQ.1) THEN
            WRITE(LUPRI,'(///A//,(A,F20.14))')
     *         ' Calculated wave function overlaps after phase fix:',
     *         ' OA1B2 = ', OA1B2, ' OA2B1 = ', OA2B1,
     *         ' OA1B1 = ', OA1B1, ' OA2B2 = ', OA2B2,
     *         ' OA1A2 = ', OA1A2, ' OB1B2 = ', OB1B2
          ELSE
            WRITE(LUPRI,'(///A//,(A,F20.14))')
     *         ' Calculated wave function overlaps after phase fix:',
     *         ' OA1B1 = ', OA1B1
          END IF
         END IF
C
C        Symmetric orthonormalization:
C        compute C11, C12, C21, C22 factors for A and B geometry
C
      IF (SYMORT .AND. NSTAT2 .EQ. 1) THEN
         X1   = DP5/SQRT(D1 + OA1A2)
         X2   = DP5/SQRT(D1 - OA1A2)
         CA11 = X1 + X2
         CA22 = CA11
         CA12 = X1 - X2
         CA21 = CA12
C
         X1   = DP5/SQRT(D1 + OB1B2)
         X2   = DP5/SQRT(D1 - OB1B2)
         CB11 = X1 + X2
         CB22 = CB11
         CB12 = X1 - X2
         CB21 = CB12
C
         IF (IPRNAD .GT. 0) THEN
            CALL HEADER('Symmetrically orthogonalized wave functions',
     *                  -1)
            WRITE (LUPRI, '(T11,A,SS,F13.8,A,SP,F14.8,A)')
     *         ' /A1''>  =',CA11,'  /A1> ',CA12,'  /A2>',
     *         ' /A2''>  =',CA21,'  /A1> ',CA22,'  /A2>',
     *         ' /B1''>  =',CB11,'  /B1> ',CB12,'  /B2>',
     *         ' /B2''>  =',CB21,'  /B1> ',CB22,'  /B2>'
         END IF
C
         XA1B2 = CA11*CB21*OA1B1 + CA12*CB21*OA2B1
     *         + CA11*CB22*OA1B2 + CA12*CB22*OA2B2
C
         XA2B1 = CA21*CB11*OA1B1 + CA22*CB11*OA2B1
     *         + CA21*CB12*OA1B2 + CA22*CB12*OA2B2
C
         XA1B1 = CA11*CB11*OA1B1 + CA12*CB11*OA2B1
     *         + CA11*CB12*OA1B2 + CA12*CB12*OA2B2
C
         XA2B2 = CA21*CB21*OA1B1 + CA12*CB22*OA1B2
     *         + CA22*CB21*OA2B1 + CA22*CB22*OA2B2
C
C        Write out overlap elements with many decimals
C
         IF (IPRNAD .GE. 5) THEN
            WRITE(LUPRI,'(///2A//,(A,F20.14))')
     *         ' Calculated wave function overlaps',
     *         ' after orthogonalization:',
     *         ' XA1B2 = ', XA1B2, ' XA2B1 = ', XA2B1,
     *         ' XA1B1 = ', XA1B1, ' XA2B2 = ', XA2B2
         END IF
C
         SYM1  = XA1B2 - XA2B1
         SYM2  = XA1B2 + XA2B1
C
         BORN11 = XA1B1
         BORN22 = XA2B2
C
         RESA1 = D1 - XA1B1*XA1B1 - XA1B2*XA1B2
         RESB1 = D1 - XA1B1*XA1B1 - XA2B1*XA2B1
         RESA2 = D1 - XA2B2*XA2B2 - XA2B1*XA2B1
         RESB2 = D1 - XA2B2*XA2B2 - XA1B2*XA1B2
C
         IF (RESA1 .LT. D0 .OR. RESA2 .LT. D0 .OR.
     *       RESB1 .LT. D0 .OR. RESB2 .LT. D0) THEN
             WRITE(LUPRI,'(///2A/A//,(A,F20.14))')
     *            ' *** WARNING *** phase error in wave function',
     *            ' as a residual squared is negative;',
     *            '                 the residuals squared:',
     *            ' RESA1 = ', RESA1, ' RESA2 = ', RESA2,
     *            ' RESB1 = ', RESB1, ' RESB2 = ', RESB2
         ELSE
            RESA1 = SQRT(RESA1)
            RESB1 = SQRT(RESB1)
            RESA2 = SQRT(RESA2)
            RESB2 = SQRT(RESB2)
         END IF
      ELSE
         BORN11 = OA1B1
C
C End of symmetric orthogonalization part
C
      END IF
C
C ********************************************************************
C
C Final part: Coupling element and ouput section
C
C
C
      CALL HEADER('<<<<< FINAL RESULTS FROM WESNAD >>>>>',-1)
      CALL TIMDAT
      IF ( NSTAT2.EQ.1 .AND. (.NOT. SYMORT .OR. IPRNAD .GE. 2) ) THEN
         CALL HEADER('Results without symmetric orthogonalization',-1)
         WRITE(LUPRI,'(//A,F12.8)') ' Displacement A to B: ',DISPL
         WRITE(LUPRI,'(/2(/A,F20.12)/)')
     *      ' 1: <a1 / b2> : ',OA1B2,
     *      ' 2: <a2 / b1> : ',OA2B1
         IF (DISPL .NE. D0) THEN
            COUPL1 = DP5 * (OA1B2 - OA2B1) / DISPL
            COUPL2 = - (OA1B2 + OA2B1) / (DISPL*DISPL)
            WRITE (LUPRI,'(/A,F20.8)')
     *         '  <    WF(1)   /  d/dQ WF(2) > =',COUPL1,
     *         '  < d/dQ WF(1) /  d/dQ WF(2) > =',COUPL2
         END IF
      END IF
      IF (SYMORT .AND. NSTAT2.EQ.1) THEN
         CALL HEADER('RESULTS WITH SYMMETRIC ORTHOGONALIZATION',-1)
         WRITE(LUPRI,'(/A/)') ' Wave function overlaps:'
         WRITE(LUPRI,'(/T22,A,F9.6,A/T22,A/T22,A,F9.6,A/T22,A/
     *             /T22,F11.6,A,F9.6/
     *             /T22,A/T22,A,F9.6,A/T22,A/T22,A,F9.6,A/)')
     *   '     /A2> ----',OA2B2,' ---- /B2>',
     *   '      |   \\                /  |',
     *   '      |  ',OA2B1,'       /    |',
     *   '      |       \\        /      |',
     *        OA1A2,'       *       ',OB1B2,
     *   '      |       /        \\      |',
     *   '      |  ',OA1B2,'       \\    |',
     *   '      |   /                \\  |',
     *   '     /A1> ----',OA1B1,' ---- /B1>'
         WRITE(LUPRI,'(/A/)') ' Overlap between '//
     *      'symmetrically orthogonalized wave functions:'
         WRITE(LUPRI,'(/T22,A,F9.6,A/T22,A/T22,A,F9.6,A/T22,A/
     *             /T22,F11.6,A,F9.6/
     *             /T21,A/T22,A,F9.6,A/T22,A/T22,A,F9.6,A/)')
     *   '    /A2''> ----',XA2B2,' ---- /B2''>',
     *   '      |   \\                /  |',
     *   '      |  ',XA2B1,'       /    |',
     *   '      |       \\        /      |',
     *         D0  ,'       *       ', D0,
     *   '      |       /        \\      |',
     *   '      |  ',XA1B2,'       \\    |',
     *   '      |   /                \\  |',
     *   '    /A1''> ----',XA1B1,' ---- /B1''>'
         IF (IPRNAD .GE. 1) WRITE(LUPRI,'(/,3(/A),4(/A,F20.8)/)')
     *      ' Coupling to other states',
     *      ' ( i.e. overlap to other states but /1> and /2> )',
     *      ' ------------------------------------------------',
     *      ' Residual coupling to /A1''> =',RESA1,
     *      ' Residual coupling to /A2''> =',RESA2,
     *      ' Residual coupling to /B1''> =',RESB1,
     *      ' Residual coupling to /B2''> =',RESB2
         WRITE(LUPRI,'(//A,F12.8)') ' Displacement A to B: ',DISPL
         IF (DISPL .NE. D0) THEN
            COUPL1 = DP5 * SYM1 / DISPL
            COUPL2 = (-SYM2) / (DISPL*DISPL)
            BORN11 = 2.00D0*(D1 - BORN11) / (DISPL*DISPL)
            BORN22 = 2.00D0*(D1 - BORN22) / (DISPL*DISPL)
            WRITE (LUPRI,'(/A,F20.8)')
     *         '  <    WF(1)   /  d/dQ WF(2) > =',COUPL1,
     *         '  < d/dQ WF(1) /  d/dQ WF(2) > =',COUPL2,
     *         '  < d/dQ WF(1) /  d/dQ WF(1) > =',BORN11,
     *         '  < d/dQ WF(2) /  d/dQ WF(2) > =',BORN22
         END IF
      ELSE IF (NSTAT2 .EQ. 0) THEN
         CALL HEADER('ADIABATIC CORRECTION',-1)
         WRITE(LUPRI,'(//A,F12.8)') ' Displacement A to B: ',DISPL
         IF (DISPL .NE. D0) THEN
            BORN11 = 2.0D0*(D1 - BORN11) / (DISPL*DISPL)
            WRITE (LUPRI,'(/A,F20.8,A)')
     *         '  < d/dQ WF(1) /  d/dQ WF(1) > =',BORN11,' au m_e'
            BORN11 = BORN11*XTKAYS
            WRITE (LUPRI,'(/A,F20.3,A)')
     *         '  < d/dQ WF(1) /  d/dQ WF(1) > =',BORN11,' cm-1 m_e'
            BORN11 = BORN11/XFAMU
            WRITE (LUPRI,'(/A,F20.3,A)')
     *         '  < d/dQ WF(1) /  d/dQ WF(1) > =',BORN11,' cm-1 u'
            BORN11A = BORN11/(2.0D0*NUCMASS(1))
            BORN11B = BORN11/(2.0D0*NUCMASS(2))
            BORN11  = BORN11A + BORN11B
            WRITE (LUPRI,'(/A,F20.3,A)')
     *         '  < d/dQ WF(1) /  d/dQ WF(1) >/2m_A =',BORN11A,' cm-1',
     *         '  < d/dQ WF(1) /  d/dQ WF(1) >/2m_B =',BORN11B,' cm-1',
     *         '  < d/dQ WF(1) /  d/dQ WF(1) >/2mu  =',BORN11 ,' cm-1'
         END IF
      END IF
C
C     END OF WESTA.WESNAD.
C
      RETURN
      END
#  /* Deck wesgmo */
      SUBROUTINE WESGMO(CMOA1,CMOA2,CMOB1,CMOB2,WRK,LWRK)
C
C Get A and B MO's
C
#include "implicit.h"
      DIMENSION CMOA1(*),CMOA2(*),CMOB1(*),CMOB2(*), WRK(*)
C
C Used from common blocks:
C  WESINF: IRDMOW
C  WESTAP: LUIT1V(2,2)
C
#include "wesinf.h"
#include "westap.h"
C
C
C
      CALL READMW(IRDMOW,LUIT1V(1,1),CMOA1)
      CALL READMW(IRDMOW,LUIT1V(2,1),CMOB1)
      IF (NSTAT2 .EQ. 1) THEN
         CALL READMW(IRDMOW,LUIT1V(1,2),CMOA2)
         CALL READMW(IRDMOW,LUIT1V(2,2),CMOB2)
      END IF
C
C
      RETURN
      END
      SUBROUTINE NADOVL(SAOAA,SAOAB,SAOBB,WRK,LWRK)
C
C   (c) Hans Joergen Aa. Jensen, Sep. 2004
C
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION SAOAA(NBAST,*), SAOAB(NBAST,*), SAOBB(NBAST,*),
     &          WRK(N2BASX,*)
      CHARACTER*8 SQHLABEL
C
C  WESINF : DISPL, ?
C  INFORB : NBAST, ?
C  WESTAP : LUPROPV(2)
C
#include "wesinf.h"
#include "inforb.h"
#include "westap.h"

C
      IF (LWRK .LT. 3*N2BASX) CALL QUIT('Too little memory NADOVL')

C     Overlap AA
      CALL WESPROP('OVERLAP ',LUPROPV(1),WRK,IANTSYM)
      IF (IANTSYM .NE. 1) CALL QUIT( 'NADOVL AA' )
      CALL DSPTSI(NBAST,WRK,SAOAA)

C     Overlap <A|B> = <A|A> - DISPL * <A|-dA/dQ>
      CALL DCOPY(N2BASX,SAOAA,1,SAOAB,1)
      K = 0
      DO JATOM = 1,NMAXD
      DO ICOOR = 1,3
         K = K + 1
         XDELQ = DELQ(ICOOR,JATOM)
      IF (XDELQ .NE. 0.0D0) THEN
         WRITE(SQHLABEL,'(A5,I3.3)') 'SQHDR',K
         write (lupri,'(2A)')'SQHDOR  label: ', SQHLABEL
         CALL WESPROP(SQHLABEL,LUPROPV(1),WRK,IANTSYM)
         IF (IANTSYM .NE. 0) CALL QUIT( 'NADOVL AB' )
         IF (IPRNAD .GE. 8) THEN
            WRITE (LUPRI,'(//2A)') '<A|-dA/dQ> matrix ',SQHLABEL
            CALL OUTPUT(WRK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         CALL DAXPY(N2BASX,-XDELQ,WRK,1,SAOAB,1)
         FAC = 0.5D0*XDELQ*XDELQ
#ifdef oldd2o
         CALL DCOPY(N2BASX,WRK,1,WRK(1,2),1)
         CALL DCOPY(N2BASX,SAOAA,1,WRK(1,3),1)
         CALL DGESOL(NBAST,NBAST,NBAST,NBAST,
     &               WRK(1,3),WRK(1,2),WRK(1,4),INFO)
C                 A       ,B       ,KPVT    ,INFO)
C         B := A(inv) B
C
         IF (INFO .ne. 0) CALL QUIT( 'NADOVL INFO .ne. 0' )
!         CALL XMPAB(FAC,
!     1              WRK,NBAST,NBAST,NBAST,NBAST,
!     2              WRK(1,2),NBAST,NBAST,NBAST,NBAST,
!     3              SAOAB,NBAST,NBAST)
          CALL MPAB(
      1             WRK,NBAST,NBAST,NBAST,NBAST,
      2             WRK(1,2),NBAST,NBAST,NBAST,NBAST,
      3             WRK(1,3),NBAST,NBAST)
         write (lupri,*) '2. derivative as sqhdo squared'
         call output(wrk(1,3),1,nbast,1,nbast,nbast,nbast,1,lupri)
         write (lupri,*) '2. derivative from sqhd2o'
         call output(wrk(1,2),1,nbast,1,nbast,nbast,nbast,1,lupri)
         call daxpy (n2basx,-1.0d0,wrk(1,2),1,wrk(1,3),1)
         write (lupri,*) '2. derivative difference'
         call output(wrk(1,3),1,nbast,1,nbast,nbast,nbast,1,lupri)
#else
         WRITE(SQHLABEL,'(A5,I3.3)') 'HD2OR',K
         write (lupri,'(2A)')'SQHD2OR label: ', SQHLABEL
         CALL WESPROP(SQHLABEL,LUPROPV(1),WRK(1,2),IANTSYM)
         IF (IPRNAD .GE. 8) THEN
            WRITE (LUPRI,'(//2A)') '<A|-d2A/dQ2> matrix ',SQHLABEL
            CALL OUTPUT(WRK(1,2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (IANTSYM .NE. 0) CALL QUIT( 'NADOVL2 AB' )
         CALL DAXPY(N2BASX,-FAC,WRK(1,2),1,SAOAB,1)
#endif
      END IF
      END DO
      END DO

C     Overlap BB
      CALL WESPROP('OVERLAP ',LUPROPV(2),WRK,IANTSYM)
      IF (IANTSYM .NE. 1) CALL QUIT( 'NADOVL BB' )
      CALL DSPTSI(NBAST,WRK,SAOBB)

C     Overlap BA, for test
!     CALL WESPROP(SQHLABEL,LUPROPV(2),WRK,IANTSYM)
!     IF (IANTSYM .NE. 0) CALL QUIT( 'NADOVL BA' )
!     IF (IPRNAD .GE. 8) THEN
!        WRITE (LUPRI,'(//A)') '<B|-dB/dQ> matrix ',SQHLABEL
!        CALL OUTPUT(WRK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
!     END IF
!     CALL DSCAL(N2BASX,DISPL,WRK,1)
!     DO I = 1, NBAST
!        VNRM = DDOT(NBAST,WRK(1,I),1,WRK(1,I),1)
!        WRK(I,I) = SQRT(SAOBB(I,I)-VNRM) - SAOBB(I,I)
!     END DO
!     CALL DAXPY(N2BASX,1.0D0,SAOBB,1,WRK,1)
!     CALL DGETRN(WRK,NBAST,NBAST)
!     CALL DAXPY(N2BASX,-1.0D0,SAOAB,1,WRK,1)
!     ABBANORM = DNRM2(N2BASX,WRK,1)
      ABBANORM = 0.0D0

!     IF (IPRNAD .GE. 8 .OR. ABBANORM .GT. 1.0D-5) THEN
      IF (IPRNAD .GE. 8) THEN
         WRITE (LUPRI,'(//A)') ' SAOAA: AO (A-A square) overlap matrix'
         CALL OUTPUT(SAOAA,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE (LUPRI,'(//A)') ' SAOAB: AO (A-B square) overlap matrix'
         CALL OUTPUT(SAOAB,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE (LUPRI,'(//A)') ' SAOBB: AO (B-B square) overlap matrix'
         CALL OUTPUT(SAOBB,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
!        WRITE (LUPRI,'(//A,1P,D10.2)')
!    &      ' Test: SAOAB - SAOBA(transp), matrix norm is',ABBANORM
!        CALL OUTPUT(WRK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
!     ELSE
!        WRITE (LUPRI,'(//A,1P,D10.2)')
!    &      ' Test: SAOAB - SAOBA(transp), matrix norm is',ABBANORM
      END IF

      RETURN
      END
