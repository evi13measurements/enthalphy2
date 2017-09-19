C
C   File: wesgp.F
C
C  /* Deck readmw */
      SUBROUTINE READMW(IRDMO,LUMO,CMO)
C
C 17-May-1995 hjaaj
C
#include "implicit.h"
      DIMENSION CMO(*)
      PARAMETER (LWORK = 1)
C
C  PRIUNIT : LUPRI
C  INFTAP : LUIT1, LUSIFC
C
#include "priunit.h"
#include "inftap.h"
C
      IF (IRDMO .NE. 7 .AND. IRDMO .NE. 9 .AND. IRDMO .NE. 10) THEN
         WRITE (LUPRI,'(//A,I6)')
     &   'READMW ERROR: specified IRDMO not implemented, IRDMO =',IRDMO
         STOP 'READMW ERROR: specified IRDMO not implemented.'
      END IF
      LUSIFC = LUMO
      LUIT1  = LUMO
      CALL READMO(CMO,IRDMO,WORK,LWORK)
      RETURN
      END
C  /* Deck wesgci */
      SUBROUTINE WESGCI(IGEOM1,ISTAT1,NCVEC1,CVEC1,
     *                  IGEOM2,ISTAT2,NCVEC2,CVEC2)
C
C  Purpose: get CI-vectors
C           mc energies are read (in EMCWES(2,2)
C
#include "implicit.h"
      DIMENSION CVEC1(2),CVEC2(2)
C
      PARAMETER ( D1 = 1.0D0 )
#include "thrzer.h"
C
C Used from common blocks:
C  PRIUNIT: LUPRI
C  WESINF: EMCWES(2,2)
C  WESTAP: LUIT1V(2,2)
C
#include "priunit.h"
#include "wesinf.h"
#include "westap.h"
#include "inftap.h"
C
C Externals:
C
      LOGICAL FNDLAB
C
      LUIT1 = LUIT1V(IGEOM1,ISTAT1)
      CALL GETREF(CVEC1,NCVEC1)
      CVNRM = DNRM2(NCVEC1,CVEC1,1)
      IF (ABS(D1-CVNRM) .GT. 1.0D-10) THEN
         WRITE (LUPRI,*) 'WESGCI ERROR: 1 - norm of CVEC1 = ',D1-CVNRM
         CALL QUIT('WESGCI ERROR, CVEC1 not normalized')
      ELSE IF (ABS(D1-CVNRM) .GT. THRZER) THEN
         WRITE (LUPRI,*) 'WESGCI WARNING: 1 - norm of CVEC1 = ',D1-CVNRM
      END IF
C...energy
      REWIND (LUIT1)
      READ (LUIT1)
      READ (LUIT1) DUM, DUM, DUM, EMCWES(IGEOM1,ISTAT1)
C
      IF (ISTAT2 .GT. 0) THEN
         LUIT1 = LUIT1V(IGEOM2,ISTAT2)
         CALL GETREF(CVEC2,NCVEC2)
         CVNRM = DNRM2(NCVEC2,CVEC2,1)
         IF (ABS(D1-CVNRM) .GT. 1.0D-10) THEN
            WRITE (LUPRI,*)
     &      'WESGCI ERROR: 1 - norm of CVEC2 = ',D1-CVNRM
            CALL QUIT('WESGCI ERROR, CVEC2 not normalized')
         ELSE IF (ABS(D1-CVNRM) .GT. THRZER) THEN
            WRITE (LUPRI,*)
     &      'WESGCI WARNING: 1 - norm of CVEC2 = ',D1-CVNRM
         END IF
C...energy
         REWIND (LUIT1)
         READ (LUIT1)
         READ (LUIT1) DUM, DUM, DUM, EMCWES(IGEOM2,ISTAT2)
      END IF
C
      RETURN
      END
C  /* Deck symmo */
      SUBROUTINE SYMMO(IPAR,CMO,VC)
C
C  WESTA.SYMMO /ha
C
C  (Revision 17-May-1986)
C
C Purpose:
C    To transform MO:s from symmetry-adapted (CMO) to
C    non-symmetry-adapted (VC) basis, where VC runs over
C    NBAST basis functions. CMO sorted in symmetry first order
C    (normal order), VC is sorted inactive-active-secondary order
C    i.e. "ISW" order.
C
C    IPAR = 1: Transform CMO to VC
C    IPAR = 2: Transform VC to CMO
C
C  (Revision 17-Sep-1989)
C
C   Using ISWR and ISXR to assure correct trasnformation also
C   Considering RAS division of active space
C
#include "implicit.h"
      DIMENSION CMO(2),VC(NBAST,2)
C
C  PRIUNIT : LUPRI
C  INFORB : NSYM,NORBT,NBAST,NORB(*),NBAS(*),IORB(*),IBAS(*)
C  INFIND : ISW(*)
C
#include "priunit.h"
#include "maxash.h"
#include "maxorb.h"
#include "wcborb.h"
#include "inforb.h"
#include "infind.h"
C
C
      IF (IPAR.NE.1) GO TO 400
C
C---- Transform CMO to VC
C
      CALL DZERO(VC,NBAST*NORBT)
C
      JCMO = 0
      DO 390 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI .EQ. 0) GO TO 390
         NBASI = NBAS(ISYM)
         IORBI = IORB(ISYM)
         IBASI = IBAS(ISYM)
         DO 380 IMO = IORBI+1,IORBI+NORBI
            IMONEW = ISWR(IMO)
            DO 370 IAO = 1,NBASI
               VC(IBASI+IAO,IMONEW) = CMO(JCMO+IAO)
  370       CONTINUE
            JCMO = JCMO + NBASI
  380    CONTINUE
  390 CONTINUE
C
      GO TO 800
C
C---- Transform VC to CMO
C
  400 CONTINUE
C
      IF (IPAR.NE.2) THEN
         WRITE(LUPRI,'(/A,I5)') ' SYMMO Parameter in error, IPAR =',IPAR
         STOP 'WESTA.SYMMO: invalid control parameter'
      ENDIF
C
      CALL DZERO(CMO,NCMOT)
C
      JCMO = 0
      DO 490 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI .EQ. 0) GO TO 490
         NBASI = NBAS(ISYM)
         IORBI = IORB(ISYM)
         IBASI = IBAS(ISYM)
         DO 480 IMO = IORBI+1,IORBI+NORBI
            IMONEW = ISWR(IMO)
            DO 470 IAO = 1,NBASI
               CMO(JCMO+IAO) = VC(IBASI+IAO,IMONEW)
  470       CONTINUE
            JCMO = JCMO + NBASI
  480    CONTINUE
  490 CONTINUE
C
C
  800 CONTINUE
      RETURN
      END
C  /* Deck wrspec */
      SUBROUTINE WRSPEC(CVEC,CMO,ESPEC)
C
C  21-MAY-1986
C
C  Purpose: To write a set of CI- and orbitals vectors (read from
C           LUSPEC) consequtively on a big temporary file (LUSPCW).
C           To retrieve spectral energies ESPEC.
C
#include "implicit.h"
      DIMENSION CVEC(2),CMO(2),ESPEC(2)
C
      PARAMETER ( D0 = 0.0D0 )
C
C
C Used from common blocks:
C  PRIUNIT : LUPRI
C  INFORB : NCMOT..
C  WESINF : NCONF2,NSTAT2.
C  WESTAP : MSTAT2, LUSPEC(MSTAT2), LUSPCW
C
#include "priunit.h"
#include "inforb.h"
#include "wesinf.h"
#include "westap.h"
C
C Externals:
C
      LOGICAL FNDLAB
C
C
      CHARACTER*8 RTNLBL(2)
      CHARACTER*7 FILENM(MSTAT2)
      DATA FILENM /'fort.61','fort.62',
     *             'fort.63','fort.64',
     *             'fort.65','fort.66',
     *             'fort.67','fort.68',
     *             'fort.69','fort.70',
     *             'fort.71','fort.72',
     *             'fort.73','fort.74',
     *             'fort.75','fort.76',
     *             'fort.77','fort.78',
     *             'fort.79','fort.80',
     *             'fort.81','fort.82',
     *             'fort.83','fort.84',
     *             'fort.85','fort.86',
     *             'fort.87','fort.88',
     *             'fort.89','fort.90',
     *             'fort.91','fort.92',
     *             'fort.93','fort.94',
     *             'fort.95','fort.96',
     *             'fort.97','fort.98',
     *             'fort.99','fort100'/
C
C     open LUSPCW
C
      OPEN(LUSPCW,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *     FILE='fort.51')
C
C    start reading vectors form LUSPEC. We assume that "NEWORB" is
C    positioned after "STARTVEC" or "RESTART " not available.
C
      REWIND LUSPCW
      DO 100 ISTAT2 = 1,NSTAT2
C
         LUITX = LUSPEC(ISTAT2)
         OPEN(LUITX,STATUS='OLD',FORM='UNFORMATTED',FILE=FILENM(ISTAT2))
         REWIND LUITX
         CALL MOLLB2 ('STARTVEC',RTNLBL,LUITX,LUPRI)
         READ(RTNLBL(1),'(2I4)') NRS, JSTA
         IF (ABS(JSTA) .NE. ICIORD(ISTAT2)) THEN
            WRITE(6,*) 'WRSPEC error: incorrect MCSCF state.'
            WRITE(6,*) 'ISTAT2,ICIORD(ISTAT2)',ISTAT2,ICIORD(ISTAT2)
            WRITE(6,*) 'NRS   ,"istate"      ',NRS,JSTA
            STOP 'WRSPEC error: incorrect MCSCF state'
         END IF
         IF (JSTA .GT. 0) THEN
            DO 200 ISKIP = 1, (ICIORD(ISTAT2) - 1)
               READ(LUITX)
  200       CONTINUE
         END IF
         CALL READT(LUITX,NCONF2,CVEC)
         REWIND LUITX
         IF ( FNDLAB('RESTART ',LUITX) ) THEN
            READ (LUITX) ESPEC(ISTAT2)
         ELSE
            REWIND LUITX
            ESPEC(ISTAT2) = D0
         END IF
         CALL MOLLAB('NEWORB  ',LUITX,LUPRI)
         CALL READT(LUITX,NCMOT,CMO)
C
         CALL WRITT(LUSPCW,NCONF2,CVEC)
         CALL WRITT(LUSPCW,NCMOT,CMO)
C
         CLOSE(UNIT=LUITX,STATUS='KEEP')
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck westst */
      SUBROUTINE WESTST(SAOAA,SAOBB,CMO,WRK,LWRK)
C
C 3-Feb-1985 hjaaj
C
C Input test, get MO's and check if they are orthonormal
C
#include "implicit.h"
      DIMENSION SAOAA(*), SAOBB(*), CMO(*), WRK(*)
C
      PARAMETER ( D1 = 1.0D0, SMOMAX = 1.D-8, SMOWRN = 1.D-12 )
C
C Used from common blocks:
C  WESINF: IRDMOW,IPRDEF,NGEOM
C  WESTAP: LUIT1V(2,2)
C  INFORB: NORBT,NBAST,N2ORBX
C  PRIUNIT: LUPRI
C
#include "priunit.h"
#include "wesinf.h"
#include "westap.h"
#include "inforb.h"
C
      INPERR = 0
C
      KVC    = 1
      KSMO   = KVC   + NORBT*NBAST
      KUTHU  = KSMO  + N2ORBX
      KWRK1  = KUTHU + NBAST
      LWRK1  = KWRK1 + N2BASX - 1
      IF (LWRK1 .GT. LWRK) CALL ERRWRK('WESTST',LWRK1,LWRK)
C
C     get and check CMOA1
C
      CALL READMW(IRDMOW,LUIT1V(1,1),CMO)
      CALL SYMMO(1,CMO,WRK(KVC))
C
C     If NGEOM = 1  CALL UTHU since overlap matrix then is packed
C
      IF (NGEOM .EQ. 1) THEN
         CALL UTHU(SAOAA,WRK(KSMO),WRK(KVC),WRK(KUTHU),NBAST,NORBT)
         II = KSMO - 1
         DO 1101 I = 1,NORBT
            II = II + I
            WRK(II) = WRK(II) - D1
 1101    CONTINUE
         IMX = IDAMAX(NNORBX,WRK(KSMO),1)
         JMX = KSMO-1 + IMX
         WMX = ABS(WRK(JMX))
C
      ELSE
C
         CALL DGEMM('N','N',NBAST,NORBT,NBAST,1.D0,
     &              SAOAA,     NBAST,
     &              WRK(KVC),  NBAST,0.D0,
     &              WRK(KWRK1),NBAST)
         CALL DGEMM('T','N',NORBT,NORBT,NBAST,1.D0,
     &              WRK(KVC),  NBAST,
     &              WRK(KWRK1),NBAST,0.D0,
     &              WRK(KSMO), NORBT)
C
         II = KSMO
         DO 1100 I = 1,NORBT
            WRK(II) = WRK(II) - D1
            II = II + NORBT + 1
 1100    CONTINUE
         IMX = IDAMAX(N2ORBX,WRK(KSMO),1)
         JMX = KSMO-1 + IMX
         WMX = ABS(WRK(JMX))
C
      ENDIF
C
      JPR = IPRDEF
      IF (JPR .GE. 4) WRITE (LUPRI,'(//A,1PD10.2)')
     *      ' Largest deviation from orthonormality, A1 orbitals:',
     *      WRK(JMX)
      IF (WMX .GT. SMOMAX) THEN
         WRITE (LUPRI,'(//A/A,I5,A,1PD10.2)')
     *      ' ERROR, A1 molecular orbitals not orthonormal',
     *      '        largest element SMO(',IMX,') =',WRK(JMX)
         INPERR = INPERR + 1
         JPR = 100
      ELSE IF (WMX .GT. SMOWRN) THEN
         WRITE (LUPRI,'(//A/A,I5,A,1PD10.2)')
     *      ' WARNING, A1 molecular orbitals not quite orthonormal',
     *      '          largest element SMO(',IMX,') =',WRK(JMX)
      END IF
      IF (JPR .GT. 10) THEN
         WRITE (LUPRI,'(//A/)') ' SMO - 1 for A1 molecular orbitals:'
         IF (NGEOM .EQ. 1) THEN
            CALL OUTPAK(WRK(KSMO),NORBT,1,LUPRI)
         ELSE
            CALL OUTPUT(WRK(KSMO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         ENDIF
      END IF
C
C     get and check CMOA2  (If WESNAD, i.e. NGEOM .gt. 1)
C
      IF (NGEOM .EQ. 1) GO TO 8000
C
      IF (NSTAT2 .EQ. 1) THEN
      CALL READMW(IRDMOW,LUIT1V(1,2),CMO)
      CALL SYMMO(1,CMO,WRK(KVC))
C
C
         CALL DGEMM('N','N',NBAST,NORBT,NBAST,1.D0,
     &              SAOAA,     NBAST,
     &              WRK(KVC),  NBAST,0.D0,
     &              WRK(KWRK1),NBAST)
         CALL DGEMM('T','N',NORBT,NORBT,NBAST,1.D0,
     &              WRK(KVC),  NBAST,
     &              WRK(KWRK1),NBAST,0.D0,
     &              WRK(KSMO), NORBT)
      II = KSMO
      DO I = 1,NORBT
         WRK(II) = WRK(II) - D1
         II = II + NORBT + 1
      END DO
      IMX = IDAMAX(N2ORBX,WRK(KSMO),1)
      JMX = KSMO-1 + IMX
      WMX = ABS(WRK(JMX))
C
C
      JPR = IPRDEF
      IF (JPR .GE. 5) WRITE (LUPRI,'(//A,1PD10.2)')
     *      ' Largest deviation from orthonormality, A2 orbitals:',
     *      WRK(JMX)
      IF (WMX .GT. SMOMAX) THEN
         WRITE (LUPRI,'(//A/A,I5,A,1PD10.2)')
     *      ' ERROR, A2 molecular orbitals not orthonormal',
     *      '        largest element SMO(',IMX,') =',WRK(JMX)
         INPERR = INPERR + 1
         JPR = 100
      ELSE IF (WMX .GT. SMOWRN) THEN
         WRITE (LUPRI,'(//A/A,I5,A,1PD10.2)')
     *      ' WARNING, A2 molecular orbitals not quite orthonormal',
     *      '          largest element SMO(',IMX,') =',WRK(JMX)
      END IF
      IF (JPR .GT. 10) THEN
         WRITE (LUPRI,'(//A/)') ' SMO - 1 for A2 molecular orbitals:'
         IF (NGEOM .EQ. 1) THEN
            CALL OUTPAK(WRK(KSMO),NORBT,1,LUPRI)
         ELSE
            CALL OUTPUT(WRK(KSMO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         ENDIF
      END IF
      END IF  ! (NSTAT2 .EQ. 1)
C
C
C
C
C     get and check CMOB1
C
      CALL READMW(IRDMOW,LUIT1V(2,1),CMO)
      CALL SYMMO(1,CMO,WRK(KVC))
         CALL DGEMM('N','N',NBAST,NORBT,NBAST,1.D0,
     &              SAOBB,     NBAST,
     &              WRK(KVC),  NBAST,0.D0,
     &              WRK(KWRK1),NBAST)
         CALL DGEMM('T','N',NORBT,NORBT,NBAST,1.D0,
     &              WRK(KVC),  NBAST,
     &              WRK(KWRK1),NBAST,0.D0,
     &              WRK(KSMO), NORBT)
      II = KSMO
      DO I = 1,NORBT
         WRK(II) = WRK(II) - D1
         II = II + (NORBT + 1)
      END DO
      IMX = IDAMAX(N2ORBX,WRK(KSMO),1)
      JMX = KSMO-1 + IMX
      WMX = ABS(WRK(JMX))
      JPR = IPRDEF
      IF (JPR .GE. 5) WRITE (LUPRI,'(//A,1PD10.2)')
     *      ' Largest deviation from orthonormality, B1 orbitals:',
     *      WRK(JMX)
      IF (WMX .GT. SMOMAX) THEN
         WRITE (LUPRI,'(//A/A,I5,A,1PD10.2)')
     *      ' ERROR, B1 molecular orbitals not orthonormal',
     *      '        largest element SMO(',IMX,') =',WRK(JMX)
         INPERR = INPERR + 1
         JPR = 100
      ELSE IF (WMX .GT. SMOWRN) THEN
         WRITE (LUPRI,'(//A/A,I5,A,1PD10.2)')
     *      ' WARNING, B1 molecular orbitals not quite orthonormal',
     *      '          largest element SMO(',IMX,') =',WRK(JMX)
      END IF
      IF (JPR .GT. 10) THEN
         WRITE (LUPRI,'(//A/)') ' SMO - 1 for B1 molecular orbitals:'
         CALL OUTPUT(WRK(KSMO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     get and check CMOB2
C
      IF (NSTAT2 .EQ. 1) THEN
      CALL READMW(IRDMOW,LUIT1V(2,2),CMO)
      CALL SYMMO(1,CMO,WRK(KVC))
         CALL DGEMM('N','N',NBAST,NORBT,NBAST,1.D0,
     &              SAOBB,     NBAST,
     &              WRK(KVC),  NBAST,0.D0,
     &              WRK(KWRK1),NBAST)
         CALL DGEMM('T','N',NORBT,NORBT,NBAST,1.D0,
     &              WRK(KVC),  NBAST,
     &              WRK(KWRK1),NBAST,0.D0,
     &              WRK(KSMO), NORBT)
      II = KSMO
      DO 4100 I = 1,NORBT
         WRK(II) = WRK(II) - D1
         II = II + (NORBT + 1)
 4100 CONTINUE
      IMX = IDAMAX(N2ORBX,WRK(KSMO),1)
      JMX = KSMO-1 + IMX
      WMX = ABS(WRK(JMX))
      JPR = IPRDEF
      IF (JPR .GE. 5) WRITE (LUPRI,'(//A,1PD10.2)')
     *      ' Largest deviation from orthonormality, B2 orbitals:',
     *      WRK(JMX)
      IF (WMX .GT. SMOMAX) THEN
         WRITE (LUPRI,'(//A/A,I5,A,1PD10.2)')
     *      ' ERROR, B2 molecular orbitals not orthonormal',
     *      '        largest element SMO(',IMX,') =',WRK(JMX)
         INPERR = INPERR + 1
         JPR = 100
      ELSE IF (WMX .GT. SMOWRN) THEN
         WRITE (LUPRI,'(//A/A,I5,A,1PD10.2)')
     *      ' WARNING, B2 molecular orbitals not quite orthonormal',
     *      '          largest element SMO(',IMX,') =',WRK(JMX)
      END IF
      IF (JPR .GT. 10) THEN
         WRITE (LUPRI,'(//A/)') ' SMO - 1 for B2 molecular orbitals:'
         CALL OUTPUT(WRK(KSMO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      END IF ! (NSTAT2 .EQ. 1)
C
C     end of WESTST
C
 8000 CONTINUE
      IF (INPERR .GT. 0) THEN
         WRITE (LUPRI,'(//A)') ' INPUT ERRORS FOUND BY WESTST.'
         STOP 'INPUT ERRORS FOUND BY WESTST.'
      END IF
      RETURN
      END
C  /* Deck timdat */
      SUBROUTINE TIMDAT
#include "implicit.h"
#include "priunit.h"
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      INTEGER CDATE, CTIME
      CDATE = DATE()
      CTIME = CLOCK()
      WRITE(LUPRI,100) CDATE, CTIME
  100 FORMAT(/,'  Date and time (CRAY) : ',A8,2X,A8)
#endif
#if defined (SYS_AIX)
C     AIX XL FORTRAN version 2.3+
      CHARACTER*(24) fdate
      WRITE (LUPRI,100) fdate()
  100 FORMAT(/,' Date and time (IBM-AIX): ',A)
#endif
#if defined (SYS_HPUX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (HP-UX)  : ',A)
#endif
#if defined (SYS_DEC)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (DEC)  : ',A)
#endif
#if defined (SYS_IRIX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (IRIX)   : ',A)
#endif
#if defined (SYS_LINUX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Linux) : ',A)
#endif
#if defined (SYS_DARWIN)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Darwin) : ',A)
#endif
#if defined (SYS_SX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Nec-SX) : ',A)
#endif
#if defined (SYS_SUN)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (SUN)   : ',A)
#endif
      RETURN
      END
      SUBROUTINE WESPROP(WORD,LUPROP,PRPAO,IANTSYM)
C
C Sep 2004 HJAaJ, based on RDPROP in sirius/sirgp.F
C
C Read AO property integrals with label WORD from LUPROP
C Output:
C  PRPAO   : AO property integrals
C  IANTSYM : -1 if integral matrix is antisymmetric
C            +1 if matrix is symmetric
C             0 if matrix is square
C
#include "implicit.h"
#include "priunit.h"
      CHARACTER*8 WORD
      DIMENSION PRPAO(*)
C
C -- common blocks:
C  INFORB : NNBASX
C
#include "inforb.h"
C
      CHARACTER*8 RTNLBL(2)
      LOGICAL FNDLB2
C
      IF (LUPROP .LE. 0) STOP 'WESPROP ERROR: LUPROP .le. 0'
C
      REWIND LUPROP
      IF ( FNDLB2(WORD,RTNLBL,LUPROP)) THEN
         IF (RTNLBL(2).EQ.'SYMMETRI') THEN
            IANTSYM = 1
            CALL READT(LUPROP,NNBASX,PRPAO)
         ELSE IF (RTNLBL(2).EQ.'ANTISYMM') THEN
            IANTSYM = -1
            CALL READT(LUPROP,NNBASX,PRPAO)
         ELSE IF (RTNLBL(2).EQ.'SQUARE  ') THEN
            IANTSYM = 0
            CALL READT(LUPROP,N2BASX,PRPAO)
         ELSE
            WRITE (LUPRI,'(//2A)')
     *         ' WESPROP error - unknown RTNLBL: ',RTNLBL(2)
           STOP 'WESPROP error - unknown RTNLBL'
         END IF
      ELSE
         WRITE (LUPRI,'(//3A,I5)') ' >>> WESPROP error: property "',
     *      WORD,'" not found on AOPROPER.',LUPROP
         STOP 'WESPROP error: property not found on AOPROPER.'
      END IF
      RETURN
      END
