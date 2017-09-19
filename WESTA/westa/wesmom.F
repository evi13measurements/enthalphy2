#ifdef TODO
!/* Comdeck todo-may95 */
!(2) remove WESPAK ??
#endif
#  /* Deck wesmom */
      SUBROUTINE WESMOM(WRK,LWRK)
C
C  Purpose:
C    To control the calculation of the transition moment
C    < WF(1) | T | WF(2) > where
C    for ICASE = 1, T =  overlap monopole operator
C    for ICASE = 2, T =  dipole length operator   (x,y or z component)
C    for ICASE = 3, T =  dipole velocity operator (x,y or z component)
C    for ICASE = 4, T =  magnetic dipole operator (x,y or z component)
C
C Changed for determinant codes , November 1988
C
C  Declarations
C
#include "implicit.h"
C
      DIMENSION WRK(2)
#include "iratdef.h"
      PARAMETER (D0  = 0.0D0, D1  =  1.0D0, D2 = 2.0D0,
     1           DP5 = 0.5D0, DM1 = -1.0D0)
      PARAMETER ( KERROR = - 1 000 000 000 )
      PARAMETER ( THRDC = 1.D-10 )
      PARAMETER ( FACEV = 27.2117D0 )
C
C  PRIUNIT : LUPRI
C  WESINF : ..., MOMTYP(*), ...
C  WESTAP : LUSYMV(*),LUIT1V(1,*),LUSPCW, LUTMPW
C  INFORB : MULD2H(8,8), NORBT,NBAST,...
C  INFVAR : NCONF...
C  INFTAP : LUSYMB
C  INFDIM : SETCI parameters
C  INFINP : NACTEL.
C
#include "priunit.h"
#include "maxorb.h"
#include "wesinf.h"
#include "westap.h"
#include "inforb.h"
#include "infvar.h"
#include "inftap.h"
#include "infdim.h"
#include "infinp.h"
C
      LOGICAL EXPAND
C
      IPRECV = IPRTRM-10
C
      WRITE (LUPRI,'(A///)') '1'
      CALL HEADER('>>> OUTPUT FROM WESMOM <<<',-1)
      IF (NGEOM .NE. 1) THEN
         WRITE (LUPRI, '(///A/A//)')
     *      ' *** ERROR (WESMOM), transition moments elements',
     *              ' require one geometry only. Nothing done.'
         RETURN
      END IF
C
      NCONFX = NCONF
      LSYMP = MULD2H(LSYM1, LSYM2)
      WRITE (LUPRI,'(/2(/A,I3,F5.1),/)')
     1   ' Symmetry and spin multiplicity of state 1 :',LSYM1,SPIN1,
     2   ' Symmetry and spin multiplicity of state 2 :',LSYM2,SPIN2
C
C
C****************************************************************
C    Do we need to expand CI vectors to get a common CAS guga ?
C****************************************************************
C
         EXPAND = .TRUE.
C
C     -- We always need ao overlap matrix
C        (for calculating biorthogonal orbitals)
      MOMTYP(1) = .TRUE.
C
C****************************************************************
C Core allocation
C****************************************************************
C
C...OVLMOM...READMW...BIORAS
      KSAOAA = 1
      KDVEKT = KSAOAA + N2BASX
      IF (MOMTYP(1)) THEN
         KSAO = KDVEKT
         KDLX = KSAO + NNBASX
         KE1  = KDLX
      ELSE
         KDLX = KDVEKT
         KSAO = KERROR
      END IF
      IF (MOMTYP(2)) THEN
         KDLY = KDLX + NNBASX
         KDLZ = KDLY + NNBASX
         KDVX = KDLZ + NNBASX
         KE1  = KDVX
      ELSE
         KDVX = KDLX
         KDLX = KERROR
         KDLY = KERROR
         KDLZ = KERROR
      END IF
      IF (MOMTYP(3)) THEN
         KDVY = KDVX + NNBASX
         KDVZ = KDVY + NNBASX
         KDMX = KDVZ + NNBASX
         KE1  = KDMX
      ELSE
         KDMX = KDVX
         KDVX = KERROR
         KDVY = KERROR
         KDVZ = KERROR
      END IF
      IF (MOMTYP(4)) THEN
         KDMY = KDMX + NNBASX
         KDMZ = KDMY + NNBASX
         KE1  = KDMZ + NNBASX
      ELSE
         KE1  = KDMX
         KDMX = KERROR
         KDMY = KERROR
         KDMZ = KERROR
      END IF
C
C...Rawmom,Finmom,Energies
C
      KEMCCH = KE1
      KRWMOM = KEMCCH + NSTAT2
      KFIMOM = KRWMOM + 9*NSTAT2
      KE11   = KFIMOM + 9*NSTAT2
C
C...WRSPEC
C
      KCICH  = KE11
      KCMOCH = KCICH  + NCONF2
      KE12   = KCMOCH + NCMOT
C
C...READMW...BIORAS
C
      KCMOA1 = KE11
      KCMOA2 = KCMOA1 + NCMOT
      KVCA1  = KCMOA2 + NCMOT
      KVCA2  = KVCA1  + NORBT*NBAST
      KTRA1A = KVCA2  + NORBT*NBAST
      KTRA2A = KTRA1A + NASHT * NASHT
      KTIN1  = KTRA2A + NASHT * NASHT
      KTIN2  = KTIN1  + NISHT
      KW1    = KTIN2  + NISHT
      LW1    = LWRK   - KW1
C
C...CI vectors
C
      KCI1EX = KW1
      KCI2EX = KCI1EX + NCONFX
      KW2    = KCI2EX + NCONFX
      LW2    = LWRK   - KW2
C
C...SCALCI...WESTRACI/WESDM
C
      KCINDX = KW2
      KW30   = KCINDX + LCINDX
C
      KDTV   = KW30
      KW3    = KDTV   + N2ASHX
      LW3    = LWRK   - KW3
      IF (KW3.GT.LWRK) CALL ERRWRK('WESMOM.WESDM',-KW3,LWRK)
C
      KDAO   = KW3
      KDSPAO = KDAO   + N2BASX
      KDAPAO = KDSPAO + NNBASX
      KW4    = KDAPAO + NNBASX
      LW4    = LWRK   - KW4
      IF (LW4 .LT. N2BASX) CALL ERRWRK('WESMOM.DAOtrans',N2BASX,LW4)
C
C ... WESECV
C
      KW5    = KW3
      LW5    = LWRK - KW5
      IF (KW5.GT.LWRK) CALL ERRWRK('WESMOM.WESECV',KW5,LWRK)
C
C...WESSWF (get symmetric orthonormalization of wave functions)
C
      KSWF   = KE11
      KTSWF  = KSWF   + NSTAT2*(NSTAT2+1)/2
      KW44   = KTSWF  + NSTAT2*NSTAT2
      LW44   = LWRK   - KW44
C     ==========================================
C
C****************************************************************
C *** CALCULATE REQUESTED GEOMETRY NON-DISPLACED MOMENT INTEGRALS ***
C****************************************************************
C
      TIMINT = SECOND()
      CALL OVLMOM(MOMTYP,WRK(KDVEKT),1)
      TIMINT = SECOND() - TIMINT
      WRITE (LUPRI,'(///A/)') ' The following properties are evaluated:'
      IF(MOMTYP(1)) WRITE(LUPRI,'(A)')' Overlap (electric monopole)'
      IF(MOMTYP(2)) WRITE(LUPRI,'(A)')' Electric dipole (length form)'
      IF(MOMTYP(3)) WRITE(LUPRI,'(A)')' Electric dipole (velocity form)'
      IF(MOMTYP(4)) WRITE(LUPRI,'(A)')' Magnetic dipole'
      IF(IPRTRM .GT. 0) WRITE (LUPRI,'(//A,F10.2,A)')
     1   ' >>> TIME USED FOR AO MOMENT INTEGRALS:',TIMINT,' SECONDS.'
C
C Sort out A-A square block SAOAA out of packed SAO
C
      CALL DSPTSI(NBAST,WRK(KSAO),WRK(KSAOAA))
C
      IF (IPRTRM .GE. 8) THEN
         WRITE (LUPRI,'(//A)') ' AO overlap matrix'
         CALL OUTPAK(WRK(KSAO),NBAST,1,LUPRI)
         IF (IPRTRM .GE. 10) THEN
            WRITE (LUPRI,'(//A)') ' AO unpacked overlap matrix'
            CALL OUTPUT(WRK(KSAOAA),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (MOMTYP(2)) THEN
            WRITE(LUPRI,'(//A)')' AO electric x dipole (length) matrix'
            CALL OUTPAK(WRK(KDLX),NBAST,1,LUPRI)
            WRITE(LUPRI,'(//A)')' AO electric y dipole (length) matrix'
            CALL OUTPAK(WRK(KDLY),NBAST,1,LUPRI)
            WRITE(LUPRI,'(//A)')' AO electric z dipole (length) matrix'
            CALL OUTPAK(WRK(KDLZ),NBAST,1,LUPRI)
         END IF
         IF (MOMTYP(3)) THEN
            WRITE (LUPRI,'(//A)')
     1         ' AO electric dipole (velocity) matrix - x'
            CALL OUTPAK(WRK(KDVX),NBAST,1,LUPRI)
            WRITE (LUPRI,'(//A)')
     1         ' AO electric dipole (velocity) matrix - y'
            CALL OUTPAK(WRK(KDVY),NBAST,1,LUPRI)
            WRITE (LUPRI,'(//A)')
     1         ' AO electric dipole (velocity) matrix - z'
            CALL OUTPAK(WRK(KDVZ),NBAST,1,LUPRI)
         END IF
         IF (MOMTYP(4)) THEN
            WRITE (LUPRI,'(//A)') ' AO magnetic dipole matrix - x'
            CALL OUTPAK(WRK(KDMX),NBAST,1,LUPRI)
            WRITE (LUPRI,'(//A)') ' AO magnetic dipole matrix - y'
            CALL OUTPAK(WRK(KDMY),NBAST,1,LUPRI)
            WRITE (LUPRI,'(//A)') ' AO magnetic dipole matrix - z'
            CALL OUTPAK(WRK(KDMZ),NBAST,1,LUPRI)
         END IF
      END IF
C
C****************************************************************
C *** CONTROL OF TRANSITION MOMENT CALCULATION ***
C****************************************************************
C
C
C****************************************************************
C *** PREPARATION ***
C****************************************************************
C
C *** TEST INPUT:
C     TEST IF CMOA1 CONSIST OF ORTHONORMAL MO´S
C
      CALL WESTST(WRK(KSAO),DUMMY,WRK(KCMOA1),WRK(KW1),LW1)
C
C
C  Move final state ci vectors and orbitals to single file
C  and get final state energies.
C
      CALL WRSPEC(WRK(KCICH),WRK(KCMOCH),WRK(KEMCCH))
C
C  Get initial state ci vector, expand if required and then write
C  it on unit LUTMPW
C
      IF (LSYM1.NE.LSYM2) THEN
         LSYM1X = 1
         LSYM2X = 1
         NOSYM  = 1
      ELSE
         LSYM1X = LSYM1
         LSYM2X = LSYM2
         NOSYM  = 0
      END IF
      IF (NASHT .GT. 0) THEN
C Obtain string information for expanded graph
         CALL GETCIX(WRK(KCINDX),LSYM,LSYM,WRK(KW5),LW5,NOSYM)
         IF (EXPAND) THEN
            CALL WESGCI(1,ISTATV(1,1),NCONF1,WRK(KCI2EX),0,0,0,DUMMY)
C           (KCI2EX used as scratch)
            CALL EXPCID(WRK(KCI2EX),WRK(KCI1EX),LSYM1,LSYM1X,LUSYMV(1),
     &                  WRK(KCINDX),WRK(KW5),LW5,IPRECV)
C
            IF (IPRTRM .GE. 35) THEN
               WRITE (LUPRI,'(//A,I8)')
     *         ' CI vector no. 1 before WESECV, length =',NCONF1
               WRITE( LUPRI,'(2(1X,5F15.12/))')
     *         (WRK(KCI2EX-1+I),I=1,NCONF1)
               WRITE (LUPRI,'(//A,I8)')
     *         ' CI vector no. 1 after WESECV, length =',NCONFX
               WRITE( LUPRI,'(2(1X,5F15.12/))')
     *         (WRK(KCI1EX-1+I),I=1,NCONFX)
            END IF
         ELSE
            CALL WESGCI(1,ISTATV(1,1),NCONF1,WRK(KCI1EX),0,0,0,DUMMY)
         END IF
         EMCINI = EMCWES(1,1)
C
C        open LUTMPW
C
#if defined (SYS_CRAY)
         OPEN(LUTMPW,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *        FILE='fort.52')
#endif
#if defined (SYS_VAX)
         L0 = NCONFX/128 + 1
         OPEN(UNIT=LUTMPW,STATUS='NEW',FORM='UNFORMATTED',
     *        BLOCKSIZE=L0,INITIALSIZE=L0)
#endif
C
         REWIND LUTMPW
         CALL WRITT(LUTMPW,NCONFX,WRK(KCI1EX))
      END IF
C
C
C****************************************************************
C  *** START LOOP OVER SEVERAL FINAL STATES (OF SYMMETRY LSYM2)
C****************************************************************
C
      REWIND LUSPCW
      DO 300 ISTAT2 = 1,NSTAT2
C
C
C
C        Get ci vectors and orbitals for final state ISTAT2
C
         IF (NASHT .EQ. 0) THEN
            WRK(KCI2EX) = D1
         ELSE IF (EXPAND) THEN
            CALL READT(LUSPCW,NCONF2,WRK(KCI1EX))
            CALL EXPCID(WRK(KCI1EX),WRK(KCI2EX),LSYM2,LSYM2X,LUSYMV(2),
     &                  WRK(KCINDX),WRK(KW5),LW5,IPRECV)
            IF (IPRTRM .GE. 35) THEN
               WRITE (LUPRI,'(//A,I8)')
     *         ' CI vector no. 2 before WESECV, length =',NCONF2
               WRITE( LUPRI,'(2(1X,5F15.12/))')
     *         (WRK(KCI1EX-1+I),I=1,NCONF2)
               WRITE (LUPRI,'(//A,I8)')
     *         ' CI vector no. 2 after WESECV, length =',NCONFX
               WRITE (LUPRI,'(2(1X,5F15.12/))')
     *         (WRK(KCI2EX-1+I),I=1,NCONFX)
            END IF
         ELSE
            CALL READT(LUSPCW,NCONF2,WRK(KCI2EX))
         END IF
C
         IF (EXPAND) THEN
            NCONF  = NCONFX
         ELSE
            NCONF  = NCONF1
         ENDIF
C
         CALL READT(LUSPCW,NCMOT,WRK(KCMOA2))
C
C        Get ci vector and orbitals for initial state
C
         IF (NASHT .EQ. 0) THEN
            WRK(KCI1EX) = D1
         ELSE
            REWIND LUTMPW
            CALL READT(LUTMPW,NCONF,WRK(KCI1EX))
         END IF
         CALL READMW(IRDMOW,LUIT1V(1,1),WRK(KCMOA1))
C
         IF (IPRTRM .GE. 30) THEN
            WRITE (LUPRI,'(//A,I8)')
     *      ' CI vector no. 1 from SIRIUS, length =',NCONF
            WRITE (LUPRI,'(2(1X,5F15.12/))') (WRK(KCI1EX-1+I),I=1,NCONF)
            WRITE (LUPRI,'(//A,I8)')
     *      ' CI vector no. 2 from SIRIUS, length =',NCONF
            WRITE (LUPRI,'(2(1X,5F15.12/))') (WRK(KCI2EX-1+I),I=1,NCONF)
            WRITE (LUPRI,'(//A,I8)')
     *      ' MO:S for state no. 2 from SIRIUS, length =',NCMOT
            CALL PRORB(WRK(KCMOA2),.FALSE.,LUPRI)
         END IF

C        1.First: obtain biorthogonal orbitals
C                 obtain OVLA12 element needed for "SYMORT" correction.
C
         CALL HEADER('Construction of biorthogonal orbitals',1)
         WRITE (LUPRI,'(/A,I5//)') ' EXCITED STATE NO. ',ISTAT2
C
C        CALL BIORAS(CMO1,CMO2,VC1,VC2,SAO12,TIN1,TIN2,TRA1A,
C    *               TRA2A,WRK,LWRK)
C
C  We take only active part of TRA-matrices out of BIORAS +
C  the inactive diagonal parts for scaling.
C
         CALL BIORAS(WRK(KCMOA1),WRK(KCMOA2),WRK(KVCA1),WRK(KVCA2),
     *               WRK(KSAOAA),WRK(KTIN1),WRK(KTIN2),WRK(KTRA1A),
     *               WRK(KTRA2A),WRK(KW5),LW5)
C
C        CALL WESTRACI(NCVEC,CVEC,LCVEC,UMO,XNDXCI,WRK,LFREE)
C
         IF (NASHT .GT. 0) THEN
            CALL WESTRACI(1,WRK(KCI1EX),NCONF,WRK(KTRA1A),WRK(KCINDX),
     *                 WRK(KW5),LW5)
            CALL WESTRACI(1,WRK(KCI2EX),NCONF,WRK(KTRA2A),WRK(KCINDX),
     *                 WRK(KW5),LW5)
         END IF
C
C Scale CI vectors by (Tin1)**2 resp. (Tin2)**2
C
      FAC1 = D1
      FAC2 = D1
      DO 1110 I =1, NISHT
        FAC1 = FAC1 * WRK(KTIN1+I-1) * WRK(KTIN1+I-1)
        FAC2 = FAC2 * WRK(KTIN2+I-1) * WRK(KTIN2+I-1)
 1110 CONTINUE
C
         CALL DSCAL(NCONF,FAC1,WRK(KCI1EX),1)
         CALL DSCAL(NCONF,FAC2,WRK(KCI2EX),1)
C
C
         IF (LSYMP .EQ. 1 .AND. SPIN1 .EQ. SPIN2) THEN
            OVLA12 = DDOT(NCONF,WRK(KCI1EX),1,WRK(KCI2EX),1)
         ELSE
            OVLA12 = D0
         ENDIF
C
         IF (IPRTRM .GE. 30) THEN
            WRITE (LUPRI,'(//A,I8)')
     *      ' CI vector no. 1 after WESTRACI/SCALCI, length =',NCONFX
            WRITE (LUPRI,'(2(1X,5F15.12/))')(WRK(KCI1EX-1+I),I=1,NCONFX)
            WRITE (LUPRI,'(//A,I8)')
     *      ' CI vector no. 2 after WESTRACI, length =',NCONFX
            WRITE (LUPRI,'(2(1X,5F15.12/))')(WRK(KCI2EX-1+I),I=1,NCONFX)
         END IF
C
C
C        2.Second:
C           obtain transition density matrix. If the symmetries of
C           the two states are equal call WESDM directly, otherwise
C           expand both initial and final state CI vectors in nosymmetry
C           format (2 calls to WESECV) before calling WESDM
C
         CALL HEADER('Transition density matrix',1)
C
         IF (SPIN1 .NE. SPIN2) THEN
            WRITE (LUPRI,'(//A,2F5.1/A)')
     *      ' *** SPIN OF STATES 1 AND 2 :',SPIN1,SPIN2,
     *      ' *** NO SPIN MOMENTS IMPLEMENTED IN WESMOM, NOTHING DONE.'
           RETURN
         END IF
C
         IF (NASHT .GT. 0) THEN
C
C
C           CALL WESDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,
C    &                 ISPIN1,TDM,XINDX,WORK,KFREE,LFREE)
            KFREE = 1
            LFREE = LW3
            ISPIN1 = 0
            CALL WESDM(LSYM1X,LSYM2X,NCONFX,NCONFX,
     &                 WRK(KCI1EX),WRK(KCI2EX),WRK(KDTV),
     &                 ISPIN1,.TRUE.,WRK(KCINDX),
     &                 WRK(KW3),KFREE,LFREE)
            WRITE (LUPRI,'(//A)')
     1         ' TRANSITION DENSITY MATRIX OVER ACTIVE MO:S'
            CALL OUTPUT(WRK(KDTV),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
C
C        3.Third:
C           Backtransform transition d.m. to ao basis (contravariant)
C
C
C        ACTIVE PART OF CONTRAVARIANT AO TRANSITION DENSITY MATRIX
C
         IF (NASHT .GT. 0) THEN
            JVCA1 = KVCA1 + NISHT*NBAST
            JVCA2 = KVCA2 + NISHT*NBAST
            CALL DGEMM('N','N',NBAST,NASHT,NASHT,1.D0,
     &                 WRK(JVCA1),NBAST,
     &                 WRK(KDTV), NASHT,0.D0,
     &                 WRK(KW4),  NBAST)
            CALL DGEMM('N','T',NBAST,NBAST,NASHT,1.D0,
     &                 WRK(KW4),  NBAST,
     &                 WRK(JVCA2),NBAST,0.D0,
     &                 WRK(KDAO), NBAST)
            IF (IPRTRM .GE. 6) THEN
               WRITE (LUPRI,'(//2A)')
     1         ' TRANSITION DENSITY MATRIX OVER ACTIVE MO:S',
     2         ' IN CONTRAVARIANT AO BASIS'
              CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            END IF
         ELSE
            CALL DZERO(WRK(KDAO), N2BASX)
         END IF
C
C        Add INACTIVE PART OF CONTRAVARIANT AO TRANSITION DENSITY MATRIX
C
         IF ( NISHT .GT. 0 .AND. ABS(OVLA12) .GT. THRDC ) THEN
            KT1 = KDAO + N2BASX
            KT2 = KT1  + N2BASX
            KTE = KT2  + N2BASX
            IF (KTE .GT. LWRK) CALL ERRWRK('WESMOM.DCAO',KTE,LWRK)
            CALL DUNIT(WRK(KT1),NISHT)
            CALL DSCAL(NISHT,(D2*OVLA12),WRK(KT1),(NISHT+1))
            CALL DGEMM('N','N',NBAST,NISHT,NISHT,1.D0,
     &                 WRK(KVCA1),NBAST,
     &                 WRK(KT1),  NISHT,0.D0,
     &                 WRK(KT2),  NBAST)
            CALL DGEMM('N','T',NBAST,NBAST,NISHT,1.D0,
     &                 WRK(KT2),  NBAST,
     &                 WRK(KVCA2),NBAST,1.D0,
     &                 WRK(KDAO), NBAST)
C
         END IF
         IF (IPRTRM .GE. 6) THEN
            WRITE (LUPRI,'(//2A)') ' FULL TRANSITION DENSITY MATRIX',
     1         ' IN CONTRAVARIANT AO BASIS'
            IF (NISHT .GT. 0) THEN
              CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            ELSE
               WRITE (LUPRI,'(/A/A/)')
     1         ' IS EQUAL TO ABOVE MATRIX FROM ACTIVE MO:S',
     2         ' (NO INACTIVE ORBITALS IN THIS CALCULATION).'
            END IF
         END IF
C
C    4.Fourth: Construct symmetrized and anti-symmetrized
C              DAO matrices in packed form
C
         CALL WESPAK(WRK(KDAO),WRK(KDSPAO),WRK(KDAPAO),NBAST)
         IF (IPRTRM .GE. 7) THEN
            WRITE (LUPRI,'(//A/A)') ' FULL TRANSITION DENSITY MATRIX',
     1         ' IN CONTRAVARIANT AO BASIS, SYMMETRIC PART'
            CALL OUTPAK(WRK(KDSPAO),NBAST,1,LUPRI)
            WRITE (LUPRI,'(//A/A)') ' FULL TRANSITION DENSITY MATRIX',
     1         ' IN CONTRAVARIANT AO BASIS, ANTISYMMETRIC PART'
            CALL OUTPAK(WRK(KDAPAO),NBAST,1,LUPRI)
         END IF
C
C    5.Fifth: take the trace of DTV and Transformed property (MO) matrix
C             to obtain final transition moment
C
         CALL HEADER('Transition moments',1)
C
         WRITE (LUPRI,'(//A,I5)') ' Final state No. ', ISTAT2
C
C
         WRITE(LUPRI,'(//A,F20.10)')' < WF(1) /        WF(2) > =',OVLA12
         IF (MOMTYP(2)) THEN
            DLXA12 = DDOT(NNBASX,WRK(KDSPAO),1,WRK(KDLX),1)
            DLYA12 = DDOT(NNBASX,WRK(KDSPAO),1,WRK(KDLY),1)
            DLZA12 = DDOT(NNBASX,WRK(KDSPAO),1,WRK(KDLZ),1)
            WRK(KRWMOM + (ISTAT2-1))            = DLXA12
            WRK(KRWMOM + (NSTAT2 + ISTAT2-1))   = DLYA12
            WRK(KRWMOM + (2*NSTAT2 + ISTAT2-1)) = DLZA12
            WRITE (LUPRI,'(//A//,(A,F20.10))')
     9         ' Electric dipole (length) :',
     1         ' < WF(1) /   X  / WF(2) > =',DLXA12,
     2         ' < WF(1) /   Y  / WF(2) > =',DLYA12,
     3         ' < WF(1) /   Z  / WF(2) > =',DLZA12
         END IF
         IF (MOMTYP(3)) THEN
            DVXA12 = DDOT(NNBASX,WRK(KDAPAO),1,WRK(KDVX),1)
            DVYA12 = DDOT(NNBASX,WRK(KDAPAO),1,WRK(KDVY),1)
            DVZA12 = DDOT(NNBASX,WRK(KDAPAO),1,WRK(KDVZ),1)
            WRK(KRWMOM + (3*NSTAT2 + ISTAT2-1)) = DVXA12
            WRK(KRWMOM + (4*NSTAT2 + ISTAT2-1)) = DVYA12
            WRK(KRWMOM + (5*NSTAT2 + ISTAT2-1)) = DVZA12
            WRITE (LUPRI,'(//A//,(A,F20.10))')
     9         ' Electric dipole (velocity) :',
     1         ' < WF(1) / d/dX / WF(2) > =',DVXA12,
     2         ' < WF(1) / d/dY / WF(2) > =',DVYA12,
     3         ' < WF(1) / d/dZ / WF(2) > =',DVZA12
            OMEA12 = WRK(KEMCCH-1+ISTAT2) - EMCINI
            WRITE (LUPRI,'(/A,F20.10)')
     1         ' Excitation energy  OMEGA =',OMEA12
            IF (OMEA12 .NE. D0) THEN
               TVXA12 = DVXA12 / OMEA12
               TVYA12 = DVYA12 / OMEA12
               TVZA12 = DVZA12 / OMEA12
               WRITE (LUPRI,'(/A//,(A,F20.10))')
     9         ' Electric dipole (velocity) divided by exc. energy :',
     1         ' < 1 / d/dX / 2 > / OMEGA =',TVXA12,
     2         ' < 1 / d/dY / 2 > / OMEGA =',TVYA12,
     1         ' < 1 / d/dZ / 2 > / OMEGA =',TVZA12
            END IF
         END IF
         IF (MOMTYP(4)) THEN
            DMXA12 = DDOT(NNBASX,WRK(KDAPAO),1,WRK(KDMX),1)
            DMYA12 = DDOT(NNBASX,WRK(KDAPAO),1,WRK(KDMY),1)
            DMZA12 = DDOT(NNBASX,WRK(KDAPAO),1,WRK(KDMZ),1)
            WRK(KRWMOM + (6*NSTAT2 + ISTAT2-1)) = DMXA12
            WRK(KRWMOM + (7*NSTAT2 + ISTAT2-1)) = DMYA12
            WRK(KRWMOM + (8*NSTAT2 + ISTAT2-1)) = DMZA12
            WRITE (LUPRI,'(//A//,(A,F20.10))')
     9         ' Magnetic dipole / angular non-adiabatic coupling :',
     1         ' < WF(1) / L(X) / WF(2) > =',DMXA12,
     2         ' < WF(1) / L(Y) / WF(2) > =',DMYA12,
     3         ' < WF(1) / L(Z) / WF(2) > =',DMZA12
         END IF
C
C End of ISTAT2 loop
C
  300 CONTINUE
C
C     Close and delete temporary file LUTMPW
C
      IF (NASHT .GT. 0) CLOSE(UNIT=LUTMPW,STATUS='DELETE')
C
C************************************************
C   Symmetric orthonormalization of final states
C************************************************
C
      IF (SYMORT) THEN
C
        LUSYMB = LUSYMV(2)
        NCONF  = NCONF2
C
C
        CALL WESSWF(WRK(KSAOAA),WRK(KSWF),WRK(KTSWF),WRK(KW44),LW44)
C       CALL WESSWF(SAOAA,SWF,TSWF,WRK,LWRK)
        LUSYMB = -1
C
C       Close and delete temporary file LUSPCW
C
        CLOSE(UNIT=LUSPCW,STATUS='DELETE')
C
        CALL HEADER('<<<<< FINAL RESULTS FROM WESTA MOMENTS >>>>>',-1)
        CALL TIMDAT
C************************************************
C       Printout of results
C************************************************
C
C       Printout of transition energies
C
        WRITE (LUPRI,'(/A,F12.6)')
     *     ' GROUND STATE ENERGY (AU)                :',EMCINI
        IF (WRK(KEMCCH) .NE. D0)
     *     WRITE (LUPRI,'(/A,F12.6)')
     *     ' TRANSITION ENERGY FOR LOWEST STATE (AU) :',
     *     (WRK(KEMCCH)-EMCINI)
        WRITE (LUPRI,'(/A,F12.6)')
     *     ' TRANSITION ENERGY FOR LOWEST STATE (EV) :',
     *     FACEV * (WRK(KEMCCH)-EMCINI)
C
C       Printout of final state overlap matrix.
C
        WRITE (LUPRI,'(//A)') ' Spectral overlap matrix:'
        CALL OUTPAK (WRK(KSWF),NSTAT2,1,LUPRI)
        WRITE (LUPRI,'(//A)') ' S-1/2 tranformation:'
        CALL OUTPUT (WRK(KTSWF),1,NSTAT2,1,NSTAT2,NSTAT2,NSTAT2,1,LUPRI)
C
C       Transform moment to symmetrically orthonomalized
C       wave function basis and write out final results.
C
        IF (MOMTYP(2)) THEN
C
           CALL DGEMM('N','N',NSTAT2,1,NSTAT2,1.D0,
     &                 WRK(KTSWF), NSTAT2,
     &                 WRK(KRWMOM),NSTAT2,0.D0,
     &                 WRK(KFIMOM),NSTAT2)
           CALL DGEMM('N','N',NSTAT2,1,NSTAT2,1.D0,
     &                 WRK(KTSWF), NSTAT2,
     &                 WRK(KRWMOM + NSTAT2),NSTAT2,0.D0,
     &                 WRK(KFIMOM + NSTAT2),NSTAT2)
           CALL DGEMM('N','N',NSTAT2,1,NSTAT2,1.D0,
     &                 WRK(KTSWF), NSTAT2,
     &                 WRK(KRWMOM + 2*NSTAT2),NSTAT2,0.D0,
     &                 WRK(KFIMOM + 2*NSTAT2),NSTAT2)
C
         CALL HEADER('<<<<< RESULTS: DIPOLE LENGTH MOMENTS >>>>>',-1)
C
         WRITE (LUPRI,'(//A,T10,A/T10,A/)') ' STATE',
     *     'X-TRANS ENERGY    X-TRANS INTENSITY    X-TRANS INTENSITY',
     *     '    ( eV )         ( non-ort )           ( sym-ort )'
          DO 1001 ISTAT2 =1,NSTAT2
             XINT1 = WRK(KRWMOM-1+ISTAT2)
             XINT2 = WRK(KFIMOM-1+ISTAT2)
             YINT1 = WRK(KRWMOM-1+ISTAT2 + NSTAT2)
             YINT2 = WRK(KFIMOM-1+ISTAT2 + NSTAT2)
             ZINT1 = WRK(KRWMOM-1+ISTAT2 + 2*NSTAT2)
             ZINT2 = WRK(KFIMOM-1+ISTAT2 + 2*NSTAT2)
             XINT1 = XINT1 * XINT1
             XINT2 = XINT2 * XINT2
             YINT1 = YINT1 * YINT1
             YINT2 = YINT2 * YINT2
             ZINT1 = ZINT1 * ZINT1
             ZINT2 = ZINT2 * ZINT2
             OMEGA = FACEV*(WRK(KEMCCH+ISTAT2-1) - EMCINI)
             WRITE (LUPRI,1010) ISTAT2,'-X',OMEGA,XINT1,XINT2
             WRITE (LUPRI,1010) ISTAT2,'-Y',OMEGA,YINT1,YINT2
             WRITE (LUPRI,1010) ISTAT2,'-Z',OMEGA,ZINT1,ZINT2
 1001     CONTINUE
C
 1010     FORMAT(I5,A,F14.6,F19.6,F21.6)
C
        ENDIF
C
        IF (MOMTYP(3)) THEN
C
           CALL DGEMM('N','N',NSTAT2,1,NSTAT2,1.D0,
     &                 WRK(KTSWF), NSTAT2,
     &                 WRK(KRWMOM + 3*NSTAT2),NSTAT2,0.D0,
     &                 WRK(KFIMOM + 3*NSTAT2),NSTAT2)
           CALL DGEMM('N','N',NSTAT2,1,NSTAT2,1.D0,
     &                 WRK(KTSWF), NSTAT2,
     &                 WRK(KRWMOM + 4*NSTAT2),NSTAT2,0.D0,
     &                 WRK(KFIMOM + 4*NSTAT2),NSTAT2)
           CALL DGEMM('N','N',NSTAT2,1,NSTAT2,1.D0,
     &                 WRK(KTSWF), NSTAT2,
     &                 WRK(KRWMOM + 5*NSTAT2),NSTAT2,0.D0,
     &                 WRK(KFIMOM + 5*NSTAT2),NSTAT2)
C
         CALL HEADER('<<<<< RESULTS: DIPOLE VELOCITY MOMENTS >>>>>',-1)
C
         WRITE (LUPRI,'(//A,T10,A/T10,A/)') ' STATE',
     *     'X-TRANS ENERGY    X-TRANS INTENSITY    X-TRANS INTENSITY',
     *     '    ( eV )         ( non-ort )           ( sym-ort )'
          DO 1002 ISTAT2 = 1,NSTAT2
             XINT1 = WRK(KRWMOM-1+ISTAT2 + 3*NSTAT2)
             XINT2 = WRK(KFIMOM-1+ISTAT2 + 3*NSTAT2)
             YINT1 = WRK(KRWMOM-1+ISTAT2 + 4*NSTAT2)
             YINT2 = WRK(KFIMOM-1+ISTAT2 + 4*NSTAT2)
             ZINT1 = WRK(KRWMOM-1+ISTAT2 + 5*NSTAT2)
             ZINT2 = WRK(KFIMOM-1+ISTAT2 + 5*NSTAT2)
             OMEGA = WRK(KEMCCH+ISTAT2-1) - EMCINI
          IF (OMEGA .EQ. D0) GO TO 1002
             FAC   = D1 / (OMEGA * OMEGA)
             XINT1 = FAC * XINT1 * XINT1
             XINT2 = FAC * XINT2 * XINT2
             YINT1 = FAC * YINT1 * YINT1
             YINT2 = FAC * YINT2 * YINT2
             ZINT1 = FAC * ZINT1 * ZINT1
             ZINT2 = FAC * ZINT2 * ZINT2
             OMEGA = FACEV * OMEGA
             WRITE (LUPRI,1010) ISTAT2,'-X',OMEGA,XINT1,XINT2
             WRITE (LUPRI,1010) ISTAT2,'-Y',OMEGA,YINT1,YINT2
             WRITE (LUPRI,1010) ISTAT2,'-Z',OMEGA,ZINT1,ZINT2
 1002     CONTINUE
C
        ENDIF
C
        IF (MOMTYP(4)) THEN
C
           CALL DGEMM('N','N',NSTAT2,1,NSTAT2,1.D0,
     &                 WRK(KTSWF), NSTAT2,
     &                 WRK(KRWMOM + 6*NSTAT2),NSTAT2,0.D0,
     &                 WRK(KFIMOM + 6*NSTAT2),NSTAT2)
           CALL DGEMM('N','N',NSTAT2,1,NSTAT2,1.D0,
     &                 WRK(KTSWF), NSTAT2,
     &                 WRK(KRWMOM + 7*NSTAT2),NSTAT2,0.D0,
     &                 WRK(KFIMOM + 7*NSTAT2),NSTAT2)
           CALL DGEMM('N','N',NSTAT2,1,NSTAT2,1.D0,
     &                 WRK(KTSWF), NSTAT2,
     &                 WRK(KRWMOM + 8*NSTAT2),NSTAT2,0.D0,
     &                 WRK(KFIMOM + 8*NSTAT2),NSTAT2)
C
         CALL HEADER('<<<<< RESULTS: MAGNETIC DIPOLE MOMENTS >>>>>',-1)
C
         WRITE (LUPRI,'(//A,T10,A/T10,A/)') ' STATE',
     *     'X-TRANS ENERGY    X-TRANS INTENSITY    X-TRANS INTENSITY',
     *     '    ( eV )         ( non-ort )           ( sym-ort )'
          DO 1003 ISTAT2 =1,NSTAT2
             XINT1 = WRK(KRWMOM-1+ISTAT2 + 6*NSTAT2)
             XINT2 = WRK(KFIMOM-1+ISTAT2 + 6*NSTAT2)
             YINT1 = WRK(KRWMOM-1+ISTAT2 + 7*NSTAT2)
             YINT2 = WRK(KFIMOM-1+ISTAT2 + 7*NSTAT2)
             ZINT1 = WRK(KRWMOM-1+ISTAT2 + 8*NSTAT2)
             ZINT2 = WRK(KFIMOM-1+ISTAT2 + 8*NSTAT2)
             XINT1 = XINT1 * XINT1
             XINT2 = XINT2 * XINT2
             YINT1 = YINT1 * YINT1
             YINT2 = YINT2 * YINT2
             ZINT1 = ZINT1 * ZINT1
             ZINT2 = ZINT2 * ZINT2
             OMEGA = FACEV*(WRK(KEMCCH+ISTAT2-1) - EMCINI)
             WRITE (LUPRI,1010) ISTAT2,'-X',OMEGA,XINT1,XINT2
             WRITE (LUPRI,1010) ISTAT2,'-Y',OMEGA,YINT1,YINT2
             WRITE (LUPRI,1010) ISTAT2,'-Z',OMEGA,ZINT1,ZINT2
 1003     CONTINUE
C
        END IF
C
      END IF
C
C End of WESMOM
C
      RETURN
      END
#  /* Deck wesdm */
      SUBROUTINE WESDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,
     &                 ISPIN1,TDM,XINDX,WORK,KFREE,LFREE)
C
C     CONSTRUCT ONE (RHO1) PARTICLE DENSITY MATRIX
C     MATRICES TO BE USED IN WESTA CALCULATION.
C
C     (based on RSPDM 950516)
C
C  OUTPUT:
C
C  RHO1(KL)   = (L/ E(KL) /R)
C
#include "implicit.h"
C
      DIMENSION CL(*), CR(*), RHO1(NASHT,*)
      DIMENSION XINDX(*), WORK(*)
C
      PARAMETER ( D0 = 0.0D0 )
C
C Used from common blocks:
C   PRIUNIT : LUPRI
C   WESINF : IPRDEF
C   INFINP : FLAG(*)
C   INFORB : NASHT,N2ASHX,NNASHX
C
#include "priunit.h"
#include "wesinf.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
C
C
      INTEGER RHOTYP
      LOGICAL USEPSM, CSFEXP, TDM, NORHO2
      PARAMETER (RHOTYP = 2, ISPIN2 = 0,
     &    USEPSM = .FALSE., NORHO2 = .TRUE.)
C     RHOTYP = 2 for non-symmetrized density matrix
C     ISPIN2 = 0 for singlet, spin symmetry of RHO1 is then
C        determined by ISPIN1
C     USEPSM is irrelevant for RHO1 (only used for RHO2)
C
C     Treat NASHT .eq. 1 as a special case
C
      IPRDNS = IPRDEF
C
      IF (NASHT .EQ. 0) GO TO 9999
      IF (NASHT .EQ. 1) THEN
         IF ( ILSYM.EQ.IRSYM ) THEN
            RHO1(1,1) = CL(1) * CR(1)
C           RHO1 = 1.0*CL(1)*CR(1) both for ISPIN1 = 0 and ISPIN1 = 1
         ELSE
            RHO1(1,1) = D0
         END IF
         GO TO 9999
      END IF
C
C     Set CSFEXP
C
      CSFEXP = .NOT.FLAG(27)
      IF (ISPIN1 .NE. 0) THEN
C        950516-hjaaj: detcin06.su checks for IDMTYP .eq. 1 or 2
C        in DENSID, spin density corresponds to IDMTYP = -1.
C        HJMAERKE: change DENSID to determine IDMTYP from
C        ISPIN1 and ISPIN2
C
         WRITE (LUPRI,'(//A/A)')
     &   'WESDM ERROR: detcin06.su version cannot calculate'//
     &   ' spin density', '16-May-1995 hjaaj.'
         STOP 'WESDM ERROR: detcin06.su cannot calc. spin dens.'
      END IF
C
C
      IF ( IPRDNS.GT.65 ) THEN
         WRITE(LUPRI,'(/A)')' ***WESDM BEFORE CALLING DENSID'
         WRITE(LUPRI,'(/A,/I6,4I7,I8)')
     *     ' ILSYM  IRSYM  NCLDIM  NCRDIM  KFREE   LFREE :',
     *       ILSYM, IRSYM, NCLDIM, NCRDIM, KFREE,  LFREE
      END IF
      IF ( IPRDNS.GT.120 ) THEN
         WRITE(LUPRI,'(/A)')' *WESDM* LEFT REFERENCE VECTOR'
         CALL OUTPUT(CL,1,NCLDIM,1,1,NCLDIM,1,1,LUPRI)
         WRITE(LUPRI,'(/A)')' *WESDM* RIGHT REFERENCE VECTOR'
         CALL OUTPUT(CR,1,NCRDIM,1,1,NCRDIM,1,1,LUPRI)
      END IF
      CALL SETVEC(RHO1,D0,N2ASHX)
      CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
     *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     *            XINDX,WORK,KFREE,LFREE)
C     CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
C    *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
C    *            XNDXCI,WORK,KFREE,LFREE)
C
C
      IF ( IPRDNS.GT.65 ) THEN
         WRITE(LUPRI,'(/A)')' *** LEAVING WESDM'
         WRITE(LUPRI,'(/A,/I6,4I7,I8)')
     *     ' ILSYM  IRSYM  NCLDIM  NCRDIM  KFREE   LFREE :',
     *       ILSYM, IRSYM, NCLDIM, NCRDIM, KFREE,  LFREE
      END IF
      IF ( IPRDNS .GE. 10 ) THEN
         WRITE (LUPRI,1100)
         CALL OUTPUT(RHO1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      ENDIF
 1100 FORMAT(/' DV = One-el. density matrix, active part, MO-basis')
C
 9999 RETURN
C     ... end of WESDM.
      END
#  /* Deck wespak */
      SUBROUTINE WESPAK(DMAT,DSP,DAP,NBAST)
C
C  13-Feb-1986 hjaaj --
C
C  Construct symmetrized and anti-symmetrized packed versions
C  of the contravariant transition density matrix.
C
#include "implicit.h"
      DIMENSION DMAT(NBAST,2), DSP(2), DAP(2)
      PARAMETER ( D0 = 0.0D0 )
C
      IJ = 0
      DO 200 J = 1,NBAST
         DO 100 I = 1,(J-1)
            IJ = IJ + 1
            DSP(IJ) = DMAT(I,J) + DMAT(J,I)
            DAP(IJ) = DMAT(I,J) - DMAT(J,I)
  100    CONTINUE
         IJ = IJ + 1
         DSP(IJ) = DMAT(I,I)
         DAP(IJ) = D0
  200 CONTINUE
C
      RETURN
      END
#  /* Deck ovlmom */
      SUBROUTINE OVLMOM(LCASE,DVEKT,NGEOM)
C
C OVLMOM from 30-Oct-1985 /HA,HJAAJ
C 17-May-1995: modified to call RDPROP instead of calculating integrals
C 22-Sep-2004: modified to call WESPROP based on RDPROP
C
C Purpose:
C     To compute overlap, electric and magnetic dipole monent integrals
C     over atomic basis functions.
C     electric dipole integrals can be calculated in length or
C     velocity form.
C     NGEOM = 1 : one geometry
C     NGEOM = 2 : two displaced geometries
C     MOMINT constructs the integral matrix with no assumption of
C     symmetry
C     LCASE(1): overlap integrals
C     LCASE(2): dipole length integrals   ( r - C)
C     LCASE(3): dipole velocity integrals ( P )
C     LCASE(4): magnetic dipole integrals ( L = (r - C) X P )
C
#include "implicit.h"
      LOGICAL LCASE(4)
      DIMENSION DVEKT(*)
C
#include "inforb.h"
#include "westap.h"
C
      LOGICAL OVERLA,DIPLEN,DIPVEL,DIPMAG
C
C
      IF (NGEOM .NE. 1) THEN
         STOP 'OVLMOM ERROR: NGEOM .ne. 1'
      END IF
      NBT2  = NGEOM * NBAST
      NNBT2 = NBT2*(NBT2+1)/2
C
      OVERLA = LCASE(1)
      DIPLEN = LCASE(2)
      DIPVEL = LCASE(3)
      DIPMAG = LCASE(4)
      NDV = 0
      IF (OVERLA) THEN
         NDV = NDV + 1
         JDLX = 1 + NNBT2
         CALL WESPROP('OVERLAP ',LUPROPV(1),DVEKT(1),IANTSYM)

      ELSE
         JDLX = 1
      ENDIF
      IF (DIPLEN) THEN
         JDLY = JDLX + NNBT2
         JDLZ = JDLY + NNBT2
         JDVX = JDLZ + NNBT2
         CALL WESPROP('XDIPLEN ',LUPROPV(1),DVEKT(JDLX),IANTSYM)
         CALL WESPROP('YDIPLEN ',LUPROPV(1),DVEKT(JDLY),IANTSYM)
         CALL WESPROP('ZDIPLEN ',LUPROPV(1),DVEKT(JDLZ),IANTSYM)
         NDV  = NDV  + 3
      ELSE
         JDVX = JDLX
      END IF
      IF (DIPVEL) THEN
         JDVY = JDVX + NNBT2
         JDVZ = JDVY + NNBT2
         JDMX = JDVZ + NNBT2
         CALL WESPROP('XDIPVEL ',LUPROPV(1),DVEKT(JDVX),IANTSYM)
         CALL WESPROP('YDIPVEL ',LUPROPV(1),DVEKT(JDVY),IANTSYM)
         CALL WESPROP('ZDIPVEL ',LUPROPV(1),DVEKT(JDVZ),IANTSYM)
         NDV  = NDV  + 3
      ELSE
         JDMX = JDVX
      END IF
      IF (DIPMAG) THEN
         JDMY = JDMX + NNBT2
         JDMZ = JDMY + NNBT2
         CALL WESPROP('XANGMOM ',LUPROPV(1),DVEKT(JDMX),IANTSYM)
         CALL WESPROP('YANGMOM ',LUPROPV(1),DVEKT(JDMY),IANTSYM)
         CALL WESPROP('ZANGMOM ',LUPROPV(1),DVEKT(JDMZ),IANTSYM)
         NDV  = NDV + 3
      ENDIF
C
      RETURN
      END
