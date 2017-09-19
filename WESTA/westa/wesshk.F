#  /* Deck wesshk */
      SUBROUTINE WESSHK(WRK,LWRK)
C
C 21-May-1986 Hans Agren and Hans Joergen Aa. Jensen
C  4-Dec-1986
C  2-May-1989 Changes for RAS. Ha, Jo.
C
C Purpose:
C   Calculate electron annihilation (shake-up) or creation (BIS)
C   overlap amplitudes and intensities
C
C   IACFLG = 0  :  No expansion
C          = 1  :  Annihilation. Expand N el. wfn to N-1 el. wfn
C          = 2  :  Creation.     Expand N el. wfn to N+1 el. wfn
C
C   IACSYM = symmetry for annihilation/creation
C   IACORB = orbital  for annihilation/creation
C
C  Outline:
C
C  Input for N-1 el. system. The annihilated electron is active also
C  in GS String information (LUSYMV(X),X=1,2) set up for N electrons
C  (GS) and N-1/N+1 electrons (Final states).
C
C  Question: how to deal with Reordering of inactive to active orbital
C            in the ground state (the one where el. is annih/created ?)
#include "implicit.h"
      DIMENSION WRK(1)
      DIMENSION SCR(600),ISCR(600)
C
#include "iratdef.h"
      PARAMETER ( D0 = 0.0D0, DUMMY = 1.0D20,  FACEV = 27.2117D0)
C
C
C Used from common blocks:
C   PRIUNIT : LUPRI
C   WESINF : NSTAT2, NCONF2, IRDMOW, EMCWES(1,1),
C            SYMORT, FROZEN, REDUCE
C   WESTAP : LUSPCW, LUIT1V(1,1), LUSYMV(*)
C   INFORB : NBAST,N2BASX, NCMOT
C   INFVAR : NCONF
C   INFDIM : NISHDI, NASHDI
C   INFTAP : LUSYMB
C   INFINP : NACTEL
C
#include "priunit.h"
#include "maxorb.h"
#include "wesinf.h"
#include "westap.h"
#include "inforb.h"
#include "infvar.h"
#include "infdim.h"
#include "inftap.h"
#include "infinp.h"
C
c
C
C     -- We always need ao overlap matrix
C        (for calculating biorthogonal orbitals)
      MOMTYP(1) = .TRUE.
c
      IPRECV = IPRSHK-10
      IF (.NOT.SHAKE) THEN
         WRITE(6,'(A)') ' WESHK called but not SHAKE '
         STOP ' WESSHK called but not SHAKE '
      END IF
      IF (IPRSHK .GT. 0) THEN
         WRITE (LUPRI,'(A///)') '1'
         CALL HEADER('<<<<< WESTA SHAKE CALCULATION >>>>>',-1)
      END IF
      IF (NSTAT2 .EQ. 0) THEN
         WRITE (LUPRI,'(//A/A/)')
     *  ' WESTA  SHAKE CALC. FINISHED, NO FINAL STATES !',
     *  ' (NSTAT2 = 0)'
         RETURN
      END IF
C
C first CALL SETCI to get NCONF2 and LCINDX
C LSYM1 = LSYM2  but we expand with NOSYM = 1 (?)
C
c     write(6,*) ' WESSHK ; NCONF1,NCONF2 = ', NCONF1,NCONF2
      NOSYM = 1
      CALL SETCI(NCONFX,NCDETS,LSYM1,WRK,LWRK,NOSYM)
c     WRITE(6,*) ' NCONFX = ', NCONFX
C     NCONF1 = NCONFX
C
C     Core allocation:
C
C  .. 1 .. general
C        EMCCH  : core hole energies
C        RAWMOM : raw moments
C        FINMOM : final moments (to symm. orthog. core hole states)
C        SAOU   : ao overlap matric stored in symmetry triangels (OCE)
C        SAO    : ao overlap matrix (full matrix)
C
C
      KEMCCH = 1
      KRWMOM = KEMCCH + NSTAT2
      KFIMOM = KRWMOM + NSTAT2
      KSAOU  = KFIMOM + NSTAT2
      KSAO   = KSAOU  + NNBASX
      KW1    = KSAO   + N2BASX
C
C  .. 2.1 .. WRSPEC (write final CI/MO:s unexpanded on LUSPCW)
C
      KCICHH = KW1
      KCMOCH = KCICHH + NCONF2
      KW2    = KCMOCH + NCMOT
      LW2    = LWRK - KW2

C
C...GS orbitals and expanded GS and FS CI vectors.
C
      KCMOGS = KW2
      KCMGS1 = KCMOGS + NCMOT
      KCINDX = KCMGS1 + NCMOT
      KCI1EX = KCINDX + LCINDX
      KCI2EX = KCI1EX + NCONFX
      KW3    = KCI2EX + NCONFX
      LW3    = LWRK   - KW3
C
      IF (KW3 .GT. LWRK) CALL ERRWRK('WESSHK.EXPCD2',KW3,LW3)
C
C  .. 4 .. WESSWF (get symmetric orthonormalization of wave functions)
C
      KSWF   = KW1
      KTSWF  = KSWF   + NSTAT2*(NSTAT2+1)/2
      KW4    = KTSWF  + NSTAT2*NSTAT2
      LW4    = LWRK   - KW4
C
C fix AO overlap matrix ==========================================
C
C  If OCEINT then open LUONEL file
C
      IF (OCEINT) THEN
C
#if defined (SYS_VAX)
      OPEN(UNIT=LUONEL,STATUS='OLD',FORM='UNFORMATTED',
     *        READONLY,SHARED,BLOCKSIZE=7680)
#else
      OPEN(LUONEL,STATUS='OLD',FORM='UNFORMATTED',FILE='AOONEINT')
#endif
      ENDIF
C
C----------------------------------------------------------------------
C     1. Get ao overlap matrix; unpack to square matrix
C
C        For OCE calculations molecule input is useless, and the
C        overlap matrix has to be read from LUONEL
C
      IF (OCEINT) THEN
         CALL DZERO(WRK(KSAO),N2BASX)
         CALL DZERO(WRK(KSAOU),NNBAST)
         CALL DZERO(WRK(KW1),NNBASX)
         REWIND LUONEL
         CALL MOLLAB('OVERLAP ',LUONEL,LUPRI)
         CALL DZERO(WRK(KSAOU),NNBAST)
   15    READ(LUONEL) SCR,ISCR,LENGTH
         DO 17 I=1,LENGTH
            WRK(KSAOU -1 + ISCR(I)) = SCR(I)
   17    CONTINUE
         IF (LENGTH .GE. 0) GO TO 15
         REWIND LUONEL
C test of AO over...
         IF (IPRSHK .GE. 45) THEN
            DO 1103 ISYM = 1,NSYM
               WRITE(6,'(A,I5)') 'AO small-triangel, sym=',ISYM
               CALL OUTPAK(WRK(KSAOU+IIBAS(ISYM)),NBAS(ISYM),1,6)
 1103       CONTINUE
         ENDIF
C  Fill a big all-symmetry triangel
         CALL DCOPY(NNBAS(1),WRK(KSAOU),1,WRK(KW1),1)
         IF (NSYM .GE. 2) THEN
         DO 18 ISYM = 2,NSYM
            IBASI  = IBAS(ISYM)
            IIBASI = IIBAS(ISYM)
            NBASI  = NBAS(ISYM)
            IBASIX = IBASI*(IBASI+1)/2
            DO 19 I =1,NBASI
               ITRIU = I*(I-1)/2
               ITRIW = IBASIX + I*IBASI + ITRIU
               CALL DCOPY(I,WRK(KSAOU + IIBASI +ITRIU),1,
     *                    WRK(KW1 + ITRIW),1)
   19       CONTINUE
   18    CONTINUE
            IF(IPRSHK .GE. 45) THEN
               WRITE (LUPRI,'(A,A)') ' AO-overlap in triangular',
     *                            ' no-symmetry form'
               CALL OUTPAK(WRK(KW1),NBAST,1,LUPRI)
            ENDIF
         ENDIF
      ELSE
C  here if not OCE calculation
COLD     CALL OVLINT(WRK(KW1),1)
      TIMINT = SECOND()
      CALL OVLMOM(MOMTYP,WRK(KSAOU),1)
      TIMINT = SECOND() - TIMINT
      WRITE (LUPRI,'(///A/)') ' The following properties are evaluated:'
      IF (MOMTYP(1)) WRITE(LUPRI,'(A)')' Overlap '
      IF (IPRSHK .GT. 0) WRITE (LUPRI,'(//A,F10.2,A)')
     1   ' >>> TIME USED FOR AO OVERLAP INTEGRALS:',TIMINT,' SECONDS.'

      ENDIF
      CALL DSPTSI(NBAST,WRK(KSAOU),WRK(KSAO))
C
      IF (IPRSHK .GT. 45) THEN
         WRITE (LUPRI,'(/A)') ' AO overlap matrix read by WESSHK'
         CALL OUTPUT(WRK(KSAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C End Overlap matrix ==========================================
C
C We write the unexpanded final state CI vectors (and their
C MO:s) to LUSPCW file. (Read and expand 1 by 1 later)
C
      IF(.NOT. CICALC) THEN
         CALL WRSPEC(WRK(KCICHH),WRK(KCMOCH),WRK(KEMCCH))
      END IF
C
C
C     3. Calculate overlap between hole/particle final state(s)
C        and the ground state with beta/alfa orbital IACORB
C        annihilated/created
C
C Obtain string information and expand graph for Ground State
C (KCI2EX used as scratch)
C
C In general the two symmetries are different
         LSYM1X = 1
         LSYM2X = 1
         NOSYM  = 1
         IACORB = IORB(IACSYM) + IACORB
C
      CALL GETCIX(WRK(KCINDX),LSYM1,LSYM1,WRK(KW2),LW2,NOSYM)
C
c     write(6,*) ' NCONF1 = ', NCONF1
      CALL WESGCI(1,1,NCONF1,WRK(KCI2EX),0,0,0,DUMMY)
c     write(6,*) ' vector read in from WESGCI '
c     call wrtmat(wrk(kci2ex),1,nconf1,1,nconf1)
c. Remove/add an electron from ground state and expand
c. into expanded graph
C
      CALL EXPCD2(WRK(KCI2EX),WRK(KCI1EX),LSYM1,LSYM1X,LUSYMV(1),
     &            WRK(KCINDX),WRK(KW3),LW3,IACFLG,IACORB)
C
C    SUBROUTINE EXPCD2(CII,CIE,ISMI,ISME,LUSTRI,XNDXCI,
C   &                  IWRK,LWORK,IAC,IACORB)
C
      IF (IPRSHK .GE. 35) THEN
          WRITE (LUPRI,'(//A,I8)')
     *    ' CI vector no. 1 after SHAKE:EXPCD2, length =',NCONFX
            WRITE( LUPRI,'(2(1X,5F15.12/))')
     *      (WRK(KCI1EX-1+I),I=1,NCONFX)
      END IF
C
      EMCGS = EMCWES(1,1)
C
      CALL READMW(IRDMOW,LUIT1V(1,1),WRK(KCMOGS))
      IF (IPRSHK .GT. 0) WRITE (LUPRI,'(/A,F20.12)')
     *   ' GROUND STATE ENERGY           (AU) :',EMCGS
C
C write expanded GS CI vector on temporary file
C        open LUTMPW
C
#if defined (SYS_VAX)
         L0 = NCONFX/128 + 1
         OPEN(UNIT=LUTMPW,STATUS='NEW',FORM='UNFORMATTED',
     *        BLOCKSIZE=L0,INITIALSIZE=L0)
#else
         OPEN(LUTMPW,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *        FILE='fort.52')
#endif
C
      REWIND LUTMPW
         CALL WRITT(LUTMPW,NCONFX,WRK(KCI1EX))
C
C     PREPARE for ISTAT2 loop and WESSWF
C     ----------------------------------
C     define LUSYMB and NCONF
C     prepare for BIORAS in WFOVLP and WESSWF.WFOVLP
C
         NCONF  = NCONF2
C
C   calculate raw moments
C
      IF (CICALC) THEN
C get MO:s if CI
         REWIND LUSPEC(1)
            CALL MOLLAB('NEWORB  ',LUSPEC(1),LUPRI)
            CALL READT(LUSPEC(1),NCMOT,WRK(KCMOCH))
         REWIND LUSPEC(1)
            CALL MOLLAB('STARTVEC',LUSPEC(1),LUPRI)
      ELSE
         REWIND LUSPCW
      ENDIF
C
      DO 300 ISTAT2 = 1,NSTAT2
         IF (IPRSHK .GE. 5) THEN
             WRITE(LUPRI,'(//A,I3/A/)')
     *       ' Raw moment calculation for core hole state no.',ISTAT2,
     *       ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         END IF
         IF (CICALC) THEN
            CALL READT(LUSPEC(1),NCONF2,WRK(KCI1EX))
         ELSE
            CALL READT(LUSPCW,NCONF2,WRK(KCI1EX))
            CALL READT(LUSPCW,NCMOT,WRK(KCMOCH))
         ENDIF
            CALL EXPCID(WRK(KCI1EX),WRK(KCI2EX),LSYM2,1,LUSYMV(2),
     &                  WRK(KCINDX),WRK(KW3),LW3,IPRECV)
C  do we know  NCONFX here ???
            IF (IPRSHK .GE. 35) THEN
               WRITE (LUPRI,'(//A,I8)')
     *         ' CI vector no. 2 after WESECV, length =',NCONFX
               WRITE (LUPRI,'(2(1X,5F15.12/))')
     *         (WRK(KCI2EX-1+I),I=1,NCONFX)
            END IF
C
            NCONF  = NCONFX
C
C We need to have the GS molecular orbital in which 1 electron is
C annihilated/created into the active space.
C How should we procced with this in RAS version? by REORMO or REORD ?
C
CTEMP     IF (SHAKE .AND. ISTAT2 .EQ. 1) THEN
CTEMP       CALL REORMO(WRK(KCMOGS),WRK(KCMOCH),WRK(KSAO),WRK(KW3),LW3)
CTEMP     END IF
C
C        retain expanded CI vector and  make copy of GS orbitals
C         since these are destroyed in WFOVLP.BIORAS
         REWIND LUTMPW
         CALL READT(LUTMPW,NCONF,WRK(KCI1EX))
         CALL DCOPY(NCMOT,WRK(KCMOGS),1,WRK(KCMGS1),1)
C
         RAWMOM = WFOVLP(WRK(KCI1EX),WRK(KCMGS1),WRK(KCI2EX),
     *                   WRK(KCMOCH),WRK(KSAO),WRK(KCINDX),
     *                   WRK(KW3),LW3)
C
C  FUNCTION WFOVLP(CVECA,CMOA,CVECB,CMOB,SAO,INDXCI,WRK,LWRK)
C
         WRK((KRWMOM-1)+ISTAT2) = RAWMOM
C
         IF (IPRSHK .GT. 0) THEN
            WRITE (LUPRI,'(/A,I3,A,F20.12/A,F16.8)')
     *      ' HOLE/PARTICLE STATE NO.',ISTAT2,' ENERGY (AU) :',
     *      WRK(ISTAT2),' OVERLAP  MOMENT FROM GROUND STATE   ',RAWMOM
         END IF
C        get transition energies (we have used KEMCCH = 1)
         IF (WRK(ISTAT2) .NE. D0) WRK(ISTAT2) = WRK(ISTAT2) - EMCGS
  300 CONTINUE
C
C     Close and delete temporary file LUTMPW
C
      IF (NASHT .GT. 0) CLOSE(UNIT=LUTMPW,STATUS='DELETE')
C
C     4. Symmetric orthonormalization of core hole states
C
      IF (SYMORT) THEN
         IF (IPRSHK .GE. 5) THEN
             WRITE(LUPRI,'(//A/A/)')
     *       ' Symmetrical orthonormalization of core hole states',
     *       ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         END IF
C
C
         CALL WESSWF(WRK(KSAO),WRK(KSWF),WRK(KTSWF),WRK(KW4),LW4)
C        CALL WESSWF(SAO,SWF,TSWF,WRK,LWRK)
C
C        Close and delete temporary file LUSPCW
C
         CLOSE(UNIT=LUSPCW,STATUS='DELETE')
C
C
C        Transform moment to symmetrically orthonomalized
C        wave function basis.
C
         CALL MPAB(WRK(KTSWF), NSTAT2,NSTAT2,NSTAT2,NSTAT2,
     *              WRK(KRWMOM),NSTAT2,1,     NSTAT2,1,
     *              WRK(KFIMOM),NSTAT2,1)
      ELSE
         CALL DZERO(WRK(KFIMOM),NSTAT2)
      END IF
C
C
C     Printout of final shake-up energies and intensities.
C
C
      CALL HEADER('<<<<< FINAL RESULTS FROM WESTA SHAKE-UP >>>>>',-1)
      CALL TIMDAT
      WRITE (LUPRI,'(/A,F12.6)')
     *   ' GROUND STATE ENERGY (AU)           :',EMCGS
      IF (WRK(KEMCCH) .NE. D0) THEN
         WRITE (LUPRI,'(/A,F12.6)')
     *   ' LOWEST HOLE/PART. ENERGY (AU)      :',(WRK(KEMCCH)+EMCGS)
         WRITE (LUPRI,'(/A,F12.6)')
     *   ' IONIZATION/AFFINITY POTENTIAL (EV) :',FACEV * WRK(KEMCCH)
      END IF
C
      IF (SYMORT) THEN
         WRITE (LUPRI,'(//A)') ' Spectral overlap matrix:'
         CALL OUTPAK (WRK(KSWF),NSTAT2,1,LUPRI)
         WRITE (LUPRI,'(//A)') ' S-1/2 tranformation:'
         CALL OUTPUT (WRK(KTSWF),1,NSTAT2,1,NSTAT2,NSTAT2,
     *   NSTAT2,1,LUPRI)
      ENDIF
C
      CALL HEADER('<<<<< RESULTS FOR SHAKE-UP SPECTRUM >>>>>',-1)
      WRITE (LUPRI,'(//A,T12,A/T12,A/)') ' STATE',
     *     'SHAKE ENERGY    SHAKE INTENSITY    SHAKE INTENSITY',
     *     '    ( eV )        ( non-ort )        ( sym-ort )'
      DO 1000 ISTAT2 =1,NSTAT2
         SHINT1 = WRK(KRWMOM-1+ISTAT2)
         SHINT1 = SHINT1 * SHINT1
         SHINT2 = WRK(KFIMOM-1+ISTAT2)
         SHINT2 = SHINT2 * SHINT2
         WRITE (LUPRI,'(I5,F16.6,F19.6,F20.6)')
     *          ISTAT2,(FACEV*(WRK(KEMCCH+ISTAT2-1) - WRK(KEMCCH))),
     *          SHINT1,SHINT2
 1000 CONTINUE
C
C  if shake restore to gs values
C
      RETURN
C     END OF WESSHK.
      END
