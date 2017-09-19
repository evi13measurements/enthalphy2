#  /* Deck westa */
      PROGRAM WESTA
C
C (c) Copyright 1985  Hans Agren and Hans Joergen Aa. Jensen
C (c) Copyright 1995  Hans Agren, Hans Joergen Aa. Jensen, Jeppe Olsen
C
C *** This is WESTA, a program for calculation of non-adiabatic
C     coupling elements and transition properties for MCSCF states
C     with non-orthogonal orbitals.
C
C   RELEASE: May-1995 version 2.0
C     November-1985 Version 1.0: Complete Active space wave functions.
C
C   Written by:
C
C     Hans Agren                 <agren@ifm.liu.se>
C
C   Institute of Physics and Measurement Technology
C   University of Linkoeping
C   S-58183 Linkoeping
C   Sweden
C
C     Hans Joergen Aa. Jensen    <hjj@dou.dk>
C
C   Department of Chemistry
C   University of Odense
C   DK-5230 Odense M
C   Denmark
C
C     Jeppe Olsen
C
C   Department of Theoretical Chemistry
C   University of Lund
C   Chemical Center
C   P.O.B. 124
C   S-221 00 LUND
C   SWEDEN
C
C
C
#include "implicit.h"
      PARAMETER (LWORK = 200 000 000)
C
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "wesinf.h"
C
C..   Include block data definition of MULD2H (in SBDORB)
C
      EXTERNAL SBDORB
C
C     GPIO_INI_DALTON: Initialze all unit numbers to neg. number,
C     define default file names and labels
C     SIRINI: Initalize Sirius unit numbers
C
      CALL GPIO_INI_DALTON
      LUPRI = 6
      CALL SIRINI(.TRUE.)
C     CALL SIRINI(NOOPEN)
C
      CALL TIMER('START ',TIME1,DUM)
      TIMIN = TIME1
      CALL WESINP(WORK,LWORK)
      CALL TIMER('WESINP',TIMIN,TIMOUT)
C
      IF (NONADI) THEN
         CALL WESNAD(WORK,LWORK)
         CALL TIMER('WESNAD',TIMIN,TIMOUT)
      END IF
C
      IF (TRAMOM) THEN
         CALL WESMOM(WORK,LWORK)
         CALL TIMER('WESMOM',TIMIN,TIMOUT)
      END IF
C
      IF (SHAKE ) THEN
         CALL WESSHK(WORK,LWORK)
         CALL TIMER('WESSHK',TIMIN,TIMOUT)
      END IF
C
      CALL TIMER('WESTA ',TIME1,TIME2)
C
C TIME AND DATE STAMP.
C
      CALL TIMDAT
      STOP 'End of WESTA calculation.'
      END
#  /* Deck wesinp */
      SUBROUTINE WESINP(WORK,LWORK)
C
C (4-Nov-1985 hjaaj,ha)
C
#include "implicit.h"
      DIMENSION WORK(LWORK)
C
      PARAMETER (D2 = 2.0D0)
      PARAMETER (NDIR = 6, NTABLE = 27)
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     *          REWORD*12, RWORD*6, TITLIN*80
C
      CHARACTER*72 STARS,RELEAS
      INTEGER LUMLC(2)
      LOGICAL INPTES, RESTAR, RDDELQ
C
C
C Used from common:
C   INFINP : LSYM, ISPIN, SPIN, NELMN*, NAS*
C   INFORB : NSYM,NISH(8),NASH(8),NSSH(8),NBAS(8)
C   WESINF : NSTAT2, NMAXD,DELQ(3,100)
C   WESTAP : MSTAT2
C
#include "priunit.h"
#include "maxorb.h"
#include "wesinf.h"
#include "westap.h"
#include "infinp.h"
#include "inforb.h"
#include "infpri.h"
#include "codata.h"
C
      DATA TABDIR /'*END OF', '*READIN', '*NUCREP',
     *             '*ONEINT', '*DIPCTL', '**END O'/
      DATA TABLE  /'.PRINT ', '.NONADI', '.DIPLEN',
     *             '.DIPVEL', '.DIPMAG', '.MOLECU',
     *             '.INPTES', '.READMO', '.RESTAR',
     *             '.SYMORT', '.WFOVER', '.DISPLA',
     *             '.STATES', '.CENTOP', '.SYMDIF',
     *             '.SHAKE ', '.NSTATE', '.REDUCE',
     *             '.FROZEN', '.OCEINT', '.NACTEL',
     *             '.CICALC', '.CIORD ', '.TITLE ',
     &             '.ADIABC', '.MASSES', '.AADISP'/
C
      DATA ISTAT1/1/, ISTAT2/2/, RDDELQ/.TRUE./,
     *     LUMLC/5,0/, INPTES/.FALSE./, RESTAR/.FALSE./
C
      DATA STARS(1:36)/'************************************'/
      DATA STARS(37:72)/'************************************'/
      DATA RELEAS(1:36)/'*WESTA* a program for transition pro'/
      DATA RELEAS(37:72)/'perties and non-adiabatic coupling '/
C
C     Initialize common block WESINP with default values
C     for the various options in block data WBDINF.
C
      EXTERNAL WBDINF
C
C     ***** WRITE PROGRAM IDENTIFICATION *****
C
      WRITE(6,2200) STARS,RELEAS,STARS
 2200 FORMAT(1X,A)
      WRITE(6,400)
  400 FORMAT(//,T2,'Authors:',
     * T12,'Hans Agren,              University of Linkoeping, Sweden'/
     * T12,'Hans Joergen Aa. Jensen, University of Odense,     Denmark'/
     * T12,'Jeppe Olsen,             University of Lund,       Sweden'/)
C
      CALL TSTAMP(' ',6)
C
      CALL SETUPW(1,WORK,LWORK)
      CALL HEADER('OUTPUT FROM WESTA INPUT PROCESSING',-1)
C
C     ***** PROCESS INPUT FOR COMMON  /WESINF/  *****
C
   10 READ (5, '(A7)',END=8000) WORD
      IF (WORD .NE. '*WESTA ' .AND. WORD .NE. '**WESTA') THEN
         GO TO 10
      ELSE
  100    READ (5, '(A7)') WORD
 1100    PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '.') THEN
            DO 1000 I = 1, NTABLE
               IF (TABLE(I) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111,112,113,114,115,116,117,118,119,120,
     *                   121,122,123,124,125,126,127), I
               END IF
 1000       CONTINUE
               WRITE (6,'(/,3A,/)') ' KEYWORD ',WORD,' NOT RECOGNIZED',
     *                              ' IN WESINP.'
               STOP ' ILLEGAL KEYWORD IN WESINP '
  101       CONTINUE
C           '.PRINT ' - print levels
C                      (IPRDEF will be default for levels not specified)
               READ (5, '(8I5)') IPRDEF,IPRNAD,IPRTRM,IPRBIO,
     *                           IPRMLC,IPRSHK
               GO TO 100
  102       CONTINUE
C           '.NONADI' - calculate non-adiabatic coupling
               NONADI = .TRUE.
               NGEOM = 2
               NSTAT2 = 1
C.. istate og lignende ..
               GO TO 100
  103       CONTINUE
C           '.DIPLEN' - calculate dipole length moments
               TRAMOM    = .TRUE.
               MOMTYP(2) = .TRUE.
               GO TO 100
  104       CONTINUE
C           '.DIPVEL' - calculate dipole velocity moments
               TRAMOM    = .TRUE.
               MOMTYP(3) = .TRUE.
               GO TO 100
  105       CONTINUE
C           '.DIPMAG' - calculate magnetic dipole moments
               TRAMOM    = .TRUE.
               MOMTYP(4) = .TRUE.
               GO TO 100
  106       CONTINUE
C           '.MOLECU' - units for molecule input
               READ (5, '(2I5)') LUMLC
               GO TO 100
  107       CONTINUE
C           '.INPTES' - only input test, no calculation
               INPTES = .TRUE.
               GO TO 100
  108       CONTINUE
C           '.READMO' - input source for molecular orbitals,
C                       using READMW
              READ (5, '(I5)') IRDMOW
              IF (IRDMOW .GT. 2) THEN
                READ (5, '(8I5)') NSYM, NELMN1,NELMX1,NELMN3, NELMX3
                READ (5, '(8I5)') (NISH(I),I=1,NSYM)
                IF (NELMN1 .EQ. 0 .AND. NELMX3 .EQ. 0) THEN
C                    ... CAS calculation
                  MCTYPE = 1
                  READ (5, '(8I5)') (NASH(I),I=1,NSYM)
                  READ (5, '(8I5)') (NSSH(I),I=1,NSYM)
                  READ (5, '(8I5)') (NBAS(I),I=1,NSYM)
                  NASHT = 0
                  DO 401 I = 1,NSYM
                     NAS1(I) = 0
                     NAS2(I) = NASH(I)
                     NAS3(I) = 0
                     NASHT   = NASHT + NASH(I)
  401             CONTINUE
                ELSE
C                 ... RAS calculation
                  MCTYPE = 2
                  READ (5, '(8I5)') (NAS1(I),I=1,NSYM)
                  READ (5, '(8I5)') (NAS2(I),I=1,NSYM)
                  READ (5, '(8I5)') (NAS3(I),I=1,NSYM)
                  READ (5, '(8I5)') (NSSH(I),I=1,NSYM)
                  READ (5, '(8I5)') (NBAS(I),I=1,NSYM)
                  NASHT = 0
                  DO 411 I = 1,NSYM
                     NASH(I) = NAS1(I) + NAS2(I) + NAS3(I)
                     NASHT   = NASHT   + NASH(I)
  411             CONTINUE
                  IF (NELMN1 .LT. 0) THEN
C                    ... special test option for CAS in RAS2
                     NELMN1 = 0
                     NELMX3 = 0
                  END IF
               END IF
             END IF
               GO TO 100
  109       CONTINUE
C           '.RESTAR' - restart, not implemented.
                  RESTAR = .TRUE.
C                 READ (5,'(A12)') REWORD
C                 IF (REWORD(1:6) .NE. '.FROM ') THEN
C                    WRITE(6,'(4A)') '"',REWORD(1:6),'" IS INCORRECT ',
C    *                     'INPUT, ".FROM " EXPECTED AFTER ".RESTAR"'
C                    STOP ' INPUT ERROR IN WESINP'
C                 ELSE
C                    RWORD = REWORD(7:12)
C                 END IF
               GO TO 100
  110       CONTINUE
C           '.SYMORT' - symmetric orthogonalization of wf(1) and wf(2)
C                       for the two geometries when non-adiabatic.
               SYMORT = .TRUE.
               WAVEOV = .TRUE.
               GO TO 100
  111       CONTINUE
C           '.WFOVER' - calculate overlap between wf(1) and wf(2)
               WAVEOV    = .TRUE.
               MOMTYP(1) = .TRUE.
               GO TO 100
  112       CONTINUE
C           '.DISPLA' - geometry displacement
               RDDELQ = .FALSE.
c              READ (5, '(I5)') NMAXD
               NMAXD = 2
               READ (5,*) ((DELQ(J,I), J = 1, 3), I = 1,NMAXD)
               DISPL = DNRM2(3*NMAXD,DELQ,1)
               GO TO 100
  113       CONTINUE
C           '.STATES' - states
               READ(5,'(2I5)') ISTAT1,ISTAT2
               GO TO 100
  114       CONTINUE
C           '.CENTOP' - Specification of origin for transition moments
               READ (5, '(3F10.6)') (CENTOP(I),I=1,3)
               GO TO 100
  115       CONTINUE
C           '.SYMDIF' - Two states of different state
C                        and/or spin symmetry
               READ (5, '(2I5)') LSYM1,LSYM2
               READ (5, '(2I5)') LSPIN1,LSPIN2
               GO TO 100
  116       CONTINUE
C           '.SHAKE ' - Shake-up calculation between a singlet
C                       ground state and doublet core hole states
C                       IACFLG, IACSYM, IACORB give type, symmetry and
C                       orbital # for annihilation/creation.
               READ (5, '(I5)') NSTAT2
               READ (5, '(3I5)') IACFLG, IACSYM, IACORB
               SHAKE  = .TRUE.
               LSPIN1 = 1
               LSPIN2 = 2
               NGEOM  = 1
               DO I = 1, NSTAT2
                 ICIORD(I) = I
               END DO
               GO TO 100
  117       CONTINUE
C           '.NSTATE' - Number of states for which we simultaneously
C                       want to obtain intensities.
               READ (5, '(I5)') NSTAT2
               DO 1116 I =1, NSTAT2
                   ICIORD(I) = I
 1116          CONTINUE
               GO TO 100
  118       CONTINUE
C           '.REDUCE' - Transition moment between wave functions with
C                       reduced graphs (e.g. X-ray moments)
               REDUCE = .TRUE.
               GO TO 100
  119       CONTINUE
C           '.FROZEN' - Core orbital is frozen
               FROZEN = .TRUE.
               GO TO 100
  120       CONTINUE
C           '.OCEINT' - OCE-integrals used, read AO-overlap from file
               OCEINT = .TRUE.
               GO TO 100
  121       CONTINUE
C           '.NACTEL' - Number of active electrons
               READ (5, '(I5)') NACTEL
               GO TO 100
  122       CONTINUE
C           '.CICALC' - CI-vectors (and MO:s) from CI calculation.
C                       Assign CI-F21 file to LUSPEC(1)
               CICALC = .TRUE.
               GO TO 100
  123       CONTINUE
C           '.CIORD ' - Order of "ISTATE" in each NSTATE CI- f21 file
C                       Assign CI-F21 file to LUSPEC(1)
               READ (5, '(12I5)') (ICIORD(I), I = 1, NSTAT2)
               GO TO 100
  124       CONTINUE
C           '.TITLE ' - Title lines
               WRITE (6,'(///A/A/)')
     1            ' WESTA JOB TITLE:',' ----------------'
 1115          READ (5, '(A)') TITLIN
               WORD = TITLIN(1:7)
               IF (WORD(1:1) .EQ. '*' .OR. WORD(1:1) .EQ. '.') THEN
                  WRITE (6,'(/)')
                  GO TO 1100
               END IF
               WRITE (6,'(1X,A)') TITLIN
               GO TO 1115
  125       CONTINUE
C           '.ADIABC - Adiabatic correction to Born-Oppenheimer
               NONADI = .TRUE.
               NGEOM = 2
               NSTAT2 = 0
               GO TO 100
  126       CONTINUE
C           '.MASSES - enter masses of the two nuclei for .ADIABC
               READ(5,*) NUCMASS(1), NUCMASS(2)
               GO TO 100
  127       CONTINUE
C           '.AADISPL' - geometry displacement in Aangstroem
               RDDELQ = .FALSE.
c              READ (5, '(I5)') NMAXD
               NMAXD = 2
               READ (5,*) ((DELQ(J,I), J = 1, 3), I = 1,NMAXD)
               FAC = 1.0D0/XTANG
               CALL DSCAL(3*NMAXD, FAC, DELQ, 1)
               DISPL = DNRM2(3*NMAXD,DELQ,1)
               GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 130
         ELSE
            WRITE (6,'(/,3A,/)') ' PROMPTER "',PROMPT,'" ILLEGAL'
            STOP ' ILLEGAL PROMPT IN WESINP '
         END IF
      END IF
  130 CONTINUE
C
C     ===============
C
      IF (LSYM2  .EQ. 0) LSYM2  = LSYM1
      IF (LSPIN2 .EQ. 0) LSPIN2 = LSPIN1
      SPIN1  = (LSPIN1 -1)/D2
      SPIN2  = (LSPIN2 -1)/D2
      LSYM   = LSYM1
      ISPIN  = LSPIN1
      SPIN   = SPIN1
C
      IF (WAVEOV .AND. .NOT. NONADI .AND. .NOT. SHAKE) TRAMOM = .TRUE.
C
C
C
      IF (INPTES) 
     &WRITE (LUPRI,'(//A/)') '** THIS IS AN INPUT TEST RUN ONLY. **'
c     WRITE (LUPRI,'(2A/)') ' THE FOLLOWING PROPERTIES WILL BE',
c    *                   ' CALCULATED IN THIS RUN:'
      CALL HEADER(
     &  'THE FOLLOWING PROPERTIES WILL BE CALCULATED IN THIS RUN:',-1)
      IF (NONADI) THEN
       IF (NSTAT2 .EQ. 0) THEN
          WRITE (LUPRI,'(6X,A)')
     &    'Adiabatic correction to reference state.'
       ELSE IF (NSTAT2 .EQ. 1) THEN
         WRITE (LUPRI,'(6X,A)') 'Non-adiabatic coupling.'
         IF (SYMORT) WRITE (LUPRI,'(9X,A)')
     *      ' The two MCSCF wave functions are symmetrically'//
     *      ' orthogonalized.'
       ELSE
         WRITE (LUPRI,'(6X,A,I5,A//A)')
     &      'Non-adiabatic coupling with',NSTAT2,' states',
     &      'ERROR : THIS IS NOT IMPLEMENTED.'
       END IF
         WRITE(LUPRI,'(/6X,A,2F12.6)') 'Nuclear masses (u):',
     &     NUCMASS(1), NUCMASS(2)
         IF (RDDELQ) THEN
            WRITE(LUPRI,'(/A)') ' INPUT ERROR: Displacement '//
     &      'must be specified in input with .DISPLA'
            STOP 'Displacement must be specified in input'
         END IF
         WRITE(LUPRI,'(//A,F12.8,A,F12.8,A)')
     &     ' Total displacement A to B:',DISPL,' au,',DISPL*XTANG,' AA'
         WRITE(LUPRI,'(/A,2(/I5,3F12.8))')
     &     ' Displacement coordinates (au):',
     &     (J, (DELQ(I,J), I=1,3), J=1,NMAXD)
      END IF
      IF (TRAMOM) THEN
         WRITE (LUPRI,'(/6X,A)') 'Transition properties:'
         IF (MOMTYP(1)) WRITE (LUPRI,'(9X,A)') 'Wave function overlap.'
         IF (MOMTYP(2)) WRITE (LUPRI,'(9X,A)')
     *      'Dipole length transition moments.'
         IF (MOMTYP(3)) WRITE (LUPRI,'(9X,A)')
     *      'Dipole velocity transition moments.'
         IF (MOMTYP(4)) WRITE (LUPRI,'(9X,A)')
     *      'Magnetic dipole transition moments.'
      END IF
      IF (SHAKE) THEN
         WRITE (LUPRI, '(/A/,A,I3,A)')
     *    ' Annihilation/creation amplitudes between the ground state',
     *    ' and the',NSTAT2,' lowest hole/particle states.'
         IF (IACFLG .EQ. 0 ) THEN
            WRITE(LUPRI,'(/A,/)')
     *      ' No Annihilation/creation, IACFLG =0'
         ELSE IF (IACFLG .EQ. 1 ) THEN
            WRITE(LUPRI,'(/A,I3/A,I3/)')
     *      ' SHAKE-annihilation of one electron in symmetry=', IACSYM,
     *      ' orbital number ', IACORB
         ELSE IF (IACFLG .EQ. 2 ) THEN
            WRITE(LUPRI,'(/A,I3/A,I3/)')
     *      ' BIS-creation of one electron in symmetry=', IACSYM,
     *      ' orbital number ', IACORB
         ELSE
            WRITE(LUPRI,'(/A/)')
     *      ' IACFLG out of RANGE (0,1,2) '
            STOP 'IACFLG outof Range (0,1,2) '
         END IF
      END IF
C
         IF (SYMORT) WRITE (LUPRI,'(/A)')
     *     ' The hole/particle MCSCF wave functions are symmetrically'//
     *     ' orthogonalized.'
         IF (NSTAT2 .GT. MSTAT2) THEN
            WRITE (LUPRI, '(//A,I3)') ' ERROR, TOO MANY FINAL STATES,'//
     *         ' CURRENT MAXIMUM IS',MSTAT2
            STOP 'WESINP: too many final states'
         END IF
          WRITE(LUPRI,'(/10X,A)') 'Order of CI ISTATE vectors per state'
          WRITE(LUPRI,'(9X,(12I5))') (ICIORD(I), I = 1, NSTAT2)
C
      IF (IPRNAD .EQ. -1) IPRNAD = IPRDEF
      IF (IPRTRM .EQ. -1) IPRTRM = IPRDEF
      IF (IPRBIO .EQ. -1) IPRBIO = IPRDEF
      IF (IPRMLC .EQ. -1) IPRMLC = IPRDEF
      IF (IPRSHK .EQ. -1) IPRSHK = IPRDEF
      WRITE (LUPRI,'(//,(A,I5))')
     *   ' GENERAL PRINT LEVEL: ',IPRDEF,
     *   ' WESNAD  PRINT LEVEL: ',IPRNAD,
     *   ' WESMOM  PRINT LEVEL: ',IPRTRM,
     *   ' BIORAS  PRINT LEVEL: ',IPRBIO,
     *   ' READIW  PRINT LEVEL: ',IPRMLC,
     *   ' WESSHK  PRINT LEVEL: ',IPRSHK
C
C     **** PROCESS INPUT FOR VARIOUS PROGRAM SECTIONS  *****
C
  200 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
           IF (WORD .EQ. TABDIR(I)) THEN
C             GO TO (1,2,3,4,5), I
              GO TO (1,2,2,2,2,1), I
           END IF
  210    CONTINUE
         WRITE (LUPRI,'(/,3A,/)') ' DIRECTORY ',WORD,' NONEXISTENT.'
         STOP ' ILLEGAL DIRECTORY IN WESINP '
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' PROMPTER "',PROMPT,'" ILLEGAL OR',
     *                        ' OUT OF ORDER.'
         STOP ' PROGRAM STOPPED IN WESINP '
      END IF
    2   CONTINUE
C   2   CALL REAINP(WORD)
      GO TO 200
C   3   CALL NUCINP(WORD)
C     GO TO 200
C     etc. 4 to 14
C
    1 CONTINUE
C
#if defined (VAR_READIW)
#if !defined (SYS_VAX)
C
C
      IF (LUMLC(1) .NE. LUINP) THEN
         LUOPEN = LUMLC(1)
         OPEN(LUMLC(1),STATUS='OLD',FORM='FORMATTED',
     *        FILE='MOLGEOA.IN',ERR=9000)
      END IF
      IF (LUMLC(2) .GT. 0 .AND. LUMLC(2) .NE. LUINP .AND.
     *    LUMLC(2) .NE. LUMLC(1)) THEN
         LUOPEN = LUMLC(2)
         OPEN(LUMLC(2),STATUS='OLD',FORM='FORMATTED',
     *        FILE='MOLGEOB.IN',ERR=9000)
      END IF
C
#endif
C
      CALL READIW(LUMLC,RDDELQ)
#endif
      IF (.NOT. INPTES) CALL SETUPW(2,WORK,LWORK)
C
C
      ISTATV(1,1) = ISTAT1
      ISTATV(1,2) = ISTAT2
      ISTATV(2,1) = ISTAT1
      ISTATV(2,2) = ISTAT2
C
      IF (INPTES) STOP 'WESTA INPUT TEST FINISHED.'
C
      CLOSE(UNIT=5)
#if defined (VAR_READIW)
      IF (LUMLC(1) .NE. 5) CLOSE(UNIT=LUMLC(1))
      IF (LUMLC(2) .GT. 0 .AND. LUMLCL(2) .NE. 5)
     &   CLOSE(UNIT=LUMLC(2))
#endif
      RETURN
C
 8000 CONTINUE
      WRITE (LUPRI, '(///A/A)')
     *   ' *** ERROR, WESTA input not found',
     *   ' *** WESTA INPUT MUST BEGIN WITH LABEL "**WESTA"'
      STOP 'WESINP ERROR, input label "**WESTA" not found.'
C
#if defined (VAR_READIW)
 9000 CONTINUE
      WRITE (LUPRI,'(///A,I5)') ' WESINP: ERROR OPENING UNIT',LUOPEN
      STOP 'ERROR OPENING MOLECULE INPUT FILE IN WESINP.'
#endif
C
C END OF WESINP.
      END
#  /* Deck setupw */
      SUBROUTINE SETUPW(ISETUP,WORK,LWORK)
C
C
C (Based on SIRIUS-SETUP Hans Agren Dec. 1983)
C
C  PART 1: Current program limits
C  PART 2: Input/output unit numbers
C  PART 3: Orbital and basis info
C  PART 4: Opening of files
C  PART 5: Offset info
C  PART 6: Index arrays
C  PART 7: CI parameters
c  PART 8: Miscellaneous
C
C CONVENTIONS FOR ARRAYS ARE IDENTICAL TO SIRIUS'
C
C
#include "implicit.h"
      DIMENSION WORK(LWORK)
#include "iratdef.h"
C -- local constants
      PARAMETER (D1=1.0D0, D2=2.0D0)
C
C  WESINF  ? NGEOM,IRDMOW,NCONF1,NCONF2,LSYM1,LSYM2,LSPIN1,LSPIN2,NSTAT2
C  WESTAP  ? LUSYMV,...
C  INFTAP  : LUSIFC
C  INFINP  : FLAG(26)
C  INFORB  ? N2BASX...
C  INFIND  ? IROW,ISAO,ISW...
C  INFVAR  ? NCONF
C  INFDIM  ? SETCI parameters
C
#include "priunit.h"
#include "maxash.h"
#include "maxorb.h"
#include "wesinf.h"
#include "westap.h"
#include "inftap.h"
#include "infinp.h"
#include "r12int.h"
#include "inforb.h"
#include "infind.h"
#include "infvar.h"
#include "infpri.h"
#include "infdim.h"
C
C-- local variables
      LOGICAL NXTLAB
      CHARACTER FORIPV*11, RTNLBL(2)*8
C
      IF (ISETUP .NE. 1) GO TO 20000
C
C  WESTAP:
C  +++++++
C  LUSYMV(1) 10   RAS CI information state 1 (written by SIRIUS)
C  LUSYMV(2) 11   RAS CI information state 2 (written by SIRIUS)
C  LUIT1V(1,1)   41   MSCSF LUIT1 file for state 1, geometry A
C  LUIT1V(2,1)   42       - - " - -              1           B
C  LUIT1V(1,2)   43       - - " - -              2           A
C  LUIT1V(2,2)   44       - - " - -              2           B
C  LUPROPV(1)    45   AOPROPER file for geometry A
C  LUPROPV(2)    46   AOPROPER file for geometry B
C
C  LUSPCW        51  Big file collecting vectors from spectrum
C  LUTMPW        52  Temporary file in WESMOM.
C  LUSPEC(I), I=1,MSTAT2 (12) units 61-72.  files containing
C             vectors for spectrum (shake and tramom)
C
      LUSYMV(1) = 10
      LUSYMV(2) = 11
C
      LUIT1V(1,1) = 41
      LUIT1V(2,1) = 42
      LUIT1V(1,2) = 43
      LUIT1V(2,2) = 44
C
      LUPROPV(1)  = 45
      LUPROPV(2)  = 46
C
      LUSPCW    = 51
      LUTMPW    = 52
      DO 1510 I = 1,MSTAT2
         LUSPEC(I) = 60 + I
 1510 CONTINUE
C
C     OPEN (UNIT=5,STATUS='OLD',ERR=9000)
      LUERR = 0
      LUW4 = 6
      LUPRI = 6
C
C     Read info from LUSIFC for geom A, state 1
C     (to define NBAS, NASH, NISH, etc.)
C
      OPEN(LUIT1V(1,1),STATUS='OLD',FORM='UNFORMATTED',
     *     FILE='SIRIFC.1A')
      LUSIFC = LUIT1V(1,1)
      IF (.NOT. NXTLAB(LBSIFC,RTNLBL,LUSIFC) ) THEN
         CALL QUIT('No labels on SIRIFC.1A')
      ELSE
         WRITE (LUPRI,*) 'Main label on SIRIFC.1A : ',LBSIFC
      END IF
      CALL SETSIR(WORK,LWORK)
      FLAG(26) = .FALSE.
C     ... do not write SIRIUS.STRINGINFO for WESTA (!)
C
C
C ********************************************************
C   SETUPW MAJOR SECTION 2.
C ********************************************************
C
C
20000 IF (ISETUP .NE. 2) GO TO 30000
      FORIPV = 'FORMATTED  '
      IF (IRDMOW .NE. 3 .AND. IRDMOW .NE. 4) FORIPV = 'UNFORMATTED'
      DST = SECOND( )
C
C PART 3: ***** SET UP ORBITAL DATA *********************************
C
C NBAS(*) and NORB(*) are input constructed in subroutine READIN.
C
      SUPSYM = .FALSE.

      NCMOT=0
      NISHT=0
      NASHT=0
      NSSHT=0
      NBAST=0
      NNORBT=0
      NNBAST=0
      N2ORBT=0
      N2BAST=0
      DO 100 ISYM=1,NSYM
         NRXR12(ISYM) = 0
         NOCC(ISYM)=NISH(ISYM)+NASH(ISYM)
         NORB(ISYM)=NSSH(ISYM)+NOCC(ISYM)
         NISHT=NISHT+NISH(ISYM)
         NASHT=NASHT+NASH(ISYM)
         NSSHT=NSSHT+NSSH(ISYM)
         NBAST=NBAST+NBAS(ISYM)
         NNORB(ISYM)=NORB(ISYM)*(NORB(ISYM)+1)/2
         NNORBT=NNORBT+NNORB(ISYM)
         N2ORB(ISYM)=NORB(ISYM)*NORB(ISYM)
         N2ORBT=N2ORBT+N2ORB(ISYM)
         NNBAS(ISYM)=NBAS(ISYM)*(NBAS(ISYM)+1)/2
         NNBAST=NNBAST+NNBAS(ISYM)
         N2BAS(ISYM)=NBAS(ISYM)*NBAS(ISYM)
         N2BAST=N2BAST+N2BAS(ISYM)
         NCMOT=NCMOT+NORB(ISYM)*NBAS(ISYM)
  100 CONTINUE
      NOCCT = NISHT + NASHT
      NORBT = NOCCT + NSSHT
      N2ISHX = NISHT*NISHT
      NNASHX = NASHT*(NASHT+1)/2
      N2ASHX = NASHT*NASHT
      NNASHY = NNASHX*(NNASHX+1)/2
      NNOCCX = NOCCT*(NOCCT+1)/2
      NNORBX=NORBT*(NORBT+1)/2
      N2ORBX=NORBT*NORBT
      NNBASX=NBAST*(NBAST+1)/2
      N2BASX=NBAST*NBAST
      N2OCCX=NOCCT*NOCCT
C
      NISHMA=0
      NASHMA=0
      NSSHMA=0
      NORBMA=0
      NBASMA=0
      NNBASM=0
      N2BASM=0
      DO 200 ISYM=1,NSYM
        NISHMA=MAX(NISHMA,NISH(ISYM))
        NASHMA=MAX(NASHMA,NASH(ISYM))
        NSSHMA=MAX(NSSHMA,NSSH(ISYM))
        NORBMA=MAX(NORBMA,NORB(ISYM))
        NBASMA=MAX(NBASMA,NBAS(ISYM))
        NNBASM=MAX(NNBASM,NNBAS(ISYM))
        N2BASM=MAX(N2BASM,N2BAS(ISYM))
  200 CONTINUE
C
C PART 4: *** Open files 
C
C
      IF (NASHT.GT.1 .AND. (TRAMOM .OR. SHAKE)) THEN
         OPEN(LUSYMV(1),STATUS='OLD',FORM='UNFORMATTED',
     *        FILE='SIRIUS.STRINGINFO.1')
         OPEN(LUSYMV(2),STATUS='OLD',FORM='UNFORMATTED',
     *        FILE='SIRIUS.STRINGINFO.2')
      ELSE
         LUSYMV(1) = -1
         LUSYMV(2) = -1
      END IF
C
C PART 5: *****  OFFSETS *********************************************
C
      IBAS(1)=0
      IIBAS(1)=0
      IORB(1)=0
      IIORB(1)=0
      I2ORB(1)=0
      IOCC(1)=0
      ICMO(1)=0
      IIISH(1)=0
      IIASH(1)=0
      IISH(1)=0
      IASH(1)=0
      ISSH(1)=0
      DO 300 ISYM=1,NSYM-1
         IBAS(ISYM+1)=IBAS(ISYM)+NBAS(ISYM)
         IIBAS(ISYM+1)=IIBAS(ISYM)+NBAS(ISYM)*(NBAS(ISYM)+1)/2
         IORB(ISYM+1)=IORB(ISYM)+NORB(ISYM)
         IIORB(ISYM+1)=IIORB(ISYM)+NORB(ISYM)*(NORB(ISYM)+1)/2
         I2ORB(ISYM+1)=I2ORB(ISYM)+NORB(ISYM)*NORB(ISYM)
         IOCC(ISYM+1)=IOCC(ISYM)+NOCC(ISYM)
         ICMO(ISYM+1)=ICMO(ISYM)+NORB(ISYM)*NBAS(ISYM)
         IIISH(ISYM+1)=IIISH(ISYM)+NISH(ISYM)*(NISH(ISYM)+1)/2
         IIASH(ISYM+1)=IIASH(ISYM)+NASH(ISYM)*(NASH(ISYM)+1)/2
         IISH(ISYM+1)=IISH(ISYM)+NISH(ISYM)
         IASH(ISYM+1)=IASH(ISYM)+NASH(ISYM)
         ISSH(ISYM+1)=ISSH(ISYM)+NSSH(ISYM)
  300 CONTINUE
C
C PART 6: **** INDEX ARRAYS ********************************************
C
C They are:
C  IROW ISW ISX ISMO ISAO NSM ICH  IOBTYP (but not IWOP)
C
C loop to set up triangular indexing array ...
      IROW(1) = 0
         DO 400 I = 1,LIROW
            IROW(I+1) = IROW(I) + I
  400    CONTINUE
C
C setup of IOBTYP (1=frozen, 2=inactive, 3=active, 4=secondary)
C
      INDX=0
      DO 500 ISYM=1,NSYM
C .. MAERKE no frozen in this version /950518-hjaaj
       NFRO(ISYM) = 0
       DO 600 NF=1,NFRO(ISYM)
       INDX=INDX+1
  600  IOBTYP(INDX)=JTFRO
       DO 700 NI=NFRO(ISYM)+1,NISH(ISYM)
       INDX=INDX+1
  700  IOBTYP(INDX)=JTINAC
       DO 800 NT=NISH(ISYM)+1,NOCC(ISYM)
       INDX=INDX+1
  800  IOBTYP(INDX)=JTACT
       DO 900 NA=NOCC(ISYM)+1,NORB(ISYM)
       INDX=INDX+1
  900  IOBTYP(INDX)=JTSEC
  500 CONTINUE
C setup of ISAO and ISMO, symmetry index of basis fu. and orbitals
      II=0
      JJ=0
         DO 1000 ISYM=1,NSYM
         NBASI=NBAS(ISYM)
         NORBI=NORB(ISYM)
         IF(NBASI.EQ.0) GO TO 1000
            DO 1010 I=1,NBASI
            II=II+1
            ISAO(II)=ISYM
            IF(I.GT.NORBI)GO TO 1010
            JJ=JJ+1
            ISMO(JJ)=ISYM
 1010       CONTINUE
 1000    CONTINUE
C
C  Compute reordering indices for molecular orbitals.
C  ISW reorders to inactive-active-secondary ordering
C  ISX reorders back
C  ISWR, ISXR the same thing for RAS wave functions
C
      CALL MAKISW
      CALL MKISWR
C
C ICH and NSM from CICTL:
C
      IN=0
      IU=0
      IR=0
      DO 1130 I=1,NSYM
      NISHI=NISH(I)
      NASHI=NASH(I)
      NSSHI=NSSH(I)
      IF(NISHI.EQ.0)GO TO 1140
      DO 1150 J=1,NISHI
      IN=IN+1
      IR=IR-1
      ICH(IN)=IR
 1150 CONTINUE
 1140 IF(NASHI.EQ.0)GO TO 1160
      DO 1170 J=1,NASHI
      IN=IN+1
      IU=IU+1
      ICH(IN)=IU
      NSM(IU)=I
 1170 CONTINUE
 1160 IF(NSSHI.EQ.0)GO TO 1130
      DO 1180 J=1,NSSHI
      IN=IN+1
      ICH(IN)=0
 1180 CONTINUE
 1130 CONTINUE
C
C >>> no orbital rotations for WESTA
C
      NWOPT  = 0
C
C FOR DIMENSIONING:
C
      NISHDI = MAX(NISHT,1)
      NASHDI = MAX(NASHT,1)
      NWOPMA = 0
      NWOPDI = 1
C
C PART 7: **** GUGA PARAMETERS **************************************
C
C read NCONF1,2 from LUSYMV(1) and (2)
C
      IF (LUSYMV(1) .GT. 0) THEN
         REWIND LUSYMV(1)
         CALL MOLLAB('CISTRING',LUSYMV(1),LUPRI)
         READ (LUSYMV(1)) NCONF1
C
         IF (LSYM1 .EQ. LSYM2)  THEN
           LSYM = LSYM1
           NOSYM = 0
         ELSE
           LSYM = 1
           NOSYM = 1
         ENDIF
         CALL SETCI(NCONFX,NCDETS,LSYM,WORK,LWORK,NOSYM)
         NCONF = NCONFX
C
         REWIND LUSYMV(2)
         CALL MOLLAB('CISTRING',LUSYMV(2),LUPRI)
         READ (LUSYMV(2)) NCONF2

         IF (NCONF .LE. NCONF1 .OR. NCONF .LE. NCONF2) THEN
            WRITE (LUPRI,'(/A/A,I12/A,2I12/A,4I5/A,8I5)') 'FATAL ERROR',
     &      ' NCONF from SETCI',NCONF,
     &      ' NCONF1/2 from STRINGINFO_1/2',NCONF1,NCONF2,
     &      ' LSYM1, LSYM2, LSYM, NOSYM',LSYM1,LSYM2,LSYNM,NOSYM,
     &      ' NASH : ', (NASH(I), I = 1, NSYM)
            CALL QUIT('Fatal error in sirius-westa interface')
         END IF
C
C     ELSE: NCONF from SETSIR call
      ENDIF
C
C
C
      L0 = 8*MAX(NCMOT,NCONF,1024)
C
C     SIRIFC.1A opened above, must always be present
C     OPEN(LUIT1V(1,1),STATUS='OLD',FORM='UNFORMATTED',
C    *     FILE='SIRIFC.1A')
C
      OPEN(LUPROPV(1),STATUS='OLD',FORM='UNFORMATTED',
     *     FILE='AOPROPER.A')
      IF (NSTAT2 .EQ. 1) THEN
         OPEN(LUIT1V(1,2),STATUS='OLD',FORM='UNFORMATTED',
     *        FILE='SIRIFC.2A')
      END IF
      IF (NGEOM .GT. 1) THEN
         OPEN(LUPROPV(2),STATUS='OLD',FORM='UNFORMATTED',
     *        FILE='AOPROPER.B')
         OPEN(LUIT1V(2,1),STATUS='OLD',FORM='UNFORMATTED',
     *        FILE='SIRIFC.1B')
         IF (NSTAT2 .EQ. 1) THEN
            OPEN(LUIT1V(2,2),STATUS='OLD',FORM='UNFORMATTED',
     *           FILE='SIRIFC.2B')
         END IF
      END IF
C
C PART 8: **** MISCELLANEOUS ******************************************
C
C
C Test for inconsistencies:
C
      NUMERR = 0
      IF (NCONF.LE.0) THEN
        WRITE (LUPRI,1610) NCONF
        NUMERR = NUMERR + 1
      ENDIF
      IF (NUMERR .GT. 0) THEN
        CALL TRACE
        STOP 'WESTA error: no configurations'
      ENDIF
 1610 FORMAT(/' SETUPW-ERROR, NCONF =',I10)
C
C *** DEFINE SOME VARIABLES NORMALLY DEFINED IN SIRINP
C     P4FLAG(3),P6FLAG(10) for CICTL (print residuals and energies)
C
      P4FLAG(3)  = .TRUE.
C     P6FLAG(10) = .TRUE.
      DTIM = SECOND() - DST
      IF (IPRDEF.GT.2) WRITE (LUPRI,1423) DTIM
 1423 FORMAT (/' *** SETUPW-INFO, time in SETUP: ',F10.2,
     X   ' seconds.'/)
C
C
C ********************************************************
C   END OF SETUPW.
C ********************************************************
C
C
30000 CONTINUE
      RETURN
C
 9000 WRITE (6,'(//A)') ' *** ERROR *** unable to open input file'
      STOP '*** ERROR (WESINP) unable to open input file (unit 5)'
      END
#  /* Deck readiw */
#if defined (VAR_READIW)
      SUBROUTINE READIW(LUMLC,RDDELQ)
C
C     NO SUBROUTINE CALLS
C
#include "implicit.h"
      INTEGER LUMLC(2)
      LOGICAL RDDELQ
C
      PARAMETER ( PI = 3.14159 26535 89793 D0 )
      PARAMETER ( MXNUC = 40 )
C
#include "wesinf.h"
#include "wesgto.h"
C
C DO NOT CALL INFIND (OR CHANGE "IS" VARIABLE)
C
#if defined (VAR_STAR2)
      INTEGER*2 ICONT
#endif
      LOGICAL IFIF, BINPUT
#if defined (SYS_VAX)
      CHARACTER*4 SPD(5),KWO(100),NAKO,IPRX(400),JPRX(400)
      CHARACTER*1 KAKO(8,3),KO(3)
      CHARACTER*2 KOKO(3)
#endif
#if !defined (SYS_VAX)
      CHARACTER*8 SPD(5),KWO(100),NAKO,IPRX(400),JPRX(400)
      CHARACTER*8 KAKO(8,3),KO(3)
      CHARACTER*8 KOKO(3)
#endif
      DIMENSION A(50),C(50),CORR(100,3),PT(8),TITLE(24),
     1   CCONT(400,8),EXPO(999),KA(8,3),
     +   CC(50),RAAP(3),RAA(3),XNO(8,151,2),
     2   CCO(260),CEX(999),CORD(100,3),CORDP(100,3)
      DIMENSION JCO(4),NUC(4,MXNUC),MULNUC(100),
     1   NY(99),NXYZ(20,3),NAMEA(20),NAMEB(20),NAMN(100),ICONT(400,8),
     2   JCONT(400),LAMN(400),MAMN(400),IAMN(400)
C
C
      DATA SPD,KO,KOKO,NXYZ,(KWO(I),I=1,20)
     1    /'S','P','D','F','G','X','Y','Z','YZ','XZ','XY',
     2     0,1,0,0,2,1,1,0,0,0,3,2,2,1,1,1,0,0,0,0,
     3     0,0,1,0,0,1,0,2,1,0,0,1,0,2,1,0,3,2,1,0,
     4     0,0,0,1,0,0,1,0,1,2,0,0,1,0,1,2,0,1,2,3,
     5     'S  ','X  ','Y  ','Z  ','XX ','XY ','XZ ','YY ','YZ ','ZZ ',
     6     'XXX','XXY','XXZ','XYY','XYZ','XZZ','YYY','YYZ','YZZ','ZZZ'/
C
      DATA ZERO,ONE,TWO,FOUR,EIGHT,ON2,ON4
     1     /0.D0,1.D0,2.D0,4.D0,8.D0,.5D0,.25D0/
      PIPPI=ONE/(EIGHT*PI*PI*PI)**ON4
C
      LUT    = 6
      LUMLCA = LUMLC(1)
      LUMLCB = LUMLC(2)
      IF (NGEOM .EQ. 1) THEN
         BINPUT = .FALSE.
         RDDELQ = .FALSE.
      ELSE IF (LUMLCB .GT. 0 .AND. LUMLCB .NE. LUMLCA) THEN
         BINPUT = .TRUE.
      ELSE
         BINPUT = .FALSE.
      END IF
C
      READ(LUMLCA,10)
      READ(LUMLCA,10) TITLE
   10 FORMAT (12A6)
      WRITE(LUT,20) TITLE
   20 FORMAT(////,5X,'MOLECULE INPUT READ BY WESTA',
     *         //,5X,'TITLE GEOMETRY A:',
     *          /,5X,'-----------------',
     *          /,2(/,1X,12A6))
      IF (BINPUT) THEN
         READ(LUMLCB,10)
         READ(LUMLCB,10) TITLE
         WRITE(LUT,21) TITLE
   21    FORMAT(//,5X,'TITLE GEOMETRY B:',
     *           /,5X,'-----------------',
     *           /,2(/,1X,12A6))
      END IF
C
      READ(LUMLCA,30) NONTYP,NSYMOP,((KAKO(J,I),I=1,3),J=1,3),ID3,
     1                TLA,TLC,IRSTRT,IPUNCH,NJAN
      IF (BINPUT) READ (LUMLCB,'()')
   30 FORMAT (2I5,10A1,2E10.2,3I5)
      READ(LUMLCA,10)
      IF (BINPUT) READ (LUMLCB,'()')
      IF (IPRMLC .GE. 1) WRITE(LUT,80) NONTYP,NSYMOP
   80 FORMAT (///,I11,' DIFFERENT TYPES OF ATOMS',
     1          /,I11,' SYMMETRY OPERATIONS')
C
      MAXLOP    = 2**NSYMOP
      IFXYZ(1)  = 0
      IFXYZ(2)  = 0
      IFXYZ(3)  = 0
      ISYTYP(1) = 1
      ISYTYP(2) = 1
      ISYTYP(3) = 1
      IF (NSYMOP.EQ.0) GO TO 200
      DO 190 J = 1,NSYMOP
        NBL = 1
        DO 90 I = 1,3
          IF (KAKO(J,I).EQ.' ') NBL = NBL + 1
          IF (KAKO(J,I).EQ.KO(1)) ISYTYP(1)=ISYTYP(1)+2**(J-1)
          IF (KAKO(J,I).EQ.KO(2)) ISYTYP(2)=ISYTYP(2)+2**(J-1)
          IF (KAKO(J,I).EQ.KO(3)) ISYTYP(3)=ISYTYP(3)+2**(J-1)
   90   CONTINUE
        GO TO (100,120,150), NBL
  100   CONTINUE
          IF (IPRMLC .GE. 1) WRITE(LUT,'(10X,A)') 'INVERSION CENTER'
          IFXYZ(1)=1
          IFXYZ(2)=1
          IFXYZ(3)=1
        GO TO 190
  120   CONTINUE
          DO 140 K=1,3
            IF (KAKO(J,1).EQ.KO(K) .OR. KAKO(J,2).EQ.KO(K) .OR.
     1          KAKO(J,3).EQ.KO(K)) GO TO 140
            IF (IPRMLC .GE. 1) WRITE(LUT,'(T11,A,A1,A)')
     1         ' ROTATION ABOUT THE ',KO(K),'-AXIS'
            IF (K.NE.1) IFXYZ(1)=1
            IF (K.NE.2) IFXYZ(2)=1
            IF (K.NE.3) IFXYZ(3)=1
  140     CONTINUE
        GO TO 190
  150   CONTINUE
          DO 160 K=1,3
            IF (KAKO(J,1).EQ.KO(K) .OR. KAKO(J,2).EQ.KO(K) .OR.
     1          KAKO(J,3).EQ.KO(K)) GO TO 170
  160     CONTINUE
        GO TO 190
  170     IF (IPRMLC .GE. 1) WRITE(LUT,'(T11,A,A2,A)')
     1       ' REFLECTION IN THE ',KOKO(K),'-PLANE'
          IFXYZ(K)=1
  190 CONTINUE
C
  200 DO 220 I=1,8
        IS=0
        DO 210 K=1,3
          KA(I,K)=(I-1)/(2**(K-1))-2*((I-1)/(2**K))
  210     IS=IS+KA(I,K)
        PT(I)=(-ONE)**IS
        MULT(I)=2**MAX(0,NSYMOP-IS)
        FMULT(I)=MULT(I)
  220 CONTINUE
      DO 250 I=1,8
        DO 240 J=1,8
          IAN=1
          IOR=1
          DO 230 K=1,3
            IAN=IAN+MIN(KA(I,K),KA(J,K))*2**(K-1)
  230       IOR=IOR+MAX(KA(I,K),KA(J,K))*2**(K-1)
          AND88(I,J)=IAN
          OR88 (I,J)=IOR
          EOR88(I,J)=IOR-IAN+1
          PTAND(I,J)=PT(IAN)
  240   CONTINUE
  250 CONTINUE
      ITYPEX(1,1)=1
      NFO=0
      NDO=0
      DO 280 JQM=1,3
         ITYPEX(2,JQM)=ISYTYP(JQM)
         DO 270 KQM=JQM,3
            NDO=NDO+1
            IVARA=ISYTYP(JQM)
            IVARB=ISYTYP(KQM)
            IVARC=EOR88(IVARA,IVARB)
            ITYPEX(3,NDO)=IVARC
            DO 260 LQM=KQM,3
               NFO=NFO+1
               IVARA=ISYTYP(LQM)
               ITYPEX(4,NFO)=EOR88(IVARA,IVARC)
  260       CONTINUE
  270    CONTINUE
  280 CONTINUE
      KHKT(1) =  1
      KHKT(2) =  3
      KHKT(3) =  6
      KHKT(4) = 10
      KMAX = 0
      LMAX = 0
      NMAX = 0
      MAXLO= 8
      NHTYP= 0
C
      DO 500 I = 1,NONTYP
         READ(LUMLCA,290) KASET,Q,NONT,IQM,(JCO(J),J=1,IQM)
         IF (BINPUT) READ (LUMLCB, '()')
#if !defined (SYS_CRAY)
  290    FORMAT (I6,F4.1,6I5)
#endif
#if defined (SYS_CRAY)
  290    FORMAT (I6,F4.1,10X,6I5)
#endif
         IF (IPRMLC .GE. 3) THEN
            IF (IQM.GT.0) THEN
               WRITE(LUT,300) I,Q,NONT,SPD(IQM)
               WRITE(LUT,310) (JCO(J),SPD(J),J=1,IQM)
            ELSE
               WRITE(LUT,301) I,Q,NONT
            END IF
         END IF
  300    FORMAT (//'     ATOMIC TYPE NUMBER',I2,
     1           /,9X,'NUCLEAR CHARGE=',F5.2,
     2           /,I10,' SYMMETRY INDEPENDENT ATOMS',
     3           /,9X,'HIGHEST TYPE IS ',A1,/)
  301    FORMAT (//'     ATOMIC TYPE NUMBER',I2,
     1           /,9X,'NUCLEAR CHARGE=',F5.2,
     2           /,I10,' SYMMETRY INDEPENDENT ATOMS',
     3           /,9X,'POINT CHARGE, NO BASIS FUNCTIONS')
  310    FORMAT (I18,' CGTO''S OF ',A1,' TYPE')
C
         READ(LUMLCA,320) (NAMEA(N),(CORR(N,J),J=1,3),N=1,NONT)
  320    FORMAT (A4,3F10.5)
C.........................................................
C        READ DELTA Q FOR CHANGING INTERNUCLEAR DISTANCE
C
         IF (BINPUT) THEN
            READ(LUMLCB,320) (NAMEB(N),(DELQ(J,NMAX+N),J=1,3),N=1,NONT)
            DO 325 N = 1,NONT
               DELQ(1,NMAX+N) = DELQ(1,NMAX+N) - CORR(N,1)
               DELQ(2,NMAX+N) = DELQ(2,NMAX+N) - CORR(N,2)
               DELQ(3,NMAX+N) = DELQ(3,NMAX+N) - CORR(N,3)
  325       CONTINUE
         ELSE IF (RDDELQ) THEN
            READ(LUMLCA,320) (NAMEB(N),(DELQ(J,NMAX+N),J=1,3),N=1,NONT)
         END IF
C
C......................................................
         NMIN = NMAX + 1
         DO 330 N=1,NONT
            NMAX = NMAX + 1
            NAMN(NMAX)  =NAMEA(N)
            CORD(NMAX,1)=CORR(N,1)
            CORD(NMAX,2)=CORR(N,2)
            CORD(NMAX,3)=CORR(N,3)
C.........................................................
            IF (NGEOM .GT. 1) THEN
               CORDP(NMAX,1)=CORD(NMAX,1)+DELQ(1,NMAX)
               CORDP(NMAX,2)=CORD(NMAX,2)+DELQ(2,NMAX)
               CORDP(NMAX,3)=CORD(NMAX,3)+DELQ(3,NMAX)
            END IF
C.........................................................
  330    CONTINUE
C........................................................
         IF (IPRMLC .GE. 2) THEN
            WRITE(LUT,340) ((CORD(N,J),J=1,3),N=NMIN,NMAX)
            IF (NGEOM .GT. 1)
     1      WRITE(LUT,341) ((DELQ(J,N),J=1,3),N=NMIN,NMAX)
         END IF
  340    FORMAT (//,9X,'NUCLEAR COORDINATES',/,(/8X,3F18.6))
  341    FORMAT (/,9X,'DISPLACEMENT OF NUCLEAR COORDINATES',
     *           /,(/8X,3F18.6))
C........................................................
         IMAX = 0
         DO 400 J=1,IQM
            KM=JCO(J)
            IF (KM .GT. MXNUC) THEN
               WRITE (LUT,'(/A,I5)') ' TOO MANY CGTO:S, MAX IS',MXNUC
               STOP 'WESTA.MOLECULE.READIN, TOO MANY CGTO:S'
            END IF
            READ(LUMLCA,40) (NUC(J,K),K=1,KM)
            IF (BINPUT) READ (LUMLCB, '()')
   40       FORMAT(10I5)
            DO 390 K=1,KM
               IMIN=IMAX+1
               IMAX=IMAX+NUC(J,K)
               READ(LUMLCA,*) (A(L),C(L),L=IMIN,IMAX)
               IF (BINPUT) READ (LUMLCB,*) (DUM,DUM,L=IMIN,IMAX)
C              READ(LUMLCA,350) (A(L),C(L),L=IMIN,IMAX)
C              IF (BINPUT) READ (LUMLCB,350) (DUM,DUM,L=IMIN,IMAX)
C 350          FORMAT (8F9.4)
C              normalization
               SUM = ZERO
               DO 370 L=IMIN,IMAX
                  DO 360 M=IMIN,L
                     T = C(L)*C(M)*
     *                   (TWO*SQRT(A(L)*A(M))/(A(L)+A(M)))**(J+ON2)
                     IF (L.NE.M) SUM=SUM+T
                     SUM=SUM+T
  360             CONTINUE
  370          CONTINUE
               SUM=ONE/SQRT(SUM)
               DO 380 L=IMIN,IMAX
                  CC(L)=C(L)
                  C(L) =C(L)*SUM*(FOUR*A(L))**(ON2*J+ON4)*PIPPI
  380          CONTINUE
C
  390       CONTINUE
  400    CONTINUE
         DO 490 N = NMIN,NMAX
            DO 420 L = 1,IMAX
               LMAX = LMAX + 1
               ALPHA(LMAX)=A(L)
               CCO(LMAX)  =CC(L)
               CONT(LMAX) =C(L)
  420       CONTINUE
            MULK=1
            LL  =1
            DO 460 L=1,NSYMOP
               DO 440 M=1,3
                  IVARA=ISYTYP(M)
                  IF (AND88(LL+1,IVARA).EQ.1) GO TO 440
                  IF (ABS(CORD(N,M)).GE.1.D-06) GO TO 450
C....................................................
C                 IF (ABS(CORDP(N,M)).GE.1.D-06) GO TO 450
C....................................................
  440          CONTINUE
               MULK=MULK+LL
  450          LL  =2*LL
  460       CONTINUE
            MAXLO    =AND88(MULK,MAXLO)
            MULNUC(N)=MULK
            NHTYP    =MAX(NHTYP,IQM)
            DO 480 J=1,IQM
               JCOJ=JCO(J)
               DO 470 K=1,JCOJ
                  KMAX=KMAX+1
                  NUCO(KMAX)  =NUC(J,K)
                  LAMN(KMAX)  =NAMN(N)
                  MULBAS(KMAX)=MULK
                  NHKT(KMAX)  =J
                  NY(KMAX)    =0
#if defined (SYS_VAX)
                  IF (NAMN(N).EQ.' NEW'.OR.NAMN(N).EQ.'NEW ') NY(KMAX)=1
#endif
                  CENT(KMAX,1,1)=CORD(N,1)
                  CENT(KMAX,2,1)=CORD(N,2)
                  CENT(KMAX,3,1)=CORD(N,3)
C..................................................
                  IF (NGEOM .GT. 1) THEN
                     CENT(KMAX,1,2)=CORDP(N,1)
                     CENT(KMAX,2,2)=CORDP(N,2)
                     CENT(KMAX,3,2)=CORDP(N,3)
                  END IF
C..................................................
  470          CONTINUE
  480       CONTINUE
  490    CONTINUE
  500 CONTINUE
C
      IF (NHTYP .GT. 4) THEN
         WRITE (LUPRI,'(/A/A,I5)')
     &   'Sorry, WESTA only implemented for spdf orbitals',
     &   'Highest ang.mom. in basis is: ',NHTYP-1
        CALL QUIT(
     &   'Sorry, WESTA only implemented for spdf orbitals')
      END IF
C
C
C     TRANS-MOMENT INTEGRATION POINT
C
C     CCCX=0.D0
C     CCCY=0.D0
C     CCCZ=0.D0
C     WRITE(LUT,3099)CCCX,CCCY,CCCZ
C3099 FORMAT(/,1X,'INTEGRATION-POINT',3F9.6)
C
      NORB =1
      IABAS=0
      IBBAS=0
C
      DO 640 LA=1,MAXLOP
        KABAS=0
        NUMR =0
        NBASSS(LA)=0
        JSTRT(1)  =0
        NSTRT(1)  =0
C
        DO 630 IA=1,KMAX
          MULA =MULBAS(IA)
          NUCA =NUCO(IA)
          JSTA =JSTRT(IA)
          NHKTA=NHKT(IA)
          KHKTA=KHKT(NHKTA)
          JSTRT(IA+1)=JSTA+NUCA
          NSTRT(IA+1)=NSTRT(IA)+KHKTA
C
          DO 620 NA=1,KHKTA
            MXYZ  =((NHKTA+1)*NHKTA*(NHKTA-1))/6+NA
            IVARA =ITYPEX(NHKTA,NA)
            IVARB =EOR88(LA,IVARA)
            IFIF  =AND88(MULA,IVARB).EQ.1
            NUMR  =NUMR+1
            OVLAP1=ZERO
            OVLAP2=ZERO
            JKB   =0
C
            DO 600 KB=1,MAXLOP
              IF (AND88(KB,MULA).NE.1) GO TO 600
              DISTA1=ZERO
              DISTA2=ZERO
              DO 510 I=1,3
                IVA   =ISYTYP(I)
                RAA(I)=CENT(IA,I,1)*(ONE-PTAND(IVA,KB))
                DISTA1=DISTA1+RAA(I)*RAA(I)
C......................................................
                IF (NGEOM .GT. 1) THEN
                  RAAP(I)=CENT(IA,I,2)*(ONE-PTAND(IVA,KB))
                  DISTA2 =DISTA2+RAAP(I)*RAAP(I)
                END IF
C......................................................
  510         CONTINUE
              IF (LA .EQ. 1) THEN
                 IABAS = IABAS + 1
                 JCONT(IABAS) = NUCA
                 IAMN(IABAS)  = LAMN(IA)
                 IPRX(IABAS)  = KWO(MXYZ)
              END IF
              JKB  =JKB  +1
              KABAS=KABAS+1
              CTRAN(NORB,JKB)=PTAND(KB,IVARB)
              ITRAN(NORB,JKB)=KABAS
C
              DO 590 JA=1,NUCA
                APA=ALPHA(JSTA+JA)
                CAA=CCO(JSTA+JA)
                CA =CONT(JSTA+JA)
                IF (LA .EQ. 1) THEN
                   IBBAS = IBBAS + 1
                   CCONT(IABAS,JA) = CA
                   ICONT(IABAS,JA) = IBBAS
                   CEX(IBBAS) =CAA
                   EXPO(IBBAS)=APA
                END IF
                IF (.NOT.IFIF) GO TO 590
C
                DO 580 JB=1,NUCA
                  B    =ALPHA(JSTA+JB)
                  CPC  =APA+B
                  SAAB1=CA*CONT(JSTA+JB)*EXP(-APA*B*DISTA1/CPC)*
     1                  PI*SQRT(PI/CPC)/CPC
C.......................................
                  IF (NGEOM.GT.1)
     1            SAAB2=CA*CONT(JSTA+JB)*EXP(-APA*B*DISTA2/CPC)*
     1                  PI*SQRT(PI/CPC)/CPC
C.......................................
                  IF (NHKTA.EQ.1) GO TO 570
                  CPC2I = ON2 / CPC
                  DO 560 IXYZ=1,3
                    LOV  =NXYZ(MXYZ,IXYZ)
                    IVA  =ISYTYP(IXYZ)
C.........................................
                    PXCS1=CENT(IA,IXYZ,1)*(APA+B*PTAND(IVA,KB))/CPC
                    AAXS =APA*RAA(IXYZ)/CPC
                    BAXS =-B*AAXS/APA
                    ANS1 =GTOMOM(LOV,LOV,0,BAXS,AAXS,PXCS1,CPC2I)
                    SAAB1=SAAB1*ANS1
C.........................................
                    IF (NGEOM .GT. 1) THEN
                      PXCS2=CENT(IA,IXYZ,2)*(APA+B*PTAND(IVA,KB))/CPC
                      AAXSP=APA*RAAP(IXYZ)/CPC
                      BAXSP=-B*AAXSP/APA
                      ANS2 =GTOMOM(LOV,LOV,0,BAXSP,AAXSP,PXCS2,CPC2I)
                      SAAB2=SAAB2*ANS2
                    END IF
C.........................................
  560             CONTINUE
C
  570             OVLAP1=OVLAP1+SAAB1*PTAND(KB,IVARB)
C.........................................
                  IF (NGEOM.GT.1) OVLAP2=OVLAP2+SAAB2*PTAND(KB,IVARB)
C.........................................
  580             CONTINUE
  590           CONTINUE
  600         CONTINUE
C
            IF (.NOT.IFIF) GO TO 620
            XNO(LA,NUMR,1) = ONE / SQRT(OVLAP1*FMULT(MULA))
C......................................
            IF (NGEOM.GT.1)
     1         XNO(LA,NUMR,2) = ONE / SQRT(OVLAP2*FMULT(MULA))
C......................................
            JTRAN(NORB)=JKB
            JPRX(NORB) =KWO(MXYZ)
            MAMN(NORB) =LAMN(IA)
CHJ-850526: keep CTRAN as actual transformation used,
C           and not as the transformation below which would
C           give normalized symmetry orbitals.
C           DO 610 IK=1,JKB
C             CTRAN(NORB,IK,1)=CTRAN(NORB,IK,1)*XNO(LA,NUMR,1)
C.............................................
C             CTRAN(NORB,IK,2)=CTRAN(NORB,IK,2)*XNO(LA,NUMR,2)
C.............................................
C 610       CONTINUE
CHJ-E
            NORB = NORB + 1
            NBASSS(LA)=NBASSS(LA)+1
  620     CONTINUE
  630   CONTINUE
  640 CONTINUE
      NNUC=0
C
      IF (NGEOM .GT. 1) THEN
         WRITE(LUT,650)
         DO 760 N=1,NMAX
            WRITE(LUT,660) NAMN(N),(CORD(N,I),I=1,3),(CORDP(N,I),I=1,3)
  760    CONTINUE
      ELSE
         WRITE(LUT,651)
         DO 761 N=1,NMAX
            WRITE(LUT,660) NAMN(N),(CORD(N,I),I=1,3)
  761    CONTINUE
      END IF
  650 FORMAT ('1',///,9X,'NUCLEAR COORDINATES IN BOHRS',
     *            //,T23,'GEOMETRY A',T57,'GEOMETRY B',/)
  651 FORMAT ('1',///,9X,'NUCLEAR COORDINATES IN BOHRS',
     *            //,T23,'GEOMETRY A',/)
  660 FORMAT (/,4X,A4,4X,3F10.6,4X,3F10.6)
C
      IF (IPRMLC .GE. 5) THEN
         WRITE(LUT,770)
         DO 780 I=1,MAXLOP
  780       WRITE(LUT,790) (EOR88(I,J),J=1,MAXLOP)
      END IF
  770 FORMAT (///,25X,' GROUP MULTIPLICATION TABLE',/)
  790 FORMAT (30X,8I4)
C
      IF (IPRMLC .GE. 4) WRITE(LUT,800) IBBAS
  800 FORMAT ('1',//,15X,'GAUSSIAN BASIS INFO',
     1        /,I18,' UNCONTRACTED BASIS FUNCTIONS',
     2        /,8X,'ATOM',15X,'EXPONENT',6X,'COEFFICIENT')
      MAKO=0
      IAB =0
C
      DO 880 IA=1,KMAX
        NHKTA=NHKT(IA)
        NAK  =(NHKTA*(NHKTA+1)*(NHKTA-1))/6
        KHKTA=KHKT(NHKTA)
        MULA =MULBAS(IA)
        LAM  =LAMN(IA)
        DO 870 NA=1,KHKTA
          NUCA=NUCO(IA)
          NAKNA=NAK+NA
          NAKO=KWO(NAKNA)
          MAKO=MAX(MAKO,NAKNA)
          DO 860 KB=1,MAXLOP
            IF (AND88(KB,MULA).NE.1) GO TO 860
            IF (IPRMLC .GE. 4) WRITE(LUT,820) LAM,KB,NAKO
  820       FORMAT (/,5X,A4,I3,3X,A3)
            DO 850 JA=1,NUCA
              IAB=IAB+1
              IF (IPRMLC .GE. 4) WRITE(LUT,830) IAB,EXPO(IAB),CEX(IAB)
  830         FORMAT (I21,1P2D15.4)
C 840         FORMAT (I6,A4,I2,A3,3X,F18.8,43X,1X)
  850       CONTINUE
  860     CONTINUE
  870   CONTINUE
  880 CONTINUE
C
C
      IF (IPRMLC .GE. 4) THEN
         WRITE(LUT,960) IABAS
         DO 990 I=1,IABAS
            J=JCONT(I)
            WRITE(LUT,980)
     *         I,IAMN(I),IPRX(I),J,(ICONT(I,K),CCONT(I,K),K=1,J)
  990    CONTINUE
      END IF
  960 FORMAT ('1',10X,' CONTRACTION INFO',
     *        /,I14,' CONTRACTED BASIS FUNCTIONS',/)
  980 FORMAT (I5,3X,A6,3X,A3,I4,4(I4,F11.5),/,(24X,4(I4,F11.5)))
C
      NORB=NORB-1
      IF (IPRMLC .GE. 5) WRITE(LUT,1010) NORB
 1010 FORMAT ('1',//12X,' SYMMETRY TRANSFORMATION INFO',
     1             /10X,I5,' SYMMETRY ADAPTED BASIS FUNCTIONS')
      I     =0
      NBT   =0
      NSYM  =0
      NIRRED=0
C
      DO 1040 LA=1,MAXLOP
        NIRRED=NIRRED+1
        NBI   =NBASSS(LA)
      IF (NBI.EQ.0) GO TO 1040
        NBT   =NBT+NBI
        NSYM  =NSYM+1
        IF (IPRMLC .GE. 5) WRITE(LUT,1020) NIRRED
 1020   FORMAT (//,10X,'IRREDUCIBLE REPRESENTATION NUMBER',I3,/)
        DO 1030 ILA=1,NBI
          I=I+1
          J=JTRAN(I)
          IF (IPRMLC .GE .5) WRITE(LUT,980)
     *       I,MAMN(I),JPRX(I),J,(ITRAN(I,K),CTRAN(I,K),K=1,J)
 1030     CONTINUE
 1040   CONTINUE
C
      MAXLOT = 1
      II = 1
      DO 1060 I = 1,NSYMOP
         II = 2*II
         IF (AND88(II+1,MAXLO).EQ.1) MAXLOT = II
 1060 CONTINUE
C
C
      IF (NGEOM.GT.1 .AND. .NOT.RDDELQ .AND. .NOT.BINPUT) THEN
         IF (NMAXD .NE. NMAX) THEN
            WRITE (LUT,'(///A,I4,A/A,I4,A)')
     *        ' ERROR, ".DISPLA" input specified',NMAXD,' atoms,',
     *        '        MOLECULE input file gives',NMAX,' atoms.'
            STOP 'READIW ERROR, WRONG NUMBER OF ATOMS IN ".DISPLA"'
         END IF
      END IF
C
      DISPL = ZERO
      IF (NGEOM .GT. 1) THEN
         DO 2200 N = 1,NMAX
            DISPA = DELQ(1,N)*DELQ(1,N) + DELQ(2,N)*DELQ(2,N)
     *            + DELQ(3,N)*DELQ(3,N)
            MULK  = MULNUC(N)
            NUMNUC = 0
            DO 2100 KB = 1,MAXLOT
               IF (AND88(MULK,KB).EQ.1) NUMNUC = NUMNUC + 1
 2100       CONTINUE
            DISPL = DISPL + NUMNUC*SQRT(DISPA)
 2200    CONTINUE
         WRITE (LUT,'(//A,F12.8)') ' Displacement A to B: ',DISPL
      END IF
C
      RETURN
C END OF READIW.
      END
#endif
#  /* Deck makisw */
      SUBROUTINE MAKISW()
C
C     8-Dec-1986
C
C Purpose: construct ISW and ISX index vectors for current
C          orbital distribution in NISH(*)/NASH(*)/NSSH(*)
C
C          ISW orders from symmetry order to inactive/active/secondary
C          order, and ISX orders the opposit way.
C
#include "implicit.h"
C
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "infind.h"
C
      IIN  = 0
      IAC  = NISHT
      ISEC = NOCCT
      I    = 0
      DO 1030 ISYM=1,NSYM
         NISHI=NISH(ISYM)
         IF (NISHI.GT.0) THEN
            DO 1040 KISH=1,NISHI
            IIN = IIN + 1
            I=I+1
            ISW(I) = IIN
            ISX(IIN) = I
 1040       CONTINUE
         ENDIF
         NASHI=NASH(ISYM)
         IF (NASHI.GT.0) THEN
            DO 1050 KASH=1,NASHI
            IAC = IAC + 1
            I=I+1
            ISW(I)=IAC
            ISX(IAC)=I
 1050       CONTINUE
         ENDIF
         NSECI=NSSH(ISYM)
         IF (NSECI.GT.0) THEN
            DO 1060 KSEC=1,NSECI
            ISEC=ISEC+1
            I=I+1
            ISW(I)=ISEC
            ISX(ISEC)=I
 1060       CONTINUE
         ENDIF
 1030 CONTINUE
C
      RETURN
      END
#  /* Deck mkiswr */
      SUBROUTINE MKISWR()
C
C     17-Sep-1989,ha
C
C based on MAKISW
C
C Purpose: construct ISWR and ISXR index vectors for current
C          orbital distribution in NISH(*)/NASH(*)/NSSH(*),
C          whwre RAS division of active space is accounted for
C
C          ISWR orders from symmetry order to inactive/RAS1,RAS2,RAS3/secondary
C          order, and ISXR orders the opposit way.
C
C
#include "implicit.h"
C
#include "maxash.h"
#include "maxorb.h"
#include "wcborb.h"
#include "inforb.h"
#include "infind.h"
C
C dimensions:
      NAS1T = 0
      NAS2T = 0
      NAS3T = 0
      DO 1000 ISYM = 1, NSYM
        NAS1T = NAS1T + NAS1(ISYM)
        NAS2T = NAS2T + NAS2(ISYM)
        NAS3T = NAS3T + NAS3(ISYM)
 1000 CONTINUE
C
      IIN   = 0
      IAC1  = NISHT
      IAC2  = NISHT + NAS1T
      IAC3  = IAC2  + NAS2T
      ISEC  = NOCCT
      I     = 0
      DO 1030 ISYM=1,NSYM
         NISHI=NISH(ISYM)
         IF (NISHI.GT.0) THEN
            DO 1040 KISH=1,NISHI
            IIN = IIN + 1
            I=I+1
            ISWR(I) = IIN
            ISXR(IIN) = I
 1040       CONTINUE
         ENDIF
C
         NAS1I=NAS1(ISYM)
         IF (NAS1I.GT.0) THEN
            DO 1051 KASH=1,NAS1I
            IAC1 = IAC1 + 1
            I=I+1
            ISWR(I)=IAC1
            ISXR(IAC1)=I
 1051       CONTINUE
         ENDIF
C
         NAS2I=NAS2(ISYM)
         IF (NAS2I.GT.0) THEN
            DO 1052 KASH=1,NAS2I
            IAC2 = IAC2 + 1
            I=I+1
            ISWR(I)=IAC2
            ISXR(IAC2)=I
 1052       CONTINUE
         ENDIF
C
         NAS3I=NAS3(ISYM)
         IF (NAS3I.GT.0) THEN
            DO 1053 KASH=1,NAS3I
            IAC3 = IAC3 + 1
            I=I+1
            ISWR(I)=IAC3
            ISXR(IAC3)=I
 1053       CONTINUE
         ENDIF
C
         NSECI=NSSH(ISYM)
         IF (NSECI.GT.0) THEN
            DO 1060 KSEC=1,NSECI
            ISEC=ISEC+1
            I=I+1
            ISWR(I)=ISEC
            ISXR(ISEC)=I
 1060       CONTINUE
         ENDIF
 1030 CONTINUE
C
      RETURN
      END
#  /* Deck wbdinf */
      BLOCK DATA WBDINF
C
C  18-May-1995
C
C  Initialize common block WESINF
C
#include "implicit.h"
#include "wesinf.h"
      DATA IPRDEF/0/,  IPRNAD/-1/, IPRTRM/-1/, IPRBIO/-1/,
     *     IPRMLC/-1/, IPRSHK/-1/, IACFLG/-1/,
     *     NONADI/.FALSE./, TRAMOM/.FALSE./, NSTAT2/1/,
     *     MOMTYP/10*.FALSE./, NGEOM/1/, IRDMOW/10/,
     *     SYMORT/.FALSE./, WAVEOV/.FALSE./,
     *     CENTOP/3*0.0D0/,LSYM1/1/,LSYM2/0/,LSPIN1/1/,LSPIN2/0/,
     *     SHAKE/.FALSE./,REDUCE/.FALSE./,FROZEN/.FALSE./,
     *     OCEINT/.FALSE./,CICALC/.FALSE./,
     *     NUCMASS/2*1.007276/
C  2*proton mass, value not double-checked /hjaaj sep 04 TODO
      END
C  /* Deck sbdorb */
      BLOCK DATA SBDORB
C     Initialize MULD2H in common block /INFORB/
#include <implicit.h>
#include <inforb.h>
C
C     MULTIPLICATION TABLE FOR SYMMETRIES (MULD2H is in /INFORB/)
C
      DATA MULD2H/1,2,3,4,5,6,7,8,
     *            2,1,4,3,6,5,8,7,
     *            3,4,1,2,7,8,5,6,
     *            4,3,2,1,8,7,6,5,
     *            5,6,7,8,1,2,3,4,
     *            6,5,8,7,2,1,4,3,
     *            7,8,5,6,3,4,1,2,
     *            8,7,6,5,4,3,2,1/
      END

