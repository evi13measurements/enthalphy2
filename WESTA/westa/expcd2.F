C
C   FILE: expcd2.F
C   Purpose: CI files for WESSHK
C
      SUBROUTINE EXPCD2(CII,CIE,ISMI,ISME,LUSTRI,XNDXCI,
     &                  IWRK,LWORK,IAC,IACORB)
C
C a CI-vector CII is given in some space, I.
C Add(IAC=1) or remove(IAC=2) orbital IACORB from this vector
C and expand into extended CI space E.
C
C The number of inactive electrons is assumed constant
C
C This is an updated version of EXPCID , to be used
C in connection with electron addition/electrin annihilation .
C
C
C Driver routine
C
#include "implicit.h"
#include "detbas.h"
#include "ciinfo.h"
#include "mxpdim.h"
#include "strnum.h"
#include "inforb.h"

      DIMENSION CII(*),CIE(*)
      DIMENSION XNDXCI(*),IWRK(LWORK)
C
C
C Jeppe Olsen November 1988
C
C SPACE FOR REORDERING OF ORBITALS
C
      LWORKG = LWORK
      KLFREE = 1
      KLREOB = 1
      KLFREE = KLFREE + NASHT
      LWORK = LWORK - NASHT
      CALL ISTVC2(IWRK(KLREOB),0,1,NASHT)
C. Set offsets to symmetry-type blocks
      CALL CIOFF(ISME,1,XNDXCI,NTEST)
C
      NTEST = 0
      IF( IAC .EQ. 0 ) THEN
        CALL EXPCI(XNDXCI(KIASTR),XNDXCI(KIBSTR),XNDXCI(KSTASA),
     &           XNDXCI(KSTASB),XNDXCI(KSTBAA),XNDXCI(KSTBAB),
     &           XNDXCI(KNSSOA),XNDXCI(KNSSOB),XNDXCI(KCOFF),
     &           XNDXCI(KICOOS),XNDXCI(KIOCOC),NASTR,NBSTR,
     &           NAEL,NBEL,MAXSYM,MULD2H,NOCTPA,NOCTPB,
     &           LUSTRI,IWRK(KLFREE),
     &           EQUAL,IWRK(KLREOB),ISMI,ISME,CII,CIE,
     &           XNDXCI(KTPFSA),XNDXCI(KTPFSB),LWORK,
     &           NDTASM,NASHT,XNDXCI(KISSOA),XNDXCI(KISSOB),NTEST)
      ELSE
        CALL EXPCI2(XNDXCI(KIASTR),XNDXCI(KIBSTR),XNDXCI(KSTASA),
     &           XNDXCI(KSTASB),XNDXCI(KSTBAA),XNDXCI(KSTBAB),
     &           XNDXCI(KNSSOA),XNDXCI(KNSSOB),XNDXCI(KISSOA),
     &           XNDXCI(KISSOB),XNDXCI(KCOFF),
     &           XNDXCI(KICOOS),XNDXCI(KIOCOC),NASTR,NBSTR,
     &           NAEL,NBEL,MAXSYM,MULD2H,
     &           NOCTPA,NOCTPB,LUSTRI,IWRK(KLFREE),
     &           EQUAL,IWRK(KLREOB),ISMI,ISME,CII,CIE,
     &           XNDXCI(KTPFSA),XNDXCI(KTPFSB),LWORK,
     &           NDTASM,NASHT,IAC,IACORB)
      END IF
      RETURN
      END
#  /* Deck expci2 */
      SUBROUTINE EXPCI2(IASTR,IBSTR,NSTASA,NSTASB,ISTBAA,ISTBAB,NSSOA,
     &           NSSOB,ISSOA,ISSOB,
     &           ICS,ICOOS,IOCOC,NASTR,NBSTR,NAEL,NBEL,
     &           MAXSME,SYMPRO,NOCTPA,NOCTPB,LUSTRI,IWORK,
     &           EQUAL,IREOB,IRFSMI,IRFSME,CII,CIE,
     &           ITPFSA,ITPFSB,LWORK,NDTASM,NORB,IAC,IACORB)
C
C Extended version of EXPCI
C
C IAC = 1 : annihilate orbital IACORB before expanding to Extended graph
C IAC = 2 : add        orbital IACORB before expanding to extended graph
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EQUAL,EQUALI
      INTEGER CALTPE,CALTPI
      DIMENSION IASTR(NAEL,NASTR),IBSTR(NBEL,NBSTR)
      DIMENSION NSTASA(MAXSME),NSTASB(MAXSME)
      DIMENSION ISTBAA(MAXSME),ISTBAB(MAXSME)
      DIMENSION NSSOA(*),NSSOB(*), ISSOA(*),ISSOB(*)
      DIMENSION ICS(MAXSME),ICOOS(NOCTPB,NOCTPA,MAXSME)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ITPFSA(*),ITPFSB(*)
      INTEGER SYMPRO(8,8),NDTASM(MAXSME)
      DIMENSION IWORK(*)
C
      NTEST = 0
C
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' **************** '
        WRITE(6,*) ' Entering EXPCI2 '
        WRITE(6,*) ' **************** '
      END IF
      NDET = NDTASM(IRFSME)
C DIMENSION OF IWORK SHOULD AT LEAST BE ???
C
C
C 1 : Read in information about dimensions of previous ras space I
C     ( Done at this place for dynamic allocation )
C
      REWIND LUSTRI
      READ(LUSTRI)
      READ(LUSTRI) MMDET
C MAXIMAL NUMBER OF SYMMETRIES
      READ(LUSTRI) MXSMI,JTYP
C. New
      READ(LUSTRI) MAEL,MBEL
C NUMBER OF OCCUPATION TYPES FOR ALPHA STRINGS AND BETA strings
      READ(LUSTRI) MOCTPA,MOCTPB
C TOTAL NUMBER OF ALPHA STRINGS AND BETA STRINGS
      READ(LUSTRI) MASTR,MBSTR
C
      WRITE(6,*) ' Number of dets from LUSTRI   ',MMDET
      WRITE(6,*) ' MXSMI AND JTYP FROM LUSTRI   ',MXSMI,JTYP
      WRITE(6,*) ' MOCTPA AND MOCTPB FROM LUSTRI',MOCTPA,MOCTPB
      WRITE(6,*) ' MAEL AND MBEL FROM LUSTRI',MAEL, MBEL
      WRITE(6,*) ' MASTR AND MBSTR FROM LUSTRI',MASTR, MBSTR
C
      IF(MAEL.EQ.MBEL) THEN
        EQUALI = .TRUE.
      ELSE
        EQUALI = .FALSE.
      END IF
C
      MEL = MAEL + MBEL
      NEL = NAEL + NBEL
      IF(.NOT. ((IAC .EQ. 0 .AND. NEL .EQ. MEL  ) .OR.
     &          (IAC .EQ. 1 .AND. NEL .EQ. MEL-1) .OR.
     &          (IAC .EQ. 2 .AND. NEL .EQ. MEL+1))) THEN
        WRITE(6,*) ' Inconsistent specs for EXPCI2 !!'
        WRITE(6,*) ' MEL NEL IAC ', MEL,NEL,IAC
        STOP ' Inconsistent specs for  EXPCI2 '
      END IF
C
C** 2 : MEMORY ALLOCATION ALLOWING INFORMATION ABOUT OLD
C       STRINGS TO BE READ IN , AND REGENERATION OF
C       NUMBERING INFORMATION ABOUT CURRENT SET OF STRINGS
C
C* 2.1 : FOR CURRENT STRINGS
      KLFREE = 1
      KLZA   = KLFREE
      KLFREE = KLFREE + NORB*NAEL
      IF( EQUAL ) THEN
        KLZB = KLZA
      ELSE
        KLZB = KLFREE
        KLFREE = KLFREE + NORB * NBEL
      END IF
C ORDERING OF E STRINGS
COLD  MAXSTA = IBION(NORB,NAEL)
COLD  MAXSTB = IBION(NORB,NBEL)
COLD  MAXSTR = MAX(MAXSTA,MAXSTB)
      MAXNSTR = MAX(NASTR,NBSTR)
      KLORDA = KLFREE
      KLFREE = KLORDA + NASTR
      IF( EQUAL ) THEN
        KLORDB = KLORDA
       ELSE
        KLORDB = KLFREE
        KLFREE = KLORDB + NBSTR
       END IF
       KLORDS = KLFREE
       KLFREE = KLORDS + MAXNSTR
C SYMMETRY OF STRINGS
       KISMSA = KLFREE
       KLFREE = KLFREE + NASTR
       IF( EQUAL ) THEN
         KISMSB = KISMSA
       ELSE
         KISMSB = KLFREE
         KLFREE = KLFREE + NBSTR
       END IF
C
C I NUMBERING TO E NUMBERING
C. Must have dimension of input strings so :
       KIEA   = KLFREE
       KLFREE = KLFREE + MASTR
       KIESA  = KLFREE
       KLFREE = KLFREE + MASTR
C
       KIEB   = KLFREE
       KLFREE = KLFREE + MBSTR
       KIESB  = KLFREE
       KLFREE = KLFREE + MBSTR
C
C. 2.2 :  ARRAYS FOR INFORMATION ABOUT I STRINGS
C
       KJASTR = KLFREE
       KLFREE = KJASTR + MAEL * MASTR
       IF ( EQUALI ) THEN
         KJBSTR = KJASTR
       ELSE
         KJBSTR = KLFREE
         KLFREE = KJBSTR + MBEL * MBSTR
       END IF
C
       KMSASA = KLFREE
       KLFREE = KMSASA + MXSMI
       IF ( EQUALI) THEN
         KMSASB = KMSASA
       ELSE
         KMSASB = KLFREE
         KLFREE = KMSASB + MXSMI
       END IF
C
       KJSBAA = KLFREE
       KLFREE = KJSBAA + MXSMI
       IF ( EQUALI ) THEN
        KJSBAB = KJSBAA
       ELSE
        KJSBAB = KLFREE
        KLFREE = KJSBAB + MXSMI
       END IF
C
       KJCS   = KLFREE
       KLFREE = KJCS + MXSMI
C
       KJCOOS = KLFREE
       KLFREE = KJCOOS + MXSMI*MOCTPA*MOCTPB
C
       KJOCOC = KLFREE
       KLFREE = KJOCOC + MOCTPA*MOCTPB
C
       IF( JTYP .EQ. 2 ) THEN
         KMSOSA = KLFREE
         KLFREE = KMSOSA + MOCTPA * MXSMI
         KJSOSA = KLFREE
         KLFREE = KLFREE + MOCTPA * MXSMI
         IF ( EQUALI) THEN
           KMSOSB = KMSOSA
           KJSOSB = KJSOSA
         ELSE
           KMSOSB = KLFREE
           KLFREE = KMSOSB + MOCTPB * MXSMI
           KJSOSB = KLFREE
           KLFREE = KLFREE + MOCTPB * MXSMI
         END IF
       END IF
C
C .. 3 : Regenerate numbering for expanded set of strings
C
C.. For alpha strings ..
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Alpha expanded strings '
         WRITE(6,*) ' ======================='
      END IF
C
      CALL WEIGHT(IWORK(KLZA),1,IWORK(KLFREE),NTEST)
C. Actual address to reverse lexical order map
      CALL ISTVC2(IWORK(KLORDA),0,1,NASTR)
      DO IA = 1, NASTR
        IWORK(KLORDS-1+IA) =
     &  IZNUM(IASTR(1,IA),NAEL,IWORK(KLZA),NAEL,NORB,IWORK(KLORDA),
     &        NTEST )
      END DO
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Reverse lexical number for alpha strings '
        CALL IWRTMA(IWORK(KLORDS),1,NASTR,1,NASTR)
      END IF
C. Reverse lexical ordering to actual address
      CALL ISETVC(IWORK(KLORDA),0,NASTR)
      DO IA = 1, NASTR
        IWORK(KLORDA-1+IWORK(KLORDS-1+IA) ) = IA
      END DO
      IF( NTEST .GE.100 ) THEN
        WRITE(6,*)
     &  ' Reverse lexicial order => actual order ( alpha strings)'
        CALL IWRTMA(IWORK(KLORDA),1,NASTR,1,NASTR)
      END IF
C.. For beta strings ..
      IF( .NOT. EQUAL ) THEN
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' beta expanded strings '
         WRITE(6,*) ' ======================='
      END IF
        CALL WEIGHT(IWORK(KLZB),2,IWORK(KLFREE),NTEST)
C. Actual address to reverse lexical order map
        CALL ISTVC2(IWORK(KLORDB),0,1,NBSTR)
        DO IB = 1, NBSTR
          IWORK(KLORDS-1+IB) =
     &    IZNUM(IBSTR(1,IB),NBEL,IWORK(KLZB),NBEL,NORB,IWORK(KLORDB),
     &          NTEST )
        END DO
        IF( NTEST .GE. 100 ) THEN
          WRITE(6,*) ' Reverse lexical number for beta strings '
          CALL IWRTMA(IWORK(KLORDS),1,NBSTR,1,NBSTR)
        END IF
C. Reverse lexical ordering to actual address
        CALL ISETVC(IWORK(KLORDB),0,NBSTR)
        DO IB = 1, NBSTR
          IWORK(KLORDB-1+IWORK(KLORDS-1+IB) ) = IB
        END DO
        IF( NTEST .NE. 0 ) THEN
          WRITE(6,*)
     &    ' Reverse lexcial order => actual order ( beta strings)'
          CALL IWRTMA(IWORK(KLORDB),1,NBSTR,1,NBSTR)
        END IF
      END IF
C. Symmetry array for E strings
      DO 200 ISM = 1,MAXSME
        IF(ISM.EQ.1)THEN
          IOFF = 1
        ELSE
           IOFF = IOFF + NSTASA(ISM-1)
        END IF
        CALL ISETVC(IWORK(KISMSA-1+IOFF),ISM,NSTASA(ISM) )
  200 CONTINUE
      IF(.NOT. EQUAL) THEN
        DO 201 ISM = 1,MAXSME
          IF(ISM.EQ.1)THEN
            IOFF = 1
          ELSE
            IOFF = IOFF + NSTASB(ISM-1)
          END IF
          CALL ISETVC(IWORK(KISMSB-1+IOFF),ISM,NSTASB(ISM) )
  201   CONTINUE
      END IF

C
C... 4 : Retrieve information about previous sets of strings
C
      READ(LUSTRI) ( IWORK(KMSASA-1+I),I=1,MXSMI)
      IF(.NOT.EQUALI)
     &READ(LUSTRI) ( IWORK(KMSASB-1+I),I=1,MXSMI)
      READ(LUSTRI) ( IWORK(KJSBAA-1+I),I=1,MXSMI)
      IF(.NOT.EQUALI)
     &READ(LUSTRI) ( IWORK(KJSBAB-1+I),I=1,MXSMI)
C
      READ(LUSTRI)( IWORK(KJASTR-1+I),I=1,MAEL*MASTR)
      IF(.NOT.EQUALI)
     &READ(LUSTRI)( IWORK(KJBSTR-1+I),I=1,MBEL*MBSTR)
C
      IF ( JTYP .EQ. 2 ) THEN
        READ(LUSTRI) (IWORK(KMSOSA-1+I),I=1,MXSMI*MOCTPA)
        READ(LUSTRI) (IWORK(KJSOSA-1+I),I=1,MXSMI*MOCTPA)
        IF( .NOT. EQUALI ) THEN
          READ(LUSTRI) (IWORK(KMSOSB-1+I),I=1,MXSMI*MOCTPB)
          READ(LUSTRI) (IWORK(KJSOSB-1+I),I=1,MXSMI*MOCTPB)
        END IF
      END IF
C
      READ(LUSTRI) (  IWORK(KJCS-1+I),I=1,MXSMI)
      READ(LUSTRI) (  IWORK(KJCOOS-1+I),I=1,MXSMI*MOCTPA*MOCTPB)
      READ(LUSTRI) (  IWORK(KJOCOC-1+I),I=1,MOCTPA*MOCTPB)
c     write(6,*) ' string information read in ! '
C
C .. 4 : Expanded number for string read in
C
      IF( IAC. EQ. 1 ) THEN
        IF(NAEL .EQ. MAEL - 1 ) THEN
          IACA = 1
          IACB = 0
        ELSE
          IACA = 0
          IACB = 1
        END IF
      ELSE IF ( IAC .EQ. 2 ) THEN
        IF(NAEL .EQ. MAEL + 1 ) THEN
          IACA = 2
          IACB = 0
        ELSE
          IACA = 0
          IACB = 2
        END IF
      END IF
      CALL EXPST2(IWORK(KJASTR),MASTR,NAEL,MAEL,NORB,IWORK(KLORDA),
     &            IWORK(KLZA),IWORK(KIEA),IWORK(KIESA),NASTR,
     &            IWORK(KLFREE),IREOB,IACA,IACORB)

      CALL EXPST2(IWORK(KJBSTR),MBSTR,NBEL,MBEL,NORB,IWORK(KLORDB),
     &              IWORK(KLZB),IWORK(KIEB),IWORK(KIESB),NBSTR,
     &              IWORK(KLFREE),IREOB,IACB,IACORB)
      IF(NTEST.GE.10) WRITE(6,*) ' finished with EXPST2 '
C
C .. 5 : EXPAND OLD CIVECTOR INTO NEW CI SPACE
C
      CALL EXPCV2(NSTASA,NSTASB,ISTBAA,ISTBAB,
     &            IOCOC,ICS,ICOOS,NOCTPA,NOCTPB,
     &            NSSOA,ISSOA,NSSOB,ISSOB,
     &            IWORK(KMSASA),IWORK(KMSASB),IWORK(KJSBAA),
     &            IWORK(KJSBAB),JTYP,
     &            IWORK(KJOCOC),IWORK(KJCS),
     &            IWORK(KJCOOS),MOCTPA,MOCTPB,IWORK(KMSOSA),
     &            IWORK(KJSOSA),IWORK(KMSOSB),IWORK(KJSOSB),
     &            IWORK(KIEA),IWORK(KIESA),CII,CIE,
     &            1,MXSMI,MAXSME,SYMPRO,IWORK(KIEB),IWORK(KIESB),
     &            IRFSME,IRFSMI,NDET,ITPFSA,ITPFSB,
     &            IWORK(KISMSA),IWORK(KISMSB))


C     EXPCV2(NSTASA,NSTASB,ISTBAA,ISTBAB,
C    &                 IOCOC,ICS,ICOOS,NOCTPA,NOCTPB,
C    &                 NSTAOS,ISTAOS,NSTBOS,ISTBOS,
C    &                 MSTASA,MSTASB,JSTBAA,JSTBAB,JTYP,
C    &                 JOCOC,JCS,JCOOS,MOCTPA,MOCTPB,
C    &                 MSTAOS,JSTAOS,MSTBOS,JSTBOS,
C    &                 ITOEA,ITOESA,CII,CIE,IZERO,MXSMI,MXSME,SYMPRO,
C    &                 ITOEB,ITOESB,
C    &                 IREFSM,JREFSM,NDET,ITPSTA,ITPSTB,
C    &                 ISMASA,ISMASB)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' **************** '
        WRITE(6,*) ' Leaving EXPCI2 '
        WRITE(6,*) ' **************** '
      END IF
      RETURN
      END
#  /* Deck expst2 */
      SUBROUTINE EXPST2(JSTR,MSTR,NEL,MEL,NORB,
     &                  IHTOE,Z,IITOE,IITOES,NSTR,
     &                  IWORK,IREORB,IAC,IACORB)
C
C Two sets of strings are given, an expanded basis and an
C initial basis . The initial set of strings is
C mapped to the expanded set of strings by
C    annihilating an electron ( IAC = 1)
C    creating an electron     ( IAC = 2 )
C    nothing                  ( IAC = 3 )
C
C The electron to be removed/added is IACORB
C Obtain number and phase shift  going
C from initial to expanded form
C
C FOR IAC = 1,2 some strings can vanish . this is indicated by
C           zeroes in IITOE and IITOES
C
C   INPUT :
C  ========
C
C EXPANDED STRINGS :
C     IHTOE : HANDY NUMBER TO ACTUAL STRING NUMBER
C     Z     : HANDY'S Z MATRIX FOR STRINGS
C
C INITIAL STRINGS
C    JSTR : OCCUPATION OF STRINGS
C    MSTR : NUMBER OF STRINGS
C    IREORB : REORDERING OF ORBITALS , GOING
C             FROM INITIAL TO EXPANDED BASIS
C    NSTR   : NUMBER OF EXPANDED STRINGS
C
C  OUTPUT :
C ==========
C
C    IITOE : IITOE(I) IS NUMBER IN EXPANDED GRAPH FOR
C            INPUT STRING I
C    IITOES : SIGN ARRAY , FROM REORDERING OF ORBITALS
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c. Input
      DIMENSION JSTR(MEL,MSTR)
      DIMENSION IHTOE(NSTR)
      DIMENSION IREORB(NORB)
      INTEGER*4 Z(*)
c. Output
      DIMENSION IITOE(MSTR)
      DIMENSION IITOES(MSTR)
c. Scratch
      DIMENSION IWORK(*)
c
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' ==============='
       WRITE(6,*) ' EXPST2 speaking '
       WRITE(6,*) ' ==============='
       WRITE(6,*) ' input strings '
       CALL IWRTMA(JSTR,MEL,MSTR,MEL,MSTR)
       WRITE(6,*) ' IAC and IACORB ', IAC,IACORB
      END IF
c
      KL1 = 1
      KL2 = KL1 + MAX(MEL,NEL)
c
      DO 100 JA = 1, MSTR
C?      write(6,*) ' Input string '
C?      CALL IWRTMA(JSTR(1,JA),MEL,1,MEL,1)
c Reorder input string
        DO IEL = 1, MEL
          IWORK(IEL) = IREORB(JSTR(IEL,JA) )
        END DO
c Change to ascending order
C            ORDSTR(IINST,IOUTST,NELMNT,ISIGN )
        CALL ORDSTR(IWORK,IWORK,MEL,ISIGN )
        IF(IAC .NE. 0 ) THEN
C              ASTRI(IAC,IORB,ISIGN,INSTR,UTSTR,NEL)
          CALL ASTRI(IAC,IACORB,ISIGN2,IWORK(KL1),
     &                IWORK(KL2),MEL)
          CALL ICOPVE(IWORK(KL2),IWORK(KL1),NEL)
        ELSE
          ISIGN2 = 1
        END IF
C?      write(6,*) ' Output string '
C?      CALL IWRTMA(IWORK(KL1),NEL,1,NEL,1)

c Address of this string in expanded basis
C       IZNUM(IOCC,NEL,Z,MAXEL,MAXORB,NEWORD,NTEST)
        IF(ISIGN2 .NE. 0 )
     &  IENUM= IZNUM(IWORK,NEL,Z,NEL,NORB,IHTOE,NTEST)
        IF( IENUM .GT. 0 .AND. ISIGN2 .NE. 0 ) THEN
          IITOE(JA) = IENUM
          IITOES(JA) = ISIGN*ISIGN2
        ELSE
          IITOE(JA) = 0
          IITOES(JA) = 0
        END IF
  100 CONTINUE
C
      NTEST = 0
      IF ( NTEST .GE. 10 ) THEN
        WRITE(6,*)
     &   ' EXPST2 : Place of input strings in expanded graph '
        CALL IWRTMA(IITOE,1,MSTR,1,MSTR)
        WRITE(6,*) 'EXPST2 :  Sign array '
        CALL IWRTMA(IITOES,1,MSTR,1,MSTR)
      END IF
C
      RETURN
      END
      SUBROUTINE ASTRI(IAC,IORB,ISIGN,INSTR,UTSTR,NEL)
C
C IAC = 1 : Annihilate orbital IORB from strin INSTR
C           to give string UTSTR
C
C IAC = 2 : Add orbital ( create ) IORB to string INSTR
C          to give string UTSTR
C
C Input string is assumed to be in ascending order , output
C string will also be delivered in ascending order .
C
C ISIGN = 0 indicates resulting string is a zero string
C ISIGN = +/- 1 gives sign to be multiplied with final string
C
C Jeppe Olsen , March 1988
C
      DIMENSION INSTR(*)
      INTEGER   UTSTR(*)
C23456
      IORBIN = 0
      ISIGN = 0
C
      DO I = 1, NEL
        IF( INSTR(I) .EQ. IORB) IORBIN = 1
      END DO
C. annihilate
      IF(IAC.EQ.1 .AND. IORBIN .EQ. 1 ) THEN
        IPLAC = 1
        DO IEL = 1,NEL
          IF(INSTR(IEL) .NE. IORB ) THEN
            UTSTR(IPLAC) = INSTR(IEL)
            IPLAC = IPLAC + 1
          ELSE
            ISIGN = (-1) ** (IEL - 1 )
          END IF
        END DO
C. Create
      ELSE IF ( IAC .EQ. 2 .AND. IORBIN .EQ. 0 ) THEN
        IPLAC = 0
        IADDED = 0
        DO IEL = 1,NEL
         IPLAC = IPLAC + 1
         IF(IORB.LT. INSTR(IEL).AND. IADDED.EQ.0 ) THEN
           UTSTR(IPLAC) = IORB
           ISIGN = (-1) ** ( IPLAC -1 )
           IPLAC = IPLAC + 1
           IADDED = 1
         END IF
         UTSTR(IPLAC) = INSTR(IEL)
        END DO
C
        IF(IADDED .EQ. 0 ) THEN
          UTSTR(NEL + 1 ) = IORB
          ISIGN = (-1) ** (NEL)
        END IF
C
      END IF
C
      RETURN
      END
#  /* Deck expcv2 */
      SUBROUTINE EXPCV2(NSTASA,NSTASB,ISTBAA,ISTBAB,
     &                 IOCOC,ICS,ICOOS,NOCTPA,NOCTPB,
     &                 NSTAOS,ISTAOS,NSTBOS,ISTBOS,
     &                 MSTASA,MSTASB,JSTBAA,JSTBAB,JTYP,
     &                 JOCOC,JCS,JCOOS,MOCTPA,MOCTPB,
     &                 MSTAOS,JSTAOS,MSTBOS,JSTBOS,
     &                 ITOEA,ITOESA,CII,CIE,IZERO,MXSMI,MXSME,SYMPRO,
     &                 ITOEB,ITOESB,
     &                 IREFSM,JREFSM,NDET,ITPSTA,ITPSTB,
     &                 ISMASA,ISMASB)
C
C A CI EXPANSION IS GIVEN IN A RAS SPACE I
C OBTAIN THE SAME VECTOR SPACE IN AN EXPANDED RASSPACE
C
C slightly modified version of EXPCIV
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NSTASA(MXSME),NSTASB(MXSME)
      DIMENSION ISTBAA(MXSME),ISTBAB(MXSME)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ICS(MXSME),ICOOS(NOCTPB,NOCTPA,MXSME)
      DIMENSION NSTAOS(NOCTPA,MXSME),NSTBOS(NOCTPB,MXSME)
      DIMENSION ISTAOS(NOCTPA,MXSME),ISTBOS(NOCTPB,MXSME)
      DIMENSION ISMASA(*),ISMASB(*)
      DIMENSION ITPSTA(*),ITPSTB(*)
C
      DIMENSION MSTASA(MXSMI),MSTASB(MXSMI)
      DIMENSION JSTBAA(MXSMI),JSTBAB(MXSMI)
      DIMENSION JOCOC(MOCTPA,MOCTPB)
      DIMENSION JCS(MXSME),JCOOS(MOCTPB,MOCTPA,MXSMI)
      DIMENSION MSTAOS(MOCTPA,MXSMI),MSTBOS(MOCTPB,MXSMI)
      DIMENSION JSTAOS(MOCTPA,MXSMI),JSTBOS(MOCTPB,MXSMI)
C
      INTEGER SYMPRO(8,8)
      DIMENSION ITOEA(*),ITOESA(*)
      DIMENSION ITOEB(*),ITOESB(*)
      DIMENSION CII(*),CIE(*)
C
      NTEST = 00
c
      IF(NTEST.GE.5) THEN
       WRITE(6,*) ' ****************'
       WRITE(6,*) ' Entering EXPCV2 '
       WRITE(6,*) ' ****************'
       WRITE(6,*)
       WRITE(6,*) ' JREFSM = ', JREFSM
C?     write(6,*) ' first elements '
C?     write(6,*) NSTAOS(1,1),ISTAOS(1,1),NSTBOS(1,1),ISTBOS(1,1)
      END IF
c
      IF( IZERO .NE. 0 ) THEN
       CALL SETVEC(CIE,0.0D0,NDET)
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' CIE zeroed '
C LOOP OVER CII
      JDET = 0
      DO 1000 JASM = 1, MXSMI
        JBSM = SYMPRO(JASM,JREFSM)
C?      write(6,*) ' jasm, jbsm ', jasm,jbsm
        DO 900 JATP = 1, MOCTPA
          NJA =   MSTAOS(JATP,JASM)
          JASTR = JSTAOS(JATP,JASM)
C?        write(6,*) ' jatp nja  jastr ', jatp,nja,jastr
          DO 800 JBTP = 1, MOCTPB
            NJB =   MSTBOS(JBTP,JBSM)
            JBSTR = JSTBOS(JBTP,JBSM)
C?          write(6,*) ' jbtp njb  jbstr ', jbtp,njb,jbstr
            IF( JOCOC(JATP,JBTP) .EQ. 0 ) GOTO 800
C
            DO 700 JB = JBSTR,JBSTR+NJB-1
C
              IB = ITOEB(JB)
C?            write(6,*) ' jb ib ', jb,ib
              IF( IB .NE. 0 ) THEN
                IBSM = ISMASB(IB)
                SIGN1= ITOESB(JB)
                IBTP = ITPSTB(IB)
C?              write(6,*) ' ibsm ibtp ', ibsm,ibtp
                IBSTR = ISTBOS(IBTP,IBSM)
                NIB =   NSTBOS(IBTP,IBSM)
                IBREL = IB - IBSTR + 1
C?              write(6,*) ' ibrel ', ibrel
              END IF
C
              DO 600 JA = JASTR,JASTR+NJA-1
                JDET = JDET + 1
                IA = ITOEA(JA)
C?              write(6,*) ' ja ia ', ja, ia
                IF( IA .EQ. 0.OR. IB.EQ. 0 ) GOTO 600
                IASM = ISMASA(IA)
                SIGN = SIGN1* ITOESA(JA)
                IATP = ITPSTA(IA)
                IASTR = ISTAOS(IATP,IASM)
                NIA = NSTAOS(IATP,IASM)
                IF(IOCOC(IATP,IBTP) .NE. 1 ) GOTO 600
                IAREL = IA - IASTR + 1
C AND THEN THE NEW ADRESS
                IOFF = ICOOS(IBTP,IATP,IASM)
                IDET = IOFF - 1 + (IBREL-1)*NIA+IAREL
C?              write(6,*) ' ibtp iatp iasm ',IBTP,IATP,IASM
C?              write(6,*) ' ioff ibrel nia iarel'
C?              write(6,*)  ioff,ibrel,nia,iarel

C
                CIE(IDET) = CII(JDET)*SIGN
  600         CONTINUE
  700       CONTINUE
  800     CONTINUE
  900   CONTINUE
 1000 CONTINUE
C
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' INPUT AND OUTPUT VECTORS FROM EXPCIV : '
        CALL WRTMAT(CII,1,JDET,1,JDET,0)
        CALL WRTMAT(CIE,1,NDET,1,NDET,0)
      END IF
C
      RETURN
      END
