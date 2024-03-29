C
C  /* Deck so_trmom */
      SUBROUTINE SO_TRMOM(MODEL,ISYMTR,NLBTOT,NEXCI,
     &                    T2AM,LT2AM,DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,
CSPAS:23/5-11: second and third moment sum rules
C    &                    TRLEN,TRVEL,TRLON,TRMAG,BSRLON,EXENG,
C    &                    WORK,LWORK)
     &                    TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,TRLON,TRMAG,
     &                    BSRLON,EXENG,
CClark:11/01/2016
     &                    BETHE,STOPP,
CClark:end
     &                    WORK,LWORK)
CKeinSPASmehr
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate RPA, RPA(D), SOPPA and SOPPA(CCSD) transition
C              moments.
C
      use so_info, only: FN_RDENSE, FN_RDENSD
C
#include "implicit.h"
#include "priunit.h"
C
#include "cbiexc.h"
CSPAS:15/03-2006: merge with Dalton-2.0
C#include "cbilrs.h"
CKeinSPASmehr
#include "ccsdsym.h"
#include "ccorb.h"
#include "soppinf.h"
#include "mxcent.h"
#include "nuclei.h"
CClark:7/1/2016
#include "dummy.h"
#include "orgcom.h"
#include "iratdef.h"
CClark:end
C
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION TRLEN(3,NSYM,MXNEXI), TRVEL(3,NSYM,MXNEXI)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TQLEN(3,3,NSYM,MXNEXI), TQVEL(3,3,NSYM,MXNEXI)
      DIMENSION TTLEN(10,NSYM,MXNEXI)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI), TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI), EXENG(NSYM,MXNEXI)
CClark:11/01/2016
      DIMENSION BETHE(3,LQ),STOPP(3,LVEL)
CClark:end
      DIMENSION T2AM(LT2AM), WORK(LWORK)
C
      PARAMETER( ONE=1.0D0 )
C
      CHARACTER*8 LABEL, PDENS_LABEL
      CHARACTER*5 MODEL
      LOGICAL  IMAGPROP
C
CClark:7/1/2016
      REAL*8      ORSAVE(3)
      REAL*8      TRGOS(3,MXNEXI)
      CHARACTER*8 RTNLBL(2)
      CHARACTER*8 LABAPPG(6)
      INTEGER     LABSYMG(6)
      LOGICAL     FNDLAB
CClark:end
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TRMOM')
C
      WRITE(PDENS_LABEL,'(A7,I1)') 'EXCITA ', ISYMTR
      LURDENSE = -1
      LURDENSD = -1
      CALL GPOPEN(LURDENSE,FN_RDENSE,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      CALL GPOPEN(LURDENSD,FN_RDENSD,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
C
C-------------------------------
C     Loop over property labels.
C-------------------------------
C
      DO 200 IPRLBL = 1, NLBTOT
C
         LABEL  = LABAPP(IPRLBL)
         KSYMOP = LABSYM(IPRLBL)
C
CClark:7/1/2016
C--------------------------------------------------------
C     For this label, read integrals and process them.
C--------------------------------------------------------
C
         CALL SO_TRMOM_1(MODEL,LABEL,KSYMOP,ISYMTR,NEXCI,
     &                    T2AM,LT2AM,DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                    TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,TRLON,TRMAG,
     &                    BSRLON,EXENG,DUMMY,WORK,LWORK,
     &                    LURDENSE,LURDENSD)
CClark:end
C
  200 CONTINUE
C
CClark:7/1/2016
C---------------------------------------
C     For calculation of stopping power.
C---------------------------------------
C
      IF (STOPPW) THEN
C
C---------------------------------------
C     Open file AOGOS(integrals) and GOS
C---------------------------------------
C
         LUAOGOS = -1
         LUGOS   = -1
C
         CALL GPOPEN(LUAOGOS,'AOGOS   ','UNKNOWN',' ',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
         CALL GPOPEN(LUGOS,'GOS     ','UNKNOWN',' ',
     &               'FORMATTED',IDUMMY,.FALSE.)
         CALL GETDAT(RTNLBL(1),RTNLBL(2))
         RTNLBL(2)='SYMMETRI'
         CALL NEWLB2('EOFLABEL',RTNLBL,LUAOGOS,LUPRI)
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
         LPAOI  = NNBASX*6
C
         KDLAB  = 1
         KIDSYM = KDLAB  + 9*MXCENT
         KIDADR = KIDSYM + 9*MXCENT
         KPAOI  = KIDADR + 9*MXCENT
         KEND1  = KPAOI + LPAOI
         LWORK1 = LWORK - KEND1
C
         CALL SO_MEMMAX ('RP_TRMOM',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT ('RP_TRMOM',' ',KEND1,LWORK)
C
C---------------------------------------------------
C     Save origins and calculate domain of q values.
C---------------------------------------------------
C
         CALL DCOPY(3,ORIGIN,1,ORSAVE,1)
C
C--------------------------------------------
C     Loop over x,y,z components of q values.
C--------------------------------------------
C
         DO IQ = 1, LQ
C
            ORIGIN(1)=QMIN+(IQ-1)*QSTEP
            ORIGIN(2)=QMIN+(IQ-1)*QSTEP
            ORIGIN(3)=QMIN+(IQ-1)*QSTEP
C
C-------------------------------
CDeterminte Cos & Sin integrals.
C-------------------------------
C
            NCOMP   = 6
            NPATOM  = 0
            NLBTOTG = 0
C
            CALL GET1IN(WORK(KPAOI),'EXPIKR ',NCOMP,WORK(KEND1),
     &                  LWORK1,WORK(KDLAB),WORK(KIDSYM),
     &                  WORK(KIDADR),IDUMMY,.FALSE.,NPATOM,.TRUE.,
     &                  DUMMY,.FALSE.,DUMMY,IPRINT)
C
            IF (IPRSOP .GT. 5) THEN
               CALL AROUND('COS KX/K Integrals in atomic basis')
               CALL OUTPAK(WORK(KPAOI),NBAST,1,LUPRI)
               CALL AROUND('COS KY/K Integrals in atomic basis')
               CALL OUTPAK(WORK(KPAOI+NNBASX),NBAST,1,LUPRI)
               CALL AROUND('COS KZ/K Integrals in atomic basis')
               CALL OUTPAK(WORK(KPAOI+2*NNBASX),NBAST,1,LUPRI)
               CALL AROUND('SIN KX/K Integrals in atomic basis')
               CALL OUTPAK(WORK(KPAOI+3*NNBASX),NBAST,1,LUPRI)
               CALL AROUND('SIN KY/K Integrals in atomic basis')
               CALL OUTPAK(WORK(KPAOI+4*NNBASX),NBAST,1,LUPRI)
               CALL AROUND('SIN KZ/K Integrals in atomic basis')
               CALL OUTPAK(WORK(KPAOI+5*NNBASX),NBAST,1,LUPRI)
            END IF
C
C------------------------------------------
C     Set labels and symmetry of integrals.
C------------------------------------------
C
            NLAB = 6
C
            CALL LABCOP(NLAB,NLBTOTG,WORK(KDLAB),WORK(KIDSYM),
     &                  LABAPPG,LABSYMG)
C
C-------------------
C     Write integals
C-------------------
C                 CALL GETDAT(RTNLBL(1),RTNLBL(2))
C         Replace time information with symmetry information$
C                 RTNLBL(2)='SYMMETRI'
C
C         Place LUAOGOS in correct position for writing
            REWIND LUAOGOS
            IF (.NOT.FNDLAB('EOFLABEL',LUAOGOS)) THEN
               WRITE (LUPRI,'(/A)')
     &                ' End of file not found in AOGOS?'
               CALL QUIT('Internal error, EOF not found in AOGOS')
            END IF
            BACKSPACE LUAOGOS

C
            LEN = MAX (4,NNBASX)
            CALL NEWLB2(LABAPPG(1),RTNLBL,LUAOGOS,LUPRI)
            CALL WRITI (LUAOGOS,IRAT*LEN,WORK(KPAOI))
            CALL NEWLB2(LABAPPG(2),RTNLBL,LUAOGOS,LUPRI)
            CALL WRITI (LUAOGOS,IRAT*LEN,WORK(KPAOI+NNBASX*1))
            CALL NEWLB2(LABAPPG(3),RTNLBL,LUAOGOS,LUPRI)
            CALL WRITI (LUAOGOS,IRAT*LEN,WORK(KPAOI+NNBASX*2))
            CALL NEWLB2(LABAPPG(4),RTNLBL,LUAOGOS,LUPRI)
            CALL WRITI (LUAOGOS,IRAT*LEN,WORK(KPAOI+NNBASX*3))
            CALL NEWLB2(LABAPPG(5),RTNLBL,LUAOGOS,LUPRI)
            CALL WRITI (LUAOGOS,IRAT*LEN,WORK(KPAOI+NNBASX*4))
            CALL NEWLB2(LABAPPG(6),RTNLBL,LUAOGOS,LUPRI)
            CALL WRITI (LUAOGOS,IRAT*LEN,WORK(KPAOI+NNBASX*5))
C
C-----------------------------------------
C     We add an extra label to signify EOF
C-----------------------------------------
C
            CALL NEWLB2('EOFLABEL',RTNLBL,LUAOGOS,LUPRI)
C
C-----------------------------------
C     Set trgos into zeor before use
C-----------------------------------
C
            CALL DZERO(TRGOS, 3*MXNEXI)
C
C----------------------------------
C     Loop over cos/sin components.
C----------------------------------
C
C
            DO L=1,6
C
C--------------------------------------------------------------------
C     Process these integrals (but the workspace has to be rewritten?
C--------------------------------------------------------------------

               LABEL  = LABAPPG(L)
               KSYMOP = LABSYMG(L)
C
C-----------------------------------------------------------------
C        For this label, read integrals and process them.
C-----------------------------------------------------------------
C
               CALL SO_TRMOM_1(MODEL,LABEL,KSYMOP,ISYMTR,NEXCI,
     &                         T2AM,LT2AM,DENSIJ,LDENSIJ,
     &                         DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                         TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,TRLON,
     &                         TRMAG,BSRLON,EXENG,TRGOS,WORK,LWORK,
     &                         LURDENSE,LURDENSD)
C
            END DO

            CALL DSCAL(NEXCI*3,2.0/ORIGIN(1)**2.0,TRGOS,1)
C
C----------------------------------
C        Write GOS into LUGOS file.
C----------------------------------
C
            DO IEXCI = 1, NEXCI
C
               WRITE (LUGOS,'(/A,D16.13)')
     &                                'Number of q valus:', ORIGIN(1)
               WRITE (LUGOS,'(/A,I4,A,I4,A,I4)') 'ISYMTR=',ISYMTR,
     &                                     'NEXCI=',NEXCI,'IEXCI=',IEXCI
               WRITE (LUGOS,'(D20.13)') TRGOS(1,IEXCI)
               WRITE (LUGOS,'(D20.13)') TRGOS(2,IEXCI)
               WRITE (LUGOS,'(D20.13)') TRGOS(3,IEXCI)
C
C-------------------------------------------
C        Calculate Bethe and Stopping Power.
C-------------------------------------------
C
               BETHE(1,IQ) = BETHE(1,IQ) + TRGOS(1,IEXCI)
               BETHE(2,IQ) = BETHE(2,IQ) + TRGOS(2,IEXCI)
               BETHE(3,IQ) = BETHE(3,IQ) + TRGOS(3,IEXCI)
C
               CALL SO_STOPPW(STOPP,TRGOS,ISYMTR,IEXCI,EXENG,IQ)
C
C               STOPP(1)    = STOPP(1) + TRGOS(1,IEXCI)*QSTEP
C               STOPP(2)    = STOPP(2) + TRGOS(2,IEXCI)*QSTEP
C               STOPP(3)    = STOPP(3) + TRGOS(3,IEXCI)*QSTEP
C
            END DO
C
            IF ( IPRSOP .GE. 1) THEN
C
               DO IEXCI = 1,NEXCI
C
                  WRITE(LUPRI,'(/,1X,A,I3)')
     &            ' RPA generalized oscilator strength',
     &            ISYMTR
                  WRITE(LUPRI,9001)
                  WRITE(LUPRI,'(A)')
     &            ' Exci. component  Total'
                  WRITE(LUPRI,9002)
C
                  WRITE(LUPRI,'(1X,I3,3X,A,10X,D20.13)') IEXCI,'X',
     &                  TRGOS(1,IEXCI)
                  WRITE(LUPRI,'(1X,I3,3X,A,10X,D20.13)') IEXCI,'Y',
     &                  TRGOS(2,IEXCI)
                  WRITE(LUPRI,'(1X,I3,3X,A,10X,D20.13)') IEXCI,'Z',
     &                  TRGOS(3,IEXCI)
C
               END DO
C
            ENDIF
C
         END DO
C
C------------------------------------------
C     Load Origins back and close the file.
C------------------------------------------
C
         CALL DCOPY(3,ORSAVE,1,ORIGIN,1)
C
         CALL GPCLOSE(LUAOGOS,'KEEP')
         CALL GPCLOSE(LUGOS,'KEEP')
C
      ENDIF
C
      CALL GPCLOSE(LURDENSE,'KEEP')
      CALL GPCLOSE(LURDENSD,'KEEP')
CClark:end
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TRMOM')
C
      RETURN
C
 9001 FORMAT(1X,'================================================',
     &       '===================')
 9002 FORMAT(1X,'------------------------------------------------',
     &       '-------------------')
C
      END
