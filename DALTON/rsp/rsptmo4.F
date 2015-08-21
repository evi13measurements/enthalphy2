!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
C
      SUBROUTINE CRTMO(VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD,
     *                 RESVEC,CMO,UDV,PV,FOCK,FC,FV,
     *                 XINDX,MJWOP,WRK,LWRK)
C
#include "implicit.h"
#include "iratdef.h"
#include "dummy.h"
C
C PURPOSE:
C CALCULATION OF THIRD ORDER TRANSITION MOMENTS
C
      LOGICAL DOMOM, DIPLEN, ATEST
C
      CHARACTER*8 ALAB,BLAB,CLAB,DLAB
C
      DIMENSION VECA(*),VECB(*),VECC(*),VECD(*)
      DIMENSION VECBC(*),VECCD(*),VECBD(*)
      DIMENSION RESVEC(*)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*)
      DIMENSION XINDX(*),WRK(*)
C
      PARAMETER ( D0 = 0.0D0, ZEROTHR = 1.0D-10 )
      DIMENSION IEQTO(5)
C
#include "priunit.h"
#include "infrsp.h"
#include "maxorb.h"
#include "infvar.h"
#include "inforb.h"
#include "infpri.h"
#include "infspi.h"
#include "wrkrsp.h"
#include "tstjep.h"
#include "infhso.h"
#include "qrinf.h"
      DIMENSION MJWOP(2,MAXWOP,8)
#include "rspprp.h"
#include "indcr.h"
      REAL*8, ALLOCATABLE  ::  RESTOM(:,:,:,:,:), RESFRE(:,:)
      REAL*8, ALLOCATABLE  ::  RESELE(:,:,:,:,:), RESMAG(:,:,:,:,:)
      REAL*8, ALLOCATABLE  ::  RESQUA(:,:,:,:,:)
#include "infcr.h"
#include "inftmo.h"
#include "codata.h"
C
      CALL QENTER('CRTMO')
      ATEST = .FALSE.
C
C     We write the results to a complementary output file. This will then
C     both serve as a file for getting a summary of the results, but more
C     importantly, it will serve as a way of avoiding already calculated
C     gamma-components in case of a crashed calculation. Check for calculated
C     components are done in BCDCHK.
C
      LURSPRES = -1
      CALL GPOPEN(LURSPRES,'RESULTS.RSP','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      WRITE(LUPRI,'(//A,A)')
     *' <<<<< CALCULATING CONTRIBUTIONS TO THIRD ORDER',
     *' TRANSITION MOMENT >>>>>'
C
      IF (THREEPHOTON) THEN
         ALLOCATE ( RESTOM(3,3,3,MXEXCR,8), RESFRE(MXEXCR,8) )
         CALL DZERO(RESTOM,3*3*3*MXEXCR*8)
         CALL DZERO(RESFRE,MXEXCR*8)
      END IF
      IF (THREEPCD) THEN
         IF(.NOT.ALLOCATED(RESFRE)) ALLOCATE(RESFRE(MXEXCR,8))
         ALLOCATE ( RESELE(3,3,3,MXEXCR,8))
         ALLOCATE ( RESMAG(3,3,3,MXEXCR,8), RESQUA(6,3,3,MXEXCR,8))
         CALL DZERO(RESELE,3*3*3*MXEXCR*8)
         CALL DZERO(RESMAG,3*3*3*MXEXCR*8)
         CALL DZERO(RESQUA,6*3*3*MXEXCR*8)
         CALL DZERO(RESFRE,MXEXCR*8)
      END IF

C
      DO 200 ISYMD = 1,NSYM
      DO 300 ISYMC = 1,NSYM
      DO 400 ISYMB = 1,NSYM
C
      ISYMDX = MULD2H(IREFSY,ISYMD)
      ISYMA = MULD2H(ISYMD,MULD2H(ISYMC,ISYMB))
      IF ( (NTMCNV(ISYMD).GT.0) .AND. (NCTMOP(ISYMC).GT.0) .AND.
     *     (NBTMOP(ISYMB).GT.0) .AND. (NATMOP(ISYMA).GT.0) ) THEN
C
      DO 500 ID   = 1,NTMCNV(ISYMD)
C
      DO 600 ICOP = 1,NCTMOP(ISYMC)
      DO 650 ICFR = 1,NCTMFR
C
      DO 700 IBOP = 1,NBTMOP(ISYMB)
      DO 750 IBFR = 1,NBTMFR
C
C     If three photon absorption or harmonic generation calculation 
C     only certain frequencies are considered.
C
      IF ( (THREEPHOTON.OR.THREEPCD) .AND.
     &     ((ABS(EXCIT2(ISYMD,ID)-3*BTMFR(IBFR)).GT.ZEROTHR) .OR.
     &      (ABS(EXCIT2(ISYMD,ID)-3*CTMFR(ICFR)).GT.ZEROTHR)) ) GOTO 750
C
      IF (CTMOHG) THEN
         IF (IBFR.NE.ICFR) GO TO 750
      END IF
C
      DO 800 IAOP = 1,NATMOP(ISYMA)
C
C     Initialize variables.
C     Check if an equivalent moment calculation already has been done,
C     DOMOM indicates the result.
C     Read response vectors and eigen vectors from disk.
C     Check if some of the response vectors are equal or zero,
C     IBCDEQ indicates the result
C
      CALL BCDCHK(DOMOM,IBCDEQ,LURSPRES,DIPLEN,DUMMY,
     *            ISYMA,ISYMB,ISYMC,ISYMD,ISYMBC,ISYMBD,ISYMCD,
     *            ALAB,BLAB,CLAB,DLAB,IAOP,IBOP,ICOP,0,
     *            IBFR,ICFR,ID,FREQA,FREQB,FREQC,FREQD,
     *            KZYVA,KZYVB,KZYVC,KZYVD,KZYVBC,KZYVBD,KZYVCD,
     *            VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD,IEQTO)
C
      IF (THREEPHOTON .AND. (.NOT. DOMOM)) THEN
        CALL DIPLAB(ALAB,IA)
        CALL DIPLAB(BLAB,IB)
        CALL DIPLAB(CLAB,IC)
        RESFRE(ID,ISYMD)=EXCIT2(ISYMD,ID)
        RESTOM(IA,IB,IC,ID,ISYMD)=
     &  RESTOM(IEQTO(1),IEQTO(2),IEQTO(3),IEQTO(4),IEQTO(5))
      END IF
      IF (THREEPCD .AND. (.NOT. DOMOM)) THEN
        IF ((ALAB(2:4).EQ.'DIP').OR.(ALAB(2:7).EQ.'ANGMOM')) THEN
          CALL DIPLAB(ALAB,IA)
          CALL DIPLAB(BLAB,IB)
          CALL DIPLAB(CLAB,IC)
          RESFRE(ID,ISYMD)=EXCIT2(ISYMD,ID)
          IF(ALAB(2:4).EQ.'DIP') THEN
             RESELE(IA,IB,IC,ID,ISYMD)=
     &       RESELE(IEQTO(1),IEQTO(2),IEQTO(3),IEQTO(4),IEQTO(5))
          ELSE IF(ALAB(2:7).EQ.'ANGMOM') THEN
             RESMAG(IA,IB,IC,ID,ISYMD)=
     &       RESMAG(IEQTO(1),IEQTO(2),IEQTO(3),IEQTO(4),IEQTO(5))
          END IF
        END IF
        IF (ALAB(3:8).EQ.'ROTSTR'.OR.ALAB(3:8).EQ.'SECMOM') THEN
          CALL QUADLAB(ALAB,IA)
          CALL DIPLAB(BLAB,IB)
          CALL DIPLAB(CLAB,IC)
          RESFRE(ID,ISYMD)=EXCIT2(ISYMD,ID)
          RESQUA(IA,IB,IC,ID,ISYMD)=
     &    RESQUA(IEQTO(1),IEQTO(2),IEQTO(3),IEQTO(4),IEQTO(5))
        END IF

      END IF
      IF (.NOT.DOMOM) GOTO 800
C
C    Initialize third order transition moment
C
      TMOM = 0
C
      IF (IPRRSP.GT.0) WRITE(LUPRI,'(//A15,2A20,/A)')
     *   'Contribution','Term','Accumulated',
     *   ' ------------------------------------------------------'
C
C
C     Calculate Na T[4] Nb Nc Nd
C
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL T4DRV(IBCDEQ,ISYMA,ISYMB,ISYMC,ISYMD,VECA,VECB,VECC,VECD,
     *           -FREQB,-FREQC,-FREQD,XINDX,UDV,PV,MJWOP,
     *           WRK,LWRK,CMO,FC)
      VAL = -DDOT(KZYVA,WRK,1,VECA,1)
      TMOM = TMOM + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A17,F18.8,F20.8)')' Na T[4] Nb Nc Nd',VAL,TMOM
C
      IF (IPRRSP.GT.5) CALL TIMER('T4DRV ',TIMSTR,TIMEND)
C
C
C     Calculate Na T[3] Nb Ncd type terms (three permutations)
C
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL T3DRV(1,ISYMA,ISYMB,ISYMCD,VECB,VECCD,ATEST,VECA,
     *           -FREQB,-FREQC-FREQD,XINDX,UDV,PV,MJWOP,
     &           WRK,LWRK,CMO,FC,FV)
      VAL = DDOT(KZYVA,WRK,1,VECA,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      IF (IBCDEQ.EQ.2) THEN
         TMPVAL = TMPVAL + VAL
         TMOM = TMOM + VAL
         CALL T3DRV(1,ISYMA,ISYMD,ISYMBC,VECD,VECBC,ATEST,VECA,
     *              -FREQD,-FREQB-FREQC,XINDX,UDV,PV,MJWOP,
     &              WRK,LWRK,CMO,FC,FV)
         VAL = DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         TMOM = TMOM + VAL
      ELSE
         CALL T3DRV(1,ISYMA,ISYMC,ISYMBD,VECC,VECBD,ATEST,VECA,
     *              -FREQC,-FREQB-FREQD,XINDX,UDV,PV,MJWOP,
     *              WRK,LWRK,CMO,FC,FV)
         VAL = DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         TMOM = TMOM + VAL
         CALL T3DRV(1,ISYMA,ISYMD,ISYMBC,VECD,VECBC,ATEST,VECA,
     *              -FREQD,-FREQB-FREQC,XINDX,UDV,PV,MJWOP,
     *              WRK,LWRK,CMO,FC,FV)
         VAL = DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         TMOM = TMOM + VAL
      END IF
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)')' Na T[3] Nx Nyz',TMPVAL,TMOM
C
      IF (IPRRSP.GT.5) CALL TIMER('T3DRV ',TIMSTR,TIMEND)

C
C     Calculate Na B[3] Nc Nd type terms 
C     (two of four permutations in each call)
C
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL X3INIT(KZYVA,KZYVC,KZYVD,ISYMA,ISYMC,ISYMD,BLAB,
     *            ISYMB,VECC,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      CALL X3INIT(KZYVA,KZYVB,KZYVD,ISYMA,ISYMB,ISYMD,CLAB,
     *            ISYMC,VECB,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVA,RESVEC,1,VECA,1)
      TMOM = TMOM + VAL
      TMPVAL = TMPVAL + VAL
C
      IF (IPRRSP.GT.0)
     *  WRITE(LUPRI,'(A15,3F20.8)')' Na X[3] Ny Nz ',TMPVAL,TMOM
C
C
C     Calculate Nb A[3] Nc Nd type terms
C     (two of six permutations in each call)
C
C
      CALL A3INIT(KZYVB,KZYVC,KZYVD,ISYMB,ISYMC,ISYMD,ALAB,
     *            ISYMA,VECC,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVB,RESVEC,1,VECB,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      CALL A3INIT(KZYVC,KZYVB,KZYVD,ISYMC,ISYMB,ISYMD,ALAB,
     *            ISYMA,VECB,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVC,RESVEC,1,VECC,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A3INIT(KZYVD,KZYVB,KZYVC,ISYMD,ISYMB,ISYMC,ALAB,
     *            ISYMA,VECB,VECC,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVD,RESVEC,1,VECD,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      IF (IPRRSP.GT.0)
     *  WRITE(LUPRI,'(A15,2F20.8)')' Nx A[3] Ny Nz ',TMPVAL,TMOM
C
C
C     Calculate Na B[2] Ncd type terms (two permutations)
C
C
      CALL X2INIT(1,KZYVA,KZYVCD,ISYMA,ISPINA,ISYMCD,0,WRK(1),VECCD,
     *            RESVEC,XINDX,UDV,PV,BLAB,ISYMB,0,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      CALL X2INIT(1,KZYVA,KZYVBD,ISYMA,ISPINA,ISYMBD,0,1,VECBD,
     *            RESVEC,XINDX,UDV,PV,CLAB,ISYMC,ISPINC,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)') ' Na X[2] Nyz   ',TMPVAL,TMOM
C
C
C     Calculate Nb A[2] Ncd type terms (six permutations)
C
C
      CALL A2INIT(1,KZYVB,KZYVCD,ISYMB,ISPINB,ISYMCD,0,1,VECCD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVB,RESVEC,1,VECB,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVCD,KZYVB,ISYMCD,0,ISYMB,ISPINB,1,VECB,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVCD,RESVEC,1,VECCD,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVC,KZYVBD,ISYMC,ISPINC,ISYMBD,0,1,VECBD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVC,RESVEC,1,VECC,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVBD,KZYVC,ISYMBD,0,ISYMC,ISPINC,1,VECC,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVBD,RESVEC,1,VECBD,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVD,KZYVBC,ISYMD,ISPIND,ISYMBC,0,1,VECBC,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVD,RESVEC,1,VECD,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVBC,KZYVD,ISYMBC,0,ISYMD,ISPIND,1,VECD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVBC,RESVEC,1,VECBC,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)') ' Nx A[2] Nyz   ',TMPVAL,TMOM
C
      IF (IPRRSP.GT.5) CALL TIMER('OTHERS',TIMSTR,TIMEND)
C
      WRITE(LUPRI,'(/A,3(/A,A10,I4,F10.6))')
     *     '@ Third order transition moment in a.u. for',
     *     '@ A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *     '@ B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB,
     *     '@ C operator, symmetry, frequency: ',CLAB,ISYMC,-FREQC
      WRITE(LUPRI,'(/A,2I4,F10.6//A,F20.8/)')
     *     '@ State no., symmetry, excitation energy:',ID,ISYMD,FREQD,
     &     '@ < 0 | ABC | f >  = ', TMOM
C
C     Write out to result file
C
      WRITE(LURSPRES,'(/A,3(/A,A10,I4,F10.6))')
     *     ' Third order transition moment in a.u. for',
     *     ' A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *     ' B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB,
     *     ' C operator, symmetry, frequency: ',CLAB,ISYMC,-FREQC
      WRITE(LURSPRES,'(/A,2I4,F10.6)')
     *     ' State no., symmetry, excitation energy:',ID,ISYMD,FREQD
      WRITE(LURSPRES,'(/A,F20.8)') ' < 0 | ABC | f >  = ', TMOM
C
      IF (THREEPHOTON) THEN
        CALL DIPLAB(ALAB,IA)
        CALL DIPLAB(BLAB,IB)
        CALL DIPLAB(CLAB,IC)
        RESTOM(IA,IB,IC,ID,ISYMD)=TMOM
        RESFRE(ID,ISYMD)=EXCIT2(ISYMD,ID)
      END IF
      IF (THREEPCD) THEN
        IF ((ALAB(2:4).EQ.'DIP').OR.(ALAB(2:7).EQ.'ANGMOM')) THEN
          CALL DIPLAB(ALAB,IA)
          CALL DIPLAB(BLAB,IB)
          CALL DIPLAB(CLAB,IC)
          RESFRE(ID,ISYMD)=EXCIT2(ISYMD,ID)
          IF(ALAB(2:4).EQ.'DIP') THEN
             RESELE(IA,IB,IC,ID,ISYMD)=TMOM
          ELSE IF(ALAB(2:7).EQ.'ANGMOM') THEN
             RESMAG(IA,IB,IC,ID,ISYMD)=TMOM
          END IF
        END IF
        IF ((ALAB(3:8).EQ.'ROTSTR').OR.(ALAB(3:8).EQ.'SECMOM')) THEN
          CALL QUADLAB(ALAB,IA)
          CALL DIPLAB(BLAB,IB)
          CALL DIPLAB(CLAB,IC)
          RESFRE(ID,ISYMD)=EXCIT2(ISYMD,ID)
          RESQUA(IA,IB,IC,ID,ISYMD)=TMOM
        END IF
       
      END IF
C
 800  CONTINUE
 750  CONTINUE
 700  CONTINUE
 650  CONTINUE
 600  CONTINUE
 500  CONTINUE
C
      END IF
C
 400  CONTINUE
 300  CONTINUE
 200  CONTINUE
      CALL GPCLOSE(LURSPRES,'KEEP')
C
C
      IF (THREEPHOTON) THEN
         CALL PRINT_THREEPHOTON(RESTOM,RESFRE,'DIPLEN',3)   
         DEALLOCATE (RESTOM)
         IF(.NOT.THREEPCD) DEALLOCATE(RESFRE)
      END IF
      IF (THREEPCD) THEN
        IF(.NOT.THREEPHOTON) THEN
         CALL PRINT_THREEPHOTON(RESELE,RESFRE,'DIPVEL',3)
         CALL PRINT_THREEPHOTON(RESMAG,RESFRE,'ANGMOM',3)
         CALL PRINT_THREEPHOTON(RESQUA,RESFRE,'ROTSTR',6)

         CALL TITLER('CALCULATING 3PCD IN VELOCITY GAUGE',
     &       '*',112)

         CALL CALC_3PCD(RESELE,RESMAG,RESQUA,NTMCNV,NSYM,MXEXCR,LUPRI)
        ELSE IF(THREEPHOTON) THEN
         CALL PRINT_THREEPHOTON(RESMAG,RESFRE,'ANGMOM',3)
         CALL PRINT_THREEPHOTON(RESQUA,RESFRE,'SECMOM',6)

         CALL TITLER('CALCULATING 3PCD IN LENGTH GAUGE',
     &       '*',112)

         CALL CALC_3PCD(RESELE,RESMAG,RESQUA,NTMCNV,NSYM,MXEXCR,LUPRI)
        END IF
         DEALLOCATE(RESELE,RESMAG,RESQUA,RESFRE)
      END IF
C
C    End of subroutine CRTMO
C
      CALL QEXIT('CRTMO')
      RETURN
      END
      SUBROUTINE PRINT_THREEPHOTON(RESTOM,RESFRE,TYPE,NASIZE)
C
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "rspprp.h"
#include "indcr.h"
#include "infcr.h"
#include "inftmo.h"
#include "codata.h"
C
      INTEGER NASIZE
      PARAMETER ( D0 = 0.0D0, ZERO = 1.0D-10 )
      DIMENSION RESTOM(NASIZE,3,3,MXEXCR,8), RESFRE(MXEXCR,8)
      CHARACTER*8 ALAB,BLAB,CLAB,DLAB,TYPE
      CHARACTER*1 ITINDX,JTINDX
      CHARACTER*2 ITINDX2
      CHARACTER*5 PRINTSTRING
      INTEGER     ILENGTHA

      CALL TITLER('FINAL RESULTS FROM THREE-PHOTON CALCULATION',
     &       '*',112)

      ILENGTHA=3

      IF(TYPE(1:6).EQ.'DIPLEN') THEN
      WRITE(LUPRI,'(A64,4(/A64))')
     & ' The three-photon absorption strength for an average molecular',
     & ' orientation is computed according to formulas given by       ',
     & ' W.M. McClain in J. Chem. Phys. 57:2264, 1972. The absorption ',
     & ' depends on the light polarization. A monochromatic light     ',
     & ' source is assumed.                                           '
      WRITE(LUPRI,'(2(/A64))')
     & ' All results are presented in atomic units, except the        ',
     & ' excitation energy which is given in eV.                      '
      WRITE(LUPRI,'(/,3(/A52))')
     &        '+----------------------------------+',
     &        '| Three-photon transition tensor T |',
     &        '+----------------------------------+'
      WRITE(LUPRI,'(2A)')
     & ' ----------------------------------------------------',
     & '--------------'
      PRINTSTRING='  T_{'
      ELSE IF(TYPE(1:6).EQ.'DIPVEL') THEN
            WRITE(LUPRI,'(A64,4(/A64))')
     & ' The three-photon transition moment tensor in velocity gauge, ',
     & ' necessary for the calculation of three-photon circular       ',
     & ' dichroism (experimental prelimibary code so far)             ',
     & '                                                              ',
     & ' Daniel H. Friese, 4. 8. 2015                                 '
      WRITE(LUPRI,'(/,3(/A52))')
     &        '+----------------------------------+',
     &        '| Three-photon transition tensor P |',
     &        '+----------------------------------+'
      WRITE(LUPRI,'(2A)')
     & ' ----------------------------------------------------',
     & '--------------'
      PRINTSTRING='  P_{'
      ELSE IF(TYPE(1:6).EQ.'ANGMOM') THEN
            WRITE(LUPRI,'(A64,4(/A64))')
     & ' The three-photon transition moment tensor                    ',
     & ' constructed from the magnetic dipole moment operator and     ',
     & ' twice the electric dipole operator in velocity gauge.        ',
     & ' Necessary for the calculation of three-photon circular       ',
     & ' dichroism (experimental prelimibary code so far)             ',
     & '                                                              ',
     & ' Daniel H. Friese, 4. 8. 2015                                 '
      WRITE(LUPRI,'(/,3(/A52))')
     &        '+----------------------------------+',
     &        '| Three-photon transition tensor M |',
     &        '+----------------------------------+'
      WRITE(LUPRI,'(2A)')
     & ' ----------------------------------------------------',
     & '--------------------'
      PRINTSTRING='  M_{'
      ELSE IF(TYPE(1:6).EQ.'ROTSTR') THEN
            WRITE(LUPRI,'(A64,4(/A64))')
     & ' The three-photon transition moment tensor in velocity gauge, ',
     & ' constructed from the electric quadrupole moment operator and ',
     & ' twice the electric dipole operator in velocity gauge.        ',
     & ' Necessary for the calculation of three-photon circular       ',
     & ' dichroism (experimental prelimibary code so far)             ',
     & '                                                              ',
     & ' Daniel H. Friese, 4. 8. 2015                                 '
      WRITE(LUPRI,'(/,3(/A52))')
     &        '+----------------------------------+',
     &        '| Three-photon transition tensor T |',
     &        '+----------------------------------+'
      WRITE(LUPRI,'(2A)')
     & ' ----------------------------------------------------',
     & '--------------------'
      PRINTSTRING=' QP_{'
      ILENGTHA=6
      ELSE IF(TYPE(1:6).EQ.'SECMOM') THEN
            WRITE(LUPRI,'(A64,4(/A64))')
     & ' The three-photon transition moment tensor in length gauge,   ',
     & ' constructed from the electric quadrupole moment operator and ',
     & ' twice the electric dipole operator in velocity gauge.        ',
     & ' Necessary for the calculation of three-photon circular       ',
     & ' dichroism (experimental prelimibary code so far)             ',
     & '                                                              ',
     & ' Daniel H. Friese, 4. 8. 2015                                 '
      WRITE(LUPRI,'(/,3(/A52))')
     &        '+----------------------------------+',
     &        '| Three-photon transition tensor Q |',
     &        '+----------------------------------+'
      WRITE(LUPRI,'(2A)')
     & ' ----------------------------------------------------',
     & '--------------------'
      PRINTSTRING=' QL_{'
      ILENGTHA=6
      END IF



C
      DO ISYMD=1,NSYM
      DO ID=1,NTMCNV(ISYMD)
C
        DF=D0
        DG=D0
        WRITE(LUPRI,'(A12,I2,A9,I2,A11,F8.3,A3)') '   Symmetry ',ISYMD,
     &  '   State ',ID,'   Energy: ',RESFRE(ID,ISYMD)*XTEV,' eV'
        WRITE(LUPRI,'(A1)')
        DO I=1,ILENGTHA
        DO J=1,3
        IF(ILENGTHA.EQ.3) THEN
            CALL LABDIP(ITINDX,I)
            CALL LABDIP(JTINDX,J)
          WRITE(LUPRI,'(3(A5,2A1,A3,F12.5))')
     &   PRINTSTRING,ITINDX,JTINDX,'x}=',RESTOM(I,J,1,ID,ISYMD),
     &   PRINTSTRING,ITINDX,JTINDX,'y}=',RESTOM(I,J,2,ID,ISYMD),
     &   PRINTSTRING,ITINDX,JTINDX,'z}=',RESTOM(I,J,3,ID,ISYMD)
        ELSE IF(ILENGTHA.EQ.6) THEN
            CALL LABQUA(ITINDX2,I)
            CALL LABDIP(JTINDX,J)
          WRITE(LUPRI,'(3(A5,A2,A1,A3,F12.5))')
     &   PRINTSTRING,ITINDX2,JTINDX,'x}=',RESTOM(I,J,1,ID,ISYMD),
     &   PRINTSTRING,ITINDX2,JTINDX,'y}=',RESTOM(I,J,2,ID,ISYMD),
     &   PRINTSTRING,ITINDX2,JTINDX,'z}=',RESTOM(I,J,3,ID,ISYMD)
        ELSE
         CALL QUIT('Weird value of ilengtha in print_threephoton')
        END IF
        END DO
        END DO
      WRITE(LUPRI,'(2A)')
     & ' ----------------------------------------------------',
     & '--------------------'
      END DO
      END DO

      IF(TYPE(1:6).EQ.'DIPLEN') THEN
      WRITE(LUPRI,'(/,6(/A60))')
     &        ' Transition probabilities (a.u.)         ',
     &        '--------------------------------------------------',
     &        ' D  =  (3*Df + 2*Dg)/35, Linear   polarization',
     &        ' D  = (-3*Df + 5*Dg)/35, Circular polarization',
     &        ' Df = sum(i,j,k){ T_iij * T_kkj }             ',
     &        ' Dg = sum(i,j,k){ T_ijk * T_ijk }             '
      WRITE(LUPRI,'(3(/A53))')
     &        '       Polarization ratio      ',
     &        '-------------------------------',
     &        '    R  = (-3*Df+5*Dg)/(3*Df+2*Dg)  '
      WRITE(LUPRI,'(/,3(/A56))')
     &        '+-------------------------------------+',
     &        '| Three-photon transition probability |',
     &        '+-------------------------------------+'
      WRITE(LUPRI,'(A62,A)')
     &'-----------------------------------------------------------',
     &'-----------'
      WRITE(LUPRI,'(A6,A4,A8,A14,3A11,A8/A62,A)') 'Sym','No',
     &        'Energy','Polarization','Df','Dg','D','R',
     &'-----------------------------------------------------------',
     &'-----------'
      DO ISYMD=1,NSYM
         DO ID = 1,NTMCNV(ISYMD)
            DF=0.0D0
            DG=0.0D0
            DO I=1,3
            DO J=1,3
            DO K=1,3
               DF=DF+RESTOM(I,I,J,ID,ISYMD)*RESTOM(K,K,J,ID,ISYMD)
               DG=DG+RESTOM(I,J,K,ID,ISYMD)**2
            END DO
            END DO
            END DO
            R=(-3*DF+5*DG)/(3*DF+2*DG)
            D=(3*DF+2*DG)/35
            WRITE(LUPRI,'(A2,2I4,F8.2,A14,3E11.3,F8.2)')' ',ISYMD,ID,
     &           RESFRE(ID,ISYMD)*XTEV,'Linear     ',DF,DG,D,R
            D=(-3*DF+5*DG)/35
            WRITE(LUPRI,'(A2,2I4,F8.2,A14,3E11.3,F8.2)')' ',ISYMD,ID,
     &              RESFRE(ID,ISYMD)*XTEV,'Circular   ',DF,DG,D,R
         END DO
      END DO
      WRITE(LUPRI,'(A62,A)')
     &'-----------------------------------------------------------',
     &'-----------'
      END IF
      RETURN
      END
      SUBROUTINE LABDIP(LAB,I)
C
#include "implicit.h"
C
      CHARACTER*1 LAB
C
C Map integer numbers to dipole operators
C 1 => x, 2 => y, and 3 => z
C
      IF (I.EQ.1) LAB(1:1)='x'
      IF (I.EQ.2) LAB(1:1)='y'
      IF (I.EQ.3) LAB(1:1)='z'
C
      RETURN
      END

      SUBROUTINE LABQUA(LAB,I)
C
#include "implicit.h"
C
      CHARACTER*2 LAB
C
C Map integer numbers to dipole operators
C 1 => x, 2 => y, and 3 => z
C
      IF (I.EQ.1) LAB(1:2)='xx'
      IF (I.EQ.2) LAB(1:2)='yy'
      IF (I.EQ.3) LAB(1:2)='zz'
      IF (I.EQ.4) LAB(1:2)='xy'
      IF (I.EQ.5) LAB(1:2)='xz'   
      IF (I.EQ.6) LAB(1:2)='yz'
C
      RETURN
      END

      SUBROUTINE CALC_3PCD(RESELE,RESMAG,RESQUA,NSTATE,NSYM,MXEXCR,
     &                     LUPRI)
C
      implicit none
C
      INTEGER NSYM,MXEXCR,LUPRI
      INTEGER NSTATE(NSYM)
      INTEGER I,J,K,L,M, ISYM, ISTATE
      REAL*8  MOM3PCD,MPABC,MPAAB1,MPAAB2,MPAAA,QPABC,QPAAB
      REAL*8  RESELE(3,3,3,MXEXCR,NSYM),
     &        RESMAG(3,3,3,MXEXCR,NSYM),
     &        RESQUA(6,3,3,MXEXCR,NSYM)
      REAL*8,ALLOCATABLE :: TPLUSTENSOR(:,:,:)
C
        WRITE(LUPRI,'(A64,18(/A64))')
     & ' Preliminary output for three-photon circular dichroism.      ',
     & '                                                              ',
     & ' The following intermediates in a.u. are being printed:       ',
     & '                                                              ',
     & ' B1: Magnetic contribution, M_abc * P_abc                     ',
     & ' B2: Magnetic contribution, M_a(p)a(p)m * P_bcc               ',
     & '       (electric dipole indices on magnetic tensor are paired)',
     & ' B3: Magnetic contribution, M_a(m)a(p)b * P_bcc               ',
     & '       (electric and magnetic index paired in magnetic tensor)',
     & ' B4: Magnetic contribution, M_aaa * P_bbb                     ',
     & ' B5: Quadrupole contribution, Q_abc * P_abc                   ',
     & ' B6: Quadrupole contribution, Q_{a(q)a(p)b * P_bcc            ',
     & '       (electric and quadrupole index paired in Q-tensor)     ',
     & '                                                              ',
     & ' These intermediates can be combined to 3PCD in manner        ',
     & ' depending on polarization and propagation of the incident    ',
     & ' photons.                                                     ',
     & '                                                              ',
     & ' Daniel H. Friese, 7. 8. 2015                                 ' 
          WRITE(LUPRI,'(/,X,60A)')
     &'+----------------------------------------------------------+'
          WRITE(LUPRI,'(X,60A)')
     &'|                    Contributions to 3PCD                 |'
          WRITE(LUPRI,'(X,60A)')
     &'+----------------------------------------------------------+'

      DO ISYM=1, NSYM
        DO ISTATE=1, NSTATE(ISYM)
          ! MAGNETIC TERMS
          MPABC =0.0D0
          MPAAB1=0.0D0
          MPAAB2=0.0D0
          MPAAA =0.0D0
          QPABC =0.0D0
          QPAAB =0.0D0
          ALLOCATE(TPLUSTENSOR(3,3,3))
          CALL CONTRACT_QUADRUTENS_RANK_3(RESQUA(1,1,1,ISTATE,ISYM),
     &            TPLUSTENSOR)
          DO I=1, 3
            DO J=1, 3
              ! B4: Magnetic, diagnoal elements
              MPAAA=MPAAA+RESMAG(I,I,I,ISTATE,ISYM)*
     &              RESELE(J,J,J,ISTATE,ISYM)
              DO K=1, 3
                 ! B1: Magnetic, all indices different
                 MPABC  = MPABC+RESMAG(I,J,K,ISTATE,ISYM)* 
     &                    RESELE(I,J,K,ISTATE,ISYM)
                 ! B2: Magnetic, two dipole indices paired
                 MPAAB1 = MPAAB1+RESMAG(J,I,I,ISTATE,ISYM)*
     &                    RESELE(J,K,K,ISTATE,ISYM)
                 ! B3: Magnetic, magnetic and dipole index paired
                 MPAAB2 = MPAAB2+RESMAG(I,I,J,ISTATE,ISYM)*
     &                    RESELE(J,K,K,ISTATE,ISYM)
                 ! B5: Quadrupole, all indices different
                 QPABC  = QPABC+TPLUSTENSOR(I,J,K)*
     &                    RESELE(I,J,K,ISTATE,ISYM)
                 ! B6: Quadrupole, quad. and dipole index paired
                 QPAAB  = QPAAB+TPLUSTENSOR(I,J,I)*
     &                    RESELE(J,K,K,ISTATE,ISYM)
              END DO
            END DO
          END DO
          DEALLOCATE(TPLUSTENSOR)
          ! PRINT
!          WRITE(LUPRI,*)'MPABC',MPABC,'MPAAB1:',MPAAB1,'MPAAB2:',
!     &        MPAAB2,'MPAAA:',MPAAA
          WRITE(LUPRI,'(X,2(A9,I4),A34)')
     &'| Irrep: ',ISYM,'  State: ',ISTATE,'                   |'
          WRITE(LUPRI,'(X,A60)')
     &'| Magnetic terms:                                          |'
!          WRITE(LUPRI,'(A7,F11.7)')'| B1: ',MPABC
           WRITE(LUPRI,'(X,2(A6,F23.7),A2)')
     &'| B1: ',MPABC,'  B2: ',MPAAB1,' |'
           WRITE(LUPRI,'(X,2(A6,F23.7),A2)')
     &'| B3: ',MPAAB2,'  B4: ',MPAAA,' |'
           WRITE(LUPRI,'(X,A60)')
     &'| Quadrupole terms:                                        |'
           WRITE(LUPRI,'(X,2(A6,F23.7),A2)')
     &'| B5: ',QPABC,' B6: ',QPAAB,' |'
           WRITE(LUPRI,'(X,A60)')
     &'+----------------------------------------------------------+'
        END DO
      END DO

      RETURN
      END SUBROUTINE CALC_3PCD

      SUBROUTINE CONTRACT_QUADRUTENS_RANK_3(ORIGINALTENSOR,SEEKEDTENSOR)

      implicit none

      INTEGER IDX1,IDX2,IDX3 ! Indices of the seeked tensor
      INTEGER IPLUS(3),IMINUS(3),QUADINDEXPLUS,QUADINDEXMINUS
      CHARACTER*2 QUADLABELPLUS,QUADLABELMINUS
      CHARACTER*1 LABMINUS1,LABMINUS2,LABPLUS1,LABPLUS2
      REAL*8 ORIGINALTENSOR(6,3,3), SEEKEDTENSOR(3,3,3)

      DO IDX1 = 1, 3
        DO IDX2 = 1, 3
          DO IDX3 = 1,3
            IF(IDX2.EQ.1) THEN
              IPLUS(1)=1
              IPLUS(2)=2
              IPLUS(3)=3
              IMINUS(1)=1
              IMINUS(2)=3
              IMINUS(3)=2
            ELSE IF(IDX2.EQ.2) THEN
              IPLUS(1)=2
              IPLUS(2)=3
              IPLUS(3)=1
              IMINUS(1)=2
              IMINUS(2)=1
              IMINUS(3)=3
            ELSE IF(IDX2.EQ.3) THEN
              IPLUS(1)=3
              IPLUS(2)=1
              IPLUS(3)=2
              IMINUS(1)=3
              IMINUS(2)=2
              IMINUS(3)=1    
            END IF
            CALL LABDIP(LABMINUS1,IDX1)
            CALL LABDIP(LABMINUS2,IMINUS(2))    
            CALL LABDIP(LABPLUS1,IDX1)
            CALL LABDIP(LABPLUS2,IPLUS(2))
            IF(IDX1.GT.IMINUS(2)) THEN
               QUADLABELMINUS(1:1)=LABMINUS2
               QUADLABELMINUS(2:2)=LABMINUS1
            ELSE
               QUADLABELMINUS(1:1)=LABMINUS1
               QUADLABELMINUS(2:2)=LABMINUS2
            END IF       
            IF(IDX1.GT.IPLUS(2)) THEN
               QUADLABELPLUS(1:1)=LABPLUS2
               QUADLABELPLUS(2:2)=LABPLUS1
            ELSE
               QUADLABELPLUS(1:1)=LABPLUS1
               QUADLABELPLUS(2:2)=LABPLUS2
            END IF  
            CALL QUADLABLOWER(QUADLABELPLUS,QUADINDEXPLUS)
            CALL QUADLABLOWER(QUADLABELMINUS,QUADINDEXMINUS)
            SEEKEDTENSOR(IDX1,IDX2,IDX3) = 
     &      ORIGINALTENSOR(QUADINDEXPLUS,IPLUS(3),IDX3) -
     &      ORIGINALTENSOR(QUADINDEXMINUS,IMINUS(3),IDX3)
          END DO
        END DO
      END DO
 
      END SUBROUTINE
