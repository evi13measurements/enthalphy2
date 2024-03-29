!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
C File: sirius/sirpop.F
C
C
C  /* Deck sirpop */
      SUBROUTINE SIRPOP(INP_POPKEY,DV,WRK,LFREE_INPUT)
C
C Major revision March 2011 Hans Joergen Aa. Jensen
C Old version: Nov-10-85 Hans Agren, revised 910708-hjaaj/0710-hjaaj
C
C Purpose :
C    To control calculations of certain properties from SIRIUS wave function
C
C MOTECC-90: This module, SIRPOP, is described in the input/ouput
C            Documentation of MOTECC-90.
C
#include "implicit.h"
      CHARACTER*(*)  INP_POPKEY
      DIMENSION      DV(*), WRK(*)
C
C   INFINP : FLAG(*)
C   INFIND : ISX()
C   INFORB ? NBAST,
C   INFPOP ? a lot
C   INFTAP : LUIT1
C
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "priunit.h"
#include "infinp.h"
#include "infind.h"
#include "inforb.h"
#include "infpop.h"
#include "inftap.h"
#include "infpri.h"
! nuclei: NATOMS
#include "nuclei.h"
C
C
      LOGICAL DIPOLM, QUADRM, MULLKM
      LOGICAL FNDLAB
      CHARACTER*(6) POPKEY

#ifdef SIRPOP_DEBUG
      logical first_call
      data first_call/.true./
#endif
C
C
      CALL QENTER('SIRPOP')
C 
      POPKEY = INP_POPKEY
      IPRMUL_SAVE  = IPRMUL
      IPOPTYP_SAVE = IPOPTYP
      IF (IPOPTYP .LE. 0) IPOPTYP = 1
      IF (POPKEY(1:4) .EQ. 'DIIS') THEN
         iprmul = -lim_poppri
         DIPOLM = .false.
         QUADRM = .false.
         MULLKM = .true.
      ELSE IF (POPKEY(1:6) .EQ. 'MCITER') THEN
         iprmul = -lim_poppri
         DIPOLM = .false.
         QUADRM = .false.
         MULLKM = .true.
      ELSE IF (POPKEY(1:5) .EQ. 'FINAL') THEN
!        WRITE (LUPRI,'(//A/)')
!    &   ' --- SIRIUS Population Analysis. ---'
         CALL TITLER('--- SIRIUS Population Analysis ---',' ',200)
         DIPOLM = .NOT.FLAG(71)
         QUADRM = .NOT.FLAG(72)
         MULLKM = .NOT.FLAG(73)
      ELSE
         CALL QUIT('Unrecognized POPKEY : '//POPKEY)
      END IF

#ifdef SIRPOP_DEBUG
      if (first_call) then
         write (lupri,'(//A)')
     &   '  *** TEST DEBUG OUTPUT SIRPOP FIRST CALL ***'
         iprmul = 21 ! DEBUG !!!!
         first_call = .false.
      end if
#endif

C
C  *****   core allocations  *****
C
C     ***---VIRIAL---POPTRA---***
C
C     The structure of the working area will be
C
C       1) CMO      MO coefficients
C       2) OCCUP    occupation numbers
C       3) T or S   kinetic energy matrix, overlap matrix
C       4) CMOCGT   MO coefficients in cgto basis
C       5) SS       overlap matrix in cgto basis
C
      KFREE_INPUT = 1
      KFREE = KFREE_INPUT
      LFREE = LFREE_INPUT
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KOCCUP,NORBT ,WRK,KFREE,LFREE)
C
      CALL MEMGET('REAL',KOVLP ,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMOCG,NBAST*NOCCT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSS   ,N2BASX,WRK,KFREE,LFREE)
C
C     ***---MULLIKEN---***
C
C     *)  POP***   matrices for accumulation of contributions
C                  to various forms of population analysis
C
      CALL MEMGET('REAL',KPOPNAB,NATOMS*NATOMS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPOPABI,NATOMS*NATOMS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPOPCG ,NBAST ,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPOPC  ,NBAST*NOCCT,WRK,KFREE,LFREE)
!     reuse allocation for "OVLP" for VIRIAL and MULLIKEN to "POPNET" in MULLIKEN
      KPOPNET = KOVLP
C
C     ***** Get the MO coefficients and NO occupation numbers
C
      REWIND LUIT1
      CALL READMO (WRK(KCMO),9)
      REWIND LUIT1
      IF ( FNDLAB('NATOCC  ',LUIT1) ) THEN
         CALL READT(LUIT1,NORBT,WRK(KOCCUP))
!        IF (INDEX(WFTYPE,'MP2') .GT. 0) THEN
! hjaaj: NOCCT should be NORBT for MP2 !!! TODO /aug2013
!        END IF
      ELSE IF (HSROHF .OR. NASHT .LE. 1) THEN
         CALL DZERO(WRK(KOCCUP),NORBT)
         DO, I = 1, NISHT
            IX = ISX(I)
            WRK(KOCCUP-1+IX) = 2.0D0
         END DO
         IF (HSROHF) THEN
            DO, I = NISHT+1, NOCCT
               IX = ISX(I)
               WRK(KOCCUP-1+IX) = 1.0D0
            END DO
         ELSE IF (NASHT .EQ. 1) THEN
            IX = ISX(NISHT+1)
            WRK(KOCCUP-1+IX) = NACTEL
         END IF
      ELSE
         NWARN = NWARN + 1
         WRITE (LUPRI,'(//A//A/A/A/)')
     &      ' *** Label "NATOCC  " not found on SIRIUS.RST in SIRPOP',
     &      ' WARNING, occupation of active orbitals not available',
     &      '          population analysis aborted.',
     &      '          Use ".TRACI" option.'
         GO TO 9999
      END IF
      IF (MCTYPE .EQ. 2) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A/A/A)')
     &       ' WARNING, population analysis for RAS wave functions'//
     &       ' uses pseudo-natural orbitals,',
     &       '          that is, the off-diagonal elements of the'//
     &       ' one-electron density matrix which',
     &       '          couple the different RAS spaces are ignored.'
      END IF
C
C     *****     Perform energy analysis
C
      CALL VIRIAL (WRK(KCMO),WRK(KOCCUP),WRK(KOVLP))
C
C     *****     Transform from symmetry orbital basis to cgto basis
C
      CALL POPTRA (WRK(KCMO),WRK(KOVLP),WRK(KCMOCG),WRK(KSS))
      IF (IPRMUL .GT. 20) THEN
         WRITE (LUPRI,'(/A)') 'SIRPOP: occuppied MOs in AO basis:'
         CALL OUTPUT(WRK(KCMOCG),1,NBAST,1,NOCCT,NBAST,NOCCT,-1,LUPRI)
         WRITE (LUPRI,'(/A)') 'SIRPOP: Overlap matrix in SO basis:'
         CALL OUTPUT(WRK(KOVLP),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         WRITE (LUPRI,'(/A)') 'SIRPOP: Overlap matrix in AO basis:'
         CALL OUTPUT(WRK(KSS),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
      call memchk('after virial',wrk,1)
C
C     *****     Perform Mulliken-type population analysis
C
      IF (MULLKM) THEN
         CALL MULLIKEN (WRK(KOCCUP), WRK(KPOPNAB),WRK(KPOPABI),
     *                  WRK(KPOPCG),WRK(KPOPC),WRK(KCMOCG),WRK(KSS),
     *                  WRK(KPOPNET))
      END IF
C
C     *****     Dipole moments
C
      IF (DIPOLM) THEN
         write(lupri,'(/A)')
     &     'INFO: dipole moment not implemented in this version'
c        KDIPM  = KFREE
c        KWRK3D = KDIPM  + NNBAST*3
c        CALL POP_DIPMOM (NNBAST,WRK(KCMO),WRK(KOCCUP),WRK(KDIPM))
      END IF
C
C     *****     Quadrupole moments  *****
C
      IF (QUADRM) THEN
         write(lupri,'(/A)')
     &     'INFO: quadrupole moment not implemented in this version'
c        KQUMAT = KFREE
c        KWRK3Q = KQUMAT + NNBAST*6
c        CALL QUADRP (NNBAST,WRK(KCMO),WRK(KOCCUP),WRK(KQUMAT))
      END IF
C
 9999 CONTINUE
      CALL MEMREL('SIRPOP',WRK,1,1,KFREE,LFREE)
      CALL QEXIT('SIRPOP')
      IPRMUL  = IPRMUL_SAVE
      IPOPTYP = IPOPTYP_SAVE
      RETURN
      END
C  /* Deck virial */
      SUBROUTINE VIRIAL (CMO,OCCUP,T)
C
C Purpose:
C      Form the kinetic energy and print an energy analysis
C
#include "implicit.h"
#include "dummy.h"
C
      DIMENSION CMO(NCMOT),OCCUP(NORBT),T(NNBAST)
C
C   INFINP : POTNUC,
C   INFIND : ISX()
C   INFORB : NBAST,
C   INFOPT : EMCSCF,
C
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "priunit.h"
#include "infpop.h"
#include "infinp.h"
#include "infind.h"
#include "inforb.h"
#include "infopt.h"
#include "infpri.h"
C
      LOGICAL   FOUND
C
C     **********     Kinetic energy
C     Read kinetic energy matrix from AOONEINT (unit LUONEL)
C
      FOUND = .FALSE.
      CALL RDONEL('KINETINT',FOUND,T,NNBAST)
      IF ( .NOT. FOUND) THEN
         WRITE (LUPRI,'(//A/A/)')
     *     ' *** Label "KINETINT" not found on AOONEINT,',
     *     '     basic information for VIRIAL routine is missing.'
         RETURN
      END IF
C
C
C     **********     Compute kinetic energy     **********
C
C              TK = Sum N  Sum C  * C  * <j/T/k>
C                    p   p j,k  pj   pk
C
C     where p goes over occupied mo's and j and k over symmetry
C     orbitals. note that only occupied orbitals are included in
C     the matrix cmo, so no skip over secondary mo's is needed.
C
C     MATBLK indicates where the current symmetry block in T matrix
C            starts.
C     IST    tells where the appropriate line in the T matrix starts.
C
      TK = 0.0D0
      IP = 0
      MATBLK = 0
C
      DO 50 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         NOCCI = NOCC(ISYM)
         IPSTRT = ICMO(ISYM)
C
         DO 40 ICC = 1,NOCCI
            IP  = IP + 1
            IST = MATBLK
C
            DO 30 J = 1,NBASI
               JORB = IPSTRT + J
C
               IPX = ISX(IP)
               FAC = 2.0D0*CMO(JORB)*OCCUP(IPX)
               DO 20 K = 1,J-1
                  TK = TK + FAC*CMO(IPSTRT+K)*T(IST+K)
   20          CONTINUE
C
C              Take care of the case j = k
C
               TK = TK + CMO(JORB)*CMO(JORB)*OCCUP(IPX)*T(IST+J)
               IST = IST + J
   30       CONTINUE
C
            IPSTRT = IPSTRT + NBASI
   40    CONTINUE
C
         MATBLK = MATBLK + NBASI*(NBASI + 1)/2
   50 CONTINUE
C
C     **********     Energy analysis and virial theorem     **********
C
      V   = EMCSCF - TK
      VIR = -V/TK
      IF (IPRMUL .GE. 0 .OR. IPRMUL .LT. -30) THEN
         WRITE(LUPRI,1001) EMCSCF,TK,V,POTNUC,VIR
      ELSE IF (IPRMUL .NE. 0) THEN
         WRITE(LUPRI,'(T7,A,F14.6)') 'Virial theorem: -V/T =',VIR
      END IF
 1001 FORMAT(/ 8X,'               Energy analysis'
     &       / 8X,58('-')
     &       /11X,'Total energy',T40,F18.10
     &       /11X,'Kinetic energy',T40,F18.10
     &       /11X,'Potential energy',T40,F18.10
     &       /11X,'Nuclear repulsion energy',T40,F18.10
     &       /11X,'-V/T',T40,F14.6
     &       / 8X,58('-'))
      RETURN
      END
C  /* Deck poptra */
      SUBROUTINE POPTRA (CMO,S,CMOCGT,SS)
C
C Purpose: transform the mo coefficients and overlap matrix
C          from symmetry orbital basis to cgto basis.
C          Called from SIRPOP
C
#include "implicit.h"
C
      DIMENSION CMO(NCMOT)         , S(NBAST,NBAST)
      DIMENSION CMOCGT(NBAST,NOCCT), SS(NBAST,NBAST)
C
C   INFORB ? NBAST,NNBAST,
C
#include "priunit.h"
#include "mxcent.h"
#include "inforb.h"
#include "infpop.h"
#include "infpri.h"
! aosotr: CTRAN, ITRAN, JTRAN
#include "maxorb.h"
#include "aosotr.h"
C
      LOGICAL FOUND
C
C     Read AO integrals for population analysis
C
      FOUND = .FALSE.
      IF (IPOPTYP .EQ. 1) THEN
C     *****   Read in overlap integrals from AOONEINT (unit LUONEL)
         CALL RDONEL('OVERLAP ',FOUND,S,NNBAST)
         CALL PKSYM1(SS,S,NBAS,NSYM,-1)
         CALL DSPTGE(NBAST,SS,S)
         IF ( .NOT. FOUND) THEN
            WRITE (LUPRI,'(//A/A/)')
     *     ' *** Label "OVERLAP" not found on AOONEINT,',
     *     '     basic information for POPTRA is missing.'
            RETURN
         END IF
         POPANA_TYPE = '   MULPOP   '
      ELSE IF (IPOPTYP .EQ. 2) THEN
C     *****   Read in overlap integrals from AOPROPER (unit LUPROP)
         CALL RDPROP('OVERLAP ',SS,FOUND)
         CALL DSPTGE(NBAST,SS,S)
         POPANA_TYPE = '   MULPOP x '
      ELSE IF (IPOPTYP .EQ. 3) THEN
         CALL RDPROP('HJPOPOVL',S,FOUND)
         POPANA_TYPE = '   HJAAJPOP '
      ELSE
         call quit('POPTRA: unknown IPOPTYP')
      END IF
C

      if (iprmul .ge. 25) then
         write (lupri,*) 'Population overlap matrix type ',
     &      popana_type, ipoptyp
         call output(S,1,nbast,1,nbast,nbast,nbast,-1,lupri)
      end if
C
C
C     *****     TRANSFORM THE MO COEFFICIENTS
C
C     for each mo you have nbas(isym) symmetry coefficients (cmo).
C     for each cmo you have a number of cgto coefficients (cmocgt)
C     so the transformation reads
C
C               C'  = C  * D
C                pa    pk   ka
C
C     where p indicates the mo and k the symmetry orbital. a indexes
C     the cgto's. For each k, all values of a as found in ITRAN are
C     examined.
C
      KSTART = 0
      IPSTC  = 0
C
      CALL DZERO(CMOCGT,NOCCT*NBAST)

      DO 60 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         NOCCI = NOCC(ISYM)
         IPSTS = ICMO(ISYM)
C
         IOCCI = IOCC(ISYM)
         DO 50 ICC = IOCCI+1,IOCCI+NOCCI
C
            DO 40 K = 1,NBASI
               KK = KSTART + K
               J  = JTRAN(KK)
C
               DO 30 JJ = 1,J
                  IA = ITRAN(KK,JJ)
                  CMOCGT(IA,ICC) =
     1               CMO(IPSTS+K)*CTRAN(KK,JJ)
   30          CONTINUE
   40       CONTINUE
            IPSTS = IPSTS + NBASI
   50    CONTINUE
         KSTART = KSTART + NBASI
   60 CONTINUE
C
C     *****     Transform the integrals
C
C     Here we have
C
C             LOOP <a/b> = SUM D  * <k/l> * D
C              a,b         k,l  ka           lb
C
C     The order of the loops is k,a,l,b.
C
      CALL DZERO(SS,N2BASX)
C
      KLSTRT = 0
C
      DO 120 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
C
         DO 110 K=1,NBASI
            KK   = KLSTRT + K
            J    = JTRAN(KK)
            FACKA= J
            FACKA= 1.0D0/FACKA
C
            DO 100 JJ = 1,J
               IA  = ITRAN(KK,JJ)
C
               DO 90 L = 1,NBASI
                  LL = KLSTRT + L
                  M  = JTRAN(LL)
                  FACLB = M
                  FACLB = 1.0D0/FACLB
C
                  DO 80 MM = 1,M
                     IB  = ITRAN(LL,MM)
                     FAC = SIGN(FACKA,CTRAN(KK,JJ))
     &                   * SIGN(FACLB,CTRAN(LL,MM))
                     SS(IA,IB) = SS(IA,IB) + S(KK,LL)*FAC
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
         KLSTRT = KLSTRT + NBASI
  120 CONTINUE
C
      RETURN
      END
C  /* Deck mullik */
      SUBROUTINE MULLIKEN (OCCUP,
     &                     POPNAB,POPABI,POPCGT,POPC,CMOCGT,SS,POPNET)
C
C  Purpose:
C         Given the mo coefficients and overlap integrals in cgto basis,
C         to compute and print out a population analysis.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxash.h"
#include "maxorb.h"
C
! aosotr: AOINFO(:,3)
! nuclei: NAMDEP(:), NATOMS, CHARGE(:), NUCDEG(:)
! symmet: FMULT(), ISYMAX()
#include "aosotr.h"
#include "nuclei.h"
#include "symmet.h"

      DIMENSION OCCUP(NORBT),
     &          POPNAB(NATOMS,NATOMS),POPABI(NATOMS,NATOMS),
     &          POPCGT(NBAST),POPC(NBAST,NOCCT),CMOCGT(NBAST,NOCCT),
     &          SS(NBAST,NBAST),POPNET(NBAST,NBAST)
C
C Used from common blocks:
C   codata.h : DEBYE, DIPSI
C   INFPOP ? MAXTYP, ...
C   INFIND : ISX(:), ISMO(:)
C   INFORB ? NBAST,
C   INFTAP : LUTEMP
C   CCOM   : GTOTYP(:)
C   SHELLS : NCENT(:)
C
#include "codata.h"
#include "infpop.h"
#include "infind.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "ccom.h"
#include "shells.h"
C
      REAL*8      POPTSM(MXCENT,2), DIPTOT(3)
      CHARACTER*1 SPDCAR

      IF (IPRMUL.GE.20) THEN
         WRITE(LUPRI,'(//A/)')
     &      ' IAOINFO(1:NBAST,1:3) (shell, center, type)'
         DO IA = 1,NBAST
            WRITE(LUPRI,'(3I20)') IAOINFO(IA,1:3)
         END DO
      END IF
C
C   Initialize the matrices where populations will be accumulated
C
      N2_POPNAB = NATOMS*NATOMS
      CALL DZERO(POPNAB,N2_POPNAB)
      CALL DZERO(POPNET,N2BASX)
      CALL DZERO(POPCGT,NBAST)
      CALL DZERO(POPC,  NBAST*NOCCT )
C
C     *****     Prepare the raw material
C
C     at this stage
C
C         POPNAB will contain the net and overlap populations between
C                nuclei A and B, summed over all mo's
C         POPABI will be the same as popnab, but resolved for each
C                mo separately (these are written out on tape)
C         POPC   will be the same as previous, but resolved for each
C                mo separately. one line of length nbast is reserved
C                for each mo in core
C         POPNET will contain net and overlap populations between
C                cgto basis functions, summed over all MO's
C
      IF (IPRMUL .GE. 5) THEN
         LUTEMP = -99
         CALL GPOPEN(LUTEMP,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUTEMP
      END IF
      DO 80 IP = 1,NOCCT
         IPX = ISX(IP)
         CALL DZERO(POPABI,N2_POPNAB)
C
         DO 75 IA = 1,NBAST
            POPAB_FAC = CMOCGT(IA,IP)*OCCUP(IPX)
            NUCA = IAOINFO(IA,2)
         IF (NUCA .GT. NATOMS) GO TO 75 ! skip auxilliary basis centers
         IF (ABS(POPAB_FAC) .LT. 1.0D-14) GO TO 75
            DO 70 IB = 1,NBAST
C
               NUCB = IAOINFO(IB,2)
            IF (NUCB .GT. NATOMS) GO TO 70 ! skip auxilliary basis centers
C
C              hjaaj Mar 2006: The factor in her1int:POPOVLINT corresponds
C                              to using SS(IB,IA) here.
               POPAB = POPAB_FAC * CMOCGT(IB,IP) * SS(IB,IA)
C
C              Distribute
C
               POPABI(NUCA,NUCB) = POPABI(NUCA,NUCB) + POPAB
               POPNAB(NUCA,NUCB) = POPNAB(NUCA,NUCB) + POPAB
               POPNET(IA,IB)     = POPNET(IA,IB)     + POPAB
C
C              (gross populations follow)
C
               POPC(IA,IP) = POPC(IA,IP) + POPAB
   70       CONTINUE ! DO 70 IB = 1,NBAST
   75    CONTINUE ! DO 75 IA = 1,NBAST
C
#ifdef SIRPOP_DEBUG
         write (lupri,*) 'popabi for orbital ',IP
         call output (popabi, 1,NATOMS,1,NATOMS,NATOMS,NATOMS,-1,lupri)
         write (lupri,*) 'popabn after orbital ',IP
         call output (popnab, 1,NATOMS,1,NATOMS,NATOMS,NATOMS,-1,lupri)
#endif
         IF (IPRMUL .GE. 5) WRITE (LUTEMP) POPABI
C
   80 CONTINUE ! DO 80 IP = 1,NOCCT
C
C     *****     compute some refined results
C
C     Now   POPTYP will contain gross populations on the cgto's
C                  of type t on center a, summed over all mo's
C           POPCGT will contain the gross populations on each cgto,
C                  summed over all mo's
C           POPNUC will contain the gross populations on each atom,
C                  summed over all mo's
C
      CALL DZERO(POPTYP, NATOMS*MAXTYP)
C
      MAX_ITYP = 0
      DO IA = 1,NBAST
         NUC  = IAOINFO(IA,2)
         ITYP = IAOINFO(IA,3)
         MAX_ITYP = MAX(MAX_ITYP, ITYP)
         DO IP = 1,NOCCT
            POPTYP(ITYP,NUC) = POPTYP(ITYP,NUC) + POPC(IA,IP)
            POPCGT(IA)       = POPCGT(IA)       + POPC(IA,IP)
         END DO
      END DO
C
      DO NUCA = 1,NATOMS
         POPA_TOT= 0.0D0
         DO NUCB = 1,NATOMS
            POPA_TOT= POPA_TOT+ POPNAB(NUCA,NUCB)
         END DO
         POPNUC(NUCA) = POPA_TOT
      END DO
C
C     *****     Now just write everything out in neat order
C
C
C     Populations between centers
C
      if (iprmul .ge. 1) then
      WRITE (LUPRI,1004) POPANA_TYPE
 1004 FORMAT(//6X,'Total net, overlap, and gross populations,',
     &            ' summed over all MOs',/6X,'Type:',A)
C
      IEND = 0
    5    IST  = IEND + 1
         IEND = MIN(NATOMS, IEND + 10)
         IF (IEND .EQ. NATOMS) THEN
            WRITE (LUPRI,1009) (NAMDEP(NUC),NUC=IST,IEND),'Gross '
         ELSE
            WRITE (LUPRI,1009) (NAMDEP(NUC),NUC=IST,IEND)
         END IF
         WRITE (LUPRI,'()')
         DO 140 NUCA = 1,NATOMS
            POPNAB_MAX = MAXVAL(POPNAB(NUCA,IST:IEND))
            POPNAB_MIN = MINVAL(POPNAB(NUCA,IST:IEND))
         IF (POPNAB_MAX .LT. 1.0D-4 .AND. POPNAB_MIN .GT. -1.0D-4)
     &      GO TO 140
            IF (IEND .EQ. NATOMS) THEN
               WRITE (LUPRI,1012)
     &         NAMDEP(NUCA),(POPNAB(NUCA,NUCB),NUCB=IST,IEND),
     &         POPNUC(NUCA)
            ELSE
               WRITE (LUPRI,1012)
     &         NAMDEP(NUCA),(POPNAB(NUCA,NUCB),NUCB=IST,IEND)
            END IF
  140    CONTINUE
      IF (IEND.LT.NATOMS)   GO TO 5
      end if
C
      IF ( IPRMUL .GE. 5 ) THEN
         WRITE (LUPRI,1001) POPANA_TYPE
 1001    FORMAT(//6X,'Net, overlap, and gross populations, ' 
     &              ,' contributions from individual MOs'
     &           /6X,'Type:',A)
         REWIND LUTEMP
         ISYM = 1
         DO 130 IP = 1,NOCCT
            IPX = ISX(IP)
            IP_SYM = ISMO( IPX )
            WRITE (LUPRI,1002) IP, OCCUP(IPX), IP_SYM
 1002       FORMAT (/6X,'Molecular orbital number',I4,
     &      ', occupation',F8.5,', symmetry',I3)
C
            READ (LUTEMP) POPABI
            IEND = 0
    3          IST  = IEND + 1
               IEND = MIN(NATOMS, IEND + 10)
               IF (IEND .EQ. NATOMS) THEN
                  WRITE (LUPRI,1009) (NAMDEP(NUC),NUC=IST,IEND),'Gross '
               ELSE
                  WRITE (LUPRI,1009) (NAMDEP(NUC),NUC=IST,IEND)
               END IF
               WRITE (LUPRI,'()')
               DO NUCA = 1,NATOMS
                  POPABI_MAX = MAXVAL(POPABI(NUCA,IST:IEND))
                  POPABI_MIN = MINVAL(POPABI(NUCA,IST:IEND))
               IF (POPABI_MAX .LT. 1.0D-4 .AND. POPABI_MIN .GT. -1.0D-4)
     &            CYCLE
                  IF (IEND .EQ. NATOMS) THEN
                     PAI = DSUM(NATOMS,POPABI(NUCA,1),NATOMS)
                     WRITE (LUPRI,1012)
     &               NAMDEP(NUCA),(POPABI(NUCA,NUCB),NUCB=IST,IEND),PAI
                  ELSE
                     WRITE (LUPRI,1012)
     &               NAMDEP(NUCA),(POPABI(NUCA,NUCB),NUCB=IST,IEND)
                  END IF
               END DO
            IF (IEND.LT.NATOMS)   GO TO 3
  130    CONTINUE
         CALL GPCLOSE(LUTEMP,'DELETE')
      END IF
C
C     Gross populations of cgto's
C
      IF (IPRMUL.GE.10) THEN
         WRITE (LUPRI,1005) POPANA_TYPE
 1005    FORMAT (//6X,'Partial orbital gross populations:',
     &   /6X,'Gross populations of AO basis functions for each MO',
     &   /6X,'Type:',A)
         DO 160 IP = 1,NOCCT
            IPX = ISX(IP)
            IP_SYM = ISMO( IPX )
            WRITE (LUPRI,1002) IP, OCCUP(IPX), IP_SYM
            WRITE (LUPRI,1006) (IA,POPC(IA,IP),IA=1,NBAST)
  160    CONTINUE
 1006    FORMAT (10(I4,F8.4))
      END IF
C
      if (iprmul .ge. 4) then
         WRITE (LUPRI,1007) POPANA_TYPE
         !old! WRITE (LUPRI,1006) (IA,POPCGT(IA),IA = 1,NBAST)
         N_ZERO = 0
         DO IA = 1, NBAST
            NUC  = IAOINFO(IA,2)
            ITYP = IAOINFO(IA,3)
            IF (ABS(POPCGT(IA)) .LT. 1.0D-5) THEN
               N_ZERO = N_ZERO + 1
               WRITE(LUPRI,'(I5,2X,A,2X,2A)')
     &         IA, NAMDEP(NUC), GTOTYP(ITYP),' Less than 1.0e-6'
            ELSE
               WRITE(LUPRI,'(I5,2X,A,2X,A,F10.4)')
     &         IA, NAMDEP(NUC), GTOTYP(ITYP),POPCGT(IA)
            END IF
         END DO
         POPTOT = DSUM(NBAST,POPCGT,1)
         WRITE (LUPRI,'(/A,F10.4/A,I5,A,I5)')
     &   ' Sum = total electronic charge  :',POPTOT,
     &   ' No. of AOs with gross population < 1.0e-5:',N_ZERO,
     &   ' out of',NBAST
      end if
 1007 FORMAT(//6X,'Total gross populations of AO basis functions,',
     &            ' summed over all MOs'/6X,'Type:',A)
C
C     GROSS POPULATIONS OF ORBITAL TYPES
C
      IF (IPRMUL.GE.5) THEN
         WRITE (LUPRI,1008) POPANA_TYPE
 1008    FORMAT(//6X,'Partial gross population of AO basis',
     &              ' functions for each MO',
     &           /6X,'distributed over AO types for each center',
     &           /6X,'Type:',A)
         DO 200 IP = 1,NOCCT
            IPX = ISX(IP)
            IP_SYM = ISMO( IPX )
            WRITE (LUPRI,1002) IP,OCCUP(IPX),IP_SYM
C
            CALL DZERO(POPTY, MAXTYP*NATOMS)
C
            NTYP = 0
            DO IA = 1,NBAST
               NUC  = IAOINFO(IA,2)
               ITYP = IAOINFO(IA,3)
               NTYP = MAX(NTYP,ITYP)
               POPTY(ITYP,NUC) = POPTY(ITYP,NUC) + POPC(IA,IP)
            END DO
C
            IEND = 0
    7       IST  = IEND + 1
            IEND = MIN(NATOMS,IEND + 10)
            WRITE (LUPRI,1009) (NAMDEP(NUC),NUC = IST,IEND)
            WRITE (LUPRI,'()')
            DO 190 ITYP = 1,NTYP
               POPTY_MAX = MAXVAL(POPTY(ITYP,IST:IEND))
               POPTY_MIN = MINVAL(POPTY(ITYP,IST:IEND))
            IF (POPTY_MAX .LT. 1.0D-4 .AND. POPTY_MIN .GT. -1.0D-4)
     &         GO TO 190
               WRITE (LUPRI,1010)
     &         GTOTYP(ITYP),(POPTY(ITYP,NUC), NUC=IST,IEND)
  190       CONTINUE
            IF (IEND.LT.NATOMS)   GO TO 7
  200    CONTINUE
      END IF

 1009 FORMAT (/11X,A6,(11(4X,A6)))
 1010 FORMAT (3X,A4,10F10.4)
 1012 FORMAT (1X,A6,10F10.4)
C
      if (iprmul .gt. 0) then
      WRITE (LUPRI,1013) POPANA_TYPE
 1013 FORMAT(/6X,'Total gross population of AO basis',
     &           ' functions, summed over all MOs',
     &       /6X,'distributed over AO types for each center',
     &       /6X,'Type:',A)

      IEND = 0
    9 IST  = IEND + 1
      IEND = MIN(NATOMS, IEND + 10)
      WRITE (LUPRI,1009) (NAMDEP(NUC),NUC = IST,IEND)
      WRITE (LUPRI,'()')
      LQ = 0
      POPTSM(IST:IEND,1:2) = 0.0D0
      DO 210 ITYP = 1,MAX_ITYP
         DO NUC = IST,IEND
            POPTSM(NUC,2) = POPTSM(NUC,2) + POPTYP(ITYP,NUC)
         END DO
         IF (LQ .GT. 0) THEN
            POPTYP_MAX = MAXVAL(POPTYP(ITYP,IST:IEND))
            POPTYP_MIN = MINVAL(POPTYP(ITYP,IST:IEND))
            IF (POPTYP_MAX .GE. 1.0D-4 .OR. POPTYP_MIN .LE. -1.0D-4)
     &         WRITE (LUPRI,1010)
     &         GTOTYP(ITYP),(POPTYP(ITYP,NUC),NUC=IST,IEND)
         END IF

      IF (GTOTYP(ITYP+1)(1:1) .EQ. GTOTYP(ITYP)(1:1)) GO TO 210
C        for cartesian GTOTYP(*)(1:1) will be 's', 'p', etc.
C        for spherical GTOTYP(*)(1:1) will be '1', '2', etc.
C        thus: when char 1 changes so does l quantum number and
C        using GTOTYP to decide works both for cartesian and spherical

         WRITE (LUPRI,1010) SPDCAR(LQ),(POPTSM(NUC,2), NUC=IST,IEND)
         WRITE (LUPRI,1099)
         POPTSM(IST:IEND,1) = POPTSM(IST:IEND,1) + POPTSM(IST:IEND,2)
         POPTSM(IST:IEND,2) = 0.0D0
         LQ = LQ + 1
  210 CONTINUE
      WRITE (LUPRI,1010) 'All',(POPTSM(NUC,1), NUC=IST,IEND)
 1099 FORMAT (1X)
      IF (IEND.LT.NATOMS)   GO TO 9
      end if
C
C     POPULATION at EACH center
C
      IF (IPRMUL.GE.6) THEN
        WRITE (LUPRI,1014) POPANA_TYPE
 1014 FORMAT (/6X,'Partial gross population of the centers',
     *            ' from each MO'/6X,'Type:',A)
C
        DO 240 IP = 1,NOCCT
          POPABI(1:NATOMS,1) = 0.0D0
          DO IA = 1,NBAST
            NUC = IAOINFO(IA,2)
            POPABI(NUC,1) = POPABI(NUC,1) + POPC(IA,IP)
          ENDDO
C
          IPX = ISX(IP)
          IP_SYM = ISMO( IPX )
          WRITE (LUPRI,1002) IP, OCCUP(IPX), IP_SYM
          WRITE (LUPRI,1015) (NAMDEP(NUC),POPABI(NUC,1),NUC=1,NATOMS)
c1015     FORMAT (4(1X,A6,F10.4))
 1015     FORMAT (5(1X,A6,F8.2,' ;'))
  240 CONTINUE
      END IF
C
      POPTOT_E = DSUM(NATOMS,POPNUC,1)
C
C     Add effective nuclear charges (nuclear charge minus ECP charge)
C     Calculate dipole moment based on atomic charges
C
      DIPTOT(1:3) = 0.0D0
      NUC = 0
      DO IATOM = 1, NUCIND

         Q_EFF = CHARGE(IATOM)
         N_DEG = NUCDEG(IATOM)
         DO I = 1, N_DEG
            NUC = NUC + 1
            POPNUC(NUC) = Q_EFF - POPNUC(NUC)
         END DO

         FAC = FMULT(ISTBNU(IATOM))*POPNUC(NUC)
         DO ICOOR = 1,3
            IF (ISYMAX(ICOOR,1) .EQ. 0) THEN
               DIPTOT(ICOOR) = DIPTOT(ICOOR) + FAC*CORD(ICOOR,IATOM)
            END IF
         END DO

      END DO

      POPTOT = DSUM(NATOMS,POPNUC,1)
      DIPMOM = SQRT(DIPTOT(1)**2 + DIPTOT(2)**2 + DIPTOT(3)**2)

 2015 FORMAT ('@ ',A12,10(A6,F6.2,'; ')/,('@',T15,10(A6,F6.2,'; ')))
      IF (iprmul .lt. 0) then
C     ... only print max -iprmul populations in each iteration

         mnuc = min(NATOMS,-iprmul)
         WRITE (LUPRI,2015) POPANA_TYPE,
     &      (NAMDEP(NUC),POPNUC(NUC),NUC=1,mnuc)
         IF (mnuc .lt. NATOMS) WRITE (LUPRI,'(3A,F12.4,A,3F12.4)')
     &      '@ ',POPANA_TYPE,' dipole moment & components:',
     &      DIPMOM, ';  ',DIPTOT(1:3)

      else

         CALL HEADER('Total atomic gross populations',10)
         WRITE (LUPRI,2015) POPANA_TYPE,
     &      (NAMDEP(NUC),POPNUC(NUC),NUC=1,NATOMS)
         WRITE (LUPRI,'(/A,F10.4)')
     &   '@ Sum = total charge:',POPTOT

       IF (IPRMUL .GE. 2) THEN
       ! note: DIPMOM cannot be expected to be accurate, it will only
       ! be accurate if the local dipole moments on each atom are zero.
       ! /hjaaj Mar 2017
         WRITE (LUPRI,'(/A)') ' Total molecular dipole '//
     &   ' moment calculated from gross populations'
         WRITE (LUPRI,'(17X,A,15X,A,10X,A/3X,3F19.6)')
     &             'au','Debye','C m (/(10**-30)',
     &             DIPMOM, DEBYE*DIPMOM, DIPSI*DIPMOM
         WRITE (LUPRI,'(/A)') ' x,y,z components:'
         CALL DP0PRI(DIPTOT)
      END IF

      end if
C
C     NET and overlap populations between cgto basis functions - sum

      IF (IPRMUL.GE.8) THEN
         WRITE (LUPRI,1017) POPANA_TYPE
 1017    FORMAT(//6X,'Total net and overlap populations between the AO',
     &          ' basis functions'
     &          /6X,'Summed over all MOs'/6X,'Type:',A)
C
         IEND = 0
   12    IST  = IEND + 1
         IEND = MIN(NBAST, IEND + 10)
         WRITE (LUPRI,1018) (IB,IB=IST,IEND)
 1018    FORMAT (/8X,'AO',10I9)
         WRITE (LUPRI,'(A)') '     AO'
         DO IA=1,NBAST
            POPNET_MAX = MAXVAL(POPNET(IA,IST:IEND))
            POPNET_MIN = MINVAL(POPNET(IA,IST:IEND))
         IF (POPNET_MAX .LT. 1.0D-4 .AND. POPNET_MIN .GT. -1.0D-4)
     &      CYCLE
            WRITE(LUPRI,'(I7,5X,10F9.4)') IA,(POPNET(IA,IB),IB=IST,IEND)
         END DO
         IF (IEND.LT.NBAST)   GO TO 12
      END IF
C
C
      RETURN
      END
C  /* Deck pop_dipmom */
      SUBROUTINE POP_DIPMOM (MATDIM,CMO,OCCUP,DIPMAT)
C
C PURPOSE:
C     Compute the electronic and total dipole moments
C
#include "implicit.h"
C
      DIMENSION CMO(NCMOT),DIPMAT(MATDIM,3),OCCUP(NORBT)
C
C   INFORB ? NBAST,
C
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "priunit.h"
#include "infpop.h"
#include "infind.h"
#include "inforb.h"
#include "infpri.h"
C
      LOGICAL FOUND
C
      WRITE (LUPRI,1001)
 1001 FORMAT (//,22X,'Dipole moment analysis')
C
C     *****     READ IN THE DIPOLE MOMENT INTEGRALS
C
      FOUND = .FALSE.
      CALL RDONEL('DIPOLMOM',FOUND,DIPMAT,MATDIM)
      IF ( .NOT. FOUND) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(//A/A/)')
     *     ' *** INFO: Label "DIPOLMOM" not found on AOONEINT,',
     *     '     basic information for POP_DIPMOM is missing.'
         RETURN
      END IF
C
C     *****     Compute the electronic contribution to dipole moment
C
      WRITE (LUPRI,1002)
 1002 FORMAT(//25X,'Orbital dipole moments'//28X,'Vector component'
     *       //10X,'Orbital',12X,'X',11X,'Y',11X,'Z',5X,'Occ.number'/)
      IP = 0
      IPSTRT = 0
      MATBLK = 0
C
      DO 50 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         NOCCI = NISH(ISYM) + NASH(ISYM)
C
         DO 40 ICC = 1,NOCCI
            IP = IP + 1
            IPX = ISX(IP)
            IST = MATBLK
            DX = 0.0D0
            DY = 0.0D0
            DZ = 0.0D0
C
            DO 30 J=1,NBASI
               JORB = IPSTRT + J
C
               DO 20 K = 1,(J-1)
                  DENS = 2.0D0*CMO(JORB)*CMO(IPSTRT+K)*OCCUP(IPX)
                  DX = DX + DENS*DIPMAT(IST+K,1)
                  DY = DY + DENS*DIPMAT(IST+K,2)
                  DZ = DZ + DENS*DIPMAT(IST+K,3)
   20          CONTINUE
C
C              Case j = k, divide by two
C              because no (j,k) permutation exist
C
               DENS = CMO(JORB)*CMO(JORB)*OCCUP(IPX)
               DX = DX + DENS*DIPMAT(IST+J,1)
               DY = DY + DENS*DIPMAT(IST+J,2)
               DZ = DZ + DENS*DIPMAT(IST+J,3)
               IST= IST + J
   30       CONTINUE
C
            IPSTRT = IPSTRT + NBASI
C
            IF (IFXYZ(1).EQ.1)   DX = 0.0D0
            IF (IFXYZ(2).EQ.1)   DY = 0.0D0
            IF (IFXYZ(3).EQ.1)   DZ = 0.0D0
            DIPOL(1) = DIPOL(1) - DX
            DIPOL(2) = DIPOL(2) - DY
            DIPOL(3) = DIPOL(3) - DZ
            IF (IPRMUL .GE. 3) THEN
               WRITE (LUPRI,1003) IP,DX,DY,DZ,OCCUP(IPX)
 1003          FORMAT(10X,I4,7X,4F12.6)
            END IF
C
   40    CONTINUE
C
         MATBLK = MATBLK + NBASI*(NBASI + 1)/2
   50 CONTINUE
C
C     To obtain exactly zero results in cases where symmetry
C     dictates that the dipole moment must vanish, explicitly
C     set the components to zero
C
      IF (IFXYZ(1).EQ.1)   DIPOL(1) = 0.0D0
      IF (IFXYZ(2).EQ.1)   DIPOL(2) = 0.0D0
      IF (IFXYZ(3).EQ.1)   DIPOL(3) = 0.0D0
C
C     *****     Print out the total dipole moments
C
C
      WRITE(LUPRI,1005) (DIPOL(I),I=1,3)
 1005 FORMAT(/11X,'TOTAL',5X,3F12.6)
C
C     MAGNITUDE
C
      DIPITO = DIPOL(1)**2 + DIPOL(2)**2 + DIPOL(3)**2
      DIPITO = SQRT(DIPITO)
      DIPITD = DIPITO * 2.5415D0
      WRITE (LUPRI,1006) DIPITO,DIPITD
 1006 FORMAT (/6X,'TOTAL DIPOLE MOMENT (A.U.) = ',F10.6,
     *       //6X,'TOTAL DIPOLE MOMENT (D)    = ',F10.6)
C
C
      RETURN
      END
C  /* Deck quadrp */
      SUBROUTINE QUADRP (MATDIM,CMO,OCCUP,QUMAT)
C
C  Purpose: Compute the quadrupole moments
C
C based on release 82 09 01
C 900709-hjaaj: corrected error in diagonalization of quadrupole tensor
C    relative to mass center (off-diagonal elements were from tensor
C    relative to origo).  Print <r2>.
C
#include "implicit.h"
C
      DIMENSION CMO(NCMOT),QUMAT(MATDIM,6),OCCUP(NORBT)
C
C   INFIND : ISX()
C   INFORB ? NBAST,
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infind.h"
#include "inforb.h"
#include "infpop.h"
#include "infpri.h"
C
      DIMENSION QUPOL(6),A(6),EIVR(9),ELQPOL(6),QOPOL(6)
C
      LOGICAL FOUND
C
      WRITE (LUPRI,1001)
 1001 FORMAT (///,21X,' QUADRUPOLE MOMENT ANALYSIS')
C
C     *****     Read in the quadrupole moment integrals
C
      FOUND = .FALSE.
      CALL RDONEL('QUADRP  ',FOUND,QUMAT,MATDIM)
      IF ( .NOT. FOUND) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(//A/A/)')
     *     ' *** INFO: Label "QUADRP" not found on AOONEINT,',
     *     '     basic information for QUADRP  is missing.'
         RETURN
      END IF
C
C     *****     Compute the electronic contribution to quadrupole moment
C
      IF (IPRMUL .GE. 3) THEN
         WRITE(LUPRI,1002)
      ELSE
         WRITE(LUPRI,2002)
      END IF
 1002 FORMAT(//23X,'Occupation weighted orbital quadrupole moments',
     *        /23X,'(except <r2> which shows spatial extent of orbital)'
     *       //' Orbital',5X,'XX',7X,'YY',7X,'ZZ',7X,'XY',7X,'XZ',
     *         7X,'YZ',7X,'RR     <r2>',/)
 2002 FORMAT(//23X,'Occupation weighted orbital quadrupole moments',
     *       //' Orbital',5X,'XX',7X,'YY',7X,'ZZ',7X,'XY',7X,'XZ',
     *         7X,'YZ',7X,'RR',/)
C
      IP = 0
      IPSTRT = 0
      MATBLK = 0
      DO 20 ICOMP = 1,6
         ELQPOL(ICOMP) = 0.0D0
   20 CONTINUE
C
      DO 90 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         NOCCI = NISH(ISYM) + NASH(ISYM)
C
         DO 80 ICC = 1,NOCCI
            IP = IP + 1
            IPX = ISX(IP)
            DO 30 ICOMP = 1,6
               QUPOL(ICOMP) = 0.0D0
   30       CONTINUE
            IST = MATBLK
C
            DO 60 J=1,NBASI
               JORB = IPSTRT + J
C
               DO 40 K = 1,(J-1)
                  DENS = 2.0D0*CMO(JORB)*CMO(IPSTRT+K)
                  DO 40 ICOMP = 1,6
                     QUPOL(ICOMP) = QUPOL(ICOMP)
     1                            + DENS*QUMAT(IST+K,ICOMP)
   40          CONTINUE
C
C              CASE J = K, DIVIDE BY TWO BECAUSE
C              NO (J,K) PERMUTATION EXIST
C
               DENS = CMO(JORB)*CMO(JORB)
               DO 50 ICOMP = 1,6
                  QUPOL(ICOMP) = QUPOL(ICOMP)
     1                         + DENS*QUMAT(IST+J,ICOMP)
   50          CONTINUE
               IST= IST + J
   60       CONTINUE
C
            IPSTRT = IPSTRT + NBASI
C           R2 = < r**2 > for this orbital
            R2 = QUPOL(1) + QUPOL(2) + QUPOL(3)
C
            DO 70 ICOMP = 1,6
               QUPOL(ICOMP)  = QUPOL(ICOMP)*OCCUP(IPX)
               ELQPOL(ICOMP) = ELQPOL(ICOMP) - QUPOL(ICOMP)
   70       CONTINUE
            RR = QUPOL(1) + QUPOL(2) + QUPOL(3)
            IF (JFXYZ(1).EQ.1)   QUPOL(4) = 0.0D0
            IF (JFXYZ(2).EQ.1)   QUPOL(5) = 0.0D0
            IF (JFXYZ(3).EQ.1)   QUPOL(6) = 0.0D0
            IF (IPRMUL .GE. 2) THEN
               WRITE (LUPRI,1003) IP,(QUPOL(ICOMP),ICOMP=1,6),RR,R2
 1003          FORMAT(I5,3X,7F9.4,F7.2)
            END IF
C
   80    CONTINUE
C
         MATBLK = MATBLK + NBASI*(NBASI + 1)/2
   90 CONTINUE
C
C     *****     Print out the total electronic quadrupole moments
C
      SUM = ELQPOL(1) + ELQPOL(2) + ELQPOL(3)
      IF (JFXYZ(1).EQ.1) ELQPOL(4) = 0.0D0
      IF (JFXYZ(2).EQ.1) ELQPOL(5) = 0.0D0
      IF (JFXYZ(3).EQ.1) ELQPOL(6) = 0.0D0
      WRITE(LUPRI,1004) ELQPOL,SUM
 1004 FORMAT(/' Summed ',7F9.4)
C
C     *****     Print out the nuclear quadrupole moment
C
      WRITE (LUPRI,1005) QPOL
 1005 FORMAT(/' Nuclear',6F9.4)
C
C     *****     Print out the total quadrupole moment
C
      DO 100 I=1,6
100      QPOL(I) = QPOL(I) + ELQPOL(I)
      SUM = QPOL(1) + QPOL(2) + QPOL(3)
      WRITE (LUPRI,1006) QPOL,SUM
 1006 FORMAT(/'  Total ',7F9.4)
C
C     *****     Quadrupole moments relative to mass center
C
      QOPOL(1)=QPOL(1)-2*QQFAC(1)*DIPOL(1)
      QOPOL(2)=QPOL(2)-2*QQFAC(2)*DIPOL(2)
      QOPOL(3)=QPOL(3)-2*QQFAC(3)*DIPOL(3)
      QOPOL(4)=QPOL(4)-  QQFAC(1)*DIPOL(2)-QQFAC(2)*DIPOL(1)
      QOPOL(5)=QPOL(5)-  QQFAC(1)*DIPOL(3)-QQFAC(3)*DIPOL(1)
      QOPOL(6)=QPOL(6)-  QQFAC(2)*DIPOL(3)-QQFAC(3)*DIPOL(2)
      SUMR    =QOPOL(1) + QOPOL(2) + QOPOL(3)
      WRITE(LUPRI,1007) QOPOL,SUMR
 1007 FORMAT(//16X,'Quadrupole moments relative to mass center',
     *       //8X,7F9.4)
C
C     *****     Diagonalize quadrupole moment tensor
C
      ITURN=0
      WRITE(LUPRI,1008)
 1008 FORMAT(//6X,
     *   'Diagonalization of quadrupole tensor relative ORIGO')
      A(1)=0.5D0*(3.0D0*QPOL(1)-SUM)
      A(3)=0.5D0*(3.0D0*QPOL(2)-SUM)
      A(6)=0.5D0*(3.0D0*QPOL(3)-SUM)
      A(2)=QPOL(4)*1.5D0
      A(4)=QPOL(5)*1.5D0
      A(5)=QPOL(6)*1.5D0
C
  145 CALL DUNIT(EIVR,3)
      CALL JACO_THR(A,EIVR,3,3,3,0.0D0)
      A(2)=A(3)
      A(3)=A(6)
      WRITE(LUPRI,140)
  140 FORMAT(//,18X,'Eigenvalues',16X,'Eigenvectors')
      IEND=0
      DO 141 I=1,3
         IST =IEND+1
         IEND=IEND+3
         WRITE(LUPRI,142)A(I),(EIVR(J),J=IST,IEND)
  142    FORMAT(/,20X,F10.6,10X,3F10.6)
  141 CONTINUE
      IF (ITURN.NE.1) THEN
         ITURN = 1
C
         WRITE(LUPRI,1009)
 1009    FORMAT(//'   Diagonalization of quadrupole tensor relative',
     *          ' to mass center')
         A(1)=0.5D0*(3.0D0*QOPOL(1)-SUM)
         A(3)=0.5D0*(3.0D0*QOPOL(2)-SUM)
         A(6)=0.5D0*(3.0D0*QOPOL(3)-SUM)
         A(2)=QOPOL(4)*1.5D0
         A(4)=QOPOL(5)*1.5D0
         A(5)=QOPOL(6)*1.5D0
         GO TO 145
C     ^-----------
      END IF
C
      RETURN
      END
C -- end of sirpop.F --
