!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C  /* Deck sosdrv */
      SUBROUTINE SOSDRV(WORK,LWORK,PASS)
C
C      Stephan P. A. Sauer 13/11-1997
C
C      This routine is the driver for the calculation of contributions
C      to the indirect nuclear spin-spin coupling constants from 
C      individual molecular orbitals.
C
#include "implicit.h"
C
C MXCOOR used from include file mxcent.h
C LUPRI used from include file priunit.h
#include "mxcent.h"
#include "priunit.h"
C
C NUCDEP used from COMMON /NUCLEI/
C NSYM, NORBT, NOCCT, NVIRT, NCMOT used from COMMON /INFORB/
C ISPPRI used from COMMON /SPNOUT/
C NVARPT used from COMMON /INFLIN/
#include "nuclei.h"
#include "inforb.h"
#include "spnout.h"
#include "inflin.h"
C
#include "maxorb.h"
#include "locinf.h"
      LOGICAL PASS
      DIMENSION WORK(LWORK)
      CALL QENTER('SOSDRV')
C
C--------------------
C     Set print level
C--------------------
C
      IPRSOS = ISPPRI
C
C     IF (SKIP) RETURN
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRSOS .GT. 0) WRITE (LUPRI,'(A,/)')
     &     '1 ---------- Output from SOSDRV ----------'
C
C--------------------------
C     Check for no symmetry
C--------------------------
C
      IF (NSYM .NE. 1) THEN
         WRITE(LUPRI,'(/2(A))') 'SOSDRV : Sum Over States calculation',
     &                          ' only possible without symmetry'
         CALL QUIT('SOSDRV : SOS only possible without symmetry')
      END IF
C
C-------------------------------
C     Initialize some dimensions
C-------------------------------
C
      NLOCC  = NOCCT
      NLVIR  = NVIRT
      NLVIRS = NLVIR
C
c      SINGNUL = .TRUE.
C
      IF (SOSOCC .OR. PIPEKL) THEN
         NLVIRS = 1
         NLOCCS = NLOCC
      ELSE
         NLOCCS = 1
      END IF
C
      IF (FBOVIR) THEN
         NLVIR  = NVSET*NOCCT
         NLVIRS = NLVIR
      ELSE IF (FBSETV) THEN
         NLVIR  = NVSET*NOCCT
         NLVIRS = NLVIR
      ELSE IF (FBSTVO) THEN
         NLVIR  = NVSET*NV2LOC
         NLVIRS = NLVIR
      END IF
C
      NLVAR2 = NLOCCS*NLOCCS*NLVIRS*NLVIRS
C
C------------------------
C     1. Work allocation:
C------------------------
C
      KCLOCO = 1
      KCLOCV = KCLOCO + NOCCT*NLOCC
      KWORK1 = KCLOCV + NVIRT*NLVIR
      LWORK1 = LWORK  - KWORK1
C
      IF (LWORK1 .LT. 0) THEN
         WRITE(LUPRI,*)'SOSDRV.1: Need :',KWORK1,'  Available :',LWORK1
         CALL QUIT('Insufficient memory in SOSDRV.1')
      ENDIF
C
C------------------------
C     2. Work allocation:
C------------------------
C
      KPRPMO = KWORK1
      KWORK2 = KPRPMO + NORBT * NORBT * 3
      LWORK2 = LWORK  - KWORK2
C
      IF (LWORK2 .LT. 0) THEN
         WRITE(LUPRI,*)'SOSDRV.2: Need :',KWORK2,'  Available :',LWORK2
         CALL QUIT('Insufficient memory in SOSDRV.2')
      ENDIF
C
C---------------------------------------------------------
C     Generate transformation matrix to localized orbitals
C---------------------------------------------------------
C
      IF (SOSOCC) THEN
C
Cspas suggestio: put IF (FOSTERBOY) THEN
C
        CALL GETLOC(IPRSOS,NLOCC,NLVIR,WORK(KCLOCO),WORK(KCLOCV),
     &              WORK(KPRPMO),WORK(KWORK2),LWORK2)
C
c        CALL GETLOC_FB(IPRSOS,NLOCC,NLVIR,WORK(KCLOCO),WORK(KCLOCV),
c     &              WORK(KPRPMO),WORK(KWORK2),LWORK2)
C
      ELSE IF (PIPEKL) THEN
C
c      CALL GETLOC_PM(IPRSOS,NLOCC,NLVIR,WORK(KCLOCO),WORK(KCLOCV),
c     &            WORK(KPRPMO),WORK(KWORK2),LWORK2)
C
      END IF
C
C------------------------
C     3. Work allocation:
C------------------------
C
      KSPDSO = KWORK2
      KWORK3 = KSPDSO + 3*NUCDEP*3*NUCDEP*NLOCC*NLOCC
      LWORK3 = LWORK  - KWORK3
C
      IF (LWORK3 .LT. 0) THEN
         WRITE(LUPRI,*)'SOSDRV.3: Need :',KWORK3,'  Available :',LWORK3
         CALL QUIT('Insufficient memory in SOSDRV.3')
      ENDIF
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
      CALL DZERO(WORK(KSPDSO),9*NUCDEP*NUCDEP*NLOCC*NLOCC)
C
C=======================
C     Calculate DSO term
C=======================
C
      CALL SOODSO(IPRSOS,NLOCC,WORK(KCLOCO),WORK(KSPDSO),
     &            WORK(KWORK3),LWORK3)
C
      NEXCIS = NVARPT
      NEXCIT = NVARPT
      IF (NSTATS .NE. 0) NEXCIS = NSTATS
      IF (NSTATT .NE. 0) NEXCIT = NSTATT
C
C--------------------------------
C     4. Triplet Work allocation:
C--------------------------------
C
      KEXCTR = KWORK3
      KSD    = KEXCTR + NEXCIT
      KFC    = KSD    + 9*9*NUCDEP*NUCDEP*NLVAR2
      KSDF   = KFC    + NUCDEP*NUCDEP*NLVAR2
      KSPSD  = KSDF   + 9*NUCDEP*NUCDEP*NLVAR2
      KSPFC  = KSPSD  + 9*NUCDEP*NUCDEP*NLVAR2
      KSPSDF = KSPFC  + 9*NUCDEP*NUCDEP*NLVAR2
      KWORK4 = KSPSDF + 9*NUCDEP*NUCDEP*NLVAR2
      LWORK4 = LWORK  - KWORK4
C
      IF (LWORK4 .LT. 0) THEN
         WRITE(LUPRI,*)'SOSDRV.4: Need :',KWORK4,'  Available :',LWORK4
         CALL QUIT('Insufficient memory in SOSDRV.4')
      ENDIF
C
C--------------------------------
C     5. Singlet work allocation:
C--------------------------------
C
      KEXCSN = KWORK4
      KSPPSO = KEXCSN + NEXCIS
      KWORK5 = KSPPSO + 9*NUCDEP*NUCDEP*NLVAR2
      LWORK5 = LWORK  - KWORK5
C
      IF (LWORK5 .LT. 0) THEN
         WRITE(LUPRI,*)'SOSDRV.5: Need :',KWORK5,'  Available :',LWORK5
         CALL QUIT('Insufficient memory in SOSDRV.5')
      ENDIF
C
C----------------------------------
C     6. Impresion work allocation:
C----------------------------------
C
      KSPTOT = KWORK5
      KWORK6 = KSPTOT + 3*NUCDEP*3*NUCDEP*NLVAR2
      LWORK6 = LWORK  - KWORK6
C
      IF (LWORK6 .LT. 0) THEN
         WRITE(LUPRI,*)'SOSDRV.6: Need :',KWORK6,'  Available :',LWORK6
         CALL QUIT('Insufficient memory in SOSDRV.6')
      ENDIF
C
      IF (SOSOCC) THEN
C
C
         IF (SOSOCS) THEN
C
C--------------------------------------
C     7. Loop over excitation energies:
C--------------------------------------
C
            NEXCII = 0
            NEXCIF = NVARPT
            NJUMPZ = NVARPT
            IF (NSTATI .NE. 0) NEXCII = ABS(NSTATI)
            IF ((NSTATF .NE. 0) .AND. (NSTATF .LE. NVARPT))
     &                         NEXCIF = ABS(NSTATF)
            IF (NITRST .NE. 0) NJUMPZ = ABS(NITRST)
CPFP Removed
C            IF (NSTATI .NE. 0) NEXCII = NSTATI
C            IF (NSTATF .NE. 0) NEXCIF = NSTATF
CCPFP I think this is wrong
CC            IF (NJUMPZ .NE. 0) NJUMPZ = NITRST
C            IF (NITRST .NE. 0) NJUMPZ = NITRST
C
            DO IEE = NEXCII,NEXCIF,NJUMPZ
C
               NEXCIS = IEE
               NEXCIT = IEE
CPFP Unnecessary
C               IF (IEE .GT. NVARPT) THEN
C                  NEXCIF = NVARPT
C                  NEXCIS = NVARPT
C                  NEXCIT = NVARPT
C               END IF
C
C--------------------------------
C     Loop over occupied orbitals
C--------------------------------
C
               JLOCCF = 1
               JLOCCL = NLOCC
C
               ILOCCF = 1
               ILOCCL = NLOCC
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
               CALL DZERO(WORK(KSD),9*9*NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KFC),NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KSDF),9*NUCDEP*NUCDEP*NLVAR2)
C
               CALL DZERO(WORK(KSPSD),9*NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KSPFC),9*NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KSPSDF),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate triplet contributions to coupling constants
C==========================================================
C
               CALL SOSJTR(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,
     &                     ILOCCL,JLOCCF,JLOCCL,NEXCIT,WORK(KCLOCO),
     &                     WORK(KCLOCV),WORK(KEXCTR),WORK(KSD),
     &                     WORK(KFC),WORK(KSDF),WORK(KSPSD),WORK(KSPFC),
     &                     WORK(KSPSDF),WORK(KWORK6),LWORK6)
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
              CALL DZERO(WORK(KSPPSO),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate singlet contributions to coupling constants
C==========================================================
C
               CALL SOSJSN(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,
     &                     ILOCCL,JLOCCF,JLOCCL,NEXCIS,WORK(KCLOCO),
     &                     WORK(KCLOCV),WORK(KEXCSN),WORK(KSPPSO),
     &                     WORK(KWORK6),LWORK6)
C
C=================================================
C     Print the results for the coupling constants
C=================================================
C
               CALL LSPRES(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,
     &                     JLOCCL,WORK(KSPDSO),WORK(KSPPSO),WORK(KSPSD),
     &                     WORK(KSPFC),WORK(KSPSDF),WORK(KSPTOT),
     &                     WORK(KWORK6),LWORK6)
C
               CALL TIMER ('SOSDRV',TIMEIN,TIMOUT)
               PASS   = .TRUE.
C
            END DO
C
CPFP 03-12-2015 adding a sort of automatic convergence
C
         ELSE IF (SOSOSC) THEN
C
CPFP 03-12-2015 It is necessary a vertor containing the J's for each step given
C               with the excitation energies.
C
C--------------------------------
C     Loop over occupied orbitals
C--------------------------------
C
c            JLOCCF = 1
c            JLOCCL = NLOCCS
cC
c            ILOCCF = 1
c            ILOCCL = NLOCCS
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
c            CALL DZERO(WORK(KSD),9*9*NUCDEP*NUCDEP*NLVAR2)
c            CALL DZERO(WORK(KFC),NUCDEP*NUCDEP*NLVAR2)
c            CALL DZERO(WORK(KSDF),9*NUCDEP*NUCDEP*NLVAR2)
cC
c            CALL DZERO(WORK(KSPSD),9*NUCDEP*NUCDEP*NLVAR2)
c            CALL DZERO(WORK(KSPFC),9*NUCDEP*NUCDEP*NLVAR2)
c            CALL DZERO(WORK(KSPSDF),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate triplet contributions to coupling constants
C==========================================================
C
c            CALL SOSJTR(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,
c     &                  ILOCCL,JLOCCF,JLOCCL,NEXCIT,WORK(KCLOCO),
c     &                  WORK(KCLOCV),WORK(KEXCTR),WORK(KSD),WORK(KFC),
c     &                  WORK(KSDF),WORK(KSPSD),WORK(KSPFC),WORK(KSPSDF),
c     &                  WORK(KWORK6),LWORK6)
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
c            CALL DZERO(WORK(KSPPSO),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate singlet contributions to coupling constants
C==========================================================
C
c            CALL SOSJSN(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,ILOCCL,
c     &                  JLOCCF,JLOCCL,NEXCIS,WORK(KCLOCO),WORK(KCLOCV),
c     &                  WORK(KEXCSN),WORK(KSPPSO),WORK(KWORK6),LWORK6)
C
C=================================================
C     Print the results for the coupling constants
C=================================================
C
c            CALL LSPRED(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,
c     &                  JLOCCL,WORK(KSPDSO),WORK(KSPPSO),WORK(KSPSD),
c     &                  WORK(KSPFC),WORK(KSPSDF),WORK(KSPTOT),
c     &                  WORK(KWORK6),LWORK6)
cC
c            CALL TIMER ('SOSDRV',TIMEIN,TIMOUT)
c            PASS   = .TRUE.
C
CPFP-TEST! AFTER THE FIRST Calculation.
C            WRITE(LUPRI,*)' After 1st. Calc. :',
C     &                    ' Need to select the wright coupling!'
Cend-Test
C
C    04-12-2015 It is needed to analyze the vector containing the previous J's
C               to see whether it converge.
C
C
            KNEETO = KWORK6
            KWORK7 = KNEETO + 5*NVARPT
            LWORK7 = LWORK  - KWORK7
C
            IF (LWORK7 .LT. 0) THEN
              WRITE(LUPRI,*)'SOSDRV.7: Need :',KWORK7,
     &                      '  Available :',LWORK7
              CALL QUIT('Insufficient memory in SOSDRV.7')
            ENDIF
C
            CALL DZERO(WORK(KWORK7),NVARPT)
C
CPFP Varables introduced in the input
C
            NOCCTA = 0
            NOMTHD = 0
            NOTHRS = 0.0D0
            IF (NCSSCC .NE. 0) NOCCTA = ABS(NCSSCC)
            IF (NCMTHD .NE. 0) NOMTHD = ABS(NCMTHD)
            IF (THRSLO .LT. 1.0D0) NOTHRS = ABS(THRSLO)
C
Cend-PFP
C
C--------------------------------------
C     7. Loop over excitation energies:
C--------------------------------------
C
            NEXCII = 0
            NEXCIF = NVARPT
            NJUMPZ = NVARPT
            IF (NSTATI .NE. 0) NEXCII = ABS(NSTATI)
            IF ((NSTATF .NE. 0) .AND. (NSTATF .LE. NVARPT))
     &                         NEXCIF = ABS(NSTATF)
            IF (NITRST .NE. 0) NJUMPZ = ABS(NITRST)
C
            DO IEE = NEXCII,NEXCIF,NJUMPZ
C
CPFP 17-01-2016 Needed for SOS-module
            ISINGUL = 0
Cend-PFP
               NEXCIS = IEE
               NEXCIT = IEE
CPFP
               INDEXT = IEE
Cend-PFP
CPFP  Unnecessary now with the "IF" above
C               IF (IEE .GT. NVARPT) THEN
C                  NEXCIF = NVARPT
C                  NEXCIS = NVARPT
C                  NEXCIT = NVARPT
C               END IF
C
C--------------------------------
C     Loop over occupied orbitals
C--------------------------------
C
               JLOCCF = 1
               JLOCCL = NLOCC
C
               ILOCCF = 1
               ILOCCL = NLOCC
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
               CALL DZERO(WORK(KSD),9*9*NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KFC),NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KSDF),9*NUCDEP*NUCDEP*NLVAR2)
C
               CALL DZERO(WORK(KSPSD),9*NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KSPFC),9*NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KSPSDF),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate triplet contributions to coupling constants
C==========================================================
C
               CALL SOSJTR(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,
     &                     ILOCCL,JLOCCF,JLOCCL,NEXCIT,WORK(KCLOCO),
     &                     WORK(KCLOCV),WORK(KEXCTR),WORK(KSD),
     &                     WORK(KFC),WORK(KSDF),WORK(KSPSD),WORK(KSPFC),
     &                     WORK(KSPSDF),WORK(KWORK6),LWORK6)
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
              CALL DZERO(WORK(KSPPSO),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate singlet contributions to coupling constants
C==========================================================
C
               CALL SOSJSN(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,
     &                     ILOCCL,JLOCCF,JLOCCL,NEXCIS,WORK(KCLOCO),
     &                     WORK(KCLOCV),WORK(KEXCSN),WORK(KSPPSO),
     &                     WORK(KWORK6),LWORK6)
C
C=================================================
C     Print the results for the coupling constants
C=================================================
C
               CALL LSPRES(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,
     &                     JLOCCL,WORK(KSPDSO),WORK(KSPPSO),WORK(KSPSD),
     &                     WORK(KSPFC),WORK(KSPSDF),WORK(KSPTOT),
     &                     WORK(KWORK6),LWORK6)
C
CPFP Here it must be called the subroutine which must analyze the vector containing the J's
C             INDEXT : Actual number of excitation energies
C             NEXCII : Initial number of excitation energies
C             NEXCIF : Final number of excitation energies
C             NJUMPZ : Step (integer)
C             NCSSCC=NOCCTA : Nuember of the SSCC of interest
C             NCMTHD=NOMTHD : Method to extrapolate (integer)
C             THRSLO=NOTHRS : Threshold for convergence (real*8)
C             RESSLT : Result in % (real*8)
C
c               WRITE(LUPRI,*) ' 1 - BEFORE CALLING SERCON !!!'
C
               CALL SERCON(NOCCTA,NOMTHD,NEXCII,NEXCIF,NJUMPZ,
     &                     NOTHRS,RESSLT,WORK(KWORK7),LWORK7)
               IF (DABS(RESSLT) .LT. DABS(THRSLO))  THEN
                 WRITE(LUPRI,'(/2(A,F8.4),A)')
     &           ' Convergence reached :',RESSLT,
     &           ' <',THRSLO,' %'
                 WRITE(LUPRI,'(/2(A,I5))')
     &           ' Number of steps', NCSSCC,' Method :', NCMTHD
                 WRITE(LUPRI,'(/(A))')
     &           ' Quit the loop over Excit. Ener.'
C
                 CALL TIMER ('SOSDRV',TIMEIN,TIMOUT)
                 PASS   = .TRUE.
                 GOTO 1010
C
               END IF
C
               CALL TIMER ('SOSDRV',TIMEIN,TIMOUT)
               PASS   = .TRUE.
C
            END DO
C
 1010       CONTINUE
Cend-PFP
         ELSE
C
C--------------------------------
C     Loop over occupied orbitals
C--------------------------------
C
            JLOCCF = 1
            JLOCCL = NLOCC
C
            ILOCCF = 1
            ILOCCL = NLOCC
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
            CALL DZERO(WORK(KSD),9*9*NUCDEP*NUCDEP*NLVAR2)
            CALL DZERO(WORK(KFC),NUCDEP*NUCDEP*NLVAR2)
            CALL DZERO(WORK(KSDF),9*NUCDEP*NUCDEP*NLVAR2)
C
            CALL DZERO(WORK(KSPSD),9*NUCDEP*NUCDEP*NLVAR2)
            CALL DZERO(WORK(KSPFC),9*NUCDEP*NUCDEP*NLVAR2)
            CALL DZERO(WORK(KSPSDF),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate triplet contributions to coupling constants
C==========================================================
C
            CALL SOSJTR(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,
     &                  ILOCCL,JLOCCF,JLOCCL,NEXCIT,WORK(KCLOCO),
     &                  WORK(KCLOCV),WORK(KEXCTR),WORK(KSD),WORK(KFC),
     &                  WORK(KSDF),WORK(KSPSD),WORK(KSPFC),WORK(KSPSDF),
     &                  WORK(KWORK6),LWORK6)
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
            CALL DZERO(WORK(KSPPSO),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate singlet contributions to coupling constants
C==========================================================
C
            CALL SOSJSN(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                  JLOCCF,JLOCCL,NEXCIS,WORK(KCLOCO),WORK(KCLOCV),
     &                  WORK(KEXCSN),WORK(KSPPSO),WORK(KWORK6),LWORK6)
C
C=================================================
C     Print the results for the coupling constants
C=================================================
C
            CALL LSPRES(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,
     &                  JLOCCL,WORK(KSPDSO),WORK(KSPPSO),WORK(KSPSD),
     &                  WORK(KSPFC),WORK(KSPSDF),WORK(KSPTOT),
     &                  WORK(KWORK6),LWORK6)
C
            CALL TIMER ('SOSDRV',TIMEIN,TIMOUT)
            PASS   = .TRUE.
C
         END IF
C
      ELSE IF (FBOVIR .OR. FBSETV .OR. FBSTVO) THEN
C
C--------------------------------
C     Loop over occupied orbitals
C--------------------------------
C
         DO  J = 1, NLOCC
C
            JLOCCF = J
            JLOCCL = J
C
            DO  I = 1, NLOCC
C
               ILOCCF = I
               ILOCCL = I
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
               CALL DZERO(WORK(KSD),9*9*NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KFC),NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KSDF),9*NUCDEP*NUCDEP*NLVAR2)
C
               CALL DZERO(WORK(KSPSD),9*NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KSPFC),9*NUCDEP*NUCDEP*NLVAR2)
               CALL DZERO(WORK(KSPSDF),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate triplet contributions to coupling constants
C==========================================================
C
               CALL SOSJTR(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,
     &                     ILOCCL,JLOCCF,JLOCCL,NEXCIT,WORK(KCLOCO),
     &                     WORK(KCLOCV),WORK(KEXCTR),WORK(KSD),
     &                     WORK(KFC),WORK(KSDF),WORK(KSPSD),WORK(KSPFC),
     &                     WORK(KSPSDF),WORK(KWORK6),LWORK6)
C
C-------------------------------------------------------
C     Initialize temporary arrays for coupling constants
C-------------------------------------------------------
C
               CALL DZERO(WORK(KSPPSO),9*NUCDEP*NUCDEP*NLVAR2)
C
C==========================================================
C     Calculate singlet contributions to coupling constants
C==========================================================
C
               CALL SOSJSN(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,
     &                     ILOCCL,JLOCCF,JLOCCL,NEXCIS,WORK(KCLOCO),
     &                     WORK(KCLOCV),WORK(KEXCSN),WORK(KSPPSO),
     &                     WORK(KWORK6),LWORK6)
C
C=================================================
C     Print the results for the coupling constants
C=================================================
C
               CALL LSPRES(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,
     &                     JLOCCL,WORK(KSPDSO),WORK(KSPPSO),WORK(KSPSD),
     &                     WORK(KSPFC),WORK(KSPSDF),WORK(KSPTOT),
     &                     WORK(KWORK6),LWORK6)
C
               CALL TIMER ('SOSDRV',TIMEIN,TIMOUT)
               PASS   = .TRUE.
C
            END DO
         END DO
C
      END IF
C
      CALL QEXIT('SOSDRV')
C
      RETURN
      END
C**********************************************************************
C  /* Deck soodso */
      SUBROUTINE SOODSO(IPRSOS,NLOCC,CLOCO,SPNDSO,WORK,LWORK)
C
C      Stephan P. A. Sauer 22/12-1999
C
C      This routine allocates some non-real arrays and calls the SOODSO_1,
C      the routine which calculates contributions from individual occupied 
C      orbitals to the diamagnetic spin orbit part of the 
C      spin-spin coupling constants.
C
#include "implicit.h"
C
C IRAT used from include file iratdef.h
C MXCENT used from include file mxcent.h
C LUPRI used from include file priunit.h
#include "iratdef.h"
#include "mxcent.h"
#include "priunit.h"
C
      DIMENSION   CLOCO(NOCCT,NLOCC)
      DIMENSION   SPNDSO(NLOCC,NLOCC,3*NUCDEP,3*NUCDEP)
      DIMENSION   WORK(LWORK)
C
C NOCCT used from COMMON /INFORB/
C NUCDEP used from COMMON /NUCLEI/
C
#include "inforb.h"
#include "nuclei.h"
C
      CALL QENTER('SOOSDO')
C
C------------------------------
C     Allocation of work memory
C------------------------------
C
      KIDSYM = 1
      KIDADR = KIDSYM + (9*MXCENT + 1)/IRAT
      KWORK1 = KIDADR + (9*MXCENT + 1)/IRAT
      LWORK1 = LWORK  - KWORK1
C
      IF (LWORK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'SOODSO.1: Need :',KWORK1,'  Available :',LWORK
         CALL QUIT('Insufficient memory in SOODSO.1')
      ENDIF
C
      CALL SOODSO_1(IPRSOS,NLOCC,CLOCO,SPNDSO,WORK(KIDSYM),
     &            WORK(KIDADR),WORK(KWORK1),LWORK1)
C
      CALL QEXIT('SOOSDO')
C
      RETURN
C
      END
C**********************************************************************
C  /* Deck soodso_1 */
      SUBROUTINE SOODSO_1(IPRSOS,NLOCC,CLOCO,SPNDSO,INTREP,INTADR,
     &                  WORK,LWORK)
C
C      Stephan P. A. Sauer 16/2-1999
C
C      This routine calculates contributions from individual occupied 
C      orbitals to the diamagnetic spin orbit part of the 
C      spin-spin coupling constants.
C
#include "implicit.h"
C
C MAXORB is needed in COMMON /INFIND/
C MAXASH is needed in COMMON /INFIND/
C MXCENT used from include file mxcent.h
C LUPRI used from include file priunit.h
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "priunit.h"
C
      PARAMETER (D0 = 0.0D+00, D1 = 1.0D+00, D2 = 2.0D+00)
      DIMENSION   CLOCO(NOCCT,NLOCC)
      DIMENSION   SPNDSO(NLOCC,NLOCC,3*NUCDEP,3*NUCDEP)
      DIMENSION   INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION   WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
      CHARACTER*8 LABEL
      CHARACTER*4 DSONUM
C
C NBAST, NNBASX, N2BASX, NBAS(8), NORB(8), NORBT, NCMOT, NOCCT, NSYM, 
C ICMO(8) used from COMMON /INFORB/
C ISX(MAXORB) used from COMMON /INFIND/
C MPQUAD used from COMMON /MAGONE/
C NUCDEP used from COMMON /NUCLEI/
C DODSO used from COMMON /SPNOUT/
C
#include "inforb.h"
#include "infind.h"
#include "magone.h"
#include "nuclei.h"
#include "spnout.h"
#include "inflin.h"
#include "inftap.h"
C
      CALL QENTER('SOODSO_1')
C
C==========================
C     Get the DSO integrals
C==========================
C
      NPATOM = 1
      NCOMP = (3*NUCDEP*(3*NUCDEP + 1)/2)
C
C------------------------------
C     Allocation of work memory
C------------------------------
C
      KPRSAO = 1
      KCMO   = KPRSAO + NNBASX * NCOMP
      KPRPAO = KCMO   + NCMOT
      KPRPMO = KPRPAO + N2BASX
      KPROCC = KPRPMO + NORBT * NORBT
      KPRLOC = KPROCC + NOCCT * NOCCT
      KDENMO = KPRLOC + NLOCC * NLOCC
      KDENLM = KDENMO + NOCCT * NOCCT
      KTEMP  = KDENLM + NLOCC * NLOCC
      KWORK1 = KTEMP  + NLOCC * NOCCT
      LWORK1 = LWORK  - KWORK1
C
      IF (LWORK1 .LT. 0) THEN
         WRITE(LUPRI,*)'SOODSO_1.1: Need :',KWORK1,'  Available :',LWORK
         CALL QUIT('Insufficient work memory in SOODSO.1')
      ENDIF
C
      IPRINT = IPRSOS - 20
C
      CALL GET1IN(WORK(KPRSAO),'DSO    ',NCOMP,WORK(KWORK1),LWORK1,
     &            LABINT,INTREP,INTADR,MPQUAD,.FALSE.,NPATOM,.TRUE.,
     &            DUMMY,.FALSE.,DUMMY,IPRINT)
C
C-------------------------
C     Read MO coefficients
C-------------------------
C
CSPAS: 27/5-11: adjusting to changes in the GPOPEN routine
C     LUSIFC = 0
CKeinSPASmehr
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &                                      IDUMMY,.FALSE.)
C
      REWIND LUSIFC
C
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
C
      READ (LUSIFC)
      READ (LUSIFC)
C
      CALL READT(LUSIFC,NCMOT,WORK(KCMO))
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C=============================
C     loop over all components
C=============================
C
      DO ICOMP = 1, NCOMP
C
         KCOMP  = KPRSAO+(ICOMP-1)*NNBASX
         LABEL  = LABINT(ICOMP)
         KSYM   = INTREP(ICOMP) + 1
C
C-----------------------------------
C        symmetrize AO integrals
C-----------------------------------
C
         CALL DSPTSI(NBAST,WORK(KCOMP),WORK(KPRPAO))
C
         IF (IPRSOS .GT. 15) THEN 
            CALL AROUND('Integrals of operator: '//LABEL)
            WRITE (LUPRI,'(A,I2)') ' Symmetry of operator:',KSYM
            CALL OUTPAK(WORK(KCOMP),NBAST,1,LUPRI)
            WRITE(LUPRI,'(/3A)') ' PROPERTY MATRIX ',LABEL,
     &                           ' IN AO BASIS AFTER SYMMETRIZATION '
            CALL OUTPUT(WORK(KPRPAO),1,NBAST,1,NBAST,NBAST,NBAST,1,
     &                  LUPRI)
         END IF
C
C-----------------------------
C        transform to MO basis
C-----------------------------
C
         CALL DZERO(WORK(KPRPMO),NORBT*NORBT) 
C
         DO ISYM = 1, NSYM
C
            JSYM = MULD2H(ISYM,KSYM)
C
            IF ((ISYM .GE. JSYM) .AND. 
     &          (NORB(ISYM) .GT. 0) .AND. (NORB(JSYM) .GT. 0)) THEN
C
               CALL UTHV(WORK(KCMO+ICMO(ISYM)),WORK(KPRPAO),
     &                   WORK(KCMO+ICMO(JSYM)),ISYM,JSYM,NBAS(ISYM),
     &                   NBAS(JSYM),WORK(KPRPMO),WORK(KWORK1))
C
               IF (IPRSOS. GT. 15) THEN 
                  WRITE(LUPR I,'(/,A,I5,A,I5)') ' ISYM= ',ISYM,
     &                  ' JSYM= ',JSYM
                  WRITE(LUPRI,'(/4A)') 
     &                  ' PROPERTY: ',LABEL,' IN MO. BASIS',
     &                  ' BEFORE (ANTI)SYMMETRIZATION '
                  CALL OUTPUT(WORK(KPRPMO),1,NORBT,1,NORBT,NORBT,NORBT,
     &                        1,LUPRI)
               END IF
C
               IF (IPRSOS .GT. 15) THEN 
                  WRITE(LUPRI,'(/A,I5,A)')
     &                  ' MO. COEFFICIENTS FOR SYMMETRY',ISYM
                  CALL OUTPUT(WORK(KCMO+ICMO(ISYM)),1,NBAS(ISYM),1,
     &                        NORB(ISYM),NBAS(ISYM),NORB(ISYM),1,LUPRI)
                  IF (ISYM .NE. JSYM) THEN
                     WRITE(LUPRI,'(/A,I5,A)')
     &                     ' MO. COEFFICIENTS FOR SYMMETRY',JSYM
                     CALL OUTPUT(WORK(KCMO+ICMO(JSYM)),1,NBAS(JSYM),1,
     &                           NORB(JSYM),NBAS(JSYM),NORB(JSYM),1,
     &                           LUPRI)
                  END IF
               END IF
C
            END IF
C
         END DO
C
         IF (KSYM .GT. 1) THEN
            ANTSYM = D1
CSPAS:27/3-06:TRANSA is changed to TRANSX
C           CALL TRANSA(WORK(KPRPMO),WORK(KPRPMO),NORBT,NORBT,ANTSYM)
            CALL TRANSX(WORK(KPRPMO),WORK(KPRPMO),NORBT,NORBT,ANTSYM,
     &                  IPRSOS)
CKeinSPASmehr
         ENDIF
C
         IF (IPRSOS .GE. 10) THEN
            WRITE(LUPRI,'(/4A)')' PROPERTY: ',LABEL,' IN MO. BASIS'
            CALL OUTPUT(WORK(KPRPMO),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                  LUPRI)
         END IF
C
C-------------------------------------------
C        get the integrals over occupied MOs
C-------------------------------------------
C
         CALL DZERO(WORK(KPROCC),NOCCT*NOCCT)
C
         DO JOCCT = 1, NOCCT
            JORBT = ISX(JOCCT)
            IOCOF = KPROCC - 1 + (JOCCT-1)*NOCCT
            IOROF = KPRPMO - 1 + (JORBT-1)*NORBT
C
            DO IOCCT = 1, NOCCT
               IORBT = ISX(IOCCT)
C
               WORK(IOCOF+IOCCT) = WORK(IOROF+IORBT)
C
            END DO
         END DO
C
         IF (IPRSOS .GE. 10) THEN
            WRITE(LUPRI,'(/4A)')' PROPERTY: ',LABEL,
     &                          ' in occupied MO basis'
            CALL OUTPUT(WORK(KPROCC),1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,
     &                  LUPRI)
         END IF
C
C------------------------------------------------
C        transform to localized occupied orbitals
C------------------------------------------------
C
         CALL DGEMM('N','N',NOCCT,NLOCC,NOCCT,D1,WORK(KPROCC),NOCCT,
     &              CLOCO,NOCCT,D0,WORK(KTEMP),NOCCT)
C
         CALL DGEMM('T','N',NLOCC,NLOCC,NOCCT,D1,CLOCO,NOCCT,
     &              WORK(KTEMP),NOCCT,D0,WORK(KPRLOC),NLOCC)
C
         IF (IPRSOS .GE. 10) THEN
            WRITE(LUPRI,'(/4A)')' PROPERTY: ',LABEL,
     &                          ' in localized occupied MO basis'
            CALL OUTPUT(WORK(KPRLOC),1,NLOCC,1,NLOCC,NLOCC,NLOCC,1,
     &                  LUPRI)
         END IF
C
C-----------------------------------------------
C        Hartree-Fock density matrix in AO basis
C-----------------------------------------------
C
         CALL DZERO(WORK(KDENMO),NOCCT*NOCCT)
C
         DO IOCCT = 1, NOCCT
C
            IOCOF = KDENMO - 1 + (IOCCT-1)*NOCCT
            WORK(IOCOF+IOCCT) = D2
C
         END DO
C
         IF (IPRSOS .GE. 10) THEN
            WRITE(LUPRI,'(/2A)')' density matrix in the basis of ',
     &                          ' occupied MOs'
            CALL OUTPUT(WORK(KDENMO),1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,
     &                  LUPRI)
         END IF
C
C------------------------------------------------
C        transform to localized occupied orbitals
C------------------------------------------------
C
         CALL DGEMM('N','N',NOCCT,NLOCC,NOCCT,D1,WORK(KDENMO),NOCCT,
     &              CLOCO,NOCCT,D0,WORK(KTEMP),NOCCT)
C
         CALL DGEMM('T','N',NLOCC,NLOCC,NOCCT,D1,CLOCO,NOCCT,
     &              WORK(KTEMP),NOCCT,D0,WORK(KDENLM),NLOCC)
C
         IF (IPRSOS .GE. 10) THEN
            WRITE(LUPRI,'(/2A)')' density matrix in the basis of ',
     &                          'localized occupied MOs'
            CALL OUTPUT(WORK(KDENLM),1,NLOCC,1,NLOCC,NLOCC,NLOCC,1,
     &                  LUPRI)
         END IF
C
C
C------------------------------------
C        multiply with density matrix
C------------------------------------
C
         DSONUM = LABEL(5:)
         CALL SCOMTR(DSONUM,ICOOR,JCOOR)
C
         IJ = -1
         DO J = 1, NLOCC
            DO I = 1, NLOCC
               IJ = IJ + 1
               SPNDSO(I,J,ICOOR,JCOOR) = - WORK(KDENLM+IJ) 
     &                                   * WORK(KPRLOC+IJ)
            END DO
         END DO
C
         CALL DCOPY(NLOCC*NLOCC,SPNDSO(1,1,ICOOR,JCOOR),1,
     &                          SPNDSO(1,1,JCOOR,ICOOR),1)

C
         IF (IPRSOS .GE. 5) THEN
            WRITE(LUPRI,'(/3A)') LABEL,' in localized occupied ',
     &                          'MO basis'
            CALL OUTPUT(SPNDSO(1,1,ICOOR,JCOOR),1,NLOCC,1,NLOCC,NLOCC,
     &                  NLOCC,1,LUPRI)
            SUM = 0
            DO I = 1, NLOCC
               SUM = SUM + SPNDSO(I,I,ICOOR,JCOOR)
            END DO
            WRITE(LUPRI,'(A,F12.4)') 'Total =',SUM
         END IF
C
      END DO
C
      CALL QEXIT('SOODSO_1')
C
      RETURN
C
      END
C
C**********************************************************************
C  /* Deck sosjsn */
C
      SUBROUTINE SOSJSN(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                  JLOCCF,JLOCCL,NEXCI,CLOCO,CLOCV,SNEXE,SPNPSO,
     &                  WORK,LWORK)
C
C      Stephan P. A. Sauer 13/11-1997
C
C      This routine calculates all singlet excitation energies 
C      at the RPA level and uses them to calculate contributions
C      from individual orbitals to the singlet parts of the 
C      spin-spin coupling constants as a sum over all excited states.
C
#include "implicit.h"
C
C
C IRAT used from include file iratdef.h
C MAXASH used from include file maxash.h
C MAXORB, MAXOCC used from include file maxorb.h
C MXCOOR used from include file mxcent.h
C LUPRI used from include file priunit.h
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "dummy.h"
C
      PARAMETER (ZERO = 0.0D+00, D05 = 0.5D+00, D1 = 1.0D+00)
C
#include "codata.h"
C
      DIMENSION CLOCO(NOCCT,NLOCC),CLOCV(NVIRT,NLVIR)
      DIMENSION SNEXE(NEXCI)
      DIMENSION SPNPSO(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION WORK(LWORK)
C
      LOGICAL HFCLC, CICLC, TRIPLE, EXECLC
      CHARACTER*8 LABEL
      LOGICAL     OLDDX
C
C MXRM, MXPHP, OOTV used from COMMON /CBIEXC/
C THRCLC, IPRCLC, MAXCLC used from COMMON /CBILRS1/
C NTRVEC used from COMMON /GDVEC/
C NVARMA used from COMMON /INFDIM/
C JTINAC, JTSEC, IOBTYP used from COMMON /INFIND/
C NVARPT used from COMMON /INFLIN/
C NASHT, NSYM, NOCCT, NVIRT used from COMMON /INFORB/
C NCONF modified in COMMON /INFVAR/
C JWOP used from COMMON /INFVAR/
C NUCDEP used from COMMON /NUCLEI/
C LUGDT, ABAGDT, LRGDT, NBGDT used from COMMON /NUCTAP/
C DOSD, DOSDFC, DOFC used from COMMON /SPNOUT/
C
#include "cbilrs.h"
#include "gdvec.h"
#include "infdim.h"
#include "infind.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "nuclei.h"
#include "spnout.h"
#include "inftap.h"
C
      CALL QENTER('SOSJSN')
C
CSPAS: 27/5-11: adjusting to changes in the GPOPEN routine
C     LUGDVES = 0
C     LUSOVES = 0
C     LUREVES = 0
      LUGDVES = -1
      LUSOVES = -1
      LUREVES = -1
CKeinSPASmehr
C
      IPRRSP = IPRCLC
      THRSOS = THRCLC
      MAXSOS = MAXCLC
C
C----------------------------------------
C     Check for Hartree-Fock wavefunction
C----------------------------------------
C
      HFCLC = NASHT .LE. 1
      IF (.NOT. HFCLC) THEN
         WRITE(LUPRI,'(/3(A))') 'SOSJSN : Sum Over States calculation',
     &                          ' only possible for HF or DFT',
     &                          ' wavefunctions'
         CALL QUIT('SOSJSN : SOS only possible for HF/DFT wavefunct.')
      ELSE
         NCONF = 1
      END IF
C
C------------------------------------------------------------
C     Open direct access file with the property gradients,
C     which were already calculated in the normal calculation 
C     of the spin-spin coupling constants
C------------------------------------------------------------
C
      CALL GPOPEN(LUGDI,ABAGDI,'UNKNOWN','DIRECT',' ',
     &            IRAT*NVARMA,OLDDX)
C
      CALL GPOPEN(LUSOVES,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVES,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVES,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
      IF (OLDDX) THEN
         WRITE (LUPRI,'(/A)') ' Old LUGDI file opened in SOSJSN.'
      ELSE
         WRITE (LUPRI,'(/A)') 'SOSJSN : Old LUGDI file not found'
         CALL QUIT('SOSJSN : Old LUGDI file not found')
      END IF
C
C*************************
C     Loop over symmetries
C*************************
C
      DO ISYM = 1, NSYM
         NREC = NGDVEC(ISYM,2)
         IF (NREC .GT. 0) THEN
C
            TRIPLE = .FALSE.
C
           CALL ABAVAR(ISYM,TRIPLE,IPRSOS,WORK,LWORK)
C
            IF (NVARPT .GT. 0) THEN
               IF (IPRSOS .GT. 5) THEN
                  WRITE (LUPRI,'(/,1X,A,I10)') 'NCONST ', NCONST
                  WRITE (LUPRI,'(1X,A,I10)')   'NVARPT ', NVARPT
                  WRITE (LUPRI,'(1X,A,I10)')   'for symmetry ',ISYM
               END IF
C
C--------------------------------
C              1. Work allocation
C--------------------------------
C
              IF ((ILOCCF .EQ. 1) .AND. (JLOCCF .EQ. 1)) THEN 
C
C=======================================================
C              Calculate all singlet excitation energies
C=======================================================
C
               IF (NEXCI .GT. NVARMA) THEN
                  WRITE(LUPRI,*) 
     &                 'SOSJSN: Symmetry ',ISYM,' Number of singlet ',
     &                 'excitations ',NEXCI,' exceeds maximum ',
     &                 'dimension ',NVARMA

                  CALL QUIT('Maximum dimension for excitations//
     &                       exceeded in SOSJSN')
               END IF
C
               CICLC  = .FALSE.
               EXECLC = .TRUE.
               NABATY = 1
               NABAOP = 1
               LABEL  = 'EXCITLAB' 
C
               CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     &                     SNEXE,NEXCI,NABATY,NABAOP,LABEL,
     &                     LUGDVES,LUSOVES,LUREVES,THRSOS,MAXSOS,IPRRSP,
     &                     MXRM,MXPHP,WORK,LWORK)
C
C-----------------------------------------------
C              Print singlet excitation energies
C-----------------------------------------------
C
               CALL HEADER('Singlet electronic excitation energies',15)
               WRITE (LUPRI,'(4X,A,/,4X,A,/,4X,A)')
     & ' Sym.      Mode              Frequency               Frequency',
     & 'ex. st.    No.               (au)                    (eV)',
     & '---------------------------------------------------------------'
               DO IEXVAL = 1, NEXCI
                  WRITE (LUPRI,'(5X,I2,5X,I5,2X,F22.5,3X,F22.5)')
     &                ISYM,IEXVAL,SNEXE(IEXVAL),SNEXE(IEXVAL)*XTEV
               END DO
C
               END IF
C
C--------------------------------
C              2. Work allocation
C--------------------------------
C
               NLVAR  = NLOCC*NLVIR
               NLVARS = NLOCCS*NLVIRS
C
               KTEMP  = 1
               KSLVZ  = KTEMP  + 2*NVARPT
               KSLVY  = KSLVZ  + NLVAR
               KGD1   = KSLVY  + NLVAR
               KGD2   = KGD1   + NLVAR
               KJLCON = KGD2   + NLVAR
               KWORK2 = KJLCON + NLVARS * NLVARS
               LWORK2 = LWORK  - KWORK2
C
               IF (LWORK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'SOSJSN.2: Need :',KWORK2,
     &                            '  Available :',LWORK
                  CALL QUIT('Insufficient memory in SOSJSN.2')
               ENDIF
C
C=====================================================
C              1. Loop over operators in this symmetry
C=====================================================
C
               DO IOP = 1, NREC
                  IREC   = IGDREC(IOP,ISYM,2)
                  ICOOR  = IGDCOR(IOP,ISYM,2)
C
                  CALL READDX (LUGDI,IREC,IRAT*NVARPT,WORK(KTEMP))
C
                  IF (IPRSOS .GT. 10) THEN
                     WRITE (LUPRI,'(/,A,I5)') ' 1. Property Gradient'
                     WRITE (LUPRI,'(A,I5)') ' Coordinate:', ICOOR
                     WRITE (LUPRI,'(A,I5)') ' Record:    ', IREC
                     CALL OUTPUT(WORK(KTEMP),1,NVARPT,1,1,NVARPT,1,1,
     &                           LUPRI)
                  END IF
C
C------------------------------------------------------------------
C                 transform property gradient to localized orbitals
C------------------------------------------------------------------
C
                  CALL LOCAI(NLOCC,NLVIR,NVARPT,CLOCO,CLOCV,
     &                       WORK(KTEMP),WORK(KGD1))
C
                  IF (IPRSOS .GT. 10) THEN
                     WRITE (LUPRI,'(/,A,I5)') 
     &                     ' localized 1. Property Gradient'
                     WRITE (LUPRI,'(A,I5)') ' Coordinate:', ICOOR
                     WRITE (LUPRI,'(A,I5)') ' Record:    ', IREC
                     CALL OUTPUT(WORK(KGD1),1,NLOCC*NLVIR,1,1,
     &                           NLOCC*NLVIR,1,1,LUPRI)
                  END IF
C
C========================================================
C                 2. Loop over operators in this symmetry
C========================================================
C
                  DO JOP = 1, IOP
                     JREC   = IGDREC(JOP,ISYM,2)
                     JCOOR  = IGDCOR(JOP,ISYM,2)
C
                     CALL READDX (LUGDI,JREC,IRAT*NVARPT,WORK(KTEMP))
C
                     IF (IPRSOS .GT. 10) THEN
                        WRITE (LUPRI,'(/,A,I5)') 
     &                        ' 2. Property Gradient'
                        WRITE (LUPRI,'(A,I5)') ' Coordinate:', JCOOR
                        WRITE (LUPRI,'(A,I5)') ' Record:    ', JREC
                        CALL OUTPUT(WORK(KTEMP),1,NVARPT,1,1,NVARPT,1,
     &                              1,LUPRI)
                     END IF
C
C---------------------------------------------------------------------
C                    transform property gradient to localized orbitals
C---------------------------------------------------------------------
C
                     CALL LOCAI(NLOCC,NLVIR,NVARPT,CLOCO,CLOCV,
     &                          WORK(KTEMP),WORK(KGD2))
C
                     IF (IPRSOS .GT. 10) THEN
                        WRITE (LUPRI,'(/,A,I5)') 
     &                         ' localized 2. Property Gradient'
                        WRITE (LUPRI,'(A,I5)') ' Coordinate:', JCOOR
                        WRITE (LUPRI,'(A,I5)') ' Record:    ', JREC
                        CALL OUTPUT(WORK(KGD2),1,NLOCC*NLVIR,1,1,
     &                              NLOCC*NLVIR,1,1,LUPRI)
                     END IF
C
C==========================================================
C                    Loop over excitations in this symmetry
C==========================================================
C
                    REWIND LUSOVES
C
                     CALL DZERO (WORK(KJLCON),NLVARS*NLVARS)
C
                     DO IEXVAL = 1, NEXCI
C
                        CALL READT(LUSOVES,2*NVARPT,WORK(KTEMP))
C
                        IF (IPRSOS .GT. 10) THEN
                           WRITE (LUPRI,'(/,A,I5)') ' Eigenvector'
                           CALL OUTPUT(WORK(KTEMP),1,NVARPT,1,2,NVARPT,
     &                                 2,1,LUPRI)
                        END IF
C
C------------------------------------------------------------------
C                       transform eigenvector to localized orbitals
C------------------------------------------------------------------
C
                        CALL LOCAI(NLOCC,NLVIR,NVARPT,CLOCO,CLOCV,
     &                             WORK(KTEMP),WORK(KSLVZ))
                        CALL LOCAI(NLOCC,NLVIR,NVARPT,CLOCO,CLOCV,
     &                             WORK(KTEMP+NVARPT),WORK(KSLVY))
C
                        IF (IPRSOS .GT. 10) THEN
                           WRITE (LUPRI,'(/,A,I5)') 
     &                            ' localized Eigenvector'
                           CALL OUTPUT(WORK(KSLVZ),1,NLOCC*NLVIR,1,1,
     &                         NLOCC*NLVIR,1,1,LUPRI)
                           CALL OUTPUT(WORK(KSLVY),1,NLOCC*NLVIR,1,1,
     &                         NLOCC*NLVIR,1,1,LUPRI)
                        END IF
C
C-------------------------------------------------------------------
C                       analyse only the contributions from occupied
C                       orbitals
C-------------------------------------------------------------------
C
                        IF (SOSOCC) THEN
C
C---------------------------------------------------------
C                       loop over localized ph excitations
C---------------------------------------------------------
C
                           IJLVAR = 0
C
                           DO JLOCC = JLOCCF, JLOCCL
C
                              JLOFF  = NLVIR * (JLOCC - 1)
C
                              IJLOFF = NLOCC * (JLOCC - 1)
C
                              DO ILOCC = ILOCCF, ILOCCL
C
                                 ILOFF  = NLVIR * (ILOCC - 1)
C
                                 IJLVAR = IJLOFF + ILOCC
                                 IJLCON = KJLCON - 1 + IJLVAR
C
                                 DO JLVIR = 1, NLVIR
C
                                    JLVAR = JLOFF + JLVIR
C
C                                   calculate contribution from
C                                   localized orbitals to the
C                                   2. transition moment
C
                                    TRMOM2 = WORK(KGD2-1+JLVAR)
     &                                     *  ( WORK(KSLVZ-1+JLVAR)
     &                                        + WORK(KSLVY-1+JLVAR))
C
                                    DO ILVIR = 1, NLVIR
C
                                       ILVAR  = ILOFF + ILVIR
C
C
C                                      calculate contribution from
C                                      localized orbitals to the
C                                      1. transition moment
C
                                       TRMOM1 = WORK(KGD1-1+ILVAR)
     &                                        * ( WORK(KSLVZ-1+ILVAR)
     &                                          + WORK(KSLVY-1+ILVAR))
C
C                                      calculate contribution from
C                                      localized orbitals to the
C                                      coupling constant
C
                                       TERM = TRMOM1 * TRMOM2 * D05
     &                                      / SNEXE(IEXVAL)
                                       WORK(IJLCON) = WORK(IJLCON)
     &                                              - TERM
C
                                       IF (IPRSOS .GT. 20) THEN
                                         WRITE (LUPRI,'(2X,A,I3,A,I3)')
     &                                      'For operator pair ',ICOOR,
     &                                      ' and',JCOOR
                                         WRITE (LUPRI,'(2X,A)')
     &                                      '========================'
                                         WRITE (LUPRI,'(5(A,D12.5))')
     &                                        'TRMOM1 :',TRMOM1,
     &                                        ' TRMOM2 :',TRMOM2,
     &                                        ' ENER :',
     &                                        SNEXE(IEXVAL),
     &                                        ' J :',term,' J :',
     &                                        WORK(IJLCON)
                                       END IF
C
                                    END DO
                                 END DO
C
                              END DO
                           END DO
C
C------------------------------------------------------------------
C                       analyse the contributions from occupied and
C                       virtual orbitals
C------------------------------------------------------------------
C
                        ELSE
C
C--------------------------------------------------------------
C                          loop over 2. localized ph excitation
C--------------------------------------------------------------
C
                           JLVAR  = 0
                           IJLVAR = 0
C
                           DO JLOCC = JLOCCF, JLOCCL
C
                              DO JLVIR = 1, NLVIR
C
                                 JLVAR = (JLOCC-1)*NLVIR + JLVIR
C
C                                calculate contribution from localized
C                                orbitals to the 2. transition moment
C
                                 TRMOM2 = WORK(KGD2-1+JLVAR)
     &                                  *  ( WORK(KSLVZ-1+JLVAR)
     &                                     + WORK(KSLVY-1+JLVAR))
C
C--------------------------------------------------------------------
C                                loop over 1. localized ph excitation
C--------------------------------------------------------------------
C
                                 ILVAR = 0
C
                                 DO ILOCC = ILOCCF, ILOCCL
C
                                    DO ILVIR = 1, NLVIR
C
                                       ILVAR = (ILOCC-1)*NLVIR + ILVIR
C
                                       JLVIRF = (JLVIR-1)*NLVIRS*NLOCCS
                                       IJLVAR =  JLVIRF + ILVIR
C
                                       IJLCON = KJLCON - 1 + IJLVAR
C
C                                      calculate contribution from 
C                                      localized orbitals to the 
C                                      1. transition moment 
C
                                       TRMOM1 = WORK(KGD1-1+ILVAR)
     &                                        * ( WORK(KSLVZ-1+ILVAR)
     &                                          + WORK(KSLVY-1+ILVAR))
C
C                                      calculate contribution from 
C                                      localized orbitals to the
C                                      coupling constant
C
                                       TERM = TRMOM1 * TRMOM2 * D05
     &                                      / SNEXE(IEXVAL)
                                       WORK(IJLCON) = WORK(IJLCON) 
     &                                              - TERM
C
                                       IF (IPRSOS .GT. 20) THEN
                                         WRITE (LUPRI,'(2X,A,I3,A,I3)')
     &                                      'For operator pair ',ICOOR,
     &                                      ' and',JCOOR
                                         WRITE (LUPRI,'(2X,A)')
     &                                      '========================'
                                         WRITE (LUPRI,'(5(A,D12.5))')
     &                                        'TRMOM1 :',TRMOM1,
     &                                        ' TRMOM2 :',TRMOM2,
     &                                        ' ENER :',
     &                                        SNEXE(IEXVAL),
     &                                        ' J :',term,' J :',
     &                                        WORK(IJLCON)
                                       END IF
C
                                    END DO
                                 END DO
C
                              END DO
                           END DO
C
                        END IF
C
                     END DO
C
C=================================================================
C                    End of loop over excitations in this symmetry
C=================================================================
C
                     IF (IPRSOS .GT. 10) THEN
                        CALL HEADER('spin-spin coupling tensors',-1)
                        WRITE (LUPRI,'(2X,A,I3,A,I3,/,2X,A)') 
     &                      'For operator pair ',ICOOR,' and',JCOOR,
     &                      '============================'
C
                        CALL PLOCPR(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                              JLOCCF,JLOCCL,WORK(KJLCON))
C
                     END IF
C
C=====================================================================
C                    Identify the operators and copy the contributions
C                    from localized orbitals to the coupling constant
C                    in the appropriate array
C=====================================================================
C
                     IF ((ICOOR .GT. 0) .AND. (JCOOR .GT. 0)) THEN
                        CALL DCOPY(NLVARS*NLVARS,WORK(KJLCON),1,
     &                             SPNPSO(1,ICOOR,JCOOR),1)
                        CALL DCOPY(NLVARS*NLVARS,WORK(KJLCON),1,
     &                             SPNPSO(1,JCOOR,ICOOR),1)
                     END IF
C
C
                  END DO
               END DO
C
C================================================================
C              End of double loop over operators in this symmetry
C================================================================
C
            END IF
         END IF
C
      END DO
C
C================================
C     End of loop over symmetries
C================================
C
      CALL GPCLOSE(LUGDI,'KEEP')
C
      CALL GPCLOSE(LUSOVES,'KEEP')
      CALL GPCLOSE(LUGDVES,'KEEP')
      CALL GPCLOSE(LUREVES,'KEEP')
C
      IF (IPRSOS .GT. 10) THEN
         CALL HEADER('PSO part of spin-spin coupling tensors',-1)
         DO IOP = 1, 3*NUCDEP
            DO JOP = 1, 3*NUCDEP
               WRITE (LUPRI,'(/,2X,A,I3,A,I3,/,2X,A)') 
     &                'For operator pair ',IOP,' and',JOP,
     &                '============================'
C
               CALL PLOCPR(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                     JLOCCF,JLOCCL,SPNPSO(1,IOP,JOP))
C
            END DO
         END DO
      END IF
C
C--------------------------------------------------
C     Print contributions to the coupling constants
C--------------------------------------------------
C
      IF (IPRSOS .GT. 5) THEN
         DO IATOM = 1, NUCDEP
            DO JATOM = 1, IATOM
C
               CALL TITLER(
     &              'Indirect spin-spin-coupling between atoms '//
     &             NAMDEP(IATOM)//' and '//NAMDEP(JATOM)//':','=',-1)
C
               CALL HEADER('Symmetry adapted PSO part of spin-spin '//
     &                     'coupling tensors',-1)
C
               CALL PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                     3*NUCDEP,IATOM,JATOM,SPNPSO,D1)
C
C
            END DO
         END DO
      END IF
C
      CALL QEXIT('SOSJSN')
C
      RETURN
C
      END
C*********************************************************************
C  /* Deck sosjtr */
C
      SUBROUTINE SOSJTR(IPRSOS,NLOCC,NLVIR,NLOCCS,NLVIRS,
     &                  ILOCCF,ILOCCL,JLOCCF,JLOCCL,NEXCI,
     &                  CLOCO,CLOCV,TREXE,SD,FC,SDF,SPNSD,
     &                  SPNFC,SPSDFC,WORK,LWORK)
C
C      Stephan P. A. Sauer 30/9-1997
C
C      This routine calculates all triplet excitation energies 
C      at the RPA level and uses them to calculate contributions
C      from individual orbitals to the triplet parts of the 
C      spin-spin coupling constants as a sum over all excited states.
C
#include "implicit.h"
C
C
C IRAT used from include file iratdef.h
C MAXASH used from include file maxash.h
C MAXORB, MAXOCC used from include file maxorb.h
C MXCOOR used from include file mxcent.h
C LUPRI used from include file priunit.h
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "dummy.h"
C
      PARAMETER (ZERO = 0.0D+00, D05 = 0.5D+00)
C
#include "codata.h"
C
      DIMENSION CLOCO(NOCCT,NLOCC),CLOCV(NVIRT,NLVIR)
      DIMENSION TREXE(NEXCI)
      DIMENSION SD(NLOCCS*NLOCCS*NLVIRS*NLVIRS,9*NUCDEP,9*NUCDEP)
      DIMENSION FC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,NUCDEP,NUCDEP)
      DIMENSION SDF(NLOCCS*NLOCCS*NLVIRS*NLVIRS,9*NUCDEP,NUCDEP)
      DIMENSION SPNSD(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPSDFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION WORK(LWORK)
C
      LOGICAL HFCLC, CICLC, TRIPLE, EXECLC
      LOGICAL OLDDX
      CHARACTER*8 LABEL
C
C THRTRP, IPRTRP, MAXTRP, MXRM, MXPHP, OOTV 
C used from COMMON /CBITRP/
C NTRVEC used from COMMON /GDVEC/
C NVARMA used from COMMON /INFDIM/
C JTINAC, JTSEC, IOBTYP used from COMMON /INFIND/
C NVARPT used from COMMON /INFLIN/
C NASHT, NSYM, NOCCT, NVIRT used from COMMON /INFORB/
C NCONF modified in COMMON /INFVAR/
C JWOP used from COMMON /INFVAR/
C NUCDEP used from COMMON /NUCLEI/
C LUGDT, ABAGDT, LRGDT, NBGDT used from COMMON /NUCTAP/
C DOSD, DOSDFC, DOFC used from COMMON /SPNOUT/
C
#include "cbitrp.h"
#include "gdvec.h"
#include "infdim.h"
#include "infind.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "nuclei.h"
#include "spnout.h"
#include "inftap.h"
C
      CALL QENTER('SOSJTR')
C
CSPAS: 27/5-11: adjusting to changes in the GPOPEN routine
C     LUGDVET = 0
C     LUSOVET = 0
C     LUREVET = 0
      LUGDVET = -1
      LUSOVET = -1
      LUREVET = -1
CKeinSPASmehr
C
      IPRRSP = IPRTRP
      THRSOS = THRTRP
      MAXSOS = MAXTRP
C
C----------------------------------------
C     Check for Hartree-Fock wavefunction
C----------------------------------------
C
      HFCLC = NASHT .LE. 1
      IF (.NOT. HFCLC) THEN
         WRITE(LUPRI,'(/3(A))') 'SOSJTR : Sum Over States calculation',
     &                          ' only possible for HF or DFT',
     &                          ' wavefunctions'
         CALL QUIT('SOSJTR : SOS only possible for HF/DFT wavefunct.')
      ELSE
         NCONF = 1
      END IF
C
C------------------------------------------------------------
C     Open direct access file with the property gradients,
C     which were already calculated in the normal calculation 
C     of the spin-spin coupling constants
C------------------------------------------------------------
C
      NREC = 0
      IF (DOSD .OR. DOSDFC) NREC = NREC + 9*NUCDEP
      IF (DOFC) NREC = NREC + NUCDEP
C
      CALL GPOPEN(LUGDT,ABAGDT,'UNKNOWN','DIRECT',' ',
     &            IRAT*NVARMA,OLDDX)
C
      CALL GPOPEN(LUSOVET,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVET,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVET,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
      IF (OLDDX) THEN
         WRITE (LUPRI,'(/A)') ' Old LUGDT file opened in SOSJTR.'
      ELSE
         WRITE (LUPRI,'(/A)') 'SOSJTR : Old LUGDT file not found'
         CALL QUIT('SOSJTR : Old LUGDT file not found')
      END IF
C
C*************************
C     Loop over symmetries
C*************************
C
      DO ISYM = 1, NSYM
         NREC = NTRVEC(ISYM)
         IF (NREC .GT. 0) THEN
C
            TRIPLE = .TRUE.
C
           CALL ABAVAR(ISYM,TRIPLE,IPRSOS,WORK,LWORK)
C
            IF (NVARPT .GT. 0) THEN
               IF (IPRSOS .GT. 5) THEN
                  WRITE (LUPRI,'(/,1X,A,I10)') 'NCONST ', NCONST
                  WRITE (LUPRI,'(1X,A,I10)')   'NVARPT ', NVARPT
                  WRITE (LUPRI,'(1X,A,I10)')   'for symmetry ',ISYM
               END IF
C
C--------------------------------
C              1. Work allocation
C--------------------------------
C
              IF ((ILOCCF .EQ. 1) .AND. (JLOCCF .EQ. 1)) THEN
C
C=======================================================
C              Calculate all triplet excitation energies
C=======================================================
C
               IF (NEXCI .GT. NVARMA) THEN
                  WRITE(LUPRI,*) 
     &                 'SOSJTR: Symmetry ',ISYM,' Number of triplet ',
     &                 'excitations ',NEXCI,' exceeds maximum ',
     &                 'dimension ',NVARMA

                  CALL QUIT('Maximum dimension for excitations//
     &                       exceeded in SOSJTR')
               END IF
C
               CICLC  = .FALSE.
               EXECLC = .TRUE.
               NABATY = 1
               NABAOP = 1
               LABEL  = 'EXCITLAB'
C
               CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     &                     TREXE,NEXCI,NABATY,NABAOP,LABEL,
     &                     LUGDVET,LUSOVET,LUREVET,THRSOS,MAXSOS,
     &                     IPRRSP,MXRM,MXPHP,WORK,LWORK)
C
C-----------------------------------------------
C              Print triplet excitation energies
C-----------------------------------------------
C
               CALL HEADER('Triplet electronic excitation energies',15)
               WRITE (LUPRI,'(4X,A,/,4X,A,/,4X,A)')
     & ' Sym.      Mode              Frequency               Frequency',
     & 'ex. st.    No.               (au)                    (eV)',
     & '---------------------------------------------------------------'
               DO IEXVAL = 1, NEXCI
                  WRITE (LUPRI,'(5X,I2,5X,I5,2X,F22.5,3X,F22.5)')
     &                ISYM,IEXVAL,TREXE(IEXVAL),TREXE(IEXVAL)*XTEV
               END DO
C
               END IF
C
C--------------------------------
C              2. Work allocation
C--------------------------------
C
               NLVAR  = NLOCC*NLVIR
               NLVARS = NLOCCS*NLVIRS
C
               KTEMP  = 1
               KSLVZ  = KTEMP  + 2*NVARPT
               KSLVY  = KSLVZ  + NLVAR
               KGD1   = KSLVY  + NLVAR
               KGD2   = KGD1   + NLVAR
               KJLCON = KGD2   + NLVAR
               KWORK2 = KJLCON + NLVARS * NLVARS
               LWORK2 = LWORK  - KWORK2
C
               IF (LWORK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'SOSJTR.2: Need :',KWORK2,
     &                            '  Available :',LWORK
                  CALL QUIT('Insufficient memory in SOSJTR.2')
               ENDIF
C
C=====================================================
C              1. Loop over operators in this symmetry
C=====================================================
C
               DO IOP = 1, NREC
                  IREC  = ITRREC(IOP,ISYM)
                  ICOOR = ITRCOR(IOP,ISYM)
C
CPFP  According to the GFORTRAN compiler :
C     Fortran runtime error: Cannot REWIND a file opened for DIRECT access
C     Then I commented it out
c                  REWIND LUGDT
C     For G77 there is not a problem
Cend-PFP
                  CALL READDX (LUGDT,IREC,IRAT*NVARPT,WORK(KTEMP))
C
                  IF (IPRSOS .GT. 10) THEN
                     WRITE (LUPRI,'(/,A,I5)') ' 1. Property Gradient'
                     WRITE (LUPRI,'(A,I5)') ' Coordinate:', ICOOR
                     WRITE (LUPRI,'(A,I5)') ' Record:    ', IREC
                     CALL OUTPUT(WORK(KTEMP),1,NVARPT,1,1,NVARPT,1,1,
     &                           LUPRI)
                  END IF
C
C------------------------------------------------------------------
C                 transform property gradient to localized orbitals
C------------------------------------------------------------------
C
                  CALL LOCAI(NLOCC,NLVIR,NVARPT,CLOCO,CLOCV,
     &                       WORK(KTEMP),WORK(KGD1))
C
                  IF (IPRSOS .GT. 10) THEN
                     WRITE (LUPRI,'(/,A,I5)') 
     &                     ' localized 1. Property Gradient'
                     WRITE (LUPRI,'(A,I5)') ' Coordinate:', ICOOR
                     WRITE (LUPRI,'(A,I5)') ' Record:    ', IREC
                     CALL OUTPUT(WORK(KGD1),1,NLOCC*NLVIR,1,1,
     &                           NLOCC*NLVIR,1,1,LUPRI)
                  END IF
C
C========================================================
C                 2. Loop over operators in this symmetry
C========================================================
C
                  DO JOP = 1, IOP
                     JREC  = ITRREC(JOP,ISYM)
                     JCOOR = ITRCOR(JOP,ISYM)
C
                     CALL READDX (LUGDT,JREC,IRAT*NVARPT,WORK(KTEMP))
C
                     IF (IPRSOS .GT. 10) THEN
                        WRITE (LUPRI,'(/,A,I5)') 
     &                        ' 2. Property Gradient'
                        WRITE (LUPRI,'(A,I5)') ' Coordinate:', JCOOR
                        WRITE (LUPRI,'(A,I5)') ' Record:    ', JREC
                        CALL OUTPUT(WORK(KTEMP),1,NVARPT,1,1,NVARPT,1,
     &                              1,LUPRI)
                     END IF
C
C---------------------------------------------------------------------
C                    transform property gradient to localized orbitals
C---------------------------------------------------------------------
C
                     CALL LOCAI(NLOCC,NLVIR,NVARPT,CLOCO,CLOCV,
     &                          WORK(KTEMP),WORK(KGD2))
C
                     IF (IPRSOS .GT. 10) THEN
                        WRITE (LUPRI,'(/,A,I5)') 
     &                         ' localized 2. Property Gradient'
                        WRITE (LUPRI,'(A,I5)') ' Coordinate:', JCOOR
                        WRITE (LUPRI,'(A,I5)') ' Record:    ', JREC
                        CALL OUTPUT(WORK(KGD2),1,NLOCC*NLVIR,1,1,
     &                              NLOCC*NLVIR,1,1,LUPRI)
                     END IF
C
C==========================================================
C                    Loop over excitations in this symmetry
C==========================================================
C
                     REWIND LUSOVET
C
                     CALL DZERO (WORK(KJLCON),NLVARS*NLVARS)
C
                     DO IEXVAL = 1, NEXCI
C
                        CALL READT(LUSOVET,2*NVARPT,WORK(KTEMP))
C
                        IF (IPRSOS .GT. 10) THEN
                           WRITE (LUPRI,'(/,A,I5,I5,I5)') 
     &                           ' Eigenvector no. :',IOP,JOP,IEXVAL
                           CALL OUTPUT(WORK(KTEMP),1,NVARPT,1,2,NVARPT,
     &                                 2,1,LUPRI)
                        END IF
C
C------------------------------------------------------------------
C                       transform eigenvector to localized orbitals
C------------------------------------------------------------------
C
                        CALL LOCAI(NLOCC,NLVIR,NVARPT,CLOCO,CLOCV,
     &                             WORK(KTEMP),WORK(KSLVZ))
                        CALL LOCAI(NLOCC,NLVIR,NVARPT,CLOCO,CLOCV,
     &                             WORK(KTEMP+NVARPT),WORK(KSLVY))
C
                        IF (IPRSOS .GT. 10) THEN
                           WRITE (LUPRI,'(/,A,I5,I5,I5)') 
     &                       ' localized Eigenvector',IOP,JOP,IEXVAL
                           CALL OUTPUT(WORK(KSLVZ),1,NLOCC*NLVIR,1,1,
     &                         NLOCC*NLVIR,1,1,LUPRI)
                           CALL OUTPUT(WORK(KSLVY),1,NLOCC*NLVIR,1,1,
     &                         NLOCC*NLVIR,1,1,LUPRI)
                        END IF
C
C-------------------------------------------------------------------
C                       analyse only the contributions from occupied
C                       orbitals
C-------------------------------------------------------------------
C
                        IF (SOSOCC) THEN
C
C---------------------------------------------------------
C                       loop over localized ph excitations
C---------------------------------------------------------
C
                           IJLVAR = 0
C
                           DO JLOCC = JLOCCF, JLOCCL
C
                              JLOFF  = NLVIR * (JLOCC - 1)
C
                              IJLOFF = NLOCC * (JLOCC - 1)
C
                              DO ILOCC = ILOCCF, ILOCCL
C
                                 ILOFF  = NLVIR * (ILOCC - 1)
C
                                 IJLVAR = IJLOFF + ILOCC
                                 IJLCON = KJLCON - 1 + IJLVAR
C
                                 DO JLVIR = 1, NLVIR
C
                                    JLVAR = JLOFF + JLVIR
C
C                                   calculate contribution from
C                                   localized orbitals to the
C                                   2. transition moment
C
                                    TRMOM2 = WORK(KGD2-1+JLVAR)
     &                                     *  (- WORK(KSLVZ-1+JLVAR)
     &                                         + WORK(KSLVY-1+JLVAR))
C
                                    DO ILVIR = 1, NLVIR
C
                                       ILVAR  = ILOFF + ILVIR
C
C
C                                      calculate contribution from
C                                      localized orbitals to the
C                                      1. transition moment
C
                                       TRMOM1 = WORK(KGD1-1+ILVAR)
     &                                        * (- WORK(KSLVZ-1+ILVAR)
     &                                           + WORK(KSLVY-1+ILVAR))
C
C                                      calculate contribution from 
C                                      localized orbitals to the
C                                      coupling constant
C
                                       TERM = TRMOM1 * TRMOM2 * D05
     &                                      / TREXE(IEXVAL)
                                       WORK(IJLCON) = WORK(IJLCON) 
     &                                              - TERM
C
                                       IF (IPRSOS .GT. 20) THEN
                                         WRITE (LUPRI,'(2X,A,I3,A,I3)')
     &                                      'For operator pair ',ICOOR,
     &                                      ' and',JCOOR
                                         WRITE (LUPRI,'(2X,A)')
     &                                      '========================'
                                         WRITE (LUPRI,'(5(A,D12.5))')
     &                                        'TRMOM1 :',TRMOM1,
     &                                        ' TRMOM2 :',TRMOM2,
     &                                        ' ENER :',
     &                                        TREXE(IEXVAL),
     &                                        ' J :',term,' J :',
     &                                        WORK(IJLCON)
                                       END IF
C
                                    END DO
                                 END DO
C
                              END DO
                           END DO
C
C------------------------------------------------------------------
C                       analyse the contributions from occupied and
C                       virtual orbitals
C------------------------------------------------------------------
C
                        ELSE
C
C--------------------------------------------------------------
C                          loop over 2. localized ph excitation
C--------------------------------------------------------------
C
                           JLVAR  = 0
                           IJLVAR = 0
C
                           DO JLOCC = JLOCCF, JLOCCL
C
                              DO JLVIR = 1, NLVIR
C
                                 JLVAR = (JLOCC-1)*NLVIR + JLVIR
C
C                                calculate contribution from localized
C                                orbitals to the 2. transition moment
C
                                 TRMOM2 = WORK(KGD2-1+JLVAR)
     &                                    *( - WORK(KSLVZ-1+JLVAR)
     &                                       + WORK(KSLVY-1+JLVAR) )
C
C--------------------------------------------------------------------
C                                loop over 1. localized ph excitation
C--------------------------------------------------------------------
C
                                 ILVAR = 0
C
                                 DO ILOCC = ILOCCF, ILOCCL
C
                                    DO ILVIR = 1, NLVIR
C
                                       ILVAR = (ILOCC-1)*NLVIR + ILVIR
C
                                       JLVIRF = (JLVIR-1)*NLVIRS*NLOCCS
                                       IJLVAR =  JLVIRF + ILVIR
C
                                       IJLCON = KJLCON - 1 + IJLVAR
C
C                                      calculate contribution from
C                                      localized orbitals to the
C                                      1. transition moment
C
                                       TRMOM1 = WORK(KGD1-1+ILVAR)
     &                                        * (- WORK(KSLVZ-1+ILVAR)
     &                                           + WORK(KSLVY-1+ILVAR))
C
C                                      calculate contribution from 
C                                      localized orbitals to the
C                                      coupling constant
C
                                       TERM = TRMOM1 * TRMOM2 * D05
     &                                      / TREXE(IEXVAL)
                                       WORK(IJLCON) = WORK(IJLCON) 
     &                                              - TERM
C
                                       IF (IPRSOS .GT. 20) THEN
                                         WRITE (LUPRI,'(2X,A,I3,A,I3)')
     &                                      'For operator pair ',ICOOR,
     &                                      ' and',JCOOR
                                         WRITE (LUPRI,'(2X,A)')
     &                                      '========================'
                                         WRITE (LUPRI,'(5(A,D12.5))')
     &                                      'TRMOM1 :',TRMOM1,
     &                                      ' TRMOM2 :',TRMOM2,
     &                                      ' ENER :',
     &                                      TREXE(IEXVAL),
     &                                      ' J :',term,' J :',
     &                                      WORK(IJLCON)
                                       END IF
C
                                    END DO
                                 END DO
C
                              END DO
                           END DO
C
C------------------------------------------------------------------
C                       End of loops over ph
C------------------------------------------------------------------
C
                        END IF
                     END DO
C
C=================================================================
C                    End of loop over excitations in this symmetry
C=================================================================
C
                     IF (IPRSOS .GT. 10) THEN
                        CALL HEADER('spin-spin coupling tensors',-1)
                        WRITE (LUPRI,'(2X,A,I3,A,I3,/,2X,A)') 
     &                      'For operator pair ',ICOOR,' and',JCOOR,
     &                      '============================'
C
                        CALL PLOCPR(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                              JLOCCF,JLOCCL,WORK(KJLCON))
C
                     END IF
C
C=====================================================================
C                    Identify the operators and copy the contributions
C                    from localized orbitals to the coupling constant
C                    in the appropriate array
C=====================================================================
C
                     IF ((ICOOR .GT. 0) .AND. (JCOOR .GT. 0)) THEN
                        CALL DCOPY(NLVARS*NLVARS,WORK(KJLCON),1,
     &                             SD(1,ICOOR,JCOOR),1)
                        CALL DCOPY(NLVARS*NLVARS,WORK(KJLCON),1,
     &                             SD(1,JCOOR,ICOOR),1)
                     ELSE IF ((ICOOR .LT. 0) .AND. (JCOOR .LT. 0)) THEN
                        CALL DCOPY(NLVARS*NLVARS,WORK(KJLCON),1,
     &                             FC(1,-ICOOR,-JCOOR),1)
                        CALL DCOPY(NLVARS*NLVARS,WORK(KJLCON),1,
     &                             FC(1,-JCOOR,-ICOOR),1)
                     ELSE IF ((ICOOR .LT. 0) .AND. (JCOOR .GT. 0)) THEN
                        CALL DCOPY(NLVARS*NLVARS,WORK(KJLCON),1,
     &                             SDF(1,JCOOR,-ICOOR),1)
                     END IF
C
C
                  END DO
               END DO
C
C================================================================
C              End of double loop over operators in this symmetry
C================================================================
C
            END IF
         END IF
C
      END DO
C
C================================
C     End of loop over symmetries
C================================
C
      CALL GPCLOSE(LUGDT,'KEEP')
C
      CALL GPCLOSE(LUSOVET,'KEEP')
      CALL GPCLOSE(LUGDVET,'KEEP')
      CALL GPCLOSE(LUREVET,'KEEP')
C
      IF (IPRSOS .GT. 10) THEN
         CALL HEADER('FC part of spin-spin coupling tensors',-1)
         DO IOP = 1, NUCDEP
            DO JOP = 1, NUCDEP
               WRITE (LUPRI,'(/,2X,A,I3,A,I3,/,2X,A)')
     &                'For operator pair ',IOP,' and',JOP,
     &                '============================'
C
               CALL PLOCPR(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                     JLOCCF,JLOCCL,FC(1,IOP,JOP))
C
            END DO
         END DO
C
C
         CALL HEADER('SD(+FC) part of spin-spin coupling tensors',-1)
         DO IOP = 1, 9*NUCDEP
            DO JOP = 1, 9*NUCDEP
               WRITE (LUPRI,'(/,2X,A,I3,A,I3,/,2X,A)') 
     &                'For operator pair ',IOP,' and',JOP,
     &                '============================'
C
               CALL PLOCPR(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                     JLOCCF,JLOCCL,SD(1,IOP,JOP))
C
            END DO
         END DO
C
C
         CALL HEADER('SD-FC part of spin-spin coupling tensors',-1)
         DO IOP = 1, 9*NUCDEP
            DO JOP = 1, NUCDEP
               WRITE (LUPRI,'(/,2X,A,I3,A,I3,/,2X,A)') 
     &                'For operator pair ',IOP,' and',JOP,
     &                '============================'
C
               CALL PLOCPR(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                     JLOCCF,JLOCCL,SDF(1,IOP,JOP))
C
            END DO
         END DO
      END IF
C
C===========================================
C     Collect triplet contributions 
C     to symmety adapated coupling constants
C===========================================
C
      CALL LTRPMV(NLOCCS,NLVIRS,IPRSOS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &            SPNSD,SPNFC,SPSDFC,SD,FC,SDF)
C
C====================================
C     Transform to non-symmetry basis
C====================================
C This was commented by somebody else (PFP)
c      CALL LSPSYM(NLOCC,NLVIRS,IPRSOS,SPNSD,SPNFC,SPSDFC,WORK,LWORK)
C
C
      CALL QEXIT('SOSJTR')
C
      RETURN
C
      END
C***********************************************************************
C  /* Deck ltrpmv */
C
      SUBROUTINE LTRPMV(NLOCCS,NLVIRS,IPRSOS,ILOCCF,ILOCCL,
     &                  JLOCCF,JLOCCL,SPNSD,SPNFC,SPSDFC,WRKSD,
     &                  WRKFC,WRKSDF)
C
C      Stephan P. A. Sauer 29/10-1997
C
C      This routine collects the different contributions to the elements
C      of the spin-dipolar, Fermi contact and Fermi contact spin-dipolar
C      cross term tensors
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
CPFP
C#include "mxorb.h"
#include "maxorb.h"
Cend-PFP
      PARAMETER (ZERO = 0.0D+00, D1 = 1.0D+00)
C
      DIMENSION WRKSD(NLOCCS*NLOCCS*NLVIRS*NLVIRS,9*NUCDEP,9*NUCDEP)
      DIMENSION WRKFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,NUCDEP,NUCDEP)
      DIMENSION WRKSDF(NLOCCS*NLOCCS*NLVIRS*NLVIRS,9*NUCDEP,NUCDEP)
      DIMENSION SPNSD(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPSDFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
C      DIMENSION TJ(3,3)
C
#include "nuclei.h"
#include "cbitrp.h"
#include "spnout.h"
#include "symmet.h"

C
      CALL QENTER('LTRPMV')
C
C=======================================
C     Sum the spin-dipolar term over the
C     components of the electronic spin
C=======================================
C
C no need! it was initialized (zeroed) in sosdrv
c      CALL DZERO(SPNSD,MXCOOR*MXCOOR*NLOCCS*NLOCCS*NLVIRS*NLVIRS)
      IF (DOSD .OR. DOSDFC) THEN
         JATOM1 = 0
         DO IREP = 0, MAXREP
         DO IATOM1 = 1, NUCIND
            IF (IAND(IREP,ISTBNU(IATOM1)).EQ.0) THEN
               JATOM1 = JATOM1 + 1
               JATOM2 = 0
               DO  IREP1 = 0, MAXREP
               DO  IATOM2 = 1, NUCIND
                  IF (IAND(IREP1,ISTBNU(IATOM2)).EQ.0) THEN
                     JATOM2 = JATOM2 + 1
                     DO ICOOR = 1, 3
                        IREP2 = IEOR(ISYMAX(ICOOR,2),IREP)
                        ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR,IREP2,2)
                        DO JCOOR = 1, 3
                           IREP3 = IEOR(ISYMAX(JCOOR,2),IREP1)
                           ISCOR2 = IPTCNT(3*(IATOM2-1)+JCOOR,IREP3,2)
                           IF (ISCOR1 .GT. 0 .AND. ISCOR2 .GT. 0) THEN
                              DO KCOOR = 1, 3
                                 IF (KCOOR .GT. ICOOR) THEN
                                   IREP4 = IEOR(ISYMAX(KCOOR,2),IREP)
                                   ISCOR3 = IPTCNT(3*(IATOM1-1)+KCOOR,
     &                                             IREP4,2)
                                   IADR1 = 3*(ISCOR3 - 1) + ICOOR
                                 ELSE
                                    IADR1 = 3*(ISCOR1 - 1) + KCOOR
                                 END IF
                                 IF (KCOOR .GT. JCOOR) THEN
                                   IREP4 = IEOR(ISYMAX(KCOOR,2),IREP1)
                                   ISCOR3 = IPTCNT(3*(IATOM2-1)+KCOOR,
     &                                             IREP4,2)
                                   IADR2 = 3*(ISCOR3 - 1) + JCOOR
                                 ELSE
                                    IADR2 = 3*(ISCOR2 - 1) + KCOOR
                                 END IF
                             DO ILVAR = 1, NLOCCS*NLOCCS*NLVIRS*NLVIRS
                                SPNSD(ILVAR,ISCOR1,ISCOR2) =
     &                                    SPNSD(ILVAR,ISCOR1,ISCOR2)
     &                                  + WRKSD(ILVAR,IADR1,IADR2)
                             END DO
                              END DO
                           END IF
                  END DO
               END DO
                  END IF
               END DO
            END DO
            END IF
            END DO
         END DO
      END IF
C
C
C====================================
C     Collect the Fermi contact terms
C====================================
C
      IF (DOFC) THEN
         JATOM1 = 0
         DO IREP = 0, MAXREP
            KATOM1 = JATOM1
            DO IATOM1 = 1, NUCIND
               IF (IAND(IREP,ISTBNU(IATOM1)).EQ.0) THEN
               JATOM1 = JATOM1 + 1
               JATOM2 = KATOM1
               DO IATOM2 = 1, NUCIND
                  IF (IAND(IREP,ISTBNU(IATOM2)).EQ.0) THEN
                  JATOM2 = JATOM2 + 1
                  DO ICOOR = 1, 3
                     IREP2 = IEOR(ISYMAX(ICOOR,2),IREP)
                     ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR,IREP2,2)
                     ISCOR2 = IPTCNT(3*(IATOM2 - 1) + ICOOR,IREP2,2)
                     IF (ISCOR1.GT.0 .AND. ISCOR2.GT.0) THEN
C
                        DO ILVAR = 1, NLOCCS*NLOCCS*NLVIRS*NLVIRS
C
                           SPNFC(ILVAR,ISCOR1,ISCOR2) = 
     &                          WRKFC(ILVAR,JATOM1,JATOM2)
                        END DO
                     END IF
                  END DO
                  END IF
               END DO
               END IF
            END DO
         END DO
      END IF
C
C====================================================
C     Sum the Fermi-contact - spin-dipolar cross term
C     over the  components of the electronic spin
C====================================================
C
      IF (DOFC .AND. DOSD) THEN
         JATOM1 = 0
         DO IREP = 0, MAXREP
            DO IATOM1 = 1, NUCIND
               IF (IAND(IREP,ISTBNU(IATOM1)).EQ.0) THEN
               JATOM1 = JATOM1 + 1
               JATOM2 = 0
               DO IREP1 = 0, MAXREP
                  DO IATOM2 = 1, NUCIND
                     IF (IAND(IREP1,ISTBNU(IATOM2)).EQ.0) THEN
                     JATOM2 = JATOM2 + 1
C
                     DO ILVAR = 1, NLOCCS*NLOCCS*NLVIRS*NLVIRS
C
                        DO ICOOR = 1, 3
                           IREP2 = IEOR(ISYMAX(ICOOR,2),IREP)
                           ISCOR1 = IPTCNT(3*(IATOM1-1)+ICOOR,IREP2,2)
                           DO JCOOR = 1, 3
                              IREP3  = IEOR(ISYMAX(JCOOR,2),IREP1)
                              ISCOR2 =
     &                            IPTCNT(3*(IATOM2-1)+JCOOR,IREP3,2)
                              IF (ISCOR1.GT.0 .AND. ISCOR2.GT.0) THEN
                                 IADR1  = 3*(ISCOR1 - 1) + JCOOR
                                 IREP4  = IEOR(ISYMAX(ICOOR,2),IREP1)
                                 ISCOR3 = 
     &                               IPTCNT(3*(IATOM2-1)+ICOOR,IREP4,2)
                                 IADR2 = 3*(ISCOR3 - 1) + JCOOR
                                 SPSDFC(ILVAR,ISCOR1,ISCOR2) = 
     &                                 WRKSDF(ILVAR,IADR1,JATOM2)
     &                               + WRKSDF(ILVAR,IADR2,JATOM1)
                              END IF
                           END DO
                        END DO
                     END DO
                     END IF
                  END DO
               END DO
               END IF
            END DO
         END DO
C
      END IF
C
C
C-------------------------------------------------------------------
C     Print contributions to the symmetry adapted coupling constants
C-------------------------------------------------------------------
C
      IF (IPRSOS .GT. 5) THEN
         DO IATOM = 1, NUCDEP
            DO JATOM = 1, IATOM
C
               CALL TITLER(
     &              'Indirect spin-spin-coupling between atoms '//
     &             NAMDEP(IATOM)//' and '//NAMDEP(JATOM)//':','=',-1)
C 
              IF (DOFC) THEN
                  CALL HEADER('Symmetry adapted FC part of spin-spin'//
     &                        ' coupling tensors',-1)
C
                  CALL PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                        3*NUCDEP,IATOM,JATOM,SPNFC,D1)
C
               END IF
C
               IF (DOSD .OR. DOSDFC) THEN
                  CALL HEADER(
     &                 'Symmetry adapted SD(+FC) part of spin-spin '//
     &                 'coupling tensors',-1)
C
                  CALL PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                        3*NUCDEP,IATOM,JATOM,SPNSD,D1)
C
               END IF
C
               IF (DOFC .AND. DOSD) THEN
                  CALL HEADER('Symmetry adapted SD-FC part of '//
     &                        'spin-spin coupling tensors',-1)
C
                  CALL PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                        3*NUCDEP,IATOM,JATOM,SPSDFC,D1)
C
               END IF
C
            END DO
         END DO
      END IF
C
      CALL QEXIT('LTRPMV')
C
      RETURN
C
      END
C*********************************************************************
C  /* Deck lspsym */
C
      SUBROUTINE LSPSYM(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                  JLOCCF,JLOCCL,SPNDSO,SPNPSO,SPNSD,
     &                  SPNFC,SPSDFC,WORK,LWORK)
C
C      Stephan P. A. Sauer 13/11-1997
C
C      This routine removes the symmetry adaptation 
C      of the nuclear magnetic moments.
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
#include "priunit.h"
C
C MXCENT, MXCOOR used from include file mxcent.h
#include "mxcent.h"
C
C NUCDEP used from COMMON /NUCLEI/
#include "nuclei.h"
C
      PARAMETER (D1 = 1.0D+00)
C
      DIMENSION SPNDSO(NLOCCS*NLOCCS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNPSO(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNSD(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPSDFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION WORK(LWORK)
C
      DIMENSION TJ(3,3)
C
      CALL QENTER('LSPSYM')
C
C--------------------------------
C     1. Work allocation
C--------------------------------
C
      KSDSO  = 1
      KSPSO  = KSDSO  + 9*NUCDEP*NUCDEP
      KSFC   = KSPSO  + 9*NUCDEP*NUCDEP
      KSSD   = KSFC   + 9*NUCDEP*NUCDEP
      KSSDFC = KSSD   + 9*NUCDEP*NUCDEP
      KCDSO  = KSSDFC + 9*NUCDEP*NUCDEP
      KCPSO  = KCDSO  + 9*NUCDEP*NUCDEP
      KCFC   = KCPSO  + 9*NUCDEP*NUCDEP
      KCSD   = KCFC   + 9*NUCDEP*NUCDEP
      KCSDFC = KCSD   + 9*NUCDEP*NUCDEP
      KCSTRA = KCSDFC + 9*NUCDEP*NUCDEP
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KWORK1 = KSCTRA + 9*NUCDEP*NUCDEP
      LWORK1 = LWORK  - KWORK1
C
      IF (LWORK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'LSPSYM.1: Need :',KWORK1,'  Available :',LWORK
         CALL QUIT('Insufficient memory in LSPSYM.1')
      ENDIF
C
C-------------------------------------------------------------
C     Transform individual contributions to non-symmetry basis
C-------------------------------------------------------------
C
         LSPRR = NLOCCS*NLOCCS
C
      DO IJLVAR = 1, LSPRR
C
         IJC = 0
         DO JCOOR = 1, 3*NUCDEP
            DO ICOOR = 1, 3*NUCDEP
               IJC = IJC + 1
               WORK(KSDSO-1+IJC)  = SPNDSO(IJLVAR,ICOOR,JCOOR)
            END DO
         END DO
C
         CALL DZERO(WORK(KCDSO),9*NUCDEP*NUCDEP)
C
         CALL TRAHES(WORK(KSDSO),3*NUCDEP,WORK(KCDSO),WORK(KCSTRA),
     &               WORK(KSCTRA),3*NUCDEP,3*NUCDEP,2)
C
         IJC = 0
         DO JCOOR = 1, 3*NUCDEP
            DO ICOOR = 1, 3*NUCDEP
               IJC = IJC + 1
               SPNDSO(IJLVAR,ICOOR,JCOOR) = WORK(KCDSO-1+IJC)
            END DO
         END DO
C
      END DO
C
C
C
      DO IJLVAR = 1, NLOCCS*NLOCCS*NLVIRS*NLVIRS
C
         IJC = 0
         DO JCOOR = 1, 3*NUCDEP
            DO ICOOR = 1, 3*NUCDEP
               IJC = IJC + 1
               WORK(KSPSO-1+IJC)  = SPNPSO(IJLVAR,ICOOR,JCOOR)
               WORK(KSFC-1+IJC)   = SPNFC(IJLVAR,ICOOR,JCOOR)
               WORK(KSSD-1+IJC)   = SPNSD(IJLVAR,ICOOR,JCOOR)
               WORK(KSSDFC-1+IJC) = SPSDFC(IJLVAR,ICOOR,JCOOR)
            END DO
         END DO
C
         CALL DZERO(WORK(KCPSO),9*NUCDEP*NUCDEP)
         CALL DZERO(WORK(KCFC),9*NUCDEP*NUCDEP)
         CALL DZERO(WORK(KCSD),9*NUCDEP*NUCDEP)
         CALL DZERO(WORK(KCSDFC),9*NUCDEP*NUCDEP)
C
         CALL TRAHES(WORK(KSPSO),3*NUCDEP,WORK(KCPSO),WORK(KCSTRA),
     &               WORK(KSCTRA),3*NUCDEP,3*NUCDEP,2)
         CALL TRAHES(WORK(KSFC),3*NUCDEP,WORK(KCFC),WORK(KCSTRA),
     &               WORK(KSCTRA),3*NUCDEP,3*NUCDEP,2)
         CALL TRAHES(WORK(KSSD),3*NUCDEP,WORK(KCSD),WORK(KCSTRA),
     &               WORK(KSCTRA),3*NUCDEP,3*NUCDEP,2)
         CALL TRAHES(WORK(KSSDFC),3*NUCDEP,WORK(KCSDFC),WORK(KCSTRA),
     &               WORK(KSCTRA),3*NUCDEP,3*NUCDEP,2)
C
         IJC = 0
         DO JCOOR = 1, 3*NUCDEP
            DO ICOOR = 1, 3*NUCDEP
               IJC = IJC + 1
               SPNPSO(IJLVAR,ICOOR,JCOOR) = WORK(KCPSO-1+IJC)
               SPNFC(IJLVAR,ICOOR,JCOOR)  = WORK(KCFC-1+IJC)
               SPNSD(IJLVAR,ICOOR,JCOOR)  = WORK(KCSD-1+IJC)
               SPSDFC(IJLVAR,ICOOR,JCOOR) = WORK(KCSDFC-1+IJC)
            END DO
         END DO
C
      END DO
C
C--------------------------------------------------
C     Print contributions to the coupling constants
C--------------------------------------------------
C
      IF (IPRSOS .GT. 5) THEN
         DO IATOM = 1, NUCDEP
            DO JATOM = 1, IATOM
C
               CALL TITLER(
     &              'Indirect spin-spin-coupling between atoms '//
     &             NAMDEP(IATOM)//' and '//NAMDEP(JATOM)//':','=',-1)
C
               CALL HEADER('DSO part of spin-spin coupling tensors',-1)
C
               CALL PLOPRT(NLOCCS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                     3*NUCDEP,IATOM,JATOM,SPNDSO,D1)
C
               CALL HEADER('PSO part of spin-spin coupling tensors',-1)
C
               CALL PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                     3*NUCDEP,IATOM,JATOM,SPNPSO,D1)
C
               CALL HEADER(
     &              'SD(+FC) part of spin-spin coupling tensors',-1)
C
               CALL PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                     3*NUCDEP,IATOM,JATOM,SPNSD,D1)
C
               CALL HEADER('FC part of spin-spin coupling tensors',-1)
C
               CALL PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                     3*NUCDEP,IATOM,JATOM,SPNFC,D1)
C
C
            END DO
         END DO
      END IF
C
      CALL QEXIT('LSPSYM')
C
      RETURN
      END
C*********************************************************************
C  /* Deck lspred */
C
      SUBROUTINE LSPRED(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                  JLOCCF,JLOCCL,SPNDSO,SPNPSO,SPNSD,SPNFC,
     &                  SPSDFC,SPNTOT,WORK,LWORK)
C
C      Stephan P. A. Sauer 17-01-2016
C
C      This routine prints the contributions to the indirect nuclear
C      spin-spin coupling constants from individual molecular orbitals
C      after the symmetry adaptation was removed and the proper 
C      g-factors were included.
C
C
#include "implicit.h"
C
C ALPHAC used from include file codata.h
C LUPRI used from include file priunit.h
C XFAMU used from include file codata.h
C PMASS used from include file pmass.h
C
#include "priunit.h"
#include "codata.h"
C#include "pmass.h"
C
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0, D1O3 = 1.0D0/D3, D0 = 0.0D0)
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2/
     &                    (4*XFAMU*XFAMU*PMASS*PMASS)
     &                   *6.5796838999D15)
C
C MXQN, MXAQN used from include file maxaqn.h
C MXCENT, MXCOOR used from include file mxcent.h
C MXCORB used from include file mxorb.h
#include "maxaqn.h"
#include "mxcent.h"
CPFP
C#include "mxorb.h"
#include "maxorb.h"
Cend-PFP
C
C SOLVNT used from COMMON /CBISOL/
C DOPERT used from COMMN /DORPS/
C CHARGE, ISTBNU, NUCDEP, NUCIND used from COMMON /NUCLEI/
C MAXOPR used from COMMON /SYMMET/
C ABUND used from COMMON /SPNOUT/
#include "cbisol.h"
#include "dorps.h"
#include "nuclei.h"
#include "spnout.h"
#include "locinf.h"
#include "symmet.h"
C
      DIMENSION SPNDSO(NLOCCS*NLOCCS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNPSO(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNSD(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPSDFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNTOT(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION WORK(LWORK)
C
      LOGICAL TEST
C
C
      CALL QENTER('LSPRED')
C
C====================================
C     Transform to non-symmetry basis
C====================================
C
      CALL LSPSYM(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &            SPNDSO,SPNPSO,SPNSD,SPNFC,SPSDFC,WORK,LWORK)
C
C----------------------------------------------------------------------
C     Calculate the contributions to the total coupling constant tensor
C----------------------------------------------------------------------
C
      NTOT = NLOCCS*NLOCCS*NLVIRS*NLVIRS*9*NUCDEP*NUCDEP
C
      CALL DZERO (SPNTOT,NTOT)
      IF (SOSOCC) CALL DAXPY (NTOT,D1,SPNDSO,1,SPNTOT,1)
      CALL DAXPY (NTOT,D1,SPNPSO,1,SPNTOT,1)
      CALL DAXPY (NTOT,D1,SPNFC,1,SPNTOT,1)
      CALL DAXPY (NTOT,D1,SPNSD,1,SPNTOT,1)
      CALL DAXPY (NTOT,D1,SPSDFC,1,SPNTOT,1)
C
C-----------------------------------------------------------
C     Calculate the isotropic contributions and print it out
C-----------------------------------------------------------
C
      CALL AROUND
     & ('Response SSCC to be estimated by SOS calculation')
C
      CALL DZERO(SSCCNI,INSPN*NCONT)
      IATOM1 = 0
CPFP
      IATIJC = 0
      INCONT = 0
Cend-PFP
      NUCINS = NUCIND
      IF (SOLVNT) NUCINS = NUCINS - 1
      DO I1 = 1, NUCINS
         DO ISYM1 = 0, MAXOPR
            IF (IAND(ISTBNU(I1),ISYM1) .EQ. 0) THEN
            IATOM1 = IATOM1 + 1
            IF (DOPERT(I1,2)) THEN
            IATOM2 = 0
            DO I2 = 1, I1
               IF (I2 .EQ. I1) THEN
                  MAXSYM = ISYM1 - 1
               ELSE
                  MAXSYM = MAXOPR
               END IF
               DO ISYM2 = 0, MAXSYM
                  IF (IAND(ISTBNU(I2),ISYM2) .EQ. 0) THEN
                  IATOM2 = IATOM2 + 1
                  IF (DOPERT(I2,2)) THEN
CPFP 04-12-2015  IATIJC is usefull as index for the automatic convergence.
C                Moreover the first equation is not used at all here.
C                  IATIJ = IATOM1*(IATOM1 -1 )/2 + IATOM2
                   IATIJC = IATIJC + 1
Cend-PFP
                  NZ1 = NINT(CHARGE(I1))
                  NZ2 = NINT(CHARGE(I2))
C                  CALL TITLER('Indirect spin-spin-coupling between '//
C     &                        NAMDEP(IATOM1)//' and '//
C     &                        NAMDEP(IATOM2)//':','=',-1)
                  TEST = .FALSE.
                  DO ISO1 = 1, 5
                     GVAL1 = DISOTP(NZ1,ISO1,'GVAL')
                     IF (GVAL1 .NE. 0) THEN
                     IF (IATOM1 .EQ. IATOM2) THEN
                        ISOMAX = ISO1
                     ELSE
                        ISOMAX = 5
                     END IF
                     DO ISO2 = 1, ISOMAX
                        GVAL2 = DISOTP(NZ2,ISO2,'GVAL')
                        IF (GVAL2 .NE. 0) THEN
                        ABUND1 = DISOTP(NZ1,ISO1,'ABUNDANCE')
                        ABUND2 = DISOTP(NZ2,ISO2,'ABUNDANCE')
                        IF (    ((ABUND1 .GE. ABUND)
     &                     .AND. (ABUND2 .GE. ABUND))
     &                     .OR.  (.NOT. TEST)) THEN
                          TEST = .TRUE.
                          NA1    = NINT(DISOTP(NZ1,ISO1,'A'))
                          NA2    = NINT(DISOTP(NZ2,ISO2,'A'))
                          FACTOR = AUTOHZ*GVAL1*GVAL2
C
                          INCONT = 1
C
                          CALL PLOPRID(IATIJC,INCONT,NLOCCS,ILOCCF,
     &                                 ILOCCL,JLOCCF,JLOCCL,3*NUCDEP,
     &                                 IATOM1,IATOM2,SPNDSO,FACTOR)
C
                          INCONT = 2
C
                          CALL PLOPRID(IATIJC,INCONT,NLOCCS,ILOCCF,
     &                                 ILOCCL,JLOCCF,JLOCCL,3*NUCDEP,
     &                                 IATOM1,IATOM2,SPNPSO,FACTOR)
C
                          INCONT = 3
C
                          CALL PLOPRID(IATIJC,INCONT,NLOCCS,ILOCCF,
     &                                 ILOCCL,JLOCCF,JLOCCL,3*NUCDEP,
     &                                 IATOM1,IATOM2,SPNSD,FACTOR)
C
                          INCONT = 4
C
                          CALL PLOPRID(IATIJC,INCONT,NLOCCS,ILOCCF,
     &                                 ILOCCL,JLOCCF,JLOCCL,3*NUCDEP,
     &                                 IATOM1,IATOM2,SPNFC,FACTOR)
C
                          INCONT = 5
C
                          CALL PLOPRID(IATIJC,INCONT,NLOCCS,ILOCCF,
     &                                 ILOCCL,JLOCCF,JLOCCL,3*NUCDEP,
     &                                 IATOM1,IATOM2,SPNTOT,FACTOR)
C
                        END IF
                        END IF
C
                     END DO
C
                     END IF
C
                  END DO
C
                  END IF
                  END IF
C
               END DO
            END DO
C
            END IF
            END IF
C
         END DO
      END DO
C
          CALL HEADER('Estimated Spin-spin coupling constant',-1)
C
      DO ITES = 1, IATIJC
        WRITE(LUPRI,'(/,2X,A,I3,5(3X,F8.4),/)')
     &           'Coupling Nº: ',ITES, (SSCCNI(ITES,JTES), JTES = 1, 5)
      END DO
C
      CALL QEXIT('LSPRED')
C
      RETURN
      END
C
C*********************************************************************
C  /* Deck lspres */
C
      SUBROUTINE LSPRES(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                  JLOCCF,JLOCCL,SPNDSO,SPNPSO,SPNSD,SPNFC,
     &                  SPSDFC,SPNTOT,WORK,LWORK)
C
C      Stephan P. A. Sauer 13/11-1997
C
C      This routine prints the contributions to the indirect nuclear
C      spin-spin coupling constants from individual molecular orbitals
C      after the symmetry adaptation was removed and the proper 
C      g-factors were included.
C
C
#include "implicit.h"
C
C ALPHAC used from include file codata.h
C LUPRI used from include file priunit.h
C XFAMU used from include file codata.h
C PMASS used from include file pmass.h
C
#include "priunit.h"
#include "codata.h"
C#include "pmass.h"
C
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0, D1O3 = 1.0D0/D3, D0 = 0.0D0)
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2/
     &                    (4*XFAMU*XFAMU*PMASS*PMASS)
     &                   *6.5796838999D15)
C
C MXQN, MXAQN used from include file maxaqn.h
C MXCENT, MXCOOR used from include file mxcent.h
C MXCORB used from include file mxorb.h
#include "maxaqn.h"
#include "mxcent.h"
CPFP
C#include "mxorb.h"
#include "maxorb.h"
Cend-PFP
C
C SOLVNT used from COMMON /CBISOL/
C DOPERT used from COMMON /DORPS/
C CHARGE, ISTBNU, NUCDEP, NUCIND used from COMMON /NUCLEI/
C MAXOPR used from COMMON /SYMMET/
C ABUND used from COMMON /SPNOUT/
#include "cbisol.h"
#include "dorps.h"
#include "nuclei.h"
#include "spnout.h"
#include "locinf.h"
#include "symmet.h"
C
      DIMENSION SPNDSO(NLOCCS*NLOCCS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNPSO(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNSD(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPSDFC(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION SPNTOT(NLOCCS*NLOCCS*NLVIRS*NLVIRS,3*NUCDEP,3*NUCDEP)
      DIMENSION WORK(LWORK)
C
      LOGICAL TEST
C
C
      CALL QENTER('LSPRES')
C
C====================================
C     Transform to non-symmetry basis
C====================================
C
      CALL LSPSYM(IPRSOS,NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &            SPNDSO,SPNPSO,SPNSD,SPNFC,SPSDFC,WORK,LWORK)
C
C----------------------------------------------------------------------
C     Calculate the contributions to the total coupling constant tensor
C----------------------------------------------------------------------
C
      NTOT = NLOCCS*NLOCCS*NLVIRS*NLVIRS*9*NUCDEP*NUCDEP
C
      CALL DZERO (SPNTOT,NTOT)
      IF (SOSOCC) CALL DAXPY (NTOT,D1,SPNDSO,1,SPNTOT,1)
      CALL DAXPY (NTOT,D1,SPNPSO,1,SPNTOT,1)
      CALL DAXPY (NTOT,D1,SPNFC,1,SPNTOT,1)
      CALL DAXPY (NTOT,D1,SPNSD,1,SPNTOT,1)
      CALL DAXPY (NTOT,D1,SPSDFC,1,SPNTOT,1)
C
C-----------------------------------------------------------
C     Calculate the isotropic contributions and print it out
C-----------------------------------------------------------
C
      CALL AROUND
     & ('Contributions to the indirect spin-spin coupling constant')
C
      IATOM1 = 0
CPFP
      IATIJC = 0
      INCONT = 0
Cend-PFP
      NUCINS = NUCIND
      IF (SOLVNT) NUCINS = NUCINS - 1
      DO I1 = 1, NUCINS
         DO ISYM1 = 0, MAXOPR
            IF (IAND(ISTBNU(I1),ISYM1) .EQ. 0) THEN
            IATOM1 = IATOM1 + 1
            IF (DOPERT(I1,2)) THEN
            IATOM2 = 0
            DO I2 = 1, I1
               IF (I2 .EQ. I1) THEN
                  MAXSYM = ISYM1 - 1
               ELSE
                  MAXSYM = MAXOPR
               END IF
               DO ISYM2 = 0, MAXSYM
                  IF (IAND(ISTBNU(I2),ISYM2) .EQ. 0) THEN
                  IATOM2 = IATOM2 + 1
                  IF (DOPERT(I2,2)) THEN
CPFP 04-12-2015  IATIJC is usefull as index for the automatic convergence.
C                Moreover the first equation is not used at all here.
C                  IATIJ = IATOM1*(IATOM1 -1 )/2 + IATOM2
                   IATIJC = IATIJC + 1
Cend-PFP
                  NZ1 = NINT(CHARGE(I1))
                  NZ2 = NINT(CHARGE(I2))
                  CALL TITLER('Indirect spin-spin-coupling between '//
     &                        NAMDEP(IATOM1)//' and '//
     &                        NAMDEP(IATOM2)//':','=',-1)
                  TEST = .FALSE.
                  DO ISO1 = 1, 5
                     GVAL1 = DISOTP(NZ1,ISO1,'GVAL')
                     IF (GVAL1 .NE. 0) THEN
                     IF (IATOM1 .EQ. IATOM2) THEN
                        ISOMAX = ISO1
                     ELSE
                        ISOMAX = 5
                     END IF
                     DO ISO2 = 1, ISOMAX
                        GVAL2 = DISOTP(NZ2,ISO2,'GVAL')
                        IF (GVAL2 .NE. 0) THEN
                        ABUND1 = DISOTP(NZ1,ISO1,'ABUNDANCE')
                        ABUND2 = DISOTP(NZ2,ISO2,'ABUNDANCE')
                        IF (     ((ABUND1 .GE. ABUND)
     &                      .AND. (ABUND2 .GE. ABUND))
     &                      .OR.  (.NOT. TEST)) THEN
                           TEST = .TRUE.
                           NA1    = NINT(DISOTP(NZ1,ISO1,'A'))
                           NA2    = NINT(DISOTP(NZ2,ISO2,'A'))
                           FACTOR = AUTOHZ*GVAL1*GVAL2
                           WRITE(LUPRI,'(/,2X,A,I3,10X,A,F8.3)')
     &                        'Mass number atom 1: ',NA1,'Abundance: ',
     &                         ABUND1
                           WRITE(LUPRI,'(2X,A,I3,10X,A,F8.3)')
     &                        'Mass number atom 2: ',NA2,'Abundance: ',
     &                         ABUND2
C
                           IF (IPRSOS .GE. 2) THEN
                              IF (SOSOCC) THEN
C
                                 CALL HEADER('DSO part of spin-spin '//
     &                                       'coupling tensor',-1)
C
                                 CALL PLOPRT(NLOCCS,ILOCCF,ILOCCL,
     &                                       JLOCCF,JLOCCL,3*NUCDEP,
     &                                       IATOM1,IATOM2,SPNDSO,
     &                                       FACTOR)
C
C
                                 CALL HEADER('PSO part of spin-spin '//
     &                                       'coupling tensor',-1)
C
                                 CALL PLOPRT(NLOCCS,ILOCCF,ILOCCL,
     &                                       JLOCCF,JLOCCL,3*NUCDEP,
     &                                       IATOM1,IATOM2,SPNPSO,
     &                                       FACTOR)
C
C
                                 CALL HEADER('FC part of spin-spin '//
     &                                       'coupling tensor',-1)
C
                                 CALL PLOPRT(NLOCCS,ILOCCF,ILOCCL,
     &                                       JLOCCF,JLOCCL,3*NUCDEP,
     &                                       IATOM1,IATOM2,SPNFC,
     &                                       FACTOR)
C
C
                                 CALL HEADER('SD part of spin-'//
     &                                       ' spin coupling tensor',-1)
C
                                 CALL PLOPRT(NLOCCS,ILOCCF,ILOCCL,
     &                                       JLOCCF,JLOCCL,3*NUCDEP,
     &                                       IATOM1,IATOM2,SPNSD,
     &                                       FACTOR)
C
C
                              ELSE
C
C
                                 CALL HEADER('DSO part of spin-spin '//
     &                                       'coupling tensor',-1)
C
                                 CALL PLOPRT(NLOCCS,ILOCCF,ILOCCL,
     &                                       JLOCCF,JLOCCL,3*NUCDEP,
     &                                       IATOM1,IATOM2,SPNDSO,
     &                                       FACTOR)
C
C
                                 CALL HEADER('PSO part of spin-spin '//
     &                                       'coupling tensor',-1)
C
                                 CALL PLCPRT(NLOCCS,NLVIRS,
     &                                      ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                                      3*NUCDEP,IATOM1,IATOM2,
     &                                      SPNPSO,FACTOR)
C
C
                                 CALL HEADER('FC part of spin-spin '//
     &                                       'coupling tensor',-1)
C
                                 CALL PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL
     &                                      ,JLOCCF,JLOCCL,3*NUCDEP
     &                                      ,IATOM1,IATOM2,SPNFC,FACTOR)
C
C
                                 CALL HEADER('SD part of spin-'//
     &                                       ' spin coupling tensor',-1)
C
                                 CALL PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL
     &                                      ,JLOCCF,JLOCCL,3*NUCDEP
     &                                      ,IATOM1,IATOM2,SPNSD,FACTOR)
C
C
                              END IF
                           END IF
C
                           IF (SOSOCC) THEN
C
                              INCONT = 1
C
                              CALL HEADER('DSO part of spin-spin '//
     &                                    'coupling constant',-1)
C
                              CALL PLOPRI(NLOCCS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNDSO,FACTOR)
C
                              INCONT = 2
C
                              CALL HEADER('PSO part of spin-spin '//
     &                                    'coupling constant',-1)
C
                              CALL PLOPRI(NLOCCS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNPSO,FACTOR)
C
                              INCONT = 3
C
                              CALL HEADER('SD part of spin-spin'//
     &                                    ' coupling constant',-1)
C
                              CALL PLOPRI(NLOCCS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNSD,FACTOR)
C
                              INCONT = 4
C
                              CALL HEADER('FC part of spin-spin '//
     &                                    'coupling constant',-1)
C
                              CALL PLOPRI(NLOCCS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNFC,FACTOR)
C
                              INCONT = 5
C
                              CALL HEADER('Total spin-spin coupling '//
     &                                    'constant',-1)
C
                              CALL PLOPRI(NLOCCS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNTOT,FACTOR)
C
C
                           ELSE
C
                              INCONT = 1
C
                              CALL HEADER('DSO part of spin-spin '//
     &                                    'coupling constant',-1)
C
                              CALL PLOPRI(NLOCCS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNDSO,FACTOR)
C
                              INCONT = 2
C
                              CALL HEADER('PSO part of spin-spin '//
     &                                    'coupling constant',-1)
C
                              CALL PLCPRI(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNPSO,FACTOR)
C
                              INCONT = 3
C
                              CALL HEADER('SD part of spin-spin'//
     &                                    ' coupling constant',-1)
C
                              CALL PLCPRI(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNSD,FACTOR)
C
                              INCONT = 4
C
                              CALL HEADER('FC part of spin-spin '//
     &                                    'coupling constant',-1)
C
                              CALL PLCPRI(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNFC,FACTOR)
C
                              INCONT = 5
C
                              CALL HEADER('Total spin-spin '//
     &                                    'coupling constant',-1)
C
                              CALL PLCPRI(NLOCCS,NLVIRS,ILOCCF,ILOCCL,
     &                                    JLOCCF,JLOCCL,3*NUCDEP,
     &                                    IATOM1,IATOM2,SPNTOT,FACTOR)
C
                           END IF
C
C
                        END IF
                        END IF
C
                     END DO
C
                     END IF
C
                  END DO
C
                  END IF
                  END IF
C
               END DO
            END DO
C
            END IF
            END IF
C
         END DO
      END DO
C
      CALL QEXIT('LSPRES')
C
      RETURN
      END

C***********************************************************************
C  /* Deck plocpr */
C
      SUBROUTINE PLOCPR(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,PROP)
C
C
C      Stephan P. A. Sauer 12/11-1997
C
C      This routine prints the contributions from pairs of orbitals to
C      a property.
C
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C TABOCL(i), TABVIL(i) used from locinf.h
#include "priunit.h"
#include "mxcent.h"
#include "spnout.h"
C#include "inflin.h"
#include "locinf.h"
C
      PARAMETER (D0 = 0.0D+00, THRPR = 1.0D-06)
C
      DIMENSION PROP(NLOCCS*NLOCCS*NLVIRS*NLVIRS)
C
      CALL QENTER('PLOCPR')
C
      WRITE (LUPRI,'(/,5X,A,A)') '    label    A  label    I',
     &                             '  label    B  label    J'
C
      IJLVAR = 0
      TOT    = D0
C
      DO JLOCC = JLOCCF,JLOCCL
C
         DO JLVIR = 1, NLVIRS
C
            DO ILOCC = ILOCCF,ILOCCL
C
               DO ILVIR = 1, NLVIRS
C
                  JLVIRF = (JLVIR - 1)*NLVIRS*NLOCCS
                  IJLVAR =  JLVIRF + ILVIR
C
                  IF (DABS(PROP(IJLVAR)) .GT. THRPR)
     &                WRITE (LUPRI,'(1X,I6,1X,4(A8,I3,1X),F12.4)')
     &                       IJLVAR,TABVIL(ILVIR),ILVIR,
     &                              TABOCL(ILOCC),ILOCC,
     &                              TABVIL(JLVIR),JLVIR,
     &                              TABOCL(JLOCC),JLOCC,PROP(IJLVAR)
                  TOT = TOT + PROP(IJLVAR)
               END DO
            END DO
         END DO
      END DO
      WRITE (LUPRI,'(10X,A,41X,D10.3,2X)') 'Total ',TOT
C
      CALL QEXIT('PLOCPR')
C
      RETURN
C
      END
C
C********************************************************************
C  /* Deck ploprid */
C
      SUBROUTINE PLOPRID(IATIJD,INCOND,NLOCCS,ILOCCF,ILOCCL,
     &                   JLOCCF,JLOCCL,NCOOR,IATOM,JATOM,PROP,FACTOR)
C
C      Stephan P. A. Sauer 17-01-2016
C
C      This routine prints the contributions from pairs of occupied
C      orbitals to the isotropic part of a property
C      and scales it with a factor.
C
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C TABOCL(i), TABVIL(i) used from locinf.h
#include "priunit.h"
#include "mxcent.h"
#include "spnout.h"
C#include "inflin.h"
#include "locinf.h"
C
      PARAMETER (D0 = 0.0D+00, D3 = 3.0D0, D1O3 = 1.0D0/D3)
      PARAMETER (THRPR = 1.0D-02)
C
      DIMENSION PROP(NLOCCS*NLOCCS,NCOOR,NCOOR)
C
      CALL QENTER('PLOPRID')
C
      IJLVAR = 0
      TOT    = D0
C
      DO JLOCC = JLOCCF,JLOCCL
         DO ILOCC = ILOCCF,ILOCCL
C
            JLOOCF = (JLOCC - 1)*NLOCCS
            IJLVAR = JLOOCF + ILOCC
C
            IC = 3*(IATOM-1)
            JC = 3*(JATOM-1)
C
            AVEISO = ( PROP(IJLVAR,IC+1,JC+1)
     &               + PROP(IJLVAR,IC+2,JC+2)
     &               + PROP(IJLVAR,IC+3,JC+3))
     &             * FACTOR * D1O3
C
            TOT = TOT + AVEISO
C
         END DO
      END DO
C
      SSCCNI(IATIJD,INCOND) = TOT
C
      CALL QEXIT('PLOPRID')
C
      RETURN
C
      END
C*********************************************************************
C  /* Deck plcpri */
C
      SUBROUTINE PLCPRI(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                  NCOOR,IATOM,JATOM,PROP,FACTOR)
C
C      Stephan P. A. Sauer 13/11-1997
C
C      This routine prints the contributions from pairs of orbitals to
C      the isotropic part of a property and scales it with a factor.
C
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C TABOCL(i), TABVIL(i) used from locinf.h
#include "priunit.h"
CPFP added to use NVARPT in locinf.h
#include "mxcent.h"
#include "spnout.h"
C#include "inflin.h"
#include "locinf.h"
C
      PARAMETER (D0 = 0.0D+00, D3 = 3.0D0, D1O3 = 1.0D0/D3)
      PARAMETER (THRPR = 1.0D-02)
C
      DIMENSION PROP(NLOCCS*NLOCCS*NLVIRS*NLVIRS,NCOOR,NCOOR)
C
      CALL QENTER('PLCPRI')
C
      WRITE (LUPRI,'(/,5X,A,A)') '    label    A  label    I',
     &                             '  label    B  label    J '
C
      IJLVAR = 0
      TOT    = D0
C
      DO JLOCC = JLOCCF,JLOCCL
C
         DO JLVIR = 1, NLVIRS
C
            DO ILOCC = ILOCCF,ILOCCL
C
               DO ILVIR = 1, NLVIRS
C
                  JLVIRF = (JLVIR - 1)*NLVIRS*NLOCCS
                  IJLVAR =  JLVIRF + ILVIR
C
                  IC = 3*(IATOM-1)
                  JC = 3*(JATOM-1)
C
                  AVEISO = ( PROP(IJLVAR,IC+1,JC+1)
     &                     + PROP(IJLVAR,IC+2,JC+2)
     &                     + PROP(IJLVAR,IC+3,JC+3))
     &                   * FACTOR * D1O3
C
                  TOT = TOT + AVEISO
C
                  IF (DABS(AVEISO) .GT. THRPR)
     &                WRITE (LUPRI,'(1X,I6,1X,4(A8,I3,1X),F12.4)')
     &                       IJLVAR,TABVIL(ILVIR),ILVIR,
     &                              TABOCL(ILOCC),ILOCC,
     &                              TABVIL(JLVIR),JLVIR,
     &                              TABOCL(JLOCC),JLOCC,AVEISO
               END DO
            END DO
         END DO
      END DO
C
      WRITE (LUPRI,'(10X,A,41X,F12.4)') 'TOTAL',TOT
CPFP 04-12-2015 for the automatic convergence
C     Keep the wright SSCC with the wright index & no value for Singularities
C          WRITE (LUPRI,*) ' IN PLCPRI :',IATIJC, NCSSCC, NCMTHD, THRSLO
      IF (SOSOSC) THEN
        IF (ISINGUL .EQ. 0) THEN
          IF ((IATIJC .EQ. NCSSCC) .AND. (INCONT .EQ. NCTERM)) THEN
            SSCCPE(INDEXT,NCTERM) = TOT
            WRITE (LUPRI,'(/,1X,A)')'Automatic convergence of the SOS'
            WRITE (LUPRI,'(1X,A)')  '--------------------------------'
            WRITE (LUPRI,'(/,1X,3(A,I4),A,F4.2)')
     &           ' SSCC Nº:', NCSSCC,'  Term Nº: ', NCTERM,
     &            '  Method:', NCMTHD,'  Threshold [%]:', THRSLO
            WRITE (LUPRI,'(1X,A,A)')
     &            '-----------------------------------',
     &            '-----------------------------------'
            WRITE (LUPRI,'(1X,A,I4,2(2X,A,F8.4))')
     &            'Pseudo Exc. Energ. Nº:', INDEXT,
     &            ' App. value :', SSCCPE(INDEXT,INCONT),
     &            ' Resp. value:', SSCCNI(IATIJC,NCTERM)
            WRITE (LUPRI,'(1X,A,A)')
     &            '-----------------------------------',
     &            '-----------------------------------'
          END IF
        ELSE
          SSCCPE(INDEXT,NCTERM) = 0.0D0
          WRITE (LUPRI,'(/,1X,A)')
     &           ' The program skips this value due to the singularity'
        END IF
      END IF
Cend-PFP
C
      CALL QEXIT('PLCPRI')
C
      RETURN
C
      END
C********************************************************************
C  /* Deck plopri */
C
      SUBROUTINE PLOPRI(NLOCCS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                  NCOOR,IATOM,JATOM,PROP,FACTOR)
C
C      Stephan P. A. Sauer 15/3-1999
C
C      This routine prints the contributions from pairs of occupied
C      orbitals to the isotropic part of a property
C      and scales it with a factor.
C
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C TABOCL(i), TABVIL(i) used from locinf.h
#include "priunit.h"
#include "mxcent.h"
#include "spnout.h"
C#include "inflin.h"
#include "locinf.h"
C
      PARAMETER (D0 = 0.0D+00, D3 = 3.0D0, D1O3 = 1.0D0/D3)
      PARAMETER (THRPR = 1.0D-02)
C
      DIMENSION PROP(NLOCCS*NLOCCS,NCOOR,NCOOR)
C
      CALL QENTER('PLOPRI')
C
      WRITE (LUPRI,'(/,5X,A)') '           label     I  label     J'
C
      IJLVAR = 0
      TOT    = D0
C
      DO JLOCC = JLOCCF,JLOCCL
         DO ILOCC = ILOCCF,ILOCCL
C
            JLOOCF = (JLOCC - 1)*NLOCCS
            IJLVAR = JLOOCF + ILOCC
C
            IC = 3*(IATOM-1)
            JC = 3*(JATOM-1)
C
            AVEISO = ( PROP(IJLVAR,IC+1,JC+1)
     &               + PROP(IJLVAR,IC+2,JC+2)
     &               + PROP(IJLVAR,IC+3,JC+3))
     &             * FACTOR * D1O3
C 
            TOT = TOT + AVEISO
C
            IF (DABS(AVEISO) .GT. THRPR)
     &           WRITE (LUPRI,'(1X,I6,1X,2(2X,A8,I3),F12.4)')
     &                  IJLVAR,TABOCL(ILOCC),ILOCC,
     &                         TABOCL(JLOCC),JLOCC,AVEISO
         END DO
      END DO
C
      WRITE (LUPRI,'(10X,A,20X,F12.4)') 'TOTAL',TOT
CPFP 04-12-2015 for the automatic convergence
C     Keep the wright SSCC with the wright index & no value for Singularities
C          WRITE (LUPRI,*) ' IN PLCPRI :',IATIJC, NCSSCC, NCMTHD, THRSLO
      IF (SOSOSC) THEN
        IF (ISINGUL .EQ. 0) THEN
          IF ((IATIJC .EQ. NCSSCC) .AND. (INCONT .EQ. NCTERM)) THEN
            SSCCPE(INDEXT,NCTERM) = TOT
            WRITE (LUPRI,'(/,1X,A)')'Automatic convergence of the SOS'
            WRITE (LUPRI,'(1X,A)')  '--------------------------------'
            WRITE (LUPRI,'(/,1X,3(A,I4),A,F4.2)')
     &            ' SSCC Nº:', NCSSCC,'  Term Nº: ', NCTERM,
     &            '  Method:', NCMTHD,'  Threshold [%]:', THRSLO
            WRITE (LUPRI,'(1X,A,A)')
     &            '-----------------------------------',
     &            '-----------------------------------'
            WRITE (LUPRI,'(1X,A,I4,2(2X,A,F8.4))')
     &            'Pseudo Exc. Energ. Nº:', INDEXT,
     &            ' App. value :', SSCCPE(INDEXT,INCONT),
     &            ' Resp. value:', SSCCNI(IATIJC,NCTERM)
            WRITE (LUPRI,'(1X,A,A)')
     &            '-----------------------------------',
     &            '-----------------------------------'
          END IF
        ELSE
          SSCCPE(INDEXT,NCTERM) = 0.0D0
          WRITE (LUPRI,'(/,1X,A)')
     &           ' The program skips this value due to the singularity'
        END IF
      END IF
Cend-PFP
C
      CALL QEXIT('PLOPRI')
C
      RETURN
C
      END

C*********************************************************************
C  /* Deck plcprt */
C
      SUBROUTINE PLCPRT(NLOCCS,NLVIRS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                  NCOOR,IATOM,JATOM,PROP,FACTOR)
C
C      Stephan P. A. Sauer 13/11-1997
C
C      This routine prints the contributions from pairs of orbitals to
C      a property tensor
C
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C TABOCL(i), TABVIL(i) used from locinf.h
#include "priunit.h"
#include "locinf.h"
C
      PARAMETER (D0 = 0.0D+00)
      PARAMETER (THRPR = 1.0D-02)
C
      DIMENSION PROP(NLOCCS*NLOCCS*NLVIRS*NLVIRS,NCOOR,NCOOR)
      DIMENSION TJ(3,3)
C
      CALL QENTER('PLCPRT')
C
      WRITE (LUPRI,'(/,5X,2A,9(5X,A,5X))') '    label    A  label    I',
     &                                       '  label    B  label    J',
     &                     'xx','xy','xz','yx','yy','yz','zx','zy','zz'
C
      IJLVAR = 0
      CALL DZERO(TJ,9)
C
      DO JLOCC = JLOCCF, JLOCCL
C
         DO JLVIR = 1, NLVIRS
C
            DO ILOCC = ILOCCF, ILOCCL
C
               DO ILVIR = 1, NLVIRS
C
                  JLVIRF = (JLVIR - 1)*NLVIRS*NLOCCS
                  IJLVAR =  JLVIRF + ILVIR
C
                  IC = 3*(IATOM-1)
                  JC = 3*(JATOM-1)
C
                  TOT = D0
                  DO I = 1, 3
                     DO J = 1, 3
                        TOT     = TOT + DABS(PROP(IJLVAR,IC+I,JC+J))
                        TJ(I,J) = TJ(I,J) + PROP(IJLVAR,IC+I,JC+J)
                     END DO
                  END DO
C
                  IF (DABS(TOT*FACTOR) .GT. THRPR)
     &                WRITE (LUPRI,'(1X,I6,1X,4(A8,I3,1X),9(D10.3,2X))')
     &                       IJLVAR,TABVIL(ILVIR),ILVIR,
     &                              TABOCL(ILOCC),ILOCC,
     &                              TABVIL(JLVIR),JLVIR,
     &                              TABOCL(JLOCC),JLOCC,
     &                     ((PROP(IJLVAR,IC+I,JC+J)*FACTOR,J=1,3),I=1,3)
C
               END DO
            END DO
         END DO
      END DO
C
      WRITE (LUPRI,'(7X,A,7X,9(D10.3,2X))') 'TOTAL',
     &                                    ((TJ(I,J)*FACTOR,J=1,3),I=1,3)
C
      CALL QEXIT('PLCPRT')
C
      RETURN
C
      END
C*******************************************************************
C  /* Deck ploprt */
C
      SUBROUTINE PLOPRT(NLOCCS,ILOCCF,ILOCCL,JLOCCF,JLOCCL,
     &                  NCOOR,IATOM,JATOM,PROP,FACTOR)
C
C      Stephan P. A. Sauer 16/3-1999
C
C      This routine prints the contributions from pairs of occupied
C      orbitals to a property tensor
C
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C TABOCL(i), TABVIL(i) used from locinf.h
#include "priunit.h"
#include "locinf.h"
C
      PARAMETER (D0 = 0.0D+00)
      PARAMETER (THRPR = 1.0D-02)
C
      DIMENSION PROP(NLOCCS*NLOCCS,NCOOR,NCOOR)
      DIMENSION TJ(3,3)
C
      CALL QENTER('PLOPRT')
C
      WRITE (LUPRI,'(/,5X,2A,9(5X,A,5X))') '    label    I  label    J',
     &                     'xx','xy','xz','yx','yy','yz','zx','zy','zz'
C
      IJLVAR = 0
      CALL DZERO(TJ,9)
C
      DO JLOCC = JLOCCF,JLOCCL
         DO ILOCC = ILOCCF,ILOCCL
C
            JLOOCF = (JLOCC - 1)*NLOCCS
            IJLVAR = JLOOCF + ILOCC
C
            IC = 3*(IATOM-1)
            JC = 3*(JATOM-1)
C
            TOT = D0
            DO I = 1, 3
               DO J = 1, 3
                  TOT     = TOT + DABS(PROP(IJLVAR,IC+I,JC+J))
                  TJ(I,J) = TJ(I,J) + PROP(IJLVAR,IC+I,JC+J)
               END DO
            END DO
C
            IF (DABS(TOT*FACTOR) .GT. THRPR)
     &           WRITE (LUPRI,'(1X,I6,1X,2(A8,I3,1X),9(D10.3,2X))')
     &                 IJLVAR,TABOCL(ILOCC),ILOCC,
     &                        TABOCL(JLOCC),JLOCC,
     &                 ((PROP(IJLVAR,IC+I,JC+J)*FACTOR,J=1,3),I=1,3)
C
         END DO
      END DO
C
      WRITE (LUPRI,'(7X,A,7X,9(D10.3,2X))') 'TOTAL',
     &                                    ((TJ(I,J)*FACTOR,J=1,3),I=1,3)
C
      CALL QEXIT('PLOPRT')
C
      RETURN
C
      END
C*********************************************************************
C  /* Deck sercon */
C
      SUBROUTINE SERCON(NOCCTA,NOMTHD,NEXCII,NEXCIF,NJUMPZ,
     &                  NOTHRS,RESSLT,WORK,LWORK)
C
C
C      Stephan P. A. Sauer 13/12-2015
C
C      This routine is the driver for the automatic convergence of the
C      to the SOS indirect nuclear spin-spin coupling constants from
C      individual molecular orbitals.
C
#include "implicit.h"
C
C MXCOOR used from include file mxcent.h
C LUPRI used from include file priunit.h
#include "mxcent.h"
#include "priunit.h"
C
C NUCDEP used from COMMON /NUCLEI/
C NSYM, NORBT, NOCCT, NVIRT, NCMOT used from COMMON /INFORB/
C ISPPRI used from COMMON /SPNOUT/
C NVARPT used from COMMON /INFLIN/
#include "nuclei.h"
#include "inforb.h"
#include "spnout.h"
#include "inflin.h"
C   used from COMMON /LOCINF/
#include "maxorb.h"
#include "locinf.h"
C
      INTEGER IV(100)
      REAL*8 CHISQ, ALAMBDA
      DIMENSION WORK(LWORK)
      DIMENSION ACOEF(3),DYESTI(3),ACOVAR(3,3),ALPHA(3,3),NA(3),V(5000)
      CALL QENTER('SERCON')
C
C   If there is a "singularity in reduced response equation"
C   then avoid any further calculation
C
      IF (ISINGUL .EQ. 0) THEN
C
        INDEX1 = INT((INDEXT - NEXCII)/NJUMPZ + 1)
C
        KWORK1 = 1
        KWORK2 = KWORK1 + INDEX1
        KWORK3 = KWORK2 + INDEX1
        KWORK4 = KWORK3 + INDEX1
        KWORK5 = KWORK4 + INDEX1
        LWORK5 = LWORK  - KWORK5
C
        IF (LWORK7 .LT. 0) THEN
          WRITE(LUPRI,*)'SERCON.1: Need :',KWORK5,
     &                '  Available :',LWORK5
          CALL QUIT('Insufficient memory in SOSDRV.1')
        END IF
C
        CALL DZERO(WORK(KWORK1),INDEX1)
        CALL DZERO(WORK(KWORK2),INDEX1)
        CALL DZERO(WORK(KWORK3),INDEX1)
        CALL DZERO(WORK(KWORK4),INDEX1)
        CALL DZERO(WORK(KWORK5),INDEX1)
C
        DESTER  = 0.0D0
        XESTIM  = 0.0D0
        YESTIM  = 0.0D0
        XESTIM1 = 0.0D0
        YESTIM1 = 0.0D0
        XESTIM2 = 0.0D0
        YESTIM2 = 0.0D0
        NCOUNT  = 0
        NNEXCF  = NEXCIF - NJUMPZ
C
        DO IEE = NEXCII,INDEXT,NJUMPZ
           IF (SSCCPE(IEE,NCTERM) .NE. 0.0D0) THEN
             NCOUNT = NCOUNT + 1
             WORK(KWORK1-1+NCOUNT) = DFLOAT(NCOUNT)
             WORK(KWORK2-1+NCOUNT) = DFLOAT(IEE)
             WORK(KWORK3-1+NCOUNT) = SSCCPE(IEE,NCTERM)
             WORK(KWORK4-1+NCOUNT) = 0.01D0*SSCCPE(IEE,NCTERM)
           END IF
        END DO
C
C      CALL DCOPY(INDEX1,WORK(KWORK1),1,SSCCPE(INDEX1,5),1)
        WRITE(LUPRI,'(/,1X,A,2I4))')'Series to be analyzed :',
     &                              INDEX1+1,INDEXT
        WRITE(LUPRI,'(/,A)')' In SERCON :    Step,     X,       Y'
        DO I = 1, NCOUNT
           WRITE(LUPRI,'(12X,3(F9.4))') WORK(KWORK1-1+I),
     &          WORK(KWORK2-1+I),WORK(KWORK3-1+I)
        END DO
C
        ILAST = NCOUNT
        IPREV = NCOUNT - 1
C
C--------------------------------------------------
C     Method 1: analyze E = |J(i) - J(i-1)| / J(i)
C--------------------------------------------------
C
        IF (NCMTHD .EQ. 1) THEN
C
          PNL = WORK(KWORK3-1+ILAST)
          IF (IPREV .LE. 0) THEN
            PNP = 0.0D0
          ELSE
            PNP = WORK(KWORK3-1+IPREV)
          END IF
C
          RESSLT = 100.0D0*DABS((PNL - PNP)/PNL)
C
          WRITE(LUPRI,'(/,1X,2(A,F6.2))')
     &              '1- The Er = |J(i) - J(i-1)| / J(i) gives :',
     &                  RESSLT,'  Required :',THRSLO
          WRITE(LUPRI,'(/,1X,2(A,F6.2))')' Resp. value :', PNL,
     &                                 '  App. value :', PNP
C
C--------------------------------------------------------
C     Method 2: analyze E = |J(final) - J(i)| / J(final)
C--------------------------------------------------------
C
        ELSE IF (NCMTHD .EQ. 2) THEN
C
          PNL = SSCCNI(NCSSCC,NCTERM)
          IF (IPREV .LE. 0) THEN
            PNP = 0.0D0
          ELSE
            PNP = WORK(KWORK3-1+ILAST)
          END IF
C
          RESSLT = 100.0D0*DABS((PNL - PNP)/PNL)
C
          WRITE(LUPRI,'(/,1X,2(A,F6.2))')
     &          '2- The Er = |J(final) - J(i)| / J(final)  gives :',
     &                  RESSLT,'  Required :',THRSLO
          WRITE(LUPRI,'(/,1X,2(A,F6.2))')' Resp. value :', PNL,
     &                                   '  App. value :', PNP
C
C--------------------------------------------------
C     Method 3: combine together Er's from 1 and 2
C--------------------------------------------------
C
        ELSE IF (NCMTHD .EQ. 3) THEN
C
          PNL1 = WORK(KWORK3-1+ILAST)
          PNL2 = SSCCNI(NCSSCC,NCTERM)
          IF (IPREV .LE. 0) THEN
            PNP1 = 0.0D0
            PNP2 = 0.0D0
          ELSE
            PNP1 = WORK(KWORK3-1+IPREV)
            PNP2 = WORK(KWORK3-1+ILAST)
          END IF
C
          RESSLT1 = 100.0D0*DABS((PNL1 - PNP1)/PNL1)
          RESSLT2 = 100.0D0*DABS((PNL2 - PNP2)/PNL2)
C
          RESSLT = MAX(RESSLT1,RESSLT2)
C
          WRITE(LUPRI,'(/,1X,2(A,F10.2),A,F10.2)')' 3- The Ers give :',
     &                  RESSLT1,' &',RESSLT2,'  Required :',THRSLO
          WRITE(LUPRI,'(/,1X,2(A,F10.2,A,F10.2))')
     &                  ' Resp. value :', PNL1,'  &', PNL2,
     &                  '  App. value :', PNP1,'  &', PNP2
          WRITE(LUPRI,'(/,1X,A)')
     &          'The largest is used for comparizon with the required'
C
C----------------------------------------------------------------------------
C     Method 4: Using polinomial fitting we
C               get an i value and analyze E = |J(f) - J(i)| / J(f)
C               infere the i=f [symb. i_(f)] value and analyze 
C               E = |J(f) - J(i_(f))| / J(f)
C----------------------------------------------------------------------------
C
C     Watch out the following three methods.
C     It would be better if one only uses the last 20 approximatelly for the vectors X & Y ... !!!
C
        ELSE IF (NCMTHD .EQ. 4) THEN
C
          YESTIM1 = 0.0D0
          YESTIM2 = 0.0D0
C
C          PNL = WORK(KWORK3-1+ILAST)
          PNL = SSCCNI(NCSSCC,NCTERM)
C
          XESTIM1 = WORK(KWORK2-1+NCOUNT) + DFLOAT(NJUMPZ)
          IF (XESTIM1 .GT. NVARPT) XESTIM1 = WORK(KWORK2-1+NCOUNT)
          XESTIM2 = DFLOAT(NVARPT)
C
          IF (ILAST .GE. 2)  THEN
C
            NLAST5 = MIN(NCOUNT,5)
C
c            CALL polint(WORK(KWORK2),WORK(KWORK3),
c     &                  NCOUNT,XESTIM1,YESTIM1,DESTER1,LUPRI)
            CALL polint(WORK(KWORK2+NCOUNT-NLAST5),
     &                  WORK(KWORK3+NCOUNT-NLAST5),
     &                  NLAST5,XESTIM1,YESTIM1,DESTER1,LUPRI)
C
            WRITE(LUPRI,'(/,1X,A,3(F14.2))')
     &           ' The polint1 method gives:',XESTIM1,YESTIM1,DESTER1
C
c            CALL polint(WORK(KWORK2),WORK(KWORK3),
c     &                  NCOUNT,XESTIM2,YESTIM2,DESTER2,LUPRI)
            CALL polint(WORK(KWORK2+NCOUNT-NLAST5),
     &                  WORK(KWORK3+NCOUNT-NLAST5),
     &                  NLAST5,XESTIM2,YESTIM2,DESTER2,LUPRI)
C
            WRITE(LUPRI,'(/,1X,A,3(F14.2))')
     &            ' The polint2 method gives:',XESTIM2,YESTIM2,DESTER2
C
          END IF
C
          PNP1 = YESTIM1
          PNP2 = YESTIM2
          RESSLT1 = 100.0D0*DABS((PNP1 - PNL)/PNL)
          RESSLT2 = 100.0D0*DABS((PNL - PNP2)/PNL)
C
          RESSLT = MAX(RESSLT1,RESSLT2)
C
          WRITE(LUPRI,'(/,1X,2(A,F14.2),A,F14.2)')' 4- The Ers give :',
     &                  RESSLT1,' &',RESSLT2,'  Required :',THRSLO
          WRITE(LUPRI,'(/,1X,3(A,F10.2))')
     &                ' Resp. value :', PNL,
     &                  '  App. values :', PNP1,'  &', PNP2
          WRITE(LUPRI,'(/,1X,A)')
     &   'The largest error is used for comparizon with the required!!!'
C
C----------------------------------------------------------------------------
C     Method 5: Using rational function fitting we
C               get an i value and analyze E = |J(f) - J(i)| / J(f)
C               get the i_(f) value and analyze E = |J(f) - J(i_(f))| / J(f)
C----------------------------------------------------------------------------
C
        ELSE IF (NCMTHD .EQ. 5) THEN
C
          YESTIM1 = 0.0D0
          YESTIM2 = 0.0D0
C          PNL = WORK(KWORK3-1+ILAST)
          PNL = SSCCNI(NCSSCC,NCTERM)
C
          XESTIM1 = WORK(KWORK2-1+NCOUNT) + DFLOAT(NJUMPZ)
          IF (XESTIM1 .GT. NVARPT) XESTIM1 = WORK(KWORK2-1+NCOUNT)
          XESTIM2 = DFLOAT(NVARPT)
C
          IF (ILAST .GE. 2)  THEN
C
            NLAST5 = MIN(NCOUNT,10)
C
c            CALL ratint(WORK(KWORK2),WORK(KWORK3),
c     &                  NCOUNT,XESTIM1,YESTIM1,DESTER1,LUPRI)
            CALL ratint(WORK(KWORK2+NCOUNT-NLAST5),
     &                  WORK(KWORK3+NCOUNT-NLAST5),
     &                  NLAST5,XESTIM1,YESTIM1,DESTER1,LUPRI)
C
            WRITE(LUPRI,'(/,1X,A,3(F14.2))')
     &            ' The ratint1 method gives:',XESTIM1,YESTIM1,DESTER1
C
c            CALL ratint(WORK(KWORK2),WORK(KWORK3),
c     &                  NCOUNT,XESTIM2,YESTIM2,DESTER2,LUPRI)
            CALL ratint(WORK(KWORK2+NCOUNT-NLAST5),
     &                  WORK(KWORK3+NCOUNT-NLAST5),
     &                  NLAST5,XESTIM2,YESTIM2,DESTER2,LUPRI)
C
            WRITE(LUPRI,'(/,1X,A,3(F14.2))')
     &            ' The ratint2 method gives:',XESTIM2,YESTIM2,DESTER2
C
          END IF
C
          PNP1 = YESTIM1
          PNP2 = YESTIM2
          RESSLT1 = 100.0D0*DABS((PNP1 - PNL)/PNL)
          RESSLT2 = 100.0D0*DABS((PNL - PNP2)/PNL)
C
          RESSLT = MAX(RESSLT1,RESSLT2)
C
          WRITE(LUPRI,'(/,1X,2(A,F14.2),A,F14.2)')
     &            ' 5- The Ers give :',RESSLT1,
     &            ' &',RESSLT2,'  Required :',THRSLO
          WRITE(LUPRI,'(/,1X,3(A,F10.2))')
     &                  ' Resp. value :', PNL,
     &                  '  App. values :', PNP1,'  &', PNP2
          WRITE(LUPRI,'(/,1X,A)')
     &   'The largest error is used for comparizon with the required!!!'
C
C----------------------------------------------------------------------------
C     Method 6: Using spline function fitting we
C               get an i value and analyze E = |J(f) - J(i)| / J(f)
C               get the i_(f) value and analyze E = |J(f) - J(i_(f))| / J(f)
C----------------------------------------------------------------------------
C
        ELSE IF (NCMTHD .EQ. 6) THEN
C
          YESTIM1 = 0.0D0
          YESTIM2 = 0.0D0
C          PNL = WORK(KWORK3-1+ILAST)
          PNL = SSCCNI(NCSSCC,NCTERM)
C
          XESTIM1 = WORK(KWORK2-1+NCOUNT) + DFLOAT(NJUMPZ)
          IF (XESTIM1 .GT. NVARPT) XESTIM1 = WORK(KWORK2-1+NCOUNT)
          XESTIM2 = DFLOAT(NVARPT)
C
          IF (ILAST .GE. 10)  THEN
C
            NLAST5 = MIN(NCOUNT,100)
C
c            CALL spline(WORK(KWORK2),WORK(KWORK3),
c     &                  NCOUNT,XESTIM1,YESTIM1,WORK(KWORK4))
c            CALL splint(WORK(KWORK2),WORK(KWORK3),WORK(KWORK4),
c     &                 NCOUNT,XESTIM1,YESTIM1)
            CALL spline(WORK(KWORK2+NCOUNT-NLAST5),
     &                  WORK(KWORK3+NCOUNT-NLAST5),
     &                  NLAST5,XESTIM1,YESTIM1,
     &                  WORK(KWORK4+NCOUNT-NLAST5))
            CALL splint(WORK(KWORK2+NCOUNT-NLAST5),
     &                  WORK(KWORK3+NCOUNT-NLAST5),
     &                  WORK(KWORK4+NCOUNT-NLAST5),
     &                  NLAST5,XESTIM1,YESTIM1)
C
            WRITE(LUPRI,'(/,1X,A,2(F14.2))')
     &          ' The spline1 method gives:',XESTIM1,YESTIM1
C
c            CALL spline(WORK(KWORK2),WORK(KWORK3),
c     &                  NCOUNT,XESTIM2,YESTIM2,WORK(KWORK4))
c            CALL splint(WORK(KWORK2),WORK(KWORK3),WORK(KWORK4),
c     &                  NCOUNT,XESTIM2,YESTIM2)
            CALL spline(WORK(KWORK2+NCOUNT-NLAST5),
     &                  WORK(KWORK3+NCOUNT-NLAST5),
     &                  NLAST5,XESTIM2,YESTIM2,
     &                  WORK(KWORK4+NCOUNT-NLAST5))
            CALL splint(WORK(KWORK2+NCOUNT-NLAST5),
     &                  WORK(KWORK3+NCOUNT-NLAST5),
     &                  WORK(KWORK4+NCOUNT-NLAST5),
     &                  NLAST5,XESTIM2,YESTIM2)
C
            WRITE(LUPRI,'(/,1X,A,2(F14.2))')
     &           ' The spline2 method gives:',XESTIM2,YESTIM2
C
          END IF
C
          PNP1 = YESTIM1
          PNP2 = YESTIM2
          RESSLT1 = 100.0D0*DABS((PNP1 - PNL)/PNL)
          RESSLT2 = 100.0D0*DABS((PNL - PNP2)/PNL)
C
          RESSLT = MAX(RESSLT1,RESSLT2)
C
          WRITE(LUPRI,'(/,1X,2(A,F14.2),A,F14.2)')
     &          ' 6- The Ers give :',RESSLT1,
     &          ' &',RESSLT2,'  Required :',THRSLO
          WRITE(LUPRI,'(/,1X,3(A,F14.2))')
     &                  ' Resp. value :', PNL,
     &                  '  App. values :', PNP1,'  &', PNP2
          WRITE(LUPRI,'(/,1X,A,/)')
     &   'The largest error is used for comparizon with the required!!!'
C
C----------------------------------------------------------------------------
C     Method 7: Using the Y = A + B*Tanh(C*N) function fitting we
C               get an i value and analyze E = |J(f) - J(i)| / J(f)
C               get the i_(f) value and analyze E = |J(f) - J(i_(f))| / J(f)
C               Here is used the Numerical Recipes'
C----------------------------------------------------------------------------
C
        ELSE IF (NCMTHD .EQ. 7) THEN
C
          CHISQ =  0.0D0
          ALAMBDA = -1.0D0
          THRS = THRSLO*0.01D0
C
          DO I = 1, 3
            NA(I) = 1
          END DO
C
c            ACOEF(1) = 0.1D0*SSCCNI(NCSSCC,NCTERM)
c            ACOEF(2) = -2.0D0*SSCCNI(NCSSCC,NCTERM)
c            ACOEF(3) = ABS(0.001D0*SSCCNI(NCSSCC,NCTERM))
C
            ACOEF(1) = 1.0D0
            ACOEF(2) = 1.0D0
            ACOEF(3) = 1.0D0
C
          CALL DZERO(ACOVAR,3*3)
          CALL DZERO(ALPHA,3*3)
C
          WRITE(LUPRI,'(/,A,/)')
     &      ' Use function "Y = A + B*Tanh(C*N)" to approach the value'
C
          IF (ILAST .GE. 2)  THEN
C
            NLAST5 = MIN(NCOUNT,100)
C
            WRITE(LUPRI,'(/,A)')
     &      ' Estimation by Levenberg-Marquart method'
            WRITE(LUPRI,'(/,1X,A,3(F10.2),/)')
     &             ' Initial guess for A, B & C :',ACOEF(1),
     &                                    ACOEF(2),ACOEF(3)
C
            DO ITER = 1, 100
              ALAMBDA = -1.0D0
C
              CALL mrqmin(WORK(KWORK2),WORK(KWORK3),WORK(KWORK4),
     &                    NLAST5,ACOEF,NA,3,ACOVAR,ALPHA,3,
     &                    CHISQ,ALAMBDA,LUPRI)
C
              ITSTP = 0
              ISTP  = 1
C
 314          CONTINUE
C
              ACHISQ = CHISQ
C
              CALL mrqmin(WORK(KWORK2),WORK(KWORK3),WORK(KWORK4),
     &                    NLAST5,ACOEF,NA,3,ACOVAR,ALPHA,3,
     &                    CHISQ,ALAMBDA,LUPRI)
C
              DLTCHI = ABS(ACHISQ-CHISQ)
              ISTP = ISTP + 1
C
              IF (CHISQ .GT. ACHISQ) THEN
                ITSTP = 0
              ELSE IF(DLTCHI .LT. 0.1D0*THRS) THEN
                ITSTP = ITSTP + 1
              END IF
C Fusible
              IF (ISTP .GT. 1000) GOTO 315
C End Fusible
              IF (ITSTP .LT. 4) GOTO 314
C
c              WRITE(LUPRI,'(/,1X,A,2(F10.2),I5)')
c     &             ' The methods seems to be okay !!!'
C
 315          CONTINUE
C
              IF (ISTP .GT. 1000)  WRITE(LUPRI,'(/,1X,A,I5,/)')
     &             ' Something does not work very well:',ISTP
C
              ALAMBDA = 0.0D0
C
              CALL mrqmin(WORK(KWORK2),WORK(KWORK3),WORK(KWORK4),
     &                    NLAST5,ACOEF,NA,3,ACOVAR,ALPHA,3,
     &                    CHISQ,ALAMBDA,LUPRI)
C
            END DO
C
            WRITE(LUPRI,'(/,1X,A,3(F10.2),/)')
     &             ' Final guess for A, B & C :',ACOEF(1),
     &                                  ACOEF(2),ACOEF(3)
C
          END IF
C
          YESTIM1 = 0.0D0
          YESTIM2 = 0.0D0
          PNL = SSCCNI(NCSSCC,NCTERM)
C
          XESTIM1 = WORK(KWORK2-1+NCOUNT) + DFLOAT(NJUMPZ)
          IF (XESTIM1 .GT. NVARPT) XESTIM1 = WORK(KWORK2-1+NCOUNT)
          XESTIM2 = DFLOAT(NVARPT)
C
          CALL FUNCS(XESTIM1,ACOEF,YESTIM1,DYESTI,3,LUPRI)
          CALL FUNCS(XESTIM2,ACOEF,YESTIM2,DYESTI,3,LUPRI)
C
          PNP1 = YESTIM1
          PNP2 = YESTIM2
          RESSLT1 = 100.0D0*DABS((PNP1 - PNL)/PNL)
          RESSLT2 = 100.0D0*DABS((PNL - PNP2)/PNL)
C
          RESSLT = MAX(RESSLT1,RESSLT2)
C
          WRITE(LUPRI,'(/,1X,2(A,F14.2),A,F14.2)')' 7- The Ers give :',
     &                  RESSLT1,' &',RESSLT2,'  Required :',THRSLO
          WRITE(LUPRI,'(/,1X,3(A,F10.2))')
     &                  ' Resp. value :', PNL,
     &                  '  App. value :', PNP1, ' & :',PNP2
          WRITE(LUPRI,'(/,1X,A)')
     &          'The largest is used for comparizon with the required'
C
C----------------------------------------------------------------------------
C     Method 8: Using the Y = A + B*Tanh(C*N) function fitting we
C               get an i value and analyze E = |J(f) - J(i)| / J(f)
C               get the i_(f) value and analyze E = |J(f) - J(i_(f))| / J(f)
C               Here it is used NL2SOL tool.
C----------------------------------------------------------------------------
C
        ELSE IF (NCMTHD .EQ. 8) THEN
C
c          CALL DFAULT(IV,V)
C
          THRS = THRSLO*0.01D0
C
          DO I = 1, 3
            NA(I) = 1
          END DO
C
            ACOEF(1) = 0.1D0*SSCCNI(NCSSCC,NCTERM)
            ACOEF(2) = -2.0D0*SSCCNI(NCSSCC,NCTERM)
            ACOEF(3) = ABS(0.001D0*SSCCNI(NCSSCC,NCTERM))
C
C       For analitical Jacobian
C
c          CALL NL2SOL(n,p,x,testr,testj,iv,v,uiparm,urparm,testr)
C
C       For numerical Jacobian
C
C          call NL2SNO(n,p,x,testr,iv,v,uiparm,urparm,testr)
C
C
          WRITE(LUPRI,'(/,A,/)')
     &    '8- Use function "Y = A + B*Tanh(C*N)" to approach the value'
C
          IF (ILAST .GE. 2)  THEN
C
            NLAST5 = MIN(NCOUNT,100)
C
            WRITE(LUPRI,'(/,A)')
     &      ' Estimation by ASSESS-NL2SOL method'
            WRITE(LUPRI,'(/,1X,A,3(F10.2),/)')
     &             ' Initial guess for A, B & C :',ACOEF(1),
     &                                    ACOEF(2),ACOEF(3)
C
            WRITE(LUPRI,'(/,1X,A,3(F10.2),/)')
     &             ' Final guess for A, B & C :',ACOEF(1),
     &                                  ACOEF(2),ACOEF(3)
C
          END IF
C
          YESTIM1 = 0.0D0
          YESTIM2 = 0.0D0
          PNL = SSCCNI(NCSSCC,NCTERM)
C
          XESTIM1 = WORK(KWORK2-1+NCOUNT) + DFLOAT(NJUMPZ)
          IF (XESTIM1 .GT. NVARPT) XESTIM1 = WORK(KWORK2-1+NCOUNT)
          XESTIM2 = DFLOAT(NVARPT)
C
          CALL FUNCS(XESTIM1,ACOEF,YESTIM1,DYESTI,3,LUPRI)
          CALL FUNCS(XESTIM2,ACOEF,YESTIM2,DYESTI,3,LUPRI)
C
          PNP1 = YESTIM1
          PNP2 = YESTIM2
          RESSLT1 = 100.0D0*DABS((PNP1 - PNL)/PNL)
          RESSLT2 = 100.0D0*DABS((PNL - PNP2)/PNL)
C
          RESSLT = MAX(RESSLT1,RESSLT2)
C
          WRITE(LUPRI,'(/,1X,2(A,F14.2),A,F14.2)')' 8- The Ers give :',
     &                  RESSLT1,' &',RESSLT2,'  Required :',THRSLO
          WRITE(LUPRI,'(/,1X,3(A,F10.2))')
     &                  ' Resp. value :', PNL,
     &                  '  App. value :', PNP1, ' & :',PNP2
          WRITE(LUPRI,'(/,1X,A)')
     &          'The largest is used for comparizon with the required'
C
        ELSE IF (NCMTHD .GT. 7) THEN
C
          WRITE(LUPRI,*)
     &           'THE REQUIRED METHOD IS NOT AVAILABLE :',NCMTHD,
     &              '  IT SHOULD BE A NUMBER FROM 1 TO 7',NOMTHD
          CALL QUIT('WRONG METHOD !!! SERCON.1')
C
        END IF
C
      ELSE IF (ISINGUL .GT. 0) THEN
C
        XESTIM1 = WORK(KWORK2-1+NCOUNT) + DFLOAT(NJUMPZ)
        XESTIM2 = DFLOAT(NVARPT)
        PNL     = WORK(KWORK3-1+ILAST)
        PNP1    = 0.0D0
        PNP2    = 0.0D0
        RESSLT  = 1.0D2
        RESSLT1 = 1.0D2
        RESSLT2 = 1.0D2
C
        WRITE(LUPRI,'(/,1X,A,I6,2A,I3,/,3(A,F14.2),//)')
     &    'There are', ISINGUL,
     &    ' "singularities in reduced response equation" so',
     &    ' the approx. values are set to 0.0D0 for the chosen method',
     &     NCMTHD, ' produce  Ers. give :',
     &             RESSLT1,' &',RESSLT2,'  Required :',THRSLO
        WRITE(LUPRI,'(/,1X,3(A,F10.2),/)')
     &                ' Resp. value :', PNL,
     &                '  App. value :', PNP1, ' & :',PNP2
C
      END IF
C
      CALL QEXIT('SERCON')
C
      RETURN
C
      END
C
C*********************************************************************
C  /* Deck polint */
C
      SUBROUTINE polint(xa,ya,n,x,y,dy,LUPRI)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
C      PARAMETER (NMAX=1000)
      INTEGER i,m,ns
C      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      REAL*8 den,dif,dift,ho,hp,w,c(n),d(n)
C
c        WRITE(LUPRI,'(/,A)')' In POLINT :    Step,     X,       Y'
c        DO I = 1, n
c           WRITE(LUPRI,'(12X,2(F9.4))') xa(I),ya(I)
c        END DO
C
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.0d0) WRITE(LUPRI,*)'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C
C*********************************************************************
C  /* Deck ratint */
C
      SUBROUTINE ratint(xa,ya,n,x,y,dy,LUPRI)
      INTEGER n,NMAX, LUPRI
      REAL*8 dy,x,y,xa(n),ya(n),TINY
C      PARAMETER (NMAX=1000,TINY=1.e-25)
      PARAMETER (TINY=1.e-25)
      INTEGER i,m,ns
C      REAL*8 dd,h,hh,t,w,c(NMAX),d(NMAX)
      REAL*8 dd,h,hh,t,w,c(n),d(n)
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.) WRITE(LUPRI,*) 'failure in ratint'
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C
C*********************************************************************
C  /* Deck spline */
C
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
C      PARAMETER (NMAX=1000)
      INTEGER i,k
C      REAL*8 p,qn,sig,un,u(NMAX)
      REAL*8 p,qn,sig,un,u(n)
      if (yp1.gt..99e30) then
        y2(1)=0.0D0
        u(1)=0.0D0
      else
        y2(1)=-0.5D0
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt.0.99D30) then
        qn=0.0D0
        un=0.0D0
      else
        qn=0.5D0
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C
C*********************************************************************
C  /* Deck splint */
C
      SUBROUTINE splint(xa,ya,y2a,n,x,y,LUPRI)
      INTEGER n,LUPRI
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) WRITE(LUPRI,*) 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0D0
      return
      END
C
C*********************************************************************
C  /* Deck mrqmin */
C
      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *alamda,LUPRI)
C
C INPUT PARAMETERS
C
C       X,Y,SIG - X, Y AND SIGMA VECTORS (REAL*8)
C       NDATA - NUMBER OF POINTS IN X,Y AND SIG (INTEGER)
C       A - COEF. OF FIT (REAL*8)
C       MA - NUMBER OF COEFFICIENTS IN A (INTEGER)
C       LISTA - VECTOR LIST OF COEF. NUMBER TO FIT (INTEGER)
C       MFIT - NUMBER OF COEF. TO FIT (LENGTH OF LISTA) (INTEGER)
C       NCA - DIMENSIONS OF COVAR AND ALPHA
C       FUNCS - FUNCTION TO FIT (SPECIFED BY)
C               SUBROUTINE FUNCS(X,A,YFIT,DYDA,MA)
C               THIS EVALUATES THE FUNCTION YFIT AT SCALAR
C                       VALUE X USING THE MA COEF. IN A.  DYDA IS A
C                       VECTOR OF PARTIAL DIREVATIVES WITH RESPECT TO
C                       EACH COEF.
C       ALAMDA - REAL*8 SCALAR
C               ON FIRST CALL SET ALAMDA<0 FOR INITIALIZATION
C               SUBSEQUENT CALLS USE ALAMDA RETURNED BY PREVIOUS
C                       ITERATION
C               ON LAST CALL SET ALAMBDA TO 0.0 TO COMPUTE COVARIANCE
C                       AND CURVATURE MATRIX
C
C OUTPUT PARAMETERS
C
C       A - UPDATED COEFFICIENTS FOR FIT
C       COVAR - COVARIANCE MATRIX (REAL*8)
C       ALPHA - CURVATURE MATRIX (REAL*8)
C       CHISQ - CHI SQUARED OF FIT (REAL*8)
C       ALAMDA - UPDATED VALUE
C-----------------------------------------------------------------------
C      INTEGER ma,nca,ndata,ia(ma),MMAX,LUPRI
      INTEGER ma,nca,ndata,ia(ma),LUPRI
      REAL*8 alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),x(ndata),y(ndata)
CPFP
C      REAL*8 alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),
C     *sig(ndata),x(ndata),y(ndata)
CPFP remove NMAX
      PARAMETER (MMAX=500)
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,mfit
      REAL*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
c      REAL*8 ochisq,atry(ma),beta(ma),da(ma)
C
      SAVE ochisq,atry,beta,da,mfit
Cend-PFP
      if(alamda.lt.0.0D0)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001D0
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,LUPRI)
CPFP
C        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,
C     *funcs,LUPRI)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      do 14 j=1,mfit
        do 13 k=1,mfit
          covar(j,k)=alpha(j,k)
13      continue
        covar(j,j)=alpha(j,j)*(1.0D0+alamda)
        da(j)=beta(j)
14    continue
      call gaussj(covar,mfit,nca,da,1,1,LUPRI)
      if(alamda.eq.0.0D0)then
        call covsrt(covar,nca,ma,ia,mfit)
        call covsrt(alpha,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
C
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,LUPRI)
C
      if(chisq.lt.ochisq)then
        alamda=0.1D0*alamda
        ochisq=chisq
        do 17 j=1,mfit
          do 16 k=1,mfit
            alpha(j,k)=covar(j,k)
16        continue
          beta(j)=da(j)
17      continue
        do 18 l=1,ma
          a(l)=atry(l)
18      continue
      else
        alamda=1.0D1*alamda
        chisq=ochisq
      endif
C
      return
      END
C
C*********************************************************************
C  /* Deck mrqcof */
C
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
     *LUPRI)
CPFP
C      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
C     *funcs,LUPRI)
      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     *y(ndata)
      EXTERNAL funcs
CPFP  removed NMAX
      PARAMETER (MMAX=500)
      INTEGER mfit,i,j,k,l,m
      REAL*8 dy,sig2i,wt,ymod,dyda(MMAX)
c      REAL*8 dy,sig2i,wt,ymod,dyda(ma)
Cend-PFP
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.0D0
12      continue
        beta(j)=0.0D0
13    continue
      chisq=0.0D0
c        WRITE(LUPRI,'(/,1X,A,4(F10.2))')
c     &         ' Chi squared IN MRQCOF 1:',chisq,a(1),a(2),a(3)
      do 16 i=1,ndata
        call funcs(x(i),a,ymod,dyda(i),ma,LUPRI)
        sig2i=1.0D0/(sig(i)*sig(i))
        dy=y(i)-ymod
        j=0
        do 15 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            wt=dyda(l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
c        WRITE(LUPRI,'(/,1X,A,I5,3(F10.2))')
c     *         ' Chi squared IN MRQCOF 2:',ndata,chisq,ymod,dyda
c        WRITE(LUPRI,'(1X,A,I5,3(F10.2))')
c     *         ' Chi squared IN MRQCOF 2-A:',ndata,chisq,sig(i),sig2i
16    continue
c        WRITE(LUPRI,'(/,1X,A,(F10.2))')
c     &         ' Chi squared IN MRQCOF 3:',chisq
c      stop
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END
C
C*********************************************************************
C  /* Deck gaussj */
C
      SUBROUTINE gaussj(a,n,np,b,m,mp,LUPRI)
      INTEGER m,mp,n,np,NMAX,LUPRI
      REAL*8 a(np,np),b(np,mp)
CPFP removed NMAX
      PARAMETER (NMAX=500)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
c      INTEGER i,icol,irow,j,k,l,ll,indxc(n),indxr(n),ipiv(n)
Cend-PFP
      REAL*8 big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.0D0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.0D0)
     *      WRITE(LUPRI,*) ' iter. n.',i,'  singular matrix in gaussj'
        pivinv=1.0D0/a(icol,icol)
        a(icol,icol)=1.0D0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.0D0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
C
C*********************************************************************
C  /* Deck covsrt */
C
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      INTEGER ma,mfit,npc,ia(ma)
      REAL*8 covar(npc,npc)
      INTEGER i,j,k
      REAL*8 swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.0D0
          covar(j,i)=0.0D0
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END
C
C*********************************************************************
C  /* Deck FUNCS */
C
      SUBROUTINE funcs(x,a,ymod,dyda,ma,LUPRI)
C
      INTEGER ma
      REAL*8 a(ma),x,dyda(ma),ymod
C
c        WRITE(LUPRI,'(/,1X,A,(F10.2))')
c     &         ' Chi squared IN FUNCS 1:',x
C
      ymod = a(1) + a(2)*TANH(a(3)*x)
C
      aux1 = (COSH(a(3)*x))**2
      dyda(1) = 1.0D0
      dyda(2) = TANH(a(3)*x)
      dyda(3) = (a(2)*x)/aux1
C
c        WRITE(LUPRI,'(/,1X,A,2(F10.2))')
c     &   ' Chi squared IN FUNCS 2:',ymod,dyda(1),dyda(2),dyda(3)
C
      RETURN
      END
C
C*********************************************************************
