!
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
      SUBROUTINE SOODSO_1(IPRSOS,NLOCC,CLOCO,SPNDSO,AINTREP,AINTADR,
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
      DIMENSION   AINTREP(9*MXCENT), AINTADR(9*MXCENT)
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
     &            LABINT,AINTREP,AINTADR,MPQUAD,.FALSE.,NPATOM,.TRUE.,
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
         KSYM   = AINTREP(ICOMP) + 1
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
C TO DO:
C a. 3 is the number of parameters to be found in methods 7 and 8.
C    Eventualy it has to be moved to locing.h
C b. UIPARM, URPARM, UFPARM has to desapear.
C c. FCN, FUNCS, FTANH2,FTANH3, CALCR,CALCJ should converge to only 
C    one unique function/subroutine
C
      PARAMETER (NP = 3)
      INTEGER NSING,INFO,IPVT(NP),UIPARM(1),IV(60+NP)
      REAL*8 CHISQ, ALAMBDA
      DIMENSION WORK(LWORK)
      DIMENSION ACOEF(NP),ACOVAR(NP,NP),ALPHA(NP,NP),DYESTI(NP)
      DIMENSION NA(NP),URPARM(1)
      EXTERNAL FUNCS,FTANH2,FTANH3,CALCJ,CALCR,UFPARM
      CALL QENTER('SERCON')
C
C   If there is a "singularity in reduced response equation"
C   then avoid any further calculation
C
      IF (ISINGUL .EQ. 0) THEN
C
C  INDEX1 should coincide with NCOUNT maximum !
C
        INDEX1 = INT((INDEXT - NEXCII)/NJUMPZ + 1)
C
        KWORK1 = 1
        KWORK2 = KWORK1 + INDEX1
        KWORK3 = KWORK2 + INDEX1
        KWORK4 = KWORK3 + INDEX1
        KWORK5 = KWORK4 + NP*INDEX1
        KWORK6 = KWORK5 + 93 + 2*NP*INDEX1 + NP*(3*NP+33)
        LWORK1 = LWORK  - KWORK6
C
        IF (LWORK7 .LT. 0) THEN
          WRITE(LUPRI,*)'SERCON.1: Need :',KWORK6,
     &                '  Available :',LWORK1
          CALL QUIT('Insufficient memory in SOSDRV.1')
        END IF
C
        CALL DZERO(WORK(KWORK1),INDEX1)
        CALL DZERO(WORK(KWORK2),INDEX1)
        CALL DZERO(WORK(KWORK3),INDEX1)
        CALL DZERO(WORK(KWORK4),INDEX1)
        CALL DZERO(WORK(KWORK5),NP*INDEX1)
        CALL DZERO(WORK(KWORK6),(93 + 2*NP*INDEX1 + NP*(3*NP+33)))
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
           END IF
        END DO
C
C       CALL DCOPY(INDEX1,WORK(KWORK1),1,SSCCPE(INDEX1,5),1)
        WRITE(LUPRI,'(/,1X,A,2(I4))')
     &    ' Series to be analyzed :',INDEX1+1,INDEXT
C
        WRITE(LUPRI,'(/,A)')' In SERCON :    Step,      X,        Y'
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
     &                                   '  App. value :', PNP
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
C               infer a J(i+1) value and analyze E = |J(r) - J(i+1)| / J(r)
C               infer the J(f) value and analyze E = |J(r) - J(f)| / J(r)
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
C               infer a J(i+1) value and analyze E = |J(r) - J(i+1)| / J(r)
C               infer the J(f) value and analyze E = |J(r) - J(f)| / J(r)
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
C               infer a J(i+1) value and analyze E = |J(r) - J(i+1)| / J(r)
C               infer the J(f) value and analyze E = |J(r) - J(f)| / J(r)
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
     &                  NLAST5,XESTIM1,YESTIM1,LUPRI)
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
     &                  NLAST5,XESTIM2,YESTIM2,LUPRI)
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
C               infer a J(i+1) value and analyze E = |J(r) - J(i+1)| / J(r)
C               infer the J(f) value and analyze E = |J(r) - J(f)| / J(r)
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
            ACOEF(1) =     -0.5D0*SSCCNI(NCSSCC,NCTERM)
            ACOEF(2) =      2.0D0*SSCCNI(NCSSCC,NCTERM)
            ACOEF(3) = ABS(0.01D0*SSCCNI(NCSSCC,NCTERM))
C
c            ACOEF(1) = 1.0D0
c            ACOEF(2) = 1.0D0
c            ACOEF(3) = 1.0D0
C
          CALL DZERO(ACOVAR,3*3)
          CALL DZERO(ALPHA,3*3)
C
          WRITE(LUPRI,'(/,A,/)')
     &      ' Use function "Y = A + B*Tanh(C*N)" to approach the value'
C
          IF (ILAST .GE. 10)  THEN
C
            NLAST5 = MIN(NCOUNT,100)
C
            WRITE(LUPRI,'(/,A)')
     &      ' Estimation by Levenberg-Marquart method'
            WRITE(LUPRI,'(/,1X,A,3(F10.2),/)')
     &             ' Initial guess for A, B & C :',
     &               ACOEF(1),ACOEF(2),ACOEF(3)
C
C-----------------------------------------------------------------
C      Calculates the mean value, variance and standard deviation
C      of the y_i values
C-----------------------------------------------------------------
C
          AMEAN = 0.0D0
C
          DO I = 1, NCOUNT
            AMEAN = AMEAN + WORK(KWORK3-1+I)
          END DO
C
C       Mean value : sum_i [WORK(KWORK3-1+I)]/NCOUNT
C
          AMEAN = AMEAN/NCOUNT
C
C       variance:                [(WORK(KWORK3-1+I)-AMEAN)**2)]/NCOUNT
C       Standard deviation: sqrt(((WORK(KWORK3-1+I)-AMEAN)**2)/NCOUNT)
C
          DO I = 1, NCOUNT
            WORK(KWORK4-1+I) =
     &                     sqrt(((WORK(KWORK3-1+I)-AMEAN)**2)/NCOUNT)
          END DO
C
            WRITE(LUPRI,'(/,A)')
     &      '   X  ,   Y    and Standard deviation'
C
          DO I = 1, NCOUNT
             WRITE(LUPRI,'(12X,3(F9.4))') WORK(KWORK2-1+I),
     &                   WORK(KWORK3-1+I),WORK(KWORK4-1+I)
          END DO
C
            WRITE(LUPRI,'(/,1X,A,(F10.2),/)')
     &             ' Mean Value :', AMEAN
C
            DO ITER = 1, 4
              ALAMBDA = -1.0D0
C
              CALL mrqmin(WORK(KWORK2),WORK(KWORK3),WORK(KWORK4),
     &                    NLAST5,ACOEF,NA,3,ACOVAR,ALPHA,3,
     &                    CHISQ,FUNCS,ALAMBDA,LUPRI,NSING)
C
              IF (NSING .NE. 0) THEN
                WRITE(NPRINT,'(/,1X,A,I4/)')
     &              ' 1. Singular matrix in gaussj',NSING
                GOTO 316
              END IF
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
     &                    CHISQ,FUNCS,ALAMBDA,LUPRI,NSING)
C
              IF (NSING .NE. 0) THEN
                WRITE(NPRINT,'(/,1X,A,I4/)')
     &              ' 2. Singular matrix in gaussj',NSING
                GOTO 316
              END IF
C
              DLTCHI = ABS(ACHISQ-CHISQ)
              ISTP = ISTP + 1
C
              IF (CHISQ .GT. ACHISQ) THEN
                ITSTP = 0
              ELSE IF(DLTCHI .LT. 0.1D0) THEN
                ITSTP = ITSTP + 1
              END IF
C Fusible
              IF (ISTP .GT. 100) THEN
                WRITE(LUPRI,'(/,1X,A,I5,/)')
     &             ' Something does not work very well:',ISTP
                GOTO 315
              ENDIF
C End Fusible
              IF (ITSTP .LT. 10) GOTO 314
C
              WRITE(LUPRI,'(/,1X,A,2(F10.2),I5)')
     &             ' The methods seems to be okay !!!'
C
 315          CONTINUE
C
              ALAMBDA = 0.0D0
C
              CALL mrqmin(WORK(KWORK2),WORK(KWORK3),WORK(KWORK4),
     &                    NLAST5,ACOEF,NA,3,ACOVAR,ALPHA,3,
     &                    CHISQ,FUNCS,ALAMBDA,LUPRI,NSING)
C
              IF (NSING .NE. 0) THEN
                WRITE(NPRINT,'(/,1X,A,I4/)')
     &              ' 3. Singular matrix in gaussj',NSING
                GOTO 316
              END IF
C
            END DO
C
 316        CONTINUE
C
            WRITE(LUPRI,'(/,1X,A,I5,A,I5,A,/)')
     &     ' Finally it seems that the program reached somewhere in :',
     &               ISTP,' and ',ITSTP,' steps'
C
            WRITE(LUPRI,'(/,1X,A,3(F10.2),/)')
     &             ' Final guess for A, B & C :',
     &                  ACOEF(1),ACOEF(2),ACOEF(3)
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
C               infer a J(i+1) value and analyze E = |J(r) - J(i+1)| / J(r)
C               infer the J(f) value and analyze E = |J(r) - J(f)| / J(r)
C               Here it is used MinPack tool.
C----------------------------------------------------------------------------
C
        ELSE IF (NCMTHD .EQ. 8) THEN
C
CPFP 2017 I still do not know how to change the precision within the method ...
C
c          CALL DFAULT(IV,V)
C
          THRS  = THRSLO*0.01D0
          INFO  = 0
          IWRK6 = 5*3+NCOUNT
C
          DO I = 1, 3
            NA(I) = 1
          END DO
C
            ACOEF(1) =     -0.5D0*SSCCNI(NCSSCC,NCTERM)
            ACOEF(2) =      2.0D0*SSCCNI(NCSSCC,NCTERM)
            ACOEF(3) = ABS(0.01D0*SSCCNI(NCSSCC,NCTERM))
c            ACOEF(1) = 0.1D0*SSCCNI(NCSSCC,NCTERM)
c            ACOEF(2) = -2.0D0*SSCCNI(NCSSCC,NCTERM)
c            ACOEF(3) = ABS(0.001D0*SSCCNI(NCSSCC,NCTERM))
C
          WRITE(LUPRI,'(/,A,/)')
     &    '8- Use function "Y = A + B*Tanh(C*N)" to approach the value'
C
          IF (ILAST .GE. 10)  THEN
C
            NLAST5 = MIN(NCOUNT,100)
C
            WRITE(LUPRI,'(/,A)')
     &      ' Estimation by MinPack - Levenberg-Marquart method :'
            WRITE(LUPRI,'(/,1X,A,3(F10.2),/)')
     &      ' Initial guess for A, B & C :',ACOEF(1),
     &                                      ACOEF(2),ACOEF(3)
C
C       For analitical Jacobian
C
c          call lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)
           call lmder1(FTANH2,NCOUNT,NP,WORK(KWORK2),WORK(KWORK3),
     &                 ACOEF,WORK(KWORK4),WORK(KWORK5),NCOUNT,THRS,
     &                 INFO,IPVT,WORK(KWORK6),IWRK6,LUPRI)
C
C       For numerical Jacobian
C           To be implemented !!!
C          call lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)
C
          WRITE(LUPRI,'(/,1X,A,3(F10.2),/)')
     &             ' Final guess for A, B & C :',ACOEF(1),
     &                                  ACOEF(2),ACOEF(3)
c
c            WRITE(LUPRI,'(/,1X,23A,I5,/)')
cfor lmder1
c     &    ' info = 0  improper input parameters.',
c     &    '',
c     &    ' info = 1  algorithm estimates that the relative error',
c     &    '           in the sum of squares is at most tol.',
c     &    '',
c     &    ' info = 2  algorithm estimates that the relative error',
c     &    '           between x and the solution is at most tol.',
c     &    '',
c     &    ' info = 3  conditions for info = 1 and info = 2 both hold.',
c     &    '',
c     &    ' info = 4  fvec is orthogonal to the columns of the',
c     &    '           jacobian to machine precision.',
c     &    '',
c     &    ' info = 5  number of calls to fcn with iflag = 1 has',
c     &    '           reached 100*(n+1).',
c     &    '',
c     &    ' info = 6  tol is too small. no further reduction in',
c     &    '           the sum of squares is possible.',
c     &    '',
c     &    ' info = 7  tol is too small. no further improvement in',
c     &    '           the approximate solution x is possible.',
c     &    '',
c     &    ' Info variable from MinPack : ',INFO
C
            WRITE(LUPRI,'(/,1X,A,I5,/)')
     &    ' Info variable from MinPack : ',INFO
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
C----------------------------------------------------------------------------
C     Method 9: Using the Y = A + B*Tanh(C*N) function fitting we
C               infer a J(i+1) value and analyze E = |J(r) - J(i+1)| / J(r)
C               infer the J(f) value and analyze E = |J(r) - J(f)| / J(r)
C               Here it is used NL2SOL tool.
C----------------------------------------------------------------------------
C
        ELSE IF (NCMTHD .EQ. 9) THEN
C
CPFP 2017 I still do not know how to change the precision within the method ...
C
          THRS  = THRSLO*0.01D0
          INFO  = 0
          IWRK6 = 5*3+NCOUNT
C
          DO I = 1, 3
            NA(I) = 1
          END DO
C
            ACOEF(1) =     -0.5D0*SSCCNI(NCSSCC,NCTERM)
            ACOEF(2) =      2.0D0*SSCCNI(NCSSCC,NCTERM)
            ACOEF(3) = ABS(0.01D0*SSCCNI(NCSSCC,NCTERM))
c            ACOEF(1) = 0.1D0*SSCCNI(NCSSCC,NCTERM)
c            ACOEF(2) = -2.0D0*SSCCNI(NCSSCC,NCTERM)
c            ACOEF(3) = ABS(0.001D0*SSCCNI(NCSSCC,NCTERM))
C
          WRITE(LUPRI,'(/,A,/)')
     &    '9- Use function "Y = A + B*Tanh(C*N)" to approach the value'
C
          IF (ILAST .GE. 10)  THEN
C
            CALL DFAULT(IV,WORK(KWORK6))
            IV(1) = 0
C
            NLAST5 = MIN(NCOUNT,100)
C
            WRITE(LUPRI,'(/,A)')
     &      ' Estimation by NL2SOL - Levenberg-Marquart method :'
            WRITE(LUPRI,'(/,1X,A,3(F10.2),/)')
     &      ' Initial guess for A, B & C :',ACOEF(1),
     &                                      ACOEF(2),ACOEF(3)
C
C       For analitical Jacobian
C
C      SUBROUTINE NL2SOL(N, P, X, XSTP, YSTP, CALCR, CALCJ, IV, V,
C     +                  UIPARM, URPARM, UFPARM, NLUPRI)
C
           call NL2SOL(NCOUNT,NP,ACOEF,WORK(KWORK2),WORK(KWORK3),
     &                 CALCR,CALCJ,IV,WORK(KWORK6),
     &                 UIPARM,URPARM,UFPARM,LUPRI)
C
C       For numerical Jacobian
C
C        To be implemented.
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
          WRITE(LUPRI,'(/,1X,2(A,F14.2),A,F14.2)')' 9- The Ers give :',
     &                  RESSLT1,' &',RESSLT2,'  Required :',THRSLO
          WRITE(LUPRI,'(/,1X,3(A,F10.2))')
     &                  ' Resp. value :', PNL,
     &                  '  App. value :', PNP1, ' & :',PNP2
          WRITE(LUPRI,'(/,1X,A)')
     &          'The largest is used for comparizon with the required'
C

        ELSE IF (NCMTHD .GT. 9) THEN
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
     *funcs,alamda,LUPRI,NSING)
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
      INTEGER ma,nca,ndata,ia(ma),LUPRI,NSING
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
      EXTERNAL FUNCS
C
      SAVE ochisq,atry,beta,da,mfit
Cend-PFP
      if(alamda.lt.0.0D0)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001D0
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,
     &              funcs,LUPRI)
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
      call gaussj(covar,mfit,nca,da,1,1,LUPRI,NSING)
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
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,
     &            funcs,LUPRI)
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
c      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
c     *LUPRI)
CPFP
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
     *funcs,LUPRI)
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
      SUBROUTINE gaussj(a,n,np,b,m,mp,LUPRI,NSINGMG)
      INTEGER m,mp,n,np,NMAX,LUPRI
      REAL*8 a(np,np),b(np,mp)
CPFP removed NMAX
      PARAMETER (NMAX=500)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      INTEGER NSINGMG
c      INTEGER i,icol,irow,j,k,l,ll,indxc(n),indxr(n),ipiv(n)
Cend-PFP
      REAL*8 big,dum,pivinv
      NSINGMG = 0
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
        if (a(icol,icol).eq.0.0D0) then
          WRITE(LUPRI,*) ' GAUSSJ: singular matrix'
          NSINGMG =  NSINGMG + 1
        endif
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
**********************************************************************
*  /* Deck FTANH2 */
*
      subroutine ftanh2(m,n,xstp,ystp,x,fvec,fjac,ldfjac,iflag,NLUPRI)
c
c   If iflag = 1
c     Calculate the residual of the function  Y(n) = A + B*tanh(C*X(n))
c     this implies the function   0 = A + B*tanh(C*X(n)) - Y(n) with this
c     residual saved in fvec(m) which has to be minimized.
c
c   If iflag = 2
c     Calculate the Jacobian matrix (3x1) in the fjac(ldfjac,n) matriz
c     dY/dA = 1 ; dY/dB = tanh(C*X) and dY/dC = (B*X) / (COSH(C*X))**2
c
c   Variables:
c    x(n)   : Input/Output Initial guess and new approximation of the
c             coeficients to be found out.
c    n      : number of unknowns (3 in the present example).
c    m      : number of points to fit the curve.
c    ldfjac : coincide with m.
c    fvec(m): Output residuals.
c    fjac(m): Output jacobian of the function to be fit
c            (whch coinside with the jacobian of the residuals).
c    xstp(m): Input dots for which there excist a value of the function.
c    ystp(m): Input values of the function.
c
      integer m,n,ldfjac,iflag,NLUPRI
      double precision x(n),fvec(m),fjac(ldfjac,n),xstp(m),ystp(m)
c
      if ( iflag .eq. 1 ) then
*
        do i = 1, m
          fvec(i) = x(1) + x(2)*TANH(x(3)*xstp(i)) - ystp(i)
        end do
*
      else if ( iflag .eq. 2 ) then
*
        do i = 1, m
          aux1 = (COSH(x(3)*xstp(i)))**2
          fjac(i,1) = 1.0D0
          fjac(i,2) = TANH(x(3)*xstp(i))
          fjac(i,3) = (x(2)*xstp(i))/aux1
        end do
*
      end if
c
      return
      end
*
**********************************************************************
*  /* Deck FTANH3 */
*
      subroutine ftanh3(m,n,xstp,ystp,x,fvec,iflag,NLUPRI)
c
c   If iflag = 1
c     Calculate the residual of the function  Y(n) = A + B*tanh(C*X(n))
c     this implies the function   0 = A + B*tanh(C*X(n)) - Y(n) with this
c     residual saved in fvec(m) which has to be minimized.
c
c   If iflag = 2
c     Calculate the Jacobian matrix (3x1) in the fjac(ldfjac,n) matriz
c     dY/dA = 1 ; dY/dB = tanh(C*X) and dY/dC = (B*X) / (COSH(C*X))**2
c
c   Variables:
c    x(n)   : Input/Output Initial guess and new approximation of the
c             coeficients to be found out.
c    n      : number of unknowns (3 in the present example).
c    m      : number of points to fit the curve.
c    ldfjac : coincide with m.
c    fvec(m): Output residuals.
c    fjac(m): Output jacobian of the function to be fit
c            (whch coinside with the jacobian of the residuals).
c    xstp(m): Input dots for which there excist a value of the function.
c    ystp(m): Input values of the function.
c
      integer m,n,iflag,NLUPRI
      double precision x(n),fvec(m),xstp(m),ystp(m)
c
c      if ( iflag .eq. 1 ) then
*
        do i = 1, m
          fvec(i) = x(1) + x(2)*TANH(x(3)*xstp(i)) - ystp(i)
        end do
*
c      else if ( iflag .eq. 2 ) then
c*
c        do i = 1, m
c          aux1 = (COSH(x(3)*xstp(i)))**2
c          fjac(i,1) = 1.0D0
c          fjac(i,2) = TANH(x(3)*xstp(i))
c          fjac(i,3) = (x(2)*xstp(i))/aux1
c        end do
c*
c      end if
c
      return
      end
*
**********************************************************************
*  /* Deck CALCR */
*
c           CALL CALCR(N,P,X,NF,R,UIPARM,URPARM,UFPARM)
      SUBROUTINE CALCR(m,n,xstp,ystp,x,iflag,fvec,NLUPRI)
*
      integer m,n,iflag,NLUPRI,NF
      double precision x(n),fvec(m),xstp(m),ystp(m)
*
        do i = 1, m
          fvec(i) = x(1) + x(2)*TANH(x(3)*xstp(i)) - ystp(i)
        end do
*
      RETURN
      END
*
**********************************************************************
*  /* Deck CALCJ */
*
c           CALL CALCJ(N,P,X,NF,J,UIPARM,URPARM,UFPARM)
      SUBROUTINE CALCJ(m,n,xstp,ystp,x,iflag,fjac,ldfjac,NLUPRI)
*
      integer m,n,ldfjac,iflag,NLUPRI,NF
      double precision x(n),fjac(m,ldfjac),xstp(m),ystp(m)
*
        do i = 1, m
          aux1 = (COSH(x(3)*xstp(i)))**2
          fjac(i,1) = 1.0D0
          fjac(i,2) = TANH(x(3)*xstp(i))
          fjac(i,3) = (x(2)*xstp(i))/aux1
        end do
*
c      do i = 1, m
c        WRITE(NLUPRI,'(1X,A,3(F12.6))')
c     &   ' IN CALCJ 1:', fjac(i,1),fjac(i,2),fjac(i,3)
c      end do
*
      RETURN
      END
*
**********************************************************************
*  PFP 2017: TAKEN FROM MinPack and slightly modified
**********************************************************************
*  /* Deck LMDER1 */
*
      subroutine lmder1(ftanh2,m,n,xstp,ystp,x,fvec,fjac,ldfjac,tol,
     *                  info,ipvt,wa,lwa,NLUPRI)
      integer m,n,ldfjac,info,lwa,NLUPRI
      integer ipvt(n)
      double precision tol,xstp(m),ystp(m)
      double precision x(n),fvec(m),fjac(ldfjac,n),wa(lwa)
      external ftanh2
c     **********
c
c     subroutine lmder1
c
c     the purpose of lmder1 is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of the
c     levenberg-marquardt algorithm. this is done by using the more
c     general least-squares solver lmder. the user must provide a
c     subroutine which calculates the functions and the jacobian.
c
c     the subroutine statement is
c
c       subroutine lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,
c                         ipvt,wa,lwa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions and the jacobian. fcn must
c         be declared in an external statement in the user
c         calling program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
c         integer m,n,ldfjac,iflag
c         double precision x(n),fvec(m),fjac(ldfjac,n)
c         ----------
c         if iflag = 1 calculate the functions at x and
c         return this vector in fvec. do not alter fjac.
c         if iflag = 2 calculate the jacobian at x and
c         return this matrix in fjac. do not alter fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmder1.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
c       fjac is an output m by n array. the upper n by n submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       tol is a nonnegative input variable. termination occurs
c         when the algorithm estimates either that the relative
c         error in the sum of squares is at most tol or that
c         the relative error between x and the solution is at
c         most tol.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  algorithm estimates that the relative error
c                   in the sum of squares is at most tol.
c
c         info = 2  algorithm estimates that the relative error
c                   between x and the solution is at most tol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  fvec is orthogonal to the columns of the
c                   jacobian to machine precision.
c
c         info = 5  number of calls to fcn with iflag = 1 has
c                   reached 100*(n+1).
c
c         info = 6  tol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  tol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c       ipvt is an integer output array of length n. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       wa is a work array of length lwa.
c
c       lwa is a positive integer input variable not less than 5*n+m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... lmder
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer maxfev,mode,nfev,njev,nprint
      double precision factor,ftol,gtol,xtol,zero
      data factor,zero /1.0d2,0.0d0/
      info = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m .or. tol .lt. zero
     *    .or. lwa .lt. 5*n + m) go to 10
c
c     call lmder.
c
      maxfev = 100*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      mode = 1
      nprint = 0
C
C     PRINT THE PARAMETERS
C
c      WRITE(NLUPRI,'(/,1X,A,F14.2)')' IN LMDER1 1: ', tol
C
      call lmder(ftanh2,m,n,xstp,ystp,x,fvec,fjac,ldfjac,ftol,xtol,
     *           gtol,maxfev,wa(1),mode,factor,nprint,info,nfev,njev,
     *           ipvt,wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1),
     *           NLUPRI)
      if (info .eq. 8) info = 4
c
c      WRITE(NLUPRI,'(/,1X,A,F14.2)')' IN LMDER1 2: ', tol
c
   10 continue
      return
c
c     last card of subroutine lmder1.
c
      end
*
**********************************************************************
*  /* Deck LMDER */
*
      subroutine lmder(ftanh2,m,n,xstp,ystp,x,fvec,fjac,ldfjac,ftol,
     *                 xtol,gtol,maxfev,diag,mode,factor,nprint,info,
     *                 nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4,NLUPRI)
      integer m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev,NLUPRI
      integer ipvt(n)
      double precision ftol,xtol,gtol,factor,xstp(m),ystp(m)
      double precision x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n),
     *                 wa1(n),wa2(n),wa3(n),wa4(m)
      external ftanh2
c     **********
c
c     subroutine lmder
c
c     the purpose of lmder is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of the
c     levenberg-marquardt algorithm. the user must provide a
c     subroutine which calculates the functions and the jacobian.
c
c     the subroutine statement is
c
c       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
c                        maxfev,diag,mode,factor,nprint,info,nfev,
c                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions and the jacobian. fcn must
c         be declared in an external statement in the user
c         calling program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
c         integer m,n,ldfjac,iflag
c         double precision x(n),fvec(m),fjac(ldfjac,n)
c         ----------
c         if iflag = 1 calculate the functions at x and
c         return this vector in fvec. do not alter fjac.
c         if iflag = 2 calculate the jacobian at x and
c         return this matrix in fjac. do not alter fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmder.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
c       fjac is an output m by n array. the upper n by n submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       ftol is a nonnegative input variable. termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most ftol.
c         therefore, ftol measures the relative error desired
c         in the sum of squares.
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol. therefore, xtol measures the
c         relative error desired in the approximate solution.
c
c       gtol is a nonnegative input variable. termination
c         occurs when the cosine of the angle between fvec and
c         any column of the jacobian is at most gtol in absolute
c         value. therefore, gtol measures the orthogonality
c         desired between the function vector and the columns
c         of the jacobian.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn with iflag = 1
c         has reached maxfev.
c
c       diag is an array of length n. if mode = 1 (see
c         below), diag is internally set. if mode = 2, diag
c         must contain positive entries that serve as
c         multiplicative scale factors for the variables.
c
c       mode is an integer input variable. if mode = 1, the
c         variables will be scaled internally. if mode = 2,
c         the scaling is specified by the input diag. other
c         values of mode are equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.).100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x, fvec, and fjac
c         available for printing. fvec and fjac should not be
c         altered. if nprint is not positive, no special calls
c         of fcn with iflag = 0 are made.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  both actual and predicted relative reductions
c                   in the sum of squares are at most ftol.
c
c         info = 2  relative error between two consecutive iterates
c                   is at most xtol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  the cosine of the angle between fvec and any
c                   column of the jacobian is at most gtol in
c                   absolute value.
c
c         info = 5  number of calls to fcn with iflag = 1 has
c                   reached maxfev.
c
c         info = 6  ftol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  xtol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c         info = 8  gtol is too small. fvec is orthogonal to the
c                   columns of the jacobian to machine precision.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn with iflag = 1.
c
c       njev is an integer output variable set to the number of
c         calls to fcn with iflag = 2.
c
c       ipvt is an integer output array of length n. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       qtf is an output array of length n which contains
c         the first n elements of the vector (q transpose)*fvec.
c
c       wa1, wa2, and wa3 are work arrays of length n.
c
c       wa4 is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
c
c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,iflag,iter,j,l
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,
     *                 one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,
     *                 sum,temp,temp1,temp2,xnorm,zero
      double precision dpmpar,enorm
      data one,p1,p5,p25,p75,p0001,zero
     *     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
c
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In LMDER 1: A, B & C :',x(1),x(2),x(3)
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m
     *    .or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero
     *    .or. maxfev .le. 0 .or. factor .le. zero) go to 300
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
c
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In LMDER 2: A, B & C :',x(1),x(2),x(3)
c
c     evaluate the function at the starting point
c     and calculate its norm.
c
      iflag = 1
      call ftanh2(m,n,xstp,ystp,x,fvec,fjac,ldfjac,iflag,NLUPRI)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(m,fvec)
c
c     initialize levenberg-marquardt parameter and iteration counter.
c
      par = zero
      iter = 1
c
c     beginning of the outer loop.
c
   30 continue
c
c        calculate the jacobian matrix.
c
         iflag = 2

c      WRITE(NLUPRI,'(/,1X,A)')' IN LMDER 3: '

         call ftanh2(m,n,xstp,ystp,x,fvec,fjac,ldfjac,iflag,NLUPRI)
         njev = njev + 1
         if (iflag .lt. 0) go to 300
c
c        if requested, call fcn to enable printing of iterates.
c
         if (nprint .le. 0) go to 40
         iflag = 0
         if (mod(iter-1,nprint) .eq. 0)
     *      call ftanh2(m,n,xstp,ystp,x,fvec,fjac,ldfjac,iflag,NLUPRI)
         if (iflag .lt. 0) go to 300
   40    continue
c
c        compute the qr factorization of the jacobian.
c
         call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
c
c        on the first iteration and if mode is 1, scale according
c        to the norms of the columns of the initial jacobian.
c
         if (iter .ne. 1) go to 80
         if (mode .eq. 2) go to 60
         do 50 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   50       continue
   60    continue
c
c        on the first iteration, calculate the norm of the scaled x
c        and initialize the step bound delta.
c
         do 70 j = 1, n
            wa3(j) = diag(j)*x(j)
   70       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   80    continue
c
c        form (q transpose)*fvec and store the first n components in
c        qtf.
c
         do 90 i = 1, m
            wa4(i) = fvec(i)
   90       continue
         do 130 j = 1, n
            if (fjac(j,j) .eq. zero) go to 120
            sum = zero
            do 100 i = j, m
               sum = sum + fjac(i,j)*wa4(i)
  100          continue
            temp = -sum/fjac(j,j)
            do 110 i = j, m
               wa4(i) = wa4(i) + fjac(i,j)*temp
  110          continue
  120       continue
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
  130       continue
c
c        compute the norm of the scaled gradient.
c
         gnorm = zero
         if (fnorm .eq. zero) go to 170
         do 160 j = 1, n
            l = ipvt(j)
            if (wa2(l) .eq. zero) go to 150
            sum = zero
            do 140 i = 1, j
               sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  140          continue
            gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
  150       continue
  160       continue
  170    continue
c
c        test for convergence of the gradient norm.
c
         if (gnorm .le. gtol) info = 4
         if (info .ne. 0) go to 300
c
c        rescale if necessary.
c
         if (mode .eq. 2) go to 190
         do 180 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  180       continue
  190    continue
c
c        beginning of the inner loop.
c
  200    continue
c
c           determine the levenberg-marquardt parameter.
c
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,
     *                 wa3,wa4)
c
c           store the direction p and x + p. calculate the norm of p.
c
            do 210 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  210          continue
            pnorm = enorm(n,wa3)
c
c           on the first iteration, adjust the initial step bound.
c
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
c
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In LMDER 3: A, B & C :',wa2(1),wa2(2),wa2(3)
c
c           evaluate the function at x + p and calculate its norm.
c
            iflag = 1
            call ftanh2(m,n,xstp,ystp,wa2,wa4,fjac,ldfjac,iflag,NLUPRI)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(m,wa4)
c
c           compute the scaled actual reduction.
c
            actred = -one
            if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
c           compute the scaled predicted reduction and
c           the scaled directional derivative.
c
            do 230 j = 1, n
               wa3(j) = zero
               l = ipvt(j)
               temp = wa1(l)
               do 220 i = 1, j
                  wa3(i) = wa3(i) + fjac(i,j)*temp
  220             continue
  230          continue
            temp1 = enorm(n,wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
c
c           compute the ratio of the actual to the predicted
c           reduction.
c
            ratio = zero
            if (prered .ne. zero) ratio = actred/prered
c
c           update the step bound.
c
            if (ratio .gt. p25) go to 240
               if (actred .ge. zero) temp = p5
               if (actred .lt. zero)
     *            temp = p5*dirder/(dirder + p5*actred)
               if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
               delta = temp*dmin1(delta,pnorm/p1)
               par = par/temp
               go to 260
  240       continue
               if (par .ne. zero .and. ratio .lt. p75) go to 250
               delta = pnorm/p5
               par = p5*par
  250          continue
  260       continue
c
c           test for successful iteration.
c
            if (ratio .lt. p0001) go to 290
c
c           successful iteration. update x, fvec, and their norms.
c
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In LMDER 4a: A, B & C :',wa2(1),wa2(2),wa2(3)
c
            do 270 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
  270          continue
            do 280 i = 1, m
               fvec(i) = wa4(i)
  280          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  290       continue
c
c           tests for convergence.
c
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In LMDER 4b: A, B & C :',x(1),x(2),x(3)
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In LMDER 4c: A, B & C :',wa2(1),wa2(2),wa2(3)
c
            if (dabs(actred) .le. ftol .and. prered .le. ftol
     *          .and. p5*ratio .le. one) info = 1
            if (delta .le. xtol*xnorm) info = 2
            if (dabs(actred) .le. ftol .and. prered .le. ftol
     *          .and. p5*ratio .le. one .and. info .eq. 2) info = 3
c
CPFP here it does not work ...
c      do 310 j = 1, n
c               x(j) = wa2(j)
c  310          continue
Cend-PFP
c
            if (info .ne. 0) go to 300
c
c           tests for termination and stringent tolerances.
c
            if (nfev .ge. maxfev) info = 5
            if (dabs(actred) .le. epsmch .and. prered .le. epsmch
     *          .and. p5*ratio .le. one) info = 6
            if (delta .le. epsmch*xnorm) info = 7
            if (gnorm .le. epsmch) info = 8
            if (info .ne. 0) go to 300
c
c           end of the inner loop. repeat if iteration unsuccessful.
c
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In LMDER 5: A, B & C :',x(1),x(2),x(3)
c
            if (ratio .lt. p0001) go to 200
c
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In LMDER Fa: A, B & C :',x(1),x(2),x(3)
c
c        end of the outer loop.
c
         go to 30
  300 continue
c
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In LMDER Fb: A, B & C :',x(1),x(2),x(3)
c        WRITE(NLUPRI,'(/,1X,A,I3,/)')
c     &      ' In LMDER Fc: info :',info
c
c     termination, either normal or user imposed.
c
CPFP
c    For some reason this is not done before and is let to the user
c    together with some convergence condition ... !!!
c
c      do 310 j = 1, n
c               x(j) = wa2(j)
c  310          continue
Cend-PFP
c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) then
        NLUPRI = nprint
           call ftanh2(m,n,xstp,ystp,x,fvec,fjac,ldfjac,iflag,NLUPRI)
      end if
c
      return
c
c     last card of subroutine lmder.
c
      end
*
**********************************************************************
*  /* Deck DPMPAR */
*
      double precision function dpmpar(i)
      integer i
c     **********
c
c     Function dpmpar
c
c     This function provides double precision machine parameters
c     when the appropriate set of data statements is activated (by
c     removing the c from column 1) and all other data statements are
c     rendered inactive. Most of the parameter values were obtained
c     from the corresponding Bell Laboratories Port Library function.
c
c     The function statement is
c
c       double precision function dpmpar(i)
c
c     where
c
c       i is an integer input variable set to 1, 2, or 3 which
c         selects the desired machine parameter. If the machine has
c         t base b digits and its smallest and largest exponents are
c         emin and emax, respectively, then these parameters are
c
c         dpmpar(1) = b**(1 - t), the machine precision,
c
c         dpmpar(2) = b**(emin - 1), the smallest magnitude,
c
c         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c     Argonne National Laboratory. MINPACK Project. November 1996.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
c
c     **********
      integer mcheps(4)
      integer minmag(4)
      integer maxmag(4)
      double precision dmach(3)
      equivalence (dmach(1),mcheps(1))
      equivalence (dmach(2),minmag(1))
      equivalence (dmach(3),maxmag(1))
c
c     Machine constants for the IBM 360/370 series,
c     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
c     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
c
c     data mcheps(1),mcheps(2) / z34100000, z00000000 /
c     data minmag(1),minmag(2) / z00100000, z00000000 /
c     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
c
c     Machine constants for the Honeywell 600/6000 series.
c
c     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
c     data minmag(1),minmag(2) / o402400000000, o000000000000 /
c     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
c
c     Machine constants for the CDC 6000/7000 series.
c
c     data mcheps(1) / 15614000000000000000b /
c     data mcheps(2) / 15010000000000000000b /
c
c     data minmag(1) / 00604000000000000000b /
c     data minmag(2) / 00000000000000000000b /
c
c     data maxmag(1) / 37767777777777777777b /
c     data maxmag(2) / 37167777777777777777b /
c
c     Machine constants for the PDP-10 (KA processor).
c
c     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
c     data minmag(1),minmag(2) / "033400000000, "000000000000 /
c     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
c
c     Machine constants for the PDP-10 (KI processor).
c
c     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
c     data minmag(1),minmag(2) / "000400000000, "000000000000 /
c     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
c
c     Machine constants for the PDP-11. 
c
c     data mcheps(1),mcheps(2) /   9472,      0 /
c     data mcheps(3),mcheps(4) /      0,      0 /
c
c     data minmag(1),minmag(2) /    128,      0 /
c     data minmag(3),minmag(4) /      0,      0 /
c
c     data maxmag(1),maxmag(2) /  32767,     -1 /
c     data maxmag(3),maxmag(4) /     -1,     -1 /
c
c     Machine constants for the Burroughs 6700/7700 systems.
c
c     data mcheps(1) / o1451000000000000 /
c     data mcheps(2) / o0000000000000000 /
c
c     data minmag(1) / o1771000000000000 /
c     data minmag(2) / o7770000000000000 /
c
c     data maxmag(1) / o0777777777777777 /
c     data maxmag(2) / o7777777777777777 /
c
c     Machine constants for the Burroughs 5700 system.
c
c     data mcheps(1) / o1451000000000000 /
c     data mcheps(2) / o0000000000000000 /
c
c     data minmag(1) / o1771000000000000 /
c     data minmag(2) / o0000000000000000 /
c
c     data maxmag(1) / o0777777777777777 /
c     data maxmag(2) / o0007777777777777 /
c
c     Machine constants for the Burroughs 1700 system.
c
c     data mcheps(1) / zcc6800000 /
c     data mcheps(2) / z000000000 /
c
c     data minmag(1) / zc00800000 /
c     data minmag(2) / z000000000 /
c
c     data maxmag(1) / zdffffffff /
c     data maxmag(2) / zfffffffff /
c
c     Machine constants for the Univac 1100 series.
c
c     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
c     data minmag(1),minmag(2) / o000040000000, o000000000000 /
c     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
c
c     Machine constants for the Data General Eclipse S/200.
c
c     Note - it may be appropriate to include the following card -
c     static dmach(3)
c
c     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
c     data mcheps/32020k,3*0/
c
c     Machine constants for the Harris 220.
c
c     data mcheps(1),mcheps(2) / '20000000, '00000334 /
c     data minmag(1),minmag(2) / '20000000, '00000201 /
c     data maxmag(1),maxmag(2) / '37777777, '37777577 /
c
c     Machine constants for the Cray-1.
c
c     data mcheps(1) / 0376424000000000000000b /
c     data mcheps(2) / 0000000000000000000000b /
c
c     data minmag(1) / 0200034000000000000000b /
c     data minmag(2) / 0000000000000000000000b /
c
c     data maxmag(1) / 0577777777777777777777b /
c     data maxmag(2) / 0000007777777777777776b /
c
c     Machine constants for the Prime 400.
c
c     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
c     data minmag(1),minmag(2) / :10000000000, :00000100000 /
c     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
c
c     Machine constants for the VAX-11.
c
c     data mcheps(1),mcheps(2) /   9472,  0 /
c     data minmag(1),minmag(2) /    128,  0 /
c     data maxmag(1),maxmag(2) / -32769, -1 /
c
c     Machine constants for IEEE machines.
c
      data dmach(1) /2.22044604926d-16/
      data dmach(2) /2.22507385852d-308/
      data dmach(3) /1.79769313485d+308/
c
      dpmpar = dmach(i)
      return
c
c     Last card of function dpmpar.
c
      end
*
**********************************************************************
*  /*Deck ENORM */
*
      double precision function enorm(n,x)
      integer n
      double precision x(n)
c     **********
c
c     function enorm
c
c     given an n-vector x, this function calculates the
c     euclidean norm of x.
c
c     the euclidean norm is computed by accumulating the sum of
c     squares in three different sums. the sums of squares for the
c     small and large components are scaled so that no overflows
c     occur. non-destructive underflows are permitted. underflows
c     and overflows do not occur in the computation of the unscaled
c     sum of squares for the intermediate components.
c     the definitions of small, intermediate and large components
c     depend on two constants, rdwarf and rgiant. the main
c     restrictions on these constants are that rdwarf**2 not
c     underflow and rgiant**2 not overflow. the constants
c     given here are suitable for every known computer.
c
c     the function statement is
c
c       double precision function enorm(n,x)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c     subprograms called
c
c       fortran-supplied ... dabs,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,
     *                 x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do 90 i = 1, n
         xabs = dabs(x(i))
         if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70
            if (xabs .le. rdwarf) go to 30
c
c              sum for large components.
c
               if (xabs .le. x1max) go to 10
                  s1 = one + s1*(x1max/xabs)**2
                  x1max = xabs
                  go to 20
   10          continue
                  s1 = s1 + (xabs/x1max)**2
   20          continue
               go to 60
   30       continue
c
c              sum for small components.
c
               if (xabs .le. x3max) go to 40
                  s3 = one + s3*(x3max/xabs)**2
                  x3max = xabs
                  go to 50
   40          continue
                  if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
   50          continue
   60       continue
            go to 80
   70    continue
c
c           sum for intermediate components.
c
            s2 = s2 + xabs**2
   80    continue
   90    continue
c
c     calculation of norm.
c
      if (s1 .eq. zero) go to 100
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
         go to 130
  100 continue
         if (s2 .eq. zero) go to 110
            if (s2 .ge. x3max)
     *         enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
            if (s2 .lt. x3max)
     *         enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
            go to 120
  110    continue
            enorm = x3max*dsqrt(s3)
  120    continue
  130 continue
      return
c
c     last card of function enorm.
c
      end
*
**********************************************************************
*  /* Deck  LMPAR */
*
      subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,
     *                 wa2)
      integer n,ldr
      integer ipvt(n)
      double precision delta,par
      double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n),
     *                 wa2(n)
c     **********
c
c     subroutine lmpar
c
c     given an m by n matrix a, an n by n nonsingular diagonal
c     matrix d, an m-vector b, and a positive number delta,
c     the problem is to determine a value for the parameter
c     par such that if x solves the system
c
c           a*x = b ,     sqrt(par)*d*x = 0 ,
c
c     in the least squares sense, and dxnorm is the euclidean
c     norm of d*x, then either par is zero and
c
c           (dxnorm-delta) .le. 0.1*delta ,
c
c     or par is positive and
c
c           abs(dxnorm-delta) .le. 0.1*delta .
c
c     this subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     qr factorization, with column pivoting, of a. that is, if
c     a*p = q*r, where p is a permutation matrix, q has orthogonal
c     columns, and r is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, then lmpar expects
c     the full upper triangle of r, the permutation matrix p,
c     and the first n components of (q transpose)*b. on output
c     lmpar also provides an upper triangular matrix s such that
c
c            t   t                   t
c           p *(a *a + par*d*d)*p = s *s .
c
c     s is employed within lmpar and may be of separate interest.
c
c     only a few iterations are generally needed for convergence
c     of the algorithm. if, however, the limit of 10 iterations
c     is reached, then the output par will contain the best
c     value obtained so far.
c
c     the subroutine statement is
c
c       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
c                        wa1,wa2)
c
c     where
c
c       n is a positive integer input variable set to the order of r.
c
c       r is an n by n array. on input the full upper triangle
c         must contain the full upper triangle of the matrix r.
c         on output the full upper triangle is unaltered, and the
c         strict lower triangle contains the strict upper triangle
c         (transposed) of the upper triangular matrix s.
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       ipvt is an integer input array of length n which defines the
c         permutation matrix p such that a*p = q*r. column j of p
c         is column ipvt(j) of the identity matrix.
c
c       diag is an input array of length n which must contain the
c         diagonal elements of the matrix d.
c
c       qtb is an input array of length n which must contain the first
c         n elements of the vector (q transpose)*b.
c
c       delta is a positive input variable which specifies an upper
c         bound on the euclidean norm of d*x.
c
c       par is a nonnegative variable. on input par contains an
c         initial estimate of the levenberg-marquardt parameter.
c         on output par contains the final estimate.
c
c       x is an output array of length n which contains the least
c         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
c         for the output par.
c
c       sdiag is an output array of length n which contains the
c         diagonal elements of the upper triangular matrix s.
c
c       wa1 and wa2 are work arrays of length n.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar,enorm,qrsolv
c
c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,iter,j,jm1,jp1,k,l,nsing
      double precision dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001,
     *                 sum,temp,zero
      double precision dpmpar,enorm
      data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/
c
c     dwarf is the smallest positive magnitude.
c
      dwarf = dpmpar(2)
c
c     compute and store in x the gauss-newton direction. if the
c     jacobian is rank-deficient, obtain a least squares solution.
c
      nsing = n
      do 10 j = 1, n
         wa1(j) = qtb(j)
         if (r(j,j) .eq. zero .and. nsing .eq. n) nsing = j - 1
         if (nsing .lt. n) wa1(j) = zero
   10    continue
      if (nsing .lt. 1) go to 50
      do 40 k = 1, nsing
         j = nsing - k + 1
         wa1(j) = wa1(j)/r(j,j)
         temp = wa1(j)
         jm1 = j - 1
         if (jm1 .lt. 1) go to 30
         do 20 i = 1, jm1
            wa1(i) = wa1(i) - r(i,j)*temp
   20       continue
   30    continue
   40    continue
   50 continue
      do 60 j = 1, n
         l = ipvt(j)
         x(l) = wa1(j)
   60    continue
c
c     initialize the iteration counter.
c     evaluate the function at the origin, and test
c     for acceptance of the gauss-newton direction.
c
      iter = 0
      do 70 j = 1, n
         wa2(j) = diag(j)*x(j)
   70    continue
      dxnorm = enorm(n,wa2)
      fp = dxnorm - delta
      if (fp .le. p1*delta) go to 220
c
c     if the jacobian is not rank deficient, the newton
c     step provides a lower bound, parl, for the zero of
c     the function. otherwise set this bound to zero.
c
      parl = zero
      if (nsing .lt. n) go to 120
      do 80 j = 1, n
         l = ipvt(j)
         wa1(j) = diag(l)*(wa2(l)/dxnorm)
   80    continue
      do 110 j = 1, n
         sum = zero
         jm1 = j - 1
         if (jm1 .lt. 1) go to 100
         do 90 i = 1, jm1
            sum = sum + r(i,j)*wa1(i)
   90       continue
  100    continue
         wa1(j) = (wa1(j) - sum)/r(j,j)
  110    continue
      temp = enorm(n,wa1)
      parl = ((fp/delta)/temp)/temp
  120 continue
c
c     calculate an upper bound, paru, for the zero of the function.
c
      do 140 j = 1, n
         sum = zero
         do 130 i = 1, j
            sum = sum + r(i,j)*qtb(i)
  130       continue
         l = ipvt(j)
         wa1(j) = sum/diag(l)
  140    continue
      gnorm = enorm(n,wa1)
      paru = gnorm/delta
      if (paru .eq. zero) paru = dwarf/dmin1(delta,p1)
c
c     if the input par lies outside of the interval (parl,paru),
c     set par to the closer endpoint.
c
      par = dmax1(par,parl)
      par = dmin1(par,paru)
      if (par .eq. zero) par = gnorm/dxnorm
c
c     beginning of an iteration.
c
  150 continue
         iter = iter + 1
c
c        evaluate the function at the current value of par.
c
         if (par .eq. zero) par = dmax1(dwarf,p001*paru)
         temp = dsqrt(par)
         do 160 j = 1, n
            wa1(j) = temp*diag(j)
  160       continue
         call qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2)
         do 170 j = 1, n
            wa2(j) = diag(j)*x(j)
  170       continue
         dxnorm = enorm(n,wa2)
         temp = fp
         fp = dxnorm - delta
c
c        if the function is small enough, accept the current value
c        of par. also test for the exceptional cases where parl
c        is zero or the number of iterations has reached 10.
c
         if (dabs(fp) .le. p1*delta
     *       .or. parl .eq. zero .and. fp .le. temp
     *            .and. temp .lt. zero .or. iter .eq. 10) go to 220
c
c        compute the newton correction.
c
         do 180 j = 1, n
            l = ipvt(j)
            wa1(j) = diag(l)*(wa2(l)/dxnorm)
  180       continue
         do 210 j = 1, n
            wa1(j) = wa1(j)/sdiag(j)
            temp = wa1(j)
            jp1 = j + 1
            if (n .lt. jp1) go to 200
            do 190 i = jp1, n
               wa1(i) = wa1(i) - r(i,j)*temp
  190          continue
  200       continue
  210       continue
         temp = enorm(n,wa1)
         parc = ((fp/delta)/temp)/temp
c
c        depending on the sign of the function, update parl or paru.
c
         if (fp .gt. zero) parl = dmax1(parl,par)
         if (fp .lt. zero) paru = dmin1(paru,par)
c
c        compute an improved estimate for par.
c
         par = dmax1(parl,par+parc)
c
c        end of an iteration.
c
         go to 150
  220 continue
c
c     termination.
c
      if (iter .eq. 0) par = zero
      return
c
c     last card of subroutine lmpar.
c
      end
*
**********************************************************************
*  /* Deck  QRFAC */
*
      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      integer m,n,lda,lipvt
      integer ipvt(lipvt)
      logical pivot
      double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
c     **********
c
c     subroutine qrfac
c
c     this subroutine uses householder transformations with column
c     pivoting (optional) to compute a qr factorization of the
c     m by n matrix a. that is, qrfac determines an orthogonal
c     matrix q, a permutation matrix p, and an upper trapezoidal
c     matrix r with diagonal elements of nonincreasing magnitude,
c     such that a*p = q*r. the householder transformation for
c     column k, k = 1,2,...,min(m,n), is of the form
c
c                           t
c           i - (1/u(k))*u*u
c
c     where u has zeros in the first k-1 positions. the form of
c     this transformation and the method of pivoting first
c     appeared in the corresponding linpack subroutine.
c
c     the subroutine statement is
c
c       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       a is an m by n array. on input a contains the matrix for
c         which the qr factorization is to be computed. on output
c         the strict upper trapezoidal part of a contains the strict
c         upper trapezoidal part of r, and the lower trapezoidal
c         part of a contains a factored form of q (the non-trivial
c         elements of the u vectors described above).
c
c       lda is a positive integer input variable not less than m
c         which specifies the leading dimension of the array a.
c
c       pivot is a logical input variable. if pivot is set true,
c         then column pivoting is enforced. if pivot is set false,
c         then no column pivoting is done.
c
c       ipvt is an integer output array of length lipvt. ipvt
c         defines the permutation matrix p such that a*p = q*r.
c         column j of p is column ipvt(j) of the identity matrix.
c         if pivot is false, ipvt is not referenced.
c
c       lipvt is a positive integer input variable. if pivot is false,
c         then lipvt may be as small as 1. if pivot is true, then
c         lipvt must be at least n.
c
c       rdiag is an output array of length n which contains the
c         diagonal elements of r.
c
c       acnorm is an output array of length n which contains the
c         norms of the corresponding columns of the input matrix a.
c         if this information is not needed, then acnorm can coincide
c         with rdiag.
c
c       wa is a work array of length n. if pivot is false, then wa
c         can coincide with rdiag.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar,enorm
c
c       fortran-supplied ... dmax1,dsqrt,min0
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,j,jp1,k,kmax,minmn
      double precision ajnorm,epsmch,one,p05,sum,temp,zero
      double precision dpmpar,enorm
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
c     compute the initial column norms and initialize several arrays.
c
      do 10 j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
   10    continue
c
c     reduce a to r with householder transformations.
c
      minmn = min0(m,n)
      do 110 j = 1, minmn
         if (.not.pivot) go to 40
c
c        bring the column of largest norm into the pivot position.
c
         kmax = j
         do 20 k = j, n
            if (rdiag(k) .gt. rdiag(kmax)) kmax = k
   20       continue
         if (kmax .eq. j) go to 40
         do 30 i = 1, m
            temp = a(i,j)
            a(i,j) = a(i,kmax)
            a(i,kmax) = temp
   30       continue
         rdiag(kmax) = rdiag(j)
         wa(kmax) = wa(j)
         k = ipvt(j)
         ipvt(j) = ipvt(kmax)
         ipvt(kmax) = k
   40    continue
c
c        compute the householder transformation to reduce the
c        j-th column of a to a multiple of the j-th unit vector.
c
         ajnorm = enorm(m-j+1,a(j,j))
         if (ajnorm .eq. zero) go to 100
         if (a(j,j) .lt. zero) ajnorm = -ajnorm
         do 50 i = j, m
            a(i,j) = a(i,j)/ajnorm
   50       continue
         a(j,j) = a(j,j) + one
c
c        apply the transformation to the remaining columns
c        and update the norms.
c
         jp1 = j + 1
         if (n .lt. jp1) go to 100
         do 90 k = jp1, n
            sum = zero
            do 60 i = j, m
               sum = sum + a(i,j)*a(i,k)
   60          continue
            temp = sum/a(j,j)
            do 70 i = j, m
               a(i,k) = a(i,k) - temp*a(i,j)
   70          continue
            if (.not.pivot .or. rdiag(k) .eq. zero) go to 80
            temp = a(j,k)/rdiag(k)
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
            if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80
            rdiag(k) = enorm(m-j,a(jp1,k))
            wa(k) = rdiag(k)
   80       continue
   90       continue
  100    continue
         rdiag(j) = -ajnorm
  110    continue
      return
c
c     last card of subroutine qrfac.
c
      end
*
**********************************************************************
*  /* Deck QRSOLV */
*
      subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
      integer n,ldr
      integer ipvt(n)
      double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
c     **********
c
c     subroutine qrsolv
c
c     given an m by n matrix a, an n by n diagonal matrix d,
c     and an m-vector b, the problem is to determine an x which
c     solves the system
c
c           a*x = b ,     d*x = 0 ,
c
c     in the least squares sense.
c
c     this subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     qr factorization, with column pivoting, of a. that is, if
c     a*p = q*r, where p is a permutation matrix, q has orthogonal
c     columns, and r is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, then qrsolv expects
c     the full upper triangle of r, the permutation matrix p,
c     and the first n components of (q transpose)*b. the system
c     a*x = b, d*x = 0, is then equivalent to
c
c                  t       t
c           r*z = q *b ,  p *d*p*z = 0 ,
c
c     where x = p*z. if this system does not have full rank,
c     then a least squares solution is obtained. on output qrsolv
c     also provides an upper triangular matrix s such that
c
c            t   t               t
c           p *(a *a + d*d)*p = s *s .
c
c     s is computed within qrsolv and may be of separate interest.
c
c     the subroutine statement is
c
c       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
c
c     where
c
c       n is a positive integer input variable set to the order of r.
c
c       r is an n by n array. on input the full upper triangle
c         must contain the full upper triangle of the matrix r.
c         on output the full upper triangle is unaltered, and the
c         strict lower triangle contains the strict upper triangle
c         (transposed) of the upper triangular matrix s.
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       ipvt is an integer input array of length n which defines the
c         permutation matrix p such that a*p = q*r. column j of p
c         is column ipvt(j) of the identity matrix.
c
c       diag is an input array of length n which must contain the
c         diagonal elements of the matrix d.
c
c       qtb is an input array of length n which must contain the first
c         n elements of the vector (q transpose)*b.
c
c       x is an output array of length n which contains the least
c         squares solution of the system a*x = b, d*x = 0.
c
c       sdiag is an output array of length n which contains the
c         diagonal elements of the upper triangular matrix s.
c
c       wa is a work array of length n.
c
c     subprograms called
c
c       fortran-supplied ... dabs,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,j,jp1,k,kp1,l,nsing
      double precision cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero
      data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/
c
c     copy r and (q transpose)*b to preserve input and initialize s.
c     in particular, save the diagonal elements of r in x.
c
      do 20 j = 1, n
         do 10 i = j, n
            r(i,j) = r(j,i)
   10       continue
         x(j) = r(j,j)
         wa(j) = qtb(j)
   20    continue
c
c     eliminate the diagonal matrix d using a givens rotation.
c
      do 100 j = 1, n
c
c        prepare the row of d to be eliminated, locating the
c        diagonal element using p from the qr factorization.
c
         l = ipvt(j)
         if (diag(l) .eq. zero) go to 90
         do 30 k = j, n
            sdiag(k) = zero
   30       continue
         sdiag(j) = diag(l)
c
c        the transformations to eliminate the row of d
c        modify only a single element of (q transpose)*b
c        beyond the first n, which is initially zero.
c
         qtbpj = zero
         do 80 k = j, n
c
c           determine a givens rotation which eliminates the
c           appropriate element in the current row of d.
c
            if (sdiag(k) .eq. zero) go to 70
            if (dabs(r(k,k)) .ge. dabs(sdiag(k))) go to 40
               cotan = r(k,k)/sdiag(k)
               sin = p5/dsqrt(p25+p25*cotan**2)
               cos = sin*cotan
               go to 50
   40       continue
               tan = sdiag(k)/r(k,k)
               cos = p5/dsqrt(p25+p25*tan**2)
               sin = cos*tan
   50       continue
c
c           compute the modified diagonal element of r and
c           the modified element of ((q transpose)*b,0).
c
            r(k,k) = cos*r(k,k) + sin*sdiag(k)
            temp = cos*wa(k) + sin*qtbpj
            qtbpj = -sin*wa(k) + cos*qtbpj
            wa(k) = temp
c
c           accumulate the tranformation in the row of s.
c
            kp1 = k + 1
            if (n .lt. kp1) go to 70
            do 60 i = kp1, n
               temp = cos*r(i,k) + sin*sdiag(i)
               sdiag(i) = -sin*r(i,k) + cos*sdiag(i)
               r(i,k) = temp
   60          continue
   70       continue
   80       continue
   90    continue
c
c        store the diagonal element of s and restore
c        the corresponding diagonal element of r.
c
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
  100    continue
c
c     solve the triangular system for z. if the system is
c     singular, then obtain a least squares solution.
c
      nsing = n
      do 110 j = 1, n
         if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1
         if (nsing .lt. n) wa(j) = zero
  110    continue
      if (nsing .lt. 1) go to 150
      do 140 k = 1, nsing
         j = nsing - k + 1
         sum = zero
         jp1 = j + 1
         if (nsing .lt. jp1) go to 130
         do 120 i = jp1, nsing
            sum = sum + r(i,j)*wa(i)
  120       continue
  130    continue
         wa(j) = (wa(j) - sum)/sdiag(j)
  140    continue
  150 continue
c
c     permute the components of z back to components of x.
c
      do 160 j = 1, n
         l = ipvt(j)
         x(l) = wa(j)
  160    continue
      return
c
c     last card of subroutine qrsolv.
c
      end
*
**********************************************************************
*  /* Deck LMDIF1 */
*
      subroutine lmdif1(ftanh3,m,n,xstp,ystp,x,fvec,tol,info,iwa,
     *                  wa,lwa,NLUPRI)
      integer m,n,info,lwa,NLUPRI
      integer iwa(n)
      double precision tol
      double precision x(n),fvec(m),wa(lwa),xstp(m),ystp(m)
      external ftanh3
c     **********
c
c     subroutine lmdif1
c
c     the purpose of lmdif1 is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of the
c     levenberg-marquardt algorithm. this is done by using the more
c     general least-squares solver lmdif. the user must provide a
c     subroutine which calculates the functions. the jacobian is
c     then calculated by a forward-difference approximation.
c
c     the subroutine statement is
c
c       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,iflag)
c         integer m,n,iflag
c         double precision x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmdif1.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
c       tol is a nonnegative input variable. termination occurs
c         when the algorithm estimates either that the relative
c         error in the sum of squares is at most tol or that
c         the relative error between x and the solution is at
c         most tol.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  algorithm estimates that the relative error
c                   in the sum of squares is at most tol.
c
c         info = 2  algorithm estimates that the relative error
c                   between x and the solution is at most tol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  fvec is orthogonal to the columns of the
c                   jacobian to machine precision.
c
c         info = 5  number of calls to fcn has reached or
c                   exceeded 200*(n+1).
c
c         info = 6  tol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  tol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c       iwa is an integer work array of length n.
c
c       wa is a work array of length lwa.
c
c       lwa is a positive integer input variable not less than
c         m*n+5*n+m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... lmdif
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer maxfev,mode,mp5n,nfev,nprint
      double precision epsfcn,factor,ftol,gtol,xtol,zero
      data factor,zero /1.0d2,0.0d0/
      info = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. m .lt. n .or. tol .lt. zero
     *    .or. lwa .lt. m*n + 5*n + m) go to 10
c
c     call lmdif.
c
      maxfev = 200*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      epsfcn = zero
      mode = 1
      nprint = 0
      mp5n = m + 5*n
      call lmdif(ftanh3,m,n,xstp,ystp,x,fvec,ftol,xtol,gtol,
     *           maxfev,epsfcn,wa(1),mode,factor,nprint,info,
     *           nfev,wa(mp5n+1),m,iwa,wa(n+1),wa(2*n+1),
     *           wa(3*n+1),wa(4*n+1),wa(5*n+1),NLUPRI)
      if (info .eq. 8) info = 4
   10 continue
      return
c
c     last card of subroutine lmdif1.
c
      end
*
**********************************************************************
*  /* Deck LMDIF */
*
      subroutine lmdif(ftanh3,m,n,xstp,ystp,x,fvec,ftol,xtol,gtol,
     *                 maxfev,epsfcn,diag,mode,factor,nprint,info,
     *                 nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4,
     *                 NLUPRI)
      integer m,n,maxfev,mode,nprint,info,nfev,ldfjac,NLUPRI
      integer ipvt(n)
      double precision ftol,xtol,gtol,epsfcn,factor,xstp(m),ystp(m)
      double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n),
     *                 wa1(n),wa2(n),wa3(n),wa4(m)
      external ftanh3
c     **********
c
c     subroutine lmdif
c
c     the purpose of lmdif is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of
c     the levenberg-marquardt algorithm. the user must provide a
c     subroutine which calculates the functions. the jacobian is
c     then calculated by a forward-difference approximation.
c
c     the subroutine statement is
c
c       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
c                        diag,mode,factor,nprint,info,nfev,fjac,
c                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,iflag)
c         integer m,n,iflag
c         double precision x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmdif.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
c       ftol is a nonnegative input variable. termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most ftol.
c         therefore, ftol measures the relative error desired
c         in the sum of squares.
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol. therefore, xtol measures the
c         relative error desired in the approximate solution.
c
c       gtol is a nonnegative input variable. termination
c         occurs when the cosine of the angle between fvec and
c         any column of the jacobian is at most gtol in absolute
c         value. therefore, gtol measures the orthogonality
c         desired between the function vector and the columns
c         of the jacobian.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn is at least
c         maxfev by the end of an iteration.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       diag is an array of length n. if mode = 1 (see
c         below), diag is internally set. if mode = 2, diag
c         must contain positive entries that serve as
c         multiplicative scale factors for the variables.
c
c       mode is an integer input variable. if mode = 1, the
c         variables will be scaled internally. if mode = 2,
c         the scaling is specified by the input diag. other
c         values of mode are equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.). 100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x and fvec available
c         for printing. if nprint is not positive, no special calls
c         of fcn with iflag = 0 are made.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  both actual and predicted relative reductions
c                   in the sum of squares are at most ftol.
c
c         info = 2  relative error between two consecutive iterates
c                   is at most xtol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  the cosine of the angle between fvec and any
c                   column of the jacobian is at most gtol in
c                   absolute value.
c
c         info = 5  number of calls to fcn has reached or
c                   exceeded maxfev.
c
c         info = 6  ftol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  xtol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c         info = 8  gtol is too small. fvec is orthogonal to the
c                   columns of the jacobian to machine precision.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn.
c
c       fjac is an output m by n array. the upper n by n submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       ipvt is an integer output array of length n. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       qtf is an output array of length n which contains
c         the first n elements of the vector (q transpose)*fvec.
c
c       wa1, wa2, and wa3 are work arrays of length n.
c
c       wa4 is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
c
c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,iflag,iter,j,l
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,
     *                 one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,
     *                 sum,temp,temp1,temp2,xnorm,zero
      double precision dpmpar,enorm
      data one,p1,p5,p25,p75,p0001,zero
     *     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      info = 0
      iflag = 0
      nfev = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m
     *    .or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero
     *    .or. maxfev .le. 0 .or. factor .le. zero) go to 300
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
c
c     evaluate the function at the starting point
c     and calculate its norm.
c
      iflag = 1
      call ftanh3(m,n,xstp,ystp,x,fvec,iflag,NLUPRI)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(m,fvec)
c
c     initialize levenberg-marquardt parameter and iteration counter.
c
      par = zero
      iter = 1
c
c     beginning of the outer loop.
c
   30 continue
c
c        calculate the jacobian matrix.
c
         iflag = 2
         call fdjac2(ftanh3,m,n,xstp,ystp,x,fvec,fjac,ldfjac,iflag,
     *               epsfcn,wa4,NLUPRI)
         nfev = nfev + n
         if (iflag .lt. 0) go to 300
c
c        if requested, call fcn to enable printing of iterates.
c
         if (nprint .le. 0) go to 40
         iflag = 0
         if (mod(iter-1,nprint) .eq. 0) 
     *              call ftanh3(m,n,xstp,ystp,x,fvec,iflag,NLUPRI)
         if (iflag .lt. 0) go to 300
   40    continue
c
c        compute the qr factorization of the jacobian.
c
         call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
c
c        on the first iteration and if mode is 1, scale according
c        to the norms of the columns of the initial jacobian.
c
         if (iter .ne. 1) go to 80
         if (mode .eq. 2) go to 60
         do 50 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   50       continue
   60    continue
c
c        on the first iteration, calculate the norm of the scaled x
c        and initialize the step bound delta.
c
         do 70 j = 1, n
            wa3(j) = diag(j)*x(j)
   70       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   80    continue
c
c        form (q transpose)*fvec and store the first n components in
c        qtf.
c
         do 90 i = 1, m
            wa4(i) = fvec(i)
   90       continue
         do 130 j = 1, n
            if (fjac(j,j) .eq. zero) go to 120
            sum = zero
            do 100 i = j, m
               sum = sum + fjac(i,j)*wa4(i)
  100          continue
            temp = -sum/fjac(j,j)
            do 110 i = j, m
               wa4(i) = wa4(i) + fjac(i,j)*temp
  110          continue
  120       continue
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
  130       continue
c
c        compute the norm of the scaled gradient.
c
         gnorm = zero
         if (fnorm .eq. zero) go to 170
         do 160 j = 1, n
            l = ipvt(j)
            if (wa2(l) .eq. zero) go to 150
            sum = zero
            do 140 i = 1, j
               sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  140          continue
            gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
  150       continue
  160       continue
  170    continue
c
c        test for convergence of the gradient norm.
c
         if (gnorm .le. gtol) info = 4
         if (info .ne. 0) go to 300
c
c        rescale if necessary.
c
         if (mode .eq. 2) go to 190
         do 180 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  180       continue
  190    continue
c
c        beginning of the inner loop.
c
  200    continue
c
c           determine the levenberg-marquardt parameter.
c
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,
     *                 wa3,wa4)
c
c           store the direction p and x + p. calculate the norm of p.
c
            do 210 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  210          continue
            pnorm = enorm(n,wa3)
c
c           on the first iteration, adjust the initial step bound.
c
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
c
c           evaluate the function at x + p and calculate its norm.
c
            iflag = 1
            call ftanh3(m,n,xstp,ystp,wa2,wa4,iflag,NLUPRI)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(m,wa4)
c
c           compute the scaled actual reduction.
c
            actred = -one
            if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
c           compute the scaled predicted reduction and
c           the scaled directional derivative.
c
            do 230 j = 1, n
               wa3(j) = zero
               l = ipvt(j)
               temp = wa1(l)
               do 220 i = 1, j
                  wa3(i) = wa3(i) + fjac(i,j)*temp
  220             continue
  230          continue
            temp1 = enorm(n,wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
c
c           compute the ratio of the actual to the predicted
c           reduction.
c
            ratio = zero
            if (prered .ne. zero) ratio = actred/prered
c
c           update the step bound.
c
            if (ratio .gt. p25) go to 240
               if (actred .ge. zero) temp = p5
               if (actred .lt. zero)
     *            temp = p5*dirder/(dirder + p5*actred)
               if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
               delta = temp*dmin1(delta,pnorm/p1)
               par = par/temp
               go to 260
  240       continue
               if (par .ne. zero .and. ratio .lt. p75) go to 250
               delta = pnorm/p5
               par = p5*par
  250          continue
  260       continue
c
c           test for successful iteration.
c
            if (ratio .lt. p0001) go to 290
c
c           successful iteration. update x, fvec, and their norms.
c
            do 270 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
  270          continue
            do 280 i = 1, m
               fvec(i) = wa4(i)
  280          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  290       continue
c
c           tests for convergence.
c
            if (dabs(actred) .le. ftol .and. prered .le. ftol
     *          .and. p5*ratio .le. one) info = 1
            if (delta .le. xtol*xnorm) info = 2
            if (dabs(actred) .le. ftol .and. prered .le. ftol
     *          .and. p5*ratio .le. one .and. info .eq. 2) info = 3
            if (info .ne. 0) go to 300
c
c           tests for termination and stringent tolerances.
c
            if (nfev .ge. maxfev) info = 5
            if (dabs(actred) .le. epsmch .and. prered .le. epsmch
     *          .and. p5*ratio .le. one) info = 6
            if (delta .le. epsmch*xnorm) info = 7
            if (gnorm .le. epsmch) info = 8
            if (info .ne. 0) go to 300
c
c           end of the inner loop. repeat if iteration unsuccessful.
c
            if (ratio .lt. p0001) go to 200
c
c        end of the outer loop.
c
         go to 30
  300 continue
c
c     termination, either normal or user imposed.
c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0)
     &    call ftanh3(m,n,xstp,ystp,x,fvec,iflag,NLUPRI)
      return
c
c     last card of subroutine lmdif.
c
      end
*
**********************************************************************
*  /* Deck FDJAC2 */
*
      subroutine fdjac2(ftanh3,m,n,xstp,ystp,x,fvec,fjac,ldfjac,iflag,
     *                  epsfcn,wa,NLUPRI)
      integer m,n,ldfjac,iflag, NLUPRI
      double precision epsfcn
      double precision x(n),fvec(m),fjac(ldfjac,n),wa(m),
     *                 xstp(m),ystp(m)
      external ftanh3
c     **********
c
c     subroutine fdjac2
c
c     this subroutine computes a forward-difference approximation
c     to the m by n jacobian matrix associated with a specified
c     problem of m functions in n variables.
c
c     the subroutine statement is
c
c       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,iflag)
c         integer m,n,iflag
c         double precision x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of fdjac2.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an input array of length n.
c
c       fvec is an input array of length m which must contain the
c         functions evaluated at x.
c
c       fjac is an output m by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       iflag is an integer variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dpmpar
c
c       fortran-supplied ... dabs,dmax1,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,j
      double precision eps,epsmch,h,temp,zero
      double precision dpmpar
      data zero /0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = dsqrt(dmax1(epsfcn,epsmch))
      do 20 j = 1, n
         temp = x(j)
         h = eps*dabs(temp)
         if (h .eq. zero) h = eps
         x(j) = temp + h
         call ftanh3(m,n,xstp,ystp,x,wa,iflag,NLUPRI)
         if (iflag .lt. 0) go to 30
         x(j) = temp
         do 10 i = 1, m
            fjac(i,j) = (wa(i) - fvec(i))/h
   10       continue
   20    continue
   30 continue
      return
c
c     last card of subroutine fdjac2.
c
      end
*
**********************************************************************
* PFP 2017: TAKEN FROM ASSESS (NL2SOL) and slightly modified
**********************************************************************
*
      SUBROUTINE ASSESS (D, IV, P, STEP, STLSTG, V, X, X0)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  ASSESS CANDIDATE STEP (NL2SOL VERSION 2.2)  ***
C
C  ***  PURPOSE  ***
C
C        THIS SUBROUTINE IS CALLED BY AN UNCONSTRAINED MINIMIZATION
C     ROUTINE TO ASSESS THE NEXT CANDIDATE STEP.  IT MAY RECOMMEND ONE
C     OF SEVERAL COURSES OF ACTION, SUCH AS ACCEPTING THE STEP, RECOM-
C     PUTING IT USING THE SAME OR A NEW QUADRATIC MODEL, OR HALTING DUE
C     TO CONVERGENCE OR FALSE CONVERGENCE.  SEE THE RETURN CODE LISTING
C     BELOW.
C
      implicit none
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   D(P),STEP(P),STLSTG(P),V(*),X(P),X0(P)
      INTEGER
     +   IV(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   EMAX,GTS,HALF,ONE,RELDX1,RFAC1,TEMP,TWO,XMAX,ZERO
      INTEGER
     +   AFCTOL,DECFAC,DST0,DSTNRM,DSTSAV,F,F0,FDIF,FLSTGD,GTSLST,
     +   GTSTEP,I,INCFAC,IRC,LMAX0,MLSTGD,MODEL,NFC,NFCALL,NFGCAL,
     +   NREDUC,PLSTGD,PREDUC,RADFAC,RADINC,RDFCMN,RDFCMX,RELDX,
     +   RESTOR,RFCTOL,STAGE,STGLIM,STPPAR,SWITCH,TOOBIG,TUNER1,
     +   TUNER2,TUNER3,XCTOL,XFTOL,XIRC
      LOGICAL
     +   GOODX
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   D1MACH,RELDST
      EXTERNAL D1MACH,RELDST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL VCOPY
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX
C
C--------------------------  PARAMETER USAGE  --------------------------
C
C     IV (I/O) INTEGER PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION
C             BELOW OF IV VALUES REFERENCED.
C      D (IN)  SCALE VECTOR USED IN COMPUTING V(RELDX) -- SEE BELOW.
C      P (IN)  NUMBER OF PARAMETERS BEING OPTIMIZED.
C   STEP (I/O) ON INPUT, STEP IS THE STEP TO BE ASSESSED.  IT IS UN-
C             CHANGED ON OUTPUT UNLESS A PREVIOUS STEP ACHIEVED A
C             BETTER OBJECTIVE FUNCTION REDUCTION, IN WHICH CASE STLSTG
C             WILL HAVE BEEN COPIED TO STEP.
C STLSTG (I/O) WHEN ASSESS RECOMMENDS RECOMPUTING STEP EVEN THOUGH THE
C             CURRENT (OR A PREVIOUS) STEP YIELDS AN OBJECTIVE FUNC-
C             TION DECREASE, IT SAVES IN STLSTG THE STEP THAT GAVE THE
C             BEST FUNCTION REDUCTION SEEN SO FAR (IN THE CURRENT ITERA-
C             TION).  IF THE RECOMPUTED STEP YIELDS A LARGER FUNCTION
C             VALUE, THEN STEP IS RESTORED FROM STLSTG AND
C             X = X0 + STEP IS RECOMPUTED.
C      V (I/O) DOUBLE PRECISION PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION
C             BELOW OF V VALUES REFERENCED.
C      X (I/O) ON INPUT, X = X0 + STEP IS THE POINT AT WHICH THE OBJEC-
C             TIVE FUNCTION HAS JUST BEEN EVALUATED.  IF AN EARLIER
C             STEP YIELDED A BIGGER FUNCTION DECREASE, THEN X IS
C             RESTORED TO THE CORRESPONDING EARLIER VALUE.  OTHERWISE,
C             IF THE CURRENT STEP DOES NOT GIVE ANY FUNCTION DECREASE,
C             THEN X IS RESTORED TO X0.
C     X0 (IN)  INITIAL OBJECTIVE FUNCTION PARAMETER VECTOR (AT THE
C             START OF THE CURRENT ITERATION).
C
C  ***  IV VALUES REFERENCED  ***
C
C    IV(IRC) (I/O) ON INPUT FOR THE FIRST STEP TRIED IN A NEW ITERATION,
C             IV(IRC) SHOULD BE SET TO 3 OR 4 (THE VALUE TO WHICH IT IS
C             SET WHEN STEP IS DEFINITELY TO BE ACCEPTED).  ON INPUT
C             AFTER STEP HAS BEEN RECOMPUTED, IV(IRC) SHOULD BE
C             UNCHANGED SINCE THE PREVIOUS RETURN OF ASSESS.
C                ON OUTPUT, IV(IRC) IS A RETURN CODE HAVING ONE OF THE
C             FOLLOWING VALUES...
C                  1 = SWITCH MODELS OR TRY SMALLER STEP.
C                  2 = SWITCH MODELS OR ACCEPT STEP.
C                  3 = ACCEPT STEP AND DETERMINE V(RADFAC) BY GRADIENT
C                       TESTS.
C                  4 = ACCEPT STEP, V(RADFAC) HAS BEEN DETERMINED.
C                  5 = RECOMPUTE STEP (USING THE SAME MODEL).
C                  6 = RECOMPUTE STEP WITH RADIUS = V(LMAX0) BUT DO NOT
C                       EVAULATE THE OBJECTIVE FUNCTION.
C                  7 = X-CONVERGENCE (SEE V(XCTOL)).
C                  8 = RELATIVE FUNCTION CONVERGENCE (SEE V(RFCTOL)).
C                  9 = BOTH X- AND RELATIVE FUNCTION CONVERGENCE.
C                 10 = ABSOLUTE FUNCTION CONVERGENCE (SEE V(AFCTOL)).
C                 11 = SINGULAR CONVERGENCE (SEE V(LMAX0)).
C                 12 = FALSE CONVERGENCE (SEE V(XFTOL)).
C                 13 = IV(IRC) WAS OUT OF RANGE ON INPUT.
C             RETURN CODE I HAS PRECDENCE OVER I+1 FOR I = 9, 10, 11.
C IV(MLSTGD) (I/O) SAVED VALUE OF IV(MODEL).
C  IV(MODEL) (I/O) ON INPUT, IV(MODEL) SHOULD BE AN INTEGER IDENTIFYING
C             THE CURRENT QUADRATIC MODEL OF THE OBJECTIVE FUNCTION.
C             IF A PREVIOUS STEP YIELDED A BETTER FUNCTION REDUCTION,
C             THEN IV(MODEL) WILL BE SET TO IV(MLSTGD) ON OUTPUT.
C IV(NFCALL) (IN)  INVOCATION COUNT FOR THE OBJECTIVE FUNCTION.
C IV(NFGCAL) (I/O) VALUE OF IV(NFCALL) AT STEP THAT GAVE THE BIGGEST
C             FUNCTION REDUCTION THIS ITERATION.  IV(NFGCAL) REMAINS
C             UNCHANGED UNTIL A FUNCTION REDUCTION IS OBTAINED.
C IV(RADINC) (I/O) THE NUMBER OF RADIUS INCREASES (OR MINUS THE NUMBER
C             OF DECREASES) SO FAR THIS ITERATION.
C IV(RESTOR) (OUT) SET TO 0 UNLESS X AND V(F) HAVE BEEN RESTORED, IN
C             WHICH CASE ASSESS SETS IV(RESTOR) = 1.
C  IV(STAGE) (I/O) COUNT OF THE NUMBER OF MODELS TRIED SO FAR IN THE
C             CURRENT ITERATION.
C IV(STGLIM) (IN)  MAXIMUM NUMBER OF MODELS TO CONSIDER.
C IV(SWITCH) (OUT) SET TO 0 UNLESS A NEW MODEL IS BEING TRIED AND IT
C             GIVES A SMALLER FUNCTION VALUE THAN THE PREVIOUS MODEL,
C             IN WHICH CASE ASSESS SETS IV(SWITCH) = 1.
C IV(TOOBIG) (IN)  IS NONZERO IF STEP WAS TOO BIG (E.G. IF IT CAUSED
C             OVERFLOW).
C   IV(XIRC) (I/O) VALUE THAT IV(IRC) WOULD HAVE IN THE ABSENCE OF
C             CONVERGENCE, FALSE CONVERGENCE, AND OVERSIZED STEPS.
C
C  ***  V VALUES REFERENCED  ***
C
C V(AFCTOL) (IN)  ABSOLUTE FUNCTION CONVERGENCE TOLERANCE.  IF THE
C             ABSOLUTE VALUE OF THE CURRENT FUNCTION VALUE V(F) IS LESS
C             THAN V(AFCTOL), THEN ASSESS RETURNS WITH IV(IRC) = 10.
C V(DECFAC) (IN)  FACTOR BY WHICH TO DECREASE RADIUS WHEN IV(TOOBIG) IS
C             NONZERO.
C V(DSTNRM) (IN)  THE 2-NORM OF D*STEP.
C V(DSTSAV) (I/O) VALUE OF V(DSTNRM) ON SAVED STEP.
C   V(DST0) (IN)  THE 2-NORM OF D TIMES THE NEWTON STEP (WHEN DEFINED,
C             I.E., FOR V(NREDUC) .GE. 0).
C      V(F) (I/O) ON BOTH INPUT AND OUTPUT, V(F) IS THE OBJECTIVE FUNC-
C             TION VALUE AT X.  IF X IS RESTORED TO A PREVIOUS VALUE,
C             THEN V(F) IS RESTORED TO THE CORRESPONDING VALUE.
C   V(FDIF) (OUT) THE FUNCTION REDUCTION V(F0) - V(F) (FOR THE OUTPUT
C             VALUE OF V(F) IF AN EARLIER STEP GAVE A BIGGER FUNCTION
C             DECREASE, AND FOR THE INPUT VALUE OF V(F) OTHERWISE).
C V(FLSTGD) (I/O) SAVED VALUE OF V(F).
C     V(F0) (IN)  OBJECTIVE FUNCTION VALUE AT START OF ITERATION.
C V(GTSLST) (I/O) VALUE OF V(GTSTEP) ON SAVED STEP.
C V(GTSTEP) (IN)  INNER PRODUCT BETWEEN STEP AND GRADIENT.
C V(INCFAC) (IN)  MINIMUM FACTOR BY WHICH TO INCREASE RADIUS.
C  V(LMAX0) (IN)  MAXIMUM REASONABLE STEP SIZE (AND INITIAL STEP BOUND).
C             IF THE ACTUAL FUNCTION DECREASE IS NO MORE THAN TWICE
C             WHAT WAS PREDICTED, IF A RETURN WITH IV(IRC) = 7, 8, 9,
C             OR 10 DOES NOT OCCUR, IF V(DSTNRM) .GT. V(LMAX0), AND IF
C             V(PREDUC) .LE. V(RFCTOL) * ABS(V(F0)), THEN ASSESS
C             RETURNS WITH IV(IRC) = 11.  IF SO DOING APPEARS WORTHWHILE,
C             THEN ASSESS REPEATS THIS TEST WITH V(PREDUC) COMPUTED FOR
C             A STEP OF LENGTH V(LMAX0) (BY A RETURN WITH IV(IRC) = 6).
C V(NREDUC) (I/O)  FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR
C             NEWTON STEP.  IF ASSESS IS CALLED WITH IV(IRC) = 6, I.E.,
C             IF V(PREDUC) HAS BEEN COMPUTED WITH RADIUS = V(LMAX0) FOR
C             USE IN THE SINGULAR CONVERVENCE TEST, THEN V(NREDUC) IS
C             SET TO -V(PREDUC) BEFORE THE LATTER IS RESTORED.
C V(PLSTGD) (I/O) VALUE OF V(PREDUC) ON SAVED STEP.
C V(PREDUC) (I/O) FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR
C             CURRENT STEP.
C V(RADFAC) (OUT) FACTOR TO BE USED IN DETERMINING THE NEW RADIUS,
C             WHICH SHOULD BE V(RADFAC)*DST, WHERE  DST  IS EITHER THE
C             OUTPUT VALUE OF V(DSTNRM) OR THE 2-NORM OF
C             DIAG(NEWD)*STEP  FOR THE OUTPUT VALUE OF STEP AND THE
C             UPDATED VERSION, NEWD, OF THE SCALE VECTOR D.  FOR
C             IV(IRC) = 3, V(RADFAC) = 1.0 IS RETURNED.
C V(RDFCMN) (IN)  MINIMUM VALUE FOR V(RADFAC) IN TERMS OF THE INPUT
C             VALUE OF V(DSTNRM) -- SUGGESTED VALUE = 0.1.
C V(RDFCMX) (IN)  MAXIMUM VALUE FOR V(RADFAC) -- SUGGESTED VALUE = 4.0.
C  V(RELDX) (OUT) SCALED RELATIVE CHANGE IN X CAUSED BY STEP, COMPUTED
C             BY FUNCTION  RELDST  AS
C                 MAX (D(I)*ABS(X(I)-X0(I)), 1 .LE. I .LE. P) /
C                    MAX (D(I)*(ABS(X(I))+ABS(X0(I))), 1 .LE. I .LE. P).
C             IF AN ACCEPTABLE STEP IS RETURNED, THEN V(RELDX) IS COM-
C             PUTED USING THE OUTPUT (POSSIBLY RESTORED) VALUES OF X
C             AND STEP.  OTHERWISE IT IS COMPUTED USING THE INPUT
C             VALUES.
C V(RFCTOL) (IN)  RELATIVE FUNCTION CONVERGENCE TOLERANCE.  IF THE
C             ACTUAL FUNCTION REDUCTION IS AT MOST TWICE WHAT WAS 
C             PREDICTED AND  V(NREDUC) .LE. V(RFCTOL)*ABS(V(F0)),  THEN
C             ASSESS RETURNS WITH IV(IRC) = 8 OR 9.  SEE ALSO V(LMAX0).
C V(STPPAR) (IN)  MARQUARDT PARAMETER -- 0 MEANS FULL NEWTON STEP.
C V(TUNER1) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION
C             REDUCTION WAS MUCH LESS THAN EXPECTED.  SUGGESTED
C             VALUE = 0.1.
C V(TUNER2) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION
C             REDUCTION WAS LARGE ENOUGH TO ACCEPT STEP.  SUGGESTED
C             VALUE = 10**-4.
C V(TUNER3) (IN)  TUNING CONSTANT USED TO DECIDE IF THE RADIUS
C             SHOULD BE INCREASED.  SUGGESTED VALUE = 0.75.
C  V(XCTOL) (IN)  X-CONVERGENCE CRITERION.  IF STEP IS A NEWTON STEP
C             (V(STPPAR) = 0) HAVING V(RELDX) .LE. V(XCTOL) AND GIVING
C             AT MOST TWICE THE PREDICTED FUNCTION DECREASE, THEN
C             ASSESS RETURNS IV(IRC) = 7 OR 9.
C  V(XFTOL) (IN)  FALSE CONVERGENCE TOLERANCE.  IF STEP GAVE NO OR ONLY
C             A SMALL FUNCTION DECREASE AND V(RELDX) .LE. V(XFTOL),
C             THEN ASSESS RETURNS WITH IV(IRC) = 12.
C
C-------------------------------  NOTES  -------------------------------
C
C  ***  APPLICATION AND USAGE RESTRICTIONS  ***
C
C        THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR
C     LEAST-SQUARES) PACKAGE.  IT MAY BE USED IN ANY UNCONSTRAINED
C     MINIMIZATION SOLVER THAT USES DOGLEG, GOLDFELD-QUANDT-TROTTER,
C     OR LEVENBERG-MARQUARDT STEPS.
C
C  ***  ALGORITHM NOTES  ***
C
C        SEE (1) FOR FURTHER DISCUSSION OF THE ASSESSING AND MODEL
C     SWITCHING STRATEGIES.  WHILE NL2SOL CONSIDERS ONLY TWO MODELS,
C     ASSESS IS DESIGNED TO HANDLE ANY NUMBER OF MODELS.
C
C  ***  USAGE NOTES  ***
C
C        ON THE FIRST CALL OF AN ITERATION, ONLY THE I/O VARIABLES
C     STEP, X, IV(IRC), IV(MODEL), V(F), V(DSTNRM), V(GTSTEP), AND
C     V(PREDUC) NEED HAVE BEEN INITIALIZED.  BETWEEN CALLS, NO I/O
C     VALUES EXECPT STEP, X, IV(MODEL), V(F) AND THE STOPPING TOLER-
C     ANCES SHOULD BE CHANGED.
C        AFTER A RETURN FOR CONVERGENCE OR FALSE CONVERGENCE, ONE CAN
C     CHANGE THE STOPPING TOLERANCES AND CALL ASSESS AGAIN, IN WHICH
C     CASE THE STOPPING TESTS WILL BE REPEATED.
C
C  ***  REFERENCES  ***
C
C     (1) DENNIS, J.E., JR., GAY, D.M., AND WELSCH, R.E. (1980),
C        AN ADAPTIVE NONLINEAR LEAST-SQUARES ALGORITHM,
C        SUBMITTED TO ACM TRANS. MATH. SOFTWARE.
C
C     (2) POWELL, M.J.D. (1970)  A FORTRAN SUBROUTINE FOR SOLVING
C        SYSTEMS OF NONLINEAR ALGEBRAIC EQUATIONS, IN NUMERICAL
C        METHODS FOR NONLINEAR ALGEBRAIC EQUATIONS, EDITED BY
C        P. RABINOWITZ, GORDON AND BREACH, LONDON.
C
C  ***  HISTORY  ***
C
C        JOHN DENNIS DESIGNED MUCH OF THIS ROUTINE, STARTING WITH
C     IDEAS IN (2). ROY WELSCH SUGGESTED THE MODEL SWITCHING STRATEGY.
C        DAVID GAY AND STEPHEN PETERS CAST THIS SUBROUTINE INTO A MORE
C     PORTABLE FORM (WINTER 1977), AND DAVID GAY CAST IT INTO ITS
C     PRESENT FORM (FALL 1978).
C
C  ***  GENERAL  ***
C
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
C     MCS-7906671.
C
C------------------------  EXTERNAL QUANTITIES  ------------------------
C
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL RELDST, VCOPY
C     DOUBLE PRECISION D1MACH, RELDST
C
C VCOPY.... COPIES ONE VECTOR TO ANOTHER.
C
C/
C  ***  NO COMMON BLOCKS  ***
C
C--------------------------  LOCAL VARIABLES  --------------------------
C
C     LOGICAL GOODX
C     INTEGER I, NFC
C     DOUBLE PRECISION EMAX, GTS, HALF, ONE, RELDX1, RFAC1,
C    +                 TEMP, TWO, XMAX, ZERO
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER AFCTOL, DECFAC, DSTNRM, DSTSAV, DST0, F, FDIF, FLSTGD, F0,
C    1        GTSLST, GTSTEP, INCFAC, IRC, LMAX0, MLSTGD, MODEL, NFCALL,
C    2        NFGCAL, NREDUC, PLSTGD, PREDUC, RADFAC, RADINC, RDFCMN,
C    3        RDFCMX, RELDX, RESTOR, RFCTOL, STAGE, STGLIM, STPPAR,
C    4        SWITCH, TOOBIG, TUNER1, TUNER2, TUNER3, XCTOL, XFTOL,
C    5        XIRC
C
C  ***  DATA INITIALIZATIONS  ***
C
      DATA HALF/0.5D+00/, ONE/1.0D+00/, TWO/2.0D+00/, ZERO/0.0D+00/
C
      DATA IRC/3/, MLSTGD/4/, MODEL/5/, NFCALL/6/,
     +     NFGCAL/7/, RADINC/8/, RESTOR/9/, STAGE/10/,
     +     STGLIM/11/, SWITCH/12/, TOOBIG/2/, XIRC/13/
      DATA AFCTOL/31/, DECFAC/22/, DSTNRM/2/, DST0/3/,
     +     DSTSAV/18/, F/10/, FDIF/11/, FLSTGD/12/, F0/13/,
     +     GTSLST/14/, GTSTEP/4/, INCFAC/23/,
     +     LMAX0/35/, NREDUC/6/, PLSTGD/15/, PREDUC/7/,
     +     RADFAC/16/, RDFCMN/24/, RDFCMX/25/,
     +     RELDX/17/, RFCTOL/32/, STPPAR/5/, TUNER1/26/,
     +     TUNER2/27/, TUNER3/28/, XCTOL/33/, XFTOL/34/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      NFC = IV(NFCALL)
      IV(SWITCH) = 0
      IV(RESTOR) = 0
      RFAC1 = ONE
      GOODX = .TRUE.
      I = IV(IRC)
      IF (I .GE. 1 .AND. I .LE. 12)
     +             GO TO (20,30,10,10,40,360,290,290,290,290,290,140), I
         IV(IRC) = 13
         GO TO 999
C
C  ***  INITIALIZE FOR NEW ITERATION  ***
C
 10   IV(STAGE) = 1
      IV(RADINC) = 0
      V(FLSTGD) = V(F0)
      IF (IV(TOOBIG) .EQ. 0) GO TO 90
         IV(STAGE) = -1
         IV(XIRC) = I
         GO TO 60
C
C  ***  STEP WAS RECOMPUTED WITH NEW MODEL OR SMALLER RADIUS  ***
C  ***  FIRST DECIDE WHICH  ***
C
 20   IF (IV(MODEL) .NE. IV(MLSTGD)) GO TO 30
C        ***  OLD MODEL RETAINED, SMALLER RADIUS TRIED  ***
C        ***  DO NOT CONSIDER ANY MORE NEW MODELS THIS ITERATION  ***
         IV(STAGE) = IV(STGLIM)
         IV(RADINC) = -1
         GO TO 90
C
C  ***  A NEW MODEL IS BEING TRIED.  DECIDE WHETHER TO KEEP IT.  ***
C
 30   IV(STAGE) = IV(STAGE) + 1
C
C     ***  NOW WE ADD THE POSSIBILTIY THAT STEP WAS RECOMPUTED WITH  ***
C     ***  THE SAME MODEL, PERHAPS BECAUSE OF AN OVERSIZED STEP.     ***
C
 40   IF (IV(STAGE) .GT. 0) GO TO 50
C
C        ***  STEP WAS RECOMPUTED BECAUSE IT WAS TOO BIG.  ***
C
         IF (IV(TOOBIG) .NE. 0) GO TO 60
C
C        ***  RESTORE IV(STAGE) AND PICK UP WHERE WE LEFT OFF.  ***
C
         IV(STAGE) = -IV(STAGE)
         I = IV(XIRC)
         GO TO (20, 30, 90, 90, 70), I
C
 50   IF (IV(TOOBIG) .EQ. 0) GO TO 70
C
C  ***  HANDLE OVERSIZE STEP  ***
C
      IF (IV(RADINC) .GT. 0) GO TO 80
         IV(STAGE) = -IV(STAGE)
         IV(XIRC) = IV(IRC)
C
 60      V(RADFAC) = V(DECFAC)
         IV(RADINC) = IV(RADINC) - 1
         IV(IRC) = 5
         GO TO 999
C
 70   IF (V(F) .LT. V(FLSTGD)) GO TO 90
C
C     *** THE NEW STEP IS A LOSER.  RESTORE OLD MODEL.  ***
C
      IF (IV(MODEL) .EQ. IV(MLSTGD)) GO TO 80
         IV(MODEL) = IV(MLSTGD)
         IV(SWITCH) = 1
C
C     ***  RESTORE STEP, ETC. ONLY IF A PREVIOUS STEP DECREASED V(F).
C
 80   IF (V(FLSTGD) .GE. V(F0)) GO TO 90
         IV(RESTOR) = 1
         V(F) = V(FLSTGD)
         V(PREDUC) = V(PLSTGD)
         V(GTSTEP) = V(GTSLST)
         IF (IV(SWITCH) .EQ. 0) RFAC1 = V(DSTNRM) / V(DSTSAV)
         V(DSTNRM) = V(DSTSAV)
         NFC = IV(NFGCAL)
         GOODX = .FALSE.
C
C
C  ***  COMPUTE RELATIVE CHANGE IN X BY CURRENT STEP  ***
C
 90   RELDX1 = RELDST(P, D, X, X0)
C
C  ***  RESTORE X AND STEP IF NECESSARY  ***
C
      IF (GOODX) GO TO 105
      DO 100 I = 1, P
         STEP(I) = STLSTG(I)
         X(I) = X0(I) + STLSTG(I)
 100     CONTINUE
C
 105  V(FDIF) = V(F0) - V(F)
      TEMP = 0.0
      IF (V(PREDUC).GT.D1MACH(1)/V(TUNER2)) TEMP = V(TUNER2) * V(PREDUC)
      IF (V(FDIF).GT.TEMP) GO TO 120
C
C        ***  NO (OR ONLY A TRIVIAL) FUNCTION DECREASE
C        ***  -- SO TRY NEW MODEL OR SMALLER RADIUS
C
         V(RELDX) = RELDX1
         IF (V(F) .LT. V(F0)) GO TO 110
              IV(MLSTGD) = IV(MODEL)
              V(FLSTGD) = V(F)
              V(F) = V(F0)
              CALL VCOPY(P, X, X0)
              IV(RESTOR) = 1
              GO TO 115
 110     IV(NFGCAL) = NFC
 115     IV(IRC) = 1
         IF (IV(STAGE) .LT. IV(STGLIM)) GO TO 130
              IV(IRC) = 5
              IV(RADINC) = IV(RADINC) - 1
              GO TO 130
C
C  ***  NONTRIVIAL FUNCTION DECREASE ACHIEVED  ***
C
 120  IV(NFGCAL) = NFC
      RFAC1 = ONE
      IF (GOODX) V(RELDX) = RELDX1
      V(DSTSAV) = V(DSTNRM)
      IF (V(FDIF) .GT. V(PREDUC)*V(TUNER1)) GO TO 200
C
C  ***  DECREASE WAS MUCH LESS THAN PREDICTED -- EITHER CHANGE MODELS
C  ***  OR ACCEPT STEP WITH DECREASED RADIUS.
C
      IF (IV(STAGE) .GE. IV(STGLIM)) GO TO 125
C        ***  CONSIDER SWITCHING MODELS  ***
         IV(IRC) = 2
         GO TO 130
C
C     ***  ACCEPT STEP WITH DECREASED RADIUS  ***
C
 125  IV(IRC) = 4
C
C  ***  SET V(RADFAC) TO FLETCHER*S DECREASE FACTOR  ***
C
 130  IV(XIRC) = IV(IRC)
      EMAX = V(GTSTEP) + V(FDIF)
      V(RADFAC) = HALF * RFAC1
      IF (EMAX .LT. V(GTSTEP)) V(RADFAC) = RFAC1 * MAX(V(RDFCMN),
     +                                           HALF * V(GTSTEP)/EMAX)
C
C  ***  DO FALSE CONVERGENCE TEST  ***
C
 140  IF (V(RELDX) .LE. V(XFTOL)) GO TO 160
         IV(IRC) = IV(XIRC)
         IF (V(F) .LT. V(F0)) GO TO 230
              GO TO 300
C
 160  IV(IRC) = 12
      GO TO 310
C
C  ***  HANDLE GOOD FUNCTION DECREASE  ***
C
 200  IF (V(FDIF) .LT. (-V(TUNER3) * V(GTSTEP))) GO TO 260
C
C     ***  INCREASING RADIUS LOOKS WORTHWHILE.  SEE IF WE JUST
C     ***  RECOMPUTED STEP WITH A DECREASED RADIUS OR RESTORED STEP
C     ***  AFTER RECOMPUTING IT WITH A LARGER RADIUS.
C
      IF (IV(RADINC) .LT. 0) GO TO 260
      IF (IV(RESTOR) .EQ. 1) GO TO 260
C
C        ***  WE DID NOT.  TRY A LONGER STEP UNLESS THIS WAS A NEWTON
C        ***  STEP.
C
         V(RADFAC) = V(RDFCMX)
         GTS = V(GTSTEP)
         IF (V(FDIF) .LT. (HALF/V(RADFAC) - ONE) * GTS)
     +            V(RADFAC) = MAX(V(INCFAC), HALF*GTS/(GTS + V(FDIF)))
         IV(IRC) = 4
         IF (V(STPPAR) .EQ. ZERO) GO TO 300
C             ***  STEP WAS NOT A NEWTON STEP.  RECOMPUTE IT WITH
C             ***  A LARGER RADIUS.
              IV(IRC) = 5
              IV(RADINC) = IV(RADINC) + 1
C
C  ***  SAVE VALUES CORRESPONDING TO GOOD STEP  ***
C
 230  V(FLSTGD) = V(F)
      IV(MLSTGD) = IV(MODEL)
      CALL VCOPY(P, STLSTG, STEP)
      V(DSTSAV) = V(DSTNRM)
      IV(NFGCAL) = NFC
      V(PLSTGD) = V(PREDUC)
      V(GTSLST) = V(GTSTEP)
      GO TO 300
C
C  ***  ACCEPT STEP WITH RADIUS UNCHANGED  ***
C
 260  V(RADFAC) = ONE
      IV(IRC) = 3
      GO TO 300
C
C  ***  COME HERE FOR A RESTART AFTER CONVERGENCE  ***
C
 290  IV(IRC) = IV(XIRC)
      IF (V(DSTSAV) .GE. ZERO) GO TO 310
         IV(IRC) = 12
         GO TO 310
C
C  ***  PERFORM CONVERGENCE TESTS  ***
C
 300  IV(XIRC) = IV(IRC)
 310  IF (ABS(V(F)) .LT. V(AFCTOL)) IV(IRC) = 10
      IF (HALF * V(FDIF) .GT. V(PREDUC)) GO TO 999
      EMAX = 0.0
      IF (ABS(V(F0)).GT.D1MACH(1)/V(RFCTOL))
     +   EMAX = V(RFCTOL) * ABS(V(F0))
      IF (V(DSTNRM) .GT. V(LMAX0) .AND. V(PREDUC) .LE. EMAX)
     +                       IV(IRC) = 11
      IF (V(DST0) .LT. ZERO) GO TO 320
      I = 0
      IF ((V(NREDUC) .GT. ZERO .AND. V(NREDUC) .LE. EMAX) .OR.
     +    (V(NREDUC) .EQ. ZERO. AND. V(PREDUC) .EQ. ZERO))  I = 2
      IF (V(STPPAR) .EQ. ZERO .AND. V(RELDX) .LE. V(XCTOL)) I = I + 1
      IF (I .GT. 0) IV(IRC) = I + 6
C
C  ***  CONSIDER RECOMPUTING STEP OF LENGTH V(LMAX0) FOR SINGULAR
C  ***  CONVERGENCE TEST.
C
 320  IF (ABS(IV(IRC)-3) .GT. 1 .AND. IV(IRC) .NE. 12) GO TO 999
      IF (V(DSTNRM) .GT. V(LMAX0)) GO TO 330
         IF (V(PREDUC) .GE. EMAX) GO TO 999
              IF (V(DST0) .LT. ZERO) GO TO 340
                   IF (HALF * V(DST0) .LE. V(LMAX0)) GO TO 999
                        GO TO 340
 330  IF (HALF * V(DSTNRM) .LE. V(LMAX0)) GO TO 999
      XMAX = V(LMAX0) / V(DSTNRM)
      IF (XMAX * (TWO - XMAX) * V(PREDUC) .GE. EMAX) GO TO 999
 340  IF (V(NREDUC) .LT. ZERO) GO TO 370
C
C  ***  RECOMPUTE V(PREDUC) FOR USE IN SINGULAR CONVERGENCE TEST  ***
C
      V(GTSLST) = V(GTSTEP)
      V(DSTSAV) = V(DSTNRM)
      IF (IV(IRC) .EQ. 12) V(DSTSAV) = -V(DSTSAV)
      V(PLSTGD) = V(PREDUC)
      IV(IRC) = 6
      CALL VCOPY(P, STLSTG, STEP)
      GO TO 999
C
C  ***  PERFORM SINGULAR CONVERGENCE TEST WITH RECOMPUTED V(PREDUC)  ***
C
 360  V(GTSTEP) = V(GTSLST)
      V(DSTNRM) = ABS(V(DSTSAV))
      CALL VCOPY(P, STEP, STLSTG)
      IV(IRC) = IV(XIRC)
      IF (V(DSTSAV) .LE. ZERO) IV(IRC) = 12
      V(NREDUC) = -V(PREDUC)
      V(PREDUC) = V(PLSTGD)
 370  IF (-V(NREDUC) .LE. V(RFCTOL) * ABS(V(F0))) IV(IRC) = 11
C
 999  RETURN
      END
      SUBROUTINE COVCLC(COVIRC, D, IV, J, N, NN, P, R, V, X)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  COMPUTE COVARIANCE MATRIX FOR NL2ITR (NL2SOL VERSION 2.2)  ***
C
C  ***  LET K = ABS(IV(COVREQ).  FOR K .LE. 2, A FINITE DIFFERENCE
C  ***  HESSIAN H IS COMPUTED (USING FUNC. AND GRAD. VALUES IF
C  ***  IV(COVREQ) IS NONNEGATIVE, AND USING ONLY FUNC. VALUES IF
C  ***  IV(COVREQ) IS NEGATIVE).  FOR SCALE = 2*F(X) / MAX(1, N-P),
C  ***  WHERE 2*F(X) IS THE RESIDUAL SUM OF SQUARES, COVCLC COMPUTES...
C  ***             K = 0 OR 1...  SCALE * H**-1 * (J**T * J) * H**-1.
C  ***             K = 2...  SCALE * H**-1.
C  ***             K .GE. 3...  SCALE * (J**T * J)**-1.
C
      implicit none
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   COVIRC,N,NN,P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   D(P),J(NN,P),R(N),V(1),X(P)
      INTEGER
     +   IV(1)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   DEL,HALF,NEGPT5,ONE,T,TWO,WK,ZERO
      INTEGER
     +   COV,COVMAT,COVREQ,DELTA,DELTA0,DLTFDC,F,FX,G,G1,GP,GSAVE1,
     +   H,HC,HMI,HPI,HPM,I,IERR,IP1,IPIV0,IPIVI,IPIVK,IPIVOT,IRC,
     +   K,KAGQT,KALM,KIND,KL,L,LMAT,M,MM1,MM1O2,MODE,NFGCAL,PP1O2,
     +   QTR,QTR1,RD,RD1,RSAVE,SAVEI,STP0,STPI,STPM,SWITCH,TOOBIG,
     +   W,W0,W1,WL,XMSAVE
      LOGICAL
     +   HAVEJ
C
C  EXTERNAL SUBROUTINES
      EXTERNAL LINVRT,LITVMU,LIVMUL,LSQRT,LTSQAR,QRFACT,VCOPY,VSCOPY
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER COVIRC, IV(1), N, NN, P
C     DOUBLE PRECISION D(P), J(NN,P), R(N), V(1), X(P)
C     DIMENSION IV(*), V(*)
C
C  ***  LOCAL VARIABLES  ***
C
C     LOGICAL HAVEJ
C     INTEGER COV, GP, GSAVE1, G1, HC, HMI, HPI, HPM, I, IPIVI, IPIVK,
C    1        IP1, IRC, K, KIND, KL, L, M, MM1, MM1O2, PP1O2, QTR1,
C    2        RD1, STPI, STPM, STP0, WL, W0, W1
C     DOUBLE PRECISION DEL, HALF, NEGPT5, ONE, T, TWO, WK, ZERO
C
C/
C  ***  EXTERNAL SUBROUTINES  ***
C
C     EXTERNAL LINVRT, LITVMU, LIVMUL, LSQRT, LTSQAR, QRFACT,
C    1         VCOPY, VSCOPY
C
C LINVRT... INVERT LOWER TRIANGULAR MATRIX.
C LITVMU... APPLY INVERSE TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX.
C LIVMUL... APPLY INVERSE OF COMPACT LOWER TRIANG. MATRIX.
C LSQRT.... COMPUTE CHOLESKY FACTOR OF (LOWER TRINAG. OF) A SYM. MATRIX.
C LTSQAR... GIVEN LOWER TRIANG. MATRIX L, COMPUTE (L**T)*L.
C QRFACT... COMPUTE QR DECOMPOSITION OF A MATRIX.
C VCOPY.... COPY ONE VECTOR TO ANOTHER.
C VSCOPY... SET ALL ELEMENTS OF A VECTOR TO A SCALAR.
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER COVMAT, COVREQ, DELTA, DELTA0, DLTFDC, F, FX, G, H, IERR,
C    1        IPIVOT, IPIV0, KAGQT, KALM, LMAT, MODE, NFGCAL, QTR,
C    2        RD, RSAVE, SAVEI, SWITCH, TOOBIG, W, XMSAVE
C
      DATA HALF/0.5D+00/, NEGPT5/-0.5D+00/, ONE/1.0D+00/, TWO/2.0D+00/,
     +     ZERO/0.0D+00/
C
      DATA COVMAT/26/, COVREQ/15/, DELTA/50/, DELTA0/44/,
     +     DLTFDC/40/, F/10/, FX/46/, G/28/, H/44/, IERR/32/,
     +     IPIVOT/61/, IPIV0/60/, KAGQT/35/, KALM/36/,
     +     LMAT/58/, MODE/38/, NFGCAL/7/, QTR/49/,
     +     RD/51/, RSAVE/52/, SAVEI/54/, SWITCH/12/,
     +     TOOBIG/2/, W/59/, XMSAVE/49/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      COV = IV(LMAT)
C
      COVIRC = 4
      KIND = IV(COVREQ)
      M = IV(MODE)
      IF (M .GT. 0) GO TO 10
         IV(KAGQT) = -1
         IF (IV(KALM) .GT. 0) IV(KALM) = 0
         IF (ABS(KIND) .GE. 3) GO TO 300
         V(FX) = V(F)
         K = IV(RSAVE)
         CALL VCOPY(N, V(K), R)
 10   IF (M .GT. P) GO TO 200
      IF (KIND .LT. 0) GO TO 100
C
C  ***  COMPUTE FINITE DIFFERENCE HESSIAN USING BOTH FUNCTION AND
C  ***  GRADIENT VALUES.
C
      GSAVE1 = IV(W) + P
      G1 = IV(G)
      IF (M .GT. 0) GO TO 15
C        ***  FIRST CALL ON COVCLC.  SET GSAVE = G, TAKE FIRST STEP  ***
         CALL VCOPY(P, V(GSAVE1), V(G1))
         IV(SWITCH) = IV(NFGCAL)
         GO TO 80
C
 15   DEL = V(DELTA)
      X(M) = V(XMSAVE)
      IF (IV(TOOBIG) .EQ. 0) GO TO 30
C
C     ***  HANDLE OVERSIZE V(DELTA)  ***
C
         IF (DEL*X(M) .GT. ZERO) GO TO 20
C             ***  WE ALREADY TRIED SHRINKING V(DELTA), SO QUIT  ***
              IV(COVMAT) = -2
              GO TO 190
C
C        ***  TRY SHRINKING V(DELTA)  ***
 20      DEL = NEGPT5 * DEL
         GO TO 90
C
 30   COV = IV(LMAT)
      GP = G1 + P - 1
C
C  ***  SET  G = (G - GSAVE)/DEL  ***
C
      DO 40 I = G1, GP
         V(I) = (V(I) - V(GSAVE1)) / DEL
         GSAVE1 = GSAVE1 + 1
 40      CONTINUE
C
C  ***  ADD G AS NEW COL. TO FINITE DIFF. HESSIAN MATRIX  ***
C
      K = COV + M*(M-1)/2
      L = K + M - 2
      IF ( M .EQ. 1) GO TO 60
C
C  ***  SET  H(I,M) = 0.5 * (H(I,M) + G(I))  FOR I = 1 TO M-1  ***
C
      DO 50 I = K, L
         V(I) = HALF * (V(I) + V(G1))
         G1 = G1 + 1
 50      CONTINUE
C
C  ***  ADD  H(I,M) = G(I)  FOR I = M TO P  ***
C
 60   L = L + 1
      DO 70 I = M, P
         V(L) = V(G1)
         L = L + I
         G1 = G1 + 1
 70      CONTINUE
C
 80   M = M + 1
      IV(MODE) = M
      IF (M .GT. P) GO TO 190
C
C  ***  CHOOSE NEXT FINITE DIFFERENCE STEP, RETURN TO GET G THERE  ***
C
      DEL = V(DELTA0) * MAX(ONE/D(M), ABS(X(M)))
      IF (X(M) .LT. ZERO) DEL = -DEL
      V(XMSAVE) = X(M)
 90   X(M) = X(M) + DEL
      V(DELTA) = DEL
      COVIRC = 2
      GO TO 999
C
C  ***  COMPUTE FINITE DIFFERENCE HESSIAN USING FUNCTION VALUES ONLY.
C
 100  STP0 = IV(W) + P - 1
      MM1 = M - 1
      MM1O2 = M*MM1/2
      IF (M .GT. 0) GO TO 105
C        ***  FIRST CALL ON COVCLC.  ***
         IV(SAVEI) = 0
         GO TO 180
C
 105  I = IV(SAVEI)
      IF (I .GT. 0) GO TO 160
      IF (IV(TOOBIG) .EQ. 0) GO TO 120
C
C     ***  HANDLE OVERSIZE STEP  ***
C
         STPM = STP0 + M
         DEL = V(STPM)
         IF (DEL*V(XMSAVE) .GT. ZERO) GO TO 110
C             ***  WE ALREADY TRIED SHRINKING THE STEP, SO QUIT  ***
              IV(COVMAT) = -2
              GO TO 999
C
C        ***  TRY SHRINKING THE STEP  ***
 110     DEL = NEGPT5 * DEL
         X(M) = V(XMSAVE) + DEL
         V(STPM) = DEL
         COVIRC = 1
         GO TO 999
C
C  ***  SAVE F(X + STP(M)*E(M)) IN H(P,M)  ***
C
 120  PP1O2 = P * (P-1) / 2
      COV = IV(LMAT)
      HPM = COV + PP1O2 + MM1
      V(HPM) = V(F)
C
C  ***  START COMPUTING ROW M OF THE FINITE DIFFERENCE HESSIAN H.  ***
C
      HMI = COV + MM1O2
      IF (MM1 .EQ. 0) GO TO 140
      HPI = COV + PP1O2
      DO 130 I = 1, MM1
         V(HMI) = V(FX) - (V(F) + V(HPI))
         HMI = HMI + 1
         HPI = HPI + 1
 130     CONTINUE
 140  V(HMI) = V(F) - TWO*V(FX)
C
C  ***  COMPUTE FUNCTION VALUES NEEDED TO COMPLETE ROW M OF H.  ***
C
      I = 1
C
 150  IV(SAVEI) = I
      STPI = STP0 + I
      V(DELTA) = X(I)
      X(I) = X(I) + V(STPI)
      IF (I .EQ. M) X(I) = V(XMSAVE) - V(STPI)
      COVIRC = 1
      GO TO 999
C
 160  X(I) = V(DELTA)
      IF (IV(TOOBIG) .EQ. 0) GO TO 170
C        ***  PUNT IN THE EVENT OF AN OVERSIZE STEP  ***
         IV(COVMAT) = -2
         GO TO 999
C
C  ***  FINISH COMPUTING H(M,I)  ***
C
 170  STPI = STP0 + I
      HMI = COV + MM1O2 + I - 1
      STPM = STP0 + M
      V(HMI) = (V(HMI) + V(F)) / (V(STPI)*V(STPM))
      I = I + 1
      IF (I .LE. M) GO TO 150
      IV(SAVEI) = 0
      X(M) = V(XMSAVE)
C
 180  M = M + 1
      IV(MODE) = M
      IF (M .GT. P) GO TO 190
C
C  ***  PREPARE TO COMPUTE ROW M OF THE FINITE DIFFERENCE HESSIAN H.
C  ***  COMPUTE M-TH STEP SIZE STP(M), THEN RETURN TO OBTAIN
C  ***  F(X + STP(M)*E(M)), WHERE E(M) = M-TH STD. UNIT VECTOR.
C
      DEL = V(DLTFDC) * MAX(ONE/D(M), ABS(X(M)))
      IF (X(M) .LT. ZERO) DEL = -DEL
      V(XMSAVE) = X(M)
      X(M) = X(M) + DEL
      STPM = STP0 + M
      V(STPM) = DEL
      COVIRC = 1
      GO TO 999
C
C  ***  RESTORE R, V(F), ETC.  ***
C
 190  K = IV(RSAVE)
      CALL VCOPY(N, R, V(K))
      V(F) = V(FX)
      IF (KIND .LT. 0) GO TO 200
         IV(NFGCAL) = IV(SWITCH)
         QTR1 = IV(QTR)
         CALL VCOPY(N, V(QTR1), R)
         IF (IV(COVMAT) .LT. 0) GO TO 999
         COVIRC = 3
         GO TO 999
C
 200  COV = IV(LMAT)
C
C  ***  THE COMPLETE FINITE DIFF. HESSIAN IS NOW STORED AT V(COV).   ***
C  ***  USE IT TO COMPUTE THE REQUESTED COVARIANCE MATRIX.           ***
C
C     ***  COMPUTE CHOLESKY FACTOR C OF H = C*(C**T)  ***
C     ***  AND STORE IT AT V(HC).  ***
C
      HC = COV
      IF (ABS(KIND) .EQ. 2) GO TO 210
         HC = ABS(IV(H))
         IV(H) = -HC
 210  CALL LSQRT(1, P, V(HC), V(COV), IRC)
      IV(COVMAT) = -1
      IF (IRC .NE. 0) GO TO 999
C
      W1 = IV(W) + P
      IF (ABS(KIND) .GT. 1) GO TO 350
C
C  ***  COVARIANCE = SCALE * H**-1 * (J**T * J) * H**-1  ***
C
      CALL VSCOPY(P*(P+1)/2, V(COV), ZERO)
      HAVEJ = IV(KALM) .EQ. (-1)
C     ***  HAVEJ = .TRUE. MEANS J IS IN ITS ORIGINAL FORM, WHILE
C     ***  HAVEJ = .FALSE. MEANS QRFACT HAS BEEN APPLIED TO J.
C
      M = P
      IF (HAVEJ) M = N
      W0 = W1 - 1
      RD1 = IV(RD)
      DO 290 I = 1, M
         IF (HAVEJ) GO TO 240
C
C        ***  SET W = IPIVOT * (ROW I OF R MATRIX FROM QRFACT).  ***
C
              CALL VSCOPY(P, V(W1), ZERO)
              IPIVI = IPIV0 + I
              L = W0 + IV(IPIVI)
              V(L) = V(RD1)
              RD1 = RD1 + 1
              IF (I .EQ. P) GO TO 260
              IP1 = I + 1
              DO 230 K = IP1, P
                   IPIVK = IPIV0 + K
                   L = W0 + IV(IPIVK)
                   V(L) = J(I,K)
 230               CONTINUE
              GO TO 260
C
C        ***  SET W = (ROW I OF J).  ***
C
 240     L = W0
         DO 250 K = 1, P
              L = L + 1
              V(L) = J(I,K)
 250          CONTINUE
C
C        ***  SET W = H**-1 * W.  ***
C
 260     CALL LIVMUL(P, V(W1), V(HC), V(W1))
         CALL LITVMU(P, V(W1), V(HC), V(W1))
C
C        ***  ADD  W * W**T  TO COVARIANCE MATRIX.  ***
C
         KL = COV
         DO 280 K = 1, P
              L = W0 + K
              WK = V(L)
              DO 270 L = 1, K
                   WL = W0 + L
                   V(KL) = V(KL)  +  WK * V(WL)
                   KL = KL + 1
 270               CONTINUE
 280          CONTINUE
 290     CONTINUE
      GO TO 380
C
C  ***  COVARIANCE = SCALE * (J**T * J)**-1.  ***
C
 300  RD1 = IV(RD)
      IF (IV(KALM) .NE. (-1)) GO TO 310
C
C        ***  APPLY QRFACT TO J  ***
C
         QTR1 = IV(QTR)
         CALL VCOPY(N, V(QTR1), R)
         W1 = IV(W) + P
         CALL QRFACT(NN, N, P, J, V(RD1), IV(IPIVOT), IV(IERR), 0,
     +               V(W1))
         IV(KALM) = -2
 310  IV(COVMAT) = -1
      IF (IV(IERR) .NE. 0) GO TO 999
      COV = IV(LMAT)
      HC = ABS(IV(H))
      IV(H) = -HC
C
C     ***  SET HC = (R MATRIX FROM QRFACT).  ***
C
      L = HC
      DO 340 I = 1, P
         IF (I .GT. 1) CALL VCOPY(I-1, V(L), J(1,I))
         L = L + I - 1
         V(L) = V(RD1)
         L = L + 1
         RD1 = RD1 + 1
 340     CONTINUE
C
C  ***  THE CHOLESKY FACTOR C OF THE UNSCALED INVERSE COVARIANCE MATRIX
C  ***  (OR PERMUTATION THEREOF) IS STORED AT V(HC).
C
C  ***  SET C = C**-1.
C
 350  CALL LINVRT(P, V(HC), V(HC))
C
C  ***  SET C = C**T * C.
C
      CALL LTSQAR(P, V(HC), V(HC))
C
      IF (HC .EQ. COV) GO TO 380
C
C     ***  C = PERMUTED, UNSCALED COVARIANCE.
C     ***  SET COV = IPIVOT * C * IPIVOT**T.
C
         DO 370 I = 1, P
              M = IPIV0 + I
              IPIVI = IV(M)
              KL = COV-1 + IPIVI*(IPIVI-1)/2
              DO 360 K = 1, I
                   M = IPIV0 + K
                   IPIVK = IV(M)
                   L = KL + IPIVK
                   IF (IPIVK .GT. IPIVI)
     +                       L = L + (IPIVK-IPIVI)*(IPIVK+IPIVI-3)/2
                   V(L) = V(HC)
                   HC = HC + 1
 360               CONTINUE
 370          CONTINUE
C
 380  IV(COVMAT) = COV
C
C  ***  APPLY SCALE FACTOR = (RESID. SUM OF SQUARES) / MAX(1,N-P).
C
      T = V(F) / (HALF * MAX(1,N-P))
      K = COV - 1 + P*(P+1)/2
      DO 390 I = COV, K
 390     V(I) = T * V(I)
C
 999  RETURN
      END
      SUBROUTINE DFAULT(IV, V)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C
C  VARIABLE DECLARATIONS
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   V(45)
      INTEGER
     +   IV(25)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   MACHEP,MEPCRT,ONE,SQTEPS,THREE
      INTEGER
     +   AFCTOL,COSMIN,COVPRT,COVREQ,D0INIT,DECFAC,DELTA0,DFAC,
     +   DINIT,DLTFDC,DLTFDJ,DTYPE,EPSLON,FUZZ,INCFAC,INITS,JTINIT,
     +   LMAX0,MXFCAL,MXITER,OUTLEV,PARPRT,PHMNFC,PHMXFC,PRUNIT,
     +   RDFCMN,RDFCMX,RFCTOL,RLIMIT,SOLPRT,STATPR,TUNER1,TUNER2,
     +   TUNER3,TUNER4,TUNER5,X0PRT,XCTOL,XFTOL
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   RMDCON
      INTEGER
     +   IMDCON
      EXTERNAL RMDCON,IMDCON
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX
C
C  ***  SUPPLY NL2SOL (VERSION 2.2) DEFAULT VALUES TO IV AND V  ***
C
C     INTEGER IV(25)
C     DOUBLE PRECISION V(45)
C/+
C     DOUBLE PRECISION MAX
C/
C     EXTERNAL IMDCON, RMDCON
C     INTEGER IMDCON
C     DOUBLE PRECISION RMDCON
C
C     DOUBLE PRECISION MACHEP, MEPCRT, ONE, SQTEPS, THREE
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER AFCTOL, COSMIN, COVPRT, COVREQ, DECFAC, DELTA0, DFAC,
C    1        DINIT, DLTFDC, DLTFDJ, DTYPE, D0INIT, EPSLON, FUZZ,
C    2        INCFAC, INITS, JTINIT, LMAX0, MXFCAL, MXITER, OUTLEV,
C    3        PARPRT, PHMNFC, PHMXFC, PRUNIT, RDFCMN, RDFCMX,
C    4        RFCTOL, RLIMIT, SOLPRT, STATPR, TUNER1, TUNER2, TUNER3,
C    5        TUNER4, TUNER5, XCTOL, XFTOL, X0PRT
C
      DATA ONE/1.0D+00/, THREE/3.0D+00/
C
C  ***  IV SUBSCRIPT VALUES  ***
C
      DATA COVPRT/14/, COVREQ/15/, DTYPE/16/, INITS/25/,
     +     MXFCAL/17/, MXITER/18/, OUTLEV/19/,
     +     PARPRT/20/, PRUNIT/21/, SOLPRT/22/,
     +     STATPR/23/, X0PRT/24/
C
C  ***  V SUBSCRIPT VALUES  ***
C
      DATA AFCTOL/31/, COSMIN/43/, DECFAC/22/,
     +     DELTA0/44/, DFAC/41/, DINIT/38/, DLTFDC/40/,
     +     DLTFDJ/36/, D0INIT/37/, EPSLON/19/, FUZZ/45/,
     +     INCFAC/23/, JTINIT/39/, LMAX0/35/, PHMNFC/20/,
     +     PHMXFC/21/, RDFCMN/24/, RDFCMX/25/,
     +     RFCTOL/32/, RLIMIT/42/, TUNER1/26/,
     +     TUNER2/27/, TUNER3/28/, TUNER4/29/,
     +     TUNER5/30/, XCTOL/33/, XFTOL/34/
C
C-----------------------------------------------------------------------
C
      IV(1) = 12
      IV(COVPRT) = 1
      IV(COVREQ) = 1
      IV(DTYPE) = 1
      IV(INITS) = 0
      IV(MXFCAL) = 200
      IV(MXITER) = 150
      IV(OUTLEV) = -1
      IV(PARPRT) = 1
      IV(PRUNIT) = IMDCON(1)
      IV(SOLPRT) = 1
      IV(STATPR) = 1
      IV(X0PRT) = 1
C
      MACHEP = RMDCON(3)
      V(AFCTOL) = 1.0D-20
      IF (MACHEP .GT. 1.0D-10) V(AFCTOL) = MACHEP**2
      V(COSMIN) = MAX(1.0D-06, 1.0D+02 * MACHEP)
      V(DECFAC) = 0.5D+00
      SQTEPS = sqrt ( machep )
      V(DELTA0) = SQTEPS
      V(DFAC) = 0.6D+00
      V(DINIT) = 0.0D+00
      MEPCRT = MACHEP ** (ONE/THREE)
      V(DLTFDC) = MEPCRT
      V(DLTFDJ) = SQTEPS
      V(D0INIT) = 1.0D+00
      V(EPSLON) = 0.1D+00
      V(FUZZ) = 1.5D+00
      V(INCFAC) = 2.0D+00
      V(JTINIT) = 1.0D-06
      V(LMAX0) = 100.0D+00
      V(PHMNFC) = -0.1D+00
      V(PHMXFC) = 0.1D+00
      V(RDFCMN) = 0.1D+00
      V(RDFCMX) = 4.0D+00
      V(RFCTOL) = MAX(1.0D-10, MEPCRT**2)
      V(RLIMIT) = RMDCON(5)
      V(TUNER1) = 0.1D+00
      V(TUNER2) = 1.0D-04
      V(TUNER3) = 0.75D+00
      V(TUNER4) = 0.5D+00
      V(TUNER5) = 0.75D+00
      V(XCTOL) = SQTEPS
      V(XFTOL) = 1.0D+02 * MACHEP
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DOTPRD(P, X, Y)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  RETURN THE INNER PRODUCT OF THE P-VECTORS X AND Y.  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   X(*),Y(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ONE,SQTETA,T,ZERO
      INTEGER
     +   I
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   RMDCON
      EXTERNAL RMDCON
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX
C
C     INTEGER P
C     DOUBLE PRECISION X(*), Y(*)
C
C     INTEGER I
C     DOUBLE PRECISION ONE, SQTETA, T, ZERO
C/+
C     DOUBLE PRECISION MAX, ABS
C/
C     EXTERNAL RMDCON
C     DOUBLE PRECISION RMDCON
C
C  ***  RMDCON(2) RETURNS A MACHINE DEPENDENT CONSTANT, SQTETA, WHICH
C  ***  IS SLIGHTLY LARGER THAN THE SMALLEST POSITIVE NUMBER THAT
C  ***  CAN BE SQUARED WITHOUT UNDERFLOWING.
C
      DATA ONE/1.0D+00/, SQTETA/0.0D+00/, ZERO/0.0D+00/
C
      DOTPRD = ZERO
      IF (P .LE. 0) GO TO 999
      IF (SQTETA .EQ. ZERO) SQTETA = RMDCON(2)
      DO 20 I = 1, P
         T = MAX(ABS(X(I)), ABS(Y(I)))
         IF (T .GT. ONE) GO TO 10
         IF (T .LT. SQTETA) GO TO 20
         T = (X(I)/SQTETA)*Y(I)
         IF (ABS(T) .LT. SQTETA) GO TO 20
 10      DOTPRD = DOTPRD + X(I)*Y(I)
 20   CONTINUE
C
 999  RETURN
      END
      SUBROUTINE DUPDAT(D, IV, J, N, NN, P, V)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  UPDATE SCALE VECTOR D FOR NL2ITR (NL2SOL VERSION 2.2)  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,NN,P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   D(P),J(NN,P),V(1)
      INTEGER
     +   IV(1)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   SII,T,VDFAC,ZERO
      INTEGER
     +   D0,DFAC,DTYPE,I,JTOL0,JTOLI,NITER,S,S1
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   V2NORM
      EXTERNAL V2NORM
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX,SQRT
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER IV(1), N, NN, P
C     DOUBLE PRECISION D(P), J(NN,P), V(1)
C     DIMENSION IV(*), V(*)
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER D0, I, JTOLI, S1
C     DOUBLE PRECISION SII, T, VDFAC
C
C     ***  CONSTANTS  ***
C     DOUBLE PRECISION ZERO
C
C/
C  ***  EXTERNAL FUNCTION  ***
C
C     EXTERNAL V2NORM
C     DOUBLE PRECISION V2NORM
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER DFAC, DTYPE, JTOL0, NITER, S
      DATA DFAC/41/, DTYPE/16/, JTOL0/86/, NITER/31/, S/53/
C
      DATA ZERO/0.0D+00/

      I = IV(DTYPE)
      IF (I .EQ. 1) GO TO 20
         IF (IV(NITER) .GT. 0) GO TO 999
C
 20   VDFAC = V(DFAC)
      D0 = JTOL0 + P
      S1 = IV(S) - 1
      DO 30 I = 1, P
         S1 = S1 + I
         SII = V(S1)
         T = V2NORM(N, J(1,I))
         IF (SII .GT. ZERO) T = SQRT(T*T + SII)
         JTOLI = JTOL0 + I
         D0 = D0 + 1
         IF (T .LT. V(JTOLI)) T = MAX(V(D0), V(JTOLI))
         D(I) = MAX(VDFAC*D(I), T)
 30      CONTINUE
C
 999  RETURN
      END
      SUBROUTINE GQTSTP(D, DIG, DIHDI, KA, L, P, STEP, V, W)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  *** COMPUTE GOLDFELD-QUANDT-TROTTER STEP BY MORE HEBDEN TECHNIQUE ***
C  ***  (NL2SOL VERSION 2.2)  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   KA,P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   D(P),DIG(P),DIHDI(1),L(1),STEP(P),V(21),W(1)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   AKI,AKK,ALPHAK,DELTA,DGXFAC,DST,EPSFAC,EPSO6,FOUR,HALF,KAPPA,
     +   LK,NEGONE,OLDPHI,ONE,P001,PHI,PHIMAX,PHIMIN,PSIFAC,RAD,ROOT,
     +   SI,SIX,SK,SW,T,T1,THREE,TWO,TWOPSI,UK,WI,ZERO
      INTEGER
     +   DGGDMX,DGNORM,DIAG,DIAG0,DST0,DSTNRM,DSTSAV,EMAX,EMIN,
     +   EPSLON,GTSTEP,I,IM1,INC,IRC,J,K,K1,KALIM,LK0,NREDUC,
     +   PHIPIN,PHMNFC,PHMXFC,PREDUC,Q,Q0,RAD0,RADIUS,STPPAR,UK0,X,
     +   X0
      LOGICAL
     +   RESTRT
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD,LSVMIN,RMDCON,V2NORM
      EXTERNAL DOTPRD,LSVMIN,RMDCON,V2NORM
C
C  EXTERNAL SUBROUTINES
      EXTERNAL LITVMU,LIVMUL,LSQRT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX,MIN,SQRT
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER KA, P
C     DOUBLE PRECISION D(P), DIG(P), DIHDI(1), L(1), V(21), STEP(P),
C    1                 W(1)
C     DIMENSION DIHDI(P*(P+1)/2), L(P*(P+1)/2), W(4*P+7)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  PURPOSE  ***
C
C        GIVEN THE (COMPACTLY STORED) LOWER TRIANGLE OF A SCALED
C     HESSIAN (APPROXIMATION) AND A NONZERO SCALED GRADIENT VECTOR,
C     THIS SUBROUTINE COMPUTES A GOLDFELD-QUANDT-TROTTER STEP OF
C     APPROXIMATE LENGTH V(RADIUS) BY THE MORE HEBDEN TECHNIQUE.  IN
C     OTHER WORDS, STEP IS COMPUTED TO (APPROXIMATELY) MINIMIZE
C     PSI(STEP) = (G**T)*STEP + 0.5*(STEP**T)*H*STEP  SUCH THAT THE
C     2-NORM OF D*STEP IS AT MOST (APPROXIMATELY) V(RADIUS), WHERE
C     G  IS THE GRADIENT,  H  IS THE HESSIAN, AND  D  IS A DIAGONAL
C     SCALE MATRIX WHOSE DIAGONAL IS STORED IN THE PARAMETER D.
C     (GQTSTP ASSUMES  DIG = D**-1 * G  AND  DIHDI = D**-1 * H * D**-1.)
C     IF G = 0, HOWEVER, STEP = 0 IS RETURNED (EVEN AT A SADDLE POINT).
C
C  ***  PARAMETER DESCRIPTION  ***
C
C     D (IN)  = THE SCALE VECTOR, I.E. THE DIAGONAL OF THE SCALE
C              MATRIX  D  MENTIONED ABOVE UNDER PURPOSE.
C   DIG (IN)  = THE SCALED GRADIENT VECTOR, D**-1 * G.  IF G = 0, THEN
C              STEP = 0  AND  V(STPPAR) = 0  ARE RETURNED.
C DIHDI (IN)  = LOWER TRIANGLE OF THE SCALED HESSIAN (APPROXIMATION),
C              I.E., D**-1 * H * D**-1, STORED COMPACTLY BY ROWS., I.E.,
C              IN THE ORDER (1,1), (2,1), (2,2), (3,1), (3,2), ETC.
C    KA (I/O) = THE NUMBER OF HEBDEN ITERATIONS (SO FAR) TAKEN TO DETER-
C              MINE STEP.  KA .LT. 0 ON INPUT MEANS THIS IS THE FIRST
C              ATTEMPT TO DETERMINE STEP (FOR THE PRESENT DIG AND DIHDI)
C              -- KA IS INITIALIZED TO 0 IN THIS CASE.  OUTPUT WITH
C              KA = 0  (OR V(STPPAR) = 0)  MEANS  STEP = -(H**-1)*G.
C     L (I/O) = WORKSPACE OF LENGTH P*(P+1)/2 FOR CHOLESKY FACTORS.
C     P (IN)  = NUMBER OF PARAMETERS -- THE HESSIAN IS A  P X P  MATRIX.
C  STEP (I/O) = THE STEP COMPUTED.
C     V (I/O) CONTAINS VARIOUS CONSTANTS AND VARIABLES DESCRIBED BELOW.
C     W (I/O) = WORKSPACE OF LENGTH 4*P + 6.
C
C  ***  ENTRIES IN V  ***
C
C V(DGNORM) (I/O) = 2-NORM OF (D**-1)*G.
C V(DSTNRM) (OUTPUT) = 2-NORM OF D*STEP.
C V(DST0)   (I/O) = 2-NORM OF D*(H**-1)*G (FOR POS. DEF. H ONLY), OR
C             OVERESTIMATE OF SMALLEST EIGENVALUE OF (D**-1)*H*(D**-1).
C V(EPSLON) (IN)  = MAX. REL. ERROR ALLOWED FOR PSI(STEP).  FOR THE
C             STEP RETURNED, PSI(STEP) WILL EXCEED ITS OPTIMAL VALUE
C             BY LESS THAN -V(EPSLON)*PSI(STEP).  SUGGESTED VALUE = 0.1.
C V(GTSTEP) (OUT) = INNER PRODUCT BETWEEN G AND STEP.
C V(NREDUC) (OUT) = PSI(-(H**-1)*G) = PSI(NEWTON STEP)  (FOR POS. DEF.
C             H ONLY -- V(NREDUC) IS SET TO ZERO OTHERWISE).
C V(PHMNFC) (IN)  = TOL. (TOGETHER WITH V(PHMXFC)) FOR ACCEPTING STEP
C             (MORE*S SIGMA).  THE ERROR V(DSTNRM) - V(RADIUS) MUST LIE
C             BETWEEN V(PHMNFC)*V(RADIUS) AND V(PHMXFC)*V(RADIUS).
C V(PHMXFC) (IN)  (SEE V(PHMNFC).)
C             SUGGESTED VALUES -- V(PHMNFC) = -0.25, V(PHMXFC) = 0.5.
C V(PREDUC) (OUT) = PSI(STEP) = PREDICTED OBJ. FUNC. REDUCTION FOR STEP.
C V(RADIUS) (IN)  = RADIUS OF CURRENT (SCALED) TRUST REGION.
C V(RAD0)   (I/O) = VALUE OF V(RADIUS) FROM PREVIOUS CALL.
C V(STPPAR) (I/O) IS NORMALLY THE MARQUARDT PARAMETER, I.E. THE ALPHA
C             DESCRIBED BELOW UNDER ALGORITHM NOTES.  IF H + ALPHA*D**2
C             (SEE ALGORITHM NOTES) IS (NEARLY) SINGULAR, HOWEVER,
C             THEN V(STPPAR) = -ALPHA.
C
C  ***  USAGE NOTES  ***
C
C     IF IT IS DESIRED TO RECOMPUTE STEP USING A DIFFERENT VALUE OF
C     V(RADIUS), THEN THIS ROUTINE MAY BE RESTARTED BY CALLING IT
C     WITH ALL PARAMETERS UNCHANGED EXCEPT V(RADIUS).  (THIS EXPLAINS
C     WHY STEP AND W ARE LISTED AS I/O).  ON AN INTIIAL CALL (ONE WITH
C     KA .LT. 0), STEP AND W NEED NOT BE INITIALIZED AND ONLY COMPO-
C     NENTS V(EPSLON), V(STPPAR), V(PHMNFC), V(PHMXFC), V(RADIUS), AND
C     V(RAD0) OF V MUST BE INITIALIZED.  TO COMPUTE STEP FROM A SADDLE
C     POINT (WHERE THE TRUE GRADIENT VANISHES AND H HAS A NEGATIVE
C     EIGENVALUE), A NONZERO G WITH SMALL COMPONENTS SHOULD BE PASSED.
C
C  ***  APPLICATION AND USAGE RESTRICTIONS  ***
C
C     THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR LEAST-
C     SQUARES) PACKAGE (REF. 1), BUT IT COULD BE USED IN SOLVING ANY
C     UNCONSTRAINED MINIMIZATION PROBLEM.
C
C  ***  ALGORITHM NOTES  ***
C
C        THE DESIRED G-Q-T STEP (REF. 2, 3, 4) SATISFIES
C     (H + ALPHA*D**2)*STEP = -G  FOR SOME NONNEGATIVE ALPHA SUCH THAT
C     H + ALPHA*D**2 IS POSITIVE SEMIDEFINITE.  ALPHA AND STEP ARE
C     COMPUTED BY A SCHEME ANALOGOUS TO THE ONE DESCRIBED IN REF. 5.
C     ESTIMATES OF THE SMALLEST AND LARGEST EIGENVALUES OF THE HESSIAN
C     ARE OBTAINED FROM THE GERSCHGORIN CIRCLE THEOREM ENHANCED BY A
C     SIMPLE FORM OF THE SCALING DESCRIBED IN REF. 6.  CASES IN WHICH
C     H + ALPHA*D**2 IS NEARLY (OR EXACTLY) SINGULAR ARE HANDLED BY
C     THE TECHNIQUE DISCUSSED IN REF. 2.  IN THESE CASES, A STEP OF
C     (EXACT) LENGTH V(RADIUS) IS RETURNED FOR WHICH PSI(STEP) EXCEEDS
C     ITS OPTIMAL VALUE BY LESS THAN -V(EPSLON)*PSI(STEP).
C
C  ***  FUNCTIONS AND SUBROUTINES CALLED  ***
C
C DOTPRD - RETURNS INNER PRODUCT OF TWO VECTORS.
C LITVMU - APPLIES INVERSE TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX.
C LIVMUL - APPLIES INVERSE OF COMPACT LOWER TRIANG. MATRIX.
C LSQRT  - FINDS CHOLESKY FACTOR (OF COMPACTLY STORED LOWER TRIANG.).
C LSVMIN - RETURNS APPROX. TO MIN. SING. VALUE OF LOWER TRIANG. MATRIX.
C RMDCON - RETURNS MACHINE DEPENDENT CONSTANTS.
C V2NORM - RETURNS 2-NORM OF A VECTOR.
C
C  ***  REFERENCES  ***
C
C 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1980), AN ADAPTIVE
C             NONLINEAR LEAST-SQUARES ALGORITHM, (SUBMITTED TO ACM
C             TRANS. MATH. SOFTWARE).
C 2.  GAY, D.M. (1979), COMPUTING OPTIMAL ELLIPTICALLY CONSTRAINED
C             STEPS, MRC TECH. SUMMARY REPORT NO. 2013, MATH. RESEARCH
C             CENTER, UNIV. OF WISCONSIN-MADISON.
C 3.  GOLDFELD, S.M., QUANDT, R.E., AND TROTTER, H.F. (1966),
C             MAXIMIZATION BY QUADRATIC HILL-CLIMBING, ECONOMETRICA 34,
C             PP. 541-551.
C 4.  HEBDEN, M.D. (1973), AN ALGORITHM FOR MINIMIZATION USING EXACT
C             SECOND DERIVATIVES, REPORT T.P. 515, THEORETICAL PHYSICS
C             DIV., A.E.R.E. HARWELL, OXON., ENGLAND.
C 5.  MORE, J.J. (1978), THE LEVENBERG-MARQUARDT ALGORITHM, IMPLEMEN-
C             TATION AND THEORY, PP.105-116 OF SPRINGER LECTURE NOTES
C             IN MATHEMATICS NO. 630, EDITED BY G.A. WATSON, SPRINGER-
C             VERLAG, BERLIN AND NEW YORK.
C 6.  VARGA, R.S. (1965), MINIMAL GERSCHGORIN SETS, PACIFIC J. MATH. 15,
C             PP. 719-729.
C
C  ***  GENERAL  ***
C
C     CODED BY DAVID M. GAY.
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
C     MCS-7906671.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  LOCAL VARIABLES  ***
C
C     LOGICAL RESTRT
C     INTEGER DGGDMX, DIAG, DIAG0, DSTSAV, EMAX, EMIN, I, IM1, INC, IRC,
C    1        J, K, KALIM, K1, LK0, PHIPIN, Q, Q0, UK0, X, X0
C     DOUBLE PRECISION ALPHAK, AKI, AKK, DELTA, DST, EPSO6, LK,
C    1                 OLDPHI, PHI, PHIMAX, PHIMIN, PSIFAC, RAD,
C    2                 ROOT, SI, SK, SW, T, TWOPSI, T1, UK, WI
C
C     ***  CONSTANTS  ***
C     DOUBLE PRECISION DGXFAC, EPSFAC, FOUR, HALF, KAPPA, NEGONE, ONE,
C    1                 P001, SIX, THREE, TWO, ZERO
C
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL DOTPRD, LITVMU, LIVMUL, LSQRT, LSVMIN, RMDCON, V2NORM
C     DOUBLE PRECISION DOTPRD, LSVMIN, RMDCON, V2NORM
C
C  ***  SUBSCRIPTS FOR V  ***
C
C     INTEGER DGNORM, DSTNRM, DST0, EPSLON, GTSTEP, STPPAR, NREDUC,
C    1        PHMNFC, PHMXFC, PREDUC, RADIUS, RAD0
      DATA DGNORM/1/, DSTNRM/2/, DST0/3/, EPSLON/19/,
     +     GTSTEP/4/, NREDUC/6/, PHMNFC/20/,
     +     PHMXFC/21/, PREDUC/7/, RADIUS/8/,
     +     RAD0/9/, STPPAR/5/
C
      DATA DGXFAC/0.0D+00/
      DATA EPSFAC/50.0D+00/
      DATA FOUR/4.0D+00/
      DATA HALF/0.5D+00/
      DATA KAPPA/2.0D+00/
      DATA NEGONE/-1.0D+00/, ONE/1.0D+00/, P001/1.0D-03/,
     +     SIX/6.0D+00/, THREE/3.0D+00/, TWO/2.0D+00/, ZERO/0.0D+00/
C
C  ***  BODY  ***
C
C     ***  STORE LARGEST ABS. ENTRY IN (D**-1)*H*(D**-1) AT W(DGGDMX).
      DGGDMX = P + 1
C     ***  STORE GERSCHGORIN OVER- AND UNDERESTIMATES OF THE LARGEST
C     ***  AND SMALLEST EIGENVALUES OF (D**-1)*H*(D**-1) AT W(EMAX)
C     ***  AND W(EMIN) RESPECTIVELY.
      EMAX = DGGDMX + 1
      EMIN = EMAX + 1
C     ***  FOR USE IN RECOMPUTING STEP, THE FINAL VALUES OF LK, UK, DST,
C     ***  AND THE INVERSE DERIVATIVE OF MORE*S PHI AT 0 (FOR POS. DEF.
C     ***  H) ARE STORED IN W(LK0), W(UK0), W(DSTSAV), AND W(PHIPIN)
C     ***  RESPECTIVELY.
      UK = 0.0D+00
      PHI = 0.0D+00
      DST = 0.0D+00
      ALPHAK = 0.0D+00
      LK0 = EMIN + 1
      PHIPIN = LK0 + 1
      UK0 = PHIPIN + 1
      DSTSAV = UK0 + 1
C     ***  STORE DIAG OF (D**-1)*H*(D**-1) IN W(DIAG),...,W(DIAG0+P).
      DIAG0 = DSTSAV
      DIAG = DIAG0 + 1
C     ***  STORE -D*STEP IN W(Q),...,W(Q0+P).
      Q0 = DIAG0 + P
      Q = Q0 + 1
      RAD = V(RADIUS)
C     ***  PHITOL = MAX. ERROR ALLOWED IN DST = V(DSTNRM) = 2-NORM OF
C     ***  D*STEP.
      PHIMAX = V(PHMXFC) * RAD
      PHIMIN = V(PHMNFC) * RAD
C     ***  EPSO6 AND PSIFAC ARE USED IN CHECKING FOR THE SPECIAL CASE
C     ***  OF (NEARLY) SINGULAR H + ALPHA*D**2 (SEE REF. 2).
      PSIFAC = TWO * V(EPSLON) / (THREE * (FOUR * (V(PHMNFC) + ONE) *
     +                       (KAPPA + ONE)  +  KAPPA  +  TWO) * RAD**2)
C     ***  OLDPHI IS USED TO DETECT LIMITS OF NUMERICAL ACCURACY.  IF
C     ***  WE RECOMPUTE STEP AND IT DOES NOT CHANGE, THEN WE ACCEPT IT.
      OLDPHI = ZERO
      EPSO6 = V(EPSLON)/SIX
      IRC = 0
      RESTRT = .FALSE.
      KALIM = KA + 50
C
C  ***  START OR RESTART, DEPENDING ON KA  ***
C
      IF (KA .GE. 0) GO TO 310
C
C  ***  FRESH START  ***
C
      K = 0
      UK = NEGONE
      KA = 0
      KALIM = 50
C
C     ***  STORE DIAG(DIHDI) IN W(DIAG0+1),...,W(DIAG0+P)  ***
C
      J = 0
      DO 20 I = 1, P
         J = J + I
         K1 = DIAG0 + I
         W(K1) = DIHDI(J)
 20      CONTINUE
C
C     ***  DETERMINE W(DGGDMX), THE LARGEST ELEMENT OF DIHDI  ***
C
      T1 = ZERO
      J = P * (P + 1) / 2
      DO 30 I = 1, J
         T = ABS(DIHDI(I))
         IF (T1 .LT. T) T1 = T
 30      CONTINUE
      W(DGGDMX) = T1
C
C  ***  TRY ALPHA = 0  ***
C
 40   CALL LSQRT(1, P, L, DIHDI, IRC)
      IF (IRC .EQ. 0) GO TO 60
C        ***  INDEF. H -- UNDERESTIMATE SMALLEST EIGENVALUE, USE THIS
C        ***  ESTIMATE TO INITIALIZE LOWER BOUND LK ON ALPHA.
         J = IRC*(IRC+1)/2
         T = L(J)
         L(J) = ONE
         DO 50 I = 1, IRC
 50           W(I) = ZERO
         W(IRC) = ONE
         CALL LITVMU(IRC, W, L, W)
         T1 = V2NORM(IRC, W)
         LK = -T / T1 / T1
         V(DST0) = -LK
         IF (RESTRT) GO TO 210
         V(NREDUC) = ZERO
         GO TO 70
C
C     ***  POSITIVE DEFINITE H -- COMPUTE UNMODIFIED NEWTON STEP.  ***
 60   LK = ZERO
      CALL LIVMUL(P, W(Q), L, DIG)
      V(NREDUC) = HALF * DOTPRD(P, W(Q), W(Q))
      CALL LITVMU(P, W(Q), L, W(Q))
      DST = V2NORM(P, W(Q))
      V(DST0) = DST
      PHI = DST - RAD
      IF (PHI .LE. PHIMAX) GO TO 280
      IF (RESTRT) GO TO 210
C
C  ***  PREPARE TO COMPUTE GERSCHGORIN ESTIMATES OF LARGEST (AND
C  ***  SMALLEST) EIGENVALUES.  ***
C
 70   V(DGNORM) = V2NORM(P, DIG)
      IF (V(DGNORM) .EQ. ZERO) GO TO 450
      K = 0
      DO 100 I = 1, P
         WI = ZERO
         IF (I .EQ. 1) GO TO 90
         IM1 = I - 1
         DO 80 J = 1, IM1
              K = K + 1
              T = ABS(DIHDI(K))
              WI = WI + T
              W(J) = W(J) + T
 80           CONTINUE
 90      W(I) = WI
         K = K + 1
 100     CONTINUE
C
C  ***  (UNDER-)ESTIMATE SMALLEST EIGENVALUE OF (D**-1)*H*(D**-1)  ***
C
      K = 1
      T1 = W(DIAG) - W(1)
      IF (P .LE. 1) GO TO 120
      DO 110 I = 2, P
         J = DIAG0 + I
         T = W(J) - W(I)
         IF (T .GE. T1) GO TO 110
              T1 = T
              K = I
 110     CONTINUE
C
 120  SK = W(K)
      J = DIAG0 + K
      AKK = W(J)
      K1 = K*(K-1)/2 + 1
      INC = 1
      T = ZERO
      DO 150 I = 1, P
         IF (I .EQ. K) GO TO 130
         AKI = ABS(DIHDI(K1))
         SI = W(I)
         J = DIAG0 + I
         T1 = HALF * (AKK - W(J) + SI - AKI)
         T1 = T1 + SQRT(T1*T1 + SK*AKI)
         IF (T .LT. T1) T = T1
         IF (I .LT. K) GO TO 140
 130     INC = I
 140     K1 = K1 + INC
 150     CONTINUE
C
      W(EMIN) = AKK - T
      UK = V(DGNORM)/RAD - W(EMIN)
C
C  ***  COMPUTE GERSCHGORIN (OVER-)ESTIMATE OF LARGEST EIGENVALUE  ***
C
      K = 1
      T1 = W(DIAG) + W(1)
      IF (P .LE. 1) GO TO 170
      DO 160 I = 2, P
         J = DIAG0 + I
         T = W(J) + W(I)
         IF (T .LE. T1) GO TO 160
              T1 = T
              K = I
 160     CONTINUE
C
 170  SK = W(K)
      J = DIAG0 + K
      AKK = W(J)
      K1 = K*(K-1)/2 + 1
      INC = 1
      T = ZERO
      DO 200 I = 1, P
         IF (I .EQ. K) GO TO 180
         AKI = ABS(DIHDI(K1))
         SI = W(I)
         J = DIAG0 + I
         T1 = HALF * (W(J) + SI - AKI - AKK)
         T1 = T1 + SQRT(T1*T1 + SK*AKI)
         IF (T .LT. T1) T = T1
         IF (I .LT. K) GO TO 190
 180     INC = I
 190     K1 = K1 + INC
 200     CONTINUE
C
      W(EMAX) = AKK + T
      LK = MAX(LK, V(DGNORM)/RAD - W(EMAX))
C
C     ***  ALPHAK = CURRENT VALUE OF ALPHA (SEE ALG. NOTES ABOVE).  WE
C     ***  USE MORE*S SCHEME FOR INITIALIZING IT.
      ALPHAK = ABS(V(STPPAR)) * V(RAD0)/RAD
C
      IF (IRC .NE. 0) GO TO 210
C
C  ***  COMPUTE L0 FOR POSITIVE DEFINITE H  ***
C
      CALL LIVMUL(P, W, L, W(Q))
      T = V2NORM(P, W)
      W(PHIPIN) = DST / T / T
      LK = MAX(LK, PHI*W(PHIPIN))
C
C  ***  SAFEGUARD ALPHAK AND ADD ALPHAK*I TO (D**-1)*H*(D**-1)  ***
C
 210  KA = KA + 1
      IF (-V(DST0) .GE. ALPHAK .OR. ALPHAK .LT. LK .OR. ALPHAK .GE. UK)
     +                      ALPHAK = UK * MAX(P001, SQRT(LK/UK))
      K = 0
      DO 220 I = 1, P
         K = K + I
         J = DIAG0 + I
         DIHDI(K) = W(J) + ALPHAK
 220     CONTINUE
C
C  ***  TRY COMPUTING CHOLESKY DECOMPOSITION  ***
C
      CALL LSQRT(1, P, L, DIHDI, IRC)
      IF (IRC .EQ. 0) GO TO 250
C
C  ***  (D**-1)*H*(D**-1) + ALPHAK*I  IS INDEFINITE -- OVERESTIMATE
C  ***  SMALLEST EIGENVALUE FOR USE IN UPDATING LK  ***
C
      J = (IRC*(IRC+1))/2
      T = L(J)
      L(J) = ONE
      DO 230 I = 1, IRC
 230     W(I) = ZERO
      W(IRC) = ONE
      CALL LITVMU(IRC, W, L, W)
      T1 = V2NORM(IRC, W)
      LK = ALPHAK - T/T1/T1
      V(DST0) = -LK
      GO TO 210
C
C  ***  ALPHAK MAKES (D**-1)*H*(D**-1) POSITIVE DEFINITE.
C  ***  COMPUTE Q = -D*STEP, CHECK FOR CONVERGENCE.  ***
C
 250  CALL LIVMUL(P, W(Q), L, DIG)
      CALL LITVMU(P, W(Q), L, W(Q))
      DST = V2NORM(P, W(Q))
      PHI = DST - RAD
      IF (PHI .LE. PHIMAX .AND. PHI .GE. PHIMIN) GO TO 290
      IF (PHI .EQ. OLDPHI) GO TO 290
      OLDPHI = PHI
      IF (PHI .GT. ZERO) GO TO 260
C        ***  CHECK FOR THE SPECIAL CASE OF  H + ALPHA*D**2  (NEARLY)
C        ***  SINGULAR.  DELTA IS .GE. THE SMALLEST EIGENVALUE OF
C        ***  (D**-1)*H*(D**-1) + ALPHAK*I.
         IF (V(DST0) .GT. ZERO) GO TO 260
         DELTA = ALPHAK + V(DST0)
         TWOPSI = ALPHAK*DST*DST + DOTPRD(P, DIG, W(Q))
         IF (DELTA .LT. PSIFAC*TWOPSI) GO TO 270
C
C  ***  UNACCEPTABLE ALPHAK -- UPDATE LK, UK, ALPHAK  ***
C
 260  IF (KA .GE. KALIM) GO TO 290
      CALL LIVMUL(P, W, L, W(Q))
      T1 = V2NORM(P, W)
C     ***  THE FOLLOWING MIN IS NECESSARY BECAUSE OF RESTARTS  ***
      IF (PHI .LT. ZERO) UK = MIN(UK, ALPHAK)
      ALPHAK = ALPHAK  +  (PHI/T1) * (DST/T1) * (DST/RAD)
      LK = MAX(LK, ALPHAK)
      GO TO 210
C
C  ***  DECIDE HOW TO HANDLE (NEARLY) SINGULAR H + ALPHA*D**2  ***
C
C     ***  IF NOT YET AVAILABLE, OBTAIN MACHINE DEPENDENT VALUE DGXFAC.
 270  IF (DGXFAC .EQ. ZERO) DGXFAC = EPSFAC * RMDCON(3)
C
C     ***  NOW DECIDE.  ***
      IF (DELTA .GT. DGXFAC*W(DGGDMX)) GO TO 350
C        ***  DELTA IS SO SMALL WE CANNOT HANDLE THE SPECIAL CASE IN
C        ***  THE AVAILABLE ARITHMETIC.  ACCEPT STEP AS IT IS.
         GO TO 290
C
C  ***  ACCEPTABLE STEP ON FIRST TRY  ***
C
 280  ALPHAK = ZERO
C
C  ***  SUCCESSFUL STEP IN GENERAL.  COMPUTE STEP = -(D**-1)*Q  ***
C
 290  DO 300 I = 1, P
         J = Q0 + I
         STEP(I) = -W(J)/D(I)
 300     CONTINUE
      V(GTSTEP) = -DOTPRD(P, DIG, W(Q))
      V(PREDUC) = HALF * (ABS(ALPHAK)*DST*DST - V(GTSTEP))
      GO TO 430
C
C
C  ***  RESTART WITH NEW RADIUS  ***
C
 310  IF (V(DST0) .LE. ZERO .OR. V(DST0) - RAD .GT. PHIMAX) GO TO 330
C
C     ***  PREPARE TO RETURN NEWTON STEP  ***
C
         RESTRT = .TRUE.
         KA = KA + 1
         K = 0
         DO 320 I = 1, P
              K = K + I
              J = DIAG0 + I
              DIHDI(K) = W(J)
 320          CONTINUE
         UK = NEGONE
         GO TO 40
C
 330  IF (KA .EQ. 0) GO TO 60
C
      DST = W(DSTSAV)
      ALPHAK = ABS(V(STPPAR))
      PHI = DST - RAD
      T = V(DGNORM)/RAD
      IF (RAD .GT. V(RAD0)) GO TO 340
C
C        ***  SMALLER RADIUS  ***
         UK = T - W(EMIN)
         LK = ZERO
         IF (ALPHAK .GT. ZERO) LK = W(LK0)
         LK = MAX(LK, T - W(EMAX))
         IF (V(DST0) .GT. ZERO) LK = MAX(LK, (V(DST0)-RAD)*W(PHIPIN))
         GO TO 260
C
C     ***  BIGGER RADIUS  ***
 340  UK = T - W(EMIN)
      IF (ALPHAK .GT. ZERO) UK = MIN(UK, W(UK0))
      LK = MAX(ZERO, -V(DST0), T - W(EMAX))
      IF (V(DST0) .GT. ZERO) LK = MAX(LK, (V(DST0)-RAD)*W(PHIPIN))
      GO TO 260
C
C  ***  HANDLE (NEARLY) SINGULAR H + ALPHA*D**2  ***
C
C     ***  NEGATE ALPHAK TO INDICATE SPECIAL CASE  ***
 350  ALPHAK = -ALPHAK
C     ***  ALLOCATE STORAGE FOR SCRATCH VECTOR X  ***
      X0 = Q0 + P
      X = X0 + 1
C
C  ***  USE INVERSE POWER METHOD WITH START FROM LSVMIN TO OBTAIN
C  ***  APPROXIMATE EIGENVECTOR CORRESPONDING TO SMALLEST EIGENVALUE
C  ***  OF (D**-1)*H*(D**-1).
C
      DELTA = KAPPA*DELTA
      T = LSVMIN(P, L, W(X), W)
C
      K = 0
C     ***  NORMALIZE W  ***
 360  DO 370 I = 1, P
 370     W(I) = T*W(I)
C     ***  COMPLETE CURRENT INV. POWER ITER. -- REPLACE W BY (L**-T)*W.
      CALL LITVMU(P, W, L, W)
      T1 = ONE/V2NORM(P, W)
      T = T1*T
      IF (T .LE. DELTA) GO TO 390
      IF (K .GT. 30) GO TO 290
      K = K + 1
C     ***  START NEXT INV. POWER ITER. BY STORING NORMALIZED W IN X.
      DO 380 I = 1, P
         J = X0 + I
         W(J) = T1*W(I)
 380     CONTINUE
C     ***  COMPUTE W = (L**-1)*X.
      CALL LIVMUL(P, W, L, W(X))
      T = ONE/V2NORM(P, W)
      GO TO 360
C
 390  DO 400 I = 1, P
 400     W(I) = T1*W(I)
C
C  ***  NOW W IS THE DESIRED APPROXIMATE (UNIT) EIGENVECTOR AND
C  ***  T*X = ((D**-1)*H*(D**-1) + ALPHAK*I)*W.
C
      SW = DOTPRD(P, W(Q), W)
      T1 = (RAD + DST) * (RAD - DST)
      ROOT = SQRT(SW*SW + T1)
      IF (SW .LT. ZERO) ROOT = -ROOT
      SI = T1 / (SW + ROOT)
C     ***  ACCEPT CURRENT STEP IF ADDING SI*W WOULD LEAD TO A
C     ***  FURTHER RELATIVE REDUCTION IN PSI OF LESS THAN V(EPSLON)/3.
      V(PREDUC) = HALF*TWOPSI
      T1 = ZERO
      T = SI*(ALPHAK*SW - HALF*SI*(ALPHAK + T*DOTPRD(P,W(X),W)))
      IF (T .LT. EPSO6*TWOPSI) GO TO 410
         V(PREDUC) = V(PREDUC) + T
         DST = RAD
         T1 = -SI
 410  DO 420 I = 1, P
         J = Q0 + I
         W(J) = T1*W(I) - W(J)
         STEP(I) = W(J) / D(I)
 420     CONTINUE
      V(GTSTEP) = DOTPRD(P, DIG, W(Q))
C
C  ***  SAVE VALUES FOR USE IN A POSSIBLE RESTART  ***
C
 430  V(DSTNRM) = DST
      V(STPPAR) = ALPHAK
      W(LK0) = LK
      W(UK0) = UK
      V(RAD0) = RAD
      W(DSTSAV) = DST
C
C     ***  RESTORE DIAGONAL OF DIHDI  ***
C
      J = 0
      DO 440 I = 1, P
         J = J + I
         K = DIAG0 + I
         DIHDI(J) = W(K)
 440     CONTINUE
      GO TO 999
C
C  ***  SPECIAL CASE -- G = 0  ***
C
 450  V(STPPAR) = ZERO
      V(PREDUC) = ZERO
      V(DSTNRM) = ZERO
      V(GTSTEP) = ZERO
      DO 460 I = 1, P
 460     STEP(I) = ZERO
C
 999  RETURN
      END
      INTEGER FUNCTION IMDCON(K)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  RETURN INTEGER MACHINE DEPENDENT CONSTANTS  ***
C
C     ***  K = 1 MEANS RETURN STANDARD OUTPUT UNIT NUMBER.   ***
C     ***  K = 2 MEANS RETURN ALTERNATE OUTPUT UNIT NUMBER.  ***
C     ***  K = 3 MEANS RETURN  INPUT UNIT NUMBER.            ***
C          (NOTE -- K = 2, 3 ARE USED ONLY BY TEST PROGRAMS.)
C
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   K
C
C  LOCAL ARRAYS
      INTEGER
     +   MDCON(3)
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   I1MACH
      EXTERNAL I1MACH
C
      MDCON(1) = I1MACH(2)
      MDCON(2) = I1MACH(3)
      MDCON(3) = I1MACH(1)
C
      IMDCON = MDCON(K)
      RETURN
      END
      SUBROUTINE ITSMRY(D, IV, P, V, X, NLUPRI)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  PRINT NL2SOL (VERSION 2.2) ITERATION SUMMARY  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   D(P),V(*),X(P)
      INTEGER
     +   IV(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   NRELDF,OLDF,PRELDF,RELDF,ZERO
      INTEGER
     +   COV1,COVMAT,COVPRT,COVREQ,DSTNRM,F,F0,FDIF,G,G1,I,I1,ICH,
     +   II,IV1,J,M,NEEDHD,NF,NFCALL,NFCOV,NG,NGCALL,NGCOV,NITER,
     +   NREDUC,OL,OUTLEV,PREDUC,PRNTIT,PRUNIT,PU,RELDX,SIZE,
     +   SOLPRT,STATPR,STPPAR,SUSED,X0PRT,NLUPRI
C
C  LOCAL ARRAYS
      CHARACTER
     +   MODEL1(3,6)*1,MODEL2(4,6)*1
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER IV(1), P
C     DOUBLE PRECISION D(P), V(1), X(P)
C     DIMENSION IV(*), V(*)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER COV1, G1, I, II, IV1, I1, J, M, NF, NG, OL, PU
C     CHARACTER*1 MODEL1(3, 6), MODEL2(4, 6)
C     DOUBLE PRECISION NRELDF, OLDF, PRELDF, RELDF, ZERO
C
C/
C  ***  NO EXTERNAL FUNCTIONS OR SUBROUTINES  ***
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER COVMAT, COVPRT, COVREQ, DSTNRM, F, FDIF, F0, G,
C    1        NEEDHD, NFCALL, NFCOV, NGCOV, NGCALL, NITER, NREDUC,
C    2        OUTLEV, PREDUC, PRNTIT, PRUNIT, RELDX, SIZE, SOLPRT,
C    3        STATPR, STPPAR, SUSED, X0PRT
C
C  ***  IV SUBSCRIPT VALUES  ***
C
      DATA COVMAT/26/, COVPRT/14/, G/28/, COVREQ/15/,
     +     NEEDHD/39/, NFCALL/6/, NFCOV/40/, NGCOV/41/,
     +     NGCALL/30/, NITER/31/, OUTLEV/19/, PRNTIT/48/,
     +     PRUNIT/21/, SOLPRT/22/, STATPR/23/, SUSED/57/,
     +     X0PRT/24/
C
C  ***  V SUBSCRIPT VALUES  ***
C
      DATA DSTNRM/2/, F/10/, F0/13/, FDIF/11/, NREDUC/6/,
     +     PREDUC/7/, RELDX/17/, SIZE/47/, STPPAR/5/
C
      DATA MODEL1(1, 1), MODEL1(2, 1), MODEL1(3, 1)
     +   /        ' ',          ' ',          ' '  /
      DATA MODEL1(1, 2), MODEL1(2, 2), MODEL1(3, 2)
     +   /        ' ',          ' ',          ' '  /
      DATA MODEL1(1, 3), MODEL1(2, 3), MODEL1(3, 3)
     +   /        ' ',          ' ',          ' '  /
      DATA MODEL1(1, 4), MODEL1(2, 4), MODEL1(3, 4)
     +   /        ' ',          ' ',          ' '  /
      DATA MODEL1(1, 5), MODEL1(2, 5), MODEL1(3, 5)
     +   /        ' ',          'G',          ' '  /
      DATA MODEL1(1, 6), MODEL1(2, 6), MODEL1(3, 6)
     +   /        ' ',          'S',          ' '  /
      DATA MODEL2(1, 1), MODEL2(2, 1), MODEL2(3, 1), MODEL2(4, 1)
     +    /       ' ',          'G',          ' ',          ' '  /
      DATA MODEL2(1, 2), MODEL2(2, 2), MODEL2(3, 2), MODEL2(4, 2)
     +   /        ' ',          'S',          ' ',          ' '  /
      DATA MODEL2(1, 3), MODEL2(2, 3), MODEL2(3, 3), MODEL2(4, 3)
     +   /        'G',          '-',          'S',          ' '  /
      DATA MODEL2(1, 4), MODEL2(2, 4), MODEL2(3, 4), MODEL2(4, 4)
     +   /        'S',          '-',          'G',          ' '  /
      DATA MODEL2(1, 5), MODEL2(2, 5), MODEL2(3, 5), MODEL2(4, 5)
     +   /        '-',          'S',          '-',          'G'  /
      DATA MODEL2(1, 6), MODEL2(2, 6), MODEL2(3, 6), MODEL2(4, 6)
     +   /        '-',          'G',          '-',          'S'  /
      DATA ZERO/0.0D+00/
C
C-----------------------------------------------------------------------
C
CPFP
      PU = IV(PRUNIT)
      IF ((NLUPRI .EQ. 0) .OR. (NLUPRI .NE. IV(PRUNIT))) PU = NLUPRI
Cend-PFP
C
      IF (PU .EQ. 0) GO TO 999
      IV1 = IV(1)
      OL = IV(OUTLEV)
      IF (IV1 .LT. 2 .OR. IV1 .GT. 15) GO TO 140
      IF (OL .EQ. 0) GO TO 20
      IF (IV1 .GE. 12) GO TO 20
      IF (IV1 .GE. 10 .AND. IV(PRNTIT) .EQ. 0) GO TO 20
      IF (IV1 .GT. 2) GO TO 10
         IV(PRNTIT) = IV(PRNTIT) + 1
         IF (IV(PRNTIT) .LT. ABS(OL)) GO TO 999
 10   NF = IV(NFCALL) - ABS(IV(NFCOV))
      IV(PRNTIT) = 0
      RELDF = ZERO
      PRELDF = ZERO
      OLDF = V(F0)
      IF (OLDF .LE. ZERO) GO TO 12
         RELDF = V(FDIF) / OLDF
         PRELDF = V(PREDUC) / OLDF
 12   IF (OL .GT. 0) GO TO 15
C
C        ***  PRINT SHORT SUMMARY LINE  ***
C
         IF (IV(NEEDHD) .EQ. 1) WRITE(PU, 1010)
 1010 FORMAT(12H0   IT    NF,6X,'F',8X,5HRELDF,6X,6HPRELDF,5X,5HRELDX)
         IV(NEEDHD) = 0
         WRITE(PU,1017) IV(NITER), NF, V(F), RELDF, PRELDF, V(RELDX)
         GO TO 20
C
C     ***  PRINT LONG SUMMARY LINE  ***
C
 15   IF (IV(NEEDHD) .EQ. 1) WRITE(PU,1015)
 1015 FORMAT(12H0   IT    NF,6X,'F',8X,5HRELDF,6X,6HPRELDF,5X,5HRELDX,
     +       4X,15HMODEL    STPPAR,6X,4HSIZE,6X,6HD*STEP,5X,7HNPRELDF)
      IV(NEEDHD) = 0
      M = IV(SUSED)
      NRELDF = ZERO
      IF (OLDF .GT. ZERO) NRELDF = V(NREDUC) / OLDF
      WRITE(PU,1017) IV(NITER), NF, V(F), RELDF, PRELDF, V(RELDX),
     +               (MODEL1(ICH, M), ICH = 1, 3),
     +               (MODEL2(ICH, M), ICH = 1, 4),
     +               V(STPPAR), V(SIZE), V(DSTNRM), NRELDF
 1017 FORMAT(1X,I5,I6,4E11.3,7A1,4E11.3)
C
 20   GO TO (999,999,30,35,40,45,50,60,70,80,90,150,110,120,130), IV1
C
 30   WRITE(PU,1030)
 1030 FORMAT(26H0***** X-CONVERGENCE *****)
      GO TO 180
C
 35   WRITE(PU,1035)
 1035 FORMAT(42H0***** RELATIVE FUNCTION CONVERGENCE *****)
      GO TO 180
C
 40   WRITE(PU,1040)
 1040 FORMAT(49H0***** X- AND RELATIVE FUNCTION CONVERGENCE *****)
      GO TO 180
C
 45   WRITE(PU,1045)
 1045 FORMAT(42H0***** ABSOLUTE FUNCTION CONVERGENCE *****)
      GO TO 180
C
 50   WRITE(PU,1050)
 1050 FORMAT(33H0***** SINGULAR CONVERGENCE *****)
      GO TO 180
C
 60   WRITE(PU,1060)
 1060 FORMAT(30H0***** FALSE CONVERGENCE *****)
      GO TO 180
C
 70   WRITE(PU,1070)
 1070 FORMAT(38H0***** FUNCTION EVALUATION LIMIT *****)
      GO TO 180
C
 80   WRITE(PU,1080)
 1080 FORMAT(28H0***** ITERATION LIMIT *****)
      GO TO 180
C
 90   WRITE(PU,1090)
 1090 FORMAT(18H0***** STOPX *****)
      GO TO 180
C
 110  WRITE(PU,1100)
 1100 FORMAT(45H0***** INITIAL SUM OF SQUARES OVERFLOWS *****)
C
      GO TO 150
C
 120  WRITE(PU,1120)
 1120 FORMAT(37H0***** BAD PARAMETERS TO ASSESS *****)
      GO TO 999
C
 130  continue
      write ( pu, * ) ' '
      WRITE(PU,1130)
 1130 FORMAT('***** J COULD NOT BE COMPUTED *****')
      IF (IV(NITER) .GT. 0) GO TO 190
      GO TO 150
C
 140  write ( pu, * ) ' '
      WRITE(PU,1140) IV1
 1140 FORMAT('***** IV(1) =',I5,' *****')
      GO TO 999
C
C  ***  INITIAL CALL ON ITSMRY  ***
C
 150  IF (IV(X0PRT) .NE. 0) then
        write ( pu, * ) ' '
        WRITE(PU,1150) (I, X(I), D(I), I = 1, P)
      end if
 1150 FORMAT('     I     INITIAL X(I)       D(I)'//(1X,I5,E17.6,E14.3))
      IF (IV1 .GE. 13) GO TO 999
      IV(NEEDHD) = 0
      IV(PRNTIT) = 0
      IF (OL .EQ. 0) GO TO 999
      IF (OL .LT. 0) WRITE(PU,1010)
      IF (OL .GT. 0) WRITE(PU,1015)
      write ( pu, * ) ' '
      WRITE(PU,1160) V(F)
 1160 FORMAT('     0     1',E11.3,11X,E11.3)
      GO TO 999
C
C  ***  PRINT VARIOUS INFORMATION REQUESTED ON SOLUTION  ***
C
 180  IV(NEEDHD) = 1
      IF (IV(STATPR) .EQ. 0) GO TO 190
         OLDF = V(F0)
         PRELDF = ZERO
         NRELDF = ZERO
         IF (OLDF .LE. ZERO) GO TO 185
              PRELDF = V(PREDUC) / OLDF
              NRELDF = V(NREDUC) / OLDF
 185     NF = IV(NFCALL) - IV(NFCOV)
         NG = IV(NGCALL) - IV(NGCOV)
         WRITE(PU,1180) V(F), V(RELDX), NF, NG, PRELDF, NRELDF
 1180 FORMAT(9H0FUNCTION,E17.6,8H   RELDX,E20.6/12H FUNC. EVALS,
     +   I8,9X,'GRAD. EVALS',I8/' PRELDF',E19.6,3X,'NPRELDF',E18.6)
C
         IF (IV(NFCOV) .GT. 0) then
           write ( pu, * ) ' '
           WRITE(PU,1185) IV(NFCOV)
         end if
 1185    FORMAT(1x,I4,' EXTRA FUNC. EVALS FOR COVARIANCE.')
         IF (IV(NGCOV) .GT. 0) WRITE(PU,1186) IV(NGCOV)
 1186    FORMAT(1X,I4,' EXTRA GRAD. EVALS FOR COVARIANCE.')

 190  IF (IV(SOLPRT) .EQ. 0) GO TO 210
         IV(NEEDHD) = 1
         G1 = IV(G)
         write ( pu, * ) ' '
         WRITE(PU,1190)
 1190 FORMAT('     I      FINAL X(I)',8X,'D(I)',10X,'G(I)'/)
         DO 200 I = 1, P
              WRITE(PU,1200) I, X(I), D(I), V(G1)
              G1 = G1 + 1
 200          CONTINUE
 1200    FORMAT(1X,I5,E17.6,2E14.3)
C
 210  IF (IV(COVPRT) .EQ. 0) GO TO 999
      COV1 = IV(COVMAT)
      IV(NEEDHD) = 1
      IF (COV1) 220, 230, 240
 220  IF (-1 .EQ. COV1) then
        write ( pu, * ) ' '
        WRITE(PU,1220)
      end if

 1220 FORMAT('++++++ INDEFINITE COVARIANCE MATRIX ++++++')
      IF (-2 .EQ. COV1) then
        write ( pu, * ) ' '
        WRITE(PU,1225)
      end if
 1225 FORMAT('++++++ OVERSIZE STEPS IN COMPUTING COVARIANCE +++++')
      GO TO 999

 230  continue
      write ( pu, * ) ' '
      WRITE(PU,1230)
 1230 FORMAT('++++++ COVARIANCE MATRIX NOT COMPUTED ++++++')
      GO TO 999

 240  I = ABS(IV(COVREQ))

      IF (I .LE. 1) then
        write ( pu, * ) ' '
        WRITE(PU,1241)
        write ( pu, * ) ' '
      end if

 1241 FORMAT('COVARIANCE = SCALE * H**-1 * (J**T * J) * H**-1')

      IF (I .EQ. 2) then
        write ( pu, * ) ' '
        WRITE(PU,1242)
        write ( pu, * ) ' '
      end if

 1242 FORMAT('COVARIANCE = SCALE * H**-1')

      IF (I .GE. 3) then
        write ( pu, * ) ' '
        WRITE(PU,1243)
        write ( pu, * ) ' '
      end if

 1243 FORMAT('COVARIANCE = SCALE * (J**T * J)**-1')

      II = COV1 - 1
      IF (OL .LE. 0) GO TO 260
      DO 250 I = 1, P
         I1 = II + 1
         II = II + I
         WRITE(PU,1250) I, (V(J), J = I1, II)
 250     CONTINUE
 1250 FORMAT(4H ROW,I3,2X,9E12.4/(9X,9E12.4))
      GO TO 999
C
 260  DO 270 I = 1, P
         I1 = II + 1
         II = II + I
         WRITE(PU,1270) I, (V(J), J = I1, II)
 270     CONTINUE
 1270 FORMAT(4H ROW,I3,2X,5E12.4/(9X,5E12.4))
C
 999  RETURN
      END
      SUBROUTINE LINVRT(N, LIN, L)

c*********************************************************************72
c
C  ***  COMPUTE  LIN = L**-1,  BOTH  N X N  LOWER TRIANG. STORED   ***
C  ***  COMPACTLY BY ROWS.  LIN AND L MAY SHARE THE SAME STORAGE.  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   L(*),LIN(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ONE,T,ZERO
      INTEGER
     +   I,II,IM1,J0,J1,JJ,K,K0,NP1
C
C  ***  PARAMETERS  ***
C
C     INTEGER N
C     DOUBLE PRECISION L(*), LIN(*)
C     DIMENSION L(N*(N+1)/2), LIN(N*(N+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER I, II, IM1, JJ, J0, J1, K, K0, NP1
C     DOUBLE PRECISION ONE, T, ZERO
      DATA ONE/1.0D+00/, ZERO/0.0D+00/
C
C  ***  BODY  ***
C
      NP1 = N + 1
      J0 = N*(NP1)/2
      DO 30 II = 1, N
         I = NP1 - II
         LIN(J0) = ONE/L(J0)
         IF (I .LE. 1) GO TO 999
         J1 = J0
         IM1 = I - 1
         DO 20 JJ = 1, IM1
              T = ZERO
              J0 = J1
              K0 = J1 - JJ
              DO 10 K = 1, JJ
                   T = T - L(K0)*LIN(J0)
                   J0 = J0 - 1
                   K0 = K0 + K - I
 10                CONTINUE
              LIN(J0) = T/L(K0)
 20           CONTINUE
         J0 = J0 - 1
 30      CONTINUE
 999  RETURN
      END
      SUBROUTINE LITVMU(N, X, L, Y)

c*********************************************************************72
c
C  ***  SOLVE  (L**T)*X = Y,  WHERE  L  IS AN  N X N  LOWER TRIANGULAR
C  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
C  ***  STORAGE.  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   L(1),X(N),Y(N)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   XI,ZERO
      INTEGER
     +   I,I0,II,IJ,IM1,J,NP1
C
      DATA ZERO/0.0D+00/
C
      DO 10 I = 1, N
 10      X(I) = Y(I)
      NP1 = N + 1
      I0 = N*(N+1)/2
      DO 30 II = 1, N
         I = NP1 - II
         XI = X(I)/L(I0)
         X(I) = XI
         IF (I .LE. 1) GO TO 999
         I0 = I0 - I
         IF (XI .EQ. ZERO) GO TO 30
         IM1 = I - 1
         DO 20 J = 1, IM1
              IJ = I0 + J
              X(J) = X(J) - XI*L(IJ)
 20           CONTINUE
 30      CONTINUE
 999  RETURN
      END
      SUBROUTINE LIVMUL(N, X, L, Y)

c*********************************************************************72
c
C  ***  SOLVE  L*X = Y, WHERE  L  IS AN  N X N  LOWER TRIANGULAR
C  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
C  ***  STORAGE.  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   L(1),X(N),Y(N)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   T,ZERO
      INTEGER
     +   I,J,K
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD
      EXTERNAL DOTPRD
C
      DATA ZERO/0.0D+00/
C
      DO 10 K = 1, N
         IF (Y(K) .NE. ZERO) GO TO 20
         X(K) = ZERO
 10      CONTINUE
      GO TO 999
 20   J = K*(K+1)/2
      X(K) = Y(K) / L(J)
      IF (K .GE. N) GO TO 999
      K = K + 1
      DO 30 I = K, N
         T = DOTPRD(I-1, L(J+1), X)
         J = J + I
         X(I) = (Y(I) - T)/L(J)
 30      CONTINUE
 999  RETURN
      END
      SUBROUTINE LMSTEP(D, G, IERR, IPIVOT, KA, P, QTR, R, STEP, V, W)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  COMPUTE LEVENBERG-MARQUARDT STEP USING MORE HEBDEN TECHNIQUE  **
C  ***  NL2SOL VERSION 2.2.  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IERR,KA,P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   D(P),G(P),QTR(P),R(1),STEP(P),V(21),W(1)
      INTEGER
     +   IPIVOT(P)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   A,ADI,ALPHAK,B,D1,D2,DFAC,DFACSQ,DST,DTOL,EIGHT,HALF,LK,
     +   NEGONE,OLDPHI,ONE,P001,PHI,PHIMAX,PHIMIN,PSIFAC,RAD,SI,SJ,
     +   SQRTAK,T,THREE,TTOL,TWOPSI,UK,WL,ZERO
      INTEGER
     +   DGNORM,DST0,DSTNRM,DSTSAV,EPSLON,GTSTEP,I,I1,IP1,J1,K,
     +   KALIM,L,LK0,NREDUC,PHIPIN,PHMNFC,PHMXFC,PP1O2,PREDUC,RAD0,
     +   RADIUS,RES,RES0,RMAT,RMAT0,STPPAR,UK0
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD,V2NORM
      EXTERNAL DOTPRD,V2NORM
C
C  EXTERNAL SUBROUTINES
      EXTERNAL LITVMU,LIVMUL,VCOPY
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX,MIN,SQRT
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER IERR, KA, P
C     INTEGER IPIVOT(P)
C     DOUBLE PRECISION D(P), G(P), QTR(P), R(1), STEP(P), V(21), W(1)
C     DIMENSION W(P*(P+5)/2 + 4)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  PURPOSE  ***
C
C        GIVEN THE R MATRIX FROM THE QR DECOMPOSITION OF A JACOBIAN
C     MATRIX, J, AS WELL AS Q-TRANSPOSE TIMES THE CORRESPONDING
C     RESIDUAL VECTOR, RESID, THIS SUBROUTINE COMPUTES A LEVENBERG
C     MARQUARDT STEP OF APPROXIMATE LENGTH V(RADIUS) BY THE MORE
C     TECHNIQUE.
C
C  ***  PARAMETER DESCRIPTION  ***
C
C      D (IN)  = THE SCALE VECTOR.
C      G (IN)  = THE GRADIENT VECTOR (J**T)*R.
C   IERR (I/O) = RETURN CODE FROM QRFACT OR QRFGS -- 0 MEANS R HAS
C             FULL RANK.
C IPIVOT (I/O) = PERMUTATION ARRAY FROM QRFACT OR QRFGS, WHICH COMPUTE
C             QR DECOMPOSITIONS WITH COLUMN PIVOTING.
C     KA (I/O).  KA .LT. 0 ON INPUT MEANS THIS IS THE FIRST CALL ON
C             LMSTEP FOR THE CURRENT R AND QTR.  ON OUTPUT KA CON-
C             TAINS THE NUMBER OF HEBDEN ITERATIONS NEEDED TO DETERMINE
C             STEP.  KA = 0 MEANS A GAUSS-NEWTON STEP.
C      P (IN)  = NUMBER OF PARAMETERS.
C    QTR (IN)  = (Q**T)*RESID = Q-TRANSPOSE TIMES THE RESIDUAL VECTOR.
C      R (IN)  = THE R MATRIX, STORED COMPACTLY BY COLUMNS.
C   STEP (OUT) = THE LEVENBERG-MARQUARDT STEP COMPUTED.
C      V (I/O) CONTAINS VARIOUS CONSTANTS AND VARIABLES DESCRIBED BELOW.
C      W (I/O) = WORKSPACE OF LENGTH P*(P+5)/2 + 4.
C
C  ***  ENTRIES IN V  ***
C
C V(DGNORM) (I/O) = 2-NORM OF (D**-1)*G.
C V(DSTNRM) (I/O) = 2-NORM OF D*STEP.
C V(DST0)   (I/O) = 2-NORM OF GAUSS-NEWTON STEP (FOR NONSING. J).
C V(EPSLON) (IN) = MAX. REL. ERROR ALLOWED IN TWONORM(R)**2 MINUS
C             TWONORM(R - J*STEP)**2.  (SEE ALGORITHM NOTES BELOW.)
C V(GTSTEP) (OUT) = INNER PRODUCT BETWEEN G AND STEP.
C V(NREDUC) (OUT) = HALF THE REDUCTION IN THE SUM OF SQUARES PREDICTED
C             FOR A GAUSS-NEWTON STEP.
C V(PHMNFC) (IN)  = TOL. (TOGETHER WITH V(PHMXFC)) FOR ACCEPTING STEP
C             (MORE*S SIGMA).  THE ERROR V(DSTNRM) - V(RADIUS) MUST LIE
C             BETWEEN V(PHMNFC)*V(RADIUS) AND V(PHMXFC)*V(RADIUS).
C V(PHMXFC) (IN)  (SEE V(PHMNFC).)
C V(PREDUC) (OUT) = HALF THE REDUCTION IN THE SUM OF SQUARES PREDICTED
C             BY THE STEP RETURNED.
C V(RADIUS) (IN)  = RADIUS OF CURRENT (SCALED) TRUST REGION.
C V(RAD0)   (I/O) = VALUE OF V(RADIUS) FROM PREVIOUS CALL.
C V(STPPAR) (I/O) = MARQUARDT PARAMETER (OR ITS NEGATIVE IF THE SPECIAL
C             CASE MENTIONED BELOW IN THE ALGORITHM NOTES OCCURS).
C
C NOTE -- SEE DATA STATEMENT BELOW FOR VALUES OF ABOVE SUBSCRIPTS.
C
C  ***  USAGE NOTES  ***
C
C     IF IT IS DESIRED TO RECOMPUTE STEP USING A DIFFERENT VALUE OF
C     V(RADIUS), THEN THIS ROUTINE MAY BE RESTARTED BY CALLING IT
C     WITH ALL PARAMETERS UNCHANGED EXCEPT V(RADIUS).  (THIS EXPLAINS
C     WHY MANY PARAMETERS ARE LISTED AS I/O).  ON AN INTIIAL CALL (ONE
C     WITH KA = -1), THE CALLER NEED ONLY HAVE INITIALIZED D, G, KA, P,
C     QTR, R, V(EPSLON), V(PHMNFC), V(PHMXFC), V(RADIUS), AND V(RAD0).
C
C  ***  APPLICATION AND USAGE RESTRICTIONS  ***
C
C     THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR LEAST-
C     SQUARES) PACKAGE (REF. 1).
C
C  ***  ALGORITHM NOTES  ***
C
C     THIS CODE IMPLEMENTS THE STEP COMPUTATION SCHEME DESCRIBED IN
C     REFS. 2 AND 4.  FAST GIVENS TRANSFORMATIONS (SEE REF. 3, PP. 60-
C     62) ARE USED TO COMPUTE STEP WITH A NONZERO MARQUARDT PARAMETER.
C        A SPECIAL CASE OCCURS IF J IS (NEARLY) SINGULAR AND V(RADIUS)
C     IS SUFFICIENTLY LARGE.  IN THIS CASE THE STEP RETURNED IS SUCH
C     THAT  TWONORM(R)**2 - TWONORM(R - J*STEP)**2  DIFFERS FROM ITS
C     OPTIMAL VALUE BY LESS THAN V(EPSLON) TIMES THIS OPTIMAL VALUE,
C     WHERE J AND R DENOTE THE ORIGINAL JACOBIAN AND RESIDUAL.  (SEE
C     REF. 2 FOR MORE DETAILS.)
C
C  ***  FUNCTIONS AND SUBROUTINES CALLED  ***
C
C DOTPRD - RETURNS INNER PRODUCT OF TWO VECTORS.
C LITVMU - APPLY INVERSE TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX.
C LIVMUL - APPLY INVERSE OF COMPACT LOWER TRIANG. MATRIX.
C VCOPY  - COPIES ONE VECTOR TO ANOTHER.
C V2NORM - RETURNS 2-NORM OF A VECTOR.
C
C  ***  REFERENCES  ***
C
C 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1980), AN ADAPTIVE
C             NONLINEAR LEAST-SQUARES ALGORITHM, (SUBMITTED TO ACM
C             TRANS. MATH. SOFTWARE).
C 2.  GAY, D.M. (1979), COMPUTING OPTIMAL ELLIPTICALLY CONSTRAINED
C             STEPS, MRC TECH. SUMMARY REPORT NO. 2013, MATH RESEARCH
C             CENTER, UNIV. OF WISCONSIN-MADISON.
C 3.  LAWSON, C.L., AND HANSON, R.J. (1974), SOLVING LEAST SQUARES
C             PROBLEMS, PRENTICE HALL, ENGLEWOOD CLIFFS, N.J.
C 4.  MORE, J.J. (1978), THE LEVENBERG-MARQUARDT ALGORITHM, IMPLEMEN-
C             TATION AND THEORY, PP.105-116 OF SPRINGER LECTURE NOTES
C             IN MATHEMATICS NO. 630, EDITED BY G.A. WATSON, SPRINGER-
C             VERLAG, BERLIN AND NEW YORK.
C
C  ***  GENERAL  ***
C
C     CODED BY DAVID M. GAY.
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
C     MCS-7906671.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER DSTSAV, I, IP1, I1, J1, K, KALIM, L, LK0, PHIPIN,
C    1        PP1O2, RES, RES0, RMAT, RMAT0, UK0
C     DOUBLE PRECISION A, ADI, ALPHAK, B, DFACSQ, DST, DTOL, D1, D2,
C    1                 LK, OLDPHI, PHI, PHIMAX, PHIMIN, PSIFAC, RAD,
C    2                 SI, SJ, SQRTAK, T, TWOPSI, UK, WL
C
C     ***  CONSTANTS  ***
C     DOUBLE PRECISION DFAC, EIGHT, HALF, NEGONE, ONE, P001, THREE,
C    1                 TTOL, ZERO
C
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL DOTPRD, LITVMU, LIVMUL, VCOPY, V2NORM
C     DOUBLE PRECISION DOTPRD, V2NORM
C
C  ***  SUBSCRIPTS FOR V  ***
C
C     INTEGER DGNORM, DSTNRM, DST0, EPSLON, GTSTEP, NREDUC, PHMNFC,
C    1        PHMXFC, PREDUC, RADIUS, RAD0, STPPAR
      DATA DGNORM/1/, DSTNRM/2/, DST0/3/, EPSLON/19/,
     +     GTSTEP/4/, NREDUC/6/, PHMNFC/20/,
     +     PHMXFC/21/, PREDUC/7/, RADIUS/8/,
     +     RAD0/9/, STPPAR/5/
C
      DATA DFAC/256.0D+00/
      DATA EIGHT/8.0D+00/
      DATA HALF/0.5D+00/
      DATA NEGONE/-1.0D+00/
      DATA ONE/1.0D+00/
      DATA P001/1.0D-03/
      DATA THREE/3.0D+00/
      DATA TTOL/2.5D+00/
      DATA ZERO/0.0D+00/
C
C  ***  BODY  ***
C
C     ***  FOR USE IN RECOMPUTING STEP, THE FINAL VALUES OF LK AND UK,
C     ***  THE INVERSE DERIVATIVE OF MORE*S PHI AT 0 (FOR NONSING. J)
C     ***  AND THE VALUE RETURNED AS V(DSTNRM) ARE STORED AT W(LK0),
C     ***  W(UK0), W(PHIPIN), AND W(DSTSAV) RESPECTIVELY.
      ALPHAK = 0.0D+00
      PSIFAC = 0.0D+00
      LK0 = P + 1
      PHIPIN = LK0 + 1
      UK0 = PHIPIN + 1
      DSTSAV = UK0 + 1
      RMAT0 = DSTSAV
C     ***  A COPY OF THE R-MATRIX FROM THE QR DECOMPOSITION OF J IS
C     ***  STORED IN W STARTING AT W(RMAT), AND A COPY OF THE RESIDUAL
C     ***  VECTOR IS STORED IN W STARTING AT W(RES).  THE LOOPS BELOW
C     ***  THAT UPDATE THE QR DECOMP. FOR A NONZERO MARQUARDT PARAMETER
C     ***  WORK ON THESE COPIES.
      RMAT = RMAT0 + 1
      PP1O2 = P * (P + 1) / 2
      RES0 = PP1O2 + RMAT0
      RES = RES0 + 1
      RAD = V(RADIUS)
      IF (RAD .GT. ZERO)
     +   PSIFAC = V(EPSLON)/((EIGHT*(V(PHMNFC) + ONE) + THREE) * RAD**2)
      PHIMAX = V(PHMXFC) * RAD
      PHIMIN = V(PHMNFC) * RAD
C     ***  DTOL, DFAC, AND DFACSQ ARE USED IN RESCALING THE FAST GIVENS
C     ***  REPRESENTATION OF THE UPDATED QR DECOMPOSITION.
      DTOL = ONE/DFAC
      DFACSQ = DFAC*DFAC
C     ***  OLDPHI IS USED TO DETECT LIMITS OF NUMERICAL ACCURACY.  IF
C     ***  WE RECOMPUTE STEP AND IT DOES NOT CHANGE, THEN WE ACCEPT IT.
      OLDPHI = ZERO
      LK = ZERO
      UK = ZERO
      KALIM = KA + 12
C
C  ***  START OR RESTART, DEPENDING ON KA  ***
C
      IF (KA) 10, 20, 370
C
C  ***  FRESH START -- COMPUTE V(NREDUC)  ***
C
 10   KA = 0
      KALIM = 12
      K = P
      IF (IERR .NE. 0) K = ABS(IERR) - 1
      V(NREDUC) = HALF*DOTPRD(K, QTR, QTR)
C
C  ***  SET UP TO TRY INITIAL GAUSS-NEWTON STEP  ***
C
 20   V(DST0) = NEGONE
      IF (IERR .NE. 0) GO TO 90
C
C  ***  COMPUTE GAUSS-NEWTON STEP  ***
C
C     ***  NOTE -- THE R-MATRIX IS STORED COMPACTLY BY COLUMNS IN
C     ***  R(1), R(2), R(3), ...  IT IS THE TRANSPOSE OF A
C     ***  LOWER TRIANGULAR MATRIX STORED COMPACTLY BY ROWS, AND WE
C     ***  TREAT IT AS SUCH WHEN USING LITVMU AND LIVMUL.
      CALL LITVMU(P, W, R, QTR)
C     ***  TEMPORARILY STORE PERMUTED -D*STEP IN STEP.
      DO 60 I = 1, P
         J1 = IPIVOT(I)
         STEP(I) = D(J1)*W(I)
 60      CONTINUE
      DST = V2NORM(P, STEP)
      V(DST0) = DST
      PHI = DST - RAD
      IF (PHI .LE. PHIMAX) GO TO 410
C     ***  IF THIS IS A RESTART, GO TO 110  ***
      IF (KA .GT. 0) GO TO 110
C
C  ***  GAUSS-NEWTON STEP WAS UNACCEPTABLE.  COMPUTE L0  ***
C
      DO 70 I = 1, P
         J1 = IPIVOT(I)
         STEP(I) = D(J1)*(STEP(I)/DST)
 70      CONTINUE
      CALL LIVMUL(P, STEP, R, STEP)
      T = ONE / V2NORM(P, STEP)
      W(PHIPIN) = (T/DST)*T
      LK = PHI*W(PHIPIN)
C
C  ***  COMPUTE U0  ***
C
 90   DO 100 I = 1, P
 100     W(I) = G(I)/D(I)
      V(DGNORM) = V2NORM(P, W)
      UK = V(DGNORM)/RAD
      IF (UK .LE. ZERO) GO TO 390
C
C     ***  ALPHAK WILL BE USED AS THE CURRENT MARQUARDT PARAMETER.  WE
C     ***  USE MORE*S SCHEME FOR INITIALIZING IT.
      ALPHAK = ABS(V(STPPAR)) * V(RAD0)/RAD
C
C
C  ***  TOP OF LOOP -- INCREMENT KA, COPY R TO RMAT, QTR TO RES  ***
C
 110  KA = KA + 1
      CALL VCOPY(PP1O2, W(RMAT), R)
      CALL VCOPY(P, W(RES), QTR)
C
C  ***  SAFEGUARD ALPHAK AND INITIALIZE FAST GIVENS SCALE VECTOR.  ***
C
      IF (ALPHAK .LE. ZERO .OR. ALPHAK .LT. LK .OR. ALPHAK .GE. UK)
     +             ALPHAK = UK * MAX(P001, SQRT(LK/UK))
      SQRTAK = SQRT(ALPHAK)
      DO 120 I = 1, P
 120     W(I) = ONE
C
C  ***  ADD ALPHAK*D AND UPDATE QR DECOMP. USING FAST GIVENS TRANS.  ***
C
      DO 270 I = 1, P
C        ***  GENERATE, APPLY 1ST GIVENS TRANS. FOR ROW I OF ALPHAK*D.
C        ***  (USE STEP TO STORE TEMPORARY ROW)  ***
         L = I*(I+1)/2 + RMAT0
         WL = W(L)
         D2 = ONE
         D1 = W(I)
         J1 = IPIVOT(I)
         ADI = SQRTAK*D(J1)
         IF (ADI .GE. ABS(WL)) GO TO 150
 130     A = ADI/WL
         B = D2*A/D1
         T = A*B + ONE
         IF (T .GT. TTOL) GO TO 150
         W(I) = D1/T
         D2 = D2/T
         W(L) = T*WL
         A = -A
         DO 140 J1 = I, P
              L = L + J1
              STEP(J1) = A*W(L)
 140          CONTINUE
         GO TO 170
C
 150     B = WL/ADI
         A = D1*B/D2
         T = A*B + ONE
         IF (T .GT. TTOL) GO TO 130
         W(I) = D2/T
         D2 = D1/T
         W(L) = T*ADI
         DO 160 J1 = I, P
              L = L + J1
              WL = W(L)
              STEP(J1) = -WL
              W(L) = A*WL
 160          CONTINUE
C
 170     IF (I .EQ. P) GO TO 280
C
C        ***  NOW USE GIVENS TRANS. TO ZERO ELEMENTS OF TEMP. ROW  ***
C
         IP1 = I + 1
         DO 260 I1 = IP1, P
              L = I1*(I1+1)/2 + RMAT0
              WL = W(L)
              SI = STEP(I1-1)
              D1 = W(I1)
C
C             ***  RESCALE ROW I1 IF NECESSARY  ***
C
              IF (D1 .GE. DTOL) GO TO 190
                   D1 = D1*DFACSQ
                   WL = WL/DFAC
                   K = L
                   DO 180 J1 = I1, P
                        K = K + J1
                        W(K) = W(K)/DFAC
 180                    CONTINUE
C
C             ***  USE GIVENS TRANS. TO ZERO NEXT ELEMENT OF TEMP. ROW
C
 190          IF (ABS(SI) .GT. ABS(WL)) GO TO 220
              IF (SI .EQ. ZERO) GO TO 260
 200          A = SI/WL
              B = D2*A/D1
              T = A*B + ONE
              IF (T .GT. TTOL) GO TO 220
              W(L) = T*WL
              W(I1) = D1/T
              D2 = D2/T
              DO 210 J1 = I1, P
                   L = L + J1
                   WL = W(L)
                   SJ = STEP(J1)
                   W(L) = WL + B*SJ
                   STEP(J1) = SJ - A*WL
 210               CONTINUE
              GO TO 240
C
 220          B = WL/SI
              A = D1*B/D2
              T = A*B + ONE
              IF (T .GT. TTOL) GO TO 200
              W(I1) = D2/T
              D2 = D1/T
              W(L) = T*SI
              DO 230 J1 = I1, P
                   L = L + J1
                   WL = W(L)
                   SJ = STEP(J1)
                   W(L) = A*WL + SJ
                   STEP(J1) = B*SJ - WL
 230               CONTINUE
C
C             ***  RESCALE TEMP. ROW IF NECESSARY  ***
C
 240          IF (D2 .GE. DTOL) GO TO 260
                   D2 = D2*DFACSQ
                   DO 250 K = I1, P
 250                    STEP(K) = STEP(K)/DFAC
 260          CONTINUE
 270     CONTINUE
C
C  ***  COMPUTE STEP  ***
C
 280  CALL LITVMU(P, W(RES), W(RMAT), W(RES))
C     ***  RECOVER STEP AND STORE PERMUTED -D*STEP AT W(RES)  ***
      DO 290 I = 1, P
         J1 = IPIVOT(I)
         K = RES0 + I
         T = W(K)
         STEP(J1) = -T
         W(K) = T*D(J1)
 290     CONTINUE
      DST = V2NORM(P, W(RES))
      PHI = DST - RAD
      IF (PHI .LE. PHIMAX .AND. PHI .GE. PHIMIN) GO TO 430
      IF (OLDPHI .EQ. PHI) GO TO 430
      OLDPHI = PHI
C
C  ***  CHECK FOR (AND HANDLE) SPECIAL CASE  ***
C
      IF (PHI .GT. ZERO) GO TO 310
         IF (KA .GE. KALIM) GO TO 430
              TWOPSI = ALPHAK*DST*DST - DOTPRD(P, STEP, G)
              IF (ALPHAK .GE. TWOPSI*PSIFAC) GO TO 310
                   V(STPPAR) = -ALPHAK
                   GO TO 440
C
C  ***  UNACCEPTABLE STEP -- UPDATE LK, UK, ALPHAK, AND TRY AGAIN  ***
C
 300  IF (PHI .LT. ZERO) UK = MIN(UK, ALPHAK)
      GO TO 320
 310  IF (PHI .LT. ZERO) UK = ALPHAK
 320  DO 330 I = 1, P
         J1 = IPIVOT(I)
         K = RES0 + I
         STEP(I) = D(J1) * (W(K)/DST)
 330     CONTINUE
      CALL LIVMUL(P, STEP, W(RMAT), STEP)
      DO 340 I = 1, P
 340     STEP(I) = STEP(I) / SQRT(W(I))
      T = ONE / V2NORM(P, STEP)
      ALPHAK = ALPHAK + T*PHI*T/RAD
      LK = MAX(LK, ALPHAK)
      GO TO 110
C
C  ***  RESTART  ***
C
 370  LK = W(LK0)
      UK = W(UK0)
      IF (V(DST0) .GT. ZERO .AND. V(DST0) - RAD .LE. PHIMAX) GO TO 20
      ALPHAK = ABS(V(STPPAR))
      DST = W(DSTSAV)
      PHI = DST - RAD
      T = V(DGNORM)/RAD
      IF (RAD .GT. V(RAD0)) GO TO 380
C
C        ***  SMALLER RADIUS  ***
         UK = T
         IF (ALPHAK .LE. ZERO) LK = ZERO
         IF (V(DST0) .GT. ZERO) LK = MAX(LK, (V(DST0)-RAD)*W(PHIPIN))
         GO TO 300
C
C     ***  BIGGER RADIUS  ***
 380  IF (ALPHAK .LE. ZERO .OR. UK .GT. T) UK = T
      LK = ZERO
      IF (V(DST0) .GT. ZERO) LK = MAX(LK, (V(DST0)-RAD)*W(PHIPIN))
      GO TO 300
C
C  ***  SPECIAL CASE -- RAD .LE. 0 OR (G = 0 AND J IS SINGULAR)  ***
C
 390  V(STPPAR) = ZERO
      DST = ZERO
      LK = ZERO
      UK = ZERO
      V(GTSTEP) = ZERO
      V(PREDUC) = ZERO
      DO 400 I = 1, P
 400     STEP(I) = ZERO
      GO TO 450
C
C  ***  ACCEPTABLE GAUSS-NEWTON STEP -- RECOVER STEP FROM W  ***
C
 410  ALPHAK = ZERO
      DO 420 I = 1, P
         J1 = IPIVOT(I)
         STEP(J1) = -W(I)
 420     CONTINUE
C
C  ***  SAVE VALUES FOR USE IN A POSSIBLE RESTART  ***
C
 430  V(STPPAR) = ALPHAK
 440  V(GTSTEP) = DOTPRD(P, STEP, G)
      V(PREDUC) = HALF * (ALPHAK*DST*DST - V(GTSTEP))
 450  V(DSTNRM) = DST
      W(DSTSAV) = DST
      W(LK0) = LK
      W(UK0) = UK
      V(RAD0) = RAD
C
      RETURN
      END
      SUBROUTINE LSQRT(N1, N, L, A, IRC)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  COMPUTE ROWS N1 THROUGH N OF THE CHOLESKY FACTOR  L  OF
C  ***  A = L*(L**T),  WHERE  L  AND THE LOWER TRIANGLE OF  A  ARE BOTH
C  ***  STORED COMPACTLY BY ROWS (AND MAY OCCUPY THE SAME STORAGE).
C  ***  IRC = 0 MEANS ALL WENT WELL.  IRC = J MEANS THE LEADING
C  ***  PRINCIPAL  J X J  SUBMATRIX OF  A  IS NOT POSITIVE DEFINITE --
C  ***  AND  L(J*(J+1)/2)  CONTAINS THE (NONPOS.) REDUCED J-TH DIAGONAL.
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IRC,N,N1
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   A(1),L(1)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   T,TD,ZERO
      INTEGER
     +   I,I0,IJ,IK,IM1,J,J0,JK,JM1,K
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C  ***  PARAMETERS  ***
C
C     INTEGER N1, N, IRC
C     DOUBLE PRECISION L(1), A(1)
C     DIMENSION L(N*(N+1)/2), A(N*(N+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER I, IJ, IK, IM1, I0, J, JK, JM1, J0, K
C     DOUBLE PRECISION T, TD, ZERO
C
C/
      DATA ZERO/0.0D+00/
C
C  ***  BODY  ***
C
      I0 = N1 * (N1 - 1) / 2
      DO 50 I = N1, N
         TD = ZERO
         IF (I .EQ. 1) GO TO 40
         J0 = 0
         IM1 = I - 1
         DO 30 J = 1, IM1
              T = ZERO
              IF (J .EQ. 1) GO TO 20
              JM1 = J - 1
              DO 10 K = 1, JM1
                   IK = I0 + K
                   JK = J0 + K
                   T = T + L(IK)*L(JK)
 10                CONTINUE
 20           IJ = I0 + J
              J0 = J0 + J
              T = (A(IJ) - T) / L(J0)
              L(IJ) = T
              TD = TD + T*T
 30           CONTINUE
 40      I0 = I0 + I
         T = A(I0) - TD
         IF (T .LE. ZERO) GO TO 60
         L(I0) = SQRT(T)
 50      CONTINUE
C
      IRC = 0
      GO TO 999
C
 60   L(I0) = T
      IRC = I
C
 999  RETURN
      END
      DOUBLE PRECISION FUNCTION LSVMIN(P, L, X, Y)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  ESTIMATE SMALLEST SING. VALUE OF PACKED LOWER TRIANG. MATRIX L
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   L(1),X(P),Y(P)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   B,HALF,ONE,PSJ,R9973,SMINUS,SPLUS,T,XMINUS,XPLUS,ZERO
      INTEGER
     +   I,II,IX,J,J0,JI,JJ,JJJ,JM1,PPLUS1
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   V2NORM
      EXTERNAL V2NORM
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MOD
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER P
C     DOUBLE PRECISION L(1), X(P), Y(P)
C     DIMENSION L(P*(P+1)/2)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  PURPOSE  ***
C
C     THIS FUNCTION RETURNS A GOOD OVER-ESTIMATE OF THE SMALLEST
C     SINGULAR VALUE OF THE PACKED LOWER TRIANGULAR MATRIX L.
C
C  ***  PARAMETER DESCRIPTION  ***
C
C  P (IN)  = THE ORDER OF L.  L IS A  P X P  LOWER TRIANGULAR MATRIX.
C  L (IN)  = ARRAY HOLDING THE ELEMENTS OF  L  IN ROW ORDER, I.E.
C             L(1,1), L(2,1), L(2,2), L(3,1), L(3,2), L(3,3), ETC.
C  X (OUT) IF LSVMIN RETURNS A POSITIVE VALUE, THEN X IS A NORMALIZED
C             APPROXIMATE LEFT SINGULAR VECTOR CORRESPONDING TO THE
C             SMALLEST SINGULAR VALUE.  THIS APPROXIMATION MAY BE VERY
C             CRUDE.  IF LSVMIN RETURNS ZERO, THEN SOME COMPONENTS OF X
C             ARE ZERO AND THE REST RETAIN THEIR INPUT VALUES.
C  Y (OUT) IF LSVMIN RETURNS A POSITIVE VALUE, THEN Y = (L**-1)*X IS AN
C             UNNORMALIZED APPROXIMATE RIGHT SINGULAR VECTOR CORRESPOND-
C             ING TO THE SMALLEST SINGULAR VALUE.  THIS APPROXIMATION
C             MAY BE CRUDE.  IF LSVMIN RETURNS ZERO, THEN Y RETAINS ITS
C             INPUT VALUE.  THE CALLER MAY PASS THE SAME VECTOR FOR X
C             AND Y (NONSTANDARD FORTRAN USAGE), IN WHICH CASE Y OVER-
C             WRITES X (FOR NONZERO LSVMIN RETURNS).
C
C  ***  APPLICATION AND USAGE RESTRICTIONS  ***
C
C     THERE ARE NO USAGE RESTRICTIONS.
C
C  ***  ALGORITHM NOTES  ***
C
C     THE ALGORITHM IS BASED ON (1), WITH THE ADDITIONAL PROVISION THAT
C     LSVMIN = 0 IS RETURNED IF THE SMALLEST DIAGONAL ELEMENT OF L
C     (IN MAGNITUDE) IS NOT MORE THAN THE UNIT ROUNDOFF TIMES THE
C     LARGEST.  THE ALGORITHM USES A RANDOM NUMBER GENERATOR PROPOSED
C     IN (4), WHICH PASSES THE SPECTRAL TEST WITH FLYING COLORS -- SEE
C     (2) AND (3).
C
C  ***  SUBROUTINES AND FUNCTIONS CALLED  ***
C
C        V2NORM - FUNCTION, RETURNS THE 2-NORM OF A VECTOR.
C
C  ***  REFERENCES  ***
C
C     (1) CLINE, A., MOLER, C., STEWART, G., AND WILKINSON, J.H.(1977),
C         AN ESTIMATE FOR THE CONDITION NUMBER OF A MATRIX, REPORT
C         TM-310, APPLIED MATH. DIV., ARGONNE NATIONAL LABORATORY.
C
C     (2) HOAGLIN, D.C. (1976), THEORETICAL PROPERTIES OF CONGRUENTIAL
C         RANDOM-NUMBER GENERATORS --  AN EMPIRICAL VIEW,
C         MEMORANDUM NS-340, DEPT. OF STATISTICS, HARVARD UNIV.
C
C     (3) KNUTH, D.E. (1969), THE ART OF COMPUTER PROGRAMMING, VOL. 2
C         (SEMINUMERICAL ALGORITHMS), ADDISON-WESLEY, READING, MASS.
C
C     (4) SMITH, C.S. (1971), MULTIPLICATIVE PSEUDO-RANDOM NUMBER
C         GENERATORS WITH PRIME MODULUS, J. ASSOC. COMPUT. MACH. 18,
C         PP. 586-593.
C
C  ***  HISTORY  ***
C
C     DESIGNED AND CODED BY DAVID M GAY (WINTER 1977/SUMMER 1978).
C
C  ***  GENERAL  ***
C
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, AND MCS76-11989.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER I, II, IX, J, JI, JJ, JJJ, JM1, J0, PPLUS1
C     DOUBLE PRECISION B, PSJ, SMINUS, SPLUS, T, XMINUS, XPLUS
C
C  ***  CONSTANTS  ***
C
C     DOUBLE PRECISION HALF, ONE, R9973, ZERO
C
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL V2NORM
C     DOUBLE PRECISION V2NORM
C
      DATA IX/2/
      DATA HALF/0.5D+00/, ONE/1.0D+00/, R9973/9973.0D+00/, ZERO/0.0D+00/
C
C  ***  BODY  ***
C
C  ***  FIRST CHECK WHETHER TO RETURN LSVMIN = 0 AND INITIALIZE X  ***
C
      II = 0
      DO 10 I = 1, P
         X(I) = ZERO
         II = II + I
         IF (L(II) .EQ. ZERO) GO TO 300
 10      CONTINUE
      IF (MOD(IX, 9973) .EQ. 0) IX = 2
      PPLUS1 = P + 1
C
C  ***  SOLVE (L**T)*X = B, WHERE THE COMPONENTS OF B HAVE RANDOMLY
C  ***  CHOSEN MAGNITUDES IN (.5,1) WITH SIGNS CHOSEN TO MAKE X LARGE.
C
C     DO J = P TO 1 BY -1...
      DO 100 JJJ = 1, P
         J = PPLUS1 - JJJ
C       ***  DETERMINE X(J) IN THIS ITERATION. NOTE FOR I = 1,2,...,J
C       ***  THAT X(I) HOLDS THE CURRENT PARTIAL SUM FOR ROW I.
         IX = MOD(3432*IX, 9973)
         B = HALF*(ONE + IX/R9973)
         XPLUS = (B - X(J))
         XMINUS = (-B - X(J))
         SPLUS = ABS(XPLUS)
         SMINUS = ABS(XMINUS)
         JM1 = J - 1
         J0 = J*JM1/2
         JJ = J0 + J
         XPLUS = XPLUS/L(JJ)
         XMINUS = XMINUS/L(JJ)
         IF (JM1 .EQ. 0) GO TO 30
         DO 20 I = 1, JM1
              JI = J0 + I
              SPLUS = SPLUS + ABS(X(I) + L(JI)*XPLUS)
              SMINUS = SMINUS + ABS(X(I) + L(JI)*XMINUS)
 20           CONTINUE
 30      IF (SMINUS .GT. SPLUS) XPLUS = XMINUS
         X(J) = XPLUS
C       ***  UPDATE PARTIAL SUMS  ***
         IF (JM1 .EQ. 0) GO TO 100
         DO 40 I = 1, JM1
              JI = J0 + I
              X(I) = X(I) + L(JI)*XPLUS
 40           CONTINUE
 100     CONTINUE
C
C  ***  NORMALIZE X  ***
C
      T = ONE/V2NORM(P, X)
      DO 110 I = 1, P
 110     X(I) = T*X(I)
C
C  ***  SOLVE L*Y = X AND RETURN SVMIN = 1/TWONORM(Y)  ***
C
      DO 200 J = 1, P
         PSJ = ZERO
         JM1 = J - 1
         J0 = J*JM1/2
         IF (JM1 .EQ. 0) GO TO 130
         DO 120 I = 1, JM1
              JI = J0 + I
              PSJ = PSJ + L(JI)*Y(I)
 120          CONTINUE
 130     JJ = J0 + J
         Y(J) = (X(J) - PSJ)/L(JJ)
 200     CONTINUE
C
      LSVMIN = ONE/V2NORM(P, Y)
      GO TO 999
C
 300  LSVMIN = ZERO
 999  RETURN
      END
      SUBROUTINE LTSQAR(N, A, L)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  SET A TO LOWER TRIANGLE OF (L**T) * L  ***
C
C  ***  L = N X N LOWER TRIANG. MATRIX STORED ROWWISE.  ***
C  ***  A IS ALSO STORED ROWWISE AND MAY SHARE STORAGE WITH L.  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   A(*),L(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   LII,LJ
      INTEGER
     +   I,I1,II,IIM1,J,K,M
C
C     INTEGER N
C     DOUBLE PRECISION A(1), L(1)
C     DIMENSION A(N*(N+1)/2), L(N*(N+1)/2)
C
C     INTEGER I, II, IIM1, I1, J, K, M
C     DOUBLE PRECISION LII, LJ
C
      II = 0
      DO 50 I = 1, N
         I1 = II + 1
         II = II + I
         M = 1
         IF (I .EQ. 1) GO TO 30
         IIM1 = II - 1
         DO 20 J = I1, IIM1
              LJ = L(J)
              DO 10 K = I1, J
                   A(M) = A(M) + LJ*L(K)
                   M = M + 1
 10                CONTINUE
 20           CONTINUE
 30      LII = L(II)
         DO 40 J = I1, II
 40           A(J) = LII * L(J)
 50      CONTINUE
C
      RETURN
      END
      SUBROUTINE NL2ITR (D, IV, J, N, NN, P, R, V, X, NLUPRI)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  CARRY OUT NL2SOL (NONLINEAR LEAST-SQUARES) ITERATIONS  ***
C  ***  (NL2SOL VERSION 2.2)  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,NN,P,NLUPRI
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   D(P),J(NN,P),R(N),V(*),X(P)
      INTEGER
     +   IV(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   E,HALF,NEGONE,ONE,RDOF1,STTSST,T,T1,ZERO
      INTEGER
     +   CNVCOD,COSMIN,COVMAT,COVPRT,COVREQ,D0INIT,DGNORM,DIG,DIG1,
     +   DINIT,DSTNRM,DTYPE,DUMMY,F,F0,FDIF,FUZZ,G,G01,G1,GTSTEP,H,
     +   H0,H1,I,IERR,IM1,INCFAC,INITS,IPIV0,IPIV1,IPIVI,IPIVK,
     +   IPIVOT,IPK,IRC,JTINIT,JTOL1,K,KAGQT,KALM,KM1,L,LKY,LKY1,
     +   LMAT,LMAT1,LMAX0,LSTGST,M,MODE,MODEL,MXFCAL,MXITER,NFCALL,
     +   NFCOV,NFGCAL,NGCALL,NGCOV,NITER,NVSAVE,PHMXFC,PP1O2,
     +   PREDUC,QTR,QTR1,RAD0,RADFAC,RADINC,RADIUS,RD,RD0,RD1,RDK,
     +   RESTOR,RLIMIT,RSAVE,RSAVE1,S,S1,SIZE,SMH,SSTEP,STEP,STEP1,
     +   STGLIM,STLSTG,STPMOD,STPPAR,SUSED,SWITCH,TEMP1,TEMP2,
     +   TOOBIG,TUNER4,TUNER5,VSAVE1,W,W1,WSCALE,X0,X01,XIRC
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD,D1MACH,V2NORM
      LOGICAL
     +   STOPX
      EXTERNAL DOTPRD,D1MACH,V2NORM,STOPX
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ASSESS,COVCLC,DUPDAT,GQTSTP,ITSMRY,LMSTEP,PARCHK,QAPPLY,
     +   QRFACT,RPTMUL,SLUPDT,SLVMUL,VAXPY,VCOPY,VSCOPY
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,SQRT
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER IV(1), N, NN, P
C     DOUBLE PRECISION D(P), J(NN,P), R(N), V(1), X(P)
C     DIMENSION IV(60+P), V(93 + 2*N + P*(3*P+31)/2)
C
C
C--------------------------  PARAMETER USAGE  --------------------------
C
C D.... SCALE VECTOR.
C IV... INTEGER VALUE ARRAY.
C J.... N BY P JACOBIAN MATRIX (LEAD DIMENSION NN).
C N.... NUMBER OF OBSERVATIONS (COMPONENTS IN R).
C NN... LEAD DIMENSION OF J.
C P.... NUMBER OF PARAMETERS (COMPONENTS IN X).
C R.... RESIDUAL VECTOR.
C V.... FLOATING-POINT VALUE ARRAY.
C X.... PARAMETER VECTOR.
C
C  ***  DISCUSSION  ***
C
C        PARAMETERS IV, N, P, V, AND X ARE THE SAME AS THE CORRESPOND-
C     ING ONES TO NL2SOL (WHICH SEE), EXCEPT THAT V CAN BE SHORTER
C     (SINCE THE PART OF V THAT NL2SOL USES FOR STORING D, J, AND R IS
C     NOT NEEDED).  MOREOVER, COMPARED WITH NL2SOL, IV(1) MAY HAVE THE
C     TWO ADDITIONAL OUTPUT VALUES 1 AND 2, WHICH ARE EXPLAINED BELOW,
C     AS IS THE USE OF IV(TOOBIG) AND IV(NFGCAL).  THE VALUES IV(D),
C     IV(J), AND IV(R), WHICH ARE OUTPUT VALUES FROM NL2SOL (AND
C     NL2SNO), ARE NOT REFERENCED BY NL2ITR OR THE SUBROUTINES IT CALLS.
C        ON A FRESH START, I.E., A CALL ON NL2ITR WITH IV(1) = 0 OR 12,
C     NL2ITR ASSUMES THAT R = R(X), THE RESIDUAL AT X, AND J = J(X),
C     THE CORRESPONDING JACOBIAN MATRIX OF R AT X.
C
C IV(1) = 1 MEANS THE CALLER SHOULD SET R TO R(X), THE RESIDUAL AT X,
C             AND CALL NL2ITR AGAIN, HAVING CHANGED NONE OF THE OTHER
C             PARAMETERS.  AN EXCEPTION OCCURS IF R CANNOT BE EVALUATED
C             AT X (E.G. IF R WOULD OVERFLOW), WHICH MAY HAPPEN BECAUSE
C             OF AN OVERSIZED STEP.  IN THIS CASE THE CALLER SHOULD SET
C             IV(TOOBIG) = IV(2) TO 1, WHICH WILL CAUSE NL2ITR TO IG-
C             NORE R AND TRY A SMALLER STEP.  THE PARAMETER NF THAT
C             NL2SOL PASSES TO CALCR (FOR POSSIBLE USE BY CALCJ) IS A
C             COPY OF IV(NFCALL) = IV(6).
C IV(1) = 2 MEANS THE CALLER SHOULD SET J TO J(X), THE JACOBIAN MATRIX
C             OF R AT X, AND CALL NL2ITR AGAIN.  THE CALLER MAY CHANGE
C             D AT THIS TIME, BUT SHOULD NOT CHANGE ANY OF THE OTHER
C             PARAMETERS.  THE PARAMETER NF THAT NL2SOL PASSES TO
C             CALCJ IS IV(NFGCAL) = IV(7).  IF J CANNOT BE EVALUATED
C             AT X, THEN THE CALLER MAY SET IV(NFGCAL) TO 0, IN WHICH
C             CASE NL2ITR WILL RETURN WITH IV(1) = 15.
C
C  ***  GENERAL  ***
C
C     CODED BY DAVID M. GAY.
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C
C     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
C     MCS-7906671.
C        (SEE NL2SOL FOR REFERENCES.)
C
C+++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++++
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER DUMMY, DIG1, G1, G01, H0, H1, I, IM1, IPIVI, IPIVK, IPIV1,
C    1        IPK, K, KM1, L, LKY1, LMAT1, LSTGST, M, PP1O2, QTR1,
C    2        RDK, RD0, RD1, RSAVE1, SMH, SSTEP, STEP1, STPMOD, S1,
C    3        TEMP1, TEMP2, W1, X01
C     DOUBLE PRECISION E, RDOF1, STTSST, T, T1
C
C     ***  CONSTANTS  ***
C
C     DOUBLE PRECISION HALF, NEGONE, ONE, ZERO
C
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL ASSESS, COVCLC, DOTPRD, DUPDAT, GQTSTP, ITSMRY, LMSTEP,
C    1         PARCHK, QAPPLY, QRFACT, RPTMUL, SLUPDT, SLVMUL, STOPX,
C    2         VAXPY, VCOPY, VSCOPY, V2NORM
C     LOGICAL STOPX
C     DOUBLE PRECISION DOTPRD, D1MACH, V2NORM
C
C ASSESS... ASSESSES CANDIDATE STEP.
C COVCLC... COMPUTES COVARIANCE MATRIX.
C DOTPRD... RETURNS INNER PRODUCT OF TWO VECTORS.
C DUPDAT... UPDATES SCALE VECTOR D.
C GQTSTP... COMPUTES GOLDFELD-QUANDT-TROTTER STEP (AUGMENTED MODEL).
C ITSMRY... PRINTS ITERATION SUMMARY AND INFO ABOUT INITIAL AND FINAL X.
C LMSTEP... COMPUTES LEVENBERG-MARQUARDT STEP (GAUSS-NEWTON MODEL).
C PARCHK... CHECKS VALIDITY OF INPUT IV AND V VALUES.
C QAPPLY... APPLIES ORTHOGONAL MATRIX Q FROM QRFACT TO A VECTOR.
C QRFACT... COMPUTES QR DECOMPOSITION OF A MATRIX VIA HOUSEHOLDER TRANS.
C RPTMUL... MULTIPLIES VECTOR BY THE R MATRIX (AND/OR ITS TRANSPOSE)
C             STORED BY QRFACT.
C SLUPDT... PERFORMS QUASI-NEWTON UPDATE ON COMPACTLY STORED LOWER TRI-
C             ANGLE OF A SYMMETRIC MATRIX.
C STOPX.... RETURNS .TRUE. IF THE BREAK KEY HAS BEEN PRESSED.
C VAXPY.... COMPUTES SCALAR TIMES ONE VECTOR PLUS ANOTHER.
C VCOPY.... COPIES ONE VECTOR TO ANOTHER.
C VSCOPY... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR.
C V2NORM... RETURNS THE 2-NORM OF A VECTOR.
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER CNVCOD, COSMIN, COVMAT, COVPRT, COVREQ, DGNORM, DIG,
C    1        DINIT, DSTNRM, DTYPE, D0INIT, F, FDIF, FUZZ,
C    2        F0, G, GTSTEP, H, IERR, INCFAC, INITS, IPIVOT, IPIV0, IRC,
C    3        JTINIT, JTOL1, KAGQT, KALM, LKY, LMAT, LMAX0, MODE, MODEL,
C    4        MXFCAL, MXITER, NFCALL, NFGCAL, NFCOV, NGCOV, NGCALL,
C    5        NITER, NVSAVE, PHMXFC, PREDUC, QTR, RADFAC, RADINC,
C    6        RADIUS, RAD0, RD, RESTOR, RLIMIT, RSAVE, S, SIZE, STEP,
C    7        STGLIM, STLSTG, STPPAR, SUSED, SWITCH, TOOBIG, TUNER4,
C    8        TUNER5, VSAVE1, W, WSCALE, XIRC, X0
C
C  ***  IV SUBSCRIPT VALUES  ***
C
      DATA CNVCOD/34/, COVMAT/26/, COVPRT/14/,
     +     COVREQ/15/, DIG/43/, DTYPE/16/, G/28/, H/44/,
     +     IERR/32/, INITS/25/, IPIVOT/61/, IPIV0/60/,
     +     IRC/3/, KAGQT/35/, KALM/36/, LKY/37/, LMAT/58/,
     +     MODE/38/, MODEL/5/, MXFCAL/17/, MXITER/18/,
     +     NFCALL/6/, NFGCAL/7/, NFCOV/40/, NGCOV/41/,
     +     NGCALL/30/, NITER/31/, QTR/49/,
     +     RADINC/8/, RD/51/, RESTOR/9/, RSAVE/52/, S/53/,
     +     STEP/55/, STGLIM/11/, STLSTG/56/, SUSED/57/,
     +     SWITCH/12/, TOOBIG/2/, W/59/, XIRC/13/, X0/60/
C
C  ***  V SUBSCRIPT VALUES  ***
C
      DATA COSMIN/43/, DGNORM/1/, DINIT/38/, DSTNRM/2/,
     +     D0INIT/37/, F/10/, FDIF/11/, FUZZ/45/,
     +     F0/13/, GTSTEP/4/, INCFAC/23/,
     +     JTINIT/39/, JTOL1/87/, LMAX0/35/,
     +     NVSAVE/9/, PHMXFC/21/, PREDUC/7/,
     +     RADFAC/16/, RADIUS/8/, RAD0/9/, RLIMIT/42/,
     +     SIZE/47/, STPPAR/5/, TUNER4/29/, TUNER5/30/,
     +     VSAVE1/78/, WSCALE/48/
C
C
      DATA HALF/0.5D+00/, NEGONE/-1.0D+00/, ONE/1.0D+00/, ZERO/0.0D+00/
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      I = IV(1)
      IF (I .EQ. 1) GO TO 20
      IF (I .EQ. 2) GO TO 50
C
C  ***  CHECK VALIDITY OF IV AND V INPUT VALUES  ***
C
C     ***  NOTE -- IF IV(1) = 0, THEN PARCHK CALLS DFAULT(IV, V)  ***
      CALL PARCHK(IV, N, NN, P, V)
      I = IV(1) - 2
      IF (I .GT. 10) GO TO 999
      GO TO (350, 350, 350, 350, 350, 350, 195, 160, 195, 10), I
C
C  ***  INITIALIZATION AND STORAGE ALLOCATION  ***
C
 10   IV(NITER) = 0
      IV(NFCALL) = 1
      IV(NGCALL) = 1
      IV(NFGCAL) = 1
      IV(MODE) = -1
      IV(STGLIM) = 2
      IV(TOOBIG) = 0
      IV(CNVCOD) = 0
      IV(COVMAT) = 0
      IV(NFCOV) = 0
      IV(NGCOV) = 0
      IV(KALM) = -1
      IV(RADINC) = 0
      IV(S) = JTOL1 + 2*P
      PP1O2 = P * (P + 1) / 2
      IV(X0) = IV(S) + PP1O2
      IV(STEP) = IV(X0) + P
      IV(STLSTG) = IV(STEP) + P
      IV(DIG) = IV(STLSTG) + P
      IV(G) = IV(DIG) + P
      IV(LKY) = IV(G) + P
      IV(RD) = IV(LKY) + P
      IV(RSAVE) = IV(RD) + P
      IV(QTR) = IV(RSAVE) + N
      IV(H) = IV(QTR) + N
      IV(W) = IV(H) + PP1O2
      IV(LMAT) = IV(W) + 4*P + 7
C     +++ LENGTH OF W = P*(P+9)/2 + 7.  LMAT IS CONTAINED IN W.
      IF (V(DINIT) .GE. ZERO) CALL VSCOPY(P, D, V(DINIT))
      IF (V(JTINIT) .GT. ZERO) CALL VSCOPY(P, V(JTOL1), V(JTINIT))
      I = JTOL1 + P
      IF (V(D0INIT) .GT. ZERO) CALL VSCOPY(P, V(I), V(D0INIT))
      V(RAD0) = ZERO
      V(STPPAR) = ZERO
      V(RADIUS) = V(LMAX0) / (ONE + V(PHMXFC))
C
C  ***  SET INITIAL MODEL AND S MATRIX  ***
C
      IV(MODEL) = 1
      IF (IV(INITS) .EQ. 2) IV(MODEL) = 2
      S1 = IV(S)
      IF (IV(INITS) .EQ. 0) CALL VSCOPY(PP1O2, V(S1), ZERO)
C
C  ***  COMPUTE FUNCTION VALUE (HALF THE SUM OF SQUARES)  ***
C
 20   T = V2NORM(N, R)
      IF (T .GT. V(RLIMIT)) IV(TOOBIG) = 1
      IF (IV(TOOBIG) .NE. 0) GO TO 30
      V(F) = 0.0
      IF (T.GT.SQRT(D1MACH(1))) V(F) = HALF * T**2
 30   IF (IV(MODE)) 40, 350, 730
C
 40   IF (IV(TOOBIG) .EQ. 0) GO TO 60
         IV(1) = 13
         GO TO 900
C
C  ***  MAKE SURE JACOBIAN COULD BE COMPUTED  ***
C
 50   IF (IV(NFGCAL) .NE. 0) GO TO 60
         IV(1) = 15
         GO TO 900
C
C  ***  COMPUTE GRADIENT  ***
C
 60   IV(KALM) = -1
      G1 = IV(G)
      DO 70 I = 1, P
         V(G1) = DOTPRD(N, R, J(1,I))
         G1 = G1 + 1
 70      CONTINUE
      IF (IV(MODE) .GT. 0) GO TO 710
C
C  ***  UPDATE D AND MAKE COPIES OF R FOR POSSIBLE USE LATER  ***
C
      IF (IV(DTYPE) .GT. 0) CALL DUPDAT(D, IV, J, N, NN, P, V)
      RSAVE1 = IV(RSAVE)
      CALL VCOPY(N, V(RSAVE1), R)
      QTR1 = IV(QTR)
      CALL VCOPY(N, V(QTR1), R)
C
C  ***  COMPUTE  D**-1 * GRADIENT  ***
C
      G1 = IV(G)
      DIG1 = IV(DIG)
      K = DIG1
      DO 80 I = 1, P
         V(K) = V(G1) / D(I)
         K = K + 1
         G1 = G1 + 1
 80      CONTINUE
      V(DGNORM) = V2NORM(P, V(DIG1))
C
      IF (IV(CNVCOD) .NE. 0) GO TO 700
      IF (IV(MODE) .EQ. 0) GO TO 570
      IV(MODE) = 0
C
C
C-----------------------------  MAIN LOOP  -----------------------------
C
C
C  ***  PRINT ITERATION SUMMARY, CHECK ITERATION LIMIT  ***
C
 150  CALL ITSMRY(D, IV, P, V, X, NLUPRI)
 160  K = IV(NITER)
      IF (K .LT. IV(MXITER)) GO TO 170
         IV(1) = 10
         GO TO 900
 170  IV(NITER) = K + 1
C
C  ***  UPDATE RADIUS  ***
C
      IF (K .EQ. 0) GO TO 185
      STEP1 = IV(STEP)
      DO 180 I = 1, P
         V(STEP1) = D(I) * V(STEP1)
         STEP1 = STEP1 + 1
 180     CONTINUE
      STEP1 = IV(STEP)
      V(RADIUS) = V(RADFAC) * V2NORM(P, V(STEP1))
C
C  ***  INITIALIZE FOR START OF NEXT ITERATION  ***
C
 185  X01 = IV(X0)
      V(F0) = V(F)
      IV(KAGQT) = -1
      IV(IRC) = 4
      IV(H) = -ABS(IV(H))
      IV(SUSED) = IV(MODEL)
C
C     ***  COPY X TO X0  ***
C
      CALL VCOPY(P, V(X01), X)
C
C  ***  CHECK STOPX AND FUNCTION EVALUATION LIMIT  ***
C
 190  IF (.NOT. STOPX(DUMMY)) GO TO 200
         IV(1) = 11
         GO TO 205
C
C     ***  COME HERE WHEN RESTARTING AFTER FUNC. EVAL. LIMIT OR STOPX.
C
 195  IF (V(F) .GE. V(F0)) GO TO 200
         V(RADFAC) = ONE
         K = IV(NITER)
         GO TO 170
C
 200  IF (IV(NFCALL) .LT. IV(MXFCAL) + IV(NFCOV)) GO TO 210
         IV(1) = 9
 205     IF (V(F) .GE. V(F0)) GO TO 900
C
C        ***  IN CASE OF STOPX OR FUNCTION EVALUATION LIMIT WITH
C        ***  IMPROVED V(F), EVALUATE THE GRADIENT AT X.
C
              IV(CNVCOD) = IV(1)
              GO TO 560
C
C. . . . . . . . . . . . .  COMPUTE CANDIDATE STEP  . . . . . . . . . .
C
 210  STEP1 = IV(STEP)
      W1 = IV(W)
      IF (IV(MODEL) .EQ. 2) GO TO 240
C
C  ***  COMPUTE LEVENBERG-MARQUARDT STEP  ***
C
         QTR1 = IV(QTR)
         IF (IV(KALM) .GE. 0) GO TO 215
              RD1 = IV(RD)
              IF (-1 .EQ. IV(KALM)) CALL QRFACT(NN, N, P, J, V(RD1),
     +                                   IV(IPIVOT), IV(IERR), 0, V(W1))
              CALL QAPPLY(NN, N, P, J, V(QTR1), IV(IERR))
 215     H1 = IV(H)
         IF (H1 .GT. 0) GO TO 230
C
C        ***  COPY R MATRIX TO H  ***
C
              H1 = -H1
              IV(H) = H1
              K = H1
              RD1 = IV(RD)
              V(K) = V(RD1)
              IF (P .EQ. 1) GO TO 230
              DO 220 I = 2, P
                   CALL VCOPY(I-1, V(K+1), J(1,I))
                   K = K + I
                   RD1 = RD1 + 1
                   V(K) = V(RD1)
 220               CONTINUE
C
 230     G1 = IV(G)
         CALL LMSTEP(D, V(G1), IV(IERR), IV(IPIVOT), IV(KALM), P,
     +               V(QTR1), V(H1), V(STEP1), V, V(W1))
         GO TO 310
C
C  ***  COMPUTE GOLDFELD-QUANDT-TROTTER STEP (AUGMENTED MODEL)  ***
C
 240  IF (IV(H) .GT. 0) GO TO 300
C
C     ***  SET H TO  D**-1 * ( (J**T)*J + S) ) * D**-1.  ***
C
         H1 = -IV(H)
         IV(H) = H1
         S1 = IV(S)
         IF (IV(KALM) .GE. 0) GO TO 270
C
C        ***  J IS IN ITS ORIGINAL FORM  ***
C
              DO 260 I = 1, P
                   T = ONE / D(I)
                   DO 250 K = 1, I
                        V(H1) = T*(DOTPRD(N,J(1,I),J(1,K))+V(S1)) / D(K)
                        H1 = H1 + 1
                        S1 = S1 + 1
 250                    CONTINUE
 260               CONTINUE
              GO TO 300
C
C  ***  LMSTEP HAS APPLIED QRFACT TO J  ***
C
 270     SMH = S1 - H1
         H0 = H1 - 1
         IPIV1 = IV(IPIVOT)
         T1 = ONE / D(IPIV1)
         RD0 = IV(RD) - 1
         RDOF1 = V(RD0 + 1)
         DO 290 I = 1, P
              L = IPIV0 + I
              IPIVI = IV(L)
              H1 = H0 + IPIVI*(IPIVI-1)/2
              L = H1 + IPIVI
              M = L + SMH
C             ***  V(L) = H(IPIVOT(I), IPIVOT(I))  ***
C             ***  V(M) = S(IPIVOT(I), IPIVOT(I))  ***
              T = ONE / D(IPIVI)
              RDK = RD0 + I
              E = V(RDK)**2
              IF (I .GT. 1) E = E + DOTPRD(I-1, J(1,I), J(1,I))
              V(L) = (E + V(M)) * T**2
              IF (I .EQ. 1) GO TO 290
              L = H1 + IPIV1
              IF (IPIVI .LT. IPIV1) L = L +
     +                               ((IPIV1-IPIVI)*(IPIV1+IPIVI-3))/2
              M = L + SMH
C             ***  V(L) = H(IPIVOT(I), IPIVOT(1))  ***
C             ***  V(M) = S(IPIVOT(I), IPIVOT(1))  ***
              V(L) = T * (RDOF1 * J(1,I)  +  V(M)) * T1
              IF (I .EQ. 2) GO TO 290
              IM1 = I - 1
              DO 280 K = 2, IM1
                   IPK = IPIV0 + K
                   IPIVK = IV(IPK)
                   L = H1 + IPIVK
                   IF (IPIVI .LT. IPIVK) L = L +
     +                               ((IPIVK-IPIVI)*(IPIVK+IPIVI-3))/2
                   M = L + SMH
C                  ***  V(L) = H(IPIVOT(I), IPIVOT(K))  ***
C                  ***  V(M) = S(IPIVOT(I), IPIVOT(K))  ***
                   KM1 = K - 1
                   RDK = RD0 + K
                   V(L) = T * (DOTPRD(KM1, J(1,I), J(1,K)) +
     +                            V(RDK)*J(K,I) + V(M)) / D(IPIVK)
 280               CONTINUE
 290          CONTINUE
C
C  ***  COMPUTE ACTUAL GOLDFELD-QUANDT-TROTTER STEP  ***
C
 300  H1 = IV(H)
      DIG1 = IV(DIG)
      LMAT1 = IV(LMAT)
      CALL GQTSTP(D, V(DIG1), V(H1), IV(KAGQT), V(LMAT1), P, V(STEP1),
     +            V, V(W1))
C
C
C  ***  COMPUTE R(X0 + STEP)  ***
C
 310  IF (IV(IRC) .EQ. 6) GO TO 350
      X01 = IV(X0)
      STEP1 = IV(STEP)
      CALL VAXPY(P, X, ONE, V(STEP1), V(X01))
      IV(NFCALL) = IV(NFCALL) + 1
      IV(1) = 1
      IV(TOOBIG) = 0
      GO TO 999
C
C. . . . . . . . . . . . .  ASSESS CANDIDATE STEP  . . . . . . . . . . .
C
 350  STEP1 = IV(STEP)
      LSTGST = IV(STLSTG)
      X01 = IV(X0)
      CALL ASSESS(D, IV, P, V(STEP1), V(LSTGST), V, X, V(X01))
C
C  ***  IF NECESSARY, SWITCH MODELS AND/OR RESTORE R  ***
C
      IF (IV(SWITCH) .EQ. 0) GO TO 360
         IV(H) = -ABS(IV(H))
         IV(SUSED) = IV(SUSED) + 2
         CALL VCOPY(NVSAVE, V, V(VSAVE1))
 360  IF (IV(RESTOR) .EQ. 0) GO TO 390
         RSAVE1 = IV(RSAVE)
         CALL VCOPY(N, R, V(RSAVE1))
 390  L = IV(IRC) - 4
      STPMOD = IV(MODEL)
      IF (L .GT. 0) GO TO (410,440,450,450,450,450,450,450,640,570), L
C
C  ***  DECIDE WHETHER TO CHANGE MODELS  ***
C
      E = V(PREDUC) - V(FDIF)
      SSTEP = IV(LKY)
      S1 = IV(S)
      CALL SLVMUL(P, V(SSTEP), V(S1), V(STEP1))
      STTSST = HALF * DOTPRD(P, V(STEP1), V(SSTEP))
      IF (IV(MODEL) .EQ. 1) STTSST = -STTSST
      IF (ABS(E + STTSST) * V(FUZZ) .GE. ABS(E)) GO TO 400
C
C     ***  SWITCH MODELS  ***
C
         IV(MODEL) = 3 - IV(MODEL)
         IF (IV(MODEL) .EQ. 1) IV(KAGQT) = -1
         IF (IV(MODEL) .EQ. 2 .AND. IV(KALM) .GT. 0) IV(KALM) = 0
         IF (-2 .LT. L) GO TO 480
              IV(H) = -ABS(IV(H))
              IV(SUSED) = IV(SUSED) + 2
              CALL VCOPY(NVSAVE, V(VSAVE1), V)
              GO TO 420
C
 400  IF (-3 .LT. L) GO TO 480
C
C     ***  RECOMPUTE STEP WITH DECREASED RADIUS  ***
C
         V(RADIUS) = V(RADFAC) * V(DSTNRM)
         GO TO 190
C
C  ***  RECOMPUTE STEP, SAVING V VALUES AND R IF NECESSARY  ***
C
 410  V(RADIUS) = V(RADFAC) * V(DSTNRM)
 420  IF (V(F) .GE. V(F0)) GO TO 190
      RSAVE1 = IV(RSAVE)
      CALL VCOPY(N, V(RSAVE1), R)
      GO TO 190
C
C  ***  COMPUTE STEP OF LENGTH V(LMAX0) FOR SINGULAR CONVERGENCE TEST
C
 440  V(RADIUS) = V(LMAX0)
      GO TO 210
C
C  ***  CONVERGENCE OR FALSE CONVERGENCE  ***
C
 450  IV(CNVCOD) = L
      IF (V(F) .GE. V(F0)) GO TO 700
         IF (IV(XIRC) .EQ. 14) GO TO 700
              IV(XIRC) = 14
C
C. . . . . . . . . . . .  PROCESS ACCEPTABLE STEP  . . . . . . . . . . .
C
 480  IV(COVMAT) = 0
C
C  ***  SET  LKY = (J(X0)**T) * R(X)  ***
C
      LKY1 = IV(LKY)
      IF (IV(KALM) .GE. 0) GO TO 500
C
C     ***  JACOBIAN HAS NOT BEEN MODIFIED  ***
C
         DO 490 I = 1, P
              V(LKY1) = DOTPRD(N, J(1,I), R)
              LKY1 = LKY1 + 1
 490          CONTINUE
         GO TO 510
C
C  ***  QRFACT HAS BEEN APPLIED TO J.  STORE COPY OF R IN QTR AND  ***
C  ***  APPLY Q TO IT.                                             ***
C
 500  QTR1 = IV(QTR)
      CALL VCOPY(N, V(QTR1), R)
      CALL QAPPLY(NN, N, P, J, V(QTR1), IV(IERR))
C
C  ***  MULTIPLY TOP P-VECTOR IN QTR BY PERMUTED UPPER TRIANGLE    ***
C  ***  STORED BY QRFACT IN J AND RD.                              ***
C
      RD1 = IV(RD)
      TEMP1 = IV(STLSTG)
      CALL RPTMUL(3, IV(IPIVOT), J, NN, P, V(RD1), V(QTR1), V(LKY1),
     +            V(TEMP1))
C
C  ***  SEE WHETHER TO SET V(RADFAC) BY GRADIENT TESTS  ***
C
 510  IF (IV(IRC) .NE. 3) GO TO 560
         STEP1 = IV(STEP)
         TEMP1 = IV(STLSTG)
         TEMP2 = IV(X0)
C
C     ***  SET  TEMP1 = HESSIAN * STEP  FOR USE IN GRADIENT TESTS  ***
C
         IF (STPMOD .EQ. 2) GO TO 530
C
C        ***  STEP COMPUTED USING GAUSS-NEWTON MODEL  ***
C        ***  -- QRFACT HAS BEEN APPLIED TO J         ***
C
              RD1 = IV(RD)
              CALL RPTMUL(2, IV(IPIVOT), J, NN, P, V(RD1),
     +                    V(STEP1), V(TEMP1), V(TEMP2))
              GO TO 560
C
C     ***  STEP COMPUTED USING AUGMENTED MODEL  ***
C
 530     H1 = IV(H)
         K = TEMP2
         DO 540 I = 1, P
              V(K) = D(I) * V(STEP1)
              K = K + 1
              STEP1 = STEP1 + 1
 540          CONTINUE
         CALL SLVMUL(P, V(TEMP1), V(H1), V(TEMP2))
         DO 550 I = 1, P
              V(TEMP1) = D(I) * V(TEMP1)
              TEMP1 = TEMP1 + 1
 550          CONTINUE
C
C  ***  SAVE OLD GRADIENT AND COMPUTE NEW ONE  ***
C
 560  IV(NGCALL) = IV(NGCALL) + 1
      G1 = IV(G)
      G01 = IV(W)
      CALL VCOPY(P, V(G01), V(G1))
      IV(1) = 2
      GO TO 999
C
C  ***  INITIALIZATIONS -- G0 = G - G0, ETC.  ***
C
 570  G01 = IV(W)
      G1 = IV(G)
      CALL VAXPY(P, V(G01), NEGONE, V(G01), V(G1))
      STEP1 = IV(STEP)
      TEMP1 = IV(STLSTG)
      TEMP2 = IV(X0)
      IF (IV(IRC) .NE. 3) GO TO 600
C
C  ***  SET V(RADFAC) BY GRADIENT TESTS  ***
C
C     ***  SET  TEMP1 = D**-1 * (HESSIAN * STEP  +  (G(X0) - G(X)))  ***
C
         K = TEMP1
         L = G01
         DO 580 I = 1, P
              V(K) = (V(K) - V(L)) / D(I)
              K = K + 1
              L = L + 1
 580          CONTINUE
C
C        ***  DO GRADIENT TESTS  ***
C
         IF (V2NORM(P, V(TEMP1)) .LE. V(DGNORM) * V(TUNER4))  GO TO 590
              IF (DOTPRD(P, V(G1), V(STEP1))
     +                  .GE. V(GTSTEP) * V(TUNER5))  GO TO 600
 590               V(RADFAC) = V(INCFAC)
C
C  ***  FINISH COMPUTING LKY = ((J(X) - J(X0))**T) * R  ***
C
C     ***  CURRENTLY LKY = (J(X0)**T) * R  ***
C
 600  LKY1 = IV(LKY)
      CALL VAXPY(P, V(LKY1), NEGONE, V(LKY1), V(G1))
C
C  ***  DETERMINE SIZING FACTOR V(SIZE)  ***
C
C     ***  SET TEMP1 = S * STEP  ***
      S1 = IV(S)
      CALL SLVMUL(P, V(TEMP1), V(S1), V(STEP1))
C
      T1 = ABS(DOTPRD(P, V(STEP1), V(TEMP1)))
      T = ABS(DOTPRD(P, V(STEP1), V(LKY1)))
      V(SIZE) = ONE
      IF (T .LT. T1) V(SIZE) = T / T1
C
C  ***  UPDATE S  ***
C
      CALL SLUPDT(V(S1), V(COSMIN), P, V(SIZE), V(STEP1), V(TEMP1),
     +            V(TEMP2), V(G01), V(WSCALE), V(LKY1))
      IV(1) = 2
      GO TO 150
C
C. . . . . . . . . . . . . .  MISC. DETAILS  . . . . . . . . . . . . . .
C
C  ***  BAD PARAMETERS TO ASSESS  ***
C
 640  IV(1) = 14
      GO TO 900
C
C  ***  CONVERGENCE OBTAINED -- COMPUTE COVARIANCE MATRIX IF DESIRED ***
C
 700  IF (IV(COVREQ) .EQ. 0 .AND. IV(COVPRT) .EQ. 0) GO TO 760
      IF (IV(COVMAT) .NE. 0) GO TO 760
      IF (IV(CNVCOD) .GE. 7) GO TO 760
      IV(MODE) = 0
 710  CALL COVCLC(I, D, IV, J, N, NN, P, R, V, X)
      GO TO (720, 720, 740, 750), I
 720  IV(NFCOV) = IV(NFCOV) + 1
      IV(NFCALL) = IV(NFCALL) + 1
      IV(RESTOR) = I
      IV(1) = 1
      GO TO 999
C
 730  IF (IV(RESTOR) .EQ. 1 .OR. IV(TOOBIG) .NE. 0) GO TO 710
      IV(NFGCAL) = IV(NFCALL)
 740  IV(NGCOV) = IV(NGCOV) + 1
      IV(NGCALL) = IV(NGCALL) + 1
      IV(1) = 2
      GO TO 999
C
 750  IV(MODE) = 0
      IF (IV(NITER) .EQ. 0) IV(MODE) = -1
C
 760  IV(1) = IV(CNVCOD)
      IV(CNVCOD) = 0
C
C  ***  PRINT SUMMARY OF FINAL ITERATION AND OTHER REQUESTED ITEMS  ***
C
 900  CALL ITSMRY(D, IV, P, V, X, NLUPRI)
C
 999  RETURN
      END
      SUBROUTINE NL2SNO(N, P, X, XSTP, YSTP, CALCR, IV, V,
     +                  UIPARM, URPARM, UFPARM, NLUPRI)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  LIKE NL2SOL, BUT WITHOUT CALCJ -- MINIMIZE NONLINEAR SUM OF  ***
C  ***  SQUARES USING FINITE DIFFERENCE JACOBIAN APPROXIMATIONS      ***
C  ***  (NL2SOL VERSION 2.2)  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,P,NLUPRI
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   URPARM(*),V(*),X(P),XSTP(N),YSTP(N)
      INTEGER
     +   IV(*),UIPARM(*)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL CALCR,UFPARM
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   H,HFAC,HLIM,NEGPT5,ONE,XK,ZERO
      INTEGER
     +   COVPRT,COVREQ,D,D1,DK,DLTFDJ,DTYPE,I,J,J1,J1K,K,NF,NFCALL,
     +   NFGCAL,R,R1,RN,TOOBIG
      LOGICAL
     +   STRTED
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   RMDCON
      EXTERNAL RMDCON
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DFAULT,ITSMRY,NL2ITR,VSCOPY
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX
C
C     INTEGER N, P, IV(1), UIPARM(1)
C     DOUBLE PRECISION X(P), V(1), URPARM(1)
C     DIMENSION IV(60+P),  V(93 + N*P + 3*N + P*(3*P+33)/2)
C     EXTERNAL CALCR, UFPARM
C
C-----------------------------  DISCUSSION  ----------------------------
C
C        THE PARAMETERS FOR NL2SNO ARE THE SAME AS THOSE FOR NL2SOL
C     (WHICH SEE), EXCEPT THAT CALCJ IS OMITTED.  INSTEAD OF CALLING
C     CALCJ TO OBTAIN THE JACOBIAN MATRIX OF R AT X, NL2SNO COMPUTES
C     AN APPROXIMATION TO IT BY FINITE (FORWARD) DIFFERENCES -- SEE
C     V(DLTFDJ) BELOW.  NL2SNO USES FUNCTION VALUES ONLY WHEN COMPUT-
C     THE COVARIANCE MATRIX (RATHER THAN THE FUNCTIONS AND GRADIENTS
C     THAT NL2SOL MAY USE).  TO DO SO, NL2SNO SETS IV(COVREQ) TO -1 IF
C     IV(COVPRT) = 1 WITH IV(COVREQ) = 0 AND TO MINUS ITS ABSOLUTE
C     VALUE OTHERWISE.  THUS V(DELTA0) IS NEVER REFERENCED AND ONLY
C     V(DLTFDC) MATTERS -- SEE NL2SOL FOR A DESCRIPTION OF V(DLTFDC).
C        THE NUMBER OF EXTRA CALLS ON CALCR USED IN COMPUTING THE JACO-
C     BIAN APPROXIMATION ARE NOT INCLUDED IN THE FUNCTION EVALUATION
C     COUNT IV(NFCALL) AND ARE NOT OTHERWISE REPORTED.
C
C V(DLTFDJ)... V(36) HELPS CHOOSE THE STEP SIZE USED WHEN COMPUTING THE
C             FINITE DIFFERENCE JACOBIAN MATRIX.  FOR DIFFERENCES IN-
C             VOLVING X(I), THE STEP SIZE FIRST TRIED IS
C                       V(DLTFDJ) * MAX(ABS(X(I)), 1/D(I)),
C             WHERE D IS THE CURRENT SCALE VECTOR (SEE REF. 1).  (IF
C             THIS STEP IS TOO BIG, I.E., IF CALCR SETS NF TO 0, THEN
C             SMALLER STEPS ARE TRIED UNTIL THE STEP SIZE IS SHRUNK 
C             BELOW 1000 * MACHEP, WHERE MACHEP IS THE UNIT ROUNDOFF.
C             DEFAULT = MACHEP**0.5.
C
C  ***  REFERENCES  ***
C
C 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1980), AN ADAPTIVE
C             NONLINEAR LEAST-SQUARES ALGORITHM, SUBMITTED TO ACM TRANS.
C             MATH. SOFTWARE.
C
C  ***  GENERAL  ***
C
C     CODED BY DAVID M. GAY.
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
C     MCS-7906671.
C
C+++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++++
C
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL DFAULT, ITSMRY, NL2ITR, RMDCON, VSCOPY
C     DOUBLE PRECISION RMDCON
C
C DFAULT... SUPPLIES DEFAULT PARAMETER VALUES.
C ITSMRY... PRINTS ITERATION SUMMARY AND INFO ABOUT INITIAL AND FINAL X.
C NL2ITR... REVERSE COMMUNICATION ROUTINE THAT CARRIES OUT NL2SOL ALGO-
C             RITHM.
C RMDCON... RETURNS MACHINE DEPENDENT CONSTANTS.
C VSCOPY... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR.
C
C     LOGICAL STRTED
C     INTEGER DK, D1, I, J1, J1K, K, NF, RN, R1
C     DOUBLE PRECISION H, HFAC, HLIM, NEGPT5, ONE, XK, ZERO
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER COVPRT, COVREQ, D, DLTFDJ, DTYPE, J, NFCALL, NFGCAL, R,
C    1        TOOBIG
C
      DATA HFAC/1.0D+03/, HLIM/0.0D+00/, NEGPT5/-0.5D+00/,
     +     ONE/1.0D+00/, ZERO/0.0D+00/
C
C  ***  IV SUBSCRIPT VALUES  ***
C
      DATA COVPRT/14/, COVREQ/15/, D/27/, DTYPE/16/, J/33/,
     +     NFCALL/6/, NFGCAL/7/, R/50/, TOOBIG/2/
C
C  ***  V SUBSCRIPT VALUES  ***
C
      DATA DLTFDJ/36/

      D1 = 94 + 2*N + P*(3*P + 31)/2
      IV(D) = D1
      R1 = D1 + P
      IV(R) = R1
      J1 = R1 + N
      IV(J) = J1
      RN = J1 - 1
      IF (IV(1) .EQ. 0) CALL DFAULT(IV, V)
      IV(COVREQ) = -ABS(IV(COVREQ))
      IF (IV(COVPRT) .NE. 0 .AND. IV(COVREQ) .EQ. 0) IV(COVREQ) = -1
      STRTED = .TRUE.
      IF (IV(1) .NE. 12) GO TO 80
         STRTED = .FALSE.
         IV(NFCALL) = 1
         IV(NFGCAL) = 1
C        ***  INITIALIZE SCALE VECTOR D TO ONES FOR COMPUTING
C        ***  INITIAL JACOBIAN.
         IF (IV(DTYPE) .GT. 0) CALL VSCOPY(P, V(D1), ONE)
C
 10   NF = IV(NFCALL)
CPFP 1
C SUBROUTINE CALCR(m, n, xstp, ystp, x, iflag, fvec, NLUPRI)
        CALL CALCR(N, P, XSTP, YSTP, X,  NF,  V(R1), NLUPRI)
      IF (STRTED) GO TO 20
         IF (NF .GT. 0) GO TO 30
              IV(1) = 13
              GO TO 90
C
 20   IF (NF .LE. 0) IV(TOOBIG) = 1
      GO TO 80
C
C  ***  COMPUTE FINITE DIFFERENCE JACOBIAN  ***
C
30    J1K = J1
      DK = D1
      DO 70 K = 1, P
         XK = X(K)
         H = V(DLTFDJ) * MAX(ABS(XK), ONE/V(DK))
         DK = DK + 1
 40      X(K) = XK + H
         NF = IV(NFGCAL)
CPFP 2
C SUBROUTINE CALCR(m, n, xstp, ystp, x, iflag, fvec, NLUPRI)
        CALL CALCR(N, P, XSTP, YSTP, X,  NF, V(J1K), NLUPRI)
         IF (NF .GT. 0) GO TO 50
              IF (HLIM .EQ. ZERO) HLIM = HFAC * RMDCON(3)
C             ***  HLIM = HFAC TIMES THE UNIT ROUNDOFF  ***
              H = NEGPT5 * H
              IF (ABS(H) .GE. HLIM) GO TO 40
                   IV(1) = 15
                   GO TO 90
 50      X(K) = XK
         DO 60 I = R1, RN
              V(J1K) = (V(J1K) - V(I)) / H
              J1K = J1K + 1
 60           CONTINUE
 70      CONTINUE
C
      STRTED = .TRUE.
C
 80   CALL NL2ITR(V(D1), IV, V(J1), N, N, P, V(R1), V, X, NLUPRI)
      IF (IV(1) - 2) 10, 30, 999
C
 90   CALL ITSMRY(V(D1), IV, P, V, X, NLUPRI)
C
 999  RETURN
      END
C
c           call NL2SOL(N, NA, ACOEF, AX, AY, CALCR, CALCJ, IV, V,
c     &                 UIPARM, URPARM, UFPARM, NPRINT)
C
      SUBROUTINE NL2SOL(N, P, X, XSTP, YSTP, CALCR, CALCJ, IV, V,
     +                  UIPARM, URPARM, UFPARM, NLUPRI)
c
c*********************************************************************72
c
C  ***  MINIMIZE NONLINEAR SUM OF SQUARES USING ANALYTIC JACOBIAN  ***
C  ***  (NL2SOL VERSION 2.2)  ***
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,P, NLUPRI
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   URPARM(*),V(*),X(P),XSTP(N),YSTP(N)
      INTEGER
     +   IV(*),UIPARM(*)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL CALCJ,CALCR,UFPARM
C
C  LOCAL SCALARS
      INTEGER
     +   D,D1,J,J1,NF,NFCALL,NFGCAL,R,R1,TOOBIG
      LOGICAL
     +   STRTED
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ITSMRY,NL2ITR
C
C
C     INTEGER N, P, IV(1), UIPARM(1)
C     DOUBLE PRECISION X(P), V(1), URPARM(1)
C     DIMENSION IV(60+P),  V(93 + N*P + 3*N + P*(3*P+33)/2)
C     DIMENSION UIPARM(*), URPARM(*)
C     EXTERNAL CALCR, CALCJ, UFPARM
C
C  ***  PURPOSE  ***
C
C        GIVEN A P-VECTOR X OF PARAMETERS, CALCR COMPUTES AN N-VECTOR
C     R = R(X) OF RESIDUALS CORRESPONDING TO X.  (R(X) PROBABLY ARISES
C     FROM A NONLINEAR MODEL INVOLVING P PARAMETERS AND N OBSERVATIONS.)
C     THIS ROUTINE INTERACTS WITH NL2ITR TO SEEK A PARAMETER VECTOR X
C     THAT MINIMIZES THE SUM OF THE SQUARES OF (THE COMPONENTS OF) R(X),
C     I.E., THAT MINIMIZES THE SUM-OF-SQUARES FUNCTION
C     F(X) = (R(X)**T) * R(X) / 2.  R(X) IS ASSUMED TO BE A TWICE CON-
C     TINUOUSLY DIFFERENTIABLE FUNCTION OF X.
C
C PARAMETER USAGE
C
C N........ (INPUT) THE NUMBER OF OBSERVATIONS, I.E., THE NUMBER OF
C                  COMPONENTS IN R(X).  N MUST BE .GE. P.
C P........ (INPUT) THE NUMBER OF PARAMETERS (COMPONENTS IN X).  P MUST
C                  BE POSITIVE.
C X........ (INPUT/OUTPUT).  ON INPUT, X IS AN INITIAL GUESS AT THE
C                  DESIRED PARAMETER ESTIMATE.  ON OUTPUT, X CONTAINS
C                  THE BEST PARAMETER ESTIMATE FOUND.
C CALCR.... (INPUT) A SUBROUTINE WHICH, GIVEN X, COMPUTES R(X).  CALCR
C                  MUST BE DECLARED EXTERNAL IN THE CALLING PROGRAM.
C                  IT IS INVOKED BY
C                       CALL CALCR(N,P,X,NF,R,UIPARM,URPARM,UFPARM)
C                  WHEN CALCR IS CALLED, NF IS THE INVOCATION COUNT
C                  FOR CALCR.  IT IS INCLUDED FOR POSSIBLE USE WITH
C                  CALCJ.  IF X IS OUT OF BOUNDS (E.G. IF IT WOULD
C                  CAUSE OVERFLOW IN COMPUTING R(X)), THEN CALCR SHOULD
C                  SET NF TO 0.  THIS WILL CAUSE A SHORTER STEP TO BE
C                  ATTEMPTED.  THE OTHER PARAMETERS ARE AS DESCRIBED
C                  ABOVE AND BELOW.  CALCR SHOULD NOT CHANGE N, P, OR X.
C CALCJ.... (INPUT) A SUBROUTINE WHICH, GIVEN X, COMPUTES THE JACOBIAN
C                  MATRIX J OF R AT X, I.E., THE N BY P MATRIX WHOSE
C                  (I,K) ENTRY IS THE PARTIAL DERIVATIVE OF THE I-TH
C                  COMPONENT OF R WITH RESPECT TO X(K).  CALCJ MUST BE
C                  DECLARED EXTERNAL IN THE CALLING PROGRAM.  IT IS
C                  INVOKED BY
C                       CALL CALCJ(N,P,X,NF,J,UIPARM,URPARM,UFPARM)
C                  NF IS THE INVOCATION COUNT FOR CALCR AT THE TIME
C                  R(X) WAS EVALUATED.  THE X PASSED TO CALCJ IS
C                  USUALLY THE ONE PASSED TO CALCR ON EITHER ITS MOST
C                  RECENT INVOCATION OR THE ONE PRIOR TO IT.  IF CALCR
C                  SAVES INTERMEDIATE RESULTS FOR USE BY CALCJ, THEN IT
C                  IS POSSIBLE TO TELL FROM NF WHETHER THEY ARE VALID
C                  FOR THE CURRENT X (OR WHICH COPY IS VALID IF TWO
C                  COPIES ARE KEPT).  IF J CANNOT BE COMPUTED AT X,
C                  THEN CALCJ SHOULD SET NF TO 0.  IN THIS CASE, NL2SOL
C                  WILL RETURN WITH IV(1) = 15.  THE OTHER PARAMETERS
C                  TO CALCJ ARE AS DESCRIBED ABOVE AND BELOW.  CALCJ
C                  SHOULD NOT CHANGE N, P, OR X.
C IV....... (INPUT/OUTPUT) AN INTEGER VALUE ARRAY OF LENGTH AT LEAST
C                  60 + P THAT HELPS CONTROL THE NL2SOL ALGORITHM AND
C                  THAT IS USED TO STORE VARIOUS INTERMEDIATE QUANTI-
C                  TIES.  OF PARTICULAR INTEREST ARE THE INITIALIZATION/
C                  RETURN CODE IV(1) AND THE ENTRIES IN IV THAT CONTROL
C                  PRINTING AND LIMIT THE NUMBER OF ITERATIONS AND FUNC-
C                  TION EVALUATIONS.  SEE THE SECTION ON IV INPUT
C                  VALUES BELOW.
C V........ (INPUT/OUTPUT) A FLOATING-POINT VALUE ARRAY OF LENGTH AT
C                  LEAST 93 + N*P + 3*N + P*(3*P+33) THAT HELPS CON-
C                  TROL THE NL2SOL ALGORITHM AND THAT IS USED TO STORE
C                  VARIOUS INTERMEDIATE QUANTITIES.  OF PARTICULAR IN-
C                  TEREST ARE THE ENTRIES IN V THAT LIMIT THE LENGTH OF
C                  THE FIRST STEP ATTEMPTED (LMAX0), SPECIFY CONVER-
C                  GENCE TOLERANCES (AFCTOL, RFCTOL, XCTOL, XFTOL),
C                  AND HELP CHOOSE THE STEP SIZE USED IN COMPUTING THE
C                  COVARIANCE MATRIX (DELTA0).  SEE THE SECTION ON
C                  (SELECTED) V INPUT VALUES BELOW.
C UIPARM... (INPUT) USER INTEGER PARAMETER ARRAY PASSED WITHOUT CHANGE
C                  TO CALCR AND CALCJ.
C URPARM... (INPUT) USER FLOATING-POINT PARAMETER ARRAY PASSED WITHOUT
C                  CHANGE TO CALCR AND CALCJ.
C UFPARM... (INPUT) USER EXTERNAL SUBROUTINE OR FUNCTION PASSED WITHOUT
C                  CHANGE TO CALCR AND CALCJ.
C
C  ***  IV INPUT VALUES (FROM SUBROUTINE DFAULT)  ***
C
C IV(1)...  ON INPUT, IV(1) SHOULD HAVE A VALUE BETWEEN 0 AND 12......
C             0 AND 12 MEAN THIS IS A FRESH START.  0 MEANS THAT
C             DFAULT(IV, V) IS TO BE CALLED TO PROVIDE ALL DEFAULT
C             VALUES TO IV AND V.  12 (THE VALUE THAT DFAULT ASSIGNS TO
C             IV(1)) MEANS THE CALLER HAS ALREADY CALLED DFAULT(IV, V)
C             AND HAS POSSIBLY CHANGED SOME IV AND/OR V ENTRIES TO NON-
C             DEFAULT VALUES.  DEFAULT = 12.
C IV(COVPRT)... IV(14) = 1 MEANS PRINT A COVARIANCE MATRIX AT THE SOLU-
C             TION.  (THIS MATRIX IS COMPUTED JUST BEFORE A RETURN WITH
C             IV(1) = 3, 4, 5, 6.)
C             IV(COVPRT) = 0 MEANS SKIP THIS PRINTING.  DEFAULT = 1.
C IV(COVREQ)... IV(15) = NONZERO MEANS COMPUTE A COVARIANCE MATRIX
C             JUST BEFORE A RETURN WITH IV(1) = 3, 4, 5, 6.  IN
C             THIS CASE, AN APPROXIMATE COVARIANCE MATRIX IS OBTAINED
C             IN ONE OF SEVERAL WAYS.  LET K = ABS(IV(COVREQ)) AND LET
C             SCALE = 2*F(X)/MAX(1,N-P),  WHERE 2*F(X) IS THE RESIDUAL
C             SUM OF SQUARES.  IF K = 1 OR 2, THEN A FINITE DIFFERENCE
C             HESSIAN APPROXIMATION H IS OBTAINED.  IF H IS POSITIVE
C             DEFINITE (OR, FOR K = 3, IF THE JACOBIAN MATRIX J AT X
C             IS NONSINGULAR), THEN ONE OF THE FOLLOWING IS COMPUTED...
C                  K = 1....  SCALE * H**-1 * (J**T * J) * H**-1.
C                  K = 2....  SCALE * H**-1.
C                  K = 3....  SCALE * (J**T * J)**-1.
C             (J**T IS THE TRANSPOSE OF J, WHILE **-1 MEANS INVERSE.)
C             IF IV(COVREQ) IS POSITIVE, THEN BOTH FUNCTION AND GRAD-
C             IENT VALUES (CALLS ON CALCR AND CALCJ) ARE USED IN COM-
C             PUTING H (WITH STEP SIZES DETERMINED USING V(DELTA0) --
C             SEE BELOW), WHILE IF IV(COVREQ) IS NEGATIVE, THEN ONLY
C             FUNCTION VALUES (CALLS ON CALCR) ARE USED (WITH STEP
C             SIZES DETERMINED USING V(DLTFDC) -- SEE BELOW).  IF
C             IV(COVREQ) = 0, THEN NO ATTEMPT IS MADE TO COMPUTE A CO-
C             VARIANCE MATRIX (UNLESS IV(COVPRT) = 1, IN WHICH CASE
C             IV(COVREQ) = 1 IS ASSUMED).  SEE IV(COVMAT) BELOW.
C             DEFAULT = 1.
C IV(DTYPE).... IV(16) TELLS HOW THE SCALE VECTOR D (SEE REF. 1) SHOULD
C             BE CHOSEN.  IV(DTYPE) .GE. 1 MEANS CHOOSE D AS DESCRIBED
C             BELOW WITH V(DFAC).  IV(DTYPE) .LE. 0 MEANS THE CALLER
C             HAS CHOSEN D AND HAS STORED IT IN V STARTING AT
C             V(94 + 2*N + P*(3*P + 31)/2).  DEFAULT = 1.
C IV(INITS).... IV(25) TELLS HOW THE S MATRIX (SEE REF. 1) SHOULD BE
C             INITIALIZED.  0 MEANS INITIALIZE S TO 0 (AND START WITH
C             THE GAUSS-NEWTON MODEL).  1 AND 2 MEAN THAT THE CALLER
C             HAS STORED THE LOWER TRIANGLE OF THE INITIAL S ROWWISE IN
C             V STARTING AT V(87+2*P).  IV(INITS) = 1 MEANS START WITH
C             THE GAUSS-NEWTON MODEL, WHILE IV(INITS) = 2 MEANS START
C             WITH THE AUGMENTED MODEL (SEE REF. 1).  DEFAULT = 0.
C IV(MXFCAL)... IV(17) GIVES THE MAXIMUM NUMBER OF FUNCTION EVALUATIONS
C             (CALLS ON CALCR, EXCLUDING THOSE USED TO COMPUTE THE CO-
C             VARIANCE MATRIX) ALLOWED.  IF THIS NUMBER DOES NOT SUF-
C             FICE, THEN NL2SOL RETURNS WITH IV(1) = 9.  DEFAULT = 200.
C IV(MXITER)... IV(18) GIVES THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C             IT ALSO INDIRECTLY LIMITS THE NUMBER OF GRADIENT EVALUA-
C             TIONS (CALLS ON CALCJ, EXCLUDING THOSE USED TO COMPUTE
C             THE COVARIANCE MATRIX) TO IV(MXITER) + 1.  IF IV(MXITER)
C             ITERATIONS DO NOT SUFFICE, THEN NL2SOL RETURNS WITH
C             IV(1) = 10.  DEFAULT = 150.
C IV(OUTLEV)... IV(19) CONTROLS THE NUMBER AND LENGTH OF ITERATION SUM-
C             MARY LINES PRINTED (BY ITSMRY).  IV(OUTLEV) = 0 MEANS DO
C             NOT PRINT ANY SUMMARY LINES.  OTHERWISE, PRINT A SUMMARY
C             LINE AFTER EACH ABS(IV(OUTLEV)) ITERATIONS.  IF IV(OUTLEV)
C             IS POSITIVE, THEN SUMMARY LINES OF LENGTH 117 (PLUS CARRI-
C             AGE CONTROL) ARE PRINTED, INCLUDING THE FOLLOWING...  THE
C             ITERATION AND FUNCTION EVALUATION COUNTS, CURRENT FUNC-
C             TION VALUE (V(F) = HALF THE SUM OF SQUARES), RELATIVE
C             DIFFERENCE IN FUNCTION VALUES ACHIEVED BY THE LATEST STEP
C             (I.E., RELDF = (F0-V(F))/F0, WHERE F0 IS THE FUNCTION
C             VALUE FROM THE PREVIOUS ITERATION), THE RELATIVE FUNCTION
C             REDUCTION PREDICTED FOR THE STEP JUST TAKEN (I.E.,
C             PRELDF = V(PREDUC) / F0, WHERE V(PREDUC) IS DESCRIBED
C             BELOW), THE SCALED RELATIVE CHANGE IN X (SEE V(RELDX)
C             BELOW), THE MODELS USED IN THE CURRENT ITERATION (G =
C             GAUSS-NEWTON, S=AUGMENTED), THE MARQUARDT PARAMETER
C             STPPAR USED IN COMPUTING THE LAST STEP, THE SIZING FACTOR
C             USED IN UPDATING S, THE 2-NORM OF THE SCALE VECTOR D
C             TIMES THE STEP JUST TAKEN (SEE REF. 1), AND NPRELDF, I.E.,
C             V(NREDUC)/F0, WHERE V(NREDUC) IS DESCRIBED BELOW -- IF
C             NPRELDF IS POSITIVE, THEN IT IS THE RELATIVE FUNCTION
C             REDUCTION PREDICTED FOR A NEWTON STEP (ONE WITH
C             STPPAR = 0).  IF NPRELDF IS ZERO, EITHER THE GRADIENT
C             VANISHES (AS DOES PRELDF) OR ELSE THE AUGMENTED MODEL
C             IS BEING USED AND ITS HESSIAN IS INDEFINITE (WITH PRELDF
C             POSITIVE).  IF NPRELDF IS NEGATIVE, THEN IT IS THE NEGA-
C             OF THE RELATIVE FUNCTION REDUCTION PREDICTED FOR A STEP
C             COMPUTED WITH STEP BOUND V(LMAX0) FOR USE IN TESTING FOR
C             SINGULAR CONVERGENCE.
C                  IF IV(OUTLEV) IS NEGATIVE, THEN LINES OF MAXIMUM
C             LENGTH 79 (OR 55 IS IV(COVPRT) = 0) ARE PRINTED, INCLUD-
C             ING ONLY THE FIRST 6 ITEMS LISTED ABOVE (THROUGH RELDX).
C             DEFAULT = 1.
C IV(PARPRT)... IV(20) = 1 MEANS PRINT ANY NONDEFAULT V VALUES ON A
C             FRESH START OR ANY CHANGED V VALUES ON A RESTART.
C             IV(PARPRT) = 0 MEANS SKIP THIS PRINTING.  DEFAULT = 1.
C IV(PRUNIT)... IV(21) IS THE OUTPUT UNIT NUMBER ON WHICH ALL PRINTING
C             IS DONE.  IV(PRUNIT) = 0 MEANS SUPPRESS ALL PRINTING.
C             (SETTING IV(PRUNIT) TO 0 IS THE ONLY WAY TO SUPPRESS THE
C             ONE LINE TERMINATION REASON MESSAGE PRINTED BY ITSMRY.)
C             DEFAULT = STANDARD OUTPUT UNIT (UNIT 6 ON MOST SYSTEMS).
C IV(SOLPRT)... IV(22) = 1 MEANS PRINT OUT THE VALUE OF X RETURNED (AS
C             WELL AS THE CORRESPONDING GRADIENT AND SCALE VECTOR D).
C             IV(SOLPRT) = 0 MEANS SKIP THIS PRINTING.  DEFAULT = 1.
C IV(STATPR)... IV(23) = 1 MEANS PRINT SUMMARY STATISTICS UPON RETURN-
C             ING.  THESE CONSIST OF THE FUNCTION VALUE (HALF THE SUM
C             OF SQUARES) AT X, V(RELDX) (SEE BELOW), THE NUMBER OF
C             FUNCTION AND GRADIENT EVALUATIONS (CALLS ON CALCR AND
C             CALCJ RESPECTIVELY, EXCLUDING ANY CALLS USED TO COMPUTE
C             THE COVARIANCE), THE RELATIVE FUNCTION REDUCTIONS PREDICT-
C             ED FOR THE LAST STEP TAKEN AND FOR A NEWTON STEP (OR PER-
C             HAPS A STEP BOUNDED BY V(LMAX0) -- SEE THE DESCRIPTIONS
C             OF PRELDF AND NPRELDF UNDER IV(OUTLEV) ABOVE), AND (IF AN
C             ATTEMPT WAS MADE TO COMPUTE THE COVARIANCE) THE NUMBER OF
C             CALLS ON CALCR AND CALCJ USED IN TRYING TO COMPUTE THE
C             COVARIANCE.  IV(STATPR) = 0 MEANS SKIP THIS PRINTING.
C             DEFAULT = 1.
C IV(X0PRT).... IV(24) = 1 MEANS PRINT THE INITIAL X AND SCALE VECTOR D
C             (ON A FRESH START ONLY).  IV(X0PRT) = 0 MEANS SKIP THIS
C             PRINTING.  DEFAULT = 1.
C
C  ***  (SELECTED) IV OUTPUT VALUES  ***
C
C IV(1)........ ON OUTPUT, IV(1) IS A RETURN CODE....
C             3 = X-CONVERGENCE.  THE SCALED RELATIVE DIFFERENCE 
C                  BETWEEN THE CURRENT PARAMETER VECTOR X AND A LOCALLY
C                  OPTIMAL PARAMETER VECTOR IS VERY LIKELY AT MOST
C                  V(XCTOL).
C             4 = RELATIVE FUNCTION CONVERGENCE.  THE RELATIVE DIFFER-
C                  ENCE BETWEEN THE CURRENT FUNCTION VALUE AND ITS LO-
C                  CALLY OPTIMAL VALUE IS VERY LIKELY AT MOST V(RFCTOL).
C             5 = BOTH X- AND RELATIVE FUNCTION CONVERGENCE (I.E., THE
C                  CONDITIONS FOR IV(1) = 3 AND IV(1) = 4 BOTH HOLD).
C             6 = ABSOLUTE FUNCTION CONVERGENCE.  THE CURRENT FUNCTION
C                  VALUE IS AT MOST V(AFCTOL) IN ABSOLUTE VALUE.
C             7 = SINGULAR CONVERGENCE.  THE HESSIAN NEAR THE CURRENT
C                  ITERATE APPEARS TO BE SINGULAR OR NEARLY SO, AND A
C                  STEP OF LENGTH AT MOST V(LMAX0) IS UNLIKELY TO YIELD
C                  A RELATIVE FUNCTION DECREASE OF MORE THAN V(RFCTOL).
C             8 = FALSE CONVERGENCE.  THE ITERATES APPEAR TO BE CONVERG-
C                  ING TO A NONCRITICAL POINT.  THIS MAY MEAN THAT THE
C                  CONVERGENCE TOLERANCES (V(AFCTOL), V(RFCTOL),
C                  V(XCTOL)) ARE TOO SMALL FOR THE ACCURACY TO WHICH
C                  THE FUNCTION AND GRADIENT ARE BEING COMPUTED, THAT
C                  THERE IS AN ERROR IN COMPUTING THE GRADIENT, OR THAT
C                  THE FUNCTION OR GRADIENT IS DISCONTINUOUS NEAR X.
C             9 = FUNCTION EVALUATION LIMIT REACHED WITHOUT OTHER CON-
C                  VERGENCE (SEE IV(MXFCAL)).
C            10 = ITERATION LIMIT REACHED WITHOUT OTHER CONVERGENCE
C                  (SEE IV(MXITER)).
C            11 = STOPX RETURNED .TRUE. (EXTERNAL INTERUPT).  SEE THE
C                  USAGE NOTES BELOW.
C            13 = F(X) CANNOT BE COMPUTED AT THE INITIAL X.
C            14 = BAD PARAMETERS PASSED TO ASSESS (WHICH SHOULD NOT
C                  OCCUR).
C            15 = THE JACOBIAN COULD NOT BE COMPUTED AT X (SEE CALCJ
C                  ABOVE).
C            16 = N OR P (OR PARAMETER NN TO NL2ITR) OUT OF RANGE --
C                  P .LE. 0 OR N .LT. P OR NN .LT. N.
C            17 = RESTART ATTEMPTED WITH N OR P (OR PAR. NN TO NL2ITR)
C                  CHANGED.
C            18 = IV(INITS) IS OUT OF RANGE.
C            19...45 = V(IV(1)-18) IS OUT OF RANGE.
C            50 = IV(1) WAS OUT OF RANGE.
C            87...(86+P) = JTOL(IV(1)-86) (I.E., V(IV(1)) IS NOT
C                  POSITIVE (SEE V(DFAC) BELOW).
C IV(COVMAT)... IV(26) TELLS WHETHER A COVARIANCE MATRIX WAS COMPUTED.
C             IF (IV(COVMAT) IS POSITIVE, THEN THE LOWER TRIANGLE OF
C             THE COVARIANCE MATRIX IS STORED ROWWISE IN V STARTING AT
C             V(IV(COVMAT)).  IF IV(COVMAT) = 0, THEN NO ATTEMPT WAS
C             MADE TO COMPUTE THE COVARIANCE.  IF IV(COVMAT) = -1,
C             THEN THE FINITE DIFFERENCE HESSIAN WAS INDEFINITE.  AND
C             AND IF IV(COVMAT) = -2, THEN A SUCCESSFUL FINITE DIFFER-
C             ENCING STEP COULD NOT BE FOUND FOR SOME COMPONENT OF X
C             (I.E., CALCR SET NF TO 0 FOR EACH OF TWO TRIAL STEPS).
C             NOTE THAT IV(COVMAT) IS RESET TO 0 AFTER EACH SUCCESSFUL
C             STEP, SO IF SUCH A STEP IS TAKEN AFTER A RESTART, THEN
C             THE COVARIANCE MATRIX WILL BE RECOMPUTED.
C IV(D)........ IV(27) IS THE STARTING SUBSCRIPT IN V OF THE CURRENT
C             SCALE VECTOR D.
C IV(G)........ IV(28) IS THE STARTING SUBSCRIPT IN V OF THE CURRENT
C             LEAST-SQUARES GRADIENT VECTOR (J**T)*R.
C IV(NFCALL)... IV(6) IS THE NUMBER OF CALLS SO FAR MADE ON CALCR (I.E.,
C             FUNCTION EVALUATIONS, INCLUDING THOSE USED IN COMPUTING
C             THE COVARIANCE).
C IV(NFCOV).... IV(40) IS THE NUMBER OF CALLS MADE ON CALCR WHEN
C             TRYING TO COMPUTE COVARIANCE MATRICES.
C IV(NGCALL)... IV(30) IS THE NUMBER OF GRADIENT EVALUATIONS (CALLS ON
C             CALCJ) SO FAR DONE (INCLUDING THOSE USED FOR COMPUTING
C             THE COVARIANCE).
C IV(NGCOV).... IV(41) IS THE NUMBER OF CALLS MADE ON CALCJ WHEN
C             TRYING TO COMPUTE COVARIANCE MATRICES.
C IV(NITER).... IV(31) IS THE NUMBER OF ITERATIONS PERFORMED.
C IV(R)........ IV(50) IS THE STARTING SUBSCRIPT IN V OF THE RESIDUAL
C             VECTOR R CORRESPONDING TO X.
C
C  ***  (SELECTED) V INPUT VALUES (FROM SUBROUTINE DFAULT)  ***
C
C V(AFCTOL)... V(31) IS THE ABSOLUTE FUNCTION CONVERGENCE TOLERANCE.
C             IF NL2SOL FINDS A POINT WHERE THE FUNCTION VALUE (HALF
C             THE SUM OF SQUARES) IS LESS THAN V(AFCTOL), AND IF NL2SOL
C             DOES NOT RETURN WITH IV(1) = 3, 4, OR 5, THEN IT RETURNS
C             WITH IV(1) = 6.  DEFAULT = MAX(10**-20, MACHEP**2), WHERE
C             MACHEP IS THE UNIT ROUNDOFF.
C V(DELTA0)... V(44) IS A FACTOR USED IN CHOOSING THE FINITE DIFFERENCE
C             STEP SIZE USED IN COMPUTING THE COVARIANCE MATRIX WHEN
C             IV(COVREQ) = 1 OR 2.  FOR COMPONENT I, STEP SIZE
C                  V(DELTA0) * MAX(ABS(X(I)), 1/D(I)) * SIGN(X(I))
C             IS USED, WHERE D IS THE CURRENT SCALE VECTOR (SEE REF. 1).
C             (IF THIS STEP RESULTS IN CALCR SETTING NF TO 0, THEN -0.5
C             TIMES THIS STEP IS ALSO TRIED.)  DEFAULT = MACHEP**0.5,
C             WHERE MACHEP IS THE UNIT ROUNDOFF.
C V(DFAC)..... V(41) AND THE D0 AND JTOL ARRAYS (SEE V(D0INIT) AND
C             V(JTINIT)) ARE USED IN UPDATING THE SCALE VECTOR D WHEN
C             IV(DTYPE) .GT. 0.  (D IS INITIALIZED ACCORDING TO
C             V(DINIT).)  LET D1(I) =
C               MAX(SQRT(JCNORM(I)**2 + MAX(S(I,I),0)), V(DFAC)*D(I)),
C             WHERE JCNORM(I) IS THE 2-NORM OF THE I-TH COLUMN OF THE
C             CURRENT JACOBIAN MATRIX AND S IS THE S MATRIX OF REF. 1.
C             IF IV(DTYPE) = 1, THEN D(I) IS SET TO D1(I) UNLESS
C             D1(I) .LT. JTOL(I), IN WHICH CASE D(I) IS SET TO
C                                MAX(D0(I), JTOL(I)).
C             IF IV(DTYPE) .GE. 2, THEN D IS UPDATED DURING THE FIRST
C             ITERATION AS FOR IV(DTYPE) = 1 (AFTER ANY INITIALIZATION
C             DUE TO V(DINIT)) AND IS LEFT UNCHANGED THEREAFTER.
C             DEFAULT = 0.6.
C V(DINIT).... V(38), IF NONNEGATIVE, IS THE VALUE TO WHICH THE SCALE
C             VECTOR D IS INITIALIZED.  DEFAULT = 0.
C V(DLTFDC)... V(40) HELPS CHOOSE THE STEP SIZE USED WHEN COMPUTING THE
C             COVARIANCE MATRIX WHEN IV(COVREQ) = -1 OR -2.  FOR
C             DIFFERENCES INVOLVING X(I), THE STEP SIZE FIRST TRIED IS
C                       V(DLTFDC) * MAX(ABS(X(I)), 1/D(I)),
C             WHERE D IS THE CURRENT SCALE VECTOR (SEE REF. 1).  (IF
C             THIS STEP IS TOO BIG THE FIRST TIME IT IS TRIED, I.E., IF
C             CALCR SETS NF TO 0, THEN -0.5 TIMES THIS STEP IS ALSO
C             TRIED.)  DEFAULT = MACHEP**(1/3), WHERE MACHEP IS THE
C             UNIT ROUNDOFF.
C V(D0INIT)... V(37), IF POSITIVE, IS THE VALUE TO WHICH ALL COMPONENTS
C             OF THE D0 VECTOR (SEE V(DFAC)) ARE INITIALIZED.  IF
C             V(DFAC) = 0, THEN IT IS ASSUMED THAT THE CALLER HAS
C             STORED D0 IN V STARTING AT V(P+87).  DEFAULT = 1.0.
C V(JTINIT)... V(39), IF POSITIVE, IS THE VALUE TO WHICH ALL COMPONENTS
C             OF THE JTOL ARRAY (SEE V(DFAC)) ARE INITIALIZED.  IF
C             V(JTINIT) = 0, THEN IT IS ASSUMED THAT THE CALLER HAS
C             STORED JTOL IN V STARTING AT V(87).  DEFAULT = 10**-6.
C V(LMAX0).... V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
C             VERY FIRST STEP THAT NL2SOL ATTEMPTS.  IT IS ALSO USED
C             IN TESTING FOR SINGULAR CONVERGENCE -- IF THE FUNCTION
C             REDUCTION PREDICTED FOR A STEP OF LENGTH BOUNDED BY
C             V(LMAX0) IS AT MOST V(RFCTOL) * ABS(F0), WHERE  F0  IS
C             THE FUNCTION VALUE AT THE START OF THE CURRENT ITERATION,
C             AND IF NL2SOL DOES NOT RETURN WITH IV(1) = 3, 4, 5, OR 6,
C             THEN IT RETURNS WITH IV(1) = 7.    DEFAULT = 100.
C V(RFCTOL)... V(32) IS THE RELATIVE FUNCTION CONVERGENCE TOLERANCE.
C             IF THE CURRENT MODEL PREDICTS A MAXIMUM POSSIBLE FUNCTION
C             REDUCTION (SEE V(NREDUC)) OF AT MOST V(RFCTOL)*ABS(F0) AT
C             THE START OF THE CURRENT ITERATION, WHERE  F0  IS THE
C             THEN CURRENT FUNCTION VALUE, AND IF THE LAST STEP ATTEMPT-
C             ED ACHIEVED NO MORE THAN TWICE THE PREDICTED FUNCTION
C             DECREASE, THEN NL2SOL RETURNS WITH IV(1) = 4 (OR 5).
C             DEFAULT = MAX(10**-10, MACHEP**(2/3)), WHERE MACHEP IS
C             THE UNIT ROUNDOFF.
C V(TUNER1)... V(26) HELPS DECIDE WHEN TO CHECK FOR FALSE CONVERGENCE
C             AND TO CONSIDER SWITCHING MODELS.  THIS IS DONE IF THE
C             ACTUAL FUNCTION DECREASE FROM THE CURRENT STEP IS NO MORE
C             THAN V(TUNER1) TIMES ITS PREDICTED VALUE.  DEFAULT = 0.1.
C V(XCTOL).... V(33) IS THE X-CONVERGENCE TOLERANCE.  IF A NEWTON STEP
C             (SEE V(NREDUC)) IS TRIED THAT HAS V(RELDX) .LE. V(XCTOL)
C             AND IF THIS STEP YIELDS AT MOST TWICE THE PREDICTED FUNC-
C             TION DECREASE, THEN NL2SOL RETURNS WITH IV(1) = 3 (OR 5).
C             (SEE THE DESCRIPTION OF V(RELDX) BELOW.)
C             DEFAULT = MACHEP**0.5, WHERE MACHEP IS THE UNIT ROUNDOFF.
C V(XFTOL).... V(34) IS THE FALSE CONVERGENCE TOLERANCE.  IF A STEP IS
C             TRIED THAT GIVES NO MORE THAN V(TUNER1) TIMES THE PREDICT-
C             ED FUNCTION DECREASE AND THAT HAS V(RELDX) .LE. V(XFTOL),
C             AND IF NL2SOL DOES NOT RETURN WITH IV(1) = 3, 4, 5, 6, OR
C             7, THEN IT RETURNS WITH IV(1) = 8.  (SEE THE DESCRIPTION
C             OF V(RELDX) BELOW.)  DEFAULT = 100*MACHEP, WHERE
C             MACHEP IS THE UNIT ROUNDOFF.
C V(*)........ DFAULT SUPPLIES TO V A NUMBER OF TUNING CONSTANTS, WITH
C             WHICH IT SHOULD ORDINARILY BE UNNECESSARY TO TINKER.  SEE
C             VERSION 2.2 OF THE NL2SOL USAGE SUMMARY (WHICH IS AN
C             APPENDIX TO REF. 1).
C
C  ***  (SELECTED) V OUTPUT VALUES  ***
C
C V(DGNORM)... V(1) IS THE 2-NORM OF (D**-1)*G, WHERE G IS THE MOST 
C             RECENTLY COMPUTED GRADIENT AND D IS THE CORRESPONDING SCALE
C             VECTOR.
C V(DSTNRM)... V(2) IS THE 2-NORM OF D*STEP, WHERE STEP IS THE MOST 
C             RECENTLY COMPUTED STEP AND D IS THE CURRENT SCALE VECTOR.
C V(F)........ V(10) IS THE CURRENT FUNCTION VALUE (HALF THE SUM OF
C             SQUARES).
C V(F0)....... V(13) IS THE FUNCTION VALUE AT THE START OF THE CURRENT
C             ITERATION.
C V(NREDUC)... V(6), IF POSITIVE, IS THE MAXIMUM FUNCTION REDUCTION
C             POSSIBLE ACCORDING TO THE CURRENT MODEL, I.E., THE FUNC-
C             TION REDUCTION PREDICTED FOR A NEWTON STEP (I.E.,
C             STEP = -H**-1 * G,  WHERE  G = (J**T) * R  IS THE CURRENT
C             GRADIENT AND H IS THE CURRENT HESSIAN APPROXIMATION --
C             H = (J**T)*J  FOR THE GAUSS-NEWTON MODEL AND
C             H = (J**T)*J + S  FOR THE AUGMENTED MODEL).
C                  V(NREDUC) = ZERO MEANS H IS NOT POSITIVE DEFINITE.
C                  IF V(NREDUC) IS NEGATIVE, THEN IT IS THE NEGATIVE OF
C             THE FUNCTION REDUCTION PREDICTED FOR A STEP COMPUTED WITH
C             A STEP BOUND OF V(LMAX0) FOR USE IN TESTING FOR SINGULAR
C             CONVERGENCE.
C V(PREDUC)... V(7) IS THE FUNCTION REDUCTION PREDICTED (BY THE CURRENT
C             QUADRATIC MODEL) FOR THE CURRENT STEP.  THIS (DIVIDED BY
C             V(F0)) IS USED IN TESTING FOR RELATIVE FUNCTION
C             CONVERGENCE.
C V(RELDX).... V(17) IS THE SCALED RELATIVE CHANGE IN X CAUSED BY THE
C             CURRENT STEP, COMPUTED AS
C                  MAX(ABS(D(I)*(X(I)-X0(I)), 1 .LE. I .LE. P) /
C                     MAX(D(I)*(ABS(X(I))+ABS(X0(I))), 1 .LE. I .LE. P),
C             WHERE X = X0 + STEP.
C
C  NOTES
C
C  ***  ALGORITHM NOTES  ***
C
C        SEE REF. 1 FOR A DESCRIPTION OF THE ALGORITHM USED.
C        ON PROBLEMS WHICH ARE NATURALLY WELL SCALED, BETTER PERFORM-
C     ANCE MAY BE OBTAINED BY SETTING V(D0INIT) = 1.0 AND IV(DTYPE) = 0,
C     WHICH WILL CAUSE THE SCALE VECTOR D TO BE SET TO ALL ONES.
C
C  ***  USAGE NOTES  ***
C
C        AFTER A RETURN WITH IV(1) .LE. 11, IT IS POSSIBLE TO RESTART,
C     I.E., TO CHANGE SOME OF THE IV AND V INPUT VALUES DESCRIBED ABOVE
C     AND CONTINUE THE ALGORITHM FROM THE POINT WHERE IT WAS INTERRUPT-
C     ED.  IV(1) SHOULD NOT BE CHANGED, NOR SHOULD ANY ENTRIES OF IV
C     AND V OTHER THAN THE INPUT VALUES (THOSE SUPPLIED BY DFAULT).
C        THOSE WHO DO NOT WISH TO WRITE A CALCJ WHICH COMPUTES THE JA-
C     COBIAN MATRIX ANALYTICALLY SHOULD CALL NL2SNO RATHER THAN NL2SOL.
C     NL2SNO USES FINITE DIFFERENCES TO COMPUTE AN APPROXIMATE JACOBIAN.
C        THOSE WHO WOULD PREFER TO PROVIDE R AND J (THE RESIDUAL AND
C     JACOBIAN) BY REVERSE COMMUNICATION RATHER THAN BY WRITING SUBROU-
C     TINES CALCR AND CALCJ MAY CALL ON NL2ITR DIRECTLY.  SEE THE COM-
C     MENTS AT THE BEGINNING OF NL2ITR.
C        THOSE WHO USE NL2SOL INTERACTIVELY MAY WISH TO SUPPLY THEIR
C     OWN STOPX FUNCTION, WHICH SHOULD RETURN .TRUE. IF THE BREAK KEY
C     HAS BEEN PRESSED SINCE STOPX WAS LAST INVOKED.  THIS MAKES IT POS-
C     SIBLE TO EXTERNALLY INTERRUPT NL2SOL (WHICH WILL RETURN WITH
C     IV(1) = 11 IF STOPX RETURNS .TRUE.).
C        STORAGE FOR J IS ALLOCATED AT THE END OF V.  THUS THE CALLER
C     MAY MAKE V LONGER THAN SPECIFIED ABOVE AND MAY ALLOW CALCJ TO USE
C     ELEMENTS OF J BEYOND THE FIRST N*P AS SCRATCH STORAGE.
C
C  ***  PORTABILITY NOTES  ***
C
C        THE NL2SOL DISTRIBUTION TAPE CONTAINS BOTH SINGLE AND DOUBLE
C     PRECISION VERSIONS OF THE NL2SOL SOURCE CODE, SO IT SHOULD BE UN-
C     NECESSARY TO CHANGE PRECISIONS.
C        ONLY THE FUNCTIONS IMDCON AND RMDCON CONTAIN MACHINE DEPENDENT
C     CONSTANTS.  TO CHANGE FROM ONE MACHINE TO ANOTHER, IT SHOULD
C     SUFFICE TO CHANGE THE (FEW) RELEVANT LINES IN THESE FUNCTIONS.
C
C  ***  REFERENCES  ***
C
C 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1980), AN ADAPTIVE
C             NONLINEAR LEAST-SQUARES ALGORITHM, SUBMITTED TO ACM TRANS.
C             MATH. SOFTWARE.
C
C  ***  GENERAL  ***
C
C     CODED BY DAVID M. GAY (WINTER 1979 - WINTER 1980).
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
C     MCS-7906671.
C
C  DECLARATIONS
C
C     EXTERNAL ITSMRY, NL2ITR
C ITSMRY... PRINTS ITERATION SUMMARY AND INFO ABOUT INITIAL AND FINAL X.
C NL2ITR... REVERSE COMMUNICATION ROUTINE THAT CARRIES OUT NL2SOL ALGO-
C             RITHM.
C
C     LOGICAL STRTED
C     INTEGER D1, J1, NF, R1
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER D, J, NFCALL, NFGCAL, R, TOOBIG
C
C  ***  IV SUBSCRIPT VALUES  ***
C
      DATA NFCALL/6/, NFGCAL/7/, TOOBIG/2/
C
C  ***  V SUBSCRIPT VALUES  ***
C
      DATA D/27/, J/33/, R/50/
C
      WRITE(NLUPRI,'(/,1X,A)')
     &   '---------------------------------'
      WRITE(NLUPRI,'(1X,A)')
     &   ' Summery of the NL2SOL package'
      WRITE(NLUPRI,'(1X,A,/)')
     &   '---------------------------------'
C
      D1 = 94 + 2*N + P*(3*P + 31)/2
      IV(D) = D1
      R1 = D1 + P
      IV(R) = R1
      J1 = R1 + N
      IV(J) = J1
      STRTED = .TRUE.
      IF (IV(1) .NE. 0 .AND. IV(1) .NE. 12) GO TO 40
         STRTED = .FALSE.
         IV(NFCALL) = 1
         IV(NFGCAL) = 1
C
 10   NF = IV(NFCALL)
C
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In NL2SOL 1: A, B & C :',x(1),x(2),x(3)
CPFP 3
c SUBROUTINE CALCR(m, n, xstp, ystp, x, iflag, fvec, NLUPRI)
        CALL CALCR(N, P, XSTP, YSTP, X,   NF, V(R1), NLUPRI)
c
c      do i = 1, N
c        WRITE(NLUPRI,'(1X,A,3(F10.2))')
c     &      ' In NL2SOL 1a: CALCR :',XSTP(i),YSTP(i),V(R1+i)
c      end do
c
      IF (STRTED) GO TO 20
         IF (NF .GT. 0) GO TO 30
              IV(1) = 13
              GO TO 60
C
 20   IF (NF .LE. 0) IV(TOOBIG) = 1
      GO TO 40
CPFP 1
C SUBROUTINE CALCJ(m, n, xstp, ystp, x, iflag, fjac, ldfjac, NLUPRI)
 30     CALL CALCJ(N, P, XSTP, YSTP, X, IV(NFGCAL), V(J1), P, NLUPRI)
c
c        WRITE(NLUPRI,'(/,1X,A,I5,/)')
c     &      ' In NL2SOL 1b: IV(NFGCAL):',IV(NFGCAL)
c
c      do i = 1, N
c        WRITE(NLUPRI,'(1X,A,3(F10.2))')
c     &      ' In NL2SOL 1c CALCJ:',V(J1+i),V(J1+N+i),V(J1+2*N+i)
c      end do
c
      IF (IV(NFGCAL) .EQ. 0) GO TO 50
      STRTED = .TRUE.
C
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In NL2SOL 2: A, B & C :',x(1),x(2),x(3)
C
 40   CALL NL2ITR(V(D1), IV, V(J1), N, N, P, V(R1), V, X, NLUPRI)
C
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In NL2SOL 1d CALCJ:',V(D1+1),V(D1+2),V(D1+3)
c
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In NL2SOL 1d CALCJ:',X(1),X(2),X(3)
c
c      do i = 1, N
c        WRITE(NLUPRI,'(1X,A,4(F10.2))')
c     &    ' In NL2SOL 1e CALCJ:',V(R1+i),V(J1+i),V(J1+N+i),V(J1+2*N+i)
c      end do
c
      IF (IV(1) - 2) 10, 30, 999
C
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In NL2SOL 4: A, B & C :',x(1),x(2),x(3)
C
 50   IV(1) = 15
 60   CALL ITSMRY(V(D1), IV, P, V, X, NLUPRI)
C
c        WRITE(NLUPRI,'(/,1X,A,3(F10.2),/)')
c     &      ' In NL2SOL 5: A, B & C :',x(1),x(2),x(3)
C
 999  WRITE(NLUPRI,'(/,1X,A)')
     &   '-------------------------------------------------------'
      WRITE(NLUPRI,'(1X,2(A,I8))')
     &   ' Leave NL2SOL after',IV(31),' iter. with INFO equal to',IV(1)
      WRITE(NLUPRI,'(1X,A,/)')
     &   '-------------------------------------------------------'
C
      RETURN
      END
c      SUBROUTINE NL2X
c
cc*********************************************************************72
cc
cC     ***  TEST NL2SOL AND NL2SNO ON MADSEN EXAMPLE  ***
cC
cC  VARIABLE DECLARATIONS
cC
cC  LOCAL ARRAYS
c      DOUBLE PRECISION
c     +   URPARM(1),V(147),X(2)
c      INTEGER
c     +   IV(62),UIPARM(1)
cC
cC  EXTERNAL SUBROUTINES
c      EXTERNAL MADJ,MADR,NL2SNO,NL2SOL,UFPARM
cC      EXTERNAL NL2SNO,NL2SOL
cC
c      X(1) = 3.0D+00
c      X(2) = 1.0D+00
c      IV(1) = 0
c      CALL NL2SOL(3, 2, X, MADR, MADJ, IV, V, UIPARM, URPARM, UFPARM)
c      IV(1) = 12
c      X(1) = 3.0D+00
c      X(2) = 1.0D+00
c      CALL NL2SNO(3, 2, X, MADR, IV, V, UIPARM, URPARM, UFPARM)
c      RETURN
c      END
      SUBROUTINE PARCHK(IV, N, NN, P, V)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  CHECK NL2SOL (VERSION 2.2) PARAMETERS, PRINT CHANGED VALUES  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,NN,P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   V(*)
      INTEGER
     +   IV(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   BIG,MACHEP,TINY,VK,ZERO
      INTEGER
     +   D0INIT,DTYPE,DTYPE0,EPSLON,I,ICH,INITS,IV1,JTINIT,JTOL0,
     +   JTOL1,JTOLP,K,L,M,NVDFLT,OLDN,OLDNN,OLDP,PARPRT,PARSV1,
     +   PRUNIT,PU
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   VM(27),VX(27)
      CHARACTER
     +   CNGD(12)*1,DFLT(12)*1,VN(8,27)*1,WHICH(12)*1
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   RMDCON
      EXTERNAL RMDCON
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DFAULT,VCOPY
C
C     INTEGER IV(21), N, NN, P
C     DOUBLE PRECISION V(33)
C     DIMENSION IV(*), V(*)
C
C     EXTERNAL DFAULT, RMDCON, VCOPY
C     DOUBLE PRECISION RMDCON
C DFAULT -- SUPPLIES DFAULT PARAMETER VALUES.
C RMDCON -- RETURNS MACHINE DEPENDENT CONSTANTS.
C VCOPY  -- COPIES ONE VECTOR TO ANOTHER.
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER I, IV1, JTOLP, K, L, M, NVDFLT, PU
C     CHARACTER*1 CNGD(12), WHICH(12)
C     CHARACTER*1 DFLT(12), VN(8,27)
C     DOUBLE PRECISION BIG, MACHEP, TINY, VK, VM(27), VX(27), ZERO
C
C  ***  IV AND V SUBSCRIPTS  ***
C
C     INTEGER DTYPE, DTYPE0, D0INIT, EPSLON, INITS, JTINIT, JTOL0,
C    1        JTOL1, OLDN, OLDNN, OLDP, PARPRT, PARSV1, PRUNIT
C
      DATA BIG/0.0D+00/, NVDFLT/27/, TINY/1.0D+00/, ZERO/0.0D+00/
C
      DATA DTYPE/16/, DTYPE0/29/, D0INIT/37/, EPSLON/19/,
     +     INITS/25/, JTINIT/39/, JTOL0/86/, JTOL1/87/,
     +     OLDN/45/, OLDNN/46/, OLDP/47/, PARPRT/20/,
     +     PARSV1/51/, PRUNIT/21/
C
      DATA
     + VN(1,1),VN(2,1),VN(3,1),VN(4,1),VN(5,1),VN(6,1),VN(7,1),VN(8,1)
     +   /'E',    'P',    'S',    'L',    'O',    'N',    '.',    '.'/
      DATA
     + VN(1,2),VN(2,2),VN(3,2),VN(4,2),VN(5,2),VN(6,2),VN(7,2),VN(8,2)
     +   /'P',    'H',    'M',    'N',    'F',    'C',    '.',    '.'/
      DATA
     + VN(1,3),VN(2,3),VN(3,3),VN(4,3),VN(5,3),VN(6,3),VN(7,3),VN(8,3)
     +   /'P',    'H',    'M',    'X',    'F',    'C',    '.',    '.'/
      DATA
     + VN(1,4),VN(2,4),VN(3,4),VN(4,4),VN(5,4),VN(6,4),VN(7,4),VN(8,4)
     +   /'D',    'E',    'C',    'F',    'A',    'C',    '.',    '.'/
      DATA
     + VN(1,5),VN(2,5),VN(3,5),VN(4,5),VN(5,5),VN(6,5),VN(7,5),VN(8,5)
     +   /'I',    'N',    'C',    'F',    'A',    'C',    '.',    '.'/
      DATA
     + VN(1,6),VN(2,6),VN(3,6),VN(4,6),VN(5,6),VN(6,6),VN(7,6),VN(8,6)
     +   /'R',    'D',    'F',    'C',    'M',    'N',    '.',    '.'/
      DATA
     + VN(1,7),VN(2,7),VN(3,7),VN(4,7),VN(5,7),VN(6,7),VN(7,7),VN(8,7)
     +   /'R',    'D',    'F',    'C',    'M',    'X',    '.',    '.'/
      DATA
     + VN(1,8),VN(2,8),VN(3,8),VN(4,8),VN(5,8),VN(6,8),VN(7,8),VN(8,8)
     +   /'T',    'U',    'N',    'E',    'R',    '1',    '.',    '.'/
      DATA
     + VN(1,9),VN(2,9),VN(3,9),VN(4,9),VN(5,9),VN(6,9),VN(7,9),VN(8,9)
     +   /'T',    'U',    'N',    'E',    'R',    '2',    '.',    '.'/
      DATA
     + VN(1,10),VN(2,10),VN(3,10),VN(4,10),VN(5,10),VN(6,10),VN(7,10),
     + VN(8,10)
     +   /'T',    'U',    'N',    'E',    'R',    '3',    '.',    '.'/
      DATA
     + VN(1,11),VN(2,11),VN(3,11),VN(4,11),VN(5,11),VN(6,11),VN(7,11),
     + VN(8,11)
     +   /'T',    'U',    'N',    'E',    'R',    '4',    '.',    '.'/
      DATA
     + VN(1,12),VN(2,12),VN(3,12),VN(4,12),VN(5,12),VN(6,12),VN(7,12),
     + VN(8,12)
     +   /'T',    'U',    'N',    'E',    'R',    '5',    '.',    '.'/
      DATA
     + VN(1,13),VN(2,13),VN(3,13),VN(4,13),VN(5,13),VN(6,13),VN(7,13),
     + VN(8,13)
     +   /'A',    'F',    'C',    'T',    'O',    'L',    '.',    '.'/
      DATA
     + VN(1,14),VN(2,14),VN(3,14),VN(4,14),VN(5,14),VN(6,14),VN(7,14),
     + VN(8,14)
     +   /'R',    'F',    'C',    'T',    'O',    'L',    '.',    '.'/
      DATA
     + VN(1,15),VN(2,15),VN(3,15),VN(4,15),VN(5,15),VN(6,15),VN(7,15),
     + VN(8,15)
     +   /'X',    'C',    'T',    'O',    'L',    '.',    '.',    '.'/
      DATA
     + VN(1,16),VN(2,16),VN(3,16),VN(4,16),VN(5,16),VN(6,16),VN(7,16),
     + VN(8,16)
     +   /'X',    'F',    'T',    'O',    'L',    '.',    '.',    '.'/
      DATA
     + VN(1,17),VN(2,17),VN(3,17),VN(4,17),VN(5,17),VN(6,17),VN(7,17),
     + VN(8,17)
     +   /'L',    'M',    'A',    'X',    '0',    '.',    '.',    '.'/
      DATA
     + VN(1,18),VN(2,18),VN(3,18),VN(4,18),VN(5,18),VN(6,18),VN(7,18),
     + VN(8,18)
     +   /'D',    'L',    'T',    'F',    'D',    'J',    '.',    '.'/
      DATA
     + VN(1,19),VN(2,19),VN(3,19),VN(4,19),VN(5,19),VN(6,19),VN(7,19),
     + VN(8,19)
     +   /'D',    '0',    'I',    'N',    'I',    'T',    '.',    '.'/
      DATA
     + VN(1,20),VN(2,20),VN(3,20),VN(4,20),VN(5,20),VN(6,20),VN(7,20),
     + VN(8,20)
     +   /'D',    'I',    'N',    'I',    'T',    '.',    '.',    '.'/
      DATA
     + VN(1,21),VN(2,21),VN(3,21),VN(4,21),VN(5,21),VN(6,21),VN(7,21),
     + VN(8,21)
     +   /'J',    'T',    'I',    'N',    'I',    'T',    '.',    '.'/
      DATA
     + VN(1,22),VN(2,22),VN(3,22),VN(4,22),VN(5,22),VN(6,22),VN(7,22),
     + VN(8,22)
     +   /'D',    'L',    'T',    'F',    'D',    'C',    '.',    '.'/
      DATA
     + VN(1,23),VN(2,23),VN(3,23),VN(4,23),VN(5,23),VN(6,23),VN(7,23),
     + VN(8,23)
     +   /'D',    'F',    'A',    'C',    '.',    '.',    '.',    '.'/
      DATA
     + VN(1,24),VN(2,24),VN(3,24),VN(4,24),VN(5,24),VN(6,24),VN(7,24),
     + VN(8,24)
     +   /'R',    'L',    'I',    'M',    'I',    'T',    '.',    '.'/
      DATA
     + VN(1,25),VN(2,25),VN(3,25),VN(4,25),VN(5,25),VN(6,25),VN(7,25),
     + VN(8,25)
     +   /'C',    'O',    'S',    'M',    'I',    'N',    '.',    '.'/
      DATA
     + VN(1,26),VN(2,26),VN(3,26),VN(4,26),VN(5,26),VN(6,26),VN(7,26),
     + VN(8,26)
     +   /'D',    'E',    'L',    'T',    'A',    '0',    '.',    '.'/
      DATA
     + VN(1,27),VN(2,27),VN(3,27),VN(4,27),VN(5,27),VN(6,27),VN(7,27),
     + VN(8,27)
     +   /'F',    'U',    'Z',    'Z',    '.',    '.',    '.',    '.'/
C
      DATA VM(1)/1.0D-03/
      DATA VM(2)/-0.99D+00/
      DATA VM(3)/1.0D-03/
      DATA VM(4)/1.0D-02/
      DATA VM(5)/1.2D+00/
      DATA VM(6)/1.0D-02/
      DATA VM(7)/1.2D+00/
      DATA VM(8)/0.0D+00/
      DATA VM(9)/0.0D+00/
      DATA VM(10)/1.0D-03/
      DATA VM(11)/-1.0D+00/
      DATA VM(15)/0.0D+00/
      DATA VM(16)/0.0D+00/
      DATA VM(19)/0.0D+00/
      DATA VM(20)/-10.0D+00/
      DATA VM(21)/0.0D+00/
      DATA VM(23)/0.0D+00/
      DATA VM(24)/1.0D+10/
      DATA VM(27)/1.01D+00/

      DATA VX(1)/0.9D+00/
      DATA VX(2)/-1.0D-03/
      DATA VX(3)/1.0D+01/
      DATA VX(4)/0.8D+00/
      DATA VX(5)/1.0D+02/
      DATA VX(6)/0.8D+00/
      DATA VX(7)/1.0D+02/
      DATA VX(8)/0.5D+00/
      DATA VX(9)/0.5D+00/
      DATA VX(10)/1.0D+00/
      DATA VX(11)/1.0D+00/
      DATA VX(14)/0.1D+00/
      DATA VX(15)/1.0D+00/
      DATA VX(16)/1.0D+00/
      DATA VX(18)/1.0D+00/
      DATA VX(22)/1.0D+00/
      DATA VX(23)/1.0D+00/
      DATA VX(25)/1.0D+00/
      DATA VX(26)/1.0D+00/
      DATA VX(27)/1.0D+02/
C
      DATA CNGD(1), CNGD(2), CNGD(3), CNGD(4), CNGD(5), CNGD(6)
     +   /     '-',     '-',     '-',     'C',     'H',     'A'/
      DATA CNGD(7), CNGD(8), CNGD(9), CNGD(10), CNGD(11), CNGD(12)
     +   /     'N',     'G',     'E',     'D',     ' ',     'V'/
      DATA DFLT(1), DFLT(2), DFLT(3), DFLT(4), DFLT(5), DFLT(6)
     +   /     'N',     'O',     'N',     'D',     'E',     'F'/
      DATA DFLT(7), DFLT(8), DFLT(9), DFLT(10), DFLT(11), DFLT(12)
     +   /     'A',     'U',     'L',     'T',     ' ',     'V'/

      IF (IV(1) .EQ. 0) CALL DFAULT(IV, V)
      PU = IV(PRUNIT)
      IV1 = IV(1)
      IF (IV1 .NE. 12) GO TO 30
         IF (NN .GE. N .AND. N .GE. P .AND. P .GE. 1) GO TO 20
              IV(1) = 16
              IF (PU .NE. 0) WRITE(PU,10) NN, N, P
 10           FORMAT(30H0///// BAD NN, N, OR P... NN =,I5,5H, N =,I5,
     +               5H, P =,I5)
              GO TO 999
 20      K = IV(21)
         CALL DFAULT(IV(21), V(33))
         IV(21) = K
         IV(DTYPE0) = IV(DTYPE+20)
         IV(OLDN) = N
         IV(OLDNN) = NN
         IV(OLDP) = P
         DO 25 ICH = 1, 12
            WHICH(ICH) = DFLT(ICH)
 25      CONTINUE
         GO TO 80
 30   IF (N .EQ. IV(OLDN) .AND. NN .EQ. IV(OLDNN) .AND. P .EQ. IV(OLDP))
     +                       GO TO 50
         IV(1) = 17
         IF (PU .NE. 0) WRITE(PU,40) IV(OLDNN), IV(OLDN), IV(OLDP), NN,
     +                               N, P
 40      FORMAT('0///// (NN,N,P) CHANGED FROM (',I5,',',I5,',',I3,
     +          ') TO (',I5,',',I5,',',I3,').')
         GO TO 999
C
 50   IF (IV1 .LE. 11 .AND. IV1 .GE. 1) GO TO 70
         IV(1) = 50
         IF (PU .NE. 0) WRITE(PU,60) IV1
 60      FORMAT('0/////  IV(1) =',I5,' SHOULD BE BETWEEN 0 AND 12.')
         GO TO 999
C
 70   DO 75 ICH = 1, 12
         WHICH(ICH) = CNGD(ICH)
 75   CONTINUE
C
 80   IF (BIG .GT. TINY) GO TO 90
         TINY = RMDCON(1)
         MACHEP = RMDCON(3)
         BIG = RMDCON(6)
         VM(12) = MACHEP
         VX(12) = BIG
         VM(13) = TINY
         VX(13) = BIG
         VM(14) = MACHEP
         VM(17) = TINY
         VX(17) = BIG
         VM(18) = MACHEP
         VX(19) = BIG
         VX(20) = BIG
         VX(21) = BIG
         VM(22) = MACHEP
         VX(24) = RMDCON(5)
         VM(25) = MACHEP
         VM(26) = MACHEP
 90   M = 0
      IF (IV(INITS) .GE. 0 .AND. IV(INITS) .LE. 2) GO TO 110
         M = 18
         IF (PU .NE. 0) WRITE(PU,100) IV(INITS)
 100     FORMAT(25H0/////  INITS... IV(25) =,I4,20H SHOULD BE BETWEEN 0,
     +          7H AND 2.)
 110  K = EPSLON
      DO 140 I = 1, NVDFLT
         VK = V(K)
         IF (VK .GE. VM(I) .AND. VK .LE. VX(I)) GO TO 130
              M = K
           IF (PU .NE. 0) WRITE(PU,120) (VN(ICH, I), ICH=1, 8),
     +                                  (VN(ICH, I), ICH=1, 8),
     +                                  K, VK, VM(I), VX(I)
 120          FORMAT(8H0/////  ,8A1,5H.. V(,I2,3H) =,E11.3,7H SHOULD,
     +               ' BE BETWEEN',E11.3,4H AND,E11.3)
 130     K = K + 1
 140     CONTINUE
C
      IF (IV1 .EQ. 12 .AND. V(JTINIT) .GT. ZERO) GO TO 170
C
C  ***  CHECK JTOL VALUES  ***
C
      JTOLP = JTOL0 + P
      DO 160 I = JTOL1, JTOLP
         IF (V(I) .GT. ZERO) GO TO 160
         K = I - JTOL0
         IF (PU .NE. 0) WRITE(PU,150) K, I, V(I)
 150     FORMAT(12H0///// JTOL(,I3,6H) = V(,I3,3H) =,E11.3,
     +          20H SHOULD BE POSITIVE.)
         M = I
 160     CONTINUE
C
 170  IF (M .EQ. 0) GO TO 180
         IV(1) = M
         GO TO 999
C
 180  IF (PU .EQ. 0 .OR. IV(PARPRT) .EQ. 0) GO TO 999
      IF (IV1 .NE. 12 .OR. IV(INITS) .EQ. 0) GO TO 200
         M = 1
         WRITE(PU,190) IV(INITS)
 190     FORMAT(22H0NONDEFAULT VALUES..../20H INITS..... IV(25) =,I3)
 200  IF (IV(DTYPE) .EQ. IV(DTYPE0)) GO TO  210
         IF (M .EQ. 0) WRITE(PU,215) (WHICH(ICH), ICH=1, 12)
         M = 1
         WRITE(PU,205) IV(DTYPE)
 205     FORMAT(20H DTYPE..... IV(16) =,I3)
 210  K = EPSLON
      L = PARSV1
      DO 240 I = 1, NVDFLT
         IF (V(K) .EQ. V(L)) GO TO 230
              IF (M .EQ. 0) WRITE(PU,215) (WHICH(ICH), ICH = 1, 12)
 215          FORMAT ('0',12A1,'ALUES....'/)
              M = 1
              WRITE (PU,220) (VN(ICH, I), ICH = 1, 8), K, V(K)
 220          FORMAT (1X, 8A1, 5H.. V(, I2, 3H) =, E15.7)
 230     K = K + 1
         L = L + 1
 240     CONTINUE
      IV(DTYPE0) = IV(DTYPE)
      CALL VCOPY(NVDFLT, V(PARSV1), V(EPSLON))
      IF (IV1 .NE. 12) GO TO 999
         IF (V(JTINIT) .GT. ZERO) GO TO 260
              JTOLP = JTOL0 + P
              WRITE(PU,250) (V(I), I = JTOL1, JTOLP)
 250          FORMAT(24H0(INITIAL) JTOL ARRAY.../(1X,6E12.3))
 260     IF (V(D0INIT) .GT. ZERO) GO TO 999
              K = JTOL1 + P
              L = K + P - 1
              WRITE(PU,270) (V(I), I = K, L)
 270          FORMAT(22H0(INITIAL) D0 ARRAY.../1X,6E12.3)
C
 999  RETURN
      END
      SUBROUTINE QAPPLY(NN, N, P, J, R, IERR)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  VARIABLE DECLARATIONS
C
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IERR,N,NN,P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   J(NN,P),R(N)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   T
      INTEGER
     +   I,K,L,NL1
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD
      EXTERNAL DOTPRD
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS
C
C     *****PARAMETERS.
C     INTEGER NN, N, P, IERR
C     DOUBLE PRECISION J(NN,P), R(N)
C
C     *****PURPOSE.
C     THIS SUBROUTINE APPLIES TO R THE ORTHOGONAL TRANSFORMATIONS
C     STORED IN J BY QRFACT
C
C     *****PARAMETER DESCRIPTION.
C     ON INPUT.
C
C        NN IS THE ROW DIMENSION OF THE MATRIX J AS DECLARED IN
C             THE CALLING PROGRAM DIMENSION STATEMENT
C
C        N IS THE NUMBER OF ROWS OF J AND THE SIZE OF THE VECTOR R
C
C        P IS THE NUMBER OF COLUMNS OF J AND THE SIZE OF SIGMA
C
C        J CONTAINS ON AND BELOW ITS DIAGONAL THE COLUMN VECTORS
C             U WHICH DETERMINE THE HOUSEHOLDER TRANSFORMATIONS
C             IDENT - U*U.TRANSPOSE
C
C        R IS THE RIGHT HAND SIDE VECTOR TO WHICH THE ORTHOGONAL
C             TRANSFORMATIONS WILL BE APPLIED
C
C        IERR IF NON-ZERO INDICATES THAT NOT ALL THE TRANSFORMATIONS
C             WERE SUCCESSFULLY DETERMINED AND ONLY THE FIRST
C             ABS(IERR) - 1 TRANSFORMATIONS WILL BE USED
C
C     ON OUTPUT.
C
C        R HAS BEEN OVERWRITTEN BY ITS TRANSFORMED IMAGE
C
C     *****APPLICATION AND USAGE RESTRICTIONS.
C     NONE
C
C     *****ALGORITHM NOTES.
C     THE VECTORS U WHICH DETERMINE THE HOUSEHOLDER TRANSFORMATIONS
C     ARE NORMALIZED SO THAT THEIR 2-NORM SQUARED IS 2.  THE USE OF
C     THESE TRANSFORMATIONS HERE IS IN THE SPIRIT OF (1).
C
C     *****SUBROUTINES AND FUNCTIONS CALLED.
C
C     DOTPRD - FUNCTION, RETURNS THE INNER PRODUCT OF VECTORS
C
C     *****REFERENCES.
C     (1) BUSINGER, P. A., AND GOLUB, G. H. (1965), LINEAR LEAST SQUARES
C        SOLUTIONS BY HOUSEHOLDER TRANSFORMATIONS, NUMER. MATH. 7,
C        PP. 269-276.
C
C     *****HISTORY.
C     DESIGNED BY DAVID M. GAY, CODED BY STEPHEN C. PETERS (WINTER 1977)
C
C     *****GENERAL.
C
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, AND MCS76-11989.
C
C     *****LOCAL VARIABLES.
C     INTEGER I, K, L, NL1
C     DOUBLE PRECISION T
C/
C     *****FUNCTIONS.
C     EXTERNAL DOTPRD
C     DOUBLE PRECISION DOTPRD
C
      K = P
      IF (IERR .NE. 0) K = ABS(IERR) - 1
      IF ( K .EQ. 0) GO TO 999
C
      DO 20 L = 1, K
         NL1 = N - L + 1
         T = -DOTPRD(NL1, J(L,L), R(L))
C
         DO 10 I = L, N
            R(I) = R(I) + T*J(I,L)
 10      CONTINUE
 20   CONTINUE
 999  RETURN
      END
      SUBROUTINE QRFACT(NM,M,N,QR,ALPHA,IPIVOT,IERR,NOPIVK,SUM)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  COMPUTE THE QR DECOMPOSITION OF THE MATRIX STORED IN QR  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IERR,M,N,NM,NOPIVK
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   ALPHA(N),QR(NM,N),SUM(N)
      INTEGER
     +   IPIVOT(N)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ALPHAK,BETA,ONE,P01,P99,QRKK,QRKMAX,RKTOL,RKTOL1,SIGMA,SUMJ,
     +   TEMP,UFETA,ZERO
      INTEGER
     +   I,J,JBAR,K,K1,MINUM,MK1
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD,RMDCON,V2NORM
      EXTERNAL DOTPRD,RMDCON,V2NORM
C
C  EXTERNAL SUBROUTINES
      EXTERNAL VAXPY,VSCOPY
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MIN,SQRT
C
C     *****PARAMETERS.
C     INTEGER NM,M,N,IPIVOT(N),IERR,NOPIVK
C     DOUBLE PRECISION  QR(NM,N),ALPHA(N),SUM(N)
C     *****LOCAL VARIABLES.
C     INTEGER I,J,JBAR,K,K1,MINUM,MK1
C     DOUBLE PRECISION  ALPHAK,BETA,QRKK,QRKMAX,SIGMA,TEMP,UFETA,RKTOL,
C    1        RKTOL1,SUMJ
C     *****FUNCTIONS.
C/+
C     INTEGER MIN
C     DOUBLE PRECISION              ABS,SQRT
C/
C     EXTERNAL DOTPRD, RMDCON, VAXPY, VSCOPY, V2NORM
C     DOUBLE PRECISION DOTPRD, RMDCON, V2NORM
C DOTPRD... RETURNS INNER PRODUCT OF TWO VECTORS.
C RMDCON... RETURNS MACHINE DEPENDENT CONSTANTS.
C VAXPY... COMPUTES SCALAR TIMES ONE VECTOR PLUS ANOTHER.
C VSCOPY... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR.
C V2NORM... RETURNS THE 2-NORM OF A VECTOR.
C
C     *****CONSTANTS.
C     DOUBLE PRECISION ONE, P01, P99, ZERO
      DATA ONE/1.0D+00/, P01/0.01D+00/, P99/0.99D+00/, ZERO/0.0D+00/
C
C
C
C     *****PURPOSE.
C
C     THIS SUBROUTINE DOES A QR-DECOMPOSITION ON THE M X N MATRIX QR,
C        WITH AN OPTIONALLY MODIFIED COLUMN PIVOTING, AND RETURNS THE
C        UPPER TRIANGULAR R-MATRIX, AS WELL AS THE ORTHOGONAL VECTORS
C        USED IN THE TRANSFORMATIONS.
C
C     *****PARAMETER DESCRIPTION.
C     ON INPUT.
C
C        NM MUST BE SET TO THE ROW DIMENSION OF THE TWO DIMENSIONAL
C             ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C             DIMENSION STATEMENT.
C
C        M MUST BE SET TO THE NUMBER OF ROWS IN THE MATRIX.
C
C        N MUST BE SET TO THE NUMBER OF COLUMNS IN THE MATRIX.
C
C        QR CONTAINS THE DOUBLE PRECISION RECTANGULAR MATRIX TO BE DECOMPOSED.
C
C     NOPIVK IS USED TO CONTROL PIVOTTING.  COLUMNS 1 THROUGH
C        NOPIVK WILL REMAIN FIXED IN POSITION.
C
C        SUM IS USED FOR TEMPORARY STORAGE FOR THE SUBROUTINE.
C
C     ON OUTPUT.
C
C        QR CONTAINS THE NON-DIAGONAL ELEMENTS OF THE R-MATRIX
C             IN THE STRICT UPPER TRIANGLE. THE VECTORS U, WHICH
C             DEFINE THE HOUSEHOLDER TRANSFORMATIONS   I - U*U-TRANSP,
C             ARE IN THE COLUMNS OF THE LOWER TRIANGLE. THESE VECTORS U
C             ARE SCALED SO THAT THE SQUARE OF THEIR 2-NORM IS 2.0.
C
C        ALPHA CONTAINS THE DIAGONAL ELEMENTS OF THE R-MATRIX.
C
C        IPIVOT REFLECTS THE COLUMN PIVOTING PERFORMED ON THE INPUT
C             MATRIX TO ACCOMPLISH THE DECOMPOSITION. THE J-TH
C             ELEMENT OF IPIVOT GIVES THE COLUMN OF THE ORIGINAL
C             MATRIX WHICH WAS PIVOTED INTO COLUMN J DURING THE
C             DECOMPOSITION.
C
C        IERR IS SET TO.
C             0 FOR NORMAL RETURN,
C             K IF NO NON-ZERO PIVOT COULD BE FOUND FOR THE K-TH
C                  TRANSFORMATION, OR
C             -K FOR AN ERROR EXIT ON THE K-TH THANSFORMATION.
C             IF AN ERROR EXIT WAS TAKEN, THE FIRST (K - 1)
C             TRANSFORMATIONS ARE CORRECT.
C
C
C     *****APPLICATIONS AND USAGE RESTRICTIONS.
C     THIS MAY BE USED WHEN SOLVING LINEAR LEAST-SQUARES PROBLEMS --
C     SEE SUBROUTINE QR1 OF ROSEPACK.  IT IS CALLED FOR THIS PURPOSE
C     BY LLSQST IN THE NL2SOL (NONLINEAR LEAST-SQUARES) PACKAGE.
C
C     *****ALGORITHM NOTES.
C     THIS VERSION OF QRFACT TRIES TO ELIMINATE THE OCCURRENCE OF
C     UNDERFLOWS DURING THE ACCUMULATION OF INNER PRODUCTS.  RKTOL1
C     IS CHOSEN BELOW SO AS TO INSURE THAT DISCARDED TERMS HAVE NO
C     EFFECT ON THE COMPUTED TWO-NORMS.
C
C     ADAPTED FROM THE ALGOL ROUTINE SOLVE (1).
C
C     *****REFERENCES.
C     (1)     BUSINGER,P. AND GOLUB,G.H., LINEAR LEAST SQUARES
C     SOLUTIONS BY HOUSHOLDER TRANSFORMATIONS, IN WILKINSON,J.H.
C     AND REINSCH,C.(EDS.), HANDBOOK FOR AUTOMATIC COMPUTATION,
C     VOLUME II. LINEAR ALGEBRA, SPRINGER-VERLAG, 111-118 (1971).
C     PREPUBLISHED IN NUMER.MATH. 7, 269-276 (1965).
C
C     *****HISTORY.
C     THIS AMOUNTS TO THE SUBROUTINE QR1 OF ROSEPACK WITH RKTOL1 USED
C     IN PLACE OF RKTOL BELOW, WITH V2NORM USED TO INITIALIZE (AND
C     SOMETIMES UPDATE) THE SUM ARRAY, AND WITH CALLS ON DOTPRD AND
C     VAXPY IN PLACE OF SOME LOOPS.
C
C     *****GENERAL.
C
C     DEVELOPMENT OF THIS PROGRAM SUPPORTED IN PART BY
C     NATIONAL SCIENCE FOUNDATION GRANT GJ-1154X3 AND
C     NATIONAL SCIENCE FOUNDATION GRANT DCR75-08802
C     TO NATIONAL BUREAU OF ECONOMIC RESEARCH, INC.
C
C
C     ==========  UFETA IS THE SMALLEST POSITIVE FLOATING POINT NUMBER
C        S.T. UFETA AND -UFETA CAN BOTH BE REPRESENTED.
C
C     ==========  RKTOL IS THE SQUARE ROOT OF THE RELATIVE PRECISION
C        OF FLOATING POINT ARITHMETIC (MACHEP).
      DATA RKTOL/0.0D+00/, UFETA/0.0D+00/
C     *****BODY OF PROGRAM.
      IF (UFETA .GT. ZERO) GO TO 10
         UFETA = RMDCON(1)
         RKTOL = RMDCON(4)
   10 IERR = 0
      RKTOL1 = P01 * RKTOL
C
      DO 20 J=1,N
         SUM(J) = V2NORM(M, QR(1,J))
         IPIVOT(J) = J
   20 CONTINUE
C
      MINUM = MIN(M,N)
C
      DO 120 K=1,MINUM
         MK1 = M - K + 1
C        ==========K-TH HOUSEHOLDER TRANSFORMATION==========
         SIGMA = ZERO
         JBAR = 0
C        ==========FIND LARGEST COLUMN SUM==========
      IF (K .LE. NOPIVK) GO TO 50
         DO 30 J=K,N
              IF (SIGMA .GE. SUM(J))  GO TO 30
              SIGMA = SUM(J)
              JBAR = J
   30    CONTINUE
C
         IF (JBAR .EQ. 0)  GO TO 220
         IF (JBAR .EQ. K)  GO TO 50
C        ==========COLUMN INTERCHANGE==========
         I = IPIVOT(K)
         IPIVOT(K) = IPIVOT(JBAR)
         IPIVOT(JBAR) = I
         SUM(JBAR) = SUM(K)
         SUM(K) = SIGMA
C
         DO 40 I=1,M
              SIGMA = QR(I,K)
              QR(I,K) = QR(I,JBAR)
              QR(I,JBAR) = SIGMA
   40    CONTINUE
C        ==========END OF COLUMN INTERCHANGE==========
   50    CONTINUE
C        ==========  SECOND INNER PRODUCT  ==========
         QRKMAX = ZERO
C
         DO 60 I=K,M
              IF (ABS( QR(I,K) ) .GT. QRKMAX)  QRKMAX = ABS( QR(I,K) )
   60    CONTINUE
C
         IF (QRKMAX .LT. UFETA)  GO TO 210
         ALPHAK = V2NORM(MK1, QR(K,K)) / QRKMAX
         SIGMA = ALPHAK**2
C
C        ==========  END SECOND INNER PRODUCT  ==========
         QRKK = QR(K,K)
         IF (QRKK .GE. ZERO)  ALPHAK = -ALPHAK
         ALPHA(K) = ALPHAK * QRKMAX
         BETA = QRKMAX * SQRT(SIGMA - (QRKK*ALPHAK/QRKMAX) )
         QR(K,K) = QRKK - ALPHA(K)
         DO 65 I=K,M
   65         QR(I,K) =  QR(I,K) / BETA
         K1 = K + 1
         IF (K1 .GT. N) GO TO 120
C
         DO 110 J = K1, N
              TEMP = -DOTPRD(MK1, QR(K,K), QR(K,J))
C
C             ***  SET QR(I,J) = QR(I,J) + TEMP*QR(I,K), I = K,...,M.
C
              CALL VAXPY(MK1, QR(K,J), TEMP, QR(K,K), QR(K,J))
C
              IF (K1 .GT. M) GO TO 110
              SUMJ = SUM(J)
              IF (SUMJ .LT. UFETA) GO TO 110
              TEMP = ABS(QR(K,J)/SUMJ)
              IF (TEMP .LT. RKTOL1) GO TO 110
              IF (TEMP .GE. P99) GO TO 90
                   SUM(J) = SUMJ * SQRT(ONE - TEMP**2)
                   GO TO 110
   90         SUM(J) = V2NORM(M-K, QR(K1,J))
  110    CONTINUE
C        ==========END OF K-TH HOUSEHOLDER TRANSFORMATION==========
  120 CONTINUE
C
      GO TO 999
C     ==========ERROR EXIT ON K-TH TRANSFORMATION==========
  210 IERR = -K
      GO TO 230
C     ==========NO NON-ZERO ACCEPTABLE PIVOT FOUND==========
  220 IERR = K
  230 DO 240 I = K, N
         ALPHA(I) = ZERO
         IF (I .GT. K) CALL VSCOPY(I-K, QR(K,I), ZERO)
 240     CONTINUE
C     ==========RETURN TO CALLER==========
  999 RETURN
      END
      DOUBLE PRECISION FUNCTION RELDST(P, D, X, X0)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  COMPUTE AND RETURN RELATIVE DIFFERENCE BETWEEN X AND X0  ***
C  ***  NL2SOL VERSION 2.2  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   D(P),X(P),X0(P)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   EMAX,T,XMAX,ZERO
      INTEGER
     +   I
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS
C
      DATA ZERO/0.0D+00/
C
      EMAX = ZERO
      XMAX = ZERO
      DO 10 I = 1, P
         T = ABS(D(I) * (X(I) - X0(I)))
         IF (EMAX .LT. T) EMAX = T
         T = D(I) * (ABS(X(I)) + ABS(X0(I)))
         IF (XMAX .LT. T) XMAX = T
 10      CONTINUE
      RELDST = ZERO
      IF (XMAX .GT. ZERO) RELDST = EMAX / XMAX
      RETURN
      END
      DOUBLE PRECISION FUNCTION RMDCON(K)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   K
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   BIG,ETA,MACHEP,ONE001,PT999
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   D1MACH
      EXTERNAL D1MACH
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C
C  ***  RETURN MACHINE DEPENDENT CONSTANTS USED BY NL2SOL  ***
C
C +++  COMMENTS BELOW CONTAIN DATA STATEMENTS FOR VARIOUS MACHINES.  +++
C +++  TO CONVERT TO ANOTHER MACHINE, PLACE A C IN COLUMN 1 OF THE   +++
C +++  DATA STATEMENT LINE(S) THAT CORRESPOND TO THE CURRENT MACHINE +++
C +++  AND REMOVE THE C FROM COLUMN 1 OF THE DATA STATEMENT LINE(S)  +++
C +++  THAT CORRESPOND TO THE NEW MACHINE.                           +++
C
C     INTEGER K
C
C  ***  THE CONSTANT RETURNED DEPENDS ON K...
C
C  ***        K = 1... SMALLEST POS. ETA SUCH THAT -ETA EXISTS.
C  ***        K = 2... SQUARE ROOT OF 1.001*ETA.
C  ***        K = 3... UNIT ROUNDOFF = SMALLEST POS. NO. MACHEP SUCH
C  ***                 THAT 1 + MACHEP .GT. 1 .AND. 1 - MACHEP .LT. 1.
C  ***        K = 4... SQUARE ROOT OF 0.999*MACHEP.
C  ***        K = 5... SQUARE ROOT OF 0.999*BIG (SEE K = 6).
C  ***        K = 6... LARGEST MACHINE NO. BIG SUCH THAT -BIG EXISTS.
C
C
      DATA ONE001/1.001D+00/, PT999/0.999D+00/
C
      BIG = D1MACH(2)
      ETA = D1MACH(1)
      MACHEP = D1MACH(4)

      GO TO (10, 20, 30, 40, 50, 60), K
C
 10   RMDCON = ETA
      GO TO 999
C
 20   RMDCON = SQRT(ONE001*ETA)
      GO TO 999
C
 30   RMDCON = MACHEP
      GO TO 999
C
 40   RMDCON = SQRT(PT999*MACHEP)
      GO TO 999
C
 50   RMDCON = SQRT(PT999*BIG)
      GO TO 999
C
 60   RMDCON = BIG
C
 999  RETURN
      END
      SUBROUTINE RPTMUL(FUNC, IPIVOT, J, NN, P, RD, X, Y, Z)

c*********************************************************************72
c
C  ***  FUNC = 1... SET  Y = RMAT * (PERM**T) * X.
C  ***  FUNC = 2... SET  Y = PERM * (RMAT**T) * RMAT * (PERM**T) * X.
C  ***  FUNC = 3... SET  Y = PERM * (RMAT**T) X.
C
C
C  ***  PERM = MATRIX WHOSE I-TH COL. IS THE IPIVOT(I)-TH UNIT VECTOR.
C  ***  RMAT IS THE UPPER TRIANGULAR MATRIX WHOSE STRICT UPPER TRIANGLE
C  ***       IS STORED IN  J  AND WHOSE DIAGONAL IS STORED IN RD.
C  ***  Z IS A SCRATCH VECTOR.
C  ***  X AND Y MAY SHARE STORAGE.
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   FUNC,NN,P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   J(NN,P),RD(P),X(P),Y(P),Z(P)
      INTEGER
     +   IPIVOT(P)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ZK
      INTEGER
     +   I,IM1,K,KM1
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD
      EXTERNAL DOTPRD
C
C
C-----------------------------------------------------------------------
C
      IF (FUNC .GT. 2) GO TO 50
C
C  ***  FIRST SET  Z = (PERM**T) * X  ***
C
      DO 10 I = 1, P
         K = IPIVOT(I)
         Z(I) = X(K)
 10      CONTINUE
C
C  ***  NOW SET  Y = RMAT * Z  ***
C
      Y(1) = Z(1) * RD(1)
      IF (P .LE. 1) GO TO 40
      DO 30 K = 2, P
         KM1 = K - 1
         ZK = Z(K)
         DO 20 I = 1, KM1
 20           Y(I) = Y(I) + J(I,K)*ZK
         Y(K) = ZK*RD(K)
 30      CONTINUE
C
 40   IF (FUNC .LE. 1) GO TO 999
      GO TO 70
C
 50   DO 60 I = 1, P
 60      Y(I) = X(I)
C
C  ***  SET  Z = (RMAT**T) * Y  ***
C
 70   Z(1) = Y(1) * RD(1)
      IF (P .EQ. 1) GO TO 90
      DO 80 I = 2, P
         IM1 = I - 1
         Z(I) = Y(I)*RD(I) + DOTPRD(IM1, J(1,I), Y)
 80      CONTINUE
C
C  ***  NOW SET  Y = PERM * Z  ***
C
 90   DO 100 I = 1, P
         K = IPIVOT(I)
         Y(K) = Z(I)
 100     CONTINUE
C
 999  RETURN
      END
      SUBROUTINE SLUPDT(A, COSMIN, P, SIZE, STEP, U, W, WCHMTD, WSCALE,
     +                  Y)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  UPDATE SYMMETRIC  A  SO THAT  A * STEP = Y  ***
C  ***  (LOWER TRIANGLE OF  A  STORED ROWWISE       ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   COSMIN,SIZE,WSCALE
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   A(1),STEP(P),U(P),W(P),WCHMTD(P),Y(P)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   DENMIN,HALF,ONE,SDOTWM,T,UI,WI,ZERO
      INTEGER
     +   I,J,K
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD,V2NORM
      EXTERNAL DOTPRD,V2NORM
C
C  EXTERNAL SUBROUTINES
      EXTERNAL SLVMUL
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MIN
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER P
C     DOUBLE PRECISION A(1), COSMIN, SIZE, STEP(P), U(P), W(P),
C    1                 WCHMTD(P), WSCALE, Y(P)
C     DIMENSION A(P*(P+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER I, J, K
C     DOUBLE PRECISION DENMIN, SDOTWM, T, UI, WI
C
C     ***  CONSTANTS  ***
C     DOUBLE PRECISION HALF, ONE, ZERO
C
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL DOTPRD, SLVMUL, V2NORM
C     DOUBLE PRECISION DOTPRD, V2NORM
C
      DATA HALF/0.5D+00/, ONE/1.0D+00/, ZERO/0.0D+00/

      SDOTWM = DOTPRD(P, STEP, WCHMTD)
      DENMIN = COSMIN * V2NORM(P,STEP) * V2NORM(P,WCHMTD)
      WSCALE = ONE
      IF (DENMIN .NE. ZERO) WSCALE = MIN(ONE, ABS(SDOTWM/DENMIN))
      T = ZERO
      IF (SDOTWM .NE. ZERO) T = WSCALE / SDOTWM
      DO 10 I = 1, P
 10      W(I) = T * WCHMTD(I)
      CALL SLVMUL(P, U, A, STEP)
      T = HALF * (SIZE * DOTPRD(P, STEP, U)  -  DOTPRD(P, STEP, Y))
      DO 20 I = 1, P
 20      U(I) = T*W(I) + Y(I) - SIZE*U(I)
C
C  ***  SET  A = A + U*(W**T) + W*(U**T)  ***
C
      K = 1
      DO 40 I = 1, P
         UI = U(I)
         WI = W(I)
         DO 30 J = 1, I
              A(K) = SIZE*A(K) + UI*W(J) + WI*U(J)
              K = K + 1
 30           CONTINUE
 40      CONTINUE
C
      RETURN
      END
      SUBROUTINE SLVMUL(P, Y, S, X)

c*********************************************************************72
c
C  ***  SET  Y = S * X,  S = P X P SYMMETRIC MATRIX.  ***
C  ***  LOWER TRIANGLE OF  S  STORED ROWWISE.         ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   S(1),X(P),Y(P)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   XI
      INTEGER
     +   I,IM1,J,K
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD
      EXTERNAL DOTPRD
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER P
C     DOUBLE PRECISION S(1), X(P), Y(P)
C     DIMENSION S(P*(P+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER I, IM1, J, K
C     DOUBLE PRECISION XI
C
C  ***  EXTERNAL FUNCTION  ***
C
C     EXTERNAL DOTPRD
C     DOUBLE PRECISION DOTPRD
C
      J = 1
      DO 10 I = 1, P
         Y(I) = DOTPRD(I, S(J), X)
         J = J + I
 10      CONTINUE
C
      IF (P .LE. 1) GO TO 999
      J = 1
      DO 40 I = 2, P
         XI = X(I)
         IM1 = I - 1
         J = J + 1
         DO 30 K = 1, IM1
              Y(K) = Y(K) + S(J)*XI
              J = J + 1
 30           CONTINUE
 40      CONTINUE
C
 999  RETURN
      END
      LOGICAL FUNCTION STOPX(IDUMMY)

c*********************************************************************72
c
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IDUMMY
C
C     *****PURPOSE...
C     THIS FUNCTION MAY SERVE AS THE STOPX (ASYNCHRONOUS INTERRUPTION)
C     FUNCTION FOR THE NL2SOL (NONLINEAR LEAST-SQUARES) PACKAGE AT
C     THOSE INSTALLATIONS WHICH DO NOT WISH TO IMPLEMENT A
C     DYNAMIC STOPX.
C
C     *****ALGORITHM NOTES...
C     AT INSTALLATIONS WHERE THE NL2SOL SYSTEM IS USED
C     INTERACTIVELY, THIS DUMMY STOPX SHOULD BE REPLACED BY A
C     FUNCTION THAT RETURNS .TRUE. IF AND ONLY IF THE INTERRUPT
C     (BREAK) KEY HAS BEEN PRESSED SINCE THE LAST CALL ON STOPX.
C
C
      STOPX = .FALSE.
      RETURN
      END
      DOUBLE PRECISION FUNCTION V2NORM(P, X)

c*********************************************************************72
c
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  RETURN THE 2-NORM OF THE P-VECTOR X, TAKING  ***
C  ***  CARE TO AVOID THE MOST LIKELY UNDERFLOWS.    ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   X(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ONE,R,SCALE,SQTETA,T,XI,ZERO
      INTEGER
     +   I,J
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   RMDCON
      EXTERNAL RMDCON
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,SQRT
C
C
      DATA ONE/1.0D+00/, SQTETA/0.0D+00/, ZERO/0.0D+00/
C
      IF (P .GT. 0) GO TO 10
         V2NORM = ZERO
         GO TO 999
 10   DO 20 I = 1, P
         IF (X(I) .NE. ZERO) GO TO 30
 20      CONTINUE
      V2NORM = ZERO
      GO TO 999
C
 30   SCALE = ABS(X(I))
      IF (I .LT. P) GO TO 40
         V2NORM = SCALE
         GO TO 999
 40   T = ONE
      IF (SQTETA .EQ. ZERO) SQTETA = RMDCON(2)
C
C     ***  SQTETA IS (SLIGHTLY LARGER THAN) THE SQUARE ROOT OF THE
C     ***  SMALLEST POSITIVE FLOATING POINT NUMBER ON THE MACHINE.
C     ***  THE TESTS INVOLVING SQTETA ARE DONE TO PREVENT UNDERFLOWS.
C
      J = I + 1
      DO 60 I = J, P
         XI = ABS(X(I))
         IF (XI .GT. SCALE) GO TO 50
              R = XI / SCALE
              IF (R .GT. SQTETA) T = T + R*R
              GO TO 60
 50           R = SCALE / XI
              IF (R .LE. SQTETA) R = ZERO
              T = ONE  +  T * R*R
         SCALE = XI
 60      CONTINUE
C
      V2NORM = SCALE * SQRT(T)
 999  RETURN
      END
      SUBROUTINE VAXPY(P, W, A, X, Y)

c*********************************************************************72
c
C  ***  SET W = A*X + Y  --  W, X, Y = P-VECTORS, A = SCALAR  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   A
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   W(*),X(*),Y(*)
C
C  LOCAL SCALARS
      INTEGER
     +   I
C
C
      DO 10 I = 1, P
 10      W(I) = A*X(I) + Y(I)
      RETURN
      END
      SUBROUTINE VCOPY(P, Y, X)

c*********************************************************************72
c
C  ***  SET Y = X, WHERE X AND Y ARE P-VECTORS  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   X(*),Y(*)
C
C  LOCAL SCALARS
      INTEGER
     +   I
C
C
      DO 10 I = 1, P
 10      Y(I) = X(I)
      RETURN
      END
      SUBROUTINE VSCOPY(P, Y, S)

c*********************************************************************72
c
C  ***  SET P-VECTOR Y TO SCALAR S  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   S
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   Y(*)
C
C  LOCAL SCALARS
      INTEGER
     +   I
C
C
      DO 10 I = 1, P
 10      Y(I) = S
      RETURN
      END
      function i1mach ( i )

c*********************************************************************72
c
cc I1MACH returns integer machine dependent constants.
c
c  Discussion:
c
c    Input/output unit numbers.
c
c      I1MACH(1) = the standard input unit.
c      I1MACH(2) = the standard output unit.
c      I1MACH(3) = the standard punch unit.
c      I1MACH(4) = the standard error message unit.
c
c    Words.
c
c      I1MACH(5) = the number of bits per integer storage unit.
c      I1MACH(6) = the number of characters per integer storage unit.
c
c    Integers.
c
c    Assume integers are represented in the S digit base A form:
c
c      Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
c
c    where 0 <= X(1:S-1) < A.
c
c      I1MACH(7) = A, the base.
c      I1MACH(8) = S, the number of base A digits.
c      I1MACH(9) = A**S-1, the largest integer.
c
c    Floating point numbers
c
c    Assume floating point numbers are represented in the T digit
c    base B form:
c
c      Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
c
c    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
c
c      I1MACH(10) = B, the base.
c
c    Single precision
c
c      I1MACH(11) = T, the number of base B digits.
c      I1MACH(12) = EMIN, the smallest exponent E.
c      I1MACH(13) = EMAX, the largest exponent E.
c
c    Double precision
c
c      I1MACH(14) = T, the number of base B digits.
c      I1MACH(15) = EMIN, the smallest exponent E.
c      I1MACH(16) = EMAX, the largest exponent E.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528,
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 16.
c
c    Output, integer I1MACH, the value of the chosen parameter.
c
      implicit none

      integer i
      integer i1mach

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
        write ( *, '(a,i12)' ) '  I = ', i
        i1mach = 0
        stop
      else if ( i == 1 ) then
        i1mach = 5
      else if ( i == 2 ) then
        i1mach = 6
      else if ( i == 3 ) then
        i1mach = 7
      else if ( i == 4 ) then
        i1mach = 6
      else if ( i == 5 ) then
        i1mach = 32
      else if ( i == 6 ) then
        i1mach = 4
      else if ( i == 7 ) then
        i1mach = 2
      else if ( i == 8 ) then
        i1mach = 31
      else if ( i == 9 ) then
        i1mach = 2147483647
      else if ( i == 10 ) then
        i1mach = 2
      else if ( i == 11 ) then
        i1mach = 24
      else if ( i == 12 ) then
        i1mach = -125
      else if ( i == 13 ) then
        i1mach = 128
      else if ( i == 14 ) then
        i1mach = 53
      else if ( i == 15 ) then
        i1mach = -1021
      else if ( i == 16 ) then
        i1mach = 1024
      else if ( 16 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
        write ( *, '(a,i12)' ) '  I = ', i
        i1mach = 0
        stop
      end if

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      function d1mach ( i )

c*********************************************************************72
c
cc D1MACH returns double precision machine dependent constants.
c
c  Discussion:
c
c    D1MACH can be used to obtain machine dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = D1MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    D1MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
c    D1MACH ( 2) = B^EMAX*(1 - B^(-T)), the largest magnitude.
c    D1MACH ( 3) = B^(-T), the smallest relative spacing.
c    D1MACH ( 4) = B^(1-T), the largest relative spacing.
c    D1MACH ( 5) = LOG10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision D1MACH, the value of the constant.
c
      implicit none

      double precision d1mach
      integer i

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      else if ( i == 1 ) then
        d1mach = 4.450147717014403D-308
      else if ( i == 2 ) then
        d1mach = 8.988465674311579D+307
      else if ( i == 3 ) then
        d1mach = 1.110223024625157D-016
      else if ( i == 4 ) then
        d1mach = 2.220446049250313D-016
      else if ( i == 5 ) then
        d1mach = 0.301029995663981D+000
      else if ( 5 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      end if

      return
      end
*
**********************************************************************
*
      subroutine UFPARM ( meqn, nvar, x )

c*********************************************************************72
c
cc UFPARM is a user-supplied external routine.
c
c  Discussion:
c
c    The name of the routine, the argument list, and even whether
c    it is a function or subroutine, are left to the user.
c
c    NL2SOL simply passes the external reference from the calling
c    program through to the residual and jacobian routines.
c
c    If the user has no need for this facility, then a dummy
c    routine like this one may be used.
c
c  Parameters:
c
c    Input, integer MEQN, the number of functions.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the current value of the variables.
c
      implicit none

      integer meqn
      integer nvar

      double precision x(nvar)

      return
      end
*
**********************************************************************
*
