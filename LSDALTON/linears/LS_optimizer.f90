!===========================================================! 
!    The main driver for geometry optimization in LSDALTON  !
!===========================================================!
! Written by Vladimir Rybkin in 11/2010
!
SUBROUTINE LS_RUNOPT(E,config,optinfo,H1,F,D,S,ls,lupri,luerr)
use precision
use ls_util 
use files
use lstiming
use configuration
use lsdalton_fock_module
use matrix_module
use dal_interface
use optimization_input
use memory_handling
use Readmolefile
Implicit Real(realk) (A-H,O-Z)
!  All these general entities needed to get energy and gradient
Type(lsitem) :: ls   ! General information,used only to get E and gradient
Type(Matrix), intent(inout) :: F,D,S ! Fock,density,overlap matrices
Type(Matrix), intent(inout) :: H1   ! One electron matrix
Type(ConfigItem), intent(inout) :: Config ! General information
!
Type(opt_setting) :: optinfo
LOGICAL ls_minend, INDXOK, STATPO, ACTIVE
LOGICAL REJGEO, TRU, FAL, TMPLOG, NEWSTP, NEWBMT
CHARACTER TMPLIN*80, WRDRSP*7
Integer :: lupri, luerr   ! File units
Real(realk) :: E   ! Energy
Real(realk), pointer ::  EGRAD(:), CSTEP(:)
Real(realk), pointer ::  GRDOLD(:), GRDMAT(:,:)
Real(realk), pointer ::  STPMAT(:,:), HESOLD(:,:)
Real(realk), pointer ::  GRDARR(:,:), STPARR(:,:)
!
!     The array geinfo contains optimization information for each
!     iteration. The first index is the iteration, the second gives
!     the property:   1  -  Energy
!                     2  -  Gradient norm
!                     3  -  Index of Hessian
!                           (a negative index indicates symmetry break)
!                     4  -  Step length
!                     5  -  Trust radius
!                     6  -  # rejected steps
!
Real(realk), pointer :: GEINFO(:,:)
Real(realk), pointer :: WILBMT(:,:), BMTRAN(:,:)
Real(realk), pointer :: HESINT(:,:), VECMOD(:)
Real(realk) :: TE,TS ! CPU time
PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
Real(realk), pointer :: KEHESS(:),KALHES(:)
Real(realk), pointer :: KBMINV(:),KPJINM(:),KEVEC(:), KCONMT(:)
Real(realk), pointer :: KTEMP1(:),KTEMP2(:),KTEMP3(:),KTEMP4(:), &
      &   KTEMP5(:),KTEMP6(:),KTEMP7(:),KTEMP8(:),KTEMP9(:)
Real(realk), pointer :: TEMP1(:,:),TEMP2(:,:)
! 
CHARACTER(len=70) :: WORD
LOGICAL :: ReadWord
Integer :: luinfo,fileStatus
Integer :: NRedint
! Quit statement for failed print level
If ((optinfo%RedInt .OR. optinfo%DelInt) .AND. optinfo%IPrint .GE. 15) then
   Call LSQuit('The print level in **OPTIMI is too high for optimization &
    redundant internals. Try to reduce it', lupri)
Endif
!
! First, grabbing number of atoms
  NAtoms = config%Molecule%nAtoms
! Second, determining numbers of various coordinates
! needed for memory allocation
  MXCENT = NAtoms
  MXCOOR = 3*NAtoms
  optinfo%ICartCoord = 3*NAtoms
  optinfo%NTempMat = 6*optinfo%ICartCoord
  optinfo%NCoordTot = optinfo%ICartCoord
  optinfo%energy = E
! Knowing number of cartesians,allocating cartesian coordinates vector
  Call mem_alloc(optinfo%Coordinates,3,MXCENT)
  Call ls_dzero(optinfo%Coordinates,3*MXCENT)
!
! Third, grabbing the initial coordinates
!
      Do i = 1, NAtoms
         optinfo%Coordinates(:,i) = config%Molecule%Atom(i)%Center(:)  
      Enddo
! Fourth, finding the number of internals if needed
  If (optinfo%RedInt .OR. optinfo%DelInt .OR. optinfo%InmdHess) then  
      ! Allocate some memory for finding number of redundant internals
      Call mem_alloc(TEMP1,NAtoms,8)
      Call mem_alloc(TEMP2,NAtoms,NAtoms)
      ! We allocate INTCRD with memory in excess, once we get
      ! the exact number of redundant internals we will allocate
      ! as much memory as needed
      Call mem_alloc(optinfo%INTCRD,NAtoms*3*8,6)
      ! Find redundant internals
      IPrint = optinfo%IPrint   ! To avoid a printout if LS_FNDRED
      optinfo%IPrint = 0
      Call LS_FNDRED(TEMP1,TEMP2,NAtoms,Config%Molecule,lupri,optinfo,.TRUE.)
      optinfo%IPrint = IPrint   ! Print level back
      ! Maximum number of coordinates set
      MXCOOR = MAX(optinfo%ICartCoord,optinfo%NIntCoord)
      ! Deallocate memory
      Call mem_dealloc(TEMP1)
      Call mem_dealloc(TEMP2)
      Call mem_dealloc(optinfo%INTCRD)
! If .Findre no optimization is carried out, but only the redundant internals are found
      If (optinfo%FindRe) then
          Call lsheader(lupri,'Redundant internals found')
          Call mem_dealloc(optinfo%Coordinates)
          Return
      Endif
  Endif
! Allocating gradient and Hessian
      Call mem_alloc(optinfo%GradMol,optinfo%ICartCoord)
      Call mem_alloc(optinfo%HessMol,optinfo%ICartCoord,optinfo%ICartCoord)
!Initializing them
      Call ls_dzero(optinfo%GradMol,optinfo%ICartCoord)
      Call ls_dzero(optinfo%HessMol,optinfo%ICartCoord*optinfo%ICartCoord)
! Allocating more
      Call mem_alloc(optinfo%STPDIA,MXCOOR)
      Call mem_alloc(optinfo%STPSYM,MXCOOR)
      Call mem_alloc(optinfo%GRDDIA,MXCOOR)
      Call mem_alloc(optinfo%EVAL,  MXCOOR)
      Call mem_alloc(optinfo%EVALOL,MXCOOR)
      Call mem_alloc(GRDMAT,MXCOOR,25)
      Call mem_alloc(STPMAT,MXCOOR,25)
      Call mem_alloc(GRDARR,MXCOOR,25)
      Call mem_alloc(STPARR,MXCOOR,25)
      Call mem_alloc(EGRAD,MXCOOR)
      Call mem_alloc(GRDOLD,MXCOOR)
      Call mem_alloc(HESOLD,MXCOOR,MXCOOR)
  If (optinfo%RedInt .OR. optinfo%DelInt .OR. optinfo%InmdHess) then 
      Call mem_alloc(optinfo%GRDINT,MXCOOR)
      Call mem_alloc(optinfo%STPINT,MXCOOR)
      Call mem_alloc(optinfo%CoordInt,MXCOOR)
      Call mem_alloc(optinfo%INTCRD,8*MXCOOR,6)
  Endif
  Call mem_alloc(HESINT,MXCOOR,MXCOOR)
  Call mem_alloc(BMTRAN,MXCOOR,MXCOOR)
  Call mem_alloc(KBMINV,optinfo%NIntCoord*MXCOOR)
  Call mem_alloc(WILBMT,optinfo%NIntCoord,MXCOOR)
  Call mem_alloc(VECMOD,MXCOOR)
!
!     Allocate GEINFO 
!
      call mem_alloc(GEINFO,optinfo%MaxIter,6,.TRUE.,.FALSE.)
      call LSTIMER('START ',TS,TE,lupri)
!
!     Initialization of variables.
!
      THRLDP = 1.0D-4
      THRIND = 5.0D-4
      TOLST  = 1.0D-5
      call ls_DZERO(GEINFO,(optinfo%MaxIter)*6)
      call ls_DZERO(GRDARR,25*optinfo%NIntCoord)
      call ls_DZERO(STPARR,25*optinfo%NIntCoord)
      ACTIVE = .FALSE.
      optinfo%KEPTIT = 0
      GEINFO(0,5) = optinfo%TrustRad
      IWOFF = 0
!
!     Perform preoptimization if requested.
!
      IF (optinfo%DoPre) call lsquit('DOPRE not an option in lsdalton',lupri)
!
!     Allocating some memory
!
      call mem_alloc(KEHESS,MXCOOR*MXCOOR)
      call mem_alloc(KALHES,MXCOOR*MXCOOR)
!
!     Calculate gradient and Hessian for second order method and
!     first order method with initial Hessian.
!
If (.NOT. optinfo%Findre ) then
!
!     First order methods only require the energy and the gradient.
!
         call Get_Gradient(lupri,NAtoms,F,D,ls,optinfo)
Endif
!
!     Make optinfo%VRLM-file of initial geometry if requested.
!
      IF (optinfo%VRLM) CALL LSQUIT('No VRLM implemented!',LUPRI)
!
!     Save initial geometry and energy to MOLDEN file
!
!      IF (MOLDEN) call LSQUIT('No MOLDEN implemented!',LUPRI)

!
!     We allocate more 
!
      Call mem_alloc(CSTEP,MXCOOR)
      Call mem_alloc(KPJINM,optinfo%NIntCoord*optinfo%NIntCoord)
      Call mem_alloc(KEVEC,MXCOOR*MXCOOR)
      Call mem_alloc(KCONMT,optinfo%NIntCoord*optinfo%NIntCoord)
      Call mem_alloc(KTEMP1,MXCOOR*MXCOOR)
      Call mem_alloc(KTEMP2,MXCOOR*MXCOOR)
      Call mem_alloc(KTEMP3,MXCOOR*MXCOOR)
      Call mem_alloc(KTEMP4,MXCOOR*MXCOOR)
      Call mem_alloc(KTEMP5,MXCOOR*MXCOOR)
      Call mem_alloc(KTEMP6,MXCOOR*MXCOOR)
      Call mem_alloc(KTEMP7,MXCOOR*MXCOOR)
      Call mem_alloc(KTEMP8,MXCOOR*MXCOOR)
      Call mem_alloc(KTEMP9,MXCOOR*MXCOOR)
!     Set cartesian step equal to zero
      call ls_dzero(CSTEP,MXCOOR)
!
!     Check if redundant internal coordinates should be used.
!
      IF (optinfo%DelInt .OR. optinfo%RedInt .OR. optinfo%InrdHess .OR. optinfo%InmdHess) THEN
         call ls_INIRED(optinfo%NIntCoord,MXCOOR,WILBMT,BMTRAN,KBMINV, &
     &        KPJINM,KTEMP1,KTEMP2,KTEMP3,KTEMP4,KTEMP5,KTEMP6, &
              Config%Molecule,NAtoms,lupri,optinfo)
      END IF
      IF (optinfo%DelInt .OR. optinfo%RedInt) THEN
         NCRDHS = optinfo%NIntCoord
      ELSE
         NCRDHS = optinfo%ICartCoord
      END IF
      IF (optinfo%RatFun) NCRDHS = NCRDHS + 1
!
!     Initialize Hessian if first order method is used.
!
 7    CONTINUE
      IF (.NOT. optinfo%Newton) THEN
        call ls_INIHES(config%Molecule, &
         &  MXCOOR,MXCOOR,GRDOLD,HESOLD, &
         &  KTEMP1,KTEMP2,KTEMP3,KTEMP4,WILBMT,BMTRAN, & 
         &  KBMINV,HESINT,lupri,optinfo)
      ENDIF
!
      call ls_DZERO(EGRAD,MXCOOR)
      call ls_DZERO(KEHESS,MXCOOR*MXCOOR)
      call ls_DZERO(KALHES,MXCOOR*MXCOOR)
      DO i = 1, optinfo%ICartCoord
         EGRAD(i) = optinfo%GradMol(i)
      ENDDO  
      JI = 1
      DO I = 1, optinfo%ICartCoord
         DO J = 1, optinfo%ICartCoord
            KEHESS(JI) = optinfo%HessMol(J,I) 
            KALHES(JI) = optinfo%HessMol(J,I) 
            JI = JI + 1
         ENDDO
      ENDDO
!      Call ls_copyGH(EGRAD,KEHESS,KALHES,optinfo%GradMol,optinfo%HessMol,optinfo%ICartCoord)
!     Construct projection operator and use it.
!     Then diagonalize Hessian.
!
      IF (optinfo%RedInt .OR. optinfo%DelInt) THEN

         IF (optinfo%Newton) call ls_CGHINT(config%Molecule, &
     &        optinfo%NIntCoord,MXCOOR,KTEMP1, &
     &        KTEMP2,KTEMP3,KTEMP4,KTEMP5, &
     &        WILBMT,KBMINV,BMTRAN,HESINT,lupri,optinfo)
         call ls_PRJINT(optinfo%NIntCoord,optinfo%NIntCoord,KPJINM,KCONMT, &
     &        HESINT,KTEMP1,KTEMP2,KTEMP3,KTEMP4,lupri,optinfo)
!
!     Note that the contents of KTEMP7 is passed on
!     from LINSRC to FNSTIN below.
!
         IF (optinfo%LnSearch .AND. optinfo%RatFun .AND. (optinfo%ItrNmr .GT. 0)) &
     &        call ls_LINSRC(optinfo%NIntCoord,optinfo%NIntCoord,optinfo%GRDINT,GRDARR(1,1), &
     &        KTEMP7,STPARR(1,1),KTEMP3,KTEMP4, &
     &        ACTIVE,EMOD,lupri,optinfo)
         IF (optinfo%RatFun .AND. optinfo%Saddle) NCRDHS = NCRDHS - 1
         call ls_DIAINT(optinfo%NIntCoord,MXCOOR,NCRDHS,KEVEC,KTEMP1, &
     &        KTEMP2,KTEMP3,KTEMP4,THRIND,HESINT,KTEMP5,lupri,optinfo)
         IF (optinfo%RatFun .AND. optinfo%Saddle) NCRDHS = NCRDHS + 1
      ELSE
!
!     Note that the contents of KTEMP1 is passed on
!     from PROJGH to DIAHES below.
!
         call ls_PROJGH(EGRAD,KEHESS,KALHES,KTEMP1, &
     &   KTEMP2,KTEMP3,KTEMP4,optinfo%NCoordTot, &
     &   optinfo%NTempMat,lupri,optinfo)
         IF (optinfo%LnSearch .AND. optinfo%RatFun .AND. (optinfo%ItrNmr .GT. 0)) &
     &      call ls_LINSRC(optinfo%ICartCoord,MXCOOR,EGRAD,GRDARR(1,1),CSTEP, &
     &        STPARR(1,1),KTEMP3,KTEMP4,ACTIVE,EMOD,lupri,optinfo)
            call ls_DIAHES(optinfo%NIntCoord,MXCOOR,NCRDHS,EGRAD,KEHESS,KALHES, &
              & KTEMP1,THRIND, &
              & KEVEC,KTEMP2,KTEMP3,KTEMP4,optinfo%NCoordTot,optinfo%NTempMat, &
              & lupri,optinfo)
      END IF
      GEINFO(0,1) = E
      GEINFO(0,3) = optinfo%IndHes*1.0D0
!!!!!! Vladimir:: currently disabled
!     Write Hessian to file (for 1st order restarts).
!
!      IF (.NOT. optinfo%NoHessianWrite) &
!     &     call ls_PNCHES(optinfo%NIntCoord,MXCOOR,HESINT,WILBMT,BMTRAN,KTEMP1, &
!     &     KTEMP2,KTEMP3,KTEMP4,WORK(KWRK2),LWRK2,lupri,optinfo)

!
!     Determine step, check for convergence, print output and
!     and update geometry.
!
      IREJ = 0
 755  CONTINUE
      IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
         call ls_FNSTIN(config%Molecule,optinfo%NIntCoord, &
     &        MXCOOR,NCRDHS,HESINT,KEVEC, &
     &        KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
     &        KTEMP5,CSTEP,WILBMT,BMTRAN,KBMINV,GRDARR, &
     &        STPARR,ACTIVE,EMOD,VECMOD,KTEMP7,lupri,optinfo)
      ELSE
         call ls_FNDSTP(MXCOOR,MXCOOR,NCRDHS,EGRAD,KEHESS, &
     &        KEVEC,KTEMP1,KTEMP2,KTEMP3, &
     &        KTEMP4,KTEMP5,CSTEP,GRDARR,STPARR, &
     &        ACTIVE,EMOD,VECMOD,lupri,optinfo)
      END IF
      optinfo%GeConv = ls_minend(optinfo%NIntCoord,BMTRAN,KTEMP1,KTEMP2,lupri,optinfo)
!
!     If there has been a completely failed step, the geometry has
!     by default not converged.
!
      IF (ABS(GEINFO(0,6)) .GT. 1.0D-3) optinfo%GeConv = .FALSE.
      call ls_PRIALL(config%Molecule,NAtoms,CSTEP,KTEMP1,lupri,optinfo)
!
!     Save this geometry and energy to MOLDEN file
!
      NEWSTP = .FALSE.
      NEWBMT = .FALSE.
!
!     To allow reinitialization
!
      optinfo%InitHess = .FALSE.
!
      IF (.NOT. optinfo%GeConv) then 
           call Find_Geometry(E,CSTEP,EGRAD,KTEMP1, &
     &     KTEMP2,IREJ,GEINFO,NEWSTP,NEWBMT,NAtoms,lupri, &
           luerr,config,ls,H1,F,D,S,optinfo)
      !    Grabbing coordinates
           Do i = 1, NAtoms
              optinfo%Coordinates(:,i) = KTEMP1(3*i-2:3*i)
           Enddo
      !    Renovating Config%Molecule
           Do i = 1,NAtoms
              Config%Molecule%Atom(i)%Center(:)=optinfo%Coordinates(:,i)
           Enddo
      ENDIF
      !
      IF (NEWSTP) GOTO 755
      GEINFO(0,2) = optinfo%GradNorm
      GEINFO(0,4) = optinfo%StepNorm
      IF (optinfo%ItrNmr .LT. optinfo%IterMax) GEINFO(1,5) = optinfo%TrustRad
      IF (ABS(GEINFO(0,6)) .LT. 1.0D-3) THEN
         GEINFO(0,6) = IREJ*1.0D0
      ELSE
         GEINFO(0,6) = -(ABS(GEINFO(0,6))+ABS(IREJ)*1.0D0)
      END IF
      optinfo%TotRj = optinfo%TotRj + ABS(IREJ)
!
!     Determine value of the various coordinates
!
      IF (optinfo%RedInt .AND. (optinfo%IPrint .GE. 1)) THEN
         call ATOM_INI(KTEMP1,Config%Molecule,optinfo,NAtoms,.TRUE.,lupri)
         call ls_GETINT(NAtoms,optinfo%NIntCoord,KTEMP1,optinfo%CoordInt,lupri,optinfo)
         call lsheader(lupri,'New internal coordinates') 
         call output(optinfo%CoordInt,1,1,1,optinfo%NIntCoord,1,optinfo%NIntCoord,1,LUPRI)
         WRITE(LUPRI,'(//)')
      END IF
!
!     If the step has failed
!
      IF (IREJ .LT. 0) THEN
         GOTO 7
!      ELSE IF (optinfo%RejIni .AND. optinfo%RedInt .AND. (optinfo%TotRj .GE. 3)) THEN
!         WRITE(LUPRI,*)'***** NOTE! *****'
!         WRITE(LUPRI,*)
!     &        'The number of dihedral angles will be reduced!'
!         call ls_RREDUN
!         optinfo%TotRj = 0
      END IF
!
      IF (optinfo%GeConv .AND. optinfo%DoPre .AND. (.NOT. optinfo%FinPre)) THEN
         optinfo%KeepHessian = .FALSE.
      END IF
!
!     Switch the iteration counter
!
!      optinfo%ItrNmr = optinfo%ItrNmr + 1
!
!     DO WHILE-loop that runs until geometry has converged or
!     maximum number of iterations is reached.
!
 10   CONTINUE
      call LSTIMER('Geom. opt.',TS,TE,lupri)

      IF ((optinfo%ItrNmr .LT. optinfo%IterMax) .AND. (.NOT. optinfo%GeConv)) THEN
         optinfo%ItrNmr = optinfo%ItrNmr + 1
         NCRD = optinfo%NCoordTot
         IF (optinfo%RedInt .OR. optinfo%DelInt) NCRD = optinfo%NIntCoord
         DO 20 I = 1, NCRD
            optinfo%EVALOL(I) = optinfo%EVAL(I)
 20      CONTINUE
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call ls_UPGDST(optinfo%NIntCoord,optinfo%NIntCoord,GRDARR,STPARR, &
            & optinfo%GRDINT,optinfo%STPINT,lupri,optinfo)
         ELSE
            call ls_UPGDST(optinfo%ICartCoord,MXCOOR,GRDARR,STPARR,& 
            & EGRAD,optinfo%STPSYM,lupri,optinfo)
         END IF
!
!     We go through the same procedure as for the first iteration,
!     but here we call Get_Energy. We called the optimizer after 
!     calculating the first energy in LSDALTON

!!!!! Vladimir:: No Newton yet!     
!         IF (optinfo%Newton) THEN
!            call ls_GTHESS(EGRAD,KEHESS,KALHES, &
!     &           EXHER,EXSIR,EXABA,WORK(KWRK1),LWRK1,IWOFF, &
!     &           WRKDLM,optinfo%NCoordTot,optinfo)
!         ELSE
!            call Get_Energy(E,config,optinfo,H1,F,D,S,ls,NAtoms,lupri,luerr)
            call Get_Gradient(lupri,NAtoms,F,D,ls,optinfo)
!         END IF
!
!     If redundant internal coordinates are used, Wilson's B matrix,
!     its derivative and its inverse must be updated.
!
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call ls_GETWIL(config%Molecule,optinfo%NIntCoord,MXCOOR, &
     &           KTEMP1,WILBMT, &
     &           BMTRAN,KTEMP2,lupri,NAtoms,optinfo)
            IF (optinfo%IPrint .GE. IPRMAX) &
     &           call ls_GETDWL(config%Molecule,optinfo%NIntCoord, &
     &                KTEMP1,KTEMP2, &
     &                KTEMP3,WILBMT,lupri,NAtoms,optinfo)
            call ls_GTBINV(optinfo%NIntCoord,KTEMP1,KTEMP2,KTEMP3, &
     &           KTEMP4,WILBMT,BMTRAN,KBMINV,KPJINM, &
     &           KTEMP5,KTEMP6,lupri,optinfo)
         END IF
!
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            NCRDHS = optinfo%NIntCoord
         ELSE
            NCRDHS = optinfo%ICartCoord
         END IF
         IF (optinfo%RatFun) NCRDHS = NCRDHS + 1
         IF (.NOT. optinfo%Newton) THEN
            IF (optinfo%Rebid) THEN
               call ls_INIHES(config%Molecule,MXCOOR, &
     &              MXCOOR,GRDOLD,HESOLD, &
     &              KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
     &              WILBMT,BMTRAN,KBMINV,HESINT,lupri,optinfo)
               optinfo%Rebid = .FALSE.
            ELSE
               call ls_UPDHES(config%Molecule,MXCOOR,MXCOOR, &
     &              GRDOLD,GRDMAT,STPMAT, &
     &              HESOLD,KTEMP1,KTEMP2,KTEMP3,KTEMP4,KTEMP5,KTEMP6, &
     &              KTEMP7,KTEMP8,KTEMP9,WILBMT,BMTRAN,KBMINV, &
     &              HESINT,NINT(ABS(GEINFO(optinfo%ItrNmr-1,6))), &
     &              NINT(ABS(GEINFO(optinfo%ItrNmr,6))),lupri,optinfo)
            END IF
!
      call ls_DZERO(EGRAD,MXCOOR)
      call ls_DZERO(KEHESS,MXCOOR*MXCOOR)
      call ls_DZERO(KALHES,MXCOOR*MXCOOR)
            DO i = 1, optinfo%ICartCoord
               EGRAD(i) = optinfo%GradMol(i)
            ENDDO  
            JI = 1
            DO I = 1, optinfo%ICartCoord
               DO J = 1, optinfo%ICartCoord
                  KEHESS(JI) = optinfo%HessMol(J,I) 
                  KALHES(JI) = optinfo%HessMol(J,I) 
                  JI = JI + 1
                ENDDO
            ENDDO
         END IF
 33      CONTINUE
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            IF (optinfo%Newton) call ls_CGHINT(config%Molecule, &
     &           optinfo%NIntCoord,MXCOOR,KTEMP1, &
     &           KTEMP2,KTEMP3,KTEMP4,KTEMP5, &
     &           WILBMT,KBMINV,BMTRAN,HESINT,lupri,optinfo)
            call ls_PRJINT(optinfo%NIntCoord,optinfo%NIntCoord,KPJINM,KCONMT, &
     &           HESINT,KTEMP1,KTEMP2,KTEMP3,KTEMP4,lupri,optinfo)
            IF (optinfo%LnSearch .AND. optinfo%RatFun .AND. (optinfo%ItrNmr .GT. 0)) &
     &           call ls_LINSRC(optinfo%NIntCoord,optinfo%NIntCoord,optinfo%GRDINT,GRDARR(1,1), &
     &           KTEMP7,STPARR(1,1),KTEMP3,KTEMP4, &
     &           ACTIVE,EMOD,lupri,optinfo)
            IF (optinfo%RatFun .AND. optinfo%Saddle) NCRDHS = NCRDHS - 1
            call ls_DIAINT(optinfo%NIntCoord,MXCOOR,NCRDHS,KEVEC,KTEMP1, &
     &           KTEMP2,KTEMP3,KTEMP4,THRIND,HESINT, &
     &           KTEMP5,lupri,optinfo)
            IF (optinfo%RatFun .AND. optinfo%Saddle) NCRDHS = NCRDHS + 1
         ELSE
            call ls_PROJGH(EGRAD,KEHESS,KALHES,KTEMP1, &
     &      KTEMP2,KTEMP3,KTEMP4,optinfo%NCoordTot, &
     &      optinfo%NTempMat,lupri,optinfo)
            IF (optinfo%LnSearch .AND. optinfo%RatFun .AND. (optinfo%ItrNmr .GT. 0)) &
     &           call ls_LINSRC(optinfo%ICartCoord,MXCOOR,EGRAD,GRDARR(1,1),CSTEP, &
     &           STPARR(1,1),KTEMP3,KTEMP4,ACTIVE,EMOD,lupri,optinfo)
            call ls_DIAHES(optinfo%NIntCoord,MXCOOR,NCRDHS,EGRAD,KEHESS, &
     &           KALHES,KTEMP1,THRIND,KEVEC,KTEMP2,KTEMP3,KTEMP4, &
     &           optinfo%NCoordTot,optinfo%NTempMat,lupri,optinfo)
         END IF
!
!     Update information for this iteration
!
         GEINFO(optinfo%ItrNmr,3) = optinfo%IndHes*1.0D0
         IREJ = 0
!!!!!! Vladimir:: currently disabled
!     Write Hessian to file
!
!         IF (.NOT. optinfo%NoHessianWrite) &
!     &      call ls_PNCHES(optinfo%NIntCoord,MXCOOR,HESINT,WILBMT,BMTRAN,KTEMP1, &
!     &      KTEMP2,KTEMP3,KTEMP4,WORK(KWRK2),LWRK2, &
!     &      lupri,optinfo)

!
 756     CONTINUE
         call ls_dzero(CSTEP,optinfo%ICartCoord)
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call ls_FNSTIN(config%Molecule,optinfo%NIntCoord, &
     &           MXCOOR,NCRDHS,HESINT,KEVEC, &
     &           KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
     &           KTEMP5,CSTEP,WILBMT,BMTRAN,KBMINV,GRDARR, &
     &           STPARR,ACTIVE,EMOD,VECMOD,KTEMP7,lupri,optinfo)
         ELSE
            call ls_FNDSTP(MXCOOR,MXCOOR,NCRDHS,EGRAD,KEHESS, &
     &           KEVEC,KTEMP1,KTEMP2,KTEMP3, &
     &           KTEMP4,KTEMP5,CSTEP,GRDARR,STPARR, &
     &           ACTIVE,EMOD,VECMOD,lupri,optinfo)
         END IF
         optinfo%GeConv = ls_minend(optinfo%NIntCoord,BMTRAN,KTEMP1,KTEMP2,lupri,optinfo)
         IF (ABS(GEINFO(optinfo%ItrNmr,6)) .GT. 1.0D-3) optinfo%GeConv = .FALSE.
         call ls_PRIALL(Config%Molecule,NAtoms,CSTEP,KTEMP1,lupri,optinfo)
!
!
         IF (optinfo%RedInt .AND. (optinfo%IPrint .GE. 1)) THEN
            call ATOM_INI(KTEMP1,Config%Molecule,optinfo,NAtoms,.TRUE.,lupri)
            call ls_GETINT(NAtoms,optinfo%NIntCoord,KTEMP1,optinfo%CoordInt,lupri,optinfo)
            call lsheader(lupri,'New internal coordinates')
            call output(optinfo%CoordInt,1,1,1,optinfo%NIntCoord,1,optinfo%NIntCoord,1,LUPRI)
            WRITE(LUPRI,'(//)')
         END IF
         NEWSTP = .FALSE.
         IF (.NOT. optinfo%GeConv) then 
              call Find_Geometry(E,CSTEP,EGRAD,KTEMP1,KTEMP2,&
     &        IREJ,GEINFO,NEWSTP,NEWBMT,NAtoms,lupri,luerr, &
     &        config,ls,H1,F,D,S,optinfo)
              !  Grabbing coordinates
              Do i = 1, NAtoms
                 optinfo%Coordinates(:,i) = KTEMP1(3*i-2:3*i)
              Enddo
              !    Renovating Config%Molecule
              Do i = 1,NAtoms
                 Config%Molecule%Atom(i)%Center(:)=optinfo%Coordinates(:,i)
              Enddo
         ENDIF 
!            
         IF (NEWSTP) GOTO 756
         GEINFO(optinfo%ItrNmr,2) = optinfo%GradNorm
         GEINFO(optinfo%ItrNmr,4) = optinfo%StepNorm
         IF (optinfo%ItrNmr .LT. optinfo%IterMax) GEINFO(optinfo%ItrNmr+1,5) = optinfo%TrustRad
         IF (ABS(GEINFO(optinfo%ItrNmr,6)) .LT. 1.0D-3) THEN
            GEINFO(optinfo%ItrNmr,6) = IREJ*1.0D0
         ELSE
            GEINFO(optinfo%ItrNmr,6) = -(ABS(GEINFO(optinfo%ItrNmr,6))+ABS(IREJ)*1.0D0)
         END IF
         optinfo%TotRj = optinfo%TotRj + ABS(IREJ)
         IF (optinfo%Rebid) THEN
            call dcopy(optinfo%NIntCoord,optinfo%STPINT,1,KTEMP7,1)
            call ls_DZERO(optinfo%STPINT,optinfo%NIntCoord)
            call ls_DZERO(GRDOLD,optinfo%NIntCoord)
            DO 605 I = 1, optinfo%NIntCoord
               DO 607 J = 1, optinfo%NIntCoord
                  optinfo%STPINT(I) = optinfo%STPINT(I) + BMTRAN(I,J)*KTEMP7(J-1)
                  GRDOLD(I) = GRDOLD(I) + BMTRAN(I,J)*optinfo%GRDINT(J)
 607           CONTINUE
 605        CONTINUE
         END IF
!
!     If the step has failed
!
         IF (IREJ .LT. 0) THEN
            IF (.NOT. optinfo%Newton) THEN
               IF (NEWBMT) THEN
                  NCRDHS = optinfo%NIntCoord
                  IF (optinfo%RatFun) NCRDHS = NCRDHS + 1
                  call ls_GETWIL(config%Molecule,optinfo%NIntCoord, &
     &                 MXCOOR,KTEMP1,WILBMT, &
     &                 BMTRAN,KTEMP2,lupri,NAtoms,optinfo)
                  call ls_GTBINV(optinfo%NIntCoord,KTEMP1,KTEMP2, &
     &                 KTEMP3,KTEMP4,WILBMT,BMTRAN, &
     &                 KBMINV,KPJINM,KTEMP5,KTEMP6,lupri,optinfo)
               END IF
               call ls_INIHES(config%Molecule,MXCOOR, &
     &              MXCOOR,GRDOLD,HESOLD, &
     &              KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
     &              WILBMT,BMTRAN,KBMINV,HESINT,lupri,optinfo)
!
               call ls_DZERO(EGRAD,MXCOOR)
               call ls_DZERO(KEHESS,MXCOOR*MXCOOR)
               call ls_DZERO(KALHES,MXCOOR*MXCOOR)
               DO i = 1, optinfo%ICartCoord
                  EGRAD(i) = optinfo%GradMol(i)
               ENDDO  
               JI = 1
               DO I = 1, optinfo%ICartCoord
                  DO J = 1, optinfo%ICartCoord
                     KEHESS(JI) = optinfo%HessMol(J,I) 
                     KALHES(JI) = optinfo%HessMol(J,I) 
                     JI = JI + 1
                  ENDDO
               ENDDO
!      Call ls_copyGH(EGRAD,KEHESS,KALHES,optinfo%GradMol,optinfo%HessMol,optinfo%ICartCoord)
            END IF
            GOTO 33
!         ELSE IF (optinfo%RejIni .AND. optinfo%RedInt .AND. (optinfo%TotRj .GE. 5)) THEN
!            WRITE(LUPRI,*)'***** NOTE! *****'
!            WRITE(LUPRI,*)
!     &           'The number of dihedral angles will be reduced!'
!            call ls_RREDUN
!            optinfo%TotRj = 0
         END IF
!
!     Check if rejected steps should cause reinitialization of Hessian.
!
         IF ((.NOT. optinfo%Newton) .AND. (optinfo%RejIni .AND. (IREJ .GE. 1))) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)'***** NOTE! *****'
            WRITE(LUPRI,*) &
     &           'Due to rejected step, Hessian is reinitialized.'
            WRITE(LUPRI,*)
            call ls_INIHES(config%Molecule,MXCOOR, &
     &           MXCOOR,GRDOLD,HESOLD, &
     &           KTEMP1,KTEMP2,KTEMP3,KTEMP4,WILBMT,BMTRAN, &
     &           KBMINV,HESINT,lupri,optinfo)
            optinfo%TrustRad = GEINFO(0,5)
            GEINFO(optinfo%ItrNmr+1,5) = optinfo%TrustRad
            optinfo%Restart = .TRUE.
         END IF
!
!     Check if increase of gradient norm should cause reinitialization
!     of Hessian. Reinitialization occurs when the norm of the gradient
!     is larger than the norm of the gradient two iterations earlier.
!
         IF (.NOT.optinfo%Newton .AND. optinfo%GradIni .AND. (optinfo%ItrNmr .GE. 2)) THEN
            IF (GEINFO(optinfo%ItrNmr,2) .GE. GEINFO(optinfo%ItrNmr-2,2)) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*)'***** NOTE! *****'
               WRITE(LUPRI,*)'Due to increasing gradient norm,   &
     &              Hessian is reinitialized.'
               WRITE(LUPRI,*)
               call ls_INIHES(config%Molecule,MXCOOR, &
     &              MXCOOR,GRDOLD,HESOLD, &
     &              KTEMP1,KTEMP2,KTEMP3,KTEMP4,WILBMT,BMTRAN, &
     &              KBMINV,HESINT,lupri,optinfo)
               optinfo%TrustRad = GEINFO(0,5)
               GEINFO(optinfo%ItrNmr+1,5) = optinfo%TrustRad
            END IF
         END IF
!
         IF (optinfo%GeConv .AND. optinfo%DoPre .AND. (.NOT. optinfo%FinPre)) THEN
            optinfo%KeepHessian = .FALSE.
         END IF
!
         GOTO 10
!
!     Finished case 1: Geometry has converged.
!
      ELSE IF (optinfo%GeConv) THEN
!
!     Final results are printed, partially through PRIINF.
!
         call lsheader(lupri,' End of Optimization ')
         call ls_PRIINF(Config%Molecule,NAtoms,GEINFO,lupri,optinfo)
         WRITE(LUPRI,*)
         IF (optinfo%ConOpt) THEN
            WRITE(LUPRI,*) 'Constrained optimization converged in ', &
     &           optinfo%ItrNmr+1, ' iterations!'
            IF (optinfo%GradNorm .GT. optinfo%GradThr) THEN
               WRITE(LUPRI,*) 'Removing the  &
     &            constraint(s) might decrease the energy further.'
            ELSE
               WRITE(LUPRI,*) 'A saddle point might have been reached.'
            END IF
         ELSE
            WRITE(LUPRI,*) 'Geometry converged in ', optinfo%ItrNmr+1, &
     &           ' iterations!'
         END IF
         IF (optinfo%Newton .AND. optinfo%Saddle .AND. (optinfo%IndHes .NE. 1)) THEN
            WRITE(LUPRI,'(/A/A)') &
     &         ' Please note that Hessian index does not correspond', &
     &         ' to a first order saddle point (transition state).'
         END IF
         E = GEINFO(optinfo%ItrNmr,1)
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A,F14.6,A)') &
     &          ' Energy at final geometry is       : ',E,' a.u.'
         ERGDIF = optinfo%energy - GEINFO (0,1)
         WRITE(LUPRI,'(A,F14.6,A)') &
     &          ' Energy change during optimization : ',ERGDIF,' a.u.'
         ERGDIF = ERGDIF * XKJMOL
         WRITE(LUPRI,'(A,F14.3,A)') &
     &        '                                     ',ERGDIF,' kJ/mol'
         WRITE(LUPRI,*)
         IF (optinfo%DoPre) THEN
            WRITE(LUPRI,'(A)') ' Preoptimization was performed using'// &
     &           ' the basis set(s):'
            DO 111 I = 1, optinfo%Pre-1
               WRITE(LUPRI,'(A,A60)') '     ',optinfo%PreText(I)
 111        CONTINUE
            WRITE(LUPRI,*)
         END IF
         IF (optinfo%DoSpE) THEN
            E = GEINFO(optinfo%ItrNmr+1,1)
            WRITE(LUPRI,'(A,A60)') ' Using the basis ',optinfo%SpBText
            WRITE(LUPRI,'(A,F14.6,A)') &
     &          ' single point energy was calculated: ',E,' a.u.'
            WRITE(LUPRI,*)
         END IF
!
!     Finished case 2: Exceeded maximum number of iterations.
!
      ELSE
!     No single point energy has been calculated.
         TMPLOG = optinfo%DoSpE
         optinfo%DoSpE = .FALSE.
         call lsheader(lupri,'Optimization Control Center')
         call ls_PRIINF(Config%Molecule,NAtoms,GEINFO,lupri,optinfo)
         optinfo%DoSpE = TMPLOG
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Geometry has NOT converged!'
         WRITE(LUPRI,*) 'Maximum number of iterations (', optinfo%IterMax, &
     &                                        ') has been reached and'
         WRITE(LUPRI,*) 'optimization halted. Increase number or ', &
     &                                   'restart from last geometry.'
         IF (optinfo%DoSpE) WRITE(LUPRI,*) 'No single point energy has been  &
     &                          calculated.'
         WRITE(LUPRI,*)
      END IF
!
!  Deallocate everything
!      
      Call mem_dealloc(GEINFO)
      Call mem_dealloc(optinfo%IConstr)
      Call mem_dealloc(KEHESS)
      Call mem_dealloc(KALHES)
      Call mem_dealloc(CSTEP)
      Call mem_dealloc(KPJINM)
      Call mem_dealloc(KEVEC)
      Call mem_dealloc(KCONMT)
      Call mem_dealloc(KTEMP1)
      Call mem_dealloc(KTEMP2)
      Call mem_dealloc(KTEMP3)
      Call mem_dealloc(KTEMP4)
      Call mem_dealloc(KTEMP5)
      Call mem_dealloc(KTEMP6)
      Call mem_dealloc(KTEMP7)
      Call mem_dealloc(KTEMP8)
      Call mem_dealloc(KTEMP9)
      Call mem_dealloc(optinfo%Coordinates)
      Call mem_dealloc(optinfo%GradMol)
      Call mem_dealloc(optinfo%HessMol)
      Call mem_dealloc(optinfo%STPDIA)
      Call mem_dealloc(optinfo%STPSYM)
      Call mem_dealloc(optinfo%GRDDIA)
      Call mem_dealloc(optinfo%EVAL)
      Call mem_dealloc(optinfo%EVALOL)
      Call mem_dealloc(GRDMAT)
      Call mem_dealloc(STPMAT)
      Call mem_dealloc(GRDARR)
      Call mem_dealloc(EGRAD)
      Call mem_dealloc(GRDOLD)
      Call mem_dealloc(STPARR)
      Call mem_dealloc(HESOLD)
  If (optinfo%RedInt .OR. optinfo%DelInt .OR. optinfo%InmdHess) then 
      Call mem_dealloc(optinfo%GRDINT)
      Call mem_dealloc(optinfo%STPINT)
      Call mem_dealloc(optinfo%CoordInt)
      Call mem_dealloc(optinfo%INTCRD)
  Endif
  Call mem_dealloc(HESINT)
  Call mem_dealloc(BMTRAN)
  Call mem_dealloc(KBMINV)
  Call mem_dealloc(WILBMT)
  Call mem_dealloc(VECMOD)
  IF (optinfo%NFreeze.GT.0) call mem_dealloc(optinfo%FreezeArray)
  IF (optinfo%NAdd.GT.0) call mem_dealloc(optinfo%AddCoordArray)
  IF (optinfo%NumPre.GT.0) call mem_dealloc(optinfo%PreText)
!
      RETURN
      END
!=========================!
! Find_Geometry           !  
!=========================!
SUBROUTINE Find_Geometry(E,CSTEP,EGRAD,COONEW,COOOLD, &
     &     IREJ,GEINFO,NEWSTP,NEWBMT,NAtoms,lupri,luerr, &
           config,ls,H1,F,D,S,optinfo)
use ls_util 
use optimization_input
use files
use molecule_module
use configuration
!
!     If the step is acceptable, the geometry is updated
!     and written to file.
!
Implicit Real(realk) (A-H,O-Z)
!
Type(ConfigItem), intent(inout) :: Config ! General information
Type(lsitem) :: ls   ! General information,used only to get E and gradient
Type(Matrix), intent(inout) :: F,D,S ! Fock,density,overlap matrices
Type(Matrix), intent(inout) :: H1   ! One electron matrix
!
Integer :: lupri,luerr,NAtoms
TYPE(opt_setting) :: optinfo
Real(realk) :: CSTEP(MXCOOR), EGRAD(MXCOOR)
Real(realk) :: COONEW(3,MXCENT), COOOLD(3,MXCENT)
Integer     :: ICRD(3)
Real(realk) :: GEINFO(0:optinfo%IterMax,6)
Real(realk) :: E
CHARACTER*10 FILENM
LOGICAL REJGEO,NEWSTP,NEWBMT
LOGICAL FAILED
SAVE FAILED, IFAILD
DATA FAILED, IFAILD /.FALSE.,0/
!
      NEWSTP = .FALSE.
      REJGEO = .TRUE.
      optinfo%energyOld = optinfo%energy
      IJ = 1
      DO 10 J = 1, NAtoms
         DO 20 I = 1, 3
            COOOLD(I,J) = optinfo%Coordinates(I,J)
            COONEW(I,J) = optinfo%Coordinates(I,J) + CSTEP(IJ)
            IJ = IJ + 1
 20      CONTINUE
 10   CONTINUE
!
!     Here we start a loop to obtain acceptable step, note that REJGEO
!     is initially set TRUE to enter the loop.
!
 50   CONTINUE
      IF ((IREJ .LE. optinfo%MaxRej) .AND. REJGEO) THEN
!
!     If geometry stabilization has been requested (experimental feature!!!),
!     all coordinates with a difference less than the limit are set equal.
!
         IF (optinfo%Stblz .GT. 0) THEN
            IF (optinfo%IPrint .GT. 6) THEN
               call lsheader(lupri,'Non-stabilized geometry')
               Do i = 1,NAtoms
                  Config%Molecule%Atom(i)%Center(:)=COONEW(:,i)
               Enddo
               call Print_Geometry(config%Molecule,lupri)
            END IF
            DO 33 J = 1, nAtoms - 1
               ICRD(1) = NINT(COONEW(1,J)*10**optinfo%Stblz)
               ICRD(2) = NINT(COONEW(2,J)*10**optinfo%Stblz)
               ICRD(3) = NINT(COONEW(3,J)*10**optinfo%Stblz)
               DO 34 JJ = J, nAtoms
                  DO 35 I = 1, 3
                     IF (ABS(ICRD(I) - NINT(COONEW(I,JJ)*10**optinfo%Stblz)) &
     &                    .LE. 1) COONEW(I,JJ) = COONEW(I,J)
 35               CONTINUE
 34            CONTINUE
 33         CONTINUE
            IF (optinfo%IPrint .GT. 6) THEN
               call lsheader(lupri,'Stabilized geometry')
               Do i = 1,NAtoms
                  Config%Molecule%Atom(i)%Center(:)=COONEW(:,i)
               Enddo
               call Print_Geometry(config%Molecule,lupri)
            END IF
         END IF
!
         IF (optinfo%IPrint .GT. 2) THEN
            call lsheader(lupri,'New geometry')
               Do i = 1,NAtoms
                  Config%Molecule%Atom(i)%Center(:)=COONEW(:,i)
               Enddo
               call Print_Geometry(config%Molecule,lupri)
         END IF
!
!     Calculate energy at new geometry, which is compared to predicted
!     energy(change) in UPTRAD.
!     The temporary update of optinfo%ItrNmr is in case the molecule input is
!     provided in the DALTON input file
!

!!! Vladimir: disabled         optinfo%ItrNmr = optinfo%ItrNmr + 1
         !     Renovating geometry first
         optinfo%Coordinates = COONEW
         Call Get_Energy(E,config,optinfo,H1,F,D,S,ls,NAtoms,lupri,luerr)
         optinfo%energy = E
         GEINFO(optinfo%ItrNmr+1,1) = E
         IF (IREJ .EQ. 0) THEN
            call ls_UPTRAD(REJGEO,lupri,optinfo)
            IF (.NOT. REJGEO) THEN
               FAILED = .FALSE.
               IFAILD = 0
            END IF
!
!     After the first failure, we are satisfied if the new energy is below
!     the last. No comparison with predicted energy is done.
!
         ELSE IF ((optinfo%energy .GT. optinfo%energyOld)) THEN
            IF (FAILED .AND. (IFAILD .LE. optinfo%MaxRej) .AND. &
     &           (ABS(optinfo%energy-optinfo%energyOld) .LT. 1.0D-5)) THEN
               WRITE(LUPRI,'(/A)') 'Trouble determining step,  &
     &              accepting small energy increase.'
               IFAILD = IFAILD + 1
               REJGEO = .FALSE.
            ELSE
               WRITE(LUPRI,'(/A)') &
     &              'Step rejected because energy is increasing.'
               WRITE(LUPRI,'(A,F10.5)')' Updated trust radius', optinfo%TrustRad
               REJGEO = .TRUE.
            END IF
         ELSE
            WRITE(LUPRI,'(/A)') 'Acceptable step has been found.'
            REJGEO = .FALSE.
            FAILED = .FALSE.
         END IF
!
         IF (REJGEO) THEN
            IREJ = IREJ + 1
!
!     Line search based on quadratic model
!
            GRADDI = 0.0D0
            DO 60 I = 1, optinfo%ICartCoord
               GRADDI = GRADDI + optinfo%GRDDIA(I)*optinfo%STPDIA(I)
 60         CONTINUE
            GRADDI = GRADDI/optinfo%StepNorm
!
            IF (optinfo%IPrint .GE. 12) THEN
               call lsheader(lupri,'Line search based on quadratic model')
               WRITE(LUPRI,'(A,F12.6)') &
     &              ' Energy at last geometry     : ', optinfo%energyOld
               WRITE(LUPRI,'(A,F12.6)') &
     &              ' Energy at rejected geometry : ', optinfo%energy
               WRITE(LUPRI,'(A,F12.6)') &
     &              ' Norm of rejected step       : ', optinfo%StepNorm
               WRITE(LUPRI,'(A,F12.6)') &
     &              ' Norm of gradient            : ', optinfo%GradNorm
               WRITE(LUPRI,'(A,F12.6)') &
     &              ' Gradient along step         : ', GRADDI
            END IF
!
!     The minimum for a quadratic model is calculated with the formula
!                    -f'(0)
!     x     =  -------------------
!      min     2*(f(1)-f(0)-f'(0))
!
            FAC = -0.5D0*GRADDI/(optinfo%energy-optinfo%energyOld-GRADDI)
!
!     If the factor found is very small or very large, we don't trust
!     it. The factor is replaced by "safer" (but rather atbitrary) numbers.
!
            IF (FAC .LT. 0.1D0) FAC = 0.25D0
            IF (FAC .GT. 0.9D0) FAC = 0.75D0
!
!     We have to update both steps and their norm.
!
            DO 70 I = 1, optinfo%NIntCoord
               optinfo%STPINT(I) = optinfo%STPINT(I)*FAC
 70         CONTINUE
            DO 75 I = 1, optinfo%ICartCoord
               optinfo%STPDIA(I) = optinfo%STPDIA(I)*FAC
               optinfo%STPSYM(I) = optinfo%STPSYM(I)*FAC
 75         CONTINUE
            optinfo%StepNorm = optinfo%StepNorm*FAC
!
!     We also set the trust radius equal to the new norm
!
            optinfo%TrustRad = optinfo%StepNorm
!
            WRITE(LUPRI,'(A,F12.6)') &
     &           ' Minimum for quadratic model : ', FAC
            WRITE(LUPRI,'(A,F12.6)') &
     &           ' Norm of new step            : ', optinfo%StepNorm
!
!     Finally we construct a new geometry based on the factor found
!
            DO 80 J = 1, nAtoms
               DO 85 I = 1, 3
                  COONEW(I,J)=FAC*COONEW(I,J)+(1.0D0-FAC)*COOOLD(I,J)
 85            CONTINUE
 80         CONTINUE
         END IF
         GOTO 50
      ELSE IF (REJGEO) THEN
!
!     Maximum number of allowed rejections reached
!
         GEINFO(optinfo%ItrNmr,4) = optinfo%StepNorm
         IF (optinfo%ItrNmr .LT. optinfo%IterMax) GEINFO(optinfo%ItrNmr+1,5) = optinfo%TrustRad
         GEINFO(optinfo%ItrNmr,6) = IREJ*1.0D0
!
!     If redundant internal coordinates are used, we try reducing the
!     number of dihedral angles to one third the original number (high
!     redundancy might cause problems). We only allow this once before
!     we give up (this should be viewed as an emergency solution!).
!
         IF ((optinfo%RedInt .AND. (.NOT. FAILED)) .AND. (.NOT. optinfo%ConOpt)) THEN
            FAILED = .TRUE.
            IREJ = -IREJ
            GEINFO(optinfo%ItrNmr,6) = 0.0D0
            WRITE(LUPRI,*) 'Maximum number of rejected steps (',optinfo%MaxRej, &
     &           ') reached.'
            WRITE(LUPRI,'(A)') 'No acceptable step found.'
            WRITE(LUPRI,'(/A)')'***** NOTE! *****'
            WRITE(LUPRI,'(A)')'As an emergency solution,  &
     &           the number of dihedral angles will be reduced!'
            call ls_RREDUN(lupri,optinfo)
            NEWBMT = .TRUE.
            optinfo%TrustRad = 0.5D0
            RETURN
         ELSE IF (((.NOT. optinfo%Newton) .AND. (.NOT. FAILED)) .AND. &
     &           (.NOT. optinfo%ConOpt)) THEN
            FAILED = .TRUE.
            IREJ = -IREJ
            GEINFO(optinfo%ItrNmr,6) = 0.0D0
            WRITE(LUPRI,*) 'Maximum number of rejected steps (',optinfo%MaxRej, &
     &           ') reached.'
            WRITE(LUPRI,'(A)') 'No acceptable step found.'
            WRITE(LUPRI,'(/A)')'***** NOTE! *****'
            WRITE(LUPRI,'(A)')'As a last resort,  &
     &           the Hessian is initialized to unity!'
            optinfo%EvLini = 1.0D0
            optinfo%TrustRad = 0.5D0
            RETURN
!
!     Otherwise we give up...
!
         ELSE
            call ls_PRIINF(config%Molecule,GEINFO,lupri,optinfo)
            WRITE(LUPRI,*) 'Maximum number of rejected steps (',optinfo%MaxRej, &
     &        ') reached.'
            WRITE(LUPRI,'(A)') 'No acceptable step found. Aborting.'
            call lsquit('*** FNDGEO *** No acceptable step found.',lupri)
         END IF
      END IF
End subroutine Find_Geometry
!======================!  
! Get_energy           !
!======================!
! Call LSDALTON to calculate energy
Subroutine Get_Energy(E,config,optinfo,H1,F,D,S,ls,NAtoms,lupri,luerr)
!
use precision
use lstiming
use configuration
use lsdalton_fock_module
use matrix_module
use dal_interface
use optimization_input
use ks_settings
use direct_dens_util
use initial_guess
Implicit none
Type(lsitem),target :: ls   ! General information,used only to get E and gradient
Type(Matrix), intent(inout) :: F,D,S     ! Fock,density,overlap matrices
Type(Matrix), intent(inout),target :: H1 ! One electron matrix
Type(ConfigItem), intent(inout) :: Config ! General information
Type(opt_setting) :: optinfo
Real(realk) :: E, PotNuc   ! Electronic energy and nuclear potential
Integer :: NAtoms,i,lupri,luerr,nbast
Logical :: do_decomp
nbast = D%nrow
  if (config%opt%cfg_incremental) call ks_free_incremental_fock()
  if (config%opt%cfg_incremental) then
    call ks_init_incremental_fock(nbast)
  endif
  do_decomp =(config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
            & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
            & config%decomp%cfg_check_converged_solution .or. &
            & config%decomp%cfg_rsp_nexcit > 0 .or. config%integral%locallink) 
  if (do_decomp) then
     call decomp_shutdown(config%decomp)
  endif

!call mat_zero(F)
!call mat_zero(H1) 
!call mat_zero(S)
!call mat_zero(D)
!
! First we renovate coordinates
!
call II_free_setting(ls%setting)
call II_init_setting(ls%setting)
! eventually empirical dispersion correction in case of dft
!CALL II_DFTDISP(LS%SETTING,DUMMY,1,1,0,LUPRI,1)

!
Do i = 1,NAtoms
      ls%input%Molecule%Atom(i)%Center(:)=optinfo%Coordinates(:,i)
      Config%Molecule%Atom(i)%Center(:)=optinfo%Coordinates(:,i)
Enddo
!
call II_set_default_setting(ls%setting,ls%input)
!
!   Setting DFT grid equal to zero in order to recalculate
!   it at new geometry
!
If (ls%input%do_dft) then
   ls%input%dalton%grdone = 0
Endif

    if ((config%opt%cfg_start_guess == 'TRILEVEL')&
      &.or.(config%opt%cfg_start_guess == 'ATOMS')&
      &.or.config%decomp%cfg_gcbasis) then
!     Not working properly for geometry-optimization
      config%diag%cfg_restart = .FALSE.
      call trilevel_basis(config%opt,ls)
    endif
!
! New nuclear repulsion
  PotNuc = 0.D0
  CALL II_get_nucpot(lupri,luerr,ls%setting,PotNuc)
  config%opt%potnuc = POTNUC
  ls%input%potnuc = POTNUC
! New energy
Write(*,*)'CALLING FOR NEW ENERGY!'
Call II_get_overlap(lupri,luerr,ls%setting,S)
Call II_get_h1(lupri,luerr,ls%setting,H1)
lsint_fock_data%ls => ls
lsint_fock_data%H1 => H1
lsint_fock_data%lupri = lupri
lsint_fock_data%luerr = luerr

Call get_initial_dens(H1,S,D,ls,config)

  if (do_decomp) then
     call decomp_init(nbast,config%decomp)
     call dd_mat_eigenvalues_to_aux(config%decomp%cfg_unres,S)
     call decomposition(config%decomp)
  else if (config%opt%cfg_start_guess == 'TRILEVEL') then
     call mat_free(config%decomp%lcv_CMO)
  endif

  if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then !FIXME: put this somewhere else!
     call mat_init(config%av%Fprev,nbast,nbast)
     call mat_init(config%av%Dprev,nbast,nbast)
  endif

Call scfloop(H1,F,D,S,E,ls,config)
!
End subroutine Get_Energy
!==================!
! Get_Gradient     !
!==================!
Subroutine Get_Gradient(lupri,NAtoms,F,D,ls,optinfo)
!
! Calls II_get_molecular_gradient  
!
use precision
use lstiming
use configuration
use lsdalton_fock_module
use matrix_module
use dal_interface
use optimization_input
use memory_handling
Implicit none
Integer :: lupri,NAtoms,i
Type(opt_setting) :: optinfo
Type(Matrix), intent(in) :: F,D   ! Fock and density matrix
Type(lsitem) :: ls
Real(realk), pointer :: Gradient(:,:)
! Allocate gradient first
Call mem_alloc(Gradient,3,nAtoms)
Gradient = 0.D0
! Calculate gradient
Call II_get_molecular_gradient(Gradient,lupri,F,D,ls%setting,ls%input%do_dft,.TRUE.)
! Expand gradient to optinfo%GradMol
Do i = 1,NAtoms
   optinfo%GradMol(3*i-2:3*i) = Gradient(:,i)
Enddo
! Deallocate gradient
Call mem_dealloc(Gradient)
!
End subroutine Get_Gradient






