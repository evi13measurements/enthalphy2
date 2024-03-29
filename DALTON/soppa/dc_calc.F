C
C  /* Deck dc_calc */
      SUBROUTINE DC_CALC(ISYMTR,NEXCI,EXVAL,LEXVAL,
     &                   DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                   DENSAI,LDENSAI,
     &                   T2MP,LT2MP,FOCKD,LFOCKD,
     &                   WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, June 1997
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Determine "double corrected RPA" excitation energies.
C
#ifdef VAR_MPI
      use so_parutils, only: parsoppa_do_eres, my_mpi_integer,
     &                       soppa_nint
#endif
      use so_info, only: sop_excita
C
#include "implicit.h"
#ifdef VAR_MPI
#include "mpif.h"
C  IRAT in order to assign space for load-balancing
#include "iratdef.h"
#endif
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      LOGICAL   NONEWT
      CHARACTER*8 PDENS_LABEL
C
      DIMENSION EXVAL(LEXVAL)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2MP(LT2MP),     FOCKD(LFOCKD),   WORK(LWORK)
#ifdef VAR_MPI
      INTEGER   CP_ISYMTR
      INTEGER   MAXNUMJOBS
C     This array is only there to ensure that the four above variables
C     are allocated consecutively, so that it can be send together. Only
C     use it for this purpose.
C     The definition must match that in soppa_nodedriver
      INTEGER   INFO_ARRAY(4)
      EQUIVALENCE (info_array(1), cp_isymtr), (info_array(2),nit),
     &            (info_array(3), nnewtr),    (info_array(4),noldtr)
      INTEGER(MPI_INTEGER_KIND) :: IERR

C      IF (NUMPROCS .GT. 1) CALL PARQUIT('RPA(D)')
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('DC_CALC')

      KEND1 = 1
      LWORK1 = LWORK
C
C--------------------------------------------------------------
C     Initialize iteration counter, number of old trialvectors,
C     and number of new trialvectors.
C--------------------------------------------------------------
C
      NIT    = 1
C
      NOLDTR = 0
C
      NNEWTR = NEXCI
C
#ifdef VAR_MPI
C------------------------------------------------------------------
C     For MPI, we need some space in which to store the indices each
C     process is to work with in so_eres.
C------------------------------------------------------------------
C
      call mpi_comm_size( mpi_comm_world, numprocs, ierr)
      maxnumjobs = soppa_nint - min(soppa_nint, numprocs) + 1
      if ( numprocs .eq. 1 ) then
C Not a real parallel job, don't bother
         lAssignedIndices = 1
         kAssignedIndices = 0
      else
         lAssignedIndices = (maxnumjobs + 1) /IRAT
         kAssignedIndices = KEND1
         KEND1 = kAssignedIndices + lAssignedIndices
         LWORK1 = LWORK - KEND1
         CALL SO_MEMMAX ('DC_CALC.1A',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('DC_CALC.1A',' ',KEND1,LWORK)
      endif
#endif
C
C---------------------------------------------------------
C     Make E[2] linear transformation of RPA eigenvectors.
C---------------------------------------------------------
C
#ifdef VAR_MPI
C In parallel, send slaves to so_eres
C
         call mpi_bcast( parsoppa_do_eres, 1, my_mpi_integer, 0,
     &                   mpi_comm_world, ierr )
C ISYMTR is a non-local parameter, we need to copy it to the info-array
         CP_ISYMTR = ISYMTR
         CALL MPI_BCAST( INFO_ARRAY, 4, MY_MPI_INTEGER, 0,
     &                   MPI_COMM_WORLD, IERR)
#endif
CPi 18.08.16: Copy from so_lrsolv
      CALL GETTIM (DUMMY,WTIMES)
      DTIME      = SECOND()
      CALL SO_ERES('DCRPA',0,NEXCI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &             T2MP,LT2MP,FOCKD,LFOCKD,DENSAI,LDENSAI,1,ISYMTR,
#ifdef VAR_MPI
     &             WORK(kAssignedIndices),maxnumjobs,
#endif
     &             WORK(KEND1),LWORK1)
      DTIME      = SECOND()   - DTIME
      SOTIME(35) = SOTIME(35) + DTIME
      CALL GETTIM (DUMMY,WTIMEE)
      SOWTIM(1)  = SOWTIM(1)  + WTIMEE - WTIMES
C
C-----------------------------------------------------------
C     Make S[2] linear transformation of trialvectors giving
C     resultvectors.
C-----------------------------------------------------------
C
      DTIME      = SECOND()
      CALL DC_SRES(0,NEXCI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &             ISYMTR,WORK,LWORK)
      DTIME      = SECOND()   - DTIME
      SOTIME(40) = SOTIME(40) + DTIME
C
C-----------------------------------------------------------------
C     Calculate diagonal parts of E[2] and S[2] and write to disk.
C-----------------------------------------------------------------
C
      DTIME      = SECOND()
      CALL SO_DIAG('DCRPA',FOCKD,LFOCKD,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &             ISYMTR,WORK,LWORK)
      DTIME      = SECOND()   - DTIME
      SOTIME(31) = SOTIME(31) + DTIME
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LTR1E = NT1AM(ISYMTR)
      LTR1D = NT1AM(ISYMTR)
C      LTR2E = NT2AM(ISYMTR)
C      LTR2D = NT2AM(ISYMTR)
CPi 15.06.16
      LTR2E = N2P2HOP(ISYMTR)
      LTR2D = N2P2HOP(ISYMTR)
Cend-Pi
C
      KTR1   = 1
      KRES1  = KTR1   + LTR1E
      KRESO1 = KRES1  + LTR1E
      KRES2  = KRESO1 + LTR1E
      KEND1  = KRES2  + LTR2E
      LWORK1 = LWORK  - KEND1
C
      CALL SO_MEMMAX ('DC_CALC',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('DC_CALC',' ',KEND1,LWORK)
C
C----------------
C     Open files.
C----------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
      CALL SO_OPEN(LURS1E,FNRS1E,LTR1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LTR1D)
      CALL SO_OPEN(LURO1E,FNRO1E,LTR1E)
      CALL SO_OPEN(LURO1D,FNRO1D,LTR1D)
      CALL SO_OPEN(LURS2E,FNRS2E,LTR2E)
      CALL SO_OPEN(LURS2D,FNRS2D,LTR2D)
C
C---------------------------
C     Loop over excitations.
C---------------------------
C
      DO 100 IEXCI = 1,NEXCI
C
C-----------------------------------------------------------------------
C        Calculate 1p-1h part of double corrected RPA excitation energy.
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KTR1),  LTR1E,LUTR1E,FNTR1E,IEXCI)
         CALL SO_READ(WORK(KRES1), LTR1E,LURS1E,FNRS1E,IEXCI)
         CALL SO_READ(WORK(KRESO1),LTR1E,LURO1E,FNRO1E,IEXCI)
C
         CALL DAXPY(LTR1E,-EXVAL(IEXCI),WORK(KRESO1),1,WORK(KRES1),1)
C
         EX1P1H = DDOT(LTR1E,WORK(KRES1),1,WORK(KTR1),1)
C
C-----------------------------------------------------------------------
C        Calculate 1h-1p part of double corrected RPA excitation energy.
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KTR1),  LTR1D,LUTR1D,FNTR1D,IEXCI)
         CALL SO_READ(WORK(KRES1), LTR1D,LURS1D,FNRS1D,IEXCI)
         CALL SO_READ(WORK(KRESO1),LTR1D,LURO1D,FNRO1D,IEXCI)
C
         CALL DAXPY(LTR1D,-EXVAL(IEXCI),WORK(KRESO1),1,WORK(KRES1),1)
C
         EX1H1P = DDOT(LTR1D,WORK(KRES1),1,WORK(KTR1),1)
C
C-----------------------------------------------------------------------
C        Calculate 2p-2h part of double corrected RPA excitation energy.
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KRES2), LTR2E,LURS2E,FNRS2E,IEXCI)
C
         IF (TRIPLET) THEN
C
            CALL DC_OMEC(EX2P2H,WORK(KRES2),LTR2E,EXVAL(IEXCI),
     &                ISYMTR,WORK(KEND1),LWORK1)
C
         ELSE
C
            CALL CC_OMEC(EX2P2H,WORK(KRES2),EXVAL(IEXCI),
     &                WORK(KEND1),LWORK1,ISYMTR)
C
C--------------------------------------------------------------------
C        Calculate the first order 2p-2h eigenvector in RPA(D) theory
C        and write to output.
C--------------------------------------------------------------------
C
            CALL DC_R1VEC(WORK(KRES2), LTR2E,EXVAL(IEXCI),ISYMTR,
     &                    WORK(KEND1),LWORK1)
C
         END IF
C
         CALL SO_WRITE(WORK(KRES2), LTR2E, LUTR2E,FNTR2E, IEXCI)
C
C-----------------------------------------------------------------------
C        Calculate 2h-2p part of double corrected RPA excitation energy.
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KRES2), LTR2D,LURS2D,FNRS2D,IEXCI)
C
         IF (TRIPLET) THEN
C
            CALL DC_OMEC(EX2H2P,WORK(KRES2),LTR2D,EXVAL(IEXCI),
     &                ISYMTR,WORK(KEND1),LWORK1)
C
         ELSE
C
            CALL CC_OMEC(EX2H2P,WORK(KRES2),EXVAL(IEXCI),
     &                WORK(KEND1),LWORK1,ISYMTR)
C
C--------------------------------------------------------------------
C        Calculate the first order 2h-2p eigenvector in RPA(D) theory
C        and write to output.
C--------------------------------------------------------------------
C
            CALL DC_R1VEC(WORK(KRES2), LTR2D,EXVAL(IEXCI),ISYMTR,
     &                    WORK(KEND1),LWORK1)
C
         END IF
C
         CALL SO_WRITE(WORK(KRES2), LTR2D, LUTR2D,FNTR2D, IEXCI)
C
C-----------------------------------------------------------
C        Add contributions to give RPA(D) excitation energy.
C-----------------------------------------------------------
C
         IF (TRIPLET) THEN
C
            EXVAL(IEXCI) = EX1P1H + EX1H1P + EX2P2H + EX2H2P
C
         ELSE
C
            EXVAL(IEXCI) = EX1P1H + EX1H1P + EX2P2H/TWO + EX2H2P/TWO
C
         END IF
C
  100 CONTINUE
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
      CALL SO_CLOSE(LURO1E,FNRO1E,'KEEP')
      CALL SO_CLOSE(LURO1D,FNRO1D,'KEEP')
      CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
      CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
C
C------------------------------------
C     Flush the standard output file.
C------------------------------------
C
      CALL FLSHFO(LUPRI)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('DC_CALC')
C
      RETURN
      END
