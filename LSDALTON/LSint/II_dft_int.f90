!> @file
!> Module containing main exchange-correlation integral driver, and routines to evaluate AOs and electron-densities
!> \brief 
!> \author T. Kjaergaard
!> \date 2009 
MODULE IIDFTINT
use memory_handling
use precision
use TYPEDEF
use dft_type
CONTAINS
!> \brief wrapper exchange-correlation integral routine that build basinf.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
     & CB, DFTDATA)
use BUILDAOBATCH
IMPLICIT NONE
INTEGER       :: LUPRI,IPRINT,NBAST,NDMAT,NGEODRV
REAL(REALK)   :: DMAT(NBAST,NBAST,NDMAT)
LOGICAL       :: DFT_ISGGA,DOGGA,DOLND !do london
EXTERNAL CB !NAME OF SUBROUTINE TO CALL 
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
TYPE(BASINF)  :: BAS
TYPE(LSSETTING) :: SETTING
!
INTEGER            :: NDER,NTYPSO,IGEODRV,NSOB
EXTERNAL DFT_ISGGA

CALL BUILD_BASINF(LUPRI,IPRINT,BAS,SETTING,SETTING%scheme%GRDONE)

DOGGA = DFT_ISGGA() ! C code
IGEODRV = NGEODRV
IF (DOGGA) IGEODRV = IGEODRV + 1
CALL II_SETUPSOS(IGEODRV,DOLND,NBAST,NDER,NTYPSO,NSOB)
CALL II_DFTINT1(LUPRI,IPRINT,DMAT,NBAST,BAS%KMAX,BAS%NATOMS,BAS%NHTYP,&
     &SETTING%SCHEME%GRDONE,BAS%MXPRIM,NDMAT,NGEODRV,DOLND,CB,DFTDATA,&
     &DFTDATA%ELECTRONS,&
     &BAS%X,BAS%Y,BAS%Z,BAS%CHARGE,BAS%NCENT,BAS%NHKT,BAS%NUCO,&
     &BAS%PRIEXP,BAS%JSTRT,BAS%RSHEL,BAS%CENT,BAS%NSTART,BAS%CC,&
     &BAS%Ushells,BAS%CCSTART,BAS%CCINDEX,DFTDATA%ENERGY,&
     &SETTING%SCHEME%RADINT,SETTING%SCHEME%ANGMIN,&
     &SETTING%SCHEME%ANGINT,SETTING%SCHEME%HRDNES,SETTING%SCHEME%NOPRUN,&
     &SETTING%SCHEME%ITERATIONS,NDER,NTYPSO,NSOB,DOGGA,&
     &SETTING%SCHEME%DFTELS,SETTING%MOLECULE(1)%p%nELECTRONS,SETTING%SCHEME%DFTHRI,&
     &SETTING%SCHEME%RHOTHR,SETTING%SCHEME%TURBO)

CALL FREE_BASINF(BAS)

END SUBROUTINE II_DFTINT

!> \brief main exchange-correlation integral routine
!> \author T. Kjaergaard
!> \date 2010
!>
!>  The routine first generates a grid (II_OPNQUA) which means that it
!>  generates a number of gridpoints/coordinates with an associated weight
!>  then follows a loop over all gridpoints - the actual implementation is 
!>  a loop where in each iteration a batch of gridpoints+weights is 
!>  read from file (REAQUA). In addition to the coordinates and weights there is
!>  also read which shells contribute to this batch of gridpoints. 
!>  The shells are transformed into orbitalindexes (BLOCKS) which state
!>  which orbitals contribute to the gridpoints (which orbitals are active). 
!>  the number of these active orbitals are determined (NactBas) 
!>  INXACT is then created which for a given active index give the orbitalindex
!>  An active Dmat is then created from the full Dmat.
!>  Next the (active) gaussian atomic orbitals are created (GAO) for each gridpoint 
!>  The GAO have dimensions (ngridpoints,NactBAS,NTYPSO)
!>  NTYPSO is 1 for LDA and 4 for GGA (because we both need the atomic orbitals and 
!>  the first geometrical derivatives). In general  
!>
!>  GAO: evaluated orbitals for a batch of grid points.
!>  GAO(:,:,1) contains orbital values.
!>  GAO(:,:,2:4) contains first geom. derivatives. - if requested.
!>  GAO(:,:,5:10) contains second derivatives - if requested.
!>  GAO(:,:,11:20) contains third derivatives - if requested.
!>  After requested geometric derivatives, london related derivatives are placed.
!>  
!>  From the GAO the density (RHO) is calculated for each grid point and the 
!>  gradient of the density (GRAD) is calculated. 
!>  The number of electrons is calculated (an estimate of the error of the grid)
!>  and a generic function call CB is called. CB is an input argument, so when 
!>  you call II_DFTINT1 you can call it with for instance II_DFT_KSMGGA, and 
!>  II_DFTINT1 therefore calls II_DFT_KSMGGA. II_DFT_KSMGGA and all other 
!>  'worker' routines are for now in the II_dft_ksm.f90 file. 
SUBROUTINE II_DFTINT1(LUPRI,IPRINT,DMAT,NBAST,KMAX,NATOMS,NHTYP,GRDONE,&
     &MXPRIM,NDMAT,NGEODRV,DOLND,CB,DFTDATA,ELECTRONS,&
     &X,Y,Z,CHARGE,NCENT,NHKT,NUCO,PRIEXP,JSTRT,RSHEL,CENT,NSTART,CC,ushells,&
     &CCSTART,CCINDEX,ENERGY,RADINT,ANGMIN,ANGINT,IHARDNESS,LNOPRUNE,ITERATIONS,&
     &NDER,NTYPSO,NSOB,DOGGA,DFTELS,NELECTRONS,DFTHRI,RHOTHR,TURBO) 
use ks_settings
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER  :: IPRINT
!> the density matrix
REAL(REALK) :: DMAT(NBAST,NBAST,NDMAT)
!> number of basis functions
INTEGER  :: NBAST
!> maximum number of shells
INTEGER  :: KMAX
!> number of atoms
INTEGER  :: natoms
!> the maximum angular momentum
INTEGER  :: NHTYP
!> if the grid is done GRDONE=1 else GRDONE=0
INTEGER :: GRDONE
!> MXPRIM is the total number of (unique) primitive orbitals
INTEGER  :: MXPRIM
!> number of density matrices
INTEGER  :: NDMAT
!> the order of geometrical derivative
INTEGER  :: NGEODRV
!> do london derivative GAOs?
LOGICAL  :: DOLND 
!> the NAME OF external SUBROUTINE TO CALL 
EXTERNAL CB 
!> contains the data that must be given to the CB routine
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> number of electrons from nummerical integration
REAL(REALK) :: ELECTRONS
!> X coordinate for each atom
REAL(REALK):: X(natoms)
!> Y coordinate for each atom
REAL(REALK):: Y(natoms)
!> Z coordinate for each atom
REAL(REALK):: Z(natoms)
!> charge for each atom used in grid-generation
INTEGER  :: CHARGE(natoms)
!> which atomic center the shell is attached to
INTEGER    :: NCENT(KMAX)
!> the angular momentum for each shell
INTEGER    :: NHKT(KMAX)
!> the number of primitives for each shell
INTEGER    :: NUCO(KMAX)
!> the unique primitve exponents
REAL(REALK):: PRIEXP(MXPRIM)
!> the index to start in PRIEXP(MXPRIM) for a given shell index 
INTEGER :: JSTRT(KMAX)
!> the radius of each shell (used in grid-generation)
REAL(REALK):: RSHEL(KMAX)
!> X,Y,Z coordinate for each shell
REAL(REALK):: CENT(3,KMAX)
!> for a given shell index it gives the corresponding starting orbital
INTEGER    :: NSTART(KMAX)
!> the number of unique shells (also different contraction coefficients matrices) 
INTEGER  :: ushells
!> a contraction coefficient matrix for all the unique shells
TYPE(LSMATRIX):: CC(ushells)
!> hmmm I do not remember
INTEGER :: CCSTART(ushells)
!> contraction index for given shell, for shellindex => unique shellindex
INTEGER :: CCINDEX(KMAX)
!> the energy
real(realk) :: ENERGY
!> some grid generation specification
REAL(REALK) :: radint
!> some grid generation specification
INTEGER :: angmin
!> some grid generation specification
INTEGER :: angint
!> some grid generation specification
INTEGER :: IHARDNESS
!> grid point pruning 
LOGICAL  :: LNOPRUNE 
!> How many iterations we do
INTEGER  :: ITERATIONS
!> level of geometrical derivatives
INTEGER  :: NDER
!> number of GAOs and geomderiv derivativ GAOs
INTEGER  :: NTYPSO
!> place in GAO to place LONDON contrib
INTEGER  :: NSOB
!> is it a GGA calc or a LDA calc
LOGICAl  :: DOGGA
!> threshold for electrons
REAL(REALK) :: DFTELS
!> the number of electrons we should have
INTEGER  :: NELECTRONS
!> threshold for value of GAOs
REAL(REALK) :: DFTHRI
!> threshold for value of RHO
REAL(REALK) :: RHOTHR
!> some grid generation specification
INTEGER  :: TURBO
!
! choose reasonably large. Exceeding this limit means that boxes are too large.
Integer, parameter :: MXBLLEN=128
Integer, parameter :: NBUFLEN=20000
REAL(REALK),pointer :: COOR(:,:), WEIGHT(:), COOR_pt(:,:)
REAL(REALK) :: TELECTRONS,ERROR
INTEGER     :: NSHLBLCK(2,KMAX)
INTEGER :: BLOCKS(2,KMAX)
LOGICAL     :: CHECKELS,SETIT
INTEGER :: NPOINTS,NLEN,NSHELL,NROW,NCOL
INTEGER :: iprune,L_prev,L_curr
INTEGER :: IPT,spSIZE,L,NCURLEN,I,J
REAL(REALK),pointer :: SPHMAT(:)
INTEGER,pointer     :: SPINDEX(:)
REAL(REALK) :: CPUTIME,CPU2,CPU1,WALLTIME,WALL2,WALL1,DMAX,GAOMAX
INTEGER   :: KCKTA,KHKTA,NHKTA,XX,IT,IBUF,IBUF_PREV,IDUM,ILEN,K,NactBAS,NRED,NRED2,IDMAT
LOGICAL   :: LDUM
REAL(REALK),pointer :: GAO(:,:,:)
REAL(REALK),pointer :: RHOA(:,:), GRADA(:,:,:)
REAL(REALK),pointer :: ACTIVE_DMAT(:,:,:)
INTEGER,pointer :: INXACT(:)

TELECTRONS = 0.D0
IT=0
SETIT=.FALSE.
IF(GRDONE .EQ. 0)SETIT=.TRUE.
!pruning: per default on
IPRUNE = 1
IF (LNOPRUNE) IPRUNE = 0
CALL LS_GETTIM(CPU1,WALL1)
CALL II_OPNQUA(NBAST,radint,angmin,angint,ihardness,iprune,natoms,& 
     &X,Y,Z,Charge,GRDONE,NCENT,NHKT,NUCO,NHTYP,KMAX,MXPRIM,PRIEXP,    & 
     &JSTRT,RSHEL,IT,TURBO,LUPRI) ! C code in grid-gen.c
CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL LS_TIMTXT('>>>  WALL Time used in gridgeneration is  ',WALLTIME,LUPRI)
CALL LS_TIMTXT('>>>  CPU  Time used in gridgeneration is  ',CPUTIME,LUPRI)

!If lsopen/close are to be commented back in, remember 'use files'! /Stinne
!!$IF (incremental_scheme) THEN
!!$   IBUF   = 0
!!$   L_curr = -1
!!$   CALL LSOPEN(L_curr,rho_curr2,'UNKNOWN','UNFORMATTED')
!!$   L_prev = -1
!!$   CALL LSOPEN(L_prev,rho_prev2,'UNKNOWN','UNFORMATTED')
!!$   CALL svap_strings(rho_curr2,rho_prev2,4)
!!$ENDIF
! I HAD TO COMMENT OUT THE INCREMENTAL SCHEME BECAUSE IT CONFLICTS WITH OPENMP
! AND I AM NOT SURE HOW TO HAVE BOTH ON THE SAME TIME. TK 

IF(SETIT)THEN
   ITERATIONS=IT
ELSE
   IT=ITERATIONS
ENDIF

IF(NHTYP .GE. 3)THEN
   spSIZE=0
   DO L=2,NHTYP-1
      spSIZE=spSIZE+(2*L+1)*(L+1)*(L+2)/2
   ENDDO
   call mem_alloc(SPHMAT,spSIZE)
   !ALLOCATE(SPHMAT(spSIZE))
   CALL LS_DZERO(SPHMAT,spSIZE)
   call mem_alloc(SPINDEX,NHTYP)
   !ALLOCATE(SPINDEX(NHTYP))
   CALL Build_PRECALCULATED_SPHMAT(LUPRI,NHTYP-1,spSIZE,SPHMAT,SPINDEX)
ELSE
   spSIZE = 9
   call mem_alloc(SPHMAT,spSIZE)
   !ALLOCATE(SPHMAT(spSIZE))
   call mem_alloc(SPINDEX,2)
   !ALLOCATE(SPINDEX(2))
   SPINDEX(1)=1
   SPINDEX(2)=1
ENDIF

!$OMP PARALLEL PRIVATE(XX,NSHELL,NSHLBLCK,COOR,COOR_pt,WEIGHT,IPT,NCURLEN&
!$OMP ,ACTIVE_DMAT,BLOCKS,RHOA,GRADA,GAO,NLEN,NactBAS,INXACT) &
!$OMP REDUCTION(+:TELECTRONS) REDUCTION(+:ENERGY) 
call mem_alloc(RHOA,MXBLLEN,NDMAT)
call mem_alloc(GRADA,3,MXBLLEN,NDMAT)
call mem_alloc(WEIGHT,NBUFLEN)
call mem_alloc(COOR,3,NBUFLEN)
!$OMP DO SCHEDULE(DYNAMIC,1) 
DO XX=1,IT
!$OMP CRITICAL
   CALL REAQUA(NSHELL,NSHLBLCK,NBUFLEN,COOR,WEIGHT,NLEN)
!$OMP END CRITICAL
!  CAlls grid_getchunk_blocked() which reads grid data from the grid file also
!  with screening information if only nblocks and shlblocks are
!  provided. The data is read to coor array that will contain
!  the grid point coordinates and to weight that will contain
!  the associated weights. These arrays must be preallocated and have
!  length maxlen. Now, the information about the basis
!  function shells relevant for this chunk of grid points is returned
!  via nshell and nshlblck arguments. We do not
!  return the active shell numbers each other separately. They can be
!  packed instead into blocks of consecutive active shells. Argument
!  nblock will be set to the number of these blocks and shlblocks
!  entries will be filled with the beginnings (nshlblck[1][:]) and
!  ends (nshlblck[2][:]) of the blocks of active shells.

   IF(NLEN .LE. 0) THEN
      CALL LSQUIT('SOMETHING WRONG IN THE DFT LOOP',lupri)
   ENDIF
   CALL SHELL_TO_ORB(LUPRI,NSHELL,NSHLBLCK,BLOCKS,NSTART,KMAX,NBAST) !OUTPUT: BLOCKS
   CALL DETERMINE_NACTIVEORB(NactBAS,NSHELL,BLOCKS,KMAX) !OUTPUT: NactBAS 
   call mem_alloc(INXACT,NactBAS)
   CALL BUILD_INXACT(NactBAS,NSHELL,BLOCKS,INXACT,KMAX) !OUTPUT: INXACT
   call mem_alloc(ACTIVE_DMAT,NactBAS,NactBAS,NDMAT)
   CALL CONSTRUCT_ACTIVE_DMAT(NSHELL,KMAX,BLOCKS,DMAT,NBAST,ACTIVE_DMAT,NactBAS,NDMAT,INXACT,DMAX)
   CALL SHELL_TO_ACTORB(NSHELL,BLOCKS,KMAX) !CHANGE ORBBLOCKS
   DO IPT = 1, NLEN, MXBLLEN
      NCURLEN=MIN(MXBLLEN,NLEN-IPT+1)
      call mem_alloc(GAO,NCURLEN,NactBAS,NTYPSO) 
      COOR_pt => COOR(:,IPT:(IPT+NCURLEN-1))
      CALL II_BLGETSOS(LUPRI,NCURLEN,GAO,COOR_pt,&
      &                NSHELL,NSHLBLCK,NactBAS,DOLND,DOGGA,DFTHRI,0,&
      &                NTYPSO,NSOB,NHTYP,spSIZE,SPHMAT,SPINDEX,KMAX,NSTART,NHKT&
      &                ,CENT,JSTRT,CC,ushells,CCSTART,CCINDEX,MXPRIM,&
      &                PRIEXP,NDER)
      ! RETURNS: GAO: evaluated orbitals for a batch of grid points.
      !     GAO(:,:,1) contains orbital values.
      !     GAO(:,:,2:4) contains first geom. derivatives. - if requested.
      !     GAO(:,:,5:10) contains second derivatives - if requested.
      !     GAO(:,:,11:20) contains third derivatives - if requested.
      !     After requested geometric derivatives, london related derivatives
      !     are placed.
      IF(DOGGA) THEN
         CALL II_GETRHO_BLOCKED_GGA(LUPRI,ACTIVE_DMAT,NactBAS,GAO,NTYPSO,&
              &NSHELL,BLOCKS,NCURLEN,NDMAT,RHOA,GRADA,RHOTHR,MXBLLEN)
         ! computes rho_a and gradients of rho, Dmat is a density matrix
         ! (it can be a total density and then one will get total density,
         ! or it can be an alpha/beta density    assert(NTYPSO>=NRHO)
      ELSE
         CALL II_GETRHO_BLOCKED_LDA(LUPRI,ACTIVE_DMAT,NactBAS,GAO,NTYPSO,&
              &NSHELL,BLOCKS,NCURLEN,NDMAT,RHOA,RHOTHR,MXBLLEN)
      END IF
      !Incremental scheme (rho and grad calculated with difference density
      !instead of the full density-matrix)
!!$      IF (incremental_scheme) THEN
!!$         IF (do_increment) THEN
!!$            IBUF_PREV = IBUF
!!$            DO I=1,NCURLEN
!!$               IF (MOD(IBUF_PREV,NBUFLEN).EQ.0) THEN
!!$                  READ(L_prev) ILEN
!!$                  READ(L_prev) (RHOA_BUF_PREV(J),J=1,ILEN)
!!$                  IF (DOGGA) THEN
!!$                     READ(L_prev) ((GRADA_BUF_PREV(J,K),J=1,3),K=1,ILEN)
!!$                  ENDIF
!!$                  IBUF_PREV = 0
!!$               ENDIF
!!$               IBUF_PREV = IBUF_PREV + 1
!!$               RHOA(I) = RHOA(I) + RHOA_BUF_PREV(IBUF_PREV)
!!$               IF (DOGGA) THEN
!!$                  GRADA(1,I) = GRADA(1,I) + GRADA_BUF_PREV(1,IBUF_PREV)
!!$                  GRADA(2,I) = GRADA(2,I) + GRADA_BUF_PREV(2,IBUF_PREV)
!!$                  GRADA(3,I) = GRADA(3,I) + GRADA_BUF_PREV(3,IBUF_PREV)
!!$               ENDIF
!!$            ENDDO
!!$         ENDIF
!!$         DO I=1,NCURLEN
!!$            IF (IBUF.EQ.NBUFLEN) THEN
!!$               WRITE(L_curr) NBUFLEN
!!$               WRITE(L_curr) (RHOA_BUF(J),J=1,NBUFLEN)
!!$               IF (DOGGA) THEN
!!$                  WRITE(L_curr) ((GRADA_BUF(J,K),J=1,3),K=1,NBUFLEN)
!!$               ENDIF
!!$               IBUF = 0
!!$            ENDIF
!!$            IBUF = IBUF + 1
!!$            RHOA_BUF(IBUF) = RHOA(I)
!!$            IF (DOGGA) THEN
!!$               GRADA_BUF(1,IBUF) = GRADA(1,I)
!!$               GRADA_BUF(2,IBUF) = GRADA(2,I)
!!$               GRADA_BUF(3,IBUF) = GRADA(3,I)
!!$            ENDIF
!!$         ENDDO
!!$      ENDIF
! I HAD TO COMMENT OUT THE INCREMENTAL SCHEME BECAUSE IT CONFLICTS WITH OPENMP
! AND I AM NOT SURE HOW TO HAVE BOTH ON THE SAME TIME. TK 
      DO IDMAT = 1,NDMAT
         DO I = 1, NCURLEN
            TELECTRONS = TELECTRONS + WEIGHT(IPT+I-1)*RHOA(I,IDMAT) 
         END DO
      ENDDO
      CALL CB(LUPRI,NCURLEN,NSHELL,BLOCKS(:,:),INXACT(:),NactBas,NBAST,NDMAT,ACTIVE_DMAT(:,:,:),NTYPSO,GAO(:,:,:),&
           &RHOA(:,:),GRADA(:,:,:),MXBLLEN,COOR(:,IPT:IPT+NCURLEN-1),WEIGHT(IPT:IPT+NCURLEN-1),ENERGY,DFTDATA,RHOTHR,DFTHRI)
      call mem_dealloc(GAO)
   ENDDO
   call mem_dealloc(INXACT)
   call mem_dealloc(ACTIVE_DMAT)
ENDDO
!$OMP END DO
   call mem_dealloc(RHOA)
   call mem_dealloc(GRADA)
   call mem_dealloc(WEIGHT)
   call mem_dealloc(COOR)
!$OMP END PARALLEL

call mem_dealloc(SPHMAT)
call mem_dealloc(SPINDEX)

CALL CLSQUA 
!
! Test on the number of electrons
!   
ELECTRONS = TELECTRONS
ERROR  = ELECTRONS - REAL(NELECTRONS)
IF (DABS(ERROR) .GT. DFTELS*REAL(NELECTRONS)) THEN
WRITE (LUPRI,'(4(/2X,A,F14.6),/2X,A)')                              &
     &' Number of electrons from numerical integration:',ELECTRONS, &
     &' Number of electrons from orbital occupations:  ',REAL(NELECTRONS),&
     &' Error in the number of electrons:              ',ERROR,     &
     &' Error larger than DFTELS (set input):          ',DFTELS,    &
     &' Calculation aborted.'
CALL LSQUIT                                                           &
     &    ('Wrong number of electrons in DFTINT. Calculation aborted.',lupri)
END IF

!!$IF (incremental_scheme) THEN
!!$   IF (IBUF.GT.0) THEN
!!$      WRITE(L_curr) IBUF
!!$      WRITE(L_curr) (RHOA_BUF(J),J=1,IBUF)
!!$      IF (DOGGA) THEN
!!$         WRITE(L_curr) ((GRADA_BUF(J,K),J=1,3),K=1,IBUF)
!!$      ENDIF
!!$   ENDIF
!!$   CALL LSCLOSE(L_prev,'KEEP')
!!$   CALL LSCLOSE(L_curr,'KEEP')
!!$ENDIF
! I HAD TO COMMENT OUT THE INCREMENTAL SCHEME BECAUSE IT CONFLICTS WITH OPENMP
! AND I AM NOT SURE HOW TO HAVE BOTH ON THE SAME TIME. TK 
END SUBROUTINE II_DFTINT1

!> \brief determine the number of active orbitals
!> \author T. Kjaergaard
!> \date 2010
!>
!>  orbitals that contribute to the current gridpoints
SUBROUTINE DETERMINE_NACTIVEORB(NactBAS,NSHELL,BLOCKS,KMAX)
IMPLICIT NONE
INTEGER,intent(in) :: NSHELL,KMAX
INTEGER,intent(in) :: BLOCKS(2,KMAX)
INTEGER,intent(out):: NactBAS
!
INTEGER            :: IBL,IORB

NactBAS = 0
DO IBL = 1, NSHELL
   DO IORB = BLOCKS(1,IBL),BLOCKS(2,IBL) !ORBITALINDEX
      NactBAS = NactBAS + 1
   ENDDO
ENDDO

END SUBROUTINE DETERMINE_NACTIVEORB

!> \brief build indexhandling from active index to orbitalindex
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE BUILD_INXACT(NactBAS,NSHELL,BLOCKS,INXACT,KMAX) 
IMPLICIT NONE
INTEGER,intent(in) :: NSHELL,KMAX,NactBAS
INTEGER,intent(in) :: BLOCKS(2,KMAX)
INTEGER,intent(out):: INXACT(NactBAS)
!
INTEGER            :: IBL,IORB,IBASIS

IBASIS = 0
DO IBL = 1, NSHELL
   DO IORB = BLOCKS(1,IBL),BLOCKS(2,IBL) !ORBITALINDEX
      IBASIS = IBASIS + 1
      INXACT(IBASIS) = IORB !FOR A GIVEN ACTIVEINDEX - THE CORRESPONDING ORBITALINDEX
   ENDDO
ENDDO

END SUBROUTINE BUILD_INXACT

!> \brief construct active Dmat
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE CONSTRUCT_ACTIVE_DMAT(NSHELL,KMAX,BLOCKS,DMAT,NBAST,ACTIVE_DMAT,NactBAS,NDMAT,INXACT,DMAX)
IMPLICIT NONE 
INTEGER,intent(in)     :: NSHELL,KMAX,NBAST,NactBAS,NDMAT
INTEGER,intent(in)     :: BLOCKS(2,KMAX)
INTEGER,intent(in)    :: INXACT(NactBAS)
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
REAL(REALK),intent(out):: ACTIVE_DMAT(NactBAS,NactBAS,NDMAT),DMAX
!
INTEGER               :: IMAT,IACT,JACT,I,J
REAL(REALK)           :: DTMP
DMAX=0.D0
DO IMAT = 1,NDMAT
   DO JACT = 1,NactBAS
      J = INXACT(JACT)
      DO IACT = 1, NactBAS
         I = INXACT(IACT)
         DTMP = DMAT(I,J,IMAT)
         ACTIVE_DMAT(IACT,JACT,IMAT) = DTMP
         DMAX = MAX(DMAX,DABS(DTMP))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE CONSTRUCT_ACTIVE_DMAT

!> \brief determine ntypso and NSOB
!> \author T. Kjaergaard
!> \date 2010
!>
!>  NTYPSO determine how many GAOs, geometrical derivatives and 
!>  london derivatives there is needed. NSOB is the place in GAO
!>  to place mag derivative GAOs 
!>
SUBROUTINE II_SETUPSOS(GEODRV,DOLND,NBAST,NDER,NTYPSO,NSOB)
IMPLICIT NONE
  INTEGER :: NDER,NTYPSO,NBAST
  INTEGER :: NSOB
! index to keep track of where london derivatives fit in GAO
  INTEGER :: GEODRV
!     GEODRV - 0 for just orbital values, 1 for first order orbital
!     derivatives, 2 for laplacian.
  LOGICAL :: DOLND
  !Compute derivatives wrt magnetic field ("london" derivatives).
  !     the computed orbitals are ordered as follows:
  !     (orbital values=O, [dO/dx,dO/dy,dO/dz, [d0/dxx, ..]],
  !      [dO/dBx, dO/dBy, dO/dBz])
  NDER = GEODRV
  IF (NDER.EQ.0) NTYPSO =  1
  IF (NDER.EQ.1) NTYPSO =  4
  IF (NDER.EQ.2) NTYPSO = 10
  IF (NDER.EQ.3) NTYPSO = 20
  IF (DOLND) THEN
     NTYPSO = NTYPSO + 3
     NSOB   = NTYPSO - 2 
     IF (NDER.GT.0) THEN
        NTYPSO = NTYPSO + 9
     END IF
  ELSE
     NSOB  = 0
  END IF

END SUBROUTINE II_SETUPSOS

!> \brief evaluates the gaussian atomic orbitals
!> \author T. Kjaergaard
!> \date 2010
!>
!>  RETURNS: GSO: evaluated orbitals for a batch of grid points.
!>     GSO(:,:,1) contains orbital values.
!>     GSO(:,:,2:4) contains first geom. derivatives. - if requested.
!>     GSO(:,:,5:10) contains second derivatives - if requested.
!>     GSO(:,:,11:20) contains third derivatives - if requested.
!>     After requested geometric derivatives, London related derivatives
!>     are placed.
!>  REWRITTEN BY T.KJAERGAARD ORIGINALLY BY  T. Helgaker sep 99, P. Salek 03
!>
SUBROUTINE II_BLGETSOS(LUPRI,NVCLEN,GAO,COOR,NBLCNT,IBLCKS,&
     &        NactBAS,DOLND,DOGGA,DFTHRI,IPRINT,NTYPSO,NSOB,NHTYP,SIZE,&
     &        SPHMAT,SPINDEX,KMAX,NSTART,NHKT,CENT,&
     &        JSTRT,CC,ushells,CCSTART,CCINDEX,MXPRIM,PRIEXP,NDER)
IMPLICIT NONE
INTEGER     :: LUPRI,NVCLEN,NBLCNT,NactBAS,IPRINT,NTYPSO,NSOB,SIZE,NHTYP,KMAX
INTEGER     :: IBLCKS(2,NBLCNT),SPINDEX(NHTYP),NSTART(KMAX),NHKT(KMAX),ushells
INTEGER     :: JSTRT(KMAX),MXPRIM,NDER,CCSTART(ushells)
INTEGER     :: CCINDEX(KMAX)
REAL(REALK) :: PRIEXP(MXPRIM)
LOGICAL     :: DOLND, DOGGA
REAL(REALK) :: GAO(NVCLEN,NactBAS,NTYPSO), COOR(3,NVCLEN),CENT(3,KMAX)
REAL(REALK) :: GAOMAX(NactBAS)
REAL(REALK) :: orig(3),SPHMAT(SIZE),CENX,CENY,CENZ
REAL(REALK),PARAMETER :: D0 = 0.0D0, DTHRS = 20.0D0
!
REAL(REALK) :: PA(3,NVCLEN), PA2(NVCLEN),DFTHRI
INTEGER     :: I,IADR,JSTA,ISHELA,KHKTA,KCKTA,NHKTA,IBL,J,oooo,IJ,SPVAL,K
INTEGER,pointer  :: LVALUE(:),MVALUE(:),NVALUE(:)
!REAL(REALK) :: CC(MXPRIM,KMAX)
TYPE(LSMATRIX) :: CC(ushells)
!Orig(1) = 0D0
!Orig(2) = 0D0
!Orig(3) = 0D0

call mem_alloc(LVALUE,1)
call mem_alloc(MVALUE,1)
call mem_alloc(NVALUE,1)
LVALUE(1)=0.d0
MVALUE(1)=0.d0
NVALUE(1)=0.d0

SPVAL=1
!!#if 1
IF(NDER.GT.1) THEN
   ! fixme - this won't scale linearly 
   CALL LS_DZERO(GAO(1,1,5),6*NactBAS*NVCLEN)
   IF(NDER.GT.2)THEN
      CALL LS_DZERO(GAO(1,1,11),10*NactBAS*NVCLEN)
   END IF
END IF
!!!#else
!! this should be activated one I have modified the II_BLGETGA2 routine
!!IF(NDER.GT.2) THEN
!!   ! fixme - this won't scale linearly 
!!   CALL LS_DZERO(GAO(1,1,11),10*NactBAS*NVCLEN)
!!END IF
!!#endif
IADR = 0
DO IBL = 1, NBLCNT
   DO ISHELA = IBLCKS(1,IBL),IBLCKS(2,IBL)
!      IORB =   NSTART(ISHELA)+1  !ORBITALINDEX NOT USED
      NHKTA = NHKT(ISHELA)
      KHKTA  = 2*(NHKTA-1)+1  
      KCKTA  = NHKTA*(NHKTA+1)/2  
      JSTA = JSTRT(ISHELA)
      CENX = CENT(1,ISHELA)!+ ORIG(1)
      CENY = CENT(2,ISHELA)!+ ORIG(2)
      CENZ = CENT(3,ISHELA)!+ ORIG(3)
      IF (NDER.GT.1.OR.NHKTA .GE. 3.OR.DOLND) THEN !d orbitals or deriv
         call mem_dealloc(LVALUE)
         call mem_dealloc(MVALUE)
         call mem_dealloc(NVALUE)
         call mem_alloc(LVALUE,KCKTA*KCKTA)
         call mem_alloc(MVALUE,KCKTA*KCKTA)
         call mem_alloc(NVALUE,KCKTA*KCKTA)
         SPVAL=KCKTA*KCKTA
         IJ=0
         DO I = 1,KCKTA
            DO J = 1,I
               IJ=IJ+1
               LVALUE(IJ)=NHKTA-I
               MVALUE(IJ)=I-J
               NVALUE(IJ)=J-1
            ENDDO
         ENDDO
      ENDIF
      DO I=1,NVCLEN !gridpoints
         PA(1,i) = COOR(1,i)-CENX 
         PA(2,i) = COOR(2,i)-CENY
         PA(3,i) = COOR(3,i)-CENZ
      END DO
      DO I=1,NVCLEN !gridpoints
         PA2(i) = PA(1,i)**2 + PA(2,i)**2 + PA(3,i)**2
      END DO
      IF (NDER.EQ.0) THEN
         CALL II_BLGETGAO(LUPRI,NVCLEN,NactBAS,NTYPSO,NHKTA,KHKTA,KCKTA,MXPRIM,JSTA,&
              & SPHMAT(SPINDEX(NHKTA):SPINDEX(NHKTA)+KHKTA*KCKTA-1),&
              & GAO,IADR,PA,PA2,DFTHRI,CC(CCINDEX(ISHELA)),&
              & CCSTART(CCINDEX(ISHELA)),PRIEXP,LVALUE,MVALUE,NVALUE,SPVAL,.TRUE.)
      ELSE IF (NDER.GT.0) THEN
         CALL II_BLGETGA1(LUPRI,NVCLEN,NactBAS,NTYPSO,NHKTA,KHKTA,KCKTA,MXPRIM,JSTA,&
              & SPHMAT(SPINDEX(NHKTA):SPINDEX(NHKTA)+KHKTA*KCKTA-1),&
              & GAO,IADR,PA,PA2,DFTHRI,CC(CCINDEX(ISHELA)),&
              & CCSTART(CCINDEX(ISHELA)),PRIEXP,LVALUE,MVALUE,NVALUE,SPVAL,.TRUE.)
         IF (NDER.GT.1) THEN
            CALL II_BLGETGA2(LUPRI,NVCLEN,NactBAS,NHKTA,KHKTA,KCKTA,MXPRIM,JSTA,&
              & SPHMAT(SPINDEX(NHKTA):SPINDEX(NHKTA)+KHKTA*KCKTA-1),&
              & GAO(:,:,5),GAO(:,:,6),&
              & GAO(:,:,7),GAO(:,:,8),&
              & GAO(:,:,9),GAO(:,:,10),IADR,&
              & PA,PA2,DFTHRI,CC(CCINDEX(ISHELA)),CCSTART(CCINDEX(ISHELA)),PRIEXP,&
              & LVALUE,MVALUE,NVALUE,SPVAL,.TRUE.)
            IF (NDER.GT.2) THEN
!               CALL II_BLGETGA3(NVCLEN,&
!               &       GAO(1,IADR,11),GAO(1,IADR,12),GAO(1,IADR,13),&
!               &       GAO(1,IADR,14),GAO(1,IADR,15),GAO(1,IADR,16),&
!               &       GAO(1,IADR,17),GAO(1,IADR,18),GAO(1,IADR,19),&
!               &       GAO(1,IADR,20),SPHMAT(SPINDEX(NHKTA)),PA,PA2,DFTHRI)
            END IF
         END IF
      END IF
      IF (DOLND) THEN
         CALL II_BLGETGB1(NVCLEN,NactBAS,IADR,KHKTA,GAO(:,:,1),PA,&
         &              GAO(:,:,NSOB),GAO(:,:,NSOB+1),GAO(:,:,NSOB+2))
         IF (DOGGA) THEN
            DO I = 1, 3
               CALL II_BLGETGB1(NVCLEN,NactBAS,IADR,KHKTA,GAO(:,:,1+I),PA,&
               & GAO(:,:,NSOB+2+I),GAO(:,:,NSOB+5+I),GAO(:,:,NSOB+8+I))
            END DO
         END IF
      END IF
      IF (NDER.GT.1.OR.NHKTA .GE. 3.OR.DOLND) THEN !d orbitals or deriv      
         call mem_dealloc(LVALUE)
         call mem_dealloc(MVALUE)
         call mem_dealloc(NVALUE)
         call mem_alloc(LVALUE,1)
         call mem_alloc(MVALUE,1)
         call mem_alloc(NVALUE,1)
         SPVAL=1
         LVALUE(1)=0.d0
         MVALUE(1)=0.d0
         NVALUE(1)=0.d0
      ENDIF
      IADR = IADR + KHKTA !UPDATE ACTIVEORBITALINDEX
   END DO
END DO

call mem_dealloc(LVALUE)
call mem_dealloc(MVALUE)
call mem_dealloc(NVALUE)

!write (lupri,*) 'integrals from BLDFTAOS '
!call output(gao(1,1,1),1,nvclen,1,nbast,nvclen,NactBAS,1,lupri)
!write (lupri,*) 'x integrals from BLDFTAOS '
!call output(gao(1,1,2),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!write (lupri,*) 'y integrals from BLDFTAOS '
!call output(gao(1,1,3),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!write (lupri,*) 'z integrals from BLDFTAOS '
!call output(gao(1,1,4),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!if (nder.eq.2) then
!   WRITE(lupri,*)'COOR(1,:)',(COOR(1,i),I=1,NVCLEN)
!   WRITE(lupri,*)'COOR(2,:)',(COOR(2,i),I=1,NVCLEN)
!   WRITE(lupri,*)'COOR(3,:)',(COOR(3,i),I=1,NVCLEN)
!   write (lupri,*) ' xx integrals from II_BLGETSOS '
!   call output(gao(1,1,5),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' xy integrals from II_BLGETSOS '
!   call output(gao(1,1,6),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' xz integrals from II_BLGETSOS '
!   call output(gao(1,1,7),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' yy integrals from II_BLGETSOS '
!   call output(gao(1,1,8),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' yz integrals from II_BLGETSOS '
!   call output(gao(1,1,9),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' zz integrals from II_BLGETSOS '
!   call output(gao(1,1,10),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!end if

END SUBROUTINE II_BLGETSOS

!> \brief determines the rho for a LDA type calc
!> \author T. Kjaergaard
!> \date 2010
!>
!>     computes  <o|dmat|o'>  i.e rho_a where dmat is a density matrix.
!>
SUBROUTINE II_GETRHO_BLOCKED_LDA(LUPRI,DMAT,NactBAS,GAO,NTYPSO,NBLOCKS,BLOCKS,NVCLEN,NDMAT,RHO,RHOTHR,MXBLLEN)

IMPLICIT NONE
INTEGER     :: NactBAS,NVCLEN,NBLOCKS,NTYPSO,LUPRI,ndmat,MXBLLEN
REAL(REALK) :: DMAT(NactBAS,NactBAS,ndmat), GAO(NVCLEN,NactBAS,NTYPSO)
INTEGER     :: BLOCKS(2,NBLOCKS)
REAL(REALK) :: TMP(NVCLEN,NactBAS), RHO(MXBLLEN,NDMAT)
INTEGER     :: IBL,ISTART,IBLEN,JBL,JSTART,JBLEN
INTEGER     :: IDX,JTOP,JDX,K,I,J
REAL(REALK) :: GAOGMX(NactBAS),RHOTHR
REAL(REALK) :: GAORED(NVCLEN,NactBAS),GAOMAX,DMAX
INTEGER     :: INXRED(NactBAS),NRED,JRED,IRED,idmat
REAL(REALK),pointer :: DRED(:,:)

! Set up maximum Gaussian AO elements
GAOMAX = 0.0d0
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      GAOGMX(I) = 0.0d0
      DO K = 1, NVCLEN
         GAOMAX = MAX(GAOMAX,ABS(GAO(K,I,1)))
         GAOGMX(I) = MAX(GAOGMX(I),ABS(GAO(K,I,1)))
      ENDDO
   ENDDO
ENDDO
DO IDMAT=1,NDMAT
   !       Set up maximum density-matrix elements
   DMAX = 0.0d0
   DO IBL=1, NBLOCKS
      DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
         DO JBL=1, NBLOCKS
            DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
               DMAX = MAX(DMAX,ABS(DMAT(I,J,IDMAT)))
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !       Set up reduced Gaussian AO's
   NRED = 0
   DO IBL=1, NBLOCKS
      DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
         IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
            NRED = NRED + 1
            INXRED(NRED) = I
            DO K = 1, NVCLEN
               GAORED(K,NRED) = GAO(K,I,1)
            ENDDO
         ENDIF
      ENDDO
   ENDDO

   IF (NRED.GT.0) THEN
      !         Set up reduced density-matrix
      ALLOCATE(DRED(NRED,NRED))
      DO JRED=1,NRED
         J = INXRED(JRED)
         DO IRED=1,NRED
            I = INXRED(IRED)
            DRED(IRED,JRED) = DMAT(I,J,IDMAT)
         ENDDO
      ENDDO
      !         First half-contraction of Gaussian AO with density-matrix
      CALL DGEMM("N","N",NVCLEN,NRED,NRED,1D0,GAORED,NVCLEN,&
           &     DRED,NRED,0.0d0,TMP,NVCLEN)
      !         Second half-contraction
      DO K = 1, NVCLEN
         RHO(K,IDMAT)= GAORED(K,1)*TMP(K,1)
      END DO
      DO I = 2, NRED
         DO K = 1, NVCLEN
            RHO(K,IDMAT)    = RHO(K,IDMAT) + GAORED(K,I)*TMP(K,I)
         END DO
      END DO
      DEALLOCATE(DRED)
      !Hack Severeal functionals does not handle a zero density - so
      !     we set these values explicitly to some small value.
      !     Should instead skip these contributions.
      DO K = 1, NVCLEN
         IF (ABS(RHO(K,IDMAT)).LE.1.0d-20) RHO(K,IDMAT) = 1.0d-20
      END DO
   ELSE
      DO K = 1, NVCLEN
         RHO(K,IDMAT) = 1.0d-20
      END DO
   ENDIF
ENDDO

END SUBROUTINE II_GETRHO_BLOCKED_LDA

!> \brief determines the rho for a GGA type calc
!> \author T. Kjaergaard
!> \date 2010
!>
!>     computes  <o|dmat|o'>
!>     i.e rho_a where dmat is a density matrix (it can be a total
!>     density end then one will get total density, or it can be an
!>     alpha/beta density.
!>     assert(NTYPSO>=NRHO)
!>
SUBROUTINE II_GETRHO_BLOCKED_GGA(LUPRI,DMAT,NactBAS,GAO,NTYPSO,NBLOCKS,BLOCKS,NVCLEN,NDMAT,RHO,GRAD,RHOTHR,MXBLLEN)
IMPLICIT NONE
INTEGER     :: NactBAS,NVCLEN,NBLOCKS,NTYPSO,LUPRI,ndmat,MXBLLEN
REAL(REALK) :: DMAT(NactBAS,NactBAS,ndmat), GAO(NVCLEN,NactBAS,NTYPSO)
INTEGER     :: BLOCKS(2,NBLOCKS)
REAL(REALK) :: RHO(MXBLLEN,NDMAT), GRAD(3,MXBLLEN,NDMAT),TMP(NVCLEN,NactBAS)
INTEGER     :: IBL,ISTART,IBLEN,JBL,JSTART,JBLEN
INTEGER     :: IDX,JTOP,JDX,K,I,J
REAL(REALK) :: GAOGMX(NactBAS),GAOMAX,DMAX,RHOTHR
REAL(REALK) :: GAORED(NVCLEN,NactBAS,4)
INTEGER     :: INXRED(NactBAS),NRED,JRED,IRED,idmat
REAL(REALK),ALLOCATABLE :: DRED(:,:)

! Set up maximum Gaussian AO elements
GAOMAX = 0.0d0
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      GAOGMX(I) = 0.0d0
      DO K = 1, NVCLEN
         GAOMAX = MAX(GAOMAX,ABS(GAO(K,I,1)))
         DO J=1,4
            GAOGMX(I) = MAX(GAOGMX(I),ABS(GAO(K,I,J)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
DO IDMAT=1,NDMAT
   !       Set up maximum density-matrix elements
   DMAX = 0.0d0
   DO IBL=1, NBLOCKS
      DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
         DO JBL=1, NBLOCKS
            DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
               DMAX = MAX(DMAX,ABS(DMAT(I,J,IDMAT)))
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !       Set up reduced Gaussian AO's
   NRED = 0
   DO IBL=1, NBLOCKS
      DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
         IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
            NRED = NRED + 1
            INXRED(NRED) = I
            DO K = 1, NVCLEN
               GAORED(K,NRED,1) = GAO(K,I,1)
               GAORED(K,NRED,2) = GAO(K,I,2)
               GAORED(K,NRED,3) = GAO(K,I,3)
               GAORED(K,NRED,4) = GAO(K,I,4)
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   
   IF (NRED.GT.0) THEN
      !         Set up reduced density-matrix
      ALLOCATE(DRED(NRED,NRED))
      DO IRED=1,NRED
         I = INXRED(IRED)
         DO JRED=1,NRED
            J = INXRED(JRED)
            DRED(IRED,JRED) = DMAT(I,J,IDMAT)
         ENDDO
      ENDDO
      !         First half-contraction of Gaussian AO with density-matrix
      CALL DGEMM("N","N",NVCLEN,NRED,NRED,1D0,GAORED,NVCLEN,&
           &     DRED,NRED,0.0d0,TMP,NVCLEN)

      !         Second half-contraction
      DO K = 1, NVCLEN
         RHO(K,IDMAT)    = GAORED(K,1,1)*TMP(K,1)
         GRAD(1,K,IDMAT) = 2*GAORED(K,1,2)*TMP(K,1)
         GRAD(2,K,IDMAT) = 2*GAORED(K,1,3)*TMP(K,1)
         GRAD(3,K,IDMAT) = 2*GAORED(K,1,4)*TMP(K,1)
      END DO
      DO I = 2, NRED
         DO K = 1, NVCLEN
            RHO(K,IDMAT)    = RHO(K,IDMAT)    +   GAORED(K,I,1)*TMP(K,I)
            GRAD(1,K,IDMAT) = GRAD(1,K,IDMAT) + 2*GAORED(K,I,2)*TMP(K,I)
            GRAD(2,K,IDMAT) = GRAD(2,K,IDMAT) + 2*GAORED(K,I,3)*TMP(K,I)
            GRAD(3,K,IDMAT) = GRAD(3,K,IDMAT) + 2*GAORED(K,I,4)*TMP(K,I)
         END DO
      END DO
      DEALLOCATE(DRED)
      !Hack Severeal functionals does not handle a zero density - so
      !     we set these values explicitly to some small value.
      !     Should instead skip these contributions.
      DO K = 1, NVCLEN
         IF (ABS(RHO(K,IDMAT)).LE.1.0d-20) RHO(K,IDMAT) = 1.0d-20
         IF (ABS(GRAD(1,K,IDMAT)).LE.1.0d-20) GRAD(1,K,IDMAT) = 1.0d-20
         IF (ABS(GRAD(2,K,IDMAT)).LE.1.0d-20) GRAD(2,K,IDMAT) = 1.0d-20
         IF (ABS(GRAD(3,K,IDMAT)).LE.1.0d-20) GRAD(3,K,IDMAT) = 1.0d-20
      END DO
   ELSE
      !Hack Severeal functionals does not handle a zero density - so
      !     we set these values explicitly to some small value.
      !     Should instead skip these contributions.
      DO K = 1, NVCLEN
         RHO(K,IDMAT) = 1.0d-20
         GRAD(1,K,IDMAT) = 1.0d-20
         GRAD(2,K,IDMAT) = 1.0d-20
         GRAD(3,K,IDMAT) = 1.0d-20
      END DO
   ENDIF
ENDDO

END SUBROUTINE II_GETRHO_BLOCKED_GGA

!> \brief evaluates the pure GAO(gaussian atomic orbitals) , so no derivatives
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_BLGETGAO(LUPRI,NVCLEN,NactBAS,NTYPSO,NHKTA,KHKTA,KCKTA,MXPRIM,JSTA,&
                     & CSP,GAO,IADR,PA,PA2,DFTHRI,CC,CCSTART,PRIEXP,LVALUE,MVALUE,&
                     &NVALUE,SPVAL,SPHRA)
IMPLICIT NONE
INTEGER                   :: NVCLEN,NactBAS,MXPRIM,KHKTA,JSTA,K,KCKTA,LUPRI,NHKTA
INTEGER                   :: SPVAL,CCSTART,IADR,ISTART,NTYPSO
REAL(REALK),PARAMETER     :: D0 = 0.0D0, D1 = 1.0D0
REAL(REALK),intent(inout) :: GAO(NVCLEN,NactBAS,NTYPSO)
REAL(REALK),intent(in)    :: PA(3,NVCLEN),PA2(NVCLEN)
REAL(REALK)               :: CSP(KHKTA,KCKTA),PRIEXP(MXPRIM)
REAL(REALK)               :: GA(NVCLEN),CINT(NVCLEN),VA(NVCLEN),VE(NVCLEN)
REAL(REALK)               :: DFTHRI,GAZ,GAXX,GAX,GAY,GAYY,GMAX,PA_1,PA_2,PA_3
REAL(REALK)               :: PRICCFVAL,PRIEXPVAL
INTEGER                   :: I,J,LVALUE(SPVAL),MVALUE(SPVAL),NVALUE(SPVAL),TEMPI
INTEGER                   :: LVALJ, MVALJ, NVALJ, LVALI, MVALI, NVALI
INTEGER                   :: ISTART1,ISTART2,ISTART3,ISTART4,ISTART5,ISTART6
LOGICAl                   :: SPHRA
TYPE(LSMATRIX)              :: CC

!     loop over primitives
ISTART=IADR
ISTART1=ISTART+1
ISTART2=ISTART+2
ISTART3=ISTART+3
ISTART4=ISTART+4
ISTART5=ISTART+5
ISTART6=ISTART+6
!CALL LS_DZERO(GA,NVCLEN)
I=1
J = JSTA+CCSTART-1+I
PRICCFVAL = CC%elms(I)
PRIEXPVAL = -PRIEXP(J)
DO K = 1, NVCLEN
   GA(K) = PRICCFVAL*DEXP(PRIEXPVAL*PA2(K))
END DO
DO I = 2,CC%nrow
   J = JSTA+CCSTART-1+I
   PRICCFVAL = CC%elms(I)
   PRIEXPVAL = -PRIEXP(J)
   DO K = 1, NVCLEN
      GA(K) = GA(K) + PRICCFVAL*DEXP(PRIEXPVAL*PA2(K))
   END DO
ENDDO

! screening based on the maximum GA value for the baych of points: GMAX \Andreas Krapp
GMAX = D0
DO K = 1, NVCLEN
   GMAX = MAX(GMAX,DABS(GA(K)))
END DO

!     
!     contracted orbitals
!
IF (NHKTA .EQ. 1) THEN      !s orbitals
   IF (GMAX .GT. DFTHRI) THEN 
      DO K = 1, NVCLEN
         GAO(K,ISTART1,1) = GA(K)
      END DO
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1,1),NVCLEN)
   END IF
ELSEIF (NHKTA .EQ. 2) THEN !p orbitals
   IF (GMAX .GT. DFTHRI) THEN
      DO K = 1, NVCLEN
         GAO(K,ISTART1,1) = PA(1,K)*GA(K)
         GAO(K,ISTART2,1) = PA(2,K)*GA(K)
         GAO(K,ISTART3,1) = PA(3,K)*GA(K)
      END DO
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1,1),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,1),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,1),NVCLEN)
   END IF
ELSEIF (NHKTA .EQ. 3) THEN !d orbitals
   IF (SPHRA) THEN
      IF (GMAX .GT. DFTHRI) THEN
         DO K = 1, NVCLEN
            PA_1 = PA(1,K)
            PA_2 = PA(2,K)
            PA_3 = PA(3,K)
            GAX  = PA_1*GA(K)
            GAY  = PA_2*GA(K)
            GAZ  = PA_3*GA(K)
            GAXX = PA_1*GAX
            GAYY = PA_2*GAY
            GAO(K,ISTART1,1) = CSP(1,2)*PA_2*GAX
            GAO(K,ISTART2,1) = CSP(2,5)*PA_2*GAZ
            GAO(K,ISTART3,1) = CSP(3,1)*GAXX + CSP(3,4)*GAYY&
            &                  + CSP(3,6)*PA_3*GAZ
            GAO(K,ISTART4,1) = CSP(4,3)*PA_1*GAZ
            GAO(K,ISTART5,1) = CSP(5,1)*GAXX + CSP(5,4)*GAYY
         END DO
      ELSE 
         CALL LS_DZERO(GAO(1,ISTART1,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART2,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART3,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART4,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART5,1),NVCLEN)
      END IF
   ELSE
      IF (GMAX .GT. DFTHRI) THEN
         DO K = 1, NVCLEN
            PA_1 = PA(1,K)
            PA_2 = PA(2,K)
            PA_3 = PA(3,K)
            GAX  = PA_1*GA(K)
            GAY  = PA_2*GA(K)
            GAZ  = PA_3*GA(K)
            GAO(K,ISTART1,1) = PA_1*GAX
            GAO(K,ISTART2,1) = PA_2*GAX
            GAO(K,ISTART3,1) = PA_3*GAX
            GAO(K,ISTART4,1) = PA_2*GAY
            GAO(K,ISTART5,1) = PA_3*GAY 
            GAO(K,ISTART6,1) = PA_3*GAZ
         END DO
      ELSE
         CALL LS_DZERO(GAO(1,ISTART1,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART2,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART3,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART4,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART5,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART6,1),NVCLEN)
      END IF
   END IF
ELSE !higher than d orbitals
IF (SPHRA) THEN
   DO I = 1, KHKTA
      CALL LS_DZERO(GAO(1,ISTART+I,1),NVCLEN)
   END DO
   IF (GMAX .GT. DFTHRI) THEN
      DO J = 1, KCKTA
         LVALJ = LVALUE(J)
         MVALJ = MVALUE(J)
         NVALJ = NVALUE(J)
         DO K = 1, NVCLEN
            CINT(K) = (PA(1,K)**LVALJ)*(PA(2,K)**MVALJ)&
            &                 *(PA(3,K)**NVALJ)*GA(K)
         END DO
         !              do a dgemm here?
         DO I = 1, KHKTA
            TEMPI=ISTART+I
            DO K = 1, NVCLEN
               GAO(K,TEMPI,1) = GAO(K,TEMPI,1) + CSP(I,J)*CINT(K)
            END DO
         END DO
      END DO
   END IF
ELSE
   IF (GMAX .GT. DFTHRI) THEN
      DO I = 1, KHKTA
         LVALI = LVALUE(I)
         MVALI = MVALUE(I)
         NVALI = NVALUE(I)
         TEMPI=ISTART+I
         DO K = 1, NVCLEN
            GAO(K,TEMPI,1) = (PA(1,K)**LVALI)*(PA(2,K)**MVALI)&
            &                 *(PA(3,K)**NVALI)*GA(K)
         END DO
      END DO
   END IF
END IF
END IF

END SUBROUTINE II_BLGETGAO

!> \brief evaluates the GAO(gaussian atomic orbitals) + first geo derivatives
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_BLGETGA1(LUPRI,NVCLEN,NactBAS,NTYPSO,NHKTA,KHKTA,KCKTA,MXPRIM,JSTA,&
                     & CSP,GAO,IADR,PA,PA2,DFTHRI,CC,CCSTART,PRIEXP,LVALUE,MVALUE,&
                     & NVALUE,SPVAL,SPHRA)
 
IMPLICIT NONE
INTEGER                   :: NVCLEN,NactBAS,MXPRIM,KHKTA,JSTA,K,KCKTA,LUPRI,NHKTA
INTEGER                   :: SPVAL,IADR,ISTART,NTYPSO
REAL(REALK),PARAMETER     :: D0 = 0.0D0, D1 = 1.0D0,D2 = 2.0D0
REAL(REALK),intent(inout) :: GAO(NVCLEN,NactBAS,NTYPSO) !NTYPSO=4
REAL(REALK),intent(in)    :: PA(3,NVCLEN),PA2(NVCLEN)
REAL(REALK)               :: CSP(KHKTA,KCKTA),PRIEXP(MXPRIM)
REAL(REALK)               :: GA(NVCLEN),VA(NVCLEN),VE(NVCLEN)
REAL(REALK)               :: DFTHRI
INTEGER                   :: I,J,LVALUE(SPVAL),MVALUE(SPVAL),NVALUE(SPVAL)
INTEGER,PARAMETER         :: I0=1, IX=2, IY=3, IZ=4
REAL(REALK)  :: GU(NVCLEN),GAX(NVCLEN),GAY(NVCLEN),GAZ(NVCLEN), P0(NVCLEN)
REAL(REALK)  :: FX(NVCLEN), FY(NVCLEN), FZ(NVCLEN),SPH0(NVCLEN), SPHX(NVCLEN)
REAL(REALK)  :: SPHY(NVCLEN), SPHZ(NVCLEN),FAC,TGX,TGY,TGZ
LOGICAL      :: SPHRA
REAL(REALK), pointer  :: CAO(:,:), CAOX(:,:)
REAL(REALK), pointer  :: CAOY(:,:),CAOZ(:,:)
REAL(REALK)  :: SPHFAC
REAL(REALK)  :: PRICCFVAL,PRIEXPVAL,GMAX,PA_1,PA_2,PA_3,GAVAL,DL,DN,DM
INTEGER      :: ICOMPA,L,M,N,CCSTART,TEMPI
INTEGER      :: ISTART1,ISTART2,ISTART3,ISTART4,ISTART5
TYPE(LSMATRIX) :: CC

!     loop over primitives
ISTART=IADR
ISTART1=ISTART+1
ISTART2=ISTART+2
ISTART3=ISTART+3
ISTART4=ISTART+4
ISTART5=ISTART+5
!CALL LS_DZERO(GA,NVCLEN)
!CALL LS_DZERO(GU,NVCLEN)
I=1
J = JSTA+CCSTART
PRICCFVAL = CC%elms(I)
PRIEXPVAL = PRIEXP(J)
DO K = 1, NVCLEN
   FAC = PRICCFVAL*DEXP(-PRIEXPVAL*PA2(K))
   GA(K) = FAC
   GU(K) = - D2*PRIEXPVAL*FAC
END DO
DO I = 2,CC%nrow
   J = JSTA+CCSTART-1+I
   PRICCFVAL = CC%elms(I)
   PRIEXPVAL = PRIEXP(J)
   DO K = 1, NVCLEN
      FAC   = PRICCFVAL*DEXP(-PRIEXPVAL*PA2(K))
      GA(K) = GA(K) + FAC
      GU(K) = GU(K) - D2*PRIEXPVAL*FAC
   END DO
END DO

IF (SPHRA) THEN
   call mem_alloc(CAO,NVCLEN,KCKTA)
   call mem_alloc(CAOX,NVCLEN,KCKTA)
   call mem_alloc(CAOY,NVCLEN,KCKTA)
   call mem_alloc(CAOZ,NVCLEN,KCKTA)
END IF

!screening based on the maximum GA value for the batch of points: GMAX
GMAX = D0
DO K = 1, NVCLEN
   GMAX = MAX(GMAX,DABS(GA(K)))
END DO

!s orbitals
IF (NHKTA .EQ. 1) THEN
   IF (GMAX .GT. DFTHRI) THEN
      DO K = 1, NVCLEN
         GAO(K,ISTART1,I0) = GA(K)
      END DO
      DO K = 1, NVCLEN
         GAO(K,ISTART1,IX) = PA(1,K)*GU(K) !CHANGE TO PA(K,1-3)
      END DO
      DO K = 1, NVCLEN
         GAO(K,ISTART1,IY) = PA(2,K)*GU(K)
      END DO
      DO K = 1, NVCLEN
         GAO(K,ISTART1,IZ) = PA(3,K)*GU(K)
      END DO  
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1,I0),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IX),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IY),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IZ),NVCLEN)
   END IF
!p orbitals
ELSE IF (NHKTA .EQ. 2) THEN
   IF (GMAX .GT. DFTHRI) THEN
      DO K = 1, NVCLEN
         PA_1 = PA(1,K)
         PA_2 = PA(2,K)
         PA_3 = PA(3,K)
         GAVAL = GA(K)
         TGX = PA_1*GU(K)     !CHANGE TO TGX(K)
         TGY = PA_2*GU(K)
         TGZ = PA_3*GU(K)
         GAO(K,ISTART1,I0) = PA_1*GAVAL
         GAO(K,ISTART2,I0) = PA_2*GAVAL
         GAO(K,ISTART3,I0) = PA_3*GAVAL
         GAO(K,ISTART1,IX) = PA_1*TGX + GAVAL
         GAO(K,ISTART2,IX) = PA_2*TGX
         GAO(K,ISTART3,IX) = PA_3*TGX
         GAO(K,ISTART1,IY) = PA_1*TGY
         GAO(K,ISTART2,IY) = PA_2*TGY + GAVAL
         GAO(K,ISTART3,IY) = PA_3*TGY
         GAO(K,ISTART1,IZ) = PA_1*TGZ
         GAO(K,ISTART2,IZ) = PA_2*TGZ
         GAO(K,ISTART3,IZ) = PA_3*TGZ + GAVAL
      END DO
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1,I0),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,I0),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,I0),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IX),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,IX),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,IX),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IY),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,IY),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,IY),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IZ),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,IZ),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,IZ),NVCLEN)
   END IF
! d and higher orbitals
ELSE
   IF (GMAX.GT.DFTHRI) THEN
      DO K = 1, NVCLEN
         FX(K) = PA(1,K)*GU(K)!CHANGE TO PA(K,1-3)
         FY(K) = PA(2,K)*GU(K)
         FZ(K) = PA(3,K)*GU(K)
      END DO
      DO ICOMPA = 1,KCKTA
         L = LVALUE(ICOMPA)
         M = MVALUE(ICOMPA)
         N = NVALUE(ICOMPA)
         DO K = 1, NVCLEN
            P0(K)  = (PA(1,K)**L)*(PA(2,K)**M)*(PA(3,K)**N) ! make different loops for L=0 or M=0 or N=0 etc ??
         END DO
         DO K = 1, NVCLEN
            GAX(K) = FX(K)*P0(K)
            GAY(K) = FY(K)*P0(K)
            GAZ(K) = FZ(K)*P0(K)
            GAVAL = GA(K)
            PA_1  = PA(1,K)
            PA_2  = PA(2,K)
            PA_3  = PA(3,K)
            IF(L.GT.0) THEN
               GAX(K) = GAX(K)+L*(PA_1**(L-1))*(PA_2**M)*&
                    &  (PA_3**N)*GAVAL
            ENDIF
            IF(M.GT.0) THEN
               GAY(K) = GAY(K)+M*(PA_1**L)*(PA_2**(M-1))*&
                    &   (PA_3**N)*GAVAL
            ENDIF
            IF(N.GT.0) THEN
               GAZ(K) = GAZ(K)+N*(PA_1**L)*(PA_2**M)*&
                    &  (PA_3**(N-1))*GAVAL
            ENDIF
         END DO
         IF (SPHRA) THEN
            DO K = 1, NVCLEN
               CAO (K, ICOMPA) = GA(K)*P0(K)
               CAOX(K, ICOMPA) = GAX(K)
               CAOY(K, ICOMPA) = GAY(K)
               CAOZ(K, ICOMPA) = GAZ(K)
            END DO
         ELSE
            TEMPI = ISTART+ICOMPA
            DO K = 1, NVCLEN
               GAO(K,TEMPI,I0) = GA(K)*P0(K)
            END DO
            DO K = 1, NVCLEN
               GAO(K,TEMPI,IX) = GAX(K)
            END DO
            DO K = 1, NVCLEN
               GAO(K,TEMPI,IY) = GAY(K)
            END DO
            DO K = 1, NVCLEN
               GAO(K,TEMPI,IZ) = GAZ(K)
            END DO
         END IF
      END DO
      IF(SPHRA) THEN
         IF (NHKTA.EQ.3) THEN !d orbitals
            DO K = 1, NVCLEN
               GAO(K,ISTART1,I0) = CSP(1,2)*CAO(K,2)
               GAO(K,ISTART2,I0) = CSP(2,5)*CAO(K,5)
               GAO(K,ISTART3,I0) = CSP(3,1)*CAO(K,1) + CSP(3,4)*CAO(K,4)& 
               &                               + CSP(3,6)*CAO(K,6)
               GAO(K,ISTART4,I0) = CSP(4,3)*CAO(K,3)
               GAO(K,ISTART5,I0) = CSP(5,1)*CAO(K,1) + CSP(5,4)*CAO(K,4)
               GAO(K,ISTART1,IX) = CSP(1,2)*CAOX(K,2)
               GAO(K,ISTART2,IX) = CSP(2,5)*CAOX(K,5)
               GAO(K,ISTART3,IX) = CSP(3,1)*CAOX(K,1) + CSP(3,4)*CAOX(K,4)& 
               &                                + CSP(3,6)*CAOX(K,6)
               GAO(K,ISTART4,IX) = CSP(4,3)*CAOX(K,3)
               GAO(K,ISTART5,IX) = CSP(5,1)*CAOX(K,1) + CSP(5,4)*CAOX(K,4)
               GAO(K,ISTART1,IY) = CSP(1,2)*CAOY(K,2)
               GAO(K,ISTART2,IY) = CSP(2,5)*CAOY(K,5)
               GAO(K,ISTART3,IY) = CSP(3,1)*CAOY(K,1) + CSP(3,4)*CAOY(K,4)&
               &                                + CSP(3,6)*CAOY(K,6)
               GAO(K,ISTART4,IY) = CSP(4,3)*CAOY(K,3)
               GAO(K,ISTART5,IY) = CSP(5,1)*CAOY(K,1) + CSP(5,4)*CAOY(K,4)
               GAO(K,ISTART1,IZ) = CSP(1,2)*CAOZ(K,2)
               GAO(K,ISTART2,IZ) = CSP(2,5)*CAOZ(K,5)
               GAO(K,ISTART3,IZ) = CSP(3,1)*CAOZ(K,1) + CSP(3,4)*CAOZ(K,4)& 
               &                                + CSP(3,6)*CAOZ(K,6)
               GAO(K,ISTART4,IZ) = CSP(4,3)*CAOZ(K,3)
               GAO(K,ISTART5,IZ) = CSP(5,1)*CAOZ(K,1) + CSP(5,4)*CAOZ(K,4)
            END DO
         ELSE
            ! FIXME: take timings without use of the temporary arrays
            DO I = 1, KHKTA
               DO K = 1, NVCLEN
                  SPH0(K) = D0 
                  SPHX(K) = D0 
                  SPHY(K) = D0 
                  SPHZ(K) = D0 
               END DO
               DO J = 1, KCKTA
                  SPHFAC = CSP(I,J)
                  IF (DABS(SPHFAC).GT.D0) THEN
                     DO K = 1, NVCLEN
                        SPH0(K) = SPH0(K) + SPHFAC*CAO (K,J)
                        SPHX(K) = SPHX(K) + SPHFAC*CAOX(K,J)
                        SPHY(K) = SPHY(K) + SPHFAC*CAOY(K,J)
                        SPHZ(K) = SPHZ(K) + SPHFAC*CAOZ(K,J)
                     END DO
                  END IF
               END DO
               TEMPI = ISTART+I
               DO K = 1, NVCLEN
                  GAO(K,TEMPI,I0) = SPH0(K)
                  GAO(K,TEMPI,IX) = SPHX(K)
                  GAO(K,TEMPI,IY) = SPHY(K)
                  GAO(K,TEMPI,IZ) = SPHZ(K)
               END DO
            END DO
         END IF
      END IF
   ELSE
     IF (.NOT.SPHRA) THEN
         DO I = 1, KCKTA
            TEMPI = ISTART+I
            CALL LS_DZERO(GAO(1,TEMPI,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,TEMPI,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,TEMPI,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,TEMPI,IZ),NVCLEN)
         END DO
      ELSE
         IF (NHKTA.EQ.3) THEN
            CALL LS_DZERO(GAO(1,ISTART1,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART2,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART3,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART4,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART5,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART1,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART2,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART3,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART4,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART5,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART1,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART2,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART3,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART4,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART5,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART1,IZ),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART2,IZ),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART3,IZ),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART4,IZ),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART5,IZ),NVCLEN)
         ELSE
            DO I = 1, KHKTA
               TEMPI = ISTART+I
               CALL LS_DZERO(GAO(1,TEMPI,I0),NVCLEN)
               CALL LS_DZERO(GAO(1,TEMPI,IX),NVCLEN)
               CALL LS_DZERO(GAO(1,TEMPI,IY),NVCLEN)
               CALL LS_DZERO(GAO(1,TEMPI,IZ),NVCLEN)
            END DO
         END IF
      END IF
   END IF
END IF

IF (SPHRA) THEN
   call mem_dealloc(CAO)
   call mem_dealloc(CAOX)
   call mem_dealloc(CAOY)
   call mem_dealloc(CAOZ)
END IF

END SUBROUTINE II_BLGETGA1

!> \brief evaluates the GAO(gaussian atomic orbitals) + first and second geo derivatives
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_BLGETGA2(LUPRI,NVCLEN,NactBAS,NHKTA,KHKTA,KCKTA,MXPRIM,JSTA,&
                     & CSP,GAOXX,GAOXY,GAOXZ,GAOYY,GAOYZ,GAOZZ,IADR,PA,PA2,&
                     &DFTHRI,CC,CCSTART,PRIEXP,LVALUETK,MVALUETK,NVALUETK,SPVAL,SPHRA)
IMPLICIT NONE
INTEGER                   :: NVCLEN,MXPRIM,KHKTA,JSTA,K,KCKTA,LUPRI,NHKTA,NactBAS
INTEGER                   :: SPVAL,CCSTART,IADR,ISTART
REAL(REALK),PARAMETER     :: D0 = 0.0D0, D1 = 1.0D0,D2 = 2.0D0
REAL(REALK),intent(inout) :: GAOXX(NVCLEN,NactBAS), GAOXY(NVCLEN,NactBAS)
REAL(REALK),intent(inout) :: GAOXZ(NVCLEN,NactBAS), GAOYY(NVCLEN,NactBAS)
REAL(REALK),intent(inout) :: GAOYZ(NVCLEN,NactBAS), GAOZZ(NVCLEN,NactBAS)
REAL(REALK),intent(in)    :: PA(3,NVCLEN),PA2(NVCLEN)
REAL(REALK)               :: CSP(KHKTA,KCKTA),PRIEXP(MXPRIM)
REAL(REALK)               :: GA,VA(NVCLEN),VE(NVCLEN)
REAL(REALK)               :: DFTHRI,GAZ,GAX,GAY
INTEGER                   :: I,J,LVALUETK(SPVAL),MVALUETK(SPVAL),NVALUETK(SPVAL)
INTEGER,PARAMETER         :: I0=1, IX=2, IY=3, IZ=4
LOGICAL      :: SPHRA
REAL(REALK)  :: CAOXX(KCKTA),CAOXY(KCKTA),CAOXZ(KCKTA)
REAL(REALK)  :: CAOYY(KCKTA),CAOYZ(KCKTA),CAOZZ(KCKTA)
REAL(REALK)  :: ALPHA,TALPH,TAPAX,TAPAY,TAPAZ,PXD,PYD,PZD,PXM,PYM,PZM
REAL(REALK)  :: PX0,PY0,PZ0,PXP,PYP,PZP,P000,GAXX,GAXY,GAXZ,GAYY,GAYZ,GAZZ
REAL(REALK)  :: SPHXX,SPHXY,SPHXZ,SPHYY,SPHYZ,SPHZZ,SPHFAC,PA_1,PA_2,PA_3,PA2_B
INTEGER      :: L,M,N,ICOMPA,IPRIMA,IV,TEMPI
TYPE(LSMATRIX) :: CC
!CHANGE THIS ACCORDING TO NEW BLGETGA2 CHANGED BY ANDREAS
ISTART=IADR
DO IV = 1, NVCLEN
   PA_1=PA(1,IV)
   PA_2=PA(2,IV)
   PA_3=PA(3,IV)
   PA2_B=PA2(IV)
   IF (SPHRA) THEN
      DO I=1,KCKTA
         CAOXX(I) = D0
         CAOXY(I) = D0
         CAOXZ(I) = D0
         CAOYY(I) = D0
         CAOYZ(I) = D0
         CAOZZ(I) = D0
      END DO
   END IF
   DO I = 1,CC%nrow
      IPRIMA = JSTA+CCSTART-1+I
      ALPHA = PRIEXP(IPRIMA)
      TALPH = -D2*ALPHA
      TAPAX = TALPH*PA_1
      TAPAY = TALPH*PA_2
      TAPAZ = TALPH*PA_3
      GA = CC%elms(I)*DEXP(-ALPHA*PA2(IV))!FAC
      IF (ABS(GA).GT.DFTHRI) THEN
         DO ICOMPA = 1, KCKTA
            L = LVALUETK(ICOMPA)
            M = MVALUETK(ICOMPA)
            N = NVALUETK(ICOMPA)
!                            
            PXD = D0
            PYD = D0
            PZD = D0
            IF (L.GT.1) PXD = (L*(L-1))*(PA_1**(L-2))
            IF (M.GT.1) PYD = (M*(M-1))*(PA_2**(M-2))
            IF (N.GT.1) PZD = (N*(N-1))*(PA_3**(N-2))
            PXM = D0
            PYM = D0
            PZM = D0
            IF (L.GT.0) PXM = (L)*(PA_1**(L-1))
            IF (M.GT.0) PYM = (M)*(PA_2**(M-1))
            IF (N.GT.0) PZM = (N)*(PA_3**(N-1))
            PX0 = PA_1**L
            PY0 = PA_2**M
            PZ0 = PA_3**N
            PXP = TAPAX*PX0
            PYP = TAPAY*PY0
            PZP = TAPAZ*PZ0
            P000 = PX0*PY0*PZ0 
            IF (NHKTA.EQ.1) THEN
               ! s orbitals
               GAXX = TAPAX**2 + TALPH
               GAYY = TAPAY**2 + TALPH
               GAZZ = TAPAZ**2 + TALPH
               GAXY = TAPAX*TAPAY
               GAXZ = TAPAX*TAPAZ
               GAYZ = TAPAY*TAPAZ
            ELSE IF (NHKTA.EQ.2) THEN
               ! p orbitals
               GAXX = (TAPAX**2 + TALPH*(2*L+1))*P000
               GAYY = (TAPAY**2 + TALPH*(2*M+1))*P000
               GAZZ = (TAPAZ**2 + TALPH*(2*N+1))*P000
               GAXY = TAPAX*TAPAY*P000+(PXP*PYM+PXM*PYP)*PZ0
               GAXZ = TAPAX*TAPAZ*P000+(PXP*PZM+PXM*PZP)*PY0
               GAYZ = TAPAY*TAPAZ*P000+(PYP*PZM+PYM*PZP)*PX0
            ELSE 
               ! d and higher orbitals
               GAXX = (TAPAX**2 + TALPH*(2*L+1))*P000 + PXD*PY0*PZ0
               GAYY = (TAPAY**2 + TALPH*(2*M+1))*P000 + PX0*PYD*PZ0
               GAZZ = (TAPAZ**2 + TALPH*(2*N+1))*P000 + PX0*PY0*PZD
               GAXY = TAPAX*TAPAY*P000 + (PXP*PYM+PXM*PYP+PXM*PYM)*PZ0
               GAXZ = TAPAX*TAPAZ*P000 + (PXP*PZM+PXM*PZP+PXM*PZM)*PY0
               GAYZ = TAPAY*TAPAZ*P000 + (PYP*PZM+PYM*PZP+PYM*PZM)*PX0
            END IF
            IF (SPHRA.AND.NHKTA.GT.2) THEN
               CAOXX(ICOMPA) = CAOXX(ICOMPA) + GAXX*GA 
               CAOXY(ICOMPA) = CAOXY(ICOMPA) + GAXY*GA 
               CAOXZ(ICOMPA) = CAOXZ(ICOMPA) + GAXZ*GA
               CAOYY(ICOMPA) = CAOYY(ICOMPA) + GAYY*GA
               CAOYZ(ICOMPA) = CAOYZ(ICOMPA) + GAYZ*GA
               CAOZZ(ICOMPA) = CAOZZ(ICOMPA) + GAZZ*GA
            ELSE
               TEMPI = ISTART+ICOMPA
               GAOXX(IV,TEMPI) = GAOXX(IV,TEMPI) + GAXX*GA 
               GAOXY(IV,TEMPI) = GAOXY(IV,TEMPI) + GAXY*GA 
               GAOXZ(IV,TEMPI) = GAOXZ(IV,TEMPI) + GAXZ*GA
               GAOYY(IV,TEMPI) = GAOYY(IV,TEMPI) + GAYY*GA
               GAOYZ(IV,TEMPI) = GAOYZ(IV,TEMPI) + GAYZ*GA
               GAOZZ(IV,TEMPI) = GAOZZ(IV,TEMPI) + GAZZ*GA
            ENDIF
         END DO
      END IF
   END DO
   IF (SPHRA.AND.NHKTA.GT.2) THEN
      DO I = 1,KHKTA
         SPHXX = D0 
         SPHXY = D0 
         SPHXZ = D0 
         SPHYY = D0 
         SPHYZ = D0 
         SPHZZ = D0 
         DO J = 1, KCKTA
            SPHFAC = CSP(I,J)
            IF (ABS(SPHFAC).GT.D0) THEN
               SPHXX = SPHXX + SPHFAC*CAOXX(J)
               SPHXY = SPHXY + SPHFAC*CAOXY(J)
               SPHXZ = SPHXZ + SPHFAC*CAOXZ(J)
               SPHYY = SPHYY + SPHFAC*CAOYY(J)
               SPHYZ = SPHYZ + SPHFAC*CAOYZ(J)
               SPHZZ = SPHZZ + SPHFAC*CAOZZ(J)
            END IF
         END DO
         TEMPI = ISTART+I
         GAOXX(IV,TEMPI) = SPHXX
         GAOXY(IV,TEMPI) = SPHXY
         GAOXZ(IV,TEMPI) = SPHXZ
         GAOYY(IV,TEMPI) = SPHYY
         GAOYZ(IV,TEMPI) = SPHYZ
         GAOZZ(IV,TEMPI) = SPHZZ
      END DO
   END IF
END DO

END SUBROUTINE II_BLGETGA2

!> \brief evaluates the london derivative GAOs (gaussian atomic orbitals)
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_BLGETGB1(NVCLEN,NactBAS,IADR,KHKTA,GAO,PA,GABX,GABY,GABZ)
IMPLICIT NONE
INTEGER :: NVCLEN,NactBAS,IADR,KHKTA
REAL(REALK) :: GAO(NVCLEN,NactBAS),GABX(NVCLEN,NactBAS),GABY(NVCLEN,NactBAS)
REAL(REALK) :: GABZ(NVCLEN,NactBAS),PA(3,NVCLEN)
!
REAL(REALK),PARAMETER :: D05 = 0.5d0
INTEGER :: I,J,K
REAL(REALK) :: GA

DO I = 1,KHKTA
 J = IADR+I
 DO K = 1, NVCLEN
    GA = D05*GAO(K,J)
    GABX(K,J) = GA*PA(1,K)
    GABY(K,J) = GA*PA(2,K)
    GABZ(K,J) = GA*PA(3,K)
 END DO
END DO

END SUBROUTINE II_BLGETGB1

!> \brief make blocks that contain the orbitals which contribute to each gridbatch
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE SHELL_TO_ORB(LUPRI,NSHELL,SHELBLOCK,ORBBLOCKS,NSTART,KMAX,NBAST)
!SHELL_TO_ORB transforms shell block indices to orbital block indices.
IMPLICIT NONE
INTEGER  :: NSHELL,KMAX,NBAST,LUPRI
INTEGER  :: SHELBLOCK(2,KMAX),ORBBLOCKS(2,KMAX),NSTART(KMAX)
!
INTEGER  :: ISHELL,I

DO I = 1,NSHELL
   ORBBLOCKS(1,I) = NSTART(SHELBLOCK(1,I))+1
ENDDO

DO I = 1,NSHELL
   IF(SHELBLOCK(2,I) .LT. KMAX)THEN
      ORBBLOCKS(2,I) = NSTART(SHELBLOCK(2,I)+1)
   ELSE
      ORBBLOCKS(2,I) = NBAST
   ENDIF
!   WRITE(LUPRI,*) '("shell ",2I4," translated to orbital ",2I4)'&
!    &,SHELBLOCK(1,I),SHELBLOCK(2,I),ORBBLOCKS(1,I),ORBBLOCKS(2,I)
ENDDO

END SUBROUTINE SHELL_TO_ORB

!> \brief make blocks that contain the active orbitals which contribute to each gridbatch
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE SHELL_TO_ACTORB(NSHELL,BLOCKS,KMAX)
!SHELL_TO_ORB transforms shell block indices to orbital block indices.
IMPLICIT NONE
INTEGER  :: NSHELL,KMAX
INTEGER  :: BLOCKS(2,KMAX)
!
INTEGER  :: IBASIS,IBL,IORB1,IORB2

IBASIS = 0
DO IBL = 1, NSHELL
   IBASIS = IBASIS + 1
   IORB1 = BLOCKS(1,IBL)
   BLOCKS(1,IBL) = IBASIS
   IORB2 = BLOCKS(2,IBL) !ORBITALINDEX
   IBASIS = IBASIS + (IORB2-IORB1)
   BLOCKS(2,IBL) = IBASIS
ENDDO

END SUBROUTINE SHELL_TO_ACTORB

!> \brief build spherical transformation matrices
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE Build_PRECALCULATED_SPHMAT(LUPRI,MAXANGMOM,SIZE,SPHMAT,SPINDEX) 
use math_fun
IMPLICIT NONE
INTEGER          :: MAXANGMOM,nAngmom,SIZE,LUPRI,SPINDEX(MAXANGMOM+1)
REAL(REALK)      :: SPHMAT(SIZE)
INTEGER          :: L,I
Real(realk), parameter :: DM1 = -1.0d0, DO = 0.0d0, D1 = 1.0d0, D2 = 2.0d0
INTEGER  :: M1,MADR,MABS,V0, NDER, IOFF
INTEGER  :: EE,FF,GG,BX,BY,BZ,II,JJ,KK,PX,PY,PZ
REAL(realk)  :: FACNRM,FAC1,FAC2, FAC3,FAC4, FACTOR
INTEGER  :: T,U,V,A,B,C,P,Q,R,X,Y,Z,TOT,IADR,AX0,AY0,AZ0,NSIZE
INTEGER  :: M,Ncol,Nrow,INDEX,NDIM,STARTINDEX

IF(MAXANGMOM .LE. 1)CALL LSQUIT('ERROR IN Build_PRECALCULATED_SPHMAT',lupri)
! CALL LS_DZERO(SPHMAT,SIZE) SHOULD BE DONE OUTSIDE

nANGMOM=MAXANGMOM+1 
NSIZE=0
STARTINDEX=1
SPINDEX(1)=1
SPINDEX(2)=1
DO I=3,nANGMOM
   SPINDEX(I)=STARTINDEX
   L = I-1 !angmom
   NRow = 2*L+1
   NCol = (L+1)*(L+2)/2
DO M1 = 0, 2*L 
   M = M1 - L
   IF (L.EQ.1) THEN
      IF (M .EQ. -1) MADR =  0  
      IF (M .EQ.  0) MADR =  1 
      IF (M .EQ.  1) MADR = -1 
   ELSE
      MADR = M
   END IF
   MABS = ABS(M)
   V0 = 0
   IF (M .LT. 0) V0 = 1 
   FACNRM = D1
   IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(LUPRI,L+MABS)*&
        &FACULT(LUPRI,L-MABS))/(FACULT(LUPRI,L)*(D2**MABS))
   FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
   FACNRM = FACNRM/SQRT(FACUL2(LUPRI,2*L-1))
   DO T = 0, L - MABS, 2
   DO U = 0, T, 2
   DO V = V0, MABS, 2
      !        almost 6.4.48 in the book
      FAC3 = FACNRM*BINOM(LUPRI,L,T/2)*BINOM(LUPRI,L-T/2,MABS+T/2)&
           &                    *BINOM(LUPRI,T/2,U/2)*BINOM(LUPRI,MABS,V)
      DO A = 0, MIN(0,T+MABS-U-V) 
      DO B = 0, MIN(0,U+V)
      DO C = 0, MIN(0,L-T-MABS)
         !           6.4.47 in the book
         DO P = 0, - A, 2
         DO Q = 0, - B, 2
         DO R = 0, - C, 2
            FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
                 &   D2**(-A-B-C-P-Q-R-T)*FAC3
            X = T+MABS-U-V-2*A-P
            Y = U+V-2*B-Q
            Z = L-T-MABS-2*C-R
            TOT = X + Y + Z
            IADR = 1 + (2*L+1)*(NCRT(X,Y,Z)-1) + L + MADR+NSIZE
            SPHMAT(IADR) = SPHMAT(IADR) + FACTOR 
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
   ENDDO
   ENDDO
   ENDDO
ENDDO
   NSIZE= NSIZE+Nrow*Ncol 
   STARTINDEX=STARTINDEX+Nrow*Ncol
ENDDO
END SUBROUTINE BUILD_PRECALCULATED_SPHMAT

!
!     DFTDISP calculates the empirical dispersion correction to the energy and the gradient 
!     as formulated by Grimme in J. Comp. Chem. 2004, 25, 1463. and  J. Comp. Chem. 2006, 27, 1787.
!
!     EDISP = -S6 * sum_i^{N-1}sum_{j=i+1}^N(C_6(ij)/R(ij)^6 * f(R_ij))
!
!     with the damping function f(R_ij) = 1/(1+exp(-d*[R_ij/Rr-1]))
!
!     where N is the number of atoms
!           R_ij is the interatomic distance
!           Rr is the sum of the van-der-Waals radii
!           d is a fixed damping parameter (d=20.0)
!           C_6(ij) = sqrt(C_6(i)*C_6(j)
!           C_6 are fixed, atomic C_6 parameters
!           S_6 is a fixed, functional dependent global scaling factor (defined for BP86, BLYP, PBE, B3LYP, TPSS)
!
!     -------------------------------------------------------------
!     EDISP:  contains the energy correction at output
!     NDERIV: order of derivative wrt nuclei 
!             0: energy only
!             1: first derivative
!             higher: not defined yet
!     -------------------------------------------------------------
!
!     03.2010 Andreas Krapp 
!
SUBROUTINE II_DFTDISP(SETTING,GRAD,DIM1,DIM2,NDERIV,LUPRI,IPRINT)
use ls_util
use ks_settings
IMPLICIT NONE
! external integer
INTEGER, INTENT(IN) :: NDERIV, LUPRI, DIM1, DIM2, IPRINT
! internal integer
INTEGER :: NCENTA, NCENTB, ISCOOA, ISCOOB, ISCOOR, IATOM, IOFF, J, NATOMS
! external real
REAL(REALK), INTENT(INOUT) :: GRAD(DIM1,DIM2)
! internal real
REAL(REALK) :: R0(54), C6(54)
REAL(REALK) :: E, EADD, S6
REAL(REALK) :: CHARGA, CORDAX, CORDAY, CORDAZ, C6A, RvdWA
REAL(REALK) :: CHARGB, CORDBX, CORDBY, CORDBZ, C6B, RvdWB
REAL(REALK) :: RX, RY, RZ, R2, R, RR, R6FAC, C6FAC
REAL(REALK) :: ALPHA, EXPOA, FDMP
REAL(REALK) :: DFAC, GRADX, GRADY, GRADZ
REAL(REALK), ALLOCATABLE :: GRDFT(:)
REAL, PARAMETER ::  CONV1=1.88972612   ! convert angstrom to bohr
REAL, PARAMETER ::  CONV2=17.3452771   ! convert Joule*nm^6/mol to Bohr^6*hartree
                                       !    1 Bohr = 52.9177 * 10^-3 nm
                                       !    1 Hartree = 2.6255*10^6 Joule/mol
!  external function
!EXTERNAL DISP_FUNCFAC
!  external types
TYPE(LSSETTING),   INTENT(INOUT) :: SETTING
!
! van der Waals radii for the elements H-Xe in Angstrom
! taken from JCC 2006, 27, 1787
!
     DATA &
!         H     He
     & R0/1.001,1.012, &
!         Li    Be    B     C     N     O     F     Ne
     &    0.825,1.408,1.485,1.452,1.397,1.342,1.287,1.243, &
!         Na    Mg    Al    Si    P     S     Cl    Ar
     &    1.144,1.364,1.639,1.716,1.705,1.683,1.639,1.595, &
!         K     Ca
     &    1.485,1.474, &
!         Sc-Zn,
     &    1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562, &
!         Ga    Ge    As    Se    Br    Kr
     &    1.650,1.727,1.760,1.771,1.749,1.727, &
!         Rb    Sr    
     &    1.628,1.606, &
!         Y-Cd,
     &    1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.639, &
!         In    Sn    Sb    Te    I     Xe
     &    1.672,1.804,1.881,1.892,1.892,1.881/ 

!
!     C6 parameters for the elements H-Xe in Joule*nm^6/mol
!     taken from JCC 2006, 27, 1787
!
      DATA &
!         H     He
     & C6/0.14 ,0.08 , &
!         Li    Be    B     C     N     O     F     Ne
     &    1.61 ,1.61 ,3.13 ,1.75 ,1.23 ,0.70 ,0.75 ,0.63 , &
!         Na    Mg    Al    Si    P     S     Cl    Ar
     &    5.71 ,5.71 ,10.79,9.23 ,7.84 ,5.57 ,5.07 ,4.61 , &
!         K     Ca
     &    10.80,10.80, &
!         Sc-Zn,
     &    10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80, &
!         Ga    Ge    As    Se    Br    Kr
     &    16.99,17.10,16.37,12.64,12.47,12.01, &
!         Rb    Sr    
     &    24.67,24.67, &
!         Y-Cd,
     &    24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67, &
!         In    Sn    Sb    Te    I     Xe
     &    37.32,38.71,38.44,31.74,31.50,29.99/


      ! only for DFT 
      IF (.NOT. SETTING%DO_DFT) RETURN

      ! we do not want dispersion correction -> return
      IF (.NOT. SETTING%SCHEME%DODISP) RETURN

      IF (NDERIV.EQ.1) THEN
         SETTING%SCHEME%DISPDONE = .FALSE.
      ELSE
         ! if we have calculated the dispersion correction we can return
         IF (SETTING%SCHEME%DISPDONE) THEN 
            IF (IPRINT .GE. 1) THEN
               WRITE(LUPRI,'(1X,A39,F24.10)') 'dispersion correction to the KS-energy: ', SETTING%EDISP
               WRITE(LUPRI,*)''
            END IF
            RETURN
         ENDIF
         IF (SETTING%MOLECULE(1)%p%NATOMS.LE.1) THEN 
            SETTING%SCHEME%DISPDONE = .FALSE.
            SETTING%EDISP = 0.0D0
            RETURN
         ELSE
            SETTING%SCHEME%DISPDONE = .TRUE.
         END IF
      END IF

!     error check
      IF (NDERIV.GT.1 .OR. NDERIV .LT. 0) THEN
         WRITE(LUPRI,'(4X,A62)') &
     &  ' WARNING: dispersion correction only for energies and gradients'
         WRITE(*,*) &
     &  'WARNING: dispersion correction only for energies and gradients'
      END IF

!     Print section
      WRITE(LUPRI,*)''
      WRITE(LUPRI,'(1X,A56)') &
     &   'Add empirical dispersion corr. to the XC-energy/gradient'
      WRITE(LUPRI,'(1X,A25)')'  following S. Grimme,   '
      WRITE(LUPRI,'(1X,A25)')'  JCC 2004, 25, 1463. and'
      WRITE(LUPRI,'(1X,A25)')'  JCC 2006, 27, 1787.    '

!     initialisations and memory for gradient
      NATOMS = SETTING%MOLECULE(1)%p%NATOMS
      E = 0.0D0
      SETTING%EDISP = 0.0D0
      ALLOCATE(GRDFT(DIM1*DIM2))
      IF (NDERIV.EQ.1) THEN
         CALL LS_DZERO(GRDFT,DIM1*DIM2)
      END IF

!     calculate correction

!
!     Run over nuclei A
!
      DO NCENTA = 1, NATOMS-1
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE
         IF (ABS(CHARGA-NINT(CHARGA)).GT.1d-10) THEN
            CALL LSQUIT('Error in DFTDISP. Not implemented for&
                       & non-integer charge!',-1)
         ENDIF
         IF (NINT(CHARGA) .LE. 54 .AND. NINT(CHARGA) .GT. 0) THEN
            ISCOOA = (NCENTA-1)*3
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(3)
            C6A    = C6(NINT(CHARGA))*CONV2
            RvdWA  = R0(NINT(CHARGA))*CONV1
!
!           Run over nuclei B
!
            DO NCENTB =  NCENTA+1, NATOMS
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE
               IF (ABS(CHARGB-NINT(CHARGB)).GT.1d-10) THEN
                  CALL LSQUIT('Error in DFTDISP. Not implemented for&
                             & non-integer charge!',-1)
               ENDIF
               IF (NINT(CHARGB).LE. 54 .AND.NINT(CHARGB).GT. 0) THEN
                  ISCOOB = (NCENTB-1)*3
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(3)
                  C6B    = C6(NINT(CHARGB))*CONV2
                  RvdWB  = R0(NINT(CHARGB))*CONV1

!                 C6 factor for atom pair A-B
                  C6FAC = SQRT(C6A*C6B)

!                 distance R between atoms A and B
!                 and R^6 
                  RX    = CORDAX-CORDBX
                  RY    = CORDAY-CORDBY
                  RZ    = CORDAZ-CORDBZ
                  R2    = RX**2 + RY**2 + RZ**2
                  R6FAC = R2*R2*R2
                  R     = SQRT(R2)

!                 sum of the van der Waals radii RR  
                  RR   = RvdWA + RvdWB

!                 damping function FDMP
                  ALPHA = -20.0*R/RR+20.0D0
                  EXPOA = EXP(ALPHA)
                  FDMP  = 1.0D0/(1.0D0 + EXPOA)

!                 dispersion correction contribution
                  EADD  = C6FAC/R6FAC * FDMP
                  E = E + EADD

                  IF (NDERIV.EQ.1) THEN
!                    derivative wrt to nuclear positions
                     DFAC = (-6.0/R + 20.0/RR * EXPOA * FDMP)*EADD / R
                     GRADX = DFAC * RX
                     GRADY = DFAC * RY
                     GRADZ = DFAC * RZ
                     GRDFT(ISCOOA + 1) = GRDFT( ISCOOA + 1) + GRADX
                     GRDFT(ISCOOA + 2) = GRDFT( ISCOOA + 2) + GRADY
                     GRDFT(ISCOOA + 3) = GRDFT( ISCOOA + 3) + GRADZ
                     GRDFT(ISCOOB + 1) = GRDFT( ISCOOB + 1) - GRADX
                     GRDFT(ISCOOB + 2) = GRDFT( ISCOOB + 2) - GRADY
                     GRDFT(ISCOOB + 3) = GRDFT( ISCOOB + 3) - GRADZ
                  END IF

               ELSE
                  WRITE(LUPRI,'(4X,A42)') 'DISPERSION CORRECTION ONLY FOR ATOMS 1-54.'
                  WRITE(LUPRI,'(4X,A30,I6,A4,I6)') 'ACTUAL CHARGE FOR NUCLEUS ',NCENTB,' IS ',NINT(CHARGB)
                  CALL LSQUIT('DISPERSION CORRECTION ONLY FOR ATOMS 1-54.',-1)
               END IF
            END DO
         ELSE 
             WRITE(LUPRI,'(4X,A42)') 'DISPERSION CORRECTION ONLY FOR ATOMS 1-54.'
             WRITE(LUPRI,'(4X,A30,I6,A4,I6)') 'ACTUAL CHARGE FOR NUCLEUS ',NCENTA,' IS ',NINT(CHARGA)
             CALL LSQUIT('DISPERSION CORRECTION ONLY FOR ATOMS 1-54.',-1)
         END IF
      END DO

!     final dispersion correction
!     S6 factor is functional dependant
      CALL DISP_FUNCFAC(S6)
      SETTING%EDISP = -S6 * E

!     Add derivative
      IF (NDERIV.EQ.1) THEN
         IF( (DIM1.NE.3) .OR. (DIM2.NE.NATOMS)) THEN 
            CALL LSQUIT('ERROR IN DFTDISP WITH GRADIENT',-1)
         ENDIF
         DO IATOM = 1, NATOMS
            ISCOOR = (IATOM-1)*3 
            GRAD(1,IATOM) = GRAD(1,IATOM) - GRDFT(ISCOOR+1)*S6
            GRAD(2,IATOM) = GRAD(2,IATOM) - GRDFT(ISCOOR+2)*S6
            GRAD(3,IATOM) = GRAD(3,IATOM) - GRDFT(ISCOOR+3)*S6
         END DO
      END IF

      DEALLOCATE(GRDFT)

!     Print section
      WRITE(LUPRI,'(1X,A13,F7.2)')'    S6 factor',S6
      WRITE(LUPRI,'(1X,A15)')     '    C6 factors:'
      DO IATOM = 1, NATOMS
         WRITE (LUPRI, '(2X,A6,F7.3)') &
     &          SETTING%MOLECULE(1)%p%ATOM(IATOM)%NAME, C6(NINT(SETTING%MOLECULE(1)%p%ATOM(IATOM)%CHARGE))
      END DO

      IF (IPRINT .GE. 1) THEN
         WRITE(LUPRI,'(1X,A39,F24.10)') 'dispersion correction to the KS-energy: ', SETTING%EDISP
         IF (NDERIV.EQ.1) THEN
            CALL LSHEADER(LUPRI,'XC-gradient including empir. disp. corr.')
            DO IATOM = 1, NATOMS
               WRITE (LUPRI, '(1X,A6,F17.10,2F24.10)') &
     &            SETTING%MOLECULE(1)%p%ATOM(IATOM)%NAME, (GRAD(J,IATOM),J=1,3)
            END DO
         END IF
         WRITE(LUPRI,*)''
      END IF

   RETURN
END SUBROUTINE II_DFTDISP

END MODULE IIDFTINT
