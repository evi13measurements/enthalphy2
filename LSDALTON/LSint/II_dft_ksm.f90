!> @file
!> Module containing soubroutine for calculation of the exchange-correlation contribution to KS-matrix
MODULE IIDFTKSM
use precision
use TYPEDEF
use IIDFTINT
use dft_type
CONTAINS
!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Main kohn-sham matrix driver - calls the nummerical integrater with the 
!> name of a worker routine (II_DFT_KSMGGA,II_DFT_KSMLDA,II_DFT_KSMGGAUNRES,
!> II_DFT_KSMLDAUNRES) which does the work for each grid point
!>
SUBROUTINE II_DFT_KSM(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,ENERGY)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL          :: DOGGA
INTEGER          :: GRDONE,NHTYP,IDMAT
REAL(REALK)      :: SUM,NELE

DOGGA = DFT_ISGGA()
DFTDATA%ENERGY = 0D0
IF(DOGGA) THEN
   IF(NDMAT.EQ.2)THEN !UNRESTRICTED
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_KSMGGAUNRES,DFTDATA)
   ELSE
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_KSMGGA,DFTDATA)
   ENDIF
   DO IDMAT = 1,NDMAT
      DO I = 1, NBAST
         DO J = 1, I-1
            AVERAG = 0.5*(DFTDATA%FKSM(J,I,IDMAT) + DFTDATA%FKSM(I,J,IDMAT))
            DFTDATA%FKSM(J,I,IDMAT) = AVERAG
            DFTDATA%FKSM(I,J,IDMAT) = AVERAG
         END DO
      END DO
   END DO
ELSE 
   IF(NDMAT.EQ.2)THEN !UNRESTRICTED
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_KSMLDAUNRES,DFTDATA)
   ELSE !DEFAULT (RESTRICTED)
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_KSMLDA,DFTDATA)
   END IF
ENDIF
NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
ENERGY = DFTDATA%ENERGY
IF(IPRINT.GE.0) WRITE(LUPRI,'(A,2F20.14,A,E9.2)')&
     &     'KS electrons/energy:', DFTDATA%ELECTRONS, ENERGY,&
     &     ' rel.err:', (DFTDATA%ELECTRONS-NELE)/(NELE)

END SUBROUTINE II_DFT_KSM

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Main molecular gradient driver - calls the nummerical integrater with the 
!> name of a worker routine (II_geoderiv_molgrad_worker) which 
!> does work for each grid point
!>
SUBROUTINE II_dft_geoderiv_molgrad(SETTING,LUPRI,IPRINT,nbast,DMAT,DFTDATA)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!
REAL(REALK)      :: ELE
INTEGER          :: I,J,KMAX,MXPRIM
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL          :: DOGGA
INTEGER          :: GRDONE,NHTYP
REAL(REALK)      :: SUM,NELE
Type(DaltonInput):: DALTON

  CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,1,1,.FALSE.,&
       & II_geoderiv_molgrad_worker,DFTDATA)
  NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
  IF(IPRINT.GE.0) WRITE(LUPRI,'(A,F20.14,A,E9.2)')&
       &     'KS electrons:', DFTDATA%ELECTRONS,' rel.err:', (DFTDATA%ELECTRONS-NELE)/(NELE)

  ! add eventually empirical dispersion correction \Andreas Krapp
  CALL II_DFTDISP(SETTING,DFTDATA%GRAD,3,SETTING%MOLECULE(1)%p%NATOMS,1,LUPRI,IPRINT)

END SUBROUTINE II_DFT_GEODERIV_MOLGRAD

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Worker routine that for a batch of gridpoints build the LDA kohn-sham matrix
!>
SUBROUTINE II_DFT_KSMLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK) :: VXC(NBLEN),VX(5),DFTENE
INTEGER     :: IPNT,I,J
Real(realk), parameter :: D2 = 2.0D0,DUMMY = 0D0
EXTERNAL DFTENE

! LDA Exchange-correlation contribution to Kohn-Sham energy
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000   = \frac{\partial  f}{\partial \rho_{\alpha}}
      CALL dft_funcderiv1(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
      ENERGY = ENERGY + DFTENE(RHO(IPNT,1),DUMMY)*WGHT(IPNT)
      !WARNING the factor 2 is due to F=F_alpha + F_beta = 2 F_alpha 
      VXC(IPNT) = D2*VX(1) 
   ELSE
      VXC(IPNT) = 0.0d0
   ENDIF
END DO

! LDA Exchange-correlation contribution to Kohn-Sham matrix
CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,1,NBLEN,VXC,GAO(:,:,1),DFTDATA%FKSM(:,:,1),DFTHRI)

END SUBROUTINE II_DFT_KSMLDA

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Worker routine that for a batch of gridpoints build the unrestricted LDA kohn-sham matrix
!>
SUBROUTINE II_DFT_KSMLDAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
Real(realk), parameter :: D2 = 2.0D0,DUMMY = 0D0
INTEGER     :: IPNT,I,J,IDMAT
REAL(REALK) :: VXC(NBLEN,NDMAT),VX(5),DFTENEUNRES
EXTERNAL DFTENEUNRES

! LDA Exchange-correlation contribution to Kohn-Sham energy
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR .OR. RHO(IPNT,2) .GT. RHOTHR)THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000   = \frac{\partial  f}{\partial \rho_{\alpha}} 
      !and same for beta
      CALL dft_funcderiv1unres(RHO(IPNT,1),RHO(IPNT,2),DUMMY,DUMMY,WGHT(IPNT),VX)
      ENERGY = ENERGY + DFTENEUNRES(RHO(IPNT,1),RHO(IPNT,2),DUMMY,DUMMY)*WGHT(IPNT)
      VXC(IPNT,1) = VX(1)
      VXC(IPNT,2) = VX(2)
   ELSE
      VXC(IPNT,1) = 0.0d0
      VXC(IPNT,2) = 0.0d0
   ENDIF
END DO
! LDA Exchange-correlation contribution to Kohn-Sham matrix
CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,1,NBLEN,&
     &VXC(:,1),GAO(:,:,1),DFTDATA%FKSM(:,:,1),DFTHRI)
CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,1,NBLEN,&
     &VXC(:,2),GAO(:,:,1),DFTDATA%FKSM(:,:,2),DFTHRI)

END SUBROUTINE II_DFT_KSMLDAUNRES

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> routine that does the actual work see II_dft_ksm.tex
!>
SUBROUTINE II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,IBLSTART,IBLEND,&
     &                    COEF,GAOS,EXCMAT,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in)  :: IBLSTART
!> end gridpoint
INTEGER,intent(in)  :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK) :: TMP(NBLEN,NACTBAST),GAOMAX
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST)
REAL(REALK),allocatable :: GAORED(:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED
REAL(REALK),pointer :: EXCRED(:,:)

ALLOCATE(GAORED(NBLEN,NACTBAST))
NK=IBLEND-IBLSTART+1
NRED = 0 
GAOMAX = 0.0d0
!        Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0d0
      DO K = IBLSTART, IBLEND
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J)))
         GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J)))
      ENDDO
   ENDDO
ENDDO
!        Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J) 
         DO K = IBLSTART, IBLEND
            GAORED(K,NRED) = GAOS(K,J)
         ENDDO
      ENDIF
   ENDDO
ENDDO
IF (NRED.GT.0) THEN
   !  First half-contraction of GAO's with potential
   DO J=1,NRED
      DO K=IBLSTART, IBLEND
         TMP(K,J) =  coef(K)* GAORED(K,J)
      ENDDO
   ENDDO
   !  Second half-contraction of GAO's with potential
   call mem_alloc(EXCRED,NRED,NRED)
   CALL DGEMM('T','N',NRED,NRED,NK,1.d0,&
        &                GAORED,NK,TMP,NK,0.0d0,&
        &                EXCRED,NRED)
   !  Distribute contributions to KS-matrix
!$OMP CRITICAL
   DO JRED=1,NRED         !Jred is reduced index
      J = INXRED(JRED)    !J is orbitalindex
      DO IRED=1,NRED      !Ired is reduced index
         I = INXRED(IRED) !I is orbitalindex
         EXCMAT(I,J) = EXCMAT(I,J) + EXCRED(IRED,JRED)
      ENDDO
   ENDDO
!$OMP END CRITICAL
   call mem_DEALLOC(EXCRED)
ENDIF
DEALLOCATE(GAORED)

END SUBROUTINE II_DFT_DIST_LDA

!> \brief main closed shell GGA kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DFT_KSMGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
Real(realk), parameter :: D2 = 2.0D0,DUMMY = 0D0,D3 = 3.0D0,D05 = 0.5D0
INTEGER     :: IPNT,I,J
REAL(REALK) :: VXC(4,NBLEN),VX(5),DFTENE,GRD,GRDA,A
EXTERNAL DFTENE
!     GGA Exchange-correlation contribution to Kohn-Sham matrix
DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000 = \frac{\partial f}{\partial \rho_{\alpha}}        
      !vx(2) = drvs.df0010 = \frac{\partial f}{\partial |\nabla \rho_{\alpha}|^{2}}
      !vx(3) = drvs.df00001= \frac{\partial f}{\partial \nabla \rho_{\alpha} \nabla \rho_{\beta}}  
      CALL dft_funcderiv1(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
      ENERGY = ENERGY + DFTENE(RHO(IPNT,1),GRD)*WGHT(IPNT)
      VXC(1,IPNT) = D2*VX(1) !this is 2 times greater than the UNRES VXC BECAUSE we need alpha + beta
      IF(GRD.GT.1D-40) THEN
         GRDA = D05*GRD
!         VXC(2,IPNT) = D2*(VX(2)/GRDA + VX(3))  
         A = D2*(VX(2)/GRDA + VX(3))
         VXC(2,IPNT) = A*GRAD(1,IPNT,1)
         VXC(3,IPNT) = A*GRAD(2,IPNT,1)
         VXC(4,IPNT) = A*GRAD(3,IPNT,1)
         !WARNING: this is the same as for the unrestricted case, but since we use \nabla rho 
         !and not \nabla rho_{\alpha} in II_DISTGGA it is the same as having a factor 2 on this 
         !meaning that we take alpha + beta.
      ELSE
         VXC(2,IPNT) = 0D0
         VXC(3,IPNT) = 0D0
         VXC(4,IPNT) = 0D0
      ENDIF
   ELSE
      VXC(:,IPNT) = 0D0
   END IF
END DO
CALL II_DISTGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,NactBast,NBAST,1,NBLEN,&
     &           VXC(:,:),GAO(:,:,1:4),DFTDATA%FKSM(:,:,1),DFTHRI)

END SUBROUTINE II_DFT_KSMGGA

!> \brief main unrestricted GGA kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DFT_KSMGGAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
Real(realk), parameter :: D2 = 2.0D0,DUMMY = 0D0
INTEGER     :: IPNT,I,J
REAL(REALK) :: VXC(NBLEN,2,NDMAT),VXM(NBLEN),VX(5),DFTENEUNRES,GRDA,GRDB
EXTERNAL DFTENEUNRES

!     GGA Exchange-correlation contribution to Kohn-Sham matrix
DO IPNT = 1, NBLEN
   GRDA = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)+&
        &GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   GRDB = SQRT(GRAD(1,IPNT,2)*GRAD(1,IPNT,2)+GRAD(2,IPNT,2)*GRAD(2,IPNT,2)+&
        &GRAD(3,IPNT,2)*GRAD(3,IPNT,2))
   IF((GRDA .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR).OR.&
        &(GRDB .GT. RHOTHR .OR. RHO(IPNT,2).GT.RHOTHR))then 
      !get the functional derivatives 
      !vx(1) = drvs.df1000 = \frac{\partial f}{\partial \rho_{\alpha}}        
      !vx(2) = drvs.df0100 = \frac{\partial f}{\partial \rho_{\beta }}        
      !vx(3) = drvs.df0010 = \frac{\partial f}{\partial |\nabla \rho_{\alpha}|}
      !vx(4) = drvs.df0001 = \frac{\partial f}{\partial |\nabla \rho_{\alpha}|}
      !vx(5) = drvs.df00001= \frac{\partial f}{\partial \nabla \rho_{\alpha} \nabla \rho_{\beta}}  
      IF(GRDA.LT.1D-40) GRDA = 1D-40
      IF(GRDB.LT.1D-40) GRDB = 1D-40
      CALL dft_funcderiv1unres(RHO(IPNT,1),RHO(IPNT,2),GRDA,GRDB,WGHT(IPNT),VX)
      ENERGY = ENERGY + DFTENEUNRES(RHO(IPNT,1),RHO(IPNT,2),GRDA,GRDB)*WGHT(IPNT)
      VXC(IPNT,1,1) = VX(1)
      VXC(IPNT,1,2) = VX(2)
      !WARNING: The factor of 2 in these next coeffficients comes from the fact 
      !that we only build half of these terms and use a symmetrization at the end to get the full contribution
      VXC(IPNT,2,1) = D2*VX(3)/GRDA  
      VXC(IPNT,2,2) = D2*VX(4)/GRDB
      VXM(IPNT) = D2*VX(5) !mixed derivate
   ELSE
      VXC(IPNT,1,1) = 0D0
      VXC(IPNT,2,1) = 0D0
      VXC(IPNT,1,2) = 0D0
      VXC(IPNT,2,2) = 0D0
      VXM(IPNT) = 0D0
   END IF
END DO

! call with drho_alpha dgradrho_alpha dmixed and gradA, gradB
CALL II_DISTGGABUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,NactBast,NBAST,1,NBLEN,&
     &  VXC(:,1,1),VXC(:,2,1),VXM,GAO(:,:,1:4),GRAD(:,:,1),GRAD(:,:,2),DFTDATA%FKSM(:,:,1),DFTHRI)
! call with drho_beta dmixed dgradrho_beta and gradA, gradB
CALL II_DISTGGABUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,NactBast,NBAST,1,NBLEN,&
     &  VXC(:,1,2),VXM,VXC(:,2,2),GAO(:,:,1:4),GRAD(:,:,1),GRAD(:,:,2),DFTDATA%FKSM(:,:,2),DFTHRI)

END SUBROUTINE II_DFT_KSMGGAUNRES

!> \brief a distribution routine 
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DISTGGABUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,IBLSTART,IBLEND,&
     &                    COEF1,COEFA,COEFB,GAOS,GRADA,GRADB,EXCMAT,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in)  :: IBLSTART
!> end gridpoint
INTEGER,intent(in)  :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEFA(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEFB(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,4)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK) :: TMP(NBLEN,NACTBAST),GRADA(:,:),GRADB(:,:),GAOMAX
REAL(REALK) :: GAORED(NBLEN,NACTBAST,4),GAOGMX(NACTBAST)
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,INXRED(NACTBAST)
INTEGER     :: IRED,JRED,NRED
REAL(REALK),pointer :: EXCRED(:,:)

NK=IBLEND-IBLSTART+1
NRED = 0
GAOMAX = 0.0d0
! Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0d0
      DO K = IBLSTART, IBLEND
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,1)))
         DO I=1,4
            GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,I)))
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO K = IBLSTART, IBLEND
            GAORED(K,NRED,1) = GAOS(K,J,1)
            GAORED(K,NRED,2) = GAOS(K,J,2)
            GAORED(K,NRED,3) = GAOS(K,J,3)
            GAORED(K,NRED,4) = GAOS(K,J,4)
         ENDDO
      ENDIF
   ENDDO
ENDDO
IF (NRED.GT.0) THEN
!  First half-contraction of GAO's with potential
   DO J=1,NRED !J is reduced index
      DO K=IBLSTART, IBLEND
         TMP(K,J) =  coef1(K)* GAORED(K,J,1)&
              &    + coefA(K)*(GAORED(K,J,2)*GRADA(1,K)+&
              &                GAORED(K,J,3)*GRADA(2,K)+&
              &                GAORED(K,J,4)*GRADA(3,K))&
              &    + coefB(K)*(GAORED(K,J,2)*GRADB(1,K)+&
              &                GAORED(K,J,3)*GRADB(2,K)+&
              &                GAORED(K,J,4)*GRADB(3,K))
      ENDDO
   ENDDO
!  Second half-contraction of GAO's with potential
   CALL MEM_ALLOC(EXCRED,NRED,NRED)
   CALL DGEMM('T','N',NRED,NRED,NK,1.0d0,&
        &                GAORED,NK,TMP,NK,0.0d0,&
        &                EXCRED,NRED)
!  Distribute contributions to KS-matrix
!$OMP CRITICAL
   DO JRED=1,NRED         !Jred is reduced index
      J = INXRED(JRED)    !J is orbital index
      DO IRED=1,NRED      !Ired is reduced index
         I = INXRED(IRED) !I is orbital index
         EXCMAT(I,J) = EXCMAT(I,J) + EXCRED(IRED,JRED)
      ENDDO
   ENDDO
!$OMP END CRITICAL
   CALL MEM_DEALLOC(EXCRED)
ENDIF

END SUBROUTINE II_DISTGGABUNRES

!> \brief the main GGA distribution routine 
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DISTGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,IBLSTART,IBLEND,&
     &                    COEF,GAOS,EXCMAT,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in)  :: IBLSTART
!> end gridpoint
INTEGER,intent(in)  :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(4,NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,4)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK) :: TMP(NBLEN,NACTBAST),GAOMAX
REAL(REALK) :: GAORED(NBLEN,NACTBAST,4),GAOGMX(NACTBAST)
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,INXRED(NACTBAST)
INTEGER     :: IRED,JRED,NRED
REAL(REALK),pointer :: EXCRED(:,:)

NK=IBLEND-IBLSTART+1
NRED = 0
GAOMAX = 0.0d0
! Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0d0
      DO K = IBLSTART, IBLEND
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,1)))
         GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,1)))
      ENDDO
      DO I=2,4
         DO K = IBLSTART, IBLEND
            GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,I)))
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO I=1,4
            DO K = IBLSTART, IBLEND
               GAORED(K,NRED,I) = GAOS(K,J,I)
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDDO
IF (NRED.GT.0) THEN
!  First half-contraction of GAO's with potential
   DO J=1,NRED
      DO K=IBLSTART, IBLEND
         TMP(K,J) =  coef(1,K)*GAORED(K,J,1)&
              &   + coef(2,K)*GAORED(K,J,2)&
              &   + coef(3,K)*GAORED(K,J,3)&
              &   + coef(4,K)*GAORED(K,J,4)
      ENDDO
   ENDDO
!  Second half-contraction of GAO's with potential
   CALL MEM_ALLOC(EXCRED,NRED,NRED)
   CALL DGEMM('T','N',NRED,NRED,NK,1.0d0,&
        &                GAORED,NK,TMP,NK,0.0d0,&
        &                EXCRED,NRED)
!  Distribute contributions to KS-matrix
!$OMP CRITICAL
   DO JRED=1,NRED          !Jred is reduced index
      J = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED       !Ired is reduced index
         I = INXRED(IRED)  !I is orbital index
         EXCMAT(I,J) = EXCMAT(I,J) + EXCRED(IRED,JRED)
      ENDDO
   ENDDO
!$OMP END CRITICAL
   CALL MEM_DEALLOC(EXCRED)
ENDIF

END SUBROUTINE II_DISTGGA

!> \brief GGA Exchange-correlation contribution to molecular gradient
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_geoderiv_molgrad_worker(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAOS,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
  IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
Real(realk), parameter :: D2 = 2.0D0,DUMMY = 0D0,D05 = 0.5D0
REAL(REALK) :: VXC(2,NBLEN),VX(5),GRDA,GRD,GAOMAX,DMAX
INTEGER  :: IPNT,I,J,JBL,IBL,K
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL :: DOGGA
INTEGER     :: KVALS(3,3)
REAL(REALK) :: GAOGMX(NACTBAST)
REAL(REALK),allocatable :: GAORED(:,:,:),GDRED(:,:,:)
REAL(REALK),allocatable :: DRED(:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,NRED,orb2atom(nbast)
INTEGER     :: atom(NACTBAST),iatom,IX,K1,K2,K3,KA,ik,jk
REAL(REALK) :: FRC,GA,GA2,GFS

orb2atom = DFTDATA%orb2atom
KVALS(1:3,1) = (/1, 2, 3/)
KVALS(1:3,2) = (/2, 4, 5/)
KVALS(1:3,3) = (/3, 5, 6/)
DOGGA = DFT_ISGGA()
IF (DOGGA) THEN
   DO IPNT = 1, NBLEN
      GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
           &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
      IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
         !get the functional derivatives 
         CALL dft_funcderiv1(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
         VXC(1,IPNT) = D2*VX(1) 
         IF(GRD.GT.1D-40) THEN
            GRDA = D05*GRD
            VXC(2,IPNT) = (VX(2)/GRDA + VX(3))  
         ELSE
            VXC(2,IPNT) = 0D0
         ENDIF
      ELSE
         VXC(1,IPNT) = 0D0
         VXC(2,IPNT) = 0D0
      END IF
   END DO
ELSE
   DO IPNT = 1, NBLEN
      IF(RHO(IPNT,1).GT.RHOTHR) THEN
         !get the functional derivatives 
         !vx(1) = drvs.df1000 = \frac{\partial f}{\partial \rho_{\alpha}}        
         CALL dft_funcderiv1(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
         VXC(1,IPNT) = D2*VX(1) 
      ELSE
         VXC(1,IPNT) = 0D0
      END IF
   END DO
ENDIF

GAOMAX = 0.0d0
! Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0d0
      DO I=1,NTYPSO
         DO K = 1,NBLEN
            GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,I)))
            GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,I)))
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Set up maximum density-matrix elements
DMAX = 0.0d0
DO JBL=1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)        !J is active index
      DO IBL=1, NBLOCKS
         DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)  !I is active index
            DMAX = MAX(DMAX,ABS(DMAT(I,J,1)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
! Count reduced number of AO's
NRED = 0
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
         NRED = NRED + 1
      ENDIF
   ENDDO
ENDDO

IF (NRED.GT.0) THEN
   ALLOCATE(DRED(NRED,NRED))
   ALLOCATE(GAORED(NBLEN,NRED,NTYPSO))
   ALLOCATE(GDRED(NBLEN,NRED,NTYPSO))
   ! Set up reduced Gaussian AO's
   IRED = 0
   DO IBL=1, NBLOCKS
      DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
         IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
            IRED = IRED + 1
            INXRED(IRED) = I
            DO J=1,NTYPSO
               DO K = 1, NBLEN
                  GAORED(K,IRED,J)  = GAOS(K,I,J)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ! Set up reduced density-matrix
   DO JRED=1,NRED            !Jred is reduced index
      J = INXRED(JRED)       !J is active index
      DO IRED=1,NRED         !Ired is reduced index
         I = INXRED(IRED)    !I is active index
         DRED(IRED,JRED) = DMAT(I,J,1)
      ENDDO
   ENDDO
   ! Set up reduced coordinate index in gradient
   DO IRED=1,NRED
      I = INXACT(INXRED(IRED)) !I is orbital index
      atom(IRED) = orb2atom(I)
   ENDDO
   
   IF (DOGGA) THEN
      ! Density-matrix contraction  
      ! \chi_{\mu}_{A,B} D_{\mu \nu}
      ! \frac{\partial \chi_{\mu}}{\frac \partial x} D_{\mu \nu}
      ! \frac{\partial \chi_{\mu}}{\frac \partial y} D_{\mu \nu}
      ! \frac{\partial \chi_{\mu}}{\frac \partial z} D_{\mu \nu}
      DO J=1,4
         CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0d0,GAORED(1,1,J),&
         &                  NBLEN,DRED,NRED,0.0d0,GDRED(1,1,J),NBLEN )
      ENDDO
      DO IRED=1,NRED
         iatom = atom(IRED)
         KA = INXRED(IRED) !KA is active index
         DO IX=1,3
            K1 = KVALS(1,IX) + 4
            K2 = KVALS(2,IX) + 4
            K3 = KVALS(3,IX) + 4
            FRC = 0.D0
            ! Assuming E_{\xc}=\int f[\rho ,\nabla \rho] d\textbf{r}
            ! \frac{ \partial E^{xc}[\rho]}{\partial R} =  
            ! \int \frac{\partial f }{\partial \rho} \frac{\partial \rho(\textbf{r})}{\partial R} d\textbf{r} + 
            ! \int \frac{\partial f }{\partial |\nabla \rho_{\alpha}|} \frac{\nabla \rho(\textbf{r})}{|\nabla \rho_{\alpha}|}  
            ! \frac{\partial \nabla \rho(\textbf{r})}{\partial R} d\textbf{r}
            ! VXC(1,I) = \frac{\partial f }{\partial \rho}
            ! VXC(2,I) = \int \frac{\partial f }{\partial |\nabla \rho_{\alpha}|}
            DO I = 1, NBLEN
               !\frac{\partial \chi_{\mu}}{\partial R_{\gamma}}
               GA  = GAOS(I,KA,IX+1)
               !\nabla \rho \nabla \chi_{\mu} D_{\mu \nu}
               GFS = GRAD(1,I,1)*GDRED(I,IRED,2)+GRAD(2,I,1)*GDRED(I,IRED,3)+GRAD(3,I,1)*GDRED(I,IRED,4)
               !\nabla \rho \frac{\partial \nabla \chi_{\mu}}{\partial R_{\gamma}}
               GA2 = GRAD(1,I,1)*GAOS(I,KA,K1)+GRAD(2,I,1)*GAOS(I,KA,K2)+GRAD(3,I,1)*GAOS(I,KA,K3)
               !\int VXC(1)\frac{\partial \chi_{\mu}}{\partial R_{\gamma}}\chi_{\nu}D_{\mu \nu}
               ! + VXC(2) \nabla \rho \frac{\partial \nabla \chi_{\mu}}{\partial R_{\gamma}} \chi_{\nu}D_{\mu \nu}
               ! + VXC(2) \nabla \rho \nabla \chi_{\mu} \frac{\partial \chi_{\mu}}{\partial R_{\gamma}} D_{\mu \nu}
               FRC = FRC + VXC(1,I)*GDRED(I,IRED,1)*GA + VXC(2,I)*(GDRED(I,IRED,1)*GA2 + GFS*GA)
            END DO
            !$OMP CRITICAL
            DFTDATA%GRAD(IX,iatom) = DFTDATA%GRAD(IX,iatom) - FRC
            !$OMP END CRITICAL
         ENDDO ! IX
      ENDDO ! IRED
   ELSE
   ! Density-matrix contraction
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0d0,GAORED,&
      &                 NBLEN,DRED,NRED,0.0d0,GDRED,NBLEN    )
      DO IRED=1,NRED
         iatom = atom(IRED)
         KA = INXRED(IRED)  !KA is active index
         DO IX=1,3
            FRC = 0.D0
            DO I = 1, NBLEN
               FRC = FRC + VXC(1,I)*GDRED(I,IRED,1)*GAOS(I,KA,IX+1)
            END DO
            !$OMP CRITICAL
            DFTDATA%GRAD(IX,iatom) = DFTDATA%GRAD(IX,iatom) - FRC
            !$OMP END CRITICAL
         ENDDO ! IX
      ENDDO ! IA
   ENDIF ! IF DOGGA
   DEALLOCATE(DRED)
   DEALLOCATE(GAORED)
   DEALLOCATE(GDRED)
ENDIF !NRED GT 0

END SUBROUTINE II_GEODERIV_MOLGRAD_WORKER

!> \brief Main linear response driver
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE ii_dft_linrsp(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL          :: DOGGA
INTEGER          :: GRDONE,NHTYP,IDMAT,IBMAT

  DOGGA = DFT_ISGGA()
  DFTDATA%ENERGY = 0D0
  IF(DOGGA) THEN
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
             & II_DFT_LINRSPGGAUNRES,DFTDATA)
     ELSE
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
             & II_DFT_LINRSPGGA,DFTDATA)
     ENDIF
     DO IBMAT = 1,DFTDATA%NBMAT
        DO I = 1, NBAST
           DO J = 1, I-1
              AVERAG = 0.5*(DFTDATA%FKSM(J,I,IBMAT) + DFTDATA%FKSM(I,J,IBMAT))
              DFTDATA%FKSM(J,I,IBMAT) = AVERAG
              DFTDATA%FKSM(I,J,IBMAT) = AVERAG
           END DO
        END DO
     END DO
  ELSE 
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
             & II_DFT_LINRSPLDAUNRES,DFTDATA)
     ELSE !DEFAULT (RESTRICTED)
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
             & II_DFT_LINRSPLDA,DFTDATA)
     END IF
  ENDIF
  NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
  IF(IPRINT.GE.0) WRITE(LUPRI,'(A,F20.14,A,E9.2)')&
       &     'KS electrons:', DFTDATA%ELECTRONS,&
       &     ' rel.err:', (DFTDATA%ELECTRONS-NELE)/(NELE)

END SUBROUTINE II_DFT_LINRSP

!> \brief Main LDA linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ii_dft_linrsplda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK),allocatable :: EXPVAL(:,:),VXC(:,:)
!EXPVAL(NBLEN,DFTDATA%NBMAT),VXC(NBLEN,DFTDATA%NBMAT),VX(9),DFTENE
REAL(REALK) :: VX(14),DFTENE
INTEGER     :: I,J,NBMAT,IPNT,IBMAT,nred
LOGICAL     :: DOCALC
Real(realk), parameter :: D2 = 2.0D0,DUMMY = 0D0
REAL(REALK) :: fRR

allocate(EXPVAL(NBLEN,DFTDATA%NBMAT))
allocate(VXC(NBLEN,DFTDATA%NBMAT))
NBMAT = DFTDATA%NBMAT
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
   !get expectation value of BMAT = \sum_{\mu \nu} \chi_{\mu} \chi_{\nu} BMAT_{\mu \nu}
 call II_get_expectationval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
        &Nactbast,NBAST,GAO,EXPVAL,DFTDATA%BMAT,DFTDATA%nBMAT,DFTHRI,NRED)
 IF(NRED.GT.0)THEN
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
    !get the functional derivatives 
    CALL dft_funcderiv2(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
    !fRR = 0.5*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
    fRR = VX(6) + VX(7)     
    DO IBMAT = 1,NBMAT
     !THE factor 0.5 IN fRR cancels a factor 2.0 in VXC so  
     VXC(IPNT,IBMAT) = D2*fRR*EXPVAL(IPNT,IBMAT) ! D4*fRR*EXPVAL(IPNT,IBMAT)
     !WARNING the factor 4 is due to 
     ! G_{\rho \sigma} &=& \frac{\delta (F^{\alpha} + F^{\beta})}{\delta D^{\alpha}_{\nu \mu}}\kappa^{\alpha}_{\nu \mu}\\
     ! &+& \frac{\delta (F^{\alpha} + F^{\beta})}{\delta D^{\beta}_{\nu \mu}}\kappa^{\beta}_{\nu \mu}\\
     ! G_{\rho \sigma} &=& 4 \frac{\delta (F^{\alpha}{\delta D^{\alpha}_{\nu \mu}} \kappa^{\alpha}_{\nu \mu}\\
     ! G_{\rho \sigma} &=& 4 \int \frac{\partial^{2} f }{\partial \rho^{2}_{\alpha}} \Omega_{\rho \sigma}
     ! \Omega_{\mu \nu} \kappa^{\alpha}_{\mu \nu} d\textbf{r}
    ENDDO
   ELSE
    VXC(IPNT,:) = 0.0d0
   ENDIF
  END DO
  DO IBMAT = 1,NBMAT
   CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
        & 1,NBLEN,VXC(:,IBMAT),GAO(:,:,1),DFTDATA%FKSM(:,:,IBMAT),DFTHRI)
  ENDDO
 ENDIF
ENDIF
deallocate(EXPVAL)
deallocate(VXC)

END SUBROUTINE II_DFT_LINRSPLDA

!> \brief computes the expectation value of a matrix BMAT
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_expectationval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,GAOS,EXPVAL,BMAT,nBmat,DFTHRI,NRED)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,1)
!> The expectation value of the B matrix
REAL(REALK),intent(inout) :: EXPVAL(NBLEN,NBMAT)
!> The B matrix (some perturbed density matrix) 
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in) :: NBMAT
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of reduced orbitals
INTEGER,intent(inout) :: NRED
!
REAL(REALK) :: TMP(NBLEN,NACTBAST),GAOMAX,BMAX
REAL(REALK) :: GAORED(NBLEN,NACTBAST),GAOGMX(NACTBAST)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,INXRED(NACTBAST)
INTEGER     :: IRED,JRED,IORB,JORB,IBMAT
REAL(REALK),allocatable :: BRED(:,:,:)

!IF(NBMAT .GT.1)CALL LSQUIT('II_get_expectationval_lda does not work correctly for more than 1 nbmat',lupri)
NRED = 0
GAOMAX = 0.0d0
! Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0d0
      DO K = 1, NBLEN
         GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,1)))
      ENDDO
      GAOMAX = MAX(GAOMAX,GAOGMX(J))
   ENDDO
ENDDO

! Set up maximum response-vector elements
BMAX = 0.0d0
DO IBMAT=1,NBMAT
DO JBL=1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL) !J is active index
      JORB = INXACT(J)                 !JORB is orbital index 
      DO IBL=1, NBLOCKS
         DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)  !I is active index
            IORB = INXACT(I)                  !IORB is orbital index 
            BMAX = MAX(BMAX,ABS(BMAT(IORB,JORB,IBMAT)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
ENDDO

! Set up reduced Gaussian AO's
NRED=0
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL) !J is active index
      IF (GAOGMX(J)*GAOMAX*BMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO K = 1, NBLEN 
            GAORED(K,NRED) = GAOS(K,J,1)
         ENDDO
      ENDIF
   ENDDO
ENDDO

IF (NRED.GT.0) THEN
   ! Set up reduced density-matrix
   ALLOCATE(BRED(NRED,NRED,NBMAT))
   DO IBMAT=1,NBMAT
    DO JRED=1,NRED         !Jred is reduced index 
      J = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED       !Ired is reduced index 
         I = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED,IBMAT) = BMAT(I,J,IBMAT)
      ENDDO
    ENDDO
   ENDDO
   
   ! First half-contraction of Gaussian AO with density-matrix
   DO IBMAT=1,NBMAT
     CALL DGEMM("N","N",NBLEN,NRED,NRED,1D0,GAORED,NBLEN,&
        &     BRED(:,:,IBMAT),NRED,0.0d0,TMP,NBLEN)
     ! Second half-contraction 
     DO K = 1, NBLEN
        EXPVAL(K,IBMAT) = GAORED(K,1)*TMP(K,1)
     END DO
     DO I = 2, NRED  !I is reduced index
        DO K = 1, NBLEN
           EXPVAL(K,IBMAT) = EXPVAL(K,IBMAT) + GAORED(K,I)*TMP(K,I)
        END DO
     END DO
   ENDDO
   DEALLOCATE(BRED)
ENDIF

END SUBROUTINE II_GET_EXPECTATIONVAL_LDA

!> \brief main unrestricted LDA linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_LINRSPLDAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
CALL LSQUIT('not implemented yet',lupri)

END SUBROUTINE II_DFT_LINRSPLDAUNRES

!> \brief main GGA linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_LINRSPGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK),allocatable :: EXPVAL(:,:),VXC(:,:,:),EXPGRAD(:,:,:)
REAL(REALK) :: VX(14),DFTENE,GRD,GRDA,MIXEDGRDA
INTEGER     :: I,J,NBMAT,IPNT,IBMAT,nred
LOGICAL     :: DOCALC
Real(realk), parameter :: D4 = 4.0D0, DUMMY = 0D0,D05 = 0.5D0
Real(realk), parameter :: D2 = 2.0D0, D8 = 8.0D0,D025 = 0.25D0
REAL(REALK) :: fR,fZ,fRR,fRZ,fZZ,fRG,fZG,fGG,fG,A,B
NBMAT = DFTDATA%NBMAT
allocate(EXPVAL(NBLEN,NBMAT))
allocate(EXPGRAD(3,NBLEN,NBMAT))
allocate(VXC(4,NBLEN,NBMAT))
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
   !CALC
   !EXPVAL : the expectation value of BMAT = \sum_{\mu \nu} \chi_{\mu} \chi_{\nu} BMAT_{\mu \nu}
   !EXPGRAD: the gradient components of BMAT = \sum_{\mu \nu} \nabla (\chi_{\mu} \chi_{\nu}) BMAT_{\mu \nu}
 call II_get_expectationval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
        &Nactbast,NBAST,GAO,EXPVAL,EXPGRAD,DFTDATA%BMAT,DFTDATA%nBMAT,DFTHRI,NRED)
 IF(NRED.GT.0)THEN
  DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD.LT.1D-40) GRD = 1D-40
   GRDA = D05*GRD
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
    !get the functional derivatives 
    CALL dft_funcderiv2(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
    fR  = D05*(VX(1) + VX(2))   !0.5*(drvs.df1000 + drvs.df0100);
    fZ  = VX(3)                    !drvs.df0010;
    fRR = D05*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
    fRZ = D05*(VX(8) + VX(9))   !0.5*(drvs.df1010 + drvs.df1001);
    fZZ = D05*(VX(11) + VX(12)) !0.5*(drvs.df0020 + drvs.df0011);
    fRG = D05*VX(10)             !0.5*drvs.df10001;   
    fZG = D05*VX(13)             !0.5*drvs.df00101; 
    fGG = D025*VX(14)            !0.25*drvs.df00002; 
    fG  = D05*VX(5)               !0.5*drvs.df00001;  
    DO IBMAT = 1,NBMAT
     MIXEDGRDA = (EXPGRAD(1,IPNT,IBMAT)*GRAD(1,IPNT,1)&
          &+EXPGRAD(2,IPNT,IBMAT)*GRAD(2,IPNT,1)&
          &+EXPGRAD(3,IPNT,IBMAT)*GRAD(3,IPNT,1))
     !the LDA part
     VXC(1,IPNT,IBMAT) =D4*fRR*EXPVAL(IPNT,IBMAT)+D4*(fRZ/GRD+fRG)*MIXEDGRDA
     !the non LDA parts
     A = D8*((fRZ/GRD + fRG)*EXPVAL(IPNT,IBMAT)&
          & + (((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)*MIXEDGRDA)
     B= D8*(fZ/GRD + fG)
     VXC(2,IPNT,IBMAT) = A*GRAD(1,IPNT,1)+B*EXPGRAD(1,IPNT,IBMAT)
     VXC(3,IPNT,IBMAT) = A*GRAD(2,IPNT,1)+B*EXPGRAD(2,IPNT,IBMAT)
     VXC(4,IPNT,IBMAT) = A*GRAD(3,IPNT,1)+B*EXPGRAD(3,IPNT,IBMAT)
    ENDDO
   ELSE
    VXC(:,IPNT,:) = 0.0d0
   ENDIF
  END DO
  DO IBMAT = 1,NBMAT
     CALL II_DISTGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,1,NBLEN,&
     &                    VXC(:,:,IBMAT),GAO,DFTDATA%FKSM(:,:,IBMAT),DFTHRI)
  ENDDO
 ENDIF
ENDIF
deallocate(EXPVAL)
deallocate(EXPGRAD)
deallocate(VXC)

END SUBROUTINE II_DFT_LINRSPGGA

!> \brief computes the expectation value and and gradient of it
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_expectationval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,GAOS,EXPVAL,EXPGRAD,BMAT,nBMAT,DFTHRI,NRED)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,4)
!> The expectation value of the B MATRIX
REAL(REALK),intent(inout) :: EXPVAL(NBLEN,NBMAT)
!> The gradient of the expectation value of the B MATRIX
REAL(REALK),intent(inout) :: EXPGRAD(3,NBLEN,NBMAT)
!> The B matrix (some perturbed density matrix) 
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in) :: NBMAT
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of reduced orbitals
INTEGER,intent(inout) :: NRED
!
REAL(REALK) :: TMP(NBLEN,NACTBAST),GAOMAX,BREDIJ,BREDJI
REAL(REALK) :: GAORED2(NBLEN,NACTBAST,4),GAOGMX(NACTBAST),BMAX
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,INXRED(NACTBAST)
INTEGER     :: IRED,JRED,IORB,JORB,NK,IBMAT
REAL(REALK),parameter :: D05 = 0.5D0
REAL(REALK),allocatable :: BRED(:,:)

GAOMAX = 0.0d0
! Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0d0
      DO K = 1,NBLEN
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,1)))
         DO I=1,4
            GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,I)))
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Set up maximum response-vector elements
BMAX = 0.0d0
DO IBMAT=1,NBMAT
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)  !I is active index
      IORB = INXACT(I)                  !IORB is orbital index
      DO JBL=1, NBLOCKS
         DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)  !J is active index
            JORB = INXACT(J)                  !JORB is orbital index
            BMAX = MAX(BMAX,ABS(BMAT(IORB,JORB,IBMAT)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
ENDDO
! Set up reduced Gaussian AO's
NRED=0
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX*BMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO K = 1, NBLEN
            GAORED2(K,NRED,1) = GAOS(K,J,1)
            GAORED2(K,NRED,2) = GAOS(K,J,2)
            GAORED2(K,NRED,3) = GAOS(K,J,3)
            GAORED2(K,NRED,4) = GAOS(K,J,4)
         ENDDO
      ENDIF
   ENDDO
ENDDO
   
IF (NRED.GT.0) THEN
   ! Set up reduced density-matrix
   ALLOCATE(BRED(NRED,NRED))
   DO IBMAT=1,NBMAT
    DO JRED=1,NRED        !JRED is reduced index
       J = INXRED(JRED)   !J is orbital index
       DO IRED=1,NRED     !IRED is reduced index
          I = INXRED(IRED)!I is orbital index
          BRED(IRED,JRED) = BMAT(I,J,IBMAT)
       ENDDO
    ENDDO
   
    ! First half-contraction of Gaussian AO with density-matrix
    DO J = 1, NRED
       BREDIJ = BRED(1,J)
       BREDJI = BRED(J,1)
       DO K = 1, NBLEN
          TMP(K,J) = GAORED2(K,1,1)*BREDIJ + GAORED2(K,1,1)*BREDJI
       END DO
       DO I = 2, NRED
          BREDIJ = BRED(I,J)
          BREDJI = BRED(J,I)
          DO K = 1, NBLEN
             TMP(K,J) = TMP(K,J) + GAORED2(K,I,1)*BREDIJ + GAORED2(K,I,1)*BREDJI
          END DO
       ENDDO
    ENDDO
      ! Second half-contraction
    DO K = 1, NBLEN
      EXPVAL(K,IBMAT) = D05*GAORED2(K,1,1)*TMP(K,1)
    END DO
    DO I = 2, NRED
      DO K = 1, NBLEN
         EXPVAL(K,IBMAT) = EXPVAL(K,IBMAT) + D05*GAORED2(K,I,1)*TMP(K,I)
      END DO
    END DO
    DO K = 1, NBLEN
      EXPGRAD(1,K,IBMAT) = GAORED2(K,1,2)*TMP(K,1)
      EXPGRAD(2,K,IBMAT) = GAORED2(K,1,3)*TMP(K,1)
      EXPGRAD(3,K,IBMAT) = GAORED2(K,1,4)*TMP(K,1)
    ENDDO
    DO I = 2, NRED
      DO K = 1, NBLEN
         EXPGRAD(1,K,IBMAT) = EXPGRAD(1,K,IBMAT) + GAORED2(K,I,2)*TMP(K,I)
         EXPGRAD(2,K,IBMAT) = EXPGRAD(2,K,IBMAT) + GAORED2(K,I,3)*TMP(K,I)
         EXPGRAD(3,K,IBMAT) = EXPGRAD(3,K,IBMAT) + GAORED2(K,I,4)*TMP(K,I)
      ENDDO
    ENDDO
   ENDDO
   DEALLOCATE(BRED)
ENDIF
END SUBROUTINE II_GET_EXPECTATIONVAL_GGA

!> \brief main unrestricted GGA linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_LINRSPGGAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
CALL LSQUIT('not implemented yet',lupri)

END SUBROUTINE II_DFT_LINRSPGGAUNRES

!> \brief Main quadratic response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ii_dft_quadrsp(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL          :: DOGGA
INTEGER          :: GRDONE,NHTYP,IDMAT,IBMAT

  DOGGA = DFT_ISGGA()
  DFTDATA%ENERGY = 0D0
  IF(DOGGA) THEN
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
!             & II_DFT_QUADRSPGGAUNRES,DFTDATA)
        CALL LSQUIT('II_DFT_QUADRSPGGAUNRES NOT DONE YET',lupri)
     ELSE
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
             & II_DFT_QUADRSPGGA,DFTDATA)
!        CALL LSQUIT('II_DFT_QUADRSPGGA NOT DONE YET',lupri)
     ENDIF
     DO IDMAT = 1,DFTDATA%NDMAT
        DO I = 1, NBAST
           DO J = 1, I-1
              AVERAG = 0.5*(DFTDATA%FKSM(J,I,IDMAT) + DFTDATA%FKSM(I,J,IDMAT))
              DFTDATA%FKSM(J,I,IDMAT) = AVERAG
              DFTDATA%FKSM(I,J,IDMAT) = AVERAG
           END DO
        END DO
     END DO
  ELSE 
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
!             & II_DFT_QUADRSPLDAUNRES,DFTDATA)
        CALL LSQUIT('II_DFT_QUADRSPLDAUNRES NOT DONE YET',lupri)
     ELSE !DEFAULT (RESTRICTED)
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
             & II_DFT_QUADRSPLDA,DFTDATA)
     END IF
  ENDIF
  NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
  IF(IPRINT.GE.0) WRITE(LUPRI,'(A,F20.14,A,E9.2)')&
       &     'KS electrons:', DFTDATA%ELECTRONS,&
       &     ' rel.err:', (DFTDATA%ELECTRONS-NELE)/(NELE)

END SUBROUTINE II_DFT_QUADRSP

!> \brief LDA Exchange-correlation contribution to quadratic response 
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ii_dft_quadrsplda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK),allocatable :: EXPVAL(:,:),VXC(:)
!EXPVAL(NBLEN,DFTDATA%NBMAT),VXC(NBLEN,DFTDATA%NBMAT),VX(9),DFTENE
REAL(REALK) :: VX(27),DFTENE
INTEGER     :: I,J,NBMAT,IPNT,IBMAT,nred
LOGICAL     :: DOCALC
Real(realk), parameter :: D2 = 2.0D0,DUMMY = 0D0,D3 = 3.0D0
REAL(REALK) :: fRRR

allocate(EXPVAL(NBLEN,DFTDATA%NBMAT))
allocate(VXC(NBLEN))
NBMAT = DFTDATA%NBMAT
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
!get expectation value of BMAT=\sum_{\mu \nu}\chi_{\mu}\chi_{\nu}BMAT_{\mu \nu}
 call II_get_expectationval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
        &Nactbast,NBAST,GAO,EXPVAL,DFTDATA%BMAT,DFTDATA%nBMAT,DFTHRI,NRED)
 IF(NRED.GT.0)THEN
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
    !get the functional derivatives 
    CALL dft_funcderiv3(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
!    fR  = VX(1)              !drvs.df1000;
!radovan: factor 0.5 "missing" here compared to "true" derivative wrt RR i'm 
!sure it's compensated elsewhere but can be confusing when comparing with 
!other codes
!     fRR = (VX(6) + VX(7))   !(drvs.df2000 + drvs.df1100);
!radovan: factor 0.25 "missing" here compared to "true" derivative wrt RRR
!        (i'm sure it's compensated elsewhere)
!         but can be confusing when comparing with other codes
    fRRR = (VX(15)+D3*VX(16))     !(drvs.df3000 + 3*drvs.df2100);
    print*,'fRRR',fRRR,'EXPVAL(IPNT,1)*EXPVAL(IPNT,2)',EXPVAL(IPNT,1)*EXPVAL(IPNT,2)  
!    VXCB(IPNT) = fRR*EXPVALB(IPNT) 
!    VXCC(IPNT) = fRR*EXPVALC(IPNT) 
    VXC(IPNT) = D2*fRRR*EXPVAL(IPNT,1)*EXPVAL(IPNT,2)  
   ELSE
    VXC(IPNT) = 0.0d0
   ENDIF
  END DO
   CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
        & 1,NBLEN,VXC(:),GAO(:,:,1),DFTDATA%FKSM(:,:,1),DFTHRI)
 ENDIF
ENDIF
deallocate(EXPVAL)
deallocate(VXC)

END SUBROUTINE II_DFT_QUADRSPLDA

!> \brief GGA Exchange-correlation contribution to quadratic response 
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_QUADRSPGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,MXBLLEN,COORD,WGHT,ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK),allocatable :: EXPVAL(:,:),VXC(:,:),EXPGRAD(:,:,:)
REAL(REALK) :: VX(27),DFTENE,GRD,GRDA,MIXEDGRDA
INTEGER     :: I,J,NBMAT,IPNT,IBMAT,nred
LOGICAL     :: DOCALC
Real(realk), parameter :: D4 = 4.0D0, DUMMY = 0D0,D05 = 0.5D0,D3 = 3.d0
Real(realk), parameter :: D2 = 2.0D0, D8 = 8.0D0,D025 = 0.25D0
REAL(REALK) :: GRD2,GRDA2,GRDA3
REAL(REALK) :: fRZ,fRG,fZZ,fRRR,fRRZ,fRRG,fRRGX,fRZZ,fZZZ,gradY,gradZ,gradYZ,A,B,C
NBMAT = DFTDATA%NBMAT
allocate(EXPVAL(NBLEN,NBMAT))
allocate(EXPGRAD(3,NBLEN,NBMAT))
allocate(VXC(4,NBLEN))
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
 call II_get_expectationval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
        &Nactbast,NBAST,GAO,EXPVAL,EXPGRAD,DFTDATA%BMAT,DFTDATA%nBMAT,DFTHRI,NRED)
 IF(NRED.GT.0)THEN
  DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD.LT.1D-40) GRD = 1D-40
   GRDA = D05*GRD
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
    !get the functional derivatives 
    GRD2 = GRD*GRD
    GRDA2 = GRDA*GRDA
    GRDA3 = GRDA2*GRDA
    CALL dft_funcderiv3(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
    fRZ = (VX(8)+VX(9))/GRD              !(drvs.df1010 + drvs.df1001)/(2*grada)
    fRG = D2*VX(10)                      !2*drvs.df10001   
    fZZ = (VX(11)+VX(12))/GRD2-VX(3)/(GRD2*GRDA) !(drvs.df0020 + drvs.df0011)/(4*grada2)-drvs.df0010/(4*grada3)
    fRRR = VX(15)+D3*VX(16)              !(drvs.df3000 + 3*drvs.df2100) 
    fRRZ = (VX(17)+VX(18)+D2*VX(20))/GRD !(drvs.df2010+drvs.df2001+2*drvs.df1110)/(2*grada)
    fRRG = VX(19)+VX(21)                 !drvs.df20001+drvs.df11001
    fRRGX = D2*(VX(19)+VX(21))           !2*(drvs.df20001+drvs.df11001)
    fRZZ = (VX(22)+VX(24)+D2*VX(23))/GRD2 - (VX(8)+VX(9))/(GRD2*GRDA)  !(drvs.df1020+drvs.df0120+2*drvs.df1011)/(4*grada2)-(drvs.df1010+drvs.df1001)/(4*grada3)
    fZZZ = ((VX(25)+D3*VX(26))/(GRDA3)-D3*(VX(11)+VX(12))/(GRDA2*GRDA2)+D3*VX(3)/(GRDA3*GRDA2))/D8
!((drvs.df0030 + 3*drvs.df0021)/grada3& 
!         &-3*(drvs.df0020 + drvs.df0011)/(grada2*grada2)&
!         &+3*drvs.df0010/(grada3*grada2))/8.0
    gradY = D05*(EXPGRAD(1,IPNT,1)*GRAD(1,IPNT,1) &
         &+EXPGRAD(2,IPNT,1)*GRAD(2,IPNT,1) &
         &+EXPGRAD(3,IPNT,1)*GRAD(3,IPNT,1))
    gradZ = D05*(EXPGRAD(1,IPNT,2)*GRAD(1,IPNT,1) &
         &+EXPGRAD(2,IPNT,2)*GRAD(2,IPNT,1) &
         &+EXPGRAD(3,IPNT,2)*GRAD(3,IPNT,1))
    gradYZ = (EXPGRAD(1,IPNT,2)*EXPGRAD(1,IPNT,1) &
         &+EXPGRAD(2,IPNT,2)*EXPGRAD(2,IPNT,1) &
         &+EXPGRAD(3,IPNT,2)*EXPGRAD(3,IPNT,1))
    VXC(1,IPNT) = D2*fRRR*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) &!OK
         &+D4*(fRRZ*EXPVAL(IPNT,1)*gradZ+fRRZ*EXPVAL(IPNT,2)*gradY) & !OK
         &+D8*gradZ*gradY*fRZZ &! OK
         &+D4*gradYZ*fRZ &! OK  
         &+D4*fRRG*EXPVAL(IPNT,1)*gradZ &! OK
         &+D4*fRRG*EXPVAL(IPNT,2)*gradY+D2*fRG*gradYZ !OK

    A = D8*fZZZ*gradY*gradZ &
         & + D4*(fRZZ*EXPVAL(IPNT,1)*gradZ + fRZZ*EXPVAL(IPNT,2)*gradY) &
         & + D2*fRRZ*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) &
         & + fRRGX*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) + D4*fZZ*gradYZ
    B = D8*fZZ*gradY + D4*fRZ*EXPVAL(IPNT,1) + D2*fRG*EXPVAL(IPNT,1)
    C = D8*fZZ*gradZ + D4*fRZ*EXPVAL(IPNT,2) + D2*fRG*EXPVAL(IPNT,2)
    
    VXC(2,IPNT) = D2*A*GRAD(1,IPNT,1) + D2*B*EXPGRAD(1,IPNT,2) + D2*C*EXPGRAD(1,IPNT,1)
    VXC(3,IPNT) = D2*A*GRAD(2,IPNT,1) + D2*B*EXPGRAD(2,IPNT,2) + D2*C*EXPGRAD(2,IPNT,1)
    VXC(4,IPNT) = D2*A*GRAD(3,IPNT,1) + D2*B*EXPGRAD(3,IPNT,2) + D2*C*EXPGRAD(3,IPNT,1)
   ELSE
    VXC(IPNT,:) = 0.0d0
   ENDIF
  END DO
  CALL II_DISTGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
       & 1,NBLEN,VXC(:,:),GAO(:,:,:),DFTDATA%FKSM(:,:,1),DFTHRI)
 ENDIF
ENDIF
deallocate(EXPVAL)
deallocate(EXPGRAD)
deallocate(VXC)

END SUBROUTINE II_DFT_QUADRSPGGA

!> \brief main magnetic derivative kohn-sham matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_magderiv_kohnsham_mat(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL          :: DOGGA,DOLONDON
INTEGER          :: GRDONE,NHTYP,IDMAT

  DOGGA = DFT_ISGGA()
  DOLONDON = .TRUE.
  DFTDATA%ENERGY = 0D0
  IF(DOGGA) THEN
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL LSQUIT('II_DFT_magderiv_kohnshamGGAUNRES not implemented',lupri)
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
!             & II_DFT_magderiv_kohnshamGGAUNRES,DFTDATA)
     ELSE
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
             & II_DFT_magderiv_kohnshamGGA,DFTDATA)
     ENDIF
  ELSE 
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL LSQUIT('II_DFT_magderiv_kohnshamLDAUNRES not implemented',lupri)
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
!             & II_DFT_magderiv_kohnshamLDAUNRES,DFTDATA)
     ELSE !DEFAULT (RESTRICTED)
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
             & II_DFT_magderiv_kohnshamLDA,DFTDATA)
     END IF
  ENDIF
  !ANTISYMMETRIZE
  DO IDMAT = 1,NDMAT*3
     DO I = 1, NBAST
        DO J = 1, I-1
           AVERAG = 0.5*(DFTDATA%FKSM(J,I,IDMAT) - DFTDATA%FKSM(I,J,IDMAT))
           DFTDATA%FKSM(J,I,IDMAT) = AVERAG
           DFTDATA%FKSM(I,J,IDMAT) = -AVERAG
        END DO
        DFTDATA%FKSM(I,I,IDMAT) = 0.d0  
     END DO
  END DO
  NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
  IF(IPRINT.GE.0) WRITE(LUPRI,'(A,F20.14,A,E9.2)')&
       &     'KS electrons/energy:', DFTDATA%ELECTRONS, &
       &     ' rel.err:', (DFTDATA%ELECTRONS-NELE)/(NELE)

END SUBROUTINE II_DFT_MAGDERIV_KOHNSHAM_MAT

!> \brief magnetic derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_magderiv_kohnshamLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,MXBLLEN,COORD,WGHT,&
     & ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
Real(realk), parameter :: D4 = 4.0D0,DUMMY = 0D0
REAL(REALK) :: VXC(NBLEN),VX(5),DFTENE
INTEGER     :: I,J,IPNT

! LDA Exchange-correlation contribution to Kohn-Sham energy
DO IPNT = 1, NBLEN
!   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000   = \frac{\partial  f}{\partial \rho_{\alpha}}
      CALL dft_funcderiv1(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
      !WARNING the factor 2 is due to F=F_alpha + F_beta = 2 F_alpha 
      VXC(IPNT) = VX(1)*D4 
!   ELSE
!      VXC(IPNT) = 0.0d0
!   ENDIF
END DO
!ntypso should be 4
CALL II_DFT_distmagderiv_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     &1,NBLEN,VXC,GAO(:,:,:),DFTDATA%FKSM(:,:,1:3),COORD,DFTHRI,NTYPSO)

END SUBROUTINE II_DFT_MAGDERIV_KOHNSHAMLDA

!> \brief distribution routine for magnetic derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTMAGDERIV_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,IBLSTART,IBLEND,COEF,GAOS,EXCMAT,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in)  :: IBLSTART
!> end gridpoint
INTEGER,intent(in)  :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,3)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
REAL(REALK) :: TMP(NBLEN,NACTBAST,3),GAOMAX
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST)
REAL(REALK),allocatable :: GAORED(:,:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
REAL(REALK),pointer :: EXCRED(:,:,:)
REAL(REALK) :: Rbeta,Rgamma
real(realk),parameter :: D2=2.D0

ALLOCATE(GAORED(NBLEN,NACTBAST,4))
NK=IBLEND-IBLSTART+1
NRED = 0 
GAOMAX = 0.0d0
!        Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0d0
      DO K = IBLSTART, IBLEND
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,N)))
         GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
      ENDDO
   ENDDO
 ENDDO
ENDDO
!        Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = INXACT(J) 
   DO N=1,NTYPSO
    DO K = IBLSTART, IBLEND
       GAORED(K,NRED,N) = GAOS(K,J,N)
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT.0) THEN
   !  First half-contraction of GAO's with potential
 DO ALPHA = 1,3
  beta = betaList(ALPHA)
  gamma= gammaList(ALPHA)
  DO J=1,NRED
   DO K=IBLSTART, IBLEND
    Rgamma = COORD(gamma,K) !- Origin(gamma)       
    Rbeta = COORD(beta,K) !- Origin(beta)
    TMP(K,J,ALPHA) = coef(K)*GAORED(K,J,1+beta)*Rgamma & 
         & - coef(K)*GAORED(K,J,1+gamma)*Rbeta
   ENDDO
  ENDDO
 ENDDO
 !  Second half-contraction of GAO's
 call mem_alloc(EXCRED,NRED,NRED,3)
 DO ALPHA = 1,3
  CALL DGEMM('T','N',NRED,NRED,NK,1D0,GAORED(:,:,1),NK,&
       &     TMP(:,:,ALPHA),NK,0.0d0,EXCRED(:,:,ALPHA),NRED)
 ENDDO
!$OMP CRITICAL
 DO ALPHA=1,3
  DO JRED=1,NRED    !Jred is reduced index
   J = INXRED(JRED) !J is orbital index
   DO IRED=1,NRED   !Ired is reduced index
    I = INXRED(IRED)!I is orbital index 
    EXCMAT(I,J,ALPHA) = EXCMAT(I,J,ALPHA) + EXCRED(IRED,JRED,ALPHA)
   ENDDO
  ENDDO
 ENDDO
!$OMP END CRITICAL
 call mem_DEALLOC(EXCRED)
ENDIF
DEALLOCATE(GAORED)

END SUBROUTINE II_DFT_DISTMAGDERIV_LDA

!> \brief  magnetic derivative Kohn-sham matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_magderiv_kohnshamGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,MXBLLEN,COORD,WGHT,&
     & ENERGY,DFTDATA,RHOTHR,DFTHRI)
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
Real(realk), parameter :: D2 = 2.0D0,DUMMY = 0D0,D3 = 3.0D0,D05 = 0.5D0
REAL(REALK) :: VXC(4,NBLEN),VX(5),DFTENE,GRD,GRDA,A
INTEGER     :: I,J
DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
      CALL dft_funcderiv1(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
      VXC(1,IPNT) = D2*VX(1) 
      IF(GRD.GT.1D-40) THEN
         GRDA = D05*GRD
         A = D2*(VX(2)/GRDA + VX(3))
         VXC(2,IPNT) = A*GRAD(1,IPNT,1)
         VXC(3,IPNT) = A*GRAD(2,IPNT,1)
         VXC(4,IPNT) = A*GRAD(3,IPNT,1)
      ELSE
         VXC(2,IPNT) = 0D0
         VXC(3,IPNT) = 0D0
         VXC(4,IPNT) = 0D0
      ENDIF
   ELSE
      VXC(:,IPNT) = 0D0
   END IF
END DO
CALL II_DFT_distmagderiv_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     &1,NBLEN,VXC,GAO(:,:,:),DFTDATA%FKSM(:,:,1:3),COORD,DFTHRI,NTYPSO)

END SUBROUTINE II_DFT_MAGDERIV_KOHNSHAMGGA

!> \brief distribution routine for magnetic derivative Kohn-sham matrix GGA
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTMAGDERIV_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,IBLSTART,IBLEND,COEF,GAOS,EXCMAT,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in)  :: IBLSTART
!> end gridpoint
INTEGER,intent(in)  :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(4,NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,3)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
REAL(REALK) :: TMP(NBLEN,NACTBAST,3),GAOMAX
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
INTEGER     :: beta1,gamma1,beta2X,beta2Y,beta2Z,gamma2X,gamma2Y,gamma2Z
REAL(REALK) :: GAOGMX(NACTBAST)
REAL(REALK),allocatable :: GAORED(:,:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,COORDINATE
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
REAL(REALK),pointer :: EXCRED(:,:,:)
REAL(REALK) :: Rbeta,Rgamma
real(realk),parameter :: D2=2.D0,D05=0.5D0

ALLOCATE(GAORED(NBLEN,NACTBAST,NTYPSO))
NK=IBLEND-IBLSTART+1
NRED = 0 
GAOMAX = 0.0d0
!        Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOGMX(J) = 0.0d0
  DO N=1,NTYPSO
   DO K = IBLSTART, IBLEND
    GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,N)))
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO
!        Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = INXACT(J) 
   DO N=1,NTYPSO
    DO K = IBLSTART, IBLEND
       GAORED(K,NRED,N) = GAOS(K,J,N)
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT.0) THEN
 DO ALPHA = 1,3
  beta = betaList(ALPHA)
  gamma= gammaList(ALPHA)
  beta1 = 4+betaList(ALPHA)
  gamma1= 4+gammaList(ALPHA)
  beta2X = 7+(betaList(ALPHA)-1)*3+1
  beta2Y = 7+(betaList(ALPHA)-1)*3+2
  beta2Z = 7+(betaList(ALPHA)-1)*3+3
  gamma2X= 7+(gammaList(ALPHA)-1)*3+1
  gamma2Y= 7+(gammaList(ALPHA)-1)*3+2
  gamma2Z= 7+(gammaList(ALPHA)-1)*3+3
  DO J=1,NRED
   DO K=IBLSTART, IBLEND
    Rgamma = COORD(gamma,K) !- Origin(gamma)       
    Rbeta = COORD(beta,K) !- Origin(beta)
    TMP(K,J,ALPHA) = &
         &+(D2*coef(1,K)*Rgamma+coef(1+gamma,K))*GAORED(K,J,beta1)&
         &-(D2*coef(1,K)*Rbeta+coef(1+beta,K))*GAORED(K,J,gamma1)&
         &-(GAORED(K,J,gamma2X)*Rbeta-GAORED(K,J,beta2X)*Rgamma)*coef(2,K)&
         &-(GAORED(K,J,gamma2Y)*Rbeta-GAORED(K,J,beta2Y)*Rgamma)*coef(3,K)&
         &-(GAORED(K,J,gamma2Z)*Rbeta-GAORED(K,J,beta2Z)*Rgamma)*coef(4,K)
   ENDDO
  ENDDO
 ENDDO
 call mem_alloc(EXCRED,NRED,NRED,3)
 DO ALPHA = 1,3
    CALL DGEMM('T','N',NRED,NRED,NK,1D0,GAORED(:,:,1),NK,&
         &     TMP(:,:,ALPHA),NK,0.0d0,EXCRED(:,:,ALPHA),NRED)
 ENDDO
 DO COORDINATE=2,4
    DO ALPHA = 1,3
     beta = betaList(ALPHA)
     gamma= gammaList(ALPHA)
     beta1 = 4+betaList(ALPHA)
     gamma1= 4+gammaList(ALPHA)
     DO J=1,NRED
      DO K=IBLSTART, IBLEND
       Rgamma = COORD(gamma,K) !- Origin(gamma)       
       Rbeta = COORD(beta,K) !- Origin(beta)
       TMP(K,J,ALPHA) = (-GAORED(K,J,gamma1)*Rbeta+GAORED(K,J,beta1)*Rgamma)&
            &*coef(COORDINATE,K)
      ENDDO
     ENDDO
    ENDDO
    DO ALPHA = 1,3
       CALL DGEMM('T','N',NRED,NRED,NK,1D0,GAORED(:,:,COORDINATE),NK,&
            &     TMP(:,:,ALPHA),NK,1.0d0,EXCRED(:,:,ALPHA),NRED)
    ENDDO
 ENDDO
!$OMP CRITICAL
 DO ALPHA=1,3
  DO JRED=1,NRED     !Jred is reduced index
   J = INXRED(JRED)  !J is orbital index
   DO IRED=1,NRED    !Ired is reduced index
    I = INXRED(IRED) !I is orbital index
    EXCMAT(I,J,ALPHA) = EXCMAT(I,J,ALPHA) + EXCRED(IRED,JRED,ALPHA)
   ENDDO
  ENDDO
 ENDDO
!$OMP END CRITICAL
 call mem_DEALLOC(EXCRED)
ENDIF
DEALLOCATE(GAORED)

END SUBROUTINE II_DFT_DISTMAGDERIV_GGA

!> \brief Main driver for the calculation of magnetic derivative of the linear response
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ii_dft_magderiv_linrsp(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL          :: DOGGA,DOLONDON
INTEGER          :: GRDONE,NHTYP,IDMAT,IBMAT,nbmat

  NBMAT = DFTDATA%NBMAT
  DOGGA = DFT_ISGGA()
  DOLONDON = .TRUE.
  DFTDATA%ENERGY = 0D0
  IF(DOGGA) THEN
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL LSQUIT('II_DFT_MAGDERIV_LINRSPGGAUNRES not implemented',lupri)
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
!             & II_DFT_MAGDERIV_LINRSPGGAUNRES,DFTDATA)
     ELSE
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
             & II_DFT_MAGDERIV_LINRSPGGA,DFTDATA)
     ENDIF
  ELSE 
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL LSQUIT('II_DFT_MAGDERIV_LINRSPLDAUNRES not implemented',lupri)
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
!             & II_DFT_MAGDERIV_LINRSPLDAUNRES,DFTDATA)
     ELSE !DEFAULT (RESTRICTED)
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
             & II_DFT_MAGDERIV_LINRSPLDA,DFTDATA)
     END IF
  ENDIF
  !ANTISYMMETRIZE FKSM
  DO IBMAT = 1,NBMAT*3
     DO I = 1, NBAST
        DO J = 1, I-1
           AVERAG = 0.5*(DFTDATA%FKSM(J,I,IBMAT) - DFTDATA%FKSM(I,J,IBMAT))
           DFTDATA%FKSM(J,I,IBMAT) =  AVERAG
           DFTDATA%FKSM(I,J,IBMAT) = -AVERAG
        END DO
        DFTDATA%FKSM(I,I,IBMAT) = 0.d0  
     END DO
  END DO
  !SYMMETRIZE FKSMS
  DO IBMAT = 1,NBMAT*3
     DO I = 1, NBAST
        DO J = 1, I-1
           AVERAG = 0.5*(DFTDATA%FKSMS(J,I,IBMAT) + DFTDATA%FKSMS(I,J,IBMAT))
           DFTDATA%FKSMS(J,I,IBMAT) =  AVERAG
           DFTDATA%FKSMS(I,J,IBMAT) =  AVERAG
        END DO
     END DO
  END DO
  NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
  IF(IPRINT.GE.0) WRITE(LUPRI,'(A,F20.14,A,E9.2)')&
       &     'KS electrons:', DFTDATA%ELECTRONS,&
       &     ' rel.err:', (DFTDATA%ELECTRONS-NELE)/(NELE)

END SUBROUTINE II_DFT_MAGDERIV_LINRSP

!> \brief magnetic derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_MAGDERIV_LINRSPLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,MXBLLEN,COORD,WGHT,&
     & ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
Real(realk), parameter :: D4 = 4.0D0,DUMMY = 0D0
REAL(REALK) :: fRR,VX(14)
INTEGER     :: I,J,NBMAT,ibmat,N,NRED,IPNT
LOGICAL     :: DOCALC,dosympart
REAL(REALK),parameter :: D2=2.d0
REAL(REALK),allocatable :: EXPVAL(:,:),VXC(:,:)
REAL(REALK),allocatable :: EXPGRAD(:,:,:),VXC2(:,:,:)

NBMAT = DFTDATA%NBMAT
dosympart = DFTDATA%dosympart
allocate(EXPVAL(NBLEN,NBMAT))
allocate(EXPGRAD(NBLEN,3,NBMAT))!magnetic gradients of EXPVAL=BMAT*chi(r)*chi(r)
allocate(VXC(NBLEN,NBMAT))
allocate(VXC2(NBLEN,NBMAT,3))
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
!get expectation value of BMAT
 call II_get_magderiv_expectationval_lda(LUPRI,NTYPSO,NBLEN,NBLOCKS,BLOCKS,&
        & INXACT,Nactbast,NBAST,GAO,COORD,EXPVAL,EXPGRAD,DFTDATA%BMAT,nBMAT,&
        & NRED,DFTHRI,dosympart)
 IF(NRED.GT.0)THEN
  VXC2 = 0.0d0
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
    !get the functional derivatives 
    CALL dft_funcderiv2(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
    !fRR = 0.5*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
    fRR = VX(6) + VX(7)     
     !THE factor 0.5 IN fRR cancels a factor 2.0 in VXC so  
    DO IBMAT = 1,NBMAT
     VXC(IPNT,IBMAT) = D2*fRR*EXPVAL(IPNT,IBMAT)!D4*fRR*EXPVAL(IPNT,IBMAT)
    ENDDO
    DO N=1,3
     DO IBMAT = 1,NBMAT
      VXC2(IPNT,IBMAT,N)=D2*fRR*EXPGRAD(IPNT,N,IBMAT)
     ENDDO
    ENDDO
   ELSE
    VXC(IPNT,:) = 0.0d0
   ENDIF
  END DO
  CALL II_DFT_DISTMAGDERIV_linrsp_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
    &Nactbast,NBAST,1,NBLEN,VXC,VXC2,NBMAT,GAO,DFTDATA%FKSM,DFTDATA%FKSMS,COORD,DFTHRI,NTYPSO,dosympart)
 ENDIF
ENDIF
deallocate(EXPVAL)
deallocate(EXPGRAD)
deallocate(VXC)
deallocate(VXC2)

END SUBROUTINE II_DFT_MAGDERIV_LINRSPLDA

!> \brief computes the expectation value and the magnetic derivative of it
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_magderiv_expectationval_lda(LUPRI,NTYPSO,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,GAOS,COORD,EXPVAL,EXPGRAD,BMAT,nBmat,NRED,DFTHRI,DOSYMPART)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> The expectationvalue
REAL(REALK),intent(inout) :: EXPVAL(NBLEN,NBMAT)
!> The magnetic gradient of expectationvalue
REAL(REALK),intent(inout) :: EXPGRAD(NBLEN,3,NBMAT)
!> B matrix some perturbed density matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> Number of reduced orbitals
INTEGER,intent(inout)  :: NRED
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> do the symmetric part 
LOGICAL,intent(in) :: DOSYMPART
!
REAL(REALK) :: TMP(NBLEN,NACTBAST,4),GAOMAX,BMAX
REAL(REALK) :: GAORED(NBLEN,NACTBAST,NTYPSO),GAOGMX(NACTBAST),X,Y,Z
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,INXRED(NACTBAST)
INTEGER     :: IRED,JRED,IORB,JORB,gabX,gabY,gabZ,N,IBMAT
REAL(REALK),allocatable :: BRED(:,:,:)

!IF(NBMAT .GT.1)CALL LSQUIT('II_get_expectationval_lda does not work correctly for more than 1 nbmat',lupri)
NRED = 0
GAOMAX = 0.0d0
! Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0d0
      DO K = 1,NBLEN
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,N)))
         GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
      ENDDO
   ENDDO
 ENDDO
ENDDO

! Set up maximum response-vector elements
BMAX = 0.0d0
DO IBMAT=1,NBMAT
DO JBL=1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL) !J is active index 
      JORB = INXACT(J) !JORB is orbital index 
      DO IBL=1, NBLOCKS  
         DO I = BLOCKS(1,IBL), BLOCKS(2,IBL) !I is active index 
            IORB = INXACT(I) !IORB is orbital index 
            BMAX = MAX(BMAX,ABS(BMAT(IORB,JORB,IBMAT)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
ENDDO

! Set up reduced Gaussian AO's
NRED=0
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX*BMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO N=1,NTYPSO
          DO K = 1, NBLEN
             GAORED(K,NRED,N) = GAOS(K,J,N)
          ENDDO 
         ENDDO
      ENDIF
   ENDDO
ENDDO

IF (NRED.GT.0) THEN
   ! Set up reduced density-matrix
 ALLOCATE(BRED(NRED,NRED,NBMAT))
 DO IBMAT=1,NBMAT
  DO JRED=1,NRED     !Jred is reduced index 
   J = INXRED(JRED)  !J is orbital index 
   DO IRED=1,NRED    !Ired is reduced index 
    I = INXRED(IRED) !I is orbital index 
    BRED(IRED,JRED,IBMAT) = BMAT(I,J,IBMAT)
   ENDDO
  ENDDO
 ENDDO

   ! First half-contraction of Gaussian AO with density-matrix
 DO IBMAT=1,NBMAT
    CALL DGEMM("N","N",NBLEN,NRED,NRED,1D0,GAORED(:,:,1),NBLEN,&
        &     BRED(:,:,IBMAT),NRED,0.0d0,TMP(:,:,1),NBLEN)
    ! Second half-contraction 
    DO K = 1, NBLEN
       EXPVAL(K,IBMAT) = GAORED(K,1,1)*TMP(K,1,1)
    END DO
    DO I = 2, NRED
       DO K = 1, NBLEN
          EXPVAL(K,IBMAT) = EXPVAL(K,IBMAT) + GAORED(K,I,1)*TMP(K,I,1)
       END DO
    END DO
 ENDDO

 IF(DOSYMPART)THEN 
    !IF MATRIX B IS SYMMETRIC THIS WILL BE ZERO
    !SO WE ONLY CALCULATE THIS IF THE B MATRIX IS NONSYMMETRIC
    !RESULTING IN A SYMMETRIC FINAL CONTRIBUTION 
   gabX = 2
   gabY = 3
   gabZ = 4
   DO IBMAT=1,NBMAT
   ! First half-contraction of Gaussian AO with density-matrix
    !TMP(K,J,N) = GAORED(K,I,N)*BRED(I,J) - TMP(K,J,1) is already built
    DO N=2,4
     CALL DGEMM("N","N",NBLEN,NRED,NRED,1D0,GAORED(:,:,N),NBLEN,&
          &     BRED(:,:,IBMAT),NRED,0.0d0,TMP(:,:,N),NBLEN)
    ENDDO
   ! Second half-contraction of Gaussian AOs
    J=1
    DO K = 1, NBLEN
     Y = COORD(2,K) !- Origin(2)      
     Z = COORD(3,K) !- Origin(3)      
     EXPGRAD(K,1,IBMAT) = &
          & +GAORED(K,J,1)*(Z*TMP(K,J,gabY)-Y*TMP(K,J,gabZ))&
          & +TMP(K,J,1)*(Y*GAORED(K,J,gabZ)-Z*GAORED(K,J,gabY))
    END DO
    DO J = 2, NRED
     DO K = 1, NBLEN
        Y = COORD(2,K) !- Origin(2)      
        Z = COORD(3,K) !- Origin(3)      
        EXPGRAD(K,1,IBMAT) = EXPGRAD(K,1,IBMAT)&
             & +GAORED(K,J,1)*(Z*TMP(K,J,gabY)-Y*TMP(K,J,gabZ))&
             & +TMP(K,J,1)*(Y*GAORED(K,J,gabZ)-Z*GAORED(K,J,gabY))
     END DO
    ENDDO
    J=1
    DO K = 1, NBLEN
     X = COORD(1,K) !- Origin(1)      
     Z = COORD(3,K) !- Origin(3)      
     EXPGRAD(K,2,IBMAT) = &
          & +GAORED(K,J,1)*(X*TMP(K,J,gabZ)-Z*TMP(K,J,gabX))&
          & +TMP(K,J,1)*(Z*GAORED(K,J,gabX)-X*GAORED(K,J,gabZ))
    END DO
    DO J = 2, NRED
     DO K = 1, NBLEN
        X = COORD(1,K) !- Origin(1)      
        Z = COORD(3,K) !- Origin(3)      
        EXPGRAD(K,2,IBMAT) = EXPGRAD(K,2,IBMAT)&
             & +GAORED(K,J,1)*(X*TMP(K,J,gabZ)-Z*TMP(K,J,gabX))&
             & +TMP(K,J,1)*(Z*GAORED(K,J,gabX)-X*GAORED(K,J,gabZ))
     END DO
    ENDDO
    J=1
    DO K = 1, NBLEN
     X = COORD(1,K) !- Origin(1)      
     Y = COORD(2,K) !- Origin(2)      
     EXPGRAD(K,3,IBMAT) = &
          & +GAORED(K,J,1)*(Y*TMP(K,J,gabX)-X*TMP(K,J,gabY))&
          & +TMP(K,J,1)*(X*GAORED(K,J,gabY)-Y*GAORED(K,J,gabX))
    END DO
    DO J = 2, NRED
     DO K = 1, NBLEN
        X = COORD(1,K) !- Origin(1)      
        Y = COORD(2,K) !- Origin(2)      
        EXPGRAD(K,3,IBMAT) = EXPGRAD(K,3,IBMAT)&
             & +GAORED(K,J,1)*(Y*TMP(K,J,gabX)-X*TMP(K,J,gabY))&
             & +TMP(K,J,1)*(X*GAORED(K,J,gabY)-Y*GAORED(K,J,gabX))
     END DO
    END DO
   ENDDO
 ENDIF
 DEALLOCATE(BRED)
ENDIF

END SUBROUTINE II_GET_MAGDERIV_EXPECTATIONVAL_LDA

!> \brief A distribution routine for the magnetic derivative of the linear response
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTMAGDERIV_linrsp_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,IBLSTART,IBLEND,COEF,COEF2,NBMAT,GAOS,EXCMAT,EXCMATSYM,COORD,DFTHRI,NTYPSO,dosympart)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in) :: IBLSTART
!> end gridpoint
INTEGER,intent(in) :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(NBLEN,NBMAT)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(NBLEN,NBMAT,3)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The magnetic derivative LDA Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,3*NBMAT)
!> The sym part magnetic derivative LDA Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMATSYM(NBAST,NBAST,3*NBMAT)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!> do the symmetric part 
LOGICAL,intent(in) :: DOSYMPART
!
REAL(REALK) :: TMP(NBLEN,NACTBAST,3),GAOMAX
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST)
REAL(REALK),allocatable :: GAORED(:,:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,IBMAT
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
REAL(REALK),pointer :: EXCRED(:,:,:)
REAL(REALK) :: Rbeta,Rgamma
real(realk),parameter :: D2=2.D0,D05=0.5d0

ALLOCATE(GAORED(NBLEN,NACTBAST,4))
NK=IBLEND-IBLSTART+1
GAOMAX = 0.0d0
!        Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOGMX(J) = 0.0d0
  DO K = IBLSTART, IBLEND
     GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,1)))
     DO N=1,NTYPSO
        GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
     ENDDO
  ENDDO
 ENDDO
ENDDO
!        Set up reduced Gaussian AO's
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = INXACT(J) 
   DO N=1,NTYPSO
    DO K = IBLSTART, IBLEND
       GAORED(K,NRED,N) = GAOS(K,J,N)
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT.0) THEN
   !  First half-contraction of GAO's with potential
 call mem_alloc(EXCRED,NRED,NRED,3)
 DO IBMAT = 1,NBMAT
  DO ALPHA = 1,3
   beta = betaList(ALPHA)
   gamma= gammaList(ALPHA)
   DO J=1,NRED
    DO K=IBLSTART, IBLEND
     Rgamma = COORD(gamma,K) !- Origin(gamma)       
     Rbeta = COORD(beta,K) !- Origin(beta)
     TMP(K,J,ALPHA) = &
          &-coef(K,IBMAT)*GAORED(K,J,1+beta)*Rgamma &
          &+coef(K,IBMAT)*GAORED(K,J,1+gamma)*Rbeta 
    ENDDO
   ENDDO
  ENDDO
  !  Second half-contraction of GAO's
  DO ALPHA = 1,3
   CALL DGEMM('T','N',NRED,NRED,NK,1D0,GAORED(:,:,1),NK,&
        &     TMP(:,:,ALPHA),NK,0.0d0,EXCRED(:,:,ALPHA),NRED)
  ENDDO
!$OMP CRITICAL
  DO ALPHA=1,3
   N = ALPHA+(IBMAT-1)*NBMAT
   DO JRED=1,NRED      !Jred is reduced index 
    J = INXRED(JRED)   !J is orbital index 
    DO IRED=1,NRED     !Ired is reduced index 
     I = INXRED(IRED)  !I is orbital index 
     EXCMAT(I,J,N) = EXCMAT(I,J,N) + EXCRED(IRED,JRED,ALPHA)
    ENDDO
   ENDDO
  ENDDO
!$OMP END CRITICAL
 ! second part magnetic differentiation on Bmat part give a sym mat and will only contribute if BMAT is non symmetric
  IF(DOSYMPART)THEN
   DO ALPHA = 1,3
    DO J=1,NRED
     DO K=IBLSTART, IBLEND
      TMP(K,J,ALPHA) = -GAORED(K,J,1)*coef2(K,IBMAT,ALPHA)
     ENDDO
    ENDDO
   ENDDO
   !  Second half-contraction of GAO's
   DO ALPHA = 1,3
     CALL DGEMM('T','N',NRED,NRED,NK,1D0,GAORED(:,:,1),NK,&
          &     TMP(:,:,ALPHA),NK,0.0d0,EXCRED(:,:,ALPHA),NRED)
   ENDDO
!$OMP CRITICAL
   DO ALPHA=1,3
    N = ALPHA+(IBMAT-1)*NBMAT
    DO JRED=1,NRED     !Jred is reduced index 
     J = INXRED(JRED)  !J is orbital index  
     DO IRED=1,NRED    !Ired is reduced index  
      I = INXRED(IRED) !I is orbital index 
      EXCMATSYM(I,J,N) = EXCMATSYM(I,J,N) + EXCRED(IRED,JRED,ALPHA)
     ENDDO
    ENDDO
   ENDDO
!$OMP END CRITICAL
  ENDIF
 ENDDO
 call mem_DEALLOC(EXCRED)
ENDIF
DEALLOCATE(GAORED)

END SUBROUTINE II_DFT_DISTMAGDERIV_linrsp_LDA

!> \brief magnetic derivative linrsp matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_MAGDERIV_LINRSPGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,MXBLLEN,COORD,WGHT,&
     & ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
INTEGER     :: I,J,NBMAT,ibmat,N,NRED,IPNT
REAL(REALK) :: VX(14),GRD,GRDA,MIXEDGRDA,MIXEDGRDAMAG(3)
LOGICAL     :: DOCALC,dosympart
REAL(REALK),allocatable :: VXC1(:,:),VXC1MAG(:,:,:)
REAL(REALK),allocatable :: EXPGRAD(:,:,:,:),VXC2(:,:,:),VXC2MAG(:,:,:,:)
Real(realk), parameter :: D4 = 4.0D0, DUMMY = 0D0,D05 = 0.5D0
Real(realk), parameter :: D2 = 2.0D0, D8 = 8.0D0,D025 = 0.25D0
REAL(REALK) :: fR,fZ,fRR,fRZ,fZZ,fRG,fZG,fGG,fG,A,B,AMAG(3)
NBMAT = DFTDATA%NBMAT
dosympart = DFTDATA%dosympart
allocate(EXPGRAD(NBLEN,4,4,NBMAT))!mixed geo,magn gradients of EXPVAL
allocate(VXC1(NBLEN,NBMAT))
allocate(VXC1MAG(NBLEN,NBMAT,3))
allocate(VXC2(NBLEN,NBMAT,3))
allocate(VXC2MAG(NBLEN,NBMAT,3,3))
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
!get expectation value of BMAT
 call II_get_magderiv_expectationval_gga(LUPRI,NTYPSO,NBLEN,NBLOCKS,BLOCKS,&
        & INXACT,Nactbast,NBAST,GAO,COORD,EXPGRAD,DFTDATA%BMAT,nBMAT,&
        & NRED,DFTHRI,dosympart)
 IF(NRED.GT.0)THEN
  VXC2 = 0.0d0
  DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
    IF(GRD.LT.1D-40) GRD = 1D-40
    GRDA = D05*GRD
    !get the functional derivatives 
    CALL dft_funcderiv2(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
    fR  = D05*(VX(1) + VX(2))   !0.5*(drvs.df1000 + drvs.df0100);
    fZ  = VX(3)                    !drvs.df0010;
    fRR = D05*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
    fRZ = D05*(VX(8) + VX(9))   !0.5*(drvs.df1010 + drvs.df1001);
    fZZ = D05*(VX(11) + VX(12)) !0.5*(drvs.df0020 + drvs.df0011);
    fRG = D05*VX(10)             !0.5*drvs.df10001;   
    fZG = D05*VX(13)             !0.5*drvs.df00101; 
    fGG = D025*VX(14)            !0.25*drvs.df00002; 
    fG  = D05*VX(5)               !0.5*drvs.df00001;  
    DO IBMAT = 1,NBMAT
     MIXEDGRDA = (EXPGRAD(IPNT,2,1,IBMAT)*GRAD(1,IPNT,1)&
          &+EXPGRAD(IPNT,3,1,IBMAT)*GRAD(2,IPNT,1)&
          &+EXPGRAD(IPNT,4,1,IBMAT)*GRAD(3,IPNT,1))
     VXC1(IPNT,IBMAT) =D4*fRR*EXPGRAD(IPNT,1,1,IBMAT)+D4*(fRZ/GRD+fRG)*MIXEDGRDA
     !magnetic differentiation of VXC1
     DO N=1,3
        MIXEDGRDAMAG(N) = (EXPGRAD(IPNT,2,N+1,IBMAT)*GRAD(1,IPNT,1)&
             &+EXPGRAD(IPNT,3,N+1,IBMAT)*GRAD(2,IPNT,1)&
             &+EXPGRAD(IPNT,4,N+1,IBMAT)*GRAD(3,IPNT,1))
     ENDDO
     DO N=1,3
      VXC1MAG(IPNT,IBMAT,N)=D4*(fRZ/GRD+fRG)*MIXEDGRDAMAG(N)+D4*fRR*EXPGRAD(IPNT,1,1+N,IBMAT)
     ENDDO
     A = D8*((fRZ/GRD + fRG)*EXPGRAD(IPNT,1,1,IBMAT)&
          & + (((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)*MIXEDGRDA)
     B= D8*(fZ/GRD + fG)
     VXC2(IPNT,IBMAT,1) = A*GRAD(1,IPNT,1)+B*EXPGRAD(IPNT,2,1,IBMAT)
     VXC2(IPNT,IBMAT,2) = A*GRAD(2,IPNT,1)+B*EXPGRAD(IPNT,3,1,IBMAT)
     VXC2(IPNT,IBMAT,3) = A*GRAD(3,IPNT,1)+B*EXPGRAD(IPNT,4,1,IBMAT)
     !magnetic differentiation of VXC2
     DO N=1,3
        AMAG(N) = D8*((fRZ/GRD + fRG)*EXPGRAD(IPNT,1,N+1,IBMAT)&
             & + (((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)*MIXEDGRDAMAG(N))
     ENDDO
     DO N=1,3
      VXC2MAG(IPNT,IBMAT,1,N)=AMAG(N)*GRAD(1,IPNT,1)+B*EXPGRAD(IPNT,2,1+N,IBMAT)
      VXC2MAG(IPNT,IBMAT,2,N)=AMAG(N)*GRAD(2,IPNT,1)+B*EXPGRAD(IPNT,3,1+N,IBMAT)
      VXC2MAG(IPNT,IBMAT,3,N)=AMAG(N)*GRAD(3,IPNT,1)+B*EXPGRAD(IPNT,4,1+N,IBMAT)
      ENDDO
    ENDDO
   ELSE
    VXC1(IPNT,:) = 0.0d0
    VXC1MAG(IPNT,:,:) = 0.0d0
    VXC2(IPNT,:,:) = 0.0d0
    VXC2MAG(IPNT,:,:,:) = 0.0d0
   ENDIF
  END DO
  CALL II_DFT_DISTMAGDERIV_linrsp_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
    &Nactbast,NBAST,1,NBLEN,VXC1,VXC1MAG,VXC2,VXC2MAG,NBMAT,GAO,&
    &DFTDATA%FKSM,DFTDATA%FKSMS,COORD,DFTHRI,NTYPSO,dosympart)
 ENDIF
ENDIF
deallocate(EXPGRAD)
deallocate(VXC1)
deallocate(VXC2)
deallocate(VXC1MAG)
deallocate(VXC2MAG)

END SUBROUTINE II_DFT_MAGDERIV_LINRSPGGA

!> \brief distribution routine for magnetic derivative linrsp matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTMAGDERIV_linrsp_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,IBLSTART,IBLEND,COEF1,COEF1MAG,COEF2,COEF2MAG,NBMAT,GAOS,&
     &EXCMAT,EXCMATSYM,COORD,DFTHRI,NTYPSO,dosympart)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in) :: IBLSTART
!> end gridpoint
INTEGER,intent(in) :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(NBLEN,NBMAT)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(NBLEN,NBMAT,3)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1MAG(NBLEN,NBMAT,3)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2MAG(NBLEN,NBMAT,3,3)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The magnetic derivative LDA Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,3*NBMAT)
!> The sym part magnetic derivative LDA Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMATSYM(NBAST,NBAST,3*NBMAT)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!> do the symmetric part 
LOGICAL,intent(in) :: DOSYMPART
!
REAL(REALK) :: TMP(NBLEN,NACTBAST,3),GAOMAX
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED,beta1,gamma1
REAL(REALK) :: GAOGMX(NACTBAST)
REAL(REALK),allocatable :: GAORED(:,:,:)
INTEGER     :: gab1beta,gab1gamma,gab2Xbeta,gab2Xgamma,gab2Ybeta,gab2Ygamma
INTEGER     :: gab2Zbeta,gab2Zgamma,COORDINATE
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,IBMAT
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
REAL(REALK),pointer :: EXCRED(:,:,:)
REAL(REALK) :: Rbeta,Rgamma
real(realk),parameter :: D2=2.D0,D05=0.5d0

ALLOCATE(GAORED(NBLEN,NACTBAST,NTYPSO))
NK=IBLEND-IBLSTART+1
GAOMAX = 0.0d0
!        Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOGMX(J) = 0.0d0
  DO K = IBLSTART, IBLEND
     GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,1)))
     DO N=1,NTYPSO
        GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
     ENDDO
  ENDDO
 ENDDO
ENDDO
!        Set up reduced Gaussian AO's
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
!  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = INXACT(J) 
   DO N=1,NTYPSO
    DO K = IBLSTART, IBLEND
       GAORED(K,NRED,N) = GAOS(K,J,N)
    ENDDO
   ENDDO
!  ENDIF
 ENDDO
ENDDO

IF (NRED.GT.0) THEN
!--------------------------------------------------
!   First term containing differentiation on Omega
!-------------------------------------------------
   !  First half-contraction of GAO's with potential
 call mem_alloc(EXCRED,NRED,NRED,3)
 DO IBMAT = 1,NBMAT
  DO ALPHA = 1,3
   beta = betaList(ALPHA)
   gamma= gammaList(ALPHA)
   beta1 = 4+beta
   gamma1= 4+gamma
   DO J=1,NRED
    DO K=IBLSTART, IBLEND
     Rgamma = COORD(gamma,K) !- Origin(gamma)       
     Rbeta = COORD(beta,K) !- Origin(beta)
     TMP(K,J,ALPHA) = &
          &-coef1(K,IBMAT)*GAORED(K,J,beta1)*Rgamma &
          &+coef1(K,IBMAT)*GAORED(K,J,gamma1)*Rbeta 
    ENDDO
   ENDDO
  ENDDO
  !  Second half-contraction of GAO's
  DO ALPHA = 1,3
   CALL DGEMM('T','N',NRED,NRED,NK,1D0,GAORED(:,:,1),NK,&
        &     TMP(:,:,ALPHA),NK,0.0d0,EXCRED(:,:,ALPHA),NRED)
  ENDDO
!$OMP CRITICAL
  DO ALPHA=1,3
   N = ALPHA+(IBMAT-1)*3
   DO JRED=1,NRED      !Jred is reduced index 
    J = INXRED(JRED)   !J is orbital index 
    DO IRED=1,NRED     !Ired is reduced index 
     I = INXRED(IRED)  !I is orbital index 
     EXCMAT(I,J,N) = EXCMAT(I,J,N) + EXCRED(IRED,JRED,ALPHA)
    ENDDO
   ENDDO
  ENDDO
!$OMP END CRITICAL
!------------------------------------------------------------
!   Second term containing differentiation on Bmat(LDA part) 
!------------------------------------------------------------
 ! second part magnetic differentiation on Bmat part give a sym mat and will only contribute if BMAT is non symmetric
  IF(DOSYMPART)THEN
   DO ALPHA = 1,3
    DO J=1,NRED
     DO K=IBLSTART, IBLEND
      TMP(K,J,ALPHA) = -GAORED(K,J,1)*coef1MAG(K,IBMAT,ALPHA)
     ENDDO
    ENDDO
   ENDDO
   !  Second half-contraction of GAO's
   DO ALPHA = 1,3
     CALL DGEMM('T','N',NRED,NRED,NK,1D0,GAORED(:,:,1),NK,&
          &     TMP(:,:,ALPHA),NK,0.0d0,EXCRED(:,:,ALPHA),NRED)
   ENDDO
!$OMP CRITICAL
   DO ALPHA=1,3
    N = ALPHA+(IBMAT-1)*NBMAT
    DO JRED=1,NRED       !Jred is reduced index 
     J = INXRED(JRED)    !J is orbital index 
     DO IRED=1,NRED      !Ired is reduced index 
      I = INXRED(IRED)   !I is orbital index 
      EXCMATSYM(I,J,N) = EXCMATSYM(I,J,N) + EXCRED(IRED,JRED,ALPHA)
     ENDDO
    ENDDO
   ENDDO
!$OMP END CRITICAL
  ENDIF
!-------------------------------------------------------
!   Third term containing differentiation on Nabla*Omega
!-------------------------------------------------------
  DO ALPHA = 1,3
   beta = betaList(ALPHA)
   gamma= gammaList(ALPHA)
   gab1beta = 4+beta
   gab1gamma= 4+gamma
   gab2Xbeta = 5+3*beta
   gab2Xgamma = 5+3*gamma
   gab2Ybeta = 6+3*beta
   gab2Ygamma = 6+3*gamma
   gab2Zbeta = 7+3*beta
   gab2Zgamma = 7+3*gamma
   DO J=1,NRED
    DO K=IBLSTART, IBLEND
     Rgamma = COORD(gamma,K) !- Origin(gamma)       
     Rbeta = COORD(beta,K) !- Origin(beta)
     TMP(K,J,ALPHA) = &
         &-D05*coef2(K,IBMAT,gamma)*GAORED(K,J,gab1beta)+D05*coef2(K,IBMAT,beta)*GAORED(K,J,gab1gamma)&
         &+D05*(GAORED(K,J,gab2Xgamma)*Rbeta-GAORED(K,J,gab2Xbeta)*Rgamma)*coef2(K,IBMAT,1)&
         &+D05*(GAORED(K,J,gab2Ygamma)*Rbeta-GAORED(K,J,gab2Ybeta)*Rgamma)*coef2(K,IBMAT,2)&
         &+D05*(GAORED(K,J,gab2Zgamma)*Rbeta-GAORED(K,J,gab2Zbeta)*Rgamma)*coef2(K,IBMAT,3)
    ENDDO
   ENDDO
  ENDDO
  DO ALPHA = 1,3
     CALL DGEMM('T','N',NRED,NRED,NK,1D0,GAORED(:,:,1),NK,&
          &     TMP(:,:,ALPHA),NK,0.0d0,EXCRED(:,:,ALPHA),NRED)
  ENDDO
  DO COORDINATE=1,3
   DO ALPHA = 1,3
     beta = betaList(ALPHA)
     gamma= gammaList(ALPHA)
     gab1beta = 4+beta
     gab1gamma= 4+gamma
     DO J=1,NRED
      DO K=IBLSTART, IBLEND
       Rgamma = COORD(gamma,K) !- Origin(gamma)       
       Rbeta = COORD(beta,K) !- Origin(beta)
       TMP(K,J,ALPHA) = -D05*(-GAORED(K,J,gab1gamma)*Rbeta+GAORED(K,J,gab1beta)*Rgamma)&
            &*coef2(K,IBMAT,COORDINATE)
      ENDDO
     ENDDO
   ENDDO
   DO ALPHA = 1,3
      CALL DGEMM('T','N',NRED,NRED,NK,1D0,GAORED(:,:,1+COORDINATE),NK,&
           &     TMP(:,:,ALPHA),NK,1.0d0,EXCRED(:,:,ALPHA),NRED)
   ENDDO
  ENDDO
!$OMP CRITICAL
  DO ALPHA=1,3
   DO JRED=1,NRED      !Jred is reduced index 
    J = INXRED(JRED)   !J is orbital index    
    DO IRED=1,NRED     !Ired is reduced index 
     I = INXRED(IRED)  !I is orbital index   
     EXCMAT(I,J,ALPHA) = EXCMAT(I,J,ALPHA) + EXCRED(IRED,JRED,ALPHA)
    ENDDO
   ENDDO
  ENDDO
!$OMP END CRITICAL
!-------------------------------------------------------
!   Fourth term containing differentiation on Bmat(GGA)
!-------------------------------------------------------
!  First half-contraction of GAO's with potential
  DO ALPHA=1,3
   DO J=1,NRED
      DO K=IBLSTART, IBLEND
         TMP(K,J,ALPHA) =  -coef2MAG(K,IBMAT,1,ALPHA)*GAORED(K,J,2)&
              &   - coef2MAG(K,IBMAT,2,ALPHA)*GAORED(K,J,3)&
              &   - coef2MAG(K,IBMAT,3,ALPHA)*GAORED(K,J,4)
      ENDDO
   ENDDO
  ENDDO
!  Second half-contraction of GAO's with potential
  DO ALPHA=1,3
   CALL DGEMM('T','N',NRED,NRED,NK,1.0d0,&
        &                GAORED(:,:,1),NK,TMP(:,:,ALPHA),NK,0.0d0,&
        &                EXCRED(:,:,ALPHA),NRED)
  ENDDO
!  Distribute contributions to KS-matrix
!$OMP CRITICAL
  DO ALPHA=1,3
   DO JRED=1,NRED          !Jred is reduced index   
      J = INXRED(JRED)     !J is orbital index   
      DO IRED=1,NRED       !Ired is reduced index   
         I = INXRED(IRED)  !I is orbital index   
         EXCMATSYM(I,J,ALPHA) = EXCMATSYM(I,J,ALPHA) + EXCRED(IRED,JRED,ALPHA)
      ENDDO
   ENDDO
  ENDDO
!$OMP END CRITICAL
 ENDDO
call mem_DEALLOC(EXCRED)
ENDIF
DEALLOCATE(GAORED)

END SUBROUTINE II_DFT_DISTMAGDERIV_linrsp_GGA

!> \brief computes the expectation value(expval), geometrical gradient of expval and magnetic derivative of both.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_magderiv_expectationval_gga(LUPRI,NTYPSO,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,GAOS,COORD,EXPGRAD,BMAT,nBmat,NRED,DFTHRI,DOSYMPART)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> expectation value and geometrical and magnetic gradient
REAL(REALK),intent(inout) :: EXPGRAD(NBLEN,4,4,NBMAT)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> Number of reduced orbitals
INTEGER,intent(inout)  :: NRED
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> do the symmetric part 
LOGICAL,intent(in) :: DOSYMPART
!
REAL(REALK) :: TMP(NBLEN,NACTBAST,4),GAOMAX,BMAX,BREDIJ,BREDJI
REAL(REALK) :: GAORED(NBLEN,NACTBAST,NTYPSO),GAOGMX(NACTBAST),X,Y,Z,Rgamma,Rbeta
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,INXRED(NACTBAST)
INTEGER     :: IRED,JRED,IORB,JORB,gabX,gabY,gabZ,N,IBMAT,gabx1,gaby1,gabz1,M
INTEGER     :: alpha,beta,gamma,gab1beta,gab1gamma,gab2Xbeta,gab2Xgamma
INTEGER     :: gab2Ybeta,gab2Ygamma,gab2Zbeta,gab2Zgamma
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
REAL(REALK),allocatable :: BRED(:,:,:)
REAL(REALK),parameter :: zero=0.d0,D05=0.5d0

!IF(NBMAT .GT.1)CALL LSQUIT('II_get_expectationval_lda does not work correctly for more than 1 nbmat',lupri)
GAOMAX = 0.0d0
! Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOGMX(J) = 0.0d0
  DO K = 1,NBLEN
   GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,1)))
   DO N=1,NTYPSO
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO

! Set up maximum response-vector elements
BMAX = 0.0d0
DO IBMAT=1,NBMAT
DO JBL=1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL) !J is active index
      JORB = INXACT(J) !Jorb is orbital index
      DO IBL=1, NBLOCKS
         DO I = BLOCKS(1,IBL), BLOCKS(2,IBL) !I is active index
            IORB = INXACT(I) !Iorb is orbital index
            BMAX = MAX(BMAX,ABS(BMAT(IORB,JORB,IBMAT)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
ENDDO

! Set up reduced Gaussian AO's
NRED=0
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
!      IF (GAOGMX(J)*GAOMAX*BMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO N=1,NTYPSO
          DO K = 1, NBLEN
             GAORED(K,NRED,N) = GAOS(K,J,N)
          ENDDO 
         ENDDO
!      ENDIF
   ENDDO
ENDDO

IF (NRED.GT.0) THEN
   ! Set up reduced density-matrix
 ALLOCATE(BRED(NRED,NRED,NBMAT))
 DO IBMAT=1,NBMAT
  DO JRED=1,NRED     !Jred is reduced index
   J = INXRED(JRED)  !J is orbital index
   DO IRED=1,NRED    !Ired is reduced index
    I = INXRED(IRED) !I is orbital index
    BRED(IRED,JRED,IBMAT) = BMAT(I,J,IBMAT)
   ENDDO
  ENDDO
 ENDDO

!-------------------------------------------------------------
!  The expectation value
!-------------------------------------------------------------
   ! First half-contraction of Gaussian AO with density-matrix
 DO IBMAT=1,NBMAT
    CALL DGEMM("N","N",NBLEN,NRED,NRED,1D0,GAORED(:,:,1),NBLEN,&
        &     BRED(:,:,IBMAT),NRED,0.0d0,TMP(:,:,1),NBLEN)
    ! Second half-contraction 
    DO K = 1, NBLEN
       EXPGRAD(K,1,1,IBMAT) = GAORED(K,1,1)*TMP(K,1,1)
    END DO
    DO I = 2, NRED
       DO K = 1, NBLEN
          EXPGRAD(K,1,1,IBMAT) = EXPGRAD(K,1,1,IBMAT) + GAORED(K,I,1)*TMP(K,I,1)
       END DO
    END DO
 ENDDO

!-------------------------------------------------------------
!  The magnetic derivative of the expectation value
!-------------------------------------------------------------
 IF(DOSYMPART)THEN 
    !IF MATRIX B IS SYMMETRIC THIS WILL BE ZERO
    !SO WE ONLY CALCULATE THIS IF THE B MATRIX IS NONSYMMETRIC
    !RESULTING IN A SYMMETRIC FINAL CONTRIBUTION 
   gabX = 5
   gabY = 6
   gabZ = 7
   gabX1 = 2
   gabY1 = 3
   gabZ1 = 4
   DO IBMAT=1,NBMAT
   ! First half-contraction of Gaussian AO with density-matrix
    !TMP(K,J,N) = GAORED(K,I,N)*BRED(I,J) - TMP(K,J,1) is already built
    DO N=2,4
       M=N+3
     CALL DGEMM("N","N",NBLEN,NRED,NRED,1D0,GAORED(:,:,M),NBLEN,&
          &     BRED(:,:,IBMAT),NRED,0.0d0,TMP(:,:,N),NBLEN)
    ENDDO
   ! Second half-contraction of Gaussian AOs
    J=1
    DO K = 1, NBLEN
     Y = COORD(2,K) !- Origin(2)      
     Z = COORD(3,K) !- Origin(3)      
     EXPGRAD(K,1,2,IBMAT) = &
          & +GAORED(K,J,1)*(Z*TMP(K,J,gabY1)-Y*TMP(K,J,gabZ1))&
          & +TMP(K,J,1)*(Y*GAORED(K,J,gabZ)-Z*GAORED(K,J,gabY))
    END DO
    DO J = 2, NRED
     DO K = 1, NBLEN
        Y = COORD(2,K) !- Origin(2)      
        Z = COORD(3,K) !- Origin(3)      
        EXPGRAD(K,1,2,IBMAT) = EXPGRAD(K,1,2,IBMAT)&
             & +GAORED(K,J,1)*(Z*TMP(K,J,gabY1)-Y*TMP(K,J,gabZ1))&
             & +TMP(K,J,1)*(Y*GAORED(K,J,gabZ)-Z*GAORED(K,J,gabY))
     END DO
    ENDDO
    J=1
    DO K = 1, NBLEN
     X = COORD(1,K) !- Origin(1)      
     Z = COORD(3,K) !- Origin(3)      
     EXPGRAD(K,1,3,IBMAT) = &
          & +GAORED(K,J,1)*(X*TMP(K,J,gabZ1)-Z*TMP(K,J,gabX1))&
          & +TMP(K,J,1)*(Z*GAORED(K,J,gabX)-X*GAORED(K,J,gabZ))
    END DO
    DO J = 2, NRED
     DO K = 1, NBLEN
        X = COORD(1,K) !- Origin(1)      
        Z = COORD(3,K) !- Origin(3)      
        EXPGRAD(K,1,3,IBMAT) = EXPGRAD(K,1,3,IBMAT)&
             & +GAORED(K,J,1)*(X*TMP(K,J,gabZ1)-Z*TMP(K,J,gabX1))&
             & +TMP(K,J,1)*(Z*GAORED(K,J,gabX)-X*GAORED(K,J,gabZ))
     END DO
    ENDDO
    J=1
    DO K = 1, NBLEN
     X = COORD(1,K) !- Origin(1)      
     Y = COORD(2,K) !- Origin(2)      
     EXPGRAD(K,1,4,IBMAT) = &
          & +GAORED(K,J,1)*(Y*TMP(K,J,gabX1)-X*TMP(K,J,gabY1))&
          & +TMP(K,J,1)*(X*GAORED(K,J,gabY)-Y*GAORED(K,J,gabX))
    END DO
    DO J = 2, NRED
     DO K = 1, NBLEN
        X = COORD(1,K) !- Origin(1)      
        Y = COORD(2,K) !- Origin(2)      
        EXPGRAD(K,1,4,IBMAT) = EXPGRAD(K,1,4,IBMAT)&
             & +GAORED(K,J,1)*(Y*TMP(K,J,gabX1)-X*TMP(K,J,gabY1))&
             & +TMP(K,J,1)*(X*GAORED(K,J,gabY)-Y*GAORED(K,J,gabX))
     END DO
    END DO
   ENDDO
 ELSE
   DO IBMAT=1,NBMAT
    DO N=1,3
     DO K = 1, NBLEN
      EXPGRAD(K,1,N,IBMAT)=zero
     ENDDO
    ENDDO
   ENDDO
 ENDIF
!-------------------------------------------------------------
!  The geometrical derivative of the expectation value
!-------------------------------------------------------------
! I should make these as subroutines also called from DISTGGA
! WARNING THIS COULD BE IMPROVED MAYBE SOME DGEMM
 DO IBMAT=1,NBMAT
    DO J = 1, NRED
     BREDIJ = BRED(1,J,IBMAT)
     BREDJI = BRED(J,1,IBMAT)
     DO K = 1, NBLEN
      TMP(K,J,1) = GAORED(K,1,1)*(BREDIJ + BREDJI)
     END DO
     DO I = 2, NRED
      BREDIJ = BRED(I,J,IBMAT)
      BREDJI = BRED(J,I,IBMAT)
      DO K = 1, NBLEN
       TMP(K,J,1) = TMP(K,J,1) + GAORED(K,I,1)*(BREDIJ + BREDJI)
      END DO
     ENDDO
    ENDDO
      ! Second half-contraction
    DO K = 1, NBLEN
      EXPGRAD(K,2,1,IBMAT) = GAORED(K,1,2)*TMP(K,1,1)
      EXPGRAD(K,3,1,IBMAT) = GAORED(K,1,3)*TMP(K,1,1)
      EXPGRAD(K,4,1,IBMAT) = GAORED(K,1,4)*TMP(K,1,1)
    ENDDO
    DO I = 2, NRED
      DO K = 1, NBLEN
         EXPGRAD(K,2,1,IBMAT) = EXPGRAD(K,2,1,IBMAT) + GAORED(K,I,2)*TMP(K,I,1)
         EXPGRAD(K,3,1,IBMAT) = EXPGRAD(K,3,1,IBMAT) + GAORED(K,I,3)*TMP(K,I,1)
         EXPGRAD(K,4,1,IBMAT) = EXPGRAD(K,4,1,IBMAT) + GAORED(K,I,4)*TMP(K,I,1)
      ENDDO
    ENDDO
 ENDDO
!-------------------------------------------------------------
!  The mixed geometrical and magnetic derivative of the expectation value
!-------------------------------------------------------------
 do IBMAT=1,NBMAT
  do I=2,4
   do J=2,4
    do K=1,NBLEN
      EXPGRAD(K,J,I,IBMAT) = 0.d0
    enddo
   enddo
  enddo
 enddo
 DO IBMAT=1,NBMAT
  DO J = 1, NRED
   DO I = 1, NRED
    DO ALPHA =1,3
     beta = betalist(ALPHA)
     gamma = gammalist(ALPHA)
     gab1beta = 4+beta
     gab1gamma = 4+gamma
     gab2Xbeta = 5+3*beta
     gab2Xgamma = 5+3*gamma
     gab2Ybeta = 6+3*beta
     gab2Ygamma = 6+3*gamma
     gab2Zbeta = 7+3*beta
     gab2Zgamma = 7+3*gamma
     DO K = 1, NBLEN
        Rgamma = COORD(gamma,K) !- Origin(gamma)       
        Rbeta = COORD(beta,K) !- Origin(beta)
        EXPGRAD(K,gamma+1,alpha+1,IBMAT) = EXPGRAD(K,gamma+1,alpha+1,IBMAT)&
             &+GAORED(K,I,gab1beta)*BRED(I,J,IBMAT)*GAORED(K,J,1)&
             &-GAORED(K,I,1)*BRED(I,J,IBMAT)*GAORED(K,J,gab1beta)

        EXPGRAD(K,beta+1,alpha+1,IBMAT) = EXPGRAD(K,beta+1,alpha+1,IBMAT)&
             &-GAORED(K,I,gab1gamma)*BRED(I,J,IBMAT)*GAORED(K,J,1)&
             &+GAORED(K,I,1)*BRED(I,J,IBMAT)*GAORED(K,J,gab1gamma)

        EXPGRAD(K,2,alpha+1,IBMAT) = EXPGRAD(K,2,alpha+1,IBMAT)&
             &+GAORED(K,I,gab2Xbeta)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,1)&
             &-GAORED(K,I,gab2Xgamma)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,1)&
             &-GAORED(K,I,1)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,gab2Xbeta)&
             &+GAORED(K,I,1)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,gab2Xgamma)&
             &+GAORED(K,I,gab1beta)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,2)&
             &-GAORED(K,I,gab1gamma)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,2)&
             &-GAORED(K,I,2)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,gab1beta)&
             &+GAORED(K,I,2)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,gab1gamma)
        
        EXPGRAD(K,3,alpha+1,IBMAT) = EXPGRAD(K,3,alpha+1,IBMAT)&
             &+GAORED(K,I,gab2Ybeta)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,1)&
             &-GAORED(K,I,gab2Ygamma)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,1)&
             &-GAORED(K,I,1)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,gab2Ybeta)&
             &+GAORED(K,I,1)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,gab2Ygamma)&
             &+GAORED(K,I,gab1beta)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,3)&
             &-GAORED(K,I,gab1gamma)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,3)&
             &-GAORED(K,I,3)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,gab1beta)&
             &+GAORED(K,I,3)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,gab1gamma)
        
        EXPGRAD(K,4,alpha+1,IBMAT) = EXPGRAD(K,4,alpha+1,IBMAT)&
             &+GAORED(K,I,gab2Zbeta)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,1)&
             &-GAORED(K,I,gab2Zgamma)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,1)&
             &-GAORED(K,I,1)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,gab2Zbeta)&
             &+GAORED(K,I,1)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,gab2Zgamma)&
             &+GAORED(K,I,gab1beta)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,4)&
             &-GAORED(K,I,gab1gamma)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,4)&
             &-GAORED(K,I,4)*Rgamma*BRED(I,J,IBMAT)*GAORED(K,J,gab1beta)&
             &+GAORED(K,I,4)*Rbeta*BRED(I,J,IBMAT)*GAORED(K,J,gab1gamma)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
 DEALLOCATE(BRED)
ENDIF

END SUBROUTINE II_GET_MAGDERIV_EXPECTATIONVAL_GGA

!> \brief main geometrical derivative of the kohn-sham matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_kohnsham_mat(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!
INTEGER          :: I,J,KMAX,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL          :: DOGGA,DOLONDON
INTEGER          :: GRDONE,NHTYP,IDMAT,NGEODERIV

  DOGGA = DFT_ISGGA()
  DOLONDON = .FALSE.
  NGEODERIV = 1
  DFTDATA%ENERGY = 0D0
  IF(DOGGA) THEN
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL LSQUIT('II_DFT_geoderiv_kohnshamGGAUNRES not implemented',lupri)
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
!             & II_DFT_geoderiv_kohnshamGGAUNRES,DFTDATA)
     ELSE
!        CALL LSQUIT('II_DFT_geoderiv_kohnshamGGA not implemented',lupri)
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
             & II_DFT_geoderiv_kohnshamGGA,DFTDATA)
     ENDIF
  ELSE 
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL LSQUIT('II_DFT_geoderiv_kohnshamLDAUNRES not implemented',lupri)
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
!             & II_DFT_geoderiv_kohnshamLDAUNRES,DFTDATA)
     ELSE !DEFAULT (RESTRICTED)
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
             & II_DFT_geoderiv_kohnshamLDA,DFTDATA)
     END IF
  ENDIF
  NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
  IF(IPRINT.GE.0) WRITE(LUPRI,'(A,2F20.14,A,E9.2)')&
       &     'KS electrons/energy:', DFTDATA%ELECTRONS, DFTDATA%ENERGY,&
       &     ' rel.err:', (DFTDATA%ELECTRONS-NELE)/(NELE)

END SUBROUTINE II_DFT_GEODERIV_KOHNSHAM_MAT

!> \brief geometrical derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_kohnshamLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,MXBLLEN,COORD,WGHT,&
     & ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK) :: VXC(NBLEN),VX(14),DFTENE,EXPVAL(NBLEN,1)
REAL(REALK) :: fR(NBLEN),fRR(NBLEN)
REAL(REALK),parameter :: D2=2.d0,D4=4.d0,DUMMY = 0D0 
INTEGER     :: I,J,nred,nbmat,ipnt
logical :: DOCALC
nbmat=1
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
 call II_get_expectationval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
!      &Nactbast,NBAST,GAO,EXPVAL,DFTDATA%BMAT,DFTDATA%nBMAT,DFTHRI,NRED)
      &Nactbast,NBAST,GAO,EXPVAL,DFTDATA%BMAT,1,DFTHRI,NRED)
 IF(NRED.GT.0)THEN
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      CALL dft_funcderiv2(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
      fR(IPNT) = VX(1)*D4 
      fRR(IPNT) = D2*(VX(6) + VX(7))*EXPVAL(IPNT,1) 
   ELSE
      VXC(IPNT) = 0.0d0
   ENDIF
  END DO
   CALL II_DFT_distgeoderiv_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,&
        &NBAST,1,NBLEN,fR,fRR,GAO(:,:,:),DFTDATA%GRAD,DFTDATA%orb2atom,&
        &DFTDATA%BMAT,nbmat,DMAT,ndmat,DFTDATA%natoms,COORD,DFTHRI,NTYPSO)
 ENDIF
ENDIF
END SUBROUTINE II_DFT_GEODERIV_KOHNSHAMLDA

!> \brief distribution routine for geometrical derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTGEODERIV_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,IBLSTART,IBLEND,COEF1,COEF2,GAOS,GRAD,orb2atom,&
     &BMAT,nbmat,DMAT,ndmat,natoms,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in) :: IBLSTART
!> end gridpoint
INTEGER,intent(in) :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> molecular gradient
REAL(REALK),intent(inout) :: GRAD(3,NATOMS)
!> for given orbital, the atom it is attached to
INTEGER,intent(in) :: orb2atom(nbast)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> The density matrix
REAL(REALK),intent(in) :: DMAT(NactBAST,NactBAST,NDMAT)
!> Number of D matrices
INTEGER,intent(in)  :: NDMAT
!> Number of atoms
INTEGER,intent(in)  :: NATOMS
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
INTEGER     :: atom(Nactbast)
REAL(REALK) :: TMPB(NBLEN,NACTBAST),TMPD(NBLEN,NACTBAST),GAOMAX
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST),TMP2,BT
REAL(REALK) :: GAORED(NBLEN,NACTBAST,4)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,Jact,Iact,iatom
INTEGER     :: COORDINATE,IORB,JORB
real(realk),parameter :: D2=2.D0,D05=0.5D0
REAL(REALK),allocatable :: DRED(:,:),BRED(:,:)
NK=IBLEND-IBLSTART+1
GAOMAX = 0.0d0
!        Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
  DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
   GAOGMX(J) = 0.0d0
   DO K = IBLSTART, IBLEND
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOMAX = MAX(GAOMAX,GAOGMX(J))
 ENDDO
ENDDO
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = J 
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT.0) THEN
! Set up reduced Gaussian AO's
   ALLOCATE(DRED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JACT = INXRED(JRED)     !J is active index
      DO IRED=1,NRED          !IRED is reduced index
         IACT = INXRED(IRED)  !I is active index
         DRED(IRED,JRED) = DMAT(IACT,JACT,1)
      ENDDO
   ENDDO

   NRED = 0 
   DO JBL = 1, NBLOCKS
    DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
     IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
      NRED = NRED + 1
      !WARNING WE CHANGE THE USE OF INXRED
      INXRED(NRED) = INXACT(J)
      !inxred for a given redindex the corresponding orbitalindex
      DO N=1,NTYPSO
       DO K = IBLSTART, IBLEND
        GAORED(K,NRED,N) = GAOS(K,J,N)
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DO IRED=1,NRED
      I = INXRED(IRED) !oribtal index
      atom(IRED) = orb2atom(I)
   ENDDO

   ALLOCATE(BRED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED) = BMAT(IORB,JORB,1)
      ENDDO
   ENDDO

   CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0d0,GAORED,&
        &          NBLEN,DRED,NRED,0.0d0,TMPD,NBLEN)
!   CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0d0,GAORED,&
!        &          NBLEN,BRED,NRED,0.0d0,TMPB,NBLEN)
   DO JRED = 1,NRED
      IRED =1
      BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
      DO K=IBLSTART, IBLEND
         TMPB(K,JRED) = GAORED(K,IRED,1)*BT
      ENDDO
      DO IRED =2,NRED
         BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
         DO K=IBLSTART, IBLEND
            TMPB(K,JRED) = TMPB(K,JRED)+GAORED(K,IRED,1)*BT
         ENDDO
      ENDDO
   ENDDO
!$OMP CRITICAL
 DO IRED = 1,NRED
   iatom = atom(IRED)
   DO COORDINATE = 1,3
      DO K=IBLSTART, IBLEND
         GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
              & - coef1(K)*GAORED(K,IRED,1+COORDINATE)*TMPB(K,IRED)&
              & - coef2(K)*GAORED(K,IRED,1+COORDINATE)*TMPD(K,IRED)
      ENDDO
   ENDDO
 ENDDO
!$OMP END CRITICAL
 DEALLOCATE(DRED)
 DEALLOCATE(BRED)
ENDIF

END SUBROUTINE II_DFT_DISTGEODERIV_LDA

!> \brief geometrical derivative Kohn-sham matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_kohnshamGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,MXBLLEN,COORD,WGHT,&
     & ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK) :: VX(14),DFTENE,EXPVAL(NBLEN,1)
REAL(REALK) :: EXPGRAD(3,NBLEN,1),GRD,GRDA,MIXEDGRDA
REAL(REALK) :: fR,fRR,fZ,fRZ,fZZ,fRG,fZG,fGG,fG,A,B
REAL(REALK) :: VXC1(NBLEN),VXC2(3,NBLEN),VXC3(NBLEN),VXC4(3,NBLEN)
REAL(REALK),parameter :: D2=2.d0, D4=4.d0, DUMMY = 0D0, D05=0.5d0, D025=0.25d0 
REAL(REALK),parameter :: D8=8.d0
INTEGER     :: I,J,nred,nbmat,IPNT
logical :: DOCALC
nbmat=1
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
 call II_get_expectationval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
!      &Nactbast,NBAST,GAO,EXPVAL,DFTDATA%BMAT,DFTDATA%nBMAT,DFTHRI,NRED)
      &Nactbast,NBAST,GAO,EXPVAL,EXPGRAD,DFTDATA%BMAT,1,DFTHRI,NRED)
 IF(NRED.GT.0)THEN
  DO IPNT = 1, NBLEN
     GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
          &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
     GRDA = D05*GRD
     IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
       CALL dft_funcderiv2(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
       fR  = D05*(VX(1) + VX(2))   !0.5*(drvs.df1000 + drvs.df0100);
       fZ  = VX(3)                    !drvs.df0010;
       fRR = D05*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
       fRZ = D05*(VX(8) + VX(9))   !0.5*(drvs.df1010 + drvs.df1001);
       fZZ = D05*(VX(11) + VX(12)) !0.5*(drvs.df0020 + drvs.df0011);
       fRG = D05*VX(10)             !0.5*drvs.df10001;   
       fZG = D05*VX(13)             !0.5*drvs.df00101; 
       fGG = D025*VX(14)            !0.25*drvs.df00002; 
       fG  = D05*VX(5)               !0.5*drvs.df00001;  
       MIXEDGRDA = (EXPGRAD(1,IPNT,1)*GRAD(1,IPNT,1)&
            &+EXPGRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
            &+EXPGRAD(3,IPNT,1)*GRAD(3,IPNT,1))

       VXC1(IPNT) = D4*VX(1)
       A = D2*(VX(3)/GRDA + VX(5))
       VXC2(1,IPNT) = A*GRAD(1,IPNT,1)
       VXC2(2,IPNT) = A*GRAD(2,IPNT,1)
       VXC2(3,IPNT) = A*GRAD(3,IPNT,1)
       !the LDA part
       VXC3(IPNT) =D4*fRR*EXPVAL(IPNT,1)+D4*(fRZ/GRD+fRG)*MIXEDGRDA
       !the non LDA parts
       A = D4*((fRZ/GRD + fRG)*EXPVAL(IPNT,1)&
            & + (((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)*MIXEDGRDA)
       B= D4*(fZ/GRD + fG)
       VXC4(1,IPNT) = A*GRAD(1,IPNT,1)+B*EXPGRAD(1,IPNT,1)
       VXC4(2,IPNT) = A*GRAD(2,IPNT,1)+B*EXPGRAD(2,IPNT,1)
       VXC4(3,IPNT) = A*GRAD(3,IPNT,1)+B*EXPGRAD(3,IPNT,1)
    ELSE
       VXC1(IPNT) = 0.0d0
       VXC2(:,IPNT) = 0.0d0
       VXC3(IPNT) = 0.0d0
       VXC4(:,IPNT) = 0.0d0
    ENDIF
  END DO
   CALL II_DFT_distgeoderiv_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,&
        &NBAST,1,NBLEN,VXC1,VXC2,VXC3,VXC4,GAO(:,:,:),DFTDATA%GRAD,&
        &DFTDATA%orb2atom,DFTDATA%BMAT,nbmat,DMAT,ndmat,DFTDATA%natoms,COORD,DFTHRI,NTYPSO)
 ENDIF
ENDIF
END SUBROUTINE II_DFT_GEODERIV_KOHNSHAMGGA

!> \brief distribution routine for geometrical derivative Kohn-sham matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTGEODERIV_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,IBLSTART,IBLEND,COEF1,COEF2,COEF3,COEF4,GAOS,GRAD,orb2atom,&
     &BMAT,nbmat,DMAT,ndmat,natoms,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in) :: IBLSTART
!> end gridpoint
INTEGER,intent(in) :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(3,NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF3(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF4(3,NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> molecular gradient
REAL(REALK),intent(inout) :: GRAD(3,NATOMS)
!> for given orbital, the atom it is attached to
INTEGER,intent(in) :: orb2atom(NBAST)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> The density matrix
REAL(REALK),intent(in) :: DMAT(NactBAST,NactBAST,NDMAT)
!> Number of D matrices
INTEGER,intent(in)  :: NDMAT
!> Number of atoms
INTEGER,intent(in)  :: NATOMS
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
INTEGER     :: atom(Nactbast)
REAL(REALK) :: GAOMAX
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST),TMP2,BT,DT
REAL(REALK) :: GAORED(NBLEN,NACTBAST,NTYPSO)
REAL(REALK) :: TMPB(NBLEN,NactBast),TMPBX(NBLEN,NactBast)
REAL(REALK) :: TMPBY(NBLEN,NactBast),TMPBZ(NBLEN,NactBast)
REAL(REALK) :: TMPD(NBLEN,NactBast),TMPDX(NBLEN,NactBast)
REAL(REALK) :: TMPDY(NBLEN,NactBast),TMPDZ(NBLEN,NactBast)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,Jact,Iact,iatom
INTEGER     :: COORDINATE,IORB,JORB,R,RX,RY,RZ
integer,parameter     :: xlist(3)=(/5,6,7/),ylist(3)=(/6,8,9/)
integer,parameter     :: zlist(3)=(/7,9,10/)
real(realk),parameter :: D2=2.D0,D05=0.5D0,D8=8.D0
REAL(REALK),allocatable :: DRED(:,:),BRED(:,:)
NK=IBLEND-IBLSTART+1
GAOMAX = 0.0d0
!        Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
  DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
   GAOGMX(J) = 0.0d0
   DO K = IBLSTART, IBLEND
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOMAX = MAX(GAOMAX,GAOGMX(J))
 ENDDO
ENDDO
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = J 
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT.0) THEN
! Set up reduced Gaussian AO's
   ALLOCATE(DRED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JACT = INXRED(JRED)     !J is active index
      DO IRED=1,NRED          !IRED is reduced index
         IACT = INXRED(IRED)  !I is active index
         DRED(IRED,JRED) = DMAT(IACT,JACT,1)
      ENDDO
   ENDDO

   NRED = 0 
   DO JBL = 1, NBLOCKS
    DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
     IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
      NRED = NRED + 1
      !WARNING WE CHANGE THE USE OF INXRED
      INXRED(NRED) = INXACT(J)
      !inxred for a given redindex the corresponding orbitalindex
      DO N=1,NTYPSO
       DO K = IBLSTART, IBLEND
        GAORED(K,NRED,N) = GAOS(K,J,N)
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DO IRED=1,NRED
      I = INXRED(IRED) !oribtal index
      atom(IRED) = orb2atom(I)
   ENDDO

   ALLOCATE(BRED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED) = BMAT(IORB,JORB,1)
      ENDDO
   ENDDO

   DO JRED = 1,NRED
      IRED =1
      BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
      DT = D05*(DRED(IRED,JRED)+DRED(JRED,IRED))
      DO K=IBLSTART, IBLEND
         TMPB(K,JRED) = GAORED(K,IRED,1)*BT
         TMPBX(K,JRED) = GAORED(K,IRED,2)*BT
         TMPBY(K,JRED) = GAORED(K,IRED,3)*BT
         TMPBZ(K,JRED) = GAORED(K,IRED,4)*BT
         TMPD(K,JRED) = GAORED(K,IRED,1)*DT
         TMPDX(K,JRED) = GAORED(K,IRED,2)*DT
         TMPDY(K,JRED) = GAORED(K,IRED,3)*DT
         TMPDZ(K,JRED) = GAORED(K,IRED,4)*DT
      ENDDO
      DO IRED =2,NRED
         !SYMMETRIZE THESE TWO MATRICES BEFORE START SO THAT WE CAN USE DGEMM
         BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
         DT = D05*(DRED(IRED,JRED)+DRED(JRED,IRED))
         DO K=IBLSTART, IBLEND
            TMPB(K,JRED) = TMPB(K,JRED)+GAORED(K,IRED,1)*BT
            TMPBX(K,JRED) = TMPBX(K,JRED)+GAORED(K,IRED,2)*BT
            TMPBY(K,JRED) = TMPBY(K,JRED)+GAORED(K,IRED,3)*BT
            TMPBZ(K,JRED) = TMPBZ(K,JRED)+GAORED(K,IRED,4)*BT
            TMPD(K,JRED) = TMPD(K,JRED)+GAORED(K,IRED,1)*DT
            TMPDX(K,JRED) = TMPDX(K,JRED)+GAORED(K,IRED,2)*DT
            TMPDY(K,JRED) = TMPDY(K,JRED)+GAORED(K,IRED,3)*DT
            TMPDZ(K,JRED) = TMPDZ(K,JRED)+GAORED(K,IRED,4)*DT
         ENDDO
      ENDDO
   ENDDO

!$OMP CRITICAL
   DO COORDINATE = 1,3
    R=1+COORDINATE
    RX=xlist(COORDINATE) 
    RY=ylist(COORDINATE) 
    RZ=zlist(COORDINATE) 
    DO IRED = 1,NRED
     iatom = atom(IRED)
!QQQQ
     DO K=IBLSTART, IBLEND
      GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
        &-coef1(K)*GAORED(K,IRED,R)*TMPB(K,IRED)
      GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
        &-coef2(1,K)*(GAORED(K,IRED,RX)*TMPB(K,IRED)+GAORED(K,IRED,R)*TMPBX(K,IRED))&
        &-coef2(2,K)*(GAORED(K,IRED,RY)*TMPB(K,IRED)+GAORED(K,IRED,R)*TMPBY(K,IRED))&
        &-coef2(3,K)*(GAORED(K,IRED,RZ)*TMPB(K,IRED)+GAORED(K,IRED,R)*TMPBZ(K,IRED))
      GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
        &-coef3(K)*GAORED(K,IRED,R)*TMPD(K,IRED)
      GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
        &-coef4(1,K)*(GAORED(K,IRED,R)*TMPDX(K,IRED)+GAORED(K,IRED,RX)*TMPD(K,IRED))&
        &-coef4(2,K)*(GAORED(K,IRED,R)*TMPDY(K,IRED)+GAORED(K,IRED,RY)*TMPD(K,IRED))&
        &-coef4(3,K)*(GAORED(K,IRED,R)*TMPDZ(K,IRED)+GAORED(K,IRED,RZ)*TMPD(K,IRED))
     ENDDO
    ENDDO
   ENDDO
!$OMP END CRITICAL
 DEALLOCATE(DRED)
 DEALLOCATE(BRED)
ENDIF

END SUBROUTINE II_DFT_DISTGEODERIV_GGA

!> \brief main geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_linrspgrad(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!
INTEGER          :: I,J,KMAX,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL          :: DOGGA,DOLONDON
INTEGER          :: GRDONE,NHTYP,IDMAT,NGEODERIV

  DOGGA = DFT_ISGGA()
  DOLONDON = .FALSE.
  NGEODERIV = 1
  DFTDATA%ENERGY = 0D0
  IF(DOGGA) THEN
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL LSQUIT('II_DFT_geoderiv_linrspGGAUNRES not implemented',lupri)
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
!             & II_DFT_geoderiv_linrspGGAUNRES,DFTDATA)
     ELSE
!        CALL LSQUIT('II_DFT_geoderiv_linrspGGA not implemented',lupri)
 CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
             & II_DFT_geoderiv_linrspGGA,DFTDATA)
     ENDIF
  ELSE 
     IF(NDMAT.EQ.2)THEN !UNRESTRICTED
        CALL LSQUIT('II_DFT_geoderiv_linrspLDAUNRES not implemented',lupri)
!        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
!             & II_DFT_geoderiv_linrspLDAUNRES,DFTDATA)
     ELSE !DEFAULT (RESTRICTED)
        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
             & II_DFT_geoderiv_linrspLDA,DFTDATA)
     END IF
  ENDIF
  NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
  IF(IPRINT.GE.0) WRITE(LUPRI,'(A,2F20.14,A,E9.2)')&
       &     'KS electrons/energy:', DFTDATA%ELECTRONS, DFTDATA%ENERGY,&
       &     ' rel.err:', (DFTDATA%ELECTRONS-NELE)/(NELE)

END SUBROUTINE II_DFT_GEODERIV_LINRSPGRAD

!> \brief LDA geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_linrspLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,MXBLLEN,COORD,WGHT,&
     & ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK) :: VXC(NBLEN),VX(27),DFTENE,EXPVAL(NBLEN,2)
REAL(REALK) :: fRRR(NBLEN),fRR(2,NBLEN),A
REAL(REALK),parameter :: D2=2.d0,D4=4.d0,D3=3.d0,DUMMY = 0D0 
INTEGER     :: I,J,nred,nbmat,IPNT
logical :: DOCALC
nbmat=DFTDATA%nBMAT
IF(nbmat.NE.2)call LSQUIT('II_DFT_geoderiv_linrspLDA requires 2 matrices',lupri)
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
 call II_get_expectationval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
      &Nactbast,NBAST,GAO,EXPVAL,DFTDATA%BMAT,nbmat,DFTHRI,NRED)
 IF(NRED.GT.0)THEN
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      CALL dft_funcderiv3(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
      A = D4*(VX(6) + VX(7))
      fRR(1,IPNT) = A*EXPVAL(IPNT,1) 
      fRR(2,IPNT) = A*EXPVAL(IPNT,2) 
      A = D2*(VX(15) + D3*VX(16))
      fRRR(IPNT) = A*EXPVAL(IPNT,1)*EXPVAL(IPNT,2)
   ELSE
      VXC(IPNT) = 0.0d0
   ENDIF
  END DO
   CALL II_DFT_distgeoderiv2_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,&
        &NBAST,1,NBLEN,fRR,fRRR,GAO(:,:,:),DFTDATA%GRAD,DFTDATA%orb2atom,&
        &DFTDATA%BMAT,nbmat,DMAT,ndmat,DFTDATA%natoms,COORD,DFTHRI,NTYPSO)
 ENDIF
ENDIF
END SUBROUTINE II_DFT_GEODERIV_LINRSPLDA

!> \brief distribution routine LDA geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTGEODERIV2_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,IBLSTART,IBLEND,COEF1,COEF2,GAOS,GRAD,orb2atom,&
     &BMAT,nbmat,DMAT,ndmat,natoms,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in) :: IBLSTART
!> end gridpoint
INTEGER,intent(in) :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(2,NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> molecular gradient
REAL(REALK),intent(inout) :: GRAD(3,NATOMS)
!> for given orbital, the atom it is attached to
INTEGER,intent(in) :: orb2atom(NBAST)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> The density matrix
REAL(REALK),intent(in) :: DMAT(NactBAST,NactBAST,NDMAT)
!> Number of D matrices
INTEGER,intent(in)  :: NDMAT
!> Number of atoms
INTEGER,intent(in)  :: NATOMS
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
INTEGER     :: atom(Nactbast)
REAL(REALK) :: TMPB(NBLEN,NACTBAST),TMPD(NBLEN,NACTBAST),GAOMAX
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: TMPA(NBLEN,NACTBAST),GAOGMX(NACTBAST),TMP2,BT,AT
REAL(REALK) :: GAORED(NBLEN,NACTBAST,4)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,Jact,Iact,iatom
INTEGER     :: COORDINATE,IORB,JORB
real(realk),parameter :: D2=2.D0,D05=0.5D0
REAL(REALK),allocatable :: DRED(:,:),BRED(:,:),ARED(:,:)
NK=IBLEND-IBLSTART+1
GAOMAX = 0.0d0
!        Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
  DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
   GAOGMX(J) = 0.0d0
   DO K = IBLSTART, IBLEND
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOMAX = MAX(GAOMAX,GAOGMX(J))
 ENDDO
ENDDO
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = J 
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT.0) THEN
! Set up reduced Gaussian AO's
   ALLOCATE(DRED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JACT = INXRED(JRED)     !J is active index
      DO IRED=1,NRED          !IRED is reduced index
         IACT = INXRED(IRED)  !I is active index
         DRED(IRED,JRED) = DMAT(IACT,JACT,1)
      ENDDO
   ENDDO

   NRED = 0 
   DO JBL = 1, NBLOCKS
    DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
     IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
      NRED = NRED + 1
      !WARNING WE CHANGE THE USE OF INXRED
      INXRED(NRED) = INXACT(J)
      !inxred for a given redindex the corresponding orbitalindex
      DO N=1,NTYPSO
       DO K = IBLSTART, IBLEND
        GAORED(K,NRED,N) = GAOS(K,J,N)
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DO IRED=1,NRED
      I = INXRED(IRED) !oribtal index
      atom(IRED) = orb2atom(I)
   ENDDO

   ALLOCATE(ARED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         ARED(IRED,JRED) = BMAT(IORB,JORB,1)
      ENDDO
   ENDDO
   ALLOCATE(BRED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED) = BMAT(IORB,JORB,2)
      ENDDO
   ENDDO

   CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0d0,GAORED,&
        &          NBLEN,DRED,NRED,0.0d0,TMPD,NBLEN)
   DO JRED = 1,NRED
      IRED =1
      AT = D05*(ARED(IRED,JRED)+ARED(JRED,IRED))
      BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
      DO K=IBLSTART, IBLEND
         TMPB(K,JRED) = GAORED(K,IRED,1)*BT
      ENDDO
      DO IRED =2,NRED
         BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
         DO K=IBLSTART, IBLEND
            TMPB(K,JRED) = TMPB(K,JRED)+GAORED(K,IRED,1)*BT
         ENDDO
      ENDDO
      IRED =1
      DO K=IBLSTART, IBLEND
         TMPA(K,JRED) = GAORED(K,IRED,1)*AT
      ENDDO
      DO IRED =2,NRED
         AT = D05*(ARED(IRED,JRED)+ARED(JRED,IRED))
         DO K=IBLSTART, IBLEND
            TMPA(K,JRED) = TMPA(K,JRED)+GAORED(K,IRED,1)*AT
         ENDDO
      ENDDO
   ENDDO
!$OMP CRITICAL
 DO IRED = 1,NRED
   iatom = atom(IRED)
   DO COORDINATE = 1,3
      DO K=IBLSTART, IBLEND
         GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
              & - coef1(1,K)*GAORED(K,IRED,1+COORDINATE)*TMPB(K,IRED)&
              & - coef1(2,K)*GAORED(K,IRED,1+COORDINATE)*TMPA(K,IRED)&
              & - coef2(K)*GAORED(K,IRED,1+COORDINATE)*TMPD(K,IRED)
      ENDDO
   ENDDO
 ENDDO
!$OMP END CRITICAL
 DEALLOCATE(DRED)
 DEALLOCATE(BRED)
 DEALLOCATE(ARED)
ENDIF

END SUBROUTINE II_DFT_DISTGEODERIV2_LDA

!> \brief the GGA geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_linrspGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,MXBLLEN,COORD,WGHT,&
     & ENERGY,DFTDATA,RHOTHR,DFTHRI)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER,intent(in)     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK) :: VX(27),DFTENE,EXPVAL(NBLEN,2)
REAL(REALK) :: EXPGRAD(3,NBLEN,2)
REAL(REALK) :: VXC1(2,NBLEN),VXC2(3,NBLEN,2),VXC3(4,NBLEN)
REAL(REALK) :: ARHOGRAD,BRHOGRAD,ABRHOGRAD,GRDA,GRD2,GRD,GRDA2,GRDA3
REAL(REALK),parameter :: D2=2.d0,D4=4.d0,DUMMY = 0D0,D05=0.5D0,D025=0.25D0 
REAL(REALK),parameter :: D8=8.d0,D3=3.d0,D16=16.D0
REAL(REALK) :: fR,fZ,fRZ,fZZ,fRG,fZG,fGG,fG,B,facW,factorRZ,A,fRR,fRRR,fRRZ,fRRG
REAL(REALK) :: fRRGX,fRZZ,fZZZ,gradA,gradB,gradAB,C
INTEGER     :: I,J,nred,nbmat,IPNT
logical :: DOCALC
nbmat=DFTDATA%nBMAT
IF(nbmat.NE.2)call LSQUIT('II_DFT_geoderiv_linrspGGA requires 2 matrices',lupri)
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
   call II_get_expectationval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
        &Nactbast,NBAST,GAO,EXPVAL,EXPGRAD,DFTDATA%BMAT,DFTDATA%nBMAT,&
        &DFTHRI,NRED)
 IF(NRED.GT.0)THEN
  DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   GRDA = D05*GRD
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
    ARHOGRAD = (EXPGRAD(1,IPNT,1)*GRAD(1,IPNT,1)&
         &+EXPGRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
         &+EXPGRAD(3,IPNT,1)*GRAD(3,IPNT,1))
    BRHOGRAD = (EXPGRAD(1,IPNT,2)*GRAD(1,IPNT,1)&
         &+EXPGRAD(2,IPNT,2)*GRAD(2,IPNT,1)&
         &+EXPGRAD(3,IPNT,2)*GRAD(3,IPNT,1))
    ABRHOGRAD = (EXPGRAD(1,IPNT,1)*EXPGRAD(1,IPNT,2)&
         &+EXPGRAD(2,IPNT,1)*EXPGRAD(2,IPNT,2)&
         &+EXPGRAD(3,IPNT,1)*EXPGRAD(3,IPNT,2))     
    GRD2 = GRD*GRD
    GRDA2 = GRDA*GRDA
    GRDA3 = GRDA2*GRDA
    CALL dft_funcderiv3(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
    !THE NON DIFFERENTIATED PART
    fR  = D05*(VX(1) + VX(2))    !0.5*(drvs.df1000 + drvs.df0100);
    fZ  = VX(3)                  !drvs.df0010;
    fRR = D05*(VX(6) + VX(7))    !0.5*(drvs.df2000 + drvs.df1100);
    fRZ = D05*(VX(8) + VX(9))    !0.5*(drvs.df1010 + drvs.df1001);
    fZZ = D05*(VX(11) + VX(12))  !0.5*(drvs.df0020 + drvs.df0011);
    fRG = D05*VX(10)             !0.5*drvs.df10001;   
    fZG = D05*VX(13)             !0.5*drvs.df00101; 
    fGG = D025*VX(14)            !0.25*drvs.df00002; 
    fG  = D05*VX(5)              !0.5*drvs.df00001;  
    VXC1(1,IPNT) = D8*fRR*EXPVAL(IPNT,1) + D8*(fRZ/GRD+fRG)*ARHOGRAD
    VXC1(2,IPNT) = D8*fRR*EXPVAL(IPNT,2) + D8*(fRZ/GRD+fRG)*BRHOGRAD
    B= D8*(fZ/GRD + fG)
    FacW = D8*(((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)
    FactorRZ = D8*(fRZ/GRD + fRG)
    A = FactorRZ*EXPVAL(IPNT,1)+FacW*ARHOGRAD
    VXC2(1,IPNT,1) = B*EXPGRAD(1,IPNT,1)+A*GRAD(1,IPNT,1)
    VXC2(2,IPNT,1) = B*EXPGRAD(2,IPNT,1)+A*GRAD(2,IPNT,1)
    VXC2(3,IPNT,1) = B*EXPGRAD(3,IPNT,1)+A*GRAD(3,IPNT,1)
    A = FactorRZ*EXPVAL(IPNT,2)+FacW*BRHOGRAD
    VXC2(1,IPNT,2) = B*EXPGRAD(1,IPNT,2)+A*GRAD(1,IPNT,1)
    VXC2(2,IPNT,2) = B*EXPGRAD(2,IPNT,2)+A*GRAD(2,IPNT,1)
    VXC2(3,IPNT,2) = B*EXPGRAD(3,IPNT,2)+A*GRAD(3,IPNT,1)
    !THE DIFFERENTIATED PART
    fRRR = (VX(15) + D3*VX(16))
    fRZ = (VX(8)+VX(9))/GRD  !warning change definition of fRZ
    fRG = D4*fRG                 !D2*VX(10)                      
    fZZ = (VX(11)+VX(12))/GRD2-VX(3)/(GRD2*GRDA)
    fRRZ = (VX(17)+VX(18)+D2*VX(20))/GRD 
    fRRG = VX(19)+VX(21)                 
    fRRGX = D2*(VX(19)+VX(21))           
    fRZZ = (VX(22)+VX(24)+D2*VX(23))/GRD2 - (VX(8)+VX(9))/(GRD2*GRDA)  
    fZZZ = ((VX(25)+D3*VX(26))/(GRDA3)-D3*(VX(11)+VX(12))/(GRDA2*GRDA2)+D3*VX(3)/(GRDA3*GRDA2))/D8

    VXC3(1,IPNT) = D2*fRRR*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) & 
         &+D2*fRRZ*(EXPVAL(IPNT,1)*BRHOGRAD+EXPVAL(IPNT,2)*ARHOGRAD) & 
         &+D2*BRHOGRAD*ARHOGRAD*fRZZ& 
         &+D4*ABRHOGRAD*fRZ &
         &+D2*fRRG*(EXPVAL(IPNT,1)*BRHOGRAD+EXPVAL(IPNT,2)*ARHOGRAD)&
         &+D2*fRG*ABRHOGRAD 
 
    A = D2*fZZZ*ARHOGRAD*BRHOGRAD &
         & + D2*(fRZZ*EXPVAL(IPNT,1)*BRHOGRAD + fRZZ*EXPVAL(IPNT,2)*ARHOGRAD) &
         & + D2*fRRZ*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) &
         & + fRRGX*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) + D4*fZZ*ABRHOGRAD
    B = D4*fZZ*ARHOGRAD + D4*fRZ*EXPVAL(IPNT,1) + D2*fRG*EXPVAL(IPNT,1)
    C = D4*fZZ*BRHOGRAD + D4*fRZ*EXPVAL(IPNT,2) + D2*fRG*EXPVAL(IPNT,2)
    
    VXC3(2,IPNT) = A*GRAD(1,IPNT,1) + B*EXPGRAD(1,IPNT,2) + C*EXPGRAD(1,IPNT,1)
    VXC3(3,IPNT) = A*GRAD(2,IPNT,1) + B*EXPGRAD(2,IPNT,2) + C*EXPGRAD(2,IPNT,1)
    VXC3(4,IPNT) = A*GRAD(3,IPNT,1) + B*EXPGRAD(3,IPNT,2) + C*EXPGRAD(3,IPNT,1)
   ELSE
      VXC1(:,IPNT) = 0.0d0
      VXC2(:,IPNT,:) = 0.0d0
      VXC3(:,IPNT) = 0.0d0
   ENDIF
  END DO
  CALL II_DFT_distgeoderiv2_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,&
       &NBAST,1,NBLEN,VXC1,VXC2,VXC3,GAO(:,:,:),DFTDATA%GRAD,DFTDATA%orb2atom,&
       &DFTDATA%BMAT,nbmat,DMAT,ndmat,DFTDATA%natoms,COORD,DFTHRI,NTYPSO)
 ENDIF
ENDIF
END SUBROUTINE II_DFT_GEODERIV_LINRSPGGA

!> \brief distribution routine GGA geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTGEODERIV2_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,IBLSTART,IBLEND,COEF1,COEF2,COEF3,GAOS,GRAD,&
     &orb2atom,BMAT,nbmat,DMAT,ndmat,natoms,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> start gridpoint
INTEGER,intent(in) :: IBLSTART
!> end gridpoint
INTEGER,intent(in) :: IBLEND
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(2,NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(3,NBLEN,2)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF3(4,NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> molecular gradient
REAL(REALK),intent(inout) :: GRAD(3,NATOMS)
!> for given orbital, the atom it is attached to
INTEGER,intent(in) :: orb2atom(NBAST)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> The density matrix
REAL(REALK),intent(in) :: DMAT(NactBAST,NactBAST,NDMAT)
!> Number of D matrices
INTEGER,intent(in)  :: NDMAT
!> Number of atoms
INTEGER,intent(in)  :: NATOMS
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
INTEGER     :: atom(Nactbast)
REAL(REALK) :: TMPB(NBLEN,NACTBAST),TMPD(NBLEN,NACTBAST,4),GAOMAX
INTEGER     :: NK,ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: TMPA(NBLEN,NACTBAST),GAOGMX(NACTBAST),TMP2,BT,AT
REAL(REALK) :: TMPAX(NBLEN,NACTBAST),TMPAY(NBLEN,NACTBAST),TMPAZ(NBLEN,NACTBAST)
REAL(REALK) :: TMPBX(NBLEN,NACTBAST),TMPBY(NBLEN,NACTBAST),TMPBZ(NBLEN,NACTBAST)
REAL(REALK) :: GAORED(NBLEN,NACTBAST,NTYPSO)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,Jact,Iact,iatom
INTEGER     :: COORDINATE,IORB,JORB,R,RX,RY,RZ
real(realk),parameter :: D2=2.D0,D05=0.5D0
REAL(REALK),allocatable :: DRED(:,:),BRED(:,:),ARED(:,:)
integer,parameter     :: xlist(3)=(/5,6,7/),ylist(3)=(/6,8,9/)
integer,parameter     :: zlist(3)=(/7,9,10/)


NK=IBLEND-IBLSTART+1
GAOMAX = 0.0d0
!        Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
  DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
   GAOGMX(J) = 0.0d0
   DO K = IBLSTART, IBLEND
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOMAX = MAX(GAOMAX,GAOGMX(J))
 ENDDO
ENDDO
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = J 
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT.0) THEN
! Set up reduced Gaussian AO's
   ALLOCATE(DRED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JACT = INXRED(JRED)     !J is active index
      DO IRED=1,NRED          !IRED is reduced index
         IACT = INXRED(IRED)  !I is active index
         DRED(IRED,JRED) = DMAT(IACT,JACT,1)
      ENDDO
   ENDDO

   NRED = 0 
   DO JBL = 1, NBLOCKS
    DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
     IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
      NRED = NRED + 1
      !WARNING WE CHANGE THE USE OF INXRED
      INXRED(NRED) = INXACT(J)
      !inxred for a given redindex the corresponding orbitalindex
      DO N=1,NTYPSO
       DO K = IBLSTART, IBLEND
        GAORED(K,NRED,N) = GAOS(K,J,N)
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DO IRED=1,NRED
      I = INXRED(IRED) !oribtal index
      atom(IRED) = orb2atom(I)
   ENDDO

   ALLOCATE(ARED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         ARED(IRED,JRED) = BMAT(IORB,JORB,1)
      ENDDO
   ENDDO
   ALLOCATE(BRED(NRED,NRED))
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED) = BMAT(IORB,JORB,2)
      ENDDO
   ENDDO
   DO I=1,4
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0d0,GAORED(:,:,I),&
           &          NBLEN,DRED,NRED,0.0d0,TMPD(:,:,I),NBLEN)
   ENDDO
   DO JRED = 1,NRED
      IRED =1
      AT = D05*(ARED(IRED,JRED)+ARED(JRED,IRED))
      BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
      DO K=IBLSTART, IBLEND
         TMPB(K,JRED) = GAORED(K,IRED,1)*BT
         TMPBX(K,JRED) = GAORED(K,IRED,2)*BT
         TMPBY(K,JRED) = GAORED(K,IRED,3)*BT
         TMPBZ(K,JRED) = GAORED(K,IRED,4)*BT
      ENDDO
      DO IRED =2,NRED
         BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
         DO K=IBLSTART, IBLEND
            TMPB(K,JRED) = TMPB(K,JRED)+GAORED(K,IRED,1)*BT
            TMPBX(K,JRED) = TMPBX(K,JRED)+GAORED(K,IRED,2)*BT
            TMPBY(K,JRED) = TMPBY(K,JRED)+GAORED(K,IRED,3)*BT
            TMPBZ(K,JRED) = TMPBZ(K,JRED)+GAORED(K,IRED,4)*BT
         ENDDO
      ENDDO
      IRED =1
      DO K=IBLSTART, IBLEND
         TMPA(K,JRED) = GAORED(K,IRED,1)*AT
         TMPAX(K,JRED) = GAORED(K,IRED,2)*AT
         TMPAY(K,JRED) = GAORED(K,IRED,3)*AT
         TMPAZ(K,JRED) = GAORED(K,IRED,4)*AT
      ENDDO
      DO IRED =2,NRED
         AT = D05*(ARED(IRED,JRED)+ARED(JRED,IRED))
         DO K=IBLSTART, IBLEND
            TMPA(K,JRED) = TMPA(K,JRED)+GAORED(K,IRED,1)*AT
            TMPAX(K,JRED) = TMPAX(K,JRED)+GAORED(K,IRED,2)*AT
            TMPAY(K,JRED) = TMPAY(K,JRED)+GAORED(K,IRED,3)*AT
            TMPAZ(K,JRED) = TMPAZ(K,JRED)+GAORED(K,IRED,4)*AT
         ENDDO
      ENDDO
   ENDDO
!$OMP CRITICAL
 DO IRED = 1,NRED
   iatom = atom(IRED)
   DO COORDINATE = 1,3
      R=1+COORDINATE
      RX=xlist(COORDINATE) 
      RY=ylist(COORDINATE) 
      RZ=zlist(COORDINATE) 
      DO K=IBLSTART, IBLEND
         GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
              & - coef1(1,K)*GAORED(K,IRED,R)*TMPB(K,IRED)&
              & - coef1(2,K)*GAORED(K,IRED,R)*TMPA(K,IRED)&
              & - coef2(1,K,1)*(TMPB(K,IRED)*GAORED(K,IRED,RX)+TMPBX(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(2,K,1)*(TMPB(K,IRED)*GAORED(K,IRED,RY)+TMPBY(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(3,K,1)*(TMPB(K,IRED)*GAORED(K,IRED,RZ)+TMPBZ(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(1,K,2)*(TMPA(K,IRED)*GAORED(K,IRED,RX)+TMPAX(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(2,K,2)*(TMPA(K,IRED)*GAORED(K,IRED,RY)+TMPAY(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(3,K,2)*(TMPA(K,IRED)*GAORED(K,IRED,RZ)+TMPAZ(K,IRED)*GAORED(K,IRED,R))&
              & - coef3(1,K)*GAORED(K,IRED,R)*TMPD(K,IRED,1)&
              & - coef3(2,K)*(TMPD(K,IRED,1)*GAORED(K,IRED,RX)+TMPD(K,IRED,2)*GAORED(K,IRED,R))&
              & - coef3(3,K)*(TMPD(K,IRED,1)*GAORED(K,IRED,RY)+TMPD(K,IRED,3)*GAORED(K,IRED,R))&
              & - coef3(4,K)*(TMPD(K,IRED,1)*GAORED(K,IRED,RZ)+TMPD(K,IRED,4)*GAORED(K,IRED,R))

      ENDDO
   ENDDO
 ENDDO
!$OMP END CRITICAL
 DEALLOCATE(DRED)
 DEALLOCATE(BRED)
 DEALLOCATE(ARED)
ENDIF

END SUBROUTINE II_DFT_DISTGEODERIV2_GGA

END MODULE IIDFTKSM
