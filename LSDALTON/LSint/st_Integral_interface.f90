!> @file
!> Contains obsolete set of subroutines for interfacing integral interface to Thermie driver
MODULE st_Integral_Intf
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use memory_handling
  use Integralinfo
  use integraldriver
  use BUILDAOBATCH
  use lstiming
  use ls_Integral_Interface

CONTAINS
!> \brief calculates the 3 center 2 electron repulsion integral matrix (ab|0alpha) for a given auxilliary basis function alpha
!> \author T. Kjaergaard
!> \date 2010
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
!> \param SETTING Integral evalualtion settings
!> \param alpha the given auxilliary basis function alpha
!> \param MAT the output matrix
!>
!> Calculates the 3 center 2 electron repulsion integral matrix (ab|0alpha) 
!> for a given auxilliary basis function alpha
!>
SUBROUTINE st_get_3center_eri_mat(LUPRI,LUERR,SETTING,ALPHA,MAT)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
TYPE(INTEGRALINPUT)   :: INT_INPUT
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
TYPE(MATRIX)          :: MAT
INTEGER               :: LUPRIOLD,LUPRI,LUERR,i,j,k,nbast,nbastAux,dim,alpha,a,b
INTEGER               :: c,ORB
TYPE(AOITEM),target   :: AO,emptyAO,AOaux,AOaux2

nbast=MAT%nrow

call LSHEADER(lupri,'st_get_3center_mat')

call init_integral_input(INT_INPUT,SETTING)

CALL BUILD_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
     &SETTING%MOLECULE(1)%p,SETTING%BASIS(1)%p%REGULAR,&
     &AO,'Default',.FALSE.,.FALSE.)
CALL BUILD_SINGLE_ORBBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
     &SETTING%MOLECULE(1)%p,SETTING%BASIS(1)%p%AUXILIARY,&
     &AOaux,'Default',alpha,dim,ORB)
! dim is the degeneracy so 1 for s orbital, 3 for p, 5 for d,....
! ORB is the orbitalcomponent 
CALL BUILD_EMPTY_AO(emptyAO,LUPRI)

     INT_INPUT%operator = 'Coulomb'
     INT_INPUT%AO(1)%p => AO
     INT_INPUT%AO(2)%p => AO
     INT_INPUT%sameLHSaos = .TRUE.
     INT_INPUT%AO(3)%p => emptyAO
     INT_INPUT%AO(4)%p => AOaux

NULLIFY(Int_Output%ResultMat)
ALLOCATE(Int_Output%ResultMat(nbast,nbast,1,dim,1))
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,&
     &INT_INPUT,INT_OUTPUT)

CALL DCOPY(nbast*nbast,Int_Output%ResultMat(:,:,1,ORB,1),1,MAT%elms,1)

DEALLOCATE(Int_Output%ResultMat)
CALL free_aoitem(lupri,AOaux)
CALL free_aoitem(lupri,AO)
CALL free_aoitem(lupri,emptyAO)

END SUBROUTINE st_get_3center_eri_mat

!> \brief Calculates the 3 center 2 electron repulsion integral matrix (ab|0alpha) for a given batch of auxilliary basis function alpha
!> \author T. Kjaergaard
!> \date 2010
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
!> \param SETTING Integral evalualtion settings
!> \param alpha the given auxilliary basis function alpha
!> \param dim is the dimension of the basisfunction alpha (1 for s orb, ..., 5 for d orb, ..)
!> \param MAT the output matrix
!>
!> Calculates the 3 center 2 electron repulsion integral matrix (ab|0alpha) 
!> for a given set of auxilliary basis function alpha
!>
SUBROUTINE st_get_3center_eri_matbatch(LUPRI,LUERR,SETTING,ALPHA,DIM,MAT)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
TYPE(INTEGRALINPUT)   :: INT_INPUT
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
INTEGER               :: LUPRIOLD,LUPRI,LUERR,i,j,k,nbast,nbastAux,dim,alpha,a,b
TYPE(MATRIX)          :: MAT(DIM)
INTEGER               :: c,ORB,dim2
TYPE(AOITEM),target   :: AO,emptyAO,AOaux,AOaux2

nbast=MAT(1)%nrow

call LSHEADER(lupri,'st_get_3center_mat')

call init_integral_input(INT_INPUT,SETTING)

CALL BUILD_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
     &SETTING%MOLECULE(1)%p,SETTING%BASIS(1)%p%REGULAR,&
     &AO,'Default',.FALSE.,.FALSE.)
CALL BUILD_SINGLE_ORBBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
     &SETTING%MOLECULE(1)%p,SETTING%BASIS(1)%p%AUXILIARY,&
     &AOaux,'Default',alpha,dim2,ORB)
! dim is the degeneracy so 1 for s orbital, 3 for p, 5 for d,....
! ORB is the orbitalcomponent 
IF(dim .NE. dim2 .OR. ORB .NE. 1)THEN
   WRITE(LUPRI,*)'SOMETHING WRONG IN st_get_3center_eri_matbatch '
   IF(dim .NE. dim2)THEN
      WRITE(LUPRI,*)'THE SPECIFIED DIMENSION AND THE ACTUAL DIMENSION&
           & ARE DIFFERENT'
   ENDIF
   IF(ORB .NE. 1)THEN
      WRITE(LUPRI,*)'THE SPECIFIED BASISFUNCTION SHOULD BE THE FIRST&
           & OF A GIVEN BATCH OF FOR INSTANCE THE 3 P ORBITALS'
   ENDIF
   CALL LSQUIT('ERROR IN st_get_3center_eri_matbatch',lupri)
ENDIF
CALL BUILD_EMPTY_AO(emptyAO,LUPRI)

     INT_INPUT%operator = 'Coulomb'
     INT_INPUT%AO(1)%p => AO
     INT_INPUT%AO(2)%p => AO
     INT_INPUT%sameLHSaos = .TRUE.
     INT_INPUT%AO(3)%p => emptyAO
     INT_INPUT%AO(4)%p => AOaux

NULLIFY(Int_Output%ResultMat)
ALLOCATE(Int_Output%ResultMat(nbast,nbast,1,dim,1))
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,&
     &INT_INPUT,INT_OUTPUT)

DO a=1,dim
   CALL DCOPY(nbast*nbast,Int_Output%ResultMat(:,:,1,a,1),1,MAT(a)%elms,1)
ENDDO

DEALLOCATE(Int_Output%ResultMat)
CALL free_aoitem(lupri,AOaux)
CALL free_aoitem(lupri,AO)
CALL free_aoitem(lupri,emptyAO)

END SUBROUTINE st_get_3center_eri_matbatch

END MODULE ST_INTEGRAL_INTF

!> \brief Does the SPHERICAL TRANSFORMATION of a 5 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TEMP the input 5 dim array i cartesian 
!> \param SPHINT the input 5 dim array i spherical
!> \param ndim1 the size of dimension 1
!> \param ndim2 the size of dimension 2
!> \param nMAT the number of input matrices
!> \param nsphmat the number og output matrices
!> \param nangmom the number of angular moments
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param LUPRI the logical unit number for the output file
SUBROUTINE SPHERICAL_TRANSFORMATION(TEMP,SPHINT,ndim1,ndim2,nMAT,nSPHMAT,nANGMOM,iprint,lupri)
use precision
IMPLICIT NONE
INTEGER     :: ANGMOM,nANGMOM,nLM,nXYZ,lupri,ELM,nMAT,nSPHMAT,iprint,ndim1,ndim2
REAL(REALK) :: TEMP(ndim1,ndim2,1,1,nMAT),SPHINT(ndim1,ndim2,1,1,nSPHMAT),COEF
REAL(REALK),pointer :: TRANSMAT(:)
REAL(REALK),PARAMETER :: D0 = 0.D0, D1 = 1.D0, D2 = 2.D0, D3 = 3.0D0
INTEGER     :: STARTI,STARTJ,I,J,K1,K2,TSIZE
nLM  = 2*nANGMOM + 1
nXYZ = (nANGMOM + 1)*(nANGMOM + 2)/2
!call mem_alloc(TRANSMAT,nXYZ,nLM)
TSIZE=nXYZ*nLM
ALLOCATE(TRANSMAT(nXYZ*nLM))

STARTJ=0
STARTI=0
DO ANGMOM=0,nANGMOM
   nLM  = 2*ANGMOM + 1
   nXYZ = (ANGMOM + 1)*(ANGMOM + 2)/2
   CALL BUILD_CART_TO_SPH_MAT(ANGMOM,TRANSMAT,TSIZE,nLM,nXYZ,lupri,iprint)
   DO I = 1, nLM
      DO J = 1, nXYZ
         COEF = TRANSMAT(J+(I-1)*nXYZ)
         IF (ABS(COEF) .GT. D0) THEN
            DO K2 = 1, ndim2
               DO K1 = 1, ndim1
                  SPHINT(K1,K2,1,1,STARTI+I) = SPHINT(K1,K2,1,1,STARTI+I) + COEF*TEMP(K1,K2,1,1,STARTJ+J)
               ENDDO
            ENDDO
         END IF
      ENDDO
   ENDDO
   STARTI=STARTI+nLM
   STARTJ=STARTJ+nXYZ
ENDDO
!call mem_dealloc(TRANSMAT)
DEALLOCATE(TRANSMAT)

END SUBROUTINE SPHERICAL_TRANSFORMATION

!> \brief Does the SPHERICAL TRANSFORMATION of a 2 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TEMP the input 2 dim array i cartesian 
!> \param SPHINT the input 2 dim array i spherical
!> \param ELM the size of dimension 1
!> \param nMAT the size of input dimension 2 (number of cartesian comp) 
!> \param nSPHMAT the size of output dimension 2 (number of spherical comp) 
!> \param nangmom the number of angular moments
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param LUPRI the logical unit number for the output file
SUBROUTINE SPHERICAL_TRANSFORMATION2(TEMP,SPHINT,ELM,nMAT,nSPHMAT,nANGMOM,iprint,lupri)
  use precision
IMPLICIT NONE
INTEGER     :: ANGMOM,nANGMOM,nLM,nXYZ,lupri,ELM,nMAT,nSPHMAT,iprint,TSIZE
REAL(REALK) :: TEMP(ELM,nMAT),SPHINT(ELM,nSPHMAT),COEF
REAL(REALK),pointer :: TRANSMAT(:)
REAL(REALK),PARAMETER :: D0 = 0.D0, D1 = 1.D0, D2 = 2.D0, D3 = 3.0D0
INTEGER     :: STARTI,STARTJ,I,J,K
nLM  = 2*nANGMOM + 1
nXYZ = (nANGMOM + 1)*(nANGMOM + 2)/2
!CALL MEM_ALLOC(TRANSMAT,nXYZ*nLM)
TSIZE=nXYZ*nLM
ALLOCATE(TRANSMAT(TSIZE))
STARTJ=0
STARTI=0
DO ANGMOM=0,nANGMOM
   nLM  = 2*ANGMOM + 1
   nXYZ = (ANGMOM + 1)*(ANGMOM + 2)/2
   CALL BUILD_CART_TO_SPH_MAT(ANGMOM,TRANSMAT,TSIZE,nLM,nXYZ,lupri,iprint)
   DO I = 1, nLM
      DO J = 1, nXYZ
         COEF = TRANSMAT(J+(I-1)*nXYZ)
         IF (ABS(COEF) .GT. D0) THEN
            DO K = 1, ELM
               SPHINT(K,STARTI+I) = SPHINT(K,STARTI+I) + COEF*TEMP(K,STARTJ+J)
            ENDDO
         END IF
      ENDDO
   ENDDO
   STARTI=STARTI+nLM
   STARTJ=STARTJ+nXYZ
ENDDO
!call mem_dealloc(TRANSMAT)
DEALLOCATE(TRANSMAT)
END SUBROUTINE SPHERICAL_TRANSFORMATION2

!> \brief build cartesian to spherical matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param ANGMOM the angular momentum
!> \param TRAMAT the output matrix
!> \param TSIZE the dim of TRAMAT
!> \param nLM the number of spherical components
!> \param nXYZ the number of cartesian components
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE BUILD_CART_TO_SPH_MAT(ANGMOM,TRAMAT,TSIZE,nLM,nXYZ,lupri,iprint)
use precision
use math_fun
IMPLICIT NONE
INTEGER               :: LVAL1,ANGMOM,nLM,nXYZ,IPRINT,lupri,TSIZE
REAL(REALK),PARAMETER :: D0 = 0.D0, D1 = 1.D0, D2 = 2.D0, D3 = 3.0D0
REAL(REALK)           :: TRAMAT(TSIZE), COSSIN(0:ANGMOM,0:ANGMOM), PL(0:ANGMOM)
INTEGER               :: K,M,IX,I,J,L,ILM,N,IXYZ
REAL(REALK)           :: CM,CMK,CMKIJ,CMKI,FAC
LVAL1 = ANGMOM+1
CALL LS_DZERO(PL(0),LVAL1)
DO K = 0, ANGMOM/2
   FAC = (-1)**K
   PL(ANGMOM-2*K) = (FAC/(2**ANGMOM))&
   &  *BINOM(lupri,ANGMOM,K)*BINOM(lupri,2*(ANGMOM-K),ANGMOM)
ENDDO
!IF(IPRINT .GT. 10)THEN
!   CALL LSHEADER(lupri,'Legendre polynomial')
!   CALL OUTPUT(PL(0:LVAL1),1,1,1,LVAL1,1,LVAL1,1,LUPRI)
!ENDIF
CALL LS_DZERO(COSSIN(0,0),LVAL1*LVAL1)
DO M = 0, ANGMOM
   COSSIN(M,0) = D1
   DO K = 1, M
      COSSIN(M,K) = COSSIN(M-1,K-1)*((-1)**(K-1))
      IF (M .GT. K) COSSIN(M,K) = COSSIN(M,K) + COSSIN(M-1,K)
   ENDDO
ENDDO
!IF(IPRINT .GT. 10)THEN
!   CALL LSHEADER(lupri,'Cosine and sine factors')
!   CALL OUTPUT(COSSIN(0:LVAL1,0:LVAL1),1,LVAL1,1,LVAL1,LVAL1,LVAL1,1,LUPRI)
!ENDIF

! Transformation coefficients
! ---------------------------
!

CALL LS_DZERO(TRAMAT,nXYZ*nLM)
DO M = 0, ANGMOM
   CM = SQRT(D2*FACULT(LUPRI,ANGMOM-M)/FACULT(LUPRI,ANGMOM+M))
   IF (M .EQ. 0) CM = D1
!   IF (MINTEG.EQ.2) CM = CM/SQRT(FACUL2(2*ANGMOM-1))
   DO K = MOD(ANGMOM - M,2), ANGMOM - M, 2
      IF (M .GT. 0) PL(K) = (K+1)*PL(K+1)
      CMK = CM*PL(K)
      DO I = 0, (ANGMOM - K - M)/2
         CMKI = CMK*BINOM(lupri,(ANGMOM - K - M)/2,I)
         DO J = 0, I
            CMKIJ = CMKI*BINOM(lupri,I,J)
            DO N = 0, M
               IX = ANGMOM - 2*J - M + N
               IX = IX*(IX + 1)/2 + LVAL1 - M - 2*I
!               IF (MORDER .EQ. 0) THEN
                  ILM = MAX(1,2*M + MOD(N,2))
!               ELSE
!                  IF (MOD(N,2) .EQ. 1) THEN
!                     ILM = 1 + ANGMOM - M
!                  ELSE
!                     ILM = 1 + ANGMOM + M
!                  END IF
!               END IF
               TRAMAT(IX+(ILM-1)*nXYZ) = TRAMAT(IX+(ILM-1)*nXYZ) + CMKIJ*COSSIN(M,N) 
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

!IF (IPRINT .GT. 4) THEN
!   CALL LSHEADER(lupri,'Cartesian to spherical transformation matrix')
!   WRITE (LUPRI,'(29X,A,I2)') ' Moment order:',ANGMOM
!   IXYZ = (ANGMOM+1)*(ANGMOM+2)/2
!   ILM  = 2*ANGMOM + 1
!   CALL OUTPUT(TRAMAT,1,IXYZ,1,ILM,NXYZ,NLM,1,LUPRI)
!END IF

END SUBROUTINE BUILD_CART_TO_SPH_MAT
