!> @file
!> Module that contains information about angular indeces
MODULE SphCart_Matrices
use precision
use math_fun

TYPE SPHMAT
REAL(REALK),pointer  :: elms(:)
END TYPE SPHMAT

TYPE SPHMATPOINTER
TYPE(SPHMAT),pointer :: p
END TYPE SPHMATPOINTER

CONTAINS
!> \brief calculate spherical transformation matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param SPH_MAT the spherical matrices 
!> \param nSPHMAT the number of spherical matrices
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE setPrecalculatedSphmat(SPH_MAT,nSPHMAT,LUPRI,IPRINT)
implicit none
TYPE(SPHMAT),pointer     :: SPH_MAT(:)
INTEGER                  :: nSPHMAT,LUPRI,IPRINT
!
INTEGER                  :: I,L,nLM,nXYZ

NULLIFY(SPH_MAT)
ALLOCATE(SPH_MAT(nSPHMAT))
DO I=1,nSPHMAT
  L = I-1
  nLM  = 2*L+1
  nXYZ = (L+1)*(L+2)/2
  NULLIFY(SPH_MAT(I)%elms)
  ALLOCATE(SPH_MAT(I)%elms(nLM*nXYZ))
  CALL Sph_to_Cart_matrix(L,SPH_MAT(I)%elms,nLM,nXYZ,LUPRI,IPRINT)
ENDDO

END SUBROUTINE setPrecalculatedSphmat

!> \brief free the calculated spherical transformation matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param SPH_MAT the spherical matrices 
!> \param nSPHMAT the number of spherical matrices
SUBROUTINE freePrecalculatedSphmat(SPH_MAT,nSPHMAT)
implicit none
TYPE(SPHMAT),pointer     :: SPH_MAT(:)
INTEGER                  :: nSPHMAT
!
INTEGER                  :: I

DO I=1,nSPHMAT
  DEALLOCATE(SPH_MAT(I)%elms)
  NULLIFY(SPH_MAT(I)%elms)
ENDDO
DEALLOCATE(SPH_MAT)
NULLIFY(SPH_MAT)

END SUBROUTINE freePrecalculatedSphmat

!> \brief calculate spherical to cartesian matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param L angular moment
!> \param SCMAT the transformation matrix
!> \param nLM the number of spherical components
!> \param nXYZ the number of cartesian components
!> \param lupri the logical unit number for the output file
!> \param iprint the printlevel, determining how much output should be generated
SUBROUTINE Sph_to_Cart_matrix(L,SCMAT,nLM,nXYZ,LUPRI,IPRINT)
implicit none
REAL(REALK) :: SCMAT(nLM,nXYZ)
INTEGER     :: L,nLM,nXYZ,LUPRI,IPRINT
!
Real(realk), parameter :: DM1 = -1.0d0, D1 = 1.0d0, D2 = 2.0d0
INTEGER     :: M1,MADR,MABS,V0
REAL(realk) :: FACNRM,FAC3,FACTOR
INTEGER     :: T,U,V,A,B,C,P,Q,R,X,Y,Z,M
INTEGER     :: iLM,iXYZ

IF(L .EQ. 0)THEN ! S
 SCMAT(1,1)=1.d0
ELSEIF(L .EQ. 1)THEN ! P
 CALL LS_DZERO(SCMAT,9)
 SCMAT(1,1)=1.d0
 SCMAT(2,2)=1.d0
 SCMAT(3,3)=1.d0
ELSEIF(L .GT. 1)THEN
 CALL LS_DZERO(SCMAT,nLM*nXYZ)
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
           iLM = 1 + L + MADR
           iXYZ = NCRT(X,Y,Z)
           SCMAT(iLM,iXYZ) = SCMAT(iLM,iXYZ) + FACTOR 
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
ENDIF

END SUBROUTINE Sph_to_Cart_matrix

!> \brief calculate cartesian to spherical matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param L angular moment
!> \param CSMAT the transformation matrix
!> \param nXYZ the number of cartesian components
!> \param nLM the number of spherical components
!> \param lupri the logical unit number for the output file
!> \param iprint the printlevel, determining how much output should be generated
SUBROUTINE Cart_to_Sph_matrix(L,CSMAT,nXYZ,nLM,LUPRI,IPRINT)
implicit none
REAL(REALK) :: CSMAT(nXYZ,nLM)
INTEGER     :: L,nXYZ,nLM,LUPRI,IPRINT
!
Real(realk), parameter :: DM1 = -1.0d0, D1 = 1.0d0, D2 = 2.0d0
INTEGER     :: M1,MADR,MABS,V0
REAL(realk) :: FACNRM,FAC3,FACTOR
INTEGER     :: T,U,V,A,B,C,P,Q,R,X,Y,Z,M
INTEGER     :: iXYZ,iLM

IF(L .EQ. 0)THEN ! S
  CSMAT(1,1)=1.d0
ELSEIF(L .EQ. 1)THEN ! P
  CALL LS_DZERO(CSMAT,9)
  CSMAT(1,1)=1.d0
  CSMAT(2,2)=1.d0
  CSMAT(3,3)=1.d0
ELSEIF(L .GT. 1)THEN
  CALL LS_DZERO(CSMAT,nXYZ*nLM)
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
	  iXYZ = NCRT(X,Y,Z)
	  iLM = 1 + L + MADR
          CSMAT(iXYZ,iLM) = CSMAT(iXYZ,iLM) + FACTOR 
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
ENDIF

END SUBROUTINE Cart_to_Sph_matrix

END MODULE SphCart_Matrices
