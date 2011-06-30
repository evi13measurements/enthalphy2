MODULE ls_util
Integer,save :: LSIUNIT = 6

CONTAINS
SUBROUTINE LSHEADER(LUPRI,HEAD)
IMPLICIT NONE
  CHARACTER HEAD*(*)
  INTEGER :: LUPRI,LENGTH,INDENT,I
  
  LENGTH = LEN(HEAD)
  INDENT = (72 - LENGTH)/2 + 1
  
  WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
  WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
  WRITE (LUPRI, '()')
  
END SUBROUTINE LSHEADER

SUBROUTINE LS_PRINT_GRADIENT(lupri,molecule,GRDMOL,natoms,TEXT)
  use molecule_type
  use precision  
  implicit none
  integer :: lupri,natoms
  type(moleculeinfo) :: molecule
  real(realk)       :: GRDMOL(3,natoms)
  CHARACTER*(*)     :: TEXT
  CHARACTER(len=15) :: PRINTTEXT
  CHARACTER(len=23) :: PRINTTEXT2
  CHARACTER(len=40) :: PRINTTEXT3
  INTEGER :: IUNIT,length,ioff,iatom,j
  IF (LUPRI.EQ.-1) THEN
     IUNIT = LSIUNIT
  ELSE
     IUNIT = LUPRI
  ENDIF

length = LEN(TEXT)
IF(length .GT. 15) CALL LSQUIT('TEXTLENGTH PROVIDED TO LS_PRINT_GRADIENT IS LIMITED TO 15',lupri)
IF (TEXT(1:5) .EQ. 'TOTAL') THEN
   PRINTTEXT2 = 'Molecular gradient (au)'
  CALL LSHEADER(iunit,PRINTTEXT2)
ELSE
   PRINTTEXT3 = TEXT(1:length)//' contribution to gradient'
   length = length+25
  CALL LSHEADER(iunit,PRINTTEXT3(1:length))
ENDIF

DO IATOM = 1, natoms
   WRITE (IUNIT,'(1X,A6,F17.10,2F24.10)') &
        &molecule%ATOM(IATOM)%Name, (GRDMOL(J,IATOM),J=1,3)
END DO

end SUBROUTINE LS_PRINT_GRADIENT

!> \brief Compute the 'almost' root-mean square of A-B, diff = RMS(A-B), where A,B are matrices
!> \author P. Merlot
!> \date 12-05-2010
!> \param diff The 'almost' RMS norm of the difference between A and B (almost, since divided by sqrt(nrow*ncol) and not by nrow*ncol) 
SUBROUTINE rms_Diff(A, B, nrow,ncol,diff)
  use precision
  implicit none
  INTEGER, intent(IN)      :: nrow,ncol
  REAL(realk), intent(IN)  :: A(nrow,ncol), B(nrow,ncol)
  REAL(realk), intent(OUT) :: diff
!
  REAL(realk)              :: tempRms(nrow,ncol)    
  REAL(realk),pointer      :: WORK
  REAL(realk), external    :: dlange

  tempRms = 0.0d0
  tempRms = abs(A)-abs(B)
  diff = dlange('F',nrow,ncol,tempRms,1,WORK) * sqrt(1.d0/(nrow*ncol))
END SUBROUTINE rms_Diff

END MODULE LS_UTIL

MODULE math_fun
use typedef
use precision

contains
FUNCTION FACULT(LUPRI,N)
IMPLICIT NONE
real(realk), PARAMETER :: D1=1.D0
integer        :: N,I,LUPRI
real(realk)    :: FACULT
IF (N .LT. 0) THEN
   WRITE (LUPRI,'(/,A,I10,/A)')&
   &         ' Argument less than zero in FACULT:',N,&
   &         ' Program cannot continue.'
   CALL LSQUIT('Illegal argument in FACULT',lupri)
ELSE
   FACULT = D1
   DO I = 1, N
      FACULT = FACULT*I
   ENDDO
END IF
END FUNCTION FACULT

FUNCTION FACUL2(LUPRI,N)
IMPLICIT NONE
real(realk), PARAMETER :: D1=1.D0
real(realk)    :: FACUL2
integer :: N,I,LUPRI
IF (N .LT. 0) THEN
   FACUL2 = DFLOAT(N + 2)
   DO I = N + 4, 1, 2
      FACUL2 = FACUL2*I
   END DO
   IF (FACUL2 .EQ. 0.d0) THEN
      WRITE (LUPRI,'(/,A,I10,/A)')&
      &            ' Double factorial undefined for ',N,&
      &            ' Program cannot continue.'
      CALL LSQUIT('Illegal argument in FACUL2',lupri)
   ELSE
      FACUL2 = D1/FACUL2
   END IF
ELSE IF (N.EQ.0) THEN
   FACUL2 = D1
ELSE ! N > 0
   FACUL2 = DFLOAT(N)
   DO I = N - 2, 1, -2
      FACUL2 = FACUL2*I
   END DO
END IF
END FUNCTION FACUL2

FUNCTION BINOM(LUPRI,I,J)
real(realk), PARAMETER  :: D1=1.D0
INTEGER :: I,J,LUPRI
real(realk)    :: BINOM

IF (I .LT. J) THEN
   WRITE (LUPRI,'(/,A,2I5,/A)')&
   &         ' Second argument larger than first argument in BINOM:',&
   &         I,J,' Program cannot continue.'
   CALL LSQUIT('Illegal arguments in BINOM',lupri)
ELSE
   BINOM = FACULT(LUPRI,I)/(FACULT(LUPRI,I-J)*FACULT(LUPRI,J))
END IF
END FUNCTION BINOM

FUNCTION NCRT(I,J,K)
IMPLICIT NONE
INTEGER  :: I,J,K,NCRT
NCRT = 1 + J + 2*K + (J + K)*(J + K - 1)/2
END FUNCTION NCRT


END MODULE math_fun

