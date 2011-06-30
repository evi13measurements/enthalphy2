!> @file 
!> timing module
MODULE LSTIMING
use precision
!use typedef
Integer,save :: LSIUNIT = 6
contains
!> \brief take time
!> \author T. Kjaergaard
!> \date 2010
!> \param text label to print along with timings
!> \param CPUTIME the cpu time
!> \param WALLTIME the wall time
!> \param lupri the logical unit number 
SUBROUTINE LSTIMER(TEXT,CPUTIME,WALLTIME,LUPRI)
implicit none
INTEGER           :: LUPRI,length
CHARACTER*(*)     :: TEXT
CHARACTER(len=15) :: PRINTTEXT
REAL(REALK)       :: TIME1,TIME2,DELTAWALL,CPUTIME,WALLTIME,DELTACPU

INTEGER :: IUNIT

IF (LUPRI.EQ.-1) THEN
  IUNIT = LSIUNIT
ELSE
  IUNIT = LUPRI
ENDIF

length = LEN(TEXT)

IF(length .GT. 15) CALL LSQUIT('TEXTLENGTH PROVIDED TO LSTIMER IS LIMITED TO 15',lupri)
IF (TEXT(1:5) .EQ. 'START') THEN
   CALL LS_GETTIM(CPUTIME,WALLTIME)
ELSE
   PRINTTEXT='               '
   CALL LS_GETTIM(TIME1,TIME2)
   DELTACPU=TIME1-CPUTIME
   PRINTTEXT(1:length) =  TEXT(1:length) 
   CALL LS_TIMTXT('>>>  CPU Time used in '//PRINTTEXT//' is',DELTACPU,IUNIT)
#ifdef VAR_WTIME
   DELTAWALL=TIME2-WALLTIME
   CALL LS_TIMTXT('>>> wall Time used in '//PRINTTEXT//' is',DELTAWALL,IUNIT)
   WALLTIME=TIME2
#endif
   CPUTIME=TIME1
ENDIF

END SUBROUTINE LSTIMER

!> \brief set the lsiunit integer
!> \author T. Kjaergaard
!> \date 2010
!> \param iunit
SUBROUTINE SET_LSIUNIT(IUNIT)
implicit none
INTEGER :: IUNIT
  LSIUNIT = IUNIT
END SUBROUTINE SET_LSIUNIT

END MODULE LSTIMING
