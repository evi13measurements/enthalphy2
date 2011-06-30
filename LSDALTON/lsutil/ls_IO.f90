!> @file 
!> Contain the IOITEM which keeps track of files written to disk
MODULE io
use files
use precision
!Integer,parameter :: maxFiles = 30000
Integer,parameter :: maxRecord = 2**28
Integer,parameter :: increment = 50

TYPE IOITEM
Integer       :: nallocFiles
Integer       :: numFiles
Character(80),pointer :: filename(:)
Integer,pointer       :: IUNIT(:)
Logical,pointer       :: isOpen(:)
END TYPE IOITEM

CONTAINS
!> \brief initialise the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
SUBROUTINE io_init(IO)
implicit none
TYPE(IOITEM)  :: IO
IO%numFiles = 0
call io_alloc(IO,increment)
END SUBROUTINE io_init

!> \brief free the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
SUBROUTINE io_free(IO)
implicit none
TYPE(IOITEM)  :: IO
IO%numFiles = 0
call io_dealloc(IO)
END SUBROUTINE io_free

!> \brief allocate the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
!> \param nsize the number of files that should be allocated
SUBROUTINE io_alloc(IO,nsize)
implicit none
integer       :: nsize
TYPE(IOITEM)  :: IO
IO%nallocFiles = nsize
nullify(IO%filename)
nullify(IO%IUNIT)
nullify(IO%isopen)
allocate(IO%filename(nsize))
allocate(IO%IUNIT(nsize))
allocate(IO%isopen(nsize))
END SUBROUTINE io_alloc

!> \brief deallocate the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
SUBROUTINE io_dealloc(IO)
implicit none
TYPE(IOITEM)  :: IO
IO%nallocFiles = 0
deallocate(IO%filename)
deallocate(IO%IUNIT)
deallocate(IO%isopen)
nullify(IO%filename)
nullify(IO%IUNIT)
nullify(IO%isopen)
END SUBROUTINE io_dealloc

!> \brief copy the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM, to be copied
!> \param IONEW the new copied IOITEM
SUBROUTINE COPY_IOITEM(IO,IONEW)
IMPLICIT NONE
TYPE(IOITEM) :: IO,IONEW
INTEGER      :: I
IOnew%numfiles=IO%numfiles
IF(IOnew%numfiles.GT.IOnew%nallocFiles)&
     &call lsquit('Programming Error in copy_ioitem, not initialized IONEW')
DO I=1,IO%numfiles
 IOnew%filename(I)=IO%filename(I)
 IOnew%IUNIT(I)=IO%IUNIT(I)
 IOnew%isopen(I)=IO%isopen(I)
ENDDO

END SUBROUTINE COPY_IOITEM

!> \brief add a filename to the IOITEM structue
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
!> \param Filename the filename to be added
!> \param LUPRI the logical unit number for the output file
SUBROUTINE io_add_filename(IO,Filename,LUPRI)
implicit none
TYPE(IOITEM)  :: IO
Character(80) :: Filename
Integer       :: LUPRI
!
Logical :: fileFound
Integer :: iFile,oldnallocFiles
TYPE(IOITEM)  :: TEMPIO

IF (io_file_exist(Filename,IO)) THEN
  WRITE(LUPRI,'(1X,3A)') 'Error in io_add_filename. File ',TRIM(Filename),&
   &                      ' allready exist in io-list!'
  CALL lsQUIT('Error in io_add_filename. Trying to add an existing filename.',lupri)
ENDIF

IF(IO%numFiles+1.GT.IO%nallocFiles)THEN
   oldnallocFiles = IO%nallocFiles
   call IO_alloc(TEMPIO,oldnallocFiles)
   CALL COPY_IOITEM(IO,TEMPIO)
   call IO_dealloc(IO)
   call IO_alloc(IO,oldnallocFiles+increment)
   CALL COPY_IOITEM(TEMPIO,IO)
   call IO_dealloc(TEMPIO)
   IF(IO%numFiles+1.GT.IO%nallocFiles)THEN
      WRITE(LUPRI,'(1X,2A)') 'Error in io_add_filename. something &
           & strange happend, last add file =',TRIM(Filename)
      CALL lsQUIT('Error in io_add_filename',lupri)
   ENDIF
endif

IO%numFiles = IO%numFiles + 1
IO%filename(IO%numFiles) = Filename
IO%IUNIT(IO%numFiles) = -1
IO%isOpen(IO%numFiles) = .FALSE.

END SUBROUTINE io_add_filename

!> \brief determines if the file exist in the IOITEM
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Filename the filename to be added
!> \param IO the IOITEM
LOGICAL FUNCTION io_file_exist(Filename,IO)
implicit none
Character(80) :: Filename
TYPE(IOITEM)  :: IO
!
Integer :: iFile
io_file_exist = .FALSE.
DO iFile=1,IO%numFiles
  IF (Filename.EQ.IO%Filename(iFile)) THEN
    io_file_exist = .TRUE.
    RETURN
  ENDIF
ENDDO
END FUNCTION io_file_exist

!> \brief determines logical unit number of filename
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
!> \param Filename the filename
INTEGER FUNCTION io_iunit(IO,Filename)
implicit none
Character(80) :: Filename
TYPE(IOITEM)  :: IO
!
Integer :: iFile

IF (.not.io_file_exist(Filename,IO)) THEN
  CALL lsQUIT('Error in io_iunit, file does not exist!',-1)
ENDIF

DO iFile=1,IO%numFiles
  IF (Filename.EQ.IO%Filename(iFile)) EXIT
ENDDO
io_iunit = IO%iUnit(iFile)
END FUNCTION io_iunit

!> \brief open a filename
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_open(Filename,IO,LUPRI,LUERR)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR
TYPE(IOITEM)  :: IO
!
Integer :: IDUM,IUNIT,iFile
Logical :: fileFound,isOpen,LDUM

fileFound = .FALSE.
DO iFile=1,IO%numFiles
  IF (Filename.EQ.IO%Filename(iFile)) THEN
    fileFound = .TRUE.
    IUNIT     = IO%IUNIT(iFile)
    isOpen    = IO%isOpen(iFile)
    EXIT
  ENDIF
ENDDO

IF (.NOT. fileFound) THEN
  WRITE(LUPRI,'(1X,2A)') 'Error in io_open. Could not find file: ',TRIM(Filename)
  CALL lsQUIT('Error in io_open. File not found!',lupri)
ENDIF

IF (isOpen) THEN
  WRITE(LUPRI,'(1X,3A)') 'Error in io_open. File: ',TRIM(Filename),' allready opened.'
  CALL lsQUIT('Error in io_open. Trying to open an allready opened file!',lupri)
ENDIF

CALL LSOPEN(IUNIT,FILENAME,'UNKNOWN','UNFORMATTED')

IO%isOpen(iFile) = .TRUE.
IO%IUNIT(iFile)  = IUNIT

END SUBROUTINE io_open

!> \brief close filename
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_close(Filename,IO,LUPRI,LUERR)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR
TYPE(IOITEM)  :: IO
!
Integer :: IDUM,IUNIT,iFile
Logical :: fileFound,isOpen,LDUM,fileExist

fileFound = .FALSE.
DO iFile=1,IO%numFiles
  IF (Filename.EQ.IO%Filename(iFile)) THEN
    fileFound = .TRUE.
    IUNIT     = IO%IUNIT(iFile)
    isOpen    = IO%isOpen(iFile)
    EXIT
  ENDIF
ENDDO

IF (.NOT. fileFound) THEN
  WRITE(LUPRI,'(1X,2A)') 'Error in io_close. Could not find file: ',TRIM(Filename)
  CALL lsQUIT('Error in io_close. File not found!',lupri)
ENDIF

IF (.NOT.isOpen) THEN
  WRITE(LUPRI,'(1X,3A)') 'Error in io_close. File: ',TRIM(Filename),' not opened.'
  CALL lsQUIT('Error in io_close. Trying to close a file that is not open!',lupri)
ENDIF

INQUIRE(file=FILENAME,exist=fileExist)

IF (.NOT.fileExist) THEN
  WRITE(LUPRI,'(1X,3A)') 'Error in io_close. File: ',TRIM(Filename),' does not exist on disk!'
  CALL lsQUIT('Error in io_close. File missing!',lupri)
ENDIF

CALL LSCLOSE(IUNIT,'KEEP')

IO%isOpen(iFile) = .FALSE.
IO%IUNIT(iFile)  = -1

END SUBROUTINE io_close

!> \brief obtain a filename based on input values used for screening
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Filename the filename
!> \param Identifier Character string to identify the filename
!> \param AO1 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 1
!> \param AO2 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 2
!> \param AO3 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 3
!> \param AO4 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 4
!> \param start1 start index for AO1
!> \param start2 start index for AO2
!> \param Oper operator label
!> \param intType 'Contracted or Primitive
!> \param FRAGMENT is this a fragment calculation
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_get_filename(Filename,Identifier,AO1,AO2,AO3,AO4,start1,start2,Oper,intType,FRAGMENT,LUPRI,LUERR)
implicit none
Character*(*)        :: Identifier,AO1,AO2,AO3,AO4,Oper,intType
Integer              :: LUPRI,LUERR,start1,start2
Character(80)        :: Filename
logical              :: FRAGMENT
!
Character(80) :: AOstring(4)
Integer       :: iLen,iFilename,IAO,i,tmp(2)
Character(len=1)  :: STRING1
Character(len=2)  :: STRING2
Character(len=3)  :: STRING3
Character(len=4)  :: STRING4
Character(len=5)  :: STRING5
Character(len=6)  :: STRING6

AOstring(1) = AO1
AOstring(2) = AO2
AOstring(3) = AO3
AOstring(4) = AO4

iLen = LEN(TRIM(Identifier))
Filename(1:iLen) = Identifier(1:iLen)
iFilename = iLen + 1
iLen = LEN(TRIM(intType))
iLen = min(iLen,3)
Filename(iFilename:iFilename+iLen-1) = intType(1:iLen)
iFilename = iFilename + iLen
iLen = LEN(TRIM(Oper))
iLen = min(iLen,3)
Filename(iFilename:iFilename+iLen-1) = Oper(1:iLen)
iFilename = iFilename + iLen
DO IAO=1,4
  IF (AOstring(IAO).EQ.'Regular') THEN
    Filename(iFilename:iFilename) = 'r'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.'DF-Aux') THEN
    Filename(iFilename:iFilename) = 'd'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.'Empty') THEN
    Filename(iFilename:iFilename) = 'e'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.'Nuclear') THEN
    Filename(iFilename:iFilename) = 'n'
    iFilename = iFilename + 1
  ELSE
    WRITE(LUPRI,'(1X,A,I1,2A)') 'Error! Wrong AOstring(',IAO,')=',AOstring(IAO),' in io_get_filename.'
    CALL lsQUIT('Error in io_get_filename',lupri)
  ENDIF
ENDDO

IF(FRAGMENT)THEN
   Filename(iFilename:iFilename) = 'F'
   iFilename = iFilename + 1
ENDIF

tmp(1)=start1
tmp(2)=start2
do I=1,2
   Filename(iFilename:iFilename) = 'S'
   iFilename = iFilename + 1
   IF(tmp(I) .LT. 10)THEN
      WRITE(STRING1,'(I1)') tmp(I)
      Filename(iFilename:iFilename) = STRING1
      iFilename = iFilename + 1
   ELSEIF(tmp(I) .LT. 100)THEN
      WRITE(STRING2,'(I2)') tmp(I)
      Filename(iFilename:iFilename+1) = STRING2
      iFilename = iFilename + 2
   ELSEIF(tmp(I) .LT. 1000)THEN
      WRITE(STRING3,'(I3)') tmp(I)
      Filename(iFilename:iFilename+2) = STRING3
      iFilename = iFilename + 3
   ELSEIF(tmp(I) .LT. 10000)THEN
      WRITE(STRING4,'(I4)') tmp(I)
      Filename(iFilename:iFilename+3) = STRING4
      iFilename = iFilename + 4
   ELSEIF(tmp(I) .LT. 100000)THEN
      WRITE(STRING5,'(I5)') tmp(I)
      Filename(iFilename:iFilename+4) = STRING5
      iFilename = iFilename + 5
   ELSEIF(tmp(I) .LT. 1000000)THEN
      WRITE(STRING6,'(I6)') tmp(I)
      Filename(iFilename:iFilename+5) = STRING6
      iFilename = iFilename + 6
   ELSE
      WRITE(LUPRI,'(1X,A,I1,2A)') 'startbasisfunction in io_get_filename is larger than 100000 not implemented'
      CALL lsQUIT('Error in io_get_filename',lupri)
   ENDIF
enddo

iFilename=iFilename-1
IF (iFilename.GT.80) THEN
  CALL lsQUIT('Error: iFilename > 80 in io_get_filename.',lupri)
ELSE
  DO i=iFilename+1,80
    Filename(i:i) = ' '
  ENDDO
ENDIF
END SUBROUTINE io_get_filename

!> \brief write mat to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat to be written to disk
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_write_mat(Mat,Filename,IO,LUPRI,LUERR)
use matrix_operations
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR,n1,n2,n3,n4,n5
TYPE(matrix)  :: Mat
TYPE(IOITEM)  :: IO
CALL io_open(Filename,IO,LUPRI,LUERR)
CALL mat_write_to_disk(io_iunit(IO,Filename),Mat)
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_write_mat

!> \brief write tensor to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat 5 dim array to be written to disk
!> \param n1 the size of dimension 1
!> \param n2 the size of dimension 2
!> \param n3 the size of dimension 3
!> \param n4 the size of dimension 4
!> \param n5 the size of dimension 5
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_write(Mat,n1,n2,n3,n4,n5,Filename,IO,LUPRI,LUERR)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR,n1,n2,n3,n4,n5
Real(realk)   :: Mat(n1,n2,n3,n4,n5)
TYPE(IOITEM)  :: IO
CALL io_open(Filename,IO,LUPRI,LUERR)
CALL io_write_tensor(Mat,n1,n2,1,1,1,io_iunit(IO,Filename))
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_write

!> \brief write tensor to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat 5 dim array to be written to disk
!> \param n1 the size of dimension 1
!> \param n2 the size of dimension 2
!> \param n3 the size of dimension 3
!> \param n4 the size of dimension 4
!> \param n5 the size of dimension 5
!> \param IUNIT the logical unit number for the file to write to
SUBROUTINE io_write_tensor(Mat,n1,n2,n3,n4,n5,IUNIT)
implicit none
Real(realk)   :: Mat(n1,n2,n3,n4,n5)
Integer       :: n1,n2,n3,n4,n5,IUNIT
!
Integer :: i1,i2,i3,i4,i5
Integer :: IDUM
Logical :: LDUM
WRITE(IUNIT) n1,n2,n3,n4,n5
DO i5=1,n5
  DO i4=1,n4
    DO i3=1,n3
      DO i2=1,n2
        WRITE(IUNIT) (Mat(i1,i2,i3,i4,i5),i1=1,n1)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE io_write_tensor

!> \brief write lstensor to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param tensor the lstensor to be written to disk  
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_write_lstensor(tensor,Filename,IO,LUPRI,LUERR)
use LSTENSOR_OPERATIONSMOD, only: write_lstensor_to_disk, LSTENSOR
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR
type(lstensor) :: tensor
TYPE(IOITEM)  :: IO
CALL io_open(Filename,IO,LUPRI,LUERR)
CALL write_lstensor_to_disk(tensor,io_iunit(IO,Filename),lupri)
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_write_lstensor

!> \brief read mat to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat to be read from file 
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_read_mat(Mat,Filename,IO,LUPRI,LUERR)
use matrix_operations
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR,n1,n2,n3,n4,n5
TYPE(matrix)  :: Mat
TYPE(IOITEM)  :: IO
CALL io_open(Filename,IO,LUPRI,LUERR)
CALL mat_read_from_disk(io_iunit(IO,Filename),Mat)
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_read_mat

!> \brief read 5 dim array to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat to be read from file 
!> \param n1 the size of dimension 1
!> \param n2 the size of dimension 2
!> \param n3 the size of dimension 3
!> \param n4 the size of dimension 4
!> \param n5 the size of dimension 5
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_read(Mat,n1,n2,n3,n4,n5,Filename,IO,LUPRI,LUERR)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR,n1,n2,n3,n4,n5
Real(realk)   :: Mat(n1,n2,n3,n4,n5)
TYPE(IOITEM)  :: IO
CALL io_open(Filename,IO,LUPRI,LUERR)
CALL io_read_tensor(Mat,n1,n2,1,1,1,io_iunit(IO,Filename))
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_read

!> \brief read 5 dim array from disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat to be read from file 
!> \param n1 the size of dimension 1
!> \param n2 the size of dimension 2
!> \param n3 the size of dimension 3
!> \param n4 the size of dimension 4
!> \param n5 the size of dimension 5
!> \param IUNIT the logical unit number for file from which to read
SUBROUTINE io_read_tensor(Mat,n1,n2,n3,n4,n5,IUNIT)
implicit none
Real(realk)   :: Mat(n1,n2,n3,n4,n5)
Integer       :: n1,n2,n3,n4,n5,IUNIT
!
Integer :: i1,i2,i3,i4,i5
Integer :: IDUM
Logical :: LDUM
READ(IUNIT) i1,i2,i3,i4,i5
IF ((i1.NE.n1).OR.(i2.NE.n2).OR.(i3.NE.n3).OR.(i4.NE.n4).OR.(i5.NE.n5)) THEN
  CALL lsQUIT('Dimension error in io_read_tensor!',-1)
ENDIF
DO i5=1,n5
  DO i4=1,n4
    DO i3=1,n3
      DO i2=1,n2
        READ(IUNIT) (Mat(i1,i2,i3,i4,i5),i1=1,n1)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE io_read_tensor

!> \brief read lstensor from disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param tensor the lstensor to be written to disk  
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_read_lstensor(tensor,Filename,IO,LUPRI,LUERR)
use LSTENSOR_OPERATIONSMOD, only: read_lstensor_from_disk, LSTENSOR
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR
type(lstensor) :: tensor
TYPE(IOITEM)  :: IO
CALL io_open(Filename,IO,LUPRI,LUERR)
CALL read_lstensor_from_disk(tensor,io_iunit(IO,Filename),lupri)
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_read_lstensor

!Type(IOitem)         :: IO
!IUNIT = -1
!CALL LSOPEN(IUNIT,FILENAME,'UNKNOWN','SEQUENTIAL','UNFORMATTED',IDUM,LDUM)

  
END MODULE io

!> \brief write vector 
!> \author S. Reine
!> \date 2010
!> \param iunit the logical unit number for the file, from which to read
!> \param vector the vectro th write
!> \param N the size of the vector
SUBROUTINE ls_write(IUNIT,vector,N)
use io, only: maxRecord
use precision
implicit none
Integer,intent(IN)     :: IUNIT,N
Real(realk),intent(IN) :: vector(N)
!
Integer :: nBUF,startBUF,endBUF,I
nBUF     = n/maxRecord + 1
startBUF = 1
DO I=1,nBUF
  endBUF = MIN(n,startBUF+maxRecord-1)
  WRITE(IUNIT) vector(startBUF:endBUF)
  startBUF = startBUF + maxRecord
ENDDO
END SUBROUTINE ls_write

!> \brief read vector
!> \author S. Reine
!> \date 2010
!> \param iunit the logical unit number for the file, from which to read
!> \param vector the vectro th write
!> \param N the size of the vector
SUBROUTINE ls_read(IUNIT,vector,N)
use io, only: maxRecord
use precision
implicit none
Integer,intent(IN)      :: IUNIT,N
Real(realk),intent(OUT) :: vector(N)
!
Integer :: nBUF,startBUF,endBUF,I
nBUF     = n/maxRecord + 1
startBUF = 1
DO I=1,nBUF
  endBUF = MIN(n,startBUF+maxRecord-1)
  READ(IUNIT) vector(startBUF:endBUF)
  startBUF = startBUF + maxRecord
ENDDO
END SUBROUTINE ls_read

