MODULE memory_handling_nonp
!Monitor memory for integral code and possibly other parts!
   use memory_handling

!Interfaces for allocating/deallocating non-pointers, i.e. allocatables
!Use pointers whereever possible!
INTERFACE mem_alloc_nonp
  MODULE PROCEDURE real_allocate_1dim_nonp, real_allocate_2dim_nonp, &
     &             real_allocate_3dim_nonp,                          &
     &             real_allocate_5dim_nonp,                          &
     &             int_allocate_1dim_nonp,  int_allocate_2dim_nonp,  &
     &                                      int_allocate_4dim_nonp,  &
     &             logic_allocate_1dim_nonp
END INTERFACE

INTERFACE mem_dealloc_nonp
  MODULE PROCEDURE real_deallocate_1dim_nonp, real_deallocate_2dim_nonp, &
     &             real_deallocate_3dim_nonp,                            &  
     &             real_deallocate_5dim_nonp,                            &  
     &             int_deallocate_1dim_nonp,  int_deallocate_2dim_nonp,  &
     &                                        int_deallocate_4dim_nonp,  &
     &             logic_deallocate_1dim_nonp
END INTERFACE

CONTAINS

!----- ALLOCATE REAL ALLOCATABLES -----! 

SUBROUTINE real_allocate_1dim_nonp(A,n)
implicit none
integer,intent(in)         :: n
REAL(REALK),allocatable,intent(inout) :: A(:)
integer :: IERR
integer (kind=long) :: nsize
   ALLOCATE(A(n),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_1dim_nonp',IERR,n
     STOP
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_1dim_nonp

SUBROUTINE real_allocate_2dim_nonp(A,n1,n2)
implicit none
integer,intent(in)         :: n1, n2
REAL(REALK),allocatable,intent(inout) :: A(:,:)
integer :: IERR
integer (kind=long) :: nsize
   ALLOCATE(A(n1,n2),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_2dim_nonp',IERR,n1,n2
     STOP
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_2dim_nonp

SUBROUTINE real_allocate_3dim_nonp(A,n1,n2,n3)
implicit none
integer,intent(in)         :: n1, n2, n3
REAL(REALK),allocatable,intent(inout) :: A(:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   ALLOCATE(A(n1,n2,n3),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_3dim_nonp',IERR,n1,n2,n3
     STOP
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_3dim_nonp

SUBROUTINE real_allocate_5dim_nonp(A,n1,n2,n3,n4,n5)
implicit none
integer,intent(in)         :: n1, n2, n3, n4, n5
REAL(REALK),allocatable,intent(inout) :: A(:,:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   ALLOCATE(A(n1,n2,n3,n4,n5),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_5dim_nonp',IERR,n1,n2,n3,n4,n5
     STOP
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_5dim_nonp

!----- DEALLOCATE REAL ALLOCATABLES -----!

SUBROUTINE real_deallocate_1dim_nonp(A)
implicit none
REAL(REALK),allocatable,intent(inout) :: A(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_deallocate_1dim_nonp',IERR
     STOP
   ENDIF
END SUBROUTINE real_deallocate_1dim_nonp

SUBROUTINE real_deallocate_2dim_nonp(A)
implicit none
REAL(REALK),allocatable,intent(inout) :: A(:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_deallocate_2dim_nonp',IERR
     STOP
   ENDIF
END SUBROUTINE real_deallocate_2dim_nonp

SUBROUTINE real_deallocate_3dim_nonp(A)
implicit none
REAL(REALK),allocatable,intent(inout) :: A(:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_deallocate_3dim_nonp',IERR
     STOP
   ENDIF
END SUBROUTINE real_deallocate_3dim_nonp

SUBROUTINE real_deallocate_5dim_nonp(A)
implicit none
REAL(REALK),allocatable,intent(inout) :: A(:,:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_deallocate_5dim_nonp',IERR
     STOP
   ENDIF
END SUBROUTINE real_deallocate_5dim_nonp

!----- ALLOCATE INTEGER ALLOCATABLES -----!

SUBROUTINE int_allocate_1dim_nonp(I,n)
implicit none
integer,intent(in)         :: n
INTEGER,allocatable,intent(inout) :: I(:)
integer :: IERR
integer (kind=long) :: nsize
   ALLOCATE(I(n),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in int_allocate_1dim_nonp',IERR,n
     call lsquit('Error in int_allocate_1dim_nonp',-1)
   ENDIF
   nsize = size(I)*mem_intsize
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_1dim_nonp

SUBROUTINE int_allocate_2dim_nonp(I,n1,n2)
implicit none
integer,intent(in)         :: n1,n2
INTEGER,allocatable,intent(inout) :: I(:,:)
integer :: IERR
integer (kind=long) :: nsize
   ALLOCATE(I(n1,n2),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in int_allocate_2dim_nonp',IERR,n1,n2
     call lsquit('Error in int_allocate_2dim_nonp',-1)
   ENDIF
   nsize = size(I)*mem_intsize
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_2dim_nonp

SUBROUTINE int_allocate_4dim_nonp(I,n1,n2,n3,n4)
implicit none
integer,intent(in)         :: n1,n2,n3,n4
INTEGER,allocatable,intent(inout) :: I(:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   ALLOCATE(I(n1,n2,n3,n4),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in int_allocate_4dim_nonp',IERR,n1,n2,n3,n4
     call lsquit('Error in int_allocate_4dim_nonp',-1)
   ENDIF
   nsize = size(I)*mem_intsize
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_4dim_nonp

!----- DEALLOCATE INTEGER ALLOCATABLES -----!

SUBROUTINE int_deallocate_1dim_nonp(I)
implicit none
INTEGER,allocatable,intent(inout) :: I(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(I)*mem_intsize
   call mem_deallocated_mem_integer(nsize)
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE.0) THEN
      write(*,*) 'Error in int_deallocate_1dim_nonp',IERR
      STOP
   ENDIF
END SUBROUTINE int_deallocate_1dim_nonp

SUBROUTINE int_deallocate_2dim_nonp(I)
implicit none
INTEGER,allocatable,intent(inout) :: I(:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(I)*mem_intsize
   call mem_deallocated_mem_integer(nsize)
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE.0) THEN
      write(*,*) 'Error in int_deallocate_2dim_nonp',IERR
      STOP
   ENDIF
END SUBROUTINE int_deallocate_2dim_nonp

SUBROUTINE int_deallocate_4dim_nonp(I)
implicit none
INTEGER,allocatable,intent(inout) :: I(:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(I)*mem_intsize
   call mem_deallocated_mem_integer(nsize)
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE.0) THEN
      write(*,*) 'Error in int_deallocate_4dim_nonp',IERR
      STOP
   ENDIF
END SUBROUTINE int_deallocate_4dim_nonp

!----- ALLOCATE LOGICAL ALLOCATABLES -----!

SUBROUTINE logic_allocate_1dim_nonp(L,n)
implicit none
integer,intent(in)         :: n
LOGICAL,allocatable,intent(inout) :: L(:)
integer :: IERR
integer (kind=long) :: nsize
   ALLOCATE(L(n),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in logic_allocate_1dim_nonp',IERR,n
     STOP
   ENDIF
   nsize = size(L)*mem_logicalsize
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic_allocate_1dim_nonp

!----- DEALLOCATE LOGICAL ALLOCATABLES -----!

SUBROUTINE logic_deallocate_1dim_nonp(L)
implicit none
LOGICAL,allocatable,intent(inout) :: L(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(L)*mem_logicalsize
   call mem_deallocated_mem_logical(nsize)
   DEALLOCATE(L,STAT = IERR)
   IF (IERR.NE.0) THEN
      write(*,*) 'Error in logic_deallocate_1dim_nonp',IERR
      STOP
   ENDIF
END SUBROUTINE logic_deallocate_1dim_nonp

END MODULE memory_handling_nonp

