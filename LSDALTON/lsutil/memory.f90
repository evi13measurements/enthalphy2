MODULE memory_handling
!Monitor memory for integral code and possibly other parts!
   use precision
   integer(KIND=long),save :: mem_allocated_global = 0, max_mem_used_global = 0         !Count all memory
   integer(KIND=long),save :: mem_allocated_type_matrix = 0, max_mem_used_type_matrix = 0 !Count memory, density opt code
   integer(KIND=long),save :: mem_allocated_real = 0, max_mem_used_real = 0             !Count 'real' memory, integral code
   integer(KIND=long),save :: mem_allocated_integer = 0, max_mem_used_integer = 0       !Count 'integer' memory, integral code
   integer(KIND=long),save :: mem_allocated_character = 0, max_mem_used_character = 0   !Count 'character' memory, integral code
   integer(KIND=long),save :: mem_allocated_logical = 0, max_mem_used_logical = 0       !Count 'logical' memory, integral code
!Memory distributed on types:
   integer(KIND=long),save :: mem_allocated_linkshell = 0, max_mem_used_linkshell = 0         !Count memory, type linkshell
   integer(KIND=long),save :: mem_allocated_integralitem = 0, max_mem_used_integralitem = 0   !Count memory, type integralitem
   integer(KIND=long),save :: mem_allocated_integrand = 0, max_mem_used_integrand = 0         !Count memory, type integrand
   integer(KIND=long),save :: mem_allocated_overlap = 0, max_mem_used_overlap = 0             !Count memory, type Overlap
   integer(KIND=long),save :: mem_allocated_ODitem = 0, max_mem_used_ODitem = 0               !Count memory, type ODitem
   integer(KIND=long),save :: mem_allocated_DFT = 0, max_mem_used_DFT = 0  !Count memory in DFT
              !Count memory, type FMM
   integer(KIND=long),save :: mem_allocated_FMM = 0, max_mem_used_FMM = 0  
   integer(KIND=long),save :: mem_allocated_lstensor = 0, max_mem_used_lstensor = 0  
   integer(KIND=long),parameter :: mem_realsize=8
   integer(KIND=long),parameter :: mem_logicalsize=4
   integer(KIND=long),parameter :: mem_complexsize=16
   !integer(KIND=long),parameter :: mem_charsize=4 It is not obvious how we
   !should count character memory!
   !integer(KIND=long),parameter :: mem_pointersize=8
#if VAR_64int
   integer(KIND=long),parameter :: mem_intsize=8
#else
   integer(KIND=long),parameter :: mem_intsize=4
#endif
!Interfaces for allocating/deallocating pointers
INTERFACE mem_alloc
  MODULE PROCEDURE real_allocate_1dim, real_allocate_2dim, &
     &             real_allocate_2dim_zero,                & 
     &             real_allocate_3dim, real_allocate_4dim, &
     &             real_allocate_5dim, real_allocate_5dim_zero, &
     &             int_allocate_1dim,  int_allocate_2dim,  &
     &             char_allocate_1dim,                     &
     &             logic_allocate_1dim
END INTERFACE

INTERFACE mem_dealloc
  MODULE PROCEDURE real_deallocate_1dim, real_deallocate_2dim, &
     &             real_deallocate_3dim, real_deallocate_4dim, &
     &             real_deallocate_5dim,                       &
     &             int_deallocate_1dim,  int_deallocate_2dim,  &
     &             char_deallocate_1dim,                       &
     &             logic_deallocate_1dim
END INTERFACE

CONTAINS
!> \brief Print current and max. amount of memory allocated for different data types. 
!> \author S. Host
!> \date 2009
  subroutine stats_mem(lupri)
    implicit none
    !> Logical unit number for output file.
    integer,intent(in) :: lupri
    integer(KIND=long) :: tmp

    WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
    WRITE(LUPRI,'("                  Memory statistics          ")')
    WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
    WRITE(LUPRI,'("  Allocated memory (TOTAL):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_global
    WRITE(LUPRI,'("  Allocated memory (type(matrix)):  ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_type_matrix
    WRITE(LUPRI,'("  Allocated memory (real):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_real
    WRITE(LUPRI,'("  Allocated memory (integer):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_integer
    WRITE(LUPRI,'("  Allocated memory (logical):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_logical
    WRITE(LUPRI,'("  Allocated memory (linkshell):     ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_linkshell
    WRITE(LUPRI,'("  Allocated memory (integrand):     ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_integrand
    WRITE(LUPRI,'("  Allocated memory (integralitem):  ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_integralitem
    WRITE(LUPRI,'("  Allocated memory (overlap):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_overlap
    WRITE(LUPRI,'("  Allocated memory (ODitem):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ODitem
    WRITE(LUPRI,'("  Allocated memory (lstensor):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_lstensor
    WRITE(LUPRI,'("  Allocated memory (FMM   ):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_FMM
    WRITE(LUPRI,'("  Allocated memory (DFT   ):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_DFT

    call print_maxmem(lupri,max_mem_used_global,'TOTAL')
    call print_maxmem(lupri,max_mem_used_type_matrix,'type(matrix)')
    CALL print_maxmem(lupri,max_mem_used_real,'real(realk)')
    CALL print_maxmem(lupri,max_mem_used_integer,'integer')
    CALL print_maxmem(lupri,max_mem_used_logical,'logical')
    CALL print_maxmem(lupri,max_mem_used_linkshell,'linkshell')
    CALL print_maxmem(lupri,max_mem_used_integrand,'integrand')
    CALL print_maxmem(lupri,max_mem_used_integralitem,'integralitem')
    CALL print_maxmem(lupri,max_mem_used_overlap,'Overlap')
    CALL print_maxmem(lupri,max_mem_used_ODitem,'ODitem')
    CALL print_maxmem(lupri,max_mem_used_lstensor,'LStensor')
    CALL print_maxmem(lupri,max_mem_used_FMM,'FMM')
    CALL print_maxmem(lupri,max_mem_used_DFT,'DFT')

    WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
    WRITE(LUPRI,*)
  end subroutine stats_mem

   !> \brief For given SCF it, print amount of memory allocated for different data types. 
   !> \author S. Host
   !> \date 2009
  subroutine scf_stats_debug_mem(lupri,it)
    implicit none
    !> Current SCF iteration
    integer, intent(in) :: it
    !> Logical unit number for output file
    integer,intent(in) :: lupri

    WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
    WRITE(LUPRI,'("                  Memory statistics, iteration", i4)') it
    WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
    call print_maxmem(lupri,max_mem_used_global,'TOTAL')
    if (max_mem_used_type_matrix > 1.0d0) call print_maxmem(lupri,max_mem_used_type_matrix,'type(matrix)')
    if (max_mem_used_real > 1.0d0) call print_maxmem(lupri,max_mem_used_real,'real')
    if (max_mem_used_integer > 1.0d0) call print_maxmem(lupri,max_mem_used_integer,'integer')
    if (max_mem_used_logical > 1.0d0) call print_maxmem(lupri,max_mem_used_logical,'logical')
    if (max_mem_used_linkshell > 1.0d0) call print_maxmem(lupri,max_mem_used_linkshell,'linkshell')
    if (max_mem_used_integrand > 1.0d0) call print_maxmem(lupri,max_mem_used_integrand,'integrand')
    if (max_mem_used_integralitem > 1.0d0) call print_maxmem(lupri,max_mem_used_integralitem,'integralitem')
    if (max_mem_used_overlap > 1.0d0) call print_maxmem(lupri,max_mem_used_overlap,'overlap')
    if (max_mem_used_ODitem > 1.0d0) call print_maxmem(lupri,max_mem_used_ODitem,'ODitem')
    if (max_mem_used_lstensor > 1.0d0) call print_maxmem(lupri,max_mem_used_lstensor,'lstensor')
    if (max_mem_used_FMM > 1.0d0) call print_maxmem(lupri,max_mem_used_FMM,'FMM    ')
    if (max_mem_used_DFT > 1.0d0) call print_maxmem(lupri,max_mem_used_DFT,'DFT    ')
    WRITE(LUPRI,*)
    call print_mem_alloc(lupri,mem_allocated_global,'TOTAL')
    if (mem_allocated_type_matrix > 1.0d0) call print_mem_alloc(lupri,mem_allocated_type_matrix,'type(matrix)')
    if (mem_allocated_real > 1.0d0) call print_mem_alloc(lupri,mem_allocated_real,'real')
    if (mem_allocated_integer > 1.0d0) call print_mem_alloc(lupri,mem_allocated_integer,'integer')
    if (mem_allocated_logical > 1.0d0) call print_mem_alloc(lupri,mem_allocated_logical,'logical')
    if (mem_allocated_linkshell > 1.0d0) call print_mem_alloc(lupri,mem_allocated_linkshell,'linkshell')
    if (mem_allocated_integrand > 1.0d0) call print_mem_alloc(lupri,mem_allocated_integrand,'integrand')
    if (mem_allocated_integralitem > 1.0d0) call print_mem_alloc(lupri,mem_allocated_integralitem,'integralitem')
    if (mem_allocated_overlap > 1.0d0) call print_mem_alloc(lupri,mem_allocated_overlap,'overlap')
    if (mem_allocated_ODitem > 1.0d0) call print_mem_alloc(lupri,mem_allocated_ODitem,'ODitem')
    if (mem_allocated_lstensor > 1.0d0) call print_mem_alloc(lupri,mem_allocated_lstensor,'lstensor')
    if (mem_allocated_FMM > 1.0d0) call print_mem_alloc(lupri,mem_allocated_FMM,'FMM   ')
    if (mem_allocated_DFT > 1.0d0) call print_mem_alloc(lupri,mem_allocated_DFT,'DFT   ')
    WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
    WRITE(LUPRI,*)
  end subroutine scf_stats_debug_mem

   !> \brief status information printout. Print amount of memory allocated for different data types. 
   !> \author T. Kjaergaard
   !> \date 2009
  subroutine debug_mem_stats(lupri)
    implicit none
    !> Logical unit number for output file
    integer,intent(in) :: lupri

    WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
    WRITE(LUPRI,'("                  Debug Memory Statistics              ")')
    WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
    call print_maxmem(lupri,max_mem_used_global,'TOTAL')
    if (max_mem_used_type_matrix > 1.0d0) call print_maxmem(lupri,max_mem_used_type_matrix,'type(matrix)')
    if (max_mem_used_real > 1.0d0) call print_maxmem(lupri,max_mem_used_real,'real')
    if (max_mem_used_integer > 1.0d0) call print_maxmem(lupri,max_mem_used_integer,'integer')
    if (max_mem_used_logical > 1.0d0) call print_maxmem(lupri,max_mem_used_logical,'logical')
    if (max_mem_used_linkshell > 1.0d0) call print_maxmem(lupri,max_mem_used_linkshell,'linkshell')
    if (max_mem_used_integrand > 1.0d0) call print_maxmem(lupri,max_mem_used_integrand,'integrand')
    if (max_mem_used_integralitem > 1.0d0) call print_maxmem(lupri,max_mem_used_integralitem,'integralitem')
    if (max_mem_used_overlap > 1.0d0) call print_maxmem(lupri,max_mem_used_overlap,'overlap')
    if (max_mem_used_ODitem > 1.0d0) call print_maxmem(lupri,max_mem_used_ODitem,'ODitem')
    if (max_mem_used_lstensor > 1.0d0) call print_maxmem(lupri,max_mem_used_lstensor,'lstensor')
    if (max_mem_used_FMM > 1.0d0) call print_maxmem(lupri,max_mem_used_FMM,'FMM    ')
    if (max_mem_used_DFT > 1.0d0) call print_maxmem(lupri,max_mem_used_DFT,'DFT    ')
    WRITE(LUPRI,*)
    call print_mem_alloc(lupri,mem_allocated_global,'TOTAL')
    if (mem_allocated_type_matrix > 1.0d0) call print_mem_alloc(lupri,mem_allocated_type_matrix,'type(matrix)')
    if (mem_allocated_real > 1.0d0) call print_mem_alloc(lupri,mem_allocated_real,'real')
    if (mem_allocated_integer > 1.0d0) call print_mem_alloc(lupri,mem_allocated_integer,'integer')
    if (mem_allocated_logical > 1.0d0) call print_mem_alloc(lupri,mem_allocated_logical,'logical')
    if (mem_allocated_linkshell > 1.0d0) call print_mem_alloc(lupri,mem_allocated_linkshell,'linkshell')
    if (mem_allocated_integrand > 1.0d0) call print_mem_alloc(lupri,mem_allocated_integrand,'integrand')
    if (mem_allocated_integralitem > 1.0d0) call print_mem_alloc(lupri,mem_allocated_integralitem,'integralitem')
    if (mem_allocated_overlap > 1.0d0) call print_mem_alloc(lupri,mem_allocated_overlap,'overlap')
    if (mem_allocated_ODitem > 1.0d0) call print_mem_alloc(lupri,mem_allocated_ODitem,'ODitem')
    if (mem_allocated_lstensor > 1.0d0) call print_mem_alloc(lupri,mem_allocated_lstensor,'lstensor')
    if (mem_allocated_FMM > 1.0d0) call print_mem_alloc(lupri,mem_allocated_FMM,'FMM   ')
    if (mem_allocated_DFT > 1.0d0) call print_mem_alloc(lupri,mem_allocated_DFT,'DFT   ')
    WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
    WRITE(LUPRI,*)

  end subroutine debug_mem_stats

   !> \brief Print routine for memory statistics (max. allocated memory). 
   !> \author S. Host
   !> \date 2009
  subroutine print_maxmem(lupri,max_mem_used,STRING)
    implicit none
    !> Logical unit number for output file
    integer :: lupri
    !> Amount of memory used
    integer(KIND=long) :: max_mem_used
    !> Name of data type for which memory usage is printed
    Character*(*)    ::  STRING
    character(len=20) :: trimstring
    !Character(len=3) ::  FORMATSTRING1
    !Character(len=23) ::  FORMATSTRING2

    trimstring = adjustl(string)
    if (max_mem_used < 100) then !Divide by 1 to typecast real
       WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " Byte")') trimstring, max_mem_used/1.0d0
    else if (max_mem_used < 1000000) then
       WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " kB")') trimstring, max_mem_used/1024.0d0
    else if (max_mem_used < 1000000000) then
       WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " MB")') trimstring, max_mem_used/(1024.0d0*1024.0d0)
    else
       WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " GB")') trimstring, max_mem_used/(1024.0d0*1024.0d0*1024.0d0)
    endif
  end subroutine print_maxmem

   !> \brief Print routine for memory statistics (current allocated memory). 
   !> \author S. Host
   !> \date 2009
  subroutine print_mem_alloc(lupri,mem_used,STRING)
    implicit none
    !> Logical unit number for output file
    integer :: lupri
    !> Amount of memory used
    integer(KIND=long) :: mem_used
    !> Name of data type for which memory usage is printed
    Character*(*)    ::  STRING
    character(len=20) :: trimstring

    trimstring = adjustl(string)
    if (mem_used < 100) then !Divide by 1 to typecast real
       WRITE(LUPRI,'("  Allocated memory, ", a24, f9.3, " Byte")') trimstring, mem_used/1.0d0
    else if (mem_used < 1000000) then
       WRITE(LUPRI,'("  Allocated memory, ", a24, f9.3, " kB")') trimstring, mem_used/1024.0d0
    else if (mem_used < 1000000000) then
       WRITE(LUPRI,'("  Allocated memory, ", a24, f9.3, " MB")') trimstring, mem_used/(1024.0d0*1024.0d0)
    else
       WRITE(LUPRI,'("  Allocated memory, ", a24, f9.3, " GB")') trimstring, mem_used/(1024.0d0*1024.0d0*1024.0d0)
    endif
  end subroutine print_mem_alloc

!SUBROUTINE mem_print(iunit)
!implicit none
!integer,intent(in) :: iunit
!  call print_mem('Global',mem_allocated_global,iunit)
!END SUBROUTINE mem_print
!
!SUBROUTINE mem_max_print(iunit)
!implicit none
!integer,intent(in) :: iunit
!  call print_mem('Mamimum',max_mem_used_global,iunit)
!END SUBROUTINE mem_max_print
!
!SUBROUTINE print_mem(id,mem,iunit)
!implicit none
!character*(*),intent(in) :: id
!integer,intent(in) :: iunit
!integer(kind=long), intent(in) :: mem
!!
!character(80) :: memory_string
!
!IF (mem.GT.(1024**3)) THEN
!  WRITE(memory_string,'(1X,I10,A1,I3,1X,A6)') int(mem/1024.d0**3),'.', &
!     &   int((mem-int(mem/1024.d0**3)*1024.d0**3)/1024.d0**3*1000),'Gbytes'
!ELSE IF (mem.GT.(1024**2)) THEN
!  WRITE(memory_string,'(1X,I10,A1,I3,1X,A6)') int(mem/1024.d0**2),'.', &
!     &   int((mem-int(mem/1024.d0**2)*1024.d0**2)/1024.d0**2*1000),'Mbytes'
!ELSE IF (mem.GT.(1024)) THEN
!  WRITE(memory_string,'(1X,I10,A1,I3,1X,A6)') int(mem/1024.d0),'.', &
!     &   int((mem-int(mem/1024.d0)*1024.d0)/1024.d0*1000),'kbytes'
!ELSE
!  WRITE(memory_string,'(1X,I14,1X,A6)') mem,' bytes'
!ENDIF
!WRITE(IUNIT,'(3X,A,X,A)') id,memory_string
!WRITE(*,'(3X,A10,A14,A40)') id,' memory usage:',memory_string
!END SUBROUTINE print_mem

!----- ALLOCATE REAL POINTERS -----! 

SUBROUTINE real_allocate_1dim(A,n)
implicit none
integer,intent(in)  :: n
REAL(REALK),pointer :: A(:)
integer :: IERR
integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_1dim',IERR,n
     write(*,'(1X,A)') IERR/0
!    CALL lsQUIT('Error in real_allocate_1dim')
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_1dim

SUBROUTINE real_allocate_2dim(A,n1,n2)
implicit none
integer,intent(in)  :: n1, n2
REAL(REALK),pointer :: A(:,:)
integer :: IERR
integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_2dim',IERR,n1,n2
     STOP
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_2dim

SUBROUTINE real_allocate_2dim_zero(A,n1,n2,First,Second)
! Allocates 2d arrays starting from zero index
! for first,second or both dimensions
implicit none
integer,intent(in)  :: n1, n2
REAL(REALK),pointer :: A(:,:)
integer :: IERR
integer (kind=long) :: nsize
Logical :: First, Second
!
   nullify(A)
! Both need zeroth element
   If (First .AND. Second) then
      ALLOCATE(A(0:n1,0:n2),STAT = IERR)
      IF (IERR.NE.0) THEN
          write(*,*) 'Error in real_allocate_2dim_zero',IERR,n1,n2
          STOP
      ENDIF
   Else 
! Only one 
        If (First) then
           ALLOCATE(A(0:n1,n2),STAT = IERR)
           IF (IERR.NE.0) THEN
               write(*,*) 'Error in real_allocate_2dim_zero',IERR,n1,n2
               STOP
           ENDIF
        Else 
           If (Second) then
              ALLOCATE(A(n1,0:n2),STAT = IERR)
              IF (IERR.NE.0) THEN
                  write(*,*) 'Error in real_allocate_2dim_zero',IERR,n1,n2
                  STOP
              ENDIF
           Else
! None :: an error, should be at least one.
              write(*,*) 'Error in real_allocate_2dim_zero'
              STOP
           Endif
        Endif
   Endif 
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_2dim_zero

SUBROUTINE real_allocate_3dim(A,n1,n2,n3)
implicit none
integer,intent(in)  :: n1, n2, n3
REAL(REALK),pointer :: A(:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2,n3),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_3dim',IERR,n1,n2,n3
     STOP
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_3dim

SUBROUTINE real_allocate_4dim(A,n1,n2,n3,n4)
implicit none
integer,intent(in)  :: n1,n2,n3,n4
REAL(REALK),pointer :: A(:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2,n3,n4),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_4dim',IERR,n1,n2,n3,n4
     STOP
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_4dim

SUBROUTINE real_allocate_5dim(A,n1,n2,n3,n4,n5)
implicit none
integer,intent(in)  :: n1,n2,n3,n4,n5
REAL(REALK),pointer :: A(:,:,:,:,:)
integer :: IERR
integer (kind=long) :: itest, nsize
   nullify(A)
   ALLOCATE(A(n1,n2,n3,n4,n5),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_5dim',IERR,n1,n2,n3,n4,n5
     STOP
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_5dim

SUBROUTINE real_allocate_5dim_zero(A,n1,n2,n3,n4,n5,z1,z2,z3,z4,z5)
! Allocates 5d arrays starting from zero index
! for some or all dimensions
implicit none
integer,intent(in)  :: n1,n2,n3,n4,n5
logical,intent(in)  :: z1,z2,z3,z4,z5
REAL(REALK),pointer :: A(:,:,:,:,:)
integer :: IERR
integer (kind=long) :: itest, nsize
integer             :: i1,i2,i3,i4,i5
   nullify(A)
   i1=1
   IF(z1)i1=0
   i2=1
   IF(z2)i2=0
   i3=1
   IF(z3)i3=0
   i4=1
   IF(z4)i4=0
   i5=1
   IF(z5)i5=0
   ALLOCATE(A(i1:n1,i2:n2,i3:n3,i4:n4,i5:n5),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_allocate_5dim_zero',IERR,n1,n2,n3,n4,n5
     STOP
   ENDIF
   nsize = size(A)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_5dim_zero

!----- DEALLOCATE REAL POINTERS -----!

SUBROUTINE real_deallocate_1dim(A)
implicit none
REAL(REALK),pointer :: A(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call lsQUIT('Error in real_deallocate_1dim - memory previously released',-1)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_deallocate_1dim',IERR
     STOP
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_1dim

SUBROUTINE real_deallocate_2dim(A)
implicit none
REAL(REALK),pointer :: A(:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call lsQUIT('Error in real_deallocate_2dim - memory previously released',-1)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_deallocate_2dim',IERR
     STOP
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_2dim

SUBROUTINE real_deallocate_3dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call lsQUIT('Error in real_deallocate_3dim - memory previously released',-1)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_deallocate_3dim',IERR
     STOP
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_3dim

SUBROUTINE real_deallocate_4dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call lsQUIT('Error in real_deallocate_4dim - memory previously released',-1)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_deallocate_4dim',IERR
     STOP
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_4dim

SUBROUTINE real_deallocate_5dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call lsQUIT('Error in real_deallocate_5dim - memory previously released',-1)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in real_deallocate_5dim',IERR
     STOP
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_5dim

!----- ALLOCATE INTEGER POINTERS -----!

SUBROUTINE int_allocate_1dim(I,n)
implicit none
integer,intent(in)  :: n
INTEGER,pointer     :: I(:)
integer :: IERR
integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in int_allocate_1dim',IERR,n
     call lsquit('Error in int_allocate_1dim',-1)
   ENDIF
   nsize = size(I)*mem_intsize
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_1dim

SUBROUTINE int_allocate_2dim(I,n1,n2)
implicit none
integer,intent(in) :: n1,n2
INTEGER,pointer    :: I(:,:)
integer :: IERR
integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n1,n2),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in int_allocate_2dim',IERR,n1,n2
     call lsquit('Error in int_allocate_2dim',-1)
   ENDIF
   nsize = size(I)*mem_intsize
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_2dim

!----- DEALLOCATE INTEGER POINTERS -----!

SUBROUTINE int_deallocate_1dim(I)
implicit none
INTEGER,pointer :: I(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(I)*mem_intsize
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call lsQUIT('Error in int_deallocate_1dim - memory previously released',-1)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE.0) THEN
      write(*,*) 'Error in int_deallocate_1dim',IERR
      STOP
   ENDIF
   nullify(I)
END SUBROUTINE int_deallocate_1dim

SUBROUTINE int_deallocate_2dim(I)
implicit none
INTEGER,pointer :: I(:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(I)*mem_intsize
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call lsQUIT('Error in int_deallocate_2dim - memory previously released',-1)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE.0) THEN
      write(*,*) 'Error in int_deallocate_2dim',IERR
      STOP
   ENDIF
   nullify(I)
END SUBROUTINE int_deallocate_2dim

!----- ALLOCATE CHARACTER POINTERS -----!

SUBROUTINE char_allocate_1dim(C,n)
implicit none
integer,intent(in)         :: n
CHARACTER(LEN=*),pointer :: C(:)
integer :: IERR
integer (kind=long) :: nsize
   nullify(C)
   ALLOCATE(C(n),STAT = IERR)
   IF (IERR.NE.0) THEN
      write(*,*) 'Error in char_allocate_1dim',IERR,n
      STOP
   ENDIF
   nsize = mem_complexsize*size(C)
   call mem_allocated_mem_character(nsize)
END SUBROUTINE char_allocate_1dim

!----- DEALLOCATE CHARACTER POINTERS -----!

SUBROUTINE char_deallocate_1dim(C)
implicit none
CHARACTER(LEN=*),pointer :: C(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = mem_complexsize*size(C)
   call mem_deallocated_mem_character(nsize)
   if (.not.ASSOCIATED(C)) then
      print *,'Memory previously released!!'
      call lsQUIT('Error in char_deallocate_1dim - memory previously released',-1)
   endif
   DEALLOCATE(C,STAT = IERR)
   IF (IERR.NE.0) THEN
      write(*,*) 'Error in char_deallocate_1dim',IERR
      STOP
   ENDIF
END SUBROUTINE char_deallocate_1dim

!----- ALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic_allocate_1dim(L,n)
implicit none
integer,intent(in) :: n
LOGICAL,pointer    :: L(:)
integer :: IERR
integer (kind=long) :: nsize
   nullify(L)
   ALLOCATE(L(n),STAT = IERR)
   IF (IERR.NE.0) THEN
     write(*,*) 'Error in logic_allocate_1dim',IERR,n
     STOP
   ENDIF
   nsize = size(L)*mem_logicalsize 
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic_allocate_1dim

!----- DEALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic_deallocate_1dim(L)
implicit none
LOGICAL,pointer :: L(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(L)*mem_logicalsize
   call mem_deallocated_mem_logical(nsize)
   if (.not.ASSOCIATED(L)) then
      print *,'Memory previously released!!'
      call lsQUIT('Error in logic_deallocate_1dim - memory previously released',-1)
   endif
   DEALLOCATE(L,STAT = IERR)
   IF (IERR.NE.0) THEN
      write(*,*) 'Error in logic_deallocate_1dim',IERR
      STOP
   ENDIF
   NULLIFY(L)
END SUBROUTINE logic_deallocate_1dim

!----- MEMORY HANDLING -----!

  subroutine mem_allocated_mem_real(nsize)
     implicit none
     integer (kind=long), intent(in) :: nsize
!$OMP CRITICAL (memory)
     mem_allocated_real = mem_allocated_real + nsize
     if (mem_allocated_real < 0) then
        write(*,*) 'Real memory negative! mem_allocated_real =', mem_allocated_real
        write(*,*) 'Real memory negative! nsize =', nsize
        call lsQUIT('Error in mem_allocated_mem_real - probably integer overflow!',-1)
     endif
     max_mem_used_real = MAX(max_mem_used_real,mem_allocated_real)
     !Count also the total memory:
     mem_allocated_global = mem_allocated_global  + nsize
     if (mem_allocated_global < 0) then
        write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
        call lsQUIT('Error in mem_allocated_mem_real - probably integer overflow!',-1)
     endif
     max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
!$OMP END CRITICAL (memory)
  end subroutine mem_allocated_mem_real

  subroutine mem_deallocated_mem_real(nsize)
     implicit none
     integer (kind=long), intent(in) :: nsize
!$OMP CRITICAL (memory)
     mem_allocated_real = mem_allocated_real - nsize
     if (mem_allocated_real < 0) then
        write(*,*) 'Real memory negative! mem_allocated_real =', mem_allocated_real
        call lsQUIT('Error in mem_deallocated_mem_real - something wrong with deallocation!',-1)
     endif
     !Count also the total memory:
     mem_allocated_global = mem_allocated_global - nsize
     if (mem_allocated_global < 0) then
        write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
        call lsQUIT('Error in mem_deallocated_mem_real - something wrong with deallocation!',-1)
     endif
!$OMP END CRITICAL (memory)
  end subroutine mem_deallocated_mem_real

  subroutine mem_allocated_mem_integer(nsize)
     implicit none
     integer(kind=long), intent(in) :: nsize
!$OMP CRITICAL (memory)
     mem_allocated_integer = mem_allocated_integer + nsize
     max_mem_used_integer = MAX(max_mem_used_integer,mem_allocated_integer)
     !Count also the total memory:
     mem_allocated_global = mem_allocated_global  + nsize
     max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
!$OMP END CRITICAL (memory)
  end subroutine mem_allocated_mem_integer

  subroutine mem_deallocated_mem_integer(nsize)
     implicit none
     integer (kind=long), intent(in) :: nsize
!$OMP CRITICAL (memory)
     mem_allocated_integer = mem_allocated_integer - nsize
     if (mem_allocated_integer < 0) then
        call lsQUIT('Error in mem_deallocated_mem_integer - probably integer overflow!',-1)
     endif
     !Count also the total memory:
     mem_allocated_global = mem_allocated_global - nsize
     if (mem_allocated_global < 0) then
        call lsQUIT('Error in mem_deallocated_mem_integer - probably integer overflow!',-1)
     endif
!$OMP END CRITICAL (memory)
  end subroutine mem_deallocated_mem_integer

  subroutine mem_allocated_mem_character(nsize)
     implicit none
     integer (kind=long), intent(in) :: nsize
!$OMP CRITICAL (memory)
     mem_allocated_character = mem_allocated_character + nsize
     max_mem_used_character = MAX(max_mem_used_character,mem_allocated_character)
     !Count also the total memory:
     mem_allocated_global = mem_allocated_global  + nsize
     max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
!$OMP END CRITICAL (memory)
  end subroutine mem_allocated_mem_character

  subroutine mem_deallocated_mem_character(nsize)
     implicit none
     integer (kind=long), intent(in) :: nsize
!$OMP CRITICAL (memory)
     mem_allocated_character = mem_allocated_character - nsize
     if (mem_allocated_character < 0) then
        call lsQUIT('Error in mem_deallocated_mem_character - probably integer overflow!',-1)
     endif
     !Count also the total memory:
     mem_allocated_global = mem_allocated_global - nsize
     if (mem_allocated_global < 0) then
        call lsQUIT('Error in mem_deallocated_mem_character - probably integer overflow!',-1)
     endif
!$OMP END CRITICAL (memory)
  end subroutine mem_deallocated_mem_character

  subroutine mem_allocated_mem_logical(nsize)
     implicit none
     integer (kind=long), intent(in) :: nsize
!$OMP CRITICAL (memory)
     mem_allocated_logical = mem_allocated_logical + nsize
     max_mem_used_logical = MAX(max_mem_used_logical,mem_allocated_logical)
     !Count also the total memory:
     mem_allocated_global = mem_allocated_global  + nsize
     max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
!$OMP END CRITICAL (memory)
  end subroutine mem_allocated_mem_logical

  subroutine mem_deallocated_mem_logical(nsize)
     implicit none
     integer (kind=long), intent(in) :: nsize
!$OMP CRITICAL (memory)
     mem_allocated_logical = mem_allocated_logical - nsize
     if (mem_allocated_logical < 0) then
        call lsQUIT('Error in mem_deallocated_mem_logical - probably integer overflow!',-1)
     endif
     !Count also the total memory:
     mem_allocated_global = mem_allocated_global - nsize
     if (mem_allocated_global < 0) then
        call lsQUIT('Error in mem_deallocated_mem_logical - probably integer overflow!',-1)
     endif
!$OMP END CRITICAL (memory)
  end subroutine mem_deallocated_mem_logical
END MODULE memory_handling

