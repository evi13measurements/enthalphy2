!> @file 
!> Contains lstensor structure and associated subroutines
MODULE LSTENSOR_OPERATIONSMOD
  use precision
  use basis_type
  use Matrix_module
  use matrix_operations
  use memory_handling
  use OD_Type
!  INTERFACE lstensor_init   
!     MODULE PROCEDURE lstensor_init_1dim, lstensor_init_2dim, &
!          &           lstensor_init_3dim, lstensor_init_4dim
!     MODULE PROCEDURE lstensor_init_4dim
!  END INTERFACE
  
  INTERFACE Build_lstensor_from_dense_mat
     MODULE PROCEDURE Build_lstensor_from_dense_mat_single, &
          &           Build_lstensor_from_dense_mat_array
  END INTERFACE

!  INTERFACE Build_dense_mat_from_lstensor
!     MODULE PROCEDURE Build_dense_mat_from_lstensor_single, &
!          &           Build_dense_mat_from_lstensor_array
!  END INTERFACE


  INTERFACE Build_matrix_from_lstensor
     MODULE PROCEDURE Build_singlematrix_from_lstensor, &
          &           Build_matrixarray_from_lstensor
  END INTERFACE

TYPE LSBATCHTENSOR
INTEGER                 :: nAngmomA
INTEGER                 :: nAngmomB
INTEGER                 :: nAngmomC
INTEGER                 :: nAngmomD
INTEGER                 :: nContA(maxAOangmom)
INTEGER                 :: nOrbCompA(maxAOangmom)
INTEGER                 :: startOrbA(maxAOangmom)
INTEGER                 :: nContB(maxAOangmom)
INTEGER                 :: nOrbCompB(maxAOangmom)
INTEGER                 :: startOrbB(maxAOangmom)
INTEGER                 :: nContC(maxAOangmom)
INTEGER                 :: nOrbCompC(maxAOangmom)
INTEGER                 :: startOrbC(maxAOangmom)
INTEGER                 :: nContD(maxAOangmom)
INTEGER                 :: nOrbCompD(maxAOangmom)
INTEGER                 :: startOrbD(maxAOangmom)
INTEGER                 :: nelmE
REAL(REALK),pointer     :: elms(:) !nContA*nOrbCompA*...nContD*nOrbCompD*nmat
INTEGER                 :: nelms
END TYPE LSBATCHTENSOR

TYPE LSAOTENSOR
TYPE(LSBATCHTENSOR),pointer :: BATCH(:,:,:,:)
INTEGER :: ATOM1
INTEGER :: ATOM2
INTEGER :: ATOM3
INTEGER :: ATOM4
INTEGER :: nAOBATCH1
INTEGER :: nAOBATCH2
INTEGER :: nAOBATCH3
INTEGER :: nAOBATCH4
LOGICAL :: ALLOC
END TYPE LSAOTENSOR

TYPE LSTENSOR
TYPE(LSAOTENSOR),pointer :: LSAO(:)
INTEGER :: nAtom1
INTEGER :: nAtom2
INTEGER :: nAtom3
INTEGER :: nAtom4
INTEGER :: nmat
INTEGER :: nbast1
INTEGER :: nbast2
INTEGER :: nbast3
INTEGER :: nbast4
INTEGER :: nLSAO
INTEGER,pointer :: INDEX(:,:,:,:) !GIVES THE INDEX IN LSAO (or zero)
logical :: gradienttensor
logical :: primCStensor
END TYPE LSTENSOR

CONTAINS
!> \brief determines the memory requirement for lstensor structure
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor 
!> \param LUDISK the logical unit number of output file
SUBROUTINE Determine_lstensor_memory(tensor,mem_alloc2)
implicit none
type(lstensor),intent(in) :: tensor
integer(kind=long),intent(out) :: mem_alloc2
integer(kind=long) :: sizebatchtensor,sizelsaotensor
integer(kind=long) :: sizelstensor
integer :: I,Ibat,Jbat,Kbat,Lbat

sizebatchtensor = maxAOangmom*mem_intsize*12+6*mem_intsize
sizelsaotensor = 8*mem_intsize+mem_logicalsize
sizelstensor =  10*mem_intsize+2*mem_logicalsize+size(tensor%INDEX)*mem_intsize
mem_alloc2 = sizelstensor
do I=1,TENSOR%nLSAO
   mem_alloc2 = mem_alloc2 + sizelsaotensor
   IF(TENSOR%LSAO(I)%ALLOC)THEN
      DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
         DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
            DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
               DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
                  mem_alloc2=mem_alloc2+sizebatchtensor
                  mem_alloc2=mem_alloc2+TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms*mem_realsize
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF
enddo
!if (mem_alloc2 < 100) then !Divide by 1 to typecast real
!   WRITE(6,'(A,A3,A,f9.3,A)')' Mem required for lstensor:',mem_alloc2/1.0d0,' Byte'
!else if (mem_alloc2 < 1000000) then
!   WRITE(6,'(A,A3,A,f9.3,A)')' Mem required for lstensor:',mem_alloc2/1024.0d0,' kB'
!else if (mem_alloc2 < 1000000000) then
!   WRITE(6,'(A,A3,A,f9.3,A)')' Mem required for lstensor:',mem_alloc2/(1024.0d0*1024.0d0),' MB'
!else
!   WRITE(6,'(A,A3,A,f9.3,A)')' Mem required for lstensor:',mem_alloc2/(1024.0d0*1024.0d0*1024.0d0),' GB'
!endif

END SUBROUTINE DETERMINE_LSTENSOR_MEMORY

!> \brief add memory used for lstensor to global statistics
!> \author T. Kjaergaard
!> \date 2010
!> \param nsize size of memory
SUBROUTINE ADD_MEM_TO_GLOBAL(nsize)
implicit none
integer(KIND=long),intent(in) :: nsize

!print*,'ADD_MEM_TO_GLOBAL1:                       SIZE',nsize
!print*,'ADD_MEM_TO_GLOBAL1: the mem_allocated_global  ',mem_allocated_global
!print*,'ADD_MEM_TO_GLOBAL1: the mem_allocated_lstensor',mem_allocated_lstensor
mem_allocated_lstensor = mem_allocated_lstensor + nsize 
max_mem_used_lstensor = MAX(mem_allocated_lstensor,max_mem_used_lstensor)
mem_allocated_global = mem_allocated_global + nsize  
max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
!print*,'ADD_MEM_TO_GLOBAL2: the mem_allocated_global  ',mem_allocated_global
!print*,'ADD_MEM_TO_GLOBAL2: the mem_allocated_lstensor',mem_allocated_lstensor

END SUBROUTINE ADD_MEM_TO_GLOBAL

!> \brief remove memory used for lstensor from global statistics
!> \author T. Kjaergaard
!> \date 2010
!> \param nsize size of memory
SUBROUTINE REMOVE_MEM_FROM_GLOBAL(nsize)
implicit none
integer(KIND=long),intent(in) :: nsize

!print*,'REMOVE_MEM_FROM_GLOBAL1:                       SIZE',nsize
!print*,'REMOVE_MEM_FROM_GLOBAL1: the mem_allocated_global  ',mem_allocated_global
!print*,'REMOVE_MEM_FROM_GLOBAL1: the mem_allocated_lstensor',mem_allocated_lstensor
mem_allocated_lstensor = mem_allocated_lstensor - nsize 
mem_allocated_global = mem_allocated_global - nsize  
!print*,'REMOVE_MEM_FROM_GLOBAL2: the mem_allocated_global  ',mem_allocated_global
!print*,'REMOVE_MEM_FROM_GLOBAL2: the mem_allocated_lstensor',mem_allocated_lstensor

END SUBROUTINE REMOVE_MEM_FROM_GLOBAL

!> \brief write the lstensor structure to disk
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor to be written to disk
!> \param LUDISK the logical unit number of file to write to
SUBROUTINE Write_lstensor_to_disk(TENSOR,LUDISK,lupri)
implicit none
type(lstensor) :: tensor
integer :: ludisk,lupri
!
integer :: iatom,j,imat2,i,ibat,jbat,np1,np2,ip
integer :: kbat,lbat,iangmom,jangmom,kangmom,langmom
integer :: Ielm,IORB,JORB,KORB,LORB,nContA,nOrbCompA,nmat
integer :: nContB,nOrbCompB,nContC,nOrbCompC,nContD,nOrbCompD
integer(kind=long) :: nmemsize

WRITE(LUDISK)TENSOR%primCStensor 
WRITE(LUDISK)TENSOR%gradienttensor 
WRITE(LUDISK)TENSOR%natom1
WRITE(LUDISK)TENSOR%natom2
WRITE(LUDISK)TENSOR%natom3
WRITE(LUDISK)TENSOR%natom4
WRITE(LUDISK)TENSOR%nbast1
WRITE(LUDISK)TENSOR%nbast2
WRITE(LUDISK)TENSOR%nbast3
WRITE(LUDISK)TENSOR%nbast4
WRITE(LUDISK)TENSOR%nmat
WRITE(LUDISK)TENSOR%nLSAO
nmat = TENSOR%nmat
IF(TENSOR%gradienttensor)THEN
   DO Iatom = 1,TENSOR%natom1
      DO J=1,3
         DO IMAT2 = 1,nmat
            WRITE(LUDISK) TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms(J+(IMAT2-1)*3)
         ENDDO
      ENDDO
   ENDDO
ELSEIF(TENSOR%primCStensor)THEN
 DO I = 1,TENSOR%nLSAO
  WRITE(LUDISK)TENSOR%LSAO(I)%ALLOC
  IF(TENSOR%LSAO(I)%ALLOC)THEN
   WRITE(LUDISK)TENSOR%LSAO(I)%ATOM1
   WRITE(LUDISK)TENSOR%LSAO(I)%ATOM2
   WRITE(LUDISK)TENSOR%LSAO(I)%nAOBATCH1
   WRITE(LUDISK)TENSOR%LSAO(I)%nAOBATCH2
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     nP1 = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompA(1)
     nP2 = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompB(1)
     WRITE(LUDISK)nP1
     WRITE(LUDISK)nP2
     DO IP = 1,nP1*nP2
        WRITE(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(IP)
     ENDDO
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ELSE
 DO I = 1,TENSOR%nLSAO
  WRITE(LUDISK)TENSOR%LSAO(I)%ALLOC
  IF(TENSOR%LSAO(I)%ALLOC)THEN
   WRITE(LUDISK)TENSOR%LSAO(I)%ATOM1
   WRITE(LUDISK)TENSOR%LSAO(I)%ATOM2
   WRITE(LUDISK)TENSOR%LSAO(I)%ATOM3
   WRITE(LUDISK)TENSOR%LSAO(I)%ATOM4
   WRITE(LUDISK)TENSOR%LSAO(I)%nAOBATCH1
   WRITE(LUDISK)TENSOR%LSAO(I)%nAOBATCH2
   WRITE(LUDISK)TENSOR%LSAO(I)%nAOBATCH3
   WRITE(LUDISK)TENSOR%LSAO(I)%nAOBATCH4
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
      DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
       WRITE(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
       WRITE(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
       WRITE(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
       WRITE(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
       DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
       ENDDO
       DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
       ENDDO
       DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)
       ENDDO
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        WRITE(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
       ENDDO
       WRITE(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE
       WRITE(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms      
       do IELM = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms*nmat
          WRITE(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ENDIF
call Determine_lstensor_memory(tensor,nmemsize)
END SUBROUTINE WRITE_LSTENSOR_TO_DISK

!> \brief read the lstensor structure from disk
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor to be written to disk
!> \param LUDISK the logical unit number of file to write to
SUBROUTINE read_lstensor_from_disk(TENSOR,LUDISK,LUPRI)
implicit none
type(lstensor) :: TENSOR 
integer :: LUDISK,lupri
!
integer :: iatom,I,J,IMAT2,ibat,jbat,DIM,nP1,nP2,IP,nbat3,nbat4
integer :: kbat,lbat,iangmom,jangmom,kangmom,langmom,nbat1,nbat2
integer :: Ielm,IORB,JORB,KORB,LORB,nContA,nOrbCompA,nmat
integer :: atom1,atom2,atom3,atom4
integer :: nContB,nOrbCompB,nContC,nOrbCompC,nContD,nOrbCompD
integer(kind=long) :: nmemsize

READ(LUDISK)TENSOR%primCStensor 
READ(LUDISK)TENSOR%gradienttensor 
READ(LUDISK)TENSOR%natom1
READ(LUDISK)TENSOR%natom2
READ(LUDISK)TENSOR%natom3
READ(LUDISK)TENSOR%natom4
READ(LUDISK)TENSOR%nbast1
READ(LUDISK)TENSOR%nbast2
READ(LUDISK)TENSOR%nbast3
READ(LUDISK)TENSOR%nbast4
READ(LUDISK)TENSOR%nmat
READ(LUDISK)TENSOR%nLSAO
nmat = TENSOR%nmat
IF(TENSOR%gradienttensor)THEN
  CALL read_lstensor_from_disk_gradient(TENSOR,nmat,LUDISK)
ELSEIF(TENSOR%primCStensor)THEN
  call read_lstensor_from_disk_primcs(TENSOR,nmat,LUDISK)
ELSE
 NULLIFY(TENSOR%LSAO)
 NULLIFY(TENSOR%INDEX)
 ALLOCATE(TENSOR%LSAO(TENSOR%natom1*TENSOR%natom2*TENSOR%natom3*TENSOR%natom4))
 ALLOCATE(TENSOR%INDEX(TENSOR%natom1,TENSOR%natom2,TENSOR%natom3,TENSOR%natom4))
 TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
 DO I = 1,TENSOR%nLSAO
  READ(LUDISK)TENSOR%LSAO(I)%ALLOC
  IF(TENSOR%LSAO(I)%ALLOC)THEN
   READ(LUDISK)TENSOR%LSAO(I)%ATOM1
   atom1=TENSOR%LSAO(I)%ATOM1
   READ(LUDISK)TENSOR%LSAO(I)%ATOM2
   atom2=TENSOR%LSAO(I)%ATOM2
   READ(LUDISK)TENSOR%LSAO(I)%ATOM3
   atom3=TENSOR%LSAO(I)%ATOM3
   READ(LUDISK)TENSOR%LSAO(I)%ATOM4
   atom4=TENSOR%LSAO(I)%ATOM4
   TENSOR%INDEX(atom1,atom2,atom3,atom4)=I
   READ(LUDISK)TENSOR%LSAO(I)%nAOBATCH1
   nbat1=TENSOR%LSAO(I)%nAOBATCH1
   READ(LUDISK)TENSOR%LSAO(I)%nAOBATCH2
   nbat2=TENSOR%LSAO(I)%nAOBATCH2
   READ(LUDISK)TENSOR%LSAO(I)%nAOBATCH3
   nbat3=TENSOR%LSAO(I)%nAOBATCH3
   READ(LUDISK)TENSOR%LSAO(I)%nAOBATCH4
   nbat4=TENSOR%LSAO(I)%nAOBATCH4
   NULLIFY(TENSOR%LSAO(I)%BATCH)
   ALLOCATE(TENSOR%LSAO(I)%BATCH(nbat1,nbat2,nbat3,nbat4))
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
      DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
       READ(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
       READ(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
       READ(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
       READ(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
       DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
       ENDDO
       DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
       ENDDO
       DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)
       ENDDO
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        READ(LUDISK)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
       ENDDO
       READ(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE
       READ(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms
       DIM = TENSOR%nmat*TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms
       NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms)
       ALLOCATE(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(DIM))
       DO IELM = 1,DIM*nmat
          READ(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ELSE
   TENSOR%LSAO(I)%ATOM1 = 0
   TENSOR%LSAO(I)%ATOM2 = 0
   TENSOR%LSAO(I)%ATOM3 = 0
   TENSOR%LSAO(I)%ATOM4 = 0
   TENSOR%LSAO(I)%nAOBATCH1 = 0
   TENSOR%LSAO(I)%nAOBATCH2 = 0
   TENSOR%LSAO(I)%nAOBATCH3 = 0
   TENSOR%LSAO(I)%nAOBATCH4 = 0
  ENDIF
 ENDDO
ENDIF

Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE read_lstensor_from_disk

!> \brief read the lstensor structure from disk, in case of primitvescreening tensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor to be written to disk
!> \param nmat the number of matrices
!> \param LUDISK the logical unit number of file to write to
SUBROUTINE read_lstensor_from_disk_primcs(TENSOR,nmat,LUDISK)
implicit none
type(lstensor) :: TENSOR 
integer :: LUDISK
!
integer :: iatom,I,J,IMAT2,ibat,jbat,DIM,nP1,nP2,IP,nbat3,nbat4
integer :: kbat,lbat,iangmom,jangmom,kangmom,langmom,nbat1,nbat2
integer :: Ielm,IORB,JORB,KORB,LORB,nContA,nOrbCompA,nmat
integer :: atom1,atom2,atom3,atom4
integer :: nContB,nOrbCompB,nContC,nOrbCompC,nContD,nOrbCompD
integer(kind=long) :: nmemsize
!
 NULLIFY(TENSOR%LSAO)
 NULLIFY(TENSOR%INDEX)
 ALLOCATE(TENSOR%LSAO(TENSOR%natom1*TENSOR%natom2))
 ALLOCATE(TENSOR%INDEX(TENSOR%natom1,TENSOR%natom2,1,1))
 TENSOR%INDEX(TENSOR%natom1,TENSOR%natom2,1,1) = 0
 DO I = 1,TENSOR%nLSAO
  READ(LUDISK)TENSOR%LSAO(I)%ALLOC
  IF(TENSOR%LSAO(I)%ALLOC)THEN
   READ(LUDISK)TENSOR%LSAO(I)%ATOM1
   READ(LUDISK)TENSOR%LSAO(I)%ATOM2
   TENSOR%INDEX(TENSOR%LSAO(I)%ATOM1,TENSOR%LSAO(I)%ATOM2,1,1) = I
   TENSOR%LSAO(I)%ATOM3 = 1
   TENSOR%LSAO(I)%ATOM4 = 1
   READ(LUDISK)TENSOR%LSAO(I)%nAOBATCH1
   READ(LUDISK)TENSOR%LSAO(I)%nAOBATCH2
   TENSOR%LSAO(I)%nAOBATCH3 = 1
   TENSOR%LSAO(I)%nAOBATCH4 = 1
   nbat1=TENSOR%LSAO(I)%nAOBATCH1
   nbat2=TENSOR%LSAO(I)%nAOBATCH2
   NULLIFY(TENSOR%LSAO(I)%BATCH)
   ALLOCATE(TENSOR%LSAO(I)%BATCH(nbat1,nbat2,1,1))
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     READ(LUDISK)nP1
     READ(LUDISK)nP2
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomA = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomB = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomC = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomD = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContA(1) = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompA(1) = nP1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbA(1) = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContB(1) = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompB(1) = nP2
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbB(1) = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContC(1) = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompC(1) = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbC(1) = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContD(1) = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompD(1) = 1 
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbD(1) = 1
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nelmE = 1
     DIM = nP1*nP2
     NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms)
     ALLOCATE(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(DIM))
     TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nelms = DIM
     DO IP = 1,DIM
        READ(LUDISK) TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(IP)
     ENDDO
    ENDDO
   ENDDO
  ELSE
   TENSOR%LSAO(I)%ATOM1 = 0
   TENSOR%LSAO(I)%ATOM2 = 0
   TENSOR%LSAO(I)%ATOM3 = 0
   TENSOR%LSAO(I)%ATOM4 = 0
   TENSOR%LSAO(I)%nAOBATCH1 = 0
   TENSOR%LSAO(I)%nAOBATCH2 = 0
   TENSOR%LSAO(I)%nAOBATCH3 = 0
   TENSOR%LSAO(I)%nAOBATCH4 = 0
  ENDIF
 ENDDO

END SUBROUTINE read_lstensor_from_disk_primcs

!> \brief read the lstensor structure from disk, in case of gradient tensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor to be written to disk
!> \param nmat the number of matrices
!> \param LUDISK the logical unit number of file to write to
SUBROUTINE read_lstensor_from_disk_gradient(TENSOR,nmat,LUDISK)
implicit none
type(lstensor) :: TENSOR 
integer :: LUDISK
!
integer :: iatom,I,J,IMAT2,ibat,jbat,DIM,nP1,nP2,IP,nbat3,nbat4
integer :: kbat,lbat,iangmom,jangmom,kangmom,langmom,nbat1,nbat2
integer :: Ielm,IORB,JORB,KORB,LORB,nContA,nOrbCompA,nmat
integer :: atom1,atom2,atom3,atom4
integer :: nContB,nOrbCompB,nContC,nOrbCompC,nContD,nOrbCompD
integer(kind=long) :: nmemsize
!
NULLIFY(TENSOR%LSAO)
NULLIFY(TENSOR%INDEX)
ALLOCATE(TENSOR%LSAO(TENSOR%natom1))
ALLOCATE(TENSOR%INDEX(TENSOR%natom1,1,1,1))
TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
DO I= 1,TENSOR%natom1
   TENSOR%INDEX(I,1,1,1) = I
   NULLIFY(TENSOR%LSAO(I)%BATCH)
   ALLOCATE(TENSOR%LSAO(I)%BATCH(1,1,1,1))
   TENSOR%LSAO(I)%ALLOC = .TRUE.
   TENSOR%LSAO(I)%ATOM1 = I
   TENSOR%LSAO(I)%ATOM2 = 0
   TENSOR%LSAO(I)%ATOM3 = 0
   TENSOR%LSAO(I)%ATOM4 = 0
   TENSOR%LSAO(I)%nAOBATCH1 = 0
   TENSOR%LSAO(I)%nAOBATCH2 = 0
   TENSOR%LSAO(I)%nAOBATCH3 = 0
   TENSOR%LSAO(I)%nAOBATCH4 = 0
   TENSOR%LSAO(I)%BATCH(1,1,1,1)%nelmE = 3*nmat
   NULLIFY(TENSOR%LSAO(I)%BATCH(1,1,1,1)%elms)
   ALLOCATE(TENSOR%LSAO(I)%BATCH(1,1,1,1)%elms(3*nmat))
   DO J=1,3
      DO IMAT2 = 1,nmat
         READ(LUDISK) TENSOR%LSAO(I)%BATCH(1,1,1,1)%elms(J+(IMAT2-1)*3)
      ENDDO
   ENDDO
ENDDO

Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE read_lstensor_from_disk_gradient

!> \brief make an empty AOITEM
!> \author T. Kjaergaard
!> \date 2010
!> \param AO the AOITEM to be built
SUBROUTINE SET_EMPTY_AO(AO)
IMPLICIT NONE
TYPE(AOITEM)  :: AO

AO%natoms = 1                       
ALLOCATE(AO%ATOMICnORB(AO%natoms))  
ALLOCATE(AO%ATOMICnBATCH(AO%natoms))
AO%ATOMICnORB(1)=1                  
AO%ATOMICnBATCH(1)=1                
AO%nbast = 1                        
AO%empty=.TRUE.
AO%nbatches=1
AO%nCC=1
AO%nExp=1
NULLIFY(AO%BATCH)
ALLOCATE(AO%BATCH(1))   
AO%BATCH(1)%batch=1
AO%BATCH(1)%nAngmom=1
AO%BATCH(1)%nContracted(1)=1
AO%BATCH(1)%startOrbital(1)=1
AO%BATCH(1)%nOrbComp(1)=1
AO%BATCH(1)%nPrimitives=1

END SUBROUTINE SET_EMPTY_AO

!> \brief free an empty AOITEM
!> \author T. Kjaergaard
!> \date 2010
!> \param AO free the AOITEM
SUBROUTINE FREE_EMPTY_AO(AO)
IMPLICIT NONE
TYPE(AOITEM)  :: AO

DEALLOCATE(AO%ATOMICnORB)  
NULLIFY(AO%ATOMICnORB)
DEALLOCATE(AO%ATOMICnBATCH)
NULLIFY(AO%ATOMICnBATCH)
DEALLOCATE(AO%BATCH)   
NULLIFY(AO%BATCH)

END SUBROUTINE FREE_EMPTY_AO

!> \brief print the lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor to be written to disk
!> \param doall print all - all elements, compared to just overall info
!> \param LUPRI logical unit number for the default output file
SUBROUTINE Print_lstensor(TENSOR,LUPRI,doall)
implicit none
TYPE(LSTENSOR)     :: TENSOR
INTEGER            :: LUPRI
logical            :: doall

INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,iatom,nmat,jatom,katom,latom
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT,i1,np1,np2,istart,iend

WRITE(LUPRI,'(A)') 'Printing lstensor'
WRITE(LUPRI,'(A,L)') 'TENSOR%primCStensor',TENSOR%primCStensor 
WRITE(LUPRI,'(A,L)') 'TENSOR%gradienttensor',TENSOR%gradienttensor 
WRITE(LUPRI,'(A,I7)') 'TENSOR%natom1',TENSOR%natom1
WRITE(LUPRI,'(A,I7)') 'TENSOR%natom2',TENSOR%natom2
WRITE(LUPRI,'(A,I7)') 'TENSOR%natom3',TENSOR%natom3
WRITE(LUPRI,'(A,I7)') 'TENSOR%natom4',TENSOR%natom4
WRITE(LUPRI,'(A,I7)') 'TENSOR%nbast1',TENSOR%nbast1
WRITE(LUPRI,'(A,I7)') 'TENSOR%nbast2',TENSOR%nbast2
WRITE(LUPRI,'(A,I7)') 'TENSOR%nbast3',TENSOR%nbast3
WRITE(LUPRI,'(A,I7)') 'TENSOR%nbast4',TENSOR%nbast4
WRITE(LUPRI,'(A,I7)') 'TENSOR%nmat  ',TENSOR%nmat
WRITE(LUPRI,'(A,I7)') 'TENSOR%nLSAO ',TENSOR%nLSAO
IF(TENSOR%gradienttensor)THEN
 DO Iatom=1,TENSOR%natom1
  DO Jatom=1,TENSOR%natom2
   DO Katom=1,TENSOR%natom3
    DO Latom=1,TENSOR%natom4
     WRITE(LUPRI,'(A,I3,A,I3,A,I3,A,I3,A,I4)')'INDEX(',Iatom,',',Jatom,',',Katom,',',Latom,')=',&
          &TENSOR%INDEX(Iatom,Jatom,Katom,Latom)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
 WRITE(LUPRI,'(A4,1X,A1,20X,A1,21X,A1,21X,A5)')'Atom ','X','Y','Z','ndmat'
 nmat = TENSOR%nmat
 DO Iatom = 1,TENSOR%natom1
    DO IMAT = 1,nmat
       WRITE(lupri,'(I4,3F22.10,I5)')Iatom,(TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms(J+(IMAT-1)*3),J=1,3),IMAT
    ENDDO
 ENDDO
ELSEIF(TENSOR%primCStensor)THEN
 DO Iatom=1,TENSOR%natom1
  DO Jatom=1,TENSOR%natom2
   DO Katom=1,TENSOR%natom3
    DO Latom=1,TENSOR%natom4
     WRITE(LUPRI,'(A,I3,A,I3,A,I3,A,I3,A,I4)')'INDEX(',Iatom,',',Jatom,',',Katom,',',Latom,')=',&
          &TENSOR%INDEX(Iatom,Jatom,Katom,Latom)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
 DO I = 1,TENSOR%nLSAO
  IF(TENSOR%LSAO(I)%ALLOC)THEN
   WRITE(LUPRI,'(A,5I3)') 'ATOMS:',TENSOR%LSAO(I)%ATOM1,TENSOR%LSAO(I)%ATOM2
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     WRITE(LUPRI,'(A,5I3)') 'I,Ibat,Jbat,Kbat,Lbat',I,Ibat,Jbat
     nP1 = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompA(1)
     nP2 = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompB(1)
     DO IMAT = 1,TENSOR%nmat
      ISTART = 1+(IMAT-1)*nP1*nP2
      IEND = IMAT*nP1*nP2
      WRITE(lupri,*)'The primGAB matrix for IMAT=',IMAT,'nP1,nP2',nP1,nP2
!      print*,'TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms',TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms
      call output(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(ISTART:IEND),1,nP1,1,nP2,nP1,nP2,1,lupri)
     ENDDO
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ELSE
 DO Iatom=1,TENSOR%natom1
  DO Jatom=1,TENSOR%natom2
   DO Katom=1,TENSOR%natom3
    DO Latom=1,TENSOR%natom4
     WRITE(LUPRI,'(A,I3,A,I3,A,I3,A,I3,A,I4)')'INDEX(',Iatom,',',Jatom,',',Katom,',',Latom,')=',&
          &TENSOR%INDEX(Iatom,Jatom,Katom,Latom)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  WRITE(LUPRI,'(A,5I3)') 'ATOM  :',TENSOR%LSAO(I)%ATOM1,TENSOR%LSAO(I)%ATOM2,&
       &TENSOR%LSAO(I)%ATOM3,TENSOR%LSAO(I)%ATOM4
  WRITE(LUPRI,'(A,5I3)') 'nBATCH:',TENSOR%LSAO(I)%nAOBATCH1,TENSOR%LSAO(I)%nAOBATCH2,&
       &TENSOR%LSAO(I)%nAOBATCH3,TENSOR%LSAO(I)%nAOBATCH4
  if(doall)then
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
      DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
       WRITE(LUPRI,'(A,5I3)') 'I,Ibat,Jbat,Kbat,Lbat',I,Ibat,Jbat,Kbat,Lbat
       IELM = 0
       DO IMAT = 1,TENSOR%nmat                    
        DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
         nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
         nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
         sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
         DO LORB = 0,nContD*nOrbCompD-1    
          DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
           nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
           nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
           sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
           DO KORB = 0,nContC*nOrbCompC-1
                
            DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
             nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
             nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
             sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
             DO JORB = 0,nContB*nOrbCompB-1
                  
              DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
               nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
               nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
               sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
               WRITE(LUPRI,'(3X,A,5I3)') 'orbitalindex:',sA,sB+JORB,sC+KORB,sD+LORB,IMAT
               do i1=iELM+1,iElM+nContA*nOrbCompA
                WRITE(LUPRI,*) 'Ielms=',i1
                if(i1 .GT. TENSOR%nmat*TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms)call lsquit('error in print lstensor',-1)
                  !                         WRITE(LUPRI,'(5F16.8)')TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(i1)
                WRITE(LUPRI,*)TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(i1)
               enddo
                 !                      WRITE(LUPRI,'(5X,5F12.8,/(5X,5F12.8))')(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%&
                 !                           &elms(i1),i1=iELM+1,iElM+nContA*nOrbCompA)
               IELM = IELM+nContA*nOrbCompA
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
    ENDDO
   ENDDO
  endif
 ENDIF
ENDDO
ENDIF

END SUBROUTINE PRINT_LSTENSOR

!> \brief copy an lstensor to a new lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR1 the original lstensor
!> \param TENSOR2 the new lstensor
SUBROUTINE copy_lstensor_to_lstensor(TENSOR1,TENSOR2)
  implicit none
  TYPE(LSTENSOR)     :: TENSOR1,TENSOR2

  INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,nmat,dim,ielms
  INTEGER    :: nAngmomA,nAngmomB,nAngmomC,nAngmomD
  integer(kind=long) :: nmemsize

  TENSOR2%primCStensor = TENSOR1%primCStensor
  TENSOR2%gradienttensor = TENSOR1%gradienttensor
  TENSOR2%natom1 = TENSOR1%natom1 
  TENSOR2%natom2 = TENSOR1%natom2 
  TENSOR2%natom3 = TENSOR1%natom3 
  TENSOR2%natom4 = TENSOR1%natom4 
  TENSOR2%nbast1 = TENSOR1%nbast1 
  TENSOR2%nbast2 = TENSOR1%nbast2 
  TENSOR2%nbast3 = TENSOR1%nbast3 
  TENSOR2%nbast4 = TENSOR1%nbast4 
  TENSOR2%nmat   = TENSOR1%nmat 
  TENSOR2%nLSAO  = TENSOR1%nLSAO

  NULLIFY(TENSOR2%LSAO)
  NULLIFY(TENSOR2%INDEX)
  ALLOCATE(TENSOR2%LSAO(TENSOR2%nLSAO))
  ALLOCATE(TENSOR2%INDEX(TENSOR2%natom1,TENSOR2%natom2,TENSOR2%natom3,TENSOR2%natom4))

  DO Ibat = 1,TENSOR1%natom1 
     DO Jbat = 1,TENSOR1%natom2
        DO Kbat = 1,TENSOR1%natom3
           DO Lbat = 1,TENSOR1%natom4 
              TENSOR2%INDEX(Ibat,Jbat,Kbat,Lbat) = TENSOR1%INDEX(Ibat,Jbat,Kbat,Lbat)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  nmat = TENSOR2%nmat
  DO I = 1,TENSOR1%nLSAO
   TENSOR2%LSAO(I)%ALLOC = TENSOR1%LSAO(I)%ALLOC
   IF(TENSOR1%LSAO(I)%ALLOC)THEN    
    TENSOR2%LSAO(I)%ATOM1 = TENSOR1%LSAO(I)%ATOM1 
    TENSOR2%LSAO(I)%ATOM2 = TENSOR1%LSAO(I)%ATOM2 
    TENSOR2%LSAO(I)%ATOM3 = TENSOR1%LSAO(I)%ATOM3 
    TENSOR2%LSAO(I)%ATOM4 = TENSOR1%LSAO(I)%ATOM4 
    TENSOR2%LSAO(I)%nAOBATCH1 = TENSOR1%LSAO(I)%nAOBATCH1 
    TENSOR2%LSAO(I)%nAOBATCH2 = TENSOR1%LSAO(I)%nAOBATCH2 
    TENSOR2%LSAO(I)%nAOBATCH3 = TENSOR1%LSAO(I)%nAOBATCH3
    TENSOR2%LSAO(I)%nAOBATCH4 = TENSOR1%LSAO(I)%nAOBATCH4 
    NULLIFY(TENSOR2%LSAO(I)%BATCH)
    ALLOCATE(TENSOR2%LSAO(I)%BATCH(TENSOR2%LSAO(I)%nAOBATCH1 ,TENSOR2%LSAO(I)%nAOBATCH2 ,&
         &TENSOR2%LSAO(I)%nAOBATCH3 ,TENSOR2%LSAO(I)%nAOBATCH4))
    DO Ibat = 1,TENSOR1%LSAO(I)%nAOBATCH1
     DO Jbat = 1,TENSOR1%LSAO(I)%nAOBATCH2
      DO Kbat = 1,TENSOR1%LSAO(I)%nAOBATCH3
       DO Lbat = 1,TENSOR1%LSAO(I)%nAOBATCH4
        dim = nmat*TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms
        NULLIFY(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms)
        ALLOCATE(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(dim))
        
        do ielms = 1,dim
           TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(ielms) = &
                &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(ielms)         
        enddo
        
        nAngmomA = TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
        nAngmomB = TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
        nAngmomC = TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
        nAngmomD = TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA = nAngmomA 
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB = nAngmomB
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC = nAngmomC
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD = nAngmomD
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(1:nAngmomA) = &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(1:nAngmomA)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(1:nAngmomA)= &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(1:nAngmomA)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(1:nAngmomA)= &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(1:nAngmomA)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(1:nAngmomB)   = &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(1:nAngmomB)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(1:nAngmomB)= &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(1:nAngmomB)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(1:nAngmomB)= &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(1:nAngmomB)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(1:nAngmomC)   = &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(1:nAngmomC)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(1:nAngmomC)= &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(1:nAngmomC)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(1:nAngmomC)= &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(1:nAngmomC)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(1:nAngmomD)   = &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(1:nAngmomD)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(1:nAngmomD)= &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(1:nAngmomD)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(1:nAngmomD)= &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(1:nAngmomD)
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE = &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE
        TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms = &
             &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms
       enddo
      enddo
     enddo
    enddo
   else
      TENSOR2%LSAO(I)%ATOM1 = 0
      TENSOR2%LSAO(I)%ATOM2 = 0
      TENSOR2%LSAO(I)%ATOM3 = 0
      TENSOR2%LSAO(I)%ATOM4 = 0
      TENSOR2%LSAO(I)%nAOBATCH1 = 0
      TENSOR2%LSAO(I)%nAOBATCH2 = 0
      TENSOR2%LSAO(I)%nAOBATCH3 = 0
      TENSOR2%LSAO(I)%nAOBATCH4 = 0
   endif
  enddo
  Call Determine_lstensor_memory(tensor2,nmemsize)
  call add_mem_to_global(nmemsize)

end SUBROUTINE copy_lstensor_to_lstensor

!> \brief copy an lstensor to a new lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR1 the original lstensor
!> \param TENSOR2 the new lstensor
SUBROUTINE lstensor_compare(TENSOR1,TENSOR2)
implicit none
TYPE(LSTENSOR)     :: TENSOR1,TENSOR2

INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,nmat,dim,ielms
INTEGER    :: nAngmomA,nAngmomB,nAngmomC,nAngmomD
INTEGER    :: iAngmomA,iAngmomB,iAngmomC,iAngmomD
integer(kind=long) :: nmemsize
print*,'lstensor_compare'
IF(TENSOR2%primCStensor.NEQV.TENSOR1%primCStensor)call lsquit('EW2',6)
IF(TENSOR2%gradienttensor.NEQV.TENSOR1%gradienttensor)call lsquit('EW3',6)
IF(TENSOR2%natom1.NE.TENSOR1%natom1)call lsquit('EW4',6) 
IF(TENSOR2%natom2.NE.TENSOR1%natom2)call lsquit('EW5',6)
IF(TENSOR2%natom3.NE.TENSOR1%natom3)call lsquit('EW6',6) 
IF(TENSOR2%natom4.NE.TENSOR1%natom4)call lsquit('EW7',6) 
IF(TENSOR2%nbast1.NE.TENSOR1%nbast1)call lsquit('EW8',6) 
IF(TENSOR2%nbast2.NE.TENSOR1%nbast2)call lsquit('EW9',6) 
IF(TENSOR2%nbast3.NE.TENSOR1%nbast3)call lsquit('EW10',6) 
IF(TENSOR2%nbast4.NE.TENSOR1%nbast4)call lsquit('EW11',6) 
IF(TENSOR2%nmat.NE.TENSOR1%nmat)call lsquit('EW12',6) 
IF(TENSOR2%nLSAO .NE.TENSOR1%nLSAO)call lsquit('EW13',6)

DO Ibat = 1,TENSOR1%natom1 
   DO Jbat = 1,TENSOR1%natom2
      DO Kbat = 1,TENSOR1%natom3
         DO Lbat = 1,TENSOR1%natom4 
            IF(TENSOR2%INDEX(Ibat,Jbat,Kbat,Lbat).NE.TENSOR1%INDEX(Ibat,Jbat,Kbat,Lbat))call lsquit('EW14',6)
         ENDDO
      ENDDO
   ENDDO
ENDDO

DO I = 1,TENSOR1%nLSAO
 IF(TENSOR2%LSAO(I)%ALLOC.NEQV.TENSOR1%LSAO(I)%ALLOC)call lsquit('EW15',6)
 IF(TENSOR1%LSAO(I)%ALLOC)THEN    
  IF(TENSOR2%LSAO(I)%ATOM1.NE.TENSOR1%LSAO(I)%ATOM1)call lsquit('EW16',6) 
  IF(TENSOR2%LSAO(I)%ATOM2.NE.TENSOR1%LSAO(I)%ATOM2)call lsquit('EW17',6) 
  IF(TENSOR2%LSAO(I)%ATOM3.NE.TENSOR1%LSAO(I)%ATOM3)call lsquit('EW18',6) 
  IF(TENSOR2%LSAO(I)%ATOM4.NE.TENSOR1%LSAO(I)%ATOM4)call lsquit('EW19',6) 
  IF(TENSOR2%LSAO(I)%nAOBATCH1.NE.TENSOR1%LSAO(I)%nAOBATCH1)call lsquit('EW20',6) 
  IF(TENSOR2%LSAO(I)%nAOBATCH2.NE.TENSOR1%LSAO(I)%nAOBATCH2)call lsquit('EW21',6) 
  IF(TENSOR2%LSAO(I)%nAOBATCH3.NE.TENSOR1%LSAO(I)%nAOBATCH3)call lsquit('EW22',6)
  IF(TENSOR2%LSAO(I)%nAOBATCH4.NE.TENSOR1%LSAO(I)%nAOBATCH4)call lsquit('EW23',6) 
  DO Ibat = 1,TENSOR1%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR1%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR1%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR1%LSAO(I)%nAOBATCH4
      IF(TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms.NE.&
           &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms)call lsquit('EW24',6)
      do ielms = 1,TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms
         IF(ABS(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(ielms)-&
              &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(ielms)).GT.1E-10)call lsquit('EW25',6)         
      enddo

nAngmomA = TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
nAngmomB = TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
nAngmomC = TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
nAngmomD = TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD

IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA.NE.TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA)&
     &call lsquit('EW26',6)
IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB.NE.TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB)&
     &call lsquit('EW27',6)  
IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC.NE.TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC)&
     &call lsquit('EW28',6)  
IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD.NE.TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD)&
     &call lsquit('EW29',6)   
do IangmomA=1,nangmomA
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(iAngmomA) .NE. TENSOR1%&
        &LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(iAngmomA))call lsquit('EW30',6)
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(iAngmomA) .NE. TENSOR1%&
        &LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(iAngmomA))call lsquit('EW31',6)
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(iAngmomA) .NE. TENSOR1%&
        &LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(iAngmomA))call lsquit('EW32',6)
enddo
do IangmomB=1,nangmomB
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(iAngmomB) .NE. TENSOR1%LSAO(I)%&
        &BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(iAngmomB))call lsquit('EW33',6)
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(iAngmomB) .NE. TENSOR1%LSAO(I)%&
        &BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(iAngmomB))call lsquit('EW34',6)
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(iAngmomB) .NE. TENSOR1%LSAO(I)%&
        &BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(iAngmomB))call lsquit('EW35',6)
enddo
do IangmomC=1,nangmomC
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(iAngmomC) .NE. TENSOR1%LSAO(I)%&
        &BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(iAngmomC))call lsquit('EW36',6)
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(iAngmomC) .NE. TENSOR1%LSAO(I)%&
        &BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(iAngmomC))call lsquit('EW37',6)
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(iAngmomC) .NE. TENSOR1%LSAO(I)%&
        &BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(iAngmomC))call lsquit('EW38',6)
enddo
do IangmomD=1,nangmomD
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(iAngmomD) .NE. TENSOR1%LSAO(I)%&
        &BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(iAngmomD))call lsquit('EW39',6)
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(iAngmomD) .NE. TENSOR1%LSAO(I)%&
        &BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(iAngmomD))call lsquit('EW40',6)
   IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(iAngmomD) .NE. TENSOR1%LSAO(I)%&
        &BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(iAngmomD))call lsquit('EW41',6)
enddo
IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE .NE. &
     &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE)call lsquit('EW42',6)
IF(TENSOR2%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms .NE. &
     &TENSOR1%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms)call lsquit('EW43',6)
     enddo
    enddo
   enddo
  enddo
 endif
enddo

end SUBROUTINE lstensor_compare

!> \brief free lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
SUBROUTINE lstensor_free(TENSOR)
implicit none
TYPE(LSTENSOR)     :: TENSOR
INTEGER :: I,Ibat,Jbat,Kbat,Lbat
integer(kind=long) :: nmemsize

Call Determine_lstensor_memory(tensor,nmemsize)
call remove_mem_from_global(nmemsize)

IF(TENSOR%gradienttensor)THEN
 DO I = 1,TENSOR%nLSAO
  IF(TENSOR%LSAO(I)%ALLOC)THEN
    IF(ASSOCIATED(TENSOR%LSAO(I)%BATCH(1,1,1,1)%elms))THEN
       DEALLOCATE(TENSOR%LSAO(I)%BATCH(1,1,1,1)%elms)
       NULLIFY(TENSOR%LSAO(I)%BATCH(1,1,1,1)%elms)
    ENDIF
    DEALLOCATE(TENSOR%LSAO(I)%BATCH)
    NULLIFY(TENSOR%LSAO(I)%BATCH)
  ENDIF
 ENDDO
 DEALLOCATE(TENSOR%LSAO)
 DEALLOCATE(TENSOR%INDEX)
 NULLIFY(TENSOR%LSAO)
 NULLIFY(TENSOR%INDEX)
ELSE
 DO I = 1,TENSOR%nLSAO
  IF(TENSOR%LSAO(I)%ALLOC)THEN
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
      DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
       IF(ASSOCIATED(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms))THEN
        DEALLOCATE(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms)
        NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms)
       ENDIF
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   DEALLOCATE(TENSOR%LSAO(I)%BATCH)
   NULLIFY(TENSOR%LSAO(I)%BATCH)
  ENDIF
 ENDDO
 DEALLOCATE(TENSOR%LSAO)
 DEALLOCATE(TENSOR%INDEX)
 NULLIFY(TENSOR%LSAO)
 NULLIFY(TENSOR%INDEX)
ENDIF

END SUBROUTINE lstensor_free

!> \brief set the lstensor to zero
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
SUBROUTINE lstensor_zero(TENSOR)
implicit none
TYPE(LSTENSOR)     :: TENSOR
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT,i1

IF(TENSOR%primCStensor)CALL lsQUIT('this type of lstensor should not be used in this fashion2',-1)
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      IELM = 0
      DO IMAT = 1,TENSOR%nmat                     
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
          DO KORB = 0,nContC*nOrbCompC-1
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1                           
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
              do i1=iELM+1,iElM+nContA*nOrbCompA
               TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(i1) = 0.d0
              enddo
              IELM = IELM+nContA*nOrbCompA
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
   ENDDO
  ENDDO
 ENDIF
ENDDO

END SUBROUTINE LSTENSOR_ZERO

!> \brief build full 5 dimensional array from an lstensor 
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT134 the 5 dimensional array
!> \param dim1 the size of the 1. dimension
!> \param dim2 the size of the 2. dimension
!> \param dim3 the size of the 3. dimension
!> \param dim4 the size of the 4. dimension
!> \param dim5 the size of the 5. dimension
SUBROUTINE Build_full_5dim_from_lstensor(TENSOR,MAT134,dim1,dim2,dim3,dim4,dim5)
implicit none
integer          :: dim1,dim2,dim3,dim4,dim5
TYPE(LSTENSOR)   :: TENSOR
REAL(REALK)      :: MAT134(dim1,dim2,dim3,dim4,dim5)
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT
LOGICAL    :: ERROR

IF(TENSOR%primCStensor)CALL lsQUIT('this type of lstensor should not be used in this fashion3',-1)
ERROR=.FALSE.
IF(TENSOR%nbast1 .NE. dim1) ERROR = .TRUE.
IF(TENSOR%nbast2 .NE. dim2) ERROR = .TRUE.
IF(TENSOR%nbast3 .NE. dim3) ERROR = .TRUE.
IF(TENSOR%nbast4 .NE. dim4) ERROR = .TRUE.
IF(ERROR)THEN
   print*,'LSTENSOR%nbast1',TENSOR%nbast1,'dim1',dim1
   print*,'LSTENSOR%nbast2',TENSOR%nbast2,'dim2',dim2
   print*,'LSTENSOR%nbast3',TENSOR%nbast3,'dim3',dim3
   print*,'LSTENSOR%nbast4',TENSOR%nbast4,'dim4',dim4
!   print*,'LSTENSOR%nbast5',TENSOR%nbast5,'dim5',dim5
   CALL lsQUIT('ERROR: Build_full_5dim_from_lstensor dim do not match',-1)
ENDIF
MAT134 = 0.d0
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      IELM = 0
      DO IMAT = 1,TENSOR%nmat
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
          DO KORB = 0,nContC*nOrbCompC-1
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1                           
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               IELM = IELM+1
               MAT134(sA+IORB,sB+JORB,sC+KORB,sD+LORB,IMAT) &
                    & = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
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
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDDO

END SUBROUTINE Build_full_5dim_from_lstensor

!> \brief build full 2 dimensional array from an lstensor 
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the 2 dimensional array
!> \param dim1 the size of the 1. dimension
!> \param dim2 the size of the 2. dimension
SUBROUTINE Build_full_2dim_from_lstensor(TENSOR,MAT2,dim1,dim2)
implicit none
INTEGER            :: dim1,dim2
TYPE(LSTENSOR)     :: TENSOR
real(realk)        :: MAT2(dim1,dim2)
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT
INTEGER    :: nbast(4),MATINDEX(4)

IF(TENSOR%primCStensor)CALL lsQUIT('this type of lstensor should not be used in this fashion4',-1)
nbast(1)=TENSOR%nbast1
nbast(2)=TENSOR%nbast2
nbast(3)=TENSOR%nbast3
nbast(4)=TENSOR%nbast4
IF(TENSOR%nmat.NE.1)CALL lsQUIT('called Build_full_2dim_from_lstensor with TENSOR%nmat.NE.1',-1)
DO I=1,4
   IF(dim1 .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(dim2 .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=dim1+dim2
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_matrix_from_lstensor dim do not match',-1)
do I=1,dim1
   do J=1,dim2
      MAT2(I,J) = 0.d0
   enddo
enddo
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      IELM = 0
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         MATINDEX(4)=sD+LORB
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
          DO KORB = 0,nContC*nOrbCompC-1
           MATINDEX(3)=sC+KORB
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1
             MATINDEX(2)=sB+JORB
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               MATINDEX(1)=sA+IORB
               IELM = IELM+1
               MAT2(MATINDEX(DIMI),MATINDEX(DIMJ)) = &
                    & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
               !test just until this routine is testet thouroghly
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
   ENDDO
  ENDDO
 ENDIF
ENDDO

END SUBROUTINE Build_full_2dim_from_lstensor

!> \brief build full 3 dimensional array from an lstensor 
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the 3 dimensional array
!> \param dim1 the size of the 1. dimension
!> \param dim2 the size of the 2. dimension
!> \param dim3 the size of the 3. dimension
SUBROUTINE Build_full_3dim_from_lstensor(TENSOR,MAT2,dim1,dim2,dim3)
implicit none
INTEGER            :: dim1,dim2,dim3
TYPE(LSTENSOR)     :: TENSOR
real(realk)        :: MAT2(dim1,dim2,dim3)
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT
INTEGER    :: nbast(4),MATINDEX(4)

IF(TENSOR%primCStensor)CALL lsQUIT('this type of lstensor should not be used in this fashion4',-1)
nbast(1)=TENSOR%nbast1
nbast(2)=TENSOR%nbast2
nbast(3)=TENSOR%nbast3
nbast(4)=TENSOR%nbast4
IF(TENSOR%nmat.NE.dim3)CALL lsQUIT('called Build_full_3dim_from_lstensor with TENSOR%nmat.NE.dim3',-1)
DO I=1,4
   IF(dim1 .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(dim2 .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=dim1+dim2
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_matrix_from_lstensor dim do not match',-1)
do IMAT=1,dim3
 do J=1,dim2
  do I=1,dim1
   MAT2(I,J,IMAT) = 0.d0
  enddo
 enddo
enddo
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      IELM = 0
      DO IMAT = 1,TENSOR%nmat
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         MATINDEX(4)=sD+LORB
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
          DO KORB = 0,nContC*nOrbCompC-1
           MATINDEX(3)=sC+KORB
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1
             MATINDEX(2)=sB+JORB
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               MATINDEX(1)=sA+IORB
               IELM = IELM+1
               MAT2(MATINDEX(DIMI),MATINDEX(DIMJ),IMAT) = &
                    & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
               !test just until this routine is testet thouroghly
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
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDDO

END SUBROUTINE Build_full_3dim_from_lstensor

!> \brief build a reduced lstensor(number of AObatches instead of orbitals) from a full tensor
!> \author T. Kjaergaard
!> \date 2010
!> \param REDTENSOR the reduced lstensor
!> \param TENSOR the original lstensor
!> \param dim1 the size of the 1. dimension
!> \param dim2 the size of the 2. dimension
!> \param dim3 the size of the 3. dimension
!> \param dim4 the size of the 4. dimension
SUBROUTINE Build_reduce_lstensor(redtensor,TENSOR,dim1,dim2,dim3,dim4)
implicit none
INTEGEr          :: dim1,dim2,dim3,dim4
TYPE(LSTENSOR)   :: TENSOR
TYPE(LSTENSOR)   :: redTENSOR
!
INTEGER    :: I,J,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,testdim1,testdim2,testdim3,testdim4
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT,oldatom1,oldatom2,oldatom3,oldatom4
INTEGER :: startOrb1(TENSOR%natom1),startOrb2(TENSOR%natom2),startOrb3(TENSOR%natom3),startOrb4(TENSOR%natom4)
LOGICAL    :: ERROR
real(realk) :: maxElm
integer(kind=long) :: nmemsize

IF(TENSOR%primCStensor)CALL lsQUIT('this type of lstensor should not be used in this fashion5',-1)
REDTENSOR%primCStensor = TENSOR%primCStensor
REDTENSOR%gradienttensor = TENSOR%gradienttensor
redTENSOR%natom1 = TENSOR%natom1 
redTENSOR%natom2 = TENSOR%natom2 
redTENSOR%natom3 = TENSOR%natom3 
redTENSOR%natom4 = TENSOR%natom4 
redTENSOR%nbast1 = dim1
redTENSOR%nbast2 = dim2
redTENSOR%nbast3 = dim3
redTENSOR%nbast4 = dim4
redTENSOR%nmat   = 1
redTENSOR%nLSAO  = TENSOR%nLSAO

NULLIFY(redTENSOR%LSAO)
NULLIFY(redTENSOR%INDEX)
ALLOCATE(redTENSOR%LSAO(TENSOR%nLSAO))
ALLOCATE(redTENSOR%INDEX(TENSOR%natom1,TENSOR%natom2,TENSOR%natom3,TENSOR%natom4))
redTENSOR%INDEX = TENSOR%INDEX

testdim1 = 0
DO I=1,TENSOR%natom1 
   startOrb1(I) = testdim1
   J=TENSOR%INDEX(I,1,1,1)
   testdim1 = testdim1 + TENSOR%LSAO(J)%nAOBATCH1 
ENDDO
testdim2 = 0
DO I=1,TENSOR%natom2
   startOrb2(I) = testdim2
   J=TENSOR%INDEX(1,I,1,1)
   testdim2 = testdim2 + TENSOR%LSAO(J)%nAOBATCH2 
ENDDO
testdim3 = 0
DO I=1,TENSOR%natom3
   startOrb3(I) = testdim3
   J=TENSOR%INDEX(1,1,I,1)
   testdim3 = testdim3 + TENSOR%LSAO(J)%nAOBATCH3 
ENDDO
testdim4 = 0
DO I=1,TENSOR%natom4
   startOrb4(I) = testdim4
   J=TENSOR%INDEX(1,1,1,I)
   testdim4 = testdim4 + TENSOR%LSAO(J)%nAOBATCH4 
ENDDO
if(testdim1 .NE. dim1)call lsquit('error dim mismatch in Build_reduce_lstensor',-1)
if(testdim2 .NE. dim2)call lsquit('error dim mismatch in Build_reduce_lstensor',-1)
if(testdim3 .NE. dim3)call lsquit('error dim mismatch in Build_reduce_lstensor',-1)
if(testdim4 .NE. dim4)call lsquit('error dim mismatch in Build_reduce_lstensor',-1)

DO I = 1,TENSOR%nLSAO
 redTENSOR%LSAO(I)%ALLOC = TENSOR%LSAO(I)%ALLOC
 redTENSOR%LSAO(I)%ATOM1 = TENSOR%LSAO(I)%ATOM1 
 redTENSOR%LSAO(I)%ATOM2 = TENSOR%LSAO(I)%ATOM2 
 redTENSOR%LSAO(I)%ATOM3 = TENSOR%LSAO(I)%ATOM3 
 redTENSOR%LSAO(I)%ATOM4 = TENSOR%LSAO(I)%ATOM4 
 redTENSOR%LSAO(I)%nAOBATCH1 = TENSOR%LSAO(I)%nAOBATCH1 
 redTENSOR%LSAO(I)%nAOBATCH2 = TENSOR%LSAO(I)%nAOBATCH2 
 redTENSOR%LSAO(I)%nAOBATCH3 = TENSOR%LSAO(I)%nAOBATCH3
 redTENSOR%LSAO(I)%nAOBATCH4 = TENSOR%LSAO(I)%nAOBATCH4 
 testdim1 = startOrb1(TENSOR%LSAO(I)%ATOM1)
 testdim2 = startOrb2(TENSOR%LSAO(I)%ATOM2)
 testdim3 = startOrb3(TENSOR%LSAO(I)%ATOM3)
 testdim4 = startOrb4(TENSOR%LSAO(I)%ATOM4)

 IF(redTENSOR%LSAO(I)%ALLOC)THEN
    
  NULLIFY(redTENSOR%LSAO(I)%BATCH)
  ALLOCATE(redTENSOR%LSAO(I)%BATCH(TENSOR%LSAO(I)%nAOBATCH1 ,TENSOR%LSAO(I)%nAOBATCH2 ,&
       &TENSOR%LSAO(I)%nAOBATCH3 ,TENSOR%LSAO(I)%nAOBATCH4))

  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4

      maxElm = 0.d0
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(1) = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(1) = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(1) = testdim1+Ibat
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(1) = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(1) = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(1) = testdim2+Jbat
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(1) = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(1) = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(1) = testdim3+Kbat
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(1) = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(1) = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(1) = testdim4+Lbat
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE = 1
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms = 1
      NULLIFY(redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms)
      ALLOCATE(redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(1))
      
      IELM = 0
      DO IMAT = 1,TENSOR%nmat
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          DO KORB = 0,nContC*nOrbCompC-1
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1                           
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               IELM = IELM+1
               maxElm = MAX(maxElm,ABS(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)))
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      redTENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(1) = maxElm
     enddo
    enddo
   enddo
  enddo
 endif
enddo

Call Determine_lstensor_memory(redtensor,nmemsize)
call add_mem_to_global(nmemsize)

end SUBROUTINE Build_reduce_lstensor

!> \brief initiate the lstensor structure
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param AO1 the Atomic orbital(AO) item for center 1 
!> \param AO2 the AO item for center 2 
!> \param AO3 the AO item for center 3 
!> \param AO4 the AO item for center 4 
!> \param nbast1 the size of the 1. dimension
!> \param nbast2 the size of the 2. dimension
!> \param nbast3 the size of the 3. dimension
!> \param nbast4 the size of the 4. dimension
!> \param nmat the number of density matrices or size of the 5. dimension
!> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO3 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO4 flag to describe if the AO1 item should be used or an empty AO 
!> \param ODscreen flag: use the overlap distribution screening 
!> \param lupri the logical unit number for the output
SUBROUTINE init_lstensor_5dim(TENSOR,AO1,AO2,AO3,AO4,nbast1,nbast2,nbast3,nbast4,nmat,useAO1,useAO2,useAO3,useAO4,ODscreen,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM),target  :: AO1,AO2,AO3,AO4
TYPE(AOITEM),pointer :: AOT1,AOT2,AOT3,AOT4
INTEGER            :: nbast1,nbast2,nbast3,nbast4,nmat,lupri
logical   :: useAO1,useAO2,useAO3,useAO4,ODscreen
!
TYPE(AOITEM),target :: AOT
!REAL(REALK) :: MAXELM
INTEGER :: nElms
! BELONGINING TO AOT1 
INTEGER :: natom1,nbastI,nbatI,Ibat,Iangmom,Iatom,IORB
INTEGER :: nContA,nOrbCompA,sA,nAngA,AOT1batch,batAOT1
! BELONGINING TO AOT2 
INTEGER :: natom2,nbastJ,nbatJ,Jbat,Jangmom,Jatom,JORB
INTEGER :: nContB,nOrbCompB,sB,nAngB,AOT2batch,batAOT2
! BELONGINING TO AOT3 
INTEGER :: natom3,nbastK,nbatK,Kbat,Kangmom,Katom,KORB
INTEGER :: nContC,nOrbCompC,sC,nAngC,AOT3batch,batAOT3
! BELONGINING TO AOT4
INTEGER :: natom4,nbastL,nbatL,Lbat,Langmom,Latom,LORB
INTEGER :: nContD,nOrbCompD,sD,nAngD,AOT4batch,batAOT4
! COMMON
INTEGER :: DIM,IELM,IMAT2,I
LOGICAL :: ERROR,ScreenAB,ScreenCD,ScreenAtom
integer(kind=long) :: nmemsize
TENSOR%primCStensor = .FALSE.
TENSOR%gradienttensor = .FALSE.
call SET_EMPTY_AO(AOT)
AOT1 => AO1
AOT2 => AO2
AOT3 => AO3
AOT4 => AO4
IF(.NOT.useAO1) AOT1 => AOT
IF(.NOT.useAO2) AOT2 => AOT
IF(.NOT.useAO3) AOT3 => AOT
IF(.NOT.useAO4) AOT4 => AOT

!ERROR = .FALSE.
!IF(nbast1 .NE. AOT1%nbast)ERROR = .TRUE.
!IF(nbast2 .NE. AOT2%nbast)ERROR = .TRUE.
!IF(nbast3 .NE. AOT3%nbast)ERROR = .TRUE.
!IF(nbast4 .NE. AOT4%nbast)ERROR = .TRUE.
!IF(ERROR)THEN
!   print*,'input dim',nbast1,nbast2,nbast3,nbast4
!   print*,'AObatch  ',AOT1%nbast,AOT2%nbast,AOT3%nbast,AOT4%nbast
!   write(lupri,*)'input dim',nbast1,nbast2,nbast3,nbast4
!   write(lupri,*)'AObatch  ',AOT1%nbast,AOT2%nbast,AOT3%nbast,AOT4%nbast
!   CALL lsQUIT('input dimension and input AObatch mismatch in init_lstensor_5dim')
!ENDIF
natom1 = AOT1%natoms
natom2 = AOT2%natoms
natom3 = AOT3%natoms
natom4 = AOT4%natoms
NULLIFY(TENSOR%LSAO)
NULLIFY(TENSOR%INDEX)
ALLOCATE(TENSOR%LSAO(natom1*natom2*natom3*natom4))
ALLOCATE(TENSOR%INDEX(natom1,natom2,natom3,natom4))
TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
TENSOR%natom1 = natom1
TENSOR%natom2 = natom2
TENSOR%natom3 = natom3
TENSOR%natom4 = natom4
TENSOR%nbast1 = nbast1
TENSOR%nbast2 = nbast2
TENSOR%nbast3 = nbast3
TENSOR%nbast4 = nbast4
TENSOR%nmat = nmat
AOT1batch=0
I = 0
DO Iatom = 1,natom1
 nbatI = AOT1%ATOMICnBatch(IATOM)
 AOT2batch=0
 DO Jatom = 1,natom2
  nbatJ = AOT2%ATOMICnBatch(JATOM)
  AOT3batch=0
  DO Katom = 1,natom3
   nbatK = AOT3%ATOMICnBatch(KATOM)
   AOT4batch=0
   DO Latom = 1,natom4
    nbatL = AOT4%ATOMICnBatch(LATOM)
    I=I+1
    TENSOR%INDEX(IATOM,JATOM,KATOM,Latom) = I
    ScreenAtom = .FALSE.
    IF(ODscreen)call determineODscreening(nbatI,nbatJ,nbatK,nbatL,AOT1batch,&
         & AOT2batch,AOT3batch,AOT4batch,useAO1,useAO2,useAO3,&
         & useAO4,AOT1%BATCH,AOT2%BATCH,AOT3%BATCH,AOT4%BATCH,ScreenAtom)
    IF(ScreenAtom)THEN
     NULLIFY(TENSOR%LSAO(I)%BATCH)
     TENSOR%LSAO(I)%ALLOC = .FALSE.
     TENSOR%LSAO(I)%ATOM1 = Iatom
     TENSOR%LSAO(I)%ATOM2 = Jatom
     TENSOR%LSAO(I)%ATOM3 = Katom
     TENSOR%LSAO(I)%ATOM4 = Latom
     TENSOR%LSAO(I)%nAOBATCH1 = nbatI
     TENSOR%LSAO(I)%nAOBATCH2 = nbatJ
     TENSOR%LSAO(I)%nAOBATCH3 = nbatK
     TENSOR%LSAO(I)%nAOBATCH4 = nbatL
    ELSE
    NULLIFY(TENSOR%LSAO(I)%BATCH)
    ALLOCATE(TENSOR%LSAO(I)%BATCH(nbatI,nbatJ,nbatK,nbatL))
    TENSOR%LSAO(I)%ALLOC = .TRUE.
    TENSOR%LSAO(I)%ATOM1 = Iatom
    TENSOR%LSAO(I)%ATOM2 = Jatom
    TENSOR%LSAO(I)%ATOM3 = Katom
    TENSOR%LSAO(I)%ATOM4 = Latom
    TENSOR%LSAO(I)%nAOBATCH1 = nbatI
    TENSOR%LSAO(I)%nAOBATCH2 = nbatJ
    TENSOR%LSAO(I)%nAOBATCH3 = nbatK
    TENSOR%LSAO(I)%nAOBATCH4 = nbatL
    batAOT1 = AOT1batch
    DO Ibat = 1,nbatI
     batAOT1 = batAOT1+1 
     batAOT2 = AOT2batch
     DO Jbat = 1,nbatJ
      batAOT2 = batAOT2+1 
      screenAB=.FALSE.
      IF(ODscreen)call getODscreening(AOT1%BATCH(batAOT1),AOT2%BATCH(batAOT2),screenAB)!empty?
      IF(screenAB)THEN
       batAOT3 = AOT3batch
       DO Kbat = 1,nbatK
        batAOT3 = batAOT3+1 
        batAOT4 = AOT4batch
        DO Lbat = 1,nbatL
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE=0
         TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms=0
         NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms) 
        ENDDO
       ENDDO
      ELSE
       batAOT3 = AOT3batch
       DO Kbat = 1,nbatK
        batAOT3 = batAOT3+1 
        batAOT4 = AOT4batch
        DO Lbat = 1,nbatL
         batAOT4 = batAOT4+1
         screenCD = .FALSE.
         IF(ODscreen)call getODscreening(AOT3%BATCH(batAOT3),AOT4%BATCH(batAOT4),screenCD)!empty?
         IF(screenCD)THEN
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE=0
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms=0
            NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms) 
         ELSE
          nAngA = AOT1%BATCH(batAOT1)%nAngmom 
          nAngB = AOT2%BATCH(batAOT2)%nAngmom 
          nAngC = AOT3%BATCH(batAOT3)%nAngmom 
          nAngD = AOT4%BATCH(batAOT4)%nAngmom 
          TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA = nAngA
          TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB = nAngB
          TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC = nAngC
          TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD = nAngD
          DO Iangmom = 1,nAngA
           nContA = AOT1%BATCH(batAOT1)%nContracted(Iangmom)
           nOrbCompA = AOT1%BATCH(batAOT1)%nOrbComp(Iangmom)
           sA = AOT1%BATCH(batAOT1)%startOrbital(Iangmom)
           DO Jangmom = 1,nAngB
            nContB = AOT2%BATCH(batAOT2)%nContracted(Jangmom)
            nOrbCompB = AOT2%BATCH(batAOT2)%nOrbComp(Jangmom)
            sB = AOT2%BATCH(batAOT2)%startOrbital(Jangmom) 
            DO Kangmom = 1,nAngC
             nContC = AOT3%BATCH(batAOT3)%nContracted(Kangmom)
             nOrbCompC = AOT3%BATCH(batAOT3)%nOrbComp(Kangmom)
             sC = AOT3%BATCH(batAOT3)%startOrbital(Kangmom)
             DO Langmom = 1,nAngD
              nContD = AOT4%BATCH(batAOT4)%nContracted(Langmom)
              nOrbCompD = AOT4%BATCH(batAOT4)%nOrbComp(Langmom)
              sD = AOT4%BATCH(batAOT4)%startOrbital(Langmom)      
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom) = nContA
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom) = nOrbCompA
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom) = sA
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom) = nContB
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom) = nOrbCompB
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom) = sB
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom) = nContC
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom) = nOrbCompC
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom) = sC
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom) = nContD
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom) = nOrbCompD
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom) = sD
              TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE = 1
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          DIM = 0
          DO IMAT2 = 1,nmat
           DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
            nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
            nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
            DO LORB = 0,nContD*nOrbCompD-1
             DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
             nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
             nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
              DO KORB = 0,nContC*nOrbCompC-1
               DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
                nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
                nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
                DO JORB = 0,nContB*nOrbCompB-1
                 DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA 
                  nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
                  nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
                  DO IORB = 0,nContA*nOrbCompA-1
                   DIM = DIM+1
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms)
          ALLOCATE(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(DIM))
          TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms = DIM/nmat
          IELM = 0
          DO IMAT2 = 1,nmat
           DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
            nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
            nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
            sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
            DO LORB = 0,nContD*nOrbCompD-1
             DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
              nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
              nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
              sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)        
              DO KORB = 0,nContC*nOrbCompC-1
               DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
                nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
                nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
                sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
                DO JORB = 0,nContB*nOrbCompB-1                               
                 DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA 
                  nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
                  nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
                  sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
                  DO IORB = 0,nContA*nOrbCompA-1
                   IELM = IELM+1
                   IF(IELM .GT. DIM)CALL lsQUIT('Build_lstensor_from_full_5dim',-1)
                   TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM) = 0.d0
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
        ENDDO
       ENDDO
      ENDIF
     ENDDO
    ENDDO
    ENDIF
    AOT4batch = AOT4batch + nbatL
   ENDDO
   AOT3batch = AOT3batch + nbatK
  ENDDO
  AOT2batch = AOT2batch + nbatJ
 ENDDO
 AOT1batch = AOT1batch + nbatI
ENDDO
TENSOR%nLSAO = I

call FREE_EMPTY_AO(AOT)
Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE Init_lstensor_5dim

!> \brief build nuclear lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param fullmat the full nuclear screening matrix 
!> \param TENSOR the lstensor
!> \param natoms the number of atoms
SUBROUTINE build_Nuclearlstensor(fullmat,TENSOR,natoms)
implicit none
TYPE(LSTENSOR)     :: TENSOR
INTEGER            :: nAtoms
real(realk)        :: fullmat(natoms,1)
!
Integer :: Iatom
integer(kind=long) :: nmemsize

TENSOR%primCStensor = .FALSE.
TENSOR%gradienttensor = .FALSE.
NULLIFY(TENSOR%LSAO)
NULLIFY(TENSOR%INDEX)
ALLOCATE(TENSOR%LSAO(natoms))
ALLOCATE(TENSOR%INDEX(natoms,1,1,1))
TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
TENSOR%natom1 = natoms
TENSOR%natom2 = 1
TENSOR%natom3 = 1
TENSOR%natom4 = 1
TENSOR%nbast1 = 1
TENSOR%nbast2 = 1
TENSOR%nbast3 = 1
TENSOR%nbast4 = 1
TENSOR%nmat = 1
DO Iatom = 1,natoms
   TENSOR%INDEX(IATOM,1,1,1) = Iatom
   NULLIFY(TENSOR%LSAO(Iatom)%BATCH)
   ALLOCATE(TENSOR%LSAO(Iatom)%BATCH(1,1,1,1))
   TENSOR%LSAO(Iatom)%ALLOC = .TRUE.
   TENSOR%LSAO(Iatom)%ATOM1 = Iatom
   TENSOR%LSAO(Iatom)%ATOM2 = 1
   TENSOR%LSAO(Iatom)%ATOM3 = 1
   TENSOR%LSAO(Iatom)%ATOM4 = 1
   TENSOR%LSAO(Iatom)%nAOBATCH1 = 1
   TENSOR%LSAO(Iatom)%nAOBATCH2 = 1
   TENSOR%LSAO(Iatom)%nAOBATCH3 = 1
   TENSOR%LSAO(Iatom)%nAOBATCH4 = 1

   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nAngmomA = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nAngmomB = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nAngmomC = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nAngmomD = 1
            
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nContA(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nOrbCompA(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%startOrbA(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nContB(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nOrbCompB(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%startOrbB(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nContC(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nOrbCompC(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%startOrbC(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nContD(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nOrbCompD(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%startOrbD(1) = 1
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nelmE = 1
   NULLIFY(TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms)
   ALLOCATE(TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms(1))
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms(1) = fullmat(Iatom,1)
   TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%nelms = 1
ENDDO
TENSOR%nLSAO = natoms
Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE 

!> \brief initiate the gradient lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param natoms the number of atoms
!> \param nmat the number of matrices
!> \param lupri the logical unit number
SUBROUTINE init_gradientlstensor(TENSOR,natoms,nmat,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
INTEGER            :: nmat,lupri,natoms
!
INTEGER :: nElms
INTEGER :: natom1,IELM,IMAT2,I,J,Iatom
integer(kind=long) :: nmemsize

TENSOR%primCStensor = .FALSE.
TENSOR%gradienttensor = .TRUE.
natom1 = natoms
NULLIFY(TENSOR%LSAO)
NULLIFY(TENSOR%INDEX)
ALLOCATE(TENSOR%LSAO(natom1))
ALLOCATE(TENSOR%INDEX(natom1,1,1,1))
TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
TENSOR%natom1 = natom1
TENSOR%natom2 = 0
TENSOR%natom3 = 0
TENSOR%natom4 = 0
TENSOR%nbast1 = 0
TENSOR%nbast2 = 0
TENSOR%nbast3 = 0
TENSOR%nbast4 = 0
TENSOR%nmat = nmat
DO Iatom = 1,natom1
   I=iatom
   TENSOR%INDEX(I,1,1,1) = Iatom
   NULLIFY(TENSOR%LSAO(I)%BATCH)
   ALLOCATE(TENSOR%LSAO(I)%BATCH(1,1,1,1))
   TENSOR%LSAO(I)%ALLOC = .TRUE.
   TENSOR%LSAO(I)%ATOM1 = Iatom
   TENSOR%LSAO(I)%ATOM2 = 0
   TENSOR%LSAO(I)%ATOM3 = 0
   TENSOR%LSAO(I)%ATOM4 = 0
   TENSOR%LSAO(I)%nAOBATCH1 = 0
   TENSOR%LSAO(I)%nAOBATCH2 = 0
   TENSOR%LSAO(I)%nAOBATCH3 = 0
   TENSOR%LSAO(I)%nAOBATCH4 = 0

   TENSOR%LSAO(I)%BATCH(1,1,1,1)%nelmE = 3*nmat
   NULLIFY(TENSOR%LSAO(I)%BATCH(1,1,1,1)%elms)
   ALLOCATE(TENSOR%LSAO(I)%BATCH(1,1,1,1)%elms(3*nmat))
   DO J=1,3
      DO IMAT2 = 1,nmat
         TENSOR%LSAO(I)%BATCH(1,1,1,1)%elms(J+(IMAT2-1)*3) = 0.d0
      ENDDO
   ENDDO
ENDDO
TENSOR%nLSAO = natom1
Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE Init_gradientlstensor

!> \brief build a full gradient 2 array from gradient lstensor 
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param GRAD the 2dim gradient array
!> \param natom the number of atoms
!> \param nmat the number of matrices
!> \param lupri the logical unit number
SUBROUTINE build_gradient_from_gradientlstensor(TENSOR,GRAD,natom,nmat,lupri)
implicit none
TYPE(LSTENSOR)    :: TENSOR
REAL(REALK)       :: GRAD(3,natom*nmat) 
INTEGER           :: lupri
!
INTEGER :: nmat,nElms
INTEGER :: natom,IELM,IMAT,I,J,Iatom

natom = TENSOR%natom1
IF(TENSOR%nmat.EQ.1)THEN
   DO Iatom = 1,natom
      DO J=1,3
         GRAD(J,Iatom) = TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms(J)
      ENDDO
   ENDDO
ELSE
   DO Iatom = 1,natom
      DO J=1,3
         DO IMAT = 1,nmat
            GRAD(J,Iatom+(IMAT-1)*natom) = TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms(J+(IMAT-1)*3)
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE Build_gradient_from_gradientlstensor

!> \brief build a lstensor from full 5 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the full 5dim array
!> \param AO1 the Atomic orbital(AO) item for center 1 
!> \param AO2 the AO item for center 2 
!> \param AO3 the AO item for center 3 
!> \param AO4 the AO item for center 4 
!> \param nbast1 the size of the 1. dimension
!> \param nbast2 the size of the 2. dimension
!> \param nbast3 the size of the 3. dimension
!> \param nbast4 the size of the 4. dimension
!> \param nmat the number of density matrices or size of the 5. dimension
!> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO3 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO4 flag to describe if the AO1 item should be used or an empty AO 
!> \param lupri the logical unit number for the output
SUBROUTINE Build_lstensor_from_full_5dim(TENSOR,MAT2,AO1,AO2,AO3,AO4,nbast1,nbast2,nbast3,nbast4,&
     &nmat,useAO1,useAO2,useAO3,useAO4,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM),target  :: AO1,AO2,AO3,AO4
TYPE(AOITEM),pointer :: AOT1,AOT2,AOT3,AOT4
INTEGER            :: nbast1,nbast2,nbast3,nbast4,nmat,lupri
REAL(REALK)        :: MAT2(nbast1,nbast2,nbast3,nbast4,nmat)
logical   :: useAO1,useAO2,useAO3,useAO4
!
TYPE(AOITEM),target :: AOT
!REAL(REALK) :: MAXELM
INTEGER :: nElms
! BELONGINING TO AOT1 
INTEGER :: natom1,nbastI,nbatI,Ibat,Iangmom,Iatom,IORB
INTEGER :: nContA,nOrbCompA,sA,nAngA,AOT1batch,batAOT1
! BELONGINING TO AOT2 
INTEGER :: natom2,nbastJ,nbatJ,Jbat,Jangmom,Jatom,JORB
INTEGER :: nContB,nOrbCompB,sB,nAngB,AOT2batch,batAOT2
! BELONGINING TO AOT3 
INTEGER :: natom3,nbastK,nbatK,Kbat,Kangmom,Katom,KORB
INTEGER :: nContC,nOrbCompC,sC,nAngC,AOT3batch,batAOT3
! BELONGINING TO AOT4
INTEGER :: natom4,nbastL,nbatL,Lbat,Langmom,Latom,LORB
INTEGER :: nContD,nOrbCompD,sD,nAngD,AOT4batch,batAOT4
! COMMON
INTEGER :: DIM,IELM,IMAT2,I
integer(kind=long) :: nmemsize

TENSOR%primCStensor = .FALSE.
TENSOR%gradienttensor = .FALSE.

call SET_EMPTY_AO(AOT)
AOT1 => AO1
AOT2 => AO2
AOT3 => AO3
AOT4 => AO4
IF(.NOT.useAO1) AOT1 => AOT
IF(.NOT.useAO2) AOT2 => AOT
IF(.NOT.useAO3) AOT3 => AOT
IF(.NOT.useAO4) AOT4 => AOT

natom1 = AOT1%natoms
natom2 = AOT2%natoms
natom3 = AOT3%natoms
natom4 = AOT4%natoms
NULLIFY(TENSOR%LSAO)
NULLIFY(TENSOR%INDEX)
ALLOCATE(TENSOR%LSAO(natom1*natom2*natom3*natom4))
ALLOCATE(TENSOR%INDEX(natom1,natom2,natom3,natom4))
TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
TENSOR%natom1 = natom1
TENSOR%natom2 = natom2
TENSOR%natom3 = natom3
TENSOR%natom4 = natom4
TENSOR%nbast1 = nbast1
TENSOR%nbast2 = nbast2
TENSOR%nbast3 = nbast3
TENSOR%nbast4 = nbast4
TENSOR%nmat = nmat
AOT1batch=0
I = 0
DO Iatom = 1,natom1
 nbatI = AOT1%ATOMICnBatch(IATOM)
 AOT2batch=0
 DO Jatom = 1,natom2
  nbatJ = AOT2%ATOMICnBatch(JATOM)
  AOT3batch=0
  DO Katom = 1,natom3
   nbatK = AOT3%ATOMICnBatch(KATOM)
   AOT4batch=0
   DO Latom = 1,natom4
    nbatL = AOT4%ATOMICnBatch(LATOM)
    I=I+1
    TENSOR%INDEX(IATOM,JATOM,KATOM,Latom) = I
    NULLIFY(TENSOR%LSAO(I)%BATCH)
    ALLOCATE(TENSOR%LSAO(I)%BATCH(nbatI,nbatJ,nbatK,nbatL))
    TENSOR%LSAO(I)%ALLOC = .TRUE.
    TENSOR%LSAO(I)%ATOM1 = Iatom
    TENSOR%LSAO(I)%ATOM2 = Jatom
    TENSOR%LSAO(I)%ATOM3 = Katom
    TENSOR%LSAO(I)%ATOM4 = Latom
    TENSOR%LSAO(I)%nAOBATCH1 = nbatI
    TENSOR%LSAO(I)%nAOBATCH2 = nbatJ
    TENSOR%LSAO(I)%nAOBATCH3 = nbatK
    TENSOR%LSAO(I)%nAOBATCH4 = nbatL
    batAOT1 = AOT1batch
    DO Ibat = 1,nbatI
     batAOT1 = batAOT1+1 
     batAOT2 = AOT2batch
     DO Jbat = 1,nbatJ
      batAOT2 = batAOT2+1 
      batAOT3 = AOT3batch
      DO Kbat = 1,nbatK
       batAOT3 = batAOT3+1 
       batAOT4 = AOT4batch
       DO Lbat = 1,nbatL
        batAOT4 = batAOT4+1 
        nAngA = AOT1%BATCH(batAOT1)%nAngmom 
        nAngB = AOT2%BATCH(batAOT2)%nAngmom 
        nAngC = AOT3%BATCH(batAOT3)%nAngmom 
        nAngD = AOT4%BATCH(batAOT4)%nAngmom 
        TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA = nAngA
        TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB = nAngB
        TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC = nAngC
        TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD = nAngD
        DO Iangmom = 1,nAngA
         nContA = AOT1%BATCH(batAOT1)%nContracted(Iangmom)
         nOrbCompA = AOT1%BATCH(batAOT1)%nOrbComp(Iangmom)
         sA = AOT1%BATCH(batAOT1)%startOrbital(Iangmom)
         DO Jangmom = 1,nAngB
          nContB = AOT2%BATCH(batAOT2)%nContracted(Jangmom)
          nOrbCompB = AOT2%BATCH(batAOT2)%nOrbComp(Jangmom)
          sB = AOT2%BATCH(batAOT2)%startOrbital(Jangmom) 
          DO Kangmom = 1,nAngC
           nContC = AOT3%BATCH(batAOT3)%nContracted(Kangmom)
           nOrbCompC = AOT3%BATCH(batAOT3)%nOrbComp(Kangmom)
           sC = AOT3%BATCH(batAOT3)%startOrbital(Kangmom)
           DO Langmom = 1,nAngD
            nContD = AOT4%BATCH(batAOT4)%nContracted(Langmom)
            nOrbCompD = AOT4%BATCH(batAOT4)%nOrbComp(Langmom)
            sD = AOT4%BATCH(batAOT4)%startOrbital(Langmom)
            
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom) = nContA
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom) = nOrbCompA
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom) = sA
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom) = nContB
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom) = nOrbCompB
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom) = sB
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom) = nContC
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom) = nOrbCompC
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom) = sC
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom) = nContD
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom) = nOrbCompD
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom) = sD
            TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelmE = 1
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        DIM = 0
        DO IMAT2 = 1,nmat
         DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
          nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
          nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
          DO LORB = 0,nContD*nOrbCompD-1
           DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
            nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
            nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
            DO KORB = 0,nContC*nOrbCompC-1
             DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
              nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
              nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
              DO JORB = 0,nContB*nOrbCompB-1
               DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
                nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
                nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
                DO IORB = 0,nContA*nOrbCompA-1
                 DIM = DIM+1
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms)
        ALLOCATE(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(DIM))
        TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms = DIM/nmat
        IELM = 0
        DO IMAT2 = 1,nmat
         DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
          nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
          nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
          sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
          DO LORB = 0,nContD*nOrbCompD-1
           DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
            nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
            nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
            sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)        
            DO KORB = 0,nContC*nOrbCompC-1                           
             DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
              nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
              nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
              sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
              DO JORB = 0,nContB*nOrbCompB-1
               DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
                nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
                nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
                sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
                DO IORB = 0,nContA*nOrbCompA-1
                 IELM = IELM+1
                 IF(IELM .GT. DIM)CALL lsQUIT('Build_lstensor_from_full_5dim',-1)
                 TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM) &
                      & = MAT2(sA+IORB,sB+JORB,sC+KORB,sD+LORB,IMAT2)
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
      ENDDO
     ENDDO
    ENDDO
    AOT4batch = AOT4batch + nbatL
   ENDDO
   AOT3batch = AOT3batch + nbatK
  ENDDO
  AOT2batch = AOT2batch + nbatJ
 ENDDO
 AOT1batch = AOT1batch + nbatI
ENDDO
TENSOR%nLSAO = I

call FREE_EMPTY_AO(AOT)
Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE Build_lstensor_from_full_5dim

!> \brief build a lstensor from full 3 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the full 3dim array
!> \param AO1 the Atomic orbital(AO) item for center 1 
!> \param AO2 the AO item for center 2 
!> \param nbast1 the size of the 1. dimension
!> \param nbast2 the size of the 2. dimension
!> \param nmat the number of density matrices or size of the 5. dimension
!> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
!> \param lupri the logical unit number for the output
SUBROUTINE Build_lstensor_from_full_3dim(TENSOR,MAT2,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM),target :: AO1,AO2
INTEGER            :: nbast1,nbast2,nmat,lupri
REAL(REALK)        :: MAT2(nbast1,nbast2,nmat)
logical   :: useAO1,useAO2
!
logical   :: useAO3,useAO4
TYPE(AOITEM),pointer :: AOT1,AOT2,AOT3,AOT4
TYPE(AOITEM),target :: AOT

call SET_EMPTY_AO(AOT)
AOT1 => AO1
AOT2 => AO2
AOT3 => AOT
AOT4 => AOT
IF(.NOT.useAO1) AOT1 => AOT
IF(.NOT.useAO2) AOT2 => AOT
useAO3=.TRUE.
useAO4=.TRUE.
call Build_lstensor_from_full_5dim(TENSOR,MAT2,AOT1,AOT2,AOT3,AOT4,nbast1,nbast2,1,1,nmat,useAO1,useAO2,useAO3,useAO4,lupri)

call FREE_EMPTY_AO(AOT)

END SUBROUTINE Build_lstensor_from_full_3dim

!> \brief build a lstensor from full 2 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the full 3dim array
!> \param AO1 the Atomic orbital(AO) item for center 1 
!> \param AO2 the AO item for center 2 
!> \param nbast1 the size of the 1. dimension
!> \param nbast2 the size of the 2. dimension
!> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
!> \param lupri the logical unit number for the output
SUBROUTINE Build_lstensor_from_full_2dim(TENSOR,MAT2,AO1,AO2,nbast1,nbast2,useAO1,useAO2,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM),target :: AO1,AO2
INTEGER            :: nbast1,nbast2,lupri
REAL(REALK)        :: MAT2(nbast1,nbast2)
logical   :: useAO1,useAO2
!
logical   :: useAO3,useAO4
TYPE(AOITEM),pointer :: AOT1,AOT2,AOT3,AOT4
TYPE(AOITEM),target :: AOT
integer :: nmat
nmat=1
call SET_EMPTY_AO(AOT)
AOT1 => AO1
AOT2 => AO2
AOT3 => AOT
AOT4 => AOT
IF(.NOT.useAO1) AOT1 => AOT
IF(.NOT.useAO2) AOT2 => AOT
useAO3=.TRUE.
useAO4=.TRUE.
call Build_lstensor_from_full_5dim(TENSOR,MAT2,AOT1,AOT2,AOT3,AOT4,nbast1,nbast2,1,1,nmat,useAO1,useAO2,useAO3,useAO4,lupri)

call FREE_EMPTY_AO(AOT)

END SUBROUTINE Build_lstensor_from_full_2dim

!> \brief wrapper routine for building a single type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output
!> \param TENSOR the lstensor
!> \param MAT the type matrix
SUBROUTINE Build_singlematrix_from_lstensor(lupri,TENSOR,MAT)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT
integer :: lupri
!
real(realk),pointer  :: fullMAT(:,:,:,:,:)
integer              :: I,nbast(4)
!print*,'build_single matrix'
select case(matrix_type)
case(mtype_dense)
   call Build_single_dense_matrix_from_lstensor(TENSOR,MAT)
case(mtype_sparse1)
   call mem_alloc(fullMAT,TENSOR%nbast1,TENSOR%nbast2,TENSOR%nbast3,TENSOR%nbast4,1)
   call Build_full_5dim_from_lstensor(TENSOR,fullMAT,TENSOR%nbast1,TENSOR%nbast2,TENSOR%nbast3,TENSOR%nbast4,1)
   IF(TENSOR%nbast2 .EQ.1) call mat_set_from_full(fullmat(:,1,:,1,1),1.d0, MAT)
   IF(TENSOR%nbast3 .EQ.1) call mat_set_from_full(fullmat(:,:,1,1,1),1.d0, MAT)
   call mem_dealloc(fullMAT)
case(mtype_unres_dense)
   call Build_single_unres_matrix_from_lstensor(TENSOR,MAT)
case(mtype_sparse_block)
   call mem_alloc(fullMAT,TENSOR%nbast1,TENSOR%nbast2,TENSOR%nbast3,TENSOR%nbast4,1)
   call Build_full_5dim_from_lstensor(TENSOR,fullMAT,TENSOR%nbast1,TENSOR%nbast2,TENSOR%nbast3,TENSOR%nbast4,1)
   IF(TENSOR%nbast2 .EQ.1) call mat_set_from_full(fullmat(:,1,:,1,1),1.d0, MAT)
   IF(TENSOR%nbast3 .EQ.1) call mat_set_from_full(fullmat(:,:,1,1,1),1.d0, MAT)
   call mem_dealloc(fullMAT)
case(mtype_csr)
   call Build_single_csr_matrix_from_lstensor(TENSOR,MAT)
case default
   stop "BUILD_SINGLEMATRIX_FROM_LSTENSOR not implemented for this type of matrix"
end select
!print*,'done build_single matrix'
END SUBROUTINE BUILD_SINGLEMATRIX_FROM_LSTENSOR

!> \brief build a single type dense matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the type matrix
SUBROUTINE Build_single_dense_matrix_from_lstensor(TENSOR,MAT)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT
INTEGER    :: nbast(4),MATINDEX(4)

nbast(1)=TENSOR%nbast1
nbast(2)=TENSOR%nbast2
nbast(3)=TENSOR%nbast3
nbast(4)=TENSOR%nbast4

IF(TENSOR%nmat .NE.1)CALL lsQUIT('ERROR: Build_single_dense_matrix_from_lstensor nmat > 1',-1)
DO I=1,4
   IF(MAT%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT%nrow+MAT%ncol
!print*,'DIMI',DIMI,'DIMJ',DIMJ
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_single_dense_matrix_from_lstensor dim do not match',-1)

MAT%elms = 0.d0
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
!   print*,'Ibat=',Ibat
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      IELM = 0
      DO IMAT = 1,TENSOR%nmat
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         MATINDEX(4)=sD+LORB
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
          DO KORB = 0,nContC*nOrbCompC-1
           MATINDEX(3)=sC+KORB
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1
             MATINDEX(2)=sB+JORB                          
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               MATINDEX(1)=sA+IORB
               IELM = IELM+1
!               print*,'MATINDEX(DIMI)',MATINDEX(DIMI),'DIMI',DIMI
!               print*,'MATINDEX(DIMJ)',MATINDEX(DIMJ),'DIMJ',DIMJ
               MAT%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT%nrow) = &
                    & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
               !test just until this routine is testet thouroghly
               DO J=1,4
                IF((J.NE.DIMI).AND.(J.NE.DIMJ))THEN
                   IF(MATINDEX(J).NE.1)&
  &CALL lsQUIT('in Build_single_dense_matrix_from_lstensor lstensor seem to have more than 2 dim',-1)
                ENDIF
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
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDDO

END SUBROUTINE Build_single_dense_matrix_from_lstensor

!> \brief build a single csr (compress sparse row) type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT the type matrix
subroutine Build_single_csr_matrix_from_lstensor(TENSOR,MAT)
use matrix_operations_csr
!include 'mkl_spblas.fi'
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IMAT,IELMs
INTEGER    :: nbast(4),MATINDEX(4),job(8),n,nnz,info,l1,olda1,l2
INTEGER,pointer    :: ROW(:),COL(:),IELM(:)
real(realk),pointer :: VAL(:)

if(mat%ncol.NE.mat%nrow)then
   call lsquit('Build_single_csr_matrix_from_lstensor requires a quadratic matrix',-1)
endif

nbast(1)=TENSOR%nbast1
nbast(2)=TENSOR%nbast2
nbast(3)=TENSOR%nbast3
nbast(4)=TENSOR%nbast4

IF(TENSOR%nmat .NE.1)CALL lsQUIT('ERROR: Build_single_csr_matrix_from_lstensor nmat > 1',-1)
DO I=1,4
   IF(MAT%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT%nrow+MAT%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_single_csr_matrix_from_lstensor dim do not match',-1)

! PROPER WAY OF DOING THIS :
!
! ATTACH TO BASIS AND THEN TO AOBATCH A LOOPORDER FOR EACH TYPE
! ALSO ATTACH A TYPE FOR ALL ATOMS
! SO THAT I CAN DO
! do iatom = 1,tensor%natom1
!    I = tensor%index(Iatom,1,1,1)
!    ID = atomtype(iatom)
!    ib=1,length(order(ID)%batch(:))
!    ibat = order(ID)%batch(ib)
!    iang = order(ID)%ang(ib)
!    calc sA 
!    sA = 0
!    DO SAA=2,iang
!      sA = sA+nOrbComp(SAA-1)*TENSOR%LSAO(I)%BATCH(Ibat,1,1,1)%nContA(SAA-1)
!    ENDDO
!    nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,1,1,1)%nOrbCompA(iang)
!    startA = TENSOR%LSAO(I)%BATCH(Ibat,1,1,1)%startOrbA(Iang)
!    do iContA = 1,TENSOR%LSAO(I)%BATCH(Ibat,1,1,1)%nContA(iang)
!     do iOrbCompA = 1,
!      ielm = sA+iOrbCompA+(icontA)*nOrbCompA
!      orbitalindex = startA + iOrbCompA+(icontA)*nOrbCompA
!     enddo
!    enddo
! enddo
! this should mean that orbitalindex is now consequtive

! FOR NOW WE JUST DO THIS
NNZ = 0
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      IELMs = 0
      DO IMAT = 1,TENSOR%nmat
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          DO KORB = 0,nContC*nOrbCompC-1
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               IELMs = IELMs+1
               IF(ABS(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELMs)).GT.zeroCSR)NNZ=NNZ+1
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
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDDO

!call mat_csr_allocate(MAT,NNZ,MAT%nrow) Stinne, Rasmus 2/8-2010
call mat_csr_allocate(MAT,NNZ)

call mem_alloc(VAL,NNZ)
call mem_alloc(ROW,NNZ)
call mem_alloc(COL,NNZ)
NNZ = 0
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
   DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
      IELMs = 0
      DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
       nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
       nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
       sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
       DO LORB = 0,nContD*nOrbCompD-1
        MATINDEX(4)=sD+LORB
        DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
         nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
         nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
         sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
         DO KORB = 0,nContC*nOrbCompC-1
          MATINDEX(3)=sC+KORB
          DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
           nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
           nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
           sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
           DO JORB = 0,nContB*nOrbCompB-1
            MATINDEX(2)=sB+JORB                          
            DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
             nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
             nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
             sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
             DO IORB = 0,nContA*nOrbCompA-1
              MATINDEX(1)=sA+IORB
              IELMs = IELMs+1
              IF(ABS(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELMs)).GT.zeroCSR)THEN
                 NNZ=NNZ+1
                 VAL(NNZ) = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELMs) 
                 ROW(NNZ) = MATINDEX(DIMI) 
                 COL(NNZ) = MATINDEX(DIMJ) 
              ENDIF
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
   ENDDO
  ENDDO
 ENDIF
ENDDO

!WRITE(6,*)'THE COORDINATE FORM   nnz=',nnz
!WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(row(j),j=1,nnz)

call QSORT_csr(ROW,COL,VAL)

!WRITE(6,*)'AFTER INITIAL SORT AFTER FIRST nnz=',nnz
!WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(row(j),j=1,nnz)

L1=1
olda1 = ROW(L1)
L2=1
do I=2,NNZ
   IF(olda1.EQ.ROW(I))THEN
      L2=L2+1
   ELSE
      call QSORT_csr(COL(L1:L2),ROW(L1:L2),VAL(L1:L2))
!      WRITE(6,*)'AFTER SORT for ROW(I)=',olda1
!      WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
!      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
!      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(ROW(j),j=1,nnz)
      L1=L2+1
      L2=L2+1
      olda1 = ROW(L1)
   ENDIF
   IF(I.EQ.NNZ)THEN
      call QSORT_csr(COL(L1:NNZ),ROW(L1:NNZ),VAL(L1:NNZ))
!      WRITE(6,*)'AFTER FINAL SORT'
!      WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
!      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
!      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(ROW(j),j=1,nnz)
   ENDIF
enddo

job(1)=1
job(2)=1
job(3)=1
job(4)=1
job(5)=nnz
job(6)=0
job(7)=1
job(8)=1
n = mat%nrow
!print*,'MKL: n',n
!WRITE(6,*)'THE COORDINATE FORM   nnz=',nnz,'n',n
!WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(row(j),j=1,nnz)
#ifdef VAR_MKL
call mkl_dcsrcoo(job,n,mat%val,mat%col,mat%row,nnz,VAL,ROW,COL,info)
#endif
!HACK HACK
!!$MAT%elms = 0.d0
!!$DO I = 1,TENSOR%nLSAO
!!$ IF(TENSOR%LSAO(I)%ALLOC)THEN
!!$  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
!!$   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
!!$    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
!!$     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
!!$      IELMs = 0
!!$      DO IMAT = 1,TENSOR%nmat
!!$       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
!!$        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
!!$        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
!!$        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
!!$        DO LORB = 0,nContD*nOrbCompD-1
!!$         MATINDEX(4)=sD+LORB
!!$         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
!!$          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
!!$          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
!!$          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
!!$          DO KORB = 0,nContC*nOrbCompC-1
!!$           MATINDEX(3)=sC+KORB
!!$           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
!!$            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
!!$            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
!!$            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
!!$            DO JORB = 0,nContB*nOrbCompB-1
!!$             MATINDEX(2)=sB+JORB                          
!!$             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
!!$              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
!!$              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
!!$              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
!!$              DO IORB = 0,nContA*nOrbCompA-1
!!$               MATINDEX(1)=sA+IORB
!!$               IELMs = IELMs+1
!!$               MAT%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT%nrow) = &
!!$                    & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELMs)
!!$               !test just until this routine is testet thouroghly
!!$              ENDDO
!!$             ENDDO
!!$            ENDDO
!!$           ENDDO
!!$          ENDDO
!!$         ENDDO
!!$        ENDDO
!!$       ENDDO
!!$      ENDDO
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$  ENDDO
!!$ ENDIF
!!$ENDDO
!HACK HACK
!print*,'THE CSR from MKL  '
!call mat_print(mat,1,mat%nrow,1,mat%ncol,6)
call mem_dealloc(VAL)
call mem_dealloc(ROW)
call mem_dealloc(COL)

END subroutine Build_single_csr_matrix_from_lstensor

!> \brief build an array of csr (compress sparse row) type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT the array of type matrix
subroutine Build_array_csr_matrix_from_lstensor(TENSOR,MAT)
use matrix_operations_csr
!include 'mkl_spblas.fi'
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT(:)
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA,nmat
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IMAT,IELMs
INTEGER    :: nbast(4),MATINDEX(4),job(8),n,info,l1,olda1,l2,maxnnz
INTEGER,pointer    :: ROW(:,:),COL(:,:),IELM(:,:),nnz(:)
real(realk),pointer :: VAL(:,:)

if(mat(1)%ncol .NE. mat(1)%nrow)then
   call lsquit('Build_single_csr_matrix_from_lstensor requires a quadratic matrix',-1)
endif
nbast(1)=TENSOR%nbast1
nbast(2)=TENSOR%nbast2
nbast(3)=TENSOR%nbast3
nbast(4)=TENSOR%nbast4
nmat = TENSOR%nmat
allocate(nnz(nmat))
DO I=1,4
   IF(MAT(1)%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT(1)%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT(1)%nrow+MAT(1)%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_single_csr_matrix_from_lstensor dim do not match',-1)

! PROPER WAY OF DOING THIS :
!
! ATTACH TO BASIS AND THEN TO AOBATCH A LOOPORDER FOR EACH TYPE
! ALSO ATTACH A TYPE FOR ALL ATOMS
! SO THAT I CAN DO
! do iatom = 1,tensor%natom1
!    I = tensor%index(Iatom,1,1,1)
!    ID = atomtype(iatom)
!    ib=1,length(order(ID)%batch(:))
!    ibat = order(ID)%batch(ib)
!    iang = order(ID)%ang(ib)
!    calc sA 
!    sA = 0
!    DO SAA=2,iang
!      sA = sA+nOrbComp(SAA-1)*TENSOR%LSAO(I)%BATCH(Ibat,1,1,1)%nContA(SAA-1)
!    ENDDO
!    nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,1,1,1)%nOrbCompA(iang)
!    startA = TENSOR%LSAO(I)%BATCH(Ibat,1,1,1)%startOrbA(Iang)
!    do iContA = 1,TENSOR%LSAO(I)%BATCH(Ibat,1,1,1)%nContA(iang)
!     do iOrbCompA = 1,
!      ielm = sA+iOrbCompA+(icontA)*nOrbCompA
!      orbitalindex = startA + iOrbCompA+(icontA)*nOrbCompA
!     enddo
!    enddo
! enddo
! this should mean that orbitalindex is now consequtive

! FOR NOW WE JUST DO THIS
NNZ = 0
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      IELMs = 0
      DO IMAT = 1,TENSOR%nmat
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          DO KORB = 0,nContC*nOrbCompC-1
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               IELMs = IELMs+1
               IF(ABS(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELMs)).GT.zeroCSR)NNZ(IMAT)=NNZ(IMAT)+1
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
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDDO

!call mat_csr_allocate(MAT,NNZ,MAT%nrow) Stinne, Rasmus 2/8-2010
do IMAT=1,nmat
   call mat_csr_allocate(MAT(IMAT),NNZ(IMAT))
enddo
maxnnz = 0
do IMAT=1,nmat
   maxnnz = max(maxnnz,nnz(IMAT))
enddo

call mem_alloc(VAL,maxNNZ,NMAT)
call mem_alloc(ROW,maxNNZ,NMAT)
call mem_alloc(COL,maxNNZ,NMAT)
NNZ = 0
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
   DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
      IELMs = 0
      DO IMAT = 1,TENSOR%nmat
      DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
       nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
       nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
       sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
       DO LORB = 0,nContD*nOrbCompD-1
        MATINDEX(4)=sD+LORB
        DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
         nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
         nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
         sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
         DO KORB = 0,nContC*nOrbCompC-1
          MATINDEX(3)=sC+KORB
          DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
           nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
           nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
           sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
           DO JORB = 0,nContB*nOrbCompB-1
            MATINDEX(2)=sB+JORB                          
            DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
             nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
             nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
             sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
             DO IORB = 0,nContA*nOrbCompA-1
              MATINDEX(1)=sA+IORB
              IELMs = IELMs+1
              IF(ABS(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELMs)).GT.zeroCSR)THEN
                 NNZ(IMAT)=NNZ(IMAT)+1
                 VAL(NNZ(IMAT),IMAT) = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELMs) 
                 ROW(NNZ(IMAT),IMAT) = MATINDEX(DIMI) 
                 COL(NNZ(IMAT),IMAT) = MATINDEX(DIMJ) 
              ENDIF
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
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDDO

!DO IMAT = 1,TENSOR%nmat
!WRITE(6,*)'THE COORDINATE FORM   nnz=',nnz(IMAT),'MATRIX NR=',I
!WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j,IMAT),j=1,nnz(IMAT))
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j,IMAT),j=1,nnz(IMAT))
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(row(j,IMAT),j=1,nnz(IMAT))
!ENDDO

DO IMAT = 1,TENSOR%nmat
   call QSORT_csr(ROW(1:NNZ(IMAT),IMAT),COL(1:NNZ(IMAT),IMAT),VAL(1:NNZ(IMAT),IMAT))

!WRITE(6,*)'AFTER INITIAL SORT AFTER FIRST nnz=',nnz
!WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(row(j),j=1,nnz)

   L1=1
   olda1 = ROW(L1,IMAT)
   L2=1
   do I=2,NNZ(IMAT)
      IF(olda1.EQ.ROW(I,IMAT))THEN
         L2=L2+1
      ELSE
         call QSORT_csr(COL(L1:L2,IMAT),ROW(L1:L2,IMAT),VAL(L1:L2,IMAT))
         !      WRITE(6,*)'AFTER SORT for ROW(I)=',olda1
         !      WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
         !      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
         !      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(ROW(j),j=1,nnz)
         L1=L2+1
         L2=L2+1
         olda1 = ROW(L1,IMAT)
      ENDIF
      IF(I.EQ.NNZ(IMAT))THEN
         call QSORT_csr(COL(L1:NNZ(IMAT),IMAT),ROW(L1:NNZ(IMAT),IMAT),VAL(L1:NNZ(IMAT),IMAT))
         !      WRITE(6,*)'AFTER FINAL SORT'
         !      WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
         !      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
         !      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(ROW(j),j=1,nnz)
      ENDIF
   enddo

   job(1)=1
   job(2)=1
   job(3)=1
   job(4)=1
   job(5)=nnz(IMAT)
   job(6)=0
   job(7)=1
   job(8)=1
   n = mat(IMAT)%nrow
   !call mat_print(mat,1,mat%nrow,1,mat%nrow,6)
#ifdef VAR_MKL
   call mkl_dcsrcoo(job,n,mat(IMAT)%val,mat(IMAT)%col,mat(IMAT)%row,nnz(IMAT),VAL(1:NNZ(IMAT),IMAT),ROW(1:NNZ(IMAT),IMAT),COL(1:NNZ(IMAT),IMAT),info)
#endif

enddo
!HACK HACK
!!$MAT%elms = 0.d0
!!$DO I = 1,TENSOR%nLSAO
!!$ IF(TENSOR%LSAO(I)%ALLOC)THEN
!!$  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
!!$   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
!!$    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
!!$     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
!!$      IELMs = 0
!!$      DO IMAT = 1,TENSOR%nmat
!!$       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
!!$        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
!!$        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
!!$        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
!!$        DO LORB = 0,nContD*nOrbCompD-1
!!$         MATINDEX(4)=sD+LORB
!!$         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
!!$          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
!!$          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
!!$          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
!!$          DO KORB = 0,nContC*nOrbCompC-1
!!$           MATINDEX(3)=sC+KORB
!!$           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
!!$            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
!!$            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
!!$            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
!!$            DO JORB = 0,nContB*nOrbCompB-1
!!$             MATINDEX(2)=sB+JORB                          
!!$             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
!!$              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
!!$              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
!!$              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
!!$              DO IORB = 0,nContA*nOrbCompA-1
!!$               MATINDEX(1)=sA+IORB
!!$               IELMs = IELMs+1
!!$               MAT%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT%nrow) = &
!!$                    & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELMs)
!!$               !test just until this routine is testet thouroghly
!!$              ENDDO
!!$             ENDDO
!!$            ENDDO
!!$           ENDDO
!!$          ENDDO
!!$         ENDDO
!!$        ENDDO
!!$       ENDDO
!!$      ENDDO
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$  ENDDO
!!$ ENDIF
!!$ENDDO
!HACK HACK
!print*,'THE CSR from MKL  '
!call mat_print(mat(IMAT),1,mat(IMAT)%nrow,1,mat(IMAT)%nrow,6)
deallocate(nnz)
call mem_dealloc(VAL)
call mem_dealloc(ROW)
call mem_dealloc(COL)

END subroutine Build_array_csr_matrix_from_lstensor

!TEMP
!!$ THIS WILL SORT BIGGEST FIRSYT
!!$RECURSIVE SUBROUTINE Qsort_csr(a1,a2,val)
!!$IMPLICIT NONE
!!$INTEGER, INTENT(INOUT) :: a1(:),a2(:)
!!$REAL(REALK) :: VAL(:)
!!$INTEGER :: split
!!$INTEGER :: I
!!$
!!$  IF(size(a1) > 7) THEN 
!!$     CALL Partition_csr(a1,a2,VAL,split)
!!$     CALL Qsort_csr(a1(:split-1),a2(:split-1),VAL(:split-1))
!!$     CALL Qsort_csr(a1(split:),a2(split:),VAL(split:))
!!$  ELSEIF(size(a1) > 1)THEN
!!$     CALL INSERTION_csr(a1,a2,VAL)
!!$  END IF
!!$ 
!!$END SUBROUTINE Qsort_csr
!!$
!!$SUBROUTINE INSERTION_csr(a,a2,INDEXES)
!!$IMPLICIT NONE
!!$INTEGER, INTENT(INOUT) :: a(:),a2(:)
!!$INTEGER                    :: I,J,temp,tempI2
!!$real(realk) :: INDEXES(:),tempI
!!$
!!$Do I = 2, Size(a)
!!$   temp = a(I)
!!$   tempI = INDEXES(I)
!!$   tempI2 = a2(I)
!!$   DO J = I-1, 1, -1
!!$      IF(temp .GT. a(J)) Then
!!$         a(J+1) = a(J)
!!$         INDEXES(J+1) = INDEXES(J)
!!$         a2(J+1) = a2(J)
!!$      ELSE
!!$         Exit
!!$      ENDIF
!!$   ENDDO
!!$   a(J+1) = temp
!!$   INDEXES(J+1) = tempI
!!$   a2(J+1) = tempI2
!!$End Do
!!$
!!$END SUBROUTINE INSERTION_csr
!!$
!!$SUBROUTINE Partition_csr(a, a2,INDEXES, marker)
!!$IMPLICIT NONE
!!$INTEGER, INTENT(INOUT) :: a(:),a2(:)
!!$INTEGER, INTENT(OUT)       :: marker
!!$INTEGER                    :: left, right, I, temp,pivot,tempI2
!!$REAL(REALK)                :: INDEXES(:),tempI
!!$ 
!!$pivot = (a(1) + a(size(a))) / 2 ! Average of first and last elements 
!!$                                !to prevent quadratic behavior with
!!$                                ! sorted or reverse sorted data
!!$IF(a(1) .EQ. pivot)THEN
!!$   DO I =size(a)-1,2,-1
!!$      IF(a(I) .NE. pivot)THEN
!!$         pivot = (a(1) + a(I)) / 2
!!$         EXIT
!!$      ENDIF
!!$   ENDDO
!!$ENDIF
!!$
!!$left = 0                         ! behavior with sorted or reverse sorted data
!!$right = size(a) + 1
!!$
!!$DO WHILE (left < right)
!!$   right = right - 1
!!$   DO WHILE (a(right) < pivot)
!!$      right = right-1
!!$   END DO
!!$   left = left + 1
!!$   DO WHILE (a(left) > pivot)
!!$      left = left + 1
!!$   END DO
!!$   IF (left < right) THEN 
!!$      temp = a(right)
!!$      a(right) = a(left)
!!$      a(left) = temp
!!$      tempI = INDEXES(right)
!!$      tempI2 = a2(right)
!!$      INDEXES(right) = INDEXES(left)
!!$      a2(right) = a2(left)
!!$      INDEXES(left) = tempI
!!$      a2(left) = tempI2
!!$   END IF
!!$END DO
!!$
!!$IF (left == right) THEN
!!$   marker = left + 1
!!$ELSE
!!$   marker = left
!!$END IF
!!$
!!$END SUBROUTINE Partition_csr

!> \brief recursive quick sort algorithm, smallest first
!> \author T. Kjaergaard
!> \date 2010
!> \param a1 the values to be sorted
!> \param a2 and array that should be sorted like a1
!> \param val and array that should be sorted like a1
RECURSIVE SUBROUTINE Qsort_csr(a1,a2,val)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: a1(:),a2(:)
REAL(REALK) :: VAL(:)
INTEGER :: split
INTEGER :: I

  IF(size(a1) > 7) THEN 
     CALL Partition_csr(a1,a2,VAL,split)
     CALL Qsort_csr(a1(:split-1),a2(:split-1),VAL(:split-1))
     CALL Qsort_csr(a1(split:),a2(split:),VAL(split:))
  ELSEIF(size(a1) > 1)THEN
     CALL INSERTION_csr(a1,a2,VAL)
  END IF
 
END SUBROUTINE Qsort_csr

!> \brief insertion sort algorithm, smallest first
!> \author T. Kjaergaard
!> \date 2010
!> \param a the values to be sorted
!> \param a2 and array that should be sorted like a1
!> \param INDEXES and array that should be sorted like a1
SUBROUTINE INSERTION_csr(a,a2,INDEXES)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: a(:),a2(:)
INTEGER                    :: I,J,temp,tempI2
real(realk) :: INDEXES(:),tempI

Do I = 2, Size(a)
   temp = a(I)
   tempI = INDEXES(I)
   tempI2 = a2(I)
   DO J = I-1, 1, -1
      IF(temp .LT. a(J)) Then
         a(J+1) = a(J)
         INDEXES(J+1) = INDEXES(J)
         a2(J+1) = a2(J)
      ELSE
         Exit
      ENDIF
   ENDDO
   a(J+1) = temp
   INDEXES(J+1) = tempI
   a2(J+1) = tempI2
End Do

END SUBROUTINE INSERTION_csr

!> \brief the worker routine used by the qucik sort (qsort_csr) routine, smallest first
!> \author T. Kjaergaard
!> \date 2010
!> \param a the values to be sorted
!> \param a2 and array that should be sorted like a1
!> \param INDEXES and array that should be sorted like a1
!> \param marker counter index
SUBROUTINE Partition_csr(a, a2,INDEXES, marker)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: a(:),a2(:)
INTEGER, INTENT(OUT)       :: marker
INTEGER                    :: left, right, I, temp,pivot,tempI2
REAL(REALK)                :: INDEXES(:),tempI
 
pivot = (a(1) + a(size(a))) / 2 ! Average of first and last elements 
                                !to prevent quadratic behavior with
                                ! sorted or reverse sorted data
IF(a(1) .EQ. pivot)THEN
   DO I =size(a)-1,2,-1
      IF(a(I) .NE. pivot)THEN
         pivot = (a(1) + a(I)) / 2
         EXIT
      ENDIF
   ENDDO
ENDIF

left = 0                         ! behavior with sorted or reverse sorted data
right = size(a) + 1

DO WHILE (left < right)
   right = right - 1
   DO WHILE (a(right) > pivot)
      right = right-1
   END DO
   left = left + 1
   DO WHILE (a(left) < pivot)
      left = left + 1
   END DO
   IF (left < right) THEN 
      temp = a(left)
      a(left) = a(right)
      a(right) = temp
      tempI = INDEXES(left)
      tempI2 = a2(left)
      INDEXES(left) = INDEXES(right)
      a2(left) = a2(right)
      INDEXES(right) = tempI
      a2(right) = tempI2
   END IF
END DO

IF (left == right) THEN
   marker = left + 1
ELSE
   marker = left
END IF

END SUBROUTINE Partition_csr

!> \brief build a single unres_dense type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT the array of type matrix
SUBROUTINE Build_single_unres_matrix_from_lstensor(TENSOR,MAT)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT
INTEGER    :: nbast(4),MATINDEX(4),nElms
logical    :: option1

nbast(1)=TENSOR%nbast1
nbast(2)=TENSOR%nbast2
nbast(3)=TENSOR%nbast3
nbast(4)=TENSOR%nbast4
option1 = .true. !same alpha beta part
IF(TENSOR%nmat .EQ.2)option1 = .false. !different alpha beta part
IF(TENSOR%nmat .GT.2)call lsQUIT('Error: Build_single_unres_matrix_from_lstensor nmat > 2',-1)
DO I=1,4
   IF(MAT%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT%nrow+MAT%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_single_unres_matrix_from_lstensor dim do not match',-1)

MAT%elms = 0.d0
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      nELms = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms
      IELM = 0
      DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
       nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
       nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
       sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
       DO LORB = 0,nContD*nOrbCompD-1
        MATINDEX(4)=sD+LORB
        DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
         nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
         nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
         sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
         DO KORB = 0,nContC*nOrbCompC-1
          MATINDEX(3)=sC+KORB
          DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
           nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
           nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
           sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
           DO JORB = 0,nContB*nOrbCompB-1
            MATINDEX(2)=sB+JORB                          
            DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
             nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
             nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
             sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
             DO IORB = 0,nContA*nOrbCompA-1
              MATINDEX(1)=sA+IORB
              IELM = IELM+1
              if(option1)then
                 MAT%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT%nrow) = &
                      & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
                 MAT%elmsb(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT%nrow) = &
                      & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
              else
                 MAT%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT%nrow) = &
                      & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
                 MAT%elmsb(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT%nrow) = &
                      & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM+nElms)
              endif
              !test just until this routine is testet thouroghly
              DO J=1,4
               IF((J.NE.DIMI).AND.(J.NE.DIMJ))THEN
                IF(MATINDEX(J).NE.1)&
                     &CALL lsQUIT('in Build_single_unres_matrix_from_lstensor lstensor seem to have more than 2 dim',-1)
               ENDIF
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
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDDO

END SUBROUTINE Build_single_unres_matrix_from_lstensor

!> \brief wrapper routine for building an array of type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number of the output file
!> \param TENSOR the lstensor
!> \param MAT the array of type matrix
SUBROUTINE Build_matrixarray_from_lstensor(lupri,TENSOR,MAT)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT(:)
integer            :: lupri
!
real(realk),pointer  :: fullMAT(:,:,:,:,:)
integer              :: I
!print*,'build_matrixarray matrix'

select case(matrix_type)
case(mtype_dense)
   call Build_dense_matrixarray_from_lstensor(TENSOR,MAT)
case(mtype_sparse1)
   call mem_alloc(fullMAT,TENSOR%nbast1,TENSOR%nbast2,TENSOR%nbast3,TENSOR%nbast4,TENSOR%nmat)
   call Build_full_5dim_from_lstensor(TENSOR,fullMAT,TENSOR%nbast1,TENSOR%nbast2,TENSOR%nbast3,TENSOR%nbast4,TENSOR%nmat)
   IF(TENSOR%nbast2 .EQ.1) THEN
      do I=1,TENSOR%nmat
         call mat_set_from_full(fullmat(:,1,:,1,I),1.d0, MAT(I))
      enddo
   ENDIF
   IF(TENSOR%nbast3 .EQ.1) THEN
      do I=1,TENSOR%nmat
         call mat_set_from_full(fullmat(:,:,1,1,I),1.d0, MAT(I))
      enddo
   ENDIF
   call mem_dealloc(fullMAT)
case(mtype_unres_dense)
   call Build_unres_matrixarray_from_lstensor(TENSOR,MAT)
case(mtype_sparse_block)
   call mem_alloc(fullMAT,TENSOR%nbast1,TENSOR%nbast2,TENSOR%nbast3,TENSOR%nbast4,TENSOR%nmat)
   call Build_full_5dim_from_lstensor(TENSOR,fullMAT,TENSOR%nbast1,TENSOR%nbast2,TENSOR%nbast3,TENSOR%nbast4,TENSOR%nmat)
   IF(TENSOR%nbast2 .EQ.1) THEN
      do I=1,TENSOR%nmat
         call mat_set_from_full(fullmat(:,1,:,1,I),1.d0, MAT(I))
      enddo
   ENDIF
   IF(TENSOR%nbast3 .EQ.1) THEN
      do I=1,TENSOR%nmat
         call mat_set_from_full(fullmat(:,:,1,1,I),1.d0, MAT(I))
      enddo
   ENDIF
   call mem_dealloc(fullMAT)
case(mtype_csr)
   call Build_array_csr_matrix_from_lstensor(TENSOR,MAT)
case default
   stop "BUILD_MATRIXARRAY_FROM_LSTENSOR not implemented for this type of matrix"
end select
!print*,'done build_matrixarray matrix'

END SUBROUTINE BUILD_MATRIXARRAY_FROM_LSTENSOR

!> \brief build an array of dense type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT273 the array of type matrix
SUBROUTINE Build_dense_matrixarray_from_lstensor(TENSOR,MAT273)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT273(:)
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT
INTEGER    :: nbast(4),MATINDEX(4)

nbast(1)=TENSOR%nbast1
nbast(2)=TENSOR%nbast2
nbast(3)=TENSOR%nbast3
nbast(4)=TENSOR%nbast4

DO I=1,4
   IF(MAT273(1)%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT273(1)%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT273(1)%nrow+MAT273(1)%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_dense_matrixarray_from_lstensor dim do not match',-1)

DO I = 1,TENSOR%nmat
   MAT273(I)%elms = 0.d0
ENDDO
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      IELM = 0
      DO IMAT = 1,TENSOR%nmat
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         MATINDEX(4)=sD+LORB
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
          DO KORB = 0,nContC*nOrbCompC-1
           MATINDEX(3)=sC+KORB
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1
             MATINDEX(2)=sB+JORB
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               MATINDEX(1)=sA+IORB
               IELM = IELM+1
               MAT273(IMAT)%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT273%nrow) = &
                    & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
               !test just until this routine is testet thouroghly
               DO J=1,4
                  IF((J.NE.DIMI).AND.(J.NE.DIMJ))THEN
 IF(MATINDEX(J).NE.1)CALL lsQUIT('in Build_dense_matrixarray_from_lstensor lstensor seem to have more than 2 dim',-1)
                  ENDIF
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
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDDO

END SUBROUTINE Build_dense_matrixarray_from_lstensor

!> \brief build an array of unres_dense type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT273 the array of type matrix
SUBROUTINE Build_unres_matrixarray_from_lstensor(TENSOR,MAT273)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT273(:)
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT
INTEGER    :: nbast(4),MATINDEX(4),nElms,nmat
logical    :: option1

nbast(1)=TENSOR%nbast1
nbast(2)=TENSOR%nbast2
nbast(3)=TENSOR%nbast3
nbast(4)=TENSOR%nbast4
option1 = .true. !same alpha beta part
IF(TENSOR%nmat .EQ. (2*size(MAT273)))option1 = .false. !different alpha beta part
!IF(TENSOR%nmat .GT.2) call lsQUIT('Error: Build_unres_matrixarray_from_lstensor nmat > 2',-1)

IF (.NOT.option1) THEN
  nmat = TENSOR%nmat/2
ELSE
  nmat = TENSOR%nmat
ENDIF

DO I=1,4
   IF(MAT273(1)%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT273(1)%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT273(1)%nrow+MAT273(1)%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_unres_matrixarray_from_lstensor dim do not match',-1)

DO I = 1,nmat
   MAT273(I)%elms = 0.d0
   MAT273(I)%elmsb = 0.d0
ENDDO
DO I = 1,TENSOR%nLSAO
 IF(TENSOR%LSAO(I)%ALLOC)THEN
  DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
   DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
    DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
     DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
      nELms = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms
      IELM = 0
      DO IMAT = 1,nmat
       IF (.NOT.option1) IELM = 0
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         MATINDEX(4)=sD+LORB
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
          DO KORB = 0,nContC*nOrbCompC-1
           MATINDEX(3)=sC+KORB
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1
             MATINDEX(2)=sB+JORB                          
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               MATINDEX(1)=sA+IORB
               IELM = IELM+1
               if(option1)then
               ! Coulomb type contributions
                 MAT273(IMAT)%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT273%nrow) = &
                      & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
                 MAT273(IMAT)%elmsb(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT273%nrow) = &
                      & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM)
               else
               ! Exchange type contributions
                 MAT273(IMAT)%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT273%nrow) = &
                      & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM+nElms*2*(IMAT-1))
                 MAT273(IMAT)%elmsb(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*MAT273%nrow) = &
                      & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM+nElms*(2*IMAT-1))
               endif
               !test just until this routine is testet thouroghly
               DO J=1,4
                  IF((J.NE.DIMI).AND.(J.NE.DIMJ))THEN
 IF(MATINDEX(J).NE.1)CALL lsQUIT('in Build_unres_matrixarray_from_lstensor lstensor seem to have more than 2 dim',-1)
                  ENDIF
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
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDDO

END SUBROUTINE Build_unres_matrixarray_from_lstensor

!> \brief add a full 2 dimensional array to an lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param fullmat the array to be added to lstensor
!> \param dim1 the size of the 1. dimension 
!> \param dim2 the size of the 2. dimension 
!> \param imat not used 
SUBROUTINE add_full_2dim_to_lstensor(TENSOR,fullMAT,dim1,dim2,imat)
implicit none
TYPE(LSTENSOR)     :: TENSOR
real(realk)        :: fullMAT(:,:)
integer            :: dim1,dim2,imat
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,iatom,natom,nmat
INTEGER    :: nbast(4),MATINDEX(4)

nbast(1)=TENSOR%nbast1
nbast(2)=TENSOR%nbast2
nbast(3)=TENSOR%nbast3
nbast(4)=TENSOR%nbast4

DO I=1,4
   IF(dim1 .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(dim2 .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=dim1+dim2
IF(SUM1.NE.SUM2)THEN
   print*,'dim1',DIMI
   print*,'dim2',DIMJ
   print*,'dim1',dim1
   print*,'dim2',dim2
   print*,'nbast1',nbast(1)
   print*,'nbast2',nbast(2)
   print*,'nbast3',nbast(3)
   print*,'nbast4',nbast(4)
   print*,'nbast',nbast(DIMI)
   print*,'nbast',nbast(DIMJ)
   CALL lsQUIT('ERROR: add_full_2dim_to_lstensor dim do not match',-1)
ENDIF
IF(TENSOR%gradienttensor)THEN
 natom = TENSOR%natom1
 nmat = TENSOR%nmat
 DO Iatom = 1,natom
  DO IMAT = 1,nmat
   DO J=1,3
    TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms(J+(IMAT-1)*3) = &
         &TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms(J+(IMAT-1)*3) + fullmat(J,Iatom+(IMAT-1)*natom)
   ENDDO
  ENDDO
 ENDDO
ELSE
 DO I = 1,TENSOR%nLSAO
  IF(TENSOR%LSAO(I)%ALLOC)THEN
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
      DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
       IELM = 0                     
       DO Langmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomD
        nContD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContD(Langmom)
        nOrbCompD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompD(Langmom)
        sD = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbD(Langmom)
        DO LORB = 0,nContD*nOrbCompD-1
         MATINDEX(4)=sD+LORB
         DO Kangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomC
          nContC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContC(Kangmom)
          nOrbCompC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompC(Kangmom)
          sC = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbC(Kangmom)    
          DO KORB = 0,nContC*nOrbCompC-1
           MATINDEX(3)=sC+KORB
           DO Jangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomB
            nContB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContB(Jangmom)
            nOrbCompB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompB(Jangmom)
            sB = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbB(Jangmom)
            DO JORB = 0,nContB*nOrbCompB-1
             MATINDEX(2)=sB+JORB
             DO Iangmom = 1,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nAngmomA
              nContA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nContA(Iangmom)
              nOrbCompA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nOrbCompA(Iangmom)
              sA = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%startOrbA(Iangmom)
              DO IORB = 0,nContA*nOrbCompA-1
               MATINDEX(1)=sA+IORB
               IELM = IELM+1
               TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM) = &
                    &TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(IELM) + fullmat(MATINDEX(DIMI),MATINDEX(DIMJ))
               !test just until this routine is testet thouroghly
               DO J=1,4
                  IF((J.NE.DIMI).AND.(J.NE.DIMJ))THEN
                     IF(MATINDEX(J).NE.1)&
                          &CALL lsQUIT('in add_full_2dim_to_lstensor lstensor seem to have more than 2 dim',-1)
                  ENDIF
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
     ENDDO
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ENDIF

END SUBROUTINE add_full_2dim_to_lstensor

!> \brief build lstensor from single dense matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param AOT1 the 1. Atomic orbital
!> \param AOT2 the 2. Atomic orbital
!> \param mat the matrix type
SUBROUTINE Build_lstensor_from_dense_mat_single(TENSOR,AOT1,AOT2,MAT)
  use Matrix_module
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM)       :: AOT1,AOT2
TYPE(MATRIX)       :: MAT
!
REAL(REALK) :: MAXELM
INTEGER :: ATOMICINDEX,nElms
INTEGER :: nbast1,nbast2,natom1,natom2,Iatom,nbastI,nbatI
INTEGER :: nbastJ,nbatJ,Ibat,Jbat,iangmom,jangmom,jatom
INTEGER :: nContA,nOrbCompA,nContB,nOrbCompB,sA,sB,I,J 
integer(kind=long) :: nmemsize

TENSOR%primCStensor=.FALSE.
nbast1 = MAT%nrow
nbast2 = MAT%ncol
natom1 = AOT1%natoms
natom2 = AOT2%natoms
NULLIFY(TENSOR%LSAO)
NULLIFY(TENSOR%INDEX)
ALLOCATE(TENSOR%LSAO(natom1*natom2))
ALLOCATE(TENSOR%INDEX(natom1,natom2,1,1))
TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
TENSOR%natom1 = natom1
TENSOR%natom2 = natom2
TENSOR%natom3 = 1
TENSOR%natom4 = 1
TENSOR%nbast1 = nbast1
TENSOR%nbast2 = nbast2
TENSOR%nbast3 = 1
TENSOR%nbast4 = 1
TENSOR%nmat = 1
ATOMICINDEX=0
DO Iatom = 1,natom1
   nbastI = AOT1%ATOMICnORB(Iatom)
   nbatI = AOT1%ATOMICnBatch(IATOM)
   DO Jatom = 1,natom2
      nbastJ = AOT2%ATOMICnOrb(Iatom)
      nbatJ = AOT2%ATOMICnBatch(JATOM)
      !DETERMINE MAXELM FOR THIS ATOM-PAIR
      MAXELM = 0.d0
      DO Ibat = 1,nbatI
         DO Jbat = 1,nbatJ
            DO Iangmom = 1,AOT1%BATCH(Ibat)%nAngmom 
               DO Jangmom = 1,AOT2%BATCH(Jbat)%nAngmom 
                  nContA = AOT1%BATCH(Ibat)%nContracted(Iangmom) 
                  nOrbCompA = AOT1%BATCH(Ibat)%nOrbComp(Iangmom)
                  nContB = AOT2%BATCH(Jbat)%nContracted(Jangmom)
                  nOrbCompB = AOT2%BATCH(Jbat)%nOrbComp(Jangmom)
                  sA = AOT1%BATCH(Ibat)%startOrbital(Iangmom)
                  sB = AOT2%BATCH(Jbat)%startOrbital(Jangmom)                     
                  DO I = 1,nContA*nOrbCompA
                     DO J = 1,nContB*nOrbCompB
                        MAXELM = MAX(MAXELM,ABS(MAT%elms(sA+I-1+(sA+J-2)*nbast1)))
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      IF(MAXELM .GT. 1.0d-10)THEN
         ATOMICINDEX = ATOMICINDEX + 1 
         TENSOR%INDEX(IATOM,JATOM,1,1) = ATOMICINDEX  
         NULLIFY(TENSOR%LSAO(ATOMICINDEX)%BATCH)
         ALLOCATE(TENSOR%LSAO(ATOMICINDEX)%BATCH(nbatI,nbatJ,1,1))
         TENSOR%LSAO(ATOMICINDEX)%ALLOC = .TRUE.
         TENSOR%LSAO(ATOMICINDEX)%ATOM1 = Iatom
         TENSOR%LSAO(ATOMICINDEX)%ATOM2 = Jatom
         TENSOR%LSAO(ATOMICINDEX)%ATOM3 = 1
         TENSOR%LSAO(ATOMICINDEX)%ATOM4 = 1
         TENSOR%LSAO(ATOMICINDEX)%nAOBATCH1 = nbatI
         TENSOR%LSAO(ATOMICINDEX)%nAOBATCH2 = nbatJ
         TENSOR%LSAO(ATOMICINDEX)%nAOBATCH3 = 1
         TENSOR%LSAO(ATOMICINDEX)%nAOBATCH4 = 1
         DO Ibat = 1,nbatI
            DO Jbat = 1,nbatJ
               TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nAngmomA = AOT1%BATCH(Ibat)%nAngmom 
               TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nAngmomB = AOT2%BATCH(Jbat)%nAngmom 
               TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nAngmomC = 1
               TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nAngmomD = 1
               DO Iangmom = 1,AOT1%BATCH(Ibat)%nAngmom 
                  DO Jangmom = 1,AOT2%BATCH(Jbat)%nAngmom 
                     nContA = AOT1%BATCH(Ibat)%nContracted(Iangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nContA(Iangmom) = nContA
                     nOrbCompA = AOT1%BATCH(Ibat)%nOrbComp(Iangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nOrbCompA(Iangmom) = nOrbCompA
                     sA = AOT1%BATCH(Ibat)%startOrbital(Iangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%startOrbA(Iangmom) = sA
                     nContB = AOT2%BATCH(Jbat)%nContracted(Jangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nContB(Jangmom) = nContB
                     nOrbCompB = AOT2%BATCH(Jbat)%nOrbComp(Jangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nOrbCompB(Jangmom) = nOrbCompB
                     sB = AOT2%BATCH(Jbat)%startOrbital(Jangmom)                     
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%startOrbB(Jangmom) = sB
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nContC(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nOrbCompC(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%startOrbC(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nContD(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nOrbCompD(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%startOrbD(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nelmE = 1
                     nElms = nContA*nContB*nOrbCompA*nOrbCompB
                     NULLIFY(TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%elms)
                     ALLOCATE(TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%elms(nElms))
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nelms = nElms
                     DO J = 1,nContB*nOrbCompB
                        DO I = 1,nContA*nOrbCompA
                           TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%elms(I+(J-1)*nContA*nOrbCompA) &
     & = MAT%elms(sA+I-1+(sA+J-2)*nbast1)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDDO
TENSOR%nLSAO = ATOMICINDEX
Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE Build_lstensor_from_dense_mat_single

!> \brief build lstensor from an array of dense matrix type
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param AOT1 the 1. Atomic orbital
!> \param AOT2 the 2. Atomic orbital
!> \param mat the matrix type array
!> \param nmat the number of matrices
SUBROUTINE Build_lstensor_from_dense_mat_array(TENSOR,AOT1,AOT2,MAT,nmat)
  use Matrix_module
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM)       :: AOT1,AOT2
INTEGER            :: nmat
TYPE(MATRIX)       :: MAT(nmat)
!
REAL(REALK) :: MAXELM
INTEGER :: ATOMICINDEX,nElms
INTEGER :: nbast1,nbast2,natom1,natom2,Iatom,nbastI,nbatI
INTEGER :: nbastJ,nbatJ,Ibat,Jbat,iangmom,jangmom,jatom
INTEGER :: nContA,nOrbCompA,nContB,nOrbCompB,sA,sB,I,J,Imat
integer(kind=long) :: nmemsize

TENSOR%primCStensor=.FALSE.
nbast1 = MAT(1)%nrow
nbast2 = MAT(1)%ncol
natom1 = AOT1%natoms
natom2 = AOT2%natoms
nullify(TENSOR%LSAO)
nullify(TENSOR%INDEX)
ALLOCATE(TENSOR%LSAO(natom1*natom2))
ALLOCATE(TENSOR%INDEX(natom1,natom2,1,1))
TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
TENSOR%natom1 = natom1
TENSOR%natom2 = natom2
TENSOR%natom3 = 1
TENSOR%natom4 = 1
TENSOR%nbast1 = nbast1
TENSOR%nbast2 = nbast2
TENSOR%nbast3 = 1
TENSOR%nbast4 = 1
TENSOR%nmat = nmat
ATOMICINDEX=0
DO Iatom = 1,natom1
   nbastI = AOT1%ATOMICnOrb(Iatom)
   nbatI = AOT1%ATOMICnBatch(IATOM)
   DO Jatom = 1,natom2
      nbastJ = AOT2%ATOMICnOrb(Iatom)
      nbatJ = AOT2%ATOMICnBatch(JATOM)
      !DETERMINE MAXELM FOR THIS ATOM-PAIR
      MAXELM = 0.d0
      DO Ibat = 1,nbatI
         DO Jbat = 1,nbatJ
            DO Iangmom = 1,AOT1%BATCH(Ibat)%nAngmom 
               DO Jangmom = 1,AOT2%BATCH(Jbat)%nAngmom 
                  nContA = AOT1%BATCH(Ibat)%nContracted(Iangmom) 
                  nOrbCompA = AOT1%BATCH(Ibat)%nOrbComp(Iangmom)
                  nContB = AOT2%BATCH(Jbat)%nContracted(Jangmom)
                  nOrbCompB = AOT2%BATCH(Jbat)%nOrbComp(Jangmom)
                  sA = AOT1%BATCH(Ibat)%startOrbital(Iangmom)
                  sB = AOT2%BATCH(Jbat)%startOrbital(Jangmom)  
                  DO IMAT = 1,NMAT
                     DO I = 1,nContA*nOrbCompA
                        DO J = 1,nContB*nOrbCompB
                           MAXELM = MAX(MAXELM,ABS(MAT(IMAT)%elms(sA+I-1+(sA+J-2)*nbast1)))
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      IF(MAXELM .GT. 1.0d-10)THEN
         ATOMICINDEX = ATOMICINDEX + 1 
         TENSOR%INDEX(IATOM,JATOM,1,1) = ATOMICINDEX  
         nullify(TENSOR%LSAO(ATOMICINDEX)%BATCH)
         ALLOCATE(TENSOR%LSAO(ATOMICINDEX)%BATCH(nbatI,nbatJ,1,1))
         TENSOR%LSAO(ATOMICINDEX)%ALLOC = .TRUE.
         TENSOR%LSAO(ATOMICINDEX)%ATOM1 = Iatom
         TENSOR%LSAO(ATOMICINDEX)%ATOM2 = Jatom
         TENSOR%LSAO(ATOMICINDEX)%ATOM3 = 1
         TENSOR%LSAO(ATOMICINDEX)%ATOM4 = 1
         TENSOR%LSAO(ATOMICINDEX)%nAOBATCH1 = nbatI
         TENSOR%LSAO(ATOMICINDEX)%nAOBATCH2 = nbatJ
         TENSOR%LSAO(ATOMICINDEX)%nAOBATCH3 = 1
         TENSOR%LSAO(ATOMICINDEX)%nAOBATCH4 = 1
         DO Ibat = 1,nbatI
            DO Jbat = 1,nbatJ
               TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nAngmomA = AOT1%BATCH(Ibat)%nAngmom 
               TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nAngmomB = AOT2%BATCH(Jbat)%nAngmom 
               TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nAngmomC = 1
               TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nAngmomD = 1
               DO Iangmom = 1,AOT1%BATCH(Ibat)%nAngmom 
                  DO Jangmom = 1,AOT2%BATCH(Jbat)%nAngmom 
                     nContA = AOT1%BATCH(Ibat)%nContracted(Iangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nContA(Iangmom) = nContA
                     nOrbCompA = AOT1%BATCH(Ibat)%nOrbComp(Iangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nOrbCompA(Iangmom) = nOrbCompA
                     sA = AOT1%BATCH(Ibat)%startOrbital(Iangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%startOrbA(Iangmom) = sA
                     nContB = AOT2%BATCH(Jbat)%nContracted(Jangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nContB(Jangmom) = nContB
                     nOrbCompB = AOT2%BATCH(Jbat)%nOrbComp(Jangmom)
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nOrbCompB(Jangmom) = nOrbCompB
                     sB = AOT2%BATCH(Jbat)%startOrbital(Jangmom)                     
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%startOrbB(Jangmom) = sB
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nContC(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nOrbCompC(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%startOrbC(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nContD(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nOrbCompD(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%startOrbD(1) = 1
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nelmE = 1
                     nElms = nContA*nContB*nOrbCompA*nOrbCompB
                     TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%nelms = nElms
                     nullify(TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%elms)
                     ALLOCATE(TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%elms(nElms*nmat))
                     DO IMAT = 1,NMAT
                        DO J = 1,nContB*nOrbCompB
                           DO I = 1,nContA*nOrbCompA
                              TENSOR%LSAO(ATOMICINDEX)%BATCH(Ibat,Jbat,1,1)%elms(I+(J-1)*nContA*nOrbCompA+(IMAT-1)*nElms) &
     & = MAT(IMAT)%elms(sA+I-1+(sA+J-2)*nbast1)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDDO
TENSOR%nLSAO = ATOMICINDEX
Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE Build_lstensor_from_dense_mat_array

!> \brief restructure the primitive screening to fit with standard lstensor structure
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the output lstensor
!> \param pTENSOR the primitive lstensor
!> \param cAO1 the 1. Atomic orbital
!> \param cAO2 the 2. Atomic orbital
!> \param lupri the logical unit number of the output file
SUBROUTINE restructure_primGab_lstensor(TENSOR,pTENSOR,cAO1,cAO2,lupri)!,RefScreenMat,n1,n2,lupri)
implicit none
TYPE(LSTENSOR)       :: TENSOR,pTENSOR
TYPE(AOITEM),target  :: cAO1,cAO2
TYPE(AOITEM),pointer :: AOT1,AOT2
INTEGER            :: nbast1,nbast2,lupri!,n1,n2
logical   :: useAO1,useAO2
!
TYPE(AOITEM),target :: AOT
!REAL(REALK) :: MAXELM
INTEGER :: nElms
! BELONGINING TO AOT1 
INTEGER :: natom1,nbastI,nbatI,Ibat,Iangmom,Iatom,IORB
INTEGER :: nContA,nOrbCompA,sA,nAngA,AOT1batch,batAOT1,nmat
! BELONGINING TO AOT2 
INTEGER :: natom2,nbastJ,nbatJ,Jbat,Jangmom,Jatom,JORB
INTEGER :: nContB,nOrbCompB,sB,nAngB,AOT2batch,batAOT2
! COMMON
INTEGER :: DIM,IELM,IMAT2,I,npbatI,npbatJ,npbatK,npbatL,mark_pIbat,mark_pJbat,mark_pKbat
INTEGER :: mark_pLbat,pIbat,pJbat,pKbat,pLbat,ip,jp,kp,lp,dima,dimb,dimc,dimd,na,nb,nc,nd
INTEGER :: lock_startD,lock_startC,lock_startB,lock_startA,pIELM,startD,startA,startB,startC
INTEGER :: IA,IB,IC,ID,nP1,nP2,pIbat_start,pJbat_start,oldIbat,oldJbat,cIbat,cJbat
LOGICAL :: ERROR
real(realk) :: maxelm,DIFF
!real(realk) :: RefscreenMat(n1,n2)
integer  :: atomA,atomB,index,iOrbCompA,iOrbCompB,J,Abatch,Bbatch,iP1S,iP2S,ip1,ip2,pbatAOT1,pbatAOT2
real(realk),pointer :: temp(:,:)
integer,pointer  :: pIbat2Ibat(:),pJbat2Jbat(:),pIELMs2(:)
integer(kind=long) :: nmemsize
TENSOR%primCStensor=.TRUE.
IF(pTENSOR%gradienttensor)CALL lsQUIT('wrong input to restructure_primGab_lstensor  gradinet tensor',-1)
IF(pTENSOR%natom3.NE.1)CALL lsQUIT('wrong input to restructure_primGab_lstensor 1',-1)
IF(pTENSOR%natom4.NE.1)CALL lsQUIT('wrong input to restructure_primGab_lstensor 2',-1)
IF(pTENSOR%nbast3.NE.1)CALL lsQUIT('wrong input to restructure_primGab_lstensor 3',-1)
IF(pTENSOR%nbast4.NE.1)CALL lsQUIT('wrong input to restructure_primGab_lstensor 4',-1)
IF(pTENSOR%nmat.NE.1)CALL lsQUIT('wrong input to restructure_primGab_lstensor 5',-1)
TENSOR%gradienttensor = .FALSE.
call SET_EMPTY_AO(AOT)

IF(cAO1%empty)THEN
   AOT1 => AOT
ELSE
   AOT1 => cAO1
ENDIF

IF(cAO2%empty)THEN
   AOT2 => AOT
ELSE
   AOT2 => cAO2
ENDIF

nbast1 = AOT1%nbast
nbast2 = AOT2%nbast
natom1 = AOT1%natoms
natom2 = AOT2%natoms
NULLIFY(TENSOR%LSAO)
NULLIFY(TENSOR%INDEX)
ALLOCATE(TENSOR%LSAO(natom1*natom2))
ALLOCATE(TENSOR%INDEX(natom1,natom2,1,1))
TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
TENSOR%natom1 = natom1
TENSOR%natom2 = natom2
TENSOR%natom3 = 1
TENSOR%natom4 = 1
TENSOR%nbast1 = nbast1
TENSOR%nbast2 = nbast2
TENSOR%nbast3 = 1
TENSOR%nbast4 = 1
TENSOR%nmat = pTENSOR%nmat 
nmat = pTENSOR%nmat 
AOT1batch=0
I = 0
DO Iatom = 1,natom1
 nbatI = AOT1%ATOMICnBatch(IATOM)
 AOT2batch=0
 DO Jatom = 1,natom2
  nbatJ = AOT2%ATOMICnBatch(JATOM)
  I=I+1
  TENSOR%INDEX(IATOM,JATOM,1,1) = I
  NULLIFY(TENSOR%LSAO(I)%BATCH)
  ALLOCATE(TENSOR%LSAO(I)%BATCH(nbatI,nbatJ,1,1))
  TENSOR%LSAO(I)%ALLOC = .TRUE.
  TENSOR%LSAO(I)%ATOM1 = Iatom
  TENSOR%LSAO(I)%ATOM2 = Jatom
  TENSOR%LSAO(I)%ATOM3 = 1
  TENSOR%LSAO(I)%ATOM4 = 1
  TENSOR%LSAO(I)%nAOBATCH1 = nbatI
  TENSOR%LSAO(I)%nAOBATCH2 = nbatJ
  TENSOR%LSAO(I)%nAOBATCH3 = 1
  TENSOR%LSAO(I)%nAOBATCH4 = 1

  INDEX = pTENSOR%INDEX(Iatom,Jatom,1,1)
  IF(INDEX .NE. I)CALL lsQUIT('Mismatch between The Contract AOs and primLstensor screening',-1)
  npbatI = pTENSOR%LSAO(INDEX)%nAOBATCH1
  npbatJ = pTENSOR%LSAO(INDEX)%nAOBATCH2  

  batAOT1 = AOT1batch
  DO Ibat = 1,nbatI
   batAOT1 = batAOT1+1 
   nP1 = AOT1%BATCH(batAOT1)%nPrimitives
   batAOT2 = AOT2batch
   DO Jbat = 1,nbatJ
    batAOT2 = batAOT2+1 
    nP2 = AOT2%BATCH(batAOT2)%nPrimitives
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomA = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomB = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomC = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomD = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContA(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompA(1) = nP1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbA(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContB(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompB(1) = nP2
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbB(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContC(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompC(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbC(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContD(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompD(1) = 1 
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbD(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nelmE = 1
    DIM = nP1*nP2
    NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms)
    ALLOCATE(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(DIM))
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nelms = DIM
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms = 0.d0
   ENDDO 
  ENDDO

!  print*,'INDEX',INDEX,'I=',I,'Iatom',Iatom,'Jatom',Jatom,'npbatI,npbatJ',npbatI,npbatJ
  write(lupri,*)'IATOM,JATOM',IATOM,JATOM

  batAOT2 = AOT2batch
  pbatAOT2 = 0
  DO Jbat = 1,nbatJ
   batAOT2 = batAOT2+1 
   nP2 = AOT2%BATCH(batAOT2)%nPrimitives
   do iP2=1,nP2
    pbatAOT2 = pbatAOT2+1
    call mem_alloc(pIELMs2,npbatI)
    pIELMs2 = 0
    DO Jangmom = 1,AOT2%BATCH(batAOT2)%nAngmom
     nOrbCompB = AOT2%BATCH(batAOT2)%nOrbComp(Jangmom)
     DO JORB = 1,nOrbCompB
      batAOT1 = AOT1batch
      pbatAOT1 = 0
      DO Ibat = 1,nbatI
       batAOT1 = batAOT1+1 
       nP1 = AOT1%BATCH(batAOT1)%nPrimitives
       do iP1=1,nP1
        pbatAOT1 = pbatAOT1+1
        IELM = iP1 + (iP2-1)*nP1
        DO Iangmom = 1,pTENSOR%LSAO(INDEX)%BATCH(pbatAOT1,pbatAOT2,1,1)%nAngmomA
         nOrbCompA = pTENSOR%LSAO(INDEX)%BATCH(pbatAOT1,pbatAOT2,1,1)%nOrbCompA(Iangmom)
         DO IORB = 1,nOrbCompA
!            write(lupri,*)'pTENSOR(',pbatAOT1,',',pbatAOT2,') -> cTENSOR(',Ibat,',',Jbat,')'
          pIELMs2(pbatAOT1) = pIELMs2(pbatAOT1)+1 
          TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(IELM) = MAX(&
               & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(IELM),&
               & pTENSOR%LSAO(INDEX)%BATCH(pbatAOT1,pbatAOT2,1,1)%elms(pIELMs2(pbatAOT1)))
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      IF(pbatAOT1 .NE.npbatI )call lsQUIT('dimension mismatch Restructure_primGab_lstensor 1',-1)
     ENDDO
    ENDDO
    call mem_dealloc(pIELMs2)
   enddo
  ENDDO
  IF(pbatAOT2 .NE.npbatJ )call lsQUIT('dimension mismatch Restructure_primGab_lstensor 2',-1)
!!$  iP2 = 0
!!$  oldJbat = 0
!!$  batAOT2 = AOT2batch
!!$  DO pJbat = 1,npbatJ
!!$   iP2 = iP2 + 1
!!$   Jbat = pJbat2Jbat(pJbat)
!!$   IF(Jbat .NE. oldJbat)THEN
!!$      iP2 = 1
!!$      batAOT2 = batAOT2+1 
!!$   ENDIF
!!$   oldJbat = Jbat
!!$   call mem_alloc(pIELMs2,npbatI)
!!$   pIELMs2 = 0
!!$   DO Jangmom = 1,AOT2%BATCH(batAOT2)%nAngmom
!!$    nOrbCompB = AOT2%BATCH(batAOT2)%nOrbComp(Jangmom)
!!$    DO JORB = 1,nOrbCompB
!!$     iP1 = 0
!!$     oldIbat = 0
!!$     batAOT1 = AOT1batch
!!$     DO pIbat = 1,npbatI 
!!$      iP1 = iP1 + 1
!!$      Ibat = pIbat2Ibat(pIbat)
!!$      IF(Ibat .NE. oldIbat)THEN
!!$         iP1 = 1
!!$         batAOT1 = batAOT1+1 
!!$         nP1 = AOT1%BATCH(batAOT1)%nPrimitives
!!$         IF(AOT1%empty)nP1=1
!!$      ENDIF
!!$      oldIbat = Ibat
!!$      IELM = iP1 + (iP2-1)*nP1
!!$      DO Iangmom = 1,pTENSOR%LSAO(INDEX)%BATCH(pIbat,pJbat,1,1)%nAngmomA
!!$       nOrbCompA = pTENSOR%LSAO(INDEX)%BATCH(pIbat,pJbat,1,1)%nOrbCompA(Iangmom)
!!$       DO IORB = 1,nOrbCompA
!!$        pIELMs2(pIbat) = pIELMs2(pIbat)+1 
!!$        TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(IELM) = MAX(&
!!$            & TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(IELM),&
!!$            & pTENSOR%LSAO(INDEX)%BATCH(pIbat,pJbat,1,1)%elms(pIELMs2(pIbat)))
!!$!        WRITE(lupri,'(A,I2,A,I2,A,2X,E11.5,2X,A,2X,E11.5,2X,A,I2,A,I2)')'temp(',iP1,',',iP2,')= ',pTENSOR%LSAO(INDEX)%BATCH(pIbat,pJbat,1,1)%elms(pIELMs2(pIbat)),&
!!$!             & '  max  ',TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(IELM),'  globalIbat=',batAOT1,'  globalJbat=',batAOT2
!!$       ENDDO
!!$      ENDDO
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$   call mem_dealloc(pIELMs2)
!!$  ENDDO

!  call mem_dealloc(pIbat2Ibat)
!  call mem_dealloc(pJbat2Jbat)
  AOT2batch = AOT2batch + nbatJ
 ENDDO
 AOT1batch = AOT1batch + nbatI
ENDDO
TENSOR%nLSAO = I
call FREE_EMPTY_AO(AOT)

!!$WRITE(LUPRI,*)'REF screenmat FULL',AOT1%nbatches,AOT2%nbatches
!!$call output(refScreenMat,1,n1,1,n2,n1,n2,1,lupri)      
!!$DO Ibat = 1,AOT1%nbatches
!!$ write(lupri,*),'Ibat',Ibat,'of',AOT1%nbatches
!!$ atomA =  AOT1%BATCH(Ibat)%atom
!!$ nP1 = AOT1%BATCH(Ibat)%nPrimitives
!!$ DO Jbat = 1,AOT2%nbatches
!!$  write(lupri,*),'Jbat',Jbat,'of',AOT2%nbatches
!!$  atomB =  AOT2%BATCH(Jbat)%atom
!!$  nP2 = AOT2%BATCH(Jbat)%nPrimitives
!!$  INDEX = TENSOR%INDEX(atomA,atomB,1,1)
!!$  maxelm = 0      
!!$  call mem_alloc(temp,nP1,nP2)
!!$  temp = 0.d0
!!$  write(lupri,*),'nangmomA',AOT1%BATCH(Ibat)%nAngmom
!!$  DO Iangmom = 1,AOT1%BATCH(Ibat)%nAngmom
!!$   write(lupri,*),'IangmomA',iangmom
!!$   sA = AOT1%BATCH(Ibat)%startPrimOrbital(Iangmom)         
!!$   nOrbCompA = AOT1%BATCH(Ibat)%nOrbComp(Iangmom)
!!$   nA = sA-1+nOrbCompA*nP1
!!$   write(lupri,*),'nangmomB',AOT2%BATCH(Jbat)%nAngmom
!!$   DO Jangmom = 1,AOT2%BATCH(Jbat)%nAngmom
!!$    write(lupri,*),'JangmomB',Jangmom
!!$    sB = AOT2%BATCH(Jbat)%startPrimOrbital(Jangmom)                     
!!$    nOrbCompB = AOT2%BATCH(Jbat)%nOrbComp(Jangmom)
!!$    nB = sB-1+nOrbCompB*nP2
!!$    WRITE(LUPRI,*)'REF screenmat Iangmom,Jangmom,sA,nA,sB,nB,n1,n2',Iangmom,Jangmom,sA,nA,sB,nB,n1,n2
!!$    call output(refScreenMat,sA,nA,sB,nB,n1,n2,1,lupri)      
!!$    do iP = 1,nP1
!!$     do iOrbCompA = 1,nOrbCompA
!!$      I = sA - 1 + iOrbCompA + (iP-1)*nOrbCompA
!!$      if(iP.EQ.2) write(lupri,*)'for iP=2 I=',I
!!$      do jP = 1,nP2
!!$       do iOrbCompB = 1,nOrbCompB
!!$        J = sB - 1 + iOrbCompB + (jP-1)*nOrbCompB
!!$        if((iP.EQ.2).AND.(jP.EQ.3)) write(lupri,*)'for jP=3 J=',J,'elm',refScreenMat(I,J)
!!$        temp(iP,jP) = MAX(temp(iP,jP),refScreenMat(I,J))
!!$       enddo
!!$      enddo
!!$     enddo
!!$    enddo            
!!$   ENDDO
!!$  ENDDO
!!$  WRITE(LUPRI,'(A,I4,I4)')'REF screenmat'
!!$  call output(temp,1,nP1,1,nP2,nP1,nP2,1,lupri)      
!!$
!!$  batAOT1 = 0
!!$  DO Iatom = 1,natom1
!!$     nbatI = AOT1%ATOMICnBatch(IATOM)
!!$     DO cIbat = 1,nbatI
!!$        batAOT1 = batAOT1+1 
!!$        IF(batAOT1.EQ.Ibat)Abatch = cIbat
!!$     ENDDO
!!$  ENDDO
!!$  batAOT2 = 0
!!$  DO Jatom = 1,natom2
!!$     nbatJ = AOT2%ATOMICnBatch(JATOM)
!!$     DO cJbat = 1,nbatJ
!!$        batAOT2 = batAOT2+1 
!!$        IF(batAOT2.EQ.Jbat)Bbatch = cJbat
!!$     ENDDO
!!$  ENDDO
!!$  
!!$  WRITE(LUPRI,'(A,I4,I4)')'MY VERSION of prim'
!!$  call output(TENSOR%LSAO(index)%BATCH(Abatch,Bbatch,1,1)%elms,1,nP1,1,nP2,nP1,nP2,1,lupri)      
!!$  DIFF = 0.d0
!!$  do iP = 1,nP1
!!$     do jP = 1,nP2
!!$        IELM = iP + (jP-1)*nP1
!!$        DIFF = DIFF + ABS(temp(iP,jP) - TENSOR%LSAO(index)%BATCH(Abatch,Bbatch,1,1)%elms(ielm))
!!$     enddo
!!$  enddo
!!$  print*,'DIFF',DIFF
!!$  write(lupri,*)'DIFF',DIFF
!!$  !      IF(DIFF.GT.1.d-16)CALL lsQUIT('somethikfjnhbbbbbgfgfg')
!!$  call mem_dealloc(temp)
!!$ ENDDO
!!$ENDDO
Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE Restructure_primGab_lstensor
!-----------------------------------------------------
!!$!> \brief init the lstensor for the primitive gab matrix
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param Output contains output specs for integral storage
!!$!> \param TENSOR the output lstensor
!!$!> \param pAO1 the 1. primitive Atomic orbital
!!$!> \param pAO2 the 2. primitive Atomic orbital
!!$!> \param cAO1 the 1. contracted Atomic orbital
!!$!> \param cAO2 the 2. contracted Atomic orbital
!!$!> \param lupri the logical unit number of the output file
!!$SUBROUTINE init_primGablstensor(Output,TENSOR,psAO1,psAO2,csAO1,csAO2,lupri)
!!$  use integraloutput_type
!!$implicit none
!!$Type(IntegralOutput) :: Output
!!$TYPE(LSTENSOR)       :: TENSOR
!!$TYPE(AOITEM),target  :: cAO1,cAO2
!!$TYPE(AOITEM),pointer :: AOT1,AOT2
!!$INTEGER            :: nbast1,nbast2,lupri!,n1,n2
!!$logical   :: useAO1,useAO2
!!$!
!!$TYPE(AOITEM),target :: AOT
!!$!REAL(REALK) :: MAXELM
!!$INTEGER :: nElms
!!$! BELONGINING TO AOT1 
!!$INTEGER :: natom1,nbastI,nbatI,Ibat,Iangmom,Iatom,IORB
!!$INTEGER :: nContA,nOrbCompA,sA,nAngA,AOT1batch,batAOT1,nmat
!!$! BELONGINING TO AOT2 
!!$INTEGER :: natom2,nbastJ,nbatJ,Jbat,Jangmom,Jatom,JORB
!!$INTEGER :: nContB,nOrbCompB,sB,nAngB,AOT2batch,batAOT2
!!$! COMMON
!!$INTEGER :: DIM,IELM,IMAT2,I,npbatI,npbatJ,npbatK,npbatL,mark_pIbat,mark_pJbat,mark_pKbat
!!$INTEGER :: mark_pLbat,pIbat,pJbat,pKbat,pLbat,ip,jp,kp,lp,dima,dimb,dimc,dimd,na,nb,nc,nd
!!$INTEGER :: lock_startD,lock_startC,lock_startB,lock_startA,pIELM,startD,startA,startB,startC
!!$INTEGER :: IA,IB,IC,ID,nP1,nP2,pIbat_start,pJbat_start,oldIbat,oldJbat,cIbat,cJbat
!!$LOGICAL :: ERROR
!!$real(realk) :: maxelm,DIFF
!!$!real(realk) :: RefscreenMat(n1,n2)
!!$integer  :: atomA,atomB,index,iOrbCompA,iOrbCompB,J,Abatch,Bbatch,iP1S,iP2S,ip1,ip2,pbatAOT1,pbatAOT2
!!$real(realk),pointer :: temp(:,:)
!!$integer,pointer  :: pIbat2Ibat(:),pJbat2Jbat(:),pIELMs2(:)
!!$TENSOR%primCStensor=.TRUE.
!!$TENSOR%gradienttensor = .FALSE.
!!$call SET_EMPTY_AO(AOT)
!!$
!!$IF(cAO1%empty)THEN
!!$   cAOT1 => AOT
!!$   pAOT1 => AOT
!!$ELSE
!!$   cAOT1 => csAO1
!!$   pAOT1 => psAO1
!!$ENDIF
!!$
!!$IF(cAO2%empty)THEN
!!$   cAOT2 => AOT
!!$   pAOT2 => AOT
!!$ELSE
!!$   cAOT2 => csAO2
!!$   pAOT2 => psAO2
!!$ENDIF
!!$
!!$nbast1 = cAOT1%nbast
!!$nbast2 = cAOT2%nbast
!!$natom1 = cAOT1%natoms
!!$natom2 = cAOT2%natoms
!!$NULLIFY(TENSOR%LSAO)
!!$NULLIFY(TENSOR%INDEX)
!!$ALLOCATE(TENSOR%LSAO(natom1*natom2))
!!$ALLOCATE(TENSOR%INDEX(natom1,natom2,1,1))
!!$TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
!!$TENSOR%natom1 = natom1
!!$TENSOR%natom2 = natom2
!!$TENSOR%natom3 = 1
!!$TENSOR%natom4 = 1
!!$TENSOR%nbast1 = nbast1
!!$TENSOR%nbast2 = nbast2
!!$TENSOR%nbast3 = 1
!!$TENSOR%nbast4 = 1
!!$TENSOR%nmat = nmat 
!!$AOT1batch=0
!!$I = 0
!!$DO Iatom = 1,natom1
!!$ nbatI = cAOT1%ATOMICnBatch(IATOM)
!!$ AOT2batch=0
!!$ DO Jatom = 1,natom2
!!$  nbatJ = cAOT2%ATOMICnBatch(JATOM)
!!$  I=I+1
!!$  TENSOR%INDEX(IATOM,JATOM,1,1) = I
!!$  NULLIFY(TENSOR%LSAO(I)%BATCH)
!!$  ALLOCATE(TENSOR%LSAO(I)%BATCH(nbatI,nbatJ,1,1))
!!$  TENSOR%LSAO(I)%ALLOC = .TRUE.
!!$  TENSOR%LSAO(I)%ATOM1 = Iatom
!!$  TENSOR%LSAO(I)%ATOM2 = Jatom
!!$  TENSOR%LSAO(I)%ATOM3 = 1
!!$  TENSOR%LSAO(I)%ATOM4 = 1
!!$  TENSOR%LSAO(I)%nAOBATCH1 = nbatI
!!$  TENSOR%LSAO(I)%nAOBATCH2 = nbatJ
!!$  TENSOR%LSAO(I)%nAOBATCH3 = 1
!!$  TENSOR%LSAO(I)%nAOBATCH4 = 1
!!$
!!$  npbatI = pAOT1%ATOMICnBatch(IATOM)
!!$  npbatJ = pAOT2%ATOMICnBatch(JATOM)
!!$
!!$  batAOT1 = AOT1batch
!!$  DO Ibat = 1,nbatI
!!$   batAOT1 = batAOT1+1 
!!$   nP1 = cAOT1%BATCH(batAOT1)%nPrimitives
!!$   batAOT2 = AOT2batch
!!$   DO Jbat = 1,nbatJ
!!$    batAOT2 = batAOT2+1 
!!$    nP2 = cAOT2%BATCH(batAOT2)%nPrimitives
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomA = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomB = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomC = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomD = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContA(1) = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompA(1) = nP1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbA(1) = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContB(1) = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompB(1) = nP2
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbB(1) = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContC(1) = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompC(1) = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbC(1) = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContD(1) = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompD(1) = 1 
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbD(1) = 1
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nelmE = 1
!!$    DIM = nP1*nP2
!!$    NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms)
!!$    ALLOCATE(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(DIM))
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nelms = DIM
!!$    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms = 0.d0
!!$   ENDDO 
!!$  ENDDO
!!$  AOT2batch = AOT2batch + nbatJ
!!$ ENDDO
!!$ AOT1batch = AOT1batch + nbatI
!!$ENDDO
!!$TENSOR%nLSAO = I
!!$
!!$CALL build_primassistarrays(OUTPUT%CnPrimitivesA,Output%startBatchA,Output%ContractedBatchA,natom1,cAOT1)
!!$
!!$CALL build_primassistarrays(OUTPUT%CnPrimitivesB,Output%startBatchB,Output%ContractedBatchB,natom2,cAOT2)
!!$
!!$DO Iatom = 1,natom
!!$   WRITE(lupri,*)'IATOM:',IATOM
!!$   pnbatI = pAOT%ATOMICnBatch(IATOM)
!!$   DO pIbat = 1,pnbatI
!!$      WRITE(lupri,*)'THE ContractedBatchA, StartBatchA     nr pbatches',pnbatI
!!$      write(lupri,*)Output%ContractedBatchA(IATOM)%BATCH(pIbat),&
!!$           &Output%startBatchA(IATOM)%BATCH(pIbat)
!!$   enddo
!!$   WRITE(lupri,*)'THE nPrimitivesA                      nr batches',cAOT%ATOMICnBatch(IATOM)
!!$   do Ibat=1,cAOT%ATOMICnBatch(IATOM)
!!$      write(lupri,*)Output%CnPrimitivesA(Ibat)
!!$   enddo
!!$enddo
!!$call FREE_EMPTY_AO(AOT)
!!$
!!$END SUBROUTINE Init_primGablstensor
!!$
!!$!> \brief build som arrays needed to place primitive screening integrals in proper place in lstensor
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param CnPrimitivesA number of primitives in the contracted basis
!!$!> \param startBatchA starting index in contracted basis
!!$!> \param ContractedBatchA contracted batch in contracted basis
!!$!> \param natom the number of atoms in the atomic orbitals
!!$!> \param cAOT a contracted Atomic orbital
!!$SUBROUTINE build_primassistarrays(CnPrimitivesA,startBatchA,ContractedBatchA,&
!!$     &natom,cAOT)
!!$implicit none
!!$type(intbatch) :: CnPrimitivesA(:),startBatchA(:),ContractedBatchA(:)
!!$TYPE(AOITEM)   :: cAOT
!!$integer :: natom
!!$!
!!$integer :: AOTbatch, nbatI,pIbat,iP1,nP1,Ibat
!!$
!!$AOTbatch=0
!!$ALLOCATE(CnPrimitivesA(natom1))
!!$ALLOCATE(startBatchA(natom1))
!!$ALLOCATE(ContractedBatchA(natom1))
!!$DO Iatom = 1,natom
!!$ nbatI = cAOT%ATOMICnBatch(IATOM)
!!$ ALLOCATE(CnPrimitivesA(IATOM)%BATCH(nbatI))
!!$ pIbat = 0
!!$ DO Ibat = 1,nbatI
!!$    nP1 = cAOT%BATCH(AOTbatch+Ibat)%nPrimitives
!!$    do iP1=1,nP1
!!$       pIbat = pIbat+1
!!$    enddo
!!$ ENDDO
!!$ ALLOCATE(startBatchA(IATOM)%BATCH(pIbat))
!!$ ALLOCATE(ContractedBatchA(IATOM)%BATCH(pIbat))
!!$ pIbat = 0
!!$ DO Ibat = 1,nbatI
!!$    nP1 = cAOT%BATCH(AOTbatch+Ibat)%nPrimitives
!!$    CnPrimitivesA(IATOM)%BATCH(Ibat)=nP1
!!$    do iP1=1,nP1
!!$       pIbat = pIbat+1
!!$       startBatchA(IATOM)%BATCH(pIbat) = pIbat-1
!!$       ContractedBatchA(IATOM)%BATCH(pIbat) = Jbat
!!$    enddo
!!$ ENDDO
!!$ AOTbatch = AOTbatch + nbatI
!!$ENDDO
!!$END SUBROUTINE BUILD_PRIMASSISTARRAYS

!> \brief determine maximum element of the lstensor structure 
!> \author T. Kjaergaard
!> \date 2010
!> \param maxelm the maximum element of the lstensor
!> \param TENSOR the output lstensor
SUBROUTINE determine_maxelm(maxelm,TENSOR)
implicit none
TYPE(LSTENSOR)     :: TENSOR
real(realk) :: maxelm

INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nContA,nOrbCompA
INTEGER    :: sA,nContB,nOrbCompB,sB,nContC,nOrbCompC,sC,J,iatom,nmat,ip1,ip2,dim
INTEGER    :: nContD,nOrbCompD,sD,IORB,JORB,KORB,LORB,IELM,IMAT,i1,np1,np2,istart,iend

maxelm = 0.d0
IF(TENSOR%gradienttensor)THEN
 nmat = TENSOR%nmat
 DO Iatom = 1,TENSOR%natom1
  DO IMAT = 1,nmat
   do J = 1,3
    maxelm = MAX(maxelm,TENSOR%LSAO(Iatom)%BATCH(1,1,1,1)%elms(J+(IMAT-1)*3))
   enddo
  ENDDO
 ENDDO
ELSEIF(TENSOR%primCStensor)THEN
 DO I = 1,TENSOR%nLSAO
  IF(TENSOR%LSAO(I)%ALLOC)THEN
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     nP1 = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompA(1)
     nP2 = TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompB(1)
     DO IMAT = 1,TENSOR%nmat
      ielm = 0
      do iP1 = 1,nP1
       do iP2 = 1,nP2
        ielm = ielm+1
        maxelm = MAX(maxelm,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(ielm))
       enddo
      enddo
     ENDDO
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ELSE
 DO I = 1,TENSOR%nLSAO
  IF(TENSOR%LSAO(I)%ALLOC)THEN
   DO Ibat = 1,TENSOR%LSAO(I)%nAOBATCH1
    DO Jbat = 1,TENSOR%LSAO(I)%nAOBATCH2
     DO Kbat = 1,TENSOR%LSAO(I)%nAOBATCH3
      DO Lbat = 1,TENSOR%LSAO(I)%nAOBATCH4
       DIM = TENSOR%nmat*TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%nelms
       do i1=1,DIM
        maxelm = MAX(maxelm,TENSOR%LSAO(I)%BATCH(Ibat,Jbat,Kbat,Lbat)%elms(i1))
       enddo
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ENDIF

END SUBROUTINE determine_maxelm

!> \brief determine if the batch should be screened away  
!> \author T. Kjaergaard
!> \date 2010
!> \param nbatI number of batch for batch1
!> \param nbatJ number of batch for batch2
!> \param nbatK number of batch for batch3
!> \param nbatL number of batch for batch4
!> \param AOT1batch AO batch index 1
!> \param AOT2batch AO batch index 2
!> \param AOT3batch AO batch index 3
!> \param AOT4batch AO batch index 4
!> \param useAO1 use the AObatch 1 in contrast to use an empty 
!> \param useAO2 use the AObatch 2 in contrast to use an empty 
!> \param useAO3 use the AObatch 3 in contrast to use an empty 
!> \param useAO4 use the AObatch 4 in contrast to use an empty 
!> \param BATCHA3 the 1. AO batch 
!> \param BATCHB3 the 2. AO batch 
!> \param BATCHC3 the 3. AO batch 
!> \param BATCHD3 the 4. AO batch 
!> \param OVERALLscreen should the batch be screened
subroutine determineODscreening(nbatI,nbatJ,nbatK,nbatL,&
     & AOT1batch,AOT2batch,AOT3batch,AOT4batch,useAO1,useAO2,&
     & useAO3,useAO4,BATCHA3,BATCHB3,BATCHC3,BATCHD3,OVERALLscreen)
  use OD_Type, only: getODscreening
implicit none
integer,intent(in) :: nbatI,nbatJ,nbatK,nbatL,AOT1batch,AOT2batch,AOT3batch,AOT4batch
TYPE(AOBATCH),intent(in)  :: BATCHA3(:)
TYPE(AOBATCH),intent(in)  :: BATCHB3(:)
TYPE(AOBATCH),intent(in)  :: BATCHC3(:)
TYPE(AOBATCH),intent(in)  :: BATCHD3(:)
logical,intent(in) :: useAO1,useAO2,useAO3,useAO4
logical,intent(out) :: overallscreen
!
integer :: batAOT1,batAOT2,batAOT3,batAOT4
integer :: Ibat,Jbat,Kbat,Lbat
logical :: ODscreenLHS,ODscreenRHS

OVERALLscreen = .TRUE.

IF(useAO1 .OR. useAO2)THEN
 batAOT1 = AOT1batch
 DO Ibat = 1,nbatI
  batAOT1 = batAOT1+1 
  batAOT2 = AOT2batch
  DO Jbat = 1,nbatJ
   batAOT2 = batAOT2+1 
   ODscreenLHS = .FALSE.
   call getODscreening(BATCHA3(batAOT1),BATCHB3(batAOT2),ODscreenLHS) 
   IF(.NOT.ODscreenLHS)THEN
 !   print*,'ODscreenLHS',ODscreenLHS,'is false so we should calc'
    !if just one element is not screen - we do not screen any
    OVERALLscreen = .FALSE.
    EXIT
   ENDIF
  ENDDO
 ENDDO
ENDIF

IF(useAO3 .OR. useAO4)THEN
 IF(OVERALLscreen)THEN
  batAOT3 = AOT3batch
  DO Kbat = 1,nbatK
   batAOT3 = batAOT3 + 1 
   batAOT4 = AOT4batch
   DO Lbat = 1,nbatL
    batAOT4 = batAOT4 + 1  
    ODscreenRHS = .FALSE.
    call getODscreening(BATCHC3(batAOT3),BATCHD3(batAOT4),ODscreenRHS) 
    IF(.NOT.ODscreenRHS)THEN
     !if just one element is not screen - we do not screen any
     OVERALLscreen = .FALSE.
     EXIT
    ENDIF
   ENDDO
  ENDDO
 ENDIF
ENDIF

end subroutine determineODscreening

END MODULE LSTENSOR_OPERATIONSMOD
