module lsmpi_mod
#ifdef VAR_LSMPI
  include 'mpif.h'
#endif

  INTERFACE ls_mpibcast
     MODULE PROCEDURE ls_mpibcast_integer, ls_mpibcast_integerV,&
          &           ls_mpibcast_realk, ls_mpibcast_realkV,&
          &           ls_mpibcast_realkM, ls_mpibcast_realkT,&
          &           ls_mpibcast_logical, ls_mpibcast_logicalV,&
          &           ls_mpibcast_logicalM,&
          &           ls_mpibcast_charac, ls_mpibcast_characV
  END INTERFACE

contains
    subroutine ls_mpibcast_integer(buffer,master)
      implicit none
      integer :: buffer
      integer :: master,ierr,count,datatype

#ifdef VAR_LSMPI
      DATATYPE = MPI_INTEGER
      COUNT = 1
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_integer

    subroutine ls_mpibcast_integerV(buffer,nbuf,master)
      implicit none
      integer :: nbuf,master,ierr,count,datatype
      integer :: buffer(nbuf)

#ifdef VAR_LSMPI
      DATATYPE = MPI_INTEGER
      COUNT = nbuf
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_integerV

    subroutine ls_mpibcast_realk(buffer,master)
      use precision
      implicit none
      real(realk) :: buffer
      integer :: master,ierr,count,datatype

#ifdef VAR_LSMPI
      DATATYPE = MPI_DOUBLE_PRECISION
      COUNT = 1
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_realk

    subroutine ls_mpibcast_realkV(buffer,nbuf,master)
      use precision
      implicit none
      integer :: nbuf,master,ierr,count,datatype
      real(realk) :: buffer(nbuf)

#ifdef VAR_LSMPI
      DATATYPE = MPI_DOUBLE_PRECISION
      COUNT = nbuf
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_realkV

    subroutine ls_mpibcast_realkM(buffer,nbuf1,nbuf2,master)
      use precision
      implicit none
      integer :: nbuf1,nbuf2,master,ierr,count,datatype
      real(realk) :: buffer(nbuf1,nbuf2)

#ifdef VAR_LSMPI
      DATATYPE = MPI_DOUBLE_PRECISION
      COUNT = nbuf1*nbuf2
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_realkM

    subroutine ls_mpibcast_realkT(buffer,nbuf1,nbuf2,nbuf3,master)
      use precision
      implicit none
      integer :: master,nbuf1,nbuf2,nbuf3,ierr,count,datatype
      real(realk) :: buffer(nbuf1,nbuf2,nbuf3)

#ifdef VAR_LSMPI
      DATATYPE = MPI_DOUBLE_PRECISION
      COUNT = nbuf1*nbuf2*nbuf3
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_realkT

    subroutine ls_mpibcast_logical(buffer,master)
      implicit none
      logical :: buffer
      integer :: master,ierr,count,datatype

#ifdef VAR_LSMPI
      DATATYPE = MPI_LOGICAL
      COUNT = 1
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_logical

    subroutine ls_mpibcast_logicalV(buffer,nbuf,master)
      implicit none
      integer :: nbuf,master,ierr,count,datatype
      logical :: buffer(nbuf)

#ifdef VAR_LSMPI
      DATATYPE = MPI_LOGICAL
      COUNT = nbuf
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_logicalV

    subroutine ls_mpibcast_logicalM(buffer,nbuf1,nbuf2,master)
      implicit none
      integer :: nbuf1,nbuf2,master,ierr,count,datatype
      logical :: buffer(nbuf1,nbuf2)

#ifdef VAR_LSMPI
      DATATYPE = MPI_LOGICAL
      COUNT = nbuf1*nbuf2
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_logicalM

    subroutine ls_mpibcast_charac(buffer,master)
      implicit none
      character :: buffer
      integer :: master,ierr,count,datatype

#ifdef VAR_LSMPI
      DATATYPE = MPI_CHARACTER
      COUNT = 1
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_charac

    subroutine ls_mpibcast_characV(buffer,nbuf,master)
      implicit none
      integer :: nbuf,master,ierr,count,datatype
      character*(*) :: buffer

#ifdef VAR_LSMPI
      DATATYPE = MPI_CHARACTER
      COUNT = nbuf
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_characV

    subroutine LSMPI_MYFAIL(IERR)
      implicit none
      integer :: IERR,IERR2,IERRCL
      CHARACTER(len=40) :: ERRBUF
#ifdef VAR_LSMPI
      CALL MPI_ERROR_CLASS(IERR,IERRCL,IERR2)
      IF (IERRCL.EQ.MPI_SUCCESS) THEN
         ERRBUF = 'No error'
      ELSE IF (IERRCL.EQ.MPI_ERR_BUFFER) THEN
         ERRBUF = 'Invalid buffer pointer'
      ELSE IF (IERRCL.EQ.MPI_ERR_COUNT) THEN
         ERRBUF = 'Invalid count argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_TYPE) THEN
         ERRBUF = 'Invalid datatype argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_TAG) THEN
         ERRBUF = 'Invalid tag argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_COMM) THEN
         ERRBUF = 'Invalid communicator'
      ELSE IF (IERRCL.EQ.MPI_ERR_RANK) THEN
         ERRBUF = 'Invalid rank'
      ELSE IF (IERRCL.EQ.MPI_ERR_REQUEST) THEN
         ERRBUF = 'Invalid request (handle)'
      ELSE IF (IERRCL.EQ.MPI_ERR_ROOT) THEN
         ERRBUF = 'Invalid root'
      ELSE IF (IERRCL.EQ.MPI_ERR_GROUP) THEN
         ERRBUF = 'Invalid group'
      ELSE IF (IERRCL.EQ.MPI_ERR_OP) THEN
         ERRBUF = 'Invalid operation'
      ELSE IF (IERRCL.EQ.MPI_ERR_TOPOLOGY) THEN
         ERRBUF = 'Invalid topology'
      ELSE IF (IERRCL.EQ.MPI_ERR_DIMS) THEN
         ERRBUF = 'Invalid dimension argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_ARG) THEN
         ERRBUF = 'Invalid argument of some other kind'
      ELSE IF (IERRCL.EQ.MPI_ERR_UNKNOWN) THEN
         ERRBUF = 'Unknown error'
      ELSE IF (IERRCL.EQ.MPI_ERR_TRUNCATE) THEN
         ERRBUF = 'Message truncated on receive'
      ELSE IF (IERRCL.EQ.MPI_ERR_OTHER) THEN
         ERRBUF = 'Known error not in this list'
      ELSE IF (IERRCL.EQ.MPI_ERR_INTERN) THEN
         ERRBUF = 'Internal MPI (implementation) error'
      ELSE IF (IERRCL.EQ.MPI_ERR_IN_STATUS) THEN
         ERRBUF = 'Error code is in status'
      ELSE IF (IERRCL.EQ.MPI_ERR_PENDING) THEN
         ERRBUF = 'Pending request'
      ELSE IF (IERRCL.EQ.MPI_ERR_LASTCODE) THEN
         ERRBUF = 'Last error code'
      ELSE
!        something we didn't know ...
         WRITE(6,'(A,I4)') 'MPI error class',IERRCL
      END IF
!
      WRITE(6,'(/A)') ' ERROR in MPI : '//ERRBUF
!
      CALL LSQUIT('Error detected in MPI. Please consult dalton output!',-1)
!
#endif
    end subroutine LSMPI_MYFAIL

  end module lsmpi_mod

    subroutine lsmpi_init
#ifdef VAR_LSMPI
    use lsmpi_mod
    use infpar_module
    integer :: ierr
    
    write(*,*) 'debug:entered lsmpi_init'
    call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, infpar%mynum, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, infpar%nodtot, ierr )

    infpar%master = 0;

    write(*,*) 'debug:mynum',infpar%mynum,infpar%master
    if (infpar%mynum.ne.infpar%master) then

      write(*,*) 'debug: slave entered'
      call lsmpi_slave

    endif
#endif
    end subroutine lsmpi_init 

#ifdef VAR_LSMPI

    subroutine lsmpi_slave
    use infpar_module
    use lsmpi_mod
    integer :: ierr
    character*8 :: job

write(*,*) 'debug:pending job',infpar%mynum
100  call ls_mpibcast(job,8,infpar%master)
!100  call MPI_BCAST(job,8,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
write(*,*) 'debug:job requested',job,infpar%mynum

     select case(job)
       case('LSGETINT');
         call lsmpi_getIntegrals_Slave
       case('QUIT');  
write(*,*) 'quit reached - ending job',infpar%mynum
         call lsmpi_finalize
       case default
         write(*,*) 'SLAVE: Job not recognized. Quit!'
         call lsmpi_finalize
     end select

write(*,*) 'job finished - pending job',infpar%mynum

     goto 100

    end subroutine lsmpi_slave
#endif

    subroutine lsmpi_finalize
    use lsmpi_mod
    use infpar_module
#ifdef VAR_LSMPI
    integer :: ierr
    
     if (infpar%mynum.eq.infpar%master) &
    &    call MPI_BCAST('QUIT    ',8,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

     call MPI_FINALIZE(ierr)

     !stop all slaves
     if (infpar%mynum.ne.infpar%master) STOP
#endif
 
    end subroutine lsmpi_finalize

    

