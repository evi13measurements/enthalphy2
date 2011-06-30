!> @file
!> Contains module for MPI integral evaluation settings
module infpar_module
#if defined (VAR_MPI) || defined (VAR_LSMPI)

 type infpar_struct
   integer :: nodtot, nodeid(129), ncode,  iprpar, mtottk, ntask,&
             &nfmat,  ndegdi, master, mynum,  mytid,  timing,   slave,&
             &debug
   character(20) :: nodnam(128), myname
 end type infpar_struct

 type(infpar_struct) :: infpar

#endif

#ifdef VAR_MPI

contains
    
    subroutine infpar_init
    integer :: iprtyp = 10, iprint=0

    !initalize infpar structure
    call infpar90_init(infpar)

    !tell slaves to run infpar_slave_init()
    call mpixbcast(iprtyp,1,'INTEGER',infpar%master)
    call mpixbcast(iprint,1,'INTEGER',infpar%master)

    end subroutine


#endif

end module infpar_module

#ifdef VAR_MPI
!   These routines heve to be in scope of F77, C and C++
    subroutine infpar_slave_init()
      use infpar_module
      call infpar90_init(infpar)
    end subroutine 

!#if 0
!   subroutine mat_test_slave
!     use infpar_module
!     use matrix_operations
!
!     Type(Matrix) :: A
!     integer      :: nrc
!
!    !if (infpar%mynum.ne.1) return
!
!     call mat_select_type(mtype_sparse_block)
!     call mpixbcast(nrc,1,'INTEGER',infpar%master)
!    !call mpixrecv(nrc,1,'INTEGER',infpar%master,0)
!
!     call mat_init(A,nrc,nrc);
!
!     call mat_mpixbcast(A, infpar%master)
!    !call mat_mpixrecv(A,infpar%master,1)
!
!     print *, "Matrix norm on slave:", infpar%mynum, mat_sqnorm2(A)
! 
!   end subroutine
!#endif
#endif
