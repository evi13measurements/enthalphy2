!> @file
!> Interface subroutines for exchange-correlation contributions

!> \brief Calculates the xc contribution to the Kohn-Sham matrix
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_get_xc_Fock_mat(LUPRI,LUERR,SETTING,nbast,D,F,EDFT)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use memory_handling
  use lstiming
  use IIDFTKSM
  use DFT_type
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D
!> The Kohn-Sham matrix
TYPE(MATRIX)          :: F
!> The xc contribution to the energy
REAL(REALK)           :: EDFT
!
TYPE(MATRIX)          :: temp
INTEGER               :: i,j,ndmat
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK)           :: TS,TE
REAL(REALK)           :: CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
REAL(REALK)   :: DUMMY(1,1)

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
ndmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
call mem_alloc(Dmat,nbast,nbast,ndmat)
call mem_alloc(DFTDATA%FKSM,nbast,nbast,ndmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*ndmat)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1D0,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2.d0,Dmat(:,:,1),1)
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_KSM(SETTING,LUPRI,1,nbast,ndmat,Dmat,DFTDATA,EDFT)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

IF(SETTING%SCHEME%DODISP) THEN 
    ! add empirical dispersion correction \Andreas Krapp
   CALL II_DFTDISP(SETTING,DUMMY,1,1,0,LUPRI,1)
   DFTDATA%ENERGY = DFTDATA%ENERGY + SETTING%EDISP
   EDFT = EDFT + SETTING%EDISP
ENDIF

!WARNING: 
! For a closed shell molecule calculated using an unrestricted and a 
! closed shell aproach  would fullfill
! F(from closed shell) = F_alpha(from onres) + F_beta(from onres) 
! but for some reason this is not what they want in SCF-loop so we 
! multiply the unrestriced result with 2

call mat_init(temp,F%nrow,F%ncol)
IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,1),1,temp%elms,1)
   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,2),1,temp%elmsb,1)
   CALL mat_DAXPY(1.0d0,temp,F)
ELSE !CLOSED_SHELL
   CALL mat_set_from_full(DFTDATA%FKSM(:,:,1),1.D0,temp,'XCmat')
   CALL mat_DAXPY(0.5d0,temp,F)
ENDIF

call mat_free(temp)
CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_Fock_mat is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_Fock_mat is  ',WALLTIME,LUPRI)

call mem_dealloc(Dmat)
call mem_dealloc(DFTDATA%FKSM)

END SUBROUTINE II_get_xc_Fock_mat

!> \brief Calculates the xc contribution to the molecular gradient
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_get_xc_geoderiv_molgrad(LUPRI,LUERR,SETTING,nbast,D,grad,natoms)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use lstiming
  use IIDFTKSM
  use DFT_type
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D
!> The xc contribution to the molecular grad
REAL(REALK)           :: GRAD(3,natoms)
!> Number of atoms
INTEGER               :: natoms
!
INTEGER               :: i,j
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK)           :: TS,TE
REAL(REALK)           :: CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ndmat,norb,norbitals

ndmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
call mem_alloc(Dmat,nbast,nbast,ndmat)
call mem_alloc(DFTDATA%grad,3,natoms)
DFTDATA%grad = 0.d0

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
DFTDATA%natoms = natoms

call mem_alloc(DFTDATA%orb2atom,nbast)
nOrbitals = 0
DO I=1,nAtoms
   nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbREG
   DO J = nOrbitals+1,nOrbitals+nOrb
      DFTDATA%orb2atom(J)=I
   ENDDO
   nOrbitals = nOrbitals + nOrb
ENDDO
IF(nOrbitals .NE. nbast)&
     &CALL LSQUIT('mismatch in orbital dimension in II_get_xc_geoderiv_molgrad',-1)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('not working - not implemented',-1)
   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1D0,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2.d0,Dmat(:,:,1),1)
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_dft_geoderiv_molgrad(setting,LUPRI,1,nbast,Dmat,DFTDATA)
CALL LSTIMER('geoderiv_molgrd',TS,TE,LUPRI)
GRAD = DFTDATA%grad
CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_geoderiv_molgrad is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_geoderiv_molgrad is  ',WALLTIME,LUPRI)

call mem_dealloc(DFTDATA%orb2atom)
call mem_dealloc(Dmat)
call mem_dealloc(DFTDATA%grad)

END SUBROUTINE II_get_xc_geoderiv_molgrad

!> \brief Calculates the xc contribution to the linear response
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_get_xc_linrsp(LUPRI,LUERR,SETTING,nbast,b,D,G,nbmat)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use memory_handling
  use lstiming
  use IIDFTKSM
  use DFT_type
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The b matrix G(b)
TYPE(MATRIX)          :: b(nbmat)
!> The density matrix
TYPE(MATRIX)          :: D
!> The xc cont to the linear response
TYPE(MATRIX)          :: G(nbmat)
!> Number of b mat
INTEGER               :: nbmat
!
INTEGER               :: i,j,ndmat
TYPE(MATRIX)          :: temp
TYPE(DFTDATATYPE)  :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:)
REAL(REALK)           :: TS,TE
REAL(REALK)           :: CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ibmat

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
ndmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
DFTDATA%nbmat = nbmat
call mem_alloc(Dmat,nbast,nbast)
call mem_alloc(DFTDATA%FKSM,nbast,nbast,nbmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*nbmat)

!WARNING: the densitymatrix for an unres calc fullfill Tr(DS)= N   
!while the closed shell calculation fullfill Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
!   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)   
!   call mem_alloc(DFTDATA%BMAT,nbast,nbast,2*nbmat)
!   CALL DCOPY(b%nrow*b%ncol,b%elms,1,DFTDATA%bmat(:,:,1),1)
!   CALL DCOPY(b%nrow*b%ncol,b%elmsb,1,DFTDATA%bmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1D0,Dmat(:,:))
   CALL DSCAL(nbast*nbast,2.d0,Dmat(:,:),1)
   call mem_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   DO IBMAT=1,nbmat
      call mat_to_full(b(IBMAT),1D0,DFTDATA%bmat(:,:,IBMAT))
   ENDDO
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_LINRSP(SETTING,LUPRI,1,nbast,ndmat,Dmat,DFTDATA)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,1),1,G%elms,1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,2),1,G%elmsb,1)
ELSE !CLOSED_SHELL
   DO IBMAT=1,nbmat
      CALL mat_set_from_full(DFTDATA%FKSM(:,:,IBMAT),1.D0,G(IBMAT),'XCmat')
      CALL mat_scal(0.5d0,G(IBMAT))
   ENDDO
ENDIF

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_linrsp is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_linrsp is  ',WALLTIME,LUPRI)

call mem_dealloc(Dmat)
call mem_dealloc(DFTDATA%FKSM)
call mem_dealloc(DFTDATA%BMAT)

END SUBROUTINE II_get_xc_linrsp

!> \brief Calculates the xc contribution to the quadratic response
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_quadrsp(LUPRI,LUERR,SETTING,nbast,b,c,D,T)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use memory_handling
  use lstiming
  use IIDFTKSM
  use DFT_type
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The b matrix T(b,c)
TYPE(MATRIX)          :: b
!> The b matrix T(b,c)
TYPE(MATRIX)          :: c
!> The density matrix
TYPE(MATRIX)          :: D
!> The xc cont to the quadratic response
TYPE(MATRIX)          :: T
!
INTEGER               :: i,j,ndmat,nbmat
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK)           :: TS,TE
REAL(REALK)           :: CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
ndmat = 1
nbmat = 2
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
DFTDATA%nbmat = nbmat
call mem_alloc(Dmat,nbast,nbast,ndmat)
call mem_alloc(DFTDATA%FKSM,nbast,nbast,ndmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*ndmat)

!WARNING: the densitymatrix for an unres calc fullfill Tr(DS)= N   
!while the closed shell calculation fullfill Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
!   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)   
!   call mem_alloc(DFTDATA%BMAT,nbast,nbast,2*nbmat)
!   CALL DCOPY(b%nrow*b%ncol,b%elms,1,DFTDATA%bmat(:,:,1),1)
!   CALL DCOPY(b%nrow*b%ncol,b%elmsb,1,DFTDATA%bmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1D0,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2.d0,Dmat(:,:,1),1)
   call mem_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   call mat_to_full(b,1D0,DFTDATA%bmat(:,:,1))
   call mat_to_full(c,1D0,DFTDATA%bmat(:,:,2))
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_QUADRSP(SETTING,LUPRI,1,nbast,ndmat,Dmat,DFTDATA)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,1),1,T%elms,1)
   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,2),1,T%elmsb,1)
ELSE !CLOSED_SHELL
   CALL mat_set_from_full(DFTDATA%FKSM(:,:,1),1.D0,T,'XCmat')
   CALL mat_scal(0.5d0,T)
ENDIF

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_quadrsp is ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_quadrsp is ',WALLTIME,LUPRI)

call mem_dealloc(Dmat)
call mem_dealloc(DFTDATA%FKSM)
call mem_dealloc(DFTDATA%BMAT)

END SUBROUTINE II_get_xc_quadrsp

!> \brief Calculates the xc contribution to the magnetic derivative kohn-sham matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,SETTING,nbast,D,F)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use memory_handling
  use lstiming
  use IIDFTKSM
  use DFT_type
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D
!> The xc cont to the magnetic deriv of F
TYPE(MATRIX)          :: F(3) !x,y and z components
INTEGER               :: i,j,ndmat
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
ndmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
call mem_alloc(Dmat,nbast,nbast,ndmat)
call mem_alloc(DFTDATA%FKSM,nbast,nbast,3*ndmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*3*ndmat)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('II_get_xc_magderiv_kohnsham_mat not implemented for unrestricted yet',lupri)
!   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
!   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1D0,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2.d0,Dmat(:,:,1),1)
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_dft_magderiv_kohnsham_mat(SETTING,LUPRI,1,nbast,ndmat,Dmat,DFTDATA)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

!WARNING: 
! For a closed shell molecule calculated using an unrestricted and a 
! closed shell aproach  would fullfill
! F(from closed shell) = F_alpha(from onres) + F_beta(from onres) 
! but for some reason this is not what they want in SCF-loop so we 
! multiply the unrestriced result with 2

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('II_get_xc_magderiv_kohnsham_mat not implemented for unrestricted yet',lupri)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,1),1,temp%elms,1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,2),1,temp%elmsb,1)
!   CALL mat_DAXPY(1.0d0,temp,F)
ELSE !CLOSED_SHELL
   do I=1,3*ndmat
      CALL mat_set_from_full(DFTDATA%FKSM(:,:,I),1.D0,F(I),'XCmat')
      CALL mat_scal(0.5d0,F(I))
   enddo
ENDIF

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_magderiv_kohnsham_mat is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_magderiv_kohnsham_mat is  ',WALLTIME,LUPRI)

call mem_dealloc(Dmat)
call mem_dealloc(DFTDATA%FKSM)

END SUBROUTINE II_get_xc_magderiv_kohnsham_mat

!> \brief Calculates the xc contribution to the mag derivative of the linear rsp
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_magderiv_linrsp(LUPRI,LUERR,SETTING,nbast,b,D,G,nbmat)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use memory_handling
  use lstiming
  use IIDFTKSM
  use DFT_type
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> the B mat, perturb density
TYPE(MATRIX)          :: b(nbmat)
!> The density matrix
TYPE(MATRIX)          :: D
!> The xc cont to the magnetic deriv of G(b)
TYPE(MATRIX)          :: G(3*nbmat)
!> number of B matrices
INTEGER               :: nbmat
INTEGER               :: i,j,ndmat
TYPE(MATRIX)          :: temp
TYPE(DFTDATATYPE)  :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ibmat

WRITE(lupri,*)'STARTING II_get_xc_magderiv_linrsp'
CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
ndmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
DFTDATA%nbmat = nbmat
call mem_alloc(Dmat,nbast,nbast)
call mem_alloc(DFTDATA%FKSM,nbast,nbast,3*nbmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*3*nbmat)
DFTDATA%dosympart = .TRUE. !should be a test
IF(DFTDATA%dosympart)THEN
   call mem_alloc(DFTDATA%FKSMS,nbast,nbast,3*nbmat)
   CALL LS_DZERO(DFTDATA%FKSMS,nbast*nbast*3*nbmat)
ENDIF
!WARNING: the densitymatrix for an unres calc fullfill Tr(DS)= N   
!while the closed shell calculation fullfill Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
!   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)   
!   call mem_alloc(DFTDATA%BMAT,nbast,nbast,2*nbmat)
!   CALL DCOPY(b%nrow*b%ncol,b%elms,1,DFTDATA%bmat(:,:,1),1)
!   CALL DCOPY(b%nrow*b%ncol,b%elmsb,1,DFTDATA%bmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1D0,Dmat(:,:))
   CALL DSCAL(nbast*nbast,2.d0,Dmat(:,:),1)
   call mem_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   DO IBMAT=1,nbmat
      call mat_to_full(b(IBMAT),1D0,DFTDATA%bmat(:,:,IBMAT))
   ENDDO
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_MAGDERIV_LINRSP(SETTING,LUPRI,1,nbast,ndmat,Dmat,DFTDATA)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,1),1,G%elms,1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,2),1,G%elmsb,1)
ELSE !CLOSED_SHELL
   DO IBMAT=1,nbmat*3
      CALL mat_set_from_full(DFTDATA%FKSM(:,:,IBMAT),1.D0,G(IBMAT),'XCmat')
      CALL mat_scal(1.d0,G(IBMAT))
   ENDDO
   do I=1,nbmat
      WRITE(LUPRI,*)'ASYM PART II_X'
      call mat_print(G(1+(I-1)*3),1,nbast,1,nbast,lupri)
      WRITE(LUPRI,*)'ASYM PART II_Y'
      call mat_print(G(2+(I-1)*3),1,nbast,1,nbast,lupri)
      WRITE(LUPRI,*)'ASYM PART II_Z'
      call mat_print(G(3+(I-1)*3),1,nbast,1,nbast,lupri)
   enddo
   IF(DFTDATA%dosympart)THEN
      call mat_init(temp,nbast,nbast)
      DO IBMAT=1,nbmat*3
         CALL mat_set_from_full(DFTDATA%FKSMS(:,:,IBMAT),1.D0,temp,'XCmat')
         WRITE(LUPRI,*)'SYMMETRIC (0.5) PART II_',IBMAT
         call mat_print(temp,1,nbast,1,nbast,lupri)
         CALL mat_daxpy(0.5d0,temp,G(IBMAT))
      ENDDO
      call mat_free(temp)
   ENDIF
ENDIF

do I=1,nbmat
   WRITE(LUPRI,*)'THE II_X magderiv linrsp component'
   call mat_print(G(1+(I-1)*3),1,nbast,1,nbast,lupri)
   WRITE(LUPRI,*)'THE II_Y magderiv linrsp component'
   call mat_print(G(2+(I-1)*3),1,nbast,1,nbast,lupri)
   WRITE(LUPRI,*)'THE II_Z magderiv linrsp component'
   call mat_print(G(3+(I-1)*3),1,nbast,1,nbast,lupri)
enddo

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_linrsp is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_linrsp is  ',WALLTIME,LUPRI)

call mem_dealloc(Dmat)
call mem_dealloc(DFTDATA%FKSM)
IF(DFTDATA%dosympart)call mem_dealloc(DFTDATA%FKSMS)
call mem_dealloc(DFTDATA%BMAT)
WRITE(lupri,*)'END OF II_get_xc_magderiv_linrsp'

END SUBROUTINE II_get_xc_magderiv_linrsp

!> \brief Calculates the xc contribution to the geo derivative of the KohnSham matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_geoderiv_FxDgrad(LUPRI,LUERR,SETTING,nbast,D,b,grad,natoms)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use lstiming
  use IIDFTKSM
  use DFT_type
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D
!> the B mat, perturb density
TYPE(MATRIX)          :: b
!> The xc cont to the geomderiv of F 
REAL(REALK)           :: GRAD(3,natoms)
!> number of atoms
INTEGER               :: natoms
!
INTEGER               :: i,j
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ndmat,norb,norbitals,nbmat,ibmat

ndmat = 1
nbmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
call mem_alloc(Dmat,nbast,nbast,ndmat)
call mem_alloc(DFTDATA%grad,3,natoms)
DFTDATA%grad = 0.d0

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
DFTDATA%natoms = natoms

call mem_alloc(DFTDATA%orb2atom,nbast)
nOrbitals = 0
DO I=1,nAtoms
   nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbREG
   DO J = nOrbitals+1,nOrbitals+nOrb
      DFTDATA%orb2atom(J)=I
   ENDDO
   nOrbitals = nOrbitals + nOrb
ENDDO
IF(nOrbitals .NE. nbast)&
     &CALL LSQUIT('mismatch in orbital dimension in II_get_xc_geoderiv_FxDgrad',-1)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('not working - not implemented',-1)
   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1D0,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2.d0,Dmat(:,:,1),1)
   call mem_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   DO IBMAT=1,nbmat
      call mat_to_full(b,1D0,DFTDATA%bmat(:,:,IBMAT))
   ENDDO
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_geoderiv_kohnsham_mat(setting,LUPRI,1,nbast,1,Dmat,DFTDATA)
CALL LSTIMER('geoderiv_FxDgrd',TS,TE,LUPRI)
GRAD = DFTDATA%grad
CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_geoderiv_FxDgrad is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_geoderiv_FxDgrad is  ',WALLTIME,LUPRI)

call mem_dealloc(DFTDATA%orb2atom)
call mem_dealloc(Dmat)
call mem_dealloc(DFTDATA%grad)
call mem_dealloc(DFTDATA%BMAT)

END SUBROUTINE II_get_xc_geoderiv_FxDgrad

!> \brief Calculates the xc contribution to the geo derivative of the linear response
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_geoderiv_GxDgrad(LUPRI,LUERR,SETTING,nbast,D,a,b,grad,natoms)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use lstiming
  use IIDFTKSM
  use DFT_type
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D
!> the A mat, perturb density A*G(B)
TYPE(MATRIX)          :: a
!> the B mat, perturb density A*G(B)
TYPE(MATRIX)          :: b
!> The xc cont to the geomderiv of G 
REAL(REALK)           :: GRAD(3,natoms)
!> number of atoms
INTEGER               :: natoms
!
INTEGER               :: i,j
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ndmat,norb,norbitals,nbmat,ibmat

ndmat = 1
nbmat = 2
DFTDATA%nbmat = nbmat
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
call mem_alloc(Dmat,nbast,nbast,ndmat)
call mem_alloc(DFTDATA%grad,3,natoms)
DFTDATA%grad = 0.d0

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
DFTDATA%natoms = natoms

call mem_alloc(DFTDATA%orb2atom,nbast)
nOrbitals = 0
DO I=1,nAtoms
   nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbREG
   DO J = nOrbitals+1,nOrbitals+nOrb
      DFTDATA%orb2atom(J)=I
   ENDDO
   nOrbitals = nOrbitals + nOrb
ENDDO
IF(nOrbitals .NE. nbast)&
     &CALL LSQUIT('mismatch in orbital dimension in II_get_xc_geoderiv_FxDgrad',-1)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('not working - not implemented',-1)
   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1D0,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2.d0,Dmat(:,:,1),1)
   call mem_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   call mat_to_full(a,1D0,DFTDATA%bmat(:,:,1))
   call mat_to_full(b,1D0,DFTDATA%bmat(:,:,2))
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_geoderiv_linrspgrad(setting,LUPRI,1,nbast,1,Dmat,DFTDATA)
CALL LSTIMER('geoderiv_FxDgrd',TS,TE,LUPRI)
GRAD = DFTDATA%grad
CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_geoderiv_FxDgrad is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_geoderiv_FxDgrad is  ',WALLTIME,LUPRI)

call mem_dealloc(DFTDATA%orb2atom)
call mem_dealloc(Dmat)
call mem_dealloc(DFTDATA%grad)
call mem_dealloc(DFTDATA%BMAT)

END SUBROUTINE II_get_xc_geoderiv_GxDgrad
