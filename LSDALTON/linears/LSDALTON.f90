!> @file 
!> Contains main SCF driver, some module wrappers and miscellaneous.

!> \brief Driver for stand-alone f90 linear scaling SCF.
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2008-10-26
SUBROUTINE lsdalton
  use configuration
  use direct_dens_util
  use files
  use initial_guess
  use linsca_debug
  use lstiming
  USE scf_stats
  use daltoninfo       
  uSE density_optimization
  use lsdalton_fock_module
  use ks_settings
  use IIDFTINT
! DEC 
!  use DEC_dalton_mod
! PROPERTIES SECTION
  use lsdalton_rsp_mod
! GEOMETRY OPTIMIZATION
  use optimization_input
  implicit none
  TYPE(lsitem),target :: ls
  real(realk)         :: Tstart,Tend,t1,t2 !,ten,tstr,E,gradnrm
  integer             :: nbast !,lu_error,lu_output
!  TYPE(util_HistoryStore) :: queue
!  TYPE(modFIFO)       :: fifoqueue
  TYPE(Matrix)         :: F,D, CMO
  Type(Matrix), target :: H1,S
!  CHARACTER(len=80)   :: BASISSETNAME
!  CHARACTER(len=9)    :: BASISLABEL
  LOGICAL             :: isdft,restart_from_dens,isgga, do_decomp!,isminimum
  integer             :: matmul1, matmul2, matmultot,idum,ldum,restart_lun, lun
  REAL(REALK)         :: exchange_fac,Epotnuc, mx
! Energy
  REAL(REALK)         :: E
!
  integer,allocatable :: degeneracy(:),orbitals(:)
  integer             :: nshell,natoms,I,J,bast,idummy, lupri, luerr
  type(configItem)    :: config
  logical             :: doDFT,mem_monitor
  real(realk), allocatable :: eival(:)

  REAL(REALK) :: TIMSTR,TIMEND
  real(realk) :: DUMMY(1,1)

! Entered LSDALTON
!  cfg_lsdalton = .true.

! Initializations 
!

! Read DALTON.INP FILE,MOLECULE.INP AND BASISSET FILES TO SET UP
! THE  DALTONINPUT STRUCTURE
  LUPRI=-1
  LUERR=-1
  CALL LSOPEN(LUPRI,'LSDALTON.OUT','NEW','FORMATTED')
  CALL LSOPEN(LUERR,'LSDALTON.ERR','UNKNOWN','FORMATTED')

  call LSTIMER('START',t1,t2,LUPRI)

  CALL PRINT_INTRO(LUPRI)

  CALL LSTIMER('START ',TIMSTR,TIMEND,lupri)

  call mat_no_of_matmuls(no_of_matmuls)
  ! READ *LINSCA SECTION UNDER **WAVEFUNCTION AND INIT SCF_CONFIG
  call config_set_default_config(config)
  ! Setting default starting-guess to ATOMS instead of HUCKEL
  config%opt%cfg_start_guess = 'ATOMS'

  call config_read_input(config,lupri,luerr,.FALSE.)
  ls%input%dalton = config%integral

  doDFT = config%opt%calctype.EQ.config%opt%dftcalc
  call ls_init(ls,lupri,luerr,nbast,config%integral,doDFT,.true.)
  call set_final_config_and_print(lupri,config,ls)

  ! eventually empirical dispersion correction in case of dft
  CALL II_DFTDISP(LS%SETTING,DUMMY,1,1,0,LUPRI,1)

  call mat_pass_info(LUPRI,config%opt%info_matop,mem_monitor)
  config%av%dsm_history_size = config%av%cfg_settings(config%av%CFG_SET_type)%max_history_size
  config%av%diis_history_size = config%av%dsm_history_size
  config%av%ediis_history_size = config%av%dsm_history_size

  CALL LSTIMER('*INPUT',TIMSTR,TIMEND,lupri)
  !lupri is now a global variable corresponding to lupri

#if defined(VAR_INT64)
  if (bit_size(nbast)==64) then
     write(lupri,*) 'Using 64-bit integers!'
  endif
#endif

! Kasper K, avoid Hartree-Fock related calculations for DEC calculation if requested

    if ((config%opt%cfg_start_guess == 'TRILEVEL')&
      &.or.(config%opt%cfg_start_guess == 'ATOMS')&
      &.or.config%decomp%cfg_gcbasis) then
!      allocate(ls)
      call trilevel_basis(config%opt,ls)
      CALL LSTIMER('*ATOM ',TIMSTR,TIMEND,lupri)
    endif


  if (config%opt%cfg_prefer_bsm) then
#ifdef HAVE_BSM
     CALL ls_initbsm(ls%input%BASIS%REGULAR,ls%input)
#else
     call lsquit('.BLOCK requested but BSM is not compiled',lupri)
#endif
  endif

  if (config%sparsetest) then
    call sparsetest(ls%setting, lupri)
  endif

  CALL mat_init(D,nbast,nbast)
  CALL mat_init(F,nbast,nbast)
  CALL mat_init(H1,nbast,nbast)
  CALL mat_init(S,nbast,nbast)

    CALL II_get_overlap(lupri,luerr,ls%setting,S)

!  CALL LSHEADER(6,'Overlap matrix - Test print')
!  CALL mat_print(S,1,S%nrow,1,S%ncol,6)
!  just for test and then compute overlap integrals
!  Call set_pbc_molecules(ls%input,ls%setting,lupri,luerr)
  CALL LSTIMER('*S    ',TIMSTR,TIMEND,lupri)
  !write(lupri,*) 'QQQ New  S:',mat_trab(S,S)
  CALL II_get_h1(lupri,luerr,ls%setting,H1)
  CALL LSTIMER('*H1   ',TIMSTR,TIMEND,lupri)
  !write(lupri,*) 'QQQ New  H1:',mat_trab(H1,H1)

  !data to pass down to fck_get_fock subroutine
  lsint_fock_data%ls => ls
  lsint_fock_data%H1 => H1
  lsint_fock_data%lupri = lupri
  lsint_fock_data%luerr = luerr

  call get_initial_dens(H1,S,D,ls,config)

  if (mem_monitor) then
     write(lupri,*)
     WRITE(LUPRI,'("Max no. of matrices allocated in Level 2 / get_initial_dens: ",I10)') max_no_of_matrices
     max_no_of_matrices = no_of_matrices
  endif
  !write (lupri,*) 'Start density:'
  !call MAT_PRINT(D, 1, D%nrow, 1, D%ncol, LUPRI)

  CALL LSTIMER('*START',TIMSTR,TIMEND,lupri)

  if (config%opt%cfg_incremental) then
    call ks_init_incremental_fock(nbast)
  endif

  do_decomp =(config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
            & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
            & config%decomp%cfg_check_converged_solution .or. &
            & config%decomp%cfg_rsp_nexcit > 0) 

  if (do_decomp) then
     call decomp_init(nbast,config%decomp)
     call dd_mat_eigenvalues_to_aux(config%decomp%cfg_unres,S)
     call decomposition(config%decomp)
  else if (config%opt%cfg_start_guess == 'TRILEVEL') then
     call mat_free(config%decomp%lcv_CMO)
  endif

  if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then !FIXME: put this somewhere else!
     call mat_init(config%av%Fprev,nbast,nbast)
     call mat_init(config%av%Dprev,nbast,nbast)
  endif

  call scfloop(H1,F,D,S,E,ls,config)

  !lcm basis
  if (config%decomp%cfg_lcm) then
  ! get orbitals
   call mat_init(Cmo,nbast,nbast)
   allocate(eival(nbast))
   call mat_diag_f(F,S,eival,Cmo)
   deallocate(eival)

  ! localize orbitals
   call leastchange_lcm(config%decomp,Cmo,config%decomp%nocc,ls)

  ! write gc basis
   call write_lsitem_to_disk(ls)

  ! write orbitals
   lun = -1
   CALL LSOPEN(lun,'lcm_orbitals.u','new','UNFORMATTED')
   call mat_write_to_disk(lun,Cmo)
   call LSclose(LUN,'KEEP')

   call leastchangeOrbspreadStandalone(mx,ls,Cmo,config%decomp%lupri,config%decomp%luerr)
   write(*,*) 'Orbspread standalone LCM: ', mx

   ! set basis to CMO
   call mat_assign(config%decomp%U_inv,Cmo)
   call mat_mul(config%decomp%U_inv,S,'t','n',1d0,0d0,config%decomp%U)


  ! free Cmo
   call mat_free(Cmo)
  endif

!
! Optimization
!
   if (config%optinfo%optimize .EQV. .TRUE.) then
      CALL LS_runopt(E,config,config%optinfo,H1,F,D,S,ls,lupri,luerr)
   endif

  call config_shutdown(config)

!write(lupri,*) 'mem_allocated_integer, max_mem_used_integer', mem_allocated_integer, max_mem_used_integer

  call lsdalton_response(lupri,ls,config,F,D,S,&
       &config%opt%calctype .EQ. config%opt%dftcalc)

! PROPERTIES SECTION
!
   if (config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. & 
     & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
     & config%decomp%cfg_check_converged_solution .or. config%decomp%cfg_rsp_nexcit > 0) then   
      call get_oao_transformed_matrices(config%decomp,F,D)
   endif
  if (do_decomp) then
     call decomp_shutdown(config%decomp)
     !call dd_shutdown(config%decomp%cfg_unres)
  endif

  CALL LSTIMER('*SCF  ',TIMSTR,TIMEND,lupri)
  WRITE(lupri,*)

  !call mat_no_of_matmuls(matmultot)
  !WRITE(lupri,'("Total no. of matmuls in SCF optimization: ",I10)') matmultot

  if ((config%opt%cfg_start_guess.eq.'TRILEVEL')&
    &.or.(config%opt%cfg_start_guess.eq.'ATOMS')&
    &.or.config%decomp%cfg_gcbasis) then
  !  CALL trilevel_shutdown
   !DEALLOCATE(ls)
  endif

  if (config%opt%cfg_incremental) call ks_free_incremental_fock()
!
! FINALIZE SECTION - releases memory n'stuff
!
  
  CALL mat_free(H1)
  CALL mat_free(F)
  CALL mat_free(D)
  CALL mat_free(S)



  call mat_no_of_matmuls(no_of_matmuls)
  WRITE(LUPRI,'("Total no. of matmuls used:                ",I10)') no_of_matmuls
  WRITE(LUPRI,'("Total no. of Fock/KS matrix evaluations:  ",I10)') ls%input%nfock
  if (mem_monitor) WRITE(LUPRI,'("Max no. of matrices allocated in Level 3: ",I10)') max_no_of_matrices
  call stats_mem(config%lupri)

  call LSTIMER('LSDALTON',t1,t2,LUPRI)
  !call CPU_TIME(tend)      redundant! /Stinne 13/10-2010
  !WRITE(lupri,*) "*****************************************************"
  !WRITE(lupri,*) "**     CPU-TIME USED IN LSDALTON: ",tend-tstart,"   **"
  !WRITE(lupri,*) "*****************************************************"
  call ls_free(ls)
  !CALL QEXIT('LSDALTON')

END SUBROUTINE LSDALTON


!> \brief Print author list for stand-alone f90 linear scaling SCF to LSDALTON.OUT
!> \author T. Kjaergaard, S. Reine
!> \date 2008-10-26
SUBROUTINE print_intro(lupri)
implicit none
!> Logical unit number for file DALTON.OUT
integer,intent(in)        :: lupri
integer                   :: I

WRITE(LUPRI,*)' '
WRITE(LUPRI,*)'    ******************************************************************    '
WRITE(LUPRI,*)'    **********  LSDALTON - An electronic structure program  **********    '
WRITE(LUPRI,*)'    ******************************************************************    '
WRITE(LUPRI,*)' '


         WRITE (LUPRI,'(5X,A)')&
     &  ' This is output from LSDALTON (Release Dalton2011)',&
     &  '       Authors:                                  ',&
     &  ' Sonia Coriani,         University of Trieste,        Italy  ',&
     &  ' Trygve Helgaker,       University of Oslo,           Norway ',&
     &  ' Stinne Hoest,          Aarhus University,            Denmark',&
     &  ' Branislav Jansik,      Aarhus University,            Denmark',&
     &  ' Poul Joergensen,       Aarhus University,            Denmark',&
     &  ' Joanna Kauczor,        Aarhus University,            Denmark',&
     &  ' Thomas Kjaergaard,     Aarhus University,            Denmark',&
     &  ' Kasper Kristensen,     Aarhus University,            Denmark',&
     &  ' Jeppe Olsen,           Aarhus University,            Denmark',&
     &  ' Simen Reine,           University of Oslo,           Norway ',&
     &  ' Pawel Salek,           KTH Stockholm,                Sweden ',&
     &  ' Andreas Thorvaldsen,   University of Tromsoe,        Norway ',&
     &  ' Lea Thoegersen,        Aarhus University,            Denmark',&
     &  ' Vladimir Rybkin,       University of Oslo,           Norway ',&
     &  ' Vebjoern Bakken,       University of Oslo,           Norway ',&
     &  ' Mark Watson,           University of Oslo,           Norway ',&
     &  ' Andreas Krapp,         University of Oslo,           Norway '
!         WRITE (LUPRI,'(/1X,69A1/)') ('-',I=1,69)

         WRITE (LUPRI,'(5X,A)')&
     &'NOTE:',&
     &' ',&
     &'This is an experimental code for the evaluation of molecular',&
     &'properties using SCF wave functions. The authors',&
     &'accept no responsibility for the performance of the code or',&
     &'for the correctness of the results.',&
     &' ',&
     &'The code (in whole or part) is provided under a licence and',&
     &'is not to be reproduced for further distribution without',&
     &'the written permission of the authors or their representatives.',&
     &' ',&
     &'See the home page "http://daltonprogram.org"',&
     &'for further information.',&
     &' ',&
     &'If results obtained with this code are published,',&
     &'an appropriate citation would be:',&
     &' ',&
     &'"LSDalton, a molecular electronic structure program, Release 2011,',&
     &' written by <INSERT AUTHOR LIST>"'
        write(lupri,*)

END SUBROUTINE

#ifdef HAVE_BSM 
!> \brief Initialize the Block-Sparse Matrix module
!> \author B. Jansik
!> \date 2008-10-26
!> \param basis contains information about the basis set
!> \param input contains information read from file DALTON.INP
!>
!> initbsm is a routine to intialize the blocked sparse library.  The
!> library needs to know the position of all atoms and the indexes of
!> the related basis function blocks in the matrices to create proper
!> permutations of the blocks. Epsilon is the truncation
!> threshold. MAXBLOCKSIZE is the block size optimal for used blas
!> implementation. It should usually be around 100-150.  TRACETHR is
!> the trace convergence threshold when the eigenvalues are supposed
!> to be well separated and after which the algorithm switches to a
!> "tightening" mode.  NTIGHT is the number of double-iterations made
!>  in the tightening mode.
!>
SUBROUTINE ls_initbsm(basis,input)
  use typedef
  implicit none
  type(basissetinfo) :: basis
  type(daltoninput),intent(in)   :: input
  integer             :: natoms, i
  real(realk),PARAMETER :: EPSILON = 1D-7
  real(realk),PARAMETER :: TRACETHR = 5D-3
  INTEGER,PARAMETER     :: MAXBLOCKSIZE = 110
  INTEGER,PARAMETER     :: NTIGHT = 3
  INTEGER, pointer      :: NATMST(:)
  real(realk),allocatable :: X(:), Y(:), Z(:)

  NATOMS = input%MOLECULE%nAtoms

  NATMST => initbsm_setlist(basis,input)

  allocate(X(NATOMS),Y(NATOMS),Z(NATOMS))

  DO I=1,NATOMS
     X(I) = input%MOLECULE%ATOM(I)%CENTER(1)
     Y(I) = input%MOLECULE%ATOM(I)%CENTER(2)
     Z(I) = input%MOLECULE%ATOM(I)%CENTER(3)
  ENDDO

  CALL bsm_lib_init(X,Y,Z,NATMST,nAtoms,MAXBLOCKSIZE,EPSILON,TRACETHR,NTIGHT)

  deallocate (NATMST,X,Y,Z)

contains

!> \brief Initialize the list for the Block-Sparse Matrix module (?) Not really sure what happens here!!
!> \author B. Jansik
!> \date 2008-10-26
!> \param basis Contains information about the basis set
!> \param dalton_inp Contains information read from DALTON.INP
 function initbsm_setlist(basis,input)
  use BUILDAOBATCH
  use typedef
implicit none
type(basissetinfo),intent(in)  :: basis
type(daltoninput),intent(in)    :: input
integer :: nAtoms
integer :: itype, iset, r,icharge
integer, pointer :: initbsm_setlist(:)

 nAtoms = input%MOLECULE%nAtoms

 allocate (initbsm_setlist(nAtoms+1))
 
 initbsm_setlist = 0
 R = basis%labelindex
 do i=1, nAtoms
    IF(R .EQ.0)THEN
       icharge = INT(input%MOLECULE%ATOM(i)%charge) 
       itype = basis%chargeindex(icharge)
    ELSE
       itype = input%MOLECULE%ATOM(i)%IDtype(1)
    ENDIF

    initbsm_setlist(i+1)= &
         &initbsm_setlist(i)  + basis%ATOMTYPE(itype)%ToTnorb
 enddo

 return 
 end function initbsm_setlist

END SUBROUTINE LS_INITBSM
#endif

!> \brief This subroutine calculates number of electron for neutral molecule defined in ls%setting structure
!> \author T. Kjaergaard, S. Reine
!> \date 2008-10-26
!> \param ls Contains information read from DALTON.INP 
!> \param nel Number of electrons
SUBROUTINE get_num_electrons_neutral(nel,ls)
  use typedef
  implicit none
  TYPE(lsitem) :: ls
  integer, intent(out)        :: nel
  integer :: i

  nel=0
  do i=1,ls%setting%MOLECULE(1)%p%nAtoms
    nel = nel + ls%setting%MOLECULE(1)%p%Atom(i)%Charge
  enddo
END SUBROUTINE get_num_electrons_neutral

!> \brief Contains the main SCF loop for LSDALTON
!> \author T. Kjaergaard, S. Reine
!> \date 2008-10-26
!> \param H1 One-electron Hamiltonian
!> \param F Fock/Kohn-Sham matrix
!> \param D Density matrix
!> \param S Overlap matrix
!> \param ls Contains information read from DALTON.INP 
!> \param config Contains all info about configuration/settings for SCF calculation
SUBROUTINE scfloop(H1,F,D,S,E,ls,config)
  use configuration
  use Matrix_Operations
  use direct_dens_util
  USE scf_stats
  use daltoninfo       
  uSE density_optimization
  use density_subspace_min
  use lstiming
  use rsp_util
   implicit none
   
   type(matrix),intent(in)         :: H1
   type(matrix),intent(inout)      :: D, F, S
   type(configItem),intent(inout)  :: config
   integer                         :: nbast
   TYPE(util_HistoryStore)         :: queue
   TYPE(modFIFO)                   :: fifoqueue
   TYPE(Matrix)                    :: grad,tempm1,tempm2,tempm3,tempm4,tempm5 
   real(realk)                     :: E, gradnrm, hessian_eigenval, ehomo, elumo
   integer                         :: iteration, matmul1, matmul2, matmultot,number_atoms, queuesize,nnz
   REAL                            :: tstart, tend, t0, MFLOPS, norm
   logical                         :: energy_converged
   logical     :: file_exists 
   real(realk) :: TSTR, TEN, TIMSTR, TIMEND, t1, t2
   type(matrix) :: x, Dchol !For debug
   type(matrix) :: unitmat
   type(matrix), pointer :: Dpointer, Fpointer
   type(lsitem) :: ls
   LOGICAL :: LSint, dalink, incremental
   real(realk) :: acceptratio, limitratio

   CALL LSTIMER('START',TSTR,TEN,config%LUPRI)
  !INITIALISE SCF CYCLE THINGS 
  nbast = H1%nrow
  acceptratio = config%av%cfg_settings(config%av%CFG_set_type)%min_density_overlap
  limitratio = config%av%cfg_settings(config%av%CFG_set_type)%max_dorth_ratio
  call DOPT_get_density_init(nbast,acceptratio,limitratio,config%diag)
  !call direct_dens_init
  if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
     queuesize = config%av%cfg_settings(config%av%cfg_set_type)%max_history_size
     !FIXME: this queue cannot be put on disk, get_from_modFIFO_disk won't work! 
     call modfifo_init(fifoqueue,queuesize,nbast,.false.)
  else
     call queue_init(config%av,queue)
  endif
  !if (DEBUG_DCHANGE) call debug_dchange_init

  call scf_stats_init(config%opt)

!!!!Moved outside scfloop, needed also for lcv localization
!#if 0
!  if (cfg_density_method == cfg_f2d_direct_dens .or. cfg_density_method == cfg_f2d_arh .or. &
!       & cfg_check_converged_solution .or. cfg_rsp_nexcit > 0) then   
!     call decomp_init(nbast,decomp)
!     call dd_mat_eigenvalues_to_aux(S)
!     call decomposition(decomp)
!  endif
!#endif
  
!
! Starting guess for density matrix - or restart from old
!
!   call get_initial_dens(H1,S,D)
!   call II_get_Fock_mat(newlupri,newluerr,ls,D,F)
    call mat_report_sparsity(S,'S    ',nnz,config%lupri)
    call mat_report_sparsity(D,'AO D ',nnz,config%lupri)
    IF(associated(F%elms))THEN
       call mat_report_sparsity(F,'AO F ',nnz,config%lupri)
    ENDIF
   write(config%lupri,*) 'Relative convergence threshold for solver:', config%solver%cfg_micro_thresh

   CALL LSTIMER('INIT SCF',TSTR,TEN,config%LUPRI)

   call mat_init(grad,nbast,nbast)

!
! SCF iterations
!
   DO iteration = 1, config%opt%cfg_max_linscf_iterations
      CALL LSTIMER('START ',t1,t2,config%lupri)
!     Incremental scheme set for the density-fitting gradient contribution. /SR 2010-10-19
      incremental = config%opt%cfg_incremental .AND. iteration.NE.1
      call II_setIncremental(ls%setting%scheme,incremental)
      if (iteration == 1) then
         !We cannot use DaLink in the 1st iteration - this gives the right energy, but the Fock matrix
         ! is so bad that ARH cannot converge. If DaLink is requested, turn it off and then back on after
         ! the 1st iteration. /Stinne, Thomas, Brano 23/11-2009
         dalink = .false.
         if (ls%setting%scheme%DALINK) then
            ls%setting%scheme%DALINK = .FALSE.
            dalink = .true.
         endif
      endif

      WRITE(config%LUPRI,'("** Get Fock matrix number ",i3)') iteration
      CALL LSTIMER('START ',TIMSTR,TIMEND,config%lupri)
      CALL get_fock(config, fifoqueue, queue, iteration, D, H1, F, E,LSint,&
           &config%lupri,config%opt%luerr,ls)
      if (iteration == 1) then
         !Turn DaLink back on, if requested:
         if (dalink) then
            ls%setting%scheme%DALINK = .true.
         endif
      endif
      CALL LSTIMER('FCK_FO ',TIMSTR,TIMEND,config%LUPRI)

      if (config%solver%step_accepted) CALL get_AO_gradient(F, D, S, grad)

      CALL LSTIMER('G_GRAD',TIMSTR,TIMEND,config%LUPRI)
      gradnrm = sqrt(mat_sqnorm2(grad))

      ! Statistic stuff
      call scf_stats_update(iteration,gradnrm,E,config%opt)
      call scf_stats_debug_mem(config%lupri,iteration)
      !if (DEBUG_DCHANGE) call debug_dchange_update(D,S)

      ! Test for convergence
      IF(gradnrm < config%opt%set_convergence_threshold) then 
         EXIT
      else if ((config%opt%cfg_hesonly .or. config%opt%cfg_diaghesonly) .and. &
            & iteration == 2) then
         write(config%lupri,*) 'You have requested only calculation of first SCF energy'
         write(config%lupri,*) '- now moving on to calculate lowest Hessian eigenvalue for this point.'
         exit
      ENDIF

      WRITE(config%LUPRI,'("** Make average of the last F and D matrices")')
      CALL Density_subspace_minimization(config, fifoqueue, queue, E, S, H1, grad, F, D, iteration)
      CALL LSTIMER('AVERAG',TIMSTR,TIMEND,config%LUPRI)

      WRITE(config%LUPRI,'("** Get new density ")')
      call mat_no_of_matmuls(matmul1)
      CALL DOPT_get_density(config, fifoqueue, queue, F, H1, D, iteration)
      call mat_no_of_matmuls(matmul2)
      WRITE(config%LUPRI,'("No. of matmuls in get_density: ",I5)') matmul2-matmul1
      CALL LSTIMER('G_DENS',TIMSTR,TIMEND,config%LUPRI)

      CALL LSTIMER('SCF iteration',t1,t2,config%lupri)
   END DO
   CALL LSTIMER('**ITER',TSTR,TEN,config%LUPRI)
   print *, "loop done"
   config%solver%set_do_2nd_order = .false.
   call mat_no_of_matmuls(matmultot)
   write(config%lupri,*)
   WRITE(config%LUPRI,'("Total no. of matmuls in SCF optimization: ",I10)') matmultot
   call scf_afterplay(config,H1,S,D,E,gradnrm,F)
   if (config%solver%cfg_2nd_order_all) then
      config%solver%set_do_2nd_order = config%solver%cfg_do_2nd_order
   endif
!
! FREE SCF OPTIMIZATION STUFF (USEUALLY DONE AFTER SCF_AFTERPLAY)
!
  call DOPT_get_density_SHUTDOWN(config%diag)
  if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
     call modfifo_free(fifoqueue)
  !endif
  else
     call queue_free(config%av,queue)
  endif
  CALL mat_free(grad)
  call scf_stats_shutdown
!!!!Moved outside scfloop, needed also for lcv localization
!#if 0
!   if (cfg_density_method == cfg_f2d_direct_dens .or. cfg_density_method == cfg_f2d_arh .or. &
!       & cfg_check_converged_solution .or. cfg_rsp_nexcit > 0) then   
!     call decomp_shutdown
!     call dd_shutdown
!  endif
!#endif

END SUBROUTINE scfloop

!> \brief Wrapper for Fock/KS matrix.
!> \author L. Thogersen
!> \date 2002
!>
!> Get the new fock-matrix F(D). If we already evaluated one in densopt because
!> of a configuration-shift test, this one is used.
!>
subroutine get_fock(config,fifoqueue,queue,iteration,D,H1,F,Etotal,LSint,newlupri,newluerr,ls)
   use configuration
   use av_utilities
   use dal_interface
   use fock_evaluator,only: fck_get_fock_LSINT
   use ks_settings
   use TYPEDEF

   IMPLICIT NONE
   !> Contains all info about configuration/settings for SCF calculation
   type(configItem),intent(inout)         :: config
   !> New queue type: Contains Fock/KS and density matrices from previous SCF iterations (if ARH)
   type(modFIFO), intent(in)              :: fifoqueue
   !> Old queue type: Contains Fock/KS and density matrices from previous SCF iterations (if DIIS or DSM)
   type(util_HistoryStore), intent(inout) :: queue
   !> Current SCF iteration
   integer, intent(in)                    :: iteration
   !> Current density matrix
   TYPE(Matrix),intent(inout),target      :: D
   !> One-electron Hamiltonian
   TYPE(Matrix),intent(inout)             :: H1
   !> Fock/KS matrix
   TYPE(Matrix),intent(inout)             :: F
   !> SCF energy corresponding to constructed Fock/KS matrix 
   real(realk), INTENT(OUT)               :: Etotal
   !> True if both new and old integral code should be used (for debugging) (?)
   LOGICAL,intent(in)                     :: LSint
   !> Alternative output file if both new and old integral code should be used
   integer,intent(in)                     :: newlupri
   !> Alternative error file if both new and old integral code should be used
   integer,intent(in)                     :: newluerr
   !> Contains settings for integral code
   type(lsitem),intent(inout)             :: ls
   type(matrix)                :: wrk
   integer :: queue_lu, ndim
   logical, external :: do_dft
   integer :: previous

   ndim = F%nrow
     if (config%diag%cfg_no_confs_checked_in_rh) then 
!       Incremental scheme - build KS-matrix based on density-difference rather
!       than the density
        if (config%opt%cfg_queue_on_disk) then
           !Possibility to dump queue to disk while contruction Fock matrix for saving memory
           if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
              call fifoqueue_on_disk(fifoqueue,queue_lu,ndim)
           else
              call queue_on_disk(queue,queue_lu,ndim)
           endif
        endif
        CALL fck_get_fock_LSINT(D,H1,F,Etotal,LSint,newlupri,newluerr,ls)
        if (config%opt%cfg_queue_on_disk) then
           !Restore queue if it has been dumped to disk
           if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
              call fifoqueue_from_disk(fifoqueue,queue_lu,ndim)
           else
              call queue_from_disk(queue,queue_lu,ndim)
           endif
        endif
        if ((config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
            config%solver%cfg_do_2nd_order) .and. iteration > 1) then
           call update_trustradius(config%solver, ls, iteration, Etotal, fifoqueue%offset)
           if (.not. config%solver%step_accepted) Etotal = config%solver%old_energy
        endif
     else
        WRITE(config%lupri,'("** Fock matrix was already found in RH-step '// &
             & 'exploring a configuration shift ")')
        F   = queue%F(queue%current_position)
        D   = queue%D(queue%current_position)
        Etotal = queue%Energy(queue%current_position)
      endif
      !WRITE(LUPRI,*) 'E_SCF right after evaluation: ',Etotal
   !endif
end subroutine get_fock

!> @file
!> Contains main SCF driver and some wrappers to starting guess, Fock/KS matrix etc.



!> \brief After SCF opt., calculate HOMO-LUMO gap, lowest Hes. eigenvalue, print statistic 'n'stuff.
!> \author L. Thogersen
!> \date 2003
subroutine scf_afterplay(config,H1,S,D,E,gradnrm,F)
   use arh_debugging
   use configuration
   use direct_dens_util_unres
   !use debug_dchange_module
   USE Fock_evaluator, only: fck_unscale_virt
   use scf_stats
   implicit none
   !> Contains all info about configuration/settings for SCF calculation
   type(configItem), intent(inout) :: config
   !> One-electron Hamiltonian
   type(matrix), intent(in) :: H1
   !> Overlap matrix
   type(matrix), intent(in) :: S
   !> Converged AO density matrix
   type(matrix), intent(in) :: D
   !> Converged SCF energy
   real(realk), intent(in) :: E
   !> Final SCF gradient norm
   real(realk), intent(in) :: gradnrm
   !> Converged Fock/KS matrix. Output only if config%opt%cfg_scale_virt = .true.
   type(Matrix), intent(inout) :: F
   type(debugItem) :: debug
   type(DDitem)    :: DD

   !write (lupri,*) 'Incoming D, scf_afterplay:'
   !call MAT_PRINT(D, 1, D%nrow, 1, D%ncol, LUPRI)
   !write (lupri,*) 'Incoming F, scf_afterplay:'
   !call MAT_PRINT(F, 1, D%nrow, 1, D%ncol, LUPRI)
   if (config%decomp%cfg_rsp_nexcit > 0) then
      write(config%lupri,*) 'Postponing calculation of HOMO-LUMO gap to response part...'
   else if (config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
          & config%decomp%cfg_check_converged_solution .or. &
          & config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
      !debug_arh_hessian = .false.
      call get_oao_transformed_matrices(config%decomp,F,D)
      call dd_init(config%decomp,DD)
      if (config%diag%cfg_unres) then
         call dd_homolumo_and_heseigen_unres(DD,config%decomp,debug,.false.,0)
      else
         call dd_homolumo_and_heseigen(DD,config%decomp,debug,.false.,0)
      endif
      if (config%opt%cfg_diaghesonly .or. config%opt%debug_diag_hessian) then 
         call debug_diag_full_hessian(config%solver,config%decomp)
      endif
      call dd_shutdown(config%decomp,DD)
   endif
   if (config%opt%cfg_scale_virt) then
     !Fock matrix is modified - unmodify it
     config%opt%cfg_scale_virt = .false.
     call fck_unscale_virt(H1,S,D,F,config%decomp%nocc)
   endif
   CALL scf_stats_end_print(config%opt)

end subroutine scf_afterplay

!> \brief Update of trust-radius - see Molecular Electronic Structure Theory p. 615
!> \author S. Host
!> \date 2005
!>
!> To be used for arh and 2nd order optimization only.
!>
   subroutine update_trustradius(arh, ls, SCF_it, E, ndens)
   use arhDensity
   use decompMod
!   use fock_evaluator
   use TYPEDEF
   implicit none
      !> Contains setting for ARH solver
      type(SolverItem)               :: arh
      !> Contains settings for integrals
      type(lsitem),intent(inout)     :: ls
      !> Current SCF iteration
      integer, intent(in)            :: SCF_it 
      !> Current SCF energy
      real(realk)                    :: E
      !> Number of densities/Fock/KS matrices in subspace
      integer, intent(in)            :: ndens
      real(realk)                    :: r, tr_ratio
      character*8                    :: tr_criterion

   arh%debug_arh_scfit = SCF_it
   if (arh%set_optxelm) then
      tr_criterion = ' (xmax) '
   else
      tr_criterion = ' (xnorm)' 
   endif
   
   arh%step_accepted = .true.

   !The predicted SCF energy change is
   !Q = E(SCF)_n + (E1_red)T * c + 1/2*cT * (E2_red) * c
   !and the ratio between actual and predicted energy change is
   !    E(SCF)_n+1 - E(SCF)_n
   !r = _____________________
   !        Q - E(SCF)_n
 
   !The ratio denominator
   !Q - E(SCF)_n = (E1_red)T * c + 1/2*cT * (E2_red) * c
   !has been calculated in ARH module and is found in arh%denom  

   !write(arh%lupri,*) 'Denominator is:', arh%denom
   !write(arh%lupri,*) 'Numerator is:',  E - arh%old_energy
   !Now calculate the ratio:
   r = (E - arh%old_energy)/arh%denom
   if (r < 0) then
      arh%Nrejections = arh%Nrejections + 1
      if (arh%info_lineq) write(arh%lupri,*) 'Number of steps rejected:', arh%Nrejections
   endif
   !write(lupri,*) '% Old and new SCF energy:', ESCF_old, ESCF_new
   !write(lupri,*) '% SCF energy change:', ESCF_new - ESCF_old 
   !write(lupri,*) '% Predicted energy change:', denom 
   !write(lupri,*) '% Ratio between observed and predicted changes:', r 

  !if (.not. cfg_do_2nd_order) then
      if (arh%set_optxelm) then
      write (arh%lupri, "(F8.5, A8, F11.5, F11.5, F8.2, F10.4, F10.4, i10, i11, '    %%%')") &
           & arh%set_max_element, tr_criterion, arh%maxelm, arh%xnorm, -arh%current_mu, r, arh%D_para_D_tot, ndens, SCF_it-1
      else
      write (arh%lupri, "(F8.5, A8, F11.5, F11.5, F8.2, F10.4, F10.4, i10, i11, '    %%%')") &
           & arh%set_max_step, tr_criterion, arh%maxelm, arh%xnorm, -arh%current_mu, r, &
           & arh%D_para_D_tot, ndens, SCF_it-1
      endif
   !else
   !   if (arh%set_optxelm) then
   !   write (arh%lupri, "(F8.5, A8, F11.5, F11.5, F8.2, F10.4, F10.4, D14.4, D13.4, i5, '    %%%')") &
   !        & arh%set_max_element, tr_criterion, arh%maxelm, arh%xnorm, -arh%current_mu, r, &
   !        & arh%D_para_D_tot, arh%gradcontrib, arh%hescontrib, SCF_it-1
   !   else
   !   write (arh%lupri, "(F8.5, A8, F11.5, F11.5, F8.2, F10.4, F10.4, D14.4, D13.4, i5, '    %%%')") &
   !        & arh%set_max_step, tr_criterion, arh%maxelm, arh%xnorm, -arh%current_mu, r, &
   !        & arh%D_para_D_tot, arh%gradcontrib, arh%hescontrib, SCF_it-1
   !   endif
   !endif
 
!16/9 - 2010: Brano uses the solver for localization of orbital, and sometimes
!the ratio gets really large, which seems to cause problems. We therefore
!introduce:
! r > 0.95 double TR -> changed to 1.5d0 > r > 0.95
! 1.5 < r < 2.0 leave TR unchanged
! r > 2.0 contract TR
! /Stinne
!Removed, didn't work well! /Stinne 27/10-10

   if (arh%set_optxelm) then
      if (arh%info_lineq) write(arh%lupri,*) 'Trust radius based on max element!'
      tr_ratio = arh%maxelm/arh%set_max_element
      if (arh%maxelm > 5.0d-2 .and. arh%set_max_element > 1.0d-1 .and. &
          & tr_ratio < 0.9d0 .and. r > 0.0d0) then ! .and. r < 1.0d0 .and. abs(arh%current_mu) > 1.0d-2) then
         if (arh%info_lineq) write(arh%lupri,*) '% Too large difference between trust-radius and actual step. &
                        & Reduce trust-radius'
         if (arh%set_max_element*arh%cfg_arh_contract > arh%maxelm) then
            if (arh%info_lineq) write(arh%lupri,*) '% h(old), h(new):', &
                & arh%set_max_element, arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_element = arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
         else
            if (arh%info_lineq) write(arh%lupri,*) '% h(old), h(new):', arh%set_max_element, arh%maxelm
            arh%set_max_element = arh%maxelm
            arh%set_max_step = arh%xnorm 
         endif
      else
      !Now take action - change the trust-radius in accordance with r:
         if (r > arh%cfg_arh_expand_crit .and. r < 1.5d0) then  !we can take larger steps - expansion
            !arh%set_max_step = xnorm
            if (arh%trustradius_decreased) then
               if (arh%info_lineq) write(arh%lupri,*) "Trust-radius was decreased to obtain convergence, don't expand xmax!"
            !else if (abs(arh%current_mu) < 0.001 .and. arh%set_max_element >= cfg_trmaxelm) then
            !   write(lupri,*) "Mu ~ 0, trust region not expanded!"
            else if (r > 0.95d0) then ! .and. r < 1.5d0) then
               if (arh%info_lineq) write(arh%lupri,*) 'Large ratio, double trust-radius'
               arh%set_max_element = arh%set_max_element*2.0d0
               arh%set_max_step = arh%set_max_step*2.0d0
            !else if (tr_ratio < 0.9d0) then
            !   write(lupri,*) '% Too large difference between trust-radius and actual step - &
            !            & do not expand trust-radius.'
            else
               if (arh%info_lineq) write(arh%lupri,*) '% Expand trust-radius by a factor',  arh%cfg_arh_expand, 'h(old), h(new):', &
                              & arh%set_max_element, arh%set_max_element*arh%cfg_arh_expand
               arh%set_max_element = arh%set_max_element*arh%cfg_arh_expand
               arh%set_max_step = arh%xnorm*arh%cfg_arh_expand
               if (arh%info_lineq) write(arh%lupri,*) '% Max ||X|| updated to', arh%set_max_step
               if (arh%info_lineq) write(arh%lupri,*) '% Max ||X|| updated, old, new:', arh%set_max_step, arh%xnorm !Keep it safe...
            endif
         else if (arh%cfg_arh_contract_crit < r .and. r < arh%cfg_arh_expand_crit) then ! .or. &
                !& (r > 1.5d0 .and. r < 2.0d0)) then
            if (arh%info_lineq) write(arh%lupri,*) '% trust-radius ok, h = ', arh%set_max_element
            arh%set_max_step = arh%xnorm
            !if (r > 1.5d0 .and. r < 2.0d0 .and. arh%info_lineq) write(arh%lupri,*) '% 1.5 < r < 2.0, trust-radius unchanged'  
         else if (0.0d0 < r .and. r < arh%cfg_arh_contract_crit) then
         !else if ((0.0d0 < r .and. r < arh%cfg_arh_contract_crit) .or. &
         !        & r > 2.0d0) then
            if (arh%info_lineq) write(arh%lupri,*) '% Reduce trust-radius, h(old), h(new):', arh%set_max_element, &
                 &arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_element = arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_step = arh%xnorm
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
            !if (r > 2.0d0 .and. arh%info_lineq) write(arh%lupri,*) '% r > 2.0, trust-radius contracted'
         else if (r < 0.0d0) then
            arh%set_max_element = arh%maxelm
            arh%set_max_step = arh%xnorm
            if (arh%info_lineq) WRITE(arh%lupri,*) &
                 &'Reject step and reduce trust-radius, h(old), h(new):',&
                 & arh%set_max_element, arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_element = arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
            arh%step_accepted = .false.
         endif
      endif
   else
      if (arh%info_lineq) write(arh%lupri,*) 'Trust radius based on max norm!'
      tr_ratio = arh%xnorm/arh%set_max_step
      if (arh%xnorm > 5.0d-2 .and. arh%set_max_step > 1.0d-1 .and. &
       & tr_ratio < 0.9d0 .and. r > 0 .and. abs(arh%current_mu) > 1.0d-2) then
         if (arh%info_lineq) write(arh%lupri,*) '% Too large difference between trust-radius and actual step. Reduce trust-radius'
         if (arh%set_max_step*arh%cfg_arh_contract > arh%xnorm) then
            if (arh%info_lineq) write(arh%lupri,*) '% h(old), h(new):', arh%set_max_step, arh%set_max_step*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
         else
            if (arh%info_lineq) write(arh%lupri,*) '% h(old), h(new):', arh%set_max_step, arh%xnorm
            arh%set_max_step = arh%xnorm
         endif
      else
         if (r > arh%cfg_arh_expand_crit) then ! .and. r < 1.5d0) then  !we can take larger steps - expansion
            if (arh%trustradius_decreased) then
               if (arh%info_lineq) write(arh%lupri,*) "Trust-radius was decreased to obtain convergence, don't expand xnorm!"
            !else if (abs(arh%current_mu) < 0.001 .and. arh%set_max_step >= cfg_trmaxstep) then
            !   write(lupri,*) "Mu ~ 0, trust region not expanded!"
            else if (r > 0.95d0) then ! .and. r < 1.5d0) then
               if (arh%info_lineq) write(arh%lupri,*) 'Large ratio, double trust-radius'
               arh%set_max_step = arh%set_max_step*2.0d0
            else
               if (arh%info_lineq) write(arh%lupri,*) '% Expand trust-radius by a factor',  arh%cfg_arh_expand, 'h(old), h(new):', &
                              & arh%set_max_step, arh%set_max_step*arh%cfg_arh_expand
               arh%set_max_step = arh%set_max_step*arh%cfg_arh_expand
            endif
         else if (arh%cfg_arh_contract_crit < r .and. r < arh%cfg_arh_expand_crit) then !  .or. &
               ! & (r > 1.5d0 .and. r < 2.0d0)) then
            if (arh%info_lineq) then           
               !if (r > 1.5d0 .and. r < 2.0d0) then
               !   write(arh%lupri,*) '% 1.5 < r < 2.0, trust-radius unchanged'
               !else
                  write(arh%lupri,*) '% trust-radius ok, h = ', arh%set_max_step
               !endif
            endif
         else if (0.0d0 < r .and. r < arh%cfg_arh_contract_crit) then
         !else if (0.0d0 < r .and. r < arh%cfg_arh_contract_crit .or. &
         !        & r > 2.0d0) then
            if (arh%info_lineq) write(arh%lupri,*) '% Reduce trust-radius, h(old), h(new):', arh%set_max_step, &
                 &arh%set_max_step*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
            !if (r > 2.0d0 .and. arh%info_lineq) write(arh%lupri,*) '% r > 2.0, trust-radius contracted'
         else if (r < 0.0d0) then
            arh%set_max_step = arh%xnorm
            if (arh%info_lineq) write(arh%lupri,*) '% Reject step and reduce trust-radius, h(old), h(new):', &
                           & arh%set_max_step, arh%set_max_step*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
            arh%step_accepted = .false.
         endif
      endif
   endif

   if (arh%set_max_element > arh%cfg_max_element) then
      if (arh%info_lineq) write (arh%lupri,*) 'Absolute max element allowed is', arh%cfg_max_element,', resetting trust radius!'
      arh%set_max_element = arh%cfg_max_element
   endif
   if (arh%set_max_step > arh%cfg_max_step) then
      if (arh%info_lineq) write (arh%lupri,*) 'Absolute max step allowed is', arh%cfg_max_step,' , resetting trust radius!'
      arh%set_max_step = arh%cfg_max_step
   endif

   !Removed 13/10/2010, it doesn't work properly. /Stinne
   !if (.not. arh%step_accepted .and. arh%Nrejections > 1 .and. arh%set_local) then
   !   !If there is more than one rejection in the semilocal/local area,
   !   !the problem is often related to integral accuracy. Set a keyword
   !   !to get higher integral accuracy: 
   !   write(arh%lupri,*) 'Warning: Too many rejections. Integral accuracy will be increased'
   !   ls%setting%scheme%threshold = ls%setting%scheme%threshold*1.0D-1
   !   !This keyword is reset to false outside arh_solver
   !endif

   end subroutine update_trustradius

