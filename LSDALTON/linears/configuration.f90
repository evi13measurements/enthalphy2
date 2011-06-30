!> @file
!> Contains configuration module. 

!> \brief Type definitions for configuration structures and routines for reading input.
!> \author S. Host and T. Kjaergaard
!> \date March 2010
module configuration
use arhDensity
use av_utilities
use diagonalization
use decompMod
use files
use typedef
use opttype
use response_wrapper_type_module
use lsdalton_response_type_mod
use optimization_input
!> \brief Contains info, settings and data for response part (defaults or read from input file).
!> \author T. Kjaergaard
!> \date April 2010
type responseitem
   !> Used to store info about MCD calculation
   type(MCDinputitem)  :: MCDinput
   !> Used to store info about polarizability (alpha) calculation
   type(ALPHAinputitem)  :: ALPHAinput
   !> Used to store info about 1st hyperpolarizability (beta) calculation
   type(BETAinputitem)  :: BETAinput
   !> Used to store info about 2nd hyperpolarizability (gamma) calculation
   type(GAMMAinputitem)  :: GAMMAinput
   !> Used to store info about standard TPA calculation.
   type(TPAinputitem)  :: TPAinput
   !> Used to store info about damped TPA calculation.
   type(DTPAinputitem)  :: DTPAinput
   !> Used to store info about excited state gradient calculation.
   type(ESGinputitem)  :: ESGinput
   !> Used to store info about excited state dipole moment calculation.
   type(ESDinputitem)  :: ESDinput
   !> Used to store info about solver that is used for calculation.
   type(RSPSOLVERinputitem) :: RSPSOLVERinput
   type(rsp_tasksitem) :: tasks
end type responseitem

!> \brief Contains info, settings and data for entire calculation (defaults or read from input file).
!> \author S. Host
!> \date March 2010
type ConfigItem
   !> Logical unit number for LSDALTON.OUT
   integer              :: lupri
   !> Used for Augmented Roothaan-Hall, direct density optimization etc.
   type(SolverItem)     :: solver
   !> Used to store OAO decomposition of overlap matrix + OAO transformed matrices
   type(DecompItem)     :: decomp
   !> Used to store info about integrals
   type(integralconfig) :: integral
   !> Used to store info about density optimization type
   type(OptItem)        :: opt
   !> Used to store info about SCF averaging
   type(AvItem)         :: av
   !> Used to store info about diagonalization
   type(DiagItem)       :: diag
   !> Used to store info about molecule
   type(moleculeinfo)   :: molecule
   !> Used to store info about which atoms have which basisset
   type(basissetlibraryitem):: lib
   !> Used to store info about response calculation
   type(responseitem) :: response
   !> Used to store info about optimization
   type(opt_setting)  :: optinfo
   !Only for testing new sparse matrix library, should be removed afterwards!
   logical            :: sparsetest
end type ConfigItem

contains

!> \brief Call routines to set default values for different structures.
!> \author S. Host
!> \date March 2010
subroutine config_set_default_config(config)
implicit none
   !> Contains info, settings and data for entire calculation
   type(ConfigItem), intent(inout) :: config

  call arh_set_default_config(config%solver)
  call decomp_set_default_config(config%decomp)
  call integral_set_default_config(config%integral)
  call opt_set_default_config(config%opt)
  call av_set_default_config(config%av)
  call diag_set_default_config(config%diag)
  !RESPONSE
  ! Polarizability
  call ALPHAinputitem_set_default_config(config%response%alphainput)
  ! 1st hyperpolarizability
  call BETAinputitem_set_default_config(config%response%betainput)
  ! 2nd hyperpolarizability
  call GAMMAinputitem_set_default_config(config%response%gammainput)
  ! Standard TPA
  call TPAinputitem_set_default_config(config%response%tpainput)
  ! Damped TPA
  call DTPAinputitem_set_default_config(config%response%dtpainput)
  ! ESG
  call ESGinputitem_set_default_config(config%response%esginput)
  ! ESD
  call ESDinputitem_set_default_config(config%response%esdinput)
  ! MCD
  call MCDinputitem_set_default_config(config%response%MCDinput)
  ! RSP solver
  call RSPSOLVERiputitem_set_default_config(config%response%RSPSOLVERinput)
  call rsp_tasks_set_default_config(config%response%tasks)
  ! geometry optimization
  call optimization_set_default_config(config%optinfo)
  !Only for testing new sparse matrix library, should be removed afterwards!
  config%sparsetest = .false.
end subroutine config_set_default_config

!> \brief Wrapper to routines for read input files DALTON.INP and MOLECULE.INP.
!> \author S. Host and T. Kjaergaard
!> \date March 2010
subroutine config_read_input(config,lupri,luerr,nonlsdaltonrun)
  use readmolefile
implicit none
   !> Contains info, settings and data for entire calculation
   type(ConfigItem), intent(inout) :: config
   !> Logical unit number for LSDALTON.OUT
   integer, intent(in)             :: lupri
   !> Logical unit number for LSDALTON.ERR
   integer, intent(in)             :: luerr
   LOGICAL :: LSint,nonlsdaltonrun

   config%lupri          = lupri
   config%solver%lupri   = lupri
   config%decomp%lupri   = lupri
   config%decomp%luerr   = luerr
   config%av%lupri       = lupri
   config%opt%lupri      = lupri
   config%opt%luerr      = luerr
   config%diag%lupri     = lupri

   !read the MOLECULE.INP and set input
   call read_mol_file_and_build_molecule(lupri,config%molecule,config%LIB,&
        &.FALSE.,0,config%integral%DoSpherical,&
        &config%integral%IntegralThreshold,config%integral%Auxbasis)
   config%integral%nelectrons = config%molecule%nelectrons 
   config%integral%molcharge = INT(config%molecule%charge)
   !read the DALTON.INP and set input
   CALL read_dalton_input(LUPRI,config,nonlsdaltonrun) 

end subroutine config_read_input

!> \brief If data has been allocated by the default setting, it should be free'd here.
!> \author S. Host
!> \date March 2010
subroutine config_shutdown(config)
implicit none
   !> Contains info, settings and data for entire calculation
   type(ConfigItem), intent(inout) :: config

   call av_shutdown(config%av)
   call free_Moleculeinfo(config%Molecule)
   call free_MCDinputitem(config%response%MCDinput)

end subroutine config_shutdown

!> \brief Read input file DALTON.INP and set configuration structure accordingly.
!> \author T. Kjaergaard
!> \date March 2010
SUBROUTINE read_dalton_input(LUPRI,config,nonlsdaltonrun)
! READ THE INPUT FOR THE INTEGRAL 
implicit none
!> Logical unit number for LSDALTON.OUT
INTEGER            :: LUPRI
!> Contains info, settings and data for entire calculation
type(ConfigItem), intent(inout) :: config
INTEGER            :: LUCMD !Logical unit number for the daltoninput
INTEGER            :: IDUMMY,IPOS,IPOS2,COUNTER
character(len=70)  :: WORD
character(len=2)   :: PROMPT
LOGICAL            :: DONE,file_exists,READWORD,LSDALTON,STARTGUESS,nonlsdaltonrun
!LINSCA variables:
real(realk)        :: shift, min_density_overlap, maxratio, zero
integer            :: nvec, i,inperr
Real(realk)  :: hfweight 

STARTGUESS = .FALSE.
Config%integral%cfg_lsdalton = .TRUE.
COUNTER = 0

INQUIRE(file='DALTON.INP',EXIST=file_exists) 
IF(file_exists)THEN
   LUCMD=-1
   CALL lsOPEN(LUCMD,'DALTON.INP','OLD','FORMATTED')
ELSE
   CALL lsQUIT('DALTON.INP does not exist',lupri)
ENDIF
READWORD=.TRUE.
DONE=.FALSE.
rewind(LUCMD)
DO
   IF(DONE)EXIT
   IF(READWORD) THEN
      READ (LUCMD, '(A40)') WORD
      READWORD=.TRUE.
      COUNTER = 0
   ELSE
      IF (COUNTER.GT.1) THEN
        WRITE(LUPRI,'(1X,2A)') 'Infinite loop for input line:',WORD
        CALL lsQUIT('Infinite loop in read_dalton_input,due to wrong input',lupri)
      ENDIF
      COUNTER = COUNTER + 1
   ENDIF
!
Write(*,*)'WORD=',WORD
   PROMPT = WORD(1:2)
   IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#'))CYCLE
!   IF (WORD(1:14) == '**DALTON INPUT') CYCLE
   IF (WORD(1:10) == '**INTEGRAL') THEN
      READWORD = .TRUE.
      CALL INTEGRAL_INPUT(config%integral,readword,word,lucmd,lupri)
   ENDIF
   IF ((WORD(1:10) == '**WAVE FUN').OR.(WORD(1:10) == '**WAVEFUNC')) THEN
      READWORD=.TRUE.
      DO   
         IF(READWORD) THEN
            READ (LUCMD, '(A40)') WORD
            READWORD=.TRUE.
         ENDIF
         PROMPT = WORD(1:2)
         IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
         IF(PROMPT(1:1) .EQ. '.') THEN
            SELECT CASE(WORD(1:7)) 
            CASE ('.HF');  config%opt%calctype = config%opt%hfcalc !Hartree-Fock calc
                     config%integral%exchangeFactor = 1.d0

            CASE ('.DFT'); config%opt%calctype = config%opt%dftcalc !DFT calc
                           config%av%CFG_SET_type = config%av%CFG_THR_dft
                           config%solver%do_dft = .true.
               DO 
                  READ (LUCMD, '(A40)') WORD
                  IF ((WORD(1:1) .EQ. '!') .OR. (WORD(1:1) .EQ. '#')) CYCLE   
                  IF (WORD(1:1) .EQ. '.' .OR. WORD(1:1) .EQ. '*') THEN
                     WRITE (LUPRI,'(/A/A//A)')&
                          & '--> Input error for line following .DFT',&
                          & '    expected functional specification but read:',&
                          & WORD
                  ELSE
!                     IF(WORD(1:3) .EQ. 'LDA') 
                     IPOS = INDEX(WORD,'CAM')
                     IPOS2 = INDEX(WORD,'cam')                     
                     IF(IPOS .NE. 0 .OR. IPOS2 .NE. 0)THEN !CAM
                        config%integral%CAM=.TRUE.
                        IPOS = INDEX(WORD,'alpha')
                        IF (IPOS .NE. 0) THEN
                           IPOS2 = INDEX(WORD(IPOS:),'=')
                           IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
                              WRITE (LUPRI,'(2X,A40)') 'Incorrect input for CAM parameters'
                              WRITE (LUPRI,'(2X,A40)') 'Format is "alpha=?  beta=? mu=?"'
                              CALL lsQUIT('Incorrect input for alpha parameter',lupri)
                           ELSE
                              READ (WORD((IPOS+IPOS2):70),*) config%integral%CAMalpha
                              IPOS = INDEX(WORD,'beta')
                              IF (IPOS .NE. 0) THEN
                                 IPOS2 = INDEX(WORD(IPOS:),'=')
                                 IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 5)) THEN
                                    WRITE (LUPRI,'(2X,A40)') 'Incorrect input for CAM parameters'
                                    WRITE (LUPRI,'(2X,A40)') 'Format is "alpha=?  beta=? mu=?"'
                                    CALL lsQUIT('Incorrect input for alpha parameter',lupri)
                                 ELSE
                                    READ (WORD((IPOS+IPOS2):70),*) config%integral%CAMbeta
                                    IPOS = INDEX(WORD,'mu')
                                    IF (IPOS .NE. 0) THEN
                                       IPOS2 = INDEX(WORD(IPOS:),'=')
                                       IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 3)) THEN
                                          WRITE (LUPRI,'(2X,A40)') 'Incorrect input for CAM parameters'
                                          WRITE (LUPRI,'(2X,A40)') 'Format is "alpha=?  beta=? mu=?"'
                                          CALL lsQUIT('Incorrect input for alpha parameter',lupri)
                                       ELSE
                                          READ (WORD((IPOS+IPOS2):70),*) config%integral%CAMmu
                                       ENDIF
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDIF
                        ELSE
                           config%integral%CAMalpha=0.19D0
                           config%integral%CAMbeta=0.46D0
                           config%integral%CAMmu=0.33D0
                        ENDIF
                        WRITE(LUPRI,*) 'This is a CAM functional with'
                        WRITE(LUPRI,*) 'config%integral%CAMalpha',config%integral%CAMalpha
                        WRITE(LUPRI,*) 'config%integral%CAMbeta',config%integral%CAMbeta
                        WRITE(LUPRI,*) 'config%integral%CAMmu',config%integral%CAMmu
                     END IF
                     INPERR = 1
                     hfweight=0.D0 
                     !it is assumed that hfweight is set to zero and only  
                     !changed if the functional require a HF weight  
                     !different from zero. 
                     print*,'WORD(1:40)',WORD(1:40)
                     CALL DFTsetFunc(WORD(1:40),INPERR,40,hfweight)
                     config%integral%exchangeFactor = hfweight
                  ENDIF
                  EXIT
               ENDDO
            CASE DEFAULT
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                    & '" not recognized in **WAVE FUNCTION'
               CALL lsQUIT('Illegal keyword in **WAVE FUNCTION',lupri)
            END SELECT
         ELSE
            WRITE(LUPRI,'(2X,A)')'Requires .DFT or .HF after **WAVEFUNCTION when .LINSCA specified'
            CALL lsQUIT('Requires .DFT or .HF after **WAVEFUNCTION when .LINSCA specified',lupri)

         ENDIF
         EXIT
      ENDDO
   ENDIF
   IF (WORD(1:10) == '*DFT INPUT') THEN
      CALL READ_WAVE_DFTINPUT(LUPRI,LUCMD,config%integral,WORD)
      READWORD=.FALSE.
!      EXIT !done with what we are interested in
   ENDIF
   IF (WORD(1:7) == '*LINSCA') THEN
      write(lupri,*) '*LINSCA under **WAVE FUNCTION has been replaced by *DENSOPT'
      call lsquit('*LINSCA under **WAVE FUNCTION has been replaced by *DENSOPT',-1)
   else if (WORD(1:8) == '*DENSOPT') THEN
      READWORD=.TRUE.
      DO   
         IF(READWORD) THEN
            READ (LUCMD, '(A40)') WORD
            READWORD=.TRUE.
         ENDIF
         PROMPT = WORD(1:2)
         IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
         IF(PROMPT(1:1) .EQ. '.') THEN
            SELECT CASE(WORD) 
            CASE('.2ND_ALL');    config%solver%cfg_2nd_order_all = .true.
                                 config%solver%cfg_do_2nd_order = .true.
                                 config%solver%set_do_2nd_order = .true.
                                 config%opt%CFG_density_method =  config%opt%CFG_F2D_DIRECT_DENS
                                 config%solver%cfg_arh_truncate = .false.
            CASE('.2ND_LOC');    config%solver%cfg_2nd_order_local = .true.
                                 config%solver%cfg_do_2nd_order = .true.
                                 config%opt%CFG_density_method =  config%opt%CFG_F2D_DIRECT_DENS
            CASE('.ARH');        config%solver%cfg_arh_truncate = .true.
                                 config%solver%cfg_arh_crop = .true.
                                 config%opt%CFG_density_method = config%opt%CFG_F2D_ARH
            CASE('.ARH FULL');   config%solver%cfg_arh_crop = .true.
                                 config%solver%cfg_arh_truncate = .false.
                                 config%opt%CFG_density_method = config%opt%CFG_F2D_ARH
            CASE('.ASYM');       config%opt%cfg_asym = .true.
            CASE('.BLOCK');      CALL mat_select_type(mtype_sparse_block)
                                 config%opt%cfg_prefer_BSM = .true.
            CASE('.CHOLESKY');   config%decomp%lowdin_diagonalize = .false.; config%decomp%cholesky_decomp   = .true.
            CASE('.CONFSHIFT');  config%diag%cfg_no_conf_shift = .false.
            CASE('.CONTFAC');    READ(LUCMD,*) config%solver%cfg_arh_contract
            CASE('.CONTRAC');    READ(LUCMD,*) config%solver%cfg_arh_contract_crit
            CASE('.CONVDYN');    READ(LUCMD,*) config%opt%cfg_convdyn_type ; config%opt%cfg_convdyn = .true.
            CASE('.CONVTHR');    READ(LUCMD,*) config%opt%cfg_convergence_threshold
                                 config%opt%set_convergence_threshold = config%opt%cfg_convergence_threshold
            CASE('.CSR');        config%opt%cfg_prefer_CSR = .true.
            CASE('.DIAGHESONLY'); config%opt%cfg_diaghesonly = .true.
            CASE('.DIIS');       config%av%CFG_averaging = config%av%CFG_AVG_DIIS
            CASE('.DISK');       config%opt%cfg_queue_on_disk = .true.
            CASE('.DISKSOLVER'); config%solver%cfg_arh_disk_micro = .true.
            !CASE('.DISKQUEUE') ; config%solver%cfg_arh_disk_macro = .true. !Not active - get_from_modFIFO_disk won't work!
            CASE('.DORTH');      config%diag%CFG_lshift = config%diag%cfg_lshift_dorth
                                 config%av%CFG_lshift = config%av%cfg_lshift_dorth
            CASE('.DSM');        config%av%CFG_averaging = config%av%CFG_AVG_DSM
            CASE('.DSMONE');     config%av%cfg_averaging = config%av%cfg_avg_dsm
                                 config%av%cfg_dsm_app = config%av%cfg_dsm_one
            CASE('.DSMXTRA');    config%av%cfg_averaging = config%av%cfg_avg_dsm
                                 config%av%cfg_dsm_app = config%av%cfg_dsm_xtra_term
            CASE('.DUMPMAT');    config%opt%dumpmatrices = .true.
            CASE('.EDIIS');      config%av%CFG_averaging = config%av%CFG_AVG_EDIIS
            CASE('.EXPAND');     READ(LUCMD,*) config%solver%cfg_arh_expand_crit
            CASE('.EXPFAC');     READ(LUCMD,*) config%solver%cfg_arh_expand 
            CASE('.FIXSHIFT');   READ(LUCMD,*) shift 
                                 config%solver%cfg_fixed_shift_param = shift ; config%solver%cfg_fixed_shift = .true.
                                 config%diag%cfg_fixed_shift_param = shift   ; config%diag%cfg_fixed_shift = .true.
            CASE('.FLUSH');      config%av%cfg_flush_vec = .true.
            CASE('.GCBASIS');    config%integral%NOSEGMENT = .TRUE. ; config%integral%TRILEVEL = .TRUE.
                                 config%decomp%cfg_gcbasis = .true.
            CASE('.HESONLY');    config%opt%cfg_hesonly = .true.
                                 config%decomp%cfg_check_converged_solution = .true.
                                 config%decomp%cfg_hessian_nvec = 1
            CASE('.HESVEC');     READ(LUCMD,*) config%decomp%cfg_hessian_nvec 
                                 config%decomp%cfg_check_converged_solution = .true.
            CASE('.HLMAXIT');    READ(LUCMD,*) config%decomp%cfg_homolumo_maxit
            CASE('.L2THR');      READ(LUCMD,*) config%opt%cfg_level2_convfactor
            CASE('.LCV');        config%decomp%cfg_lcv = .true.
            CASE('.LCM');        config%decomp%cfg_lcv = .true. ; config%decomp%cfg_lcm=.true.
            CASE('.LEVELSH');    ALLOCATE(config%diag%cfg_levelshifts(100)) ; config%diag%cfg_levelshifts = 0.0d0
                                 READ(LUCMD,*) config%diag%cfg_nshifts,(config%diag%cfg_levelshifts(i),i=1,config%diag%cfg_nshifts)
                                 config%diag%cfg_fixed_shift = .true. ; config%diag%cfg_custom_shift = .true.
            CASE('.LINCOMB');    READ(LUCMD,*) config%opt%cfg_weight_param 
            CASE('.LWITER') ;    config%decomp%lowdin_diagonalize = .false.; config%decomp%lowdin_iterative  = .true. 
            CASE('.LWQITER') ;   config%decomp%lowdin_diagonalize = .false.; config%decomp%lowdin_qiterative = .true.
            CASE('.MAXELM');     READ(LUCMD,*) config%solver%cfg_max_element
                                               config%solver%set_max_element = config%solver%cfg_max_element
            CASE('.MAXIT');      READ(LUCMD,*) config%opt%cfg_max_linscf_iterations
            CASE('.MAXRATI');    READ(LUCMD,*) maxratio
                                               config%av%cfg_settings%max_dorth_ratio = maxratio
            CASE('.MAXSTEP');    READ(LUCMD,*) config%solver%cfg_max_step 
                                               config%solver%set_max_step = config%solver%cfg_max_step
            CASE('.MICTHRS');    READ(LUCMD,*) config%solver%cfg_micro_thresh
                                 config%decomp%cfg_micro_thresh = config%solver%cfg_micro_thresh
            CASE('.MICROVECS');  READ(LUCMD,*) config%solver%cfg_arh_microvecs
            CASE('.MINDAMP');    READ(LUCMD,*) config%solver%cfg_min_lshift
                                               config%diag%cfg_min_lshift = config%solver%cfg_min_lshift
            CASE('.MOCHANGE');   config%diag%cfg_lshift = config%diag%cfg_lshift_MOchange
                                 config%av%cfg_lshift = config%av%cfg_lshift_MOchange
            CASE('.MUOPT');      config%diag%CFG_lshift = config%diag%CFG_lshift_search
                                 config%av%CFG_lshift = config%av%CFG_lshift_search
            CASE('.NALPHA');     read(LUCMD,*) config%decomp%nocca ; config%decomp%alpha_specified = .true.
            CASE('.NBETA');      read(LUCMD,*) config%decomp%noccb ; config%decomp%beta_specified = .true.
            CASE('.NOAV');       config%av%CFG_averaging =   config%av%CFG_AVG_none  
            CASE('.NEWDAMP');    config%solver%cfg_arh_newdamp = .true.
            CASE('.NVEC');       READ(LUCMD,*) NVEC; config%av%cfg_settings%max_history_size = NVEC
                                 config%av%dsm_history_size = NVEC
                                 config%av%diis_history_size = NVEC
                                 config%av%ediis_history_size = NVEC
            CASE('.NVECDSM');    READ(LUCMD,*) config%av%dsm_history_size
            CASE('.NVECDII');    READ(LUCMD,*) NVEC
                                 config%av%diis_history_size = NVEC
                                 config%av%ediis_history_size = NVEC
            CASE('.NOPREC');     config%solver%cfg_NOPREC = .true.
                                 config%decomp%cfg_NOPREC = .true.
            CASE('.INCREMENT');  config%opt%cfg_incremental = .true.
            CASE('.NOSHIFT');    config%diag%cfg_lshift = config%diag%cfg_lshift_none
                                 config%av%cfg_lshift = config%av%cfg_lshift_none
                                 config%solver%cfg_nodamp = .true.
            CASE('.OAO');        config%decomp%cfg_do_in_oao = .true.; config%opt%cfg_start_guess = 'H1OAO'
            CASE('.OVERLAP');    READ(LUCMD,*) min_density_overlap
                                 config%av%cfg_settings%min_density_overlap &
                                     & = min_density_overlap
            !CASE('.PURIFY');     config%opt%cfg_density_method = config%opt%cfg_f2d_purification - NO LONGER SUPPORTED! /Stinne 16-08-2010
            !                     read(LUCMD,*) config%opt%cfg_purification_method
            CASE('.RESTART');    config%diag%CFG_restart =  .TRUE. 
            CASE('.RH');         config%opt%CFG_density_method =  config%opt%CFG_F2D_ROOTHAAN
            CASE('.SAFE');       config%av%CFG_safe = .true.
            CASE('.SCALVIR');    config%opt%cfg_scale_virt = .true.
            CASE('.SOEO');       config%opt%cfg_soeo = .true.
            CASE('.SPARSE');     CALL mat_select_type(mtype_sparse1)
            CASE('.SPARSETEST'); config%sparsetest = .true.
            CASE('.SPIN');       READ(LUCMD,*) config%decomp%spin
            CASE('.STABILITY');  config%decomp%cfg_check_converged_solution = .true.
                                 config%decomp%cfg_hessian_nvec = 1
            CASE('.STAB MAXIT'); READ(LUCMD,*) config%decomp%cfg_check_maxit
            CASE('.START');      READ(LUCMD,*) config%opt%cfg_start_guess 
                                 STARTGUESS = .TRUE.
                                 if (config%opt%cfg_start_guess == 'TRILEVEL' .or. &
                                 &   config%opt%cfg_start_guess == 'ATOMS') then
                                    !this is default 
                                    config%integral%nosegment = .TRUE.
                                    config%integral%NOFAMILY = .TRUE.
                                    config%integral%trilevel = .TRUE.
                                 else !h1diag
                                    config%integral%NOSEGMENT = .FALSE.
                                    config%integral%TRILEVEL = .FALSE.
                                    config%integral%NOFAMILY = .FAlSE.
                                 endif
                                 !READ (LUCMD, '(A40)') WORD
                                 !IPOS = INDEX(WORD,'TRILEVEL')
                                 !IF(IPOS.NE.0)THEN
                                 !   config%integral%nosegment = .TRUE.
                                 !   config%integral%trilevel = .TRUE.
                                 !ENDIF
                                 !IPOS = INDEX(WORD,'ATOMS')
                                 !IF(IPOS.NE.0)THEN
                                 !   config%integral%nosegment = .TRUE.
                                 !   config%integral%trilevel = .TRUE.
                                 !ENDIF
            CASE('.TRSCF');      config%opt%CFG_density_method = config%opt%CFG_F2D_ROOTHAAN
                                 config%diag%cfg_lshift = config%diag%cfg_lshift_dorth
                                 config%av%cfg_lshift = config%av%cfg_lshift_dorth
                                 config%av%CFG_averaging = config%av%CFG_AVG_DSM
            CASE('.TrFD');       config%opt%CFG_density_method =  config%opt%CFG_F2D_DIRECT_DENS
            CASE('.TrFD FULL');  config%opt%CFG_density_method =  config%opt%CFG_F2D_DIRECT_DENS
                                 config%solver%cfg_arh_truncate = .false.
            CASE('.UNREST');     config%decomp%cfg_unres=.true.
                                 config%integral%unres=.true.
                                 config%diag%cfg_unres=.true.
                                 config%opt%cfg_unres=.true.
                                 config%response%RSPsolverinput%cfg_unres = .true.
            CASE('.UNSAFE');     config%solver%cfg_arh_crop_safe = .false.
            CASE('.VanLenthe');  config%opt%CFG_density_method =  config%opt%CFG_F2D_ROOTHAAN !Diagonalization
                                 config%av%CFG_averaging = config%av%CFG_AVG_van_lenthe
                                 config%av%CFG_lshift = config%av%cfg_lshift_vanlenthe
                                 config%diag%CFG_lshift = config%diag%cfg_lshift_vanlenthe
            CASE('.ZERO');       READ(LUCMD,*) zero
                                 call mat_zero_cutoff(zero)
            CASE DEFAULT
               WRITE(config%LUPRI,*) ' Keyword ',WORD,' not recognized in read_dalton_input'
               !CALL lsQUIT('Illegal keyword in read_dalton_input.',config%lupri)
            END SELECT
         ELSE IF (PROMPT(1:1) .EQ. '$') THEN
            IF (WORD == '$INFO') THEN
              call config_info_input(config,lucmd)
              cycle
            ELSE
              WRITE(LUPRI,*) ' Keyword ',WORD,' not recognized in read_dalton_input'
              CALL lsQUIT('Illegal keyword in read_dalton_input.',config%lupri)
              cycle
            ENDIF
         ENDIF
         IF(PROMPT .EQ. '**') THEN
            DONE=.TRUE.
            EXIT
            READWORD=.FALSE.
         ENDIF
         IF(PROMPT(1:1) .EQ. '*') THEN
            EXIT
            READWORD=.FALSE.
         ENDIF
      ENDDO
   ENDIF


   ! KK, change from $RESPONS to **RESPONS to be consistent with other input structure.
   ResponseInput: IF (WORD(1:9) == '**RESPONS') THEN
      READWORD=.TRUE.
      call config_rsp_input(config,lucmd,readword)
   END IF ResponseInput
!   
! Find optimization input section
!
   IF (WORD(1:7) .EQ. '**OPTIM') THEN  
      config%optinfo%optimize = .TRUE.
      CALL LS_optimization_input(config%optinfo,readword,word,lucmd, &
      lupri,config%molecule%nAtoms)
   ENDIF
!
   IF (WORD == '*END OF INPUT') THEN
      DONE=.TRUE.
   ENDIF
ENDDO
!IF(LSDALTON .AND. (.NOT. STARTGUESS)) config%integral%TRILEVEL = .TRUE.
!IF(.NOT.config%integral%LINSCA)THEN
!   WRITE (LUPRI,'(/,3A,/)') ' The Keyword .LINSCA was not used in&
!   & **INTEGRAL section.'
!   CALL lsQUIT('.LINSCA WAS NOT USED.')
!ENDIF
CALL lsCLOSE(LUCMD,'KEEP')
END SUBROUTINE read_dalton_input

subroutine INTEGRAL_INPUT(integral,readword,word,lucmd,lupri)
  implicit none
  LOGICAL,intent(inout)                :: READWORD
  TYPE(integralconfig),intent(inout)   :: integral
  character(len=70)  :: WORD
  INTEGER,intent(in) :: LUCMD !Logical unit number for the daltoninput
  INTEGER,intent(in) :: LUPRI !Logical unit number for the daltonoutput file
!
  INTEGER            :: IDUMMY
  character(len=2)   :: PROMPT
  integral%LINSCA = .TRUE. !should be obsolete

  DO   
     IF(READWORD) THEN
        READ (LUCMD, '(A40)') WORD
        READWORD=.TRUE.
     ENDIF
     PROMPT = WORD(1:2)
     IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
     IF(PROMPT .EQ. '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(PROMPT(1:1) .EQ. '.') THEN
        SELECT CASE(WORD) 
        CASE ('.FRAGMENT'); READ(LUCMD,*) INTEGRAL%numAtomsPerFragment; INTEGRAL%FRAGMENT = .TRUE.
        CASE ('.2CENTERERI'); INTEGRAL%DO2CENTERERI = .TRUE.
        CASE ('.3CENTEROVL'); INTEGRAL%DO3CENTEROVL = .TRUE.
        CASE ('.4CENTERERI');  INTEGRAL%DO4CENTERERI = .TRUE.
        CASE ('.AOPRINT');  READ(LUCMD,*) INTEGRAL%AOPRINT
        CASE ('.TEST_NDMAT_COULOMB'); INTEGRAL%TEST_NDMAT_COULOMB = .TRUE.
        CASE ('.TEST_NDMAT_EXCHANGE'); INTEGRAL%TEST_NDMAT_EXCHANGE = .TRUE.
        CASE ('.BASPRINT');  READ(LUCMD,*) INTEGRAL%BASPRINT
        CASE ('.DEBUG PS');  INTEGRAL%PS_DEBUG = .TRUE.
        CASE ('.DEBUGNUCREP'); INTEGRAL%DEBUGNUCREP = .TRUE.
        CASE ('.DEBUGKINETIC'); INTEGRAL%DEBUGKINETIC = .TRUE.
        CASE ('.DEBUGOVERLAP');  INTEGRAL%DEBUGOVERLAP = .TRUE.
        CASE ('.DEBUGDECPACKED');  INTEGRAL%DEBUGDECPACKED = .TRUE.
        CASE ('.DEBUG4CENTER');  INTEGRAL%DEBUG4CENTER = .TRUE.
        CASE ('.DEBUG4CENTERERI');  INTEGRAL%DEBUG4CENTER_ERI = .TRUE.
        CASE ('.DEBUGCCFRAGMENT');  INTEGRAL%DEBUGCCFRAGMENT = .TRUE.
        CASE ('.CARMOM');  READ(LUCMD,*) INTEGRAL%CARMOM
        CASE ('.SPHMOM');  READ(LUCMD,*) INTEGRAL%SPHMOM
        CASE ('.CART-E'); INTEGRAL%DoCartesian = .TRUE.
        CASE ('.INTPRINT');  READ(LUCMD,*) INTEGRAL%INTPRINT
        CASE ('.NOJENGINE'); INTEGRAL%JENGINE = .FALSE.
        CASE ('.MAXPASSES'); READ(LUCMD,*) INTEGRAL%maxpasses
        CASE ('.NOLINK'); INTEGRAL%LINK = .FALSE.
        CASE ('.DALINK'); INTEGRAL%LINK = .TRUE.
           INTEGRAL%DALINK = .TRUE.
        CASE ('.NSETUV'); INTEGRAL%nonSphericalETUV = .TRUE.
        CASE ('.ETUVIL'); INTEGRAL%ETUVinsideLOOP = .TRUE.
        CASE('.LOW RJ000 ACCURACY'); INTEGRAL%HIGH_RJ000_ACCURACY = .FALSE.
        CASE('.PrimAng'); INTEGRAL%OrderAngPrim = .FALSE.
        CASE('.FTUVMAXPRIM'); READ(LUCMD,*)INTEGRAL%FTUVmaxprim
        CASE ('.LINSCAPRINT')
           READ(LUCMD,*) INTEGRAL%LINSCAPRINT
           INTEGRAL%MOLPRINT=INTEGRAL%LINSCAPRINT
           INTEGRAL%BASPRINT=INTEGRAL%LINSCAPRINT
           INTEGRAL%AOPRINT=INTEGRAL%LINSCAPRINT
           INTEGRAL%INTPRINT=INTEGRAL%LINSCAPRINT
        CASE ('.PRINTATOMCOORD');  INTEGRAL%PRINTATOMCOORD = .TRUE.
        CASE ('.MOLPRINT');  READ(LUCMD,*) INTEGRAL%MOLPRINT
        CASE ('.NO SCREEN');  
           INTEGRAL%CS_SCREEN = .FALSE. 
           INTEGRAL%PS_SCREEN = .FALSE. 
           INTEGRAL%OD_SCREEN = .FALSE.
        CASE ('.THRESH')
           READ(LUCMD,*) INTEGRAL%THRESHOLD
        CASE ('.NO PASS');  INTEGRAL%DOPASS = .FALSE. 
        CASE ('.NO CS');  INTEGRAL%CS_SCREEN = .FALSE. 
        CASE ('.NO PS');  INTEGRAL%PS_SCREEN = .FALSE. 
        CASE ('.NOFAMILY'); INTEGRAL%NOFAMILY = .TRUE.
        CASE ('.OVERLAP-DF-J'); INTEGRAL%OVERLAP_DF_J=.TRUE.
        CASE ('.THR_CS');   READ(LUCMD,*) INTEGRAL%CS_THRESHOLD
           !SET IN RELATION TO THE ONE THRESHOLD
           INTEGRAL%CS_THRESHOLD = INTEGRAL%CS_THRESHOLD/INTEGRAL%THRESHOLD 
        CASE ('.THR_PS');   READ(LUCMD,*) INTEGRAL%PS_THRESHOLD
           !SET IN RELATION TO THE ONE THRESHOLD
           INTEGRAL%PS_THRESHOLD = INTEGRAL%PS_THRESHOLD/INTEGRAL%THRESHOLD 
        CASE ('.TIMINGS'); INTEGRAL%TIMINGS = .TRUE.
        CASE ('.UNCONT'); INTEGRAL%UNCONT = .TRUE.
        CASE ('.NOSEGMENT'); INTEGRAL%NOSEGMENT = .TRUE.
        CASE ('.MIXEDOVERLAP'); INTEGRAL%MIXEDOVERLAP = .TRUE.
        CASE ('.NO MM_FILES'); INTEGRAL%NO_MMFILES = .TRUE.
        CASE ('.DENSFIT'); integral%DENSFIT = .TRUE.                 
        CASE ('.RUNMM');  integral%FMM = .TRUE.
        CASE ('.DECGRA'); integral%run_dec_gradient_test=.true.
        CASE DEFAULT
           WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                & '" not recognized in **INTEGRALS readin.'
           print*,'Keyword ',WORD
           !CALL lsQUIT('Illegal keyword in **INTEGRAL.',lupri)
        END SELECT
     ENDIF
     IF (WORD(1:2) == '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(PROMPT(1:1) .EQ. '*')THEN
        SELECT CASE(WORD(1:7))
        CASE('*READIN'); WRITE(LUPRI,*)'*READIN is not yet implemented in new int driver'
        CASE('*ONEINT'); WRITE(LUPRI,*)'*ONEINT is not yet implemented in new int driver'
        CASE('*TWOINT'); WRITE(LUPRI,*)'*TWOINT is not yet implemented in new int driver'
        CASE('*SUPINT'); WRITE(LUPRI,*)'*SUPINT is not yet implemented in new int driver'
        CASE('*ER2INT'); WRITE(LUPRI,*)'*ER2INT is not yet implemented in new int driver'
        CASE('*SORINT'); WRITE(LUPRI,*)'*SORINT is not yet implemented in new int driver'
        CASE('*DENFIT') 
           CALL READ_INTEGRALS_DENFIT_INPUT(LUPRI,LUCMD,integral,WORD)
           READWORD=.FALSE.
        CASE('*FMM   ')
           CALL READ_INTEGRALS_FMM_INPUT(LUPRI,LUCMD,integral,WORD)
           READWORD=.FALSE.
        CASE('*FCK3  '); 
           CALL READ_INTEGRALS_FCK3_INPUT(LUPRI,LUCMD,integral,WORD)
           READWORD=.FALSE.
           WRITE(LUPRI,*)'WORD LEAVING FCK3 ',WORD
        CASE('*PERIOD'); WRITE(LUPRI,*)'*PERIOD is not yet implemented in new int driver'
        CASE('*RIFOCK'); WRITE(LUPRI,*)'*RIFOCK is not yet implemented in new int driver'
        CASE('*END OF'); WRITE(LUPRI,*)'FOUND *END OF:'           
        CASE DEFAULT
           WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                & '" not recognized in **INTEGRALS'
                          CALL lsQUIT('Illegal keyword in INTEGRALS.',lupri)
        END SELECT
     ENDIF
  ENDDO
END subroutine INTEGRAL_INPUT

!> \brief Read the $INFO section under *LINSCA in input file DALTON.INP and set configuration structure accordingly.
!> \author S. Host
!> \date March 2010
SUBROUTINE config_info_input(config,lucmd)
  implicit none
  !> Contains info, settings and data for entire calculation
  type(configItem),intent(inout) :: config
  !> Logical unit number for DALTON.INP
  integer,intent(in) :: lucmd
  character(len=40) :: word
  integer :: i

  READ (LUCMD, '(A40)') word
  DO
     if (WORD(1:12) == '$END INFO') exit
     IF (WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#') THEN
         READ (LUCMD, '(A40)') word
         CYCLE
     endif
     SELECT CASE(WORD)
     CASE('DEBUG_ARH_LINTRA')
          config%solver%DEBUG_ARH_LINTRA = .true.
     CASE('DEBUG_ARH_PRECOND')
          config%solver%DEBUG_ARH_PRECOND = .true.
     CASE('DEBUG_CONVERT')
           config%opt%DEBUG_CONVERT = .true.
           READ(LUCMD,*) config%opt%cfg_which_conversion 
     !CASE('DEBUG_DCHANGE')
     !     DEBUG_DCHANGE = .true.
     CASE('DEBUG_DD')
          READ(LUCMD,*) config%solver%cfg_nits_debug
          config%decomp%cfg_check_converged_solution = .true.
          config%solver%DEBUG_DD = .true.
     CASE('DEBUG_DD_LINTRA')
          config%solver%DEBUG_DD_LINTRA = .true.
     CASE('DEBUG_DD_HOMOLUMO')
          config%solver%DEBUG_DD_HOMOLUMO = .true.
     CASE('DEBUG_DIAG_REDSPACE')
          config%solver%DEBUG_DIAG_REDSPACE = .true.
     CASE('DEBUG_DIAG_HESSIAN')
          config%opt%DEBUG_DIAG_HESSIAN = .true.
     !CASE('DEBUG_DSM_DCHANGE')
     !     DEBUG_DSM_DCHANGE = .true.
     !CASE('DEBUG_DSM_EMODEL')
     !     DEBUG_DSM_EMODEL = .true.
     !CASE('DEBUG_DSM_METRIC')
     !     DEBUG_DSM_METRIC = .true.
     !CASE('DEBUG_EMODEL_CHANGE')
     !     DEBUG_EMODEL_CHANGE = .true.
     CASE('DEBUG_HESSIAN')
          config%solver%DEBUG_HESSIAN = .true.
     CASE('DEBUG_HESSIAN_EXACT')
          config%solver%DEBUG_HESSIAN_EXACT = .true. ; config%solver%DEBUG_HESSIAN = .true.
     CASE('DEBUG_IDEMPOTENCY')
          config%diag%DEBUG_IDEMPOTENCY = .true.
     !CASE('DEBUG_OAO_GRADIENT')
     !     DEBUG_OAO_GRADIENT = .true.
     CASE('DEBUG_RH_DSM_ECHANGE')
          config%diag%DEBUG_RH_DSM_ECHANGE = .true.
     !CASE('DEBUG_RH_MU_E')
     !     DEBUG_RH_MU_E = .true.
     !     READ(LUCMD,*) cfg_nits_debug,  cfg_mu_max_debug 
     !     READ(LUCMD,*) (cfg_its_debug(i),i=1,cfg_nits_debug)
     !CASE('DEBUG_RSP')
     !     DEBUG_RSP = .true.
     !CASE('DEBUG_RSP_LINSCA')
     !     DEBUG_RSP_LINSCA = .true.
     CASE('INFO_CROP')
          config%solver%INFO_CROP = .true.
     CASE('INFO_LEVELSHIFT')
          config%solver%INFO_LEVELSHIFT    = .true.
     CASE('INFO_STABILITY')
          config%decomp%INFO_STABILITY      = .true.
     CASE('INFO_DIIS')
          config%av%INFO_DIIS = .true.
          config%av%INFO_WEIGHT_FINAL = .true. 
     CASE('INFO_DSM')
          config%av%INFO_DSM_EIGENVAL     = .true.
          config%av%INFO_DSM_ENERGY       = .true.
          config%av%INFO_DSM_EXIT         = .true.
          config%av%INFO_DSM_PROJ         = .true.
          config%av%INFO_DSM_STEP_TOTAL   = .true.
          config%av%INFO_WEIGHT_FINAL     = .true.
     CASE('INFO_DSM_DETAIL')
          config%av%INFO_D_PROJ           = .true.
          config%av%INFO_DSM_CNORM_MU_FIG = .true.
          !config%av%INFO_DSM_DELTA        = .true.
          !config%av%INFO_DSM_DERIVATIVES  = .true.
          config%av%INFO_DSM_EIGENVAL     = .TRUE.
          config%av%INFO_DSM_ENERGY       = .TRUE.
          config%av%INFO_DSM_EQ           = .true.
          config%av%INFO_DSM_EXIT         = .TRUE.
          config%av%INFO_DSM_GRAD         = .true.
          config%av%INFO_DSM_METRIC       = .true.
          config%av%INFO_DSM_NIT          = .true.
          config%av%INFO_DSM_PROJ         = .TRUE.
          config%av%INFO_DSM_RATIO        = .true.
          config%av%INFO_DSM_STEP         = .true.
          config%av%INFO_DSM_STEP_BRACKET = .true.
          config%av%INFO_DSM_STEP_TOTAL   = .TRUE.
          config%av%INFO_DSM_TRUSTR       = .true.
          config%av%INFO_WEIGHTS          = .true.
     CASE('INFO_LINEQ')
          config%solver%INFO_LINEQ = .true.
     CASE('INFO_MATOP')
          config%opt%INFO_MATOP = .true.
     !CASE('INFO_ORB_E')
     !     INFO_ORB_E = .true.
     CASE('INFO_STABILITY_REDSPACE')
           config%decomp%info_stability_redspace = .true.
     CASE('INFO_RH')
          config%diag%INFO_RH_ITERATIONS    = .true.
          config%diag%INFO_RH_MU            = .true. 
          config%diag%INFO_RH_GAP           = .true.
     CASE('INFO_RH_DETAIL')
          !config%diag%INFO_RH_EPRED         = .true.
          config%diag%INFO_RH_GAP           = .true.
          config%diag%INFO_RH_GRADIENT      = .true.
          config%diag%INFO_RH_ITERATIONS    = .true.
          config%diag%INFO_RH_MU            = .true.
          config%diag%INFO_RH_MOCHANGE      = .true.
     CASE('INFO_RSP')
        config%response%rspsolverinput%INFO_RSP = .true.
     CASE('INFO_RSP_REDSPACE')
          config%response%rspsolverinput%INFO_RSP_REDSPACE = .true.
     !CASE('INFO_TIME_MAT_OPERATIONS')
     !     INFO_TIME_MAT_OPERATIONS = .true.
     !     call mat_timings
     CASE('INFO_WEIGHT')
          config%av%INFO_WEIGHT_FINAL = .true.
     CASE DEFAULT
          WRITE(config%LUPRI,*) ' Keyword ',WORD,' not recognized in config_info_input'
          CALL lsQUIT('Illegal keyword in config_info_input.',config%lupri)
     END SELECT
     READ (LUCMD, '(A40)') word
  ENDDO

  !if (INFO_TIME_MAT_OPERATIONS) call mat_timings
  
END SUBROUTINE config_info_input



!> \brief Read the **RESPONS section under *LINSCA in input file DALTON.INP and set configuration structure accordingly.
!> \author S. Host
!> \date March 2010
SUBROUTINE config_rsp_input(config,lucmd,readword)
  USE RSP_type_module
  implicit none
  !> Logical for keeping track of when to read
  LOGICAL,intent(inout)                :: READWORD
  !> Contains info, settings and data for entire calculation, including response
  type(configItem),intent(inout) :: config
  !> Logical unit number for DALTON.INP
  integer,intent(in) :: lucmd
  character(len=40) :: word
  character(len=8) :: labels(2)
  character(len=8) :: QRlabels(3)
  integer :: i,j,k,n,nops,nlabel
  logical :: xcomp,ycomp,zcomp
!temporary variables, fix to make test cases run!
  logical :: cfg_run_pdbs, cfg_rsp_run_quadratic
  integer :: Bterm_nr=0
  integer,allocatable     :: Bterm_index(:)
  integer :: STEPS
  real(realk) :: TMP1,TMP2,STEP,TMP3

  cfg_run_pdbs = .false.
  cfg_rsp_run_quadratic = .false.
  nops = 0
  DO
     IF(READWORD) THEN
        READ (LUCMD, '(A40)') WORD
        READWORD=.TRUE.
     ENDIF
     IF ((WORD(1:1) .EQ. '!') .OR. (WORD(1:1) .EQ. '#')) CYCLE
     IF(WORD(1:2) .EQ. '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(WORD(1:13) == '*END OF INPUT') THEN
        backspace(LUCMD)
        EXIT
     END IF

     if (WORD(1:1) == '*') then
       !which type of response is wanted??
       SELECT CASE(WORD)
       ! Kasper K
       CASE('*ALPHA')
           config%response%tasks%doALPHA=.true.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.BFREQ')
                 config%response%alphainput%real_frequencies_in_input = .true.
                 READ(LUCMD,*) config%response%alphainput%nfreq 
                 allocate(config%response%alphainput%bfreq(config%response%alphainput%nfreq))
                 read(LUCMD,*) config%response%alphainput%bfreq(1:config%response%alphainput%nfreq)
              CASE('.IMBFREQ')
                 config%response%alphainput%imag_frequencies_in_input = .true.
                 READ(LUCMD,*) config%response%alphainput%nimfreq 
                 allocate(config%response%alphainput%imbfreq(config%response%alphainput%nimfreq))
                 read(LUCMD,*) config%response%alphainput%imbfreq(1:config%response%alphainput%nimfreq)
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *ALPHA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*BETA')
           config%response%tasks%doBETA=.true.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.BFREQ')
                 config%response%betainput%real_bfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%betainput%nbfreq
                 allocate(config%response%betainput%bfreq(config%response%betainput%nbfreq))
                 read(LUCMD,*) config%response%betainput%bfreq(1:config%response%betainput%nbfreq)
              CASE('.IMBFREQ')
                 config%response%betainput%imag_bfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%betainput%nimbfreq 
                 allocate(config%response%betainput%imbfreq(config%response%betainput%nimbfreq))
                 read(LUCMD,*) config%response%betainput%imbfreq(1:config%response%betainput%nimbfreq)
              CASE('.CFREQ')
                 config%response%betainput%real_cfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%betainput%ncfreq
                 allocate(config%response%betainput%cfreq(config%response%betainput%ncfreq))
                 read(LUCMD,*) config%response%betainput%cfreq(1:config%response%betainput%ncfreq)
              CASE('.IMCFREQ')
                 config%response%betainput%imag_cfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%betainput%nimcfreq 
                 allocate(config%response%betainput%imcfreq(config%response%betainput%nimcfreq))
                 read(LUCMD,*) config%response%betainput%imcfreq(1:config%response%betainput%nimcfreq)
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *BETA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*GAMMA')
           config%response%tasks%doGAMMA=.true.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.BFREQ')
                 config%response%gammainput%real_bfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%nbfreq
                 allocate(config%response%gammainput%bfreq(config%response%gammainput%nbfreq))
                 read(LUCMD,*) config%response%gammainput%bfreq(1:config%response%gammainput%nbfreq)
              CASE('.IMBFREQ')
                 config%response%gammainput%imag_bfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%nimbfreq 
                 allocate(config%response%gammainput%imbfreq(config%response%gammainput%nimbfreq))
                 read(LUCMD,*) config%response%gammainput%imbfreq(1:config%response%gammainput%nimbfreq)
              CASE('.CFREQ')
                 config%response%gammainput%real_cfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%ncfreq
                 allocate(config%response%gammainput%cfreq(config%response%gammainput%ncfreq))
                 read(LUCMD,*) config%response%gammainput%cfreq(1:config%response%gammainput%ncfreq)
              CASE('.IMCFREQ')
                 config%response%gammainput%imag_cfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%nimcfreq 
                 allocate(config%response%gammainput%imcfreq(config%response%gammainput%nimcfreq))
                 read(LUCMD,*) config%response%gammainput%imcfreq(1:config%response%gammainput%nimcfreq)
              CASE('.DFREQ')
                 config%response%gammainput%real_dfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%ndfreq
                 allocate(config%response%gammainput%dfreq(config%response%gammainput%ndfreq))
                 read(LUCMD,*) config%response%gammainput%dfreq(1:config%response%gammainput%ndfreq)
              CASE('.IMDFREQ')
                 config%response%gammainput%imag_dfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%nimdfreq 
                 allocate(config%response%gammainput%imdfreq(config%response%gammainput%nimdfreq))
                 read(LUCMD,*) config%response%gammainput%imdfreq(1:config%response%gammainput%nimdfreq)
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *GAMMA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*TPA')
           config%response%tasks%doTPA=.true.
           ! Sanity check 1. We only run TPA calculations if .NEXCIT have been specified!
           TPASanityCheck1: if(config%decomp%cfg_rsp_nexcit == 0) then
              write(config%lupri,*) 'Error in TPA input: To run a TPA calculation you must&
                   & define the maximum number of excited states using the .NEXCIT&
                   & keyword under **RESPONSE before specifying *TPA!'
              write(config%lupri,*) 'Example of TPA response input for the first five excited states:'
              write(config%lupri,*) '**RESPONSE'
              write(config%lupri,*) '.NEXCIT'
              write(config%lupri,*) '5'
              write(config%lupri,*) '*TPA'
                  call lsQUIT('Error in TPA input: To run a TPA calculation you must&
                 & define the maximum number of excited states using the .NEXCIT&
                 & keyword under **RESPONSE before specifying *TPA!&
                 & See end of output file for input example!', config%lupri)
           end if TPASanityCheck1
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.EXSTATES')
                 config%response%tpainput%specific_states_in_input = .true.
                 READ(LUCMD,*) config%response%tpainput%tpa_nexci
                 allocate(config%response%tpainput%ExStates(config%response%tpainput%tpa_nexci))
                 read(LUCMD,*) config%response%tpainput%ExStates(1:config%response%tpainput%tpa_nexci)
                 ! Sanity check 2 
                 do i=1,config%response%tpainput%tpa_nexci
                      TPASanityCheck2: if(config%response%tpainput%ExStates(i) > &
                                &config%decomp%cfg_rsp_nexcit) then
                   write(config%lupri,*) 'Error in .EXSTATES in *TPA input: The maximum excited state&
                          & defined by .EXSTATES must be less than or equal to the number of&
                          & excited states defined by .NEXCIT under **RESPONSE.'
                   write(config%lupri,*) 'Example of TPA response input for the four excited states 1,3, 5, and 8:'
                   write(config%lupri,*) '**RESPONSE'
                   write(config%lupri,*) '.NEXCIT'
                   write(config%lupri,*) '8'
                   write(config%lupri,*) '*TPA'
                   write(config%lupri,*) '.EXSTATES'
                   write(config%lupri,*) '4'
                   write(config%lupri,*) '1 3 5 8'
                  CALL lsQUIT('Error in .EXSTATES in *TPA input: The maximum excited state&
                     & defined by .EXSTATES must be less than or equal the number of excited states defined by&
                     & .NEXCIT under **RESPONSE!&
                     & See end of output file for input example!', config%lupri)
               end if TPASanityCheck2
                end do
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *TPA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*DAMPED_TPA')
           config%response%tasks%doDTPA=.true.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.OPFREQ') ! One-photon frequencies
                 READ(LUCMD,*) config%response%dtpainput%nfreq
                 allocate(config%response%dtpainput%freq(config%response%dtpainput%nfreq))
                 READ(LUCMD,*) config%response%dtpainput%freq(1:config%response%dtpainput%nfreq)
              CASE('.GAMMA') ! One-photon frequencies
                 config%response%dtpainput%gamma_specified=.true.
                 READ(LUCMD,*) config%response%dtpainput%gamma
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *DAMPED_TPA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*ESDIPOLE')
           config%response%tasks%doESD=.true.
           ! Sanity check 1. We only run ESD calculations if .NEXCIT have been specified!
           ESDSanityCheck1: if(config%decomp%cfg_rsp_nexcit == 0) then
              write(config%lupri,*) 'Error in excited state gradient input: &
                   & To run an ESD calculation you must&
                   & define the maximum number of excited states using the .NEXCIT&
                   & keyword under **RESPONSE before specifying *ESDIPOLE!'
              write(config%lupri,*) 'Example of ESD response input for the first five excited states:'
              write(config%lupri,*) '**RESPONSE'
              write(config%lupri,*) '.NEXCIT'
              write(config%lupri,*) '5'
              write(config%lupri,*) '*ESDIPOLE'
                  call lsQUIT('Error in excited state gradient input: To run an ESD calculation you must&
                 & define the maximum number of excited states using the .NEXCIT&
                 & keyword under **RESPONSE before specifying *ESDIPOLE!&
                 & See end of output file for input example!', config%lupri)
               end if ESDSanityCheck1
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.EXSTATES')
                 config%response%esdinput%specific_states_in_input = .true.
                 READ(LUCMD,*) config%response%esdinput%esd_nexci
                 allocate(config%response%esdinput%ExStates(config%response%esdinput%esd_nexci))
                 read(LUCMD,*) config%response%esdinput%ExStates(1:config%response%esdinput%esd_nexci)
                 ! Sanity check 2 
                 do i=1,config%response%esdinput%esd_nexci
                      ESDSanityCheck2: if(config%response%esdinput%ExStates(i) > &
                              &config%decomp%cfg_rsp_nexcit) then
                   write(config%lupri,*) 'Error in .EXSTATES in * input: The maximum excited state&
                          & defined by .EXSTATES must be less than or equal to the number of&
                          & excited states defined by .NEXCIT under **RESPONSE.'
                   write(config%lupri,*) 'Example of excited state gradient response input &
                                &for the four excited states 1,3, 5, and 8:'
                   write(config%lupri,*) '**RESPONSE'
                   write(config%lupri,*) '.NEXCIT'
                   write(config%lupri,*) '8'
                   write(config%lupri,*) '*ESDIPOLE'
                   write(config%lupri,*) '.EXSTATES'
                   write(config%lupri,*) '4'
                   write(config%lupri,*) '1 3 5 8'
                  CALL lsQUIT('Error in .EXSTATES in *ESDIPOLE input: The maximum excited state&
                     & defined by .EXSTATES must be less than or equal the number of excited states defined by&
                     & .NEXCIT under **RESPONSE!&
                     & See end of output file for input example!', config%lupri)
               end if ESDSanityCheck2
                end do
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *ESDIPOLE input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*ESGRAD')
           config%response%tasks%doESG=.true.
           ! Sanity check 1. We only run ESG calculations if .NEXCIT have been specified!
           ESGSanityCheck1: if(config%decomp%cfg_rsp_nexcit == 0) then
              write(config%lupri,*) 'Error in excited state gradient input: &
                   & To run an ESG calculation you must&
                   & define the maximum number of excited states using the .NEXCIT&
                   & keyword under **RESPONSE before specifying *ESGRAD!'
              write(config%lupri,*) 'Example of ESG response input for the first five excited states:'
              write(config%lupri,*) '**RESPONSE'
              write(config%lupri,*) '.NEXCIT'
              write(config%lupri,*) '5'
              write(config%lupri,*) '*ESGRAD'
                  call lsQUIT('Error in excited state gradient input: To run an ESG calculation you must&
                 & define the maximum number of excited states using the .NEXCIT&
                 & keyword under **RESPONSE before specifying *ESGRAD!&
                 & See end of output file for input example!', config%lupri)
               end if ESGSanityCheck1
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.EXSTATES')
                 config%response%esginput%specific_states_in_input = .true.
                 READ(LUCMD,*) config%response%esginput%esg_nexci
                 allocate(config%response%esginput%ExStates(config%response%esginput%esg_nexci))
                 read(LUCMD,*) config%response%esginput%ExStates(1:config%response%esginput%esg_nexci)
                 ! Sanity check 2 
                 do i=1,config%response%esginput%esg_nexci
                      ESGSanityCheck2: if(config%response%esginput%ExStates(i) > &
                              &config%decomp%cfg_rsp_nexcit) then
                   write(config%lupri,*) 'Error in .EXSTATES in * input: The maximum excited state&
                          & defined by .EXSTATES must be less than or equal to the number of&
                          & excited states defined by .NEXCIT under **RESPONSE.'
                   write(config%lupri,*) 'Example of excited state gradient response input &
                                &for the four excited states 1,3, 5, and 8:'
                   write(config%lupri,*) '**RESPONSE'
                   write(config%lupri,*) '.NEXCIT'
                   write(config%lupri,*) '8'
                   write(config%lupri,*) '*ESGRAD'
                   write(config%lupri,*) '.EXSTATES'
                   write(config%lupri,*) '4'
                   write(config%lupri,*) '1 3 5 8'
                  CALL lsQUIT('Error in .EXSTATES in *ESGRAD input: The maximum excited state&
                     & defined by .EXSTATES must be less than or equal the number of excited states defined by&
                     & .NEXCIT under **RESPONSE!&
                     & See end of output file for input example!', config%lupri)
               end if ESGSanityCheck2
                end do
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *ESGRAD input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Thomas K and Kasper K
       CASE('*MOLGRA')
                    WRITE(config%LUPRI,*) 'Ground state molecular gradient calculations are carried out.'
                    config%response%tasks%dograd = .True.
             ! Joanna K
        CASE('*SOLVER')
            do
               READ(LUCMD,'(A40)') word
               if(word(1:1) == '!' .or. word(1:1) == '#') cycle
               if(word(1:1) == '*') then ! New property or *END OF INPUT
                  backspace(LUCMD)
                  exit
               end if
               SELECT CASE(word)
               CASE('.COMPLEX')
                  config%response%rspsolverinput%rsp_complex = .true.
                  READ(LUCMD,*) config%response%rspsolverinput%rsp_gamma
               CASE('.CONVTHR') 
                  READ(LUCMD,*) config%response%rspsolverinput%rsp_thresh
               CASE('.NEW_SOLVER')
                  config%response%rspsolverinput%rsp_new_solver = .true.
               CASE('.MAXIT')
                  READ(LUCMD,*) config%response%rspsolverinput%rsp_maxit 
                 config%response%rspsolverinput%rsp_maxred=2*config%response%rspsolverinput%rsp_maxit 
               CASE('.S_NORM')
                  config%response%rspsolverinput%rsp_single_norm =.true.
               CASE('.CONVDYN')
                  READ(LUCMD,*) config%response%rspsolverinput%rsp_convdyn_type
                  config%response%rspsolverinput%rsp_convdyn =.true.
                  SELECT CASE(config%response%rspsolverinput%rsp_convdyn_type)
                  CASE('TIGHT'); config%response%rspsolverinput%rsp_conv_factor = 1.0d-3
                  CASE('STAND'); config%response%rspsolverinput%rsp_conv_factor = 1.0d-2
                  CASE('SLOPP'); config%response%rspsolverinput%rsp_conv_factor = 1.0d-1
                  CASE DEFAULT
                  WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',config%response%rspsolverinput%rsp_convdyn_type,&
                                    & '" not recognized with .CONVDYN'
                  WRITE (config%LUPRI,*) 'Options are TIGHT, STANDARD, and SLOPPY.'
                  CALL lsQUIT('Illegal keyword with .CONVDYN',config%lupri)
                  END SELECT
               CASE('.OLSEN')
                  config%response%rspsolverinput%rsp_olsen = .true.
               CASE('.QUIET')
                  config%response%rspsolverinput%rsp_quiet = .true.
               CASE('.AOPREC')
                  config%response%rspsolverinput%rsp_mo_precond = .false.
               CASE('.AOSTART')
                  config%response%rspsolverinput%rsp_mostart = .false.
               CASE('.NOPREC')
                  config%response%rspsolverinput%rsp_no_precond = .true.
               CASE('.RESTEXC')
                  READ (LUCMD,*) config%response%rspsolverinput%rsp_restart_nexci
                  config%response%rspsolverinput%rsp_restart_exci = .true.
               CASE ('.NSTART');   READ(LUCMD,*) config%response%rspsolverinput%rsp_no_of_startvectors
                  config%response%rspsolverinput%rsp_startvectors = .true.  
                  config%decomp%cfg_startvectors = .TRUE.
               CASE('.TWOSTART')
                  config%response%rspsolverinput%rsp_damp_2start=.true.
               CASE DEFAULT
                  WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                       & '" not recognized in RESPONSE *SOLVER input.'
                  CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
                 
               END SELECT
            enddo


! KASPER K: IMPORTANT NOTE!!!
! As far as I know all response input keywords below 
! (except .NEXCIT and perhaps some solver thresholds) will be
! obsolote in LSDALTON and should be removed for the release.
! If anyone has objections please let me know!



       ! kasperk
       CASE('*AORESPONSE')
                    WRITE(config%LUPRI,*) 'AO-based response function &
                                   &calculations are carried out...'
                    !cfg_run_AOresponse=.true.
                    !allocate(AOR_input_temp(100))
                    !i=1
                    !do
                    !  READ(LUCMD,'(A40)') AOR_input_i
                    !    if(AOR_input_i(1:1) == '!' .or. AOR_input_i(1:1) == '#') cycle
                    !    if(AOR_input_i(1:15) == '*END AORESPONSE') exit
                    !    AOR_input_temp(i) = AOR_input_i
                    !    if(i>1) then
                    !      if(AOR_input_temp(i-1)(1:7) == '.NEXCIT') then
                    !        read(AOR_input_i, '(I4)') nexci_temp
                    !        if(config%decomp%cfg_rsp_nexcit < nexci_temp) config%decomp%cfg_rsp_nexcit = nexci_temp
                    !      endif
                    !    endif
                    !    i=i+1
                    !enddo
                    !AOR_stacksize=i-1
                    !allocate(AOR_input(AOR_stacksize))
                    !AOR_input(1:AOR_stacksize) = AOR_input_temp(1:AOR_stacksize)
                    !deallocate(AOR_input_temp)
       CASE('*QUASIMCD')
           config%response%tasks%doMCD=.true.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:15) == '*END QUASIMCD') exit
              SELECT CASE(word)
              CASE('.DAMPEDXCOOR')
                 READ(LUCMD,*) config%response%MCDinput%nXcoor
                 call mem_alloc(config%response%MCDinput%Xcoor,config%response%MCDinput%nXcoor)
                 DO I=1,config%response%MCDinput%nXcoor
                    READ(LUCMD,*) config%response%MCDinput%Xcoor(I)
                 ENDDO
              CASE('.DAMPEDRANGE')
                 READ(LUCMD,*) TMP1
                 READ(LUCMD,*) TMP2
                 READ(LUCMD,*) STEP
                 STEPS=1
                 TMP3=TMP1
                 DO WHILE(TMP3 < TMP2)
                    STEPS=STEPS+1
                    TMP3=TMP3+STEP
                 ENDDO
                 config%response%MCDinput%nXcoor = STEPS
                 call mem_alloc(config%response%MCDinput%Xcoor,STEPS)
                 config%response%MCDinput%Xcoor(1) = TMP1
                 TMP3=TMP1
                 DO I=2,config%response%MCDinput%nXcoor
                    TMP3=TMP3+STEP
                    config%response%MCDinput%Xcoor(I) = TMP3
                 ENDDO
              CASE('.MCDEXCIT')
                 READ(LUCMD,*) config%response%MCDinput%nexci
              CASE('.NO LONDON')
                 config%response%MCDinput%london=.FALSE.
              CASE('.NO NONLONDON')
                 config%response%MCDinput%nolondon=.FALSE.                 
              CASE('.NO SIMULATE')
                 config%response%MCDinput%simulate=.FALSE.    
              CASE('.NO ATERM')
                 config%response%MCDinput%doAterms =.FALSE.    
              CASE('.NO BTERM')
                 config%response%MCDinput%doBterms =.FALSE.    
              CASE('.NO DAMPED')
                 config%response%MCDinput%dampedMCD =.FALSE.   
              CASE('.GAUSSIAN')
                 config%response%MCDinput%lorentz=.FALSE.                 
              CASE('.LINESHAPEPARAM')
                 config%response%MCDinput%useinputgamma=.TRUE.                 
                 READ(LUCMD,*) config%response%MCDinput%Gamma 
              CASE('.NVECFORPEAK')
                 READ(LUCMD,*) config%response%MCDinput%nVecForPeak
              CASE('.NSTEPS')
                 READ(LUCMD,*) config%response%MCDinput%Nsteps
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *QUASIMCD input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       CASE('*LINRSP')
                    config%response%tasks%dolinrsp = .True.
                    WRITE(config%LUPRI,*) 'Linear response calculations are carried out.'
                    nlabel = 2
       !SONIANEW
       CASE('*ESGRSP')
                    WRITE(config%LUPRI,'(/A)') 'Excited state gradient is available in LSDALTON&
                         & only using the *ESGRAD keyword!'
                    CALL lsQUIT('Use *ESGRAD to get the excited state gradient!',config%lupri)
                    !WRITE(LUPRI,*) 'Excited state gradient calculations are carried out'
                    !cfg_rsp_run_exgr = .true.
                    !cfg_rsp_grad_purify = .true.
                    !cfg_start_guess = 'H1DIAG'    
                    !IF(cfg_do_in_oao) cfg_start_guess = 'H1OAO'    
                    !nlabel = 1
                    !WRITE(LUPRI,'(A)')'The Calculation uses the density matrix from the previous gometry step'
                    !WRITE(LUPRI,'(A)')'as an initial guess in the new geometry step.'
                    !WRITE(LUPRI,'(A,A,A)')'If this fails the fallback is ',cfg_start_guess,' As Huckel has not been implemented'
       !SONIANEW
       CASE('*HESSMAG')
                    WRITE(config%LUPRI,*) 'Hessian/Magnetizability calcs carried out'
                    nlabel = 2
                    !cfg_rsp_run_mag = .true.
                    !!Sonia: Make separate set of options: use_eq_79 etc..
                    !!cfg_rsp_run_hes = .true.
       !THOMAS_NEW
       CASE('*SHIELD')
                    WRITE(config%LUPRI,*)'You have requeste a calculation of the &
                    & Nuclear Magnetic Shielding tensor, unfortunately the code&
                    & has been corrupted which is what happens when you do not&
                    & make a test case. Let that be a warning. Contact me if&
                    & this is really important. Thomas Kjaergaard'
                    CALL FLUSH(config%LUPRI)
                    CALL lsQUIT('Code Corrupted *SHIELD nolonger valid keyword, contact Thomas Kjaergaard',-1)
                    !WRITE(LUPRI,*) 'Nuclear Magnetic Shielding&
                    ! & calculations are carried out'
                    !cfg_rsp_run_shield = .true.
       CASE('*PDBS')
                    WRITE(config%LUPRI,*) 'Pertubation dependent basis set &
                    & calculations are carried out'
                    cfg_run_PDBS = .true.
                    if (cfg_rsp_run_quadratic) then
                       WRITE(config%LUPRI,*) 'You cannot run Quadratic&
                       & Response Calculation and a Perturbation dependent &
                       & basis set calculation at the same time' 
                       CALL lsQUIT('REMOVE *PDBS or *QUADRSP',-1)
                    endif
       CASE('*QUIET')
                    !cfg_rsp_quiet=.TRUE.
       CASE DEFAULT
          WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
             & '" not recognized in config_rsp_input.'
          !CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri) 
       END SELECT

     else ! keywords not starting with *
       SELECT CASE(WORD)
       CASE('.NEXCIT') 
                       READ (LUCMD,*) config%decomp%cfg_rsp_nexcit
! KK, in LSDALTON_ONLY calculate one-photon absorption strengths for all
! excitation energies by default when .NEXCIT is specified.
                       config%response%tasks%doOPA=.true.
       CASE('.FREQ')
                       !ALLOCATE(cfg_rsp_freq(10))
                       !READ(LUCMD,*) cfg_rsp_nfreq,(cfg_rsp_freq(i),i=1,cfg_rsp_nfreq)
                       !cfg_rsp_freqs_specified = .true.
       CASE('.IMFREQ')
                       !ALLOCATE(cfg_rsp_imfreq(10))
                       !READ(LUCMD,*) cfg_rsp_nimfreq,(cfg_rsp_imfreq(i),i=1,cfg_rsp_nimfreq)
                       !cfg_rsp_imfreqs_specified = .true.
       CASE('.DAMP_FREQ')
                       !READ(LUCMD,*) cfg_rsp_nfreq,cfg_rsp_minfreq,cfg_rsp_maxfreq
                       !ALLOCATE(cfg_rsp_freq(cfg_rsp_nfreq))
                       !if( (cfg_rsp_minfreq/=cfg_rsp_maxfreq) .and.&
                       !& (cfg_rsp_nfreq==1)) then
                       !   WRITE (LUPRI,'(/,3A,/)') ' Error. Wrong number of frequencies. '
                       !   CALL lsQUIT('Error in damp_freq')
                       !else
                       !   cfg_rsp_freq(1)=cfg_rsp_minfreq
                       !   cfg_rsp_freq(cfg_rsp_nfreq)=cfg_rsp_maxfreq
                       !   if (cfg_rsp_nfreq>2) then
                       !      k=0
                       !      do i=2,cfg_rsp_nfreq-1
                       !        k=k+1
                       !       cfg_rsp_freq(i)=cfg_rsp_minfreq+k*((cfg_rsp_maxfreq-cfg_rsp_minfreq)/(cfg_rsp_nfreq-1))
                       !        enddo
                       !        endif
                       !cfg_rsp_freqs_specified = .true.
                       !endif
       CASE('.DAMP_TPA')   !kasperk
! 1. Determine isotropically averaged damped TPA. Input:
! .DAMP_TPA
! 0
! 
! 2. Determine only specific TPA component(s), e.g. XXXX and XXYY. Input example:
! .DAMP_TPA
! 2                 ! Number of components
! XXXX              ! stored in tp_info(1)
! XXYY              ! stored in tp_info(2)
                       !cfg_rsp_run_damped_tpa = .true.
                       !READ(LUCMD,*) Ntp_comp
                       !if(Ntp_comp /= 0) then
                       !  allocate(tp_info(Ntp_comp))
                       !    do i=1,Ntp_comp
                       !      READ(LUCMD,*) tp_info(i)
                       !    enddo
                       !endif
       CASE('.OPERATOR')
                       !WRITE(LUPRI,*)
                       !WRITE(LUPRI,*) 'Response functions/transition moments will be calculated'
                       !WRITE(LUPRI,*) 'for the following (pairs of) operators:'
                       !WRITE(LUPRI,*) 
                       !do 
                       !  READ (LUCMD, '(a40)') word
                       !  if (word(1:1) == '.' .or. word(1:1) == '$') exit
                       !  if (word(1:1) == '*') exit
                       !  if (word(1:1) == '!' .or. word(1:1) == '#') cycle
                       !  nops = nops + 1
                       !  !TODO: check that it is within ..operators boundaries
                       !  !SONIANEW: initialized within given section (LINRSP or ESGRSP)
                       !  !nlabel = 2
                       !  call trim_string(word,nlabel,labels)
                       !  do i = 1,nlabel
                       !    cfg_rsp_operators(i,nops) = CFG_RSP_LABEL_TO_INT(labels(i))
                       !  enddo
                       !enddo
                       !cycle
       CASE('.DIPOLE')
                       !cfg_rsp_dipole = .true. 
       CASE('.QRSP OAO')
                       !cfg_QRSP_OAO=.true.
       CASE('.QRSP LRF')
                       !cfg_QR_do_LRF=.true.
       CASE('.ATERM')
                       !if(cfg_run_PDBS)then
                       !  cfg_run_Aterm=.true.
                       !  READ (LUCMD,*) Aterm_nr
                       !  !PDBSitems=PDBSitems+Aterm_nr
                       !  !PDBSitems set when it is examined how many degenerate
                       !  !states there is 
                       !  allocate(Aterm_index(Aterm_nr))
                       !  do i=1,Aterm_nr 
                       !     READ(LUCMD,*) Aterm_index(i)
                       !  enddo
                       !  IF(config%decomp%cfg_rsp_nexcit .EQ. 0)THEN
                       !     config%decomp%cfg_rsp_nexcit = Aterm_nr
                       !     DO i=1,Aterm_nr 
                       !        config%decomp%cfg_rsp_nexcit = MAX(config%decomp%cfg_rsp_nexcit,Aterm_index(i))
                       !     ENDDO
                       !  ENDIF
                       !endif
       CASE('.BTERM')
                       if(cfg_run_PDBS)then
                         !cfg_run_Bterm=.true.
                         READ (LUCMD,*) Bterm_nr
                         !PDBSitems=PDBSitems+Bterm_nr
                         allocate(Bterm_index(Bterm_nr))
                         do i=1,Bterm_nr 
                            READ(LUCMD,*) Bterm_index(i)
                         enddo
                         IF(config%decomp%cfg_rsp_nexcit .EQ. 0)THEN
                            config%decomp%cfg_rsp_nexcit = Bterm_nr
                            DO i=1,Bterm_nr 
                               config%decomp%cfg_rsp_nexcit = MAX(config%decomp%cfg_rsp_nexcit,Bterm_index(i))
                            ENDDO
                         ENDIF
                       endif
       CASE('.VERDET')
                       !if(cfg_run_PDBS)then
                       !  cfg_run_VERDET=.true.
                       !  READ (LUCMD,*) verdet_nr
                       !  PDBSitems=PDBSitems+verdet_nr
                       !  allocate(Verdet_freq(verdet_nr))
                       !  do i=1,verdet_nr 
                       !    READ(LUCMD,*) Verdet_freq(i)
                       !  enddo
                       !endif
       CASE('.PDBSINFO')
                       !cfg_PDBS_INFO=.true.
       CASE('.ESGRSP')
                       !cfg_run_Exgrad = .true. 
                       !READ (LUCMD,*) Exg_nr
                       !PDBSitems=PDBSitems+Exg_nr
                       !allocate(Exg_index(Exg_nr))
                       !do i=1,Exg_nr 
                       !  READ(LUCMD,*) Exg_index(i)
                       !enddo
       CASE('.RAMAN')
                       !cfg_run_RAMAN=.true.
                       !READ (LUCMD,*) raman_nr
                       !PDBSitems=PDBSitems+raman_nr
                       !allocate(Raman_freq(raman_nr))
                       !do i=1,raman_nr 
                       !  READ(LUCMD,*) Raman_freq(i)
                       !enddo
       CASE('.HBT')
                       !cfg_run_HerzbergTeller=.true.
                       !READ (LUCMD,*) HBT_nr
                       !PDBSitems=PDBSitems+HBT_nr
                       !allocate(HBT_index(HBT_nr))
                       !do i=1,HBT_nr 
                       !  READ(LUCMD,*) HBT_index(i)
                       !enddo
       CASE('.HOTFCHT')
                       !cfg_hotfcht = .true.
       !SONIANEW = select the excited state for the exc. state gradient
       CASE('.WHICHEXST') 
                       !READ (LUCMD,*) cfg_rsp_whichexst
       CASE('.RESTEXC')
                       !READ (LUCMD,*) cfg_restart_nexci; cfg_restart_exci = .true.
       !Joanna - rsp solver using conjugate gradient with optimal vectors
       CASE('.CGOP')
                      !cfg_rsp_cgop = .true.
                      !cfg_cg_truncate = .true.
                      !no_pairing=.true.
       CASE('.CGOP_NT')
                      !cfg_rsp_cgop = .true.
                      !cfg_cg_truncate = .false.
                      !no_pairing=.true.
       CASE('.CGOP_PAIR')      
                      !cfg_rsp_cgop = .true.       
       CASE('.CG_VEC')
                      !READ (LUCMD,*) cfg_red_truc
       CASE('.NO_CGOP_PREC')                 
               !cfg_cgop_prec=.false.
       !Joanna- complex response solver
       CASE('.COMPLEX')
                   READ (LUCMD,*) config%response%rspsolverinput%rsp_gamma
                         config%response%rspsolverinput%rsp_complex=.true.
                     !cfg_rsp_moprec= .true.
                     !cfg_rsp_mostart=.true.
       CASE('.COMP_OR')
                      !READ (LUCMD,*) cfg_rsp_gamma
                      !cfg_rsp_complex_or = .true.
       !FILIP:If requested, turn off the use of the density matrix from the previous geometry step
       !      as an initial guess for the new geometry step.
       !      gradients.
       CASE('.NODRESTART')
                      ! if (.not.cfg_rsp_run_grad) then
                      !    write(lupri,*)'.NODRESTART option requested, but *MOLGRA not turned on!'
                      !    call lsquit('.NODRESTART option requested, but *MOLGRA not turned on!')
                      ! endif
                      ! cfg_rsp_grad_drestart = .false.
       !FILIP:If requested, turn off McWeeny purification of the density matrix for the molecular
       !      gradients.
       CASE('.NOPURIFY')
                      ! if (.not.cfg_rsp_run_grad) then
                      !    write(lupri,*)'.NOPURIFY option requested, but *MOLGRA not turned on!'
                      !    call lsquit('.NOPURIFY option requested, but *MOLGRA not turned on!')
                      ! endif
                      ! cfg_rsp_grad_purify = .false.
       !FILIP:Change the threshold for the difference abs(2.0*trace-Nelectrons), below which 
       !      we disregard the McWeeny purified D and fall back to a standard initial guess.
       CASE('.THRNEL')
                      ! if (.not.cfg_rsp_run_grad) then
                      !    write(lupri,*)'.THRNEL option requested, but *MOLGRA not turned on!'
                      !    call lsquit('THRNEL option requested, but *MOLGRA not turned on!')
                      ! endif
                      ! READ (LUCMD,*) THRNEL
       CASE('.POLARIZ')
                      ! cfg_rsp_polariz = .true. !polarizability
       CASE('.HYPOLAR2')
                      ! cfg_rsp_hypolar2 = .true. !1st hyperpol., n+1
       CASE('.HYPOLAR')
                      ! cfg_rsp_hypolar = .true. !1st hyperpol., 2n+1
       CASE('.SECHYP3')
                      ! cfg_rsp_sechyp3 = .true. !2nd hyperpol., n+1
       CASE('.SECHYP')
                      ! cfg_rsp_sechyp = .true. !2nd hyperpol., 2n+1 (1+2+1)
       CASE('.SECHYP1')
                      ! cfg_rsp_sechyp1 = .true. !2nd hyperpol., 2n+1 (2+1+1)
       CASE('.MAGNET')
                      ! cfg_rsp_magnet = .true.   !ajt
       CASE('.EFGB')
                      ! cfg_rsp_efgb = .true.     !ajt
       CASE('.CME')
                      ! cfg_rsp_cme = .true.      !ajt
       CASE('.ROA')
                      ! cfg_rsp_roa = .true.      !ajt
       CASE('.CARS')
                      ! cfg_rsp_cars = .true.     !ajt
       CASE('.JONES')
                      ! cfg_rsp_jones = .true.    !ajt
       CASE('.VIBBETA')
                      ! cfg_rsp_vibbeta = .true.  !ajt
       CASE('.PROPTEST')
                      ! cfg_rsp_proptest = .true. !ajt
       CASE DEFAULT
          WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
               & '" not recognized in config_rsp_input'
          !CALL lsQUIT('Illegal keyword in config_rsp_input',config%lupri)
       END SELECT
     endif
  ENDDO

  if (config%decomp%cfg_rsp_nexcit > 0) then !Stinne
      !cfg_rsp_maxred = config%decomp%cfg_rsp_nexcit*cfg_rsp_maxred
  else
     !if (cfg_rsp_freqs_specified) then
     !   WRITE(LUPRI,*)
     !   WRITE(LUPRI, "('The response functions will be calculated for each of the ', i3, &
     !               & ' chosen frequencies.')") cfg_rsp_nfreq
     !  !kk and ajt: default allocation of zero cfg_rsp_imfreq
     !   if (.not.cfg_rsp_imfreqs_specified) then
     !      cfg_rsp_nimfreq = cfg_rsp_nfreq
     !      allocate(cfg_rsp_imfreq(cfg_rsp_nimfreq))
     !      cfg_rsp_imfreq = 0.0d0
     !   endif
     !else
     !   WRITE(LUPRI,*)
     !   WRITE(LUPRI, "('No frequencies specified in input. Only static properties will be &
     !               & calculated.')") 
     !endif
  endif

  if (config%decomp%cfg_rsp_nexcit > 0) then
    !if (cfg_rsp_dipole) then
    !   do i = 1,3  !run over dip X,Y,Z
    !     nops = nops + 1
    !     cfg_rsp_operators(1,nops) = i
    !   enddo
    !endif
    !!SONIANEW: chosen excited state must be within required solutions
    !if (cfg_rsp_run_exgr) then
    !   if (cfg_rsp_whichexst > config%decomp%cfg_rsp_nexcit) then
    !      WRITE(LUPRI,*) 'You requested an excited state which has not been determined'
    !      WRITE(LUPRI,*) 'Calculation stops. Change NEXCIT or WHICHEXST and come back'
    !      STOP '*** ERROR IN ESGRSP INPUT ( WHICHEXST > NEXCIT) '
    !   elseif (cfg_rsp_whichexst == 0 ) then
    !      WRITE(LUPRI,*) 'No excited state number specified. Use default value'
    !      WRITE(LUPRI,*) 'cfg_rsp_whichexst = config%decomp%cfg_rsp_nexcit (highest state)'
    !      cfg_rsp_whichexst = config%decomp%cfg_rsp_nexcit
    !   end if
    !endif
  else
    !if (cfg_rsp_dipole) then
    !   do j = 1,3  !run over dipole_X,Y,Z - see top of file
    !     do i = j,3  !run over dipole_X,Y,Z
    !       nops = nops + 1
    !       cfg_rsp_operators(1,nops) = i
    !       cfg_rsp_operators(2,nops) = j
    !     enddo
    !   enddo
    !endif
  endif
  !number of "vectors" = number of operator pair to evaluate
  !cfg_rsp_nvec = nops  

END SUBROUTINE config_rsp_input



SUBROUTINE READ_WAVE_DFTINPUT(LUPRI,LUCMD,DALTON,WORD)
implicit none
TYPE(integralconfig)   :: DALTON
INTEGER            :: LUCMD !Logical unit number for the daltoninput
character(len=40),intent(out)  :: WORD
character(len=1)   :: PROMPT
INTEGER            :: LUPRI
character(len=80)  :: LINE

DO
   READ (LUCMD, '(A40)') WORD
   PROMPT = WORD(1:1)
   IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
   IF (PROMPT .EQ. '.') THEN
      SELECT CASE(WORD) 
      !PLEASE KEEP IN ALPHABETICAL ORDER
      CASE ('.DFTELS'); READ(LUCMD,*) DALTON%DFTELS
      CASE ('.DFTTHR'); READ(LUCMD,*) DALTON%DFTHR0,DALTON%DFTHRL, DALTON%DFTHRI, DALTON%RHOTHR
      CASE ('.RADINT'); READ(LUCMD,*) DALTON%RADINT
      CASE ('.ANGINT'); 
         READ(LUCMD,*) DALTON%ANGINT
         IF(DALTON%ANGINT .GT. 64) CALL lsquit('ANGINT > 64 not implemented',-1)
      CASE ('.ANGMIN'); READ(LUCMD,*) DALTON%ANGMIN
      CASE ('.NOPRUN'); DALTON%NOPRUN = .TRUE.
      CASE ('.GRID TYPE'); 
         READ(LUCMD,'(A80)') LINE
         CALL DFTGRIDINPUT(LINE,DALTON%TURBO)
      CASE ('.DFTAC');
         CALL lsQUIT(' Not implemented in the new code - if you know what this&
              & you could tell me and I might implement it. TK',lupri)
         DALTON%DFTASC = .TRUE.
         DALTON%DFTPOT = .TRUE.
         READ (LUCMD,*) DALTON%DFTIPT, DALTON%DFTBR1, DALTON%DFTBR2
      CASE('.CARPAR');
         call lsquit('Cannot call DFTCARTESIANINPUT, only new integral code is compiled',-1)
      CASE ('.COARSE'); DALTON%RADINT = 1.D-11; DALTON%ANGINT = 35; 
      CASE ('.NORMAL'); DALTON%RADINT = 1.D-13; DALTON%ANGINT = 35; 
      CASE ('.FINE'  ); DALTON%RADINT = 1.D-13; DALTON%ANGINT = 42; 
      CASE ('.ULTRAF'); DALTON%RADINT = 1.D-15; DALTON%ANGINT = 64;
      CASE ('HARTRE'); 
         DALTON%DFTADD = .FALSE.
      CASE ('.GRID1' ) 
         CALL DFTGRIDINPUT("TURBO BLOCK",DALTON%TURBO);
         DALTON%RADINT = 1.D-5; DALTON%ANGINT = 17;
      CASE ('.GRID2' ) 
         CALL DFTGRIDINPUT("TURBO BLOCK",DALTON%TURBO);
         DALTON%RADINT = 2.15447D-7; DALTON%ANGINT = 23;
      CASE ('.GRID3' ) 
         CALL DFTGRIDINPUT("TURBO BLOCK",DALTON%TURBO);
         DALTON%RADINT = 4.64159D-9; DALTON%ANGINT = 29;
      CASE ('.GRID4' ) 
         CALL DFTGRIDINPUT("TURBO BLOCK",DALTON%TURBO);
         DALTON%RADINT = 5.01187D-14; DALTON%ANGINT = 35;
      CASE ('.GRID5' ) 
         CALL DFTGRIDINPUT("TURBO BLOCK",DALTON%TURBO);
         DALTON%RADINT = 2.15443D-17; DALTON%ANGINT = 47;
      CASE ('.AOSAVE' ) 
         CALL lsQUIT('Not implemented in the new code, because I am not&
                    &convinced about the performance, convince me and I will. TK',lupri);
      CASE ('.HARDNESS' ); READ(LUCMD,*) DALTON%HRDNES
      CASE ('.DISPER' )
         DALTON%DODISP = .TRUE.
         CALL DFTDISPCHECK()
      CASE DEFAULT
         WRITE (LUPRI,'(/,3A,/)') ' Keyword ',WORD,&
              & ' not recognized in *DFT INPUT'
         CALL lsQUIT('Illegal keyword in *DFT INPUT',lupri)
      END SELECT
   ENDIF
   IF (PROMPT .EQ. '*') EXIT
ENDDO

!WRITE (LUPRI,'(4X,A,17X,3F12.2)')&
!&        ' DFT LSint thresholds:', DALTON%DFTHR0, DALTON%DFTHRL, DALTON%DFTHRI 
!WRITE (LUPRI,'(4X,A,17X,3D12.2)')&
!&        ' DFT LSint threshold for number of electrons: ', DALTON%DFTELS 
!WRITE (LUPRI,'(4X,A,F8.4,I4)')' DFT LSint radial quadrature accuracy/ang. &
!     & expansion order:',DALTON%RADINT,DALTON%ANGINT
!
END SUBROUTINE READ_WAVE_DFTINPUT

SUBROUTINE READ_INTEGRALS_DENFIT_INPUT(LUPRI,LUCMD,DALTON,word)
implicit none
TYPE(integralconfig)   :: DALTON
INTEGER            :: LUCMD !Logical unit number for the daltoninput
character(len=40),intent(out)  :: WORD
character(len=1)   :: PROMPT
INTEGER            :: LUPRI

! the old version is SUBROUTINE DFIINP

DO
  READ (LUCMD, '(A40)') WORD
  PROMPT = WORD(1:1)
  IF (PROMPT(1:1) .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
  IF (PROMPT .EQ. '.') THEN
    SELECT CASE(WORD) 
      CASE ('.DIATOM'); WRITE(LUPRI,*)'.DIATOM NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.DOBOX '); WRITE(LUPRI,*)'.DOBOX NOT IMPLEMENTED IN NEW DRIVER'  
      CASE ('.BOXORB'); WRITE(LUPRI,*)'.BOXORB NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.BUFSIZ'); WRITE(LUPRI,*)'.BUFSIZ NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.CONSTR'); WRITE(LUPRI,*)'.CONSTR NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOROBU'); WRITE(LUPRI,*)'.NOROBU NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.ATTENU'); WRITE(LUPRI,*)'.ATTENU NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.GAUDAM'); WRITE(LUPRI,*)'.GAUDAM NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NONLIN'); WRITE(LUPRI,*)'.NONLIN NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.HERMIT'); WRITE(LUPRI,*)'.HERMIT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.3CENT '); WRITE(LUPRI,*)'.3CENT  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.OAOINT'); WRITE(LUPRI,*)'.OAOINT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NO_J_M'); WRITE(LUPRI,*)'.NO_J_M NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NO_K_M'); WRITE(LUPRI,*)'.NO_K_M NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.RI-LIN'); WRITE(LUPRI,*)'.RI-LIN NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NONRMA'); WRITE(LUPRI,*)'.NONRMA NOT IMPLEMENTED IN NEW DRIVER' 
      CASE DEFAULT
                  WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                  & '" not recognized in DENFIT INTEGRAL.'
                 CALL lsQUIT('Illegal keyword in DENFIT INTEGRAL.',lupri)
    END SELECT
  ENDIF
  IF (PROMPT .EQ. '*') EXIT
ENDDO
END SUBROUTINE READ_INTEGRALS_DENFIT_INPUT

SUBROUTINE READ_INTEGRALS_FMM_INPUT(LUPRI,LUCMD,DALTON,word)
implicit none
TYPE(integralconfig)   :: DALTON
INTEGER            :: LUCMD !Logical unit number for the daltoninput
character(len=40),intent(out)  :: WORD
character(len=1)   :: PROMPT
INTEGER            :: LUPRI

! the old version is SUBROUTINE DFIINP

DO
  READ (LUCMD, '(A40)') WORD
  PROMPT = WORD(1:1)
  IF (PROMPT(1:1) .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
  IF (PROMPT .EQ. '.') THEN
    SELECT CASE(WORD(1:7)) 
      CASE ('.SKIP  '); WRITE(LUPRI,*)'.SKIP   NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.PRINT '); WRITE(LUPRI,*)'.PRINT  NOT IMPLEMENTED IN NEW DRIVER'  
      CASE ('.STOP  '); WRITE(LUPRI,*)'.STOP   NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.LMAX  '); READ(LUCMD,*) DALTON%MM_LMAX
      CASE ('.NLEVEL'); WRITE(LUPRI,*)'.NLEVEL NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.FQUAD '); WRITE(LUPRI,*)'.FQUAD  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NNQUAD'); WRITE(LUPRI,*)'.NNQUAD NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.BQUAD '); WRITE(LUPRI,*)'.BQUAD  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NLOGN '); WRITE(LUPRI,*)'.NLOGN  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.FMM   '); WRITE(LUPRI,*)'.FMM    NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.GRSRCH'); WRITE(LUPRI,*)'.GRSRCH NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.TLMAX '); READ(LUCMD,*) DALTON%MM_TLMAX
      CASE ('.UMAT  '); WRITE(LUPRI,*)'.UMAT   NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.BRFREE'); WRITE(LUPRI,*)'.BRFREE NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.GRAIN '); WRITE(LUPRI,*)'.GRAIN  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.DYNLMA'); WRITE(LUPRI,*)'.DYNLMA NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.ALLSQR'); WRITE(LUPRI,*)'.ALLSQR NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.TRSRCH'); WRITE(LUPRI,*)'.TRSRCH NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOBOXP'); WRITE(LUPRI,*)'.NOBOXP NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.PRSTAT'); WRITE(LUPRI,*)'.PRSTAT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.TCON  '); WRITE(LUPRI,*)'.TCON   NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.SCREEN'); READ(LUCMD,*) DALTON%MM_SCREEN
         DALTON%MM_SCREEN = DALTON%MM_SCREEN/DALTON%THRESHOLD
      CASE ('.SKIPNN'); WRITE(LUPRI,*)'.SKIPNN NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.MMSAVE'); WRITE(LUPRI,*)'.MMSAVE NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOONE '); DALTON%MM_NO_ONE = .TRUE.
      CASE ('.NOMMBU'); DALTON%USEBUFMM  = .FALSE.
      CASE DEFAULT
                  WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                  & '" not recognized in FMM INTEGRAL.'
                 CALL lsQUIT('Illegal keyword in FMM INTEGRAL.',lupri)
    END SELECT
  ENDIF
  IF (PROMPT .EQ. '*') EXIT
ENDDO
END SUBROUTINE READ_INTEGRALS_FMM_INPUT

SUBROUTINE READ_INTEGRALS_FCK3_INPUT(LUPRI,LUCMD,DALTON,word)
implicit none
TYPE(integralconfig)   :: DALTON
INTEGER            :: LUCMD !Logical unit number for the daltoninput
character(len=40),intent(out)  :: WORD
character(len=1)   :: PROMPT
INTEGER            :: LUPRI

! the old version is SUBROUTINE FCK3INP

DO
  READ (LUCMD, '(A40)') WORD
  PROMPT = WORD(1:1)
  IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
  IF (PROMPT .EQ. '.') THEN
    SELECT CASE(WORD) 
      CASE ('.FIXTHR'); WRITE(LUPRI,*)'.FIXTHR NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.THRFCK'); WRITE(LUPRI,*)'.THRFCK NOT IMPLEMENTED IN NEW DRIVER'
         READ(LUCMD, '(A40)') WORD
         WRITE(LUPRI,*)'next word',WORD
      CASE ('.PRINT'); WRITE(LUPRI,*)'.PRINT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSPEX'); WRITE(LUPRI,*)'.NOSPEX NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSCRP'); WRITE(LUPRI,*)'.NOSCRP NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSKCL'); WRITE(LUPRI,*)'.NOSKCL NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSKCS'); WRITE(LUPRI,*)'.NOSKCS NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.SCRFBT'); WRITE(LUPRI,*)'.SCRFBT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOPACK'); WRITE(LUPRI,*)'.NOPACK NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSKAL'); WRITE(LUPRI,*)'.NOSKAL NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.FULL_J'); WRITE(LUPRI,*)'.FULL_J  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.UNCONT'); WRITE(LUPRI,*)'.UNCONT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.BRFREE'); WRITE(LUPRI,*)'.BRFREE NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.CHCKNA'); WRITE(LUPRI,*)'.CHCKNA NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NCSODB'); WRITE(LUPRI,*)'.NCSODB NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NCSTUV'); WRITE(LUPRI,*)'.NCSTUV NOT IMPLEMENTED IN NEW DRIVER' 

      CASE ('.TWOGRA'); WRITE(LUPRI,*)'.TWOGRA NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.HERMIT'); WRITE(LUPRI,*)'.HERMIT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.TWOSHI'); WRITE(LUPRI,*)'.TWOSHI NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.PRPSCR'); WRITE(LUPRI,*)'.PRPSCR NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.PRPSTH'); WRITE(LUPRI,*)'.PRPSTH NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NUCATR'); WRITE(LUPRI,*)'.NUCATR NOT IMPLEMENTED IN NEW DRIVER'
      CASE DEFAULT
                  WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                  & '" not recognized in *FCK INTEGRAL.'
                 CALL lsQUIT('Illegal keyword in FCK3 INTEGRAL.',lupri)
    END SELECT
  ENDIF
  IF (PROMPT .EQ. '*') EXIT
ENDDO

END SUBROUTINE READ_INTEGRALS_FCK3_INPUT

!> \brief Check if keywords conform and print configuration.
!> \author S. Host
!> \date March 2010
!>
!> If keywords specified in DALTON.INP do not conform, there are two options: \n
!> 1. Clean up, i.e. change the settings specified by the user to something
!>    meaningful. Remember to clarify this in output! E.g. \n
!>    'H1DIAG does not work well with only few saved microvectors for ARH. 
!>    Resetting max. size of subspace in ARH linear equations to', <something meaningful>. \n
!> 2. Quit, if there is no logical way to recover. \n
!> After deciding what the final configuration should be, print selected details
!> which may be useful for the user.
!>
subroutine set_final_config_and_print(lupri,config,ls)
use scf_stats
implicit none
   !> Contains info, settings and data for entire calculation
   type(configItem), intent(inout) :: config
   !> Logical unit number for LSDALTON.OUT
   integer, intent(in)             :: lupri
   !> Object containing integral settings and molecule
   type(lsitem), intent(inout)     :: ls
!
   integer                         :: i
   logical                         :: file_exists
   real(realk)                     :: conv_factor, potnuc
   CHARACTER*24, PARAMETER :: AVG_NAMES(5) = &
        &  (/ 'None                    ', &
        &     'DSM                     ', &
        &     'DIIS                    ', &
        &     'EDIIS                   ', &
        &     'Van Lenthe modified DIIS' /)
   CHARACTER*49, PARAMETER :: dsm_names(4) = &
        &  (/ 'Standard DSM                                     ', &
        &     'Only one iteration in DSM                        ', &
        &     'Line search in the steplength after one iteration', &
        &     'Extra accurate DSM energy model                  '/) 
   CHARACTER*35, PARAMETER :: F2D_NAMES(4) = (/ &
        & 'Diagonalization            ',&
        & 'Direct density optimization',&
        & 'Purification scheme        ',&
        & 'Augmented RH optimization  ' /)
   CHARACTER*35, PARAMETER :: SHIFT_NAMES(5) = (/ &
        & 'Level shifting by MO overlap       ',&
        & 'Level shifting by line search in mu',&
        & 'Level shifting by ||Dorth|| ratio  ',&
        & 'No level shifting                  ',&
        & 'Van Lenthe fixed level shifts      '/)

   write(config%lupri,*) 'Configuration - as obtained from new input structure:'

   if (config%opt%cfg_prefer_BSM) then
#if !defined(HAVE_BSM)
      CALL lsQUIT('.BLOCK requested but BSM not there',config%lupri)
#endif
   endif

   config%av%dsm_history_size   = config%av%cfg_settings(config%av%CFG_SET_type)%max_history_size
   config%av%diis_history_size  = config%av%dsm_history_size
   config%av%ediis_history_size = config%av%dsm_history_size

!Printing the configuration for the calculation:
!===============================================

   if(ls%input%DO_DFT)CALL DFTREPORT(lupri) !print the functional

   WRITE (config%LUPRI,'(4X,A,17X,3D12.2)')&
   &        ' DFT LSint thresholds:', config%integral%DFTHR0, config%integral%DFTHRL, config%integral%DFTHRI 
   WRITE (config%LUPRI,'(4X,A,17X,3D12.2)')&
   &        ' DFT LSint threshold for number of electrons: ', config%integral%DFTELS 
   WRITE (config%LUPRI,'(4X,A,D12.4,I4)')' DFT LSint radial quadrature accuracy/ang. &
        & expansion order:',config%integral%RADINT,config%integral%ANGINT

   WRITE(config%LUPRI,*)
   if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
      write (config%lupri,*) '  You have requested Augmented Roothaan-Hall optimization'
      write (config%lupri,*) '  => explicit averaging is turned off!'
      config%av%cfg_averaging = config%av%cfg_avg_none
      !if (config%solver%debug_arh) then
      !   write (config%lupri,*) '  Warning: Debugging Augmented Roothaan-Hall'
      !   write (config%lupri,*) '  => ARH terms turned off!'
      !   config%solver%cfg_arhterms = .false.
      !endif
      write(config%lupri,"('Expand trust radius if ratio is larger than:          ',F6.2)") config%solver%cfg_arh_expand_crit
      write(config%lupri,"('Contract trust radius if ratio is smaller than:       ',F6.2)") config%solver%cfg_arh_contract_crit
      write(config%lupri,"('On expansion, trust radius is expanded by a factor    ',F6.2)") config%solver%cfg_arh_expand
      write(config%lupri,"('On contraction, trust radius is contracted by a factor',F6.2)") config%solver%cfg_arh_contract
      WRITE(config%LUPRI,*)
      if (config%solver%cfg_arh_truncate) then
         if (config%opt%cfg_start_guess == 'H1DIAG' .and. config%solver%cfg_arh_microvecs < 5) then
            WRITE(config%LUPRI,*) 'H1DIAG does not work well with only few saved microvectors for ARH...'
            config%solver%cfg_arh_microvecs = 5
            WRITE(config%LUPRI,*) 'Resetting max. size of subspace in ARH linear equations to:', config%solver%cfg_arh_microvecs
         else        
            WRITE(config%LUPRI,*) 'Maximum size of subspace in ARH linear equations:', config%solver%cfg_arh_microvecs
         endif
         WRITE(config%LUPRI,*)
      endif
   endif

   if (config%solver%cfg_2nd_order_all) then
      config%av%CFG_averaging = config%av%CFG_AVG_none 
      write (config%lupri,*) 'You have requested 2nd order optimization => no averaging (DIIS or DSM)!'
      config%solver%set_do_2nd_order = .true.
      WRITE(config%LUPRI,*)
   endif

   WRITE(config%LUPRI,"('LINSCF configuration :')")
   WRITE(config%LUPRI,"('Density subspace min. method    : ',A)") AVG_NAMES(config%av%CFG_averaging)
   if (config%av%cfg_averaging == config%av%cfg_avg_dsm) then
     WRITE(config%LUPRI,"('  dsm approach: ',A)") dsm_names(config%av%cfg_dsm_app)
   endif

   WRITE(config%LUPRI,"('Density optimization : ',A)") F2D_NAMES(config%opt%CFG_density_method)
   !if (config%opt%CFG_density_method == config%opt%CFG_F2D_ROOTHAAN) then
   !   cfg_rsp_mostart = .true.
   !   cfg_rsp_moprec = .true.
   !endif

!Settings concerning diagonalization and averaging:
!==================================================

   if (config%diag%cfg_lshift .ne. config%diag%cfg_lshift_none) then
     WRITE(config%LUPRI,"('with ',A)") SHIFT_NAMES(config%diag%CFG_lshift)
     if (config%diag%cfg_lshift == config%diag%cfg_lshift_Mochange) then
       WRITE(config%LUPRI,"('  where ',f6.3,' as smallest accepted overlap  grep')") &
                  & config%av%cfg_settings(config%av%CFG_SET_TYPE)%min_density_overlap
     elseif (config%diag%cfg_lshift ==config%diag%cfg_lshift_dorth) then
       WRITE(config%LUPRI,"('  where ',f6.3,' is largest accepted Dorth ratio  grep')") &
                   & config%av%cfg_settings(config%av%CFG_SET_TYPE)%max_dorth_ratio
     endif
   endif
   WRITE(config%LUPRI,*)
   if (config%diag%cfg_custom_shift) then
      write(config%lupri,"(i3,' levelshifts with the values')") config%diag%cfg_nshifts
      do i = 1, config%diag%cfg_nshifts
         write(config%lupri,"(f5.2)") config%diag%cfg_levelshifts(i)
      enddo
      write(config%lupri, "('have been requested.')")
   endif

   WRITE(config%LUPRI,*)
   !find the maximum number of stored vectors
   config%av%cfg_settings%max_history_size = MAX(config%av%dsm_history_size,config%av%diis_history_size)

   WRITE(config%LUPRI,*) 'Maximum size of Fock/density queue in averaging:', &
      &  config%av%cfg_settings(config%av%CFG_SET_type)%max_history_size

   if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then
      ALLOCATE(config%diag%cfg_levelshifts(100))

      config%diag%cfg_fixed_shift = .true.
      config%diag%cfg_custom_shift = .true.

      if (config%opt%calctype == config%opt%dftcalc) then
         config%diag%cfg_levelshifts = 1.0d0
      else
         config%diag%cfg_levelshifts = 0.3d0
         do i = 1, 5
            config%diag%cfg_levelshifts(i) = 1.0d0
         enddo
      endif
   endif

!Settings for unrestricted:
!==========================

   IF(MOD(config%integral%nelectrons,2) == 0)THEN
      !Even number of electrons
!      config%decomp%nocc = (config%integral%nelectrons - config%integral%molcharge)/2
      !Cecilie change 07/07 2010: Why subtract molcharge?
      ! Shouldn't it just be nelectrons / 2 ??
      ! This way it gives nocc = 3 for four electrons in H2
      config%decomp%nocc = config%integral%nelectrons/2
      config%decomp%nactive = 0

   ELSE
      !Odd number of electrons
      !Stinne change 23/4-2010: why subtract one here???
      !config%decomp%nocc = (config%integral%nelectrons - 1 - config%integral%molcharge)/2
!      config%decomp%nocc = (config%integral%nelectrons - config%integral%molcharge)/2
      !Cecilie change 07/07 2010: Same here
      config%decomp%nocc = config%integral%nelectrons/2
      config%decomp%nactive = 1
   ENDIF
   config%diag%nocc = config%decomp%nocc

   if (config%decomp%alpha_specified .or. config%decomp%beta_specified) then
      config%integral%unres =.TRUE.
      config%decomp%cfg_unres =.TRUE.
      config%diag%cfg_unres =.TRUE.
      config%opt%cfg_unres =.TRUE.
      !write(lupri,*) 'alpha_specified, beta_specified', alpha_specified, beta_specified
      if (config%decomp%alpha_specified .and. config%decomp%beta_specified) then
         if (config%decomp%nocca + config%decomp%noccb /= 2*config%decomp%nocc + config%decomp%nactive) then
            call lsquit('Nalpha + Nbeta differs from number of electrons!',config%lupri)
         endif
      else if (config%decomp%alpha_specified) then
         config%decomp%NOCCB = 2*config%decomp%NOCC + config%decomp%nactive - config%decomp%NOCCA
      else if (config%decomp%beta_specified) then
         config%decomp%NOCCA = 2*config%decomp%NOCC + config%decomp%nactive - config%decomp%NOCCB
      endif 
      config%integral%nelectrons = config%decomp%NOCCA + config%decomp%NOCCB
      write(LUPRI,'(/,1x,a)') '--------------------------'
      write(LUPRI,'(1x,a)')   '<Unrestricted calculation>'
      write(LUPRI,'(1x,a)')   '--------------------------'
      write(LUPRI,'(1x,a,i6)')   'ALPHA spin occupancy =',config%decomp%nocca
      write(LUPRI,'(1x,a,i6,/)') 'BETA  spin occupancy =',config%decomp%noccb
      call mat_select_type(mtype_unres_dense)
   else IF(config%decomp%nactive /= 0 .or. config%decomp%cfg_unres) THEN
      !unrestricted SCF if Nelec uneven or if cfg_unres=.true.

      config%integral%unres = .true.
      config%decomp%cfg_unres = .true.
      config%diag%cfg_unres = .true.
      config%opt%cfg_unres = .true.

      config%decomp%NOCCA = config%decomp%NOCC
      config%decomp%NOCCB = config%decomp%NOCC + config%decomp%nactive

      config%diag%nocca = config%decomp%NOCCA
      config%diag%noccb = config%decomp%NOCCb
      if(config%integral%nelectrons /= 0) then
write(config%lupri,*) 'WARNING WARNING WARNING spin check commented out!!! /Stinne'

!FIXME: What should be here?? What does spin mean here??? What do we actually
!want to check???
         !if (config%decomp%spin == 0) then
         !   if (MOD(config%integral%nelectrons, 2) /= 0) then
         !      write(config%lupri,*) 'Spin:', config%decomp%spin
         !      write(config%lupri,*) 'Nelectrons:', config%integral%nelectrons
         !      write(config%lupri,*) 'MOD(Nelectrons,2):', MOD(config%integral%nelectrons, 2)
         !      call lsquit('Spin and number of electrons do not conform!')
         !   endif
         !   config%decomp%NOCCA = config%integral%nelectrons/2
         !   config%decomp%NOCCB = config%integral%nelectrons/2
         !else if (config%decomp%spin == 1) then
         !   if (MOD(config%integral%nelectrons, 2) /= 1) then
         !      write(config%lupri,*) 'Spin:', config%decomp%spin
         !      write(config%lupri,*) 'Nelectrons:', config%integral%nelectrons
         !      write(config%lupri,*) 'MOD(Nelectrons,2):', MOD(config%integral%nelectrons, 2)
         !      call lsquit('Spin and number of electrons do not conform!') 
         !   endif
         !   config%decomp%NOCCA = config%integral%nelectrons/2
         !   config%decomp%NOCCB = config%integral%nelectrons - config%decomp%NOCCA
         !else if(config%decomp%spin == 2) then
         !   if (MOD(config%integral%nelectrons,2) /= 0) then
         !      write(config%lupri,*) 'Spin:', config%decomp%spin
         !      write(config%lupri,*) 'Nelectrons:', config%integral%nelectrons
         !      write(config%lupri,*) 'MOD(Nelectrons,2):', MOD(config%integral%nelectrons, 2)
         !      call lsquit('Spin and number of electrons do not conform!') 
         !   endif
         !   config%decomp%NOCCA = config%integral%nelectrons/2 - 1
         !   config%decomp%NOCCB = config%integral%nelectrons - config%decomp%NOCCA
         !else
         !   call lsquit('Only spin = 0, 1, and 2 implemented!') 
         !endif
      else
         config%decomp%NOCCA = config%integral%nelectrons/2 
         config%decomp%NOCCB = config%integral%nelectrons/2
         if (config%decomp%spin == 0) then
            !occupations as above  
         else if (config%decomp%spin == 1) then
            if(MOD(config%integral%nelectrons, 2) /= 1) &
                 & call lsquit('Spin and Nelec do not conform!',config%lupri) 
            config%decomp%NOCCA = config%integral%nelectrons/2
            config%decomp%NOCCB = config%integral%nelectrons - config%decomp%NOCCA
         else if(config%decomp%spin == 2) then
            if (MOD(config%integral%nelectrons, 2) /= 0) &
                 & call lsquit('Spin and Nelec do not conform!',config%lupri) 
            config%decomp%NOCCA = config%integral%nelectrons/2 - 1
            config%decomp%NOCCB = config%integral%nelectrons - config%decomp%NOCCA
         else
            call lsquit('Only spin = 0, 1, and 2 implemented!',config%lupri) 
         endif
      endif
      write(LUPRI,'(/,1x,a)') '--------------------------'
      write(LUPRI,'(1x,a)')   '<Unrestricted calculation>'
      write(LUPRI,'(1x,a)')   '--------------------------'
      if(config%decomp%spin == 2) &
           & write(LUPRI,'(1x,a)') 'Spin symmetry = Triplet'
      write(LUPRI,'(1x,a,i6)')   'ALPHA spin occupancy =', config%decomp%nocca
      write(LUPRI,'(1x,a,i6,/)') 'BETA  spin occupancy =', config%decomp%noccb
      !fixme: should be available for other matrix types as well
      call mat_select_type(mtype_unres_dense)
   ENDIF

!Settings concerning SCF gradient convergence threshold:
!=======================================================

   if (config%opt%cfg_convdyn) then
      !Until someone writes a better XC grid generation code, we have to use different thresholds for HF and DFT
      if (config%solver%do_dft) then
         SELECT CASE(config%opt%cfg_convdyn_type)
         CASE('TIGHT'); conv_factor = 1.0d-5
         CASE('STAND'); conv_factor = 1.0d-4
         CASE('SLOPP'); conv_factor = 1.0d-3
         CASE DEFAULT
         WRITE (LUPRI,'(/,3A,/)') ' Keyword "',config%opt%cfg_convdyn_type,&
              & '" not recognized with .CONVDYN'
         WRITE (LUPRI,*) 'Options are TIGHT, STANDARD, and SLOPPY.'
         CALL lsQUIT('Illegal keyword with .CONVDYN',config%lupri)
         END SELECT
      else
         SELECT CASE(config%opt%cfg_convdyn_type)
         CASE('TIGHT'); conv_factor = 1.0d-6
         CASE('STAND'); conv_factor = 1.0d-5
         CASE('SLOPP'); conv_factor = 1.0d-4
         CASE DEFAULT
         WRITE (LUPRI,'(/,3A,/)') ' Keyword "',config%opt%cfg_convdyn_type,&
              & '" not recognized with .CONVDYN'
         WRITE (LUPRI,*) 'Options are TIGHT, STANDARD, and SLOPPY.'
         CALL lsQUIT('Illegal keyword with .CONVDYN',config%lupri)
         END SELECT
      endif

      if (config%decomp%cfg_unres) then
         config%opt%cfg_convergence_threshold = conv_factor*sqrt((config%decomp%nocca+config%decomp%noccb)*1.0d0)
         config%opt%set_convergence_threshold = conv_factor*sqrt((config%decomp%nocca+config%decomp%noccb)*1.0d0)
      else
         config%opt%cfg_convergence_threshold = conv_factor*sqrt(config%decomp%nocc*2.0d0)
         config%opt%set_convergence_threshold = conv_factor*sqrt(config%decomp%nocc*2.0d0)
      endif

      !WRITE(lupri,*)
      WRITE(config%LUPRI,"('Dynamic convergence threshold for gradient: ',E10.2)") &
           & config%opt%set_convergence_threshold
      WRITE(config%lupri,*)
   else
      WRITE(config%LUPRI,"('Convergence threshold for gradient: ',E10.2)") &
           & config%opt%set_convergence_threshold
   endif

!Settings for HOMO-LUMO gap, Hessian eigenvalues and rsp starting guess:
!=======================================================================

   if (config%decomp%cfg_startvectors) then
      if (config%decomp%cfg_no_of_startvectors < config%decomp%cfg_rsp_nexcit) then
         WRITE(config%LUPRI,*)
         WRITE(config%LUPRI,"('WARNING: Number of start vectors smaller than number of excitation energies!')")
         WRITE(config%LUPRI,"('Resetting .NSTART to match .NEXCI')")
         config%decomp%cfg_no_of_startvectors = config%decomp%cfg_rsp_nexcit
      else if (config%decomp%cfg_no_of_startvectors < config%decomp%cfg_hessian_nvec) then
         WRITE(config%LUPRI,*)
         WRITE(config%LUPRI,"('WARNING: Number of start vectors smaller than number of requested Hessian eigenvalues!')")
         WRITE(config%LUPRI,"('Resetting .NSTART to match .HESVEC')")
         config%decomp%cfg_no_of_startvectors = config%decomp%cfg_hessian_nvec
      endif
      if (config%decomp%cfg_homolumo_maxit < 100*config%decomp%cfg_no_of_startvectors) then
          WRITE(config%LUPRI,"('Increasing max. no. of HOMO-LUMO iterations to ', i6)") 100*config%decomp%cfg_no_of_startvectors
          config%decomp%cfg_homolumo_maxit = 100*config%decomp%cfg_no_of_startvectors
      endif
      if (config%decomp%cfg_check_maxit < 50*config%decomp%cfg_hessian_nvec) then
          WRITE(config%LUPRI,"('Increasing max. no. of Hessian iterations to ', i6)") 50*config%decomp%cfg_hessian_nvec
          config%decomp%cfg_check_maxit = 50*config%decomp%cfg_hessian_nvec
      endif
   else if (config%decomp%cfg_rsp_nexcit > 0) then
      if (config%decomp%cfg_homolumo_maxit < 100*config%decomp%cfg_rsp_nexcit) then
          WRITE(config%LUPRI,"('Increasing max. no. of HOMO-LUMO iterations to ', i6)") 100*config%decomp%cfg_rsp_nexcit
          config%decomp%cfg_homolumo_maxit = 100*config%decomp%cfg_rsp_nexcit
      endif
   else if (config%decomp%cfg_hessian_nvec > 1) then
      if (config%decomp%cfg_check_maxit < 50*config%decomp%cfg_hessian_nvec) then
          WRITE(config%LUPRI,"('Increasing max. no. of Hessian iterations to ', i6)") 50*config%decomp%cfg_hessian_nvec
          config%decomp%cfg_check_maxit = 50*config%decomp%cfg_hessian_nvec
      endif
   endif

   !Check if HOMO-LUMO is really necessary for further calculations (Hessian eival or exci energies)
   ! - if not, we can print a warning instead of quitting if convergence fails
   if (config%decomp%cfg_check_converged_solution .or. config%decomp%cfg_rsp_nexcit > 0) then
      config%decomp%cfg_hlgap_needed = .true.
   endif

! Check integral input:
!======================

   if(config%integral%densfit .AND. (.NOT. config%integral%auxbasis))then
      WRITE(config%LUPRI,'(/A)') &
           &     'You have specified .DENSFIT in the dalton input but not supplied a fitting basis set'
      CALL lsQUIT('Density fitting input inconsitensy: add fitting basis set',config%lupri)
   endif

   if(config%integral%DALINK .AND. config%opt%cfg_incremental)THEN
      WRITE(config%lupri,*)'DalinK does currently not work with the incremental Fock matrix scheme'
      WRITE(config%lupri,*)'either use .LINK or .NONCREM to turn off incremental Fock matrix scheme'
      print*,'DalinK does currently not work with the incremental Fock matrix scheme'
      print*,'either use .LINK or .NONCREM to turn off incremental Fock matrix scheme'
      CALL lsQUIT('DalinK does currently not work with the incremental Fock matrix scheme',config%lupri)
   endif

   if((.NOT.config%integral%LINK).AND.(.NOT.config%integral%LINK) )THEN
      WRITE(config%lupri,*)'You have chosen to run without the use of Jengine or LinK'
      WRITE(config%lupri,*)'We therefor deactivate the incremental Fock matrix scheme'
      WRITE(config%lupri,*)'as it will provide no speedup - only reduced accuracy'
      config%opt%cfg_incremental = .FALSE.
   endif
   
   !Note that the placement of this if statement is important, as the config%opt%cfg_incremental is subject to change
   if(config%opt%cfg_incremental)THEN
      IF(ABS(1.0D-10 - config%integral%threshold).LT.1.0D-12)THEN
         !incremental and no change in threshold
         WRITE(config%lupri,*)'You have chosen to run with the incremental Fock matrix scheme (default using Jengine or LinK)'
         WRITE(config%lupri,*)'Due to the accumulated error in the incremental scheme we thighten the integral threshold'
         WRITE(config%lupri,*)'This new threshold is the same that you will obtain by adding'
         WRITE(config%lupri,*)'.THRESH'
         WRITE(config%lupri,*)'1.D-11'
         WRITE(config%lupri,*)'To the DALTON.INP file'
         config%integral%threshold = 1.0D-11
      ELSE
         !incremental, but manual change in threshold
         WRITE(config%lupri,*)'You have chosen to run with the incremental Fock matrix scheme (default using Jengine or LinK)'
         WRITE(config%lupri,*)'Due to the accumulated error in the incremental scheme we would usually thighten the'
         WRITE(config%lupri,*)'integral threshold to 1.D-11, but since you have specified a non default threshold we assume'
         WRITE(config%lupri,*)'that you know what you are doing' 
      ENDIF
   endif

!Check setting for linear equations iterative solver:
!====================================================

  config%opt%do_trustregion = &
      & (config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
      & (config%solver%cfg_2nd_order_all .or. config%solver%cfg_2nd_order_local))

   !Currently turned off, get_from_modFIFO_disk won't work!
   !if (config%solver%cfg_arh_disk_macro .and. config%opt%cfg_queue_on_disk) then
   !   write(config%lupri,*) 'It does not make sense to use .DUMP (dump previous'
   !   write(config%lupri,*) 'F, D matrices to disk while constructing new F), since'
   !   write(config%lupri,*) 'you have also requested that these should always'
   !   write(config%lupri,*) 'be stored on disk! Ignoring .DUMP !!'
   !   config%opt%cfg_queue_on_disk = .false.
   !endif

   if (config%opt%cfg_queue_on_disk) then
      write(config%lupri,*) 'Dump previous F, D matrices to disk while constructing new F'
      write(config%lupri,*) '- memory saving, but time consuming!'
   endif

   !if (config%opt%cfg_arh_disk_macro) then !Currently turned off, get_from_modFIFO_disk won't work!
   !   write(config%lupri,*) 'Keep queue of previous F, D matrices on disk instead of in core'
   !   write(config%lupri,*) '- memory saving, but time consuming!'
   !endif

   if (config%solver%cfg_arh_disk_micro) then
      write(config%lupri,*) 'ARH solver: Keep trial and sigma vectors on disk instead of in core'
      write(config%lupri,*) '- memory saving, but time consuming!'
   endif

   if (config%solver%cfg_arh_newdamp .and. .not. config%solver%cfg_arh_truncate) then
      WRITE(config%LUPRI,'(/A)') &
      &     'Can only use new damping scheme with .TRUNCATE'
       CALL lsQUIT('Can only use new damping scheme with .TRUNCATE',config%lupri)
   endif 

   if (config%opt%CFG_density_method == config%opt%CFG_F2D_DIRECT_DENS) then
      config%solver%cfg_arhterms = .false.
      config%solver%set_arhterms = .false.
      config%solver%cfg_arh_crop = .true.
      config%solver%cfg_arh_truncate = .false.
   endif

!Check for stuff not implemented for unrestricted:
!=================================================

   if (config%opt%cfg_start_guess=='TRILEVEL' .and. config%decomp%cfg_unres) then
      call lsquit('Sorry, trilevel starting guess not implemented for unrestricted!',config%lupri)
   endif

   if (config%solver%debug_dd) then
      if (config%decomp%cfg_unres) call lsquit('DEBUG_DD not implemented for unrestricted',config%lupri)
      call scf_stats_arh_header(config%lupri)
   endif

!Grand canonical basis stuff:
!============================

   if (config%opt%cfg_start_guess == 'ATOMS' .or. config%opt%cfg_start_guess == 'TRILEVEL') then
      config%decomp%cfg_gcbasis = .true.
   endif

!MKL sanity check:
!==================

   if (config%opt%cfg_prefer_CSR) then
      if (matrix_type == mtype_unres_dense) then
         call lsquit('Compressed Sparse Row (CSR) not implemented for unrestricted!',config%lupri)
      else
#ifdef VAR_MKL
         CALL mat_select_type(mtype_csr)
#else
         call lsquit('.CSR requires MKL library and -DVAR_MKL precompiler flag',config%lupri)
#endif
      endif
   endif

   if (matrix_type == mtype_csr) then
      if (config%opt%CFG_density_method == config%opt%CFG_F2D_roothaan) then
         write(config%lupri,*)
         write(config%lupri,*) 'The combination of diagonalization and CSR is very inefficient!'
         write(config%lupri,*) 'Please choose another type of density optimization' 
         write(config%lupri,*) '(e.g. .ARH, .TrFD, 2ND_ALL, 2ND_LOC) and welcome back...' 
         call lsquit('Combining diagonalization and CSR is inefficient!',-1)
      endif
   endif

!Local LINK check:
!=================

   write(config%lupri,*) 'End of configuration from new input structure!'

   ls%input%dalton%unres = config%decomp%cfg_unres

   CALL II_get_nucpot(lupri,lupri,ls%setting,POTNUC)
   config%opt%potnuc = POTNUC
   ls%input%potnuc = POTNUC

end subroutine set_final_config_and_print

subroutine init_integralconfig(integral,lupri)
implicit none
TYPE(integralconfig):: integral
INTEGER             :: LUCMD !Logical unit number for the daltoninput
INTEGER             :: IDUMMY,LUPRI,IPOS,IPOS2,COUNTER
character(len=70)   :: WORD
character(len=2)    :: PROMPT
LOGICAL             :: DONE,file_exists,READWORD,LSDALTON,STARTGUESS
character(len=16)   ::  start_guess

call integral_set_default_config(integral)

STARTGUESS = .FALSE.
integral%cfg_lsdalton = .TRUE.
COUNTER = 0

INQUIRE(file='DALTON.INP',EXIST=file_exists) 
IF(file_exists)THEN
   LUCMD=-1
   CALL lsOPEN(LUCMD,'DALTON.INP','OLD','FORMATTED')
ELSE
   CALL lsQUIT('DALTON.INP does not exist',lupri)
ENDIF
READWORD=.TRUE.
DONE=.FALSE.
rewind(LUCMD)
DO
   IF(DONE)EXIT
   IF(READWORD) THEN
      READ (LUCMD, '(A40)') WORD
      READWORD=.TRUE.
      COUNTER = 0
   ELSE
      IF (COUNTER.GT.1) THEN
        WRITE(LUPRI,'(1X,2A)') 'Infinite loop for input line:',WORD
        CALL lsQUIT('Infinite loop in init_integralconfig, due to input error',lupri)
      ENDIF
      COUNTER = COUNTER + 1
   ENDIF
   PROMPT = WORD(1:2)
   IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#'))CYCLE
   IF (WORD(1:10) == '**INTEGRAL') THEN
      READWORD = .TRUE.
      CALL INTEGRAL_INPUT(integral,readword,word,lucmd,lupri)
   ENDIF
   IF ((WORD(1:10) == '**WAVE FUN').OR.(WORD(1:10) == '**WAVEFUNC')) THEN
      READWORD=.TRUE.
      DO   
         IF(READWORD) THEN
            READ (LUCMD, '(A40)') WORD
            READWORD=.TRUE.
         ENDIF
         PROMPT = WORD(1:2)
         IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
         IF(PROMPT(1:1) .EQ. '.') EXIT
      ENDDO
   ENDIF
   IF (WORD(1:7) == '*LINSCA') THEN
      READWORD=.TRUE.
      DO   
         IF(READWORD) THEN
            READ (LUCMD, '(A40)') WORD
            READWORD=.TRUE.
         ENDIF
         PROMPT = WORD(1:2)
         IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
         IF(PROMPT(1:1) .EQ. '.') THEN
            SELECT CASE(WORD) 
            CASE('.GCBASIS');    integral%NOSEGMENT = .TRUE.; 
               integral%TRILEVEL = .TRUE.;
               integral%NOFAMILY = .TRUE.
               !all this is now default
            CASE('.START');      READ(LUCMD,*) start_guess 
                                 STARTGUESS = .TRUE.
                                 if (start_guess(1:8) == 'TRILEVEL' .or. &
                                 &   start_guess(1:5) == 'ATOMS') then
                                    integral%nosegment = .TRUE.
                                    integral%NOFAMILY = .TRUE.
                                    integral%trilevel = .TRUE.
                                    !all this is now default
                                 else !h1diag
                                    integral%nosegment = .FALSE.
                                    integral%NOFAMILY = .FALSE.
                                    integral%trilevel = .FALSE.                                    
                                 endif
            CASE DEFAULT
               !do nothing
            END SELECT
         ELSE IF (PROMPT(1:1) .EQ. '$') THEN
            CYCLE
         ENDIF
         IF(PROMPT .EQ. '**') THEN
            DONE=.TRUE.
            EXIT
            READWORD=.FALSE.
         ENDIF
         IF(PROMPT(1:1) .EQ. '*') THEN
            EXIT
            READWORD=.FALSE.
         ENDIF
      ENDDO
   ENDIF
   IF (WORD == '*END OF INPUT') THEN
      DONE=.TRUE.
   ENDIF
ENDDO

CALL lsCLOSE(LUCMD,'KEEP')

END SUBROUTINE init_integralconfig

end module configuration
