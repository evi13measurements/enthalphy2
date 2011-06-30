!> @file 
!> Contains response driver.

MODULE lsdalton_rsp_mod
!> \brief Driver for stand-alone f90 linear scaling SCF.
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2008-10-26
!>
!> PROPERTIES SECTION
!>
use response_wrapper_module
contains
  SUBROUTINE lsdalton_response(lu_pri,ls,config,F,D,S,dodft)
    use configuration
    use lstiming
    use prop_contribs, only: prop_oneave
    use rspsolver,     only: prop_molcfg
    use rsp_equations, only: rsp_eq_sol_empty
    use rsp_util,      only: util_save_MOinfo,util_free_MOstuff
    use TYPEDEF,       only: LSSETTING
    use molecule_type, only: MOLECULE_PT,MOLECULEINFO
    use matrix_defop,  only: operator(*)
    implicit none
    integer,intent(in)      :: lu_pri
    TYPE(lsitem),target     :: ls
    type(configItem),target :: config
    logical                 :: dodft
    real(realk)             :: Tstart,Tend,t1,t2 !,ten,tstr,E,gradnrm
    real(realk)             :: DipoleMoment(3)
    TYPE(Matrix)            :: F,D,S
    type(prop_molcfg)       :: molcfg
! Molecular gradient
    real(realk), pointer   :: Grad(:,:)

    !ajt This call generates Cmo/orbe in rsp_util, which rsp_solver needs
    if(config%response%RSPSOLVERinput%rsp_mo_precond) then
       call util_save_MOinfo(F,S,config%decomp%nocc) 
    endif
    call CPU_TIME(tstart)
    call LSTIMER('START',t1,t2,LU_PRI)

    !create config struct to be passed to prop_contribs / rsp_equations
    molcfg = prop_molcfg(0d0*S,ls%setting%MOLECULE(1)%p%Natoms, &
         & config%decomp%lupri,config%decomp%luerr, &
         & ls%setting,config%decomp,config%response%rspsolverinput)

    ! Kasper K, we ALWAYS calculate the permanent electric dipole for LSDALTON -
    ! also if no response properties have been requested.
    call Get_dipole_moment(molcfg,F,D,S,.true.,DipoleMoment)

    ! Determine excited states and one-photon absorption
    if(config%response%tasks%doOPA)then
       ! Determine transition density matrices and corresponding excitation energies.
       ! The transition density matrices and corresponding excitation energies
       ! are stored in solved_eqs(:)%D and solved_eqs(:)%freq(1) inside rsp_equations.
       ! IMPORTANT!!! This must be done BEFORE any response calculations involving excited states,
       ! i.e. residues of response functions etc.
       ! Therefore, this call should be the first thing in lsdalton_response!
       call calculate_and_store_transition_density_matrices(molcfg,F,D,S)
       ! Calculate one-photon absorption strengths for the excited states determined above
       call OPAresponse_driver(molcfg,F,D,S)
    endif

    ! KK: Place all response properties involving excited states here:
    ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! TPA
    if(config%response%tasks%doTPA)then
       call TPAresponse_driver(molcfg,F,D,S,config%response%tpainput)
    endif

    ! Excited state gradient
    if(config%response%tasks%doESG)then
       call ESGresponse_driver(molcfg,F,D,S,config%response%ESGinput)
    endif

    ! Excited state dipole moment
    if(config%response%tasks%doESD)then
       call ESDresponse_driver(molcfg,F,D,S,config%response%ESDinput,DipoleMoment)
    endif


    ! MCD, three-photon absorption etc.
    ! etc.

    ! Now all response properties involving excited states have been determined
    ! and we may delete the transition moment density matrices
    ! stored in solved_eqs(:)%D inside rsp_equations.
    ! IMPORTANT!!! This call has to be placed AFTER all response properties
    ! involving excited states!!!
    if(config%response%tasks%doOPA)then
       call free_transition_density_matrices(molcfg)
    endif

    ! Molecular gradient
    ! 29/10/11 Vladimir Rybkin: need Grad as input argument to 
    ! ii_get_molecular gradient elsewhere, therefore it's declared
    ! and allocated/deallocated here as well
    if(config%response%tasks%doGrad)then
       call mem_alloc(Grad,3,config%Molecule%NAtoms)
       call ii_get_molecular_gradient(Grad,lu_pri,F,D, &
       & ls%setting,dodft,.true.)
       call mem_dealloc(Grad)
    endif

    ! Polarizability
    if(config%response%tasks%doALPHA)then
       call ALPHAresponse_driver(molcfg,F,D,S,config%response%alphainput)
    endif

    ! 1st hyperpolarizability
    if(config%response%tasks%doBETA)then
       call BETAresponse_driver(molcfg,F,D,S,config%response%betainput,DipoleMoment)
    endif

    ! 2nd hyperpolarizability
    if(config%response%tasks%doGAMMA)then
       call GAMMAresponse_driver(molcfg,F,D,S,config%response%gammainput)
    endif

    ! Damped two-photon absorption
    if(config%response%tasks%doDTPA)then
       call DTPAresponse_driver(molcfg,F,D,S,config%response%dtpainput)
    endif

    call LSTIMER('LSDALTON RSP',t1,t2,LU_PRI)
    call CPU_TIME(tend)
    WRITE(lu_pri,*) "*****************************************************"
    Write(lu_pri,*) "**     CPU-TIME USED IN LSDALTON RESPONSE: ",tend-tstart,"   **"
    WRITE(lu_pri,*) "*****************************************************"

    ! Clear saved solutions of response equations, stored in
    call rsp_eq_sol_empty()  !rsp_eq_sol in module rsp_equations

    ! Free any allocated Cmo/orbe in rsp_util
    if(config%response%RSPSOLVERinput%rsp_mo_precond)then
       call util_free_MOstuff()
    endif

  END SUBROUTINE LSDALTON_RESPONSE

end MODULE lsdalton_rsp_mod
