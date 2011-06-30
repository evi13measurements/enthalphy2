#ifndef LSDALTON_ONLY
!=========================================!
! Main driver for linear scaling dynamics !
!=========================================!
Module Dynamics_driver
Contains
!===================!
! LS_dyn_run        !
!===================!
Subroutine LS_dyn_run(E,config,H1,F,D,S,ls,lupri,luerr)
Use Precision
Use ls_dynamics
Use configuration
Use matrix_module
Use lsdalton_fock_module
Use dalton_interface
!Use LSTiming
IMPLICIT NONE
Real(realk) :: E
Real(realk) :: TStep
Type(lsitem) :: ls
Type(Matrix), intent(inout) :: F,D,S
Type(Matrix), intent(in) :: H1
Type(ConfigItem), intent(inout) :: Config
Type(trajtype) :: Traj
Integer :: lupri, luerr
Integer :: NAtoms,err,i,j
!
Call Allocate_traj(Config%Molecule%nAtoms,Traj)
!
! Grabbing information from Config%Molecule to Trajectory
!
NAtoms = Config%Molecule%nAtoms
Write(*,*)'NATOMS',NAtoms
Write(*,*)'Mass',Config%Molecule%Atom(1)%Mass
Do i = 1,nAtoms
   Write(*,*)'i=',i
!   Write(*,*) Config%Molecule
   Traj%Mass(i) = Config%Molecule%Atom(i)%Mass
   ! Coordinates
   Do j = 1,3
      Traj%Coordinates(3*(i-1)+j)= &
      Config%Molecule%Atom(i)%Center(j) 
   Enddo   
Enddo
! Velocities
Traj%Velocities = Config%dynamics%Initial_Velocities
TStep = config%dynamics%TimeStep
!
Write(*,*)'Initial potential',E
Write(*,*)'Initial coordinates',Traj%Coordinates
Write(*,*)'Initial velocities',Traj%Velocities
!
! Get gradient
Traj%gradient = 0.D0
!
   Write(*,*)'Here!'
! Getting accelerations
Do i = 1,NAtoms
   Do j = 1,3
    Traj%Accel(3*(i-1)+j)=-Traj%Gradient(3*(i-1)+j)/Traj%Mass(i)
   Enddo
Enddo
!
! Taking first Verlet step
!
Traj%Velocities = &
Traj%Velocities + Traj%Accel*TStep*0.5D0
Traj%Coordinates = &
Traj%Coordinates + Traj%Velocities*TStep
! Updating congig
Call Pack_coordinates(config%Molecule,Traj%Coordinates,NAtoms)
! Reset setting for integral-code and update molecular coordinates
! Simen: ToDo - remove ls%input%Molecule and use only config%Molecule
call II_free_setting(ls%setting)
call II_init_setting(ls%setting)
Call Pack_coordinates(ls%input%Molecule,Traj%Coordinates,NAtoms)
call II_set_default_setting(ls%setting,ls%input)
! New energy
Write(*,*)'CALLING FOR NEW ENERGY!'
Write(*,*)'new coordinates'
Write(*,*)Config%Molecule%Atom(1)%Center
Write(*,*)'new ls-setting coordinates'
Write(*,*)ls%setting%Molecule(1)%p%Atom(1)%Center
Call II_get_overlap(lupri,luerr,ls%setting,S)
Call II_get_h1(lupri,luerr,ls%setting,H1)
Call get_initial_dens(H1,S,D,ls,config)
Call scfloop(H1,F,D,S,E,ls,config)
! New gradient
Traj%Gradient = 1.D0
! New accelerations
Do i = 1,NAtoms
   Do j = 1,3
    Traj%Accel(3*(i-1)+j)=-Traj%Gradient(3*(i-1)+j)/Traj%Mass(i)
   Enddo
Enddo
! Second velocity half-step
Traj%Velocities = Traj%Velocities + &
Traj%Accel*TStep*0.5D0
! End of Verlet step
Write(*,*)'New potential',E
Write(*,*)'New coordinates',Traj%Coordinates
Call Deallocate_traj(traj)
Deallocate (Config%dynamics%Initial_velocities)
!
End subroutine LS_dyn_run
!
End module Dynamics_driver
#endif
