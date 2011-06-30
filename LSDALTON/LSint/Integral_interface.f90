!> @file
!> Library-like integral-interface routines 

!> \brief Calculates overlap integral matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param S the overlap matrix
SUBROUTINE II_get_overlap(LUPRI,LUERR,SETTING,S)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(MATRIX)          :: S
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: nbast
real(realk)         :: TS,TE

!CALL LSHEADER(lupri,'II_get_overlap')
CALL LSTIMER('START ',TS,TE,LUPRI)
nbast = S%nrow
IF(SETTING%SCHEME%DEBUGOVERLAP)THEN
  call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,1)
  CALL ls_getIntegrals('Regular','Regular','Empty','Empty',&
       &1,1,'Overlap','Regular','Contracted',SETTING,LUPRI,LUERR)
ELSE
  call initIntegralOutputDims(setting%Output,nbast,1,nbast,1,1)
  CALL ls_getIntegrals('Regular','Empty','Regular','Empty',&
       &1,1,'Overlap','Regular','Contracted',SETTING,LUPRI,LUERR)
ENDIF
CALL retrieve_Output(lupri,setting,S)

CALL LSTIMER('OVERLAP',TS,TE,LUPRI)
END SUBROUTINE II_get_overlap

!> \brief Calculates overlap integral matrix between 2 different AO basis's
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
SUBROUTINE II_get_mixed_overlap(LUPRI,LUERR,SETTING)
  use files
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,IPRINT
!
TYPE(MATRIX)        :: S
Integer             :: nbast,i,j,LU,nbastaux

IPRINT=SETTING%SCHEME%INTPRINT
nbast = SETTING%MOLECULE(1)%p%nbastREG
nbastaux = SETTING%MOLECULE(1)%p%nbastAUX
call mat_init(S,nbast,nbastaux)
call initIntegralOutputDims(setting%output,nbast,nbastaux,1,1,1)
CALL ls_getIntegrals('Regular','DF-Aux','Empty','Empty',&
     &1,1,'Overlap','Regular','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,S)
WRITE(LUPRI,'(A,2X,F16.8)')'Mixed Overlap',mat_dotproduct(S,S)
IF(IPRINT.GT.1000)THEN
   call mat_print(S,1,nbast,1,nbastaux,lupri)
ENDIF
LU=-1
CALL LSOPEN(LU,'dualbasisoverlap','NEW','UNFORMATTED')
call mat_write_to_disk(LU,S)
CALL LSCLOSE(LU,'KEEP')
call mat_free(S)

END SUBROUTINE II_get_mixed_overlap

!> \brief Calculates one electron fock matrix contribution
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the one electron fock matrix contribution
SUBROUTINE II_get_h1(LUPRI,LUERR,SETTING,h)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(MATRIX),target   :: h
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: nbast
TYPE(MATRIX),target :: tmp
Real(realk)         :: OLDTHRESH
CALL II_get_nucel_mat(LUPRI,LUERR,SETTING,h)

nbast = h%nrow
CALL mat_init(tmp,nbast,nbast)

CALL II_get_kinetic(LUPRI,LUERR,SETTING,tmp)

call mat_daxpy(1.D0,tmp,h)
CALL mat_free(tmp)

END SUBROUTINE II_get_h1

!> \brief Calculates the kinetic energy fock matrix contribution
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param T the kinetic energy fock matrix contribution
SUBROUTINE II_get_kinetic(LUPRI,LUERR,SETTING,T)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(MATRIX),target   :: T
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer               :: nbast
real(realk)           :: OLDTHRESH
CALL mat_zero(T)
OLDTHRESH = SETTING%SCHEME%THRESHOLD
SETTING%SCHEME%THRESHOLD = SETTING%SCHEME%THRESHOLD*1.0D-5
nbast = T%nrow
call initIntegralOutputDims(setting%output,nbast,1,nbast,1,1)
CALL ls_getIntegrals('Regular','Empty','Regular','Empty',&
     &1,1,'Kinetic','Regular','Contracted',SETTING,LUPRI,LUERR)
SETTING%SCHEME%THRESHOLD = OLDTHRESH
CALL retrieve_Output(lupri,setting,T)

END SUBROUTINE II_get_kinetic

!> \brief Calculates the nuclear attraction fock matrix contribution
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the nuclear attraction fock matrix contribution
SUBROUTINE II_get_nucel_mat(LUPRI,LUERR,SETTING,h)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(MATRIX),target   :: h
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: nbast
TYPE(MATRIX)        :: temp
real(realk)         :: OLDTHRESH

CALL mat_zero(h)
OLDTHRESH = SETTING%SCHEME%THRESHOLD
SETTING%SCHEME%THRESHOLD = SETTING%SCHEME%THRESHOLD*1.0D-5
nbast = h%nrow
! Caclculate multipole moments when using FMM
IF(SETTING%SCHEME%FMM) THEN
  call mat_init(temp,nbast,nbast)
  CALL mat_zero(temp)
  CALL ls_attachDmatToSetting(temp,1,setting,'LHS',1,2,lupri)
  CALL ls_attachDmatToSetting(temp,1,setting,'RHS',3,4,lupri)
  call ls_multipolemoment(LUPRI,LUERR,SETTING,nbast,0,nbast,nbast,nbast,nbast,&
     &'Regular','Regular','Regular','Regular','Coulomb','Contracted',.FALSE.)
  CALL ls_freeDmatFromSetting(setting)
  call mat_free(temp)
ENDIF 
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,1)
CALL ls_getIntegrals('Regular','Regular','Nuclear','Empty',&
     &1,1,'Nucrep','Regular','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,h)
SETTING%SCHEME%THRESHOLD = OLDTHRESH

END SUBROUTINE II_get_nucel_mat

!> \brief Calculates the full molecular gradient
!> \author T. Kjaergaard
!> \date 2010-04-26
!> \param lupri Default print unit
!> \param F The Fock/Kohn-Sham matrix
!> \param D density matrix
!> \param setting Integral evalualtion settings
!> \param dodft is it a DFT or HF calc 
!> \param doprint should we print the gradient
subroutine II_get_molecular_gradient(GRAD,lupri,F,D,setting,dodft,doprint)
  use precision
  use Matrix_module
  use Matrix_Operations
  use TYPEDEF  
  use ls_util
  use LSTIMING
  implicit none
  integer,intent(in)        :: lupri
  type(Matrix),intent(in)   :: F,D
  TYPE(LSSETTING)           :: SETTING
  logical                   :: dodft,doprint
!
  type(matrixp)               :: Dmat(1)
  type(Matrix),target         :: tempm1,tempm3
  type(Matrix)                :: tempm2
  integer   :: nbast,natom,ndmat,ix,iatom,luerr
  real(realk),pointer  :: tmpGRAD(:,:)
  real(realk), intent(out) :: GRAD(3,SETTING%MOLECULE(1)%p%Natoms)
  real(realk) :: nrm,ts,te

  call lstimer('START ',ts,te,lupri)

 if(doprint) then
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(8X,A)') '*************************************************************'
    write(lupri,'(8X,A)') '*            MOLECULAR GRADIENT RESULTS (in a.u.)           *'
    write(lupri,'(8X,A)') '*************************************************************'
    write(lupri,*) 
    write(lupri,*) 
 end if

  nbast = D%nrow    
  natom = setting%molecule(1)%p%natoms
  ndmat = 1
  luerr = 0
  call mat_init(tempm1,nbast,nbast)
  tempm1=D
  call mat_scal(2.d0,tempm1)
  call II_get_nn_gradient(Grad,setting,lupri,luerr)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,GRAD,natom,'nuclear rep')
  call lstimer('nn-grad ',ts,te,lupri)
  Dmat(1)%p => tempm1
  call mem_alloc(tmpGrad,3,nAtom)
  CALL II_get_twoElectron_gradient(tmpGrad,natom,Dmat,Dmat,ndmat,ndmat,setting,lupri,luerr)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'twoElectron')
  call lstimer('two-grad',ts,te,lupri)
  CALL DAXPY(3*natom,1.d0,tmpGrad,1,Grad,1)
  CALL II_get_ne_gradient(tmpGrad,Dmat,ndmat,setting,lupri,luerr)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'nuclear att')
  call lstimer('ne-grad ',ts,te,lupri)
  CALL DAXPY(3*natom,1.d0,tmpGrad,1,Grad,1)
  if (dodft) THEN
     CALL II_get_xc_geoderiv_molgrad(lupri,luerr,setting,nbast,D,tmpGrad,natom)
     if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'exchange-corr')
     CALL DAXPY(3*natom,1.d0,tmpGrad,1,Grad,1)
     call lstimer('xc-grad ',ts,te,lupri)
  endif
  CALL II_get_kinetic_gradient(tmpGrad,Dmat,ndmat,setting,lupri,luerr)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'kinetic')
  CALL DAXPY(3*natom,1.d0,tmpGrad,1,Grad,1)
  call lstimer('kin-grad',ts,te,lupri)
  call mat_init(tempm2,nbast,nbast)
  call mat_init(tempm3,nbast,nbast)
  call mat_mul(tempm1,F,'n','n',1.d0,0.d0,tempm2)
  call mat_mul(tempm2,D,'n','n',-1.d0,0.d0,tempm3)
  call mat_free(tempm1)
  call mat_free(tempm2)

  Dmat(1)%p => tempm3 !the DFD mat
  CALL II_get_reorthoNormalization(tmpGrad,Dmat,ndmat,setting,lupri,luerr)
  call mat_free(tempm3)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'reorthonomal')
  call lstimer('reorth-g',ts,te,lupri)
  CALL DAXPY(3*natom,1.d0,tmpGrad,1,Grad,1)
  CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,GRAD,natom,'TOTAL')

 if(doprint) then
    nrm = 0.d0
    DO iatom = 1,natom
      DO ix=1,3
        nrm = nrm + Grad(ix,iatom)*Grad(ix,iatom)
      ENDDO
    ENDDO
    nrm = sqrt(nrm/3/natom)

    write(lupri,*) 
    write(lupri,*)
    write(lupri,'(1X,A24,F23.16)') 'RMS gradient norm (au): ',nrm
    write(lupri,*) 
    write(lupri,*) 
 end if
 


  call mem_dealloc(tmpGrad)

end subroutine II_GET_MOLECULAR_GRADIENT

!> \brief Calculates the two-electron contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param twoGrad The two-electron gradient
!> \param natoms The number of atoms
!> \param DmatLHS The left-hand-side density matrix
!> \param DmatRHS The reft-hand-side density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_twoElectron_gradient(twoGrad,natoms,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
IMPLICIT NONE
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs,natoms
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)     :: twoGrad(3,natoms)
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
Real(realk),pointer :: exchangeGrad(:,:)

!nAtoms = setting%molecule(1)%p%nAtoms
CALL mem_alloc(exchangeGrad,3,nAtoms)

CALL II_get_J_gradient(twoGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,twoGRAD,natoms,'Coulomb')

exchangeGrad = 0d0

CALL II_get_K_gradient(exchangeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,exchangeGRAD,natoms,'exchange')

CALL DAXPY(3*nAtoms,1.d0,exchangeGrad,1,twoGrad,1)
CALL mem_dealloc(exchangeGrad)

END SUBROUTINE II_get_twoElectron_gradient

!> \brief Calculates the exchange contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param kGrad The exchange gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_K_gradient(kGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(OUT)       :: kGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
integer               :: nAtoms,iDmat,I,J
real(realk)           :: Factor
CHARACTER(len=7)      :: Oper
logical               :: Dsym

Integer                :: symLHS(ndlhs),symRHS(ndrhs)
Integer                :: nlhs,nrhs
Type(matrixp)          :: DLHS(2*ndlhs),DRHS(2*ndrhs)

! Check first if the exchange-contribution should be calculated. If not exit 
! this subroutine
IF (SETTING%SCHEME%exchangeFactor.EQ.0.0d0) RETURN

! Check symetry. Split non-symmetric matrices to symmetric and anti symmetric parts. 
! Symmetric, anti-symmetric and anti-symmetric, symmetric paris vanish.
call II_split_dmats(DmatLHS,DmatRHS,ndlhs,ndrhs,DLHS,DRHS,symLHS,symRHS,nlhs,nrhs)

IF (nlhs.EQ.0) RETURN

! Actual calculation

nAtoms = setting%molecule(1)%p%nAtoms

CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',2,4,lupri)
CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',1,3,lupri)

IF (SETTING%SCHEME%CAM) THEN
  Oper = 'Erfc   '   !Coulomb attenuated method
ELSE
  Oper = 'Coulomb'   !Regular Coulomb metric 
ENDIF

!Calculates the HF-exchange contribution

call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
Dsym = .TRUE.
call ls_get_exchange_mat(Dsym,'Regular','Regular','Regular','Regular',&
     &                   Oper,'Gradient','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,kGrad)

Factor = SETTING%SCHEME%exchangeFactor
do I=1,3
   do J=1,natoms
      kgrad(I,J)=Factor*kGrad(I,J)
   enddo
enddo

CALL ls_freeDfull(setting)
setting%RHSdmat = .FALSE.
setting%LHSdmat = .FALSE.
!call mem_dealloc(exchangeGradient)

! Free allocated memory
call II_free_split_dmats(DLHS,DRHS,symLHS,symRHS,ndlhs,ndrhs)


END SUBROUTINE II_get_K_gradient

!> \brief Driver for calculating the electron-electron repulsion contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param eeGrad The electron-electron-repulsion gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_J_gradient(eeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)     :: eeGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
Integer                :: symLHS(ndlhs),symRHS(ndrhs)
Integer                :: nlhs,nrhs
Type(matrixp)          :: DLHS(ndlhs),DRHS(ndrhs)

! Check symetry and remove all anti-symmetric components (which will not contribute) 
! to the Coulomb ERI gradient
call II_symmetrize_dmats(DmatLHS,DmatRHS,ndlhs,ndrhs,DLHS,DRHS,symLHS,symRHS,nlhs,nrhs)

IF (nlhs.EQ.0) RETURN

! Actual calculation
IF (setting%scheme%densfit) THEN
  CALL II_get_df_J_gradient(eeGrad,DLHS,DRHS,nlhs,nrhs,setting,lupri,luerr)
ELSE
  CALL II_get_J_gradient_regular(eeGrad,DLHS,DRHS,nlhs,nrhs,setting,lupri,luerr)
ENDIF

! Free allocated memory
call II_free_symmetrized_dmats(DLHS,DRHS,symLHS,symRHS,ndlhs,ndrhs)

END SUBROUTINE II_get_J_gradient

!> \brief Calculates the density-fitted electron-electron repulsion contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param eeGrad The electron-electron-repulsion gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_df_J_gradient(eeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)     :: eeGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
!
integer                   :: nAtoms,iDmat,nlhs,nrhs
type(matrixp)             :: eeGradMat(1)
type(matrix),target       :: eeGradTarget
Type(matrixp),pointer     :: DLHS(:),DRHS(:)
logical                   :: same
TYPE(matrix),target       :: calpha ! TEMPORARY SOLUTION
Character(80)             :: Filename
Real(realk)               :: eeGradtmp(3,setting%molecule(1)%p%nAtoms)

integer :: nbast,naux

! Calculate multipole moments when using FMM
IF(SETTING%SCHEME%FMM) THEN
  CALL LSQUIT('Electron-electron gradient not yet implemented with FMM',-1)
ENDIF 

nbast = DmatLHS(1)%p%nrow
naux  = setting%molecule(1)%p%nbastAUX

same = ls_same_mats(DmatLHS,DmatRHS,ndlhs,ndrhs)
IF (.NOT.same) CALL lsQUIT('Error in II_get_df_J_gradient. Only working for equal density matrices!',-1)
IF ((ndlhs.GT.1).OR.(ndrhs.GT.1)) CALL lsQUIT('Error in II_get_df_J_gradient. Only working for ndmat = 1!',-1)
Filename = 'LSCALPHA'
IF (.not.io_file_exist(Filename,SETTING%IO)) call lsquit('Error in II_get_df_J_gradient. CALPHA does not exsit!',-1)


CALL mat_init(calpha,naux,1)
call io_read_mat(calpha,Filename,Setting%IO,lupri,luerr)
!call mat_print(CALPHA,1,naux,1,1,lupri)
! ADD reading of fitting coefficients from file CALPHA

NULLIFY(DLHS)
NULLIFY(DRHS)
nlhs = ndlhs
nrhs = ndrhs
ALLOCATE(DLHS(nlhs))
ALLOCATE(DRHS(nrhs))

eeGrad = 0.d0

!****************** Calculate (rho^e|tilde rho) ***********
DO idmat = 1,nlhs
DLHS(idmat)%p => DmatLHS(idmat)%p
ENDDO
DO idmat = 1,nrhs
  DRHS(idmat)%p => calpha
ENDDO

CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',1,2,lupri)
CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',3,4,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

eeGradtmp = 0.d0
call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
CALL ls_jengine('Regular','Regular','DF-Aux','Empty',&
     &          'Coulomb','Gradient','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,eeGRADtmp)
CALL ls_freeDmatFromSetting(setting)
eeGrad = eeGrad + 2.d0*eeGradtmp
!**************** End calculate (rho^e|tilde rho) *********

!****************** Calculate (tilde rho^e|rho) ***********
DO idmat = 1,nlhs
DLHS(idmat)%p => calpha
ENDDO
DO idmat = 1,nrhs
  DRHS(idmat)%p => DmatRHS(idmat)%p
ENDDO

CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',1,2,lupri)
CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',3,4,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

eeGradtmp = 0.d0
call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
CALL ls_jengine('DF-Aux','Empty','Regular','Regular',&
     &          'Coulomb','Gradient','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,eeGRADtmp)
CALL ls_freeDmatFromSetting(setting)
eeGrad = eeGrad + eeGradtmp
!**************** End calculate (tilde rho^e|rho) *********

!**************** Calculate (tilde rho^e|tilde rho) ***********
DO idmat = 1,nlhs
DLHS(idmat)%p => calpha
ENDDO
DO idmat = 1,nrhs
  DRHS(idmat)%p => calpha
ENDDO

CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',1,2,lupri)
CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',3,4,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

eeGRADtmp = 0.d0
call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
CALL ls_jengine('DF-Aux','Empty','DF-Aux','Empty',&
     &          'Coulomb','Gradient','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,eeGRADtmp)
eeGrad = eeGrad - 2.d0*eeGRADtmp
CALL ls_freeDmatFromSetting(setting)
!************** End calculate (tilde rho^e|tilde rho) *********

DEALLOCATE(DLHS)
DEALLOCATE(DRHS)

CALL mat_free(calpha)

END SUBROUTINE II_get_df_J_gradient

!> \brief Calculates the (regular) electron-electron repulsion contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param eeGrad The electron-electron-repulsion gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_J_gradient_regular(eeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)     :: eeGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
!
integer                   :: nAtoms,iDmat,nlhs,nrhs
type(matrixp)             :: eeGradMat(1)
type(matrix),target       :: eeGradTarget
Type(matrixp),pointer     :: DLHS(:),DRHS(:)
logical                   :: same

integer :: nbast

! Calculate multipole moments when using FMM
IF(SETTING%SCHEME%FMM) THEN
  CALL LSQUIT('Electron-electron gradient not yet implemented with FMM',-1)
ENDIF 

nbast = DmatLHS(1)%p%nrow

same = ls_same_mats(DmatLHS,DmatRHS,ndlhs,ndrhs)
NULLIFY(DLHS)
NULLIFY(DRHS)
IF (.NOT.same.AND.ndlhs.EQ.ndrhs) THEN
! If the density-matrices are not the same we make a dual set of matrices -
! for two matrices DLHS = D and DRHS = P we have the two contributions
! D_ab(ab'|cd)P_cd and P_ab(ab'|cd)D_cd  (since we only differentiate on the LHS)
  nlhs = 2*ndlhs
  nrhs = 2*ndrhs
  ALLOCATE(DLHS(nlhs))
  ALLOCATE(DRHS(nrhs))
  DO idmat = 1,ndrhs
    DLHS(idmat)%p       => DmatLHS(idmat)%p
    DLHS(ndlhs+idmat)%p => DmatRHS(idmat)%p
    DRHS(idmat)%p       => DmatRHS(idmat)%p
    DRHS(ndlhs+idmat)%p => DmatLHS(idmat)%p
  ENDDO
ELSE
  nlhs = ndlhs
  nrhs = ndrhs
  ALLOCATE(DLHS(nlhs))
  ALLOCATE(DRHS(nrhs))
  DO idmat = 1,nlhs
    DLHS(idmat)%p => DmatLHS(idmat)%p
  ENDDO
  DO idmat = 1,nrhs
    DRHS(idmat)%p => DmatRHS(idmat)%p
  ENDDO
ENDIF

CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',1,2,lupri)
CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',3,4,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
CALL ls_jengine('Regular','Regular','Regular','Regular',&
     &          'Coulomb','Gradient','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,eeGRAD)

CALL ls_freeDmatFromSetting(setting)

DEALLOCATE(DLHS)
DEALLOCATE(DRHS)

IF (.NOT.same) call dscal(3*nAtoms,0.5d0,eeGrad,1)

END SUBROUTINE II_get_J_gradient_regular


!> \brief Calculates the one-electron contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param oneGrad The one-electron gradient contribution
!> \param Dmat The density matrices
!> \param ndmat The number of density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_oneElectron_gradient(oneGrad,Dmat,ndmat,setting,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(OUT)       :: oneGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrixp),intent(IN)      :: Dmat(ndmat)
!
Real(realk),pointer :: gradCont(:,:)
Integer :: nAtoms
!
! Nuclear-electron attraction gradient
CALL II_get_ne_gradient(oneGrad,Dmat,ndmat,setting,lupri,luerr)

nAtoms = setting%molecule(1)%p%nAtoms
CALL mem_alloc(gradCont,3,nAtoms)
! Kinetic energy gradient
CALL II_get_kinetic_gradient(gradCont,Dmat,ndmat,setting,lupri,luerr)
!
CALL DAXPY(3*nAtoms,1.d0,gradCont,1,oneGrad,1)
CALL mem_dealloc(gradCont)
!
END SUBROUTINE II_get_oneElectron_gradient


!> \brief Calculates the nuclear-electron attraction contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param neGrad The nuclear-electron attraction gradient
!> \param Dmat The density matrices
!> \param ndmat The number of density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_ne_gradient(neGrad,Dmat,ndmat,setting,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(OUT)       :: neGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrixp),intent(IN)      :: Dmat(ndmat)
!
integer             :: nAtoms
real(realk)         :: OLDTHRESH

OLDTHRESH = SETTING%SCHEME%THRESHOLD
SETTING%SCHEME%THRESHOLD = SETTING%SCHEME%THRESHOLD*1.0D-5

! Caclculate multipole moments when using FMM
IF(SETTING%SCHEME%FMM) THEN
  CALL LSQUIT('Nuclear-electronic gradient not yet implemented with FMM',lupri)
ENDIF 

CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'LHS',1,2,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

call initIntegralOutputDims(setting%output,3,nAtoms,1,1,1)
CALL ls_getIntegrals('Regular','Regular','Nuclear','Empty',&
     &ndmat,1,'Nucrep','Gradient','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,neGrad)

CALL ls_freeDmatFromSetting(setting)
SETTING%SCHEME%THRESHOLD = OLDTHRESH

END SUBROUTINE II_get_ne_gradient

!> \brief Calculates the kinetic energy contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param kinGrad The kinetic energy gradient
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_kinetic_gradient(kinGrad,Dmat,ndmat,setting,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(OUT)       :: kinGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrixp),intent(IN)      :: Dmat(ndmat)
!
integer             :: nAtoms
real(realk)         :: OLDTHRESH

OLDTHRESH = SETTING%SCHEME%THRESHOLD
SETTING%SCHEME%THRESHOLD = SETTING%SCHEME%THRESHOLD*1.0D-5

CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'LHS',1,3,lupri)

nAtoms = setting%molecule(1)%p%nAtoms
call initIntegralOutputDims(setting%output,3,nAtoms,1,1,1)
CALL ls_getIntegrals('Regular','Empty','Regular','Empty',ndmat,1,&
     &'Kinetic','Gradient','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,kinGrad)
CALL ls_freeDmatFromSetting(setting)

SETTING%SCHEME%THRESHOLD = OLDTHRESH

END SUBROUTINE II_get_kinetic_gradient

!> \brief Calculates the nuclear repulsion energy contribution
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nucpot the nuclear repulsion energy contribution
SUBROUTINE II_get_nucpot(LUPRI,LUERR,SETTING,NUCPOT)
  use precision
  use TYPEDEF  
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
integer             :: usemat
INTEGER               :: LUPRI,LUERR
REAL(realk)           :: nucpot
Integer               :: I,J
real(realk)           :: pq(3),distance

NUCPOT=0
DO I=1,SETTING%MOLECULE(1)%p%Natoms
  DO J=I+1,SETTING%MOLECULE(1)%p%Natoms
    pq(1) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(1)
    pq(2) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(2)
    pq(3) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(3)
    Distance = sqrt(pq(1)*pq(1)+pq(2)*pq(2)+pq(3)*pq(3))
    NUCPOT=NUCPOT+SETTING%MOLECULE(1)%p%ATOM(I)%Charge*SETTING%MOLECULE(1)%p%ATOM(J)%Charge/Distance
 ENDDO
ENDDO
END SUBROUTINE II_get_nucpot

!> \brief Calculates the nuclear-nuclear contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-02-29
!> \param nucGrad The nuclear-nuclear molecular gradient vector
!> \param settings The settings for integral evaluation
SUBROUTINE II_get_nn_gradient(nucGrad,SETTING,lupri,luerr)
  use precision
  use TYPEDEF  
IMPLICIT NONE
TYPE(LSSETTING),intent(IN) :: SETTING
REAL(realk),intent(OUT)    :: nucGrad(3,SETTING%MOLECULE(1)%p%Natoms)
Integer,intent(IN)         :: lupri,luerr
!
Integer     :: I,J
real(realk) :: pq(3),distance,temp

nucGrad=0.d0
DO I=1,SETTING%MOLECULE(1)%p%Natoms
  DO J=I+1,SETTING%MOLECULE(1)%p%Natoms
    pq(1) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(1)
    pq(2) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(2)
    pq(3) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(3)
    Distance = sqrt(pq(1)*pq(1)+pq(2)*pq(2)+pq(3)*pq(3))
    temp = SETTING%MOLECULE(1)%p%ATOM(I)%Charge*&
     &     SETTING%MOLECULE(1)%p%ATOM(J)%Charge/Distance**3
    nucGrad(1,I) = nucGrad(1,I) - pq(1)*temp
    nucGrad(2,I) = nucGrad(2,I) - pq(2)*temp
    nucGrad(3,I) = nucGrad(3,I) - pq(3)*temp
    nucGrad(3,J) = nucGrad(3,J) + pq(3)*temp
    nucGrad(1,J) = nucGrad(1,J) + pq(1)*temp
    nucGrad(2,J) = nucGrad(2,J) + pq(2)*temp
 ENDDO
ENDDO
END SUBROUTINE II_get_nn_gradient

!> \brief Calculates the cartesian moments 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param carmom the cartesian moments 
!> \param nmat the number of matrices 
!> \param nderiv The order of cartesian moments
SUBROUTINE II_get_carmom(LUPRI,LUERR,SETTING,carmom,nmat,nderiv)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,nderiv,nmat
integer             :: usemat
TYPE(MATRIX),target   :: carmom(nmat)
TYPE(MATRIX)          :: TMP
TYPE(LSSETTING)       :: SETTING
Real(realk),pointer   :: integrals(:,:,:,:,:)
INTEGER               :: nbast,I
type(matrixp)         :: intmat(nmat)

nbast = carmom(1)%nrow
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nmat)
CALL ls_getIntegrals('Regular','Regular','Empty','Empty',&
                      &nmat,nderiv,'Carmom','Regular','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,carmom)

END SUBROUTINE II_get_carmom

!> \brief Calculates the dipole moment contribution from the nuclei
!> \param setting Integral evalualtion settings
!> \param nucdip Nuclear contrib to dipole moment
SUBROUTINE II_get_nucdip(setting,nucdip)
  use precision
  use TYPEDEF
IMPLICIT NONE
TYPE(LSSETTING),intent(in) :: setting
REAL(realk),intent(out)    :: nucdip(3)
INTEGER            :: i
TYPE(ATOM),pointer :: mol(:)
nucdip = 0
mol => setting%MOLECULE(1)%p%ATOM
DO i=1,setting%MOLECULE(1)%p%Natoms
  nucdip = nucdip + mol(i)%Charge * mol(i)%CENTER
ENDDO
END SUBROUTINE II_get_nucdip

!> \brief Calculates property integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param mat the matrix of property integrals
!> \param nmat the number of matrices 
!> \param STRING the label of property integral
!>
!> IMPLEMENTED SO FAR 
!> DIPLEN :3  XDIPLEN ,YDIPLEN ,ZDIPLEN
!> SECMOM :6  XXSECMOM,XYSECMOM,XZSECMOM,YYSECMOM,YZSECMOM,ZZSECMOM
!> WILL IMPLEMT AT SOME POINT
!> ANGMOM :3  XANGMOM,YANGMOM,ZANGMOM
!>
SUBROUTINE II_get_integral(LUPRI,LUERR,SETTING,mat,nmat,STRING)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
character(len=7)     :: STRING
TYPE(LSSETTING)      :: SETTING
INTEGER              :: LUPRI,LUERR,nderiv,nmat
TYPE(MATRIX),target  :: mat(nmat)
!
Character(len=6)     :: Oper
TYPE(MATRIX),pointer :: TMP(:)
INTEGER              :: nbast,I,nmat2

nbast = mat(1)%nrow
Oper = 'Carmom'
SELECT CASE(STRING)
CASE('DIPLEN ')
   IF(nmat.NE.3)call LSQUIT('II_get_integral:DIPLEN INPUT REQUIRES A DIM = 3',lupri)
   nmat2=4
   nderiv = 1
CASE('SECMOM ')
   IF(nmat.NE.6)call LSQUIT('II_get_integral:SECMOM INPUT REQUIRES A DIM = 6',lupri)
   nmat2=10
   nderiv = 2
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in II_get_integral =',STRING
   CALL LSQUIT('Wrong case in II_get_integral',lupri)
END SELECT

call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nmat2)
CALL ls_getIntegrals('Regular','Regular','Empty','Empty',&
     &nmat2,nderiv,Oper,'Regular','Contracted',SETTING,LUPRI,LUERR)

ALLOCATE(TMP(nmat2))
DO I=1,nmat2
   call mat_init(TMP(I),nbast,nbast)
ENDDO
CALL retrieve_Output(lupri,setting,TMP)

SELECT CASE(STRING)
CASE('DIPLEN ')
   DO I=1,3
      call mat_assign(mat(I),tmp(I+1))
   ENDDO
CASE('SECMOM ')
   DO I=1,6
      call mat_assign(mat(I),tmp(I+4))
   ENDDO
END SELECT

DO I=1,nmat2
   call mat_free(TMP(I))
ENDDO
DEALLOCATE(TMP)

END SUBROUTINE II_get_integral

!> \brief Calculates the spherical moments
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param sphmat the matrices of the spherical moments
!> \param nsphmat the number of matrices 
!> \param nderiv the order of moments
SUBROUTINE II_get_sphmom(LUPRI,LUERR,SETTING,sphmom,nsphmat,nderiv)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,nderiv,nsphmat,usemat,nbast,I
TYPE(MATRIX),target   :: sphmom(nsphmat)
TYPE(LSSETTING)       :: SETTING

nbast = sphmom(1)%nrow
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nsphmat)
CALL ls_getIntegrals('Regular','Regular','Empty','Empty',nsphmat,nderiv,'Sphmom','Regular','Contracted',&
     &                SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,sphmom)

END SUBROUTINE II_get_sphmom

!> \brief Calculates the 3 center overlap integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupriold old Default print unit
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nbast The number of basis functions
!> \param auxoption If the auxiliary basis should be used
SUBROUTINE II_get_3center_overlap(LUPRIOLD,LUPRI,LUERR,SETTING,nbast,auxoption)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use tco_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast,lupriold,i,j,k,nbastaux,iprint
TYPE(MATRIX),target   :: TMP
integer             :: usemat
LOGICAL               :: auxoption
Real(realk),pointer   :: integrals(:,:,:,:,:)
Real(realk)           :: SUM
Real(realk),pointer   :: tco_integrals(:,:,:)
Real(realk)           :: tco_SUM,diff
type(matrixp)         :: intmat(1)

iprint=SETTING%SCHEME%INTPRINT
IF(auxoption)THEN
   nbastaux = SETTING%MOLECULE(1)%p%nbastAUX
   call initIntegralOutputDims(setting%output,nbast,nbast,1,nbastaux,1)
   CALL ls_getIntegrals('Regular','Regular','Empty','DF-Aux',&
        &1,1,'Overlap','Regular','Contracted',SETTING,LUPRI,LUERR)
   call mem_alloc(integrals,nbast,nbast,1,nbastaux,1)
   CALL retrieve_Output(lupri,setting,integrals)

   call mem_alloc(tco_integrals,nbast,nbast,nbastaux)
   CALL ls_DZERO(tco_integrals,nbast*nbast*nbastaux)
   CALL tco_getIntegrals(tco_integrals,'Regular','Regular','DF-Aux','Contracted',&
                         SETTING,LUPRI,LUERR)
ELSE
   nbastaux = nbast
   call initIntegralOutputDims(setting%output,nbast,nbast,1,nbast,1)
   CALL ls_getIntegrals('Regular','Regular','Empty','Regular',&
        &1,1,'Overlap','Regular','Contracted',SETTING,LUPRI,LUERR)
   call mem_alloc(integrals,nbast,nbast,1,nbast,1)
   CALL retrieve_Output(lupri,setting,integrals)

   call mem_alloc(tco_integrals,nbast,nbast,nbast)
   CALL ls_DZERO(tco_integrals,nbast*nbast*nbast)
   CALL tco_getIntegrals(tco_integrals,'Regular','Regular','Regular','Contracted',&
                         SETTING,LUPRI,LUERR)
ENDIF

SUM=0
tco_SUM=0
diff=0
IF(IPRINT.GT.100)THEN
   WRITE(LUPRIOLD,'(2X,A)')'3 Center Overlap (a|b|c)'
   WRITE(LUPRIOLD,'(2X,A)')'FORMAT IS (a,b) c=1, c=2 , ..'
ENDIF
DO i=1,nbast
   DO j=1,nbast
      IF(IPRINT.GT.100)THEN
         WRITE(LUPRIOLD,'(2X,A1,I2,A1,I2,A2,5F16.6/,(10X,5F16.6))')&
              &'(',i,',',j,')=',&
              &(Integrals(i,j,1,k,1),k=1,nbastaux)
      ENDIF
      DO k=1,nbastaux
         SUM=SUM+Integrals(i,j,1,k,1)
         tco_SUM=tco_SUM+tco_integrals(i,j,k)
         diff=diff+abs(integrals(i,j,1,k,1)-tco_integrals(i,j,k))
      ENDDO
   ENDDO
ENDDO

WRITE(LUPRIOLD,'(2X,A,F20.12)')'SUM OF ALL 3 CENTER OVERLAPS '
WRITE(LUPRIOLD,'(2X,A,F20.12)')'FROM THE STANDARD DRIVER ',SUM
WRITE(LUPRIOLD,'(2X,A,F20.12)')'FROM THE TCO DRIVER      ',tco_SUM
WRITE(LUPRIOLD,'(2X,A,D20.9)') 'ABSOLUTE DIFFERENCE      ',diff

call mem_dealloc(integrals)
call mem_dealloc(tco_integrals)

END SUBROUTINE II_get_3center_overlap

!> \brief Calculates the 2 center electron repulsion integrals (eri)
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param F the matrix containing the 2 center eri
SUBROUTINE II_get_2center_eri(LUPRI,LUERR,SETTING,F)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(MATRIX),target   :: F
integer               :: usemat
TYPE(LSSETTING)       :: SETTING
TYPE(INTEGRALINPUT)   :: INTINPUT
TYPE(INTEGRALOUTPUT)  :: INTOUTPUT
INTEGER               :: LUPRI,LUERR
!
Real(realk),pointer :: integrals(:,:,:,:,:)
Integer             :: nbast
type(matrixp)         :: intmat(1)

nbast = F%nrow
call initIntegralOutputDims(setting%output,nbast,1,nbast,1,1)
CALL ls_getIntegrals('Regular','Empty','Regular','Empty',1,1,'Coulomb','Regular','Contracted',&
      &              SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F)
END SUBROUTINE II_get_2center_eri

!
! THESE SUBROUTINES ARE NEVER CALLED AND THEREFORE NEVER TESTED
! AND THEY CANNOT BE TRUSTED
!
!!$SUBROUTINE II_get_2center_aux_eri(LUPRI,LUERR,SETTING,F)
!!$  use precision
!!$  use TYPEDEF  
!!$  use Matrix_module
!!$  use Matrix_Operations
!!$  use ls_Integral_Interface
!!$IMPLICIT NONE
!!$TYPE(MATRIX),target   :: F
!!$TYPE(MATRIX)          :: TMP
!!$integer             :: usemat
!!$TYPE(LSSETTING)       :: SETTING
!!$TYPE(INTEGRALINPUT)   :: INTINPUT
!!$TYPE(INTEGRALOUTPUT)  :: INTOUTPUT
!!$INTEGER               :: LUPRI,LUERR
!!$!
!!$Real(realk),pointer :: integrals(:,:,:,:,:)
!!$Integer             :: nbastaux
!!$type(matrixp)         :: intmat(1)
!!$
!!$nbastaux = SETTING%MOLECULE(1)%p%nbastAUX
!!$call mem_alloc(integrals,1,1,1,1,1)
!!$!NULLIFY(integrals)
!!$!ALLOCATE(integrals(1,1,1,1,1))
!!$call ls_dzero(integrals,nbastaux*nbastaux)
!!$usemat=2
!!$intmat(1)%p => F
!!$CALL ls_getIntegrals(integrals,intmat,usemat,'DF-Aux','Empty','DF-Aux','Empty',1,1,'Coulomb','Regular','Contracted',&
!!$      &              SETTING,LUPRI,LUERR)
!!$call mem_dealloc(integrals)
!!$!DEALLOCATE(integrals)
!!$END SUBROUTINE II_get_2center_aux_eri
!!$
!!$SUBROUTINE II_get_3center_eri(LUPRIOLD,LUPRI,LUERR,SETTING,auxoption)
!!$  use precision
!!$  use TYPEDEF  
!!$  use Matrix_module
!!$  use Matrix_Operations
!!$  use ls_Integral_Interface
!!$IMPLICIT NONE
!!$TYPE(LSSETTING)       :: SETTING
!!$TYPE(MATRIX)          :: TMP(1)
!!$integer             :: usemat
!!$INTEGER               :: LUPRI,LUERR,nbast,lupriold
!!$LOGICAL               :: auxoption
!!$!
!!$Real(realk),pointer :: integrals(:,:,:,:,:)
!!$Integer :: nbastaux,i,j,k
!!$type(matrixp)         :: intmat(1)
!!$
!!$!NULLIFY(integrals)
!!$IF(auxoption)THEN
!!$!simen Remove nbast from calling list
!!$  nbast = SETTING%MOLECULE(1)%p%nbastREG
!!$  nbastaux = SETTING%MOLECULE(1)%p%nbastAUX
!!$  call mem_alloc(integrals,nbast,nbast,1,nbastaux,1)
!!$  !NULLIFY(integrals)
!!$  !ALLOCATE(integrals(nbast,nbast,1,nbastaux,1))
!!$  call ls_dzero(integrals,nbast*nbast*nbastaux)
!!$! CALL MAT_INIT(TMP(1),1,1)
!!$ usemat=0
!!$ CALL ls_getIntegrals(integrals,intmat,usemat,'Regular','Regular','Empty','DF-Aux',1,1,'Coulomb','Regular','Contracted',&
!!$      &                SETTING,LUPRI,LUERR)
!!$ CALL MAT_free(TMP(1))
!!$  WRITE(LUPRIOLD,'(2X,A)')'QQ 3 Center eri (a,b|alpha)'
!!$  DO i=1,nbast         
!!$     DO j=1,nbast
!!$        WRITE(LUPRIOLD,'(3X,A2,1X,A2,10X,A8,8X,A8,8X,A8)')' A',' B',' alpha=1',' alpha=2','  ...   '
!!$         WRITE(LUPRIOLD,'(2X,A1,I2,A1,I2,A2,5F12.6/,(10X,5F12.6))')&
!!$                &'(',i,',',j,')=',(integrals(i,j,1,k,1),k=1,nbastAUX)
!!$     ENDDO
!!$  ENDDO
!!$ELSE
!!$  nbast = SETTING%MOLECULE(1)%p%nbastREG
!!$  call mem_alloc(integrals,nbast,nbast,1,nbast,1)
!!$  !NULLIFY(integrals)
!!$  !ALLOCATE(integrals(nbast,nbast,1,nbast,1))
!!$  call ls_dzero(integrals,nbast*nbast*nbast)
!!$!  CALL MAT_INIT(TMP(1),1,1)
!!$  usemat=0
!!$  CALL ls_getIntegrals(integrals,intmat,usemat,'Regular','Regular','Empty','Regular',1,1,'Coulomb','Regular','Contracted',&
!!$      &                SETTING,LUPRI,LUERR)
!!$  CALL MAT_FREE(TMP(1))
!!$  WRITE(LUPRIOLD,'(2X,A)')'QQ 3 Center eri (a,b|c)'
!!$  DO i=1,nbast         
!!$     DO j=1,nbast
!!$        WRITE(LUPRIOLD,'(3X,A2,1X,A2,10X,A8,8X,A8,8X,A8)')' A',' B',' alpha=1',' alpha=2','  ...   '
!!$         WRITE(LUPRIOLD,'(2X,A1,I2,A1,I2,A2,5F12.6/,(10X,5F12.6))')&
!!$                &'(',i,',',j,')=',(integrals(i,j,1,k,1),k=1,nbast)
!!$     ENDDO
!!$  ENDDO
!!$ENDIF
!!$call mem_dealloc(integrals)
!!$!DEALLOCATE(integrals)
!!$
!!$END SUBROUTINE II_get_3center_eri
!!$
!!$SUBROUTINE II_get_3center_erfc_eri(LUPRIOLD,LUPRI,LUERR,SETTING,auxoption)
!!$  use precision
!!$  use TYPEDEF  
!!$  use Matrix_module
!!$  use Matrix_Operations
!!$  use ls_Integral_Interface
!!$IMPLICIT NONE
!!$TYPE(LSSETTING)       :: SETTING
!!$TYPE(MATRIX)          :: TMP(1)
!!$integer             :: usemat
!!$INTEGER               :: LUPRI,LUERR,lupriold
!!$LOGICAL               :: auxoption
!!$!
!!$Real(realk),pointer :: integrals(:,:,:,:,:)
!!$Integer :: nbastaux,i,j,k,nbast
!!$type(matrixp)         :: intmat(1)
!!$
!!$!CALL LSHEADER(lupri,'II_get_3center_erfc_eri')
!!$
!!$!NULLIFY(integrals)
!!$IF(auxoption)THEN
!!$  nbast = SETTING%MOLECULE(1)%p%nbastREG
!!$  nbastaux = SETTING%MOLECULE(1)%p%nbastAUX
!!$  call mem_alloc(integrals,nbast,nbast,1,nbastaux,1)
!!$  !NULLIFY(integrals)
!!$  !ALLOCATE(integrals(nbast,nbast,1,nbastaux,1))
!!$  call ls_dzero(integrals,nbast*nbast*nbastaux)
!!$!  CALL MAT_INIT(TMP(1),1,1)
!!$  usemat=0
!!$  CALL ls_getIntegrals(integrals,intmat,usemat,'Regular','Regular','Empty','DF-Aux',1,1,'Erfc   ',&
!!$       &'Regular','Contracted', SETTING,LUPRI,LUERR)
!!$!  CALL MAT_FREE(TMP(1))
!!$ELSE
!!$  nbast = SETTING%MOLECULE(1)%p%nbastREG
!!$  call mem_alloc(integrals,nbast,nbast,1,nbast,1)
!!$  !NULLIFY(integrals)
!!$  !ALLOCATE(integrals(nbast,nbast,1,nbast,1))
!!$  call ls_dzero(integrals,nbast*nbast*nbast)
!!$!  CALL MAT_INIT(TMP(1),1,1)
!!$  usemat=0
!!$  CALL ls_getIntegrals(integrals,intmat,usemat,'Regular','Regular','Empty','Regular',1,1,'Erfc',&
!!$       &'Regular','Contracted', SETTING,LUPRI,LUERR)
!!$!  CALL MAT_FREE(TMP(1))
!!$ENDIF
!!$call mem_dealloc(integrals)
!!$!DEALLOCATE(integrals)
!!$
!!$END SUBROUTINE II_get_3center_erfc_eri
!!$
!!$SUBROUTINE II_get_3center_erf_eri(LUPRIOLD,LUPRI,LUERR,SETTING,auxoption)
!!$  use precision
!!$  use TYPEDEF  
!!$  use Matrix_module
!!$  use Matrix_Operations
!!$  use ls_Integral_Interface
!!$IMPLICIT NONE
!!$TYPE(LSSETTING)       :: SETTING
!!$INTEGER               :: LUPRI,LUERR,lupriold
!!$TYPE(MATRIX)          :: TMP(1)
!!$integer             :: usemat
!!$LOGICAL               :: auxoption
!!$!
!!$Real(realk),pointer :: integrals(:,:,:,:,:)
!!$Integer :: nbastaux,i,j,k,nbast
!!$TYPE(MATRIXP)          :: intmat(1)
!!$
!!$!CALL LSHEADER(lupri,'II_get_3center_erf_eri')
!!$
!!$!NULLIFY(integrals)
!!$IF(auxoption)THEN
!!$  nbast = SETTING%MOLECULE(1)%p%nbastREG
!!$  nbastaux = SETTING%MOLECULE(1)%p%nbastAUX
!!$  call mem_alloc(integrals,nbast,nbast,1,nbastaux,1)
!!$  !NULLIFY(integrals)
!!$  !ALLOCATE(integrals(nbast,nbast,1,nbastaux,1))
!!$  call ls_dzero(integrals,nbast*nbast*nbastaux)
!!$!  CALL MAT_INIT(TMP(1),1,1)
!!$  usemat=0
!!$  CALL ls_getIntegrals(integrals,intmat,usemat,'Regular','Regular','Empty','DF-Aux',1,1,'Erf    ',&
!!$       &'Regular','Contracted', SETTING,LUPRI,LUERR)
!!$!  CALL MAT_FREE(TMP(1))
!!$ELSE
!!$  nbast = SETTING%MOLECULE(1)%p%nbastREG
!!$  call mem_alloc(integrals,nbast,nbast,1,nbast,1)
!!$  !NULLIFY(integrals)
!!$  !ALLOCATE(integrals(nbast,nbast,1,nbast,1))
!!$  call ls_dzero(integrals,nbast*nbast*nbast)
!!$!  CALL MAT_INIT(TMP(1),1,1)
!!$  usemat=0
!!$  CALL ls_getIntegrals(integrals,intmat,usemat,'Regular','Regular','Empty','Regular',1,1,'Erf ',&
!!$       &'Regular','Contracted', SETTING,LUPRI,LUERR)
!!$!  CALL MAT_FREE(TMP(1))
!!$ENDIF
!!$call mem_dealloc(integrals)
!!$!DEALLOCATE(integrals)
!!$
!!$END SUBROUTINE II_get_3center_erf_eri

!> \brief Calculates the coulomb matrix
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use lstiming
IMPLICIT NONE
Integer               :: ndmat
TYPE(MATRIX)          :: D(ndmat),F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
!Real(realk),pointer :: Dfull(:,:,:)
!Real(realk),pointer :: Ffull(:,:,:,:,:)
Real(realk)         :: TS,TE,fac
integer :: idmat

CALL LSTIMER('START ',TS,TE,LUPRI)
IF(SETTING%SCHEME%JENGINE)THEN
!   IF(SETTING%SCHEME%FMM)THEN
!     call II_calc_and_write_MMfile_for_FMM(LUPRI,LUERR,SETTING,D)
!     call write_MM_DENS_file(LUPRI,LUERR,SETTING,D)
!     CALL LSTIMER('CalcMM',TS,TE,LUPRI)
!   ENDIF
   call II_get_jengine_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
ELSE
  call II_get_default_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
ENDIF

fac = 2.d0
IF(matrix_type .EQ. mtype_unres_dense)fac = 1.d0
DO IDMAT=1,ndmat
   WRITE(lupri,*)'The Coulomb energy contribution ',fac*0.5d0*mat_dotproduct(D(IDMAT),F(IDMAT))
ENDDO

END SUBROUTINE II_get_coulomb_mat

!> \brief Calculates the coulomb matrix using the jengine method
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_jengine_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
Integer             :: ndmat
TYPE(MATRIX),target   :: D(ndmat),F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Real(realk),pointer :: Dfull(:,:,:)
Real(realk),pointer :: Ffull(:,:,:,:,:)
Real(realk)         :: TS,TE
integer             :: I
type(matrixp)       :: Jmat(1),Dmat(1)

!CALL LSHEADER(lupri,'II_get_coulomb_mat')
!CALL LSTIMER('START ',TS,TE,LUPRI)
!ndmat = 1
DO I=1,ndmat
   CALL MAT_ZERO(F(I))
ENDDO
CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,lupri)
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,lupri)
IF(SETTING%SCHEME%FMM)call ls_multipolemoment(LUPRI,LUERR,SETTING,F(1)%nrow,0,&
     & D(1)%nrow,D(1)%ncol,D(1)%nrow,D(1)%ncol,&
     & 'Regular','Regular','Regular','Regular','Coulomb','Contracted',.TRUE.)
CALL ls_freeDmatFromSetting(setting)
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,F(1)%nrow,F(1)%ncol,1,1,ndmat)
call ls_jengine('Regular','Regular','Regular','Regular','Coulomb','Regular','Contracted',&
     &          SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F)
CALL ls_freeDmatFromSetting(setting)
!CALL LSTIMER('Jengin',TS,TE,LUPRI)

END SUBROUTINE II_get_jengine_mat

!> \brief Calculates the coulomb matrix using default explicit 4 center int
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_default_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use ls_Integral_Interface
IMPLICIT NONE
Integer               :: ndmat
TYPE(MATRIX)          :: D(ndmat)
TYPE(MATRIX)          :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Real(realk),pointer :: Dfull(:,:,:)
Real(realk),pointer :: Ffull(:,:,:,:,:)
Real(realk)         :: TS,TE
integer             :: I
type(matrixp)       :: Jmat(1),Dmat(1)

!CALL LSHEADER(lupri,'II_get_coulomb_mat')

!ndmat = 1
CALL LSTIMER('START ',TS,TE,LUPRI)
!if UHF we add the densities together and then do J(the coulomb matrix)
!as the alpha and beta part of J is the same.
call mem_alloc(Dfull,D(1)%nrow,D(1)%ncol,ndmat)
IF(matrix_type .EQ. mtype_unres_dense)THEN
  DO I=1,ndmat
    CALL DCOPY(D(1)%nrow*D(1)%ncol,D(I)%elms,1,Dfull(:,:,I),1)
    CALL DAXPY(D(1)%nrow*D(1)%ncol,1.d0,D(I)%elmsb,1,Dfull(:,:,I),1)
    CALL DSCAL(D(1)%nrow*D(1)%ncol,0.5d0,Dfull(:,:,I),1)
  ENDDO
ELSE !CLOSED_SHELL
  DO I=1,ndmat
   call mat_to_full(D(I),1D0,Dfull(:,:,I))
  ENDDO
ENDIF

call initIntegralOutputDims(setting%Output,F(1)%nrow,F(1)%ncol,1,1,ndmat)
CALL ls_attachDmatToSetting(Dfull,ndmat,setting,'RHS',3,4,lupri)
call ls_get_coulomb_mat('Regular','Regular','Regular','Regular',&
     &                   'Coulomb','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F)
DO I=1,ndmat
   call mat_scal(2.d0,F(I))
ENDDO
CALL ls_freeDmatFromSetting(setting)

call mem_dealloc(Dfull)
CALL LSTIMER('reg-J ',TS,TE,LUPRI)

END SUBROUTINE II_get_default_coulomb_mat

!> \brief Calculates the coulomb matrix using density fitting
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
SUBROUTINE II_get_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
IMPLICIT NONE
TYPE(MATRIX),target   :: D,F
TYPE(MATRIX),target   :: galpha,calpha
integer               :: usemat
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: ndmat,nbasis,naux,info
Real(realk),pointer :: alphabetafull(:,:,:,:,:)
Real(realk),pointer :: galphafull(:,:,:,:,:)
Real(realk),pointer :: calphafull(:,:,:)
Real(realk)         :: TSTART,TEND
Character(80)       :: Filename
type(matrixp)       :: Jmat(1),Dmat(1)
type(matrixp)       :: Intmat(1)

IF (SETTING%SCHEME%OVERLAP_DF_J) THEN
  CALL II_get_overlap_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F)
ELSE
  CALL II_get_regular_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F)
ENDIF

END SUBROUTINE II_get_df_coulomb_mat

!> \brief Calculates the coulomb matrix using regular density fitting
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
SUBROUTINE II_get_regular_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use ls_Integral_Interface
!use linsolvdf
IMPLICIT NONE
TYPE(MATRIX),target   :: D,F
TYPE(MATRIX),target   :: galpha,calpha,calpha_old
integer               :: usemat
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: ndmat,nbasis,naux,info
Real(realk),pointer :: alphabetafull(:,:,:,:,:)
Real(realk),pointer :: galphafull(:,:,:,:,:)
Real(realk),pointer :: calphafull(:,:,:)
Real(realk)         :: TSTART,TEND
Character(80)       :: Filename
type(matrixp)       :: Jmat(1),Dmat(1)
type(matrixp)       :: Intmat(1)

!call LSHEADER(lupri,'II_get_regular_df_coulomb_mat')
call LSTIMER('START ',TSTART,TEND,LUPRI)
nbasis = SETTING%MOLECULE(1)%p%nbastREG
naux   = SETTING%MOLECULE(1)%p%nbastAUX
ndmat = 1
!(alpha|rho)
CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,lupri)
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,lupri)
IF(SETTING%SCHEME%FMM)call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
     & D%nrow,D%ncol,D%nrow,D%ncol,&
     &'DF-Aux','Empty','Regular','Regular','Coulomb','Contracted',.TRUE.)
CALL ls_freeDmatFromSetting(setting)
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,naux,1,1,1,ndmat)
call ls_jengine('DF-Aux','Empty','Regular','Regular','Coulomb','Regular','Contracted',&
     &          SETTING,LUPRI,LUERR)
call mem_alloc(galphafull,naux,1,1,1,ndmat)
IF(ndmat.EQ.1)then
   CALL retrieve_Output(lupri,setting,galphafull(:,:,1,1,1))
ELSE
   CALL retrieve_Output(lupri,setting,galphafull(:,:,1,1,:))
ENDIF
CALL ls_freeDmatFromSetting(setting)
call LSTIMER('GALPHA',TSTART,TEND,LUPRI)
!(alpha|beta)
call mem_alloc(alphabetafull,naux,1,naux,1,1)
call io_get_filename(Filename,'ALBE','DF-Aux','Empty','DF-Aux','Empty',0,0,&
     &'Coulomb','Contracted',.FALSE.,LUPRI,LUERR)
!Refine filename for FRAGMENT = .TRUE.
IF (io_file_exist(Filename,SETTING%IO).AND..NOT.(SETTING%numFragments.GT.1)) THEN
   call io_read(alphabetafull,naux,naux,1,1,1,Filename,SETTING%IO,LUPRI,LUERR)
ELSE
   call ls_DZERO(alphabetafull,naux*naux)
   call initIntegralOutputDims(setting%output,naux,1,naux,1,1)
   call ls_getIntegrals('DF-Aux','Empty',&
        &'DF-Aux','Empty',1,1,'Coulomb','Regular','Contracted',SETTING,LUPRI,LUERR)
   call retrieve_Output(lupri,setting,alphabetafull)
!Simen FRAGMENT
   IF (.NOT.(SETTING%numFragments.GT.1)) THEN
      call io_add_filename(SETTING%IO,Filename,LUPRI)
      call io_write(alphabetafull,naux,naux,1,1,1,Filename,SETTING%IO,LUPRI,LUERR)
   ENDIF
ENDIF
call LSTIMER('ALBE  ',TSTART,TEND,LUPRI)
!c_alpha = (alpha|beta)^-1 (beta|rho)
call mem_alloc(calphafull,naux,1,ndmat)
call DCOPY(naux*ndmat,galphafull,1,calphafull,1)
call DPOSV('U',naux,ndmat,alphabetafull,naux,calphafull,naux,info)
IF (info.ne.0) THEN
   WRITE(LUPRI,'(1X,A,I5)') 'DPOSV error in II_get_regular_df_coulomb_mat. Info =',info
   call LSQUIT('DPOSV error in II_get_regular_df_coulomb_mat',lupri)
ENDIF
call mem_dealloc(alphabetafull)
call mem_dealloc(galphafull)

! **** Write calpha to file in matrix format for storage (used for gradients)
call mat_init(calpha,naux,1)
call mat_set_from_full(calphafull(:,:,1),1.d0,calpha)
Filename = 'LSCALPHA'
IF (.NOT.io_file_exist(Filename,SETTING%IO)) THEN
  call io_add_filename(SETTING%IO,Filename,LUPRI)
ELSE
  IF (setting%scheme%incremental) THEN
    call mat_init(calpha_old,naux,1)
    call io_read_mat(calpha_old,Filename,Setting%IO,lupri,luerr)
    call mat_daxpy(1.d0,calpha_old,calpha)
    call mat_free(calpha_old)
  ENDIF
ENDIF
call io_write_mat(calpha,Filename,Setting%IO,lupri,luerr)
!call mat_print(CALPHA,1,naux,1,1,lupri) Stinne comment 14/9-2010!
call mat_free(calpha)
! **** End write

call LSTIMER('LINSOL',TSTART,TEND,LUPRI)
!Jfit_ab = (ab|rho_fit)
CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,lupri)
CALL ls_attachDmatToSetting(calphafull,ndmat,setting,'RHS',3,4,lupri)
IF(SETTING%SCHEME%FMM)call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
     & D%nrow,D%ncol,naux,1,&
     &'Regular','Regular','DF-Aux','Empty','Coulomb','Contracted',.TRUE.)
CALL ls_freeDmatFromSetting(setting)
CALL ls_attachDmatToSetting(calphafull,ndmat,setting,'RHS',3,4,lupri)
call mat_zero(F)              
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,F%nrow,F%ncol,1,1,1)
call ls_jengine('Regular','Regular','DF-Aux','Empty','Coulomb','Regular','Contracted',&
     &          SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F)
CALL ls_freeDmatFromSetting(setting)
call mem_dealloc(calphafull)
call LSTIMER('FIT-J ',TSTART,TEND,LUPRI)
END SUBROUTINE II_get_regular_df_coulomb_mat

!> \brief Calculates the coulomb matrix using overlap density fitting
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
SUBROUTINE II_get_overlap_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use ls_Integral_Interface
!use tco_Integral_Interface
use linsolvdf
IMPLICIT NONE
TYPE(MATRIX),target   :: D,F
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: ndmat,nbasis,naux
Real(realk),pointer :: g1alphafull(:,:,:,:,:)
!Real(realk),pointer :: g2alpha(:,:,:,:,:)
Real(realk),pointer :: calphafull(:,:,:)
TYPE(matrix),target :: g1alpha,calpha,g2alpha,F2
Real(realk)         :: TSTART,TEND
Character(10)       :: oper
Logical             :: Coulomb
Logical             :: Frag
type(matrixp)       :: Jmat(1),Dmat(1)
type(matrix),target :: tmpF,tmpD

!CALL LSHEADER(lupri,'II_get_overlap_df_coulomb_mat')
CALL LSTIMER('START ',TSTART,TEND,LUPRI)
oper='Overlap'
Coulomb = oper.EQ.'Coulomb'
ndmat = 1
nbasis = SETTING%MOLECULE(1)%p%nbastREG
naux   = SETTING%MOLECULE(1)%p%nbastAUX

!<alpha|w|rho>
CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,lupri)
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,lupri)
IF(SETTING%SCHEME%FMM.AND.Coulomb) call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
     & D%nrow,D%ncol,D%nrow,D%ncol,&
     & 'DF-Aux','Empty','Regular','Regular','Coulomb','Contracted',.TRUE.)
CALL ls_freeDmatFromSetting(setting)
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,lupri)
!CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,naux,1,1,1,1)
call ls_jengine('DF-Aux','Empty','Regular','Regular',oper,'Regular','Contracted',&
     &          SETTING,LUPRI,LUERR)
call mat_init(g1alpha,naux,1) 
CALL mat_ZERO(g1alpha)
CALL retrieve_Output(lupri,setting,g1alpha)
CALL ls_freeDmatFromSetting(setting)

call mem_alloc(g1alphafull,naux,1,1,1,ndmat)
call mat_to_full(g1alpha,1D0,g1alphafull(:,:,1,1,1))
call mat_free(g1alpha)
call mem_alloc(calphafull,naux,1,ndmat)
!|rho_fit) = <alpha|w|beta>^-1 <beta|w|rho>
call linsolv_df(calphafull,g1alphafull,'DF-Aux',oper,naux,ndmat,SETTING,LUPRI,LUERR)
call mem_dealloc(g1alphafull)
call mat_init(calpha,naux,1)
CALL mat_set_from_full(calphafull(:,:,1),1.D0,calpha)
call mem_dealloc(calphafull)
!Calculate multipole moments (lor read from file) if using FMM
CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,lupri)
CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,lupri)
IF(SETTING%SCHEME%FMM)call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
     & D%nrow,D%ncol,calpha%nrow,calpha%ncol,&
     &'Regular','Regular','DF-Aux','Empty','Coulomb','Contracted',.TRUE.)
CALL ls_freeDmatFromSetting(setting)
CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,lupri)
!Jfit_ab(1) = (ab|rho_fit)
CALL MAT_ZERO(F)
!Jmat(1)%p => F
!Dmat(1)%p => calpha
!CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,F%nrow,F%ncol,1,1,1)
call ls_jengine('Regular','Regular','DF-Aux','Empty','Coulomb','Regular','Contracted',&
     &          SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F)
CALL ls_freeDmatFromSetting(setting)

if(.NOT.Coulomb) then
   !(alpha|rho)
   CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,lupri)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,lupri)
   IF(SETTING%SCHEME%FMM) call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
     & D%nrow,D%ncol,calpha%nrow,calpha%ncol,&
     &'DF-Aux','Empty','Regular','Regular','Coulomb','Contracted',.TRUE.)
   CALL ls_freeDmatFromSetting(setting)
   CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,lupri)
   !Jmat,ndmat,
   call initIntegralOutputDims(setting%Output,naux,1,1,1,1)
   call ls_jengine('DF-Aux','Empty','Regular','Regular','Coulomb','Regular','Contracted',&
        &          SETTING,LUPRI,LUERR)
   call mat_init(g1alpha,naux,1) 
   call mat_zero(g1alpha) 
   CALL retrieve_Output(lupri,setting,g1alpha)
   CALL ls_freeDmatFromSetting(setting)
   !(alpha|rho_fit)
   CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,lupri)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,lupri)
   IF(SETTING%SCHEME%FMM)call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
        & D%nrow,D%ncol,calpha%nrow,calpha%ncol,&
        & 'DF-Aux','Empty','DF-Aux','Empty','Coulomb','Contracted',.TRUE.)
   CALL ls_freeDmatFromSetting(setting)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,lupri)
   call initIntegralOutputDims(setting%Output,naux,1,1,1,1)
   call ls_jengine('DF-Aux','Empty','DF-Aux','Empty','Coulomb','Regular','Contracted',&
        &          SETTING,LUPRI,LUERR)
   call mat_init(g2alpha,naux,1)
   CALL mat_ZERO(g2alpha)
   CALL retrieve_Output(lupri,setting,g2alpha)
   CALL ls_freeDmatFromSetting(setting)
   !(alpha|delta_rho)
   call mat_daxpy(-1.D0,g2alpha,g1alpha)
   call mat_free(g2alpha)
   !c_alpha = <alpha|w|beta>^-1 (beta|delta_rho)
   call mem_alloc(g1alphafull,naux,1,1,1,ndmat)
   call mat_to_full(g1alpha,1D0,g1alphafull(:,:,1,1,1))
   call mat_free(g1alpha)
   call mem_alloc(calphafull,naux,1,ndmat)
   call linsolv_df(calphafull,g1alphafull,'DF-Aux',oper,naux,ndmat,SETTING,LUPRI,LUERR)
   call mem_dealloc(g1alphafull)
   CALL mat_set_from_full(calphafull(:,:,1),1.D0,calpha)
   call mem_dealloc(calphafull)

   !Jfit_ab(2,3) = <ab|w|delta_rho>
   CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,lupri)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,lupri)
   IF(SETTING%SCHEME%FMM.AND.Coulomb)call ls_multipolemoment(LUPRI,LUERR,SETTING,&
        & nbasis,naux,D%nrow,D%ncol,calpha%nrow,calpha%ncol,&
        &'Regular','Regular','DF-Aux','Empty','Coulomb','Contracted',.TRUE.)
   CALL ls_freeDmatFromSetting(setting)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,lupri)
   call initIntegralOutputDims(setting%Output,nbasis,nbasis,1,1,1)
   call ls_jengine('Regular','Regular','DF-Aux','Empty',oper,'Regular','Contracted',&
        &          SETTING,LUPRI,LUERR)
   call mat_init(F2,nbasis,nbasis)
   CALL MAT_ZERO(F2)
   CALL retrieve_Output(lupri,setting,F2)
   CALL ls_freeDmatFromSetting(setting)
   call mat_daxpy(1.D0,F2,F)
   call mat_free(F2)
endif
call mat_free(calpha)
CALL LSTIMER('FIT-JO',TSTART,TEND,LUPRI)

END SUBROUTINE II_get_overlap_df_coulomb_mat

!> \brief Calculates the exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
SUBROUTINE II_get_exchange_mat(LUPRI,LUERR,SETTING,D,ndmat,Dsym,F)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use ls_Integral_Interface
IMPLICIT NONE
Integer               :: ndmat
TYPE(MATRIX)          :: D(ndmat)
TYPE(MATRIX)          :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
LOGICAL               :: Dsym
!
IF (SETTING%SCHEME%exchangeFactor.EQ.0.0d0) RETURN

CALL II_get_exchange_mat_regular(LUPRI,LUERR,SETTING,D,ndmat,Dsym,F)

END SUBROUTINE II_get_exchange_mat

!> \brief Calculates the exchange matrix using explicit 4 center 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
SUBROUTINE II_get_exchange_mat_regular(LUPRI,LUERR,SETTING,D,ndmat,Dsym,F)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use ls_Integral_Interface
IMPLICIT NONE
Integer               :: ndmat
TYPE(MATRIX),target   :: D(ndmat)
TYPE(MATRIX)          :: F(ndmat),K(ndmat),temp,temp2,K2
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
LOGICAL               :: Dsym
!
Integer             :: idmat,incdmat,nrow,ncol,ndmats
Real(realk),pointer :: Kfull(:,:,:,:,:)
Real(realk),pointer :: DfullLHS(:,:,:)
Real(realk),pointer :: DfullRHS(:,:,:)
Real(realk)         :: TS,TE,fac
CHARACTER(len=7)    :: Oper
integer :: nCalcInt,nCalcIntZero,nCalcIntZeroContrib

CALL LSTIMER('START ',TS,TE,LUPRI)
IF(matrix_type .EQ. mtype_unres_dense) THEN
  ndmats = 2*ndmat
ELSE
  ndmats = ndmat
ENDIF

!CALL LSHEADER(lupri,'II_get_exchange_mat')
setting%RHSdmat = .TRUE.
setting%nDmatRHS = ndmats

nrow = D(1)%nrow
ncol = D(1)%ncol
DO idmat=2,ndmat
  IF ((D(idmat)%nrow.NE.nrow).OR.(D(idmat)%ncol.NE.ncol)) THEN
     CALL LSQUIT('Error in calling II_get_exchange_mat. Mismatching D dimensions!',lupri)
  ENDIF
ENDDO

call mem_alloc(DfullRHS,nrow,ncol,ndmats)

DO idmat=1,ndmat
  IF(matrix_type .EQ. mtype_unres_dense)THEN
     CALL DCOPY(nrow*ncol,D(idmat)%elms, 1,DfullRHS(:,:,2*idmat-1),  1)
     CALL DCOPY(nrow*ncol,D(idmat)%elmsb,1,DfullRHS(:,:,2*idmat),1)
  ELSE !CLOSED_SHELL
     call mat_to_full(D(idmat),1D0,DfullRHS(:,:,idmat))
  ENDIF
ENDDO

IF (setting%scheme%daLinK) THEN
  call mem_alloc(DfullLHS,nrow,ncol,ndmats)
  setting%LHSdmat = .TRUE.
  setting%nDmatLHS = ndmats
  DO idmat=1,ndmat
    IF(matrix_type .EQ. mtype_unres_dense)THEN
       CALL DCOPY(nrow*ncol,D(idmat)%elms, 1,DfullLHS(:,:,2*idmat-1),  1)
       CALL DCOPY(nrow*ncol,D(idmat)%elmsb,1,DfullLHS(:,:,2*idmat),1)
    ELSE !CLOSED_SHELL
       call mat_to_full(D(idmat),1D0,DfullLHS(:,:,idmat))
    ENDIF
  ENDDO
  CALL ls_attachDmatToSetting(DfullLHS,ndmats,setting,'LHS',2,4,lupri)
ENDIF

CALL ls_attachDmatToSetting(DfullRHS,ndmats,setting,'RHS',1,3,lupri)

nrow = F(1)%nrow
ncol = F(1)%ncol
DO idmat=2,ndmat
  IF ((F(idmat)%nrow.NE.nrow).OR.(F(idmat)%ncol.NE.ncol)) THEN
     CALL LSQUIT('Error in calling II_get_exchange_mat. Mismatching F dimensions!',lupri)
  ENDIF
ENDDO

IF (SETTING%SCHEME%CAM) THEN
  Oper = 'Erfc   '   !Coulomb attenuated method
ELSE
  Oper = 'Coulomb'   !Regular Coulomb metric 
ENDIF
!Calculates the HF-exchange contribution
call initIntegralOutputDims(setting%Output,nrow,ncol,1,1,ndmats)
call ls_get_exchange_mat(Dsym,'Regular','Regular','Regular','Regular',&
     &                   Oper,'Regular','Contracted',SETTING,LUPRI,LUERR)

DO idmat=1,ndmat
  call mat_init(K(idmat),nrow,ncol)
ENDDO
CALL retrieve_Output(lupri,setting,K)

fac = 2.d0
IF(matrix_type .EQ. mtype_unres_dense)fac = 1.d0
DO IDMAT=1,ndmat
   WRITE(lupri,*)'The Exchange energy contribution ',-SETTING%SCHEME%exchangeFactor*fac*0.5d0*mat_dotproduct(D(idmat),K(idmat))
ENDDO
DO idmat=1,ndmat
!  call mat_scal(-SETTING%SCHEME%exchangeFactor,K(idmat))
  call mat_daxpy(-SETTING%SCHEME%exchangeFactor,K(idmat),F(idmat))
  call mat_free(K(idmat))
ENDDO

CALL ls_freeDfull(setting)
setting%RHSdmat = .FALSE.
setting%LHSdmat = .FALSE.
call mem_dealloc(DfullRHS)
IF (setting%scheme%daLinK)call mem_dealloc(DfullLHS)

IF(SETTING%SCHEME%DALINK .AND.SETTING%SCHEME%LINK)THEN
   CALL LSTIMER('DaLINK-Kbuild',TS,TE,LUPRI)
ELSEIF(SETTING%SCHEME%LINK)THEN
   CALL LSTIMER('LINK-Kbuild',TS,TE,LUPRI)
ELSE
   CALL LSTIMER('st-Kbuild',TS,TE,LUPRI)
ENDIF
END SUBROUTINE II_get_exchange_mat_regular

!> \brief Calculates the coulomb and exchange matrix using explicit 4 center 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the exchange matrix
!> \param ndmat the number of density matrices
SUBROUTINE II_get_coulomb_and_exchange_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use ls_Integral_Interface
IMPLICIT NONE
Integer               :: ndmat
TYPE(MATRIX)          :: D(ndmat),F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Real(realk)         :: TS,TE

CALL LSTIMER('START ',TS,TE,LUPRI)

!CALL LSHEADER(lupri,'II_get_coulomb_and_exchange_mat')
!ndmat = 1
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,lupri)
call initIntegralOutputDims(setting%Output,F(1)%nrow,F(1)%ncol,1,1,ndmat)
call ls_get_coulomb_and_exchange_mat('Regular','Regular','Regular','Regular',&
     &                   'Coulomb','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F)
CALL LSTIMER('reg-Fbuild',TS,TE,LUPRI)
END SUBROUTINE II_get_coulomb_and_exchange_mat

!> \brief Calculates the Fock matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
!> \param ndmat the number of density matrices
SUBROUTINE II_get_Fock_mat(LUPRI,LUERR,SETTING,D,Dsym,F,ndmat)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use lstiming
IMPLICIT NONE
integer               :: ndmat
TYPE(MATRIX)          :: D(ndmat),F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
LOGICAL               :: Dsym
!
Real(realk)  :: TS,TE
integer :: I
!CALL LSHEADER(lupri,'II_get_Fock_mat')
IF (SETTING%SCHEME%DENSFIT) THEN
   IF(matrix_type .EQ. mtype_unres_dense)&
        &CALL LSQUIT('Density fitting not implemented for unrestricted - spam Simen Reine',lupri)
   IF(ndmat.NE.1)CALL LSQUIT('For the time being the II_get_df_coulomb_mat is &
        &not implemted for more than 1 density matrix - spam Simen Reine',lupri)
   call II_get_df_coulomb_mat(LUPRI,LUERR,SETTING,D(1),F(1))
   call II_get_exchange_mat(LUPRI,LUERR,SETTING,D(1:1),1,Dsym,F(1:1))
ELSEIF (SETTING%SCHEME%JENGINE .OR. SETTING%SCHEME%LinK) THEN
   call II_get_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
   
   call II_get_exchange_mat(LUPRI,LUERR,SETTING,D,ndmat,Dsym,F)
   
ELSE
   CALL LSTIMER('START',TS,TE,LUPRI)
   !  call II_get_coulomb_and_exchange_mat(LUPRI,LUERR,SETTING,D,F)
   IF((matrix_type .EQ. mtype_unres_dense).OR. SETTING%SCHEME%CAM)THEN
      call II_get_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
      call II_get_exchange_mat(LUPRI,LUERR,SETTING,D,ndmat,Dsym,F)
   ELSE 
      call II_get_coulomb_and_exchange_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
   ENDIF
   CALL LSTIMER('reg-F  ',TS,TE,LUPRI)
ENDIF

END SUBROUTINE II_get_Fock_mat

!> \brief Calculates the 3 center eri matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param ALPHA the orbital index (ab|0alpha)
!> \param MAT the output matrix
!> 
!> Calculates the 3 center 2 electron repulsion 
!> integral matrix (ab|0alpha) for a given auxilliary 
!> basis function alpha
!>
SUBROUTINE II_GET_3CENTER_ERI_MAT(LUPRI,LUERR,SETTING,ALPHA,MAT)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use st_Integral_Intf
IMPLICIT NONE
TYPE(MATRIX)          :: MAT
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
Integer               :: alpha
!
!CALL LSHEADER(lupri,'II_get_3center_eri_mat')
call st_get_3Center_eri_mat(LUPRI,LUERR,SETTING,ALPHA,MAT)

END SUBROUTINE II_GET_3CENTER_ERI_MAT

!> \brief Calculates the 3 center eri matrix 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param ALPHA the orbital index (ab|0alpha)
!> \param DIM the dimension of orbital alpha
!> \param MAT the output matrix
!> 
!> Calculates the 3 center 2 electron repulsion 
!> integral matrix (ab|0alpha) for a given auxilliary 
!> basis function alpha
!>
SUBROUTINE II_GET_3CENTER_ERI_MATBATCH(LUPRI,LUERR,SETTING,ALPHA,DIM,MAT)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use st_Integral_Intf
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,alpha,DIM
TYPE(MATRIX)          :: MAT(DIM)
TYPE(LSSETTING)       :: SETTING

call st_get_3Center_eri_matbatch(LUPRI,LUERR,SETTING,ALPHA,DIM,MAT)

END SUBROUTINE II_GET_3CENTER_ERI_MATBATCH

!> \brief Calculates the 4 center shell eri dimensions
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param ALPHA the orbital index (alpha beta|delta gamma)
!> \param BETA the orbital index (alpha beta|delta gamma)
!> \param GAMMA the orbital index (alpha beta|delta gamma)
!> \param DELTA the orbital index (alpha beta|delta gamma)
!> \param DIM1 the dimension of orbital alpha
!> \param DIM2 the dimension of orbital beta
!> \param DIM3 the dimension of orbital delta
!> \param DIM4 the dimension of orbital gamma
!> 
!> Calculates the dimensions of the 4 center 2 electron 
!> repulsion integral matrix (alpha beta|delta gamma) for 
!> a given shell indeces alpha,beta,delta,gamma
!>
SUBROUTINE II_GET_4CENTERSHELL_ERI_DIM(LUPRI,LUERR,SETTING,ALPHA,BETA,GAMMA,DELTA,dim1,dim2,dim3,dim4)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use BUILDAOBATCH

IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
Integer               :: alpha,beta,gamma,delta
Integer               :: dim1,dim2,dim3,dim4
TYPE(AOITEM)          :: AO
Integer               :: ORB,iprint
write(lupri,*) 'II_GET_4CENTERSHELL_ERI_DIM'
IPRINT = SETTING%SCHEME%AOPRINT
CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,IPRINT,&
     &SETTING%MOLECULE(1)%p,&
     &SETTING%BASIS(1)%p%REGULAR,&
     &AO,'Default',.FALSE.,.FALSE.,alpha,dim1)
CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,IPRINT,&
     &SETTING%MOLECULE(1)%p,&
     &SETTING%BASIS(1)%p%REGULAR,&
     &AO,'Default',.FALSE.,.FALSE.,beta,dim2)
CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,IPRINT,&
     &SETTING%MOLECULE(1)%p,&
     &SETTING%BASIS(1)%p%REGULAR,&
     &AO,'Default',.FALSE.,.FALSE.,gamma,dim3)
CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,IPRINT,&
     &SETTING%MOLECULE(1)%p,&
     &SETTING%BASIS(1)%p%REGULAR,&
     &AO,'Default',.FALSE.,.FALSE.,delta,dim4)

END SUBROUTINE II_GET_4CENTERSHELL_ERI_DIM

!> \brief Calculates the 4 center shell eri tensor
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param ALPHA the orbital index (alpha beta|delta gamma)
!> \param BETA the orbital index (alpha beta|delta gamma)
!> \param GAMMA the orbital index (alpha beta|delta gamma)
!> \param DELTA the orbital index (alpha beta|delta gamma)
!> \param DIM1 the dimension of orbital alpha
!> \param DIM2 the dimension of orbital beta
!> \param DIM3 the dimension of orbital delta
!> \param DIM4 the dimension of orbital gamma
!> \param TENSOR the output tensor
!> 
!> Calculates the 4 center 2 electron repulsion integral matrix (ab|cd) so
!> collect a and b indices to one index and c and d to one index
!>
SUBROUTINE II_GET_4CENTERSHELL_ERI(LUPRI,LUERR,SETTING,ALPHA,BETA,GAMMA,DELTA,dim1,dim2,dim3,dim4,TENSOR)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use BUILDAOBATCH
  use integraldriver
  use Integralinfo
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,nbast
Integer               :: alpha,beta,gamma,delta
INTEGER               :: dim1,dim2,dim3,dim4
INTEGER               :: ORB,iprint
TYPE(INTEGRALINPUT)   :: INT_INPUT
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
REAL(REALK),target    :: TENSOR(dim1,dim2,dim3,dim4,1)
TYPE(LSSETTING)       :: SETTING
TYPE(AOITEM),target   :: AO1,AO2,AO3,AO4

call init_integral_input(INT_INPUT,SETTING)

CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
     &SETTING%MOLECULE(1)%p,&
     &SETTING%BASIS(1)%p%REGULAR,&
     &AO1,'Default',.FALSE.,.FALSE.,alpha,dim1)
CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
     &SETTING%MOLECULE(1)%p,&
     &SETTING%BASIS(1)%p%REGULAR,&
     &AO2,'Default',.FALSE.,.FALSE.,beta,dim2)
CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
     &SETTING%MOLECULE(1)%p,&
     &SETTING%BASIS(1)%p%REGULAR,&
     &AO3,'Default',.FALSE.,.FALSE.,gamma,dim3)
CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
     &SETTING%MOLECULE(1)%p,&
     &SETTING%BASIS(1)%p%REGULAR,&
     &AO4,'Default',.FALSE.,.FALSE.,delta,dim4)

INT_INPUT%operator = 'Coulomb'
INT_INPUT%AO(1)%p => AO1
INT_INPUT%AO(2)%p => AO2
INT_INPUT%AO(3)%p => AO3
INT_INPUT%AO(4)%p => AO4

Int_Output%ResultMat => TENSOR
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,&
     &INT_INPUT,INT_OUTPUT)

END SUBROUTINE II_GET_4CENTERSHELL_ERI

!> \brief Calculates the explicit 4 center eri tensor
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the explicit 4 center eri
!> \param DIM1 the dimension of orbital alpha
!> \param DIM2 the dimension of orbital beta
!> \param DIM3 the dimension of orbital delta
!> \param DIM4 the dimension of orbital gamma
SUBROUTINE II_get_4center_eri(LUPRI,LUERR,SETTING,outputintegral,dim1,dim2,dim3,dim4)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,dim1,dim2,dim3,dim4
REAL(REALK),target    :: outputintegral(dim1,dim2,dim3,dim4,1)
!
integer               :: usemat
Real(realk),pointer   :: integrals(:,:,:,:,:)
type(matrixp)         :: intmat(1)

call initIntegralOutputDims(setting%output,dim1,dim2,dim3,dim4,1)
CALL ls_getIntegrals('Regular','Regular','Regular','Regular',&
     &1,1,'Coulomb','Regular','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,outputintegral)

END SUBROUTINE II_get_4center_eri

!> \brief Calculates the dimensions of the shell indeces alpha
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param alpha the shell index 
!> \param DIM1 the dimension of shell index alpha
SUBROUTINE II_get_shellindex_dim(LUPRI,LUERR,SETTING,alpha,dim1)
use precision
use TYPEDEF  
use Matrix_module
use Matrix_Operations
use BUILDAOBATCH

IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
Integer               :: alpha,beta
Integer               :: dim1,dim2
TYPE(AOITEM)          :: AO
Integer               :: ORB,iprint

IPRINT = SETTING%SCHEME%AOPRINT
CALL BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,IPRINT,&
     &SETTING%MOLECULE(1)%p,&
     &SETTING%BASIS(1)%p%REGULAR,&
     &AO,'Default',.FALSE.,.FALSE.,alpha,dim1)

END SUBROUTINE II_get_shellindex_dim

!> \brief Calculates the decpacked explicit 4 center eri
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the output (full,full,batchB,batchC)
!> \param batchB batch index 
!> \param batchC batch index 
!> \param nbast full orbital dimension
!> \param dim2 the dimension of batch index 
!> \param dim3 the dimension of batch index 
SUBROUTINE II_GET_DECPACKED4CENTER_K_ERI(LUPRI,LUERR,SETTING,&
     &outputintegral,batchB,batchC,nbast,dim2,dim3)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast,dim2,dim3,batchB,batchC
REAL(REALK),target    :: outputintegral(nbast,nbast,dim2,dim3)
!
integer               :: I,J
type(matrixp)         :: intmat(1)
integer               :: nAO
logical,pointer       :: OLDsameMOLE(:,:)
nAO = setting%nAO
nullify(OLDsameMOLE)
allocate(OLDsameMOLE(nAO,nAO))
OLDsameMOLE = setting%sameMol
setting%batchindex(2)=batchB
setting%batchindex(3)=batchC
setting%batchdim(2)=dim2
setting%batchdim(3)=dim3
setting%sameMol(2,3)=batchB .EQ.batchC
setting%sameMol(3,2)=batchB .EQ.batchC
DO I=2,3
   J=1
   setting%sameMol(I,J)=.FALSE.
   setting%sameMol(J,I)=.FALSE.
   J=4
   setting%sameMol(I,J)=.FALSE.
   setting%sameMol(J,I)=.FALSE.
ENDDO
setting%scheme%decpacked = .true.
CALL II_get_4center_eri(LUPRI,LUERR,SETTING,outputintegral,nbast,nbast,dim2,dim3)
!back to normal
setting%batchindex(2)=0
setting%batchindex(3)=0
setting%batchdim(2)=0
setting%batchdim(3)=0
setting%sameFrag=.TRUE.
setting%scheme%decpacked = .false.
setting%sameMol = OLDsameMOLE
deallocate(OLDsameMOLE) 
nullify(OLDsameMOLE) 
END SUBROUTINE II_GET_DECPACKED4CENTER_K_ERI

!> \brief Calculates the decpacked explicit 4 center eri
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the output (full,full,batchC,batchD)
!> \param batchC batch index 
!> \param batchD batch index 
!> \param nbast full orbital dimension
!> \param dim3 the dimension of batch index 
!> \param dim4 the dimension of batch index 
SUBROUTINE II_GET_DECPACKED4CENTER_J_ERI(LUPRI,LUERR,SETTING,&
     &outputintegral,batchC,batchD,nbast,dim3,dim4)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast,dim3,dim4,batchC,batchD
REAL(REALK),target    :: outputintegral(nbast,nbast,dim3,dim4)
!
integer               :: I,J
type(matrixp)         :: intmat(1)
logical               :: psscreen
integer               :: nAO
logical,pointer       :: OLDsameMOLE(:,:)
nAO = setting%nAO
nullify(OLDsameMOLE)
allocate(OLDsameMOLE(nAO,nAO))
OLDsameMOLE = setting%sameMol

setting%batchindex(3)=batchC
setting%batchindex(4)=batchD
setting%batchdim(3)=dim3
setting%batchdim(4)=dim4
setting%sameMol(3,4)=batchC .EQ.batchD
setting%sameMol(4,3)=batchC .EQ.batchD
DO I=1,2
   DO J=3,4
      setting%sameMol(I,J)=.FALSE.
      setting%sameMol(J,I)=.FALSE.
   ENDDO
ENDDO
CALL II_get_4center_eri(LUPRI,LUERR,SETTING,outputintegral,nbast,nbast,dim3,dim4)
!back to normal
setting%batchindex(3)=0
setting%batchindex(4)=0
setting%batchdim(3)=0
setting%batchdim(4)=0
setting%sameFrag=.TRUE.
setting%sameMol = OLDsameMOLE
deallocate(OLDsameMOLE) 
nullify(OLDsameMOLE) 

END SUBROUTINE II_GET_DECPACKED4CENTER_J_ERI

!> \brief Calculates the diagonal explicit 4 center 2 eri
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param intType either primitive or contracted
!> \param AO1 options 'Regular','DF-Aux' or 'Empty'
!> \param AO2 options 'Regular','DF-Aux' or 'Empty'
!> \param setting Integral evalualtion settings
!> \param Gab the output matrix
SUBROUTINE II_get_2int_diag(LUPRI,LUERR,intType,AO1,AO2,SETTING,GAB)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(MATRIX),target   :: GAB
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
Character*(*)         :: AO1,AO2,intType
!
Real(realk),pointer :: integrals(:,:)
Integer             :: nbast,I,J
real(realk)         :: TS,TE
nbast = GAB%nrow
call mat_zero(GAB)
call mem_alloc(integrals,nbast,nbast)
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,1)
CALL ls_getScreenIntegrals('II_get_2int_diag',AO1,AO2,'Coulomb',intType,SETTING,LUPRI,LUERR,.TRUE.)
CALL retrieve_Output(lupri,setting,integrals)
DO J=1,nbast
   DO I=1,nbast
      integrals(I,J)=integrals(I,J)*integrals(I,J)
   ENDDO
ENDDO
CALL DCOPY(nbast*nbast,integrals,1,Gab%elms,1)
call mem_dealloc(integrals)
END SUBROUTINE II_get_2int_diag

!> \brief Calculates the 4 center 2 eri screening mat
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param Gab the output matrix
SUBROUTINE II_get_2int_ScreenMat(LUPRI,LUERR,SETTING,GAB)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(MATRIX),target   :: GAB
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: nbast
nbast = GAB%nrow
call mat_zero(GAB)
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,1)
CALL ls_getScreenIntegrals('II_get_2int_ScreenMat','Regular','Regular',&
     &'Coulomb','Contracted',SETTING,LUPRI,LUERR,.TRUE.)
CALL retrieve_Output(lupri,setting,GAB)

END SUBROUTINE II_get_2int_ScreenMat

!> \brief set the Incremental scheme
!> \author T. Kjaergaard
!> \date 2010
!> \param scheme the integral scheme
!> \param increm the increm logical
SUBROUTINE II_setIncremental(scheme,increm)
use precision
use TYPEDEF  
implicit none
TYPE(LSINTSCHEME),INTENT(INOUT) :: scheme
LOGICAL,INTENT(IN)              :: increm
scheme%incremental = increm
END SUBROUTINE II_setIncremental

!> \brief build BatchOrbitalInfo
!> \author T. Kjaergaard
!> \date 2010
!> \param setting Integral evalualtion settings
!> \param AO options 'Regular','DF-Aux' or 'Empty'
!> \param intType options contracted or primitive
!> \param nbast the number of basis functions
!> \param orbtobast for a given basis function it provides the batch index
!> \param nBatches the number of batches
!> \param lupri Default print unit
!> \param luerr Default error print unit
SUBROUTINE II_getBatchOrbitalInfo(Setting,AO,intType,nBast,orbToBatch,nBatches,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
  use ao_type
implicit none
TYPE(LSSETTING),intent(IN) :: Setting
Character*(*),intent(IN)   :: AO,intType
Integer,intent(IN)         :: nBast,lupri,luerr
Integer,intent(OUT)        :: orbToBatch(nBast)
Integer,intent(OUT)        :: nBatches
!
TYPE(BATCHORBITALINFO) :: BO
TYPE(AOITEM)           :: AObuild
Integer                :: nDim
!
CALL initBatchOrbitalInfo(BO,nBast)
CALL setAObatch(AObuild,0,nDim,AO,intType,'Default',Setting%scheme,Setting%fragment(1)%p,&
     &          setting%basis(1)%p,lupri,luerr)
CALL setBatchOrbitalInfo(BO,AObuild,lupri)
orbToBatch = BO%orbToBatch
nBatches   = BO%nBatches
CALL freeBatchOrbitalInfo(BO)
CALL free_aoitem(lupri,AObuild)
END SUBROUTINE II_getBatchOrbitalInfo

!> \brief Calculates the reorthonormalization gradient term
!> \author S. Reine
!> \date 2010-03-22
!> \param reOrtho The reorthonormalization gradient term
!> \param DFD The matrix product DFD = D F D
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_reorthoNormalization(reOrtho,DFDmat,ndmat,setting,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(OUT)       :: reOrtho(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrixp),intent(IN)      :: DFDmat(ndmat)
!
integer             :: nAtoms
real(realk)         :: OLDTHRESH

OLDTHRESH = SETTING%SCHEME%THRESHOLD
SETTING%SCHEME%THRESHOLD = SETTING%SCHEME%THRESHOLD*1.0D-5

CALL ls_attachDmatToSetting(DFDmat,ndmat,setting,'LHS',1,2,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

call initIntegralOutputDims(setting%output,3,nAtoms,1,1,1)
CALL ls_getIntegrals('Regular','Regular','Empty','Empty',ndmat,1,&
     &'Overlap','Gradient','Contracted',SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,reOrtho)
CALL ls_freeDmatFromSetting(setting)

SETTING%SCHEME%THRESHOLD = OLDTHRESH

END SUBROUTINE II_get_reorthoNormalization

!> \brief Check symetry and remove all anti-symmetric components (which will not contribute) to the Coulomb ERI gradient
!> \author S. Reine
!> \date 2010
!> \param DmatLHS The left-hand-side density matrix
!> \param DmatRHS The right-hand-side density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param DLHS The output left-hand-side density matrix
!> \param DRHS The output right-hand-side density matrix
!> \param symLHS The symmetry of LHS density matrices
!> \param symRHS The symmetry of RHS density matrices
!> \param nlhs The output number of LHS density matrices
!> \param nrhs The output number of RHS density matrices
SUBROUTINE II_symmetrize_dmats(DmatLHS,DmatRHS,ndlhs,ndrhs,DLHS,DRHS,symLHS,symRHS,nlhs,nrhs)
  use precision
  use TYPEDEF
  use Matrix_module
  use matrix_util, only : mat_get_isym,util_get_symm_part
implicit none
Integer,intent(in)          :: ndlhs,ndrhs
Type(matrixp),intent(IN)    :: DmatLHS(ndlhs),DmatRHS(ndrhs)
Integer,intent(inout)       :: symLHS(ndlhs),symRHS(ndrhs)
Integer,intent(inout)       :: nlhs,nrhs
Type(matrixp),intent(inout) :: DLHS(ndlhs),DRHS(ndrhs)
!
Real(realk), parameter :: thresh = 1.0d-14
Integer                :: idmat,isym

IF (ndlhs.NE.ndrhs) CALL lsquit('Erorr in II_symmetrize_dmats: ndlhs different from ndrhs.',-1)

nlhs = 0
DO idmat=1,ndlhs
  isym = mat_get_isym(DmatLHS(idmat)%p,thresh)
  symLHS(idmat) = isym
  If (isym.EQ.1) THEN
    nlhs = nlhs + 1
    DLHS(nlhs)%p => DmatLHS(idmat)%p
  ELSE IF  (isym.EQ.3) THEN
    nlhs = nlhs + 1
    ALLOCATE (DLHS(nlhs)%p)
    call mat_init(DLHS(nlhs)%p,DmatLHS(idmat)%p%nrow,DmatLHS(idmat)%p%ncol)
    call mat_assign(DLHS(nlhs)%p,DmatLHS(idmat)%p)
    call util_get_symm_part(DLHS(nlhs)%p)
  ELSE IF ((isym.EQ.2).OR.(isym.EQ.4)) THEN
    ! Do nothing (i.e. remove the component from the calculation)
  ELSE
    call lsquit('Error in II_symmetrize_dmats. isym wrong in lhs-loop!',-1)
  ENDIF
ENDDO

nrhs = 0
DO idmat=1,ndrhs
  isym = mat_get_isym(DmatRHS(idmat)%p,thresh)
  symRHS(idmat) = isym
  If (isym.EQ.1) THEN
    nrhs = nrhs + 1
    DRHS(nrhs)%p => DmatRHS(idmat)%p
  ELSE IF  (isym.EQ.3) THEN
    nrhs = nrhs + 1
    ALLOCATE (DRHS(nrhs)%p)
    call mat_init(DRHS(nrhs)%p,DmatRHS(idmat)%p%nrow,DmatRHS(idmat)%p%ncol)
    call mat_assign(DRHS(nrhs)%p,DmatRHS(idmat)%p)
    call util_get_symm_part(DRHS(nrhs)%p)
  ELSE IF ((isym.EQ.2).OR.(isym.EQ.4)) THEN
    ! Do nothing (i.e. remove the component from the calculation)
  ELSE
    call lsquit('Error in II_symmetrize_dmats. isym wrong in rhs-loop!',-1)
  ENDIF
ENDDO

DO idmat=1,ndrhs
  IF (symLHS(idmat).NE.symRHS(idmat)) &
 &  call lsquit('Erorr in II_symmetrize_dmats. Different matrix symmetries not implemented.',-1)
ENDDO

END SUBROUTINE II_symmetrize_dmats

!> \brief free the symmetrized dmats
!> \author S. Reine
!> \date 2010
!> \param DLHS The left-hand-side density matrix
!> \param DRHS The right-hand-side density matrix
!> \param symLHS The symmetry of LHS density matrices
!> \param symRHS The symmetry of RHS density matrices
!> \param nlhs The output number of LHS density matrices
!> \param nrhs The output number of RHS density matrices
SUBROUTINE II_free_symmetrized_dmats(DLHS,DRHS,symLHS,symRHS,ndlhs,ndrhs)
  use precision
  use TYPEDEF
  use Matrix_module
implicit none
Integer,intent(IN)          :: ndlhs,ndrhs
Integer,intent(IN)          :: symLHS(ndlhs),symRHS(ndrhs)
Type(matrixp),intent(inout) :: DLHS(ndlhs),DRHS(ndrhs)
!
Integer                :: nlhs,nrhs,idmat,isym

nlhs = 0
DO idmat=1,ndlhs
  isym = symLHS(idmat)
  If (isym.EQ.1) THEN
    nlhs = nlhs + 1
  ELSE IF  (isym.EQ.3) THEN
    nlhs = nlhs + 1
    call mat_free(DLHS(nlhs)%p)
    DEALLOCATE(DLHS(nlhs)%p)
  ENDIF
ENDDO

nrhs = 0
DO idmat=1,ndrhs
  isym = symRHS(idmat) 
  If (isym.EQ.1) THEN
    nrhs = nrhs + 1
  ELSE IF  (isym.EQ.3) THEN
    nrhs = nrhs + 1
    call mat_free(DRHS(nrhs)%p)
    DEALLOCATE (DRHS(nrhs)%p)
  ENDIF
ENDDO

END SUBROUTINE II_free_symmetrized_dmats

!> \brief split the dmats
!> \author S. Reine
!> \date 2010
!> \param DmatLHS The left-hand-side density matrix
!> \param DmatRHS The right-hand-side density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param DLHS The output left-hand-side density matrix
!> \param DRHS The output right-hand-side density matrix
!> \param symLHS The symmetry of LHS density matrices
!> \param symRHS The symmetry of RHS density matrices
!> \param nlhs The output number of LHS density matrices
!> \param nrhs The output number of RHS density matrices
SUBROUTINE II_split_dmats(DmatLHS,DmatRHS,ndlhs,ndrhs,DLHS,DRHS,symLHS,symRHS,nlhs,nrhs)
  use precision
  use TYPEDEF
  use Matrix_module
  use matrix_util, only : mat_get_isym,util_get_symm_part,util_get_antisymm_part
implicit none
Integer,intent(in)          :: ndlhs,ndrhs
Type(matrixp),intent(IN)    :: DmatLHS(ndlhs),DmatRHS(ndrhs)
Integer,intent(inout)       :: symLHS(ndlhs),symRHS(ndrhs)
Integer,intent(inout)       :: nlhs,nrhs
Type(matrixp),intent(inout) :: DLHS(2*ndlhs),DRHS(2*ndrhs)
!
Real(realk), parameter :: thresh = 1.0d-14
Integer                :: idmat,isym

IF (ndlhs.NE.ndrhs) CALL lsquit('Error in II_split_dmats: ndlhs different from ndrhs.',-1)

nlhs = 0
DO idmat=1,ndlhs
  isym = mat_get_isym(DmatLHS(idmat)%p,thresh)
  symLHS(idmat) = isym
  If (isym.EQ.1) THEN
    nlhs = nlhs + 1
    DLHS(nlhs)%p => DmatLHS(idmat)%p
  ELSE IF  (isym.EQ.3) THEN
    nlhs = nlhs + 2
    ALLOCATE (DLHS(nlhs-1)%p)
    ALLOCATE (DLHS(nlhs)%p)
    call mat_init(DLHS(nlhs-1)%p,DmatLHS(idmat)%p%nrow,DmatLHS(idmat)%p%ncol)
    call mat_init(DLHS(nlhs)%p,DmatLHS(idmat)%p%nrow,DmatLHS(idmat)%p%ncol)
    call mat_assign(DLHS(nlhs-1)%p,DmatLHS(idmat)%p)
    call util_get_symm_part(DLHS(nlhs-1)%p)
    call util_get_antisymm_part(DmatLHS(idmat)%p,DLHS(nlhs)%p)
  ELSE IF ((isym.EQ.2).OR.(isym.EQ.4)) THEN
    ! Do nothing (i.e. remove the component from the calculation)
  ELSE
    call lsquit('Error in II_split_dmats. isym wrong in lhs-loop!',-1)
  ENDIF
ENDDO

nrhs = 0
DO idmat=1,ndrhs
  isym = mat_get_isym(DmatRHS(idmat)%p,thresh)
  symRHS(idmat) = isym
  If (isym.EQ.1) THEN
    nrhs = nrhs + 1
    DRHS(nrhs)%p => DmatRHS(idmat)%p
  ELSE IF  (isym.EQ.3) THEN
    nrhs = nrhs + 2
    ALLOCATE (DRHS(nrhs-1)%p)
    ALLOCATE (DRHS(nrhs)%p)
    call mat_init(DRHS(nrhs-1)%p,DmatRHS(idmat)%p%nrow,DmatRHS(idmat)%p%ncol)
    call mat_init(DRHS(nrhs)%p,DmatRHS(idmat)%p%nrow,DmatRHS(idmat)%p%ncol)
    call mat_assign(DRHS(nrhs-1)%p,DmatRHS(idmat)%p)
    call util_get_symm_part(DRHS(nrhs-1)%p)
    call util_get_antisymm_part(DmatRHS(idmat)%p,DRHS(nrhs)%p)
  ELSE IF ((isym.EQ.2).OR.(isym.EQ.4)) THEN
    ! Do nothing (i.e. remove the component from the calculation)
  ELSE
    call lsquit('II_split_dmats in II_symmetrize_dmats. isym wrong in rhs-loop!',-1)
  ENDIF
ENDDO

DO idmat=1,ndrhs
  IF (symLHS(idmat).NE.symRHS(idmat)) &
 &  call lsquit('Error in II_split_dmats. Different matrix symmetries not implemented.',-1)
ENDDO

END SUBROUTINE II_split_dmats

!> \brief free the split dmats
!> \author S. Reine
!> \date 2010
!> \param DLHS The left-hand-side density matrix
!> \param DRHS The right-hand-side density matrix
!> \param symLHS The symmetry of LHS density matrices
!> \param symRHS The symmetry of RHS density matrices
!> \param nlhs The output number of LHS density matrices
!> \param nrhs The output number of RHS density matrices
SUBROUTINE II_free_split_dmats(DLHS,DRHS,symLHS,symRHS,ndlhs,ndrhs)
  use precision
  use TYPEDEF
  use Matrix_module
implicit none
Integer,intent(IN)          :: ndlhs,ndrhs
Integer,intent(IN)          :: symLHS(ndlhs),symRHS(ndrhs)
Type(matrixp),intent(inout) :: DLHS(2*ndlhs),DRHS(2*ndrhs)
!
Integer                :: nlhs,nrhs,idmat,isym

nlhs = 0
DO idmat=1,ndlhs
  isym = symLHS(idmat)
  If (isym.EQ.1) THEN
    nlhs = nlhs + 1
  ELSE IF  (isym.EQ.3) THEN
    nlhs = nlhs + 2
    call mat_free(DLHS(nlhs-1)%p)
    call mat_free(DLHS(nlhs)%p)
    DEALLOCATE(DLHS(nlhs-1)%p)
    DEALLOCATE(DLHS(nlhs)%p)
  ENDIF
ENDDO

nrhs = 0
DO idmat=1,ndrhs
  isym = symRHS(idmat) 
  If (isym.EQ.1) THEN
    nrhs = nrhs + 1
  ELSE IF  (isym.EQ.3) THEN
    nrhs = nrhs + 2
    call mat_free(DRHS(nrhs-1)%p)
    call mat_free(DRHS(nrhs)%p)
    DEALLOCATE (DRHS(nrhs-1)%p)
    DEALLOCATE (DRHS(nrhs)%p)
  ENDIF
ENDDO

END SUBROUTINE II_free_split_dmats

