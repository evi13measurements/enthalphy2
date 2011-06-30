!> @file
!> Contains general library routines for integral evaluation

!> \brief Calculates the one-electron AO-matrix
!> \author S. Reine
!> \date 2010-02-26
!> \param h1 The one-electron AO-matrix
!> \param nbast The number of oribtals
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_h1(h1,nbast,lupri,luerr)
  use configuration, only: init_integralconfig
!  use configuration, only: configitem, config_set_default_config, config_read_input
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
IMPLICIT NONE
Integer,intent(in)      :: nbast,lupri,luerr
Real(realk),intent(out) :: h1(nbast,nbast)
!
TYPE(MATRIX),target :: h
type(lsitem)        :: ls
Integer             :: mtype_save, nbas
!type(configItem)    :: config
TYPE(integralconfig)   :: integral

call init_integralconfig(integral,lupri)
!Initialize the ls-item from DALTON.INP and MOLECULE.INP
call ls_init(ls,lupri,luerr,nbas,integral,.false.,.false.)

mtype_save = matrix_type
matrix_type = mtype_dense
CALL mat_init(h,nbast,nbast)

CALL II_get_h1(lupri,luerr,ls%setting,h)

call dcopy(nbast*nbast,h%elms,1,h1,1)

CALL mat_free(h)
matrix_type = mtype_save

call ls_free(ls)

END SUBROUTINE LSlib_get_h1

!> \brief Calculates the gradient of the nuclear potential
!> \author S. Reine
!> \date 2010-02-26
!> \param nucGrad The nuclear potential gradient
!> \param nAtoms The number of atoms
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_nn_gradient(nucGrad,nAtoms,lupri,luerr)
use configuration, only: init_integralconfig
use TYPEDEF
use daltonInfo
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr
Real(realk),intent(OUT) :: nucGrad(3,nAtoms)
!
type(lsitem) :: ls
integer      :: nbast
!type(configItem)    :: config
TYPE(integralconfig)   :: integral

call init_integralconfig(integral,lupri)

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
!call config_set_default_config(config)
!call config_read_input(config,lupri,luerr,.FALSE.)
ls%input%dalton = integral!config%integral
call ls_init(ls,lupri,luerr,nbast,integral,.false.,.false.)

call II_get_nn_gradient(nucGrad,ls%setting,lupri,luerr)

call ls_free(ls)

END SUBROUTINE LSlib_get_nn_gradient

!> \brief Calculates the one-electron gradient contribution
!> \author S. Reine
!> \date 2010-03-22
!> \param oneGrad Nuclear-electronic-attraction gradient
!> \param D The density matrix
!> \param nbast The number of orbitals/basis functions
!> \param ndmat The number of density matrices
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_oneElectron_gradient(oneGrad,D,nbast,ndmat,nAtoms,lupri,luerr)
use configuration, only: init_integralconfig
use TYPEDEF
use daltonInfo
use matrix_operations
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast,ndmat
Real(realk),intent(OUT) :: oneGrad(3,nAtoms)
Real(realk),intent(IN)  :: D(nbast,nbast,ndmat)
!
type(lsitem)        :: ls
type(matrixp)       :: Dmat(ndmat)
type(matrix),target :: Dtarget(ndmat)
integer             :: idmat,nbasis
TYPE(integralconfig)   :: integral

call init_integralconfig(integral,lupri)
call ls_init(ls,lupri,luerr,nbasis,integral,.false.,.false.)

IF (ndmat.lt.1) CALL lsQUIT('ndmat<1 in LSlib_get_oneElectron_gradient',lupri)
IF (nAtoms.lt.ls%setting%molecule(1)%p%nAtoms) CALL lsQUIT('nAtoms inconsistency in LSlib_get_oneElectron_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_oneElectron_gradient',lupri)

!Copy density-matrices from D to Dmat (using Dtarget)
DO idmat=1,ndmat
  Dmat(idmat)%p => Dtarget(idmat)
  call mat_init(Dmat(idmat)%p,nbast,nbast)
  call mat_set_from_full(D(:,:,idmat),1.D0,Dmat(idmat)%p)
ENDDO

!Calculate the nuclear-electron attraction gradient
CALL II_get_oneElectron_gradient(oneGrad,Dmat,ndmat,ls%setting,lupri,luerr)

!Free density-matrices Dmat
DO idmat=1,ndmat
  call mat_free(Dmat(idmat)%p)
ENDDO

call ls_free(ls)

END SUBROUTINE LSlib_get_oneElectron_gradient

!> \brief Calculates the nuclear-electronic-attraction gradient
!> \author S. Reine
!> \date 2010-03-22
!> \param neGrad Nuclear-electronic-attraction gradient
!> \param D The density matrix
!> \param nbast The number of orbitals/basis functions
!> \param ndmat The number of density matrices
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_ne_gradient(neGrad,D,nbast,ndmat,nAtoms,lupri,luerr)
use configuration, only: init_integralconfig
use TYPEDEF
use daltonInfo
use matrix_operations
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast,ndmat
Real(realk),intent(OUT) :: neGrad(3,nAtoms)
Real(realk),intent(IN)  :: D(nbast,nbast,ndmat)
!
type(lsitem)        :: ls
type(matrixp)       :: Dmat(ndmat)
type(matrix),target :: Dtarget(ndmat)
integer             :: idmat,nbasis
!type(configItem)    :: config
TYPE(integralconfig)   :: integral

call init_integralconfig(integral,lupri)

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
!call config_set_default_config(config)
!call config_read_input(config,lupri,luerr,.FALSE.)
!ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbasis,integral,.false.,.false.)

IF (ndmat.lt.1) CALL lsQUIT('ndmat<1 in LSlib_get_ne_gradient',lupri)
IF (nbast.lt.1) CALL lsQUIT('nbast<1 in LSlib_get_ne_gradient',lupri)
IF (nAtoms.lt.1) CALL lsQUIT('nAtoms<1 in LSlib_get_ne_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_ne_gradient',lupri)

!Copy density-matrices from D to Dmat (using Dtarget)
DO idmat=1,ndmat
  Dmat(idmat)%p => Dtarget(idmat)
  call mat_init(Dmat(idmat)%p,nbast,nbast)
  call mat_set_from_full(D(:,:,idmat),1.D0,Dmat(idmat)%p)
ENDDO

!Calculate the nuclear-electron attraction gradient
CALL II_get_ne_gradient(neGrad,Dmat,ndmat,ls%setting,lupri,luerr)

!Free density-matrices Dmat
DO idmat=1,ndmat
  call mat_free(Dmat(idmat)%p)
ENDDO

call ls_free(ls)

END SUBROUTINE LSlib_get_ne_gradient

!> \brief Calculates the electron-electron repulsion gradient
!> \author S. Reine
!> \date 2010-04-06
!> \param eeGrad The exchange gradient contribution
!> \param DLHS The density matrix of the first electron (or left-hand-side)
!> \param DRHS The density matrix of the second electron (or right-hand-side)
!> \param nbast The number of orbitals/basis functions
!> \param ndlhs The number of density LHS matrices
!> \param ndrhs The number of density RHS matrices
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_twoElectron_gradient(eeGrad,DLHS,DRHS,nbast,ndlhs,ndrhs,nAtoms,lupri,luerr)
use configuration, only: configitem, config_set_default_config, config_read_input
use TYPEDEF
use daltonInfo
use matrix_operations
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast,ndlhs,ndrhs
Real(realk),intent(OUT) :: eeGrad(3,nAtoms)
Real(realk),intent(IN)  :: DLHS(nbast,nbast,ndlhs)
Real(realk),intent(IN)  :: DRHS(nbast,nbast,ndrhs)
!
type(lsitem)        :: ls
type(matrixp)       :: DmatLHS(ndlhs),DmatRHS(ndrhs)
type(matrix),target :: DtargetLHS(ndlhs),DtargetRHS(ndrhs)
integer             :: idmat,nbasis
type(configItem)    :: config
logical             :: coeff
Character(80)       :: Filename

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
call config_set_default_config(config)
call config_read_input(config,lupri,luerr,.FALSE.)
ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbasis,config%integral,.false.,.false.)

IF (ndlhs.lt.1) CALL lsQUIT('ndlhs<1 in LSlib_get_twoElectron_gradient',lupri)
IF (ndrhs.lt.1) CALL lsQUIT('ndrhs<1 in LSlib_get_twoElectron_gradient',lupri)
IF (nbast.lt.1) CALL lsQUIT('nbast<1 in LSlib_get_twoElectron_gradient',lupri)
IF (nAtoms.lt.1) CALL lsQUIT('nAtoms<1 in LSlib_get_twoElectron_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_twoElectron_gradient',lupri)

!Copy density-matrices from DLHS to DmatLHS (using DtargetLHS)
DO idmat=1,ndlhs
  DmatLHS(idmat)%p => DtargetLHS(idmat)
  call mat_init(DmatLHS(idmat)%p,nbast,nbast)
  call mat_set_from_full(DLHS(:,:,idmat),1.D0,DmatLHS(idmat)%p)
ENDDO
!Copy density-matrices from DRHS to DmatRHS (using DtargetRHS)
DO idmat=1,ndrhs
  DmatRHS(idmat)%p => DtargetRHS(idmat)
  call mat_init(DmatRHS(idmat)%p,nbast,nbast)
  call mat_set_from_full(DRHS(:,:,idmat),1.D0,DmatRHS(idmat)%p)
ENDDO

IF (ls%setting%scheme%densfit) THEN
  Filename = 'LSCALPHA'
  INQUIRE(file=Filename,exist=coeff)
  IF (.not.coeff) &
  &CALL lsQUIT('Error in LSlib_get_twoElectron_gradient. No fitting coefficients CALPHA on file')
  call io_add_filename(ls%setting%io,Filename,lupri)
ENDIF
!Calculate the nuclear-electron attraction gradient
CALL II_get_twoElectron_gradient(eeGrad,nAtoms,DmatLHS,DmatRHS,ndlhs,ndrhs,ls%setting,lupri,luerr)

!Free density-matrices DmatLHS and DmatRHS
DO idmat=1,ndlhs
  call mat_free(DmatLHS(idmat)%p)
ENDDO
DO idmat=1,ndrhs
  call mat_free(DmatRHS(idmat)%p)
ENDDO

call ls_free(ls)

END SUBROUTINE LSlib_get_twoElectron_gradient

!> \brief Calculates the exchange-correlation contribution to the molecular gradient
!> \author T. Kjaergaard
!> \date 2010-04-26
!> \param exGrad The exchange-correlation gradient contribution
!> \param DMAT The density matrix
!> \param nbast The number of orbitals/basis functions
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_xc_gradient(exGrad,Dmat,nbast,nAtoms,lupri,luerr)
use configuration, only: configitem, config_set_default_config, config_read_input
use TYPEDEF
use daltonInfo
use matrix_operations
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast
Real(realk),intent(OUT) :: exGrad(3,nAtoms)
Real(realk),intent(IN)  :: DMAT(nbast,nbast)
!
type(lsitem)        :: ls
type(matrix)        :: D
type(configItem)    :: config
integer             :: nbasis
!Initialize the ls-item from DALTON.INP and MOLECULE.INP
call config_set_default_config(config)
call config_read_input(config,lupri,luerr,.FALSE.)
ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbasis,config%integral,.TRUE.,.false.)

IF (nbast.lt.1) CALL lsQUIT('nbast<1 in LSlib_get_twoElectron_gradient',lupri)
IF (nAtoms.lt.1) CALL lsQUIT('nAtoms<1 in LSlib_get_twoElectron_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_twoElectron_gradient',lupri)

call mat_init(D,nbast,nbast)
call mat_set_from_full(DMAT,1.D0,D)

!Calculate the exchange-correlation gradient
CALL II_get_xc_geoderiv_molgrad(lupri,luerr,ls%setting,nbast,D,exGrad,natoms)

call mat_free(D)
call ls_free(ls)

END SUBROUTINE LSlib_get_xc_gradient

!> \brief Calculates the Coulomb contribution to the gradient
!> \author S. Reine
!> \date 2010-04-06
!> \param coulombGrad The Coulomb gradient contribution
!> \param DLHS The density matrix of the first electron (or left-hand-side)
!> \param DRHS The density matrix of the second electron (or right-hand-side)
!> \param nbast The number of orbitals/basis functions
!> \param ndlhs The number of density LHS matrices
!> \param ndrhs The number of density RHS matrices
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_J_gradient(coulombGrad,DLHS,DRHS,nbast,ndlhs,ndrhs,nAtoms,lupri,luerr)
use configuration, only: configitem, config_set_default_config, config_read_input
use TYPEDEF
use daltonInfo
use matrix_operations
implicit none
Integer,intent(IN)        :: nAtoms,lupri,luerr,nbast,ndlhs,ndrhs
Real(realk),intent(INOUT) :: coulombGrad(3,nAtoms)
Real(realk),intent(IN)    :: DLHS(nbast,nbast,ndlhs)
Real(realk),intent(IN)    :: DRHS(nbast,nbast,ndrhs)
!
type(lsitem)        :: ls
type(matrixp)       :: DmatLHS(ndlhs),DmatRHS(ndrhs)
type(matrix),target :: DtargetLHS(ndlhs),DtargetRHS(ndrhs)
integer             :: idmat,nbasis
type(configItem)    :: config
Character(80)       :: Filename
Logical             :: coeff

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
call config_set_default_config(config)
call config_read_input(config,lupri,luerr,.FALSE.)
ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbasis,config%integral,.false.,.false.)

IF (ndlhs.lt.1) CALL lsQUIT('ndlhs<1 in LSlib_get_J_gradient',lupri)
IF (ndrhs.lt.1) CALL lsQUIT('ndrhs<1 in LSlib_get_J_gradient',lupri)
IF (nbast.lt.1) CALL lsQUIT('nbast<1 in LSlib_get_J_gradient',lupri)
IF (nAtoms.lt.1) CALL lsQUIT('nAtoms<1 in LSlib_get_J_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_J_gradient',lupri)

!Copy density-matrices from DLHS to DmatLHS (using DtargetLHS)
DO idmat=1,ndlhs
  call mat_init(DtargetLHS(idmat),nbast,nbast)
  call mat_set_from_full(DLHS(:,:,idmat),1.D0,DtargetLHS(idmat))
  DmatLHS(idmat)%p => DtargetLHS(idmat)
ENDDO
!Copy density-matrices from DRHS to DmatRHS (using DtargetRHS)
DO idmat=1,ndrhs
  call mat_init(DtargetRHS(idmat),nbast,nbast)
  call mat_set_from_full(DRHS(:,:,idmat),1.D0,DtargetRHS(idmat))
  DmatRHS(idmat)%p => DtargetRHS(idmat)
ENDDO

IF (ls%setting%scheme%densfit) THEN
  Filename = 'LSCALPHA'
  INQUIRE(file=Filename,exist=coeff)
  IF (.not.coeff) &
  &CALL lsQUIT('Error in LSlib_get_twoElectron_gradient. No fitting coefficients CALPHA on file')
  call io_add_filename(ls%setting%io,Filename,lupri)
ENDIF
!Calculate the nuclear-electron attraction gradient
CALL II_get_J_gradient(coulombGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,ls%setting,lupri,luerr)

!Free density-matrices DmatLHS and DmatRHS
DO idmat=1,ndlhs
  call mat_free(DtargetLHS(idmat))
ENDDO
DO idmat=1,ndrhs
  call mat_free(DtargetRHS(idmat))
ENDDO

call ls_free(ls)

END SUBROUTINE LSlib_get_J_gradient

!> \brief Calculates the exchange contribution to the gradient
!> \author S. Reine
!> \date 2010-04-06
!> \param exchangeGrad The exchange gradient contribution
!> \param DLHS The density matrix of the first electron (or left-hand-side)
!> \param DRHS The density matrix of the second electron (or right-hand-side)
!> \param nbast The number of orbitals/basis functions
!> \param ndlhs The number of density LHS matrices
!> \param ndrhs The number of density RHS matrices
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_K_gradient(exchangeGrad,DLHS,DRHS,nbast,ndlhs,ndrhs,nAtoms,lupri,luerr)
use configuration, only: configitem, config_set_default_config, config_read_input
use TYPEDEF
use daltonInfo
use matrix_operations
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast,ndlhs,ndrhs
Real(realk),intent(OUT) :: exchangeGrad(3,nAtoms)
Real(realk),intent(IN)  :: DLHS(nbast,nbast,ndlhs)
Real(realk),intent(IN)  :: DRHS(nbast,nbast,ndrhs)
!
type(lsitem)        :: ls
type(matrixp)       :: DmatLHS(ndlhs),DmatRHS(ndrhs)
type(matrix),target :: DtargetLHS(ndlhs),DtargetRHS(ndrhs)
integer             :: idmat,nbasis
type(configItem)    :: config

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
call config_set_default_config(config)
call config_read_input(config,lupri,luerr,.FALSE.)
ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbasis,config%integral,.false.,.false.)

IF (nbast.lt.1) CALL lsQUIT('nbast<1 in LSlib_get_K_gradient',lupri)
IF (nAtoms.lt.1) CALL lsQUIT('nAtoms<1 in LSlib_get_K_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_K_gradient',lupri)

!Copy density-matrices from DLHS to DmatLHS (using DtargetLHS)
DO idmat=1,ndlhs
  DmatLHS(idmat)%p => DtargetLHS(idmat)
  call mat_init(DmatLHS(idmat)%p,nbast,nbast)
  call mat_set_from_full(DLHS(:,:,idmat),1.D0,DmatLHS(idmat)%p)
ENDDO
!Copy density-matrices from DRHS to DmatRHS (using DtargetRHS)
DO idmat=1,ndrhs
  DmatRHS(idmat)%p => DtargetRHS(idmat)
  call mat_init(DmatRHS(idmat)%p,nbast,nbast)
  call mat_set_from_full(DRHS(:,:,idmat),1.D0,DmatRHS(idmat)%p)
ENDDO

!Calculate the nuclear-electron attraction gradient
CALL II_get_K_gradient(exchangeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,ls%setting,lupri,luerr)

!Free density-matrices DmatLHS and DmatRHS
DO idmat=1,ndlhs
  call mat_free(DmatLHS(idmat)%p)
ENDDO
DO idmat=1,ndrhs
  call mat_free(DmatRHS(idmat)%p)
ENDDO

call ls_free(ls)

END SUBROUTINE LSlib_get_K_gradient

!> \brief Calculates the kinetic energy gradient
!> \author S. Reine
!> \date 2010-03-22
!> \param kinGrad Kinetic energy gradient
!> \param D The density matrix
!> \param nbast The number of orbitals/basis functions
!> \param ndmat The number of density matrices
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_kinetic_gradient(kinGrad,D,nbast,ndmat,nAtoms,lupri,luerr)
use configuration, only: init_integralconfig
use TYPEDEF
use daltonInfo
use matrix_operations
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast,ndmat
Real(realk),intent(OUT) :: kinGrad(3,nAtoms)
Real(realk),intent(IN)  :: D(nbast,nbast,ndmat)
!
type(lsitem)        :: ls
type(matrixp)       :: Dmat(ndmat)
type(matrix),target :: Dtarget(ndmat)
integer             :: idmat,nbasis
!type(configItem)    :: config
TYPE(integralconfig)   :: integral

call init_integralconfig(integral,lupri)

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
!call config_set_default_config(config)
!call config_read_input(config,lupri,luerr,.FALSE.)
!ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbasis,integral,.false.,.false.)

IF (ndmat.lt.1) CALL lsQUIT('ndmat<1 in LSlib_get_kinetic_gradient',lupri)
IF (nbast.lt.1) CALL lsQUIT('nbast<1 in LSlib_get_kinetic_gradient',lupri)
IF (nAtoms.lt.1) CALL lsQUIT('nAtoms<1 in LSlib_get_kinetic_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_kinetic_gradient',lupri)

!Copy density-matrices from D to Dmat (using Dtarget)
DO idmat=1,ndmat
  Dmat(idmat)%p => Dtarget(idmat)
  call mat_init(Dmat(idmat)%p,nbast,nbast)
  call mat_set_from_full(D(:,:,idmat),1.D0,Dmat(idmat)%p)
ENDDO

!Calculate the kinetic energy gradient
CALL II_get_kinetic_gradient(kinGrad,Dmat,ndmat,ls%setting,lupri,luerr)

!Free density-matrices Dmat
DO idmat=1,ndmat
  call mat_free(Dmat(idmat)%p)
ENDDO

call ls_free(ls)

END SUBROUTINE LSlib_get_kinetic_gradient

!> \brief Calculates the reorthonormalization gradient term
!> \author S. Reine
!> \date 2010-03-22
!> \param reOrtho The reorthonormalization gradient term
!> \param DFD The matrix product DFD = - D F D
!> \param nbast The number of orbital (basis functions)
!> \param ndmat The number of density matrices
!> \param nAtom The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error print unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_reorthoNormalization(reOrtho,DFD,nbast,ndmat,nAtoms,lupri,luerr)
use configuration, only: init_integralconfig
use TYPEDEF
use daltonInfo
use matrix_operations
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast,ndmat
Real(realk),intent(OUT) :: reOrtho(3,nAtoms)
Real(realk),intent(IN)  :: DFD(nbast,nbast,ndmat)
!
type(lsitem)        :: ls
type(matrixp)       :: DFDmat(ndmat)
type(matrix),target :: DFDtarget(ndmat)
integer             :: idmat,nbasis
!type(configItem)    :: config
TYPE(integralconfig)   :: integral

call init_integralconfig(integral,lupri)

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
!call config_set_default_config(config)
!call config_read_input(config,lupri,luerr,.FALSE.)
!ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbasis,integral,.false.,.false.)

IF (ndmat.lt.1) CALL lsQUIT('ndmat<1 in LSlib_get_reorthoNormalization',lupri)
IF (nbast.lt.1) CALL lsQUIT('nbast<1 in LSlib_get_reorthoNormalization',lupri)
IF (nAtoms.lt.1) CALL lsQUIT('nAtoms<1 in LSlib_get_reorthoNormalization',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_reorthoNormalization',lupri)

!Copy density-matrices from D to Dmat (using Dtarget)
DO idmat=1,ndmat
  call mat_init(DFDtarget(idmat),nbast,nbast)
  call mat_set_from_full(DFD(:,:,idmat),1.D0,DFDtarget(idmat))
  DFDmat(idmat)%p => DFDtarget(idmat)
ENDDO

!Calculate the nuclear-electron attraction gradient
CALL II_get_reorthoNormalization(reOrtho,DFDmat,ndmat,ls%setting,lupri,luerr)

!Free density-matrices Dmat
DO idmat=1,ndmat
  call mat_free(DFDtarget(idmat))
ENDDO

call ls_free(ls)

END SUBROUTINE LSlib_get_reorthoNormalization

!> \brief Sets up a list indicating which AO-batch each orbital is accosiated with
!> \author S. Reine
!> \date 18-03-2010
!> \param orbToBatch For each orbital index this vector returns the AO-batch it belongs to
!> \param nBast The number of orbitals/basis functions
!> \param nBatches The number of AO-batches (use in the integral code)
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> The AO-batches are use by the integral code, and resembles the traditional shells.
!> In a shell a set of primitive exponents gives rise to a contracted functions and a 
!> set of angular components for a fixed angular component. So say we have a p-orbital
!> block in basis-set file according to
!>     e1 c11 c21
!>     e2 c12 c22
!>     e3 c13 c23
!> The three primitive exponents e1-e3 give rise to two contracted functions, both with
!> p_x, p_y and p_z orital components - overall six orbitals. A batch is in most cases
!> the same as a shell, but may consist of different angular component block - all sharing
!> the same set of primitives.
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_getBatchOrbitalInfo(orbToBatch,nBast,nBatches,lupri,luerr)
  use precision
  use TYPEDEF  
  use ls_Integral_Interface
  use ao_type
  use configuration, only: init_integralconfig
  use daltonInfo
implicit none
Integer,intent(IN)         :: nBast,lupri,luerr
Integer,intent(OUT)        :: orbToBatch(nBast)
Integer,intent(OUT)        :: nBatches
!
TYPE(lsitem)     :: ls
!type(configItem) :: config
Integer :: nbasis
TYPE(integralconfig)   :: integral

call init_integralconfig(integral,lupri)

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
!call config_set_default_config(config)
!call config_read_input(config,lupri,luerr,.FALSE.)
!ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbasis,integral,.false.,.false.)
!
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_getBatchOrbitalInfo. Basis-function mismatch',lupri)
!
CALL II_getBatchOrbitalInfo(ls%setting,'Regular','Contracted',nBast,orbToBatch,nBatches,lupri,luerr)
!
call ls_free(ls)
!
END SUBROUTINE LSlib_getBatchOrbitalInfo

!> \brief Calculates the diagonal elements of the two-electron 4 center integrals - same as the Gab screening matrix squared. (ab|ab)
!> \author T. Kjaergaard
!> \date 2010-03-17
!> \param Gab the diagonal elements of the two-electron 4 center integrals
!> \param intType : Options 'Contracted' or 'Primitive'
!> \param AO1 : Options 'Regular','DF-Aux','Empty' denoting regular/aux/empty basis for center 1
!> \param AO2 : Options 'Regular','DF-Aux','Empty' denoting regular/aux/empty basis for center 2
!> \param nbast The number of oribtals (or basis functions)
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_2int_diag(Gab,intType,AO1,AO2,nbast,lupri,luerr)
  use configuration, only: init_integralconfig
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
IMPLICIT NONE
Integer,intent(in)      :: nbast,lupri,luerr
Real(realk),intent(out) :: gab(nbast,nbast)
Character*(*)           :: AO1,AO2,intType
!
TYPE(MATRIX),target :: tmp_gab
type(lsitem)        :: ls
Integer             :: mtype_save, nbas
!type(configItem)    :: config
TYPE(integralconfig)   :: integral

call init_integralconfig(integral,lupri)

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
!call config_set_default_config(config)
!call config_read_input(config,lupri,luerr,.FALSE.)
!ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbas,integral,.false.,.false.)

mtype_save = matrix_type
matrix_type = mtype_dense
CALL mat_init(tmp_Gab,nbast,nbast)

CALL II_get_2int_diag(lupri,luerr,intType,AO1,AO2,ls%setting,tmp_Gab)

call dcopy(nbast*nbast,tmp_Gab%elms,1,Gab,1)
CALL mat_free(tmp_Gab)
matrix_type = mtype_save
call ls_free(ls)

END SUBROUTINE LSlib_get_2int_diag

!> \brief Calculates the (ab|cd) with fixed a and b batchindexes so that the output would be a 4dim tensor with dim (dimAbatch,dimBbatch,fulldimC,fulldimD)
!> \author T. Kjaergaard
!> \date 2010-03-17
!> \param integrals the (ab|cd) tensor
!> \param batchA the requested A batchindex
!> \param batchB the requested B batchindex
!> \param dimA dimension of the requested  A batchindex
!> \param dimB dimension of the requested  B batchindex
!> \param nbast The number of oribtals (or basis functions) for c and d
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_ABresctricted_4CenterEri(integrals,batchA,batchB,dimA,dimB,nbast,lupri,luerr)
  use configuration, only: init_integralconfig
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
IMPLICIT NONE
Integer,intent(in)      :: nbast,dimA,dimB,lupri,luerr,batchA,batchB
Real(realk),intent(out) :: integrals(dimA,dimB,nbast,nbast)
!
TYPE(MATRIX),target :: tmp_gab
type(lsitem)        :: ls
Integer             :: mtype_save, nbas,I,J
!type(configItem)    :: config
TYPE(integralconfig)   :: integral

call init_integralconfig(integral,lupri)

!Initialize the ls-item from DALTON.INP and MOLECULE.INP
!call config_set_default_config(config)
!call config_read_input(config,lupri,luerr,.FALSE.)
!ls%input%dalton = config%integral
call ls_init(ls,lupri,luerr,nbas,integral,.false.,.false.)
ls%setting%batchindex(1)=batchA
ls%setting%batchindex(2)=batchB
ls%setting%batchdim(1)=dimA
ls%setting%batchdim(2)=dimB
ls%setting%sameMol(1,2)=batchA .EQ.batchB
ls%setting%sameMol(2,1)=batchA .EQ.batchB
DO I=1,2
   DO J=3,4
      ls%setting%sameMol(I,J)=.FALSE.
      ls%setting%sameMol(J,I)=.FALSE.
   ENDDO
ENDDO
CALL II_get_4center_eri(LUPRI,LUERR,ls%SETTING,integrals,dimA,dimB,nbast,nbast)
call ls_free(ls)

END SUBROUTINE LSlib_get_ABresctricted_4CenterEri

