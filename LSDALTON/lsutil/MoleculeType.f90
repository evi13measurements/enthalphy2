!> @file
!> molecule type module, contains also standard operation subroutines for this type
!> \author S. Reine and T.Kjaergaard
!> \date 2010-02-21
MODULE molecule_type
use precision
!*****************************************
!*
!* OBJECT CONTAINING INFORMATION ABOUT THE MOLECULE
!*
!*****************************************
TYPE ATOM
Integer           :: Isotope !Isotope number acc.to abundancy 
Character(len=4)  :: Name    !Name of atom
real(realk)       :: Mass    !Atomic mass 
real(realk)       :: CovRad  ! Covalent radius
real(realk)       :: Frag    ! Assigned fragment
real(realk)       :: CENTER(3)
Integer           :: Atomic_number !Atomic number
real(realk)       :: Charge        !Atomic Charge
integer           :: nbasis
Character(len=9)  :: basislabel(2) !REGULAR or AUXILIARY
INTEGER           :: Basisindex(2) !which set in the BASISSETLIBRARY
INTEGER           :: IDtype(2) !A unique identifier - identifying the type
LOGICAL           :: GHOST !true if basisset but no atoms 
! THE FOLLOWING ARE ADDED AFTER BUILD BASIS TO DO PARALLELIZATION
!INTEGER           :: SphOrbREG !Spherical orbitals for Regular basis
!INTEGER           :: CAROrbREG !Cartesian orbitals
INTEGER           :: nContOrbREG !# contracted orbitals
INTEGER           :: nPrimOrbREG !# primitives orbitals
!INTEGER           :: SphOrbAUX !Spherical orbitals for Aux bas
!INTEGER           :: CAROrbAUX !Cartesian orbitals
INTEGER           :: nContOrbAUX !# contracted orbitals
INTEGER           :: nPrimOrbAUX !# primitives orbitals
INTEGER           :: nContOrbHUC !# contracted orbitals for huckel
INTEGER           :: nPrimOrbHUC !# primitives orbitals for huckel
END TYPE ATOM

TYPE MOLECULEINFO
Character(len=22)    :: label
TYPE(ATOM), pointer  :: ATOM(:) !length = nAtomtypes
INTEGER              :: nAtoms
INTEGER              :: nelectrons
real(realk)          :: charge !molecular charge
! THE FOLLOWING ARE ADDED AFTER BUILD BASIS
INTEGER              :: nbastREG
INTEGER              :: nbastAUX
INTEGER              :: nbastHUC
INTEGER              :: nprimbastREG
INTEGER              :: nprimbastAUX
INTEGER              :: nprimbastHUC
END TYPE MOLECULEINFO

TYPE MOLECULE_PT
TYPE(MOLECULEINFO),pointer :: p
END TYPE MOLECULE_PT

TYPE MOLECULARORBITALINFO
INTEGER         :: nAtoms
INTEGER         :: nBastReg
INTEGER         :: nBastAux
INTEGER,POINTER :: numAtomicOrbitalsReg(:)
INTEGER,POINTER :: startAtomicOrbitalsReg(:)
INTEGER,POINTER :: endAtomicOrbitalsReg(:)
INTEGER,POINTER :: numAtomicOrbitalsAux(:)
INTEGER,POINTER :: startAtomicOrbitalsAux(:)
INTEGER,POINTER :: endAtomicOrbitalsAux(:)
END TYPE MOLECULARORBITALINFO

CONTAINS
!#################################################################
!#
!# SUBROUTINES THAT ONLY AFFECT THE TYPES DEFINED IN THIS FILE
!# init_moleculeinfo
!# free_Moleculeinfo
!# build_atomicmolecule
!# COPY_ATOM
!# COPY_MOLECULE
!# MPICOPY_ATOM
!# MPICOPY_MOLECULE
!# write_moleculeinfo_to_disk
!# read_moleculeinfo_from_disk
!# 
!################################################################

!> \brief write the atom structure to disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param iatom Contains the atom structure to be written
!> \param lun logic unit number of file to write to
SUBROUTINE write_atom_to_disk(lun,IATOM)
implicit none
TYPE(ATOM),intent(in)  :: IATOM
INTEGER,intent(in)             :: lun

write(lun)IATOM%Isotope
write(lun)IATOM%Name
write(lun)IATOM%CENTER(1)
write(lun)IATOM%CENTER(2)
write(lun)IATOM%CENTER(3)
write(lun)IATOM%Charge     
write(lun)IATOM%nbasis
write(lun)IATOM%basislabel(1)
write(lun)IATOM%basislabel(2)
write(lun)IATOM%Basisindex(1)
write(lun)IATOM%Basisindex(2)
write(lun)IATOM%IDtype(1) 
write(lun)IATOM%IDtype(2) 
write(lun)IATOM%GHOST
write(lun)IATOM%nContOrbREG
write(lun)IATOM%nPrimOrbREG
write(lun)IATOM%nContOrbAUX
write(lun)IATOM%nPrimOrbAUX 
write(lun)IATOM%nContOrbHUC
write(lun)IATOM%nPrimOrbHUC

END SUBROUTINE write_atom_to_disk

!> \brief read the atom structure to disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param iatom Contains the atom structure to be written
!> \param lun logic unit number of file to read from
SUBROUTINE read_atom_from_disk(lun,IATOM)
implicit none
TYPE(ATOM),intent(inout)  :: IATOM
INTEGER,intent(in)             :: lun

read(lun)IATOM%Isotope
read(lun)IATOM%Name
read(lun)IATOM%CENTER(1)
read(lun)IATOM%CENTER(2)
read(lun)IATOM%CENTER(3)
read(lun)IATOM%Charge     
read(lun)IATOM%nbasis
read(lun)IATOM%basislabel(1)
read(lun)IATOM%basislabel(2)
read(lun)IATOM%Basisindex(1)
read(lun)IATOM%Basisindex(2)
read(lun)IATOM%IDtype(1) 
read(lun)IATOM%IDtype(2) 
read(lun)IATOM%GHOST
read(lun)IATOM%nContOrbREG
read(lun)IATOM%nPrimOrbREG
read(lun)IATOM%nContOrbAUX
read(lun)IATOM%nPrimOrbAUX 
read(lun)IATOM%nContOrbHUC
read(lun)IATOM%nPrimOrbHUC

END SUBROUTINE read_atom_from_disk

!> \brief write the molecule structure to disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param Molecule Contains the information about the original molecule
!> \param lun logic unit number of file to write to
SUBROUTINE write_moleculeinfo_to_disk(lun,MOLECULE)
implicit none
TYPE(MOLECULEINFO),intent(in)  :: MOLECULE
INTEGER,intent(in)             :: lun
!
INTEGER                        :: I

write(lun)MOLECULE%nAtoms
write(lun)MOLECULE%label
DO I=1,MOLECULE%nAtoms
   call write_atom_to_disk(lun,MOLECULE%ATOM(I))
ENDDO
write(lun)MOLECULE%nelectrons
write(lun)MOLECULE%charge
write(lun)MOLECULE%nbastREG
write(lun)MOLECULE%nbastAUX
write(lun)MOLECULE%nbastHUC
write(lun)MOLECULE%nprimbastREG
write(lun)MOLECULE%nprimbastAUX
write(lun)MOLECULE%nprimbastHUC

END SUBROUTINE write_moleculeinfo_to_disk

!> \brief read the molecule structure from disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param Molecule Contains the information about the original molecule
!> \param lun logic unit number of file to write to
SUBROUTINE read_moleculeinfo_from_disk(lun,MOLECULE)
implicit none
TYPE(MOLECULEINFO),intent(inout)  :: MOLECULE
INTEGER,intent(in)                :: lun
!
INTEGER                           :: I

read(lun)MOLECULE%nAtoms
read(lun)MOLECULE%label
ALLOCATE(MOLECULE%ATOM(MOLECULE%nAtoms))
DO I=1,MOLECULE%nAtoms
   call read_atom_from_disk(lun,MOLECULE%ATOM(I))
ENDDO
read(lun)MOLECULE%nelectrons
read(lun)MOLECULE%charge
read(lun)MOLECULE%nbastREG
read(lun)MOLECULE%nbastAUX
read(lun)MOLECULE%nbastHUC
read(lun)MOLECULE%nprimbastREG
read(lun)MOLECULE%nprimbastAUX
read(lun)MOLECULE%nprimbastHUC

END SUBROUTINE read_moleculeinfo_from_disk

!> \brief Initializes the molecular info
!> \author S. Host
!> \date 2010-02-21
!> \param molecule Contains the information about the molecule
!> \param nAtoms The number of atoms
!> \param label  Label to uniqely specify the molecule (use for instance for screening matrix storage)
SUBROUTINE init_Moleculeinfo(Molecule,nAtoms,label)
implicit none
TYPE(MOLECULEINFO),intent(inout) :: MOLECULE
INTEGER,intent(IN)               :: nAtoms
Character(len=22),intent(IN)     :: label

NULLIFY(MOLECULE%ATOM)
MOLECULE%nAtoms = nAtoms
IF (nAtoms.GT.0) ALLOCATE(MOLECULE%ATOM(nAtoms))
MOLECULE%label = label

MOLECULE%nbastREG = 0
MOLECULE%nbastAUX = 0
MOLECULE%nbastHUC = 0
MOLECULE%nprimbastREG = 0
MOLECULE%nprimbastAUX = 0
MOLECULE%nprimbastHUC = 0
! call dens_stat_allocated_memory(nsize)

END SUBROUTINE init_Moleculeinfo

!> \brief Frees the molecular info
!> \author S. Host
!> \date 2010-02-21
!> \param molecule Contains the information about the molecule
SUBROUTINE free_Moleculeinfo(Molecule)
implicit none
TYPE(MOLECULEINFO) :: MOLECULE

if (MOLECULE%nAtoms.GT.0) then
  if (.not.ASSOCIATED(MOLECULE%ATOM)) then
    print*,'memory previously released!!'
    STOP 'Error in Molecule_free - memory previously released'
  endif
  DEALLOCATE(MOLECULE%ATOM)
  NULLIFY(MOLECULE%ATOM)
endif

END SUBROUTINE free_Moleculeinfo

!> \brief build moleculeinfo containing only 1 atom from full molecule
!> \author T. Kjaergaard
!> \date 2010-03-23
!> \param Molecule Contains the information about the original molecule
!> \param atomicMolecule the moleculeinfo to be built
!> \param iatom atom in full molecule
SUBROUTINE build_atomicmolecule(molecule,atomicmolecule,iatom,lupri)
IMPLICIT NONE
TYPE(MOLECULEINFO),intent(INOUT) :: atomicmolecule
TYPE(MOLECULEINFO),intent(IN)  :: MOLECULE
integer,intent(in) :: iatom,lupri
character(len=22) :: label

nullify(atomicmolecule%ATOM)
allocate(atomicmolecule%ATOM(1))
atomicmolecule%nAtoms=1
write(label,'(A6,I16)') 'GCATOM',iatom
atomicmolecule%label=label

atomicmolecule%nbastREG     = 0
atomicmolecule%nPrimbastREG = 0
atomicmolecule%nbastAUX     = 0
atomicmolecule%nPrimbastAUX = 0
atomicmolecule%nbastHUC     = 0
atomicmolecule%nPrimbastHUC = 0

call copy_atom(MOLECULE,Iatom,atomicmolecule,1,lupri)
atomicmolecule%nelectrons=INT(atomicmolecule%ATOM(1)%Charge)
atomicmolecule%charge=0

END SUBROUTINE build_atomicmolecule

!> \brief Copies a oldMOLECULE to newMOLECULE
!> \author T. Kjaergaard
!> \date 2010-05-31
!> \param oldMolecule Contains the information about the original molecule
!> \param newMolecule new molecule - same as oldMolecule afterwards
!> \param LUPRI LOGICAL UNIT NUMBER FOR OUTPUT
subroutine copy_molecule(oldMolecule,newMolecule,lupri)
implicit none
TYPE(MOLECULEINFO),intent(INOUT) :: newMOLECULE
TYPE(MOLECULEINFO),intent(IN)    :: oldMOLECULE
integer :: lupri
!
integer :: I

newMOLECULE%nAtoms = oldMOLECULE%nAtoms
newMOLECULE%label = oldMOLECULE%label
newMOLECULE%nbastREG = 0
newMOLECULE%nbastAUX = 0
newMOLECULE%nbastHUC = 0
newMOLECULE%nprimbastREG = 0
newMOLECULE%nprimbastAUX = 0
newMOLECULE%nprimbastHUC = 0

nullify(newMOLECULE%ATOM)
allocate(newMOLECULE%ATOM(newMOLECULE%nAtoms))
do I = 1, newMOLECULE%nAtoms
   call copy_atom(oldMOLECULE,I,newMOLECULE,I,LUPRI)
enddo
newMOLECULE%nelectrons = oldMOLECULE%nelectrons
newMOLECULE%charge = oldMOLECULE%charge

end subroutine copy_molecule

!> \brief Copies an atom from MOLECULE to FRAGMENT
!> \author S. Reine
!> \date 2010-02-21
!> \param Molecule Contains the information about the original molecule
!> \param I Atomic number of atom to be copied
!> \param Fragment Contains the information about the molecule copy
!> \param J Atomic number to be given for the copied atom
!> \param LUPRI LOGICAL UNIT NUMBER FOR OUTPUT
SUBROUTINE copy_atom(MOLECULE,I,FRAGMENT,J,LUPRI)
IMPLICIT NONE
INTEGER,intent(IN)               :: I,J,LUPRI
TYPE(MOLECULEINFO),intent(INOUT) :: FRAGMENT
TYPE(MOLECULEINFO),intent(IN)    :: MOLECULE
integer :: K

IF ((J.GT.FRAGMENT%nAtoms).OR.(J.LT.0)) THEN
  WRITE(LUPRI,'(1X,A,I4,A,I4)') 'Error in COPY_ATOM. Fragment atomic number ',J,&
     &                          ' is wrong - natoms = ',FRAGMENT%nAtoms
  CALL LSQUIT('Input in-consistency in COPY_ATOM',lupri)
ELSEIF ((I.GT.MOLECULE%nAtoms).OR.(I.LT.0)) THEN
  WRITE(LUPRI,'(1X,A,I4,A,I4)') 'Error in COPY_ATOM. Molecular atomic number ',I,&
     &                          ' is wrong - natoms = ',MOLECULE%nAtoms
  CALL LSQUIT('Input in-consistency in COPY_ATOM',lupri)
ENDIF

FRAGMENT%ATOM(J)%Isotope=MOLECULE%ATOM(I)%Isotope
FRAGMENT%ATOM(J)%Name=MOLECULE%ATOM(I)%Name
FRAGMENT%ATOM(J)%CENTER(1)=MOLECULE%ATOM(I)%CENTER(1)
FRAGMENT%ATOM(J)%CENTER(2)=MOLECULE%ATOM(I)%CENTER(2)
FRAGMENT%ATOM(J)%CENTER(3)=MOLECULE%ATOM(I)%CENTER(3)
FRAGMENT%ATOM(J)%Charge=MOLECULE%ATOM(I)%Charge
FRAGMENT%ATOM(J)%nbasis=MOLECULE%ATOM(I)%nbasis
do K = 1,MOLECULE%ATOM(I)%nbasis
   FRAGMENT%ATOM(J)%basislabel(K)=MOLECULE%ATOM(I)%basislabel(K)
   FRAGMENT%ATOM(J)%Basisindex(K)=MOLECULE%ATOM(I)%Basisindex(K)
   FRAGMENT%ATOM(J)%IDtype(K)=MOLECULE%ATOM(I)%IDtype(K)
enddo
FRAGMENT%ATOM(J)%GHOST=MOLECULE%ATOM(I)%GHOST
FRAGMENT%ATOM(J)%nContOrbREG=MOLECULE%ATOM(I)%nContOrbREG
FRAGMENT%ATOM(J)%nPrimOrbREG=MOLECULE%ATOM(I)%nPrimOrbREG
FRAGMENT%ATOM(J)%nContOrbAUX=MOLECULE%ATOM(I)%nContOrbAUX
FRAGMENT%ATOM(J)%nPrimOrbAUX=MOLECULE%ATOM(I)%nPrimOrbAUX
FRAGMENT%ATOM(J)%nContOrbHUC=MOLECULE%ATOM(I)%nContOrbHUC
FRAGMENT%ATOM(J)%nPrimOrbHUC=MOLECULE%ATOM(I)%nPrimOrbHUC

FRAGMENT%nbastREG     = FRAGMENT%nbastREG     + MOLECULE%ATOM(I)%nContOrbREG
FRAGMENT%nPrimbastREG = FRAGMENT%nPrimbastREG + MOLECULE%ATOM(I)%nPrimOrbREG
FRAGMENT%nbastAUX     = FRAGMENT%nbastAUX     + MOLECULE%ATOM(I)%nContOrbAUX
FRAGMENT%nPrimbastAUX = FRAGMENT%nPrimbastAUX + MOLECULE%ATOM(I)%nPrimOrbAUX
FRAGMENT%nbastHUC     = FRAGMENT%nbastHUC     + MOLECULE%ATOM(I)%nContOrbHUC
FRAGMENT%nPrimbastHUC = FRAGMENT%nPrimbastHUC + MOLECULE%ATOM(I)%nPrimOrbHUC

END SUBROUTINE copy_atom

!> \brief MPI Copies(Broadcasts) a MOLECULE
!> \author T. Kjaergaard
!> \date 2010-05-31
!> \param Molecule Contains the information about the molecule
!> \param Slave if this processor is a slave 
!> \param Master the integer for the master process
subroutine mpicopy_molecule(Molecule,SLAVE,Master)
use lsmpi_mod
implicit none
TYPE(MOLECULEINFO),intent(INOUT) :: MOLECULE
logical :: Slave
integer :: master
!
integer :: I

call LS_MPIBCAST(MOLECULE%nAtoms,Master)

IF(SLAVE)THEN
   MOLECULE%nbastREG = 0
   MOLECULE%nbastAUX = 0
   MOLECULE%nbastHUC = 0
   MOLECULE%nprimbastREG = 0
   MOLECULE%nprimbastAUX = 0
   MOLECULE%nprimbastHUC = 0
   nullify(MOLECULE%ATOM)
   allocate(MOLECULE%ATOM(MOLECULE%nAtoms))
ENDIF
do I = 1, MOLECULE%nAtoms
   call mpicopy_atom(MOLECULE,I,Slave,Master)
enddo

call LS_MPIBCAST(MOLECULE%nelectrons,Master)
call LS_MPIBCAST(MOLECULE%charge,Master)

end subroutine mpicopy_molecule

!> \brief MPI Copies(Broadcasts) an atom from MOLECULE
!> \author S. Reine
!> \date 2010-02-21
!> \param Molecule Contains the information about the original molecule
!> \param I Atomic number of atom to be copied
!> \param Slave if this processor is a slave 
!> \param Master the integer for the master process
SUBROUTINE mpicopy_atom(MOLECULE,I,Slave,Master)
use lsmpi_mod
IMPLICIT NONE
INTEGER,intent(IN)               :: I,Master
TYPE(MOLECULEINFO),intent(INOUT)    :: MOLECULE
Logical :: Slave
!
integer :: K

call LS_MPIBCAST(MOLECULE%ATOM(I)%Isotope,Master)
call LS_MPIBCAST(MOLECULE%ATOM(I)%Name,len(MOLECULE%ATOM(I)%Name),Master)

call LS_MPIBCAST(MOLECULE%ATOM(I)%CENTER,3,Master)
call LS_MPIBCAST(MOLECULE%ATOM(I)%Charge,Master)
call LS_MPIBCAST(MOLECULE%ATOM(I)%nbasis,Master)
do K = 1,MOLECULE%ATOM(I)%nbasis
   call LS_MPIBCAST(MOLECULE%ATOM(I)%basislabel(K),len(MOLECULE%ATOM(I)%basislabel(K)),Master)
   call LS_MPIBCAST(MOLECULE%ATOM(I)%basisindex(K),Master)
   call LS_MPIBCAST(MOLECULE%ATOM(I)%IDtype(K),Master)
enddo
call LS_MPIBCAST(MOLECULE%ATOM(I)%GHOST,Master)

call LS_MPIBCAST(MOLECULE%ATOM(I)%nContOrbREG,Master)
call LS_MPIBCAST(MOLECULE%ATOM(I)%nPrimOrbREG,Master)
call LS_MPIBCAST(MOLECULE%ATOM(I)%nContOrbAUX,Master)
call LS_MPIBCAST(MOLECULE%ATOM(I)%nPrimOrbAUX,Master)
call LS_MPIBCAST(MOLECULE%ATOM(I)%nContOrbHUC,Master)
call LS_MPIBCAST(MOLECULE%ATOM(I)%nPrimOrbHUC,Master)
IF(SLAVE)THEN
   MOLECULE%nbastREG     = MOLECULE%nbastREG     + MOLECULE%ATOM(I)%nContOrbREG
   MOLECULE%nPrimbastREG = MOLECULE%nPrimbastREG + MOLECULE%ATOM(I)%nPrimOrbREG
   MOLECULE%nbastAUX     = MOLECULE%nbastAUX     + MOLECULE%ATOM(I)%nContOrbAUX
   MOLECULE%nPrimbastAUX = MOLECULE%nPrimbastAUX + MOLECULE%ATOM(I)%nPrimOrbAUX
   MOLECULE%nbastHUC     = MOLECULE%nbastHUC     + MOLECULE%ATOM(I)%nContOrbHUC
   MOLECULE%nPrimbastHUC = MOLECULE%nPrimbastHUC + MOLECULE%ATOM(I)%nPrimOrbHUC
ENDIF
END SUBROUTINE mpicopy_atom

END MODULE molecule_type
