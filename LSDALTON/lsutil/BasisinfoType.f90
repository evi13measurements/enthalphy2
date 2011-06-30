!> @file
!> basisinfo type module, contains also standard operation subroutines for this type
!> \brief basisinfo type module and associated subroutines for this type 
!> \author T. Kjaergaard
!> \date 2010
MODULE basis_type
 use lsmatrix_operations_dense
 use precision
 use AO_type
 INTEGER, PARAMETER      :: maxBASISsegment=25
 INTEGER, PARAMETER      :: maxBasisSetInLIB=10
 INTEGER, PARAMETER      :: maxNumberOfChargesinLIB=10

TYPE segment
INTEGER                 :: nrow
INTEGER                 :: ncol
real(realk),pointer     :: elms(:)     
real(realk),pointer     :: Exponents(:)  !length=nrow
END TYPE segment

TYPE SHELL
INTEGER                 :: nprim  !=#primitives
INTEGER                 :: norb   !=#orbitals
!TYPE(segment),pointer :: segment(:) 
TYPE(segment)           :: segment(maxBASISsegment) 
INTEGER                 :: nsegments
END TYPE SHELL

TYPE ATOMTYPEITEM
LOGICAL                  :: FAMILY !TYPE BASISSET
INTEGER                  :: nAngmom
INTEGER                  :: ToTnorb !total number of orbitals
INTEGER                  :: ToTnprim!total number of primitives
INTEGER                  :: Charge
!TYPE(ANGMOMITEM),pointer :: ANGMOM(:)
TYPE(SHELL)              :: SHELL(maxAOangmom)
Character(len=80)        :: NAME   
END TYPE ATOMTYPEITEM

TYPE BASISSETINFO
INTEGER                    :: natomtypes
TYPE(ATOMTYPEITEM),pointer :: ATOMTYPE(:)
INTEGER                    :: labelindex 
!  Labelindex indicate the ordering of the different types in ATOMTYPEITEM
!
!  labelindex=n for n>0 indicate MoleculeSpecific ordering which means that 
!  the molecule%ATOM(iatom)%IDtype(n) determines which ATOMTYPE the given atom has
!  labelindex=1 is for regularMoleculeSpecific ordering 
!  labelindex=2 is for auxiliaryMoleculeSpecific ordering
!  This is the case when ATOMBASIS is used in MOLECULE.INP
!   
!  labelindex=0 indicate charge based indexing which means that all atoms 
!  share the same basisset like 6-31G or cc-pVTZ - this is the case when 
!  BASIS is used in MOLECULE.INP, This means that the chargeindex is used
!  to acces the given ATOMTYPE
INTEGER,pointer            :: Chargeindex(:)!length nChargeindex=maxcharge
!only allocated when labelindex=0
INTEGER                    :: nChargeindex
INTEGER                    :: nbast,nprimbast
Character(len=9)           :: label !REGULAR for a regular basis
END TYPE BASISSETINFO

TYPE BASISINFO
TYPE(BASISSETINFO)        :: REGULAR
TYPE(BASISSETINFO)        :: HUCKEL
TYPE(BASISSETINFO)        :: AUXILIARY
END TYPE BASISINFO

TYPE BASIS_PT
TYPE(BASISINFO),pointer :: p
END TYPE BASIS_PT

TYPE BASINF
INTEGER             :: KMAX !number of shells for regular basis
! so an atom with 10s8p5d1f1g would have 25 shells (75 basisfunctions)
INTEGER,pointer     :: NHKT(:) !size kmax, angmom +1 for this shell
INTEGER,pointer     :: NUCO(:) !size kmax, Nr. primitives
INTEGER,pointer     :: NCENT(:) !size kmax, which atom the shell belongs to
INTEGER,pointer     :: NSTART(:)!size kmax, the shells start contracted Orbital 
INTEGER,pointer     :: JSTRT(:) !size kmax, the accumulated number of primitives
INTEGER,pointer     :: CCSTART(:)!size ushells, where the CCmatrix starts 
INTEGER,pointer     :: CCINDEX(:)!size kmax, where the CCmatrix starts 
real(realk),pointer :: X(:) !size natoms
real(realk),pointer :: Y(:) !size natoms
real(realk),pointer :: Z(:) !size natoms
real(realk),pointer :: PRIEXP(:) !size mxprim, primitive exponents 
TYPE(LSMatrix),pointer  :: CC(:) !size ushells, ContractionCoefficients for shell
real(realk),pointer :: RSHEL(:) !size kmax, radii of shells
real(realk),pointer :: CENT(:,:) !size 3,kmax x,y,z coord of shell
INTEGER             :: NHTYP,nAtoms,GRDONE,MXPRIM,Ushells 
INTEGER,pointer     :: CHARGE(:)!size natoms
END TYPE BASINF

CONTAINS
!#################################################################
!#
!# SUBROUTINES THAT ONLY AFFECT THE TYPES DEFINED IN THIS FILE
!# copy_basissetinfo
!# write_atomtypeitem_to_disk
!# write_basissetinfo_to_disk
!# write_basisinfo_to_disk
!# read_atomtypeitem_from_disk
!# read_basissetinfo_from_disk
!# read_basisinfo_from_disk
!# INIT_BASISSETINFO_Types
!# INIT_BASISSETINFO_ContractionMatrix
!# INIT_BASISSETINFO_elms
!# ALLOC_AND_TAKE_SUBBASISSETINFO
!# FREE_BASISSETINFO
!# LSMPI_ALLOC_BASISSETINFO
!# PRINT_BASISSETINFO
!################################################################

!> \brief write the atomtype structure to disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param Atomtype Contains the basis information for the atomtype, and type specifiers
!> \param lun logic unit number of file to write to
SUBROUTINE write_atomtypeitem_to_disk(lun,ATOMTYPE)
implicit none
TYPE(ATOMTYPEITEM),intent(in)    :: ATOMTYPE
INTEGER,intent(in)               :: lun
!
Integer :: I,J,K,nrow,ncol

write(lun)ATOMTYPE%nAngmom
write(lun)ATOMTYPE%FAMILY
write(lun)ATOMTYPE%ToTnorb
write(lun)ATOMTYPE%ToTnprim
write(lun)ATOMTYPE%Charge
DO I=1,ATOMTYPE%nAngmom
   WRITE(lun)ATOMTYPE%SHELL(I)%nprim
   WRITE(lun)ATOMTYPE%SHELL(I)%norb
   WRITE(lun)ATOMTYPE%SHELL(I)%nsegments
   DO J=1,ATOMTYPE%SHELL(I)%nsegments
      nrow = ATOMTYPE%SHELL(I)%segment(J)%nrow
      ncol = ATOMTYPE%SHELL(I)%segment(J)%ncol
      WRITE(lun)nrow
      WRITE(lun)ncol
      WRITE(lun)(ATOMTYPE%SHELL(I)%segment(J)%elms(K),K=1,nrow*ncol)
      WRITE(lun)(ATOMTYPE%SHELL(I)%segment(J)%Exponents(K),K=1,nrow)
   ENDDO
ENDDO
write(lun)ATOMTYPE%Name

END SUBROUTINE write_atomtypeitem_to_disk

!> \brief write the basissetinfo structure to disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basissetinfo Contains the basis information for given set
!> \param lun logic unit number of file to write to
SUBROUTINE write_basissetinfo_to_disk(lun,BASISSET)
implicit none
TYPE(BASISSETINFO),intent(in)    :: BASISSET
INTEGER,intent(in)               :: lun
!
INTEGER :: I

write(lun)BASISSET%natomtypes
DO I=1,BASISSET%natomtypes
   call write_atomtypeitem_to_disk(lun,BASISSET%ATOMTYPE(I))
ENDDO
write(lun)BASISSET%labelindex
write(lun)BASISSET%nChargeindex
DO I=1,BASISSET%nChargeindex
   write(lun)BASISSET%Chargeindex(I)
ENDDO
write(lun)BASISSET%nbast
write(lun)BASISSET%nprimbast
write(lun)BASISSET%label

END SUBROUTINE write_basissetinfo_to_disk

!> \brief write the basisinfo structure to disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information
!> \param lun logic unit number of file to write to
SUBROUTINE write_basisinfo_to_disk(lun,BASIS)
implicit none
TYPE(BASISINFO),intent(in)    :: BASIS
INTEGER,intent(in)            :: lun
!
INTEGER :: ZERO
ZERO=0
IF(BASIS%REGULAR%natomtypes.NE.0)THEN
   CALL WRITE_BASISSETINFO_TO_DISK(lun,BASIS%REGULAR)
ELSE
   write(lun)ZERO
ENDIF
IF(BASIS%HUCKEL%natomtypes.NE.0)THEN
   CALL WRITE_BASISSETINFO_TO_DISK(lun,BASIS%HUCKEL)
ELSE
   write(lun)ZERO
ENDIF
IF(BASIS%AUXILIARY%natomtypes.NE.0)THEN
   CALL WRITE_BASISSETINFO_TO_DISK(lun,BASIS%AUXILIARY)
ELSE
   write(lun)ZERO
ENDIF
END SUBROUTINE write_basisinfo_to_disk

!> \brief read the atomtype structure from disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param Atomtype Contains the basis information for the atomtype, and type specifiers
!> \param lun logic unit number of file to write to
SUBROUTINE read_atomtypeitem_from_disk(lun,ATOMTYPE)
implicit none
TYPE(ATOMTYPEITEM),intent(inout)    :: ATOMTYPE
INTEGER,intent(in)               :: lun
!
Integer :: I,J,nsize,K,nrow,ncol

read(lun)ATOMTYPE%nAngmom
read(lun)ATOMTYPE%FAMILY
read(lun)ATOMTYPE%ToTnorb
read(lun)ATOMTYPE%ToTnprim
read(lun)ATOMTYPE%Charge
DO I=1,ATOMTYPE%nAngmom
   READ(lun)ATOMTYPE%SHELL(I)%nprim
   READ(lun)ATOMTYPE%SHELL(I)%norb
   READ(lun)ATOMTYPE%SHELL(I)%nsegments
   DO J=1,ATOMTYPE%SHELL(I)%nsegments
      READ(lun)nrow
      READ(lun)ncol
      ATOMTYPE%SHELL(I)%segment(J)%nrow = nrow
      ATOMTYPE%SHELL(I)%segment(J)%ncol = ncol
      ALLOCATE(ATOMTYPE%SHELL(I)%segment(J)%elms(nrow*ncol))
      READ(lun)(ATOMTYPE%SHELL(I)%segment(J)%elms(K),K=1,nrow*ncol)
      ALLOCATE(ATOMTYPE%SHELL(I)%segment(J)%Exponents(nrow))
      READ(lun)(ATOMTYPE%SHELL(I)%segment(J)%Exponents(K),K=1,nrow)
   ENDDO
ENDDO
read(lun)ATOMTYPE%Name

END SUBROUTINE read_atomtypeitem_from_disk

!> \brief read the basissetinfo structure from disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basissetinfo Contains the basis information for given set
!> \param lun logic unit number of file to read from
SUBROUTINE read_basissetinfo_from_disk(lun,BASISSET)
implicit none
TYPE(BASISSETINFO),intent(inout)    :: BASISSET
INTEGER,intent(in)               :: lun
!
INTEGER :: I

read(lun)BASISSET%natomtypes
IF(BASISSET%natomtypes.GT.0)THEN
   NULLIFY(BASISSET%ATOMTYPE)
   ALLOCATE(BASISSET%ATOMTYPE(BASISSET%natomtypes))
   DO I=1,BASISSET%natomtypes
      call read_atomtypeitem_from_disk(lun,BASISSET%ATOMTYPE(I))
   ENDDO
   read(lun)BASISSET%labelindex
   read(lun)BASISSET%nChargeindex
   NULLIFY(BASISSET%Chargeindex)
   ALLOCATE(BASISSET%Chargeindex(BASISSET%nChargeindex))
   DO I=1,BASISSET%nChargeindex
      read(lun)BASISSET%Chargeindex(I)
   ENDDO
   read(lun)BASISSET%nbast
   read(lun)BASISSET%nprimbast
   read(lun)BASISSET%label
ENDIF

END SUBROUTINE read_basissetinfo_from_disk

!> \brief read the basisinfo structure from disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information
!> \param lun logic unit number of file to read from
SUBROUTINE read_basisinfo_from_disk(lun,BASIS)
implicit none
TYPE(BASISINFO),intent(inout)    :: BASIS
INTEGER,intent(in)            :: lun

CALL READ_BASISSETINFO_FROM_DISK(lun,BASIS%REGULAR)
CALL READ_BASISSETINFO_FROM_DISK(lun,BASIS%HUCKEL)
CALL READ_BASISSETINFO_FROM_DISK(lun,BASIS%AUXILIARY)

END SUBROUTINE read_basisinfo_from_disk

!> \brief initialise basissetinfo structure
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information for given set
!> \param natomtypes is number of atomtypes
SUBROUTINE init_basissetinfo_types(BasisInfo,natomtypes)
implicit none
TYPE(BASISSETINFO),intent(inout) :: BasisInfo
INTEGER,intent(in)               :: natomtypes

NULLIFY(BasisInfo%AtomType)
BasisInfo%natomtypes=natomtypes
ALLOCATE(BasisInfo%AtomType(natomtypes))

END SUBROUTINE init_basissetinfo_types

!> \brief initialise shell in basissetinfo structure
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information for given set
!> \param natomtype is number of atomtypes
!> \param nAngmom is the current angularmoment
!> \param segments is number of segments in the shell for this angular moment
SUBROUTINE init_basissetinfo_ContractionMatrix(BasisInfo,natomtype,&
                                                         &nAngmom,segments)
implicit none
TYPE(BASISSETINFO),intent(inout) :: BasisInfo
INTEGER,intent(in)            :: natomtype,nAngmom,segments

BasisInfo%AtomType(natomtype)%nAngmom=nAngmom
BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%nsegments=segments
IF(segments .GT. maxBASISsegment) THEN
   print*,'You have alot of segments you probably use a very&
        & large basisset which is not generally contracted. This&
        & is okay, but you need to increase the parameter&
        & maxBASISsegment in the lsutil/BasisinfoType.f90 file and recompile.&
        & Currently the parameter is set to ',maxBASISsegment
   print*,'You must increase the number in TYPE-DEF.f90 to at least ',segments
   print*,'Thomas Kjaergaard'
   CALL LSQUIT('Increase maxBASISsegment in TYPE-DEF.f90 file',-1)
ENDIF

END SUBROUTINE init_basissetinfo_ContractionMatrix

!> \brief initialise exponents and contractionmatrix in basissetinfo structure
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information for given set
!> \param natomtype is number of atomtypes
!> \param nAngmom is the current angularmoment
!> \param nsegments is number of segments in the shell for this angular moment
!> \param nrow is the number of exponents and the number of rows for the contraction matrix
!> \param ncol is the number of collums for the contraction matrix
SUBROUTINE init_basissetinfo_elms(BasisInfo,natomtype,nAngmom,&
                                                       &segments,nrow,ncol)
implicit none
TYPE(BASISSETINFO),intent(inout) :: BasisInfo
INTEGER,intent(in)            :: natomtype,nAngmom,segments,nrow,ncol
!
INTEGER            :: nsize,i
nsize=nrow*ncol
BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                                     &segment(segments)%nrow=nrow
BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                                      &segment(segments)%ncol=ncol
NULLIFY(BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                                          &segment(segments)%elms)
ALLOCATE(BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                                    &segment(segments)%elms(nsize))
do i = 1,nsize
  BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                   &segment(segments)%elms(i) = 0.0d0
enddo
NULLIFY(BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                       &segment(segments)%Exponents)
ALLOCATE(BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                       &segment(segments)%Exponents(nrow))

END SUBROUTINE init_basissetinfo_elms

!> \brief allocate and build a basissetinfo for a given type of a full basissetinfo  
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param oldbas is the old full BASISSETINFO
!> \param itype is the type in the full basissetinfo that is requested 
!> \param newbas is the new small BASISSETINFO which only contain 1 type, the itype in OLDBAS
SUBROUTINE alloc_and_take_subbasissetinfo(OLDBAS,itype,NEWBAS)
  implicit none
  TYPE(BASISSETINFO),intent(in)    :: OLDBAS
  TYPE(BASISSETINFO),intent(inout) :: NEWBAS
  INTEGER,intent(in)  :: itype
!
  INTEGER            :: nsegments,isegment,nbast,nprimbast
  INTEGER            :: nset,iset,ntype,iang,nrow,ncol,icharge
  
  nbast = 0
  nprimbast = 0
  NEWBAS%labelindex = OLDBAS%labelindex
  ALLOCATE(NEWBAS%ATOMTYPE(1))
  NEWBAS%ATOMTYPE(1)%name = OLDBAS%ATOMTYPE(itype)%name
  NEWBAS%nATOMTYPES = 1
  NEWBAS%ATOMTYPE(1) = OLDBAS%ATOMTYPE(itype)
  DO iang = 1,OLDBAS%ATOMTYPE(itype)%nangmom
     NEWBAS%ATOMTYPE(1)%SHELL(iang) = OLDBAS%ATOMTYPE(itype)%SHELL(iang)
     nsegments = OLDBAS%ATOMTYPE(itype)%SHELL(iang)%nsegments
     DO isegment = 1,nsegments
        nrow = OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)%nrow
        ncol = OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)%ncol
        nbast = nbast+ncol
        nprimbast = nprimbast+nrow
        NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment) = OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)
        ALLOCATE(NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%elms(nrow*ncol))
        ALLOCATE(NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%Exponents(nrow))
        NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%elms = &
             &OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)%elms
        NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%Exponents = &
             &OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)%Exponents
     ENDDO
  ENDDO
  NEWBAS%nbast = nbast
  NEWBAS%nprimbast = nprimbast
  IF(OLDBAS%labelindex .EQ. 0)THEN
     icharge = OLDBAS%ATOMTYPE(itype)%Charge
     ALLOCATE(NEWBAS%Chargeindex(icharge))
     NEWBAS%Chargeindex = 0
     NEWBAS%Chargeindex(icharge) = 1
  ENDIF
 
END SUBROUTINE alloc_and_take_subbasissetinfo

!> \brief copy basissetinfo to another basissetinfo  
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param oldbas is the original BASISSETINFO
!> \param newbas is the copied BASISSETINFO
SUBROUTINE copy_basissetinfo(OLDBAS,NEWBAS)
  implicit none
  TYPE(BASISSETINFO),intent(in)    :: OLDBAS
  TYPE(BASISSETINFO),intent(inout) :: NEWBAS
!
  INTEGER            :: I,J,K,nrow,ncol
  
  NEWBAS%natomtypes = OLDBAS%natomtypes
  NEWBAS%labelindex = OLDBAS%labelindex
  NEWBAS%nChargeindex = OLDBAS%nChargeindex
  NEWBAS%nbast = OLDBAS%nbast
  NEWBAS%nprimbast = OLDBAS%nprimbast
  NEWBAS%label = OLDBAS%label
  NULLIFY(NEWBAS%ATOMTYPE)
  ALLOCATE(NEWBAS%ATOMTYPE(NEWBAS%natomtypes))
  DO I = 1,NEWBAS%natomtypes
     NEWBAS%ATOMTYPE(I)%nAngmom = OLDBAS%ATOMTYPE(I)%nAngmom 
     NEWBAS%ATOMTYPE(I)%family = OLDBAS%ATOMTYPE(I)%family 
     NEWBAS%ATOMTYPE(I)%ToTnorb = OLDBAS%ATOMTYPE(I)%ToTnorb
     NEWBAS%ATOMTYPE(I)%ToTnprim = OLDBAS%ATOMTYPE(I)%ToTnprim
     NEWBAS%ATOMTYPE(I)%Charge = OLDBAS%ATOMTYPE(I)%Charge
     NEWBAS%ATOMTYPE(I)%NAME  = OLDBAS%ATOMTYPE(I)%NAME 
     DO J=1,NEWBAS%ATOMTYPE(I)%nAngmom
        NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim  = OLDBAS%ATOMTYPE(I)%SHELL(J)%nprim
        NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = OLDBAS%ATOMTYPE(I)%SHELL(J)%norb
        NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments  = OLDBAS%ATOMTYPE(I)%SHELL(J)%nsegments
        DO K=1,NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments
           nrow = OLDBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow
           ncol = OLDBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol
           NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow  = nrow
           NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol  = ncol
           NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms)
           ALLOCATE(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms(nrow*ncol))
           NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents)
           ALLOCATE(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents(nrow))
           NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms = &
                &OLDBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms
           NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents = &
                &OLDBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents
        ENDDO
     ENDDO
  ENDDO
  IF(NEWBAS%nChargeindex .NE. 0)THEN
     NULLIFY(NEWBAS%Chargeindex)
     ALLOCATE(NEWBAS%Chargeindex(NEWBAS%nChargeindex))
     DO I = 1,NEWBAS%nChargeindex
        NEWBAS%Chargeindex(I) = OLDBAS%Chargeindex(I)  
     ENDDO
  ENDIF
 
END SUBROUTINE copy_basissetinfo

!> \brief MPI braodcasts basissetinfo
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param bas is the BASISSETINFO
!> \param slave logical to tell if processor is slave
!> \param master integer of master process
SUBROUTINE mpicopy_basissetinfo(BAS,slave,master)
  use lsmpi_mod
  implicit none
  TYPE(BASISSETINFO),intent(inout)    :: BAS
  Integer            :: master
  logical            :: slave
!
  INTEGER            :: I,J,K,nrow,ncol
  
  call LS_MPIBCAST(BAS%natomtypes,Master)

  call LS_MPIBCAST(BAS%labelindex,Master)
  call LS_MPIBCAST(BAS%nChargeindex,Master)
  call LS_MPIBCAST(BAS%nbast,Master)
  call LS_MPIBCAST(BAS%nprimbast,Master)
  call LS_MPIBCAST(BAS%label,len(BAS%label),Master)

  IF(slave)THEN
     NULLIFY(BAS%ATOMTYPE)
     ALLOCATE(BAS%ATOMTYPE(BAS%natomtypes))
  ENDIF
  DO I = 1,BAS%natomtypes
     call LS_MPIBCAST(BAS%ATOMTYPE(I)%nAngmom,Master)
     call LS_MPIBCAST(BAS%ATOMTYPE(I)%family,Master)
     call LS_MPIBCAST(BAS%ATOMTYPE(I)%ToTnorb,Master)
     call LS_MPIBCAST(BAS%ATOMTYPE(I)%ToTnprim,Master)
     call LS_MPIBCAST(BAS%ATOMTYPE(I)%Charge,Master)
     call LS_MPIBCAST(BAS%ATOMTYPE(I)%NAME,LEN(BAS%ATOMTYPE(I)%NAME),Master)
     DO J=1,BAS%ATOMTYPE(I)%nAngmom
        call LS_MPIBCAST(BAS%ATOMTYPE(I)%SHELL(J)%nprim,Master)
        call LS_MPIBCAST(BAS%ATOMTYPE(I)%SHELL(J)%norb,Master)
        call LS_MPIBCAST(BAS%ATOMTYPE(I)%SHELL(J)%nsegments,Master)
        DO K=1,BAS%ATOMTYPE(I)%SHELL(J)%nsegments
           call LS_MPIBCAST(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow,Master)
           call LS_MPIBCAST(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol,Master)
           nrow = BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow
           ncol = BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol
           IF(SLAVE)THEN
              NULLIFY(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms)
              ALLOCATE(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms(nrow*ncol))
              NULLIFY(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents)
              ALLOCATE(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents(nrow))
           ENDIF
           call LS_MPIBCAST(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms,nrow*ncol,Master)
           call LS_MPIBCAST(BAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents,nrow,Master)
        ENDDO
     ENDDO
  ENDDO
  IF(BAS%nChargeindex .NE. 0)THEN
     IF(SLAVE)THEN
        NULLIFY(BAS%Chargeindex)
        ALLOCATE(BAS%Chargeindex(BAS%nChargeindex))
     ENDIF
     DO I = 1,BAS%nChargeindex
        call LS_MPIBCAST(BAS%Chargeindex(I),Master)
     ENDDO
  ENDIF
 
END SUBROUTINE mpicopy_basissetinfo

!> \brief deallocate BASISSETINFO
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param BasisInfo the BasissetInfo to be deallocated
SUBROUTINE free_basissetinfo(BasisInfo)
implicit none
TYPE(BASISSETINFO),intent(inout) :: BasisInfo
!
INTEGER  :: I,J,K,L,nsets,natomtypes,nAngmom,nsegments!,nsize,nrow,ncol
natomtypes=BASISINFO%natomtypes
IF(natomtypes.GT.0)THEN
 DO J=1,natomtypes
   nangmom=BASISINFO%ATOMTYPE(J)%nAngmom
   DO K=1,nAngmom
      nsegments=BASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments
      IF(nsegments .NE. 0)THEN
         DO L=1,nsegments
            if (.not.ASSOCIATED(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms)) then
               print*,'memory previously released!!'
               STOP 'Error in FREE_BASISSETINFO1 - memory previously released'
            endif
            DEALLOCATE(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms)
            NULLIFY(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms)  
            if (.not.ASSOCIATED(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents)) then
               print*,'memory previously released!!'
               STOP 'Error in FREE_BASISSETINFO2 - memory previously released'
            endif
            DEALLOCATE(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents)
            NULLIFY(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents)
         ENDDO
      ENDIF
   ENDDO
 ENDDO
 if ((natomtypes.NE.0).AND.(.not.ASSOCIATED(BASISINFO%ATOMTYPE))) then
   print*,'memory previously released!!'
   STOP 'Error in FREE_BASISSETINFO4 - memory previously released'
 endif
 DEALLOCATE(BASISINFO%ATOMTYPE)
 NULLIFY(BASISINFO%ATOMTYPE)
 IF(BASISINFO%Labelindex.EQ.0)THEN
   IF(.not.ASSOCIATED(BASISINFO%Chargeindex))THEN
     print*,'memory previously released!!'
     STOP 'Error in FREE_BASISSETINFO5  - memory previously released'
   ENDIF
   DEALLOCATE(BASISINFO%Chargeindex)   
   NULLIFY(BASISINFO%Chargeindex)   
 ENDIF
 BASISINFO%natomtypes=0
ENDIF
END SUBROUTINE free_basissetinfo

!> \brief allocate BASISSETINFO, from already set values in BASISSETINFO, used in MPI parallelization
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param BasisInfo the BasissetInfo to be allocated
SUBROUTINE lsmpi_alloc_basissetinfo(BASISINFO)
IMPLICIT NONE
TYPE(BASISSETINFO),intent(inout) :: BASISINFO
!
INTEGER            :: I,J,K,L,nrow,nsize

NULLIFY(BASISINFO%ATOMTYPE)
ALLOCATE(BASISINFO%ATOMTYPE(BASISINFO%natomtypes))   
DO J=1,BASISINFO%nAtomtypes
   !NO need to allocate SHELL
   DO K=1,BASISINFO%ATOMTYPE(J)%nAngmom
      !NO need to allocate segments
      DO L=1,BASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments
         nrow=BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%nrow
         nsize=nrow*BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%ncol
         NULLIFY(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms)
         NULLIFY(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents)
         ALLOCATE(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms(nSIZE))
         ALLOCATE(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents(nrow))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE lsmpi_alloc_basissetinfo

!> \brief print BASISSETINFO routine
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param lupri the logical unit number of the output file
!> \param BasisInfo the BasissetInfo to be printed
SUBROUTINE print_basissetinfo(LUPRI,BASISINFO)
implicit none
TYPE(BASISSETINFO),intent(in)  :: BASISINFO
INTEGER,intent(in)      :: LUPRI
!
INTEGER             :: nSets,nAngmom,natomtypes,I,J,K,L,Charge,segments
INTEGER             :: nprim,nContractions,nPrimitives,ncol
CHARACTER(len=1)    :: SPDFGH(10)=(/'S','P','D','F','G','H','I','J','K','L'/) 

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'BASISSETINFORMATION'
WRITE(LUPRI,*)'nbast',BASISINFO%nbast
natomtypes=BASISINFO%natomtypes
WRITE(LUPRI,*)'Number of different types of Atoms with this basisset:',natomtypes
DO J=1,natomtypes
   nAngmom=BASISINFO%ATOMTYPE(J)%nAngmom
   Charge=BASISINFO%ATOMTYPE(J)%Charge
   WRITE(LUPRI,*)'------------------------------------------------------------'
   WRITE(LUPRI,'(A8,I4,A14,I4)')' TYPE   :',J,' has charge ',Charge
   IF(BASISINFO%ATOMTYPE(J)%FAMILY) WRITE(LUPRI,'(A)')' This atomtype is a Family basis set type'
   WRITE(LUPRI,*)'------------------------------------------------------------'
   DO K=1,nAngmom
      segments=BASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments
      WRITE(LUPRI,*)'Number of ContractionCoefficient blocks =',segments
      IF(segments .NE. 0)THEN
         DO L=1,segments
            nPrimitives=BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%nrow
            ncol=BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%ncol
            IF (K.LE.10) THEN
               WRITE(LUPRI,'(2X,A1,A)')SPDFGH(K),'-TYPE FUNCTIONS'
            ELSE 
               WRITE(LUPRI,'(2X,I2,A)')K,'-TYPE FUNCTIONS'
            ENDIF
            WRITE(LUPRI,'(A22,I3,A22,I3)')'Exponents BLOCK=',L
            call OUTPUT(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%&
                 &Exponents, 1, nPrimitives, 1, 1, nPrimitives, 1, 1, LUPRI)
            WRITE(LUPRI,*)'segment BLOCK=',L
            CALL OUTPUT(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms,1,nPrimitives,1,ncol,nPrimitives,ncol,1,LUPRI)
         ENDDO
      ELSE
         IF (K.LE.10) THEN
            WRITE(LUPRI,'(2X,A3,A1,A)')'NO ',SPDFGH(K),'-TYPE FUNCTIONS'
         ELSE 
            WRITE(LUPRI,'(2X,A3,I2,A)')'NO ',K,'-TYPE FUNCTIONS'
         ENDIF
      ENDIF
   ENDDO
ENDDO

WRITE(LUPRI,*)' '

END SUBROUTINE print_basissetinfo

END MODULE basis_type

