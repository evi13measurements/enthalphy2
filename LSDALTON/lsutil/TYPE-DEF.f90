!> @file 
!> contains many structure and associated subroutine
MODULE TYPEDEF
 use files
 use precision
! use matrix_module
 use molecule_type
 use basis_type
! use aobatch_type
 use io
 use integral_type
 use ao_type
 use lsmatrix_type
 use matrix_module
 use integralOutput_type

INTERFACE retrieve_output
   MODULE PROCEDURE retrieve_output_mat_single, &
        & retrieve_output_mat_array,&
        & retrieve_output_5dim, retrieve_output_2dim,&
        & retrieve_output_3dim, retrieve_output_lstensor
END INTERFACE

INTERFACE II_setMolecules
  MODULE PROCEDURE II_setMolecules_4
  MODULE PROCEDURE II_setMolecules_2
  MODULE PROCEDURE II_setMolecules_2_1
  MODULE PROCEDURE II_setMolecules_1_1_1
  MODULE PROCEDURE II_setMolecules_1_1_1_1
END INTERFACE II_setMolecules

TYPE INTEGERP
INTEGER,pointer :: p
END TYPE INTEGERP

!*****************************************
!*
!* OBJECT CONTAINING INFORMATION ABOUT THE CALCULATION
!* THE DALTON INPUT FILE
!*
!*****************************************
TYPE integralconfig
!PARAMETERS FROM **INTEGRALS   DECLERATION
LOGICAL  :: UNRES
LOGICAL  :: CFG_LSDALTON !.RUN LINSCA SPECIFIED 
LOGICAL  :: TRILEVEL
LOGICAL  :: DOPASS
LOGICAL  :: DENSFIT
LOGICAL  :: LINSCA
!*DENSFIT PARAMETERS

!*LINSCA PRINT PARAMETERS
INTEGER  :: LINSCAPRINT
INTEGER  :: AOPRINT
INTEGER  :: MOLPRINT
INTEGER  :: INTPRINT
INTEGER  :: BASPRINT
LOGICAL  :: PRINTATOMCOORD
!*MAIN INTEGRAL PARAMETERS
LOGICAL  :: JENGINE
LOGICAL  :: LOCALLINK
REAL(REALK)  :: LOCALLINKmulthr
LOGICAL  :: LOCALLINKDcont
REAL(REALK)  :: LOCALLINKDthr
LOGICAL  :: LOCALLINKsimmul
INTEGER  :: LOCALLINKoption
LOGICAL  :: LOCALLINKincrem
INTEGER  :: FTUVmaxprim
INTEGER  :: maxpasses
LOGICAL  :: FMM
LOGICAL  :: LINK
LOGICAL  :: DALINK
LOGICAL  :: DEBUGOVERLAP
LOGICAL  :: DEBUG4CENTER
LOGICAL  :: DEBUG4CENTER_ERI
LOGICAL  :: DEBUGDECPACKED
LOGICAL  :: DEBUGCCFRAGMENT
LOGICAL  :: DEBUGKINETIC
LOGICAL  :: DEBUGNUCREP
LOGICAL  :: DO4CENTERERI
LOGICAL  :: OVERLAP_DF_J
LOGICAL  :: TIMINGS
LOGICAL  :: nonSphericalETUV
LOGICAL  :: ETUVinsideLOOP
LOGICAL  :: HIGH_RJ000_ACCURACY
LOGICAL  :: OrderAngPrim
!*FMM PARAMETERS
Integer     :: MM_LMAX
Integer     :: MM_TLMAX
REAL(realk) :: MM_SCREEN
LOGICAL     :: NO_MMFILES
LOGICAL     :: MM_NO_ONE
LOGICAL     :: CREATED_MMFILES
LOGICAL     :: USEBUFMM
Integer     :: MMunique_ID1
!*BASIS PARAMETERS
LOGICAL  :: ATOMBASIS
LOGICAL  :: BASIS
LOGICAL  :: AUXBASIS
LOGICAL  :: NOFAMILY
LOGICAL  :: DoCartesian
LOGICAL  :: DoSpherical
LOGICAL  :: UNCONT !FORCE UNCONTRACTED BASIS
LOGICAL  :: NOSEGMENT !DISABLE SEGMENTS 

!* JOB REQUESTS
LOGICAL  :: DO3CENTEROVL
LOGICAL  :: DO2CENTERERI
INTEGER  :: CARMOM
INTEGER  :: SPHMOM
LOGICAL  :: MIXEDOVERLAP
LOGICAL  :: TEST_NDMAT_COULOMB
LOGICAL  :: TEST_NDMAT_EXCHANGE

!*CAUCHY-SCHWARZ INTEGRAL PARAMETERS
!THE ONE THRESHOLD TO RULE THEM ALL
REAL(REALK):: THRESHOLD  
!THESE THRESHOLDS TELL HOW THEY SHOULD BE SET COMPARED TO THE ONE THRESHOLD
REAL(REALK) :: CS_THRESHOLD
REAL(REALK) :: OE_THRESHOLD
REAL(REALK) :: PS_THRESHOLD
Real(realk) :: OD_THRESHOLD
REAL(REALK) :: Integralthreshold
!OTHER CAUCHY-SCHWARZ INTEGRAL PARAMETERS
LOGICAL  :: CS_SCREEN
LOGICAL  :: OE_SCREEN
LOGICAL  :: GAB_ON_FILE
LOGICAL  :: WRITE_GAB_TO_FILE
!*PRIMITIVE INTEGRAL PARAMETERS
LOGICAL  :: PS_SCREEN
LOGICAL  :: PRIM_GAB_ON_FILE
LOGICAL  :: WRITE_PRIM_GAB_TO_FILE
LOGICAL  :: PS_DEBUG
!Screen OD-batches by AO-batch extent
LOGICAL     :: OD_SCREEN 
!Fragment molecule into to distinct parts, and construct matrices block by block
LOGICAL     :: FRAGMENT
!Approximate number of atoms per fragment
Integer     :: numAtomsPerFragment

LOGICAL     :: PUREFMM
!FMM
INTEGER     :: LU_LUINTM
INTEGER     :: LU_LUINTR
!Coulomb attenuated method CAM parameters
LOGICAL     :: CAM
REAL(REALK) :: CAMalpha
REAL(REALK) :: CAMbeta
REAL(REALK) :: CAMmu

!DFT PARAMETERS
INTEGER     :: GRDONE !IF GRID HAS BEEN CREATED 
INTEGER     :: ITERATIONS !batches of dft-grid-points
REAL(REALK) :: RADINT
INTEGER     :: ANGMIN
INTEGER     :: ANGINT
INTEGER     :: HRDNES ! hardness of the partition function in the becke schemes
!
REAL(REALK) :: DFTELS !       DEFAULT: 1.0D-3
REAL(REALK) :: DFTHR0 !                1.0D-9
REAL(REALK) :: DFTHRI !                2.0D-12
REAL(REALK) :: DFTHRL !                1.0D-10
REAL(REALK) :: RHOTHR !                2.0D-15
LOGICAL     :: NOPRUN !                .FALSE.
LOGICAL     :: DFTASC !                .FALSE.
LOGICAL     :: DFTPOT !                .FALSE.
LOGICAL     :: DODISP !                .FALSE. ; empirical dispersion correction following Grimme
REAL(REALK) :: DFTIPT !                1.0D-20      
REAL(REALK) :: DFTBR1 !                1.0D-20
REAL(REALK) :: DFTBR2 !                1.0D-20
LOGICAL     :: DFTADD !                .TRUE.
LOGICAL     :: DISPDONE !              .FALSE.
!EXCHANGE FACTOR
REAL(REALK) :: exchangeFactor
INTEGER     :: TURBO
!MOLECULE INFO
INTEGER     :: nelectrons
INTEGER     :: molcharge

! TESTING FUNCTIONALITIES FOR DEC
LOGICAL     :: run_dec_gradient_test
END TYPE integralconfig

TYPE LSINTSCHEME
!PARAMETERS FROM **INTEGRALS   DECLERATION
LOGICAL  :: CFG_LSDALTON
LOGICAL  :: DOPASS
LOGICAL  :: DENSFIT
INTEGER  :: AOPRINT
INTEGER  :: INTPRINT
LOGICAL  :: JENGINE
INTEGER  :: FTUVmaxprim
INTEGER  :: maxpasses
LOGICAL  :: FMM
LOGICAL  :: LINK
LOGICAL  :: DALINK
LOGICAL  :: DEBUGOVERLAP
LOGICAL  :: DEBUG4CENTER
LOGICAL  :: DEBUG4CENTER_ERI
LOGICAL  :: DEBUGCCFRAGMENT
LOGICAL  :: DEBUGKINETIC
LOGICAL  :: DEBUGNUCREP
LOGICAL  :: DO4CENTERERI
LOGICAL  :: OVERLAP_DF_J
LOGICAL  :: TIMINGS
LOGICAL  :: nonSphericalETUV
LOGICAL  :: ETUVinsideLOOP
LOGICAL  :: HIGH_RJ000_ACCURACY
LOGICAL  :: OrderAngPrim
!*FMM PARAMETERS
Integer     :: MM_LMAX
Integer     :: MM_TLMAX
REAL(realk) :: MM_SCREEN
LOGICAL     :: NO_MMFILES
LOGICAL     :: MM_NO_ONE
LOGICAL     :: CREATED_MMFILES
LOGICAL     :: USEBUFMM
Integer     :: MMunique_ID1
!*BASIS PARAMETERS
LOGICAL  :: AUXBASIS
LOGICAL  :: NOFAMILY
LOGICAL  :: DoCartesian
LOGICAL  :: DoSpherical
LOGICAL  :: UNCONT !FORCE UNCONTRACTED BASIS
LOGICAL  :: NOSEGMENT !DISABLE SEGMENTS 

!* JOB REQUESTS
LOGICAL  :: DO3CENTEROVL
LOGICAL  :: DO2CENTERERI
INTEGER  :: CARMOM
INTEGER  :: SPHMOM
LOGICAL  :: MIXEDOVERLAP
LOGICAL  :: TEST_NDMAT_COULOMB
LOGICAL  :: TEST_NDMAT_EXCHANGE

!*CAUCHY-SCHWARZ INTEGRAL PARAMETERS
!THE ONE THRESHOLD TO RULE THEM ALL
REAL(REALK):: THRESHOLD  
!THESE THRESHOLDS TELL HOW THEY SHOULD BE SET COMPARED TO THE ONE THRESHOLD
REAL(REALK) :: CS_THRESHOLD
REAL(REALK) :: OE_THRESHOLD
REAL(REALK) :: PS_THRESHOLD
Real(realk) :: OD_THRESHOLD
REAL(REALK) :: Integralthreshold
!OTHER CAUCHY-SCHWARZ INTEGRAL PARAMETERS
LOGICAL  :: CS_SCREEN
LOGICAL  :: OE_SCREEN
LOGICAL  :: GAB_ON_FILE
LOGICAL  :: WRITE_GAB_TO_FILE
!*PRIMITIVE INTEGRAL PARAMETERS
LOGICAL  :: PS_SCREEN
LOGICAL  :: PRIM_GAB_ON_FILE
LOGICAL  :: WRITE_PRIM_GAB_TO_FILE
LOGICAL  :: PS_DEBUG
!Screen OD-batches by AO-batch extent
LOGICAL     :: OD_SCREEN 
!Fragment molecule into to distinct parts, and construct matrices block by block
LOGICAL     :: FRAGMENT
!Approximate number of atoms per fragment
Integer     :: numAtomsPerFragment

LOGICAL     :: DECPACKED !special packing for the dec-program
!FMM
LOGICAL     :: PUREFMM
INTEGER     :: LU_LUINTM
INTEGER     :: LU_LUINTR
!Coulomb attenuated method CAM parameters
LOGICAL     :: CAM
REAL(REALK) :: CAMalpha
REAL(REALK) :: CAMbeta
REAL(REALK) :: CAMmu
REAL(REALK) :: exchangeFactor !EXCHANGE FACTOR
!DFT PARAMETERS
INTEGER     :: GRDONE !IF GRID HAS BEEN CREATED 
INTEGER     :: ITERATIONS !batches of dft-grid-points
REAL(REALK) :: RADINT
INTEGER     :: ANGMIN
INTEGER     :: ANGINT
INTEGER     :: HRDNES ! hardness of the partition function in the becke schemes
REAL(REALK) :: DFTELS !       DEFAULT: 1.0D-3
REAL(REALK) :: DFTHR0 !                1.0D-9
REAL(REALK) :: DFTHRI !                2.0D-12
REAL(REALK) :: DFTHRL !                1.0D-10
REAL(REALK) :: RHOTHR !                2.0D-15
LOGICAL     :: NOPRUN !                .FALSE.
LOGICAL     :: DFTASC !                .FALSE.
LOGICAL     :: DFTPOT !                .FALSE.
LOGICAL     :: DODISP !                .FALSE. ; empirical disperision correction
REAL(REALK) :: DFTIPT !                1.0D-20      
REAL(REALK) :: DFTBR1 !                1.0D-20
REAL(REALK) :: DFTBR2 !                1.0D-20
LOGICAL     :: DFTADD !                .TRUE.
LOGICAL     :: DISPDONE !              .FALSE.
INTEGER     :: TURBO

LOGICAL :: INCREMENTAL !Use incremental scheme (density-difference KS-matrix build)
  
END TYPE LSINTSCHEME

!*****************************************
!*
!* OBJECT CONTAINING BASISSETLIBRARY
!*
!*****************************************
TYPE BASISSETLIBRARYITEM
Character(len=80)         :: BASISSETNAME(maxBasisSetInLIB)
!'GHOST' if no basissets  
integer                   :: nbasissets   
integer                   :: nCharges(maxBasisSetInLIB)
integer                   :: Charges(maxBasisSetInLIB,maxNumberOfChargesinLIB)
END TYPE BASISSETLIBRARYITEM

TYPE MAT3D
real(realk), pointer :: elements(:,:,:)
INTEGER              :: dim1,dim2,dim3
END TYPE MAT3D

TYPE BLOCK
Integer :: fragment1
Integer :: fragment2
Logical :: sameFragments
Integer :: nbast1
Integer :: nbast2
Integer :: nprimbast1
Integer :: nprimbast2
integer :: startOrb1
integer :: startOrb2
integer :: startprimOrb1
integer :: startprimOrb2
Integer :: node
Integer :: nAtoms1
Integer :: nAtoms2
END TYPE BLOCK

!******** BLOCKINFO ********
TYPE BLOCKINFO
Integer :: numBlocks
Logical :: sameAOs
TYPE(BLOCK),pointer :: blocks(:)
END TYPE BLOCKINFO


!******** FRAGMENTINFO ********
TYPE FRAGMENTINFO
Integer :: numFragments
Logical :: numberOrbialsSet
Integer :: atomsInMolecule
Integer,pointer :: fragmentIndex(:)  !Index giving the fragment of each atom
Integer,pointer :: nAtoms(:) !atoms in each fragment
! First dimension numFragments, second dimension for different basis sets:
!    1: Regular, 2: DF-Aux, 3: Huckel
Integer,pointer :: nContOrb(:,:)
Integer,pointer :: nPrimOrb(:,:)
Integer,pointer :: nStartContOrb(:,:)
Integer,pointer :: nStartPrimOrb(:,:)
END TYPE FRAGMENTINFO

TYPE FRAGMENTINFO_PT
TYPE(FRAGMENTINFO),pointer :: p
END TYPE FRAGMENTINFO_PT

!One for each of the four AO's
TYPE FRAGMENTITEM
!TYPE(MOLECULE_PT) :: MOLECULE(4)
INTEGER               :: numFragments(4) !Number of fragments to partition each molecule into
TYPE(FRAGMENTINFO_PT) :: INFO(4)
Logical               :: infoAllocated(4) = .FALSE.
Logical               :: identical(4,4)
TYPE(BLOCKINFO)       :: LHSblock
integer               :: iLHSBlock
TYPE(BLOCKINFO)       :: RHSblock
integer               :: iRHSBlock
END TYPE FRAGMENTITEM

TYPE FRAGMENT_PT
TYPE(FRAGMENTITEM),pointer :: p
END TYPE FRAGMENT_PT

!*****************************************
!*
!* OBJECT CONTAINING INFORMATION ABOUT 
!* THE INPUT TO THE DALTON PROGRAM
!*
!*****************************************
TYPE DALTONINPUT
LOGICAL                    :: DO_DFT
REAL(REALK)                :: POTNUC
integer                    :: nfock ! number of fock matrix build
TYPE(MOLECULEINFO),pointer :: MOLECULE
TYPE(integralconfig)       :: DALTON
TYPE(BASISINFO),pointer    :: BASIS
TYPE(IOITEM)               :: IO
INTEGER                    :: numFragments !Number of fragments to partition molecule into
Integer                    :: numNodes     !Number of MPI nodes
Integer                    :: node         !Integer value defining the node number
!if you add a structure to this type remember to add it to MPI_ALLOC_DALTONINPUT
END TYPE DALTONINPUT

!*****************************************
!*
!* OBJECT CONTAINING INFORMATION ABOUT 
!* THE INTEGRAL SETTINGS
!*  - used for interfacing with LSint (IntegralInterface.f90)
!*
!*****************************************
!> \brief Contrains the settings for the ls-integral routines
TYPE LSSETTING
LOGICAL                    :: DO_DFT
REAL(REALK)                :: EDISP  !empiricial dispersion correction
INTEGER                    :: nAO = 4!Number of AOs, e.g. the 4 for four-center
                                     !two-electron Coulomb-repulsion integrals (ab|cd)
TYPE(MOLECULE_PT),pointer  :: MOLECULE(:) !One for each AO
TYPE(BASIS_PT),pointer     :: BASIS(:)    !One for each AO
TYPE(MOLECULE_PT),pointer  :: FRAGMENT(:) !One for each AO
TYPE(LSINTSCHEME)          :: SCHEME !The specifications on how to run the integrals
TYPE(IOITEM)               :: IO !Keeps track of the different files that are stored on disk
INTEGER,pointer            :: Batchindex(:) !One for each AO, zero if full AObatch is requested
INTEGER,pointer            :: Batchdim(:) !dim for teach AO Batchindex, zero if full AObatch is requested
INTEGER,pointer            :: molID(:) !unique identifier for different molecules, used for screening matrices, default 0
LOGICAL,pointer            :: sameMOL(:,:)  !Specifies if the different MOLECULES are identical
LOGICAL,pointer            :: sameBAS(:,:)  !Specifies if the different BASIS are identical
LOGICAL,pointer            :: sameFRAG(:,:) !Specifies if the different FRAGMENTS are identical
LOGICAL,pointer            :: molBuild(:) !Specifies if the MOLECULEs have been built 
                                          !(i.e. not set to point to a molecule)
LOGICAL,pointer            :: basBuild(:) !Specifies if the BASISes have been built
LOGICAL,pointer            :: fragBuild(:) !Specifies if the FRAGMENTs have been built
!Density-matrix information
TYPE(matrixp),pointer      :: DmatLHS(:) !Used for passing LHS density
TYPE(matrixp),pointer      :: DmatRHS(:) !Used for passing RHS density
Real(realk),pointer        :: DfullLHS(:,:,:) !Used for passing LHS density
Real(realk),pointer        :: DfullRHS(:,:,:) !Used for passing RHS density
Integer                    :: nDmatLHS   !Number of LHS densities
Integer                    :: nDmatRHS   !Number of RHS densities
Logical                    :: LHSdmat    !Specifying whether the LHS density-matrix has been assigned
Logical                    :: RHSdmat    !Specifying whether the RHS density-matrix has been assigned
Logical                    :: LHSdfull    !Specifying whether the LHS density-matrix has been assigned
Logical                    :: RHSdfull    !Specifying whether the RHS density-matrix has been assigned
Logical                    :: LHSdalloc   !Specifying if LHSdfull alloced 
Logical                    :: RHSdalloc   !Specifying if RHSdfull alloced 
!> \param Specifying the symmetry of the LHS density-matrix for each nDmatLHS (0 not set, 1=sym, 2=anti-sym, 3=no-sym, 4=zero)
Integer,pointer            :: DsymLHS(:) 
!> \param Specifying the symmetry of the RHS density-matrix for each nDmatRHS (0 not set, 1=sym, 2=anti-sym, 3=no-sym, 4=zero)
Integer,pointer            :: DsymRHS(:) 
Integer                    :: LHSdmatAOindex1
Integer                    :: LHSdmatAOindex2
Integer                    :: RHSdmatAOindex1
Integer                    :: RHSdmatAOindex2
! Fragment and node info
INTEGER                    :: numFragments !Number of fragments to partition molecule into
TYPE(FRAGMENTITEM)         :: FRAGMENTS    !Information about the fragments
Integer                    :: numNodes     !Number of MPI nodes
Integer                    :: node         !Integer value defining the node number
!TYPE(OUTPUTSPEC)           :: OUTSPEC !The specifications on how the output should be given
TYPE(INTEGRALOUTPUT)       :: OUTPUT  !The structure containing the output 
!if you add a structure to this type remember to add it to MPI_ALLOC_DALTONINPUT
END TYPE LSSETTING

!*****************************************
!*
!* OBJECT CONTAINING INFORMATION ABOUT 
!* THE INTEGRAL SETTINGS OBJECT
!*   - used for IntegralInterface.f90
!*
!*****************************************
TYPE LSITEM
TYPE(DALTONINPUT) :: INPUT   !Input handling (of DALTON.INP and MOLECULE.INP)
TYPE(LSSETTING)   :: SETTING !Settings for integral evaluation
INTEGER           :: LUPRI = -1  !Output-file unit number
INTEGER           :: LUERR = -1  !Error-file unit number
LOGICAL           :: fopen = .false. !Determines wether the LSDALTON.OUT and LSDALTON.ERR has been opened
END TYPE LSITEM

Contains
!> \brief write lsitem to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param LS lsitem to be written
SUBROUTINE write_lsitem_to_disk(LS)
implicit none
TYPE(LSITEM)    :: LS
!
LOGICAL         :: fileexist  
INTEGER         :: lun,i

INQUIRE(file='lsitem',EXIST=fileexist)
IF(fileexist)THEN
   lun = -1
   CALL LSOPEN(LUN,'lsitem','old','UNFORMATTED')
   call LSclose(LUN,'DELETE')
ENDIF
lun = -1
CALL LSOPEN(lun,'lsitem','new','UNFORMATTED')

CALL WRITE_DALTONINPUT_TO_DISK(lun,LS%INPUT)
! IT DOES NOT MAKE SENSE TO WRITE THE LUPRI,LUERR
   call LSclose(LUN,'KEEP')

END SUBROUTINE write_lsitem_to_disk

!> \brief read lsitem from disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param LS lsitem to be written
SUBROUTINE read_lsitem_from_disk(LS)
implicit none
TYPE(LSITEM)    :: LS
!
LOGICAL         :: fileexist  
INTEGER         :: lun,i

INQUIRE(file='lsitem',EXIST=fileexist)
IF(fileexist)THEN
   lun = -1
   CALL LSOPEN(lun,'lsitem','old','UNFORMATTED')
   REWIND(lun)
   CALL READ_DALTONINPUT_FROM_DISK(lun,LS%INPUT)
   call LSclose(LUN,'KEEP')
   CALL II_init_setting(ls%setting)
   CALL II_set_default_setting(ls%setting,ls%input)
ELSE
   CALL LSQUIT('In call to read_lsitem_from_disk FILE: lsitem do not exit ',-1)
ENDIF

END SUBROUTINE read_lsitem_from_disk

!> \brief write daltoninput to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lun logical unit number to write to
!> \param input the daltoninput structure
SUBROUTINE write_daltoninput_to_disk(lun,INPUT)
implicit none
TYPE(DALTONINPUT),intent(in)  :: INPUT
INTEGER,intent(in)            :: lun

CALL WRITE_MOLECULEINFO_TO_DISK(lun,INPUT%MOLECULE)
CALL WRITE_DALTONITEM_TO_DISK(lun,INPUT%DALTON)
CALL WRITE_BASISINFO_TO_DISK(lun,INPUT%BASIS)
! IT DOES NOT MAKE SENSE TO WRITE THE IOITEM
! IT DOES NOT MAKE SENSE TO WRITE numFragments,numNodes,node

END SUBROUTINE write_daltoninput_to_disk

!> \brief read daltoninput from disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lun logical unit number to write to
!> \param input the daltoninput structure
SUBROUTINE read_daltoninput_from_disk(lun,INPUT)
implicit none
TYPE(DALTONINPUT),intent(inout)  :: INPUT
INTEGER,intent(in)            :: lun

NULLIFY(INPUT%MOLECULE)
ALLOCATE(INPUT%MOLECULE)
CALL READ_MOLECULEINFO_FROM_DISK(lun,INPUT%MOLECULE)
CALL READ_DALTONITEM_FROM_DISK(lun,INPUT%DALTON)
NULLIFY(INPUT%BASIS)
ALLOCATE(INPUT%BASIS)
CALL READ_BASISINFO_FROM_DISK(lun,INPUT%BASIS)
! IT DOES NOT MAKE SENSE TO READ THE IOITEM
CALL io_init(INPUT%IO)
! IT DOES NOT MAKE SENSE TO READ numFragments,numNodes,node
INPUT%numFragments=0
INPUT%numNodes=0
INPUT%node=0

END SUBROUTINE read_daltoninput_from_disk

!> \brief write daltonitem to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lun logical unit number to write to
!> \param dalton the daltonitem structure
SUBROUTINE write_daltonitem_to_disk(lun,DALTON)
implicit none
TYPE(integralconfig),intent(in)   :: DALTON
INTEGER,intent(in)            :: lun

WRITE(LUN) DALTON%UNRES
WRITE(LUN) DALTON%CFG_LSDALTON
WRITE(LUN) DALTON%TRILEVEL
WRITE(LUN) DALTON%DOPASS
WRITE(LUN) DALTON%DENSFIT
WRITE(LUN) DALTON%LINSCA
WRITE(LUN) DALTON%PRINTATOMCOORD
WRITE(LUN) DALTON%JENGINE
WRITE(LUN) DALTON%LOCALLINK
WRITE(LUN) DALTON%LOCALLINKmulthr
WRITE(LUN) DALTON%LOCALLINKsimmul
WRITE(LUN) DALTON%LOCALLINKoption
WRITE(LUN) DALTON%LOCALLINKincrem
WRITE(LUN) DALTON%LOCALLINKDcont
WRITE(LUN) DALTON%LOCALLINKDthr
WRITE(LUN) DALTON%FMM
WRITE(LUN) DALTON%LINK
WRITE(LUN) DALTON%DALINK
WRITE(LUN) DALTON%DEBUGOVERLAP
WRITE(LUN) DALTON%DEBUG4CENTER
WRITE(LUN) DALTON%DEBUGDECPACKED
WRITE(LUN) DALTON%DEBUG4CENTER_ERI
WRITE(LUN) DALTON%DEBUGCCFRAGMENT
WRITE(LUN) DALTON%DEBUGKINETIC
WRITE(LUN) DALTON%DEBUGNUCREP
WRITE(LUN) DALTON%DO4CENTERERI
WRITE(LUN) DALTON%OVERLAP_DF_J
WRITE(LUN) DALTON%TIMINGS
WRITE(LUN) DALTON%nonSphericalETUV
WRITE(LUN) DALTON%ETUVinsideLOOP
WRITE(LUN) DALTON%HIGH_RJ000_ACCURACY
WRITE(LUN) DALTON%OrderAngPrim
WRITE(LUN) DALTON%NO_MMFILES
WRITE(LUN) DALTON%MM_NO_ONE
WRITE(LUN) DALTON%CREATED_MMFILES
WRITE(LUN) DALTON%USEBUFMM
WRITE(LUN) DALTON%ATOMBASIS
WRITE(LUN) DALTON%BASIS
WRITE(LUN) DALTON%AUXBASIS
WRITE(LUN) DALTON%NOFAMILY
WRITE(LUN) DALTON%DoCartesian
WRITE(LUN) DALTON%DoSpherical
WRITE(LUN) DALTON%UNCONT
WRITE(LUN) DALTON%NOSEGMENT
WRITE(LUN) DALTON%DO3CENTEROVL
WRITE(LUN) DALTON%DO2CENTERERI
WRITE(LUN) DALTON%MIXEDOVERLAP
WRITE(LUN) DALTON%TEST_NDMAT_COULOMB
WRITE(LUN) DALTON%TEST_NDMAT_EXCHANGE
WRITE(LUN) DALTON%CS_SCREEN
WRITE(LUN) DALTON%OE_SCREEN
WRITE(LUN) DALTON%GAB_ON_FILE
WRITE(LUN) DALTON%WRITE_GAB_TO_FILE
WRITE(LUN) DALTON%PS_SCREEN
WRITE(LUN) DALTON%PRIM_GAB_ON_FILE
WRITE(LUN) DALTON%WRITE_PRIM_GAB_TO_FILE
WRITE(LUN) DALTON%PS_DEBUG
WRITE(LUN) DALTON%OD_SCREEN 
WRITE(LUN) DALTON%FRAGMENT
WRITE(LUN) DALTON%PUREFMM
WRITE(LUN) DALTON%CAM
WRITE(LUN) DALTON%NOPRUN
WRITE(LUN) DALTON%DFTASC
WRITE(LUN) DALTON%DFTPOT
WRITE(LUN) DALTON%DODISP
WRITE(LUN) DALTON%DFTADD

WRITE(LUN) DALTON%LINSCAPRINT
WRITE(LUN) DALTON%AOPRINT
WRITE(LUN) DALTON%MOLPRINT
WRITE(LUN) DALTON%INTPRINT
WRITE(LUN) DALTON%BASPRINT
WRITE(LUN) DALTON%FTUVmaxprim
WRITE(LUN) DALTON%maxpasses
WRITE(LUN) DALTON%MM_LMAX
WRITE(LUN) DALTON%MM_TLMAX
WRITE(LUN) DALTON%MMunique_ID1
WRITE(LUN) DALTON%CARMOM
WRITE(LUN) DALTON%SPHMOM
WRITE(LUN) DALTON%numAtomsPerFragment
WRITE(LUN) DALTON%LU_LUINTM
WRITE(LUN) DALTON%LU_LUINTR
WRITE(LUN) DALTON%GRDONE 
WRITE(LUN) DALTON%DISPDONE
WRITE(LUN) DALTON%ITERATIONS 
WRITE(LUN) DALTON%ANGMIN
WRITE(LUN) DALTON%ANGINT
WRITE(LUN) DALTON%HRDNES 
WRITE(LUN) DALTON%TURBO

WRITE(LUN) DALTON%MM_SCREEN
WRITE(LUN) DALTON%THRESHOLD  
WRITE(LUN) DALTON%CS_THRESHOLD
WRITE(LUN) DALTON%OE_THRESHOLD
WRITE(LUN) DALTON%PS_THRESHOLD
WRITE(LUN) DALTON%OD_THRESHOLD
WRITE(LUN) DALTON%Integralthreshold
WRITE(LUN) DALTON%CAMalpha
WRITE(LUN) DALTON%CAMbeta
WRITE(LUN) DALTON%CAMmu
WRITE(LUN) DALTON%RADINT
WRITE(LUN) DALTON%DFTELS 
WRITE(LUN) DALTON%DFTHR0 
WRITE(LUN) DALTON%DFTHRI 
WRITE(LUN) DALTON%DFTHRL 
WRITE(LUN) DALTON%RHOTHR 
WRITE(LUN) DALTON%DFTIPT 
WRITE(LUN) DALTON%DFTBR1 
WRITE(LUN) DALTON%DFTBR2 
WRITE(LUN) DALTON%exchangeFactor

END SUBROUTINE write_daltonitem_to_disk

!> \brief read daltonitem from disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lun logical unit number to write to
!> \param dalton the daltonitem structure
SUBROUTINE read_daltonitem_from_disk(lun,DALTON)
implicit none
TYPE(integralconfig),intent(inout)   :: DALTON
INTEGER,intent(in)            :: lun

READ(LUN) DALTON%UNRES
READ(LUN) DALTON%CFG_LSDALTON
READ(LUN) DALTON%TRILEVEL
READ(LUN) DALTON%DOPASS
READ(LUN) DALTON%DENSFIT
READ(LUN) DALTON%LINSCA
READ(LUN) DALTON%PRINTATOMCOORD
READ(LUN) DALTON%JENGINE
READ(LUN) DALTON%LOCALLINK
READ(LUN) DALTON%LOCALLINKmulthr
READ(LUN) DALTON%LOCALLINKsimmul
READ(LUN) DALTON%LOCALLINKoption
READ(LUN) DALTON%LOCALLINKincrem
READ(LUN) DALTON%LOCALLINKDcont
READ(LUN) DALTON%LOCALLINKDthr
READ(LUN) DALTON%FMM
READ(LUN) DALTON%LINK
READ(LUN) DALTON%DALINK
READ(LUN) DALTON%DEBUGOVERLAP
READ(LUN) DALTON%DEBUG4CENTER
READ(LUN) DALTON%DEBUGDECPACKED
READ(LUN) DALTON%DEBUG4CENTER_ERI
READ(LUN) DALTON%DEBUGCCFRAGMENT
READ(LUN) DALTON%DEBUGKINETIC
READ(LUN) DALTON%DEBUGNUCREP
READ(LUN) DALTON%DO4CENTERERI
READ(LUN) DALTON%OVERLAP_DF_J
READ(LUN) DALTON%TIMINGS
READ(LUN) DALTON%nonSphericalETUV
READ(LUN) DALTON%ETUVinsideLOOP
READ(LUN) DALTON%HIGH_RJ000_ACCURACY
READ(LUN) DALTON%OrderAngPrim
READ(LUN) DALTON%NO_MMFILES
READ(LUN) DALTON%MM_NO_ONE
READ(LUN) DALTON%CREATED_MMFILES
READ(LUN) DALTON%USEBUFMM
READ(LUN) DALTON%ATOMBASIS
READ(LUN) DALTON%BASIS
READ(LUN) DALTON%AUXBASIS
READ(LUN) DALTON%NOFAMILY
READ(LUN) DALTON%DoCartesian
READ(LUN) DALTON%DoSpherical
READ(LUN) DALTON%UNCONT
READ(LUN) DALTON%NOSEGMENT
READ(LUN) DALTON%DO3CENTEROVL
READ(LUN) DALTON%DO2CENTERERI
READ(LUN) DALTON%MIXEDOVERLAP
READ(LUN) DALTON%TEST_NDMAT_COULOMB
READ(LUN) DALTON%TEST_NDMAT_EXCHANGE
READ(LUN) DALTON%CS_SCREEN
READ(LUN) DALTON%OE_SCREEN
READ(LUN) DALTON%GAB_ON_FILE
READ(LUN) DALTON%WRITE_GAB_TO_FILE
READ(LUN) DALTON%PS_SCREEN
READ(LUN) DALTON%PRIM_GAB_ON_FILE
READ(LUN) DALTON%WRITE_PRIM_GAB_TO_FILE
READ(LUN) DALTON%PS_DEBUG
READ(LUN) DALTON%OD_SCREEN 
READ(LUN) DALTON%FRAGMENT
READ(LUN) DALTON%PUREFMM
READ(LUN) DALTON%CAM
READ(LUN) DALTON%NOPRUN
READ(LUN) DALTON%DFTASC
READ(LUN) DALTON%DFTPOT
READ(LUN) DALTON%DODISP
READ(LUN) DALTON%DFTADD

READ(LUN) DALTON%LINSCAPRINT
READ(LUN) DALTON%AOPRINT
READ(LUN) DALTON%MOLPRINT
READ(LUN) DALTON%INTPRINT
READ(LUN) DALTON%BASPRINT
READ(LUN) DALTON%FTUVmaxprim
READ(LUN) DALTON%maxpasses
READ(LUN) DALTON%MM_LMAX
READ(LUN) DALTON%MM_TLMAX
READ(LUN) DALTON%MMunique_ID1
READ(LUN) DALTON%CARMOM
READ(LUN) DALTON%SPHMOM
READ(LUN) DALTON%numAtomsPerFragment
READ(LUN) DALTON%LU_LUINTM
READ(LUN) DALTON%LU_LUINTR
READ(LUN) DALTON%GRDONE 
READ(LUN) DALTON%DISPDONE
READ(LUN) DALTON%ITERATIONS 
READ(LUN) DALTON%ANGMIN
READ(LUN) DALTON%ANGINT
READ(LUN) DALTON%HRDNES 
READ(LUN) DALTON%TURBO

READ(LUN) DALTON%MM_SCREEN
READ(LUN) DALTON%THRESHOLD  
READ(LUN) DALTON%CS_THRESHOLD
READ(LUN) DALTON%OE_THRESHOLD
READ(LUN) DALTON%PS_THRESHOLD
READ(LUN) DALTON%OD_THRESHOLD
READ(LUN) DALTON%Integralthreshold
READ(LUN) DALTON%CAMalpha
READ(LUN) DALTON%CAMbeta
READ(LUN) DALTON%CAMmu
READ(LUN) DALTON%RADINT
READ(LUN) DALTON%DFTELS 
READ(LUN) DALTON%DFTHR0 
READ(LUN) DALTON%DFTHRI 
READ(LUN) DALTON%DFTHRL 
READ(LUN) DALTON%RHOTHR 
READ(LUN) DALTON%DFTIPT 
READ(LUN) DALTON%DFTBR1 
READ(LUN) DALTON%DFTBR2 
READ(LUN) DALTON%exchangeFactor

END SUBROUTINE read_daltonitem_from_disk

!> \brief set the integralconfig to default values
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param dalton the integralconfig structure
SUBROUTINE integral_set_default_config(DALTON)
IMPLICIT NONE
TYPE(integralconfig)   :: DALTON

DALTON%unres = .FALSE.
DALTON%cfg_lsdalton = .FALSE.
DALTON%TRILEVEL = .TRUE.
DALTON%DOPASS = .TRUE.
DALTON%DENSFIT = .FALSE.
DALTON%LINSCA = .FALSE.
!PRINTING KEYWORDS
DALTON%LINSCAPRINT = 0
DALTON%BASPRINT = 0
DALTON%AOPRINT = 0
DALTON%MOLPRINT = 0
DALTON%INTPRINT = 0
DALTON%PRINTATOMCOORD = .FALSE.
DALTON%BASIS = .FALSE.
DALTON%ATOMBASIS = .FALSE.
DALTON%AUXBASIS = .FALSE.
DALTON%NOFAMILY = .TRUE. 
DALTON%DoCartesian = .FALSE.
DALTON%JENGINE = .TRUE.
DALTON%LOCALLINK = .FALSE.
DALTON%LOCALLINKmulthr = 1.D-5
DALTON%LOCALLINKsimmul = .FALSE.
DALTON%LOCALLINKoption = 1
DALTON%LOCALLINKincrem = .FALSE.
DALTON%LOCALLINKDcont = .FALSE.
DALTON%LOCALLINKDthr = 1.D-10
DALTON%LINK = .TRUE.
DALTON%DALINK = .FALSE.
DALTON%HIGH_RJ000_ACCURACY = .TRUE.
DALTON%OrderAngPrim = .TRUE.
!FMM PARAMETERS
DALTON%FMM = .FALSE.
DALTON%FTUVmaxprim = 64
DALTON%maxpasses = 40
DALTON%PUREFMM = .FALSE.
DALTON%MM_LMAX   = 8
DALTON%MM_TLMAX  = 20
DALTON%MM_NO_ONE = .FALSE.
DALTON%NO_MMFILES = .FALSE.
DALTON%CREATED_MMFILES = .FALSE.
DALTON%USEBUFMM = .TRUE.
DALTON%nonSphericalETUV = .FALSE.
DALTON%ETUVinsideLOOP = .FALSE.
DALTON%DEBUGOVERLAP = .FALSE.
DALTON%DEBUG4CENTER = .FALSE.
DALTON%DEBUGDECPACKED = .FALSE.
DALTON%DEBUG4CENTER_ERI = .FALSE.
DALTON%DEBUGCCFRAGMENT = .FALSE.
DALTON%DEBUGNUCREP = .FALSE.
DALTON%DO3CENTEROVL = .FALSE.
DALTON%DO2CENTERERI = .FALSE.
DALTON%DO4CENTERERI = .FALSE.
DALTON%OVERLAP_DF_J = .FALSE.
DALTON%TEST_NDMAT_COULOMB = .FALSE.
DALTON%TEST_NDMAT_EXCHANGE = .FALSE.
DALTON%TIMINGS = .FALSE.
DALTON%UNCONT = .FALSE.
DALTON%NOSEGMENT = .TRUE.
!THE ONE THRESHOLD
DALTON%THRESHOLD = 1.0D-10
DALTON%CS_THRESHOLD = 1.0D+0 ! so 1.0D-10 }   
DALTON%OE_THRESHOLD = 1.0D-1 ! so 1.0D-11 }
DALTON%PS_THRESHOLD = 1.0D-1 ! so 1.0D-11 }
DALTON%OD_THRESHOLD = 1.0D-1 ! so 1.0D-11 } in integralinput
DALTON%MM_SCREEN    = 1.0D-1 ! so 1.0D-11 }
DALTON%CS_SCREEN = .TRUE.
DALTON%OE_SCREEN = .TRUE.
DALTON%GAB_ON_FILE = .FALSE.
DALTON%WRITE_GAB_TO_FILE = .TRUE.
DALTON%PS_SCREEN = .TRUE.
DALTON%OD_SCREEN = .TRUE.
DALTON%PRIM_GAB_ON_FILE = .FALSE.
DALTON%WRITE_PRIM_GAB_TO_FILE = .TRUE.
DALTON%PS_DEBUG = .FALSE.
DALTON%DEBUGKINETIC = .FALSE.
DALTON%CARMOM = 0
DALTON%SPHMOM = 0
!DALTON%FRAGMENT = .TRUE.
DALTON%FRAGMENT = .FALSE.
!Default is to make the number of atoms so large fragmentation is not used
DALTON%numAtomsPerFragment = 10000000
DALTON%MIXEDOVERLAP = .FALSE.

!CAM PARAMETERS
DALTON%CAM = .FALSE.
DALTON%CAMalpha=0.19D0
DALTON%CAMbeta=0.46D0
DALTON%CAMmu=0.33D0

!DFT PARAMETERS
DALTON%GRDONE = 0 !FALSE = GRID NOT DONE YET  
DALTON%RADINT = 1.0D-11
DALTON%ANGMIN = 5
DALTON%ANGINT = 31
DALTON%HRDNES = 3 ! hardness of the becke partioning function
DALTON%ITERATIONS = 0 
DALTON%DFTHR0 = 1.0D-9  !not used
DALTON%DFTHRI = 2.0D-12
DALTON%DFTELS = 1.0D-3
DALTON%DFTHRL = 2.0D-10 !not used - obsolete
DALTON%RHOTHR = 2.0D-15
DALTON%NOPRUN = .FALSE. !not used at all - now it is used ! \Andreas
DALTON%DFTASC = .FALSE.
DALTON%DFTPOT = .FALSE.
DALTON%DODISP = .FALSE.
DALTON%DFTIPT = 1.0D-20      
DALTON%DFTBR1 = 1.0D-20
DALTON%DFTBR2 = 1.0D-20
DALTON%DFTADD = .TRUE.  !disable DFT D
DALTON%DISPDONE = .FALSE.
DALTON%TURBO = 0 !TURBOMOLE grid

! DEC TEST PARAMETERS
DALTON%run_dec_gradient_test=.false.
END SUBROUTINE integral_set_default_config

!> \brief attach dmat to integral input structure
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param integralinput the integral input structure
!> \param Dmat to attach
!> \param ndim1 size of dimension 1
!> \param ndim2 size of dimension 2
!> \param ndmat number of density matrices
!> \param Dmatside the side of the matrix LHS or RHS
SUBROUTINE attachDmatToInput(Input,Dmat,ndim1,ndim2,ndmat,DmatSide)
implicit none
Type(IntegralInput) :: Input
Integer             :: ndim1,ndim2,ndmat
Real(realk),target  :: Dmat(ndim1,ndim2,ndmat)
Character(3)        :: DmatSide
IF (DmatSide.EQ.'LHS') THEN
  Input%LHS_DMAT     = .TRUE.
  Input%DMAT_LHS     => Dmat
  Input%NDMAT_LHS    = ndmat
  Input%NDIM_LHS(1)  = ndim1
  Input%NDIM_LHS(2)  = ndim2
ELSE IF (DmatSide.EQ.'RHS') THEN
  Input%RHS_DMAT     = .TRUE.
  Input%DMAT_RHS     => Dmat
  Input%NDMAT_RHS    = ndmat
  Input%NDIM_RHS(1)  = ndim1
  Input%NDIM_RHS(2)  = ndim2
ELSE
  CALL LSQUIT('Programming error: attachDmatToInput called with wrong argument',-1)
ENDIF
END SUBROUTINE attachDmatToInput

!*****************************************
!*
!*  AOBATCH INITIATION ROUTINES
!*
!*****************************************

!> \brief mpi allocate the daltoninput structure
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param dalton the dalton input structure
SUBROUTINE LSMPI_ALLOC_DALTONINPUT(DALTON)
IMPLICIT NONE
TYPE(DALTONINPUT) :: DALTON
! THE MOLECULE
NULLIFY(DALTON%MOLECULE%ATOM)
ALLOCATE(DALTON%MOLECULE%ATOM(DALTON%MOLECULE%nAtoms))
!THE BASISSET
CALL LSMPI_ALLOC_BASISSETINFO(DALTON%BASIS%REGULAR)
CALL LSMPI_ALLOC_BASISSETINFO(DALTON%BASIS%HUCKEL)
CALL LSMPI_ALLOC_BASISSETINFO(DALTON%BASIS%AUXILIARY)

END SUBROUTINE LSMPI_ALLOC_DALTONINPUT

!> \brief print lsitem 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param LS the lsitem
!> \param lupri the logical unit number to write to
SUBROUTINE PRINT_LSITEM(LS,LUPRI)
IMPLICIT NONE
TYPE(LSITEM) :: LS
INTEGER      :: LUPRI

CALL PRINT_DALTONINPUT(LS%INPUT,LUPRI)
CALL PRINT_LSSETTING(LS%SETTING,LUPRI)

END SUBROUTINE PRINT_LSITEM

!> \brief print the dalton input structure
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param dalton the dalton input
!> \param lupri the logical unit number to write to
SUBROUTINE PRINT_DALTONINPUT(DALTON,LUPRI)
IMPLICIT NONE
TYPE(DALTONINPUT) :: DALTON
INTEGER           :: LUPRI

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'THE DALTON INPUT STUCTUR'
WRITE(LUPRI,*) '                     '
CALL PRINT_MOLECULEINFO(LUPRI,DALTON%MOLECULE,DALTON%BASIS)
CALL PRINT_MOLECULE_AND_BASIS(LUPRI,DALTON%MOLECULE,DALTON%BASIS%REGULAR)
CALL PRINT_BASISSETINFO(LUPRI,DALTON%BASIS%REGULAR)
CALL PRINT_DALTONITEM(LUPRI,DALTON%DALTON)
CALL PRINT_IOITEM(DALTON%IO,LUPRI)

IF(DALTON%DALTON%AUXBASIS)THEN
   write(lupri,*)'THE DALTON%BASIS%AUXILIARY'
   CALL PRINT_MOLECULE_AND_BASIS(LUPRI,DALTON%MOLECULE,DALTON%BASIS%AUXILIARY)
   CALL PRINT_BASISSETINFO(LUPRI,DALTON%BASIS%AUXILIARY)
ENDIF

END SUBROUTINE PRINT_DALTONINPUT

!> \brief print the lssetting 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param set the lssetting structure
!> \param lupri the logical unit number to write to
SUBROUTINE PRINT_LSSETTING(SET,LUPRI)
IMPLICIT NONE
TYPE(LSSETTING) :: SET
INTEGER         :: LUPRI
!
INTEGER :: I,J,dim1,dim2,nrow,ncol,imat

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'THE SETTING STUCTUR'
WRITE(LUPRI,*) '                     '
WRITE(LUPRI,*)'DO_DFT',SET%DO_DFT
DO I=1,SET%nAO
   WRITE(LUPRI,*)'AO NUMBER',I
   WRITE(LUPRI,*)'MOLECULE AND BASIS'
   CALL PRINT_MOLECULEINFO(LUPRI,SET%MOLECULE(I)%p,SET%BASIS(I)%p)
   CALL PRINT_MOLECULE_AND_BASIS(LUPRI,SET%MOLECULE(I)%p,SET%BASIS(I)%p%REGULAR)   
   CALL PRINT_BASISSETINFO(LUPRI,SET%BASIS(I)%p%REGULAR)
   WRITE(LUPRI,*)'FRAGMENT AND BASIS'
   CALL PRINT_MOLECULEINFO(LUPRI,SET%FRAGMENT(I)%p,SET%BASIS(I)%p)
   CALL PRINT_MOLECULE_AND_BASIS(LUPRI,SET%FRAGMENT(I)%p,SET%BASIS(I)%p%REGULAR)   
   WRITE(LUPRI,*)'Batchindex(',I,')=',SET%Batchindex(I)
   WRITE(LUPRI,*)'Batchdim  (',I,')=',SET%Batchdim(I)
   WRITE(LUPRI,*)'molID     (',I,')=',SET%molID(I)
   WRITE(LUPRI,*)'molBuild  (',I,')=',SET%molBUILD(I)
   WRITE(LUPRI,*)'basBuild  (',I,')=',SET%basBUILD(I)
   WRITE(LUPRI,*)'fragBuild  (',I,')=',SET%fragBUILD(I)
ENDDO
IF(SET%RHSdmat)THEN
   DO imat = 1,SET%nDmatRHS
      nrow = SET%DmatRHS(Imat)%p%nrow
      ncol = SET%DmatRHS(Imat)%p%ncol
      WRITE(LUPRI,*)'THE RHS DENSITY TYPE(MATRIX) NR=',Imat,'Dim=',nrow,ncol,'Sym=',SET%DsymRHS(Imat)
      call mat_print(SET%DmatRHS(Imat)%p,1,nrow,1,ncol,lupri)
   enddo
ENDIF
IF(SET%RHSdfull)THEN
   dim1 = SIZE(SET%DfullRHS, 1)  
   dim2 = SIZE(SET%DfullRHS, 2)  
   DO imat = 1,SET%nDmatRHS
      WRITE(LUPRI,*)'THE RHS DENSITY REAL(REALK) NR=',Imat,'Dim=',dim1,dim2,'Sym=',SET%DsymRHS(Imat)
      call output(SET%DfullRHS(:,:,imat),1,dim1,1,dim2,dim1,dim2,1,lupri)
   ENDDO
ENDIF
IF(SET%LHSdmat)THEN
   DO imat = 1,SET%nDmatLHS
      nrow = SET%DmatLHS(Imat)%p%nrow
      ncol = SET%DmatLHS(Imat)%p%ncol
      WRITE(LUPRI,*)'THE RHS DENSITY TYPE(MATRIX) NR=',Imat,'Dim=',nrow,ncol,'Sym=',SET%DsymLHS(Imat)
      call mat_print(SET%DmatLHS(Imat)%p,1,nrow,1,ncol,lupri)
   enddo
ENDIF
IF(SET%LHSdfull)THEN
   dim1 = SIZE(SET%DfullLHS, 1)  
   dim2 = SIZE(SET%DfullLHS, 2)  
   DO imat = 1,SET%nDmatLHS
      WRITE(LUPRI,*)'THE RHS DENSITY REAL(REALK) NR=',Imat,'Dim=',dim1,dim2,'Sym=',SET%DsymLHS(Imat)
      call output(SET%DfullLHS(:,:,imat),1,dim1,1,dim2,dim1,dim2,1,lupri)
   ENDDO
ENDIF

WRITE(LUPRI,*)'THE IO'
CALL PRINT_IOITEM(SET%IO,LUPRI)
WRITE(LUPRI,*)'THE LSINTSCHEM'
call printLsintScheme(SET%SCHEME,LUPRI)

DO I=1,SET%nAO
   WRITE(LUPRI,*)'SameMol ',(SET%SameMol(I,J),J=1,SET%nAO)
ENDDO
DO I=1,SET%nAO
   WRITE(LUPRI,*)'SameBas ',(SET%Samebas(I,J),J=1,SET%nAO)
ENDDO
DO I=1,SET%nAO
   WRITE(LUPRI,*)'SameFrag',(SET%SameFrag(I,J),J=1,SET%nAO)
ENDDO

WRITE(LUPRI,*)'LHSdmatAOindex1',SET%LHSdmatAOindex1
WRITE(LUPRI,*)'LHSdmatAOindex2',SET%LHSdmatAOindex2
WRITE(LUPRI,*)'RHSdmatAOindex1',SET%RHSdmatAOindex1
WRITE(LUPRI,*)'RHSdmatAOindex2',SET%RHSdmatAOindex2

WRITE(LUPRI,*)'numFragments',SET%numFragments
WRITE(LUPRI,*)'numNodes    ',SET%numNodes
WRITE(LUPRI,*)'node        ',SET%node

END SUBROUTINE PRINT_LSSETTING

!> \brief print the moleculeinfo structure
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number to write to
!> \param molecule the molecule info structure
!> \param basis the basis info structure
SUBROUTINE PRINT_MOLECULEINFO(LUPRI,MOLECULE,BASIS)
implicit none
TYPE(MOLECULEINFO) :: MOLECULE
TYPE(BASISINFO)    :: BASIS
INTEGER            :: I,J
INTEGER            :: LUPRI,ICHARGE,ITYPE1,ITYPE2

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'THE MOLECULE'

WRITE(LUPRI,*) '--------------------------------------------------------------------'
WRITE(LUPRI,'(A38,2X,F8.4)')'Molecular Charge                 :',MOLECULE%charge
WRITE(LUPRI,'(2X,A38,2X,I7)')'Regular basisfunctions             :',MOLECULE%nbastREG
WRITE(LUPRI,'(2X,A38,2X,I7)')'Auxiliary basisfunctions           :',MOLECULE%nbastAUX
WRITE(LUPRI,'(2X,A38,2X,I7)')'Huckel basisfunctions              :',MOLECULE%nbastHUC
WRITE(LUPRI,'(2X,A38,2X,I7)')'Primitive Regular basisfunctions   :',MOLECULE%nprimbastREG
WRITE(LUPRI,'(2X,A38,2X,I7)')'Primitive Auxiliary basisfunctions :',MOLECULE%nprimbastAUX
WRITE(LUPRI,'(2X,A38,2X,I7)')'Primitive Huckel basisfunctions    :',MOLECULE%nprimbastHUC
WRITE(LUPRI,*) '--------------------------------------------------------------------'
WRITE(LUPRI,*) '                     '

WRITE(LUPRI,*) '                     '
IF(MOLECULE%ATOM(1)%nbasis == 2) THEN
   WRITE(LUPRI,*) '--------------------------------------------------------------------'
   WRITE(LUPRI,'(2X,A4,2X,A6,2X,A12,2X,A20,2X,A7,2X,A8,2X,A8)')'atom',&
        &'charge','Atomicbasis ','Auxiliarybasisset',' GHOST ','nPrimREG','nContREG'
   WRITE(LUPRI,*) '--------------------------------------------------------------------'
ELSE
   WRITE(LUPRI,*) '--------------------------------------------------------------------'
   WRITE(LUPRI,'(2X,A4,2X,A6,2X,A12,2X,A7,2X,A8,2X,A8)')'atom',&
        &'charge','Atomicbasis ',' GHOST ','nPrimREG','nContREG'
   WRITE(LUPRI,*) '--------------------------------------------------------------------'
ENDIF

IF(MOLECULE%nAtoms .GT. 30)THEN
   IF(BASIS%REGULAR%Labelindex .EQ. 0)THEN
      DO I=1,30
         IF(MOLECULE%ATOM(I)%nbasis == 2) THEN
            ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
            ITYPE1 = BASIS%REGULAR%CHARGEINDEX(ICHARGE)
            ITYPE2 = BASIS%AUXILIARY%CHARGEINDEX(ICHARGE)
            WRITE(LUPRI,'(2X,I4,2X,F6.3,2X,A12,2X,A20,4X,L1,10X,I5,7X,I5)') I,MOLECULE%ATOM(I)%Charge,&
                 &BASIS%REGULAR%ATOMTYPE(ITYPE1)%NAME,&
                 &BASIS%AUXILIARY%ATOMTYPE(ITYPE2)%NAME,&
                 &MOLECULE%ATOM(I)%GHOST,&
                 &MOLECULE%ATOM(I)%nPrimOrbREG,MOLECULE%ATOM(I)%nContOrbREG
         ELSE
            ICHARGE = MOLECULE%ATOM(I)%CHARGE
            ITYPE1 = BASIS%REGULAR%CHARGEINDEX(ICHARGE)
            WRITE(LUPRI,'(2X,I4,2X,F6.3,2X,A12,4X,L1,10X,I5,7X,I5)') I,MOLECULE%ATOM(I)%Charge,&
                 &BASIS%REGULAR%ATOMTYPE(ITYPE1)%NAME,&
                 &MOLECULE%ATOM(I)%GHOST,&
                 &MOLECULE%ATOM(I)%nPrimOrbREG,MOLECULE%ATOM(I)%nContOrbREG
         ENDIF
      ENDDO
      WRITE(LUPRI,'(2X,A)')'Since you have more than 30 atoms only the first 30'
      WRITE(LUPRI,'(2X,A)')'are printed in order to limit output'
   ELSE
      DO I=1,30
         IF(MOLECULE%ATOM(I)%nbasis == 2) THEN
            ITYPE1 = MOLECULE%ATOM(I)%IDtype(1)
            ITYPE2 = MOLECULE%ATOM(I)%IDtype(2)
            WRITE(LUPRI,'(2X,I4,2X,F6.3,2X,A12,2X,A20,4X,L1,10X,I5,7X,I5)') I,MOLECULE%ATOM(I)%Charge,&
                 &BASIS%REGULAR%ATOMTYPE(ITYPE1)%NAME,&
                 &BASIS%AUXILIARY%ATOMTYPE(ITYPE2)%NAME,&
                 &MOLECULE%ATOM(I)%GHOST,&
                 &MOLECULE%ATOM(I)%nPrimOrbREG,MOLECULE%ATOM(I)%nContOrbREG
         ELSE
            ITYPE1 = MOLECULE%ATOM(I)%IDtype(1)
            WRITE(LUPRI,'(2X,I4,2X,F6.3,2X,A12,4X,L1,10X,I5,7X,I5)') I,MOLECULE%ATOM(I)%Charge,&
                 &BASIS%REGULAR%ATOMTYPE(ITYPE1)%NAME,&
                 &MOLECULE%ATOM(I)%GHOST,&
                 &MOLECULE%ATOM(I)%nPrimOrbREG,MOLECULE%ATOM(I)%nContOrbREG
         ENDIF
      ENDDO
      WRITE(LUPRI,'(2X,A)')'Since you have more than 30 atoms only the first 30'
      WRITE(LUPRI,'(2X,A)')'are printed in order to limit output'
   ENDIF
ELSE
   IF(BASIS%REGULAR%Labelindex .EQ. 0)THEN
      DO I=1,MOLECULE%nAtoms
         IF(MOLECULE%ATOM(I)%nbasis == 2) THEN
            ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
            ITYPE1 = BASIS%REGULAR%CHARGEINDEX(ICHARGE)
            ITYPE2 = BASIS%AUXILIARY%CHARGEINDEX(ICHARGE)
            WRITE(LUPRI,'(2X,I4,2X,F6.3,2X,A12,2X,A20,4X,L1,10X,I5,7X,I5)') I,MOLECULE%ATOM(I)%Charge,&
                 &BASIS%REGULAR%ATOMTYPE(itype1)%NAME,&
                 &BASIS%AUXILIARY%ATOMTYPE(itype2)%NAME,&
                 &MOLECULE%ATOM(I)%GHOST,&
                 &MOLECULE%ATOM(I)%nPrimOrbREG,MOLECULE%ATOM(I)%nContOrbREG
         ELSE
            ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
            ITYPE1 = BASIS%REGULAR%CHARGEINDEX(ICHARGE)
            WRITE(LUPRI,'(2X,I4,2X,F6.3,2X,A12,4X,L1,10X,I5,7X,I5)') I,MOLECULE%ATOM(I)%Charge,&
                 &BASIS%REGULAR%ATOMTYPE(itype1)%NAME,&
                 &MOLECULE%ATOM(I)%GHOST,&
                 &MOLECULE%ATOM(I)%nPrimOrbREG,MOLECULE%ATOM(I)%nContOrbREG
         ENDIF
      ENDDO
   ELSE
      DO I=1,MOLECULE%nAtoms
         IF(MOLECULE%ATOM(I)%nbasis == 2) THEN
            ITYPE1 = MOLECULE%ATOM(I)%IDtype(1)
            ITYPE2 = MOLECULE%ATOM(I)%IDtype(2)
            WRITE(LUPRI,'(2X,I4,2X,F6.3,2X,A12,2X,A20,4X,L1,10X,I5,7X,I5)') I,MOLECULE%ATOM(I)%Charge,&
                 &BASIS%REGULAR%ATOMTYPE(itype1)%NAME,&
                 &BASIS%AUXILIARY%ATOMTYPE(itype2)%NAME,&
                 &MOLECULE%ATOM(I)%GHOST,&
                 &MOLECULE%ATOM(I)%nPrimOrbREG,MOLECULE%ATOM(I)%nContOrbREG
         ELSE
            ITYPE1 = MOLECULE%ATOM(I)%IDtype(1)
            WRITE(LUPRI,'(2X,I4,2X,F6.3,2X,A12,4X,L1,10X,I5,7X,I5)') I,MOLECULE%ATOM(I)%Charge,&
                 &BASIS%REGULAR%ATOMTYPE(itype1)%NAME,&
                 &MOLECULE%ATOM(I)%GHOST,&
                 &MOLECULE%ATOM(I)%nPrimOrbREG,MOLECULE%ATOM(I)%nContOrbREG
         ENDIF
      ENDDO
   ENDIF
ENDIF

WRITE(LUPRI,*) '                     '

    WRITE(LUPRI,'(2X,A4,2X,A4,2X,A7,2X,A16,2X,A16,2X,A16)')'ATOM','NAME',&
    &'ISOTOPE','      X     ','      Y     ','      Z     '
IF(MOLECULE%nAtoms .GT. 30)THEN
   DO I=1,30
      WRITE(LUPRI,'(2X,I4,2X,A4,2X,I7,2X,F16.8,2X,F16.8,2X,F16.8)') I,&
           & MOLECULE%ATOM(I)%Name,&
           & MOLECULE%ATOM(I)%Isotope,&
           & MOLECULE%ATOM(I)%CENTER(1),&
           & MOLECULE%ATOM(I)%CENTER(2),&
           & MOLECULE%ATOM(I)%CENTER(3)
   ENDDO
   WRITE(LUPRI,'(2X,A)')'Since you have more than 30 atoms only the first 30'
   WRITE(LUPRI,'(2X,A)')'are printed in order to limit output'
ELSE
   DO I=1,MOLECULE%nAtoms
      WRITE(LUPRI,'(2X,I4,2X,A4,2X,I7,2X,F16.8,2X,F16.8,2X,F16.8)') I,&
           & MOLECULE%ATOM(I)%Name,&
           & MOLECULE%ATOM(I)%Isotope,&
           & MOLECULE%ATOM(I)%CENTER(1),&
           & MOLECULE%ATOM(I)%CENTER(2),&
           & MOLECULE%ATOM(I)%CENTER(3)
   ENDDO
ENDIF

WRITE(LUPRI,*) '                     '

END SUBROUTINE PRINT_MOLECULEINFO

!> \brief PRINT BASISSETLIBRARY
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number to write to
!> \param basissetlibrary
SUBROUTINE PRINT_BASISSETLIBRARY(LUPRI,BASISSETLIBRARY)
implicit none
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
INTEGER            :: I,J
INTEGER            :: LUPRI
CHARACTER(len=9)   :: STRINGFORMAT
WRITE(LUPRI,*) '  '
WRITE(LUPRI,'(A)')'BASISSETLIBRARY'
WRITE(LUPRI,*)'Number of Basisset',BASISSETLIBRARY%nbasissets
DO I=1,BASISSETLIBRARY%nbasissets
   WRITE(LUPRI,'(A10,2X,A50)')'BASISSET:',BASISSETLIBRARY%BASISSETNAME(I)
   IF(BASISSETLIBRARY%nCharges(I) < 10)THEN
      WRITE(StringFormat,'(A5,I1,A3)') '(A10,',BASISSETLIBRARY%nCharges(I),'I4)'
   ELSE
      WRITE(StringFormat,'(A5,I2,A3)') '(A10,',BASISSETLIBRARY%nCharges(I),'I4)'
   ENDIF
   WRITE(LUPRI,StringFormat)'CHARGES:',(BASISSETLIBRARY%Charges(I,J)&
        &,J=1,BASISSETLIBRARY%nCharges(I))
ENDDO
WRITE(LUPRI,*) '                     '

END SUBROUTINE PRINT_BASISSETLIBRARY

!> \brief PRINT the integral config structure
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number to write to
!> \param dalton the integral config structure to be printet
SUBROUTINE PRINT_DALTONITEM(LUPRI,DALTON)
implicit none
TYPE(integralconfig)  :: DALTON
INTEGER           :: LUPRI

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'THE DALTONITEM'
WRITE(LUPRI,*)' '
WRITE(LUPRI,'(2X,A35,7X,L1)')'DENSFIT',DALTON%DENSFIT
WRITE(LUPRI,'(2X,A35,7X,L1)')'LINSCA',DALTON%LINSCA
WRITE(LUPRI,'(2X,A35,I8)')'LINSCAPRINT',DALTON%LINSCAPRINT
WRITE(LUPRI,'(2X,A35,I8)')'AOPRINT',DALTON%AOPRINT
WRITE(LUPRI,'(2X,A35,I8)')'MOLPRINT',DALTON%MOLPRINT
WRITE(LUPRI,'(2X,A35,7X,L1)')'JENGINE',DALTON%JENGINE
WRITE(LUPRI,'(2X,A35,7X,L1)')'LOCALLINK',DALTON%LOCALLINK
WRITE(LUPRI,'(2X,A35,7X,F16.8)')'LOCALLINKmulthr',DALTON%LOCALLINKmulthr
WRITE(LUPRI,'(2X,A35,7X,L1)')'LOCALLINKsimmul',DALTON%LOCALLINKsimmul
WRITE(LUPRI,'(2X,A35,7X,I8)')'LOCALLINKoption',DALTON%LOCALLINKoption
WRITE(LUPRI,'(2X,A35,7X,L1)')'LOCALLINKincrem',DALTON%LOCALLINKincrem
WRITE(LUPRI,'(2X,A35,7X,L1)')'LOCALLINKcont',DALTON%LOCALLINKDcont
WRITE(LUPRI,'(2X,A35,7X,F16.8)')'LOCALLINKDthr',DALTON%LOCALLINKDthr
WRITE(LUPRI,'(2X,A35,7X,L1)')'FMM',DALTON%FMM
WRITE(LUPRI,'(2X,A35,7X,L1)')'LINK',DALTON%LINK
WRITE(LUPRI,'(2X,A35,7X,L1)')'DALINK',DALTON%DALINK
WRITE(LUPRI,'(2X,A35,I8)')'INTPRINT',DALTON%INTPRINT
WRITE(LUPRI,'(2X,A35,7X,L1)')'DEBUGOVERLAP',DALTON%DEBUGOVERLAP
WRITE(LUPRI,'(2X,A35,7X,L1)')'DEBUG4CENTER',DALTON%DEBUG4CENTER
WRITE(LUPRI,'(2X,A35,7X,L1)')'DEBUGDECPACKED',DALTON%DEBUGDECPACKED
WRITE(LUPRI,'(2X,A35,7X,L1)')'DEBUG4CENTER_ERI',DALTON%DEBUG4CENTER_ERI
WRITE(LUPRI,'(2X,A35,7X,L1)')'DEBUGCCFRAGMENT',DALTON%DEBUGCCFRAGMENT
WRITE(LUPRI,'(2X,A35,7X,L1)')'DEBUGKINETIC',DALTON%DEBUGKINETIC
WRITE(LUPRI,'(2X,A35,7X,L1)')'DEBUGNUCREP',DALTON%DEBUGNUCREP
WRITE(LUPRI,'(2X,A35,7X,L1)')'DO4CENTERERI',DALTON%DO4CENTERERI
WRITE(LUPRI,'(2X,A35,7X,L1)')'OVERLAP_DF_J',DALTON%OVERLAP_DF_J
WRITE(LUPRI,'(2X,A35,7X,L1)')'TIMINGS',DALTON%TIMINGS
WRITE(LUPRI,'(2X,A35,7X,L1)')'nonSphericalETUV',DALTON%nonSphericalETUV
WRITE(LUPRI,'(2X,A35,7X,L1)')'ETUVinsideLOOP',DALTON%ETUVinsideLOOP

WRITE(LUPRI,'(2X,A35,I8)')'BASPRINT',DALTON%BASPRINT
WRITE(LUPRI,'(2X,A35,7X,L1)')'ATOMBASIS',DALTON%ATOMBASIS
WRITE(LUPRI,'(2X,A35,7X,L1)')'BASIS',DALTON%BASIS
WRITE(LUPRI,'(2X,A35,7X,L1)')'AUXBASIS',DALTON%AUXBASIS
WRITE(LUPRI,'(2X,A35,7X,L1)')'NOFAMILY',DALTON%NOFAMILY
WRITE(LUPRI,'(2X,A35,7X,L1)')'DoCartesian',DALTON%DoCartesian
WRITE(LUPRI,'(2X,A35,7X,L1)')'DoSpherical',DALTON%DoSpherical
WRITE(LUPRI,'(2X,A35,7X,L1)')'UNCONT',DALTON%UNCONT
WRITE(LUPRI,'(2X,A35,7X,L1)')'NOSEGMENT',DALTON%NOSEGMENT

WRITE(LUPRI,'(2X,A35,7X,L1)')'DO3CENTEROVL',DALTON%DO3CENTEROVL
WRITE(LUPRI,'(2X,A35,7X,L1)')'DO2CENTERERI',DALTON%DO2CENTERERI
WRITE(LUPRI,'(2X,A35,I8)')'CARMOM',DALTON%CARMOM
WRITE(LUPRI,'(2X,A35,7X,L1)')'MIXEDOVERLAP',DALTON%MIXEDOVERLAP

!*CAUCHY-SCHWARZ INTEGRAL PARAMETERS
WRITE(LUPRI,'(2X,A35,F16.8)')'CS_THRESHOLD',DALTON%CS_THRESHOLD
WRITE(LUPRI,'(2X,A35,7X,L1)')'CS_SCREEN',DALTON%CS_SCREEN
WRITE(LUPRI,'(2X,A35,7X,L1)')'GAB_ON_FILE',DALTON%GAB_ON_FILE
WRITE(LUPRI,'(2X,A35,7X,L1)')'WRITE_GAB_TO_FILE',DALTON%WRITE_GAB_TO_FILE
!*PRIMITIVE INTEGRAL PARAMETERS
WRITE(LUPRI,'(2X,A35,F16.8)')'PS_THRESHOLD',DALTON%PS_THRESHOLD
WRITE(LUPRI,'(2X,A35,7X,L1)')'PS_SCREEN',DALTON%PS_SCREEN
WRITE(LUPRI,'(2X,A35,7X,L1)')'PRIM_GAB_ON_FILE',DALTON%PRIM_GAB_ON_FILE
WRITE(LUPRI,'(2X,A35,7X,L1)')'WRITE_PRIM_GAB_TO_FILE',DALTON%WRITE_PRIM_GAB_TO_FILE
WRITE(LUPRI,'(2X,A35,7X,L1)')'PS_DEBUG',DALTON%PS_DEBUG
WRITE(LUPRI,'(2X,A35,7X,L1)')'FRAGMENT',DALTON%FRAGMENT
WRITE(LUPRI,'(2X,A35,7X,L1)')'numAtomsPerFragment',DALTON%numAtomsPerFragment

!Coulomb attenuated method CAM parameters
WRITE(LUPRI,'(2X,A35,7X,L1)')'CAM',DALTON%CAM
WRITE(LUPRI,'(2X,A35,F16.8)') 'CAMalpha',DALTON%CAMalpha
WRITE(LUPRI,'(2X,A35,F16.8)') 'CAMbeta',DALTON%CAMbeta
WRITE(LUPRI,'(2X,A35,F16.8)') 'CAMmu',DALTON%CAMmu

!DFT PARAMETERS
WRITE(LUPRI,'(2X,A35,I8)')'GRDONE',DALTON%GRDONE 
WRITE(LUPRI,'(2X,A35,I8)')'ITERATIONS ',DALTON%ITERATIONS 
WRITE(LUPRI,'(2X,A35,F16.8)') 'RADINT',DALTON%RADINT
WRITE(LUPRI,'(2X,A35,I8)')'ANGMIN',DALTON%ANGMIN
WRITE(LUPRI,'(2X,A35,I8)')'ANGINT',DALTON%ANGINT
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTELS',DALTON%DFTELS
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTHR0',DALTON%DFTHR0
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTHRI',DALTON%DFTHRI
WRITE(LUPRI,'(2X,A35,L)') 'DISPER', DALTON%DODISP
WRITE(LUPRI,'(2X,A35,L)') 'DISPDONE',DALTON%DISPDONE
!EXCHANGE FACTOR
WRITE(LUPRI,'(2X,A35,F16.8)') 'exchangeFactor',DALTON%exchangeFactor

END SUBROUTINE PRINT_DALTONITEM

!> \brief print the IO item structure
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the ioitem to be printet
!> \param lupri the logical unit number to write to
SUBROUTINE PRINT_IOITEM(IO,LUPRI)
implicit none
TYPE(IOITEM)  :: IO
INTEGER       :: LUPRI
!
INTEGER :: I

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'THE IOITEM'
WRITE(LUPRI,*)' '
WRITE(LUPRI,'(2X,A35,7X,I4)')'numfiles',IO%numFiles
DO I=1,IO%numFiles
   WRITE(LUPRI,'(2X,A15,A30)')'FILENAME',IO%filename(I)
   WRITE(LUPRI,'(2X,A15,I10)')'IUNIT   ',IO%IUNIT(I)
   WRITE(LUPRI,'(2X,A15,L1)')'isopen  ',IO%isopen(I)
ENDDO

END SUBROUTINE PRINT_IOITEM

!> \brief print the molceule and basis info
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number to write to
!> \param molecule the molecule structure
!> \param basinfo the basis info structure
SUBROUTINE PRINT_MOLECULE_AND_BASIS(LUPRI,MOLECULE,basInfo)
IMPLICIT NONE
TYPE(BASISSETINFO),intent(in) :: basInfo
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
INTEGER,intent(in)            :: LUPRI
!
INTEGER             :: I,TOTCHARGE,TOTprim,TOTcont,icharge,R,set,type
CHARACTER(len=35)   :: CC 
LOGICAl             :: printbas,printed_dots
INTEGER             :: oldset,oldtype,K,J,atom_of_type,L
WRITE(LUPRI,*) '                     '

   WRITE(LUPRI,'(A)')'Atoms and basis sets'
   WRITE (LUPRI,'(A,I7)') '  Total number of atoms        :',MOLECULE%natoms
   TOTCHARGE=0
   TOTprim=0
   TOTcont=0

R=basInfo%labelindex
WRITE(LUPRI,'(2X,A3,2X,A9,A10,I4)')'THE',basInfo%label,' is on R =',R
   WRITE (LUPRI,'(A)')'--------------------------------------------------------&
        &-------------'
   WRITE (LUPRI,'(1X,A4,1X,A5,1X,A6,2X,A8,12X,A4,1X,A4,1X,A5)')'atom','label','charge'&
        &,'basisset','prim','cont','basis'
   WRITE (LUPRI,'(A)')'--------------------------------------------------------&
                                                                &-------------'
oldset=0
oldtype=0
atom_of_type = 0
printed_dots=.FALSE.
PRINTBAS = .TRUE.
DO I=1,MOLECULE%nAtoms
   IF(basInfo%labelindex .EQ.0)THEN
      ICHARGE = INT(MOLECULE%ATOM(I)%charge) 
      type= basInfo%Chargeindex(ICHARGE)
   ELSE
      type=MOLECULE%ATOM(I)%IDtype(R)
   ENDIF
   IF(oldtype .EQ. type)THEN !OLDTYPE
      atom_of_type = atom_of_type+1
      IF(atom_of_type.GT.10)THEN
         PRINTBAS = .FALSE.
         IF(.NOT.printed_dots)THEN
            WRITE (LUPRI,'(A7,1X,A4,1X,A7,1X,A20,1X,A7,2X,A7,1X,A25)') &
                 & '      :',':   ' ,'    :  ','    :               ',&
                 & '      :','      :','       :                 '
            printed_dots=.TRUE.
         ENDIF
      ENDIF
   ELSE
      atom_of_type = 0 
      printed_dots=.FALSE.
   ENDIF

   TOTCHARGE=TOTCHARGE+basInfo%ATOMTYPE(type)%Charge
   TOTprim=TOTprim+basInfo%ATOMTYPE(type)%Totnprim
   TOTcont=TOTcont+basInfo%ATOMTYPE(type)%Totnorb
   IF(PRINTBAS)THEN
      CALL BASTYP(LUPRI,MOLECULE,basInfo,I,type,CC)         
      WRITE (LUPRI,'(I7,1X,A4,1X,F7.3,1X,A20,1X,I7,2X,I7,1X,A25)') &
           & I,MOLECULE%ATOM(I)%NAME,&
           & MOLECULE%ATOM(I)%Charge,&
           & basInfo%ATOMTYPE(type)%Name,&
           & basInfo%ATOMTYPE(type)%Totnprim,&
           & basInfo%ATOMTYPE(type)%Totnorb, CC
   ENDIF
   oldtype = type
   PRINTBAS = .TRUE.
ENDDO

WRITE (LUPRI,'(A)')'--------------------------------------------------------&
     &-------------'
WRITE (LUPRI,'(A9,I7,25X,I7,1X,I7)')'total         ',TOTCHARGE,TOTprim,TOTcont
WRITE (LUPRI,'(A)')'--------------------------------------------------------&
     &-------------'
WRITE(LUPRI,*) '                     '
END SUBROUTINE PRINT_MOLECULE_AND_BASIS

!> \brief determine the label CC used in print molecule and basis
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number to write to
!> \param molecule the molecule structure
!> \param basisinfo the basis info structure
!> \param I not used
!> \param type the atomtype in the BASISINFO structure 
!> \param CC the label used in print molecule and basis
SUBROUTINE BASTYP(LUPRI,MOLECULE,BASISINFO,I,type,CC)
implicit none
TYPE(BASISSETINFO)  :: BASISINFO
TYPE(MOLECULEINFO)  :: MOLECULE
CHARACTER(len=35)   :: CC 
CHARACTER(len=1)    :: CC1
CHARACTER(len=1)    :: spdfgh(10)=(/'s','p','d','f','g','h','i','j','k','l'/) 
CHARACTER(len=2)    :: CC2 
INTEGER             :: I,J,L,IND,K,type,set,LUPRI

CC='                                   '
CC(1:1) = '['
IND=2

DO J=1,BASISINFO%ATOMTYPE(type)%nANGMOM
  L=BASISINFO%ATOMTYPE(type)%SHELL(J)%nprim
  IF(L .LT. 100 .AND. L .GT. 0)THEN
    IF(L<10)THEN
      CC1=Char(L+48)
      CC(IND:IND)=CC1
      IND=IND+1
    ELSE
      CC2=Char(L/10+48)//Char(mod(L,10)+48)
      CC(IND:IND+1)=CC2
      IND=IND+2
    ENDIF
  ELSE
!    IF(L .NE. 0) PRINT*,'ERROR DO YOU REALLY HAVE MORE THAN 99 primitives?'
  ENDIF
  IF(L .NE. 0) THEN
     IF (J.GT.10) CALL LSQUIT('Need to modify BASTYP for nANGMOM.GT.10',-1)
     CC(IND:IND)=spdfgh(J)
     IND=IND+1
  ENDIF
ENDDO

CC(IND:IND) = '|'
IND = IND + 1

DO J=1,BASISINFO%ATOMTYPE(type)%nANGMOM
  L=BASISINFO%ATOMTYPE(type)%SHELL(J)%norb
  IF(L .LT. 100 .AND. L .GT. 0)THEN
    IF(L<10)THEN
      CC1=Char(L+48)
      CC(IND:IND)=CC1
      IND=IND+1
    ELSE
      CC2=Char(L/10+48)//Char(mod(L,10)+48)
      CC(IND:IND+1)=CC2
      IND=IND+2
    ENDIF
  ELSE
!    IF(L .NE. 0) PRINT*,'ERROR DO YOU REALLY HAVE MORE THAN 99 primitives?'
  ENDIF
  IF(L .NE. 0) THEN
     CC(IND:IND)=spdfgh(J)
     IND=IND+1
  ENDIF
ENDDO

CC(IND:IND) = ']'

END SUBROUTINE BASTYP

!> \brief print the fragmentinfo
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number to write to
!> \param fragment the fragmentinfo to be printet
!> \param Label 
SUBROUTINE PRINT_FRAGMENTINFO(LUPRI,FRAGMENT,LABEL)
implicit none
TYPE(FRAGMENTINFO) :: FRAGMENT
INTEGER            :: I,LABEL
INTEGER            :: LUPRI
CHARACTER(len=7)   :: STRING(3)
WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'THE FRAGMENTINFO'
WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(2X,A,2X,I3)') 'Number of fragments    = ',FRAGMENT%numFragments
WRITE(LUPRI,'(2X,A,2X,L1)') 'Number of Orbital Sets = ',FRAGMENT%numberOrbialsSet

WRITE(LUPRI,*) '--------------------------------------------------------------------------'
WRITE(LUPRI,'(2X,A,2X,A,2X,A,2X,A,2X,A,2X,A,2X,A)')'Fragment','Basis','nPrimOrb',&
     &'nContOrb','nStartContOrb','nStartPrimOrb','Atoms'
WRITE(LUPRI,*) '--------------------------------------------------------------------------'
STRING(1) = 'Regular'
STRING(2) = 'DF-Aux '
STRING(3) = 'Huckel '
DO I=1,FRAGMENT%numFragments
   WRITE(LUPRI,'(2X,I3,6X,A7,6X,I3,7X,I3,12X,I3,12X,I3,4X,I3)') I, STRING(LABEL), &
        & FRAGMENT%nPrimOrb(I,LABEL),& 
        & FRAGMENT%nContOrb(I,LABEL), FRAGMENT%nStartContOrb(I,LABEL), &
        & FRAGMENT%nStartPrimOrb(I,LABEL), FRAGMENT%nATOMS(I)
   WRITE(LUPRI,*) I, STRING(LABEL), &
        & FRAGMENT%nPrimOrb(I,LABEL),& 
        & FRAGMENT%nContOrb(I,LABEL), FRAGMENT%nStartContOrb(I,LABEL), &
        & FRAGMENT%nStartPrimOrb(I,LABEL), FRAGMENT%nATOMS(I)
   print*,I, STRING(LABEL), &
        & FRAGMENT%nPrimOrb(I,LABEL),& 
        & FRAGMENT%nContOrb(I,LABEL), FRAGMENT%nStartContOrb(I,LABEL), &
        & FRAGMENT%nStartPrimOrb(I,LABEL), FRAGMENT%nATOMS(I)
ENDDO
END SUBROUTINE PRINT_FRAGMENTINFO

!> \brief print the blockinfo
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number to write to
!> \param blocks the blockinfo to be printet
!> \param Label 
SUBROUTINE PRINT_BLOCKINFO(LUPRI,BLOCKS,LABEL)
implicit none
TYPE(BLOCKINFO) :: BLOCKS
INTEGER         :: I
CHARACTER LABEL*(*)
INTEGER         :: LUPRI

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'THE '//LABEL//'BLOCKINFO'
WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(1X,A,I3)') 'Number of blocks    = ',BLOCKS%numBlocks

   WRITE(LUPRI,*) '--------------------------------------------------------------------'
   WRITE(LUPRI,'(1X,A5,1X,A5,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A4)')&
        &'Frag1','Frag2','nbast1','nbast2','nAtoms1','nAtoms2',&
        &'start1','start2','node'
   WRITE(LUPRI,*) '--------------------------------------------------------------------'

DO I=1,BLOCKS%numBlocks
   WRITE(LUPRI,'(1X,I5,1X,I5,1X,I6,1X,I6,1X,I6,1X,I6,2X,I6,2X,I6,1X,I4)') &
        & BLOCKS%blocks(I)%fragment1,&
        & BLOCKS%blocks(I)%fragment2, BLOCKS%blocks(I)%nbast1,& 
        & BLOCKS%blocks(I)%nbast2, BLOCKS%blocks(I)%nAtoms1, &
        & BLOCKS%blocks(I)%nAtoms2, BLOCKS%blocks(I)%startOrb1,& 
        & BLOCKS%blocks(I)%startOrb2, BLOCKS%blocks(I)%node
ENDDO

WRITE(LUPRI,*) '                     '

END SUBROUTINE PRINT_BLOCKINFO

!> \brief this subroutine is needed in the multipole moment calculation when using a buffer for writing it puts information into the real and integer buffer
!> \author S. Reine
!> \date 2010
SUBROUTINE LS_FILLBUFFER(OUTPUT,JE,M,IB,IA,J1,J2,IBTCH,IND,EXT,OD1,OD2,OD3,SPHINT,DENS)
   IMPLICIT NONE
   TYPE(INTEGRALOUTPUT) :: OUTPUT
   INTEGER      :: M, IB, IA, J1,J2,IBTCH,IND,JE
   INTEGER      :: I
   REAL(REALK)  :: EXT,OD1,OD2,OD3,SPHINT,DENS
   I = OUTPUT%IBUFI
   OUTPUT%IBUF(1,I)=JE 
   OUTPUT%IBUF(2,I)=M
   OUTPUT%IBUF(3,I)=iB
   OUTPUT%IBUF(4,I)=iA
   OUTPUT%IBUF(5,I)=J1
   OUTPUT%IBUF(6,I)=J2
   OUTPUT%IBUF(7,I)=IBTCH
   OUTPUT%IBUF(8,I)=IND
!  FIXME: if  we have the gradient we have to save more elements 
!   IF (DO_GRADIENT) THEN
!   END IF
! the real buffer
   OUTPUT%RBUF(1,I)=EXT
   OUTPUT%RBUF(2,I)=OD1
   OUTPUT%RBUF(3,I)=OD2
   OUTPUT%RBUF(4,I)=OD3
   OUTPUT%RBUF(5,I)=SPHINT
   OUTPUT%RBUF(6,I)=DENS
!  FIXME: if  we have the gradient we have to save more elements 
!   IF (DO_GRADIENT) THEN
!   END IF
   OUTPUT%IBUFI= OUTPUT%IBUFI+ 1
END SUBROUTINE LS_FILLBUFFER

!> \brief this subroutine is needed in the multipole moment calculation when using a buffer for writing it empties the integer buffer
!> \author S. Reine
!> \date 2010
SUBROUTINE LS_EMPTYIBUF(OUTPUT,IBUFFER,LUINTM)
! ATTENTION: counter is not reset !
   IMPLICIT NONE
   TYPE(INTEGRALOUTPUT) :: OUTPUT
   INTEGER     :: IBUFFER((OUTPUT%IBUFI-1)*OUTPUT%MAXBUFI)
   INTEGER     :: LUINTM, II, I
   II = OUTPUT%MAXBUFI*(OUTPUT%IBUFI-1)
   WRITE(LUINTM) II
   WRITE(LUINTM) (IBUFFER(I),I=1,II)
END SUBROUTINE LS_EMPTYIBUF

!> \brief this subroutine is needed in the multipole moment calculation when using a buffer for writing it empties the real buffer
!> \author S. Reine
!> \date 2010
SUBROUTINE LS_EMPTYRBUF(OUTPUT,RBUFFER,LUINTR)
! ATTENTION: counter is reset
   IMPLICIT NONE
   TYPE(INTEGRALOUTPUT) :: OUTPUT
   REAL(REALK) :: RBUFFER((OUTPUT%IBUFI-1)*OUTPUT%MAXBUFR)
   INTEGER     :: LUINTR, II, I
   II = OUTPUT%MAXBUFR*(OUTPUT%IBUFI-1)
   WRITE(LUINTR) II
   WRITE(LUINTR) (RBUFFER(I),I=1,II)
   OUTPUT%IBUFI = 1
END SUBROUTINE LS_EMPTYRBUF

!SUBROUTINE LS_EMPTYBUF(OUTPUT,IBUFFER,RBUFFER,LUINTM,LUINTMR)
! this subroutine is needed in the multipole moment calculation when using a buffer for writing
! it empties the integer and the real buffer
!   IMPLICIT NONE
!   TYPE(INTEGRALOUTPUT) :: OUTPUT
!   INTEGER     :: IBUFFER(OUTPUT%MAXBUFI*(OUTPUT%IBUFI-1))
!!   INTEGER     :: LUINTM, LUINTMR, IBUFI, IBUFR, I
!   REAL(REALK) :: RBUFFER(OUTPUT%MAXBUFR*(OUTPUT%IBUFI-1))
!   IBUFI = OUTPUT%MAXBUFI*(OUTPUT%IBUFI-1)
!   IBUFR = OUTPUT%MAXBUFR*(OUTPUT%IBUFI-1)
!   WRITE(LUINTM) IBUFI
!   WRITE(LUINTM) (IBUFFER(I),I=1,IBUFI)
!   WRITE(LUINTMR) IBUFR
!   WRITE(LUINTMR) (RBUFFER(I),I=1,IBUFR)
!END SUBROUTINE LS_EMPTYBUF

!> \brief this subroutine is needed in the multipole moment calculation when using a buffer for writing it puts nuclear information into the buffer
!> \author S. Reine
!> \date 2010
SUBROUTINE LS_FILLNUCBUF(OUTPUT,CHARGE,X,Y,Z)
   IMPLICIT NONE 
   TYPE(INTEGRALOUTPUT) :: OUTPUT
   REAL(REALK)          :: CHARGE, X,Y,Z
   INTEGER              :: I
   I = OUTPUT%IBUFN
   OUTPUT%NBUF(1,I) = CHARGE
   OUTPUT%NBUF(2,I) = X
   OUTPUT%NBUF(3,I) = Y
   OUTPUT%NBUF(4,I) = Z
   OUTPUT%IBUFN = OUTPUT%IBUFN+1
END SUBROUTINE LS_FILLNUCBUF

!> \brief this subroutine is needed in the multipole moment calculation when using a buffer for writing it empties the nuclear position buffer
!> \author S. Reine
!> \date 2010
SUBROUTINE LS_EMPTYNUCBUF(OUTPUT,RBUFFER,LUINTMR)
   IMPLICIT NONE
   TYPE(INTEGRALOUTPUT) :: OUTPUT
   REAL(REALK) :: RBUFFER((OUTPUT%IBUFN-1)*OUTPUT%MAXBUFN)
   INTEGER :: LUINTMR, II,I
   II = (OUTPUT%IBUFN-1)*OUTPUT%MAXBUFN
   WRITE(LUINTMR) II
   WRITE(LUINTMR) (RBUFFER(I),I=1,II)   
   OUTPUT%IBUFN = 1
END SUBROUTINE LS_EMPTYNUCBUF

!> \brief this subroutine is needed in the multipole moment calculation when using a buffer for writing it initializes the buffer arrays
!> \author S. Reine
!> \date 2010
!>
!> and ATTENTION: passes the MMBUFLEN, MAXBUFN, MAXBUFR, MAXBUFI parameters 
!>                to the COMMON block variables in cbifmm.h
!>                which are needed in the fmm routines
!>
SUBROUTINE LS_INITMMBUF(OUTPUT,NDER)
   IMPLICIT NONE
   TYPE(INTEGRALOUTPUT) :: OUTPUT
   INTEGER :: NDER
   OUTPUT%MMBUFLEN  = 2000
   OUTPUT%IBUFI  = 1
   OUTPUT%IBUFN  = 1
   OUTPUT%MAXBUFN = 4
   IF (NDER .EQ. 1)  THEN
      OUTPUT%MAXBUFR   = 6
      OUTPUT%MAXBUFI   = 8
      ALLOCATE(OUTPUT%IBUF(OUTPUT%MAXBUFI,OUTPUT%MMBUFLEN))
      ALLOCATE(OUTPUT%RBUF(OUTPUT%MAXBUFR,OUTPUT%MMBUFLEN))
      ALLOCATE(OUTPUT%NBUF(OUTPUT%MAXBUFN,OUTPUT%MMBUFLEN))
      CALL INITBUFMEM(OUTPUT%IBUF,OUTPUT%RBUF,OUTPUT%NBUF,OUTPUT%MAXBUFI,OUTPUT%MAXBUFR,OUTPUT%MAXBUFN,OUTPUT%MMBUFLEN)
   ELSEIF (NDER .EQ. 2) THEN
      OUTPUT%MAXBUFR   = 12
      OUTPUT%MAXBUFI   = 16
      WRITE(*,*) 'USE OF BUFFER NOT IMPLEMENTED YET FOR DERIVATIVES'
      WRITE(*,*) 'SET USEBUFMM = .FALSE.'
      OUTPUT%USEBUFMM=.FALSE.
      !ALLOCATE(OUTPUT%IBUF(OUTPUT%MAXBUFI,OUTPUT%MMBUFLEN))
      !ALLOCATE(OUTPUT%RBUF(OUTPUT%MAXBUFR,OUTPUT%MMBUFLEN))
      !ALLOCATE(OUTPUT%NBUF(OUTPUT%MAXBUFN,OUTPUT%MMBUFLEN))
   ELSE 
      CALL LSQUIT('UNDEFINED BUFFERLENGTH FOR MM-BUFFER FOR HIGHER DERIVATIVES?',-1)
   ENDIF
   CALL SETMMBUFINFO(OUTPUT%USEBUFMM,OUTPUT%MMBUFLEN,OUTPUT%MAXBUFI,OUTPUT%MAXBUFR,OUTPUT%MAXBUFN)

END SUBROUTINE LS_INITMMBUF

!> \brief this subroutine is needed in the multipole moment calculation when using a buffer for writing it empties the nuclear position buffe
!> \author S. Reine
!> \date 2010
SUBROUTINE LS_FREEMMBUF(OUTPUT)
   IMPLICIT NONE
   TYPE(INTEGRALOUTPUT) :: OUTPUT
   DEALLOCATE(OUTPUT%IBUF)
   DEALLOCATE(OUTPUT%RBUF)
   DEALLOCATE(OUTPUT%NBUF)
END SUBROUTINE LS_FREEMMBUF

!> \brief 
!> \author S. Reine
!> \date 2010
SUBROUTINE INITBUFMEM(IBUFFER,RBUFFER,NBUFFER,MAXI,MAXR,MAXN,MAXL)
   IMPLICIT NONE
   INTEGER IBUFFER(MAXI*MAXL), I,MAXI,MAXR,MAXN,MAXL
   REAL(REALK) RBUFFER(MAXR*MAXL), NBUFFER(MAXN*MAXL)
   DO I=1, MAXR*MAXL
      RBUFFER(I) = 0.0D0
   ENDDO
   DO I=1, MAXI*MAXL
      IBUFFER(I) = 0
   ENDDO
   DO I=1, MAXN*MAXL
      NBUFFER(I) = 0.0D0
   ENDDO
END SUBROUTINE INITBUFMEM

!> \brief retrieve the output from the setting and put it into a single matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output
!> \param setting contains the output in lstensor format 
!> \param singleMAT the output matrix 
SUBROUTINE retrieve_output_mat_single(lupri,setting,singleMAT)
implicit none
TYPE(LSSETTING) :: setting
TYPE(MATRIX)    :: singleMAT 
integer   :: lupri

call Build_matrix_from_lstensor(lupri,setting%Output%ResultMat2,singleMAT)
call lstensor_free(setting%Output%ResultMat2)

END SUBROUTINE retrieve_output_mat_single

!> \brief retrieve the output from the setting and put it into a single lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output
!> \param setting contains the output in lstensor format 
!> \param tensor the output lstensor
SUBROUTINE retrieve_output_lstensor(lupri,setting,tensor)
implicit none
TYPE(LSSETTING) :: setting
TYPE(lstensor)    :: tensor
integer   :: lupri

call copy_lstensor_to_lstensor(setting%Output%ResultMat2,tensor)
call lstensor_free(setting%Output%ResultMat2)

END SUBROUTINE retrieve_output_lstensor

!> \brief retrieve the output from the setting and put it into an array matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output
!> \param setting contains the output in lstensor format 
!> \param arrayMAT the output matrix 
SUBROUTINE retrieve_output_mat_array(lupri,setting,arrayMAT)
implicit none
TYPE(LSSETTING) :: setting
TYPE(MATRIX)    :: arrayMAT(:)
integer   :: lupri

call Build_matrix_from_lstensor(lupri,setting%Output%ResultMat2,arrayMAT)
call lstensor_free(setting%Output%ResultMat2)

END SUBROUTINE retrieve_output_mat_array

!> \brief retrieve the output from the setting and put it into an 5 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output
!> \param setting contains the output in lstensor format 
!> \param fullmat the output 5 dim array
SUBROUTINE retrieve_output_5dim(lupri,setting,fullMAT)
implicit none
TYPE(LSSETTING) :: setting
REAL(REALK)      :: fullMAT(:,:,:,:,:)
integer   :: lupri,ndim1,ndim2,ndim3,ndim4,ndim5

ndim1 = setting%Output%ResultMat2%nbast1
ndim2 = setting%Output%ResultMat2%nbast2
ndim3 = setting%Output%ResultMat2%nbast3
ndim4 = setting%Output%ResultMat2%nbast4
ndim5 = setting%Output%ResultMat2%nmat
Call Build_full_5dim_from_lstensor(setting%Output%ResultMat2,&
     &fullMAT,ndim1,ndim2,ndim3,ndim4,ndim5)
call lstensor_free(setting%Output%ResultMat2)

END SUBROUTINE retrieve_output_5dim

!> \brief retrieve the output from the setting and put it into an 3 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output
!> \param setting contains the output in lstensor format 
!> \param fullmat the output 3 dim array
SUBROUTINE retrieve_output_3dim(lupri,setting,fullMAT)
implicit none
TYPE(LSSETTING) :: setting
REAL(REALK)      :: fullMAT(:,:,:)
integer   :: lupri,ndim1,ndim2,ndim3,ndim4,ndim5,n1,n2,n3
REAL(REALK),pointer      :: fullMATtmp(:,:,:,:,:)

ndim1 = setting%Output%ResultMat2%nbast1
ndim2 = setting%Output%ResultMat2%nbast2
ndim3 = setting%Output%ResultMat2%nbast3
ndim4 = setting%Output%ResultMat2%nbast4
ndim5 = setting%Output%ResultMat2%nmat
n1 = size(fullMAT(:,1,1))
n2 = size(fullMAT(1,:,1))
n3 = size(fullMAT(1,1,:))
IF (ndim1*ndim2*ndim3*ndim4*ndim5.NE.n1*n2*n3) THEN
  WRITE(LUPRI,'(1X,A)') 'Error in retrieve_output_3dim. Mismatching dimensions'
  CALL LSQUIT('Error in retrieve_output_3dim. Mismatching dimensions',-1)
ENDIF
call mem_alloc(fullMATtmp,ndim1,ndim2,ndim3,ndim4,ndim5)
Call Build_full_5dim_from_lstensor(setting%Output%ResultMat2,&
     &fullMATtmp,ndim1,ndim2,ndim3,ndim4,ndim5)
call lstensor_free(setting%Output%ResultMat2)
call dcopy(ndim1*ndim2*ndim3*ndim4*ndim5,fullMATtmp,1,fullMAT,1)
call mem_dealloc(fullMATtmp)

END SUBROUTINE retrieve_output_3dim
#if 0
SUBROUTINE retrieve_output_3dim(lupri,setting,fullMAT)
implicit none
TYPE(LSSETTING) :: setting
REAL(REALK)      :: fullMAT(:,:,:)
integer   :: lupri,I,J

call Build_full_3dim_from_lstensor(setting%Output%ResultMat2,&
     & fullMAT,setting%Output%ndim(1),setting%Output%ndim(2),&
     & setting%Output%ndim(5))
call lstensor_free(setting%Output%ResultMat2)

END SUBROUTINE retrieve_output_3dim
#endif

!> \brief retrieve the output from the setting and put it into an 2 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output
!> \param setting contains the output in lstensor format 
!> \param fullmat the output 2 dim array
SUBROUTINE retrieve_output_2dim(lupri,setting,fullMAT)
implicit none
TYPE(LSSETTING) :: setting
REAL(REALK)     :: fullMAT(:,:)
integer         :: lupri
!
integer             :: I,J,ndim1,ndim2,ndim3,ndim4,ndim5,n1,n2
Real(realk),pointer :: fullMATtmp(:,:,:,:,:)

IF(setting%output%resultmat2%gradienttensor)THEN
   call build_gradient_from_gradientlstensor(setting%output%resultmat2,fullmat,&
        &setting%output%resultmat2%natom1,setting%output%resultmat2%nmat,lupri)
   call lstensor_free(setting%Output%ResultMat2)
ELSE
  ndim1 = setting%Output%ResultMat2%nbast1
  ndim2 = setting%Output%ResultMat2%nbast2
  ndim3 = setting%Output%ResultMat2%nbast3
  ndim4 = setting%Output%ResultMat2%nbast4
  ndim5 = setting%Output%ResultMat2%nmat
  n1 = size(fullMAT(:,1))
  n2 = size(fullMAT(1,:))
  IF (ndim1*ndim2*ndim3*ndim4*ndim5.NE.n1*n2) THEN
    WRITE(LUPRI,'(1X,A)') 'Error in retrieve_output_2dim. Mismatching dimensions'
    CALL LSQUIT('Error in retrieve_output_2dim. Mismatching dimensions',-1)
  ENDIF
  call mem_alloc(fullMATtmp,ndim1,ndim2,ndim3,ndim4,ndim5)
  Call Build_full_5dim_from_lstensor(setting%Output%ResultMat2,&
       &fullMATtmp,ndim1,ndim2,ndim3,ndim4,ndim5)
  call lstensor_free(setting%Output%ResultMat2)
  call dcopy(ndim1*ndim2*ndim3*ndim4*ndim5,fullMATtmp,1,fullMAT,1)
  call mem_dealloc(fullMATtmp)
ENDIF

!IF(setting%output%resultmat2%gradienttensor)THEN
!   call build_gradient_from_gradientlstensor(setting%output%resultmat2,fullmat,&
!        &setting%output%resultmat2%natom1,setting%output%resultmat2%nmat,lupri)
!   call lstensor_free(setting%Output%ResultMat2)
!ELSE
!   call Build_full_2dim_from_lstensor(setting%Output%ResultMat2,&
!        &fullMAT,setting%Output%ndim(1),setting%Output%ndim(2))
!   call lstensor_free(setting%Output%ResultMat2)
!ENDIF

END SUBROUTINE retrieve_output_2dim

!> \brief print the lsintscheme structure
!> \author T. Kjaergaard
!> \date 2010
!> \param scheme the lsintschem to be printet
!> \param iunit the logical unit number for the output
SUBROUTINE printLsintScheme(scheme,IUNIT)
implicit none
TYPE(LSINTSCHEME),INTENT(IN) :: scheme
INTEGER,INTENT(IN)           :: IUNIT

WRITE(IUNIT,'(3X,A22,L7)') 'CFG_LSDALTON          ', scheme%CFG_LSDALTON          
WRITE(IUNIT,'(3X,A22,L7)') 'DOPASS                ', scheme%DOPASS                
WRITE(IUNIT,'(3X,A22,L7)') 'DENSFIT               ', scheme%DENSFIT               
WRITE(IUNIT,'(3X,A22,I7)') 'AOPRINT               ', scheme%AOPRINT              
WRITE(IUNIT,'(3X,A22,I7)') 'INTPRINT              ', scheme%INTPRINT              
WRITE(IUNIT,'(3X,A22,L7)') 'JENGINE               ', scheme%JENGINE               
WRITE(IUNIT,'(3X,A22,I7)') 'FTUVmaxprim           ', scheme%FTUVmaxprim           
WRITE(IUNIT,'(3X,A22,I7)') 'maxpasses             ', scheme%maxpasses             
WRITE(IUNIT,'(3X,A22,L7)') 'FMM                   ', scheme%FMM                   
WRITE(IUNIT,'(3X,A22,L7)') 'LINK                  ', scheme%LINK                  
WRITE(IUNIT,'(3X,A22,L7)') 'DALINK                ', scheme%DALINK                
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUGOVERLAP          ', scheme%DEBUGOVERLAP
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUG4CENTER          ', scheme%DEBUG4CENTER
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUG4CENTER_ERI      ', scheme%DEBUG4CENTER_ERI
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUGCCFRAGMENT       ', scheme%DEBUGCCFRAGMENT
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUGKINETIC          ', scheme%DEBUGKINETIC          
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUGNUCREP           ', scheme%DEBUGNUCREP           
WRITE(IUNIT,'(3X,A22,L7)') 'DO4CENTERERI          ', scheme%DO4CENTERERI          
WRITE(IUNIT,'(3X,A22,L7)') 'OVERLAP_DF_J          ', scheme%OVERLAP_DF_J          
WRITE(IUNIT,'(3X,A22,L7)') 'TIMINGS               ', scheme%TIMINGS               
WRITE(IUNIT,'(3X,A22,L7)') 'nonSphericalETUV      ', scheme%nonSphericalETUV      
WRITE(IUNIT,'(3X,A22,L7)') 'ETUVinsideLOOP        ', scheme%ETUVinsideLOOP        
WRITE(IUNIT,'(3X,A22,L7)') 'HIGH_RJ000_ACCURACY   ', scheme%HIGH_RJ000_ACCURACY   
WRITE(IUNIT,'(3X,A22,L7)') 'OrderAngPrim          ', scheme%OrderAngPrim          
WRITE(IUNIT,'(3X,A22,I7)') 'MM_LMAX               ', scheme%MM_LMAX               
WRITE(IUNIT,'(3X,A22,I7)') 'MM_TLMAX              ', scheme%MM_TLMAX              
WRITE(IUNIT,'(3X,A22,G14.2)') 'MM_SCREEN             ', scheme%MM_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'NO_MMFILES            ', scheme%NO_MMFILES            
WRITE(IUNIT,'(3X,A22,L7)') 'MM_NO_ONE             ', scheme%MM_NO_ONE             
WRITE(IUNIT,'(3X,A22,L7)') 'CREATED_MMFILES       ', scheme%CREATED_MMFILES       
WRITE(IUNIT,'(3X,A22,L7)') 'USEBUFMM              ', scheme%USEBUFMM              
WRITE(IUNIT,'(3X,A22,I7)') 'MMunique_ID1          ', scheme%MMunique_ID1          
WRITE(IUNIT,'(3X,A22,L7)') 'AUXBASIS              ', scheme%AUXBASIS              
WRITE(IUNIT,'(3X,A22,L7)') 'NOFAMILY              ', scheme%NOFAMILY              
WRITE(IUNIT,'(3X,A22,L7)') 'DoCartesian           ', scheme%DoCartesian           
WRITE(IUNIT,'(3X,A22,L7)') 'DoSpherical           ', scheme%DoSpherical           
WRITE(IUNIT,'(3X,A22,L7)') 'UNCONT                ', scheme%UNCONT                
WRITE(IUNIT,'(3X,A22,L7)') 'NOSEGMENT             ', scheme%NOSEGMENT             
WRITE(IUNIT,'(3X,A22,L7)') 'DO3CENTEROVL          ', scheme%DO3CENTEROVL          
WRITE(IUNIT,'(3X,A22,L7)') 'DO2CENTERERI          ', scheme%DO2CENTERERI          
WRITE(IUNIT,'(3X,A22,I7)') 'CARMOM                ', scheme%CARMOM                
WRITE(IUNIT,'(3X,A22,I7)') 'SPHMOM                ', scheme%SPHMOM                
WRITE(IUNIT,'(3X,A22,L7)') 'MIXEDOVERLAP          ', scheme%MIXEDOVERLAP          
WRITE(IUNIT,'(3X,A22,L7)') 'TEST_NDMAT_COULOMB    ', scheme%TEST_NDMAT_COULOMB    
WRITE(IUNIT,'(3X,A22,L7)') 'TEST_NDMAT_EXCHANGE   ', scheme%TEST_NDMAT_EXCHANGE   
WRITE(IUNIT,'(3X,A22,G14.2)') 'THRESHOLD             ', scheme%THRESHOLD             
WRITE(IUNIT,'(3X,A22,G14.2)') 'CS_THRESHOLD          ', scheme%CS_THRESHOLD          
WRITE(IUNIT,'(3X,A22,G14.2)') 'OE_THRESHOLD          ', scheme%OE_THRESHOLD          
WRITE(IUNIT,'(3X,A22,G14.2)') 'PS_THRESHOLD          ', scheme%PS_THRESHOLD          
WRITE(IUNIT,'(3X,A22,G14.2)') 'OD_THRESHOLD          ', scheme%OD_THRESHOLD          
WRITE(IUNIT,'(3X,A22,G14.2)') 'Integralthreshold     ', scheme%Integralthreshold     
WRITE(IUNIT,'(3X,A22,L7)') 'CS_SCREEN             ', scheme%CS_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'OE_SCREEN             ', scheme%OE_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'GAB_ON_FILE           ', scheme%GAB_ON_FILE           
WRITE(IUNIT,'(3X,A22,L7)') 'WRITE_GAB_TO_FILE     ', scheme%WRITE_GAB_TO_FILE     
WRITE(IUNIT,'(3X,A22,L7)') 'PS_SCREEN             ', scheme%PS_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'PRIM_GAB_ON_FILE      ', scheme%PRIM_GAB_ON_FILE      
WRITE(IUNIT,'(3X,A22,L7)') 'WRITE_PRIM_GAB_TO_FILE', scheme%WRITE_PRIM_GAB_TO_FILE
WRITE(IUNIT,'(3X,A22,L7)') 'PS_DEBUG              ', scheme%PS_DEBUG              
WRITE(IUNIT,'(3X,A22,L7)') 'OD_SCREEN             ', scheme%OD_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'FRAGMENT              ', scheme%FRAGMENT              
WRITE(IUNIT,'(3X,A22,I7)') 'numAtomsPerFragment   ', scheme%numAtomsPerFragment   
WRITE(IUNIT,'(3X,A22,L7)') 'PUREFMM               ', scheme%PUREFMM               
WRITE(IUNIT,'(3X,A22,I7)') 'LU_LUINTM             ', scheme%LU_LUINTM             
WRITE(IUNIT,'(3X,A22,I7)') 'LU_LUINTR             ', scheme%LU_LUINTR             
WRITE(IUNIT,'(3X,A22,L7)')'DECPACKED                ', scheme%DECPACKED

WRITE(IUNIT,'(3X,A22,L7)') 'CAM                   ', scheme%CAM                   
WRITE(IUNIT,'(3X,A22,G14.2)') 'CAMalpha              ', scheme%CAMalpha              
WRITE(IUNIT,'(3X,A22,G14.2)') 'CAMbeta               ', scheme%CAMbeta               
WRITE(IUNIT,'(3X,A22,G14.2)') 'CAMmu                 ', scheme%CAMmu                 
WRITE(IUNIT,'(3X,A22,G14.2)') 'exchangeFactor        ', scheme%exchangeFactor        
WRITE(IUNIT,'(3X,A22,L7)')'GRDONE                ', scheme%GRDONE 
WRITE(IUNIT,'(3X,A22,I7)')'ITERATIONS            ', scheme%ITERATIONS
WRITE(IUNIT,'(3X,A22,I7)')'ANGMIN                ', scheme%ANGMIN 
WRITE(IUNIT,'(3X,A22,I7)')'HRDNES                ', scheme%HRDNES 
WRITE(IUNIT,'(3X,A22,I7)')'TURBO                 ', scheme%TURBO   
WRITE(IUNIT,'(3X,A22,G14.2)') 'RADINT                ', scheme%RADINT
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTELS                ', scheme%DFTELS
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTHR0                ', scheme%DFTHR0
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTHRI                ', scheme%DFTHRI
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTHRL                ', scheme%DFTHRL
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTHRL                ', scheme%DFTHRL
WRITE(IUNIT,'(3X,A22,G14.2)') 'RHOTHR                ', scheme%RHOTHR
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTIPT                ', scheme%DFTIPT
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTBR1                ', scheme%DFTBR1
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTBR2                ', scheme%DFTBR2
WRITE(IUNIT,'(3X,A22,L7)')'NOPRUN                ', scheme%NOPRUN 
WRITE(IUNIT,'(3X,A22,L7)')'DFTASC                ', scheme%DFTASC
WRITE(IUNIT,'(3X,A22,L7)')'DFTPOT                ', scheme%DFTPOT
WRITE(IUNIT,'(3X,A22,L7)')'DISPER                ', scheme%DODISP
WRITE(IUNIT,'(3X,A22,L7)')'DISPDONE              ', scheme%DISPDONE
WRITE(IUNIT,'(3X,A22,L7)')'DFTADD                ', scheme%DFTADD

END SUBROUTINE printLsintScheme

!> \brief copy the lssetting structure
!> \author T. Kjaergaard
!> \date 2010
!> \param newsetting the new setting structure
!> \param oldsetting the original setting structure
!> \param lupri the logical unit number for the output
SUBROUTINE copy_setting(newsetting,oldsetting,lupri)
implicit none
type(lssetting),intent(in)    :: oldsetting
type(lssetting),intent(inout) :: newsetting
integer :: lupri
!
integer :: I,nAO,ndmat,dim1,dim2,dim3

newsetting = oldsetting 
newsetting%DO_DFT = oldsetting%DO_DFT
newsetting%nAO = oldsetting%nAO
nAO = newsetting%nAO
NULLIFY(newsetting%Molecule)
NULLIFY(newsetting%Basis)
NULLIFY(newsetting%Fragment)
NULLIFY(newsetting%Batchindex)
NULLIFY(newsetting%Batchdim)
NULLIFY(newsetting%molID)

ALLOCATE(newsetting%Molecule(nAO))
ALLOCATE(newsetting%Fragment(nAO))
ALLOCATE(newsetting%Basis(nAO))
ALLOCATE(newsetting%Batchindex(nAO))
ALLOCATE(newsetting%Batchdim(nAO))
ALLOCATE(newsetting%molID(nAO))

do I = 1,nAO
   nullify(newsetting%Molecule(I)%p)
   allocate(newsetting%Molecule(I)%p)
   call copy_molecule(oldsetting%Molecule(I)%p,newsetting%Molecule(I)%p,lupri)
   nullify(newsetting%Fragment(I)%p)
   allocate(newsetting%Fragment(I)%p)
   call copy_molecule(oldsetting%Fragment(I)%p,newsetting%Fragment(I)%p,lupri)
   nullify(newsetting%Basis(I)%p)
   allocate(newsetting%Basis(I)%p)
   call copy_basissetinfo(oldsetting%basis(I)%p%REGULAR,newsetting%basis(I)%p%REGULAR) 
   call copy_basissetinfo(oldsetting%basis(I)%p%AUXILIARY,newsetting%basis(I)%p%AUXILIARY) 
   newsetting%Batchindex(I) = oldsetting%Batchindex(I)
   newsetting%Batchdim(I) = oldsetting%Batchdim(I)
   newsetting%molID(I) = oldsetting%molID(I) 
enddo

NULLIFY(newSETTING%molBuild)
NULLIFY(newSETTING%basBuild)
NULLIFY(newSETTING%fragBuild)
ALLOCATE(newSETTING%molBuild(nAO))
ALLOCATE(newSETTING%basBuild(nAO))
ALLOCATE(newSETTING%fragBuild(nAO))

newsetting%molBuild = .TRUE.!oldSETTING%molBuild
newsetting%basBuild = .TRUE.!oldSETTING%basBuild
newsetting%fragBuild = .TRUE.!oldSETTING%fragBuild

NULLIFY(newSETTING%sameMOL)
NULLIFY(newSETTING%sameBAS)
NULLIFY(newSETTING%sameFRAG)
ALLOCATE(newSETTING%sameMOL(nAO,nAO))
ALLOCATE(newSETTING%sameBAS(nAO,nAO))
ALLOCATE(newSETTING%sameFRAG(nAO,nAO))

newsetting%sameMOL = oldSETTING%sameMOL
newsetting%sameBAS = oldSETTING%sameBAS
newsetting%sameFRAG = oldSETTING%sameFRAG

NULLIFY(newSETTING%DmatLHS)
NULLIFY(newSETTING%DmatRHS)
NULLIFY(newSETTING%DfullLHS)
NULLIFY(newSETTING%DfullRHS)
newSETTING%nDmatLHS = oldSETTING%nDmatLHS
newSETTING%nDmatRHS = oldSETTING%nDmatRHS

newSETTING%LHSdmat = oldSETTING%LHSdmat
newSETTING%RHSdmat = oldSETTING%RHSdmat
NULLIFY(newSETTING%DsymLHS)
NULLIFY(newSETTING%DsymRHS)
 
IF(newSETTING%LHSdmat)THEN
  ndmat = newSETTING%nDmatLHS
  nullify(newsetting%DmatLHS)
  allocate(newsetting%DmatLHS(ndmat))
  do I=1,ndmat
     dim1 = oldsetting%DmatLHS(I)%p%nrow
     dim2 = oldsetting%DmatLHS(I)%p%ncol
     call mat_print(oldsetting%DmatLHS(I)%p,1,dim1,1,dim2,lupri)
     nullify(newsetting%DmatLHS(I)%p)
     newsetting%DmatLHS(I)%p = oldsetting%DmatLHS(I)%p
     call mat_print(newsetting%DmatLHS(I)%p,1,dim1,1,dim2,lupri)
  enddo
  call mem_alloc(newsetting%DsymLHS,ndmat)
  newsetting%DsymLHS = oldsetting%DsymLHS
ENDIF
IF(newSETTING%RHSdmat)THEN
  ndmat = newSETTING%nDmatRHS
  nullify(newsetting%DmatRHS)
  allocate(newsetting%DmatRHS(newSETTING%nDmatRHS))
  do I=1,newSETTING%nDmatRHS
     dim1 = oldsetting%DmatRHS(I)%p%nrow
     dim2 = oldsetting%DmatRHS(I)%P%ncol
     call mat_print(oldsetting%DmatRHS(I)%p,1,dim1,1,dim2,lupri)
     nullify(newsetting%DmatRHS(I)%p)
     newsetting%DmatRHS(I)%p = oldsetting%DmatRHS(I)%p
     call mat_print(newsetting%DmatRHS(I)%p,1,dim1,1,dim2,lupri)
  enddo
  call mem_alloc(newsetting%DsymRHS,ndmat)
  newsetting%DsymRHS = oldsetting%DsymRHS 
ENDIF

newSETTING%LHSdalloc = .FALSE.
newSETTING%RHSdalloc = .FALSE.

newSETTING%LHSdfull = oldSETTING%LHSdfull 
IF(newSETTING%LHSdfull)THEN
   ndmat = newSETTING%nDmatLHS
   nullify(newsetting%DfullLHS)
   dim1 = SIZE(oldSETTING%DfullLHS, 1)  
   dim2 = SIZE(oldSETTING%DfullLHS, 2)  
   dim3 = SIZE(oldSETTING%DfullLHS, 3)  
   allocate(newSETTING%DfullLHS(dim1,dim2,dim3))
   call output(oldSETTING%DfullLHS(:,:,1),1,dim1,1,dim2,dim1,dim2,1,6)
   newSETTING%DfullLHS = oldSETTING%DfullLHS 
   newSETTING%LHSdalloc = .TRUE.
   IF(.NOT.ASSOCIATED(newsetting%DsymLHS))THEN
      call mem_alloc(newsetting%DsymLHS,ndmat)
      newsetting%DsymLHS = oldsetting%DsymLHS
   ENDIF
ENDIF 

newSETTING%RHSdfull = oldSETTING%RHSdfull
IF(newSETTING%RHSdfull)THEN
   ndmat = newSETTING%nDmatRHS
   nullify(newsetting%DfullRHS)
   dim1 = SIZE(oldSETTING%DfullRHS, 1)  
   dim2 = SIZE(oldSETTING%DfullRHS, 2)  
   dim3 = SIZE(oldSETTING%DfullRHS, 3)  
   allocate(newSETTING%DfullRHS(dim1,dim2,dim3))
   call output(oldSETTING%DfullLHS(:,:,1),1,dim1,1,dim2,dim1,dim2,1,6)
   newSETTING%DfullRHS = oldSETTING%DfullRHS 
   newSETTING%RHSdalloc = .TRUE.
   IF(.NOT.ASSOCIATED(newsetting%DsymRHS))THEN
      call mem_alloc(newsetting%DsymRHS,ndmat)
      newsetting%DsymRHS = oldsetting%DsymRHS 
   ENDIF
ENDIF 

newSETTING%LHSdmatAOindex1 = oldSETTING%LHSdmatAOindex1 
newSETTING%LHSdmatAOindex2 = oldSETTING%LHSdmatAOindex2 
newSETTING%RHSdmatAOindex1 = oldSETTING%RHSdmatAOindex1
newSETTING%RHSdmatAOindex2 = oldSETTING%RHSdmatAOindex2
newSETTING%numFragments = oldSETTING%numFragments
newSETTING%numNodes = oldSETTING%numNodes
newSETTING%node  = oldSETTING%node 
newsetting%SCHEME = oldsetting%SCHEME

call io_init(newsetting%IO)

END SUBROUTINE copy_setting

!> \brief mpi copy the lssetting structure
!> \author T. Kjaergaard
!> \date 2010
!> \param setting the setting structure to broadcast
SUBROUTINE mpicopy_setting(setting)
use lsmpi_mod
use infpar_module
implicit none
#ifdef VAR_LSMPI
!include 'mpif.h'
#endif
type(lssetting),intent(inout) :: setting
integer :: lupri
!
#ifdef VAR_LSMPI
integer :: I,nAO,ndmat,dim1,dim2,dim3,ierr,master,iAO
logical :: SLAVE

 call  MPI_COMM_RANK( MPI_COMM_WORLD,SETTING%node, ierr )
 call  MPI_COMM_SIZE( MPI_COMM_WORLD,SETTING%numNodes, ierr )

Master       = infpar%master
Setting%node = infpar%mynum
SLAVE = SETTING%node.ne.0

CALL LS_MPIBCAST(setting%numFragments,Master)
CALL LS_MPIBCAST(setting%nAO,Master)
nAO = setting%nAO

IF(SLAVE)THEN
   NULLIFY(setting%Molecule)
   NULLIFY(setting%Basis)
   NULLIFY(setting%Fragment)
   NULLIFY(setting%Batchindex)
   NULLIFY(setting%Batchdim)
   NULLIFY(setting%molID)
   ALLOCATE(setting%Molecule(nAO))
   ALLOCATE(setting%Fragment(nAO))
   ALLOCATE(setting%Basis(nAO))
   ALLOCATE(setting%Batchindex(nAO))
   ALLOCATE(setting%Batchdim(nAO))
   ALLOCATE(setting%molID(nAO))
ENDIF

do I = 1,nAO
   IF(SLAVE)THEN
      nullify(setting%Molecule(I)%p)
      allocate(setting%Molecule(I)%p)
      nullify(setting%Fragment(I)%p)
      allocate(setting%Fragment(I)%p)
      nullify(setting%Basis(I)%p)
      allocate(setting%Basis(I)%p)
   ENDIF
   call mpicopy_molecule(setting%Molecule(I)%p,Slave,Master)
!   call mpicopy_molecule(setting%Fragment(I)%p,Slave,Master)
   call mpicopy_basissetinfo(setting%basis(I)%p%REGULAR,Slave,Master) 
   call mpicopy_basissetinfo(setting%basis(I)%p%AUXILIARY,Slave,Master) 
   call LS_MPIBCAST(setting%Batchindex,nAO,Master)
   call LS_MPIBCAST(setting%Batchdim,nAO,Master)
   call LS_MPIBCAST(setting%molID,nAO,Master)
enddo

IF(SLAVE)THEN
   NULLIFY(SETTING%molBuild)
   NULLIFY(SETTING%basBuild)
   NULLIFY(SETTING%fragBuild)
   ALLOCATE(SETTING%molBuild(nAO))
   ALLOCATE(SETTING%basBuild(nAO))
   ALLOCATE(SETTING%fragBuild(nAO))
   setting%molBuild = .TRUE.
   setting%basBuild = .TRUE.
   setting%fragBuild = .FALSE.
   DO iAO=1,nAO
     setting%Fragment(iAO)%p => setting%Molecule(iAO)%p
   ENDDO
   NULLIFY(SETTING%sameMOL)
   NULLIFY(SETTING%sameBAS)
   NULLIFY(SETTING%sameFRAG)
   ALLOCATE(SETTING%sameMOL(nAO,nAO))
   ALLOCATE(SETTING%sameBAS(nAO,nAO))
   ALLOCATE(SETTING%sameFRAG(nAO,nAO))
ENDIF

call LS_MPIBCAST(setting%sameMOL,nAO,nAO,Master)
call LS_MPIBCAST(setting%sameBas,nAO,nAO,Master)
call LS_MPIBCAST(setting%sameFrag,nAO,nAO,Master)

IF(SLAVE)THEN
   NULLIFY(setting%DmatLHS)
   NULLIFY(setting%DmatRHS)
   NULLIFY(setting%DfullLHS)
   NULLIFY(setting%DfullRHS)
   NULLIFY(setting%DsymLHS)
   NULLIFY(setting%DsymRHS)
ENDIF

call LS_MPIBCAST(setting%nDmatLHS,Master)
call LS_MPIBCAST(setting%nDmatRHS,Master)
call LS_MPIBCAST(setting%LHSdmat,Master)
call LS_MPIBCAST(setting%RHSdmat,Master)
 
IF(SETTING%LHSdmat)THEN
  ndmat = SETTING%nDmatLHS
  IF(SLAVE)THEN
     nullify(setting%DmatLHS)
     allocate(setting%DmatLHS(ndmat))
  ENDIF
  do I=1,ndmat
     IF(SLAVE)NULLIFY(setting%DmatLHS(I)%p)
     call mat_mpicopy(setting%DmatLHS(I)%p,Slave,Master)
  enddo
  IF(SLAVE)call mem_alloc(setting%DsymLHS,ndmat)
  call LS_MPIBCAST(setting%DsymLHS,ndmat,Master)
ENDIF
IF(SETTING%RHSdmat)THEN
  ndmat = SETTING%nDmatRHS
  IF(SLAVE)THEN
     nullify(setting%DmatRHS)
     allocate(setting%DmatRHS(SETTING%nDmatRHS))
  ENDIF
  do I=1,SETTING%nDmatRHS
     IF(SLAVE)NULLIFY(setting%DmatRHS(I)%p)
     call mat_mpicopy(setting%DmatRHS(I)%p,Slave,Master)
  enddo
  IF(SLAVE)call mem_alloc(setting%DsymRHS,ndmat)
  call LS_MPIBCAST(setting%DsymRHS,ndmat,Master) 
ENDIF

IF(SLAVE)THEN
   SETTING%LHSdalloc = .FALSE.
   SETTING%RHSdalloc = .FALSE.
ENDIF

call LS_MPIBCAST(SETTING%LHSdfull,Master) 
IF(SETTING%LHSdfull)THEN
   ndmat = SETTING%nDmatLHS
   IF(.NOT.SLAVE)THEN
      dim1 = SIZE(SETTING%DfullLHS, 1)  
      dim2 = SIZE(SETTING%DfullLHS, 2)  
      dim3 = SIZE(SETTING%DfullLHS, 3)  
      call LS_MPIBCAST(dim1,Master) 
      call LS_MPIBCAST(dim2,Master) 
      call LS_MPIBCAST(dim3,Master) 
   ENDIF
   IF(SLAVE)THEN
      nullify(setting%DfullLHS)
      allocate(SETTING%DfullLHS(dim1,dim2,dim3))
      SETTING%LHSdalloc = .TRUE.
   ENDIF
   call LS_MPIBCAST(SETTING%DfullLHS,dim1,dim2,dim3,Master) 
   IF(.NOT.ASSOCIATED(setting%DsymLHS))THEN
      IF(SLAVE)call mem_alloc(setting%DsymLHS,ndmat)
      call LS_MPIBCAST(setting%DsymLHS,ndmat,Master)
   ENDIF
ENDIF 

call LS_MPIBCAST(SETTING%RHSdfull,Master) 
IF(SETTING%RHSdfull)THEN
   ndmat = SETTING%nDmatRHS
   IF(.NOT.SLAVE)THEN
      dim1 = SIZE(SETTING%DfullRHS, 1)  
      dim2 = SIZE(SETTING%DfullRHS, 2)  
      dim3 = SIZE(SETTING%DfullRHS, 3)  
      call LS_MPIBCAST(dim1,Master) 
      call LS_MPIBCAST(dim2,Master) 
      call LS_MPIBCAST(dim3,Master) 
   ENDIF
   IF(SLAVE)THEN
      nullify(setting%DfullRHS)
      allocate(SETTING%DfullRHS(dim1,dim2,dim3))
      SETTING%RHSdalloc = .TRUE.
   ENDIF
   call LS_MPIBCAST(SETTING%DfullRHS,dim1,dim2,dim3,Master) 
   IF(.NOT.ASSOCIATED(setting%DsymRHS))THEN
      IF(SLAVE)call mem_alloc(setting%DsymRHS,ndmat)
      call LS_MPIBCAST(setting%DsymRHS,ndmat,Master)
   ENDIF
ENDIF 

IF(SLAVE)THEN
   CALL io_init(SETTING%IO)
ENDIF

call LS_MPIBCAST(SETTING%LHSdmatAOindex1,Master)
call LS_MPIBCAST(SETTING%LHSdmatAOindex2,Master)
call LS_MPIBCAST(SETTING%RHSdmatAOindex1,Master)
call LS_MPIBCAST(SETTING%RHSdmatAOindex2,Master)

call mpicopy_scheme(setting%SCHEME,Slave,Master)

call mpicopy_output(setting%output,Slave,Master)
#endif

END SUBROUTINE mpicopy_setting

#ifdef VAR_LSMPI
SUBROUTINE mpicopy_output(output,Slave,Master)
  use lsmpi_mod
implicit none
type(integralOutput) :: output
logical :: slave
integer :: master
call LS_MPIBCAST(output%ndim,5,Master)
call LS_MPIBCAST(output%doGrad,Master)
call LS_MPIBCAST(output%USEBUFMM,Master)
call LS_MPIBCAST(output%MMBUFLEN,Master)
call LS_MPIBCAST(output%MAXBUFI,Master)
call LS_MPIBCAST(output%MAXBUFR,Master)
call LS_MPIBCAST(output%MAXBUFN,Master)
call LS_MPIBCAST(output%IBUFI,Master)
call LS_MPIBCAST(output%IBUFN,Master)
call LS_MPIBCAST(output%LUITNM,Master)
call LS_MPIBCAST(output%LUITNMR,Master)
call LS_MPIBCAST(output%docontrule,Master)
call LS_MPIBCAST(output%ncontrule,Master)
END SUBROUTINE mpicopy_output

SUBROUTINE mpicopy_scheme(scheme,slave,master)
  use lsmpi_mod
implicit none
type(lsintscheme) :: scheme
logical :: slave
integer :: master

!PARAMETERS FROM **INTEGRALS   DECLERATION
call LS_MPIBCAST(scheme%CFG_LSDALTON,Master)
call LS_MPIBCAST(scheme%DOPASS,Master)
call LS_MPIBCAST(scheme%DENSFIT,Master)
call LS_MPIBCAST(scheme%AOPRINT,Master)
call LS_MPIBCAST(scheme%INTPRINT,Master)
call LS_MPIBCAST(scheme%JENGINE,Master)
call LS_MPIBCAST(scheme%FTUVmaxprim,Master)
call LS_MPIBCAST(scheme%maxpasses,Master)
call LS_MPIBCAST(scheme%FMM,Master)
call LS_MPIBCAST(scheme%LINK,Master)
call LS_MPIBCAST(scheme%DALINK,Master)
call LS_MPIBCAST(scheme%DEBUGOVERLAP,Master)
call LS_MPIBCAST(scheme%DEBUG4CENTER,Master)
call LS_MPIBCAST(scheme%DEBUG4CENTER_ERI,Master)
call LS_MPIBCAST(scheme%DEBUGCCFRAGMENT,Master)
call LS_MPIBCAST(scheme%DEBUGKINETIC,Master)
call LS_MPIBCAST(scheme%DEBUGNUCREP,Master)
call LS_MPIBCAST(scheme%DO4CENTERERI,Master)
call LS_MPIBCAST(scheme%OVERLAP_DF_J,Master)
call LS_MPIBCAST(scheme%TIMINGS,Master)
call LS_MPIBCAST(scheme%nonSphericalETUV,Master)
call LS_MPIBCAST(scheme%ETUVinsideLOOP,Master)
call LS_MPIBCAST(scheme%HIGH_RJ000_ACCURACY,Master)
call LS_MPIBCAST(scheme%OrderAngPrim,Master)
!*FMM PARAMETERS
call LS_MPIBCAST(scheme%MM_LMAX,Master)
call LS_MPIBCAST(scheme%MM_TLMAX,Master)
call LS_MPIBCAST(scheme%MM_SCREEN,Master)
call LS_MPIBCAST(scheme%NO_MMFILES,Master)
call LS_MPIBCAST(scheme%MM_NO_ONE,Master)
call LS_MPIBCAST(scheme%CREATED_MMFILES,Master)
call LS_MPIBCAST(scheme%USEBUFMM,Master)
call LS_MPIBCAST(scheme%MMunique_ID1,Master)
!*BASIS PARAMETERS
call LS_MPIBCAST(scheme%AUXBASIS,Master)
call LS_MPIBCAST(scheme%NOFAMILY,Master)
call LS_MPIBCAST(scheme%DoCartesian,Master)
call LS_MPIBCAST(scheme%DoSpherical,Master)
call LS_MPIBCAST(scheme%UNCONT,Master) !FORCE UNCONTRACTED BASIS
call LS_MPIBCAST(scheme%NOSEGMENT,Master) !DISABLE SEGMENTS 
!* JOB REQUESTS
call LS_MPIBCAST(scheme%DO3CENTEROVL,Master)
call LS_MPIBCAST(scheme%DO2CENTERERI,Master)
call LS_MPIBCAST(scheme%CARMOM,Master)
call LS_MPIBCAST(scheme%SPHMOM,Master)
call LS_MPIBCAST(scheme%MIXEDOVERLAP,Master)
call LS_MPIBCAST(scheme%TEST_NDMAT_COULOMB,Master)
call LS_MPIBCAST(scheme%TEST_NDMAT_EXCHANGE,Master)
!*CAUCHY-SCHWARZ INTEGRAL PARAMETERS
!THE ONE THRESHOLD TO RULE THEM ALL
call LS_MPIBCAST(scheme%THRESHOLD,Master)
!THESE THRESHOLDS TELL HOW THEY SHOULD BE SET COMPARED TO THE ONE THRESHOLD
call LS_MPIBCAST(scheme%CS_THRESHOLD,Master)
call LS_MPIBCAST(scheme%OE_THRESHOLD,Master)
call LS_MPIBCAST(scheme%PS_THRESHOLD,Master)
call LS_MPIBCAST(scheme%OD_THRESHOLD,Master)
call LS_MPIBCAST(scheme%Integralthreshold,Master)
!OTHER CAUCHY-SCHWARZ INTEGRAL PARAMETERS
call LS_MPIBCAST(scheme%CS_SCREEN,Master)
call LS_MPIBCAST(scheme%OE_SCREEN,Master)
call LS_MPIBCAST(scheme%GAB_ON_FILE,Master)
call LS_MPIBCAST(scheme%WRITE_GAB_TO_FILE,Master)
!*PRIMITIVE INTEGRAL PARAMETERS
call LS_MPIBCAST(scheme%PS_SCREEN,Master)
call LS_MPIBCAST(scheme%PRIM_GAB_ON_FILE,Master)
call LS_MPIBCAST(scheme%WRITE_PRIM_GAB_TO_FILE,Master)
call LS_MPIBCAST(scheme%PS_DEBUG,Master)
!Screen OD-batches by AO-batch extent
call LS_MPIBCAST(scheme%OD_SCREEN,Master)
!Fragment molecule into to distinct parts, and construct matrices block by block
call LS_MPIBCAST(scheme%FRAGMENT,Master)
!Approximate number of atoms per fragment
call LS_MPIBCAST(scheme%numAtomsPerFragment,Master)
call LS_MPIBCAST(scheme%DECPACKED,Master)
!FMM
call LS_MPIBCAST(scheme%PUREFMM,Master)
call LS_MPIBCAST(scheme%LU_LUINTM,Master)
call LS_MPIBCAST(scheme%LU_LUINTR,Master)
!Coulomb attenuated method CAM parameters
call LS_MPIBCAST(scheme%CAM,Master)
call LS_MPIBCAST(scheme%CAMalpha,Master)
call LS_MPIBCAST(scheme%CAMbeta,Master)
call LS_MPIBCAST(scheme%CAMmu,Master)
call LS_MPIBCAST(scheme%exchangeFactor,Master)
!DFT PARAMETERS
call LS_MPIBCAST(scheme%GRDONE,Master)
call LS_MPIBCAST(scheme%ITERATIONS,Master)
call LS_MPIBCAST(scheme%RADINT,Master)
call LS_MPIBCAST(scheme%ANGMIN,Master)
call LS_MPIBCAST(scheme%ANGINT,Master)
call LS_MPIBCAST(scheme%HRDNES,Master)
call LS_MPIBCAST(scheme%DFTELS,Master)
call LS_MPIBCAST(scheme%DFTHR0,Master)
call LS_MPIBCAST(scheme%DFTHRI,Master)
call LS_MPIBCAST(scheme%DFTHRL,Master)
call LS_MPIBCAST(scheme%RHOTHR,Master)
call LS_MPIBCAST(scheme%NOPRUN,Master)
call LS_MPIBCAST(scheme%DFTASC,Master)
call LS_MPIBCAST(scheme%DFTPOT,Master)
call LS_MPIBCAST(scheme%DODISP,Master)
call LS_MPIBCAST(scheme%DFTIPT,Master)
call LS_MPIBCAST(scheme%DFTBR1,Master)
call LS_MPIBCAST(scheme%DFTBR2,Master)
call LS_MPIBCAST(scheme%DFTADD,Master)
call LS_MPIBCAST(scheme%DISPDONE,Master)
call LS_MPIBCAST(scheme%TURBO,Master)
  
END SUBROUTINE mpicopy_scheme
#endif

SUBROUTINE II_setMolecules_2_1(setting,mol1,ao1,ao2,mol2,ao3)
implicit none
TYPE(LSSETTING),intent(inout)        :: setting
TYPE(MOLECULEINFO),intent(in),target :: mol1,mol2
Integer, intent(in)                  :: ao1,ao2,ao3
!
Integer :: iao

setting%molecule(ao1)%p => mol1
setting%molecule(ao2)%p => mol1
setting%molecule(ao3)%p => mol2

setting%sameMol = .false.
DO iao=1,4
  setting%sameMol(iao,iao) = .true.
ENDDO
setting%sameMol(ao1,ao2) = .true.
setting%sameMol(ao2,ao1) = .true.

END SUBROUTINE II_setMolecules_2_1

SUBROUTINE II_setMolecules_2(setting,mol1,ao1,ao2)
implicit none
TYPE(LSSETTING),intent(inout)        :: setting
TYPE(MOLECULEINFO),intent(in),target :: mol1
Integer, intent(in)                  :: ao1,ao2
!
Integer :: iao

setting%molecule(ao1)%p => mol1
setting%molecule(ao2)%p => mol1

setting%sameMol = .false.
DO iao=1,4
  setting%sameMol(iao,iao) = .true.
ENDDO
setting%sameMol(ao1,ao2) = .true.
setting%sameMol(ao2,ao1) = .true.

END SUBROUTINE II_setMolecules_2

SUBROUTINE II_setMolecules_1_1_1(setting,mol1,ao1,mol2,ao2,mol3,ao3)
implicit none
TYPE(LSSETTING),intent(inout)        :: setting
TYPE(MOLECULEINFO),intent(in),target :: mol1,mol2,mol3
Integer, intent(in)                  :: ao1,ao2,ao3
!
Integer :: iao

setting%molecule(ao1)%p => mol1
setting%molecule(ao2)%p => mol2
setting%molecule(ao3)%p => mol3

setting%sameMol = .false.
DO iao=1,4
  setting%sameMol(iao,iao) = .true.
ENDDO

END SUBROUTINE II_setMolecules_1_1_1

SUBROUTINE II_setMolecules_1_1_1_1(setting,mol1,ao1,mol2,ao2,mol3,ao3,mol4,ao4)
implicit none
TYPE(LSSETTING),intent(inout)        :: setting
TYPE(MOLECULEINFO),intent(in),target :: mol1,mol2,mol3,mol4
Integer, intent(in)                  :: ao1,ao2,ao3,ao4
!
Integer :: iao

setting%molecule(ao1)%p => mol1
setting%molecule(ao2)%p => mol2
setting%molecule(ao3)%p => mol3
setting%molecule(ao4)%p => mol4

setting%sameMol = .false.
DO iao=1,4
  setting%sameMol(iao,iao) = .true.
ENDDO

END SUBROUTINE II_setMolecules_1_1_1_1

SUBROUTINE II_setMolecules_4(setting,mol,ao1,ao2,ao3,ao4)
implicit none
TYPE(LSSETTING),intent(inout)        :: setting
TYPE(MOLECULEINFO),intent(in),target :: mol
Integer, intent(in)                  :: ao1,ao2,ao3,ao4

setting%molecule(ao1)%p => mol
setting%molecule(ao2)%p => mol
setting%molecule(ao3)%p => mol
setting%molecule(ao4)%p => mol

setting%sameMol = .true.

END SUBROUTINE II_setMolecules_4

!> \brief set default values for the lssetting structure
!> \author T. Kjaergaard
!> \date 2010
!> \param setting the setting structure
!> \param input the daltoninput structure
SUBROUTINE II_set_default_setting(SETTING,INPUT)
implicit none
TYPE(LSSETTING)   :: SETTING
TYPE(DALTONINPUT) :: INPUT
!
Integer :: iAO
SETTING%DO_DFT = INPUT%DO_DFT
DO iAO=1,SETTING%nAO
  SETTING%MOLECULE(iAO)%p => INPUT%MOLECULE
  SETTING%BASIS(iAO)%p    => INPUT%BASIS
  SETTING%FRAGMENT(iAO)%p => INPUT%MOLECULE
  SETTING%BATCHINDEX(iAO)= 0
  SETTING%BATCHDIM(iAO)= 0
  SETTING%molID(iAO)= 0
ENDDO
SETTING%sameMOL  = .TRUE.
SETTING%sameBAS  = .TRUE.
SETTING%sameFRAG = .TRUE.

!Simen ToDo Should be removed from input
SETTING%node = INPUT%node
SETTING%numNodes = INPUT%numNodes
SETTING%numFragments = INPUT%numFragments
CALL II_setIntegralSchemeFromInput(INPUT%DALTON,SETTING%SCHEME)

END SUBROUTINE II_set_default_setting

!> \brief initialise the lssetting structure
!> \author T. Kjaergaard
!> \date 2010
!> \param setting the setting structure
SUBROUTINE II_init_setting(SETTING)
implicit none
TYPE(LSSETTING)  :: SETTING
INTEGER          :: nAO

nAO = SETTING%nAO

NULLIFY(SETTING%MOLECULE)
NULLIFY(SETTING%BASIS)
NULLIFY(SETTING%FRAGMENT)
NULLIFY(SETTING%BATCHINDEX)
NULLIFY(SETTING%BATCHDIM)
NULLIFY(SETTING%molID)
ALLOCATE(SETTING%MOLECULE(nAO))
ALLOCATE(SETTING%BASIS(nAO))
ALLOCATE(SETTING%FRAGMENT(nAO))
ALLOCATE(SETTING%BATCHINDEX(nAO))
ALLOCATE(SETTING%BATCHDIM(nAO))
ALLOCATE(SETTING%molID(nAO))
NULLIFY(SETTING%molBuild)
NULLIFY(SETTING%basBuild)
NULLIFY(SETTING%fragBuild)
ALLOCATE(SETTING%molBuild(nAO))
ALLOCATE(SETTING%basBuild(nAO))
ALLOCATE(SETTING%fragBuild(nAO))
SETTING%molBuild = .FALSE.
SETTING%basBuild = .FALSE.
SETTING%fragBuild = .FALSE.
NULLIFY(SETTING%sameMOL)
NULLIFY(SETTING%sameBAS)
NULLIFY(SETTING%sameFRAG)
ALLOCATE(SETTING%sameMOL(nAO,nAO))
ALLOCATE(SETTING%sameBAS(nAO,nAO))
ALLOCATE(SETTING%sameFRAG(nAO,nAO))
NULLIFY(SETTING%DmatLHS)
NULLIFY(SETTING%DmatRHS)
SETTING%nDmatLHS = 1
SETTING%nDmatRHS = 1
SETTING%LHSdmat = .FALSE.
SETTING%RHSdmat = .FALSE.
SETTING%LHSdfull = .FALSE.
SETTING%RHSdfull = .FALSE.
SETTING%LHSdalloc = .FALSE.
SETTING%RHSdalloc = .FALSE.
NULLIFY(SETTING%DsymLHS)
NULLIFY(SETTING%DsymRHS)
CALL io_init(SETTING%IO)

END SUBROUTINE II_init_setting

!> \brief free the lssetting structure
!> \author T. Kjaergaard
!> \date 2010
!> \param setting the setting structure
SUBROUTINE II_free_setting(SETTING)
!use molecule_module
implicit none
TYPE(LSSETTING)  :: SETTING
INTEGER          :: nAO,iAO

nAO = SETTING%nAO

DO iAO=1,nAO
  IF (SETTING%molBuild(iAO)) call free_Moleculeinfo(SETTING%MOLECULE(iAO)%p)
  IF (SETTING%basBuild(iAO)) THEN
!Simen ToDo make a freeBasis
    call free_basissetinfo(SETTING%BASIS(iAO)%p%REGULAR)
    IF (SETTING%SCHEME%AUXBASIS) call free_basissetinfo(SETTING%BASIS(iAO)%p%AUXILIARY)
  ENDIF
  IF (SETTING%fragBuild(iAO)) call free_Moleculeinfo(SETTING%FRAGMENT(iAO)%p)
ENDDO
DEALLOCATE(SETTING%MOLECULE)
DEALLOCATE(SETTING%BASIS)
DEALLOCATE(SETTING%FRAGMENT)
DEALLOCATE(SETTING%BATCHINDEX)
DEALLOCATE(SETTING%BATCHDIM)
DEALLOCATE(SETTING%molID)
NULLIFY(SETTING%MOLECULE)
NULLIFY(SETTING%BASIS)
NULLIFY(SETTING%FRAGMENT)
NULLIFY(SETTING%Batchindex)
NULLIFY(SETTING%Batchdim)
NULLIFY(SETTING%molID)
DEALLOCATE(SETTING%molBuild)
DEALLOCATE(SETTING%basBuild)
DEALLOCATE(SETTING%fragBuild)
NULLIFY(SETTING%molBuild)
NULLIFY(SETTING%basBuild)
NULLIFY(SETTING%fragBuild)
DEALLOCATE(SETTING%sameMOL)
DEALLOCATE(SETTING%sameBAS)
DEALLOCATE(SETTING%sameFRAG)
NULLIFY(SETTING%sameMOL)
NULLIFY(SETTING%sameBAS)
NULLIFY(SETTING%sameFRAG)
CALL io_free(SETTING%IO)

END SUBROUTINE II_free_setting

!> \brief set up the lsint scheme from the integralconfig 
!> \author T. Kjaergaard
!> \date 2010
!> \param dalton_inp the integralconfig structure
!> \param scheme the lsintscheme structure
SUBROUTINE II_setIntegralSchemeFromInput(dalton_inp,scheme)
implicit none
TYPE(integralconfig), INTENT(IN) :: dalton_inp
TYPE(LSINTSCHEME),INTENT(OUT) :: scheme

scheme%CFG_LSDALTON          = dalton_inp%CFG_LSDALTON
scheme%DOPASS                = dalton_inp%DOPASS
scheme%DENSFIT               = dalton_inp%DENSFIT
scheme%AOPRINT               = dalton_inp%AOPRINT
scheme%INTPRINT              = dalton_inp%INTPRINT
scheme%JENGINE               = dalton_inp%JENGINE
scheme%FTUVmaxprim           = dalton_inp%FTUVmaxprim
scheme%maxpasses             = dalton_inp%maxpasses
scheme%FMM                   = dalton_inp%FMM
scheme%LINK                  = dalton_inp%LINK
scheme%DALINK                = dalton_inp%DALINK
scheme%DEBUGOVERLAP          = dalton_inp%DEBUGOVERLAP
scheme%DEBUG4CENTER          = dalton_inp%DEBUG4CENTER
scheme%DEBUG4CENTER_ERI      = dalton_inp%DEBUG4CENTER_ERI
scheme%DEBUGCCFRAGMENT       = dalton_inp%DEBUGCCFRAGMENT
scheme%DEBUGKINETIC          = dalton_inp%DEBUGKINETIC
scheme%DEBUGNUCREP           = dalton_inp%DEBUGNUCREP
scheme%DO4CENTERERI          = dalton_inp%DO4CENTERERI
scheme%OVERLAP_DF_J          = dalton_inp%OVERLAP_DF_J
scheme%TIMINGS               = dalton_inp%TIMINGS
scheme%nonSphericalETUV      = dalton_inp%nonSphericalETUV
scheme%ETUVinsideLOOP        = dalton_inp%ETUVinsideLOOP
scheme%HIGH_RJ000_ACCURACY   = dalton_inp%HIGH_RJ000_ACCURACY
scheme%OrderAngPrim          = dalton_inp%OrderAngPrim
scheme%MM_LMAX               = dalton_inp%MM_LMAX
scheme%MM_TLMAX              = dalton_inp%MM_TLMAX
scheme%MM_SCREEN             = dalton_inp%MM_SCREEN
scheme%NO_MMFILES            = dalton_inp%NO_MMFILES
scheme%MM_NO_ONE             = dalton_inp%MM_NO_ONE
scheme%CREATED_MMFILES       = dalton_inp%CREATED_MMFILES
scheme%USEBUFMM              = dalton_inp%USEBUFMM
scheme%MMunique_ID1          = dalton_inp%MMunique_ID1
scheme%AUXBASIS              = dalton_inp%AUXBASIS
scheme%NOFAMILY              = dalton_inp%NOFAMILY
scheme%DoCartesian           = dalton_inp%DoCartesian
scheme%DoSpherical           = dalton_inp%DoSpherical
scheme%UNCONT                = dalton_inp%UNCONT
scheme%NOSEGMENT             = dalton_inp%NOSEGMENT
scheme%DO3CENTEROVL          = dalton_inp%DO3CENTEROVL
scheme%DO2CENTERERI          = dalton_inp%DO2CENTERERI
scheme%CARMOM                = dalton_inp%CARMOM
scheme%SPHMOM                = dalton_inp%SPHMOM
scheme%MIXEDOVERLAP          = dalton_inp%MIXEDOVERLAP
scheme%TEST_NDMAT_COULOMB    = dalton_inp%TEST_NDMAT_COULOMB
scheme%TEST_NDMAT_EXCHANGE   = dalton_inp%TEST_NDMAT_EXCHANGE
scheme%THRESHOLD             = dalton_inp%THRESHOLD
scheme%CS_THRESHOLD          = dalton_inp%CS_THRESHOLD
scheme%OE_THRESHOLD          = dalton_inp%OE_THRESHOLD
scheme%PS_THRESHOLD          = dalton_inp%PS_THRESHOLD
scheme%OD_THRESHOLD          = dalton_inp%OD_THRESHOLD
scheme%Integralthreshold     = dalton_inp%Integralthreshold
scheme%CS_SCREEN             = dalton_inp%CS_SCREEN
scheme%OE_SCREEN             = dalton_inp%OE_SCREEN
scheme%GAB_ON_FILE           = dalton_inp%GAB_ON_FILE
scheme%WRITE_GAB_TO_FILE     = dalton_inp%WRITE_GAB_TO_FILE
scheme%PS_SCREEN             = dalton_inp%PS_SCREEN
scheme%PRIM_GAB_ON_FILE      = dalton_inp%PRIM_GAB_ON_FILE
scheme%WRITE_PRIM_GAB_TO_FILE= dalton_inp%WRITE_PRIM_GAB_TO_FILE
scheme%PS_DEBUG              = dalton_inp%PS_DEBUG
scheme%OD_SCREEN             = dalton_inp%OD_SCREEN
scheme%FRAGMENT              = dalton_inp%FRAGMENT
scheme%numAtomsPerFragment   = dalton_inp%numAtomsPerFragment
scheme%PUREFMM               = dalton_inp%PUREFMM
scheme%LU_LUINTM             = dalton_inp%LU_LUINTM
scheme%LU_LUINTR             = dalton_inp%LU_LUINTR
scheme%CAM                   = dalton_inp%CAM
scheme%CAMalpha              = dalton_inp%CAMalpha
scheme%CAMbeta               = dalton_inp%CAMbeta
scheme%CAMmu                 = dalton_inp%CAMmu
scheme%exchangeFactor        = dalton_inp%exchangeFactor

scheme%GRDONE  = dalton_inp%GRDONE
scheme%ITERATIONS   = dalton_inp%ITERATIONS
scheme%RADINT  = dalton_inp%RADINT
scheme%ANGMIN  = dalton_inp%ANGMIN
scheme%ANGINT  = dalton_inp%ANGINT
scheme%HRDNES   = dalton_inp%HRDNES
scheme%DFTELS   = dalton_inp%DFTELS
scheme%DFTHR0  = dalton_inp%DFTHR0
scheme%DFTHRI  = dalton_inp%DFTHRI
scheme%DFTHRL  = dalton_inp%DFTHRL
scheme%RHOTHR  = dalton_inp%RHOTHR
scheme%NOPRUN  = dalton_inp%NOPRUN
scheme%DFTASC  = dalton_inp%DFTASC
scheme%DFTPOT  = dalton_inp%DFTPOT 
scheme%DODISP  = dalton_inp%DODISP
scheme%DFTIPT  = dalton_inp%DFTIPT
scheme%DFTBR1  = dalton_inp%DFTBR1
scheme%DFTBR2  = dalton_inp%DFTBR2
scheme%DFTADD  = dalton_inp%DFTADD
scheme%DISPDONE = dalton_inp%DISPDONE
scheme%TURBO    = dalton_inp%TURBO
scheme%DECPACKED    = .FALSE.
scheme%INCREMENTAL  = .FALSE.

END SUBROUTINE II_setIntegralSchemeFromInput

!> \brief print the lsint scheme
!> \author T. Kjaergaard
!> \date 2010
!> \param scheme the lsintscheme structure to be printet
!> \param iunit the logical unit number of output 
SUBROUTINE II_printScheme(scheme,IUNIT)
implicit none
TYPE(LSINTSCHEME),INTENT(IN) :: scheme
INTEGER,INTENT(IN)           :: IUNIT

WRITE(IUNIT,'(3X,A22,L7)') 'CFG_LSDALTON          ', scheme%CFG_LSDALTON          
WRITE(IUNIT,'(3X,A22,L7)') 'DOPASS                ', scheme%DOPASS                
WRITE(IUNIT,'(3X,A22,L7)') 'DENSFIT               ', scheme%DENSFIT               
WRITE(IUNIT,'(3X,A22,I7)') 'AOPRINT               ', scheme%AOPRINT              
WRITE(IUNIT,'(3X,A22,I7)') 'INTPRINT              ', scheme%INTPRINT              
WRITE(IUNIT,'(3X,A22,L7)') 'JENGINE               ', scheme%JENGINE               
WRITE(IUNIT,'(3X,A22,I7)') 'FTUVmaxprim           ', scheme%FTUVmaxprim           
WRITE(IUNIT,'(3X,A22,I7)') 'maxpasses             ', scheme%maxpasses             
WRITE(IUNIT,'(3X,A22,L7)') 'FMM                   ', scheme%FMM                   
WRITE(IUNIT,'(3X,A22,L7)') 'LINK                  ', scheme%LINK                  
WRITE(IUNIT,'(3X,A22,L7)') 'DALINK                ', scheme%DALINK                
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUGOVERLAP          ', scheme%DEBUGOVERLAP
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUG4CENTER          ', scheme%DEBUG4CENTER
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUG4CENTER_ERI      ', scheme%DEBUG4CENTER_ERI
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUGCCFRAGMENT       ', scheme%DEBUGCCFRAGMENT
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUGKINETIC          ', scheme%DEBUGKINETIC          
WRITE(IUNIT,'(3X,A22,L7)') 'DEBUGNUCREP           ', scheme%DEBUGNUCREP           
WRITE(IUNIT,'(3X,A22,L7)') 'DO4CENTERERI          ', scheme%DO4CENTERERI          
WRITE(IUNIT,'(3X,A22,L7)') 'OVERLAP_DF_J          ', scheme%OVERLAP_DF_J          
WRITE(IUNIT,'(3X,A22,L7)') 'TIMINGS               ', scheme%TIMINGS               
WRITE(IUNIT,'(3X,A22,L7)') 'nonSphericalETUV      ', scheme%nonSphericalETUV      
WRITE(IUNIT,'(3X,A22,L7)') 'ETUVinsideLOOP        ', scheme%ETUVinsideLOOP        
WRITE(IUNIT,'(3X,A22,L7)') 'HIGH_RJ000_ACCURACY   ', scheme%HIGH_RJ000_ACCURACY   
WRITE(IUNIT,'(3X,A22,L7)') 'OrderAngPrim          ', scheme%OrderAngPrim          
WRITE(IUNIT,'(3X,A22,I7)') 'MM_LMAX               ', scheme%MM_LMAX               
WRITE(IUNIT,'(3X,A22,I7)') 'MM_TLMAX              ', scheme%MM_TLMAX              
WRITE(IUNIT,'(3X,A22,G14.2)') 'MM_SCREEN             ', scheme%MM_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'NO_MMFILES            ', scheme%NO_MMFILES            
WRITE(IUNIT,'(3X,A22,L7)') 'MM_NO_ONE             ', scheme%MM_NO_ONE             
WRITE(IUNIT,'(3X,A22,L7)') 'CREATED_MMFILES       ', scheme%CREATED_MMFILES       
WRITE(IUNIT,'(3X,A22,L7)') 'USEBUFMM              ', scheme%USEBUFMM              
WRITE(IUNIT,'(3X,A22,I7)') 'MMunique_ID1          ', scheme%MMunique_ID1          
WRITE(IUNIT,'(3X,A22,L7)') 'AUXBASIS              ', scheme%AUXBASIS              
WRITE(IUNIT,'(3X,A22,L7)') 'NOFAMILY              ', scheme%NOFAMILY              
WRITE(IUNIT,'(3X,A22,L7)') 'DoCartesian           ', scheme%DoCartesian           
WRITE(IUNIT,'(3X,A22,L7)') 'DoSpherical           ', scheme%DoSpherical           
WRITE(IUNIT,'(3X,A22,L7)') 'UNCONT                ', scheme%UNCONT                
WRITE(IUNIT,'(3X,A22,L7)') 'NOSEGMENT             ', scheme%NOSEGMENT             
WRITE(IUNIT,'(3X,A22,L7)') 'DO3CENTEROVL          ', scheme%DO3CENTEROVL          
WRITE(IUNIT,'(3X,A22,L7)') 'DO2CENTERERI          ', scheme%DO2CENTERERI          
WRITE(IUNIT,'(3X,A22,I7)') 'CARMOM                ', scheme%CARMOM                
WRITE(IUNIT,'(3X,A22,I7)') 'SPHMOM                ', scheme%SPHMOM                
WRITE(IUNIT,'(3X,A22,L7)') 'MIXEDOVERLAP          ', scheme%MIXEDOVERLAP          
WRITE(IUNIT,'(3X,A22,L7)') 'TEST_NDMAT_COULOMB    ', scheme%TEST_NDMAT_COULOMB    
WRITE(IUNIT,'(3X,A22,L7)') 'TEST_NDMAT_EXCHANGE   ', scheme%TEST_NDMAT_EXCHANGE    
WRITE(IUNIT,'(3X,A22,G14.2)') 'THRESHOLD             ', scheme%THRESHOLD             
WRITE(IUNIT,'(3X,A22,G14.2)') 'CS_THRESHOLD          ', scheme%CS_THRESHOLD          
WRITE(IUNIT,'(3X,A22,G14.2)') 'OE_THRESHOLD          ', scheme%OE_THRESHOLD          
WRITE(IUNIT,'(3X,A22,G14.2)') 'PS_THRESHOLD          ', scheme%PS_THRESHOLD          
WRITE(IUNIT,'(3X,A22,G14.2)') 'OD_THRESHOLD          ', scheme%OD_THRESHOLD          
WRITE(IUNIT,'(3X,A22,G14.2)') 'Integralthreshold     ', scheme%Integralthreshold     
WRITE(IUNIT,'(3X,A22,L7)') 'CS_SCREEN             ', scheme%CS_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'OE_SCREEN             ', scheme%OE_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'GAB_ON_FILE           ', scheme%GAB_ON_FILE           
WRITE(IUNIT,'(3X,A22,L7)') 'WRITE_GAB_TO_FILE     ', scheme%WRITE_GAB_TO_FILE     
WRITE(IUNIT,'(3X,A22,L7)') 'PS_SCREEN             ', scheme%PS_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'PRIM_GAB_ON_FILE      ', scheme%PRIM_GAB_ON_FILE      
WRITE(IUNIT,'(3X,A22,L7)') 'WRITE_PRIM_GAB_TO_FILE', scheme%WRITE_PRIM_GAB_TO_FILE
WRITE(IUNIT,'(3X,A22,L7)') 'PS_DEBUG              ', scheme%PS_DEBUG              
WRITE(IUNIT,'(3X,A22,L7)') 'OD_SCREEN             ', scheme%OD_SCREEN             
WRITE(IUNIT,'(3X,A22,L7)') 'FRAGMENT              ', scheme%FRAGMENT              
WRITE(IUNIT,'(3X,A22,I7)') 'numAtomsPerFragment   ', scheme%numAtomsPerFragment   
WRITE(IUNIT,'(3X,A22,L7)') 'PUREFMM               ', scheme%PUREFMM               
WRITE(IUNIT,'(3X,A22,I7)') 'LU_LUINTM             ', scheme%LU_LUINTM             
WRITE(IUNIT,'(3X,A22,I7)') 'LU_LUINTR             ', scheme%LU_LUINTR             
WRITE(IUNIT,'(3X,A22,L7)') 'CAM                   ', scheme%CAM                   
WRITE(IUNIT,'(3X,A22,G14.2)') 'CAMalpha              ', scheme%CAMalpha              
WRITE(IUNIT,'(3X,A22,G14.2)') 'CAMbeta               ', scheme%CAMbeta               
WRITE(IUNIT,'(3X,A22,G14.2)') 'CAMmu                 ', scheme%CAMmu                 
WRITE(IUNIT,'(3X,A22,G14.2)') 'exchangeFactor        ', scheme%exchangeFactor        

WRITE(IUNIT,'(3X,A22,L7)')'GRDONE                ', scheme%GRDONE 
WRITE(IUNIT,'(3X,A22,I7)')'ITERATIONS            ', scheme%ITERATIONS
WRITE(IUNIT,'(3X,A22,I7)')'ANGMIN                ', scheme%ANGMIN 
WRITE(IUNIT,'(3X,A22,I7)')'HRDNES                ', scheme%HRDNES 
WRITE(IUNIT,'(3X,A22,I7)')'TURBO                 ', scheme%TURBO   
WRITE(IUNIT,'(3X,A22,G14.2)') 'RADINT                ', scheme%RADINT
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTELS                ', scheme%DFTELS
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTHR0                ', scheme%DFTHR0
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTHRI                ', scheme%DFTHRI
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTHRL                ', scheme%DFTHRL
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTHRL                ', scheme%DFTHRL
WRITE(IUNIT,'(3X,A22,G14.2)') 'RHOTHR                ', scheme%RHOTHR
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTIPT                ', scheme%DFTIPT
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTBR1                ', scheme%DFTBR1
WRITE(IUNIT,'(3X,A22,G14.2)') 'DFTBR2                ', scheme%DFTBR2
WRITE(IUNIT,'(3X,A22,L7)')'NOPRUN                ', scheme%NOPRUN 
WRITE(IUNIT,'(3X,A22,L7)')'DFTASC                ', scheme%DFTASC
WRITE(IUNIT,'(3X,A22,L7)')'DFTPOT                ', scheme%DFTPOT
WRITE(IUNIT,'(3X,A22,L7)')'DISPER                ', scheme%DODISP
WRITE(IUNIT,'(3X,A22,L7)')'DISPDONE              ', scheme%DISPDONE
WRITE(IUNIT,'(3X,A22,L7)')'DFTADD                ', scheme%DFTADD
WRITE(IUNIT,'(3X,A22,L7)')'DECPACKED             ', scheme%DECPACKED
WRITE(IUNIT,'(3X,A22,L7)')'INCREMENTAL           ', scheme%INCREMENTAL

END SUBROUTINE II_printScheme

END MODULE TYPEDEF

