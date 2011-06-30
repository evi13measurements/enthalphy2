MODULE mm_global_paras_mod

   IMPLICIT NONE
   PUBLIC   

!===============================================================================

   ! We select here the REAL and INTEGER precisions
   ! For consistency with the rest of DALTON we opt for DOUBLE PRECISION (1D0)
#if defined (VAR_INT64)
   INTEGER, PARAMETER :: INTK  = SELECTED_INT_KIND(17),           &
#else
   INTEGER, PARAMETER :: INTK  = SELECTED_INT_KIND(9),           &
#endif
                         REALK = KIND(1D0) 

!===============================================================================
! Derived types !
!================

    ! Structure for counting raw multipole moments by type
   TYPE mm_counters
      INTEGER(INTK) :: tot, nuc, elec, nbast, nauxbas
   END TYPE mm_counters
 
!-------------------------------------------------------------------------------

    ! Complete prescription of MM execution held in "scheme"
   TYPE T_contract_schm 
      INTEGER(INTK) :: ID 
      INTEGER(INTK) :: NN_ID    ! allow possibility for different NN contractor
      INTEGER(INTK) :: T_buffer
      INTEGER(INTK) :: NN_T_buffer
      INTEGER(INTK) :: sort_para
      INTEGER(INTK) :: LHS_mm_type, RHS_mm_type 
   END TYPE T_contract_schm 

   TYPE W_contract_schm 
      INTEGER(INTK) :: ID 
      INTEGER(INTK) :: W_buffer
      INTEGER(INTK) :: sort_para
   END TYPE W_contract_schm 

    ! T-pair interactions search algorithm controls
   TYPE T_searcher_type
      INTEGER(INTK) :: algorithm, shape
      LOGICAL :: all_square
   END TYPE T_searcher_type

!FIXME: this has probably not been updated correctly!!!
 !NOTE: when changing the definition, change also MPI distribution of the type.
   TYPE scheme_paras
      INTEGER(INTK)         :: algorithm 
      INTEGER(INTK)         :: NLEVEL
      LOGICAL               :: branch_free
      LOGICAL               :: inc_NN
      LOGICAL               :: dynamic_levels
      LOGICAL               :: dynamic_LMAX_on
      LOGICAL               :: contracted
      INTEGER(INTK)         :: raw_LMAX, trans_LMAX
      INTEGER(INTK)         :: LEXTRA 
      INTEGER(INTK)         :: LHS_mm_range, RHS_mm_range 
      LOGICAL               :: pack_LHS, pack_RHS
      LOGICAL               :: LHS_dens, RHS_dens
      INTEGER(INTK)         :: J_maker
      REAL(REALK)           :: cls_radius 
      REAL(REALK)           :: system_size
      REAL(REALK)           :: grain_input 
      REAL(REALK)           :: dens_screen_thr
      TYPE(T_searcher_type) :: T_searcher(5)
      LOGICAL               :: NN_box_pretesting
      TYPE(T_contract_schm) :: T_con 
      TYPE(W_contract_schm) :: W_con 
   END TYPE scheme_paras
 
!-------------------------------------------------------------------------------

    ! raw (from file) multipole moment data structure
   TYPE raw_mm_paras
      REAL(REALK)   :: cntr(3), ext
      INTEGER(INTK) :: batch            ! batch id of unique centres/extents
      INTEGER(INTK) :: id               ! map to raw moments array
      INTEGER(INTK) :: Lmin             ! OD order (s=>0, p=>1 etc.)
      INTEGER(INTK) :: map_up           ! map to boxed multipole moment
      INTEGER(INTK) :: box(3), bra
      REAL(REALK)   :: box_cntr(3)
   END TYPE raw_mm_paras
   TYPE J_index_type
      INTEGER(INTK) :: i_indx, j_indx   ! map to J_matrix elements
   END TYPE J_index_type
   TYPE raw_mm_data 
      TYPE(raw_mm_paras), POINTER :: paras(:)
      REAL(REALK),        POINTER :: dens(:)  
      REAL(REALK),        POINTER :: qlm(:,:)  
      REAL(REALK),        POINTER :: qlm_T(:,:)  
      REAL(REALK),        POINTER :: qlm_W(:,:)  
      TYPE(J_index_type), POINTER :: J_indices(:)
      TYPE(id_list),      POINTER :: batch_map(:) ! maps raw paras in one batch
!ANDREAS
      REAL(REALK),        POINTER :: qlm_der(:,:,:)
      REAL(REALK),        POINTER :: qlm_der_T(:,:,:)
      INTEGER(INTK),      POINTER :: mom2atom(:,:)
      INTEGER(INTK) :: endnuc
      INTEGER(INTK) :: startnuc
!
   END TYPE raw_mm_data 
 
    ! boxed multipole moment data structure
   TYPE box_mm_paras
      INTEGER(INTK) :: box(3)
      REAL(REALK)   :: cntr(3)
      INTEGER(INTK) :: bra, level
      INTEGER(INTK) :: map_up            ! map to moment at next level 
      REAL(REALK)   :: cntr_up(3)        ! centre of parent box
      INTEGER(INTK) :: id                ! map to moment at this level
   END TYPE box_mm_paras
   TYPE box_mm_data
      TYPE(box_mm_paras), POINTER :: LHS_paras(:)
      TYPE(box_mm_paras), POINTER :: RHS_paras(:)
      ! just RHS moments (only RHS is built up in box hierarchy)
      REAL(REALK),        POINTER :: qlm_T(:,:)  
      REAL(REALK),        POINTER :: qlm_W(:,:)  
   END TYPE box_mm_data
 
    ! map structure to link packed parameters with all raw moment members
   TYPE id_node
      INTEGER(INTK) :: id
! we add these two here as J-matrix indices
!FIXME: rethink how we can go about passing J-matrix indices in general!
      INTEGER(INTK) :: i_indx, j_indx
      TYPE(id_node), POINTER :: next
   END TYPE id_node
   TYPE id_list
      INTEGER(INTK) :: occ 
      TYPE(id_node), POINTER :: head
   END TYPE id_list

    ! generalised multipole moment parameter structure
   TYPE gen_mm_paras
      TYPE(raw_mm_paras), POINTER :: raw_paras(:)
      TYPE(box_mm_paras), POINTER :: box_paras(:)
      TYPE(id_list),      POINTER :: box_map(:)   ! maps raw paras in one box
   END TYPE gen_mm_paras
 
!-------------------------------------------------------------------------------

   TYPE LHS_RHS_type 
      INTEGER(INTK) :: LHS, RHS
   END TYPE LHS_RHS_type

   TYPE mm_range
      INTEGER(INTK) :: hi, lo, tot
   END TYPE mm_range
 
!-------------------------------------------------------------------------------

!FIXME: we now describe the pair types fundamental to this code
! Orginally, we had separate descriptions for the T-pair interactions
! and the W-pair translations, but we are now in the process of combining
! both descriptions into a single, unified framework...
! However, the names currently remain as T_xxxxx
! We now also generalise to include "J-pairs" describing the a contribution
! to an element of the J-matrix, requiring 2 more indices.
! Can we break this up/ make things more elegant, since these 2 integers
! are never used when getting energies.

   TYPE T_paras
!FIXME: try other types to reduce memory use here??
      INTEGER(INTK) :: LHS_LMAX, LHS_id
      INTEGER(INTK) :: RHS_LMAX, RHS_id
      INTEGER(INTK) :: weight
!FIXME: try changing to SINGLE PRECISION here??
      REAL(REALK)   :: ratio
   END TYPE T_paras

   TYPE T_pair_single
      TYPE(T_paras) :: paras
      REAL(REALK)   :: r_ab(3)
      INTEGER(INTK) :: LMAX, lm_max
      ! used only for W_pairs when translating to distinguish qlm and Vff modes
      CHARACTER(1)  :: N_or_T
   END TYPE T_pair_single

   TYPE T_pair_list
      TYPE(T_paras), POINTER :: paras(:) 
      REAL(REALK)   :: r_ab(3)
      INTEGER(INTK) :: LMAX, lm_max
      INTEGER(INTK) :: LHS_LMAX, RHS_LMAX
      ! used only for W_pairs when translating to distinguish qlm and Vff modes
      CHARACTER(1)  :: N_or_T
   END TYPE T_pair_list

   TYPE T_pair_batch
      TYPE(T_pair_single), POINTER :: items(:)
      INTEGER(INTK) :: ndim
   END TYPE T_pair_batch
 
    ! index type for OLD and NEW W_translations
   TYPE old_new 
      INTEGER(INTK) :: old, new 
   END TYPE old_new 

!===============================================================================
! Global variables  !
!===================!

    ! Unit number for output file writing
   INTEGER(INTK), SAVE :: LUPRI = 10
    ! Flag to say if we link with main DALTON code, or run "stand alone"
   LOGICAL, SAVE :: MM_STANDALONE = .FALSE.

   ! no KIND is specified here - intentionally! - to achieve interoperability
   ! with MPI.
   ! serial calculation has 1 node with number 0.
   INTEGER :: NNODES=1, MYNUM=0
   ! whether extra statistics should be printed (only initialised here)
   LOGICAL, SAVE :: mm_stats_printed = .FALSE.
 
!===============================================================================
! Global parameters !
!===================!
 
   REAL(REALK), PARAMETER :: zero = 0.0_REALK,                   &
                             one  = 1.0_REALK,                   &
                             two  = 2.0_REALK,                   &
                             half = 0.5_REALK
 
!-------------------------------------------------------------------------------
!
! Parameters for default scheme
!
   REAL(REALK),   PARAMETER :: GRAIN_DF = 2.0_REALK 
   REAL(REALK),   PARAMETER :: RPQMIN_DF = 0.1_REALK 
    ! default parameter for density based screening. (see also DENS_SCREEN_CUT)
   REAL(REALK),   PARAMETER :: DENS_SCREEN_DF = 1e-12_REALK
   INTEGER(INTK), PARAMETER :: ALGORITHM_DF = 5
   ! default NLEVEL is -1 because the default is actually to use dynamic levels
   INTEGER(INTK), PARAMETER :: NLEVEL_DF = -1
   INTEGER(INTK), PARAMETER :: TLMAX_DF = 20 
   INTEGER(INTK), PARAMETER :: LEXTRA_DF = 2
   INTEGER(INTK), PARAMETER :: T_CONTRACTOR_DF = 2
   INTEGER(INTK), PARAMETER :: TMATM_DF = 25
   LOGICAL,       PARAMETER :: INC_NN_DF = .TRUE.
   LOGICAL,       PARAMETER :: USEUMAT_DF = .FALSE.
   LOGICAL,       PARAMETER :: BRFREE_DF = .FALSE.
   LOGICAL,       PARAMETER :: DYNLMAX_DF = .FALSE.
!   LOGICAL,       PARAMETER :: ALLSQR_DF = .FALSE.
!FIXME: changed conservatively since need to make FULL_J work
   LOGICAL,       PARAMETER :: ALLSQR_DF = .TRUE.
   LOGICAL,       PARAMETER :: TRSRCH_DF = .FALSE.
   LOGICAL,       PARAMETER :: GRSRCH_DF = .FALSE.
   LOGICAL,       PARAMETER :: NOBOXP_DF = .FALSE.
   LOGICAL,       PARAMETER :: NOCONT_DF = .FALSE.
   LOGICAL,       PARAMETER :: TRBUFF_DF = .TRUE.
   LOGICAL,       PARAMETER :: DYNLEV_DF = .TRUE.
   LOGICAL,       PARAMETER :: PACK_RHS_DF = .TRUE.
   LOGICAL,       PARAMETER :: PACK_LHS_DF = .TRUE.

!-------------------------------------------------------------------------------

    ! Threshold parameter for density based screening.
    ! Screening will only be performed if threshold is greater than this.
    ! (see also DENS_SCREEN_THR)
   REAL(REALK),   PARAMETER :: DENS_SCREEN_CUT = 1e-20_REALK

    ! Parameters for branch and box structure
   INTEGER(INTK), PARAMETER :: WS_MIN = 2
   REAL(REALK),   PARAMETER :: GRAIN_BUFFER  = 1.00001_REALK 
   LOGICAL,       PARAMETER :: JOIN_BRANCHES = .TRUE. 
 
    ! Parameters for T-list tree sorting
   INTEGER(INTK), PARAMETER :: START_LEN   = 8
   INTEGER(INTK), PARAMETER :: TREE_LENGTH = 50000
   ! MAX_AVG_PER_NODE limits the maximal packing ratio.
   ! there will be never more interactions in the tree than
   ! tree size times MAX_AVG_PER_NODE.
   INTEGER(INTK), PARAMETER :: MAX_AVG_PER_NODE = 15
 
    ! Noise buffer when testing for classical interactions
 !FIXME: have this passed from clsfmm.h ... used for nothing else!!
   REAL(REALK), PARAMETER :: ZERO_DIST_TOL = 1e-14_REALK
 
    ! Parameters for contractions and translations
 !FIXME: clarify the use of these tolerance parameters!!
   REAL(REALK), PARAMETER :: DISTINCT_T_TOL  = 1e-15_REALK
   REAL(REALK), PARAMETER :: EXTENT_TEST_TOL = 1e-15_REALK
   REAL(REALK), PARAMETER :: ZERO_VECT_TOL   = 1e-10_REALK
 
    ! Flag for using UNSCALED solid harmonic formulation
   LOGICAL, PARAMETER :: USE_UNSCALED_HARMONICS = .FALSE. 
 
    ! maxmimum LMAX that is used for constructing FULL T-matrix
   INTEGER(INTK), PARAMETER :: FULL_T_LMAX = 4

!-------------------------------------------------------------------------------
!
! error checks
!
   INTEGER(INTK), PARAMETER :: TLMAX_HIGH = 25 
   INTEGER(INTK), PARAMETER :: LMAX_HIGH = 15 
    ! Parameter to define level depth used in box hierarchy
   INTEGER(INTK), PARAMETER :: TOP_LEVEL = 2
   INTEGER(INTK), PARAMETER :: MAX_LEVEL = 25

!-------------------------------------------------------------------------------
!
! Named parameters for coding use (avoidance of "magic" numbers etc.)
!
    ! Interface file names (contains moments and co-ordinate data)
   CHARACTER(7),  PARAMETER :: INPUT_FILE0 = 'MM_CNT0'
   CHARACTER(7),  PARAMETER :: INPUT_FILE1 = 'MM_DATA' 
   CHARACTER(7),  PARAMETER :: INPUT_FILE2 = 'MM_CNTS' 
   CHARACTER(7),  PARAMETER :: INPUT_FILE3 = 'MM_DENS'
   CHARACTER(7),  PARAMETER :: INPUT_FILE4 = 'MM_DENL'
   CHARACTER(7),  PARAMETER :: INPUT_FILE5 = 'MM_DENR'
   CHARACTER(7),  PARAMETER :: INPUT_FILE6 = 'MM_FLDN'
   CHARACTER(7),  PARAMETER :: INPUT_FILE7 = 'MM_FRDN'
   CHARACTER(7),  PARAMETER :: INPUT_FILE8 = 'MM_DATR'

    ! Named parameters for energy or J-matrix request or gradient
   LOGICAL, PARAMETER :: GET_ENERGY    = .TRUE.,   &
                         GET_JMAT      = .FALSE.,  &
                         DO_GRADIENT   = .TRUE.,   & 
                         DO_NOGRADIENT = .FALSE. 

    ! Named parameters for different MM methods implemented
   INTEGER(INTK), PARAMETER :: DO_FQ    = 1 ,     &
                               DO_NN    = 2 ,     &
                               DO_BQ    = 3 ,     &
                               DO_NlogN = 4 ,     &
                               DO_FMM   = 5 

    ! Named parameters for different J_matrix call types
   INTEGER(INTK), PARAMETER :: INTTYPE_ONE_EL=1,   &
                               INTTYPE_TWO_EL=2,   &
                               INTTYPE_FULL_J=3,   &
                               INTTYPE_NUC_AT=4,   &
                               INTTYPE_NUC_EL=5
 
    ! Named parameters for T-pair interaction types
   INTEGER(INTK), PARAMETER :: LHS_raw_RHS_raw = 1,     &
                               LHS_raw_RHS_box = 2,     &
                               LHS_box_RHS_raw = 3,     &
                               LHS_box_RHS_box = 4

    ! Named parameters for T-pair interaction search algorithms
   INTEGER(INTK), PARAMETER :: FQ_LOOPS    = 1,     &
                               GRID_SEARCH = 2,     &
                               TREE_SEARCH = 3
    ! and sub-type
   INTEGER(INTK), PARAMETER :: SHAPE_SQUARE     = 1,     &
                               SHAPE_TRIANGULAR = 2

    ! Named parameters for different contractors implemented
   INTEGER(INTK), PARAMETER :: T_CONTRACTOR_DIRECT = 100,        &
                               T_CONTRACTOR_TREE = 101,          &
                               T_CONTRACTOR_DYN = 103,           &
                               T_CONTRACTOR_SCALE = 104,         &
                               T_CONTRACTOR_SCALE_TREE = 1041,   &
                               T_CONTRACTOR_MULTI = 105,         &
                               T_CONTRACTOR_FULL = 106,          &
                               W_CONTRACTOR_DIRECT = 206,        &
                               W_CONTRACTOR_B = 207,             &
                               W_CONTRACTOR_FAST = 208
 
    ! Named parameters for different J-builders (contractors) implemented
   INTEGER(INTK), PARAMETER :: J_BUILDER_DIRECT = 1,        &
                               J_BUILDER_TREE = 2

    ! Named parameters for moments used in contraction
   INTEGER(INTK), PARAMETER :: ELECTRONIC_ONLY = 1,        &
                               NUCLEAR_ONLY    = 2,        &
                               ALL_MOMENTS     = 3,        &
                               ELE_AND_NUC     = 4,        &
                               REG_MIN_AUX     = 5
 
    ! Named parameters for T-pair interaction tests
   INTEGER(INTK), PARAMETER :: TEST_EXT     = 1,        &
                               TEST_R_IJ    = 2,        &
                               TEST_NN_EXT  = 3,        &
                               TEST_NN_R_IJ = 4,        &
                               TEST_FF      = 5,        &
                               TEST_LFF     = 6

    ! Named parameters for auxiliary moments
   INTEGER(INTK), PARAMETER :: USE_RAW_QLM   = 0 ,        &
                               USE_T_SYM_QLM = 1 
 
    ! Named parameters for T and W evaluation schemes
   INTEGER(INTK), PARAMETER :: NON_DIRECT = 0 ,        &
                               RUN_DIRECT = 1 
 
    ! Named parameters for lists available 
   INTEGER(INTK), PARAMETER :: NULL_T_BUFFER  = 1 ,        &
                               NULL_W_BUFFER  = 2 ,        &
                               TREE_T_BUFFER  = 3 ,        &
                               TREE_W_BUFFER  = 4 ,        &
                               SKIP_T_BUFFER  = 5 ,        &
                               SKIP_W_BUFFER  = 6 ,        &
                               MULTI_T_BUFFER = 7 ,        &
                               SCALE_T_BUFFER = 8

    ! Sort orders expected by interaction evaluators.
   INTEGER(INTK), PARAMETER :: NO_SORT         = 0 ,        &
                               SORT_BY_SCALE   = 1,        &
                               SORT_BY_RHS_MMS = 2
 
    ! Named parameters for memory management
   INTEGER(INTK), PARAMETER :: NMEMDIVS = 5
   CHARACTER(7),  PARAMETER :: NSPACE(NMEMDIVS) = (/ 'raw_qlm',   &
                                                     'raw_Vff',   &            
                                                     'box_Vff',   &            
                                                     'Vff_tmp',   &            
                                                     'box_qlm' /)
   INTEGER(INTK), PARAMETER :: MEM_RAW_QLM = 1,        &
                               MEM_RAW_VFF = 2,        &
                               MEM_BOX_VFF = 3,        &
                               MEM_VFF_TMP = 4,        &
                               MEM_BOX_QLM = 5

!===============================================================================

END MODULE mm_global_paras_mod

