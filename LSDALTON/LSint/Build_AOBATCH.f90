!> @file
!> Module containing subroutines for building AOitem (information about AO-shells/batches)
!> \author T. Kjaergaard
!> \date 2008
MODULE BUILDAOBATCH
  use precision
  use TYPEDEF
  use lsmatrix_operations_dense
  use memory_handling

TYPE ORGANIZE
   INTEGER                    :: type,Segment(maxAOangmom),atom
   INTEGER                    :: exponentindex,angmoms,ANGMOM(maxAOangmom)
   INTEGER                    :: CC(maxAOangmom)
END TYPE ORGANIZE
  
TYPE AOorganizer
   INTEGER                    :: nbatches
   TYPE(ORGANIZE),pointer     :: ORG(:)
END TYPE AOORGANIZER

CONTAINS
!> \brief builds the AOitem
!> \author T. Kjaergaard
!> \date 2008
!>
!> The most general routine which builds the AOitem from the BASISSET
!> in the BASISINFO and from the molecule 
!>
SUBROUTINE BUILD_AO(LUPRI,SCHEME,IPRINT,MOLECULE,BASISINFO,AO&
     &,CARTESIAN,UNCONTRACTED,INTNRM)
use molecule_type
implicit none
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> contains all info about the integralscheme requested (thresholds,use cartesian,..)
TYPE(LSINTSCHEME)         :: SCHEME
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASISINFO
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: IPRINT
!> label (Hermite,Cartesian or Default)
Character*(*)             :: Cartesian
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!> if INTERMIDIATE NORMALIZATION shoule be used, USED FOR SCREENING
LOGICAL                   :: INTNRM 
!
TYPE(AOITEM),pointer      :: AOmodel(:)
TYPE(AOorganizer)         :: AOorganize 
INTEGER                   :: nsegments,nAngmom,nbatches,L,A,B
INTEGER                   :: UNIQUEMATRIXES,UNIQUEEXPONENTS,SUM
INTEGER                   :: I,J,K,type,nPrimitives,nContracted
INTEGER                   :: R,aobatches,UNIQUE2,AOsum
INTEGER                   :: orbitalIndex,tmporbitalindex
INTEGER                   :: GHOSTFUNCTIONS,AOmodelbat
INTEGER,pointer           :: MODELTYPES(:)
INTEGER                   :: tmpprimorbitalindex,primorbitalindex,icharge
INTEGER                   :: nbat,SUM1,nMODELEXP
INTEGER,pointer           :: MODELEXP(:)
LOGICAL                   :: NOFAMILY,OLDATOM
Character(len=80)         :: Cartesian2

IF(BASISINFO%natomtypes.EQ.0)CALL LSQUIT('Error BUILD_AO called with empty basis')
AO%natoms = MOLECULE%natoms

NULLIFY(AO%ATOMICnORB)
NULLIFY(AO%ATOMICnBATCH)
ALLOCATE(AO%ATOMICnORB(AO%natoms))
ALLOCATE(AO%ATOMICnBATCH(AO%natoms))
IF(INTNRM)THEN
   IF(.NOT. UNCONTRACTED)THEN
      WRITE(LUPRI,'(2X,A)')'INTNRM REQUIRE UNCONTRACTED TO BE TRUE IN BUILD_AO'
      CALL LSQUIT('INTNRM REQUIRE UNCONTRACTED TO BE TRUE IN BUILD_AO',lupri)
   ENDIF
ENDIF
AO%empty=.FALSE.
NOFAMILY=SCHEME%NOFAMILY
!print*,'BUILD_AOBATCH:  SCHEME%NOFAMILY',SCHEME%NOFAMILY
SELECT CASE(Cartesian(1:4))
CASE('Cart'); Cartesian2 = 'Cartesian'
CASE('Herm'); Cartesian2 = 'Hermite'
CASE('Defa'); !DEFAULT IS TAKEN FROM INPUT 
   IF(SCHEME%DoCartesian)THEN
      Cartesian2 = 'Cartesian'
   ELSE
      Cartesian2 = 'Hermite'
   ENDIF
CASE DEFAULT      
   WRITE (LUPRI,'(/,3A,/)') ' Keyword "',Cartesian,&
        & '" not recognized in Build_AOBATCH'
   print*,'Keyword ',Cartesian,' not recognized'
   CALL LSQUIT('Illegal keyword in Build_AOBATCH',lupri)
END SELECT

call mem_alloc(MODELTYPES,BASISINFO%natomtypes)
SUM=0
L=0
IF(UNCONTRACTED)THEN
   DO I=1,BASISINFO%natomtypes    
      SUM1=0
      L=L+1
      MODELTYPES(I)=L
      DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
         SUM1=SUM1+BASISINFO%ATOMTYPE(I)%SHELL(K)%nprim
      ENDDO
      SUM=MAX(SUM1,SUM)         
   ENDDO
ELSE !DEFAULT
   DO I=1,BASISINFO%natomtypes
      SUM1=0
      L=L+1
      MODELTYPES(I)=L
      DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
         SUM1=SUM1+BASISINFO%ATOMTYPE(I)%SHELL(K)%nsegments
      ENDDO
      SUM=MAX(SUM1,SUM)
   ENDDO
ENDIF
IF(SUM .EQ. 0) THEN
   WRITE(lupri,*)'BUILD_AO CALLED BUT NO BATCHES - SOMETHING WRONG&
        & This may be because you want to use densityfitting but&
        & have not specified a auxillary basis set in the molecule input file'
   CALL FLUSH(LUPRI)
   print*,'BUILD_AO CALLED BUT NO BATCHES - SOMETHING WRONG&
        & This may be because you want to use densityfitting but&
        & have not specified a auxillary basis set in the molecule input file'
   CALL LSQUIT('BUILD_AO CALLED BUT NO BATCHES - SOMETHING WRONG&
        & This may be because you want to use densityfitting but&
        & have not specified a auxillary basis set in the molecule input file',lupri)
ENDIF
!SUM IS NOW THE MAXIMUM NUMBER OF BATCHES A SINGLE ATOM CAN HAVE 
AOmodelbat=L
GHOSTFUNCTIONS=0 !SHOULD BE CHANGED TO ACCOUNT FOR GHOST FUNCTIONS
aobatches=(MOLECULE%natoms-GHOSTFUNCTIONS)*SUM
NULLIFY(AO%BATCH)
ALLOCATE(AO%BATCH(aobatches))
IF(IPRINT .GT. 15)WRITE(lupri,*)aobatches,' aobatches should be more than sufficient'
AOsum=SUM
ALLOCATE(AOmodel(AOmodelbat))
DO I=1,AOmodelbat
   AOmodel(I)%natoms=1
   NULLIFY(AOmodel(I)%ATOMICnORB)
   ALLOCATE(AOmodel(I)%ATOMICnORB(1))
   AOmodel(I)%ATOMICnORB(1)=1
   NULLIFY(AOmodel(I)%ATOMICnBATCH)
   ALLOCATE(AOmodel(I)%ATOMICnBATCH(1))
   NULLIFY(AOmodel(I)%BATCH)
   ALLOCATE(AOmodel(I)%BATCH(SUM))
   AOmodel(I)%nbatches = 0
   NULLIFY(AOmodel(I)%CC)
   NULLIFY(AOmodel(I)%Exponents)
   ALLOCATE(AOmodel(I)%CC(1)) !not used
   AOmodel(I)%nCC = 1 !not used
   AOmodel(I)%nExp = 0
   ALLOCATE(AOmodel(I)%Exponents(SUM))
   CALL lsmat_dense_init(AOmodel(I)%CC(1),1,1)
ENDDO
IF(SUM .EQ. 0) CALL LSQUIT('SUM EQ ZERO SOMETHINGS WRONG',lupri)
nMODELEXP=SUM
call mem_alloc(MODELEXP,SUM)

NULLIFY(AOorganize%ORG)
ALLOCATE(AOorganize%ORG(aobatches))
DO I=1,aobatches
   AOorganize%ORG(I)%exponentindex=0
   DO J=1,maxAOangmom
      AOorganize%ORG(I)%CC(J)=0
   ENDDO
   AOorganize%ORG(I)%atom=0
ENDDO

NULLIFY(AO%CC)
ALLOCATE(AO%CC(SUM*AOmodelbat))
NULLIFY(AO%Exponents)
ALLOCATE(AO%Exponents(SUM*AOmodelbat))

nbatches=0
UNIQUEMATRIXES=0
UNIQUEEXPONENTS=0
AO%nbatches = 0 
AOorganize%nbatches =  0 
orbitalIndex = 1
primorbitalIndex = 1
R = BASISINFO%Labelindex
DO I=1,MOLECULE%natoms   
   IF(IPRINT .GT. 2) THEN
      WRITE(LUPRI,'(2X,A6,I4,2X,A8,F8.5)')'Atom: ',I,'Charge: ',&
           &MOLECULE%ATOM(I)%Charge
   ENDIF
   IF(R.EQ.0)THEN
      ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
      type = BASISINFO%CHARGEINDEX(ICHARGE)
   ELSE
      type = MOLECULE%ATOM(I)%IDtype(R)
   ENDIF
   IF(IPRINT .GT. 2) THEN
      WRITE(LUPRI,'(2X,A10,A20,2X,A6,I4)')'Basisset: ',BASISINFO%ATOMTYPE(type)%NAME(1:20),'Type: ',type
      WRITE(LUPRI,'(2X,A26)')'-------------------------------------'
   ENDIF
   
   !TEST IF THIS TYPE OF ATOM HAS ALREADY BEEN PROCESSED -> OLDATOM=TRUE
   OLDATOM=.FALSE.
   IF(MODELTYPES(type) .NE. 0)THEN
      IF(AOmodel(MODELTYPES(type))%nbatches .NE. 0) OLDATOM=.TRUE.
   ENDIF
   
   IF(OLDATOM) THEN
      IF(IPRINT .GT. 2) WRITE(LUPRI,*)'AS THIS IS AN OLD ATOM WE COPY FROM MODELAOBATCH'
      L=MODELTYPES(type)
      CALL COPY_FROM_MODEL_AO(AOmodel(L),AO,I,MOLECULE,orbitalindex,primorbitalindex,lupri)
      
   ELSE
      IF(IPRINT .GT. 2) WRITE(LUPRI,*)'THIS IS A NEW ATOM SO WE BUILD A MODELAO BATCH'
      L=MODELTYPES(type)
      
      AOmodel(L)%nbatches = 0
      TMPorbitalIndex = 0
      TMPprimorbitalIndex = 0
      
      nAngmom=BASISINFO%ATOMTYPE(type)%nAngmom
      nsegments=0
      DO B=1,nAngmom
         nsegments=nsegments+BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
      ENDDO
      
      DO B=1,nAngmom
         IF(IPRINT .GT. 2) WRITE(LUPRI,'(2X,A,I3,A,I3)')'Angmom:',B,&
              &' of ',nAngmom
         IF(BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments .EQ. 0)THEN
            !NON OF THESE ORBITALS - DO NOT ADD BATCH            
         ELSE
            DO J=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
               IF(IPRINT .GT. 2)  WRITE(LUPRI,'(2X,A,I3,A,I3)')&
                    &'Segment:',J,' of ',&
                    &BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
               
               !TEST IF EXPONENTS is already in AO%Exponents(:) 
               !OR PUT IN ANOTHER WAY 
               !TEST IF IT SHARES EXPONENTS WITH SOME OTHER SEGMENT
               IF(AOmodel(L)%nExp .EQ. 0 )THEN
                  UNIQUE2=0
               ELSE
                  IF(BASISINFO%ATOMTYPE(type)%FAMILY)THEN
                     CALL SHARES_EXPONENTS(LUPRI,IPRINT,AOmodel(L),AOmodel(L)%nExp,BASISINFO,type,B,J,&
                          &UNIQUE2,UNCONTRACTED,MODELEXP,nMODELEXP)
                  ELSE
                     UNIQUE2 = 0
                  ENDIF
               ENDIF
               IF(UNIQUE2 .EQ. 0 .OR. NOFAMILY) THEN
                  !***********************************************************************
                  !*
                  !*   NEW SEGMENT SO NEW BATCH IS ADDED 
                  !*
                  !***********************************************************************
                  
                  IF(IPRINT .GT. 2) WRITE(LUPRI,'(2X,A)')'This is a new segment so we &
                       &increase the number of batches'
                  
                  CALL ADD_BATCH(SCHEME,MOLECULE,BASISINFO,AOmodel(L),AOorganize,I,&
                       &type,B,J,lupri,iprint,nPrimitives,nContracted,CARTESIAN2,&
                       &nbatches,UNCONTRACTED)
                  CALL ADD_SEGMENT_AND_CC(AO,AOmodel(L),AOorganize,BASISINFO,nbatches,nPrimitives,&
                       &nContracted,type,B,J,lupri,iprint,UNIQUEMATRIXES,UNIQUEEXPONENTS,UNCONTRACTED,&
                       &INTNRM,MODELEXP,nMODELEXP)
                  
                  CALL DIRECT_POINTERS(AO,AOmodel(L),AOorganize,nContracted,nPrimitives,nbatches,B,&
                       &lupri,iprint,TMPorbitalIndex,TMPprimorbitalIndex,UNCONTRACTED)
                  
               ELSE
                  !***********************************************************************
                  !*
                  !*   THIS SEGMENT SHARES EXPONENTS WITH AN EARLY SEGMENT
                  !*
                  !*   SO A.FAMILY BASISSETIS IS EMPLOYED
                  !*
                  !***********************************************************************
                  IF(IPRINT .GT. 2) WRITE(LUPRI,'(2X,A,I3)')&
                       &'SHARES EXPONENTS WITH ',UNIQUE2
                  CALL DETERMINE_nbathces(AOorganize,I,UNIQUE2,nbat)
                  CALL EXPAND_BATCH(AOmodel(L),AOorganize,basisinfo,type,B,J,nbat,&
                       &A,lupri,UNCONTRACTED)
                  
                  CALL ADD_CC(AO,AOorganize,nbat,nPrimitives,nContracted,lupri,iprint,&
                       &BASISINFO,type,B,J,A,UNIQUEMATRIXES,UNCONTRACTED,INTNRM)
                  
                  CALL DIRECT_CC_POINTERS(AO,AOmodel(L),AOorganize,nContracted,&
                       &nPrimitives,A,nbat,B,lupri,iprint,TMPorbitalIndex,TMPprimorbitalindex,UNCONTRACTED)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      CALL SET_AO_EXTENT(AOmodel(L),SCHEME)
      AOmodel(L)%ATOMICnBATCH(1)=AOmodel(L)%nbatches
      IF(IPRINT .GT. 40)THEN
         WRITE(LUPRI,'(2X,A)')'The Model AO for this atom of this type and basisset'
         CALL PRINT_AO(LUPRI,AOmodel(L))
      ENDIF
      CALL COPY_FROM_MODEL_AO(AOmodel(L),AO,I,MOLECULE,orbitalindex,primorbitalindex,lupri)
   ENDIF
ENDDO

if (.not.ASSOCIATED(AOorganize%ORG)) then
  print*,'memory previously released!!'
  STOP 'Error in BUILD_AO - memory previously released'
endif



DO I =1,AOmodelbat
   call free_AOitem(lupri,AOmodel(I))
ENDDO
DEALLOCATE(AOmodel)
DEALLOCATE(AOorganize%ORG)
NULLIFY(AOorganize%ORG)
call mem_dealloc(MODELTYPES)
call mem_dealloc(MODELEXP)
!DEALLOCATE(MODELTYPES)

!AO%nbatches=nbatches
AO%nCC=UNIQUEMATRIXES
AO%nExp=UNIQUEEXPONENTS
AO%nbast = orbitalindex-1
IF(IPRINT .GT. 20) THEN
   WRITE(LUPRI,*)'BUILD_AOBATCH:  PRINTING FINAL AO BATCH  UNCONTRACTED',UNCONTRACTED
   CALL PRINT_AO(LUPRI,AO)
ENDIF

!CALL DETERMINE_AOBATCH_MEM(AO)

END SUBROUTINE BUILD_AO

!!$SUBROUTINE DETERMINE_AOBATCH_MEM(AO)
!!$IMPLICIT NONE
!!$TYPE(AOITEM)              :: AO
!!$INTEGER(KIND=long)        :: alloc_memory
!!$INTEGER                    :: I
!!$alloc_memory = 0
!!$alloc_memory = alloc_memory + sizeof(AO%EMPTY) 
!!$alloc_memory = alloc_memory + sizeof(AO%nbatches) 
!!$alloc_memory = alloc_memory + sizeof(AO%nCC) 
!!$alloc_memory = alloc_memory + sizeof(AO%nExp) 
!!$DO I = 1,AO%nbatches
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%type)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%spherical)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%CENTER)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nPrimitives)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%maxContracted)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%maxAngmom)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%pExponents)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nAngmom)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%extent)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%ANGMOM)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nContracted)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%startOrbital)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%startprimOrbital)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nOrbComp)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nPrimOrbComp)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nOrbitals)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%pCC)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%CCindex)
!!$ENDDO
!!$DO I = 1,AO%nCC
!!$   alloc_memory = alloc_memory+sizeof(AO%CC(I)%elms)+sizeof(AO%CC(I)%nrow)+sizeof(AO%CC(I)%ncol)+sizeof(AO%CC(I)%complex)
!!$ENDDO
!!$DO I = 1,AO%nEXP
!!$   alloc_memory = alloc_memory+sizeof(AO%Exponents(I)%elms)+sizeof(AO%Exponents(I)%nrow)+sizeof(AO%Exponents(I)%ncol)+sizeof(AO%Exponents(I)%complex)
!!$ENDDO
!!$END SUBROUTINE DETERMINE_AOBATCH_MEM

!> \brief builds an empty AOitem
!> \author T. Kjaergaard
!> \date 2008
!>
!> build an empty AOitem used, for when calculating non 4center integrals
!>
SUBROUTINE BUILD_EMPTY_AO(AO,LUPRI)
IMPLICIT NONE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the logical unit number for the output file
INTEGER                   :: LUPRI
AO%natoms = 1                       
NULLIFY(AO%ATOMICnORB)
NULLIFY(AO%ATOMICnBATCH)
ALLOCATE(AO%ATOMICnORB(AO%natoms))  
ALLOCATE(AO%ATOMICnBATCH(AO%natoms))
AO%ATOMICnORB(1)=1                  
AO%ATOMICnBATCH(1)=1                
AO%nbast = 1                

AO%empty=.TRUE.
AO%nbatches=1
AO%nCC=1
AO%nExp=1
NULLIFY(AO%BATCH)
ALLOCATE(AO%BATCH(1))   
NULLIFY(AO%CC)
ALLOCATE(AO%CC(1))
NULLIFY(AO%Exponents)
ALLOCATE(AO%Exponents(1))
CALL lsmat_dense_init(AO%CC(1),1,1)
CALL lsmat_dense_init(AO%Exponents(1),1,1)
AO%CC(1)%elms(1)=1
AO%Exponents(1)%elms(1)=0
AO%BATCH(1)%type_Empty = .TRUE.
AO%BATCH(1)%type_Hermite = .FALSE.
AO%BATCH(1)%type_Cartesian = .FALSE.
AO%BATCH(1)%type_Nucleus = .FALSE.
AO%BATCH(1)%spherical=.true.
AO%BATCH(1)%atom=1
AO%BATCH(1)%batch=1
AO%BATCH(1)%CENTER(1)=0.d0
AO%BATCH(1)%CENTER(2)=0.d0
AO%BATCH(1)%CENTER(3)=0.d0
AO%BATCH(1)%nPrimitives=1
AO%BATCH(1)%maxContracted=1
AO%BATCH(1)%maxAngmom=0
AO%BATCH(1)%pExponents => AO%Exponents(1)
AO%BATCH(1)%nAngmom=1
AO%BATCH(1)%extent=0.d0
AO%BATCH(1)%ANGMOM(1)=0
AO%BATCH(1)%nContracted(1)=1
AO%BATCH(1)%startOrbital(1)=1
AO%BATCH(1)%startprimOrbital(1)=1
AO%BATCH(1)%nOrbComp(1)=1
AO%BATCH(1)%nPrimOrbComp(1)=1
AO%BATCH(1)%nOrbitals(1)=1
AO%BATCH(1)%pCC(1)%p => AO%CC(1)
AO%BATCH(1)%CCindex(1) = 0

!WRITE(LUPRI,*)'BUILD_AOBATCH:  PRINTING FINAL AO BATCH'
!CALL PRINT_AO(LUPRI,AO)
!CALL DETERMINE_AOBATCH_MEM(AO)

END SUBROUTINE BUILD_EMPTY_AO

!> \brief builds an empty nuclear AOitem
!> \author T. Kjaergaard
!> \date 2008
!>
!> build an empty nuclear AOitem, used for nuclear attraction integrals
!>
SUBROUTINE BUILD_EMPTY_NUCLEAR_AO(AO,MOLECULE,LUPRI)
use molecule_type
implicit none
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!
INTEGER                   :: aobatches,I,ncharges,J,NEWCHARGE
REAL(REALK),pointer       :: CHARGES(:)

AO%natoms = 1                            
NULLIFY(AO%ATOMICnORB)       
NULLIFY(AO%ATOMICnBATCH)     
ALLOCATE(AO%ATOMICnORB(AO%natoms))       
ALLOCATE(AO%ATOMICnBATCH(AO%natoms))     
AO%ATOMICnORB(1)=0                       
AO%ATOMICnBATCH(1)=1                     
AO%nbast = 1                     

AO%empty=.TRUE.
aobatches=MOLECULE%natoms
AO%nbatches=aobatches
!AO%nCC=1
AO%nExp=1

NULLIFY(AO%BATCH)
ALLOCATE(AO%BATCH(aobatches))

NULLIFY(AO%Exponents)
ALLOCATE(AO%Exponents(1))

CALL lsmat_dense_init(AO%Exponents(1),1,1)

call mem_alloc(CHARGES,MOLECULE%natoms)
!ALLOCATE(CHARGES(MOLECULE%natoms))
CHARGES(1)=MOLECULE%ATOM(1)%Charge
ncharges=1
DO I=2,MOLECULE%natoms
   NEWCHARGE=0
   DO J=1,ncharges
      IF(ABS(MOLECULE%ATOM(I)%Charge-charges(J) ).LE.1.d-10 ) NEWCHARGE=J
   ENDDO
   IF(NEWCHARGE .EQ. 0)THEN
      ncharges=ncharges+1
      CHARGES(ncharges)=MOLECULE%ATOM(I)%Charge
   ENDIF
ENDDO

AO%nCC=nCharges

NULLIFY(AO%CC)
ALLOCATE(AO%CC(nCharges))

DO I=1,nCharges
   CALL lsmat_dense_init(AO%CC(I),1,1)
   AO%CC(I)%elms(1)=-CHARGES(I)
ENDDO

DO I=1,MOLECULE%natoms
   DO J=1,ncharges
      IF(ABS(MOLECULE%ATOM(I)%Charge-charges(J)).LE.1.d-10) NEWCHARGE=J
   ENDDO
   AO%BATCH(I)%pCC(1)%p => AO%CC(NEWCHARGE)   
   AO%BATCH(I)%CCindex(1) = NEWCHARGE   
ENDDO

call mem_dealloc(CHARGES)
!DEALLOCATE(CHARGES)

AO%Exponents(1)%elms(1)=0.d0

DO I=1,MOLECULE%natoms
   AO%BATCH(I)%type_Nucleus = .TRUE.
   AO%BATCH(I)%TYPE_Empty = .FALSE.
   AO%BATCH(I)%TYPE_Hermite = .FALSE.
   AO%BATCH(I)%TYPE_Cartesian = .FALSE.
   AO%BATCH(I)%spherical=.true.
   AO%BATCH(I)%atom=I
   AO%BATCH(I)%batch=1
   AO%BATCH(I)%CENTER(1)=MOLECULE%ATOM(I)%CENTER(1)
   AO%BATCH(I)%CENTER(2)=MOLECULE%ATOM(I)%CENTER(2)
   AO%BATCH(I)%CENTER(3)=MOLECULE%ATOM(I)%CENTER(3)
   AO%BATCH(I)%nPrimitives=1
   AO%BATCH(I)%maxContracted=1
   AO%BATCH(I)%maxAngmom=0
   AO%BATCH(I)%pExponents => AO%Exponents(1)
   AO%BATCH(I)%nAngmom=1
   AO%BATCH(I)%extent=0.d0
   AO%BATCH(I)%ANGMOM(1)=0
   AO%BATCH(I)%nContracted(1)=1
   AO%BATCH(I)%startOrbital(1)=1
   AO%BATCH(I)%startprimOrbital(1)=I
   AO%BATCH(I)%nOrbComp(1)=1
   AO%BATCH(I)%nPrimOrbComp(1)=1
   AO%BATCH(I)%nOrbitals(1)=1
ENDDO

!WRITE(LUPRI,*)'BUILD_AOBATCH:  PRINTING FINAL AO BATCH'
!CALL PRINT_AO(LUPRI,AO)
!CALL DETERMINE_AOBATCH_MEM(AO)

END SUBROUTINE BUILD_EMPTY_NUCLEAR_AO

!> \brief builds an AOitem with a single Orbital
!> \author T. Kjaergaard
!> \date 2008
!>
!> builds an AOitem with a single Orbital used for coupled-cluster
!> calculations
!>
SUBROUTINE BUILD_SINGLE_ORBBATCH_AO(LUPRI,SCHEME,IPRINT,MOLECULE,&
     &BASISINFO,AO,CARTESIAN,basfunc,dim,OrbitalComponent)
use molecule_type
implicit none
!> the logical unit number for the output file
INTEGER,intent(in)                   :: LUPRI
!> contains all info about the integralscheme requested (thresholds,use cartesian,..)
TYPE(LSINTSCHEME),intent(in) :: SCHEME
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in)           :: IPRINT
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO),intent(in):: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in):: BASISINFO
!> the AOitem to be build
TYPE(AOITEM),intent(inout):: AO
!> label (Hermite,Cartesian or Default)
Character*(*)             :: Cartesian
!> which basisfunction the AOitem should be build from.
Integer,intent(in)        :: basfunc
!> dimension of the orbital batch (1 for s, 3 for p, 5 for d,etc)
Integer,intent(out)       :: dim
! the orbital component requested (for instance 2 could be the Py component)
Integer,intent(out)       :: OrbitalComponent
!
Character(len=80)         :: Cartesian2
Integer  :: natoms,type,orbitalindex,ncont,nprim
Integer  :: I,K,KA,AA,L,ANGMOM,ORBITAL,iATOM,ATOMtype2,R,ATOMnprim
Integer  :: SEGMENT,norb,icharge
Logical  :: DONE
INTEGER,pointer :: ATOMtype(:)

SELECT CASE(Cartesian(1:4))
CASE('Cart'); Cartesian2 = 'Cartesian'
CASE('Herm'); Cartesian2 = 'Hermite'
CASE('Defa'); !DEFAULT IS TAKEN FROM INPUT 
   IF(SCHEME%DoCartesian)THEN
      Cartesian2 = 'Cartesian'
   ELSE
      Cartesian2 = 'Hermite'
   ENDIF
CASE DEFAULT      
   WRITE (LUPRI,'(/,3A,/)') ' Keyword "',Cartesian,&
        & '" not recognized in Build_AOBATCH'
   print*,'Keyword ',Cartesian,' not recognized'
   CALL LSQUIT('Illegal keyword in Build_AOBATCH',lupri)
END SELECT

AO%empty=.FALSE.

natoms = MOLECULE%nAtoms
call mem_alloc(ATOMtype,natoms)
R = BASISINFO%LabelIndex
IF(R.EQ.0)THEN
   DO I=1,nAtoms
      ICHARGE = INT(MOLECULE%ATOM(I)%Charge)
      ATOMtype(I) = BASISINFO%ChargeIndex(ICHARGE)
   ENDDO
ELSE
   DO I=1,nAtoms
      ATOMtype(I)=MOLECULE%ATOM(I)%IDtype(R)
   ENDDO
ENDIF
orbitalindex = 0
DONE=.FALSE.
DO I=1,nAtoms
   type=ATOMtype(I)
   DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
      norb=BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
      DO KA=1,BASISINFO%ATOMTYPE(type)%SHELL(K)%nsegments
         ncont=BASISINFO%ATOMTYPE(type)%SHELL(K)%SEGMENT(KA)%ncol
         nprim=BASISINFO%ATOMTYPE(type)%SHELL(K)%SEGMENT(KA)%nrow
         DO L=1,ncont
            DO AA=1,(2*(K-1)+1)
               orbitalindex=orbitalindex+1
!               primorbitalindex=primorbitalindex+nprim
               IF(orbitalindex .EQ. BASFUNC) THEN
                  ANGMOM=K
                  ORBITAL=L
                  iATOM=I
                  ATOMnprim=nprim
                  ATOMtype2=type
                  DONE=.TRUE.
                  OrbitalComponent=AA
                  SEGMENT=KA
               ENDIF
               IF(DONE)EXIT
            ENDDO
            IF(DONE)EXIT
         ENDDO
         IF(DONE)EXIT
      ENDDO
      IF(DONE)EXIT
!      primorbitalIndex = primorbitalIndex + nprim*(2*(K-1)+1)
   ENDDO
   IF(DONE)EXIT
ENDDO

AO%nbatches=1

AO%nCC=1
AO%nExp=1
NULLIFY(AO%BATCH)
ALLOCATE(AO%BATCH(1))   
NULLIFY(AO%CC)
ALLOCATE(AO%CC(1))
NULLIFY(AO%Exponents)
ALLOCATE(AO%Exponents(1))

CALL lsmat_dense_init(AO%CC(1),ATOMnPrim,1)
CALL lsmat_dense_init(AO%Exponents(1),ATOMnPrim,1)

DO I=1,ATOMnPrim
   AO%CC(1)%elms(I)=BASISINFO%ATOMTYPE(ATOMtype2)&
        &%SHELL(ANGMOM)%SEGMENT(SEGMENT)%elms(I+(ORBITAL-1)*ATOMnPrim)
ENDDO
call dcopy (nPrim,BASISINFO%ATOMTYPE(ATOMtype2)&
     &%SHELL(ANGMOM)%SEGMENT(SEGMENT)%Exponents,1,&
     &AO%Exponents(1)%elms,1)


!AO%BATCH(1)%type=Cartesian2
AO%BATCH(1)%type_cartesian  = Cartesian2(1:4).EQ.'Cart'
AO%BATCH(1)%type_hermite  = Cartesian2(1:4).EQ.'Herm'
AO%BATCH(1)%type_Empty  = .FALSE. 
AO%BATCH(1)%type_Nucleus  = .FALSE.

AO%BATCH(1)%atom=Iatom
AO%BATCH(1)%batch=1

AO%BATCH(1)%spherical=SCHEME%DoSpherical

AO%BATCH(1)%CENTER(1)=MOLECULE%ATOM(iATOM)%CENTER(1)
AO%BATCH(1)%CENTER(2)=MOLECULE%ATOM(iATOM)%CENTER(2)
AO%BATCH(1)%CENTER(3)=MOLECULE%ATOM(iATOM)%CENTER(3)

AO%BATCH(1)%nPrimitives=nPrim
AO%BATCH(1)%maxContracted=nCont
AO%BATCH(1)%maxAngmom=ANGMOM-1

AO%BATCH(1)%pExponents => AO%Exponents(1)
AO%BATCH(1)%pCC(1)%p => AO%CC(1)
AO%BATCH(1)%CCindex(1) = 1
AO%BATCH(1)%nAngmom=1
AO%BATCH(1)%extent = getExtent(AO%BATCH(1)%pExponents%elms,AO%BATCH(1)%pCC(1)%p%elms,&
     &                         nPrim,nCont,ANGMOM-1,SCHEME%OD_SCREEN,SCHEME%THRESHOLD*SCHEME%OD_THRESHOLD)
AO%BATCH(1)%ANGMOM(1)=ANGMOM-1
AO%BATCH(1)%nContracted(1)=1
AO%BATCH(1)%startOrbital(1)=1

IF (AO%BATCH(1)%spherical) THEN
   AO%BATCH(1)%nOrbComp(1) = 2*ANGMOM-1
ELSE
   AO%BATCH(1)%nOrbComp(1) = ANGMOM*(ANGMOM+1)/2
ENDIF
dim=AO%BATCH(1)%nOrbComp(1)
AO%BATCH(1)%nOrbitals(1)=AO%BATCH(1)%nOrbComp(1)
AO%BATCH(1)%nPrimOrbComp(1)= ANGMOM*(ANGMOM+1)/2
AO%BATCH(1)%startprimOrbital(1)=1

AO%natoms = 1         
NULLIFY(AO%ATOMICnORB)       
NULLIFY(AO%ATOMICnBATCH)                  
ALLOCATE(AO%ATOMICnORB(AO%natoms)) 
ALLOCATE(AO%ATOMICnBATCH(AO%natoms))    
AO%ATOMICnORB(1)=AO%BATCH(1)%nOrbComp(1)
AO%nbast = AO%BATCH(1)%nOrbComp(1)
AO%ATOMICnBATCH(1)=1

IF(IPRINT .GT. 20) THEN
   WRITE(LUPRI,*)'BUILD_SINGLE_ORB_AO: PRINTING AO BATCH'
   CALL PRINT_AO(LUPRI,AO)
ENDIF

call mem_dealloc(ATOMtype)
!CALL DETERMINE_AOBATCH_MEM(AO)

END SUBROUTINE BUILD_SINGLE_ORBBATCH_AO

!> \brief builds an AOitem with a single shellbatch (an S, or PxPyPz, etc.)
!> \author T. Kjaergaard
!> \date 2008
!>
!> builds an AOitem with a single Orbitalbatch used for coupled-cluster
!> calculations, where the 4 center 2 electron integrals are calculated 
!> in batches to avoid memory problems 
!>
SUBROUTINE BUILD_SINGLE_SHELLBATCH_AO(LUPRI,SCHEME,IPRINT,MOLECULE,&
     &BASISINFO,AO_output,CARTESIAN,UNCONTRACTED,INTNRM,RequestedBatchIndex,dim)
use molecule_type
implicit none
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> contains all info about the integralscheme requested (thresholds,use cartesian,..)
TYPE(LSINTSCHEME)         :: SCHEME
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: IPRINT
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASISINFO
!> the AOitem to be build
TYPE(AOITEM)              :: AO_output
!> label (Hermite,Cartesian or Default)
Character*(*)             :: Cartesian
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!> if INTERMIDIATE NORMALIZATION shoule be used, USED FOR SCREENING
LOGICAL                   :: INTNRM 
!> the requested batch index to build the AOitem from
INTEGER                   :: RequestedBatchIndex
!> dimension of the orbital batch (1 for s, 3 for p, 5 for d,etc)
INTEGER                   :: dim
!
TYPE(AOITEM),pointer      :: AOmodel(:)
TYPE(AOorganizer)         :: AOorganize 
INTEGER                   :: nbatchesDUM
INTEGER                   :: nsegments,nAngmom,nbatches,L,A,B
INTEGER                   :: UNIQUEMATRIXES,UNIQUEEXPONENTS,SUM,IBAT
INTEGER                   :: I,I2,J,K,type,nPrimitives,nContracted
INTEGER                   :: R,aobatches,UNIQUE2,AOsum,iangmom
INTEGER                   :: orbitalIndex, nOrbitals,tmporbitalindex
INTEGER                   :: GHOSTFUNCTIONS,AOmodelbat
TYPE(AOITEM)              :: AO
INTEGER,pointer           :: MODELTYPES(:)
INTEGER                   :: tmpprimorbitalindex,primorbitalindex,icharge
INTEGER                   :: nbat,SUM1,nMODELEXP,iBB,jBB
INTEGER,pointer           :: MODELEXP(:)
LOGICAL                   :: NOFAMILY,OLDATOM
Character(len=80)         :: Cartesian2

AO_output%empty=.FALSE.
NOFAMILY=SCHEME%NOFAMILY
IF(INTNRM)THEN
   IF(.NOT. UNCONTRACTED)THEN
      WRITE(LUPRI,'(2X,A)')'INTNRM REQUIRE UNCONTRACTED TO BE TRUE IN BUILD_AO'
      CALL LSQUIT('INTNRM REQUIRE UNCONTRACTED TO BE TRUE IN BUILD_AO',lupri)
   ENDIF
ENDIF
SELECT CASE(Cartesian(1:4))
CASE('Cart'); Cartesian2 = 'Cartesian'
CASE('Herm'); Cartesian2 = 'Hermite'
CASE('Defa'); !DEFAULT IS TAKEN FROM INPUT 
   IF(SCHEME%DoCartesian)THEN
      Cartesian2 = 'Cartesian'
   ELSE
      Cartesian2 = 'Hermite'
   ENDIF
CASE DEFAULT      
   WRITE (LUPRI,'(/,3A,/)') ' Keyword "',Cartesian,&
        & '" not recognized in Build_AOBATCH'
   print*,'Keyword ',Cartesian,' not recognized'
   CALL LSQUIT('Illegal keyword in Build_AOBATCH',lupri)
END SELECT

call mem_alloc(MODELTYPES,BASISINFO%natomtypes)
SUM=0
L=0
IF(UNCONTRACTED)THEN
   DO I=1,BASISINFO%natomtypes    
      SUM1=0
      L=L+1
      MODELTYPES(I)=L
      DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
         SUM1=SUM1+BASISINFO%ATOMTYPE(I)%SHELL(K)%nprim
      ENDDO
      SUM=MAX(SUM1,SUM)         
   ENDDO
ELSE !DEFAULT
   DO I=1,BASISINFO%natomtypes
      SUM1=0
      L=L+1
      MODELTYPES(I)=L
      DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
         SUM1=SUM1+BASISINFO%ATOMTYPE(I)%SHELL(K)%nsegments
      ENDDO
      SUM=MAX(SUM1,SUM)
   ENDDO
ENDIF
IF(SUM .EQ. 0) THEN
   WRITE(lupri,*)'BUILD_AO CALLED BUT NO BATCHES - SOMETHING WRONG&
        & This may be because you want to use densityfitting but&
        & have not specified a auxillary basis set in the molecule input file'
   CALL FLUSH(LUPRI)
   print*,'BUILD_AO CALLED BUT NO BATCHES - SOMETHING WRONG&
        & This may be because you want to use densityfitting but&
        & have not specified a auxillary basis set in the molecule input file'
   CALL LSQUIT('BUILD_AO CALLED BUT NO BATCHES - SOMETHING WRONG&
        & This may be because you want to use densityfitting but&
        & have not specified a auxillary basis set in the molecule input file',lupri)
ENDIF
!SUM IS NOW THE MAXIMUM NUMBER OF BATCHES A SINGLE ATOM CAN HAVE 
AOmodelbat=L
GHOSTFUNCTIONS=0 !SHOULD BE CHANGED TO ACCOUNT FOR GHOST FUNCTIONS
aobatches=(MOLECULE%natoms-GHOSTFUNCTIONS)*SUM
AOsum=SUM
ALLOCATE(AOmodel(AOmodelbat))
DO I=1,AOmodelbat
   NULLIFY(AOmodel(I)%ATOMICnORB)
   ALLOCATE(AOmodel(I)%ATOMICnORB(1))
   AOmodel(I)%ATOMICnORB(1)=1
   NULLIFY(AOmodel(I)%ATOMICnBATCH)
   ALLOCATE(AOmodel(I)%ATOMICnBATCH(1))
   AOmodel(I)%ATOMICnBATCH(1)=1
   NULLIFY(AOmodel(I)%BATCH)
   ALLOCATE(AOmodel(I)%BATCH(SUM))
   AOmodel(I)%nbatches = 0
   NULLIFY(AOmodel(I)%CC)
   NULLIFY(AOmodel(I)%Exponents)
   ALLOCATE(AOmodel(I)%CC(1)) !not used
   AOmodel(I)%nCC = 1 !not used
   AOmodel(I)%nExp = 0
   ALLOCATE(AOmodel(I)%Exponents(SUM))
   CALL lsmat_dense_init(AOmodel(I)%CC(1),1,1)
ENDDO
IF(SUM .EQ. 0) CALL LSQUIT('SUM EQ ZERO SOMETHINGS WRONG',lupri)
nMODELEXP=SUM
call mem_alloc(MODELEXP,SUM)

NULLIFY(AOorganize%ORG)
ALLOCATE(AOorganize%ORG(aobatches))
DO I=1,aobatches
   AOorganize%ORG(I)%exponentindex=0
   DO J=1,maxAOangmom
      AOorganize%ORG(I)%CC(J)=0
   ENDDO
   AOorganize%ORG(I)%atom=0
ENDDO

NULLIFY(AO%CC)
ALLOCATE(AO%CC(SUM*AOmodelbat))
NULLIFY(AO%Exponents)
ALLOCATE(AO%Exponents(SUM*AOmodelbat))

nbatches=0
UNIQUEMATRIXES=0
UNIQUEEXPONENTS=0
AOorganize%nbatches =  0 
orbitalIndex = 1
primorbitalIndex = 1
R = BASISINFO%Labelindex
DO I=1,MOLECULE%natoms   
   IF(R.EQ.0)THEN
      ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
      type = BASISINFO%CHARGEINDEX(ICHARGE)
   ELSE
      type = MOLECULE%ATOM(I)%IDtype(R)
   ENDIF
   !TEST IF THIS TYPE OF ATOM HAS ALREADY BEEN PROCESSED -> OLDATOM=TRUE
   OLDATOM=.FALSE.
   IF(MODELTYPES(type) .NE. 0)THEN
      IF(AOmodel(MODELTYPES(type))%nbatches .NE. 0) OLDATOM=.TRUE.
   ENDIF
   IF(OLDATOM) THEN
      L=MODELTYPES(type)
      DO IBAT = 1,AOmodel(L)%nbatches
         nbatches = nbatches + 1
         IF(nbatches .EQ. RequestedBatchIndex)THEN
            call build_aoitem_from_aobatch(LUPRI,AOmodel(L)%BATCH(IBAT),MOLECULE,I,AO_OUTPUT)
            dim = 0
            DO iangmom = 1,AO_output%BATCH(1)%nangmom
               dim = dim + AO_output%BATCH(1)%nOrbitals(iangmom)
            ENDDO
            DO I2 =1,AOmodelbat
               call free_AOitem(lupri,AOmodel(I2))
            ENDDO
            DEALLOCATE(AOmodel)
            DEALLOCATE(AOorganize%ORG)
            NULLIFY(AOorganize%ORG)
            call mem_dealloc(MODELTYPES)
            call mem_dealloc(MODELEXP)
            DO I2=1,UNIQUEMATRIXES
               CALL LSMAT_DENSE_FREE(AO%CC(I2))
            ENDDO
            DO I2=1,UNIQUEEXPONENTS
               CALL LSMAT_DENSE_FREE(AO%Exponents(I2))
            ENDDO 
            DEALLOCATE(AO%CC)
            NULLIFY(AO%CC)
            DEALLOCATE(AO%Exponents)
            NULLIFY(AO%Exponents)
            AO_OUTPUT%natoms = 1         
            NULLIFY(AO_OUTPUT%ATOMICnORB)       
            NULLIFY(AO_OUTPUT%ATOMICnBATCH)             
            ALLOCATE(AO_OUTPUT%ATOMICnORB(AO_OUTPUT%natoms)) 
            ALLOCATE(AO_OUTPUT%ATOMICnBATCH(AO_OUTPUT%natoms))
            norbitals = 0
            IBB=1
            DO JBB = 1,AO_OUTPUT%BATCH(IBB)%nangmom
               norbitals = norbitals + AO_OUTPUT%BATCH(IBB)%nOrbitals(JBB)
            ENDDO
            AO_OUTPUT%ATOMICnORB(1)=norbitals
            AO_OUTPUT%nbast = norbitals
            AO_OUTPUT%ATOMICnBATCH(1)=1
            IF(iprint.GT.20)then
               write(lupri,*)'write from single shell  1'
               CALL PRINT_AO(LUPRI,AO_OUTPUT)
            endif
            return
         ENDIF
      ENDDO
   ELSE
      L=MODELTYPES(type)      
      AOmodel(L)%nbatches = 0
      TMPorbitalIndex = 0
      TMPprimorbitalIndex = 0
      nAngmom=BASISINFO%ATOMTYPE(type)%nAngmom
      nsegments=0
      DO B=1,nAngmom
         nsegments=nsegments+BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
      ENDDO
      DO B=1,nAngmom
         IF(IPRINT .GT. 2) WRITE(LUPRI,'(2X,A,I3,A,I3)')'Angmom:',B,&
              &' of ',nAngmom
         IF(BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments .EQ. 0)THEN
            !NON OF THESE ORBITALS - DO NOT ADD BATCH            
         ELSE
            DO J=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
               !TEST IF EXPONENTS is already in AO%Exponents(:) 
               !OR PUT IN ANOTHER WAY 
               !TEST IF IT SHARES EXPONENTS WITH SOME OTHER SEGMENT
               IF(AOmodel(L)%nExp .EQ. 0 )THEN
                  UNIQUE2=0
               ELSE
                  IF(BASISINFO%ATOMTYPE(type)%FAMILY)THEN
                     CALL SHARES_EXPONENTS(LUPRI,IPRINT,AOmodel(L),AOmodel(L)%nExp,BASISINFO,type,B,J,&
                     &UNIQUE2,UNCONTRACTED,MODELEXP,nMODELEXP)
                  ELSE
                     UNIQUE2 = 0
                  ENDIF
               ENDIF
               IF(UNIQUE2 .EQ. 0 .OR. NOFAMILY) THEN
                  !***********************************************************************
                  !*   NEW SEGMENT SO NEW BATCH IS ADDED TO AOMODEL
                  !***********************************************************************
                  CALL ADD_BATCH(SCHEME,MOLECULE,BASISINFO,AOmodel(L),AOorganize,I,&
                       &type,B,J,lupri,iprint,nPrimitives,nContracted,CARTESIAN2,&
                       &nbatchesDUM,UNCONTRACTED)
                  IF(UNCONTRACTED)THEN
                     nbatches = nbatches + BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
                  ELSE
                     nbatches = nbatches + 1
                  ENDIF
                  CALL ADD_SEGMENT_AND_CC(AO,AOmodel(L),AOorganize,BASISINFO,nbatchesDUM,nPrimitives,&
                       &nContracted,type,B,J,lupri,iprint,UNIQUEMATRIXES,UNIQUEEXPONENTS,UNCONTRACTED,&
                       &INTNRM,MODELEXP,nMODELEXP)
                  CALL DIRECT_POINTERS(AO,AOmodel(L),AOorganize,nContracted,nPrimitives,nbatchesDUM,B,&
                       &lupri,iprint,TMPorbitalIndex,TMPprimorbitalIndex,UNCONTRACTED)
               ELSE
                  !***********************************************************************
                  !*   THIS SEGMENT SHARES EXPONENTS WITH AN EARLY SEGMENT
                  !*   SO A.FAMILY BASISSETIS IS EMPLOYED
                  !***********************************************************************
                  CALL DETERMINE_nbathces(AOorganize,I,UNIQUE2,nbat)
                  CALL EXPAND_BATCH(AOmodel(L),AOorganize,basisinfo,type,B,J,nbat,&
                       &A,lupri,UNCONTRACTED)                  
                  CALL ADD_CC(AO,AOorganize,nbat,nPrimitives,nContracted,lupri,iprint,&
                       &BASISINFO,type,B,J,A,UNIQUEMATRIXES,UNCONTRACTED,INTNRM)
                  CALL DIRECT_CC_POINTERS(AO,AOmodel(L),AOorganize,nContracted,&
                       &nPrimitives,A,nbat,B,lupri,iprint,TMPorbitalIndex,TMPprimorbitalindex,UNCONTRACTED)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      CALL SET_AO_EXTENT(AOmodel(L),SCHEME)
      DO IBAT = 1,AOmodel(L)%nbatches
         IF(nbatches-AOmodel(L)%nbatches+IBAT .EQ. RequestedBatchIndex)THEN
            call build_aoitem_from_aobatch(LUPRI,AOmodel(L)%BATCH(IBAT),MOLECULE,I,AO_OUTPUT)
            dim = 0
            DO iangmom = 1,AO_output%BATCH(1)%nangmom
               dim = dim + AO_output%BATCH(1)%nOrbitals(iangmom)
            ENDDO
            DO I2 =1,AOmodelbat
               call free_AOitem(lupri,AOmodel(I2))
            ENDDO
            DEALLOCATE(AOmodel)
            DEALLOCATE(AOorganize%ORG)
            NULLIFY(AOorganize%ORG)
            call mem_dealloc(MODELTYPES)
            call mem_dealloc(MODELEXP)
            DO I2=1,UNIQUEMATRIXES
               CALL LSMAT_DENSE_FREE(AO%CC(I2))
            ENDDO
            DO I2=1,UNIQUEEXPONENTS
               CALL LSMAT_DENSE_FREE(AO%Exponents(I2))
            ENDDO 
            DEALLOCATE(AO%CC)
            NULLIFY(AO%CC)
            DEALLOCATE(AO%Exponents)
            NULLIFY(AO%Exponents)
            AO_OUTPUT%natoms = 1         
            NULLIFY(AO_OUTPUT%ATOMICnORB)       
            NULLIFY(AO_OUTPUT%ATOMICnBATCH)             
            ALLOCATE(AO_OUTPUT%ATOMICnORB(AO_OUTPUT%natoms)) 
            ALLOCATE(AO_OUTPUT%ATOMICnBATCH(AO_OUTPUT%natoms))
            norbitals = 0
            IBB=1
            DO JBB = 1,AO_OUTPUT%BATCH(IBB)%nangmom
               norbitals = norbitals + AO_OUTPUT%BATCH(IBB)%nOrbitals(JBB)
            ENDDO
            AO_OUTPUT%ATOMICnORB(1)=norbitals
            AO_OUTPUT%nbast = norbitals
            AO_OUTPUT%ATOMICnBATCH(1)=1
            IF(iprint.GT.20)then
               write(lupri,*)'write from single shell  1'
               CALL PRINT_AO(LUPRI,AO_OUTPUT)
            endif
            RETURN
         ENDIF
      ENDDO
   ENDIF
ENDDO

! IF YOU GET TO THIS STAGE THEN SOMETHING IS WRONG
WRITE(LUPRI,*)'THE REQUESTEDBATCHINDEX GIVEN AS INPUT WAS NOT FOUND'
WRITE(LUPRI,*)'REQUESTEDBATCHINDEX GIVEN AS INPUT',RequestedBatchIndex
WRITE(LUPRI,*)'number of batches',nbatches
CALL LSQUIT('THE REQUESTEDBATCHINDEX GIVEN AS INPUT WAS NOT FOUND',lupri)

END SUBROUTINE BUILD_SINGLE_SHELLBATCH_AO

!> \brief builds an AOitem from a single AObatch
!> \author T. Kjaergaard
!> \date 2008
!>
!> builds an AOitem from a single AObatch
!>
SUBROUTINE build_aoitem_from_aobatch(LUPRI,AOBAT,MOLECULE,IATOM,AO)
implicit none
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the AObatch from which to build the AOitem
TYPE(AOBATCH) :: AOBAT
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> the atom that the AObatch bnelongs to
INTEGER  :: IATOM
!> the AOitem to be build
TYPE(AOITEM) :: AO
!
INTEGER       :: iangmom,nprim,ncont,nangmom
INTEGER       :: nOrbital,nprimOrb

AO%nbatches = 1
nangmom = AOBAT%nAngmom 
NULLIFY(AO%BATCH)
ALLOCATE(AO%BATCH(1))
AO%BATCH(1)%TYPE_Empty       = AOBAT%TYPE_Empty
AO%BATCH(1)%TYPE_Hermite     = AOBAT%TYPE_Hermite
AO%BATCH(1)%TYPE_Cartesian   = AOBAT%TYPE_Cartesian
AO%BATCH(1)%TYPE_Nucleus     = AOBAT%TYPE_Nucleus
AO%BATCH(1)%spherical        = AOBAT%spherical
AO%BATCH(1)%atom             = 1
AO%BATCH(1)%batch            = 1
AO%BATCH(1)%CENTER(1)        = MOLECULE%ATOM(iATOM)%CENTER(1)
AO%BATCH(1)%CENTER(2)        = MOLECULE%ATOM(iATOM)%CENTER(2)
AO%BATCH(1)%CENTER(3)        = MOLECULE%ATOM(iATOM)%CENTER(3)
AO%BATCH(1)%nPrimitives      = AOBAT%nPrimitives
AO%BATCH(1)%maxContracted    = AOBAT%maxContracted
AO%BATCH(1)%maxAngmom        = AOBAT%maxAngmom

AO%nExp =1
NULLIFY(AO%Exponents)
ALLOCATE(AO%Exponents(1))
nPrim = AOBAT%nPrimitives
CALL lsmat_dense_init(AO%Exponents(1),nPrim,1)
call dcopy(nPrim,AOBAT%pExponents%elms,1,&
        &AO%Exponents(1)%elms,1)
NULLIFY(AO%BATCH(1)%pExponents)
AO%BATCH(1)%pExponents  =>  AO%Exponents(1)

AO%BATCH(1)%nAngmom          = AOBAT%nAngmom     
AO%BATCH(1)%extent           = AOBAT%extent
AO%BATCH(1)%ANGMOM           = AOBAT%ANGMOM
AO%BATCH(1)%nContracted      = AOBAT%nContracted
nOrbital = 1
DO iangmom = 1,AOBAT%nAngmom     
   AO%BATCH(1)%startOrbital(iangmom) = nOrbital
   nOrbital=nOrbital+AOBAT%nContracted(iangmom)*AOBAT%nOrbComp(iangmom)
ENDDO
AO%BATCH(1)%startprimOrbital(1) = 1
nprimOrb = 1
DO iangmom = 1,AOBAT%nAngmom     
   AO%BATCH(1)%startprimOrbital(iangmom) = nprimOrb
   nprimOrb = nprim+AOBAT%nPrimitives*AOBAT%nOrbComp(iangmom)
ENDDO
DO iangmom = 1,AOBAT%nAngmom     
   AO%BATCH(1)%CCindex(iangmom)  = iangmom 
ENDDO
AO%BATCH(1)%nOrbComp         = AOBAT%nOrbComp
AO%BATCH(1)%nPrimOrbComp     = AOBAT%nPrimOrbComp
AO%BATCH(1)%nOrbitals        = AOBAT%nOrbitals

AO%nCC = nangmom
NULLIFY(AO%CC)
ALLOCATE(AO%CC(nangmom))
nPrim = AOBAT%nPrimitives
DO iangmom = 1,nangmom
   nCont = AOBAT%pCC(iangmom)%p%ncol
   CALL lsmat_dense_init(AO%CC(iangmom),nPrim,nCont)
   call dcopy(nPrim*nCont,AOBAT%pCC(iangmom)%p%elms,1,AO%CC(iangmom)%elms,1)
   NULLIFY(AO%BATCH(1)%pCC(iangmom)%p)
   AO%BATCH(1)%pCC(iangmom)%p => AO%CC(iangmom)
ENDDO

END SUBROUTINE BUILD_AOITEM_FROM_AOBATCH

!> \brief set the AO extent
!> \author S. Reine
!> \date 2010
SUBROUTINE SET_AO_EXTENT(AO,SCHEME)
implicit none
!> the AOitem to be build
TYPE(AOITEM)      :: AO
!> contains all info about the integralscheme requested (thresholds)
TYPE(LSINTSCHEME) :: SCHEME
!
Integer :: iBatch

DO iBatch=1,AO%nBatches
  AO%BATCH(iBatch)%extent = getAObatchExtent(AO%BATCH(iBatch),SCHEME)
ENDDO

END SUBROUTINE SET_AO_EXTENT

!> \brief calculate the AObatch extent
!> \author S. Reine
!> \date 2010
REAL(REALK) FUNCTION getAObatchExtent(AOB,SCHEME)
implicit none
!> the AObatch from which to calc extent
TYPE(AOBATCH)     :: AOB
!> contains all info about the integralscheme requested (thresholds)
TYPE(LSINTSCHEME) :: SCHEME
!
Integer     :: iAngmom,angmom,nPrim,nCont
Real(realk) :: extent

extent = 0.d0
nPrim = AOB%nPrimitives

DO iAngmom=1,AOB%nAngmom
  angmom = AOB%Angmom(iAngmom)
  nCont  = AOB%pCC(iAngmom)%p%ncol
  extent = max(extent,getExtent(AOB%pExponents%elms,AOB%pCC(iAngmom)%p%elms,&
     &           nPrim,nCont,angmom,SCHEME%OD_SCREEN,SCHEME%THRESHOLD*SCHEME%OD_THRESHOLD))
ENDDO

getAObatchExtent = extent

END FUNCTION getAObatchExtent

!> \brief calculate the extent of a single angular momentum for a single AObatch
!> \author S. Reine documented by T. Kjaergaard
!> \date 2010
REAL(REALK) FUNCTION getExtent(Exponents,CC,nPrim,nCont,angmom,Screen,Threshold)
implicit none
!> the exponents 
Real(realk) :: Exponents(nPrim)
!> the contractioncoefficients 
Real(realk) :: CC(nPrim,nCont)
!> the number of primitive orbitals
Integer     :: nPrim
!> the number of contractedorbitals
Integer     :: nCont
!> the angular momentum
Integer     :: angmom
!> if Overlap distribution screening ODscreening is used in this calc
Logical     :: Screen
!> the threshold for the Overlap distribution screening
Real(realk) :: Threshold
!
Integer     :: iPrim,iCont
Real(realk) :: extent,maxContraction,extent2,r2,fun,funD,rold,rnew

extent = 0.d0
IF (Screen) THEN
  extent2 = 0.d0
  DO iPrim=1,nPrim
    maxContraction = 0.d0
    DO iCont=1,nCont
      maxContraction = max(maxContraction,abs(CC(iPrim,iCont)))
    ENDDO
    IF ((Threshold.LT.tiny(1d0)).OR.(nCont.LT.1).OR.(Exponents(iPrim).LT.tiny(1d0))) THEN
      WRITE(*,*) 'Error in FUNCTION getExtent', Threshold, nCont,Exponents(iPrim)
      CALL LSQUIT('Error in FUNCTION getExtent',-1)
    ENDIF
!    r2 = (-log(Threshold)+log(maxContraction)+log(real(nCont)))/Exponents(iPrim)
    r2 = (-log(Threshold)+log(maxContraction))/Exponents(iPrim)
!   Take the above expression to be a constant A. We should then solve 
!      r^2 = A + l/a*ln(r)       (i)
!   to get a proper extent. We instead take A to be an estimate of r^2
!   and make the correction r^2 = A + l/a*ln(sqrt(A)). The equation (i)
!   can of course instead be solved iteratively.
    IF (r2.GT.0.d0) THEN
       IF(angmom .GT. 0)THEN
          Rold = sqrt(r2 + angmom*log(sqrt(r2))/Exponents(iPrim))
          fun = ABS(maxContraction*(rold**angmom)*exp(-Exponents(iPrim)*Rold**2))
          IF(threshold .LE. fun)THEN
             DO
                funD = ABS( maxContraction*angmom*(rold**(angmom-1))*&
                     &exp(-Exponents(iPrim)*Rold**2)+maxContraction*(rold**angmom)&
                     &*(-2*Exponents(iPrim)*Rold)*exp(-Exponents(iPrim)*Rold**2))
                Rnew = ROLD - (Threshold-fun)/funD
                IF(ABS(Rnew-ROLD) .LE. 1.0D-24)THEN
                   write(*,*)'calculation stalled'
                   EXIT
                ENDIF
                fun = ABS(maxContraction*(Rnew)**angmom*exp(-Exponents(iPrim)*Rnew**2))
                Rold=Rnew
                IF(ABS(fun-Threshold).LE. Threshold*1.0D-11)EXIT
                IF(Rold.NE.Rold)THEN
                   write(*,*)'found NaN in aobatch iteratvie extent'
                   EXIT
                ENDIF
             ENDDO
          ENDIF
          r2=Rold*Rold
       ENDIF
    ELSE
    ENDIF
    extent2 = max(extent2,r2)
  ENDDO
  IF (extent2.LT.0.d0) THEN
      WRITE(*,*) 'Negative squared distance in FUNCTION getExtent',extent2
      CALL LSQUIT('Negative squared distance in FUNCTION getExtent',-1)
  ENDIF
  extent = sqrt(extent2)
ENDIF
getExtent = extent

END FUNCTION getExtent

!> \brief determines if this segment of exponents are identical to earlier stored exponentes - used for family basis set
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE SHARES_EXPONENTS(LUPRI,IPRINT,AO,UNIQUEEXPONENTS,BASISINFO,type,B,J,UNIQUE,&
                                &UNCONTRACTED,MODELEXP,nMODELEXP)
implicit none
      !> the logical unit number for the output file
      INTEGER                   :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in)           :: IPRINT
!> The AOitem that contains already processed exponents
TYPE(AOITEM)              :: AO
!> the number of unique exponents saved in the AO%Exponents 
INTEGER                   :: UNIQUEEXPONENTS
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASISINFO
!> Atomtype for the current atom being processed
INTEGER                   :: type
!> Angular moment for the current atom being processed
INTEGER                   :: B
!> The segment for the current atom being processed
INTEGER                   :: J
!> Output: either 0 if the exponents are unique or the index of the place in AO%Exponents where these exponents are already stored.
INTEGER :: UNIQUE
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!> the number of differents segments of exponents saved in MODELEXP
INTEGER :: nMODELEXP
!> contains the index of where the differents segments of exponents are placed in AO%Exponents(I)
INTEGER :: MODELEXP(nMODELEXP)
!
INTEGER                   :: I,K,nPrim
real(realk)               :: SUM
REAL(REALK),PARAMETER     :: NULL=0.d0
IF(UNIQUEEXPONENTS .GT. 0) THEN
   IF(UNCONTRACTED)THEN
      UNIQUE=0
      DO I=1,UNIQUEEXPONENTS
         IF(ABS(BASISINFO%ATOMTYPE(type)%SHELL(B)%segment(J)%Exponents(1)-AO%Exponents(I)%elms(1)) < ExpThr)THEN
            UNIQUE=MODELEXP(I)       
         ENDIF
      ENDDO
   ELSE !DEFAULT
      UNIQUE=0
      nPrim=BASISINFO%ATOMTYPE(type)%SHELL(B)%segment(J)%nrow
      DO I=1,UNIQUEEXPONENTS
         IF(AO%Exponents(I)%nrow .EQ. nPrim)THEN
            SUM=NULL
            DO K=1,nPrim
               SUM=SUM+ABS(BASISINFO%ATOMTYPE(type)%SHELL(B)%segment(J)%Exponents(K)&
                    &-AO%Exponents(I)%elms(K))
            ENDDO
            IF (SUM<ExpThr) THEN
               UNIQUE=MODELEXP(I)       
            ENDIF
         ENDIF
      ENDDO
   ENDIF
ELSE
UNIQUE=0
ENDIF
IF(IPRINT .GT. 2)THEN
   IF(UNIQUE/=0) WRITE(LUPRI,*)'IT WAS DETERMINED THAT THIS SEGMENT SHARES'
   IF(UNIQUE/=0) WRITE(LUPRI,*)'EXPONENTS WITH SEGMENT:',UNIQUE
ENDIF

END SUBROUTINE SHARES_EXPONENTS

!> \brief add batch to aoitem
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE ADD_BATCH(SCHEME,MOLECULE,BASISINFO,AO,AOorganize,iATOM,type,B,J,&
     &lupri,iprint,nPrim,nCont,CARTESIANTYPE,nbatches,UNCONTRACTED)
use molecule_type
implicit none
!> contains all info about the integralscheme requested (thresholds,use cartesian,..)
TYPE(LSINTSCHEME)         :: SCHEME
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASISINFO 
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> contains info about the full AOorganize
TYPE(AOorganizer)         :: AOorganize
!> the atomtype for the current atom
INTEGER                   :: type
!> the angular momentum for the current atom
INTEGER                   :: B
!> the segment for the current atom
INTEGER                   :: J
!> the logical unit number for the output file
INTEGER                   :: lupri
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: iprint
!> number of contracted functions
INTEGER,intent(out)       :: nCont
!> number of primitive functions
INTEGER,intent(out)       :: nPrim
!> label (Hermite,Cartesian or Default)
Character(len=80)         :: CartesianTYPE
!> nbatches = old number of batches + new number of batches
INTEGER                   :: nbatches
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!
INTEGER                   :: iATOM,I

IF(UNCONTRACTED)THEN
   !WE ADD ONE BATCH FOR EACH PRIMITIVES
   nPrim=BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
   nCont=nPrim
   DO I=1,nPrim
      nbatches = AO%nbatches + I
      AOorganize%ORG(nbatches)%type=type
      AOorganize%ORG(nbatches)%angmoms=1
      AOorganize%ORG(nbatches)%ANGMOM(1)=B-1
      AOorganize%ORG(nbatches)%Segment(1)=J
      AOorganize%ORG(nbatches)%atom=iATOM
      AO%BATCH(nbatches)%nAngmom = 1
      AO%BATCH(nbatches)%type_cartesian  = CartesianTYPE(1:4).EQ.'Cart'
      AO%BATCH(nbatches)%type_hermite  = CartesianTYPE(1:4).EQ.'Herm'
      AO%BATCH(nbatches)%type_Empty  = .FALSE. 
      AO%BATCH(nbatches)%type_Nucleus  = .FALSE.
      AO%BATCH(nbatches)%spherical= SCHEME%DoSpherical
      AO%BATCH(nbatches)%Angmom(1)= B - 1
      AO%BATCH(nbatches)%maxAngmom= B - 1
      AO%BATCH(nbatches)%CENTER(1)=MOLECULE%ATOM(iATOM)%CENTER(1)
      AO%BATCH(nbatches)%CENTER(2)=MOLECULE%ATOM(iATOM)%CENTER(2)
      AO%BATCH(nbatches)%CENTER(3)=MOLECULE%ATOM(iATOM)%CENTER(3)
   ENDDO
   AO%nbatches = nbatches
   AOorganize%nbatches = nbatches
ELSE !DEFAULT
   nbatches = AO%nbatches + 1
   AO%nbatches = nbatches 
   AOorganize%nbatches = nbatches
   AOorganize%ORG(nbatches)%type=type
   AOorganize%ORG(nbatches)%angmoms=1
   AOorganize%ORG(nbatches)%ANGMOM(1)=B-1
   AOorganize%ORG(nbatches)%Segment(1)=J
   AOorganize%ORG(nbatches)%atom=iATOM
   AO%BATCH(nbatches)%nAngmom = 1
   AO%BATCH(nbatches)%type_cartesian  = CartesianTYPE(1:4).EQ.'Cart'
   AO%BATCH(nbatches)%type_hermite  = CartesianTYPE(1:4).EQ.'Herm'
   AO%BATCH(nbatches)%type_Empty  = .FALSE. 
   AO%BATCH(nbatches)%type_Nucleus  = .FALSE.
   AO%BATCH(nbatches)%spherical= SCHEME%DoSpherical
   AO%BATCH(nbatches)%Angmom(1)= B - 1
   AO%BATCH(nbatches)%maxAngmom= B - 1
   AO%BATCH(nbatches)%CENTER(1)=MOLECULE%ATOM(iATOM)%CENTER(1)
   AO%BATCH(nbatches)%CENTER(2)=MOLECULE%ATOM(iATOM)%CENTER(2)
   AO%BATCH(nbatches)%CENTER(3)=MOLECULE%ATOM(iATOM)%CENTER(3)
   nPrim=BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
   nCont=BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
ENDIF

END SUBROUTINE ADD_BATCH

!> \brief Does different initializations before get_density is called.
!> \author T. Kjaergaard
!> \date 2008
!>
!> Add a segment of exponents and the corresponding contraction coefficients  
!> to the unique exponents and contraction coefficients
!>
SUBROUTINE ADD_SEGMENT_AND_CC(AO,AOmodel,AOorganize,BASISINFO,nbatches,nPrim,nCont,&
     & type,B,J,lupri,iprint,UNIQUEMATRIXES,UNIQUEEXPONENTS,UNCONTRACTED,&
     & INTNRM,MODELEXP,nMODELEXP)
IMPLICIT NONE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the model AOitem for the atom 
TYPE(AOITEM)              :: AOmodel
!> the AOorganizer that contain the unique exponents
TYPE(AOorganizer)         :: AOorganize 
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASISINFO 
!> the number of batches already in the full AO
INTEGER                   :: nbatches
!> the number of primitive functions
INTEGER                   :: nPrim
!> the number of contracted functions
INTEGER                   :: nCont
!> the atomtype for the current atom
INTEGER                   :: type
!> the angular momentum for the current atom
INTEGER                   :: B
!> the segment for the current atom
INTEGER                   :: J
!> the logical unit number for the output file
INTEGER                   :: lupri
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: iprint
!> the number of unique contractioncoefficientmatrices saved
INTEGER                   :: UNIQUEMATRIXES
!> the number of unique exponents saved
INTEGER                   :: UNIQUEEXPONENTS
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!> if INTERMIDIATE NORMALIZATION shoule be used, USED FOR SCREENING
LOGICAL                   :: INTNRM 
!> the number of differents segments of exponents saved in MODELEXP
INTEGER :: nMODELEXP
!> contains the index of where the differents segments of exponents are placed in AO%Exponents(I)
INTEGER :: MODELEXP(nMODELEXP)
INTEGER                   :: nbat,I,icont,iprim
REAL(REALK)               :: PI,EXPO,PIPPI,maxcont,element
PI=3.14159265358979323846D0

IF(UNCONTRACTED)THEN
 PIPPI = (0.5D0/PI)**(0.75D0)
 nbat = nbatches-nPrim
 IF(INTNRM)THEN
   !DETERMINE MAX CONTRACTION ELEMENT
   MAXcont=0.d0
   DO iprim=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
     DO iCont=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
       element=ABS(BASISINFO%ATOMTYPE(type)&
       &%SHELL(B)%SEGMENT(J)%elms(iprim+(iCont-1)*nPrim))
       MAXcont=MAX(MAXcont,element)
     ENDDO
   ENDDO
 ENDIF

 DO I=1,nPrim
    nbatches = nbat+I
    UNIQUEMATRIXES=UNIQUEMATRIXES+1
    UNIQUEEXPONENTS=UNIQUEEXPONENTS+1
    AOmodel%nExp=AOmodel%nExp+1
    MODELEXP(AOmodel%nExp)=UNIQUEEXPONENTS
    AOorganize%ORG(nbatches)%exponentindex=UNIQUEEXPONENTS
    AOorganize%ORG(nbatches)%CC(1)=UNIQUEMATRIXES
    CALL lsmat_dense_init(AO%CC(UNIQUEMATRIXES),1,1)
    CALL lsmat_dense_init(AO%Exponents(UNIQUEEXPONENTS),1,1)
    CALL lsmat_dense_init(AOmodel%Exponents(AOmodel%nExp),1,1)
    
    Expo = BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%Exponents(I)

    AO%Exponents(UNIQUEEXPONENTS)%elms(1)= Expo
    AOmodel%Exponents(AOmodel%nExp)%elms(1)= Expo

    IF(IPRINT .GT. 2)THEN
       WRITE(LUPRI,*)'THE NEW EXPONENTS',UNIQUEEXPONENTS
       call lsmat_dense_print(AO%Exponents(UNIQUEEXPONENTS),&
            &1,1,1,1,lupri)
    ENDIF
    
    IF(INTNRM)THEN
       AO%CC(UNIQUEMATRIXES)%elms(1)= MAXcont*&
            & (4*Expo)**(0.5D0*B+0.25d0)*PIPPI 
    ELSE
       AO%CC(UNIQUEMATRIXES)%elms(1)= &
            & (4*Expo)**(0.5D0*B+0.25d0)*PIPPI 
    ENDIF
    
    IF(IPRINT .GT. 2)THEN
       WRITE(LUPRI,*)'THE NEW CC',UNIQUEMATRIXES
       call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
            &1,1,1,1,lupri)
    ENDIF
 ENDDO
ELSE !DEFAULT
   UNIQUEMATRIXES=UNIQUEMATRIXES+1
   UNIQUEEXPONENTS=UNIQUEEXPONENTS+1
   AOmodel%nExp=AOmodel%nExp+1
   MODELEXP(AOmodel%nExp)=UNIQUEEXPONENTS
   AOorganize%ORG(nbatches)%exponentindex=UNIQUEEXPONENTS
   AOorganize%ORG(nbatches)%CC(1)=UNIQUEMATRIXES
   CALL lsmat_dense_init(AO%CC(UNIQUEMATRIXES),nPrim,nCont)
   CALL lsmat_dense_init(AO%Exponents(UNIQUEEXPONENTS),nPrim,1)
   CALL lsmat_dense_init(AOmodel%Exponents(AOmodel%nExp),nPrim,1)
   call dcopy (nPrim*nCont,BASISINFO%ATOMTYPE(type)%SHELL(B)%&
        &SEGMENT(J)%elms,1,AO%CC(UNIQUEMATRIXES)%elms,1)
   
   IF(IPRINT .GT. 2)THEN
      WRITE(LUPRI,*)'THE NEW CC',UNIQUEMATRIXES
      call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
           &1,nPrim,1,nCont,lupri)
   ENDIF
   
   call dcopy (nPrim,BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%Exponents,1,&
        &AO%Exponents(UNIQUEEXPONENTS)%elms,1)
   call dcopy (nPrim,BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%Exponents,1,&
        &AOmodel%Exponents(AOmodel%nExp)%elms,1)
   
   IF(IPRINT .GT. 2)THEN
      WRITE(LUPRI,*)'THE NEW EXPONENTS',UNIQUEEXPONENTS
      call lsmat_dense_print(AO%Exponents(UNIQUEEXPONENTS),&
           &1,nPrim,1,1,lupri)
   ENDIF
ENDIF

END SUBROuTINE ADD_SEGMENT_AND_CC

!> \brief Direct pointers to point to the unique exponents and contraction coefficients. 
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE DIRECT_POINTERS(AO,AOmodel,AOorganize,nCont,nPrim,nbatches,B,lupri,iprint,orbitalIndex,&
     &primorbitalindex,UNCONTRACTED)
IMPLICIT NONE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the AOitem for the current atom
TYPE(AOITEM)              :: AOmodel
!> the AOorganizer
TYPE(AOorganizer)         :: AOorganize 
!> the number of contracted functions
INTEGER                   :: nCont
!> the number of primitive functions
INTEGER                   :: nPrim
!> the number of batches already in the full AO
INTEGER                   :: nbatches
!> the angular momentum for the current atom
INTEGER                   :: B
!> the logical unit number for the output file
INTEGER                   :: lupri
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: iprint
!> the orbital index
INTEGER                   :: orbitalIndex
!> the primitive orbital index
INTEGER                   :: primorbitalindex
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
INTEGER                   :: nOrbitals,I,nbat,A

IF(UNCONTRACTED)THEN
   nbat = nbatches-nPrim
   DO I=1,nPrim
      nbatches = nbat+I
      A=AOmodel%BATCH(nbatches)%nangmom
      IF(A .GT. maxAOangmom)THEN
         WRITE(LUPRI,*)'You have at least one atom with many different orbitals&
              & - more than just the normal s, p, d, f, g, h, i, j, k and l&
              & orbitals so you need to increase the maxAOangmom parameter&
              & in TYPE-DEF.f90 file and recompile. It should be even to or&
              & greater than the number of different orbitals you need.&
              & Thomas Kjaergaard'
         CALL LSQUIT('Increase maxAOangmom in TYPE-DEF.f90',lupri)
      ENDIF
      AOmodel%BATCH(nbatches)%nContracted(A) = 1
      AOmodel%BATCH(nbatches)%maxContracted  = 1
      AOmodel%BATCH(nbatches)%nPrimitives=1
   
      IF(IPRINT .GT. 2) WRITE(LUPRI,*)'CC pointer to ',AOorganize%ORG(nbatches)%CC(1)
   
      AOmodel%BATCH(nbatches)%pCC(A)%p=>&
        &AO%CC(AOorganize%ORG(nbatches)%CC(1))
      AOmodel%BATCH(nbatches)%CCindex(A)=AOorganize%ORG(nbatches)%CC(1)
      IF(IPRINT .GT. 30) THEN
         WRITE(LUPRI,*)'Which is'
         CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbatches)%pCC(A)%p,1,1,1,1,lupri)
      ENDIF
   
      IF(IPRINT .GT. 2)WRITE(LUPRI,*)'Exp pointer to ',AOorganize%ORG(nbatches)%exponentindex

      AOmodel%BATCH(nbatches)%pExponents=>&
           &AO%Exponents(AOorganize%ORG(nbatches)%exponentindex)

      IF(IPRINT .GT. 30) THEN
         WRITE(LUPRI,*)'Which is'
!         CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbatches)%pExponents,1,1,1,1,lupri)
         CALL LSMAT_DENSE_PRINT(AO%Exponents(AOorganize%ORG(nbatches)%exponentindex),1,1,1,1,lupri)
      ENDIF
   
      ! Orbital index (used in for example density-matrix, Fock/KS-matrix,
      ! fitting coefficient etc.)
      AOmodel%BATCH(nbatches)%startOrbital(A) = orbitalIndex
      AOmodel%BATCH(nbatches)%startprimOrbital(A) = primorbitalIndex
      AOmodel%BATCH(nbatches)%nPrimOrbComp(A) = B*(B+1)/2
      IF (AOmodel%BATCH(nbatches)%spherical) THEN
         AOmodel%BATCH(nbatches)%nOrbComp(A) = 2*B-1
      ELSE
         AOmodel%BATCH(nbatches)%nOrbComp(A) = B*(B+1)/2
      ENDIF
      nOrbitals = AOmodel%BATCH(nbatches)%nOrbComp(A)
      orbitalIndex = orbitalIndex + nOrbitals
      primorbitalIndex = primorbitalIndex + nOrbitals
      AOmodel%BATCH(nbatches)%nOrbitals(A)    = nOrbitals
   ENDDO
ELSE
   A=AOmodel%BATCH(nbatches)%nangmom
   IF(A .GT. maxAOangmom)THEN
      WRITE(LUPRI,*)'You have at least one atom with many different orbitals&
           & - more than just the normal s, p, d, f, g, h, i, j, k and l&
           & orbitals so you need to increase the maxAOangmom parameter&
           & in TYPE-DEF.f90 file and recompile. It should be even to or&
           & greater than the number of different orbitals you need.&
           & Thomas Kjaergaard'
      CALL LSQUIT('Increase maxAOangmom in TYPE-DEF.f90',lupri)
   ENDIF
   AOmodel%BATCH(nbatches)%nContracted(A) = nCont
   AOmodel%BATCH(nbatches)%maxContracted  = nCont
   AOmodel%BATCH(nbatches)%nPrimitives=nPrim
   
   IF(IPRINT .GT. 2) WRITE(LUPRI,*)'CC pointer to ',AOorganize%ORG(nbatches)%CC(1)
   
   AOmodel%BATCH(nbatches)%pCC(A)%p=>&
        &AO%CC(AOorganize%ORG(nbatches)%CC(1))
   AOmodel%BATCH(nbatches)%CCindex(A)=AOorganize%ORG(nbatches)%CC(1)
   
   IF(IPRINT .GT. 30) THEN
      WRITE(LUPRI,*)'Which is'
      CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbatches)%pCC(A)%p,1,nPrim,1,nCont,lupri)
   ENDIF
   
   
   ! Orbital index (used in for example density-matrix, Fock/KS-matrix,
   ! fitting coefficient etc.)
   AOmodel%BATCH(nbatches)%startOrbital(A) = orbitalIndex
   AOmodel%BATCH(nbatches)%startprimOrbital(A) = primorbitalIndex
   AOmodel%BATCH(nbatches)%nPrimOrbComp(A) = B*(B+1)/2
   IF (AOmodel%BATCH(nbatches)%spherical) THEN
      AOmodel%BATCH(nbatches)%nOrbComp(A) = 2*B-1
   ELSE
      AOmodel%BATCH(nbatches)%nOrbComp(A) = B*(B+1)/2
   ENDIF
   nOrbitals = AOmodel%BATCH(nbatches)%nOrbComp(A)*nCont
   orbitalIndex = orbitalIndex + nOrbitals
   primorbitalIndex = primorbitalIndex &
        &+ AOmodel%BATCH(nbatches)%nOrbComp(A)*nprim
   AOmodel%BATCH(nbatches)%nOrbitals(A)    = nOrbitals
   IF(IPRINT .GT. 2)WRITE(LUPRI,*)'Exp pointer to ',AOorganize%ORG(nbatches)%exponentindex
   AOmodel%BATCH(nbatches)%pExponents=>&
        &AO%Exponents(AOorganize%ORG(nbatches)%exponentindex)
   
   IF(IPRINT .GT. 30) THEN
      WRITE(LUPRI,*)'Which is'
      CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbatches)%pExponents,1,nPrim,1,1,lupri)
   ENDIF
ENDIF
END SUBROUTINE DIRECT_POINTERS

!> \brief Determine the number of batches for the atom index atom 
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE DETERMINE_nbathces(AOorganize,atom,UNIQUE,nbatches)
implicit none
TYPE(AOorganizer) :: AOorganize
INTEGER           :: atom,nbatches,I,UNIQUE
nbatches=0
DO I=AOorganize%nbatches,1,-1
   IF(AOorganize%ORG(I)%atom==atom)THEN
      IF(AOorganize%ORG(I)%exponentindex==UNIQUE)THEN
         nbatches=I
         RETURN
      ENDIF
   ENDIF
ENDDO

END SUBROUTINE DETERMINE_nbathces

!> \brief adding extra angular moments to the already existing BATCH(family basisset)
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE EXPAND_BATCH(AO,AOorganize,BASISINFO,type,B,J,nbat,A,lupri,&
     &UNCONTRACTED)
IMPLICIT NONE
INTEGER                   :: A,lupri,I,nprim
INTEGER,intent(in)        :: B,nbat,J,type
TYPE(BASISSETINFO)        :: BASISINFO 
TYPE(AOITEM)              :: AO
TYPE(AOorganizer)         :: AOorganize 
LOGICAL                   :: UNCONTRACTED

IF(UNCONTRACTED)THEN
   nPrim=BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%nrow
   DO I=0,nPrim-1
      A=AO%BATCH(nbat+I)%nAngmom+1
      IF(A .GT. maxAOangmom)THEN
        WRITE(LUPRI,*)'You have at least one atom with many different orbitals&
             & - more than just the normal s, p, d, f, g, h, i, j, k and l&
             & orbitals so you need to increase the maxAOangmom parameter&
             & in TYPE-DEF.f90 file and recompile. It should be even to or&
             & greater than the number of different orbitals you need.&
             & Thomas Kjaergaard'
        CALL LSQUIT('Increase maxAOangmom in TYPE-DEF.f90',lupri)
      ENDIF
      AO%BATCH(nbat+I)%nAngmom = A
      AO%BATCH(nbat+I)%Angmom(A)   = B - 1
      AO%BATCH(nbat+I)%maxAngmom = B - 1!MAX(AO%BATCH(UNIQUE)%maxAngmom,B)
      AOorganize%ORG(nbat+I)%angmoms= A
      AOorganize%ORG(nbat+I)%ANGMOM(A)=B-1
   ENDDO
ELSE
   A=AO%BATCH(nbat)%nAngmom+1
   IF(A .GT. maxAOangmom)THEN
      WRITE(LUPRI,*)'You have at least one atom with many different orbitals&
           & - more than just the normal s, p, d, f, g, h, i, j, k and l&
           & orbitals so you need to increase the maxAOangmom parameter&
           & in TYPE-DEF.f90 file and recompile. It should be even to or&
           & greater than the number of different orbitals you need.&
           & Thomas Kjaergaard'
      CALL LSQUIT('Increase maxAOangmom in TYPE-DEF.f90',lupri)
   ENDIF
   AO%BATCH(nbat)%nAngmom = A
   AO%BATCH(nbat)%Angmom(A)   = B - 1
   AO%BATCH(nbat)%maxAngmom = B - 1!MAX(AO%BATCH(UNIQUE)%maxAngmom,B)
   AOorganize%ORG(nbat)%angmoms= A
   AOorganize%ORG(nbat)%ANGMOM(A)=B-1
ENDIF

END SUBROUTINE EXPAND_BATCH

!> \brief add a contraction coefficient matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ADD_CC(AO,AOorganize,nbatches,nPrim,nCont,lupri,iprint,&
     &BASISINFO,type,B,J,A,UNIQUEMATRIXES,UNCONTRACTED,INTNRM)
IMPLICIT NONE
TYPE(AOITEM)              :: AO
TYPE(AOorganizer)         :: AOorganize
TYPE(BASISSETINFO)        :: BASISINFO 
INTEGER                   :: UNIQUEMATRIXES,nbatches,nPrim
INTEGER                   :: lupri,nCont,iprint
INTEGER                   :: type,B,J,A,I,icont,iprim
REAL(REALK)               :: PI,EXP,PIPPI,maxcont
LOGICAL                   :: UNCONTRACTED,INTNRM
PI=3.14159265358979323846D0

   nPrim=BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%nrow
   nCont=BASISINFO%ATOMTYPE(type)%SHELL(B)&
        &%SEGMENT(J)%ncol

IF(UNCONTRACTED)THEN  
   PIPPI = (0.5D0/PI)**(0.75D0)
   IF(INTNRM)THEN
      !DETERMINE MAX CONTRACTION ELEMENT
      MAXcont=0.d0
      DO iprim=1,nPrim
         DO iCont=1,nCont
            MAXcont=MAX(MAXcont,ABS(BASISINFO%ATOMTYPE(type)&
                 &%SHELL(B)%SEGMENT(J)%elms(iprim+(iCont-1)*nPrim)))
!            WRITE(LUPRI,*)'MAXcont',MAXcont

         ENDDO
      ENDDO
!      WRITE(LUPRI,*)'FINAL MAXcont',MAXcont
   ENDIF
   DO I=0,nPrim-1
      UNIQUEMATRIXES=UNIQUEMATRIXES+1
      AOorganize%ORG(nbatches+I)%CC(A)=UNIQUEMATRIXES
      CALL lsmat_dense_init(AO%CC(UNIQUEMATRIXES),1,1)
      Exp = BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%Exponents(I+1)
      IF(INTNRM)THEN
         AO%CC(UNIQUEMATRIXES)%elms(1)= MAXcont* &
              & (4*Exp)**(0.5D0*B+0.25d0)*PIPPI 
!         WRITE(LUPRI,*)'CCMAT',MAXcont* &
!              & (4*Exp)**(0.5D0*B+0.25d0)*PIPPI
!         WRITE(LUPRI,*)'CCMAT-MAXCONT',&
!              & (4*Exp)**(0.5D0*B+0.25d0)*PIPPI

      ELSE
         AO%CC(UNIQUEMATRIXES)%elms(1)=&
              & (4*Exp)**(0.5D0*B+0.25d0)*PIPPI 
      ENDIF
      IF(IPRINT .GT. 2)THEN
         WRITE(LUPRI,*)'THE NEW CC',UNIQUEMATRIXES
         call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
              &1,1,1,1,lupri)
      ENDIF
   ENDDO
ELSE
   UNIQUEMATRIXES=UNIQUEMATRIXES+1
   AOorganize%ORG(nbatches)%CC(A)=UNIQUEMATRIXES
   CALL lsmat_dense_init(AO%CC(UNIQUEMATRIXES),nPrim,nCont)
   call dcopy (nPrim*nCont,BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%elms,1,&
        &AO%CC(UNIQUEMATRIXES)%elms,1)
   IF(IPRINT .GT. 2)THEN
      WRITE(LUPRI,*)'THE NEW CC',UNIQUEMATRIXES
      call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
           &1,nPrim,1,nCont,lupri)
   ENDIF
ENDIF

END SUBROuTINE ADD_CC

!> \brief direct contraction coefficient matrix pointers to the contraction coefficients in the AOorganize
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE DIRECT_CC_POINTERS(AO,AOmodel,AOorganize,nCont,nPrim,A,nbat,B,&
     &lupri,iprint,orbitalIndex,primorbitalIndex,UNCONTRACTED)
IMPLICIT NONE
INTEGER,intent(in)        :: A
TYPE(AOITEM)              :: AO,AOmodel
TYPE(AOorganizer)         :: AOorganize 
INTEGER                   :: nOrbitals,orbitalIndex,lupri,iprint,B
INTEGER                   :: nCont,nPrim,nbat,I,primorbitalindex
LOGICAL                   :: UNCONTRACTED
IF(UNCONTRACTED)THEN
   DO I=0,nPrim-1
      AOmodel%BATCH(nbat+I)%pCC(A)%p=>&
           &AO%CC(AOorganize%ORG(nbat+I)%CC(A))
      AOmodel%BATCH(nbat+I)%CCindex(A)=AOorganize%ORG(nbat+I)%CC(A)
      IF(IPRINT .GT. 2) WRITE(LUPRI,*)'CC pointer to        ',AOorganize%ORG(nbat)%CC(A)
      IF(IPRINT .GT. 5)THEN
         WRITE(LUPRI,*)'Exp already points to '
         CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbat+I)%pExponents,1,1,1,1,lupri)
      ENDIF
      AOmodel%BATCH(nbat+I)%nContracted(A) = 1
      AOmodel%BATCH(nbat+I)%maxContracted= 1
      AOmodel%BATCH(nbat+I)%startOrbital(A) = orbitalIndex
      AOmodel%BATCH(nbat+I)%startprimOrbital(A) = primorbitalIndex
      AOmodel%BATCH(nbat+I)%nPrimOrbComp(A) = B*(B+1)/2
      IF (AOmodel%BATCH(nbat+I)%spherical) THEN
         AOmodel%BATCH(nbat+I)%nOrbComp(A) = 2*B-1
      ELSE
         AOmodel%BATCH(nbat+I)%nOrbComp(A) = B*(B+1)/2
      ENDIF
      nOrbitals = AOmodel%BATCH(nbat+I)%nOrbComp(A)
      orbitalIndex = orbitalIndex + nOrbitals
      primorbitalIndex = primorbitalIndex + nOrbitals
   ENDDO
ELSE !DEFAULT
   AOmodel%BATCH(nbat)%pCC(A)%p=>&
        &AO%CC(AOorganize%ORG(nbat)%CC(A))
   AOmodel%BATCH(nbat)%CCindex(A)=AOorganize%ORG(nbat)%CC(A)

   IF(IPRINT .GT. 2) WRITE(LUPRI,*)'CC pointer to        ',AOorganize%ORG(nbat)%CC(A)
   IF(IPRINT .GT. 5)THEN
      WRITE(LUPRI,*)'Exp already points to '
      CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbat)%pExponents,1,&
           &AOmodel%BATCH(nbat)%pExponents%nrow,1,&
           &AOmodel%BATCH(nbat)%pExponents%ncol,lupri)
   ENDIF
   AOmodel%BATCH(nbat)%nContracted(A) = nCont
   AOmodel%BATCH(nbat)%maxContracted=MAX(AOmodel%BATCH(nbat)%maxContracted,nCont)
   AOmodel%BATCH(nbat)%startOrbital(A) = orbitalIndex
   AOmodel%BATCH(nbat)%startprimOrbital(A) = primorbitalIndex
   AOmodel%BATCH(nbat)%nPrimOrbComp(A) = B*(B+1)/2
   IF (AOmodel%BATCH(nbat)%spherical) THEN
      AOmodel%BATCH(nbat)%nOrbComp(A) = 2*B-1
   ELSE
      AOmodel%BATCH(nbat)%nOrbComp(A) = B*(B+1)/2
   ENDIF
   nOrbitals = AOmodel%BATCH(nbat)%nOrbComp(A)*nCont
   orbitalIndex = orbitalIndex + nOrbitals
   primorbitalIndex = primorbitalIndex&
        & + AOmodel%BATCH(nbat)%nOrbComp(A)*nprim
   AOmodel%BATCH(nbat)%nOrbitals(A)    = nOrbitals      
ENDIF

END SUBROUTINE DIRECT_CC_POINTERS

!> \brief in case of old atoms we copy AObatches from a modelAO to the AO
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE COPY_FROM_MODEL_AO(AOmodel,AO,iATOM,MOLECULE,orbitalIndex,primorbitalindex,lupri)
use molecule_type
IMPLICIT NONE
TYPE(MOLECULEINFO)        :: MOLECULE
TYPE(AOITEM)              :: AO,AOmodel
INTEGER                   :: I,J,bat,iATOM,orbitalIndex,lupri,SUM
INTEGER                   :: primorbitalindex,pSUM,norbitals

SUM=0
pSUM=0
AO%ATOMICnBatch(IATOM) = AOmodel%nbatches 
DO I=1,AOmodel%nbatches
   bat=AO%nbatches+1
   AO%nbatches = bat
   AO%BATCH(bat)%atom = IATOM
   AO%BATCH(bat)%batch = I
   AO%BATCH(bat)%type_empty = AOmodel%BATCH(I)%type_empty
   AO%BATCH(bat)%type_hermite = AOmodel%BATCH(I)%type_hermite
   AO%BATCH(bat)%type_cartesian = AOmodel%BATCH(I)%type_cartesian
   AO%BATCH(bat)%type_nucleus = AOmodel%BATCH(I)%type_nucleus
   AO%BATCH(bat)%spherical = AOmodel%BATCH(I)%spherical
   AO%BATCH(bat)%CENTER(1) = MOLECULE%ATOM(iATOM)%CENTER(1)
   AO%BATCH(bat)%CENTER(2) = MOLECULE%ATOM(iATOM)%CENTER(2)
   AO%BATCH(bat)%CENTER(3) = MOLECULE%ATOM(iATOM)%CENTER(3)
   AO%BATCH(bat)%nPrimitives = AOmodel%BATCH(I)%nPrimitives
   AO%BATCH(bat)%maxContracted = AOmodel%BATCH(I)%maxContracted
   AO%BATCH(bat)%maxAngmom = AOmodel%BATCH(I)%maxAngmom
   AO%BATCH(bat)%pExponents => AOmodel%BATCH(I)%pExponents
   AO%BATCH(bat)%nAngmom = AOmodel%BATCH(I)%nAngmom
   AO%BATCH(bat)%extent  = AOmodel%BATCH(I)%extent 
   DO J=1,AOmodel%BATCH(I)%nAngmom
      AO%BATCH(bat)%Angmom(J) = AOmodel%BATCH(I)%Angmom(J)
      AO%BATCH(bat)%nContracted(J) = AOmodel%BATCH(I)%nContracted(J)
      AO%BATCH(bat)%nOrbComp(J) = AOmodel%BATCH(I)%nOrbComp(J)
      AO%BATCH(bat)%nPrimOrbComp(J) = AOmodel%BATCH(I)%nPrimOrbComp(J)
      AO%BATCH(bat)%pCC(J)%p =>  AOmodel%BATCH(I)%pCC(J)%p
      AO%BATCH(bat)%CCindex(J) =  AOmodel%BATCH(I)%CCindex(J)

      AO%BATCH(bat)%startOrbital(J) = orbitalIndex &
           &+ AOmodel%BATCH(I)%startOrbital(J)  
      AO%BATCH(bat)%startprimOrbital(J) = primorbitalIndex &
           &+ AOmodel%BATCH(I)%startprimOrbital(J)  

      nOrbitals = AOmodel%BATCH(I)%nOrbComp(J)*AOmodel%BATCH(I)%nContracted(J)
      AO%BATCH(bat)%nOrbitals(J) = nOrbitals
      SUM=SUM+nOrbitals
      pSUM=pSUM+AOmodel%BATCH(I)%nOrbComp(J)*AOmodel%BATCH(I)%nPrimitives
   ENDDO
ENDDO
orbitalIndex = orbitalIndex + SUM
AO%ATOMICnORB(IATOM) = SUM
primorbitalIndex = primorbitalIndex + pSUM

END SUBROUTINE COPY_FROM_MODEL_AO

!> \brief free the AOitem structure
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE FREE_AOITEM(LUPRI,AO)
IMPLICIT NONE
TYPE(AOITEM)              :: AO
INTEGER                   :: I,J,LUPRI

DEALLOCATE(AO%ATOMICnORB)
NULLIFY(AO%ATOMICnORB)
DEALLOCATE(AO%ATOMICnBATCH)
NULLIFY(AO%ATOMICnBATCH)

DO I=1,AO%nCC
   CALL LSMAT_DENSE_FREE(AO%CC(I))
ENDDO
DO I=1,AO%nExp
   CALL LSMAT_DENSE_FREE(AO%Exponents(I))
ENDDO
DO I=1,AO%nbatches
   DO J=1,AO%BATCH(I)%nAngmom
      NULLIFY(AO%BATCH(I)%pCC(J)%p)
   ENDDO
   NULLIFY(AO%BATCH(I)%pExponents)
ENDDO
DEALLOCATE(AO%BATCH)
NULLIFY(AO%BATCH)
DEALLOCATE(AO%CC)
NULLIFY(AO%CC)
DEALLOCATE(AO%Exponents)
NULLIFY(AO%Exponents)

END SUBROUTINE FREE_AOITEM

!> \brief build basinf structure used in exchange-correlation calculations
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE BUILD_BASINF(LUPRI,IPRINT,BAS,SETTING,GRDONE)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER           :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER           :: IPRINT
!> the basinf structure to be build
TYPE(BASINF)      :: BAS
!> logical for grid. 1 if grid have been build and 0 if not 
INTEGER           :: GRDONE 
!> the setting structure containing info about molecule,basisset,...
TYPE(LSSETTING)   :: SETTING
!
INTEGER           :: KMAX !number of shells see type def of BASINF
INTEGER           :: NHTYP
INTEGER           :: nprim,TOTPRIM,nrow,ncol,N,OLDnrow
INTEGER           :: TYP,NR,USHELL,nshells,ushell2
INTEGER           :: natoms,type,SHELL,orbitalindex,L,I,K,norb,R,M,MXPRIM
INTEGER,pointer :: ATOMtype(:)
INTEGER,pointer :: MODELTYPES(:)
LOGICAL           :: NEWATOM
REAL(REALK)       :: FACL(10),R2,THLOG,EXP
!TYPE(CCMODEL),ALLOCATABLE  :: UniqueAtom(:)
REAL(REALK),PARAMETER     :: DMIN=1.D-13
INTEGER           :: irow,nrow2,nrow3,NRSIZE,ICHARGE
LOGICAL           :: END
INTEGER,pointer   :: NSHELLINDEX(:),USHELLINDEX(:)
INTEGER,pointer   :: UATOM(:)

natoms = SETTING%MOLECULE(1)%p%nAtoms
call mem_alloc(ATOMtype,natoms)
BAS%nAtoms = nAtoms
NULLIFY(BAS%X)
ALLOCATE(BAS%X(nAtoms))
NULLIFY(BAS%Y)
ALLOCATE(BAS%Y(nAtoms))
NULLIFY(BAS%Z)
ALLOCATE(BAS%Z(nAtoms))
NULLIFY(BAS%CHARGE)
ALLOCATE(BAS%CHARGE(nAtoms))

DO I=1,nAtoms
     BAS%X(I) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1)
     BAS%Y(I) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2)
     BAS%Z(I) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3)
     BAS%CHARGE(I) = INT(SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE)
ENDDO

KMAX=0
MXPRIM=0
R = SETTING%BASIS(1)%p%REGULAR%Labelindex
IF(R.EQ.0)THEN
   DO I=1,nAtoms
      ICHARGE = INT(SETTING%MOLECULE(1)%p%ATOM(I)%Charge)
      type = SETTING%BASIS(1)%p%REGULAR%Chargeindex(ICHARGE)
      ATOMtype(I)=type
      !determine KMAX = number of shells 
      DO K=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%nAngmom
         KMAX=KMAX+SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%norb
         MXPRIM=MXPRIM+SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%nprim
      ENDDO
   ENDDO
ELSE
   DO I=1,nAtoms
      ICHARGE = INT(SETTING%MOLECULE(1)%p%ATOM(I)%Charge)
      type = SETTING%MOLECULE(1)%p%ATOM(I)%IDtype(R)
      ATOMtype(I)=type
      !determine KMAX = number of shells 
      DO K=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%nAngmom
         KMAX=KMAX+SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%norb
         MXPRIM=MXPRIM+SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%nprim
      ENDDO
   ENDDO
ENDIF
BAS%KMAX=KMAX

NULLIFY(BAS%NCENT)
ALLOCATE(BAS%NCENT(KMAX))
NULLIFY(BAS%NSTART)
ALLOCATE(BAS%NSTART(KMAX))
NULLIFY(BAS%NHKT)
ALLOCATE(BAS%NHKT(KMAX))
NULLIFY(BAS%NUCO)
ALLOCATE(BAS%NUCO(KMAX))
NULLIFY(BAS%JSTRT)
ALLOCATE(BAS%JSTRT(KMAX))
NULLIFY(BAS%CENT)
ALLOCATE(BAS%CENT(3,KMAX))

orbitalindex = 0
SHELL = 1
NHTYP=0
TOTPRIM=0
DO I=1,nAtoms
   type=ATOMtype(I)
   DO K=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%nAngmom
      norb=SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%norb
      nprim=SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%nprim
      DO L=1,norb
      BAS%NCENT(SHELL) = I
      BAS%CENT(1,SHELL) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1) 
      BAS%CENT(2,SHELL) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2) 
      BAS%CENT(3,SHELL) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3) 
      BAS%NSTART(SHELL) = orbitalindex
      BAS%JSTRT(SHELL) = TOTPRIM !the accumulated number of primitives
      BAS%NHKT(SHELL) = K
      BAS%NUCO(SHELL) = nprim
      orbitalindex=orbitalindex+(2*(K-1)+1)
      SHELL = SHELL+1
      ENDDO
      TOTPRIM=TOTPRIM+nprim
      NHTYP = MAX(K,NHTYP)
   ENDDO
ENDDO

NULLIFY(BAS%PRIEXP)
ALLOCATE(BAS%PRIEXP(MXPRIM))
MXPRIM = 1
DO I=1,nAtoms
   type=ATOMtype(I)
   DO K=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%nAngmom
      DO L=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%nsegments
         DO M= 1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%segment(L)%nrow
            BAS%PRIEXP(MXPRIM) = SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%segment(L)%Exponents(M)
            MXPRIM = MXPRIM+1
         ENDDO
      ENDDO
   ENDDO
ENDDO
BAS%GRDONE = 0
BAS%NHTYP=NHTYP
BAS%MXPRIM=MXPRIM-1
MXPRIM = MXPRIM-1

!-------------------------------------------------------
!THIS SHOULD NOT BE CONSTUCTED EVERY TIME - JUST SAVE CC,CCINDEX,CCSTART in DALTON
!------------------------------------------------------------
TYP=SETTING%BASIS(1)%p%REGULAR%nAtomtypes
call mem_alloc(MODELTYPES,TYP)
NR=0
DO I=1,SETTING%BASIS(1)%p%REGULAR%natomtypes    
   NR=NR+1
   MODELTYPES(I)=NR
ENDDO
NRSIZE=NR
call mem_alloc(UATOM,NR)
call mem_alloc(NSHELLINDEX,NR)
call mem_alloc(USHELLINDEX,NR)
DO I=1,NR
   UATOM(I) = 0 
ENDDO

USHELL=1
DO I=1,nAtoms
   type=ATOMtype(I)  
   ! TEST IF THIS TYPE OF ATOM HAS ALREADY BEEN PROCESSED -> OLDATOM=TRUE
   NEWATOM=.FALSE.
   NR=MODELTYPES(type)
   IF(UATOM(NR) .EQ. 0) NEWATOM=.TRUE.
   IF(NEWATOM) THEN
     UATOM(NR)=1
     DO K=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%nAngmom
      DO L=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%nsegments
         ncol=SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%ncol
         DO N=1,ncol
            USHELL = USHELL+1
         ENDDO
      ENDDO
     ENDDO
   ENDIF
ENDDO

DO I=1,NRSIZE
   UATOM(I) = 0 
ENDDO

BAS%Ushells=USHELL-1
NULLIFY(BAS%CC)
ALLOCATE(BAS%CC(USHELL))
NULLIFY(BAS%CCSTART)
ALLOCATE(BAS%CCSTART(USHELL))
NULLIFY(BAS%CCINDEX)
ALLOCATE(BAS%CCINDEX(KMAX))

SHELL=1
USHELL=1
DO I=1,nAtoms
   type=ATOMtype(I)  
   ! TEST IF THIS TYPE OF ATOM HAS ALREADY BEEN PROCESSED -> OLDATOM=TRUE
   NEWATOM=.FALSE.
   NR=MODELTYPES(type)
   IF(UATOM(NR) .EQ. 0) NEWATOM=.TRUE.
   IF(NEWATOM) THEN
      UATOM(NR)=I
      USHELLINDEX(NR)=USHELL
      NSHELLINDEX(NR)=0
     DO K=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%nAngmom
         OLDnrow=1
      DO L=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%nsegments
         nrow=SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%nrow
         ncol=SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%ncol
         DO N=1,ncol
            nrow2=0
            DO irow = 1,nrow
               IF(ABS(SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%&
                    &SHELL(K)%SEGMENT(L)%elms(irow+(N-1)*nrow)) > DMIN)THEN
                  nrow2=nrow2+1
               ENDIF
            ENDDO
            if(nrow2 .EQ. nrow)then !all element nonzero!      
               call lsmat_dense_init(BAS%CC(USHELL),nrow,1)
               call dcopy(nrow,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%&
                    &SHELL(K)%SEGMENT(L)%elms(1+(N-1)*nrow:N*nrow),&
                    &1,BAS%CC(USHELL)%elms,1)
               BAS%CCSTART(USHELL)=OLDnrow
               BAS%CCINDEX(SHELL)=USHELL
               SHELL = SHELL+1
               USHELL = USHELL+1
               NSHELLINDEX(NR)=NSHELLINDEX(NR)+1
            ELSE !some zeros
               nrow3=nrow
               END=.FALSE.
               DO irow = 1,nrow
                  IF(END)EXIT
                  IF(ABS(SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%&
                       &SHELL(K)%SEGMENT(L)%elms(irow+(N-1)*nrow)) < DMIN)THEN
                     nrow3=nrow3-1
                  ELSE
                     call lsmat_dense_init(BAS%CC(USHELL),nrow3,1)
                     call dcopy(nrow3,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%&
                          &SHELL(K)%SEGMENT(L)%elms(irow+(N-1)*nrow:N*nrow),&
                          &1,BAS%CC(USHELL)%elms,1)
                     BAS%CCSTART(USHELL)=OLDnrow+irow-1
                     END=.TRUE.
                  ENDIF
               ENDDO
               BAS%CCINDEX(SHELL)=USHELL
               SHELL = SHELL+1
               USHELL = USHELL+1
               NSHELLINDEX(NR)=NSHELLINDEX(NR)+1
            ENDIF
         ENDDO
         OLDnrow=nrow+OLDnrow
      ENDDO
     ENDDO
   ELSE ! OLD ATOM TYPE -> COPY FROM Uniqueatom
      NR=MODELTYPES(type)
      nshells = NSHELLINDEX(NR)
      USHELL2 = USHELLINDEX(NR) 
      DO L=1,nSHELLS
         BAS%CCINDEX(SHELL)=USHELL2
         SHELL = SHELL+1
         USHELL2 = USHELL2+1
      ENDDO
   ENDIF
ENDDO

call mem_dealloc(NSHELLINDEX)
call mem_dealloc(USHELLINDEX)
call mem_dealloc(UATOM)
!DEALLOCATE(NSHELLINDEX)
!DEALLOCATE(USHELLINDEX)
!DEALLOCATE(UATOM)

!---------------------------------------------------------------------------

IF(SHELL-1 .NE. KMAX) CALL LSQUIT('number of shell is wrong in Build_basinf',lupri) 


IF(GRDONE .EQ. 0)THEN
   ! get radii of all shells as defined by specified threshold DFTHRI.
   FACL(1)= 1D0
   FACL(2)= 1.3333D0
   FACL(3)= 1.6D0
   FACL(4)= 1.83D0
   FACL(5)= 2.03D0
   FACL(6)= 2.22D0
   FACL(7)= 2.39D0
   FACL(8)= 2.55D0
   FACL(9)= 2.70D0
   FACL(10)= 2.84D0

   NULLIFY(BAS%RSHEL)
   ALLOCATE(BAS%RSHEL(KMAX))
   CALL LS_DZERO(BAS%RSHEL,KMAX)
   THLOG = LOG(SETTING%SCHEME%DFTHRI)
   SHELL = 1
   DO I=1,nAtoms
      type=ATOMtype(I)
      DO K=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%nAngmom
         DO L=1,SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%nsegments
            nrow=SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%nrow
            ncol=SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%ncol
            DO N=1,ncol
               DO M=1,nrow
                  EXP=SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%SHELL(K)%&
                       &SEGMENT(L)%Exponents(M)
                  R2 = (LOG(ABS(SETTING%BASIS(1)%p%REGULAR%ATOMTYPE(type)%&
                       &SHELL(K)%SEGMENT(L)%elms(M+(N-1)*nrow) ) ) -THLOG)/EXP
                  IF(BAS%RSHEL(SHELL).LT.R2) BAS%RSHEL(SHELL) = R2
               ENDDO
               SHELL = SHELL+1
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   DO SHELL=1,KMAX
      BAS%RSHEL(SHELL) = SQRT(BAS%RSHEL(SHELL))*FACL(BAS%NHKT(SHELL))
   END DO
ELSE
   !not used but allocated to satisfy some compilers
   NULLIFY(BAS%RSHEL)
   ALLOCATE(BAS%RSHEL(1))
ENDIF

call mem_dealloc(ATOMtype)
call mem_dealloc(MODELTYPES)

END SUBROUTINE BUILD_BASINF

!> \brief free the basinf structure used in exchange-correlation calculations
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE FREE_BASINF(BAS)
TYPE(BASINF)      :: BAS
INTEGER           :: ushells,I
ushells = BAS%ushells

DEALLOCATE(BAS%X)
DEALLOCATE(BAS%Y)
DEALLOCATE(BAS%Z)
DEALLOCATE(BAS%CHARGE)
DEALLOCATE(BAS%NCENT)
DEALLOCATE(BAS%NSTART)
DEALLOCATE(BAS%NHKT)
DEALLOCATE(BAS%NUCO)
DEALLOCATE(BAS%JSTRT)
DEALLOCATE(BAS%CENT)
DEALLOCATE(BAS%PRIEXP)
DO I=1,ushells
   call lsmat_dense_free(BAS%CC(I))
ENDDO   
DEALLOCATE(BAS%CC)
DEALLOCATE(BAS%CCSTART)
DEALLOCATE(BAS%CCINDEX)

NULLIFY(BAS%X)
NULLIFY(BAS%Y)
NULLIFY(BAS%Z)
NULLIFY(BAS%CHARGE)
NULLIFY(BAS%NCENT)
NULLIFY(BAS%NSTART)
NULLIFY(BAS%NHKT)
NULLIFY(BAS%NUCO)
NULLIFY(BAS%JSTRT)
NULLIFY(BAS%CENT)
NULLIFY(BAS%PRIEXP)
NULLIFY(BAS%CC)
NULLIFY(BAS%CCSTART)

!RSHEL SHALL NOT ALLOCATED BECAUSE
!THAT HAPPENS INSIDE THE C CODE
DEALLOCATE(BAS%RSHEL)
NULLIFY(BAS%RSHEL)

END SUBROUTINE FREE_BASINF

SUBROUTINE ADD_TO_CC(CC,STARTPRIM,NROW,TMP,N,MXPRIM,KMAX)
IMPLICIT NONE
INTEGER     :: M,STARTPRIM,NROW,N,KMAX,MXPRIM
REAL(REALK) :: CC(MXPRIM*KMAX),TMP(:)

DO M=1,nrow
   CC(STARTPRIM+M) = TMP(M+(N-1)*nrow)  
ENDDO

END SUBROUTINE ADD_TO_CC

SUBROUTINE SHELL_NUMBER(LUPRI,IPRINT,DALTON,BASISINFO,KMAX)
IMPLICIT NONE
INTEGER           :: LUPRI,IPRINT
TYPE(DALTONINPUT) :: DALTON
TYPE(BASISSETINFO)        :: BASISINFO
INTEGER           :: KMAX !number of shells see type def of BASINF
!
INTEGER           :: I,R,type,K,icharge

KMAX=0
R = BASISINFO%Labelindex

IF(R.EQ.0)THEN
   DO I=1,DALTON%MOLECULE%nAtoms
      ICHARGE = INT(DALTON%MOLECULE%ATOM(I)%Charge)
      type = BASISINFO%Chargeindex(ICHARGE)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         KMAX=KMAX+BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
      ENDDO
   ENDDO
ELSE
   DO I=1,DALTON%MOLECULE%nAtoms
      type = DALTON%MOLECULE%ATOM(I)%IDtype(R)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         KMAX=KMAX+BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE SHELL_NUMBER

SUBROUTINE ORBITAL_DEGENERACY(LUPRI,IPRINT,DALTON,BASISINFO,KMAX,DEGENERACY)
IMPLICIT NONE
INTEGER           :: LUPRI,IPRINT,KMAX!number of shells
TYPE(DALTONINPUT) :: DALTON
TYPE(BASISSETINFO)        :: BASISINFO
INTEGER           :: DEGENERACY(KMAX) 
!
INTEGER           :: I,R,type,K,shell,norb,L,ICHARGE
R = BASISINFO%Labelindex

IF(R.EQ.0)THEN
   SHELL = 1
   DO I=1,DALTON%MOLECULE%nAtoms
      ICHARGE = INT(DALTON%MOLECULE%ATOM(I)%Charge)
      type = BASISINFO%Chargeindex(ICHARGE)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         norb=BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
         DO L=1,norb
            DEGENERACY(SHELL)=(2*(K-1)+1)
            SHELL = SHELL+1
         ENDDO
      ENDDO
   ENDDO
ELSE
   SHELL = 1
   DO I=1,DALTON%MOLECULE%nAtoms
      type = DALTON%MOLECULE%ATOM(I)%IDtype(R)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         norb=BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
         DO L=1,norb
            DEGENERACY(SHELL)=(2*(K-1)+1)
            SHELL = SHELL+1
         ENDDO
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE ORBITAL_DEGENERACY

SUBROUTINE ATOMIC_ORBITALS(LUPRI,IPRINT,DALTON,BASISINFO,natoms,orbitals)
IMPLICIT NONE
INTEGER           :: LUPRI,IPRINT,natoms!number of shells
TYPE(DALTONINPUT) :: DALTON
TYPE(BASISSETINFO)        :: BASISINFO
INTEGER           :: orbitals(natoms) 
!
INTEGER           :: I,R,type,K,norb,ICHARGE

CALL LS_IZERO(ORBITALS,natoms)

R = BASISINFO%Labelindex
IF(R.EQ.0)THEN
   DO I=1,DALTON%MOLECULE%nAtoms
      ICHARGE = INT(DALTON%MOLECULE%ATOM(I)%Charge)
      type = BASISINFO%Chargeindex(ICHARGE)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         norb=BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
         ORBITALS(I)=ORBITALS(I)+(2*(K-1)+1)*norb
      ENDDO
   ENDDO
ELSE
   DO I=1,DALTON%MOLECULE%nAtoms
      type = DALTON%MOLECULE%ATOM(I)%IDtype(R)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         norb=BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
         ORBITALS(I)=ORBITALS(I)+(2*(K-1)+1)*norb
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE ATOMIC_ORBITALS

END MODULE BUILDAOBATCH


