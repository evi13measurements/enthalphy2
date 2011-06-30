!> @file
!> Module containing information about the OD-batches (combination of zero, one or two AO-batches)
MODULE ODbatches
 use precision
 use typedef
 use memory_handling
 use ls_util
 use OD_type
 !> An OD-batch is a set of:
 !>  i)   product overlaps between two batches of basis-functions -
 !>       where each AO-batch belong to a given center, and
 !>       share a common set of primitive basis-functions
 !>  ii)  AO basis-functions (beloning to one
 !>       center and with a set of primitive functions).
 !>  iii) an empty batch (used for debugging purposes only -
 !>       useful when testing contruction of E-coefficients).
 !>
CONTAINS
!> \brief Creates OD-batches from two sets of AO-batches
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param OD ODbatch to be built
!> \param INPUT integral input which contain all input to the integral eval.
!> \param SIDE label LHS or RHS 
!> \param IUNIT the logical unit number for the output file
SUBROUTINE Create_ODbatches(OD,INPUT,SIDE,IUNIT)
IMPLICIT NONE
 TYPE(ODITEM)      :: OD
 TYPE(INTEGRALINPUT),target :: INPUT
 LOGICAL           :: sameAOs
! Character(len=80) :: IDENTIFIER
 Character*(*)     :: SIDE
 Integer           :: IOD,IA,IB,IBSTART
 Integer           :: SA
 Integer           :: AOA,AOB,IUNIT
 Logical           :: screen
 Real(realk)       :: maxGab,maxGabElm
 Integer           :: i2,atomA,atomB,Gindex
 Integer           :: batchA,batchB
! real(realk),parameter    :: R2PI52 = 5.91496717279561287782D0
 type(lstensor),pointer :: GAB
!                     R2PI52 = sqrt(2 * sqrt(PI^5) )
 SELECT CASE(SIDE)
 CASE('LHS')
    sameAOs=INPUT%sameLHSaos
    AOA = 1
    AOB = 2
    NULLIFY(GAB)
    IF(INPUT%CS_SCREEN)THEN
       GAB => Input%LST_GAB_LHS
    ENDIF
    CALL INIT_OD(OD,INPUT%AO(AOA)%p,INPUT%AO(AOB)%p,sameAOs,GAB,&
     &           Input%CS_Screen,Input%OD_Screen,Input%CS_threshold,Input%CS_MAXELM_RHS)
 CASE('RHS')
    sameAOs=INPUT%sameRHSaos
    AOA = 3
    AOB = 4
    IF(INPUT%CS_SCREEN)GAB => Input%LST_GAB_RHS
    CALL INIT_OD(OD,INPUT%AO(AOA)%p,INPUT%AO(AOB)%p,sameAOs,GAB,&
     &           Input%CS_Screen,Input%OD_Screen,Input%CS_threshold,Input%CS_MAXELM_LHS)
 CASE DEFAULT
    WRITE(IUNIT,'(1X,2A)') 'Wrong case in Create_ODbatches =',SIDE
    CALL LSQUIT('Wrong case in Create_ODbatches',-1)
 END SELECT

! OD%identifier = IDENTIFIER  -  not used
 OD%sameAOs    = sameAOs
 IOD = 0
 DO IA=1,INPUT%AO(AOA)%p%nbatches
   IBSTART = 1
   IF (sameAOs) IBSTART = IA
   DO IB=IBSTART,INPUT%AO(AOB)%p%nbatches
!    Determine screening criteria
     screen = .FALSE.
!    Cauchy-Scwartz screening |(a|b)| <= sqrt{(a|a)} sqrt{(b|b)}
     IF(INPUT%CS_SCREEN)THEN
        atomA = INPUT%AO(AOA)%p%BATCH(IA)%atom
        atomB = INPUT%AO(AOB)%p%BATCH(IB)%atom
        batchA = INPUT%AO(AOA)%p%BATCH(IA)%batch
        batchB = INPUT%AO(AOB)%p%BATCH(IB)%batch
        GINDEX = GAB%INDEX(atomA,atomB,1,1)
        !THE MAXIMUM ELEMENT IN THE SCHWARZ(GAB) MATRIX FROM THE 2 AOBATCHES
        IF(SIDE(1:3) .EQ. 'LHS')THEN 
           CALL CALC_MAXGAB(maxGab,&
                &GAB%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%elms,&
                &GAB%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%nelms)
           maxGabElm = INPUT%CS_MAXELM_RHS
        ELSE
           CALL CALC_MAXGAB(maxGab,&
                &GAB%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%elms,&
                &GAB%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%nelms)
           maxGabElm = INPUT%CS_MAXELM_LHS
        ENDIF
        IF (maxGab.LE.Input%CS_Threshold/maxGabElm) screen = .TRUE.
     ENDIF
!    Extent screening
     IF (INPUT%OD_SCREEN) THEN
       CALL getODscreening(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB),screen) 
     ENDIF

     IF (.NOT.screen) THEN
!      Create new OD-batch
       IOD         = IOD + 1
       OD%BATCH(IOD)%IA            =  IA
       OD%BATCH(IOD)%IB            =  IB
       OD%BATCH(IOD)%AO(1)%p       => INPUT%AO(AOA)%p%BATCH(IA)
       OD%BATCH(IOD)%AO(2)%p       => INPUT%AO(AOB)%p%BATCH(IB)
       OD%BATCH(IOD)%nAngmom       =  INPUT%AO(AOA)%p%BATCH(IA)%nAngmom*INPUT%AO(AOB)%p%BATCH(IB)%nAngmom
       OD%BATCH(IOD)%nPrimitives   =  INPUT%AO(AOA)%p%BATCH(IA)%nPrimitives*INPUT%AO(AOB)%p%BATCH(IB)%nPrimitives
       OD%BATCH(IOD)%maxContracted =  INPUT%AO(AOA)%p%BATCH(IA)%maxContracted*INPUT%AO(AOB)%p%BATCH(IB)%maxContracted
       OD%BATCH(IOD)%maxGAB = 0
       OD%BATCH(IOD)%sameAO = (sameAOs .AND. IA .EQ. IB)
!
       IF (OD%BATCH(IOD)%sameAO) OD%BATCH(IOD)%nAngmom &
     &       = INPUT%AO(AOA)%p%BATCH(IA)%nAngmom*(INPUT%AO(AOA)%p%BATCH(IA)%nAngmom+1)/2
!     
       IF(INPUT%CS_SCREEN) OD%BATCH(IOD)%maxGAB  = maxGab

!      Non-classical screening screen away all contributions that can be calculated 
!      classically to within some threshold. For this screning we therefore set up
!      an (effective) center as well as the non-classical extent of the OD-batches
       IF(INPUT%DO_MULMOM.OR.INPUT%NonClassical_SCREEN)THEN
         CALL getODcenter(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB), &
     &                    OD%BATCH(IOD),INPUT%MM_SCREENTHR,IUNIT)
         CALL getODextentNonClassical(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB),&
     &                                OD%BATCH(IOD),INPUT%MM_SCREENTHR,IUNIT)
!      For overlap-integrals we screen if the distance between the LHS and RHS OB-batches exceeds
!      the sum of their extents. The OD-batch centers are the same weigthed center as for the 
!      non-classical screening, the primitive extents are determined by when the (gaussian product)
!      OD-function becomes smaller than a given threshold. As for non-classical screening, the
!      maximum primitive extent, augmented by the distance to the weigthed OD-center, is 
!      taken as the extent.
       ELSEIF (INPUT%OE_Screen) THEN
         CALL getODcenter(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB), &
     &                    OD%BATCH(IOD),INPUT%OE_THRESHOLD,IUNIT)
         CALL getODextentOverlap(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB),&
     &                           OD%BATCH(IOD),INPUT%OE_THRESHOLD,IUNIT)
       ENDIF

    ENDIF
 ENDDO
ENDDO

 IF (IOD.NE.OD%nbatches) THEN
   WRITE(IUNIT,'(1X,1A,I5)') 'Number of batches from Create_ODbatches', IOD
   WRITE(IUNIT,'(1X,1A,I5)') 'Number of batches from Get_Nodbatches  ', OD%nbatches
   WRITE(IUNIT,'(1X,1A)')    'Programming error: Mismatch between number of OD-batches'
   CALL LSQUIT('Mismatch between number of OD-batches in Create_ODbatches',-1)
 ENDIF
! CALL DETERMINE_ODITEM_MEM(OD)

END SUBROUTINE Create_ODbatches

!!$SUBROUTINE DETERMINE_ODITEM_MEM(OD)
!!$IMPLICIT NONE
!!$TYPE(ODitem)        :: OD
!!$INTEGER(KIND=long)  :: mem_alloc
!!$INTEGER             :: I
!!$
!!$mem_alloc = 0
!!$mem_alloc = mem_alloc+sizeof(OD%sameAOs)
!!$DO I = 1,OD%nbatches
!!$!mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%identifier)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%AO)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%nAngmom)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%nPrimitives)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%maxContracted)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%maxGAB)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%sameAO)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%IA)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%IB)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%ODcenter)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%ODextent)
!!$ENDDO
!!$
!!$END SUBROUTINE DETERMINE_ODITEM_MEM

!> \brief Determine the center of the overlap distributions
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param ODB the overlap distribution 
!> \param THRESH threshold for zero contribution   
!> \param IUNIT the logical unit number for the output file
SUBROUTINE getODcenter(AOA,AOB,ODB,THRESH,IUNIT)
IMPLICIT NONE
 TYPE(AOBATCH)       :: AOA,AOB
 TYPE(ODBATCH)       :: ODB
 Integer           :: IOD
 Integer           :: IUNIT
 real(realk)       :: X,Y,Z,DIST12,THRESH,SM,AM,BM
 !real(realk),allocatable  :: EXTPMAX(:),primcenter(:,:) 
 real(realk)       :: e1,e2,p,mu,MAXCC1SUM,CCSUM,MAXCC2SUM,pm1,WT
 real(realk)       :: FAC,DIST
 real(realk),pointer :: CC1SUM(:),CC2SUM(:)
 Integer           :: nP1,nP2,i2,i1,i12,iA1,iA2,iC1,iC2,n12
 real(realk)       :: R2PI52 = 5.91496717279561287782D0
 logical           :: nonzero,nucleus

nucleus = AOA%type_Nucleus .OR. AOB%type_Nucleus
IF (nucleus) THEN
  ODB%ODextent = 0.d0
!Simen     Set up an artificial extent of 1 because the error in the electron-nuclear attraction
!Simen     is larger than the electronic repulsion contribution. Should be analyzed!
!Simen     Matched by an exual extent set up in mm_read_in_raw_data (mm_interface.f90)
  ODB%ODextent = 1.d0
  IF (AOA%type_Nucleus) THEN
    ODB%ODcenter(1) = AOA%center(1)
    ODB%ODcenter(2) = AOA%center(2)
    ODB%ODcenter(3) = AOA%center(3)
  ELSE IF (AOB%type_Nucleus) THEN
    ODB%ODcenter(1) = AOB%center(1)
    ODB%ODcenter(2) = AOB%center(2)
    ODB%ODcenter(3) = AOB%center(3)
  ELSE
    CALL LSQUIT('Error in getODcenter, both AO-batches are nuclei!',-1)
  ENDIF
ELSE
  X = AOA%center(1)-AOB%center(1)
  Y = AOA%center(2)-AOB%center(2)
  Z = AOA%center(3)-AOB%center(3)
  DIST12 = X*X+Y*Y+Z*Z
  !NOT YET IMPLEMENTED HARDCODED FOR NOW
  SM=0.0D0
  AM=0.0D0
  BM=0.0D0
  NONZERO = .FALSE.
  call mem_alloc(CC1SUM,AOA%nAngmom)
  call mem_alloc(CC2SUM,AOB%nAngmom)
  !ALLOCATE(CC1SUM(AOA%nAngmom))
  !ALLOCATE(CC2SUM(AOB%nAngmom))
  nP2=AOB%nPrimitives
  nP1=AOA%nPrimitives
  i12=0
  DO i2=1,nP2
     e2  = ODB%AO(2)%p%pexponents%elms(i2)
     DO i1=1,nP1
        e1  = ODB%AO(1)%p%pexponents%elms(i1)
        p = e1+e2
        pm1 = 1/p
        mu = e1*e2*pm1
        CC1SUM = 0.D0
        CC2SUM = 0.D0
        DO iA1 = 1,AOA%nAngmom
           DO iC1=1,AOA%nContracted(iA1)
              CC1SUM(iA1) = CC1SUM(iA1) + ABS(AOA%pCC(iA1)%p%elms(i1+(iC1-1)*nP1))
           ENDDO
        ENDDO
        MAXCC1SUM = 0.d0
        DO iA1 = 1,AOA%nAngmom
           MAXCC1SUM = MAX(MAXCC1SUM,CC1SUM(iA1)) 
        ENDDO
        
        DO iA2 = 1,AOB%nAngmom
           DO iC2=1,AOB%nContracted(iA2)
              CC2SUM(iA2) = CC2SUM(iA2) + ABS(AOB%pCC(iA2)%p%elms(i2+(iC2-1)*nP2))
           ENDDO
        ENDDO
        MAXCC2SUM = 0.d0
        DO iA2 = 1,AOB%nAngmom
           MAXCC2SUM = MAX(MAXCC2SUM,CC2SUM(iA2)) 
        ENDDO
        CCSUM = MAXCC1SUM*MAXCC2SUM
        
        FAC = R2PI52*EXP(-MU*DIST12)*pm1*CCSUM
        IF(FAC .GT. THRESH)THEN
           NONZERO = .TRUE.
           MAXCC1SUM = 0.d0
           DO iA1 = 1,AOA%nAngmom
              MAXCC1SUM = MAXCC1SUM+CC1SUM(iA1) 
           ENDDO
           MAXCC2SUM = 0.d0
           DO iA2 = 1,AOB%nAngmom
              MAXCC2SUM = MAXCC2SUM+CC2SUM(iA2) 
           ENDDO
           
           WT = MAXCC1SUM*MAXCC2SUM
           
           i12 = i12 + 1
           
           WT = WT * EXP(-MU*DIST12)
           SM = SM + WT
           AM = AM + WT * e1 * pm1
           BM = BM + WT * e2 * pm1
        ENDIF
     ENDDO
  ENDDO
  n12=i12
  call mem_dealloc(CC1SUM)
  call mem_dealloc(CC2SUM)
  !DEALLOCATE(CC1SUM)
  !DEALLOCATE(CC2SUM)
  
  IF(NONZERO)THEN
    AM = AM/SM
    BM = BM/SM
    ODB%ODcenter(1) = AM*AOA%center(1) + BM*AOB%center(1)
    ODB%ODcenter(2) = AM*AOA%center(2) + BM*AOB%center(2)
    ODB%ODcenter(3) = AM*AOA%center(3) + BM*AOB%center(3)
  ELSE
    ODB%ODcenter(1) = (AOA%center(1) + AOB%center(1))/2.D0
    ODB%ODcenter(2) = (AOA%center(2) + AOB%center(2))/2.D0
    ODB%ODcenter(3) = (AOA%center(3) + AOB%center(3))/2.D0
  ENDIF

ENDIF
END SUBROUTINE getODcenter

!> \brief Determine the non classical overlap distribution extent
!> \author S. Reine
!> \date 2010
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param ODB the overlap distribution 
!> \param THRESH threshold for extent (when the OD falls below THRESH)
!> \param IUNIT the logical unit number for the output file
SUBROUTINE getODextentNonClassical(AOA,AOB,ODB,THRESH,IUNIT)
IMPLICIT NONE
 TYPE(AOBATCH)       :: AOA,AOB
 TYPE(ODBATCH)       :: ODB
 Integer           :: IOD
 Integer           :: IUNIT
 real(realk)       :: LS_ERFCIV, X,Y,Z,DIST12,THRESH,EXTPMAX2,FACCLS
 real(realk),pointer  :: EXTPMAX(:),primcenter(:,:) 
 real(realk)       :: e1,e2,p,mu,MAXCC1SUM,CCSUM,MAXCC2SUM,pm1,EXTP
 real(realk)       :: FAC,dist
 real(realk),pointer :: CC1SUM(:),CC2SUM(:)
 Integer           :: nP1,nP2,i2,i1,i12,iA1,iA2,iC1,iC2,n12
 real(realk),parameter       :: R2PI52 = 5.91496717279561287782D0
 logical           :: nucleus

nucleus = AOA%type_Nucleus .OR. AOB%type_Nucleus
IF (nucleus) THEN
  ODB%ODextent = 0.d0
!Simen     Set up an artificial extent of 1 because the error in the electron-nuclear attraction
!Simen     is larger than the electronic repulsion contribution. Should be analyzed!
!Simen     Matched by an exual extent set up in mm_read_in_raw_data (mm_interface.f90)
  ODB%ODextent = 1.d0
ELSE
  X = AOA%center(1)-AOB%center(1)
  Y = AOA%center(2)-AOB%center(2)
  Z = AOA%center(3)-AOB%center(3)
  DIST12 = X*X+Y*Y+Z*Z
  !NOT YET IMPLEMENTED HARDCODED FOR NOW
  FACCLS = LS_ERFCIV(THRESH)
  EXTPMAX2=0.0D0
  call mem_alloc(CC1SUM,AOA%nAngmom)
  call mem_alloc(CC2SUM,AOB%nAngmom)
  !ALLOCATE(CC1SUM(AOA%nAngmom))
  !ALLOCATE(CC2SUM(AOB%nAngmom))
  nP2=AOB%nPrimitives
  nP1=AOA%nPrimitives
  call mem_alloc(EXTPMAX,nP1*nP2)
  call mem_alloc(Primcenter,3,nP1*nP2)
  !ALLOCATE(EXTPMAX(nP1*nP2))
  !ALLOCATE(Primcenter(3,nP1*nP2))
  i12=0
  DO i2=1,nP2
     e2  = ODB%AO(2)%p%pexponents%elms(i2)
     DO i1=1,nP1
        e1  = ODB%AO(1)%p%pexponents%elms(i1)
        p = e1+e2
        pm1 = 1/p
        mu = e1*e2*pm1
        CC1SUM = 0.D0
        CC2SUM = 0.D0
        DO iA1 = 1,AOA%nAngmom
           DO iC1=1,AOA%nContracted(iA1)
              CC1SUM(iA1) = CC1SUM(iA1) + ABS(AOA%pCC(iA1)%p%elms(i1+(iC1-1)*nP1))
           ENDDO
        ENDDO
        MAXCC1SUM = 0.d0
        DO iA1 = 1,AOA%nAngmom
           MAXCC1SUM = MAX(MAXCC1SUM,CC1SUM(iA1)) 
        ENDDO
        
        DO iA2 = 1,AOB%nAngmom
           DO iC2=1,AOB%nContracted(iA2)
              CC2SUM(iA2) = CC2SUM(iA2) + ABS(AOB%pCC(iA2)%p%elms(i2+(iC2-1)*nP2))
           ENDDO
        ENDDO
        MAXCC2SUM = 0.d0
        DO iA2 = 1,AOB%nAngmom
           MAXCC2SUM = MAX(MAXCC2SUM,CC2SUM(iA2)) 
        ENDDO
        CCSUM = MAXCC1SUM*MAXCC2SUM
        
        FAC = R2PI52*EXP(-MU*DIST12)*pm1*CCSUM
        IF(FAC .GT. THRESH)THEN
           i12 = i12 + 1
           
           EXTP = SQRT(pm1) * FACCLS
           EXTPMAX(i12) = EXTP

           Primcenter(1,i12) = (e1*AOA%center(1) + e2*AOB%center(1))/(e1+e2)
           Primcenter(2,i12) = (e1*AOA%center(2) + e2*AOB%center(2))/(e1+e2)
           Primcenter(3,i12) = (e1*AOA%center(3) + e2*AOB%center(3))/(e1+e2)
        ENDIF
     ENDDO
  ENDDO
  n12=i12
  call mem_dealloc(CC1SUM)
  call mem_dealloc(CC2SUM)
  !DEALLOCATE(CC1SUM)
  !DEALLOCATE(CC2SUM)
  
  EXTPMAX2 = 0.D0
  DO I12=1,n12
     X = ODB%ODcenter(1) - Primcenter(1,i12)
     Y = ODB%ODcenter(2) - Primcenter(2,i12)
     Z = ODB%ODcenter(3) - Primcenter(3,i12)
     DIST = SQRT(X*X+Y*Y+Z*Z)
     EXTPMAX2 = MAX(EXTPMAX2,EXTPMAX(i12)+DIST)
  ENDDO
  ODB%ODextent = EXTPMAX2
  call mem_dealloc(EXTPMAX)
  call mem_dealloc(Primcenter)
  !DEALLOCATE(EXTPMAX)
  !DEALLOCATE(Primcenter)
ENDIF
END SUBROUTINE getODextentNonClassical

!> \brief Determine the overlap distribution extent
!> \author S. Reine
!> \date 2010
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param ODB the overlap distribution 
!> \param THRESH threshold for extent (when the OD falls below THRESH)
!> \param IUNIT the logical unit number for the output file
SUBROUTINE getODextentOverlap(AOA,AOB,ODB,THRESH,IUNIT)
IMPLICIT NONE
 TYPE(AOBATCH)       :: AOA,AOB
 TYPE(ODBATCH)       :: ODB
 Integer           :: IOD
 Integer           :: IUNIT
 real(realk)       :: X,Y,Z,DIST12,THRESH,EXTPMAX2
 real(realk),pointer  :: EXTPMAX(:),primcenter(:,:) 
 real(realk)       :: e1,e2,p,mu,MAXCC1SUM,CCSUM,MAXCC2SUM,pm1,EXTP
 real(realk)       :: FACSCREEN,FAC,DIST
 real(realk),pointer :: CC1SUM(:),CC2SUM(:)
 Integer           :: nP1,nP2,i2,i1,i12,iA1,iA2,iC1,iC2,n12
 real(realk)       :: R2PI52 = 5.91496717279561287782D0
 logical           :: nucleus

nucleus = AOA%type_Nucleus .OR. AOB%type_Nucleus
IF (nucleus) THEN
  CALL LSQUIT('Not an option getODextentOverlap with nuclei',-1)
ELSE
  X = AOA%center(1)-AOB%center(1)
  Y = AOA%center(2)-AOB%center(2)
  Z = AOA%center(3)-AOB%center(3)
  DIST12 = X*X+Y*Y+Z*Z
  EXTPMAX2=0.0D0
  call mem_alloc(CC1SUM,AOA%nAngmom)
  call mem_alloc(CC2SUM,AOB%nAngmom)
  !ALLOCATE(CC1SUM(AOA%nAngmom))
  !ALLOCATE(CC2SUM(AOB%nAngmom))
  nP2=AOB%nPrimitives
  nP1=AOA%nPrimitives
  call mem_alloc(EXTPMAX,nP1*nP2)
  call mem_alloc(Primcenter,3,nP1*nP2)
  !ALLOCATE(EXTPMAX(nP1*nP2))
  !ALLOCATE(Primcenter(3,nP1*nP2))
  i12=0
  DO i2=1,nP2
     e2  = ODB%AO(2)%p%pexponents%elms(i2)
     DO i1=1,nP1
        e1  = ODB%AO(1)%p%pexponents%elms(i1)
        p = e1+e2
        pm1 = 1/p
        mu = e1*e2*pm1
        CC1SUM = 0.D0
        CC2SUM = 0.D0
        DO iA1 = 1,AOA%nAngmom
           DO iC1=1,AOA%nContracted(iA1)
              CC1SUM(iA1) = CC1SUM(iA1) + ABS(AOA%pCC(iA1)%p%elms(i1+(iC1-1)*nP1))
           ENDDO
        ENDDO
        MAXCC1SUM = 0.d0
        DO iA1 = 1,AOA%nAngmom
           MAXCC1SUM = MAX(MAXCC1SUM,CC1SUM(iA1)) 
        ENDDO
        
        DO iA2 = 1,AOB%nAngmom
           DO iC2=1,AOB%nContracted(iA2)
              CC2SUM(iA2) = CC2SUM(iA2) + ABS(AOB%pCC(iA2)%p%elms(i2+(iC2-1)*nP2))
           ENDDO
        ENDDO
        MAXCC2SUM = 0.d0
        DO iA2 = 1,AOB%nAngmom
           MAXCC2SUM = MAX(MAXCC2SUM,CC2SUM(iA2)) 
        ENDDO
        CCSUM = MAXCC1SUM*MAXCC2SUM
        
!Simen  Not proper for overlap integrals, refine!
        FAC = R2PI52*EXP(-MU*DIST12)*pm1*CCSUM
        IF(FAC .GT. THRESH)THEN
           i12 = i12 + 1
           
           EXTP = getODbatchOverlapExtent(p,AOA%maxAngmom+AOB%maxAngmom,&
     &                                    EXP(-MU*DIST12)*CCSUM,THRESH)
           
           EXTPMAX(i12) = EXTP

           Primcenter(1,i12) = (e1*AOA%center(1) + e2*AOB%center(1))/(e1+e2)
           Primcenter(2,i12) = (e1*AOA%center(2) + e2*AOB%center(2))/(e1+e2)
           Primcenter(3,i12) = (e1*AOA%center(3) + e2*AOB%center(3))/(e1+e2)
        ENDIF
     ENDDO
  ENDDO
  n12=i12
  call mem_dealloc(CC1SUM)
  call mem_dealloc(CC2SUM)
  !DEALLOCATE(CC1SUM)
  !DEALLOCATE(CC2SUM)
  
  EXTPMAX2 = 0.D0
  DO I12=1,n12
     X = ODB%ODcenter(1) - Primcenter(1,i12)
     Y = ODB%ODcenter(2) - Primcenter(2,i12)
     Z = ODB%ODcenter(3) - Primcenter(3,i12)
     DIST = SQRT(X*X+Y*Y+Z*Z)
     EXTPMAX2 = MAX(EXTPMAX2,EXTPMAX(i12)+DIST)
  ENDDO
  ODB%ODextent = EXTPMAX2
  call mem_dealloc(EXTPMAX)
  call mem_dealloc(Primcenter)
  !DEALLOCATE(EXTPMAX)
  !DEALLOCATE(Primcenter)
ENDIF
END SUBROUTINE getODextentOverlap

!> \brief Determine the overlap distribution extent
!> \author S. Reine
!> \date 2010
!> \param exponent the p=p1+p2 exponent of the OD 
!> \param angmax the maximum angular momentum
!> \param prefactor the prefactor exp(-mu *Rab), mu is reduced exponent
!> \param Threshold the threshold
!>
!>Function to return the extent of gaussian solid harmonic function;
!>   prefactor*S^{angmax,m}(r_P)*exp(-p r_R^2) 
!>      .LE. prefactor*r_P^l * exp(-p r_R^2) .LE. threshold
!>
REAL(REALK) FUNCTION getODbatchOverlapExtent(exponent,angmax,prefactor,threshold)
implicit none
REAL(REALK) :: exponent,prefactor,threshold
INTEGER     :: angmax
!
REAL(REALK) :: r2

r2 = -(1/exponent)*log(threshold/prefactor)
!Simen: Not exact for l>0 refine!
IF (r2.GT.1.d0) r2=r2+angmax*log(sqrt(r2))
getODbatchOverlapExtent = sqrt(r2)

END FUNCTION getODbatchOverlapExtent

!> \brief Determine the number of OD-batches
!> \author S. Reine
!> \date 2010
!> \param nbatches the number of batches to be determined
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param sameAOs if the AOs on AOA and AOB is the same
!> \param GAB the screening lstensor
!> \param CS_screen if Cauchy-Schwarz screening is used
!> \param OD_screen if the overlap distribution screening is used
!> \param CS_threshold the Cauchy-Schwarz screening threshold
!> \param maxgabelm the maximum element of the gab matrix
SUBROUTINE GET_NODBATCHES(nbatches,AOA,AOB,sameAOs,GAB,CS_Screen,OD_Screen,CS_threshold,maxGabElm)
IMPLICIT NONE
 Integer      :: nbatches,IA,IB,IBSTART
 TYPE(AOITEM) :: AOA, AOB
 LOGICAL      :: sameAOs
! Real(realk)  :: GAB(:,:)
type(lstensor),pointer :: GAB
 Logical      :: CS_Screen, OD_Screen
 Real(realk)  :: CS_threshold,maxGabElm
!
 Integer      :: atomA,atomB,Gindex,batcha,BatchB
 Real(realk)  :: maxGab
 Logical      :: screen

 nbatches=0
 DO IA=1,AOA%nbatches
   IBSTART = 1
   IF (sameAOs) IBSTART = IA
   DO IB=IBSTART,AOB%nbatches
     screen = .FALSE.
     IF(CS_SCREEN)THEN
        atomA = AOA%BATCH(IA)%atom
        BatchA = AOA%BATCH(IA)%batch
        atomB = AOB%BATCH(IB)%atom
        BatchB = AOB%BATCH(IB)%batch
        GINDEX = GAB%INDEX(atomA,atomB,1,1)
        CALL CALC_MAXGAB(maxGab,&
             &GAB%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%elms,&
             &GAB%LSAO(Gindex)%BATCH(batchA,batchB,1,1)%nelms)
        IF (maxGab.LE.CS_Threshold/maxGabElm) screen = .TRUE.
     ENDIF
     IF (OD_SCREEN) THEN
        CALL getODscreening(AOA%BATCH(IA),AOB%BATCH(IB),screen) 
     ENDIF
     IF (.NOT.screen) THEN
       nbatches = nbatches + 1
     ENDIF
   ENDDO
 ENDDO
END SUBROUTINE GET_NODBATCHES

!> \brief Initialize OD-item
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param OD the Overlap distribution 
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param sameAOs if the AOs on AOA and AOB is the same
!> \param GAB the screening lstensor
!> \param CS_screen if Cauchy-Schwarz screening is used
!> \param OD_screen if the overlap distribution screening is used
!> \param CS_threshold the Cauchy-Schwarz screening threshold
!> \param maxgabelm the maximum element of the gab matrix
SUBROUTINE Init_OD(OD,AOA,AOB,sameAOs,GAB,CS_Screen,OD_Screen,CS_threshold,maxGabElm)
  use memory_handling
IMPLICIT NONE
 TYPE(ODITEM)   :: OD
 TYPE(AOITEM) :: AOA, AOB
 Integer      :: nbatches,SIZE
 LOGICAL      :: sameAOs
 type(lstensor),pointer :: GAB
! Real(realk)  :: GAB(:,:) 
 Logical      :: OD_Screen,CS_Screen
 Real(realk)  :: CS_threshold,maxGabElm
!
 CALL GET_NODBATCHES(nbatches,AOA,AOB,sameAOs,GAB,CS_Screen,OD_Screen,CS_threshold,maxGabElm)
OD%nbatches = nbatches
 IF(nbatches.GT. 0)THEN
    Nullify(OD%BATCH)
    Allocate(OD%BATCH(nbatches))
 ENDIF
 !SIZE = SIZE OF ODBATCH
 SIZE = 5*mem_realsize+5*mem_intsize+1*mem_logicalsize
 mem_allocated_ODitem = mem_allocated_ODitem + SIZE*nbatches 
 max_mem_used_ODitem = MAX(max_mem_used_ODitem,mem_allocated_ODitem)
END SUBROUTINE Init_OD
!
!> \brief calculate the maximum gab element
!> \author T. Kjaergaard
!> \date 2010
!> \param maxgab the maximum gab element
!> \param GAB an array/block/part of the screening tensor 
!> \param ndim the dimension of the small array 
SUBROUTINE CALC_MAXGAB(maxGab,GAB,ndim)
IMPLICIT NONE
INTEGER       :: I,ndim
real(realk)   :: maxGab
real(realk)   :: GAB(ndim)

maxGab = 0.0d0
DO I=1,ndim
   maxGAB = MAX(maxGAB,GAB(I)) 
ENDDO

END SUBROUTINE CALC_MAXGAB

END MODULE ODbatches


!> \brief 
!> \author 
!> \date 
!> \param THR the threshold 
FUNCTION LS_ERFCIV(THR)
 use precision
implicit none
real(realk),PARAMETER :: D0 = 0.0D0,THRMAX = 1.0D-15
real(realk)           :: ERFCI(0:15),LS_ERFCIV,THR
REAL(realk)           :: ARG

  ERFCI = (/0.00000D0, 1.16309D0, 1.82139D0, 2.32675D0, &
     &                       2.75106D0, 3.12341D0, 3.45891D0, 3.76656D0,& 
     &                       4.05224D0, 4.32001D0, 4.57282D0, 4.81292D0,&
     &                       5.04203D0, 5.26151D0, 5.47248D0, 5.67585D0/)
  ARG = MAX(MAX(ABS(THR),THRMAX),D0)
  LS_ERFCIV = ERFCI(-INT(DLOG10(ARG)))
END FUNCTION
