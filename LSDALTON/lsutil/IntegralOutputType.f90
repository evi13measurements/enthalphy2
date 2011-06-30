!> @file 
!> Contains OBJECT CONTAINING INFORMATION ABOUT THE INTEGRAL OUTPUT
MODULE integraloutput_type
  use precision
  use lstensor_operationsmod

TYPE contractionrule
integer :: OutIntCont(5)
integer :: OutRhsCont(5)
integer :: RhsIntCont(5)
logical :: outerCont
real(realk) :: factor
logical :: RHScontraction
END TYPE contractionrule

type intbatch
integer,pointer :: BATCH(:)
end type intbatch


TYPE INTEGRALOUTPUT
REAL(REALK)          :: Result
REAL(REALK),pointer  :: ResultMat(:,:,:,:,:)
type(lstensor)       :: ResultMat2
type(lstensor)       :: RHScont
REAL(REALK),pointer  :: GRAD(:,:)
LOGICAL              :: doGRAD
Integer              :: ndim(5)
! buffer 
LOGICAL              :: USEBUFMM    ! flag for using/not using a buffer 
                                    !to write multipole moments to file
INTEGER,pointer      :: IBUF(:,:)   ! integer buffer
REAL(REALK),pointer  :: RBUF(:,:)   ! real buffer
REAL(REALK),pointer  :: NBUF(:,:)   ! real buffer for nuclear position information
INTEGER              :: MMBUFLEN    ! length of the integer buffer is : MMBUFLEN*MAXBUFI;  
                                    !        of the real buffer MMBUFLEN*MAXBUFR;
                                    !        of the nuclear position buffer  MMBUFLEN*MAXBUFN
INTEGER              :: MAXBUFI,MAXBUFR,MAXBUFN
INTEGER              :: IBUFI       ! counter for the real and integer buffer
INTEGER              :: IBUFN       ! counter for the nuclear buffer
INTEGER              :: LUITNM, LUITNMR ! logical units for the integer and the real buffer files
! end buffer
logical              :: docontrule
type(contractionrule),pointer :: contrule(:)
integer              :: ncontrule
type(intbatch),pointer ::  CnPrimitivesA(:)
type(intbatch),pointer ::  CnPrimitivesB(:)
type(intbatch),pointer ::  startBatchA(:)
type(intbatch),pointer ::  startBatchB(:)
type(intbatch),pointer ::  ContractedBatchA(:)
type(intbatch),pointer ::  ContractedBatchB(:)
END TYPE INTEGRALOUTPUT

!!$TYPE OUTPUTSPEC
!!$LOGICAL :: FULLoutput
!!$LOGICAL :: MATRIXoutput
!!$LOGICAL :: LSTENSORoutput
!!$!----------------------------------------------------
!!$REAL(REALK),pointer :: FULLresult(:,:,:,:,:)
!!$INTEGER             :: FULLdim1
!!$INTEGER             :: FULLdim2
!!$INTEGER             :: FULLdim3
!!$INTEGER             :: FULLdim4
!!$INTEGER             :: FULLdim5
!!$!----------------------------------------------------
!!$TYPE(MATRIX),pointer :: MATRIXresult(:)
!!$INTEGER             :: nRmat
!!$!----------------------------------------------------
!!$!----------------------------------------------------
!!$integer :: lsdim(5)
!!$END TYPE OUTPUTSPEC

CONTAINS
!> \brief initialise the integral output structure
!> \author T. Kjaergaard
!> \date 2010
!> \param IntOut the integraloutput structure to be initialised
!> \param dim1 size og dimension 1
!> \param dim2 size og dimension 2
!> \param dim3 size og dimension 3
!> \param dim4 size og dimension 4
!> \param dim5 size og dimension 5
SUBROUTINE initIntegralOutput(IntOut,dim1,dim2,dim3,dim4,dim5)
implicit none
TYPE(INTEGRALOUTPUT) :: IntOut
INTEGER              :: dim1,dim2,dim3,dim4,dim5

CALL initIntegralOutputDims(IntOut,dim1,dim2,dim3,dim4,dim5)
IntOut%docontrule = .false.
IntOut%ncontrule = 0
NULLIFY(IntOut%ResultMat)
ALLOCATE(IntOut%ResultMat(dim1,dim2,dim3,dim4,dim5))
CALL LS_DZERO(IntOut%ResultMat,dim1*dim2*dim3*dim4*dim5)

END SUBROUTINE initIntegralOutput

!> \brief set the dimensions of the integral output structure
!> \author T. Kjaergaard
!> \date 2010
!> \param IntOut the integraloutput structure to be initialised
!> \param dim1 size og dimension 1
!> \param dim2 size og dimension 2
!> \param dim3 size og dimension 3
!> \param dim4 size og dimension 4
!> \param dim5 size og dimension 5
SUBROUTINE initIntegralOutputDims(IntOut,dim1,dim2,dim3,dim4,dim5)
implicit none
TYPE(INTEGRALOUTPUT) :: IntOut
INTEGER              :: dim1,dim2,dim3,dim4,dim5

IntOut%docontrule = .false.
IntOut%ncontrule = 0
IntOut%ndim(1) = dim1
IntOut%ndim(2) = dim2
IntOut%ndim(3) = dim3
IntOut%ndim(4) = dim4
IntOut%ndim(5) = dim5

END SUBROUTINE initIntegralOutputDims

!!$SUBROUTINE initOutputSpec(Out)
!!$implicit none
!!$TYPE(OUTPUTSPEC) :: Out
!!$
!!$Out%FULLoutput = .FALSE.
!!$Out%MATRIXoutput = .FALSE.
!!$Out%LSTENSORoutput = .FALSE.
!!$NULLIFY(FULLresult)
!!$Out%FULLdim1 = 1
!!$Out%FULLdim2 = 1
!!$Out%FULLdim3 = 1
!!$Out%FULLdim4 = 1
!!$Out%FULLdim5 = 1
!!$NULLIFY(MATRIXresult)
!!$Out%nRmat = 1
!!$
!!$END SUBROUTINE initOutputSpec

!> \brief add a contraction rule
!> \author T. Kjaergaard
!> \date 2010
!>
!>  call add_contrule(setting%output,1,2,3,4,5,0,0,0,0,0,0,0,0,0,0,A,.false.)
!>  would give the OUTPUT(n1,n2,n3,n4,n5) as output
!>
!>  call add_contrule(setting%output,1,4,2,3,5,0,0,0,0,0,0,0,0,0,0,A,.false.)
!>  will change the order of the output and produce
!>  OUTPUT(n1,n4,n2,n3,n5)
!>
!>  call add_contrule(output,2,4,0,0,5,0,0,0,0,0,1,3,0,0,0,A,.true.)
!>  will contract the integral with a rhs density
!>  OUTPUT(n2,n4,1,1,n5) = sum_{n1,n3} A*INTEGRAL(n1,n2,n3,n4,n5) * DMAT(n1,n3) 
!>
!>  call add_contrule(setting%output,1,2,0,0,5,0,0,0,0,0,3,4,0,0,0,A,.true.)
!>  will contract the integral with a rhs density
!>  OUTPUT(n1,n2,1,1,n5) = sum_{n3,n4} A*INTEGRAL(n1,n2,n3,n4,n5) * DMAT(n3,n4) 
!>
!> \param Output the integraloutput structure
!> \param OI1 the index for the first output integral dim
!> \param OI2 the index for the first output integral dim
!> \param OI3 the index for the first output integral dim
!> \param OI4 the index for the first output integral dim
!> \param OI5 the index for the first output integral dim
!> \param OR1 the index for the RHS outer contraction
!> \param OR2 the index for the RHS outer contraction
!> \param OR3 the index for the RHS outer contraction
!> \param OR4 the index for the RHS outer contraction
!> \param OR5 the index for the RHS outer contraction
!> \param RI1 the index for the RHS inner contraction
!> \param RI2 the index for the RHS inner contraction
!> \param RI3 the index for the RHS inner contraction
!> \param RI4 the index for the RHS inner contraction
!> \param RI5 the index for the RHS inner contraction
!> \param factor a factor to include in the DMAT contraction
!> \param RHScontraction flag to do a contraction with rhs density
subroutine add_contrule(output,OI1,OI2,OI3,OI4,OI5,OR1,OR2,OR3,OR4,OR5,RI1,RI2,RI3,RI4,RI5,factor,RHScontraction)
implicit none
type(INTEGRALOUTPUT),intent(inout)    :: output
integer,intent(in) :: OI1,OI2,OI3,OI4,OI5,OR1,OR2,OR3,OR4,OR5,RI1,RI2,RI3,RI4,RI5
logical,intent(in)     :: RHScontraction
real(realk),intent(in) :: factor
!
type(contractionrule),pointer  :: tmpcontrule(:)
Integer :: I

Output%docontrule = .true.
if(output%ncontrule.GT.0)then
   nullify(tmpcontrule)
   allocate(tmpcontrule(output%ncontrule))
   do I=1,output%ncontrule
      tmpcontrule(I)=output%contrule(I)
   enddo
   deallocate(output%contrule)
   nullify(output%contrule)
   allocate(output%contrule(output%ncontrule+1))
   do I=1,output%ncontrule
      output%contrule(I)=tmpcontrule(I)
   enddo
   deallocate(tmpcontrule)
   nullify(tmpcontrule)
else
   nullify(output%contrule)
   allocate(output%contrule(output%ncontrule+1))
endif
output%ncontrule = output%ncontrule + 1
I = output%ncontrule

output%contrule(I)%RHScontraction = RHScontraction
output%contrule(I)%OutIntCont(1) = OI1
output%contrule(I)%OutIntCont(2) = OI2
output%contrule(I)%OutIntCont(3) = OI3
output%contrule(I)%OutIntCont(4) = OI4
output%contrule(I)%OutIntCont(5) = OI5

output%contrule(I)%OutRhsCont(1) = OR1
output%contrule(I)%OutRhsCont(2) = OR2
output%contrule(I)%OutRhsCont(3) = OR3
output%contrule(I)%OutRhsCont(4) = OR4
output%contrule(I)%OutRhsCont(5) = OR5

output%contrule(I)%RhsIntCont(1) = RI1
output%contrule(I)%RhsIntCont(2) = RI2
output%contrule(I)%RhsIntCont(3) = RI3
output%contrule(I)%RhsIntCont(4) = RI4
output%contrule(I)%RhsIntCont(5) = RI5

output%contrule(I)%OuterCont = .False.
output%contrule(I)%factor = Factor
!print*,'Factor chosen=',output%contrule(I)%factor
IF(OR1+OR2+OR3+OR4+OR5 .GT. 0 ) output%contrule(I)%OuterCont = .True.

end subroutine add_contrule

!> \brief free the contraction rule
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param Output the integraloutput structure
subroutine free_contrule(output)
implicit none
type(INTEGRALOUTPUT),intent(inout)    :: output

if(output%ncontrule.GT.0)then
   deallocate(output%contrule)
   nullify(output%contrule)
endif
output%ncontrule=0

end subroutine free_contrule

!> \brief init the lstensor for the primitive gab matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param Output contains output specs for integral storage
!> \param TENSOR the output lstensor
!> \param pAO1 the 1. primitive Atomic orbital
!> \param pAO2 the 2. primitive Atomic orbital
!> \param cAO1 the 1. contracted Atomic orbital
!> \param cAO2 the 2. contracted Atomic orbital
!> \param lupri the logical unit number of the output file
SUBROUTINE init_primGablstensor(Output,TENSOR,psAO1,psAO2,csAO1,csAO2,lupri)
implicit none
Type(IntegralOutput) :: Output
TYPE(LSTENSOR)       :: TENSOR
TYPE(AOITEM),target  :: csAO1,csAO2,psAO1,psAO2
TYPE(AOITEM),pointer :: cAOT1,cAOT2,pAOT1,pAOT2
INTEGER            :: nbast1,nbast2,lupri!,n1,n2
logical   :: useAO1,useAO2
!
TYPE(AOITEM),target :: AOT
!REAL(REALK) :: MAXELM
INTEGER :: nElms
! BELONGINING TO AOT1 
INTEGER :: natom1,nbastI,nbatI,Ibat,Iangmom,Iatom,IORB
INTEGER :: nContA,nOrbCompA,sA,nAngA,AOT1batch,batAOT1,nmat
! BELONGINING TO AOT2 
INTEGER :: natom2,nbastJ,nbatJ,Jbat,Jangmom,Jatom,JORB
INTEGER :: nContB,nOrbCompB,sB,nAngB,AOT2batch,batAOT2
! COMMON
INTEGER :: DIM,IELM,IMAT2,I,npbatI,npbatJ,npbatK,npbatL,mark_pIbat,mark_pJbat,mark_pKbat
INTEGER :: mark_pLbat,pIbat,pJbat,pKbat,pLbat,ip,jp,kp,lp,dima,dimb,dimc,dimd,na,nb,nc,nd
INTEGER :: lock_startD,lock_startC,lock_startB,lock_startA,pIELM,startD,startA,startB,startC
INTEGER :: IA,IB,IC,ID,nP1,nP2,pIbat_start,pJbat_start,oldIbat,oldJbat,cIbat,cJbat
LOGICAL :: ERROR
real(realk) :: maxelm,DIFF
!real(realk) :: RefscreenMat(n1,n2)
integer  :: atomA,atomB,index,iOrbCompA,iOrbCompB,J,Abatch,Bbatch,iP1S,iP2S,ip1,ip2,pbatAOT1,pbatAOT2
real(realk),pointer :: temp(:,:)
integer,pointer  :: pIbat2Ibat(:),pJbat2Jbat(:),pIELMs2(:)
integer(kind=long) :: nmemsize

TENSOR%primCStensor=.TRUE.
TENSOR%gradienttensor = .FALSE.
call SET_EMPTY_AO(AOT)

IF(csAO1%empty)THEN
   cAOT1 => AOT
   pAOT1 => AOT
ELSE
   cAOT1 => csAO1
   pAOT1 => psAO1
ENDIF

IF(csAO2%empty)THEN
   cAOT2 => AOT
   pAOT2 => AOT
ELSE
   cAOT2 => csAO2
   pAOT2 => psAO2
ENDIF

nbast1 = cAOT1%nbast
nbast2 = cAOT2%nbast
natom1 = cAOT1%natoms
natom2 = cAOT2%natoms
NULLIFY(TENSOR%LSAO)
NULLIFY(TENSOR%INDEX)
ALLOCATE(TENSOR%LSAO(natom1*natom2))
ALLOCATE(TENSOR%INDEX(natom1,natom2,1,1))
TENSOR%INDEX = 0 !if 0 lsaotensor not allocated 
TENSOR%natom1 = natom1
TENSOR%natom2 = natom2
TENSOR%natom3 = 1
TENSOR%natom4 = 1
TENSOR%nbast1 = nbast1
TENSOR%nbast2 = nbast2
TENSOR%nbast3 = 1
TENSOR%nbast4 = 1
TENSOR%nmat = 1
AOT1batch=0
I = 0
DO Iatom = 1,natom1
 nbatI = cAOT1%ATOMICnBatch(IATOM)
 AOT2batch=0
 DO Jatom = 1,natom2
  nbatJ = cAOT2%ATOMICnBatch(JATOM)
  I=I+1
  TENSOR%INDEX(IATOM,JATOM,1,1) = I
  NULLIFY(TENSOR%LSAO(I)%BATCH)
  ALLOCATE(TENSOR%LSAO(I)%BATCH(nbatI,nbatJ,1,1))
  TENSOR%LSAO(I)%ALLOC = .TRUE.
  TENSOR%LSAO(I)%ATOM1 = Iatom
  TENSOR%LSAO(I)%ATOM2 = Jatom
  TENSOR%LSAO(I)%ATOM3 = 1
  TENSOR%LSAO(I)%ATOM4 = 1
  TENSOR%LSAO(I)%nAOBATCH1 = nbatI
  TENSOR%LSAO(I)%nAOBATCH2 = nbatJ
  TENSOR%LSAO(I)%nAOBATCH3 = 1
  TENSOR%LSAO(I)%nAOBATCH4 = 1

  npbatI = pAOT1%ATOMICnBatch(IATOM)
  npbatJ = pAOT2%ATOMICnBatch(JATOM)

  batAOT1 = AOT1batch
  DO Ibat = 1,nbatI
   batAOT1 = batAOT1+1 
   nP1 = cAOT1%BATCH(batAOT1)%nPrimitives
   batAOT2 = AOT2batch
   DO Jbat = 1,nbatJ
    batAOT2 = batAOT2+1 
    nP2 = cAOT2%BATCH(batAOT2)%nPrimitives
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomA = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomB = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomC = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nAngmomD = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContA(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompA(1) = nP1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbA(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContB(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompB(1) = nP2
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbB(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContC(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompC(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbC(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nContD(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nOrbCompD(1) = 1 
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%startOrbD(1) = 1
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nelmE = 1
    DIM = nP1*nP2
    NULLIFY(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms)
    ALLOCATE(TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms(DIM))
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%nelms = DIM
    TENSOR%LSAO(I)%BATCH(Ibat,Jbat,1,1)%elms = 0.d0
   ENDDO 
  ENDDO
  AOT2batch = AOT2batch + nbatJ
 ENDDO
 AOT1batch = AOT1batch + nbatI
ENDDO
TENSOR%nLSAO = I

CALL build_primassistarrays(OUTPUT%CnPrimitivesA,Output%startBatchA,Output%ContractedBatchA,natom1,cAOT1)

CALL build_primassistarrays(OUTPUT%CnPrimitivesB,Output%startBatchB,Output%ContractedBatchB,natom2,cAOT2)

!!$DO Iatom = 1,natom1
!!$   WRITE(lupri,*)'IATOM:',IATOM
!!$   DO pIbat = 1,pAOT1%ATOMICnBatch(IATOM)
!!$      WRITE(lupri,*)'THE ContractedBatchA, StartBatchA     nr pbatches',pAOT1%ATOMICnBatch(IATOM)
!!$      write(lupri,*)Output%ContractedBatchA(IATOM)%BATCH(pIbat),&
!!$           &Output%startBatchA(IATOM)%BATCH(pIbat)
!!$   enddo
!!$   WRITE(lupri,*)'THE nPrimitivesA                      nr batches',cAOT1%ATOMICnBatch(IATOM)
!!$   do Ibat=1,cAOT1%ATOMICnBatch(IATOM)
!!$      write(lupri,*)Output%CnPrimitivesA(IATOM)%BATCH(Ibat)
!!$   enddo
!!$enddo
call FREE_EMPTY_AO(AOT)
!print*,'end init_primGablstensor'

Call Determine_lstensor_memory(tensor,nmemsize)
call add_mem_to_global(nmemsize)

END SUBROUTINE Init_primGablstensor

!> \brief build som arrays needed to place primitive screening integrals in proper place in lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param CnPrimitivesA number of primitives in the contracted basis
!> \param startBatchA starting index in contracted basis
!> \param ContractedBatchA contracted batch in contracted basis
!> \param natom the number of atoms in the atomic orbitals
!> \param cAOT a contracted Atomic orbital
SUBROUTINE build_primassistarrays(CnPrimitivesA,startBatchA,ContractedBatchA,&
     &natom,cAOT)
implicit none
type(intbatch),pointer :: CnPrimitivesA(:),startBatchA(:),ContractedBatchA(:)
TYPE(AOITEM)   :: cAOT
integer :: natom
!
integer :: AOTbatch, nbatI,pIbat,iP1,nP1,Ibat,Iatom

AOTbatch=0
ALLOCATE(CnPrimitivesA(natom))
ALLOCATE(startBatchA(natom))
ALLOCATE(ContractedBatchA(natom))
DO Iatom = 1,natom
 nbatI = cAOT%ATOMICnBatch(IATOM)
 call mem_alloc(CnPrimitivesA(IATOM)%BATCH,nbatI)
 pIbat = 0
 DO Ibat = 1,nbatI
    nP1 = cAOT%BATCH(AOTbatch+Ibat)%nPrimitives
    do iP1=1,nP1
       pIbat = pIbat+1
    enddo
 ENDDO
 call mem_alloc(startBatchA(IATOM)%BATCH,pIbat)
 call mem_alloc(ContractedBatchA(IATOM)%BATCH,pIbat)
 pIbat = 0
 DO Ibat = 1,nbatI
    nP1 = cAOT%BATCH(AOTbatch+Ibat)%nPrimitives
    CnPrimitivesA(IATOM)%BATCH(Ibat)=nP1
    do iP1=1,nP1
       pIbat = pIbat+1
       startBatchA(IATOM)%BATCH(pIbat) = iP1-1!pIbat-1
       ContractedBatchA(IATOM)%BATCH(pIbat) = Ibat
    enddo
 ENDDO
 AOTbatch = AOTbatch + nbatI
ENDDO
END SUBROUTINE BUILD_PRIMASSISTARRAYS

!> \brief free the prim assist arrays
!> \author T. Kjaergaard
!> \date 2010
!> \param Output contains output specs for integral storage
SUBROUTINE free_primassistarrays(Output,TENSOR)
implicit none
Type(IntegralOutput) :: Output
TYPE(LSTENSOR)       :: TENSOR
!
integer :: iatom

IF(ASSOCIATED(OUTPUT%CnPrimitivesA))THEN
   DO IATOM=1,TENSOR%natom1
      call mem_dealloc(OUTPUT%CnPrimitivesA(IATOM)%BATCH)
      call mem_dealloc(OUTPUT%startBatchA(IATOM)%BATCH)
      call mem_dealloc(OUTPUT%ContractedBatchA(IATOM)%BATCH)
   ENDDO
   DEALLOCATE(OUTPUT%CnPrimitivesA)
   DEALLOCATE(OUTPUT%startBatchA)
   DEALLOCATE(OUTPUT%ContractedBatchA)
ENDIF
IF(ASSOCIATED(OUTPUT%CnPrimitivesB))THEN
   DO IATOM=1,TENSOR%natom2
      call mem_dealloc(OUTPUT%CnPrimitivesB(IATOM)%BATCH)
      call mem_dealloc(OUTPUT%startBatchB(IATOM)%BATCH)
      call mem_dealloc(OUTPUT%ContractedBatchB(IATOM)%BATCH)
   ENDDO
   DEALLOCATE(OUTPUT%CnPrimitivesB)
   DEALLOCATE(OUTPUT%startBatchB)
   DEALLOCATE(OUTPUT%ContractedBatchB)
ENDIF

END SUBROUTINE FREE_PRIMASSISTARRAYS

end MODULE integraloutput_type
