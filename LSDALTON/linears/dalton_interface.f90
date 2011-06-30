!===========================================================================
! dal_interface module contains all the routines that allow for
! exchanging information with the remaining part of dalton which
! may have different memory layout etc.

MODULE dal_interface
   use LStiming
   USE Matrix_operations

INTERFACE di_GET_GbDs
  MODULE PROCEDURE di_GET_GbDsSingle, di_GET_GbDsArray
END INTERFACE

INTERFACE di_GET_GbDs_LSDALTON
  MODULE PROCEDURE di_GET_GbDsSingle_LSDALTON, di_GET_GbDsArray_LSDALTON
END INTERFACE

   PUBLIC::&
        & di_debug_4center_eri, &
        & di_debug_4center, &
        & di_decpacked, &
        & di_debug_ccfragment, &
        & di_get_overlap_and_H1, &
        & di_get_fock_LSDALTON, &
        & di_GET_GbDs, &
        & di_GET_GbDs_lsdalton

   PRIVATE

   ! local pointer to H1 matrix allocated in linscf() subroutine. This allows
   ! access to H1 within this module, avoiding re-reading H1 from disk and
   ! avoiding extensive modifications to existing subroutine interfaces.
   type(MATRIX), pointer, save  :: lH1
CONTAINS
   SUBROUTINE di_debug_4center_eri(lu_pri,lu_err,ls,nbast)
     use TYPEDEF
     IMPLICIT NONE
     integer,intent(in)      :: lu_pri,lu_err,nbast
     type(lsitem),intent(inout) :: ls
     !
     real(realk),pointer   :: integrals(:,:,:,:)
     integer :: iB,iA,iC,iD
     !THIS IS A SPECIAL HARDCODED DEBUG ROUTINE FOR WATER (TESTCASE LS_DEBUG4CENTER)
     call mem_alloc(integrals,nbast,nbast,nbast,nbast)
     call II_get_4center_eri(LU_PRI,LU_ERR,ls%setting,integrals,nbast,nbast,nbast,nbast)
     WRITE(lu_pri,*)'di_debug_4center_eri: some random 4 center integrals'
     !SAME ATOM ON ALL 4
     iA = 1; iB = 1; iC = 1; iD = 1
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 2; iC = 4; iD = 4
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 3; iB = 3; iC = 9; iD = 9
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 3; iC = 6; iD = 9
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     !SAME ATOM ON 1,2,3 and new on 4.
     iA = 1; iB = 1; iC = 1; iD = 11
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 2; iC = 4; iD = 11
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 3; iB = 3; iC = 9; iD = 12
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 3; iC = 6; iD = 12
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     !SAME ATOM ON LHS and same on RHS.
     iA = 1; iB = 1; iC = 11; iD = 11
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 2; iC = 11; iD = 11
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 3; iB = 3; iC = 11; iD = 11
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 3; iC = 11; iD = 11
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     !SAME ATOM ON LHS and different on RHS.
     iA = 1; iB = 1; iC = 13; iD = 11
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 2; iC = 13; iD = 11
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 3; iB = 3; iC = 10; iD = 12
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 3; iC = 10; iD = 12
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     !SAME ATOM ON 1 and same on 2,3,4
     iA = 10; iB = 10; iC = 10; iD = 10
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 11; iB = 10; iC = 10; iD = 10
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 12; iB = 10; iC = 10; iD = 10
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 13; iB = 10; iC = 10; iD = 10
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 10; iB = 1; iC = 1; iD = 1
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 11; iB = 1; iC = 1; iD = 1
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 12; iB = 1; iC = 1; iD = 1
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 13; iB = 1; iC = 1; iD = 1
     WRITE(lu_pri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     call mem_dealloc(integrals)

   END SUBROUTINE di_debug_4center_eri

   SUBROUTINE di_debug_4center(lu_pri,lu_err,ls,nbast,D)
   use TYPEDEF
      IMPLICIT NONE
      TYPE batchTOorb 
         INTEGER,pointer :: orbindex(:)
         INTEGER :: norbindex
      END TYPE batchTOorb
      TYPE(Matrix),intent(in) :: D
      integer,intent(in)      :: lu_pri,lu_err,nbast
      type(lsitem),intent(inout) :: ls
      !
      TYPE(Matrix)  :: tempm1,tempm2,tempm3
      real(realk),pointer   :: integrals(:,:,:,:)
      real(realk) :: CoulombFactor,FAB,int
      integer,pointer :: orb2batch(:),batchdim(:)
      type(batchTOorb),pointer :: batch2orb(:)
      integer :: Abatch,Bbatch,IORB,JORB,dimA,dimB,J,nbatches
      integer :: iB,iA,iC,iD,offset,offset2,BORB,AORB,iprint

      CoulombFactor = 2.d0
      !build the coulomb matrix 
      call mat_init(tempm1,nbast,nbast)
      call mat_init(tempm2,nbast,nbast)
      call mat_init(tempm3,nbast,nbast)
      CALL II_get_coulomb_mat(lu_pri,lu_err,ls%setting,D,tempm1,1)

      allocate(orb2batch(nbast))
      call LSlib_getBatchOrbitalInfo(orb2Batch,nBast,nBatches,lu_pri,lu_err)

      allocate(batchdim(nbatches))
      batchdim = 0
      DO IORB=1,nbast
         batchdim(orb2batch(IORB)) = batchdim(orb2batch(IORB))+1 
      ENDDO

      allocate(batch2orb(nbatches))
      DO Abatch=1,nbatches
         allocate(batch2orb(Abatch)%orbindex(batchdim(Abatch)))
         batch2orb(Abatch)%orbindex = 0
         batch2orb(Abatch)%norbindex = 0
      ENDDO
      DO IORB=1,nbast
         Abatch = orb2batch(IORB)
         batch2orb(Abatch)%norbindex = batch2orb(Abatch)%norbindex+1 
         J = batch2orb(Abatch)%norbindex
         batch2orb(Abatch)%orbindex(J) = IORB
      ENDDO
      DO Bbatch=1,nbatches
         dimB = batchdim(Bbatch)
         DO Abatch=1,nbatches
            dimA = batchdim(Abatch)
            nullify(integrals)
            allocate(integrals(dimA,dimB,nbast,nbast))
!            integrals=0.d0
            CALL LSlib_get_ABresctricted_4CenterEri(integrals,Abatch,Bbatch,dimA,&
                 &dimB,nbast,lu_pri,lu_err)            
            DO iB = 1,dimB
               !translate to orbitalindex: B orbital
               BORB = batch2orb(Bbatch)%orbindex(iB)
               offset2 = (BORB-1)*nbast
               DO iA = 1,dimA
                  !translate to orbitalindex: B orbital
                  AORB = batch2orb(Abatch)%orbindex(iA) 
                  FAB =0.d0
                  DO id=1,nbast
                     offset=(id-1)*nbast
                     DO ic=1,nbast
                        FAB=FAB+D%elms(ic+offset)*CoulombFactor*integrals(ia,ib,ic,id)
                     ENDDO
                  ENDDO
                  tempm2%elms(AORB+offset2) = FAB
               ENDDO
            ENDDO
            deallocate(integrals)
            nullify(integrals)
!            call mem_dealloc(integrals)
         ENDDO
      ENDDO
      DO Abatch=1,nbatches
         deallocate(batch2orb(Abatch)%orbindex)
      ENDDO
      deallocate(batch2orb)
      deallocate(batchdim)
      deallocate(orb2batch)
      call mat_add(1.d0,tempm1,-1.d0,tempm2,tempm3)
      write(lu_pri,*) 'QQQ DI_DEBUG_4CENTER STD ',mat_trab(tempm1,tempm1)
      write(lu_pri,*) 'QQQ DI_DEBUG_4CENTER CHOL',mat_trab(tempm2,tempm2)
      IF(ABS(mat_trab(tempm3,tempm3)).LE.1.D-15)THEN
         write(lu_pri,*)'QQQ SUCCESFUL CHOLESKY TEST'
      ELSE
         CALL lsQUIT('CHOLESKY TEST TEST',lu_pri)
      ENDIF
      call mat_free(tempm1)
      call mat_free(tempm2)
      call mat_free(tempm3)

    END SUBROUTINE DI_DEBUG_4CENTER

   SUBROUTINE di_decpacked(lu_pri,lu_err,ls,nbast,D)
   use TYPEDEF
      IMPLICIT NONE
   TYPE batchTOorb 
      INTEGER,pointer :: orbindex(:)
      INTEGER :: norbindex
   END TYPE batchTOorb
      TYPE(Matrix),intent(in) :: D
      integer,intent(in)      :: lu_pri,lu_err,nbast
      type(lsitem),intent(inout) :: ls
      !
      TYPE(Matrix)  :: K,Kdec,tempm3,J, Jdec
      real(realk),pointer   :: integrals(:,:,:,:)
      real(realk) :: ExchangeFactor,FAB,int,Dmatelm,CoulombFactor
      integer,pointer :: orb2batch(:),batchdim(:)
      type(batchTOorb),pointer :: batch2orb(:)
      integer :: Cbatch,Bbatch,IORB,JORB,dimC,dimB,JK
      integer :: nbatches,Dbatch,dimD,iDbatch
      integer :: iB,iA,iC,iD,offset,offset2,BORB,CORB,iprint,ibbatch,icbatch

      ExchangeFactor = -1.d0
      !build the coulomb matrix 
      call mat_init(K,nbast,nbast)
      call mat_init(Kdec,nbast,nbast)
      call mat_init(tempm3,nbast,nbast)
      call mat_zero(K)
      CALL II_get_exchange_mat(lu_pri,lu_err,ls%setting,D,1,.TRUE.,K)

      call mem_alloc(orb2batch,nbast)
      call LSlib_getBatchOrbitalInfo(orb2Batch,nBast,nBatches,lu_pri,lu_err)

      allocate(batchdim(nbatches))
      batchdim = 0
      DO IORB=1,nbast
         batchdim(orb2batch(IORB)) = batchdim(orb2batch(IORB))+1 
      ENDDO

      allocate(batch2orb(nbatches))
      DO Cbatch=1,nbatches
         allocate(batch2orb(Cbatch)%orbindex(batchdim(Cbatch)))
         batch2orb(Cbatch)%orbindex = 0
         batch2orb(Cbatch)%norbindex = 0
      ENDDO
      DO IORB=1,nbast
         Cbatch = orb2batch(IORB)
         batch2orb(Cbatch)%norbindex = batch2orb(Cbatch)%norbindex+1 
         JK = batch2orb(Cbatch)%norbindex
         batch2orb(Cbatch)%orbindex(JK) = IORB
      ENDDO

      call mat_zero(Kdec)
      DO Cbatch=1,nbatches
         dimC = batchdim(Cbatch)
         DO Bbatch=1,nbatches
            dimB = batchdim(Bbatch)
            nullify(integrals)
            allocate(integrals(nbast,nbast,dimB,dimC))
            call II_GET_DECPACKED4CENTER_K_ERI(LU_PRI,LU_ERR,ls%SETTING,&
                 &integrals,Bbatch,Cbatch,nbast,dimB,dimC)
            DO iCbatch = 1,dimC
               !translate to orbitalindex: C orbital
               iC = batch2orb(Cbatch)%orbindex(iCbatch)
               offset=(iC-1)*nbast
               DO iBbatch = 1,dimB
                  !translate to orbitalindex: B orbital
                  iB = batch2orb(Bbatch)%orbindex(iBbatch) 
                  Dmatelm = D%elms(iB+offset)*exchangeFactor
                  DO id=1,nbast
                     offset2 = (iD-1)*nbast
                     DO ia=1,nbast
                        Kdec%elms(iA+offset2)=Kdec%elms(ia+offset2)+integrals(ia,id,iBbatch,iCbatch)*Dmatelm
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            deallocate(integrals)
            nullify(integrals)
         ENDDO
      ENDDO
      call mat_add(1.d0,Kdec,-1.d0,K,tempm3)
      write(lu_pri,*) 'QQQ DI_DEBUG_DECPACK K STD    ',mat_trab(K,K)
      write(lu_pri,*) 'QQQ DI_DEBUG_DECPACK K DECPACK',mat_trab(Kdec,Kdec)
      write(lu_pri,*) 'QQQ DIFF',ABS(mat_trab(tempm3,tempm3))
      IF(ABS(mat_trab(tempm3,tempm3)).LE.1.D-15)THEN
         write(lu_pri,*)'QQQ SUCCESFUL DECPACK K TEST'
      ELSE
         CALL lsQUIT('DECPACKED K TEST FAILED',lu_pri)
      ENDIF
      call mat_free(K)
      call mat_free(Kdec)
      CoulombFactor = 2.d0
      call mat_init(J,nbast,nbast)
      call mat_init(Jdec,nbast,nbast)
      call mat_zero(J)

      CALL II_get_coulomb_mat(lu_pri,lu_err,ls%setting,D,J,1)

      call mat_zero(Jdec)
      DO Dbatch=1,nbatches
         dimD = batchdim(Dbatch)
         DO Cbatch=1,nbatches
            dimC = batchdim(Cbatch)
            nullify(integrals)
            allocate(integrals(nbast,nbast,dimC,dimD))
            call II_GET_DECPACKED4CENTER_J_ERI(LU_PRI,LU_ERR,ls%SETTING,&
                 &integrals,Cbatch,Dbatch,nbast,dimC,dimD)
            DO iDbatch = 1,dimD
               !translate to orbitalindex: D orbital
               iD = batch2orb(Dbatch)%orbindex(iDbatch)
               offset=(iD-1)*nbast
               DO iCbatch = 1,dimC
                  !translate to orbitalindex: B orbital
                  iC = batch2orb(Cbatch)%orbindex(iCbatch) 
                  Dmatelm = D%elms(iC+offset)*CoulombFactor
                  DO ib=1,nbast
                     offset2 = (ib-1)*nbast
                     DO ia=1,nbast
                        Jdec%elms(iA+offset2)=Jdec%elms(iA+offset2)+integrals(ia,ib,iCbatch,iDbatch)*Dmatelm
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            deallocate(integrals)
            nullify(integrals)
         ENDDO
      ENDDO
      call mat_add(1.d0,Jdec,-1.d0,J,tempm3)
      write(lu_pri,*) 'QQQ DI_DEBUG_DECPACK_J  STD    ',mat_trab(J,J)
      write(lu_pri,*) 'QQQ DI_DEBUG_DECPACK_J  DECPACK',mat_trab(Jdec,Jdec)
      write(lu_pri,*) 'QQQ DIFF',ABS(mat_trab(tempm3,tempm3))
      IF(ABS(mat_trab(tempm3,tempm3)).LE.1.D-15)THEN
         write(lu_pri,*)'QQQ SUCCESFUL DECPACK J TEST'
      ELSE
         CALL lsQUIT('DECPACKED J TEST FAILED',lu_pri)
      ENDIF
      call mem_dealloc(orb2batch)
      call mat_free(J)
      call mat_free(Jdec)
      call mat_free(tempm3)

    END SUBROUTINE DI_DECPACKED

   SUBROUTINE di_debug_ccfragment(lu_pri,lu_err,ls,nbast,D,CMO,nocc)
   use TYPEDEF
   use AO_Type
   use BUILDAOBATCH
   use DALTONINFO
      IMPLICIT NONE
      TYPE(Matrix),intent(in) :: D,CMO
      integer,intent(in)      :: lu_pri,lu_err,nbast,nocc
      type(lsitem),intent(inout) :: ls
      !
      type(lsitem),pointer :: fragment_ls(:)
      type(AOITEM) :: AO
      TYPE(Matrix)  :: tempm1,tempm2,tempm3,tempm4,D2,refF
      real(realk),pointer   :: densPtmp(:,:)
      integer,pointer :: orb2batch(:),orb2atom(:),ATOMS(:)
      logical,pointer :: LATOMS(:)
      integer :: iatom,nbast1,nbast2,I,J,dim1,dim2,dim3,dim4,igamma
      integer :: delta,iorb,natoms,nBatches,natom2,iprint
      real(realk) :: MAXCMO
      iprint=0
      WRITE(LU_PRI,*)'DEBUGGING THE CCFRAGMENT'
      WRITE(*,*)'DEBUGGING THE CCFRAGMENT'

      !build the coulomb matrix 
      call mat_init(tempm1,nbast,nbast)
      call mat_init(tempm2,nbast,nbast)
      call mat_init(tempm3,nbast,nbast)
      call mat_init(refF,nbast,nbast)
      CALL II_get_coulomb_mat(LU_PRI,LU_ERR,ls%setting,D,refF,1)

      natoms = ls%input%molecule%natoms
      nullify(LATOMS)
      nullify(ATOMS)
      allocate(LATOMS(natoms))
      ALLOCATE(ATOMS(natoms))

      nullify(orb2batch)
      allocate(orb2batch(nbast))
      call LSlib_getBatchOrbitalInfo(orb2Batch,nBast,nBatches,lu_pri,lu_err)
      call BUILD_AO(LU_PRI,ls%setting%SCHEME,0,ls%setting%MOLECULE(1)%p,ls%setting%BASIS(1)%p%REGULAR,AO&
     &,'Default',.false.,.false.)

      call mem_alloc(orb2atom,nbast)
      do iorb = 1,nbast
         orb2atom(iorb) = AO%BATCH(orb2batch(iorb))%atom
      enddo

      dim1=nbast
      dim2=nbast

      nullify(fragment_ls)
      allocate(fragment_ls(nocc))
      call mat_zero(tempm3)
      MAXCMO = 0.d0
      DO iorb = 1,nocc
         do igamma=1,nbast
            MAXCMO = MAX(MAXCMO,CMO%elms(igamma+(iorb-1)*nbast))
         enddo
         LATOMS = .FALSE. 
         do igamma=1,nbast
            IF(ABS(MAXCMO*CMO%elms(igamma+(iorb-1)*nbast)).GT. 1.D-13)THEN
               LATOMS(orb2atom(igamma)) = .TRUE. 
            ENDIF
         enddo
!         print*,'LATOMS',LATOMS
         nbast1 = 0
         do igamma=1,nbast
            IF(LATOMS(orb2atom(igamma)))THEN
               nbast1 = nbast1+1
            ENDIF
         enddo
         dim3=nbast1
         dim4=nbast1
!         print*,'dim3',dim3,'dim4',dim4
         natom2=0
         ATOMS=0
         DO iatom = 1,natoms
            IF(LATOMS(iatom))THEN
               natom2=natom2+1
               ATOMS(natom2)=iatom
            ENDIF
         ENDDO
!         print*,'ATOMS',ATOMS
         nbast1 = 0
         nullify(DensPtmp)
         allocate(DensPtmp(dim3,dim4))
         DensPtmp = 0.d0
         do igamma=1,nbast
            if(LATOMS(orb2atom(igamma)))THEN
               nbast1 = nbast1+1
               nbast2 = 0
               do delta=1,nbast
                  if(LATOMS(orb2atom(delta)))THEN
                     nbast2 = nbast2+1
                     DensPtmp(nbast1,nbast2) = CMO%elms(igamma+(iorb-1)*nbast)*CMO%elms(delta+(iorb-1)*nbast)
                  endif
               enddo
            endif
         enddo
         call build_ccfragmentlsitem(LS,fragment_ls(iorb),ATOMS(1:NATOM2),NATOM2,LU_PRI,IPRINT)
         ! here we set the ab in (ab|cd) to full
         fragment_ls(iorb)%SETTING%MOLECULE(1)%p => LS%SETTING%MOLECULE(1)%p
         fragment_ls(iorb)%SETTING%BASIS(1)%p => LS%SETTING%BASIS(1)%p
         fragment_ls(iorb)%SETTING%FRAGMENT(1)%p => LS%SETTING%FRAGMENT(1)%p
         fragment_ls(iorb)%SETTING%MOLECULE(2)%p => LS%SETTING%MOLECULE(2)%p
         fragment_ls(iorb)%SETTING%BASIS(2)%p => LS%SETTING%BASIS(2)%p
         fragment_ls(iorb)%SETTING%FRAGMENT(2)%p => LS%SETTING%FRAGMENT(2)%p
         fragment_ls(iorb)%SETTING%molID(1)=0
         fragment_ls(iorb)%SETTING%molID(2)=0
         fragment_ls(iorb)%SETTING%molID(3)=iorb
         fragment_ls(iorb)%SETTING%molID(4)=iorb
         fragment_ls(iorb)%SETTING%sameMOL=.TRUE.
         DO I=1,2
            DO J=3,4          
               fragment_ls(iorb)%SETTING%sameMOL(I,J) = .FALSE.
               fragment_ls(iorb)%SETTING%sameMOL(J,I) = .FALSE.
            ENDDO
         ENDDO
         call mat_init(D2,dim3,dim4)
         CALL mat_set_from_full(DensPtmp,1.d0, D2,'DenP')
         deallocate(DensPtmp)
         call mat_zero(tempm1)
         CALL II_get_coulomb_mat(LU_PRI,LU_ERR,fragment_ls(iorb)%setting,D2,tempm1,1)
         call mat_free(D2)

!         allocate(DensPtmp(nbast,nbast))
!         DensPtmp=0.d0
!         do delta=1,nbast
!            do gamma=1,nbast
!               DensPtmp(gamma,delta) = CMO%elms(gamma+(iorb-1)*nbast)*CMO%elms(delta+(iorb-1)*nbast)
!            enddo
!         enddo
!         call mat_init(D2,nbast,nbast)
!         CALL mat_set_from_full(DensPtmp,1.d0, D2,'DenP')
!         call mat_zero(tempm2)
!         CALL II_get_coulomb_mat(LUPRI,LUERR,ls%setting,D2,tempm2)
!         call mat_init(tempm4,nbast,nbast)
!         call mat_add(1.d0,tempm1,-1.d0,tempm2,tempm4)
!         write(lupri,*) 'QQQ DI_DEBUG',iorb,'             :',mat_trab(tempm2,tempm2)
!         write(lupri,*) 'QQQ DI_DEBUG lsfragment',iorb,'  :',mat_trab(tempm1,tempm1)
!         write(lupri,*) 'QQQ DIFF                         :',mat_trab(tempm4,tempm4)
!         IF(ABS(mat_trab(tempm4,tempm4)).LE.1.D-8)THEN
!            write(lupri,*)'QQQ SUCCESFUL'
!         ELSE
!            WRITE(lupri,*)'ZZ THE FULL DMAT IORB =',iorb
!            call mat_print(D2,1,nbast,1,nbast,lupri)
!            WRITE(lupri,*)'THE FRAGMENT FOCK'
!            call mat_print(tempm1,1,nbast,1,nbast,lupri)
!            WRITE(lupri,*)'THE LS FOCK'
!            call mat_print(tempm2,1,nbast,1,nbast,lupri)
!            WRITE(lupri,*)'THE DIFF'
!            call mat_print(tempm4,1,nbast,1,nbast,lupri)            
!            WRITE(LUPRI,*)'ZZZ LS SETTING'
!            call PRINT_LSSETTING(ls%setting,LUPRI)
!            WRITE(LUPRI,*)'QQQ FRAGMNET SETTING'
!            call PRINT_LSSETTING(fragment(iorb)%setting,LUPRI)
!            
!            CALL lsQUIT('TEST TEST IORB')
!         ENDIF
!         call mat_free(D2)
!         call mat_free(tempm4)

         call ls_free(fragment_ls(iorb))
         !** F = h + G
         call mat_daxpy(1d0,tempm1,tempm3)
      ENDDO

      call mat_init(tempm4,nbast,nbast)
      call mat_add(1.d0,tempm3,-1.d0,refF,tempm4)
      write(lu_pri,*) 'QQQ DI_DEBUG     STD   ',mat_trab(refF,refF)
      write(lu_pri,*) 'QQQ DI_DEBUG lsfragment',mat_trab(tempm3,tempm3)
      write(lu_pri,*) 'QQQ DIFF               ',mat_trab(tempm4,tempm4)
      IF(ABS(mat_trab(tempm4,tempm4)).LE.1.D-8)THEN
         write(lu_pri,*)'QQQ SUCCESFUL'
      ELSE
         WRITE(lu_pri,*)'THE DIFF'
         call mat_print(tempm4,1,nbast,1,nbast,lu_pri)
         CALL lsQUIT('TEST TEST',lu_pri)
      ENDIF
      call mem_dealloc(orb2atom)
      call mat_free(tempm1)
      call mat_free(tempm2)
      call mat_free(tempm3)
      call mat_free(tempm4)
      call mat_free(refF)
      call free_aoitem(lu_pri,AO)

    END SUBROUTINE di_debug_ccfragment

   SUBROUTINE di_get_overlap_and_H1(S, H1,LSint,newlupri,newluerr,ls)
   use TYPEDEF
      IMPLICIT NONE
      TYPE(Matrix), target :: S, H1
      TYPE(Matrix) :: tempm1,tempm2,Kinetic
      integer  :: NNBAST, n2basx,newlupri,newluerr
      type(lsitem) :: ls
      logical :: LFOUND,LSint
      real(realk), ALLOCATABLE :: wrk(:), wrk1(:)
      !local pointer to H1
      LFOUND = .FALSE.
      lH1 => H1
      call II_get_overlap(newlupri,newluerr,ls%setting,S)
      call II_get_h1(newlupri,newluerr,ls%setting,H1)

   END SUBROUTINE di_get_overlap_and_H1

   SUBROUTINE di_get_fock_LSDALTON(D,h1,F,Etotal,newlupri,newluerr,ls)
   ! ===================================================================
   ! di_get_fock obtains total fock matrix and corresponding energy.
   ! WE have to go through the interface to dalton before the fock
   ! evaluator learns how to handle arbitrary-type arrays.
   ! ===================================================================
      use ks_settings
      use TYPEDEF
      IMPLICIT NONE
      TYPE(Matrix), INTENT(IN)    :: H1, D
     !TYPE(Matrix), save          :: D0, F0
      TYPE(Matrix), INTENT(INOUT) :: F
     !TYPE(Matrix)                :: Ddiff
      type(lsitem) :: ls
      real(realk), INTENT(OUT) :: Etotal
      real(realk)   :: edfty
      integer nbast, newlupri,newluerr
      integer :: error
      logical :: Dsym

      if (incremental_scheme) then
       call mat_assign(Ddiff,D)
      else
       call mat_clone(Ddiff, D)
      endif

      if (do_increment) call mat_daxpy(-1d0,D0,Ddiff)

      IF(.not.ls%input%dalton%unres) THEN !closed shell case

         Dsym = .TRUE. !symmetric Density matrix
         !IF(cfg_increase_integral_accuracy)THEN
         !   ls%setting%scheme%threshold = ls%setting%scheme%threshold*1.0D-1 
         !ENDIF
         ls%input%nfock = ls%input%nfock + 1
         call II_get_Fock_mat(newlupri,newluerr,ls%setting,Ddiff,Dsym,F,1)

         if (do_increment) call mat_daxpy(1d0,F0,F)
         !write(lupri,*) 'QQQ New  F:',mat_trab(F,F)

         if (incremental_scheme) then
          call mat_assign(F0,F); call mat_assign(D0,D)
         endif

!        Special case for incremental scheme (Fock-matrix build using density-difference
!        rather than the full density-matrix)

         Etotal = fockenergy_f(F,D,H1,ls%input%dalton%unres,ls%input%potnuc)
         IF(ls%setting%do_dft) THEN
            nbast = D%nrow
            !call printmem(lupri,'Memory before II_get_xc_fock_mat: ')
            call II_get_xc_fock_mat(newlupri,newluerr,ls%setting,nbast,D,F,Edfty)
            !write(lupri,*) 'QQQ New  Kohn-Sham:',mat_trab(F,F)
            !call printmem(lupri,'Memory after II_get_xc_fock_mat: ')
! I HAD TO REMOVE THE INCREMENTAL SCHEME AS IT CONFLICTS WITH OPENMP
!            call II_get_xc_fock_mat(newlupri,newluerr,ls%input,nbast,Ddiff,F,Edfty)
            Etotal = Etotal + Edfty
         ENDIF
      ELSE               !unrestricted open shell case
         Dsym = .TRUE. !symmetric Density matrix
         !IF(cfg_increase_integral_accuracy)THEN
         !   ls%setting%scheme%threshold = ls%setting%scheme%threshold*1.0D-1 
         !ENDIF
         ls%input%nfock = ls%input%nfock + 1
         call II_get_Fock_mat(newlupri,newluerr,ls%setting,Ddiff,Dsym,F,1)

         if (do_increment) call mat_daxpy(1d0,F0,F)

         if (incremental_scheme) then
          call mat_assign(F0,F); call mat_assign(D0,D)
         endif

!        Special case for incremental scheme (Fock-matrix build using density-difference
!        rather than the full density-matrix)

         Etotal = fockenergy_f(F,D,H1,ls%input%dalton%unres,ls%input%potnuc)

         IF(ls%setting%do_dft) THEN
            nbast = D%nrow
            call II_get_xc_fock_mat(newlupri,newluerr,ls%setting,nbast,D,F,Edfty)
            Etotal = Etotal + Edfty
         ENDIF
      ENDIF

      if (incremental_scheme) do_increment = .true. !call mat_free(Ddiff)

      !** F = h + G
      call mat_daxpy(1d0,H1,F)
      !call printmem(lupri,'Memory, end di_get_fock_LSDALTON: ')
   END SUBROUTINE di_get_fock_LSDALTON


   double precision function fockenergy_F(F,D,H1,unres,pot_nuc)
     !E = Tr(h + F)D + hnuc ! No factor since molecular orbitals
     implicit none
     TYPE(matrix), intent(in) :: F,D,H1
     double precision :: hD, FD, fac,pot_nuc
     double precision, external :: DF_E_ROB_CORR
     integer :: ndim
     logical :: unres
     
     ndim=F%nrow
     fac=2d0
     if(unres) fac=1d0
     !Get the one-electron hamiltonian
     !Tr(FD)
     !Tr(hD)
     hD =       mat_dotproduct(D,H1)
     FD = 0.5d0*mat_dotproduct(D,F)
     !E(HF) = Tr(h +F)D + hnuc
     fockenergy_F = (hd + FD)*fac + POT_NUC !Stinne: Remove call to boxed DF contribution (obsolete), 08-06-2010
     ! FIXME + potnuc
   end function fockenergy_F


   
!##########################################################
      subroutine di_GET_GbDsSingle(Dens,GbDs,nocca,noccb)
        use lsdalton_fock_module
        !*********************************************************
        ! Determine the G matrix for the 2-e contribution to sigma
        ! vector in RSP
        ! G([b,D]s) = 2-e part of Fock Matrix with a modified
        !             density [b,D]s (here called Dens)
        ! Sonia, October 2004
        ! Thomas, Feb 2010 (fixed unrestricted + added lsdalton lsint)
        !*********************************************************
        implicit none
        !Input/Output
        type(Matrix), intent(in) :: Dens
        type(Matrix), intent(inout) :: GbDs  !output
        !> Number of alpha spin electrons
        integer, intent(in) :: nocca
        !> Number of beta spin electrons
        integer, intent(in) :: noccb
        type(Matrix) :: TEMP,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9  !output
        !local
        real(realk), allocatable :: Dfull(:),Gfull(:)
        real(realk), allocatable :: DfullAB(:,:),FfullAB(:,:),FCAO(:),DCAO(:)
        integer :: nbast, nbast2,nfmat
        real(realk) :: dummy1, dummy2,fjab
        integer :: IFCTYP(100), IX, MATSYM, ISYMDM(100), ndmat,i
        logical :: direct
        integer :: IX1, IX2,IXAB,IXA,IXB
        real(realk) :: thrzer = 1.0d-14, norm_alfa, norm_beta
        logical, external  :: ISDIRECT
        logical :: alfa,Dsym,LINK
        Link = .FALSE.
        if (nocca+noccb /= 1) then !if only 1 electron -> 2 el part = 0
           Dsym = .FALSE. !NONsymmetric Density matrix
           !This should be changed to a test like the MATSYM function for full matrices
           IF(lsint_fock_data%ls%setting%scheme%link)THEN
              lsint_fock_data%ls%setting%scheme%link = .FALSE.
              LINK = .TRUE.
           ENDIF
           !            lsint_fock_data%ls%input%nfock = lsint_fock_data%ls%input%nfock + 1
           call II_get_Fock_mat(lsint_fock_data%lupri,lsint_fock_data%luerr,lsint_fock_data%ls%setting,Dens,Dsym,GbDs,1)
           IF(LINK)THEN
              lsint_fock_data%ls%setting%scheme%link = .TRUE.
           ENDIF
        endif
      end subroutine di_GET_GbDsSingle

      subroutine di_GET_GbDsArray(Dens,GbDs,nDmat,nocca,noccb)
        use lsdalton_fock_module
        !*********************************************************
        ! Determine the G matrix for the 2-e contribution to sigma
        ! vector in RSP
        ! G([b,D]s) = 2-e part of Fock Matrix with a modified
        !             density [b,D]s (here called Dens)
        ! Sonia, October 2004
        ! Thomas, Feb 2010 (fixed unrestricted + added lsdalton lsint)
        !*********************************************************
        implicit none
        !Input/Output
        type(Matrix), intent(in) :: Dens(nDmat)
        type(Matrix), intent(inout) :: GbDs(nDmat)  !output
        !> Number of alpha spin electrons
        integer, intent(in) :: nocca
        !> Number of beta spin electrons
        integer, intent(in) :: noccb
        type(Matrix),pointer :: temp(:),temp3(:)
        type(Matrix) :: temp2
        !local
        real(realk), allocatable :: Dfull(:),Gfull(:)
        real(realk), allocatable :: DfullAB(:,:),FfullAB(:,:),FCAO(:),DCAO(:)
        integer :: nbast, nbast2,nfmat
        real(realk) :: dummy1, dummy2,fjab
        integer :: IFCTYP(100), IX, MATSYM, ISYMDM(100), ndmat,i
        logical :: direct
        integer :: IX1, IX2,IXAB,IXA,IXB
        real(realk) :: thrzer = 1.0d-14, norm_alfa, norm_beta
        logical, external  :: ISDIRECT
        logical :: alfa,Dsym,LINK
        interface 
           SUBROUTINE II_get_Fock_mat(LUPRI,LUERR,SETTING,D,Dsym,F,ndmat)
             use precision
             use TYPEDEF  
             use Matrix_module
             use Matrix_Operations
             IMPLICIT NONE
             integer               :: ndmat
             TYPE(MATRIX)          :: D(ndmat),F(ndmat)
             TYPE(LSSETTING)       :: SETTING
             INTEGER               :: LUPRI,LUERR
             LOGICAL               :: Dsym
           end SUBROUTINE II_get_Fock_mat
        end interface
        Link = .FALSE.
         if (nocca+noccb /= 1) then !if only 1 electron -> 2 el part = 0
          Dsym = .FALSE. !NONsymmetric Density matrix
          !This should be changed to a test like the MATSYM function for full matrices
          IF(lsint_fock_data%ls%setting%scheme%link)THEN
           lsint_fock_data%ls%setting%scheme%link = .FALSE.
           LINK = .TRUE.
          ENDIF
!          lsint_fock_data%ls%input%nfock = lsint_fock_data%ls%input%nfock + 1
          call II_get_Fock_mat(lsint_fock_data%lupri,lsint_fock_data%luerr,&
               &lsint_fock_data%ls%setting,Dens,Dsym,GbDs,ndmat)
          IF(LINK)THEN
           lsint_fock_data%ls%setting%scheme%link = .TRUE.
          ENDIF
         endif
      end subroutine di_GET_GbDsArray
      
      subroutine di_GET_GbDsSingle_LSDALTON(lu_pri,lu_err,setting,Dens,GbDs)
        use typedef
        !*********************************************************
        ! Determine the G matrix for the 2-e contribution to sigma
        ! vector in RSP
        ! G([b,D]s) = 2-e part of Fock Matrix with a modified
        !             density [b,D]s (here called Dens)
        ! Thomas, sep 2010 
        !*********************************************************
        implicit none
        type(lssetting),intent(inout) :: setting
        integer, intent(in) :: lu_pri,lu_err
        type(Matrix), intent(in) :: Dens
        type(Matrix), intent(inout) :: GbDs  
 !       
        logical :: Dsym,LINK
        Link = .FALSE.
!        if (cfg_nocca+cfg_noccb /= 1) then !if only 1 electron -> 2 el part = 0
           Dsym = .FALSE. !NONsymmetric Density matrix
           !This should be changed to a test like the MATSYM function for full matrices
           IF(setting%scheme%link)THEN
              setting%scheme%link = .FALSE.
              LINK = .TRUE.
           ENDIF
           call II_get_Fock_mat(lu_pri,lu_err,setting,Dens,Dsym,GbDs,1)
           IF(LINK)THEN
              setting%scheme%link = .TRUE.
            ENDIF
!        endif
      end subroutine di_GET_GbDsSingle_LSDALTON

      subroutine di_GET_GbDsArray_LSDALTON(lu_pri,lu_err,setting,Dens,GbDs,nDmat)
        use lsdalton_fock_module
        !*********************************************************
        ! Determine the G matrix for the 2-e contribution to sigma
        ! vector in RSP
        ! G([b,D]s) = 2-e part of Fock Matrix with a modified
        !             density [b,D]s (here called Dens)
        ! Sonia, October 2004
        ! Thomas, Feb 2010 (fixed unrestricted + added lsdalton lsint)
        !*********************************************************
        implicit none
        !Input/Output
        type(lssetting),intent(inout) :: setting
        integer, intent(in) :: lu_pri,lu_err,ndmat
        type(Matrix), intent(in) :: Dens(nDmat)
        type(Matrix), intent(inout) :: GbDs(nDmat)  !output
        !
        integer :: i
        logical :: Dsym, LINK
        interface 
           SUBROUTINE II_get_Fock_mat(LUPRI,LUERR,SETTING,D,Dsym,F,ndmat)
             use precision
             use TYPEDEF  
             use Matrix_module
             use Matrix_Operations
             IMPLICIT NONE
             integer               :: ndmat
             TYPE(MATRIX)          :: D(ndmat),F(ndmat)
             TYPE(LSSETTING)       :: SETTING
             INTEGER               :: LUPRI,LUERR
             LOGICAL               :: Dsym
           end SUBROUTINE II_get_Fock_mat
        end interface

        Link = .FALSE.
        !I do not have acces to these things right now
!        if (cfg_nocca+cfg_noccb /= 1) then !if only 1 electron -> 2 el part = 0
           IF(setting%scheme%link)THEN
              setting%scheme%link = .FALSE.
              LINK = .TRUE.
           ENDIF
           !This should be changed to a test like the MATSYM function for full matrices   
           Dsym = .FALSE. !NONsymmetric Density matrix
           call II_get_Fock_mat(lu_pri,lu_err,setting,Dens,Dsym,GbDs,ndmat)
           IF(LINK)THEN
              setting%scheme%link = .TRUE.
           ENDIF
!        endif
      end subroutine di_GET_GbDsArray_LSDALTON

      subroutine di_lowdin(S,S_sqrt,S_minus_sqrt,n,lupri)
        use matrix_operations
        use lowdin_module
        implicit none
        real(realk), target :: S(n*n), S_sqrt(n*n), S_minus_sqrt(n*n)
        type(Matrix)        :: tS, tS_sqrt, tS_minus_sqrt
        integer             :: n, lupri
        
        call lowdin_diag(n,S,S_sqrt, S_minus_sqrt, lupri)
#if 0 && defined(HAVE_BSM)
        select case (matrix_type)
        case(mtype_sparse_block)
           
           call mat_init(tS,n,n)
           call mat_init(tS_sqrt,n,n)
           call mat_init(tS_minus_sqrt,n,n)
           
           call mat_set_from_full(S,1.0d0,tS)
           
           call lowdin_schulz(tS,tS_sqrt,tS_minus_sqrt,lupri)
           
           call mat_to_full(tS_sqrt,1.0d0,S_sqrt)
           call mat_to_full(tS_minus_sqrt,1.0d0,S_minus_sqrt)
           
           call mat_free(tS)
           call mat_free(tS_sqrt)
           call mat_free(tS_minus_sqrt)
        case default
           call lsquit('Only dense and BSM matrices are supported by di_lowdin()!',-1)
        end select
#endif
      end subroutine di_lowdin

    end MODULE dal_interface
