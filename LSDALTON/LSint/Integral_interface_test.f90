!> @file Interface file used to test implemented functionality (not covered by other tests) - 
!>       for instance integrals for CC dec machinery should be placed in this file
SUBROUTINE II_test_tco_densityContraction(LUPRIOLD,LUPRI,LUERR,SETTING,D,auxoption)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use tco_Integral_Interface
  use linsolvdf
IMPLICIT NONE
TYPE(MATRIX),target   :: D
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRIOLD,LUPRI,LUERR
LOGICAL               :: auxoption
!
Integer             :: ndmat,nbasis,naux
Real(realk),pointer :: Dfull(:,:,:)
Real(realk),pointer :: tco_Ffull(:,:,:)
Real(realk),pointer :: galphafull(:,:,:),tco_galphafull(:,:,:)
Real(realk),pointer :: calphafull(:,:,:)
type(matrix),target :: galpha,calpha,F
type(matrix)        :: tco_galpha,tco_F,TMP
REAL(REALK)         :: SUM,tco_SUM,diff
CHARACTER*80        :: BasisType
INTEGER             :: i,j,iprint
type(matrixp)       :: Jmat(1),Dmat(1)
Logical             :: FRAGMENT_SAVE

FRAGMENT_SAVE = SETTING%SCHEME%FRAGMENT
SETTING%SCHEME%FRAGMENT = .FALSE.
!CALL LSHEADER(lupri,'II_test_tco_densityContraction')
ndmat = 1
iprint = SETTING%SCHEME%INTPRINT
iprint=10000
nbasis = SETTING%MOLECULE(1)%p%nbastREG
IF(auxoption) THEN
  BasisType = 'DF-Aux'
  naux = SETTING%MOLECULE(1)%p%nbastAUX
ELSE
  BasisType = 'Regular'
  naux = SETTING%MOLECULE(1)%p%nbastREG
ENDIF

IF(matrix_type .EQ. mtype_unres_dense)CALL LSQUIT('ndmat hardcoded to 1 in II_test_tco_densityContraction',lupri)
call mem_alloc(Dfull,nbasis,nbasis,1)
call mat_to_full(D,1D0,Dfull(:,:,1))
call mat_init(galpha,naux,1)
call mat_zero(galpha)
Jmat(1)%p => galpha
Dmat(1)%p => D
CALL ls_attachDmatToSetting(Dmat,1,setting,'RHS',3,4,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,galpha%nrow,galpha%ncol,1,1,ndmat)
print*,'LS_JENGINE IN II_test_tco_densityContraction'
call ls_jengine(BasisType,'Empty','Regular','Regular','Overlap','Regular','Contracted',&
                SETTING,LUPRI,LUERR)
print*,'DONE LS_JENGINE IN II_test_tco_densityContraction'
CALL retrieve_Output(lupri,setting,galpha)

CALL ls_freeDmatFromSetting(setting)

call mem_alloc(tco_galphafull,naux,1,1)
CALL ls_DZERO(tco_galphafull,naux)
call tco_contract_matrix(tco_galphafull,BasisType,'Regular','Regular','Contracted',&
                         Dfull,nbasis,nbasis,SETTING,LUPRI,LUERR)
call mem_dealloc(Dfull)
!in ls_jengine for (AO3='Regular' AND AO4='Regular') integrals are multiplied by CoulombFactor=2

call mat_init(tco_galpha,naux,1)
call mat_set_from_full(tco_galphafull(:,:,1),2.d0,tco_galpha)
call mem_dealloc(tco_galphafull)
IF(IPRINT.GT.100)THEN
   WRITE(LUPRI,'(2X,A)')'3 Center Overlap contracted with matrix (a|b|c) x m_bc'
   call mat_print(galpha,1,naux,1,1,lupri)
   call mat_print(tco_galpha,1,naux,1,1,lupri)
   call mat_init(TMP,naux,1)
   call mat_add(1.D0,tco_galpha,-1.D0,galpha,TMP)
   call mat_print(TMP,1,naux,1,1,lupri)
   call mat_free(TMP)
ENDIF

call mat_init(TMP,naux,1)
call mat_add(1.D0,tco_galpha,-1.D0,galpha,TMP)
call mat_print(TMP,1,naux,1,1,lupri)
WRITE(LUPRIOLD,'(2X,A,F20.12)')'SUM OF ALL 3 CENTER OVERLAPS CONTRACTED WITH MATRIX'
WRITE(LUPRIOLD,'(2X,A,F20.12)')'FROM THE STANDARD DRIVER ',mat_sum(galpha)
WRITE(LUPRIOLD,'(2X,A,F20.12)')'FROM THE TCO DRIVER      ',mat_sum(tco_galpha)
WRITE(LUPRIOLD,'(2X,A,D20.9)') 'ABSOLUTE DIFFERENCE      ',mat_sum(TMP)
call mat_free(TMP)
call mat_free(tco_galpha)

call mem_alloc(galphafull,naux,1,1)
call mat_to_full(galpha,1D0,galphafull(:,:,1))
call mat_free(galpha)
call mem_alloc(calphafull,naux,1,1)
call linsolv_df(calphafull,galphafull,BasisType,'Overlap',naux,ndmat,SETTING,LUPRI,LUERR)
call mem_dealloc(galphafull)
call mat_init(calpha,naux,1)
call mat_set_from_full(calphafull(:,:,1),1.d0,calpha)
call mat_init(F,nbasis,nbasis)
CALL MAT_ZERO(F)
Jmat(1)%p => F
Dmat(1)%p => calpha
CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,F%nrow,F%ncol,1,1,1)
call ls_jengine('Regular','Regular',BasisType,'Empty','Overlap','Regular','Contracted',&
                SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F)

CALL ls_freeDmatFromSetting(setting)

call mat_free(calpha)
call mem_alloc(tco_Ffull,nbasis,nbasis,1)
CALL ls_DZERO(tco_Ffull,nbasis*nbasis)
call tco_contract_vector(tco_Ffull,'Regular','Regular',BasisType,'Contracted',&
                         calphafull,naux,SETTING,LUPRI,LUERR)
call mem_dealloc(calphafull)
call mat_init(tco_F,nbasis,nbasis)
call mat_set_from_full(tco_Ffull(:,:,1),1.d0,tco_F)

IF(IPRINT.GT.100)THEN
   WRITE(LUPRI,'(2X,A)')'3 Center Overlap contracted with vector (a|b|c) x v_c'
   call mat_init(TMP,nbasis,nbasis)
   call mat_print(F,1,nbasis,1,nbasis,lupri)
   call mat_print(tco_F,1,nbasis,1,nbasis,lupri)
   call mat_add(1.D0,tco_F,-1.D0,F,TMP)
   call mat_print(TMP,1,nbasis,1,nbasis,lupri)
   call mat_free(TMP)
ENDIF
call mat_init(TMP,nbasis,nbasis)
call mat_add(1.D0,tco_F,-1.D0,F,TMP)
WRITE(LUPRIOLD,'(2X,A,F20.12)')'SUM OF ALL 3 CENTER OVERLAPS CONTRACTED WITH VECTOR'
WRITE(LUPRIOLD,'(2X,A,F20.12)')'FROM THE STANDARD DRIVER ',mat_sum(F)
WRITE(LUPRIOLD,'(2X,A,F20.12)')'FROM THE TCO DRIVER      ',mat_sum(tco_F)
WRITE(LUPRIOLD,'(2X,A,D20.9)') 'ABSOLUTE DIFFERENCE      ',mat_sum(TMP)
call mat_free(TMP)

call mem_dealloc(tco_Ffull)
call mat_free(F)
call mat_free(tco_F)

SETTING%SCHEME%FRAGMENT = FRAGMENT_SAVE

END SUBROUTINE II_test_tco_densityContraction

SUBROUTINE II_test_coulomb_ndmat(LUPRI,LUERR,SETTING,D,F)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(MATRIX),target   :: D,F
TYPE(MATRIX),target   :: D2(2),F2(2),F3
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: ndmat
Real(realk),pointer :: Dfull(:,:,:),TMP(:,:)
Real(realk),pointer :: Ffull(:,:,:,:,:)
Real(realk)         :: TS,TE
type(matrixp)       :: Jmat(2),Dmat(2)

!THIS IS ONLY A TEST ROUTINE TO MAKE SURE THAT COULOMB
! WORKS FOR MORE THAN 1 MATRIX, AT THE MOMENT (19/11/2009) WE HAVE NO NEED FOR
! THE COULOMB TO WORK WITH MORE THAN 1 MATRIX AT A TIME BUT THIS 
! SUBROUTINE IS TO MAINTAIN THAT FEATURE - FOR LATER USEGE
ndmat = 2
call mat_init(D2(1),D%nrow,D%ncol)
call mat_init(D2(2),D%nrow,D%ncol)
D2(1) = D
D2(2) = D
call mat_scal(5.123456789123d0,D2(1))
call mat_init(F2(1),D%nrow,D%ncol)
call mat_init(F2(2),D%nrow,D%ncol)
CALL mat_ZERO(F)
CALL mat_ZERO(F2(1))
CALL mat_ZERO(F2(2))

IF (SETTING%SCHEME%JENGINE) THEN     
   Dmat(1)%p => D2(1)
   Dmat(2)%p => D2(2)
   CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,lupri)
   call initIntegralOutputDims(setting%Output,F2(1)%nrow,F2(1)%ncol,1,1,2)
   call ls_jengine('Regular','Regular','Regular','Regular','Coulomb','Regular','Contracted',&
        &          SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,F2)
   CALL ls_freeDmatFromSetting(setting)
ELSE
   call mem_alloc(Dfull,D%nrow,D%ncol,ndmat)
   IF(matrix_type .EQ. mtype_unres_dense)THEN
      CALL DCOPY(D2(1)%nrow*D2(1)%ncol,D2(1)%elms,1,Dfull(:,:,1),1)
      CALL DAXPY(D2(1)%nrow*D2(1)%ncol,1.d0,D2(1)%elmsb,1,Dfull(:,:,1),1)
      CALL DCOPY(D2(2)%nrow*D2(2)%ncol,D2(2)%elms,1,Dfull(:,:,2),1)
      CALL DAXPY(D2(2)%nrow*D2(2)%ncol,1.d0,D2(2)%elmsb,1,Dfull(:,:,2),1)
   ELSE !CLOSED_SHELL
      call mat_to_full(D2(1),1D0,Dfull(:,:,1))
      call mat_to_full(D2(2),1D0,Dfull(:,:,2))
   ENDIF
   call mem_alloc(Ffull,F%nrow,F%ncol,1,1,ndmat)
   call initIntegralOutputDims(setting%Output,F%nrow,F%ncol,1,1,ndmat)
   CALL ls_DZERO(Ffull,F%nrow*F%ncol*ndmat)
   CALL ls_attachDmatToSetting(Dfull,ndmat,setting,'RHS',3,4,lupri)
   call ls_get_coulomb_mat('Regular','Regular','Regular','Regular',&
        &                   'Coulomb','Contracted',SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,F2)
   call mem_dealloc(Dfull)
   call mem_dealloc(Ffull)
ENDIF

call mat_init(F3,F%nrow,F%ncol)
CALL mat_ZERO(F3)
call II_get_coulomb_mat(LUPRI,LUERR,SETTING,D2(1),F3,1)
call mat_add(1.D0,F3,-1.D0,F2(1),D2(1))
IF(ABS(MAT_SUM(D2(1))) .GT. 1D-10)THEN
   write(LUPRI,*)'QQQ REF  ',mat_trab(F3,F3)
   write(LUPRI,*)'QQQ TEST ',mat_trab(F2(1),F2(1))
   write(LUPRI,*)'TMP PRINT ',mat_trab(F2(2),F2(2))
   WRITE(LUPRI,*)'UNSUCCESFULL NDMAT TEST'
ELSE
   write(LUPRI,*)'QQQ REF  ',mat_trab(F3,F3)
   write(LUPRI,*)'QQQ TEST ',mat_trab(F2(1),F2(1))
   write(LUPRI,*)'TMP PRINT ',mat_trab(F2(2),F2(2))
   WRITE(LUPRI,*)'SUCCESFULL NDMAT TEST'
   F = F2(2) 
ENDIF
call mat_free(F3)
call mat_free(F2(1))
call mat_free(F2(2))
call mat_free(D2(1))
call mat_free(D2(2))

END SUBROUTINE II_test_coulomb_ndmat

SUBROUTINE II_test_exchange_ndmat(LUPRI,LUERR,SETTING,D,F)
  use precision
  use TYPEDEF  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
IMPLICIT NONE
TYPE(MATRIX),target   :: D,F
TYPE(MATRIX),target   :: D2(2),F2(2),F3
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: ndmat
Real(realk),pointer :: Dfull(:,:,:),TMP(:,:)
Real(realk),pointer :: Ffull(:,:,:,:,:)
Real(realk)         :: TS,TE
type(matrixp)       :: Jmat(2),Dmat(2)

!THIS IS ONLY A TEST ROUTINE TO MAKE SURE THAT COULOMB
! WORKS FOR MORE THAN 1 MATRIX, AT THE MOMENT (19/11/2009) WE HAVE NO NEED FOR
! THE COULOMB TO WORK WITH MORE THAN 1 MATRIX AT A TIME BUT THIS 
! SUBROUTINE IS TO MAINTAIN THAT FEATURE - FOR LATER USEGE
ndmat = 2
call mat_init(D2(1),D%nrow,D%ncol)
call mat_init(D2(2),D%nrow,D%ncol)
D2(1) = D
D2(2) = D
call mat_scal(5.123456789123d0,D2(1))
call mat_init(F2(1),D%nrow,D%ncol)
call mat_init(F2(2),D%nrow,D%ncol)
CALL mat_ZERO(F2(1))
CALL mat_ZERO(F2(2))

call mem_alloc(Ffull,F%nrow,F%ncol,1,1,ndmat)
CALL ls_DZERO(Ffull,F%nrow*F%ncol*ndmat)
call II_get_exchange_mat(LUPRI,LUERR,SETTING,D2,2,.true.,F2)
call mem_dealloc(Ffull)

call mat_init(F3,F%nrow,F%ncol)
CALL mat_ZERO(F3)
call II_get_exchange_mat(LUPRI,LUERR,SETTING,D2(1),1,.true.,F3)
call mat_add(1.D0,F3,-1.D0,F2(1),D2(1))
IF(ABS(MAT_SUM(D2(1))) .GT. 1D-10)THEN
   write(LUPRI,*)'QQQ REF  ',mat_trab(F3,F3)
   write(LUPRI,*)'QQQ TEST ',mat_trab(F2(1),F2(1))
   write(LUPRI,*)'TMP PRINT ',mat_trab(F2(2),F2(2))
   WRITE(LUPRI,*)'UNSUCCESFULL NDMAT TEST'
ELSE
   write(LUPRI,*)'QQQ REF  ',mat_trab(F3,F3)
   write(LUPRI,*)'QQQ TEST ',mat_trab(F2(1),F2(1))
   write(LUPRI,*)'TMP PRINT ',mat_trab(F2(2),F2(2))
   WRITE(LUPRI,*)'SUCCESFULL NDMAT TEST'
   call mat_daxpy(1.d0,F2(2),F)
ENDIF
call mat_free(F3)
call mat_free(F2(1))
call mat_free(F2(2))
call mat_free(D2(1))
call mat_free(D2(2))

END SUBROUTINE II_test_exchange_ndmat


