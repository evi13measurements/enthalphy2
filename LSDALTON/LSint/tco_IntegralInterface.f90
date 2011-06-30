!> @file
!> Contains module for bridging the integral interface to the toc main driver.
MODULE tco_Integral_Interface
use precision
use typedef
use buildbasisset
use buildaobatch
use tco_Integral_driver
use molecule_module

CONTAINS

SUBROUTINE tco_getIntegrals(integrals,AO1,AO2,AO3,intType,SETTING,LUPRI,LUERR)
implicit none
REAL(REALK),pointer      :: integrals(:,:,:)
CHARACTER*(*)            :: AO1,AO2,AO3,intType
TYPE(LSSETTING)          :: SETTING
INTEGER                  :: LUPRI,LUERR
!
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
TYPE(AOITEM),target      :: AObuild(3)
INTEGER                  :: nAObuilds

CALL init_TCO_integralinput(INPUT,SETTING)
CALL set_TCO_inputAO(INPUT,AO1,AO2,AO3,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR)
OUTPUT%ndim(1) = INPUT%AOdim(1)
OUTPUT%ndim(2) = INPUT%AOdim(2)
OUTPUT%ndim(3) = INPUT%AOdim(3)
OUTPUT%ResultMat => integrals

CALL TCO_INTEGRALDRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INPUT,OUTPUT)

CALL free_TCO_inputAO(AObuild,nAObuilds,LUPRI)

END SUBROUTINE tco_getIntegrals

SUBROUTINE tco_contract_vector(integrals,AO1,AO2,AO3,intType,Vect,n,&
                               SETTING,LUPRI,LUERR)
implicit none
REAL(REALK),pointer      :: integrals(:,:,:)
CHARACTER*(*)            :: AO1,AO2,AO3,intType
REAL(REALK)              :: Vect(n)
INTEGER                  :: n
TYPE(LSSETTING)          :: SETTING
INTEGER                  :: LUPRI,LUERR

CALL tco_contract_driver(integrals,AO1,AO2,AO3,intType,1,Vect,n,1,&
                         SETTING,LUPRI,LUERR)

END SUBROUTINE tco_contract_vector

SUBROUTINE tco_contract_matrix(integrals,AO1,AO2,AO3,intType,Matr,n1,n2,&
                               SETTING,LUPRI,LUERR)
implicit none
REAL(REALK),pointer      :: integrals(:,:,:)
CHARACTER*(*)            :: AO1,AO2,AO3,intType
REAL(REALK)              :: Matr(n1,n2)
INTEGER                  :: n1,n2
TYPE(LSSETTING)          :: SETTING
INTEGER                  :: LUPRI,LUERR

CALL tco_contract_driver(integrals,AO1,AO2,AO3,intType,2,Matr,n1,n2,&
                         SETTING,LUPRI,LUERR)

END SUBROUTINE tco_contract_matrix

SUBROUTINE tco_contract_driver(integrals,AO1,AO2,AO3,intType,contractionType,Matr,n1,n2,&
                               SETTING,LUPRI,LUERR)
implicit none
REAL(REALK),pointer      :: integrals(:,:,:)
CHARACTER*(*)            :: AO1,AO2,AO3,intType
REAL(REALK)              :: Matr(n1,n2)
INTEGER                  :: contractionType,n1,n2
TYPE(LSSETTING)          :: SETTING
INTEGER                  :: LUPRI,LUERR
!
TYPE(TCO_INTEGRALINPUT)  :: INPUT
TYPE(TCO_INTEGRALOUTPUT) :: OUTPUT
TYPE(AOITEM),target      :: AObuild(3)
INTEGER                  :: nAObuilds,ntemp
LOGICAL                  :: correct_dimension

CALL init_TCO_integralinput(INPUT,SETTING)
IF(contractionType.EQ.1) THEN
  INPUT%CONT_VECTOR = .TRUE.
ELSEIF(contractionType.EQ.2) THEN
  INPUT%CONT_MATRIX = .TRUE.
ELSE
  CALL LSQUIT('Programming error, wrong contraction type in TCO_CONTRACT',lupri)
ENDIF

CALL set_TCO_inputAO(INPUT,AO1,AO2,AO3,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR)
IF(INPUT%CONT_VECTOR) THEN
  correct_dimension = n1.EQ.INPUT%AOdim(3)
ELSEIF(INPUT%CONT_MATRIX) THEN
  correct_dimension = n1.EQ.INPUT%AOdim(2) .AND. n2.EQ.INPUT%AOdim(3)
ELSE
  CALL LSQUIT('Programming error, wrong contraction type in TCO_CONTRACT',lupri)
ENDIF
IF(.NOT.correct_dimension) THEN
  CALL LSQUIT('ERROR: Wrong dimensions in TCO_CONTRACT',lupri)
ENDIF
CALL attach_TCO_matrixToInput(INPUT,Matr,n1,n2)

IF(INPUT%CONT_VECTOR) THEN
  OUTPUT%ndim(1) = INPUT%AOdim(1)
  OUTPUT%ndim(2) = INPUT%AOdim(2)
  OUTPUT%ndim(3) = 1
ELSEIF(INPUT%CONT_MATRIX) THEN
  OUTPUT%ndim(1) = INPUT%AOdim(1)
  OUTPUT%ndim(2) = 1
  OUTPUT%ndim(3) = 1
ELSE
  CALL LSQUIT('Programming error, wrong contraction type in TCO_CONTRACT',lupri)
ENDIF
OUTPUT%ResultMat => integrals

CALL TCO_INTEGRALDRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INPUT,OUTPUT)

CALL deattach_TCO_matrixFromInput(INPUT)
CALL free_TCO_inputAO(AObuild,nAObuilds,LUPRI)

END SUBROUTINE tco_contract_driver

SUBROUTINE set_TCO_inputAO(INPUT,AO1,AO2,AO3,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR)
implicit none
TYPE(TCO_INTEGRALINPUT)  :: INPUT
CHARACTER*(*)            :: AO1,AO2,AO3,intType
INTEGER                  :: nAObuilds
TYPE(AOITEM),target      :: AObuild(3)
TYPE(LSSETTING)          :: SETTING
INTEGER                  :: LUPRI,LUERR
!
CHARACTER(80)            :: AOstring(3)
TYPE(BASISSETINFO)       :: AObasis
INTEGER                  :: IAO,JAO,indAO,indexUnique(3),ndim(3)
LOGICAL                  :: uncont,intnrm,uniqueAO,emptyAO

IF(intType.EQ.'Primitive') THEN
  uncont = .TRUE.
  intnrm = .TRUE.
ELSEIF(intType.EQ.'Contracted') THEN
  uncont = .FALSE.
  intnrm = .FALSE.
ELSE
  WRITE(LUPRI,'(1X,2A)') 'Wrong case in Set_TCO_InputAO, intType =',intType
  CALL LSQUIT('Error - wrong intType in Set_TCO_InputAO',lupri)
ENDIF

INPUT%same12AOs = AO1.EQ.AO2 .AND. AO1.NE.'Empty'
INPUT%same13AOs = AO1.EQ.AO3 .AND. AO1.NE.'Empty'
INPUT%same23AOs = AO2.EQ.AO3 .AND. AO2.NE.'Empty'

AOstring(1) = AO1
AOstring(2) = AO2
AOstring(3) = AO3

nAObuilds = 0
DO IAO=1,3
  uniqueAO = .TRUE.
  DO JAO=1,IAO-1
    IF(AOstring(IAO).EQ.AOstring(JAO)) THEN
      uniqueAO = .FALSE.
      indAO=indexUnique(JAO)
      EXIT
    ENDIF
  ENDDO
  IF (uniqueAO) THEN
    nAObuilds = nAObuilds+1
    indAO = nAObuilds
    indexUnique(IAO) = indAO
    emptyAO = .FALSE.
    SELECT CASE(AOstring(IAO))
    CASE ('Regular')
      AObasis = SETTING%BASIS(1)%p%REGULAR
    CASE ('Huckel')
      AObasis = SETTING%BASIS(1)%p%HUCKEL
    CASE ('DF-Aux')
      AObasis = SETTING%BASIS(1)%p%AUXILIARY
    CASE ('Empty')
      emptyAO = .TRUE.
      CALL BUILD_EMPTY_AO(AObuild(indAO),LUPRI)
      ndim(indAO) = 1
    CASE DEFAULT
      WRITE(*,*)     'case: ',AOstring(IAO)
      WRITE(lupri,*) 'case: ',AOstring(IAO)
      WRITE(luerr,*) 'case: ',AOstring(IAO)
      CALL LSQuit('Programming error: Not a case in Set_TCO_InputAO!',lupri)
    END SELECT
    IF (.NOT.emptyAO) THEN
      CALL BUILD_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
                    SETTING%MOLECULE(1)%p,AObasis,&
                    AObuild(indAO),'Default',uncont,intnrm)
      ndim(indAO) = getNbasis(AOstring(IAO),intType,SETTING%MOLECULE(1)%p,LUPRI)
    ENDIF
  ENDIF
  INPUT%AO(IAO)%p => AObuild(indAO)
  INPUT%AOdim(IAO) = ndim(indAO)
ENDDO

END SUBROUTINE set_TCO_inputAO

SUBROUTINE free_TCO_inputAO(AObuild,nAObuilds,LUPRI)
implicit none
TYPE(AOITEM),target  :: AObuild(3)
INTEGER              :: nAObuilds,LUPRI
!
INTEGER              :: iAO

DO iAO=1,nAObuilds
  CALL free_aoitem(LUPRI,AObuild(IAO))
ENDDO
 
END SUBROUTINE free_TCO_inputAO

END MODULE tco_Integral_Interface
