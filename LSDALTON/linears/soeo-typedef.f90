!> @file
!> Contains the soeo_typedef module
!
!> brief Structure with information for SOEO
!> author C. Nygaard
!> date 2010.07.07
module soeo_typedef

use matrix_util

implicit none

!Contains:
!  type soeoItem

type soeoItem

!New AO matrices:
!> Density matrix in AO basis
type(matrix)              :: Dao
!> Fock matrix in the AO basis
type(matrix)              :: Fao
!> Overlap matrix
type(matrix)              :: S
!> MO coefficient matrix
type(matrix)              :: C

!New MO matrices:
!> Density matrix in the MO basis
type(matrix)              :: Dmo
!> Fock matrix in the MO basis
type(matrix)              :: Fmo
!> First derivatives of Dmo wrt occupations
type(matrix), allocatable :: nfirst(:)
!> Second derivatives of Dmo wrt occupations
type(matrix), allocatable :: nsecond(:,:)
!> Occupations Dmo(i,i) = cos^2(theta(i))
type(matrix)              :: oldtheta

!Old matrices (matrices from last iteration is saved in case of rejection of the step)
!> Old occupations
type(matrix)              :: old_oldtheta
!> Old density in MO basis
type(matrix)              :: old_Dmo
!> Old density in AO basis
type(matrix)              :: old_Dao
!> Old Fock matrix in MO basis
type(matrix)              :: old_Fmo
!> Old Fock matrix in AO basis
type(matrix)              :: old_Fao
!> Old MO coefficient matrix
type(matrix)              :: old_C
!> Old first derivatives of Dmo wrt occupations
type(matrix), allocatable :: old_nfirst(:)
!> Old second derivatives of Dmo wrt occupations
type(matrix), allocatable :: old_nsecond(:,:)

!Space:
!> Size of all space (number o fbasis functions)
integer                   :: Nbast
!> Size of occupied space
integer                   :: Nocc
!> Size of active space
integer                   :: Nact

!Energy:
!> The total energy
real(realk)               :: Etotal
!> The predicted energy change
real(realk)               :: dEpred
!> The actual energy change
real(realk)               :: dE

!Thresholds and stuff:
!> Maximum macro-iterations
integer                   :: macromaxiter
!> Convergence threshold for macro-iterations
real(realk)               :: macrothresh
!> Maximum micro-iterations
integer                   :: micromaxiter
!> Convergence threshold for micro-iterations
real(realk)               :: microthresh
!> Trust radius
real(realk)               :: trust

end type soeoItem

end module soeo_typedef
