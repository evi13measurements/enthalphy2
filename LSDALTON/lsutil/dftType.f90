!> dft type module
!> \author T.Kjaergaard
!> \date 2010-02-21
MODULE dft_type
use precision

!> USED IN II_DFTINT TO save data like Fockmatrix and gradients
!> calculated at each gridpoint 
TYPE DFTDATATYPE
INTEGER             :: nbast
INTEGER             :: ndmat !# Density matrices
INTEGER             :: nbmat !# response vectors
real(realk)         :: energy
real(realk)         :: electrons
real(realk),pointer :: BMAT(:,:,:)!nbast,nbast,nbmat
real(realk),pointer :: FKSM(:,:,:)!nbast,nbast,ndmat
real(realk),pointer :: FKSMS(:,:,:)!nbast,nbast,ndmat
logical             :: dosympart
INTEGER             :: natoms
INTEGER,pointer     :: orb2atom(:)
real(realk),pointer :: grad(:,:)
END TYPE DFTDATATYPE

END MODULE dft_type
