module mlcc_drive
!
!  Contains input routine called from mlcc_interface_input outside
!  and stores input in mlcc_data
!
contains
!
subroutine mlcc_drv(work,lwork,lupri)
!
   use mlcc_types
   use mlcc_workspace
   use mlcc_data
   use mlcc_init
   use mlcc_utilities
   use mlcc_input_data
   use mlcc_t2_init
!
!  mlcc3 driver
!  Author Rolf H. Myhre
!  December 2016
!
   implicit none
!
   integer, intent(in)                    :: lupri !general output unit
   integer, intent(in)                    :: lwork !free space in work
!
   real(dp), intent(in), dimension(lwork) :: work !work static array
!
   integer :: i,j,a,b,nai,nbj,naibj
!
   ml_lupri = lupri 
!
   write(ml_lupri,*)
   write(ml_lupri,*) 'In mlcc_drv'
   write(ml_lupri,*)
!
   write(ml_lupri,*)
   write(ml_lupri,*) 'mlcc_active: ', mlcc_active
   write(ml_lupri,*) 'print_mlcc:  ', print_mlcc 
   write(ml_lupri,*)
!
   call work_init(mem)
   call hf_reader()
!
!  Allocate amplitudes
!
   call allocator(t2am,n_t2am_pack,1)
   t2am=zero
!  Read in IAJB integrals
!
   call mlcc_iajb(t2am)
!
!
   call mlcc_get_cholesky()
   call t2_init(t2am)
   call deallocator(t2am,n_t2am_pack,1)
end subroutine mlcc_drv
!
end module mlcc_drive
