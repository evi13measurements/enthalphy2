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
   use mlcc_omega
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
!  Initialize memory variables
!
   call work_init(mem)
!
!  Read in Fock diagonal and MO coefficients from Dalton's 
!  Hartree-Fock routine
!
   call hf_reader
!
!
!  Allocate amplitudes and the omega vector 
!
   call allocator(omega1,n_vir,n_occ) ! Omega_a,i
   call allocator(omega2,n_t2am_pack,1)
!
   call allocator(t1am,n_vir,n_occ)
   call allocator(t2am,n_t2am_pack,1)
!
   omega1 = zero
   omega2 = zero
!
   t1am = zero
   t2am = zero
!
!  Read in Cholesky vectors in AO basis, transform them to MO basis,
!  and save the MO Cholesky vectors to file (IJ,AB,IA,AI)
!
   call mlcc_get_cholesky
!
!  Set initial guess for the doubles amplitudes 
!
call allocator(mo_fock_mat,n_orbitals,n_orbitals)
call mlcc_fock
  call t2_init
!
!  Calculate the omega vector 
!
   call mlcc_omega_calc
!
   call deallocator(t2am,n_t2am_pack,1)
   call deallocator(fock_diagonal,n_orbitals,1)
   call deallocator(orb_coefficients,n_lambda,1)
   call deallocator(mo_fock_mat,n_orbitals,n_orbitals)
!
!
end subroutine mlcc_drv
!
end module mlcc_drive
