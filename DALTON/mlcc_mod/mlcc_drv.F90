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
   luprint = lupri 
!
   write(luprint,*)
   write(luprint,*) 'In mlcc_drv'
   write(luprint,*)
!
   write(luprint,*)
   write(luprint,*) 'mlcc_active: ', mlcc_active
   write(luprint,*) 'print_mlcc:  ', print_mlcc 
   write(luprint,*)
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
   call allocator(omega2,n_ov_ov_packed,1)
!
   call allocator(t1am,n_vir,n_occ)
   call allocator(t2am,n_ov_ov_packed,1)
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
!   call mlcc_fock
   call t2_init
!
!  Calculate the omega vector 
!
!   call mlcc_omega_calc
!
   call deallocator(t2am,n_ov_ov_packed,1)
   call deallocator(fock_diagonal,n_orbitals,1)
   call deallocator(orb_coefficients,n_lambda,1)
   call deallocator(mo_fock_mat,n_orbitals,n_orbitals)
   call deallocator(omega1,n_vir,n_occ) ! Omega_a,i
   call deallocator(omega2,n_ov_ov_packed,1)
!
!
end subroutine mlcc_drv
!
end module mlcc_drive
