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
   use mlcc_fock
   use mlcc_energy
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
   luprint = lupri 
!
   write(luprint,*) 'Calling oo driver from old driver'
   call flshfo(luprint)
!
   call mlcc_oo_driver
!
   write(luprint,*) 'Beginning old driver'
   call flshfo(luprint)
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
   call allocator(F_i_j,n_occ,n_occ)
   call allocator(F_i_a,n_occ,n_vir)
   call allocator(F_a_i,n_vir,n_occ)
   call allocator(F_a_b,n_vir,n_vir)
   call mlcc_get_fock
   call t2_init
!
!  Calculate the omega vector 
!
!   call mlcc_omega_calc
!
! !
! !  Start the coupled cluster solver for the ground state energy
! !
    call mlcc_energy_drv
! !
   call deallocator(t2am,n_ov_ov_packed,1)
   call deallocator(t1am,n_vir,n_occ)
   call deallocator(fock_diagonal,n_orbitals,1)
   call deallocator(orb_coefficients,n_lambda,1)
   call deallocator(omega1,n_vir,n_occ) ! Omega_a,i
   call deallocator(omega2,n_ov_ov_packed,1)
   call deallocator(F_i_j,n_occ,n_occ)
   call deallocator(F_i_a,n_occ,n_vir)
   call deallocator(F_a_i,n_vir,n_occ)
   call deallocator(F_a_b,n_vir,n_vir)
!
!
end subroutine mlcc_drv
!
end module mlcc_drive
