module mlcc_data
!
!
!  mlcc3 types
!  authors henrik koch and rolf h. myhre
!  december 2014
!
!  purpose: store variables required througout mlcc3
!
   use mlcc_types
!   
   implicit none
!
!
!  model variable from old dalton
   character(10)                    :: model ='ccsd      '
!   
!  energy or response vector
   logical                          :: resp_option = .false.
!   
!  integer variables
   integer                          :: n_orbitals, n_basis, n_lambda, n_occ, n_vir,n_J,n_reduced
   integer                          :: n_t1am, n_t2am, n_t2am_pack
   integer                          :: n_v_2, n_v_3, n_basis_2, n_basis_2_pack, n_bas_orb
   integer                          :: n_ao_ints, n_ov, n_vv, n_oo, n_oo_packed, n_oov, n_ovv
   integer                          :: n_ooo,n_vv_packed,n_ov_ov_packed
!
!  info read from file  
   real(dp)                         :: nuclear_potential,scf_energy
!   
!  basic variable
   logical                          :: mlcc3_active
   logical                          :: mlcc3_nrg_spa
   logical                          :: mlcc3_nrg_gen
   integer                          :: n_active, n_general
   integer                          :: n_occ_inp,n_vir_inp
   integer                          :: n_gen_o_inp,n_gen_v_inp
!
!  various pointers
!
!  packed omega vectors from ccsd
   real(dp), dimension(:,:), pointer  :: omega1            => null()
   real(dp), dimension(:,:), pointer  :: omega2            => null()
!   
!  t1 and t2 amplitudes, t2 packed
   real(dp), dimension(:,:), pointer  :: t1am              => null()
   real(dp), dimension(:,:), pointer  :: t2am              => null()
!
!  hf orbital coefficients and fock diagonal elements
   real(dp), dimension(:,:), pointer  :: orb_coefficients  => null()
   real(dp), dimension(:,:), pointer  :: fock_diagonal     => null()
!   
!  the mo fock matrix, standard and t1- and c1-transformed
   real(dp), dimension(:,:), pointer  :: mo_fock_mat       => null()
!
   integer                          :: luprint !general output unit
   integer                          :: mem = 10000000
end module mlcc_data

module mlcc_input_data
!
!  purpose: store input variables from interface
!
   use mlcc_types
!   
   logical  :: mlcc_active = .false.
   integer  :: print_mlcc = 1
!   
!
end module mlcc_input_data

