module ccsd_class
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!                                                                 
!      Coupled cluster singles and doubles (CCSD) class module                                 
!   Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017  
!                                                                 
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Use general tools 
!
   use types
   use utils
   use workspace
   use input_output
!
!  Use ancestor class module (CCS)
!
   use ccs_class
!
   implicit none 
!
!  -::- Definition of the CCSD class -::- 
!
   type, extends(cc_singles) :: cc_singles_doubles
!
!     Amplitude attributes
!
      integer(i15) :: n_t2am = 0                    ! Number of doubles amplitudes
      real(dp), dimension(:,:), allocatable :: t2am ! Doubles amplitude vector
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_cc_singles_doubles
      procedure :: drv  => drv_cc_singles_doubles
!
!     Initialization routine for the (singles, doubles) amplitudes
!
      procedure :: initialize_amplitudes => initialize_amplitudes_cc_singles_doubles
!
!     Routine to calculate the energy from the current amplitudes
!
      procedure :: get_energy => get_energy_cc_singles_doubles 
!
!     Routine to set a number of orbital shorthands (see below)
!
      procedure :: set_orbital_shorthands => set_orbital_shorthands_cc_singles_doubles
!
   end type cc_singles_doubles
!
!  -::- Module variables and routines not belonging to the class -::-
!
!  Shorthands for orbital information used extensively in the
!  matrix multiplications of coupled cluster theory 
!
   integer(i15) :: n_o            ! n_occ
   integer(i15) :: n_v            ! n_vir
   integer(i15) :: n_ov           ! n_occ * n_vir
   integer(i15) :: n_oo           ! n_occ^2 
   integer(i15) :: n_vv           ! n_vir^2 
   integer(i15) :: n_oo_packed    ! n_occ * (n_occ + 1) / 2
   integer(i15) :: n_vv_packed    ! n_vir * (n_vir + 1) / 2
   integer(i15) :: n_oov          ! n_occ^2 * n_vir 
   integer(i15) :: n_ovv          ! n_occ * n_vir^2 
   integer(i15) :: n_ooo          ! n_occ^3 
   integer(i15) :: n_ov_ov_packed ! n_ov * (n_ov + 1) / 2
!
!
contains
!
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!              Initialization and driver routines
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
   subroutine init_cc_singles_doubles(wfn)
!
!     Initialize CCSD object
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Performs the following tasks:
!
!     1. Sets HF orbital and energy information by reading from file (read_hf_info)
!     2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!     3. Allocates the Fock matrix and sets it to zero
!     4. Sets orbital shorthands (n_o, n_oo, n_ov, etc.)
!     5. Initializes the amplitudes (sets their initial values and associated variables)
!
      implicit none 
!
      class(cc_singles_doubles) :: wfn

!     Read Hartree-Fock info from SIRIUS
!
      call wfn % read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wfn % read_transform_cholesky 
!
!     Set orbital shorthands (n_o, n_oo, n_ov, etc.) that are 
!     members of the module (not the class)
!
      call wfn % set_orbital_shorthands
!
!     Initialize (singles and doubles) amplitudes
!
      call wfn % initialize_amplitudes
!
   end subroutine init_cc_singles_doubles
!
!
   subroutine drv_cc_singles_doubles(wfn)
!
      implicit none 
!
      class(cc_singles_doubles) :: wfn
!
! To do...
!
   end subroutine drv_cc_singles_doubles
!
!
   subroutine initialize_amplitudes_cc_singles_doubles(wfn)
!
!     Initialize Amplitudes (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Allocates the amplitudes, sets them to zero, calculates
!     the number of amplitudes, and sets the doubles amplitudes
!     to the perturbative MP2 estimate.
!
      implicit none 
!
      class(cc_singles_doubles) :: wfn
!
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_ia_jb
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0
      integer(i15) :: ai = 0, bj = 0, ia = 0, jb = 0, aibj = 0
!
!     Calculate the number of singles and doubles amplitudes
!
      wfn % n_t1am = (wfn % n_occ)*(wfn % n_vir) 
      wfn % n_t2am = (wfn % n_t1am)*(wfn % n_t1am + 1)/2
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator(wfn % t1am, wfn % n_t1am, 1)
      wfn % t1am = zero
!
!     Allocate the doubles amplitudes and set to zero
!
      call allocator (wfn % t2am, wfn % n_t2am, 1)
      wfn % t2am = zero
!
!
!     -::- Initialize the doubles amplitudes to the MP2 estimate -::-
!
!
!     Allocate L_ia_J and g_ia_jb
!
      call allocator(L_ia_J, n_ov, wfn % n_J)
      L_ia_J = zero
!
      call allocator(g_ia_jb, n_ov, n_ov)
      g_ia_jb = zero
!
!     Calculate g_ia_jb = g_iajb
!
      call dgemm('N','T',    &
                  n_ov,      &
                  n_ov,      &
                  wfn % n_J, &
                  one,       &
                  L_ia_J,    &
                  n_ov,      &
                  L_ia_J,    &
                  n_ov,      &
                  zero,      &
                  g_ia_jb,   &
                  n_ov)
!
!     Set the doubles amplitudes
!
      do a = 1, wfn % n_vir
         do b = 1, wfn % n_vir
            do i = 1, wfn % n_occ
               do j = 1, wfn % n_occ
!
!                 Get necessary indices
!
                  ai = index_two(a, i, wfn % n_vir)
                  bj = index_two(b, j, wfn % n_vir)
                  ia = index_two(i, a, wfn % n_occ)
                  jb = index_two(j, b, wfn % n_occ)
!
!                 Set the doubles indices
!
                  if (ai .le. bj) then
!
                     aibj = index_packed(ai,bj)
                     wfn % t2am(aibj, 1) = - g_ia_jb(ia,jb)/(wfn % fock_diagonal(wfn % n_occ + a, 1) + &
                                                             wfn % fock_diagonal(wfn % n_occ + b, 1) - &
                                                             wfn % fock_diagonal(i, 1) - &
                                                             wfn % fock_diagonal(j, 1))
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine initialize_amplitudes_cc_singles_doubles
!
!
   subroutine get_energy_cc_singles_doubles(wfn)
!
!     Get Energy (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the CCSD energy
!
      implicit none 
!
      class(cc_singles_doubles) :: wfn 
!
      real(dp), dimension(:,:), allocatable :: L_ia_J  ! L_ia^J
      real(dp), dimension(:,:), allocatable :: g_ia_jb ! g_iajb
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0, ai = 0
      integer(i15) :: bj = 0, aibj = 0, ia = 0, jb = 0, ib = 0, ja = 0
!
!     Allocate the Cholesky vector L_ia_J = L_ia^J and set to zero 
!
      call allocator(L_ia_J, (wfn % n_occ)*(wfn % n_vir), wfn % n_J)
      L_ia_J = zero
!
!     Get the Cholesky vector L_ia_J 
!
      call wfn % get_cholesky_ia(L_ia_J)
!
!     Allocate g_ia_jb = g_iajb and set it to zero
!
      call allocator(g_ia_jb, n_ov, n_ov)
      g_ia_jb = zero
!
!     Calculate the integrals g_ia_jb from the Cholesky vector L_ia_J 
!
      call dgemm('N','T',    &
                  n_ov,      &
                  n_ov,      &
                  wfn % n_J, &
                  one,       &
                  L_ia_J,    &
                  n_ov,      &
                  L_ia_J,    &
                  n_ov,      &
                  zero,      &
                  g_ia_jb,   &
                  n_ov)
!
!     Deallocate the Cholesky vector L_ia_J 
!
      call deallocator(L_ia_J, n_ov, wfn % n_J)
!
!     Set the initial value of the energy 
!
      wfn % energy = wfn % scf_energy
!
!     Add the correlation energy E = E + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!
      do i = 1, wfn % n_occ
         do a = 1, wfn % n_vir
            do j = 1, wfn % n_occ
               do b = 1, wfn % n_vir
!
!                 Calculate the necessary indices 
!
                  ai   = index_two(a, i, wfn % n_vir)
                  bj   = index_two(b, j, wfn % n_vir)
!
                  aibj = index_packed(ai, bj)
!
                  ia   = index_two(i, a, wfn % n_occ)
                  jb   = index_two(j, b, wfn % n_occ)
!
                  ib   = index_two(i, b, wfn % n_occ)
                  ja   = index_two(j, a, wfn % n_occ)
!
!                 Add the correlation energy 
!
                  wfn % energy = wfn % energy + & 
                                 (wfn % t2am(aibj,1) + (wfn % t1am(a,i))*(wfn % t1am(b,j)))*&
                                 (two*g_ia_jb(ia,jb) - g_ia_jb(ib,ja))
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ia_jb
!
      call deallocator(g_ia_jb, n_ov, n_ov)
!
   end subroutine get_energy_cc_singles_doubles
!
!
   subroutine set_orbital_shorthands_cc_singles_doubles(wfn)
!
!     Set Orbital Shorthands (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017   
!
!     Sets orbital shorthands useful for matrix multiplications in the
!     coupled cluster code 
!
      implicit none 
!
      class(cc_singles_doubles), intent(in) :: wfn
!
      n_o = wfn % n_occ
      n_v = wfn % n_vir 
!
      n_ov = n_o*n_v
      n_oo = n_o*n_o 
      n_vv = n_v*n_v 
!
      n_oo_packed = n_o*(n_o+1)/2
      n_vv_packed = n_v*(n_v+1)/2
!
      n_oov = n_o*n_o*n_v 
      n_ovv = n_o*n_v*n_v 
      n_ooo = n_o*n_o*n_o
!
      n_ov_ov_packed = n_ov*(n_ov+1)/2
!
   end subroutine set_orbital_shorthands_cc_singles_doubles
!
!
end module ccsd_class
