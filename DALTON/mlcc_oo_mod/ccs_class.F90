module ccs_class
!
   use input_output
   use workspace
   use mlcc_types
   use hf_class
!
   type, extends(hartree_fock) :: cc_singles
!
      integer(i15) :: n_t1am = 0
      real(dp), dimension(:,:), allocatable :: t1am ! Singles amplitudes
!
   contains 
!
      procedure :: init => init_cc_singles
      procedure :: drv  => drv_cc_singles
!
   end type cc_singles
!
contains
!
   subroutine init_cc_singles(wavefn)
!
      implicit none 
!
      class(cc_singles) :: wavefn
!
      write(unit_output,*) 'Initializing CCS...'
!
!     Initialize the Hartree-Fock-specific quantities
!
      write(unit_output,*) '- Calling HF initializer'
      call init_hartree_fock(wavefn)
!
!     Calculate the number of singles amplitudes
!
      write(unit_output,*) '- Initializing singles amplitudes'
      wavefn % n_t1am = (wavefn % n_occ)*(wavefn % n_vir) 
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator ( wavefn % t1am, wavefn % n_t1am, 1)
      wavefn % t1am = zero
!
   end subroutine init_cc_singles
!
   subroutine drv_cc_singles(wavefn)
!
      implicit none 
!
      class(cc_singles) :: wavefn
!
      write(unit_output,*) 'CCS driver has begun...'
!
   end subroutine drv_cc_singles
!
end module ccs_class
