module cholesky_integrals_class
!
   use mlcc_types
!
   type cholesky_integrals
      integer(i15)      :: n_J       = 0
      integer(i15)      :: n_ao      = 0
      integer(i15)      :: n_mo      = 0
      character(len=10) :: integral_program
   contains
!	
      procedure :: init => init_cholesky_integrals
!
   end type cholesky_integrals
!
contains
!
   subroutine init_cholesky_integrals(cholesky)
!
      implicit none
!
      class(cholesky_integrals)   :: cholesky
!
      cholesky % integral_program = 'DALTON    '
!
      if (cholesky % integral_program .eq. 'DALTON    ') then
!
         call dalton_interface_driver(cholesky)
!
      endif
!
   end subroutine init_cholesky_integrals
!
end module cholesky_integrals_class
