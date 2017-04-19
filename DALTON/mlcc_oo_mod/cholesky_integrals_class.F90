module cholesky_integrals_class
!
   type cholesky_integrals
      integer           :: n_J  = 0
      integer           :: n_ao = 0
      integer           :: n_mo = 0
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
      use dalton_interface
!
      implicit none
!
      class(cholesky_integrals)   :: cholesky
!
      cholesky % integral_program = 'DALTON    '
!
      if (cholesky % integral_program .eq. 'DALTON    ') then
!
         call dalton_inteface_drv(cholesky)
!
      endif
!
   end subroutine init_cholesky_integrals
!
end module cholesky_integrals_class
