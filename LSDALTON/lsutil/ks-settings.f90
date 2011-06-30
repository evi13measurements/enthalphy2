MODULE KS_settings
use matrix_module
! This module evaluates the Fock/Kohn-Sham matrix using chosen
! algorithm, matrix representation, etc.
!
!SAVE
  LOGICAL :: incremental_scheme = .FALSE.
  LOGICAL :: do_increment = .FALSE.
  CHARACTER(4) :: rho_curr = 'RHO1'
  CHARACTER(4) :: rho_prev = 'RHO2'
  CHARACTER(4) :: rho_curr2 = 'RHO3'
  CHARACTER(4) :: rho_prev2 = 'RHO4'
  TYPE(Matrix),save :: D0, F0, Ddiff
CONTAINS

subroutine ks_free_incremental_fock
use matrix_operations
implicit none
  incremental_scheme = .false.
  do_increment       = .false.
  call mat_free(D0)
  call mat_free(F0)
  call mat_free(Ddiff)
end subroutine ks_free_incremental_fock

subroutine ks_init_incremental_fock(nbast)
use matrix_operations
implicit none
integer :: nbast
  incremental_scheme = .true.
  do_increment       = .false.
  call mat_init(D0,nbast,nbast)
  call mat_init(F0,nbast,nbast)
  call mat_init(Ddiff,nbast,nbast)
end subroutine ks_init_incremental_fock



  subroutine svap_strings(st1,st2,len)
  implicit none
  integer :: len
  character(len) :: st1,st2,temp
  temp = st1
  st1  = st2
  st2  = temp
  end subroutine svap_strings
END MODULE KS_settings

SUBROUTINE get_incremental_settings(inc_scheme,do_inc)
use KS_settings
LOGICAL :: inc_scheme,do_inc
inc_scheme = incremental_scheme
do_inc = do_increment
END SUBROUTINE get_incremental_settings

