!> @file 
!> Contains the precision specifications
MODULE precision
  INTEGER, PARAMETER :: realk = KIND(1D0)
  integer, parameter :: long = selected_int_kind(10)
END MODULE precision
