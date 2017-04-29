module types
!
!  Types module
!  Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!  Defines the real and integer types used throughout the module
!
   implicit none
!
!  Kind numbers for real and integers
!  Usage: real(dp) :: foo, integer(i15) :: foo_int, etc.
!
   integer, parameter :: sp  = selected_real_kind(6,37)
   integer, parameter :: dp  = selected_real_kind(15,307)
   integer, parameter :: qp  = selected_real_kind(33,4931)
   integer, parameter :: i15 = selected_int_kind(15)
!   
!  Integers as reals 
!
   real(dp), parameter :: zero = 0.0D0, one = 1.0D0, two = 2.0D0, half = 0.5D0
   real(dp), parameter :: three = 3.0D0, four = 4.0D0, five = 5.0D0, six = 6.0D0
!
end module types
