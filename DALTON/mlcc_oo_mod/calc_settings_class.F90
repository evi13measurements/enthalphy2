module calc_settings_class 
!
!
!                         Calculation settings class module                                 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017         
!
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  :::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the calc_settings class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: settings 
!
      real(dp) :: energy_threshold = 1.0D-06
      real(dp) :: ampeqs_threshold = 1.0D-06
!
      integer(i15) :: ampeqs_max_iterations = 25
!
   end type settings                                                                            
!
!
contains
!
end module calc_settings_class