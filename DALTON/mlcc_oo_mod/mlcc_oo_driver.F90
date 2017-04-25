subroutine mlcc_oo_driver
!  
   use ccs_class
   use input_output
!
   implicit none
!
   type(cc_singles) :: ccs
!
!  Create & open the main output file, used throughout the program
!
   call init_output_file
   open(unit=unit_output,file='mlcc.out',status='old',form='formatted')
!
!  Run the calculation
!
   call ccs % init
   call ccs % drv
!
end subroutine mlcc_oo_driver
