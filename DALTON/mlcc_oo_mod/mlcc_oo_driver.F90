subroutine mlcc_oo_driver
!  
   use ccsd_class
   use workspace
   use input_output
!
   implicit none
!
   type(cc_singles_doubles) :: ccsd
!
!  Set up workspace controller
!
   call work_init
!
!  Create & open the main output file, used throughout the program
!
   call init_output_file
   open(unit=unit_output,file='mlcc.out',status='old',form='formatted')
!
!  Run the calculation
!
   call ccsd % init
   call ccsd % drv
!
end subroutine mlcc_oo_driver
