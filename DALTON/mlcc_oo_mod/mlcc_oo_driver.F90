subroutine mlcc_oo_driver
!  
!  Model classes
!
   use hf_class
   use ccs_class
   use ccsd_class
   use cc2_class
!
!  Memory handling and IO modules
!
   use workspace
   use input_output
!
   implicit none
!
!  Declare wavefunction
!
   type(ccsd) :: wf 
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
!  Start the calculation
!
   call wf % init
   call wf % drv
!
!  Close main output file
!
   close(unit_output)
!
end subroutine mlcc_oo_driver
