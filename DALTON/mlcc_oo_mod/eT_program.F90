program eT_program
!
!
!
!                 Coupled cluster module eT - Main program                                
!         Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017         
!                                                                           
!
   use input_output
!
   use hf_class
   use ccs_class
   use cc2_class
   use ccsd_class
!
   implicit none
!
   type(cc2), allocatable  :: cc2
   type(ccsd), allocatable :: ccsd
!
   class(hf), pointer :: wf
!
!  Open input file
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Reading method section of input file -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set method - Initialize object
!
!  call method_reader(STRING)
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Reading calculation section of input file -::- 
!  :::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set calculation specifications
!
!  call calculation_reader(wf%tasks)
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Reading settings section of input file -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set calculation settings settings
!
!  call settings_reader(wf%settings) 
!
!  Close input file
!
!  :::::::::::::::::::::::::::::
!  -::- Prepare output file -::- 
!  :::::::::::::::::::::::::::::
!
!  :::::::::::::::::::::::::
!  -::- Run Calculation -::- 
!  :::::::::::::::::::::::::
!
!  call wf%drv()
!
!  :::::::::::::::::::::::::::
!  -::- Close output file -::- 
!  :::::::::::::::::::::::::::
!
end program eT_program