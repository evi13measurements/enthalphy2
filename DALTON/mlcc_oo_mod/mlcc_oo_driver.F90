subroutine mlcc_oo_driver
!  
   use hf_class
   use printing
!
   type(hartree_fock) :: hf
!
   call init_output_file
!
   open(unit=unit_output,file='mlcc.out',status='old',form='formatted')
   write(unit_output,*) 'Called mlcc_oo_drv'
!
   call hf % init
!
end subroutine mlcc_oo_driver
