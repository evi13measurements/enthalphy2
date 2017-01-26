module mlcc_drive
!
!  Contains input routine called from mlcc_interface_input outside
!  and stores input in mlcc_data
!
contains
!
subroutine mlcc_input(input,i_count,lupri)
!
use mlcc_data
!
!  mlcc input routine
!  Author Rolf H. Myhre
!  December 2016
!
   implicit none
!
   integer, intent(in)           :: lupri, i_count
   character(len=*), intent(in)  :: input
!   
   character(len=254)            :: buffer
!   
   integer                       :: i, j, r_off, r_len
   integer                       :: trolo, kk
   integer, dimension(3)         :: kkx
   logical                       :: lupri_open
!
   inquire(unit=lupri,opened=lupri_open)
!
   if (.not. lupri_open) then
      write(*,*) 'lupri is not open'
      stop 
   end if
!
!  Say hi
   write(lupri,*) 
   write(lupri,*) 'Entered mlcc_input'
   write(lupri,*) 

!  initialize some defaults
!
   print_mlcc = 1 ! print level
!
   r_off = 1
   i = 0
   
   read(input(r_off:),*) buffer
   r_off = r_off + len(trim(buffer)) + 1
   i = i + 1
!
   if (trim(buffer) .eq. '*MLCC') then
!
      do while (i .lt. i_count)
!
!        Get the next keyword
         read(input(r_off:),*) buffer
         r_off = r_off + len(trim(buffer)) + 1
         i = i + 1

         select case (trim(buffer))
!        
            case('.MLACTIVE')
!        
               mlcc_active = .true.
!        
            case('.PRINT')
!        
               read(input(r_off:),*) buffer
               r_off = r_off + len(trim(buffer)) + 1
               read(buffer(1:),*) print_mlcc
               i = i + 1
!        
            case('.LALALA')
!        
               write(lupri,*) 'lalala'
!        
            case('.TROLOLOLO')
!        
               read(input(r_off:),*) buffer
               r_off = r_off + len(trim(buffer)) + 1
               read(buffer(1:),*) trolo
               i = i + 1
!
               write(lupri,*) 'trolo ', trolo
!        
            case('.KK')
!        
               read(input(r_off:),*) buffer
               r_off = r_off + len(trim(buffer)) + 1
               read(buffer(1:),*) kk
               i = i + 1
!
               if (kk .ge. 1) then
                  do j=1,kk
                     read(input(r_off:),*) buffer
                     r_off = r_off + len(trim(buffer)) + 1
                     i = i + 1
                     read(buffer(1:),*) kkx(j)
                  end do
               else
                  write(lupri,*) 'kk must be ge 1'
                  stop
               end if
!
               write(lupri,*) 'kk', kk
               write(lupri,*) 'kkx', kkx
!
            case default
               write(lupri,*) 'Keyword ', trim(buffer), ' not recognized in mlcc_input'
               stop
!        
         end select
!
      end do
!
   else 
!
      write(lupri,*) 'mlcc_input called with wrong input'
      stop 
!
   end if
!
!  Sanity checks
!
   if(print_mlcc .ge. 3) then
!
      write(lupri,*)
      write(lupri,*) 'Output from mlcc3_input'
      write(lupri,*) 'mlcc_active: ', mlcc_active
      write(lupri,*) 'print_mlcc : ', print_mlcc
      write(lupri,*)
!
   end if
!
end subroutine mlcc_input
!
subroutine mlcc_drv(work,lwork,lupri)
!
   use mlcc_types
   use mlcc_work
   use mlcc_data
   use mlcc_init
!
!  mlcc3 driver
!  Author Rolf H. Myhre
!  December 2016
!
   implicit none
!
   integer, intent(in)                    :: lupri !general output unit
   integer, intent(in)                    :: lwork !free space in work
!
   real(dp), intent(in), dimension(lwork) :: work !work static array
!
   ml_lupri = lupri
!
   write(lupri,*)
   write(lupri,*) 'In mlcc_drv'
   write(lupri,*)
!
   call work_init(mem,lupri)
   call hf_reader()
end subroutine mlcc_drv
!
end module mlcc_drive
