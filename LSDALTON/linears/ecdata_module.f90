module EcData
use precision
!grand canonical electronic configuration data
!for most atomic numbers 1 to 86 (H to Rn)

TYPE ElementEcData
 real(realk) :: occ_s(7)
 real(realk) :: occ_p(15)
 real(realk) :: occ_d(15)
 real(realk) :: occ_f(7)

 integer     :: a
 integer     :: nocc_s
 integer     :: nocc_p
 integer     :: nocc_d
 integer     :: nocc_f
END TYPE 



Type(ElementEcData) :: ElementTable(86)

contains

subroutine EcData_init
implicit none
integer :: i
real(realk) :: j

!Zero out
 do i=1,86
  ElementTable(i)%a=0
  ElementTable(i)%occ_s=0d0;
  ElementTable(i)%occ_p=0d0;
  ElementTable(i)%occ_d=0d0;
  ElementTable(i)%occ_f=0d0;
  ElementTable(i)%nocc_s=0;
  ElementTable(i)%nocc_p=0;
  ElementTable(i)%nocc_d=0;
  ElementTable(i)%nocc_f=0;
 enddo

!Hydrogen
 ElementTable(1)%a=1
 ElementTable(1)%occ_s(1)=0.5d0 
 ElementTable(1)%nocc_s=1


!Helium
 ElementTable(2)%a=2
 ElementTable(2)%occ_s(1)=1d0 
 ElementTable(2)%nocc_s=1

 ! Li, Be
 j=0d0
 do i=3,4
 j=j+1d0
 ElementTable(i)=ElementTable(2)
 ElementTable(i)%a=i
 ElementTable(i)%occ_s(2)=j/2d0
 ElementTable(i)%nocc_s=2
 enddo


! B C N O F Ne
 j = 0d0
 do i=5,10
 j = j+  1d0
 ElementTable(i)=ElementTable(4)
 ElementTable(i)%a=i
 ElementTable(i)%occ_p(1:3)=(/ j/6d0, j/6d0, j/6d0 /)
 ElementTable(i)%nocc_p=3
 enddo

! Na, Mg
 j=0d0
 do i=11,12
 j=j+1d0
 ElementTable(i)=ElementTable(10)
 ElementTable(i)%a=i
 ElementTable(i)%occ_s(3)=j/2d0
 ElementTable(i)%nocc_s=3
 enddo


! Al, Si, P, S, Cl, Ar
 j=0d0
 do i=13,18
 j=j+1
 ElementTable(i)=ElementTable(12)
 ElementTable(i)%a=i
 ElementTable(i)%occ_p(4:6)=(/ j/6d0, j/6d0, j/6d0 /)
 ElementTable(i)%nocc_p=6
 enddo

! K, Ca
 j=0d0
 do i=19,20
 j=j+1d0
 ElementTable(i)=ElementTable(18)
 ElementTable(i)%a=i
 ElementTable(i)%occ_s(4)=j/2d0
 ElementTable(i)%nocc_s=4
 enddo


!Sc Ti V Cr Mn Fe Co Ni Cu Zn

 j=0d0
 do i=21,30
 j = j+ 1d0
 ElementTable(i)=ElementTable(20)
 ElementTable(i)%a=i
 ElementTable(i)%occ_d(1:5)=(/ j/10d0, j/10d0, j/10d0, j/10d0, j/10d0 /)
 ElementTable(i)%nocc_d=5
 enddo

! Cr
 ElementTable(24)%occ_s(4)=0.5d0
 ElementTable(24)%occ_d(1:5)=(/ 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0 /)

! Cu
 ElementTable(29)%occ_s(4)=0.5d0
 ElementTable(29)%occ_d(1:5)=(/ 1d0, 1d0, 1d0, 1d0, 1d0 /)
 

! Ga Ge As Se Br Kr
 j=0d0
 do i=31,36
 j = j+1
 ElementTable(i)=ElementTable(30)
 ElementTable(i)%a=i
 ElementTable(i)%occ_p(7:9)=(/ j/6d0, j/6d0, j/6d0 /)
 ElementTable(i)%nocc_p=9
 enddo 

! Rb, Sr
 j=0d0
 do i=37,38
 j=j+1d0
 ElementTable(i)=ElementTable(36)
 ElementTable(i)%a=i
 ElementTable(i)%occ_s(5)=j/2d0
 ElementTable(i)%nocc_s=5
 enddo

! Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd
 j=2d0
 do i=39,48
 j=j+1d0
 ElementTable(i)=ElementTable(36)
 ElementTable(i)%a=i
 ElementTable(i)%occ_s(5)=j/12
 ElementTable(i)%occ_d(6:10)=(/ j/12d0, j/12d0, j/12d0, j/12d0, j/12d0 /)
 ElementTable(i)%nocc_d=10
 ElementTable(i)%nocc_s=5
 enddo


! In Sn Sb Te I Xe
  j=0d0
 do i=49,54
 j = j+1
 ElementTable(i)=ElementTable(48)
 ElementTable(i)%a=i
 ElementTable(i)%occ_p(10:12)=(/ j/6d0, j/6d0, j/6d0 /)
 ElementTable(i)%nocc_p=12
 enddo 

! Cs Ba
 j=0d0
 do i=55,56
 j=j+1d0
 ElementTable(i)=ElementTable(54)
 ElementTable(i)%a=i
 ElementTable(i)%occ_s(6)=j/2d0
 ElementTable(i)%nocc_s=6
 enddo

! La - Hg
 j=0d0
 do i=57,80
 j = j + 1d0
 ElementTable(i)=ElementTable(56)
 ElementTable(i)%a=i
 ElementTable(i)%occ_d(11:15)=(/ j/24d0, j/24d0, j/24d0, j/24d0, j/24d0 /)
 ElementTable(i)%occ_f(1:7) = (/ j/24d0, j/24d0, j/24d0, j/24d0, j/24d0, j/24d0, j/24d0 /)
 ElementTable(i)%nocc_d=15
 ElementTable(i)%nocc_f=7
 enddo

! Tl - Rn
  j=0d0
 do i=81,86
 j = j+1
 ElementTable(i)=ElementTable(80)
 ElementTable(i)%a=i
 ElementTable(i)%occ_p(13:15)=(/ j/6d0, j/6d0, j/6d0 /)
 ElementTable(i)%nocc_p=15
 enddo 
 
#if 0
 do i=1,86

  j=0d0
  do k=1,ElementTable(i)%nocc_s
   j = j + ElementTable(i)%occ_s(k)
  enddo

  do k=1,ElementTable(i)%nocc_p
   j = j + ElementTable(i)%occ_p(k)
  enddo

  do k=1,ElementTable(i)%nocc_d
   j = j + ElementTable(i)%occ_d(k)
  enddo

  do k=1,ElementTable(i)%nocc_f
   j = j + ElementTable(i)%occ_f(k)
  enddo


  write (*,*), ElementTable(i)%a, ElementTable(i)%a/2d0, j
  if ((abs(j - ElementTable(i)%a/2d0)).ge.1d-6) then
     stop
  endif

 enddo
#endif

end subroutine EcData_init
end module EcData
