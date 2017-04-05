module class_shape
    implicit none
!
!
    type :: shape 
        integer :: color = 1
        real    :: area = 0.0d0
    contains
        procedure :: get_area => shape_area
        procedure :: get_color
    end type shape
!
contains
    subroutine shape_area(this)
!
    implicit none
!
        class(shape), intent(inout) :: this
        write(*,*)'Area is:',this%area
!
    end subroutine shape_area
!
    subroutine get_color(this)
!
        implicit none
!
        class(shape), intent(inout) :: this
        write(*,*)'The color is:',this%color
    end subroutine get_color
end module class_shape

module class_triangle
    use class_shape
    implicit none
!
    type, extends(shape) :: triangle
!
        real :: baseline
        real :: height
    contains
        procedure :: get_area => triangle_area
    end type triangle
!
contains
    subroutine  triangle_area(this)
        implicit none
!
        class(triangle), intent(inout) :: this
!
        this%area=this%baseline*this%height/2.0d0
!
        write(*,*)'Area is:',this%area
!
    end subroutine  triangle_area
end module class_triangle

program test_shape
    use class_triangle
    implicit none
!
    type(triangle)   :: s
    s=triangle(1,0.0d0,3.0d0,1.0d0)
    call s%get_area
    call s%get_color
!
end program test_shape
