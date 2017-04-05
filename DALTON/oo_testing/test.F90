program test_shape
    use class_shape
    implicit none
!
    type(shape)   :: s
    s=shape(.false.)
    call s%check
!
end program test_shape