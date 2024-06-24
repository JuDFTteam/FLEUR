module m_map_to_unit
contains    
   function map_to_unit(x) result(y)
    implicit none 
    real, intent(in) :: x(:)
    real             :: y(size(x))

    y = x

    ! everything close to 0 or 1 get's mapped to 0 and 1
    where (abs(y - anint(y)) < 1e-6) y = anint(y)
    ! map to 0 -> 1 interval
    y = y - floor(y)
 end function map_to_unit
end module m_map_to_unit
