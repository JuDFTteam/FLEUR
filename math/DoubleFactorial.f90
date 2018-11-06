module m_DoubleFactorial

implicit none


contains

  
real function DoubleFactorial( n_upper, n_lower )
! calculates ( 2 * n_upper + 1 ) !! / ( 2 * n_lower + 1 ) !! or just ( 2 * n_upper + 1 ) !!, if n_lower is not present

  integer           :: n_upper
  integer, optional :: n_lower
  integer           :: i, i_lower

  i_lower = 1
  if( present(n_lower) ) i_lower = n_lower + 1

  DoubleFactorial = 1.
  do i = i_lower, n_upper
    DoubleFactorial = DoubleFactorial * ( 2 * i + 1 )
  end do

end function DoubleFactorial


end module m_DoubleFactorial
