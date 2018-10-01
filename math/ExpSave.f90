module m_ExpSave

  contains

  pure real function exp_save( x )
    ! replace exp by a function that does not under/overflow dw09
    implicit none
    real, intent(in)      :: x
    real, parameter       :: maxexp = log( 2.0 ) * maxexponent( 2.0 )
    real, parameter       :: minexp = log( 2.0 ) * minexponent( 2.0 )

    if ( abs( x ) > minexp .and. abs( x ) < maxexp ) then
      exp_save = exp( x )
    else
      if ( x > 0 ) then
        if ( x > minexp ) then
          exp_save = exp( maxexp )
        else
          exp_save = exp( minexp )
        endif
      else
        if ( - x > minexp ) then
          exp_save = exp( - maxexp )
        else
          exp_save = exp( - minexp )
        endif
      endif
    endif

  end function exp_save

end module m_ExpSave
