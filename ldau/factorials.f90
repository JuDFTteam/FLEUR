module m_factorials

    implicit none

    contains

    real function doubleFactorial(n)

        !Calculates the double factorial n!!

        integer, intent(in) :: n

        integer upper_bound, k

        if (mod(n,2) == 1) then
            upper_bound = n+1/2
        else
            upper_bound = n/2
        endif

        doubleFactorial = 1
        do k = 1, upper_bound
            if (mod(n,2) == 1) then
                doubleFactorial = doubleFactorial * (2*k-1)
            else
                doubleFactorial = doubleFactorial * 2*k
            endif
        enddo

    end function

    real function factorial(n)

        !Calculates the factorial n!

        integer, intent(in) :: n

        integer k

        factorial = 1
        do k = 1, n
            factorial = factorial * k
        enddo

    end function

end module