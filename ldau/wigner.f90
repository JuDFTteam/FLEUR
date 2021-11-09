module m_wigner

    !Implementations of wigner3j and wigner6j symbols
    !These implementations were translated into fortran
    !from the implementations of these symbols in the
    !sympy python package
    !https://www.sympy.org/en/index.html
    !The copyright notice of sympy is reproduced below
    !
    ! Copyright (c) 2006-2021 SymPy Development Team

    ! All rights reserved.

    ! Redistribution and use in source and binary forms, with or without
    ! modification, are permitted provided that the following conditions are met:

    ! a. Redistributions of source code must retain the above copyright notice,
    !     this list of conditions and the following disclaimer.
    ! b. Redistributions in binary form must reproduce the above copyright
    !     notice, this list of conditions and the following disclaimer in the
    !     documentation and/or other materials provided with the distribution.
    ! c. Neither the name of SymPy nor the names of its contributors
    !     may be used to endorse or promote products derived from this software
    !     without specific prior written permission.


    ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    ! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    ! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ! ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
    ! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    ! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    ! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    ! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    ! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    ! OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    ! DAMAGE.

    use m_juDFT

    implicit none

    interface wigner3j
        procedure :: wigner3jf, wigner3ji
    end interface

    interface wigner6j
        procedure :: wigner6jf, wigner6ji
    end interface

    private
    public :: wigner3j, wigner6j

    contains

    function factorial_list(maxfac)

        integer, intent(in) :: maxfac

        real factorial_list(0:maxfac)

        integer i

        factorial_list(0) = 1
        do i = 1, maxfac
            factorial_list(i) = factorial_list(i-1) * i
        enddo

    end function

    real function wigner3ji(l1,l2,l3,m1,m2,m3,fac)

        integer, intent(in) :: l1,l2,l3,m1,m2,m3
        real, optional, intent(in) :: fac(0:)

        wigner3ji = wigner3jf(real(l1),real(l2),real(l3),real(m1),real(m2),real(m3), fac=fac)

    end function

    real function wigner3jf(l1,l2,l3,m1,m2,m3,fac)

        real, intent(in) :: l1,l2,l3,m1,m2,m3
        real, optional, intent(in) :: fac(0:)

        integer :: tmin,tmax,t,f1,f2,f3,f4,f5
        integer :: maxfact
        real :: fac_list(0:int(max(l1 + l2 + l3 + 1, l1 + abs(m1), l2 + abs(m2),l3 + abs(m3)))) 
        real :: ressqrt 

        maxfact = int(max(l1 + l2 + l3 + 1, l1 + abs(m1), l2 + abs(m2),l3 + abs(m3)))
        wigner3jf = 0
  
        if (abs(int(l1*2)-l1*2)>1e-12&
            .or.abs(int(l2*2)-l2*2)>1e-12&
            .or.abs(int(l3*2)-l3*2)>1e-12) then
            call juDFT_error('l values should be integer or half-integer', calledby='wigner3j')
        endif
        if (abs(int(m1*2)-m1*2)>1e-12&
            .or.abs(int(m2*2)-m2*2)>1e-12&
            .or.abs(int(m3*2)-m3*2)>1e-12) then
            call juDFT_error('m values should be integer or half-integer', calledby='wigner3j')
        endif

        if (l1 + l2 - l3 < 0) return
        if (l1 - l2 + l3 < 0) return
        if (-l1 + l2 + l3 < 0) return

        if(abs(m1+m2+m3)>1e-12)  return
        if(abs(m1).gt.l1) return
        if(abs(m2).gt.l2) return
        if(abs(m3).gt.l3) return
        if(l3.lt.abs(l1-l2).or.l3.gt.l1+l2) return
  
        f1   = int(l3-l2+m1)
        f2   = int(l3-l1-m2)
        f3   = int(l1+l2-l3)
        f4   = int(l1-m1)
        f5   = int(l2+m2)

        tmin = int(max(0,-f1,-f2)) ! The arguments to fac (see below)
        tmax = int(min(f3,f4,f5))  ! must not be negative.

        ! The following line is only for testing and should be removed at a later time.      
        if(tmax-tmin .ne. int(min(l1+m1,l1-m1,l2+m2,l2-m2,l3+m3,l3-m3, &
                              l1+l2-l3,l1-l2+l3,-l1+l2+l3))) then
            call juDFT_error('Number of terms incorrect.', calledby='wigner3j')
        endif

        if (present(fac)) then
            if (size(fac) < maxfact) then
                call juDFT_error('Provided insufficient number of factorials', calledby='wigner3j')
            endif
            fac_list = fac(0:maxfact)
        else
            fac_list = factorial_list(maxfact)
        endif

        ressqrt = sqrt(real(fac_list(int(l1+l2-l3)))) * &
                  sqrt(real(fac_list(int(l1-l2+l3)))) * &
                  sqrt(real(fac_list(int(-l1+l2+l3)))) * &
                  sqrt(real(fac_list(int(l1-m1)))) * &
                  sqrt(real(fac_list(int(l1+m1)))) * &
                  sqrt(real(fac_list(int(l2-m2)))) * &
                  sqrt(real(fac_list(int(l2+m2)))) * &
                  sqrt(real(fac_list(int(l3-m3)))) * &
                  sqrt(real(fac_list(int(l3+m3)))) / &
                  sqrt(real(fac_list(int(l1 + l2 + l3 + 1))))

        if (tmin.le.tmax) then
            do t = tmin, tmax
                wigner3jf = wigner3jf + (-1)**t / &
                                   real(  fac_list(t)    * fac_list(f1+t) * fac_list(f2+t) &
                                         *fac_list(f3-t) * fac_list(f4-t) * fac_list(f5-t) )
            enddo
        endif
        wigner3jf = wigner3jf * ressqrt * (-1) ** int(l1 - l2 - m3)

    end function

    real function wigner6ji(l1,l2,l3,m1,m2,m3,fac)

        integer, intent(in) :: l1,l2,l3,m1,m2,m3
        real, optional, intent(inout) :: fac(0:)

        wigner6ji = wigner6jf(real(l1),real(l2),real(l3),real(m1),real(m2),real(m3), fac=fac)

    end function

    real function wigner6jf(l1,l2,l3,m1,m2,m3,fac)

        real, intent(in) :: l1,l2,l3,m1,m2,m3
        real, optional, intent(in) :: fac(0:)

        integer :: tmin,tmax,t,f1,f2,f3,f4,f5
        integer :: maxfact
        real :: prefactor
        real :: fac_list(0:max(int(min(l1 + l2 + l3 + m1, l1 + m1 + m2 + m3, l2 + l3 + m2 + m3))+ 1,&
                         int(l1 + l2 + l3 + m1), int(l1 + m1 + m2 + m3), &
                         int(l2 + l3 + m2 + m3) ))

        wigner6jf = 0.0
        prefactor = (-1) ** int(l1 + l2 + m1 + m2) * &
                    big_delta_coeff(l1,l2,m2, fac=fac) * &
                    big_delta_coeff(l3,m1,m2, fac=fac) * &
                    big_delta_coeff(l1,l3,m3, fac=fac) * &
                    big_delta_coeff(l2,m1,m3, fac=fac)

        tmin = int(max(l1 + l2 + m2, l3 + m1 + m2, l1 + l3 + m3, l2 + m1 + m3))
        tmax = int(min(l1 + l2 + l3 + m1, l1 + m1 + m2 + m3, l2 + l3 + m2 + m3))

        maxfact = max(tmax + 1, int(l1 + l2 + l3 + m1), int(l1 + m1 + m2 + m3), &
                      int(l2 + l3 + m2 + m3))
        
        if (present(fac)) then
            if (size(fac) < maxfact) then
                call juDFT_error('Provided insufficient number of factorials', calledby='wigner3j')
            endif
        else
            fac_list = factorial_list(maxfact)
        endif
        
        wigner6jf = 0.0
        do t= tmin, tmax
            wigner6jf = wigner6jf + (-1)**t * fac_list(t+1) / (&
                        fac_list(int(t-l1-l2-m2)) * &
                        fac_list(int(t-l3-m1-m2)) * &
                        fac_list(int(t-l1-l3-m3)) * &
                        fac_list(int(t-l2-m1-m3)) * &
                        fac_list(int(l1+l2+l3+m1-t)) * &
                        fac_list(int(l1+m1+m2+m3-t)) * &
                        fac_list(int(l2+l3+m2+m3-t)))
        enddo

        wigner6jf = prefactor * wigner6jf * (-1) ** int(l1 + l2 + l3 + m1)

    end function

    real function big_delta_coeff(l1,l2,l3, fac)

        real, intent(in) :: l1,l2,l3
        real,optional, intent(in) :: fac(0:)

        integer :: maxfact
        real :: fac_list(0:int(max(l1 + l2 - l3, l1 + l3 - l2, l2 + l3 - l1, l1 + l2 + l3 + 1)))

        big_delta_coeff = 0.0

        if (abs(int(l1 + l2 - l3) - (l1 + l2 - l3))>1e-12.or.&
            abs(int(l1 - l2 + l3) - (l1 - l2 + l3))>1e-12.or.&
            abs(int(-l1 + l2 + l3) - (-l1 + l2 + l3))>1e-12) then
            call juDFT_error('l values should be integer or half-integer', calledby='big_delta_coef')
        endif
        if ((l1 + l2 - l3) < 0) return
        if ((l1 - l2 + l3) < 0) return
        if ((-l1 + l2 + l3) < 0) return

        maxfact = int(max(l1 + l2 - l3, l1 + l3 - l2, l2 + l3 - l1, l1 + l2 + l3 + 1))
        if (present(fac)) then
            if (size(fac) < maxfact) then
                Call juDFT_error('Provided insufficient number of factorials', calledby='wigner3j')
            endif
        else
            fac_list = factorial_list(maxfact)
        endif

        big_delta_coeff = sqrt(fac_list(int(l1 + l2 - l3))) * &
                          sqrt(fac_list(int(l1 + l3 - l2))) * &
                          sqrt(fac_list(int(l2 + l3 - l1))) / &
                          sqrt(fac_list(int(l1 + l2 + l3 + 1)))

    end function
    
end module m_wigner