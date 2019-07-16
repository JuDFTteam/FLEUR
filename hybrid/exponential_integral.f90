! Calculate the exponential integral using the algorithm of
! [1] Tseng, Lee, Journal of Hydrology, 205 (1998) 38-51
module m_exponential_integral

   implicit none

   real, parameter :: series_laguerre = 4.0

contains

   ! Calculate the exponential integral E_1(x):
   !
   !        inf
   !          /     -t
   !         |     e
   ! E (x) = | dt -----
   !  1      |      t
   !        /
   !         x
   !
   ! Input:  arg - position at which exponential integral is evaluated (arg > 0)
   ! Output: res - E_1(arg)
   pure subroutine calculateExponentialIntegral(arg, res)

      implicit none

      real, intent(in)  :: arg
      real, intent(out) :: res

      ! For arguments smaller than 4 the series expansion is used
      if (arg < series_laguerre) then
         res = seriesExpansion(arg)

         ! otherwise a Gauss-Laguerre expansion is better
      else
         res = exp(-arg)*gauss_laguerre(arg)
      endif

   end subroutine calculateExponentialIntegral

   ! Series expansion of the exponential integral
   !
   !                          n_cut
   !                          -----     n  n
   !                           \    (-1)  x
   ! E (x) = -gamma - ln(x) -   )   --------
   !  1                        /     n * n!
   !                          -----
   !                          n = 1
   !
   ! where gamma is the Euler constant.
   ! n_cut is set to 25
   ! Input: arg - argument for which the exponential integral is approximated
   ! Return: approximation by series expansion for E_1(arg)
   pure real function seriesExpansion(arg)

      implicit none

      real, intent(in) :: arg

      real    :: res, fact  ! result of the summation, 1 / n
      integer :: i          ! counter variable

      real, parameter :: EULER_GAMMA = 0.57721566490153286060651209008241 ! Euler constant
      integer, parameter :: ITERATION = 25 ! Cutoff for series expansion

      ! initialize summation result
      res = 0.0

      ! perform the summation
      do i = ITERATION, 2, -1
         ! calculate 1/n
         fact = 1.0/i
         ! add next term of summation
         res = arg*fact*(fact - res)
      end do

      ! calculate the final result
      seriesExpansion = -EULER_GAMMA - log(arg) + arg*(1.0 - res)

   end function seriesExpansion

   ! The Gauss Laguerre expansion of the exponential integral can be written as
   !
   !              N
   ! E (arg)    -----    a
   !  1          \        n
   ! ------- =    )   --------
   !   -arg      /    x  + arg
   !  e         -----  n
   !             n=1
   !
   ! where the a_n and x_n are determined by least quadrature and are given in [1]
   ! Input: arg - point at which Gaussian Laguerre quadrature is calculated
   ! Return: E_1(arg) in this approximation
   pure real function gauss_laguerre(arg)

      implicit none

      real, intent(in) :: arg

      ! the quadrature constants a_n and x_n from [1]
      real, parameter :: a(1:15) = (/ &
                         0.2182348859400869e+00, 0.3422101779228833e+00, 0.2630275779416801e+00, &
                         0.1264258181059305e+00, 0.4020686492100091e-01, 0.8563877803611838e-02, &
                         0.1212436147214252e-02, 0.1116743923442519e-03, 0.6459926762022901e-05, &
                         0.2226316907096273e-06, 0.4227430384979365e-08, 0.3921897267041089e-10, &
                         0.1456515264073126e-12, 0.1483027051113301e-15, 0.1600594906211133e-19/)
      real, parameter :: x(1:15) = (/ &
                         0.9330781201728180e-01, 0.4926917403018839e+00, 0.1215595412070949e+01, &
                         0.2269949526203743e+01, 0.3667622721751437e+01, 0.5425336627413553e+01, &
                         0.7565916226613068e+01, 0.1012022856801911e+02, 0.1313028248217572e+02, &
                         0.1665440770832996e+02, 0.2077647889944877e+02, 0.2562389422672878e+02, &
                         0.3140751916975394e+02, 0.3853068330648601e+02, 0.4802608557268579e+02/)

      ! Calculate the summation
      gauss_laguerre = sum(a/(x + arg))

   end function gauss_laguerre

end module m_exponential_integral
