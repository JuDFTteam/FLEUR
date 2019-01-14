module m_juDFT_string
   implicit none
   character(len=3), parameter :: whitespaces = " " // achar(9) // achar(13) ! list of all whitespaces
contains
   function strip(input) result(output)
      implicit none
      character(len=*), intent(in)  :: input
      character(:), allocatable     :: output

      integer                       :: front, back

      front = 1
      do while(index(whitespaces, input(front:front)) /= 0 )
         front = front + 1
      enddo

      back  = len(input)
      do while(index(whitespaces, input(back:back)) /= 0)
         back = back - 1
      enddo

      output = input(front:back)
   end function strip
end module m_juDFT_string
