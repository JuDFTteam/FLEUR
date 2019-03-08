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

   function str2int(str) result(int)
      implicit none
      character(len=*), intent(in) :: str
      integer                      :: int
      integer                      :: stat

      read (str, *, iostat=stat) int
      if (stat /= 0) then
         write (*, *) "str reading failed", str
         stop 9
      endif
   end function str2int

   function int2str(num) result(ret_str)
      implicit none
      integer, intent(in)            :: num
      character(len=:), allocatable  :: ret_str

      allocate(character(100) :: ret_str)
      write (ret_str,*) num
      ret_str = strip(ret_str)
   end function int2str
end module m_juDFT_string
