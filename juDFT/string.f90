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

   function float2str(num) result(ret_str)
      implicit none
      real, intent(in)               :: num
      character(len=:), allocatable  :: ret_str

      allocate(character(100) :: ret_str)
      if(num >= 1e-1 .and. num <= 1e4) then
         write (ret_str,"(F10.5)") num
      else
         write (ret_str,"(ES12.4)") num
      endif
      ret_str = strip(ret_str)
   end function float2str

   function gen_filename(base, iter, kpt, ext) result(out_str)
      implicit none
      character(len=*), intent(in)           :: base
      integer, intent(in), optional          :: iter, kpt
      character(len=*), intent(in), optional :: ext
      character(len=:), allocatable          :: out_str

      out_str = base
      if(present(iter)) out_str = out_str // "_it=" // int2str(iter)
      if(present(kpt))  out_str = out_str // "_kpt=" // int2str(kpt)

      if(present(ext)) then
         if(ext(1:1) == ".") then
            out_str = out_str // ext
         else
            out_str = out_str // "." // ext
         endif
      endif
   end function gen_filename

end module m_juDFT_string
