module m_juDFT_string
   implicit none
   character(len=3), parameter :: whitespaces = " "//achar(9)//achar(13) ! list of all whitespaces
   interface int2str
      module procedure int2str_int4, int2str_int8
   end interface int2str
contains
   function strip(input) result(output)
      implicit none
      character(len=*), intent(in)  :: input
      character(:), allocatable     :: output

      integer                       :: front, back

      front = 1
      do while (index(whitespaces, input(front:front)) /= 0)
         front = front + 1
         if (front>len(input)) THEN
            front=front-1
            EXIT
         endif
      enddo

      back = len(input)
      do while (index(whitespaces, input(back:back)) /= 0)
         if (back==1) exit
         back = back - 1
      enddo

      output = input(front:back)
   end function strip

   function str2int(str) result(int)
      implicit none
      character(len=*), intent(in) :: str
      integer                      :: int
      integer                      :: stat

      read(str, *, iostat=stat) int
      if(stat /= 0) then
         write(*, *) "str reading failed", str
         stop 9
      endif
   end function str2int

   function int2str_int4(num, length) result(ret_str)
      implicit none
      integer, intent(in)            :: num
      integer, intent(in), optional  :: length
      character(len=:), allocatable  :: ret_str

      allocate (character(100) :: ret_str)
      write (ret_str, *) num
      ret_str = strip(ret_str)

      if(present(length)) then 
         do while (len(ret_str) < length) 
            ret_str = " " // ret_str 
         enddo 
      endif
   end function int2str_int4

   function int2str_int8(num, length) result(ret_str)
      implicit none
      integer(8), intent(in)         :: num
      integer, intent(in), optional  :: length
      character(len=:), allocatable  :: ret_str

      allocate (character(100) :: ret_str)
      write (ret_str, *) num
      ret_str = strip(ret_str)

      if(present(length)) then 
         do while (len(ret_str) < length) 
            ret_str = " " // ret_str 
         enddo 
      endif
   end function int2str_int8

   function float2str(num) result(ret_str)
      implicit none
      real, intent(in)               :: num
      character(len=:), allocatable  :: ret_str

      allocate (character(100) :: ret_str)
      if (num >= 1e-1 .and. num <= 1e4) then
         write (ret_str, "(F10.5)") num
      else
         write (ret_str, "(ES12.4)") num
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
      if (present(iter)) out_str = out_str//"_it="//int2str(iter)
      if (present(kpt)) out_str = out_str//"_kpt="//int2str(kpt)

      if (present(ext)) then
         if (ext(1:1) == ".") then
            out_str = out_str//ext
         else
            out_str = out_str//"."//ext
         endif
      endif
   end function gen_filename

   function replace_text(s, text, rep) RESULT(outs)
      CHARACTER(*)           :: s, text, rep
      CHARACTER(LEN(s) + 100) :: outs     ! provide outs with extra 100 char len
      INTEGER                 :: i, nt, nr

      outs = s; nt = LEN_TRIM(text); nr = LEN_TRIM(rep)
      DO
         i = INDEX(outs, text(:nt)); IF (i == 0) EXIT
         outs = outs(:i - 1)//rep(:nr)//outs(i + nt:)
      END DO
   end function Replace_Text

   function get_byte_str(byte) result(str)
      implicit none
      real, intent(in)              :: byte
      character(len=:), allocatable :: str
      integer :: l10

      l10 = floor(log10(byte))

      if (l10 >= 12) then
         str = float2str(1e-12*byte) //" TB"
      elseif (l10 >= 9) then
         str = float2str(1e-9*byte) //" GB"
      elseif (l10 >= 6) then
         str = float2str(1e-6*byte) //" MB"
      elseif (l10 >= 3) then
         str = float2str(1e-3*byte) //" kB"
      else
         str = float2str(byte) //" byte"
      endif
   end function get_byte_str
end module m_juDFT_string
