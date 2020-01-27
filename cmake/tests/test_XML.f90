      program test
      integer:: i
      interface
         function xml()bind(C,name="xmlInitParser")
              use iso_c_binding
              INTEGER(c_int) xmlInitParser
         end function xml
      end interface
  
      i=xml()
      end program test
