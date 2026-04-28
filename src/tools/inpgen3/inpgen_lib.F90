!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------


   subroutine say_hello(num) BIND(C, name="say_hello")
      USE iso_c_binding
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: num
      INTEGER :: i
      DO i = 1, num
      print *, "Hello from inpgen_lib Fortran module!"
      enddo
   end subroutine say_hello

   SUBROUTINE make_inp_py(simple_input,len_simple_input,profileName,len_profileName,nosym) BIND(C, name="make_inp")
      USE iso_c_binding
      use m_inpgen_make_kpts
      use m_inpgen_make_inp   
      IMPLICIT NONE
      type(c_ptr), INTENT(IN),value :: simple_input
      integer(c_int), INTENT(IN), value :: len_simple_input
      type(c_ptr), INTENT(IN),value :: profileName
      integer(c_int), INTENT(IN), value :: len_profileName
      LOGICAL, INTENT(IN),value :: nosym

      CHARACTER(KIND=c_char,len=len_simple_input),pointer :: simple_input_f
      CHARACTER(KIND=c_char,len=len_profileName),pointer :: profileName_f

   
      call c_f_pointer(simple_input,simple_input_f)
      call c_f_pointer(profileName,profileName_f)


      print *, "In make_inp_py with input of length:",len_trim(simple_input_f)
      print *, simple_input_f
      print *, "Profile name:", profileName_f
      print *, "No symmetry flag:", nosym
      OPEN(UNIT=97, STATUS='scratch')
      WRITE(97, *) simple_input_f
      REWIND(97)

      CALL make_inp_xml(97, "inp.xml", .TRUE., profileName_f, &
                        [.TRUE., .TRUE., .TRUE., .TRUE.], [.TRUE.,.TRUE.,.TRUE.],nosym)

   END SUBROUTINE make_inp_py

   subroutine make_kpt_py(len_kpts_str,kpts_str,len_kpts_path,kpts_path,nosym) BIND(C, name="make_kpt")
      USE iso_c_binding
      use m_inpgen_make_kpts
      IMPLICIT NONE
      INTEGER(c_int), INTENT(IN),value :: len_kpts_str,len_kpts_path
      type(c_ptr), INTENT(IN),value :: kpts_str, kpts_path
      LOGICAL, INTENT(IN),value :: nosym

    
      CHARACTER(KIND=c_char,len=len_kpts_str),pointer :: kpts_str_f
      CHARACTER(KIND=c_char,len=len_kpts_path),pointer :: kpts_path_f

      print *, "In make_kpt_py with lengths:",len_kpts_str,len_kpts_path
  
      call c_f_pointer(kpts_str,kpts_str_f)
      call c_f_pointer(kpts_path,kpts_path_f)

      print *, "In make_kpt_py with kpts_str of length:",len_trim(kpts_str_f)
      print*, "kpts_str:", kpts_str_f
      print *, "In make_kpt_py with kpts_path of length:",len_trim(kpts_path_f)


      CALL add_kpoints(kpts_str_f,kpts_path_f,nosym)

   end subroutine make_kpt_py
