!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_calculator
      use m_juDFT
      !
      !  This module implements a parser able to evaluate expressions in
      !  input files
      !
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: evaluate,ASSIGN_var,delete_vars,evaluatefirst,
     $          makenumberstring,show,evaluateFirstOnly,
     $          evaluateFirstIntOnly,evaluateFirstBoolOnly

      CHARACTER(len = 10),SAVE,ALLOCATABLE :: var_names(:)
      REAL,ALLOCATABLE,SAVE                :: var_values(:)
      INTEGER,SAVE                         :: n_vars
      INTEGER,PARAMETER                    :: num_predef = 6

      CONTAINS


      !<-- S: priv_error(m)
      SUBROUTINE show( )
!-----------------------------------------------
!
!           (last modified: 2012-05-22) pbaum
!-----------------------------------------------
      IMPLICIT NONE

      !<-- Locals
      INTEGER             :: n
      !>

      WRITE(*,*) "Defined variables:"

      DO n = 1,n_vars
         WRITE(*,"(9a)") TRIM(var_names(n)), ' = '
     $     ,makenumberstring(var_values(n))
      ENDDO

      END SUBROUTINE
      !>


      !<-- S: priv_error(m)
      SUBROUTINE priv_error(m)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(in) :: m
      !>
      !<-- Locals
      INTEGER             :: ierr,n
      !>

      WRITE(*,*) m
      WRITE(6,*) m
      WRITE(*,*) "Defined variables:"

      DO n = 1,n_vars
         WRITE(*,"(2a,f0.10)") TRIM(var_names(n)),' = ',var_values(n)
      ENDDO

      CALL juDFT_error("Error in expression",calledby="calculator")
      END SUBROUTINE
      !> 

      !<-- S: priv_increase_storage()
      SUBROUTINE priv_increase_storage()
!-----------------------------------------------
!    increase the storage for the variables if needed
!    This is a very slow procedure, so be careful to
!    adjust the buffer size if it is called very often
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<-- Locals
      CHARACTER(len = 10),ALLOCATABLE :: tmp_names(:)
      REAL   ,ALLOCATABLE             :: tmp_values(:)
      INTEGER,PARAMETER :: min_buffer = 5 
      INTEGER :: i
      !>

      IF (ALLOCATED(var_names)) THEN
         IF (n_vars+1 < SIZE(var_names)) RETURN  !nothing to be done
      ELSE
!        WRITE(*,*) "Predefine variables:", num_predef
         n_vars = num_predef
         ALLOCATE(var_values(n_vars))
         ALLOCATE(var_names(n_vars))
         var_names (1) = 'Pi'
         var_values(1) = 3.1415926535897932384626433832795
         var_names (2) = 'Deg'
         var_values(2) = 17.453292519943295E-3
         var_names (3) = 'Ang'
         var_values(3) =   1.8897261247728981
         var_names (4) = 'nm'
         var_values(4) =   18.897261247728981
         var_names (5) = 'pm'
         var_values(5) = 0.018897261247728981
         var_names (6) = 'Bohr'
         var_values(6) = 1.0
      ENDIF


      !<-- copy old data
      IF (ALLOCATED(var_names)) THEN
         ALLOCATE(tmp_names(SIZE(var_names)))
         tmp_names(:) = var_names
         ALLOCATE(tmp_values(SIZE(var_values)))
         tmp_values(:) = var_values
         DEALLOCATE(var_values,var_names)
      ENDIF
      !>

      ALLOCATE(var_values(n_vars+min_buffer))
      ALLOCATE(var_names(n_vars+min_buffer))

      !<-- Copy data back
      IF (ALLOCATED(tmp_names)) THEN
         var_names(:SIZE(tmp_names)) = tmp_names
         var_values(:SIZE(tmp_values)) = tmp_values
         DEALLOCATE(tmp_names,tmp_values)
      ENDIF
      !>
      END SUBROUTINE
      !> 

      !<-- S: delete_vars()
      SUBROUTINE delete_vars()
!-----------------------------------------------
!    
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      IF (ALLOCATED(var_names)) THEN
         DEALLOCATE(var_names)
         DEALLOCATE(var_values)
         n_vars = 0
      ENDIF
      END SUBROUTINE
      !> 

      !<-- S: assign_var(var,value)
      SUBROUTINE ASSIGN_var(var,value)
!-----------------------------------------------
!  assign a value to a variable
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len = *),INTENT(IN) :: var
      REAL   ,INTENT(IN)            :: value
      !>

      !<-- Locals
      INTEGER             :: n
      CHARACTER(len = 10) :: s
      !>

      s = TRIM(ADJUSTL(var))

      DO n = 1,n_vars
         IF (s == var_names(n)) THEN
            IF (n > num_predef) THEN
               !variable exists, new value assigned
               var_values(n) = value
               RETURN
            ELSE
               CALL priv_error(
     $           'attempt to override predefined constant!')
            END IF
         ENDIF
      ENDDO

      !this is a new variable
      CALL priv_increase_storage()
      n_vars = n_vars+1
      var_names(n_vars) = s
      var_values(n_vars) = value
      END SUBROUTINE
      !> 

      !<-- F: priv_number(string) result(number)
      FUNCTION priv_number(string) result(number)
!-----------------------------------------------
!    read the first part of string as a number
!             (last modified: 06-04-11) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<-- Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      REAL                             :: number
      !>

      !<-- Locals
      INTEGER             :: pos
      LOGICAL             :: dot
      !>

      number = 0.0
      pos = 0
      dot = .false.
      loop:DO
         pos = pos+1
         IF (pos>LEN(string)) THEN
            EXIT
         ENDIF
         SELECT CASE(string(pos:pos))
             CASE ('0':'9') 
                CYCLE loop
             CASE ('+','-')
                IF(pos /= 1) EXIT loop
             CASE ('.') 
                IF(dot) THEN
                   EXIT loop
                ELSE
                   dot = .TRUE.
                ENDIF
             CASE default
                EXIT loop
         END SELECT
      ENDDO loop

      IF (pos == 0) RETURN 
      READ(string(:pos-1),*)number
      IF (pos>LEN(string)) THEN
         string=' '
      ELSE
         string = string(pos:)
      ENDIF
      END function
      !>

      !<-- S: priv_text(string, command, number)
      SUBROUTINE  priv_text(string, command, number)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      CHARACTER(len = 10), INTENT(out) :: command
      REAL, INTENT(out)                :: number
      !>

      !<-- Locals
      INTEGER  :: n,l
      CHARACTER(len = 5),DIMENSION(13), PARAMETER :: functions = (/
     $     'cos( ','sin( ','tan( ','exp( ','log( ','abs( ','sqrt('
     $     ,'acos(','asin(','atan(','cosh(','sinh(','tanh('/)
      !>

      command = ' '
      number = 0.0

      !<--check if this is a function call
      DO n = 1,SIZE(functions)
         l = len_TRIM(functions(n))
         IF (len_TRIM(string) <= l) CYCLE
         IF (string(:l) == functions(n)(:l)) THEN
            command = functions(n)(:l)
            string=string(l:)
            RETURN
         ENDIF
      ENDDO
      !>

      !<-- it must be a variable
      !separate the name of the variable

      l = 1
      DO 
         l = l+1
         IF (l>len_TRIM(string)) EXIT
         SELECT CASE (string(l:l)) 
         CASE ('*','+','-','/',')',' ','^','%')
            EXIT
         END SELECT
      ENDDO

      l = l-1
      DO n = 1,n_vars
         IF (l /= len_TRIM(var_names(n))) CYCLE
         IF (string(:l) == var_names(n)(:l)) THEN
            command = 'variable'
            number  = var_values(n)
            IF (len_TRIM(string)>l) THEN
               string = string(l+1:)
            ELSE
               string=' '
            ENDIF
            RETURN
         ENDIF
      ENDDO
      !>
      CALL priv_error("Unknown character string found: "//TRIM(string
     $     ))
      END SUBROUTINE
      !> 

      !<-- F: priv_operator(string)
      FUNCTION priv_operator(string) result(command)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      CHARACTER                        :: command
      !>

      IF (len_TRIM(string)<2)
     $   CALL priv_error("Parsing error (operator): "//trim(string))

      SELECT CASE( string(1:1) )
         CASE( '+', '-', '/', '%', '^', '*' )
            command = string(1:1)
            !<-- see if the 2nd char is also an operator
            SELECT CASE( string(2:2) )
               CASE( '*' )
                  command = '^' ! convert ** to ^
                  string = string(2:) ! cut away the 1st '*'
               CASE( '+', '-', '/', '%', '^' )
                  CALL priv_error(
     $               'Operator following operator '//TRIM(string))
            END SELECT
         CASE default
            CALL priv_error('Unknown operator: '//TRIM(string))
      END SELECT
      string = string(2:) ! cut away the 1st char

      END FUNCTION
      !>

      !<-- F: priv_bracket(string)
      FUNCTION priv_bracket(string)RESULT(substring)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      CHARACTER(len = LEN(string)) :: substring
      !>

      !<-- Locals
      INTEGER             :: n,pos,cnt
      !>

      pos = 0
      cnt = 0

      loop:DO
         pos = pos+1
         IF (pos>len_TRIM(string)) CALL priv_error("Unbalanced brackets"
     $        )
         SELECT CASE(string(pos:pos))
            CASE('(')
              cnt = cnt+1
            CASE(')')
              cnt = cnt-1
         END SELECT
         IF (cnt == 0) EXIT loop
      ENDDO loop

      substring = TRIM(string(2:pos-1))
      IF (len_TRIM(substring) == 0) CALL priv_error
     $     ("Empty brackets found")
      IF (len_TRIM(string) < pos+1) THEN
         string = ' '
      ELSE
         string = string(pos+1:)
      END IF
      END FUNCTION
      !> 

      !<-- F: priv_calc(string,command,value)
      RECURSIVE FUNCTION priv_calc(string,value,command)
     $     RESULT(number)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len=*), INTENT(inout) :: string
      CHARACTER(len = 10), INTENT(in) :: command
      REAL, INTENT(in)                :: value
      REAL                            :: number
      !>

      !<-- Locals
      CHARACTER(len = 10) :: nextcommand
      REAL                :: nextnumber
      !>

      SELECT CASE(command)
         CASE('number','variable','end')
            number = value
         CASE('bracket')
            number = evaluate(priv_bracket(string))
         CASE('+')
            number = value + priv_evaluateblock(string, command)
         CASE('-')
            number = value - priv_evaluateblock(string, command)
         CASE('*')
            number = value * priv_evaluateblock(string, command)
         CASE('%')
            number = modulo(value,priv_evaluateblock(string, command))
         CASE('/')
            number = priv_evaluateblock(string, command)
            IF (number == 0.0) CALL priv_error("Divide by zero")
            number = value/number
         CASE('^','**')
            number = priv_evaluateblock(string, command)
            IF (number > 0.0) THEN
               number = value ** number
            ELSEIF (value < 0.0) THEN
               IF (INT(number) == number) THEN
                  number = value ** INT(number)
               ELSE
                  CALL priv_error('x^y, x<0 and y not integer')
               END IF
            ELSEIF (number /= 0.0) THEN
               number = 0.0
            ELSE
               CALL priv_error('Undefined result 0^0')
            END IF     
         CASE('cos(','sin(','exp(','log(','abs(','sqrt(','acos(','
     $           asin(','atan(','cosh(','tanh(','tan(')
            call priv_getnextatom(string, nextnumber, nextcommand)
            number = priv_calc(string, nextnumber, nextcommand)
            SELECT CASE (command)
               CASE('sin(')
                  number = SIN(number)
               CASE('cos(')
                  number = COS(number)
               CASE('tan(')
                  number = TAN(number)
               CASE('exp(')
                  number = EXP(number)
               CASE('log(')
                  IF (number <= 0) CALL priv_error("log(x),x <= 0 ")
                  number = LOG(number)
               CASE('abs(')
                  number = ABS(number)
               CASE('sqrt(')
                  IF (number < 0) CALL priv_error("sqrt(x),x <0 ")
                  number = SQRT(number)
               CASE('acos(')
                  IF (ABS(number)>1) CALL priv_error("acos(x), |x|>1")
                  number = ACOS(number)
               CASE('asin(')
                  IF (ABS(number)>1) CALL priv_error("asin(x), |x|>1")
                  number = ASIN(number)
               CASE('atan(')
                  number = ATAN(number)
               CASE('cosh(')
                  number = COSH(number)
               CASE('sinh(')
                  number = SINH(number)
               CASE('tanh(')
                  number = TANH(number)
            END SELECT
         CASE default
            CALL priv_error("Parsing error: "//command)
      END SELECT
      END FUNCTION
      !> 

      !<-- S: priv_getnextatom(func, number, command)
      SUBROUTINE priv_getnextatom(string, number, command)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      REAL, INTENT(inout)              :: number
      CHARACTER(len = 10), INTENT(out) :: command
      !>

      SELECT CASE(string(1:1))
      CASE('0':'9', '.')
         number = priv_number(string)
         command = 'number'
      CASE('+', '-', '/', '*', '%', '^')
         command = priv_OPERATOR(string)
      CASE('a':'z','A':'Z')
         CALL priv_text(string, command, number)
      CASE('(')
         command = 'bracket'
      CASE default
         IF (len_TRIM(string) > 1) CALL
     $        priv_error('Unknowninput:'//TRIM(string))
         command = 'end'
      END SELECT
      END SUBROUTINE
      !>

      !<-- F: priv_order(command) 

      FUNCTION priv_order(command) result(order)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(in) :: command
      INTEGER                       :: order
      !>

      order = 0

      SELECT CASE(TRIM(command))
         CASE('+','-')
            order = 10
         CASE('*','/','%')
            order = 100
         CASE('**','^')
            order = 1000
      END SELECT
      END FUNCTION
      !> 

      !<-- S: priv_peeknextatom(string, number, command)
      SUBROUTINE priv_peeknextatom(string, number, command)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(in) :: string
      REAL, INTENT(inout)           :: number
      CHARACTER(len = 10), INTENT(inout) :: command
      !>

      !<-- Locals
      CHARACTER(len   = LEN(string)) :: s
      !>

      s=string
      CALL priv_getnextatom(s,number,command)
      END SUBROUTINE
      !> 

      !<-- F: priv_evaluateblock(func, blockcommand) result(number)
      RECURSIVE FUNCTION priv_evaluateblock(string, blockcommand)
     $     RESULT(number)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      CHARACTER(len=*), INTENT(in)     :: blockcommand
      REAL                             :: number
      !>

      !<-- Locals
      CHARACTER(len = 10) :: command, nextcommand
      REAL                :: nextnumber
      !>

      CALL priv_peeknextatom(string, nextnumber, nextcommand)

      IF (TRIM(nextcommand) == 'end') CALL priv_error
     $     ('Premature end of function')

      DO WHILE((priv_order(nextcommand) == 0 .OR.
     $     (priv_order(nextcommand) >priv_order(blockcommand))) .AND.
     $     TRIM(nextcommand) /= 'end')
         CALL priv_getnextatom(string, number, command)
         number = priv_calc(string, number, command)
         CALL priv_peeknextatom(string, nextnumber, nextcommand)
      END DO
      END FUNCTION
      !> 

      !<-- F: evaluate(s) 
      RECURSIVE FUNCTION evaluate(s) RESULT(number)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(in  ) :: s
      REAL                            :: number
      !>

      !<-- Locals
      CHARACTER(len = 10) :: command
      CHARACTER(len = len_trim(s)) :: tmp_s
      !>

      tmp_s  = TRIM(ADJUSTL(s))
      number = 0
      command = ' '

      DO WHILE(command /= 'end')
         CALL priv_getnextatom(tmp_s, number, command)
         number = priv_calc(tmp_s, number, command)
      END DO
      END FUNCTION
      !> 

      !<-- F: makenumberstring(x)
      FUNCTION makenumberstring(x) result( str )
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!             (last modified: 2012-05-22) pbaum
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      REAL   ,INTENT(IN)     :: x
      CHARACTER(len=20)      :: str
      !>

      !<-- Locals
      INTEGER                :: n, i, m
      REAL                   :: xx, xm
      CHARACTER(len =10)     :: fn ! function string
      CHARACTER              :: p2 ! closing parenthesis
      CHARACTER              :: sg ! sign
      !>

      DO i = 0,2
         !<-- assume that xx = fn(x)
         !<-- apply the inverse function x = fn^{-1}(xx)
         SELECT CASE (i)
         CASE ( 2 )   ; fn = 'tan('   ; p2 = ')' ; xx = atan(abs(x))
         CASE ( 1 )   ; fn = 'sqrt('  ; p2 = ')' ; xx = x*x
         CASE DEFAULT ; fn = ''       ; p2 = ' ' ; xx = abs(x)
         END SELECT

         !<-- restore the sign to appear in front of the function
         sg = ' ' ; IF (x < 0.) sg = '-'

         !<-- check if xx is a simple fraction
         n = 0 ! 0: not found
         m = 0 ! init
         DO WHILE (m < 25)
           m = m+1
           xm = xx*m
           IF (ABS(xm-NINT(xm))<1E-6) THEN
             n = m
             EXIT
           END IF
         END DO

         IF (n == 1) THEN
            WRITE(str,"(2a,i0,a)")
     $         TRIM(sg), TRIM(fn), NINT(xx), TRIM(p2)
         ELSE IF (n > 1) THEN
            WRITE(str,"(2a,i0,a,i0,a)")
     $         TRIM(sg), TRIM(fn), NINT(xx*n), '/', n, TRIM(p2)
         ENDIF
         IF (n > 0) RETURN
      END DO
      !>

      !ok nothing found
      WRITE(str,"(f20.12)") x
      i = LEN(str)
      DO WHILE( str(i:i) == '0' )
         str(i:i) = ' '
         i = i-1
      END DO
      str = ADJUSTL(str)
      END FUNCTION
      !> 

      !<-- F: evaluateFirst(string)
      FUNCTION evaluateFirst(s,n)result(number)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE

      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: s
      INTEGER,OPTIONAL                 :: n
      REAL                             :: number
      !>

      !<-- Locals
      INTEGER             :: pos
      !>

      s = ADJUSTL(s)
      IF (len_TRIM(s) == 0) THEN
         number = 0
         RETURN
      ENDIF
      pos = INDEX(s," ")
      IF (pos == 0) pos = LEN(s)
      IF (PRESENT(n)) pos = MAX(pos,n)
      number  = evaluate(s(:pos))
      IF (pos<LEN_TRIM(s)) THEN
         s = s(pos:)
      ELSE
         s =" "
      ENDIF
      END FUNCTION
      !> 


      FUNCTION evaluateFirstOnly(s)result(number)

      IMPLICIT NONE

      CHARACTER(len =*), INTENT(IN)    :: s
      REAL                             :: number

      !<-- Locals
      CHARACTER(len=LEN(s)) :: tempS
      INTEGER               :: pos
      !>

      tempS = ADJUSTL(s)
      IF (len_TRIM(tempS) == 0) THEN
         number = 0
         RETURN
      ENDIF
      pos = INDEX(tempS," ")
      IF (pos == 0) pos = LEN(tempS)
      number  = evaluate(tempS(:pos))

      END FUNCTION


      FUNCTION evaluateFirstIntOnly(s)result(number)

      IMPLICIT NONE

      CHARACTER(len =*), INTENT(in) :: s
      INTEGER                          :: number

      !<-- Locals
      INTEGER               :: pos
      CHARACTER(len=LEN(s)) :: tempS
      !>

      tempS = ADJUSTL(s)
      IF (len_TRIM(tempS) == 0) THEN
         number = 0
         RETURN
      END IF
      pos = INDEX(tempS," ")
      IF (pos == 0) pos = LEN(tempS)
      number = NINT(evaluate(tempS(:pos)))

      END FUNCTION

      FUNCTION evaluateFirstBoolOnly(s)result(bool)

      IMPLICIT NONE

      CHARACTER(len =*), INTENT(in)    :: s
      LOGICAL                          :: bool

      !<-- Locals
      INTEGER               :: pos
      CHARACTER(len=LEN(s)) :: tempS
      !>

      tempS = ADJUSTL(s)
      IF (len_TRIM(tempS) == 0) THEN
         CALL juDFT_error("String is empty.",
     +                    calledby ="calculator")
         RETURN
      END IF

      pos = INDEX(tempS," ")
      IF (pos == 0) pos = LEN(tempS)
      SELECT CASE (tempS(:pos))
         CASE ('F','f','false','FALSE')
            bool = .FALSE.
         CASE ('T','t','true','TRUE')
            bool = .TRUE.
         CASE DEFAULT
            CALL juDFT_error('No valid bool at start of: ' // tempS,
     +                       calledby ="calculator")
      END SELECT

      END FUNCTION

      END MODULE m_calculator
