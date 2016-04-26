      MODULE m_relax
      use m_juDFT
      USE m_calculator
!-----------------------------------------------
! DESC: This module contains subroutines for the geometry 
!       optimization.
!       subroutines:
!       relax:              must be called from the force part of FLEUR
!                           (only public subroutine)  
!       read_relax_file:    reads the relax_inp file (see documentation
!                           later)
!       calculate_*_shift:  different strategies for calculating shift
!   
!                 Daniel Wortmann, (06-04-27)
!-----------------------------------------------
      PRIVATE

      TYPE t_force
         CHARACTER(len = 20) :: formula,coord,default
         REAL                :: pos,force,shift
      END TYPE

      INTEGER,SAVE              :: mode
      INTEGER,SAVE              :: n_force
      TYPE(t_force)             :: f(100)
      REAL,PARAMETER            :: max_shift = 0.25
      REAL                      :: min_shift = 0.01

      PUBLIC relax

      CONTAINS
 
      !<-- S: read_relax_file()
#ifndef __PGI
      SUBROUTINE READ_relax_file()
!-----------------------------------------------
!     All behavior of the relaxation is controlled by the relax_inp file
!     Example:
!     broyden min_shift=0.01
!     atom 2: z=10.0
!     coord d12=10.0 : fz1-fz2
!
!     In the first line the broyden/simple scheme can be choosen
!     min_shift gives the convergence target for the shifts
!     Two ways of specifying coordinates are possible:
!     1: give the atom keyword+ no of atom type
!        you choose the x,y,z coordinate
!        The factor for the simple relaxation must be specified
!     2: a name of a variable defined by def in the beginning of the
!        inp file can be specified by the coord keyword
!        Again the propotionality factor should be given
!        After the colon the formula to calculate the force for this coordinate
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      USE m_calculator
      IMPLICIT NONE

      !<--Locals
      CHARACTER(len = 100) :: line,s,varname,value
      INTEGER              :: pos,atomnumber,err
      !>

      OPEN(99,file ="relax_inp",status ="old",iostat= err)
      IF (err /= 0) CALL juDFT_error("No relax_inp found")
      !<-- read header of relax file

      READ(99,"(a)") line
      SELECT CASE(line(1:1))
         CASE ("s","S")
            mode = 0
         CASE ("g","G")
            mode = 1
         CASE("b","B")
            mode = 2
         CASE default
            mode = 0
      END SELECT
      IF (INDEX(line,"min_shift =")>0) THEN
         min_shift             = evaluatefirst(line(index(line
     $        ,"min_shift =")+11:))
      ENDIF

      !>
      DO
         READ(99,"(a)",END = 100) line
         line = TRIM(ADJUSTL(line))
         SELECT CASE (line(1:1))
            CASE ("a")
               line = TRIM(ADJUSTL(line(5:)))
               pos = INDEX(line,":")
               IF (pos < 2.OR.pos >= len_TRIM(line)) CALL
     +             juDFT_error("Error in reading atom line in relax_inp"
     +             ,calledby ="relax")
               atomnumber = evaluate(line(:pos-1))
               line = line(pos+1:)
               DO WHILE (len_TRIM(line)>1)
                  line = TRIM(ADJUSTL(line))
                  IF (INDEX(line,",")>0) THEN
                     s    = line(:INDEX(line,",")-1)
                     line = line(INDEX(line,",")+1:)
                  ELSE
                     s    = line
                     line =" "
                  ENDIF
                  IF (INDEX(s,"=")>1) THEN
                     varname = s(:INDEX(s,"=")-1)
                     value   = s(INDEX(s,"=")+1:)
                  ELSE
                     varname = s
                     value   = "0.0"
                  ENDIF
                  CALL add_force_atom(atomnumber,ADJUSTL(varname),value)
               ENDDO
            CASE("c")
               line = TRIM(ADJUSTL(line(6:)))
               pos = INDEX(line,":")
               IF (pos < 2.OR.pos >= len_TRIM(line)) CALL
     +            juDFT_error("Error in reading coord line in relax_inp"
     +            ,calledby ="relax")
               s = TRIM(ADJUSTL(line(:pos-1)))
               line = line(pos+1:)
               IF (INDEX(s,"=") /= 0) THEN
                  varname   = s(:INDEX(s,"=")-1)
                  value     = s(INDEX(s,"=")+1:)
               ELSE
                  varname = s
                  value   = "0.0"
               ENDIF
               CALL add_force_coord(varname,value,line)
         END SELECT
      ENDDO
 100  CLOSE(99)

      END SUBROUTINE

      !> 

      !<-- S: add_force_atom(no,var,value)

      SUBROUTINE add_force_atom(no,var,value)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)           :: no
      CHARACTER(len =*),INTENT(in) :: var,value
      !>

      n_force = n_force+1
      IF (n_force>SIZE(f)) CALL juDFT_error("Too many forces",calledby
     +     ="relax")
      WRITE(f(n_force)%coord,"(a1,a1,i0)") "#",var(1:1),no
      WRITE(f(n_force)%formula,"(a1,a1,i0)") "f",var(1:1),no
      f(n_force)%default = value
      
      END SUBROUTINE

      !> 
      !<-- S: add_force_coord(no,var,value)

      SUBROUTINE add_force_coord(varname,value,line)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len=*),INTENT(in) ::varname,value,line
      !>

      n_force = n_force+1

      IF (n_force>SIZE(f)) CALL juDFT_error("Too many forces",calledby
     +     ="relax")
      f(n_force)%coord = varname
      f(n_force)%formula = line
      f(n_force)%default = value
      
      END SUBROUTINE

      !> 

      !<-- S: set_positions(pos)
      SUBROUTINE set_positions(pos,neq)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      REAL   ,INTENT(IN)     :: pos(:,:)
      INTEGER,INTENT(IN)     :: neq(:)
      !>
      !<-- Locals
      INTEGER :: n,i
      !>

      DO n = 1,n_force 
         IF (f(n)%coord(1:1) =="#") THEN
            !<-- set all positions that are atomic coordinates
            READ(f(n)%coord(3:),"(i4)") i
            i = SUM(neq(:i-1))+1
            SELECT CASE(f(n)%coord(2:2))
               CASE("x","X")
                  f(n)%pos = pos(1,i)
               CASE("y","Y")
                  f(n)%pos = pos(3,i)
               CASE("z","Z")
                  f(n)%pos = pos(3,i)
            END SELECT
            !>
         ELSE
            !<-- the coordinate names should be known to calculator
            f(n)%pos = evaluate(f(n)%coord)
            !>
         ENDIF
      ENDDO
      

      END SUBROUTINE
      !> 
      !<-- S: set_forces(f_in)
      SUBROUTINE set_forces(f_in)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      REAL   ,INTENT(IN)     :: f_in(:,:)
      !>
      !<-- Locals
      CHARACTER(len = 6) :: s
      INTEGER            :: n
      !>
      
      !<-- put all forces into the calculator

      DO n = 1,SIZE(f_in,2)
         WRITE(s,"(a1,a1,i0)") "f","x",n
         CALL ASSIGN_var(s,f_in(1,n))
         WRITE(s,"(a1,a1,i0)") "f","y",n
         CALL ASSIGN_var(s,f_in(2,n))
         WRITE(s,"(a1,a1,i0)") "f","z",n
         CALL ASSIGN_var(s,f_in(3,n))
      ENDDO

      !>

      !<--calculate all forces
      DO n=1,n_force
         f(n)%force = evaluate(f(n)%formula)
      ENDDO
      !>


      END SUBROUTINE
      !> 

      !<-- S: calculate_simple_shift()

      SUBROUTINE calculate_simple_shift()
!-----------------------------------------------
!   calculates the shift by simply multiplying the default with the force
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<-- Locals
      INTEGER   :: n
      !>

      DO n = 1,n_force
         f(n)%shift = evaluate(f(n)%default)
         IF (f(n)%shift == 0.0) f(n)%shift  = 1.0
         f(n)%shift = f(n)%force*f(n)%shift
         f(n)%shift = MIN(max_shift,f(n)%shift) !limit shift
         f(n)%shift = MAX(-max_shift,f(n)%shift) !limit shift
      ENDDO
      CALL WRITE_history(.TRUE.)
      END SUBROUTINE

      !> 
      !<-- S: calculate_grad_shift()
      SUBROUTINE calculate_grad_shift()
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments


      !>
      !<--Locals
      INTEGER             :: n,n_old
      REAL                :: f1,f2
      REAL,POINTER        :: p_old(:,:),s_old(:,:),f_old(:,:)
      !>

      !<-- first read the old positions

      IF(READ_history(p_old,f_old,s_old) == 0) THEN
         WRITE(*,*) "No forces.hist found, perform simple relaxation"
         CALL calculate_simple_shift()
         RETURN
      ENDIF

      !>

      !<-- Here we use only the last two entries
      n_old = SIZE(p_old,2)-1
      DO n = 1,n_force
         f1 = (f_old(n,n_old+1)-f_old(n,n_old))/(p_old(n,n_old+1)
     $        -p_old(n,n_old))
         f2 = f_old(n,n_old+1)-f1*p_old(n,n_old+1)
         f(n)%shift = -1*f2/f1-f(n)%pos
      ENDDO
      !>
      DEALLOCATE(p_old,f_old,s_old)
      CALL WRITE_history(.FALSE.)
      END SUBROUTINE
      !> 
      !<-- S: calculate_broyden_shift()
      SUBROUTINE calculate_broyden_shift()
!-----------------------------------------------
!  Use the broyden method to calculate the shifts
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments


      !>
      !<--Locals
      REAL   ,POINTER :: p_old(:,:),f_old(:,:),s_old(:,:)
      INTEGER         :: n,i,j,hist_length
      REAL            :: h(n_force,n_force)
      REAL            :: p(n_force),y(n_force),v(n_force)
      REAL            :: py,yy,gamma
      LOGICAL         :: initial

      !>
      
      !<-- first read the old positions

      IF(READ_history(p_old,f_old,s_old) == 0) THEN
         WRITE(*,*) "No forces.hist found, perform simple relaxation"
         CALL calculate_simple_shift()
         RETURN
      ENDIF

      !>
      
      !<--calculate approx. Hessian
      !<--initialize H
      h = 0.0
      DO n = 1,n_force
         h(n,n) = 1.0
      ENDDO
      !>
      !<--loop over all iterations (including current)
      hist_length=size(s_old,1)
      WRITE(*,*) "Calculating Hessian with history of:",hist_length-1
      DO n = 2,hist_length
         ! differences
         p(:) = p_old(:,n)-p_old(:,n-1)
         y(:) = f_old(:,n-1)-f_old(:,n)
         ! get necessary inner products and H|y>
         py = dot_PRODUCT(p,y)
         v = MATMUL(y,h)
         yy = dot_PRODUCT(y,v)
         !check that update will leave h positive definite;
         IF (py <= 0.0) THEN
            WRITE (6,*) '  bfgs: <p|y> < 0'
            WRITE (6,*) '  check convergence of forces'
            !<-- Starting over with initial hessian
            h = 0.0
            DO j = 1,n_force
               h(j,j)     = 1.0
            ENDDO
            initial = .TRUE.
            CYCLE 
            !>
         ELSE
            !<-- update h
            IF (n == 2) THEN
               gamma = py/yy
            ELSE
               gamma = 1.0
            ENDIF
            DO j = 1,n_force
               DO i = 1,n_force
                  h(i,j) = (h(i,j) - (v(i)*p(j)+p(i)*v(j))/py)*gamma
     +                 + (1.+gamma*yy/py)*p(i)*p(j)/py
               ENDDO
            ENDDO
            !>
            initial = .FALSE.
         ENDIF
      ENDDO 
      !>
      !>
      !<-- new shift
      DEALLOCATE(f_old,s_old,p_old)
      CALL WRITE_history(.FALSE.)
      IF (initial) THEN
         call calculate_simple_shift()
      ELSE
         WRITE(*,*) "Hessian:",h,"End H"
         f(:n_force)%shift = MATMUL(f(:n_force)%force,h)
      ENDIF
      !>
      END SUBROUTINE
      !> 

      !<-- S: update_positions(pos,mrot,tau,amat,ngopr)
      SUBROUTINE update_positions(pos,neq,mrot,tau,amat,bmat,ngopr
     $     ,invtab)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      use m_constants
      IMPLICIT NONE
      !<--Arguments
      REAL   ,INTENT(INOUT)     :: pos(:,:)
      INTEGER,INTENT(IN)        :: mrot(:,:,:),neq(:)
      REAL   ,INTENT(IN)        :: tau(:,:)
      INTEGER,INTENT(IN)        :: invtab(:),ngopr(:)
      REAL   ,INTENT(IN)        :: bmat(3,3),amat(3,3)
      !>
      !<-- Locals
      INTEGER             :: n,i,na,jop
      REAL                :: p(3)
      !>
      DO n = 1,n_force 
         IF (f(n)%coord(1:1) =="#") THEN
            !<-- set all positions that are atomic coordinates
            READ(f(n)%coord(3:),"(i4)") i
            write(*,*) i
            IF (i>1) THEN 
               i = SUM(neq(:i-1))+1
            ENDIF
            write(*,*) i
            SELECT CASE(f(n)%coord(2:2))
               CASE("x","X")
                  pos(1,i) = f(n)%pos+f(n)%shift
               CASE("y","Y")
                  pos(2,i) = f(n)%pos+f(n)%shift
               CASE("z","Z")
                  pos(3,i) = f(n)%pos+f(n)%shift
            END SELECT
            !>
         ENDIF
      ENDDO
      
      write(*,*) "Calculated positions"
      DO n=1,SIZE(pos,2)
         WRITE(*,"(i5,3(1x,f0.8))") n,pos(:,n)
      ENDDO
      
      !<-- Now recalculate the equivalent positions

      na = 0
      DO n = 1,size(neq)
         na = na+1
         p = MATMUL(pos(:,na),bmat)/2/pimach()
         DO i = 2,neq(n)
            na=na+1
            jop = invtab(ngopr(na))
            WRITE(*,*) "p",p
            write(*,*) "jop",jop
            write(*,*) "mrot",mrot(:,:,jop)
            write(*,*) "tau",tau(:,jop)
            pos(:,na)=matmul(amat,matmul(mrot(:,:,jop),p)+tau(:,jop))
         ENDDO
      ENDDO

      write(*,*) matmul(amat,bmat)

      !>

      WRITE(*,*) "Equivalent positions"
      DO n=1,SIZE(pos,2)
         WRITE(*,"(i5,3(1x,f0.8))") n,pos(:,n)  
      ENDDO

      END SUBROUTINE
      !> 

      !<-- S: write_inp_new(pos,neq)
      SUBROUTINE write_inp_new(pos,neq,bmat,film)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      USE m_constants
      IMPLICIT NONE
      !<--Arguments
      REAL   ,INTENT(IN)     :: pos(:,:)
      INTEGER,INTENT(IN)     :: neq(:)
      REAL   ,INTENT(IN)     :: bmat(3,3)
      LOGICAL,INTENT(IN)     :: film
      !>
      !<-- Locals
      REAL                :: p(3)
      CHARACTER(len = 100) :: line
      INTEGER             :: n,na,i,ii,iii
      !>

      open(98,file="inp",status="old")
      open(99,file="inp_new",status="new")

      !<--Copy all up to first line starting with a "*"
      line=" "
      DO WHILE(line(1:1) /="*")
         READ(98,"(a)") line
         !<--if the line is a "def coord", replace the coordinate
         IF (line(1:4) =="def ") THEN
            DO n  = 1,n_force
               IF (f(n)%coord(1:1) =="#") CYCLE
               na = len_TRIM(f(n)%coord)
               IF (len_TRIM(line)-4 < na) CYCLE
               WRITE(*,*) "Testing:",line(5:na+5),f(n)%coord(:na)
               IF (line(5:na+4) == f(n)%coord(:na)) THEN
                  write(*,*) "MATCH"
                  line = line(:na+4)//"="
     $                 //makenumberstring(f(n)%pos+f(n)%shift)
               ENDIF
            ENDDO
            !>
         ENDIF
         WRITE(99,"(a)") line
      ENDDO
      !>
      !<-- write atomic positions
      na = 0
      DO n = 1,SIZE(neq)
         READ(98,"(a)") line; WRITE(99,"(a)") line
         READ(98,"(a)") line; WRITE(99,"(a)") line
         READ(98,"(a)") line; WRITE(99,"(a)") line
         atomloop:DO i = 1,neq(n)
            na = na+1
            READ(98,"(a)") line
            !<--Check if coordinates of this atom are modified
            DO ii = 1,n_force
               IF (f(ii)%coord(1:1) =="#") THEN
                  READ(f(n)%coord(3:),"(i4)") iii
                  IF (iii == n) THEN
                     !<--Forces found on some atomic coordinate
                     p = MATMUL(pos(:,na),bmat)/2/pimach()
                     IF (film) THEN
                        p(3) = pos(3,na)
                     ENDIF
                     WRITE(99,"(5a)") makenumberstring(p(1))," ",
     $                    makenumberstring(p(2))," ",
     $                    makenumberstring(p(3))
                     CYCLE atomloop
                     !>
                  ENDIF
               ENDIF
            ENDDO
            !>
            WRITE(99,"(a)") line
         ENDDO atomloop
      ENDDO
      !>
      !<-- copy rest of inp file
      DO 
         READ(98,"(a)",END = 110) line
         WRITE(99,"(a)") line
      ENDDO
      !>
 110  CLOSE(98)
      CLOSE(99)

      END SUBROUTINE
      !> 
#endif
      !<-- S: relax(pos,neq,mrot,tau,amat,bmat,ngopr,invtab,f_in)

      SUBROUTINE relax(film,pos_in,neq,mrot,tau,amat,bmat,ngopr,invtab
     $     ,f_in)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      REAL   ,INTENT(IN)        :: pos_in(:,:)
      INTEGER,INTENT(IN)        :: mrot(:,:,:),neq(:)
      REAL   ,INTENT(IN)        :: tau(:,:)
      INTEGER,INTENT(IN)        :: invtab(:),ngopr(:)
      REAL   ,INTENT(IN)        :: bmat(3,3),amat(3,3)
      REAL   ,INTENT(IN)        :: f_in(:,:)
      LOGICAL,INTENT(IN)        :: film

      !>
      !<-- Locals
      REAL                :: pos(SIZE(pos_in,1),SIZE(pos_in,2))

      !>
#ifndef __PGI
      pos = pos_in

      CALL READ_relax_file()
      CALL set_positions(pos,neq)
      CALL set_forces(f_in)
      
      SELECT CASE(mode)
         CASE (0)
            CALL calculate_simple_shift()
         CASE(1)
            CALL calculate_grad_shift()
         CASE(2)
            CALL calculate_broyden_shift()
      END SELECT

      CALL update_positions(pos,neq,mrot,tau,amat,bmat,ngopr
     $     ,invtab)

!      IF (ALL(ABS(f%shift)<min_shift)) THEN
!         !traditional stop
!         WRITE (6,'(a)') "Des woars!"
!          CALL juDFT_error(" RELAX Des woars ",calledby="relax")
!      ENDIF

      CALL WRITE_inp_new(pos,neq,bmat,film)

      !ok all done
      n_force = 0
#endif
      END SUBROUTINE

      !> 

      !<-- S: write_history(new)
      SUBROUTINE WRITE_history(new)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      LOGICAL,INTENT(IN)     :: new
      !>
      !<-- Locals
      INTEGER             :: n

      !>
      IF (new) THEN
         OPEN(99,file ="forces.hist",status="replace")
      ELSE
         OPEN(99,file ="forces.hist",position="append")
      ENDIF

      DO n = 1,n_force
         WRITE(99,"(a10,3f15.10)") f(n)%coord,f(n)%pos,f(n)%force,f(n
     $        )%shift
      ENDDO
      CLOSE(99)


      END SUBROUTINE
      !> 

      !<-- F: read_history(p_old,f_old,s_old)
      INTEGER FUNCTION READ_history(p_old,f_old,s_old)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      REAL,POINTER ::p_old(:,:),f_old(:,:),s_old(:,:)
      !>
      !<-- Locals
      INTEGER             :: n,n_old,err,nn
      !>
      READ_history = 0
      OPEN(99,file ="forces.hist",status ="old",iostat= err)
      IF (err /= 0) RETURN

      !<--Check how many entries can be read
      n_old=0
      DO
         DO n = 1,n_force
            READ(99,*,END = 120)
         ENDDO
         n_old = n_old+1
      ENDDO
 120  REWIND(99)
      !>
      ALLOCATE(p_old(n_force,n_old+1))
      ALLOCATE(f_old(n_force,n_old+1))
      ALLOCATE(s_old(n_force,n_old+1))
      !<-- read file
      DO nn=1,n_old
         DO n = 1,n_force
            READ(99,"(10x,3f15.10)") p_old(n,nn),f_old(n,nn),s_old(n,nn)
         ENDDO
      ENDDO
      !>
      WRITE(*,"(a,i0,a)") "History file contains ",n_old
     $     ," previous force step(s)"
      !<-- add current force and position as last entry
      p_old(:,n_old+1) = f(:n_force)%pos
      f_old(:,n_old+1) = f(:n_force)%force
      !>
      CLOSE(99)
      READ_history = n_old
      END FUNCTION
      !> 
      END
