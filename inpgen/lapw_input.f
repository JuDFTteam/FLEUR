      MODULE m_lapwinput      
      use m_juDFT
!-----------------------------------------------------------------------
! read in some lapw-specific input or set appropriate defauts
!-----------------------------------------------------------------------
      INTEGER, PARAMETER  :: dbgfh=6, errfh=6, warnfh=6

      CONTAINS
      SUBROUTINE lapw_input(
     >                      infh,nline,xl_buffer,bfh,buffer,
     <                      jspins,kcrel,ndvgrd,nkpt,div,
     <                      frcor,ctail,chng,tria,kmax,gmax,gmaxxc,
     <                      dvac,dtild,tkb,namex,relcor)
     
      USE m_readrecord
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: xl_buffer,infh,bfh
      INTEGER, INTENT (OUT) :: jspins,kcrel,ndvgrd,nkpt,div(3)
      LOGICAL, INTENT (OUT) :: frcor,ctail,tria
      REAL,    INTENT (OUT) :: kmax,gmax,gmaxxc,tkb,chng
      REAL,    INTENT (INOUT) :: dvac,dtild
      CHARACTER(len=4), INTENT (OUT) :: namex
      CHARACTER(len=12),INTENT (OUT) :: relcor
      CHARACTER(len=xl_buffer)       :: buffer
      
      INTEGER iflag,div1,div2,div3,nline,nbuffer,ios
      LOGICAL h_film,h_comp,h_exco,h_kpt,fatalerror,relxc
      CHARACTER(len=4) :: xctyp

!---> namelists for input
      NAMELIST /comp/   jspins, frcor, ctail, kcrel, 
     &                  gmax, gmaxxc, kmax
      NAMELIST /exco/   xctyp, relxc 
      NAMELIST /film/   dvac, dtild
      NAMELIST /kpt/    nkpt, div1, div2, div3, tkb, tria


      h_film=.false. ; h_comp=.false.
      h_exco=.false. ; h_kpt=.false.
      fatalerror=.false.
                          ! jspins, gmax, gmaxxc, kmax were set before
        frcor  = .false.  ! no frozen core
        ctail  = .true.   ! always core-tail correction
        kcrel = 0         ! no fully-magnetic dirac core

        relcor = 'non-relativi'
        namex = 'pbe '
        ndvgrd = 6 ; chng= -1.0e-12 
        nkpt = 0 ; div = 0 
        tkb = 0.001 ; tria = .false.
 

!===> read input

      nbuffer = len_trim(buffer)

      loop: DO

      IF (nbuffer == 0) then
        DO
          CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
          IF (ios==1) GOTO 999
          IF (ios == 2)  CALL juDFT_error
     +         ("end of file while reading a record",calledby
     +         ="lapw_input")
          IF (buffer(1:1)=='&') EXIT
          CALL err(0)
          fatalerror = .true.
        ENDDO
      ENDIF

!===> comp
      IF (buffer(1:5)=='&comp') THEN
        IF (h_comp) CALL err(1)
         
        READ (bfh,comp,err=912, end=912, iostat=ios)
        h_comp = .true.
        IF (jspins>2 .OR. jspins<1)  CALL juDFT_error
     +       ("jspins>2 .OR. jspins<",calledby ="lapw_input")
        IF (kcrel <0 .or. kcrel >1)  CALL juDFT_error
     +       ("kcrel <0 .or. kcrel >1",calledby ="lapw_input")

!===> exco
      ELSEIF (buffer(1:5)=='&exco') THEN
        IF (h_exco) CALL err(1)
        
        READ (bfh,exco,err=912, end=912, iostat=ios)
        h_exco = .true.

        iflag=-9999
       IF ((xctyp == 'l91 ').OR.(xctyp == 'xa  ').OR.(xctyp == 'wign').
     &  OR.(xctyp == 'hl  ').OR.(xctyp == 'bh  ').OR.(xctyp == 'mjw ').
     &  OR.(xctyp == 'vwn ').OR.(xctyp == 'pz  ') ) iflag=1

       IF ((xctyp == 'pw91').OR.(xctyp == 'pbe ').OR.(xctyp == 'rpbe').
     &  OR.(xctyp == 'Rpbe').OR.(xctyp == 'wc  ') ) iflag=2
        
       IF (iflag < -1 )  CALL juDFT_error("error reading lda/gga",
     +      calledby="lapw_input")

        namex = xctyp
        IF (relxc) relcor = 'relativistic'

        IF (iflag==2) THEN
        ELSE
          ndvgrd = 0 ; chng= 0.0
        ENDIF

!===> film
      ELSEIF (buffer(1:5)=='&film') THEN
        IF (h_film) CALL err(1)

        READ (bfh,film,err=912, end=912, iostat=ios)
        h_film=.true.
        IF (dvac > dtild)  CALL juDFT_error("dvac > dtild",calledby
     +       ="lapw_input")

!===> kpt
      ELSEIF (buffer(1:4)=='&kpt') THEN
        IF (h_kpt) CALL err(1)

        div1 = 1 ; div2 = 1 ; div3 = 1

        READ (bfh,kpt,err=912, end=912, iostat=ios)
        h_kpt=.true.
        div(1) = div1 ; div(2) = div2 ; div(3) = div3 

!===> end
      ELSEIF (buffer(1:4)=='&end') THEN
        WRITE (dbgfh,*) 'end of input record in line:',nline
        EXIT loop

!===> atom, allatoms
      ELSEIF (buffer(1:5)=='&atom' .OR.
     &        buffer(1:9)=='&allatoms') THEN
        WRITE(errfh,*) 'buffer ',buffer(1:9),
     &       ' out of place in line:',nline
        fatalerror = .true.

!===> unknown namelist
      ELSE
          call err(2)
      ENDIF
!===>
      nbuffer = 0

      ENDDO loop

 999  CONTINUE
      IF (fatalerror) 
     &     CALL juDFT_error
     +     ("ERROR(S) reading input. Check output for details.",calledby
     +     ="lapw_input")

      RETURN

!===> error handling

 911  CONTINUE
      WRITE (errfh,*) 'atom_input: ERROR reading input. ios  =',ios,
     &               ', line =',nline
      CALL juDFT_error("atom_input: ERROR reading input",calledby
     +     ="lapw_input")

 912  CONTINUE
      WRITE (errfh,*) 'atom_input: ERROR reading namelist.',
     &               ' ios =',ios,
     &               ' line =',nline
      WRITE (errfh,*) buffer(1:nbuffer)
      WRITE (errfh,*) 'The cause of this error may be ...'
      WRITE (errfh,*) '        a variable not defined in this namelist,'
      WRITE (errfh,*) '        wrong type of data for a variable.'
      CALL juDFT_error("atom_input: ERROR reading input",calledby
     +     ="lapw_input")

 913  CONTINUE
      WRITE (errfh,*) 'atom_input: ERROR reading record.',
     &               ' ios =',ios,
     &               ' line =',nline
      WRITE (errfh,*) buffer(1:nbuffer)
      CALL juDFT_error("atom_input: ERROR reading input",calledby
     +     ="lapw_input")

!----------------------------------------------------------------
      CONTAINS   ! INTERNAL subroutines
!----------------------------------------------------------------
      SUBROUTINE err( n )

      INTEGER, INTENT (IN) :: n

      WRITE(errfh,*)
      IF (n==1) THEN
        WRITE (errfh,*) 'atom_input: ERROR multiple namelists.',
     &               ' line =',nline
      ELSEIF (n==2) THEN
        WRITE (errfh,*) 'atom_input: ERROR unknown namelist.',
     &               ' line =',nline
      ELSEIF (n==3) THEN
        WRITE (errfh,*) 'atom_input: ERROR reading namelist.',
     &               ' line =',nline
      ELSE
        WRITE (errfh,*) 'atom_input: ERROR reading input.',
     &               ' line =',nline
      ENDIF
      WRITE (errfh,*) buffer(1:nbuffer)
      WRITE (errfh,*)
      fatalerror = .true.
      RETURN
      END SUBROUTINE err

      END SUBROUTINE lapw_input
      END MODULE m_lapwinput      
