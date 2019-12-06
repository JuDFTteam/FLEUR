      MODULE m_structinput
      use m_juDFT
!********************************************************************
!     read in lattice information and generate space group operations
!********************************************************************
      CONTAINS
      SUBROUTINE struct_input( 
     >                        infh,errfh,warnfh,symfh,symfn,bfh,
     >                        natmax,nop48,
     X                        nline,xl_buffer,buffer,
     <                        title,film,cal_symm,checkinp,symor,
     <                    cartesian,oldfleur,a1,a2,a3,dvac,aa,scale,i_c,
     <                       factor,natin,atomid,atompos,ngen,mmrot,ttr,
     <                       atomLabel,
     <                        l_hyb,l_soc,l_ss,theta,phi,qss,inistop)

      use m_calculator
      USE m_readrecord
      USE m_rwsymfile
      USE m_lattice, ONLY : lattice2
      IMPLICIT NONE

!===> Arguments
      INTEGER, INTENT (IN)    :: infh, errfh,  warnfh, symfh,bfh
      INTEGER, INTENT (IN)    :: natmax, nop48, xl_buffer
      INTEGER, INTENT (INOUT)  :: nline
      CHARACTER(len=xl_buffer) :: buffer
      LOGICAL                 :: cal_symm, checkinp, symor, film
      LOGICAL                 :: cartesian,oldfleur,inistop
      LOGICAL, INTENT (OUT)   :: l_hyb,l_soc,l_ss
      INTEGER, INTENT (OUT)   :: natin,i_c
      INTEGER, INTENT (OUT)   :: ngen
      REAL,    INTENT (OUT)   :: aa,theta,phi
      REAL,    INTENT (OUT)   :: dvac
      REAL,    INTENT (OUT)   :: a1(3),a2(3),a3(3)
      REAL,    INTENT (OUT)   :: scale(3),factor(3),qss(3)
      REAL,    INTENT (OUT)   :: atompos(3,natmax)
      REAL,    INTENT (OUT)   :: atomid(natmax)
      INTEGER, INTENT (OUT)   :: mmrot(3,3,nop48)
      REAL,    INTENT (OUT)   :: ttr(3,nop48)
      CHARACTER(len=80), INTENT (OUT) :: title
      CHARACTER(len=7),  INTENT (IN)  :: symfn
      CHARACTER(LEN=20), INTENT (OUT) :: atomLabel(natmax)

!===> data
      REAL,             PARAMETER :: eps=1.e-7
      CHARACTER(len=1), PARAMETER :: cops(-1:3)=(/'2','3','4','6','1'/)

!===> Local Variables
      INTEGER :: n, ng, op, nbuffer, ios,nop2, i, j
      REAL    :: shift(3),rdummy(3,3),z_max,z_min,mat(3,3),x(3), rest
      LOGICAL :: oldfleurset,l_symfile,l_gen,hybrid
      CHARACTER(len=10)        :: chtmp
      CHARACTER(len=3)         :: ch_test
!===> namelists
      NAMELIST /input/ film, cartesian, cal_symm, checkinp, inistop,
     &                 symor, oldfleur, hybrid

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!--->    set defaults
      film      = .false.   ! bulk calculation is default

      cartesian = .false.   ! read in atomic positions 
                            ! in either lattice units (.false.)
                            ! or scaled cartesian coordinates (.true.)

      cal_symm  = .true.    ! calculate space group symmetry    (.true.)
                            ! read in space group symmetry info (.false.)

      checkinp  = .false.   ! =T program reads input and stops
      
      inistop = .FALSE.     ! = T program stops after strho,swsp,flip

      op = -911             ! =1 => checkinp=t
                            ! =2 => inistop=t
                            ! =4 => itmax=0

      oldfleur  = .true.    ! =T fleur21 compatibility

      symor     = .false.   ! =T select the largest symmorphic subgroup

      hybrid    = .false.   ! =T create inp file for hybrid functionals, too
      l_ss      = .false.   ! =T spin-spiral calculation ... may affect
      l_soc     = .false.   ! =T spin-orbit interaction ... the symmetry

      factor(:) = 1.0
      theta = 0.0 ; phi = 0.0
      qss(:) = (/0.0,0.0,0.0/)
      mat = 0.0
!initialize the calculator
      DO
         READ (UNIT = infh,FMT = "(a3)",iostat=ios) ch_test
         if (ios.ne.0) exit
         BACKSPACE(infh)
         IF (ch_test   /="def") EXIT
         READ(unit = infh,FMT="(4x,a)") buffer
         n = INDEX(buffer,"=")
         IF (n == 0.OR.n>len_TRIM(buffer)-1) STOP
     $        "Error in variable definitions"
         CALL ASSIGN_var(buffer(:n-1),evaluate(buffer(n+1:)))
      ENDDO
      backspace(infh)


!===> start reading input

      CALL read_record(
     >                 infh,xl_buffer,bfh,
     X                 nline,
     <                 nbuffer,buffer,ios )

      READ (buffer,'(a)') title

      IF ( buffer(1:1) == '&' ) THEN         ! already read some input
        title = 'unnamed project'
      ELSE   
        CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
      ENDIF

      oldfleurset = .FALSE.
      l_hyb = .FALSE.
      IF ( buffer(1:6)=='&input' ) THEN      ! get namelist 'input'
        READ (bfh,input)
        l_hyb       = hybrid
        IF ( index(buffer,'oldfleur')>0 ) oldfleurset = .true.
        op = 0 
        IF ( op > 0 ) THEN
          IF ( btest(op,0) ) checkinp = .true.
          IF ( BTEST(op,1) ) inistop  = .TRUE.
          IF ( btest(op,2) ) WRITE (6,*) 'action N/A'
!dbg+
          IF ( btest(op,0) ) WRITE (6,*) 'bit 0 set'
          IF ( btest(op,1) ) WRITE (6,*) 'bit 1 set'
          IF ( btest(op,2) ) WRITE (6,*) 'bit 2 set'
!dbg-
        ENDIF

        CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
      ENDIF

      IF ( buffer(1:1) == '&' ) THEN

        IF ( buffer(1:8)=='&lattice' ) THEN
          CALL lattice2( 
     >                  buffer,xl_buffer,errfh,bfh,nline,
     <                  a1,a2,a3,aa,scale,mat,i_c,ios )
          dvac = 0.00
          IF ( ios.NE.0 ) THEN
            WRITE (errfh,*)
            WRITE (errfh,*) 'struct_input: ERROR! ',
     &                   'while reading &lattice in line',nline,'.'
            WRITE (errfh,*)
            CALL juDFT_error("ERROR! while reading &lattice",calledby
     +           ="struct_input")
          ENDIF
        ELSE
          WRITE(errfh,*)
          WRITE(errfh,*) 'struct_input: ERROR! line',nline,'.'
          WRITE(errfh,*)
     &     'Expecting either namelist &lattice or dircet lattice input.'
          WRITE (errfh,*)
          CALL juDFT_error("ERROR! Cannot find lattice info.",calledby
     +         ="struct_input")
        ENDIF

      ELSE

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!--->    title:
!        cannot begin with an & and should not contain an !
!
!--->    lattice vectors:
!--->
!--->    lattice vectors are input in scaled cartesian coordinates:
!--->
!--->    the overall scale is set by aa and scale(:) as follows:
!--->    assume that we want the lattice vectors to be given by
!--->      a_i = ( a_i(1) xa , a_i(2) xb , a_i(3) xc )
!--->    then choose aa, scale such that: xa = aa * scale(1), etc.
!--->    to make it easy to input sqrts, if scale(i)<0, then
!--->    scale = sqrt(|scale|)
!--->    Example: hexagonal lattice
!           a1 = ( sqrt(3)/2 a , -1/2 a , 0.      )
!           a2 = ( sqrt(3)/2 a ,  1/2 a , 0.      )
!           a3 = ( 0.          , 0.     , c=1.62a )
!
!        input:
!            0.5  -0.5  0.0     ! a1
!            0.5   0.5  0.0     ! a2
!            0.0   0.0  1.0     ! a3
!           6.2                 ! lattice constant
!          -3.0   0.0   1.62    ! scale(2) is 1 by default

!--->    read in (scaled) lattice vectors (and dvac, if present)

         !READ (buffer,*) a1
         a1(1)=evaluatefirst(buffer)
         a1(2)=evaluatefirst(buffer)
         a1(3)=evaluatefirst(buffer)
         CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
         !READ (buffer,*) a2
         a2(1)=evaluatefirst(buffer)
         a2(2)=evaluatefirst(buffer)
         a2(3)=evaluatefirst(buffer)
         CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
         !READ (buffer,*, err=811,end=811, iostat=ios) a3, dvac
         a3(1)=evaluatefirst(buffer)
         a3(2)=evaluatefirst(buffer)
         a3(3)=evaluatefirst(buffer)
         dvac=evaluatefirst(buffer)
         IF (film.AND.(dvac <= 0.00)) THEN
            WRITE(*,*)'Film calculation but no reasonable dVac provided'
            WRITE(*,*)'Setting default for dVac'
            dvac = ABS(a3(3)) ! This is later set to the real default by the chkmt result
         END IF
 811     CONTINUE              ! obviously no film calculation
         !READ(buffer,*) a3

!--->    read in overall lattice constant

         CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
         !READ (buffer,*) aa
         aa=evaluatefirst(buffer)

!--->    read in scale
         scale = 0.00
         CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
         !READ (buffer,*) scale
         scale(1)=evaluatefirst(buffer)
         scale(2)=evaluatefirst(buffer)
         scale(3)=evaluatefirst(buffer)
         
      ENDIF ! &lattice ...

!===>    program configuration
!     if oldfleur was not set in the input, set it here, dependent on 
!     film/bulk calculation

      IF ( .not.oldfleurset ) THEN
        oldfleur = .false.
        IF ( film ) oldfleur = .true.
      ENDIF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!--->    atomic positions:
!--->
!--->    atomic positions input can be either in scaled cartesian
!--->    or lattice vector units, as determined by logical cartesian.
!--->    (for supercells, sometimes more natural to input positions
!--->    in scaled cartesian.)
!--->
!--->    if ntin < 0, then the representative atoms only are given;
!--->    this requires that the space group symmetry be given as input.

!--->    read in number of atoms or types

      CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
      !READ (buffer,*) natin
      natin=evaluatefirst(buffer)
!--->    read in atomic positions
!--->    and atomic identification number (atomid)
!--->    to distinguish different atom types. 
!--->    (atomid is used later as default for atom Z value (zatom)

      DO n = 1, abs(natin)
        CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
        !READ (buffer,*) atomid(n), atompos(:,n)
        atomid(n)=evaluatefirst(buffer)
        atompos(1,n)=evaluatefirst(buffer)
        atompos(2,n)=evaluatefirst(buffer)
        atompos(3,n)=evaluatefirst(buffer)
        IF(TRIM(ADJUSTL(buffer)).NE.'') THEN
           atomLabel(n) = TRIM(ADJUSTL(buffer))
        ELSE
           WRITE(atomLabel(n),'(i0)') n
        END IF
        DO i = 2,40
           DO j = 1, 3
              rest = ABS(i*atompos(j,n) - NINT(i*atompos(j,n)))
              IF (rest.LT.(i*0.000001)) THEN
                 atompos(j,n) = NINT(i*atompos(j,n)) / real(i)
              END IF
           END DO
        END DO
      ENDDO

      CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)

      IF ( buffer(1:6)=='&shift') then
        buffer = buffer(7:len_trim(buffer)-1)
        shift = -911.0
        READ (buffer,*, err=821,end=821, iostat=ios) shift
 821    CONTINUE
        READ (buffer,*) shift(1)
        IF ( shift(3)<-900.0 ) shift(3) = shift(1)
        IF ( shift(2)<-900.0 ) shift(2) = shift(1)
        DO n = 1, 3
          atompos(n,1:abs(natin)) = atompos(n,1:abs(natin))+shift(n)
        ENDDO
        
        CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
      ENDIF

      IF ( buffer(1:7)=='&factor') THEN
        buffer = buffer(8:len_trim(buffer)-1)
        factor = -911.0 
        READ (buffer,*, err=831,end=831, iostat=ios) factor
 831    CONTINUE
        READ (buffer,*) factor(1)
        IF ( factor(3)<-900.0 ) factor(3) = factor(1)
        IF ( factor(2)<-900.0 ) factor(2) = factor(1)
        DO n = 1, 3
          atompos(n,1:abs(natin)) = atompos(n,1:abs(natin))/factor(n)
        ENDDO
        
        CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
      ENDIF

      IF ( buffer(1:1).NE.'&') THEN
        WRITE (warnfh,*)
        WRITE (warnfh,*) 'struct_input: WARNING! ',
     &       'Number of atoms to small or too many atoms in list?.'
        WRITE (warnfh,*)
      ENDIF

      IF (abs(mat(1,1)).GT.0.0000001) THEN ! transform hex->trig
        CALL recip(a1,a2,a3,rdummy)
        DO n = 1, abs(natin)
          x = matmul(mat,atompos(:,n))
          write(*,'(3f10.5)') x(1:3)
          write(*,'(3f10.5)') rdummy
          atompos(:,n) = matmul(rdummy,x)
          write(*,'(3f10.5)') atompos(1:3,n)
        ENDDO
        i_c = 1
      ENDIF

      IF (film) THEN

        z_max = MAXVAL( atompos(3,1:abs(natin)) )  ! check the outmost atomic position
        z_min = MINVAL( atompos(3,1:abs(natin)) )  
        z_max = 2 * (MAX(z_max,-z_min) + 3.0)      ! how much space do we need in z-dir.
        a3(3) = MAX( a3(3), z_max/(aa*scale(3)) )  ! adjust a3(3) so that it fits

        IF(.NOT.cartesian) THEN
           atompos(3,1:abs(natin)) =                  ! rescale to internal coordinates
     +     atompos(3,1:abs(natin))/(a3(3)*aa*scale(3))
        END IF

      ENDIF

!===> read symmetry from file or from namelist

      INQUIRE ( file=trim(symfn), exist=l_symfile )

      IF ( l_symfile ) THEN

        WRITE (6,*) 'DBG: l_symfile=',l_symfile
        CALL rw_symfile(
     >                  'r',symfh,symfn,nop48,rdummy,
     X                   mmrot,ttr,ngen,nop2,symor)
        cal_symm = .false.

      ELSE

        l_gen = .false.
        IF ( buffer(1:4)=='&gen' ) l_gen = .true.

        IF ( buffer(1:4)=='&gen' .or.
     &       buffer(1:4)=='&sym'     ) THEN

          WRITE (6,*) 'DBG: &sym=',buffer(1:4)

          buffer  = ADJUSTL(buffer(5:nbuffer))
          nbuffer = LEN_TRIM(buffer)

          READ (buffer,*,err=913, end=913, iostat=ios) ngen
          n = scan(buffer,' ')
          IF ( n>0 ) buffer = buffer(n+1:nbuffer)
          buffer  = adjustl(buffer(1:nbuffer))
          nbuffer = len_trim(buffer)
          READ (buffer,*,err=913, end=913, iostat=ios) 
     &         ( mmrot(1,1,n),mmrot(1,2,n),mmrot(1,3,n),ttr(1,n),
     &           mmrot(2,1,n),mmrot(2,2,n),mmrot(2,3,n),ttr(2,n),
     &           mmrot(3,1,n),mmrot(3,2,n),mmrot(3,3,n),ttr(3,n),
     &           n = 2, ngen+1 )

        cal_symm = .false.

        CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)

        ENDIF ! &gen or &sym

      ENDIF ! from_symfile
      IF ( buffer(1:4)=='&soc' ) THEN
         l_soc=.true. 
         buffer  = ADJUSTL(buffer(5:nbuffer))
         nbuffer = LEN_TRIM(buffer)
         READ (buffer,*,err=913, end=913, iostat=ios) theta,phi
         CALL read_record(infh,xl_buffer,bfh,nline,nbuffer,buffer,ios)
      ENDIF
      IF ( buffer(1:4)=='&qss' ) THEN
         l_ss=.true.
         buffer  = ADJUSTL(buffer(5:nbuffer))
         nbuffer = LEN_TRIM(buffer)
         READ (buffer,*,err=913,end=913,iostat=ios) qss(1),qss(2),qss(3)
      ENDIF


      IF ( .not.cal_symm ) THEN   ! &gen or &sym

!---> make sure idenity is first operation          
        mmrot(:,:,1) = reshape((/ 1,0,0, 0,1,0, 0,0,1 /),(/ 3,3 /))
        ttr(:,1) = 0.0

        DO ng = 2, ngen + 1
          IF ( all( mmrot(:,:,ng)==mmrot(:,:,1) ) .AND.          ! identity was entered 
     &         all( abs( ttr(:,ng)-ttr(:,1) ) < eps ) ) THEN     ! explicitely as matrix 'ng'
            DO n = ng, ngen
              mmrot(:,:,n) = mmrot(:,:,n+1)                      ! shift by '-1' & exit
              ttr(:,n)     = ttr(:,n+1) 
            ENDDO
            ngen = ngen - 1
            EXIT
          ENDIF
        END DO

        IF ( l_gen ) mmrot(1,1,1) = 0 ! is used later to distinguish
                                      ! between generators and full group

        WRITE (6,*) 'DBG: mmrot(1,1,1)=',mmrot(1,1,1)

      ENDIF ! .not.cal_symm

      RETURN

 912  CONTINUE
      WRITE (errfh,*) 'struct_input: ERROR reading namelist.',
     &               ' ios =',ios,
     &               ' line =',nline
      WRITE (errfh,*) buffer(1:nbuffer)
      WRITE (errfh,*) 'The cause of this error may be ...'
      WRITE (errfh,*) '        a variable not defined in this namelist,'
      WRITE (errfh,*) '        wrong type of data for a variable.'
      CALL juDFT_error("ERROR reading input",calledby ="struct_input")

 913  CONTINUE
      WRITE (errfh,*) 'struct_input: ERROR reading record.',
     &               ' ios =',ios,
     &               ' line =',nline
      WRITE (errfh,*) buffer(1:nbuffer)
       CALL juDFT_error("ERROR reading input",calledby="struct_input")

      END SUBROUTINE struct_input
!-------------------------------------
      SUBROUTINE recip(a1,a2,a3,b)

      USE m_constants, ONLY : pimach
      IMPLICIT NONE
      REAL, INTENT (IN) :: a1(3),a2(3),a3(3)
      REAL, INTENT (OUT):: b(3,3)
      REAL volume

!  volume (without scaling factor aa^3)
      volume  = a1(1)*a2(2)*a3(3) + a2(1)*a3(2)*a1(3) +
     &          a3(1)*a1(2)*a2(3) - a1(3)*a2(2)*a3(1) -
     &          a2(3)*a3(2)*a1(1) - a3(3)*a1(2)*a2(1)

!  reciprocal lattice vectors in scaled Cartesian units
      b(1,1) = (a2(2)*a3(3) - a2(3)*a3(2))
      b(1,2) = (a2(3)*a3(1) - a2(1)*a3(3))
      b(1,3) = (a2(1)*a3(2) - a2(2)*a3(1))
      b(2,1) = (a3(2)*a1(3) - a3(3)*a1(2))
      b(2,2) = (a3(3)*a1(1) - a3(1)*a1(3))
      b(2,3) = (a3(1)*a1(2) - a3(2)*a1(1))
      b(3,1) = (a1(2)*a2(3) - a1(3)*a2(2))
      b(3,2) = (a1(3)*a2(1) - a1(1)*a2(3))
      b(3,3) = (a1(1)*a2(2) - a1(2)*a2(1))
!      b = 2.0*pimach()*b/volume
      b = b/volume
      END SUBROUTINE recip

      END MODULE m_structinput
