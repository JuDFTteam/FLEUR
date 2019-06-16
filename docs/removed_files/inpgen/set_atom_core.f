      MODULE m_setatomcore
      use m_juDFT
      CONTAINS
!================================
!     setatom_bystr
!     setcore_bystr
!================================
      SUBROUTINE setatom_bystr(
     >                        l_buffer,nwdd,econfig,
     <                        natomst, ncorest, nvalst, nelec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This subroutine sets up the atomic configuration, definded by a string
! In principle all features of 'setatom' are implemented, except for the
! storage of the configuration in the arrays nprnc, kappa, occ
!
! Different states are defined by n,l and the occupation
! the syntax for this string is rather strict:
!     nloo   where n is the single digit major quantum number n 1,2,3,...
!                  l is the single digit minor quantum number l s,p,d,f
!                 oo is the double digit integer occupation 1,2,...,10,11,...
!     the entries in the string are separated by spaces
! Core electrons and valence electrons and different energy-windows 
! are separated by one of 'separators'
!
!                                                           r.g.2001
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE

      INTEGER,                 INTENT(IN)    :: l_buffer,nwdd
      CHARACTER(LEN=l_buffer), INTENT(IN)    :: econfig
      INTEGER,       OPTIONAL, INTENT(OUT)   :: natomst,ncorest,nvalst
      REAL,          OPTIONAL, INTENT(INOUT) :: nelec(0:nwdd)

      INTEGER     :: n,o,win,newst
      CHARACTER   :: l
      INTEGER     :: nat, nco, nva, nel(0:nwdd)
      CHARACTER(LEN=l_buffer) :: workconf

      IF ( TRIM(econfig).EQ."NULL".OR.              ! check, wether we're 
     &     LEN_TRIM(econfig).EQ.0 )  RETURN         ! supposed to do something

      nat = 0 ! overall number of atomic states
      nco = 0 ! number of corestates
      nva = 0 ! number of valence states
      win = 0 ! index to current window (0=core)
      nel = 0 ! initialize electron-counter

      workconf=econfig            ! working copy of econfig
      CALL expandconfig(workconf) ! check econfig and expand noble-gas shortcuts

      DO WHILE ( LEN_TRIM(workconf).GT.0 ) ! scan through character string 
                                           ! containing the electronic configuration

        CALL getconfig(workconf,           ! get next item of the configuration
     <                 n,l,o)

        DO WHILE ( n.EQ.0.AND.l.EQ."0".AND.o.EQ.0 )  ! if a window-deliminator was hit,
                                                     ! increase the window=index
          win=win+1
          CALL getconfig(workconf,
     <                   n,l,o)
        ENDDO

        SELECT CASE (l)                    ! store n,l,o to the atomic arrays and 
                                           !increase the state-counters accordingly      
          CASE("s")
                IF ( o.GT.0 ) newst=1
          CASE("p")
                IF ( o.GT.0 ) newst=1
                IF ( o.GT.2 ) newst=2
          CASE("d") 
                IF ( o.GT.0 ) newst=1
                IF ( o.GT.4 ) newst=2
          CASE("f")
                IF ( o.GT.0 ) newst=1
                IF ( o.GT.6 ) newst=2
        END SELECT
        nat=nat+newst
        IF ( win.EQ.0 ) THEN ; nco=nco+newst
                        ELSE ; nva=nva+newst
        ENDIF
        nel(win)=nel(win)+o
      ENDDO

      IF ( PRESENT(natomst) ) natomst=nat
      IF ( PRESENT(ncorest) ) ncorest=nco
      IF ( PRESENT(nvalst)  ) nvalst =nva
      IF ( PRESENT(nelec)   ) nelec  =nel

      RETURN
      END SUBROUTINE setatom_bystr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  SUBROUTINE setcore_bystr
!
! This subroutine sets up the core configuration, definded by a string
! It stores the core-quantum-numbers to coreqn and the 
! occupation-numbers to coreocc, which will be used by 'setcore'
!
! Different states are defined by n,l and the occupation
! The syntax for this string is rather strict:
!     nloo   where n is the single digit major quantum number n 1,2,3,...
!                  l is the single digit minor quantum number l s,p,d,f
!                 oo is the double digit integer occupation 1,2,...,10,11,...
!     The entries in the string are separated by spaces
! core electrons and valence electrons are separated by one of 'separators'
!
! the routine will not exit as the first window-separator is found, but 
! count up all levels. Core states reach up to nallst
!
!                                                           r.g.12.2000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE setcore_bystr(
     >                         atype,nstd,ntype,l_buffer,
     X                         econfig,nallst,ncorest,
     <                         coreqn,coreocc)

      IMPLICIT none

      INTEGER, PARAMETER :: errfh = 6
!      INTEGER          :: storeconfig

      INTEGER, INTENT(IN)     :: atype,nstd,ntype,l_buffer
      CHARACTER(LEN=l_buffer), INTENT(IN)    :: econfig(ntype)
      INTEGER, INTENT(INOUT)  :: ncorest,nallst
      INTEGER, INTENT(OUT)    :: coreqn(2,nstd,ntype) ! core states (relativistic)
      REAL,    INTENT(OUT)    :: coreocc(nstd,ntype)  ! core occupations

      INTEGER                 :: n,o,i,win,newst
      CHARACTER               :: l
      CHARACTER(LEN=l_buffer) :: workconf

! check, wether we're supposed to do something
      IF (     TRIM(econfig(atype)).EQ."NULL".OR.
     &     LEN_TRIM(econfig(atype)).EQ.0 ) RETURN

      nallst = 0                  ! number of all states
      ncorest = 0                 ! number of corestates
      workconf = econfig(atype)   ! working copy of electronic configuration
      CALL expandconfig(workconf) ! check econfig and expand noble-gas shortcuts

      win = 0
      DO WHILE ( LEN_TRIM(workconf).GT.0 )    ! scan through character string 
                                              ! containing the electronic configuration
        CALL getconfig(workconf,n,l,o)        ! get next item of the configuration

        DO WHILE ( n.EQ.0.AND.l.EQ."0".AND.o.EQ.0 )  ! if a window-deliminator was hit
          win = win + 1                              ! increase window index
          CALL getconfig(workconf,n,l,o)
        ENDDO

        newst=storeconfig(n,l,o,nstd,nallst,                        ! store n,l,o to the atomic
     &    coreqn(1,:,atype),coreqn(2,:,atype),coreocc(:,atype))     ! arrays and increase the
        nallst  = nallst  + newst                                   ! state-counters accordingly
        IF (win == 0) ncorest = ncorest + newst                                   

      ENDDO

      IF (win == 0) THEN
! actually, the routine should be exited from within the WHILE-loop
! obviously no window-sepatator was found - ABORT!
      WRITE(errfh,*) "An error was found while processing the"
      WRITE(errfh,*) "  electronic configuration for atom #",
     &               atype
      WRITE(errfh,*) "Please divide core- and valence states "
      WRITE(errfh,*) "  by inserting one of '|/\-'"
       CALL juDFT_error("setcore_bystr",calledby="set_atom_core",hint=
     +     "An error was found while processing the"//
     +     " electronic configuration. Please divide core- and "//
     +     "valence states by inserting one of '|/\-'")
      ENDIF

      END SUBROUTINE setcore_bystr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  SUBROUTINE expandconfig
!
! this subroutine check the validity of econfig (not very strict though) and
! expands noble gas shortcuts
!                                                           r.g.12.2000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE expandconfig(econfig)

      IMPLICIT none

      INTEGER, PARAMETER :: errfh = 6, l_buffer= 512, s_buffer= 20
      CHARACTER (LEN=*), PARAMETER :: 
     C    He_core="1s2",
     C    Ne_core="1s2 2s2 2p6",
     C    Ar_core="1s2 2s2 2p6 3s2 3p6",
     C    Kr_core="1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6",
     C    Xe_core="1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6",
     C    Rn_core="1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 "//
     C                                "5p6 5d10 6s2 6p6"
      CHARACTER (LEN=*), PARAMETER :: 
     C    separators="|/\-",
     C    legalchars="spdf0123456789 [HeNeArKrXeRn]"//separators

      CHARACTER(LEN=*), INTENT(INOUT) :: econfig

      INTEGER :: n,o,error,i,win
      INTEGER :: strpos,newpos
      CHARACTER(LEN=l_buffer)      :: bigbuf
      CHARACTER(LEN=s_buffer)      :: buf
      CHARACTER                    :: l

! scan for illegal charcters
      error=VERIFY(econfig,legalchars)
      IF ( error.NE.0 ) THEN
        WRITE(errfh,'(A)')  "*** An illegal character was found in "//
     &                      "the electronic configuration:"
        WRITE(errfh,'(A)')  TRIM(econfig)
        bigbuf=REPEAT("-",(error-1))//"^"
        WRITE(errfh,FMT='(A)') TRIM(bigbuf)
        CALL juDFT_error("invalid electronic configuration",calledby
     +       ="set_atom_core")
      ENDIF

      IF ( SCAN(econfig,separators).EQ.0 ) THEN
! obviously no window-sepatator was found - ABORT!
        WRITE(errfh,*) "An error was found while processing the"
        WRITE(errfh,*) "  electronic configuration ",TRIM(econfig)
        WRITE(errfh,*) "Please divide core- and valence states "
        WRITE(errfh,*) "  by inserting one of '|/\-'"
       CALL juDFT_error("expnadconfig",calledby="set_atom_core",hint=
     +     "An error was found while processing the"//
     +     " electronic configuration. Please divide core- and "//
     +     "valence states by inserting one of '|/\-'")
      ENDIF

! condense multiple spaces to single spaces
      econfig=ADJUSTL(econfig)
      DO WHILE ( INDEX(TRIM(econfig),"  ").GT.0 )
        strpos=INDEX(econfig,"  ")
        econfig=econfig(:strpos)//econfig(strpos+2:)
      END DO
 
! first, look if there's some noble gas core configuration included
! if so, replace the shortcut with the full configuration
      if ( INDEX(econfig,"[").GT.0 ) THEN
        strpos=INDEX(econfig,"[")
        newpos=INDEX(econfig,"]")
        SELECT CASE (econfig((strpos+1):(newpos-1)))
          CASE ("He")
            econfig=He_core//" "//ADJUSTL(econfig((newpos+1):))
          CASE ("Ne")
            econfig=Ne_core//" "//ADJUSTL(econfig((newpos+1):))
          CASE ("Ar")
            econfig=Ar_core//" "//ADJUSTL(econfig((newpos+1):))
          CASE ("Kr")
            econfig=Kr_core//" "//ADJUSTL(econfig((newpos+1):))
          CASE ("Xe")
            econfig=Xe_core//" "//ADJUSTL(econfig((newpos+1):))
          CASE ("Rn")
            econfig=Rn_core//" "//ADJUSTL(econfig((newpos+1):))
          CASE DEFAULT
            WRITE(*,'(2A)') "*** Unknown noble gas ",
     A                      econfig((strpos+1):(newpos-1))
            CALL juDFT_error("unknown configuration shortcut",calledby
     +           ="set_atom_core")
        END SELECT
      ENDIF

      RETURN
      END SUBROUTINE expandconfig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  SUBROUTINE getconfig
!
! this subroutine chops off the first set of electrons of econfig and returns 
! the corresponding n,l,o values and the shortened econfig
!
! the syntax for this string is rather strict:
!     nloo   where n is the single digit major quantum number n 1,2,3,...
!                  l is the single digit minor quantum number l s,p,d,f
!                 oo is the double digit integer occupation 1,2,...,10,11,...
!     the entries in the string are separated by spaces
!
! core electrons and valence electrons are separated by |,...
! if such a separator was encountered, n=l=o=0 is returned
!
! make sure, the configuration covers _all_ electrons - core and valence
!
!                                                           r.g.12.2000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getconfig(econfig,n,l,o)

      IMPLICIT NONE

      INTEGER, PARAMETER :: errfh = 6, l_buffer= 512, s_buffer= 20
      CHARACTER (LEN=*), PARAMETER :: separators="|/\-"

      CHARACTER(LEN=*), INTENT(INOUT) :: econfig
      INTEGER,          INTENT(OUT)   :: n,o
      CHARACTER,        INTENT(OUT)   :: l

      INTEGER :: error,i,win
      INTEGER :: strpos,newpos
      CHARACTER(LEN=l_buffer)      :: bigbuf
      CHARACTER(LEN=s_buffer)      :: buf

      win = 0 !gs

      econfig=ADJUSTL(econfig)
! look if a separator is in first place. if so, exit
      IF ( VERIFY(econfig(1:1),separators).EQ.0 ) THEN
        win=win+1
        econfig=ADJUSTL(econfig(2:))
        n=0
        l="0"
        o=0
        RETURN
      ENDIF

! copy first item to buf - 
! items are separated by a space
! or an separator (see separators), to mark core/valence and win/win - bounds
      strpos=SCAN(econfig,(separators//" "))
!      write(*,*) "getconfig: ",TRIM(econfig)
!      write(*,*) "           ",REPEAT("-",(strpos-1)),"^"
      buf=econfig(:(strpos-1))
!      write(*,*) "getconfig: ",TRIM(buf)
! check if these n-l values occur again within the string
! if so, abort
      IF ( INDEX(econfig(2:),buf(1:2)).GT.0 ) THEN
        WRITE(*,'(3A)') "*** There can only be one set of ",buf(1:2),
     A             "-electrons"
        CALL juDFT_error("getconfig: invalid electronic configuration"
     +       ,calledby ="set_atom_core")
      ENDIF
! split and convert description to numeric/character var's
      READ(buf,FMT='(I1,A1,I2)',IOSTAT=error) n,l,o
      IF ( error.NE.0 ) THEN
        WRITE(*,'(2A)') "*** error encountered, "//
     A                  "while processing entry ",buf," of econfig"
        WRITE(*,'(A)') "*** valid syntax is  nloo  (e.g. 3d10)"
        CALL juDFT_error("error converting configuration string"
     +          ,calledby ="set_atom_core")
      ENDIF

! remove first item from econfig
      econfig=ADJUSTL(econfig(strpos:))

      RETURN

      END SUBROUTINE getconfig


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      INTEGER FUNCTION storeconfig
!
!  This function stores the states described by n,l,o (inn,inl,ino) to the 
!  arrays coren(nstd) (n) corek(nstd) (kappa) and coreocc(nstd)
!
!  nst is the number of the current number of states in these arrays
!  it will _NOT_ be changed here. 
!
!  the return value is the number of states added to the arrays
!
!                                                          r.g.2001
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER FUNCTION storeconfig(
     >                             inn,inl,ino,nstd,nst,
     <                             coren,corek,coreo)

      IMPLICIT none

      INTEGER,       INTENT(IN)    :: inn,ino,nstd
      CHARACTER,     INTENT(IN)    :: inl
      INTEGER,       INTENT(IN)    :: nst
      INTEGER,       INTENT(INOUT) :: coren(nstd)
      INTEGER,       INTENT(INOUT) :: corek(nstd)
!      INTEGER,       INTENT(INOUT) :: coreqn(2,nstd)
      REAL,          INTENT(INOUT) :: coreo(nstd)

      INTEGER            :: naddst,n,o
      CHARACTER(LEN=1)   :: l

      naddst=0
      n=inn
      l=inl
      o=ino

      SELECT CASE (TRIM(l))
! --- s-states
        CASE("s")
          IF ( o.GT.2 ) THEN
            WRITE(*,'(A,I3,A)') "*** ",o,
     A                          " s-electrons are impossible!"
            CALL juDFT_error("storeconfig: too many s-electrons",
     +           calledby="set_atom_core")
          ENDIF
          naddst=naddst+1
          coren(nst+naddst) = n
          corek(nst+naddst) = -1
          coreo(nst+naddst)  = REAL(o)
! --- p-states
        CASE("p")
          IF ( o.GT.6 ) THEN
            WRITE(*,'(A,I3,A)') "*** ",o,
     A                          " p-electrons are impossible!"
            CALL juDFT_error("storeconfig: too many p-electrons",
     +           calledby="set_atom_core")
          ENDIF
          IF ( n.LT.2 ) THEN 
            WRITE(*,'(a,i0,a)') "*** There's nothing like ",n,"p"
            CALL juDFT_error("storeconfig: invalid specifier",
     +                      calledby="set_atom_core")
          ENDIF
          naddst=naddst+1
          coren(nst+naddst)  = n
          corek(nst+naddst)  = 1
          IF ( o.GT.2 ) THEN; coreo(nst+naddst) = 2.0 
                        ELSE; coreo(nst+naddst) = REAL(o)
          ENDIF
          o = o-2
          IF ( o.GT.0 ) THEN
            naddst=naddst+1
            coren(nst+naddst) = n
            corek(nst+naddst) = -2
            coreo(nst+naddst)  = REAL(o)
          ENDIF
! --- d-states
        CASE("d")
          IF ( o.GT.10 ) THEN
            WRITE(*,'(A,I3,A)') "*** ",o,
     A                          " d-electrons are impossible!"
            CALL juDFT_error("storeconfig: too many d-electrons",
     +           calledby="set_atom_core")
          ENDIF
          IF ( n.LT.3 ) THEN
            WRITE(*,'(a,i0,a)') "**** There's nothing like ",n,"d"
            CALL juDFT_error("storeconfig: invalid specifier",calledby
     +           ="set_atom_core")
          ENDIF
          naddst=naddst+1
          coren(nst+naddst) = n
          corek(nst+naddst) = 2
c
c         Prefer a magnetic configuration to a nonmagnetic one
c         Example: 3 d electrons are distributed as
c         2 electrons in the kappa=2 and 1 electron in the kappa=-3 state
c
          IF ( o <= 2 ) THEN
            coreo(nst+naddst) = REAL(o)
          ELSEIF ( o < 6 ) THEN
            coreo(nst+naddst) = 2.0
          ELSEIF ( o <= 7 ) THEN
            coreo(nst+naddst) = REAL(o) - 3.0
          ELSE
            coreo(nst+naddst) = 4.0
          ENDIF
          o = o - coreo(nst+naddst)
          IF ( o.GT.0 ) THEN
            naddst=naddst+1
            coren(nst+naddst) = n
            corek(nst+naddst) = -3
            coreo(nst+naddst)  = REAL(o)
          ENDIF
! --- f-states
        CASE("f")
          IF ( o.GT.14 ) THEN
            WRITE(*,'(A,I3,A)') "*** ",o,
     A                          " f-electrons are impossible!"
            CALL juDFT_error("storeconfig: too many f-electrons",
     +           calledby="set_atom_core")
          ENDIF
          IF ( n.LT.4 ) THEN
            WRITE(*,'(a,i0,a)') "*** There's nothing like ",n,"f"
            CALL juDFT_error("storeconfig: invalid specifier",calledby
     +           ="set_atom_core")
          ENDIF
          naddst=naddst+1
          coren(nst+naddst) = n
          corek(nst+naddst) = 3
c
c         Prefer a magnetic configuration to a nonmagnetic one
c         Example: 9 f electrons are distributed as
c         5 electrons in the kappa=3 and 4 electrons in the kappa=-4 state
c
          IF ( o <= 3 ) THEN
            coreo(nst+naddst) = REAL(o)
          ELSEIF ( o < 8 ) THEN
            coreo(nst+naddst) = 3.0
          ELSEIF ( o <= 10 ) THEN
            coreo(nst+naddst) = REAL(o) - 4.0
          ELSE
            coreo(nst+naddst) = 6.0
          ENDIF
          o = o - coreo(nst+naddst)
          IF ( o.GT.0 ) THEN
            naddst=naddst+1
            coren(nst+naddst) = n
            corek(nst+naddst) = -4
            coreo(nst+naddst) = REAL(o)
          ENDIF
! --- other states, abort!!
        CASE DEFAULT
          WRITE(*,'(2A)') "*** invalid l-specifier: ",l
          CALL juDFT_error("storeconfig: invalid specifier",calledby
     +         ="set_atom_core")
      END SELECT

      storeconfig=naddst

      END FUNCTION storeconfig

      END MODULE m_setatomcore
