      MODULE m_atominput
      use m_juDFT

      INTEGER, PARAMETER  :: dbgfh=6, errfh=6, bfh=93, warnfh=6
      REAL, PARAMETER     :: eps=0.00000001

      CONTAINS
!***********************************************************************
!     reads in the parameters associated with atom types from input
!     file. Part of inp-generator
!***********************************************************************
      SUBROUTINE atom_input(
     >                      infh,ntype,zatom,xl_buffer,buffer,nlod,
     >                      jspins,film,neq,idlist,xmlCoreRefOccs,
     X                      nline,jri,lmax,lnonsph,ncst,rmt,dx,bmu,
     X                     xmlCoreStates,xmlPrintCoreStates,xmlCoreOccs,
     <                      nlo,llo,nel,el0,ello0,evac0 )

      USE m_readrecord
      USE m_setatomcore, ONLY : setatom_bystr, setcore_bystr
      USE m_constants,   ONLY : namat_const

      IMPLICIT NONE

! ... Arguments ...
      INTEGER, INTENT (IN)    :: infh  ! file number of input-file
      INTEGER, INTENT (INOUT) :: nline ! current line in this file
      INTEGER, INTENT (INOUT) :: nel   ! number of valence electrons

      INTEGER, INTENT (IN)     :: ntype,xl_buffer,nlod,jspins
      LOGICAL, INTENT (IN)     :: film
      INTEGER, INTENT (IN)     :: neq(ntype)
      REAL   , INTENT (IN)     :: zatom(ntype),idlist(ntype)
      REAL   , INTENT (IN)     :: xmlCoreRefOccs(29)
      INTEGER, INTENT (INOUT)  :: jri(ntype)      ! mt radial mesh points
      INTEGER, INTENT (INOUT)  :: lmax(ntype)     ! max. l to include for density, overlap etc.
      INTEGER, INTENT (INOUT)  :: lnonsph(ntype)  ! max. l for nonspherical MT-contributions
      INTEGER, INTENT (INOUT)  :: ncst(ntype)     ! # of core levels
      INTEGER, INTENT (INOUT)  :: nlo(ntype)      ! # of local orbitals
      INTEGER, INTENT (INOUT)  :: llo(nlod,ntype) ! l-values of the LO's
      REAL, INTENT (INOUT)     :: rmt(ntype)      ! muffin-tin radius 
      REAL, INTENT (INOUT)     :: dx(ntype)       ! log. spacing
      REAL, INTENT (INOUT)     :: bmu(ntype)      ! magnetic moment
      REAL, INTENT (INOUT)     :: xmlCoreOccs(2,29,ntype)
      LOGICAL, INTENT (INOUT)  :: xmlCoreStates(29,ntype)
      LOGICAL, INTENT (INOUT)  :: xmlPrintCoreStates(29,ntype)
      REAL, INTENT (OUT)    :: el0(0:3,ntype),ello0(nlod,ntype),evac0(2)
      CHARACTER(len=xl_buffer) :: buffer

!===> data
      INTEGER, PARAMETER ::  l_buffer=512   ! maximum length of e-config string
      INTEGER, PARAMETER ::  nwdd=2         ! maximum number of windows
      INTEGER, PARAMETER ::  nstd=31        ! maximum number of core states

!===> Local Variables
      INTEGER :: nbuffer,ios,n,i,j,l,d1,d10,aoff,up,dn
      INTEGER :: lmax0_def,lnonsph0_def,jri0_def,ncst0_def
      INTEGER :: lmax0,lnonsph0,jri0,ncst0,nlod0,llod
      INTEGER :: natomst,ncorest,nvalst,z,nlo0
      INTEGER :: xmlCoreStateNumber
      REAL    :: rmt0_def,dx0_def,bmu0_def
      REAL    :: rmt0,dx0,bmu0,zat0,id
      LOGICAL :: fatalerror, h_atom, h_allatoms
      LOGICAL :: idone(ntype) 
      INTEGER :: lonqn(nlod,ntype),skiplo(ntype),z_int(ntype)
      INTEGER :: coreqn(2,nstd,ntype),lval(nstd,ntype),llo0(nlod)
      REAL    :: nelec(0:nwdd),coreocc(nstd,ntype)
      LOGICAL :: lchange(0:3,ntype),llochg(2,ntype)


      CHARACTER(len= l_buffer) :: econfig0_def,econfig0
      CHARACTER(len= l_buffer) :: econfig(ntype) ! verbose electronic config
      CHARACTER(len=80) :: lo(ntype), lo0
      CHARACTER(len=13) :: fname

      CHARACTER(len=1) :: lotype(0:3)

      DATA lotype /'s','p','d','f'/

!---> initialize some variables

      fatalerror = .false.
      h_atom=.false.;h_allatoms=.false.

      idone(1:ntype) = .false.
      z_int(1:ntype) = NINT(zatom(1:ntype))
      lo = ' ' ; nlod0 = 0 ; nlo = 0 ; llod = 0 ; lonqn = 0
      ncst = 0 ; econfig(1:ntype) = ' '
!
      lmax0_def    = -9999  
      lnonsph0_def = -9999
      rmt0_def     = -9999.9
      dx0_def      = -9999.9
      jri0_def     = -9999
      ncst0_def    = -9999
      econfig0_def = 'NONE'
      bmu0_def     = -9999.9

      WRITE(6,*)
      WRITE(6,'(a50)') '==============================================='
      WRITE(6,'(a50)') '===  modifying atomic input for &(all)atom  ==='
      WRITE(6,'(a50)') '==============================================='
      WRITE(6,*)

!===> continue reading input
      
      nbuffer = len_trim(buffer)

      IF ((buffer(1:9)=='&allatoms') .OR. 
     &    (buffer(1:5)=='&atom') .OR.
!    resetting nbuffer for &qss or &soc interferes with the lapw
!    namelist, therefore, its contributions are also checked for.
!    might interfere with other namelists, too.
!    Klueppelberg Jul 2012
     &    (buffer(1:5)=='&comp') .OR.
     &    (buffer(1:5)=='&exco') .OR.
     &    (buffer(1:5)=='&film') .OR.
     &    (buffer(1:4)=='&kpt') ) THEN
      ELSE
        nbuffer = 0 ! reset, to read in after &qss or &soc 
      ENDIF
      
      loop: DO

      IF (nbuffer == 0) THEN
        DO
          CALL read_record(infh,xl_buffer,nline,nbuffer,buffer,ios)
          IF (ios==1) GOTO 999
          IF (ios == 2)  CALL juDFT_error
     +         ("end of file while reading a record",
     +         calledby ="atom_input")
          IF (buffer(1:1)=='&') EXIT
          CALL err(0)
          fatalerror = .true.
        ENDDO
      ENDIF

!===> allatoms

      IF (buffer(1:9)=='&allatoms') THEN
        IF (h_allatoms) CALL err(1)
        h_allatoms = .true.
        IF (h_atom) then
          WRITE (errfh,*)
          WRITE (errfh,*) 'atom_input: ERROR',
     &     'namelist &allatoms must appear before namelist(s) &atom.'
          WRITE (errfh,*)
          fatalerror = .true.
        ELSE
!--->     read defaults for atom defaults
          CALL read_allatoms(
     >                       l_buffer,
     <                       rmt0_def,dx0_def,jri0_def,lmax0_def,
     <                       lnonsph0_def,ncst0_def,econfig0_def,
     <                       bmu0_def,ios)

          IF (ios.NE.0) GOTO 912
          IF (rmt0_def > -9999.8) THEN
            rmt     = rmt0_def
            WRITE (6,'(a25,f12.6)') 'globally changed rmt to',rmt0_def
          ENDIF
          IF (dx0_def  > -9999.8)   THEN
            dx      = dx0_def
            WRITE (6,'(a25,f12.6)') 'globally changed dx  to',dx0_def
          ENDIF
          IF (jri0_def > -9998  )   THEN
            jri     = jri0_def
            WRITE (6,'(a25,i12)') 'globally changed jri to',jri0_def
          ENDIF
          IF (lmax0_def > -9998 )   THEN
            lmax    = lmax0_def
            WRITE (6,'(a26,i12)') 'globally changed lmax to',
     &                                                       lmax0_def
          ENDIF
          IF (lnonsph0_def > -9998) THEN
            lnonsph = lnonsph0_def
            WRITE (6,'(a28,i12)') 'globally changed lnonsph to ',
     &                                                    lnonsph0_def
          ENDIF
          IF (ncst0_def > -9998 )   THEN
            ncst    = ncst0_def
            WRITE (6,'(a26,i12)') 'globally changed ncst to',
     &                                                       ncst0_def
          ENDIF
          IF (econfig0_def.NE.'NONE') THEN
            econfig = econfig0_def
            WRITE (6,'(a26,a80)') 'globally set econfig to ',
     &                                                    econfig0_def
          ENDIF
          IF (bmu0_def > -9999.8)   THEN
            bmu     = bmu0_def
            WRITE (6,'(a25,f12.6)') 'globally changed bmu to',bmu0_def
          ENDIF
        ENDIF

!===> atom
      ELSEIF (buffer(1:5)=='&atom') THEN
        h_atom=.true.

!--->   set atom defaults
        lmax0    = -9999  
        lnonsph0 = -9999
        rmt0     = -9999.9
        dx0      = -9999.9
        jri0     = -9999
        ncst0    = -9999
        econfig0 = 'NONE'
        bmu0     = -9999.9
        lo0      = ' '

!--->   read namelist
        CALL read_atom(
     >                 l_buffer,lotype,
     <                 id,zat0,rmt0,jri0,dx0,lmax0,lnonsph0,
     <                 ncst0,econfig0,bmu0,lo0,nlod0,llod,ios)
        IF (ios.ne.0) THEN
          CALL err(3)
        ELSE
!--->     put the data into the correct place
          DO n = 1, ntype
            IF (abs( id - idlist(n) ) > 0.001) CYCLE
            IF (idone(n)) then
              WRITE (errfh,*) 'atom_input: ERROR. did that one already'
              fatalerror=.true.
              EXIT
            ELSE
              IF (rmt0 > -9999.8) THEN
                rmt(n)  = rmt0
                WRITE (6,'(a9,i4,2a2,a16,f12.6)') 'for atom ',n,
     &                ' (',namat_const(z_int(n)),') changed rmt to',rmt0
              ENDIF
              IF (dx0 > -9999.8) THEN
                dx(n)  = dx0
                WRITE (6,'(a9,i4,2a2,a16,f12.6)') 'for atom ',n,
     &                ' (',namat_const(z_int(n)),') changed dx  to', dx0
              ENDIF
              IF (jri0 > -9998  ) THEN
                jri(n)  = jri0
                WRITE (6,'(a9,i4,2a2,a16,i12)') 'for atom ',n,
     &                ' (',namat_const(z_int(n)),') changed jri to',jri0
              ENDIF
              IF (lmax0 > -9998  ) THEN
                lmax(n)  = lmax0
                WRITE (6,'(a9,i4,2a2,a17,i12)') 'for atom ',n,
     &              ' (',namat_const(z_int(n)),') changed lmax to',lmax0
              ENDIF
              IF (lnonsph0 > -9998  ) THEN
                lnonsph(n)  = lnonsph0
                WRITE (6,'(a9,i4,2a2,a20,i12)') 'for atom ',n,
     &        ' (',namat_const(z_int(n)),') changed lnonsph to',lnonsph0
              ENDIF
              IF (bmu0 > -9999.8  ) THEN
                bmu(n)  = bmu0
                WRITE (6,'(a9,i4,2a2,a16,f12.6)') 'for atom ',n,
     &              ' (',namat_const(z_int(n)),  ') changed bmu to',bmu0
              ENDIF
              IF (ncst0 > -9998  ) THEN
                ncst(n)  = ncst0
                WRITE (6,'(a9,i4,2a2,a17,i12)') 'for atom ',n,
     &             ' (',namat_const(z_int(n)), ') changed ncst to',ncst0
              ENDIF
! ===> electronic configuration
              IF (econfig0.NE.'NONE') THEN
                 econfig(n) = econfig0
                 WRITE (6,'(a9,i4,2a2,a17,a80)') 'for atom ',n,
     &           ' (',namat_const(z_int(n)),') set econfig to ',econfig0
                 CALL setatom_bystr(
     >                              l_buffer,nwdd,econfig(n),
     <                          natomst,ncorest,nvalst,nelec)
                 WRITE (6,'("   corestates =",i3," with",f6.1,
     &                                 " electrons")')  ncorest,nelec(0)
                 WRITE (6,'("   valence st.=",i3," with",f6.1,
     &                                 " electrons")')   nvalst,nelec(1)
                 IF (nelec(2) /= 0) THEN
                 WRITE (6,'("second window found!")')
                 WRITE (6,'("   valence st.=",i3," with",f6.1,
     &                                 " electrons")')   nvalst,nelec(2)
                 ENDIF
                 IF (nelec(0)+nelec(1)+nelec(2)-zatom(n) > 0.01) THEN
                    CALL juDFT_error
     +                   ("econfig does not fit to this atom type!"
     +                   ,calledby ="atom_input")
                 ENDIF
                 IF (ncst0 > -9998  ) THEN
                   IF (ncorest /= ncst0) THEN
                     WRITE (6,'("  ==> core-states (ncst):",i3,
     &                              " =/= (econfig):",i3)') ncst,ncorest
                     CALL juDFT_error
     +                    ("econfig does not fit to the specified ncst"
     +                    ,calledby ="atom_input")
                   ENDIF
                 ELSE
                   ncst(n) = ncorest
                 ENDIF
              ENDIF
! ===> local orbitals
              IF (lo0 /= ' ') THEN
                WRITE (6,'(6a,i3,7a,i3,3a,a80)')
     &                     "nlod =",nlod0," llod =",llod," : ",lo0
                lo(n)      = lo0
                IF (nlod0 > nlod)  CALL juDFT_error("atom_input: too "
     &                              //"many lo",calledby="atom_input")

                nlo(n) = len_trim(lo(n))/2
                DO i = 1, nlo(n)
                  j = 2*i
                  DO l = 0, 3
                    IF (lo(n)(j:j) == lotype(l)) THEN
                      llo(i,n) = l
                    ENDIF
                  ENDDO
                  j = j - 1
                  READ (lo(n)(j:j),*) lonqn(i,n)
                ENDDO
                WRITE (6,'("   nlo(",i3,") = ",i2," llo = ",8i2)') n,
     &                                 nlo(n),(llo(i,n),i=1,nlo(n))
                WRITE (6,'("   lonqn = ",8i2)') (lonqn(i,n),i=1,nlo(n))
              ENDIF
              idone(n)   = .true.
            ENDIF
          ENDDO
        ENDIF

!===> not an atom related namelist, we are done
      ELSE
        exit loop
      ENDIF

      nbuffer = 0
      ENDDO loop

 999  CONTINUE
      IF (fatalerror) 
     &   CALL juDFT_error("ERROR(S) reading input. Check output for "
     &                  //"details.",calledby="atom_input")

!----------- adjust the core-levels, lo's and the energy parameters ----

      coreqn(1:2,1:nstd,1:ntype) = 0
      coreocc(1:nstd,1:ntype) = -1.0

      nel = 0;  el0 = -9999.9
      ello0 = -9999.9; evac0 = -9999.9
      atoms: DO n = 1, ntype

        CALL setcore_bystr(
     >                      n,nstd,ntype,l_buffer,
     X                      econfig,natomst,ncorest,
     <                      coreqn,coreocc)

        IF ( coreqn(1,1,n) /= 0 ) THEN
          DO i = 1, natomst
            IF (coreqn(2,i,n) < 0) THEN
               lval(i,n) = - coreqn(2,i,n) - 1
            ELSE
               lval(i,n) = coreqn(2,i,n)
            ENDIF
          ENDDO 

           d1  = mod(nint(zatom(n)),10)
           d10 = int( (nint(zatom(n)) + 0.5)/10 )
          aoff = iachar('1')-1
          fname = 'corelevels.'//achar(d10+aoff)//achar(d1+aoff)
          OPEN (27,file=fname,form='formatted')
          write(27,'(i3)') natomst

          WRITE (6,*) '----------'
          DO i = 1, ncorest
            WRITE(6,'("     core :",2i3,f6.1)') 
     &             coreqn(1,i,n),coreqn(2,i,n),coreocc(i,n)
            j = INT(coreocc(i,n) / 2)
            IF (coreocc(i,n) > 2*j) THEN
              j = - coreocc(i,n)
            ENDIF
            write(27,'(4i3)') coreqn(1,i,n),coreqn(2,i,n),j,j
            xmlCoreStateNumber = 0
            SELECT CASE(coreqn(1,i,n))
               CASE (1)
                  IF(coreqn(2,i,n).EQ.-1) xmlCoreStateNumber = 1   !(1s1/2)
               CASE (2)
                  IF(coreqn(2,i,n).EQ.-1) xmlCoreStateNumber = 2   !(2s1/2)
                  IF(coreqn(2,i,n).EQ.1) xmlCoreStateNumber = 3    !(2p1/2)
                  IF(coreqn(2,i,n).EQ.-2) xmlCoreStateNumber = 4   !(2p3/2)
               CASE (3)
                  IF(coreqn(2,i,n).EQ.-1) xmlCoreStateNumber = 5   !(3s1/2)
                  IF(coreqn(2,i,n).EQ.1) xmlCoreStateNumber = 6    !(3p1/2)
                  IF(coreqn(2,i,n).EQ.-2) xmlCoreStateNumber = 7   !(3p3/2)
                  IF(coreqn(2,i,n).EQ.2) xmlCoreStateNumber = 9    !(3d3/2)
                  IF(coreqn(2,i,n).EQ.-3) xmlCoreStateNumber = 10  !(3d5/2)
               CASE (4)
                  IF(coreqn(2,i,n).EQ.-1) xmlCoreStateNumber = 8   !(4s1/2)
                  IF(coreqn(2,i,n).EQ.1) xmlCoreStateNumber = 11   !(4p1/2)
                  IF(coreqn(2,i,n).EQ.-2) xmlCoreStateNumber = 12  !(4p3/2)
                  IF(coreqn(2,i,n).EQ.2) xmlCoreStateNumber = 14   !(4d3/2)
                  IF(coreqn(2,i,n).EQ.-3) xmlCoreStateNumber = 15  !(4d5/2)
                  IF(coreqn(2,i,n).EQ.3) xmlCoreStateNumber = 19   !(4f5/2)
                  IF(coreqn(2,i,n).EQ.-4) xmlCoreStateNumber = 20  !(4f7/2)
               CASE (5)
                  IF(coreqn(2,i,n).EQ.-1) xmlCoreStateNumber = 13  !(5s1/2)
                  IF(coreqn(2,i,n).EQ.1) xmlCoreStateNumber = 16   !(5p1/2)
                  IF(coreqn(2,i,n).EQ.-2) xmlCoreStateNumber = 17  !(5p3/2)
                  IF(coreqn(2,i,n).EQ.2) xmlCoreStateNumber = 21   !(5d3/2)
                  IF(coreqn(2,i,n).EQ.-3) xmlCoreStateNumber = 22  !(5d5/2)
                  IF(coreqn(2,i,n).EQ.3) xmlCoreStateNumber = 26   !(5f5/2)
                  IF(coreqn(2,i,n).EQ.-4) xmlCoreStateNumber = 27  !(5f7/2)
               CASE (6)
                  IF(coreqn(2,i,n).EQ.-1) xmlCoreStateNumber = 18  !(6s1/2)
                  IF(coreqn(2,i,n).EQ.1) xmlCoreStateNumber = 23   !(6p1/2)
                  IF(coreqn(2,i,n).EQ.-2) xmlCoreStateNumber = 24  !(6p3/2)
                  IF(coreqn(2,i,n).EQ.2) xmlCoreStateNumber = 28   !(6d3/2)
                  IF(coreqn(2,i,n).EQ.-3) xmlCoreStateNumber = 29  !(6d5/2)
               CASE (7)
                  IF(coreqn(2,i,n).EQ.-1) xmlCoreStateNumber = 25  !(7s1/2)
            END SELECT
            IF(xmlCoreStateNumber.EQ.0) STOP 'Invalid core state!'
            xmlCoreStates(xmlCoreStateNumber,n) = .TRUE.
            xmlPrintCoreStates(xmlCoreStateNumber,n) = 
     +         coreocc(i,n).NE.xmlCoreRefOccs(xmlCoreStateNumber)
            SELECT CASE(xmlCoreStateNumber)
               CASE (9:10,14:15,19:22,26:29)
                  up = MIN((xmlCoreRefOccs(xmlCoreStateNumber)/2),
     +                          coreocc(i,n))
                  dn = MAX(0.0,coreocc(i,n)-up)
               CASE DEFAULT
                  up = CEILING(coreocc(i,n)/2)
                  dn = FLOOR(coreocc(i,n)/2)
            END SELECT
            xmlCoreOccs(1,xmlCoreStateNumber,n) = up
            xmlCoreOccs(2,xmlCoreStateNumber,n) = dn
          ENDDO
          DO i = ncorest+1, natomst
            WRITE(6,'("  valence :",2i3,f6.1,i4,a1)') 
     &             coreqn(1,i,n),coreqn(2,i,n),coreocc(i,n),
     &                      coreqn(1,i,n),lotype(lval(i,n))
            nel = nel + coreocc(i,n) *neq(n)

c           In d and f shells a magnetic alignment of the spins
c           is preferred in the valence bands
c           Hence the up and down occupation is chosen such that
c           the total spin is maximized

            IF ( abs(coreqn(2,i,n)+0.5) > 2.499 )
     +      THEN
              IF ( coreocc(i,n) > abs(coreqn(2,i,n)) ) THEN
                up = abs(coreqn(2,i,n))
                dn = coreocc(i,n) - abs(coreqn(2,i,n))
              ELSE
                up = coreocc(i,n)
                dn = 0
              END IF

c           in s and p states equal occupation of up and down states

            ELSE
              j = INT(coreocc(i,n) / 2)
              IF (coreocc(i,n) > 2*j) THEN
                j = - coreocc(i,n)
              ENDIF
              up = j
              dn = j
            END IF
            WRITE(27,'(4i3,i4,a1)') coreqn(1,i,n),coreqn(2,i,n),up,dn,
     &                      coreqn(1,i,n),lotype(lval(i,n))
          ENDDO
          WRITE (6,*) '----------'
          CLOSE(27)

          DO i = natomst,1,-1                    ! determine valence states
            IF (el0(lval(i,n),n) < -9999.8) THEN ! not processed already
              el0(lval(i,n),n) = REAL(coreqn(1,i,n))
              IF (i <= ncorest) THEN
                el0(lval(i,n),n) = coreqn(1,i,n) + 1.0 ! was already in the core
              ENDIF
            ENDIF
          ENDDO
          DO j = 0,3
            IF (el0(j,n) < -9999.8) THEN
              el0(j,n) = REAL(j+1)
            ENDIF
          ENDDO

        ELSE  ! determine defauts  as usual

          z = NINT(zatom(n))
          nlo0 = nlo(n)
          llo0 = llo(:,n)
          CALL atom_defaults(
     >                       n,ntype,nlod,z,neq,
     X                       ncst0,nel,nlo,llo)

          IF (ncst(n) == 0) ncst(n) = ncst0
          IF (lonqn(1,n) /= 0) THEN ! already set before
            DO i = 1,nlo(n)                       ! subtract lo-charge
              nel = nel - 2*(2*llo(i,n)+1)*neq(n)
              IF (llo(i,n) == 0) ncst(n) = ncst(n) + 1
              IF (llo(i,n) >  0) ncst(n) = ncst(n) + 2
            ENDDO
            nlo(n) = nlo0                         ! set old values
            llo(:,n) = llo0 
            DO i = 1,nlo(n)                       ! add old lo-charge
              nel = nel + 2*(2*llo(i,n)+1)*neq(n)   
              IF (llo(i,n) == 0) ncst(n) = ncst(n) - 1
              IF (llo(i,n) >  0) ncst(n) = ncst(n) - 2
            ENDDO
          ELSE
             lonqn(1:nlo(n),n) = 0 !LO check below should not be needed
                                   !for default setting of enparas
          ENDIF


        ENDIF

        IF (nlo(n) /= 0) THEN                    ! check for local orbitals
          DO i = 1, nlo(n)
            ello0(i,n) = REAL(lonqn(i,n))
            IF ( lonqn(i,n) == NINT(el0(llo(i,n),n)) ) THEN  ! increase qn
              el0(llo(i,n),n) = el0(llo(i,n),n) + 1          ! in LAPW's by 1
            ENDIF
          ENDDO
        ENDIF
        skiplo(n) = 0
        DO i = 1, nlo(n)
          skiplo(n) = skiplo(n) + (2*llo(i,n) + 1)
        ENDDO

      ENDDO atoms

      DO n = 1, ntype
         z = NINT(zatom(n))
         IF (all(el0(:,n)>-9999.)) cycle !enpara was set already
         IF ( z < 3 ) THEN
            el0(:,n) = real( (/1,2,3,4/) )
         ELSEIF ( z < 11 ) THEN
            el0(:,n) = real( (/2,2,3,4/) )
         ELSEIF ( z < 19 ) THEN
            el0(:,n) = real( (/3,3,3,4/) )
         ELSEIF ( z < 31 ) THEN
            el0(:,n) = real( (/4,4,3,4/) )
         ELSEIF ( z < 37 ) THEN
            el0(:,n) = real( (/4,4,4,4/) )
         ELSEIF ( z < 49 ) THEN
            el0(:,n) = real( (/5,5,4,4/) )
         ELSEIF ( z < 55 ) THEN
            el0(:,n) = real( (/5,5,5,4/) )
         ELSEIF ( z < 72 ) THEN
            el0(:,n) = real( (/6,6,5,4/) )
         ELSEIF ( z < 81 ) THEN
            el0(:,n) = real( (/6,6,5,5/) )
         ELSEIF ( z < 87 ) THEN
            el0(:,n) = real( (/6,6,6,5/) )
         ELSE
            el0(:,n) = real( (/7,7,6,5/) )
         ENDIF
! correct valence charge
         DO i = 1,nlo(n)
            IF (llo(i,n).GT.3) THEN
               nel = nel - 2*(2*llo(i,n)+1)*neq(n)   
               IF (llo(i,n) == 0) ncst(n) = ncst(n) + 1
               IF (llo(i,n) >  0) ncst(n) = ncst(n) + 2
            ELSE IF (ello0(i,n).GE.el0(llo(i,n),n)) THEN
               nel = nel - 2*(2*llo(i,n)+1)*neq(n)   
               IF (llo(i,n) == 0) ncst(n) = ncst(n) + 1
               IF (llo(i,n) >  0) ncst(n) = ncst(n) + 2
            END IF
         ENDDO

         DO i = 1, nlo(n)
            IF (ello0(i,n).EQ.0.0) THEN
               ello0(i,n) = el0(llo(i,n),n) - 1
            END IF
         ENDDO
      ENDDO

      WRITE (6,'("Valence Electrons =",i5)') nel


      RETURN

!===> error handling

 911  CONTINUE
      WRITE (errfh,*) 'atom_input: ERROR reading input. ios  =',ios,
     &               ', line =',nline
      CALL juDFT_error("atom_input: ERROR reading input",
     &                calledby="atom_input")

 912  CONTINUE
      WRITE (errfh,*) 'atom_input: ERROR reading namelist.',
     &               ' ios =',ios,
     &               ' line =',nline
      WRITE (errfh,*) buffer(1:nbuffer)
      WRITE (errfh,*) 'The cause of this error may be ...'
      WRITE (errfh,*) '        a variable not defined in this namelist,'
      WRITE (errfh,*) '        wrong type of data for a variable.'
      CALL juDFT_error("atom_input: ERROR reading input",
     &                calledby="atom_input")

 913  CONTINUE
      WRITE (errfh,*) 'atom_input: ERROR reading record.',
     &               ' ios =',ios,
     &               ' line =',nline
      WRITE (errfh,*) buffer(1:nbuffer)
      CALL juDFT_error("atom_input: ERROR reading input",
     &                calledby="atom_input")

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
!----------------------------------------------------------------
      END SUBROUTINE atom_input
!----------------------------------------------------------------
!================================================================
      SUBROUTINE read_allatoms(
     >                         l_buffer,
     X                         rmt,dx,jri,lmax,lnonsph,ncst,econfig,
     <                         bmu,ios)
!****************************************************************
!     reads in defaults for muffin-tin radius, mesh, etc.
!****************************************************************

      IMPLICIT NONE

      INTEGER, INTENT (IN)    :: l_buffer
      INTEGER, INTENT (INOUT) :: jri     ! mt radial mesh points
      INTEGER, INTENT (INOUT) :: lmax    ! max. l to include for density, overlap etc.
      INTEGER, INTENT (INOUT) :: lnonsph ! max. l for nonspherical MT-contributions
      INTEGER, INTENT (INOUT) :: ncst    ! # of core levels
      INTEGER, INTENT (OUT)   :: ios
      REAL, INTENT (INOUT)    :: rmt, dx ! muffin-tin radius and log. spacing
      REAL, INTENT (INOUT)    :: bmu     ! magnetic moment
      CHARACTER(len=l_buffer) :: econfig ! verbose electronic config

      NAMELIST /allatoms/ rmt,dx,jri,lmax,lnonsph,ncst,econfig,
     &                    bmu

      READ (bfh,allatoms,err=911,end=911,iostat=ios)

 911  CONTINUE
      END SUBROUTINE read_allatoms
!================================================================
      SUBROUTINE read_atom(
     >                     l_buffer,lotype,
     X                     id,z,rmt,jri,dx,lmax,lnonsph,ncst,econfig,
     <                     bmu,lo,nlod,llod,ios )
!***********************************************************************
!     reads in muffin-tin radius, mesh, etc.
!***********************************************************************

      USE m_element, ONLY : z_namat
      IMPLICIT NONE

! ... arguments ...
      INTEGER, INTENT (IN)     :: l_buffer
      REAL, INTENT (OUT)       :: id,z,rmt,dx,bmu
      INTEGER                  :: lmax,lnonsph,ncst,jri,nlod,llod
      CHARACTER(len=l_buffer)  :: econfig
      CHARACTER(len=80)        :: lo
      INTEGER, INTENT (OUT)    :: ios
      CHARACTER(len=1), INTENT (IN) :: lotype(0:3)

! ... internal variables ...
      INTEGER                  :: i,j,k,l,n
      REAL                     :: zz
      CHARACTER(len=2)         :: element
      CHARACTER(len=80)        :: lo1

      CHARACTER(len=2) :: lotype2(0:3)
      DATA lotype2 /'sS','pP','dD','fF'/

      NAMELIST /atom/ id,z,rmt,dx,jri,lmax,lnonsph,ncst,
     &                econfig,bmu,lo,element

      id = -9999.9
      z  = -9999.9
      element = ' '

      READ (bfh,atom,err=911,end=911,iostat=ios)

! -> determine which atom we are concerned with ...

      IF ((z < -9999.8).AND.(element.EQ.' ')) THEN
        WRITE (errfh,*)
     &       'ERROR! No element specified  in namelist atom...'
        WRITE (errfh,*) 'use z=.. or element=.. to define it!'
        ios = 3001
        RETURN
      ENDIF
      IF (id < -9999.8) THEN     ! if no id specified
        zz = REAL(z_namat(element))
        IF (z < 0.00) THEN
          IF (zz > eps) THEN
            id = zz              ! use element name
          ELSE
            id = z               ! or use "z" for id
          ENDIF
        ELSE
          IF (zz > 0 .AND. abs(zz-z)>0.5) THEN
            WRITE (warnfh,*)
            WRITE (warnfh,*) 'atom_input: WARNING! ',
     &       'z and z of specified element differ by more than 0.5. '
            WRITE (warnfh,*) '  z = ', z, 'element =',element
            WRITE (warnfh,*)
          ENDIF
        ENDIF
      ENDIF

!---> order local orbitals determine nlod, llod
      lo1 = adjustl(lo)
      lo = ' '
      i = 0
      n = 0
      DO l = 0, 3
        DO
          j = SCAN(lo1,lotype2(l))  ! search for 's' or 'S', 'p' or 'P' etc.
          IF (j > 0) THEN
            lo1(j:j) = ' '
            n = n + 1
            i = i + 1
            IF (j > 1) THEN
              k = SCAN(lo1(j-1:j-1),'123456') ! determine principal quantum number
              IF (k > 0) THEN
                lo(i:i) = lo1(j-1:j-1)
                lo1(j-1:j-1) = ' '
              ELSE
                lo(i:i) = '0'
              ENDIf
            ELSE
              lo(i:i) = '0'
            ENDIF
            i = i + 1
            lo(i:i) = lotype(l)
            nlod = max( nlod, n )
            llod = max( llod, l )
          ELSE
            EXIT
          ENDIF
        ENDDO
      ENDDO
      IF (len_trim(lo1) > 0) then
        WRITE (errfh,*) 'ERROR reading local orbital input...',lo1
        ios = 3002
      ENDIF

 911  CONTINUE
      END SUBROUTINE read_atom
!================================================================

      SUBROUTINE atom_defaults(
     >                         n,ntype,nlod,z,neq,
     X                         ncst2,nel,nlo,llo)

      IMPLICIT NONE

      INTEGER, INTENT (IN)    :: n,ntype,nlod,z
      INTEGER, INTENT (IN)    :: neq(ntype)
      INTEGER, INTENT (INOUT) :: nel,ncst2
      INTEGER, INTENT (INOUT) :: nlo(ntype),llo(nlod,ntype)
      

      INTEGER locore,i
      INTEGER ncst1(0:103),nce(0:24),nval(0:3)
!
! number of core levels for each element; the INT(ncst1/100) number
! provides information about possible local orbitals: 1...(s,p)-LO
! 2...p-LO and 3...d-LO
!
      ncst1 =(/0,0,                                                0,  ! Va,H,He
     +     01, 01,                                  1, 1, 1, 1, 1, 1,  ! Li - Ne
     +     04, 04,                                  4, 4, 4, 4, 4, 4,  ! Na - Ar
     +    107,107,207,207, 7, 7, 7, 7, 7, 7, 7, 7,309, 9, 9, 9, 9, 9,  ! K - Kr
     +    112,112,212,212,12,12,12,12,12,12,12,12,314,14,14,14,14,14,  ! Rb - Xe
     +    117,117,217,217,17,17,17,17,17,17,17,17, 17,17,17,17,17,     ! Cs - Lu
     +                219,19,19,19,19,19,19,19,19,321,21,21,21,21,21,  ! Hf - Rn
     +    124,124,224,224,224,24,24,24,24,24,24,24, 24,24,24,24,24/)   ! Fr - Lw
!
! electrons associated with a given number of core-levels
!
      nce(0) = 0   ; nce(1) = 2   ; nce(4) = 10  ; nce(7) = 18
      nce(9) = 28  ; nce(12) = 36 ; nce(14) = 46 ; nce(17) = 54
      nce(19) = 68 ; nce(21) = 78 ; nce(24) = 86

!
!--> determine core levels
!
      ncst2 = ncst1( z )
      IF (ncst2.GT.300) THEN     ! should add d-LO
        ncst2 = ncst2 - 300 ; locore = 10
        nlo(n) = 1 ; llo(1,n) = 2
      ELSEIF (ncst2.GT.200) THEN ! should add p-LO
        ncst2 = ncst2 - 200 ; locore = 6
        nlo(n) = 1 ; llo(1,n) = 1
      ELSEIF (ncst2.GT.100) THEN ! should add (s,p)-LO
        ncst2 = ncst2 - 100 ; locore = 8
        nlo(n) = 2 ; llo(1,n) = 0 ; llo(2,n) = 1
      ELSE
        nlo(n) = 0 ; locore = 0
      ENDIF
      nel = nel + ( z - nce(ncst2) + locore ) * neq(n)
      IF ((locore == 6).OR.(locore == 10)) ncst2 = ncst2 - 2
      IF (locore == 8) ncst2 = ncst2 - 3

!     WRITE (6,9070) z,ncst2
!9070 FORMAT (i3,3i5)
!     WRITE (6,9090) neq(n),.false.,nlo(n),(llo(i,n),i=1,nlo(n))
!9090 FORMAT (i2,',force =',l1,',nlo=',i2,',llo=',20i2)
!
!--> determine valence states
!


      END SUBROUTINE atom_defaults

      END MODULE m_atominput
