      MODULE m_enpara
      use m_juDFT
!     *************************************************************
!     Module containing three subroutines
!     r_enpara: read enpara file
!     w_enpara: write enpara file
!     mix_enpara: calculate new energy parameters
!     *************************************************************
      CONTAINS
      SUBROUTINE w_enpara(   &
     &                    atoms,jspin,film,&
     &                    enpara,&
     &                    id)
!
! write enpara-file
!
      USE m_types
      IMPLICIT NONE
      
      INTEGER, INTENT (IN) :: jspin,id
      LOGICAL,INTENT(IN)   :: film
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_enpara),INTENT(IN) :: enpara

      INTEGER n,l,lo
      LOGICAL l_opened

      INQUIRE(unit=40,OPENED=l_opened)
      if (.not.l_opened) return


      WRITE (40,FMT=8035) jspin,enpara%enmix(jspin)
      WRITE (40,FMT=8036)
 8035 FORMAT (5x,'energy parameters          for spin ',i1,' mix=',f10.6)
 8036 FORMAT (t6,'atom',t15,'s',t24,'p',t33,'d',t42,'f')
      DO n = 1,atoms%ntype
         WRITE (6,FMT=8040)  n, (enpara%el0(l,n,jspin),l=0,3),&
     &                          (enpara%lchange(l,n,jspin),l=0,3),enpara%skiplo(n,jspin)
         WRITE (id,FMT=8040) n, (enpara%el0(l,n,jspin),l=0,3),&
     &                          (enpara%lchange(l,n,jspin),l=0,3),enpara%skiplo(n,jspin)
         WRITE (40,FMT=8040) n, (enpara%el0(l,n,jspin),l=0,3),&
     &                          (enpara%lchange(l,n,jspin),l=0,3),enpara%skiplo(n,jspin)
!--->    energy parameters for the local orbitals
         IF (atoms%nlo(n).GE.1) THEN
            WRITE (6,FMT=8039) (enpara%ello0(lo,n,jspin),lo=1,atoms%nlo(n))
            WRITE (6,FMT=8038) (enpara%llochg(lo,n,jspin),lo=1,atoms%nlo(n))
            WRITE (id,FMT=8039) (enpara%ello0(lo,n,jspin),lo=1,atoms%nlo(n))
            WRITE (id,FMT=8038) (enpara%llochg(lo,n,jspin),lo=1,atoms%nlo(n))
            WRITE (40,FMT=8039) (enpara%ello0(lo,n,jspin),lo=1,atoms%nlo(n))
            WRITE (40,FMT=8038) (enpara%llochg(lo,n,jspin),lo=1,atoms%nlo(n))
         END IF

      ENDDO
 8038 FORMAT (' --> change   ',60(l1,8x))
 8039 FORMAT (' --> lo ',60f9.5)
 8040 FORMAT (' -->',i3,1x,4f9.5,' change: ',4l1,' skiplo: ',i3)

      IF (film) THEN
         WRITE (40,FMT=8050) enpara%evac0(1,jspin),enpara%lchg_v(1,jspin),enpara%evac0(2,jspin)
         WRITE (6,FMT=8050)  enpara%evac0(1,jspin),enpara%lchg_v(1,jspin),enpara%evac0(2,jspin)
         WRITE (id,FMT=8050) enpara%evac0(1,jspin),enpara%lchg_v(1,jspin),enpara%evac0(2,jspin)
 8050    FORMAT ('  vacuum parameter=',f9.5,' change: ',l1,&
     &           ' second vacuum=',f9.5)
      ENDIF

      RETURN
      END SUBROUTINE w_enpara
!
!------------------------------------------------------------------
      SUBROUTINE r_enpara(&
     &                    atoms,input,jsp,&
     &                    enpara)
!------------------------------------------------------------------
      USE m_types
      IMPLICIT NONE

      INTEGER, INTENT (IN)        :: jsp
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_enpara),INTENT(INOUT):: enpara

      INTEGER n,l,lo,skip_t,io_err
      enpara%lchange(:,:,jsp) =.false.
      enpara%el0(:,:,jsp)     =0.0
      enpara%ello0(:,:,jsp)   =0.0

!-->  first line contains mixing parameter!

      enpara%enmix(jsp) = 0.0
      READ (40,FMT ='(48x,f10.6)',iostat=io_err) enpara%enmix(jsp)
      IF (io_err /= 0) THEN
         !use defaults
         enpara%lchange(:,:,jsp)=.false.
         enpara%llochg(:,:,jsp)=.false.
         enpara%evac0(:,jsp)=-0.1
         enpara%skiplo(:,jsp) = 0
         enpara%enmix(jsp) = 0.0
         enpara%lchg_v(:,jsp)=.false.
         enpara%el0(:,:,jsp)=-99999.
         CALL default_enpara(jsp,atoms,enpara)
         WRITE (6,FMT=8001) jsp
         WRITE (6,FMT = 8000)
         DO n = 1,atoms%ntype
            WRITE (6,FMT=8140) n,(enpara%el0(l,n,jsp),l=0,3),&
     &                          (enpara%lchange(l,n,jsp),l=0,3),enpara%skiplo(n,jsp)    
            IF (atoms%nlo(n)>=1) THEN
               WRITE (6,FMT = 8139)          (enpara%ello0(lo,n,jsp),lo=1,atoms%nlo(n))
               WRITE (6,FMT = 8138)         (enpara%llochg(lo,n,jsp),lo = 1,atoms%nlo(n))
            ENDIF
         ENDDO                   
         RETURN
      ENDIF
      READ (40,*)                       ! skip next line
      IF (enpara%enmix(jsp).EQ.0.0) enpara%enmix(jsp) = 1.0
      WRITE (6,FMT=8001) jsp
      WRITE (6,FMT=8000)
      skip_t = 0
      DO n = 1,atoms%ntype
         READ (40,FMT=8040,END=200) (enpara%el0(l,n,jsp),l=0,3),&
     &                          (enpara%lchange(l,n,jsp),l=0,3),enpara%skiplo(n,jsp)    
         WRITE (6,FMT=8140) n,(enpara%el0(l,n,jsp),l=0,3),&
     &                          (enpara%lchange(l,n,jsp),l=0,3),enpara%skiplo(n,jsp)    
!
!--->    energy parameters for the local orbitals
!
         IF (atoms%nlo(n).GE.1) THEN
             skip_t = skip_t + enpara%skiplo(n,jsp) * atoms%neq(n)
             READ (40,FMT=8039,END=200)  (enpara%ello0(lo,n,jsp),lo=1,atoms%nlo(n))
             READ (40,FMT=8038,END=200) (enpara%llochg(lo,n,jsp),lo=1,atoms%nlo(n))
             WRITE (6,FMT=8139)          (enpara%ello0(lo,n,jsp),lo=1,atoms%nlo(n))
             WRITE (6,FMT=8138)         (enpara%llochg(lo,n,jsp),lo=1,atoms%nlo(n))
         ELSEIF (enpara%skiplo(n,jsp).GT.0) THEN
             WRITE (6,*) "for atom",n," no LO's were specified"
             WRITE (6,*) 'but skiplo was set to',enpara%skiplo 
             CALL juDFT_error("No LO's but skiplo",calledby ="enpara",&
     &        hint="If no LO's are set skiplo must be 0 in enpara")
         END IF
!
!--->    set the energy parameters with l>3 to the value of l=3
!
         DO  l = 4,atoms%lmax(n)
             enpara%el0(l,n,jsp) = enpara%el0(3,n,jsp)
         ENDDO
      ENDDO   ! atoms%ntype
 
      IF (input%film) THEN
         enpara%lchg_v = .true.
         READ (40,FMT=8050,END=200) enpara%evac0(1,jsp),enpara%lchg_v(1,jsp),enpara%evac0(2,jsp)
         WRITE (6,FMT=8150)         enpara%evac0(1,jsp),enpara%lchg_v(1,jsp),enpara%evac0(2,jsp)
      ENDIF
      IF (atoms%nlod.GE.1) THEN               
         WRITE (6,FMT=8090) jsp,skip_t
         WRITE (6,FMT=8091) 
      END IF

! input formats

 8038 FORMAT (14x,60(l1,8x))
 8039 FORMAT (8x,60f9.5)
 8040 FORMAT (8x,4f9.5,9x,4l1,9x,i3)
 8050 FORMAT (19x,f9.5,9x,l1,15x,f9.5)

! output formats

 8138 FORMAT (' --> change   ',60(l1,8x))
 8139 FORMAT (' --> lo ',60f9.5)
 8140 FORMAT (' -->',i3,1x,4f9.5,' change: ',4l1,' skiplo: ',i3)
 8150 FORMAT ('  vacuum parameter=',f9.5,' change: ',l1,&
     &           ' second vacuum=',f9.5)
 8001 FORMAT ('READING enpara for spin: ',i1)
 8000 FORMAT (/,' energy parameters:',/,t10,'s',t20,&
     &        'p',t30,'d',t37,'higher l - - -')
 8090 FORMAT ('Spin: ',i1,' -- ',i3,'eigenvalues')
 8091 FORMAT ('will be skipped for energyparameter computation')

      RETURN

 200  WRITE (6,*) 'the end of the file enpara has been reached while'
      WRITE (6,*) 'reading the energy-parameters.'
      WRITE (6,*) 'possible reason: energy parameters have not been'
      WRITE (6,*) 'specified for all atom types.'
      WRITE (6,FMT='(a,i4)')&
     &     'the actual number of atom-types is: ntype=',atoms%ntype
      CALL juDFT_error&
     &     ("unexpected end of file enpara reached while reading"&
     &     ,calledby ="enpara")

      END SUBROUTINE r_enpara
!
!------------------------------------------------------------------
      SUBROUTINE mix_enpara(&
     &                      jsp,atoms,vacuum,&
     &                      obsolete,input,enpara,&
     &                      vr,vz,pvac,svac,&
     &                      ener,sqal,enerlo,sqlo)
!------------------------------------------------------------------
      USE m_types
      IMPLICIT NONE
      INTEGER,INTENT(IN)             :: jsp
      TYPE(t_atoms),INTENT(IN)       :: atoms
      TYPE(t_vacuum),INTENT(IN)      :: vacuum
      TYPE(t_obsolete),INTENT(IN)    :: obsolete
      TYPE(t_input),INTENT(IN)       :: input
      TYPE(t_enpara),INTENT(INOUT)   :: enpara

      REAL,    INTENT(IN) :: vr(atoms%jmtd,atoms%ntype)
      REAL,    INTENT(IN) :: ener(0:3,atoms%ntype),sqal(0:3,atoms%ntype)
      REAL,    INTENT(IN) :: enerlo(atoms%nlod,atoms%ntype),sqlo(atoms%nlod,atoms%ntype)
      REAL,    INTENT(IN) :: pvac(2),svac(2),vz(vacuum%nmzd,2)

      INTEGER ityp,j,l,lo
      REAl    vbar,maxdist
      INTEGER same(atoms%nlod)
      LOGICAL l_int_enpara

      l_int_enpara=all(enpara%el0==int(enpara%el0)) !test only enpara of lapw, no lo

      maxdist=0.0
      DO ityp = 1,atoms%ntype
!        look for LO's energy parameters equal to the LAPW (and previous LO) ones
         same = 0
         DO lo = 1,atoms%nlo(ityp)
           IF(enpara%el0(atoms%llo(lo,ityp),ityp,jsp).eq.enpara%ello0(lo,ityp,jsp)) same(lo)=-1
           DO l = 1,lo-1
             IF(atoms%llo(l,ityp).ne.atoms%llo(lo,ityp)) cycle
             IF(enpara%ello0(l,ityp,jsp).eq.enpara%ello0(lo,ityp,jsp).and.same(lo).eq.0)&
     &         same(lo)=l
           ENDDO
         ENDDO
!
!--->   change energy parameters
!
         IF ( obsolete%lepr.EQ.1) THEN
            j = atoms%jri(ityp) - (log(4.0)/atoms%dx(ityp)+1.51)
            vbar = vr(j,ityp)/( atoms%rmt(ityp)*exp(atoms%dx(ityp)*(j-atoms%jri(ityp))) )
         ELSE
            vbar = 0.0
         END IF
 
         DO l = 0,3
            IF ( enpara%lchange(l,ityp,jsp) ) THEN
               write(6,*) 'Type:',ityp,' l:',l
               write(6,FMT=777) enpara%el0(l,ityp,jsp),&
     &              (ener(l,ityp)/sqal(l,ityp) - vbar),&
     &              abs(enpara%el0(l,ityp,jsp)-(ener(l,ityp)/sqal(l,ityp) - vbar))
               maxdist=max(maxdist,&
     &              abs(enpara%el0(l,ityp,jsp)-(ener(l,ityp)/sqal(l,ityp) - vbar)))
               enpara%el0(l,ityp,jsp) =(1.0-enpara%enmix(jsp))*enpara%el0(l,ityp,jsp) + &
     &              enpara%enmix(jsp)*(ener(l,ityp)/sqal(l,ityp) - vbar)
           ENDIF
         ENDDO
         DO l = 4, atoms%lmaxd
            IF ( enpara%lchange(3,ityp,jsp) ) THEN
              enpara%el0(l,ityp,jsp) = enpara%el0(3,ityp,jsp)
            ENDIF
         ENDDO
!
!--->    determine and change local orbital energy parameters
!
         DO lo = 1,atoms%nlo(ityp)
            IF (atoms%l_dulo(lo,ityp)) THEN
               enpara%ello0(lo,ityp,jsp) =enpara%el0(atoms%llo(lo,ityp),ityp,jsp)
            ELSEIF (enpara%llochg(lo,ityp,jsp) ) THEN
               IF(same(lo).eq.-1) THEN
                 enpara%ello0(lo,ityp,jsp) = enpara%el0(atoms%llo(lo,ityp),ityp,jsp)
                 cycle
               ELSE IF(same(lo).gt.0) THEN
                 enpara%ello0(lo,ityp,jsp) = enpara%ello0(same(lo),ityp,jsp)
                 cycle
               ENDIF 
               write(6,*) 'Type:',ityp,' lo:',lo
               write(6,FMT=777) enpara%ello0(lo,ityp,jsp),&
     &           (enerlo(lo,ityp)/sqlo(lo,ityp) - vbar),&
     &          abs(enpara%ello0(lo,ityp,jsp)-(enerlo(lo,ityp)/sqlo(lo,ityp)-vbar))
               maxdist=max(maxdist,&
     &         abs(enpara%ello0(lo,ityp,jsp)-(enerlo(lo,ityp)/sqlo(lo,ityp)-vbar)))
               enpara%ello0(lo,ityp,jsp) =(1.0-enpara%enmix(jsp))*enpara%ello0(lo,ityp,jsp)+&
     &              enpara%enmix(jsp)*(enerlo(lo,ityp)/sqlo(lo,ityp) - vbar)
            END IF
         END DO

      END DO

      IF (input%film) THEN

         IF (obsolete%lepr.eq.1) THEN
           vbar = vz(1,1)
         ELSE
           vbar = 0.0
         ENDIF

         IF (enpara%lchg_v(1,jsp) ) THEN
            write(6,*) 'Vacuum:'
            write(6,FMT=777) enpara%evac0(1,jsp),(pvac(1)/svac(1) - vbar),&
     &              abs(enpara%evac0(1,jsp)-(pvac(1)/svac(1) - vbar))
            maxdist=max(maxdist,abs(enpara%evac0(1,jsp)-(pvac(1)/svac(1) - vbar)))
            enpara%evac0(1,jsp) =(1.0-enpara%enmix(jsp))*enpara%evac0(1,jsp)+&
     &              enpara%enmix(jsp)*(pvac(1)/svac(1) - vbar)
            IF (vacuum%nvac.EQ.2) THEN
               IF (obsolete%lepr.eq.1) vbar = vz(1,vacuum%nvac)
               enpara%evac0(2,jsp) = (1.0-enpara%enmix(jsp))*enpara%evac0(2,jsp)+&
     &              enpara%enmix(jsp)*(pvac(2)/svac(2) - vbar)
            ELSE
               enpara%evac0(2,jsp) = enpara%evac0(1,jsp)
            ENDIF
         ENDIF

      ENDIF
      WRITE(6,'(a36,f12.6)') 'Max. mismatch of energy parameters:',&
     &                                                       maxdist
      if (maxdist>1.0.and..not.l_int_enpara) call juDFT_warn&
     &     ("Energy parameter mismatch too large",hint&
     &     ="If any energy parameters calculated from the output "//&
     &     "differ from the input by more than 1Htr, chances are "//&
     &     "high that your initial setup was broken.")
      RETURN
 777  FORMAT('Old:',f8.5,' new:',f8.5,' diff:',f8.5)

      END SUBROUTINE mix_enpara

      subroutine default_enpara(jsp,atoms,enpara)
      !Create default values for enpara given as simple integers
      use m_types
      implicit none
      INTEGER   ,INTENT(IN)       :: jsp
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_enpara),INTENT(INOUT):: enpara

      INTEGER             :: n,i

 
      DO n = 1,atoms%ntype
         IF (all(enpara%el0(:,n,jsp)>-9999.)) cycle !enpara was set already
         IF ( atoms%nz(n) < 3 ) THEN
            enpara%el0(:,n,jsp) = real( (/1,2,3,4/) )
         ELSEIF ( atoms%nz(n) < 11 ) THEN
            enpara%el0(:,n,jsp) = real( (/2,2,3,4/) )
         ELSEIF ( atoms%nz(n) < 19 ) THEN
            enpara%el0(:,n,jsp) = real( (/3,3,3,4/) )
         ELSEIF ( atoms%nz(n) < 31 ) THEN
            enpara%el0(:,n,jsp) = real( (/4,4,3,4/) )
         ELSEIF ( atoms%nz(n) < 37 ) THEN
            enpara%el0(:,n,jsp) = real( (/4,4,4,4/) )
         ELSEIF ( atoms%nz(n) < 49 ) THEN
            enpara%el0(:,n,jsp) = real( (/5,5,4,4/) )
         ELSEIF ( atoms%nz(n) < 55 ) THEN
            enpara%el0(:,n,jsp) = real( (/5,5,5,4/) )
         ELSEIF ( atoms%nz(n) < 72 ) THEN
            enpara%el0(:,n,jsp) = real( (/6,6,5,4/) )
         ELSEIF ( atoms%nz(n) < 81 ) THEN
            enpara%el0(:,n,jsp) = real( (/6,6,5,5/) )
         ELSEIF ( atoms%nz(n) < 87 ) THEN
            enpara%el0(:,n,jsp) = real( (/6,6,6,5/) )
         ELSE
            enpara%el0(:,n,jsp) = real( (/7,7,6,5/) )
         ENDIF
         DO i = 1, atoms%nlo(n)
            enpara%ello0(i,n,jsp) = enpara%el0(atoms%llo(i,n),n,jsp) - 1
         ENDDO
      ENDDO
      END subroutine
      END MODULE m_enpara
