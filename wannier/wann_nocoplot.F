      MODULE m_wann_nocoplot
      USE m_juDFT
c     ++++++++++++++++++++++++++++++++++++++++++++++++
c     +  If noco=t, transform wave functions within  +
c     +  mu-th MT back to the global frame to plot   +
c     +  to plot WFs in a meaningful way.            +
c     +                                              +
c     ++++++++++++++++++++++++++++++++++++++++++++++++
      CONTAINS
      SUBROUTINE wann_nocoplot(atoms,slice,nnne,amat,bmat
     >                        ,nkpts,odi,film,natd,ntypd,jmtd
     >                        ,ntype,neq,pos,jri,rmsh,alph,beta
     >                        ,nqpts,qss,z1,zatom)
!    *****************************************************
      USE m_types
      USE m_xsf_io
      USE m_wann_gwf_tools, ONLY : get_index_q
      USE m_constants

      IMPLICIT NONE

      TYPE(t_atoms),INTENT(IN):: atoms

      LOGICAL, INTENT(IN) :: slice,film
      INTEGER, INTENT(IN) :: nnne,nkpts,ntypd,jmtd,natd
      INTEGER, INTENT(IN) :: ntype,neq(ntypd)
      INTEGER, INTENT(IN) :: jri(ntypd),nqpts
      REAL,    INTENT(IN) :: amat(3,3),bmat(3,3),qss(3,nqpts)
      REAL,    INTENT(IN) :: rmsh(jmtd,ntypd)
      REAL,    INTENT(IN) :: pos(3,natd),z1,zatom(:)
      REAL,    INTENT(IN) :: alph(ntypd),beta(ntypd)

      TYPE(od_inp), INTENT(in) :: odi

      LOGICAL :: twodim,xsf,cartesian
      INTEGER :: nbmin,nbmax,nplot,nplo,nbn
      INTEGER :: nslibd,nkqpts
      INTEGER :: ix,iy,iz,ikpt,jspin
      INTEGER :: ii1,ii2,ii3
      INTEGER :: i1,i2,i3,iintsp
      INTEGER :: gx,gy,gz,help_kpt
      INTEGER :: na,nt,nq,iqpt
      INTEGER :: grid(3),count_mt,count_int,count_vac
      REAL :: vec1(3),vec2(3),vec3(3),zero(3)
      REAL :: pt(3),point(3),rcc2(3)
      REAL :: s,arg,pi,u_r,u_i
      COMPLEX :: phasfac
      COMPLEX :: U(2,2)
      COMPLEX :: wf_local(2),wf_global,xdnout
      CHARACTER(len=30) :: filename
      CHARACTER(len=20) :: vandername,name1,name2,name3

      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename
    
      intrinsic real,aimag,conjg,exp,cmplx

      pi = pimach()
      nkqpts = nkpts*nqpts

      INQUIRE(file ="plot_inp",exist=twodim)
      IF(.NOT.twodim) THEN
         CALL juDFT_error("Need the plot_inp from RNK generation!"
     >                  ,calledby="wann_nocoplot")
      ENDIF

      !<-- Open the plot_inp file for input
      OPEN (18,file='plot_inp')
      READ(18,'(i2,5x,l1)') nplot,xsf

      IF (nplot.ge.2) 
     &     CALL juDFT_error
     +     ("plots one by one, please, this is not charge density"
     +     ,calledby="wann_nocoplot")

         ! the defaults
         twodim = .TRUE.;cartesian=.TRUE.;grid=(/100,100,100/)
         vec1 = (/0.,0.,0./);vec2=(/0.,0.,0./);vec3=(/0.,0.,0./)
         zero = (/0.,0.,0./);filename="default"
         READ(18,plot)
         IF (twodim.AND.ANY(grid(1:2)<1)) 
     +        CALL juDFT_error("Illegal grid size in plot",calledby
     +        ="wann_nocoplot")
         IF (.NOT.twodim.AND.ANY(grid<1)) 
     +        CALL juDFT_error("Illegal grid size in plot",calledby
     +        ="wann_nocoplot")
         IF (twodim) grid(3) = 1
         !calculate cartesian coordinates if needed
         IF (.NOT.cartesian) THEN
            vec1=matmul(amat,vec1)
            vec2=matmul(amat,vec2)
            vec3=matmul(amat,vec3)
            zero=matmul(amat,zero)
         ENDIF
         IF (filename =="default") WRITE(filename,'(a,i2)') "plot",1
      CLOSE(18)

         ! loop over k-points
         DO ikpt=1,nkqpts
            count_mt=0; count_int=0; count_vac=0
            iqpt = get_index_q(ikpt,nkpts)
            write(*,*)'kq',ikpt,' q',iqpt,' qz',qss(3,iqpt)

            ! open the old RNK.ikpt.jspin files
            DO jspin=1,2  ! local frame axis (inside MT)
               WRITE(vandername,202) ikpt,jspin
               OPEN(400+jspin,file=vandername)
               READ(400+jspin,7)gx,gy,gz,help_kpt,nslibd
            ENDDO 

            IF(.not.xsf) THEN
            ! open the new UNK.ikpt.iintsp files
            DO iintsp=1,2  ! global frame axis (INT and vacuum)              
               WRITE(vandername,201) ikpt,iintsp
               OPEN(500+iintsp,file=vandername,status='unknown')
               WRITE(500+iintsp,7)grid(1),grid(2),grid(3),ikpt,nslibd
            ENDDO
            ENDIF

         ! loop over all bands
               nbmin=1
               nbmax=nslibd
         bands:DO nbn = nbmin,nbmax

         IF (xsf) THEN
           do jspin=1,2
            write (name1,22) ikpt,nbn,jspin
   22       format (i5.5,'.',i3.3,'.real.',i1,'.xsf')
            write (name2,23) ikpt,nbn,jspin
   23       format (i5.5,'.',i3.3,'.imag.',i1,'.xsf')
            write (name3,24) ikpt,nbn,jspin
   24       format (i5.5,'.',i3.3,'.absv.',i1,'.xsf')
            OPEN(600+jspin,file=name1)
            CALL xsf_WRITE_atoms(600+jspin,atoms,film,odi%d1,amat)
            OPEN(602+jspin,file=name2)
            CALL xsf_WRITE_atoms(602+jspin,atoms,film,odi%d1,amat)
            OPEN(604+jspin,file=name3)
            CALL xsf_WRITE_atoms(604+jspin,atoms,film,odi%d1,amat)
            CALL xsf_WRITE_header(600+jspin,twodim,filename,(vec1),
     &       (vec2),(vec3),zero
     $           ,grid)
            CALL xsf_WRITE_header(602+jspin,twodim,filename,(vec1),
     &       (vec2),(vec3),zero
     $           ,grid)
            CALL xsf_WRITE_header(604+jspin,twodim,filename,(vec1),
     &       (vec2),(vec3),zero
     $           ,grid)
            enddo
         ENDIF

         ! loop over real space grid points
         DO iz = 0,grid(3)-1
          DO iy = 0,grid(2)-1
           xloop:DO ix = 0,grid(1)-1
            point = zero+vec1*REAL(ix)/grid(1)+vec2*REAL(iy)
     $                 /grid(2)
            IF (.NOT.twodim) point = point+vec3*REAL(iz)/grid(3)

            ! read old RNK information at grid point
            DO jspin=1,2
               READ(400+jspin,8)u_r,u_i
               wf_local(jspin) = cmplx(u_r,u_i)
            ENDDO

             ! is point in MT?
             ii1 = 3
             ii2 = 3
             ii3 = 3
             IF (film .AND. .NOT.odi%d1) ii3 = 0
             IF (odi%d1) THEN
                ii1 = 0 ; ii2 = 0
             END IF
             DO  i1 = -ii1,ii1
              DO  i2 = -ii2,ii2
               DO  i3 = -ii3,ii3
                pt = point+MATMUL(amat,(/i1,i2,i3/))
                na = 0
                DO nt = 1,ntype
                 DO nq = 1,neq(nt)
                  na   = na + 1
                  s  = SQRT(dot_PRODUCT(pos(:,na)-pt,pos(:,na)-pt))
                  IF (s<rmsh(jri(nt),nt)) THEN
                    count_mt=count_mt+1
                    ! we are inside the MT with alph(nt),beta(nt)
                    ! set up transformation local -> global
                   U(1,1) =  exp(-ImagUnit*alph(nt)/2.)*cos(beta(nt)/2.) 
                   U(1,2) = -exp(-ImagUnit*alph(nt)/2.)*sin(beta(nt)/2.)
                   U(2,1) =  exp( ImagUnit*alph(nt)/2.)*sin(beta(nt)/2.)
                   U(2,2) =  exp( ImagUnit*alph(nt)/2.)*cos(beta(nt)/2.)
                    
                    ! transform wfs to global frame
                    DO iintsp=1,2
                       pt = matmul(bmat,rcc2)/tpi_const
                       arg = -pi*real(2*iintsp-3)*( qss(1,iqpt)*rcc2(1)
     >                                             +qss(2,iqpt)*rcc2(2)
     >                                             +qss(3,iqpt)*rcc2(3))
                       phasfac = cmplx(cos(arg),sin(arg))

                       wf_global= phasfac*(  U(iintsp,1)*wf_local(1)
     >                                     + U(iintsp,2)*wf_local(2) )

                       if(xsf) THEN
                          xdnout=wf_global
                          WRITE(600+iintsp,*) real(xdnout)
                          WRITE(602+iintsp,*) aimag(xdnout)
                          WRITE(604+iintsp,*) real(xdnout*conjg(xdnout))
                       ELSE
                          WRITE(500+iintsp,8)real(wf_global),
     >                                       aimag(wf_global)
                       ENDIF
                    ENDDO

                   CYCLE xloop
                  ENDIF
                 ENDDO
                ENDDO !nt
               ENDDO
              ENDDO
             ENDDO !i1

             ! VACUUM region
             IF (SQRT((pt(1))**2+(pt(2))**2)>=z1)THEN
             count_vac=count_vac+1
             DO iintsp=1,2
                phasfac=cmplx(1.,0.)
                wf_global = phasfac*wf_local(iintsp)

                if(xsf) THEN
                   xdnout=wf_global
                   WRITE(600+iintsp,*) real(xdnout)
                   WRITE(602+iintsp,*) aimag(xdnout)
                   WRITE(604+iintsp,*) real(xdnout*conjg(xdnout))      
                ELSE
                   WRITE(500+iintsp,8)real(wf_global),
     >                  aimag(wf_global)
                ENDIF
             ENDDO
 
             CYCLE xloop
             ENDIF

             ! if we are here, point is in INTERSTITIAL
             ! therefore just copy wfs
             count_int = count_int+1
             DO iintsp=1,2
                phasfac=cmplx(1.,0.)
                wf_global = phasfac*wf_local(iintsp)

                if(xsf) THEN
                   xdnout=wf_global
                   WRITE(600+iintsp,*) real(xdnout)
                   WRITE(602+iintsp,*) aimag(xdnout)
                   WRITE(604+iintsp,*) real(xdnout*conjg(xdnout))      
                ELSE
                   WRITE(500+iintsp,8)real(wf_global),
     >                  aimag(wf_global)
                ENDIF
             ENDDO

            ENDDO xloop
           ENDDO
          ENDDO !z-loop

          IF (xsf) THEN    
             DO iintsp=1,2
              CALL xsf_WRITE_endblock(600+iintsp,twodim)
              CALL xsf_WRITE_endblock(602+iintsp,twodim)
              CALL xsf_WRITE_endblock(604+iintsp,twodim)
              CLOSE(600+iintsp); CLOSE(602+iintsp); CLOSE(604+iintsp)
             ENDDO
          ENDIF

               ENDDO bands

               IF(.not.xsf) THEN
               ! close new UNK.ikpt.iintsp files
               DO iintsp=1,2
                  CLOSE(500+iintsp)
               ENDDO
               ENDIF

               ! delete RNK.ikpt.jspin files
               DO jspin=1,2   
                  IF(xsf) THEN
                     CLOSE(400+jspin)
                  ELSE
                     CLOSE(400+jspin,status='delete')
                  ENDIF
               ENDDO

               write(*,*)count_mt,count_int,count_vac

         ENDDO  ! end ikpt loop

  201 FORMAT ('UNK',i5.5,'.',i1)
  202 FORMAT ('RNK',i5.5,'.',i1)
    7 FORMAT (5i4)
    8 FORMAT (f20.12,1x,f20.12)

      END SUBROUTINE wann_nocoplot
!------------------------------------------
      
      END MODULE m_wann_nocoplot
