      MODULE m_u2f
!     ****************************************************
!     write formatted density or potential onto unit '99'
!     e. wimmer   march 1985
!     ****************************************************
      CONTAINS
      SUBROUTINE u2f(&
     &               stars,input,atoms,sphhar,vacuum,&
     &               cell,sym,l_noco)
!
      USE m_loddop
      USE m_constants,ONLY:namat_const
      USE m_types
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_cell),INTENT(IN)   :: cell
      LOGICAL,INTENT(IN)        :: l_noco
!     ..
!     ..
!     .. Local Scalars ..
      INTEGER i,ivac,j,jsp,k,lh,n,izn,na,urec  
      LOGICAL n_exist
!     ..
!     .. Local Arrays ..
      INTEGER it
      CHARACTER(len=8) dop,iop,name(10)
      COMPLEX, ALLOCATABLE :: fpw(:,:),fzxy(:,:,:,:)
      COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
      REAL,    ALLOCATABLE :: fz(:,:,:),fr(:,:,:,:)
!     ..
!     ..
      !atoms%jmtd = maxval(atoms%jri(:))
      !sphhar%nlhd = maxval(sphhar%nlh(:))

      IF (l_noco) THEN
         OPEN (98,file='rhomat_inp',form='unformatted')
         OPEN (99,file='rhomat_form',form='formatted')
      ELSE
         OPEN (98,file='f_unf',form='unformatted')
         OPEN (99,file='f_form',form='formatted')
      ENDIF
!
      ALLOCATE( fpw(stars%ng3,input%jspins),fzxy(vacuum%nmzxy,stars%ng2-1,2,input%jspins) )
      ALLOCATE( cdom(stars%ng3),cdomvz(vacuum%nmz,2),cdomvxy(vacuum%nmzxy,stars%ng2-1,2) )
      ALLOCATE( fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins),fz(vacuum%nmz,2,input%jspins) )

      CALL loddop(&
     &            stars,vacuum,atoms,sphhar,&
     &            input,sym,&
     &            98,&
     &            it,fr,fpw,fz,fzxy)
      write(*,*) 'after loddop:', l_noco
      IF (l_noco) THEN
         cdom(:) = (0.0,0.0)
         cdomvz(:,:) = (0.0,0.0) ; cdomvxy(:,:,:) = (0.0,0.0)
         READ (98,END=101,ERR=101) (cdom(k),k=1,stars%ng3)
         IF (input%film) THEN
            READ (98) ((cdomvz(j,ivac),j=1,vacuum%nmz),ivac=1,vacuum%nvac)
            READ (98) (((cdomvxy(j,k-1,ivac),j=1,vacuum%nmzxy),k=2,stars%ng2),&
     &                                              ivac=1,vacuum%nvac)
         ENDIF
      ENDIF
         
  101 CLOSE (98)

      WRITE (99,FMT=8000) name
 8000 FORMAT (10a8)
      WRITE (6,FMT=8010) name
 8010 FORMAT (' loddop title:',10a8)
      WRITE (99,FMT=8020) iop,dop,it
 8020 FORMAT (a8,1x,a8,' iteration=',i3)
      DO 60 jsp = 1,input%jspins
         WRITE (99,FMT=8030) jsp
 8030    FORMAT ('spin=',i1)
         WRITE (99,FMT=8040) atoms%ntype
 8040    FORMAT ('ntype=',i3)
         na = 1
         DO 30 n = 1,atoms%ntype
            izn = atoms%zatom(n) + 0.01
            WRITE (99,FMT=8050) namat_const(izn),n,atoms%jri(n),atoms%rmt(n),atoms%dx(n)
 8050       FORMAT (a2,2x,'n=',i3,2x,'jri=',i3,2x,'rmt=',f10.6,2x,'dx=',&
     &             f10.6)
            WRITE (99,FMT=8060) atoms%ntypsy(na),sphhar%nlh(atoms%ntypsy(na))
 8060       FORMAT ('ntypsy=',i2,2x,'nlh=',i3)
            DO 20 lh = 0,sphhar%nlh(atoms%ntypsy(na))
               WRITE (99,FMT=8070) lh
 8070          FORMAT ('lh=',i3)
               WRITE (99,FMT=8080) (fr(i,lh,n,jsp),i=1,atoms%jri(n))
 8080          FORMAT (4e20.13)
   20       CONTINUE
            na = na + atoms%neq(n)
   30    CONTINUE
         WRITE (99,FMT=8090) stars%ng3
 8090    FORMAT ('nq3=',i6)
         IF (sym%invs) THEN
            WRITE (99,FMT=8080) (real(fpw(k,jsp)),k=1,stars%ng3)
         ELSE
            WRITE (99,FMT=8080) (fpw(k,jsp),k=1,stars%ng3)
         END IF
         IF (input%film) THEN
            DO 50 ivac = 1,vacuum%nvac
               WRITE (99,FMT=8100) ivac
 8100          FORMAT ('ivac=',i1)
               WRITE (99,FMT=8110) vacuum%nmz,cell%z1,vacuum%delz
 8110          FORMAT ('nmz=',i3,2x,'z1=',f20.13,2x,'delz=',f8.4)
               WRITE (99,FMT=8080) (fz(i,ivac,jsp),i=1,vacuum%nmz)
               WRITE (99,FMT=8120) stars%ng2,vacuum%nmzxy
 8120          FORMAT ('nq2=',i5,2x,'nmzxy=',i5)
               DO 40 k = 2,stars%ng2
                  IF (sym%invs2) THEN
                     WRITE (99,FMT=8080) (real(fzxy(j,k-1,ivac,jsp)),&
     &                                                    j=1,vacuum%nmzxy)
                  ELSE
                     WRITE (99,FMT=8080) (fzxy(j,k-1,ivac,jsp),j=1,&
     &                 vacuum%nmzxy)
                  END IF
   40          CONTINUE
   50       CONTINUE
         END IF
   60 CONTINUE

      IF (l_noco) THEN
         WRITE (99,8080) (cdom(k),k=1,stars%ng3)
         IF (input%film) THEN
            WRITE (99,8080) ((cdomvz(j,ivac),j=1,vacuum%nmz),ivac=1,vacuum%nvac)
            WRITE (99,8080) (((cdomvxy(j,k-1,ivac),j=1,vacuum%nmzxy),k=2,stars%ng2),&
     &                                              ivac=1,vacuum%nvac)
         ENDIF
      ENDIF
      CLOSE (99)
      DEALLOCATE( fpw,fzxy,cdom,cdomvz,cdomvxy,fr,fz)

      END SUBROUTINE u2f
      END MODULE m_u2f
