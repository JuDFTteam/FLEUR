      MODULE m_olap

      CONTAINS 

c     Calculates plane-wave overlap matrix olap defined by GPT(1:3,1:NGPT).
c     (Muffin-tin spheres are cut out.)
c     olap_pw calculates full overlap matrix

      SUBROUTINE olap_pw(olap,gpt,ngpt,ntype,natd,jmtd,neq,omtil,
     &                   rmt,taual,bmat)
      
      USE m_constants 
      IMPLICIT NONE

c     - scalars -
      INTEGER,INTENT(IN)       :: ngpt,ntype,natd,jmtd
      REAL,INTENT(IN)          :: omtil
c     - arrays -
      INTEGER,INTENT(IN)       :: gpt(3,ngpt),neq(ntype)
      REAL,INTENT(IN)          :: rmt(ntype),taual(3,natd),
     &                            bmat(3,3)

#if ( defined(CPP_INVERSION) )
      REAL,INTENT(OUT)         ::  olap(ngpt,ngpt)
#else
      COMPLEX,INTENT(OUT)      ::  olap(ngpt,ngpt)
#endif
c     - local -
      INTEGER                  :: i,j,itype,icent,ineq
      REAL                     :: g,r,fgr
      COMPLEX,PARAMETER        :: img=(0d0,1d0)
      INTEGER                  :: dg(3)
      
     
      DO i=1,ngpt
        DO j=1,i
          dg        = gpt(:,j)-gpt(:,i)
          g         = gptnorm(dg,bmat)
          olap(i,j) = 0
          IF(g.eq.0) THEN
            DO itype=1,ntype
              r          = rmt(itype)
              olap(i,j) = olap(i,j) - neq(itype) *fpi_const*r**3/3/omtil
            END DO
          ELSE
            icent = 0
            DO itype=1,ntype
              r   = g * rmt(itype)
              fgr = fpi_const* ( sin(r) - r*cos(r) ) /g**3 / omtil
              DO ineq=1,neq(itype)
                icent     = icent+1
                olap(i,j) = olap(i,j) - fgr *
     &                 exp(img*tpi_const*dot_product(dg,taual(:,icent)))
              END DO
            END DO
          END IF
          IF(i.eq.j) olap(i,j) = olap(i,j) + 1
#if ( defined(CPP_INVERSION) )
         olap(j,i) = olap(i,j)
#else
         olap(j,i) = conjg(olap(i,j))
#endif          
        END DO
      END DO

      END SUBROUTINE olap_pw

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

c     olap_pwp  calculates upper triangular part of overlap matrix

      SUBROUTINE olap_pwp(olap,gpt,ngpt,ntype,natd,jmtd,neq,omtil,
     &                    rmt,taual,bmat)
      

      USE m_constants 
      IMPLICIT NONE

c     - scalars -
      INTEGER,INTENT(IN)       :: ngpt,ntype,natd,jmtd
      REAL,INTENT(IN)          :: omtil
c     - arrays -
      INTEGER,INTENT(IN)       :: gpt(3,ngpt),neq(ntype)
      REAL,INTENT(IN)          :: rmt(ntype),taual(3,natd),
     &                            bmat(3,3)

#if ( defined(CPP_INVERSION) )
      REAL,INTENT(OUT)         ::  olap(ngpt*(ngpt+1)/2)
#else
      COMPLEX,INTENT(OUT)      ::  olap(ngpt*(ngpt+1)/2)
#endif
c     - local -
      INTEGER                  :: i,j,k,itype,icent,ineq
      REAL                     :: g,r,fgr
      COMPLEX,PARAMETER        :: img=(0d0,1d0)
      INTEGER                  :: dg(3)

      k  = 0
      DO i=1,ngpt
        DO j=1,i
          k       = k + 1
          dg      = gpt(:,i)-gpt(:,j)
          g       = gptnorm(dg,bmat)
          olap(k) = 0
          IF(g.eq.0) THEN
            DO itype=1,ntype
              r       = rmt(itype)
              olap(k) = olap(k) - neq(itype) * fpi_const*r**3/3 / omtil
            END DO
          ELSE
            icent = 0
            DO itype=1,ntype
              r   = g * rmt(itype)
              fgr = fpi_const* ( sin(r) - r*cos(r) ) /g**3 / omtil
              DO ineq=1,neq(itype)
                icent   = icent+1
                olap(k) = olap(k) - fgr * 
     &               exp( img*tpi_const*dot_product(dg,taual(:,icent)) )
              END DO
            END DO
          END IF
          IF(i.eq.j) olap(k) = olap(k) + 1
        END DO
      END DO

      END SUBROUTINE olap_pwp


      PURE FUNCTION gptnorm(gpt,bmat)
      
        IMPLICIT NONE

        REAL                :: gptnorm
        INTEGER,INTENT(IN)  :: gpt(3)
        REAL,INTENT(IN)     :: bmat(3,3)

        gptnorm = sqrt ( sum ( matmul ( gpt(:),bmat(:,:) ) **2 ) )

      END FUNCTION gptnorm

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   

      SUBROUTINE wfolap_init(olappw,olapmt,gpt,ngpt,
     &                       ntype,lmax,lmaxd,nindx,neq,taual,
     &                       jmtd,jri,dx,rmsh,rmt,bmat,omtil,
     &                       bas1,bas2,maxindx,natd)
      
      USE m_util    ,ONLY: intgrf,intgrf_init

      IMPLICIT NONE
      
c     - scalars -
      INTEGER,INTENT(IN)       :: ntype,lmaxd,jmtd,natd,ngpt,maxindx
      REAL,INTENT(IN)          :: omtil

c     - arrays -
      INTEGER,INTENT(IN)       :: lmax(ntype),neq(ntype),
     &                            nindx(0:lmaxd,ntype),jri(ntype)
      INTEGER,INTENT(IN)       :: gpt(3,ngpt)
      REAL,INTENT(IN)          :: dx(ntype),rmsh(jmtd,ntype),rmt(ntype)
      REAL,INTENT(IN)          :: taual(3,natd),bmat(3,3)
      REAL,INTENT(IN)         ::  bas1(jmtd,maxindx,0:lmaxd,ntype),
     &                            bas2(jmtd,maxindx,0:lmaxd,ntype)
      REAL,INTENT(OUT)         :: olapmt(maxindx,maxindx,0:lmaxd,ntype)
#if ( defined(CPP_INVERSION) )
      REAL,INTENT(OUT)         :: olappw(ngpt,ngpt)
#else
      COMPLEX,INTENT(OUT)      :: olappw(ngpt,ngpt)
#endif

c     - local -
      INTEGER                  :: itype,l,nn,n1,n2

      REAL,ALLOCATABLE         :: gridf(:,:)


      CALL intgrf_init(ntype,jmtd,jri,dx,rmsh,gridf)
      olapmt = 0
      DO itype = 1,ntype
        DO l = 0,lmax(itype)
          nn = nindx(l,itype)
          DO n2 = 1,nn
            DO n1 = 1,nn!n2
              !IF( n1 .gt. 2 .or. n2 .gt. 2) CYCLE
              olapmt(n1,n2,l,itype) = intgrf ( 
     &                            bas1(:,n1,l,itype)*bas1(:,n2,l,itype)
     &                           +bas2(:,n1,l,itype)*bas2(:,n2,l,itype),
     &                            jri,jmtd,rmsh,dx,ntype,itype,gridf)
!               olapmt(n2,n1,l,itype) = olapmt(n1,n2,l,itype)
            END DO
          END DO
        END DO
      END DO

      CALL olap_pw(olappw,gpt,ngpt,ntype,natd,jmtd,neq,omtil,rmt,taual,
     &             bmat)


      END SUBROUTINE wfolap_init
      
      
      
      
      FUNCTION wfolap(cmt1,cpw1,cmt2,cpw2,ngpt1,ngpt2,olappw,olapmt,
     &                ntype,neq,natd,lmax,lmaxd,nindx,maxindx,maxlmindx)
      
      USE m_wrapper

      IMPLICIT NONE


c     - scalars -      
      COMPLEX                :: wfolap
      INTEGER,INTENT(IN)     :: ntype,lmaxd,natd,maxindx,
     &                          ngpt1,ngpt2,maxlmindx
c     - arrays - 
      INTEGER,INTENT(IN)     :: neq(ntype),lmax(ntype),
     &                          nindx(0:lmaxd,ntype)
      COMPLEX,INTENT(IN)     :: cmt1(maxlmindx,natd),
     &                          cmt2(maxlmindx,natd)
#if ( defined(CPP_INVERSION) )
      REAL,INTENT(IN)        :: cpw1(ngpt1)
#else
      COMPLEX,INTENT(IN)     :: cpw1(ngpt1)
#endif
      COMPLEX,INTENT(IN)     :: cpw2(ngpt2)
#if ( defined(CPP_INVERSION) )
      REAL,INTENT(IN)        :: olappw(ngpt1,ngpt2)
#else
      COMPLEX,INTENT(IN)     :: olappw(ngpt1,ngpt2)
#endif
      REAL,INTENT(IN)        :: olapmt(maxindx,maxindx,0:lmaxd,ntype)
c     - local -
      INTEGER                :: itype,ieq,iatom,l,m,lm,nn
      COMPLEX                :: carr(ngpt1),cdum
      REAL                   :: rarr1(ngpt1),rarr2(ngpt1),rdum1,rdum2

      wfolap = 0
      iatom  = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom = iatom + 1
          lm = 0
          DO l = 0,lmax(itype)
            DO m = -l,l
              nn     = nindx(l,itype)
              wfolap = wfolap + 
     &                 dot_product ( cmt1(lm+1:lm+nn,iatom),
     &                               matmul( olapmt(:nn,:nn,l,itype),
     &                                       cmt2(lm+1:lm+nn,iatom)   ))
              lm     = lm + nn
            END DO
          END DO
        END DO
      END DO
      
       
      wfolap = wfolap + dot_product(cpw1,matmul(olappw,cpw2))

!       CALL dgemv('N',ngpt1,ngpt2,1d0,olappw,ngpt1,real(cpw2),1,0d0,rarr1,1)
!       CALL dgemv('N',ngpt1,ngpt2,1d0,olappw,ngpt1,aimag(cpw2),1,0d0,rarr2,1)
! 
!       rdum1 = dotprod(cpw1,rarr1)
!       rdum2 = dotprod(cpw1,rarr2)
!       cdum  = cmplx( rdum1, rdum2 )

!       wfolap = wfolap + cdum

      END FUNCTION wfolap


      FUNCTION wfolap1(cmt1,cpw1,cmt2,cpw2,ngpt1,ngpt2,olappw,olapmt,
     &                ntype,neq,natd,lmax,lmaxd,nindx,maxindx,maxlmindx)
     
      IMPLICIT NONE
      
c     -scalars - 
      COMPLEX                :: wfolap1
      INTEGER,INTENT(IN)     :: ntype,lmaxd,natd,maxindx,
     &                          ngpt1,ngpt2,maxlmindx
c     - arrays -
      INTEGER,INTENT(IN)     :: neq(ntype),lmax(ntype),
     &                          nindx(0:lmaxd,ntype)
      COMPLEX,INTENT(IN)     :: cmt1(maxlmindx,natd),
     &                          cmt2(maxlmindx,natd)
#if ( defined(CPP_INVERSION) )
      REAL,INTENT(IN)        :: cpw1(ngpt1),cpw2(ngpt2)
#else
      COMPLEX,INTENT(IN)     :: cpw1(ngpt1),cpw2(ngpt2)
#endif
#if ( defined(CPP_INVERSION) )
      REAL,INTENT(IN)        :: olappw(ngpt1,ngpt2)
#else
      COMPLEX,INTENT(IN)     :: olappw(ngpt1,ngpt2)
#endif
      REAL,INTENT(IN)        :: olapmt(maxindx,maxindx,0:lmaxd,ntype)
     
c     - local -
      INTEGER                :: itype,ieq,iatom,l,m,lm,nn
      
      wfolap1 = 0
      iatom   = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom = iatom + 1
          lm = 0
          DO l = 0,lmax(itype)
            DO m = -l,l
              nn      = nindx(l,itype)
              wfolap1 = wfolap1 + 
     &                  dot_product ( cmt1(lm+1:lm+nn,iatom),
     &                                matmul( olapmt(:nn,:nn,l,itype),
     &                                        cmt2(lm+1:lm+nn,iatom)  ))
              lm      = lm + nn
            END DO
          END DO
        END DO
      END DO
      
      wfolap1 = wfolap1 + dot_product(cpw1,matmul(olappw,cpw2))
      
      END FUNCTION wfolap1


      FUNCTION wfolap2(cmt1,cpw1,cmt2,cpw2,ngpt1,ngpt2,olappw,olapmt,
     &                ntype,neq,natd,lmax,lmaxd,nindx,maxindx,maxlmindx)
      
      IMPLICIT NONE


c     - scalars -      
      COMPLEX                :: wfolap2
      INTEGER,INTENT(IN)     :: ntype,lmaxd,natd,maxindx,
     &                          ngpt1,ngpt2,maxlmindx
c     - arrays - 
      INTEGER,INTENT(IN)     :: neq(ntype),lmax(ntype),
     &                          nindx(0:lmaxd,ntype)
      COMPLEX,INTENT(IN)     :: cmt1(maxlmindx,natd),
     &                          cmt2(maxlmindx,natd)
! #if ( defined(CPP_INVERSION) )
!       REAL,INTENT(IN)        :: cpw1(ngpt1)
! #else
      COMPLEX,INTENT(IN)     :: cpw1(ngpt1)
! #endif
      COMPLEX,INTENT(IN)     :: cpw2(ngpt2)
#if ( defined(CPP_INVERSION) )
      REAL,INTENT(IN)        :: olappw(ngpt1,ngpt2)
#else
      COMPLEX,INTENT(IN)     :: olappw(ngpt1,ngpt2)
#endif
      REAL,INTENT(IN)        :: olapmt(maxindx,maxindx,0:lmaxd,ntype)
c     - local -
      INTEGER                :: itype,ieq,ic,l,m,lm,nn
      
      wfolap2 = 0
      ic     = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          ic = ic + 1
          lm = 0
          DO l = 0,lmax(itype)
            DO m = -l,l
              nn     = nindx(l,itype)
              wfolap2 = wfolap2 + 
     &                 dot_product ( cmt1(lm+1:lm+nn,ic),
     &                               matmul( olapmt(:nn,:nn,l,itype),
     &                                       cmt2(lm+1:lm+nn,ic)     ) )
              lm     = lm + nn
            END DO
          END DO
        END DO
      END DO
      
      wfolap2 = wfolap2 + dot_product(cpw1,matmul(olappw,cpw2))
      
      END FUNCTION wfolap2

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 
c     Calculates overlap between core and valence wavefunctions

      SUBROUTINE olap_cv(maxindx,nkpt,maxlcutc,maxindxc,ntype,
     &                   lmaxd,lmax,lmaxc,lmaxcd,nindxc,nindx,
     &                   core1,core2,bas1,bas2,cmt,neq,jmtd,neigd,
     &                   maxlmindx,natd,jri,dx,rmsh,gridf)
      
      USE m_util    ,ONLY: intgrf,intgrf_init,chr
      
      IMPLICIT NONE

c     - scalars -
      INTEGER,INTENT(IN)    ::  maxindx,nkpt,maxlcutc,maxindxc,ntype,
     &                          lmaxd,lmaxcd,jmtd,neigd,maxlmindx,natd

c     - arrays -
      INTEGER,INTENT(IN)    ::  neq(ntype)
      INTEGER,INTENT(IN)    ::  lmaxc(ntype),lmax(ntype),jri(ntype)
      INTEGER,INTENT(IN)    ::  nindxc(0:maxlcutc,ntype),
     &                          nindx(0:lmaxd,ntype)
      REAL,INTENT(IN)       ::  dx(ntype),rmsh(jmtd,ntype)
      REAL,INTENT(IN)       ::  core1(jmtd,0:lmaxcd,maxindxc,ntype),
     &                          core2(jmtd,0:lmaxcd,maxindxc,ntype)
      REAL,INTENT(IN)       ::  bas1(jmtd,maxindx,0:lmaxd,ntype),
     &                          bas2(jmtd,maxindx,0:lmaxd,ntype)
      COMPLEX,INTENT(IN)    ::  cmt(neigd,nkpt,maxlmindx,natd)

c     - local scalars -
      INTEGER               :: itype,icent,l,m,lm,i,j

c     - local arrays -
      INTEGER,ALLOCATABLE   :: olapcv_loc(:,:,:,:,:)
      REAL,ALLOCATABLE      :: gridf(:,:)
      REAL,ALLOCATABLE      :: olapcb(:)
      REAL,ALLOCATABLE      :: olapcv_avg(:,:,:,:),olapcv_max(:,:,:,:)  
      COMPLEX,ALLOCATABLE   :: olapcv(:,:)
      CHARACTER, PARAMETER  :: lchar(0:38) =
     &          (/'s','p','d','f','g','h','i','j','k','l','m','n','o',
     &            'x','x','x','x','x','x','x','x','x','x','x','x','x',
     &            'x','x','x','x','x','x','x','x','x','x','x','x','x' /)


c      CALL intgrf_init(ntype,jmtd,jri,dx,rmsh,gridf)
      WRITE(6,'(/A)') 'Overlap <core|basis>'
      ALLOCATE ( olapcb(maxindx),olapcv(neigd,nkpt),
     &     olapcv_avg(  -maxlcutc:maxlcutc,maxindxc,0:maxlcutc,ntype),
     &     olapcv_max(  -maxlcutc:maxlcutc,maxindxc,0:maxlcutc,ntype),
     &     olapcv_loc(2,-maxlcutc:maxlcutc,maxindxc,0:maxlcutc,ntype) )
 
      DO itype=1,ntype
        IF(ntype.gt.1) WRITE(6,'(A,I3)') 'Atom type',itype
        DO l=0,lmaxc(itype)
          IF(l.gt.lmax(itype)) THEN
            WRITE(*,*) 'l greater then lmax(itype)'
            EXIT ! very improbable case
          END IF
c          WRITE(6,8001) (lchar(l),i=1,min(3,nindx(l,itype)))
          DO i=1,nindxc(l,itype)
            WRITE(6,'(I1,A,2X)',advance='no') i+l,lchar(l)
            DO j=1,nindx(l,itype)
              olapcb(j) = intgrf( core1(:,l,i,itype)*bas1(:,j,l,itype)+
     &                            core2(:,l,i,itype)*bas2(:,j,l,itype),
     &                            jri,jmtd,rmsh,dx,ntype,itype,gridf )
              WRITE(6,'(F10.6)',advance='no') olapcb(j)
            ENDDO

            lm    = sum ( (/ (nindx(j,itype)*(2*j+1),j=0,l-1) /) )
            icent = sum(neq(1:itype-1))+1 ! take first of group of equivalent atoms
            DO m=-l,l
              olapcv = 0
              DO j=1,nindx(l,itype)
                lm = lm + 1
                olapcv(:,:) = olapcv(:,:) + olapcb(j)*cmt(:,:,lm,icent)
              END DO
              olapcv_avg(  m,i,l,itype) = sqrt( sum(abs(olapcv(:,:))**2)
     &                                          /nkpt/neigd )
              olapcv_max(  m,i,l,itype) = maxval(abs(olapcv(:,:)))
              olapcv_loc(:,m,i,l,itype) = maxloc(abs(olapcv(:,:)))
            END DO
            WRITE(6,*)
              
          END DO
        END DO
      END DO
      
      WRITE(6,'(/A)') 'Average overlap <core|val>'
  
      DO itype=1,ntype
        IF(ntype.gt.1) write(6,'(A,I3)') 'Atom type',itype
        DO l=0,lmaxc(itype)
          DO i=1,nindxc(l,itype)
            WRITE(6,'(I1,A,2X)',advance='no') i+l,lchar(l)
            WRITE(6,'('//chr(2*l+1)//'F10.6)')olapcv_avg(-l:l,i,l,itype)
          END DO
        END DO
      END DO
  
      WRITE(6,'(/A)') 'Maximum overlap <core|val> at (band/kpoint)'
      DO itype=1,ntype
        IF(ntype.gt.1) write(6,'(A,I3)') 'Atom type',itype
        DO l=0,lmaxc(itype)
          DO i=1,nindxc(l,itype)
            WRITE(6,'(I1,A,2X)',advance='no') i+l,lchar(l)
            WRITE(6,'('//chr(2*l+1)//
     &               '(F10.6,'' ('',I3.3,''/'',I4.3,'')''))')
     &                      (olapcv_max(  m,i,l,itype),
     &                      olapcv_loc(:,m,i,l,itype),m=-l,l)
          END DO
        END DO
      END DO
      
      DEALLOCATE (olapcb,olapcv,olapcv_avg,olapcv_max,olapcv_loc)


      END SUBROUTINE olap_cv

      
      END MODULE m_olap
