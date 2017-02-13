      MODULE m_stepf
      USE m_juDFT
      CONTAINS
        SUBROUTINE stepf(sym,stars,atoms,oneD, input,cell, vacuum)
          !
          !*********************************************************************
          !     calculates the fourier components of the interstitial step
          !     function for the reciprocal vectors of the star list.
          !           m. weinert  1986
          !*********************************************************************
          !
          !     also set up FFT of U(G) on a (-2G:+2G) grid for convolutions
          !
          !*********************************************************************
          USE m_cfft
          USE m_constants
          USE m_od_cylbes
          USE m_types
          IMPLICIT NONE
          !     ..
          TYPE(t_sym),INTENT(IN)        :: sym
          TYPE(t_stars),INTENT(INOUT)   :: stars
          TYPE(t_atoms),INTENT(IN)      :: atoms
          TYPE(t_oneD),INTENT(IN)       :: oneD
          TYPE(t_input),INTENT(IN)      :: input
          TYPE(t_cell),INTENT(IN)       :: cell
          TYPE(t_vacuum),INTENT(IN)     :: vacuum
          !     ..
          !     .. Local Scalars ..
          COMPLEX c_c,c_phs
          REAL c,dd,gs,th,inv_omtil,r_phs
          REAL g_rmt,g_sqr,help,g_abs,fp_omtil,r_c,gr,gx,gy
          INTEGER i,k,n,n3,na,nn,i1,i2,i3,ic,ifft2d,ifftd,kk
          INTEGER ic1,ic2,ic3,icc,im1,im2,im3,loopstart
          !     ..
          !     .. Local Arrays ..
          COMPLEX sf(stars%ng3)
          REAL g(3),gm(3),fJ
          REAL,    ALLOCATABLE :: bfft(:)
          INTEGER, ALLOCATABLE :: icm(:,:,:)
          !     ..
          !     ..
          !--->    if step function on unit14, then just read it in
          !
          ifftd = 27*stars%mx1*stars%mx2*stars%mx3
          !
          OPEN (14,file='wkf2',form='unformatted',status='unknown')
          REWIND 14
          READ (14,END=10,err=10) n3,n
          IF (n3.NE.stars%ng3) GO TO 10
          IF (n.NE.ifftd) GO TO 10
          READ (14) (stars%ustep(i),i=1,stars%ng3)
          READ (14) (stars%ufft(i),i=0,ifftd-1)
          CLOSE (14)
          RETURN

10        CONTINUE

          IF (input%film) THEN
             dd = vacuum%dvac*cell%area/cell%omtil
             IF (oneD%odd%d1) dd = cell%vol/cell%omtil
          ELSE
             dd = 1.0
          END IF
          !--->    G=0 star
          c = 0.0
          DO  n = 1,atoms%ntype
             c = c + atoms%neq(n)*atoms%volmts(n)/cell%omtil
          ENDDO
          stars%ustep(1) = CMPLX(dd-c,0.0)
          !--->    G(parallel)=0  (for film)
          IF (input%film .AND. .NOT.oneD%odd%d1) THEN
             DO  k = 2,stars%ng3
                IF (stars%ig2(k).EQ.1) THEN
                   th = cell%bmat(3,3)*stars%kv3(3,k)*cell%z1
                   stars%ustep(k) = CMPLX(cell%vol*SIN(th)/th/cell%omtil,0.0)
                ELSE
                   stars%ustep(k) = CMPLX(0.0,0.0)
                END IF
             ENDDO
             !-odim
          ELSEIF (oneD%odd%d1) THEN
             DO k = 2,stars%ng3
                gr = 0.0
                IF (stars%kv3(3,k).EQ.0) THEN
                   kk = stars%ig2(k)
                   gr = stars%sk2(kk)
                   CALL od_cylbes(1,gr*cell%z1,fJ)
                   stars%ustep(k) = CMPLX(2.*dd*fJ/(gr*cell%z1),0.)
                ELSE
                   stars%ustep(k) =CMPLX(0.,0.)
                END IF

             ENDDO
             !+odim
          ELSE
             DO  k = 2,stars%ng3
                stars%ustep(k) = CMPLX(0.0,0.0)
             END DO
          END IF
          !--->    sphere contributions
          na = 0
          DO  n = 1,atoms%ntype
             c = 3.*atoms%volmts(n)/cell%omtil
             !-->     structure factors: loop over equivalent atoms
             na = na + 1
             DO  k = 2,stars%ng3
                th = -tpi_const* DOT_PRODUCT(stars%kv3(:,k),atoms%taual(:,na))
                sf(k) = CMPLX(COS(th),SIN(th))
             END DO
             DO  nn = 2,atoms%neq(n)
                na = na + 1
                DO  k = 2,stars%ng3
                   th = -tpi_const* DOT_PRODUCT(stars%kv3(:,k),atoms%taual(:,na))
                   sf(k) = sf(k) + CMPLX(COS(th),SIN(th))
                END DO
             END DO
             !--->    update step function
             DO  k = 2,stars%ng3
                gs = stars%sk3(k)*atoms%rmt(n)
                stars%ustep(k) = stars%ustep(k) - (c* (SIN(gs)/gs-COS(gs))/ (gs*gs))* sf(k)
             ENDDO
          ENDDO
          !
          ! --> set up stepfunction on fft-grid:
          !
          ALLOCATE (  bfft(0:27*stars%mx1*stars%mx2*stars%mx3-1) )
          im1=CEILING(1.5*stars%mx1); im2=CEILING(1.5*stars%mx2); im3=CEILING(1.5*stars%mx3) 
          ALLOCATE ( icm(-im1:im1,-im2:im2,-im3:im3) )
          icm = 0
          ic=1
          inv_omtil=1.0/cell%omtil
          fp_omtil=  -fpi_const*inv_omtil
          !DO first vector before loop
          stars%ufft(0)=0.0
          bfft(0)=0.0
          DO n=1,atoms%ntype
             stars%ufft(0)=stars%ufft(0)+atoms%neq(n)*atoms%volmts(n)
          ENDDO
          stars%ufft(0)=1.0-stars%ufft(0)*inv_omtil
          loopstart=1
          DO i3=0,3*stars%mx3-1
             gm(3)=REAL(i3)
             IF ( gm(3) > 1.5*stars%mx3 ) gm(3)=gm(3)-3.0*stars%mx3
             DO i2=0,3*stars%mx2-1
                gm(2)=REAL(i2)
                IF ( gm(2) > 1.5*stars%mx2 ) gm(2)=gm(2)-3.0*stars%mx2
                DO i1=loopstart,3*stars%mx1-1
                   loopstart=0 !all further loops start at i1=0
                   gm(1)=REAL(i1)
                   IF ( gm(1) > 1.5*stars%mx1 ) gm(1)=gm(1)-3.0*stars%mx1
                   !
                   !-> use inversion <-> c.c.
                   !
                   ic1 = NINT(gm(1)) ; ic2 = NINT(gm(2)) ; ic3 = NINT(gm(3))
                   IF ( gm(3) < 0.0 ) THEN  ! retreive from table icm()
                      icc = icm(-ic1,-ic2,-ic3)
                      !IF (icc.EQ.0) THEN
                      !  write(*,*) ic1,ic2,ic3,icc
                      !   CALL juDFT_error(" error in stepf! ",calledby="stepf")
                      !ENDIF
                      stars%ufft(ic) = stars%ufft(icc)
                      IF (.NOT.sym%invs) bfft(ic) = - bfft(icc)

                      ic=ic+1
                      CYCLE 
                   ELSE                         ! store number in table icm()
                      icm(ic1,ic2,ic3) = ic
                      IF (ic1 == im1) icm(-ic1,ic2,ic3) = ic
                      IF (ic2 == im2) icm(ic1,-ic2,ic3) = ic
                      IF ((ic1 == im1).AND.(ic2 == im2)) icm(-ic1,-ic2,ic3) = ic
                   ENDIF
                   g=MATMUL(TRANSPOSE(cell%bmat),gm)
                   g_sqr = DOT_PRODUCT(g,g)
                   g_abs = SQRT(g_sqr)
                   help = fp_omtil/g_sqr
                   IF (sym%invs) THEN
                      r_c = 0.0
                      !       Better no OpenMP, huge overhead! Parallel region is located in a
                      !       nested loop and is therefore created more then a billion times. 
                      !           U.Alekseeva 15.10.2015 
!!$OMP  PARALLEL DO PRIVATE(r_phs,nn,th,na,g_rmt,n) DEFAULT(SHARED) REDUCTION(+:r_c)
                      DO n=1,atoms%ntype
                         r_phs = 0.0
                         na=SUM(atoms%neq(:n-1))
                         DO nn=1,atoms%neq(n)
                            th=-tpi_const*DOT_PRODUCT(gm,atoms%taual(:,na+nn))
                            r_phs = r_phs + COS(th)
                         ENDDO
                         g_rmt = g_abs * atoms%rmt(n)
                         r_c=r_c+atoms%rmt(n)*(SIN(g_rmt)/g_rmt-COS(g_rmt))*r_phs
                      ENDDO
!!$OMP END PARALLEL DO
                      stars%ufft(ic) = help * r_c
                   ELSE
                      c_c=CMPLX(0.0,0.0)
!!$OMP  PARALLEL DO PRIVATE(c_phs,nn,th,na,g_rmt,n) DEFAULT(SHARED)REDUCTION(+:c_c)

                      DO n=1,atoms%ntype
                         c_phs = CMPLX(0.0,0.0)
                         na=SUM(atoms%neq(:n-1))
                         DO nn=1,atoms%neq(n)
                            th=-tpi_const*DOT_PRODUCT(gm,atoms%taual(:,na+nn))
                            c_phs = c_phs + EXP(CMPLX(0,th))
                         ENDDO
                         g_rmt = g_abs * atoms%rmt(n)
                         c_c=c_c+atoms%rmt(n)*(SIN(g_rmt)/g_rmt-COS(g_rmt))*c_phs
                      ENDDO
!!$OMP END PARALLEL DO
                      stars%ufft(ic) = help * REAL(c_c)
                      bfft(ic) = help * AIMAG(c_c)
                   ENDIF


                   IF (((i3.EQ.3*stars%mx3/2).OR. (i2.EQ.3*stars%mx2/2)).OR. (i1.EQ.3*stars%mx1/2)) THEN
                      stars%ufft(ic)=0.0 
                      bfft(ic)=0.0 
                   ENDIF
                   !-odim
                   IF (oneD%odd%d1) THEN
                      IF (ic.LT.9*stars%mx1*stars%mx2 .AND. ic.NE.0) THEN
                         gx = (cell%bmat(1,1)*gm(1) + cell%bmat(2,1)*gm(2))
                         gy = (cell%bmat(1,2)*gm(1) + cell%bmat(2,2)*gm(2))
                         gr = SQRT(gx**2 + gy**2)
                         CALL od_cylbes(1,gr*cell%z1,fJ)
                         stars%ufft(ic) = stars%ufft(ic) +2*cell%vol*fJ/(gr*cell%z1*cell%omtil)
                      END IF
                   END IF
                   !+odim
                   ic=ic+1
                ENDDO
             ENDDO
          ENDDO
          !
          ! --> add film-contributions
          !
          IF (input%film .AND. .NOT.oneD%odd%d1) THEN

             ifft2d=9*stars%mx1*stars%mx2
             stars%ufft(0)=stars%ufft(0)+cell%vol*inv_omtil-1.0

             DO i3=1,3*stars%mx3-1
                gm(3)=REAL(i3)
                IF ( gm(3) > 1.5*stars%mx3 ) gm(3)=gm(3)-3.0*stars%mx3
                th=cell%bmat(3,3)*gm(3)*cell%z1
                stars%ufft(i3*ifft2d)=stars%ufft(i3*ifft2d)+cell%vol*inv_omtil*SIN(th)/th
             ENDDO

          ELSEIF (oneD%odd%d1) THEN
             !-odim
             stars%ufft(0) = stars%ufft(0)+cell%vol*inv_omtil-1.0
             !+odim

          ENDIF
          !
          ! --> make fft
          !
          IF (sym%invs) bfft=0.0
          CALL cfft(stars%ufft,bfft,ifftd,3*stars%mx1,3*stars%mx1,+1)
          CALL cfft(stars%ufft,bfft,ifftd,3*stars%mx2,9*stars%mx1*stars%mx2,+1)
          CALL cfft(stars%ufft,bfft,ifftd,3*stars%mx3,ifftd,+1)

          DEALLOCATE ( bfft , icm )

          !--->    store on unit14
          REWIND 14
          WRITE (14) stars%ng3,ifftd
          WRITE (14) (stars%ustep(i),i=1,stars%ng3)
          WRITE (14) (stars%ufft(i),i=0,ifftd-1)

          CLOSE (14)

        END SUBROUTINE stepf
      END MODULE m_stepf
