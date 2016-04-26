MODULE m_od_vvacis
CONTAINS
  SUBROUTINE od_vvacis(&
       &     n2d_1,DIMENSION,vacuum,nq2_1,&
       &     kv2_1,cell,MM,stars,&
       &     nstr2_1,&
       &     oneD,&
       &     rht,rhtxy,psq,vz,sym,&
       &     vxy,vpw)

    !     generates m^2+gz^2.neq.0 coefficients of the vacuum potential
    !     and fourier coefficients of the interstitial potential
    !     in the case of 1-dimensional calculations
    !     based on the using of the green's function for the 
    !     radial equation for the potential:
    !     [1/r(d/dr(r*d/dr))-m^2/r^2-gz^2]V(r) = 4pi*rho(r)
    !     when gz.neq.0 green's function is represented in
    !     the following way: 
    !          gf = 4pi*I_m(gz*r<)*K_m(gz*r>)
    !     where I and K are modified bessel funcions
    !     in the case of gz.eq.0 the green's function is:
    !            gf = 2pi*((r</r>)^m)/m
    !                 Y.Mokrousov, autumn 2002

    !     Fully symmetrized version, Mai 2003

    USE m_constants
    USE m_od_cylbes
    USE m_modcyli
    USE m_modcylk
    USE m_vacp5_0
    USE m_vacp5_z
    USE m_visp5_0
    USE m_visp5_z
    USE m_angle
    USE m_qsf
    USE m_fft2d


    USE m_types
    IMPLICIT NONE

    TYPE(t_dimension),INTENT(IN)   :: DIMENSION

    TYPE(t_oneD),INTENT(IN)   :: oneD

    TYPE(t_vacuum),INTENT(IN)   :: vacuum

    TYPE(t_sym),INTENT(IN)   :: sym

    TYPE(t_stars),INTENT(IN)   :: stars

    TYPE(t_cell),INTENT(IN)   :: cell

    INTEGER, INTENT (IN) :: n2d_1  ,MM 
    INTEGER, INTENT (IN) :: nq2_1  
    INTEGER, INTENT (IN) :: nstr2_1(n2d_1)
    INTEGER, INTENT (IN) :: kv2_1(2,n2d_1) 
    COMPLEX, INTENT (INOUT) :: psq(stars%n3d)
    REAL,    INTENT (IN) :: vz(vacuum%nmzd,2,DIMENSION%jspd) 
    REAL,    INTENT (IN) :: rht(vacuum%nmzd,2,DIMENSION%jspd)
    COMPLEX, INTENT (IN) :: rhtxy(vacuum%nmzxyd,n2d_1-1,2,DIMENSION%jspd)
    COMPLEX, INTENT (OUT):: vxy(vacuum%nmzxyd,n2d_1-1,2,DIMENSION%jspd)
    COMPLEX, INTENT (OUT):: vpw(stars%n3d,DIMENSION%jspd)

    !     local
    INTEGER :: m
    !------> reciprocal vectors staff

    INTEGER irec2,irec3,k1,k2,gzi,gzi1,k3,gzmin
    REAL    gz,gx,gy,g,g2 

    !------> different garbage

    REAL, PARAMETER :: tol_21 = 1.0e-21
    INTEGER  imz,imz1,i,ivac,iirec1,m1,j,irc1,im,irc
    INTEGER ix,iy,iz,s,n
    REAL    x,y,r,z,b,q,zf
    REAL    mult
    REAL    ani1,ani2,rhti
    INTEGER ivfft1,ivfft2,ivfft2d
    COMPLEX ic

    !------> different staff with special functions

    REAL  IR(1:stars%k3d,0:MM)
    REAL, ALLOCATABLE :: II(:),KK(:)
    REAL  fJ,fJ2,fJ1,iJ2,IIIR,fJ3
    REAL, ALLOCATABLE :: III(:),IIII(:,:,:)
    REAL, ALLOCATABLE :: fJJ(:,:),iJJ(:,:)

    !------> used for the optimization 

    INTEGER irec1(4),l,lmin ,l1

    !------> some factors used for the construction of the ch. density

    REAL, ALLOCATABLE :: fact(:)
    COMPLEX aa,a

    !------> values of the potential on the vacuum boundary

    COMPLEX val_help
    COMPLEX, ALLOCATABLE :: val(:),val_m(:,:)

    !------> Charge density

    COMPLEX, ALLOCATABLE :: rxy(:)

    !---------------> POTENTIALS

    !-> interstitial potential by its gz,m - components and total on the 
    !-> real grid 

    COMPLEX, ALLOCATABLE :: vis(:,:,:),vis_tot(:,:)

    !-> potential in the vacuum caused by the vacuum charge density      

    COMPLEX, ALLOCATABLE :: pvac(:)

    !-> potential in the vacuum caused by the interst. density

    COMPLEX, ALLOCATABLE :: pint(:)

    !-> radial components of the potential caused by the vacuum  
    !-> density on the real grid  

    COMPLEX, ALLOCATABLE :: vis_help(:,:),vpw_help(:)

    !-> real grid vis_z,vis_0 for the fft2d

    REAL,ALLOCATABLE :: af2(:),bf2(:)

    !-> radial grids

    INTEGER, ALLOCATABLE :: rmap(:,:)
    REAL   , ALLOCATABLE :: rr(:)
    INTEGER nrd

    !--> time

    REAL    gxy0,fxy0,phi
    COMPLEX gxy(stars%n2d-1)
    COMPLEX fxy(stars%n2d-1)

    INTRINSIC REAL,aimag

    !--------- preparations ---------->

    ALLOCATE ( KK(vacuum%nmzxyd),II(vacuum%nmzxyd),III(9*stars%k1d*stars%k2d),&
         &     IIII(9*stars%k1d*stars%k3d,1:stars%k3d,0:MM),&
         &     rmap(0:3*stars%k1d-1,0:3*stars%k2d-1),rr(1:9*stars%k1d*stars%k2d),&
         &     fact(vacuum%nmzxyd),val(n2d_1),fJJ(0:MM+1,stars%n2d),&
         &     val_m(-stars%k3d:stars%k3d,-MM:MM),rxy(vacuum%nmzxyd),iJJ(0:MM+1,1:stars%k3d),&
         &     vis(0:3*stars%k1d-1,0:3*stars%k2d-1,n2d_1),&
         &     vis_tot(0:3*stars%k1d-1,0:3*stars%k2d-1),pvac(vacuum%nmzxyd),&
         &     pint(vacuum%nmzxyd),vis_help(0:3*stars%k1d-1,0:3*stars%k2d-1),&
         &     af2(0:9*stars%k1d*stars%k2d-1),bf2(0:9*stars%k1d*stars%k2d-1),vpw_help(stars%n3d) )

    ivfft2d = 9*stars%k1d*stars%k2d
    ic = CMPLX(0.,1.)
    ivfft1 = 3*stars%k1d
    ivfft2 = 3*stars%k2d
    ani1 = 1./REAL(ivfft1)
    ani2 = 1./REAL(ivfft2)

    !--------- initializations -------->

    !----> vpw in the '1st aproximation' (V - tilde)

    vpw(1,1) = CMPLX(0.,0.)

    DO irec3 = 2,stars%ng3

       g = stars%sk3(irec3)

       vpw(irec3,1) = fpi_const*psq(irec3)/(g*g)

    ENDDO

    DO irc1 = 2,nq2_1
       DO i = 1,vacuum%nmzxy
          vxy(i,irc1-1,1,1) = CMPLX(0.,0.)
       END DO
    END DO

    !----> values of the potential in the 1st approximation on the boundary
    !----> if nstr2.ne.1 then it should be changed!!!

    DO m = 0,MM+1
       DO k3 = 1,stars%k3d
          CALL modcyli(m,cell%bmat(3,3)*k3*cell%z1,iJJ(m,k3))
       END DO
       DO irec2 = 1,stars%ng2
          g2 = stars%sk2(irec2)
          IF (irec2.NE.0) THEN
             CALL od_cylbes(m,g2*cell%z1,fJJ(m,irec2))
          END IF
       END DO
    END DO

    DO irc1 = 1,nq2_1
       val(irc1) = CMPLX(0.,0.)
       m = kv2_1(2,irc1)
       IF (m.LT.0) THEN
          mult = REAL((-1)**m)
       ELSE
          mult = 1.
       END IF
       k3 = kv2_1(1,irc1)
       DO irec2 = 1,stars%ng2
          phi = stars%phi2(irec2)
          irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),k3)
          IF (irec3.NE.0) THEN
             val(irc1) = val(irc1) +&
                  &              (ic**m)*vpw(irec3,1)*EXP(-ic*&
                  &              m*phi)*fJJ(iabs(m),irec2)*&
                  &              stars%nstr2(irec2)*mult
          END IF
       END DO
    END DO

    !-----> preparing arrays for radial grids:selecting from
    !-----> x & y only the radius

    nrd = 0

    DO ix = 0,ivfft1 - 1
       iy_loop: DO  iy = 0,ivfft2 - 1
          x = ix*ani1
          IF (x.GT.0.5) x = x - 1.
          y = iy*ani2
          IF (y.GT.0.5) y = y - 1.
          r = SQRT((x*cell%amat(1,1) + y*cell%amat(1,2))**2 +&
               &               (x*cell%amat(2,1) + y*cell%amat(2,2))**2)
          DO i = 1,nrd
             IF (ABS(r-rr(i)).LE.1.e-6) THEN
                rmap(ix,iy) = i
                CYCLE iy_loop
             END IF
          END DO
          nrd = nrd + 1
          rmap(ix,iy) = nrd
          rr(nrd) = r
       ENDDO iy_loop
    END DO

    DO gzi = -stars%k3d,stars%k3d
       DO m = -MM,MM
          IF (m.LT.0) THEN
             mult = REAL((-1)**m)
          ELSE
             mult = 1.
          END IF
          val_m(gzi,m) = CMPLX(0.,0.)
          irc1 = oneD%ig1(gzi,m)
          IF (irc1.NE.0) THEN
             val_m(gzi,m) = val(irc1)
          ELSE
             DO irec2 = 1,stars%ng2
                phi = stars%phi2(irec2)
                irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),gzi)
                IF (irec3.NE.0) THEN
                   val_m(gzi,m) = val_m(gzi,m) +&
                        &                    (ic**m)*vpw(irec3,1)*EXP(-ic*&
                        &                    m*phi)*fJJ(iabs(m),irec2)*&
                        &                    stars%nstr2(irec2)*mult
                END IF
             END DO
          END IF
       END DO
    END DO

    !-------  During the following staff we miss the m=0,gz=0 -> irec1 = 1
    !-------  component of the potential which should be also added
    !-------  in order to get the total potential for the FFT's

    DO ix = 0,ivfft1 - 1
       DO iy = 0,ivfft2 - 1
          irc = rmap(ix,iy)
          r = rr(irc)
          IF (r.GT.cell%z1) THEN
             zf = (r-cell%z1)/vacuum%delz + 1.0
             im = zf
             q = zf - im
             vis(ix,iy,1) = 0.5*(q-1.)*&
                  &              (q-2.)*vz(im,1,1) -&
                  &              q*(q-2.)*vz(im+1,1,1) +&
                  &              0.5*q*(q-1.)*vz(im+2,1,1)
          ELSE
             vis(ix,iy,1) = &
                  &              vz(1,1,1) - val(1) + tpi_const*&
                  &              psq(1)*(cell%z1*cell%z1 - r*r)/2.
          END IF
          DO irc1 = 2,nq2_1
             vis(ix,iy,irc1) = CMPLX(0.,0.)
          END DO
       END DO
    END DO

    DO m = 0,MM
       DO k3 = 1,stars%k3d
          CALL modcyli(m,cell%bmat(3,3)*k3*cell%z1,IR(k3,m))
       END DO
    END DO

    DO i = 1,nrd
       DO m = 0,MM
          DO k3 = 1,stars%k3d
             CALL modcyli(m,cell%bmat(3,3)*k3*rr(i),IIII(i,k3,m))
          END DO
       END DO
    END DO

    !------- cycle by positive gz---------->

    DO gzi = 0,stars%k3d                        ! gz

       !------- cycle by positive m ---------->

       m_loop: DO m = 0,MM

          !-------------------------------------->

          IF (m.NE.0 .OR. gzi.NE.0) THEN ! m^2 + gz^2.ne.0

             irec1(1) = oneD%ig1(gzi,m)
             irec1(2) = oneD%ig1(gzi,-m)
             irec1(3) = oneD%ig1(-gzi,m)
             irec1(4) = oneD%ig1(-gzi,-m)

             DO l = 1,3
                DO l1 = l+1,4
                   IF (irec1(l).EQ.irec1(l1)) irec1(l1) = 0
                END DO
             END DO

             !---> if all the irec1 not equal to zero

             s = 0

             DO l = 1,4
                s = s + irec1(l)*irec1(l)
             END DO

             !---> preparing special functions which depend only on the 
             !---> absolute value of m and gz

             IF (s.NE.0 .AND. gzi.NE.0) THEN
                gz = cell%bmat(3,3)*gzi
                DO  imz = 1,vacuum%nmzxy
                   z = cell%z1 + vacuum%delz*(imz-1)  
                   CALL modcylk(m,gz*z,KK(imz))
                   CALL modcyli(m,gz*z,II(imz))
                END DO
                IIIR = II(1)
                DO irc = 1,nrd
                   III(irc) = IIII(irc,gzi,m)
                END DO
             ELSEIF (s.EQ.0) THEN
                CYCLE m_loop
             ELSEIF (s.NE.0 .AND. gzi.EQ.0) THEN  
                gz = 0.
             END IF

             !---> now we start the cycle by +-m,gz

             DO  l = 1,4

                IF (irec1(l).NE.0) THEN

                   !--------------------------------------->   

                   DO ix = 0,ivfft1-1
                      DO iy = 0,ivfft2-1
                         vis_help(ix,iy) = CMPLX(0.,0.)
                      END DO
                   END DO

                   DO i = 1,vacuum%nmzxy
                      fact(i) = 0.
                      rxy(i) = CMPLX(0.,0.)
                      pvac(i) = CMPLX(0.,0.)
                      pint(i) = CMPLX(0.,0.)
                   END DO

                   !-------------- gz = 0 ------------------------------------->
                   !----------------------------------------------------------->

                   IF (gzi.EQ.0) THEN

                      !----- this form of the density is just more easy to use

                      DO imz = 1,vacuum%nmzxy
                         rxy(imz) = rhtxy(imz,irec1(l)-1,1,1)
                      END DO

                      !----- vacuum potential caused by the vacuum density

                      CALL vacp5_0(&
                           &        vacuum%nmzxyd,vacuum%nmzxy,cell%z1,tpi_const,rxy,m,vacuum%delz,&
                           &        pvac,fact)

                      !----- vacuum potential caused by the interstitial density

                      aa = CMPLX(0.,0.)
                      m1 = kv2_1(2,irec1(l))
                      gzi1 = kv2_1(1,irec1(l))

                      IF (m1.LT.0) THEN
                         mult = REAL((-1)**m)
                      ELSE
                         mult = 1.
                      END IF

                      DO irec2 = 1,stars%ng2
                         irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),0)
                         IF (irec3.NE.0 .AND. irec2.NE.1) THEN
                            phi = stars%phi2(irec2)
                            g2 = stars%sk2(irec2)
                            aa = aa +    &
                                 &              (ic**m1)*psq(irec3)*&
                                 &              EXP(CMPLX(0.,-m1*phi))*&
                                 &              CMPLX(mult*fJJ(m+1,irec2)/g2,0.)*stars%nstr2(irec2)
                         END IF
                      END DO

                      !----- total vacuum potential

                      pint(:vacuum%nmzxy) =  fact(:vacuum%nmzxy)*aa 

                      vxy(:vacuum%nmzxy,irec1(l)-1,1,1) = pvac(:vacuum%nmzxy) + pint(:vacuum%nmzxy)

                      !----- array val further is a boundary values of the
                      !----- potential V- \tilde \tilde which is created to compensate 
                      !----- the influence of the outer 'noice' charge density - which 
                      !----- is just a periodical continuation of the interstitial charge
                      !----- density, V - \tilde and V - \tilde\tilde are then added in
                      !----- order to obtain the real interstitial potential           

                      val_help = vxy(1,irec1(l)-1,1,1) - val(irec1(l))

                      !----- potential \tilde\tilde{V} is a solution of the Laplase equation
                      !----- in the interstitial with the boundary conditions val_0 and val_z
                      !----- further, it is generated on the uniform grid in the unit cell
                      !----- \tilde{\Omega}, in the space between the cylindrical 
                      !----- interstitial boundary and the squre boundaries it is put to
                      !----- the vacuum potential

                      CALL visp5_0(&
                           &        vacuum%nmzxyd,vacuum%nmzxy,vacuum%delz,m,ivfft1,ivfft2,l,&
                           &        rxy,ani1,ani2,cell%z1,cell%amat,&
                           &        pvac,pint,tpi_const,DIMENSION%jspd,val_help,&
                           &        vis_help)

                      DO ix = 0,ivfft1 - 1
                         DO iy = 0,ivfft2 - 1
                            vis(ix,iy,irec1(l)) = vis_help(ix,iy)
                         END DO
                      END DO


                      !------- gz.NEQ.0--------------------------------------------->
                   ELSE   
                      !------------------------------------------------------------->

                      DO  imz = 1,vacuum%nmzxy
                         rxy(imz) = rhtxy(imz,irec1(l)-1,1,1)
                      END DO

                      !----- vacuum potential caused by the vacuum density        

                      CALL vacp5_z(&
                           &        vacuum%nmzxyd,vacuum%nmzxyd,cell%z1,vacuum%delz,fpi_const,II,KK,rxy,m,&
                           &        pvac)

                      !----- vacuum potential caused by the intst. density

                      a = CMPLX(0.,0.)
                      m1 = kv2_1(2,irec1(l))
                      gzi1 = kv2_1(1,irec1(l))

                      IF (m1.LT.0) THEN
                         mult = REAL((-1)**m)
                      ELSE
                         mult = 1.
                      END IF

                      DO irec2 = 1,stars%ng2
                         irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),gzi1)
                         IF (irec3.NE.0) THEN
                            g = stars%sk3(irec3)
                            g2 = stars%sk2(irec2)
                            phi = stars%phi2(irec2)
                            b = cell%z1*( gz*iJJ(m+1,gzi)*fJJ(m,irec2) + &
                                 &              g2*II(1)*fJJ(m+1,irec2))/(g*g)
                            a = a +    &
                                 &              (ic**m1)*mult*EXP(-ic*m1*phi)&
                                 &              *psq(irec3)*b*stars%nstr2(irec2)
                         END IF
                      END DO

                      !----- total vacuum potential ---------------

                      DO imz = 1,vacuum%nmzxy
                         pint(imz) = fpi_const*a*KK(imz) 
                         vxy(imz,irec1(l)-1,1,1) =  pint(imz) + pvac(imz)  
                      END DO

                      val_help = vxy(1,irec1(l)-1,1,1) - val(irec1(l))

                      CALL visp5_z(&
                           &        vacuum%nmzxyd,vacuum%nmzxyd,vacuum%delz,m,ivfft1,ivfft2,IIIR,&
                           &        rxy,ani1,ani2,cell%z1,cell%amat,pvac,pint,tpi_const,l,&
                           &        fpi_const,val_help,III,m1,gz,rmap,rr,&
                           &        vis_help)

                      DO ix = 0,ivfft1 - 1
                         DO iy = 0,ivfft2 - 1
                            vis(ix,iy,irec1(l)) = vis_help(ix,iy)
                         END DO
                      END DO

                      !---  end of vacuum contibution to vis ---------------------->

                   END IF                    ! gz.atoms%neq.0

                   !---- > finishing the cycle by +-m,gz

                END IF

             ENDDO
          END IF                    ! m^2+gz^2.atoms%neq.0

       ENDDO m_loop                  ! m

    ENDDO                  ! gz 

    !*************************************************************
    !-------> finding the Fourier coefficients of the potential
    !-------  in the interstitial
    !-------  the scheme (making the potential continuous) is as follows:
    !-------  1. transforming the \tilde{V} from pw-representation
    !-------     to the real grid, putting to zero outside the
    !-------     cylindrical interstitial
    !-------  2. Adding the \tilde\tilde{V} on the real grid 
    !-------  3. Transforming back to the pw-representation 
    !-------  Everything is done in \tilda {\Omega}

    gzmin = -stars%k3d

    IF (sym%zrfs.OR.sym%invs) gzmin = 0

    DO k3 = gzmin,stars%k3d          ! collect the Fourier components

       fxy0 = 0.

       rhti = 0.

       DO irec2 = 1,stars%ng2 - 1
          fxy(irec2) = CMPLX(0.,0.)
       END DO

       DO irec2 = 1,stars%ng2
          irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),k3)
          IF (irec3.NE.0) THEN
             IF (irec2.EQ.1) THEN
                fxy0 = REAL(vpw(irec3,1))
                rhti = AIMAG(vpw(irec3,1))
             ELSE
                fxy(irec2-1) = vpw(irec3,1)
             END IF
          END IF
       END DO

       af2=0.0
       bf2=0.0

       CALL fft2d(&
            &        stars,&
            &        af2(0),bf2,&
            &        fxy0,rhti,fxy,&
            &        1,+1)


       !--> although some harmonics of the vacuum potential are equal to zero
       !--> the same harmonics of the interstitial potential \Tilde{V} are not,
       !--> because this potential is determined from pw coefficients of the
       !--> charge density, which has the periodicity of the rectangular 
       !--> unit cell, hence, these harmonics should be extracted from the 
       !--> final interstitial potential

       i = 0

       DO iy = 0,3*stars%k2d-1
          DO ix = 0,3*stars%k1d-1
             x = ix*ani1
             IF (x.GT.0.5) x = x - 1.
             y = iy*ani2
             IF (y.GT.0.5) y = y - 1.
             r = SQRT((x*cell%amat(1,1) + y*cell%amat(1,2))**2 +&
                  &                  (x*cell%amat(2,1) + y*cell%amat(2,2))**2)               
             phi = angle(x*cell%amat(1,1) + y*cell%amat(1,2),&
                  &                     x*cell%amat(2,1) + y*cell%amat(2,2))
             vis_tot(ix,iy) = CMPLX(0.,0.)
             j = rmap(ix,iy)
             DO m = -MM,MM
                irc1 = oneD%ig1(k3,m)
                IF (irc1.NE.0) THEN
                   vis_tot(ix,iy) = vis_tot(ix,iy) +&
                        &                    vis(ix,iy,irc1)
                ELSE
                   IF (k3.EQ.0) THEN
                      IF (r.LE.cell%z1) THEN
                         vis_tot(ix,iy) = vis_tot(ix,iy) - &
                              &                          val_m(k3,m)*&
                              &                          EXP(CMPLX(0.,m*phi))*&
                              &                          (r**(iabs(m)))/(cell%z1**(iabs(m)))
                      END IF
                   ELSE
                      IF (r.LE.cell%z1) THEN
                         vis_tot(ix,iy) = vis_tot(ix,iy) -&
                              &                    val_m(k3,m)*IIII(j,iabs(k3),iabs(m))*&
                              &                    EXP(CMPLX(0.,m*phi))/IR(iabs(k3),iabs(m))
                      END IF
                   END IF
                END IF
             END DO

             IF (r.LE.cell%z1) THEN
                af2(i) = af2(i) + REAL(vis_tot(ix,iy))
                bf2(i) = bf2(i) + AIMAG(vis_tot(ix,iy))
                !$$$                 af2(i) = real(vis_tot(ix,iy))
                !$$$                 bf2(i) = aimag(vis_tot(ix,iy))
             ELSE
                af2(i) = REAL(vis_tot(ix,iy))
                bf2(i) = AIMAG(vis_tot(ix,iy))
             END IF
             i = i+1
          END DO
       END DO

       rhti = 0.

       CALL fft2d(&
            &        stars,&
            &        af2(0),bf2,&
            &        gxy0,rhti,gxy,&
            &        1,-1) 

       DO irec2 = 1,stars%ng2
          irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),k3)
          IF (irec3.NE.0) THEN
             IF (irec2.EQ.1) THEN
                vpw_help(irec3) = CMPLX(gxy0,rhti)
             ELSE
                vpw_help(irec3) = gxy(irec2-1)
             END IF
          END IF
       END DO

    END DO                    ! gz -> Vpw(.,.,gz)

    DO irec3 = 1,stars%ng3
       vpw(irec3,1) = vpw_help(irec3)
       !$$$         vpw(irec3,1) = vpw(irec3,1) + vpw_help(irec3)
    END DO


    DO irc1 = 2,nq2_1
       DO imz = 1,vacuum%nmzxy
          IF (ABS(vxy(imz,irc1-1,1,1)).LE.tol_21)&
               &          vxy(imz,irc1-1,1,1) = CMPLX(0.,0.)
       END DO
    END DO



    DEALLOCATE ( KK,II,III,IIII,fact,val,val_m,rxy,&
         &     vis,vis_tot,pvac,pint,vis_help,rmap,rr,&
         &     af2,bf2,vpw_help,fJJ,iJJ )

    RETURN
  END SUBROUTINE od_vvacis
END MODULE m_od_vvacis












