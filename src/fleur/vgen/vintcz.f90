!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vintcz
  !     *************************************************************
  !     z-dependent part of Coulomb potential in the interstitial   *
  !     [-d/2,d/2] region            c.l.fu, r.podloucky            *
  !     *************************************************************
  !     modified for thick films to avoid underflows gb`06
  !---------------------------------------------------------------
   USE m_constants
   USE m_cfft
   USE m_lagrange3
   USE m_types_cell
   USE m_types_field
   USE m_types_stars
   USE m_types_vacuum
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: vintcz, vintcz_sumgrid, exp_safe
CONTAINS
   SUBROUTINE vintcz_sumgrid(stars,cell,nrec2,psq,n,vsum)
      !! Sum over q_z in the interstitial part of vintcz for a star with G_par /= 0
      !! (Neumann boundary condition) on the film z-grid z_i = a3*i/n, i=0..n-1
      !! (a point with z_i > a3/2 stands for z_i - a3).
      !!
      !! exp(g(z-dh)), exp(-g(z+dh)) and the phases exp(+-i q dh) do not depend on q_z, so
      !!   sum_q c_q [ e^{iqz} + vcons2 ( (g+iq) e_p(z) e^{iq dh} + (g-iq) e_m(z) e^{-iq dh} ) ]
      !!     = F(z) + vcons2 ( A e_p(z) + B e_m(z) ),
      !! with F(z) = sum_q c_q e^{iqz}, A = sum_q c_q (g+iq) e^{iq dh}, B = sum_q c_q (g-iq) e^{-iq dh}.
      !! A and B cost O(mx3); F on the whole grid is one inverse FFT. The result is the value
      !! vintcz would sum up for |z| < z1; for |z| >= z1 vintcz does not use it.
      IMPLICIT NONE

      TYPE(t_stars), INTENT(IN)  :: stars
      TYPE(t_cell),  INTENT(IN)  :: cell
      INTEGER,       INTENT(IN)  :: nrec2, n
      COMPLEX,       INTENT(IN)  :: psq(stars%ng3)
      COMPLEX,       INTENT(OUT) :: vsum(0:n-1)

      REAL    :: a(n), b(n)
      REAL    :: dh, g, g3, q, vcons2, cos_q, sin_q, z, ani, e_m, e_p
      COMPLEX :: vcons1, ca, cb
      INTEGER :: iq, ig3n, k1, k2, i

      dh = cell%z1
      k1 = stars%kv2(1,nrec2)
      k2 = stars%kv2(2,nrec2)
      g = stars%sk2(nrec2)
      vcons2 = - 1.0 / (2.*g)
      ani = 1.0 / REAL(n)

      a = 0.0
      b = 0.0
      ca = CMPLX(0.0,0.0)
      cb = CMPLX(0.0,0.0)
      DO iq = -stars%mx3, stars%mx3
         ig3n = stars%ig(k1,k2,iq)
         IF (ig3n == 0) CYCLE
         q = iq*cell%bmat(3,3)
         g3 = stars%sk3(ig3n)
         vcons1 = fpi_const*psq(ig3n)*stars%rgphs(k1,k2,iq) / (g3*g3)
         ! the frequencies |iq| <= mx3 < n/2 do not overlap on the grid
         i = MODULO(iq,n) + 1
         a(i) = REAL(vcons1)
         b(i) = AIMAG(vcons1)
         cos_q = COS(q*dh)
         sin_q = SIN(q*dh)
         ca = ca + vcons1*(g + ImagUnit*q)*CMPLX(cos_q, sin_q)
         cb = cb + vcons1*(g - ImagUnit*q)*CMPLX(cos_q, -sin_q)
      END DO

      CALL cfft(a, b, n, n, n, 1) ! sum_k c_k exp(+2 pi i k j/n)

      DO i = 0, n-1
         z = cell%amat(3,3) * i * ani
         IF (z > cell%amat(3,3)/2.) z = z - cell%amat(3,3)
         e_m = exp_safe( -g*(z+dh) ) !exp_safe handles overflow
         e_p = exp_safe( g*(z-dh) )
         vsum(i) = CMPLX(a(i+1), b(i+1)) + vcons2*(ca*e_p + cb*e_m)
      END DO

   END SUBROUTINE vintcz_sumgrid

   COMPLEX FUNCTION vintcz(stars,vacuum,cell,field,z,nrec2,psq,vnew,rhobar,sig1dh,vz1dh,alphm,vslope,l_dfptvgen,vsum)

      IMPLICIT NONE

      TYPE(t_stars),  INTENT(IN) :: stars
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      TYPE(t_cell),   INTENT(IN) :: cell
      TYPE(t_field),  INTENT(IN) :: field

      INTEGER,        INTENT(IN) :: nrec2
      COMPLEX,        INTENT(IN) :: rhobar,vslope
      COMPLEX,        INTENT(IN) :: sig1dh,vz1dh
      REAL,           INTENT(IN) :: z

      COMPLEX,        INTENT(IN) :: psq(stars%ng3),vnew(:,:,:)
      COMPLEX,        INTENT(IN) :: alphm(stars%ng2,2)
      LOGICAL,        INTENT(IN) :: l_dfptvgen
      COMPLEX, OPTIONAL, INTENT(IN) :: vsum ! sum over q_z (G_par /= 0, Neumann) at this z, from vintcz_sumgrid

      COMPLEX                    :: sumrr,vcons1,test,c_ph,phas
      REAL                       :: bj0,dh,fit,g,g3,q,qdh,vcons2,zf
      REAL                       :: e_m,e_p,cos_q,sin_q
      INTEGER                    :: ig3n,im,iq,ivac,k1,k2,nrec2r
      
      dh = cell%z1
      sumrr = (0.,0.)
      vintcz = (0.,0.)

      
      !--->    if z is in the vacuum, use vacuum representations (m.w.)
      IF (ABS(z).GE.cell%z1) THEN
         ivac = 1
         IF (z.LT.0.0) THEN
            ivac = 2
            IF (vacuum%nvac==1) ivac = 1
         END IF
         zf = (ABS(z)-cell%z1)/vacuum%delz + 1.0
         im = zf
         q = zf - im
         ! For q/=0 in DFPT, the is no G+q=0, so all stars are treated in the G/=0 way.
         IF (nrec2.EQ.1.AND.((.NOT.l_dfptvgen).OR.norm2(stars%center)<1e-8)) THEN
            fit = lagrange3(q,REAL(vnew(im,1,ivac)),REAL(vnew(im+1,1,ivac)),REAL(vnew(im+2,1,ivac)))
            vintcz = CMPLX(fit,0.0)
         ELSE IF (im+2.LE.vacuum%nmzxy) THEN
            if (z<0) THEN
               call stars%map_2nd_vac(vacuum,nrec2,nrec2r,phas)
            else
               nrec2r=nrec2
               phas=cmplx(1.0,0.0)
            end if
            vintcz = phas*lagrange3(q,vnew(im,nrec2r,ivac),vnew(im+1,nrec2r,ivac),vnew(im+2,nrec2r,ivac))
         END IF
         RETURN
      END IF

      ! For q/=0 in DFPT, there is no G+q=0, so all stars are treated in the G/=0 way.
      IF (nrec2==1.AND.((.NOT.l_dfptvgen).OR.norm2(stars%center)<1e-8)) THEN    !     ---->    g=0 coefficient

         DO  iq = -stars%mx3,stars%mx3
            IF (iq.EQ.0) CYCLE
            ig3n = stars%ig(0,0,iq)
            !     ----> use only stars within the g_max sphere (oct.97 shz)
            IF (ig3n.NE.0) THEN
               q = iq*cell%bmat(3,3)
               qdh = q*dh
               bj0 = SIN(qdh)/qdh
               sumrr = (cmplx(cos(q*z),sin(q*z)) - cmplx(cos(qdh),sin(qdh)))/(q*q) + ImagUnit*cos(qdh)*(dh-z)/q + bj0*(z*z-dh*dh)/2.
               vintcz = vintcz + fpi_const*psq(ig3n)*sumrr
            END IF
         END DO
         !           -----> v2(z)
         vintcz = vintcz + vz1dh - fpi_const* (dh-z)*&
            &              (sig1dh-rhobar/2.* (dh-z))
         

         IF (field%efield%dirichlet .AND. vslope /= 0.0) THEN
            vintcz = vintcz + vslope * (dh-z)
         END IF
         !     ---->    (g.ne.0)  coefficients
      ELSE
         k1 = stars%kv2(1,nrec2)
         k2 = stars%kv2(2,nrec2)
         g = stars%sk2(nrec2)

         IF (PRESENT(vsum) .AND. .NOT.field%efield%dirichlet) THEN
            vintcz = vsum
         ELSE
         DO  iq = -stars%mx3,stars%mx3
            ig3n = stars%ig(k1,k2,iq)
            !     ----> use only stars within the g_max sphere (oct.97 shz)
            IF (ig3n.NE.0) THEN
               c_ph = stars%rgphs(k1,k2,iq)
               !           -----> v3(z)
               q = iq*cell%bmat(3,3)
               g = stars%sk2(nrec2)
               g3 = stars%sk3(ig3n)
               vcons1 = fpi_const*psq(ig3n)*c_ph / (g3*g3)
               IF (field%efield%dirichlet) THEN
                  e_m = 0.0
                  e_p = 0.0
                  vcons1  = vcons1/(g*SINH(g*2*(dh+field%efield%zsigma)))
                  e_m = e_m - EXP(-ImagUnit*q*dh)*(g*COSH(g*(-field%efield%zsigma))+ImagUnit*q*SINH(g*(-field%efield%zsigma)))
                  e_p = e_p + EXP(ImagUnit*q*dh)*(-g*COSH(g*(-field%efield%zsigma))+ImagUnit*q*SINH(g*(-field%efield%zsigma)))
                  sumrr = EXP(ImagUnit*q*z)
                  e_m = e_m + sumrr*(g*COSH(g*(z+dh+field%efield%zsigma))-ImagUnit*q*SINH(g*(z+dh+field%efield%zsigma)))
                  e_p = e_p + sumrr*(g*COSH(g*(z-dh-field%efield%zsigma))-ImagUnit*q*SINH(g*(z-dh-field%efield%zsigma)))
                  vintcz = vintcz+ vcons1*(e_m*SINH(g*(field%efield%zsigma+dh-z))&
                                          +e_p*SINH(g*(field%efield%zsigma+dh+z)))
               ELSE
                  sumrr = (0.0,0.0)
                  vcons2 = - 1.0 / (2.*g)
                  e_m = exp_safe( -g*(z+dh) ) !exp_safe handles overflow
                  e_p = exp_safe( g*(z-dh) )
                  !               --> sum over gz-stars
                  cos_q = COS(q*dh)
                  sin_q = SIN(q*dh)
                  sumrr = sumrr + CMPLX(COS(q*z),SIN(q*z)) + vcons2 *&
                        &                      ( (g + ImagUnit*q) * e_p * (cos_q + ImagUnit*sin_q) +&
                        &                        (g - ImagUnit*q) * e_m * (cos_q - ImagUnit*sin_q) )
                  vintcz = vintcz + vcons1*sumrr

               END IF ! Neumann (vs. Dirichlet)
            END IF ! ig3d /= 0
         END DO
         END IF ! precomputed sum over q_z
         !  ----> v4(z)
         IF (field%efield%dirichlet) THEN
            e_m = SINH(g*(field%efield%zsigma+dh - z))
            e_p = SINH(g*(field%efield%zsigma+dh + z))
            test = e_m*alphm(nrec2,2) + e_p*alphm(nrec2,1)
            test = fpi_const/(g*SINH(g*2*(field%efield%zsigma+dh))) * test
            IF ( 2.0 * test == test ) test = CMPLX(0.0,0.0)
            vintcz = vintcz + test
            IF (ALLOCATED (field%efield%C1)) THEN
               g = stars%sk2(nrec2)
               e_m = exp_safe (-g*z)
               e_p = exp_safe ( g*z)
               vintcz = vintcz + field%efield%C1(nrec2-1)*e_m &
                  &            + field%efield%C2(nrec2-1)*e_p
            END IF
         ELSE ! Neumann
            e_m = exp_safe( -g*z  )
            e_p = exp_safe( g*z  )
            test = e_m*alphm(nrec2,2) + e_p*alphm(nrec2,1)
            IF ( 2.0 * test == test ) test = CMPLX(0.0,0.0)
            vintcz = vintcz + tpi_const/g* test

         END IF
      END IF
   END FUNCTION vintcz

   PURE REAL FUNCTION exp_safe(x)
      ! replace exp by a function that does not under/overflow dw09
      IMPLICIT NONE

      REAL, INTENT(IN) :: x
      REAL, PARAMETER  :: maxexp = LOG(2.0)*MAXEXPONENT(2.0)
      REAL, PARAMETER  :: minexp = LOG(2.0)*MINEXPONENT(2.0)

      IF ( ABS(x)>minexp .AND. ABS(x)<maxexp ) THEN
         exp_safe = EXP(x)
      ELSE
         IF ( x > 0 ) THEN
            IF ( x > minexp ) THEN
               exp_safe = EXP(maxexp)
            ELSE
               exp_safe = EXP(minexp)
            END IF
         ELSE
            IF ( -x > minexp ) THEN
               exp_safe = EXP(-maxexp)
            ELSE
               exp_safe = EXP(-minexp)
            END IF
         END IF
      END IF
   END FUNCTION exp_safe
END MODULE m_vintcz
