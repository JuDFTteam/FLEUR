MODULE m_vintcz
  !     *************************************************************
  !     z-dependent part of Coulomb potential in the interstitial   *
  !     [-d/2,d/2] region            c.l.fu, r.podloucky            *
  !     *************************************************************
  !     modified for thick films to avoid underflows gb`06
  !---------------------------------------------------------------
CONTAINS
   COMPLEX FUNCTION vintcz(stars,vacuum,cell,input,field,z,nrec2,psq,vnew,rhobar,sig1dh,vz1dh,alphm,vslope,sigma_disc,l_dfptvgen,l_corr,diff_vmz1dh,sigma_disc2)
      USE m_constants
      USE m_types

      IMPLICIT NONE

      TYPE(t_stars),  INTENT(IN) :: stars
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      TYPE(t_cell),   INTENT(IN) :: cell
      TYPE(t_input),  INTENT(IN) :: input
      TYPE(t_field),  INTENT(IN) :: field

      INTEGER,        INTENT(IN) :: nrec2
      COMPLEX,        INTENT(IN) :: rhobar,vslope
      COMPLEX,        INTENT(IN) :: sig1dh,vz1dh
      REAL,           INTENT(IN) :: z

      COMPLEX,        INTENT(IN) :: psq(stars%ng3),vnew(:,:,:)!,vxy(:,:,:) !(vacuum%nmzxyd,stars%ng2-1,2)
      COMPLEX,        INTENT(IN) :: alphm(stars%ng2,2)
      complex,        intent(in) :: sigma_disc(2)
      LOGICAL,        INTENT(IN) :: l_dfptvgen, l_corr
      COMPLEX,        INTENT(IN) :: diff_vmz1dh

      complex, optional, intent(in) :: sigma_disc2(2)
      !REAL,           INTENT(IN) :: vz(:,:) !(vacuum%nmzd,2,jspins)

      COMPLEX                    :: argr,sumrr,vcons1,test,c_ph,phas
      REAL                       :: bj0,dh,fit,g,g3,q,qdh,vcons2,zf
      REAL                       :: e_m,e_p,cos_q,sin_q
      INTEGER                    :: ig3n,im,iq,ivac,k1,k2,nrec2r
      COMPLEX                    :: newdp, newdm, newdp2, newdm2
      
      dh = cell%z1
      sumrr = (0.,0.)
      vintcz = (0.,0.)

      newdp = CMPLX(0.0,0.0)
      newdm = CMPLX(0.0,0.0)
      newdp2 = CMPLX(0.0,0.0)
      newdm2 = CMPLX(0.0,0.0)
      
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
            fit = 0.5* (q-1.)* (q-2.)*REAL(vnew(im,1,ivac)) -&
               &            q* (q-2.)*REAL(vnew(im+1,1,ivac)) +&
               &        0.5*q* (q-1.)*REAL(vnew(im+2,1,ivac))
            vintcz = CMPLX(fit,0.0)
         ELSE IF (im+2.LE.vacuum%nmzxy) THEN
            if (z<0) THEN 
               call stars%map_2nd_vac(vacuum,nrec2,nrec2r,phas) ! TODO: AN TB; will this work?
            else
               nrec2r=nrec2 
               phas=cmplx(1.0,0.0)
            end if      
            vintcz = phas*0.5* (q-1.)* (q-2.)*vnew(im,nrec2r,ivac) -&
                                    q* (q-2.)*vnew(im+1,nrec2r,ivac) +&
                                0.5*q* (q-1.)*vnew(im+2,nrec2r,ivac)
         END IF
         RETURN
      END IF

      ! For q/=0 in DFPT, there is no G+q=0, so all stars are treated in the G/=0 way.
      IF (nrec2==1.AND.((.NOT.l_dfptvgen).OR.norm2(stars%center)<1e-8)) THEN    !     ---->    g=0 coefficient

         ! New discontinuity correction
         newdp = - psq(1) * dh
         newdp2 = -psq(1) * dh**2 

         DO  iq = -stars%mx3,stars%mx3
            IF (iq.EQ.0) CYCLE
            ig3n = stars%ig(0,0,iq)
            !     ----> use only stars within the g_max sphere (oct.97 shz)
            IF (ig3n.NE.0) THEN
               q = iq*cell%bmat(3,3)
               qdh = q*dh
               bj0 = SIN(qdh)/qdh
               argr = ImagUnit*q*z
               sumrr = (EXP(argr)-EXP(ImagUnit*qdh))/ (q*q) + ImagUnit*COS(qdh)* (dh-z)/q + bj0* (z*z-dh*dh)/2.
               vintcz = vintcz + fpi_const*psq(ig3n)*sumrr
               
               ! New discontinuity correction
               newdp = newdp + ImagUnit * psq(ig3n) * cmplx(cos(qdh), sin(qdh)) * dh / qdh
               newdp2 = newdp2 + psq(ig3n) * cmplx(cos(qdh), sin(qdh)) / q**2

            END IF
         END DO
         !           -----> v2(z)
         vintcz = vintcz + vz1dh - fpi_const* (dh-z)*&
            &              (sig1dh-rhobar/2.* (dh-z))
         ! Discontinuity correction
         !vintcz = vintcz - fpi_const * sigma_disc(1)*(dh-z)
         !IF (PRESENT(sigma_disc2)) vintcz = vintcz + fpi_const * sigma_disc2(1)

         ! New discontinuity correction
         IF (l_dfptvgen) vintcz = vintcz - fpi_const * newdp * (dh-z) + fpi_const * newdp2

         ! Correct the interstitial potential by the mismatch we calculated as
         ! a boundary condition for DFPT q=0.
         !IF (l_dfptvgen) vintcz = vintcz - diff_vmz1dh*z/(2*dh) + diff_vmz1dh/2

         IF (field%efield%dirichlet .AND. vslope /= 0.0) THEN
            vintcz = vintcz + vslope * (dh-z)
         END IF
         !     ---->    (g.ne.0)  coefficients
      ELSE
         k1 = stars%kv2(1,nrec2)
         k2 = stars%kv2(2,nrec2)

         newdp = cmplx(0.0,0.0)
         newdm = cmplx(0.0,0.0)
         newdp2 = cmplx(0.0,0.0)
         newdm2 = cmplx(0.0,0.0)

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

                  ! New discontinuity correction
                  IF (iq==0) newdp = newdp - psq(ig3n) * dh
                  IF (iq/=0) newdp = newdp + ImagUnit * psq(ig3n) * cmplx(cos_q, sin_q) * dh / (q*dh)
                  IF (iq==0) newdm = newdm - psq(ig3n) * dh
                  IF (iq/=0) newdm = newdm - ImagUnit * psq(ig3n) * cmplx(cos_q,-sin_q) * dh / (q*dh)
                  IF (iq==0) newdp2 = newdp2 - psq(ig3n) * dh**2
                  IF (iq/=0) newdp2 = newdp2 + psq(ig3n) * cmplx(cos_q, sin_q) / q**2
                  IF (iq==0) newdm2 = newdm2 + psq(ig3n) * dh**2
                  IF (iq/=0) newdm2 = newdm2 - psq(ig3n) * cmplx(cos_q,-sin_q) / q**2

               END IF ! Neumann (vs. Dirichlet)
            END IF ! ig3d /= 0
         END DO
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

            ! New discontinuity correction
            !IF (l_dfptvgen.AND.l_corr.AND.nrec2==1) THEN
            !   e_m = exp_safe( - g * abs(z-cell%z1) )
            !   e_p = exp_safe( - g * abs(z+cell%z1) )
            !   test = e_m * newdp + e_p * newdm
            !   if ( 2.0 * test == test ) test = cmplx( 0.0, 0.0 )
            !   vintcz = vintcz + tpi_const/g * test
            !   !IF (abs(z-cell%z1)<1e-9) e_m = 0.0
            !   !IF (abs(z+cell%z1)<1e-9) e_p = 0.0
            !   test = e_m * newdp2 * sign(1.0,z-cell%z1)+ e_p * newdm2 * sign(1.0,z+cell%z1)
            !   if ( 2.0 * test == test ) test = cmplx( 0.0, 0.0 )
            !   vintcz = vintcz - 0*tpi_const * test
            !END IF

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
