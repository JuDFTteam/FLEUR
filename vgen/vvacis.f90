MODULE m_vvacis
  !     **********************************************************
  !     g.ne.0 coefficients of vacuum coulomb potential          *
  !     due to the interstitial charge density inside slab       *
  !                                   c.l.fu, r.podloucky        *
  !     **********************************************************
  !     modified for thick films to avoid underflows gb`06
  !---------------------------------------------------------------
CONTAINS
  SUBROUTINE vvacis(&
       &                  stars,vacuum,&
       &                  sym,cell,&
       &                  psq, input,&
       &                  vxy)

    USE m_constants
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)  :: vacuum
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_cell),INTENT(IN)    :: cell
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (IN) :: psq(stars%n3d)
    COMPLEX, INTENT (OUT):: vxy(vacuum%nmzxyd,stars%n2d-1,2,input%jspins)
    !     ..
    !     .. Local Scalars ..
    COMPLEX arg,ci
    REAL dh,g,qz,sign,signz,vcons,z,e_m
    REAL arg_r,arg_i
    INTEGER i2d,ig3n,imz,imzxy,ivac,k1,k2,kz,m0,nrec2,nrz,nz
    !     ..
    !     .. Local Arrays ..
    COMPLEX sumr(2)
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC exp
    !     ..
    ci = CMPLX(0.0,1.0)

    vxy(:,:,:,1) = CMPLX(0.,0.)
    dh = cell%z1
    m0 = -stars%mx3
    IF (sym%zrfs) m0 = 0
    DO  nrec2 = 2,stars%ng2
       k1 = stars%kv2(1,nrec2)
       k2 = stars%kv2(2,nrec2)
       g = stars%sk2(nrec2)
       IF (input%efield%dirichlet) THEN
          vcons = 2.0*tpi_const/(g*SINH(g*2.0*(input%efield%zsigma+dh)))
          arg_r = g*(dh+input%efield%zsigma+dh)
       ELSE ! Neumann
          vcons = tpi_const/g
          arg_r = exp_save( - 2*dh*g )
       END IF
       DO ivac = 1,vacuum%nvac
          sumr(ivac) = (0.0,0.0)
          sign = 3. - 2.*ivac
          DO kz = m0,stars%mx3
             ig3n = stars%ig(k1,k2,kz)
             !     ----> use only stars within the g_max sphere (oct.97 shz)
             IF (ig3n.NE.0) THEN
                nz = 1
                IF (sym%zrfs) nz = stars%nstr(ig3n)/stars%nstr2(nrec2)
                qz = kz*cell%bmat(3,3)
                !     ---> sum over gz-stars
                DO  nrz = 1,nz
                   signz = 3. - 2.*nrz
                   IF (input%efield%dirichlet) THEN
                      ! prefactor
                      arg = EXP(-ci*signz*qz*dh)&
                           &                      /(2*(g**2 + qz**2)) * psq(ig3n)
                      IF (ivac == 1) THEN
                         sumr(ivac) = sumr(ivac) + EXP(-arg_r)*arg*(&
                              &                     (- EXP(2*g*(input%efield%zsigma+dh))&
                              &                      + EXP(2*(ci*signz*qz*dh+arg_r)))&
                              &                     *(g-ci*signz*qz)&
                              &                    +(- EXP(2*g*dh)&
                              &                      + EXP(2*ci*signz*qz*dh))&
                              &                     *(g+ci*signz*qz) )
                      ELSE
                         sumr(ivac) = sumr(ivac) + arg*(&
                              &                     EXP(arg_r)*(g+ci*signz*qz)&
                              &                     +(g-ci*signz*qz)*EXP(-arg_r)&
                              &                     +2*EXP(2*(ci*signz*qz*dh))&
                              &                      *(-g*COSH(g*(-input%efield%zsigma))&
                              &                        +ci*signz*qz*SINH(g*(-input%efield%zsigma))) )
                      END IF
                   ELSE
                      arg = g + sign*ci*signz*qz
                      arg_i = sign*signz*qz*dh
                      sumr(ivac) = sumr(ivac) + psq(ig3n)*(&
                           &                     COS(arg_i)*( 1 - arg_r ) +&
                           &                  ci*SIN(arg_i)*( 1 + arg_r ) ) / arg
                   END IF
                ENDDO  ! nrz 
             ENDIF
          ENDDO       ! kz 
          z = 0 ! moved cell%z1 into above equations gb`06
          DO imz = 1,vacuum%nmzxy
             IF (input%efield%dirichlet) THEN
                e_m = SINH(g*(input%efield%zsigma-z))
             ELSE ! NEUMANN
                e_m = exp_save( -g*z  )
             END IF
             vxy(imz,nrec2-1,ivac,1) = vxy(imz,nrec2-1,ivac,1) +&
                  &                                   vcons*sumr(ivac)*e_m
             z = z + vacuum%delz
          ENDDO  ! imz 
       ENDDO     ! ivac
    ENDDO        ! nrec2

  END SUBROUTINE vvacis
  !------------------------------------------------------------------
  PURE REAL FUNCTION exp_save(x)
    ! replace exp by a function that does not under/overflow dw09
    IMPLICIT NONE
    REAL   ,INTENT(IN)     :: x
    REAL, PARAMETER ::    maxexp = LOG(2.0)*MAXEXPONENT(2.0)
    REAL, PARAMETER ::    minexp = LOG(2.0)*MINEXPONENT(2.0)

    IF ( ABS(x)>minexp .AND. ABS(x)<maxexp ) THEN
       exp_SAVE = EXP(x)
    ELSE
       IF ( x > 0 ) THEN
          IF ( x > minexp ) THEN
             exp_save = EXP(maxexp)
          ELSE
             exp_save = EXP(minexp)
          ENDIF
       ELSE
          IF ( -x > minexp ) THEN
             exp_save = EXP(-maxexp)
          ELSE
             exp_save = EXP(-minexp)
          ENDIF
       ENDIF
    ENDIF
  END FUNCTION exp_save

END MODULE m_vvacis
