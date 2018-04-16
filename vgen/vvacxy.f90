MODULE m_vvacxy
  USE m_juDFT
  !     **********************************************************
  !     g.ne.0 coefficient of vacuum coulomb potential           *
  !     due to warped vacuum charged density                     *
  !                                 c.l.fu, r.podloucky          *
  !     **********************************************************
  !     modified for thick films to avoid underflows gb`06
  !---------------------------------------------------------------
CONTAINS
  SUBROUTINE vvacxy(&
       &                  stars,vacuum,cell,sym,input,&
       &                  rhtxy,&
       &                  vxy,&
       &                  alphm)

    USE m_intgr, ONLY : intgz1
    USE m_constants
    USE m_types
    USE m_qsf
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_stars),INTENT(IN)       :: stars
    TYPE(t_cell),INTENT(IN)        :: cell
    !     ..
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (IN)    :: rhtxy(vacuum%nmzxyd,stars%ng2-1,2)
    COMPLEX, INTENT (INOUT) :: vxy(vacuum%nmzxyd,stars%ng2-1,2)
    COMPLEX, INTENT (OUT)   :: alphm(stars%ng2,2)
    !     ..
    !     .. Local Scalars ..
    COMPLEX alph0 ,alph2,alph1,alphaz,betaz,test
    REAL g,vcons,z,e_m,e_p
    INTEGER imz,ip,irec2,ivac,ncsh
    LOGICAL tail
    !     ..
    !     .. Local Arrays ..
    REAL fra(vacuum%nmzxyd),frb(vacuum%nmzxyd),fia(vacuum%nmzxyd),fib(vacuum%nmzxyd)
    REAL alpha(vacuum%nmzxyd,2,2),beta(vacuum%nmzxyd,2,2)
    REAL, ALLOCATABLE :: sig_top(:), sig_bot(:)
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC aimag,cmplx,conjg,exp,REAL
    !     ..
    IF (ALLOCATED (input%efield%rhoEF) .OR. input%efield%dirichlet) THEN
       ncsh = input%efield%zsigma/vacuum%delz + 1.01
       ! If nmzxy < ncsh, the inhomogenous field cannot be represented.
       ! If nmzxy > ncsh, the boundary condition is wrong - and the
       ! potential is very wavy - try it yourself, if you don't believe.
       IF (vacuum%nmzxyd < ncsh) THEN
          WRITE (6,*) 'ERROR vvacxy.f: vacuum%nmzxyd =', vacuum%nmzxyd,&
               &                '< ncsh(zsigma) = ',ncsh
          CALL juDFT_error("ERROR: vacuum%nmzxyd < ncsh",calledby ="vvacxy")
       ELSE IF (vacuum%nmzxyd > ncsh) THEN
          WRITE (6,*) 'WARNING vvacxy.f: vacuum%nmzxyd =', vacuum%nmzxyd,&
               &                '> ncsh(zsigma) = ',ncsh
          WRITE (0,*) 'WARNING vvacxy.f: vacuum%nmzxyd =', vacuum%nmzxyd,&
               &                '> ncsh(zsigma) = ',ncsh
          CALL juDFT_warn("nmzxyd > ncsh",calledby ="vvacxy")
       END IF
    END IF

    !     ..
    !     2-dim star loop g.ne.0
    ip = vacuum%nmzxy + 1
    irec2_loop: DO irec2 = 2,stars%ng2
       g = stars%sk2(irec2)
       vcons = tpi_const/g

       ! ********** DIRICHLET ************************************
       IF (input%efield%dirichlet) THEN
          IF (ALLOCATED (input%efield%rhoEF)) THEN
             vxy(ncsh:vacuum%nmzxy,irec2-1,1) = input%efield%rhoEF(irec2-1, 1)
             IF (vacuum%nvac == 2) THEN
                vxy(ncsh:vacuum%nmzxy,irec2-1,2) = input%efield%rhoEF(irec2-1, 2)
             END IF
          ELSE
             vxy(ncsh:vacuum%nmzxy,irec2-1,1:vacuum%nvac) = 0.0
          END IF

          vcons = 2.0*vcons/SINH(g*2.0*(input%efield%zsigma+cell%z1))

          ivac_loop1: DO ivac = 1,vacuum%nvac
             z = cell%z1
             imz_loop1: DO imz = 1,ncsh-1
                ! As "z" > 0 in this subroutine, the integrand is the same
                ! for both ivac -- but the integral bounds are reversed
                e_m = SINH(g*(input%efield%zsigma+cell%z1-z))
                e_p = SINH(g*(z+input%efield%zsigma+cell%z1))
                fra(ncsh-imz) = REAL(rhtxy(imz,irec2-1,ivac))* e_m
                fia(ncsh-imz) = AIMAG(rhtxy(imz,irec2-1,ivac))*e_m
                frb(imz) = REAL(rhtxy(imz,irec2-1,ivac))* e_p
                fib(imz) = AIMAG(rhtxy(imz,irec2-1,ivac))*e_p
                z = z + vacuum%delz
             END DO imz_loop1
             CALL intgz1(fra,vacuum%delz,ncsh-1,alpha(1,ivac,1),tail=.FALSE.)
             CALL intgz1(fia,vacuum%delz,ncsh-1,alpha(1,ivac,2),tail=.FALSE.)
             CALL qsf(vacuum%delz,frb,beta(1,ivac,1),ncsh-1,1)
             CALL qsf(vacuum%delz,fib,beta(1,ivac,2),ncsh-1,1)
          END DO ivac_loop1

          IF (ivac == 2) THEN
             ! Honour reversed integral bounds
             alpha(:,ivac,1) = -alpha(:,ivac,1)
             alpha(:,ivac,2) = -alpha(:,ivac,2)
             beta(:,ivac,1) = -beta(:,ivac,1)
             beta(:,ivac,2) = -beta(:,ivac,2)
          END IF

          alph1 = CMPLX(alpha(ncsh-1,1,1),alpha(ncsh-1,1,2))
          IF (vacuum%nvac == 1) THEN
             IF (sym%invs) THEN
                alph2 = CONJG(alph1)
             ELSE
                alph2 = alph1
             END IF
          ELSE
             alph2 = CMPLX(alpha(ncsh-1,2,1),alpha(ncsh-1,2,2))
          END IF
          ivac_loop2: DO ivac = 1,vacuum%nvac
             z = cell%z1
             IF (ivac == 1) alph0 = alph2
             IF (ivac == 2) alph0 = alph1
             imz_loop2: DO imz = 1,ncsh-1
                betaz = CMPLX(beta(imz,ivac,1),beta(imz,ivac,2))
                alphaz = CMPLX(alpha(ncsh-imz,ivac,1),&
                     &                        alpha(ncsh-imz,ivac,2))
                e_m = SINH(g*(input%efield%zsigma+cell%z1-z))
                e_p = SINH(g*(z+input%efield%zsigma+cell%z1))
                test = e_m*(alph0+betaz) + e_p*alphaz
                IF ( 2.0 * test == test ) test = CMPLX(0.0,0.0)
                vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac)&
                     &                  +  vcons * test
                IF (ALLOCATED (input%efield%C1)) THEN
                   e_m = exp_save( -g*z  )
                   e_p = exp_save( g*z  )
                   IF (ivac == 1) THEN ! z > 0
                      vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac)&
                           &                                       + input%efield%C1(irec2-1)*e_p&
                           &                                       + input%efield%C2(irec2-1)*e_m
                   ELSE ! z < 0
                      vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac)&
                           &                                       + input%efield%C1(irec2-1)*e_m&
                           &                                       + input%efield%C2(irec2-1)*e_p
                   END IF
                END IF
                z = z + vacuum%delz
             END DO imz_loop2
          END DO ivac_loop2
          alphm(irec2-1,1) = alph1
          alphm(irec2-1,2) = alph2

          ! ********** NEUMANN ************************************
       ELSE
          ivac_loop3: DO ivac = 1,vacuum%nvac
             z = cell%z1
             imz_loop3: DO imz = 1,vacuum%nmzxy
                e_m = exp_save( -g*z )
                e_p = exp_save( g*z )
                fra(ip-imz) = REAL(rhtxy(imz,irec2-1,ivac))* e_m
                fia(ip-imz) = AIMAG(rhtxy(imz,irec2-1,ivac))*e_m
                frb(imz) = REAL(rhtxy(imz,irec2-1,ivac))* e_p
                fib(imz) = AIMAG(rhtxy(imz,irec2-1,ivac))*e_p
                z = z + vacuum%delz
             END DO imz_loop3

             ! Add external field, if segmented
             IF (ALLOCATED (input%efield%rhoEF)) THEN
                z = cell%z1 + input%efield%zsigma
                e_m = exp_save( -g*z )
                e_p = exp_save( g*z )

                ! The equation has a minus sign as "rhtxy" contains the electron density
                ! (a positive number representing a negative charge) while rhoEF
                ! specifies the charges in terms of the (positive) elementary charge "e".
                fra(ip-ncsh) = fra(ip-ncsh)&
                     &                        - REAL (input%efield%rhoEF(irec2-1, ivac))*e_m
                fia(ip-ncsh) = fia(ip-ncsh)&
                     &                        - AIMAG (input%efield%rhoEF(irec2-1, ivac))*e_m
                frb(ncsh)    = frb(ncsh)&
                     &                       - REAL (input%efield%rhoEF(irec2-1, ivac))*e_p
                fib(ncsh)    = fib(ncsh)&
                     &                       - AIMAG (input%efield%rhoEF(irec2-1, ivac))*e_p
             END IF
             CALL intgz1(fra,vacuum%delz,vacuum%nmzxy,alpha(1,ivac,1),tail=.TRUE.)
             CALL intgz1(fia,vacuum%delz,vacuum%nmzxy,alpha(1,ivac,2),tail=.TRUE.)
             CALL qsf(vacuum%delz,frb,beta(1,ivac,1),vacuum%nmzxy,1)
             CALL qsf(vacuum%delz,fib,beta(1,ivac,2),vacuum%nmzxy,1)
          END DO ivac_loop3

          alph1 = CMPLX(alpha(vacuum%nmzxy,1,1),alpha(vacuum%nmzxy,1,2))
          IF (vacuum%nvac.EQ.1) THEN
             IF (sym%invs) THEN
                alph2 = CONJG(alph1)
             ELSE
                alph2 = alph1
             END IF
          ELSE
             alph2 = CMPLX(alpha(vacuum%nmzxy,2,1),alpha(vacuum%nmzxy,2,2))
          END IF
          ivac_loop4: DO ivac = 1,vacuum%nvac
             z = cell%z1
             IF (ivac.EQ.1) alph0 = alph2
             IF (ivac.EQ.2) alph0 = alph1
             imz_loop4: DO imz = 1,vacuum%nmzxy
                betaz = CMPLX(beta(imz,ivac,1),beta(imz,ivac,2))
                alphaz = CMPLX(alpha(ip-imz,ivac,1),alpha(ip-imz,ivac,2))
                e_m = exp_save( -g*z  )
                e_p = exp_save( g*z  )
                test = e_m*(alph0+betaz) + e_p*alphaz
                IF ( 2.0 * test == test ) test = CMPLX(0.0,0.0)
                vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac) +&
                     &                    vcons * test
                z = z + vacuum%delz
             END DO imz_loop4
          END DO ivac_loop4
          alphm(irec2-1,1) = alph1
          alphm(irec2-1,2) = alph2
       END IF ! Neumann (vs. Dirichlet)
    END DO irec2_loop
  END SUBROUTINE vvacxy

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

END MODULE m_vvacxy
