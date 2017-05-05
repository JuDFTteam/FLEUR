MODULE m_intnv
  !     ************************************************
  !     calculates the integral of charge density 
  !     and potential in the unit cell
  !     ************************************************
CONTAINS
  SUBROUTINE int_nv(&
       &                  stars,vacuum,atoms,sphhar,&
       &                  cell,sym,input,oneD,&
       &                  qpw,vpw_w,&
       &                  rhtxy,vxy,&
       &                  rht,vz,&
       &                  rho,vr,&
       &                  RESULT)

    USE m_intgr, ONLY : intgr3,intgz0
    USE m_types
    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..
    REAL  RESULT  
    TYPE(t_stars),INTENT(IN)  :: stars
    TYPE(t_vacuum),INTENT(IN) :: vacuum
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_sphhar),INTENT(IN) :: sphhar
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_oneD),INTENT(IN)   :: oneD

    !     ..
    !     .. Array Arguments ..
    COMPLEX qpw(stars%ng3),rhtxy(vacuum%nmzxyd,oneD%odi%n2d-1,2),&
         &       vpw_w(stars%ng3),vxy(vacuum%nmzxyd,oneD%odi%n2d-1,2)
    REAL    rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype),rht(vacuum%nmzd,2),&
         &      vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype),vz(vacuum%nmzd,2)
   
    !     ..
    !     .. Local Scalars ..
    REAL dpdot,facv,tis,tmt,tvac,tvact
    INTEGER i,ip,ivac,j,k2,lh,n,npz,nat
    LOGICAL tail
    !     ..
    !     .. Local Arrays ..
    REAL dpj(atoms%jmtd),dpz(vacuum%nmzd)
    !     ..
    !     ..
    !     -----> CALCULATE DENSITY-POTENTIAL INTEGRALS
    !
    !  ******************* INTERSTITIAL REGION**********************
    !
    !  -> warping has been moved to vgen and visxc resp. ...gustav
    !
    tis = cell%omtil * REAL( DOT_PRODUCT(vpw_w,qpw))

    WRITE (6,FMT=8020) tis
    WRITE (16,FMT=8020) tis
8020 FORMAT (/,10x,'interstitial :',t40,f20.10)

    RESULT = RESULT + tis
    !
    !   ******************M.T. SPHERES*******************
    !
    tmt = 0.
    nat = 1
    DO n = 1,atoms%ntype
       DO lh = 0,sphhar%nlh(atoms%ntypsy(nat))
          DO j = 1,atoms%jri(n)
             dpj(j) = rho(j,lh,n)*vr(j,lh,n)
          ENDDO
          CALL intgr3(dpj,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),dpdot)
          tmt = tmt + dpdot*atoms%neq(n)
       ENDDO
       nat = nat + atoms%neq(n)
    ENDDO
    WRITE (6,FMT=8030) tmt
    WRITE (16,FMT=8030) tmt
8030 FORMAT (/,10x,'muffin tin spheres :',t40,f20.10)
    RESULT = RESULT + tmt
    !
    ! *********** VACUUM REGION**************
    !
    IF (input%film .AND. .NOT.oneD%odi%d1) THEN
       npz = vacuum%nmz + 1
       tail = .TRUE.
       IF (sym%zrfs .OR. sym%invs) THEN
          facv = 2.0
       ELSE
          facv = 1.0
       END IF
       tvac = 0.
       tvact = 0.
       !     set array dpz to zero
       dpz=0.0
       DO ivac = 1,vacuum%nvac
          DO ip = 1,vacuum%nmz
             dpz(npz-ip) = rht(ip,ivac)*vz(ip,ivac)
             !         --->  WARPING REGION
          ENDDO
          DO  k2 = 2,stars%ng2
             DO  ip = 1,vacuum%nmzxy
                dpz(npz-ip) = dpz(npz-ip) +&
                     &                       stars%nstr2(k2)*rhtxy(ip,k2-1,ivac)*&
                     &                          CONJG(vxy(ip,k2-1,ivac))
             ENDDO
          ENDDO
          CALL intgz0(dpz,vacuum%delz,vacuum%nmz,tvac,tail)
          tvact = tvact + cell%area*tvac*facv
       ENDDO
       WRITE (6,FMT=8040) tvact
       WRITE (16,FMT=8040) tvact
8040   FORMAT (/,10x,'vacuum :',t40,f20.10)
       RESULT = RESULT + tvact
    ELSEIF (oneD%odi%d1) THEN
       !-odim
       npz = vacuum%nmz +1
       tail = .TRUE.
       tvac = 0.
       tvact = 0.
       !     set array dpz to zero
       dpz=0.0
       DO  ip = 1,vacuum%nmz
          dpz(npz-ip) = (cell%z1+vacuum%delz*(ip-1))*&
               &                    rht(ip,vacuum%nvac)*vz(ip,vacuum%nvac)
          !          ---> WARPING REGION
       ENDDO
       DO  k2 = 2,oneD%odi%nq2
          DO  ip = 1,vacuum%nmzxy
             dpz(npz-ip) = dpz(npz-ip)+&
                  &             (cell%z1+vacuum%delz*(ip-1))*&
                  &             rhtxy(ip,k2-1,vacuum%nvac)*&
                  &             CONJG(vxy(ip,k2-1,vacuum%nvac))
          ENDDO
       ENDDO

       CALL intgz0(dpz,vacuum%delz,vacuum%nmz,tvac,tail)
       tvact = tvact + cell%area*tvac
       WRITE (6,FMT=8041) tvact
       WRITE (16,FMT=8041) tvact
8041   FORMAT (/,10x,'vacuum :',t40,f20.10)
       RESULT = RESULT + tvact
       !+odim
    END IF

  END SUBROUTINE int_nv
END MODULE m_intnv
