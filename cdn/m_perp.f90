!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_m_perp
CONTAINS
  SUBROUTINE m_perp(atoms,itype,iRepAtom,noco,nococonv,vr0, chmom,qa21)
    !***********************************************************************
    ! calculates the perpendicular part of the local moment.
    ! if l_relax is true the angle of the output local moment is calculated
    ! and mixed with the input angles using mix_b as the mixing parameter
    ! if l_constr is true the output constraint b-field is calculated and
    ! mixed with the input contraint field using mix_b
    ! Philipp Kurz 2000-02-09
    !***********************************************************************

    USE m_constants
    USE m_intgr, ONLY : intgr3
    USE m_polangle
    USE m_rotdenmat
    USE m_types
    IMPLICIT NONE
    TYPE(t_noco),INTENT(IN)          :: noco
    TYPE(t_nococonv),INTENT(INOUT)   :: nococonv
    TYPE(t_atoms),INTENT(IN)         :: atoms

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: itype, iRepAtom
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT    (IN) :: chmom(:,:)!(atoms%ntype,input%jspins)
    REAL, INTENT    (IN) :: vr0(:,:,:)!(atoms%jmtd,atoms%ntype,jspd)
    COMPLEX, INTENT (IN) :: qa21(atoms%ntype)
    !     ..
    !     .. Local Scalars ..
    INTEGER iri
    REAL b_xavh,scale,b_con_outx,b_con_outy,mx,my,mz,&
         &     alphh,betah,mz_tmp,mx_mix,my_mix,mz_mix,absmag
    REAL    rho11,rho22, alphdiff
    COMPLEX rho21
    !     ..
    !     .. Local Arrays ..
    REAL b_xc_h(atoms%jmtd),b_xav(atoms%ntype)

    ! angles in nocoinp file are (alph-alphdiff)
    IF (noco%l_ss) THEN
       alphdiff = 2.0*pi_const*(nococonv%qss(1)*atoms%taual(1,iRepAtom) + &
                                nococonv%qss(2)*atoms%taual(2,iRepAtom) + &
                                nococonv%qss(3)*atoms%taual(3,iRepAtom) )
    ELSE
       alphdiff = 0.0
    END IF

    !---> calculated the comp. of the local moment vector
    mx = 2*REAL(qa21(itype))
    my = 2*AIMAG(qa21(itype))
    mz = chmom(itype,1) - chmom(itype,2)
    absmag=SQRT(mx*mx+my*my+mz*mz)
    WRITE  (6,8025) mx,my,mz,absmag
    !---> determine the polar angles of the moment vector in the local frame
    CALL pol_angle(mx,my,mz,betah,alphh)
    WRITE  (6,8026) betah,alphh
8025 FORMAT(2x,'--> local frame: ','mx=',f9.5,' my=',f9.5,' mz=',f9.5,' |m|=',f9.5)
8026 FORMAT(2x,'-->',10x,' local beta=',f9.5,&
         &                   '  local alpha=',f9.5)

    IF(noco%l_alignMT) THEN
      WRITE  (6,8400) nococonv%beta,nococonv%alph
      8400   FORMAT(2x,'-->',10x,'nococonv%beta=',f9.5, ' nococonv%alpha=',f9.5)
    END IF
    IF (noco%l_relax(itype)) THEN
       !--->    rotate the (total (integrated) density matrix to obtain
       !--->    it in the global spin coordinate frame
       rho11 = chmom(itype,1)
       rho22 = chmom(itype,2)
       rho21 = qa21(itype)
       CALL rot_den_mat(nococonv%alph(itype),nococonv%beta(itype), rho11,rho22,rho21)
       !--->    determine the polar angles of the mom. vec. in the global frame
       mx = 2*REAL(rho21)
       my = 2*AIMAG(rho21)
       mz = rho11 - rho22
       IF ( mz .LT. 0.0 ) THEN
          mx = (-1.0) * mx_mix
          my = (-1.0) * my_mix
          mz = (-1.0) * mz_mix
       ENDIF
       CALL pol_angle(mx,my,mz,betah,alphh)
       WRITE  (6,8027) nococonv%beta(itype),nococonv%alph(itype)-alphdiff
       WRITE  (6,8028) betah,alphh-alphdiff
8027   FORMAT(2x,'-->',10x,' input nococonv%beta=',f9.5, '  input nococonv%alpha=',f9.5)
8028   FORMAT(2x,'-->',10x,'output nococonv%beta=',f9.5, ' output nococonv%alpha=',f9.5)

       !  ff    do the same for mixed density: rho21 = mix_b * rho21
       rho11 = chmom(itype,1)
       rho22 = chmom(itype,2)
       rho21 = qa21(itype)
       rho21 = noco%mix_b * rho21
       CALL rot_den_mat(nococonv%alph(itype),nococonv%beta(itype), rho11,rho22,rho21)
       !--->    determine the polar angles of the mom. vec. in the global frame
       mx_mix = 2*REAL(rho21)
       my_mix = 2*AIMAG(rho21)
       mz_mix = rho11 - rho22
       WRITE  (6,8031) mx_mix,my_mix
8031   FORMAT(2x,'--> global frame: ','mixed mx=',f9.5,' mixed my=',f9.5)
       ! if magnetic moment (in local frame!) is negative, direction of quantization
       ! has to be antiparallel!
       mz_tmp = chmom(itype,1) - chmom(itype,2)
       IF ( mz_tmp .LT. 0.0 ) THEN
          mx_mix = (-1.0) * mx_mix
          my_mix = (-1.0) * my_mix
          mz_mix = (-1.0) * mz_mix
       ENDIF
       ! calculate angles alpha and beta in global frame
       CALL pol_angle(mx_mix,my_mix,mz_mix,betah,alphh)
       WRITE  (6,8029) betah,alphh-alphdiff
8029   FORMAT(2x,'-->',10x,' new nococonv%beta  =',f9.5, '  new nococonv%alpha  =',f9.5)
       nococonv%alph(itype) = alphh
       nococonv%beta(itype) = betah
    ENDIF

    IF (noco%l_constr) THEN
       !--->    calculate the average value of B_xc (<B_xc>)
       DO iri = 1,atoms%jri(itype)
          b_xc_h(iri) = (  vr0(iri,itype,1) - vr0(iri,itype,2) )*atoms%rmsh(iri,itype)
       ENDDO
       CALL intgr3(b_xc_h,atoms%rmsh(1,itype),atoms%dx(itype),atoms%jri(itype),b_xavh)
       b_xav(itype) = fpi_const*b_xavh/atoms%volmts(itype)
       !--->    calculate the output constraint B-field (B_con)
       !        take negative of absolute value! gb`05
       scale = -ABS(b_xav(itype)/(chmom(itype,1)-chmom(itype,2)))
       b_con_outx = scale*mx
       b_con_outy = scale*my
       !--->    mix input and output constraint fields
       WRITE  (6,8100) nococonv%b_con(1,itype),nococonv%b_con(2,itype)
       WRITE  (6,8200) b_con_outx,b_con_outy
       nococonv%b_con(1,itype) = nococonv%b_con(1,itype) + noco%mix_b*b_con_outx
       nococonv%b_con(2,itype) = nococonv%b_con(2,itype) + noco%mix_b*b_con_outy
    ENDIF

8100 FORMAT (2x,'-->',10x,' input B_con_x=',f12.6,&
         &                    '  input B_con_y=',f12.6,&
         &                    ' B_xc average=',f12.6)
8200 FORMAT (2x,'-->',10x,' delta B_con_x=',f12.6,&
         &                    ' delta B_con_y=',f12.6)

  END SUBROUTINE m_perp
END MODULE m_m_perp
