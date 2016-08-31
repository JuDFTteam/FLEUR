!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_m_perp
CONTAINS 
  SUBROUTINE m_perp(atoms,itype,noco,vr0, chmom,qa21,alphdiff)
    !***********************************************************************
    ! calculates the perpendicular part of the local moment.
    ! if l_relax is true the angle of the output local moment is calculated
    ! and mixed with the input angles using mix_b as the mixing parameter
    ! if l_constr is true the output constraint b-field is calculated and
    ! mixed with the input contraint field using mix_b
    ! Philipp Kurz 2000-02-09
    !***********************************************************************

    USE m_intgr, ONLY : intgr3
    USE m_constants, ONLY : fpi_const
    USE m_polangle
    USE m_rotdenmat
    USE m_types
    IMPLICIT NONE
    TYPE(t_noco),INTENT(INOUT)   :: noco
    TYPE(t_atoms),INTENT(IN)     :: atoms

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: itype
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT    (IN) :: chmom(:,:)!(atoms%ntypd,dimension%jspd)
    REAL, INTENT    (IN) :: alphdiff(atoms%ntypd) 
    REAL, INTENT    (IN) :: vr0(:,:,:)!(atoms%jmtd,atoms%ntypd,jspd)
    COMPLEX, INTENT (IN) :: qa21(atoms%ntypd)
    !     ..
    !     .. Local Scalars ..
    INTEGER iri
    REAL b_xavh,scale,b_con_outx,b_con_outy,mx,my,mz,&
         &     alphh,betah,mz_tmp,mx_mix,my_mix,mz_mix
    REAL    rho11,rho22
    COMPLEX rho21
    !     ..
    !     .. Local Arrays ..
    REAL b_xc_h(atoms%jmtd),b_xav(atoms%ntypd)


    !---> calculated the comp. of the local moment vector
    mx = 2*REAL(qa21(itype))
    my = 2*AIMAG(qa21(itype))
    mz = chmom(itype,1) - chmom(itype,2)
    WRITE  (6,8025) mx,my
    WRITE (16,8025) mx,my
    !---> determine the polar angles of the moment vector in the local frame
    CALL pol_angle(mx,my,mz,betah,alphh)
    WRITE  (6,8026) betah,alphh
    WRITE (16,8026) betah,alphh
8025 FORMAT(2x,'--> local frame: ','mx=',f9.5,' my=',f9.5)
8026 FORMAT(2x,'-->',10x,' delta beta=',f9.5,&
         &                   '  delta alpha=',f9.5)

    IF (noco%l_relax(itype)) THEN
       !--->    rotate the (total (integrated) density matrix to obtain
       !--->    it in the global spin coordinate frame
       rho11 = chmom(itype,1)
       rho22 = chmom(itype,2)
       rho21 = qa21(itype)
       CALL rot_den_mat(noco%alph(itype),noco%beta(itype), rho11,rho22,rho21)
       !--->    determine the polar angles of the mom. vec. in the global frame
       mx = 2*REAL(rho21)
       my = 2*AIMAG(rho21)
       mz = rho11 - rho22
       CALL pol_angle(mx,my,mz,betah,alphh)
       WRITE  (6,8027) noco%beta(itype),noco%alph(itype)-alphdiff(itype)
       WRITE (16,8027) noco%beta(itype),noco%alph(itype)-alphdiff(itype)
       WRITE  (6,8028) betah,alphh-alphdiff(itype)
       WRITE (16,8028) betah,alphh-alphdiff(itype)
8027   FORMAT(2x,'-->',10x,' input noco%beta=',f9.5, '  input noco%alpha=',f9.5)
8028   FORMAT(2x,'-->',10x,'output noco%beta=',f9.5, ' output noco%alpha=',f9.5)

       !  ff    do the same for mixed density: rho21 = mix_b * rho21
       rho11 = chmom(itype,1)
       rho22 = chmom(itype,2)
       rho21 = qa21(itype)
       rho21 = noco%mix_b * rho21
       CALL rot_den_mat(noco%alph(itype),noco%beta(itype), rho11,rho22,rho21)
       !--->    determine the polar angles of the mom. vec. in the global frame
       mx_mix = 2*REAL(rho21)
       my_mix = 2*AIMAG(rho21)
       mz_mix = rho11 - rho22
       WRITE  (6,8031) mx_mix,my_mix
       WRITE (16,8031) mx_mix,my_mix 
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
       WRITE  (6,8029) betah,alphh-alphdiff(itype)
       WRITE (16,8029) betah,alphh-alphdiff(itype)
8029   FORMAT(2x,'-->',10x,' new noco%beta  =',f9.5, '  new noco%alpha  =',f9.5)
       noco%alph(itype) = alphh
       noco%beta(itype) = betah
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
       WRITE  (6,8100) noco%b_con(1,itype),noco%b_con(2,itype)
       WRITE (16,8100) noco%b_con(1,itype),noco%b_con(2,itype)
       WRITE  (6,8200) b_con_outx,b_con_outy
       WRITE (16,8200) b_con_outx,b_con_outy
       noco%b_con(1,itype) = noco%b_con(1,itype) + noco%mix_b*b_con_outx
       noco%b_con(2,itype) = noco%b_con(2,itype) + noco%mix_b*b_con_outy
    ENDIF

8100 FORMAT (2x,'-->',10x,' input B_con_x=',f12.6,&
         &                    '  input B_con_y=',f12.6,&
         &                    ' B_xc average=',f12.6)
8200 FORMAT (2x,'-->',10x,' delta B_con_x=',f12.6,&
         &                    ' delta B_con_y=',f12.6)

  END SUBROUTINE m_perp
END MODULE m_m_perp
