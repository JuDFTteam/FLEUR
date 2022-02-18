!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_m_perp
IMPLICIT NONE
CONTAINS
  SUBROUTINE m_perp(atoms,itype,noco,nococonv, chmom,qa21,vtot)
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
    INTEGER, INTENT (IN) :: itype
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT    (IN) :: chmom(:)!(atoms%ntype,input%jspins)
    COMPLEX, INTENT (IN) :: qa21
    TYPE(T_potden), INTENT    (IN),OPTIONAL :: vtot
    !     ..
    !     .. Local Scalars ..
    INTEGER iri, iRepAtom
    REAL b_xavh,scale,b_con_outx,b_con_outy,mx,my,mz,&
         &     alphh,betah,mz_tmp,mx_mix,my_mix,mz_mix,absmag
    REAL    rho11,rho22, alphdiff
    COMPLEX rho21
    real :: magmom(0:3)
    !     ..
    !     .. Local Arrays ..
    REAL b_xc_h(atoms%jmtd),b_xav(atoms%ntype)

    irepAtom=sum(atoms%neq(:itype-1))+1
    ! angles in nocoinp file are (alph-alphdiff)
     IF (noco%l_ss) THEN
       alphdiff = 2.0*pi_const*(nococonv%qss(1)*atoms%taual(1,iRepAtom) + &
                                nococonv%qss(2)*atoms%taual(2,iRepAtom) + &
                                nococonv%qss(3)*atoms%taual(3,iRepAtom) )
     ELSE
       alphdiff = 0.0
     END IF

     !---> calculated the comp. of the local moment vector
     magmom=nococonv%denmat_to_mag(chmom(1),chmom(2),qa21)
     !---> determine the polar angles of the moment vector in the local frame
     CALL pol_angle(magmom(1),magmom(2),magmom(3),betah,alphh,.true.)
     call priv_output_moments(itype,magmom(1:3),betah,alphh,global=.false.)

     rho11=chmom(1)
     rho22=chmom(2)
     rho21=qa21
     call nococonv%rotdenmat(itype, rho11,rho22,rho21, toGlobal=.true.)
     magmom=nococonv%denmat_to_mag(rho11,rho22,rho21)
     !now also give output in global frame
     !call nococonv%rot_magvec(itype,magmom,toGlobal=.true.)
     CALL pol_angle(magmom(1),magmom(2),magmom(3),betah,alphh,.true.)
     call priv_output_moments(itype,magmom(1:3),betah,alphh,global=.true.)
     



    IF(noco%l_alignMT(itype)) THEN
      WRITE  (oUnit,8400) itype,nococonv%beta(itype),nococonv%alph(itype)
      8400   FORMAT(2x,'Atom:',I9.1,' -->',10x,'nococonv%beta=',f9.5, ' nococonv%alpha=',f9.5)
    END IF

    IF (noco%l_alignMT(itype).and..not.noco%l_unrestrictMT(itype)) THEN
       !--->    rotate the (total (integrated) density matrix to obtain
       !--->    it in the global spin coordinate frame
       rho11 = chmom(1)
       rho22 = chmom(2)
       rho21 = qa21
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
       CALL pol_angle(mx,my,mz,betah,alphh,.true.)
       WRITE  (oUnit,8027) nococonv%beta(itype),nococonv%alph(itype)-alphdiff
       WRITE  (oUnit,8028) betah,alphh-alphdiff
8027   FORMAT(2x,'-->',10x,' input nococonv%beta=',f9.5, '  input nococonv%alpha=',f9.5)
8028   FORMAT(2x,'-->',10x,'output nococonv%beta=',f9.5, ' output nococonv%alpha=',f9.5)

       !  ff    do the same for mixed density: rho21 = mix_b * rho21
       rho11 = chmom(1)
       rho22 = chmom(2)
       rho21 = qa21
       rho21 = noco%mix_b * rho21
       CALL rot_den_mat(nococonv%alph(itype),nococonv%beta(itype), rho11,rho22,rho21)
       !--->    determine the polar angles of the mom. vec. in the global frame
       mx_mix = 2*REAL(rho21)
       my_mix = 2*AIMAG(rho21)
       mz_mix = rho11 - rho22
       WRITE  (oUnit,8031) mx_mix,my_mix
8031   FORMAT(2x,'--> global frame: ','mixed mx=',f9.5,' mixed my=',f9.5)
       ! if magnetic moment (in local frame!) is negative, direction of quantization
       ! has to be antiparallel!
       mz_tmp = chmom(1) - chmom(2)
       IF ( mz_tmp .LT. 0.0 ) THEN
          mx_mix = (-1.0) * mx_mix
          my_mix = (-1.0) * my_mix
          mz_mix = (-1.0) * mz_mix
       ENDIF
       ! calculate angles alpha and beta in global frame
       CALL pol_angle(mx_mix,my_mix,mz_mix,betah,alphh)
       WRITE  (oUnit,8029) betah,alphh-alphdiff
8029   FORMAT(2x,'-->',10x,' new nococonv%beta  =',f9.5, '  new nococonv%alpha  =',f9.5)
       nococonv%alph(itype) = alphh
       nococonv%beta(itype) = betah
    ENDIF

    IF (noco%l_constrained(itype).and.present(vtot)) THEN
       !--->    calculate the average value of B_xc (<B_xc>)
       DO iri = 1,atoms%jri(itype)
          b_xc_h(iri) = (  vtot%mt(iri,0,itype,1) - vtot%mt(iri,0,itype,2) )*atoms%rmsh(iri,itype)
       ENDDO
       CALL intgr3(b_xc_h,atoms%rmsh(1,itype),atoms%dx(itype),atoms%jri(itype),b_xavh)
       b_xav(itype) = fpi_const*b_xavh/atoms%volmts(itype)
       !--->    calculate the output constraint B-field (B_con)
       !        take negative of absolute value! gb`05
       scale = -ABS(b_xav(itype)/(chmom(1)-chmom(2)))
       b_con_outx = scale*mx
       b_con_outy = scale*my
       !--->    mix input and output constraint fields
       WRITE  (oUnit,8100) nococonv%b_con(1,itype),nococonv%b_con(2,itype)
       WRITE  (oUnit,8200) b_con_outx,b_con_outy
       nococonv%b_con(1,itype) = nococonv%b_con(1,itype) + noco%mix_b*b_con_outx
       nococonv%b_con(2,itype) = nococonv%b_con(2,itype) + noco%mix_b*b_con_outy
    ENDIF

8100 FORMAT (2x,'-->',10x,' input B_con_x=',f12.6,&
         &                    '  input B_con_y=',f12.6)
8200 FORMAT (2x,'-->',10x,' delta B_con_x=',f12.6,&
         &                    ' delta B_con_y=',f12.6)

  END SUBROUTINE 

  subroutine priv_output_moments(itype,magmom,beta,alpha,global)
   USE m_xmlOutput
   USE m_constants
   integer,intent(in):: itype 
   real,intent(in)   :: magmom(3),beta,alpha
   logical,intent(in) :: global

   character(len=15):: label 
   character(len=30):: attributes(2)      

   if (global) then 
      WRITE  (oUnit,8125) itype,magmom(1),magmom(2),magmom(3)
      WRITE  (oUnit,8126) itype,beta,alpha
   else
      WRITE  (oUnit,8025) itype,magmom(1),magmom(2),magmom(3),sqrt(dot_product(magmom(1:3),magmom(1:3)))
      WRITE  (oUnit,8026) itype,beta,alpha
   endif
8125 FORMAT(2x,'Atom:',I9.1,' --> global frame: ','mx=',f9.5,' my=',f9.5,' mz=',f9.5,' |m|=',f9.5)
8126 FORMAT(2x,'Atom:',I9.1,' -->',10x,' global beta=',f9.5,'  global alpha=',f9.5)
8025 FORMAT(2x,'Atom:',I9.1,' --> local frame: ','mx=',f9.5,' my=',f9.5,' mz=',f9.5,' |m|=',f9.5)
8026 FORMAT(2x,'Atom:',I9.1,' -->',10x,' local beta=',f9.5,'  local alpha=',f9.5)

   WRITE(attributes(1),'(i0)') iType
   WRITE(attributes(2),'(3(f9.7,1x))') magmom(1),magmom(2),magmom(3)
   label=trim(merge("globalMagMoment","localMagMoment ",global))
   CALL writeXMLElementFormPoly(label,(/'atomType','vec     '/),&
                             attributes,reshape((/8,3,6,30/),(/2,2/)))


   end subroutine
END MODULE m_m_perp
