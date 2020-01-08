!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vgen_finalize
  USE m_juDFT
   USE m_xcBfield
CONTAINS
  SUBROUTINE vgen_finalize(mpi,oneD,field,cell,atoms,stars,vacuum,sym,noco,input,sphhar,vTot,vCoul,denRot)
    !     ***********************************************************
    !     FLAPW potential generator                           *
    !     ***********************************************************
    !     some rescaling is done here
    !     ***********************************************************
    !     in noco case vmatgen is called to generate 2x2 int-potential
    !     **********************************************************
    USE m_constants
    USE m_vmatgen
    USE m_types
    USE m_rotate_mt_den_tofrom_local
    USE m_sfTests
    USE m_magnMomFromDen
    IMPLICIT NONE
    TYPE(t_mpi),       INTENT(IN)     :: mpi
    TYPE(t_oneD),      INTENT(IN)     :: oneD
    TYPE(t_field),                INTENT(INOUT)  :: field
    TYPE(t_cell),      INTENT(IN)     :: cell
    TYPE(t_vacuum),INTENT(IN)       :: vacuum
    TYPE(t_noco),INTENT(IN)         :: noco
    TYPE(t_sym),INTENT(IN)          :: sym
    TYPE(t_stars),INTENT(IN)        :: stars
    TYPE(t_atoms),INTENT(IN)        :: atoms
    TYPE(t_input),INTENT(IN)        :: input
    TYPE(t_sphhar),INTENT(IN)       :: sphhar
    TYPE(t_potden),INTENT(INOUT)    :: vTot,vCoul,denRot
    !     ..
    !     .. Local Scalars ..
    INTEGER i,js,n

      TYPE(t_potden) :: div, phi, checkdiv
      TYPE(t_potden), DIMENSION(3) :: cvec, corrB, bxc
      REAL                         :: b(3,atoms%ntype),dummy1(atoms%ntype),dummy2(atoms%ntype)



    ! Rescale vTot%pw_w with number of stars
    IF (.NOT.noco%l_noco) THEN
       DO js=1,SIZE(vtot%pw_w,2)
          DO i=1,stars%ng3
             vTot%pw_w(i,js)=vtot%pw_w(i,js) / stars%nstr(i)
          END DO
       END DO
    ELSEIF(noco%l_noco) THEN
       CALL vmatgen(stars,atoms,vacuum,sym,input,denRot,vTot)
       IF (noco%l_mtnocoPot) THEN
          CALL rotate_mt_den_from_local(atoms,sphhar,sym,denRot,noco,vtot)
       END IF
    ENDIF

    ! Source-free testwise
    !CALL sftest(mpi,dimension,field,stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,1,inDen,1.0)
    !CALL sftest(mpi,dimension,field,stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,1,vTot,2.0)



    ! Once it is tested:
    IF (noco%l_mtnocoPot.AND.noco%l_sourceFree) THEN ! l_sf will go here
       CALL magnMomFromDen(input,atoms,noco,vTot,b,dummy1,dummy2)
       DO i=1,atoms%ntype
          WRITE  (6,8025) i,b(1,i),b(2,i),b(3,i),SQRT(b(1,i)**2+b(2,i)**2+b(3,i)**2)
          8025 FORMAT(2x,'--> Bfield before SF (atom ',i2,': ','Bx 1=',f9.5,' By=',f9.5,' Bz=',f9.5,' |B|=',f9.5)
       END DO
       CALL makeVectorField(sym,stars,atoms,sphhar,vacuum,input,noco,vTot,2.0,bxc)
       CALL sourcefree(mpi,field,stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,bxc,div,phi,cvec,corrB,checkdiv)
       CALL div%resetPotDen()
       CALL checkdiv%resetPotDen()
       CALL phi%resetPotDen()
       DO i=1,3
          CALL bxc(i)%resetPotDen()
          CALL corrB(i)%resetPotDen()
       END DO
       CALL correctPot(vTot,cvec)
       CALL magnMomFromDen(input,atoms,noco,vTot,b,dummy1,dummy2)
       DO i=1,atoms%ntype
          WRITE  (6,8026) i,b(1,i),b(2,i),b(3,i),SQRT(b(1,i)**2+b(2,i)**2+b(3,i)**2)
          8026 FORMAT(2x,'--> Bfield after SF (atom ',i2,': ','Bx 1=',f9.5,' By=',f9.5,' Bz=',f9.5,' |B|=',f9.5)
       END DO
    END IF



  !           ---> store v(l=0) component as r*v(l=0)/sqrt(4pi)

 DO js = 1,input%jspins !Used input%jspins instead of SIZE(vtot%mt,4) since the off diag, elements of VTot%mt need no rescaling.
       DO n = 1,atoms%ntype
          vTot%mt(:atoms%jri(n),0,n,js)  = atoms%rmsh(:atoms%jri(n),n)*vTot%mt(:atoms%jri(n),0,n,js)/sfp_const
       ENDDO
    ENDDO     ! js =1,input%jspins

    ! Rescale vCoul%pw_w with number of stars
    DO js = 1, SIZE(vCoul%pw_w,2)
       DO i = 1, stars%ng3
          vcoul%pw_w(i,js) = vcoul%pw_w(i,js) / stars%nstr(i)  !this normalization is needed for gw
       END DO
    END DO

    !Copy first vacuum into second vacuum if this was not calculated before
    IF (vacuum%nvac==1) THEN
       vTot%vacz(:,2,:)  = vTot%vacz(:,1,:)
       IF (sym%invs) THEN
          vTot%vacxy(:,:,2,:)  = CMPLX(vTot%vacxy(:,:,1,:))
       ELSE
          vTot%vacxy(:,:,2,:)  = vTot%vacxy(:,:,1,:)
       ENDIF
    ENDIF

  END SUBROUTINE vgen_finalize
END MODULE m_vgen_finalize
