!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vbfield
   USE m_juDFT
CONTAINS
  SUBROUTINE dfpt_vbfield(input,starsq,noco,atoms,sym,sphhar,cell,vTot,vTotIm)
    !This subroutine calculates the Zeeman field perturbation
    USE m_types
    USE m_constants
    USE m_sphbes
    USE m_phasy1
    USE m_dfpt_gradient

    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)::input
    TYPE(t_stars),INTENT(IN) :: starsq
    TYPE(t_noco),INTENT(IN) ::noco
    TYPE(t_atoms),INTENT(IN)::atoms
    TYPE(t_sym),INTENT(IN)  ::sym
    TYPE(t_sphhar),INTENT(IN)::sphhar
    TYPE(t_cell),INTENT(IN) ::cell
    TYPE(t_potden),INTENT(INOUT)::vTot,vTotIm

    INTEGER :: iType, iSpin, l, m, lm, ll1, i, imax, lmax
    REAL    :: bsign, qabs
    LOGICAL :: l_afm

    COMPLEX, ALLOCATABLE :: vbfield_mt(:,:,:,:)
    REAL,    ALLOCATABLE :: resultreal(:,:,:,:), resultimag(:,:,:,:)
    REAL,    ALLOCATABLE :: sbf(:,:)
    COMPLEX              :: pylm((atoms%lmaxd+1)**2, atoms%ntype)

    l_afm = .FALSE.

    IF (input%jspins.NE.2) CALL judft_error("B-fields can only be used in spin-polarized calculations")
    !IF (noco%l_noco) CALL judft_error("B-fields not implemented in noco case")

    qabs = starsq%sk3(1) ! |q|, analogue to qlim for efield

    vTot%pw(:,:)     = 0.0
    vTot%mt(:,:,:,:) = 0.0
    vTotIm%mt(:,:,:,:) = 0.0

    IF (l_afm) THEN
       IF (qabs.GT.1e-8) CALL judft_error("AFM B-field only implemented for q = 0",calledby="dfpt_vbfield")
       !No IR contribution; alternate sign between sublattices in MT
       DO iType = 1, atoms%ntype
          bsign = MERGE(-1.0, 1.0, MOD(iType,2) == 1)
          vTot%mt(:atoms%jri(iType),0,iType,1) = -bsign * sfp_const/2.0
          vTot%mt(:atoms%jri(iType),0,iType,2) =  bsign * sfp_const/2.0
       END DO
    ELSE
       !interstitial region: G=0 star of starsq is e^{iq.r}
       vTot%pw(1,1) = -1/2.
       vTot%pw(1,2) = +1/2.

       !MT-region: e^{iq.r} = sum_lm 4pi i^l j_l(|q|r) e^{iq.tau} Y*_lm(q) Y_lm(r)
       !At q = 0 this reduces to sfp_const in l=0 (j_0(0)=1, j_l(0)=0 for l>0).
       ALLOCATE(vbfield_mt(atoms%jmtd,(atoms%lmaxd+1)**2,atoms%nat,1))
       vbfield_mt = CMPLX(0.0,0.0)
       ALLOCATE(resultreal(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))
       ALLOCATE(resultimag(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))

       CALL phasy1(atoms, starsq, sym, cell, 1, pylm)

       DO iType = 1, atoms%ntype
          lmax = atoms%lmax(iType)
          imax = atoms%jri(iType)
          ALLOCATE(sbf(0:lmax,imax))
          DO i = 1, imax
             CALL sphbes(lmax,qabs*atoms%rmsh(i,iType),sbf(:,i))
          END DO
          DO l = 0, lmax
             ll1 = l*(l+1)+1
             DO m = -l, l
                lm = ll1 + m
                vbfield_mt(1:imax,lm,atoms%firstAtom(iType),1) = sbf(l,1:imax)*pylm(lm,iType)
             END DO
          END DO
          DEALLOCATE(sbf)
       END DO

       !go to lattice harmonics:
       CALL sh_to_lh(sym, atoms, sphhar, 1, 3, vbfield_mt, resultreal, resultimag)

       DO iSpin = 1, 2
          bsign = MERGE(-1.0, 1.0, iSpin == 1)
          vTot%mt(:,:,:,iSpin)   = bsign * resultreal(:,:,:,1)/2.0
          vTotIm%mt(:,:,:,iSpin) = bsign * resultimag(:,:,:,1)/2.0
       END DO
    END IF



  END SUBROUTINE dfpt_vbfield
END MODULE m_dfpt_vbfield