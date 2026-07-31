!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_etabinit
  USE m_juDFT
  !     *******************************************************
  !     *****   set up etab via old core program          *****
  !     *******************************************************
  !     modified to run with core-levels as provided by setcor
  !     ntab & ltab transport this info to core.F        gb`02
  !------------------------------------------------------------
CONTAINS
  SUBROUTINE etabinit(atoms, input, iType, vr, etab, ntab, ltab, nkmust)

    USE m_constants
    USE m_differ
    USE m_types
    IMPLICIT NONE

    TYPE(t_atoms), INTENT(IN)  :: atoms
    TYPE(t_input), INTENT(IN)  :: input
    INTEGER,       INTENT(IN)  :: iType
    REAL,          INTENT(IN)  :: vr(atoms%jmtd)
    REAL,          INTENT(OUT) :: etab(100,atoms%ntype)
    INTEGER,       INTENT(OUT) :: ntab(100,atoms%ntype), ltab(100,atoms%ntype)
    INTEGER,       INTENT(OUT) :: nkmust(atoms%ntype)

    REAL    :: c,d,dxx,e,fj,fl,fn,rn,rnot,t2,z,t1,rr,weight
    INTEGER :: i,ic,ipos,j,korb,nst,ncmsh,ierr
    INTEGER :: kappa(maxval(atoms%econf%num_states)), nprnc(maxval(atoms%econf%num_states))
    REAL    :: eig(maxval(atoms%econf%num_states)), occ(maxval(atoms%econf%num_states),1)
    REAL    :: vrd(atoms%msh), a(atoms%msh), b(atoms%msh)

    c     = c_light(1.0)
    ncmsh = atoms%msh

    z    = atoms%zatom(iType)
    dxx  = atoms%dx(iType)
    rnot = atoms%rmsh(1,iType)
    d    = EXP(atoms%dx(iType))
    rn   = rnot * (d**(ncmsh-1))

    CALL atoms%econf(iType)%get_core(nst, nprnc, kappa, occ)

    WRITE (oUnit,FMT=8000) z, rnot, dxx, atoms%jri(iType)

    DO j = 1, atoms%jri(iType)
       vrd(j) = vr(j)
    END DO
    IF (input%l_core_confpot) THEN
       t1 = 0.125
       t2 = vrd(atoms%jri(iType)) / atoms%rmt(iType) - atoms%rmt(iType)*t1
       rr = atoms%rmt(iType)
    ELSE
       t2 = vrd(atoms%jri(iType)) / (atoms%jri(iType) - atoms%msh)
    END IF
    IF (atoms%jri(iType) < atoms%msh) THEN
       DO i = atoms%jri(iType) + 1, atoms%msh
          IF (input%l_core_confpot) THEN
             rr     = d*rr
             vrd(i) = rr * (t2 + rr*t1)
          ELSE
             vrd(i) = vrd(atoms%jri(iType)) + t2 * (i - atoms%jri(iType))
          END IF
       END DO
    END IF

    nst = atoms%econf(iType)%num_core_states
    DO korb = 1, nst
       fn     = nprnc(korb)
       fj     = iabs(kappa(korb)) - 0.5e0
       weight = 2*fj + 1.e0
       fl     = fj + 0.5e0 * isign(1, kappa(korb))
       e      = -2.0 * (z/(fn+fl))**2
       CALL differ(fn, fl, fj, c, z, dxx, rnot, rn, d, atoms%msh, vrd, e, a, b, ierr)
       IF (ierr /= 0) CALL juDFT_error("error in core-levels", calledby="etabinit")
       WRITE (oUnit,FMT=8010) fn, fl, fj, e, weight
       eig(korb) = e
    END DO

    ic = 0
    DO korb = 1, nst
       fn     = nprnc(korb)
       fj     = iabs(kappa(korb)) - 0.5e0
       weight = 2*fj + 1.e0
       fl     = fj + 0.5e0 * isign(1, kappa(korb))
       DO i = 1, INT(weight)
          ic = ic + 1
          IF (kappa(korb) > 0) THEN
             ipos = ic + 1 + i
          ELSE IF (kappa(korb) < -1) THEN
             ipos = ic - 2*(iabs(kappa(korb))-1) + MAX(i-2, 0)
          ELSE
             ipos = ic
          END IF
          etab(ipos,iType) = eig(korb)
          ntab(ipos,iType) = NINT(fn)
          ltab(ipos,iType) = NINT(fl)
       END DO
    END DO
    nkmust(iType) = ic

    DO i = 1, nkmust(iType)
       WRITE(oUnit,'(f12.6,2i3)') etab(i,iType), ntab(i,iType), ltab(i,iType)
    END DO

8000 FORMAT (/,/,10x,'z=',f4.0,5x,'r(1)=',e14.6,5x,'dx=',f8.6,5x,&
                'm.t.index=',i4,/,15x,'n',4x,'l',5x,'j',4x,'energy',7x,'weight')
8010 FORMAT (12x,2f5.0,f6.1,f10.4,f12.4)

  END SUBROUTINE etabinit
END MODULE m_etabinit
