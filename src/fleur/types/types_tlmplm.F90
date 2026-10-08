!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_tlmplm
   !! Local muffin-tin Hamiltonian in the unified radial basis of t_radfun.
   !! Basis index of type n: u(lm) = lm for l<=lrange(n), then udot(lm) of the LAPW channels in the
   !! order of the matching coefficients (atoms%udot_rows), udot(lm) of the APW channels, and the
   !! (2l+1) functions of each LO in atoms%llo order; see ind.
   use m_types_rsoc
   use m_types_radfun
   use m_judft
   use m_types_atoms
  IMPLICIT NONE
  PRIVATE

  TYPE t_tlmplm
     COMPLEX,ALLOCATABLE :: h(:,:,:,:,:)            !(0:nbasd-1,0:nbasd-1,ntype,ispin,jspin)
     COMPLEX,ALLOCATABLE :: h_loc_nonsph(:,:,:,:,:) !non-spherical block of u and LAPW udot (l<=lnonsph) incl. LDA+U; Cholesky factor if spin-diagonal
     INTEGER,ALLOCATABLE :: ind(:,:,:)              !position of (slot,lm) of type n in h (slot,0:lmd,ntype); -1 if not in basis
     INTEGER,ALLOCATABLE :: nbas(:)                 !size of the unified basis per type
     INTEGER,ALLOCATABLE :: lrange(:)               !largest l of the LAPW part per type
     REAL,ALLOCATABLE    :: e_shift(:,:)
     TYPE(t_radfun),ALLOCATABLE :: radfun(:)        !radial basis of h
     REAL,ALLOCATABLE    :: h_sph_extra(:,:,:,:)    !MetaGGA: lambda=0 elements (uu,ud,du,dd) for l>lnonsph, added by hsmt_sph (0:3,0:lmaxd,ntype,jspins)
     LOGICAL             :: l_sph_extra = .FALSE.
     TYPE(t_rsoc)        :: rsoc
   CONTAINS
     PROCEDURE,PASS :: init => tlmplm_init
  END TYPE t_tlmplm
  PUBLIC t_tlmplm
CONTAINS
  SUBROUTINE tlmplm_init(td,atoms,jspins,l_fulllmax)
    !! l_fulllmax: LAPW part up to lmax (forces) instead of lnonsph
    CLASS(t_tlmplm),INTENT(INOUT):: td
    TYPE(t_atoms),INTENT(IN)     :: atoms
    INTEGER,INTENT(in)           :: jspins
    LOGICAL,INTENT(IN)           :: l_fulllmax

    INTEGER :: n,l,m,lo,s,nb,sns
    INTEGER :: boff(0:atoms%lmaxd)

    IF (ALLOCATED(td%h)) DEALLOCATE(td%h,td%h_loc_nonsph,td%ind,td%nbas,td%lrange,td%e_shift,td%radfun,td%h_sph_extra)
    td%lrange = MERGE(atoms%lmax,atoms%lnonsph,l_fulllmax)
    ALLOCATE(td%nbas(atoms%ntype))
    ALLOCATE(td%ind(MAXVAL([(atoms%num_radial_functions_per_l(n),n=1,atoms%ntype)]),0:atoms%lmaxd*(atoms%lmaxd+2),atoms%ntype),source=-1)
    DO n = 1,atoms%ntype
       IF (ANY(atoms%llo(:atoms%nlo(n),n)>td%lrange(n))) &
          CALL judft_error("Local orbital with l larger than the non-spherical cutoff lnonsph",calledby="types_tlmplm")
       s = td%lrange(n)*(td%lrange(n)+2)+1
       ! udot of the LAPW channels in the order of the matching coefficients, then udot of APW channels
       boff(0:td%lrange(n)) = atoms%udot_rows(td%lrange(n),n)
       nb = atoms%num_ab_rows(td%lrange(n),n)
       DO l = 0,td%lrange(n)
          DO m = -l,l
             td%ind(1,l*(l+1)+m,n) = l*(l+1)+m
             IF (boff(l)>=0) td%ind(2,l*(l+1)+m,n) = boff(l)+l+m
          END DO
          IF (boff(l)>=0) CYCLE
          DO m = -l,l
             td%ind(2,l*(l+1)+m,n) = nb+l+m
          END DO
          nb = nb+2*l+1
       END DO
       DO lo = 1,atoms%nlo(n)
          IF (atoms%l_dulo(lo,n)) CYCLE ! APW LO lives on the udot slot
          l = atoms%llo(lo,n)
          DO m = -l,l
             td%ind(atoms%slot_of_lo(lo,n),l*(l+1)+m,n) = nb+l+m
          END DO
          nb = nb+2*l+1
       END DO
       td%nbas(n) = nb
    END DO
    sns = MAXVAL(atoms%lnonsph*(atoms%lnonsph+2)+1)
    ALLOCATE(td%h(0:MAXVAL(td%nbas)-1,0:MAXVAL(td%nbas)-1,atoms%ntype,jspins,jspins),source=CMPLX(0.0,0.0))
    ALLOCATE(td%h_loc_nonsph(0:2*sns-1,0:2*sns-1,atoms%ntype,jspins,jspins),source=CMPLX(0.0,0.0))
    ALLOCATE(td%e_shift(atoms%ntype,jspins),source=0.0)
    ALLOCATE(td%radfun(atoms%ntype))
    ALLOCATE(td%h_sph_extra(0:3,0:atoms%lmaxd,atoms%ntype,jspins),source=0.0)
    td%l_sph_extra = .FALSE.
  END SUBROUTINE tlmplm_init

END MODULE m_types_tlmplm
