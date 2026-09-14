!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_qpwtonmt
  !***************************************************************
  !     This subroutine calculates the lattice harmonic expansions
  !     for the plane wave density inside the spheres.
  !
  !             Stefan Bl"ugel  , IFF, Nov. 1997
  !***************************************************************
  IMPLICIT NONE
  PRIVATE

  INTERFACE qpw_to_nmt
     !> qpw_to_nmt_spin:   a single plane wave field added to one spin component
     !>                    of rho (collinear case)
     !> qpw_to_nmt_denmat: several plane wave fields distributed over the
     !>                    components of the muffin-tin density matrix
     MODULE PROCEDURE qpw_to_nmt_spin, qpw_to_nmt_denmat
  END INTERFACE

  PUBLIC :: qpw_to_nmt

CONTAINS

  SUBROUTINE qpw_to_nmt_spin(sphhar,atoms,stars,sym,cell,fmpi,jspin,l_cutoff,qpwc,rho)

    USE m_types

    IMPLICIT NONE

    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_mpi),INTENT(IN)      :: fmpi

    INTEGER, INTENT (IN)    :: jspin,l_cutoff
    COMPLEX, INTENT (IN)    :: qpwc(:)       !(stars%ng3)
    REAL,    INTENT (INOUT) :: rho(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,jspins)

    CALL qpw_to_nmt_denmat(sphhar,atoms,stars,sym,cell,fmpi,l_cutoff,&
                           RESHAPE(qpwc,[SIZE(qpwc),1]),rho(:,:,:,jspin:jspin))

  END SUBROUTINE qpw_to_nmt_spin

  SUBROUTINE qpw_to_nmt_denmat(sphhar,atoms,stars,sym,cell,fmpi,l_cutoff,qpwc,rho,wmt)
    !***************************************************************
    !     qpwc(:,i) holds the plane wave representation of the i-th field. Field
    !     i is added to the component j of the muffin-tin density matrix of atom
    !     type n with the weight wmt(i,j,n); without wmt the fields are added to
    !     the components of the same index.
    !
    !     This is used to add the (spherical) core tails of all atoms to all
    !     muffin-tin spheres in a non-collinear calculation: the plane wave part
    !     is given as the charge density and the magnetization vector in the
    !     GLOBAL frame,
    !        qpwc(:,1) = n , qpwc(:,2:4) = (m_x,m_y,m_z) ,
    !     while the muffin-tin density matrix is stored in the LOCAL frame of
    !     each atom type. The weights then rotate the magnetization vector into
    !     the local frame of the atom type it is added to, see cdnovlp_noco.
    !     Note that all these fields are real in real space, so each component
    !     is accumulated in exactly the same way as in the collinear case.
    !***************************************************************
    USE m_constants
    USE m_phasy1
    USE m_sphbes
    USE m_types

    IMPLICIT NONE

    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_mpi),INTENT(IN)      :: fmpi

    INTEGER, INTENT (IN)    :: l_cutoff
    COMPLEX, INTENT (IN)    :: qpwc(:,:)     !(stars%ng3,nsrc)
    REAL,    INTENT (INOUT) :: rho(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,ndest)
    REAL, OPTIONAL, INTENT (IN) :: wmt(:,:,:) !(nsrc,ndest,atoms%ntype)

    REAL,    PARAMETER :: zero = 0.0
    COMPLEX, PARAMETER :: czero = (0.0,0.0)

    INTEGER :: in,j,jl,j0,jm,k,l,lh,m,n,nd,nrm,n1,n2,na,lm,i,idest,nsrc,ndest
    REAL    :: d0,gr,r0,rr
    COMPLEX :: sm

    COMPLEX :: pylm( (atoms%lmaxd+1)**2,atoms%ntype )
    REAL    :: bsl(0:atoms%lmaxd),bsl2(atoms%jmtd,0:atoms%lmaxd)
    INTEGER :: mr(atoms%ntype),lmx(atoms%ntype),lmxx(atoms%ntype),ntypsy_o(atoms%ntype)
    COMPLEX,ALLOCATABLE :: cp(:),cp0(:)
    REAL,ALLOCATABLE    :: w(:,:,:),rho_tmp(:,:,:,:)

    CALL timestart("qpw_to_nmt")

    nsrc  = SIZE(qpwc,2)
    ndest = SIZE(rho,4)
    ALLOCATE(cp0(ndest),w(nsrc,ndest,atoms%ntype))

    IF (PRESENT(wmt)) THEN
       w = wmt
    ELSE
       w = zero
       DO i = 1, MIN(nsrc,ndest)
          w(i,i,:) = 1.0
       END DO
    END IF

    !----> cut-off l-expansion of non-spherical charge contribution
    !      from coretails of neighboring atom for l> l_cutoff
    DO n = 1,atoms%ntype
       na = atoms%firstAtom(n)
       lmx(n) = MIN( atoms%lmax(n) , l_cutoff )
       ntypsy_o(n) = sym%ntypsy(na)
    END DO
    !
    !----> identify atoms with the same radial mesh
    !
    j0 = 0
    r0 = zero
    d0 = zero
    nrm= 0
    DO n = 1 , atoms%ntype
       IF (.NOT.(atoms%jri(n).EQ.j0 .AND. atoms%rmsh(1,n).EQ.r0 &
            &                          .AND. atoms%dx(n).EQ.d0)) THEN
          j0 = atoms%jri(n)
          r0 = atoms%rmsh(1,n)
          d0 = atoms%dx(n)
          nrm= nrm + 1
          lmxx(nrm) = lmx(n)
       END IF
       mr(nrm)=n
       lmxx(nrm) = MAX( lmx(n) , lmx(nrm) )
    END DO
    !
    !=====> Loop over the g-vectors
    !
    ! ----> g=0 component
    IF (fmpi%irank == 0) THEN
       DO  n = 1 , atoms%ntype
          DO idest = 1, ndest
             cp0(idest) = stars%nstr(1)*SUM(qpwc(1,:)*w(:,idest,n))
          END DO
          DO j = 1,atoms%jri(n)
             rr = atoms%rmsh(j,n)*atoms%rmsh(j,n)
             DO idest = 1, ndest
                rho(j,0,n,idest) = rho(j,0,n,idest) + rr*sfp_const*REAL(cp0(idest))
             END DO
          ENDDO
       ENDDO
    ELSE
       rho = zero
    ENDIF

    ! ----> g.ne.0 components
    !
    !     g=|=0 terms : \sum_{g =|= 0} 4 \pi i^l \rho_int(g) r_i^{2} \times
    !                    j_{l} (gr_i) \exp{ig\xi_i} Y^*_{lm} (g)
    !
    ! ---->     calculate structure constant for each atom
    !     (4pi*i**l/nop(3)*sum(R){exp(iRG(taual-taur)*conjg(ylm(RG))
    !
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP SHARED(fmpi,stars,atoms,sym,cell,sphhar,qpwc,rho,w,ndest,lmx,lmxx,mr,nrm,ntypsy_o) &
    !$OMP PRIVATE(k,pylm,n1,in,n2,j,gr,bsl,bsl2,jl,n,nd,lh,l,sm,jm,m,lm,cp,idest,rho_tmp)
    ALLOCATE(cp(ndest))
    ALLOCATE(rho_tmp(SIZE(rho,1),0:SIZE(rho,2)-1,SIZE(rho,3),ndest))
    rho_tmp = zero
    !$OMP DO
    DO k = fmpi%irank+2,stars%ng3,fmpi%isize
       CALL phasy1(atoms,stars,sym,cell,k,pylm)
       n1 = 1
       DO in = 1 , nrm
          n2 = mr(in)
          bsl2 = zero
          DO j = 1,atoms%jri(n1)
             gr = stars%sk3(k)*atoms%rmsh(j,n1)
             CALL sphbes(lmxx(in),gr,bsl)
             DO jl = 0,lmxx(in)
                bsl2(j,jl) = bsl(jl) * atoms%rmsh(j,n1)*atoms%rmsh(j,n1)
             END DO
          END DO
          DO n = n1,n2
             nd = ntypsy_o(n)
             DO idest = 1, ndest
                cp(idest) = stars%nstr(k)*SUM(qpwc(k,:)*w(:,idest,n))
             END DO
             DO lh = 0,sphhar%nlh(nd)
                l = sphhar%llh(lh,nd)
                IF ( l.GT.lmx(n) ) CYCLE
                sm = czero
                DO jm = 1,sphhar%nmem(lh,nd)
                   m  = sphhar%mlh(jm,lh,nd)
                   lm = l*(l+1) + m + 1
                   sm = sm + CONJG(sphhar%clnu(jm,lh,nd)) *pylm(lm,n)
                END DO
                DO idest = 1, ndest
                   rho_tmp(:atoms%jri(n),lh,n,idest) = rho_tmp(:atoms%jri(n),lh,n,idest) &
                                                     + bsl2(:atoms%jri(n),l) * REAL(cp(idest)*sm)
                END DO
             END DO
          END DO
          n1 = n2 + 1
       ENDDO
    ENDDO
    !$OMP END DO
    !$OMP CRITICAL
    rho = rho + rho_tmp
    !$OMP END CRITICAL
    DEALLOCATE(rho_tmp,cp)
    !$OMP END PARALLEL
    !
    CALL timestop("qpw_to_nmt")
    !
  END SUBROUTINE qpw_to_nmt_denmat
END MODULE m_qpwtonmt
