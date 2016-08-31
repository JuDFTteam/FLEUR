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
CONTAINS
  SUBROUTINE qpw_to_nmt(&
       &                      sphhar,atoms,stars,&
       &                      sym,cell,oneD,mpi,&
       &                      jspin,l_cutoff,qpwc,&
       &                      rho)
    !
    USE m_constants
    USE m_phasy1
    USE m_sphbes
    USE m_types
    USE m_od_phasy

    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_mpi),INTENT(IN)     :: mpi

    INTEGER, INTENT (IN) :: jspin,l_cutoff    
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (IN) :: qpwc(stars%n3d)
    REAL,    INTENT (INOUT) :: rho(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)
    !-odim
    !+odim
    !     ..
    !     .. Local Scalars ..
    REAL,    PARAMETER :: zero = 0.0
    COMPLEX, PARAMETER :: czero = (0.0,0.0)
    INTEGER in,j,jl,j0,jm,k,l,lh,m,n,nd,nrm,n1,n2,na,lm
    REAL    d0,gr,r0,rr
    COMPLEX cp,sm,cprr2
    !     ..
    !     .. Local Arrays ..
    COMPLEX pylm( (atoms%lmaxd+1)**2,atoms%ntypd ) !,bsl_c(atoms%jmtd,0:lmaxd)
    REAL    bsl(0:atoms%lmaxd),bsl_r(atoms%jmtd,0:atoms%lmaxd),bsl_i(atoms%jmtd,0:atoms%lmaxd)
    INTEGER mr(atoms%ntypd),lmx(atoms%ntypd),lmxx(atoms%ntypd),ntypsy_o(atoms%ntypd)
    !     ..
    !$      REAL,ALLOCATABLE        :: rho_tmp(:,:,:)


    !----> cut-off l-expansion of non-spherical charge contribution
    !      from coretails of neighboring atom for l> l_cutoff
    !
    na = 1
    DO n = 1,atoms%ntype
       lmx(n) = MIN( atoms%lmax(n) , l_cutoff )
       ntypsy_o(n) = atoms%ntypsy(na)
       na = na + atoms%neq(n)
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
    CALL timestart("qpw_to_nmt")
    IF (mpi%irank == 0) THEN
       cp = qpwc(1)*stars%nstr(1)
       DO  n = 1 , atoms%ntype
          DO j = 1,atoms%jri(n)
             rr = atoms%rmsh(j,n)*atoms%rmsh(j,n)
             rho(j,0,n,jspin) = rho(j,0,n,jspin) + rr*sfp_const*cp
          ENDDO
       ENDDO
    ELSE
       rho(:,:,:,jspin) = 0.0
    ENDIF

    ! ----> g.ne.0 components
    !
    !     g=|=0 terms : \sum_{g =|= 0} 4 \pi i^l \rho_int(g) r_i^{2} \times
    !                    j_{l} (gr_i) \exp{ig\xi_i} Y^*_{lm} (g)
    !
    ! ---->     calculate structure constant for each atom
    !     (4pi*i**l/nop(3)*sum(R){exp(iRG(taual-taur)*conjg(ylm(RG)) 
    !
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,cp,pylm,n1,in,n2,j)      &
    !$OMP  PRIVATE(cprr2,gr,bsl,jl,bsl_r,bsl_i,n,nd,lh,l,sm,jm,m,lm)  &
    !$OMP  PRIVATE(rho_tmp)


    !$    ALLOCATE(rho_tmp(size(rho,1),0:size(rho,2)-1,size(rho,3)))
    !$    rho_tmp=0
    !$OMP DO
    DO k = mpi%irank+2,stars%ng3,mpi%isize
       cp = qpwc(k)*stars%nstr(k)
       IF (.NOT.oneD%odi%d1) THEN
          CALL phasy1(&
               &                  atoms,stars,sym,&
               &                  cell,k,&
               &                  pylm)
       ELSE
          !-odim
          CALL od_phasy(&
               &              atoms%ntypd,stars%n3d,atoms%natd,atoms%lmaxd,atoms%ntype,atoms%neq,atoms%lmax,&
               &              atoms%taual,cell%bmat,stars%kv3,k,oneD%odi,oneD%ods,&
               &              pylm) !keep
          !+odim
       END IF
       !
       n1 = 1
       DO in = 1 , nrm
          n2 = mr(in)
          DO j = 1,atoms%jri(n1)
             cprr2 = cp*atoms%rmsh(j,n1)*atoms%rmsh(j,n1)
             gr = stars%sk3(k)*atoms%rmsh(j,n1)
             CALL sphbes(lmxx(in),gr,bsl)
             DO jl=0,lmxx(in)
                bsl_r(j,jl) = bsl(jl) *  REAL(cprr2)
                bsl_i(j,jl) = bsl(jl) * AIMAG(cprr2)
             END DO
          END DO
          DO n = n1,n2 
             nd = ntypsy_o(n)
             DO lh = 0,sphhar%nlh(nd)
                l = sphhar%llh(lh,nd)
                IF ( l.LE.lmx(n) ) THEN
                   sm = czero
                   DO jm = 1,sphhar%nmem(lh,nd)
                      m  = sphhar%mlh(jm,lh,nd)
                      lm = l*(l+1) + m + 1
                      sm = sm + CONJG(sphhar%clnu(jm,lh,nd))&
                           &                              *pylm(lm,n)
                   END DO
                   !$                if (.false.) THEN
                   rho(:,lh,n,jspin) = rho(:,lh,n,jspin)&
                        &                            + bsl_r(:,l) *  REAL(sm)
                   rho(:,lh,n,jspin) = rho(:,lh,n,jspin)&
                        &                            - bsl_i(:,l) * AIMAG(sm)
                   !$                 endif
                   !$            rho_tmp(:,lh,n)=rho_tmp(:,lh,n)+bsl_r(:,l)*real(sm)
                   !$            rho_tmp(:,lh,n)=rho_tmp(:,lh,n)-bsl_i(:,l)*aimag(sm)
                END IF
             END DO
          END DO
          n1 = n2 + 1 
       ENDDO
    ENDDO
    !$OMP END DO
    !$OMP CRITICAL
    !$          rho(:,:,:,jspin)=rho(:,:,:,jspin)+rho_tmp
    !$OMP END CRITICAL
    !$      DEALLOCATE(rho_tmp)
    !$OMP END PARALLEL
    !
    CALL timestop("qpw_to_nmt")
    !
  END SUBROUTINE qpw_to_nmt
END MODULE m_qpwtonmt
