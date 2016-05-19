MODULE m_vmts
  !     *******************************************************************
  !     this subroutine calculates the lattice harmonics expansion coeffi-*
  !     cients of the coulomb potential for all atom types                *
  !                                c.l.fu, r.podloucky                    *
  !     *******************************************************************
CONTAINS
  SUBROUTINE vmts(mpi,stars,sphhar,atoms,&
       &                sym,cell,oneD,&
       &                vpw,rho,&
       &                vr)

#include"cpp_double.h"
    USE m_constants
    USE m_types
    USE m_intgr, ONLY : intgr2
    USE m_phasy1
    USE m_sphbes
    USE m_od_phasy

    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ...
    TYPE(t_mpi),INTENT(IN)     :: mpi
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_cell),INTENT(IN)    :: cell
    TYPE(t_oneD),INTENT(IN)    :: oneD
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (IN) :: vpw(:,:)!(stars%n3d,input%jspins)
    REAL,    INTENT (IN) :: rho(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,input%jspins)
    REAL,    INTENT (OUT):: vr(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,input%jspins)
    !-odim
    !+odim
    !     ..
    !     .. Local Scalars ..
    COMPLEX cp,sm
    REAL rmt2l,rmtl,ror,rr,rrlr,fpl21
    INTEGER i,jm,k,l,l21,lh ,n,nd ,lm,n1,nat,m
    !     ..
    !     .. Local Arrays ..
    COMPLEX vtl(0:sphhar%nlhd,atoms%ntypd)
    COMPLEX pylm( (atoms%lmaxd+1)**2 ,atoms%ntypd )
    REAL    f1r(atoms%jmtd),f2r(atoms%jmtd),x1r(atoms%jmtd),x2r(atoms%jmtd)
    REAL    sbf(0:atoms%lmaxd),rrl(atoms%jmtd),rrl1(atoms%jmtd)
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER ierr(3)
    COMPLEX, ALLOCATABLE :: c_b(:)
    ! ..
    ! ..  External Subroutines
    EXTERNAL MPI_REDUCE
#endif
    !     ..
    !     ..
    !     ----> calculate lattice harmonics expansion coefficients of the
    !     ----> interstitial coulomb potential on the sphere boundaries
    nat = 1
    DO n = 1,atoms%ntype
       DO lh = 0,sphhar%nlh(atoms%ntypsy(nat))
          vtl(lh,n) = 0.e0
       ENDDO
       nat = nat + atoms%neq(n) 
    ENDDO
#ifdef CPP_MPI
    vtl(:,:) = CMPLX(0.0,0.0)
    CALL MPI_BCAST(vpw,SIZE(vpw),CPP_MPI_COMPLEX,0,&
         &                          mpi,ierr)
#endif

    !           ----> g=0 component
    IF (mpi%irank == 0) THEN
       DO n = 1,atoms%ntype
          vtl(0,n) = sfp_const*vpw(1,1)
       ENDDO
    ENDIF
    !           ----> g.ne.0 components
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(cp,pylm,nat,n,sbf,nd,lh,&
    !$OMP& sm,jm,m,lm,l) REDUCTION(+:vtl)
    DO k = mpi%irank+2, stars%ng3, mpi%isize
       cp = vpw(k,1)*stars%nstr(k)
       IF (.NOT.oneD%odi%d1) THEN
          CALL phasy1(&
               &                  atoms,stars,sym,&
               &                  cell,k,&
               &                  pylm)
       ELSE
          !-odim
CALL od_phasy(&
               &           atoms%ntype,stars%n3d,atoms%nat,atoms%lmaxd,atoms%ntype,atoms%neq,atoms%lmax,&
               &           atoms%taual,cell%bmat,stars%kv3,k,oneD%odi,oneD%ods,&
               &           pylm)
          !+odim
       END IF
       !
       nat = 1
       DO n = 1,atoms%ntype
          CALL sphbes(atoms%lmax(n),stars%sk3(k)*atoms%rmt(n),sbf)
          nd = atoms%ntypsy(nat)
          DO lh = 0,sphhar%nlh(nd)
             l = sphhar%llh(lh,nd)
             sm = (0.,0.)
             DO jm = 1,sphhar%nmem(lh,nd)
                m = sphhar%mlh(jm,lh,nd)
                lm = l*(l+1) + m + 1 
                sm = sm + CONJG(sphhar%clnu(jm,lh,nd))*pylm(lm,n)
             ENDDO
             vtl(lh,n) = vtl(lh,n) + cp*sbf(l)*sm
          ENDDO
          nat = nat + atoms%neq(n)
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
#ifdef CPP_MPI
    n1 = (sphhar%nlhd+1)*atoms%ntypd
    ALLOCATE(c_b(n1))
    CALL MPI_REDUCE(vtl,c_b,n1,CPP_MPI_COMPLEX,MPI_SUM,0, mpi%mpi_comm,ierr)
    IF (mpi%irank.EQ.0) THEN
       vtl=reshape(c_b,(/sphhar%nlhd+1,atoms%ntypd/))
    ENDIF
    DEALLOCATE (c_b)
#endif

    !     ----> solution of the poisson's equation
    nat = 1
    DO  n = 1,atoms%ntype
       nd = atoms%ntypsy(nat)
       DO  lh = 0,sphhar%nlh(nd)
          l = sphhar%llh(lh,nd)
          l21 = 2*l + 1
          fpl21 = fpi_const/l21
          DO i = 1,atoms%jri(n)
             rrl(i) = atoms%rmsh(i,n)**l
             rrl1(i) = 1./( rrl(i) * atoms%rmsh(i,n) )
             x1r(i) = rrl(i)*rho(i,lh,n,1)
             x2r(i) = rrl1(i)*rho(i,lh,n,1)
          ENDDO
          CALL intgr2(x1r,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),f1r)
          CALL intgr2(x2r,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),f2r)
          rmtl = 1./atoms%rmt(n)**l
          rmt2l = 1./atoms%rmt(n)**l21
          DO  i = 1,atoms%jri(n)
             rrlr = rrl(i)*rmt2l
             ror = rrl(i)*rmtl
             vr(i,lh,n,1) = fpl21 * (rrl1(i)*f1r(i)-rrlr*f1r(atoms%jri(n))+&
                  &                   rrl(i) * (f2r(atoms%jri(n))-f2r(i))) + ror*vtl(lh,n)
          ENDDO
       ENDDO
       nat = nat + atoms%neq(n)
    ENDDO
    DO  n = 1,atoms%ntype
       DO  i = 1,atoms%jri(n)
          rr = atoms%rmsh(i,n)/atoms%rmt(n)
          vr(i,0,n,1) = vr(i,0,n,1)-sfp_const*(1.-rr)/atoms%rmsh(i,n)*atoms%zatom(n)
       ENDDO
    ENDDO
  END SUBROUTINE vmts
END MODULE m_vmts
