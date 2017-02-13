MODULE m_mpmom
  !     ***********************************************************
  !     determine the multipole moments of (original charge minus
  !     plane wave charge) for each atom type
  !                                     c.l.fu
  ! cf. M.Weinert J.Math.Phys. 22(11) (1981) p.2434 eq. (10)-(15)
  !
  !     qlmo(m,l,n) : mult.mom. of the mufftn-tin charge density
  !     qlmp(m,l,n) : mult.mom. of the plane-wave charge density
  !     qlm (m,l,n) : (output) difference of the former quantities
  !     
  !     ***********************************************************
CONTAINS
  SUBROUTINE mpmom(mpi,atoms,sphhar,stars,&
       &                 sym,cell,oneD,&
       &                 qpw,rho,&
       &                 qlm)

#include"cpp_double.h"
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)         :: mpi
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_stars),INTENT(IN)       :: stars
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_sphhar),INTENT(IN)      :: sphhar
    TYPE(t_atoms),INTENT(IN)       :: atoms
    !
    !
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: rho(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd)
    COMPLEX, INTENT (IN) :: qpw(:,:)     !(stars%ng3,dimension%jspd) 
    COMPLEX, INTENT (OUT):: qlm(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    !-odim
    !+odim
    !     ..
    !     .. Local Scalars ..

    INTEGER j,jm,lh ,mb,mem,mems,n,nd ,l,nat,m

    !     ..
    !     .. Local Arrays ..

    COMPLEX qlmo(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    COMPLEX qlmp(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)


    !
    !     multipole moments of original charge (q_{lm}^i)
    !
    IF (mpi%irank == 0) THEN
       CALL mt_moments(atoms,sphhar,&
            &              rho(:,:,:,1),qlmo)
    ENDIF ! mpi%irank == 0

    CALL pw_moments(mpi,stars,atoms,cell,&
         &                      sym,oneD,&
         &                      qpw(:,1),qlmp)
    !
    ! eq.(15): \tilde q_(lm}^i = q_{lm}^i - q_{lm}^{Ii}
    !
    IF (mpi%irank == 0) THEN
       qlm=qlmo-qlmp

       !
       ! Output section
       !
       nat = 1
       DO  n = 1,atoms%ntype
          WRITE (6,FMT=8000) n
          nd = atoms%ntypsy(nat)
          DO lh = 0,sphhar%nlh(nd)
             l = sphhar%llh(lh,nd)
             mems = sphhar%nmem(lh,nd)
             DO mem = 1,mems
                m = sphhar%mlh(mem,lh,nd)
                WRITE (6,FMT=8010) l,m,qlmo(m,l,n),qlmp(m,l,n)
                !     write(16,1002) l,m,qlmo(m,l,n),qlmp(m,l,n)
             ENDDO
          ENDDO
          nat = nat + atoms%neq(n)
       ENDDO
       !
8000   FORMAT (/,10x,'multipole moments for atom type=',i5,/,/,t3,'l',t7,&
            &       'm',t27,'original',t57,'plane wave')
8010   FORMAT (1x,i2,2x,i3,2x,2 (5x,2e15.5))
       !
    ENDIF ! mpi%irank == 0

  END SUBROUTINE mpmom


  SUBROUTINE mt_moments(atoms,sphhar,&
       &              rho,qlmo)
    !multipole moments of original charge (q_{lm}^i)
    USE m_intgr, ONLY : intgr3
    USE m_constants,ONLY:sfp_const
    USE m_types
    IMPLICIT NONE
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    REAL,INTENT(IN)   :: rho(: ,0:, :) 
    COMPLEX,INTENT(OUT)::qlmo(-atoms%lmaxd:,0:,:)

    INTEGER :: n,ns,jm,nl,l,j,mb ,m,nat
    REAL    :: fint
    REAL f(MAXVAL(atoms%jri))

    qlmo=0.0
    nat = 1
    DO n = 1, atoms%ntype
       ns = atoms%ntypsy(nat)
       jm = atoms%jri(n)
       DO nl = 0, sphhar%nlh(ns)
          l = sphhar%llh(nl,ns)
          DO j = 1, jm
             f(j) = (atoms%rmsh(j,n)**l)*rho(j,nl,n)
          ENDDO
          CALL intgr3(f,atoms%rmsh(:,n),atoms%dx(n),jm,fint)
          DO mb = 1, sphhar%nmem(nl,ns)
             m = sphhar%mlh(mb,nl,ns)
             qlmo(m,l,n) = qlmo(m,l,n) + sphhar%clnu(mb,nl,ns)*fint
          ENDDO
       ENDDO
       qlmo(0,0,n) = qlmo(0,0,n) - atoms%zatom(n)/sfp_const
       nat = nat + atoms%neq(n)
    ENDDO
  END SUBROUTINE mt_moments


  SUBROUTINE pw_moments(mpi,stars,atoms,cell,&
       &                      sym,oneD,&
       &                      qpw,qlmp)
    !multipole moments of plane wave charge inside the spheres (q_{lm}^{Ii})
    USE m_phasy1
    USE m_sphbes
    USE m_od_phasy
    USE m_constants,ONLY:sfp_const
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)   :: mpi
    TYPE(t_oneD),INTENT(IN)  :: oneD
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_stars),INTENT(IN) :: stars
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_atoms),INTENT(IN) :: atoms
    COMPLEX,INTENT(IN) :: qpw(:)
    COMPLEX,INTENT(OUT):: qlmp(-atoms%lmaxd:,0:,:)
    !locals
    INTEGER:: n,k,l,ll1 ,lm,ierr(3),m
    COMPLEX sk3i,cil,nqpw
    COMPLEX pylm( (MAXVAL(atoms%lmax)+1)**2 ,atoms%ntype )
    REAL::sk3r,rl3
    REAL aj(0:MAXVAL(atoms%lmax)+1)
    COMPLEX, ALLOCATABLE :: c_b(:)
    LOGICAL :: od
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    EXTERNAL MPI_REDUCE
#endif
    qlmp= 0.0
    IF (mpi%irank==0) THEN
       !g eq 0 term : \sqrt{4 \pi}/3 R_i^3 \rho_I(0) \delta_{l,0}

       DO n = 1,atoms%ntype
          qlmp(0,0,n) = qpw(1)*stars%nstr(1)*atoms%volmts(n)/sfp_const
       ENDDO
    ENDIF
#ifdef CPP_MPI
    CALL MPI_BCAST(qpw,SIZE(qpw),CPP_MPI_COMPLEX,0,&
         &                          mpi%mpi_comm,ierr)
#endif
    !      g ne 0 terms : \sum_{K \= 0} 4 \pi i^l \rho_I(K) R_i^{l+3} \times
    !      j_{l+1} (KR_i) / KR_i \exp{iK\xi_i} Y^*_{lm} (K)
    od=oneD%odi%d1
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pylm,nqpw,n,sk3r,aj,rl3,sk3i,&
    !$OMP& l,cil,ll1,m,lm) REDUCTION(+:qlmp)
    DO k = mpi%irank+2, stars%ng3, mpi%isize
       IF (od) THEN
          CALL od_phasy(&
               &           atoms%ntype,stars%ng3,atoms%nat,atoms%lmaxd,atoms%ntype,atoms%neq,atoms%lmax,&
               &           atoms%taual,cell%bmat,stars%kv3,k,oneD%odi,oneD%ods,&
               &           pylm)
       ELSE
          CALL phasy1(&
               &               atoms,stars,sym,&
               &               cell,k,&
               &               pylm)
       END IF
       !
       nqpw = qpw(k)*stars%nstr(k)
       DO n = 1,atoms%ntype
          sk3r = stars%sk3(k)*atoms%rmt(n)
          CALL sphbes(atoms%lmax(n)+1,sk3r,aj)
          rl3 = atoms%rmt(n)**3
          sk3i = nqpw/sk3r
          DO l = 0,atoms%lmax(n)
             cil = aj(l+1)*sk3i*rl3
             ll1 = l*(l+1) + 1
             DO m = -l,l
                lm = ll1 + m 
                qlmp(m,l,n) = qlmp(m,l,n) + cil*pylm(lm,n)
             ENDDO
             rl3 = rl3*atoms%rmt(n)
          ENDDO                 ! l = 0, atoms%lmax(n)
       ENDDO                    ! n = 1, atoms%ntype
    ENDDO                       ! k = 2, stars%ng3
    !$OMP END PARALLEL DO
#ifdef CPP_MPI
    n = SIZE(qlmp)
    ALLOCATE(c_b(n))
    CALL MPI_REDUCE(qlmp,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0,&
         &                                   mpi%mpi_comm,ierr)
    IF (mpi%irank.EQ.0) THEN
       qlmp=RESHAPE(c_b,(/SIZE(qlmp,1),SIZE(qlmp,2),SIZE(qlmp,3)/))
    ENDIF
    DEALLOCATE (c_b)
#endif

  END SUBROUTINE pw_moments
END MODULE m_mpmom
