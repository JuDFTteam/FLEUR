MODULE m_psqpw
  !     ***********************************************************
  !     generates the fourier coefficients of pseudo charge density
  !                                   c.l.fu
  !         corrected april 1990   m.w.
  !
  ! cf. M.Weinert J.Math.Phys. 22(11) (1981) p.2434 eq. (10)-(15)
  !
  !
  !     parallelized 04/08 gb
  !     ***********************************************************

CONTAINS

  SUBROUTINE psqpw( mpi, atoms, sphhar, stars, vacuum, DIMENSION, cell, input, sym, oneD, &
       &     qpw, rho, rht, l_xyav, yukawa_residual, psq )

#include"cpp_double.h"
    USE m_constants
    USE m_phasy1
    USE m_mpmom 
    USE m_sphbes
    USE m_qsf
    USE m_od_phasy
    USE m_od_cylbes
    USE m_types
    use m_DoubleFactorial
    use m_SphBessel
    IMPLICIT NONE

    TYPE(t_mpi),        INTENT(IN)  :: mpi
    TYPE(t_atoms),      INTENT(IN)  :: atoms
    TYPE(t_sphhar),     INTENT(IN)  :: sphhar
    TYPE(t_stars),      INTENT(IN)  :: stars
    TYPE(t_vacuum),     INTENT(IN)  :: vacuum
    TYPE(t_dimension),  INTENT(IN)  :: DIMENSION
    TYPE(t_cell),       INTENT(IN)  :: cell
    TYPE(t_input),      INTENT(IN)  :: input
    TYPE(t_sym),        INTENT(IN)  :: sym
    TYPE(t_oneD),       INTENT(IN)  :: oneD
    LOGICAL,            INTENT(IN)  :: l_xyav
    COMPLEX,            INTENT(IN)  :: qpw(stars%ng3) 
    REAL,               INTENT(IN)  :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype) 
    REAL,               INTENT(IN)  :: rht(vacuum%nmzd,2)
    logical,            intent(in)  :: yukawa_residual
    COMPLEX,            INTENT(OUT) :: psq(stars%ng3)

    COMPLEX                         :: psint, sa, sl, sm
    REAL                            :: f, fact, fpo, gz, p, qvac, rmtl, s, fJ, gr, g
    INTEGER                         :: ivac, k, l, n, n1, nc, ncvn, lm, ll1, nd, m, nz
    COMPLEX                         :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%ntype)
    COMPLEX                         :: qlm(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    REAL                            :: q2(vacuum%nmzd)
    real                            :: pn(0:atoms%lmaxd,atoms%ntype)
    real                            :: aj(0:atoms%lmaxd+DIMENSION%ncvd+1)
    REAL                            :: rht1(vacuum%nmz)
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER                         :: ierr(3)
    COMPLEX, ALLOCATABLE            :: c_b(:)
#endif

    ! Calculate multipole moments
    CALL mpmom( mpi, atoms, sphhar, stars, sym, cell, oneD, qpw, rho, yukawa_residual, qlm )
#ifdef CPP_MPI
    psq(:) = CMPLX(0.0,0.0)
    CALL MPI_BCAST(qpw,size(qpw),CPP_MPI_COMPLEX,0,mpi%mpi_comm,ierr)
    nd = (2*atoms%lmaxd+1)*(atoms%lmaxd+1)*atoms%ntype
    CALL MPI_BCAST(qlm,nd,CPP_MPI_COMPLEX,0,mpi%MPI_COMM,ierr)
#endif
    !
    ! pn(l,n) = (2l + 2nc(n) + 3)!! / (2l + 1)!! R^l  ;   ncv(n)=n+l in paper
    ! 
    DO n = 1,atoms%ntype
       rmtl = 1.0
       DO l = 0,atoms%lmax(n)
          IF (l.GE.atoms%ncv(n)) THEN
             pn(l,n) = 0.0
          ELSE
             p = 1.
             DO nc = l,atoms%ncv(n)
                p = p* (2*nc+3)
             ENDDO
             pn(l,n) = p/rmtl
          END IF
          rmtl = rmtl*atoms%rmt(n)
       ENDDO
    ENDDO
    !
    ! G eq 0 term (eq.29) : \tilde \rho_s (0) = \sqrt{4 pi} / \Omega \sum_i \tilde q_{00}^i
    !
    s = 0.
    DO n = 1,atoms%ntype
       s = s + atoms%neq(n)*REAL(qlm(0,0,n))
    ENDDO
    IF (mpi%irank == 0) THEN
       psq(1) = qpw(1) + (sfp_const/cell%omtil)*s
    ENDIF
    !
    ! G ne 0 term (eq.28) : \tilde \rho_s (K) = 4 pi / \Omega \sum_{lmi} (-i)^l \exp{-iK\xi_i}
    !                    (2n+3)!!/(2l+1)!! * 1/R_i^l * j_{n+1}(KR_i)/(KR_i)^{n-l+1} Y_{lm} (K)
    !
    fpo = 1./cell%omtil
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pylm,sa,n,ncvn,aj,&
    !$OMP&  sl,l,n1,ll1,sm,m,lm)
    DO k = mpi%irank+2, stars%ng3, mpi%isize
       IF (.NOT.oneD%odi%d1) THEN
          CALL phasy1(&
               &               atoms,stars,sym,&
               &               cell,k,&
               &               pylm)
       ELSE
          !-odim
          CALL od_phasy(&
               &        atoms%ntype,stars%ng3,atoms%nat,atoms%lmaxd,atoms%ntype,atoms%neq,atoms%lmax,&
               &        atoms%taual,cell%bmat,stars%kv3,k,oneD%odi,oneD%ods,&
               &        pylm)
          !+odim
       END IF
       !
       sa = 0.
       DO n = 1,atoms%ntype
          ncvn = atoms%ncv(n)
          CALL sphbes(ncvn+1,stars%sk3(k)*atoms%rmt(n),aj)
          sl = 0.
          DO l = 0,atoms%lmax(n)
             IF (l.GE.ncvn) GO TO 60
             n1 = ncvn - l + 1
             ll1 = l*(l+1) + 1
             sm = 0.
             DO m = -l,l
                lm = ll1 + m 
                sm = sm + qlm(m,l,n)*CONJG(pylm(lm,n))
             ENDDO
60           sl = sl + pn(l,n)/ ((stars%sk3(k)*atoms%rmt(n))**n1)*aj(ncvn+1)*sm
          ENDDO
          sa = sa + atoms%neq(n)*sl
       ENDDO
       psq(k) = qpw(k) + fpo*sa
    ENDDO
    !$OMP END PARALLEL DO
#ifdef CPP_MPI
    ALLOCATE(c_b(stars%ng3))
    CALL MPI_REDUCE(psq,c_b,stars%ng3,CPP_MPI_COMPLEX,MPI_SUM,0,mpi%MPI_COMM,ierr)
    IF (mpi%irank.EQ.0) THEN
       psq(:stars%ng3)=c_b(:stars%ng3)
    ENDIF
    DEALLOCATE (c_b)
#endif
    IF (mpi%irank == 0) THEN
       !
       ! Check: integral of the pseudo charge density within the slab
       !
       IF (input%film .AND. .NOT.oneD%odi%d1) THEN
          psint = psq(1)*stars%nstr(1)*vacuum%dvac
          DO k = 2,stars%ng3
             IF (stars%ig2(k).EQ.1) THEN
                gz = stars%kv3(3,k)*cell%bmat(3,3)
                f = 2.*SIN(gz*cell%z1)/gz
                psint = psint + stars%nstr(k)*psq(k)*f
             END IF
          ENDDO
          psint = cell%area*psint
       ELSEIF (input%film .AND. oneD%odi%d1) THEN
          !-odim
          psint = (0.0,0.0)
          DO k = 2,stars%ng3
             IF (stars%kv3(3,k).EQ.0) THEN
                g = (stars%kv3(1,k)*cell%bmat(1,1) + stars%kv3(2,k)*cell%bmat(2,1))**2 +&
                     &             (stars%kv3(1,k)*cell%bmat(1,2) + stars%kv3(2,k)*cell%bmat(2,2))**2
                gr = SQRT(g)
                CALL od_cylbes(1,gr*cell%z1,fJ)
                f = 2*cell%vol*fJ/(gr*cell%z1)
                psint = psint + stars%nstr(k)*psq(k)*f
             END IF
          ENDDO
          psint = psint + psq(1)*stars%nstr(1)*cell%vol
          !+odim
       ELSEIF (.NOT.input%film) THEN
          psint = psq(1)*stars%nstr(1)*cell%omtil
       ENDIF
       WRITE (6,FMT=8000) psint
       WRITE (16,FMT=8000) psint
8000   FORMAT (/,10x,'integral of pseudo charge density inside the slab='&
            &       ,5x,2f11.6)
       IF (.NOT.input%film) RETURN
       !
       ! Normalized pseudo density
       !
       IF (.NOT.oneD%odi%d1) THEN
          qvac = 0.0
          DO ivac = 1,vacuum%nvac
             CALL qsf(vacuum%delz,rht(1,ivac),q2,vacuum%nmz,0)
             q2(1) = q2(1)*cell%area
             qvac = qvac + q2(1)*2./REAL(vacuum%nvac)
          ENDDO
          qvac = qvac - 2*input%sigma
       ELSE
          !-odim
          qvac = 0.0
          DO nz = 1,vacuum%nmz
             rht1(nz) = (cell%z1+(nz-1)*vacuum%delz)*&
                  &           rht(nz,vacuum%nvac)
          ENDDO
          CALL qsf(vacuum%delz,rht1(1),q2,vacuum%nmz,0)
          qvac = cell%area*q2(1)
          !+odim
       END IF
       !      fact = abs(qvac/psint)
       !      DO k = 1,nq3
       !         psq(k) = fact*psq(k)
       !      ENDDO
       IF (l_xyav) RETURN
       fact = (qvac + psint)/(stars%nstr(1)*cell%vol)
       psq(1) = psq(1) - fact
       !-gu
       WRITE (6,FMT=8010) fact*1000
       WRITE (16,FMT=8010) fact*1000
8010   FORMAT (/,10x,'                     1000 * normalization const. ='&
            &       ,5x,2f11.6)
       !
    ENDIF ! mpi%irank == 0 

  END SUBROUTINE psqpw

END MODULE m_psqpw
