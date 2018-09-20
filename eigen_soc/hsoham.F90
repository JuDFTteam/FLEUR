MODULE m_hsoham
  !
  !*********************************************************************
  ! set up spin-orbit contribution to hamiltonian
  !*********************************************************************
  !
CONTAINS
  SUBROUTINE hsoham(&
       atoms,noco,input,nsz,chelp,rsoc,ahelp,bhelp,&
       nat_start,nat_stop,n_rank,n_size,SUB_COMM,&
       hsomtx)

#include"cpp_double.h"

    USE m_types
    IMPLICIT NONE
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER ierr(3)
#endif
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_noco),INTENT(IN)    :: noco
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_rsoc),INTENT(IN)    :: rsoc
    !     ..
    !     .. Scalar Arguments ..
    !     ..
    INTEGER, INTENT (IN) ::  nat_start,nat_stop,n_rank,n_size,SUB_COMM
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nsz(:)!(dimension%jspd)  
    COMPLEX, INTENT (IN) :: ahelp(:,:,:,:)!(lmd,nat_l,dimension%neigd,dimension%jspd)
    COMPLEX, INTENT (IN) :: bhelp(:,:,:,:)!(lmd,nat_l,dimension%neigd,dimension%jspd)
    COMPLEX, INTENT (IN) :: chelp(-atoms%llod :,:,:,:,:)!(-llod:llod ,dimension%neigd,atoms%nlod,nat_l,dimension%jspd)
    COMPLEX, INTENT (OUT):: hsomtx(:,:,:,:)!(dimension%neigd,neigd,2,2)
    !     ..
    !     .. Local Scalars ..
    COMPLEX c_1,c_2,c_3,c_4,c_5
    INTEGER i,j,jsp,jsp1,l,lwn ,m1,n,na,nn,i1,j1,ilo,ilop,m,nat_l,na_g,lm,ll1
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: c_b(:,:,:),c_a(:,:,:),c_c(:,:,:),c_buf(:)
    !     ..
    !
    !---------------------------------------------------------------------
    !  ss'  _
    ! H  = \  (xhelp(s,i,na,l,m) conjg(yhelp(s',j,na,l,m')*rsoxy(na,l,s,s')
    !           *<slm|L*S|s'lm'>
    !  ij  /_
    !       na,l,m,m'
    !                       x,y = a,b
    !---------------------------------------------------------------------
    !
    nat_l = nat_stop - nat_start + 1  ! atoms processed by this pe
    !
    !---> update hamiltonian matrices: upper triangle
    !
    DO i1 = 1,2
       jsp = i1
       IF (input%jspins.EQ.1) jsp = 1
       DO j1 = 1,2
          jsp1 = j1
          IF (input%jspins.EQ.1) jsp1 = 1
          !!$OMP PARALLEL DEFAULT(none)&
          !!$OMP PRIVATE(j,na,na_g,n,nn,l,m,m1,ilo,i,lwn,ilop)& 
          !!$OMP PRIVATE(c_a,c_b,c_c,c_1,c_2,c_3,c_4,c_5) &
          !!$OMP SHARED(hsomtx,i1,jsp,j1,jsp1,nsz,atoms)& 
          !!$OMP SHARED(ahelp,bhelp,chelp,noco,nat_start,nat_stop,nat_l)&
          !!$OMP SHARED(rsoc)

          ALLOCATE ( c_b(-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,nat_l),&
                     c_a(-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,nat_l),&
                     c_c(-atoms%llod :atoms%llod ,atoms%nlod ,nat_l) )

          !!$OMP DO 

          DO j = 1,nsz(jsp1)
             !
             ! prepare \sum_m' conjg( xhelp(m',l,na,j,jsp1) ) * soangl(l,m,i1,l,m',j1)
             !
             na = 0 ; na_g = 0
             DO n = 1,atoms%ntype
                DO nn = 1, atoms%neq(n)
                   na_g = na_g + 1
                   IF ((na_g.GE.nat_start).AND.(na_g.LE.nat_stop)) THEN
                      na = na + 1
                      !--> regular part
                      DO l = 1,atoms%lmax(n)
                         ll1 = l*(l+1) 
                         DO m = -l,l
                            c_a(m,l,na) = CMPLX(0.,0.)
                            c_b(m,l,na) = CMPLX(0.,0.)
                            DO m1 = -l,l
                               lm = ll1 + m1
                               c_a(m,l,na) = c_a(m,l,na) + rsoc%soangl(l,m,i1,l,m1,j1)&
                                    *CONJG(ahelp(lm,na,j,jsp1))
                               c_b(m,l,na) = c_b(m,l,na) + rsoc%soangl(l,m,i1,l,m1,j1)&
                                    *CONJG(bhelp(lm,na,j,jsp1))
                            ENDDO
                         ENDDO
                      ENDDO
                      !--> LO contribution
                      DO ilo = 1,atoms%nlo(n)
                         l = atoms%llo(ilo,n)
                         IF (l.GT.0) THEN
                            DO m = -l,l
                               c_c(m,ilo,na) = CMPLX(0.,0.)
                               DO m1 = -l,l
                                  c_c(m,ilo,na) = c_c(m,ilo,na) + CONJG(&
                                       chelp(m1,j,ilo,na,jsp1))*rsoc%soangl(l,m,i1,l,m1,j1)
                               ENDDO
                            ENDDO
                         ENDIF
                      ENDDO
                      ! end lo's
                   ENDIF
                ENDDO  ! nn
             ENDDO     ! n
                !
             ! continue loop structure
             !
             DO i = 1,nsz(jsp)
                hsomtx(i,j,i1,j1) = CMPLX(0.,0.)
                na = 0 ; na_g = 0
                !
                !--->    loop over each atom type
                !
                DO n = 1,atoms%ntype
                   lwn = atoms%lmax(n)
                   !
                   !--->    loop over equivalent atoms
                   !
                   DO  nn = 1,atoms%neq(n)
                      na_g = na_g + 1
                      IF ((na_g.GE.nat_start).AND.(na_g.LE.nat_stop)) THEN
                         na = na + 1
                         DO l = 1,lwn
                            ll1 = l*(l+1) 
                            DO m = -l,l
                               lm = ll1 + m
                               c_1 =   rsoc%rsopp(n,l,i1,j1) * ahelp(lm,na,i,jsp) +&
                                      rsoc%rsopdp(n,l,i1,j1) * bhelp(lm,na,i,jsp)
                               c_2 =  rsoc%rsoppd(n,l,i1,j1) * ahelp(lm,na,i,jsp) +&
                                     rsoc%rsopdpd(n,l,i1,j1) * bhelp(lm,na,i,jsp)
                               hsomtx(i,j,i1,j1) = hsomtx(i,j,i1,j1) +&
                                    c_1*c_a(m,l,na) + c_2*c_b(m,l,na)  
                            ENDDO
                            ! 
                         ENDDO
                         !--> LO contribution
                         DO ilo = 1,atoms%nlo(n)
                            l = atoms%llo(ilo,n)
                            ll1 = l*(l+1)
                            IF (l.GT.0) THEN
                               DO m = -l,l
                                  lm = ll1 + m
                                  c_3 = rsoc%rsopplo(n,ilo,i1,j1) *ahelp(lm,na,i,jsp) +&
                                       rsoc%rsopdplo(n,ilo,i1,j1) *bhelp(lm,na,i,jsp)
                                  c_4 = rsoc%rsoplop(n,ilo,i1,j1) *chelp(m,i,ilo,na,jsp)
                                  c_5 =rsoc%rsoplopd(n,ilo,i1,j1) *chelp(m,i,ilo,na,jsp)
                                  hsomtx(i,j,i1,j1) = hsomtx(i,j,i1,j1) + &
                                       c_4*c_a(m,l,na) + c_5*c_b(m,l,na) +&
                                       c_3*c_c(m,ilo,na)
                               ENDDO
                               DO ilop = 1,atoms%nlo(n)
                                  IF (atoms%llo(ilop,n).EQ.l) THEN
                                     DO m = -l,l
                                        hsomtx(i,j,i1,j1) = hsomtx(i,j,i1,j1) + &
                                             rsoc%rsoploplop(n,ilop,ilo,i1,j1) * &
                                             chelp(m,i,ilop,na,jsp) * c_c(m,ilo,na)
                                     ENDDO
                                  ENDIF
                               ENDDO
                            ENDIF
                         ENDDO
                         ! end lo's
                      ENDIF
                   ENDDO
                ENDDO ! atoms
             ENDDO
             !!i
          ENDDO
          !!j
          !!$OMP END DO
          DEALLOCATE (c_a,c_b,c_c)
          !!$OMP END PARALLEL
       ENDDO
       !!jsp1
    ENDDO
    !!jsp
    !
    !---> update hamiltonian matrices: lower triangle
    !
    DO i = 1,nsz(1)
       DO j = 1,nsz(input%jspins)
          hsomtx(j,i,2,1) = CONJG(hsomtx(i,j,1,2))
       ENDDO
    ENDDO
    !
#ifdef CPP_MPI
    CALL MPI_BARRIER(SUB_COMM,ierr)
    n = 4*nsz(1)*nsz(input%jspins)
    ALLOCATE(c_buf(n))
    CALL MPI_REDUCE(hsomtx(1,1,1,1),c_buf,n,CPP_MPI_COMPLEX,MPI_SUM,0,SUB_COMM,ierr)
    IF (n_rank.EQ.0) THEN
        CALL CPP_BLAS_ccopy(n, c_buf, 1, hsomtx(1,1,1,1), 1)
    ENDIF
    DEALLOCATE(c_buf)
#endif
    !
    RETURN
  END SUBROUTINE hsoham
END MODULE m_hsoham
