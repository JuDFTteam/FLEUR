MODULE m_hsoham
  !
  !*********************************************************************
  ! set up spin-orbit contribution to hamiltonian
  !*********************************************************************
  !
CONTAINS
  SUBROUTINE hsoham(&
       atoms,noco,input,nsz,chelp,&
       rsoc,ahelp,bhelp,&
       hsomtx)
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_noco),INTENT(IN)    :: noco
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_rsoc),INTENT(IN)    :: rsoc
    !     ..
    !     .. Scalar Arguments ..
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nsz(:)!(dimension%jspd)  
    COMPLEX, INTENT (IN) :: ahelp(-atoms%lmaxd:,:,:,:,:)!(-lmaxd:lmaxd,lmaxd,atoms%nat,dimension%neigd,dimension%jspd)
    COMPLEX, INTENT (IN) :: bhelp(-atoms%lmaxd:,:,:,:,:)!(-lmaxd:lmaxd,lmaxd,atoms%nat,dimension%neigd,dimension%jspd)
    COMPLEX, INTENT (IN) :: chelp(-atoms%llod :,:,:,:,:)!(-llod:llod ,dimension%neigd,atoms%nlod,atoms%nat ,dimension%jspd)
    COMPLEX, INTENT (OUT):: hsomtx(:,:,:,:)!(2,2,dimension%neigd,neigd)
    !     ..
    !     .. Local Scalars ..
    COMPLEX c_1,c_2,c_3,c_4,c_5
    INTEGER i,j,jsp,jsp1,l,lwn ,m1,n,na,nn,i1,j1,ilo,ilop,m
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: c_b(:,:,:),c_a(:,:,:),c_c(:,:,:)
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
    !---> update hamiltonian matrices: upper triangle
    !

    DO i1 = 1,2
       jsp = i1
       IF (input%jspins.EQ.1) jsp = 1
       DO j1 = 1,2
          jsp1 = j1
          IF (input%jspins.EQ.1) jsp1 = 1
          !$OMP PARALLEL DEFAULT(none)&
          !$OMP PRIVATE(j,na,n,nn,l,m,m1,ilo,i,lwn,ilop)& 
          !$OMP PRIVATE(c_a,c_b,c_c,c_1,c_2,c_3,c_4,c_5) &
          !$OMP SHARED(hsomtx,i1,jsp,j1,jsp1,nsz,atoms)& 
          !$OMP SHARED(ahelp,bhelp,chelp,noco)&
          !$OMP SHARED(rsoc)

          ALLOCATE ( c_b(-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,atoms%nat),&
               c_a(-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,atoms%nat),&
               c_c(-atoms%llod :atoms%llod ,atoms%nlod ,atoms%nat) )

          !$OMP DO 

          DO j = 1,nsz(jsp1)
             !
             ! prepare \sum_m' conjg( xhelp(m',l,na,j,jsp1) ) * soangl(l,m,i1,l,m',j1)
             !
             na = 0
             DO n = 1,atoms%ntype
                DO nn = 1, atoms%neq(n)
                   na = na + 1
                   !--> regular part
                   DO l = 1,atoms%lmax(n)
                      DO m = -l,l
                         c_a(m,l,na) = CMPLX(0.,0.)
                         c_b(m,l,na) = CMPLX(0.,0.)
                         DO m1 = -l,l
                            c_a(m,l,na) = c_a(m,l,na) + rsoc%soangl(l,m,i1,l,m1,j1)&
                                 *CONJG(ahelp(m1,l,na,j,jsp1))
                            c_b(m,l,na) = c_b(m,l,na) + rsoc%soangl(l,m,i1,l,m1,j1)&
                                 *CONJG(bhelp(m1,l,na,j,jsp1))
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
                ENDDO
             ENDDO
             !
             ! continue loop structure
             !
             DO i = 1,nsz(jsp)
                hsomtx(i1,j1,i,j) = CMPLX(0.,0.)
                na = 0
                !
                !--->    loop over each atom type
                !
                DO n = 1,atoms%ntype
                      lwn = atoms%lmax(n)
                      !
                      !--->    loop over equivalent atoms
                      !
                      DO  nn = 1,atoms%neq(n)
                         na = na + 1
                         DO l = 1,lwn
                            ! 
                            DO m = -l,l
                               c_1 =   rsoc%rsopp(n,l,i1,j1) * ahelp(m,l,na,i,jsp) +&
                                    rsoc%rsopdp(n,l,i1,j1) * bhelp(m,l,na,i,jsp)
                               c_2 =  rsoc%rsoppd(n,l,i1,j1) * ahelp(m,l,na,i,jsp) +&
                                    rsoc%rsopdpd(n,l,i1,j1) * bhelp(m,l,na,i,jsp)
                               hsomtx(i1,j1,i,j) = hsomtx(i1,j1,i,j) +&
                                    c_1*c_a(m,l,na) + c_2*c_b(m,l,na)  
                            ENDDO
                            ! 
                         ENDDO
                         !--> LO contribution
                         DO ilo = 1,atoms%nlo(n)
                            l = atoms%llo(ilo,n)
                            IF (l.GT.0) THEN
                               DO m = -l,l
                                  c_3 = rsoc%rsopplo(n,ilo,i1,j1) *ahelp(m,l,na,i,jsp) +&
                                       rsoc%rsopdplo(n,ilo,i1,j1) *bhelp(m,l,na,i,jsp)
                                  c_4 = rsoc%rsoplop(n,ilo,i1,j1) *chelp(m,i,ilo,na,jsp)
                                  c_5 =rsoc%rsoplopd(n,ilo,i1,j1) *chelp(m,i,ilo,na,jsp)
                                  hsomtx(i1,j1,i,j) = hsomtx(i1,j1,i,j) + &
                                       c_4*c_a(m,l,na) + c_5*c_b(m,l,na) +&
                                       c_3*c_c(m,ilo,na)
                               ENDDO
                               DO ilop = 1,atoms%nlo(n)
                                  IF (atoms%llo(ilop,n).EQ.l) THEN
                                     DO m = -l,l
                                        hsomtx(i1,j1,i,j) = hsomtx(i1,j1,i,j) + &
                                             rsoc%rsoploplop(n,ilop,ilo,i1,j1) * &
                                             chelp(m,i,ilop,na,jsp) * c_c(m,ilo,na)
                                     ENDDO
                                  ENDIF
                               ENDDO
                            ENDIF
                         ENDDO
                         ! end lo's
                      ENDDO

                   ENDDO
                !
             ENDDO
             !!i
          ENDDO
          !!j
          !$OMP END DO
          DEALLOCATE (c_a,c_b,c_c)
          !$OMP END PARALLEL
       ENDDO
       !!jsp1
    ENDDO
    !!jsp
    !
    !---> update hamiltonian matrices: lower triangle
    !
    DO i = 1,nsz(1)
       DO j = 1,nsz(input%jspins)
          hsomtx(2,1,j,i) = CONJG(hsomtx(1,2,i,j))
       ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE hsoham
END MODULE m_hsoham
