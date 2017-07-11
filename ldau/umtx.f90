MODULE m_umtx
  USE m_juDFT
  !*********************************************************************
  !* The calculation of the "U"-contribution to Hartree-Fock matrix.   *
  !*********************************************************************
CONTAINS
  SUBROUTINE umtx(&
       lmaxb,ntype,n_u,lda_u,f0,f2,f4,f6,&
       u)

    USE m_constants
    USE m_sgaunt
    IMPLICIT NONE

    INTEGER, PARAMETER   :: lmaxw=3,lmmaxw1=(2*lmaxw+2)**2
    INTEGER, INTENT (IN) :: n_u,lmaxb,ntype
    INTEGER, INTENT (IN) :: lda_u(ntype)
    REAL,    INTENT (IN) :: f0(n_u),f2(n_u),f4(n_u),f6(n_u)
    REAL,    INTENT (OUT) :: u(-lmaxb:lmaxb,-lmaxb:lmaxb,&
         -lmaxb:lmaxb,-lmaxb:lmaxb,n_u)

    INTEGER i,j,k,l,m,mk,nfk,n,itype
    INTEGER m1,m2,m3,m4,lm1,lm2,lm3,lm4,kf,l_l(n_u)
    REAL    uk,uq,avu,avj,cgk1,cgk2,tol
    REAL    fk(lmaxb+1,n_u)
    REAL,   ALLOCATABLE :: c(:,:,:)
    !
    tol = 1.0e-14
    !
    ! transformation to Hr-units:
    !
    n = 0
    DO itype = 1,ntype
       IF (lda_u(itype).GE.0) THEN
          n = n + 1
          l_l(n) = lda_u(itype)
          fk(1,n) = f0(n) / hartree_to_ev_const
          fk(2,n) = f2(n) / hartree_to_ev_const
          fk(3,n) = f4(n) / hartree_to_ev_const
          IF (l_l(n).EQ.3) THEN
             fk(4,n) = f6(n) / hartree_to_ev_const
          ELSEIF (l_l(n).GT.3) THEN
             CALL juDFT_error("LDA+U for p, d or f-states!", calledby="umtx")
          ENDIF
       ENDIF
    ENDDO
    !
    ! evaluate Gaunt parameter
    !
    ALLOCATE( c(0:2*lmaxw+1,lmmaxw1,lmmaxw1) )
    DO k = 1,lmmaxw1
       DO j = 1,lmmaxw1
          DO i = 0,2*lmaxw+1
             c(i,j,k) = 0.0
          ENDDO
       ENDDO
    ENDDO
    CALL sgaunt(lmaxw,lmmaxw1,lmaxb,&
         c)
    !
    ! lda_u(n) is here only the 'l' for atom 'n'
    !
    DO  n = 1,n_u                      !!! over d-atoms
       l = l_l(n)
       kf = 2*l
       DO m1 = -l,l
          lm1 = l*(l+1)+m1+1
          DO m2 = -l,l
             lm2 = l*(l+1)+m2+1
             DO m3 = -l,l
                lm3 = l*(l+1)+m3+1
                DO m4 = -l,l
                   lm4 = l*(l+1)+m4+1
                   uk = 0.0e0
                   DO k=0,kf,2
                      uq = 0.e0
                      DO mk=-k,k
                         IF (mk.NE.m1-m3)  CYCLE
                         cgk1 = c(k/2,lm1,lm3)
                         IF (ABS(cgk1).LT.tol) CYCLE
                         IF (mk.NE.m4-m2)  CYCLE
                         cgk2 = c(k/2,lm4,lm2)
                         IF (ABS(cgk2).LT.tol) CYCLE
                         uq = uq+cgk1*cgk2
                      ENDDO                   ! mk
                      IF (ABS(uq).LT.tol) CYCLE
                      nfk=k/2+1
                      uk=uk+uq*fk(nfk,n)*4*pi_const/(2*k+1)
                   ENDDO                     ! k
                   u(m1,m2,m3,m4,n)=uk
                ENDDO                       ! m4 etc.
             ENDDO
          ENDDO
       ENDDO
       avu=0.e0
       avj=0.e0

       DO i = -l,l
          DO j = -l,l
             avu = avu+u(i,j,i,j,n)
             avj = avj+(u(i,j,i,j,n)-u(i,j,j,i,n))
          ENDDO
       ENDDO
       avu = avu/(2*l+1)/(2*l+1)
       avj = avj/(2*l+1)/(2*l)
       avj = avu-avJ
       !        WRITE (6,*) 'U-matr:'
       !        IF (l.eq.2) WRITE (6,111) ((u(i,j,i,j,n),i=-l,l),j=-l,l)
       !        IF (l.eq.3) WRITE (6,211) ((u(i,j,i,j,n),i=-l,l),j=-l,l)
       !        WRITE (6,*) 'J-matr:'
       !        IF (l.eq.2) WRITE (6,111) ((u(i,j,j,i,n),i=-l,l),j=-l,l)
       !        IF (l.eq.3) WRITE (6,211) ((u(i,j,j,i,n),i=-l,l),j=-l,l)
       !         PRINT*,'U-av:',avu*hartree_to_ev_const
       !         PRINT*,'J-av:',avj*hartree_to_ev_const
111    FORMAT (5f8.4)
211    FORMAT (7f8.4)
112    FORMAT (10e20.10)
       !c         WRITE (9,112) ((((u(i,j,k,m,n),i=-l,l),j=-l,l),k=-l,l),m=-l,l)
    ENDDO
    DEALLOCATE (c)

  END SUBROUTINE umtx
END MODULE m_umtx
