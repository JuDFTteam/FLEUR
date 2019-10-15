MODULE m_umtx
   USE m_juDFT
   !*********************************************************************
   !* The calculation of the "U"-contribution to Hartree-Fock matrix.   *
   !*-------------------------------------------------------------------*
   !* Extension to multiple U per atom type by G.M. 2017                *
   !*********************************************************************
CONTAINS
   SUBROUTINE umtx(u_in,n_u,f0,f2,f4,f6,u)

      USE m_constants
      USE m_sgaunt
      USE m_types
      IMPLICIT NONE

      INTEGER,       INTENT(IN)  :: n_u
      TYPE(t_utype), INTENT(IN)  :: u_in(:)
      REAL,          INTENT(IN)  :: f0(:),f2(:),f4(:),f6(:)
      REAL,          INTENT(OUT) :: u(-lmaxU_const:,-lmaxU_const:,-lmaxU_const:,-lmaxU_const:,:)

      INTEGER, PARAMETER         :: lmaxw=3,lmmaxw1=(2*lmaxw+2)**2

      INTEGER i,j,k,l,m,mk,nfk,itype,i_u
      INTEGER m1,m2,m3,m4,lm1,lm2,lm3,lm4,kf
      REAL    uk,uq,avu,avj,cgk1,cgk2,tol
      REAL    fk(lmaxU_const+1,n_u)
      REAL,   ALLOCATABLE :: c(:,:,:)
      !
      tol = 1.0e-14
      !
      ! transformation to Hr-units:
      !
      DO i_u = 1, n_u
         itype = u_in(i_u)%atomType
         l = u_in(i_u)%l
         fk(1,i_u) = f0(i_u) / hartree_to_ev_const
         fk(2,i_u) = f2(i_u) / hartree_to_ev_const
         fk(3,i_u) = f4(i_u) / hartree_to_ev_const
         IF (l.EQ.3) THEN
            fk(4,i_u) = f6(i_u) / hartree_to_ev_const
         ELSE IF (l.GT.3) THEN
            CALL juDFT_error("LDA+U for p, d or f-states!", calledby="umtx")
         END IF
      END DO
      !
      ! evaluate Gaunt parameter
      !
      ALLOCATE(c(0:2*lmaxw+1,lmmaxw1,lmmaxw1))
      DO k = 1,lmmaxw1
         DO j = 1,lmmaxw1
            DO i = 0,2*lmaxw+1
               c(i,j,k) = 0.0
            END DO
         END DO
      END DO

      CALL sgaunt(lmaxw,lmmaxw1,lmaxU_const,c)

      DO i_u = 1, n_u
         l = u_in(i_u)%l
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
                        END DO                   ! mk
                        IF (ABS(uq).LT.tol) CYCLE
                        nfk=k/2+1
                        uk=uk+uq*fk(nfk,i_u)*4*pi_const/(2*k+1)
                     END DO                     ! k
                     u(m1,m2,m3,m4,i_u)=uk
                  END DO                       ! m4 etc.
               END DO
            END DO
         END DO
         avu=0.e0
         avj=0.e0

         DO i = -l,l
            DO j = -l,l
               avu = avu+u(i,j,i,j,i_u)
               avj = avj+(u(i,j,i,j,i_u)-u(i,j,j,i,i_u))
            END DO
         END DO
         avu = avu/(2*l+1)/(2*l+1)
         avj = avj/(2*l+1)/(2*l)
         avj = avu-avJ
         !        WRITE (6,*) 'U-matr:'
         !        IF (l.eq.2) WRITE (6,111) ((u(i,j,i,j,i_u),i=-l,l),j=-l,l)
         !        IF (l.eq.3) WRITE (6,211) ((u(i,j,i,j,i_u),i=-l,l),j=-l,l)
         !        WRITE (6,*) 'J-matr:'
         !        IF (l.eq.2) WRITE (6,111) ((u(i,j,j,i,i_u),i=-l,l),j=-l,l)
         !        IF (l.eq.3) WRITE (6,211) ((u(i,j,j,i,i_u),i=-l,l),j=-l,l)
         !         PRINT*,'U-av:',avu*hartree_to_ev_const
         !         PRINT*,'J-av:',avj*hartree_to_ev_const
111      FORMAT (5f8.4)
211      FORMAT (7f8.4)
112      FORMAT (10e20.10)
         !c         WRITE (9,112) ((((u(i,j,k,m,i_u),i=-l,l),j=-l,l),k=-l,l),m=-l,l)
      END DO
      DEALLOCATE (c)

   END SUBROUTINE umtx
END MODULE m_umtx
