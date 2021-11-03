MODULE m_uj2f
   USE m_juDFT
   !  *********************************************************************
   !  * The calculation of slater integrals from u&j                      *
   !  * input in eV; output in htr.                                       *
   !  *-------------------------------------------------------------------*
   !  * Extension to multiple U per atom type by G.M. 2017                *
   !  * Extension for uses beyond LDA+U by H.J 2019                       *
   !  *********************************************************************
   USE m_types

   IMPLICIT NONE

   INTERFACE uj2f
      procedure :: uj2f_simple, uj2f_spins, uj2f_single
      procedure :: uj2f_single_onelist, uj2f_multiple_onelist
   END INTERFACE

   CONTAINS

   subroutine uj2f_single_onelist(jspins,u_in,f)

      INTEGER,          INTENT(IN)  :: jspins
      TYPE(t_utype),    INTENT(IN)  :: u_in
      REAL,             INTENT(OUT) :: f(0:6)

      real :: f0, f2, f4, f6

      f = 0.0
      CALL uj2f_single(jspins,u_in,f0,f2,f4,f6)
      f(0) = f0
      f(2) = f2
      f(4) = f4
      f(6) = f6

   end subroutine

   subroutine uj2f_multiple_onelist(jspins,u_in,n_u,f)

      INTEGER,          INTENT(IN)  :: jspins
      INTEGER,          INTENT(IN)  :: n_u
      TYPE(t_utype),    INTENT(IN)  :: u_in(:)
      REAL, ALLOCATABLE,INTENT(OUT) :: f(:,:)

      real :: f0(n_u), f2(n_u), f4(n_u), f6(n_u)

      allocate(f(0:6,n_u), source=0.0)

      CALL uj2f_simple(jspins,u_in,n_u,f0,f2,f4,f6)

      f(0,:) = f0
      f(2,:) = f2
      f(4,:) = f4
      f(6,:) = f6

   end subroutine

   SUBROUTINE uj2f_single(jspins,u_in,f0,f2,f4,f6)

      INTEGER,          INTENT(IN)  :: jspins
      TYPE(t_utype),    INTENT(IN)  :: u_in
      REAL,             INTENT(OUT) :: f0,f2
      REAL,             INTENT(OUT) :: f4,f6

      REAL :: f0List(1),f2List(1)
      REAL :: f4List(1),f6List(1)

      CALL uj2f_simple(jspins,[u_in],1,f0List,f2List,f4List,f6List)

      f0 = f0List(1)
      f2 = f2List(1)
      f4 = f4List(1)
      f6 = f6List(1)

   END SUBROUTINE uj2f_single

   SUBROUTINE uj2f_simple(jspins,u_in,n_u,f0,f2,f4,f6)

      INTEGER,          INTENT(IN)  :: jspins
      INTEGER,          INTENT(IN)  :: n_u
      TYPE(t_utype),    INTENT(IN)  :: u_in(:)
      REAL,             INTENT(OUT) :: f0(:),f2(:)
      REAL,             INTENT(OUT) :: f4(:),f6(:)

      REAL :: f0Spins(n_u,jspins),f2Spins(n_u,jspins)
      REAL :: f4Spins(n_u,jspins),f6Spins(n_u,jspins)

      CALL uj2f_spins(jspins,u_in,n_u,f0Spins,f2Spins,f4Spins,f6Spins)

      f0 = (f0Spins(:,1) + f0Spins(:,jspins))/ 2.0
      f2 = (f2Spins(:,1) + f2Spins(:,jspins))/ 2.0
      f4 = (f4Spins(:,1) + f4Spins(:,jspins))/ 2.0
      f6 = (f6Spins(:,1) + f6Spins(:,jspins))/ 2.0

   END SUBROUTINE uj2f_simple

   SUBROUTINE uj2f_spins(jspins,u_in,n_u,f0,f2,f4,f6)

      INTEGER,          INTENT(IN)  :: jspins
      INTEGER,          INTENT(IN)  :: n_u
      TYPE(t_utype),    INTENT(IN)  :: u_in(:)
      REAL,             INTENT(OUT) :: f0(:,:),f2(:,:)
      REAL,             INTENT(OUT) :: f4(:,:),f6(:,:)

      INTEGER l,itype,ltest,ispin,i_u
      REAL u,j,a,ftest(4)
      LOGICAL l_exist

      l_exist=.FALSE.
      INQUIRE (file='slaterf',exist=l_exist)

      IF (l_exist) THEN
         !
         ! --> f's have been calculated in cored ; read from file
         !
         OPEN (45,file='slaterf',form='formatted',status='old')
         DO ispin = 1, jspins
            DO i_u = 1, n_u
               itype = u_in(i_u)%atomType
               l = u_in(i_u)%l
               f2(i_u,ispin)=0.0 ; f4(i_u,ispin)=0.0 ; f6(i_u,ispin)=0.0
100            READ (45,'(i3,4f20.10)') ltest,ftest(1:4)
               IF (ltest.EQ.l) THEN
                  f0(i_u,ispin) = ftest(1)
                  IF (l.GT.0) THEN
                     f2(i_u,ispin) = ftest(2)
                     IF (l.GT.1) THEN
                        f4(i_u,ispin) = ftest(3)
                        IF (l.GT.2) THEN
                           f6(i_u,ispin) = ftest(4)
                        END IF
                     END IF
                  END IF
               ELSE
                  GOTO 100
               END IF
               READ (45,'(i3,4f20.10)') ltest,ftest(1)
               !                IF (ltest.EQ.0) THEN
               !                   f0(n,ispin) = f0(n,ispin) - ftest(1)
               !                ENDIF

               !              write(*,*) n,ispin,l,f0(n,ispin),f2(n,ispin),
               !    +                              f4(n,ispin),f6(n,ispin)
            END DO ! n_u
         ENDDO
         CLOSE (45)
      ELSE
         !
         ! lda_u%l: orb.mom; lda_u%u,j: in eV
         !
         DO i_u = 1, n_u

            itype = u_in(i_u)%atomType
            l = u_in(i_u)%l
            u = u_in(i_u)%u
            j = u_in(i_u)%j
            !
            !        l.eq.0 :  f0 = u (the l=0 and l=1 case approximated g.b.`01)
            !
            IF (l.EQ.0) THEN
               f0(i_u,1) = u
               f2(i_u,1) = 0.0
               f4(i_u,1) = 0.0
               f6(i_u,1) = 0.0
               IF (j>0.00001) CALL juDFT_error("lda+u: no magnetic s-states", calledby ="uj2f")
               !
               !        l == 1 :  j = f2 / 5  (from PRL 80,5758 g.b.)
               !
            ELSE IF (l.EQ.1) THEN
               f0(i_u,1) = u
               f2(i_u,1) = 5.0*j
               f4(i_u,1) = 0.0
               f6(i_u,1) = 0.0
               !
               !        l.eq.2 : 3d: j=(f2+f4)/14; f4/f2 = 0.625
               !
            ELSE IF (l.EQ.2) THEN
               !             PRINT*, 'd-states'
               f0(i_u,1) = u
               f2(i_u,1) = 14.0*j/1.625
               f4(i_u,1) = f2(i_u,1)*0.625
               f6(i_u,1) = 0.0
               !
               !        l.eq. 3 : 4f: j=(286f2+195f4+250f6)/6435; f2/f4 = 675/451; f2/f6=2025/1001
               !
            ELSE IF (l.EQ.3) THEN
               !             PRINT*, 'f-states'
               f0(i_u,1) = u
               a= 286.0 + 195.0*451.0/675.0 + 250.0*1001.0/2025.0
               f2(i_u,1) = 6435.0*j/a
               f4(i_u,1) = 451.0/675.0*f2(i_u,1)
               f6(i_u,1) = 1001.0/2025.0*f2(i_u,1)
            ELSE
               CALL juDFT_error('lda+U is restricted to l<=3 !', calledby="uj2f")
            END IF
            IF (jspins.EQ.2) THEN
               f0(i_u,jspins) = f0(i_u,1)
               f2(i_u,jspins) = f2(i_u,1)
               f4(i_u,jspins) = f4(i_u,1)
               f6(i_u,jspins) = f6(i_u,1)
            ENDIF

         END DO ! n_u
      ENDIF

   END SUBROUTINE uj2f_spins
END MODULE m_uj2f
