!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_od_strgn1
      CONTAINS
        SUBROUTINE od_strgn1(&
             &     xcpot,cell,sym,oneD)

          !**********************************************************
          !     one-dimensional stars generator
          !                                 Y.Mokrousov, 2002-2003
          !***********************************************************
          USE m_types

          IMPLICIT NONE
          CLASS(t_xcpot),INTENT(IN) :: xcpot
          TYPE(t_cell),INTENT(IN)   :: cell
          TYPE(t_sym),INTENT(IN)    :: sym
          TYPE(t_oneD),INTENT(INOUT):: oneD



          INTEGER nfftx_1,nffty_1,kfx_1,kfy_1,kfft_1,m
          REAL gfx_1,gfy_1
          INTEGER i ,m1,z,z1

          !     odl%nn2d = odi%nn2d
          !     odg%nn2d = odi%nn2d

          !     ALLOCATE ( odi%kv(1:2,odi%n2d),
          !    &           odi%ig(-odi%k3:odi%k3,-odi%M:odi%M),
          !    &           odi%nst2(1:odi%n2d),
          !    &           odl%igf(0:odl%nn2d-1,1:2),
          !    &           odl%pgf(0:odl%nn2d-1) )
          oneD%igfft1(0:oneD%odd%nn2d-1,1:2) = 0
          oneD%pgfft1(0:oneD%odd%nn2d-1) = 0.

          IF (xcpot%needs_grad()) THEN
             !         ALLOCATE ( odg%pgfx(0:odg%nn2d-1),
             !     &              odg%pgfy(0:odg%nn2d-1),
             !     &              odg%pgfxx(0:odg%nn2d-1),
             !     &              odg%pgfxy(0:odg%nn2d-1),
             !     &              odg%pgfyy(0:odg%nn2d-1) )
             oneD%pgft1x(0:oneD%odd%nn2d-1) = 0.
             oneD%pgft1y(0:oneD%odd%nn2d-1) = 0.
             oneD%pgft1xx(0:oneD%odd%nn2d-1) = 0.
             oneD%pgft1xy(0:oneD%odd%nn2d-1) = 0.
             oneD%pgft1yy(0:oneD%odd%nn2d-1) = 0.
          END IF

          !      odi%nq2 = odi%n2d
          !      odd%kimax2 = odd%n2d - 1
          oneD%nstr1(1:oneD%odd%n2d) = 0

          !---> generating mapping arrays       

          i = 0

          DO  z1 = 0,2*oneD%odd%k3
             DO  m1 = 0,2*oneD%odd%M
                IF (m1.LE.oneD%odd%M) m = m1
                IF (m1.GT.oneD%odd%M) m = oneD%odd%M - m1
                IF (z1.LE.oneD%odd%k3) z = z1
                IF (z1.GT.oneD%odd%k3) z = oneD%odd%k3-z1
                IF (oneD%odd%chi.EQ.1) THEN
                   IF (sym%zrfs) THEN
                      IF (MOD(m,oneD%odd%rot).EQ.0) THEN
                         IF (z.GE.0) THEN
                            i = i+1
                            oneD%ig1(z,m) = i
                            oneD%kv1(1,i) = z
                            oneD%kv1(2,i) = m
                            IF (z.EQ.0) THEN
                               oneD%nstr1(i) = 1
                            ELSE
                               oneD%nstr1(i) = 2
                            END IF
                         ELSE
                            oneD%ig1(z,m) = oneD%ig1(-z,m)
                         END IF
                      ELSE
                         oneD%ig1(z,m) = 0
                      END IF
                   ELSE IF (.NOT.sym%zrfs) THEN
                      IF (MOD(m,oneD%odd%rot).EQ.0) THEN
                         i = i+1
                         oneD%ig1(z,m) = i
                         oneD%kv1(1,i) = z
                         oneD%kv1(2,i) = m
                         oneD%nstr1(i) = 1
                      ELSE
                         oneD%ig1(z,m) = 0
                      END IF
                   END IF
                ELSE
                   IF (MOD(m+(oneD%odd%rot)*z,oneD%odd%chi).EQ.0) THEN
                      i = i+1
                      oneD%ig1(z,m) = i
                      oneD%kv1(1,i) = z
                      oneD%kv1(2,i) = m
                      oneD%nstr1(i) = 1
                   ELSE
                      oneD%ig1(z,m) = 0 
                   END IF
                END IF
             ENDDO
          ENDDO

          !---> preparations for 2dim vacuum fft
          !---> at the moment we have no symmetries

          nfftx_1 = 3*oneD%odd%k3
          nffty_1 = 3*oneD%odd%M
          DO  i = 1,oneD%odd%nq2
             kfx_1 = oneD%kv1(1,i)
             kfy_1 = oneD%kv1(2,i)
             IF (kfx_1.LT.0) kfx_1 = kfx_1 + nfftx_1
             IF (kfy_1.LT.0) kfy_1 = kfy_1 + nffty_1
             kfft_1 = kfx_1 + kfy_1*nfftx_1
             oneD%igfft1(i-1,1) = i
             oneD%igfft1(i-1,2) = kfft_1
             oneD%pgfft1(i-1) = 1.
          ENDDO

          IF (xcpot%needs_grad()) THEN
             DO  i = 1,oneD%odd%nq2
                kfx_1 = oneD%kv1(1,i)
                kfy_1 = oneD%kv1(2,i)
                gfx_1 = cell%bmat(3,3)*kfx_1
                gfy_1 = kfy_1
                oneD%pgft1x(i-1)  = gfx_1
                oneD%pgft1y(i-1)  = gfy_1
                oneD%pgft1xx(i-1) = gfx_1*gfx_1
                oneD%pgft1yy(i-1) = gfy_1*gfy_1
                oneD%pgft1xy(i-1) = gfx_1*gfy_1
             ENDDO
          END IF

          !     odi%kimax2 = odi%nq2 - 1

          RETURN
        END SUBROUTINE od_strgn1
      END MODULE m_od_strgn1
 

