      MODULE m_inpeigdim
!*********************************************************************
!     inputs the necessary quantities for the eigenvalue part (energy
!     parameters, k-points, wavefunction cutoffs, etc.).
!                  m. weinert   jan. 1987
!*********************************************************************
      CONTAINS
      SUBROUTINE inpeig_dim(input,obsolete,cell,noco, oneD,jij,kpts,dimension,stars)

      USE m_constants, ONLY : pi_const,tpi_const
      USE m_types
      use m_apwsdim
      IMPLICIT NONE
      TYPE(t_input),INTENT(INOUT)     :: input
      TYPE(t_obsolete),INTENT(INOUT)  :: obsolete
      TYPE(t_cell),INTENT(INOUT)      :: cell
      TYPE(t_noco),INTENT(INOUT)      :: noco
      TYPE(t_stars),INTENT(INOUT)     :: stars
      TYPE(t_dimension),INTENT(INOUT) :: dimension
      TYPE(t_kpts),INTENT(INOUT)      :: kpts
      TYPE(t_oneD),INTENT(INOUT)      :: oneD
      TYPE(t_Jij),INTENT(INOUT)       :: Jij
!-odim
!-odim
      INTEGER nk,nq,i,nv,nv2,j,kq1,kq2,kq3
      REAL  s1,s2,scale,bk(3)
      LOGICAL xyu
   !     ..
      kpts%nkpt = 0 ; dimension%nvd = 0 ; dimension%nv2d = 0
      stars%kq1d = 0 ; stars%kq2d = 0 ; stars%kq3d = 0
      !cell%aamat=matmul(transpose(cell%amat),cell%amat)
      cell%bbmat=matmul(cell%bmat,transpose(cell%bmat))
!
      jij%nqpt=1
      IF (jij%l_J) THEN
        OPEN (113,file='qpts',form='formatted',status='old')
        READ (113,*) jij%nqpt
      ENDIF
      OPEN (41,file='kpts',form='formatted',status='old')

!--->    k-mesh: given in units of the reciprocal lattice basis vectors
!--->    scale is a factor to make input easier (default=1.0). k-pt
!--->    weights can be relative weights since they are renormalized.
!--->    input: for bulk - k1,k2,k3,wtkpt
!--->           for film - k1,k2,wtkpt
!--->    read k-points from file 41='kpts'
         IF (input%film) THEN
            READ (41,fmt=8050) kpts%nkpt,scale,xyu
         ELSE
            READ (41,*) kpts%nkpt,scale
         END IF
 8030    FORMAT (4f10.5)
 8040    FORMAT (i5,f20.10)
 8050    FORMAT (i5,f20.10,3x,l1)

         kpts%nkpt = max(kpts%nkpt,kpts%nkpt)
 8060    FORMAT (i5,f20.10)
         IF (scale.EQ.0.0) scale = 1.0
         DO nq=1,jij%nqpt
           IF(jij%l_J) THEN
             READ (113,fmt=8070) noco%qss(1),noco%qss(2),noco%qss(3)
 8070        FORMAT(2(f14.10,1x),f14.10)
           ENDIF

           DO nk = 1,kpts%nkpt
             IF(input%film.and..not.oneD%odd%d1) THEN
                READ (41,fmt=8080) (bk(i),i=1,2)
 8080        FORMAT (3f10.5)
             ELSE
                READ (41,fmt=8030) (bk(i),i=1,3)
             ENDIF
             IF (oneD%odd%d1) THEN
               bk(1) = 0.
               bk(2) = 0.
             ELSEIF (input%film .AND. .NOT.oneD%odd%d1) THEN
               bk(3) = 0.0
               IF (xyu) THEN
!            transform to cartesian coordinates
                  IF (cell%latnam.EQ.'hex') THEN
                     bk(1) = bk(1)*tpi_const/cell%amat(2,2)
                     bk(2) = bk(2)*pi_const/cell%amat(1,1)
                  ELSE
                     bk(1) = bk(1)*pi_const/cell%amat(1,1)
                     bk(2) = bk(2)*pi_const/cell%amat(2,2)
                  END IF
!            transform to internal coordinates
                  s1 = 0.0
                  s2 = 0.0
                  DO j = 1,2
                     s1 = s1 + cell%amat(j,1)*bk(j)/tpi_const
                     s2 = s2 + cell%amat(j,2)*bk(j)/tpi_const
                  ENDDO
                  bk(1) = s1
                  bk(2) = s2
               END IF
             END IF
             DO i = 1,3
               bk(i) = bk(i)/scale
             ENDDO
             CALL apws_dim(&
     &                     bk(:),cell,input,noco,oneD,&
     &                     nv,nv2,kq1,kq2,kq3)
             stars%kq1d = max(kq1,stars%kq1d)
             stars%kq2d = max(kq2,stars%kq2d)
             stars%kq3d = max(kq3,stars%kq3d)
             
             dimension%nvd = max(dimension%nvd,nv)
             dimension%nv2d = max(dimension%nv2d,nv2)

           ENDDO ! k=pts
           REWIND(41)
           READ (41,*)
        ENDDO   ! q-pts

      IF (jij%l_J) THEN
       CLOSE(113)
      ENDIF
      CLOSE (41)

      END SUBROUTINE inpeig_dim
      END MODULE m_inpeigdim
