      MODULE m_xyavden
      CONTAINS
      SUBROUTINE xy_av_den(&
     &                    stars,vacuum,cell,psq,rht)

      USE m_types
      USE m_cfft
      IMPLICIT NONE
      TYPE(t_vacuum),INTENT(IN)   :: vacuum
      TYPE(t_stars),INTENT(IN)   :: stars
      TYPE(t_cell),INTENT(IN)   :: cell
      REAL,    INTENT(IN) :: rht(vacuum%nmzd,2)
      COMPLEX, INTENT(IN) :: psq(stars%ng3)

      INTEGER  ivfft,i,j,k
      REAL     ani,z
      REAL,    ALLOCATABLE :: af1(:),bf1(:)

      ivfft =  3*stars%mx3
      ALLOCATE (af1(ivfft),bf1(ivfft))

      af1(:) = 0.0 ; bf1(:) = 0.0
      DO i = 1, stars%ng3
        IF (stars%ig2(i) == 1) THEN
          k = stars%kv3(3,i)
          IF ( k < 0 ) THEN
            k = ivfft + k + 1 
          ELSE
            k = k + 1
          ENDIF
          af1(k) = real(psq(i))
          bf1(k) = aimag(psq(i))
        ENDIF
      ENDDO

      CALL cfft(af1,bf1,ivfft,ivfft,ivfft,+1)

      OPEN(77,file='qws',status='unknown')
      j = 1
      k = 3 - 2*j
      DO i = vacuum%nmz,1,-1
        z = (vacuum%dvac/2 + (i-1)*vacuum%delz) * k 
        WRITE(77,'(2f20.10)') z,rht(i,j)*cell%area
      ENDDO
      ani = 1.0/real(ivfft)
      j = 0
      DO i = 0,ivfft - 1
        j = j + 1
        z = cell%amat(3,3)*i*ani
        IF (z > cell%amat(3,3)/2) z = z - cell%amat(3,3)
        IF ( abs(z) < vacuum%dvac/2 ) THEN
          WRITE(77,'(2f20.10)') z,af1(j)*cell%area
        ENDIF
      ENDDO
      j = 1
      k = 3 - 2*j
      DO i = 1, vacuum%nmz
        z = (vacuum%dvac/2 + (i-1)*vacuum%delz) * k 
        WRITE(77,'(2f20.10)') z,rht(i,j)*cell%area
      ENDDO
      
      CLOSE(77)
      DEALLOCATE (af1,bf1)
      STOP

      END SUBROUTINE xy_av_den
      END MODULE m_xyavden

