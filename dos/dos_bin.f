      MODULE m_dosbin
c*********************************************************
c    this subroutine generate the idos, the ldos, the partial
c    ldos in the spheres and the z-dependent dos for the
c    vacuum.
c
c*********************************************************
      CONTAINS
      SUBROUTINE dos_bin(e,neig,wtkpt,eig,qal,
     <                   g)
c
      IMPLICIT NONE

c
      INTEGER,INTENT (IN) :: neig(:,:)
      REAL,   INTENT (IN) :: wtkpt(:),e(:)
      REAL,   INTENT (IN) :: eig(:,:,:),qal(:,:,:)
      REAL,   INTENT (OUT):: g(:,:)

c------> local variables
c
      INTEGER  nl,k,j, i,js
      REAL  de,wk,emin
c     ..
      de=e(2)-e(1)
      g=0.0
      emin=minval(e)
c
c----> put weights in the right bins
c
      DO js=1,size(qal,3)
        DO k = 1 , size(qal,2)
          wk = wtkpt(k)/de
          DO j = 1 , neig(k,js)
            i = NINT((eig(j,k,js)-emin)/de) + 1
            IF ( (i.LE.size(e)) .AND. (i.GE.1) ) THEN
              g(i,js) = g(i,js) + wk*qal(j,k,js)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE dos_bin
      END MODULE m_dosbin
