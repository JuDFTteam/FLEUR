      MODULE m_dosbin
c*********************************************************
c    this subroutine generate the idos, the ldos, the partial
c    ldos in the spheres and the z-dependent dos for the
c    vacuum.
c
c*********************************************************
      CONTAINS
      SUBROUTINE dos_bin(
     >                   jspins,ndos,ne,emin,emax,neigd,
     >                   nkpt,neig,wtkpt,eig,qal,
     <                   g)
c
      IMPLICIT NONE
c
      INTEGER,INTENT (IN) :: jspins,ndos,ne,nkpt,neigd,neig(nkpt)
      REAL,   INTENT (IN) :: emin,emax,wtkpt(nkpt)
      REAL,   INTENT (IN) :: eig(neigd,nkpt),qal(ndos,neigd,nkpt)
      REAL,   INTENT (OUT):: g(ne,ndos)

c------> local variables
c
      INTEGER  nl,k,j, i
      REAL  de,wk
c     ..
      de=(emax-emin)/real(ne-1)
      g=0.0
c
c----> put weights in the right bins
c
      DO k = 1 , nkpt
         wk = 2.*wtkpt(k)/real(jspins)/de
         DO j = 1 , neig(k)
            i = NINT((eig(j,k)-emin)/de) + 1
            IF ( (i.LE.ne) .AND. (i.GE.1) ) THEN
                  DO nl = 1 , ndos
                     g(i,nl) = g(i,nl) + wk*qal(nl,j,k)
                  ENDDO
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE dos_bin
      END MODULE m_dosbin
