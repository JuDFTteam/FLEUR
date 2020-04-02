MODULE m_fertetra

   USE m_types
   USE m_juDFT
   USE m_tetrahedronInit

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE fertetra(input,noco,kpts,mpi,ne,eig,ef,w,seigv)

      TYPE(t_kpts),        INTENT(IN)     :: kpts
      TYPE(t_noco),        INTENT(IN)     :: noco
      TYPE(t_input),       INTENT(IN)     :: input
      TYPE(t_mpi),         INTENT(IN)     :: mpi
      INTEGER,             INTENT(IN)     :: ne(:,:)
      REAL,                INTENT(IN)     :: eig(:,:,:)
      REAL,                INTENT(INOUT)  :: seigv
      REAL,                INTENT(INOUT)  :: ef
      REAL,                INTENT(INOUT)  :: w(:,:,:)

      INTEGER :: jspin,jspins,ikpt,it,iBand
      REAL    :: dlow,dup,dfermi,s1,s,chmom
      REAL    :: lowBound,upperBound


      jspins = MERGE(1,input%jspins,noco%l_noco)
      !---------------------------------------------
      !Find the interval, where ef should be located
      !---------------------------------------------
      !Initial guess
      lowBound = MINVAL(eig)-0.01
      upperBound = ef + 0.2

      !First check the lower bound
      dlow = 0.0
      DO jspin = 1, jspins
         CALL tetrahedronInit(kpts,eig(:,:,jspin),MINVAL(ne(:,jspin)),&
                              lowBound,input%film,w(:,:,jspin))
         DO ikpt = 1, kpts%nkpt
            DO iBand = 1, MINVAL(ne(:,jspin))
               dlow = dlow + w(iBand,ikpt,jspin) * 2.0/input%jspins
            ENDDO
         ENDDO
      ENDDO

      IF (dlow.GT.input%zelec) THEN
         WRITE(6,9000) lowBound,dlow,input%zelec
         CALL juDFT_error("valence band too high ",calledby="fertetra")
      ENDIF
9000  FORMAT (' valence band too high ',/,&
              '  elow ',f10.5,' dlow ',f10.5,' nelec ',f10.5)

      it = 0
      DO
         !Now check the upper bound
         dup = 0.0
         DO jspin = 1, jspins
            CALL tetrahedronInit(kpts,eig(:,:,jspin),MINVAL(ne(:,jspin)),&
                                 upperBound,input%film,w(:,:,jspin))
            DO ikpt = 1, kpts%nkpt
               DO iBand = 1, MINVAL(ne(:,jspin))
                  dup = dup + w(iBand,ikpt,jspin) * 2.0/input%jspins
               ENDDO
            ENDDO
         ENDDO

         IF (dup.GT.input%zelec) THEN
            EXIT
         ELSE
            !Raise the upper bound
            upperBound = upperBound + 0.2
            it = it + 1
            IF(it.GT.10) THEN
               WRITE (6,9100) upperBound,dup,input%zelec
9100           FORMAT (' valence band too low ',/,&
                       '  eup  ',f10.5,' dup  ',f10.5,' nelec ',f10.5)
               CALL juDFT_error("valence band too low ",calledby ="fertetra")
            ENDIF
         ENDIF
      ENDDO

      !-----------------------------------------------------------------------------------
      !Now that the fermi energy is guaranteed to be in the interval [lowBound,upperBound]
      !We use the bisection method to find it
      !-----------------------------------------------------------------------------------

      DO WHILE(upperBound-lowBound.GT.1e-10)

         ef = (lowBound+upperBound)/2.0
         dfermi = 0.0
         DO jspin = 1, jspins
            !-------------------------------------------------------
            ! Compute the weights for charge density integration
            !-------------------------------------------------------
            CALL tetrahedronInit(kpts,eig(:,:,jspin),MINVAL(ne(:,jspin)),&
                                 ef,input%film,w(:,:,jspin))
            DO ikpt = 1, kpts%nkpt
               DO iBand = 1, MINVAL(ne(:,jspin))
                  dfermi = dfermi + w(iBand,ikpt,jspin) * 2.0/input%jspins
               ENDDO
            ENDDO
         ENDDO
         IF(dfermi.GT.input%zelec) THEN
            !Occupation to large -> search in the right interval
            upperBound = ef
         ELSE IF(dfermi.LE.input%zelec) THEN
            !Occupation to small -> search in the right interval
            lowBound = ef
         ENDIF
      ENDDO
      WRITE (6,9200) ef,dfermi,input%zelec
9200  FORMAT (//,'Tetrahedron method: ',//,'   fermi energy ',f10.5,&
                 ' dtot ',f10.5,' nelec ',f10.5)

      !----------------------------------------------
      ! Obtain sum of weights and valence eigenvalues
      !----------------------------------------------
      s1 = 0.
      seigv = 0.
      DO jspin = 1,jspins
         s = 0.
         DO iBand = 1,MAXVAL(ne(:,jspin))
            DO ikpt = 1,kpts%nkpt
               s = s + w(iBand,ikpt,jspin)
               seigv = seigv + w(iBand,ikpt,jspin)*eig(iBand,ikpt,jspin)
            ENDDO
         ENDDO
         s1 = s1 + s
      ENDDO
      seigv = 2.0/input%jspins*seigv
      chmom = s1 - jspins*s
      IF ( mpi%irank == 0 ) THEN
        WRITE (6,FMT=9300) seigv,s1,chmom
      END IF
9300  FORMAT (/,10x,'sum of valence eigenvalues=',f20.6,5x,&
             'sum of weights=',f10.6,/,10x,'moment=',f12.6)

   END SUBROUTINE fertetra
END MODULE m_fertetra