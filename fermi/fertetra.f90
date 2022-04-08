MODULE m_fertetra

   USE m_types
   USE m_constants
   USE m_juDFT
   USE m_tetrahedronInit
   USE m_xmlOutput

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
      REAL    :: dlow,dup,dfermi,s1,s,chmom,seigvTemp
      REAL    :: lowBound,upperBound,weightSum
      CHARACTER(LEN=20)    :: attributes(2)


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
         CALL tetrahedronInit(kpts,input,eig(:,:,jspin),MINVAL(ne(:,jspin)),&
                              lowBound,weightSum=weightSum)
         dlow = dlow + weightSum * 2.0/input%jspins
      ENDDO

      IF (dlow.GT.input%zelec) THEN
         WRITE(oUnit,9000) lowBound,dlow,input%zelec
         CALL juDFT_error("valence band too high ",calledby="fertetra")
      ENDIF
9000  FORMAT (' valence band too high ',/,&
              '  elow ',f10.5,' dlow ',f10.5,' nelec ',f10.5)

      it = 0
      DO
         !Now check the upper bound
         dup = 0.0
         DO jspin = 1, jspins
            CALL tetrahedronInit(kpts,input,eig(:,:,jspin),MINVAL(ne(:,jspin)),&
                                 upperBound,weightSum=weightSum)
            dup = dup + weightSum * 2.0/input%jspins
         ENDDO

         IF (dup.GT.input%zelec) THEN
            EXIT
         ELSE
            !Raise the upper bound
            upperBound = upperBound + 0.2
            it = it + 1
            IF(it.GT.10) THEN
               WRITE (oUnit,9100) upperBound,dup,input%zelec
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
            ! Compute the current occupation
            !-------------------------------------------------------
            CALL tetrahedronInit(kpts,input,eig(:,:,jspin),MINVAL(ne(:,jspin)),&
                                 ef,weightSum=weightSum)
            dfermi = dfermi + weightSum * 2.0/input%jspins
         ENDDO
         IF(ABS(dfermi-input%zelec).LT.1e-12) THEN
            EXIT
         ELSE IF(dfermi-input%zelec.GT.0.0) THEN
            !Occupation to large -> search in the left interval
            upperBound = ef
         ELSE
            !Occupation to small -> search in the right interval
            lowBound = ef
         ENDIF
      ENDDO

      !---------------------------------------------------------------------
      !Calculate final occupation and weights for individual kpts and bands
      !---------------------------------------------------------------------
      ef = (lowBound+upperBound)/2.0
      dfermi = 0.0
      DO jspin = 1, jspins
         !-------------------------------------------------------
         ! Compute the weights for charge density integration
         !-------------------------------------------------------
         CALL tetrahedronInit(kpts,input,eig(:,:,jspin),MINVAL(ne(:,jspin)),&
                              ef,weightSum=weightSum,weights=w(:,:,jspin))
         dfermi = dfermi + weightSum * 2.0/input%jspins
      ENDDO

      WRITE (oUnit,9200) ef,dfermi,input%zelec
9200  FORMAT (//,'Tetrahedron method: ',//,'   fermi energy =',f10.5,&
                 ' dtot ',f10.5,' nelec ',f10.5)

      !----------------------------------------------
      ! Obtain sum of weights and valence eigenvalues
      !----------------------------------------------
      s1 = 0.0
      seigv = 0.0
      DO jspin = 1,jspins
         s = 0.0
         DO ikpt = 1,kpts%nkpt
            DO iBand = 1,ne(ikpt,jspin)
               s = s + w(iBand,ikpt,jspin)
               seigv = seigv + w(iBand,ikpt,jspin)*eig(iBand,ikpt,jspin)
            ENDDO
         ENDDO
         s1 = s1 + s
      ENDDO
      seigv = 2.0/input%jspins*seigv
      chmom = s1 - jspins*s
      
      seigvTemp = seigv
      IF (noco%l_soc .AND. (.NOT. noco%l_noco)) THEN
         seigvTemp = seigvTemp / 2.0
      END IF
      
      IF ( mpi%irank == 0 ) THEN
         attributes = ''
         WRITE(attributes(1),'(f20.10)') seigvTemp
         WRITE(attributes(2),'(a)') 'Htr'
         CALL writeXMLElement('sumValenceSingleParticleEnergies',(/'value','units'/),attributes)
         WRITE (oUnit,FMT=9300) seigvTemp,s1,chmom
      END IF
9300  FORMAT (/,10x,'sum of valence eigenvalues=',f20.10,5x,&
             'sum of weights=',f10.6,/,10x,'moment=',f12.6)

   END SUBROUTINE fertetra
END MODULE m_fertetra
