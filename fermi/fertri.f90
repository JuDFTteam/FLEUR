!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fertri
   !
   !     calculates fermi energy and weights using triangular (tetrahedron) method
   !
   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_tetraef
   USE m_dosef
   USE m_dosint
   USE m_doswt
   USE m_xmlOutput

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE fertri(input,noco,kpts,irank,ne,jspins,zc,eig,sfac,&
                     ef,seigv,w)

      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_noco),    INTENT(IN)    :: noco
      TYPE(t_kpts),  INTENT(IN)    :: kpts
      INTEGER,       INTENT(IN)    :: jspins,irank
      REAL,          INTENT(IN)    :: zc,sfac
      INTEGER,       INTENT(IN)    :: ne(:,:)!(nkpt,jspins)
      REAL,          INTENT(OUT)   :: seigv
      REAL,          INTENT(INOUT) :: ef
      REAL,          INTENT(OUT)   :: w(:,:,:) !(neig,nkpt,jspins)
      REAL,          INTENT(INOUT) :: eig(:,:,:)!(neig,nkpt,jspins)

      REAL    :: chmom,ct,del,dez,ei,emax,emin,s,s1,workf
      REAL    :: lb,ub,e_set,seigvTemp
      INTEGER :: i,ic,j,jsp,k,neig,jj
      INTEGER :: nemax(2)
      CHARACTER(LEN=20)    :: attributes(2)
      REAL, PARAMETER :: de = 5.0e-3 !Step for initial search

      IF (irank == 0) THEN
         WRITE (oUnit,FMT=8000)
8000     FORMAT (/,/,10x,'linear triangular method')
      END IF
      w=0.0
      !
      !--->   clear w and set eig=-9999.9
      e_set = -9999.9
      IF (.NOT.input%film) e_set = 1.0e10
      DO jsp = 1,jspins
         nemax(jsp) = 0.0
         DO k = 1,kpts%nkpt
            nemax(jsp) = max0(nemax(jsp),ne(k,jsp))
            DO i = 1,ne(k,jsp)
               w(i,k,jsp) = 0.
            ENDDO
            DO i = ne(k,jsp)+1,size(w,1)
               w(i,k,jsp) = 0.
               eig(i,k,jsp) = e_set
            ENDDO
         ENDDO
      ENDDO
      !
      !      sfac = 2.0/real(jspins)

      IF(.not.input%film) THEN
         lb = MINVAL(eig) - 0.01
         ub = ef + 0.2
         CALL tetra_ef(kpts,jspins,lb,ub,eig,zc,sfac,ef,w)
      ELSE
         IF (irank == 0) THEN
            WRITE (oUnit,FMT=*) 'ef_hist=',ef
         END IF

         ei = ef
!jr      emin = -9999.9
         emin = +9999.9
         emax = -emin
         !Find initial boundaries for the bisection method
         ic = 1
         DO WHILE (emin.GT.emax)
            ic = ic + 1

            CALL dosint(ei,nemax,jspins,kpts,sfac,eig,ct)

            IF (irank == 0) WRITE (oUnit,FMT=*) 'ct=',ct

            IF (ct.LT.zc) THEN ! ei < ef
               emin = ei
               ei = ei + de
            ELSEIF (ct.GT.zc) THEN ! ei > ef
               emax = ei
               ei = ei - de
            ENDIF

            IF (irank==0 .AND. ic.GT.100) THEN
               WRITE (oUnit,FMT=8050) ei,ef,emin,emax,ct,zc
8050           FORMAT (/,/,10x,'error fertri: initial guess of ef off by 25 mry',&
                        ' ei,ef,emin,emax,ct,zc',/,10x,6e16.7,/,10x,&
                        'check number of bands')
               CALL juDFT_error("initial guess of ef off by 25 mry",calledby="fertri")
            END IF

         ENDDO

         IF (ct.NE.zc) THEN
            IF (irank == 0) WRITE (oUnit,FMT=*) '2nd dosint'
            !---> refine ef to a value of 5 mry * (2**-20)
            iterate : DO i = 1, 40
               ei = 0.5* (emin+emax)

               CALL dosint(ei,nemax,jspins,kpts,sfac,eig,ct)

               IF (irank == 0) WRITE (oUnit,FMT=*) 'i=',i,', ct=',ct
               IF ( ct == zc ) THEN
                  EXIT iterate
               ELSEIF ( ct > zc ) THEN
                  emax = ei
               ELSE
                  emin = ei
               ENDIF
            ENDDO iterate
         ENDIF
         ef = ei
         del = emax - emin
         dez = zc - ct
         workf = -hartree_to_ev_const*ef
         IF (irank == 0) THEN
            WRITE (oUnit,FMT=8030) ef,workf,del,dez
8030        FORMAT(/,10x,'fermi energy=',f10.5,' har',/,10x,'work function='&
                    ,f10.5,' ev',/,10x,'uncertainity in energy and weights=',&
                     2e16.6)
         END IF
         !
         !--->   obtain dos at ef
         !
         CALL dosef(ei,nemax,jspins,kpts,sfac,eig)
         !
         !--->   obtain weights needed for integration
         !
         CALL doswt(ei,nemax,jspins,kpts,eig,w)

      ENDIF ! .NOT.input%film

      !
      !--->   write weights
      !
!     DO jsp = 1,jspins
!        neig = nemax(jsp)
!        DO i = 1,neig
!           DO k = 1,kpts%nkpt
!              WRITE (oUnit,FMT=*) 'w(',i,',',k,',',jsp,')=',w(i,k,jsp)
!           ENDDO
!        ENDDO
!     ENDDO

      !find degenerate states
     DO jsp=1,jspins
         DO k=1,kpts%nkpt
            i=1   
            DO while(i<nemax(jsp))
               j=1
               do while (abs(eig(i,k,jsp)-eig(i+j,k,jsp))<1E-9)
                  j=j+1
                  if (i+j>nemax(jsp)) exit
               ENDDO
               if (j>1) THEN
                  j=j-1
                  !Make sure all weights are equal
                  w(i:i+j,k,jsp)=sum(w(i:i+j,k,jsp))/j
                  i=i+j
               endif
               i=i+1   
            enddo      
         enddo
      enddo
      !
      !--->   obtain sum of weights and valence eigenvalues
      !
      
      s1 = 0.
      seigv = 0.
      DO jsp = 1,jspins
         s = 0.
         neig = nemax(jsp)
         DO i = 1,neig
            DO k = 1,kpts%nkpt
               s = s + w(i,k,jsp)
               seigv = seigv + w(i,k,jsp)*eig(i,k,jsp)
            ENDDO
         ENDDO
         s1 = s1 + s
      ENDDO
      seigv = sfac*seigv
      chmom = s1 - jspins*s

      seigvTemp = seigv
      IF (noco%l_soc .AND. (.NOT. noco%l_noco)) THEN
         seigvTemp = seigvTemp / 2.0
      END IF

      IF (irank == 0) THEN
         attributes = ''
         WRITE(attributes(1),'(f20.10)') seigvTemp
         WRITE(attributes(2),'(a)') 'Htr'
         CALL writeXMLElement('sumValenceSingleParticleEnergies',(/'value','units'/),attributes)
         WRITE (oUnit,FMT=8040) seigvTemp,s1,chmom
8040     FORMAT (/,10x,'sum of valence eigenvalues=',f20.10,5x,&
                  'sum of weights=',f10.6,/,10x,'moment=',f12.6)
      END IF

   END SUBROUTINE fertri
END MODULE m_fertri
