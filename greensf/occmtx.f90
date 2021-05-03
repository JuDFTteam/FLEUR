MODULE m_occmtx

   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_rotMMPmat

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE occmtx(g,gfinp,input,atoms,noco,nococonv,mmpMat,spin,l_write,check,occError)

      !calculates the occupation of a orbital treated with DFT+HIA from the related greens function
      !The Greens-function should already be prepared on a energy contour ending at e_fermi
      !The occupation is calculated with:
      !
      ! n^sigma_mm' = -1/2pi int^Ef dz (G^+(z)^sigma_mm'-G^-(z)^sigma_mm')
      !
      ! If l_write is given the density matrix together with the spin up/down trace is written to the out files

      TYPE(t_greensf),                  INTENT(IN)    :: g
      TYPE(t_gfinp),                    INTENT(IN)    :: gfinp
      TYPE(t_input),                    INTENT(IN)    :: input
      TYPE(t_atoms),                    INTENT(IN)    :: atoms
      TYPE(t_noco),                     INTENT(IN)    :: noco
      TYPE(t_nococonv),                 INTENT(IN)    :: nococonv
      COMPLEX,                          INTENT(INOUT) :: mmpMat(-lmaxU_const:,-lmaxU_const:,:)
      INTEGER,                 OPTIONAL,INTENT(IN)    :: spin
      LOGICAL,                 OPTIONAL,INTENT(IN)    :: l_write !write the occupation matrix to out file
      LOGICAL,                 OPTIONAL,INTENT(IN)    :: check
      LOGICAL,                 OPTIONAL,INTENT(INOUT) :: occError

      INTEGER :: ind1,ind2,ipm,iz,ispin,l,lp,spin_ind
      INTEGER :: atomType,atomTypep,m,mp,i,j,ns,spin_start,spin_end
      REAL    :: nup,ndwn,tr
      COMPLEX :: weight, offd
      TYPE(t_mat) :: gmat
      CHARACTER(len=300) :: message
      CHARACTER(len=2) :: l_type
      CHARACTER(len=8) :: l_form
      TYPE(t_contourInp) :: contourInp

      !Get the element information
      l  = g%elem%l
      lp = g%elem%lp
      atomType  = g%elem%atomType
      atomTypep = g%elem%atomTypep
      contourInp = gfinp%contour(g%elem%iContour)

      !Check for Contours not reproducing occupations
      IF(contourInp%shape.EQ.CONTOUR_SEMICIRCLE_CONST.AND.ABS(contourInp%et).GT.1e-12) &
         WRITE(oUnit,*) "Energy contour not ending at efermi: These are not the actual occupations"
      IF(contourInp%shape.EQ.CONTOUR_DOS_CONST.AND..NOT.contourInp%l_dosfermi) &
         WRITE(oUnit,*) "Energy contour not weighted for occupations: These are not the actual occupations"

      IF(PRESENT(spin)) THEN
         mmpMat(:,:,spin) = cmplx_0
         spin_start = spin
         spin_end   = spin
      ELSE
         mmpMat = cmplx_0
         spin_start = 1
         IF(ALLOCATED(g%gmmpMat)) spin_end = SIZE(g%gmmpMat,4)
         IF(ALLOCATED(g%uu)) spin_end = SIZE(g%uu,4)
      ENDIF

      DO ispin = spin_start, spin_end
         spin_ind = MERGE(3, ispin, ispin.EQ.4)
         DO ipm = 1, 2
            !Integrate over the contour:
            DO iz = 1, g%contour%nz
               !get the corresponding gf-matrix
               weight = MERGE(g%contour%de(iz),conjg(g%contour%de(iz)),ipm.EQ.1)
               CALL g%get(atoms,iz,ipm.EQ.2,ispin,gmat)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0
                  DO mp = -lp,lp
                     ind2 = ind2 + 1
                     IF(ispin<=3) THEN
                        mmpMat(m,mp,spin_ind) = mmpMat(m,mp,spin_ind) + ImagUnit/tpi_const * (-1)**(ipm-1) * gmat%data_c(ind1,ind2) &
                                                                       * weight
                     ELSE
                        mmpMat(m,mp,spin_ind) = mmpMat(m,mp,spin_ind) - 1.0/tpi_const * (-1)**(ipm-1) * gmat%data_c(ind1,ind2) &
                                                                       * weight
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            !For the contour 3 (real Axis just shifted with sigma) we can add the tails on both ends
            IF(contourInp%shape.EQ.CONTOUR_DOS_CONST.AND.contourInp%l_anacont) THEN
               !left tail
               weight = MERGE(g%contour%de(1),conjg(g%contour%de(1)),ipm.EQ.1)
               CALL g%get(atoms,1,ipm.EQ.2,ispin,gmat)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0
                  DO mp = -lp,lp
                     ind2 = ind2 + 1
                     IF(ispin<=3) THEN
                        mmpMat(m,mp,spin_ind) = mmpMat(m,mp,spin_ind) - 1/tpi_const * gmat%data_c(ind1,ind2) * weight
                     ELSE
                        mmpMat(m,mp,spin_ind) = mmpMat(m,mp,spin_ind) - ImagUnit/tpi_const * gmat%data_c(ind1,ind2) * weight
                     ENDIF
                  ENDDO
               ENDDO
               !right tail
               weight = MERGE(g%contour%de(g%contour%nz),conjg(g%contour%de(g%contour%nz)),ipm.EQ.1)
               CALL g%get(atoms,g%contour%nz,ipm.EQ.2,ispin,gmat)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0
                  DO mp = -lp,lp
                     ind2 = ind2 + 1
                     IF(ispin<=3) THEN
                        mmpMat(m,mp,spin_ind) = mmpMat(m,mp,spin_ind) + 1/tpi_const * gmat%data_c(ind1,ind2) * weight
                     ELSE
                        mmpMat(m,mp,spin_ind) = mmpMat(m,mp,spin_ind) + ImagUnit/tpi_const * gmat%data_c(ind1,ind2) * weight
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      !Rotate the occupation matrix into the global frame in real-space
      IF(noco%l_noco) THEN
         mmpmat(:,:,spin_start:MIN(spin_end,3)) = rotMMPmat(mmpmat(:,:,spin_start:MIN(spin_end,3)),nococonv%alph(atomType),nococonv%beta(atomType),0.0,l)
      ELSE IF(noco%l_soc) THEN
         mmpmat(:,:,spin_start:MIN(spin_end,3)) = rotMMPmat(mmpmat(:,:,spin_start:MIN(spin_end,3)),nococonv%phi,nococonv%theta,0.0,l)
      ENDIF


      !Sanity check are the occupations reasonable?
      IF(PRESENT(check)) THEN
         IF(check) THEN
            IF(PRESENT(occError)) occError = .FALSE.
            DO ispin = spin_start, spin_end
               IF(ispin>input%jspins) CYCLE !Only the spin-diagonal parts
               tr = 0.0
               DO m = -l,l
                  tr = tr + REAL(mmpmat(m,m,ispin))/(3.0-input%jspins)
                  IF(REAL(mmpmat(m,m,ispin))/(3.0-input%jspins).GT. 1.05 .OR.&
                     REAL(mmpmat(m,m,ispin))/(3.0-input%jspins).LT.-0.01) THEN

                     IF(PRESENT(occError)) THEN
                        occError = .TRUE.
                     ELSE
                        WRITE(message,9100) ispin,m,REAL(mmpmat(m,m,ispin))
9100                    FORMAT("Invalid element in mmpmat (spin ",I1,",m ",I2"): ",f14.8)
                        CALL juDFT_warn(TRIM(ADJUSTL(message)),calledby="occmtx")
                     ENDIF
                  ENDIF
               ENDDO
               IF(tr.LT.-0.01.OR.tr.GT.2*l+1.1) THEN
                  IF(PRESENT(occError)) THEN
                     occError = .TRUE.
                  ELSE
                     WRITE(message,9110) ispin,tr
9110                 FORMAT("Invalid occupation for spin ",I1,": ",f14.8)
                     CALL juDFT_warn(TRIM(ADJUSTL(message)),calledby="occmtx")
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      !Io-part (ATM this subroutine is only called from rank 0)
      IF(PRESENT(l_write)) THEN
         IF(l_write) THEN
            !Write to file
            WRITE (l_type,'(i2)') 2*(2*l+1)
            l_form = '('//l_type//'f8.4)'
9000        FORMAT(/,"Occupation matrix obtained from the green's function for atom: ",I3," l: ",I3)
            WRITE(oUnit,9000) atomType, l
            WRITE(oUnit,"(A)") "In the |L,S> basis:"
            DO ispin = 1, MERGE(3, input%jspins, gfinp%l_mperp)
               WRITE(oUnit,'(A,I0)') "Spin: ", ispin
               WRITE(oUnit,l_form) ((mmpmat(i,j,ispin),i=-l,l),j=-lp,lp)
            ENDDO


            IF(l.EQ.lp) THEN
               nup = 0.0
               DO i = -l, l
                  nup = nup + REAL(mmpmat(i,i,1))
               ENDDO
               WRITE(oUnit,'(/,1x,A,I0,A,A,A,f8.4)') "l--> ",l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Up trace: ", nup

               IF(input%jspins.EQ.2) THEN
                  ndwn = 0.0
                  DO i = -l, l
                     ndwn = ndwn + REAL(mmpmat(i,i,2))
                  ENDDO
                  WRITE(oUnit,'(1x,A,I0,A,A,A,f8.4)') "l--> ",l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Down trace: ", ndwn
               ENDIF

               IF(gfinp%l_mperp) THEN
                  offd = cmplx_0
                  DO i = -l, l
                     offd = offd + mmpmat(i,i,3)
                  ENDDO
                  WRITE(oUnit,'(1x,A,I0,A,A,A,f8.4)') "l--> ",l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Offd trace (x): ", REAL(offd)
                  WRITE(oUnit,'(1x,A,I0,A,A,A,f8.4)') "l--> ",l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Offd trace (y): ", AIMAG(offd)
               ENDIF
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE occmtx

END MODULE m_occmtx
