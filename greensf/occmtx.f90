MODULE m_occmtx

   USE m_juDFT
   USE m_types
   USE m_types_scalarGF
   USE m_constants
   USE m_rotMMPmat

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE occmtx(g,gfinp,input,atoms,noco,nococonv,mmpMat,spin,usdus,denCoeffsOffDiag,scalarGF,l_write,check,occError)

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
      TYPE(t_usdus),           OPTIONAL,INTENT(IN)    :: usdus
      TYPE(t_denCoeffsOffDiag),OPTIONAL,INTENT(IN)    :: denCoeffsOffDiag
      TYPE(t_scalarGF),        OPTIONAL,INTENT(IN)    :: scalarGF
      LOGICAL,                 OPTIONAL,INTENT(IN)    :: l_write !write the occupation matrix to out file
      LOGICAL,                 OPTIONAL,INTENT(IN)    :: check
      LOGICAL,                 OPTIONAL,INTENT(INOUT) :: occError

      INTEGER :: ind1,ind2,ipm,iz,ispin,l,lp
      INTEGER :: atomType,atomTypep,m,mp,i,j,ns,spin_start,spin_end
      REAL    :: nup,ndwn,tr
      COMPLEX :: weight
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
         DO ipm = 1, 2
            !Integrate over the contour:
            DO iz = 1, g%contour%nz
               !get the corresponding gf-matrix
               weight = MERGE(g%contour%de(iz),conjg(g%contour%de(iz)),ipm.EQ.1)
               CALL g%get(atoms,iz,ipm.EQ.2,spin,gmat,usdus=usdus,&
                          denCoeffsOffDiag=denCoeffsOffDiag,scalarGF=scalarGF)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0
                  DO mp = -lp,lp
                     ind2 = ind2 + 1
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) + ImagUnit/tpi_const * (-1)**(ipm-1) * gmat%data_c(ind1,ind2) &
                                                             * weight
                  ENDDO
               ENDDO
            ENDDO
            !For the contour 3 (real Axis just shifted with sigma) we can add the tails on both ends
            IF(contourInp%shape.EQ.CONTOUR_DOS_CONST.AND.contourInp%l_anacont) THEN
               !left tail
               weight = MERGE(g%contour%de(1),conjg(g%contour%de(1)),ipm.EQ.1)
               CALL g%get(atoms,1,ipm.EQ.2,spin,gmat,usdus=usdus,&
                          denCoeffsOffDiag=denCoeffsOffDiag,scalarGF=scalarGF)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0
                  DO mp = -lp,lp
                     ind2 = ind2 + 1
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) - 1/tpi_const * gmat%data_c(ind1,ind2) * weight
                  ENDDO
               ENDDO
               !right tail
               weight = MERGE(g%contour%de(g%contour%nz),conjg(g%contour%de(g%contour%nz)),ipm.EQ.1)
               CALL g%get(atoms,g%contour%nz,ipm.EQ.2,spin,gmat,usdus=usdus,&
                          denCoeffsOffDiag=denCoeffsOffDiag,scalarGF=scalarGF)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0
                  DO mp = -lp,lp
                     ind2 = ind2 + 1
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) + 1/tpi_const * gmat%data_c(ind1,ind2) * weight
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      !Rotate the occupation matrix into the global frame in real-space
      IF(noco%l_noco) THEN
         mmpmat(:,:,spin_start:spin_end) = rotMMPmat(mmpmat(:,:,spin_start:spin_end),nococonv%alph(atomType),nococonv%beta(atomType),0.0,l)
      ELSE IF(noco%l_soc) THEN
         mmpmat(:,:,spin_start:spin_end) = rotMMPmat(mmpmat(:,:,spin_start:spin_end),nococonv%phi,nococonv%theta,0.0,l)
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
            !Construct the full matrix in the |L,ml,ms> basis (real)
            ns = 2*l+1
            CALL gmat%init(.FALSE.,2*ns,2*ns)
            DO m = -l, l
               DO mp = -l, l
                  gmat%data_c(m+l+1,mp+l+1) = mmpmat(m,mp,1)/(3-input%jspins)
                  IF(input%jspins.EQ.1) THEN
                     gmat%data_c(m+l+1+ns,mp+l+1+ns) = mmpmat(-m,-mp,1)/(3-input%jspins)
                  ELSE
                     gmat%data_c(m+l+1+ns,mp+l+1+ns) = mmpmat(m,mp,2)
                  ENDIF
               ENDDO
            ENDDO
            !spin-offdiagonal
            IF(gfinp%l_mperp) THEN
               gmat%data_c(1:ns,ns+1:2*ns) = mmpmat(-l:l,-l:l,3)
               gmat%data_c(ns+1:2*ns,1:ns) = conjg(transpose(mmpmat(-l:l,-l:l,3)))
            ENDIF
            !Calculate the spin-up/down occupation
            nup = 0.0
            DO i = 1, ns
               nup = nup + REAL(gmat%data_c(i,i))
            ENDDO
            ndwn = 0.0
            DO i = ns+1, 2*ns
               ndwn = ndwn + REAL(gmat%data_c(i,i))
            ENDDO
            !Write to file
            WRITE (l_type,'(i2)') 4*ns
            l_form = '('//l_type//'f8.4)'
9000        FORMAT(/,"Occupation matrix obtained from the green's function for atom: ",I3," l: ",I3)
            WRITE(oUnit,9000) atomType, l
            WRITE(oUnit,"(A)") "In the |L,S> basis:"
            WRITE(oUnit,l_form) ((gmat%data_c(i,j),i=1,2*ns),j=1,2*ns)
            WRITE(oUnit,'(1x,A,I0,A,A,A,f8.4)') "l--> ",l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Up trace: ", nup
            WRITE(oUnit,'(1x,A,I0,A,A,A,f8.4)') "l--> ",l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Down trace: ", ndwn
         ENDIF
      ENDIF

   END SUBROUTINE occmtx

END MODULE m_occmtx
