MODULE m_occmtx



   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_ind_greensf
   USE m_lsTOjmj

   IMPLICIT NONE


CONTAINS

   SUBROUTINE occmtx(g,l,nType,atoms,sym,input,mmpMat,lp,nTypep,l_write,check)

      !calculates the occupation of a orbital treated with DFT+HIA from the related greens function
      !The Greens-function should already be prepared on a energy contour ending at e_fermi
      !The occupation is calculated with:
      !
      ! n^sigma_mm' = -1/2pi int^Ef dz (G^+(z)^sigma_mm'-G^-(z)^sigma_mm')
      !
      ! If l_write is given the density matrix together with the spin up/down trace is written to the out files
      ! Additionally the transformation to the |J,mj> subspace is performed via the clebsch gordan coefficients
      ! And the occupations of the respective j states are given 

      TYPE(t_greensf),        INTENT(IN)  :: g
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      TYPE(t_sym),            INTENT(IN)  :: sym
      TYPE(t_input),          INTENT(IN)  :: input
      COMPLEX,                INTENT(OUT) :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,3)
      INTEGER,                INTENT(IN)  :: l
      INTEGER,                INTENT(IN)  :: nType
      INTEGER, OPTIONAL,      INTENT(IN)  :: lp
      INTEGER, OPTIONAL,      INTENT(IN)  :: nTypep
      LOGICAL, OPTIONAL,      INTENT(IN)  :: l_write !write the occupation matrix to out file in both |L,S> and |J,mj>
      LOGICAL, OPTIONAL,      INTENT(IN)  :: check



      INTEGER ind1,ind2,ipm,iz,ispin,m,mp,lp_loop,i,ns
      LOGICAL l_vertcorr
      REAL    re,imag,nup,ndwn,nhi,nlow,tr
      TYPE(t_mat) :: gmat,cmat,jmat
      CHARACTER(len=300) :: message

      l_vertcorr = .false. !Enables/Disables a correction for the vertical parts of the rectangular contour

      mmpMat = 0.0

      IF(.NOT.PRESENT(lp)) THEN
         lp_loop = l 
      ELSE 
         lp_loop = lp 
      ENDIF

      DO ispin = 1, MERGE(3,input%jspins,input%l_gfmperp)
         DO ipm = 1, 2
            !Integrate over the contour:
            DO iz = 1, g%nz
               !get the corresponding gf-matrix
               CALL g%get_gf(gmat,atoms,input,iz,l,nType,ipm.EQ.2,spin=ispin,lp=lp,nTypep=nTypep)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0 
                  DO mp = -lp_loop,lp_loop
                     ind2 = ind2 + 1 
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) - 1/(2.0*pi_const*ImagUnit) * (-1)**(ipm-1) * gmat%data_c(ind1,ind2) &
                                                             * MERGE(g%de(iz),conjg(g%de(iz)),ipm.EQ.1)
                  ENDDO
               ENDDO
               CALL gmat%free()
            ENDDO
            !For the contour 3 (real Axis just shifted with sigma) we can add the tails on both ends
            IF(g%mode.EQ.3.AND.input%gf_anacont) THEN
               !left tail
               CALL g%get_gf(gmat,atoms,input,1,l,nType,ipm.EQ.2,spin=ispin,lp=lp,nTypep=nTypep)
               ind1 = 0
               DO m = -l, l 
                  ind1 = ind1 + 1
                  ind2 = 0 
                  DO mp = -lp_loop,lp_loop
                     ind2 = ind2 + 1 
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) - 1/(2.0*pi_const*ImagUnit) * (-1)**(ipm-1) * gmat%data_c(ind1,ind2) &
                                                             * MERGE(g%de(1),conjg(g%de(1)),ipm.EQ.1)
                  ENDDO
               ENDDO
               CALL gmat%free()
               !right tail
               CALL g%get_gf(gmat,atoms,input,g%nz,l,nType,ipm.EQ.2,spin=ispin,lp=lp,nTypep=nTypep)
               ind1 = 0
               DO m = -l, l 
                  ind1 = ind1 + 1
                  ind2 = 0 
                  DO mp = -lp_loop,lp_loop
                     ind2 = ind2 + 1 
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) + 1/(2.0*pi_const*ImagUnit) * (-1)**(ipm-1) * gmat%data_c(ind1,ind2) &
                                                             * MERGE(g%de(g%nz),conjg(g%de(g%nz)),ipm.EQ.1)
                  ENDDO
               ENDDO
               CALL gmat%free()
            ENDIF
         ENDDO
      ENDDO

      !Sanity check are the occupations reasonable?
      IF(PRESENT(check)) THEN
         IF(check) THEN
            DO ispin = 1, input%jspins 
               tr = 0.0
               DO i = -l,l
                  tr = tr + REAL(mmpmat(i,i,ispin))/(3-input%jspins)
                  IF(REAL(mmpmat(i,i,ispin))/(3-input%jspins).GT.1.01&
                     .OR.REAL(mmpmat(i,i,ispin))/(3-input%jspins).LT.0.0)&
                        CALL juDFT_warn("Invalid element in mmpMat")
               ENDDO
               IF(tr.LT.0.OR.tr.GT.2*l+1.1) THEN
                  WRITE(message,9100) ispin,tr
9100              FORMAT("Invalid occupation for spin ",I1,": ",f14.8)
                  CALL juDFT_warn(TRIM(ADJUSTL(message)),calledby="occmtx") 
               ENDIF
            ENDDO
         ENDIF
      ENDIF


      !Io-part (ATM this subroutine is only called from rank 0)
      IF(PRESENT(l_write)) THEN
         IF(l_write.AND.lp_loop.EQ.l.AND.(.NOT.PRESENT(nTypep).OR.(PRESENT(nTypep).AND.nTypep.EQ.nType))) THEN
            !Construct the full matrix in the |L,ml,ms> basis (real)
            ns = 2*l+1
            CALL gmat%init(.TRUE.,2*ns,2*ns)
            CALL jmat%init(.TRUE.,2*ns,2*ns)
            DO m = -l, l 
               DO mp = -l, l
                  gmat%data_r(m+l+1,mp+l+1) = REAL(mmpmat(m,mp,1))/(3-input%jspins)
                  IF(input%jspins.EQ.1) THEN
                     gmat%data_r(m+l+1+ns,mp+l+1+ns) = REAL(mmpmat(-m,-mp,MIN(2,input%jspins)))/(3-input%jspins)
                  ELSE
                     gmat%data_r(m+l+1+ns,mp+l+1+ns) = REAL(mmpmat(m,mp,MIN(2,input%jspins)))/(3-input%jspins)
                  ENDIF
               ENDDO
            ENDDO
            !spin-offdiagonal
            WRITE(*,"(14f14.8)") mmpmat(-l:l,-l:l,3)
            IF(input%l_gfmperp) THEN
               gmat%data_r(1:ns,ns+1:2*ns) = REAL(mmpmat(-l:l,-l:l,3))
               gmat%data_r(ns+1:2*ns,1:ns) = REAL(transpose(mmpmat(-l:l,-l:l,3)))
            ENDIF
            !Calculate the spin-up/down occupation
            nup = 0.0
            DO i = 1, ns
               nup = nup + gmat%data_r(i,i)
            ENDDO
            ndwn = 0.0
            DO i = ns+1, 2*ns
               ndwn = ndwn + gmat%data_r(i,i)
            ENDDO
            !Write to file
            WRITE(6,*)
9000        FORMAT("Occupation matrix obtained from the green's function for atom: ",I3," l: ",I3)
            WRITE(6,9000) nType, l
            WRITE(6,"(A)") "In the |L,S> basis:"
            WRITE(6,"(14f8.4)") gmat%data_r
            WRITE(6,"(1x,A,f8.4)") "Spin-Up trace: ", nup
            WRITE(6,"(1x,A,f8.4)") "Spin-Down trace: ", ndwn

            !Obtain the conversion matrix to the |J,mj> basis
            CALL cmat%init(.TRUE.,2*ns,2*ns)
            CALL lsTOjmj(cmat,l)
            !Perform the transformation
            jmat%data_r = matmul(gmat%data_r,cmat%data_r)
            jmat%data_r = matmul(transpose(cmat%data_r),jmat%data_r)
            !Calculate the low/high j trace
            nlow = 0.0
            DO i = 1, ns-1
               nlow = nlow + jmat%data_r(i,i)
            ENDDO
            nhi = 0.0
            DO i = ns, 2*ns
               nhi = nhi + jmat%data_r(i,i)
            ENDDO
            
            !Write to file
            WRITE(6,"(A)") "In the |J,mj> basis:"
            WRITE(6,"(14f8.4)") jmat%data_r
            WRITE(6,"(1x,A,f8.4)") "Low J trace: ", nlow
            WRITE(6,"(1x,A,f8.4)") "High J trace: ", nhi
            WRITE(6,*)
         ENDIF
      ENDIF

   END SUBROUTINE occmtx



END MODULE m_occmtx