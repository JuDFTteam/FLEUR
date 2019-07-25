MODULE m_occmtx



   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_ind_greensf
   USE m_lsTOjmj

   IMPLICIT NONE


CONTAINS

   SUBROUTINE occmtx(g,l,nType,atoms,sym,input,mmpMat,lp,nTypep,l_write)

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
      COMPLEX,                INTENT(OUT) :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MERGE(3,input%jspins,input%l_gfmperp))
      INTEGER,                INTENT(IN)  :: l
      INTEGER,                INTENT(IN)  :: nType
      INTEGER, OPTIONAL,      INTENT(IN)  :: lp
      INTEGER, OPTIONAL,      INTENT(IN)  :: nTypep
      LOGICAL, OPTIONAL,      INTENT(IN)  :: l_write !write the occupation matrix to out file in both |L,S> and |J,mj>

      INTEGER ind1,ind2,ipm,iz,ispin,m,mp,lp_loop,i,ns
      LOGICAL l_vertcorr
      REAL    re,imag,nup,ndwn,nhi,nlow
      TYPE(t_mat) :: gmat,cmat,jmat

      l_vertcorr = .false. !Enables/Disables a correction for the vertical parts of the rectangular contour

      mmpMat(:,:,:) = CMPLX(0.0,0.0)

      IF(.NOT.PRESENT(lp)) THEN
         lp_loop = l 
      ELSE 
         lp_loop = lp 
      ENDIF

      !REPLACE: input%jspins --> MERGE(3,input%jspins,input%l_gfmperp)
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
         ENDDO
      ENDDO

      !Io-part (ATM this subroutine is only called from rank 0)
      IF(l_write.AND.lp_loop.EQ.l.AND.(.NOT.PRESENT(nTypep).OR.(PRESENT(nTypep).AND.nTypep.EQ.nType))) THEN
         !Construct the full matrix in the |L,ml,ms> basis (real)
         ns = 2*l+1
         CALL gmat%init(.TRUE.,2*ns,2*ns)
         CALL jmat%init(.TRUE.,2*ns,2*ns)
         !spin-up
         gmat%data_r(1:ns,1:ns) = REAL(mmpmat(-l:l,-l:l,1))
         !spin-down
         gmat%data_r(ns+1:2*ns,ns+1:2*ns) = REAL(mmpmat(-l:l,-l:l,2))
         !spin-offdiagonal
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
9000     FORMAT("Occupation matrix obtained from the green's function for atom: ",I3," l: ",I3)
         WRITE(6,9000) nType, l
         WRITE(6,"(A)") "In the |L,S> basis:"
         WRITE(6,"(14f8.5)") gmat%data_r
         WRITE(6,"(1x,A,f8.5)") "Spin-Up trace: ", nup
         WRITE(6,"(1x,A,f8.5)") "Spin-Down trace: ", ndwn

         !Obtain the conversion matrix to the |J,mj> basis
         CALL cmat%init(.TRUE.,2*ns,2*ns)
         CALL lsTOjmj(cmat,l)
         !Perform the transformation
         jmat%data_r = matmul(gmat%data_r,cmat%data_r)
         jmat%data_r = matmul(transpose(cmat%data_r),jmat%data_r)
         !Calculate the low/high j trace
         nlow = 0.0
         DO i = 1, ns-1
            nlow = nlow + gmat%data_r(i,i)
         ENDDO
         nhi = 0.0
         DO i = ns, 2*ns
            nhi = nhi + gmat%data_r(i,i)
         ENDDO

         WRITE(6,"(A)") "In the |J,mj> basis:"
         WRITE(6,"(14f8.5)") jmat%data_r
         WRITE(6,"(1x,A,f8.5)") "Low J trace: ", nlow
         WRITE(6,"(1x,A,f8.5)") "High J trace: ", nhi
      ENDIF

   END SUBROUTINE occmtx



END MODULE m_occmtx