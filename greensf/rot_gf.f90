MODULE m_rot_gf

!-----------------------------------------------
! Rotates the contribution from eqivalent atoms
!-----------------------------------------------

CONTAINS


   SUBROUTINE rot_projDOS(sym,atoms,input,angle,greensfCoeffs)

      USE m_types
      USE m_juDFT
      USE m_constants

      IMPLICIT NONE

      TYPE(t_sym),            INTENT(IN)     :: sym
      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_input),          INTENT(IN)     :: input
      TYPE(t_greensfCoeffs),  INTENT(INOUT)  :: greensfCoeffs
      REAL,                   INTENT(IN)     :: angle(:)

      COMPLEX, ALLOCATABLE :: curr_dos(:,:,:),calc_mat(:,:,:)
      COMPLEX d_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

      INTEGER i_gf,l,nType,nn,natom,ispin,imat,it,is,isi,m,mp,ie
      REAL fac
      COMPLEX phase

      CALL timestart("Green's function: Rotate")

      ALLOCATE(curr_dos(greensfCoeffs%ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),&
               calc_mat(greensfCoeffs%ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))

      DO i_gf = 1, atoms%n_gf

         l     = atoms%gfelem(i_gf)%l
         nType = atoms%gfelem(i_gf)%atomType

         !Loop through equivalent atoms
         DO nn = 1, atoms%neq(nType)
            natom = SUM(atoms%neq(:nType-1)) + nn
            !Rotate the eqivalent atom into the irreducible brillouin zone
            fac = 1.0/(sym%invarind(natom)*atoms%neq(nType))
            IF(sym%invarind(natom).EQ.0) CALL juDFT_error("No symmetry operations available",calledby="greensfImag")
            DO ispin = 1, MERGE(3,input%jspins,input%l_gfmperp)
               DO imat = 1, MERGE(1,5,input%l_gfsphavg)
                  IF(imat.EQ.1) THEN
                     curr_dos(:,:,:) = greensfCoeffs%projdos(:,:,:,nn,i_gf,ispin)
                  ELSE IF(imat.EQ.2) THEN
                     curr_dos(:,:,:) = greensfCoeffs%uu(:,:,:,nn,i_gf,ispin)
                  ELSE IF(imat.EQ.3) THEN
                     curr_dos(:,:,:) = greensfCoeffs%dd(:,:,:,nn,i_gf,ispin)
                  ELSE IF(imat.EQ.4) THEN
                     curr_dos(:,:,:) = greensfCoeffs%ud(:,:,:,nn,i_gf,ispin)
                  ELSE IF(imat.EQ.5) THEN
                     curr_dos(:,:,:) = greensfCoeffs%du(:,:,:,nn,i_gf,ispin)
                  ENDIF
                  DO m = -l ,l
                     IF(ANY(AIMAG(curr_dos(:,m,m)).GT.0.0).AND.ispin<3) CALL juDFT_error("curr_dos>0")
                  ENDDO
                  DO it = 1, sym%invarind(natom)
                     is = sym%invarop(natom,it)
                     isi = sym%invtab(is)
                     d_mat(:,:) = cmplx(0.0,0.0)
                     DO m = -l,l
                        DO mp = -l,l
                           d_mat(m,mp) = sym%d_wgn(m,mp,l,isi)
                        ENDDO
                     ENDDO
                     phase = MERGE(exp(ImagUnit*angle(isi)),CMPLX(1.0,0.0),ispin.EQ.3)
                     DO ie = 1, greensfCoeffs%ne
                        calc_mat(ie,:,:) = matmul( transpose( conjg(d_mat) ) , curr_dos(ie,:,:))
                        calc_mat(ie,:,:) = matmul( calc_mat(ie,:,:), d_mat )
                     ENDDO
                     DO ie = 1, greensfCoeffs%ne
                        IF(imat.EQ.1) THEN
                           greensfCoeffs%projdos(ie,:,:,0,i_gf,ispin) = greensfCoeffs%projdos(ie,:,:,0,i_gf,ispin) + phase *AIMAG(fac *  calc_mat(ie,:,:))
                        ELSE IF(imat.EQ.2) THEN
                           greensfCoeffs%uu(ie,:,:,0,i_gf,ispin) = greensfCoeffs%uu(ie,:,:,0,i_gf,ispin) + AIMAG(fac * phase * calc_mat(ie,:,:))
                        ELSE IF(imat.EQ.3) THEN
                           greensfCoeffs%dd(ie,:,:,0,i_gf,ispin) = greensfCoeffs%dd(ie,:,:,0,i_gf,ispin) + AIMAG(fac * phase * calc_mat(ie,:,:))
                        ELSE IF(imat.EQ.4) THEN
                           greensfCoeffs%ud(ie,:,:,0,i_gf,ispin) = greensfCoeffs%ud(ie,:,:,0,i_gf,ispin) + AIMAG(fac * phase * calc_mat(ie,:,:))
                        ELSE IF(imat.EQ.5) THEN
                           greensfCoeffs%du(ie,:,:,0,i_gf,ispin) = greensfCoeffs%du(ie,:,:,0,i_gf,ispin) + AIMAG(fac * phase * calc_mat(ie,:,:))
                        ENDIF
                     ENDDO!ie
                  ENDDO!it
               ENDDO!imat
            ENDDO
         ENDDO
      ENDDO

      CALL timestop("Green's function: Rotate")


   END SUBROUTINE rot_projDOS

END MODULE m_rot_gf