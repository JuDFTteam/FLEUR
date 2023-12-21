!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_vmix
   !
   ! mix the old and new density matrix for the lda+U method
   !                                                 gb.2001
   ! --------------------------------------------------------
   ! Extension to multiple U per atom type by G.M. 2017
   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_xmlOutput

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE v_mix(input,atoms,noco,nIJ_llp_mmp_in,nIJ_llp_mmp_out)

      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_noco), INTENT(IN)    :: noco
      COMPLEX,      INTENT(IN)    :: nIJ_llp_mmp_out(-lmaxU_const:,-lmaxU_const:,:,:)
      COMPLEX,      INTENT(INOUT) :: nIJ_llp_mmp_in (-lmaxU_const:,-lmaxU_const:,:,:)


      INTEGER :: jsp, i_pair, i_v, latom1, latom2, matom1, matom2, natom2
!      REAL    :: alpha,spinf,gam,del,uParam,jParam
      REAL    :: dist(SIZE(nIJ_llp_mmp_in,4))

      CHARACTER(LEN=20)   :: attributes(6)
!      COMPLEX,ALLOCATABLE :: nIJ_llp_mmp(:,:,:,:)

      !
      ! check for possible rotation of n_mmp
      !
      !zero=0.0
      !CALL nmat_rot(zero,-atoms%lda_u%theta,-atoms%lda_u%phi,3,atoms%n_u,input%jspins,atoms%lda_u%l,n_mmp_out)

      ! Write out n_mmp_out to out.xml file
      CALL openXMLElementNoAttributes('ldaVDensityMatrix')
!      DO jsp = 1, SIZE(n_mmp_out,4)
!         DO i_u = 1, atoms%n_u
!            l = atoms%lda_u(i_u)%l
!            itype = atoms%lda_u(i_u)%atomType
!            uParam = atoms%lda_u(i_u)%u
!            jParam = atoms%lda_u(i_u)%j
!            attributes = ''
!            WRITE(attributes(1),'(i0)') jsp
!            WRITE(attributes(2),'(i0)') itype
!            WRITE(attributes(3),'(i0)') i_u
!            WRITE(attributes(4),'(i0)') l
!            WRITE(attributes(5),'(f15.8)') uParam
!            WRITE(attributes(6),'(f15.8)') jParam
!            CALL writeXMLElementMatrixPoly('densityMatrixFor',&
!                                          (/'spin    ','atomType','uIndex  ','l       ','U       ','J       '/),&
!                                          attributes,n_mmp_out(-l:l,-l:l,i_u,jsp))
!         END DO
!      END DO
      CALL closeXMLElement('ldaVDensityMatrix')

      ! exit subroutine if density matrix does not exist
      IF(.NOT.ANY(ABS(nIJ_llp_mmp_in(:,:,:,:)).GT.1e-12)) RETURN

      !Calculate distance
      dist = 0.0
      DO jsp = 1, SIZE(nIJ_llp_mmp_in,4)
         i_pair=0 !counts number of pairs
         DO i_v = 1, atoms%n_v  !loop over pairs which are corrected by U+V 
            latom1 = atoms%lda_v(i_v)%thisAtomL
            Do natom2 = 1, atoms%lda_v(i_v)%numOtherAtoms
               i_pair = i_pair + 1
               latom2 = atoms%lda_v(i_v)%otherAtomL
               Do matom1 = -latom1, latom1
                  Do matom2 = -latom2, latom2
                     dist(jsp) = dist(jsp) + ABS(nIJ_llp_mmp_out(matom1,matom2,i_pair,jsp) - nIJ_llp_mmp_in(matom1,matom2,i_pair,jsp))
                  END DO
               END DO
            END DO
         END DO
      END DO

      !Write to outfile
      IF(input%jspins.EQ.1) THEN
         WRITE (oUnit,'(a,f12.6)') 'nIJ_llp_mmp distance =',dist(1)
      ELSE
         DO jsp = 1, SIZE(nIJ_llp_mmp_in,4)
!            if (jsp > 2 .and. .not.any(noco%l_spinoffd_ldau)) cycle
            WRITE (oUnit,'(a,i1,a,f12.6)') 'nIJ_llp_mmp distance spin ',jsp,' =',dist(jsp)
         ENDDO
      ENDIF

      IF (.FALSE.) THEN !(input%ldauLinMix) THEN

         ! mix here straight with given mixing factors
!         ALLOCATE (n_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(n_mmp_in,dim=3),SIZE(n_mmp_in,dim=4)))
!         n_mmp = cmplx_0

!         alpha = input%ldauMixParam
!         spinf = input%ldauSpinf

!         IF (input%jspins.EQ.1) THEN
!            DO i_u = 1, atoms%n_u
!               DO m = -lmaxU_const,lmaxU_const
!                  DO mp = -lmaxU_const,lmaxU_const

!                     n_mmp(m,mp,i_u,1) =      alpha * n_mmp_out(m,mp,i_u,1) + &
!                                        (1.0-alpha) * n_mmp_in (m,mp,i_u,1)

!                  END DO
!               END DO
!            END DO
!         ELSE
!            gam = 0.5 * alpha * (1.0 + spinf)
!            del = 0.5 * alpha * (1.0 - spinf)
!            DO i_u = 1,atoms%n_u
!               DO m = -lmaxU_const,lmaxU_const
!                  DO mp = -lmaxU_const,lmaxU_const

!                     n_mmp(m,mp,i_u,1) =       gam * n_mmp_out(m,mp,i_u,1) + &
!                                         (1.0-gam) * n_mmp_in (m,mp,i_u,1) - &
!                                               del * n_mmp_out(m,mp,i_u,2) + &
!                                               del * n_mmp_in (m,mp,i_u,2)

!                     n_mmp(m,mp,i_u,2) =       gam * n_mmp_out(m,mp,i_u,2) + &
!                                         (1.0-gam) * n_mmp_in (m,mp,i_u,2) - &
!                                               del * n_mmp_out(m,mp,i_u,1) + &
!                                               del * n_mmp_in (m,mp,i_u,1)
!                     IF(noco%l_mperp) THEN
!                        n_mmp(m,mp,i_u,3) =       alpha * n_mmp_out(m,mp,i_u,3) + &
!                                            (1.0-alpha) * n_mmp_in (m,mp,i_u,3)
!                     ENDIF

!                  END DO
!               END DO
!            END DO

!         ENDIF
!         n_mmp_in = n_mmp
!         DEALLOCATE(n_mmp)
      ENDIF

      CALL openXMLElementNoAttributes('ldaVDensityMatrixConvergence')
!      DO jsp = 1, SIZE(dist)
!         if (jsp > 2 .and. .not.any(noco%l_spinoffd_ldau)) cycle
!         attributes = ''
!         WRITE(attributes(1),'(i0)') jsp
!         WRITE(attributes(2),'(f13.6)') dist(jsp)
!         CALL writeXMLElementForm('distance',['spin    ','distance'],attributes(:2),reshape([4,8,1,13],[2,2]))
!      ENDDO
      CALL closeXMLElement('ldaVDensityMatrixConvergence')

   END SUBROUTINE v_mix
END MODULE m_vmix
