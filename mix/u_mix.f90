!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_umix
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

   SUBROUTINE u_mix(input,atoms,noco,n_mmp_in,n_mmp_out)

      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_noco), INTENT(IN)    :: noco
      COMPLEX,      INTENT(IN)    :: n_mmp_out(-lmaxU_const:,-lmaxU_const:,:,:)
      COMPLEX,      INTENT(INOUT) :: n_mmp_in (-lmaxU_const:,-lmaxU_const:,:,:)


      INTEGER :: mp,m,l,itype,i_u,jsp
      REAL    :: alpha,spinf,gam,del,uParam,jParam
      REAL    :: zero(atoms%n_u),dist(SIZE(n_mmp_in,4))

      CHARACTER(LEN=20)   :: attributes(6)
      COMPLEX,ALLOCATABLE :: n_mmp(:,:,:,:)

      !
      ! check for possible rotation of n_mmp
      !
      !zero=0.0
      !CALL nmat_rot(zero,-atoms%lda_u%theta,-atoms%lda_u%phi,3,atoms%n_u,input%jspins,atoms%lda_u%l,n_mmp_out)

      ! Write out n_mmp_out to out.xml file
      CALL openXMLElementNoAttributes('ldaUDensityMatrix')
      DO jsp = 1, SIZE(n_mmp_out,4)
         DO i_u = 1, atoms%n_u
            l = atoms%lda_u(i_u)%l
            itype = atoms%lda_u(i_u)%atomType
            uParam = atoms%lda_u(i_u)%u
            jParam = atoms%lda_u(i_u)%j
            attributes = ''
            WRITE(attributes(1),'(i0)') jsp
            WRITE(attributes(2),'(i0)') itype
            WRITE(attributes(3),'(i0)') i_u
            WRITE(attributes(4),'(i0)') l
            WRITE(attributes(5),'(f15.8)') uParam
            WRITE(attributes(6),'(f15.8)') jParam
            CALL writeXMLElementMatrixPoly('densityMatrixFor',&
                                          (/'spin    ','atomType','uIndex  ','l       ','U       ','J       '/),&
                                          attributes,n_mmp_out(-l:l,-l:l,i_u,jsp))
         END DO
      END DO
      CALL closeXMLElement('ldaUDensityMatrix')

      ! exit subroutine if density matrix does not exist
      IF(.NOT.ANY(ABS(n_mmp_in(:,:,1:atoms%n_u,:)).GT.1e-12)) RETURN

      !Calculate distance
      dist = 0.0
      DO i_u = 1, atoms%n_u
         DO m = -lmaxU_const,lmaxU_const
            DO mp = -lmaxU_const,lmaxU_const
               DO jsp = 1, SIZE(n_mmp_in,4)
                  dist(jsp) = dist(jsp) + ABS(n_mmp_out(m,mp,i_u,jsp) - n_mmp_in(m,mp,i_u,jsp))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !Write to outfile
      IF(input%jspins.EQ.1) THEN
         WRITE (oUnit,'(a,f12.6)') 'n_mmp distance =',dist(1)
      ELSE
         DO jsp = 1, SIZE(n_mmp_in,4)
            if (jsp > 2 .and. .not.any(noco%l_spinoffd_ldau)) cycle
            WRITE (oUnit,9000) 'n_mmp distance spin ',jsp,' =',dist(jsp)
9000        FORMAT(a,I1,a,f12.6)
         ENDDO
      ENDIF

      IF (input%ldauLinMix) THEN

         ! mix here straight with given mixing factors
         ALLOCATE (n_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(n_mmp_in,dim=3),SIZE(n_mmp_in,dim=4)))
         n_mmp = cmplx_0

         alpha = input%ldauMixParam
         spinf = input%ldauSpinf

         IF (input%jspins.EQ.1) THEN
            DO i_u = 1, atoms%n_u
               DO m = -lmaxU_const,lmaxU_const
                  DO mp = -lmaxU_const,lmaxU_const

                     n_mmp(m,mp,i_u,1) =      alpha * n_mmp_out(m,mp,i_u,1) + &
                                        (1.0-alpha) * n_mmp_in (m,mp,i_u,1)

                  END DO
               END DO
            END DO
         ELSE
            gam = 0.5 * alpha * (1.0 + spinf)
            del = 0.5 * alpha * (1.0 - spinf)
            DO i_u = 1,atoms%n_u
               DO m = -lmaxU_const,lmaxU_const
                  DO mp = -lmaxU_const,lmaxU_const

                     n_mmp(m,mp,i_u,1) =       gam * n_mmp_out(m,mp,i_u,1) + &
                                         (1.0-gam) * n_mmp_in (m,mp,i_u,1) - &
                                               del * n_mmp_out(m,mp,i_u,2) + &
                                               del * n_mmp_in (m,mp,i_u,2)

                     n_mmp(m,mp,i_u,2) =       gam * n_mmp_out(m,mp,i_u,2) + &
                                         (1.0-gam) * n_mmp_in (m,mp,i_u,2) - &
                                               del * n_mmp_out(m,mp,i_u,1) + &
                                               del * n_mmp_in (m,mp,i_u,1)
                     IF(noco%l_mperp) THEN
                        n_mmp(m,mp,i_u,3) =       alpha * n_mmp_out(m,mp,i_u,3) + &
                                            (1.0-alpha) * n_mmp_in (m,mp,i_u,3)
                     ENDIF

                  END DO
               END DO
            END DO

         ENDIF
         n_mmp_in = n_mmp
         DEALLOCATE(n_mmp)
      ENDIF

      CALL openXMLElementNoAttributes('ldaUDensityMatrixConvergence')
      DO jsp = 1, SIZE(dist)
         attributes = ''
         WRITE(attributes(1),'(i0)') jsp
         WRITE(attributes(2),'(f13.6)') dist(jsp)
         CALL writeXMLElementForm('distance',['spin    ','distance'],attributes(:2),reshape([4,8,1,13],[2,2]))
      ENDDO
      CALL closeXMLElement('ldaUDensityMatrixConvergence')

   END SUBROUTINE u_mix
END MODULE m_umix
