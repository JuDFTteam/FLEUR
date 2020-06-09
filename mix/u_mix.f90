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
   USE m_nmat_rot
   USE m_xmlOutput

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE u_mix(input,atoms,noco,n_mmp_in,n_mmp_out)

      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_noco), INTENT(IN)    :: noco
      COMPLEX,      INTENT(INOUT) :: n_mmp_out(-lmaxU_const:,-lmaxU_const:,:,:)
      COMPLEX,      INTENT(INOUT) :: n_mmp_in (-lmaxU_const:,-lmaxU_const:,:,:)


      INTEGER :: j,k,l,itype,i_u,jsp
      REAL    :: alpha,spinf,gam,del,sum1,sum2,sum3,uParam,jParam
      REAL    :: zero(atoms%n_u)

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

      IF (input%ldauLinMix) THEN

         ! mix here straight with given mixing factors
         ALLOCATE (n_mmp,mold=n_mmp_in)
         n_mmp = cmplx_0

         alpha = input%ldauMixParam
         spinf = input%ldauSpinf

         sum1 = 0.0
         IF (input%jspins.EQ.1) THEN
            DO i_u = 1, atoms%n_u
               DO j = -3,3
                  DO k = -3,3
                     sum1 = sum1 + ABS(n_mmp_out(k,j,i_u,1) - n_mmp_in(k,j,i_u,1))
                     n_mmp(k,j,i_u,1) = alpha * n_mmp_out(k,j,i_u,1) + (1.0-alpha) * n_mmp_in(k,j,i_u,1)
                  END DO
               END DO
            END DO
            WRITE (oUnit,'(a16,f12.6)') 'n_mmp distance =',sum1
         ELSE
            sum2 = 0.0
            sum3 = 0.0
            gam = 0.5 * alpha * (1.0 + spinf)
            del = 0.5 * alpha * (1.0 - spinf)
            DO i_u = 1,atoms%n_u
               DO j = -3,3
                  DO k = -3,3
                     sum1 = sum1 + ABS(n_mmp_out(k,j,i_u,1) - n_mmp_in(k,j,i_u,1))
                     sum2 = sum2 + ABS(n_mmp_out(k,j,i_u,2) - n_mmp_in(k,j,i_u,2))
                     IF(noco%l_mperp) sum3 = sum3 + ABS(n_mmp_out(k,j,i_u,3) - n_mmp_in(k,j,i_u,3))

                     n_mmp(k,j,i_u,1) =       gam * n_mmp_out(k,j,i_u,1) + &
                                        (1.0-gam) * n_mmp_in (k,j,i_u,1) + &
                                              del * n_mmp_out(k,j,i_u,2) - &
                                              del * n_mmp_in (k,j,i_u,2)

                     n_mmp(k,j,i_u,2) =       gam * n_mmp_out(k,j,i_u,2) + &
                                        (1.0-gam) * n_mmp_in (k,j,i_u,2) + &
                                              del * n_mmp_out(k,j,i_u,1) - &
                                              del * n_mmp_in (k,j,i_u,1)
                     IF(noco%l_mperp) THEN
                        n_mmp(k,j,i_u,3) =       alpha * n_mmp_out(k,j,i_u,3) + &
                                           (1.0-alpha) * n_mmp_in (k,j,i_u,3)
                     ENDIF
                  END DO
               END DO
            END DO
            WRITE (oUnit,'(a23,f12.6)') 'n_mmp distance spin 1 =',sum1
            WRITE (oUnit,'(a23,f12.6)') 'n_mmp distance spin 2 =',sum2
            IF(noco%l_mperp) WRITE (oUnit,'(a23,f12.6)') 'n_mmp distance spin 3 =',sum3
         ENDIF
         n_mmp_in = n_mmp
         DEALLOCATE(n_mmp)
      ELSE ! input%ldauLinMix

         ! only calculate distance

         sum1 = 0.0
         DO i_u = 1, atoms%n_u
            DO j = -3,3
               DO k = -3,3
                  sum1 = sum1 + ABS(n_mmp_out(k,j,i_u,1) - n_mmp_in(k,j,i_u,1))
               END DO
            END DO
         END DO
         IF (input%jspins.EQ.1) THEN
            WRITE (oUnit,'(a16,f12.6)') 'n_mmp distance =',sum1
         ELSE
            sum2 = 0.0
            WRITE (oUnit,'(a23,f12.6)') 'n_mmp distance spin 1 =',sum1
            DO i_u = 1, atoms%n_u
               DO j = -3,3
                  DO k = -3,3
                     sum2 = sum2 + ABS(n_mmp_out(k,j,i_u,2) - n_mmp_in(k,j,i_u,2))
                  END DO
               END DO
            END DO
            WRITE (oUnit,'(a23,f12.6)') 'n_mmp distance spin 2 =',sum2
            IF(noco%l_mperp) THEN
               !Spin off-diagonal
               sum3 = 0.0
               DO i_u = 1, atoms%n_u
                  DO j = -3,3
                     DO k = -3,3
                        sum3 = sum3 + ABS(n_mmp_out(k,j,i_u,3) - n_mmp_in(k,j,i_u,3))
                     END DO
                  END DO
               END DO
               WRITE (oUnit,'(a23,f12.6)') 'n_mmp distance spin 3 =',sum3
            ENDIF
         END IF
      END IF ! input%ldauLinMix

   END SUBROUTINE u_mix
END MODULE m_umix
