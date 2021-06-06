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
   USE m_ldau_density
   use m_types_mixvector

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE u_mix(input,atoms,noco,l_mixldau,sphhar,vacuum,enpara,vTot,fmpi,hub1data,inden,outden, ldau_outden)

      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_noco),     INTENT(IN)    :: noco
      TYPE(t_sphhar),   INTENT(IN)    :: sphhar
      TYPE(t_vacuum),   INTENT(IN)    :: vacuum
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_potden),   INTENT(IN)    :: vTot
      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_hub1data), INTENT(IN)    :: hub1data
      LOGICAL,          INTENT(IN)    :: l_mixldau
      TYPE(t_potden),   INTENT(INOUT) :: inden
      TYPE(t_potden),   INTENT(INOUT) :: outden
      TYPE(t_potden),   INTENT(OUT)   :: ldau_outden


      INTEGER :: mp,m,l,itype,i_u,jsp,ispin,nType
      REAL    :: alpha,spinf,gam,del,uParam,jParam
      REAL    :: zero(atoms%n_u),dist(SIZE(inden%mmpmat,4))
      COMPLEX :: s, rho21

      CHARACTER(LEN=20)   :: attributes(6)
      COMPLEX,ALLOCATABLE :: n_mmp(:,:,:,:)

      !
      ! check for possible rotation of n_mmp
      !
      !zero=0.0
      !CALL nmat_rot(zero,-atoms%lda_u%theta,-atoms%lda_u%phi,3,atoms%n_u,input%jspins,atoms%lda_u%l,n_mmp_out)

      CALL density_from_mmpmat_coeffs(outden, atoms, input, noco, sphhar, enpara,&
                                      vTot, hub1data, fmpi, ldau_outDen, outden%mmpmat)

      IF(fmpi%irank==0) THEN
         ! Write out outden%mmpmat to out.xml file
         CALL openXMLElementNoAttributes('ldaUDensityMatrix')
         DO jsp = 1, SIZE(outden%mmpmat,4)
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
                                             attributes,outden%mmpmat(-l:l,-l:l,i_u,jsp))
            END DO
         END DO
         CALL closeXMLElement('ldaUDensityMatrix')
      ENDIF

      ! exit subroutine if density matrix does not exist
      IF(.NOT.ANY(ABS(inden%mmpmat(:,:,1:atoms%n_u,:)).GT.1e-12)) RETURN

      IF(fmpi%irank==0) THEN
         !Calculate distance
         dist = 0.0
         DO i_u = 1, atoms%n_u
            DO m = -lmaxU_const,lmaxU_const
               DO mp = -lmaxU_const,lmaxU_const
                  DO jsp = 1, SIZE(inden%mmpmat,4)
                     dist(jsp) = dist(jsp) + ABS(outden%mmpmat(m,mp,i_u,jsp) - inden%mmpmat(m,mp,i_u,jsp))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !Write to outfile
         IF(input%jspins.EQ.1) THEN
            WRITE (oUnit,'(a,f12.6)') 'n_mmp distance =',dist(1)
         ELSE
            DO jsp = 1, SIZE(inden%mmpmat,4)
               WRITE (oUnit,9000) 'n_mmp distance spin ',jsp,' =',dist(jsp)
9000           FORMAT(a,I1,a,f12.6)
            ENDDO
         ENDIF

         IF (input%ldauLinMix) THEN

            ! mix here straight with given mixing factors
            ALLOCATE (n_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(inDen%mmpmat,dim=3),SIZE(inden%mmpmat,dim=4)))
            n_mmp = cmplx_0

            alpha = input%ldauMixParam
            spinf = input%ldauSpinf

            IF (input%jspins.EQ.1) THEN
               DO i_u = 1, atoms%n_u
                  DO m = -lmaxU_const,lmaxU_const
                     DO mp = -lmaxU_const,lmaxU_const

                        n_mmp(m,mp,i_u,1) =      alpha * outden%mmpmat(m,mp,i_u,1) + &
                                           (1.0-alpha) * inden%mmpmat (m,mp,i_u,1)

                     END DO
                  END DO
               END DO
            ELSE
               gam = 0.5 * alpha * (1.0 + spinf)
               del = 0.5 * alpha * (1.0 - spinf)
               DO i_u = 1,atoms%n_u
                  DO m = -lmaxU_const,lmaxU_const
                     DO mp = -lmaxU_const,lmaxU_const

                        n_mmp(m,mp,i_u,1) =       gam * outden%mmpmat(m,mp,i_u,1) + &
                                            (1.0-gam) * inden%mmpmat (m,mp,i_u,1) - &
                                                  del * outden%mmpmat(m,mp,i_u,2) + &
                                                  del * inden%mmpmat (m,mp,i_u,2)

                        n_mmp(m,mp,i_u,2) =       gam * outden%mmpmat(m,mp,i_u,2) + &
                                            (1.0-gam) * inden%mmpmat (m,mp,i_u,2) - &
                                                  del * outden%mmpmat(m,mp,i_u,1) + &
                                                  del * inden%mmpmat (m,mp,i_u,1)
                        IF(noco%l_mperp) THEN
                           n_mmp(m,mp,i_u,3) =       alpha * outden%mmpmat(m,mp,i_u,3) + &
                                               (1.0-alpha) * inden%mmpmat (m,mp,i_u,3)
                        ENDIF

                     END DO
                  END DO
               END DO

            ENDIF
            inden%mmpmat = n_mmp
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
      ENDIF

      CALL inDen%distribute(fmpi%mpi_comm)

   END SUBROUTINE u_mix

   SUBROUTINE precond_ldau(imix,inDen, outDen, ldau_inden,ldau_outden, l_performed,sm, fsm)

      INTEGER,          INTENT(IN)     :: imix
      TYPE(t_potden),   INTENT(IN)     :: inden
      TYPE(t_potden),   INTENT(IN)     :: outden
      TYPE(t_potden),   INTENT(IN)     :: ldau_inden
      TYPE(t_potden),   INTENT(IN)     :: ldau_outden
      LOGICAL,          INTENT(OUT)    :: l_performed
      TYPE(t_mixvector),INTENT(INOUT)  :: sm
      TYPE(t_mixvector),INTENT(INOUT)  :: fsm

      TYPE(t_potden) :: inden_subbed, outden_subbed, delta_den

      l_performed =.FALSE.
      IF(imix==0) RETURN
      IF(.NOT.ALLOCATED(ldau_inden%mt)) RETURN
      l_performed =.TRUE.
      WRITE(*,*) 'precond'

      call inden_subbed%subPotDen(inDen,ldau_inden)
      call outden_subbed%subPotDen(outDen,ldau_outden)

      inden_subbed%mmpmat_uu = inden%mmpMAT_uu
      inden_subbed%mmpmat_ud = inden%mmpMAT_ud
      inden_subbed%mmpmat_du = inden%mmpMAT_du
      inden_subbed%mmpmat_dd = inden%mmpMAT_dd

      call sm%alloc()
      call sm%from_density(inden_subbed)

      call delta_den%subPotDen(outden_subbed,inden_subbed)

      delta_den%mmpmat_uu = outden%mmpmat_uu - inden%mmpMAT_uu
      delta_den%mmpmat_ud = outden%mmpmat_ud - inden%mmpMAT_ud
      delta_den%mmpmat_du = outden%mmpmat_du - inden%mmpMAT_du
      delta_den%mmpmat_dd = outden%mmpmat_dd - inden%mmpMAT_dd

      call fsm%alloc()
      call fsm%from_density(delta_den)

   END SUBROUTINE precond_ldau

END MODULE m_umix
