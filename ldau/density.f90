MODULE m_ldau_density

   USE m_types
   USE m_juDFT
   USE m_genMTBasis

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE density_from_mmpmat_coeffs(den, atoms, input, noco, sphhar, enpara, vTot, hub1data, fmpi, ldauDen, mmpmat)

      TYPE(t_potden),   INTENT(IN) :: den
      TYPE(t_atoms),    INTENT(IN) :: atoms
      TYPE(t_input),    INTENT(IN) :: input
      TYPE(t_noco),     INTENT(IN) :: noco
      TYPE(t_sphhar),   INTENT(IN) :: sphhar
      TYPE(t_enpara),   INTENT(IN) :: enpara
      TYPE(t_potden),   INTENT(IN) :: vTot
      TYPE(t_hub1data), INTENT(IN) :: hub1data
      TYPE(t_mpi),      INTENT(IN) :: fmpi
      TYPE(t_potden),   INTENT(OUT):: ldauDen
      COMPLEX, OPTIONAL, ALLOCATABLE, INTENT(INOUT) :: mmpmat(:,:,:,:)


      INTEGER :: nType, l, m, ispin, i_u, iGrid
      COMPLEX :: s, rho21

      TYPE(t_usdus) :: usdus
      TYPE(t_denCoeffsOffdiag) :: denCoeffsOffdiag
      REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:),flo(:,:,:,:) ! radial functions

      WRITE(*,*) ALLOCATED(den%mmpmat_uu), SHAPE(den%mmpmat_uu)
      WRITE(*,*) ALLOCATED(den%mmpmat_ud), SHAPE(den%mmpmat_ud)
      WRITE(*,*) ALLOCATED(den%mmpmat_du), SHAPE(den%mmpmat_du)
      WRITE(*,*) ALLOCATED(den%mmpmat_dd), SHAPE(den%mmpmat_dd)
      WRITE(*,*) ALLOCATED(den%mmpmat), SHAPE(den%mmpmat)
      call ldauDen%copyPotDen(den)
      call ldauDen%resetPotDen()

      IF(PRESENT(mmpMat)) mmpMat = ldauDen%mmpMat

      CALL usdus%init(atoms,input%jspins)
      CALL denCoeffsOffdiag%init(atoms,noco,sphhar,.FALSE.,.FALSE.)

      ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)) ! Deallocation before mpi_col_den
      ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins))
      ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins))
      !Calculate the mmpmat from the coefficients and the corresponding
      !density if needed
      DO i_u = 1, atoms%n_u
         nType = atoms%lda_u(i_u)%atomType
         l = atoms%lda_u(i_u)%l
         DO ispin = 1, input%jspins
            CALL genMTBasis(atoms,enpara,vTot,fmpi,nType,ispin,usdus,f(:,:,0:,ispin),g(:,:,0:,ispin),flo(:,:,:,ispin),&
                            hub1data=hub1data,l_writeArg=.FALSE.)

            DO m = -l, l
               DO iGrid = 1,atoms%jri(nType)
                  s =  den%mmpmat_uu(m,m,i_u,ispin)*( f(iGrid,1,l,ispin)*f(iGrid,1,l,ispin)+f(iGrid,2,l,ispin)*f(iGrid,2,l,ispin) )&
                     + den%mmpmat_dd(m,m,i_u,ispin)*( g(iGrid,1,l,ispin)*g(iGrid,1,l,ispin)+g(iGrid,2,l,ispin)*g(iGrid,2,l,ispin) )&
                     + den%mmpmat_du(m,m,i_u,ispin)*( f(iGrid,1,l,ispin)*g(iGrid,1,l,ispin)+f(iGrid,2,l,ispin)*g(iGrid,2,l,ispin) )&
                     + den%mmpmat_ud(m,m,i_u,ispin)*( g(iGrid,1,l,ispin)*f(iGrid,1,l,ispin)+g(iGrid,2,l,ispin)*f(iGrid,2,l,ispin) )
                  ldauDen%mt(iGrid,0,nType,ispin) = ldauDen%mt(iGrid,0,nType,ispin) + REAL(s)/sfp_const
               ENDDO
            ENDDO

            IF(PRESENT(mmpMat)) THEN
               mmpmat(:,:,i_u,ispin) =  den%mmpmat_uu(:,:,i_u,ispin) &
                                      + den%mmpmat_dd(:,:,i_u,ispin) * usdus%ddn(l,nType,ispin)
            ENDIF
         END DO
         IF(noco%l_mperp) THEN

            DO m = -l, l
               DO iGrid = 1,atoms%jri(nType)
                  s =  den%mmpmat_uu(m,m,i_u,3)*( f(iGrid,1,l,2)*f(iGrid,1,l,1)+f(iGrid,2,l,2)*f(iGrid,2,l,1) )&
                     + den%mmpmat_dd(m,m,i_u,3)*( g(iGrid,1,l,2)*g(iGrid,1,l,1)+g(iGrid,2,l,2)*g(iGrid,2,l,1) )&
                     + den%mmpmat_du(m,m,i_u,3)*( f(iGrid,1,l,2)*g(iGrid,1,l,1)+f(iGrid,2,l,2)*g(iGrid,2,l,1) )&
                     + den%mmpmat_ud(m,m,i_u,3)*( g(iGrid,1,l,2)*f(iGrid,1,l,1)+g(iGrid,2,l,2)*f(iGrid,2,l,1) )
                  rho21=CONJG(s)/sfp_const

                  ldauDen%mt(iGrid,0,nType,3) = ldauDen%mt(iGrid,0,nType,3) - REAL(rho21)
                  ldauDen%mt(iGrid,0,nType,4) = ldauDen%mt(iGrid,0,nType,4) + AIMAG(rho21)
               ENDDO
            ENDDO

            IF(PRESENT(mmpMat)) THEN

               CALL denCoeffsOffdiag%addRadFunScalarProducts(atoms,f,g,flo,nType)
               mmpmat(:,:,i_u,3) =  den%mmpmat_uu(:,:,i_u,3) * denCoeffsOffdiag%uu21n(l,nType) &
                                  + den%mmpmat_ud(:,:,i_u,3) * denCoeffsOffdiag%ud21n(l,nType) &
                                  + den%mmpmat_du(:,:,i_u,3) * denCoeffsOffdiag%du21n(l,nType) &
                                  + den%mmpmat_dd(:,:,i_u,3) * denCoeffsOffdiag%dd21n(l,nType)
            ENDIF
         ENDIF

      ENDDO

   END SUBROUTINE density_from_mmpmat_coeffs

END MODULE m_ldau_density