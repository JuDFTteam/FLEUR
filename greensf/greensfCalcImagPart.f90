MODULE m_greensfCalcImagPart

   USE m_types
   USE m_constants
   USE m_tetrahedronInit
   USE m_juDFT

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfCalcImagPart(cdnvalJob,spin_ind,gfinp,atoms,input,kpts,noco,mpi,&
                                  results,greensfBZintCoeffs,greensfImagPart)

      TYPE(t_cdnvalJob),         INTENT(IN)     :: cdnvalJob
      INTEGER,                   INTENT(IN)     :: spin_ind
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_kpts),              INTENT(IN)     :: kpts
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_mpi),               INTENT(IN)     :: mpi
      TYPE(t_results),           INTENT(IN)     :: results
      TYPE(t_greensfBZintCoeffs),INTENT(IN)     :: greensfBZintCoeffs
      TYPE(t_greensfImagPart),   INTENT(INOUT)  :: greensfImagPart


      INTEGER  :: ikpt_i,ikpt,nBands,jsp,i_gf
      INTEGER  :: l,lp,m,mp,iBand,ie,j,eGrid_start,eGrid_end
      INTEGER  :: indUnique,i_elem
      LOGICAL  :: l_zero
      REAL     :: del,eb,wtkpt
      COMPLEX  :: fac,weight
      INTEGER, ALLOCATABLE :: ev_list(:)
      REAL,    ALLOCATABLE :: eig(:)
      REAL,    ALLOCATABLE :: eMesh(:)
      REAL,    ALLOCATABLE :: dosWeights(:,:)
      INTEGER, ALLOCATABLE :: indBound(:,:)


      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(results%ef,del,eb,eMesh=eMesh)

      !Spin degeneracy and additional factors
      fac = -2.0/input%jspins * ImagUnit * pi_const

      DO ikpt_i = 1, SIZE(cdnvalJob%k_list)
         ikpt    = cdnvalJob%k_list(ikpt_i)
         ev_list = cdnvaljob%compact_ev_list(ikpt_i,.TRUE.)
         nBands  = SIZE(ev_list)
         jsp     = MERGE(1,spin_ind,noco%l_noco)
         eig     = results%eig(ev_list,ikpt,jsp)

         SELECT CASE(input%bz_integration)
         CASE(0) !Histogram method
            wtkpt = kpts%wtkpt(ikpt)
         CASE(3) !Tetrahedron method
            CALL timestart("Green's Function: TetrahedronWeights")
            ALLOCATE(dosWeights(gfinp%ne,nBands),source=0.0)
            ALLOCATE(indBound(nBands,2),source=0)
            CALL tetrahedronInit(kpts,ikpt,results%eig(ev_list,:,jsp),nBands,eMesh,gfinp%ne,&
                                 input%film,dosWeights,bounds=indBound,dos=.TRUE.)
            CALL timestop("Green's Function: TetrahedronWeights")
         CASE DEFAULT
            CALL juDFT_error("Invalid Brillouin-zone integration mode for Green's Functions",&
                              hint="Choose either hist or tetra",calledby="greensfCalcImagPart")
         END SELECT

         CALL timestart("Green's Function: Imaginary Part")
         !Loop over Green's Function elements
         !$OMP PARALLEL DEFAULT(NONE) &
         !$OMP SHARED(gfinp,input,greensfBZintCoeffs,greensfImagPart) &
         !$OMP SHARED(ikpt_i,ikpt,ev_list,nBands,del,eb,eig,dosWeights,indBound,fac,wtkpt,spin_ind) &
         !$OMP PRIVATE(i_gf,l,lp,m,mp,iBand,j,eGrid_start,eGrid_end,ie,weight,l_zero)&
         !$OMP PRIVATE(indUnique,i_elem)
         !$OMP DO
         DO i_gf = 1, gfinp%n

            !Get the information about the current element
            l  = gfinp%elem(i_gf)%l
            lp = gfinp%elem(i_gf)%lp

            CALL uniqueElements_gfinp(gfinp,i_elem,ind=i_gf,indUnique=indUnique)

            IF(i_gf/=indUnique) CYCLE

            DO m = -l, l
               DO mp = -lp, lp
                  DO iBand = 1, nBands

                     !Check for a non-zero weight in the energy interval
                     l_zero = .TRUE.
                     SELECT CASE(input%bz_integration)
                     CASE(0) !Histogram Method
                        j = FLOOR((eig(iBand)-eb)/del)+1
                        IF(j.LE.gfinp%ne.AND.j.GE.1) l_zero = .FALSE.
                        eGrid_start = j
                        eGrid_end   = j
                     CASE(3) !Tetrahedron method
                        IF(ANY(ABS(dosWeights(indBound(iBand,1):indBound(iBand,2),iBand)).GT.1e-14)) l_zero = .FALSE.
                        eGrid_start = indBound(iBand,1)
                        eGrid_end   = indBound(iBand,2)
                     CASE DEFAULT
                     END SELECT

                     IF(l_zero) CYCLE !No non-zero weight for this band

                     DO ie = eGrid_start, eGrid_end

                        SELECT CASE(input%bz_integration)
                        CASE(0) !Histogram Method
                           weight = fac * wtkpt/del
                        CASE(3) !Tetrahedron method
                           weight = fac * dosWeights(ie,iBand)
                        CASE DEFAULT
                        END SELECT

                        IF(gfinp%l_sphavg) THEN
                           greensfImagPart%sphavg(ie,m,mp,i_elem,spin_ind) = greensfImagPart%sphavg(ie,m,mp,i_elem,spin_ind) &
                                                                           + AIMAG(weight * greensfBZintCoeffs%sphavg(iBand,m,mp,ikpt_i,i_elem,spin_ind))
                        ELSE
                           greensfImagPart%uu(ie,m,mp,i_elem,spin_ind) = greensfImagPart%uu(ie,m,mp,i_elem,spin_ind) &
                                                                        + AIMAG(weight * greensfBZintCoeffs%uu(iBand,m,mp,ikpt_i,i_elem,spin_ind))
                           greensfImagPart%ud(ie,m,mp,i_elem,spin_ind) = greensfImagPart%uu(ie,m,mp,i_elem,spin_ind) &
                                                                        + AIMAG(weight * greensfBZintCoeffs%ud(iBand,m,mp,ikpt_i,i_elem,spin_ind))
                           greensfImagPart%du(ie,m,mp,i_elem,spin_ind) = greensfImagPart%uu(ie,m,mp,i_elem,spin_ind) &
                                                                        + AIMAG(weight * greensfBZintCoeffs%du(iBand,m,mp,ikpt_i,i_elem,spin_ind))
                           greensfImagPart%dd(ie,m,mp,i_elem,spin_ind) = greensfImagPart%uu(ie,m,mp,i_elem,spin_ind) &
                                                                        + AIMAG(weight * greensfBZintCoeffs%dd(iBand,m,mp,ikpt_i,i_elem,spin_ind))
                        ENDIF

                     ENDDO!ie
                  ENDDO!ib
               ENDDO!mp
            ENDDO!m
         ENDDO!i_gf
         !$OMP END DO
         !$OMP END PARALLEL
         CALL timestop("Green's Function: Imaginary Part")

         IF(input%bz_integration==3) DEALLOCATE(dosWeights,indBound)
      ENDDO!k-point loop

      !Collect the results from all mpi ranks
      CALL greensfImagPart%collect(spin_ind,mpi%mpi_comm)

   END SUBROUTINE greensfCalcImagPart
END MODULE m_greensfCalcImagPart