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

#include"cpp_double.h"

      INTEGER  :: ikpt_i,ikpt,nBands,jsp,i_gf
      INTEGER  :: l,lp,m,mp,iBand,ie,j,eGrid_start,eGrid_end
      INTEGER  :: indUnique,i_elem
      LOGICAL  :: l_zero
      REAL     :: del,eb,wtkpt
      COMPLEX  :: fac,weight
      INTEGER, ALLOCATABLE :: ev_list(:)
      REAL,    ALLOCATABLE :: eig(:)
      REAL,    ALLOCATABLE :: eMesh(:)
      REAL,    ALLOCATABLE :: imagReal(:,:)
      COMPLEX, ALLOCATABLE :: weights(:,:),imag(:,:)
      REAL,    ALLOCATABLE :: dosWeights(:,:)
      REAL,    ALLOCATABLE :: resWeights(:,:)
      INTEGER, ALLOCATABLE :: indBound(:,:)


      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(results%ef,del,eb,eMesh=eMesh)

      !Spin degeneracy factors
      fac = 2.0/input%jspins

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
            ALLOCATE(dosWeights(SIZE(eMesh),nBands),source=0.0)
            ALLOCATE(resWeights(SIZE(eMesh),nBands),source=0.0)
            ALLOCATE(indBound(nBands,2),source=0)
            IF(.FALSE..AND..NOT.input%film) THEN !Here we also calculate the resolvent weights
               CALL tetrahedronInit(kpts,ikpt,results%eig(ev_list,:,jsp),nBands,eMesh,SIZE(eMesh),&
                                    input%film,dosWeights,bounds=indBound,resWeights=resWeights,dos=.TRUE.)
            ELSE
               CALL tetrahedronInit(kpts,ikpt,results%eig(ev_list,:,jsp),nBands,eMesh,SIZE(eMesh),&
                                    input%film,dosWeights,bounds=indBound,dos=.TRUE.)
            ENDIF
            ALLOCATE(weights(SIZE(eMesh),nBands),source=cmplx_0)
            weights = fac * (resWeights - ImagUnit * pi_const * dosWeights)
            DEALLOCATE(dosWeights,resWeights)
            CALL timestop("Green's Function: TetrahedronWeights")
         CASE DEFAULT
            CALL juDFT_error("Invalid Brillouin-zone integration mode for Green's Functions",&
                              hint="Choose either hist or tetra",calledby="greensfCalcImagPart")
         END SELECT

         CALL timestart("Green's Function: Imaginary Part")
         !Loop over Green's Function elements
         DO i_gf = 1, gfinp%n

            !Get the information about the current element
            l  = gfinp%elem(i_gf)%l
            lp = gfinp%elem(i_gf)%lp

            i_elem = gfinp%uniqueElements(ind=i_gf,indUnique=indUnique)

            IF(i_gf/=indUnique) CYCLE

            !$OMP parallel default(none) &
            !$OMP shared(gfinp,input,greensfBZintCoeffs,greensfImagPart) &
            !$OMP shared(i_elem,l,lp,ikpt_i,nBands,eMesh)&
            !$OMP shared(del,eb,eig,weights,indBound,fac,wtkpt,spin_ind) &
            !$OMP private(ie,m,mp,iBand,j,eGrid_start,eGrid_end,weight,imag,imagReal,l_zero)
            ALLOCATE(imag(SIZE(eMesh),MERGE(1,4,gfinp%l_sphavg)),source=cmplx_0)
            ALLOCATE(imagReal(SIZE(eMesh),MERGE(1,4,gfinp%l_sphavg)),source=0.0)
            !$OMP do collapse(2)
            DO mp = -lp, lp
               DO m = -l, l
                  imag = cmplx_0
                  DO iBand = 1, nBands

                     !Check for a non-zero weight in the energy interval
                     l_zero = .TRUE.
                     SELECT CASE(input%bz_integration)
                     CASE(0) !Histogram Method
                        j = FLOOR((eig(iBand)-eb)/del)+1
                        IF(j.LE.SIZE(eMesh).AND.j.GE.1) l_zero = .FALSE.
                        eGrid_start = j
                        eGrid_end   = j
                     CASE(3) !Tetrahedron method
                        IF(ANY(ABS(weights(indBound(iBand,1):indBound(iBand,2),iBand)).GT.1e-14)) l_zero = .FALSE.
                        eGrid_start = indBound(iBand,1)
                        eGrid_end   = indBound(iBand,2)
                     CASE DEFAULT
                     END SELECT

                     IF(l_zero) CYCLE !No non-zero weight for this band

                     IF(eGrid_start==eGrid_end) THEN
                        ie = eGrid_start

                        SELECT CASE(input%bz_integration)
                        CASE(0) !Histogram Method
                           weight = -fac * ImagUnit * pi_const * wtkpt/del
                        CASE(3) !Tetrahedron method
                           weight = weights(ie,iBand)
                        CASE DEFAULT
                        END SELECT

                        IF(gfinp%l_sphavg) THEN
                           imag(ie,1) = imag(ie,1) + weight * greensfBZintCoeffs%sphavg(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                        ELSE
                           imag(ie,1) = imag(ie,1) + weight * greensfBZintCoeffs%uu(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                           imag(ie,2) = imag(ie,2) + weight * greensfBZintCoeffs%dd(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                           imag(ie,3) = imag(ie,3) + weight * greensfBZintCoeffs%ud(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                           imag(ie,4) = imag(ie,4) + weight * greensfBZintCoeffs%du(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                        ENDIF

                     ELSE IF(eGrid_start==1 .AND. eGrid_end==SIZE(eMesh)) THEN!Here we always use the tetrahedron method
                        !We can only use the BLAS routine on the full array
                        IF(gfinp%l_sphavg) THEN
                           CALL CPP_BLAS_caxpy(SIZE(eMesh),greensfBZintCoeffs%sphavg(iBand,m,mp,ikpt_i,i_elem,spin_ind),&
                                               weights(:,iBand),1,imag(:,1),1)
                        ELSE
                           CALL CPP_BLAS_caxpy(SIZE(eMesh),greensfBZintCoeffs%uu(iBand,m,mp,ikpt_i,i_elem,spin_ind),&
                                               weights(:,iBand),1,imag(:,1),1)
                           CALL CPP_BLAS_caxpy(SIZE(eMesh),greensfBZintCoeffs%dd(iBand,m,mp,ikpt_i,i_elem,spin_ind),&
                                               weights(:,iBand),1,imag(:,2),1)
                           CALL CPP_BLAS_caxpy(SIZE(eMesh),greensfBZintCoeffs%ud(iBand,m,mp,ikpt_i,i_elem,spin_ind),&
                                               weights(:,iBand),1,imag(:,3),1)
                           CALL CPP_BLAS_caxpy(SIZE(eMesh),greensfBZintCoeffs%du(iBand,m,mp,ikpt_i,i_elem,spin_ind),&
                                               weights(:,iBand),1,imag(:,4),1)
                        ENDIF
                     ELSE
                        IF(gfinp%l_sphavg) THEN
                           imag(eGrid_start:eGrid_end,1) = imag(eGrid_start:eGrid_end,1) + weights(eGrid_start:eGrid_end,iBand)&
                                                          * greensfBZintCoeffs%sphavg(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                        ELSE
                           imag(eGrid_start:eGrid_end,1) = imag(eGrid_start:eGrid_end,1) + weights(eGrid_start:eGrid_end,iBand)&
                                                          * greensfBZintCoeffs%uu(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                           imag(eGrid_start:eGrid_end,2) = imag(eGrid_start:eGrid_end,2) + weights(eGrid_start:eGrid_end,iBand)&
                                                          * greensfBZintCoeffs%dd(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                           imag(eGrid_start:eGrid_end,3) = imag(eGrid_start:eGrid_end,3) + weights(eGrid_start:eGrid_end,iBand)&
                                                          * greensfBZintCoeffs%ud(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                           imag(eGrid_start:eGrid_end,4) = imag(eGrid_start:eGrid_end,4) + weights(eGrid_start:eGrid_end,iBand)&
                                                          * greensfBZintCoeffs%du(iBand,m,mp,ikpt_i,i_elem,spin_ind)
                        ENDIF
                     ENDIF

                  ENDDO!ib
                  imagReal = AIMAG(imag)
                  IF(gfinp%l_sphavg) THEN
                     CALL CPP_BLAS_saxpy(SIZE(eMesh),1.0,imagReal(:,1),1,greensfImagPart%sphavg(:,m,mp,i_elem,spin_ind),1)
                  ELSE
                     CALL CPP_BLAS_saxpy(SIZE(eMesh),1.0,imagReal(:,1),1,greensfImagPart%uu(:,m,mp,i_elem,spin_ind),1)
                     CALL CPP_BLAS_saxpy(SIZE(eMesh),1.0,imagReal(:,2),1,greensfImagPart%dd(:,m,mp,i_elem,spin_ind),1)
                     CALL CPP_BLAS_saxpy(SIZE(eMesh),1.0,imagReal(:,3),1,greensfImagPart%ud(:,m,mp,i_elem,spin_ind),1)
                     CALL CPP_BLAS_saxpy(SIZE(eMesh),1.0,imagReal(:,4),1,greensfImagPart%du(:,m,mp,i_elem,spin_ind),1)
                  ENDIF

               ENDDO!m
            ENDDO!mp
            !$OMP end do
            !$OMP end parallel
         ENDDO!i_gf
         CALL timestop("Green's Function: Imaginary Part")

         IF(input%bz_integration==3) DEALLOCATE(weights,indBound)
      ENDDO!k-point loop

      !Collect the results from all mpi ranks
      CALL greensfImagPart%collect(spin_ind,mpi%mpi_comm)

   END SUBROUTINE greensfCalcImagPart
END MODULE m_greensfCalcImagPart