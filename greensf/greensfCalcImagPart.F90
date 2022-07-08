MODULE m_greensfCalcImagPart

   USE m_types
   USE m_constants
   USE m_tetrahedronInit
   USE m_juDFT

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfCalcImagPart_single_kpt(ikpt,ikpt_i,ev_list,jspin,gfinp,atoms,input,kpts,noco,mpi,&
                                             results,greensfBZintCoeffs,greensfImagPart)

      INTEGER,                   INTENT(IN)     :: ikpt, ikpt_i
      INTEGER,                   intent(IN)     :: ev_list(:)
      INTEGER,                   INTENT(IN)     :: jspin
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_kpts),              INTENT(IN)     :: kpts
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_mpi),               INTENT(IN)     :: mpi
      TYPE(t_results),           INTENT(IN)     :: results
      TYPE(t_greensfBZintCoeffs),INTENT(IN)     :: greensfBZintCoeffs
      TYPE(t_greensfImagPart),   INTENT(INOUT)  :: greensfImagPart

      INTEGER  :: nBands,jsp,i_gf,nLO,imatSize, spin_start, spin_end, spin_ind
      INTEGER  :: l,lp,m,mp,iBand,ie,j,eGrid_start,eGrid_end
      INTEGER  :: indUnique,i_elem,imat,iLO,iLOp,i_elemLO,i_elem_imag,i_elemLO_imag
      LOGICAL  :: l_zero,l_sphavg,l_kresolved_int,l_kresolved
      REAL     :: del,eb,wtkpt
      COMPLEX  :: fac,weight
      REAL,    ALLOCATABLE :: eig(:)
      REAL,    ALLOCATABLE :: eMesh(:)
      COMPLEX, ALLOCATABLE :: imag(:,:)
      REAL,    ALLOCATABLE :: weights(:,:)
      INTEGER, ALLOCATABLE :: indBound(:,:)

      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(results%ef,del,eb,eMesh=eMesh)

      spin_start = MERGE(1,jspin,gfinp%l_mperp)
      spin_end   = MERGE(3,jspin,gfinp%l_mperp)

      spin_start = MERGE(1           ,spin_start,noco%l_mperp.AND..NOT.gfinp%l_mperp)
      spin_end   = MERGE(input%jspins,spin_end  ,noco%l_mperp.AND..NOT.gfinp%l_mperp)

      do spin_ind=spin_start, spin_end
         !Spin degeneracy factors
         fac = 2.0/input%jspins

         nBands  = SIZE(ev_list)
         jsp     = MERGE(1,spin_ind,noco%l_noco)
         eig     = results%eig(ev_list,ikpt,jsp)

         SELECT CASE(input%bz_integration)
         CASE(BZINT_METHOD_HIST) !Histogram method
            wtkpt = kpts%wtkpt(ikpt)
         CASE(BZINT_METHOD_TETRA) !Tetrahedron method
            CALL timestart("Green's Function: TetrahedronWeights")
            ALLOCATE(weights(SIZE(eMesh),nBands),source=0.0)
            ALLOCATE(indBound(nBands,2),source=0)
            CALL tetrahedronInit(kpts,input,ikpt,results%eig(ev_list,:,jsp),nBands,eMesh,&
                                 weights,bounds=indBound,dos=.TRUE.)
            weights = - pi_const * fac * weights
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
            l_sphavg = gfinp%elem(i_gf)%l_sphavg
            l_kresolved_int = gfinp%elem(i_gf)%l_kresolved_int
            l_kresolved = gfinp%elem(i_gf)%l_kresolved


            IF(.NOT.gfinp%isUnique(i_gf, distinct_kresolved_int=.TRUE.)) CYCLE

            i_elem = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg)
            i_elemLO = gfinp%uniqueElements(atoms,max_index=i_gf,lo=.TRUE.,l_sphavg=l_sphavg)

            i_elem_imag = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg,l_kresolved_int=l_kresolved_int)
            i_elemLO_imag = gfinp%uniqueElements(atoms,max_index=i_gf,lo=.TRUE.,l_sphavg=l_sphavg,l_kresolved_int=l_kresolved_int)

            nLO = 0
            imatSize = 1
            IF(.NOT.l_sphavg) THEN
               imatSize = 4
               nLO = gfinp%elem(i_gf)%countLOs(atoms)
               IF(nLO/=0) THEN
                  imatSize = 4+4*nLO+nLO**2
               ENDIF
            ENDIF

            !$OMP parallel default(none) &
            !$OMP shared(input,gfinp,greensfBZintCoeffs,greensfImagPart) &
            !$OMP shared(i_elem,i_elemLO,i_elem_imag,i_elemLO_imag,nLO,l,lp,nBands,eMesh,l_sphavg,imatSize,l_kresolved_int,l_kresolved)&
            !$OMP shared(del,eb,eig,weights,indBound,fac,wtkpt,spin_ind,ikpt_i) &
            !$OMP private(ie,iLO,iLOp,imat,m,mp,iBand,j,eGrid_start,eGrid_end,weight,imag,l_zero)
            ALLOCATE(imag(SIZE(eMesh),imatSize),source=cmplx_0)
            !$OMP do collapse(2)
            DO mp = -lp, lp
               DO m = -l, l
                  imag = cmplx_0
                  DO iBand = 1, nBands

                     !Check for a non-zero weight in the energy interval
                     l_zero = .TRUE.
                     SELECT CASE(input%bz_integration)
                     CASE(BZINT_METHOD_HIST) !Histogram Method
                        j = FLOOR((eig(iBand)-eb)/del)+1
                        IF(j.LE.SIZE(eMesh).AND.j.GE.1) l_zero = .FALSE.
                        eGrid_start = j
                        eGrid_end   = j
                     CASE(BZINT_METHOD_TETRA) !Tetrahedron method
                        IF(ANY(ABS(weights(indBound(iBand,1):indBound(iBand,2),iBand)).GT.1e-14)) l_zero = .FALSE.
                        eGrid_start = indBound(iBand,1)
                        eGrid_end   = indBound(iBand,2)
                     CASE DEFAULT
                     END SELECT

                     IF(l_zero) CYCLE !No non-zero weight for this band

                     IF(eGrid_start==eGrid_end) THEN
                        ie = eGrid_start

                        SELECT CASE(input%bz_integration)
                        CASE(BZINT_METHOD_HIST) !Histogram Method
                           IF(l_kresolved) THEN
                              weight = -fac * pi_const/del
                           ELSE
                              weight = -fac * pi_const * wtkpt/del
                           ENDIF
                        CASE(BZINT_METHOD_TETRA) !Tetrahedron method
                           weight = weights(ie,iBand)
                        CASE DEFAULT
                        END SELECT

                        IF(l_sphavg) THEN
                           imag(ie,1) = imag(ie,1) + weight * greensfBZintCoeffs%sphavg(iBand,m,mp,i_elem,spin_ind)
                        ELSE
                           imag(ie,1) = imag(ie,1) + weight * greensfBZintCoeffs%uu(iBand,m,mp,i_elem,spin_ind)
                           imag(ie,2) = imag(ie,2) + weight * greensfBZintCoeffs%dd(iBand,m,mp,i_elem,spin_ind)
                           imag(ie,3) = imag(ie,3) + weight * greensfBZintCoeffs%ud(iBand,m,mp,i_elem,spin_ind)
                           imag(ie,4) = imag(ie,4) + weight * greensfBZintCoeffs%du(iBand,m,mp,i_elem,spin_ind)
                           IF(nLO>0) THEN
                              imat = 0
                              DO iLO = 1, nLO
                                 imat = imat + 4
                                 imag(ie,imat+1) = imag(ie,imat+1) + weight * greensfBZintCoeffs%uulo(iBand,m,mp,iLO,i_elemLO,spin_ind)
                                 imag(ie,imat+2) = imag(ie,imat+2) + weight * greensfBZintCoeffs%ulou(iBand,m,mp,iLO,i_elemLO,spin_ind)
                                 imag(ie,imat+3) = imag(ie,imat+3) + weight * greensfBZintCoeffs%dulo(iBand,m,mp,iLO,i_elemLO,spin_ind)
                                 imag(ie,imat+4) = imag(ie,imat+4) + weight * greensfBZintCoeffs%ulod(iBand,m,mp,iLO,i_elemLO,spin_ind)
                              ENDDO
                              imat = 0
                              DO iLO = 1, nLO
                                 DO iLOp = 1, nLO
                                    imat = imat + 1
                                    imag(ie,4 + 4*nLO+imat) = imag(ie,4 + 4*nLO+imat) + weight * greensfBZintCoeffs%uloulop(iBand,m,mp,imat,i_elemLO,spin_ind)
                                 ENDDO
                              ENDDO
                           ENDIF
                        ENDIF
                     ELSE
                        IF(l_sphavg) THEN
                           imag(eGrid_start:eGrid_end,1) = imag(eGrid_start:eGrid_end,1) + weights(eGrid_start:eGrid_end,iBand)&
                                                            * greensfBZintCoeffs%sphavg(iBand,m,mp,i_elem,spin_ind)
                        ELSE
                           imag(eGrid_start:eGrid_end,1) = imag(eGrid_start:eGrid_end,1) + weights(eGrid_start:eGrid_end,iBand)&
                                                            * greensfBZintCoeffs%uu(iBand,m,mp,i_elem,spin_ind)
                           imag(eGrid_start:eGrid_end,2) = imag(eGrid_start:eGrid_end,2) + weights(eGrid_start:eGrid_end,iBand)&
                                                            * greensfBZintCoeffs%dd(iBand,m,mp,i_elem,spin_ind)
                           imag(eGrid_start:eGrid_end,3) = imag(eGrid_start:eGrid_end,3) + weights(eGrid_start:eGrid_end,iBand)&
                                                            * greensfBZintCoeffs%ud(iBand,m,mp,i_elem,spin_ind)
                           imag(eGrid_start:eGrid_end,4) = imag(eGrid_start:eGrid_end,4) + weights(eGrid_start:eGrid_end,iBand)&
                                                            * greensfBZintCoeffs%du(iBand,m,mp,i_elem,spin_ind)
                           IF(nLO>0) THEN
                              imat = 0
                              DO iLO = 1, nLO
                                 imat = imat + 4
                                 imag(eGrid_start:eGrid_end,imat+1) = imag(eGrid_start:eGrid_end,imat+1) + weights(eGrid_start:eGrid_end,iBand) &
                                                                     * greensfBZintCoeffs%uulo(iBand,m,mp,iLO,i_elemLO,spin_ind)
                                 imag(eGrid_start:eGrid_end,imat+2) = imag(eGrid_start:eGrid_end,imat+2) + weights(eGrid_start:eGrid_end,iBand) &
                                                                     * greensfBZintCoeffs%ulou(iBand,m,mp,iLO,i_elemLO,spin_ind)
                                 imag(eGrid_start:eGrid_end,imat+3) = imag(eGrid_start:eGrid_end,imat+3) + weights(eGrid_start:eGrid_end,iBand) &
                                                                     * greensfBZintCoeffs%dulo(iBand,m,mp,iLO,i_elemLO,spin_ind)
                                 imag(eGrid_start:eGrid_end,imat+4) = imag(eGrid_start:eGrid_end,imat+4) + weights(eGrid_start:eGrid_end,iBand) &
                                                                     * greensfBZintCoeffs%ulod(iBand,m,mp,iLO,i_elemLO,spin_ind)
                              ENDDO
                              imat = 0
                              DO iLO = 1, nLO
                                 DO iLOp = 1, nLO
                                    imat = imat + 1
                                    imag(eGrid_start:eGrid_end,4 + 4*nLO+imat) = imag(eGrid_start:eGrid_end,4 + 4*nLO+imat) + weights(eGrid_start:eGrid_end,iBand) &
                                                                                 * greensfBZintCoeffs%uloulop(iBand,m,mp,imat,i_elemLO,spin_ind)
                                 ENDDO
                              ENDDO
                           ENDIF
                        ENDIF
                     ENDIF

                  ENDDO!ib
                  IF(l_sphavg.AND..NOT.l_kresolved_int) THEN
                     CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,1),1,greensfImagPart%sphavg(:,m,mp,i_elem_imag,spin_ind),1)
                  ELSE IF(.NOT.l_kresolved_int) THEN
                     CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,1),1,greensfImagPart%uu(:,m,mp,i_elem_imag,spin_ind),1)
                     CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,2),1,greensfImagPart%dd(:,m,mp,i_elem_imag,spin_ind),1)
                     CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,3),1,greensfImagPart%ud(:,m,mp,i_elem_imag,spin_ind),1)
                     CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,4),1,greensfImagPart%du(:,m,mp,i_elem_imag,spin_ind),1)

                     IF(nLO>0) THEN
                        imat = 0
                        DO iLO = 1, nLO
                           imat = imat + 4
                           CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,imat+1),1,greensfImagPart%uulo(:,m,mp,iLO,i_elemLO_imag,spin_ind),1)
                           CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,imat+2),1,greensfImagPart%ulou(:,m,mp,iLO,i_elemLO_imag,spin_ind),1)
                           CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,imat+3),1,greensfImagPart%dulo(:,m,mp,iLO,i_elemLO_imag,spin_ind),1)
                           CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,imat+4),1,greensfImagPart%ulod(:,m,mp,iLO,i_elemLO_imag,spin_ind),1)
                        ENDDO
                        imat = 0
                        DO iLO = 1, nLO
                           DO iLOp = 1, nLO
                              imat = imat + 1
                              CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,4 + 4*nLO+imat),1,greensfImagPart%uloulop(:,m,mp,iLO,iLOp,i_elemLO_imag,spin_ind),1)
                           ENDDO
                        ENDDO
                     ENDIF
                  ELSE
                     CALL zaxpy(SIZE(eMesh),cmplx_1,imag(:,1),1,greensfImagPart%sphavg_k(:,m,mp,i_elem_imag,spin_ind,ikpt_i),1)
                  ENDIF

               ENDDO!m
            ENDDO!mp
            !$OMP end do
            DEALLOCATE(imag)
            !$OMP end parallel
         ENDDO!i_gf
         CALL timestop("Green's Function: Imaginary Part")
      ENDDO !spin

      IF(input%bz_integration==BZINT_METHOD_TETRA) DEALLOCATE(weights,indBound)
   end subroutine

END MODULE m_greensfCalcImagPart
