!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_greensfCoeffs

   USE m_juDFT
   USE m_types_setup
   USE m_constants

   IMPLICIT NONE

   PRIVATE

      TYPE t_greensfBZintCoeffs

         !Contains only the coefficients for each kpt and band handled by the current mpi rank

         COMPLEX, ALLOCATABLE :: sphavg(:,:,:,:,:)

         ! These arrays are only used in the case we want the green's function with radial dependence
         COMPLEX, ALLOCATABLE :: uu(:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dd(:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: du(:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ud(:,:,:,:,:)

         !LO-Valence Contribution
         COMPLEX, ALLOCATABLE :: uulo(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ulou(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dulo(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ulod(:,:,:,:,:,:)

         !LO-LO contribution
         !Here we need to compress the (lo,lop) index pair into one index because PGI allows a max of 7 indices
         COMPLEX, ALLOCATABLE :: uloulop(:,:,:,:,:,:)

         CONTAINS
            PROCEDURE, PASS :: init             => greensfBZintCoeffs_init
            PROCEDURE, PASS :: add_contribution => greensfBZintCoeffs_add_contribution
            PROCEDURE, PASS :: reset            => greensfBZintCoeffs_reset
      END TYPE t_greensfBZintCoeffs


      TYPE t_greensfImagPart

         !Contains the imaginary part of the greens function
         INTEGER, ALLOCATABLE :: kkintgr_cutoff(:,:,:)
         REAL   , ALLOCATABLE :: scalingFactorSphavg(:,:)
         REAL   , ALLOCATABLE :: scalingFactorRadial(:,:)
         REAL   , ALLOCATABLE :: scalingFactorSphavgKres(:,:)
         LOGICAL :: l_calc = .FALSE.

         COMPLEX, ALLOCATABLE :: sphavg(:,:,:,:,:)

         ! These arrays are only used in the case we want the green's function with radial dependence
         COMPLEX, ALLOCATABLE :: uu(:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dd(:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: du(:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ud(:,:,:,:,:)

         !LO-Valence Contribution
         COMPLEX, ALLOCATABLE :: uulo(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ulou(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dulo(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ulod(:,:,:,:,:,:)

         !LO-LO contribution
         !Here the (lo,lop) index pair is explicit again
         COMPLEX, ALLOCATABLE :: uloulop(:,:,:,:,:,:,:)


         !K-resolved greens functions
         !(Radially resolved not yet implemented since we hit the indices limit on PGI)
         COMPLEX, ALLOCATABLE :: sphavg_k(:,:,:,:,:,:)

         CONTAINS
            PROCEDURE, PASS :: init        =>  greensfImagPart_init
            PROCEDURE, PASS :: collect     =>  greensfImagPart_collect
            PROCEDURE, PASS :: mpi_bc      =>  greensfImagPart_mpi_bc
            PROCEDURE       :: scale       =>  greensfImagPart_scale
            PROCEDURE       :: applyCutoff =>  greensfImagPart_applyCutoff
            PROCEDURE       :: checkEmpty  =>  greensfImagPart_checkEmpty
      END TYPE t_greensfImagPart

   PUBLIC t_greensfBZintCoeffs, t_greensfImagPart

   CONTAINS

      SUBROUTINE greensfBZintCoeffs_init(this,gfinp,atoms,noco,jsp_start,jsp_end,nbands)

         CLASS(t_greensfBZintCoeffs),  INTENT(INOUT)  :: this
         TYPE(t_gfinp),                INTENT(IN)     :: gfinp
         TYPE(t_atoms),                INTENT(IN)     :: atoms
         TYPE(t_noco),                 INTENT(IN)     :: noco
         INTEGER,                      INTENT(IN)     :: jsp_start,jsp_end
         INTEGER,                      INTENT(IN)     :: nbands !number of kpts and bands handled by this rank

         INTEGER lmax, uniqueElementsSphavg,uniqueElementsRadial, maxSpin,uniqueElementsLO,maxLO

         lmax = lmaxU_const

         IF(gfinp%l_mperp.AND.jsp_end==2) THEN
            maxSpin = 3
         ELSE
            maxSpin = jsp_end
         ENDIF

         !Determine number of unique gf elements
         uniqueElementsSphavg  = gfinp%uniqueElements(atoms,l_sphavg=.TRUE.) !How many spherically averaged elements
         uniqueElementsRadial  = gfinp%uniqueElements(atoms,l_sphavg=.FALSE.) !How many elements with radial dependence

         IF(uniqueElementsSphavg>0) THEN
            ALLOCATE (this%sphavg(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsSphavg,jsp_start:maxSpin),source=cmplx_0)
         ENDIF
         IF(uniqueElementsRadial>0) THEN
            ALLOCATE (this%uu(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,jsp_start:maxSpin),source=cmplx_0)
            ALLOCATE (this%dd(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,jsp_start:maxSpin),source=cmplx_0)
            ALLOCATE (this%du(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,jsp_start:maxSpin),source=cmplx_0)
            ALLOCATE (this%ud(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,jsp_start:maxSpin),source=cmplx_0)

            uniqueElementsLO = gfinp%uniqueElements(atoms,lo=.TRUE.,l_sphavg=.FALSE.,maxLO=maxLO)

            IF(uniqueElementsLO>0) THEN
               ALLOCATE (this%uulo(nbands,-lmax:lmax,-lmax:lmax,maxLO,uniqueElementsLO,jsp_start:maxSpin),source=cmplx_0)
               ALLOCATE (this%ulou(nbands,-lmax:lmax,-lmax:lmax,maxLO,uniqueElementsLO,jsp_start:maxSpin),source=cmplx_0)
               ALLOCATE (this%dulo(nbands,-lmax:lmax,-lmax:lmax,maxLO,uniqueElementsLO,jsp_start:maxSpin),source=cmplx_0)
               ALLOCATE (this%ulod(nbands,-lmax:lmax,-lmax:lmax,maxLO,uniqueElementsLO,jsp_start:maxSpin),source=cmplx_0)

               ALLOCATE (this%uloulop(nbands,-lmax:lmax,-lmax:lmax,maxLO**2,uniqueElementsLO,jsp_start:maxSpin),source=cmplx_0)
            ENDIF
         ENDIF

      END SUBROUTINE greensfBZintCoeffs_init

      subroutine greensfBZintCoeffs_reset(this)

         class(t_greensfBZintCoeffs), intent(inout) :: this

         if (allocated(this%sphavg)) this%sphavg = cmplx_0
         if (allocated(this%uu)) then
            this%uu = cmplx_0
            this%ud = cmplx_0
            this%du = cmplx_0
            this%dd = cmplx_0
         endif
         if (allocated(this%uulo)) then
            this%uulo = cmplx_0
            this%ulou = cmplx_0
            this%dulo = cmplx_0
            this%ulod = cmplx_0
            this%uloulop = cmplx_0
         endif

      end subroutine

      SUBROUTINE greensfBZintCoeffs_add_contribution(this, i_elem, i_elemLO, iBand, ispin, nLO, imat, l_sphavg, contribution)

         CLASS(t_greensfBZintCoeffs),  INTENT(INOUT)   :: this
         INTEGER,                      INTENT(IN)      :: i_elem,i_elemLO,ispin,nLO,iBand,imat
         LOGICAL,                      INTENT(IN)      :: l_sphavg
         COMPLEX,                      INTENT(IN)      :: contribution(-lmaxU_const:,-lmaxU_const:)

         INTEGER :: iLO

         IF(l_sphavg) THEN
            !Spherically averaged (already multiplied with scalar products)
            this%sphavg(iBand,:,:,i_elem,ispin) = &
               this%sphavg(iBand,:,:,i_elem,ispin) + contribution
         ELSE IF(imat.EQ.1) THEN
            !imat 1-4: coefficients for Valence-Valence contribution
            this%uu(iBand,:,:,i_elem,ispin) = &
               this%uu(iBand,:,:,i_elem,ispin) + contribution
         ELSE IF(imat.EQ.2) THEN
            this%dd(iBand,:,:,i_elem,ispin) = &
               this%dd(iBand,:,:,i_elem,ispin) + contribution
         ELSE IF(imat.EQ.3) THEN
            this%ud(iBand,:,:,i_elem,ispin) = &
               this%ud(iBand,:,:,i_elem,ispin) + contribution
         ELSE IF(imat.EQ.4) THEN
            this%du(iBand,:,:,i_elem,ispin) = &
               this%du(iBand,:,:,i_elem,ispin) + contribution
         ELSE IF((imat-4.0)/2.0<=nLO) THEN
            !imat 5 - 4+2*numberofLOs: coefficients for Valence-LO contribution
            iLO = CEILING(REAL(imat-4.0)/2.0)
            IF(MOD(imat-4,2)==1) THEN
               this%uulo(iBand,:,:,iLO,i_elemLO,ispin) = &
                  this%uulo(iBand,:,:,iLO,i_elemLO,ispin) + contribution
            ELSE IF(MOD(imat-4,2)==0) THEN
               this%dulo(iBand,:,:,iLO,i_elemLO,ispin) = &
                  this%dulo(iBand,:,:,iLO,i_elemLO,ispin) + contribution
            ENDIF
         ELSE IF((imat-4.0)/2.0<=2.0*nLO) THEN
            !imat 4+2*numberofLOs+1 - 4+4*numberofLOs: coefficients for LO-Valence contribution
            iLO = CEILING(REAL(imat-4.0-2*nLO)/2.0)
            IF(MOD(imat-4-2*nLO,2)==1) THEN
               this%ulou(iBand,:,:,iLO,i_elemLO,ispin) = &
                  this%ulou(iBand,:,:,iLO,i_elemLO,ispin) + contribution
            ELSE IF(MOD(imat-4-2*nLO,2)==0) THEN
               this%ulod(iBand,:,:,iLO,i_elemLO,ispin) = &
                  this%ulod(iBand,:,:,iLO,i_elemLO,ispin) + contribution
            ENDIF
         ELSE
            !imat 4+4*numberofLOs+1 - 4+4*numberofLOs+numberofLOs**2: coefficients for LO-LO contribution
            iLO = imat - 4 - 4*nLO
            this%uloulop(iBand,:,:,iLO,i_elemLO,ispin) = &
                  this%uloulop(iBand,:,:,iLO,i_elemLO,ispin) + contribution
         ENDIF

      END SUBROUTINE greensfBZintCoeffs_add_contribution

      SUBROUTINE greensfImagPart_init(this,gfinp,atoms,input,noco,l_calc,nkpts)

         CLASS(t_greensfImagPart),  INTENT(INOUT)  :: this
         TYPE(t_gfinp),             INTENT(IN)     :: gfinp
         TYPE(t_atoms),             INTENT(IN)     :: atoms
         TYPE(t_input),             INTENT(IN)     :: input
         TYPE(t_noco),              INTENT(IN)     :: noco
         LOGICAL,                   INTENT(IN)     :: l_calc
         INTEGER,                   INTENT(IN)     :: nkpts

         INTEGER lmax,spin_dim,uniqueElementsSphavg,uniqueElementsRadial,uniqueElementsLO,maxLO
         INTEGER uniqueElementsSphavg_kres

         spin_dim = MERGE(3,input%jspins,gfinp%l_mperp)
         lmax = lmaxU_const

         this%l_calc = l_calc

          !Determine number of unique gf elements
         uniqueElementsSphavg  = gfinp%uniqueElements(atoms,l_sphavg=.TRUE.,l_kresolved_int=.FALSE.) !How many spherically averaged elements
         uniqueElementsRadial  = gfinp%uniqueElements(atoms,l_sphavg=.FALSE.,l_kresolved_int=.FALSE.) !How many elements with radial dependence

         uniqueElementsSphavg_kres = gfinp%uniqueElements(atoms,l_sphavg=.TRUE.,l_kresolved_int=.TRUE.) !How many spherically averaged elements with k-resolution

         ALLOCATE (this%kkintgr_cutoff(gfinp%n,input%jspins,2),source=0)
         IF(uniqueElementsSphavg>0) THEN
            ALLOCATE (this%sphavg(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsSphavg,spin_dim),source=cmplx_0)
            ALLOCATE (this%scalingFactorSphavg(uniqueElementsSphavg,input%jspins),source=1.0)
         ENDIF
         IF(uniqueElementsRadial>0) THEN
            ALLOCATE (this%uu(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,spin_dim),source=cmplx_0)
            ALLOCATE (this%dd(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,spin_dim),source=cmplx_0)
            ALLOCATE (this%du(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,spin_dim),source=cmplx_0)
            ALLOCATE (this%ud(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,spin_dim),source=cmplx_0)
            ALLOCATE (this%scalingFactorRadial(uniqueElementsRadial,input%jspins),source=1.0)

            uniqueElementsLO = gfinp%uniqueElements(atoms,lo=.TRUE.,l_sphavg=.FALSE.,maxLO=maxLO, l_kresolved_int=.FALSE.)
            IF(uniqueElementsLO>0) THEN
               ALLOCATE (this%uulo(gfinp%ne,-lmax:lmax,-lmax:lmax,maxLO,uniqueElementsLO,spin_dim),source=cmplx_0)
               ALLOCATE (this%ulou(gfinp%ne,-lmax:lmax,-lmax:lmax,maxLO,uniqueElementsLO,spin_dim),source=cmplx_0)
               ALLOCATE (this%dulo(gfinp%ne,-lmax:lmax,-lmax:lmax,maxLO,uniqueElementsLO,spin_dim),source=cmplx_0)
               ALLOCATE (this%ulod(gfinp%ne,-lmax:lmax,-lmax:lmax,maxLO,uniqueElementsLO,spin_dim),source=cmplx_0)

               ALLOCATE (this%uloulop(gfinp%ne,-lmax:lmax,-lmax:lmax,maxLO,maxLO,uniqueElementsLO,spin_dim),source=cmplx_0)
            ENDIF
         ENDIF
         IF(uniqueElementsSphavg_kres>0) THEN
            ALLOCATE (this%sphavg_k(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsSphavg_kres,spin_dim,nkpts),source=cmplx_0)
            ALLOCATE (this%scalingFactorSphavgKres(uniqueElementsSphavg_kres,input%jspins),source=1.0)
         ENDIF

      END SUBROUTINE greensfImagPart_init

      SUBROUTINE greensfImagPart_collect(this,spin_ind,mpi_communicator)

#ifdef CPP_MPI
         USE mpi
#endif

         CLASS(t_greensfImagPart),     INTENT(INOUT) :: this
         INTEGER,                      INTENT(IN)    :: spin_ind
         INTEGER,                      INTENT(IN)    :: mpi_communicator
#ifdef CPP_MPI
         INTEGER:: ierr,n
         COMPLEX,ALLOCATABLE::ctmp(:)

         IF(ALLOCATED(this%sphavg)) THEN
            n = SIZE(this%sphavg,1)*SIZE(this%sphavg,2)*SIZE(this%sphavg,3)*SIZE(this%sphavg,4)
            ALLOCATE(ctmp(n))
            CALL MPI_ALLREDUCE(this%sphavg(:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%sphavg(:,:,:,:,spin_ind),1)
            DEALLOCATE(ctmp)
         ENDIF
         IF(ALLOCATED(this%uu)) THEN
            n = SIZE(this%uu,1)*SIZE(this%uu,2)*SIZE(this%uu,3)*SIZE(this%uu,4)
            ALLOCATE(ctmp(n))
            CALL MPI_ALLREDUCE(this%uu(:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%uu(:,:,:,:,spin_ind),1)
            CALL MPI_ALLREDUCE(this%ud(:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%ud(:,:,:,:,spin_ind),1)
            CALL MPI_ALLREDUCE(this%du(:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%du(:,:,:,:,spin_ind),1)
            CALL MPI_ALLREDUCE(this%dd(:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%dd(:,:,:,:,spin_ind),1)
            DEALLOCATE(ctmp)
         ENDIF
         IF(ALLOCATED(this%uulo)) THEN
            n = SIZE(this%uulo,1)*SIZE(this%uulo,2)*SIZE(this%uulo,3)*SIZE(this%uulo,4)*SIZE(this%uulo,5)
            ALLOCATE(ctmp(n))
            CALL MPI_ALLREDUCE(this%uulo(:,:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%uulo(:,:,:,:,:,spin_ind),1)
            CALL MPI_ALLREDUCE(this%ulou(:,:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%ulou(:,:,:,:,:,spin_ind),1)
            CALL MPI_ALLREDUCE(this%dulo(:,:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%dulo(:,:,:,:,:,spin_ind),1)
            CALL MPI_ALLREDUCE(this%ulod(:,:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%ulod(:,:,:,:,:,spin_ind),1)
            DEALLOCATE(ctmp)
         ENDIF
         IF(ALLOCATED(this%uloulop)) THEN
            n = SIZE(this%uloulop,1)*SIZE(this%uloulop,2)*SIZE(this%uloulop,3)*SIZE(this%uloulop,4)&
               *SIZE(this%uloulop,5)*SIZE(this%uloulop,6)
            ALLOCATE(ctmp(n))
            CALL MPI_ALLREDUCE(this%uloulop(:,:,:,:,:,:,spin_ind),ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL zcopy(n,ctmp,1,this%uloulop(:,:,:,:,:,:,spin_ind),1)
            DEALLOCATE(ctmp)
         ENDIF
#endif

      END SUBROUTINE greensfImagPart_collect

      SUBROUTINE greensfImagPart_mpi_bc(this,mpi_comm,irank)
         USE m_mpi_bc_tool
         CLASS(t_greensfImagPart), INTENT(INOUT)::this
         INTEGER, INTENT(IN):: mpi_comm
         INTEGER, INTENT(IN), OPTIONAL::irank
         INTEGER ::rank
         IF (PRESENT(irank)) THEN
            rank = irank
         ELSE
            rank = 0
         END IF

         CALL mpi_bc(this%l_calc,rank,mpi_comm)

         IF(ALLOCATED(this%kkintgr_cutoff)) CALL mpi_bc(this%kkintgr_cutoff,rank,mpi_comm)
         IF(ALLOCATED(this%scalingFactorSphavg)) CALL mpi_bc(this%scalingFactorSphavg,rank,mpi_comm)
         IF(ALLOCATED(this%scalingFactorRadial)) CALL mpi_bc(this%scalingFactorRadial,rank,mpi_comm)
         IF(ALLOCATED(this%scalingFactorSphavgKres)) CALL mpi_bc(this%scalingFactorSphavgKres,rank,mpi_comm)
         IF(ALLOCATED(this%sphavg)) CALL mpi_bc(this%sphavg,rank,mpi_comm)
         IF(ALLOCATED(this%uu)) CALL mpi_bc(this%uu,rank,mpi_comm)
         IF(ALLOCATED(this%ud)) CALL mpi_bc(this%ud,rank,mpi_comm)
         IF(ALLOCATED(this%du)) CALL mpi_bc(this%du,rank,mpi_comm)
         IF(ALLOCATED(this%dd)) CALL mpi_bc(this%dd,rank,mpi_comm)
         IF(ALLOCATED(this%uulo)) CALL mpi_bc(this%uulo,rank,mpi_comm)
         IF(ALLOCATED(this%ulou)) CALL mpi_bc(this%ulou,rank,mpi_comm)
         IF(ALLOCATED(this%dulo)) CALL mpi_bc(this%dulo,rank,mpi_comm)
         IF(ALLOCATED(this%ulod)) CALL mpi_bc(this%ulod,rank,mpi_comm)
         IF(ALLOCATED(this%uloulop)) CALL mpi_bc(this%uloulop,rank,mpi_comm)

      END SUBROUTINE greensfImagPart_mpi_bc

      SUBROUTINE greensfImagPart_scale(this,i_elem,i_elemLO,l_sphavg,nLO,k_resolved)

         CLASS(t_greensfImagPart), INTENT(INOUT):: this
         INTEGER,                  INTENT(IN)   :: i_elem
         INTEGER,                  INTENT(IN)   :: i_elemLO
         LOGICAL,                  INTENT(IN)   :: l_sphavg
         INTEGER,                  INTENT(IN)   :: nLO
         LOGICAL, OPTIONAL,        INTENT(IN)   :: k_resolved

         INTEGER :: jspin
         LOGICAL :: k_resolved_arg

         k_resolved_arg = .FALSE.
         IF(PRESENT(k_resolved)) k_resolved_arg = k_resolved

         IF(l_sphavg.AND..NOT.k_resolved_arg) THEN
            IF(ALLOCATED(this%sphavg)) THEN
               IF(SIZE(this%sphavg,5)==2) THEN
                  DO jspin = 1, SIZE(this%sphavg,5)
                     this%sphavg(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin) = this%scalingFactorSphavg(i_elem,jspin) &
                                                                              * this%sphavg(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin)
                  ENDDO
               ENDIF
            ENDIF
         ELSE IF(.NOT.k_resolved_arg) THEN
            IF(ALLOCATED(this%uu)) THEN
               IF(SIZE(this%uu,5)==2) THEN
                  DO jspin = 1, SIZE(this%uu,5)
                     this%uu(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin) = this%scalingFactorRadial(i_elem,jspin) &
                                                                          * this%uu(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin)
                     this%dd(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin) = this%scalingFactorRadial(i_elem,jspin) &
                                                                          * this%dd(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin)
                     this%ud(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin) = this%scalingFactorRadial(i_elem,jspin) &
                                                                          * this%ud(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin)
                     this%du(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin) = this%scalingFactorRadial(i_elem,jspin) &
                                                                          * this%du(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin)
                  ENDDO
               ENDIF
            ENDIF
            IF(ALLOCATED(this%uulo)) THEN
               IF(nLO>0) THEN
                  IF(SIZE(this%uulo,6)==2) THEN
                     DO jspin = 1, SIZE(this%uulo,6)
                        this%uulo(:,-lmaxU_const:,-lmaxU_const:,:,i_elemLO,jspin) = this%scalingFactorRadial(i_elem,jspin) &
                                                                                   * this%uulo(:,-lmaxU_const:,-lmaxU_const:,:,i_elemLO,jspin)
                        this%ulou(:,-lmaxU_const:,-lmaxU_const:,:,i_elemLO,jspin) = this%scalingFactorRadial(i_elem,jspin) &
                                                                                   * this%ulou(:,-lmaxU_const:,-lmaxU_const:,:,i_elemLO,jspin)
                        this%dulo(:,-lmaxU_const:,-lmaxU_const:,:,i_elemLO,jspin) = this%scalingFactorRadial(i_elem,jspin) &
                                                                                   * this%dulo(:,-lmaxU_const:,-lmaxU_const:,:,i_elemLO,jspin)
                        this%ulod(:,-lmaxU_const:,-lmaxU_const:,:,i_elemLO,jspin) = this%scalingFactorRadial(i_elem,jspin) &
                                                                                   * this%ulod(:,-lmaxU_const:,-lmaxU_const:,:,i_elemLO,jspin)

                        this%uloulop(:,-lmaxU_const:,-lmaxU_const:,:,:,i_elemLO,jspin) = this%scalingFactorRadial(i_elem,jspin) &
                                                                                        * this%uloulop(:,-lmaxU_const:,-lmaxU_const:,:,:,i_elemLO,jspin)
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            IF(ALLOCATED(this%sphavg_k)) THEN
               IF(SIZE(this%sphavg_k,5)==2) THEN
                  DO jspin = 1, SIZE(this%sphavg_k,5)
                     this%sphavg_k(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin,:) = this%scalingFactorSphavgKres(i_elem,jspin) &
                                                                                  * this%sphavg_k(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin,:)
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

      END SUBROUTINE greensfImagPart_scale

      PURE FUNCTION greensfImagPart_applyCutoff(this,i_elem,i_gf,spin,l_sphavg,imat,iLO,iLOp,ikpt) Result(imagpartCut)

         CLASS(t_greensfImagPart), INTENT(IN)   :: this
         INTEGER,                  INTENT(IN)   :: i_elem
         INTEGER,                  INTENT(IN)   :: i_gf
         INTEGER,                  INTENT(IN)   :: spin
         LOGICAL,                  INTENT(IN)   :: l_sphavg
         INTEGER, OPTIONAL,        INTENT(IN)   :: imat !which radial dependence array
         INTEGER, OPTIONAL,        INTENT(IN)   :: iLO,iLOp !which local orbitals
         INTEGER, OPTIONAL,        INTENT(IN)   :: ikpt

         COMPLEX, ALLOCATABLE :: imagpartCut(:,:,:)

         INTEGER :: spin_ind, kkcut, ne

         ne = -1
         IF(ALLOCATED(this%sphavg)) THEN
            ne = SIZE(this%sphavg,1)
         ELSE IF(ALLOCATED(this%uu)) THEN
            ne = SIZE(this%uu,1)
         ELSE IF(ALLOCATED(this%sphavg_k)) THEN
            ne = SIZE(this%sphavg_k,1)
         ENDIF
         IF(ne<0) RETURN

         IF(.NOT.ALLOCATED(imagpartCut)) ALLOCATE(imagpartCut(ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),source=cmplx_0)


         IF(PRESENT(ikpt)) THEN
            IF(ALLOCATED(this%sphavg_k)) THEN
               imagpartCut = this%sphavg_k(:,:,:,i_elem,spin,ikpt)
            ENDIF
         ELSE IF(l_sphavg) THEN
            IF(ALLOCATED(this%sphavg)) THEN
               imagpartCut = this%sphavg(:,:,:,i_elem,spin)
            ENDIF
         ELSE IF(.NOT.PRESENT(iLO).AND..NOT.PRESENT(iLOp)) THEN
            !Valence-Valence arrays
            IF(ALLOCATED(this%uu)) THEN
               IF(PRESENT(imat)) THEN
                  IF(imat.EQ.1) THEN
                     imagpartCut = this%uu(:,:,:,i_elem,spin)
                  ELSE IF(imat.EQ.2) THEN
                     imagpartCut = this%dd(:,:,:,i_elem,spin)
                  ELSE IF(imat.EQ.3) THEN
                     imagpartCut = this%ud(:,:,:,i_elem,spin)
                  ELSE IF(imat.EQ.4) THEN
                     imagpartCut = this%du(:,:,:,i_elem,spin)
                  ENDIF
               ENDIF
            ENDIF
         ELSE IF(.NOT.PRESENT(iLOp)) THEN
            !LO-Valence arrays
            IF(ALLOCATED(this%uulo)) THEN
               IF(PRESENT(imat)) THEN
                  IF(imat.EQ.1) THEN
                     imagpartCut = this%uulo(:,:,:,iLO,i_elem,spin)
                  ELSE IF(imat.EQ.2) THEN
                     imagpartCut = this%ulou(:,:,:,iLO,i_elem,spin)
                  ELSE IF(imat.EQ.3) THEN
                     imagpartCut = this%dulo(:,:,:,iLO,i_elem,spin)
                  ELSE IF(imat.EQ.4) THEN
                     imagpartCut = this%ulod(:,:,:,iLO,i_elem,spin)
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            !LO-LO arrays
            IF(ALLOCATED(this%uloulop)) THEN
               imagpartCut = this%uloulop(:,:,:,iLO,iLOp,i_elem,spin)
            ENDIF
         ENDIF

         IF(ALLOCATED(imagpartCut)) THEN
            !Apply Cutoff
            spin_ind = MERGE(1,spin,spin>2)
            kkcut = this%kkintgr_cutoff(i_gf,spin_ind,2)
            IF(kkcut.ne.SIZE(imagpartCut,1)) imagpartCut(kkcut+1:,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const) = cmplx_0
         ENDIF

      END FUNCTION greensfImagPart_applyCutoff

      PURE FUNCTION greensfImagPart_checkEmpty(this,i_elem,i_elemLO,nLO,spin,l_sphavg,ikpt) Result(l_empty)

         CLASS(t_greensfImagPart), INTENT(IN)   :: this
         INTEGER,                  INTENT(IN)   :: i_elem
         INTEGER,                  INTENT(IN)   :: i_elemLO
         INTEGER,                  INTENT(IN)   :: nLO
         INTEGER,                  INTENT(IN)   :: spin
         LOGICAL,                  INTENT(IN)   :: l_sphavg
         INTEGER, OPTIONAL,        INTENT(IN)   :: ikpt

         LOGICAL :: l_empty

         IF(PRESENT(ikpt)) THEN
            IF(ALLOCATED(this%sphavg_k)) THEN
               l_empty = ALL(ABS(this%sphavg_k(:,:,:,i_elem,spin,ikpt)).LT.1e-12)
            ENDIF
         ELSE IF(l_sphavg) THEN
            IF(ALLOCATED(this%sphavg)) THEN
               l_empty = ALL(ABS(this%sphavg(:,:,:,i_elem,spin)).LT.1e-12)
            ENDIF
         ELSE
            IF(ALLOCATED(this%uu)) THEN
               l_empty =     ALL(ABS(this%uu(:,:,:,i_elem,spin)).LT.1e-12) &
                        .AND.ALL(ABS(this%dd(:,:,:,i_elem,spin)).LT.1e-12) &
                        .AND.ALL(ABS(this%ud(:,:,:,i_elem,spin)).LT.1e-12) &
                        .AND.ALL(ABS(this%du(:,:,:,i_elem,spin)).LT.1e-12)
               IF(ALLOCATED(this%uulo).AND.nLO>0) THEN
                  l_empty = l_empty .AND. ALL(ABS(this%uulo(:,:,:,:nLO,i_elemLO,spin)).LT.1e-12) &
                           .AND.ALL(ABS(this%ulou(:,:,:,:nLO,i_elemLO,spin)).LT.1e-12) &
                           .AND.ALL(ABS(this%dulo(:,:,:,:nLO,i_elemLO,spin)).LT.1e-12) &
                           .AND.ALL(ABS(this%dulo(:,:,:,:nLO,i_elemLO,spin)).LT.1e-12) &
                           .AND.ALL(ABS(this%uloulop(:,:,:,:nLO,:nLO,i_elemLO,spin)).LT.1e-12)
               ENDIF
            ENDIF
         ENDIF

      END FUNCTION greensfImagPart_checkEmpty

END MODULE m_types_greensfCoeffs
