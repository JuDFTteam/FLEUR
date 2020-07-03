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

         COMPLEX, ALLOCATABLE :: sphavg(:,:,:,:,:,:)

         ! These arrays are only used in the case we want the green's function with radial dependence
         COMPLEX, ALLOCATABLE :: uu(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dd(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: du(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ud(:,:,:,:,:,:)

         !LO-Valence Contribution
         COMPLEX, ALLOCATABLE :: uulo(:,:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ulou(:,:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dulo(:,:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ulod(:,:,:,:,:,:,:)

         !LO-LO contribution
         !Here we need to compress the (lo,lop) index pair into one index because PGI allows a max of 7 indices
         COMPLEX, ALLOCATABLE :: uloulop(:,:,:,:,:,:,:)

         CONTAINS
            PROCEDURE, PASS :: init    =>  greensfBZintCoeffs_init
      END TYPE t_greensfBZintCoeffs


      TYPE t_greensfImagPart

         !Contains the imaginary part of the greens function
         INTEGER, ALLOCATABLE :: kkintgr_cutoff(:,:,:)
         REAL   , ALLOCATABLE :: scalingFactorSphavg(:,:)
         REAL   , ALLOCATABLE :: scalingFactorRadial(:,:)
         LOGICAL :: l_calc = .FALSE.

         REAL, ALLOCATABLE :: sphavg(:,:,:,:,:)

         ! These arrays are only used in the case we want the green's function with radial dependence
         REAL, ALLOCATABLE :: uu(:,:,:,:,:)
         REAL, ALLOCATABLE :: dd(:,:,:,:,:)
         REAL, ALLOCATABLE :: du(:,:,:,:,:)
         REAL, ALLOCATABLE :: ud(:,:,:,:,:)

         !LO-Valence Contribution
         REAL, ALLOCATABLE :: uulo(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: ulou(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: dulo(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: ulod(:,:,:,:,:,:)

         !LO-LO contribution
         !Here the (lo,lop) index pair is explicit again
         REAL, ALLOCATABLE :: uloulop(:,:,:,:,:,:,:)

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

      SUBROUTINE greensfBZintCoeffs_init(this,gfinp,atoms,input,noco,jsp_start,jsp_end,nkpts,nbands)

         CLASS(t_greensfBZintCoeffs),  INTENT(INOUT)  :: this
         TYPE(t_gfinp),                INTENT(IN)     :: gfinp
         TYPE(t_atoms),                INTENT(IN)     :: atoms
         TYPE(t_input),                INTENT(IN)     :: input
         TYPE(t_noco),                 INTENT(IN)     :: noco
         INTEGER,                      INTENT(IN)     :: jsp_start,jsp_end
         INTEGER,                      INTENT(IN)     :: nkpts,nbands !number of kpts and bands handled by this rank

         INTEGER lmax, uniqueElementsSphavg,uniqueElementsRadial, maxSpin

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
            ALLOCATE (this%sphavg(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsSphavg,nkpts,jsp_start:maxSpin),source=cmplx_0)
         ENDIF
         IF(uniqueElementsRadial>0) THEN
            ALLOCATE (this%uu(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,nkpts,jsp_start:maxSpin),source=cmplx_0)
            ALLOCATE (this%dd(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,nkpts,jsp_start:maxSpin),source=cmplx_0)
            ALLOCATE (this%du(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,nkpts,jsp_start:maxSpin),source=cmplx_0)
            ALLOCATE (this%ud(nbands,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,nkpts,jsp_start:maxSpin),source=cmplx_0)
         ENDIF

      END SUBROUTINE greensfBZintCoeffs_init


      SUBROUTINE greensfImagPart_init(this,gfinp,atoms,input,noco,l_calc)

         CLASS(t_greensfImagPart),  INTENT(INOUT)  :: this
         TYPE(t_gfinp),             INTENT(IN)     :: gfinp
         TYPE(t_atoms),             INTENT(IN)     :: atoms
         TYPE(t_input),             INTENT(IN)     :: input
         TYPE(t_noco),              INTENT(IN)     :: noco
         LOGICAL,                   INTENT(IN)     :: l_calc

         INTEGER lmax,spin_dim,uniqueElementsSphavg,uniqueElementsRadial

         spin_dim = MERGE(3,input%jspins,gfinp%l_mperp)
         lmax = lmaxU_const

         this%l_calc = l_calc

          !Determine number of unique gf elements
         uniqueElementsSphavg  = gfinp%uniqueElements(atoms,l_sphavg=.TRUE.) !How many spherically averaged elements
         uniqueElementsRadial  = gfinp%uniqueElements(atoms,l_sphavg=.FALSE.) !How many elements with radial dependence

         ALLOCATE (this%kkintgr_cutoff(gfinp%n,input%jspins,2),source=0)
         IF(uniqueElementsSphavg>0) THEN
            ALLOCATE (this%sphavg(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsSphavg,spin_dim),source=0.0)
            ALLOCATE (this%scalingFactorSphavg(uniqueElementsSphavg,input%jspins),source=1.0)
         ENDIF
         IF(uniqueElementsRadial>0) THEN
            ALLOCATE (this%uu(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,spin_dim),source=0.0)
            ALLOCATE (this%dd(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,spin_dim),source=0.0)
            ALLOCATE (this%du(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,spin_dim),source=0.0)
            ALLOCATE (this%ud(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElementsRadial,spin_dim),source=0.0)
            ALLOCATE (this%scalingFactorRadial(uniqueElementsRadial,input%jspins),source=1.0)
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
#include"cpp_double.h"
         INTEGER:: ierr,n
         REAL,ALLOCATABLE::rtmp(:)

         IF(ALLOCATED(this%sphavg)) THEN
            n = SIZE(this%sphavg,1)*SIZE(this%sphavg,2)*SIZE(this%sphavg,3)*SIZE(this%sphavg,4)
            ALLOCATE(rtmp(n))
            CALL MPI_ALLREDUCE(this%sphavg(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_scopy(n,rtmp,1,this%sphavg(:,:,:,:,spin_ind),1)
            DEALLOCATE(rtmp)
         ENDIF
         IF(ALLOCATED(this%uu)) THEN
            n = SIZE(this%uu,1)*SIZE(this%uu,2)*SIZE(this%uu,3)*SIZE(this%uu,4)
            ALLOCATE(rtmp(n))
            CALL MPI_ALLREDUCE(this%uu(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_scopy(n,rtmp,1,this%uu(:,:,:,:,spin_ind),1)
            CALL MPI_ALLREDUCE(this%ud(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_scopy(n,rtmp,1,this%ud(:,:,:,:,spin_ind),1)
            CALL MPI_ALLREDUCE(this%du(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_scopy(n,rtmp,1,this%du(:,:,:,:,spin_ind),1)
            CALL MPI_ALLREDUCE(this%dd(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_scopy(n,rtmp,1,this%dd(:,:,:,:,spin_ind),1)
            DEALLOCATE(rtmp)
         ENDIF
#endif

      END SUBROUTINE greensfImagPart_collect

      SUBROUTINE greensfImagPart_mpi_bc(this,mpi_comm,irank)
         USE m_mpi_bc_tool
         CLASS(t_greensfImagPart), INTENT(INOUT)::this
         INTEGER, INTENT(IN):: mpi_comm
         INTEGER, INTENT(IN), OPTIONAL::irank
         INTEGER ::rank,myrank,n,ierr
         IF (PRESENT(irank)) THEN
            rank = irank
         ELSE
            rank = 0
         END IF

         CALL mpi_bc(this%l_calc,rank,mpi_comm)

         IF(ALLOCATED(this%kkintgr_cutoff)) CALL mpi_bc(this%kkintgr_cutoff,rank,mpi_comm)
         IF(ALLOCATED(this%scalingFactorSphavg)) CALL mpi_bc(this%scalingFactorSphavg,rank,mpi_comm)
         IF(ALLOCATED(this%scalingFactorRadial)) CALL mpi_bc(this%scalingFactorRadial,rank,mpi_comm)
         IF(ALLOCATED(this%sphavg)) CALL mpi_bc(this%sphavg,rank,mpi_comm)
         IF(ALLOCATED(this%uu)) CALL mpi_bc(this%uu,rank,mpi_comm)
         IF(ALLOCATED(this%ud)) CALL mpi_bc(this%ud,rank,mpi_comm)
         IF(ALLOCATED(this%du)) CALL mpi_bc(this%du,rank,mpi_comm)
         IF(ALLOCATED(this%dd)) CALL mpi_bc(this%dd,rank,mpi_comm)

      END SUBROUTINE greensfImagPart_mpi_bc

      SUBROUTINE greensfImagPart_scale(this,i_elem,l_sphavg)

         CLASS(t_greensfImagPart), INTENT(INOUT):: this
         INTEGER,                  INTENT(IN)   :: i_elem
         LOGICAL,                  INTENT(IN)   :: l_sphavg

         INTEGER :: jspin

         IF(l_sphavg) THEN
            IF(ALLOCATED(this%sphavg)) THEN
               IF(SIZE(this%sphavg,5)==2) THEN
                  DO jspin = 1, SIZE(this%sphavg,5)
                     this%sphavg(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin) = this%scalingFactorSphavg(i_elem,jspin) &
                                                                              * this%sphavg(:,-lmaxU_const:,-lmaxU_const:,i_elem,jspin)
                  ENDDO
               ENDIF
            ENDIF
         ELSE
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
         ENDIF

      END SUBROUTINE greensfImagPart_scale

      PURE FUNCTION greensfImagPart_applyCutoff(this,i_elem,i_gf,m,mp,spin,l_sphavg,imat) Result(imagpartCut)

         CLASS(t_greensfImagPart), INTENT(IN)   :: this
         INTEGER,                  INTENT(IN)   :: i_elem
         INTEGER,                  INTENT(IN)   :: i_gf
         INTEGER,                  INTENT(IN)   :: m
         INTEGER,                  INTENT(IN)   :: mp
         INTEGER,                  INTENT(IN)   :: spin
         LOGICAL,                  INTENT(IN)   :: l_sphavg
         INTEGER, OPTIONAL,        INTENT(IN)   :: imat !which radial dependence array

         REAL, ALLOCATABLE :: imagpartCut(:)

         INTEGER :: spin_ind, kkcut

         IF(l_sphavg) THEN
            IF(ALLOCATED(this%sphavg)) THEN
               IF(.NOT.ALLOCATED(imagpartCut)) ALLOCATE(imagpartCut(SIZE(this%sphavg,1)),source=0.0)
               imagpartCut = this%sphavg(:,m,mp,i_elem,spin)
            ENDIF
         ELSE
            IF(ALLOCATED(this%uu)) THEN
               IF(.NOT.ALLOCATED(imagpartCut)) ALLOCATE(imagpartCut(SIZE(this%uu,1)),source=0.0)
               IF(PRESENT(imat)) THEN
                  IF(imat.EQ.1) THEN
                     imagpartCut = this%uu(:,m,mp,i_elem,spin)
                  ELSE IF(imat.EQ.2) THEN
                     imagpartCut = this%dd(:,m,mp,i_elem,spin)
                  ELSE IF(imat.EQ.3) THEN
                     imagpartCut = this%ud(:,m,mp,i_elem,spin)
                  ELSE IF(imat.EQ.4) THEN
                     imagpartCut = this%du(:,m,mp,i_elem,spin)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF(ALLOCATED(imagpartCut)) THEN
            !Apply Cutoff
            spin_ind = MERGE(1,spin,spin>2)
            kkcut = this%kkintgr_cutoff(i_gf,spin_ind,2)

            IF(kkcut.ne.SIZE(imagpartCut)) imagpartCut(kkcut+1:) = 0.0
         ENDIF

      END FUNCTION greensfImagPart_applyCutoff

      PURE FUNCTION greensfImagPart_checkEmpty(this,i_elem,m,mp,spin,l_sphavg) Result(l_empty)

         CLASS(t_greensfImagPart), INTENT(IN)   :: this
         INTEGER,                  INTENT(IN)   :: i_elem
         INTEGER,                  INTENT(IN)   :: m
         INTEGER,                  INTENT(IN)   :: mp
         INTEGER,                  INTENT(IN)   :: spin
         LOGICAL,                  INTENT(IN)   :: l_sphavg

         LOGICAL :: l_empty

         IF(l_sphavg) THEN
            IF(ALLOCATED(this%sphavg)) THEN
               l_empty = ALL(ABS(this%sphavg(:,m,mp,i_elem,spin)).LT.1e-12)
            ENDIF
         ELSE
            IF(ALLOCATED(this%uu)) THEN
               l_empty =     ALL(ABS(this%uu(:,m,mp,i_elem,spin)).LT.1e-12) &
                        .AND.ALL(ABS(this%dd(:,m,mp,i_elem,spin)).LT.1e-12) &
                        .AND.ALL(ABS(this%ud(:,m,mp,i_elem,spin)).LT.1e-12) &
                        .AND.ALL(ABS(this%du(:,m,mp,i_elem,spin)).LT.1e-12)
            ENDIF
         ENDIF

      END FUNCTION greensfImagPart_checkEmpty

END MODULE m_types_greensfCoeffs