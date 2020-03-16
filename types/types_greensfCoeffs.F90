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

         CONTAINS
            PROCEDURE, PASS :: init    =>  greensfBZintCoeffs_init
      END TYPE t_greensfBZintCoeffs


      TYPE t_greensfImagPart

         !Contains the imaginary part of the greens function
         INTEGER, ALLOCATABLE :: kkintgr_cutoff(:,:,:)

         REAL, ALLOCATABLE :: sphavg(:,:,:,:,:)

         ! These arrays are only used in the case we want the green's function with radial dependence
         REAL, ALLOCATABLE :: uu(:,:,:,:,:)
         REAL, ALLOCATABLE :: dd(:,:,:,:,:)
         REAL, ALLOCATABLE :: du(:,:,:,:,:)
         REAL, ALLOCATABLE :: ud(:,:,:,:,:)

         CONTAINS
            PROCEDURE, PASS :: init    =>  greensfImagPart_init
            PROCEDURE, PASS :: collect =>  greensfImagPart_collect
      END TYPE t_greensfImagPart

   PUBLIC t_greensfBZintCoeffs, t_greensfImagPart

   CONTAINS

      SUBROUTINE greensfBZintCoeffs_init(this,gfinp,input,jsp_start,jsp_end,nkpts,nbands)

         CLASS(t_greensfBZintCoeffs),  INTENT(INOUT)  :: this
         TYPE(t_gfinp),                INTENT(IN)     :: gfinp
         TYPE(t_input),                INTENT(IN)     :: input
         INTEGER,                      INTENT(IN)     :: jsp_start,jsp_end
         INTEGER,                      INTENT(IN)     :: nkpts,nbands !number of kpts and bands handled by this rank

         INTEGER lmax, uniqueElements

         lmax = lmaxU_const

         !Determine number of unique gf elements
         CALL uniqueElements(gfinp,uniqueElements)

         IF(gfinp%l_sphavg) THEN
            ALLOCATE (this%sphavg(nbands,-lmax:lmax,-lmax:lmax,nkpts,uniqueElements,jsp_start:jsp_end),source=cmplx_0)
         ELSE
            ALLOCATE (this%uu(nbands,-lmax:lmax,-lmax:lmax,nkpts,uniqueElements,jsp_start:jsp_end),source=cmplx_0)
            ALLOCATE (this%dd(nbands,-lmax:lmax,-lmax:lmax,nkpts,uniqueElements,jsp_start:jsp_end),source=cmplx_0)
            ALLOCATE (this%du(nbands,-lmax:lmax,-lmax:lmax,nkpts,uniqueElements,jsp_start:jsp_end),source=cmplx_0)
            ALLOCATE (this%ud(nbands,-lmax:lmax,-lmax:lmax,nkpts,uniqueElements,jsp_start:jsp_end),source=cmplx_0)
         ENDIF

      END SUBROUTINE greensfBZintCoeffs_init


      SUBROUTINE greensfImagPart_init(this,gfinp,input,noco)

         CLASS(t_greensfImagPart),  INTENT(INOUT)  :: this
         TYPE(t_gfinp),             INTENT(IN)     :: gfinp
         TYPE(t_input),             INTENT(IN)     :: input
         TYPE(t_noco),              INTENT(IN)     :: noco

         INTEGER lmax,spin_dim, uniqueElements

         spin_dim = MERGE(3,input%jspins,gfinp%l_mperp)
         lmax = lmaxU_const

         !Determine number of unique gf elements
         CALL uniqueElements(gfinp,uniqueElements)

         ALLOCATE (this%kkintgr_cutoff(gfinp%n,spin_dim,2),source=0)
         IF(gfinp%l_sphavg) THEN
            ALLOCATE (this%sphavg(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElements,spin_dim),source=0.0)
         ELSE
            ALLOCATE (this%uu(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElements,spin_dim),source=0.0)
            ALLOCATE (this%dd(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElements,spin_dim),source=0.0)
            ALLOCATE (this%du(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElements,spin_dim),source=0.0)
            ALLOCATE (this%ud(gfinp%ne,-lmax:lmax,-lmax:lmax,uniqueElements,spin_dim),source=0.0)
         ENDIF

      END SUBROUTINE greensfImagPart_init

      SUBROUTINE greensfImagPart_collect(this,gfinp,spin_ind,mpi_comm)

         CLASS(t_greensfImagPart),     INTENT(INOUT) :: this
         TYPE(t_gfinp),                INTENT(IN)    :: gfinp
         INTEGER,                      INTENT(IN)    :: spin_ind
         INTEGER,                      INTENT(IN)    :: mpi_comm
#ifdef CPP_MPI
         include 'mpif.h'
#include"cpp_double.h"
         INTEGER:: ierr,irank,n
         REAL,ALLOCATABLE::rtmp(:)
         CALL MPI_COMM_RANK(mpi_comm,irank,ierr)

         IF(gfinp%l_sphavg) THEN
            n = SIZE(this%sphavg,1)*SIZE(this%sphavg,2)*SIZE(this%sphavg,3)*SIZE(this%sphavg,4)
            ALLOCATE(rtmp(n))
            CALL MPI_REDUCE(this%sphavg(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF(irank.EQ.0) CALL CPP_BLAS_scopy(n,rtmp,1,this%sphavg(:,:,:,:,spin_ind),1)
         ELSE
            n = SIZE(this%uu,1)*SIZE(this%uu,2)*SIZE(this%uu,3)*SIZE(this%uu,4)
            ALLOCATE(rtmp(n))
            CALL MPI_REDUCE(this%uu(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF(irank.EQ.0) CALL CPP_BLAS_scopy(n,rtmp,1,this%uu(:,:,:,:,spin_ind),1)
            CALL MPI_REDUCE(this%ud(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF(irank.EQ.0) CALL CPP_BLAS_scopy(n,rtmp,1,this%ud(:,:,:,:,spin_ind),1)
            CALL MPI_REDUCE(this%du(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF(irank.EQ.0) CALL CPP_BLAS_scopy(n,rtmp,1,this%du(:,:,:,:,spin_ind),1)
            CALL MPI_REDUCE(this%dd(:,:,:,:,spin_ind),rtmp,n,CPP_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF(irank.EQ.0) CALL CPP_BLAS_scopy(n,rtmp,1,this%dd(:,:,:,:,spin_ind),1)
         ENDIF
         DEALLOCATE(rtmp)
#endif

      END SUBROUTINE greensfImagPart_collect

END MODULE m_types_greensfCoeffs