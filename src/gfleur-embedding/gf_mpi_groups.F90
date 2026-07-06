!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_mpi_groups
#ifdef CPP_MPI
      USE mpi
#endif
      IMPLICIT NONE
!-----------------------------------------------
! DESC:determine different parallelization levels
!      The GF parallelization works on t_gfmpi: the standard FLEUR
!      t_mpi (gmpi%fmpi, configured "serial-like" for calls into the
!      FLEUR kernels) plus the k x layer x energy sub-communicators.
!                 Daniel Wortmann, (08-03-12)
!-----------------------------------------------
      INTEGER,PARAMETER :: GGT_PARALLEL=0
      INTEGER,PARAMETER :: FULL_PARALLEL=1
      INTEGER,PARAMETER :: UNBALANCED_PARALLEL=2
      CONTAINS
#ifndef CPP_TEST

      !<-- S: gf_setup_mpi_groups(gmpi,layers,kpts,n_energies)
      SUBROUTINE gf_setup_mpi_groups(gmpi,layers,kpts,n_energies)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      USE m_gf_types
      IMPLICIT NONE
      !<--Arguments
      TYPE(t_kpts),INTENT(IN)    :: kpts
      TYPE(t_layers),INTENT(IN)  :: layers
      TYPE(t_gfmpi),INTENT(INOUT):: gmpi !fmpi%mpi_comm must be set
      INTEGER                    :: n_energies
      !>
      INTEGER             :: ierr
#ifdef CPP_MPI
      CALL MPI_COMM_RANK (gmpi%fmpi%mpi_comm,gmpi%fmpi%irank,ierr)
      CALL MPI_COMM_SIZE (gmpi%fmpi%mpi_comm,gmpi%fmpi%isize,ierr)
#else
      gmpi%fmpi%irank = 0
      gmpi%fmpi%isize = 1
#endif
      gmpi%pe0 = gmpi%fmpi%irank == 0

      CALL priv_full(gmpi,kpts,layers,n_energies)

      CALL priv_subcom(gmpi)

      CALL priv_print_info(gmpi,layers%num_layers,kpts%nkpt)

      END SUBROUTINE
      !>

      !<-- S:priv_subcom(gmpi)
      SUBROUTINE priv_subcom(gmpi)
!-----------------------------------------------
!  create a subcom with all members that write on gf_pot and gf_iodop
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      USE m_gf_types
      IMPLICIT NONE
      !<--Arguments
      TYPE(t_gfmpi),INTENT(INOUT) :: gmpi
      !>
      !<-- Locals
      INTEGER             :: iogroupsize,i
      INTEGER,ALLOCATABLE :: members(:)
#ifdef CPP_MPI
      INTEGER             :: WORLD_GROUP,SUB_GROUP,ierr
      !>
      !<-- Group of PE that construct the potential
      iogroupsize = gmpi%k_PEperK
      ALLOCATE(members(iogroupsize))
      members = (/(i,i=0,gmpi%k_PEperK-1)/)


      CALL MPI_COMM_GROUP (gmpi%fmpi%mpi_comm,WORLD_GROUP,ierr)
      CALL MPI_GROUP_INCL (WORLD_GROUP,iogroupsize,members,             &
     &                                      SUB_GROUP,ierr)
      CALL MPI_COMM_CREATE (gmpi%fmpi%mpi_comm,SUB_GROUP,               &
     &     gmpi%iodop_SUBCOM,ierr)
      DEALLOCATE(members)
      !>
      !<-- Group of PE that work on same layer
      iogroupsize = gmpi%fmpi%isize/gmpi%k_PEperK
      ALLOCATE(members(iogroupsize))
      members = (/(i,i = 0,gmpi%fmpi%isize-1,gmpi%k_PEperK)/)
      members = members+MOD(gmpi%fmpi%irank,gmpi%k_PEperK)
      CALL MPI_GROUP_INCL (WORLD_GROUP,iogroupsize,members,             &
     &     SUB_GROUP,ierr)
      CALL MPI_COMM_CREATE (gmpi%fmpi%mpi_comm,SUB_GROUP               &
     &     ,gmpi%samelayer_SUBCOM,ierr)
      DEALLOCATE(members)
      !>
      gmpi%self_subcom = MPI_COMM_SELF
#else
      gmpi%self_subcom = 1
      gmpi%iodop_subcom = 1
      gmpi%samelayer_subcom = 1
#endif
      !configure the embedded standard t_mpi such that FLEUR kernels
      !called per (layer,k) run without their own eigenvalue
      !parallelization
      gmpi%fmpi%n_rank = 0
      gmpi%fmpi%n_size = 1
      gmpi%fmpi%sub_comm = gmpi%self_subcom

      END SUBROUTINE
      !>

      !<-- S: priv_full(gmpi,kpts,layers,n_energies)
      SUBROUTINE PRIV_full(gmpi,kpts,layers,n_energies)
!-----------------------------------------------
!  full parallelization over kpts+(layers or energies)
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      USE m_gf_types
      use m_juDFT
      IMPLICIT NONE
      !<--Arguments
      TYPE(t_kpts),INTENT(IN)     :: kpts
      TYPE(t_layers),INTENT(IN)   :: layers
      TYPE(t_gfmpi),INTENT(INOUT) :: gmpi
      INTEGER                     :: n_energies
      !>
      !<-- Locals
      INTEGER             :: ggt,n

      !>

      !k-point parallelization
      ggt = priv_ggt(kpts%nkpt,gmpi%fmpi%isize)

      gmpi%k_kperPE = kpts%nkpt/ggt
      gmpi%k_PEperK = gmpi%fmpi%isize/ggt

      ALLOCATE( gmpi%k_kpts(gmpi%k_kperPE))
      DO n = 1,gmpi%k_kperPE
         gmpi%k_rank =mod(gmpi%fmpi%irank,gmpi%k_PEperK)
         gmpi%k_kpts(n) = n+gmpi%k_kperPE*(gmpi%fmpi%irank/gmpi%k_PEperK)
      ENDDO


      !additional layer parallelization
      gmpi%kl_LayerPerPE = layers%num_layers/gmpi%k_PEperK

      IF (gmpi%kl_LayerPerPE == 0) CALL juDFT_error                      &
     &     ("Parallelization error 1")

      IF (MOD(layers%num_layers,gmpi%k_PEperK)>0)                        &
     &     gmpi%kl_LayerPerPE = gmpi%kl_LayerPerPE+1
      IF (gmpi%k_rank == gmpi%k_PEperK-1) THEN
         gmpi%kl_LayerPerPE = layers%num_layers-gmpi%kl_LayerPerPE       &
     &        *gmpi%k_rank
         ALLOCATE(gmpi%kl_layers(gmpi%kl_LayerPerPE))
         DO n = 1,gmpi%kl_LayerPerPE
            gmpi%kl_layers(n) = layers%num_layers-gmpi%kl_LayerPerPE+n
         ENDDO
      ELSE
         ALLOCATE(gmpi%kl_layers(gmpi%kl_LayerPerPE))
         DO n = 1,gmpi%kl_LayerPerPE
            gmpi%kl_layers(n) = n+gmpi%kl_LayerPerPE*gmpi%k_rank
         ENDDO
      ENDIF

      !additional energy parallelization

      gmpi%ke_ENPerPE = n_energies/gmpi%k_PEperK

      IF (gmpi%ke_EnPerPE == 0) THEN
         WRITE(*,"('PE:',i3,' Energies:',i8,' PE/kpt:',i8):")         &
            gmpi%fmpi%irank,n_energies,gmpi%k_PEperK
         CALL juDFT_error                           &
     &     ("Parallelization error 2")
      ENDIF
      IF (MOD(n_energies,gmpi%k_PEperK)>0)                               &
     &     gmpi%ke_ENPerPE = gmpi%ke_ENPerPE+1
      IF (gmpi%k_rank == gmpi%k_PEperK-1) THEN
         gmpi%ke_EnPerPE = n_energies-gmpi%ke_EnPerPE                    &
     &        *gmpi%k_rank
         ALLOCATE(gmpi%ke_energies(gmpi%ke_EnPerPE))
         DO n = 1,gmpi%ke_EnPerPE
            gmpi%ke_energies(n) = n_energies-gmpi%ke_EnPerPE+n
         ENDDO
      ELSE
         ALLOCATE(gmpi%ke_energies(gmpi%ke_EnPerPE))
         DO n = 1,gmpi%ke_ENPerPE
            gmpi%ke_energies(n) = n+gmpi%ke_ENPerPE*gmpi%k_rank
         ENDDO
      ENDIF


      END SUBROUTINE
      !>



      !<-- S: priv_print_info(gmpi,nl,nk)
      SUBROUTINE priv_PRINT_info(gmpi,nl,nk)
!-----------------------------------------------
!  prints the information on parallelization
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      USE m_gf_types
      USE m_constants, ONLY: oUnit
      IMPLICIT NONE
      !<--Arguments
                                      !no of layers, no of kpts
      INTEGER,INTENT(IN)     :: nl,nk
      TYPE(t_gfmpi),INTENT(IN) :: gmpi
      !>
      !<-- Locals
      INTEGER             :: n
      INTEGER             :: info_kpts(nk+1)
      INTEGER             :: info_layer(nl+1)
      INTEGER             :: info_com(2)
#ifdef CPP_MPI
                                    !MPI error+status
      INTEGER::e,s(MPI_STATUS_SIZE)

      !>
      !info on kpts for this PE
      info_kpts=0
      info_kpts(1) = gmpi%k_kperPE
      info_kpts(2:gmpi%k_kperPE+1) = gmpi%k_kpts(1:gmpi%k_kperPE)
      !info on layers in second parallelization
      info_layer=0
      info_layer(1) = gmpi%kl_layerperpe
      info_layer(2:gmpi%kl_layerperpe+1) =                               &
     &     gmpi%kl_layers(:gmpi%kl_layerperpe)
       info_com=(/gmpi%iodop_subcom,gmpi%samelayer_subcom/)


      WRITE(oUnit,*)
      WRITE(oUnit,*) "Parallelization setup:"
      IF (gmpi%fmpi%irank >0) THEN
         CALL MPI_SEND(info_kpts,SIZE(info_kpts),                       &
     &        MPI_INTEGER,0,11,gmpi%fmpi%mpi_comm,E)
         CALL MPI_SEND(info_layer,SIZE(info_layer),                     &
     &        MPI_INTEGER,0,22,gmpi%fmpi%mpi_comm,E)
         CALL MPI_SEND(info_com,SIZE(info_com),                     &
     &        MPI_INTEGER,0,33,gmpi%fmpi%mpi_comm,E)
      ELSE
         DO n = 0,gmpi%fmpi%isize-1
            IF (n>0) THEN
               CALL MPI_RECV(info_kpts,SIZE(info_kpts),                 &
     &              MPI_INTEGER,n,11,gmpi%fmpi%mpi_comm,s,E)
               CALL MPI_RECV(info_layer,SIZE(info_layer),               &
     &              MPI_INTEGER,n,22,gmpi%fmpi%mpi_comm,s,E)
               CALL MPI_RECV(info_com,SIZE(info_com),               &
     &              MPI_INTEGER,n,33,gmpi%fmpi%mpi_comm,s,E)

            ENDIF
            WRITE(oUnit,"('For PE:',i6)") n
            WRITE(oUnit,"('Kpts  :',100(i0,1x))")                       &
     &           info_kpts(2:info_kpts(1)+1)
            WRITE(oUnit,"('Layers:',100(i0,1x))")                       &
     &           info_layer(2:info_layer(1)+1)
            write(oUnit,"('Subcoms: iodop=',i0,'  samelayer=',i0)")     &
     &           info_com
         ENDDO
      ENDIF
      CALL MPI_Barrier(gmpi%fmpi%mpi_comm,E)
#endif

      END SUBROUTINE
      !>
#endif
     !<-- F: priv_ggt(n1,n2)
      PURE FUNCTION priv_ggt(n1,n2)
!-----------------------------------------------
! find the larges common denominator
!             (last modified: 08-03-12) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)     :: n1,n2
      INTEGER                :: priv_ggt
      !>
      !<-- Locals
      INTEGER             :: i,ii,iii
      !>
      IF (n1>n2) THEN
         i = n1; ii = n2
      ELSE
         i = n2;ii=n1
      ENDIF

      DO iii = ii,1,-1
         IF (MOD(i,iii) == 0.AND.MOD(ii,iii) == 0) THEN
            priv_ggt = iii
            RETURN
         ENDIF
      ENDDO
      priv_ggt = 1


      END FUNCTION
      !>
#ifdef CPP_MPI
      SUBROUTINE private_parallelize_loop(com,loopsize,strategy,local_elements,subcom)
      !Parallelize a loop of loopsize elements for the current communicator
      !Returns the sub-communicator for the parallelization and the loop-elements
      !that are calculated locally
      !Different strategies are implemented:
      !GGT_PARALLEL : in this case the largest common divisor is used to distribute the Parallize as good as possible
      !FULL_PARALLEL: Here the parallelization must be complete in the sense that either com->isize is a multiple of loopsize or
      !               loopsize a multiple of com->isize, if this fails an invalid subcom and unallocated local_elements is returned
      !UNBALANCED_PARALLEL: if isize<loop each PE does one or more loop-elements
      !                     if isize>loop each loop-element is done by one or more PE
      use m_juDFT
      IMPLICIT NONE
      INTEGER,INTENT(IN)               ::com,loopsize,strategy
      INTEGER,ALLOCATABLE,INTENT(OUT)  ::local_elements(:)
      INTEGER,INTENT(OUT)              ::subcom

      !locals
      INTEGER :: irank,isize,e(MPI_STATUS_SIZE)
      INTEGER :: ggt,n,n_local,n_region
      !Determine current dimensions
      CALL MPI_COMM_RANK(com,irank,e)
      CALL MPI_COMM_SIZE(com,isize,e)

      IF (strategy<UNBALANCED_PARALLEL) THEN
          ggt=priv_ggt(isize,loopsize) !this is the number of parallel regions
          IF (strategy==FULL_PARALLEL) THEN
               !Check if we get full parallelization
               IF ((ggt.ne.isize).and.(ggt.ne.loopsize)) THEN
                    subcom=MPI_COMM_NULL
                    return
               ENDIF
          ENDIF
          n_local=loopsize/ggt      !number of loop-elements per region
          n_region=mod(irank,ggt)   !index of parallel region
          CALL MPI_COMM_SPLIT(com,n_region,irank,subcom,e) !create communicator
          ALLOCATE(local_elements(n_local))
          local_elements=(/(n,n=1,n_local)/)
          local_elements=local_elements+n_region*n_local
          return
       ENDIF
       IF (strategy.ne.UNBALANCED_PARALLEL) CALL juDFT_error("BUG in gf_mpi_groups")

       IF (isize>loopsize) THEN
          !there will be loopsize parallel regions
          ALLOCATE(local_elements(1))
          local_elements(1)=mod(irank,loopsize)
          if (local_elements(1)==0) local_elements(1)=loopsize
          CALL MPI_COMM_SPLIT(com,local_elements(1),irank,subcom,e)
          return
       ELSE
          !there will be isize parallel regions
          n_local=count((/(mod(n,isize)==irank,n=1,loopsize)/))
          ALLOCATE(local_elements(n_local))
          n_local=0
          DO n=1,loopsize
            IF (mod(n,isize)==irank) THEN
                n_local=n_local+1
                local_elements(n_local)=n
            ENDIF
          ENDDO
          CALL MPI_COMM_SPLIT(com,irank,irank,subcom,e)
          RETURN
       ENDIF
       END SUBROUTINE
#endif
#ifdef CPP_TEST
       subroutine priv_plot_comm(com1,com2)
       implicit none
       integer,intent(in)::com1,com2
       integer           ::rank1,rank2,size2,g1,g2,e
       integer,allocatable::r1(:),r2(:)

       CALL MPI_COMM_GROUP(com1,g1,e)
       CALL MPI_COMM_GROUP(com2,g2,e)

       CALL MPI_GROUP_rank(g1,rank1,e)
       CALL MPI_GROUP_rank(g2,rank2,e)
       CALL MPI_GROUP_SIZE(g2,size2,e)
       ALLOCATE(r1(size2),r2(size2))
       r2=(/(e,e=0,size2-1)/)

       CALL MPI_GROUP_translate_ranks(g2,size2,r2,g1,r1,e)

       write(*,*) rank1,r2,'->',r1
       end subroutine


#endif


      END

#ifdef CPP_TEST
program test
use m_gf_mpi_groups
#ifdef CPP_MPI
use mpi
#endif

INTEGER :: loopsize

integer :: com,irank,isize,e,lsize,lrank
integer,allocatable:: locals(:)

!Determine current dimensions
      CALL MPI_INIT(e)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,e)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,e)

IF (irank==0) THEN
    !read(*,*) loopsize
    loopsize=6
ENDIF

CALL MPI_BCAST(loopsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,e)

IF (irank==0) WRITE(*,*) "Test:",isize,"PE for loop of",loopsize

CALL private_parallelize_loop(MPI_COMM_WORLD,loopsize,GGT_PARALLEL,locals,com)

IF (com==MPI_COMM_NULL) THEN
       write(*,*) irank,"GGT failed"
ELSE
      CALL MPI_COMM_RANK(com,lrank,e)
      CALL MPI_COMM_SIZE(com,lsize,e)
       WRITE(*,"(a,99i4)") "GGT:",irank,lrank,lsize,locals

DEALLOCATE(locals)
CALL MPI_COMM_FREE(com,e)
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,e)

CALL private_parallelize_loop(MPI_COMM_WORLD,loopsize,FULL_PARALLEL,locals,com)
IF (com==MPI_COMM_NULL) THEN
       write(*,*) irank,"FULL failed"
ELSE
      CALL MPI_COMM_RANK(com,lrank,e)
      CALL MPI_COMM_SIZE(com,lsize,e)
       WRITE(*,"(a,99i4)") "FULL:",irank,lrank,lsize,locals
      call priv_plot_comm(MPI_COMM_WORLD,com)
DEALLOCATE(locals)
CALL MPI_COMM_FREE(com,e)
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,e)

CALL private_parallelize_loop(MPI_COMM_WORLD,loopsize,UNBALANCED_PARALLEL,locals,com)
IF (com==MPI_COMM_NULL) THEN
       write(*,*) irank,"failed"
ELSE
      CALL MPI_COMM_RANK(com,lrank,e)
      CALL MPI_COMM_SIZE(com,lsize,e)
       WRITE(*,"(a,99i4)") "UNB:",irank,lrank,lsize,locals
       call priv_plot_comm(MPI_COMM_WORLD,com)
DEALLOCATE(locals)
CALL MPI_COMM_FREE(com,e)
ENDIF
CALL MPI_FINALIZE(e)
end
#endif
