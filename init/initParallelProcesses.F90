!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_InitParallelProcesses

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine distributes input data from inp.xml to all parallel
! processes
!
!                                GM'16
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initParallelProcesses(atoms,vacuum,input,stars,sliceplot,banddos,&
                                 dimension,cell,sym,xcpot,noco,oneD,hybrid,&
                                 kpts,enpara,sphhar,mpi,results,obsolete)

   USE m_types

   IMPLICIT NONE

   TYPE(t_mpi),      INTENT(INOUT) :: mpi
   TYPE(t_input),    INTENT(INOUT) :: input
   TYPE(t_sym),      INTENT(INOUT) :: sym
   TYPE(t_stars),    INTENT(INOUT) :: stars 
   TYPE(t_atoms),    INTENT(INOUT) :: atoms
   TYPE(t_vacuum),   INTENT(INOUT) :: vacuum
   TYPE(t_kpts),     INTENT(INOUT) :: kpts
   TYPE(t_oneD),     INTENT(INOUT) :: oneD
   TYPE(t_hybrid),   INTENT(INOUT) :: hybrid
   TYPE(t_cell),     INTENT(INOUT) :: cell
   TYPE(t_banddos),  INTENT(INOUT) :: banddos
   TYPE(t_sliceplot),INTENT(INOUT) :: sliceplot
   TYPE(t_xcpot),    INTENT(INOUT) :: xcpot
   TYPE(t_noco),     INTENT(INOUT) :: noco
   TYPE(t_dimension),INTENT(INOUT) :: dimension
   TYPE(t_enpara),   INTENT(INOUT) :: enpara
   TYPE(t_sphhar),   INTENT(INOUT) :: sphhar
   TYPE(t_results),  INTENT(INOUT) :: results
   TYPE(t_obsolete), INTENT(INOUT) :: obsolete
#ifdef CPP_MPI
   INCLUDE 'mpif.h'

   INTEGER ierr(3)

   EXTERNAL MPI_BCAST

   CALL MPI_BCAST(atoms%ntype,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%ntype,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%nat,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%nat,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%nlod,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%lmaxd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%llod,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%jmtd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sym%nop,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sym%nop2,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%ng3,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%ng2,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%mx1,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%mx2,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%mx3,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%kq1_fft,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%kq2_fft,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%kq3_fft,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sphhar%nlhd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sphhar%ntypsd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sphhar%memd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%jspd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%nstd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%nn3d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%nn2d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%ncvd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%nvd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%neigd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%nv2d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%msh,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%nspd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(kpts%numSpecialPoints,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(kpts%nkpt,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(kpts%nkpt,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(kpts%nkpt,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(kpts%ntet,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(input%jspins,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(vacuum%layerd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(vacuum%nmzxyd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(vacuum%nmzd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%k3,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%M,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%mb,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%n2d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%nop,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%nn2d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   !CALL MPI_BCAST(obsolete%nwdd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)

   IF (mpi%irank.NE.0) THEN
      IF(ALLOCATED(atoms%neq)) DEALLOCATE(atoms%neq)
      IF(ALLOCATED(atoms%volmts)) DEALLOCATE(atoms%volmts)
      IF(ALLOCATED(atoms%taual)) DEALLOCATE(atoms%taual)
      IF(ALLOCATED(atoms%rmt)) DEALLOCATE(atoms%rmt)
      ALLOCATE(atoms%nz(atoms%ntype),atoms%zatom(atoms%ntype)) !nz and zatom have the same content!
      ALLOCATE(atoms%jri(atoms%ntype),atoms%dx(atoms%ntype),atoms%rmt(atoms%ntype))
      ALLOCATE(atoms%lmax(atoms%ntype),atoms%nlo(atoms%ntype),atoms%lnonsph(atoms%ntype))
      ALLOCATE(atoms%ncst(atoms%ntype),atoms%lda_u(4*atoms%ntype))
      ALLOCATE(atoms%nflip(atoms%ntype),atoms%bmu(atoms%ntype),atoms%neq(atoms%ntype))
      ALLOCATE(atoms%l_geo(atoms%ntype),atoms%relax(3,atoms%ntype))
      ALLOCATE(atoms%taual(3,atoms%nat),atoms%pos(3,atoms%nat))
      ALLOCATE(atoms%numStatesProvided(atoms%ntype))
      ALLOCATE(atoms%rmsh(atoms%jmtd,atoms%ntype))
      ALLOCATE(atoms%volmts(atoms%ntype))
      ALLOCATE(atoms%vr0(atoms%ntype))  ! This should actually not be in the atoms type!

      ALLOCATE(atoms%ncv(atoms%ntype))
      ALLOCATE(atoms%ngopr(atoms%nat))
      ALLOCATE(atoms%lapw_l(atoms%ntype))
      ALLOCATE(atoms%invsat(atoms%nat))
      ALLOCATE(atoms%nlhtyp(atoms%ntype),atoms%ntypsy(atoms%nat))

      ALLOCATE(sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
      ALLOCATE(sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd))
      ALLOCATE(sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
      ALLOCATE(sphhar%nlh(sphhar%ntypsd),sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd))

      ALLOCATE(noco%l_relax(atoms%ntype),noco%b_con(2,atoms%ntype))
      ALLOCATE(noco%alphInit(atoms%ntype),noco%alph(atoms%ntype),noco%beta(atoms%ntype))

  
      ALLOCATE(kpts%specialPoints(3,kpts%numSpecialPoints))
      ALLOCATE(kpts%specialPointNames(kpts%numSpecialPoints))
      ALLOCATE(kpts%bk(3,kpts%nkpt))
      ALLOCATE(kpts%wtkpt(kpts%nkpt))
      ALLOCATE(kpts%ntetra(4,kpts%ntet))
      ALLOCATE(kpts%voltet(kpts%ntet))

      ALLOCATE(enpara%evac0(2,input%jspins))
      ALLOCATE(enpara%lchg_v(2,input%jspins),enpara%skiplo(atoms%ntype,input%jspins))
      ALLOCATE(enpara%enmix(input%jspins))

      ALLOCATE(sym%mrot(3,3,sym%nop),sym%tau(3,sym%nop))
      ALLOCATE(sym%invarop(atoms%nat,sym%nop),sym%invarind(atoms%nat))
      ALLOCATE(sym%multab(sym%nop,sym%nop),sym%invtab(sym%nop))
      ALLOCATE(sym%invsatnr(atoms%nat),sym%d_wgn(-3:3,-3:3,3,sym%nop))

      ALLOCATE(atoms%llo(atoms%nlod,atoms%ntype))
      ALLOCATE(atoms%ulo_der(atoms%nlod,atoms%ntype))
      ALLOCATE(enpara%ello0(atoms%nlod,atoms%ntype,input%jspins))
      ALLOCATE(enpara%llochg(atoms%nlod,atoms%ntype,input%jspins))
      ALLOCATE(enpara%el0(0:atoms%lmaxd,atoms%ntype,input%jspins))
      ALLOCATE(enpara%lchange(0:atoms%lmaxd,atoms%ntype,input%jspins))
      ALLOCATE(atoms%l_dulo(atoms%nlod,atoms%ntype)) ! For what is this?
      ALLOCATE(atoms%lo1l(0:atoms%llod,atoms%ntype))
      ALLOCATE(atoms%nlol(0:atoms%llod,atoms%ntype))

      ALLOCATE(atoms%coreStateOccs(dimension%nstd,2,atoms%ntype))
      ALLOCATE(atoms%coreStateNprnc(dimension%nstd,atoms%ntype))
      ALLOCATE(atoms%coreStateKappa(dimension%nstd,atoms%ntype))

      ALLOCATE(vacuum%izlay(vacuum%layerd,2))

      ALLOCATE(stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
      ALLOCATE(stars%ig2(stars%ng3))
      ALLOCATE(stars%kv2(2,stars%ng2),stars%kv3(3,stars%ng3))
      ALLOCATE(stars%nstr2(stars%ng2),stars%nstr(stars%ng3))
      ALLOCATE(stars%sk2(stars%ng2),stars%sk3(stars%ng3),stars%phi2(stars%ng2))
      ALLOCATE(stars%igfft(0:dimension%nn3d-1,2),stars%igfft2(0:dimension%nn2d-1,2))
      ALLOCATE(stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
      ALLOCATE(stars%pgfft(0:dimension%nn3d-1),stars%pgfft2(0:dimension%nn2d-1))
      IF(ALLOCATED(stars%ufft)) DEALLOCATE(stars%ufft)
      ALLOCATE(stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1),stars%ustep(stars%ng3))

      ALLOCATE(results%force(3,atoms%ntype,dimension%jspd))
      ALLOCATE(results%force_old(3,atoms%ntype))

      ALLOCATE(oneD%ig1(-oneD%odd%k3:oneD%odd%k3,-oneD%odd%M:oneD%odd%M))
      ALLOCATE(oneD%kv1(2,oneD%odd%n2d),oneD%nstr1(oneD%odd%n2d))
      ALLOCATE(oneD%ngopr1(atoms%nat),oneD%mrot1(3,3,oneD%odd%nop),oneD%tau1(3,oneD%odd%nop))
      ALLOCATE(oneD%invtab1(oneD%odd%nop),oneD%multab1(oneD%odd%nop,oneD%odd%nop))
      ALLOCATE(oneD%igfft1(0:oneD%odd%nn2d-1,2),oneD%pgfft1(0:oneD%odd%nn2d-1))

      ALLOCATE(hybrid%nindx(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE(hybrid%select1(4,atoms%ntype),hybrid%lcutm1(atoms%ntype))
      ALLOCATE(hybrid%lcutwf(atoms%ntype))

      IF (xcpot%is_gga()) THEN
         ALLOCATE (stars%ft2_gfx(0:dimension%nn2d-1),stars%ft2_gfy(0:dimension%nn2d-1))
         ALLOCATE (oneD%pgft1x(0:oneD%odd%nn2d-1),oneD%pgft1xx(0:oneD%odd%nn2d-1),&
                   oneD%pgft1xy(0:oneD%odd%nn2d-1),&
                   oneD%pgft1y(0:oneD%odd%nn2d-1),oneD%pgft1yy(0:oneD%odd%nn2d-1))
      ELSE
         ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
         ALLOCATE (oneD%pgft1x(0:1),oneD%pgft1xx(0:1),oneD%pgft1xy(0:1),&
                   oneD%pgft1y(0:1),oneD%pgft1yy(0:1))
      END IF

      ! Explicit atom-dependent xc functional
      ALLOCATE(atoms%namex(atoms%ntype))
      ALLOCATE(atoms%relcor(atoms%ntype))
      ALLOCATE(atoms%icorr(atoms%ntype))
      ALLOCATE(atoms%krla(atoms%ntype))
      atoms%namex = ''
      atoms%icorr = -99

      oneD%odd%nq2 = oneD%odd%n2d
      atoms%vr0(:)         = 0.0
      results%force(:,:,:) = 0.0
      stars%sk2(:) = 0.0
      stars%phi2(:) = 0.0
   END IF
#endif
END SUBROUTINE initParallelProcesses

END MODULE m_InitParallelProcesses
