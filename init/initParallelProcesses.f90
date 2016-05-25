MODULE m_InitParallelProcesses

CONTAINS

SUBROUTINE initParallelProcesses(atoms,vacuum,input,stars,sliceplot,banddos,&
                                 dimension,cell,sym,xcpot,noco,jij,oneD,hybrid,&
                                 kpts,enpara,sphhar,mpi,results,obsolete)

   USE m_types

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   TYPE(t_mpi),      INTENT(INOUT) :: mpi
   TYPE(t_input),    INTENT(INOUT) :: input
   TYPE(t_sym),      INTENT(INOUT) :: sym
   TYPE(t_stars),    INTENT(INOUT) :: stars 
   TYPE(t_atoms),    INTENT(INOUT) :: atoms
   TYPE(t_vacuum),   INTENT(INOUT) :: vacuum
   TYPE(t_kpts),     INTENT(INOUT) :: kpts
   TYPE(t_oneD),     INTENT(INOUT) :: oneD
   TYPE(t_hybrid),   INTENT(INOUT) :: hybrid
   TYPE(t_Jij),      INTENT(INOUT) :: Jij
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

   INTEGER ierr(3)

   EXTERNAL MPI_BCAST

   CALL MPI_BCAST(atoms%ntype,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%ntypd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%nat,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%natd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%nlod,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%lmaxd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%llod,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(atoms%jmtd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sym%nop,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sym%nop2,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%n3d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%n2d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%k1d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%k2d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(stars%k3d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sphhar%nlhd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sphhar%ntypsd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(sphhar%memd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%jspd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%nstd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%nn3d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(dimension%nn2d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(kpts%numSpecialPoints,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(kpts%nkpts,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(kpts%nkptd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(kpts%nkpt,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(input%jspins,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(vacuum%layerd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%k3,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%M,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%n2d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%nop,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(oneD%odd%nn2d,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(obsolete%nwdd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(jij%nqptd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)

   IF (mpi%irank.NE.0) THEN
      ALLOCATE(atoms%nz(atoms%ntype),atoms%zatom(atoms%ntype)) !nz and zatom have the same content!
      ALLOCATE(atoms%jri(atoms%ntype),atoms%dx(atoms%ntype),atoms%rmt(atoms%ntype))
      ALLOCATE(atoms%lmax(atoms%ntype),atoms%nlo(atoms%ntype),atoms%lnonsph(atoms%ntype))
      ALLOCATE(atoms%ncst(atoms%ntype),atoms%lda_u(atoms%ntype))
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

      ALLOCATE(noco%soc_opt(atoms%ntype+2),noco%l_relax(atoms%ntype),noco%b_con(2,atoms%ntype))
      ALLOCATE(noco%alph(atoms%ntype),noco%beta(atoms%ntype))

      ALLOCATE(Jij%alph1(atoms%ntype),Jij%l_magn(atoms%ntype),Jij%M(atoms%ntype))
      ALLOCATE(Jij%magtype(atoms%ntype),Jij%nmagtype(atoms%ntype))

      ALLOCATE(kpts%specialPoints(3,kpts%numSpecialPoints))
      ALLOCATE(kpts%specialPointNames(kpts%numSpecialPoints))
      ALLOCATE(kpts%bk(3,kpts%nkpts))
      ALLOCATE(kpts%weight(kpts%nkpts))
      ALLOCATE(kpts%wtkpt(kpts%nkpt))

      ALLOCATE(enpara%evac0(2,input%jspins))
      ALLOCATE(enpara%lchg_v(2,input%jspins),enpara%skiplo(atoms%ntypd,input%jspins))
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

      ALLOCATE(stars%ig(-stars%k1d:stars%k1d,-stars%k2d:stars%k2d,-stars%k3d:stars%k3d))
      ALLOCATE(stars%ig2(stars%n3d),stars%igz(stars%n3d))
      ALLOCATE(stars%kv2(2,stars%n2d),stars%kv3(3,stars%n3d))
      ALLOCATE(stars%nstr2(stars%n2d),stars%nstr(stars%n3d))
      ALLOCATE(stars%sk2(stars%n2d),stars%sk3(stars%n3d),stars%phi2(stars%n2d))
      ALLOCATE(stars%igfft(0:dimension%nn3d-1,2),stars%igfft2(0:dimension%nn2d-1,2))
      ALLOCATE(stars%rgphs(-stars%k1d:stars%k1d,-stars%k2d:stars%k2d,-stars%k3d:stars%k3d))
      ALLOCATE(stars%pgfft(0:dimension%nn3d-1),stars%pgfft2(0:dimension%nn2d-1))
      ALLOCATE(stars%ufft(0:27*stars%k1d*stars%k2d*stars%k3d-1),stars%ustep(stars%n3d))

      ALLOCATE(results%force(3,atoms%ntype,dimension%jspd))
      ALLOCATE(results%force_old(3,atoms%ntype))

      ALLOCATE(oneD%ig1(-oneD%odd%k3:oneD%odd%k3,-oneD%odd%M:oneD%odd%M))
      ALLOCATE(oneD%kv1(2,oneD%odd%n2d),oneD%nstr1(oneD%odd%n2d))
      ALLOCATE(oneD%ngopr1(atoms%natd),oneD%mrot1(3,3,oneD%odd%nop),oneD%tau1(3,oneD%odd%nop))
      ALLOCATE(oneD%invtab1(oneD%odd%nop),oneD%multab1(oneD%odd%nop,oneD%odd%nop))
      ALLOCATE(oneD%igfft1(0:oneD%odd%nn2d-1,2),oneD%pgfft1(0:oneD%odd%nn2d-1))

      ALLOCATE(hybrid%nindx(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE(hybrid%select1(4,atoms%ntype),hybrid%lcutm1(atoms%ntype))
      ALLOCATE(hybrid%select2(4,atoms%ntype),hybrid%lcutm2(atoms%ntype),hybrid%lcutwf(atoms%ntype))
      ALLOCATE(hybrid%ddist(dimension%jspd))

      IF (xcpot%igrd.NE.0) THEN
         ALLOCATE (stars%ft2_gfx(0:dimension%nn2d-1),stars%ft2_gfy(0:dimension%nn2d-1))
         ALLOCATE (oneD%pgft1x(0:oneD%odd%nn2d-1),oneD%pgft1xx(0:oneD%odd%nn2d-1),&
                   oneD%pgft1xy(0:oneD%odd%nn2d-1),&
                   oneD%pgft1y(0:oneD%odd%nn2d-1),oneD%pgft1yy(0:oneD%odd%nn2d-1))
      ELSE
         ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
         ALLOCATE (oneD%pgft1x(0:1),oneD%pgft1xx(0:1),oneD%pgft1xy(0:1),&
                   oneD%pgft1y(0:1),oneD%pgft1yy(0:1))
      END IF

      oneD%odd%nq2 = oneD%odd%n2d
      atoms%vr0(:)         = 0.0
      jij%M(:)             = 0.0
      jij%l_magn(:)        =.FALSE.
      results%force(:,:,:) = 0.0
      jij%l_wr=.TRUE.
      jij%nqptd=1
      jij%nmagn=1
      jij%mtypes=1
      jij%phnd=1
      hybrid%ddist     = 1.0
      stars%sk2(:) = 0.0
      stars%phi2(:) = 0.0
   END IF

END SUBROUTINE initParallelProcesses

END MODULE m_InitParallelProcesses
