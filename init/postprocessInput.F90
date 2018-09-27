!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_postprocessInput

CONTAINS

SUBROUTINE postprocessInput(mpi,input,field,sym,stars,atoms,vacuum,obsolete,kpts,&
     oneD,hybrid,cell,banddos,sliceplot,xcpot,forcetheo,&
     noco,DIMENSION,enpara,sphhar,l_opti,noel,l_kpts)

  USE m_juDFT
  USE m_types
  USE m_constants
  USE m_lapwdim
  USE m_ylm
  USE m_convndim
  USE m_chkmt
  USE m_localsym
  USE m_strgndim
  USE m_od_chisym
  USE m_dwigner
  USE m_mapatom
  USE m_cdn_io
  USE m_strgn
  USE m_od_strgn1
  USE m_prpqfft
  USE m_writegw
  USE m_prpxcfft
  USE m_stepf
  USE m_convn
  USE m_efield
  USE m_od_mapatom
  USE m_od_kptsgen
  USE m_gen_bz
  USE m_nocoInputCheck
  USE m_kpoints
  IMPLICIT NONE

  TYPE(t_mpi)      ,INTENT   (IN) :: mpi
  CLASS(t_forcetheo),ALLOCATABLE,INTENT(IN)   :: forcetheo
  TYPE(t_input),    INTENT(INOUT) :: input
  TYPE(t_sym),      INTENT(INOUT) :: sym
  TYPE(t_stars),    INTENT(INOUT) :: stars 
  TYPE(t_atoms),    INTENT(INOUT) :: atoms
  TYPE(t_vacuum),   INTENT(INOUT) :: vacuum
  TYPE(t_obsolete), INTENT(INOUT) :: obsolete
  TYPE(t_kpts),     INTENT(INOUT) :: kpts
  TYPE(t_oneD),     INTENT(INOUT) :: oneD
  TYPE(t_hybrid),   INTENT(INOUT) :: hybrid
  TYPE(t_cell),     INTENT(INOUT) :: cell
  TYPE(t_banddos),  INTENT(INOUT) :: banddos
  TYPE(t_sliceplot),INTENT(INOUT) :: sliceplot
  CLASS(t_xcpot),   INTENT(INOUT) :: xcpot
  TYPE(t_noco),     INTENT(INOUT) :: noco
  TYPE(t_dimension),INTENT(INOUT) :: dimension
  TYPE(t_enpara)   ,INTENT(INOUT) :: enpara
  TYPE(t_sphhar)   ,INTENT  (OUT) :: sphhar
  TYPE(t_field),    INTENT(INOUT) :: field
  LOGICAL,          INTENT  (OUT) :: l_opti
  LOGICAL,          INTENT   (IN) :: l_kpts
  CHARACTER(len=3), ALLOCATABLE, INTENT(IN) :: noel(:)

  INTEGER              :: i, j, n, na, n1, n2, iType, l, ilo, ikpt
  INTEGER              :: minNeigd, nv, nv2, kq1, kq2, kq3, jrc, jsp, ii
  INTEGER              :: ios, ntst, ierr
  REAL                 :: sumWeight, rmtmax, zp, radius, dr
  REAL                 :: kmax1, dtild1, dvac1
  REAL                 :: bk(3)
  LOGICAL              :: l_vca, l_test,l_gga
 
  INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)
  INTEGER, ALLOCATABLE :: jri1(:), lmax1(:)
  REAL,    ALLOCATABLE :: rmt1(:), dx1(:)

#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of input postprocessing (calculate missing parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (mpi%irank.EQ.0) THEN

  ! Check the LO stuff and call setlomap (from inped):

     ALLOCATE(atoms%lo1l(0:atoms%llod,atoms%ntype))
     ALLOCATE(atoms%nlol(0:atoms%llod,atoms%ntype))

     IF (ANY(atoms%lapw_l(:).NE.-1)) input%l_useapw = .TRUE.
     DO iType = 1, atoms%ntype
        IF (atoms%nlo(iType).GE.1) THEN
           IF (input%secvar) THEN
              CALL juDFT_error("LO + sevcar not implemented",calledby ="postprocessInput")
           END IF
           IF (input%isec1<input%itmax) THEN
              CALL juDFT_error("LO + Wu not implemented",calledby ="postprocessInput")
           END IF
           IF (atoms%nlo(iType).GT.atoms%nlod) THEN
              WRITE (6,*) 'nlo(n) =',atoms%nlo(iType),' > nlod =',atoms%nlod
              CALL juDFT_error("nlo(n)>nlod",calledby ="postprocessInput")
           END IF
           DO j=1,atoms%nlo(iType)
              IF (.NOT.input%l_useapw) THEN
                 IF (atoms%llo(j,iType).LT.0) THEN ! CALL juDFT_error("llo<0 ; compile with DCPP_APW!",calledby="inped")
                    WRITE(6,'(A)') 'Info: l_useapw not set.'
                    WRITE(6,'(A,I2,A,I2,A)') '      LO #', j, ' at atom type', iType, ' is an e-derivative.'
                 END IF
              ENDIF
              IF ( (atoms%llo(j,iType).GT.atoms%llod).OR.(mod(-atoms%llod,10)-1).GT.atoms%llod ) THEN
                 WRITE (6,*) 'llo(j,n) =',atoms%llo(j,iType),' > llod =',atoms%llod
                 CALL juDFT_error("llo(j,n)>llod",calledby ="postprocessInput")
              END IF
           END DO

           ! Replace call to setlomap with the following 3 loops (preliminary).
           ! atoms%nlol and atoms%lo1l arrays are strange. This should be solved differently.
           DO l = 0,atoms%llod
              atoms%nlol(l,iType) = 0
              atoms%lo1l(l,iType) = 0
           END DO

           DO ilo = 1,atoms%nlod
              atoms%l_dulo(ilo,iType) = .FALSE.
           END DO

           DO ilo = 1,atoms%nlo(iType)
              if (input%l_useapW) THEN
                 IF (atoms%ulo_der(ilo,iType).EQ.1) THEN
                    atoms%l_dulo(ilo,iType) = .TRUE.
                 END IF
              endif
              WRITE(6,'(A,I2,A,I2)') 'I use',atoms%ulo_der(ilo,iType),'. derivative of l =',atoms%llo(ilo,iType)
              IF (atoms%llo(ilo,iType)>atoms%llod) CALL juDFT_error(" l > llod!!!",calledby="postprocessInput")
              l = atoms%llo(ilo,iType)
              IF (ilo.EQ.1) THEN
                 atoms%lo1l(l,iType) = ilo
              ELSE
                 IF (l.NE.atoms%llo(ilo-1,iType)) THEN
                    atoms%lo1l(l,iType) = ilo
                 END IF
              END IF
              atoms%nlol(l,iType) = atoms%nlol(l,iType) + 1
           END DO
           WRITE (6,*) 'atoms%lapw_l(n) = ',atoms%lapw_l(iType)
        END IF

     END DO

     ! Check lda+u stuff (from inped)

     DO i = 1, atoms%n_u
        n = atoms%lda_u(i)%atomType
        IF (atoms%nlo(n).GE.1) THEN
           DO j = 1, atoms%nlo(n)
              IF ((ABS(atoms%llo(j,n)).EQ.atoms%lda_u(i)%l) .AND. (.NOT.atoms%l_dulo(j,n)) ) &
                 WRITE (*,*) 'LO and LDA+U for same l not implemented'
           END DO
        END IF
     END DO

     IF (atoms%n_u.GT.0) THEN
        IF (input%secvar) CALL juDFT_error("LDA+U and sevcar not implemented",calledby ="postprocessInput")
        IF (input%isec1<input%itmax) CALL juDFT_error("LDA+U and Wu not implemented",calledby ="postprocessInput")
        IF (noco%l_mperp) CALL juDFT_error("LDA+U and l_mperp not implemented",calledby ="postprocessInput")
     END IF

     ! Check DOS related stuff (from inped)

     IF ((banddos%ndir.LT.0).AND..NOT.banddos%dos) THEN
        CALL juDFT_error('STOP banddos: the inbuild dos-program  <0'//&
             ' can only be used if dos = true',calledby ="postprocessInput")
     END IF

     IF ((banddos%ndir.LT.0).AND.banddos%dos) THEN
        IF (banddos%e1_dos-banddos%e2_dos.LT.1e-3) THEN
           CALL juDFT_error("STOP banddos: no valid energy window for "//&
                "internal dos-program",calledby ="postprocessInput")
        END IF
        IF (banddos%sig_dos.LT.0) THEN
           CALL juDFT_error("STOP DOS: no valid broadening (sig_dos) for "//&
                "internal dos-PROGRAM",calledby ="postprocessInput")
        END IF
     END IF

     IF (banddos%vacdos) THEN
        IF (.NOT.banddos%dos) THEN
           CALL juDFT_error("STOP DOS: only set vacdos = .true. if dos = .true.",calledby ="postprocessInput")
        END IF
        IF (.NOT.vacuum%starcoeff.AND.(vacuum%nstars.NE.1))THEN
           CALL juDFT_error("STOP banddos: if stars = f set vacuum=1",calledby ="postprocessInput")
        END IF
        IF (vacuum%layers.LT.1) THEN
           CALL juDFT_error("STOP DOS: specify layers if vacdos = true",calledby ="postprocessInput")
        END IF
        DO i=1,vacuum%layers
           IF (vacuum%izlay(i,1).LT.1) THEN
              CALL juDFT_error("STOP DOS: all layers must be at z>0",calledby ="postprocessInput")
           END IF
        END DO
     END IF

     ! Check noco stuff and calculate missing noco parameters

     IF (noco%l_noco) THEN
        CALL nocoInputCheck(atoms,input,vacuum,noco)

        IF (noco%l_ss) THEN

           !--->    the angle beta is relative to the spiral in a spin-spiral
           !--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
           !--->    that means that the moments are "in line" with the spin-spiral
           !--->    (beta = qss * taual). note: this means that only atoms within
           !--->    a plane perpendicular to qss can be equivalent!

           na = 1
           DO iType = 1,atoms%ntype
              noco%phi = tpi_const*dot_product(noco%qss,atoms%taual(:,na))
              noco%alph(iType) = noco%alphInit(iType) + noco%phi
              na = na + atoms%neq(iType)
           END DO
        END IF
     ELSE
        IF (noco%l_ss) THEN
           CALL judft_warn("l_noco=F and l_ss=T is meaningless. Setting l_ss to F.")
           noco%l_ss = .FALSE.
        END IF
     END IF

     ! Calculate missing kpts parameters
     CALL kpoints(oneD,sym,cell,input,noco,banddos,kpts,l_kpts)
    
     ! Generate missing general parameters
     
     minNeigd = MAX(5,NINT(0.75*input%zelec) + 1)
     IF (noco%l_soc.and.(.not.noco%l_noco)) minNeigd = 2 * minNeigd
     IF (noco%l_soc.and.noco%l_ss) minNeigd=(3*minNeigd)/2
     IF (dimension%neigd.LT.minNeigd) THEN
        IF (dimension%neigd>0) THEN
           WRITE(*,*) 'numbands is too small. Setting parameter to default value.'
           WRITE(*,*) 'changed numbands (dimension%neigd) to ',minNeigd
        ENDIF
        dimension%neigd = minNeigd
     END IF

   
     !cell%aamat=matmul(transpose(cell%amat),cell%amat)
     cell%bbmat=matmul(cell%bmat,transpose(cell%bmat))

     CALL lapw_dim(kpts,cell,input,noco,oneD,forcetheo,DIMENSION)

     CALL lapw_fft_dim(cell,input,noco,stars)
     
        
     obsolete%lepr = 0

     IF (noco%l_noco) dimension%neigd = 2*dimension%neigd

     ! Generate missing parameters for atoms and calculate volume of the different regions

     cell%volint = cell%vol
     atoms%jmtd = maxval(atoms%jri(:))
     CALL ylmnorm_init(atoms%lmaxd)
     dimension%nspd=(atoms%lmaxd+1+mod(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)
     rmtmax = maxval(atoms%rmt(:))
     rmtmax = rmtmax*stars%gmax
     CALL convn_dim(rmtmax,dimension%ncvd)
     dimension%msh = 0
     ALLOCATE(atoms%rmsh(atoms%jmtd,atoms%ntype))
     ALLOCATE(atoms%volmts(atoms%ntype))
     na = 0
     DO iType = 1, atoms%ntype
        l_vca = .FALSE.
        INQUIRE (file="vca.in", exist=l_vca)
        IF (l_vca) THEN
           WRITE(*,*) 'Note: Implementation for virtual crystal approximation should be changed in postprocessInput!'
           WRITE(*,*) 'I am not sure whether the implementation actually makes any sense. It is from inped.'
           WRITE(*,*) 'We have to get rid of the file vca.in!'
           OPEN (17,file='vca.in',form='formatted')
           DO i= 1, iType
              READ (17,*,IOSTAT=ios) ntst,zp
              IF (ios /= 0) EXIT
              IF (ntst == iType) THEN
                 atoms%zatom(iType) = atoms%zatom(iType) + zp
              END IF
           END DO
           CLOSE (17)
        END IF

        ! Calculate mesh for valence states
        radius = atoms%rmt(iType)*exp(atoms%dx(iType)*(1-atoms%jri(iType)))
        dr = exp(atoms%dx(iType))
        DO i = 1, atoms%jri(iType)
           atoms%rmsh(i,iType) = radius
           radius = radius*dr
        END DO
        ! Calculate mesh dimension for core states
        radius = atoms%rmt(iType)
        jrc = atoms%jri(iType)
        DO WHILE (radius < atoms%rmt(iType) + 20.0)
           jrc = jrc + 1
           radius = radius*dr
        END DO
        dimension%msh = max(dimension%msh,jrc)

        atoms%volmts(iType) = (fpi_const/3.0)*atoms%rmt(iType)**3
        cell%volint = cell%volint - atoms%volmts(iType)*atoms%neq(iType)
     END DO


     ! Check muffin tin radii

     ALLOCATE (jri1(atoms%ntype), lmax1(atoms%ntype))
     ALLOCATE (rmt1(atoms%ntype), dx1(atoms%ntype))
     l_test = .TRUE. ! only checking, dont use new parameters
     l_gga=xcpot%needs_grad()
     CALL chkmt(atoms,input,vacuum,cell,oneD,l_gga,noel,l_test,&
                kmax1,dtild1,dvac1,lmax1,jri1,rmt1,dx1)
     DEALLOCATE (jri1,lmax1,rmt1,dx1)

     ! Dimensioning of lattice harmonics

     ALLOCATE(atoms%nlhtyp(atoms%ntype),atoms%ntypsy(atoms%nat))
     ALLOCATE(sphhar%clnu(1,1,1),sphhar%nlh(1),sphhar%llh(1,1),sphhar%nmem(1,1),sphhar%mlh(1,1,1))
     sphhar%ntypsd = 0
     IF (.NOT.oneD%odd%d1) THEN
        CALL local_sym(atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
             atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat,&
             atoms%taual,sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.true.,&
             atoms%nlhtyp,atoms%ntypsy,sphhar%nlh,sphhar%llh,&
             sphhar%nmem,sphhar%mlh,sphhar%clnu)
     ELSE IF (oneD%odd%d1) THEN
        WRITE(*,*) 'Note: I would be surprised if lattice harmonics generation works'
        WRITE(*,*) 'Dimensioning of local arrays seems to be inconsistent with routine local_sym'
        ALLOCATE (nq1(atoms%nat),lmx1(atoms%nat),nlhtp1(atoms%nat))
        ii = 1
        nq1=1
        DO i = 1,atoms%ntype
           DO j = 1,atoms%neq(i)
              lmx1(ii) = atoms%lmax(i)
              ii = ii + 1
           END DO
        END DO
        CALL local_sym(atoms%lmaxd,lmx1,sym%nop,sym%mrot,sym%tau,&
                       atoms%nat,atoms%nat,nq1,cell%amat,cell%bmat,atoms%taual,&
                       sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.true.,nlhtp1,&
                       atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem,&
                       sphhar%mlh,sphhar%clnu)        
        ii = 1
        DO i = 1,atoms%ntype
           atoms%nlhtyp(i) = nlhtp1(ii)
           ii = ii + atoms%neq(i)
        END DO
        DEALLOCATE (nq1,lmx1,nlhtp1)
     END IF
     DEALLOCATE(sphhar%clnu,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh)

     ALLOCATE(sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
     ALLOCATE(sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd))
     ALLOCATE(sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
     ALLOCATE(sphhar%nlh(sphhar%ntypsd),sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd))

     ! Dimensioning of stars

     IF (input%film.OR.(sym%namgrp.ne.'any ')) THEN
        CALL strgn1_dim(stars%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
                        sym%tau,sym%nop,sym%nop2,stars%mx1,stars%mx2,stars%mx3,&
                        stars%ng3,stars%ng2,oneD%odd)

     ELSE
        CALL strgn2_dim(stars%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
                        sym%tau,sym%nop,stars%mx1,stars%mx2,stars%mx3,&
                        stars%ng3,stars%ng2)
        oneD%odd%n2d = stars%ng2
        oneD%odd%nq2 = stars%ng2
        oneD%odd%nop = sym%nop
     END IF

     stars%kimax2= (2*stars%mx1+1)* (2*stars%mx2+1)-1
     stars%kimax = (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1)-1
     IF (oneD%odd%d1) THEN
        oneD%odd%k3 = stars%mx3
        oneD%odd%nn2d = (2*(oneD%odd%k3)+1)*(2*(oneD%odd%M)+1)
     ELSE
        oneD%odd%k3 = 0
        oneD%odd%M = 0
        oneD%odd%nn2d = 1
        oneD%odd%mb = 0
     END IF
     ALLOCATE (stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
     ALLOCATE (stars%ig2(stars%ng3))
     ALLOCATE (stars%kv2(2,stars%ng2),stars%kv3(3,stars%ng3))
     ALLOCATE (stars%nstr2(stars%ng2),stars%nstr(stars%ng3))
     ALLOCATE (stars%sk2(stars%ng2),stars%sk3(stars%ng3),stars%phi2(stars%ng2))
     ALLOCATE (stars%igfft(0:stars%kimax,2),stars%igfft2(0:stars%kimax2,2))
     ALLOCATE (stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
     ALLOCATE (stars%pgfft(0:stars%kimax),stars%pgfft2(0:stars%kimax2))
     ALLOCATE (stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1),stars%ustep(stars%ng3))

     stars%sk2(:) = 0.0
     stars%phi2(:) = 0.0

     ! Initialize xc fft box

     CALL prp_xcfft_box(xcpot%gmaxxc,cell%bmat,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft)

     ! Initialize missing 1D code arrays

     ALLOCATE (oneD%ig1(-oneD%odd%k3:oneD%odd%k3,-oneD%odd%M:oneD%odd%M))
     ALLOCATE (oneD%kv1(2,oneD%odd%n2d),oneD%nstr1(oneD%odd%n2d))
     ALLOCATE (oneD%ngopr1(atoms%nat),oneD%mrot1(3,3,oneD%odd%nop),oneD%tau1(3,oneD%odd%nop))
     ALLOCATE (oneD%invtab1(oneD%odd%nop),oneD%multab1(oneD%odd%nop,oneD%odd%nop))
     ALLOCATE (oneD%igfft1(0:oneD%odd%nn2d-1,2),oneD%pgfft1(0:oneD%odd%nn2d-1))

     ! Initialize missing hybrid functionals arrays

     ALLOCATE (hybrid%nindx(0:atoms%lmaxd,atoms%ntype))
   
     ! Generate lattice harmonics

     IF (.NOT.oneD%odd%d1) THEN
        CALL local_sym(atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
                       atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat,atoms%taual,&
                       sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.,&
                       atoms%nlhtyp,atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu)
        sym%nsymt = sphhar%ntypsd
        oneD%mrot1(:,:,:) = sym%mrot(:,:,:)
        oneD%tau1(:,:) = sym%tau(:,:)
     ELSE IF (oneD%odd%d1) THEN
        WRITE(*,*) 'Note: I would be surprised if lattice harmonics generation works'
        WRITE(*,*) 'Dimensioning of local arrays seems to be inconsistent with routine local_sym'
        CALL od_chisym(oneD%odd,oneD%mrot1,oneD%tau1,sym%zrfs,sym%invs,sym%invs2,cell%amat)
        ALLOCATE (nq1(atoms%nat),lmx1(atoms%nat),nlhtp1(atoms%nat))
        ii = 1
        DO i = 1,atoms%ntype
           DO j = 1,atoms%neq(i)
              nq1(ii) = 1
              lmx1(ii) = atoms%lmax(i)
              ii = ii + 1
           END DO
        END DO
        CALL local_sym(atoms%lmaxd,lmx1,sym%nop,sym%mrot,sym%tau,&
                       atoms%nat,atoms%nat,nq1,cell%amat,cell%bmat,atoms%taual,&
                       sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.,&
                       nlhtp1,atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu)
        sym%nsymt = sphhar%ntypsd
        ii = 1
        DO i = 1,atoms%ntype
           atoms%nlhtyp(i) = nlhtp1(ii)
           ii = ii + atoms%neq(i)
        END DO
        DEALLOCATE (lmx1,nlhtp1)
     END IF

     ! Calculate additional symmetry information

     IF (atoms%n_u.GT.0) THEN
        CALL d_wigner(sym%nop,sym%mrot,cell%bmat,3,sym%d_wgn)
     END IF
     IF (.NOT.oneD%odd%d1) THEN
        CALL mapatom(sym,atoms,cell,input,noco)
        oneD%ngopr1(1:atoms%nat) = atoms%ngopr(1:atoms%nat)
        !     DEALLOCATE ( nq1 )
     ELSE
        CALL juDFT_error("The oneD version is broken here. Compare call to mapatom with old version")
        CALL mapatom(sym,atoms,cell,input,noco)
        CALL od_mapatom(oneD,atoms,sym,cell)
     END IF

     !Calculate kpoint in the full BZ
     IF (kpts%l_gamma.and. banddos%ndir .eq. 0.and.kpts%specificationType==2) THEN
        CALL gen_bz(kpts,sym)
     ELSE
        kpts%nkptf=0
     ENDIF

     ! Missing xc functionals initializations
     IF (xcpot%needs_grad()) THEN
        ALLOCATE (stars%ft2_gfx(0:stars%kimax2),stars%ft2_gfy(0:stars%kimax2))
        ALLOCATE (oneD%pgft1x(0:oneD%odd%nn2d-1),oneD%pgft1xx(0:oneD%odd%nn2d-1),&
                  oneD%pgft1xy(0:oneD%odd%nn2d-1),&
                  oneD%pgft1y(0:oneD%odd%nn2d-1),oneD%pgft1yy(0:oneD%odd%nn2d-1))
     ELSE
        ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
        ALLOCATE (oneD%pgft1x(0:1),oneD%pgft1xx(0:1),oneD%pgft1xy(0:1),&
                  oneD%pgft1y(0:1),oneD%pgft1yy(0:1))
     END IF
     oneD%odd%nq2 = oneD%odd%n2d
     oneD%odi%nq2 = oneD%odd%nq2

     ! Store structure data

     CALL storeStructureIfNew(input, atoms, cell, vacuum, oneD, sym)

     ! Generate stars

     IF (input%film.OR.(sym%namgrp.NE.'any ')) THEN
        CALL strgn1(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
        IF (oneD%odd%d1) THEN
           CALL od_strgn1(xcpot,cell,sym,oneD)
        END IF
     ELSE
        CALL strgn2(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
     END IF

     ! Other small stuff

     input%strho = .FALSE.

     INQUIRE(file="cdn1",exist=l_opti)
     if (noco%l_noco) INQUIRE(file="rhomat_inp",exist=l_opti)
     l_opti=.not.l_opti
     IF ((sliceplot%iplot).OR.(input%strho).OR.(input%swsp).OR.&
         (input%lflip).OR.(input%l_bmt)) l_opti = .TRUE.

     IF (.NOT.l_opti) THEN
        !      The following call to inpeig should not be required.
        !      CALL inpeig(atoms,cell,input,oneD%odd%d1,kpts,enpara)
     END IF

     CALL prp_qfft(stars,cell,noco,input)

     IF (input%gw.GE.1) THEN
        CALL write_gw(atoms%ntype,sym%nop,1,input%jspins,atoms%nat,&
                      atoms%ncst,atoms%neq,atoms%lmax,sym%mrot,cell%amat,cell%bmat,input%rkmax,&
                      atoms%taual,atoms%zatom,cell%vol,1.0,DIMENSION%neigd,atoms%lmaxd,&
                      atoms%nlod,atoms%llod,atoms%nlo,atoms%llo,noco%l_soc)
     END IF

     CALL prp_xcfft(stars,input,cell,xcpot)
 
  END IF !(mpi%irank.EQ.0)
#ifdef CPP_MPI
  CALL MPI_BCAST(sliceplot%iplot,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
#endif

  CALL stepf(sym,stars,atoms,oneD,input,cell,vacuum,mpi)
  IF (.NOT.sliceplot%iplot) THEN   
     IF (mpi%irank.EQ.0) THEN
        CALL convn(DIMENSION,atoms,stars)
        CALL e_field(atoms,DIMENSION,stars,sym,vacuum,cell,input,field%efield)
     END IF !(mpi%irank.EQ.0)
  END IF

  ! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of input postprocessing (calculate missing parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE postprocessInput

END MODULE m_postprocessInput
