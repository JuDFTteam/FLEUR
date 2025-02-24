!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_metagga
   !USE m_cdnval
   PUBLIC  :: calc_EnergyDen
   PRIVATE :: calc_EnergyDen_auxillary_weights, &
              calc_kinEnergyDen_pw, &
              calc_kinEnergyDen_mt

   type t_RS_potden
      REAL, ALLOCATABLE :: is(:,:), mt(:,:)
   end type t_RS_potden

   TYPE t_kinED
      logical             :: set=.false.
      real, allocatable   :: is(:,:)   ! (nsp*jmtd, jspins)
      real, allocatable   :: mt(:,:,:) ! (nsp*jmtd, jspins, local num of types)
   contains
      procedure           :: alloc_mt => kED_alloc_mt
   END TYPE t_kinED

CONTAINS

  SUBROUTINE kED_alloc_mt(kED,nsp_x_jmtd, jspins, n_start, n_types, n_stride)
    IMPLICIT NONE
    class(t_kinED), intent(inout)   :: kED
    integer, intent(in)            :: nsp_x_jmtd, jspins, n_start, n_types, n_stride
    integer                        :: cnt, n

    if(.not. allocated(kED%mt)) then
      cnt = 0
      do n = n_start,n_types,n_stride
        cnt = cnt + 1
      enddo
      allocate(kED%mt(nsp_x_jmtd, jspins, cnt), source=0.0)
    endif
  END SUBROUTINE kED_alloc_mt


   SUBROUTINE calc_kinEnergyDen_pw(EnergyDen_rs, vTot_rs, den_rs, kinEnergyDen_RS)
      USE m_juDFT_stop
      !use m_cdngen
      IMPLICIT NONE
      REAL, INTENT(in)                 :: den_RS(:,:), EnergyDen_RS(:,:), vTot_RS(:,:)
      REAL, INTENT(inout), allocatable :: kinEnergyDen_RS(:,:)
#ifdef CPP_LIBXC
      REAL, PARAMETER                  :: eps = 1e-15

      kinEnergyDen_RS = EnergyDen_RS - vTot_RS * den_RS
#else
      CALL juDFT_error("MetaGGA require LibXC",hint="compile Fleur with LibXC (e.g. by giving '-external libxc' to ./configure")
#endif
   END SUBROUTINE calc_kinEnergyDen_pw

   SUBROUTINE calc_kinEnergyDen_mt(EnergyDen_RS, vTot_rs, vTot0_rs, core_den_rs, val_den_rs, &
                                   kinEnergyDen_RS)
      USE m_juDFT_stop
      implicit none
      REAL, INTENT(in)                 :: EnergyDen_RS(:,:), vTot_rs(:,:), vTot0_rs(:,:), core_den_rs(:,:), val_den_rs(:,:)
      REAL, INTENT(inout)              :: kinEnergyDen_RS(:,:)

#ifdef CPP_LIBXC
      ! Original 
      !kinEnergyDen_RS = EnergyDen_RS - (vTot0_rs * core_den_rs + vTot_rs * val_den_rs)
      kinEnergyDen_RS = EnergyDen_RS - ( vTot_rs * val_den_rs)
#else
      CALL juDFT_error("MetaGGA require LibXC",hint="compile Fleur with LibXC (e.g. by giving '-external libxc' to ./configure")
#endif
   END SUBROUTINE calc_kinEnergyDen_mt


   SUBROUTINE calc_EnergyDen(eig_id, fmpi, kpts, noco,nococonv, input, banddos, cell, atoms, enpara, stars, &
         vacuum,  sphhar, sym, gfinp, hub1inp, vTot,   results, EnergyDen,xcpot)
      ! calculates the energy density
      ! EnergyDen = \sum_i n_i(r) \varepsilon_i
      ! where n_i(r) is the one-particle density
      ! and \varepsilon_i are the eigenenergies

      USE m_types_setup
      USE m_types_potden
      USE m_types_kpts
      USE m_types_mpi
      USE m_types_enpara
      USE m_types_misc
      USE m_types_regionCharges
      USE m_types_dos
      USE m_types_vacdos
      USE m_types_cdnval
      use m_types_nococonv
      USE m_cdnval
      IMPLICIT NONE

      INTEGER,           INTENT(in)           :: eig_id
      TYPE(t_mpi),       INTENT(in)           :: fmpi
      TYPE(t_kpts),      INTENT(in)           :: kpts
      TYPE(t_noco),      INTENT(in)           :: noco
      TYPE(t_nococonv),  INTENT(in)           :: nococonv
      TYPE(t_input),     INTENT(in)           :: input
      TYPE(t_banddos),   INTENT(in)           :: banddos
      TYPE(t_cell),      INTENT(in)           :: cell
      TYPE(t_atoms),     INTENT(in)           :: atoms
      TYPE(t_enpara),    INTENT(in)           :: enpara
      TYPE(t_stars),     INTENT(in)           :: stars
      TYPE(t_vacuum),    INTENT(in)           :: vacuum

      TYPE(t_sphhar),    INTENT(in)           :: sphhar
      TYPE(t_sym),       INTENT(in)           :: sym
      TYPE(t_gfinp),     INTENT(in)           :: gfinp
      TYPE(t_hub1inp),   INTENT(in)           :: hub1inp
      TYPE(t_potden),    INTENT(in)           :: vTot
       
      TYPE(t_results),   INTENT(in)           :: results
      TYPE(t_potden),    INTENT(inout)        :: EnergyDen

      ! local
      INTEGER                         :: i,jspin

      TYPE(t_regionCharges)           :: regCharges
      TYPE(t_dos)                     :: dos
      TYPE(t_vacdos)                  :: vacdos

      TYPE(t_moments)                 :: moments
      TYPE(t_results)                 :: tmp_results
      TYPE(t_cdnvalJob)               :: cdnvalJob
      TYPE(t_potden)                  :: aux_den, real_den

      CLASS(t_xcpot),INTENT(IN)     :: xcpot  
      

      CALL regCharges%init(input, atoms)
      CALL dos%init(input,        atoms, kpts, banddos,results%eig)
      CALL vacdos%init(input,        atoms, kpts, banddos,results%eig)
!      CALL moments%init(input,    atoms)
      CALL moments%init(fmpi,input,sphhar,atoms)
      tmp_results = results

      DO jspin = 1,input%jspins
         CALL cdnvalJob%init(fmpi,input,kpts,noco,results,jspin)


         ! replace brillouin weights with auxillary weights
         CALL calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, cdnvalJob%weights)

         !WRITE(11111,*) "METAGGA"
         !WRITE(5432,*) 'Point A'
         CALL cdnval(xcpot,eig_id, fmpi, kpts, jspin, noco,nococonv, input, banddos, cell, atoms, &
            enpara, stars, vacuum,  sphhar, sym, vTot,   cdnvalJob, &
            EnergyDen, regCharges, dos, vacdos,tmp_results, moments, gfinp, hub1inp)
         !WRITE(5432,*) 'Point B'
         !DO i = 1, SIZE(EnergyDen%pw,1)
         !   WRITE(5432,'(2i6,2f15.8)') jspin, i, EnergyDen%pw(i,jspin)
         !END DO
      ENDDO

   END SUBROUTINE calc_EnergyDen

   SUBROUTINE calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, f_ik)
      USE m_types_kpts
      USE m_eig66_io
      IMPLICIT NONE
      ! calculates new (auxillary-)weights as
      ! f_iks = w_iks * E_iks
      !, where  f_iks are the new (auxillary-)weights
      ! w_iks are the weights used in brillouin zone integration
      ! E_iks are the eigen energies

      INTEGER,      INTENT(in)        :: eig_id
      INTEGER,      INTENT(in)        :: jspin
      TYPE(t_kpts), INTENT(in)        :: kpts
      REAL,         INTENT(inout)     :: f_ik(:,:) ! f_ik(band_idx, kpt_idx)

      ! local vars
      REAL                       :: e_i(SIZE(f_ik,dim=1))
      INTEGER                    :: i, ikpt

      DO ikpt = 1,kpts%nkpt
         CALL read_eig(eig_id,ikpt,jspin, eig=e_i)
         f_ik(:,ikpt) = f_ik(:,ikpt) * e_i
         DO i = 1, SIZE(f_ik,1)
            !WRITE(4653,'(2i6,2f15.8)') ikpt, i, e_i(i), f_ik(i,ikpt)
         END DO
      ENDDO
   END SUBROUTINE calc_EnergyDen_auxillary_weights

   subroutine set_zPrime(dim_idx, zMat, kpt, lapw, cell, zPrime)
      USE m_types
      USE m_constants
      implicit none
      INTEGER, intent(in)      :: dim_idx
      TYPE (t_mat), intent(in) :: zMat
      REAL, intent(in)         :: kpt(3)
      TYPE(t_lapw), intent(in) :: lapw
      TYPE(t_cell), intent(in) :: cell
      TYPE (t_mat)             :: zPrime

      REAL                     :: k_plus_g(3), fac
      INTEGER                  :: basis_idx

      call zPrime%free()
      call zPrime%init(zMat)

      do basis_idx = 1,size(lapw%gvec,dim=2)
         k_plus_g = kpt + lapw%gvec(:,basis_idx,1)
         k_plus_g = internal_to_rez(cell, k_plus_g)

         fac = k_plus_g(dim_idx)
         if(zPrime%l_real) then
            zPrime%data_r(basis_idx,:) =            fac * zMat%data_r(basis_idx,:)
         else
            zPrime%data_c(basis_idx,:) = ImagUnit * fac * zMat%data_c(basis_idx,:)
         endif
      enddo
   end subroutine set_zPrime

   function internal_to_rez(cell, vec) result(res)
      use m_types
      implicit none
      type(t_cell), intent(in) :: cell
      real, intent(in)      :: vec(3)
      real                  :: res(3)

      res = matmul(transpose(cell%bmat), vec)
   end function internal_to_rez

   subroutine undo_vgen_finalize(vtot, atoms, noco, stars)
      use m_types
      use m_constants
      use m_judft
      implicit none
      TYPE(t_potden), intent(inout)  :: vtot
      type(t_atoms), intent(in)      :: atoms
      type(t_noco), intent(in)       :: noco
      type(t_stars), intent(in)      :: stars

      integer                        :: js, n, st

      do js = 1,size(vtot%mt,4)
         do n = 1,atoms%ntype
            vTot%mt(:atoms%jri(n),0,n,js) = vtot%mt(:atoms%jri(n),0,n,js) &
                  / (atoms%rmsh(:atoms%jri(n),n) / sfp_const )
         enddo
      enddo

      if(.not. noco%l_noco) then
         do js=1,size(vtot%pw_w,2)
            do st=1,stars%ng3
               vTot%pw_w(st,js) = vTot%pw_w(st,js) * stars%nstr(st)
            enddo
         enddo
      else
         call juDFT_error("undo vgen_finalize not implemented for noco")
      endif
   end subroutine undo_vgen_finalize

   subroutine set_kinED(fmpi,   sphhar, atoms, sym,  xcpot, &
                        input, noco,   stars, vacuum ,cell,     den,     EnergyDen, vTot,kinED)

      USE m_pw_tofrom_grid             
      use m_types
      use m_cdn_io
      USE m_pwden

      implicit none
      TYPE(t_mpi),INTENT(IN)       :: fmpi
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sym), INTENT(IN)      :: sym
      CLASS(t_xcpot),INTENT(IN)    :: xcpot
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
       
      TYPE(t_cell),INTENT(IN)      :: cell
      !TYPE(t_potden),INTENT(IN)    :: den, EnergyDen, vTot
      TYPE(t_kinED),INTENT(OUT)    :: kinED
      TYPE(t_potden)               :: vTot_corrected


!     ------------ Yamashita_add

      TYPE(t_potden)  :: den, EnergyDen, vTot
      INTEGER      :: jspin
      TYPE (t_mat) :: zMat
      REAL  :: kpt(3)
      TYPE (t_force)             :: force  
      TYPE(t_lapw) :: lapw  
      TYPE (t_mat) :: zPrime
      INTEGER      :: eig_id
      TYPE(t_kpts) :: kpts
      TYPE(t_nococonv) :: nococonv
      TYPE(t_banddos)  :: banddos
      TYPE(t_dos)      :: dos
      TYPE(t_enpara)   :: enpara
      REAL,    ALLOCATABLE  :: we(:),eig(:)
      REAL                  :: bkpt(3)
      INTEGER, ALLOCATABLE  :: ev_list(:)
      INTEGER :: ne,ikpt,ikpt_i
      LOGICAL :: l_empty
      TYPE(t_results)      :: results
      TYPE(t_gradients)            :: tmp_grad
      INTEGER :: I,REPLACE
      TYPE(t_cdnvalJob) :: cdnvalJob

!     --------------------- 

      LOGICAL :: perform_MetaGGA
      TYPE(t_potden)    :: core_den, val_den
      real :: rdum, tempDistance
      logical :: ldum
      perform_MetaGGA = ALLOCATED(EnergyDen%mt) &
                      .AND. (xcpot%exc_is_MetaGGA() .or. xcpot%vx_is_MetaGGA())
      if(.not.perform_MetaGGA) return
#ifdef CPP_LIBXC
      core_den=den
      CALL core_den%resetPotDen
      call readDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym ,&
                          CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                           0,rdum,tempDistance,ldum,core_den,inFilename='cdnc')
      call val_den%subPotDen(den,core_den)

      call vTot_corrected%copyPotDen(vTot)
      call undo_vgen_finalize(vTot_corrected, atoms, noco, stars)

      call set_kinED_is(xcpot, input, noco, stars, sym, cell, den, EnergyDen, vTot_corrected,kinED)
      call set_kinED_mt(fmpi,   sphhar,    atoms, sym, noco,core_den, val_den, &
                           xcpot, EnergyDen, input, vTot_corrected,kinED)


!     ========= Yamashita_add
     
!      Replace=0
!      IF(Replace==1)THEN

       !allocate(den%pw(stars%ng3,input%jspins))
!       allocate(kinED%ispw(stars%ng3,input%jspins))
!       allocate(ev_list())


!       ikpt=cdnvalJob%k_list(ikpt_i) 
!       ev_list=cdnvaljob%compact_ev_list(ikpt,l_empty)

!       DO I=1,3

!       CALL set_zPrime(I, zMat, kpt, lapw, cell, zPrime)
!       CALL pwden(stars, kpts, banddos,   input, fmpi, noco, nococonv, cell, atoms, sym, &
!          ikpt, jspin, lapw, ne, ev_list, we, eig, den, results, force%f_b8, zPrime, dos)

!       kinED%ispw=kinED%ispw+den%pw

!       ENDDO

!       CALL init_pw_grid(stars,sym,cell,xcpot) 
!       CALL pw_to_grid(.FALSE., input%jspins, noco%l_noco, stars, &
!                             cell,  kinED%ispw,       tmp_grad, xcpot,   kinED%is) 
!       CALL finish_pw_grid()

!      ENDIF

!     ======== 


  #endif
   end subroutine set_kinED
#ifdef CPP_LIBXC
   subroutine set_kinED_is(xcpot, input, noco, stars, sym, cell, den, EnergyDen, vTot,kinED)
      use m_types
      use m_pw_tofrom_grid
      !USE m_cdnval
      implicit none
      CLASS(t_xcpot),INTENT(IN)    :: xcpot
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_sym), INTENT(IN)      :: sym
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_potden),INTENT(IN)    :: den, EnergyDen, vTot
      TYPE(t_kinED),INTENT(INOUT)  :: kinED

      !local arrays
      REAL, ALLOCATABLE            :: den_rs(:,:), ED_rs(:,:), vTot_rs(:,:)
      TYPE(t_gradients)            :: tmp_grad

      INTEGER :: I 

      CALL init_pw_grid(stars,sym,cell,xcpot)

      CALL pw_to_grid(xcpot%needs_grad(), input%jspins, noco%l_noco, stars, &
                      cell,  EnergyDen%pw, tmp_grad, xcpot,    ED_rs)
      CALL pw_to_grid(xcpot%needs_grad(), input%jspins, noco%l_noco, stars, &
                      cell,  vTot%pw,      tmp_grad, xcpot,    vTot_rs)
      CALL pw_to_grid(xcpot%needs_grad(), input%jspins, noco%l_noco, stars, &
                      cell,  den%pw,       tmp_grad, xcpot,   den_rs)

      CALL finish_pw_grid()

      call calc_kinEnergyDen_pw(ED_rs, vTot_rs, den_rs, kinED%is)
      !xcpot%kinED%is  = ED_RS - vTot_RS * den_RS
      kinED%set = .True.
    
       DO I=1,size(kinED%is,DIM=1)
!        WRITE(111111,"(I5,F15.9)") I,kinED%is(I,1)
!        WRITE(222222,"(I5,F15.9)") I,kinED%is(I,2)
        !WRITE(333,"(I5,5F15.9)") I,ED_rs(I-1,1),vTot_RS(I-1,1),den_RS(I-1,1),kinED%is(I,1),ED_RS(I-1,1)-vTot_RS(I-1,1)*den_RS(I-1,1)
        !WRITE(444,"(I5,5F15.9)") I,ED_rs(I-1,2),vTot_RS(I-1,2),den_RS(I-1,2),kinED%is(I,2),ED_RS(I-1,2)-vTot_RS(I-1,2)*den_RS(I-1,2)
       ENDDO
      

   end subroutine set_kinED_is

   subroutine set_kinED_mt(fmpi,   sphhar,    atoms, sym, noco,core_den, val_den, &
                           xcpot, EnergyDen, input, vTot,kinED)
      use m_types
      use m_mt_tofrom_grid
      implicit none
      TYPE(t_mpi),INTENT(IN)         :: fmpi
      TYPE(t_sphhar),INTENT(IN)      :: sphhar
      TYPE(t_atoms),INTENT(IN)       :: atoms
      TYPE(t_sym), INTENT(IN)        :: sym
      TYPE(t_noco), INTENT(IN)       :: noco
      TYPE(t_potden),INTENT(IN)      :: core_den, val_den, EnergyDen, vTot
      CLASS(t_xcpot),INTENT(IN)      :: xcpot
      TYPE(t_input),INTENT(IN)       :: input
      TYPE(t_kinED),INTENT(INOUT)    :: kinED
      INTEGER                        :: I,jr, loc_n, n, n_start, n_stride, cnt
      REAL,ALLOCATABLE               :: vTot_mt(:,:,:), ED_rs(:,:), vTot_rs(:,:), vTot0_rs(:,:),&
                                        core_den_rs(:,:), val_den_rs(:,:)
      TYPE(t_gradients)              :: tmp_grad
      TYPE(t_sphhar)                 :: tmp_sphhar

#ifdef CPP_MPI
      n_start=fmpi%irank+1
      n_stride=fmpi%isize
#else
      n_start=1
      n_stride=1
#endif
      CALL init_mt_grid(input%jspins,atoms,sphhar,xcpot%needs_grad(),sym)
      loc_n = 0
      allocate(ED_rs(atoms%nsp()*atoms%jmtd, input%jspins))
      allocate(vTot_rs, mold=ED_rs)
      allocate(vTot0_rs, mold=ED_rs)
      allocate(core_den_rs, mold=ED_rs)
      allocate(val_den_rs, mold=ED_rs)

      call kinED%alloc_mt(atoms%nsp()*atoms%jmtd, input%jspins, &
                                n_start,                atoms%ntype,  n_stride)
      loc_n = 0
      do n = n_start,atoms%ntype,n_stride
         loc_n = loc_n + 1

         if(.not. allocated(vTot_mt)) then
            allocate(vTot_mt(lbound(vTot%mt, dim=1):ubound(vTot%mt, dim=1),&
                             lbound(vTot%mt, dim=2):ubound(vTot%mt, dim=2),&
                             lbound(vTot%mt, dim=4):ubound(vTot%mt, dim=4)))
         endif

         do jr=1,atoms%jri(n)
            vTot_mt(jr,0:,:) = vTot%mt(jr,0:,n,:) * atoms%rmsh(jr,n)**2
         enddo
         CALL mt_to_grid(xcpot%needs_grad(), input%jspins, atoms,sym, sphhar,.TRUE., EnergyDen%mt(:, 0:, n, :), &
                         n,  noco,   tmp_grad,     ED_rs)
         CALL mt_to_grid(xcpot%needs_grad(), input%jspins, atoms, sym,sphhar,.TRUE., vTot_mt(:,0:,:), &
                         n,     noco,tmp_grad,     vTot_rs)

         tmp_sphhar%nlhd = sphhar%nlhd
         tmp_sphhar%nlh  = [(0, cnt=1,size(sphhar%nlh))]

         CALL mt_to_grid(xcpot%needs_grad(), input%jspins, atoms, sym,tmp_sphhar,.TRUE., vTot_mt(:,0:0,:), &
                         n,    noco, tmp_grad,     vTot0_rs)
         CALL mt_to_grid(xcpot%needs_grad(), input%jspins, atoms, sym,sphhar,.TRUE., &
                         core_den%mt(:,0:,n,:), n,noco, tmp_grad, core_den_rs)
         CALL mt_to_grid(xcpot%needs_grad(), input%jspins, atoms, sym,sphhar,.TRUE., &
                         val_den%mt(:,0:,n,:), n,noco, tmp_grad, val_den_rs)

         call calc_kinEnergyDen_mt(ED_RS, vTot_rs, vTot0_rs, core_den_rs, val_den_rs, &
                                   kinED%mt(:,:,loc_n))
         !xcpot%kinED%mt(:,:,loc_n) = ED_RS - (vTot0_rs * core_den_rs + vTot_rs * val_den_rs)
          
         enddo

!          DO I=1,size(kinED%mt,DIM=1)
!WRITE(3333,"(I10,2F30.15)") I-1,ED_rs(I-1,1),kinED%mt(I-1,1,loc_n)!vTot_RS(I-1,1),den_RS(I-1,1),kinED%is(I,1),ED_RS(I-1,1)-vTot_RS(I-1,1)*den_RS(I-1,1)
!WRITE(4444,"(I10,2F30.15)") I-1,ED_rs(I-1,2),kinED%mt(I-1,2,loc_n)!vTot_RS(I-1,1),den_RS(I-1,1),kinED%is(I,1),ED_RS(I-1,1)-vTot_RS(I-1,1)*den_RS(I-1,1)
!          ENDDO
   
   
         kinED%set = .True.
      CALL finish_mt_grid()
   end subroutine set_kinED_mt
#endif
END MODULE m_metagga
