!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_stars
  USE m_judft
  USE m_types_fleur_setup
  USE m_json_tools
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_stars
     !max-length of star
     REAL :: gmax
     REAL :: gmaxInit
     !no of 3d-stars
     INTEGER :: ng3
     !no of 2d-stars
     INTEGER :: ng2
     !dim of box
     INTEGER ::mx1
     INTEGER ::mx2
     INTEGER ::mx3
     !No of elements in FFT
     INTEGER ::kimax
     !No of elements in 2D-FFT
     INTEGER ::kimax2

     !Box for FFT in pwden
     INTEGER :: kq1_fft
     INTEGER :: kq2_fft
     INTEGER :: kq3_fft
     INTEGER :: kmxq_fft !no of g-vectors in sphere
     INTEGER, ALLOCATABLE :: igq_fft(:)
     INTEGER, ALLOCATABLE :: igq2_fft(:)

     !fft box for xc-pot
     INTEGER :: kxc1_fft
     INTEGER :: kxc2_fft
     INTEGER :: kxc3_fft

     INTEGER :: ng3_fft
     INTEGER :: kmxxc_fft !<number of g-vectors forming the nxc3_fft stars in the charge density or xc-density sphere

     INTEGER :: nxc3_fft !< number of stars in the  charge density  fft-box
     !rep. g-vector of star
     INTEGER,ALLOCATABLE ::kv3(:,:)
     !length of star
     REAL,ALLOCATABLE    ::sk3(:)
     !mapping of g-vectors to stars
     INTEGER,ALLOCATABLE ::ig(:,:,:)
     !No of g-vectors in star
     INTEGER,ALLOCATABLE ::nstr(:)
     !rep. g-vector of 2D-star
     INTEGER,ALLOCATABLE ::kv2(:,:)
     !length of 2D-star
     REAL,ALLOCATABLE    ::sk2(:)
     !No of g-vecs in 2D-star
     INTEGER,ALLOCATABLE ::nstr2(:)
     !mapping of
     INTEGER,ALLOCATABLE ::ig2(:)
     !
     REAL,ALLOCATABLE:: phi2(:) !<(n2d)
     !phase phactor of g-vector
     COMPLEX,ALLOCATABLE    ::rgphs(:,:,:)
     !mapping of stars to FFT-box
     INTEGER, ALLOCATABLE :: igfft(:,:)
     !same for 2D
     INTEGER, ALLOCATABLE :: igfft2(:,:)
     !phasefactors for mapping
     COMPLEX,ALLOCATABLE  :: pgfft(:)
     !same of 2D
     COMPLEX,ALLOCATABLE  :: pgfft2(:)
     !
     REAL,ALLOCATABLE     :: ft2_gfx(:),ft2_gfy(:)
     COMPLEX, ALLOCATABLE :: ustep(:)
     REAL, ALLOCATABLE    :: ufft(:)
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_stars
     PROCEDURE,PASS :: WRITE=>WRITE_stars
     PROCEDURE,PASS :: READ=>READ_stars
     PROCEDURE,PASS :: read_xml=>read_xml_stars
  END TYPE t_stars

CONTAINS
  SUBROUTINE broadcast_stars(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE m_bc_tools
#endif    
    IMPLICIT NONE
    CLASS(t_stars),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr

    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF

    CALL MPI_BCAST(tt%gmax,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%gmaxinit,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)

    CALL MPI_BCAST(tt%ng3,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%ng2,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%mx1,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%mx2,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%mx3,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kimax,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kimax2,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kq1_fft,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kq2_fft,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kq3_fft,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kmxq_fft ,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kxc1_fft,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kxc2_fft,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kxc3_fft,1,MPI_INTEGER,pe,mpi_comm,ierr)

    CALL MPI_BCAST(tt%ng3_fft,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kmxxc_fft,1,MPI_INTEGER,pe,mpi_comm,ierr)

    CALL MPI_BCAST(tt%nxc3_fft,1,MPI_INTEGER,pe,mpi_comm,ierr)

    CALL MPI_BC(tt%igq_fft,pe,mpi_comm)
    CALL MPI_BC(tt%igq2_fft,pe,mpi_comm)
    !rep. g-vector of star
    CALL MPI_BC(tt%kv3,pe,mpi_comm)
    !length of star
    CALL MPI_BC(tt%sk3,pe,mpi_comm)
    !mapping of g-vectors to stars
    CALL MPI_BC(tt%ig,pe,mpi_comm)
    !No of g-vectors in star
    CALL MPI_BC(tt%nstr,pe,mpi_comm)
    !rep. g-vector of 2D-star
    CALL MPI_BC(tt%kv2,pe,mpi_comm)
    !length of 2D-star
    CALL MPI_BC(tt%sk2,pe,mpi_comm)
    !No of g-vecs in 2D-star
    CALL MPI_BC(tt%nstr2,pe,mpi_comm)
    !mapping of
    CALL MPI_BC(tt%ig2,pe,mpi_comm)
    !
    CALL MPI_BC(tt%phi2,pe,mpi_comm) !<(n2d)
    !phase phactor of g-vector
    CALL MPI_BC(tt%rgphs,pe,mpi_comm)
    !mapping of stars to FFT-box
    CALL MPI_BC(tt%igfft,pe,mpi_comm)
    !same for 2D
    CALL MPI_BC(tt%igfft2,pe,mpi_comm)
    !phasefactors for mapping
    CALL MPI_BC(tt%pgfft,pe,mpi_comm)
    !same of 2D
    CALL MPI_BC(tt%pgfft2,pe,mpi_comm)
    !
    CALL MPI_BC(tt%ft2_gfx,pe,mpi_comm)
    CALL MPI_BC(tt%ft2_gfy,pe,mpi_comm)
    CALL MPI_BC(tt%ustep,pe,mpi_comm)
    CALL MPI_BC(tt%ufft,pe,mpi_comm)


#endif

  END SUBROUTINE broadcast_stars

  SUBROUTINE write_stars(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_stars),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    WRITE(unit,*,IOSTAT=iostat) '"stars":{'


    CALL json_print(unit,"gmax",tt%gmax)
    CALL JSON_PRINT(unit,"gmaxinit",tt%gmaxInit)
    WRITE(unit,*,IOSTAT=iostat) '}'

  END SUBROUTINE write_stars
  SUBROUTINE read_stars(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_stars),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    CHARACTER(len=40)::string
    REAL,ALLOCATABLE:: rtemp(:)
    CALL json_open_class("stars",unit,iostat)
    IF (iostat.NE.0)   RETURN

    CALL json_read(unit,"gmax",tt%gmax)
    CALL json_read(unit,"gmaxInit",tt%gmaxinit)

    CALL json_close_class(unit,iostat)

  END SUBROUTINE read_stars


  SUBROUTINE read_xml_stars(tt)
    USE m_xmlIntWrapFort
    USE m_calculator
    IMPLICIT NONE
    CLASS(t_stars),INTENT(OUT):: tt

    tt%gmax = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/cutoffs/@Gmax'))
    tt%gmaxinit=tt%gmax

  END SUBROUTINE read_xml_stars


  SUBROUTINE init_stars(stars,sym,input)
    USE m_types_sym
    USE m_types_input
    USE m_prpqfftmap
    IMPLICIT NONE
    CLASS(t_stars),INTENT(INOUT):: stars
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_input),INTENT(IN)    :: input
    

    
     CALL lapw_fft_dim(cell,input,noco,stars)

     IF (input%film) THEN
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

      ! Missing xc functionals initializations
     IF (xcpot%is_gga()) THEN
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

       ! Generate stars

     IF (input%film) THEN
        CALL strgn1(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
        IF (oneD%odd%d1) THEN
           CALL od_strgn1(xcpot,cell,sym,oneD)
        END IF
     ELSE
        CALL strgn2(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
     END IF

       CALL prp_qfft(stars,cell,noco,input)

    

       CALL prp_xcfft(stars,input,cell,xcpot)

       CALL stepf(sym,stars,atoms,oneD,input,cell,vacuum,mpi)
           ALLOCATE (stars%igq_fft(0:stars%kq1_fft*stars%kq2_fft*stars%kq3_fft-1))
    ALLOCATE (stars%igq2_fft(0:stars%kq1_fft*stars%kq2_fft-1))

    ! Set up pointer for backtransformation from g-vector in positive 
    ! domain of carge density fftibox into stars
    CALL prp_qfft_map(stars,sym,input)
    
   END SUBROUTINE init_stars
END MODULE m_types_stars
