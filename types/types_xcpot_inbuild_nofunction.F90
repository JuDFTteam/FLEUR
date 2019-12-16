!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_xcpot_inbuild_nofunction
   !This module contains the xcpot-type used for the in-build xc-implementations
   USE m_types_xcpot_data
   USE m_types_xcpot
   USE m_judft
   IMPLICIT NONE
   PRIVATE
   REAL, PARAMETER, PRIVATE :: hrtr_half = 0.5
   CHARACTER(len=4),PARAMETER:: xc_names(20)=[&
                                'l91 ','x-a ','wign','mjw ','hl  ','bh  ','vwn ','pz  ', &
                                'pw91','pbe ','rpbe','Rpbe','wc  ','PBEs', &
                                'pbe0','hse ','vhse','lhse','exx ','hf  ']

   LOGICAL,PARAMETER:: priv_LDA(20)=[&
                       .FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,&
                       .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                       .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.]


   LOGICAL,PARAMETER:: priv_gga(20)=[&
                       .TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                       .TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,&
                       .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,.TRUE.]

   LOGICAL,PARAMETER:: priv_hybrid(20)=[&
                       .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                       .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                       .TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.]

   REAL, PARAMETER       ::  amix_pbe0 = 0.25
   REAL, PARAMETER       ::  amix_hse  = 0.25
   REAL, PARAMETER       ::  amix_hf   = 1.00

   TYPE, EXTENDS(t_xcpot):: t_xcpot_inbuild_nf
      INTEGER          :: icorr=0
      TYPE(t_xcpot_data) :: data

      LOGICAL,ALLOCATABLE :: lda_atom(:)

   CONTAINS
      !overloading t_xcpot:
      PROCEDURE        :: vx_is_LDA => xcpot_vx_is_LDA
      PROCEDURE        :: vx_is_GGA => xcpot_vx_is_GGA

      PROCEDURE        :: vc_is_LDA => xcpot_vc_is_LDA
      PROCEDURE        :: vc_is_GGA => xcpot_vc_is_GGA

      PROCEDURE        :: exc_is_LDA => xcpot_exc_is_LDA
      PROCEDURE        :: exc_is_gga => xcpot_exc_is_gga
      PROCEDURE        :: is_hybrid  => xcpot_is_hybrid

      PROCEDURE        :: get_exchange_weight => xcpot_get_exchange_weight
      PROCEDURE        :: get_vxc             => xcpot_get_vxc
      PROCEDURE        :: get_exc             => xcpot_get_exc
      !not overloaded
      PROCEDURE        :: get_name => xcpot_get_name
      PROCEDURE        :: relativistic_correction
      PROCEDURE        :: is_name  => xcpot_is_name
      PROCEDURE        :: init     => xcpot_init
      PROCEDURE        :: mpi_bc => mpi_bc_xcpot_ib
   END TYPE t_xcpot_inbuild_nf
   PUBLIC t_xcpot_inbuild_nf
 CONTAINS

   Subroutine Mpi_bc_xcpot_ib(This,Mpi_comm,Irank)
    Use M_mpi_bc_tool
    Class(t_xcpot_inbuild_nf),Intent(Inout)::This
    Integer,Intent(In):: Mpi_comm
    Integer,Intent(In),Optional::Irank
    Integer ::Rank
    If (Present(Irank)) Then
       Rank=Irank
    Else
       Rank=0
    End If

    call this%t_xcpot%mpi_bc(mpi_comm,irank)

    CALL mpi_bc(this%icorr,rank,mpi_comm)
    CALL mpi_bc(this%data%is_rpbe,rank,mpi_comm)
    CALL mpi_bc(this%data%is_wc,rank,mpi_comm)
    CALL mpi_bc(this%data%is_hse,rank,mpi_comm)
    CALL mpi_bc(this%data%uk,rank,mpi_comm)
    CALL mpi_bc(this%data%um,rank,mpi_comm)
    CALL mpi_bc(this%data%is_pbes,rank,mpi_comm)
    CALL mpi_bc(this%data%is_pbe0,rank,mpi_comm)
    CALL mpi_bc(this%data%is_bh,rank,mpi_comm)
    CALL mpi_bc(this%data%is_mjw,rank,mpi_comm)
    CALL mpi_bc(this%data%exchange_weight,rank,mpi_comm)
    CALL mpi_bc(this%data%krla,rank,mpi_comm)
    CALL mpi_bc(this%lda_atom,rank,mpi_comm)



  END SUBROUTINE mpi_bc_xcpot_ib

   LOGICAL FUNCTION relativistic_correction(xcpot)
     IMPLICIT NONE
     CLASS(t_xcpot_inbuild_nf),INTENT(IN)    :: xcpot
     relativistic_correction=xcpot%DATA%krla==1
   END FUNCTION relativistic_correction

   CHARACTER(len=4) FUNCTION xcpot_get_name(xcpot)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild_nf),INTENT(IN)    :: xcpot
      IF (xcpot%icorr==0) CALL judft_error("xc-potential not initialized",calledby="types_xcpot.F90")
      xcpot_get_name=xc_names(xcpot%icorr)
   END FUNCTION xcpot_get_name

   SUBROUTINE xcpot_init(xcpot,ntype)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild_nf),INTENT(INOUT)    :: xcpot
      INTEGER,INTENT(IN)           :: ntype
      INTEGER:: n
      !Determine icorr from name

      IF (.NOT.xcpot%l_inbuild) CALL judft_error("Could not initialize inbuild xcpot")

      ALLOCATE(xcpot%lda_atom(ntype))
      xcpot%lda_atom=.FALSE.
      xcpot%icorr=0
      DO n=1,SIZE(xc_names)
         IF (TRIM(ADJUSTL(xcpot%inbuild_name))==TRIM(xc_names(n))) THEN
            xcpot%icorr=n
         ENDIF
      ENDDO
      if (xcpot%icorr==0) CALL judft_error("Unkown xc-potential:"//xcpot%inbuild_name,calledby="types_xcpot.F90")
      IF (xcpot%l_relativistic)THEN
         xcpot%DATA%krla=1
      ELSE
         xcpot%DATA%krla=0
      END IF

      !Code from exchpbe to speed up determination of constants
      IF (xcpot%is_name("rpbe")) THEN
         xcpot%data%uk=1.2450
      ELSE
         xcpot%data%uk=0.8040
      ENDIF
      IF (xcpot%is_name("PBEs")) THEN     ! pbe_sol
         xcpot%data%um=0.123456790123456d0
      ELSE
         xcpot%data%um=0.2195149727645171e0
      ENDIF
      xcpot%data%is_hse=xcpot%is_name("hse").OR.xcpot%is_name("lhse").OR.xcpot%is_name("vhse")
      xcpot%data%is_rpbe=xcpot%is_name("Rpbe") !Rpbe
      xcpot%data%is_wc=xcpot%is_name("wc")
      xcpot%data%is_pbes=xcpot%is_name("PBEs")
      xcpot%data%is_pbe0=xcpot%is_name("pbe0")
      xcpot%data%is_mjw=xcpot%is_name("mjw")
      xcpot%data%is_bh=xcpot%is_name("bh")
      xcpot%DATA%exchange_weight=xcpot%get_exchange_weight()

   END SUBROUTINE xcpot_init

   !! LDA
   logical function xcpot_exc_is_lda(xcpot)
      implicit none
      CLASS(t_xcpot_inbuild_nf),INTENT(IN):: xcpot
      xcpot_exc_is_lda= xcpot%vxc_is_lda()
   end function xcpot_exc_is_lda

   logical function xcpot_vx_is_lda(xcpot)
      implicit none
      CLASS(t_xcpot_inbuild_nf),INTENT(IN):: xcpot
      xcpot_vx_is_lda=(.not. xcpot%vxc_is_gga()) .and. (.not. xcpot%is_hybrid())
   end function xcpot_vx_is_lda

   logical function xcpot_vc_is_lda(xcpot)
      implicit none
      CLASS(t_xcpot_inbuild_nf),INTENT(IN):: xcpot
      xcpot_vc_is_lda=(.not. xcpot%vxc_is_gga()) .and. (.not. xcpot%is_hybrid())
   end function xcpot_vc_is_lda

   logical function xcpot_vx_is_gga(xcpot)
      implicit none
      CLASS(t_xcpot_inbuild_nf),INTENT(IN):: xcpot
      xcpot_vx_is_gga=priv_gga(xcpot%icorr)
   end function xcpot_vx_is_gga

   logical function xcpot_vc_is_gga(xcpot)
      implicit none
      CLASS(t_xcpot_inbuild_nf),INTENT(IN):: xcpot
      xcpot_vc_is_gga=priv_gga(xcpot%icorr)
   end function xcpot_vc_is_gga

   LOGICAL FUNCTION xcpot_exc_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild_nf),INTENT(IN):: xcpot
      xcpot_exc_is_gga = xcpot%vxc_is_gga()
   END FUNCTION xcpot_exc_is_gga

   LOGICAL FUNCTION xcpot_is_hybrid(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild_nf),INTENT(IN):: xcpot
      xcpot_is_hybrid=priv_hybrid(xcpot%icorr)
   END FUNCTION xcpot_is_hybrid

   FUNCTION xcpot_get_exchange_weight(xcpot) RESULT(a_ex)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild_nf),INTENT(IN):: xcpot

      REAL:: a_ex

      a_ex=-1
      IF (xcpot%is_name("pbe0")) a_ex=amix_pbe0
      IF (xcpot%is_name("hf")) a_ex=amix_hf
      IF (xcpot%is_name("hse")) a_ex=amix_hse
      IF (xcpot%is_name("vhse")) a_ex=amix_hse
   END FUNCTION xcpot_get_exchange_weight

   SUBROUTINE xcpot_get_vxc(xcpot,jspins,rh, vxc,vx, grad)
     !
      IMPLICIT NONE
!c
!c---> running mode parameters
!c
      CLASS(t_xcpot_inbuild_nf),INTENT(IN) :: xcpot
      INTEGER, INTENT (IN)     :: jspins
!c
!c---> charge density
!c
      REAL,INTENT (IN) :: rh(:,:)
!c
!c---> xc potential
!c
      REAL, INTENT (OUT) :: vx (:,:)
      REAL, INTENT (OUT) :: vxc(:,:)

      ! optional arguments for GGA
      TYPE(t_gradients),INTENT(INOUT),OPTIONAL::grad

      CALL judft_error("BUG: dummy xcxpot type is not functional and should not be called")

   END SUBROUTINE xcpot_get_vxc

!***********************************************************************
   SUBROUTINE xcpot_get_exc(xcpot,jspins,rh,exc,grad,kinEnergyDen_KS, mt_call)
!***********************************************************************
      IMPLICIT NONE

      CLASS(t_xcpot_inbuild_nf),INTENT(IN)     :: xcpot
      INTEGER, INTENT (IN)                  :: jspins
      REAL,INTENT (IN)                      :: rh(:,:)
      REAL, INTENT (OUT)                    :: exc(:)
      TYPE(t_gradients),OPTIONAL,INTENT(IN) ::grad
      LOGICAL, OPTIONAL, INTENT(IN)         :: mt_call
      REAL, INTENT(IN), OPTIONAL            :: kinEnergyDen_KS(:,:)

!c
!c ---> local scalars
      CALL judft_error("BUG: dummy xcxpot type is not functional and should not be called")

   END SUBROUTINE xcpot_get_exc

   LOGICAL FUNCTION xcpot_is_name(xcpot,name)
      CLASS(t_xcpot_inbuild_nf),INTENT(IN):: xcpot
      CHARACTER(len=*),INTENT(IN)  :: name
      xcpot_is_name=(TRIM(xc_names(xcpot%icorr))==TRIM((name)))
   END FUNCTION xcpot_is_name

 END MODULE m_types_xcpot_inbuild_nofunction
