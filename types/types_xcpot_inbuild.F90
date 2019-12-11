!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_xcpot_inbuild
   !This module contains the xcpot-type used for the in-build xc-implementations
   USE m_types_xcpot_data
   USE m_types_xcpot
   USE m_types_xcpot_inbuild_nofunction
   USE m_judft
   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=4),PARAMETER:: xc_names(20)=[&
                                'l91 ','x-a ','wign','mjw ','hl  ','bh  ','vwn ','pz  ', &
                                'pw91','pbe ','rpbe','Rpbe','wc  ','PBEs', &
                                'pbe0','hse ','vhse','lhse','exx ','hf  ']

   TYPE, EXTENDS(t_xcpot_inbuild_nf):: t_xcpot_inbuild
   CONTAINS
      !overloading t_xcpot:
      PROCEDURE        :: get_vxc             => xcpot_get_vxc
      PROCEDURE        :: get_exc             => xcpot_get_exc
  END TYPE t_xcpot_inbuild
   PUBLIC t_xcpot_inbuild
 CONTAINS


   SUBROUTINE xcpot_get_vxc(xcpot,jspins,rh, vxc,vx, grad)
!
      USE m_xcxal, ONLY : vxcxal
      USE m_xcwgn, ONLY : vxcwgn
      USE m_xcbh,  ONLY : vxcbh
      USE m_xcvwn, ONLY : vxcvwn
      USE m_xcpz,  ONLY : vxcpz
      USE m_vxcl91
      USE m_vxcwb91
      USE m_vxcpw91
      USE m_vxcepbe
      IMPLICIT NONE
!c
!c---> running mode parameters
!c
      CLASS(t_xcpot_inbuild),INTENT(IN) :: xcpot
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
!c
!c ---> local scalars
      INTEGER :: ngrid
      REAL, PARAMETER :: hrtr_half = 0.5

      !used to be dummy arguments for testing
      INTEGER,PARAMETER   :: idsprs=0,isprsv=0,iofile=6
      REAL,PARAMETER      :: sprsv=0.0
      LOGICAL,PARAMETER   :: lwbc=.false. ! l-white-bird-current (ta)
!c
!c.....------------------------------------------------------------------
!c
!c-----> determine exchange correlation potential
!c
      vx (:,:) = 0.0
      vxc(:,:) = 0.0
      ngrid=SIZE(rh,1)

      IF (xcpot%needs_grad()) THEN
         IF (.NOT.PRESENT(grad)) CALL judft_error("Bug: You called get_vxc for a GGA potential without providing derivatives")
         IF (xcpot%is_name("l91")) THEN    ! local pw91
            CALL vxcl91(jspins,ngrid,ngrid,rh,grad%agrt(:ngrid),grad%agru(:ngrid),grad%agrd(:ngrid), grad%g2rt(:ngrid),&
                 grad%g2ru(:ngrid),grad%g2rd(:ngrid),grad%gggrt(:ngrid),grad%gggru(:ngrid),grad%gggrd(:ngrid),&
                 grad%gzgr(:ngrid), vx(:ngrid,:),vxc(:ngrid,:), isprsv,sprsv)
         ELSEIF (xcpot%is_name("pw91")) THEN  ! pw91
            IF (lwbc) THEN
               CALL vxcwb91(jspins,ngrid,ngrid,rh(:ngrid,:),grad%agrt(:ngrid),grad%agru(:ngrid),grad%agrd(:ngrid),&
                 grad%g2rt(:ngrid),grad%g2ru(:ngrid),grad%g2rd(:ngrid),grad%gggrt(:ngrid),grad%gggru(:ngrid),&
                 grad%gggrd(:ngrid),grad%gzgr(:ngrid), vx(:ngrid,:),vxc(:ngrid,:), idsprs,isprsv,sprsv)
            ELSE

               CALL vxcpw91(jspins,ngrid,ngrid,rh(:ngrid,:),grad%agrt(:ngrid),grad%agru(:ngrid),grad%agrd(:ngrid),&
                 grad%g2rt(:ngrid),grad%g2ru(:ngrid),grad%g2rd(:ngrid),grad%gggrt(:ngrid),grad%gggru(:ngrid),&
                 grad%gggrd,grad%gzgr, vx(:ngrid,:),vxc(:ngrid,:), idsprs,isprsv,sprsv)

            ENDIF
         ELSE  ! pbe or similar
            CALL vxcepbe(xcpot%DATA,jspins,ngrid,ngrid,rh(:ngrid,:), grad%agrt,grad%agru,grad%agrd,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd, vx(:ngrid,:),vxc(:ngrid,:))
         ENDIF
      ELSE  !LDA potentials
         IF (xcpot%is_name("x-a"))  THEN   ! X-alpha method
            CALL vxcxal(xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
         ELSEIF (xcpot%is_name("wign")) THEN    ! Wigner interpolation formula
            CALL vxcwgn(xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
         ELSEIF (xcpot%is_name("mjw").OR.xcpot%is_name("bh")) THEN ! von Barth,Hedin correlation
            CALL vxcbh(iofile,xcpot%data,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))

         ELSEIF (xcpot%is_name("vwn")) THEN     ! Vosko,Wilk,Nusair correlation
            CALL vxcvwn(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
         ELSEIF (xcpot%is_name("pz")) THEN     ! Perdew,Zunger correlation
            CALL vxcpz(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
         ELSEIF (xcpot%is_name("hf")) THEN
            ! Hartree-Fock  calculation: X-alpha potential is added to generate a rational local potential,
            !                            later it is subtracted again
!            CALL juDFT_error('HF should now be treated as a GGA functional', calledby='xcpot_get_vxc')
            CALL vxcxal(xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
            !         vxc=0
         ELSEIF (xcpot%is_name("exx")) THEN
            ! if exact exchange calculation do nothing
            vxc = 0
         ELSE
            CALL juDFT_error("Unkown LDA potential",calledby="type xcpot")
         ENDIF
      ENDIF
!
!-----> hartree units
!
      vx  = hrtr_half*vx
      vxc = hrtr_half*vxc

   END SUBROUTINE xcpot_get_vxc

!***********************************************************************
   SUBROUTINE xcpot_get_exc(xcpot,jspins,rh,exc,grad,kinEnergyDen_KS, mt_call)
!***********************************************************************
      USE m_xcxal, ONLY : excxal
      USE m_xcwgn, ONLY : excwgn
      USE m_xcbh,  ONLY : excbh
      USE m_xcvwn, ONLY : excvwn
      USE m_xcpz,  ONLY : excpz
      USE m_excl91
      USE m_excwb91
      USE m_excpw91
      USE m_excepbe
      IMPLICIT NONE

      CLASS(t_xcpot_inbuild),INTENT(IN)     :: xcpot
      INTEGER, INTENT (IN)                  :: jspins
      REAL,INTENT (IN)                      :: rh(:,:)
      REAL, INTENT (OUT)                    :: exc(:)
      TYPE(t_gradients),OPTIONAL,INTENT(IN) ::grad
      LOGICAL, OPTIONAL, INTENT(IN)         :: mt_call
      REAL, INTENT(IN), OPTIONAL            :: kinEnergyDen_KS(:,:)

!c
!c ---> local scalars
      INTEGER :: ngrid
      REAL, PARAMETER :: hrtr_half = 0.5

      !used to be dummy arguments for testing
      INTEGER,PARAMETER   :: idsprs=0,isprsv=0,iofile=6
      REAL,PARAMETER      :: sprsv=0.0
      LOGICAL,PARAMETER   :: lwbc=.false. ! l-white-bird-current (ta)
!c
!c-----> determine exchange correlation energy density
!c
      exc(:) = 0.0
      ngrid=SIZE(rh,1)
      IF (xcpot%exc_is_gga()) THEN
         IF (.NOT.PRESENT(grad)) CALL judft_error("Bug: You called get_exc for a GGA potential without providing derivatives")
         IF (xcpot%is_name("l91")) THEN  ! local pw91
            CALL excl91(jspins,ngrid,ngrid,rh(:ngrid,:),grad%agrt,grad%agru,grad%agrd,grad%g2rt,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd,grad%gzgr, exc, isprsv,sprsv)
         ELSEIF (xcpot%is_name("pw91")) THEN     ! pw91
            IF (lwbc) THEN
               CALL excwb91(ngrid,ngrid,rh(:ngrid,1),rh(:ngrid,2),grad%agrt,grad%agru,grad%agrd, grad%g2rt,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd,grad%gzgr, exc, idsprs,isprsv,sprsv)
            ELSE
               CALL excpw91(jspins,ngrid,ngrid,rh(:ngrid,:),grad%agrt,grad%agru,grad%agrd, grad%g2rt,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd,grad%gzgr, exc, idsprs,isprsv,sprsv)
            ENDIF
         ELSE
            CALL excepbe(xcpot%data,jspins,ngrid,ngrid, rh(:ngrid,:),grad%agrt,grad%agru,grad%agrd,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd, exc)
         ENDIF
      ELSE !LDA
         IF (xcpot%is_name("x-a"))  THEN   ! X-alpha method
            CALL excxal(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("wign")) THEN    ! Wigner interpolation formula
            CALL excwgn(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("mjw").OR.xcpot%is_name("bh")) THEN ! von Barth,Hedin correlation
            CALL excbh(iofile,xcpot%data,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("vwn")) THEN     ! Vosko,Wilk,Nusair correlation
            CALL excvwn(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("pz")) THEN     ! Perdew,Zunger correlation
            CALL excpz(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("hf")) THEN
!            CALL juDFT_error('HF should now be treated as a GGA functional', calledby='xcpot_get_exc')
!            exc=0
            CALL excxal(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("exx")) THEN
            CALL juDFT_error('EXX should now be treated as a GGA functional', calledby='xcpot_get_exc')
         ELSE
            CALL juDFT_error("Unkown LDA potential",calledby="type xcpot")
         ENDIF
      ENDIF
!c-----> hartree units
      exc= hrtr_half*exc

   END SUBROUTINE xcpot_get_exc

 END MODULE m_types_xcpot_inbuild
