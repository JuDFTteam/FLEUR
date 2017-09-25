!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_xcallg
      use m_juDFT

c  vxcallg: driver subroutine for the exchange-correlation potentials
c  excallg: driver subroutine for the exchange-correlation energy density.
c
c  inputs a total bare electron density rh(r,1) for jspins=1
c         a bare electron density rh(r,1),rh(r,2) of majority and
c          minority for jspins=2.
c  outputs the effective exch-corr energy density due to this electron
c  density in hartree units.
c
c  depending on the variable xcpot
c  the contribution of exchange and correlation is determined in
c  different subroutines.
c
c
c     krla=1 means : relativistic correction of exchange energy
c                    density related to dirac kinetic energy, according
c                    to: a.h. macdonald and s.h. vosko, j. phys. c12,
c                    2977(1979)
c
c     based on a subroutine from s.bluegel
c     r.pentcheva, 22.01.96

      CONTAINS
c***********************************************************************
      SUBROUTINE vxcallg(xcpot,lwbc,jspins,mfftwk,nfftwk,rh,
     +                  agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,
     +                  gzgr,
     +                  vx,vxc)
c***********************************************************************
c
      USE m_vxcl91
      USE m_vxcwb91
      USE m_vxcpw91
      USE m_vxcepbe
      USE m_types
      IMPLICIT NONE
c
c---> running mode parameters
c
      TYPE(t_xcpot),INTENT(IN) :: xcpot
      INTEGER, INTENT (IN) :: mfftwk,nfftwk ! radial mesh,number of mesh points
      INTEGER, INTENT (IN) :: jspins
      LOGICAL, INTENT (IN) :: lwbc          ! l-white-bird-current (ta)
c
c---> charge density
c
      REAL,INTENT (IN) :: agr(mfftwk)
      REAL,INTENT (IN) :: rh(mfftwk,jspins)
      REAL,INTENT (IN) :: agru(mfftwk),agrd(mfftwk)
      REAL,INTENT (IN) :: g2r(mfftwk),g2ru(mfftwk),g2rd(mfftwk)
      REAL,INTENT (IN) :: gggr(mfftwk),gggru(mfftwk)
      REAL, INTENT (IN) :: gggrd(mfftwk),gzgr(mfftwk)
c
c---> xc potential
c
      REAL, INTENT (OUT) :: vx (mfftwk,jspins)
      REAL, INTENT (OUT) :: vxc(mfftwk,jspins)
c
c ---> local scalars
      INTEGER ir,js
      REAL, PARAMETER :: hrtr_half = 0.5e0

      !used to be dummy arguments for testing
      INTEGER,PARAMETER   :: idsprs=0,isprsv=0
      REAL,PARAMETER      :: sprsv=0.0
c
c.....------------------------------------------------------------------
c
c-----> determine exchange correlation potential
c
      vx (:,:) = 0.0
      vxc(:,:) = 0.0

      IF (xcpot%is_name("l91")) THEN    ! local pw91

       CALL vxcl91(
     >             jspins,mfftwk,nfftwk,rh,agr,agru,agrd,
     >             g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <             vx,vxc,
     >             isprsv,sprsv)

      ELSEIF (xcpot%is_name("pw91")) THEN  ! pw91

       IF (lwbc) THEN
          CALL vxcwb91(
     >                 jspins,mfftwk,nfftwk,rh,agr,agru,agrd,
     >                 g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                 vx,vxc,
     >                 idsprs,isprsv,sprsv)
       ELSE

          CALL vxcpw91(
     >                 jspins,mfftwk,nfftwk,rh,agr,agru,agrd,
     >                 g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                 vx,vxc,
     >                 idsprs,isprsv,sprsv)

       ENDIF

      ELSE  ! pbe or similar

        CALL vxcepbe(
     >               xcpot,jspins,mfftwk,nfftwk,rh,
     >               agr,agru,agrd,g2ru,g2rd,gggr,gggru,gggrd,
     <               vx,vxc)


      ENDIF

c
c-----> hartree units
c
      IF (jspins.EQ.2) THEN
          DO ir = 1,nfftwk
              vx(ir,1)      = hrtr_half*vx(ir,1)
              vx(ir,jspins) = hrtr_half*vx(ir,jspins)        
          
              vxc(ir,1)      = hrtr_half*vxc(ir,1)
              vxc(ir,jspins) = hrtr_half*vxc(ir,jspins)
          ENDDO
      ELSEIF (jspins.eq.1) THEN
          DO ir = 1,nfftwk
              vx(ir,1)  = hrtr_half*vx(ir,1)
              
              vxc(ir,1) = hrtr_half*vxc(ir,1)
          ENDDO
      ELSE
          WRITE (6,fmt='('' error in jspins, jspins ='',i2)')
     +      jspins
          CALL juDFT_error("error in jspins",calledby ="xcallg")
      ENDIF

      END SUBROUTINE vxcallg
!*********************************************************************
      SUBROUTINE excallg(
     >                  xcpot,lwbc,jspins,nfftwk,
     >                  rh,agr,agru,agrd,
     >                  g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                  exc)
c***********************************************************************
c
      USE m_excl91
      USE m_excwb91
      USE m_excpw91
      USE m_excepbe
      USE m_types
   
      IMPLICIT NONE
c
c---> running mode parameters
c
      TYPE(t_xcpot),INTENT(IN) :: xcpot
      LOGICAL, INTENT (IN) :: lwbc          ! l-white-bird-current (ta)
      INTEGER, INTENT (IN) :: nfftwk ! radial mesh,number of mesh points
      INTEGER, INTENT (IN) :: jspins
c
c---> charge density  & gradients
c
      REAL,INTENT (IN) :: agr(:)
      REAL,INTENT (IN) :: rh(:,:)
      REAL,INTENT (IN) :: agru(size(agr)),agrd(size(agr))
      REAL,INTENT (IN) :: g2r(size(agr)),g2ru(size(agr)),g2rd(size(agr))
      REAL,INTENT (IN) :: gggr(size(agr)),gggru(size(agr))
      REAL,INTENT (IN) :: gggrd(size(agr)),gzgr(size(agr))
c
c---> xc energy density
c
      REAL, INTENT (OUT) :: exc(size(agr))
c
c ---> local scalars
      INTEGER ir
      REAL, PARAMETER :: hrtr_half = 0.5e0

      !used to be dummy arguments for testing
      INTEGER,PARAMETER   :: idsprs=0,isprsv=0
      REAL,PARAMETER      :: sprsv=0.0
c
c-----> determine exchange correlation energy density
c
c     write(6,'(/'' icorr,krla,igrd,jspins,lwbc,='',5i5,l2)')
c    &  icorr,krla,igrd,jspins,lwbc

      IF (xcpot%is_name("l91")) THEN  ! local pw91

        CALL excl91(
     >              jspins,size(rh,1),nfftwk,rh,agr,agru,agrd,
     >              g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <              exc,
     >              isprsv,sprsv)

      ELSEIF (xcpot%is_name("pw91")) THEN     ! pw91

        IF (lwbc) THEN
          CALL excwb91(
     >                 size(agr),nfftwk,
     >                 rh(:,1),rh(:,2),agr,agru,agrd,
     >                 g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                 exc,
     >                 idsprs,isprsv,sprsv)
        ELSE

          CALL excpw91(
     >                 jspins,size(agr),nfftwk,rh,agr,agru,agrd,
     >                 g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                 exc,
     >                 idsprs,isprsv,sprsv)

         ENDIF

      ELSE
        CALL excepbe(
     >               xcpot,jspins,size(rh,1),nfftwk,
     >               rh,agr,agru,agrd,g2ru,g2rd,gggr,gggru,gggrd,
     <               exc)

      ENDIF
c
c-----> hartree units
c
      DO ir = 1,nfftwk
        exc(ir) = hrtr_half*exc(ir)
      ENDDO


      END SUBROUTINE excallg
      END MODULE m_xcallg
