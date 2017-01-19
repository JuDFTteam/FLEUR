!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vgen
  USE m_juDFT
CONTAINS
  SUBROUTINE vgen(reap,input,xcpot,dimension, atoms,sphhar,stars,&
       vacuum,sym, obsolete,cell,oneD,sliceplot,mpi, results,noco)
    !     ***********************************************************
    !     FLAPW potential generator                           *
    !     ***********************************************************
    !     calculates the density-potential integrals needed for the
    !     total energy
    !     TE_VCOUL  :   charge density-coulomb potential integral
    !     TE_VEFF:   charge density-effective potential integral
    !     TE_EXC :   charge density-ex-corr.energy density integral
    !     ***********************************************************
#include"cpp_double.h"
    USE m_constants
    USE m_vmts
    USE m_intnv
    USE m_vmtxcg
    USE m_vmtxc
    USE m_vvacxc
    USE m_vvacxcg
    USE m_visxc
    USE m_visxcg
    USE m_vvac
    USE m_vvacis
    USE m_vvacxy
    USE m_vintcz
    USE m_checkdop
    USE m_wrtdop
    USE m_cdn_io
    USE m_qfix
    USE m_types
    USE m_od_vvac
    USE m_od_vvacis
    USE m_cylpts
    USE m_convol
    USE m_xyavden
    USE m_psqpw
    USE m_potmod
    USE m_intgr,         ONLY : intgr3
    USE m_hybridmix
    USE m_icorrkeys
    USE m_cfft
    USE m_sphpts
    USE m_points
    USE m_fleur_vdw
    IMPLICIT NONE
    TYPE(t_results),INTENT(INOUT)   :: results
    TYPE(t_xcpot),INTENT(IN)        :: xcpot
    TYPE(t_mpi),INTENT(IN)          :: mpi
    TYPE(t_dimension),INTENT(IN)    :: dimension
    TYPE(t_oneD),INTENT(IN)         :: oneD
    TYPE(t_obsolete),INTENT(IN)     :: obsolete
    TYPE(t_sliceplot),INTENT(IN)    :: sliceplot
    TYPE(t_input),INTENT(INOUT)     :: input  !efield can be modified
    TYPE(t_vacuum),INTENT(IN)       :: vacuum
    TYPE(t_noco),INTENT(IN)         :: noco
    TYPE(t_sym),INTENT(IN)          :: sym
    TYPE(t_stars),INTENT(IN)        :: stars
    TYPE(t_cell),INTENT(IN)         :: cell
    TYPE(t_sphhar),INTENT(IN)       :: sphhar
    TYPE(t_atoms),INTENT(INOUT)     :: atoms !vr0 is updated
    !     ..
    !     .. Scalar Arguments ..
    LOGICAL, INTENT (IN) :: reap       
    !     ..
    !     .. Local Scalars ..
    COMPLEX vintcza,xint,rhobar
    INTEGER i,i3,irec2,irec3,ivac,j,js,k,k3,lh,n,nzst1
    INTEGER imz,imzxy,ichsmrg,ivfft,npd 
    INTEGER ifftd,ifftd2, ifftxc3d,iter,datend
    INTEGER itypsym,itype,jsp,l,nat,archiveType
    !      INTEGER i_sm,n_sm,i_sta,i_end
    REAL ani,g3,signum,z,rhmn,fix,mfie
    REAL sig1dh,vz1dh,zat_l(atoms%ntype),rdum,dpdot ! ,delta,deltb,corr
    LOGICAL l_pottot,l_vdw
    LOGICAL exi
    LOGICAL, PARAMETER :: l_xyav=.FALSE.
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: alphm(:,:)
    COMPLEX, ALLOCATABLE :: excpw(:),excxy(:,:,:),vpw_w(:,:),psq(:)
    REAL,    ALLOCATABLE :: vbar(:),af1(:),bf1(:),xp(:,:)
    REAL,    ALLOCATABLE :: rhoc(:),rhoc_vx(:)
    !.....density
    REAL,    ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:)
    COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
    ! ff
    COMPLEX, ALLOCATABLE :: cdom(:)
    COMPLEX, ALLOCATABLE :: cdomvz(:,:),cdomvxy(:,:,:)
    !.....potential
    REAL,    ALLOCATABLE :: vr(:,:,:,:),vz(:,:,:),vr_exx(:,:,:,:)
    REAL,    ALLOCATABLE :: vxr(:,:,:,:),veffr(:,:,:,:)
    COMPLEX, ALLOCATABLE :: vpw(:,:),vxy(:,:,:,:)
    COMPLEX, ALLOCATABLE :: vpw_exx(:,:),vpw_wexx(:,:)
    COMPLEX, ALLOCATABLE :: vxpw(:,:),vxpw_w(:,:),veffpw_w(:,:)

    !.....energy density
    REAL,    ALLOCATABLE :: excz(:,:),excr(:,:,:)
    !
    ! if you want to calculate potential gradients
    !
    !pg      COMPLEX vlm(-lmaxd:lmaxd,0:lmaxd,ntypd)
    !pg      REAL fint,f(jmtd)
    !pg      INTEGER ns, nl, l, jm, m, mb
    !     ..
    !     ..
    !     ..
    !     ----> note the following conventions:
    !     ivac=2: lower (negative z) vacuum
    !     ivac=1: upper (positive z) vacuum
    !     units: hartrees
    !
    ALLOCATE ( alphm(stars%n2d,2),excpw(stars%n3d),excxy(vacuum%nmzxyd,oneD%odi%n2d-1,2),&
         vbar(dimension%jspd),af1(3*stars%k3d),bf1(3*stars%k3d),xp(3,dimension%nspd),&
         rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd),rht(vacuum%nmzd,2,dimension%jspd),&
         qpw(stars%n3d,dimension%jspd),rhtxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,dimension%jspd),&
         vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd),vz(vacuum%nmzd,2,dimension%jspd),&
         vpw(stars%n3d,dimension%jspd),vpw_exx(stars%n3d,dimension%jspd),vpw_wexx(stars%n3d,dimension%jspd),&
         vxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,dimension%jspd),&
         excz(vacuum%nmzd,2),excr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype),&
         vpw_w(stars%n3d,dimension%jspd),psq(stars%n3d) )

    ALLOCATE( vxpw(stars%n3d,dimension%jspd),vxpw_w(stars%n3d,dimension%jspd) )
    ALLOCATE( vxr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd) )
    vxr  = 0

    IF (noco%l_noco) THEN
       ALLOCATE ( cdom(stars%n3d), cdomvz(vacuum%nmzd,2),cdomvxy(vacuum%nmzxyd,oneD%odi%n2d-1,2) )
       archiveType = CDN_ARCHIVE_TYPE_NOCO_const
    ELSE
       ALLOCATE ( cdom(1),cdomvz(1,1),cdomvxy(1,1,1) )
       archiveType = CDN_ARCHIVE_TYPE_CDN1_const
    END IF
    !

    IF (mpi%irank == 0) THEN
       !
       ! --  total = .false. and reap = .false. means, that we now calculate
       !     the potential from the output density
       !
       IF ((.NOT.input%total).AND.(.NOT.reap)) THEN
          IF (noco%l_noco) THEN
             CALL juDFT_error("vgen:1",calledby ="vgen")
          ENDIF
          CALL readDensity(stars,vacuum,atoms,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN1_const,CDN_OUTPUT_DEN_const,&
                           0,iter,rho,qpw,rht,rhtxy,cdom,cdomvz,cdomvxy)
       ELSE
          CALL readDensity(stars,vacuum,atoms,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                           0,iter,rho,qpw,rht,rhtxy,cdom,cdomvz,cdomvxy)
       END IF

       IF (.NOT.l_xyav) THEN
          CALL timestart("Qfix")
          CALL qfix(stars,atoms,sym,vacuum, sphhar,input,cell,oneD, qpw,rhtxy,rho,rht,.FALSE., fix)
          CALL timestop("Qfix")
       ENDIF

       IF (input%total.OR.reap) THEN
          CALL writeDensity(stars,vacuum,atoms,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                            iter,rho,qpw,rht,rhtxy,cdom,cdomvz,cdomvxy)
       END IF

       WRITE (6,FMT=8000)
8000   FORMAT (/,/,t10,' p o t e n t i a l   g e n e r a t o r',/)
       vpw  = CMPLX(0.,0.)
       vxpw = CMPLX(0.,0.)
       vxpw_w = CMPLX(0.,0.)
       !     ---> perform spin summation of charge densities
       !     ---> for the calculation of the coulomb potentials
       IF ( (input%jspins.EQ.2).AND.(.NOT.l_xyav) ) THEN
          rho(:,0:,:,1)=rho(:,0:,:,1)+rho(:,0:,:,2)
          qpw(:,1)=qpw(:,1)+qpw(:,2)
          IF (input%film) THEN
             rhtxy(:,:,:,1)=rhtxy(:,:,:,1)+rhtxy(:,:,:,2)
             rht(:,:,1) = rht(:,:,1) + rht(:,:,2)
          END IF
       END IF
       !
       !     ************** coulomb potential ***********************
       IF (l_xyav) THEN
          zat_l(:) = 0.          ! for xy-averaged densities do not
       ELSE                     ! include nuclear charge
          zat_l(:) = atoms%zatom(:)
       ENDIF

       !     ----> create pesudo-charge density coefficients
    ENDIF ! (mpi%irank == 0)
    CALL timestart("psqpw")      
    CALL psqpw(mpi, atoms,sphhar,stars,vacuum, dimension,cell,input,sym,oneD,&
         qpw,rho,rht,l_xyav, psq)
    CALL timestop("psqpw")
    IF (mpi%irank == 0) THEN

       IF (l_xyav) THEN        ! write out xy-averaged density & stop
          CALL xy_av_den(stars,vacuum, cell,psq,rht)
          CALL juDFT_error("xy-averaged density calculated", calledby ="vgen")
       ENDIF

       !     ------------------------------------------
       IF (oneD%odi%d1) THEN
          !-odim
          CALL timestart("coulomb potential")

          !---> generates the m=0,gz=0 component of the vacuum potential
          CALL od_vvac(stars,vacuum,cell, psq,rht, vz)

          !---> generation of the vacuum warped potential components and
          !---> interstitial pw potential
          !---> vvacxy_5.F is a symmetrized potential generator

          CALL od_vvacis(oneD%odi%n2d,dimension,vacuum,oneD%odi%nq2,&
               oneD%odi%kv,cell,oneD%odi%M,stars,oneD%odi%nst2,&
               oneD, rht,rhtxy,psq,vz,sym, vxy,vpw)
          CALL timestop("coulomb potential")

          !+odim
       ELSEIF (input%film .AND. .NOT.oneD%odi%d1) THEN
          !     ----> potential in the  vacuum  region
          !       
          CALL timestart("p vac") 
          CALL vvac(vacuum,stars, cell,sym,input, psq,rht, vz,rhobar,sig1dh,vz1dh)
          CALL vvacis(stars,vacuum, sym,cell, psq, input, vxy)

          CALL vvacxy(stars,vacuum,cell,sym,input, rhtxy, vxy, alphm)
          CALL timestop("p vac")
       END IF
       !     ------------------------------------------
       !     ----> potential in the  interstitial  region
       CALL timestart("p int")
       WRITE (6,FMT=8010)
8010   FORMAT (/,5x,'coulomb potential in the interstitial region:')
       IF (input%film .AND. .NOT.oneD%odi%d1) THEN
          !           -----> create v(z) for each 2-d reciprocal vector
          ivfft =  3*stars%k3d 
          !         ivfft = 2*mx3 - 1
          ani = 1.0/REAL(ivfft)
          DO  irec2 = 1,stars%ng2
             i = 0
             DO i3 = 0,ivfft - 1
                i = i + 1
                z = cell%amat(3,3)*i3*ani
                IF (z.GT.cell%amat(3,3)/2.) z = z - cell%amat(3,3)
                vintcza = vintcz(stars,vacuum,cell,sym,input,&
                     z,irec2, psq,vxy,vz,rhobar,sig1dh,vz1dh,alphm)
                af1(i) = REAL(vintcza)
                bf1(i) = AIMAG(vintcza)
             ENDDO
             !                z = (i_sm-1)*ani
             !                IF (z > 0.5) z = z - 1.0
             !                af1(i_sm) = af1(i_sm) + z * delta
             !                bf1(i_sm) = bf1(i_sm) + z * deltb
             !              ENDDO
             !            ENDIF
             !        --> 1-d fourier transform and store the coefficients in vpw( ,1)
             CALL cfft(af1,bf1,ivfft,ivfft,ivfft,-1)
             !            delta = ivfft * delta * 2 / fpi ! * amat(3,3)**2 * ani
             i = 0
             DO  i3 = 0,ivfft - 1
                k3 = i3
                IF (k3 > FLOOR(ivfft/2.) ) k3 = k3 - ivfft
                i = i + 1
                IF ((k3.GE.-stars%mx3).AND.(k3.LE.stars%mx3)) THEN
                   irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),k3)

                   !                 IF ( (irec2 == 1).AND.(i3 > 0) ) THEN                 ! smooth potential
                   !                   corr = 2.0*mod(abs(k3),2) - 1.0
                   !                   bf1(i) = bf1(i) + delta * corr / k3
                   !                 ENDIF

                   !       ----> only stars within g_max sphere (shz oct.97)
                   IF (irec3.NE.0) THEN
                      !
                      xint = CMPLX(af1(i),bf1(i))*ani
                      nzst1 = stars%nstr(irec3)/stars%nstr2(irec2)
                      vpw(irec3,1) = vpw(irec3,1) + xint/nzst1
                   END IF
                ENDIF
             ENDDO
          ENDDO
       ELSEIF (.NOT.input%film) THEN
          vpw(1,1) = CMPLX(0.0,0.0)
          vpw(2:stars%ng3,1)=fpi_const*psq(2:stars%ng3)/(stars%sk3(2:stars%ng3)*stars%sk3(2:stars%ng3))       
       END IF

       CALL timestop("p int")

    ENDIF ! mpi%irank == 0
    !     --------------------------------------------
    !     ---> potential in the muffin-tin spheres

    CALL timestart("p vmts")
    CALL vmts(mpi, stars,sphhar,atoms, sym,cell,oneD, vpw,rho, vr)
    !     --------------------------------------------
    CALL timestop("p vmts")
    IF (mpi%irank == 0) THEN
       !     ---> check continuity of coulomb potential

       IF (input%vchk) THEN
          CALL timestart("checking")
          !           ----> vacuum boundaries
          IF (input%film .AND. .NOT.oneD%odi%d1) THEN
             npd = MIN(dimension%nspd,25)
             CALL points(xp,npd)
             DO ivac = 1,vacuum%nvac
                signum = 3. - 2.*ivac
                xp(3,:npd) = signum*cell%z1/cell%amat(3,3)
                CALL checkdop(xp,npd,0,0,ivac,1,1,.FALSE.,dimension,atoms, sphhar,stars,sym,&
                     vacuum,cell,oneD, vpw,vr,vxy,vz)
             ENDDO
          ELSEIF (oneD%odi%d1) THEN
             !-odim
             npd = MIN(dimension%nspd,25)
             CALL cylpts(xp,npd,cell%z1)
             !           DO j = 1,npd
             !              xp(1,j) = xp(1,j)/amat(1,1)
             !              xp(2,j) = xp(2,j)/amat(2,2)
             !           ENDDO
             CALL checkdop(xp,npd,0,0,vacuum%nvac,1,1,.FALSE.,dimension,atoms,&
                  sphhar,stars,sym, vacuum,cell,oneD, vpw,vr,vxy,vz)
             !+odim
          END IF
          !           ----> m.t. boundaries
          nat = 1
          DO  n = 1,atoms%ntype
             CALL sphpts(xp,dimension%nspd,atoms%rmt(n),atoms%pos(1,nat))
             CALL checkdop(xp,dimension%nspd,n,nat,0,-1,1,.FALSE.,dimension,atoms,&
                  sphhar,stars,sym, vacuum,cell,oneD, vpw,vr,vxy,vz)
             nat = nat + atoms%neq(n)
          ENDDO
          CALL timestop("checking")
       END IF
       !
       !========TOTAL==============================================
       !
       !      IF (l_xyav) THEN        ! write out xy-averaged potential & stop
       !        CALL xy_av_den(
       !     >                 n3d,k3d,nq3,nmzd,nmz,dvac,delz,
       !     >                 area,ig2,kv3,amat,vpw,vz(1,1,1))
       !         CALL juDFT_error("xy-averaged potential calculated",calledby="vgen")
       !      ENDIF

       IF (input%total) THEN
          CALL timestart("int_nv")

          !
          !      ---> AVERAGE COULOMB POTENTIAL ON THE SPHERE 
          !          FOR CALCULATING THE MADELUNG TERM in totale.f
          !           r=Rmt
          DO n=1,atoms%ntype
             atoms%vr0(n)=vr(atoms%jri(n),0,n,1)
          ENDDO
          !
          !     CALCULATE THE INTEGRAL OF n*Vcoulomb
          !
          WRITE (6,FMT=8020)
          WRITE (16,FMT=8020)
8020      FORMAT (/,10x,'density-coulomb potential integrals',/)
          !
          !       interstitial first
          !
          !       convolute ufft and pot: F(G) = \sum_(G') U(G - G') V(G')
          !
          CALL convol(stars, vpw_w, vpw)
          !
          IF (input%jspins.EQ.2) CALL CPP_BLAS_ccopy(stars%ng3,vpw_w(1,1),1,vpw_w(1,input%jspins),1)
          !
          results%te_vcoul = 0.0
          CALL int_nv(stars,vacuum,atoms,sphhar, cell,sym,input,oneD,&
               qpw,vpw_w, rhtxy,vxy, rht,vz, rho,vr, results%te_vcoul)

          WRITE (6,FMT=8030) results%te_vcoul
          WRITE (16,FMT=8030) results%te_vcoul
8030      FORMAT (/,10x,'total density-coulomb potential integral :', t40,f20.10)

          CALL timestop("int_nv")

          INQUIRE(file='vdW_kernel_table',exist=l_vdw)
          IF (l_vdw) THEN

             CALL timestart("fleur_vdW")
             ! calculate vdW contribution to potential
             CALL fleur_vdW(mpi,atoms,sphhar,stars, input,dimension,&
                  cell,sym,oneD,vacuum, qpw(:,1),rho(:,:,:,1), vpw_w(:,1),vr(:,:,:,:))
             CALL timestop("fleur_vdW")

          ENDIF

       END IF
       !ENDIF !irank==0
       !
       !==========END TOTAL===================================================
       !
       !     ----> reload the density for calculating vxc (for spin-pol. case)
       !
       IF (input%jspins.EQ.2) THEN
          CALL readDensity(stars,vacuum,atoms,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                           0,iter,rho,qpw,rht,rhtxy,cdom,cdomvz,cdomvxy)
          vr(:,0:,:,2) = vr(:,0:,:,1)
          vpw(:,2) = vpw(:,1)
          IF (input%film) THEN
             vxy(:,:,:,2) = vxy(:,:,:,1)
             vz(:,:,2)=vz(:,:,1)
          END IF
       END IF
       IF (input%total) THEN
          OPEN (11,file='potcoul',form='unformatted',status='unknown')
          DO js = 1,input%jspins
             ! to enable a GW calculation,
             vpw_w(1:stars%ng3,js)=vpw_w(1:stars%ng3,js)/stars%nstr(1:stars%ng3)     ! the PW-coulomb part is not
             ! used otherwise anyway.
          ENDDO
          CALL wrtdop(stars,vacuum,atoms,sphhar, input,sym,&
               11, iter,vr,vpw_w,vz,vxy)
          DO js = 1,input%jspins
             DO i = 1,stars%ng3
                vpw_w(i,js)=vpw_w(i,js)*stars%nstr(i)
             ENDDO
          ENDDO
          CLOSE(11)
       END IF
       IF (sliceplot%plpot) THEN
          OPEN (11,file='potcoul_pl',form='unformatted',status='unknown')
          CALL wrtdop(stars,vacuum,atoms,sphhar, input,sym,&
               11, iter,vr,vpw,vz,vxy)
          CLOSE(11)
       END IF

       !     ******** exchange correlation potential******************
       !+ta
       !     rhmn: rho-min.
       !     ichsmrg: i-charge-small-region.
       !          0(not watched)
       !          1(in muffin-tin), 2(interstitial), 3(warped-vac),4(vacuum)

       ichsmrg=0
       rhmn=1.e+10
       !-ta
       excpw(:) = CMPLX(0.,0.)
       excz(:,:) = 0.0
       excxy(:,:,:) = CMPLX(0.,0.)
       excr(:,:,:) = 0.0

       !     ---> vacuum region
       IF (input%film) THEN

          CALL timestart("Vxc in vacuum")

          ifftd2 = 9*stars%k1d*stars%k2d
          IF (oneD%odi%d1) ifftd2 = 9*stars%k3d*oneD%odi%M

          IF ((xcpot%igrd == 0).AND.(xcpot%icorr /= -1)) THEN  ! LDA

             IF (.NOT.oneD%odi%d1) THEN

                CALL vvacxc(ifftd2,stars,vacuum,xcpot,input,noco,&
                     rhtxy,rht,cdomvxy,cdomvz, vxy,vz, excxy,excz)

             ELSE
                CALL judft_error("OneD broken")
                !           CALL vvacxc(&
                !     &                 stars,oneD%M,vacuum,odi%n2d,dimension,ifftd2,&
                !     &                 xcpot,input,odi%nq2,&
                !     &                 odi%nst2,rhtxy,rht,cdomvxy,cdomvz,noco,&
                !     &                 odi%kimax2%igf,odl%pgf,&
                !     &                 vxy,vz,&
                !     &                 excxy,excz)

             ENDIF

          ELSE      ! GGA

             IF (oneD%odi%d1) THEN
                CALL judft_error("OneD broken")

                CALL vvacxcg(ifftd2,stars,vacuum,noco,oneD,&
                     cell,xcpot,input,obsolete, ichsmrg,&
                     rhtxy,rht,cdomvxy,cdomvz, vxy,vz,rhmn, excxy,excz)

             ELSE
                CALL vvacxcg(ifftd2,stars,vacuum,noco,oneD,&
                     cell,xcpot,input,obsolete, ichsmrg,&
                     rhtxy,rht,cdomvxy,cdomvz, vxy,vz,rhmn, excxy,excz)

             END IF

          END IF
          CALL timestop("Vxc in vacuum")
          !+odim
       END IF
       !     ----------------------------------------
       !     ---> interstitial region
       CALL timestart("Vxc in interstitial")

       ifftd=27*stars%k1d*stars%k2d*stars%k3d

       IF ( (.NOT. obsolete%lwb) .OR. ( (xcpot%igrd == 0) .AND. (xcpot%icorr /= -1) ) ) THEN
          ! no White-Bird-trick

          ifftxc3d = stars%kxc1d*stars%kxc2d*stars%kxc3d

          IF ( (xcpot%igrd == 0) .AND. (xcpot%icorr /= -1) ) THEN
             ! LDA

             CALL visxc(ifftd,stars,noco,xcpot,input, qpw,cdom,&
                  vpw,vpw_w,vxpw,vxpw_w, excpw)

          ELSE ! GGA

             CALL visxcg(ifftd,stars,sym, ifftxc3d, cell, qpw,cdom, xcpot,input,&
                  obsolete,noco, rhmn,ichsmrg, vpw,vpw_w,vxpw,vxpw_w, excpw)

          END IF

       ELSE
          ! White-Bird-trick

          WRITE(6,'(a)') "W+B trick cancelled out. visxcwb uses at present common block cpgft3.",&
             "visxcwb needs to be reprogrammed according to visxcg.f"
          CALL juDFT_error("visxcwb",calledby ="vgen")
          !sb       CALL visxcwb(
          !sb  >                 qpw,kimax,igfft,pgfft,ufft,
          !sb  >                 icorr,total,krla,
          !sb  >                 igrd,ndvgrd,idsprs,isprsv,
          !sb  >                 idsprsi,chng,sprsv,lwb,rhmn,ichsmrg,
          !sb  =                 vpw,vpw_w,
          !sb  <                 excpw)

       END IF
       !
       ! --> on output vpw_w contains the warped effective potential and
       !               excpw the warped XC-energy density
       !

       !
       !     add IR EXX potential to vpw_w
       !
       INQUIRE(file='vpw_wexx',exist=exi)
       IF( exi ) THEN 
          WRITE(*,*) 'Read in vpw_wexx...'
          OPEN (351,file='vpw_wexx',form='formatted')
          DO js = 1,input%jspins
             DO i = 1,stars%ng3
                READ(351,'(2f30.15)') vpw_wexx(i,js) 
                vpw_w(i,js) = vpw_w(i,js) + vpw_wexx(i,js)*stars%nstr(i)
             END DO
          END DO
          CLOSE(351)
       END IF

       INQUIRE(file='vpw_exx',exist=exi)
       IF( exi ) THEN 
          WRITE(*,*) 'Read in vpw_exx...'
          OPEN (351,file='vpw_exx',form='formatted')
          DO js = 1,input%jspins
             DO i = 1,stars%ng3
                READ(351,'(2f30.15)') vpw_exx(i,js) 
                vpw(i,js) = vpw(i,js) + vpw_exx(i,js)
             END DO
          END DO
          CLOSE(351)
       END IF
       CALL timestop("Vxc in interstitial")

       !      ---> evaluate the interstitial average potential
       vbar(:) = vpw_w(1,:)*cell%omtil/cell%volint
       !-gu
       WRITE (6,FMT=8040) (vbar(js),js=1,input%jspins)
       WRITE (16,FMT=8040) (vbar(js),js=1,input%jspins)
8040   FORMAT (/,5x,'interstitial potential average (vbar) =',2f10.6)
       !
       !     ------------------------------------------
       !     ----> muffin tin spheres region

       CALL timestart ("Vxc in MT")
       IF ((xcpot%igrd.EQ.0).AND.(xcpot%icorr.NE.-1)) THEN
          CALL vmtxc(dimension,sphhar,atoms, rho,xcpot,input,sym, vr, excr,vxr)
       ELSEIF ((xcpot%igrd.GT.0).OR.(xcpot%icorr.EQ.-1)) THEN
          CALL vmtxcg(dimension,sphhar,atoms, rho,xcpot,input,sym,&
               obsolete, vxr,vr,rhmn,ichsmrg, excr)
       ELSE
          CALL juDFT_error("something wrong with xcpot before vmtxc" ,calledby ="vgen")
       ENDIF


       !
       ! add MT EXX potential to vr
       !

       INQUIRE(file='vr_exx',exist=exi)
       IF( exi ) THEN
          ALLOCATE( vr_exx(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd) )
          OPEN (350,file='vr_exx',form='formatted')
          DO jsp = 1,dimension%jspd
             DO  itype = 1,atoms%ntype
                itypsym = atoms%ntypsy( SUM(atoms%neq(:itype-1)) + 1 )
                DO lh = 0,sphhar%nlh(itypsym)
                   l = sphhar%llh(lh,itypsym)
                   DO i = 1,atoms%jmtd
                      READ(350,'(f30.15)') vr_exx(i,lh,itype,jsp)
                   END DO

                   IF( l .EQ. 0 ) THEN
                      vr_exx(:,lh,itype,jsp) = vr_exx(:,lh,itype,jsp)*sfp_const/atoms%rmsh(:,itype)
                   END IF

                END DO
             END DO
          END DO
          CLOSE(350)

          vr = vr + vr_exx

       END IF
       CALL timestop ("Vxc in MT")
       !     ------------------------------------------
       !     ---> check continuity of total potential

       IF (input%vchk) THEN
          !           ----> vacuum boundaries
          IF (input%film .AND. .NOT.oneD%odi%d1) THEN
             npd = MIN(dimension%nspd,25)
             CALL points(xp,npd)
             DO ivac = 1,vacuum%nvac
                signum = 3. - 2.*ivac
                xp(3,:npd) = signum*cell%z1/cell%amat(3,3)
                CALL checkdop(xp,npd,0,0,ivac,1,1,.FALSE.,dimension,atoms, sphhar,stars,sym,&
                     vacuum,cell,oneD, vpw,vr,vxy,vz)
             ENDDO ! ivac = 1,vacuum%nvac
          ELSEIF (oneD%odi%d1) THEN
             !-odim
             npd = MIN(dimension%nspd,25)
             CALL cylpts(xp,npd,cell%z1)
             !           DO j = 1,npd
             !              xp(1,j) = xp(1,j)/amat(1,1)
             !              xp(2,j) = xp(2,j)/amat(2,2)
             !           ENDDO
             CALL checkdop(xp,npd,0,0,vacuum%nvac,1,1,.FALSE.,dimension,atoms,&
                  sphhar,stars,sym, vacuum,cell,oneD, vpw,vr,vxy,vz)
             !+odim
          END IF
          !           ----> m.t. boundaries
          nat = 1
          DO n = 1, atoms%ntype
             CALL sphpts(xp,dimension%nspd,atoms%rmt(n),atoms%pos(1,nat))
             CALL checkdop(xp,dimension%nspd,n,nat,0,-1,1,.FALSE.,dimension,&
                  atoms,sphhar,stars,sym, vacuum,cell,oneD, vpw,vr,vxy,vz)
             nat = nat + atoms%neq(n)
          ENDDO ! n = 1, atoms%ntype
       END IF


       CALL pot_mod(atoms,sphhar,vacuum,stars, input, vr,vxy,vz,vpw,vpw_w)
       !
       !============TOTAL======================================
       !
       IF (input%total) THEN

          IF (noco%l_noco) THEN ! load qpw,rht,rhtxy from 'cdn'-file
             CALL readDensity(stars,vacuum,atoms,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                              0,iter,rho,qpw,rht,rhtxy,cdom,cdomvz,cdomvxy)
          ENDIF
          !
          !     CALCULATE THE INTEGRAL OF n1*Veff1 + n2*Veff2
          !     Veff = Vcoulomb + Vxc
          !
          ALLOCATE( veffpw_w(stars%n3d,dimension%jspd) )
          ALLOCATE( veffr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd) )
          IF( xcpot%icorr .EQ. icorr_pbe0 ) THEN
             veffpw_w = vpw_w - amix_pbe0 * vxpw_w
             veffr    = vr    - amix_pbe0 * vxr
          ELSE IF ( xcpot%icorr==icorr_hse .OR. xcpot%icorr==icorr_hseloc ) THEN
             veffpw_w = vpw_w - aMix_HSE * vxpw_w
             veffr    = vr    - aMix_HSE * vxr
          ELSE IF ( xcpot%icorr .EQ. icorr_vhse ) THEN
             veffpw_w = vpw_w - aMix_VHSE() * vxpw_w
             veffr    = vr    - aMix_VHSE() * vxr
          ELSE IF ( xcpot%icorr .EQ. icorr_hf ) THEN
             veffpw_w = vpw_w - amix_hf  * vxpw_w
             veffr    = vr    - amix_hf  * vxr
          ELSE
             veffpw_w = vpw_w
             veffr    = vr
          END IF

          !HF     kinetic energy correction for core states
          IF ( xcpot%icorr == icorr_pbe0 .OR. xcpot%icorr == icorr_hse .OR.&
               xcpot%icorr == icorr_hf   .OR. xcpot%icorr == icorr_vhse    ) THEN
             OPEN (17,file='cdnc',form='unformatted',status='unknown')
             REWIND 17
             ALLOCATE( rhoc(atoms%jmtd), rhoc_vx(atoms%jmtd) )
          END IF
          !HF

          results%te_veff = 0.0
          DO 370 js = 1,input%jspins
             WRITE (6,FMT=8050) js
             WRITE (16,FMT=8050) js
8050         FORMAT (/,10x,'density-effective potential integrals for spin ',i2,/)

             CALL int_nv(stars,vacuum,atoms,sphhar, cell,sym,input,oneD,&
                  qpw(:,js),veffpw_w(:,js), rhtxy(:,:,:,js),vxy(:,:,:,js),&
                  rht(:,:,js),vz(:,:,js), rho(1,0,1,js),veffr(1,0,1,js), results%te_veff)

             !HF
             IF ( xcpot%icorr == icorr_pbe0 .OR. xcpot%icorr == icorr_hse .OR.&
                  xcpot%icorr == icorr_hf   .OR. xcpot%icorr == icorr_vhse    ) THEN
                nat = 1
                DO n = 1,atoms%ntype
                   READ (17) ( rhoc(j), j = 1,atoms%jri(n) )
                   !             Skip over parts in cdnc not used
                   READ (17) rdum
                   !             calculate exchange correction to kinetic energy by core states
                   DO j = 1, atoms%jri(n)
                      rhoc_vx(j) = rhoc(j) * vxr(j,0,n,js) / sfp_const
                   END DO
                   CALL intgr3(rhoc_vx,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),dpdot)
                   IF( xcpot%icorr .EQ. icorr_pbe0 ) THEN
                      results%te_veff = results%te_veff + amix_pbe0 * dpdot*atoms%neq(n)
                   ELSE IF( xcpot%icorr .EQ. icorr_hse ) THEN
                      results%te_veff = results%te_veff + aMix_HSE * dpdot*atoms%neq(n)
                   ELSE IF( xcpot%icorr .EQ. icorr_vhse ) THEN
                      results%te_veff = results%te_veff + aMix_VHSE() * dpdot*atoms%neq(n)
                   ELSE IF( xcpot%icorr .EQ. icorr_hf  ) THEN
                      results%te_veff = results%te_veff + amix_hf * dpdot*atoms%neq(n)
                   END IF
                   nat     = nat + atoms%neq(n)
                END DO
                !           Skip over parts in cdnc not used
                READ (17) ( rdum, n = 1,atoms%ntype )
             END IF
             !HF

370          CONTINUE

             !HF
             IF ( xcpot%icorr == icorr_pbe0 .OR. xcpot%icorr == icorr_hse .OR.&
                  xcpot%icorr == icorr_hf   .OR. xcpot%icorr == icorr_vhse ) THEN
                CLOSE ( 17 )
                DEALLOCATE( rhoc, rhoc_vx )
             END IF
             !HF     end kinetic energy correction

             DEALLOCATE( veffpw_w,veffr )
             WRITE (6,FMT=8060) results%te_veff
             WRITE (16,FMT=8060) results%te_veff
8060         FORMAT (/,10x,'total density-effective potential integral :', t40,f20.10)
             !
             !     CALCULATE THE INTEGRAL OF n*exc
             !
             !     ---> perform spin summation of charge densities
             !     ---> for the calculation of Exc
             IF (input%jspins.EQ.2) THEN
                nat = 1
                DO n = 1,atoms%ntype
                   rho(:atoms%jri(n),0:sphhar%nlh(atoms%ntypsy(nat)),n,1) = rho(:atoms%jri(n),0:sphhar%nlh(atoms%ntypsy(nat)),n,1) + rho(:atoms%jri(n),0:sphhar%nlh(atoms%ntypsy(nat)),n,input%jspins)

                   nat = nat + atoms%neq(n)
                ENDDO
                qpw(:stars%ng3,1) = qpw(:stars%ng3,1) + qpw(:stars%ng3,input%jspins)
                IF (input%film) THEN
                   rhtxy(:vacuum%nmzxy,:oneD%odi%nq2 - 1,:vacuum%nvac,1) = &
                        rhtxy(:vacuum%nmzxy,:oneD%odi%nq2 - 1,:vacuum%nvac,1) + &
                        rhtxy(:vacuum%nmzxy,:oneD%odi%nq2 - 1,:vacuum%nvac,input%jspins)
                   rht(:vacuum%nmz,:vacuum%nvac,1) = rht(:vacuum%nmz,:vacuum%nvac,1) +&
                        rht(:vacuum%nmz,:vacuum%nvac,input%jspins)
                END IF
             END IF
             WRITE (6,FMT=8070)
             WRITE (16,FMT=8070)
8070         FORMAT (/,10x,'charge density-energy density integrals',/)

             results%te_exc = 0.0
             CALL int_nv(stars,vacuum,atoms,sphhar, cell,sym,input,oneD,&
                  qpw(:,1),excpw(1), rhtxy,excxy, rht,excz, rho,excr, results%te_exc)
             WRITE (6,FMT=8080) results%te_exc
             WRITE (16,FMT=8080) results%te_exc

8080         FORMAT (/,10x,'total charge density-energy density integral :', t40,f20.10)

          END IF
          !
          !==========END TOTAL============================================
          !
          !           ---> store v(l=0) component as r*v(l=0)/sqrt(4pi)

          DO js = 1,input%jspins

             l_pottot = .FALSE.                  ! adds a B-field, if file
             IF (dimension%jspd.EQ.2) THEN                       !   mfee is present
                INQUIRE (file='mfee',exist=l_pottot)
                IF (l_pottot) THEN
                   OPEN (88,file='mfee',form='formatted',status='unknown')
                   REWIND 88
                   DO n=1,atoms%ntype
                      READ (88,*) i,mfie
                      WRITE (*,*) 'type,field:',i,mfie
                      IF (i/=n)  CALL juDFT_error("wrong types in mfee", calledby="vgen")
                      IF (js.EQ.1) THEN
                         vr(:atoms%jri(n),0,n,js) = vr(:atoms%jri(n),0,n,js) - mfie/2.
                      ELSE
                         vr(:atoms%jri(n),0,n,js) = vr(:atoms%jri(n),0,n,js) + mfie/2.
                      ENDIF
                   ENDDO
                   CLOSE (88)
                ENDIF
             ENDIF



             DO n = 1,atoms%ntype
                vr(:atoms%jri(n),0,n,js)  = atoms%rmsh(:atoms%jri(n),n)*vr(:atoms%jri(n),0,n,js)/sfp_const
                vxr(:atoms%jri(n),0,n,js) = atoms%rmsh(:atoms%jri(n),n)*vxr(:atoms%jri(n),0,n,js)/sfp_const
             ENDDO

          ENDDO     ! js =1,input%jspins


          IF ((.NOT.reap).OR.(noco%l_noco)) THEN
             IF (input%total) THEN
                OPEN (9,file='nrp',form='unformatted',status='unknown')
             ELSE
                OPEN (9,file='nrp',form='unformatted',position='append')
             ENDIF
             CALL wrtdop(stars,vacuum,atoms,sphhar, input,sym,&
                  9, iter,vr,vpw,vz,vxy)
             CLOSE(9)
          ENDIF

          !     **************** reanalyze vpw *************************
          !                        call cpu_time(cp0)
          IF (input%total) THEN
             !     ----->write potential to file 8
             ! -> the following procedure is required in order to run
             !    correctly on the helga parallel cluster in Hamburg end 2005
             !                                                Paolo & YM
             l_pottot = .FALSE.
             INQUIRE (file='pottot',exist=l_pottot)
             IF (l_pottot) THEN
                OPEN (8,file='pottot',form='unformatted',status='unknown')
                CLOSE (8,status='delete')
                WRITE(6,*) 'vgen: pottot deleted'
             ENDIF
             OPEN (8,file='pottot',form='unformatted',status='unknown')
             REWIND 8
             DO js=1,input%jspins
                DO i=1,stars%ng3
                   vpw_w(i,js)=vpw_w(i,js)/stars%nstr(i)
                ENDDO
             ENDDO
             CALL wrtdop(stars,vacuum,atoms,sphhar, input,sym,&
                  8, iter,vr,vpw_w,vz,vxy) ! vpw_w
             CLOSE(8)

             OPEN (8,file='potx',form='unformatted',status='unknown')
             REWIND 8

             DO js=1,input%jspins
                DO i=1,stars%ng3
                   vxpw_w(i,js)=vxpw_w(i,js)/stars%nstr(i)
                ENDDO
             ENDDO

             CALL wrtdop(stars,vacuum,atoms,sphhar, input,sym,&
                  8, iter,vxr,vxpw_w,vz,vxy)
             CLOSE(8)


          END IF

       ENDIF ! mpi%irank == 0

       DEALLOCATE ( cdom,cdomvz,cdomvxy )
       DEALLOCATE (alphm,excpw,excxy,vbar,af1,bf1,xp,rho,rht,qpw,rhtxy,vr,vz,&
            vpw,vxy,excz,excr,vpw_w,psq)
       DEALLOCATE (vxpw,vxpw_w,vxr)

     END SUBROUTINE vgen
   END MODULE m_vgen
