MODULE m_vmtxcg
  !.....------------------------------------------------------------------
  !     fit spherical-harmonics expansion of exchange-correlation
  !     potential inside muffint-tin spheres and add it to coulomb
  !     potential
  !                                     c.l.fu and r.podloucky           *
  !     for the gradient correction. t.a. 1996.
  !.....------------------------------------------------------------------
  !     instead of vmtxcor.f: the different exchange-correlation
  !     potentials defined through the key icorr are called through
  !     the driver subroutine vxcallg.f, subroutines vectorized
  !     ** r.pentcheva 22.01.96
  !     *********************************************************
  !     angular mesh calculated on speacial gauss-legendre points
  !     in order to use orthogonality of lattice harmonics and
  !     avoid a least square fit
  !     ** r.pentcheva 04.03.96
  !     *********************************************************
  !     MPI and OpenMP parallelization
  !             U.Alekseeva, February 2017
  !     *********************************************************

CONTAINS
  SUBROUTINE vmtxcg(&
       &                dimension,mpi,sphhar,atoms,&
       &                rho,xcpot,input,sym,&
       &                obsolete,&
       &                vxr,vr,rhmn,ichsmrg,&
       &                excr)

#include"cpp_double.h"
    USE m_lhglptg
    USE m_grdchlh
    USE m_mkgylm
    USE m_gaussp
    USE m_xcallg, ONLY : vxcallg,excallg

    USE m_types
    IMPLICIT NONE

    TYPE(t_xcpot),INTENT(IN)       :: xcpot
    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_mpi),INTENT(IN)         :: mpi
    TYPE(t_obsolete),INTENT(IN)    :: obsolete
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_sphhar),INTENT(IN)      :: sphhar
    TYPE(t_atoms),INTENT(IN)       :: atoms
#ifdef CPP_MPI
    include "mpif.h"
#endif
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT(INOUT):: ichsmrg
    REAL,    INTENT(INOUT):: rhmn
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd)
    REAL,    INTENT (OUT):: vxr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd)
    REAL,    INTENT (INOUT):: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd)
    REAL,    INTENT (OUT)  :: excr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
    !     ..
    !     .. Local Scalars ..
    INTEGER jr,js,k,lh,n,nd,ist,nsp,ixpm ,i,nat
    REAL    rhmnm,d_15,elh,vlh
    LOGICAL lwbc              ! if true, white-bird trick
    !     ..
    !     .. Local Arrays ..
    REAL vx(dimension%nspd,dimension%jspd),vxc(dimension%nspd,dimension%jspd),exc(dimension%nspd),rx(3,dimension%nspd)
    REAL vxcl(DIMENSION%nspd,DIMENSION%jspd),excl(DIMENSION%nspd),divi
    TYPE(t_xcpot)::xcpot_tmp
    REAL wt(dimension%nspd),rr2(atoms%jmtd),thet(dimension%nspd)
    REAL agr(dimension%nspd),agru(dimension%nspd),agrd(dimension%nspd),g2r(dimension%nspd),g2ru(dimension%nspd)
    REAL g2rd(dimension%nspd),gggr(dimension%nspd),gggru(dimension%nspd),gggrd(dimension%nspd)
    REAL grgru(dimension%nspd),grgrd(dimension%nspd),gzgr(dimension%nspd)
    REAL, ALLOCATABLE :: ylh(:,:,:),ylht(:,:,:),ylhtt(:,:,:)
    REAL, ALLOCATABLE :: ylhf(:,:,:),ylhff(:,:,:),ylhtf(:,:,:)
    REAL, ALLOCATABLE :: chlh(:,:,:),chlhdr(:,:,:),chlhdrr(:,:,:)
    REAL, ALLOCATABLE :: ch(:,:),chdr(:,:),chdt(:,:),chdf(:,:)
    REAL, ALLOCATABLE :: chdrr(:,:),chdtt(:,:),chdff(:,:),chdtf(:,:)
    REAL, ALLOCATABLE :: chdrt(:,:),chdrf(:,:)

    !locals for mpi
    integer :: ierr
    integer:: n_start,n_stride
#ifdef CPP_MPI
    REAL :: vr_local(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd)
    REAL :: vxr_local(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd)  
    REAL :: excr_local(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
    INTEGER :: ichsmrg_local
    REAL :: rhmn_local, rhmn_reduced
#endif    
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC max,mod,min
    LOGICAL l_gga
    !ta+
    !.....------------------------------------------------------------------
    l_gga=xcpot%is_gga()

#ifdef CPP_MPI
    CALL MPI_BCAST(obsolete%ndvgrd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    vr_local = 0.d0
    vxr_local = 0.d0
    excr_local = 0.d0
    rhmn_local = rhmn
#endif
    vxr = 0.d0

    lwbc=.false.

    d_15 = 1.e-15
    ixpm = 1
    ist  = 1

    !     generates nspd points on a sherical shell with radius 1.0
    !     angular mesh equidistant in phi,
    !     theta are zeros of the legendre polynomials
    !
    CALL gaussp(&
         &            atoms%lmaxd,&
         &            rx,wt)

    nsp = dimension%nspd
    !
    !     generates the lattice harmonics on the angular mesh
    !
    ALLOCATE ( ylh(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd),ylht(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd),&
         &           ylhtt(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd),ylhf(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd),&
         &           ylhff(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd),ylhtf(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd) )

    CALL lhglptg(&
         &             sphhar,atoms,&
         &             rx,nsp,xcpot,sym,&
         &             ylh,thet,ylht,ylhtt,ylhf,ylhff,ylhtf)
    !
    !
    !     loop over topologically non-equivalent atoms
    !


#ifdef CPP_MPI
    n_start=mpi%irank+1
    n_stride=mpi%isize
#else
    n_start=1
    n_stride=1
#endif

    ALLOCATE ( chlh(atoms%jmtd,0:sphhar%nlhd,dimension%jspd),&
            &       chlhdr(atoms%jmtd,0:sphhar%nlhd,dimension%jspd),chlhdrr(atoms%jmtd,0:sphhar%nlhd,dimension%jspd))

    DO 200 n = n_start,atoms%ntype,n_stride
       IF (xcpot%lda_atom(n))THEN
          IF((.NOT.xcpot%is_name("pw91"))) CALL judft_error("Using locally LDA only possible with pw91 functional")
          CALL xcpot_tmp%init("l91",.FALSE.)
       ENDIF
       nat=sum(atoms%neq(:n-1))+1
       nd = atoms%ntypsy(nat)

       DO jr = 1,atoms%jri(n)
          rr2(jr) = 1.e0/ (atoms%rmsh(jr,n)*atoms%rmsh(jr,n))
       ENDDO

       DO lh = 0,sphhar%nlh(nd)

          !         calculates gradients of radial charge densities of l=> 0.
          !         rho*ylh/r**2 is charge density. chlh=rho/r**2.
          !         charge density=sum(chlh*ylh).
          !         chlhdr=d(chlh)/dr, chlhdrr=dd(chlh)/drr.

          DO js = 1,input%jspins

             DO jr = 1,atoms%jri(n)
                chlh(jr,lh,js) = rho(jr,lh,n,js)*rr2(jr)
             ENDDO

             IF (l_gga) THEN 
                CALL grdchlh(&
                     &                     ixpm,ist,atoms%jri(n),atoms%dx(n),atoms%rmsh(1,n),&
                     &                     chlh(1,lh,js),obsolete%ndvgrd,&
                     &                     chlhdr(1,lh,js),chlhdrr(1,lh,js))
             ENDIF

          ENDDO ! js
       ENDDO   ! lh
       
       !
       !-->    loop over radial mesh 
       !
       !$OMP PARALLEL DEFAULT(none) &
#ifdef CPP_MPI
       !$OMP& SHARED(vr_local,vxr_local,excr_local,ichsmrg_local,rhmn_local) &
#endif
       !$OMP& SHARED(vr,vxr,excr,rhmn,ichsmrg,l_gga) &
       !$OMP& SHARED(dimension,mpi,sphhar,atoms,rho,xcpot,input,sym,obsolete)&
       !$OMP& SHARED(n,nd,ist,ixpm,nsp,nat,d_15,lwbc) &
       !$OMP& SHARED(rx,wt,rr2,thet,xcpot_tmp) &
       !$OMP& SHARED(ylh,ylht,ylhtt,ylhf,ylhff,ylhtf) &
       !$OMP& SHARED(chlh,chlhdr,chlhdrr) &
       !$OMP& SHARED(ierr,n_start,n_stride) &
       !$OMP& PRIVATE(js,k,lh,i,rhmnm,elh,vlh) &
       !$OMP& PRIVATE(vx,vxc,exc,excl,vxcl,divi) &
       !$OMP& PRIVATE(agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,grgru,grgrd,gzgr) &
       !$OMP& PRIVATE(ch,chdr,chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf)
       ALLOCATE ( ch(dimension%nspd,dimension%jspd),chdr(dimension%nspd,dimension%jspd),chdt(dimension%nspd,dimension%jspd),&
            &       chdf(dimension%nspd,dimension%jspd),chdrr(dimension%nspd,dimension%jspd),chdtt(dimension%nspd,dimension%jspd),&
            &       chdff(dimension%nspd,dimension%jspd),chdtf(dimension%nspd,dimension%jspd),chdrt(dimension%nspd,dimension%jspd),&
            &       chdrf(dimension%nspd,dimension%jspd))
       !$OMP DO 
       DO 190 jr = 1,atoms%jri(n)
          !
          !         following are at points on jr-th sphere.
          !         initialize arrays
          !
          ch(:,:)    = 0.0     ! charge density
          chdr(:,:)  = 0.0     ! d(ch)/dr
          chdt(:,:)  = 0.0     ! d(ch)/dtheta
          chdf(:,:)  = 0.0     ! d(ch)/dfai
          chdrr(:,:) = 0.0     ! dd(ch)/drr
          chdtt(:,:) = 0.0     ! dd(ch)/dtt
          chdff(:,:) = 0.0     ! dd(ch)/dff
          chdtf(:,:) = 0.0     ! dd(ch)/dtf
          chdrt(:,:) = 0.0     ! d(d(ch)/dr)dt
          chdrf(:,:) = 0.0     ! d(d(ch)/dr)df

          DO js = 1,input%jspins
             DO lh = 0,sphhar%nlh(nd)

                !             generate the densities on an angular mesh

                DO k = 1,nsp
                   ch(k,js) = ch(k,js) + ylh(k,lh,nd)*chlh(jr,lh,js)
                ENDDO

                IF (l_gga) THEN
                   ! 
                   DO k = 1,nsp
                      chdr(k,js) =chdr(k,js)+ ylh(k,lh,nd)*chlhdr(jr,lh,js)
                      chdrr(k,js)=chdrr(k,js)+ylh(k,lh,nd)*chlhdrr(jr,lh,js)
                   ENDDO

                   DO k = 1,nsp
                      chdrt(k,js)=chdrt(k,js)+ylht(k,lh,nd)*chlhdr(jr,lh,js)
                      chdrf(k,js)=chdrf(k,js)+ylhf(k,lh,nd)*chlhdr(jr,lh,js)
                      chdt(k,js) =chdt(k,js) +ylht(k,lh,nd)*chlh(jr,lh,js)
                      chdf(k,js) =chdf(k,js) +ylhf(k,lh,nd)*chlh(jr,lh,js)
                      chdtt(k,js)=chdtt(k,js)+ylhtt(k,lh,nd)*chlh(jr,lh,js)
                      chdff(k,js)=chdff(k,js)+ylhff(k,lh,nd)*chlh(jr,lh,js)
                      chdtf(k,js)=chdtf(k,js)+ylhtf(k,lh,nd)*chlh(jr,lh,js)
                   ENDDO


                ENDIF

             ENDDO ! lh
          ENDDO   ! js

          IF (l_gga) THEN
             CALL mkgylm(&
                  &                  input%jspins,atoms%rmsh(jr,n),thet,nsp,dimension%nspd,dimension%jspd,ch,chdr,&
                  &                  chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf,&
                  &                  agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,&
                  &                  grgru,grgrd,gzgr)!keep
          ELSE
             agr(:)   = 0.0 ; agru(:)  = 0.0 ; agrd(:)  = 0.0 
             g2r(:)   = 0.0 ; g2ru(:)  = 0.0 ; g2rd(:)  = 0.0 
             gggr(:)  = 0.0 ; gggru(:) = 0.0 ; gggrd(:) = 0.0 
             grgru(:) = 0.0 ; grgrd(:) = 0.0 ; gzgr(:)  = 0.0 
          ENDIF

          !
          !         rhmnm: rho_minimum_muffin-tin..

          rhmnm=10.e+10

          DO js=1,input%jspins
             DO i=1,nsp
                ch(i,js) = max(ch(i,js),d_15)
                rhmnm = min(rhmnm,ch(i,js))
             ENDDO
          ENDDO

#ifdef CPP_MPI 
          IF (rhmnm.LT.rhmn_local) THEN
             !$OMP ATOMIC WRITE
             rhmn_local = rhmnm
          ENDIF
#else
          IF (rhmnm.LT.rhmn) THEN
             !$OMP ATOMIC WRITE
             rhmn = rhmnm
             ichsmrg = 1
          ENDIF
#endif

          IF (rhmn.LT.obsolete%chng) THEN
             WRITE (6,'(/'' rhmn.lt.obsolete%chng in vmtxc. rhmn,obsolete%chng='',&
                  &        2d9.2)') rhmn,obsolete%chng
             !             CALL juDFT_error("vmtxcg: rhmn.lt.chng",calledby="vmtxcg")
          ENDIF

          !
          !         calculate the ex.-cor. potential

         
          CALL vxcallg(&
               &                 xcpot,lwbc,input%jspins,nsp,nsp,ch,agr,agru,agrd,&
               &                 g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,&
               &                 vx,vxc)

          IF (xcpot%lda_atom(n)) THEN
             ! Use local part of pw91 for this atom
             CALL vxcallg(&
                  xcpot_tmp,lwbc,input%jspins,nsp,nsp,ch,agr,agru,agrd,&
                  g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,&
                  vx,vxcl)
             !Mix the potentials
             divi = 1.0 / (atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(1,n))
             vxc(:,:) = ( vxcl(:,:) * ( atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(jr,n) ) +&
                     vxc(:,:) * ( atoms%rmsh(jr,n) - atoms%rmsh(1,n) ) ) * divi
          ENDIF
 

          IF (mpi%irank == 0) THEN
          IF (mod(jr,1000).eq.0)&
               &              WRITE (6,'(/'' 999vxc''/(10d15.7))')&
               &                    ((vxc(k,js),k=1,nsp),js=1,input%jspins)
          ENDIF !irank==0


          !         now determine the corresponding potential number

          DO js = 1,input%jspins
             !
             !           multiplicate vx/vxc with the weights of the k-points
             !
             DO k = 1,nsp
                vx (k,js) = vx (k,js)*wt(k)
                vxc(k,js) = vxc(k,js)*wt(k)
             ENDDO

             IF (mpi%irank == 0) THEN
                IF (mod(jr,1500).EQ.0)&
                     &        WRITE (6,'('' 999wt''/(10d15.7))') (wt(k),k=1,nsp)
                IF (mod(jr,1500).EQ.0)&
                     &        WRITE (6,'('' 999vxc''/(10d15.7))') (vxc(k,js),k=1,nsp)
             ENDIF !irank==0

             DO lh = 0,sphhar%nlh(nd)
                !
                ! --->        determine the corresponding potential number
                !c            through gauss integration
                !
                vlh=dot_product(vxc(:nsp,js),ylh(:nsp,lh,nd))
#ifdef CPP_MPI             
                vr_local(jr,lh,n,js) = vr_local(jr,lh,n,js) + vlh
#else
                vr(jr,lh,n,js) = vr(jr,lh,n,js) + vlh
#endif

                ! --->        add to the given potential

                IF (mod(jr,1500).EQ.0) THEN
                   WRITE(6,'('' 999lh,js,jr,n,vlh='',4i5,d15.7)')&
                        &            lh,js,jr,n,vlh
                   WRITE(6,'('' 9ylh''/(10d15.7))') (ylh(k,lh,nd),k=1,nsp)
                ENDIF

                vlh=dot_product(vx(:nsp,js),ylh(:nsp,lh,nd))

#ifdef CPP_MPI             
                vxr_local(jr,lh,n,js) = vxr_local(jr,lh,n,js) + vlh
#else
                vxr(jr,lh,n,js) = vxr(jr,lh,n,js) + vlh
#endif

             ENDDO ! lh
          ENDDO   ! js

          IF (input%total) then
             !
             !           calculate the ex.-cor energy density
             !
             CALL excallg(xcpot,lwbc,input%jspins,nsp,&
                  ch,agr,agru,agrd,g2r,g2ru,g2rd,&
                  gggr,gggru,gggrd,gzgr, exc)
             
             IF (xcpot%lda_atom(n)) THEN
             ! Use local part of pw91 for this atom
                CALL excallg(xcpot_tmp,lwbc,input%jspins,nsp,&
                     ch,agr,agru,agrd,g2r,g2ru,g2rd,&
                     gggr,gggru,gggrd,gzgr, excl)
                !Mix the potentials
                exc(:) = ( excl(:) * ( atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(jr,n) ) +&
                     exc(:) * ( atoms%rmsh(jr,n) - atoms%rmsh(1,n) ) ) * divi
             ENDIF
             IF (mpi%irank == 0) THEN
             IF (mod(jr,10000).EQ.0)&
                  &        WRITE (6,'(/'' 999exc''/(10d15.7))') (exc(k),k=1,nsp)
             ENDIF !irank==0

          ENDIF

          !         now determine the corresponding energy density number
          !
          !         multiplicate exc with the weights of the k-points
          !
          DO k = 1,nsp
             exc(k) = exc(k)*wt(k)
          ENDDO

          DO lh = 0,sphhar%nlh(nd)
             !
             !           determine the corresponding potential number through gauss
             !           integration
             !
             elh = 0.0
             DO k = 1,nsp
                elh = elh + exc(k)*ylh(k,lh,nd)
             ENDDO
#ifdef CPP_MPI
             excr_local(jr,lh,n) =  elh
#else
             excr(jr,lh,n) =  elh
#endif

          ENDDO

190    ENDDO
       !$OMP END DO
       DEALLOCATE (ch,chdr,chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf)
       !$OMP END PARALLEL 

       !        WRITE(6,'(/'' n='',i3/'' 9999vr''/(10d15.7))') n,
       !     &   (((vr(jr,lh,n,js),jr=1,jri(n),100),lh=0,ntypsy(nat)),js=1,jspins)
       !        WRITE(6,'(/'' 9999excr''/(10d15.7))')
       !     &   ((excr(jr,lh,n),jr=1,jri(n),100),lh=0,ntypsy(nat))

200 ENDDO
    DEALLOCATE (chlh,chlhdr,chlhdrr)

#ifdef CPP_MPI
    CALL MPI_ALLREDUCE(vxr_local,vxr,atoms%jmtd*(1+sphhar%nlhd)*atoms%ntype*dimension%jspd,CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)     !ToDo:CPP_MPI_REAL?
    !using vxr_local as a temporal buffer
    CALL MPI_ALLREDUCE(vr_local,vxr_local,atoms%jmtd*(1+sphhar%nlhd)*atoms%ntype*dimension%jspd,CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)    
    vr = vr + vxr_local
    CALL MPI_ALLREDUCE(excr_local,excr,atoms%jmtd*(1+sphhar%nlhd)*atoms%ntype,CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)    
     CALL MPI_ALLREDUCE(rhmn_local,rhmn_reduced,1,CPP_MPI_REAL,MPI_MIN,mpi%mpi_comm,ierr)
     IF (rhmn_reduced.LT.rhmn) THEN
             rhmn = rhmn_reduced
             ichsmrg = 1
     ENDIF
#endif
    !

    RETURN
  END SUBROUTINE vmtxcg
END MODULE m_vmtxcg
