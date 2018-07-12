!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vmtxcg
  !.....------------------------------------------------------------------
  !     Calculate the GGA xc-potential in the MT-spheres
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
  SUBROUTINE vmtxcg(dimension,mpi,sphhar,atoms,&
       den,xcpot,input,sym, obsolete,vxc,vx,exc)

#include"cpp_double.h"
   
    USE m_grdchlh
    USE m_mkgylm
    USE m_gaussp
    USE m_types_xcpot_inbuild
    USE m_types
    IMPLICIT NONE

    CLASS(t_xcpot),INTENT(IN)      :: xcpot
    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_mpi),INTENT(IN)         :: mpi
    TYPE(t_obsolete),INTENT(IN)    :: obsolete
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_sphhar),INTENT(IN)      :: sphhar
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_potden),INTENT(IN)      :: den
    TYPE(t_potden),INTENT(INOUT)   :: vxc,vx,exc
#ifdef CPP_MPI
    include "mpif.h"
#endif
    !     ..
    !     .. Local Scalars ..
    TYPE(t_gradients)::grad
    INTEGER jr,js,k,lh,n,nd,ist,nsp,ixpm ,i,nat
    REAL    rhmnm,d_15,elh,vlh
    LOGICAL lwbc              ! if true, white-bird trick
    !     ..
    !     .. Local Arrays ..
    REAL v_x(dimension%nspd,dimension%jspd),v_xc(dimension%nspd,dimension%jspd),e_xc(dimension%nspd),rx(3,dimension%nspd)
    REAL vxcl(DIMENSION%nspd,DIMENSION%jspd),excl(DIMENSION%nspd),divi
    TYPE(t_xcpot_inbuild)::xcpot_tmp
    REAL wt(dimension%nspd),rr2(atoms%jmtd),thet(dimension%nspd)
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
#endif    
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC max,mod,min
    !ta+
    !.....------------------------------------------------------------------

#ifdef CPP_MPI
    CALL MPI_BCAST(obsolete%ndvgrd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    vr_local = 0.d0
    vxr_local = 0.d0
    excr_local = 0.d0
#endif
  
    lwbc=.false.

    d_15 = 1.e-15
    ixpm = 1
    ist  = 1

    !     generates nspd points on a sherical shell with radius 1.0
    !     angular mesh equidistant in phi,
    !     theta are zeros of the legendre polynomials
    !
    CALL gaussp(atoms%lmaxd, rx,wt)

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
       SELECT TYPE(xcpot)
       TYPE IS(t_xcpot_inbuild)
          IF (xcpot%lda_atom(n))THEN
             IF((.NOT.xcpot%is_name("pw91"))) &
                  CALL judft_warn("Using locally LDA only possible with pw91 functional")
             CALL xcpot_tmp%init("l91",.FALSE.,atoms%ntype)
          ENDIF
       END SELECT
       
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
                chlh(jr,lh,js) = den%mt(jr,lh,n,js)*rr2(jr)
             ENDDO

             CALL grdchlh(ixpm,ist,atoms%jri(n),atoms%dx(n),atoms%rmsh(1,n),&
                  chlh(1,lh,js),obsolete%ndvgrd, chlhdr(1,lh,js),chlhdrr(1,lh,js))
       
          ENDDO ! js
       ENDDO   ! lh
       
       !
       !-->    loop over radial mesh 
       !
       !!$OMP PARALLEL DEFAULT(none) &
#ifdef CPP_MPI
       !!$OMP& SHARED(vr_local,vxr_local,excr_local) &
#endif
       !!$OMP& SHARED(vxc,vx,exc) &
       !!$OMP& SHARED(dimension,mpi,sphhar,atoms,den,xcpot,input,sym,obsolete)&
       !!$OMP& SHARED(n,nd,ist,ixpm,nsp,nat,d_15,lwbc) &
       !!$OMP& SHARED(rx,wt,rr2,thet,xcpot_tmp) &
       !!$OMP& SHARED(ylh,ylht,ylhtt,ylhf,ylhff,ylhtf) &
       !!$OMP& SHARED(chlh,chlhdr,chlhdrr) &
       !!$OMP& SHARED(ierr,n_start,n_stride) &
       !!$OMP& PRIVATE(js,k,lh,i,rhmnm,elh,vlh) &
       !!$OMP& PRIVATE(v_x,v_xc,e_xc,excl,vxcl,divi) &
       !!$OMP& PRIVATE(grad) &
       !!$OMP& PRIVATE(ch,chdr,chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf)
       ALLOCATE ( ch(DIMENSION%nspd,DIMENSION%jspd),chdr(DIMENSION%nspd,DIMENSION%jspd),&
            chdt(DIMENSION%nspd,DIMENSION%jspd), chdf(DIMENSION%nspd,DIMENSION%jspd),&
            chdrr(DIMENSION%nspd,DIMENSION%jspd),chdtt(DIMENSION%nspd,DIMENSION%jspd), &
            chdff(DIMENSION%nspd,DIMENSION%jspd),chdtf(DIMENSION%nspd,DIMENSION%jspd),&
            chdrt(DIMENSION%nspd,DIMENSION%jspd), chdrf(DIMENSION%nspd,DIMENSION%jspd))
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


             ENDDO ! lh
          ENDDO   ! js

          CALL xcpot%alloc_gradients(nsp,input%jspins,grad)
          CALL mkgylm(input%jspins,atoms%rmsh(jr,n),thet,nsp,DIMENSION%nspd,DIMENSION%jspd,ch,chdr,&
               chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf,grad)

          !Set charge to minimum value
          ch(:,:)=MAX(ch(:,:),d_15)
          !
          !         calculate the ex.-cor. potential

          CALL xcpot%get_vxc(input%jspins,ch(:nsp,:),v_xc,v_x,grad)
          SELECT TYPE(xcpot)
          TYPE IS (t_xcpot_inbuild)
             IF (xcpot%lda_atom(n)) THEN
                ! Use local part of pw91 for this atom
                CALL xcpot_tmp%get_vxc(input%jspins,ch(:nsp,:),vxcl,v_x,grad)
                !Mix the potentials
                divi = 1.0 / (atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(1,n))
                v_xc(:,:) = ( vxcl(:,:) * ( atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(jr,n) ) +&
                     v_xc(:,:) * ( atoms%rmsh(jr,n) - atoms%rmsh(1,n) ) ) * divi
             ENDIF
          END SELECT

          IF (mpi%irank == 0) THEN
          IF (mod(jr,1000).eq.0) WRITE (6,'(/'' 999vxc''/(10d15.7))')&
               ((v_xc(k,js),k=1,nsp),js=1,input%jspins)
          ENDIF !irank==0


          !         now determine the corresponding potential number

          DO js = 1,input%jspins
             !
             !           multiplicate vx/vxc with the weights of the k-points
             !
             DO k = 1,nsp
                v_x (k,js) = v_x (k,js)*wt(k)
                v_xc(k,js) = v_xc(k,js)*wt(k)
             ENDDO

             IF (mpi%irank == 0) THEN
                IF (mod(jr,1500).EQ.0) WRITE (6,'('' 999wt''/(10d15.7))') (wt(k),k=1,nsp)
                IF (mod(jr,1500).EQ.0) WRITE (6,'('' 999vxc''/(10d15.7))') (v_xc(k,js),k=1,nsp)
             ENDIF !irank==0

             DO lh = 0,sphhar%nlh(nd)
                !
                ! --->        determine the corresponding potential number
                !c            through gauss integration
                !
                vlh=dot_product(v_xc(:nsp,js),ylh(:nsp,lh,nd))
#ifdef CPP_MPI             
                vr_local(jr,lh,n,js) = vr_local(jr,lh,n,js) + vlh
#else
                vxc%mt(jr,lh,n,js) = vxc%mt(jr,lh,n,js) + vlh
#endif

                ! --->        add to the given potential

                IF (mod(jr,1500).EQ.0) THEN
                   WRITE(6,'('' 999lh,js,jr,n,vlh='',4i5,d15.7)')&
                        &            lh,js,jr,n,vlh
                   WRITE(6,'('' 9ylh''/(10d15.7))') (ylh(k,lh,nd),k=1,nsp)
                ENDIF

                vlh=dot_product(v_x(:nsp,js),ylh(:nsp,lh,nd))

#ifdef CPP_MPI             
                vxr_local(jr,lh,n,js) = vxr_local(jr,lh,n,js) + vlh
#else
                vx%mt(jr,lh,n,js) = vx%mt(jr,lh,n,js) + vlh
#endif

             ENDDO ! lh
          ENDDO   ! js

          IF (ALLOCATED(exc%mt)) THEN
             !
             !           calculate the ex.-cor energy density
             !
             CALL xcpot%get_exc(input%jspins,ch(:nsp,:),e_xc,grad)
             SELECT TYPE(xcpot)
             TYPE IS(t_xcpot_inbuild)
                IF (xcpot%lda_atom(n)) THEN
                   ! Use local part of pw91 for this atom
                   CALL xcpot_tmp%get_exc(input%jspins,ch(:nsp,:),excl,grad)
                   !Mix the potentials
                   e_xc(:) = ( excl(:) * ( atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(jr,n) ) +&
                        e_xc(:) * ( atoms%rmsh(jr,n) - atoms%rmsh(1,n) ) ) * divi
                ENDIF
             END SELECT
             IF (mpi%irank == 0) THEN
             IF (mod(jr,10000).EQ.0) WRITE (6,'(/'' 999exc''/(10d15.7))') (e_xc(k),k=1,nsp)
             ENDIF !irank==0


          !         now determine the corresponding energy density number
          !
          !         multiplicate exc with the weights of the k-points
          !
          DO k = 1,nsp
             e_xc(k) = e_xc(k)*wt(k)
          ENDDO

          DO lh = 0,sphhar%nlh(nd)
             !
             !           determine the corresponding potential number through gauss
             !           integration
             !
             elh = 0.0
             DO k = 1,nsp
                elh = elh + e_xc(k)*ylh(k,lh,nd)
             ENDDO
#ifdef CPP_MPI
             excr_local(jr,lh,n) =  elh
#else
             exc%mt(jr,lh,n,1) =  elh
#endif

          ENDDO
          ENDIF

190    ENDDO
       !!$OMP END DO
       DEALLOCATE (ch,chdr,chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf)
       !!$OMP END PARALLEL 

       !        WRITE(6,'(/'' n='',i3/'' 9999vr''/(10d15.7))') n,
       !     &   (((vr(jr,lh,n,js),jr=1,jri(n),100),lh=0,ntypsy(nat)),js=1,jspins)
       !        WRITE(6,'(/'' 9999excr''/(10d15.7))')
       !     &   ((excr(jr,lh,n),jr=1,jri(n),100),lh=0,ntypsy(nat))

200 ENDDO
    DEALLOCATE (chlh,chlhdr,chlhdrr)

#ifdef CPP_MPI
    CALL MPI_ALLREDUCE(vxr_local,vx%mt,atoms%jmtd*(1+sphhar%nlhd)*atoms%ntype*dimension%jspd,CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)     !ToDo:CPP_MPI_REAL?
    !using vxr_local as a temporal buffer
    CALL MPI_ALLREDUCE(vr_local,vxr_local,atoms%jmtd*(1+sphhar%nlhd)*atoms%ntype*dimension%jspd,CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)    
    vxc%mt = vxc%mt + vxr_local
    CALL MPI_ALLREDUCE(excr_local,exc%mt(:,:,:,1),atoms%jmtd*(1+sphhar%nlhd)*atoms%ntype,CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)    
#endif
    !

    RETURN
  END SUBROUTINE vmtxcg
END MODULE m_vmtxcg
