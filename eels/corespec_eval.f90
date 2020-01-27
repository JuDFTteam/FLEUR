!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_corespec_eval

  USE m_types_setup
  USE m_types_usdus       
  USE m_types_cdnval, ONLY: t_eigVecCoeffs        
  USE m_constants        
  USE m_corespec

  IMPLICIT NONE

  CONTAINS

!===============================================================================
!
!  S U B R O U T I N E   C O R E S P E C _ G A U N T
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE corespec_gaunt()

!    use factorials
    
    use m_clebsch

    implicit none

!    real :: threejsymbol

    logical :: cmsum,clevn,ctiq1,ctiq2,ctiq3
    real :: twol1p1,twola1p1,twolip1

    smeno = "corespec_gaunt"

    write(*,'(/,a)') trim(smeno)//ssep

!    call init_factorials(6*(lmaxd+1)+1)

    ln = min(0,minval(csv%lc)-1)
    lx = max(csi%lx,maxval(csv%lc)+1)

    lan = 0
    lax = csi%lx+maxval(csv%lc)+1

    lin = minval(csv%lc)-1
    lix = maxval(csv%lc)+1

!!$    print*,"ln,lx,lan,lax,lin,lix"
!!$    print*,ln,lx,lan,lax,lin,lix

    if(.not.allocated(csv%gaunt)) &
         &allocate(csv%gaunt(ln:lx,-lx:lx,lan:lax,-lax:lax,lin:lix,-lix:lix))
    csv%gaunt = 0.0

! m<=l condition fulfilled by looping m within l value interval {-l,...,+l}

    csv%gaunt = 0.0
    do l1 = ln,lx
      do m1 = -l1,l1
        do la1 = lan,lax
          do mu1 = -la1,la1
            do li = lin,lix
              do mi = -li,li
                cmsum = (m1+mu1-mi).eq.0  ! sum of m q-nos. = 0
                clevn = mod((l1+la1+li),2).eq.0  ! sum of l q-nos. is even
                ctiq1 = (la1+li-l1).ge.0  ! triangle inequality 1
                ctiq2 = (l1+li-la1).ge.0  ! triangle inequality 2
                ctiq3 = (l1+la1-li).ge.0  ! triangle inequality 3
                twol1p1 = dble(2*l1+1)
                twola1p1 = dble(2*la1+1)
                twolip1 = dble(2*li+1)
                if(cmsum.and.clevn.and.ctiq1.and.ctiq2.and.ctiq3) then
                  csv%gaunt(l1,m1,la1,mu1,li,mi) = &
!                       &threejsymbol((l1),(la1),0,0,(li),0)*&
!                       &threejsymbol((l1),(la1),(m1),(mu1),(li),-(mi)))&
                       &clebsch(real(l1),real(la1),0.0,0.0,real(li),0.0)*&
                       &clebsch(real(l1),real(la1),real(m1),real(mu1),real(li),real(mi))*&
                       &sqrt(twol1p1*twola1p1/(4.0*pi_const*twolip1))*&
                       &(-1)**(mi)
                  if(csv%gaunt(l1,m1,la1,mu1,li,mi).ne.0.0) &
                       &write(53,'(6i5,f12.6)') l1,m1,la1,mu1,li,-mi,csv%gaunt(l1,m1,la1,mu1,li,mi)
!!$                  if(abs(csv%gaunt(l1,m1,la1,mu1,li,mi)).lt.1.d-6) &
!!$                       &write(*,'(6i5,f24.20)') l1,m1,la1,mu1,li,-mi,csv%gaunt(l1,m1,la1,mu1,li,mi)
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    if(csi%verb.eq.1) write(*,*) ""

  end subroutine corespec_gaunt
!
!===============================================================================
!===============================================================================
!
!  S U B R O U T I N E   C O R E S P E C _ R M E
!
!-------------------------------------------------------------------------------
!
  subroutine corespec_rme(atoms,input,itype,nstd,&
                          jspins,jspin,efermi,&
                          msh,vr,f,g)

    USE m_constants, ONLY : c_light
    USE m_setcor
    USE m_differ
    USE m_intgr, ONLY : intgr3
    USE m_dr2fdr
    USE m_sphbes
    USE m_intgr, ONLY : intgr3

    implicit none

    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_input),INTENT(IN)   :: input

    integer, intent(in) :: itype  ! call in ntype loop with itype = n
    integer, intent(in) :: nstd
    integer, intent(in) :: jspins,jspin
    real, intent(in) :: efermi
    integer, intent(in) :: msh
    real, intent    (in) :: vr(atoms%jmtd,atoms%ntype,jspins)
    real, intent (in) :: f(atoms%jmtd,2,0:atoms%lmaxd,jspin:jspin)
    real, intent (in) :: g(atoms%jmtd,2,0:atoms%lmaxd,jspin:jspin)

    integer :: nr,lx,lax,lin,lix,nqv,nen,nex

    integer :: ir,id,iljc,ic,il,ila,iqv,ie,ierr
    integer :: nst,kappa(nstd),nprnc(nstd)
    real :: nc,nlc,njc
    real :: c,bmu,t2,weight,e,d,rn,res,qr
    real :: vrd(msh),occ(nstd,jspins),a(msh),b(msh)
    real :: resd

    real, allocatable :: fpd(:)
    real, allocatable :: fp(:),fc(:),fsb(:)
    real :: sum1,sum2,sum3,sum1d,sum2d

    smeno = "corespec_rme"

    if(itype.ne.csi%atomType) return

    write(*,'(/,a)') trim(smeno)//ssep

    c = c_light(1.0)

    nr = atoms%jri(itype)

    allocate(fp(nr),fpd(nr),fc(nr))

  ! CORE functions
  ! csv%fc(ir,:,:,:) : ir = 1:nr
  ! csv%fc(:,id,:,:) : id = 1 { r*fc(r) } or 2 { r*[dfc(r)/dr] }
  ! csv%fc(:,:,iljc,:) : iljc = 1:csv%nljc
  ! csv%fc(:,:,:,ic) : ic = 1 { large component } or 2 { small component }
    if(.not.allocated(csv%fc)) allocate(csv%fc(nr,2,csv%nljc,2))
    csv%fc = 0.0

  ! core setup
    bmu = 0.0

    CALL setcor(itype,jspins,atoms,input,bmu,nst,kappa,nprnc,occ)

  ! extend core potential
    vrd(1:nr) = vr(1:nr,itype,jspin)
    t2 = vrd(nr)/(nr-msh)
    do ir = nr+1,msh
      vrd(ir) = vrd(nr)+t2*(ir-nr)
    enddo

  ! calculate core radial functions
    nc = real(csv%nc)
    do iljc = 1,csv%nljc
      njc = real(edgej(csi%edgeidx(iljc)))/2.0
      nlc = real(edgel(csi%edgeidx(iljc)))
      weight = 2*njc+1.0
      csv%eedge(iljc) = -2*(atoms%zatom(itype)/(nc+nlc))**2
      d = exp(atoms%dx(itype))
      rn = atoms%rmsh(1,itype)*(d**(msh-1))

      CALL differ(nc,nlc,njc,c,atoms%zatom(itype),atoms%dx(itype),&
                  atoms%rmsh(1,itype),rn,d,msh,vrd,&
                  e,&
                  a,b,ierr)

      csv%eedge(iljc)=dble(e)
      csv%fc(:,1,iljc,1) = a(1:nr)  ! large component
      csv%fc(:,1,iljc,2) = b(1:nr)  ! small component
      do ic = 1,2
        fp(:) = real(csv%fc(:,1,iljc,ic)*atoms%rmsh(1:nr,itype))
        CALL dr2fdr(fp,atoms%rmsh(1,itype),nr,fc)
        csv%fc(:,2,iljc,ic)=dble(fc(:)/atoms%rmsh(1:nr,itype))

        if(ic.eq.1) then
        do ir=1,nr
        write(90,'(2i5,16e12.4)') iljc,ir,atoms%rmsh(ir,itype),csv%fc(ir,1,iljc,ic),csv%fc(ir,2,iljc,ic)
        enddo
        write(90,*) ''
        write(90,*) ''
        endif

      enddo
      
      fp = csv%fc(:,1,iljc,1)**2
      CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,sum1)
      fp = csv%fc(:,1,iljc,2)**2
      CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,sum2)
      write(*,'(a,i5,3f8.4)') "ui",0,sum1,sum2,sum1+sum2

      fp = csv%fc(:,2,iljc,1)**2
      CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,sum1)
      fp = csv%fc(:,2,iljc,2)**2
      CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,sum2)
      write(*,'(a,i5,3f8.4)') "ui",0,sum1,sum2,sum1+sum2

      write(60,*) ""
      csv%occ(iljc) = dble(occ((csv%nc-1)**2+csi%edgeidx(iljc),jspin))
      write(*,"(a,2(a,i2),a,f3.1,2(a,i2),a,f16.8,a)") trim(smeno)//ssep,&
           &"core state: iljc = ",iljc,&
           &", nc = ",nint(nc),&
           &", njc = ",njc,&
           &", nlc = ",nint(nlc),&
           &", occ. csv%occ = ",nint(csv%occ(iljc)),&
           &", energy csv%eedge(iljc) = ",csv%eedge(iljc)," Ha found"
      if(efermi-csv%eedge(iljc).lt.ecoredeep) then
        write(*,csmsgsfs)  trim(smeno),&
             &"core state energy found not very deep: ",&
             &"efermi-csv%eedge(iljc) = ",&
             &(efermi-csv%eedge(iljc))*hartree_to_ev_const,"eV ; are you sure ? "//csmsgwar
      endif
    enddo

    CALL corespec_eloss_qv(efermi)  ! set-up csv%eloss and csv%qv arrays

    lx = csi%lx  ! lmax for l index
    lax = lx+maxval(csv%lc)+1  ! lmax for la index
    lin = minval(csv%lc)  ! minimum lc q-no.
    lix = maxval(csv%lc)  ! maximum lc q-no.
    nqv = csv%nqv
    nen = csv%nen
    nex = csv%nex

    allocate(fsb(0:lax))

  ! VALENCE functions
  ! csv%fv(ir,:,:,:) : ir = 1:nr
  ! csv%fv(:,il,:,:) : il = 0:csi%lx
  ! csv%fv(:,:,id,:) : id = 1 { a.u } or 2 { b.u' }
  ! csv%fv(:,:,:,ic) : ic = 1 { large component } or 2 { small component }
    if(.not.allocated(csv%fv)) allocate(csv%fv(nr,0:lx,2,2))
    csv%fv = 0.0

    do ic = 1,2
      do il = 0,lx
        csv%fv(:,il,1,ic) = f(1:nr,ic,il,jspin)
        csv%fv(:,il,2,ic) = g(1:nr,ic,il,jspin)

        if(ic.eq.1) then
        do ir=1,nr
        write(70,'(3i5,16e12.4)') ic,il,ir,atoms%rmsh(ir,itype),csv%fv(ir,il,1,ic),csv%fv(ir,il,2,ic)
        enddo
        write(70,*) ''
        write(70,*) ''
        endif

       fp(:) = csv%fv(:,il,1,ic)**2
        CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,sum1)
        fp(:) = csv%fv(:,il,2,ic)**2
        CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,sum2)
        fp(:) = csv%fv(:,il,1,ic)*csv%fv(:,il,2,ic)
        CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,sum3)
        write(*,'(a,i5,3f8.4)') "u ",il,sum1,sum2,sum3

      enddo
    enddo

  ! BESSEL functions
  ! csv%fb(ir,:,:,:,:) : ir = 1:nr
  ! csv%fb(:,il,:,:,:) : il = 0:lax
  ! csv%fb(:,,:,iljc,:,:) : iljc = 1:csv%nljc
  ! csv%fb(:,:,:,iqv,:) : iqv = 1:nqv
  ! csv%fb(:,:,:,:,ie) : ie = nen:nex
    if(.not.allocated(csv%fb)) allocate(csv%fb(nr,0:lax,csv%nljc,nqv,nen:nex))
    csv%fb = 0.0

    do ie = nen,nex
      do iqv = 1,nqv
        do iljc = 1,csv%nljc
          do ir = 1,nr
            fsb=0.0
            qr = real(csv%qv(0,iljc,iqv,ie)*atoms%rmsh(ir,itype))
            CALL sphbes(lax,qr,fsb)
            csv%fb(ir,:,iljc,iqv,ie) = dble(fsb)
!            write(70,'(4i5,16e12.4)') ie,iqv,iljc,ir,atoms%rmsh(ir,itype),fsb
          enddo
!          write(70,*) ''
        enddo
      enddo
    enddo

    if(.NOT.ALLOCATED(csv%rmeA)) THEN
       ALLOCATE(csv%rmeA(2,0:lx,0:lax,csv%nljc,2,nqv,nen:nex))
       ALLOCATE(csv%rmeB(2,0:lx,0:lax,csv%nljc,2,nqv,nen:nex))
       ALLOCATE(csv%rmeC(2,0:lx,0:lax,csv%nljc,2,nqv,nen:nex))
    END IF
    csv%rmeA = 0.0
    csv%rmeB = 0.0
    csv%rmeC = 0.0

    do ie = nen,nex
      do iqv = 1,nqv
        do ic = 1,2
          do iljc = 1,csv%nljc
            do ila = 0,lax
              do il = 0,lx
                do id = 1,2
                  fp(:)=csv%fc(1:nr,1,iljc,ic)*&
                       &csv%fv(1:nr,il,id,ic)*&
                       &csv%fb(1:nr,ila,iljc,iqv,ie)
                  CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,res)
                  csv%rmeA(id,il,ila,iljc,ic,iqv,ie)=dble(res)
                  fp(:)=fp(:)/atoms%rmsh(1:nr,itype)
                  CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,res)
                  csv%rmeC(id,il,ila,iljc,ic,iqv,ie)=dble(res)
                  fp(:)=csv%fc(1:nr,2,iljc,ic)*&
                       &csv%fv(1:nr,il,id,ic)*&
                       &csv%fb(1:nr,ila,iljc,iqv,ie)!/atoms%rmsh(1:nr,itype)
                  CALL intgr3(fp,atoms%rmsh(1,itype),atoms%dx(itype),nr,res)
                  csv%rmeB(id,il,ila,iljc,ic,iqv,ie)=dble(res)
                  write(41,'(7(a,i5),3f12.6)') 'ie=',ie,' iqv=',iqv,' ic=',ic,&
                        ' iljc=',iljc,' id=',id,' ila=',ila,' il=',il,&
                        csv%rmeA(id,il,ila,iljc,ic,iqv,ie),&
                        csv%rmeB(id,il,ila,iljc,ic,iqv,ie),&
                        csv%rmeC(id,il,ila,iljc,ic,iqv,ie)
                enddo  ! id
              enddo  ! il
            enddo  ! ila
          enddo  ! iljc
        enddo  ! ic
      enddo  ! iqv
    enddo  ! ie

    print*,size(3*csv%rmeA)
 
    deallocate(fsb,fc,fpd,fp)

    if(csi%verb.eq.1) write(*,*) ""

  end subroutine corespec_rme
!
!===============================================================================
!===============================================================================
!
!  S U B R O U T I N E   C O R E S P E C _ D O S
!
!-------------------------------------------------------------------------------
!
  subroutine corespec_dos(atoms,usdus,ispin,lmd,nkpt,ikpt,&
                          neigd,noccbd,efermi,sig_dos,eig,we,eigVecCoeffs)

    IMPLICIT NONE

    TYPE (t_atoms), INTENT(IN)      :: atoms
    TYPE (t_usdus), INTENT(IN)      :: usdus
    TYPE(t_eigVecCoeffs),INTENT(IN) :: eigVecCoeffs

!     .. Scalar Arguments ..
    integer, intent(in) :: ispin,lmd,nkpt,ikpt
    integer, intent(in) :: neigd,noccbd
    real, intent(in) :: efermi,sig_dos
!     .. Array Arguments ..
    real, intent (in) :: eig(neigd),we(noccbd)

! local variables
    integer :: lx,lmx,nen,nex
    integer :: iatom,iband,l1,m1,l2,m2,lm1,lm2,ie!,ljc,iqv
    real :: sigma,eigos(noccbd)
    real :: sum11,sum22

    smeno = "corespec_dos"

    lx = csi%lx
    lmx = lx*(lx+2)
    nen = csv%nen
    nex = csv%nex
    iatom = atoms%neq(csi%atomType)
    sigma = sqrt(2.0)*sig_dos*hartree_to_ev_const
    sigma = sig_dos*hartree_to_ev_const
    eigos(1:noccbd) = (eig(1:noccbd)-efermi)*hartree_to_ev_const/dble(sigma)

    if(ikpt.eq.1) then
      write(*,'(/,a)') trim(smeno)//ssep
      if(.not.allocated(csv%dose)) allocate(csv%dose(2,2,0:lmx,0:lmx,0:nex))
      if(.not.allocated(csv%dosb)) allocate(csv%dosb(2,2,0:lmx,0:lmx,noccbd))
      if(.not.allocated(csv%eos)) then
        allocate(csv%eos(0:nex))
        csv%eos(:) = csv%egrid(:)/dble(sigma)
      endif
      csv%dose = 0.0
    endif
    csv%dosb = 0.0

    do iband = 1,noccbd
      do l1 = 0,lx
        do m1 = -l1,l1
          lm1 = l1*(l1+1)+m1
          do l2 = 0,lx!!$
            do m2 = -l2,l2!!$
              lm2 = l2*(l2+1)+m2!!$
!!!! for dose:
!!!! order of xcof, xcof' : aa', ab', ba', bb'
!!!! is meant by            11 , 12 , 21 , 22
!!!! or, put another way, first index is unprimed (i.e. the outer loop furter down), second index is primed (i.e. the inner loop further down)

!!!! Check what we(1) is and does, if necessary, add a we(1) contribution to all acofs and bcofs
          csv%dosb(1,1,lm2,lm1,iband) = dble(eigVecCoeffs%acof(iband,lm2,iatom,ispin)*&
               &conjg(eigVecCoeffs%acof(iband,lm1,iatom,ispin)))!*we(1)
          csv%dosb(1,2,lm2,lm1,iband) = dble(eigVecCoeffs%acof(iband,lm2,iatom,ispin)*&
               &conjg(eigVecCoeffs%bcof(iband,lm1,iatom,ispin)))
          csv%dosb(2,1,lm2,lm1,iband) = dble(eigVecCoeffs%bcof(iband,lm2,iatom,ispin)*&
               &conjg(eigVecCoeffs%acof(iband,lm1,iatom,ispin)))
          csv%dosb(2,2,lm2,lm1,iband) = dble(eigVecCoeffs%bcof(iband,lm2,iatom,ispin)*&
               &conjg(eigVecCoeffs%bcof(iband,lm1,iatom,ispin)))!*we(1)*usdus%ddn(l1,csi%atomType,ispin)
!!!!! this has to be checked: is >> ddn << factor necessary !!!!!
!!!!! Check if we(iband) should be multiplied with everything
        enddo!!$
        enddo!!$
        enddo
      enddo
      if(eigos(iband)+3.0*sigma.ge.csv%eos(0).and.&
           &eigos(iband)-3.0*sigma.le.csv%eos(nex)) then
        do ie = 0,nex
          csv%dose(:,:,:,:,ie) = csv%dose(:,:,:,:,ie)+&
               &csv%dosb(:,:,:,:,iband)*exp(-(eigos(iband)-csv%eos(ie))**2)
        enddo
      endif
    enddo

    if(ikpt.eq.nkpt) then
      csv%dose = csv%dose/(sqrt(pi_const)*sigma)
      do ie=0,nex
        write(36,*) csv%egrid(ie),sum(csv%dose(1,1,:,:,ie)+csv%dose(2,2,:,:,ie))
      enddo
      write(36,*) ""
      write(*,'(10i8)') atoms%llod,noccbd,atoms%nlod,atoms%nat,neigd,atoms%ntype,atoms%lmaxd
      write(*,'(10i8)') lmd,atoms%ntype

      if(csi%verb.eq.1) write(*,*) ""
    endif

  end subroutine corespec_dos
!
!===============================================================================
!===============================================================================
!
!  S U B R O U T I N E   C O R E S P E C _ D D S C S
!
!-------------------------------------------------------------------------------
!
  subroutine corespec_ddscs(jspin,jspins)

    use m_ylm


    implicit none

    integer, intent(in) :: jspin,jspins

    integer :: lx,lmx,lan,lax,nqv,nen,nex,nor

    integer :: ic,ie,iqv,ior,it,iljc,imi,id1,id2,ip1,ip2
    integer :: l1,l2,m1,m2,lm1,lm2
    integer :: la1,la2,mu1,mu2
    integer :: li,mi
    integer :: lamu,lamu1,lamu2

    real :: gamma,beta,rho,qepref
    real, allocatable :: orvec(:,:)
!    real, allocatable :: orw(:)
    real :: ga(0:2,2)
    real :: prd(0:2,0:2)
    complex :: td(2),orfac,ila1la2
    complex, allocatable :: tdy(:,:),orpref(:),ylm(:,:)

    smeno = "corespec_ddscs"

    write(*,'(/,a)') trim(smeno)//ssep

    lx = csi%lx
    lmx = lx*(lx+2)
    lan = 0
    lax = csi%lx+maxval(csv%lc)+1
    nqv = csv%nqv
    nen = csv%nen
    nex = csv%nex

    nor = 1
!    nor = 26
    if(.not.allocated(orvec)) allocate(orvec(1:nor,3))
!    if(.not.allocated(orw)) allocate(orw(0:nor))
!    call lebedev(nor,orvec,orw)
    orvec(1,:) = (/1.0,0.0,0.0/)

    if(.not.allocated(csv%ddscs)) then
      allocate(csv%ddscs(2,0:nor,1:csv%nljc,0:nqv,0:nex))
      csv%ddscs = cmplx(0.0,0.0)
    endif
    if(.not.allocated(tdy)) allocate(tdy(0:nor,2))
    if(.not.allocated(orpref)) then
      allocate(orpref(0:nor))
      orpref(0) = 1.0
      if(nor.gt.0) then
        orpref(1:nor) = (4.0*pi_const)**2
        if(.not.allocated(ylm)) allocate(ylm(0:lax*(lax+2),nor))
        do ior = 1,nor
          CALL ylm4(lax,orvec(ior,:),ylm(:,ior))
          do la1 = lan,lax ; do mu1 = -la1,la1
            lamu = la1*(la1+1)+mu1
            write(98,'(3i5,2f12.8)') la1,mu1,lamu,ylm(lamu,ior)
          enddo; enddo
        enddo
      endif
    endif

    ic = 1
    gamma = csv%gamma
    beta = csv%beta

    rho = alpha*beta*sqrt(4.0*pi_const/3.0)
    print*,gamma,beta,rho
!    rho = 0.0

    do ie = nen,nex  ! energy
      do iqv = 1,nqv  ! q-vector
        do iljc = 1,csv%nljc  ! core levels
          li = edgel(csi%edgeidx(iljc))
          qepref = 4.0*gamma**2*csv%qv1(iljc,iqv,ie)/csv%qv0/&
               &(csv%qv(0,iljc,iqv,ie)**2-(csv%eloss(iljc,ie)*alpha)**2)**2
!!$          write(*,'(2i5,3f20.4)') ie,iljc,csv%qv(0,iljc,iqv,ie),csv%eloss(iljc,ie)*alpha,qepref
          tdy = cmplx(0.0,0.0)

          do imi = 1,(edgej(csi%edgeidx(iljc))+1)/2!min(nint(csv%occ(iljc)*jspins/2),2*li+1)
            mi = sign(jspin)*(edgej(csi%edgeidx(iljc))-4*(imi-1)-1)/2
!!            print*,jspin,ie,iljc,li,mi
            write(39,*) jspin,ie,iljc,li,mi

            do l1 = 0,lx ; do m1 = -l1,l1
              lm1 = l1*(l1+1)+m1
              do la1 = lan,lax ; do mu1 = -la1,la1
                lamu1 = la1*(la1+1)+mu1

                ga(0,1) = csv%gaunt(l1,-m1,la1,mu1,li,mi)
                ga(1,1) = csv%gaunt(li+1,-mi,li,mi,1,0)*&
                     &csv%gaunt(l1,-m1,la1,mu1,li+1,mi)+&
                     &csv%gaunt(li-1,-mi,li,mi,1,0)*&
                     &csv%gaunt(l1,-m1,la1,mu1,li-1,mi)
                ga(2,1) = csv%gaunt(li+1,-mi,li,mi+1,1,-1)*&
                     &csv%gaunt(l1,-m1,la1,mu1,li+1,mi)+&
                     &csv%gaunt(li-1,-mi,li,mi+1,1,-1)*&
                     &csv%gaunt(l1,-m1,la1,mu1,li-1,mi)*&
                     &sqrt(dble(2*(li-mi)*(li+mi+1)))+mi*ga(1,1)

                do l2 = 0,lx ; do m2 = -l2,l2
                  lm2 = l2*(l2+1)+m2
                  do la2 = lan,lax ; do mu2 = -la2,la2
                    lamu2 = la2*(la2+1)+mu2
                    
 !                   if(l1.eq.l2.and.m1.eq.m2) then

                    ga(0,2) = csv%gaunt(l2,-m2,la2,mu2,li,mi)
                    ga(1,2) = csv%gaunt(li+1,-mi,li,mi,1,0)*&
                             &csv%gaunt(l2,-m2,la2,mu2,li+1,mi)+&
                             &csv%gaunt(li-1,-mi,li,mi,1,0)*&
                             &csv%gaunt(l2,-m2,la2,mu2,li-1,mi)
                    ga(2,2) = csv%gaunt(li+1,-mi,li,mi+1,1,-1)*&
                             &csv%gaunt(l2,-m2,la2,mu2,li+1,mi)+&
                             &csv%gaunt(li-1,-mi,li,mi+1,1,-1)*&
                             &csv%gaunt(l2,-m2,la2,mu2,li-1,mi)*&
                             &sqrt(dble(2*(li-mi)*(li+mi+1)))+mi*ga(1,2)

                    prd = 0.0

                    do id1 = 1,2 ; 
                      do id2 = 1,2
                        prd(0,0) = prd(0,0)+ &
                           csv%rmeA(id1,l1,la1,iljc,ic,iqv,ie)*csv%rmeA(id2,l2,la2,iljc,ic,iqv,ie)*csv%dose(id1,id2,lm1,lm2,ie)
                        prd(0,1) = prd(0,1)+ &
                           csv%rmeA(id1,l1,la1,iljc,ic,iqv,ie)*csv%rmeB(id2,l2,la2,iljc,ic,iqv,ie)*csv%dose(id1,id2,lm1,lm2,ie)
                        prd(0,2) = prd(0,2)+ &
                           csv%rmeA(id1,l1,la1,iljc,ic,iqv,ie)*csv%rmeC(id2,l2,la2,iljc,ic,iqv,ie)*csv%dose(id1,id2,lm1,lm2,ie)
                        prd(1,0) = prd(1,0)+ &
                           csv%rmeB(id1,l1,la1,iljc,ic,iqv,ie)*csv%rmeA(id2,l2,la2,iljc,ic,iqv,ie)*csv%dose(id1,id2,lm1,lm2,ie)
                        prd(1,1) = prd(1,1)+ &
                           csv%rmeB(id1,l1,la1,iljc,ic,iqv,ie)*csv%rmeB(id2,l2,la2,iljc,ic,iqv,ie)*csv%dose(id1,id2,lm1,lm2,ie)
                        prd(1,2) = prd(1,2)+ &
                           csv%rmeB(id1,l1,la1,iljc,ic,iqv,ie)*csv%rmeC(id2,l2,la2,iljc,ic,iqv,ie)*csv%dose(id1,id2,lm1,lm2,ie)
                        prd(2,0) = prd(2,0)+ &
                           csv%rmeC(id1,l1,la1,iljc,ic,iqv,ie)*csv%rmeA(id2,l2,la2,iljc,ic,iqv,ie)*csv%dose(id1,id2,lm1,lm2,ie)
                        prd(2,1) = prd(2,1)+ &
                           csv%rmeC(id1,l1,la1,iljc,ic,iqv,ie)*csv%rmeB(id2,l2,la2,iljc,ic,iqv,ie)*csv%dose(id1,id2,lm1,lm2,ie)
                        prd(2,2) = prd(2,2)+ &
                           csv%rmeC(id1,l1,la1,iljc,ic,iqv,ie)*csv%rmeC(id2,l2,la2,iljc,ic,iqv,ie)*csv%dose(id1,id2,lm1,lm2,ie)
                      enddo 
                    enddo

                    td(1) = prd(0,0)*ga(0,1)*ga(0,2)
                    td(2) = cone*rho**2*(&
                           &prd(1,1)*ga(1,1)*ga(1,2)&
                          &+prd(2,2)*ga(2,1)*ga(2,2)&
                          &-prd(1,2)*ga(1,1)*ga(2,2)&
                          &-prd(2,1)*ga(2,1)*ga(1,2))&
                          &+cimu*rho*(-1)**(li+1)*(&
                          &-prd(0,1)*ga(0,1)*ga(1,2)&
                          &+prd(0,2)*ga(0,1)*ga(2,2)&
                          &+prd(1,0)*ga(1,1)*ga(0,2)&
                          &-prd(2,0)*ga(2,1)*ga(0,2))&
                          &+td(1)
                    
                    ila1la2 = cimu**(la1-la2)
                    
                    if(abs(real(td(1))).gt.0.0.or.abs(real(td(2))).gt.0.0.or.abs(aimag(td(2))).gt.0.0) then
                    write(39,'(2f4.0,i2,6i4,a,6i4,a,6f7.3,a,4f10.6)') ila1la2,la1-la2,l1,-m1,la1,mu1,li,mi,'  ',l2,-m2,la2,mu2,li,mi,'  ',ga(0,1),ga(0,2),ga(1,1),ga(1,2),ga(2,1),ga(2,2),'  ',1000000*td
                    endif

                    do ior = 0,nor  ! orientation
                      if(ior.eq.0) then
                        orfac = cone
                      else
                        orfac = ylm(lamu1,ior)*conjg(ylm(lamu2,ior))
                      endif
                      
                      tdy(ior,1:2) = tdy(ior,1:2)+td(1:2)*orfac*ila1la2

                    enddo  ! ior

 !                 endif

                  enddo; enddo
                enddo; enddo

              enddo; enddo
            enddo; enddo

          enddo  ! mi

          do it = 1,2
            do ior = 0,nor
              csv%ddscs(it,ior,iljc,iqv,ie) = csv%ddscs(it,ior,iljc,iqv,ie)+&
                   &qepref*orpref(ior)*tdy(ior,it)
!!        calculate the integral over all q-vectors, save the result in iqv=0
          csv%ddscs(it,ior,iljc,0,ie) = csv%ddscs(it,ior,iljc,0,ie)+&
                   &csv%ddscs(it,ior,iljc,iqv,ie)*csv%qv(4,iljc,iqv,ie)
            enddo
          enddo

        enddo  ! iljc
      enddo  ! iqv
    enddo  ! ie

    if(jspin.eq.1) then
      do ior = 0,nor
      do iljc = 1,csv%nljc
        do ie = nen,nex
!!        write(37,'(2i5,f8.3,4es16.4)') ior,iljc,csv%eloss(iljc,ie)*hartree_to_ev_const,csv%ddscs(1,ior,iljc,1,ie),csv%ddscs(2,ior,iljc,1,ie)
          write(37,'(2i5,f16.3,4es16.4)') ior,iljc,csv%eloss(iljc,ie)*hartree_to_ev_const,csv%ddscs(1,ior,iljc,0,ie),csv%ddscs(2,ior,iljc,0,ie)
        enddo
        write(37,*) ""
      enddo
      write(37,*) ""
      enddo
    endif
    if(jspin.eq.2) then
      do ior = 0,nor
      do iljc = 1,csv%nljc
        do ie = nen,nex
!!        write(38,'(2i5,f8.3,4es16.4)') ior,iljc,csv%eloss(iljc,ie)*hartree_to_ev_const,csv%ddscs(1,ior,iljc,1,ie),csv%ddscs(2,ior,iljc,1,ie)
          write(38,'(2i5,f16.3,4es16.4)') ior,iljc,csv%eloss(iljc,ie)*hartree_to_ev_const,csv%ddscs(1,ior,iljc,0,ie),csv%ddscs(2,ior,iljc,0,ie)
        enddo
        write(38,*) ""
      enddo
      write(38,*) ""
      enddo
    endif

    if(csi%verb.eq.1) write(*,*) ""

  end subroutine corespec_ddscs
!
!===============================================================================
!===============================================================================
!
!  S U B R O U T I N E   C O R E S P E C _ E L O S S _ Q V
!
!-------------------------------------------------------------------------------
!
  subroutine corespec_eloss_qv(efermi)

    implicit none

    real, intent(in) :: efermi

    integer :: ie,iljc,iqv,iphi,ir
    real :: eout,relfac,pi,ri,r,dr,p,alpha,beta,geofac,gf1,gf2,normfac
    pi = 3.141592653589793238462643
    smeno = "corespec_eloss_qv"
    normfac = 1!4*csv%nqr**2!/(pi*(r**2))
    write(*,'(/,a)') trim(smeno)//ssep
!    csv%nqphi = 12
!    csv%nqr = 20
    csv%nqv = 1+csv%nqphi*csv%nqr
!    write(*,'(2i6,3f16.7)')csv%nqr,csv%nqphi,csv%alpha_ex,csv%beta_ex,csv%I0
    if(.not.allocated(csv%eloss)) &
         &allocate(csv%eloss(csv%nljc,csv%nen:csv%nex))
    if(.not.allocated(csv%qv1)) &
         &allocate(csv%qv1(csv%nljc,csv%nqv,csv%nen:csv%nex))
    do ie = csv%nen,csv%nex
      do iljc = 1,csv%nljc
        csv%eloss(iljc,ie) = csv%egrid(ie)/hartree_to_ev_const+dble(efermi)-csv%eedge(iljc)
!!$        print*,iljc,ie,csv%egrid(ie),csv%eloss(iljc,ie)
      enddo
    enddo

    csv%qv0 = e2q(csi%ek0/hartree_to_ev_const)
    relfac = (mec2)**2/(csi%ek0+mec2)**2
    alpha=csv%qv0*csv%alpha_ex
!!$    print*,csi%ek0,csv%qv0

    if(.not.allocated(csv%qv)) &
         &allocate(csv%qv(0:4,csv%nljc,csv%nqv,csv%nen:csv%nex))
!!  qv(0) = |qv(1:3)|
!!  qv(4) = weight of qv(1:3)
    csv%qv=0.0
    if(csv%nqv.gt.1)then    
    do ie = csv%nen,csv%nex
      do iljc = 1,csv%nljc
          eout = csi%ek0/hartree_to_ev_const-csv%eloss(iljc,ie)
        do iqv = 1,csv%nqv
          csv%qv1(iljc,iqv,ie) = e2q(eout)
!!   set up circular 2D mesh with qz==(csv%qv1(iljc,iqv,ie)-csv%qv0)*relfac
!!   and (qx,qy) = r*(sin(phi),cos(phi))
!!   R = alpha + beta
!!   r_0 = 0
!!   r_i = i/N_r*R , i>0
!!   phi_i = i/N_phi*2pi
!!   Areas of each volume element, this corresponds to the weights of the point in the integral
!!   A_1 = pi/4 * (r_0+r_1)²
!!   A_i = pi/(N_phi*4)*(r_{i+1}²-r_{i-1}²+2*r_i*(r_{i+1}-r_{i-1})) (1<i<=N_r)
!!   A_{N_r+1} = pi/(N_phi)*(r_{N_r}² - 1/4*(r_{N_r}+r_{N_r-1})²)
!!   A_ges = pi *(alpha + beta)²
!!   Numbering of the nodes (the i-index above are independent of each other, now we make a 2-D grid with 1-D indexing by counting upwards around the clock):
!!   j=1 => center node (r=0, phi=0)
!!   j=2 ... N_phi+1 => nodes of r=r_1, phi=phi_{i=mod_{N_phi}(j-1)}
!!   j...  => nodes of r=r_{1+frac{j-2-mod_{n_r}(j-2)}{n_r}} and phi = phi_{i=mod_{N_phi}(j-1)}

          beta=csv%beta_ex*csv%qv1(iljc,iqv,ie)
!!  r is the radius of the q-disc which sits at z=q_min and contains all the allowed q-vectors 
          r=alpha + beta !small angle approximation: sin(a) ~ a
          dr = r/csv%nqr
          iphi = modulo(iqv-1,csv%nqphi)
          ir = 1+(iqv-2-modulo(iqv-2,csv%nqr))/csv%nqr
          ri = (ir-0.5)*dr
!          normfac=normfac/(pi*(r**2))
!!        write the weight of qv, i.e. the area it represents
          csv%qv(4,iljc,iqv,ie) = 1.
!          write(*,'(6f16.10)')dr,csv%nqr,csv%nqphi,ir,r
!!        write weights and values of q_x and q_y for the q-vectors:
          if(ir.eq.0) then
          csv%qv(1,iljc,iqv,ie) = 0 ! here is the angular dependency
          csv%qv(2,iljc,iqv,ie) = 0 ! here is the angular dependency
             csv%qv(4,iljc,iqv,ie) = pi*0.0625*dr**2
          elseif(ir.eq.1) then
          csv%qv(1,iljc,iqv,ie) = ri*SIN(iphi/csv%nqphi*2*pi) ! here is the angular dependency
          csv%qv(2,iljc,iqv,ie) = ri*COS(iphi/csv%nqphi*2*pi) ! here is the angular dependency
             csv%qv(4,iljc,iqv,ie) = pi*0.9735*dr**2!!!!pi/csv%nqphi*(r**2-0.25*(2.*r-dr)**2) (old, less sensible mesh described above)
          else
          csv%qv(1,iljc,iqv,ie) = ri*SIN(iphi/csv%nqphi*2*pi) ! here is the angular dependency
          csv%qv(2,iljc,iqv,ie) = ri*COS(iphi/csv%nqphi*2*pi) ! here is the angular dependency
             csv%qv(4,iljc,iqv,ie) = pi/csv%nqphi*(2*ir-1)*dr**2
          endif
!!        write z coordinates:
          csv%qv(3,iljc,iqv,ie) = (csv%qv0-csv%qv1(iljc,iqv,ie))*relfac ! here is no angular dependency
!!        write the length of qv
          csv%qv(0,iljc,iqv,ie) = sqrt(&
               &dot_product(csv%qv(1:3,iljc,iqv,ie),csv%qv(1:3,iljc,iqv,ie)))
!          write(*,'(f16.6)')csv%qv(4,iljc,iqv,ie)

!!        calculate the g_alpha_beta function and multiply with the weight to obtain the overall weight of the specific q-vector
!!        Step 1: calculate all the relevant point of the overlapping circles
          p=0.5*(ri**2+(alpha)**2-(beta)**2)/(ri)
!          write(*,'(f16.10)')p
           geofac=0.
!          gf1=csv%I0/(csv%alpha_ex**2)*min(alpha,beta)**2
!          gf2=csv%I0/(csv%alpha_ex**2)*(0.5*pi*(alpha**2 + beta**2)-p*sqrt(alpha**2-p**2)-(ri-p)*sqrt(beta**2-(ri-p)**2)&
!                    &-beta**2*asin((ri-p)/beta)-alpha**2*asin(p/alpha))
          if(ri.LE.abs(alpha-beta)) then
              geofac=csv%I0/(csv%alpha_ex**2)*min(alpha,beta)**2
!             geofac=csv%I0/(csv%alpha_ex**2)*min(alpha,beta)**2
!             write(*,'(f16.6)')geofac
          elseif(ri.GE.(alpha+beta)) then
             geofac=0.
             write(*,csmsgsis)'geofac is 0'
          else
              geofac=csv%I0/(csv%alpha_ex**2)*(0.5*pi*(alpha**2 + beta**2)-p*sqrt(alpha**2-p**2)-(ri-p)*sqrt(beta**2-(ri-p)**2)&
                    &-beta**2*asin((ri-p)/beta)-alpha**2*asin(p/alpha))
!             geofac=csv%I0/(csv%alpha_ex**2)*(0.5*pi*(alpha**2 + beta**2)-p*sqrt(alpha**2-p**2)-(ri-p)*sqrt(beta**2-(ri-p)**2)&
!                    &-beta**2*asin((ri-p)/beta)-alpha**2*asin(p/alpha))
!             write(*,csmsgsis)'geofac is not 0'
!             write(*,'(3f16.6)')alpha**2-p**2,beta**2-(ri-p)**2, geofac
           endif
!          write(*,'(f16.10)')geofac
          csv%qv(4,iljc,iqv,ie) = csv%qv(4,iljc,iqv,ie)*geofac *normfac/(pi*(r**2))
!          csv%qv(4,iljc,iqv,ie) = 1.
!          write(*,'(f16.6)')csv%qv(4,iljc,iqv,ie)
!!$          write(*,'(3i5,2f16.2,6f16.6)') ie,iqv,iljc,csi%ek0,eout*hartree_to_ev_const,csv%qv1(iljc,iqv,ie),csv%eloss(iljc,ie),csv%qv(:,iljc,iqv,ie)
          write(*,'(5f16.5)')alpha,beta,ri,abs(alpha-beta),alpha+beta!,csv%nqr,csv%qv(4,iljc,iqv,ie)
        enddo
      enddo
    enddo
    else !number of q-vectors == 1:
    do ie = csv%nen,csv%nex
      do iljc = 1,csv%nljc
          eout = csi%ek0/hartree_to_ev_const-csv%eloss(iljc,ie)
        do iqv = 1,csv%nqv
          csv%qv1(iljc,iqv,ie) = e2q(eout)
          beta=csv%beta_ex*csv%qv1(iljc,iqv,ie)
          r=alpha + beta !small angle approximation: sin(a) ~ a
!!        only q||z vectors are calcualted:
!!        write x, y, and z coordinates:
          csv%qv(1,iljc,iqv,ie) = 0 ! here is no angular dependency
          csv%qv(2,iljc,iqv,ie) = 0 ! here is no angular dependency
 
          csv%qv(3,iljc,iqv,ie) = (csv%qv0-csv%qv1(iljc,iqv,ie))*relfac ! here is no angular dependency
!!        write the length of qv
          csv%qv(0,iljc,iqv,ie) = sqrt(&
               &dot_product(csv%qv(1:3,iljc,iqv,ie),csv%qv(1:3,iljc,iqv,ie)))
!!        write the weight of qv, i.e. the area it represents, normalized by the
!total area, i.e. 1.
          csv%qv(4,iljc,iqv,ie) = 1.!(pi*(r**2))!!*0.25
          dr = r
          write(*,'(7f16.5)')alpha,beta,r*500,pi,r**2,pi*dr*dr!,csv%nqr,csv%qv(4,iljc,iqv,ie)
!!        calculate the g_alpha_beta function and multiply with the weight to obtain the overall weight of the specific q-vector
!!        Step 1: calculate all the relevant point of the overlapping circles
           geofac=0.
              geofac=csv%I0/(csv%alpha_ex**2)*min(alpha,beta)**2
          csv%qv(4,iljc,iqv,ie) = csv%qv(4,iljc,iqv,ie)*geofac
!          csv%qv(4,iljc,iqv,ie) = 1.
!          write(*,'(f16.6)')csv%qv(4,iljc,iqv,ie)
!          write(*,'(3i5,2f16.2,6f16.6)') ie,iqv,iljc,csi%ek0,eout*hartree_to_ev_const,csv%qv1(iljc,iqv,ie),csv%eloss(iljc,ie),csv%qv(:,iljc,iqv,ie)
        enddo
      enddo
    enddo
    
    endif

    if(csi%verb.eq.1) write(*,*) ""

!    implicit none

!    real, intent(in) :: efermi

!    integer :: ie,iljc,iqv
!    real :: eout,relfac

!    smeno = "corespec_eloss_qv"

!    write(*,'(/,a)') trim(smeno)//ssep

!    csv%nqv = 1

!    if(.not.allocated(csv%eloss)) &
!         &allocate(csv%eloss(csv%nljc,csv%nen:csv%nex))
!    if(.not.allocated(csv%qv1)) &
!         &allocate(csv%qv1(csv%nljc,csv%nqv,csv%nen:csv%nex))
!    do ie = csv%nen,csv%nex
!      do iljc = 1,csv%nljc
!        csv%eloss(iljc,ie) = csv%egrid(ie)/hartree_to_ev_const+dble(efermi)-csv%eedge(iljc)
!!$        print*,iljc,ie,csv%egrid(ie),csv%eloss(iljc,ie)
!      enddo
!    enddo

!    csv%qv0 = e2q(csi%ek0/hartree_to_ev_const)
!    relfac = (mec2)**2/(csi%ek0+mec2)**2
!!$    print*,csi%ek0,csv%qv0

!    if(.not.allocated(csv%qv)) &
!         &allocate(csv%qv(0:3,csv%nljc,csv%nqv,csv%nen:csv%nex))
!    csv%qv=0.0
!    do ie = csv%nen,csv%nex
!      do iqv = 1,csv%nqv
!        do iljc = 1,csv%nljc
!          eout = csi%ek0/hartree_to_ev_const-csv%eloss(iljc,ie)
!          csv%qv1(iljc,iqv,ie) = e2q(eout)
!          csv%qv(3,iljc,iqv,ie) = (csv%qv1(iljc,iqv,ie)-csv%qv0)*relfac
!          csv%qv(0,iljc,iqv,ie) = sqrt(&
!               &dot_product(csv%qv(1:3,iljc,iqv,ie),csv%qv(1:3,iljc,iqv,ie)))
!!$          write(*,'(3i5,2f16.2,6f16.6)') ie,iqv,iljc,csi%ek0,eout*hartree_to_ev_const,csv%qv1(iljc,iqv,ie),csv%eloss(iljc,ie),csv%qv(:,iljc,iqv,ie)
!        enddo
!      enddo
!    enddo

!    if(csi%verb.eq.1) write(*,*) ""

  end subroutine corespec_eloss_qv
!
!===============================================================================
!===============================================================================
!  F U N C T I O N   E 2 Q
!-------------------------------------------------------------------------------
!
  real function e2q(e)

    use m_corespec, only : mec2,alpha
    implicit none
    real, intent(in) :: e

    e2q=sqrt(e**2+2.0*e*mec2/hartree_to_ev_const)*alpha

  end function e2q
!
!===============================================================================



!
!===============================================================================
!===============================================================================
!
!  S U B R O U T I N E   L E B E D E V
!
!-------------------------------------------------------------------------------
!

  subroutine lebedev(nleb,r2leb,wleb)
    implicit none
    integer, intent(in) :: nleb
    double precision, intent(out) :: r2leb(nleb,3),wleb(nleb)

    integer :: ileb,ctrln
    double precision :: vec(0)

    if(nleb.eq. 0006) call LD0006(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0014) call LD0014(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0026) call LD0026(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0038) call LD0038(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0050) call LD0050(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0074) call LD0074(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0086) call LD0086(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0110) call LD0110(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0146) call LD0146(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0170) call LD0170(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0194) call LD0194(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0230) call LD0230(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0266) call LD0266(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0302) call LD0302(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0350) call LD0350(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0434) call LD0434(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0590) call LD0590(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0770) call LD0770(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 0974) call LD0974(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 1202) call LD1202(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 1454) call LD1454(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 1730) call LD1730(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 2030) call LD2030(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 2354) call LD2354(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 2702) call LD2702(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 3074) call LD3074(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 3470) call LD3470(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 3890) call LD3890(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 4334) call LD4334(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 4802) call LD4802(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 5294) call LD5294(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)
    if(nleb.eq. 5810) call LD5810(r2leb(1,1),r2leb(1,2),r2leb(1,3),wleb,ctrln)

    write(*,'(i8)') nleb
    do ileb = 1,nleb
       write(*,'(4f12.6)') r2leb(ileb,1:3),wleb(ileb)
    enddo

  end subroutine lebedev
!
!===============================================================================


       subroutine gen_oh(code, num, x, y, z, w, a, b, v)
       implicit logical(a-z)
       double precision x(*),y(*),z(*),w(*)
       double precision a,b,v
       integer code
       integer num
       double precision c
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated from C to fortran77 by hand.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
!    
!      Given a point on a sphere (specified by a and b), generate all
!      the equivalent points under Oh symmetry, making grid points with
!      weight v.
!      The variable num is increased by the number of different points
!      generated.
!    
!      Depending on code, there are 6...48 different but equivalent
!      points.
!    
!      code=1:   (0,0,1) etc                                (  6 points)
!      code=2:   (0,a,a) etc, a=1/sqrt(2)                   ( 12 points)
!      code=3:   (a,a,a) etc, a=1/sqrt(3)                   (  8 points)
!      code=4:   (a,a,b) etc, b=sqrt(1-2 a^2)               ( 24 points)
!      code=5:   (a,b,0) etc, b=sqrt(1-a^2), a input        ( 24 points)
!      code=6:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input  ( 48 points)
!    
       goto (1,2,3,4,5,6) code
       write (6,*) 'Gen_Oh: Invalid Code'
       stop 
    1  continue
       a=1.0
       x(1) =  a
       y(1) =  0.0
       z(1) =  0.0
       w(1) =  v
       x(2) = -a
       y(2) =  0.0
       z(2) =  0.0
       w(2) =  v
       x(3) =  0.0
       y(3) =  a
       z(3) =  0.0
       w(3) =  v
       x(4) =  0.0
       y(4) = -a
       z(4) =  0.0
       w(4) =  v
       x(5) =  0.0
       y(5) =  0.0
       z(5) =  a
       w(5) =  v
       x(6) =  0.0
       y(6) =  0.0
       z(6) = -a
       w(6) =  v
       num=num+6
       return
!    
    2  continue
       a=sqrt(0.5)
       x( 1) =  0.0
       y( 1) =  a
       z( 1) =  a
       w( 1) =  v
       x( 2) =  0.0
       y( 2) = -a
       z( 2) =  a
       w( 2) =  v
       x( 3) =  0.0
       y( 3) =  a
       z( 3) = -a
       w( 3) =  v
       x( 4) =  0.0
       y( 4) = -a
       z( 4) = -a
       w( 4) =  v
       x( 5) =  a
       y( 5) =  0.0
       z( 5) =  a
       w( 5) =  v
       x( 6) = -a
       y( 6) =  0.0
       z( 6) =  a
       w( 6) =  v
       x( 7) =  a
       y( 7) =  0.0
       z( 7) = -a
       w( 7) =  v
       x( 8) = -a
       y( 8) =  0.0
       z( 8) = -a
       w( 8) =  v
       x( 9) =  a
       y( 9) =  a
       z( 9) =  0.0
       w( 9) =  v
       x(10) = -a
       y(10) =  a
       z(10) =  0.0
       w(10) =  v
       x(11) =  a
       y(11) = -a
       z(11) =  0.0
       w(11) =  v
       x(12) = -a
       y(12) = -a
       z(12) =  0.0
       w(12) =  v
       num=num+12
       return
!    
    3  continue
       a = sqrt(1.0/3.0)
       x(1) =  a
       y(1) =  a
       z(1) =  a
       w(1) =  v
       x(2) = -a
       y(2) =  a
       z(2) =  a
       w(2) =  v
       x(3) =  a
       y(3) = -a
       z(3) =  a
       w(3) =  v
       x(4) = -a
       y(4) = -a
       z(4) =  a
       w(4) =  v
       x(5) =  a
       y(5) =  a
       z(5) = -a
       w(5) =  v
       x(6) = -a
       y(6) =  a
       z(6) = -a
       w(6) =  v
       x(7) =  a
       y(7) = -a
       z(7) = -a
       w(7) =  v
       x(8) = -a
       y(8) = -a
       z(8) = -a
       w(8) =  v
       num=num+8
       return
!    
    4  continue
       b = sqrt(1.0 - 2.0*a*a)
       x( 1) =  a
       y( 1) =  a
       z( 1) =  b
       w( 1) =  v
       x( 2) = -a
       y( 2) =  a
       z( 2) =  b
       w( 2) =  v
       x( 3) =  a
       y( 3) = -a
       z( 3) =  b
       w( 3) =  v
       x( 4) = -a
       y( 4) = -a
       z( 4) =  b
       w( 4) =  v
       x( 5) =  a
       y( 5) =  a
       z( 5) = -b
       w( 5) =  v
       x( 6) = -a
       y( 6) =  a
       z( 6) = -b
       w( 6) =  v
       x( 7) =  a
       y( 7) = -a
       z( 7) = -b
       w( 7) =  v
       x( 8) = -a
       y( 8) = -a
       z( 8) = -b
       w( 8) =  v
       x( 9) =  a
       y( 9) =  b
       z( 9) =  a
       w( 9) =  v
       x(10) = -a
       y(10) =  b
       z(10) =  a
       w(10) =  v
       x(11) =  a
       y(11) = -b
       z(11) =  a
       w(11) =  v
       x(12) = -a
       y(12) = -b
       z(12) =  a
       w(12) =  v
       x(13) =  a
       y(13) =  b
       z(13) = -a
       w(13) =  v
       x(14) = -a
       y(14) =  b
       z(14) = -a
       w(14) =  v
       x(15) =  a
       y(15) = -b
       z(15) = -a
       w(15) =  v
       x(16) = -a
       y(16) = -b
       z(16) = -a
       w(16) =  v
       x(17) =  b
       y(17) =  a
       z(17) =  a
       w(17) =  v
       x(18) = -b
       y(18) =  a
       z(18) =  a
       w(18) =  v
       x(19) =  b
       y(19) = -a
       z(19) =  a
       w(19) =  v
       x(20) = -b
       y(20) = -a
       z(20) =  a
       w(20) =  v
       x(21) =  b
       y(21) =  a
       z(21) = -a
       w(21) =  v
       x(22) = -b
       y(22) =  a
       z(22) = -a
       w(22) =  v
       x(23) =  b
       y(23) = -a
       z(23) = -a
       w(23) =  v
       x(24) = -b
       y(24) = -a
       z(24) = -a
       w(24) =  v
       num=num+24
       return
!    
    5  continue
       b=sqrt(1.0-a*a)
       x( 1) =  a
       y( 1) =  b
       z( 1) =  0.0
       w( 1) =  v
       x( 2) = -a
       y( 2) =  b
       z( 2) =  0.0
       w( 2) =  v
       x( 3) =  a
       y( 3) = -b
       z( 3) =  0.0
       w( 3) =  v
       x( 4) = -a
       y( 4) = -b
       z( 4) =  0.0
       w( 4) =  v
       x( 5) =  b
       y( 5) =  a
       z( 5) =  0.0
       w( 5) =  v
       x( 6) = -b
       y( 6) =  a
       z( 6) =  0.0
       w( 6) =  v
       x( 7) =  b
       y( 7) = -a
       z( 7) =  0.0
       w( 7) =  v
       x( 8) = -b
       y( 8) = -a
       z( 8) =  0.0
       w( 8) =  v
       x( 9) =  a
       y( 9) =  0.0
       z( 9) =  b
       w( 9) =  v
       x(10) = -a
       y(10) =  0.0
       z(10) =  b
       w(10) =  v
       x(11) =  a
       y(11) =  0.0
       z(11) = -b
       w(11) =  v
       x(12) = -a
       y(12) =  0.0
       z(12) = -b
       w(12) =  v
       x(13) =  b
       y(13) =  0.0
       z(13) =  a
       w(13) =  v
       x(14) = -b
       y(14) =  0.0
       z(14) =  a
       w(14) =  v
       x(15) =  b
       y(15) =  0.0
       z(15) = -a
       w(15) =  v
       x(16) = -b
       y(16) =  0.0
       z(16) = -a
       w(16) =  v
       x(17) =  0.0
       y(17) =  a
       z(17) =  b
       w(17) =  v
       x(18) =  0.0
       y(18) = -a
       z(18) =  b
       w(18) =  v
       x(19) =  0.0
       y(19) =  a
       z(19) = -b
       w(19) =  v
       x(20) =  0.0
       y(20) = -a
       z(20) = -b
       w(20) =  v
       x(21) =  0.0
       y(21) =  b
       z(21) =  a
       w(21) =  v
       x(22) =  0.0
       y(22) = -b
       z(22) =  a
       w(22) =  v
       x(23) =  0.0
       y(23) =  b
       z(23) = -a
       w(23) =  v
       x(24) =  0.0
       y(24) = -b
       z(24) = -a
       w(24) =  v
       num=num+24
       return
!    
    6  continue
       c=sqrt(1.0 - a*a - b*b)
       x( 1) =  a
       y( 1) =  b
       z( 1) =  c
       w( 1) =  v
       x( 2) = -a
       y( 2) =  b
       z( 2) =  c
       w( 2) =  v
       x( 3) =  a
       y( 3) = -b
       z( 3) =  c
       w( 3) =  v
       x( 4) = -a
       y( 4) = -b
       z( 4) =  c
       w( 4) =  v
       x( 5) =  a
       y( 5) =  b
       z( 5) = -c
       w( 5) =  v
       x( 6) = -a
       y( 6) =  b
       z( 6) = -c
       w( 6) =  v
       x( 7) =  a
       y( 7) = -b
       z( 7) = -c
       w( 7) =  v
       x( 8) = -a
       y( 8) = -b
       z( 8) = -c
       w( 8) =  v
       x( 9) =  a
       y( 9) =  c
       z( 9) =  b
       w( 9) =  v
       x(10) = -a
       y(10) =  c
       z(10) =  b
       w(10) =  v
       x(11) =  a
       y(11) = -c
       z(11) =  b
       w(11) =  v
       x(12) = -a
       y(12) = -c
       z(12) =  b
       w(12) =  v
       x(13) =  a
       y(13) =  c
       z(13) = -b
       w(13) =  v
       x(14) = -a
       y(14) =  c
       z(14) = -b
       w(14) =  v
       x(15) =  a
       y(15) = -c
       z(15) = -b
       w(15) =  v
       x(16) = -a
       y(16) = -c
       z(16) = -b
       w(16) =  v
       x(17) =  b
       y(17) =  a
       z(17) =  c
       w(17) =  v
       x(18) = -b
       y(18) =  a
       z(18) =  c
       w(18) =  v
       x(19) =  b
       y(19) = -a
       z(19) =  c
       w(19) =  v
       x(20) = -b
       y(20) = -a
       z(20) =  c
       w(20) =  v
       x(21) =  b
       y(21) =  a
       z(21) = -c
       w(21) =  v
       x(22) = -b
       y(22) =  a
       z(22) = -c
       w(22) =  v
       x(23) =  b
       y(23) = -a
       z(23) = -c
       w(23) =  v
       x(24) = -b
       y(24) = -a
       z(24) = -c
       w(24) =  v
       x(25) =  b
       y(25) =  c
       z(25) =  a
       w(25) =  v
       x(26) = -b
       y(26) =  c
       z(26) =  a
       w(26) =  v
       x(27) =  b
       y(27) = -c
       z(27) =  a
       w(27) =  v
       x(28) = -b
       y(28) = -c
       z(28) =  a
       w(28) =  v
       x(29) =  b
       y(29) =  c
       z(29) = -a
       w(29) =  v
       x(30) = -b
       y(30) =  c
       z(30) = -a
       w(30) =  v
       x(31) =  b
       y(31) = -c
       z(31) = -a
       w(31) =  v
       x(32) = -b
       y(32) = -c
       z(32) = -a
       w(32) =  v
       x(33) =  c
       y(33) =  a
       z(33) =  b
       w(33) =  v
       x(34) = -c
       y(34) =  a
       z(34) =  b
       w(34) =  v
       x(35) =  c
       y(35) = -a
       z(35) =  b
       w(35) =  v
       x(36) = -c
       y(36) = -a
       z(36) =  b
       w(36) =  v
       x(37) =  c
       y(37) =  a
       z(37) = -b
       w(37) =  v
       x(38) = -c
       y(38) =  a
       z(38) = -b
       w(38) =  v
       x(39) =  c
       y(39) = -a
       z(39) = -b
       w(39) =  v
       x(40) = -c
       y(40) = -a
       z(40) = -b
       w(40) =  v
       x(41) =  c
       y(41) =  b
       z(41) =  a
       w(41) =  v
       x(42) = -c
       y(42) =  b
       z(42) =  a
       w(42) =  v
       x(43) =  c
       y(43) = -b
       z(43) =  a
       w(43) =  v
       x(44) = -c
       y(44) = -b
       z(44) =  a
       w(44) =  v
       x(45) =  c
       y(45) =  b
       z(45) = -a
       w(45) =  v
       x(46) = -c
       y(46) =  b
       z(46) = -a
       w(46) =  v
       x(47) =  c
       y(47) = -b
       z(47) = -a
       w(47) =  v
       x(48) = -c
       y(48) = -b
       z(48) = -a
       w(48) =  v
       num=num+48
       return
       end
       SUBROUTINE LD0006(X,Y,Z,W,N)
       DOUBLE PRECISION X(   6)
       DOUBLE PRECISION Y(   6)
       DOUBLE PRECISION Z(   6)
       DOUBLE PRECISION W(   6)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV    6-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.1666666666666667
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0014(X,Y,Z,W,N)
       DOUBLE PRECISION X(  14)
       DOUBLE PRECISION Y(  14)
       DOUBLE PRECISION Z(  14)
       DOUBLE PRECISION W(  14)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV   14-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.6666666666666667e-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7500000000000000e-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0026(X,Y,Z,W,N)
       DOUBLE PRECISION X(  26)
       DOUBLE PRECISION Y(  26)
       DOUBLE PRECISION Z(  26)
       DOUBLE PRECISION W(  26)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV   26-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.4761904761904762e-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3809523809523810e-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714e-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0038(X,Y,Z,W,N)
       DOUBLE PRECISION X(  38)
       DOUBLE PRECISION Y(  38)
       DOUBLE PRECISION Z(  38)
       DOUBLE PRECISION W(  38)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV   38-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.9523809523809524e-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714e-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4597008433809831
       V=0.2857142857142857e-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0050(X,Y,Z,W,N)
       DOUBLE PRECISION X(  50)
       DOUBLE PRECISION Y(  50)
       DOUBLE PRECISION Z(  50)
       DOUBLE PRECISION W(  50)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV   50-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.1269841269841270e-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2257495590828924e-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2109375000000000e-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3015113445777636
       V=0.2017333553791887e-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0074(X,Y,Z,W,N)
       DOUBLE PRECISION X(  74)
       DOUBLE PRECISION Y(  74)
       DOUBLE PRECISION Z(  74)
       DOUBLE PRECISION W(  74)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV   74-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.5130671797338464e-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1660406956574204e-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2958603896103896e-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4803844614152614
       V=0.2657620708215946e-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3207726489807764
       V=0.1652217099371571e-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0086(X,Y,Z,W,N)
       DOUBLE PRECISION X(  86)
       DOUBLE PRECISION Y(  86)
       DOUBLE PRECISION Z(  86)
       DOUBLE PRECISION W(  86)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV   86-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.1154401154401154e-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1194390908585628e-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3696028464541502
       V=0.1111055571060340e-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6943540066026664
       V=0.1187650129453714e-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3742430390903412
       V=0.1181230374690448e-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0110(X,Y,Z,W,N)
       DOUBLE PRECISION X( 110)
       DOUBLE PRECISION Y( 110)
       DOUBLE PRECISION Z( 110)
       DOUBLE PRECISION W( 110)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  110-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.3828270494937162e-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9793737512487512e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1851156353447362
       V=0.8211737283191111e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6904210483822922
       V=0.9942814891178103e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3956894730559419
       V=0.9595471336070963e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4783690288121502
       V=0.9694996361663028e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0146(X,Y,Z,W,N)
       DOUBLE PRECISION X( 146)
       DOUBLE PRECISION Y( 146)
       DOUBLE PRECISION Z( 146)
       DOUBLE PRECISION W( 146)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  146-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.5996313688621381e-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7372999718620756e-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7210515360144488e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6764410400114264
       V=0.7116355493117555e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4174961227965453
       V=0.6753829486314477e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1574676672039082
       V=0.7574394159054034e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1403553811713183
       B=0.4493328323269557
       V=0.6991087353303262e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0170(X,Y,Z,W,N)
       DOUBLE PRECISION X( 170)
       DOUBLE PRECISION Y( 170)
       DOUBLE PRECISION Z( 170)
       DOUBLE PRECISION W( 170)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  170-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.5544842902037365e-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6071332770670752e-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6383674773515093e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2551252621114134
       V=0.5183387587747790e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6743601460362766
       V=0.6317929009813725e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318910696719410
       V=0.6201670006589077e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2613931360335988
       V=0.5477143385137348e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4990453161796037
       B=0.1446630744325115
       V=0.5968383987681156e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0194(X,Y,Z,W,N)
       DOUBLE PRECISION X( 194)
       DOUBLE PRECISION Y( 194)
       DOUBLE PRECISION Z( 194)
       DOUBLE PRECISION W( 194)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  194-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.1782340447244611e-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5716905949977102e-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5573383178848738e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6712973442695226
       V=0.5608704082587997e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2892465627575439
       V=0.5158237711805383e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4446933178717437
       V=0.5518771467273614e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1299335447650067
       V=0.4106777028169394e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3457702197611283
       V=0.5051846064614808e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1590417105383530
       B=0.8360360154824589
       V=0.5530248916233094e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0230(X,Y,Z,W,N)
       DOUBLE PRECISION X( 230)
       DOUBLE PRECISION Y( 230)
       DOUBLE PRECISION Z( 230)
       DOUBLE PRECISION W( 230)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  230-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=-0.5522639919727325e-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4450274607445226e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4492044687397611
       V=0.4496841067921404e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2520419490210201
       V=0.5049153450478750e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6981906658447242
       V=0.3976408018051883e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6587405243460960
       V=0.4401400650381014e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4038544050097660e-1
       V=0.1724544350544401e-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5823842309715585
       V=0.4231083095357343e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3545877390518688
       V=0.5198069864064399e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2272181808998187
       B=0.4864661535886647
       V=0.4695720972568883e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0266(X,Y,Z,W,N)
       DOUBLE PRECISION X( 266)
       DOUBLE PRECISION Y( 266)
       DOUBLE PRECISION Z( 266)
       DOUBLE PRECISION W( 266)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  266-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=-0.1313769127326952e-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2522728704859336e-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4186853881700583e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7039373391585475
       V=0.5315167977810885e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1012526248572414
       V=0.4047142377086219e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4647448726420539
       V=0.4112482394406990e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3277420654971629
       V=0.3595584899758782e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6620338663699974
       V=0.4256131351428158e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8506508083520399
       V=0.4229582700647240e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3233484542692899
       B=0.1153112011009701
       V=0.4080914225780505e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2314790158712601
       B=0.5244939240922365
       V=0.4071467593830964e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0302(X,Y,Z,W,N)
       DOUBLE PRECISION X( 302)
       DOUBLE PRECISION Y( 302)
       DOUBLE PRECISION Z( 302)
       DOUBLE PRECISION W( 302)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  302-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.8545911725128148e-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3599119285025571e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3515640345570105
       V=0.3449788424305883e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6566329410219612
       V=0.3604822601419882e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4729054132581005
       V=0.3576729661743367e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9618308522614784e-1
       V=0.2352101413689164e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2219645236294178
       V=0.3108953122413675e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7011766416089545
       V=0.3650045807677255e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2644152887060663
       V=0.2982344963171804e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5718955891878961
       V=0.3600820932216460e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2510034751770465
       B=0.8000727494073952
       V=0.3571540554273387e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1233548532583327
       B=0.4127724083168531
       V=0.3392312205006170e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0350(X,Y,Z,W,N)
       DOUBLE PRECISION X( 350)
       DOUBLE PRECISION Y( 350)
       DOUBLE PRECISION Z( 350)
       DOUBLE PRECISION W( 350)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  350-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.3006796749453936e-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3050627745650771e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7068965463912316
       V=0.1621104600288991e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4794682625712025
       V=0.3005701484901752e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1927533154878019
       V=0.2990992529653774e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6930357961327123
       V=0.2982170644107595e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3608302115520091
       V=0.2721564237310992e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6498486161496169
       V=0.3033513795811141e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1932945013230339
       V=0.3007949555218533e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3800494919899303
       V=0.2881964603055307e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2899558825499574
       B=0.7934537856582316
       V=0.2958357626535696e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9684121455103957e-1
       B=0.8280801506686862
       V=0.3036020026407088e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1833434647041659
       B=0.9074658265305127
       V=0.2832187403926303e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0434(X,Y,Z,W,N)
       DOUBLE PRECISION X( 434)
       DOUBLE PRECISION Y( 434)
       DOUBLE PRECISION Z( 434)
       DOUBLE PRECISION W( 434)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  434-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.5265897968224436e-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2548219972002607e-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2512317418927307e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6909346307509111
       V=0.2530403801186355e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1774836054609158
       V=0.2014279020918528e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4914342637784746
       V=0.2501725168402936e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6456664707424256
       V=0.2513267174597564e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2861289010307638
       V=0.2302694782227416e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7568084367178018e-1
       V=0.1462495621594614e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3927259763368002
       V=0.2445373437312980e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8818132877794288
       V=0.2417442375638981e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9776428111182649
       V=0.1910951282179532e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2054823696403044
       B=0.8689460322872412
       V=0.2416930044324775e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5905157048925271
       B=0.7999278543857286
       V=0.2512236854563495e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5550152361076807
       B=0.7717462626915901
       V=0.2496644054553086e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9371809858553722
       B=0.3344363145343455
       V=0.2236607760437849e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0590(X,Y,Z,W,N)
       DOUBLE PRECISION X( 590)
       DOUBLE PRECISION Y( 590)
       DOUBLE PRECISION Z( 590)
       DOUBLE PRECISION W( 590)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  590-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.3095121295306187e-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1852379698597489e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7040954938227469
       V=0.1871790639277744e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6807744066455243
       V=0.1858812585438317e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6372546939258752
       V=0.1852028828296213e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5044419707800358
       V=0.1846715956151242e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4215761784010967
       V=0.1818471778162769e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3317920736472123
       V=0.1749564657281154e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2384736701421887
       V=0.1617210647254411e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1459036449157763
       V=0.1384737234851692e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6095034115507196e-1
       V=0.9764331165051050e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6116843442009876
       V=0.1857161196774078e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3964755348199858
       V=0.1705153996395864e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1724782009907724
       V=0.1300321685886048e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5610263808622060
       B=0.3518280927733519
       V=0.1842866472905286e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4742392842551980
       B=0.2634716655937950
       V=0.1802658934377451e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5984126497885380
       B=0.1816640840360209
       V=0.1849830560443660e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3791035407695563
       B=0.1720795225656878
       V=0.1713904507106709e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2778673190586244
       B=0.8213021581932511e-1
       V=0.1555213603396808e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5033564271075117
       B=0.8999205842074875e-1
       V=0.1802239128008525e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0770(X,Y,Z,W,N)
       DOUBLE PRECISION X( 770)
       DOUBLE PRECISION Y( 770)
       DOUBLE PRECISION Z( 770)
       DOUBLE PRECISION W( 770)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  770-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.2192942088181184e-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1436433617319080e-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1421940344335877e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5087204410502360e-1
       V=0.6798123511050502e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1228198790178831
       V=0.9913184235294912e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2026890814408786
       V=0.1180207833238949e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2847745156464294
       V=0.1296599602080921e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3656719078978026
       V=0.1365871427428316e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4428264886713469
       V=0.1402988604775325e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5140619627249735
       V=0.1418645563595609e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6306401219166803
       V=0.1421376741851662e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6716883332022612
       V=0.1423996475490962e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6979792685336881
       V=0.1431554042178567e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1446865674195309
       V=0.9254401499865368e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3390263475411216
       V=0.1250239995053509e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5335804651263506
       V=0.1394365843329230e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6944024393349413e-1
       B=0.2355187894242326
       V=0.1127089094671749e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2269004109529460
       B=0.4102182474045730
       V=0.1345753760910670e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8025574607775339e-1
       B=0.6214302417481605
       V=0.1424957283316783e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1467999527896572
       B=0.3245284345717394
       V=0.1261523341237750e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1571507769824727
       B=0.5224482189696630
       V=0.1392547106052696e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2365702993157246
       B=0.6017546634089558
       V=0.1418761677877656e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7714815866765732e-1
       B=0.4346575516141163
       V=0.1338366684479554e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3062936666210730
       B=0.4908826589037616
       V=0.1393700862676131e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3822477379524787
       B=0.5648768149099500
       V=0.1415914757466932e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0974(X,Y,Z,W,N)
       DOUBLE PRECISION X( 974)
       DOUBLE PRECISION Y( 974)
       DOUBLE PRECISION Z( 974)
       DOUBLE PRECISION W( 974)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV  974-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.1438294190527431e-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1125772288287004e-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4292963545341347e-1
       V=0.4948029341949241e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1051426854086404
       V=0.7357990109125470e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1750024867623087
       V=0.8889132771304384e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2477653379650257
       V=0.9888347838921435e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3206567123955957
       V=0.1053299681709471e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3916520749849983
       V=0.1092778807014578e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4590825874187624
       V=0.1114389394063227e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5214563888415861
       V=0.1123724788051555e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6253170244654199
       V=0.1125239325243814e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6637926744523170
       V=0.1126153271815905e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6910410398498301
       V=0.1130286931123841e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7052907007457760
       V=0.1134986534363955e-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1236686762657990
       V=0.6823367927109931e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2940777114468387
       V=0.9454158160447096e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4697753849207649
       V=0.1074429975385679e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6334563241139567
       V=0.1129300086569132e-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5974048614181342e-1
       B=0.2029128752777523
       V=0.8436884500901954e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1375760408473636
       B=0.4602621942484054
       V=0.1075255720448885e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3391016526336286
       B=0.5030673999662036
       V=0.1108577236864462e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1271675191439820
       B=0.2817606422442134
       V=0.9566475323783357e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2693120740413512
       B=0.4331561291720157
       V=0.1080663250717391e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1419786452601918
       B=0.6256167358580814
       V=0.1126797131196295e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6709284600738255e-1
       B=0.3798395216859157
       V=0.1022568715358061e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7057738183256172e-1
       B=0.5517505421423520
       V=0.1108960267713108e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2783888477882155
       B=0.6029619156159187
       V=0.1122790653435766e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1979578938917407
       B=0.3589606329589096
       V=0.1032401847117460e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2087307061103274
       B=0.5348666438135476
       V=0.1107249382283854e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4055122137872836
       B=0.5674997546074373
       V=0.1121780048519972e-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD1202(X,Y,Z,W,N)
       DOUBLE PRECISION X(1202)
       DOUBLE PRECISION Y(1202)
       DOUBLE PRECISION Z(1202)
       DOUBLE PRECISION W(1202)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 1202-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.1105189233267572e-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9205232738090741e-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9133159786443561e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3712636449657089e-1
       V=0.3690421898017899e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9140060412262223e-1
       V=0.5603990928680660e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1531077852469906
       V=0.6865297629282609e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2180928891660612
       V=0.7720338551145630e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2839874532200175
       V=0.8301545958894795e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3491177600963764
       V=0.8686692550179628e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4121431461444309
       V=0.8927076285846890e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4718993627149127
       V=0.9060820238568219e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5273145452842337
       V=0.9119777254940867e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6209475332444019
       V=0.9128720138604181e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6569722711857291
       V=0.9130714935691735e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6841788309070143
       V=0.9152873784554116e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012604330123631
       V=0.9187436274321654e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1072382215478166
       V=0.5176977312965694e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2582068959496968
       V=0.7331143682101417e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4172752955306717
       V=0.8463232836379928e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5700366911792503
       V=0.9031122694253992e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9827986018263947
       B=0.1771774022615325
       V=0.6485778453163257e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9624249230326228
       B=0.2475716463426288
       V=0.7435030910982369e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9402007994128811
       B=0.3354616289066489
       V=0.7998527891839054e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9320822040143202
       B=0.3173615246611977
       V=0.8101731497468018e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9043674199393299
       B=0.4090268427085357
       V=0.8483389574594331e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8912407560074747
       B=0.3854291150669224
       V=0.8556299257311812e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8676435628462708
       B=0.4932221184851285
       V=0.8803208679738260e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8581979986041619
       B=0.4785320675922435
       V=0.8811048182425720e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8396753624049856
       B=0.4507422593157064
       V=0.8850282341265444e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8165288564022188
       B=0.5632123020762100
       V=0.9021342299040653e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8015469370783529
       B=0.5434303569693900
       V=0.9010091677105086e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7773563069070351
       B=0.5123518486419871
       V=0.9022692938426915e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7661621213900394
       B=0.6394279634749102
       V=0.9158016174693465e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7553584143533510
       B=0.6269805509024392
       V=0.9131578003189435e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7344305757559503
       B=0.6031161693096310
       V=0.9107813579482705e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7043837184021765
       B=0.5693702498468441
       V=0.9105760258970126e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD1454(X,Y,Z,W,N)
       DOUBLE PRECISION X(1454)
       DOUBLE PRECISION Y(1454)
       DOUBLE PRECISION Z(1454)
       DOUBLE PRECISION W(1454)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 1454-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.7777160743261247e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7557646413004701e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3229290663413854e-1
       V=0.2841633806090617e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8036733271462222e-1
       V=0.4374419127053555e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1354289960531653
       V=0.5417174740872172e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1938963861114426
       V=0.6148000891358593e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2537343715011275
       V=0.6664394485800705e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3135251434752570
       V=0.7025039356923220e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3721558339375338
       V=0.7268511789249627e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4286809575195696
       V=0.7422637534208629e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4822510128282994
       V=0.7509545035841214e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5320679333566263
       V=0.7548535057718401e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6172998195394274
       V=0.7554088969774001e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6510679849127481
       V=0.7553147174442808e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6777315251687360
       V=0.7564767653292297e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6963109410648741
       V=0.7587991808518730e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7058935009831749
       V=0.7608261832033027e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9955546194091857
       V=0.4021680447874916e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9734115901794209
       V=0.5804871793945964e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9275693732388626
       V=0.6792151955945159e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8568022422795103
       V=0.7336741211286294e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7623495553719372
       V=0.7581866300989608e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5707522908892223
       B=0.4387028039889501
       V=0.7538257859800743e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5196463388403083
       B=0.3858908414762617
       V=0.7483517247053123e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4646337531215351
       B=0.3301937372343854
       V=0.7371763661112059e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4063901697557691
       B=0.2725423573563777
       V=0.7183448895756934e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3456329466643087
       B=0.2139510237495250
       V=0.6895815529822191e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2831395121050332
       B=0.1555922309786647
       V=0.6480105801792886e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2197682022925330
       B=0.9892878979686097e-1
       V=0.5897558896594636e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1564696098650355
       B=0.4598642910675510e-1
       V=0.5095708849247346e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6027356673721295
       B=0.3376625140173426
       V=0.7536906428909755e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5496032320255096
       B=0.2822301309727988
       V=0.7472505965575118e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4921707755234567
       B=0.2248632342592540
       V=0.7343017132279698e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4309422998598483
       B=0.1666224723456479
       V=0.7130871582177445e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3664108182313672
       B=0.1086964901822169
       V=0.6817022032112776e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2990189057758436
       B=0.5251989784120085e-1
       V=0.6380941145604121e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6268724013144998
       B=0.2297523657550023
       V=0.7550381377920310e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5707324144834607
       B=0.1723080607093800
       V=0.7478646640144802e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5096360901960365
       B=0.1140238465390513
       V=0.7335918720601220e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4438729938312456
       B=0.5611522095882537e-1
       V=0.7110120527658118e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6419978471082389
       B=0.1164174423140873
       V=0.7571363978689501e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5817218061802611
       B=0.5797589531445219e-1
       V=0.7489908329079234e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD1730(X,Y,Z,W,N)
       DOUBLE PRECISION X(1730)
       DOUBLE PRECISION Y(1730)
       DOUBLE PRECISION Z(1730)
       DOUBLE PRECISION W(1730)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 1730-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.6309049437420976e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6398287705571748e-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6357185073530720e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2860923126194662e-1
       V=0.2221207162188168e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7142556767711522e-1
       V=0.3475784022286848e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1209199540995559
       V=0.4350742443589804e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1738673106594379
       V=0.4978569136522127e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2284645438467734
       V=0.5435036221998053e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2834807671701512
       V=0.5765913388219542e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3379680145467339
       V=0.6001200359226003e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3911355454819537
       V=0.6162178172717512e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4422860353001403
       V=0.6265218152438485e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4907781568726057
       V=0.6323987160974212e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5360006153211468
       V=0.6350767851540569e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6142105973596603
       V=0.6354362775297107e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6459300387977504
       V=0.6352302462706235e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718056125089225
       V=0.6358117881417972e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6910888533186254
       V=0.6373101590310117e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7030467416823252
       V=0.6390428961368665e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8354951166354646e-1
       V=0.3186913449946576e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2050143009099486
       V=0.4678028558591711e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3370208290706637
       V=0.5538829697598626e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4689051484233963
       V=0.6044475907190476e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5939400424557334
       V=0.6313575103509012e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1394983311832261
       B=0.4097581162050343e-1
       V=0.4078626431855630e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1967999180485014
       B=0.8851987391293348e-1
       V=0.4759933057812725e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2546183732548967
       B=0.1397680182969819
       V=0.5268151186413440e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3121281074713875
       B=0.1929452542226526
       V=0.5643048560507316e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3685981078502492
       B=0.2467898337061562
       V=0.5914501076613073e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4233760321547856
       B=0.3003104124785409
       V=0.6104561257874195e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4758671236059246
       B=0.3526684328175033
       V=0.6230252860707806e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5255178579796463
       B=0.4031134861145713
       V=0.6305618761760796e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5718025633734589
       B=0.4509426448342351
       V=0.6343092767597889e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2686927772723415
       B=0.4711322502423248e-1
       V=0.5176268945737826e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3306006819904809
       B=0.9784487303942695e-1
       V=0.5564840313313692e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3904906850594983
       B=0.1505395810025273
       V=0.5856426671038980e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4479957951904390
       B=0.2039728156296050
       V=0.6066386925777091e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5027076848919780
       B=0.2571529941121107
       V=0.6208824962234458e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5542087392260217
       B=0.3092191375815670
       V=0.6296314297822907e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6020850887375187
       B=0.3593807506130276
       V=0.6340423756791859e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4019851409179594
       B=0.5063389934378671e-1
       V=0.5829627677107342e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4635614567449800
       B=0.1032422269160612
       V=0.6048693376081110e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5215860931591575
       B=0.1566322094006254
       V=0.6202362317732461e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5758202499099271
       B=0.2098082827491099
       V=0.6299005328403779e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6259893683876795
       B=0.2618824114553391
       V=0.6347722390609353e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5313795124811891
       B=0.5263245019338556e-1
       V=0.6203778981238834e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5893317955931995
       B=0.1061059730982005
       V=0.6308414671239979e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6426246321215801
       B=0.1594171564034221
       V=0.6362706466959498e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6511904367376113
       B=0.5354789536565540e-1
       V=0.6375414170333233e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD2030(X,Y,Z,W,N)
       DOUBLE PRECISION X(2030)
       DOUBLE PRECISION Y(2030)
       DOUBLE PRECISION Z(2030)
       DOUBLE PRECISION W(2030)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 2030-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.4656031899197431e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5421549195295507e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2540835336814348e-1
       V=0.1778522133346553e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6399322800504915e-1
       V=0.2811325405682796e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1088269469804125
       V=0.3548896312631459e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1570670798818287
       V=0.4090310897173364e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2071163932282514
       V=0.4493286134169965e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2578914044450844
       V=0.4793728447962723e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3085687558169623
       V=0.5015415319164265e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3584719706267024
       V=0.5175127372677937e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4070135594428709
       V=0.5285522262081019e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4536618626222638
       V=0.5356832703713962e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4979195686463577
       V=0.5397914736175170e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5393075111126999
       V=0.5416899441599930e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6115617676843916
       V=0.5419308476889938e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6414308435160159
       V=0.5416936902030596e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6664099412721607
       V=0.5419544338703164e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6859161771214913
       V=0.5428983656630975e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6993625593503890
       V=0.5442286500098193e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7062393387719380
       V=0.5452250345057301e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7479028168349763e-1
       V=0.2568002497728530e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1848951153969366
       V=0.3827211700292145e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3059529066581305
       V=0.4579491561917824e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4285556101021362
       V=0.5042003969083574e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5468758653496526
       V=0.5312708889976025e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6565821978343439
       V=0.5438401790747117e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1253901572367117
       B=0.3681917226439641e-1
       V=0.3316041873197344e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1775721510383941
       B=0.7982487607213301e-1
       V=0.3899113567153771e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2305693358216114
       B=0.1264640966592335
       V=0.4343343327201309e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2836502845992063
       B=0.1751585683418957
       V=0.4679415262318919e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3361794746232590
       B=0.2247995907632670
       V=0.4930847981631031e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3875979172264824
       B=0.2745299257422246
       V=0.5115031867540091e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4374019316999074
       B=0.3236373482441118
       V=0.5245217148457367e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4851275843340022
       B=0.3714967859436741
       V=0.5332041499895321e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5303391803806868
       B=0.4175353646321745
       V=0.5384583126021542e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726197380596287
       B=0.4612084406355461
       V=0.5411067210798852e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2431520732564863
       B=0.4258040133043952e-1
       V=0.4259797391468714e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3002096800895869
       B=0.8869424306722721e-1
       V=0.4604931368460021e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3558554457457432
       B=0.1368811706510655
       V=0.4871814878255202e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4097782537048887
       B=0.1860739985015033
       V=0.5072242910074885e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4616337666067458
       B=0.2354235077395853
       V=0.5217069845235350e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5110707008417874
       B=0.2842074921347011
       V=0.5315785966280310e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5577415286163795
       B=0.3317784414984102
       V=0.5376833708758905e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6013060431366950
       B=0.3775299002040700
       V=0.5408032092069521e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3661596767261781
       B=0.4599367887164592e-1
       V=0.4842744917904866e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4237633153506581
       B=0.9404893773654421e-1
       V=0.5048926076188130e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4786328454658452
       B=0.1431377109091971
       V=0.5202607980478373e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5305702076789774
       B=0.1924186388843570
       V=0.5309932388325743e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5793436224231788
       B=0.2411590944775190
       V=0.5377419770895208e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6247069017094747
       B=0.2886871491583605
       V=0.5411696331677717e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4874315552535204
       B=0.4804978774953206e-1
       V=0.5197996293282420e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5427337322059053
       B=0.9716857199366665e-1
       V=0.5311120836622945e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5943493747246700
       B=0.1465205839795055
       V=0.5384309319956951e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6421314033564943
       B=0.1953579449803574
       V=0.5421859504051886e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6020628374713980
       B=0.4916375015738108e-1
       V=0.5390948355046314e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6529222529856881
       B=0.9861621540127005e-1
       V=0.5433312705027845e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD2354(X,Y,Z,W,N)
       DOUBLE PRECISION X(2354)
       DOUBLE PRECISION Y(2354)
       DOUBLE PRECISION Z(2354)
       DOUBLE PRECISION W(2354)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 2354-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.3922616270665292e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4703831750854424e-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4678202801282136e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2290024646530589e-1
       V=0.1437832228979900e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5779086652271284e-1
       V=0.2303572493577644e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9863103576375984e-1
       V=0.2933110752447454e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1428155792982185
       V=0.3402905998359838e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1888978116601463
       V=0.3759138466870372e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2359091682970210
       V=0.4030638447899798e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2831228833706171
       V=0.4236591432242211e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3299495857966693
       V=0.4390522656946746e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3758840802660796
       V=0.4502523466626247e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4204751831009480
       V=0.4580577727783541e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4633068518751051
       V=0.4631391616615899e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5039849474507313
       V=0.4660928953698676e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5421265793440747
       V=0.4674751807936953e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6092660230557310
       V=0.4676414903932920e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6374654204984869
       V=0.4674086492347870e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6615136472609892
       V=0.4674928539483207e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6809487285958127
       V=0.4680748979686447e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6952980021665196
       V=0.4690449806389040e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7041245497695400
       V=0.4699877075860818e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6744033088306065e-1
       V=0.2099942281069176e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1678684485334166
       V=0.3172269150712804e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2793559049539613
       V=0.3832051358546523e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3935264218057639
       V=0.4252193818146985e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5052629268232558
       V=0.4513807963755000e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6107905315437531
       V=0.4657797469114178e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1135081039843524
       B=0.3331954884662588e-1
       V=0.2733362800522836e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1612866626099378
       B=0.7247167465436538e-1
       V=0.3235485368463559e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2100786550168205
       B=0.1151539110849745
       V=0.3624908726013453e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2592282009459942
       B=0.1599491097143677
       V=0.3925540070712828e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3081740561320203
       B=0.2058699956028027
       V=0.4156129781116235e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3564289781578164
       B=0.2521624953502911
       V=0.4330644984623263e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4035587288240703
       B=0.2982090785797674
       V=0.4459677725921312e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4491671196373903
       B=0.3434762087235733
       V=0.4551593004456795e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4928854782917489
       B=0.3874831357203437
       V=0.4613341462749918e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5343646791958988
       B=0.4297814821746926
       V=0.4651019618269806e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5732683216530990
       B=0.4699402260943537
       V=0.4670249536100625e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2214131583218986
       B=0.3873602040643895e-1
       V=0.3549555576441708e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2741796504750071
       B=0.8089496256902013e-1
       V=0.3856108245249010e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3259797439149485
       B=0.1251732177620872
       V=0.4098622845756882e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3765441148826891
       B=0.1706260286403185
       V=0.4286328604268950e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4255773574530558
       B=0.2165115147300408
       V=0.4427802198993945e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4727795117058430
       B=0.2622089812225259
       V=0.4530473511488561e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5178546895819012
       B=0.3071721431296201
       V=0.4600805475703138e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5605141192097460
       B=0.3508998998801138
       V=0.4644599059958017e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6004763319352512
       B=0.3929160876166931
       V=0.4667274455712508e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3352842634946949
       B=0.4202563457288019e-1
       V=0.4069360518020356e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3891971629814670
       B=0.8614309758870850e-1
       V=0.4260442819919195e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4409875565542281
       B=0.1314500879380001
       V=0.4408678508029063e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4904893058592484
       B=0.1772189657383859
       V=0.4518748115548597e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5375056138769549
       B=0.2228277110050294
       V=0.4595564875375116e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5818255708669969
       B=0.2677179935014386
       V=0.4643988774315846e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6232334858144959
       B=0.3113675035544165
       V=0.4668827491646946e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4489485354492058
       B=0.4409162378368174e-1
       V=0.4400541823741973e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5015136875933150
       B=0.8939009917748489e-1
       V=0.4514512890193797e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5511300550512623
       B=0.1351806029383365
       V=0.4596198627347549e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5976720409858000
       B=0.1808370355053196
       V=0.4648659016801781e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6409956378989354
       B=0.2257852192301602
       V=0.4675502017157673e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5581222330827514
       B=0.4532173421637160e-1
       V=0.4598494476455523e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6074705984161695
       B=0.9117488031840314e-1
       V=0.4654916955152048e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6532272537379033
       B=0.1369294213140155
       V=0.4684709779505137e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6594761494500487
       B=0.4589901487275583e-1
       V=0.4691445539106986e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD2702(X,Y,Z,W,N)
       DOUBLE PRECISION X(2702)
       DOUBLE PRECISION Y(2702)
       DOUBLE PRECISION Z(2702)
       DOUBLE PRECISION W(2702)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 2702-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.2998675149888161e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4077860529495355e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2065562538818703e-1
       V=0.1185349192520667e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5250918173022379e-1
       V=0.1913408643425751e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8993480082038376e-1
       V=0.2452886577209897e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1306023924436019
       V=0.2862408183288702e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1732060388531418
       V=0.3178032258257357e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2168727084820249
       V=0.3422945667633690e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2609528309173586
       V=0.3612790520235922e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3049252927938952
       V=0.3758638229818521e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3483484138084404
       V=0.3868711798859953e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3908321549106406
       V=0.3949429933189938e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4320210071894814
       V=0.4006068107541156e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4715824795890053
       V=0.4043192149672723e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5091984794078453
       V=0.4064947495808078e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5445580145650803
       V=0.4075245619813152e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6072575796841768
       V=0.4076423540893566e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6339484505755803
       V=0.4074280862251555e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6570718257486958
       V=0.4074163756012244e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6762557330090709
       V=0.4077647795071246e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6911161696923790
       V=0.4084517552782530e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012841911659961
       V=0.4092468459224052e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7064559272410020
       V=0.4097872687240906e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6123554989894765e-1
       V=0.1738986811745028e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1533070348312393
       V=0.2659616045280191e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2563902605244206
       V=0.3240596008171533e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3629346991663361
       V=0.3621195964432943e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4683949968987538
       V=0.3868838330760539e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5694479240657952
       V=0.4018911532693111e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6634465430993955
       V=0.4089929432983252e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1033958573552305
       B=0.3034544009063584e-1
       V=0.2279907527706409e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1473521412414395
       B=0.6618803044247135e-1
       V=0.2715205490578897e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1924552158705967
       B=0.1054431128987715
       V=0.3057917896703976e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2381094362890328
       B=0.1468263551238858
       V=0.3326913052452555e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2838121707936760
       B=0.1894486108187886
       V=0.3537334711890037e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3291323133373415
       B=0.2326374238761579
       V=0.3700567500783129e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3736896978741460
       B=0.2758485808485768
       V=0.3825245372589122e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4171406040760013
       B=0.3186179331996921
       V=0.3918125171518296e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4591677985256915
       B=0.3605329796303794
       V=0.3984720419937579e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4994733831718418
       B=0.4012147253586509
       V=0.4029746003338211e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5377731830445096
       B=0.4403050025570692
       V=0.4057428632156627e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5737917830001331
       B=0.4774565904277483
       V=0.4071719274114857e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2027323586271389
       B=0.3544122504976147e-1
       V=0.2990236950664119e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2516942375187273
       B=0.7418304388646328e-1
       V=0.3262951734212878e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3000227995257181
       B=0.1150502745727186
       V=0.3482634608242413e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3474806691046342
       B=0.1571963371209364
       V=0.3656596681700892e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3938103180359209
       B=0.1999631877247100
       V=0.3791740467794218e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4387519590455703
       B=0.2428073457846535
       V=0.3894034450156905e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4820503960077787
       B=0.2852575132906155
       V=0.3968600245508371e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5234573778475101
       B=0.3268884208674639
       V=0.4019931351420050e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5627318647235282
       B=0.3673033321675939
       V=0.4052108801278599e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5996390607156954
       B=0.4061211551830290
       V=0.4068978613940934e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3084780753791947
       B=0.3860125523100059e-1
       V=0.3454275351319704e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3589988275920223
       B=0.7928938987104867e-1
       V=0.3629963537007920e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4078628415881973
       B=0.1212614643030087
       V=0.3770187233889873e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4549287258889735
       B=0.1638770827382693
       V=0.3878608613694378e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5000278512957279
       B=0.2065965798260176
       V=0.3959065270221274e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5429785044928199
       B=0.2489436378852235
       V=0.4015286975463570e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5835939850491711
       B=0.2904811368946891
       V=0.4050866785614717e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6216870353444856
       B=0.3307941957666609
       V=0.4069320185051913e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4151104662709091
       B=0.4064829146052554e-1
       V=0.3760120964062763e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4649804275009218
       B=0.8258424547294755e-1
       V=0.3870969564418064e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5124695757009662
       B=0.1251841962027289
       V=0.3955287790534055e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5574711100606224
       B=0.1679107505976331
       V=0.4015361911302668e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5998597333287227
       B=0.2102805057358715
       V=0.4053836986719548e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6395007148516600
       B=0.2518418087774107
       V=0.4073578673299117e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5188456224746252
       B=0.4194321676077518e-1
       V=0.3954628379231406e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5664190707942778
       B=0.8457661551921499e-1
       V=0.4017645508847530e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6110464353283153
       B=0.1273652932519396
       V=0.4059030348651293e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6526430302051563
       B=0.1698173239076354
       V=0.4080565809484880e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6167551880377548
       B=0.4266398851548864e-1
       V=0.4063018753664651e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6607195418355383
       B=0.8551925814238349e-1
       V=0.4087191292799671e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD3074(X,Y,Z,W,N)
       DOUBLE PRECISION X(3074)
       DOUBLE PRECISION Y(3074)
       DOUBLE PRECISION Z(3074)
       DOUBLE PRECISION W(3074)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 3074-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.2599095953754734e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3603134089687541e-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3586067974412447e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1886108518723392e-1
       V=0.9831528474385880e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4800217244625303e-1
       V=0.1605023107954450e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8244922058397242e-1
       V=0.2072200131464099e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1200408362484023
       V=0.2431297618814187e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1595773530809965
       V=0.2711819064496707e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2002635973434064
       V=0.2932762038321116e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2415127590139982
       V=0.3107032514197368e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2828584158458477
       V=0.3243808058921213e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3239091015338138
       V=0.3349899091374030e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3643225097962194
       V=0.3430580688505218e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4037897083691802
       V=0.3490124109290343e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4420247515194127
       V=0.3532148948561955e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4787572538464938
       V=0.3559862669062833e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5137265251275234
       V=0.3576224317551411e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5466764056654611
       V=0.3584050533086076e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6054859420813535
       V=0.3584903581373224e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6308106701764562
       V=0.3582991879040586e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6530369230179584
       V=0.3582371187963125e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718609524611158
       V=0.3584353631122350e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6869676499894013
       V=0.3589120166517785e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6980467077240748
       V=0.3595445704531601e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7048241721250522
       V=0.3600943557111074e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5591105222058232e-1
       V=0.1456447096742039e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1407384078513916
       V=0.2252370188283782e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2364035438976309
       V=0.2766135443474897e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3360602737818170
       V=0.3110729491500851e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4356292630054665
       V=0.3342506712303391e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5321569415256174
       V=0.3491981834026860e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6232956305040554
       V=0.3576003604348932e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9469870086838469e-1
       B=0.2778748387309470e-1
       V=0.1921921305788564e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1353170300568141
       B=0.6076569878628364e-1
       V=0.2301458216495632e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1771679481726077
       B=0.9703072762711040e-1
       V=0.2604248549522893e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2197066664231751
       B=0.1354112458524762
       V=0.2845275425870697e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2624783557374927
       B=0.1750996479744100
       V=0.3036870897974840e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3050969521214442
       B=0.2154896907449802
       V=0.3188414832298066e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3472252637196021
       B=0.2560954625740152
       V=0.3307046414722089e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3885610219026360
       B=0.2965070050624096
       V=0.3398330969031360e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4288273776062765
       B=0.3363641488734497
       V=0.3466757899705373e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4677662471302948
       B=0.3753400029836788
       V=0.3516095923230054e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5051333589553359
       B=0.4131297522144286
       V=0.3549645184048486e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5406942145810492
       B=0.4494423776081795
       V=0.3570415969441392e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5742204122576457
       B=0.4839938958841502
       V=0.3581251798496118e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1865407027225188
       B=0.3259144851070796e-1
       V=0.2543491329913348e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2321186453689432
       B=0.6835679505297343e-1
       V=0.2786711051330776e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2773159142523882
       B=0.1062284864451989
       V=0.2985552361083679e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3219200192237254
       B=0.1454404409323047
       V=0.3145867929154039e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3657032593944029
       B=0.1854018282582510
       V=0.3273290662067609e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4084376778363622
       B=0.2256297412014750
       V=0.3372705511943501e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4499004945751427
       B=0.2657104425000896
       V=0.3448274437851510e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4898758141326335
       B=0.3052755487631557
       V=0.3503592783048583e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5281547442266309
       B=0.3439863920645423
       V=0.3541854792663162e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5645346989813992
       B=0.3815229456121914
       V=0.3565995517909428e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5988181252159848
       B=0.4175752420966734
       V=0.3578802078302898e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2850425424471603
       B=0.3562149509862536e-1
       V=0.2958644592860982e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3324619433027876
       B=0.7330318886871096e-1
       V=0.3119548129116835e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3785848333076282
       B=0.1123226296008472
       V=0.3250745225005984e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4232891028562115
       B=0.1521084193337708
       V=0.3355153415935208e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4664287050829722
       B=0.1921844459223610
       V=0.3435847568549328e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5078458493735726
       B=0.2321360989678303
       V=0.3495786831622488e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5473779816204180
       B=0.2715886486360520
       V=0.3537767805534621e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5848617133811376
       B=0.3101924707571355
       V=0.3564459815421428e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6201348281584888
       B=0.3476121052890973
       V=0.3578464061225468e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3852191185387871
       B=0.3763224880035108e-1
       V=0.3239748762836212e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4325025061073423
       B=0.7659581935637135e-1
       V=0.3345491784174287e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4778486229734490
       B=0.1163381306083900
       V=0.3429126177301782e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5211663693009000
       B=0.1563890598752899
       V=0.3492420343097421e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5623469504853703
       B=0.1963320810149200
       V=0.3537399050235257e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6012718188659246
       B=0.2357847407258738
       V=0.3566209152659172e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6378179206390117
       B=0.2743846121244060
       V=0.3581084321919782e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4836936460214534
       B=0.3895902610739024e-1
       V=0.3426522117591512e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5293792562683797
       B=0.7871246819312640e-1
       V=0.3491848770121379e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726281253100033
       B=0.1187963808202981
       V=0.3539318235231476e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6133658776169068
       B=0.1587914708061787
       V=0.3570231438458694e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6515085491865307
       B=0.1983058575227646
       V=0.3586207335051714e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5778692716064976
       B=0.3977209689791542e-1
       V=0.3541196205164025e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6207904288086192
       B=0.7990157592981152e-1
       V=0.3574296911573953e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6608688171046802
       B=0.1199671308754309
       V=0.3591993279818963e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6656263089489130
       B=0.4015955957805969e-1
       V=0.3595855034661997e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD3470(X,Y,Z,W,N)
       DOUBLE PRECISION X(3470)
       DOUBLE PRECISION Y(3470)
       DOUBLE PRECISION Z(3470)
       DOUBLE PRECISION W(3470)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 3470-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.2040382730826330e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3178149703889544e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1721420832906233e-1
       V=0.8288115128076110e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4408875374981770e-1
       V=0.1360883192522954e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7594680813878681e-1
       V=0.1766854454542662e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1108335359204799
       V=0.2083153161230153e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1476517054388567
       V=0.2333279544657158e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1856731870860615
       V=0.2532809539930247e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2243634099428821
       V=0.2692472184211158e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2633006881662727
       V=0.2819949946811885e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3021340904916283
       V=0.2920953593973030e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3405594048030089
       V=0.2999889782948352e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3783044434007372
       V=0.3060292120496902e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4151194767407910
       V=0.3105109167522192e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4507705766443257
       V=0.3136902387550312e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4850346056573187
       V=0.3157984652454632e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5176950817792470
       V=0.3170516518425422e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5485384240820989
       V=0.3176568425633755e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6039117238943308
       V=0.3177198411207062e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6279956655573113
       V=0.3175519492394733e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6493636169568952
       V=0.3174654952634756e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6677644117704504
       V=0.3175676415467654e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6829368572115624
       V=0.3178923417835410e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6946195818184121
       V=0.3183788287531909e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7025711542057026
       V=0.3188755151918807e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7066004767140119
       V=0.3191916889313849e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5132537689946062e-1
       V=0.1231779611744508e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1297994661331225
       V=0.1924661373839880e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2188852049401307
       V=0.2380881867403424e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3123174824903457
       V=0.2693100663037885e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4064037620738195
       V=0.2908673382834366e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4984958396944782
       V=0.3053914619381535e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5864975046021365
       V=0.3143916684147777e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6686711634580175
       V=0.3187042244055363e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8715738780835950e-1
       B=0.2557175233367578e-1
       V=0.1635219535869790e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1248383123134007
       B=0.5604823383376681e-1
       V=0.1968109917696070e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1638062693383378
       B=0.8968568601900765e-1
       V=0.2236754342249974e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2035586203373176
       B=0.1254086651976279
       V=0.2453186687017181e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2436798975293774
       B=0.1624780150162012
       V=0.2627551791580541e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2838207507773806
       B=0.2003422342683208
       V=0.2767654860152220e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3236787502217692
       B=0.2385628026255263
       V=0.2879467027765895e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3629849554840691
       B=0.2767731148783578
       V=0.2967639918918702e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4014948081992087
       B=0.3146542308245309
       V=0.3035900684660351e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4389818379260225
       B=0.3519196415895088
       V=0.3087338237298308e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4752331143674377
       B=0.3883050984023654
       V=0.3124608838860167e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5100457318374018
       B=0.4235613423908649
       V=0.3150084294226743e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5432238388954868
       B=0.4574484717196220
       V=0.3165958398598402e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5745758685072442
       B=0.4897311639255524
       V=0.3174320440957372e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1723981437592809
       B=0.3010630597881105e-1
       V=0.2182188909812599e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2149553257844597
       B=0.6326031554204694e-1
       V=0.2399727933921445e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2573256081247422
       B=0.9848566980258631e-1
       V=0.2579796133514652e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2993163751238106
       B=0.1350835952384266
       V=0.2727114052623535e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3407238005148000
       B=0.1725184055442181
       V=0.2846327656281355e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3813454978483264
       B=0.2103559279730725
       V=0.2941491102051334e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4209848104423343
       B=0.2482278774554860
       V=0.3016049492136107e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4594519699996300
       B=0.2858099509982883
       V=0.3072949726175648e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4965640166185930
       B=0.3228075659915428
       V=0.3114768142886460e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5321441655571562
       B=0.3589459907204151
       V=0.3143823673666223e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5660208438582166
       B=0.3939630088864310
       V=0.3162269764661535e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5980264315964364
       B=0.4276029922949089
       V=0.3172164663759821e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2644215852350733
       B=0.3300939429072552e-1
       V=0.2554575398967435e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3090113743443063
       B=0.6803887650078501e-1
       V=0.2701704069135677e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3525871079197808
       B=0.1044326136206709
       V=0.2823693413468940e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3950418005354029
       B=0.1416751597517679
       V=0.2922898463214289e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4362475663430163
       B=0.1793408610504821
       V=0.3001829062162428e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4760661812145854
       B=0.2170630750175722
       V=0.3062890864542953e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5143551042512103
       B=0.2545145157815807
       V=0.3108328279264746e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5509709026935597
       B=0.2913940101706601
       V=0.3140243146201245e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5857711030329428
       B=0.3274169910910705
       V=0.3160638030977130e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6186149917404392
       B=0.3623081329317265
       V=0.3171462882206275e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3586894569557064
       B=0.3497354386450040e-1
       V=0.2812388416031796e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4035266610019441
       B=0.7129736739757095e-1
       V=0.2912137500288045e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4467775312332510
       B=0.1084758620193165
       V=0.2993241256502206e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4883638346608543
       B=0.1460915689241772
       V=0.3057101738983822e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5281908348434601
       B=0.1837790832369980
       V=0.3105319326251432e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5661542687149311
       B=0.2212075390874021
       V=0.3139565514428167e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6021450102031452
       B=0.2580682841160985
       V=0.3161543006806366e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6360520783610050
       B=0.2940656362094121
       V=0.3172985960613294e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4521611065087196
       B=0.3631055365867002e-1
       V=0.2989400336901431e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4959365651560963
       B=0.7348318468484350e-1
       V=0.3054555883947677e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5376815804038283
       B=0.1111087643812648
       V=0.3104764960807702e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5773314480243768
       B=0.1488226085145408
       V=0.3141015825977616e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6148113245575056
       B=0.1862892274135151
       V=0.3164520621159896e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6500407462842380
       B=0.2231909701714456
       V=0.3176652305912204e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5425151448707213
       B=0.3718201306118944e-1
       V=0.3105097161023939e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5841860556907931
       B=0.7483616335067346e-1
       V=0.3143014117890550e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6234632186851500
       B=0.1125990834266120
       V=0.3168172866287200e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6602934551848843
       B=0.1501303813157619
       V=0.3181401865570968e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6278573968375105
       B=0.3767559930245720e-1
       V=0.3170663659156037e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6665611711264577
       B=0.7548443301360158e-1
       V=0.3185447944625510e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD3890(X,Y,Z,W,N)
       DOUBLE PRECISION X(3890)
       DOUBLE PRECISION Y(3890)
       DOUBLE PRECISION Z(3890)
       DOUBLE PRECISION W(3890)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 3890-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.1807395252196920e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2848008782238827e-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2836065837530581e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1587876419858352e-1
       V=0.7013149266673816e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4069193593751206e-1
       V=0.1162798021956766e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7025888115257997e-1
       V=0.1518728583972105e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1027495450028704
       V=0.1798796108216934e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1371457730893426
       V=0.2022593385972785e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1727758532671953
       V=0.2203093105575464e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2091492038929037
       V=0.2349294234299855e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2458813281751915
       V=0.2467682058747003e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2826545859450066
       V=0.2563092683572224e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3191957291799622
       V=0.2639253896763318e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3552621469299578
       V=0.2699137479265108e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3906329503406230
       V=0.2745196420166739e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4251028614093031
       V=0.2779529197397593e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4584777520111870
       V=0.2803996086684265e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4905711358710193
       V=0.2820302356715842e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5212011669847385
       V=0.2830056747491068e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5501878488737995
       V=0.2834808950776839e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6025037877479342
       V=0.2835282339078929e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6254572689549016
       V=0.2833819267065800e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6460107179528248
       V=0.2832858336906784e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6639541138154251
       V=0.2833268235451244e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6790688515667495
       V=0.2835432677029253e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6911338580371512
       V=0.2839091722743049e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6999385956126490
       V=0.2843308178875841e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7053037748656896
       V=0.2846703550533846e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4732224387180115e-1
       V=0.1051193406971900e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1202100529326803
       V=0.1657871838796974e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2034304820664855
       V=0.2064648113714232e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2912285643573002
       V=0.2347942745819741e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3802361792726768
       V=0.2547775326597726e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4680598511056146
       V=0.2686876684847025e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5528151052155599
       V=0.2778665755515867e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6329386307803041
       V=0.2830996616782929e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8056516651369069e-1
       B=0.2363454684003124e-1
       V=0.1403063340168372e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1156476077139389
       B=0.5191291632545936e-1
       V=0.1696504125939477e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1520473382760421
       B=0.8322715736994519e-1
       V=0.1935787242745390e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1892986699745931
       B=0.1165855667993712
       V=0.2130614510521968e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2270194446777792
       B=0.1513077167409504
       V=0.2289381265931048e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2648908185093273
       B=0.1868882025807859
       V=0.2418630292816186e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3026389259574136
       B=0.2229277629776224
       V=0.2523400495631193e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3400220296151384
       B=0.2590951840746235
       V=0.2607623973449605e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3768217953335510
       B=0.2951047291750847
       V=0.2674441032689209e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4128372900921884
       B=0.3307019714169930
       V=0.2726432360343356e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4478807131815630
       B=0.3656544101087634
       V=0.2765787685924545e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4817742034089257
       B=0.3997448951939695
       V=0.2794428690642224e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5143472814653344
       B=0.4327667110812024
       V=0.2814099002062895e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5454346213905650
       B=0.4645196123532293
       V=0.2826429531578994e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5748739313170252
       B=0.4948063555703345
       V=0.2832983542550884e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1599598738286342
       B=0.2792357590048985e-1
       V=0.1886695565284976e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1998097412500951
       B=0.5877141038139065e-1
       V=0.2081867882748234e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2396228952566202
       B=0.9164573914691377e-1
       V=0.2245148680600796e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2792228341097746
       B=0.1259049641962687
       V=0.2380370491511872e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3184251107546741
       B=0.1610594823400863
       V=0.2491398041852455e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3570481164426244
       B=0.1967151653460898
       V=0.2581632405881230e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3949164710492144
       B=0.2325404606175168
       V=0.2653965506227417e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318617293970503
       B=0.2682461141151439
       V=0.2710857216747087e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4677221009931678
       B=0.3035720116011973
       V=0.2754434093903659e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5023417939270955
       B=0.3382781859197439
       V=0.2786579932519380e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5355701836636128
       B=0.3721383065625942
       V=0.2809011080679474e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5672608451328771
       B=0.4049346360466055
       V=0.2823336184560987e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5972704202540162
       B=0.4364538098633802
       V=0.2831101175806309e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2461687022333596
       B=0.3070423166833368e-1
       V=0.2221679970354546e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2881774566286831
       B=0.6338034669281885e-1
       V=0.2356185734270703e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3293963604116978
       B=0.9742862487067941e-1
       V=0.2469228344805590e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3697303822241377
       B=0.1323799532282290
       V=0.2562726348642046e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4090663023135127
       B=0.1678497018129336
       V=0.2638756726753028e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4472819355411712
       B=0.2035095105326114
       V=0.2699311157390862e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4842513377231437
       B=0.2390692566672091
       V=0.2746233268403837e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5198477629962928
       B=0.2742649818076149
       V=0.2781225674454771e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5539453011883145
       B=0.3088503806580094
       V=0.2805881254045684e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5864196762401251
       B=0.3425904245906614
       V=0.2821719877004913e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6171484466668390
       B=0.3752562294789468
       V=0.2830222502333124e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3350337830565727
       B=0.3261589934634747e-1
       V=0.2457995956744870e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3775773224758284
       B=0.6658438928081572e-1
       V=0.2551474407503706e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4188155229848973
       B=0.1014565797157954
       V=0.2629065335195311e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4586805892009344
       B=0.1368573320843822
       V=0.2691900449925075e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4970895714224235
       B=0.1724614851951608
       V=0.2741275485754276e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5339505133960747
       B=0.2079779381416412
       V=0.2778530970122595e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5691665792531440
       B=0.2431385788322288
       V=0.2805010567646741e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6026387682680377
       B=0.2776901883049853
       V=0.2822055834031040e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6342676150163307
       B=0.3113881356386632
       V=0.2831016901243473e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4237951119537067
       B=0.3394877848664351e-1
       V=0.2624474901131803e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4656918683234929
       B=0.6880219556291447e-1
       V=0.2688034163039377e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5058857069185980
       B=0.1041946859721635
       V=0.2738932751287636e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5443204666713996
       B=0.1398039738736393
       V=0.2777944791242523e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5809298813759742
       B=0.1753373381196155
       V=0.2806011661660987e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6156416039447128
       B=0.2105215793514010
       V=0.2824181456597460e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6483801351066604
       B=0.2450953312157051
       V=0.2833585216577828e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5103616577251688
       B=0.3485560643800719e-1
       V=0.2738165236962878e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5506738792580681
       B=0.7026308631512033e-1
       V=0.2778365208203180e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5889573040995292
       B=0.1059035061296403
       V=0.2807852940418966e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6251641589516930
       B=0.1414823925236026
       V=0.2827245949674705e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6592414921570178
       B=0.1767207908214530
       V=0.2837342344829828e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5930314017533384
       B=0.3542189339561672e-1
       V=0.2809233907610981e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6309812253390175
       B=0.7109574040369549e-1
       V=0.2829930809742694e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6666296011353230
       B=0.1067259792282730
       V=0.2841097874111479e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6703715271049922
       B=0.3569455268820809e-1
       V=0.2843455206008783e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD4334(X,Y,Z,W,N)
       DOUBLE PRECISION X(4334)
       DOUBLE PRECISION Y(4334)
       DOUBLE PRECISION Z(4334)
       DOUBLE PRECISION W(4334)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 4334-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.1449063022537883e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2546377329828424e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1462896151831013e-1
       V=0.6018432961087496e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3769840812493139e-1
       V=0.1002286583263673e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6524701904096891e-1
       V=0.1315222931028093e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9560543416134648e-1
       V=0.1564213746876724e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1278335898929198
       V=0.1765118841507736e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1613096104466031
       V=0.1928737099311080e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1955806225745371
       V=0.2062658534263270e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2302935218498028
       V=0.2172395445953787e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2651584344113027
       V=0.2262076188876047e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2999276825183209
       V=0.2334885699462397e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3343828669718798
       V=0.2393355273179203e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3683265013750518
       V=0.2439559200468863e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4015763206518108
       V=0.2475251866060002e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4339612026399770
       V=0.2501965558158773e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4653180651114582
       V=0.2521081407925925e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4954893331080803
       V=0.2533881002388081e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5243207068924930
       V=0.2541582900848261e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5516590479041704
       V=0.2545365737525860e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6012371927804176
       V=0.2545726993066799e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6231574466449819
       V=0.2544456197465555e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6429416514181271
       V=0.2543481596881064e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6604124272943595
       V=0.2543506451429194e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6753851470408250
       V=0.2544905675493763e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6876717970626160
       V=0.2547611407344429e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6970895061319234
       V=0.2551060375448869e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7034746912553310
       V=0.2554291933816039e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7067017217542295
       V=0.2556255710686343e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4382223501131123e-1
       V=0.9041339695118195e-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1117474077400006
       V=0.1438426330079022e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1897153252911440
       V=0.1802523089820518e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2724023009910331
       V=0.2060052290565496e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3567163308709902
       V=0.2245002248967466e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4404784483028087
       V=0.2377059847731150e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5219833154161411
       V=0.2468118955882525e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5998179868977553
       V=0.2525410872966528e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6727803154548222
       V=0.2553101409933397e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7476563943166086e-1
       B=0.2193168509461185e-1
       V=0.1212879733668632e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1075341482001416
       B=0.4826419281533887e-1
       V=0.1472872881270931e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1416344885203259
       B=0.7751191883575742e-1
       V=0.1686846601010828e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1766325315388586
       B=0.1087558139247680
       V=0.1862698414660208e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2121744174481514
       B=0.1413661374253096
       V=0.2007430956991861e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2479669443408145
       B=0.1748768214258880
       V=0.2126568125394796e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2837600452294113
       B=0.2089216406612073
       V=0.2224394603372113e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3193344933193984
       B=0.2431987685545972
       V=0.2304264522673135e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3544935442438745
       B=0.2774497054377770
       V=0.2368854288424087e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3890571932288154
       B=0.3114460356156915
       V=0.2420352089461772e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4228581214259090
       B=0.3449806851913012
       V=0.2460597113081295e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4557387211304052
       B=0.3778618641248256
       V=0.2491181912257687e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4875487950541643
       B=0.4099086391698978
       V=0.2513528194205857e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5181436529962997
       B=0.4409474925853973
       V=0.2528943096693220e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5473824095600661
       B=0.4708094517711291
       V=0.2538660368488136e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5751263398976174
       B=0.4993275140354637
       V=0.2543868648299022e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1489515746840028
       B=0.2599381993267017e-1
       V=0.1642595537825183e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1863656444351767
       B=0.5479286532462190e-1
       V=0.1818246659849308e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2238602880356348
       B=0.8556763251425254e-1
       V=0.1966565649492420e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2612723375728160
       B=0.1177257802267011
       V=0.2090677905657991e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2984332990206190
       B=0.1508168456192700
       V=0.2193820409510504e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3351786584663333
       B=0.1844801892177727
       V=0.2278870827661928e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3713505522209120
       B=0.2184145236087598
       V=0.2348283192282090e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4067981098954663
       B=0.2523590641486229
       V=0.2404139755581477e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4413769993687534
       B=0.2860812976901373
       V=0.2448227407760734e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4749487182516394
       B=0.3193686757808996
       V=0.2482110455592573e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5073798105075426
       B=0.3520226949547602
       V=0.2507192397774103e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5385410448878654
       B=0.3838544395667890
       V=0.2524765968534880e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5683065353670530
       B=0.4146810037640963
       V=0.2536052388539425e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5965527620663510
       B=0.4443224094681121
       V=0.2542230588033068e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2299227700856157
       B=0.2865757664057584e-1
       V=0.1944817013047896e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2695752998553267
       B=0.5923421684485993e-1
       V=0.2067862362746635e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3086178716611389
       B=0.9117817776057715e-1
       V=0.2172440734649114e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3469649871659077
       B=0.1240593814082605
       V=0.2260125991723423e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3845153566319655
       B=0.1575272058259175
       V=0.2332655008689523e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4211600033403215
       B=0.1912845163525413
       V=0.2391699681532458e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4567867834329882
       B=0.2250710177858171
       V=0.2438801528273928e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4912829319232061
       B=0.2586521303440910
       V=0.2475370504260665e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5245364793303812
       B=0.2918112242865407
       V=0.2502707235640574e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5564369788915756
       B=0.3243439239067890
       V=0.2522031701054241e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5868757697775287
       B=0.3560536787835351
       V=0.2534511269978784e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6157458853519617
       B=0.3867480821242581
       V=0.2541284914955151e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3138461110672113
       B=0.3051374637507278e-1
       V=0.2161509250688394e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3542495872050569
       B=0.6237111233730755e-1
       V=0.2248778513437852e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3935751553120181
       B=0.9516223952401907e-1
       V=0.2322388803404617e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4317634668111147
       B=0.1285467341508517
       V=0.2383265471001355e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4687413842250821
       B=0.1622318931656033
       V=0.2432476675019525e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5044274237060283
       B=0.1959581153836453
       V=0.2471122223750674e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5387354077925727
       B=0.2294888081183837
       V=0.2500291752486870e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5715768898356105
       B=0.2626031152713945
       V=0.2521055942764682e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6028627200136111
       B=0.2950904075286713
       V=0.2534472785575503e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6325039812653463
       B=0.3267458451113286
       V=0.2541599713080121e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3981986708423407
       B=0.3183291458749821e-1
       V=0.2317380975862936e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4382791182133300
       B=0.6459548193880908e-1
       V=0.2378550733719775e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4769233057218166
       B=0.9795757037087952e-1
       V=0.2428884456739118e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5140823911194238
       B=0.1316307235126655
       V=0.2469002655757292e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5496977833862983
       B=0.1653556486358704
       V=0.2499657574265851e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5837047306512727
       B=0.1988931724126510
       V=0.2521676168486082e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6160349566926879
       B=0.2320174581438950
       V=0.2535935662645334e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6466185353209440
       B=0.2645106562168662
       V=0.2543356743363214e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4810835158795404
       B=0.3275917807743992e-1
       V=0.2427353285201535e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5199925041324341
       B=0.6612546183967181e-1
       V=0.2468258039744386e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5571717692207494
       B=0.9981498331474143e-1
       V=0.2500060956440310e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5925789250836378
       B=0.1335687001410374
       V=0.2523238365420979e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6261658523859670
       B=0.1671444402896463
       V=0.2538399260252846e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6578811126669331
       B=0.2003106382156076
       V=0.2546255927268069e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5609624612998100
       B=0.3337500940231335e-1
       V=0.2500583360048449e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5979959659984670
       B=0.6708750335901803e-1
       V=0.2524777638260203e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6330523711054002
       B=0.1008792126424850
       V=0.2540951193860656e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6660960998103972
       B=0.1345050343171794
       V=0.2549524085027472e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6365384364585819
       B=0.3372799460737052e-1
       V=0.2542569507009158e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6710994302899275
       B=0.6755249309678028e-1
       V=0.2552114127580376e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD4802(X,Y,Z,W,N)
       DOUBLE PRECISION X(4802)
       DOUBLE PRECISION Y(4802)
       DOUBLE PRECISION Z(4802)
       DOUBLE PRECISION W(4802)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 4802-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.9687521879420705e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2307897895367918e-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2297310852498558e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2335728608887064e-1
       V=0.7386265944001919e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4352987836550653e-1
       V=0.8257977698542210e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6439200521088801e-1
       V=0.9706044762057630e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9003943631993181e-1
       V=0.1302393847117003e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1196706615548473
       V=0.1541957004600968e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1511715412838134
       V=0.1704459770092199e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1835982828503801
       V=0.1827374890942906e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2165081259155405
       V=0.1926360817436107e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2496208720417563
       V=0.2008010239494833e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2827200673567900
       V=0.2075635983209175e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3156190823994346
       V=0.2131306638690909e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3481476793749115
       V=0.2176562329937335e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3801466086947226
       V=0.2212682262991018e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4114652119634011
       V=0.2240799515668565e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4419598786519751
       V=0.2261959816187525e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4714925949329543
       V=0.2277156368808855e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4999293972879466
       V=0.2287351772128336e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5271387221431248
       V=0.2293490814084085e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5529896780837761
       V=0.2296505312376273e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6000856099481712
       V=0.2296793832318756e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6210562192785175
       V=0.2295785443842974e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6401165879934240
       V=0.2295017931529102e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6571144029244334
       V=0.2295059638184868e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718910821718863
       V=0.2296232343237362e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6842845591099010
       V=0.2298530178740771e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6941353476269816
       V=0.2301579790280501e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012965242212991
       V=0.2304690404996513e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7056471428242644
       V=0.2307027995907102e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4595557643585895e-1
       V=0.9312274696671092e-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1049316742435023
       V=0.1199919385876926e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1773548879549274
       V=0.1598039138877690e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2559071411236127
       V=0.1822253763574900e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3358156837985898
       V=0.1988579593655040e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4155835743763893
       V=0.2112620102533307e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4937894296167472
       V=0.2201594887699007e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5691569694793316
       V=0.2261622590895036e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6405840854894251
       V=0.2296458453435705e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7345133894143348e-1
       B=0.2177844081486067e-1
       V=0.1006006990267000e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1009859834044931
       B=0.4590362185775188e-1
       V=0.1227676689635876e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1324289619748758
       B=0.7255063095690877e-1
       V=0.1467864280270117e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1654272109607127
       B=0.1017825451960684
       V=0.1644178912101232e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1990767186776461
       B=0.1325652320980364
       V=0.1777664890718961e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2330125945523278
       B=0.1642765374496765
       V=0.1884825664516690e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2670080611108287
       B=0.1965360374337889
       V=0.1973269246453848e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3008753376294316
       B=0.2290726770542238
       V=0.2046767775855328e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3344475596167860
       B=0.2616645495370823
       V=0.2107600125918040e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3675709724070786
       B=0.2941150728843141
       V=0.2157416362266829e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4001000887587812
       B=0.3262440400919066
       V=0.2197557816920721e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318956350436028
       B=0.3578835350611916
       V=0.2229192611835437e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4628239056795531
       B=0.3888751854043678
       V=0.2253385110212775e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4927563229773636
       B=0.4190678003222840
       V=0.2271137107548774e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5215687136707969
       B=0.4483151836883852
       V=0.2283414092917525e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5491402346984905
       B=0.4764740676087880
       V=0.2291161673130077e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5753520160126075
       B=0.5034021310998277
       V=0.2295313908576598e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1388326356417754
       B=0.2435436510372806e-1
       V=0.1438204721359031e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1743686900537244
       B=0.5118897057342652e-1
       V=0.1607738025495257e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2099737037950268
       B=0.8014695048539634e-1
       V=0.1741483853528379e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2454492590908548
       B=0.1105117874155699
       V=0.1851918467519151e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2807219257864278
       B=0.1417950531570966
       V=0.1944628638070613e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3156842271975842
       B=0.1736604945719597
       V=0.2022495446275152e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3502090945177752
       B=0.2058466324693981
       V=0.2087462382438514e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3841684849519686
       B=0.2381284261195919
       V=0.2141074754818308e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4174372367906016
       B=0.2703031270422569
       V=0.2184640913748162e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4498926465011892
       B=0.3021845683091309
       V=0.2219309165220329e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4814146229807701
       B=0.3335993355165720
       V=0.2246123118340624e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5118863625734701
       B=0.3643833735518232
       V=0.2266062766915125e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5411947455119144
       B=0.3943789541958179
       V=0.2280072952230796e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5692301500357246
       B=0.4234320144403542
       V=0.2289082025202583e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5958857204139576
       B=0.4513897947419260
       V=0.2294012695120025e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2156270284785766
       B=0.2681225755444491e-1
       V=0.1722434488736947e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2532385054909710
       B=0.5557495747805614e-1
       V=0.1830237421455091e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2902564617771537
       B=0.8569368062950249e-1
       V=0.1923855349997633e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3266979823143256
       B=0.1167367450324135
       V=0.2004067861936271e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3625039627493614
       B=0.1483861994003304
       V=0.2071817297354263e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3975838937548699
       B=0.1803821503011405
       V=0.2128250834102103e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318396099009774
       B=0.2124962965666424
       V=0.2174513719440102e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4651706555732742
       B=0.2445221837805913
       V=0.2211661839150214e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4974752649620969
       B=0.2762701224322987
       V=0.2240665257813102e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5286517579627517
       B=0.3075627775211328
       V=0.2262439516632620e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5586001195731895
       B=0.3382311089826877
       V=0.2277874557231869e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5872229902021319
       B=0.3681108834741399
       V=0.2287854314454994e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6144258616235123
       B=0.3970397446872839
       V=0.2293268499615575e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2951676508064861
       B=0.2867499538750441e-1
       V=0.1912628201529828e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3335085485472725
       B=0.5867879341903510e-1
       V=0.1992499672238701e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3709561760636381
       B=0.8961099205022284e-1
       V=0.2061275533454027e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4074722861667498
       B=0.1211627927626297
       V=0.2119318215968572e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4429923648839117
       B=0.1530748903554898
       V=0.2167416581882652e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4774428052721736
       B=0.1851176436721877
       V=0.2206430730516600e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5107446539535904
       B=0.2170829107658179
       V=0.2237186938699523e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5428151370542935
       B=0.2487786689026271
       V=0.2260480075032884e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5735699292556964
       B=0.2800239952795016
       V=0.2277098884558542e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6029253794562866
       B=0.3106445702878119
       V=0.2287845715109671e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6307998987073145
       B=0.3404689500841194
       V=0.2293547268236294e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3752652273692719
       B=0.2997145098184479e-1
       V=0.2056073839852528e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4135383879344028
       B=0.6086725898678011e-1
       V=0.2114235865831876e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4506113885153907
       B=0.9238849548435643e-1
       V=0.2163175629770551e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4864401554606072
       B=0.1242786603851851
       V=0.2203392158111650e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5209708076611709
       B=0.1563086731483386
       V=0.2235473176847839e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5541422135830122
       B=0.1882696509388506
       V=0.2260024141501235e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5858880915113817
       B=0.2199672979126059
       V=0.2277675929329182e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6161399390603444
       B=0.2512165482924867
       V=0.2289102112284834e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6448296482255090
       B=0.2818368701871888
       V=0.2295027954625118e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4544796274917948
       B=0.3088970405060312e-1
       V=0.2161281589879992e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4919389072146628
       B=0.6240947677636835e-1
       V=0.2201980477395102e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5279313026985183
       B=0.9430706144280313e-1
       V=0.2234952066593166e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5624169925571135
       B=0.1263547818770374
       V=0.2260540098520838e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5953484627093287
       B=0.1583430788822594
       V=0.2279157981899988e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6266730715339185
       B=0.1900748462555988
       V=0.2291296918565571e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6563363204278871
       B=0.2213599519592567
       V=0.2297533752536649e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5314574716585696
       B=0.3152508811515374e-1
       V=0.2234927356465995e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5674614932298185
       B=0.6343865291465561e-1
       V=0.2261288012985219e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6017706004970264
       B=0.9551503504223951e-1
       V=0.2280818160923688e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6343471270264178
       B=0.1275440099801196
       V=0.2293773295180159e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6651494599127802
       B=0.1593252037671960
       V=0.2300528767338634e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6050184986005704
       B=0.3192538338496105e-1
       V=0.2281893855065666e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6390163550880400
       B=0.6402824353962306e-1
       V=0.2295720444840727e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6711199107088448
       B=0.9609805077002909e-1
       V=0.2303227649026753e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6741354429572275
       B=0.3211853196273233e-1
       V=0.2304831913227114e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD5294(X,Y,Z,W,N)
       DOUBLE PRECISION X(5294)
       DOUBLE PRECISION Y(5294)
       DOUBLE PRECISION Z(5294)
       DOUBLE PRECISION W(5294)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 5294-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.9080510764308163e-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2084824361987793e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2303261686261450e-1
       V=0.5011105657239616e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3757208620162394e-1
       V=0.5942520409683854e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5821912033821852e-1
       V=0.9564394826109721e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8403127529194872e-1
       V=0.1185530657126338e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1122927798060578
       V=0.1364510114230331e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1420125319192987
       V=0.1505828825605415e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1726396437341978
       V=0.1619298749867023e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2038170058115696
       V=0.1712450504267789e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2352849892876508
       V=0.1789891098164999e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2668363354312461
       V=0.1854474955629795e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2982941279900452
       V=0.1908148636673661e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3295002922087076
       V=0.1952377405281833e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3603094918363593
       V=0.1988349254282232e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3905857895173920
       V=0.2017079807160050e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4202005758160837
       V=0.2039473082709094e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4490310061597227
       V=0.2056360279288953e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4769586160311491
       V=0.2068525823066865e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5038679887049750
       V=0.2076724877534488e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5296454286519961
       V=0.2081694278237885e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5541776207164850
       V=0.2084157631219326e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5990467321921213
       V=0.2084381531128593e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6191467096294587
       V=0.2083476277129307e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6375251212901849
       V=0.2082686194459732e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6540514381131168
       V=0.2082475686112415e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6685899064391510
       V=0.2083139860289915e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6810013009681648
       V=0.2084745561831237e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6911469578730340
       V=0.2087091313375890e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6988956915141736
       V=0.2089718413297697e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7041335794868720
       V=0.2092003303479793e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7067754398018567
       V=0.2093336148263241e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3840368707853623e-1
       V=0.7591708117365267e-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9835485954117399e-1
       V=0.1083383968169186e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1665774947612998
       V=0.1403019395292510e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2405702335362910
       V=0.1615970179286436e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3165270770189046
       V=0.1771144187504911e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3927386145645443
       V=0.1887760022988168e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4678825918374656
       V=0.1973474670768214e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5408022024266935
       V=0.2033787661234659e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6104967445752438
       V=0.2072343626517331e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6760910702685738
       V=0.2091177834226918e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6655644120217392e-1
       B=0.1936508874588424e-1
       V=0.9316684484675566e-4
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9446246161270182e-1
       B=0.4252442002115869e-1
       V=0.1116193688682976e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1242651925452509
       B=0.6806529315354374e-1
       V=0.1298623551559414e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1553438064846751
       B=0.9560957491205369e-1
       V=0.1450236832456426e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1871137110542670
       B=0.1245931657452888
       V=0.1572719958149914e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2192612628836257
       B=0.1545385828778978
       V=0.1673234785867195e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2515682807206955
       B=0.1851004249723368
       V=0.1756860118725188e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2838535866287290
       B=0.2160182608272384
       V=0.1826776290439367e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3159578817528521
       B=0.2470799012277111
       V=0.1885116347992865e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3477370882791392
       B=0.2781014208986402
       V=0.1933457860170574e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3790576960890540
       B=0.3089172523515731
       V=0.1973060671902064e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4097938317810200
       B=0.3393750055472244
       V=0.2004987099616311e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4398256572859637
       B=0.3693322470987730
       V=0.2030170909281499e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4690384114718480
       B=0.3986541005609877
       V=0.2049461460119080e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4973216048301053
       B=0.4272112491408562
       V=0.2063653565200186e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5245681526132446
       B=0.4548781735309936
       V=0.2073507927381027e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5506733911803888
       B=0.4815315355023251
       V=0.2079764593256122e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5755339829522475
       B=0.5070486445801855
       V=0.2083150534968778e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1305472386056362
       B=0.2284970375722366e-1
       V=0.1262715121590664e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1637327908216477
       B=0.4812254338288384e-1
       V=0.1414386128545972e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1972734634149637
       B=0.7531734457511935e-1
       V=0.1538740401313898e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2308694653110130
       B=0.1039043639882017
       V=0.1642434942331432e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2643899218338160
       B=0.1334526587117626
       V=0.1729790609237496e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2977171599622171
       B=0.1636414868936382
       V=0.1803505190260828e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3307293903032310
       B=0.1942195406166568
       V=0.1865475350079657e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3633069198219073
       B=0.2249752879943753
       V=0.1917182669679069e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3953346955922727
       B=0.2557218821820032
       V=0.1959851709034382e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4267018394184914
       B=0.2862897925213193
       V=0.1994529548117882e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4573009622571704
       B=0.3165224536636518
       V=0.2022138911146548e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4870279559856109
       B=0.3462730221636496
       V=0.2043518024208592e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5157819581450322
       B=0.3754016870282835
       V=0.2059450313018110e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5434651666465393
       B=0.4037733784993613
       V=0.2070685715318472e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5699823887764627
       B=0.4312557784139123
       V=0.2077955310694373e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5952403350947741
       B=0.4577175367122110
       V=0.2081980387824712e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2025152599210369
       B=0.2520253617719557e-1
       V=0.1521318610377956e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2381066653274425
       B=0.5223254506119000e-1
       V=0.1622772720185755e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2732823383651612
       B=0.8060669688588620e-1
       V=0.1710498139420709e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3080137692611118
       B=0.1099335754081255
       V=0.1785911149448736e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3422405614587601
       B=0.1399120955959857
       V=0.1850125313687736e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3758808773890420
       B=0.1702977801651705
       V=0.1904229703933298e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4088458383438932
       B=0.2008799256601680
       V=0.1949259956121987e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4410450550841152
       B=0.2314703052180836
       V=0.1986161545363960e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4723879420561312
       B=0.2618972111375892
       V=0.2015790585641370e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5027843561874343
       B=0.2920013195600270
       V=0.2038934198707418e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5321453674452458
       B=0.3216322555190551
       V=0.2056334060538251e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5603839113834030
       B=0.3506456615934198
       V=0.2068705959462289e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5874150706875146
       B=0.3789007181306267
       V=0.2076753906106002e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6131559381660038
       B=0.4062580170572782
       V=0.2081179391734803e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2778497016394506
       B=0.2696271276876226e-1
       V=0.1700345216228943e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3143733562261912
       B=0.5523469316960465e-1
       V=0.1774906779990410e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3501485810261827
       B=0.8445193201626464e-1
       V=0.1839659377002642e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3851430322303653
       B=0.1143263119336083
       V=0.1894987462975169e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4193013979470415
       B=0.1446177898344475
       V=0.1941548809452595e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4525585960458567
       B=0.1751165438438091
       V=0.1980078427252384e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4848447779622947
       B=0.2056338306745660
       V=0.2011296284744488e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5160871208276894
       B=0.2359965487229226
       V=0.2035888456966776e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5462112185696926
       B=0.2660430223139146
       V=0.2054516325352142e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5751425068101757
       B=0.2956193664498032
       V=0.2067831033092635e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6028073872853596
       B=0.3245763905312779
       V=0.2076485320284876e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6291338275278409
       B=0.3527670026206972
       V=0.2081141439525255e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3541797528439391
       B=0.2823853479435550e-1
       V=0.1834383015469222e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3908234972074657
       B=0.5741296374713106e-1
       V=0.1889540591777677e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4264408450107590
       B=0.8724646633650199e-1
       V=0.1936677023597375e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4609949666553286
       B=0.1175034422915616
       V=0.1976176495066504e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4944389496536006
       B=0.1479755652628428
       V=0.2008536004560983e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5267194884346086
       B=0.1784740659484352
       V=0.2034280351712291e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5577787810220990
       B=0.2088245700431244
       V=0.2053944466027758e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5875563763536670
       B=0.2388628136570763
       V=0.2068077642882360e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6159910016391269
       B=0.2684308928769185
       V=0.2077250949661599e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6430219602956268
       B=0.2973740761960252
       V=0.2082062440705320e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4300647036213646
       B=0.2916399920493977e-1
       V=0.1934374486546626e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4661486308935531
       B=0.5898803024755659e-1
       V=0.1974107010484300e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5009658555287261
       B=0.8924162698525409e-1
       V=0.2007129290388658e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5344824270447704
       B=0.1197185199637321
       V=0.2033736947471293e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5666575997416371
       B=0.1502300756161382
       V=0.2054287125902493e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5974457471404752
       B=0.1806004191913564
       V=0.2069184936818894e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6267984444116886
       B=0.2106621764786252
       V=0.2078883689808782e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6546664713575417
       B=0.2402526932671914
       V=0.2083886366116359e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5042711004437253
       B=0.2982529203607657e-1
       V=0.2006593275470817e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5392127456774380
       B=0.6008728062339922e-1
       V=0.2033728426135397e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726819437668618
       B=0.9058227674571398e-1
       V=0.2055008781377608e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6046469254207278
       B=0.1211219235803400
       V=0.2070651783518502e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6350716157434952
       B=0.1515286404791580
       V=0.2080953335094320e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6639177679185454
       B=0.1816314681255552
       V=0.2086284998988521e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5757276040972253
       B=0.3026991752575440e-1
       V=0.2055549387644668e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6090265823139755
       B=0.6078402297870770e-1
       V=0.2071871850267654e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6406735344387661
       B=0.9135459984176636e-1
       V=0.2082856600431965e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6706397927793709
       B=0.1218024155966590
       V=0.2088705858819358e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6435019674426665
       B=0.3052608357660639e-1
       V=0.2083995867536322e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6747218676375681
       B=0.6112185773983089e-1
       V=0.2090509712889637e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD5810(X,Y,Z,W,N)
       DOUBLE PRECISION X(5810)
       DOUBLE PRECISION Y(5810)
       DOUBLE PRECISION Z(5810)
       DOUBLE PRECISION W(5810)
       INTEGER N
       DOUBLE PRECISION A,B,V
!  
!      LEBEDEV 5810-POINT ANGULAR GRID
!  
!    
!       This subroutine is part of a set of subroutines that generate
!       Lebedev grids [1-6] for integration on a sphere. The original 
!       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!       translated into fortran by Dr. Christoph van Wuellen.
!       This subroutine was translated using a C to fortran77 conversion
!       tool written by Dr. Christoph van Wuellen.
!    
!       Users of this code are asked to include reference [1] in their
!       publications, and in the user- and programmers-manuals 
!       describing their codes.
!    
!       This code was distributed through CCL (http://www.ccl.net/).
!    
!       [1] V.I. Lebedev, and D.N. Laikov
!           "A quadrature formula for the sphere of the 131st
!            algebraic order of accuracy"
!           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!    
!       [2] V.I. Lebedev
!           "A quadrature formula for the sphere of 59th algebraic
!            order of accuracy"
!           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!    
!       [3] V.I. Lebedev, and A.L. Skorokhodov
!           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!    
!       [4] V.I. Lebedev
!           "Spherical quadrature formulas exact to orders 25-29"
!           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!    
!       [5] V.I. Lebedev
!           "Quadratures on a sphere"
!           Computational Mathematics and Mathematical Physics, Vol. 16,
!           1976, pp. 10-24. 
!    
!       [6] V.I. Lebedev
!           "Values of the nodes and weights of ninth to seventeenth 
!            order Gauss-Markov quadrature formulae invariant under the
!            octahedron group with inversion"
!           Computational Mathematics and Mathematical Physics, Vol. 15,
!           1975, pp. 44-51.
!    
       N=1
       V=0.9735347946175486e-5
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1907581241803167e-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1901059546737578e-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1182361662400277e-1
       V=0.3926424538919212e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3062145009138958e-1
       V=0.6667905467294382e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5329794036834243e-1
       V=0.8868891315019135e-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7848165532862220e-1
       V=0.1066306000958872e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1054038157636201
       V=0.1214506743336128e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1335577797766211
       V=0.1338054681640871e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1625769955502252
       V=0.1441677023628504e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1921787193412792
       V=0.1528880200826557e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2221340534690548
       V=0.1602330623773609e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2522504912791132
       V=0.1664102653445244e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2823610860679697
       V=0.1715845854011323e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3123173966267560
       V=0.1758901000133069e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3419847036953789
       V=0.1794382485256736e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3712386456999758
       V=0.1823238106757407e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3999627649876828
       V=0.1846293252959976e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4280466458648093
       V=0.1864284079323098e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4553844360185711
       V=0.1877882694626914e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4818736094437834
       V=0.1887716321852025e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5074138709260629
       V=0.1894381638175673e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5319061304570707
       V=0.1898454899533629e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5552514978677286
       V=0.1900497929577815e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5981009025246183
       V=0.1900671501924092e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6173990192228116
       V=0.1899837555533510e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6351365239411131
       V=0.1899014113156229e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6512010228227200
       V=0.1898581257705106e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6654758363948120
       V=0.1898804756095753e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6778410414853370
       V=0.1899793610426402e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6881760887484110
       V=0.1901464554844117e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6963645267094598
       V=0.1903533246259542e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7023010617153579
       V=0.1905556158463228e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7059004636628753
       V=0.1907037155663528e-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3552470312472575e-1
       V=0.5992997844249967e-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9151176620841283e-1
       V=0.9749059382456978e-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1566197930068980
       V=0.1241680804599158e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2265467599271907
       V=0.1437626154299360e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2988242318581361
       V=0.1584200054793902e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3717482419703886
       V=0.1694436550982744e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4440094491758889
       V=0.1776617014018108e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5145337096756642
       V=0.1836132434440077e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5824053672860230
       V=0.1876494727075983e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6468283961043370
       V=0.1899906535336482e-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6095964259104373e-1
       B=0.1787828275342931e-1
       V=0.8143252820767350e-4
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8811962270959388e-1
       B=0.3953888740792096e-1
       V=0.9998859890887728e-4
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1165936722428831
       B=0.6378121797722990e-1
       V=0.1156199403068359e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1460232857031785
       B=0.8985890813745037e-1
       V=0.1287632092635513e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1761197110181755
       B=0.1172606510576162
       V=0.1398378643365139e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2066471190463718
       B=0.1456102876970995
       V=0.1491876468417391e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2374076026328152
       B=0.1746153823011775
       V=0.1570855679175456e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2682305474337051
       B=0.2040383070295584
       V=0.1637483948103775e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2989653312142369
       B=0.2336788634003698
       V=0.1693500566632843e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3294762752772209
       B=0.2633632752654219
       V=0.1740322769393633e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3596390887276086
       B=0.2929369098051601
       V=0.1779126637278296e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3893383046398812
       B=0.3222592785275512
       V=0.1810908108835412e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4184653789358347
       B=0.3512004791195743
       V=0.1836529132600190e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4469172319076166
       B=0.3796385677684537
       V=0.1856752841777379e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4745950813276976
       B=0.4074575378263879
       V=0.1872270566606832e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5014034601410262
       B=0.4345456906027828
       V=0.1883722645591307e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5272493404551239
       B=0.4607942515205134
       V=0.1891714324525297e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5520413051846366
       B=0.4860961284181720
       V=0.1896827480450146e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5756887237503077
       B=0.5103447395342790
       V=0.1899628417059528e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1225039430588352
       B=0.2136455922655793e-1
       V=0.1123301829001669e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1539113217321372
       B=0.4520926166137188e-1
       V=0.1253698826711277e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1856213098637712
       B=0.7086468177864818e-1
       V=0.1366266117678531e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2174998728035131
       B=0.9785239488772918e-1
       V=0.1462736856106918e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2494128336938330
       B=0.1258106396267210
       V=0.1545076466685412e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2812321562143480
       B=0.1544529125047001
       V=0.1615096280814007e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3128372276456111
       B=0.1835433512202753
       V=0.1674366639741759e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3441145160177973
       B=0.2128813258619585
       V=0.1724225002437900e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3749567714853510
       B=0.2422913734880829
       V=0.1765810822987288e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4052621732015610
       B=0.2716163748391453
       V=0.1800104126010751e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4349335453522385
       B=0.3007127671240280
       V=0.1827960437331284e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4638776641524965
       B=0.3294470677216479
       V=0.1850140300716308e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4920046410462687
       B=0.3576932543699155
       V=0.1867333507394938e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5192273554861704
       B=0.3853307059757764
       V=0.1880178688638289e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5454609081136522
       B=0.4122425044452694
       V=0.1889278925654758e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5706220661424140
       B=0.4383139587781027
       V=0.1895213832507346e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5946286755181518
       B=0.4634312536300553
       V=0.1898548277397420e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1905370790924295
       B=0.2371311537781979e-1
       V=0.1349105935937341e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2242518717748009
       B=0.4917878059254806e-1
       V=0.1444060068369326e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2577190808025936
       B=0.7595498960495142e-1
       V=0.1526797390930008e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2908724534927187
       B=0.1036991083191100
       V=0.1598208771406474e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3236354020056219
       B=0.1321348584450234
       V=0.1659354368615331e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3559267359304543
       B=0.1610316571314789
       V=0.1711279910946440e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3876637123676956
       B=0.1901912080395707
       V=0.1754952725601440e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4187636705218842
       B=0.2194384950137950
       V=0.1791247850802529e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4491449019883107
       B=0.2486155334763858
       V=0.1820954300877716e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4787270932425445
       B=0.2775768931812335
       V=0.1844788524548449e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5074315153055574
       B=0.3061863786591120
       V=0.1863409481706220e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5351810507738336
       B=0.3343144718152556
       V=0.1877433008795068e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5619001025975381
       B=0.3618362729028427
       V=0.1887444543705232e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5875144035268046
       B=0.3886297583620408
       V=0.1894009829375006e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6119507308734495
       B=0.4145742277792031
       V=0.1897683345035198e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2619733870119463
       B=0.2540047186389353e-1
       V=0.1517327037467653e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2968149743237949
       B=0.5208107018543989e-1
       V=0.1587740557483543e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3310451504860488
       B=0.7971828470885599e-1
       V=0.1649093382274097e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3646215567376676
       B=0.1080465999177927
       V=0.1701915216193265e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3974916785279360
       B=0.1368413849366629
       V=0.1746847753144065e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4295967403772029
       B=0.1659073184763559
       V=0.1784555512007570e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4608742854473447
       B=0.1950703730454614
       V=0.1815687562112174e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4912598858949903
       B=0.2241721144376724
       V=0.1840864370663302e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5206882758945558
       B=0.2530655255406489
       V=0.1860676785390006e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5490940914019819
       B=0.2816118409731066
       V=0.1875690583743703e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5764123302025542
       B=0.3096780504593238
       V=0.1886453236347225e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6025786004213506
       B=0.3371348366394987
       V=0.1893501123329645e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6275291964794956
       B=0.3638547827694396
       V=0.1897366184519868e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3348189479861771
       B=0.2664841935537443e-1
       V=0.1643908815152736e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3699515545855295
       B=0.5424000066843495e-1
       V=0.1696300350907768e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4042003071474669
       B=0.8251992715430854e-1
       V=0.1741553103844483e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4375320100182624
       B=0.1112695182483710
       V=0.1780015282386092e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4699054490335947
       B=0.1402964116467816
       V=0.1812116787077125e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5012739879431952
       B=0.1694275117584291
       V=0.1838323158085421e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5315874883754966
       B=0.1985038235312689
       V=0.1859113119837737e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5607937109622117
       B=0.2273765660020893
       V=0.1874969220221698e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5888393223495521
       B=0.2559041492849764
       V=0.1886375612681076e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6156705979160163
       B=0.2839497251976899
       V=0.1893819575809276e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6412338809078123
       B=0.3113791060500690
       V=0.1897794748256767e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4076051259257167
       B=0.2757792290858463e-1
       V=0.1738963926584846e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4423788125791520
       B=0.5584136834984293e-1
       V=0.1777442359873466e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4760480917328258
       B=0.8457772087727143e-1
       V=0.1810010815068719e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5085838725946297
       B=0.1135975846359248
       V=0.1836920318248129e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5399513637391218
       B=0.1427286904765053
       V=0.1858489473214328e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5701118433636380
       B=0.1718112740057635
       V=0.1875079342496592e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5990240530606021
       B=0.2006944855985351
       V=0.1887080239102310e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6266452685139695
       B=0.2292335090598907
       V=0.1894905752176822e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6529320971415942
       B=0.2572871512353714
       V=0.1898991061200695e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4791583834610126
       B=0.2826094197735932e-1
       V=0.1809065016458791e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5130373952796940
       B=0.5699871359683649e-1
       V=0.1836297121596799e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5456252429628476
       B=0.8602712528554394e-1
       V=0.1858426916241869e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5768956329682385
       B=0.1151748137221281
       V=0.1875654101134641e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6068186944699046
       B=0.1442811654136362
       V=0.1888240751833503e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6353622248024907
       B=0.1731930321657680
       V=0.1896497383866979e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6624927035731797
       B=0.2017619958756061
       V=0.1900775530219121e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5484933508028488
       B=0.2874219755907391e-1
       V=0.1858525041478814e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5810207682142106
       B=0.5778312123713695e-1
       V=0.1876248690077947e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6120955197181352
       B=0.8695262371439526e-1
       V=0.1889404439064607e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6416944284294319
       B=0.1160893767057166
       V=0.1898168539265290e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6697926391731260
       B=0.1450378826743251
       V=0.1902779940661772e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6147594390585488
       B=0.2904957622341456e-1
       V=0.1890125641731815e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6455390026356783
       B=0.5823809152617197e-1
       V=0.1899434637795751e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6747258588365477
       B=0.8740384899884715e-1
       V=0.1904520856831751e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6772135750395347
       B=0.2919946135808105e-1
       V=0.1905534498734563e-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END




end module m_corespec_eval
