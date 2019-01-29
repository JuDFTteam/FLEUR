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
    csv%gaunt = 0.d0

! m<=l condition fulfilled by looping m within l value interval {-l,...,+l}

    csv%gaunt = 0.d0
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
                       &sqrt(twol1p1*twola1p1/(4.d0*pi_const*twolip1))*&
                       &(-1)**(mi)
                  if(csv%gaunt(l1,m1,la1,mu1,li,mi).ne.0.d0) &
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
    csv%fc = 0.d0

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
    csv%fv = 0.d0

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
    csv%fb = 0.d0

    do ie = nen,nex
      do iqv = 1,nqv
        do iljc = 1,csv%nljc
          do ir = 1,nr
            fsb=0.d0
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
  subroutine corespec_dos(atoms,usdus,ispin,nkpt,ikpt,&
                          neigd,noccbd,efermi,sig_dos,eig,we,eigVecCoeffs)

    IMPLICIT NONE

    TYPE (t_atoms), INTENT(IN)      :: atoms
    TYPE (t_usdus), INTENT(IN)      :: usdus
    TYPE(t_eigVecCoeffs),INTENT(IN) :: eigVecCoeffs

!     .. Scalar Arguments ..
    integer, intent(in) :: ispin,nkpt,ikpt
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
    sigma = sqrt(2.d0)*sig_dos*hartree_to_ev_const
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
      csv%dose = 0.d0
    endif
    csv%dosb = 0.d0

    do iband = 1,noccbd
      do l1 = 0,lx
        do m1 = -l1,l1
          lm1 = l1*(l1+1)+m1
!!$          do l2 = 0,lx
!!$            do m2 = -l2,l2
!!$              lm2 = l2*(l2+1)+m2
          csv%dosb(1,1,lm1,lm1,iband) = dble(eigVecCoeffs%acof(iband,lm1,iatom,ispin)*&
               &conjg(eigVecCoeffs%acof(iband,lm1,iatom,ispin)))!*we(1)
          csv%dosb(1,2,lm1,lm1,iband) = dble(eigVecCoeffs%acof(iband,lm1,iatom,ispin)*&
               &conjg(eigVecCoeffs%bcof(iband,lm1,iatom,ispin)))
          csv%dosb(2,1,lm1,lm1,iband) = dble(eigVecCoeffs%bcof(iband,lm1,iatom,ispin)*&
               &conjg(eigVecCoeffs%acof(iband,lm1,iatom,ispin)))
          csv%dosb(2,2,lm1,lm1,iband) = dble(eigVecCoeffs%bcof(iband,lm1,iatom,ispin)*&
               &conjg(eigVecCoeffs%bcof(iband,lm1,iatom,ispin)))!*we(1)*usdus%ddn(l1,csi%atomType,ispin)
!!!!! this has to be checked: is >> ddn << factor necessary !!!!!
!!$        enddo
!!$        enddo
        enddo
      enddo
      if(eigos(iband)+3.d0*sigma.ge.csv%eos(0).and.&
           &eigos(iband)-3.d0*sigma.le.csv%eos(nex)) then
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
      WRITE(*,'(10i8)') atoms%lmaxd*(atoms%lmaxd+2),atoms%ntype

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

    real :: gamma,beta,rho,qepref,orvec(3)
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
    orvec = (/1.d0,0.d0,0.d0/)

    if(.not.allocated(csv%ddscs)) then
      allocate(csv%ddscs(2,0:nor,1:csv%nljc,nqv,0:nex))
      csv%ddscs = cmplx(0.d0,0.d0)
    endif
    if(.not.allocated(tdy)) allocate(tdy(0:nor,2))
    if(.not.allocated(orpref)) then
      allocate(orpref(0:nor))
      orpref(0) = 1.d0
      if(nor.gt.0) then
        orpref(1:nor) = (4.d0*pi_const)**2
        if(.not.allocated(ylm)) allocate(ylm(0:lax*(lax+2),nor))
        do ior = 1,nor
          CALL ylm4(lax,orvec,ylm(:,ior))
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

    rho = alpha*beta*sqrt(4.d0*pi_const/3.d0)
    print*,gamma,beta,rho
!    rho = 0.d0

    do ie = nen,nex  ! energy
      do iqv = 1,nqv  ! q-vector
        do iljc = 1,csv%nljc  ! core levels
          li = edgel(csi%edgeidx(iljc))
          qepref = 4.d0*gamma**2*csv%qv1(iljc,iqv,ie)/csv%qv0/&
               &(csv%qv(0,iljc,iqv,ie)**2-(csv%eloss(iljc,ie)*alpha)**2)**2
!!$          write(*,'(2i5,3f20.4)') ie,iljc,csv%qv(0,iljc,iqv,ie),csv%eloss(iljc,ie)*alpha,qepref
          tdy = cmplx(0.d0,0.d0)

          do imi = 1,nint(csv%occ(iljc)*jspins/2)
            mi = sign(jspin)*(edgej(csi%edgeidx(iljc))-4*(imi-1)-1)/2
            print*,jspin,ie,iljc,li,mi
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
                    
                    if(l1.eq.l2.and.m1.eq.m2) then

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

                    prd = 0.d0

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
                          &+cimu*rho*(&
                          &-prd(0,1)*ga(0,1)*ga(1,2)&
                          &+prd(0,2)*ga(0,1)*ga(2,2)&
                          &+prd(1,0)*ga(1,1)*ga(0,2)&
                          &-prd(2,0)*ga(2,1)*ga(0,2))
                    
                    ila1la2 = cimu**(la1-la2)
                    
                    if(abs(real(td(1))).gt.0.d0.or.abs(real(td(2))).gt.0.d0.or.abs(aimag(td(2))).gt.0.d0) then
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

                  endif

                  enddo; enddo
                enddo; enddo

              enddo; enddo
            enddo; enddo

          enddo  ! mi

          do it = 1,2
            do ior = 0,nor
              csv%ddscs(it,ior,iljc,iqv,ie) = csv%ddscs(it,ior,iljc,iqv,ie)+&
                   &qepref*orpref(ior)*tdy(ior,it)
            enddo
          enddo

        enddo  ! iljc
      enddo  ! iqv
    enddo  ! ie

    if(jspin.eq.1) then
      do ior = 0,nor
      do iljc = 1,csv%nljc
        do ie = nen,nex
          write(37,'(2i5,f8.3,4es16.4)') ior,iljc,csv%eloss(iljc,ie)*hartree_to_ev_const,csv%ddscs(1,ior,iljc,1,ie),csv%ddscs(2,ior,iljc,1,ie)
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
          write(38,'(2i5,f8.3,4es16.4)') ior,iljc,csv%eloss(iljc,ie)*hartree_to_ev_const,csv%ddscs(1,ior,iljc,1,ie),csv%ddscs(2,ior,iljc,1,ie)
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

    integer :: ie,iljc,iqv
    real :: eout,relfac

    smeno = "corespec_eloss_qv"

    write(*,'(/,a)') trim(smeno)//ssep

    csv%nqv = 1

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
!!$    print*,csi%ek0,csv%qv0

    if(.not.allocated(csv%qv)) &
         &allocate(csv%qv(0:3,csv%nljc,csv%nqv,csv%nen:csv%nex))
    csv%qv=0.d0
    do ie = csv%nen,csv%nex
      do iqv = 1,csv%nqv
        do iljc = 1,csv%nljc
          eout = csi%ek0/hartree_to_ev_const-csv%eloss(iljc,ie)
          csv%qv1(iljc,iqv,ie) = e2q(eout)
          csv%qv(3,iljc,iqv,ie) = (csv%qv1(iljc,iqv,ie)-csv%qv0)*relfac
          csv%qv(0,iljc,iqv,ie) = sqrt(&
               &dot_product(csv%qv(1:3,iljc,iqv,ie),csv%qv(1:3,iljc,iqv,ie)))
!!$          write(*,'(3i5,2f16.2,6f16.6)') ie,iqv,iljc,csi%ek0,eout*hartree_to_ev_const,csv%qv1(iljc,iqv,ie),csv%eloss(iljc,ie),csv%qv(:,iljc,iqv,ie)
        enddo
      enddo
    enddo

    if(csi%verb.eq.1) write(*,*) ""

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

    e2q=sqrt(e**2+2.d0*e*mec2/hartree_to_ev_const)*alpha

  end function e2q
!
!===============================================================================


end module m_corespec_eval
