!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_corespec_io

  USE m_corespec
  USE m_types
  USE m_juDFT

  IMPLICIT NONE

  CONTAINS

!===============================================================================
!
!  S U B R O U T I N E   C O R E S P E C _ I N I T
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE corespec_init(input,atoms,coreSpecInput)

    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)         :: input
    TYPE(t_atoms),INTENT(IN)         :: atoms
    TYPE(t_coreSpecInput),INTENT(IN) :: coreSpecInput

    INTEGER                    :: ui,i
    LOGICAL                    :: lexist

    namelist /csinp/ csi

    l_cs = .false.

! initialization of input parameters: type csi

    csi%verb = 1
    csi%atomType = 0
    csi%edge = ""
    csi%edgeidx = 0
    csi%lx = -1
    csi%ek0 = 0.d0
    csi%emn = -2.d0
    csi%emx = 20.d0
    csi%ein = 0.1d0

    csi%nqphi = 12
    csi%nqr = 20
    csi%alpha_ex=0.1 
    csi%beta_ex=0.08 
    csi%I0 = 10 

! reading of input parameters from 'corespec_inp' file

    call iounit(ui)
    inquire(file="corespec_inp",exist=lexist)
    if(lexist) then
      open(ui,file="corespec_inp",status='old')
      read(ui,nml=csinp)
      close(ui)
    else if(input%l_coreSpec) THEN
      csi = coreSpecInput
    else
      return
    endif

    smeno = "corespec_init"

    write(*,'(/,a)') trim(smeno)//ssep

    IF(ANY(atoms%nlo(:).NE.0)) THEN
       WRITE(*,*) 'PLEASE NOTE:'
       WRITE(*,*) 'EELS + LOs is only for experienced users.'
       WRITE(*,*) 'I hope you know what you are doing.'
!       CALL juDFT_error("EELS + LOs not available at the moment!" ,calledby ="corespec_io")
    END IF

! sanity check of the input parameters; if they are not correct, program stops
! unit conversion if necessary
! if csi%verb = 1, detailed information is written to stdout

!   csv%nqphi
    if(csi%nqphi.lt.0) then
      write(*,csmsgs)  trim(smeno),"found csi%nqphi < 0 !"//csmsgerr ; stop
    endif
    csv%nqphi = csi%nqphi
    if(csi%verb.eq.1) write(*,csmsgsis)  trim(smeno),&
           &"Mesh number of angles: ","csi%nqphi = ",csv%nqphi,"will be used"
 
! csi%nqr
    if(csi%nqr.lt.0) then
      write(*,csmsgs)  trim(smeno),"found csi%nqr < 0 !"//csmsgerr ; stop
    endif
    csv%nqr = csi%nqr
    if(csi%verb.eq.1) write(*,csmsgsis)  trim(smeno),&
           &"Mesh number of radii: ","csi%nqr = ",csv%nqr,"will be used"

! csi%alpha_ex
    if(csi%alpha_ex.lt.0) then
      write(*,csmsgs)  trim(smeno),"found csi%alpha_ex < 0 !"//csmsgerr ; stop
    endif
    csv%alpha_ex = csi%alpha_ex
    if(csi%verb.eq.1) write(*,csmsgsfs)  trim(smeno),&
           &"Experimental convergence angle of incoming e-: ","csi%alpha_ex = ",csv%alpha_ex,"will be used"

! csi%beta_ex
    if(csi%beta_ex.lt.0) then
      write(*,csmsgs)  trim(smeno),"found csi%beta_ex < 0 !"//csmsgerr ; stop
    endif
    csv%beta_ex = csi%beta_ex
    if(csi%verb.eq.1) write(*,csmsgsfs)  trim(smeno),&
           &"Experimental convergence angle of outgoing e-: ","csi%beta_ex = ",csv%beta_ex,"will be used"

! csi%I0
    if(csi%I0.le.0) then
      write(*,csmsgs)  trim(smeno),"found csi%I0 <= 0 !"//csmsgerr ; stop
    endif
    csv%I0 = csi%I0
    if(csi%verb.eq.1) write(*,csmsgsfs)  trim(smeno),&
           &"Intensity of incoming electrons: ","csi%I0 = ",csv%I0,"will be used"



! csi%atomType
    if(csi%atomType.le.0) then
      write(*,csmsgs)  trim(smeno),"found csi%atomType <= 0 !"//csmsgerr ; stop
    endif
    if(csi%atomType.gt.atoms%ntype) then
      write(*,csmsgs)  trim(smeno),"found csi%atomType > atoms%ntype!"//csmsgerr ; stop
    endif
    if(csi%verb.eq.1) write(*,csmsgsis)  trim(smeno),&
           &"atomic type: ","csi%atomType = ",csi%atomType,"will be used"

! csi%edge -> csv%nc
    csv%nc = 0
    if(len(trim(csi%edge)).ne.0) csv%nc = ichar(csi%edge)-ichar('K')+1
    if(csi%edge.eq."") then
      write(*,csmsgs) trim(smeno),"found empty csi%edge !"//csmsgerr ; stop
    endif
    if(csv%nc.lt.1.or.csv%nc.gt.6) then
      write(*,csmsgs) trim(smeno),&
           &"specify csi%edge as one of {K,L,M,N,O,P} !"//csmsgerr ; stop
    endif
    if(csi%verb.eq.1) write(*,csmsgsss)  trim(smeno),&
           &"edge: ","csi%edge = ",csi%edge,"will be used"
    if(csi%verb.eq.1) write(*,csmsgsis)  trim(smeno),&
           &"main quantum no.: ","csv%nc = ",csv%nc,"will be used"

! csi%edgeidx(:) -> csv%nljc
    if(maxval(csi%edgeidx).gt.(2*csv%nc-1)) then
      write(*,csmsgs) trim(smeno),&
         &"found csi%edgeidx > 2*csv%nc-1 !"//csmsgerr ; stop
    endif
    if(count(csi%edgeidx.gt.0).gt.(2*csv%nc-1)) then
      write(*,csmsgs) trim(smeno),&
         &"found more than 2*csv%nc-1 of csi%edgeidx > 0 !"//csmsgerr ; stop
    endif
    if((csv%nc-1)**2+maxval(csi%edgeidx).gt.atoms%ncst(csi%atomType)) then
      write(*,csmsgs) trim(smeno),&
         &"found (csv%nc-1)^2+maxval(csi%edgeidx) > atoms%ncst(csi%atomType)!"//csmsgerr
      stop
    endif
    csv%nljc = count(csi%edgeidx.gt.0)
    if(.not.allocated(csv%lc)) allocate(csv%lc(csv%nljc))
    csv%lc = (/(edgel(csi%edgeidx(i)), i = 1,csv%nljc)/)
!!$    print*,csv%lc
    if(csi%verb.eq.1) write(*,csmsgsis)  trim(smeno),&
           &"nljc edge lines: ","csv%nljc = ",csv%nljc,"will be used"

! csi%lx
    if(csi%lx.lt.0) then
      write(*,csmsgs)  trim(smeno),"found csi%lx < 0 !"//csmsgerr ; stop
    endif
    if(csi%lx.gt.atoms%lmax(csi%atomType)) then
      write(*,csmsgs)  trim(smeno),&
           &"found csi%lx > atoms%lmax(csi%atomType)!"//csmsgerr ; stop
    endif
    if(csi%verb.eq.1) write(*,csmsgsis)  trim(smeno),&
           &"maximum l: ","csi%lx = ",csi%lx,"will be used"

! csi%ek0
    if(csi%ek0.le.0.d0) then
      write(*,csmsgs)  trim(smeno),"found csi%ek0 <= 0.0 !"//csmsgerr ; stop
    endif
    csi%ek0 = csi%ek0*1000.d0 ! conversion from keV to eV
    csv%gamma = 1.d0+csi%ek0/mec2
    csv%beta = sqrt(1.d0-1.d0/(csv%gamma**2))
    if(csi%verb.eq.1) then
      write(*,csmsgses)  trim(smeno),&
           &"kinetic energy of incoming electrons: ","csi%ek0 = ",csi%ek0,&
           &"eV will be used"
      write(*,csmsgses)  trim(smeno),&
           &"Lorentz factor (gamma): ","csv%gamma = ",csv%gamma,&
           &"will be used"
      write(*,csmsgses)  trim(smeno),&
           &"v/c factor (beta): ","csv%beta = ",csv%beta,&
           &"will be used"
    endif

! csi%emn csi%emx csi%ein -> csv%nex csv%egrid(0:csv%nex)
    if(csi%emn.gt.csi%emx) then
      write(*,csmsgs)  trim(smeno),"found csi%emn > csi%emx !"//csmsgerr ; stop
    endif
    if(csi%ein.le.0.d0) then
      write(*,csmsgs)  trim(smeno),"found csi%ein <= 0.d0 !"//csmsgerr ; stop
    endif
    if(((csi%emx-csi%emn)/csi%ein)-int((csi%emx-csi%emn)/csi%ein).ne.0) then
      write(*,csmsgs)  trim(smeno),&
           &"found non-integer (csi%emx-csi%emn)/csi%ein !"//csmsgerr ; stop
    endif
    csv%nex = int((csi%emx-csi%emn)/csi%ein)
    if(.not.allocated(csv%egrid)) allocate(csv%egrid(0:csv%nex))
    csv%egrid = (/(csi%emn+csi%ein*dble(i), i = 0,csv%nex)/)
    csv%nen = 0
!!$    do i = 0,csv%nex
!!$      if(csv%egrid(i).ge.0.d0) then
!!$        csv%nen = i
!!$        exit
!!$      endif
!!$    enddo

!!$    print*,csv%egrid
!!$    print*,csv%nen
    if(csi%verb.eq.1) write(*,csmsgsfs)  trim(smeno),&
           &"energy spectrum lower bound: ","csi%emn = ",csi%emn,&
           &"eV will be used"
    if(csi%verb.eq.1) write(*,csmsgsfs)  trim(smeno),&
           &"energy spectrum upper bound: ","csi%emx = ",csi%emx,&
           &"eV will be used"
    if(csi%verb.eq.1) write(*,csmsgsfs)  trim(smeno),&
           &"energy spectrum increment: ","csi%ein = ",csi%ein,"eV will be used"
    if(csi%verb.eq.1) write(*,csmsgsis)  trim(smeno),&
           &"no. of energy spectrum grid points: ","csv%nex = 0 : ",&
           &csv%nex,"will be used"
    if(csi%verb.eq.1) write(*,csmsgsis)  trim(smeno),&
           &"minimum index for which egrid >= 0: ","csv%nen = ",&
           &csv%nen,"will be used"

    if(.not.allocated(csv%eedge)) allocate(csv%eedge(csv%nljc))
    csv%eedge = 0.d0
    if(.not.allocated(csv%occ)) allocate(csv%occ(csv%nljc))
    csv%occ = 0.d0

    l_cs = .true.

    if(csi%verb.eq.1) write(*,*) ""

  end subroutine corespec_init
!
!===============================================================================
!===============================================================================
!
!  S U B R O U T I N E   I O U N I T
!
!-------------------------------------------------------------------------------
!
  subroutine iounit(unit)
!
! Assignes unused integer to I/O unit
!
    implicit none
!
    integer, intent(out) :: unit
!
    integer :: i
    logical :: fopen
!
    i = 9
    fopen = .true.
    do while(fopen)
      i = i+1
      inquire(unit=i,opened=fopen)
    enddo
!
    unit = i
!
  end subroutine iounit
!
!===============================================================================



end module m_corespec_io
