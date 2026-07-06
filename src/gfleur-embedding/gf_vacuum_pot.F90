!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_gf_vacuum_pot
    implicit none
    private
    real,allocatable,save    :: vz(:),q_z(:,:)
    complex,allocatable,save :: vxy(:,:),q_xy(:,:,:)
	integer,parameter      :: nmz=250,nmzxy=100
	real,parameter         :: delz=0.1


    public gf_vacuum_coulpot,gf_vacuum_totalpot,gf_vacuum_totalcharge
    contains

	subroutine gf_vacuum_makepot_g0(efield,delz,rht,vz)
	use m_qsf
	use m_constants,onlY:pimach
	implicit none
	real,intent(in):: efield
    real,intent(in):: delz
    real,intent(in):: rht(:)
    real,intent(out):: vz(:)

    real  ::  sig(size(rht))
    integer:: nmz,n
    nmz=size(rht)

    CALL qsf(delz,rht,sig,nmz,1)
	sig = sig(nmz) - sig

	write(6,*) "Vacuum-Charge integral:",sig(1)

    CALL qsf(delz,sig,vz,nmz,1)
    vz=-4.*pimach()*(vz(nmz)-vz)

    !Add an electric field if specified
    if (abs(Efield)>1E-9) then
          DO n=1,nmz
          	vz(n)=vz(n)+delz*n*Efield
          enddo
    endif


    write(6,*) "g_parallel=0 of Coulomb pot in Vacuum"
    do n=1,nmz
    	write(6,"(i5,5(1x,f10.5))") n,vz(n),rht(n),sig(n)
    enddo

    end subroutine

    subroutine gf_vacuum_makepot_g(delz,g,n2,rhtxy,vxy)
	use m_qsf
	USE m_intgr, ONLY : intgz1
    USE m_constants, ONLY : pimach

	implicit none
	real,intent(in)    :: delz,g
    integer,intent(in) :: n2
    complex,intent(in) :: rhtxy(:,:)
    complex,intent(inout):: vxy(:,:)


    integer :: ip,nmzxy,imz
    real    :: z,e_m,e_p
    complex :: alph0,betaz,alphaz,test
    real,dimension(size(rhtxy,1))::fra,fia,frb,fib
    real,dimension(size(rhtxy,1),2)::alpha,beta


    nmzxy=size(rhtxy,1)
    ip = nmzxy + 1
    z = 0
    DO  imz = 1,nmzxy
        e_m = exp_save( -g*z )
        e_p = exp_save( g*z )
        fra(ip-imz) = real(rhtxy(imz,n2-1))* e_m
        fia(ip-imz) = aimag(rhtxy(imz,n2-1))*e_m
        frb(imz) = real(rhtxy(imz,n2-1))* e_p
        fib(imz) = aimag(rhtxy(imz,n2-1))*e_p
        z = z + delz
    ENDDO

    CALL intgz1(fra,delz,nmzxy,alpha(:,1),tail=.true.)
    CALL intgz1(fia,delz,nmzxy,alpha(:,2),tail=.true.)
    CALL qsf(delz,frb,beta(:,1),nmzxy,1)
    CALL qsf(delz,fib,beta(:,2),nmzxy,1)

    alph0 = cmplx(alpha(nmzxy,1),alpha(nmzxy,2))
    z = 0
    DO imz = 1,nmzxy
         betaz = cmplx(beta(imz,1),beta(imz,2))
         alphaz = cmplx(alpha(ip-imz,1),alpha(ip-imz,2))
         e_m = exp_save( -g*z  )
         e_p = exp_save( g*z  )
         test = e_m*(alph0+betaz) + e_p*alphaz
         IF ( 2.0 * test == test ) test = cmplx(0.0,0.0)
         vxy(imz,n2-1) =  2.*pimach()/g * test
         z = z + delz
    enddo

      contains
!------------------------------------------------------------------
      PURE REAL FUNCTION exp_save(x)
      ! replace exp by a function that does not under/overflow dw09
      IMPLICIT NONE
      REAL   ,INTENT(IN)     :: x
      REAL, PARAMETER ::    maxexp = LOG(2.0)*MAXEXPONENT(2.0)
      REAL, PARAMETER ::    minexp = LOG(2.0)*MINEXPONENT(2.0)

      IF ( ABS(x)>minexp .AND. ABS(x)<maxexp ) THEN
         exp_SAVE = EXP(x)
      ELSE
         IF ( x > 0 ) THEN
            IF ( x > minexp ) THEN
               exp_save = EXP(maxexp)
            ELSE
               exp_save = EXP(minexp)
            ENDIF
         ELSE
            IF ( -x > minexp ) THEN
               exp_save = EXP(-maxexp)
            ELSE
               exp_save = EXP(-minexp)
            ENDIF
         ENDIF
      ENDIF
      END FUNCTION

	end subroutine

	real function gf_vacuum_totalcharge(jspins,cell,stars)
	use m_gf_types
	use m_gf_iodop
	use m_qsf
	integer,intent(in) :: jspins
	type(t_cell),intent(in) :: cell
	type(t_stars),intent(in) ::stars
	!Local variables
	real    :: qz(nmz,jspins)
	complex :: qxy(nmzxy,stars%nq2-1,jspins)
	real    :: q(nmz)

	call gf_iodop_readvacuum(GF_CDNFILE,qxy,qz)
	!Calculate total vaccum charge
	if (jspins>1) qz(:,1)=qz(:,1)+qz(:,2)
	CALL qsf(delz,qz(1,1),q,nmz,0)
    gf_vacuum_totalcharge = q(1)*cell%area
	end function


	subroutine gf_vacuum_coulpot(jspins,efield,stars,vb,l_fixed)
    use m_gf_iodop
    use m_gf_types
    implicit none
	integer,intent(in)       :: jspins
	real,intent(in)          :: efield
	type(t_stars),intent(in) :: stars
	complex,intent(out)      :: vb(:)
	logical,intent(in)       :: l_fixed

	!locals
	integer                :: n
	real                   :: evac
	real,allocatable       :: rht(:,:)
	complex,allocatable    :: rhtxy(:,:,:)
	logical                :: l_exist


    !read the charge density
    allocate(rht(nmz,jspins),rhtxy(nmzxy,stars%nq2-1,jspins))
    if (.not.allocated(q_z)) then
    	allocate(q_z(nmz,jspins),q_xy(nmzxy,stars%nq2-1,jspins))
    endif
	call gf_iodop_readvacuum(GF_CDNFILE,rhtxy,rht)
	q_z=rht
	q_xy=rhtxy


	if (jspins>1) then
		rht(:,1)=rht(:,1)+rht(:,2)
		rhtxy(:,:,1)=rhtxy(:,:,1)+rhtxy(:,:,2)
	endif


	!Generate potential
	allocate(vz(nmz),vxy(nmzxy,stars%nq2-1))
	call gf_vacuum_makepot_g0(efield,delz,rht(:,1),vz)

	DO n=2,stars%nq2
		call gf_vacuum_makepot_g(delz,stars%sk2(n),n,rhtxy(:,:,1),vxy)
	enddo

	!return the boundary values
	vb(1)=(vz(4)/3.-11./6.*vz(1)+3.*vz(2)-3./2.*vz(3))/delz

    write(6,"(a,3(f10.8,2x))") "Slope at vacuum-interstitial boundary:",real(vb(1))
    if (l_fixed) then
         vb(1)=vz(1)
    else
	  inquire(file="slope",exist=l_exist)
	  if (l_exist) then
		open(99,file="slope")
		read(99,*) vb(1)
		close(99)
	  endif
	endif

	do n=2,stars%nq2
		vb(n)=vxy(1,n-1)
	enddo
	deallocate(rht,rhtxy)
	end subroutine

	subroutine gf_vacuum_totalpot(vb,stars,xcpot,cell,jspins,l_noco,mpi_com)
	use m_gf_types
	use m_gf_iodop
    use m_vvacxc
    use m_vvacxcg
    use m_od_types
    implicit none
	real,intent(inout)      :: vb
	type(t_stars),intent(in):: stars
	type(t_xcpot),intent(in):: xcpot
	type(t_cell),intent(in) :: cell
	integer,intent(in)      :: jspins,mpi_com
	logical,intent(in)      :: l_noco


	real        :: evac
	COMPLEX     :: cdomvz(nmz,2)
    COMPLEX     :: cdomvxy(nmzxy,stars%nq2-1,2)
    REAL        :: excz(nmz,2)
    COMPLEX     :: excxy(nmzxy,stars%nq2-1,2)
	real,allocatable    :: rht(:,:),rht_xc(:,:,:),vz_xc(:,:,:)
	complex,allocatable :: rhtxy(:,:,:),rhtxy_xc(:,:,:,:),vxy_xc(:,:,:,:)
	integer             :: ichsmrg
	REAL                :: rhmn
	type(od_inp)        :: odi

	!shift the g==0 coulomb potential
	write(6,*) "Vacuum potential boundary:",vb,vz(1)
	vb=vb-vz(1)
	vz=vz+vb
	!read density
	if (l_noco) then
	   allocate(rht(nmz,4),rhtxy(nmzxy,stars%nq2-1,3))
	else
	   allocate(rht(nmz,jspins),rhtxy(nmzxy,stars%nq2-1,jspins))
	endif
	call gf_iodop_readvacuum(GF_CDNFILE,rhtxy,rht)
	allocate(rht_xc(nmz,2,jspins),rhtxy_xc(nmzxy,stars%nq2-1,2,jspins))
	rht_xc(:,1,:jspins)=rht(:,:jspins)     !need extra dimension for nvac
	rhtxy_xc(:,:,1,:jspins)=rhtxy(:,:,:jspins)
	if (l_noco) then
	   cdomvz(:,1)=cmplx(rht(:,3),rht(:,4))
	   cdomvxy(:,:,1)=rhtxy(:,:,3)
	else
	   cdomvz=0.0
	   cdomvxy=0.0
	endif


	!double coulomb potential if needed
	allocate(vz_xc(size(vz,1),2,jspins),vxy_xc(size(vxy,1),size(vxy,2),2,jspins))
	vz_xc(:,1,1)=vz(:)
	vxy_xc(:,:,1,1)=vxy(:,:)
	if (jspins>1) then
		vz_xc(:,1,2)=vz(:)
		vxy_xc(:,:,1,2)=vxy(:,:)
	endif
	if (allocated(vz)) deallocate(vz,vxy)  !deallocate saved module variable
	!Calculate xc-potential
	IF ( (xcpot%igrd == 0) .and. (xcpot%icorr /= -1) ) THEN
		call vvacxc(                                                                 &
        	   stars%mx1,stars%mx2,nmz,nmzxy,stars%nq2,jspins,9*stars%mx1*stars%mx2, &
           	   xcpot%icorr,.false.,xcpot%krla,nmzxy,jspins,stars%nq2,nmz,stars%nstr2,&
           	   rhtxy_xc,rht_xc,cdomvxy,cdomvz,l_noco,                                &
               stars%kimax2,stars%igfft2,stars%pgfft2,1,                             &
               vxy_xc,vz_xc,                                                         &
               excxy,excz)
    ELSE

	    rhmn=1e+10
	    odi%d1=.false.
	    write(*,*) jspins
        CALL vvacxcg(                                                                                              &
               stars%mx1,stars%mx2,nmz,nmzxy,stars%nq2,jspins,size(stars%igfft2,1)+1,9*stars%mx1*stars%mx2,1,0.0,  &
               delz,xcpot%icorr,.false.,nmzxy,jspins,stars%nq2,cell%bmat,nmz,stars%nstr2,                          &
               stars%pgft2x,stars%pgft2y,stars%pgft2xx,stars%pgft2yy,stars%pgft2xy,                                &
               xcpot%igrd,xcpot%ndvgrd,xcpot%chng,ichsmrg,                  &
               rhtxy_xc,rht_xc,cdomvxy,cdomvz,l_noco,.false.,(/0.0,0.0,0.0/),                                            & !lnoco=flase
               stars%kimax2,stars%igfft2,stars%pgfft2,odi,                                                         &
               vxy_xc,vz_xc,rhmn,                                                                                 &
               excxy,excz)
    ENDIF

	rht(:,:jspins)=vz_xc(:,1,:)
	rhtxy(:,:,:jspins)=vxy_xc(:,:,1,:)
	call gf_iodop_writevacuum(GF_POTFILE,rhtxy,rht,mpi_com,vb)
	deallocate(rht,rhtxy,vz_xc,vxy_xc,rht_xc,rhtxy_xc)
	end subroutine


end module m_gf_vacuum_pot
