      MODULE DFT_D3
!
      USE m_setr0ab, ONLY: setr0ab
!
      IMPLICIT NONE
!      PRIVATE
!
      INTEGER  :: max_elem,maxc
      REAL :: k1,k2,k3
      REAL :: rs6pbe,rs6revpbe
      REAL :: alp6
      REAL :: rcov(94)
!
!      PUBLIC :: driver_DFT_D3
!
      CONTAINS
!
!===================================================================
      SUBROUTINE initialize_DFT_D3(autoang,rcov_in)
!===================================================================
!
      REAL :: autoang
      REAL :: rcov_in(94)

      maxc    = 5
      max_elem=94
!
      k1=16.0 ! 16.0
      k2= 4.0/3.0
      k3=-4.0 ! -4.0
!
      rs6pbe   =1.217
      rs6revpbe=0.923
!
      autoang= 0.52917726
      alp6   =14.0
!
! covalent radii (taken from Pyykko and Atsumi, 
!   Chem. Eur. J. 15, 2009, 188-197)
! values for metals decreased by 10 %
!
      rcov_in = (/                                                     
     +  0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67, 
     +  1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54, 
     +  1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09, 
     +  1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39, 
     +  1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26, 
     +  1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57, 
     +  1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53, 
     +  1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32, 
     +  1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58, 
     +  1.52, 1.53, 1.54, 1.55 /)
!
      END SUBROUTINE initialize_DFT_D3
!
!===================================================================
      SUBROUTINE ncoord(natoms,rcov,iz,xyz,cn)
!===================================================================
      INTEGER, INTENT(IN)  :: natoms,iz(:)
      REAL,INTENT(IN)  :: xyz(:,:),rcov(94)
      REAL,INTENT(OUT) :: cn(:)
      INTEGER          :: i1,iat
      REAL             :: damp,dx,dy,dz
      REAL             :: r,rr,rco,xn
!
! calculate the fractional coordination number CN^{A}(r_{AB})
!   (Eq. 15 from J. Chem. Phys. 132, 154104 (2010))
!
      do i1=1,natoms
        xn=0.0
        do iat=1,natoms
          if(iat.ne.i1)then
            dx=xyz(1,iat)-xyz(1,i1)
            dy=xyz(2,iat)-xyz(2,i1)
            dz=xyz(3,iat)-xyz(3,i1)
            r=sqrt(dx*dx+dy*dy+dz*dz)
! covalent distance in Bohr
            rco=rcov(iz(i1))+rcov(iz(iat))
            rr=rco/r
! counting function exponential has a better long-range behavior than 
! MHGs inverse damping
            damp=1.0/(1.0+exp(-k1*(rr-1.0)))
            xn=xn+damp
          endif
        enddo
	cn(i1)=xn
!        write(*,*) i1,xn
      enddo
!
      END SUBROUTINE ncoord
!
!===================================================================
      SUBROUTINE getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,c6)
!===================================================================
      INTEGER, INTENT(IN)  :: maxc,max_elem,mxc(max_elem)
      INTEGER, INTENT(IN)  :: iat,jat
      REAL,INTENT(IN)  :: c6ab(max_elem,max_elem,maxc,maxc,3)
      REAL,INTENT(IN)  :: nci,ncj
      REAL,INTENT(OUT) :: c6
      INTEGER          :: i1,j1
      REAL             :: c6mem,cn1,cn2
      REAL             :: csum,rsum,r,tmp1
!
! calculate the C6^{AB} coefficient for atomic pair {A,B} 
!   (Eq. 16 from J. Chem. Phys. 132, 154104 (2010))
! mxc(iat)=N_A is the number of reference molecules for atom A
! mxc(jat)=N_B is the number of reference molecules for atom B
!
      c6mem=-1.d+99
      rsum=0.0
      csum=0.0
      c6  =0.0
!
      do i1=1,mxc(iat)
        do j1=1,mxc(jat)
          c6=c6ab(iat,jat,i1,j1,1)
          if(c6 .gt. 0.0)then
            c6mem=c6
            cn1=c6ab(iat,jat,i1,j1,2)
            cn2=c6ab(iat,jat,i1,j1,3)
!
! distance dependence through nci and ncj
!
            r=(cn1-nci)**2+(cn2-ncj)**2
            tmp1=exp(k3*r)
            rsum=rsum+tmp1     
            csum=csum+tmp1*c6
          endif
        enddo
      enddo
!
      if(rsum .gt. 0.0) then
        c6=csum/rsum
      else
        print*,'Debug c6: rsum lt 0.0'
	c6=c6mem
      endif
!
      END SUBROUTINE getc6
!
!===================================================================
      SUBROUTINE copyc6(maxc,max_elem,c6ab,maxci)
!===================================================================
      INTEGER, INTENT(IN)  :: maxc,max_elem
      INTEGER, INTENT(OUT) :: maxci(max_elem)
      REAL    ,INTENT(OUT) :: c6ab(max_elem,max_elem,maxc,maxc,3)
      INTEGER              :: nlines,nn,kk,iat,jat,iadr,jadr
      INCLUDE 'pars_21Jul2010.inc'
!
      c6ab=-1
      maxci=0
!
      kk=1
      do nn=1,nlines
        iat=int(pars(kk+1))
        jat=int(pars(kk+2))
!
        call limit(iat,jat,iadr,jadr)
        maxci(iat)=max(maxci(iat),iadr)
        maxci(jat)=max(maxci(jat),jadr)
!
        c6ab(iat,jat,iadr,jadr,1)=pars(kk)  
        c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
        c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)
!
        c6ab(jat,iat,jadr,iadr,1)=pars(kk) 
        c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
        c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
!
        kk=(nn*5)+1
!
      enddo
!
      END SUBROUTINE copyc6
!
!===================================================================
      SUBROUTINE limit(iat,jat,iadr,jadr)
!===================================================================
      INTEGER :: iat,jat,iadr,jadr
!
      iadr=1
      jadr=1
!
 10   if(iat .gt. 100) then
         iat=iat-100
         iadr=iadr+1
         goto 10
      endif
!
 20   if(jat .gt. 100) then
         jat=jat-100
         jadr=jadr+1
         goto 20
      endif
!
      END SUBROUTINE limit
!
! subroutines to evaluate vdW forces
!
!===================================================================
      SUBROUTINE force_ncoord(natoms,rcov,iz,xyz,dcn)
!===================================================================
      INTEGER, INTENT(IN)  :: natoms,iz(:)
      REAL,INTENT(IN)  :: xyz(:,:),rcov(94)
      REAL,INTENT(OUT) :: dcn(:,:)
      INTEGER          :: i1,iat
      REAL             :: tmp1,damp_x,damp_y,damp_z
      REAL             :: dx,dy,dz
      REAL             :: r,rr,rco
!
! note that rcov is scaled to au and by k2
! rcov()=k2*rcov/autoang
!
      do i1=1,natoms
        damp_x=0.0
	damp_y=0.0
	damp_z=0.0
        do iat=1,natoms
          if(iat.ne.i1)then
            dx=xyz(1,iat)-xyz(1,i1)
            dy=xyz(2,iat)-xyz(2,i1)
            dz=xyz(3,iat)-xyz(3,i1)
            r=sqrt(dx*dx+dy*dy+dz*dz)
! covalent distance in Bohr
            rco=rcov(iz(i1))+rcov(iz(iat))
            rr=rco/r
! counting function exponential has a better long-range behavior than 
! MHGs inverse damping
!            damp=1.0/(1.0+exp(-k1*(rr-1.0)))
!            xn=xn+damp
!
            tmp1=-k1*rco*exp(-k1*(rr-1.0)) / 
     +           ((1.0+exp(-k1*(rr-1.0)))**2)/(r**3)
	    damp_x=damp_x+dx*tmp1
	    damp_y=damp_y+dy*tmp1
	    damp_z=damp_z+dz*tmp1
!
          endif
        enddo
	dcn(1,i1)=damp_x
	dcn(2,i1)=damp_y
	dcn(3,i1)=damp_z
      enddo
!
      END SUBROUTINE force_ncoord
!
!===================================================================
      SUBROUTINE force_getc6(maxc,max_elem,c6ab,mxc,iat,jat, 
!===================================================================
     +                            nci,ncj,dnci,dncj,c6deriv)
      INTEGER, INTENT(IN)  :: maxc,max_elem,mxc(max_elem)
      INTEGER, INTENT(IN)  :: iat,jat
      REAL,INTENT(IN)  :: c6ab(max_elem,max_elem,maxc,maxc,3)
      REAL,INTENT(IN)  :: nci,ncj
      REAL,INTENT(IN)  :: dnci(3),dncj(3)
      REAL,INTENT(OUT) :: c6deriv(3)
      INTEGER          :: i1,j1
      REAL             :: c6mem,cn1,cn2,c6
      REAL             :: r,rsum,csum
      REAL             :: tmpx1,tmpy1,tmpz1,tmp1
      REAL             :: tmpx2,tmpy2,tmpz2
!
      c6mem=-1.d+99
      tmpx1=0.0; tmpx2=0.0
      tmpy1=0.0; tmpy2=0.0
      tmpz1=0.0; tmpz2=0.0
      rsum =0.0; csum =0.0
      c6        =0.0
      c6deriv(:)=0.0
!
! (d/dx)C_{AB}=(Z'W-ZW')/W^2
!  Z=csum
!  W=rsum
!
      do i1=1,mxc(iat)
        do j1=1,mxc(jat)
          c6=c6ab(iat,jat,i1,j1,1)
          if(c6 .gt. 0.0) then
            c6mem=c6
            cn1=c6ab(iat,jat,i1,j1,2)
            cn2=c6ab(iat,jat,i1,j1,3)
!
! distance dependence through nci and ncj
!
            r=(cn1-nci)**2+(cn2-ncj)**2
!
            tmp1 =exp(k3*r)
            rsum =rsum+tmp1   
            csum =csum+tmp1*c6
!
	    tmpx1=tmpx1+2.0*c6*k3*( dnci(1)*(cn1-nci) +      
     +                                 dncj(1)*(cn2-ncj) )*tmp1
	    tmpx2=tmpx2+2.0*   k3*( dnci(1)*(cn1-nci) +      
     +                                 dncj(1)*(cn2-ncj) )*tmp1
	    tmpy1=tmpy1+2.0*c6*k3*( dnci(2)*(cn1-nci) +      
     +                                 dncj(2)*(cn2-ncj) )*tmp1
	    tmpy2=tmpy2+2.0*   k3*( dnci(2)*(cn1-nci) +      
     +                                 dncj(2)*(cn2-ncj) )*tmp1
	    tmpz1=tmpz1+2.0*c6*k3*( dnci(3)*(cn1-nci) +      
     +                                 dncj(3)*(cn2-ncj) )*tmp1
	    tmpz2=tmpz2+2.0*   k3*( dnci(3)*(cn1-nci) +      
     +                                 dncj(3)*(cn2-ncj) )*tmp1
!
          endif
        enddo
      enddo
!      write(*,'(6e13.6)') rsum,csum,tmpx1,tmpy1,tmpx2,tmpy2
!
!#if 0
!!
!! the initial implementation
!!
!      c6deriv(1)=(tmpx1*rsum-csum*tmpx2)/rsum**2
!      c6deriv(2)=(tmpy1*rsum-csum*tmpy2)/rsum**2
!      c6deriv(3)=(tmpz1*rsum-csum*tmpz2)/rsum**2
!#else
      if(rsum .gt. 1.d-14) then
        c6deriv(1)=(tmpx1*rsum-csum*tmpx2)/rsum**2
        c6deriv(2)=(tmpy1*rsum-csum*tmpy2)/rsum**2
        c6deriv(3)=(tmpz1*rsum-csum*tmpz2)/rsum**2
      else
        c6deriv(1:3)=0.0
      endif
!#endif
!
      END SUBROUTINE force_getc6
!
!===================================================================
      SUBROUTINE calc_ene_DFT_D3(cart_init,z_init,cart_large,z_large,
     >                           autoang,
     <                           ener,force_vdW)
!===================================================================
      INTEGER,INTENT(IN),DIMENSION(:)   :: z_init,z_large
      REAL,INTENT(IN),   DIMENSION(:,:) :: cart_init,cart_large
      REAL,INTENT(OUT),  DIMENSION(:,:) :: force_vdW
      REAL,INTENT(IN) :: autoang
      REAL,INTENT(OUT):: ener
!
      INTEGER          :: iat,jat
      INTEGER          :: nr_tot_atoms,nr_atoms
      INTEGER          :: mxc(max_elem)
      REAL,ALLOCATABLE :: cn_init(:),cn_large(:)
      REAL,ALLOCATABLE :: dcn_init(:,:),dcn_large(:,:)
      REAL             :: c6,ftmp(3)
      REAL             :: c6deriv(3)
      REAL             :: c6ab(max_elem,max_elem,maxc,maxc,3)
      REAL             :: r0ab(max_elem,max_elem),ri,rj
      REAL             :: dist,damp,vdW_energy
      REAL             :: rr,tmp
      REAL             :: dx,dy,dz,r0ab_tmp
!
      vdW_energy=0.0
!
      nr_tot_atoms=size(cart_large(3,:))
      nr_atoms    =size(cart_init(3,:))
!
      allocate(   cn_init(nr_atoms),   cn_large(nr_tot_atoms))
      allocate(dcn_init(3,nr_atoms),dcn_large(3,nr_tot_atoms))
!
      force_vdW(1:3,1:nr_atoms)=0.0
!
! read the c6ab coefficients
!
      call copyc6(maxc,max_elem,c6ab,mxc)
!
! calculate coordination number CN          (now done in driver)
! scale rcov(:)         and convert to au
! scale cart_init(:,:)  and convert to au
! scale cart_large(:,:) and convert to au
!
      call ncoord(    nr_atoms,rcov, z_init, cart_init, cn_init)
      call ncoord(nr_tot_atoms,rcov,z_large,cart_large,cn_large)
!
! for vdW forces
!
      call force_ncoord(    nr_atoms,rcov, z_init, cart_init, dcn_init)
      call force_ncoord(nr_tot_atoms,rcov,z_large,cart_large,dcn_large)
!
! calculate r0ab
!
      call setr0ab(max_elem,autoang,r0ab)
!
!      print'(20x,A)','Forces'
!
      do iat=1,size(cart_init(3,:))
	do jat=1,size(cart_large(3,:))
!
	  dist=sqrt( (cart_init(1,iat)-cart_large(1,jat))**2 + 
     +               (cart_init(2,iat)-cart_large(2,jat))**2 + 
     +               (cart_init(3,iat)-cart_large(3,jat))**2 )
!
! avoid self-interaction
!
          if(dist > 0.0001) then
!
            rr =r0ab(z_large(jat),z_init(iat))/dist
	    tmp=rs6pbe*rr
	    damp = 1.0/(1.0+6.0*tmp**alp6)
! get C6
            call getc6(maxc,max_elem,c6ab,mxc,     
     +        z_init(iat), z_large(jat), 
     +        cn_init(iat),cn_large(jat),c6)
!
            vdW_energy = vdW_energy - 0.5*damp*c6/dist**6
!            write(*,*) cn_init(iat),cn_large(jat),c6
!
            call force_getc6(maxc,max_elem,c6ab,mxc,               
     +        z_init(iat),      z_large(jat),      
     +        cn_init(iat),     cn_large(jat),      
     +        dcn_init(1:3,iat),dcn_large(1:3,jat),c6deriv)
!
! define some temporary variables
!
            dx=cart_init(1,iat)-cart_large(1,jat)
	    dy=cart_init(2,iat)-cart_large(2,jat)
	    dz=cart_init(3,iat)-cart_large(3,jat)
	    r0ab_tmp=rs6pbe*r0ab(z_large(jat),z_init(iat))
!
!             write(*,'(5f15.10)') dx,dy,dz,c6deriv(:)
!            force_vdW(1,iat)=force_vdW(1,iat)                         + &
             ftmp(1) = 
     +         42.0*c6*(r0ab_tmp**14)*2.0*dx /                     
     +         (1.0+6.0*(r0ab_tmp**14)/dist**14)**2 / dist**22 + 
     +         c6deriv(1) / 
     +         (1.0+6.0*(r0ab_tmp**14)/dist**14)    / dist**6  - 
     +         6.0*c6*dx / 
     +         (1.0+6.0*(r0ab_tmp**14)/dist**14)    / dist**8
!
!	    force_vdW(2,iat)=force_vdW(2,iat)                         + &
             ftmp(2) = 
     +         42.0*c6*(r0ab_tmp**14)*2.0*dy /                     
     +         (1.0+6.0*(r0ab_tmp**14)/dist**14)**2 / dist**22 + 
     +         c6deriv(2) / 
     +         (1.0+6.0*(r0ab_tmp**14)/dist**14)    / dist**6  - 
     +         6.0*c6*dy / 
     +         (1.0+6.0*(r0ab_tmp**14)/dist**14)    / dist**8
!
!	    force_vdW(3,iat)=force_vdW(3,iat)                         + &
            ftmp(3) = 
     +         42.0*c6*(r0ab_tmp**14)*2.0*dz /                     
     +         (1.0+6.0*(r0ab_tmp**14)/dist**14)**2 / dist**22 + 
     +         c6deriv(3) / 
     +         (1.0+6.0*(r0ab_tmp**14)/dist**14)    / dist**6  - 
     +         6.0*c6*dz / 
     +         (1.0+6.0*(r0ab_tmp**14)/dist**14)    / dist**8

            force_vdW(:,iat) = force_vdW(:,iat) + ftmp(:)
!            write(*,'(2i5,3f15.10)') iat,jat,ftmp(:)
!
	  endif
!
	enddo
!
!	print'(A,I3,3F14.8)','Atom ',iat,force_vdW(1:3,iat)* 
!     +                              27.21138386/0.529177249
!
      enddo
!
      if(allocated(cn_init))   deallocate(cn_init)
      if(allocated(cn_large))  deallocate(cn_large)
      if(allocated(dcn_init))  deallocate(dcn_init)
      if(allocated(dcn_large)) deallocate(dcn_large)
!
      ener = vdW_energy*27.21138386
!      print'(A47,4x,F16.12)',
!     &   'vdW energy calculated by DFT-D3 method (in eV):',ener
!
!      print'(A47,4x,F16.12)',
!     &   'vdW energy calculated by DFT-D3 method (in au):',vdW_energy
!
      END SUBROUTINE calc_ene_DFT_D3
!
!===================================================================
      SUBROUTINE calc_ene_DFT_D3_debug(cart_init,z_init,
     +                                 cart_large,z_large,autoang)
!===================================================================

      INTEGER,INTENT(IN),DIMENSION(:)   :: z_init,z_large
      REAL   ,INTENT(IN),DIMENSION(:,:) :: cart_init,cart_large
      REAL   ,INTENT(IN) :: autoang

      INTEGER          :: iat,jat
      INTEGER          :: nr_tot_atoms,nr_atoms
      INTEGER          :: mxc(max_elem)
      REAL,ALLOCATABLE :: cn_init(:),cn_large(:)
      REAL             :: c6
      REAL             :: c6ab(max_elem,max_elem,maxc,maxc,3)
      REAL             :: r0ab(max_elem,max_elem),ri,rj
      REAL             :: dist,damp,vdW_energy
      REAL             :: rr,tmp
!
! calculate the vdW energy without repeating in space the initial unit cell
!
      vdW_energy=0.0
!
      nr_tot_atoms=size(cart_large(3,:))
      nr_atoms    =size(cart_init(3,:))
!
      allocate(cn_init(nr_tot_atoms),cn_large(nr_tot_atoms))
!
! read the c6ab coefficients
!
      call copyc6(maxc,max_elem,c6ab,mxc)
!
! calculate coordination number CN
! scale rcov(:)        and convert to au
! scale cart_init(:,:) and convert to au
!
!! now be careful if calc_ene_DFT_D3_debug is called after 
!!   calc_ene_DFT_D3 since in this case rcov(:) and cart_init(:,:) are
!!   ALREADY scaled!!!
!!
!#if USE_calc_ene_DFT_D3 == 0
!!
!      rcov(:)       =k2*rcov(:)/autoang
!      cart_init(:,:)=cart_init(:,:)/autoang
!#endif
!!
      call ncoord(nr_atoms,rcov,z_init,cart_init,cn_init)
!
! calculate r0ab
!
      call setr0ab(max_elem,autoang,r0ab)
!
! calculate ONLY for the initial unit cell
!
      do iat=1,size(cart_init(3,:))-1
	do jat=iat+1,size(cart_init(3,:))
!
	  dist=sqrt( (cart_init(1,iat)-cart_init(1,jat))**2 + 
     +               (cart_init(2,iat)-cart_init(2,jat))**2 + 
     +               (cart_init(3,iat)-cart_init(3,jat))**2 )
!
! avoid self-interaction
!
          if(dist > 0.0001) then
!
	    rr =r0ab(z_init(jat),z_init(iat))/dist
	    tmp=rs6pbe*rr
	    damp = 1.0/(1.0+6.0*tmp**alp6)
! get C6
            call getc6(maxc,max_elem,c6ab,mxc,    
     +                 z_init(iat), z_init(jat), 
     +                 cn_init(iat),cn_init(jat),c6)
!
            vdW_energy = vdW_energy - damp*c6/dist**6
	    !print*,'*',z_init(iat), z_init(jat),cn_init(iat),cn_init(jat)
	    !print*,dist,c6,damp,vdW_energy
!
	  endif
!
	enddo
      enddo
!
      if(allocated(cn_init))  deallocate(cn_init)
      if(allocated(cn_large)) deallocate(cn_large)
!
      print'(A)',
     + 'vdW energy calculated by DFT-D3 (debug) method (in eV):'
      print'(A)','                 (only one cell)             '
      print'(4x,F16.12)',vdW_energy*27.21138386
!
      print'(A)',
     + 'vdW energy calculated by DFT-D3 (debug) method (in au):'
      print'(A)','                 (only one cell)             '
      print'(4x,F16.12)',vdW_energy

      END SUBROUTINE calc_ene_DFT_D3_debug
!
!===================================================================
      SUBROUTINE driver_DFT_D3(cart_init_in, z_init,
     >                         cart_large_in,z_large,
     >                         l_in_au,
     <                         ener,force_vdW)
!===================================================================
      INTEGER, DIMENSION(:),INTENT(IN)   :: z_init,z_large
      REAL    ,DIMENSION(:,:),INTENT(IN) :: cart_init_in,cart_large_in
      REAL    ,INTENT(OUT),DIMENSION(:,:):: force_vdW
      LOGICAL,INTENT(IN)   :: l_in_au
      REAL,  INTENT(OUT)   :: ener
!
      INTEGER nr_atoms,nr_tot_atoms
      REAL autoang,rcov_in(94),rcov(94)
      REAL, ALLOCATABLE :: cart_init(:,:),cart_large(:,:)

      CALL initialize_DFT_D3 (autoang,rcov_in)
!
      nr_atoms    =size(cart_init_in(3,:))
      nr_tot_atoms=size(cart_large_in(3,:))
      ALLOCATE ( cart_init(3,nr_atoms),cart_large(3,nr_tot_atoms) )

      IF (.not.l_in_au) THEN ! convert to a.u.
        rcov(:)        =     k2*rcov_in(:)/autoang
        cart_init(:,:) = cart_init_in(:,:)/autoang
        cart_large(:,:)=cart_large_in(:,:)/autoang
      ELSE
        rcov(:)        =     k2*rcov_in(:)
        cart_init(:,:) = cart_init_in(:,:)
        cart_large(:,:)=cart_large_in(:,:)
      ENDIF
!
      CALL calc_ene_DFT_D3(cart_init,z_init,cart_large,z_large,autoang,
     <                     ener,force_vdW)
!      CALL calc_ene_DFT_D3_debug(cart_init,z_init,cart_large,z_large,autoang))
!
      END SUBROUTINE driver_DFT_D3
!
      END MODULE DFT_D3
