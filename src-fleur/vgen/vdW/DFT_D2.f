      MODULE DFT_D2
!
      IMPLICIT NONE
      PRIVATE
!
      REAL,DIMENSION(94,2) :: c6r
!
      PUBLIC :: driver_DFT_D2
!
      CONTAINS
!
      SUBROUTINE initialize_c6r
!
! the first index is "corresponds" to atomic number Z, i.e., 1  to H, 
!   2 to He etc.
! the second index corresponds to 
!   1 c6 coefficient
!   2 r_vdW of that element
!
      c6r(1:94,1)=0.0
      c6r(1:94,1)=1.0
!
      c6r(1,1) =  1.451020284	! H atom
      c6r(1,2) =  1.001
      c6r(2,1) =  0.829154448	! He atom
      c6r(2,2) =  1.012
      c6r(5,1) = 32.440667778	! B atom
      c6r(5,2) =  1.485
      c6r(6,1) = 18.137753550	! C atom
      c6r(6,2) =  1.452
      c6r(7,1) = 12.748249638	! N atom
      c6r(7,2) =  1.397
      c6r(8,1) =  7.25510142		! O atom
      c6r(8,2) =  1.342
      c6r(9,1) =  7.77332295		! F atom
      c6r(9,2) =  1.287
      c6r(11,1)= 59.180898726	! Na atom
      c6r(11,2)=  1.144
      c6r(16,1)= 57.729878442	! S atom 
      c6r(16,2)=  1.683
      c6r(17,1)= 52.547663142	! Cl atom
      c6r(17,2)=  1.639
      c6r(19,1)=111.93585048		! K atom
      c6r(19,2)=  1.485
      c6r(26,1)=111.93585048		! Fe atom
      c6r(26,2)=  1.562
      c6r(27,1)=111.93585048		! Co atom
      c6r(27,2)=  1.562
      c6r(29,1)=111.93585048		! Cu atom
      c6r(29,2)=  1.562
      c6r(35,1)=129.244449582	! Br atom
      c6r(35,2)=  1.749
      c6r(47,1)=255.690502902	! Ag atom
      c6r(47,2)=  1.639
      c6r(63,1)=  0.0		! Eu atom
      c6r(63,2)=  1.0
      c6r(77,1)=130.0		! Ir atom
      c6r(77,2)=  1.848
!
      END SUBROUTINE initialize_c6r
!
      SUBROUTINE calc_ene_DFT_D2(cart_init,z_init,cart_large,z_large)
      INTEGER, INTENT(IN),DIMENSION(:)   :: z_init,z_large
      REAL,INTENT(IN),DIMENSION(:,:)     :: cart_init,cart_large
!
      INTEGER  :: i1,i2
      REAL     :: dx,dy,dz
      REAL     :: c6,c6i,c6j,ri,rj
      REAL     :: dist,damp,vdW_energy
      REAL,ALLOCATABLE,DIMENSION(:,:) :: force_vdW
!
      vdW_energy=0.0
!
!      print*,'1: ',size(cart_init(3,:)),size(cart_large(3,:))
!      print*,z_init
!      print*,'2: ',size(z_init),size(z_large)
!
      if(size(cart_init(3,:)) <= 0) then
        STOP 'size(cart_init(3,:)) <= 0'
      else
        allocate(force_vdW(3,size(cart_init(3,:))))
	force_vdW(:,:)=0.0
      endif
!
      print'(20x,A)','Forces'
!
      do i1=1,size(cart_init(3,:))
	c6i=c6r(z_init(i1),1)
	ri =c6r(z_init(i1),2)
	do i2=1,size(cart_large(3,:))
	  c6j=c6r(z_large(i2),1)
	  rj =c6r(z_large(i2),2)
!
	  c6=0.75*sqrt(c6i*c6j)
!
	  dist=sqrt( (cart_init(1,i1)-cart_large(1,i2))**2 + 
     +               (cart_init(2,i1)-cart_large(2,i2))**2 + 
     +               (cart_init(3,i1)-cart_large(3,i2))**2 )
!
! avoid self-interaction
!
          if(dist > 0.0001) then
!
            damp = 1.0/(1.0+exp(-20.0*(dist/(ri+rj)-1.0)))
!
            vdW_energy = vdW_energy - 0.5*damp*c6/dist**6
!
            dx=cart_init(1,i1)-cart_large(1,i2)
	    dy=cart_init(2,i1)-cart_large(2,i2)
	    dz=cart_init(3,i1)-cart_large(3,i2)
!
            force_vdW(1,i1)=force_vdW(1,i1)+ c6*dx*                    
     +        ((damp**2)*20.0*exp(-20.0*(dist/(ri+rj)-1.0))/ 
     +        ((ri+rj)*dist**7) - 6.0*damp/dist**8)
	    force_vdW(2,i1)=force_vdW(2,i1)+ c6*dy*
     +        ((damp**2)*20.0*exp(-20.0*(dist/(ri+rj)-1.0))/ 
     +        ((ri+rj)*dist**7) - 6.0*damp/dist**8)
	    force_vdW(3,i1)=force_vdW(3,i1)+ c6*dz*                    
     +        ((damp**2)*20.0*exp(-20.0*(dist/(ri+rj)-1.0))/ 
     +        ((ri+rj)*dist**7) - 6.0*damp/dist**8)
!
	  endif
!
	enddo
!
	print'(A5,I3,3F14.8)','Atom ',i1,force_vdW(1:3,i1)
!
      enddo
!
      print'(A47,4x,F16.12)',
     &  'vdW energy calculated by DFT-D2 method (in eV):',vdW_energy
!
      print'(A47,4x,F16.12)',
     &  'vdW energy calculated by DFT-D2 method (in au):'
     &  ,vdW_energy/27.21138386
!
      if(allocated(force_vdW)) deallocate(force_vdW)
!
      END SUBROUTINE calc_ene_DFT_D2
!
      SUBROUTINE driver_DFT_D2(cart_init,z_init,cart_large,z_large)
      INTEGER, DIMENSION(:)   :: z_init,z_large
      REAL    ,DIMENSION(:,:) :: cart_init,cart_large
!
      call initialize_c6r
      call calc_ene_DFT_D2(cart_init,z_init,cart_large,z_large)
!
      END SUBROUTINE driver_DFT_D2
!
      END MODULE DFT_D2
!
