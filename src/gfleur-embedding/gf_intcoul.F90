!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_intcoul
      use m_juDFT
	  USE m_gf_types
	  USE m_constants
      IMPLICIT NONE
      PRIVATE
      REAL,PARAMETER:: CORRECT_CUTOFF=0.5 ! cutoff determine how many g_|| to fix
      TYPE t_chargelayer
      	 SEQUENCE
         REAL                :: c,dt,d
         INTEGER             :: layer,points
         COMPLEX,ALLOCATABLE :: pseudo_charge(:,:)
         COMPLEX,ALLOCATABLE :: ps_charge_boundary(:,:)
         REAL,ALLOCATABLE    :: atpos(:,:)
         REAL,ALLOCATABLE    :: rmt2(:)
         COMPLEX,ALLOCATABLE :: potential(:)
         COMPLEX,ALLOCATABLE :: project(:,:)
         COMPLEX,ALLOCATABLE :: dproject(:,:)
         INTEGER,ALLOCATABLE :: gz_max(:)
      END TYPE
      LOGICAL,SAVE :: l_fix  !is the vacuum fixed?
      public gf_makeintcoulpot
      CONTAINS
      SUBROUTINE gf_makeintcoulpot(jspins,layers,stars,mpi,gfinp             &
     &     ,potential,vac_pot,atoms,cell,sym)
      USE m_gf_iodop
      USE m_gf_plot
      IMPLICIT NONE
      integer,intent(in)           :: jspins
      TYPE(t_layers),INTENT(IN)    :: layers
      TYPE(t_stars),INTENT(IN)     :: stars(:)
      TYPE(t_atoms),INTENT(IN)     :: atoms(:)
      TYPE(t_cell),INTENT(IN)      :: cell(:)
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_mpi),INTENT(IN)       :: mpi
      TYPE(t_embinp),INTENT(IN)     :: gfinp
      TYPE(t_potden),INTENT(INOUT) :: potential(:)
      REAL   ,INTENT(OUT)          :: vac_pot

	  !locals
	  INTEGER:: n2,layer,n
	  TYPE(t_chargelayer),ALLOCATABLE::lay(:)
	  COMPLEX,ALLOCATABLE ::left(:),right(:),pdata(:,:,:)
	  REAL:: pos_l,pos_r
	  INTEGER,PARAMETER::printn2=1

      CALL timestart("gf_intcoul:setup")
	  ! read charges and generate setup
	  ALLOCATE(lay(layers%num_layers+2))
          lay%points=0
	  CALL priv_read_charges(lay,layers,stars(1),gfinp%l_surface)
      CALL priv_setup_layers(lay,stars,atoms,gfinp%l_surface)
      CALL priv_read_boundaries(layers%num_layers,jspins,stars(1),                 &
     &                    gfinp%l_surface,gfinp%efield,left,right,pos_l,pos_r)

	  !for printing
	  ALLOCATE(pdata(maxval(lay%points),4,size(lay)-2))

      !map the boundary charges
      DO n2=1,stars(1)%ng2
          DO layer=2,size(lay)-1
          	CALL priv_make_boundary_charge(lay(layer),lay(layer-1),lay(layer+1),n2)
          ENDDO
      ENDDO

      !loop over all layers and generate joined charge
      DO layer=2,size(lay)-1
       	CALL priv_join_charges(lay(layer),lay(layer-1),lay(layer+1),stars(layer-1),cell(layer-1))
      ENDDO

      DO n2=1,stars(1)%ng2
          !loop over all layers and generate potential
          DO layer=2,size(lay)-1
          	!CALL priv_make_boundary_charge(lay(layer),lay(layer-1),lay(layer+1),n2)
          	!IF (n2==printn2) pdata(:lay(layer)%points,1,layer-1)=lay(layer)%pseudo_charge(:,n2)
          	!IF (n2==printn2) pdata(:lay(layer)%points,2,layer-1)=lay(layer)%work+lay(layer)%pseudo_charge(:,n2)
          	CALL priv_make_potential(lay(layer),stars(lay(layer)%layer),n2)
          	IF (n2==printn2) pdata(:lay(layer)%points,3,layer-1)=lay(layer)%potential
          ENDDO
          !correct for external boundary conditions , match only small g-vectors
          IF (stars(1)%sk2(n2)*minval(lay(2:size(lay)-1)%dt/lay(2:size(lay)-1)%points)<CORRECT_CUTOFF) THEN
              CALL priv_correct_boundary(lay,stars(1),left,right,n2,gfinp%l_surface,vac_pot)
          ENDIF
          !IF (n2<60)
          !CALL priv_correct_boundary(lay,stars(1),left,right,n2,gfinp%l_surface,vac_pot)
          DO layer=2,size(lay)-1
          	 IF (n2==printn2) pdata(:lay(layer)%points,4,layer-1)=lay(layer)%potential
          ENDDO
          IF (n2==printn2) CALL priv_print_data(lay,pdata,n2)
          CALL priv_collect_potentials(lay,potential,stars,n2)

      ENDDO
      DO layer=1,layers%num_layers
      	CALL gf_wrtcoul(GF_POTFILE,layer,stars(layer),potential(layer  &
     &        )%pw(:,1))
        DO n=-stars(layer)%mx3,stars(layer)%mx3
           write(44,*) n,potential(layer)%pw(max(1,stars(layer)%ig(0,0,n)),1)
        enddo
        CALL gf_plot(layer,stars(layer),cell(layer),atoms(layer),sym,1,                     &
        potential(layer)%pw(:,:),GF_PLOT_HARTREE)


      ENDDO



      CALL timestop("gf_intcoul:setup")

      END SUBROUTINE

     SUBROUTINE priv_setup_layers(lay,stars,atoms,l_surface)
     IMPLICIT NONE
     TYPE(t_chargelayer),INTENT(INOUT)::lay(:)
     TYPE(t_stars),INTENT(IN)         ::stars(:)
     TYPE(t_atoms),INTENT(IN)         :: atoms(:)
     LOGICAL,INTENT(IN)               :: l_surface

	 !locals
	 INTEGER:: l,n
	 COMPLEX:: g




     DO l=1,size(lay)
     	lay(l)%layer=l-1
     	IF (l==1) THEN
     	      lay(l)%layer=0
     	      CYCLE
     	 ELSEIF(l==size(lay)) THEN
     	       lay(l)%layer=0
     	       CYCLE
     	 ENDIF
     	 !Set the atom info
     	 allocate(lay(l)%atpos(3,atoms(l-1)%nat))
     	 allocate(lay(l)%rmt2(atoms(l-1)%nat))
     	 lay(l)%atpos=atoms(l-1)%pos
     	 lay(l)%rmt2=atoms(l-1)%rmt**2

     	 lay(l)%points=stars(l-1)%mx3*2+1
     	 allocate(lay(l)%gz_max(stars(l-1)%ng2))
     	 DO n=1,stars(l-1)%ng2
     	    lay(l)%gz_max(n)=maxval(stars(l-1)%kv3(3,:),mask=stars(l-1)%ig2==n)
     	 ENDDO
     	 ALLOCATE(lay(l)%project(lay(l)%points,2))
     	 ALLOCATE(lay(l)%dproject(lay(l)%points,2))
     	 ALLOCATE(lay(l)%potential(lay(l)%points))
     	 ALLOCATE(lay(l)%ps_charge_boundary(lay(l)%points,stars(1)%ng2))


     	 !construct the projectors
     	 DO n = 1,lay(l)%points
           g = n-1
           IF (n>lay(l)%points/2+1) g=n-lay(l)%points-1
           g=cmplx(0.0,1.0)*g*2.*pi_const/lay(l)%dt
           lay(l)%project(n,1)=exp(g*(-lay(l)%c/2))
           lay(l)%project(n,2)=exp(g*(lay(l)%c/2))
           lay(l)%dproject(n,1)=g*exp(g*(-lay(l)%c/2))
           lay(l)%dproject(n,2)=g*exp(g*(lay(l)%c/2))
         ENDDO
      ENDDO
      END SUBROUTINE


     SUBROUTINE priv_correct_boundary(lay,stars,left,right,n2,l_surface,vac_pot)
     USE m_gf_math
     USE m_gf_fft_singleton
     use m_juDFT
     IMPLICIT NONE
     TYPE(t_chargelayer),INTENT(INOUT)::lay(:)
     TYPE(t_stars),INTENT(IN)         ::stars
     COMPLEX,INTENT(IN)               ::right(:),left(:)
     INTEGER,INTENT(IN)               ::n2
     LOGICAL,INTENT(IN)               :: l_surface
     real,intent(out)                 :: vac_pot
	 !locals
	 INTEGER :: n,i,row,column
	 COMPLEX :: c(2,2),dc(2,2),v(2),dv(2)
	 COMPLEX,ALLOCATABLE:: corr(:,:,:),rs(:),a(:),mat(:,:)
     logical pr
     pr=n2==1
     pr=.true.
	 n=2*(size(lay)-2)
	ALLOCATE(corr(maxval(lay%points),2,n/2))
	ALLOCATE(mat(n,n),rs(n),a(n))
	corr=0.0
	mat=0.0
	rs=0.0
	IF (pr) WRITE(oUnit,*) "Uncorrected Boundary values n2=",n2
	IF (pr) WRITE(oUnit,*) "Left Bulk:",left(n2)

	DO n=2,size(lay)-1
		v(1)=dot_product(conjg(lay(n)%project(:,1)),lay(n)%potential)
	    v(2)=dot_product(conjg(lay(n)%project(:,2)),lay(n)%potential)
	    dv(1)=dot_product(conjg(lay(n)%dproject(:,1)),lay(n)%potential)
	    dv(2)=dot_product(conjg(lay(n)%dproject(:,2)),lay(n)%potential)
	    IF (pr) WRITE(oUnit,"(a,i4,4f10.4)") "L",n,v(1),dv(1)
	    IF (pr) WRITE(oUnit,"(a,i4,4f10.4)") "R",n,v(2),dv(2)
	 ENDDO
	 IF (pr) THEN
	   IF (l_surface) THEN
	      WRITE(oUnit,*) "Surface calculation:",right(n2)
	   ELSE
	 	  WRITE(oUnit,*) "Right Bulk:",right(n2)
	   ENDIF
	ENDIF
	!Construct the linear equation
	row=1
	column=1
	rs(1)=left(n2)
	rs(size(rs))=right(n2)
	DO n=2,size(lay)-1
	     CALL priv_correction_potentials(lay(n),stars,n2,corr(:lay(n)%points,:,n-1),c,dc)
	     v(1)=dot_product(conjg(lay(n)%project(:,1)),lay(n)%potential)
	     v(2)=dot_product(conjg(lay(n)%project(:,2)),lay(n)%potential)
	     dv(1)=dot_product(conjg(lay(n)%dproject(:,1)),lay(n)%potential)
	     dv(2)=dot_product(conjg(lay(n)%dproject(:,2)),lay(n)%potential)
	     IF (n==2) THEN
	     	mat(row,column)=c(1,1)
	     	mat(row,column+1)=c(2,1)
	     	rs(row)=rs(row)-v(1)
	     	if (n==size(lay)-1) then
	     	    !we only have a single layer
	     	    if (l_surface.and..not.l_fix) then
	     	    	mat(row+2,column)=dc(1,2)
	        		mat(row+2,column+1)=dc(2,2)
	        		rs(row+2)=rs(row+2)-dv(2)
	        		cycle
	        	endif
	     	    mat(row+1,column)=c(1,2)
	        	mat(row+1,column+1)=c(2,2)
	        	rs(row+1)=rs(row+1)-v(2)
	        	cycle
	        endif
	     	mat(row+1,column)=c(1,2)
	     	mat(row+1,column+1)=c(2,2)
	     	rs(row+1)=rs(row+1)-v(2)
	     	mat(row+2,column)=dc(1,2)
	     	mat(row+2,column+1)=dc(2,2)
	     	rs(row+2)=rs(row+2)-dv(2)
	     	row=row+1
	     	column=column+2
	     	CYCLE
	     ENDIF
	     IF (n==size(lay)-1) THEN
	        mat(row,column)=-c(1,1)
	        mat(row,column+1)=-c(2,1)
	        rs(row)=rs(row)+v(1)
	        mat(row+1,column)=-dc(1,1)
	        mat(row+1,column+1)=-dc(2,1)
	        rs(row+1)=rs(row+1)+dv(1)
	        IF (l_surface.and..not.l_fix) THEN
	        	mat(row+2,column)=dc(1,2)
	        	mat(row+2,column+1)=dc(2,2)
	        	rs(row+2)=rs(row+2)-dv(2)
	        ELSE
	        	mat(row+2,column)=c(1,2)
	        	mat(row+2,column+1)=c(2,2)
	        	rs(row+2)=rs(row+2)-v(2)
	        ENDIF
		    CYCLE
		 ENDIF
        mat(row,column)=-c(1,1)
        mat(row,column+1)=-c(2,1)
        rs(row)=rs(row)+v(1)
        mat(row+1,column)=-dc(1,1)
        mat(row+1,column+1)=-dc(2,1)
        rs(row+1)=rs(row+1)+dv(1)
     	mat(row+2,column)=c(1,2)
     	mat(row+2,column+1)=c(2,2)
     	rs(row+2)=rs(row+2)-v(2)
     	mat(row+3,column)=dc(1,2)
     	mat(row+3,column+1)=dc(2,2)
     	rs(row+3)=rs(row+3)-dv(2)
	    row=row+2
	    column=column+2
 	ENDDO
 	IF (pr) THEN
 		WRITE(oUnit,*) "Matching Matrix"
		DO row=1,size(mat,1)
	    	WRITE(oUnit,"(999f8.4)") mat(row,:)
		ENDDO
		WRITE(oUnit,*) "Right-hand side"
		WRITE(oUnit,"(999f8.4)") rs(:)
    ENDIF
	a=lin_equation(mat,rs)
    IF (pr) WRITE(oUnit,*) "Solution"
	IF (pr) WRITE(oUnit,"(999f8.4)") a(:)
	!now correct the potential in all layers
	i=1
	DO n=2,size(lay)-1
		lay(n)%potential=lay(n)%potential+a(i)*corr(:lay(n)%points,1,n-1)+a(i+1)*corr(:lay(n)%points,2,n-1)
		!print *, "Onlycorrection"
		!lay(n)%potential=a(i)*corr(:lay(n)%points,1,n-1)+a(i+1)*corr(:lay(n)%points,2,n-1)
		!if (n2.ne.40) lay(n)%potential=0
		i=i+2
	ENDDO
	!print the boundary values
	!IF (n2/=1) RETURN
	WRITE(oUnit,*) "Boundary values n2=",n2
	WRITE(oUnit,*) "Left Bulk:",left(n2)
	DO n=2,size(lay)-1
	     v(1)=dot_product(conjg(lay(n)%project(:,1)),lay(n)%potential)
	     v(2)=dot_product(conjg(lay(n)%project(:,2)),lay(n)%potential)
	     dv(1)=dot_product(conjg(lay(n)%dproject(:,1)),lay(n)%potential)
	     dv(2)=dot_product(conjg(lay(n)%dproject(:,2)),lay(n)%potential)
	    WRITE(oUnit,"(a,i4,4f10.4)") "L",n,v(1),dv(1)
	    WRITE(oUnit,"(a,i4,4f10.4)") "R",n,v(2),dv(2)
	 ENDDO
	 IF (l_surface) THEN
	   WRITE(oUnit,*) "Surface Calculation:",right(n2)
	 ELSE
	 	WRITE(oUnit,*) "Right Bulk:",right(n2)
     ENDIF
     if (pr) vac_pot=v(2)
     END SUBROUTINE


	 FUNCTION priv_1d_solution(kappa,n2,c,d,v0,vc,pos)RESULT(func)
	 IMPLICIT NONE
	 INTEGER,INTENT(IN) :: n2
	 REAL   ,INTENT(IN) :: kappa,c,d
	 COMPLEX,INTENT(IN) :: v0,vc
	 REAL   ,INTENT(IN) :: pos(:)
	 REAL               :: func(size(pos))
     !locals
	 REAL    :: maxpos,z,a1,a2,a3,a4,d1,d2,a,b,dd
	 INTEGER :: n

	 maxpos=maxval(pos)
	 maxpos=maxpos-minval(pos)
	 maxpos=maxpos
	 dd=maxpos-d

	 !determine coefficients for smooth function in aux-region
	 IF (n2==1) THEN
	 	d1=(vc-v0)/c
	 	d2=d1
	 ELSE
	    d1=kappa*v0
	    d2=kappa*vc
	 ENDIF
	 a1 = (dd*d1+d2*dd+2*vc-2*v0)/dd**3
	 a2 = -(dd*d1+2*d2*dd+3*vc-3*v0)/dd**2
	 a3 = d2
	 a4 = vc

	 !make linear interpolation in aux-region
	 a2=v0
	 a4=vc
	 a1=(v0-vc)/(maxpos-d)

	 DO n=1,size(pos)
	    IF (pos(n)>(c-d)/2.AND.pos(n)<c+(d-c)/2) THEN
	    !if (2>1) then
	       IF(n2==1) THEN !linear function between 0 and c
	             func(n)=pos(n)*vc/c+(c-pos(n))*v0/c
	       ELSE
	             a = -(-vc+v0*exp(kappa*c))/(exp(-kappa*c)-exp(kappa*c))
	             b = (exp(-kappa*c)*v0-vc)/(exp(-kappa*c)-exp(kappa*c))
	             func(n)=a*exp(-kappa*pos(n))+b*exp(kappa*pos(n))
	       ENDIF
	     ELSE
	         !in between we have to interpolate nicely
	          IF (pos(n)<0) THEN
	              func(n)=a1*(pos(n)-(c-d)/2)+a2
	          ELSE
	              func(n)=a1*(pos(n)-(c+(d-c)/2))+a4
	          ENDIF
!             if (pos(n)<0) then
!                  z=maxpos+pos(n)-c
!             else
!                  z=pos(n)-c
!             endif
!             func(n)=a1*z**3+a2*z**2+a3*z+a4
         ENDIF
     ENDDO
     END FUNCTION

     SUBROUTINE priv_collect_potentials(lay,potential,stars,n2)
     IMPLICIT NONE
     TYPE(t_chargelayer),INTENT(IN)  ::lay(:)
     TYPE(t_potden),INTENT(INOUT) ::potential(:)
     TYPE(t_stars),INTENT(IN)        ::stars(:)
     INTEGER,INTENT(IN)              ::n2

     INTEGER:: index,l,layer,n,g


     DO l=2,size(lay)-1
	    layer=lay(l)%layer
	    DO n=1,lay(l)%points
	        g=n-1
	        IF (n>lay(l)%points/2+1) g=n-lay(l)%points-1
	        IF (abs(g)>stars(layer)%mx3) CYCLE
	        index=stars(layer)%ig(stars(layer)%kv2(1,n2),stars(layer)%kv2(2,n2),g)
	        IF (index==0) CYCLE
	        potential(layer)%pw(index,1)=potential(layer)%pw(index,1)+lay(l)%potential(n)
	    ENDDO
	 ENDDO
	 END SUBROUTINE

      SUBROUTINE priv_READ_boundaries(num_layers,jspins,stars,l_surface,efield,left,right,pos_l,pos_r)
!-----------------------------------------------
!
!           (last modified:08-07-09) D. Wortmann
!-----------------------------------------------
      USE m_hdf_tools
      USE hdf5
      USE m_gf_io2dmat
      USE m_constants
      USE m_gf_energies
      USE m_gf_vacuum_pot
      IMPLICIT NONE
      !Arguments
      INTEGER,INTENT(IN)                :: num_layers,jspins
      real,intent(in)                   :: efield
      type(t_stars),intent(in)          :: stars
      COMPLEX,ALLOCATABLE,INTENT(OUT)   :: left(:),right(:)
      REAL,INTENT(OUT)                  :: pos_l,pos_r
      LOGICAL,INTENT(IN)                :: l_surface
      !<-- Locals
      INTEGER             :: n,n2,hdferr
      REAL                :: bias
      INTEGER(HID_T)      :: fid,gid
      REAL   ,ALLOCATABLE :: b(:,:)

      !left side
      fid=gf_io2dmatFID(IO2D_EMB,1,1,IO2D_READ)
      CALL io_gopen(fid,"Boundary",gid,hdferr)
      CALL io_READ_att(gid,"n2",n2)
      ALLOCATE(left(n2))
      ALLOCATE(b(n2,2))
      CALL io_READ_att(gid,"bulk_real",b(:,1))
      CALL io_READ_att(gid,"bulk_imag",b(:,2))
      left = CMPLX(b(:,1),b(:,2))
      DEALLOCATE(b)
      CALL io_READ_att(gid,"pos",pos_l)
      CALL io_gclose(gid,hdferr)


      ! right side

      IF (.NOT.l_surface) THEN
         fid = gf_io2dmatFID(IO2D_EMB,num_layers,2,IO2D_READ)
         CALL io_gopen(fid,"Boundary",gid,hdferr)
         CALL io_READ_att(gid,"n2",n2)
         ALLOCATE(right(n2))
         ALLOCATE(b(n2,2))
         CALL io_READ_att(gid,"bulk_real",b(:,1))
         CALL io_READ_att(gid,"bulk_imag",b(:,2))
         right = CMPLX(b(:,1),b(:,2))
         DEALLOCATE(b)
	     CALL io_READ_att(gid,"pos",pos_r)
         CALL io_gclose(gid,hdferr)
         n = gf_bias_layer(bias)
         right(1) = right(1)+bias
      ELSE
         allocate(right(n2))
         inquire(file="fix.vacuum",exist=l_fix)
         call gf_vacuum_coulpot(jspins,efield,stars,right,l_fix)

         if (l_fix) then
              open(99,file="fix.vacuum")
              read(99,*) bias
              close(99)
              right(1)=right(1)+bias
          endif

      ENDIF
      END SUBROUTINE




      SUBROUTINE priv_read_charges(lay,layers,stars,l_surface)
!-----------------------------------------------
!
!           (last modified:08-07-10) D. Wortmann
!-----------------------------------------------
      USE m_gf_types
      USE m_hdf_tools
      USE hdf5
      USE m_gf_io2dmat
      USE m_constants
      USE m_gf_iodop
      use m_juDFT
      IMPLICIT NONE
      !<--Arguments
      TYPE(t_chargelayer ),INTENT(INOUT) :: lay(:)
      TYPE(t_layers),INTENT(IN)          :: layers
      TYPE(t_stars),INTENT(IN)           :: stars
      LOGICAL,INTENT(IN)                 :: l_surface
      !<-- Locals
      INTEGER(HID_T)             :: fid,gid,did
      INTEGER                    :: DIM(3)
      INTEGER                    :: l,hdferr

      REAL:: D(2),ph !phase can only be real!?
      integer::n2

      !<-- read charge for left and right layer
      !left
      lay(1)%layer = 0
      fid=gf_io2dmatFID(IO2D_EMB,1,1,IO2D_READ)
      CALL io_gopen(fid,"Boundary",gid,hdferr)
      IF (hdferr /= 0) CALL juDFT_error("Could not read left Boundary charge")
      CALL io_dopen(gid,"pseudo",did,hdferr)
      CALL io_datadim(did,dim)
      ALLOCATE(lay(1)%pseudo_charge(DIM(2),stars%ng2))
      lay(1)%points=dim(2)
      lay(1)%pseudo_charge = 0.0
      IF (DIM(3)>stars%ng2) DIM(3) = stars%ng2
      IF (DIM(3) /= stars%ng2) WRITE(*,*) "WARNING         ! Insuficient coulomb charge in left embpot"
      dim(1)=1
      CALL io_READ(did,(/-1,1,1/),dim,"lay",lay(1)%pseudo_charge(:,:DIM(3)))
      CALL io_dclose(did,hdferr)
      CALL io_READ_att(gid,"c",lay(1)%c)
      CALL io_READ_att(gid,"dt",lay(1)%dt)
      !read atom info
      CALL io_dopen(gid,"rmt",did,hdferr)
      CALL io_datadim(did,dim)
      ALLOCATE(lay(1)%rmt2(dim(1)),lay(1)%atpos(3,dim(1)))
      CALL io_read(did,(/1/),dim,"lay",lay(1)%rmt2)
      lay(1)%rmt2=lay(1)%rmt2**2
      CALL io_dclose(did,hdferr)
      CALL io_dopen(gid,"atpos",did,hdferr)
      CALL io_read(did,(/1,1/),(/3,dim(1)/),"lay",lay(1)%atpos)
      CALL io_dclose(did,hdferr)

      CALL io_gclose(gid,hdferr)

      !shift the charge if needed
      !print *, "Manual shifting of psq"
      !d=(/0.5,0.5/)
      !DO n2=1,stars%ng2
      !    ph=exp(cmplx(0,2.*pi_const*dot_product(d,stars%kv2(:,n2))))
      !    lay(1)%pseudo_charge(:,n2)=lay(1)%pseudo_charge(:,n2)*ph
      !ENDDO
      !Warning
      !write(*,*) "WARNING, left pseudo-charge set to zero"
      !lay(1)%pseudo_charge=0.0
      !right
      IF (.NOT.l_surface) THEN
         lay(SIZE(lay))%layer = 0
         fid = gf_io2dmatFID(IO2D_EMB,layers%num_layers,2,IO2D_READ)
         CALL io_gopen(fid,"Boundary",gid,hdferr)
         IF (hdferr /= 0) CALL juDFT_error("Could not read right Boundary charge")
         CALL io_dopen(gid,"pseudo",did,hdferr)
         CALL io_datadim(did,dim)
         ALLOCATE(lay(SIZE(lay))%pseudo_charge(DIM(2),stars%ng2))
         lay(size(lay))%points=dim(2)
         lay(SIZE(lay))%pseudo_charge = 0.0
         IF (DIM(3)>stars%ng2) DIM(3)    = stars%ng2
         IF (DIM(3) /= stars%ng2) WRITE(*,*)  "WARNING      ! Insuficient coulomb charge in right embpot"
         DIM(1) = 1
         CALL io_READ(did,(/-1,1,1/),dim,"lay",lay(SIZE(lay))%pseudo_charge(:,:DIM(3)))
         CALL io_dclose(did,hdferr)
         CALL io_READ_att(gid,"c",lay(SIZE(lay))%c)
         CALL io_READ_att(gid,"dt",lay(SIZE(lay))%dt)
         !read atom info
         CALL io_dopen(gid,"rmt",did,hdferr)
         CALL io_datadim(did,dim)
         ALLOCATE(lay(SIZE(lay))%rmt2(dim(1)),lay(SIZE(lay))%atpos(3,dim(1)))
         CALL io_read(did,(/1/),dim,"lay",lay(SIZE(lay))%rmt2)
         lay(SIZE(lay))%rmt2=lay(SIZE(lay))%rmt2**2
         CALL io_dclose(did,hdferr)
         CALL io_dopen(gid,"atpos",did,hdferr)
         CALL io_read(did,(/1,1/),(/3,dim(1)/),"lay",lay(SIZE(lay))%atpos)
         CALL io_dclose(did,hdferr)
         CALL io_gclose(gid,hdferr)
      ELSE
         ALLOCATE(lay(SIZE(lay))%pseudo_charge(3,stars%ng2))
         lay(SIZE(lay))%pseudo_charge=0.0
         lay(SIZE(lay))%c=lay(1)%c
         lay(SIZE(lay))%dt=lay(1)%dt
         ALLOCATE(lay(SIZE(lay))%rmt2(0))
      ENDIF
     !shift the charge if needed
     ! print *, "Manual shifting of psq"
     ! d=(/0.5,0.5/)
     ! DO n2=1,stars%ng2
     !     ph=exp(cmplx(0,2.*pi_const*dot_product(d,stars%kv2(:,n2))))
     !     lay(size(lay))%pseudo_charge(:,n2)=lay(size(lay))%pseudo_charge(:,n2)*ph
     ! ENDDO


      !<-- loop over all layers
      DO l = 1,layers%num_layers
         CALL gf_iodop_readpseudo(l,lay(l+1)%pseudo_charge)
         lay(l+1)%dt = layers%dt(l)
         lay(l+1)%c  = layers%c(l)
         lay(l+1)%d  = layers%d(l)
         lay(l+1)%layer = l
      ENDDO


      END SUBROUTINE

      SUBROUTINE priv_make_boundary_charge(layer,leftlayer,rightlayer,n2)
      USE m_gf_fft_singleton
      IMPLICIT NONE
      TYPE(t_chargelayer),INTENT(INOUT)::layer
      TYPE(t_chargelayer),INTENT(IN)   ::leftlayer
      TYPE(t_chargelayer),INTENT(IN)   ::rightlayer
	  INTEGER,INTENT(IN)               ::n2

      REAL    :: pos(layer%points),pos_left(layer%points),pos_right(layer%points)
      COMPLEX :: charge_left(layer%points),charge_right(layer%points)
      INTEGER :: n,i,g



	  !positions in current layer
	  pos=(/ (layer%dt/(layer%points+1)*n,n=0,layer%points/2),                 &
       (layer%dt/(layer%points+1)*(n-layer%points),n=layer%points/2+1,layer%points-1) /)
	  pos=(/ (layer%dt/(layer%points)*n,n=0,layer%points-1)/)
	  pos=pos- (/(0.0,n=0,layer%points/2),(layer%dt,n=layer%points/2+1,layer%points-1) /)
	  !positions in ajacent layer
	  pos_left=pos+(layer%c+leftlayer%c)/2.0
	  pos_right=pos-(layer%c+rightlayer%c)/2.0

	  !evaluate the left charge on grid
	  charge_left=0.0
	  DO n=1,layer%points
	      !IF (pos_left(n)<0.0.OR.pos_left(n)>leftlayer%dt/2) CYCLE
	      IF (pos_left(n)>leftlayer%dt/2) CYCLE

	      DO i=1,leftlayer%points
	        g=i-1
	      	IF (i>leftlayer%points/2+1) g=i-leftlayer%points-1
	      	charge_left(n)=charge_left(n)+leftlayer%pseudo_charge(i,n2)*exp(cmplx(0.0,pos_left(n)*g*2*pi_const/leftlayer%dt))
	  	  ENDDO
	  ENDDO

	  !evaluate the left charge on grid
	  charge_right=0.0
	  DO n=1,layer%points
	      !IF (pos_right(n)<-rightlayer%dt/2.0.OR.pos_right(n)>0.0) CYCLE
	      IF (pos_right(n)<-rightlayer%dt/2.0) CYCLE
	      DO i=1,rightlayer%points
	        g=i-1
	      	IF (i>rightlayer%points/2+1) g=i-rightlayer%points-1
	      	charge_right(n)=charge_right(n)+rightlayer%pseudo_charge(i,n2)*exp(cmplx(0.0,pos_right(n)*g*2*pi_const/rightlayer%dt))
	  	  ENDDO
	  ENDDO

      !charge_left=fft(charge_left)/layer%points
	  !charge_right=fft(charge_right)/layer%points

	  layer%ps_charge_boundary(:,n2)=charge_left+charge_right

	  END SUBROUTINE

	  SUBROUTINE priv_join_charges(layer,leftlayer,rightlayer,stars,cell)
	  USE m_gf_fft_singleton
	  IMPLICIT NONE
	  TYPE(t_chargelayer),INTENT(INOUT) :: layer
	  TYPE(t_chargelayer),INTENT(IN)    :: leftlayer
	  TYPE(t_chargelayer),INTENT(IN)    :: rightlayer
	  TYPE(t_stars),INTENT(IN)          :: stars
	  TYPE(t_cell),INTENT(IN)           :: cell

	  COMPLEX  :: chargelayer(3*stars%mx1,3*stars%mx2,2)
	  REAL     :: pos(3)
	  INTEGER  :: n1,n2,n3,i1,i2,in
      ! For all layers
      ! Put pseudo-charge onto real-space grid

      DO n2=1,stars%ng2
          layer%pseudo_charge(:,n2)=fft(layer%pseudo_charge(:,n2),inv=.true.)
          !layer%pseudo_charge(:,n2)=0.0

      ENDDO
	  DO n3=1,layer%points
	     chargelayer=0.0
	     IF (n3<=layer%points/2+1) THEN
	         pos(3)=layer%dt/(layer%points)*(n3-1)
	     ELSE
	         pos(3)=layer%dt/(layer%points)*(n3-1)-layer%dt
	     ENDIF
	     !Put charge on FFT grid
	     DO n1=-stars%mx1,stars%mx1
	          IF (n1<0) THEN
	             i1=3*stars%mx1+1+n1
	          ELSE
	             i1=n1+1
	          ENDIF
	          DO n2=-stars%mx2,stars%mx2
	            IF (n2<0) THEN
	                 i2=3*stars%mx2+1+n2
	            ELSE
	                 i2=n2+1
	            ENDIF
	            in=stars%ig(n1,n2,0)
	            if (in<1) cycle
	            in=stars%ig2(in)
	            chargelayer(i1,i2,1)=layer%pseudo_charge(n3,in)
	            !chargelayer(i1,i2,2)=layer%pseudo_charge(n3,in)
	            chargelayer(i1,i2,2)=layer%ps_charge_boundary(n3,in)
	            !print *,"No Charge"

	          ENDDO
         ENDDO
         !FFt to real space
         chargelayer(:,:,1)=fft(chargelayer(:,:,1))
         chargelayer(:,:,2)=fft(chargelayer(:,:,2))

         !Determine which charge to use (==step function)
         DO n1=1,stars%mx1*3
             IF (n1<=stars%mx1*3/2) THEN
                pos(1)=(n1-1)/(stars%mx1*3.0)
             ELSE
                pos(1)=(n1-1)/(stars%mx1*3.0)-1.0
             ENDIF
             DO n2=1,stars%mx2*3
                IF (n2<=stars%mx2*3/2) THEN
                    pos(2)=(n2-1)/(stars%mx2*3.0)
                ELSE
                    pos(2)=(n2-1)/(stars%mx2*3.0)-1.0
                ENDIF
                if (priv_insphere(layer,cell,pos,0.0)) cycle
                if (priv_insphere(leftlayer,cell,pos,(layer%c+leftlayer%c)/2.0).or.           &
                    priv_insphere(rightlayer,cell,pos,-(layer%c+rightlayer%c)/2.0)) THEN
                       chargelayer(n1,n2,1)=chargelayer(n1,n2,2)
                       cycle
                endif
                if (abs(pos(3))<=layer%c/2) cycle
                chargelayer(n1,n2,1)=chargelayer(n1,n2,2)
            ENDDO
         ENDDO
         !Put back on 2D FFT-Grid
         DO n1=1,stars%mx1*3
           DO n2=1,stars%mx2*3
              write(56,"(3i3,2e15.5)") n3,n2,n1,chargelayer(n1,n2,1)
         enddo
         enddo

         chargelayer(:,:,1)=real(chargelayer(:,:,1))
         chargelayer(:,:,1)=fft(chargelayer(:,:,1),inv=.true.)/size(chargelayer(:,:,1))
         layer%pseudo_charge(n3,:)=0.0
         DO n1=-stars%mx1,stars%mx1
	          IF (n1<0) THEN
	             i1=3*stars%mx1+1+n1
	          ELSE
	             i1=n1+1
	          ENDIF
	          DO n2=-stars%mx2,stars%mx2
	            IF (n2<0) THEN
	                 i2=3*stars%mx2+1+n2
	            ELSE
	                 i2=n2+1
	            ENDIF
	            in=stars%ig(n1,n2,0)
	            if (in<1) cycle
	            in=stars%ig2(in)
	            layer%pseudo_charge(n3,in)=chargelayer(i1,i2,1)
              ENDDO
         ENDDO
      ENDDO !n3 loop
      DO n2=1,stars%ng2
          layer%pseudo_charge(:,n2)=fft(layer%pseudo_charge(:,n2))/layer%points
      ENDDO

      END SUBROUTINE

      LOGICAL FUNCTION priv_insphere(layer,cell,pos,shift)
      IMPLICIT NONE
      TYPE(t_chargelayer),INTENT(IN):: layer
      type(t_cell       ),INTENT(IN):: cell
      real,INTENT(IN)               :: pos(3),shift

      real    :: p(3),d(3)
      INTEGER :: n,i1,i2

      priv_insphere=.true.

      DO i1=-1,1
        DO i2=-1,1
            p=pos+(/1.*i1,1.*i2,shift/)
            p(1:2)=matmul(cell%amat(1:2,1:2),p(1:2))
            DO n=1,size(layer%rmt2)
                d=p-layer%atpos(:,n)
                IF (dot_product(d,d)<=layer%rmt2(n)) return
            ENDDO
        ENDDO
      ENDDO
      priv_insphere=.false.
      END FUNCTION

      SUBROUTINE priv_make_potential(lay,stars,n2)
	  USE m_gf_fft_singleton
      IMPLICIT NONE
      TYPE(t_chargelayer),INTENT(INOUT):: lay
      TYPE(t_stars),INTENT(IN)         :: stars
      INTEGER,INTENT(IN)               :: n2

      !locals
      INTEGER  :: n
      REAL     :: g,charge_neutrality,s
      complex     :: g2(size(lay%potential)*2)
      ! first we have to make the system charge neutral

      lay%potential=lay%pseudo_charge(:,n2)
      !print *,"Charge neutrality!!"
      IF (n2==1) THEN
          !Double the grid to construct charge neutral system
          lay%potential=fft(lay%potential,inv=.true.)/lay%points
          g2(1:lay%points)=lay%potential
          g2(lay%points+1:)=lay%potential

          !construct compensation charge
          s=lay%points/8.
          IF (s<2) CALL juDFT_warn("Too small d-tilde")
          charge_neutrality=-1.*sum(g2)
          lay%potential=0.0
          DO n=1,lay%points/2-1
              lay%potential(lay%points/2+2-n)=exp(-1.*(n/s)**2)
              lay%potential(lay%points/2+1+n)=exp(-1.*(n/s)**2)
          ENDDO
          lay%potential=lay%potential/sum(lay%potential)*charge_neutrality
          g2(lay%points/2+1:lay%points/2+lay%points)=g2(lay%points/2+1:lay%points/2+lay%points)+lay%potential


          ! Now calculate potential by FFT
          g2=fft(g2)
          DO n=1,size(g2)
                g=n-1
                IF (n>size(g2)/2+1) g=n-size(g2)-1
                g=(g*2.*pi_const/lay%dt/2.0)**2
                IF (abs(g)<1E-90) THEN
                    g2(n)=0.0
                ELSE
                    g2(n)=g2(n)*4.*pi_const/g
                 !print *, "NoPot"
                ENDIF
          ENDDO
          g2=fft(g2,inv=.true.)/size(g2)


          !Now map to half system
          DO n=1,size(g2)
             write(33,"(i3,2e15.5)") n,g2(n)
          enddo

          lay%potential(:)=g2(size(g2)-lay%points+1:)
          lay%potential(:lay%points/2+1)=g2(:lay%points/2+1)

          DO n=1,lay%points
             write(33,"(i3,2e15.5)") n,lay%potential(n)
          enddo

          lay%potential=fft(lay%potential)
          lay%potential(1)=1.0
          return

          s=(1.-lay%d/lay%dt)*lay%points/3.
          s=min(s,5.)
          IF (s<2) CALL juDFT_warn("Too small d-tilde")
          charge_neutrality=0.0
          lay%potential=0.0
          DO n=1,lay%points/2-1
              lay%potential(lay%points/2+2-n)=exp(-1.*(n/s)**2)
              lay%potential(lay%points/2+1+n)=exp(-1.*(n/s)**2)
          ENDDO
          !lay%potential(lay%points/2+1:lay%points/2+2)=charge_neutrality/2.0
          lay%potential=fft(lay%potential)
          lay%potential=lay%potential*charge_neutrality/lay%potential(1)
          lay%potential=lay%pseudo_charge(:,1)-lay%potential
      ENDIF

      !now calculate the coloumb potential
      !cut-off high-freqency components!
      IF (size(lay%potential)-lay%gz_max(n2)>lay%gz_max(n2)+2) lay%potential(lay%gz_max(n2)+2:size(lay%potential)-lay%gz_max(n2))=0.0


      DO n=1,lay%points
         g=n-1
         IF (n>lay%points/2+1) g=n-lay%points-1
         g=(g*2.*pi_const/lay%dt)**2+stars%sk2(n2)**2
         !g>0 here!
           !print *,"NoPot"
          lay%potential(n)=lay%potential(n)*4.*pi_const/g
     ENDDO


     END SUBROUTINE



     SUBROUTINE priv_correction_potentials(lay,stars,n2,corr,c,dc)
     USE m_gf_fft_singleton
     IMPLICIT NONE
     TYPE(t_chargelayer),INTENT(IN)       ::lay
     TYPE(t_stars),INTENT(IN)             ::stars
     INTEGER,INTENT(IN)                   ::n2
	 COMPLEX,INTENT(OUT)                  ::corr(lay%points,2)
	 COMPLEX,INTENT(OUT)                  ::c(2,2),dc(2,2)



     REAL   :: g
     INTEGER:: n,i
     complex::igb,gg
	 COMPLEX::test(size(corr,1)),test1(size(corr,1)*8)
	 real   ::pos(8*size(corr,1))

	 integer:: points,n_shift


     !points=lay%gz_max(n2)*2+1
     points=lay%points
     n_shift=(1.-(lay%d)/lay%dt)*points/2-2
     n_shift=max(n_shift,1)
     !Construct the two correcting potentials
     corr=0.0

	 if (n2==1) then
	   !construct potential such that first derivative is ok
	   corr(:,1)=1.0
	   DO n=1,points/2
	      corr(n+1,2)=n
	      corr(points+1-n,2)=-n
	   ENDDO
	   corr(:,1)=fft(corr(:,1))/size(corr,1)
	   corr(:,2)=fft(corr(:,2))/size(corr,1)
	 else
	    corr(points/2+1+n_shift,1)=1.0

        corr(points/2+2-n_shift,2)=1.0

	    corr(:points,1)=fft(corr(:points,1))
        corr(:points,2)=fft(corr(:points,2))
        DO n=1,points
           g=n-1
           IF (n>points/2+1) g=n-points-1
           if (abs(g)>lay%gz_max(n2)) THEN
                 corr(n,:)=0.0
                 cycle
           endif
           g=(g*2.*pi_const/lay%dt)**2+stars%sk2(n2)**2
           corr(n,1)=corr(n,1)*4.*pi_const/g
           corr(n,2)=corr(n,2)*4.*pi_const/g
        ENDDO
     endif


     !calculate boundary values
     c(1,1)=dot_product(conjg(lay%project(:,1)),corr(:,1))
     c(2,1)=dot_product(conjg(lay%project(:,1)),corr(:,2))
     c(1,2)=dot_product(conjg(lay%project(:,2)),corr(:,1))
     c(2,2)=dot_product(conjg(lay%project(:,2)),corr(:,2))
     !calculate boundary derivatives
    	dc(1,1)=dot_product(conjg(lay%dproject(:,1)),corr(:,1))
     	dc(2,1)=dot_product(conjg(lay%dproject(:,1)),corr(:,2))
     	dc(1,2)=dot_product(conjg(lay%dproject(:,2)),corr(:,1))
     	dc(2,2)=dot_product(conjg(lay%dproject(:,2)),corr(:,2))
	return
	!rest is only debugging stuff
	test=fft(corr(:,1),inv=.true.)
	 do n=1,size(test)
	    g=n-1
	    if (n>size(test)/2+1) g=n-size(test)-1
	    write(55,"(4f10.3)") g,test(n)
	 enddo
	 pos=(/((n-1.0)/(size(pos)-1.),n = 1,size(pos))/)
	 pos=pos*lay%dt-lay%dt/2.0
	 test1=0.0
	 DO n = 1,lay%points
         gg = n-1
         IF (n>lay%points/2+1) gg=n-lay%points-1
         gg=cmplx(0.0,1.0)*gg*2.*pi_const/lay%dt
         test1(:)=test1(:)+corr(n,1)*gg*exp(gg*pos)
     enddo
     do n=1,size(test1)
         write(55,"(4f10.3)") pos(n),test1(n)
     enddo

     END SUBROUTINE
     SUBROUTINE priv_correction_potentials_OLD(lay,stars,n2,corr,c,dc)
     USE m_gf_fft_singleton
     IMPLICIT NONE
     TYPE(t_chargelayer),INTENT(IN)       ::lay
     TYPE(t_stars),INTENT(IN)             ::stars
     INTEGER,INTENT(IN)                   ::n2
	 COMPLEX,INTENT(OUT)                  ::corr(lay%points,2)
	 COMPLEX,INTENT(OUT)                  ::c(2,2),dc(2,2)



     REAL   :: g
     INTEGER:: n,i
     complex::igb,gg
	 COMPLEX::test(size(corr,1)),test1(size(corr,1)*8)
	 real   ::pos(8*size(corr,1))
	 real   :: smooth(size(corr,1)/2)

	 integer:: points

     n=(1.-(lay%d+0.5)/lay%dt)*lay%points/2
     if (n<3) THEN
        n=3
        CALL juDFT_warn("ratio of d/dtilde and no of points too small in intcoul")
     ENDIF
	 !Construct the two correcting potentials
     corr=0.0
	if (n2==1) then
	   !construct potential such that first derivative is ok
	   corr(:,1)=1.0
	   smooth=0.0
	   do i=1,2*n
	   	  g= i-n-0.5
	   	  smooth(i)=exp(-g**2/real(n**2)*2.5)
	   enddo
	   smooth=smooth/sum(smooth)
	   corr(lay%points/2-n+2:lay%points/2+n+1,1)=corr(lay%points/2-n+2:lay%points/2+n+1,1)-sum(corr(:,1))*smooth(:2*n)

	   corr(:,1)=fft(corr(:,1))

	   DO i=1,lay%points
	      g=i-1
	      if (i>lay%points/2+1) g=i-1-lay%points
	      if (abs(g)<1E-99) then
	            corr(i,1)=0.0
	            corr(i,2)=1.0
	      else
	       corr(i,:)=corr(i,1)/cmplx(0,g)
	      endif
	   enddo
	 else
     points=lay%points*sqrt(stars%gmax**2-stars%sk2(n2)**2)/stars%gmax
     points=min(points,lay%points)
     points=max(points,2)
     n=(1.-(lay%d+0.5)/lay%dt)*points/2
     n=max(n,1)
     !print *,n2,lay%points,points,n
	 corr=0.0
	 corr(points/2+2,1)=1.0
     corr(points/2+2+n,1)=-1.0

     corr(points/2+1,2)=1.0
     corr(points/2-n+1,2)=-1.0

	 corr(:points,1)=fft(corr(:points,1))
     corr(:points,2)=fft(corr(:points,2))
     DO n=1,points
        g=n-1
        IF (n>points/2+1) g=n-points-1
        g=(g*2.*pi_const/lay%dt)**2+stars%sk2(n2)**2
        corr(n,1)=corr(n,1)*4.*pi_const/g
        corr(n,2)=corr(n,2)*4.*pi_const/g
     ENDDO
     endif


     !calculate boundary values
     c(1,1)=dot_product(conjg(lay%project(:,1)),corr(:,1))
     c(2,1)=dot_product(conjg(lay%project(:,1)),corr(:,2))
     c(1,2)=dot_product(conjg(lay%project(:,2)),corr(:,1))
     c(2,2)=dot_product(conjg(lay%project(:,2)),corr(:,2))
     !calculate boundary derivatives
     if (n2==0) then
        dc=0.0
     	DO n = 1,lay%points
         gg = n-1
         IF (n>lay%points/2+1) gg=n-lay%points-1
         gg=cmplx(0.0,1.0)*gg*2.*pi_const/lay%dt
         dc=dc+corr(n,1)*gg
        enddo
     else
     	dc(1,1)=dot_product(conjg(lay%dproject(:,1)),corr(:,1))
     	dc(2,1)=dot_product(conjg(lay%dproject(:,1)),corr(:,2))
     	dc(1,2)=dot_product(conjg(lay%dproject(:,2)),corr(:,1))
     	dc(2,2)=dot_product(conjg(lay%dproject(:,2)),corr(:,2))
	endif
	return
	!rest is only debugging stuff
	test=fft(corr(:,1),inv=.true.)
	 do n=1,size(test)
	    g=n-1
	    if (n>size(test)/2+1) g=n-size(test)-1
	    write(55,"(4f10.3)") g,test(n)
	 enddo
	 pos=(/((n-1.0)/(size(pos)-1.),n = 1,size(pos))/)
	 pos=pos*lay%dt-lay%dt/2.0
	 test1=0.0
	 DO n = 1,lay%points
         gg = n-1
         IF (n>lay%points/2+1) gg=n-lay%points-1
         gg=cmplx(0.0,1.0)*gg*2.*pi_const/lay%dt
         test1(:)=test1(:)+corr(n,1)*gg*exp(gg*pos)
     enddo
     do n=1,size(test1)
         write(55,"(4f10.3)") pos(n),test1(n)
     enddo

     END SUBROUTINE

	SUBROUTINE priv_print_data(lay,pdata,n2)
	USE m_gf_fft_singleton
	IMPLICIT NONE
	TYPE(t_chargelayer),INTENT(IN):: lay(:)
	COMPLEX,INTENT(INOUT)            :: pdata(:,:,:)
	INTEGER,INTENT(IN)            :: n2

	INTEGER:: n,l,i
	real   :: z1,z2

	DO l=2,size(lay)-1
	    DO i=1,4
	        pdata(:lay(l)%points,i,l-1)=fft(pdata(:lay(l)%points,i,l-1),inv=.TRUE.)
	    ENDDO
	ENDDO
    WRITE(oUnit,*) "Potential for n2=",n2
    WRITE(oUnit,*) "Layer nz charge correction_charge potential smooth_potential"
    DO l=2,size(lay)-1
       DO n=1,lay(l)%points
           if ((n-1.)<=lay(l)%points/2.) then
              z1=lay(l)%dt/lay(l)%points*(n-1)
           else
              z1=lay(l)%dt/lay(l)%points*(n-1)-lay(l)%dt
           endif
           z2= sum(lay(2:l-1)%c)+lay(l)%c/2+z1
           WRITE(oUnit,"(2i5,999f15.9)") l-1,n,z1,z2,real(pdata(n,:,l-1))
       ENDDO
   ENDDO

	END SUBROUTINE

END MODULE m_gf_intcoul
