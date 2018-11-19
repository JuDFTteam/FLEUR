      SUBROUTINE parawrite(&
     &                     sym,stars,atoms,sphhar,dimension,vacuum,obsolete,&
     &                     kpts,oneD,input)

      USE m_types
      IMPLICIT NONE
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_stars),INTENT(IN)     :: stars 
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_dimension),INTENT(IN) :: dimension
      TYPE(t_vacuum),INTENT(INOUT) :: vacuum
      TYPE(t_obsolete),INTENT(IN)  :: obsolete
      TYPE(t_kpts),INTENT(IN)      :: kpts
      TYPE(t_oneD),INTENT(IN)      :: oneD
      TYPE(t_input),INTENT(IN)     :: input
   
     

      write(6,*) "-----------fl7para file starts here-----------"

      WRITE (6,'(6x,''Symmetry elements and FFT parameters '')')
      WRITE (6,8080) sym%nop,stars%mx1,stars%mx2,stars%mx3,stars%ng3,stars%ng2

 8080 FORMAT (6x,'parameter (nop =',i2,',k1d=',i3,',k2d=',i3,',k3d=',i3,&
     &       ',n3d=',i6,',n2d=',i4,')')

!+sb(cdn_fft;Feb.97)
      WRITE (6,'(6x,''FFT-parameters for charge density'')')
      WRITE (6,8090) stars%kq1_fft,stars%kq2_fft,stars%kq3_fft

 8090 FORMAT (6x,'parameter (kq1d=',i3,',kq2d=',i3,',kq3d=',i3,')')

      WRITE (6,'(6x,''FFT-parameters for XC-potential'')')
      WRITE (6,8100) stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft

 8100 FORMAT (6x,'parameter (kxc1d=',i3,',kxc2d=',i3,',kxc3d=',i3,')')

      WRITE (6,'(6x,''(Inequivalent) atoms and radial mesh'')')
      WRITE (6,8110) atoms%ntype,atoms%nat,atoms%jmtd

 8110 FORMAT (6x,'parameter (ntypd=',i3,',natd=',i3,',jmtd=',i4,')')


      WRITE (6,'(6x,''Different lattice harmonics components'')')
      WRITE (6,8120) sphhar%ntypsd,sphhar%nlhd,sphhar%memd

 8120 FORMAT (6x,'parameter (ntypsd=',i3,',nlhd=',i3,',memd=',i2,')')

      WRITE (6,'(6x,''L-cutoff of potential, charge & wf '')')
      WRITE (6,8130) atoms%lmaxd

 8130 FORMAT (6x,'parameter (lmaxd=',i2,')')

      WRITE (6,'(6x,''Number of spins and vacua'')')
      WRITE (6,8140) input%jspins,vacuum%nvac

 8140 FORMAT (6x,'parameter (jspd=',i1,',nvacd=',i1,')')

      vacuum%nmz = 250
      vacuum%nmzxy = 100
      WRITE (6,'(6x,''Vacuum layers for G=0 and G=/=0'')')
      WRITE (6,8150) vacuum%nmz,vacuum%nmzxy

 8150 FORMAT (6x,'parameter (nmzd=',i3,',nmzxyd=',i3,')')

!+gu
      WRITE (6,'(6x,''3 & 2D planewaves, windows, k-points'')')
      WRITE (6,8180) dimension%nvd,dimension%nv2d,kpts%nkpt

 8180 FORMAT (6x,'parameter (nvd=',i5,',nv2d=',i4,',nwdd=1', ',nkptd=',i5,')')

      WRITE (6,'(6x,''Number of (occupied) bands'')')
      WRITE (6,8190) dimension%neigd,dimension%neigd

 8190 FORMAT (6x,'parameter (nobd=',i4,',neigd=',i4,')')
!-gu

      WRITE (6,'(6x,''Nuclear mesh and core levels'')')
      WRITE (6,8200) dimension%msh,dimension%nstd

 8200 FORMAT(6x,'parameter (msh=',i4,',nstd=',i2,')')

      WRITE (6,'(6x,''Max. l-value for pseudocharge exp.'')')
      WRITE (6,8210) dimension%ncvd

 8210 FORMAT (6x,'parameter (ncvd=',i3,')')

      WRITE (6,'(6x,''Layers for vacuum DOS'')')
      WRITE (6,'(6x,''parameter(layerd='',i3,'')'')') vacuum%layers

      WRITE (6,'(6x,''Local Orbital Parameters'')')
      WRITE (6,8220) atoms%nlod,atoms%llod

 8220 FORMAT(6x,'parameter (nlod=',i3,',llod=',i3,')')

      IF (oneD%odd%d1) THEN
        WRITE (6,'(6x,''One-dimensional parameters'')')
        WRITE (6,8230) oneD%odd%mb,oneD%odd%M,oneD%odd%m_cyl,oneD%odd%chi,oneD%odd%rot,&
     &                           oneD%odd%nop,oneD%odd%n2d,oneD%odd%d1
      ELSE
        WRITE (6,'(6x,''One-dimensional parameters'')')
        WRITE (6,8230) 0,0,0,0,0,sym%nop,stars%ng2,.false.
      END IF

 8230   FORMAT (6x,'parameter (vM=',i3,',MM=',i3,',m_cyl=',i3,&
     &                       ',chi=',i3,&
     &                       ',rot=',i3,',nop=',i3,',n2d=',i6,&
     &                       ',d1=',l1,')')
      WRITE(6,*) "-----------fl7para file ends here-----------"
 
      RETURN

      END
