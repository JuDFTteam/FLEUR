      SUBROUTINE parawrite(&
     &                     sym,stars,atoms,sphhar,vacuum,&
     &                     kpts ,input)

      USE m_types
      USE m_constants
      IMPLICIT NONE
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_stars),INTENT(IN)     :: stars 
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      
      TYPE(t_vacuum),INTENT(INOUT) :: vacuum
      TYPE(t_kpts),INTENT(IN)      :: kpts
       
      TYPE(t_input),INTENT(IN)     :: input
   
     

      write(oUnit,*) "-----------fl7para file starts here-----------"

      WRITE (oUnit,'(6x,''Symmetry elements and FFT parameters '')')
      WRITE (oUnit,8080) sym%nop,stars%mx1,stars%mx2,stars%mx3,stars%ng3,stars%ng2

 8080 FORMAT (6x,'parameter (nop =',i2,',k1d=',i3,',k2d=',i3,',k3d=',i3,&
     &       ',n3d=',i6,',n2d=',i4,')')

!+sb(cdn_fft;Feb.97)
      WRITE (oUnit,'(6x,''FFT-parameters for charge density'')')
      WRITE (oUnit,8090) stars%kq1_fft,stars%kq2_fft,stars%kq3_fft

 8090 FORMAT (6x,'parameter (kq1d=',i3,',kq2d=',i3,',kq3d=',i3,')')

      WRITE (oUnit,'(6x,''FFT-parameters for XC-potential'')')
      WRITE (oUnit,8100) stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft

 8100 FORMAT (6x,'parameter (kxc1d=',i3,',kxc2d=',i3,',kxc3d=',i3,')')

      WRITE (oUnit,'(6x,''(Inequivalent) atoms and radial mesh'')')
      WRITE (oUnit,8110) atoms%ntype,atoms%nat,atoms%jmtd

 8110 FORMAT (6x,'parameter (ntypd=',i3,',natd=',i3,',jmtd=',i4,')')


      WRITE (oUnit,'(6x,''Different lattice harmonics components'')')
      WRITE (oUnit,8120) sphhar%ntypsd,sphhar%nlhd,sphhar%memd

 8120 FORMAT (6x,'parameter (ntypsd=',i3,',nlhd=',i3,',memd=',i2,')')

      WRITE (oUnit,'(6x,''L-cutoff of potential, charge & wf '')')
      WRITE (oUnit,8130) atoms%lmaxd

 8130 FORMAT (6x,'parameter (lmaxd=',i2,')')

      WRITE (oUnit,'(6x,''Number of spins and vacua'')')
      WRITE (oUnit,8140) input%jspins,vacuum%nvacd

 8140 FORMAT (6x,'parameter (jspd=',i1,',nvacd=',i1,')')

      vacuum%nmz = 250
      vacuum%nmzxy = 100
      WRITE (oUnit,'(6x,''Vacuum layers for G=0 and G=/=0'')')
      WRITE (oUnit,8150) vacuum%nmz,vacuum%nmzxy

 8150 FORMAT (6x,'parameter (nmzd=',i3,',nmzxyd=',i3,')')

!+gu
      WRITE (oUnit,'(6x,''3 & 2D planewaves, windows, k-points'')')
      WRITE (oUnit,8180) lapw%dim_nvd(),lapw%dim_nv2d(),kpts%nkpt

 8180 FORMAT (6x,'parameter (nvd=',i5,',nv2d=',i4,',nwdd=1', ',nkptd=',i5,')')

      WRITE (oUnit,'(6x,''Number of (occupied) bands'')')
      WRITE (oUnit,8190) input%neig,input%neig

 8190 FORMAT (6x,'parameter (nobd=',i4,',neigd=',i4,')')
!-gu

      WRITE (oUnit,'(6x,''Nuclear mesh and core levels'')')
      WRITE (oUnit,8200) atoms%msh,29

 8200 FORMAT(6x,'parameter (msh=',i4,',nstd=',i2,')')


 8210 FORMAT (6x,'parameter (ncvd=',i3,')')

      WRITE (oUnit,'(6x,''Layers for vacuum DOS'')')
      WRITE (oUnit,'(6x,''parameter(layerd='',i3,'')'')') banddos%layers

      WRITE (oUnit,'(6x,''Local Orbital Parameters'')')
      WRITE (oUnit,8220) atoms%nlod,atoms%llod

 8220 FORMAT(6x,'parameter (nlod=',i3,',llod=',i3,')')

      
        WRITE (oUnit,'(6x,''One-dimensional parameters'')')
        WRITE (oUnit,8230) 0,0,0,0,0,sym%nop,stars%ng2,.false.
      

 8230   FORMAT (6x,'parameter (vM=',i3,',MM=',i3,',m_cyl=',i3,&
     &                       ',chi=',i3,&
     &                       ',rot=',i3,',nop=',i3,',n2d=',i6,&
     &                       ',d1=',l1,')')
      WRITE(oUnit,*) "-----------fl7para file ends here-----------"
 
      RETURN

      END
