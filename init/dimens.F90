      MODULE m_dimens
      use m_juDFT
      private
      public :: dimens
      CONTAINS
        SUBROUTINE dimens(&
             &                  mpi,input,sym,stars,&
             &                  atoms,sphhar,dimension,vacuum,&
             &                  obsolete,kpts,oneD,hybrid,jij)

          USE m_types
          USE m_dimen7
          USE m_firstglance
          USE m_constants
          IMPLICIT NONE
          TYPE(t_mpi),INTENT(INOUT) :: mpi
          TYPE(t_input),INTENT(INOUT) :: input
          TYPE(t_sym),INTENT(INOUT) :: sym
          TYPE(t_stars),INTENT(INOUT) :: stars 
          TYPE(t_atoms),INTENT(INOUT) :: atoms
          TYPE(t_sphhar),INTENT(INOUT) :: sphhar
          TYPE(t_dimension),INTENT(INOUT) :: dimension
          TYPE(t_vacuum),INTENT(INOUT) :: vacuum
          TYPE(t_obsolete),INTENT(INOUT) :: obsolete
          TYPE(t_kpts),INTENT(INOUT) :: kpts
          TYPE(t_oneD),INTENT(INOUT) :: oneD
          TYPE(t_hybrid),INTENT(INOUT) :: hybrid
          TYPE(t_Jij),INTENT(INOUT)    :: Jij

          TYPE(t_cell)     :: cell

          LOGICAL l_kpts,l_qpts,l_inpexist,ldum
          INTEGER n1,n2,n3,n4,n5,n6,n7,n8(3),n9,n10(3),i,j
          INTEGER i_vec(33)


#ifdef CPP_MPI
          INCLUDE 'mpif.h'
          INTEGER ierr(3)
#endif
          oneD%odd%d1=.TRUE.
          l_kpts=.TRUE.

          IF (mpi%irank.EQ.0) call priv_hello(version_const)

          WRITE (6,*) 'Your parameters: '


          OPEN (1,file='fl7para',form='formatted',status='old',err=200) ! La

          READ (1,*,ERR=200,END=200)     
          READ (1,901,ERR=200,END=200) sym%nop,stars%k1d,stars%k2d,stars%k3d,stars%n3d,stars%n2d
          IF (mpi%irank.EQ.0) WRITE (6,1001) sym%nop,stars%k1d,stars%k2d,stars%k3d,stars%n3d,stars%n2d
901       FORMAT (22x,i2,5x,i3,5x,i3,5x,i3,5x,i6,5x,i4)
1001      FORMAT (6x,'parameter (sym%nop= ',i2,',stars%k1d=',i3,',stars%k2d=',i3,',stars%k3d=',&
               &   i3,',stars%n3d=',i6,',stars%n2d=',i4,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,910,ERR=200,END=200) stars%kq1d,stars%kq2d,stars%kq3d
          IF (mpi%irank.EQ.0) WRITE (6,1010) stars%kq1d,stars%kq2d,stars%kq3d
910       FORMAT (22x,i3,6x,i3,6x,i3)
1010      FORMAT (6x,'parameter (stars%kq1d=',i3,',stars%kq2d=',i3,',stars%kq3d=',i3,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,914,ERR=200,END=200) stars%kxc1d,stars%kxc2d,stars%kxc3d
          IF (mpi%irank.EQ.0) WRITE (6,1014) stars%kxc1d,stars%kxc2d,stars%kxc3d
914       FORMAT (23x,i3,7x,i3,7x,i3)
1014      FORMAT (6x,'parameter (stars%kxc1d=',i3,',stars%kxc2d=',i3,',stars%kxc3d=',i3,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,902,ERR=200,END=200) atoms%ntypd,atoms%natd,atoms%jmtd
          IF (mpi%irank.EQ.0) WRITE (6,1002) atoms%ntypd,atoms%natd,atoms%jmtd
902       FORMAT (23x,i3,6x,i3,6x,i4)
1002      FORMAT (6x,'parameter (atoms%ntypd=',i3,',atoms%natd=',i3,',atoms%jmtd=',i4,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,903,ERR=200,END=200) sphhar%ntypsd,sphhar%nlhd,sphhar%memd
          IF (mpi%irank.EQ.0) WRITE (6,1003) sphhar%ntypsd,sphhar%nlhd,sphhar%memd
903       FORMAT (24x,i3,6x,i3,6x,i2)
1003      FORMAT (6x,'parameter (sphhar%ntypsd=',i3,',sphhar%nlhd=',i3,',sphhar%memd=',i2,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,904,ERR=200,END=200) atoms%lmaxd
          IF (mpi%irank.EQ.0) WRITE (6,1004) atoms%lmaxd
904       FORMAT (23x,i2)
1004      FORMAT (6x,'parameter (atoms%lmaxd=',i2,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,905,ERR=200,END=200) dimension%jspd,vacuum%nvacd
          IF (mpi%irank.EQ.0) WRITE (6,1005) dimension%jspd,vacuum%nvacd
905       FORMAT (22x,i1,7x,i1)
1005      FORMAT (6x,'parameter (dimension%jspd=',i1,',vacuum%nvacd=',i1,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,906,ERR=200,END=200) vacuum%nmzd,vacuum%nmzxyd
          IF (mpi%irank.EQ.0) WRITE (6,1006) vacuum%nmzd,vacuum%nmzxyd
906       FORMAT (22x,i3,8x,i3)
1006      FORMAT(6x,'parameter (vacuum%nmzd=',i3,',vacuum%nmzxyd=',i3,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,909,ERR=200,END=200) dimension%nvd,dimension%nv2d,kpts%nkptd
          IF (mpi%irank.EQ.0) WRITE (6,1009) dimension%nvd,dimension%nv2d,1,kpts%nkptd
909       FORMAT (21x,i5,6x,i4,6x,i1,7x,i5,6x,i4)
1009      FORMAT(6x,'parameter (nvd=',i5,',nv2d=',i4,',nwdd=1',',nkptd=',i5,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,911,ERR=200,END=200) dimension%neigd,dimension%neigd
          IF (mpi%irank.EQ.0) WRITE (6,1011) dimension%neigd,dimension%neigd
911       FORMAT (22x,i4,7x,i4)
1011      FORMAT(6x,'parameter (nobd=',i4,',dimension%neigd=',i4,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,912,ERR=200,END=200) dimension%msh,dimension%nstd
          IF (mpi%irank.EQ.0) WRITE (6,1012) dimension%msh,dimension%nstd
912       FORMAT (21x,i4,6x,i2)
1012      FORMAT(6x,'parameter (dimension%msh=',i4,',dimension%nstd=',i2,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,913,ERR=200,END=200) dimension%ncvd
          IF (mpi%irank.EQ.0) WRITE (6,1013) dimension%ncvd
913       FORMAT (22x,i3)
1013      FORMAT(6x,'parameter (dimension%ncvd=',i3,')')

          READ (1,*,ERR=200,END=200)     
          READ (1,915,ERR=200,END=200) vacuum%layerd
          IF (mpi%irank.EQ.0) WRITE (6,'(6x,''parameter(vacuum%layerd='',i3,'')'')') &
               &                       vacuum%layerd
915       FORMAT (23x,i3)

          READ (1,*,ERR=200,END=200)    
          READ (1,916,ERR=200,END=200) atoms%nlod,atoms%llod
          IF (mpi%irank.EQ.0) WRITE (6,1016) atoms%nlod,atoms%llod
916       FORMAT (22x,i3,6x,i3)
1016      FORMAT(6x,'parameter (atoms%nlod=',i3,',atoms%llod=',i3,')')
          !-odim
          IF (oneD%odd%d1) THEN
             READ (1,*,ERR=200,END=200)
             READ (1,917,ERR=200,END=200) oneD%odd%mb,oneD%odd%M,oneD%odd%m_cyl,oneD%odd%chi,&
                  &     oneD%odd%rot,oneD%odd%nop,oneD%odd%n2d,oneD%odd%d1
             IF (mpi%irank.EQ.0) WRITE (6,1017) oneD%odd%mb,oneD%odd%M,oneD%odd%m_cyl,oneD%odd%chi,&
                  &     oneD%odd%rot,oneD%odd%nop,oneD%odd%n2d,oneD%odd%d1
917          FORMAT (20x,i3,4x,i3,7x,i3,5x,i3,5x,i3,5x,i3,5x,i6,4x,l1)
1017         FORMAT (6x,'parameter (vM=',i3,',MM=',i3,',m_cyl=',i3,&
                  &     ',chi=',i3,&
                  &     ',rot=',i3,',sym%nop=',i3,',stars%n2d=',i6,',d1=',l1,')')
             !+odim
          END IF
          dimension%nspd=(atoms%lmaxd+1+mod(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)
          CLOSE (1)
          IF (mpi%irank.EQ.0) THEN

             INQUIRE (file='inp',exist=l_inpexist)
             IF (.not.l_inpexist) THEN
                INQUIRE (file='input',exist=l_inpexist)
                IF (.not.l_inpexist) THEN
                   CALL juDFT_error("no inp- or input-file found!",calledby&
                        &            ="dimens")
                ENDIF
                !         CALL inp_gen()
             ENDIF

             WRITE (6,*)
             CALL first_glance(&
                  &                    n1,n2,n3,n5,n6,input%itmax,&
                  &                    l_kpts,l_qpts,ldum,n7,n8,n9,n10)
             !
             IF (n1>atoms%ntypd)   CALL juDFT_error("atoms%ntypd  too small in fl7para" ,calledby ="dimens")
             IF (n2.LT.24) THEN
                IF (n2>sym%nop )    CALL juDFT_error("sym%nop   too small in fl7para" ,calledby ="dimens")
             ENDIF
             IF (n3>atoms%natd )   CALL juDFT_error("atoms%natd   too small in fl7para" ,calledby ="dimens")
             IF (n5>atoms%nlod )   CALL juDFT_error("atoms%nlod   too small in fl7para" ,calledby ="dimens")
             IF (n6>vacuum%layerd)  CALL juDFT_error("vacuum%layerd too small in fl7para" ,calledby ="dimens")
             IF ((.not.l_kpts).OR.(.not.l_qpts))  GOTO 201
          ENDIF

          RETURN

200       CONTINUE

          CLOSE (1)


#ifdef CPP_MPI
          CALL MPI_BARRIER(mpi%Mpi_comm,ierr)
#endif

201       IF (mpi%irank == 0) THEN
             IF (l_kpts) WRITE (6,*) ' No fl7para-file found, '
             WRITE (6,*) ' invoking dimen7... '
             !call first_glance to generate k-points
             CALL first_glance(&
                  &                    n1,n2,n3,n5,n6,input%itmax,&
                  &                    l_kpts,l_qpts,ldum,n7,n8,n9,n10)
             CALL dimen7(&
                  &              input,sym,stars,atoms,sphhar,&
                  &              dimension,vacuum,obsolete,kpts,&
                  &              oneD,hybrid,Jij,cell)

          ENDIF
          !     in case of a parallel calculation we have to broadcast
#ifdef CPP_MPI
          i_vec = (/sym%nop,stars%k1d,stars%k2d,stars%k3d,stars%n3d,stars%n2d,stars%kq1d,stars%kq2d,stars%kq3d,stars%kxc1d,stars%kxc2d,stars%kxc3d&
               &     ,atoms%ntypd,atoms%natd,atoms%jmtd,sphhar%ntypsd,sphhar%nlhd,sphhar%memd,atoms%lmaxd,dimension%jspd,vacuum%nvacd,dimension%nvd,dimension%nv2d&
               &     ,1,kpts%nkptd,dimension%nstd,dimension%neigd,dimension%msh,dimension%ncvd,vacuum%layerd,atoms%nlod,atoms%llod,input%itmax/)
          CALL MPI_BCAST(i_vec,33,MPI_INTEGER,0,mpi%Mpi_comm,ierr)
          sym%nop=i_vec(1);stars%k1d=i_vec(2);stars%k2d=i_vec(3);stars%k3d=i_vec(4);stars%n3d=i_vec(5)
          stars%n2d = i_vec(6);stars%kq1d=i_vec(7);stars%kq2d=i_vec(8);stars%kq3d=i_vec(9)
          stars%kxc1d = i_vec(10);stars%kxc2d = i_vec(11);stars%kxc3d = i_vec(12)
          atoms%ntypd = i_vec(13);atoms%natd =i_vec(14);atoms%jmtd=i_vec(15);sphhar%ntypsd=i_vec(16)
          sphhar%nlhd = i_vec(17);sphhar%memd=i_vec(18);atoms%lmaxd=i_vec(19);dimension%jspd=i_vec(20)
          vacuum%nvacd=i_vec(21);dimension%nvd=i_vec(22);dimension%nv2d=i_vec(23)
          kpts%nkptd = i_vec(25); dimension%nstd=i_vec(26);dimension%neigd=i_vec(27);dimension%msh=i_vec(28)
          dimension%ncvd=i_vec(29);vacuum%layerd=i_vec(30);atoms%nlod=i_vec(31);atoms%llod=i_vec(32)
          input%itmax=i_vec(33)
          CALL MPI_BCAST(oneD%odd%d1,1,MPI_LOGICAL,0,mpi%Mpi_comm,ierr)
          !      IF (odd%d1) THEN
          i_vec(:7) = (/oneD%odd%mb,oneD%odd%M,oneD%odd%m_cyl,oneD%odd%chi,oneD%odd%rot,oneD%odd%nop&
               &        ,oneD%odd%n2d/)
          CALL MPI_BCAST(i_vec,7,MPI_INTEGER,0,mpi%Mpi_comm,ierr)
          oneD%odd%mb = i_vec(1);oneD%odd%M = i_vec(2);oneD%odd%m_cyl=i_vec(3)
          oneD%odd%chi = i_vec(4);oneD%odd%rot = i_vec(5);oneD%odd%nop=i_vec(6)
          oneD%odd%n2d= i_vec(7)
          !      ELSE
          !         odd%nop = nop
          !      ENDIF
#endif
          dimension%nspd=(atoms%lmaxd+1+mod(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)
          vacuum%nmzd = 250
          vacuum%nmzxyd = 100


        END SUBROUTINE dimens


      SUBROUTINE priv_hello(ivers)
      IMPLICIT NONE
      CHARACTER(len=9), INTENT (IN)  :: ivers
      CHARACTER(len=9) :: cppflag(9)
      INTEGER          :: i,j
      WRITE (6,*) 'This output is generated by ',ivers,'  * * '
#if ( defined(CPP_AIX) )
      WRITE (6,*)  '                                    * \\:/ *'
#else
      WRITE (6,*)  '                                    * \:/ *'
#endif
      WRITE (6,*)  '                                    *  |  *'
      WRITE (6,*)  '                                      * *  '
      WRITE (6,*) 

        i = 0               ! First determine the architecture
#ifdef CPP_APC
        i = i + 1
        cppflag(i) = 'APC'
#endif
#ifdef CPP_DEC
        i = i + 1
        cppflag(i) = 'DEC'
#endif
#ifdef CPP_AIX
        i = i + 1
        cppflag(i) = 'AIX'
#endif
#ifdef CPP_T90
        i = i + 1
#ifdef CPP_MPI
        cppflag(i) = 'T3E'
#else
        cppflag(i) = 'T90'
#endif
#endif
        IF (i.GT.1) THEN 
          WRITE (6,*) 'You set compiler flags for more than one'
          WRITE (6,*) 'architecture: ', (cppflag(j),j=1,i)
          WRITE (6,*) 'Define only one system architecture! '
          CALL juDFT_error("Define only one system architecture! "&
     &         ,calledby ="dimens")
        ENDIF
        IF (i == 0) THEN
          WRITE (6,*) 'No system architecture specified in Makefile'
          cppflag(1) = 'GEN'
        ENDIF 
!
!       check for double precision etc.
!
#ifdef CPP_DOUBLE                
        cppflag(2) = 'DOUBLE'        
#else
        cppflag(2) = 'SINGLE'
#ifndef CPP_T90
        CALL juDFT_error(" define CPP_DOUBLE on non-Cray architectures!"&
     &       ,calledby ="dimens")
#endif
#endif
        WRITE (6,'(16a,4a,5a,7a,11a)') 'You compiled for ',&
     &        trim(cppflag(1)),' with ',trim(cppflag(2)),' precision,'
#ifdef CPP_INVERSION
        cppflag(1) = 'with'
#else
        cppflag(1) = 'without'
#endif
#ifdef CPP_SOC
        cppflag(2) = 'with'
#else
        cppflag(2) = 'without'
#endif
        WRITE (6,'(12a,7a,15a,7a,5a)') 'for systems ',&
     &  trim(cppflag(1)),' INVERSION and ',trim(cppflag(2)),' SOC.'
        i = 0
#ifdef CPP_MPI
        i = i + 1
        cppflag(i) = 'CPP_MPI'
#endif
#ifdef CPP_APW
        i = i + 1
        cppflag(i) = 'CPP_APW'
#endif
#ifdef CPP_CORE
        i = i + 1
        cppflag(i) = 'CPP_CORE'
#endif
#ifdef CPP_HTML
        i = i + 1
        cppflag(i) = 'CPP_HTML'
#endif
#ifdef CPP_HDF
        i = i + 1
        cppflag(i) = 'CPP_HDF'
#endif
#ifdef CPP_F90
        i = i + 1
        cppflag(i) = 'CPP_F90'
#endif
#ifdef CPP_WANN
        i = i + 1
        cppflag(i) = 'CPP_WANN'
#endif
#ifdef CPP_600
        i = i + 1
        cppflag(i) = 'CPP_600'
#endif
#ifdef CPP_GF
        i = i + 1
        cppflag(i) = 'CPP_GF'
#endif
#ifdef CPP_NOSPMVEC
        i = i + 1
        cppflag(i) = '+NOSPMVEC'
#endif
#ifdef CPP_IRAPPROX
        i = i + 1
        cppflag(i) = '+IRAPPROX'
#endif
        IF (i.GT.0) THEN
           WRITE (6,*) 'Additional flags are: ', (cppflag(j),j=1,i)
        ENDIF
      
      END SUBROUTINE

      END MODULE m_dimens
