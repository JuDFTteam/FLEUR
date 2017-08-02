!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_rw_eig
      use m_juDFT
#ifdef CPP_MPI
#define CPP_HYBEIG
#endif
c****************************************************************
c    write WF1.eig and determine maximum of number of bands
c      Frank Freimuth, October 2006
c****************************************************************
      CONTAINS
      SUBROUTINE wann_write_eig(
     >              ntapwf,
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              isize,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,fullnkpts,irecl,
     >              l_bzsym,l_byindex,l_bynumber,l_byenergy,
     >              irreduc,odi,band_min,band_max,numbands,
     >              e1s,e2s,ef,l_paulimag,nkpt,    
     <              nbnd,kpoints,l_gwf,iqpt)

      use m_cdnread, only:cdn_read
      use m_types,only:od_inp

      IMPLICIT NONE
      integer,intent(in) :: ntapwf
      integer,intent(in) :: lmaxd,ntypd,nlod,neigd,nvd,jspd
      integer,intent(in) :: isize,jspin,nbasfcn,nlotot
      logical,intent(in) :: l_ss,l_noco
      integer,intent(in) :: nrec,fullnkpts
#ifndef CPP_MPI
      integer(8), intent (in) :: irecl
#else
      include 'mpif.h'
      INTEGER(KIND=MPI_OFFSET_KIND),intent(in) :: irecl
#endif
      logical,intent(in) :: l_byindex,l_bynumber,l_byenergy
      integer,intent(in) :: irreduc(fullnkpts)
      type (od_inp),intent(in) :: odi
      integer,intent(in) :: band_min,band_max,numbands
      logical,intent(in) :: l_bzsym
      real,intent(in)    :: e1s,e2s,ef
      logical,intent(in) :: l_paulimag
      integer,intent(in) :: nkpt

      integer,intent(out):: nbnd
      real,intent(out)   :: kpoints(fullnkpts)

      logical          :: l_eig
      character(len=3) :: spin12(2)
      data spin12/'WF1' , 'WF2'/
      integer          :: ikpt,kptibz
      integer          :: nmat,nbands,nv(jspd)
      real             :: wk, bkpt(3),eig(neigd),cp_time(9)
      integer          :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      integer          :: nkbnd,i
      REAL             :: ello(nlod,ntypd,jspd),evdu(2,jspd)
      REAL             :: epar(0:lmaxd,ntypd,jspd)
      INTEGER          :: kveclo(nlotot)
      integer          :: n_start,n_end,co
      integer          :: bnd1,bnd2,kpt
      complex, allocatable :: paulimat(:,:,:,:) 
      real             :: sum1,sum2,sum3
      integer          :: num_bands,err
      
! BEGIN QPOINTS
      LOGICAL, INTENT(IN) :: l_gwf
      INTEGER, INTENT(IN) :: iqpt
      CHARACTER(len=12) :: fname
! END QPOINTS
      real,parameter   :: hartree=27.21138505

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX  :: z(nbasfcn,neigd)
#else
      REAL     :: z(nbasfcn,neigd)
#endif
!      WRITE(*,*)'min',band_min,'max',band_max,'num',numbands
!      WRITE(*,*)'wann_rw_eig: neigd',neigd
      n_start=1
      n_end=neigd
      IF (l_gwf) THEN
         WRITE(fname,'("_",i4.4,".eig")')iqpt
         fname=spin12(jspin)//trim(fname)
      ELSE
         fname=spin12(jspin)//'.eig'
      ENDIF

!      WRITE(*,*)'wann_rw_eig.F: writing eig file ',fname

      inquire (file=fname,exist=l_eig)
      if (l_eig) then
         open (306,file=fname,
     &                               form='formatted',status='old')
         rewind (306)
      else
         open (306,file=fname,
     &                               form='formatted',status='new')
      endif!l_eig

      if(l_paulimag.and..false.)then
        num_bands=band_max-band_min+1
        allocate( paulimat(num_bands,num_bands,3,nkpt),stat=err )
        IF (err/=0)  CALL juDFT_error
     +       ("error allocating paulimat",calledby
     +       ="wann_rw_eig")
        open(655,file='wpmat')
        do ikpt=1,nkpt
          read(655,*)
          do co=1,3
            read(655,*)
            do bnd2=1,num_bands
               do bnd1=1,num_bands
                  read(655,*)paulimat(bnd1,bnd2,co,ikpt)
               enddo
            enddo
          enddo
        enddo      
        close(655)
      endif
      sum1=0.0
      sum2=0.0
      sum3=0.0

      nbnd=0
      do ikpt = 1,fullnkpts

        kptibz=ikpt 
        if(l_bzsym) kptibz=irreduc(ikpt)
            call wann_read_eig(
     >              ntapwf,     
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              0,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,irecl,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time,66,
     >              l_gwf,iqpt)


        nkbnd = 0
!        WRITE(*,*)'wann_rw_eig: nbands',nbands
        do i = 1,nbands
         if((eig(i).ge.e1s .and. nkbnd.lt.numbands.and.l_bynumber).or.
     &      (eig(i).ge.e1s.and.eig(i).le.e2s.and.l_byenergy ).or.
     &      (i.ge.band_min.and.i.le.band_max.and.l_byindex))  then
           nkbnd = nkbnd + 1
           if(l_bynumber.or.l_byindex)then
              write (306,'(2i12,f19.13)') nkbnd,ikpt,(eig(i)-ef)*
     &                                    hartree
           endif
           if(l_paulimag.and..false.)then
              if(eig(i).le.ef)then
                 sum1=sum1+paulimat(nkbnd,nkbnd,1,kptibz)
                 sum2=sum2+paulimat(nkbnd,nkbnd,2,kptibz)
                 sum3=sum3+paulimat(nkbnd,nkbnd,3,kptibz)
              endif
           endif   
         endif 
        enddo
!        WRITE(*,*)'wann_rw_eig: nkbnd',nkbnd,'nbnd',nbnd 
        if (nkbnd.ge.nbnd) nbnd = nkbnd

        if (odi%d1) kpoints(ikpt) = bkpt(3)

      enddo !ikpt

      if(l_paulimag)then
        write(6,*)"sum1=",sum1/fullnkpts
        write(6,*)"sum2=",sum2/fullnkpts
        write(6,*)"sum3=",sum3/fullnkpts
      endif  


      if(l_byenergy)then   !now we know the maximum of bands
         do ikpt=1,fullnkpts
            kptibz=ikpt
            if(l_bzsym)kptibz=irreduc(ikpt)
            call wann_read_eig(     
     >              ntapwf,       
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              0,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,irecl,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time,66,
     >              l_gwf,iqpt)


              nkbnd = 0
              do i = 1,nbands
                if(eig(i).ge.e1s .and. nkbnd.lt.nbnd)then
                     nkbnd = nkbnd + 1
                  write (306,'(2i12,f19.13)')nkbnd,ikpt,(eig(i)-ef)*
     &                                       hartree
                endif
              enddo 
         enddo   
      endif
      close (306)
      end subroutine wann_write_eig
c*********************************************************************
c     each process writes its own eig file if hdf is not used for mpi
c     this routine finds the appropriate eig file for a given kpoint
c     Frank Freimuth, October 2006
c*********************************************************************
      subroutine wann_mpi_read_eig(     
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              irank,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,irecl,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time,funit,zMat,
     >              l_gwf,iqpt)

      use m_cdnread, only:cdn_read
      USE m_types
      implicit none
      INTEGER, INTENT (IN) :: irank,isize,kptibz,nbasfcn,neigd,nlotot
      INTEGER, INTENT (IN) :: nrec,nvd,jspd,jspin
      INTEGER, INTENT (IN) :: lmaxd,ntypd,nlod,iqpt
      LOGICAL, INTENT (IN) :: l_ss,l_noco,l_gwf
#ifndef CPP_MPI
      integer(8), intent (in) :: irecl
#else
      include 'mpif.h'
      INTEGER(KIND=MPI_OFFSET_KIND),intent(in) :: irecl
#endif
      INTEGER, INTENT (OUT) :: nbands,nmat
      REAL,    INTENT (OUT) :: wk

      INTEGER, INTENT (OUT) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (OUT) :: nv(jspd),kveclo(nlotot)
      REAL,    INTENT (OUT) :: bkpt(3),eig(neigd)
      REAL,    INTENT (INOUT) :: cp_time(9)
      REAL,    INTENT (OUT) :: ello(nlod,ntypd,jspd),evdu(2,jspd)
      REAL,    INTENT (OUT) :: epar(0:lmaxd,ntypd,jspd)

      TYPE(t_zmat), INTENT (INOUT) :: zmat !z(nbasfcn,noccbd) !can be real/complex

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, INTENT (OUT) :: z(nbasfcn,neigd)
#else
      REAL,    INTENT (OUT) :: z(nbasfcn,neigd)
#endif

      INTEGER, INTENT(IN) :: funit

      character(len=12) :: fname,fend
      integer :: aoff,d1,d10
      integer :: cpu_index
      logical :: l_eig
      integer :: n_start,n_end

      n_start=1
      n_end=neigd

      fend=''
      if(l_gwf) write(fend,'("_",i4.4)')iqpt

      do cpu_index=0,isize-1

        if (mod(kptibz-1,isize).eq.cpu_index) then



        if(cpu_index.eq.irank)then

           call cdn_read(
     >              eig_id,
     >              nvd,jspd,irank,isize,kptibz,jspin,nbasfcn,
     >              l_ss,l_noco,neigd,n_start,n_end,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eig,zMat)

        else !eigenvalues are stored in the file of another process


           d10=int((cpu_index+0.5)/10)
           d1=mod(cpu_index,10)
           aoff=iachar('1')-1
           if(d10.ge.10)aoff=iachar('7')
           fname='eig'//achar(d10+aoff)//achar(d1+iachar('1')-1)
           fname=trim(fname)//trim(fend)
           l_eig=.false.
           !write(*,*)fname
c           if(irank.eq.0) write(*,*) 'opening file=',fname

           inquire(file=fname,exist=l_eig)
           IF(.NOT.l_eig)  CALL juDFT_error("one of eig files not found"
     +          ,calledby ="wann_rw_eig")
           OPEN(666,file=fname,access='direct',form='unformatted',
     &                            recl=irecl,status='old')
           call cdn_read(
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              cpu_index,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,kptibz,666,neigd,
     >              n_start,n_end,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time)

           close(666)
        endif   !cpu_index==irank
        endif !find process
        enddo !cpu_index
      end subroutine wann_mpi_read_eig  
c*****************************************************************
c     read in eig file
c*****************************************************************
      subroutine wann_read_eig(     
     >              ntapwf,
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              irank,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,irecl,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time,funit,
     >              l_gwf,iqpt)

      use m_cdnread, only:cdn_read
      implicit none
      integer, intent (in) :: ntapwf
      INTEGER, INTENT (IN) :: irank,isize,kptibz,nbasfcn,neigd,nlotot
      INTEGER, INTENT (IN) :: nrec,nvd,jspd,jspin
      INTEGER, INTENT (IN) :: lmaxd,ntypd,nlod,iqpt
      LOGICAL, INTENT (IN) :: l_ss,l_noco,l_gwf
#ifndef CPP_MPI
      integer(8), intent (in) :: irecl
#else
      include 'mpif.h'
      INTEGER(KIND=MPI_OFFSET_KIND),intent(in) :: irecl
#endif
      INTEGER, INTENT (OUT) :: nbands,nmat
      REAL,    INTENT (OUT) :: wk

      INTEGER, INTENT (OUT) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (OUT) :: nv(jspd),kveclo(nlotot)
      REAL,    INTENT (OUT) :: bkpt(3),eig(neigd)
      REAL,    INTENT (INOUT) :: cp_time(9)
      REAL,    INTENT (OUT) :: ello(nlod,ntypd,jspd),evdu(2,jspd)
      REAL,    INTENT (OUT) :: epar(0:lmaxd,ntypd,jspd)

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, INTENT (OUT) :: z(nbasfcn,neigd)
#else
      REAL,    INTENT (OUT) :: z(nbasfcn,neigd)
#endif

      INTEGER, INTENT(IN) :: funit    !QPOINTS

      integer :: n_start,n_end

      n_start=1
      n_end=neigd

#if( defined(CPP_MPI) && !defined(CPP_HDF) )
! MPI without HDF
      call wann_mpi_read_eig(     
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              irank,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,irecl,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time,funit,
     >              l_gwf,iqpt)

#else
!MPI with HDF or no MPI
        call cdn_read(
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              0,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,kptibz,funit,
     >              neigd,n_start,n_end,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time)

#endif
      END SUBROUTINE wann_read_eig  
      END MODULE m_wann_rw_eig
