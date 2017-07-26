      MODULE m_wann_maxbnd
      USE m_fleurenv
c****************************************************************
c    determine maximum of number of bands
c****************************************************************
      CONTAINS
      SUBROUTINE wann_maxbnd(
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              isize,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,fullnkpts,irecl,
     >              l_bzsym,l_byindex,l_bynumber,l_byenergy,
     >              irreduc,odi,band_min,band_max,numbands,
     >              e1s,e2s,ef,nkpt,nbnd,l_gwf,iqpt)

      use m_cdnread, only:cdn_read
      use m_od_types,only:od_inp
      use m_wann_rw_eig

      IMPLICIT NONE
      integer,intent(in) :: lmaxd,ntypd,nlod,neigd,nvd,jspd
      integer,intent(in) :: isize,jspin,nbasfcn,nlotot,iqpt
      logical,intent(in) :: l_ss,l_noco,l_gwf
      integer,intent(in) :: nrec,fullnkpts,irecl
      logical,intent(in) :: l_byindex,l_bynumber,l_byenergy
      integer,intent(in) :: irreduc(fullnkpts)
      type (od_inp),intent(in) :: odi
      integer,intent(in) :: band_min,band_max,numbands
      logical,intent(in) :: l_bzsym
      real,intent(in)    :: e1s,e2s,ef
      integer,intent(in) :: nkpt

      integer,intent(out):: nbnd

      integer          :: ikpt,kptibz
      integer          :: nmat,nbands,nv(jspd)
      real             :: wk, bkpt(3),eig(neigd),cp_time(9)
      integer          :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      integer          :: nkbnd,i
      REAL             :: ello(nlod,ntypd,jspd),evdu(2,jspd)
      REAL             :: epar(0:lmaxd,ntypd,jspd)
      INTEGER          :: kveclo(nlotot)
      integer          :: n_start,n_end,co
      integer          :: num_bands
      
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX  :: z(nbasfcn,neigd)
#else
      REAL     :: z(nbasfcn,neigd)
#endif
      n_start=1
      n_end=neigd

      nbnd=0
      do ikpt = 1,fullnkpts

        kptibz=ikpt 
        if(l_bzsym) kptibz=irreduc(ikpt)
            call wann_read_eig(     
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              0,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,irecl,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time,66,
     >              l_gwf,iqpt)

        nkbnd = 0
        do i = 1,nbands
         if((eig(i).ge.e1s .and. nkbnd.lt.numbands.and.l_bynumber).or.
     &      (eig(i).ge.e1s.and.eig(i).le.e2s.and.l_byenergy ).or.
     &      (i.ge.band_min.and.i.le.band_max.and.l_byindex))  then
           nkbnd = nkbnd + 1
         endif 
        enddo
        if (nkbnd.ge.nbnd) nbnd = nkbnd

      enddo !ikpt

      end subroutine wann_maxbnd
      END MODULE m_wann_maxbnd
