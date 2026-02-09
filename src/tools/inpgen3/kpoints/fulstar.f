      MODULE m_fulstar
      use m_juDFT
      CONTAINS
      SUBROUTINE  fulstar(
     >                    iofile,iokpt,kpri,ktest,
     >                    ccr,nsym,
     >                    vkrep,nkstar,mkpt,mface,mdir,
     =                    nkpt,vkxyz,wghtkp)
c ====================================================================
c
c    this subroutine generates k-points in full stars
c    from representative k-points in irreducible wedge of BZ
c
c    Meaning of variables:
c    INPUT:
c
c    Symmetry elements of point group
c    nsym     : number of symmetry elements of points group
c    ccr     : rotation matrix for symmetry element
c                   in cartesian representation
c
c    representative k-points:
c    vkrep    : vkrep(ix,n); representative k-point vectors in irrBZ
c    nkstar   : number of representative k-vectors
c
c    OUTPUT: new k-point set
c    nkpt     : total number of k-points generated in full stars
c    vkxyz    : generated k-point vectors in cartesian representation
c    wghtkp   : (augmented) weight of vkxyz for BZ integration
c
c ====================================================================
      IMPLICIT NONE
C
C-----> PARAMETER STATEMENTS
C
       INTEGER, INTENT (IN) :: mkpt,mface,mdir
c
c ---> file number for read and write
c
      integer  iofile,iokpt
c
c ---> running mode parameter
c
      integer  kpri,ktest
C
C----->  Symmetry information
C
      integer  nsym,idsyst,idtype
      real     ccr(3,3,48)
C
C----->  BRAVAIS LATTICE INFORMATION
C
      real     bltv(3,3)
C
C----->  RECIPROCAL LATTICE INFORMATION
C
c     integer  ncorn,nface,nedge
c     real     xvec(3),rltv(3,3),fnorm(3,mface),fdist(mface)
C
C----->  BRILLOUINE ZONE INTEGRATION
C
      integer  nmop,nreg
      integer  nkpt,nkstar,ifstar(mkpt)
      real     vkxyz(3,mkpt),kzero(3),wghtkp(mkpt)
C
C --->  local variables
c
      character*80 blank
      integer  isumnkpt
      integer  i1,i2,i3,ii,ik,is,isym,ifac
      integer  dirmin,dirmax,ndir1,ndir2,idir,lim(3) ,nbound
      integer  iplus,iminus,iside(mface)
      integer  kpl,kpm,kpn,nstar(mdir)
      integer  ikpn(48,mkpt),irrkpn(mkpt),nirrbz,ntest
      real     vkrep(3,mkpt), vkstar(3,48), wght(mkpt)
      real     sumwght
      real     fract(mkpt),fsig(2),vktes(3)
      real     orient(mface),ortest
      real     aivnkpt,ainvnmop, sum,denom
      real     invtpi, zero,one,half, eps,eps1
C
C --->  intrinsic functions
c
      intrinsic   abs,max,real
C
C --->  save and data statements
c
      save     one,zero,half,eps,eps1,iplus,iminus
      data     zero/0.00/,one/1.00/,half/0.50/,
     +         eps/1.0e-8/,eps1/1.0e-9/,
     +         iplus/1/, iminus/-1/
c
c-----------------------------------------------------------------------
      if (kpri .ge. 3) then
        write(iofile,'(/)')
        write(iofile,'(3x,'' *<* fulstar *>* '')')
        write(iofile,'(3x,'' generate full stars of k-points'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')')
        write(iofile,'(3x,'' in 1. Brillouin zone'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~'')')
      end if
c
      if (ktest.ge. 5) then
        write(iofile,'(1x,i4,10x,''iofile'')') iofile
        write(iofile,'(1x,i4,10x,''iokpt'')')  iokpt
        write(iofile,'(1x,i4,10x,''kpri'')')  kpri
        write(iofile,'(1x,i4,10x,''ktest'')')  ktest
c       write(iofile,'(1x,3(f10.7,1x),10x,''xvec'')') (xvec(ii),ii=1,3)
c       write(iofile,'(1x,i4,10x,''ncorn'')')  ncorn
c       write(iofile,'(1x,i4,10x,''nedge'')')  nedge
c       write(iofile,'(1x,i4,10x,''nface'')')  nface
c       do 10 ifac = 1,nface
c       write(iofile,'(1x,i4,1x,3(f10.7,1x),10x,''fnorm'')')
c    +                 ifac,(fnorm(ii,ifac),ii=1,3)
c       write(iofile,'(1x,i4,1x,f10.7,1x,10x,''fdist'')')
c    +                 ifac,fdist(ifac)
c10     continue
        write(iofile,'(1x,i4,10x,''nkstar'')')  nkstar
        write(iofile,'(1x,i4,10x,''nkpt'')')  nkpt
      end if
c    printout heading
        write(iofile,'(1x,/,1x,''printout of generated stars '',
     >   ''of k-points'')')
        write(iofile,'(1x,i4,10x,''nkstar: number of stars '')') nkstar
c
c --->   save transferred wghtkp (applicable for set of vkrep in irrBZ)
c
             do 100 ik = 1,nkstar
                wght(ik) = wghtkp(ik)
100          continue
c
               nkpt = 0
c
             do 200 ik = 1,nkstar
c
c --->   assign repr k-vector vkrep(ik) to first k-point of full star
c
             do 205 ii   = 1,3
               vkstar(ii,1) = vkrep(ii,ik)
 205         continue
c
c --->   generate k-points establishing the symmetry star of vkrep(ik)
c
             do 210 isym = 2,nsym
               vkstar(1,isym) = ccr(1,1,isym)*vkrep(1,ik)
     +                        + ccr(1,2,isym)*vkrep(2,ik)
     +                        + ccr(1,3,isym)*vkrep(3,ik)
               vkstar(2,isym) = ccr(2,1,isym)*vkrep(1,ik)
     +                        + ccr(2,2,isym)*vkrep(2,ik)
     +                        + ccr(2,3,isym)*vkrep(3,ik)
               vkstar(3,isym) = ccr(3,1,isym)*vkrep(1,ik)
     +                        + ccr(3,2,isym)*vkrep(2,ik)
     +                        + ccr(3,3,isym)*vkrep(3,ik)
 210         continue
c
c
      if (ktest.ge. 5) then
               write(iofile,'(/,''star # '',i4,/)') ik
               write(iofile,'(1x,i4,3(1x,f10.7),10x,
     +                  '' vkrep(isym): represent k-point of star'')')
     +                                 1, (vkrep(i2,ik),i2=1,3)
             do 211 isym = 2,nsym
               write(iofile,'(1x,i4,3(1x,f10.7),15x,
     +              '' is, vkstar(isym): index and kpoint in star'')')
     +                                 isym, (vkstar(i2,isym),i2=1,3)
 211         continue
              write(iofile,'(/)')
       end if
c
c --->   eliminate equal k-points from symmetry star to form full star
c        generate ifstar(ik) accordingly
c
c        use an index field to assign the different k-points
c
                  ifstar(ik) = 1
                  ikpn(1,ik) = 1
               do 220 isym = 2,nsym
                  ikpn(isym,ik) = 0
 220           continue
c     scan over all generated points vkstar in current star
               do 225 isym = 2,nsym
c     compare with points already found to be distinct in current star
               do 226 is = 1,ifstar(ik)
               i1 = ikpn(is,ik)
               if   (abs(vkstar(1,i1)-vkstar(1,isym)).le.eps
     +         .and. abs(vkstar(2,i1)-vkstar(2,isym)).le.eps
     +         .and. abs(vkstar(3,i1)-vkstar(3,isym)).le.eps)
     +                                                      go to 225
 226         continue
c
c --->   we have found a distinct k-point
c
                  ifstar(ik) = ifstar(ik) +1
 
                  ikpn(ifstar(ik),ik) = isym
 225         continue
c
c --->   output of vectors in full star
c
           if (ktest.ge. 3) then
             write(iofile,'(1x,i4,1x,i4,43x,
     +           ''ik, ifstar(ik): index and order of full k-star'')')
     +                                                    ik, ifstar(ik)
             write(iofile,'(1x,''k-points in full star:'')')
            do 250 is = 1,ifstar(ik)
             write(iofile,'(1x,i4,3(1x,f10.7),15x,
     +                '' is, vkstar(is): index and kpoint in star'')')
     +                               is, (vkstar(i2,ikpn(is,ik)),i2=1,3)
 250        continue
c       output of parameters determining wghtkp
             write(iofile,'(1x,2(i4,1x),f17.14,10x,
     +                      '' ik, ifstar(ik),wght(ik)'')')
     +                                ik, ifstar(ik),wght(ik)
           end if
c
c --->   assign k-points and calculate weights
c              - assign vkxyz(ix,kpn) = vkstar(ix,ikpn(is,ik));
c                        ix=1,3; kpn=1,nkpt; ik=1,nstar; is=1,ifstar(ik)
c              - calculate wghtkp(kpn)=wghtkp_old(ik)/ifstar(ik)
c                                kpn=1,nkpt; ik=1,nstar
c
                 do 270 is = 1,ifstar(ik)
                     nkpt = nkpt + 1
                     vkxyz(1,nkpt) = vkstar(1,ikpn(is,ik))
                     vkxyz(2,nkpt) = vkstar(2,ikpn(is,ik))
                     vkxyz(3,nkpt) = vkstar(3,ikpn(is,ik))
                     wghtkp(nkpt)  = wght(ik)/real(ifstar(ik))
 270             continue
c
 200  continue
c
c
c --->   test sumrules for k-points in full stars
c
                 sumwght = zero
                isumnkpt = 0
             do 280 kpn = 1,nkpt
                 sumwght = sumwght + wghtkp(kpn)
 280             continue
c
                 do 290 ik = 1,nkstar
                    isumnkpt = isumnkpt + ifstar(ik)
 290             continue
c
      if (nkpt .ne. isumnkpt) then
         write(iofile,'(2(1x,i4),'' nkpt,isumnkpt do not coincide'')' )
     +                                         nkpt,isumnkpt
      else
         write(iofile,'(2(1x,i4),'' nkpt and isumnkpt do coincide'')' )
     +                                         nkpt,isumnkpt
      end if
c
      if (abs(sumwght-one) .gt. eps1) then
         write(iofile,'(1x,'' WARNING!!!!'')')
         write(iofile,'(1x,f17.10,1x,
     +   ''sumwght not equal one'')' ) sumwght
          CALL juDFT_error("sum wghtkp",calledby="fulstar")
      else
         write(iofile,'(1x,f12.10,1x,''abs(sumwght-one) le '',d10.3)' )
     +                                                      sumwght,eps1
      end if
c
      return
      END SUBROUTINE fulstar
      END MODULE m_fulstar
