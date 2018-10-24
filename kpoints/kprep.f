      MODULE m_kprep
      CONTAINS
      SUBROUTINE  kprep(
     >                  iofile,iokpt,kpri,ktest,
     >                  nface,fnorm,fdist,iside,
     >                  vktes,nkstar,mkpt,mface,mdir,
     =                  nkrep,vkrep)
c ====================================================================
c
c    this subroutine
c    checks, if k-point vktes lies within irreducible wedge of BZ
c    which is characterized by
c          iside(i)= sign( (xvec,fnorm(i))-fdist(i) ) ;(i=1,nface )
c
c    IF TRUE, it returns a nonzero counter nkrep
c             and the new repr k-vektor vkrep(i) (i=1,3)
c
c    Meaning of variables:
c    INPUT:
c
c    representation of the irreducible part of the BZ:
c    fnorm    : normal vector of the planes bordering the irrBZ
c    fdist    : distance vector of the planes bordering the irrBZ
c    iside    : characterizing the inner side of each face of the irrBZ
c    nface    : number of faces of the irrBZ
c
c    k-point to be tested:
c    vktes    : k-point vector to be tested
c    nkstar   : index of star to which vktes belongs
c
c    OUTPUT: representative k-point
c    nkrep    : index (for each star nkstar); set to 1 if
c               representative k-point in irrBZ has been found
c    vkrep    : representative k-point in irrBZ for current star;
c               set to vktes, if condition fulfilled.
c ====================================================================
      IMPLICIT NONE
C
C-----> PARAMETER STATEMENTS
C
      INTEGER, INTENT (IN) :: mface,mkpt,mdir
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
c     integer  nsym,idsyst,idtype
C
C----->  BRAVAIS LATTICE INFORMATION
C
c     real     bltv(3,3)
C
C----->  RECIPROCAL LATTICE INFORMATION
C
      integer  nface
      real     fnorm(3,mface),fdist(mface)
c     real     xvec(3),rltv(3,3)
C
C----->  BRILLOUINE ZONE INTEGRATION
C
      integer  nkstar
C
C --->  local variables
c
      integer  i1,ii,ifac
      integer  iside(mface)
      integer  nkrep
      real     vkrep(3)
      real     ortest,vktes(3)
      real     invtpi, zero,one,half, eps,eps1
C
C --->  intrinsic functions
c
      intrinsic   real,abs
C
C --->  save and data statements
c
      save     one,zero,half,eps,eps1
      data     zero/0.0d0/,one/1.0d0/,half/0.5d0/,
     +         eps/1.0d-8/,eps1/1.0d-5/
c
c-----------------------------------------------------------------------
      if (kpri .ge. 3) then
        write(iofile,'(/)')
        write(iofile,'(3x,'' *<* kprep *>* '')')
        write(iofile,'(3x,'' check if k-vectors'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~'')')
        write(iofile,'(3x,'' are in irreducible wedge'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~'')')
        write(iofile,'(3x,'' of 1. Brillouin zone'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~'')')
        write(iofile,'(/)')
      end if
      if (ktest.ge. 5) then
        write(iofile,'(1x,i4,10x,''iofile'')') iofile
        write(iofile,'(1x,i4,10x,''iokpt'')')  iokpt
        write(iofile,'(1x,i4,10x,''kpri'')')  kpri
        write(iofile,'(1x,i4,10x,''ktest'')')  ktest
c       write(iofile,'(1x,3(f10.7,1x),10x,''xvec'')') (xvec(ii),ii=1,3)
c       write(iofile,'(1x,i4,10x,''ncorn'')')  ncorn
c       write(iofile,'(1x,i4,10x,''nedge'')')  nedge
        write(iofile,'(1x,i4,10x,''nface'')')  nface
        do 10 ifac = 1,nface
        write(iofile,'(1x,i4,1x,3(f10.7,1x),10x,''fnorm'')')
     +                 ifac,(fnorm(ii,ifac),ii=1,3)
        write(iofile,'(6x,f10.7,1x,20x,''fdist'')')
     +                      fdist(ifac)
        write(iofile,'(6x,i4,1x,26x,''iside'')')
     +                      iside(ifac)
 10     continue
        write(iofile,'(1x,i4,10x,''nkstar'')')  nkstar
c       write(iofile,'(1x,i4,10x,''nkpt'')')  nkpt
        write(iofile,'(/)')
      end if
c ======================================================================
c    start calculation
c
c
c ---> check if vktes lies in irred wedge of BZ
c      (i.e. on the same side of all boundary faces of irr wedge of BZ
c                                                              as xvec);
c
           do 70 ifac = 1,nface
              ortest = zero
            do 71 ii = 1,3
              ortest = ortest + vktes(ii)*fnorm(ii,ifac)
  71        continue
              ortest = ortest - fdist(ifac)
      if (ktest.ge. 4)
     +      write(iofile,'(1x,2(i4,2x),f10.7,10x,''ifac,iside,ortest'',
     +                         '' for vktes'')') ifac,iside(ifac),ortest
c
            if (abs(ortest) .lt. eps) go to 70
            if (ortest*iside(ifac) .lt. zero) go to 60
c
 70        continue
c
c    we have found a k-point inside irr BZ
c
c    (a) make sure it is not yet stored previously
c
            if(abs(vktes(1)-vkrep(1)).le.eps1
     +        .and. abs(vktes(2)-vkrep(2)).le.eps1
     +          .and. abs(vktes(3)-vkrep(3)).le.eps1
     +                                     .and. nkrep.gt.0) go to 60
c
               nkrep = nkrep+1
            do 80 ii = 1,3
              vkrep(ii) = vktes(ii)
 80        continue
c
      if (ktest.ge. 3) then
              write(iofile,'(1x,2(i4,1x),3(1x,f10.7),/,1x,
     +           ''nkstar,nkrep,vkrep(nkstar): kpoint in irr BZ'')')
     +            nkstar,nkrep,(vkrep(i1),i1=1,3)
      end if
c
 60   continue
c
      RETURN
      END SUBROUTINE kprep
      END MODULE m_kprep
