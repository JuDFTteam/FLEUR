      MODULE m_ordstar
      use m_juDFT
c-----------------------------------------------------------------------
c
c --->   this program sorts
c        generated k-points in (not necessary full) stars
c        using the provided symmetry elements
c
c --->   order generated k-points in stars by applying symmetry:
c        - determine number of different stars nkstar .le. nkpt
c        - determine order of star iostar(kpn) .le. nsym
c        - assign pointer ikpn(i,ik); i=1,iostar(ik); ik=1,nkstar
c        - determine representative vector vkrep(3,ik); ik=1,nkstar
c        - for nreg=0:
c                    - assign nkpt= nkstar
c                    - assign vkxyz(ix,ik) = vkrep(ix,ik); ik=1,nkpt
c                                                          ix=1,3
c
c-----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE  ordstar(
     >                    iokpt,kpri,ktest,
     >                    fnorm,fdist,nface,iside,
     >                    nsym,ccr,rltv,mkpt,mface,mdir,
     =                    nkpt,vkxyz,
     <                    nkstar,iostar,ikpn,vkrep,nkrep)
cc
c    Meaning of variables:
c    INPUT:
c
c    Symmetry of lattice:
c    rltv     : cartesian coordinates of basis vectors for
c               reciprocal lattice rltv(ix,jn), ix=1,3; jn=1,3
c    nsym     : number of symmetry elements of points group
c    ccr     : rotation matrix for symmetry element
c                   in cartesian representation
c
c    representation of the irreducible part of the BZ:
c    fnorm    : normal vector of the planes bordering the irrBZ
c    fdist    : distance vector of the planes bordering the irrBZ
c    iside    : characterizing the inner side of each face of the irrBZ
c    nface    : number of faces of the irrBZ
c
c    k-point set:
c    nkpt     : number of k-points generated in set
c    vkxyz    : vector of kpoint generated; in cartesian representation
c
c    OUTPUT: Characteristics of k-point stars
c    nkstar   : number of stars for k-points generated by MOP
c    iostar   : number of k-points in each star
c    ikpn     : index field for the k-points in each star
c    vkrep    : representative k-vector in irrBZ for each star
c    nkrep    : index for each star;
c               1 if representative k-vector vkrep has been found
c-----------------------------------------------------------------------
      USE m_kprep
      IMPLICIT NONE
C
C-----> PARAMETER STATEMENTS
C
c
      INTEGER, INTENT (IN) :: mface,mkpt,mdir
c
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
      integer  nsym
      real     ccr(3,3,48)
C
C----->  RECIPROCAL LATTICE INFORMATION
C
      integer  nface
      real     rltv(3,3),fnorm(3,mface),fdist(mface)
C
C----->  BRILLOUINE ZONE INTEGRATION
C
      integer  nkpt,nkstar,iostar(mkpt)
      real     vkxyz(3,mkpt)
C
C --->  local variables
c
      character*80 blank
      integer  i,idim,i1,i2,i3,ii,ij,ik,is,isym,ifac, iik,iiik
      integer  ikc, i1red,nred
      integer  dirmin,dirmax,ndir1,ndir2,idir,lim(3)
      integer  kpl,kpm,kpn,nstnew
      integer  iplus,iminus,iside(mface),nkrep(mkpt),isi(7)
      real     orient,vkrep(3,mkpt),vktra(3)
      integer  ikpn(48,mkpt),irrkpn(mkpt),nfract(3),nleft,nirrbz
      real     fract(mkpt,3),fsig(2),vktes(3),vk_bzb(3)
      real     aivnkpt, sum,denom, t_bzb, t_len
      integer  isumkpt, iosub, ix, iy, iz, is1
      real     invtpi, zero,one,half, eps,eps1
      real     vkstar(3,48)
C
C --->  intrinsic functions
c
      intrinsic   abs,real
C
C --->  save and data statements
c
      save     one,zero,half,eps,eps1,iplus,iminus
      data     zero/0.0/,one/1.0/,half/0.5/,
     +         eps/1e-8/,eps1/1e-5/,iplus/1/,iminus/-1/
c
c-----------------------------------------------------------------------
c
c --->   set file numbers
c
c     iokpt  = 15
      iofile = 9
      OPEN (iofile,form='formatted',status='scratch')
      if (kpri .ge. 1) then
c       write(iofile,'(/)')
        write(iofile,'(3x,'' *<* ordstar *>* '')')
        write(iofile,'(3x,'' orders generated k-vectors'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~~~'')')
        write(iofile,'(3x,'' in (not neccessary full) stars'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')')
      end if
c
c --->   start calculation
c =====================================================================
c
c ---> set sign constants
       isi(1) = 0
       isi(2) = iminus
       isi(3) = iplus
       isi(4) = -2
       isi(5) =  2
       isi(6) = -3
       isi(7) =  3
c
          nkstar = 0
c
c --->   initialize pointers ikpn(ii,kpn) = 0
c               and          iostar(kpn) = 0
c
      do 10 kpn = 1,mkpt
       iostar(kpn) = 0
  10  continue
      do 15 ii = 1,nsym
       do 16 kpn = 1,mkpt
         ikpn(ii,kpn) = 0
  16   continue
  15  continue
c
c --->   screen over all kpoints
c
          do 100 kpn = 1,nkpt
c
c --->  check if given kpoint is assigned to the star of a previous one;
c             then skip
c
              do 105 ik = 1,nkstar
               do 106 is = 1,iostar(ik)
                 if(ikpn(is,ik) .eq. kpn) go to 100
 106           continue
 105          continue
c
c --->   new star is started for kpoint kpn: assign counters
c                            and set representative vkrep = zero
c
               nkstar = nkstar + 1
               is= 1
               iostar(nkstar)= 1
               ikpn(is,nkstar) = kpn
               nkrep(nkstar)= 0
               vkrep(1,nkstar) = zero
               vkrep(2,nkstar) = zero
               vkrep(3,nkstar) = zero
         if (ktest.ge.1) then
        write(iofile,'(/,1x,''start new star for vkxyz'')')
        write(iofile,'(1x,i4,1x,3(f10.7,1x),10x,''vkxyz'')')
     +                 is,(vkxyz(ii,kpn),ii=1,3)
         end if
c
c --->   generate kpoints establishing the star of vkxyz(kpn)
c
             do 110 isym = 1,nsym
               vktes(1) = ccr(1,1,isym)*vkxyz(1,kpn)
     +                    + ccr(1,2,isym)*vkxyz(2,kpn)
     +                      + ccr(1,3,isym)*vkxyz(3,kpn)
               vktes(2) = ccr(2,1,isym)*vkxyz(1,kpn)
     +                    + ccr(2,2,isym)*vkxyz(2,kpn)
     +                      + ccr(2,3,isym)*vkxyz(3,kpn)
               vktes(3) = ccr(3,1,isym)*vkxyz(1,kpn)
     +                    + ccr(3,2,isym)*vkxyz(2,kpn)
     +                      + ccr(3,3,isym)*vkxyz(3,kpn)
         if (ktest.ge.4) then
              write(iofile,'(1x,i4,1x,3(f10.7,1x),10x,''vktes'')')
     +                                           isym,(vktes(ii),ii=1,3)
         end if
c
c
c --->   check if any other kpoints belong to the star of vkxyz(kpn);
c              if so, assign new values to counters and pointers
c
               do 120 ik = kpn+1,nkpt
               if(abs(vktes(1)-vkxyz(1,ik)).le.eps1
     +           .and. abs(vktes(2)-vkxyz(2,ik)).le.eps1
     +             .and. abs(vktes(3)-vkxyz(3,ik)).le.eps1) then
c
c --->   -  make sure we have found a new point within current star
c
                    do 125 i1 = 1,is
                     if(ikpn(i1,nkstar) .eq. ik) go to 120
 125                continue
                     is= is+1
                     iostar(nkstar)= iostar(nkstar)+1
                     ikpn(is,nkstar) = ik
               end if
 120         continue
c
c --->   check, if vktes is representative k-point in irrBZ
c                                   for current star
            CALL kprep(
     >                 iofile,iokpt,kpri,ktest,
     >                 nface,fnorm,fdist,iside,
     >                 vktes,nkstar,mkpt,mface,mdir,
     =                 nkrep(nkstar),vkrep(1,nkstar))
c
 110         continue
 100      continue
c
c --->   ordering of kpoints into stars partly finished:
c
c          there are nkstar different k-stars of order iostar(ik);
c          the k-points belonging to the star are assigned to the
c          index-field ikpn(is,ik)
c          for the stars which lie inside the BZ a representative
c          vektor vkrep out of the irrBZ has been assigned.
c          We have to work on the points outside of the 1. BZ
c
      if (ktest.ge.2) then
c
c --->   printout of ordered k-points
c
        write(iofile,'(/,1x,i4,10x,''nkstar: number of stars'')') nkstar
      do 140 ik = 1,nkstar
        write(iofile,'(1x,i4,1x,i4,43x,
     +  ''ik, iostar(kpn): index and order of k-star'')')
     +       ik, iostar(ik)
        write(iofile,'(1x,''k-points in star:'')')
       do 150 is = 1,iostar(ik)
        write(iofile,'(1x,i4,3(1x,f10.7),1x,i4,10x,
     +         ''ikpn, vkxyz(kpn),is: kpoint and index in star'',/)')
     +          ikpn(is,ik), (vkxyz(i1,ikpn(is,ik)),i1=1,3), is
 150   continue
c
c       printout of representative vector vkrep of star
c
          if (nkrep(ik) .lt. 1)
     +      write(iofile,'(1x,''WARNING: we have found no '',
     +          ''k-point of star no '',i4, ''  inside irr BZ'')') ik
          if (nkrep(ik) .gt. 1)
     +      write(iofile,'(1x,''WARNING: we have found more than '',
     +      ''one k-point of star no '',i4, ''  inside irr BZ'')') ik
            write(iofile,'(1x,i4,3(1x,f10.7),1x,i4,10x,
     +         ''nkstar,vkrep(kpn),nkrep: repr k-point in irr BZ'',/)')
     +            ik,(vkrep(i1,ik),i1=1,3),nkrep(ik)
c
 140  continue
c
      end if
c
c --->   for those stars whose representative point in irrBZ have
c        not been found (nrep(ik)=0)
c        find representative point by shifting the first vector in star
c        by multiples of reciprocal latice vectors
c
      do 151 ik = 1,nkstar
          if (nkrep(ik) .eq. 0) then
c
           if (ktest.ge.2) then
             write(iofile,'(1x,i4,1x,i4,43x,
     +                ''ik, iostar(ik): index and order of k-star'')')
     +                                                    ik, iostar(ik)
           end if
c
      do 152 i3 = 1,7
        do 153 i2 = 1,7
          do 154 i1 = 1,7
                if (isi(i1).ne.0 .or. isi(i2).ne.0
     +                           .or. isi(i3).ne.0) then
                  vktra(1) = vkxyz(1,ikpn(1,ik)) + rltv(1,1)*isi(i1)
     +                                           + rltv(1,2)*isi(i2)
     +                                           + rltv(1,3)*isi(i3)
                  vktra(2) = vkxyz(2,ikpn(1,ik)) + rltv(2,1)*isi(i1)
     +                                           + rltv(2,2)*isi(i2)
     +                                           + rltv(2,3)*isi(i3)
                  vktra(3) = vkxyz(3,ikpn(1,ik)) + rltv(3,1)*isi(i1)
     +                                           + rltv(3,2)*isi(i2)
     +                                           + rltv(3,3)*isi(i3)
c
c --->   create star of vktra
c
         do 155 isym = 1,nsym
               vktes(1) = ccr(1,1,isym)*vktra(1)
     +                    + ccr(1,2,isym)*vktra(2)
     +                      + ccr(1,3,isym)*vktra(3)
               vktes(2) = ccr(2,1,isym)*vktra(1)
     +                    + ccr(2,2,isym)*vktra(2)
     +                      + ccr(2,3,isym)*vktra(3)
               vktes(3) = ccr(3,1,isym)*vktra(1)
     +                    + ccr(3,2,isym)*vktra(2)
     +                      + ccr(3,3,isym)*vktra(3)
              if (ktest.ge.4) then
                   write(iofile,'(1x,i4,1x,3(f10.7,1x),10x,''vktes'')')
     +                                           isym,(vktes(ii),ii=1,3)
              end if
c
c
c --->   check if vkrep of any other star coincides with vktes{vktra}
c              if so, assign new values to counters and pointers
c
         do 156 iik = 1,nkstar
            if (nkrep(iik) .ne. 0) then
               if(abs(vktes(1)-vkrep(1,iik)).le.eps1
     +              .and. abs(vktes(2)-vkrep(2,iik)).le.eps1
     +                .and. abs(vktes(3)-vkrep(3,iik)).le.eps1) then
c
c --->   -  assign k-point indices of star ik to star iik
c           and change counter
c
                    do 157 idir = 1,iostar(ik)
                     ikpn(iostar(iik)+idir,iik) = ikpn(idir,ik)
 157                continue
                     iostar(iik)= iostar(iik) + iostar(ik)
c
           if (ktest.ge.2) then
             write(iofile,'(1x,i4,1x,i4,43x,
     +              ''iik, iostar(iik): index and order of k-star'')')
     +                                                  iik, iostar(iik)
           end if
c
                go to 151
               end if
            end if
c
 156     continue
c
c --->   check, if vktes is in irrBZ
c               then vktes is representative k-point and
c                    current star remains a distinct star
c
            CALL kprep(
     >                 iofile,iokpt,kpri,ktest,
     >                 nface,fnorm,fdist,iside,
     >                 vktes,ik,mkpt,mface,mdir,
     =                 nkrep(ik),vkrep(1,ik))
c
                  if (nkrep(ik) .gt. 0) go to 151
c
 155     continue
c
                end if
 154      continue
 153    continue
 152  continue
          end if
 151  continue
c
c --->   -  reshuffle numbering of stars with nkrep > 0
c
                     nstnew = 0
                     isumkpt = 0
         do 158 iiik = 1,nkstar
              if (nkrep(iiik) .ne. 0) then
                     nstnew = nstnew + 1
                     isumkpt = isumkpt + iostar(iiik)
c
                   do 159 i1 = 1,iostar(iiik)
                     ikpn(i1,nstnew) = ikpn(i1,iiik)
 159               continue
                     iostar(nstnew)= iostar(iiik)
                     nkrep(nstnew) = nkrep(iiik)
                     vkrep(1,nstnew) = vkrep(1,iiik)
                     vkrep(2,nstnew) = vkrep(2,iiik)
                     vkrep(3,nstnew) = vkrep(3,iiik)
              end if
 158     continue
c
           if (ktest.ge.2) then
             write(iofile,'(/,'' result of ordering :'')')
             write(iofile,'(1x,i4,1x,i4,2x,
     +              ''no of stars, no of k-points contained in them'')')
     +                                                   nstnew, isumkpt
           end if
c
c --->   reduce number of stars
c
                nkstar = nstnew
c
c --->   final step:
c        check, if all representative vectors are distinct;
c        (sometimes in the previous step a different, but aquivalent vkrep
c        has been assigned) whose equivalence is found by translation 
c        and rotation the representative vectors vkrep
c
      do 1151 ik = nkstar,1,-1
ctest     if (nkrep(ik) .eq. 0) then
c
           if (ktest.ge.2) then
             write(iofile,'(1x,i4,1x,i4,43x,
     +                ''ik, iostar(ik): index and order of k-star'')')
     +                                                    ik, iostar(ik)
           end if
c
      do 1152 i3 = 1,3
        do 1153 i2 = 1,3
          do 1154 i1 = 1,3
                if (isi(i1).ne.0 .or. isi(i2).ne.0
     +                           .or. isi(i3).ne.0) then
                  vktra(1) = vkrep(1,ik) + rltv(1,1)*isi(i1)
     +                                           + rltv(1,2)*isi(i2)
     +                                           + rltv(1,3)*isi(i3)
                  vktra(2) = vkrep(2,ik) + rltv(2,1)*isi(i1)
     +                                           + rltv(2,2)*isi(i2)
     +                                           + rltv(2,3)*isi(i3)
                  vktra(3) = vkrep(3,ik) + rltv(3,1)*isi(i1)
     +                                           + rltv(3,2)*isi(i2)
     +                                           + rltv(3,3)*isi(i3)
c
c --->   create star of vktra
c
         do 1155 isym = 1,nsym
               vktes(1) = ccr(1,1,isym)*vktra(1)
     +                    + ccr(1,2,isym)*vktra(2)
     +                      + ccr(1,3,isym)*vktra(3)
               vktes(2) = ccr(2,1,isym)*vktra(1)
     +                    + ccr(2,2,isym)*vktra(2)
     +                      + ccr(2,3,isym)*vktra(3)
               vktes(3) = ccr(3,1,isym)*vktra(1)
     +                    + ccr(3,2,isym)*vktra(2)
     +                      + ccr(3,3,isym)*vktra(3)
              if (ktest.ge.4) then
                   write(iofile,'(1x,i4,1x,3(f10.7,1x),10x,''vktes'')')
     +                                           isym,(vktes(ii),ii=1,3)
              end if
c
c
c --->   check if vkrep of any other star coincides with vktes{vktra}
c              if so, assign new values to counters and pointers
c
         do 1156 iik = 1,nkstar
            if (iik .lt. ik) then
               if(abs(vktes(1)-vkrep(1,iik)).le.eps1
     +              .and. abs(vktes(2)-vkrep(2,iik)).le.eps1
     +                .and. abs(vktes(3)-vkrep(3,iik)).le.eps1) then
c
c --->   -  assign k-point indices of star ik to star iik
c           and change counter
c
                    do 1157 idir = 1,iostar(ik)
                     ikpn(iostar(iik)+idir,iik) = ikpn(idir,ik)
 1157                continue
                     iostar(iik)= iostar(iik) + iostar(ik)
c
           if (ktest.ge.2) then
             write(iofile,'(1x,i4,1x,i4,43x,
     +              ''iik, iostar(iik): index and order of k-star'')')
     +                                                  iik, iostar(iik)
           end if
c
                nkrep(ik) = 0
                go to 1151
               end if
            end if
c
 1156     continue
c
c --->   check, if vktes is in irrBZ
c               then vktes is representative k-point and
c                    current star remains a distinct star
c NOT NECESSARY ANYMORE
ctest            call  k p r e p
ctest     >                     (iofile,iokpt,kpri,ktest,
ctest     >                      nface,fnorm,fdist,iside,
ctest     >                      vktes,ik,
ctest     =                      nkrep(ik),vkrep(1,ik))
c
ctest                  if (nkrep(ik) .gt. 0) go to 1151
c
 1155     continue
c
                end if
 1154      continue
 1153    continue
 1152  continue
ctest     end if
 1151  continue
c
c --->   -  reshuffle numbering of stars with nkrep > 0
c
                     nstnew = 0
                     isumkpt = 0
         do 1158 iiik = 1,nkstar
              if (nkrep(iiik) .ne. 0) then
                     nstnew = nstnew + 1
                     isumkpt = isumkpt + iostar(iiik)
c
                   do 1159 i1 = 1,iostar(iiik)
                     ikpn(i1,nstnew) = ikpn(i1,iiik)
 1159               continue
                     iostar(nstnew)= iostar(iiik)
                     nkrep(nstnew) = nkrep(iiik)
                     vkrep(1,nstnew) = vkrep(1,iiik)
                     vkrep(2,nstnew) = vkrep(2,iiik)
                     vkrep(3,nstnew) = vkrep(3,iiik)
              end if
 1158     continue
c
           if (ktest.ge.2) then
             write(iofile,'(/,'' result of ordering :'')')
             write(iofile,'(1x,i4,1x,i4,2x,
     +              ''no of stars, no of k-points contained in them'')')
     +                                                   nstnew, isumkpt
           end if
c
c --->   reduce number of stars
c
                nkstar = nstnew
c
c --->   ordering of kpoints into stars totally finished:
c
c        - there are nkstar different k-stars of order iostar(ik);
c        - the k-points belonging to the star are assigned to the
c          index-field ikpn(is,ik)
c        - for every star a representative
c          vektor vkrep in the irrBZ has been assigned.
c
      if (ktest.ge.1) then
c
c    printout of ordered k-points
c
        write(iofile,'(/,1x,i4,10x,''nkstar: number of stars'')') nkstar
      do 1140 ik = 1,nkstar
        write(iofile,'(1x,i4,1x,i4,43x,
     +  ''ik, iostar(kpn): index and order of k-star'')')
     +       ik, iostar(ik)
        write(iofile,'(1x,''k-points in star:'')')
       do 1150 is = 1,iostar(ik)
        write(iofile,'(1x,i4,3(1x,f10.7),1x,i4,10x,
     +         ''ikpn, vkxyz(kpn),is: kpoint and index in star'',/)')
     +          ikpn(is,ik), (vkxyz(i1,ikpn(is,ik)),i1=1,3), is
1150   continue
c
c --->  printout of representative vector vkrep of star
c
          if (nkrep(ik) .lt. 1)
     +      write(iofile,'(1x,''WARNING: we have found no '',
     +          ''k-point of star no '',i4, ''  inside irr BZ; '',
     +          ''vkrep set to zero'')') ik
          if (nkrep(ik) .gt. 1)
     +      write(iofile,'(1x,''WARNING: we have found more than '',
     +      ''one k-point of star no '',i4, ''  inside irr BZ; '',
     +      ''last vkrep shown'')') ik
            write(iofile,'(1x,i4,3(1x,f10.7),1x,i4,10x,
     +           ''nkstar,vkrep(kpn),nkrep: repr k-point in irr BZ'')')
     +            ik,(vkrep(i1,ik),i1=1,3),nkrep(ik)
c
1140  continue
c
      end if
c
c --->   the index of every kpoint is contained in one of the stars.
c
c           check for ''left over'' k-points
c
      nleft = 0
      write(iofile,'(1x,''the following kpoints are not included'',
     + '' in any star'')')
      do 160 kpn = 1,nkpt
        do 165 ik = 1,nkstar
          do 166 is = 1,iostar(ik)
            if (kpn .eq. ikpn(is,ik)) go to 160
 166      continue
 165    continue
            write(iofile,'(1x,i4,3(1x,f10.7),15x,
     +                       ''kpn, vkxyz(kpn): kpoint '')')
     +                         kpn, (vkxyz(i1,kpn),i1=1,3)
             nleft = nleft+1
 160  continue
         if(nleft.eq.0) then
             write(iofile,'(/,1x,i4,20x,
     +                      ''no leftover points found'',/)')  nleft
c
         else
             write(iofile,'(/,1x,i4,20x,''WARNING: '',
     +      ''number of leftover points not equal zero'',/)')  nleft
!          CALL juDFT_error("leftover points",calledby="ordstar")
         end if
c
      CLOSE (iofile)
      
!
!-> check for bz-boundaries
!
      DO ik = 1, nkstar
!        write(*,'(a5,i3,a5,3f10.5)') 'star ',ik,' rep:',vkrep(:,ik)
        iosub = 0
        members: DO is = 1, iostar(ik)
         vk_bzb(:) = vkxyz(:,ikpn(is,ik))

         DO ix = -2,2
         DO iy = -2,2
         DO iz = -2,2
         DO is1 = is+1, iostar(ik)
         vk_bzb(:) = vkxyz(:,ikpn(is,ik)) + 
     +              ix * rltv(:,1) + iy * rltv(:,2) + iz * rltv(:,3)
         vk_bzb(:) = vk_bzb(:) - vkxyz(:,ikpn(is1,ik))  ! vkrep(:,ik)
         t_len = DOT_PRODUCT( vk_bzb, vk_bzb )
         IF (ABS(t_len) < 0.00001 ) THEN
           iosub = iosub + 1
           CYCLE members
         ENDIF
         ENDDO
         ENDDO
         ENDDO
         ENDDO

        ENDDO members
!        WRITE (*,*) 'iosub =',iosub
        iostar(ik) = iostar(ik)-iosub
      ENDDO

      RETURN
      END SUBROUTINE  ordstar
      END MODULE m_ordstar
