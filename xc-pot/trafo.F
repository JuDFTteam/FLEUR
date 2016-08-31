!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_trafo

      CONTAINS

      SUBROUTINE waveftrafo_symm( 
     &                      cmt_out,z_out,cmt,z,bandi,ndb,
     &                      nk,iop,ntype,neq,
     &                      lmax,lmaxd,nindx,nkpt,nkpt3,maxlmindx,
     &                      natd,nw,nwd,bk,nop,nsym,mrot,tau,invtab,
     &                      taual,jsp,jspd,nbasfcn,neigd,map,tvec,d_wgn,
     &                      bmat,nvd,gpt,nv )

      USE m_constants     
      USE m_util          ,ONLY:  modulo1 
      USE m_wrapper

      IMPLICIT NONE

c     - scalars -
      INTEGER,INTENT(IN)      ::  ntype,nkpt,maxlmindx
      INTEGER,INTENT(IN)      ::  natd,nk,nw,nwd,jsp,jspd,
     &                            neigd,nbasfcn,ndb
      INTEGER,INTENT(IN)      ::  lmaxd,nvd,nop,nsym
      INTEGER,INTENT(IN)      ::  bandi,iop

c     - arrays -
      INTEGER,INTENT(IN)      ::  neq(ntype),lmax(ntype),
     &                            nindx(0:lmaxd,ntype),nv(jspd)
      INTEGER,INTENT(IN)      ::  mrot(3,3,nop)
      INTEGER,INTENT(IN)      ::  invtab(nop)
      INTEGER,INTENT(IN)      ::  gpt(3,nv(jsp))
      INTEGER,INTENT(IN)      ::  nkpt3(3)
      INTEGER,INTENT(IN)      ::  map(natd,nsym),tvec(3,natd,nsym)

      REAL,INTENT(IN)         ::  taual(3,natd),bk(3,nkpt),tau(3,nop)
      REAL,INTENT(IN)         ::  bmat(3,3)
      COMPLEX,INTENT(IN)      ::  d_wgn(-lmaxd:lmaxd,-lmaxd:lmaxd,
     +                                       0:lmaxd,nsym)
      COMPLEX,INTENT(IN)      ::  cmt(neigd,maxlmindx,natd)
#if ( defined(CPP_INVERSION) )
      REAL,INTENT(IN)         ::  z(nbasfcn,neigd)
#else
      COMPLEX,INTENT(IN)      ::  z(nbasfcn,neigd)
#endif
      COMPLEX,INTENT(OUT)     ::  cmt_out(maxlmindx,natd,ndb)
      COMPLEX,INTENT(OUT)     ::  z_out(nv(jsp),ndb)

c           - local -

c     - scalars -
      INTEGER                 ::  iatom,iatom1,iiatom,itype,igpt,igpt1,
     +                            ieq,ieq1,iiop
      INTEGER                 ::  i,l,n,nn,lm0,lm1,lm2,m1,m2
      COMPLEX                 ::  cdum,tpiimg
      COMPLEX,PARAMETER       ::  img=(0d0,1d0)

c     - arrays -
      INTEGER                 ::  rrot(3,3),invrrot(3,3)
      INTEGER                 ::  g(3),g1(3)
      REAL                    ::  tau1(3),rtaual(3),rkpt(3),rkpthlp(3),
     +                            trans(3)
      COMPLEX                 ::  cmthlp(2*lmaxd+1)
      LOGICAL                 ::  trs

      tpiimg  = -tpi_const*img

#ifdef CPP_INVERSION
      rrot    = transpose( mrot(:,:,invtab(iop)) )
      invrrot = transpose( mrot(:,:,iop) )
      trans   = tau(:,iop)
#else
      IF( iop .le. nop ) THEN
        trs     = .false.
        rrot    = transpose( mrot(:,:,invtab(iop)) )
        invrrot = transpose( mrot(:,:,iop) )
        trans   = tau(:,iop)       
      ELSE
        trs    = .true.
        iiop   = iop - nop
        rrot   = -transpose( mrot(:,:,invtab(iiop)) )
        invrrot= -transpose( mrot(:,:,iiop) )
        trans  = tau(:,iiop)
      END IF
#endif


      rkpt    = matmul(rrot,bk(:,nk))
      rkpthlp = rkpt
      rkpt    = modulo1(rkpt,nkpt3)
      g1      = nint(rkpt - rkpthlp)


      ! MT coefficients
      cmt_out = 0
      iatom   = 0
      iiatom  = 0

      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom  = iatom + 1

          iatom1 = map(iatom,iop)
          tau1   = tvec(:,iatom,iop)

          cdum   = exp( tpiimg * dotprod(rkpt,tau1) ) 

          lm0 = 0
          DO l = 0,lmax(itype)
            nn = nindx(l,itype)
            DO n = 1,nn
              lm1 = lm0 + n
              lm2 = lm0 + n + 2*l*nn
              DO i = 1,ndb
#ifdef CPP_INVERSION
                cmt_out(lm1:lm2:nn,iatom1,i) = cdum * 
     &                  matmul(cmt(bandi+i-1,lm1:lm2:nn,iatom),
     &                         d_wgn(-l:l,-l:l,l,iop)          )
#else
                IF( trs ) THEN
                  cmthlp(:2*l+1) = CONJG ( 
     +                                 cmt(bandi+i-1,lm1:lm2:nn,iatom) )
                ELSE
                  cmthlp(:2*l+1) =     cmt(bandi+i-1,lm1:lm2:nn,iatom)
                ENDIF
                cmt_out(lm1:lm2:nn,iatom1,i) = 
     +              cdum * matmul(cmthlp(:2*l+1),d_wgn(-l:l,-l:l,l,iop))
#endif
              END DO
            END DO
            lm0 = lm2
          END DO
        END DO
        iiatom = iiatom + neq(itype)
      END DO

      ! PW coefficients
      z_out = 0

      DO igpt = 1,nv(jsp)
        g    = matmul( invrrot,gpt(:,igpt)+g1 )
        !determine number of g
        igpt1 = 0
        DO i = 1,nv(jsp)
          IF ( maxval( abs( g - gpt(:,i) ) ) .le. 1E-06  ) THEN
            igpt1 = i
            EXIT
          END IF
        END DO 
        IF ( igpt1 .eq. 0 ) THEN 
          STOP 'wavetrafo_symm: rotated G vector not found'
        END IF
        cdum =  exp( tpiimg*dotprod(rkpt+gpt(:,igpt),trans(:)) )
#ifdef CPP_INVERSION
        z_out(igpt,1:ndb) = cdum * z(igpt1,bandi:bandi+ndb-1)
#else
        IF(trs) THEN
          z_out(igpt,1:ndb) = cdum * CONJG( z(igpt1,bandi:bandi+ndb-1) )
        ELSE
          z_out(igpt,1:ndb) = cdum * z(igpt1,bandi:bandi+ndb-1)
        END IF
#endif
      END DO


      END SUBROUTINE waveftrafo_symm



      SUBROUTINE waveftrafo_genwavf( 
     &                 cmt_out,z_out,cmt,z,nk,iop,ntype,neq,
     &                 lmax,lmaxd,nindx,nkpt,nkpt3,maxlmindx,
     &                 natd,nw,nwd,bk,nop,nsym,mrot,tau,invtab,map,tvec,
     &                 taual,jsp,jspd,nbasfcn,neigd,nbands,d_wgn,
     &                 bmat,nvd,gpt_nk,nv_nk,gpt_rkpt,nv_rkpt,phase)

      use m_juDFT
      USE m_constants   
      USE m_util          ,ONLY:  modulo1
      USE m_wrapper

      IMPLICIT NONE

c     - scalars -
      INTEGER,INTENT(IN)      ::  ntype,nkpt,maxlmindx,natd,nk,
     &                            nw,nwd,jsp,jspd,neigd,nbasfcn,nbands
      INTEGER,INTENT(IN)      ::  lmaxd,nvd,nop,nsym
      INTEGER,INTENT(IN)      ::  iop
      LOGICAL                 ::  phase
c     - arrays -
      INTEGER,INTENT(IN)      ::  neq(ntype),lmax(ntype),
     &                            nindx(0:lmaxd,ntype)
      INTEGER,INTENT(IN)      ::  mrot(3,3,nop)
      INTEGER,INTENT(IN)      ::  invtab(nop)
      INTEGER,INTENT(IN)      ::  nv_nk  (jspd),gpt_nk  (3,nv_nk  (jsp))
      INTEGER,INTENT(IN)      ::  nv_rkpt(jspd),gpt_rkpt(3,nv_rkpt(jsp))
      INTEGER,INTENT(IN)      ::  nkpt3(3)
      INTEGER,INTENT(IN)      ::  map(natd,nsym),tvec(3,natd,nsym)
      REAL,INTENT(IN)         ::  taual(3,natd),bk(3,nkpt),tau(3,nop)
      REAL,INTENT(IN)         ::  bmat(3,3)
      COMPLEX,INTENT(IN)      ::  d_wgn(-lmaxd:lmaxd,-lmaxd:lmaxd,
     +                                       0:lmaxd,nsym)

      COMPLEX,INTENT(IN)      ::  cmt(neigd,maxlmindx,natd)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL,INTENT(IN)         ::  z(nbasfcn,neigd)
      REAL,INTENT(OUT)        ::  z_out(nbasfcn,neigd)
#else
      COMPLEX,INTENT(IN)      ::  z(nbasfcn,neigd)
      COMPLEX,INTENT(OUT)     ::  z_out(nbasfcn,neigd)
#endif

      COMPLEX, INTENT(OUT)    ::  cmt_out(neigd,maxlmindx,natd)
c        - local - 

c     - scalars -
      INTEGER                 ::  itype,iatom,iatom1,iiatom,igpt,igpt1,
     +                            ieq,ieq1,iiop
      INTEGER                 ::  i,l,n,nn,lm0,lm1,lm2,m1,m2
      COMPLEX                 ::  cdum,tpiimg
      COMPLEX,PARAMETER       ::  img=(0d0,1d0)
      LOGICAL                 ::  trs

c     - arrays -
      INTEGER                 ::  rrot(3,3),invrrot(3,3)
      INTEGER                 ::  g(3),g1(3)
      REAL                    ::  tau1(3),rkpt(3),rkpthlp(3),trans(3)
      COMPLEX                 ::  zhlp(nbasfcn,neigd)
      COMPLEX                 ::  cmthlp(2*lmaxd+1)

      tpiimg = -tpi_const*img
#ifdef CPP_INVERSION
      rrot   = transpose( mrot(:,:,invtab(iop)) )
      invrrot= transpose( mrot(:,:,iop) )
      trans  = tau(:,iop)
#else
      IF( iop .le. nop ) THEN
        trs    = .false.
        rrot   = transpose( mrot(:,:,invtab(iop)) )
        invrrot= transpose( mrot(:,:,iop) )
        trans  = tau(:,iop)
      ELSE
        ! in the case of SOC (l_soc=.true.) 
        ! time reversal symmetry is not valid anymore;
        ! nsym should thus equal nop
        trs    = .true.
        iiop   = iop - nop
        rrot   = -transpose( mrot(:,:,invtab(iiop)) )
        invrrot= -transpose( mrot(:,:,iiop) )
        trans  = tau(:,iiop)
      END IF
#endif

      rkpt    = matmul(rrot,bk(:,nk))
      rkpthlp = rkpt
      rkpt    = modulo1(rkpt,nkpt3)
      g1      = nint(rkpt - rkpthlp)

      ! MT coefficients
      cmt_out = 0
      iatom   = 0
      iiatom  = 0

      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom  = iatom + 1

          iatom1 = map(iatom,iop)
          tau1   = tvec(:,iatom,iop)

          cdum   = exp( tpiimg*dotprod(rkpt,tau1) ) 

          lm0 = 0
          DO l = 0,lmax(itype)
            nn = nindx(l,itype)
            DO n = 1,nn
              lm1 = lm0 + n
              lm2 = lm0 + n + 2*l*nn

              DO i = 1,nbands
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
                cmt_out(i,lm1:lm2:nn,iatom1) = cdum * 
     &                     matmul( cmt(i,lm1:lm2:nn,iatom),
     &                             d_wgn(-l:l,-l:l,l,iop)   )
#else
                IF( trs ) THEN
                  cmthlp(:2*l+1) = conjg ( cmt(i,lm1:lm2:nn,iatom) )
                ELSE
                  cmthlp(:2*l+1) =         cmt(i,lm1:lm2:nn,iatom)
                END IF
                cmt_out(i,lm1:lm2:nn,iatom1) = cdum * 
     +                     matmul(cmthlp(:2*l+1),d_wgn(-l:l,-l:l,l,iop))
#endif


              END DO
            END DO
            lm0 = lm2
          END DO
        END DO
        iiatom = iiatom + neq(itype)
      END DO

      ! PW coefficients

      zhlp = 0
      DO igpt = 1,nv_rkpt(jsp)
        g    = matmul( invrrot,gpt_rkpt(:,igpt)+g1 )
        !determine number of g
        igpt1 = 0
        DO i=1,nv_nk(jsp)
          IF ( maxval( abs( g - gpt_nk(:,i) ) ) .le. 1E-06  ) THEN
            igpt1 = i
            EXIT
          END IF
        END DO 
        IF ( igpt1 .eq. 0 ) CYCLE
        cdum = exp(tpiimg*dotprod(rkpt+gpt_rkpt(:,igpt),trans ) )
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
        zhlp(igpt,:nbands)   =  cdum*z(igpt1,:nbands)
#else
        IF(trs) THEN
          zhlp(igpt,:nbands) =  cdum*conjg( z(igpt1,:nbands) )
        ELSE
          zhlp(igpt,:nbands) =  cdum*z(igpt1,:nbands)
        END IF
#endif
      END DO

      ! If phase and inversion-sym. is true, 
      ! define the phase such that z_out is real.

      IF( phase ) THEN
        DO i=1,nbands
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )

          CALL commonphase(cdum,zhlp(:,i),nbasfcn)

          IF(any(abs(aimag(zhlp(:,i)/cdum)).gt.1d-8)) THEN
            WRITE(*,*) maxval(abs(aimag(zhlp(:,i)/cdum)))
            WRITE(*,*) zhlp
            STOP 'waveftrafo1: Residual imaginary part.'
          END IF
          z_out(:,i)   = zhlp(:,i)   / cdum
          cmt_out(i,:,:) = cmt_out(i,:,:) / cdum
#else
          z_out(:,i) = zhlp(:,i)
#endif
        END DO
      ELSE
        z_out = zhlp
      END IF


      END SUBROUTINE waveftrafo_genwavf


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Symmetrizes MT part of input matrix according to inversion symmetry.
      ! This is achieved by a transformation to
      !        1/sqrt(2) * ( exp(ikR) Y_lm(r-R) + (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) )
      ! and                                                                                 if R /=0 or m<0
      !        i/sqrt(2) * ( exp(ikR) Y_lm(r-R) - (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) ) .
      !
      !  or
      !        i*Y_l,0(r)                                                                   if R=0,m=0 and l odd
      ! These functions have the property f(-r)=f(r)* which makes the output matrix real symmetric.
      ! (Array mat is overwritten! )

      SUBROUTINE symmetrize( mat,dim1,dim2,imode,lreal,
     &                       ntypd,ntype,neq,lcutm,maxlcutm,
     &                       nindxm,natd,invsat,invsatnr)

      IMPLICIT NONE

c     - scalars -
      INTEGER,INTENT(IN)    ::  imode,dim1,dim2
      INTEGER,INTENT(IN)    ::  ntypd,ntype,natd,maxlcutm
      LOGICAL,INTENT(IN)    ::  lreal

c     - arrays -
      INTEGER,INTENT(IN)    ::  neq(ntypd),lcutm(ntypd)
      INTEGER,INTENT(IN)    ::  nindxm(0:maxlcutm,ntypd)
      INTEGER,INTENT(IN)    ::  invsat(natd),invsatnr(natd)
      COMPLEX,INTENT(INOUT) ::  mat(dim1,dim2)

c     -local scalars -
      INTEGER               ::  i,j,itype,ieq,ic,ic1,i1,i2,l,m,n,nn,
     +                          ifac,ishift
      REAL                  ::  rfac,rdum,rmax
      COMPLEX               ::  img = (0d0,1d0)

c     - local arrays -
      COMPLEX               ::  carr(max(dim1,dim2)),cfac

      rfac = sqrt(0.5d0)
      cfac = sqrt(0.5d0)*img
      ic   = 0
      i    = 0

      DO itype = 1,ntype
        nn  = sum( (/ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) /) )
        DO ieq = 1,neq(itype)
          ic  = ic + 1
          IF( invsat(ic) .eq. 0) THEN
            ! if the structure is inversion-symmetric, but the equivalent atom belongs to a different unit cell
            ! invsat(atom) = 0, invsatnr(atom) = 0
            ! but we need invsatnr(atom) = natom
            ic1 = ic
          ELSE 
            ic1 = invsatnr(ic)
          END IF
          !ic1 = invsatnr(ic)
          IF( ic1 .lt. ic ) THEN
            i = i + nn
            CYCLE
          END IF
!           IF( ic1 .lt. ic ) cycle
          DO l = 0,lcutm(itype)
            ifac = -1
            DO m = -l,l
              ifac   =  -ifac
              ishift = (ic1-ic)*nn - 2*m*nindxm(l,itype)
              DO n = 1,nindxm(l,itype)
                i = i + 1
                j = i + ishift
                IF( ic1 .ne. ic .or. m .lt. 0 ) THEN
                  IF( iand(imode,1).ne.0 ) THEN
                    carr(:dim2) = mat(i,:)
                    mat(i,:)    = 
     +                       ( carr(:dim2) + ifac * mat(j,:) ) * rfac
                    mat(j,:)    = 
     +                       ( carr(:dim2) - ifac * mat(j,:) ) * (-cfac)
                  END IF
                  IF( iand(imode,2).ne.0) THEN
                    carr(:dim1) = mat(:,i)
                    mat(:,i)    = 
     +                       ( carr(:dim1) + ifac * mat(:,j) ) * rfac
                    mat(:,j)    = 
     +                       ( carr(:dim1) - ifac * mat(:,j) ) * cfac
                  END IF
                ELSE IF( m .eq. 0 .and. ifac .eq. -1) THEN
                  IF( iand(imode,1) .ne. 0) THEN
                    mat(i,:) = -img * mat(i,:)
                  END IF
                  IF( iand(imode,2) .ne. 0) THEN
                    mat(:,i) =  img * mat(:,i)
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO

      IF(lreal) THEN
!         ! Determine common phase factor and devide by it to make the output matrix real.
!         rmax = 0
!         DO i = 1,dim1
!           DO j = 1,dim2 
!             rdum = abs(real(mat(i,j)))+abs(aimag(mat(i,j)))
!             IF(rdum.gt.1d-6) THEN
!               cfac = mat(i,j)/abs(mat(i,j))
!               GOTO 1
!             ELSE IF(rdum.gt.rmax) THEN
!               cfac = mat(i,j)/abs(mat(i,j))
!               rmax = rdum
!             END IF
!           END DO
!         END DO
!         IF(1-abs(cfac)   .gt.1d-8) THEN ; mat = 0 ; RETURN ; END IF
!  1      IF(abs(1-cfac**2).gt.1d-8) mat = mat/cfac
!
!       IF(any(abs(aimag(mat)).gt.1d-8)) THEN
!         WRITE(*,*) maxval(aimag(mat))
!         STOP 'symmetrize: Residual imaginary part. Symmetrization failed.'

        ! Determine common phase factor and divide by it to make the output matrix real.
        CALL commonphase(cfac,mat,dim1*dim2)
        mat = mat / cfac
        IF(any(abs(aimag(mat)).gt.1d-8)) 
     +STOP 'symmetrize: Residual imaginary part. Symmetrization failed.'
      END IF

      END SUBROUTINE symmetrize

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! Undoes symmetrization with routine symmetrize.
      SUBROUTINE desymmetrize(mat,dim1,dim2,imode,
     &                        ntypd,ntype,neq,lcutm,maxlcutm,
     &                        nindxm,natd,invsat,invsatnr)

      IMPLICIT NONE

c     - scalars -
      INTEGER,INTENT(IN)      ::  imode,dim1,dim2
      INTEGER,INTENT(IN)      ::  ntypd,ntype,maxlcutm,natd

c     - arrays -
      INTEGER,INTENT(IN)      ::  neq(ntypd),lcutm(ntypd)
      INTEGER,INTENT(IN)      ::  nindxm(0:maxlcutm,ntypd)
      INTEGER,INTENT(IN)      ::  invsat(natd),invsatnr(natd)
      COMPLEX,INTENT(INOUT)   ::  mat(dim1,dim2)

c     - local scalars -
      INTEGER                 ::  ifac,i,j,itype,ieq,ic,ic1,i1,i2,l,m,n,
     +                            nn,ishift
      REAL                    ::  rfac1,rfac2
      COMPLEX                 ::  img = (0d0,1d0)
c     - local arrays -
      COMPLEX                 ::  carr(max(dim1,dim2))

      rfac1 = sqrt(0.5d0)
      ic    = 0
      i     = 0
      DO itype = 1,ntype
        nn  = sum( (/ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) /) )
        DO ieq = 1,neq(itype)
          ic  = ic + 1
          IF( invsat(ic) .eq. 0) THEN
            ! if the structure is inversion-symmetric, but the equivalent atom belongs to a different unit cell
            ! invsat(atom) = 0, invsatnr(atom) =0
            ! but we need invsatnr(atom) = natom
            ic1 = ic
          ELSE 
            ic1 = invsatnr(ic)
          END IF
          !ic1 = invsatnr(ic)
          !IF( ic1 .lt. ic ) cycle
          IF ( ic1.lt.ic ) THEN
            i = i + nn
            CYCLE
          END IF
          DO l = 0,lcutm(itype)
            ifac = -1
            DO m = -l,l
              ifac  = -ifac
              rfac2  = rfac1 * ifac
              ishift = (ic1-ic)*nn - 2*m*nindxm(l,itype)
              DO n = 1,nindxm(l,itype)
                i = i + 1
                j = i + ishift
                IF( ic1 .ne. ic .or. m .lt. 0) THEN
                  IF( iand(imode,1) .ne.0) THEN
                    carr(:dim2) = mat(i,:)
                    mat(i,:)    = 
     +                          ( carr(:dim2) + img * mat(j,:) ) * rfac1
                    mat(j,:)    =
     +                          ( carr(:dim2) - img * mat(j,:) ) * rfac2
                  END IF
                  IF( iand(imode,2).ne.0) THEN
                    carr(:dim1) = mat(:,i)
                    mat(:,i)    = 
     +                          ( carr(:dim1) - img * mat(:,j) ) * rfac1
                    mat(:,j)    = 
     +                          ( carr(:dim1) + img * mat(:,j) ) * rfac2
                  END IF
                ELSE IF( m .eq. 0 .and. ifac .eq. -1 ) THEN
                  IF( iand(imode,1) .ne. 0 ) THEN
                    mat(i,:) =  img * mat(i,:)
                  END IF
                  IF( iand(imode,2) .ne. 0) THEN
                    mat(:,i) = -img * mat(:,i)
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO

      END SUBROUTINE desymmetrize

      ! bra_trafo1 rotates cprod at ikpt0(<=> not irreducible k-point) to cprod at ikpt1 (bkp(ikpt0)), which is the
      ! symmetrie equivalent one
      ! isym maps ikpt0 on ikpt1


      SUBROUTINE bra_trafo2 (
     &          vecout,vecin,dim,nobd,nbands,ikpt0,ikpt1,iop,nop,nsym,
     &          mrot,invtab,map,nkpt,nkpt3,bk,tau,bmat,maxlcutm,ntype,
     &          neq,lcutm,natd,taual,nindxm,maxindxm,nw,nwd,lmaxd,d_wgn,
     &          ngptmall,gptm,ngptm,pgptm,nbasp,invsat,
     &          invsatnr,phase)

      !  ikpt0  ::  parent of ikpt1
      !  iop maps ikpt0 on ikpt1

      USE m_constants     
      USE m_dwigner
      USE m_util

      IMPLICIT NONE

c     - scalars -
      INTEGER,INTENT(IN)      ::  ikpt0,ikpt1,iop,dim,nobd,nbands
      INTEGER,INTENT(IN)      ::  nop,nkpt,maxlcutm,ntype,natd,maxindxm,
     &                            nw,nwd,nbasp,nsym
      INTEGER,INTENT(IN)      ::  lmaxd
c     - arrays -
      INTEGER,INTENT(IN)      ::  mrot(3,3,nop),invtab(nop),neq(ntype),
     &                            lcutm(ntype),nindxm(0:maxlcutm,ntype)
      INTEGER,INTENT(IN)      ::  ngptmall
      INTEGER,INTENT(IN)      ::  gptm(3,ngptmall),
     &                            ngptm(nkpt),
     &                            pgptm(maxval(ngptm),nkpt)
      INTEGER,INTENT(IN)      ::  invsat(natd),invsatnr(natd)
      INTEGER,INTENT(IN)      ::  nkpt3(3)
      INTEGER,INTENT(IN)      ::  map(natd,nsym)
      REAL,INTENT(IN)         ::  bk(3,nkpt),tau(3,nop),bmat(3,3),
     &                            taual(3,natd)
#if ( defined(CPP_INVERSION) )
      REAL,INTENT(IN)         ::  vecin(dim,nobd,nbands)
      REAL,INTENT(OUT)        ::  vecout(dim,nobd,nbands)
#else
      COMPLEX,INTENT(IN)      ::  vecin(dim,nobd,nbands)
      COMPLEX,INTENT(OUT)     ::  vecout(dim,nobd,nbands)
#endif
      COMPLEX,INTENT(IN)      ::  d_wgn(-lmaxd:lmaxd,-lmaxd:lmaxd,
     +                                       0:lmaxd,nsym)
      COMPLEX,INTENT(OUT)     ::  phase(nobd,nbands)

c          - local -

c     - scalars -  
      INTEGER                 ::  nrkpt,rcent,itype,ieq,ic,l,n,i,j,nn,
     &                            i1,i2,j1,j2,m1,m2,ok
      INTEGER                 ::  igptm,igptm2,igptp,icent1,iiatom,iiop,
     +                            inviop
      COMPLEX                 ::  cexp,cdum
      COMPLEX,PARAMETER       ::  img=(0d0,1d0)
c     - arrays - 

      INTEGER                 ::  rrot(3,3),invrot(3,3)
      INTEGER                 ::  pnt(maxindxm,0:maxlcutm,natd)
      INTEGER                 ::  g(3),g1(3)
      REAL                    ::  rkpt(3),rkpthlp(3),rtaual(3),trans(3)
      REAL                    ::  arg
      COMPLEX                 ::  dwgn(-maxlcutm:maxlcutm,
     &                                 -maxlcutm:maxlcutm,0:maxlcutm)
!       COMPLEX                 ::  vecin1(dim,nobd,nbands),vecout1(dim,nobd,nbands)
      COMPLEX, ALLOCATABLE    ::  vecin1(:,:,:),vecout1(:,:,:)

      ALLOCATE ( vecin1( dim,nobd,nbands), 
     +           vecout1(dim,nobd,nbands), stat=ok )
      IF ( ok /= 0 ) 
     +             STOP 'bra_trafo2: error allocating vecin1 or vecout1'
      vecin1 = 0 ; vecout1 = 0

      IF( maxlcutm .gt. lmaxd ) STOP 'bra_trafo2: maxlcutm > lmaxd'   ! very improbable case

c     transform back to unsymmetrized product basis in case of inversion symmetry
      vecin1 = vecin
#if ( defined(CPP_INVERSION) )
      DO i=1,nbands
        DO j=1,nobd
          CALL desymmetrize(vecin1(:nbasp,j,i),nbasp,1,1,
     &                      ntype,ntype,neq,lcutm,maxlcutm,
     &                      nindxm,natd,invsat,invsatnr)
        END DO
      END DO
#endif

    
      IF( iop .le. nop ) THEN
        inviop = invtab(iop)
        rrot   = transpose( mrot(:,:,invtab(iop)) )
        invrot = mrot(:,:,invtab(iop))
        trans  = tau(:,iop)

        dwgn (-maxlcutm:maxlcutm,-maxlcutm:maxlcutm,0:maxlcutm) 
     &= d_wgn(-maxlcutm:maxlcutm,-maxlcutm:maxlcutm,0:maxlcutm,inviop)


      ELSE
        iiop   = iop -nop
        inviop = invtab(iiop)+nop
        rrot   = -transpose( mrot(:,:,invtab(iiop)) )
        invrot = mrot(:,:,invtab(iiop))
        trans  = tau(:,iiop)

        dwgn (-maxlcutm:maxlcutm,-maxlcutm:maxlcutm,0:maxlcutm)
     &= conjg( d_wgn(-maxlcutm:maxlcutm,
     +               -maxlcutm:maxlcutm,
     +                0:maxlcutm,inviop))

      END IF


      rkpt    = matmul(rrot,bk(:,ikpt0))
      rkpthlp = rkpt
      rkpt    = modulo1(rkpt,nkpt3)
      g       = nint(rkpthlp-rkpt)

#ifdef CPP_DEBUG
       !test
      DO i=1,nkpt
        IF ( maxval( abs(rkpt - bk(:,i)) ) .le. 1E-06  ) THEN
          nrkpt = i
          EXIT
        END IF
      END DO
      IF( nrkpt .ne. ikpt1 ) STOP 'bra_trafo2: rotation failed'
#endif

c     Define pointer to first mixed-basis functions (with m = -l)
      i  = 0
      ic = 0 
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          ic = ic + 1
          DO l = 0,lcutm(itype)
            DO n = 1,nindxm(l,itype)
              i           = i + 1
              pnt(n,l,ic) = i
            END DO
            i = i + nindxm(l,itype) * 2*l
          END DO
        END DO
      END DO

c     Multiplication
      ! MT
      cexp   = exp(img*tpi_const*dot_product( bk(:,ikpt1)+g,trans(:) ) )
      ic     = 0
      iiatom = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          ic   = ic + 1 

          rcent = map(ic,iop)

          cdum = cexp *exp(-img*tpi_const*dot_product(g,taual(:,rcent)))

          DO l = 0,lcutm(itype)
            nn = nindxm(l,itype)
            DO n = 1,nn

              i1 = pnt(n,l,ic)
              i2 = i1 + nn * 2*l
              j1 = pnt(n,l,rcent)
              j2 = j1 + nn * 2*l

              DO i=1,nbands
                DO j=1,nobd
                  vecout1(i1:i2:nn,j,i) = 
     +             cdum * matmul(vecin1(j1:j2:nn,j,i),dwgn(-l:l,-l:l,l))
                END DO
              END DO


            END DO
          END DO
        END DO
        iiatom = iiatom + neq(itype)
      END DO

      ! PW
      DO igptm = 1,ngptm(ikpt0)
        igptp  = pgptm(igptm,ikpt0)
        g1     = matmul(rrot,gptm(:,igptp)) + g
        igptm2 = 0
        DO i=1,ngptm(ikpt1)
          IF ( maxval(abs(g1-gptm(:,pgptm(i,ikpt1)))) .le. 1E-06 ) THEN
            igptm2 = i
            EXIT
          END IF
        END DO 
        IF(igptm2 .eq. 0) THEN
          WRITE(*,*) ikpt0,ikpt1,g1
          WRITE(*,*) ngptm(ikpt0),ngptm(ikpt1)
          WRITE(*,*)
          WRITE(*,*) gptm(:,igptp)
          WRITE(*,*) g
          STOP 'bra_trafo2: G-point not found in G-point set.'
        END IF
        cdum = exp(img*tpi_const*dot_product(bk(:,ikpt1)+g1,trans(:)))

        vecout1(nbasp+igptm,:,:)= cdum * vecin1(nbasp+igptm2,:,:)
      END DO

      DEALLOCATE ( vecin1 )

#if ( defined(CPP_INVERSION) )
      DO i=1,nbands
        DO j=1,nobd

          CALL symmetrize(vecout1(:,j,i),dim,1,1,.false.,
     &                    ntype,ntype,neq,lcutm,maxlcutm,
     &                    nindxm,natd,invsat,invsatnr)

          CALL commonphase(phase(j,i),vecout1(:,j,i),dim)
          vecout1(:,j,i) = vecout1(:,j,i) / phase(j,i)
          IF(any(abs(aimag(vecout1(:,j,i))).gt.1d-8)) THEN
            WRITE(*,*) vecout1(:,j,i)
            STOP 'bra_trafo2: Residual imaginary part.'
          END IF

        END DO
      END DO 
#else
      phase = (1d0,0d0)
#endif

      vecout = vecout1

      DEALLOCATE ( vecout1 )

      END SUBROUTINE bra_trafo2


c     This routine is not very fast at the moment.
      SUBROUTINE matrixtrafo( 
     X       matout,matin,ikpt0,isym,lsymmetrize,nwd,nw,ntype, 
     >       nkpt,nkpt3,bk,nop,nsym,neq,mrot,tau,invtab,invsat,invsatnr,
     >       map,bmat,maxindxm,nindxm,nbasm,ngptm,ngptmall,nbasp,
     >       pgptm,gptm,lcutm,maxlcutm,natd,taual)


      USE m_wrapper
      USE m_dwigner
      USE m_util      ,ONLY: modulo1
      USE m_constants  
      IMPLICIT NONE

      ! - scalars -
      INTEGER , INTENT(IN)  ::  ikpt0,isym
      INTEGER , INTENT(IN)  ::  ntype,nkpt,nop,natd,nwd,nw,nbasp,nsym
      INTEGER , INTENT(IN)  ::  maxlcutm,maxindxm
      LOGICAL , INTENT(IN)  ::  lsymmetrize

      ! - arrays -
      INTEGER , INTENT(IN)  ::  mrot(3,3,nop),invtab(nop)
      INTEGER , INTENT(IN)  ::  invsat(natd),invsatnr(natd)
      INTEGER , INTENT(IN)  ::  neq(ntype)
      INTEGER , INTENT(IN)  ::  nbasm(nkpt),ngptm(nkpt)
      INTEGER , INTENT(IN)  ::  nindxm(0:maxlcutm,ntype)
      INTEGER , INTENT(IN)  ::  pgptm(maxval(ngptm),nkpt)
      INTEGER , INTENT(IN)  ::  ngptmall
      INTEGER , INTENT(IN)  ::  gptm(3,ngptmall)
      INTEGER , INTENT(IN)  ::  lcutm(ntype)
      INTEGER , INTENT(IN)  ::  nkpt3(3)
      INTEGER , INTENT(IN)  ::  map(natd,nsym)

      REAL    , INTENT(IN)  ::  bk(3,nkpt),bmat(3,3)
      REAL    , INTENT(IN)  ::  taual(3,natd),tau(3,nop)

      COMPLEX , INTENT(IN)  ::  matin(nbasm(ikpt0),nbasm(ikpt0))
      COMPLEX , INTENT(OUT) ::  matout(nbasm(ikpt0),nbasm(ikpt0))
      ! - local scalars -
      INTEGER               ::  iatom,iatom1,iiatom,itype,ieq,ieq1,ic,
     +                          l,n,i,nn,i1,i2,j1,j2
      INTEGER               ::  igptm,igptm1,igptm2,igptp,igptp1,
     +                          ikpt1,isymi,iisym
      INTEGER               ::  m1,m2

      REAL                  ::  time1,time2
      COMPLEX               ::  cexp,cdum 
      COMPLEX , PARAMETER   ::  img=(0d0,1d0)

       ! - local arrays -
      INTEGER               ::  pnt(maxindxm,0:maxlcutm,natd),
     +                          g(3),g1(3),iarr(ngptm(ikpt0))
      INTEGER               ::  rot(3,3),invrot(3,3),
     +                          rrot(3,3),invrrot(3,3)

      REAL                  ::  rkpt(3),rkpthlp(3),rtaual(3)
      REAL                  ::  trans(3)

      COMPLEX               ::  matin1(nbasm(ikpt0),nbasm(ikpt0))
      COMPLEX               ::  matout1(nbasm(ikpt0),nbasm(ikpt0))
      COMPLEX               ::  dwgn   (-maxlcutm:maxlcutm,
     +                                  -maxlcutm:maxlcutm,
     +                                          0:maxlcutm)
      COMPLEX               ::  dwgninv(-maxlcutm:maxlcutm,
     +                                  -maxlcutm:maxlcutm,
     +                                          0:maxlcutm)
      COMPLEX               ::  carr(ngptm(ikpt0))

      CALL cpu_time(time1)

   
c     Transform back to unsymmetrized product basis in case of inversion symmetry.
      matin1 = matin      
# ifdef CPP_INVERSION
      IF(lsymmetrize) THEN
        CALL desymmetrize(matin1,nbasm(ikpt0),nbasm(ikpt0),3,
     +                    ntype,ntype,neq,lcutm,maxlcutm,
     +                    nindxm,natd,invsat,invsatnr)
      END IF
# endif

      IF( isym .le. nop ) THEN
        iisym   = isym
        rot     = mrot(:,:,iisym)
        invrot  = mrot(:,:,invtab(iisym))
        rrot    = transpose( mrot(:,:,invtab(iisym)) )
        invrrot = transpose( mrot(:,:,iisym) )
        rkpt    = matmul(rrot,bk(:,ikpt0))
        rkpthlp = modulo1(rkpt,nkpt3)
        g       = nint(rkpt - rkpthlp)


        CALL d_wigner(invrot,bmat,maxlcutm,dwgn(:,:,1:maxlcutm))
        dwgn(0,0,0) = 1

        DO l=0,maxlcutm
          dwgn(:,:,l) = transpose(dwgn(:,:,l))
        END DO 
      ELSE
        iisym   = isym - nop
        rot     = mrot(:,:,iisym)
        invrot  = mrot(:,:,invtab(iisym))
        rrot    = -transpose( mrot(:,:,invtab(iisym)) )
        invrrot = -transpose( mrot(:,:,iisym) )
        rkpt    = matmul(rrot,bk(:,ikpt0))
        rkpthlp = modulo1(rkpt,nkpt3)
        g       = nint(rkpt - rkpthlp)
        matin1  = conjg(matin1)

        CALL d_wigner(invrot,bmat,maxlcutm,dwgn(:,:,1:maxlcutm))
        dwgn(0,0,0) = 1

        DO l=0,maxlcutm
          dwgn(:,:,l) = transpose(dwgn(:,:,l))
        END DO

        DO l = 0,maxlcutm
          DO m1 = -l,l
            DO m2 = -l,-1
              cdum           = dwgn(m1, m2,l)
              dwgn(m1, m2,l) = dwgn(m1,-m2,l) * (-1)**m2
              dwgn(m1,-m2,l) = cdum           * (-1)**m2
            END DO
          END DO
        END DO

      END IF
      ! determine number of rotated k-point bk(:,ikpt) -> ikpt1
      ! 
      DO i=1,nkpt
        IF ( maxval( abs(rkpthlp - bk(:,i)) ) .le. 1E-06  ) THEN
          ikpt1 = i
          EXIT
        END IF
      END DO



      DO l = 0,maxlcutm
        dwgninv(-l:l,-l:l,l) = conjg(transpose(dwgn(-l:l,-l:l,l)))
      END DO

c     Define pointer to first mixed-basis functions (with m = -l)
      i  = 0
      ic = 0 
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          ic = ic + 1
          DO l = 0,lcutm(itype)
            DO n = 1,nindxm(l,itype)
              i           = i + 1
              pnt(n,l,ic) = i
            END DO
            i = i + nindxm(l,itype) * 2*l
          END DO
        END DO
      END DO


c     Right-multiplication
      ! MT
      cexp = exp(img*tpi_const* dot_product(bk(:,ikpt1)+g,tau(:,iisym)))
      iatom  = 0
      iiatom = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom  = iatom + 1
          cdum   =cexp*exp(-img*tpi_const*dot_product(g,taual(:,iatom)))

!           rtaual = matmul(invrot,taual(:,iatom)) + tau(:,invtab(iisym)) 
!           iatom1 = 0
!           DO ieq1 = 1,neq(itype)
!             IF( all(abs(modulo(rtaual-taual(:,iiatom+ieq1)+1d-12,1d0))
!      &                                               .lt. 1d-10) ) THEN ! The 1d-12 is a dirty fix.
!               iatom1 = iiatom + ieq1
!             END IF
!           END DO
!           IF( iatom1 .eq. 0 ) STOP 'matrixtrafo atom not found'          

          iatom1 = map(iatom,invtab(iisym))
          DO l = 0,lcutm(itype)
            nn = nindxm(l,itype)
            DO n = 1,nn

              i1 = pnt(n,l,iatom)
              i2 = i1 + nn * 2*l
              j1 = pnt(n,l,iatom1)
              j2 = j1 + nn * 2*l

              matout1(:,i1:i2:nn) = cdum * matmat( matin1(:,j1:j2:nn), 
     +                                             dwgn(-l:l,-l:l,l) )

            END DO
          END DO
        END DO
        iiatom = iiatom + neq(itype)
      END DO


      ! PW
      DO igptm = 1,ngptm(ikpt1)
        igptp  = pgptm(igptm,ikpt1)
        g1     = matmul(invrrot,gptm(:,igptp)-g)
        igptm2 = 0
        DO i=1,ngptm(ikpt0)
          IF ( maxval( abs(g1-gptm(:,pgptm(i,ikpt0) )) ) 
     &                                               .le. 1E-06  ) THEN
            igptm2 = i
            EXIT
          END IF
        END DO 
!         igptm2 = pntgptm(g1(1),g1(2),g1(3),ikpt0)
        IF(igptm2.eq.0)
     +             STOP 'matrixtrafo: G point not found in G-point set.'

        cdum = exp(img * tpi_const *
     +              dot_product(bk(:,ikpt1)+gptm(:,igptp),tau(:,iisym)))

        matout1(:,nbasp+igptm) = cdum * matin1(:,nbasp+igptm2) 

      END DO


c     Left-multiplication
      ! MT
      matin1 = matout1
      cexp   = conjg(cexp)
      iatom  = 0
      iiatom = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom = iatom + 1
          cdum = cexp * exp(img*tpi_const*dot_product(g,taual(:,iatom)))

!           rtaual = matmul(invrot,taual(:,iatom)) + tau(:,invtab(iisym)) 
!           iatom1 = 0
!           DO ieq1 = 1,neq(itype)
!             IF( all(abs(modulo(rtaual-taual(:,iiatom+ieq1)+1d-12,1d0))
!      &                                               .lt. 1d-10) ) THEN ! The 1d-12 is a dirty fix.
!               iatom1 = iiatom + ieq1
!             END IF
!           END DO
!           IF( iatom1 .eq. 0 ) STOP 'matrixtrafo atom not found'  
          iatom1 = map(iatom,invtab(iisym))     

          DO l = 0,lcutm(itype)
            nn = nindxm(l,itype)
            DO n = 1,nn

              i1 = pnt(n,l,iatom)
              i2 = i1 + nn * 2*l
              j1 = pnt(n,l,iatom1)
              j2 = j1 + nn * 2*l

              matout1(i1:i2:nn,:) = cdum * 
     +                   matmat(dwgninv(-l:l,-l:l,l),matin1(j1:j2:nn,:))

            END DO
          END DO
        END DO
        iiatom = iiatom + neq(itype)
      END DO


      ! PW
      DO igptm = 1,ngptm(ikpt1)
        igptp  = pgptm(igptm,ikpt1)
        g1     = matmul(invrrot,gptm(:,igptp)-g)
        igptm2 = 0
        DO i=1,ngptm(ikpt0)
          IF ( maxval( abs(g1-gptm(:,pgptm(i,ikpt0) )) ) 
     &                                               .le. 1E-06  ) THEN
            igptm2 = i
            EXIT
          END IF
        END DO 
!         igptm2 = pntgptm(g1(1),g1(2),g1(3),ikpt0)
        IF(igptm2.eq.0) 
     +             STOP 'matrixtrafo: G point not found in G-point set.'
        iarr(igptm) = igptm2
        carr(igptm) = exp(-img * tpi_const * 
     +              dot_product(bk(:,ikpt1)+gptm(:,igptp),tau(:,iisym)))
c        cdum  = exp(-img * 2*pi * dot_product(bk(:,ikpt1)+gptm(:,igptp),tau(:,isym)))
c        matout1(nbasp+igptm,:) = cdum * matin1(nbasp+igptm2,:)
      END DO
      DO i2 = 1,nbasm(ikpt1)
        DO i1 = 1,ngptm(ikpt1)
          matout1(nbasp+i1,i2) = carr(i1) * matin1(nbasp+iarr(i1),i2)
        END DO
      END DO

      ! If inversion symmetry is applicable, symmetrize to make the values real.
# ifdef CPP_INVERSION
      IF(lsymmetrize) THEN
        CALL symmetrize(matout1,nbasm(ikpt0),nbasm(ikpt0),3,.false.,
     +                  ntype,ntype,neq,lcutm,maxlcutm,
     +                  nindxm,natd,invsat,invsatnr)
      END IF
# endif
      matout = matout1
      END SUBROUTINE matrixtrafo


      SUBROUTINE matrixtrafo1(
     X                 matout,matin,ikpt0,isym,lsymmetrize,nwd,nw,ntype,
     >                 nkpt,nkpt3,bk,nop,nsym,neq,mrot,d_wgn,lmaxd,tau,
     >                 invtab,invsat,invsatnr,map,bmat,maxindxm,nindxm,
     >                 nbasm,ngptm,ngptmall,nbasp,pgptm,gptm,lcutm,
     >                 maxlcutm,natd,taual)


      USE m_wrapper
      USE m_dwigner
      USE m_util      ,ONLY: modulo1
      USE m_constants 

      IMPLICIT NONE

      ! - scalars -
      INTEGER , INTENT(IN)  ::  ikpt0,isym
      INTEGER , INTENT(IN)  ::  ntype,nkpt,nop,natd,nwd,nw,nbasp,nsym
      INTEGER , INTENT(IN)  ::  maxlcutm,maxindxm
      INTEGER , INTENT(IN)  ::  lmaxd
      LOGICAL , INTENT(IN)  ::  lsymmetrize


      ! - arrays -
      INTEGER , INTENT(IN)  ::  mrot(3,3,nop),invtab(nop)
      INTEGER , INTENT(IN)  ::  invsat(natd),invsatnr(natd)
      INTEGER , INTENT(IN)  ::  neq(ntype)
      INTEGER , INTENT(IN)  ::  nbasm(nkpt),ngptm(nkpt)
      INTEGER , INTENT(IN)  ::  nindxm(0:maxlcutm,ntype)
      INTEGER , INTENT(IN)  ::  pgptm(maxval(ngptm),nkpt)
      INTEGER , INTENT(IN)  ::  ngptmall
      INTEGER , INTENT(IN)  ::  gptm(3,ngptmall)
      INTEGER , INTENT(IN)  ::  lcutm(ntype)
      INTEGER , INTENT(IN)  ::  nkpt3(3)
      INTEGER , INTENT(IN)  ::  map(natd,nsym)

      REAL    , INTENT(IN)  ::  bk(3,nkpt),bmat(3,3)
      REAL    , INTENT(IN)  ::  taual(3,natd),tau(3,nop)

      COMPLEX , INTENT(IN)  ::  d_wgn(-lmaxd:lmaxd,-lmaxd:lmaxd,
     +                                     0:lmaxd,nsym)
      COMPLEX , INTENT(IN)  ::  matin(nbasm(ikpt0),nbasm(ikpt0))
      COMPLEX , INTENT(OUT) ::  matout(nbasm(ikpt0),nbasm(ikpt0))
      ! - local scalars -
      INTEGER               ::  iatom,iatom1,iiatom,itype,ieq,ieq1,ic,l,
     +                          n,i,nn,i1,i2,j1,j2
      INTEGER               ::  igptm,igptm1,igptm2,igptp,igptp1,ikpt1,
     +                          isymi,iisym
      INTEGER               ::  m1,m2

      REAL                  ::  time1,time2
      COMPLEX               ::  cexp,cdum 
      COMPLEX , PARAMETER   ::  img=(0d0,1d0)

       ! - local arrays -
      INTEGER               ::  pnt(maxindxm,0:maxlcutm,natd),g(3),
     +                          g1(3),iarr(ngptm(ikpt0))
      INTEGER               ::  rrot(3,3)

      REAL                  ::  rkpt(3),rkpthlp(3),rtaual(3)
      REAL                  ::  trans(3)

      COMPLEX               ::  matin1(nbasm(ikpt0),nbasm(ikpt0))
      COMPLEX               ::  matout1(nbasm(ikpt0),nbasm(ikpt0))
      COMPLEX               ::  dwgn   (-maxlcutm:maxlcutm,
     +                                  -maxlcutm:maxlcutm,
     +                                          0:maxlcutm)
      COMPLEX               ::  dwgninv(-maxlcutm:maxlcutm,
     +                                  -maxlcutm:maxlcutm,
     +                                          0:maxlcutm)
      COMPLEX               ::  carr(ngptm(ikpt0))


      IF( maxlcutm .gt. lmaxd )STOP 'matrixtrafo1: maxlcutm .gt. lmaxd'
      CALL cpu_time(time1)

      
c     Transform back to unsymmetrized product basis in case of inversion symmetry.
      matin1 = matin
# ifdef CPP_INVERSION
      IF(lsymmetrize) THEN
        CALL desymmetrize(matin1,nbasm(ikpt0),nbasm(ikpt0),3,
     &                    ntype,ntype,neq,lcutm,maxlcutm,
     &                    nindxm,natd,invsat,invsatnr)
      END IF
# endif

      IF( isym .le. nop ) THEN
        iisym   = isym
        rrot    = transpose( mrot(:,:,invtab(iisym)) )
      ELSE
        iisym   = isym - nop
        rrot    = -transpose( mrot(:,:,invtab(iisym)) )
      END IF

      DO l = 0,maxlcutm
        dwgn(-maxlcutm:maxlcutm,-maxlcutm:maxlcutm,l) =
     &    transpose(d_wgn(-maxlcutm:maxlcutm,-maxlcutm:maxlcutm,l,isym))
      END DO

      rkpt    = matmul(rrot,bk(:,ikpt0))
      rkpthlp = modulo1(rkpt,nkpt3)
      g       = nint(rkpt - rkpthlp)

      ! determine number of rotated k-point bk(:,ikpt) -> ikpt1
      DO i=1,nkpt
        IF ( maxval( abs(rkpthlp - bk(:,i)) ) .le. 1E-06  ) THEN
          ikpt1 = i
          EXIT
        END IF
      END DO

c     Define pointer to first mixed-basis functions (with m = -l)
      i  = 0
      ic = 0 
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          ic = ic + 1
          DO l = 0,lcutm(itype)
            DO n = 1,nindxm(l,itype)
              i           = i + 1
              pnt(n,l,ic) = i
            END DO
            i = i + nindxm(l,itype) * 2*l
          END DO
        END DO
      END DO


c     Right-multiplication
      ! MT
      cexp = exp(-img*tpi_const*dot_product(bk(:,ikpt1)+g,tau(:,iisym)))
      iatom  = 0
      iiatom = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom  = iatom + 1

          iatom1 = map(iatom,iisym)
          cdum   = cexp * 
     +              exp(img *tpi_const*dot_product(g,taual(:,iatom1)))


          DO l = 0,lcutm(itype)
            nn = nindxm(l,itype)
            DO n = 1,nn

              i1 = pnt(n,l,iatom)
              i2 = i1 + nn * 2*l
              j1 = pnt(n,l,iatom1)
              j2 = j1 + nn * 2*l

              matout1(:,i1:i2:nn) = cdum * matmat( matin1(:,j1:j2:nn) ,
     +                                             dwgn(-l:l,-l:l,l))

            END DO
          END DO
        END DO
        iiatom = iiatom + neq(itype)
      END DO

      ! PW
      DO igptm = 1,ngptm(ikpt0)
        igptp  = pgptm(igptm,ikpt0)
        g1     = matmul(rrot,gptm(:,igptp)) + g
        igptm1 = 0
        DO i=1,ngptm(ikpt1)
          IF( maxval( abs(g1-gptm(:,pgptm(i,ikpt1) )) ) <= 1E-06  ) THEN
            igptm1 = i
            igptp1 = pgptm(i,ikpt1)
            EXIT
          END IF
        END DO
        IF(igptm1.eq.0) 
     +            STOP 'matrixtrafo1: G point not found in G-point set.'

        cdum = exp(-img * tpi_const * 
     +             dot_product(bk(:,ikpt1)+gptm(:,igptp1),tau(:,iisym)))

        matout1(:,nbasp+igptm) = cdum * matin1(:,nbasp+igptm1)

      END DO 



c     Left-multiplication
      ! MT

      DO l = 0,maxlcutm
        dwgninv(-l:l,-l:l,l) = conjg(transpose(dwgn(-l:l,-l:l,l)))
      END DO

      matin1 = matout1
      cexp   = conjg(cexp)

      iatom  = 0
      iiatom = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom = iatom + 1


          iatom1 = map(iatom,iisym)
          cdum = cexp*exp(-img*tpi_const*dot_product(g,taual(:,iatom1)))

          DO l = 0,lcutm(itype)
            nn = nindxm(l,itype)
            DO n = 1,nn

              i1 = pnt(n,l,iatom)
              i2 = i1 + nn * 2*l
              j1 = pnt(n,l,iatom1)
              j2 = j1 + nn * 2*l

              matout1(i1:i2:nn,:) = cdum * matmat(dwgninv(-l:l,-l:l,l),
     +                                            matin1(j1:j2:nn,:))

            END DO
          END DO
        END DO
        iiatom = iiatom + neq(itype)
      END DO


      ! PW
       DO igptm = 1,ngptm(ikpt0)
        igptp  = pgptm(igptm,ikpt0)
        g1     = matmul(rrot,gptm(:,igptp)) + g
        igptm1 = 0
        DO i=1,ngptm(ikpt1)
          IF( maxval( abs(g1-gptm(:,pgptm(i,ikpt1) )) ) <= 1E-06 ) THEN
            igptm1 = i
            igptp1 = pgptm(i,ikpt1)
            EXIT
          END IF
        END DO
        IF(igptm1.eq.0) 
     +            STOP 'matrixtrafo1: G point not found in G-point set.'
        iarr(igptm) = igptm1
        carr(igptm) = exp(img * tpi_const * 
     +             dot_product(bk(:,ikpt1)+gptm(:,igptp1),tau(:,iisym)))
      END DO 
      DO i2 = 1,nbasm(ikpt0)
        DO i1 = 1,ngptm(ikpt0)
          matout1(nbasp+i1,i2) = carr(i1) * matin1(nbasp+iarr(i1),i2)
        END DO
      END DO


      ! If inversion symmetry is applicable, symmetrize to make the values real.
# ifdef CPP_INVERSION
      IF(lsymmetrize) THEN
        CALL symmetrize(matout1,nbasm(ikpt0),nbasm(ikpt0),3,.false.,
     &                  ntype,ntype,neq,lcutm,maxlcutm,
     &                  nindxm,natd,invsat,invsatnr)
      END IF
# endif

      IF( isym .le. nop ) THEN
        matout = matout1
      ELSE
        matout = conjg(matout1)
      END IF

      END SUBROUTINE matrixtrafo1


      SUBROUTINE ket_trafo(
     X        vecout,vecin,ikpt0,isym,lreal,lsymmetrize,
     >        nwd,nw,ntype,nkpt,nkpt3,bk,nop,nsym,neq,mrot,tau,invtab,
     >        invsat,invsatnr,map,bmat,maxindxm,nindxm,nbasm,
     >        ngptm,ngptmall,nbasp,pgptm,gptm,lcutm,maxlcutm,natd,taual)

      USE m_constants 
      USE m_util      ,ONLY: modulo1
      USE m_dwigner

      IMPLICIT NONE

      ! -scalars -
      INTEGER , INTENT(IN)  ::  ikpt0,isym
      INTEGER , INTENT(IN)  ::  ntype,nkpt,nop,natd,nwd,nw,nbasp,nsym
      INTEGER , INTENT(IN)  ::  maxlcutm,maxindxm
      LOGICAL , INTENT(IN)  ::  lreal,lsymmetrize

      ! - arrays -

      INTEGER , INTENT(IN)  ::  mrot(3,3,nop),invtab(nop)
      INTEGER , INTENT(IN)  ::  invsat(natd),invsatnr(natd)
      INTEGER , INTENT(IN)  ::  neq(ntype)
      INTEGER , INTENT(IN)  ::  nbasm(nkpt),ngptm(nkpt)
      INTEGER , INTENT(IN)  ::  nindxm(0:maxlcutm,ntype)
      INTEGER , INTENT(IN)  ::  pgptm(maxval(ngptm),nkpt)
      INTEGER , INTENT(IN)  ::  ngptmall
      INTEGER , INTENT(IN)  ::  gptm(3,ngptmall)
      INTEGER , INTENT(IN)  ::  lcutm(ntype)
      INTEGER , INTENT(IN)  ::  nkpt3(3)
      INTEGER , INTENT(IN)  ::  map(natd,nsym)

      REAL    , INTENT(IN)  ::  bk(3,nkpt),bmat(3,3)
      REAL    , INTENT(IN)  ::  taual(3,natd),tau(3,nop)

      COMPLEX , INTENT(IN)  ::  vecin(nbasm(ikpt0))
      COMPLEX , INTENT(OUT) ::  vecout(nbasm(ikpt0))

      ! -local scalars -
      INTEGER               ::  iatom,iatom1,iiatom,itype,ieq,ieq1,ic,l,
     +                          n,i,nn,i1,i2,j1,j2
      INTEGER               ::  igptm,igptm1,igptm2,igptp,igptp1,ikpt1,
     +                          isymi,iisym
      INTEGER               ::  m1,m2

      REAL                  ::  time1,time2
      COMPLEX               ::  cexp,cdum 
      COMPLEX , PARAMETER   ::  img=(0d0,1d0)

      ! - local arrays -
      INTEGER               ::  pnt(maxindxm,0:maxlcutm,natd),g(3),g1(3)
      INTEGER               ::  rot(3,3),invrot(3,3),
     +                          rrot(3,3),invrrot(3,3)

      REAL                  ::  rkpt(3),rkpthlp(3),rtaual(3)

      COMPLEX               ::  vecin1 (nbasm(ikpt0))
      COMPLEX               ::  vecout1(nbasm(ikpt0))
      COMPLEX               ::  dwgn   (-maxlcutm:maxlcutm,
     +                                  -maxlcutm:maxlcutm,
     +                                          0:maxlcutm)
      COMPLEX               ::  dwgninv(-maxlcutm:maxlcutm,
     +                                  -maxlcutm:maxlcutm,
     +                                          0:maxlcutm)


  
c     Transform back to unsymmetrized product basis in case of inversion symmetry.
      vecin1 = vecin!(:nbasp)
# ifdef CPP_INVERSION
      IF(lsymmetrize) THEN
        CALL desymmetrize(vecin1(:nbasp),1,nbasp,2,
     &                    ntype,ntype,neq,lcutm,maxlcutm,
     &                    nindxm,natd,invsat,invsatnr)
      END IF
# endif

      IF( isym .le. nop ) THEN
        iisym   = isym
        rot     = mrot(:,:,iisym)
        invrot  = mrot(:,:,invtab(iisym))
        rrot    = transpose( mrot(:,:,invtab(iisym)) )
        invrrot = transpose( mrot(:,:,iisym) )
        rkpt    = matmul(rrot,bk(:,ikpt0))
        rkpthlp = modulo1(rkpt,nkpt3)
        g       = nint(rkpt - rkpthlp)

        CALL d_wigner(invrot,bmat,maxlcutm,dwgn(:,:,1:maxlcutm))
        dwgn(0,0,0) = 1

        DO l=0,maxlcutm
          dwgn(:,:,l) = transpose(dwgn(:,:,l))
        END DO
      ELSE
        iisym   = isym - nop
        rot     = mrot(:,:,iisym)
        rrot    = -transpose( mrot(:,:,invtab(iisym)) )
        invrot  = mrot(:,:,invtab(iisym))
        invrrot = -transpose( mrot(:,:,iisym) )
        rkpt    = matmul(rrot,bk(:,ikpt0))
        rkpthlp = modulo1(rkpt,nkpt3)
        g       = nint(rkpt - rkpthlp)
        vecin1  = conjg(vecin1)

        CALL d_wigner(invrot,bmat,maxlcutm,dwgn(:,:,1:maxlcutm))
        dwgn(0,0,0) = 1

        DO l=0,maxlcutm
          dwgn(:,:,l) = transpose(dwgn(:,:,l))
        END DO

        DO l = 0,maxlcutm
          DO m1 = -l,l
            DO m2 = -l,-1
              cdum           = dwgn(m1, m2,l)
              dwgn(m1, m2,l) = dwgn(m1,-m2,l) * (-1)**m2
              dwgn(m1,-m2,l) = cdum           * (-1)**m2
            END DO
          END DO
        END DO
      END IF

      ! determine number of rotated k-point bk(:,ikpt) -> ikpt1
      ! 
      DO i=1,nkpt
        IF ( maxval( abs(rkpthlp - bk(:,i)) ) .le. 1E-06  ) THEN
          ikpt1 = i
          EXIT
        END IF
      END DO



!       DO l = 0,maxlcutm
!         dwgninv(-l:l,-l:l,l) = conjg(transpose(dwgn(-l:l,-l:l,l)))
!       END DO

c     Define pointer to first mixed-basis functions (with m = -l)
      i  = 0
      ic = 0 
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          ic = ic + 1
          DO l = 0,lcutm(itype)
            DO n = 1,nindxm(l,itype)
              i           = i + 1
              pnt(n,l,ic) = i
            END DO
            i = i + nindxm(l,itype) * 2*l
          END DO
        END DO
      END DO

c     Multiplication
      ! MT
      cexp   =exp(img*tpi_const*dot_product(bk(:,ikpt1)+g,tau(:,iisym)))
      iatom  = 0
      iiatom = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom  = iatom + 1
          cdum   =cexp*exp(-img*tpi_const*dot_product(g,taual(:,iatom)))

          rtaual = matmul(rot,taual(:,iatom)) + tau(:,iisym)
          iatom1 = 0
          DO ieq1 = 1,neq(itype)
            IF( all(abs(modulo(rtaual-taual(:,iiatom+ieq1)+1d-12,1d0))
     &                                               .lt. 1d-10) ) THEN ! The 1d-12 is a dirty fix.
              iatom1 = iiatom + ieq1
            END IF
          END DO
          IF( iatom1 .eq. 0 ) STOP 'ket_trafo: rotated atom not found'

          DO l = 0,lcutm(itype)
            nn = nindxm(l,itype)
            DO n = 1,nn

              i1 = pnt(n,l,iatom)
              i2 = i1 + nn * 2*l
              j1 = pnt(n,l,iatom1)
              j2 = j1 + nn * 2*l

              vecout1(i1:i2:nn) = cdum * matmul( vecin1(j1:j2:nn),
     +                                           dwgn(-l:l,-l:l,l) )

            END DO
          END DO
        END DO
        iiatom = iiatom + neq(itype)
      ENDDO

      ! PW
      DO igptm = 1,ngptm(ikpt1)
        igptp  = pgptm(igptm,ikpt1)
        g1     = matmul(invrrot,gptm(:,igptp)-g)
        igptm2 = 0
        DO i=1,ngptm(ikpt0)
          IF ( maxval( abs(g1-gptm(:,pgptm(i,ikpt0) )) ) <= 1E-06 ) THEN
            igptm2 = i
            EXIT
          END IF
        END DO 
        IF(igptm2.eq.0) 
     +               STOP 'ket_trafo: G point not found in G-point set.'

        cdum = exp(img * tpi_const * 
     +              dot_product(bk(:,ikpt1)+gptm(:,igptp),tau(:,iisym)))

        vecout1(nbasp+igptm) = cdum *         vecin1(nbasp+igptm2)

      END DO

      ! If inversion symmetry is applicable, define the phase of vecout and symmetrize to make the values real.
#ifdef CPP_INVERSION
       IF(lsymmetrize) CALL symmetrize(vecout1,1,nbasm(ikpt1),2,lreal,
     &                                 ntype,ntype,neq,lcutm,maxlcutm,
     &                                 nindxm,natd,invsat,invsatnr)
#endif
      vecout = vecout1

      END SUBROUTINE ket_trafo


      SUBROUTINE ket_trafo1(vecout,vecin,ikpt0,isym,lreal,lsymmetrize,
     >                     nwd,nw,ntype,nkpt,nkpt3,bk,nop,nsym,neq,mrot,
     >                     d_wgn,lmaxd,tau,invtab,invsat,invsatnr,map,
     >                     bmat,maxindxm,nindxm,nbasm,ngptm,ngptmall,
     >                     nbasp,pgptm,gptm,lcutm,maxlcutm,natd,taual)

      USE m_constants
      USE m_util      ,ONLY: modulo1
      USE m_dwigner

      IMPLICIT NONE

      ! -scalars -
      INTEGER , INTENT(IN)  ::  ikpt0,isym
      INTEGER , INTENT(IN)  ::  ntype,nkpt,nop,natd,nwd,nw,nbasp,nsym
      INTEGER , INTENT(IN)  ::  maxlcutm,maxindxm
      INTEGER , INTENT(IN)  ::  lmaxd
      LOGICAL , INTENT(IN)  ::  lreal,lsymmetrize

      ! - arrays -

      INTEGER , INTENT(IN)  ::  mrot(3,3,nop),invtab(nop)
      INTEGER , INTENT(IN)  ::  invsat(natd),invsatnr(natd)
      INTEGER , INTENT(IN)  ::  neq(ntype)
      INTEGER , INTENT(IN)  ::  nbasm(nkpt),ngptm(nkpt)
      INTEGER , INTENT(IN)  ::  nindxm(0:maxlcutm,ntype)
      INTEGER , INTENT(IN)  ::  pgptm(maxval(ngptm),nkpt)
      INTEGER , INTENT(IN)  ::  ngptmall
      INTEGER , INTENT(IN)  ::  gptm(3,ngptmall)
      INTEGER , INTENT(IN)  ::  lcutm(ntype)
      INTEGER , INTENT(IN)  ::  nkpt3(3)
      INTEGER , INTENT(IN)  ::  map(natd,nsym)

      REAL    , INTENT(IN)  ::  bk(3,nkpt),bmat(3,3)
      REAL    , INTENT(IN)  ::  taual(3,natd),tau(3,nop)

      COMPLEX , INTENT(IN)  ::  d_wgn(-lmaxd:lmaxd,-lmaxd:lmaxd,
     +                                     0:lmaxd,nsym)
      COMPLEX , INTENT(IN)  ::  vecin(nbasm(ikpt0))
      COMPLEX , INTENT(OUT) ::  vecout(nbasm(ikpt0))

      ! -local scalars -
      INTEGER               ::  iatom,iatom1,iiatom,itype,ieq,ieq1,ic,l,
     +                          n,i,nn,i1,i2,j1,j2
      INTEGER               ::  igptm,igptm1,igptm2,igptp,igptp1,ikpt1,
     +                          isymi,iisym
      INTEGER               ::  m1,m2

      REAL                  ::  time1,time2
      COMPLEX               ::  cexp,cdum 
      COMPLEX , PARAMETER   ::  img=(0d0,1d0)

      ! - local arrays -
      INTEGER               ::  pnt(maxindxm,0:maxlcutm,natd),g(3),g1(3)
      INTEGER               ::  rrot(3,3)

      REAL                  ::  rkpt(3),rkpthlp(3),rtaual(3)

      COMPLEX               ::  vecin1 (nbasm(ikpt0))
      COMPLEX               ::  vecout1(nbasm(ikpt0))
      COMPLEX               ::  dwgn(-maxlcutm:maxlcutm,
     +                               -maxlcutm:maxlcutm,
     +                                       0:maxlcutm)


      IF( maxlcutm .gt. lmaxd) STOP 'kettrafo1: maxlcutm > lmaxd'

  
c     Transform back to unsymmetrized product basis in case of inversion symmetry.
      vecin1 = vecin!(:nbasp)
# ifdef CPP_INVERSION
      IF(lsymmetrize) CALL desymmetrize(vecin1(:nbasp),1,nbasp,2,
     &                                  ntype,ntype,neq,lcutm,maxlcutm,
     &                                  nindxm,natd,invsat,invsatnr)
# endif

      IF( isym .le. nop ) THEN
        iisym   = isym
        rrot    = transpose( mrot(:,:,invtab(iisym)) )
      ELSE
        iisym   = isym - nop
        rrot    = -transpose( mrot(:,:,invtab(iisym)) )
      END IF

      rkpt    = matmul(rrot,bk(:,ikpt0))
      rkpthlp = modulo1(rkpt,nkpt3)
      g       = nint(rkpt - rkpthlp)

      DO l = 0,maxlcutm
        dwgn(-maxlcutm:maxlcutm,-maxlcutm:maxlcutm,l) =
     &    transpose(d_wgn(-maxlcutm:maxlcutm,-maxlcutm:maxlcutm,l,isym))
      END DO


      !
      ! determine number of rotated k-point bk(:,ikpt) -> ikpt1
      !
      DO i=1,nkpt
        IF ( maxval( abs(rkpthlp - bk(:,i)) ) .le. 1E-06  ) THEN
          ikpt1 = i
          EXIT
        END IF
      END DO

c     Define pointer to first mixed-basis functions (with m = -l)
      i  = 0
      ic = 0 
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          ic = ic + 1
          DO l = 0,lcutm(itype)
            DO n = 1,nindxm(l,itype)
              i           = i + 1
              pnt(n,l,ic) = i
            END DO
            i = i + nindxm(l,itype) * 2*l
          END DO
        END DO
      END DO

c     Multiplication
      ! MT
      cexp  =exp(-img*tpi_const*dot_product(bk(:,ikpt1)+g,tau(:,iisym)))
      iatom  = 0
      iiatom = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          iatom  = iatom + 1

          iatom1 = map(iatom,iisym)
          cdum   =cexp*exp(img*tpi_const*dot_product(g,taual(:,iatom1)))

          DO l = 0,lcutm(itype)
            nn = nindxm(l,itype)
            DO n = 1,nn

              i1 = pnt(n,l,iatom)
              i2 = i1 + nn * 2*l
              j1 = pnt(n,l,iatom1)
              j2 = j1 + nn * 2*l

              vecout1(i1:i2:nn) = cdum * matmul(vecin1(j1:j2:nn),
     +                                          dwgn(-l:l,-l:l,l))

            END DO
          END DO
        END DO
        iiatom = iiatom + neq(itype)
      ENDDO

      ! PW
      DO igptm = 1,ngptm(ikpt0)
        igptp  = pgptm(igptm,ikpt0)
        g1     = matmul(rrot,gptm(:,igptp)) +g
        igptm1 = 0
        DO i=1,ngptm(ikpt1)
          IF (maxval( abs(g1-gptm(:,pgptm(i,ikpt1) )) ) <= 1E-06 ) THEN
            igptm1 = i
            igptp1 = pgptm(i,ikpt1)
            EXIT
          END IF
        END DO
        IF(igptm1.eq.0)
     +               STOP 'ket_trafo: G point not found in G-point set.'

        cdum = exp(-img * tpi_const * 
     +             dot_product(bk(:,ikpt1)+gptm(:,igptp1),tau(:,iisym)))

        vecout1(nbasp+igptm) = cdum * vecin1(nbasp+igptm1)

      END DO 

      ! If inversion symmetry is applicable, define the phase of vecout and symmetrize to make the values real.
#ifdef CPP_INVERSION
       IF(lsymmetrize) CALL symmetrize(vecout1,1,nbasm(ikpt1),2,lreal,
     &                                 ntype,ntype,neq,lcutm,maxlcutm,
     &                                 nindxm,natd,invsat,invsatnr)
#endif
      IF( isym .le. nop ) THEN
        vecout = vecout1      
      ELSE
        vecout = conjg(vecout1)  
      END IF

      END SUBROUTINE ket_trafo1

      ! Determines common phase factor (with unit norm)
      SUBROUTINE commonphase(cfac,carr,n)

      IMPLICIT NONE
      INTEGER,INTENT(IN)      :: n
      COMPLEX,INTENT(IN)      :: carr(n)
      COMPLEX,INTENT(OUT)     :: cfac
      REAL                    :: rdum,rmax
      INTEGER                 :: i

!       IF( all( abs(carr) .lt. 1E-12 ) ) THEN
!         cfac = 1
!         RETURN
!       END IF

      cfac = 0
      rmax = 0
      DO i = 1,n
        rdum = abs(carr(i))
        IF     (rdum.gt.1d-6) THEN ; cfac = carr(i) / rdum ; EXIT
        ELSE IF(rdum.gt.rmax) THEN ; cfac = carr(i) / rdum ; rmax = rdum
        END IF
      END DO
      IF(cfac .eq. 0 .and. all(carr .ne. 0) ) THEN
        WRITE(999,*) carr
        STOP 'commonphase: Could not determine common phase factor. &
     &(Wrote carr to fort.999)'
      END IF
      END SUBROUTINE commonphase

      SUBROUTINE bramat_trafo(
     <      vecout,igptm_out,
     >      vecin,igptm_in,ikpt0,iop,writevec,pointer,invtab,
     >      nop,nsym,rrot,invrrot,map,nkptd,nkptf,nkpt3,bk,tau,maxlcutm,
     >      ntype,neq,lcutm,natd,taual,nindxm,maxindxm,nw,nwd,lmaxd,
     >      dwgn,ngptmall,gptm,ngptm,pgptm,nbasp,nbasm,invsat,invsatnr)

      USE m_constants     
      USE m_util

      IMPLICIT NONE

c     - scalars
      INTEGER,INTENT(IN)      ::  ikpt0,igptm_in,iop,nop,maxindxm
      INTEGER,INTENT(IN)      ::  nkptd,nkptf,maxlcutm,ntype,natd
      INTEGER,INTENT(IN)      ::  nw,nwd,nbasp,nsym,lmaxd
      INTEGER,INTENT(IN)      ::  ngptmall
      LOGICAL,INTENT(IN)      ::  writevec
      INTEGER,INTENT(OUT)     ::  igptm_out
c     - arrays -
      INTEGER,INTENT(IN)      ::  rrot(3,3),invrrot(3,3),invtab(nop)
      INTEGER,INTENT(IN)      ::  neq(ntype),lcutm(ntype),
     +                            nindxm(0:maxlcutm,ntype)
      INTEGER,INTENT(IN)      ::  gptm(3,ngptmall),
     &                            ngptm(nkptf),
     &                            pgptm(maxval(ngptm),nkptf)
      INTEGER,INTENT(IN)      ::  invsat(natd),invsatnr(natd)
      INTEGER,INTENT(IN)      ::  nkpt3(3),nbasm(nkptf)
      INTEGER,INTENT(IN)      ::  map(natd,nsym)
      INTEGER,INTENT(IN)      ::  pointer(
     &                          minval(gptm(1,:))-1:maxval(gptm(1,:))+1,
     &                          minval(gptm(2,:))-1:maxval(gptm(2,:))+1,
     &                          minval(gptm(3,:))-1:maxval(gptm(3,:))+1)

      REAL,   INTENT(IN)      ::  bk(3,nkptd),tau(3,nop),taual(3,natd)
      COMPLEX,INTENT(IN)      ::  vecin(nbasm(ikpt0))
      COMPLEX,INTENT(IN)      ::  dwgn(-maxlcutm:maxlcutm,
     +                                 -maxlcutm:maxlcutm,
     +                                         0:maxlcutm)
      COMPLEX,INTENT(OUT)     ::  vecout(nbasm(ikpt0))

c     - private scalars -
      INTEGER                 ::  itype,ieq,ic,l,n,i,nn,i1,i2,j1,j2,m1
      INTEGER                 ::  m2,igptm,igptm2,igptp,iiop,isym
      INTEGER                 ::  ikpt1,isymi,rcent
      LOGICAL                 ::  trs
      COMPLEX,PARAMETER       ::  img=(0d0,1d0)
      COMPLEX                 ::  cexp,cdum
c     - private arrays -
      INTEGER                 ::  pnt(maxindxm,0:maxlcutm,natd),g(3),
     +                            g1(3),iarr(ngptm(ikpt0))
      REAL                    ::  rkpt(3),rkpthlp(3),trans(3)
      COMPLEX                 ::  vecin1(nbasm(ikpt0))
      COMPLEX                 ::  carr(ngptm(ikpt0))

   
      IF( iop .le. nop ) THEN

        isym   = iop
        trs    = .false.
        trans  = tau(:,isym)

      ELSE

        isym   = iop - nop
        trs    = .true.
        trans  = tau(:,isym)

      END IF

      rkpthlp = matmul(rrot,bk(:,ikpt0))
      rkpt    = modulo1(rkpthlp,nkpt3)
      g       = nint(rkpthlp-rkpt)
      !
      ! determine number of rotated k-point bk(:,ikpt) -> ikpt1
      ! 
      DO i=1,nkptd
        IF ( maxval( abs(rkpt - bk(:,i)) ) .le. 1E-06  ) THEN
          ikpt1 = i
          EXIT
        END IF
      END DO

      DO igptm = 1,ngptm(ikpt1)
        igptp  = pgptm(igptm,ikpt1)
        g1     = matmul(invrrot,gptm(:,igptp)-g)
        igptm2 = pointer(g1(1),g1(2),g1(3))
        IF(igptm2.eq.igptm_in) THEN
          igptm_out = igptm
          IF(writevec) THEN
            cdum = exp(img * tpi_const * 
     +                     dot_product(bk(:,ikpt1)+gptm(:,igptp),trans))
            EXIT
          ELSE
            RETURN
          END IF
        END IF
      END DO

c     Transform back to unsymmetrized product basis in case of inversion symmetry.
      vecout = vecin
# ifdef CPP_INVERSION
      CALL desymmetrize(vecout,nbasp,1,1,
     &                  ntype,ntype,neq,lcutm,maxlcutm,
     &                  nindxm,natd,invsat,invsatnr)
# endif

c     Right-multiplication
      ! PW
      IF(trs) THEN ; vecin1 = cdum * conjg ( vecout )
      ELSE         ; vecin1 = cdum *         vecout 
      END IF

c     Define pointer to first mixed-basis functions (with m = -l)
      i  = 0
      ic = 0 
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          ic = ic + 1
          DO l = 0,lcutm(itype)
            DO n = 1,nindxm(l,itype)
              i           = i + 1
              pnt(n,l,ic) = i
            END DO
            i = i + nindxm(l,itype) * 2*l
          END DO
        END DO
      END DO

c     Left-multiplication
      ! MT
      cexp = exp( -img *tpi_const * dot_product(bk(:,ikpt1)+g,trans))
      ic   = 0
      DO itype = 1,ntype
        DO ieq = 1,neq(itype)
          ic    = ic + 1
          rcent = map(ic,invtab(isym))
          cdum  = cexp * exp(img *tpi_const*dot_product(g,taual(:,ic))) !rcent)))
          cdum = conjg(cdum)
          DO l = 0,lcutm(itype)
            nn = nindxm(l,itype)
            DO n = 1,nn

              i1 = pnt(n,l,ic)
              i2 = i1 + nn * 2*l
              j1 = pnt(n,l,rcent)
              j2 = j1 + nn * 2*l

              vecout(i1:i2:nn) = cdum * matmul(dwgn(-l:l,-l:l,l),
     +                                         vecin1(j1:j2:nn) )

            END DO
          END DO
        END DO
      END DO

      ! PW
      DO igptm = 1,ngptm(ikpt1)
        igptp       = pgptm(igptm,ikpt1)
        g1          = matmul(invrrot,gptm(:,igptp)-g)
        iarr(igptm) = pointer(g1(1),g1(2),g1(3))
        carr(igptm) = exp(-img * tpi_const * 
     +                     dot_product(bk(:,ikpt1)+gptm(:,igptp),trans))
      END DO
      DO i1 = 1,ngptm(ikpt1)
        vecout(nbasp+i1) = carr(i1) * vecin1(nbasp+iarr(i1))
      END DO

      ! If inversion symmetry is applicable, symmetrize to make the values real.
# ifdef CPP_INVERSION
      CALL symmetrize(vecout,nbasp,1,1,.false.,
     &                ntype,ntype,neq,lcutm,maxlcutm,
     &                nindxm,natd,invsat,invsatnr)
# endif

      END SUBROUTINE bramat_trafo

      END MODULE m_trafo
