      MODULE m_bs_comfort
      USE m_juDFT

      CONTAINS
      SUBROUTINE bs_comfort(
     > neigd,nkptd,jspd,nwdd,ntypd,lmaxd,nlod,
     > jspins,nwd,ntype,nkpt,ntapwf,irecl,
     > delgau,film,l_noco,l_ss,
     X l_disp,nkpt_l,eig_l,irank,isize,param,socfile)


#ifdef CPP_HDF
      USE m_eig66_hdf, ONLY : read_eig, read_kptenpar
#endif

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,jspd,nwdd,ntypd,lmaxd,nlod
      INTEGER, INTENT (IN) :: jspins,nwd,ntype,irecl,ntapwf,nkpt_l
      LOGICAL, INTENT (IN) :: film,socfile
      LOGICAL, INTENT (IN) :: l_noco,l_ss,l_disp
      REAL,    INTENT (IN) :: delgau
      INTEGER, INTENT (IN) :: irank,isize,param
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nkpt(nwdd)
      REAL,    INTENT (IN) :: eig_l(neigd+5,nkpt_l)

C     ..
C     .. Local Scalars ..
      REAL del,seigsc,seigv,spindg,ssc,weight,ws,zc,tkb_1
      INTEGER i,idummy,j,jsp,k,l,n,nbands,nstef,nv,nw,nrec,nmat,nspins
      INTEGER n_help
C     ..
C     .. Local Arrays ..
C
      REAL,    ALLOCATABLE :: eig(:,:,:)
      INTEGER ne(nkptd,jspd)
      REAL bk(3,nkptd),el(0:lmaxd,ntypd,jspd),evac(2,jspd)
      REAL wtkpt(nkptd),ello(nlod,ntypd,jspd)

      LOGICAL :: l_etotskip = .false.
      REAL    :: etotskip_val = 0.0

C     ..
c***********************************************************************
c                  ABBREVIATIONS
c
c     eig        : array of eigenvalues within all energy-windows
c     wtkpt      : list of the weights of each k-point (from inp-file)
c     e          : linear list of the eigenvalues within the highest
c                  energy-window
c     we         : list of weights of the eigenvalues in e
c     zelec      : number of electrons in a window
c     spindg     : spindegeneracy (2 in nonmagnetic calculations)
c     seigv      : weighted sum of the occupied valence eigenvalues
c     seigsc     : weighted sum of the semi-core eigenvalues
c     seigscv    : sum of seigv and seigsc
C     ts         : entropy contribution to the free energy
c     
c***********************************************************************
C     .. Data statements ..
      DATA del/1.0e-6/
C     ..
      ALLOCATE (eig(neigd,nkptd,jspd))
c
      IF (socfile) THEN
        IF (nwd/=1) THEN
           CALL juDFT_error("eig.soc and multiple windows",calledby
     +          ="bs_comfort")
        ENDIF
        OPEN(67,file='eig.soc',form='unformatted',action='read')
      ENDIF
c
c---> READ IN EIGENVALUES
c
      spindg = 2.0/real(jspins)
      n = 0
      nrec = 0
      seigsc = 0.0
      ssc = 0.0
      n_help = 0
c
c---> pk non-collinear
      IF (l_noco) THEN
         nspins = 1
      ELSE
         nspins = jspins
      ENDIF
c---> pk non-collinear
c
      DO 50 nw = 1,nwd
        DO 40 jsp = 1,nspins
          DO 30 k = 1,nkpt(nw)

 100         CONTINUE

             IF (socfile) THEN

               IF (jsp==1) THEN
                 READ(67) i
                 IF (i/=k) THEN
                    CALL juDFT_error("error in eig.soc",
     +                               calledby ="bs_comfort")
                 ENDIF
                 READ(67) bk(1,k), bk(2,k), bk(3,k)
                 READ(67) wtkpt(k)
                 READ(67) ne(k,jsp)
                 DO i= 1,ne(k,jsp)
                   READ(67) eig(i,k,jsp)
                 ENDDO
               ELSE
                 ne(k,jsp)= ne(k,1)
                 DO i= 1,ne(k,1)
                   eig(i,k,jsp)= eig(i,k,1)
                 ENDDO
               ENDIF

               nv= -1

             ELSE
               nrec = nrec + 1
#ifdef CPP_HDF
               IF (nw>1)  CALL juDFT_error("HDF and multiple windows.."
     +              ,calledby ="bs_comfort")
               CALL read_eig(k,jsp,ne(k,jsp),eig(:,k,jsp))
               CALL read_kptenpar(k,jsp,bk(1,k),wtkpt(k),el(:,:,jsp),
     +                            ello(:,:,jsp),evac(1,jsp))
               WRITE (6,'(a2,3f10.5,f12.6)') 'at',bk(:,k),wtkpt(k)
               WRITE (6,'(i5,a14)') ne(k,jsp),' eigenvalues :' 
               WRITE (6,'(8f12.6)') (eig(i,k,jsp),i=1,ne(k,jsp))
               nv= -1 
#else
                 IF (l_ss) THEN
                  READ (ntapwf,rec=nrec) el,evac,ello,
     +                               (bk(i,k),i=1,3),wtkpt(k),
     +                                ne(k,jsp),nv,idummy,nmat,
     +                               (eig(i,k,jsp),i=1,neigd)
                 ELSEIF (l_noco) THEN
                  READ (ntapwf,rec=nrec) el,evac,ello,
     +                               (bk(i,k),i=1,3),wtkpt(k),
     +                                ne(k,jsp),nv,nmat,
     +                               (eig(i,k,jsp),i=1,neigd)
                 ELSE
                  READ (ntapwf,rec=nrec) el(:,:,jsp),evac(:,jsp),
     +                               ello(:,:,jsp),
     +                               (bk(i,k),i=1,3),wtkpt(k),
     +                                ne(k,jsp),nv,nmat,
     +                               (eig(i,k,jsp),i=1,neigd)
                 ENDIF
#endif
              ENDIF!(socfile)

#ifdef CPP_MPI
               IF (.NOT.socfile) THEN 
                 n_help = n_help + ne(k,jsp)
                 IF (abs(evac(1,jsp)-999.9).LT.1.e-9) THEN
c
c                  obviously this record was only one part of a k-point
c
                   GOTO 100
                 ENDIF
                 ne(k,jsp) = n_help
                 n_help = 0
c                write(*,'(8f10.5)') (e(n-ne(k,jsp)+i),i=1,ne(k,jsp))
               ENDIF
#endif


   30       CONTINUE

            IF(nw.EQ.nwd) THEN
            DO 20 i=1,neigd
               DO 10 k=1,nkpt(nw)
                  write(776+jsp,*)param,k,eig(i,k,jsp)
 10            CONTINUE
               write(776+jsp,*)
 20         CONTINUE

            DO 21 i=1,neigd
               DO 11 k=1,nkpt(nw)
                  if(k.eq.param)write(778+jsp,*)param,k,eig(i,k,jsp)
 11            CONTINUE
               write(778+jsp,*)
 21         CONTINUE
            ENDIF

   40    CONTINUE
   50 CONTINUE

      DEALLOCATE ( eig )

      RETURN
      END SUBROUTINE bs_comfort
      END MODULE m_bs_comfort
