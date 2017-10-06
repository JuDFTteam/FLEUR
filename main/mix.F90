!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mix

  !*************************************************************************
  !  mixing of charge densities or potentials:
  !    IMIX= 0 : linear mixing                                     
  !    IMIX = 3 : BROYDEN'S FIRST METHOD                            
  !    IMIX = 5 : BROYDEN'S SECOND METHOD                           
  !    IMIX = 7 : GENERALIZED ANDERSEN METHOD                       
  !************************************************************************

CONTAINS

SUBROUTINE mix(stars,atoms,sphhar,vacuum,input,sym,cell,noco,oneD,&
               hybrid,inDen,outDen,results)

#include"cpp_double.h"

   USE m_juDFT
   USE m_constants
   USE m_cdn_io
   USE m_brysh1
   USE m_stmix
   USE m_broyden
   USE m_brysh2
   USE m_metric
   USE m_qfix
   USE m_types
   USE m_xmlOutput

   IMPLICIT NONE

   TYPE(t_oneD),INTENT(IN)       :: oneD
   TYPE(t_hybrid),INTENT(IN)     :: hybrid 
   TYPE(t_input),INTENT(IN)      :: input
   TYPE(t_vacuum),INTENT(IN)     :: vacuum
   TYPE(t_noco),INTENT(IN)       :: noco
   TYPE(t_sym),INTENT(IN)        :: sym
   TYPE(t_stars),INTENT(IN)      :: stars
   TYPE(t_cell),INTENT(IN)       :: cell
   TYPE(t_sphhar),INTENT(IN)     :: sphhar
   TYPE(t_atoms),INTENT(INOUT)   :: atoms !n_u is modified temporarily
   TYPE(t_potden),INTENT(IN)     :: inDen, outDen
   TYPE(t_results),INTENT(INOUT) :: results

   !Local type instances
   TYPE(t_potden) :: mixDen

   !Local Scalars
   REAL fix,intfac,vacfac
   INTEGER i,iter,imap,js,mit,irecl
   INTEGER mmap,mmaph,nmaph,nmap,mapmt,mapvac,mapvac2
   INTEGER iofl,archiveType,n_u_keep
   LOGICAL lexist,l_ldaU

   !Local Arrays
   REAL dist(6)
   REAL, ALLOCATABLE :: sm(:), fsm(:)
   CHARACTER(LEN=20) :: attributes(2)

   !External functions
   REAL CPP_BLAS_sdot
   EXTERNAL CPP_BLAS_sdot

   ! YM: I have exported 'vol' from outside, be aware
   !     IF (film) THEN
   !        vol = 2.0 * z1 * area
   !     ELSE
   !        vol = omtil
   !     ENDIF

   !In systems without inversions symmetry the interstitial star-
   !coefficients are complex. Thus twice as many numbers have to be
   !stored.
   IF (sym%invs) THEN
      intfac = 1.0
   ELSE
      intfac = 2.0
   END IF

   !The corresponding is true for the coeff. of the warping vacuum
   !density depending on the two dimensional inversion.
   IF (sym%invs2) THEN
      vacfac = 1.0
   ELSE
      vacfac = 2.0
   END IF

   mmaph = intfac*stars%ng3 + atoms%ntype*(sphhar%nlhd+1)*atoms%jmtd +&
           vacfac*vacuum%nmzxyd*(oneD%odi%n2d-1)*vacuum%nvac + vacuum%nmzd*vacuum%nvac
   mmap  = mmaph*input%jspins
   !in a non-collinear calculations extra space is needed for the
   !off-diag. part of the density matrix. these coeff. are generally
   !complex independ of invs and invs2.
   IF (noco%l_noco) THEN
      mmap = mmap + 2*stars%ng3 + 2*vacuum%nmzxyd*(oneD%odi%n2d-1)*vacuum%nvac + &
             2*vacuum%nmzd*vacuum%nvac
   END IF

   n_u_keep=atoms%n_u
   INQUIRE (file='n_mmp_mat',exist=l_ldaU) 
   IF (l_ldaU) THEN
      !In an LDA+U caclulation, also the density matrix is included in the
      !supervectors (sm,fsm) if no mixing factors are in the n_mmp_mat-file

      OPEN (69,file='n_mmp_mat',status='old',form='formatted')
      i = 0 
      DO
         READ (69,*,iostat=iofl)
         IF (iofl < 0) EXIT
         i = i + 1
      END DO
      CLOSE (69)

      IF ( MOD(i,14*input%jspins) == 1 ) THEN     ! was already mixed in u_mix
         atoms%n_u = 0
      ELSE IF ( MOD(i,28*input%jspins)== 0 ) THEN ! mix here 
         atoms%n_u = i / (28 * input%jspins ) ! calculate number of U parameters
         mmap = mmap + 7 * i / 2   ! add 7*7 complex numbers per atoms%n_u
      ELSE
         CALL juDFT_error("strange n_mmp_mat-file...",calledby ="mix")
      END IF
   ELSE
      atoms%n_u=0
   END IF ! l_ldaU

   ALLOCATE (sm(mmap),fsm(mmap))

   INQUIRE (file='broyd.'//CHAR(input%imix+48),exist=lexist)
   DO i = 1,6
      dist(i) = 0.0
   END DO

   !determine type of mixing:
   !imix=0:straight, imix=o broyden first, imix=5:broyden second
   !imix=:generalozed anderson mixing
   mit = 0
   IF (input%imix.EQ.0) THEN
      WRITE (16,FMT='(a,2f10.5)') 'STRAIGHT MIXING',input%alpha
   ELSE IF (input%imix.EQ.3) THEN
      IF ( .NOT.lexist) mit = 1
      WRITE (16,FMT='(a,f10.5)') 'BROYDEN FIRST MIXING',input%alpha
   ELSE IF (input%imix.EQ.5) THEN
      IF (.NOT.lexist) mit = 1
      WRITE (16,FMT='(a,f10.5)') 'BROYDEN SECOND MIXING',input%alpha
   ELSE IF (input%imix.EQ.7) THEN
      IF (.NOT.lexist) mit = 1
      WRITE (16,FMT='(a,f10.5)') 'ANDERSON GENERALIZED',input%alpha
   ELSE
      CALL juDFT_error("mix: input%imix =/= 0,3,5,7 ",calledby ="mix")
   END IF

   IF (input%jspins.EQ.2.AND.input%imix.NE.0) THEN
      WRITE(6,'(''WARNING : for QUASI-NEWTON METHODS SPINF=1'')')
   END IF

   !put input charge density into arrays sm 
   !in the spin polarized case the arrays consist of 
   !spin up and spin down densities
   CALL brysh1(input,stars,atoms,sphhar,noco,vacuum,sym,oneD,&
               intfac,vacfac,inDen%pw,inDen%mt,inDen%vacz,inDen%vacxy,inDen%cdom,&
               inDen%cdomvz,inDen%cdomvxy,inDen%mmpMat(-lmaxU_const,-lmaxU_const,1,1),&
               nmap,nmaph,mapmt,mapvac,mapvac2,sm) 

   !put output charge density into arrays fsm 
   CALL brysh1(input,stars,atoms,sphhar,noco,vacuum,sym,oneD,&
               intfac,vacfac,outDen%pw,outDen%mt,outDen%vacz,outDen%vacxy,outDen%cdom,&
               outDen%cdomvz,outDen%cdomvxy,outDen%mmpMat(-lmaxU_const,-lmaxU_const,1,1),&
               nmap,nmaph,mapmt,mapvac,mapvac2,fsm)

   !store fsm - sm the difference on fsm
   DO imap = 1,nmap
      fsm(imap) = fsm(imap) - sm(imap)
   END DO

   ! open files for broyden
   irecl=(nmap+1)*8
   IF (input%imix.GE.3) THEN
      IF (hybrid%l_calhf) THEN
         OPEN (57,file='hf_broyd',form='unformatted',status='unknown')
         OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
               recl=irecl,form='unformatted',status='unknown')
      ELSE
         OPEN (57,file='broyd',form='unformatted',status='unknown')
         OPEN (59,file='broyd.'//CHAR(input%imix+48),access='direct',&
               recl=irecl,form='unformatted',status='unknown')
      ENDIF
   END IF

   !mixing of the densities
   IF (input%imix.EQ.0) THEN
      CALL stmix(atoms,input,noco, nmap,nmaph,fsm, sm)
   ELSE
      CALL broyden(cell,stars,atoms,vacuum,sphhar,input,noco,oneD,sym,&
                   mmap,nmaph,mapmt,mapvac2,nmap,fsm,mit,sm)
   END IF

   !initiatlize mixed density and extract it with brysh2 call
   CALL mixDen%init(stars,atoms,sphhar,vacuum,oneD,input%jspins,.FALSE.)
   IF (noco%l_noco) THEN
      ALLOCATE (mixDen%cdom(stars%ng3),mixDen%cdomvz(vacuum%nmzd,2))
      ALLOCATE (mixDen%cdomvxy(vacuum%nmzxyd,oneD%odi%n2d-1,2))
      archiveType = CDN_ARCHIVE_TYPE_NOCO_const
   ELSE
      ALLOCATE (mixDen%cdom(1),mixDen%cdomvz(1,1),mixDen%cdomvxy(1,1,1))
      archiveType = CDN_ARCHIVE_TYPE_CDN1_const
   ENDIF
   IF (atoms%n_u.GT.0) THEN
      ALLOCATE (mixDen%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_u,input%jspins))
   ELSE
      ALLOCATE (mixDen%mmpMat(-lmaxU_const:-lmaxU_const,-lmaxU_const:-lmaxU_const,1,input%jspins))
   END IF
   mixDen%cdom = CMPLX(0.0,0.0)
   mixDen%cdomvz = CMPLX(0.0,0.0)
   mixDen%cdomvxy = CMPLX(0.0,0.0)
   mixDen%mmpMat = CMPLX(0.0,0.0)

   CALL brysh2(input,stars,atoms,sphhar,noco,vacuum,sym,sm,mixDen%mmpMat,oneD,&
               mixDen%pw,mixDen%mt,mixDen%vacz,mixDen%vacxy,mixDen%cdom,&
               mixDen%cdomvz,mixDen%cdomvxy) 

   !calculate the distance of charge densities...

   !induce metric in fsm use sm as an output array: |sm> = w |fsm>
   CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
               mmap,nmaph,mapmt,mapvac2,fsm, sm)

   !calculate the charge density distance for each spin
   IF(hybrid%l_calhf) THEN
      CALL openXMLElement('densityConvergence',(/'units  ','comment'/),(/'me/bohr^3','HF       '/))
   ELSE
      CALL openXMLElement('densityConvergence',(/'units'/),(/'me/bohr^3'/))
   END IF
   iter = inDen%iter
   DO js = 1,input%jspins
      dist(js) = CPP_BLAS_sdot(nmaph,fsm(nmaph*(js-1)+1),1, sm(nmaph*(js-1)+1),1)

      attributes = ''
      WRITE(attributes(1),'(i0)') js
      WRITE(attributes(2),'(f20.10)') 1000*SQRT(ABS(dist(js)/cell%vol))
      CALL writeXMLElementForm('chargeDensity',(/'spin    ','distance'/),attributes,reshape((/4,8,1,20/),(/2,2/)))
      IF( hybrid%l_calhf ) THEN
         WRITE (16,FMT=7901) js,iter,1000*SQRT(ABS(dist(js)/cell%vol))
         WRITE ( 6,FMT=7901) js,iter,1000*SQRT(ABS(dist(js)/cell%vol))
      ELSE
         WRITE (16,FMT=7900) js,iter,1000*SQRT(ABS(dist(js)/cell%vol))
         WRITE ( 6,FMT=7900) js,iter,1000*SQRT(ABS(dist(js)/cell%vol))
      END IF
   END DO
   IF (noco%l_noco) dist(6) = CPP_BLAS_sdot((nmap-2*nmaph), fsm(nmaph*2+1),1,sm(nmaph*2+1),1)
   IF (noco%l_noco) WRITE (6,FMT=7900) 3,iter,1000*SQRT(ABS(dist(6)/cell%vol))

   !calculate the distance of total charge and spin density
   !|rho/m(o) - rho/m(i)| = |rh1(o) -rh1(i)|+ |rh2(o) -rh2(i)| +/_
   !                        +/_2<rh2(o) -rh2(i)|rh1(o) -rh1(i)>
   IF (input%jspins.EQ.2) THEN
      dist(3) = CPP_BLAS_sdot(nmaph,fsm,1,sm(nmaph+1),1)
      dist(4) = dist(1) + dist(2) + 2.0e0*dist(3)
      dist(5) = dist(1) + dist(2) - 2.0e0*dist(3)
      CALL writeXMLElementFormPoly('overallChargeDensity',(/'distance'/),&
                                   (/1000*SQRT(ABS(dist(4)/cell%vol))/),reshape((/10,20/),(/1,2/)))
      CALL writeXMLElementFormPoly('spinDensity',(/'distance'/),&
                                   (/1000*SQRT(ABS(dist(5)/cell%vol))/),reshape((/19,20/),(/1,2/)))
      IF( hybrid%l_calhf ) THEN
         WRITE (16,FMT=8001) iter,1000*SQRT(ABS(dist(4)/cell%vol))
         WRITE (16,FMT=8011) iter,1000*SQRT(ABS(dist(5)/cell%vol))
         WRITE ( 6,FMT=8001) iter,1000*SQRT(ABS(dist(4)/cell%vol))
         WRITE ( 6,FMT=8011) iter,1000*SQRT(ABS(dist(5)/cell%vol))
      ELSE
         WRITE (16,FMT=8000) iter,1000*SQRT(ABS(dist(4)/cell%vol))
         WRITE (16,FMT=8010) iter,1000*SQRT(ABS(dist(5)/cell%vol))
         WRITE ( 6,FMT=8000) iter,1000*SQRT(ABS(dist(4)/cell%vol))
         WRITE ( 6,FMT=8010) iter,1000*SQRT(ABS(dist(5)/cell%vol))
      END IF

      !dist/vol should always be >= 0 ,
      !but for dist=0 numerically you might obtain dist/vol < 0
      !(e.g. when calculating non-magnetic systems with jspins=2).
   END IF
   results%last_distance=maxval(1000*SQRT(ABS(dist/cell%vol)))
   DEALLOCATE (sm,fsm)
   CALL closeXMLElement('densityConvergence')

   !fix charge of the new density
   CALL qfix(stars,atoms,sym,vacuum, sphhar,input,cell,oneD,&
             mixDen%pw,mixDen%vacxy,mixDen%mt,mixDen%vacz,.FALSE.,.false., fix)

   !write out mixed density
   CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                     1,results%last_distance,results%ef,.TRUE.,iter,mixDen%mt,mixDen%pw,mixDen%vacz,&
                     mixDen%vacxy,mixDen%cdom,mixDen%cdomvz,mixDen%cdomvxy)

   IF (atoms%n_u > 0) THEN
      OPEN (69,file='n_mmp_mat',status='replace',form='formatted')
      WRITE (69,'(7f20.13)') mixDen%mmpMat(:,:,:,:)
      CLOSE (69)
   ENDIF

   IF (input%imix.GT.0) THEN
      CLOSE (57)
      CLOSE (59)
   END IF

   atoms%n_u=n_u_keep

   7900 FORMAT (/,'---->    distance of charge densities for spin ',i2,'                 it=',i5,':',f13.6,' me/bohr**3')
   7901 FORMAT (/,'----> HF distance of charge densities for spin ',i2,'                 it=',i5,':',f13.6,' me/bohr**3')
   8000 FORMAT (/,'---->    distance of charge densities for it=',i5,':', f13.6,' me/bohr**3')
   8001 FORMAT (/,'----> HF distance of charge densities for it=',i5,':', f13.6,' me/bohr**3')
   8010 FORMAT (/,'---->    distance of spin densities for it=',i5,':', f13.6,' me/bohr**3')
   8011 FORMAT (/,'----> HF distance of spin densities for it=',i5,':', f13.6,' me/bohr**3')
   8020 FORMAT (4d25.14)
   8030 FORMAT (10i10)

END SUBROUTINE mix

END MODULE m_mix
