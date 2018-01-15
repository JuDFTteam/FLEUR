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
               hybrid,archiveType,inDen,outDen,results)

#include"cpp_double.h"

   USE m_juDFT
   USE m_constants
   USE m_cdn_io
   USE m_broyd_io
   USE m_brysh1
   USE m_stmix
   USE m_broyden
   USE m_broyden2
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
   TYPE(t_potden),INTENT(IN)     :: outDen
   TYPE(t_results),INTENT(INOUT) :: results
   TYPE(t_potden),INTENT(INOUT)  :: inDen
   INTEGER, INTENT(IN)           :: archiveType

   !Local Scalars
   REAL fix,intfac,vacfac
   INTEGER i,imap,js
   INTEGER mmap,mmaph,nmaph,nmap,mapmt,mapvac,mapvac2
   INTEGER iofl,n_u_keep
   LOGICAL l_exist,l_ldaU, l_densityMatrixPresent, l_pot

   !Local Arrays
   REAL dist(6)
   REAL, ALLOCATABLE :: sm(:), fsm(:), fmMet(:), smMet(:)
   CHARACTER(LEN=20) :: attributes(2)
   COMPLEX           :: n_mmpTemp(-3:3,-3:3,MAX(1,atoms%n_u),input%jspins)

   !External functions
   REAL CPP_BLAS_sdot
   EXTERNAL CPP_BLAS_sdot

   ! YM: I have exported 'vol' from outside, be aware
   !     IF (film) THEN
   !        vol = 2.0 * z1 * area
   !     ELSE
   !        vol = omtil
   !     ENDIF

   l_densityMatrixPresent = ANY(inDen%mmpMat(:,:,:,:).NE.0.0)

   !In systems without inversions symmetry the interstitial star-
   !coefficients are complex. Thus twice as many numbers have to be
   !stored.
   intfac = 2.0
   IF (sym%invs) intfac = 1.0

   !The corresponding is true for the coeff. of the warping vacuum
   !density depending on the two dimensional inversion.
   vacfac = 2.0
   IF (sym%invs2) vacfac = 1.0

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

   ! LDA+U (start)
   n_mmpTemp = inDen%mmpMat
   n_u_keep=atoms%n_u
   IF (l_densityMatrixPresent) THEN
      !In an LDA+U caclulation, also the density matrix is included in the
      !supervectors (sm,fsm) if no linear mixing is performed on it.
      IF(input%ldauLinMix) THEN
         atoms%n_u = 0
      ELSE
         mmap = mmap + 7 * 7 * 2 * atoms%n_u * input%jspins ! add 7*7 complex numbers per atoms%n_u and spin
      END IF
   ELSE
      atoms%n_u = 0
   END IF
   ! LDA+U (end)

   ALLOCATE (sm(mmap),fsm(mmap))
   ALLOCATE (smMet(mmap),fmMet(mmap))
   dist(:) = 0.0

   !determine type of mixing:
   !imix=0:straight, imix=o broyden first, imix=5:broyden second
   !imix=:generalozed anderson mixing
   IF (input%imix.EQ.0) THEN
      WRITE (16,FMT='(a,2f10.5)') 'STRAIGHT MIXING',input%alpha
   ELSE IF (input%imix.EQ.3) THEN
      WRITE (16,FMT='(a,f10.5)') 'BROYDEN FIRST MIXING',input%alpha
   ELSE IF (input%imix.EQ.5) THEN
      WRITE (16,FMT='(a,f10.5)') 'BROYDEN SECOND MIXING',input%alpha
   ELSE IF (input%imix.EQ.7) THEN
      WRITE (16,FMT='(a,f10.5)') 'ANDERSON GENERALIZED',input%alpha
   ELSE
      CALL juDFT_error("mix: input%imix =/= 0,3,5,7 ",calledby ="mix")
   END IF

   IF (input%jspins.EQ.2.AND.input%imix.NE.0) THEN
      WRITE(6,'(''WARNING : for QUASI-NEWTON METHODS SPINF=1'')')
   END IF

   !put input charge density into array sm

   !(in the spin polarized case the arrays sm and fsm consist of spin up and spin down densities)
   CALL brysh1(input,stars,atoms,sphhar,noco,vacuum,sym,oneD,&
               intfac,vacfac,inDen,nmap,nmaph,mapmt,mapvac,mapvac2,sm) 

   !put output charge density into array fsm
   CALL brysh1(input,stars,atoms,sphhar,noco,vacuum,sym,oneD,&
               intfac,vacfac,outDen,nmap,nmaph,mapmt,mapvac,mapvac2,fsm)

   !store the difference fsm - sm in fsm
   fsm(:nmap) = fsm(:nmap) - sm(:nmap)

   l_pot = .FALSE.
   ! Apply metric w to fsm and store in fmMet:  w |fsm>
   CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
               mmap,nmaph,mapmt,mapvac2,fsm,fmMet,l_pot)

   !mixing of the densities
   IF (input%imix.EQ.0) THEN
      CALL stmix(atoms,input,noco, nmap,nmaph,fsm, sm)
   ELSE
      CALL broyden(cell,stars,atoms,vacuum,sphhar,input,noco,oneD,sym,&
                   hybrid,mmap,nmaph,mapmt,mapvac2,nmap,fsm,sm)

!     Replace the broyden call above by the commented metric and broyden2 calls
!     below to switch on the continuous restart of the Broyden method.
!      ! Apply metric w to sm and store in smMet:  w |sm>
!      CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
!                  mmap,nmaph,mapmt,mapvac2,sm,smMet,l_pot)
!
!      CALL broyden2(cell,stars,atoms,vacuum,sphhar,input,noco,oneD,sym,&
!                    hybrid,mmap,nmaph,mapmt,mapvac2,nmap,fsm,sm,fmMet,smMet)
   END IF

   !initiatlize mixed density and extract it with brysh2 call
   inDen%cdom = CMPLX(0.0,0.0)
   inDen%cdomvz = CMPLX(0.0,0.0)
   inDen%cdomvxy = CMPLX(0.0,0.0)
   inDen%mmpMat = CMPLX(0.0,0.0)

   CALL brysh2(input,stars,atoms,sphhar,noco,vacuum,sym,sm,oneD,inDen) 

   !calculate the distance of charge densities...

   !induce metric in fsm use sm as an output array: |sm> = w |fsm>
!   CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
!               mmap,nmaph,mapmt,mapvac2,fsm, sm)

   !calculate the charge density distance for each spin
   IF(hybrid%l_calhf) THEN
      CALL openXMLElement('densityConvergence',(/'units  ','comment'/),(/'me/bohr^3','HF       '/))
   ELSE
      CALL openXMLElement('densityConvergence',(/'units'/),(/'me/bohr^3'/))
   END IF

   DO js = 1,input%jspins
      dist(js) = CPP_BLAS_sdot(nmaph,fsm(nmaph*(js-1)+1),1, fmMet(nmaph*(js-1)+1),1)

      attributes = ''
      WRITE(attributes(1),'(i0)') js
      WRITE(attributes(2),'(f20.10)') 1000*SQRT(ABS(dist(js)/cell%vol))
      CALL writeXMLElementForm('chargeDensity',(/'spin    ','distance'/),attributes,reshape((/4,8,1,20/),(/2,2/)))
      IF( hybrid%l_calhf ) THEN
         WRITE (16,FMT=7901) js,inDen%iter,1000*SQRT(ABS(dist(js)/cell%vol))
         WRITE ( 6,FMT=7901) js,inDen%iter,1000*SQRT(ABS(dist(js)/cell%vol))
      ELSE
         WRITE (16,FMT=7900) js,inDen%iter,1000*SQRT(ABS(dist(js)/cell%vol))
         WRITE ( 6,FMT=7900) js,inDen%iter,1000*SQRT(ABS(dist(js)/cell%vol))
      END IF
   END DO
   IF (noco%l_noco) dist(6) = CPP_BLAS_sdot((nmap-2*nmaph), fsm(nmaph*2+1),1,fmMet(nmaph*2+1),1)
   IF (noco%l_noco) WRITE (6,FMT=7900) 3,inDen%iter,1000*SQRT(ABS(dist(6)/cell%vol))

   !calculate the distance of total charge and spin density
   !|rho/m(o) - rho/m(i)| = |rh1(o) -rh1(i)|+ |rh2(o) -rh2(i)| +/_
   !                        +/_2<rh2(o) -rh2(i)|rh1(o) -rh1(i)>
   IF (input%jspins.EQ.2) THEN
      dist(3) = CPP_BLAS_sdot(nmaph,fsm,1,fmMet(nmaph+1),1)
      dist(4) = dist(1) + dist(2) + 2.0e0*dist(3)
      dist(5) = dist(1) + dist(2) - 2.0e0*dist(3)
      CALL writeXMLElementFormPoly('overallChargeDensity',(/'distance'/),&
                                   (/1000*SQRT(ABS(dist(4)/cell%vol))/),reshape((/10,20/),(/1,2/)))
      CALL writeXMLElementFormPoly('spinDensity',(/'distance'/),&
                                   (/1000*SQRT(ABS(dist(5)/cell%vol))/),reshape((/19,20/),(/1,2/)))
      IF( hybrid%l_calhf ) THEN
         WRITE (16,FMT=8001) inDen%iter,1000*SQRT(ABS(dist(4)/cell%vol))
         WRITE (16,FMT=8011) inDen%iter,1000*SQRT(ABS(dist(5)/cell%vol))
         WRITE ( 6,FMT=8001) inDen%iter,1000*SQRT(ABS(dist(4)/cell%vol))
         WRITE ( 6,FMT=8011) inDen%iter,1000*SQRT(ABS(dist(5)/cell%vol))
      ELSE
         WRITE (16,FMT=8000) inDen%iter,1000*SQRT(ABS(dist(4)/cell%vol))
         WRITE (16,FMT=8010) inDen%iter,1000*SQRT(ABS(dist(5)/cell%vol))
         WRITE ( 6,FMT=8000) inDen%iter,1000*SQRT(ABS(dist(4)/cell%vol))
         WRITE ( 6,FMT=8010) inDen%iter,1000*SQRT(ABS(dist(5)/cell%vol))
      END IF

      !dist/vol should always be >= 0 ,
      !but for dist=0 numerically you might obtain dist/vol < 0
      !(e.g. when calculating non-magnetic systems with jspins=2).
   END IF
   results%last_distance=maxval(1000*SQRT(ABS(dist/cell%vol)))
   DEALLOCATE (sm,fsm,smMet,fmMet)
   CALL closeXMLElement('densityConvergence')

   !fix charge of the new density
   CALL qfix(stars,atoms,sym,vacuum, sphhar,input,cell,oneD,&
             inDen%pw,inDen%vacxy,inDen%mt,inDen%vacz,.FALSE.,.false., fix)

   IF(atoms%n_u.NE.n_u_keep) THEN
      inDen%mmpMat = n_mmpTemp
   END IF

   atoms%n_u=n_u_keep

   IF(vacuum%nvac.EQ.1) THEN
      inDen%vacz(:,2,:) = inDen%vacz(:,1,:)
      IF (sym%invs) THEN
         inDen%vacxy(:,:,2,:) = CONJG(inDen%vacxy(:,:,1,:))
      ELSE
         inDen%vacxy(:,:,2,:) = inDen%vacxy(:,:,1,:)
      END IF
   END IF

   IF (atoms%n_u > 0) THEN
      IF (.NOT.l_densityMatrixPresent) THEN
         inDen%mmpMat(:,:,:,:) = outDen%mmpMat(:,:,:,:)
         CALL resetBroydenHistory()
      END IF
   ENDIF

   !write out mixed density
   CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                     1,results%last_distance,results%ef,.TRUE.,inDen)

   inDen%iter = inDen%iter + 1

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
