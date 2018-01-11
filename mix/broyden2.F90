!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! This is a continuously restartable implementation of the original
!!! broyden subroutine. It constructs the required u and v vectors for
!!! the relevant iterations on the fly in terms of linear combinations
!!! of the deltaN_i and deltaF_i vectors. It is slower than the original
!!! broyden subroutine so if the continuous restartability is not needed
!!! the old routine should be used.
!!!
!!!                                   GM'2018
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_broyden2
  USE m_juDFT
  !################################################################
  !     IMIX = 3 : BROYDEN'S FIRST METHOD
  !     IMIX = 5 : BROYDEN'S SECOND METHOD
  !     IMIX = 7 : GENERALIZED ANDERSEN METHOD
  !     sm   : input charge density of iteration m
  !            afterwards update rho(m+1)
  !     fm   : output minus input charge density of iteration m
  !################################################################
CONTAINS
  SUBROUTINE broyden2(cell,stars,atoms,vacuum,sphhar,input,noco,oneD,sym,&
                      hybrid,mmap,nmaph,mapmt,mapvac2,nmap,fm,sm,lpot)

#include"cpp_double.h"

    USE m_metric
    USE m_types
    USE m_broyd_io

    IMPLICIT NONE

    TYPE(t_oneD),INTENT(IN)    :: oneD
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)  :: vacuum
    TYPE(t_noco),INTENT(IN)    :: noco
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_cell),INTENT(IN)    :: cell
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_hybrid),INTENT(IN)  :: hybrid

    ! Scalar Arguments
    INTEGER, INTENT (IN)        :: mmap,nmap
    INTEGER, INTENT (IN)        :: mapmt,mapvac2
    LOGICAL,OPTIONAL,INTENT(IN) :: lpot

    ! Array Arguments
    REAL,    INTENT (IN)    :: fm(nmap) 
    REAL,    INTENT (INOUT) :: sm(nmap)

    ! Local Scalars
    INTEGER         :: i,j,it,k,nit,iread,nmaph, mit, historyLength
    REAL            :: vFMetProd,alphan,coeff,vNorm
    LOGICAL         :: l_pot, l_exist

    ! Local Arrays
    REAL, ALLOCATABLE :: fm1(:),sm1(:),uVec(:),vVec(:)
    REAL, ALLOCATABLE :: dNVec(:),dFVec(:),dNMet(:), dFMet(:), FMet(:)
    REAL, ALLOCATABLE :: deltaN_i(:), deltaF_i(:)
    REAL, ALLOCATABLE :: dNdNMat(:,:), dNdFMat(:,:), dFdNMat(:,:), dFdFMat(:,:)
    REAL, ALLOCATABLE :: dNdNLast(:), dNdFLast(:), dFdNLast(:), dFdFLast(:)
    REAL, ALLOCATABLE :: uDNTableau(:,:), vDNTableau(:,:), wDNTableau(:,:)
    REAL, ALLOCATABLE :: uDFTableau(:,:), vDFTableau(:,:), wDFTableau(:,:)

    ! External Functions
    REAL CPP_BLAS_sdot
    EXTERNAL CPP_BLAS_sdot

    ! External Subroutines
    EXTERNAL CPP_BLAS_saxpy,CPP_BLAS_sscal

    l_pot = .FALSE.
    IF (PRESENT(lpot)) l_pot = lpot

    ALLOCATE (fm1(mmap),sm1(mmap),dNVec(mmap),dFVec(mmap),uVec(mmap),vVec(mmap))
    ALLOCATE (deltaN_i(mmap),deltaF_i(mmap))
    ALLOCATE (dNdNLast(input%maxiter),dFdFLast(input%maxiter))
    ALLOCATE (dNdFLast(input%maxiter),dFdNLast(input%maxiter))

    ALLOCATE (dnMet(mmap), dFMet(mmap), FMet(mmap))

    fm1 = 0.0
    sm1 = 0.0
    dNVec = 0.0
    dFVec = 0.0
    vVec  = 0.0
    dNMet = 0.0
    dFMet = 0.0

    mit = 0
    l_exist = initBroydenHistory(input,hybrid,nmap) ! returns true if there already exists a Broyden history
    IF(.NOT.l_exist) mit = 1

    IF (mit.NE.1) THEN
       ! load input charge density (sm1) and difference of 
       ! in and out charge densities (fm1) from previous iteration (m-1)

       CALL readLastIterInAndDiffDen(hybrid,nmap,mit,alphan,sm1(:nmap),fm1(:nmap))
       IF (ABS(input%alpha-alphan) > 0.0001) THEN
          WRITE (6,*) 'mixing parameter has been changed; reset'
          WRITE (6,*) 'broyden algorithm or set alpha to',alphan
          CALL juDFT_error("mixing parameter (input) changed", calledby ="broyden")
       END IF

       ! generate F_m   - F_(m-1)  ... sm1
       !      and rho_m - rho_(m-1) .. fm1
       dNVec(1:nmap) = sm(1:nmap) - sm1(1:nmap)
       dFVec(1:nmap) = fm(1:nmap) - fm1(1:nmap)
    END IF

    ! save F_m and rho_m for next iteration
    nit = mit +1
    IF (nit > input%maxiter+1) nit = 1
    CALL writeLastIterInAndDiffDen(hybrid,nmap,nit,input%alpha,sm,fm)

    IF (mit.EQ.1) THEN 
       !     update for rho for mit=1 is straight mixing
       !     sm = sm + alpha*fm
       CALL CPP_BLAS_saxpy(nmap,input%alpha,fm,1,sm,1)
    ELSE

       CALL writeDeltaNVec(input,hybrid,nmap,mit,dNVec)
       CALL writeDeltaFVec(input,hybrid,nmap,mit,dFVec)

       ! Apply metric w to dNVec and store in dNMet:  w |dNVec>  
       CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
                   mmap,nmaph,mapmt,mapvac2,dNVec,dNMet,l_pot)

       ! Apply metric w to dFVec and store in dFMet:  w |dFVec>  
       CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
                   mmap,nmaph,mapmt,mapvac2,dFVec,dFMet,l_pot)

       ! Extend overlap matrices <delta n(i) | delta n(j)>, <delta F(i) | delta F(j)>,
       !                         <delta n(i) | delta F(j)>, <delta F(i) | delta n(j)> -start-

       iread = MIN(mit-1,input%maxiter+1)
       historyLength = iread

       dNdNLast = 0.0
       dFdFLast = 0.0
       dNdFLast = 0.0
       dFdNLast = 0.0

       DO it = 2, iread
          CALL readDeltaNVec(input,hybrid,nmap,it-mit,mit,deltaN_i)
          CALL readDeltaFVec(input,hybrid,nmap,it-mit,mit,deltaF_i)

          dNdNLast(it-1) = CPP_BLAS_sdot(nmap,deltaN_i,1,dNMet,1)
          dFdFLast(it-1) = CPP_BLAS_sdot(nmap,deltaF_i,1,dFMet,1)
          dNdFLast(it-1) = CPP_BLAS_sdot(nmap,deltaN_i,1,dFMet,1)
          dFdNLast(it-1) = CPP_BLAS_sdot(nmap,deltaF_i,1,dNMet,1)
       END DO

       dNdNLast(historyLength) = CPP_BLAS_sdot(nmap,dNVec,1,dNMet,1)
       dFdFLast(historyLength) = CPP_BLAS_sdot(nmap,dFVec,1,dFMet,1)
       dNdFLast(historyLength) = CPP_BLAS_sdot(nmap,dNVec,1,dFMet,1)
       dFdNLast(historyLength) = CPP_BLAS_sdot(nmap,dFVec,1,dNMet,1)

       CALL writeBroydenOverlapExt(input,hybrid,mit,historyLength,&
                                   dNdNLast,dFdFLast,dNdFLast,dFdNLast)

       ALLOCATE (dNdNMat(historyLength,historyLength))
       ALLOCATE (dFdFMat(historyLength,historyLength))
       ALLOCATE (dNdFMat(historyLength,historyLength))
       ALLOCATE (dFdNMat(historyLength,historyLength))

       CALL readBroydenOverlaps(input,hybrid,mit,historyLength,&
                                dNdNMat,dFdFMat,dNdFMat,dFdNMat)

       ! Extend overlap matrices <delta n(i) | delta n(j)>, <delta F(i) | delta F(j)>,
       !                         <delta n(i) | delta F(j)>, <delta F(i) | delta n(j)> -end-

       ! Construct u_i, v_i tableaus (prefactors for delta n(i) and delta F(i) vectors) -start-

       ALLOCATE (uDNTableau(historyLength,historyLength)) !first index: iteration of DN, second index: iteration of u
       ALLOCATE (vDNTableau(historyLength,historyLength)) !first index: iteration of DN, second index: iteration of v
       ALLOCATE (uDFTableau(historyLength,historyLength)) !first index: iteration of DF, second index: iteration of u
       ALLOCATE (vDFTableau(historyLength,historyLength)) !first index: iteration of DF, second index: iteration of v

       uDNTableau = 0.0
       vDNTableau = 0.0
       uDFTableau = 0.0
       vDFTableau = 0.0

       IF (input%imix.EQ.3) THEN ! Broyden's first method

          DO i = 1, historyLength
             uDNTableau(i,i) = 1.0
             uDFTableau(i,i) = input%alpha
             vDNTableau(i,i) = input%alpha

             DO j = 1, i-1
                ! Calculations for u_i
                coeff = 0.0
                DO k = 1, historyLength
                   coeff = coeff + vDNTableau(k,j)*dNdFMat(k,i)
                   coeff = coeff + vDFTableau(k,j)*dFdFMat(k,i)
                END DO
                DO k = 1, historyLength
                   uDNTableau(k,i) = uDNTableau(k,i) - coeff*uDNTableau(k,j)
                   uDFTableau(k,i) = uDFTableau(k,i) - coeff*uDFTableau(k,j)
                END DO

                ! Calculations for v_i (numerator)
                coeff = 0.0
                DO k = 1, historyLength
                   coeff = coeff + uDNTableau(k,j)*dNdNMat(k,i) ! This is in agreement with the old implementation.
                   coeff = coeff + uDFTableau(k,j)*dFdNMat(k,i)
!                   coeff = coeff + vDNTableau(k,j)*dNdNMat(k,i) ! This is in agreement with RP's thesis.
!                   coeff = coeff + vDFTableau(k,j)*dFdNMat(k,i)
                END DO
                DO k = 1, historyLength
                   vDNTableau(k,i) = vDNTableau(k,i) - coeff*vDNTableau(k,j) ! This is in agreement with the old implementation.
                   vDFTableau(k,i) = vDFTableau(k,i) - coeff*vDFTableau(k,j)
!                   vDNTableau(k,i) = vDNTableau(k,i) - coeff*uDNTableau(k,j) ! This is in agreement with RP's thesis.
!                   vDFTableau(k,i) = vDFTableau(k,i) - coeff*uDFTableau(k,j)
                END DO
             END DO

             ! Calculations for v_i (denominator)
             vNorm = dNdNMat(i,i)
             DO k = 1, historyLength
                vNorm = vNorm - uDNTableau(k,i)*dNdNMat(k,i)
                vNorm = vNorm - uDFTableau(k,i)*dFdNMat(k,i)
             END DO

             vNorm = -vNorm ! Note: This is in agreement with the old implementation.
                            !       The equations in RP's thesis don't have this sign.

             DO k = 1, historyLength
                vDNTableau(k,i) = vDNTableau(k,i) / vNorm
                vDFTableau(k,i) = vDFTableau(k,i) / vNorm
             END DO

          END DO
          
       ELSE IF (input%imix.EQ.5) THEN ! Broyden's second method

          DO i = 1, historyLength
             uDNTableau(i,i) = 1.0
             uDFTableau(i,i) = input%alpha
             vDFTableau(i,i) = 1.0 / dFdFMat(i,i)

             DO j = 1, i-1
                coeff = dFdFMat(j,i)*vDFTableau(j,j)
                DO k = 1, historyLength
                   uDNTableau(k,i) = uDNTableau(k,i) - coeff*uDNTableau(k,j)
                   uDFTableau(k,i) = uDFTableau(k,i) - coeff*uDFTableau(k,j)
                END DO
             END DO
          END DO

       ELSE IF (input%imix.EQ.7) THEN ! generalized anderson method

          ALLOCATE (wDNTableau(historyLength,historyLength)) !first index: iteration of DN, second index: iteration of w
          ALLOCATE (wDFTableau(historyLength,historyLength)) !first index: iteration of DF, second index: iteration of w

          wDNTableau = 0.0
          wDFTableau = 0.0

          DO i = 1, historyLength

             uDNTableau(i,i) = 1.0
             uDFTableau(i,i) = input%alpha
             wDFTableau(i,i) = 1.0

             DO j = 1, i-1
                ! Calculations for u_i and w_i
                coeff = 0.0
                DO k = 1, historyLength
                   coeff = coeff + vDNTableau(k,j)*dNdFMat(k,i)
                   coeff = coeff + vDFTableau(k,j)*dFdFMat(k,i)
                END DO
                DO k = 1, historyLength
                   ! Calculations for u_i
                   uDNTableau(k,i) = uDNTableau(k,i) - coeff*uDNTableau(k,j)
                   uDFTableau(k,i) = uDFTableau(k,i) - coeff*uDFTableau(k,j)

                   ! Calculations for w_i
                   wDNTableau(k,i) = wDNTableau(k,i) - coeff*wDNTableau(k,j)
                   wDFTableau(k,i) = wDFTableau(k,i) - coeff*wDFTableau(k,j)
                END DO

             END DO

             ! Calculations for v_i
             vNorm = 0.0
             DO k = 1, historyLength
                vNorm = vNorm + wDNTableau(k,i)*dNdFMat(k,i)
                vNorm = vNorm + wDFTableau(k,i)*dFdFMat(k,i)
             END DO

             DO k = 1, historyLength
                vDNTableau(k,i) = wDNTableau(k,i) / vNorm
                vDFTableau(k,i) = wDFTableau(k,i) / vNorm
             END DO

          END DO

          DEALLOCATE (wDNTableau, wDFTableau)

       END IF

       ! Construct u_i, v_i tableau (prefactors for delta n(i) and delta F(i) vectors) -end-

       ! Construct the final uVec and vVec -start-
       uVec = 0.0
       vVec = 0.0

       DO it = 2, iread
          CALL readDeltaNVec(input,hybrid,nmap,it-mit,mit,deltaN_i)
          CALL readDeltaFVec(input,hybrid,nmap,it-mit,mit,deltaF_i)

          DO k = 1, nmap
             uVec(k) = uVec(k) + uDNTableau(it-1,historyLength)*deltaN_i(k)
             uVec(k) = uVec(k) + uDFTableau(it-1,historyLength)*deltaF_i(k)
             vVec(k) = vVec(k) + vDNTableau(it-1,historyLength)*deltaN_i(k)
             vVec(k) = vVec(k) + vDFTableau(it-1,historyLength)*deltaF_i(k)
          END DO
       END DO

       DO k = 1, nmap
          uVec(k) = uVec(k) + uDNTableau(historyLength,historyLength)*dNVec(k)
          uVec(k) = uVec(k) + uDFTableau(historyLength,historyLength)*dFVec(k)
          vVec(k) = vVec(k) + vDNTableau(historyLength,historyLength)*dNVec(k)
          vVec(k) = vVec(k) + vDFTableau(historyLength,historyLength)*dFVec(k)
       END DO
       ! Construct the final uVec and vVec -end-

       ! Apply metric w to dFVec and store in FMet:  w |dFVec>
       CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
                   mmap,nmaph,mapmt,mapvac2,fm,FMet,l_pot)

       vFMetProd = CPP_BLAS_sdot(nmap,vVec,1,FMet,1)

       ! update rho(m+1)
       ! calculate sm(:) = (1.0-vFMetProd)*uVec(:) + sm
       CALL CPP_BLAS_saxpy(nmap,1.0-vFMetProd,uVec,1,sm,1)
    END IF

    DEALLOCATE (fm1,sm1,dNVec,dFVec,vVec)

  END SUBROUTINE broyden2
END MODULE m_broyden2
