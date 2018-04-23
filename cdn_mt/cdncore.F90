!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdncore

CONTAINS

SUBROUTINE cdncore(results,mpi,dimension,oneD,sliceplot,input,vacuum,noco,sym,&
                   stars,cell,sphhar,atoms,vTot,outDen,moments)

   USE m_constants
   USE m_cdn_io
   USE m_cdnovlp
   USE m_cored
   USE m_coredr
   USE m_types
   USE m_xmlOutput
   USE m_magMoms
   USE m_orbMagMoms
#ifdef CPP_MPI
   USE m_mpi_bc_coreden
#endif

   IMPLICIT NONE

   TYPE(t_results),INTENT(INOUT)    :: results
   TYPE(t_mpi),INTENT(IN)           :: mpi
   TYPE(t_dimension),INTENT(IN)     :: dimension
   TYPE(t_oneD),INTENT(IN)          :: oneD
   TYPE(t_sliceplot),INTENT(IN)     :: sliceplot
   TYPE(t_input),INTENT(IN)         :: input
   TYPE(t_vacuum),INTENT(IN)        :: vacuum
   TYPE(t_noco),INTENT(IN)          :: noco
   TYPE(t_sym),INTENT(IN)           :: sym
   TYPE(t_stars),INTENT(IN)         :: stars
   TYPE(t_cell),INTENT(IN)          :: cell
   TYPE(t_sphhar),INTENT(IN)        :: sphhar
   TYPE(t_atoms),INTENT(IN)         :: atoms
   TYPE(t_potden),INTENT(IN)        :: vTot
   TYPE(t_potden),INTENT(INOUT)     :: outDen
   TYPE(t_moments),INTENT(INOUT)    :: moments

   INTEGER                          :: jspin, n, iType
   REAL                             :: seig, rhoint, momint
   LOGICAL, PARAMETER               :: l_st=.FALSE.

   REAL                             :: rh(dimension%msh,atoms%ntype,dimension%jspd)
   REAL                             :: qint(atoms%ntype,dimension%jspd)
   REAL                             :: tec(atoms%ntype,DIMENSION%jspd)
   REAL                             :: rhTemp(dimension%msh,atoms%ntype,dimension%jspd)

   results%seigc = 0.0
   IF (mpi%irank.EQ.0) THEN
      DO jspin = 1,input%jspins
         DO n = 1,atoms%ntype
            moments%svdn(n,jspin) = outDen%mt(1,0,n,jspin) / (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
         END DO
      END DO
   END IF

   IF (input%kcrel.EQ.0) THEN
      ! Generate input file ecore for subsequent GW calculation
      ! 11.2.2004 Arno Schindlmayr
      IF ((input%gw.eq.1 .or. input%gw.eq.3).AND.(mpi%irank.EQ.0)) THEN
         OPEN (15,file='ecore',status='unknown', action='write',form='unformatted')
      END IF

      rh = 0.0
      tec = 0.0
      qint = 0.0
      IF (input%frcor) THEN
         IF (mpi%irank.EQ.0) THEN
            CALL readCoreDensity(input,atoms,dimension,rh,tec,qint)
         END IF
#ifdef CPP_MPI
         CALL mpi_bc_coreDen(mpi,atoms,input,dimension,rh,tec,qint)
#endif
      END IF
   END IF

   IF (.NOT.sliceplot%slice) THEN
      !add in core density
      IF (mpi%irank.EQ.0) THEN
         IF (input%kcrel.EQ.0) THEN
            DO jspin = 1,input%jspins
               CALL cored(input,jspin,atoms,outDen%mt,dimension,sphhar,vTot%mt(:,0,:,jspin), qint,rh,tec,seig)
               rhTemp(:,:,jspin) = rh(:,:,jspin)
               results%seigc = results%seigc + seig
            END DO
         ELSE
            CALL coredr(input,atoms,seig, outDen%mt,dimension,sphhar,vTot%mt(:,0,:,:),qint,rh)
            results%seigc = results%seigc + seig
         END IF
      END IF
      DO jspin = 1,input%jspins
         IF (mpi%irank.EQ.0) THEN
            DO n = 1,atoms%ntype
               moments%stdn(n,jspin) = outDen%mt(1,0,n,jspin) / (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
            END DO
         END IF
         IF ((noco%l_noco).AND.(mpi%irank.EQ.0)) THEN
            IF (jspin.EQ.2) THEN
               !pk non-collinear (start)
               !add the coretail-charge to the constant interstitial
               !charge (star 0), taking into account the direction of
               !magnetisation of this atom
               DO iType = 1,atoms%ntype
                  rhoint = (qint(iType,1) + qint(iType,2)) /cell%volint/input%jspins/2.0
                  momint = (qint(iType,1) - qint(iType,2)) /cell%volint/input%jspins/2.0
                  !rho_11
                  outDen%pw(1,1) = outDen%pw(1,1) + rhoint + momint*cos(noco%beta(iType))
                  !rho_22
                  outDen%pw(1,2) = outDen%pw(1,2) + rhoint - momint*cos(noco%beta(iType))
                  !real part rho_21
                  outDen%pw(1,3) = outDen%pw(1,3) + cmplx(0.5*momint *cos(noco%alph(iType))*sin(noco%beta(iType)),0.0)
                  !imaginary part rho_21
                  outDen%pw(1,3) = outDen%pw(1,3) + cmplx(0.0,-0.5*momint *sin(noco%alph(iType))*sin(noco%beta(iType)))
               END DO
               !pk non-collinear (end)
            END IF
         ELSE
            IF (input%ctail) THEN
               !+gu hope this works as well
               CALL cdnovlp(mpi,sphhar,stars,atoms,sym,dimension,vacuum,&
                            cell,input,oneD,l_st,jspin,rh(:,:,jspin),&
                            outDen%pw,outDen%vacxy,outDen%mt,outDen%vacz)
            ELSE IF (mpi%irank.EQ.0) THEN
               DO iType = 1,atoms%ntype
                  outDen%pw(1,jspin) = outDen%pw(1,jspin) + qint(iType,jspin)/input%jspins/cell%volint
               END DO
            END IF
         END IF
      END DO
   END IF

   IF (input%kcrel.EQ.0) THEN
      IF (mpi%irank.EQ.0) THEN
         CALL writeCoreDensity(input,atoms,dimension,rhTemp,tec,qint)
      END IF
      IF ((input%gw.eq.1 .or. input%gw.eq.3).AND.(mpi%irank.EQ.0)) CLOSE(15)
   END IF

END SUBROUTINE cdncore

END MODULE m_cdncore
