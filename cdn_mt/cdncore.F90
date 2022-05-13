!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdncore

CONTAINS

SUBROUTINE cdncore(fmpi ,input,vacuum,noco,nococonv,sym,&
                   stars,cell,sphhar,atoms,vTot,outDen,moments,results, EnergyDen)

   USE m_constants
   USE m_judft
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


   TYPE(t_mpi),        INTENT(IN)              :: fmpi

    
   TYPE(t_input),      INTENT(IN)              :: input
   TYPE(t_vacuum),     INTENT(IN)              :: vacuum
   TYPE(t_noco),       INTENT(IN)              :: noco
   TYPE(t_nococonv),   INTENT(IN)              :: nococonv
   TYPE(t_sym),        INTENT(IN)              :: sym
   TYPE(t_stars),      INTENT(IN)              :: stars
   TYPE(t_cell),       INTENT(IN)              :: cell
   TYPE(t_sphhar),     INTENT(IN)              :: sphhar
   TYPE(t_atoms),      INTENT(IN)              :: atoms
   TYPE(t_potden),     INTENT(IN)              :: vTot
   TYPE(t_potden),     INTENT(INOUT)           :: outDen
   TYPE(t_moments),    INTENT(INOUT)           :: moments
   TYPE(t_results),    INTENT(INOUT)           :: results
   TYPE(t_potden),     INTENT(INOUT), OPTIONAL :: EnergyDen

   INTEGER                          :: jspin, n, iType, ierr
   REAL                             :: seig, rhoint, momint
   LOGICAL, PARAMETER               :: l_st=.FALSE.
   LOGICAL                          :: l_coreDenPresent

   REAL                             :: rh(atoms%msh,atoms%ntype,input%jspins)
   REAL                             :: qint(atoms%ntype,input%jspins)
   REAL                             :: tec(atoms%ntype,input%jspins)
   REAL                             :: rhTemp(atoms%msh,atoms%ntype,input%jspins)

   results%seigc = 0.0
   IF (fmpi%irank==0) THEN
      DO jspin = 1,input%jspins
         DO n = 1,atoms%ntype
            moments%svdn(n,jspin) = outDen%mt(1,0,n,jspin) / (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
         END DO
      END DO
   END IF

   l_CoreDenPresent = .FALSE.
   IF (input%kcrel==0) THEN
      ! Generate input file ecore for subsequent GW calculation
      ! 11.2.2004 Arno Schindlmayr
      IF ((input%gw==1 .or. input%gw==3).AND.(fmpi%irank==0)) THEN
         OPEN (15,file='ecore',status='unknown', action='write',form='unformatted')
      END IF

      rh = 0.0
      tec = 0.0
      qint = 0.0
      IF (input%frcor) THEN
         IF (fmpi%irank==0) THEN
            IF(isCoreDensityPresent()) THEN
               CALL readCoreDensity(input,atoms,rh,tec,qint)
               l_coreDenPresent = .TRUE.
            END IF
         END IF
#ifdef CPP_MPI
         CALL MPI_BCAST(l_CoreDenPresent,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
         CALL mpi_bc_coreDen(fmpi,atoms,input,rh,tec,qint)
#endif
      END IF
   END IF

   !add in core density
   IF (fmpi%irank==0) THEN
      IF (input%kcrel==0) THEN
         DO jspin = 1,input%jspins
            IF(PRESENT(EnergyDen)) THEN
               CALL cored(input,jspin,atoms,outDen%mt,sphhar,l_CoreDenPresent,vTot%mt(:,0,:,jspin), qint,rh ,tec,seig, EnergyDen%mt)
            ELSE
               CALL cored(input,jspin,atoms,outDen%mt,sphhar,l_CoreDenPresent,vTot%mt(:,0,:,jspin), qint,rh ,tec,seig)
            ENDIF

            rhTemp(:,:,jspin) = rh(:,:,jspin)
            results%seigc = results%seigc + seig
         END DO
      ELSE
         IF(PRESENT(EnergyDen)) call juDFT_error("Energyden not implemented for relativistic core calculations")
         CALL coredr(input,atoms,seig, outDen%mt,sphhar,vTot%mt(:,0,:,:),qint,rh)
         results%seigc = results%seigc + seig
      END IF
   END IF
   DO jspin = 1,input%jspins
      IF (fmpi%irank==0) THEN
         DO n = 1,atoms%ntype
            moments%stdn(n,jspin) = outDen%mt(1,0,n,jspin) / (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
         END DO
      END IF
      IF ((noco%l_noco.and..not.input%ctail).AND.(fmpi%irank==0)) THEN
         IF (jspin==2) THEN

            IF(PRESENT(EnergyDen)) call juDFT_error("Energyden not implemented for noco")
            !pk non-collinear (start)
            !add the coretail-charge to the constant interstitial
            !charge (star 0), taking into account the direction of
            !magnetisation of this atom
            DO iType = 1,atoms%ntype
               rhoint = (qint(iType,1) + qint(iType,2)) /(cell%volint * input%jspins * 2.0)
               momint = (qint(iType,1) - qint(iType,2)) /(cell%volint * input%jspins * 2.0)
               !rho_11
               outDen%pw(1,1) = outDen%pw(1,1) + rhoint + momint*cos(nococonv%beta(iType))
               !rho_22
               outDen%pw(1,2) = outDen%pw(1,2) + rhoint - momint*cos(nococonv%beta(iType))
               !real part rho_21
               outDen%pw(1,3) = outDen%pw(1,3) + cmplx( -1*momint *cos(nococonv%alph(iType))*sin(nococonv%beta(iType)),&
               !imaginary part rho_21
                                                          momint *sin(nococonv%alph(iType))*sin(nococonv%beta(iType)))
               !TODO: Should be +,+ for no magic minus.
            END DO
            !pk non-collinear (end)
         END IF
      END IF
   END DO
   DO jspin = 1,input%jspins
      IF (input%ctail) THEN
         IF (noco%l_noco.and.jspin==1) THEN
            rh(:,:,1)=(rh(:,:,1)+rh(:,:,2))/2.
            rh(:,:,2)=rh(:,:,1)
         END IF
         IF(PRESENT(EnergyDen)) call juDFT_error("Energyden not implemented for ctail")
            !+gu hope this works as well
            CALL cdnovlp(fmpi,sphhar,stars,atoms,sym,vacuum,&
                         cell,input ,l_st,jspin,rh(:,:,jspin),&
                         outDen%pw,outDen%vacxy,outDen%mt,outDen%vacz,vTot%pw_w,vTot%mt)
      ELSE IF ((fmpi%irank==0).AND.(.NOT.noco%l_noco)) THEN
         DO iType = 1,atoms%ntype
            outDen%pw(1,jspin) = outDen%pw(1,jspin) + qint(iType,jspin) / (input%jspins * cell%volint)
         END DO
      END IF
   END DO

   IF (input%kcrel==0) THEN
      IF (fmpi%irank==0) THEN
         CALL writeCoreDensity(input,atoms,rhTemp,tec,qint)
         outDen%mtCore = rhTemp
         outDen%tec = tec
         outDen%qint = qint
      END IF
      IF ((input%gw==1 .or. input%gw==3).AND.(fmpi%irank==0)) CLOSE(15)
   END IF

END SUBROUTINE cdncore

END MODULE m_cdncore
