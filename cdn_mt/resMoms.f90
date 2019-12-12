! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_resMoms

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates and writes out intraatomic electric and magnetic dipole
! moments resolved with respect to their orbital (angular momentum) origins.
!
!                                           GM'2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE resMoms(sym,input,atoms,sphhar,noco,den,rhoLRes)

   USE m_constants
   USE m_types
   USE m_juDFT
   USE m_magDiMom

   IMPLICIT NONE
   TYPE(t_sym),           INTENT(IN)    :: sym
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_sphhar),        INTENT(IN)    :: sphhar
   TYPE(t_noco),          INTENT(IN)    :: noco
   TYPE(t_potden),        INTENT(IN)    :: den
   REAL,                  INTENT(IN)    :: rhoLRes(:,0:,0:,:,:)

   REAL,    ALLOCATABLE :: rhoTemp(:,:,:,:)

   REAL    :: t_op(3,atoms%ntype), elecDip(3,atoms%ntype)
   REAL    :: res_T_op(3,atoms%ntype,0:(atoms%lmaxd*(atoms%lmaxd+1))/2+atoms%lmaxd)
   REAL    :: resElecDip(3,atoms%ntype,0:(atoms%lmaxd*(atoms%lmaxd+1))/2+atoms%lmaxd)

   INTEGER :: iType, l, lp, llp

   IF(input%jspins.EQ.1) RETURN
   IF(.NOT.noco%l_noco) RETURN

   t_op = 0.0
   res_T_op = 0.0
   elecDip = 0.0
   resElecDip = 0.0
   ALLOCATE(rhoTemp(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,4))

   rhoTemp = 0.0

   rhoTemp(:,:,:,1) = den%mt(:,:,:,1)
   rhoTemp(:,:,:,2) = den%mt(:,:,:,2)

   IF (noco%l_mperp) THEN
      rhoTemp(:,:,:,3) = den%mt(:,:,:,3)
      rhoTemp(:,:,:,4) = den%mt(:,:,:,4)
!      WRITE(5000,'(f15.8)') den%mt(:,:,:,3)
!      WRITE(5000,'(f15.8)') den%mt(:,:,:,4)
   END IF

   CALL magDiMom(sym,input,atoms,sphhar,noco,noco%l_mperp,rhoTemp,t_op,elecDip)

   DO l = 0, atoms%lmaxd
      DO lp = 0, l
         llp = (l* (l+1))/2 + lp
         rhoTemp = 0.0
         rhoTemp(:,:,:,1) = rhoLRes(:,:,llp,:,1)
         rhoTemp(:,:,:,2) = rhoLRes(:,:,llp,:,2)
         rhoTemp(:,:,:,3) = rhoLRes(:,:,llp,:,3)
         rhoTemp(:,:,:,4) = rhoLRes(:,:,llp,:,4)
         CALL magDiMom(sym,input,atoms,sphhar,noco,noco%l_mperp,rhoTemp,res_T_op(:,:,llp),resElecDip(:,:,llp))
      END DO
   END DO

   DO iType = 1, atoms%ntype
      WRITE(6,*) 'Intraatomic electric and magnetic dipole moments for atom type ', iType,':'
      WRITE(6,'(a)')        '             lowL  largeL      p_x            p_y            p_z            t_x            t_y            t_z'
      WRITE(6,'(a,6f15.8)') 'Overall:              ', elecDip(:,iType), t_op(:,iType)
      DO l = 0, atoms%lmax(iType)
         DO lp = 0, l
            llp = (l* (l+1))/2 + lp
            IF(ALL(ABS(res_T_op(:,iType,llp)).LT.1.0e-8).AND.&
               ALL(ABS(resElecDip(:,iType,llp)).LT.1.0e-8)) CYCLE
            WRITE(6,'(a,2i6,6f15.8)') '          ', lp, l, resElecDip(:,iType,llp),res_T_op(:,iType,llp)
         END DO
      END DO
   END DO

END SUBROUTINE resMoms

END MODULE m_resMoms
