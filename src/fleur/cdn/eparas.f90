!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eparas
  !***********************************************************************
  ! Calculates qlo, enerlo and sqlo, which are needed to determine the
  ! new energy parameters.
  ! Philipp Kurz 99/04
  !***********************************************************************
  ! also the 'normal' energy parameters are now included...
  !
  ! if (l_mcd) then mcd contains mcd spectrum: first index = polarization
  ! second = core level ; third = band index                  gb.2001
  ! corrected to work also for multiple LO's of same l at the same atom
  !                                                           gb.2005
  !*************** ABBREVIATIONS *****************************************
  ! qlo     : charge density of one local orbital at the current k-point
  ! sqlo    : qlo integrated over the Brillouin zone
  ! enerlo  : qlo*energy integrated over the Brillouin zone
  !***********************************************************************
  !
CONTAINS
  SUBROUTINE eparas(jsp,atoms,banddos,noccbd,ev_list,fmpi,ikpt,ne,we,eig,skip_t,l_evp,eigVecCoeffs,&
                    usdus,dos,mcd)
    USE m_types
    use m_types_dos
    use m_types_mcd
    IMPLICIT NONE
    TYPE(t_usdus),         INTENT(IN)    :: usdus
    TYPE(t_mpi),           INTENT(IN)    :: fmpi
    TYPE(t_atoms),         INTENT(IN)    :: atoms
    TYPE(t_banddos),       INTENT(IN)    :: banddos
    TYPE(t_eigVecCoeffs),  INTENT(IN)    :: eigVecCoeffs
    TYPE(t_dos),           INTENT(INOUT) :: dos
    TYPE(t_mcd), OPTIONAL, INTENT(INOUT) :: mcd
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: noccbd,jsp
    INTEGER, INTENT (IN) :: ne,ikpt  ,skip_t
    LOGICAL, INTENT (IN) :: l_evp
    INTEGER, INTENT (IN) :: ev_list(noccbd)
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN)  :: eig(:)!(input%neig),
    REAL,    INTENT (IN)  :: we(noccbd)

    !     ..
    !     .. Local Scalars ..
    INTEGER i,l,lo,lop ,natom,nn,ntyp,m,n_dos
    INTEGER nt1,nt2,lm,n,ll1,ipol,icore,index
    REAL fac
    COMPLEX suma,sumb,sumab,sumba,sumaa,sumbb
    !     ..
    !     .. Local Arrays ..
    REAL qlo(noccbd,atoms%nlod,atoms%nlod,atoms%ntype)
    REAL qaclo(noccbd,atoms%nlod,atoms%ntype),qbclo(noccbd,atoms%nlod,atoms%ntype)
    LOGICAL atomTypeCovered(atoms%ntype)
    !     ..
    !
    !---> initialize ener, sqal, enerlo and sqlo on first call
    !

    IF ((ikpt.LE.fmpi%isize).AND..NOT.l_evp) THEN
       dos%qal(:,:,:,ikpt,jsp) = 0.0
    END IF

    atomTypeCovered(:) = .FALSE.
    !--->    l-decomposed density for each occupied state
    !
    !         DO 140 i = (skip_t+1),ne    ! this I need for all states
    DO i = 1,ne              ! skip in next loop
      DO n_dos = 1,size(banddos%dos_typelist)
        n=banddos%dos_typelist(n_dos)
        atomTypeCovered(n) = .TRUE.
        fac = 1./atoms%neq(n)
        nt1 = atoms%firstAtom(n)
        nt2 = nt1 + atoms%neq(n) - 1
        DO l = 0,3
          ll1 = l* (l+1)
          DO m = -l,l
            lm = ll1 + m
            IF (banddos%l_mcd) THEN
              sumaa = CMPLX(0.,0.) ; sumab = CMPLX(0.,0.)
              sumbb = CMPLX(0.,0.) ; sumba = CMPLX(0.,0.)
              DO natom = nt1,nt2
                sumaa = sumaa + eigVecCoeffs%abcof(i,lm,1,natom,jsp)*CONJG(eigVecCoeffs%abcof(i,lm,1,natom,jsp))
                sumbb = sumbb + eigVecCoeffs%abcof(i,lm,2,natom,jsp)*CONJG(eigVecCoeffs%abcof(i,lm,2,natom,jsp))
                sumab = sumab + eigVecCoeffs%abcof(i,lm,1,natom,jsp) *CONJG(eigVecCoeffs%abcof(i,lm,2,natom,jsp))
                sumba = sumba + eigVecCoeffs%abcof(i,lm,2,natom,jsp) *CONJG(eigVecCoeffs%abcof(i,lm,1,natom,jsp))
              ENDDO
              DO icore = 1, mcd%ncore(n)
                DO ipol = 1, 3
                  index = 3*(n_dos-1) + ipol
                  mcd%mcd(index,icore,ev_list(i),ikpt,jsp) = mcd%mcd(index,icore,ev_list(i),ikpt,jsp) + fac * &
                     (sumaa * CONJG(mcd%m_mcd(icore,lm+1,index,1))*mcd%m_mcd(icore,lm+1,index,1) + &
                      sumbb * CONJG(mcd%m_mcd(icore,lm+1,index,2))*mcd%m_mcd(icore,lm+1,index,2) + &
                      sumab * CONJG(mcd%m_mcd(icore,lm+1,index,2))*mcd%m_mcd(icore,lm+1,index,1) + &
                      sumba * CONJG(mcd%m_mcd(icore,lm+1,index,1))*mcd%m_mcd(icore,lm+1,index,2))
                ENDDO
              ENDDO
            ENDIF     ! end MCD
          ENDDO
        ENDDO


        
      ENDDO
   ENDDO



     !
    !--->    perform Brillouin zone integration and summation over the
    !--->    bands in order to determine the energy parameters for each
    !--->    atom and angular momentum
    !
    




  END SUBROUTINE eparas
END MODULE m_eparas
