!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_jDOS

   !--------------------------------------------------------------------
   ! Calculate the decomposition into the total angular momentum states
   ! characterized by j= l+-1/2 using the Clebsch Gordan coefficients
   !--------------------------------------------------------------------
   USE m_types
   USE m_clebsch

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE jDOS_comp(ikpt,noccbd,ev_list,atoms,input,usdus,denCoeffsOffdiag,eigVecCoeffs,jDOS)

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_usdus),             INTENT(IN)     :: usdus
      TYPE(t_denCoeffsOffdiag),  INTENT(IN)     :: denCoeffsOffdiag
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      INTEGER,                   INTENT(IN)     :: ikpt
      INTEGER,                   INTENT(IN)     :: noccbd
      INTEGER,                   INTENT(IN)     :: ev_list(:)
      TYPE(t_jDOS),              INTENT(INOUT)  :: jDOS

      INTEGER, PARAMETER :: lmax = 3 !Maximum l considered in j decomposition

      INTEGER :: iType,iBand,nn,natom,l,jj,lmup,lmdown,j_ind,spin
      REAL    :: j,mj,mup,mdown,facup,facdown,c(6),tot
      COMPLEX :: aup,bup,adown,bdown

      DO iType = 1, atoms%ntype
         DO nn =1, atoms%neq(iType)
            natom = SUM(atoms%neq(:iType-1)) + nn
            DO iBand = 1, noccbd
               j_ind = 0
               c = 0.0
               tot = 0.0
               DO l = 1, lmax
                  DO jj = 1, 2
                     j_ind = j_ind+1
                     ! j = l +- 1/2
                     j = l + (jj-1.5)
                     mj = -j
                     DO WHILE(mj <= j)
                        !mj = -l-+1/2, .... , l+-1/2

                        mup   = mj - 0.5
                        mdown = mj + 0.5

                        IF(input%jspins.EQ.1) THEN
                           mdown = mdown * (-1)
                           spin = 1
                        ELSE
                           spin = 2
                        ENDIF

                        IF(ABS(mup) <= l) THEN
                           lmup   = l*(l+1) + INT(mup)
                           facup = clebsch(REAL(l),0.5,mup,0.5,j,mj)
                           aup   = facup   * eigVecCoeffs%acof(iBand,lmup  ,natom,1)
                           bup   = facup   * eigVecCoeffs%bcof(iBand,lmup  ,natom,1)
                        ELSE
                           aup = 0.0
                           bup = 0.0
                        ENDIF

                        IF(ABS(mdown) <= l) THEN
                           lmdown = l*(l+1) + INT(mdown)
                           facdown = clebsch(REAL(l),0.5,mdown,-0.5,j,mj)
                           adown = facdown * eigVecCoeffs%acof(iBand,lmdown,natom,spin)
                           bdown = facdown * eigVecCoeffs%bcof(iBand,lmdown,natom,spin)
                        ELSE
                           adown = 0.0
                           bdown = 0.0
                        ENDIF

                        !c := norm of facup |up> + facdown |down>
                        !We have to write it out explicitely because
                        !of the offdiagonal scalar products that appear
                        c(j_ind) = c(j_ind) + &
                                  +        aup  *CONJG(aup)   &
                                  +        adown*CONJG(adown) &
                                  +        bup  *CONJG(bup)    * usdus%ddn(l,iType,1) &
                                  +        bdown*CONJG(bdown)  * usdus%ddn(l,iType,spin) &
                                  + 2*REAL(aup  *CONJG(adown)) * denCoeffsOffdiag%uu21n(l,iType) &
                                  + 2*REAL(bup  *CONJG(bdown)) * denCoeffsOffdiag%dd21n(l,iType) &
                                  + 2*REAL(aup  *CONJG(bdown)) * denCoeffsOffdiag%ud21n(l,iType) &
                                  + 2*REAL(adown*CONJG(bup))   * denCoeffsOffdiag%du21n(l,iType)
                        mj = mj + 1
                     ENDDO
                     tot = tot + c(j_ind)
                  ENDDO
               ENDDO
               !TODO: LOs
               j_ind=0
               DO l = 1, 3
                  DO jj = 1, 2
                     j_ind = j_ind+1
                     jDOS%comp(ev_list(iBand),l,jj,natom,ikpt) = c(j_ind)*100.0/tot
                     jDOS%qmtp(ev_list(iBand),natom,ikpt) = 100.0*tot
                  ENDDO
               ENDDO
            ENDDO

         ENDDO
      ENDDO

   END SUBROUTINE jDOS_comp
END MODULE m_jDOS