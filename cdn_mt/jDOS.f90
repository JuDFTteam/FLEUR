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

      INTEGER :: iType,iBand,nn,natom,l,jj,lmup,lmdown
      REAL    :: j,mj,mup,mdown,facup,facdown,c
      COMPLEX :: aup,bup,adown,bdown

      DO iType = 1, atoms%ntype
         DO nn =1, atoms%neq(iType)
            natom = SUM(atoms%neq(:iType-1)) + nn
            DO l = 0, lmax
               DO jj = -1, 1, 2
                  ! j = l +- 1/2
                  j = (2*l+jj)/2.0
                  mj = -j
                  DO WHILE(mj <= j)
                     !mj = -l-+1/2, .... , l+-1/2

                     mup   = mj - 0.5
                     mdown = mj + 0.5

                     lmup   = l*(l+1) + INT(mup)
                     lmdown = l*(l+1) + INT(mdown)

                     IF(ABS(mup) <= l) THEN
                        facup = clebsch(REAL(l),0.5,mup,0.5,j,mj)
                     ELSE
                        facup = 0.0
                     ENDIF

                     IF(ABS(mdown) <= l) THEN
                        facdown = clebsch(REAL(l),0.5,mdown,-0.5,j,mj)
                     ELSE
                        facdown = 0.0
                     ENDIF

                     DO iBand = 1, noccbd
                        aup   = facup   * eigVecCoeffs%acof(iBand,lmup  ,natom,1)
                        adown = facdown * eigVecCoeffs%acof(iBand,lmdown,natom,2)
                        bup   = facup   * eigVecCoeffs%bcof(iBand,lmup  ,natom,1)
                        bdown = facdown * eigVecCoeffs%bcof(iBand,lmdown,natom,2)

                        !c := norm of facup |up> + facdown |down>
                        !We have to write it out explicitely because
                        !of the offdiagonal scalar products that appear
                        c =         aup  *CONJG(aup)   &
                           +        adown*CONJG(adown) &
                           +        bup  *CONJG(bup)    * usdus%ddn(l,iType,1) &
                           +        bdown*CONJG(bdown)  * usdus%ddn(l,iType,2) &
                           + 2*REAL(aup  *CONJG(adown)) * denCoeffsOffdiag%uu21n(l,iType) &
                           + 2*REAL(bup  *CONJG(bdown)) * denCoeffsOffdiag%dd21n(l,iType) &
                           + 2*REAL(aup  *CONJG(bdown)) * denCoeffsOffdiag%ud21n(l,iType) &
                           + 2*REAL(adown*CONJG(bup))   * denCoeffsOffdiag%du21n(l,iType)

                        !TODO: LOs

                        jDOS%comp(ev_list(iBand),l,INT(jj*0.5+0.5),natom,ikpt) = &
                           jDOS%comp(ev_list(iBand),l,INT(jj*0.5+0.5),natom,ikpt) + c

                     ENDDO
                     mj = mj + 1
                  ENDDO
               ENDDO
            ENDDO

         ENDDO
      ENDDO

   END SUBROUTINE jDOS_comp
END MODULE m_jDOS