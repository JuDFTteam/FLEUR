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
   use m_types_jDOS

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE jDOS_comp(ikpt,noccbd,ev_list,we,atoms,banddos,input,usdus,&
                        denCoeffsOffdiag,eigVecCoeffs,jDOS)

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_banddos),           INTENT(IN)     :: banddos
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_usdus),             INTENT(IN)     :: usdus
      TYPE(t_denCoeffsOffdiag),  INTENT(IN)     :: denCoeffsOffdiag
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      INTEGER,                   INTENT(IN)     :: ikpt
      INTEGER,                   INTENT(IN)     :: noccbd
      INTEGER,                   INTENT(IN)     :: ev_list(:)
      REAL,                      INTENT(IN)     :: we(:)
      TYPE(t_jDOS),              INTENT(INOUT)  :: jDOS

      INTEGER, PARAMETER :: lmax = 3 !Maximum l considered in j decomposition

      INTEGER :: n_dos
      INTEGER :: iType,iBand,nn,natom,l,jj,j_ind,lmup,lmdown,spin,ilo,ilop
      REAL    :: j,mj,mup,mdown
      REAL    :: facup,facdown,summed,cf
      COMPLEX :: aup,bup,cup,adown,bdown,cdown,cupp,cdownp
      REAL    :: c(0:lmax*2)


      DO iType = 1, atoms%ntype
         DO nn =1, atoms%neq(iType)
            natom = SUM(atoms%neq(:iType-1)) + 1 !Representative atom
            if (.not.banddos%dos_atom(natom)) cycle
            !find index for dos
            DO n_dos=1,size(banddos%dos_atomlist)
               if (banddos%dos_atomlist(n_dos)==natom) exit
            ENDDO
            DO iBand = 1, noccbd
               j_ind = 0
               c = 0.0
               DO l = 0, lmax
                  IF(l == 0) THEN
                     !s-states (are not split up by SOC)
                     DO spin = 1, input%jspins
                        c(0) = c(0) + eigVecCoeffs%abcof(iBand,0,0,natom,spin)*CONJG(eigVecCoeffs%abcof(iBand,0,0,natom,spin)) &
                                    + eigVecCoeffs%abcof(iBand,0,1,natom,spin)*CONJG(eigVecCoeffs%abcof(iBand,0,1,natom,spin)) &
                                       *usdus%ddn(0,iType,spin)

                        DO ilo  = 1, atoms%nlo(iType)
                           IF(atoms%llo(ilo,iType) /= 0) CYCLE
                           c(0) = c(0) + 2*REAL(eigVecCoeffs%abcof(iBand,0,0,natom,spin)*eigVecCoeffs%ccof(0,iBand,ilo,natom,spin))* usdus%uulon(ilo,iType,spin) &
                                       + 2*REAL(eigVecCoeffs%abcof(iBand,0,1,natom,spin)*eigVecCoeffs%ccof(0,iBand,ilo,natom,spin))* usdus%dulon(ilo,iType,spin)
                           DO ilop  = 1, atoms%nlo(iType)
                              IF(atoms%llo(ilo,iType) /= 0) CYCLE
                              c(0) = c(0) + eigVecCoeffs%ccof(0,iBand,ilo,natom,spin)*CONJG(eigVecCoeffs%ccof(0,iBand,ilop,natom,spin))*usdus%uloulopn(ilo,ilop,iType,spin)
                           ENDDO
                        ENDDO
                     ENDDO
                  ELSE
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
                              aup   = facup   * eigVecCoeffs%abcof(iBand,lmup,0,natom,1)
                              bup   = facup   * eigVecCoeffs%abcof(iBand,lmup,1,natom,1)
                           ELSE
                              aup = 0.0
                              bup = 0.0
                           ENDIF

                           IF(ABS(mdown) <= l) THEN
                              lmdown = l*(l+1) + INT(mdown)
                              facdown = clebsch(REAL(l),0.5,mdown,-0.5,j,mj)
                              adown = - facdown * eigVecCoeffs%abcof(iBand,lmdown,0,natom,spin)
                              bdown = - facdown * eigVecCoeffs%abcof(iBand,lmdown,1,natom,spin)
                           ELSE
                              adown = 0.0
                              bdown = 0.0
                           ENDIF

                           !c := norm of facup |up> + facdown |down>
                           !We have to write it out explicitely because
                           !of the offdiagonal scalar products that appear
                           c(j_ind) = c(j_ind) &
                                     +        aup  *CONJG(aup)   &
                                     +        adown*CONJG(adown) &
                                     +        bup  *CONJG(bup)    * usdus%ddn(l,iType,1) &
                                     +        bdown*CONJG(bdown)  * usdus%ddn(l,iType,spin) &
                                     + 2*REAL(aup  *CONJG(adown)) * denCoeffsOffdiag%uu21n(l,iType) &
                                     + 2*REAL(bup  *CONJG(bdown)) * denCoeffsOffdiag%dd21n(l,iType) &
                                     + 2*REAL(aup  *CONJG(bdown)) * denCoeffsOffdiag%ud21n(l,iType) &
                                     + 2*REAL(adown*CONJG(bup))   * denCoeffsOffdiag%du21n(l,iType)

                           !Local orbitals
                           DO ilo = 1, atoms%nlo(iType)
                              IF(atoms%llo(ilo,iType) /= l) CYCLE

                              IF(ABS(mup) <= l) THEN
                                 cup = facup  * eigVecCoeffs%ccof(INT(mup),iBand,ilo,natom,1)
                              ELSE
                                 cup = 0.0
                              ENDIF

                              IF(ABS(mdown) <= l) THEN
                                 cdown = - facdown  * eigVecCoeffs%ccof(INT(mdown),iBand,ilo,natom,spin)
                              ELSE
                                 cdown = 0.0
                              ENDIF

                              !Local orbital times ab coeff contribution
                              c(j_ind) = c(j_ind) &
                                        + 2*REAL(aup  *CONJG(cup))   * usdus%uulon(ilo,iType,1) &
                                        + 2*REAL(adown*CONJG(cdown)) * usdus%uulon(ilo,iType,spin) &
                                        + 2*REAL(bup  *CONJG(cup))   * usdus%dulon(ilo,iType,1) &
                                        + 2*REAL(bdown*CONJG(cdown)) * usdus%dulon(ilo,iType,spin) &
                                        + 2*REAL(cup  *CONJG(adown)) * denCoeffsOffdiag%uulo21n(ilo,iType) &
                                        + 2*REAL(cdown*CONJG(aup))   * denCoeffsOffdiag%ulou21n(ilo,iType) &
                                        + 2*REAL(cup  *CONJG(bdown)) * denCoeffsOffdiag%dulo21n(ilo,iType) &
                                        + 2*REAL(cdown*CONJG(bup))   * denCoeffsOffdiag%ulod21n(ilo,iType)

                              !Local orbital times Local orbital contribution
                              DO ilop = 1, atoms%nlo(iType)
                                 IF(atoms%llo(ilop,iType) /= l) CYCLE

                                 IF(ABS(mup) <= l) THEN
                                    cupp = facup  * eigVecCoeffs%ccof(INT(mup),iBand,ilop,natom,1)
                                 ELSE
                                    cupp = 0.0
                                 ENDIF

                                 IF(ABS(mdown) <= l) THEN
                                    cdownp = - facdown  * eigVecCoeffs%ccof(INT(mdown),iBand,ilop,natom,spin)
                                 ELSE
                                    cdownp = 0.0
                                 ENDIF

                                 c(j_ind) = c(j_ind) &
                                           +        cup  *CONJG(cupp)    * usdus%uloulopn(ilo,ilop,iType,1) &
                                           +        cdown*CONJG(cdownp)  * usdus%uloulopn(ilo,ilop,iType,spin) &
                                           + 2*REAL(cup  *CONJG(cdownp)) * denCoeffsOffDiag%uloulop21n(ilo,ilop,iType)
                              ENDDO
                           ENDDO

                           mj = mj + 1
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               summed = SUM(c(0:2*lmax))
               cf = 100.0/summed
               j_ind=0
               DO l = 0, 3
                  DO jj = 1, 2
                     IF(l /= 0) j_ind = j_ind+1
                     jDOS%comp(ev_list(iBand),l,jj,n_dos,ikpt) = c(j_ind)*cf
                     jDOS%qmtp(ev_list(iBand),n_dos,ikpt) = 100.0*summed
                     jDOS%occ(l,jj,natom) = jDOS%occ(l,jj,n_dos) + we(iBand) * c(j_ind)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO
      ENDDO

   END SUBROUTINE jDOS_comp
END MODULE m_jDOS
