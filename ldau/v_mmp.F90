!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_vmmp
  !     ************************************************************
  !     This subroutine calculates the potential matrix v^{s}_{m,m'}
  !     for a given atom 'n' and l-quantum number 'l'. The l,u,j's for
  !     all atoms are stored in lda_u, the density matrix is ns_mmp,
  !     and the e-e- interaction potential is u(m1,m2,m3,m4,n).
  !     For details see Eq.(16) of Shick et al. PRB 60, 10765 (1999)
  !
  !     Additionally, the total energy contribution of LDA+U (Eq.24)
  !     is calculated (e_ldau).
  !     Part of the LDA+U package                   G.B., Oct. 2000
  !     ************************************************************
CONTAINS
  SUBROUTINE v_mmp(atoms,jspins,lmaxb,ns_mmp,u,f0,f2, vs_mmp,results)

    USE m_types
    IMPLICIT NONE

    TYPE(t_results),INTENT(INOUT)   :: results
    TYPE(t_atoms),INTENT(IN)        :: atoms
    !
    ! ..  Arguments ..
    INTEGER, INTENT (IN) :: lmaxb,jspins 
    REAL,    INTENT (IN) :: u(-lmaxb:lmaxb,-lmaxb:lmaxb, -lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u)
    REAL,    INTENT (IN) :: f0(atoms%n_u),f2(atoms%n_u)

    COMPLEX           :: ns_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,jspins)
    COMPLEX,INTENT(OUT)::vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,jspins)
    !
    ! ..  Local Variables ..
    INTEGER n,ispin,jspin,l ,mp,p,q,itype,m
    REAL rho_tot,u_htr,j_htr,e_ee,ns_sum,spin_deg,e_dc,e_dcc
    REAL rho_sig(jspins),v_diag(jspins),eta(0:jspins)
    !
    ! Use around-mean-field limit (true) of atomic limit (false)
    !
    !
    ! Loop over atoms
    !
    spin_deg = 1.0 / (3 - jspins)
    results%e_ldau = 0.0
    n = 0
    DO itype = 1, atoms%ntype
       l = atoms%lda_u(itype)%l
       IF (l.GE.0) THEN
          n = n + 1
          u_htr = atoms%lda_u(itype)%u / 27.21
          j_htr = atoms%lda_u(itype)%j / 27.21
          u_htr = f0(n)/27.21
          IF (l.EQ.1) THEN
             j_htr = f2(n)/(5*27.21)
          ELSEIF (l.EQ.2) THEN
             j_htr = 1.625*f2(n)/(14*27.21)
          ELSEIF (l.EQ.3) THEN
             j_htr = (286.+195*451/675+250*1001/2025)*f2(n)/(6435*27.21)
          ENDIF
          !
          ! calculate spin-density 'rho_sig' and total density 'rho_tot'
          !
          rho_tot = 0.0
          DO ispin = 1,jspins
             rho_sig(ispin) = 0.0
             DO m = -l,l
                rho_sig(ispin) = rho_sig(ispin) + REAL(ns_mmp(m,m,n,ispin))
             ENDDO
             rho_tot = rho_tot + rho_sig(ispin)
          ENDDO
          rho_sig(1) = rho_sig(1) * spin_deg  ! if jspins = 1, divide by 2
          IF (atoms%lda_u(itype)%l_amf) THEN
             eta(1) = rho_sig(1) / (2*l + 1) 
             eta(jspins) = rho_sig(jspins) / (2*l + 1) 
             eta(0) = (eta(1) + eta(jspins) ) / 2
          ELSE
             eta(0) = 1.0
             eta(1) = 1.0
             eta(jspins) = 1.0
          ENDIF
          !
          !--------------------------------------------------------------------------------------------+
          !  s       --                                        s'                    1        s   1    |
          ! V     =  >  ( <m,p|V|m',q> - <m,p|V|q,m'> d     ) n     + d    ( -U (n - -) + J (n  - -) ) |
          !  m,m'    --                                s,s'    p,q     m,m'          2            2    |
          !        p,q,s'                                                                              |
          !--------------------------------------------------------------------------------------------+     
          ! initialise vs_mmp
          !
#ifdef CPP_INVERSION
          vs_mmp(:,:,n,:) = ns_mmp(:,:,n,:)
          DO ispin = 1,jspins
             DO m = -l,l
                DO mp = -l,l
                   ns_mmp(m,mp,n,ispin) = vs_mmp(-m,-mp,n,ispin)
                ENDDO
             ENDDO
          ENDDO
#endif
          vs_mmp(:,:,n,:) = CMPLX(0.0,0.0)
          !
          ! outer spin loop - set up v_mmp
          !
          DO ispin = 1,jspins
             DO m = -l,l
                DO mp =-l,l 

                   DO jspin = 1,jspins
                      IF (ispin.EQ.jspin) THEN
                         DO p = -l,l
                            DO q = -l,l
                               vs_mmp(m,mp,n,ispin) = vs_mmp(m,mp,n,ispin) +  &
                                    ns_mmp(p, q,n,jspin) * ( u(m,p,mp,q,n) - u(m,p,q,mp,n) ) 
                            ENDDO
                         ENDDO
                      ENDIF
                      IF ((ispin.NE.jspin).OR.(jspins.EQ.1)) THEN
                         DO p = -l,l
                            DO q = -l,l
                               vs_mmp(m,mp,n,ispin) = vs_mmp(m,mp,n,ispin) +  &
                                    u(m,p,mp,q,n) * ns_mmp(p, q,n,jspin) 
                            ENDDO
                         ENDDO
                      ENDIF
                   ENDDO

                ENDDO ! m' loop
             ENDDO   ! m  loop
          ENDDO      ! outer spin loop
          !
          !  set diagonal terms and correct for non-spin-polarised case
          !
          DO ispin = 1,jspins
             v_diag(ispin) = - u_htr * ( rho_tot - 0.5*eta(0) ) + j_htr * ( rho_sig(ispin) - 0.5*eta(ispin) )
             DO m = -l,l
                DO mp = -l,l
                   vs_mmp(m,mp,n,ispin) = vs_mmp(m,mp,n,ispin) * spin_deg
                ENDDO
                vs_mmp(m,m,n,ispin) = vs_mmp(m,m,n,ispin) + v_diag(ispin)
             ENDDO
          ENDDO

          !----------------------------------------------------------------------+
          !              s                                                       !
          !  ee      1  ---   s        s                     1        s  1       !
          ! E  (n) = -  >    n      ( V     + d     ( U (n - -) - J (n - -) ))   !
          !          2  ---   m,m'     m,m'    m,m'          2           2       !
          !             m,m'                                                     !
          !----------------------------------------------------------------------+

          e_ee = 0.0
          DO ispin = 1,jspins
             DO m = -l,l
                DO mp =-l,l
                   e_ee=e_ee+REAL(vs_mmp(m,mp,n,ispin)*ns_mmp(m,mp,n,ispin))
                ENDDO
                e_ee = e_ee - v_diag(ispin) * REAL( ns_mmp(m,m,n,ispin) )
             ENDDO
          ENDDO

          !----------------------------------------------------------------------+
          !   dc       ee      U           J  --   s   s       1                 |
          !  E      = E  (n) - - n (n-1) + -  >   n  (n -1)  - - (U-J) n         |
          !   LDA+U            2           2  --               2                 |
          !                                    s                                 |
          !----------------------------------------------------------------------+

          ns_sum = 0.0
          DO ispin = 1,jspins
             ns_sum = ns_sum + rho_sig(ispin) * &
                  ( rho_sig(ispin) - eta(ispin) )
          ENDDO
          e_dc = u_htr * rho_tot * ( rho_tot - eta(0) ) - j_htr * ns_sum 
          e_dcc = (u_htr - j_htr) * rho_tot

          ns_sum = ns_sum / spin_deg
          !       e_ldau = e_ldau + (e_ee -  u_htr * rho_tot * ( rho_tot - 1. ) 
          !    +    + j_htr * ns_sum  - (u_htr - j_htr) * rho_tot) * neq(itype)
          !       write(*,*) e_ldau
          results%e_ldau = results%e_ldau + ( e_ee - e_dc - e_dcc) * atoms%neq(itype)
          !       write(*,*) e_ldau

       ENDIF
    ENDDO ! loop over atoms

    results%e_ldau = results%e_ldau / 2

  END SUBROUTINE v_mmp
END MODULE m_vmmp
