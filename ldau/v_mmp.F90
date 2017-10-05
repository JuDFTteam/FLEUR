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
  !     
  !     Extension to multiple U per atom type  G.M. 2017
  !     ************************************************************
CONTAINS
  SUBROUTINE v_mmp(sym,atoms,jspins,ns_mmp,u,f0,f2, vs_mmp,results)

    USE m_types
    USE m_constants
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)          :: sym
    TYPE(t_results),INTENT(INOUT)   :: results
    TYPE(t_atoms),INTENT(IN)        :: atoms
    !
    ! ..  Arguments ..
    INTEGER, INTENT(IN)    :: jspins 
    REAL,    INTENT(IN)    :: u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
                                -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_u)
    REAL,    INTENT(IN)    :: f0(atoms%n_u),f2(atoms%n_u)
    COMPLEX, INTENT(OUT)   :: vs_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_u,jspins)

    COMPLEX, INTENT(INOUT) :: ns_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_u,jspins)

    ! ..  Local Variables ..
    INTEGER ispin,jspin,l ,mp,p,q,itype,m,i_u
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

    DO i_u = 1, atoms%n_u
       iType = atoms%lda_u(i_u)%atomType
       l = atoms%lda_u(i_u)%l
       u_htr = atoms%lda_u(i_u)%u / hartree_to_ev_const
       j_htr = atoms%lda_u(i_u)%j / hartree_to_ev_const
       u_htr = f0(i_u)/hartree_to_ev_const
       IF (l.EQ.1) THEN
          j_htr = f2(i_u)/(5*hartree_to_ev_const)
       ELSE IF (l.EQ.2) THEN
          j_htr = 1.625*f2(i_u)/(14*hartree_to_ev_const)
       ELSE IF (l.EQ.3) THEN
          j_htr = (286.+195*451/675+250*1001/2025)*f2(i_u)/(6435*hartree_to_ev_const)
       END IF
       !
       ! calculate spin-density 'rho_sig' and total density 'rho_tot'
       !
       rho_tot = 0.0
       DO ispin = 1,jspins
          rho_sig(ispin) = 0.0
          DO m = -l,l
             rho_sig(ispin) = rho_sig(ispin) + REAL(ns_mmp(m,m,i_u,ispin))
          END DO
          rho_tot = rho_tot + rho_sig(ispin)
       END DO
       rho_sig(1) = rho_sig(1) * spin_deg  ! if jspins = 1, divide by 2
       IF (atoms%lda_u(i_u)%l_amf) THEN
          eta(1) = rho_sig(1) / (2*l + 1) 
          eta(jspins) = rho_sig(jspins) / (2*l + 1) 
          eta(0) = (eta(1) + eta(jspins) ) / 2
       ELSE
          eta(0) = 1.0
          eta(1) = 1.0
          eta(jspins) = 1.0
       END IF
       !
       !--------------------------------------------------------------------------------------------+
       !  s       --                                        s'                    1        s   1    |
       ! V     =  >  ( <m,p|V|m',q> - <m,p|V|q,m'> d     ) n     + d    ( -U (n - -) + J (n  - -) ) |
       !  m,m'    --                                s,s'    p,q     m,m'          2            2    |
       !        p,q,s'                                                                              |
       !--------------------------------------------------------------------------------------------+     
       ! initialise vs_mmp
       !
       IF (sym%invs) THEN
          vs_mmp(:,:,i_u,:) = ns_mmp(:,:,i_u,:)
          DO ispin = 1,jspins
             DO m = -l,l
                DO mp = -l,l
                   ns_mmp(m,mp,i_u,ispin) = vs_mmp(-m,-mp,i_u,ispin)
                END DO
             END DO
          END DO
       END IF
       vs_mmp(:,:,i_u,:) = CMPLX(0.0,0.0)
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
                            vs_mmp(m,mp,i_u,ispin) = vs_mmp(m,mp,i_u,ispin) +  &
                               ns_mmp(p, q,i_u,jspin) * ( u(m,p,mp,q,i_u) - u(m,p,q,mp,i_u) ) 
                         END DO
                      END DO
                   END IF
                   IF ((ispin.NE.jspin).OR.(jspins.EQ.1)) THEN
                      DO p = -l,l
                         DO q = -l,l
                            vs_mmp(m,mp,i_u,ispin) = vs_mmp(m,mp,i_u,ispin) +  &
                               u(m,p,mp,q,i_u) * ns_mmp(p, q,i_u,jspin) 
                         END DO
                      END DO
                   END IF
                END DO

             END DO ! m' loop
          END DO   ! m  loop
       END DO      ! outer spin loop
       !
       !  set diagonal terms and correct for non-spin-polarised case
       !
       DO ispin = 1,jspins
          v_diag(ispin) = - u_htr * ( rho_tot - 0.5*eta(0) ) + j_htr * ( rho_sig(ispin) - 0.5*eta(ispin) )
          DO m = -l,l
             DO mp = -l,l
                vs_mmp(m,mp,i_u,ispin) = vs_mmp(m,mp,i_u,ispin) * spin_deg
             END DO
             vs_mmp(m,m,i_u,ispin) = vs_mmp(m,m,i_u,ispin) + v_diag(ispin)
          END DO
       END DO

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
                e_ee=e_ee+REAL(vs_mmp(m,mp,i_u,ispin)*ns_mmp(m,mp,i_u,ispin))
             END DO
             e_ee = e_ee - v_diag(ispin) * REAL( ns_mmp(m,m,i_u,ispin) )
          END DO
       END DO

       !----------------------------------------------------------------------+
       !   dc       ee      U           J  --   s   s       1                 |
       !  E      = E  (n) - - n (n-1) + -  >   n  (n -1)  - - (U-J) n         |
       !   LDA+U            2           2  --               2                 |
       !                                    s                                 |
       !----------------------------------------------------------------------+

       ns_sum = 0.0
       DO ispin = 1,jspins
          ns_sum = ns_sum + rho_sig(ispin) * (rho_sig(ispin) - eta(ispin))
       END DO
       e_dc = u_htr * rho_tot * ( rho_tot - eta(0) ) - j_htr * ns_sum
       e_dcc = (u_htr - j_htr) * rho_tot

       ns_sum = ns_sum / spin_deg
       !       e_ldau = e_ldau + (e_ee -  u_htr * rho_tot * ( rho_tot - 1. ) 
       !    +    + j_htr * ns_sum  - (u_htr - j_htr) * rho_tot) * neq(itype)
       !       write(*,*) e_ldau
       results%e_ldau = results%e_ldau + ( e_ee - e_dc - e_dcc) * atoms%neq(itype)
       !       write(*,*) e_ldau

    END DO ! loop over U parameters

    results%e_ldau = results%e_ldau / 2

  END SUBROUTINE v_mmp
END MODULE m_vmmp
