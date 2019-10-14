MODULE m_lhglptg
  !.....------------------------------------------------------------------
  !     calculates lattice harmonics and their gradients on the
  !       gauss-legendre angular mesh - r.p. and t.a.
  !     for gradient. t.a. 1996.
  !.....------------------------------------------------------------------
CONTAINS
  SUBROUTINE lhglptg(&
       &                   sphhar,atoms,&
       &                   rx,nsp,dograds,sym,&
       &                   ylh,thet,phi,ylht1,ylht2,ylhf1,ylhf2,ylhtf)
    !
    USE m_polangle
    USE m_ylm
    USE m_dylm
    USE m_types
    IMPLICIT NONE

    LOGICAL, INTENT(IN)         :: dograds
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nsp  
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: rx(:,:)!(3,dimension%nspd)
    REAL,    INTENT (OUT):: thet(:) !nspd
    REAL,    INTENT (OUT):: phi(:) !nspd
    REAL,    INTENT (OUT):: ylh(:,0:,:)!(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd),thet(nspd)
    REAL,    INTENT (OUT):: ylht1(:,0:,:)!(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd)
    REAL,    INTENT (OUT):: ylht2(:,0:,:)!(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd)
    REAL,    INTENT (OUT):: ylhtf(:,0:,:)!(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd)
    REAL,    INTENT (OUT):: ylhf1(:,0:,:)!(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd)
    REAL,    INTENT (OUT):: ylhf2(:,0:,:)!(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd)
    !     ..
    !     .. Local Scalars ..
    REAL s,st1,st2,sf1,sf2,stf
    INTEGER k,lh,mem,nd,lm,ll1
    !     ..
    !     .. Local Arrays ..
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    COMPLEX dylmt1( (atoms%lmaxd+1)**2 ), dylmt2( (atoms%lmaxd+1)**2 )
    COMPLEX dylmf1( (atoms%lmaxd+1)**2 ), dylmf2( (atoms%lmaxd+1)**2 )
    COMPLEX dylmtf( (atoms%lmaxd+1)**2 )
    !     ..

    !.....------------------------------------------------------------------
    !     ..
    DO  nd = 1,sym%nsymt

       DO  k = 1,nsp

          CALL ylm4(&
               &                   atoms%lmaxd,rx(:,k),&
               &                   ylm)
          CALL pol_angle(&
               &                       rx(1,k),rx(2,k),rx(3,k),&
               &                       thet(k),phi(k))

          IF (dograds) THEN
             CALL dylm3(&
                  &                     atoms%lmaxd,atoms%lmaxd,rx(:,k),ylm,&
                  &                     dylmt1,dylmt2,dylmf1,dylmf2,dylmtf)
          ENDIF

          DO  lh = 0,sphhar%nlh(nd)
             s   = 0
             st1 = 0
             st2 = 0
             sf1 = 0
             sf2 = 0
             stf = 0
             ll1 = sphhar%llh(lh,nd) * ( sphhar%llh(lh,nd) + 1 ) + 1

             DO mem = 1,sphhar%nmem(lh,nd)
                lm = ll1 + sphhar%mlh(mem,lh,nd)
                s = s + REAL( sphhar%clnu(mem,lh,nd) * ylm(lm) )
             ENDDO

             ylh(k,lh,nd) = s

             IF (dograds) THEN

                DO mem = 1,sphhar%nmem(lh,nd)
                   lm = ll1 + sphhar%mlh(mem,lh,nd)
                   s   = s   + REAL( sphhar%clnu(mem,lh,nd)* ylm(lm) )
                   st1 = st1 + REAL( sphhar%clnu(mem,lh,nd)*dylmt1(lm) )
                   st2 = st2 + REAL( sphhar%clnu(mem,lh,nd)*dylmt2(lm) )
                   sf1 = sf1 + REAL( sphhar%clnu(mem,lh,nd)*dylmf1(lm) )
                   sf2 = sf2 + REAL( sphhar%clnu(mem,lh,nd)*dylmf2(lm) )
                   stf = stf + REAL( sphhar%clnu(mem,lh,nd)*dylmtf(lm) )
                ENDDO

                ylht1(k,lh,nd) = st1
                ylht2(k,lh,nd) = st2
                ylhf1(k,lh,nd) = sf1
                ylhf2(k,lh,nd) = sf2
                ylhtf(k,lh,nd) = stf

             ENDIF

          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE lhglptg
END MODULE m_lhglptg
