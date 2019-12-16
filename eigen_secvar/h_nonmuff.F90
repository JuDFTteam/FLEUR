!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hnonmuff
  !*********************************************************************
  !     updates hamiltonian by adding non-spherical matrix elements in
  !     the second-variation scheme. usage of tlmplm-nonmuff required
  !                r. p  1995
  !*********************************************************************
CONTAINS
  SUBROUTINE h_nonmuff(atoms,input,sym,cell, jsp,ne,usdus,td, bkpt,lapw, h,l_real,z_r,z_c)

    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_sphbes
    USE m_dsphbs
    USE m_ylm
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_usdus),INTENT(IN)       :: usdus
    TYPE(t_lapw),INTENT(IN)        :: lapw
    !     ..
    !     .. Scalar Arguments ..
    LOGICAL,INTENT(IN)   :: l_real
    INTEGER, INTENT (IN) :: jsp,ne
    !     ..
    TYPE(t_tlmplm),INTENT(IN)::td
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: bkpt(3)
    REAL,    INTENT (INOUT) :: h(ne*(ne+1)/2)

    REAL,    OPTIONAL,INTENT (IN) :: z_r(lapw%dim_nbasfcn(),ne)
    COMPLEX, OPTIONAL,INTENT (IN) :: z_c(lapw%dim_nbasfcn(),ne)
    !     ..
    !     .. Local Scalars ..
    COMPLEX dtd,dtu,hij,phase,sij,utd,utu
    REAL con1,ff,gg,gs,th,ws
    INTEGER l,l1,ll1,lm,lmp,lwn,invsfct
    INTEGER i,im,in,j,k,ke ,m1,n,na,nn,np,ii,ij,m
    !     ..
    !     .. Local Arrays ..
    COMPLEX a(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2)),ax(input%neig)
    COMPLEX b(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2)),bx(input%neig), ylm( (atoms%lmaxd+1)**2 )
    REAL vmult(3),vsmult(3),f(0:atoms%lmaxd,SIZE(lapw%k1,1)),g(0:atoms%lmaxd,SIZE(lapw%k1,1))
    !     ..
    !     ..

    con1 = fpi_const/SQRT(cell%omtil)
    !--->    loop over each atom type
    na = 0
    DO  n = 1,atoms%ntype
       lwn = atoms%lmax(n)
       !--->    set up wronskians for the matching conditions for each ntype
       DO k = 1,lapw%nv(jsp)
          gs = lapw%rk(k,jsp)*atoms%rmt(n)
          CALL sphbes(lwn,gs, f(0,k))
          CALL dsphbs(lwn,gs,f(0,k), g(0,k))
       ENDDO
       DO  l = 0,lwn
          ws = usdus%uds(l,n,jsp)*usdus%dus(l,n,jsp) - usdus%us(l,n,jsp)*usdus%duds(l,n,jsp)
          DO  k = 1,lapw%nv(jsp)
             ff = f(l,k)
             gg = lapw%rk(k,jsp)*g(l,k)
             f(l,k) = con1* (usdus%uds(l,n,jsp)*gg-ff*usdus%duds(l,n,jsp))/ws
             g(l,k) = con1* (usdus%dus(l,n,jsp)*ff-gg*usdus%us(l,n,jsp))/ws
          ENDDO
       ENDDO
       !--->    loop over equivalent atoms
       DO  nn = 1,atoms%neq(n)
          na = na + 1
          !+inv
          IF ((sym%invsat(na).EQ.0) .OR. (sym%invsat(na).EQ.1)) THEN
             IF (sym%invsat(na).EQ.0) invsfct = 1
             IF (sym%invsat(na).EQ.1) invsfct = 2
             np = sym%ngopr(na)
             !---> a and b
             a(:ne,:) = CMPLX(0.0,0.0)
             b(:ne,:) = CMPLX(0.0,0.0)
             DO  k = 1,lapw%nv(jsp)
                vmult=bkpt+(/lapw%k1(k,jsp),lapw%k2(k,jsp),lapw%k3(k,jsp)/)
                th = tpi_const* DOT_PRODUCT( vmult,atoms%taual(:,na))
                phase = CMPLX(COS(th),-SIN(th))
                !-->     apply the rotation that brings this atom into the
                !-->     representative for hamiltonian (this is the definition
                !-->     of ngopr(na)) and transform to cartesian coordinates
                vsmult=MATMUL(vmult,sym%mrot(:,:,np))
                vmult=MATMUL(vsmult,cell%bmat)
                CALL ylm4(lwn,vmult, ylm)
                !-->     synthesize the complex conjugates of a and b
                if (l_real) THEN
                   DO l = 0,lwn
                      ll1 = l* (l+1)
                      DO m = -l,l
                         lm = ll1 + m
                         hij = f(l,k) * ( phase * ylm(lm+1) )
                         sij = g(l,k) * ( phase * ylm(lm+1) )
                         a(:ne,lm) = a(:ne,lm) + hij*z_r(k,:ne)
                         b(:ne,lm) = b(:ne,lm) + sij*z_r(k,:ne)
                      END DO
                   END DO

                else
                   DO l = 0,lwn
                      ll1 = l* (l+1)
                      DO m = -l,l
                         lm = ll1 + m
                         hij = f(l,k) * ( phase * ylm(lm+1) )
                         sij = g(l,k) * ( phase * ylm(lm+1) )
                         a(:ne,lm) = a(:ne,lm) + hij*z_c(k,:ne)
                         b(:ne,lm) = b(:ne,lm) + sij*z_c(k,:ne)
                      END DO
                   END DO
                endif
             ENDDO
             DO  l = 0,lwn
                DO  m = -l,l
                   lmp = l* (l+1) + m
                   !--->    initialize ax and bx
                   ax = CMPLX(0.0,0.0)
                   bx = CMPLX(0.0,0.0)
                   !--->    loop over l,m
                   DO  l1 = 0,lwn
                      DO  m1 = -l1,l1
                         lm = l1* (l1+1) + m1
                         in = td%ind(lmp,lm,nn,jsp)
                         IF (in.NE.-9999) THEN
                            IF (in.GE.0) THEN
                               utu = CONJG(td%tuu(in,nn,jsp))*invsfct
                               dtu = CONJG(td%tdu(in,nn,jsp))*invsfct
                               utd = CONJG(td%tud(in,nn,jsp))*invsfct
                               dtd = CONJG(td%tdd(in,nn,jsp))*invsfct
                            ELSE
                               im = -in
                               utu = td%tuu(im,nn,jsp)*invsfct
                               dtd = td%tdd(im,nn,jsp)*invsfct
                               utd = td%tdu(im,nn,jsp)*invsfct
                               dtu = td%tud(im,nn,jsp)*invsfct
                            END IF
                            !--->    update ax, bx
                            ax(:ne) = ax(:ne) + CONJG(utu*a(:ne,lm)+utd*b(:ne,lm))
                            bx(:ne) = bx(:ne) + CONJG(dtu*a(:ne,lm)+dtd*b(:ne,lm))
                         END IF
                      ENDDO
                   ENDDO

                   !--->    update hamiltonian in upper packed storage mode
                   DO i = 1,ne
                      ii = (i-1)*i/2
                      DO j = 1,i - 1
                         ij = ii + j
                         hij = a(i,lmp)*ax(j) + b(i,lmp)*bx(j)
                         h(ij) = h(ij) + REAL(hij)
                      END DO
                      h(ii+i) = h(ii+i) + REAL(a(i,lmp)*ax(i)+&
                           &                        b(i,lmp)*bx(i))
                   END DO
                ENDDO
             ENDDO

          ENDIF
          !-inv
       ENDDO
    ENDDO
  END SUBROUTINE h_nonmuff
END MODULE m_hnonmuff
