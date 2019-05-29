!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_wann_orbcomp
  !*************************************************************************
  !     Compute matrix elements of the first 16 angular momentum projectors
  !     within the MT-spheres. This is done separately for each atom type.
  !
  !     These are:
  !     1: |s><s|
  !     2: |p_x><p_x|
  !     3: |p_y><p_y|
  !     4: |p_z><p_z|
  !     5: |d_xy><d_xy|
  !     6: |d_yz><d_yz|
  !     7: |d_zx><d_zx|
  !     8: |d_x2-y2><d_x2-y2|
  !     9: |d_z2><dz2|
  !     10: |f><f|
  !     11: |f><f|
  !     12: |f><f|
  !     13: |f><f|
  !     14: |f><f|
  !     15: |f><f|
  !     16: |f><f|
  !
  !     Based on code in orb_comp by Yury Koroteev
  !
  !     Frank Freimuth
  !***********************************************************************
CONTAINS
  SUBROUTINE wann_orbcomp(atoms,usdus,jspin,&
       acof,bcof,ccof,&
       num_orbs,&
       orbs,&
       l_oc_f,&
       mmn)
    USE m_types
    IMPLICIT NONE
    !     .. scalar arguments ..
    TYPE(t_atoms),INTENT(in) :: atoms
    TYPE(t_usdus),INTENT(IN) :: usdus
    INTEGER,INTENT(IN)       :: jspin
    !     .. array arguments ..
    COMPLEX, INTENT (in)  :: ccof(-atoms%llod:,:,:,:) !ccof(-llod:llod,noccbd,nlod,natd)
    COMPLEX, INTENT (in)  :: acof(:,0:,:)!acof(noccbd,0:lmd,natd)
    COMPLEX, INTENT (in)  :: bcof(:,0:,:)!bcof(noccbd,0:lmd,natd)
    COMPLEX, INTENT (inout) :: mmn(:,:,:,:)!mmn(16,noccbd,noccbd,ntype)

    INTEGER, INTENT (in)  :: num_orbs
    INTEGER, INTENT (in)  :: orbs(:)
    LOGICAL, INTENT (in)  :: l_oc_f

    !    ..Local Scalars 

    COMPLEX  ca00(SIZE(acof,1))
    COMPLEX  ca01(SIZE(acof,1))
    COMPLEX  ca02(SIZE(acof,1))
    COMPLEX  ca03(SIZE(acof,1))
    COMPLEX  ca04(SIZE(acof,1))
    COMPLEX  ca05(SIZE(acof,1))
    COMPLEX  ca06(SIZE(acof,1))
    COMPLEX  ca07(SIZE(acof,1))
    COMPLEX  ca08(SIZE(acof,1))
    COMPLEX  ca09(SIZE(acof,1))
    COMPLEX  ca10(SIZE(acof,1))
    COMPLEX  ca11(SIZE(acof,1))
    COMPLEX  ca12(SIZE(acof,1))
    COMPLEX  ca13(SIZE(acof,1))
    COMPLEX  ca14(SIZE(acof,1))
    COMPLEX  ca15(SIZE(acof,1))

    COMPLEX  cb00(SIZE(acof,1))
    COMPLEX  cb01(SIZE(acof,1))
    COMPLEX  cb02(SIZE(acof,1))
    COMPLEX  cb03(SIZE(acof,1))
    COMPLEX  cb04(SIZE(acof,1))
    COMPLEX  cb05(SIZE(acof,1))
    COMPLEX  cb06(SIZE(acof,1))
    COMPLEX  cb07(SIZE(acof,1))
    COMPLEX  cb08(SIZE(acof,1))
    COMPLEX  cb09(SIZE(acof,1))
    COMPLEX  cb10(SIZE(acof,1))
    COMPLEX  cb11(SIZE(acof,1))
    COMPLEX  cb12(SIZE(acof,1))
    COMPLEX  cb13(SIZE(acof,1))
    COMPLEX  cb14(SIZE(acof,1))
    COMPLEX  cb15(SIZE(acof,1))

    COMPLEX  cc00(SIZE(acof,1))
    COMPLEX  cc01(SIZE(acof,1))
    COMPLEX  cc02(SIZE(acof,1))
    COMPLEX  cc03(SIZE(acof,1))
    COMPLEX  cc04(SIZE(acof,1))
    COMPLEX  cc05(SIZE(acof,1))
    COMPLEX  cc06(SIZE(acof,1))
    COMPLEX  cc07(SIZE(acof,1))
    COMPLEX  cc08(SIZE(acof,1))
    COMPLEX  cc09(SIZE(acof,1))
    COMPLEX  cc10(SIZE(acof,1))
    COMPLEX  cc11(SIZE(acof,1))
    COMPLEX  cc12(SIZE(acof,1))
    COMPLEX  cc13(SIZE(acof,1))
    COMPLEX  cc14(SIZE(acof,1))
    COMPLEX  cc15(SIZE(acof,1))


    INTEGER ::  lo,orbi,rep
    REAL    :: ddn0,ddn1,ddn2,ddn3
    INTEGER :: mt,ityp,n,m,mt1,l
    REAL    :: h,g,c3,c5
    LOGICAL :: l_doit
    !      ..
    !      ..
    !      Intrinsic Function
    INTRINSIC sqrt,conjg
    !
    DATA h/0.50/ g/0.0625/
    !**************************************************** 
    c5 = SQRT(5.0)
    c3 = SQRT(3.0)
    !
    mt=0
    DO  ityp = 1,atoms%ntype
       DO rep=1,atoms%neq(ityp)
          mt=mt+1 

          l_doit=.FALSE.
          DO l=1,num_orbs
             IF(orbs(l)==mt)THEN
                l_doit=.TRUE.
                orbi=l
             ENDIF
          ENDDO
          IF(.NOT.l_doit)CYCLE

          ddn0 = usdus%ddn(0,ityp,jspin)
          ddn1 = usdus%ddn(1,ityp,jspin)        
          ddn2 = usdus%ddn(2,ityp,jspin)        
          ddn3 = usdus%ddn(3,ityp,jspin) 
          DO n=1,SIZE(acof,1)
             !
             !acof
             !         s-states
             ca00(n)=acof(n,0,mt)
             !        p-states
             ca01(n) = acof(n,1,mt) - acof(n,3,mt)
             ca02(n) = acof(n,1,mt) + acof(n,3,mt)
             ca03(n) = acof(n,2,mt)
             !        d-states
             ca04(n) = acof(n,4,mt) - acof(n,8,mt)
             ca05(n) = acof(n,5,mt) + acof(n,7,mt)
             ca06(n) = acof(n,5,mt) - acof(n,7,mt)
             ca07(n) = acof(n,4,mt) + acof(n,8,mt)
             ca08(n) = acof(n,6,mt)
             !
             !        f-states: a cubic set (cub) 
             ! 
             ca09(n) = ( acof(n,9,mt)  - acof(n,15,mt) )*c5 -&
                  ( acof(n,11,mt) - acof(n,13,mt) )*c3
             ca10(n) = ( acof(n,9,mt)  + acof(n,15,mt) )*c5 +&
                  ( acof(n,11,mt) + acof(n,13,mt) )*c3 
             ca11(n) =   acof(n,12,mt)
             ca12(n) = ( acof(n,9,mt)  + acof(n,15,mt) )*c3 -&
                  ( acof(n,11,mt) + acof(n,13,mt) )*c5 
             ca13(n) =   acof(n,10,mt) + acof(n,14,mt)
             ca14(n) = ( acof(n,9,mt)  - acof(n,15,mt) )*c3 +&
                  ( acof(n,11,mt) - acof(n,13,mt) )*c5 
             ca15(n) =   acof(n,10,mt) - acof(n,14,mt) 

             !
             !bcof
             !       s-states
             cb00(n) =  bcof(n,0,mt)
             !       p-states
             cb01(n) =  bcof(n,1,mt) - bcof(n,3,mt) 
             cb02(n) =  bcof(n,1,mt) + bcof(n,3,mt) 
             cb03(n) =  bcof(n,2,mt)
             !        d-states
             cb04(n) =  bcof(n,4,mt) - bcof(n,8,mt) 
             cb05(n) =  bcof(n,5,mt) + bcof(n,7,mt) 
             cb06(n) =  bcof(n,5,mt) - bcof(n,7,mt) 
             cb07(n) =  bcof(n,4,mt) + bcof(n,8,mt) 
             cb08(n) =  bcof(n,6,mt)
             !
             !       f-states: a cubic set (cub)
             !
             cb09(n) = ( bcof(n,9,mt)  - bcof(n,15,mt) )*c5 -&
                  ( bcof(n,11,mt) - bcof(n,13,mt) )*c3
             cb10(n) = ( bcof(n,9,mt)  + bcof(n,15,mt) )*c5 +&
                  ( bcof(n,11,mt) + bcof(n,13,mt) )*c3 
             cb11(n) =   bcof(n,12,mt)
             cb12(n) = ( bcof(n,9,mt)  + bcof(n,15,mt) )*c3 -&
                  ( bcof(n,11,mt) + bcof(n,13,mt) )*c5 
             cb13(n) =   bcof(n,10,mt) + bcof(n,14,mt)
             cb14(n) = ( bcof(n,9,mt)  - bcof(n,15,mt) )*c3 +&
                  ( bcof(n,11,mt) - bcof(n,13,mt) )*c5
             cb15(n) =   bcof(n,10,mt) - bcof(n,14,mt) 
          ENDDO !n   

          DO n=1,SIZE(acof,1)
             DO m=1,SIZE(acof,1)
                !  s
                mmn(1,orbi,m,n)  =   ca00(n)*CONJG(ca00(m)) &
                     + cb00(n)*CONJG(cb00(m))*ddn0 
                !  p
                mmn(2,orbi,m,n)  = ( ca01(n)*CONJG(ca01(m)) &
                     + cb01(n)*CONJG(cb01(m))*ddn1 )*h
                mmn(3,orbi,m,n)  = ( ca02(n)*CONJG(ca02(m)) &
                     + cb02(n)*CONJG(cb02(m))*ddn1 )*h
                mmn(4,orbi,m,n)  =   ca03(n)*CONJG(ca03(m)) &
                     + cb03(n)*CONJG(cb03(m))*ddn1 
                !  d
                mmn(5,orbi,m,n)  = ( ca04(n)*CONJG(ca04(m)) &
                     + cb04(n)*CONJG(cb04(m))*ddn2 )*h
                mmn(6,orbi,m,n)  = ( ca05(n)*CONJG(ca05(m)) &
                     + cb05(n)*CONJG(cb05(m))*ddn2 )*h
                mmn(7,orbi,m,n)  = ( ca06(n)*CONJG(ca06(m)) &
                     + cb06(n)*CONJG(cb06(m))*ddn2 )*h
                mmn(8,orbi,m,n)  = ( ca07(n)*CONJG(ca07(m)) &
                     + cb07(n)*CONJG(cb07(m))*ddn2 )*h
                mmn(9,orbi,m,n)  =   ca08(n)*CONJG(ca08(m)) &
                     + cb08(n)*CONJG(cb08(m))*ddn2 
                !  f: a cubic set
                IF(l_oc_f)THEN
                   mmn(10,orbi,m,n) = ( ca09(n)*CONJG(ca09(m)) &
                        + cb09(n)*CONJG(cb09(m))*ddn3 )*g       
                   mmn(11,orbi,m,n) = ( ca10(n)*CONJG(ca10(m)) &
                        + cb10(n)*CONJG(cb10(m))*ddn3 )*g
                   mmn(12,orbi,m,n) =   ca11(n)*CONJG(ca11(m)) &
                        + cb11(n)*CONJG(cb11(m))*ddn3 
                   mmn(13,orbi,m,n) = ( ca12(n)*CONJG(ca12(m)) &
                        + cb12(n)*CONJG(cb12(m))*ddn3 )*g
                   mmn(14,orbi,m,n) = ( ca13(n)*CONJG(ca13(m)) &
                        + cb13(n)*CONJG(cb13(m))*ddn3 )*h
                   mmn(15,orbi,m,n) = ( ca14(n)*CONJG(ca14(m)) &
                        + cb14(n)*CONJG(cb14(m))*ddn3 )*g
                   mmn(16,orbi,m,n) = ( ca15(n)*CONJG(ca15(m)) &
                        + cb15(n)*CONJG(cb15(m))*ddn3 )*h
                ENDIF !l_oc_f
             ENDDO !m
          ENDDO !n 

          !--------------------------------------------------------------------
          !ccof   ( contributions from local orbitals )
          !
          lo_loop:DO lo = 1,atoms%nlo(ityp)
             l = atoms%llo(lo,ityp)
             !lo-s
             IF ( l.EQ.0 )  THEN
                DO n=1,SIZE(acof,1)
                   cc00(n) = ccof(0,n,lo,mt)
                ENDDO
                !     
                DO n=1,SIZE(acof,1)
                   DO m=1,SIZE(acof,1)
                      mmn(1,orbi,m,n)=mmn(1,orbi,m,n)+&
                           ( ca00(n)*CONJG(cc00(m)) + &
                           cc00(n)*CONJG(ca00(m)) )*usdus%uulon(lo,ityp,jspin) + &
                           ( cb00(n)*CONJG(cc00(m)) + &
                           cc00(n)*CONJG(cb00(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc00(n)*CONJG(cc00(m))*usdus%uloulopn(lo,lo,ityp,jspin) 
                   ENDDO !m
                ENDDO !n
                CYCLE lo_loop
             ENDIF
             !lo-p
             IF ( l.EQ.1 )  THEN
                DO n=1,SIZE(acof,1)
                   cc01(n) = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
                   cc02(n) = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
                   cc03(n) = ccof( 0,n,lo,mt)
                ENDDO
                DO n=1,SIZE(acof,1)
                   DO m=1,SIZE(acof,1)
                      mmn(2,orbi,m,n) = mmn(2,orbi,m,n)  +&
                           (   ( ca01(n)*CONJG(cc01(m)) + &
                           cc01(n)*CONJG(ca01(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb01(n)*CONJG(cc01(m)) + &
                           cc01(n)*CONJG(cb01(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc01(n)*CONJG(cc01(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*h       

                      mmn(3,orbi,m,n) = mmn(3,orbi,m,n)  +&
                           (   ( ca02(n)*CONJG(cc02(m)) + &
                           cc02(n)*CONJG(ca02(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb02(n)*CONJG(cc02(m)) + &
                           cc02(n)*CONJG(cb02(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc02(n)*CONJG(cc02(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*h       

                      mmn(4,orbi,m,n) = mmn(4,orbi,m,n)  +&
                           ( ca03(n)*CONJG(cc03(m)) + &
                           cc03(n)*CONJG(ca03(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb03(n)*CONJG(cc03(m)) + &
                           cc03(n)*CONJG(cb03(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc03(n)*CONJG(cc03(m))*usdus%uloulopn(lo,lo,ityp,jspin) 


                   ENDDO !m
                ENDDO !n 
                CYCLE lo_loop
             ENDIF
             !lo-d
             IF ( l.EQ.2 )  THEN
                DO n=1,SIZE(acof,1)
                   cc04(n) = ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
                   cc05(n) = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
                   cc06(n) = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
                   cc07(n) = ccof(-2,n,lo,mt) + ccof(2,n,lo,mt)
                   cc08(n) = ccof( 0,n,lo,mt)
                ENDDO !n
                !
                DO n=1,SIZE(acof,1)
                   DO m=1,SIZE(acof,1)
                      mmn(5,orbi,m,n) = mmn(5,orbi,m,n)  +&
                           (( ca04(n)*CONJG(cc04(m)) + &
                           cc04(n)*CONJG(ca04(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb04(n)*CONJG(cc04(m)) + &
                           cc04(n)*CONJG(cb04(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc04(n)*CONJG(cc04(m))*usdus%uloulopn(lo,lo,ityp,jspin))*h 
                      mmn(6,orbi,m,n) = mmn(6,orbi,m,n)  +&
                           (( ca05(n)*CONJG(cc05(m)) + &
                           cc05(n)*CONJG(ca05(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb05(n)*CONJG(cc05(m)) + &
                           cc05(n)*CONJG(cb05(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc05(n)*CONJG(cc05(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*h 
                      mmn(7,orbi,m,n) = mmn(7,orbi,m,n)  +&
                           (( ca06(n)*CONJG(cc06(m)) + &
                           cc06(n)*CONJG(ca06(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb06(n)*CONJG(cc06(m)) + &
                           cc06(n)*CONJG(cb06(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc06(n)*CONJG(cc06(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*h
                      mmn(8,orbi,m,n) = mmn(8,orbi,m,n)  +&
                           (( ca07(n)*CONJG(cc07(m)) + &
                           cc07(n)*CONJG(ca07(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb07(n)*CONJG(cc07(m)) + &
                           cc07(n)*CONJG(cb07(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc07(n)*CONJG(cc07(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*h       
                      mmn(9,orbi,m,n) = mmn(9,orbi,m,n)  +&
                           ( ca08(n)*CONJG(cc08(m)) + &
                           cc08(n)*CONJG(ca08(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb08(n)*CONJG(cc08(m)) + &
                           cc08(n)*CONJG(cb08(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc08(n)*CONJG(cc08(m))*usdus%uloulopn(lo,lo,ityp,jspin)  
                   ENDDO !m
                ENDDO !n
                CYCLE lo_loop                        
             ENDIF
             !lo-f
             IF ( l.EQ.3 .AND. l_oc_f)  THEN
                !
                !  a cubic set (cub)
                !
                DO n=1,SIZE(acof,1)
                   cc09(n) = ( ccof(-3,n,lo,mt) - ccof(3,n,lo,mt) )*c5 -&
                        ( ccof(-1,n,lo,mt) - ccof(1,n,lo,mt) )*c3 
                   cc10(n) = ( ccof(-3,n,lo,mt) + ccof(3,n,lo,mt) )*c5 +&
                        ( ccof(-1,n,lo,mt) + ccof(1,n,lo,mt) )*c3 
                   cc11(n) =   ccof( 0,n,lo,mt)
                   cc12(n) = ( ccof(-3,n,lo,mt) + ccof(3,n,lo,mt) )*c3 -&
                        ( ccof(-1,n,lo,mt) + ccof(1,n,lo,mt) )*c5 
                   cc13(n) =   ccof(-2,n,lo,mt) + ccof(2,n,lo,mt) 
                   cc14(n) = ( ccof(-3,n,lo,mt) - ccof(3,n,lo,mt) )*c3 +&
                        ( ccof(-1,n,lo,mt) - ccof(1,n,lo,mt) )*c5
                   cc15(n) =   ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
                ENDDO !n    
                DO n=1,SIZE(acof,1)
                   DO m=1,SIZE(acof,1)
                      mmn(10,orbi,m,n) = mmn(10,orbi,m,n)  +&
                           (( ca09(n)*CONJG(cc09(m)) + &
                           cc09(n)*CONJG(ca09(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb09(n)*CONJG(cc09(m)) + &
                           cc09(n)*CONJG(cb09(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc09(n)*CONJG(cc09(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*g       
                      mmn(11,orbi,m,n) = mmn(11,orbi,m,n)  +&
                           (( ca10(n)*CONJG(cc10(m)) + &
                           cc10(n)*CONJG(ca10(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb10(n)*CONJG(cc10(m)) + &
                           cc10(n)*CONJG(cb10(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc10(n)*CONJG(cc10(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*g 
                      mmn(12,orbi,m,n) = mmn(12,orbi,m,n)  +&
                           ( ca11(n)*CONJG(cc11(m)) + &
                           cc11(n)*CONJG(ca11(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb11(n)*CONJG(cc11(m)) + &
                           cc11(n)*CONJG(cb11(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc11(n)*CONJG(cc11(m))*usdus%uloulopn(lo,lo,ityp,jspin) 
                      mmn(13,orbi,m,n) = mmn(13,orbi,m,n)  +&
                           (( ca12(n)*CONJG(cc12(m)) + &
                           cc12(n)*CONJG(ca12(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb12(n)*CONJG(cc12(m)) + &
                           cc12(n)*CONJG(cb12(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc12(n)*CONJG(cc12(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*g
                      mmn(14,orbi,m,n) = mmn(14,orbi,m,n)  +&
                           (( ca13(n)*CONJG(cc13(m)) + &
                           cc13(n)*CONJG(ca13(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb13(n)*CONJG(cc13(m)) + &
                           cc13(n)*CONJG(cb13(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc13(n)*CONJG(cc13(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*h 
                      mmn(15,orbi,m,n) = mmn(15,orbi,m,n)  +&
                           (( ca14(n)*CONJG(cc14(m)) + &
                           cc14(n)*CONJG(ca14(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb14(n)*CONJG(cc14(m)) + &
                           cc14(n)*CONJG(cb14(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc14(n)*CONJG(cc14(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*g
                      mmn(16,orbi,m,n) = mmn(16,orbi,m,n)  +&
                           (( ca15(n)*CONJG(cc15(m)) + &
                           cc15(n)*CONJG(ca15(m)) )*usdus%uulon(lo,ityp,jspin) +&
                           ( cb15(n)*CONJG(cc15(m)) + &
                           cc15(n)*CONJG(cb15(m)) )*usdus%dulon(lo,ityp,jspin) +&
                           cc15(n)*CONJG(cc15(m))*usdus%uloulopn(lo,lo,ityp,jspin) )*h
                   ENDDO !m   
                ENDDO !n    
             ENDIF
          ENDDO lo_loop


       ENDDO !rep   

    ENDDO     ! types (ityp)

  END SUBROUTINE wann_orbcomp
END MODULE m_wann_orbcomp
