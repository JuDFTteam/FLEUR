MODULE m_orbcomp
CONTAINS
  SUBROUTINE orb_comp(jspin,nobd,atoms,ne,usdus, acof,bcof,ccof, orbcomp,qmtp)
    !***********************************************************************
    !     Calculates an orbital composition of eigen states
    !     
    !                                   Yury  Koroteev  2003-12-24 
    !***********************************************************************
    !                     ABBREVIATIONS
    !          dimentions
    ! nobd                  : in, number of considered bands   
    ! lmd                   : in, (lmaxd + 1)**2
    ! natd                  : in, number of atoms in a film
    ! lmaxd                 : in, max of l    
    ! ntypd                 : in, number of mt-sphere types
    ! nlod                  : in, number of local orbitals in mt-sphere types
    ! llod                  : in, l max for local orbitals in mt-sphere types
    ! ---------------------------------------------------------------------- 
    ! neq(ntypd)            : in, number of mt-spheres of the same type
    ! acof(nobd,0:lmd,natd) : in, a,b  coefficients of linearized 
    ! bcof(nobd,0:lmd,natd) : in, mt-wavefunctions for each band and atom
    ! ccof(-llod:llod,nobd, :
    !     :      nobd,natd) : in, c coefficients for local orbitals
    ! ddn(16,ntypd)         : in,   
    ! uulon(16,ntypd)       : in,   
    ! dulon(16,ntypd)       : in,  
    ! uloulopn(16,ntypd)    : in,   
    ! nlo(ntypd)            : in, 
    ! llo(nlod,ntypd)       : in,
    !-----------------------------------------------------------------------
    ! orbcomp(nobd,16,natd) : out, an orbital composition of  states
    ! qmtp(nobd,natd)       : out, the portion of the state in mt-sphere
    !-----------------------------------------------------------------------
    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_usdus),INTENT(IN)   :: usdus
    !	..
    !	..Scalar Argument
    INTEGER, INTENT  (IN) :: nobd,ne,jspin
    !	..
    !	..Array Arguments
    COMPLEX, INTENT  (IN) :: acof(:,0:,:) !(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT  (IN) :: bcof(:,0:,:) !nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT  (IN) :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%nat)
    REAL,    INTENT (OUT) :: orbcomp(:,:,:)!(dimension%neigd,23,atoms%nat) 
    REAL,    INTENT (OUT) :: qmtp(:,:)!(dimension%neigd,atoms%nat) 
    !	..
    !	..Local Scalars 
    INTEGER  n,mt,ityp,imt,lm,lo
    INTEGER  l,lme,nate,lmaxe,jspe,nobc,nei
    REAL     sum,cf
    REAL     ddn0,ddn1,ddn2,ddn3,ddn12,ddn22,ddn32
    COMPLEX  ca00,ca01,ca02,ca03,ca04,ca05,ca06,ca07,ca08,ca09
    COMPLEX  ca10,ca11,ca12,ca13,ca14,ca15,ca16,ca17,ca18,ca19
    COMPLEX  ca20,ca21,ca22
    COMPLEX  cb00,cb01,cb02,cb03,cb04,cb05,cb06,cb07,cb08,cb09
    COMPLEX  cb10,cb11,cb12,cb13,cb14,cb15,cb16,cb17,cb18,cb19
    COMPLEX  cb20,cb21,cb22
    COMPLEX  cc00,cc01,cc02,cc03,cc04,cc05,cc06,cc07,cc08,cc09
    COMPLEX  cc10,cc11,cc12,cc13,cc14,cc15,cc16,cc17,cc18,cc19
    COMPLEX  cc20,cc21,cc22
    COMPLEX  ck00,ck01,ck02,ck03,ck04,ck05,ck06,ck07,ck08,ck09
    COMPLEX  ck10,ck11,ck12,ck13,ck14,ck15,ck16,ck17,ck18,ck19
    COMPLEX  ck20,ck21,ck22
    !	..
    !	..Local Arrays 
    REAL     comp(23)
    !	..
    !
    REAL,PARAMETER :: h=0.50, g=0.0625
    !**************************************************** 
    !
    mt=0
    DO   ityp = 1,atoms%ntype
       ddn0 = usdus%ddn(0,ityp,jspin)
       ddn1 = usdus%ddn(1,ityp,jspin) 	 
       ddn2 = usdus%ddn(2,ityp,jspin) 	 
       ddn3 = usdus%ddn(3,ityp,jspin) 
       DO  imt=1,atoms%neq(ityp)
          mt=mt+1
          DO  n=1,ne
             !
             ! acof
             !   s-states
             ca00 = acof(n,0,mt)
             !   p-states
             ca01 = acof(n,1,mt) - acof(n,3,mt)
             ca02 = acof(n,1,mt) + acof(n,3,mt)
             ca03 = acof(n,2,mt)
             !   d-states
             ca04 = acof(n,4,mt) - acof(n,8,mt)
             ca05 = acof(n,5,mt) + acof(n,7,mt)
             ca06 = acof(n,5,mt) - acof(n,7,mt)
             ca07 = acof(n,4,mt) + acof(n,8,mt)
             ca08 = acof(n,6,mt)
             !
             !   f-states: a cubic set (cub) 
             ! 
             ca09 = ( acof(n,9,mt)  - acof(n,15,mt) )*SQRT(5.0) - ( acof(n,11,mt) - acof(n,13,mt) )*SQRT(3.0)
             ca10 = ( acof(n,9,mt)  + acof(n,15,mt) )*SQRT(5.0) + ( acof(n,11,mt) + acof(n,13,mt) )*SQRT(3.0) 
             ca11 =   acof(n,12,mt)
             ca12 = ( acof(n,9,mt)  + acof(n,15,mt) )*SQRT(3.0) - ( acof(n,11,mt) + acof(n,13,mt) )*SQRT(5.0) 
             ca13 =   acof(n,10,mt) + acof(n,14,mt)
             ca14 = ( acof(n,9,mt)  - acof(n,15,mt) )*SQRT(3.0) + ( acof(n,11,mt) - acof(n,13,mt) )*SQRT(5.0) 
             ca15 =   acof(n,10,mt) - acof(n,14,mt) 
             !
             !   f-states:	a low symmetry set (lss)
             !
             ca16 =  acof(n,11,mt) - acof(n,13,mt)
             ca17 =  acof(n,11,mt) + acof(n,13,mt)
             ca18 =  acof(n,12,mt)
             ca19 =  acof(n,10,mt) - acof(n,14,mt)
             ca20 =  acof(n,10,mt) + acof(n,14,mt)
             ca21 =  acof(n,9,mt)  - acof(n,15,mt)
             ca22 =  acof(n,9,mt)  + acof(n,15,mt)
             !
             ! bcof
             !   s-states
             cb00 =  bcof(n,0,mt)
             !   p-states
             cb01 =  bcof(n,1,mt) - bcof(n,3,mt) 
             cb02 =  bcof(n,1,mt) + bcof(n,3,mt) 
             cb03 =  bcof(n,2,mt)
             !   d-states
             cb04 =  bcof(n,4,mt) - bcof(n,8,mt) 
             cb05 =  bcof(n,5,mt) + bcof(n,7,mt) 
             cb06 =  bcof(n,5,mt) - bcof(n,7,mt) 
             cb07 =  bcof(n,4,mt) + bcof(n,8,mt) 
             cb08 =  bcof(n,6,mt)
             !
             !   f-states: a cubic set (cub)
             !
             cb09 = ( bcof(n,9,mt)  - bcof(n,15,mt) )*SQRT(5.0) - ( bcof(n,11,mt) - bcof(n,13,mt) )*SQRT(3.0)
             cb10 = ( bcof(n,9,mt)  + bcof(n,15,mt) )*SQRT(5.0) + ( bcof(n,11,mt) + bcof(n,13,mt) )*SQRT(3.0) 
             cb11 =   bcof(n,12,mt)
             cb12 = ( bcof(n,9,mt)  + bcof(n,15,mt) )*SQRT(3.0) - ( bcof(n,11,mt) + bcof(n,13,mt) )*SQRT(5.0) 
             cb13 =   bcof(n,10,mt) + bcof(n,14,mt)
             cb14 = ( bcof(n,9,mt)  - bcof(n,15,mt) )*SQRT(3.0) + ( bcof(n,11,mt) - bcof(n,13,mt) )*SQRT(5.0)
             cb15 =   bcof(n,10,mt) - bcof(n,14,mt) 
             !
             !   f-states:	a low symmetry set (lss)
             !
             cb16 =  bcof(n,11,mt) - bcof(n,13,mt)
             cb17 =  bcof(n,11,mt) + bcof(n,13,mt)
             cb18 =  bcof(n,12,mt)
             cb19 =  bcof(n,10,mt) - bcof(n,14,mt)
             cb20 =  bcof(n,10,mt) + bcof(n,14,mt)
             cb21 =  bcof(n,9,mt)  - bcof(n,15,mt)
             cb22 =  bcof(n,9,mt)  + bcof(n,15,mt)
             !------------------------------------------------------------------
             !  s
             comp(1)  =   ca00*CONJG(ca00) + cb00*CONJG(cb00)*ddn0 
             !  p
             comp(2)  = ( ca01*CONJG(ca01) + cb01*CONJG(cb01)*ddn1 )*h
             comp(3)  = ( ca02*CONJG(ca02) + cb02*CONJG(cb02)*ddn1 )*h
             comp(4)  =   ca03*CONJG(ca03) + cb03*CONJG(cb03)*ddn1 
             !  d
             comp(5)  = ( ca04*CONJG(ca04) + cb04*CONJG(cb04)*ddn2 )*h
             comp(6)  = ( ca05*CONJG(ca05) + cb05*CONJG(cb05)*ddn2 )*h
             comp(7)  = ( ca06*CONJG(ca06) + cb06*CONJG(cb06)*ddn2 )*h
             comp(8)  = ( ca07*CONJG(ca07) + cb07*CONJG(cb07)*ddn2 )*h
             comp(9)  =   ca08*CONJG(ca08) + cb08*CONJG(cb08)*ddn2 
             !  f: a cubic set
             comp(10) = ( ca09*CONJG(ca09) + cb09*CONJG(cb09)*ddn3 )*g       
             comp(11) = ( ca10*CONJG(ca10) + cb10*CONJG(cb10)*ddn3 )*g
             comp(12) =   ca11*CONJG(ca11) + cb11*CONJG(cb11)*ddn3 
             comp(13) = ( ca12*CONJG(ca12) + cb12*CONJG(cb12)*ddn3 )*g
             comp(14) = ( ca13*CONJG(ca13) + cb13*CONJG(cb13)*ddn3 )*h
             comp(15) = ( ca14*CONJG(ca14) + cb14*CONJG(cb14)*ddn3 )*g
             comp(16) = ( ca15*CONJG(ca15) + cb15*CONJG(cb15)*ddn3 )*h
             !  f: a low symmetry set
             comp(17) = ( ca16*CONJG(ca16) + cb16*CONJG(cb16)*ddn3 )*h    
             comp(18) = ( ca17*CONJG(ca17) + cb17*CONJG(cb17)*ddn3 )*h
             comp(19) =   ca18*CONJG(ca18) + cb18*CONJG(cb18)*ddn3 
             comp(20) = ( ca19*CONJG(ca19) + cb19*CONJG(cb19)*ddn3 )*h
             comp(21) = ( ca20*CONJG(ca20) + cb20*CONJG(cb20)*ddn3 )*h
             comp(22) = ( ca21*CONJG(ca21) + cb21*CONJG(cb21)*ddn3 )*h
             comp(23) = ( ca22*CONJG(ca22) + cb22*CONJG(cb22)*ddn3 )*h
             !--------------------------------------------------------------------
             ! ccof   ( contributions from local orbitals )
             !
             DO  lo = 1,atoms%nlo(ityp)
                l = atoms%llo(lo,ityp)
                ! lo-s
                IF ( l.EQ.0 )  THEN
	           cc00 = ccof(0,n,lo,mt)
                   ck00 = CONJG(cc00)

                   comp(1)  =  comp(1)  +&
                        ( ca00*ck00 + cc00*CONJG(ca00) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb00*ck00 + cc00*CONJG(cb00) )*usdus%dulon(lo,ityp,jspin) + cc00*ck00*usdus%uloulopn(lo,lo,ityp,jspin) 
	           CYCLE
                ENDIF
                ! lo-p
                IF ( l.EQ.1 )  THEN
	           cc01 = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
	           cc02 = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
	           cc03 = ccof( 0,n,lo,mt)

                   ck01 = CONJG(cc01) 
                   ck02 = CONJG(cc02) 
                   ck03 = CONJG(cc03)
                   !
                   comp(2) = comp(2)  + (( ca01*ck01 + cc01*CONJG(ca01) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb01*ck01 + cc01*CONJG(cb01) )*usdus%dulon(lo,ityp,jspin) + cc01*ck01*usdus%uloulopn(lo,lo,ityp,jspin) )*h 	
	           comp(3) = comp(3)  + (( ca02*ck02 + cc02*CONJG(ca02) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb02*ck02 + cc02*CONJG(cb02) )*usdus%dulon(lo,ityp,jspin) + cc02*ck02*usdus%uloulopn(lo,lo,ityp,jspin) )*h 
                   comp(4) = comp(4)  + ( ca03*ck03 + cc03*CONJG(ca03) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb03*ck03 + cc03*CONJG(cb03) )*usdus%dulon(lo,ityp,jspin) + cc03*ck03*usdus%uloulopn(lo,lo,ityp,jspin) 
	           CYCLE
                ENDIF
                ! lo-d
                IF ( l.EQ.2 )  THEN
	           cc04 = ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
	           cc05 = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
	           cc06 = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
	           cc07 = ccof(-2,n,lo,mt) + ccof(2,n,lo,mt)
	           cc08 = ccof( 0,n,lo,mt)

                   ck04 = CONJG(cc04)
                   ck05 = CONJG(cc05)
                   ck06 = CONJG(cc06)
                   ck07 = CONJG(cc07)
                   ck08 = CONJG(cc08)

                   comp(5) = comp(5)  + (( ca04*ck04 + cc04*CONJG(ca04) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb04*ck04 + cc04*CONJG(cb04) )*usdus%dulon(lo,ityp,jspin) + cc04*ck04*usdus%uloulopn(lo,lo,ityp,jspin) )*h 	
	           comp(6) = comp(6)  + (( ca05*ck05 + cc05*CONJG(ca05) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb05*ck05 + cc05*CONJG(cb05) )*usdus%dulon(lo,ityp,jspin) + cc05*ck05*usdus%uloulopn(lo,lo,ityp,jspin) )*h 
	           comp(7) = comp(7)  + (( ca06*ck06 + cc06*CONJG(ca06) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb06*ck06 + cc06*CONJG(cb06) )*usdus%dulon(lo,ityp,jspin) + cc06*ck06*usdus%uloulopn(lo,lo,ityp,jspin) )*h
     		   comp(8) = comp(8)  + (( ca07*ck07 + cc07*CONJG(ca07) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb07*ck07 + cc07*CONJG(cb07) )*usdus%dulon(lo,ityp,jspin) + cc07*ck07*usdus%uloulopn(lo,lo,ityp,jspin) )*h 	
	           comp(9) = comp(9)  + ( ca08*ck08 + cc08*CONJG(ca08) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb08*ck08 + cc08*CONJG(cb08) )*usdus%dulon(lo,ityp,jspin) + cc08*ck08*usdus%uloulopn(lo,lo,ityp,jspin)  
                   CYCLE				
                ENDIF
                ! lo-f
                IF ( l.EQ.3 )  THEN
                   !
                   !  a cubic set (cub)
                   !
	           cc09 = ( ccof(-3,n,lo,mt) - ccof(3,n,lo,mt) )*SQRT(5.0) - ( ccof(-1,n,lo,mt) - ccof(1,n,lo,mt) )*SQRT(3.0) 
	           cc10 = ( ccof(-3,n,lo,mt) + ccof(3,n,lo,mt) )*SQRT(5.0) + ( ccof(-1,n,lo,mt) + ccof(1,n,lo,mt) )*SQRT(3.0) 
	           cc11 =   ccof( 0,n,lo,mt)
	           cc12 = ( ccof(-3,n,lo,mt) + ccof(3,n,lo,mt) )*SQRT(3.0) - ( ccof(-1,n,lo,mt) + ccof(1,n,lo,mt) )*SQRT(5.0) 
	           cc13 =   ccof(-2,n,lo,mt) + ccof(2,n,lo,mt) 
	           cc14 = ( ccof(-3,n,lo,mt) - ccof(3,n,lo,mt) )*SQRT(3.0) + ( ccof(-1,n,lo,mt) - ccof(1,n,lo,mt) )*SQRT(5.0)
	           cc15 =   ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
            !
                   ck09 = CONJG(cc09)
                   ck10 = CONJG(cc10)
                   ck11 = CONJG(cc11)
                   ck12 = CONJG(cc12)
                   ck13 = CONJG(cc13)
                   ck14 = CONJG(cc14)
                   ck15 = CONJG(cc15)
                   !
                   comp(10) = comp(10)  + (( ca09*ck09 + cc09*CONJG(ca09) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb09*ck09 + cc09*CONJG(cb09) )*usdus%dulon(lo,ityp,jspin) + cc09*ck09*usdus%uloulopn(lo,lo,ityp,jspin) )*g 	
	           comp(11) = comp(11)  + (( ca10*ck10 + cc10*CONJG(ca10) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb10*ck10 + cc10*CONJG(cb10) )*usdus%dulon(lo,ityp,jspin) + cc10*ck10*usdus%uloulopn(lo,lo,ityp,jspin) )*g 
	           comp(12) = comp(12)  + ( ca11*ck11 + cc11*CONJG(ca11) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb11*ck11 + cc11*CONJG(cb11) )*usdus%dulon(lo,ityp,jspin) + cc11*ck11*usdus%uloulopn(lo,lo,ityp,jspin) 
	           comp(13) = comp(13)  + (( ca12*ck12 + cc12*CONJG(ca12) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb12*ck12 + cc12*CONJG(cb12) )*usdus%dulon(lo,ityp,jspin) + cc12*ck12*usdus%uloulopn(lo,lo,ityp,jspin) )*g
	           comp(14) = comp(14)  + (( ca13*ck13 + cc13*CONJG(ca13) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb13*ck13 + cc13*CONJG(cb13) )*usdus%dulon(lo,ityp,jspin) + cc13*ck13*usdus%uloulopn(lo,lo,ityp,jspin) )*h 
	           comp(15) = comp(15)  + (( ca14*ck14 + cc14*CONJG(ca14) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb14*ck14 + cc14*CONJG(cb14) )*usdus%dulon(lo,ityp,jspin) + cc14*ck14*usdus%uloulopn(lo,lo,ityp,jspin) )*g
     		   comp(16) = comp(16)  + (( ca15*ck15 + cc15*CONJG(ca15) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb15*ck15 + cc15*CONJG(cb15) )*usdus%dulon(lo,ityp,jspin) + cc15*ck15*usdus%uloulopn(lo,lo,ityp,jspin) )*h
          !
          !  a low symmetry set (lss)
          !
	           cc16 = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
	           cc17 = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
	           cc18 = ccof( 0,n,lo,mt)
	           cc19 = ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
	           cc20 = ccof(-2,n,lo,mt) + ccof(2,n,lo,mt)
	           cc21 = ccof(-3,n,lo,mt) - ccof(3,n,lo,mt)
	           cc22 = ccof(-3,n,lo,mt) + ccof(3,n,lo,mt)
            !
                   ck16 = CONJG(cc16)
                   ck17 = CONJG(cc17)
                   ck18 = CONJG(cc18)
                   ck19 = CONJG(cc19)
                   ck20 = CONJG(cc20)
                   ck21 = CONJG(cc21)
                   ck22 = CONJG(cc22)
                   !
	           comp(17) = comp(17)  + (( ca16*ck16 + cc16*CONJG(ca16) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb16*ck16 + cc16*CONJG(cb16) )*usdus%dulon(lo,ityp,jspin) + cc16*ck16*usdus%uloulopn(lo,lo,ityp,jspin) )*h 
	           comp(18) = comp(18)  + (( ca17*ck17 + cc17*CONJG(ca17) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb17*ck17 + cc17*CONJG(cb17) )*usdus%dulon(lo,ityp,jspin) + cc17*ck17*usdus%uloulopn(lo,lo,ityp,jspin) )*h
	           comp(19) = comp(19)  + ( ca18*ck18 + cc18*CONJG(ca18) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb18*ck18 + cc18*CONJG(cb18) )*usdus%dulon(lo,ityp,jspin) + cc18*ck18*usdus%uloulopn(lo,lo,ityp,jspin)
	           comp(20) = comp(20)  + (( ca19*ck19 + cc19*CONJG(ca19) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb19*ck19 + cc19*CONJG(cb19) )*usdus%dulon(lo,ityp,jspin) + cc19*ck19*usdus%uloulopn(lo,lo,ityp,jspin) )*h
     		   comp(21) = comp(21)  + (( ca20*ck20 + cc20*CONJG(ca20) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb20*ck20 + cc20*CONJG(cb20) )*usdus%dulon(lo,ityp,jspin) + cc20*ck20*usdus%uloulopn(lo,lo,ityp,jspin) )*h 
	           comp(22) = comp(22)  + (( ca21*ck21 + cc21*CONJG(ca21) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb21*ck21 + cc21*CONJG(cb21) )*usdus%dulon(lo,ityp,jspin) + cc21*ck21*usdus%uloulopn(lo,lo,ityp,jspin) )*h
	           comp(23) = comp(23)  + (( ca22*ck22 + cc22*CONJG(ca22) )*usdus%uulon(lo,ityp,jspin) +&
                        ( cb22*ck22 + cc22*CONJG(cb22) )*usdus%dulon(lo,ityp,jspin) + cc22*ck22*usdus%uloulopn(lo,lo,ityp,jspin) )*h
                ENDIF
             ENDDO
             !-------------------------------------------------------------------
             !    calculate an orbital cnomposition in percets
             !
             sum = 0.0
             DO   lm=1,16
                sum = sum + comp(lm)
             ENDDO
             cf = 100.0/sum 
             qmtp(n,mt) = sum*100.0         
             orbcomp(n,:,mt) = comp(:)*cf
             !----------------------------------------------------
          ENDDO ! bands (n)
       ENDDO    ! atoms (imt) -> mt (=atoms%nat)
    ENDDO       ! types (ityp)
    !
  END SUBROUTINE orb_comp
END MODULE m_orbcomp
