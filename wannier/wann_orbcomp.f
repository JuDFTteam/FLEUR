      MODULE m_wann_orbcomp
c*************************************************************************
c     Compute matrix elements of the first 16 angular momentum projectors
c     within the MT-spheres. This is done separately for each atom type.
c
c     These are:
c     1: |s><s|
c     2: |p_x><p_x|
c     3: |p_y><p_y|
c     4: |p_z><p_z|
c     5: |d_xy><d_xy|
c     6: |d_yz><d_yz|
c     7: |d_zx><d_zx|
c     8: |d_x2-y2><d_x2-y2|
c     9: |d_z2><dz2|
c     10: |f><f|
c     11: |f><f|
c     12: |f><f|
c     13: |f><f|
c     14: |f><f|
c     15: |f><f|
c     16: |f><f|
c
c     Based on code in orb_comp by Yury Koroteev
c
c     Frank Freimuth
c***********************************************************************
      CONTAINS
      SUBROUTINE wann_orbcomp(
     >                  llod,noccbd,nlod,natd,ntypd,lmaxd,lmd,
     >                  ntype,neq,nlo,llo,acof,bcof,ccof,
     >                  ddn,uulon,dulon,uloulopn,
     >                  num_orbs,
     >                  orbs,
     >                  l_oc_f,
     =                  mmn)
      implicit none
c     .. scalar arguments ..
      integer, intent (in) :: llod,nlod,natd,ntypd,lmaxd,lmd
      integer, intent (in) :: ntype,noccbd
c     .. array arguments ..
      integer, intent (in)  :: neq(:)!neq(ntypd)
      integer, intent (in)  :: nlo(:)!nlo(ntypd)
      integer, intent (in)  :: llo(:,:)!llo(nlod,ntypd)
      real,    intent (in)  :: ddn(0:,:)!ddn(0:lmaxd,ntypd)
      real,    intent (in)  :: uloulopn(:,:,:)!uloulopn(nlod,nlod,ntypd)
      real,    intent (in)  :: uulon(:,:)!uulon(nlod,ntypd)
      real,    intent (in)  :: dulon(:,:)!dulon(nlod,ntypd)
      complex, intent (in)  :: ccof(-llod:,:,:,:) !ccof(-llod:llod,noccbd,nlod,natd)
      complex, intent (in)  :: acof(:,0:,:)!acof(noccbd,0:lmd,natd)
      complex, intent (in)  :: bcof(:,0:,:)!bcof(noccbd,0:lmd,natd)
      complex, intent (inout) :: mmn(:,:,:,:)!mmn(16,noccbd,noccbd,ntype)

      integer, intent (in)  :: num_orbs
      integer, intent (in)  :: orbs(:)
      logical, intent (in)  :: l_oc_f

c	..Local Scalars 

	COMPLEX  ca00(noccbd)
        complex  ca01(noccbd)
        complex  ca02(noccbd)
        complex  ca03(noccbd)
        complex  ca04(noccbd)
        complex  ca05(noccbd)
        complex  ca06(noccbd)
        complex  ca07(noccbd)
        complex  ca08(noccbd)
        complex  ca09(noccbd)
        complex  ca10(noccbd)
        complex  ca11(noccbd)
        complex  ca12(noccbd)
        complex  ca13(noccbd)
        complex  ca14(noccbd)
        complex  ca15(noccbd)

	COMPLEX  cb00(noccbd)
        complex  cb01(noccbd)
        complex  cb02(noccbd)
        complex  cb03(noccbd)
        complex  cb04(noccbd)
        complex  cb05(noccbd)
        complex  cb06(noccbd)
        complex  cb07(noccbd)
        complex  cb08(noccbd)
        complex  cb09(noccbd)
        complex  cb10(noccbd)
        complex  cb11(noccbd)
        complex  cb12(noccbd)
        complex  cb13(noccbd)
        complex  cb14(noccbd)
        complex  cb15(noccbd)

	COMPLEX  cc00(noccbd)
        complex  cc01(noccbd)
        complex  cc02(noccbd)
        complex  cc03(noccbd)
        complex  cc04(noccbd)
        complex  cc05(noccbd)
        complex  cc06(noccbd)
        complex  cc07(noccbd)
        complex  cc08(noccbd)
        complex  cc09(noccbd)
        complex  cc10(noccbd)
        complex  cc11(noccbd)
        complex  cc12(noccbd)
        complex  cc13(noccbd)
        complex  cc14(noccbd)
        complex  cc15(noccbd)


	INTEGER ::  lo,orbi,rep
	REAL    :: ddn0,ddn1,ddn2,ddn3
        integer :: mt,ityp,n,m,mt1,l
        real    :: h,g,c3,c5
        logical :: l_doit
c	..
c	..
c	Intrinsic Function
	Intrinsic sqrt,conjg
c
	DATA h/0.50/ g/0.0625/
c**************************************************** 
	c5 = sqrt(5.0)
	c3 = sqrt(3.0)
c
	mt=0
	DO 10  ityp = 1,ntypd
         do rep=1,neq(ityp)
           mt=mt+1 

           l_doit=.false.
           do l=1,num_orbs
              if(orbs(l)==mt)then
                 l_doit=.true.
                 orbi=l
              endif
           enddo
           if(.not.l_doit)cycle

	   ddn0 = ddn(0,ityp)
	   ddn1 = ddn(1,ityp) 	 
	   ddn2 = ddn(2,ityp) 	 
	   ddn3 = ddn(3,ityp) 
           DO n=1,noccbd
c
c acof
c   s-states
	         ca00(n)=acof(n,0,mt)
c   p-states
	         ca01(n) = acof(n,1,mt) - acof(n,3,mt)
	         ca02(n) = acof(n,1,mt) + acof(n,3,mt)
	         ca03(n) = acof(n,2,mt)
c   d-states
	         ca04(n) = acof(n,4,mt) - acof(n,8,mt)
	         ca05(n) = acof(n,5,mt) + acof(n,7,mt)
	         ca06(n) = acof(n,5,mt) - acof(n,7,mt)
	         ca07(n) = acof(n,4,mt) + acof(n,8,mt)
	         ca08(n) = acof(n,6,mt)
c
c   f-states: a cubic set (cub) 
c 
	         ca09(n) = ( acof(n,9,mt)  - acof(n,15,mt) )*c5 -
     -                  ( acof(n,11,mt) - acof(n,13,mt) )*c3
	         ca10(n) = ( acof(n,9,mt)  + acof(n,15,mt) )*c5 +
     +                  ( acof(n,11,mt) + acof(n,13,mt) )*c3 
	         ca11(n) =   acof(n,12,mt)
	         ca12(n) = ( acof(n,9,mt)  + acof(n,15,mt) )*c3 -
     -                  ( acof(n,11,mt) + acof(n,13,mt) )*c5 
	         ca13(n) =   acof(n,10,mt) + acof(n,14,mt)
	         ca14(n) = ( acof(n,9,mt)  - acof(n,15,mt) )*c3 +
     +                  ( acof(n,11,mt) - acof(n,13,mt) )*c5 
	         ca15(n) =   acof(n,10,mt) - acof(n,14,mt) 

c
c bcof
c   s-states
	         cb00(n) =  bcof(n,0,mt)
c   p-states
	         cb01(n) =  bcof(n,1,mt) - bcof(n,3,mt) 
	         cb02(n) =  bcof(n,1,mt) + bcof(n,3,mt) 
	         cb03(n) =  bcof(n,2,mt)
c   d-states
	         cb04(n) =  bcof(n,4,mt) - bcof(n,8,mt) 
	         cb05(n) =  bcof(n,5,mt) + bcof(n,7,mt) 
	         cb06(n) =  bcof(n,5,mt) - bcof(n,7,mt) 
	         cb07(n) =  bcof(n,4,mt) + bcof(n,8,mt) 
	         cb08(n) =  bcof(n,6,mt)
c
c   f-states: a cubic set (cub)
c
	         cb09(n) = ( bcof(n,9,mt)  - bcof(n,15,mt) )*c5 -
     -                  ( bcof(n,11,mt) - bcof(n,13,mt) )*c3
	         cb10(n) = ( bcof(n,9,mt)  + bcof(n,15,mt) )*c5 +
     +                  ( bcof(n,11,mt) + bcof(n,13,mt) )*c3 
	         cb11(n) =   bcof(n,12,mt)
	         cb12(n) = ( bcof(n,9,mt)  + bcof(n,15,mt) )*c3 -
     -                  ( bcof(n,11,mt) + bcof(n,13,mt) )*c5 
	         cb13(n) =   bcof(n,10,mt) + bcof(n,14,mt)
	         cb14(n) = ( bcof(n,9,mt)  - bcof(n,15,mt) )*c3 +
     +                  ( bcof(n,11,mt) - bcof(n,13,mt) )*c5
	         cb15(n) =   bcof(n,10,mt) - bcof(n,14,mt) 
           enddo !n   

           do n=1,noccbd
            do m=1,noccbd
c  s
	 mmn(1,orbi,m,n)  =   ca00(n)*conjg(ca00(m)) 
     +                      + cb00(n)*conjg(cb00(m))*ddn0 
c  p
	 mmn(2,orbi,m,n)  = ( ca01(n)*conjg(ca01(m)) 
     +                      + cb01(n)*conjg(cb01(m))*ddn1 )*h
	 mmn(3,orbi,m,n)  = ( ca02(n)*conjg(ca02(m)) 
     +                      + cb02(n)*conjg(cb02(m))*ddn1 )*h
	 mmn(4,orbi,m,n)  =   ca03(n)*conjg(ca03(m)) 
     +                      + cb03(n)*conjg(cb03(m))*ddn1 
c  d
	 mmn(5,orbi,m,n)  = ( ca04(n)*conjg(ca04(m)) 
     +                      + cb04(n)*conjg(cb04(m))*ddn2 )*h
	 mmn(6,orbi,m,n)  = ( ca05(n)*conjg(ca05(m)) 
     +                      + cb05(n)*conjg(cb05(m))*ddn2 )*h
	 mmn(7,orbi,m,n)  = ( ca06(n)*conjg(ca06(m)) 
     +                      + cb06(n)*conjg(cb06(m))*ddn2 )*h
	 mmn(8,orbi,m,n)  = ( ca07(n)*conjg(ca07(m)) 
     +                      + cb07(n)*conjg(cb07(m))*ddn2 )*h
	 mmn(9,orbi,m,n)  =   ca08(n)*conjg(ca08(m)) 
     +                      + cb08(n)*conjg(cb08(m))*ddn2 
c  f: a cubic set
         if(l_oc_f)then
	 mmn(10,orbi,m,n) = ( ca09(n)*conjg(ca09(m)) 
     +                      + cb09(n)*conjg(cb09(m))*ddn3 )*g       
         mmn(11,orbi,m,n) = ( ca10(n)*conjg(ca10(m)) 
     +                      + cb10(n)*conjg(cb10(m))*ddn3 )*g
	 mmn(12,orbi,m,n) =   ca11(n)*conjg(ca11(m)) 
     +                      + cb11(n)*conjg(cb11(m))*ddn3 
	 mmn(13,orbi,m,n) = ( ca12(n)*conjg(ca12(m)) 
     +                      + cb12(n)*conjg(cb12(m))*ddn3 )*g
	 mmn(14,orbi,m,n) = ( ca13(n)*conjg(ca13(m)) 
     +                      + cb13(n)*conjg(cb13(m))*ddn3 )*h
	 mmn(15,orbi,m,n) = ( ca14(n)*conjg(ca14(m)) 
     +                      + cb14(n)*conjg(cb14(m))*ddn3 )*g
	 mmn(16,orbi,m,n) = ( ca15(n)*conjg(ca15(m)) 
     +                      + cb15(n)*conjg(cb15(m))*ddn3 )*h
         endif !l_oc_f
            enddo !m
           enddo !n 

c--------------------------------------------------------------------
c ccof   ( contributions from local orbitals )
c
	 DO 60 lo = 1,nlo(ityp)
	    l = llo(lo,ityp)
c lo-s
	    IF ( l.EQ.0 )  THEN
                  do n=1,noccbd
	           cc00(n) = ccof(0,n,lo,mt)
                  enddo 
c     
                  do n=1,noccbd
                   do m=1,noccbd
	           mmn(1,orbi,m,n)=mmn(1,orbi,m,n)+
     +                 ( ca00(n)*conjg(cc00(m)) + 
     +                   cc00(n)*conjg(ca00(m)) )*uulon(lo,ityp) + 
     +                 ( cb00(n)*conjg(cc00(m)) + 
     +                   cc00(n)*conjg(cb00(m)) )*dulon(lo,ityp) +
     +             cc00(n)*conjg(cc00(m))*uloulopn(lo,lo,ityp) 
                   enddo !m
                  enddo !n
	          GOTO 60
            ENDIF
c lo-p
	    IF ( l.EQ.1 )  THEN
               do n=1,noccbd
	           cc01(n) = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
	           cc02(n) = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
	           cc03(n) = ccof( 0,n,lo,mt)
               enddo 
               do n=1,noccbd
                do m=1,noccbd
		   mmn(2,orbi,m,n) = mmn(2,orbi,m,n)  +
     +          (   ( ca01(n)*conjg(cc01(m)) + 
     +                cc01(n)*conjg(ca01(m)) )*uulon(lo,ityp) +
     +              ( cb01(n)*conjg(cc01(m)) + 
     +                cc01(n)*conjg(cb01(m)) )*dulon(lo,ityp) +
     +                cc01(n)*conjg(cc01(m))*uloulopn(lo,lo,ityp) )*h 	

		   mmn(3,orbi,m,n) = mmn(3,orbi,m,n)  +
     +          (   ( ca02(n)*conjg(cc02(m)) + 
     +                cc02(n)*conjg(ca02(m)) )*uulon(lo,ityp) +
     +              ( cb02(n)*conjg(cc02(m)) + 
     +                cc02(n)*conjg(cb02(m)) )*dulon(lo,ityp) +
     +                cc02(n)*conjg(cc02(m))*uloulopn(lo,lo,ityp) )*h 	

		   mmn(4,orbi,m,n) = mmn(4,orbi,m,n)  +
     +              ( ca03(n)*conjg(cc03(m)) + 
     +                cc03(n)*conjg(ca03(m)) )*uulon(lo,ityp) +
     +              ( cb03(n)*conjg(cc03(m)) + 
     +                cc03(n)*conjg(cb03(m)) )*dulon(lo,ityp) +
     +                cc03(n)*conjg(cc03(m))*uloulopn(lo,lo,ityp) 


                enddo !m
               enddo !n 
               GOTO 60
	    ENDIF
c lo-d
	    IF ( l.EQ.2 )  THEN
               do n=1,noccbd
	           cc04(n) = ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
	           cc05(n) = ccof(-1,n,lo,mt) + ccof(1,n,lo,mt)
	           cc06(n) = ccof(-1,n,lo,mt) - ccof(1,n,lo,mt)
	           cc07(n) = ccof(-2,n,lo,mt) + ccof(2,n,lo,mt)
	           cc08(n) = ccof( 0,n,lo,mt)
               enddo !n
c
               do n=1,noccbd
                do m=1,noccbd
		   mmn(5,orbi,m,n) = mmn(5,orbi,m,n)  +
     +             (( ca04(n)*conjg(cc04(m)) + 
     +                cc04(n)*conjg(ca04(m)) )*uulon(lo,ityp) +
     +              ( cb04(n)*conjg(cc04(m)) + 
     +                cc04(n)*conjg(cb04(m)) )*dulon(lo,ityp) +
     +                cc04(n)*conjg(cc04(m))*uloulopn(lo,lo,ityp) )*h 
	           mmn(6,orbi,m,n) = mmn(6,orbi,m,n)  +
     +		   (( ca05(n)*conjg(cc05(m)) + 
     +                cc05(n)*conjg(ca05(m)) )*uulon(lo,ityp) +
     +              ( cb05(n)*conjg(cc05(m)) + 
     +                cc05(n)*conjg(cb05(m)) )*dulon(lo,ityp) +
     +                cc05(n)*conjg(cc05(m))*uloulopn(lo,lo,ityp) )*h 
	           mmn(7,orbi,m,n) = mmn(7,orbi,m,n)  +
     +		   (( ca06(n)*conjg(cc06(m)) + 
     +                cc06(n)*conjg(ca06(m)) )*uulon(lo,ityp) +
     +              ( cb06(n)*conjg(cc06(m)) + 
     +                cc06(n)*conjg(cb06(m)) )*dulon(lo,ityp) +
     +                cc06(n)*conjg(cc06(m))*uloulopn(lo,lo,ityp) )*h
     		   mmn(8,orbi,m,n) = mmn(8,orbi,m,n)  +
     +             (( ca07(n)*conjg(cc07(m)) + 
     +                cc07(n)*conjg(ca07(m)) )*uulon(lo,ityp) +
     +              ( cb07(n)*conjg(cc07(m)) + 
     +                cc07(n)*conjg(cb07(m)) )*dulon(lo,ityp) +
     +                cc07(n)*conjg(cc07(m))*uloulopn(lo,lo,ityp) )*h 	
	           mmn(9,orbi,m,n) = mmn(9,orbi,m,n)  +
     +		    ( ca08(n)*conjg(cc08(m)) + 
     +                cc08(n)*conjg(ca08(m)) )*uulon(lo,ityp) +
     +              ( cb08(n)*conjg(cc08(m)) + 
     +                cc08(n)*conjg(cb08(m)) )*dulon(lo,ityp) +
     +                cc08(n)*conjg(cc08(m))*uloulopn(lo,lo,ityp)  
                enddo !m
               enddo !n
               GOTO 60				
	    ENDIF
c lo-f
	    IF ( l.EQ.3 .and. l_oc_f)  THEN
c
c  a cubic set (cub)
c
               do n=1,noccbd
	           cc09(n) = ( ccof(-3,n,lo,mt) - ccof(3,n,lo,mt) )*c5 -
     -                    ( ccof(-1,n,lo,mt) - ccof(1,n,lo,mt) )*c3 
	           cc10(n) = ( ccof(-3,n,lo,mt) + ccof(3,n,lo,mt) )*c5 +
     +                    ( ccof(-1,n,lo,mt) + ccof(1,n,lo,mt) )*c3 
	           cc11(n) =   ccof( 0,n,lo,mt)
	           cc12(n) = ( ccof(-3,n,lo,mt) + ccof(3,n,lo,mt) )*c3 -
     -                    ( ccof(-1,n,lo,mt) + ccof(1,n,lo,mt) )*c5 
	           cc13(n) =   ccof(-2,n,lo,mt) + ccof(2,n,lo,mt) 
	           cc14(n) = ( ccof(-3,n,lo,mt) - ccof(3,n,lo,mt) )*c3 +
     +                    ( ccof(-1,n,lo,mt) - ccof(1,n,lo,mt) )*c5
	           cc15(n) =   ccof(-2,n,lo,mt) - ccof(2,n,lo,mt)
               enddo !n    
               do n=1,noccbd
                do m=1,noccbd
		   mmn(10,orbi,m,n) = mmn(10,orbi,m,n)  +
     +             (( ca09(n)*conjg(cc09(m)) + 
     +                cc09(n)*conjg(ca09(m)) )*uulon(lo,ityp) +
     +              ( cb09(n)*conjg(cc09(m)) + 
     +                cc09(n)*conjg(cb09(m)) )*dulon(lo,ityp) +
     +                cc09(n)*conjg(cc09(m))*uloulopn(lo,lo,ityp) )*g 	
	           mmn(11,orbi,m,n) = mmn(11,orbi,m,n)  +
     +		   (( ca10(n)*conjg(cc10(m)) + 
     +                cc10(n)*conjg(ca10(m)) )*uulon(lo,ityp) +
     +              ( cb10(n)*conjg(cc10(m)) + 
     +                cc10(n)*conjg(cb10(m)) )*dulon(lo,ityp) +
     +                cc10(n)*conjg(cc10(m))*uloulopn(lo,lo,ityp) )*g 
	           mmn(12,orbi,m,n) = mmn(12,orbi,m,n)  +
     +		    ( ca11(n)*conjg(cc11(m)) + 
     +                cc11(n)*conjg(ca11(m)) )*uulon(lo,ityp) +
     +              ( cb11(n)*conjg(cc11(m)) + 
     +                cc11(n)*conjg(cb11(m)) )*dulon(lo,ityp) +
     +                cc11(n)*conjg(cc11(m))*uloulopn(lo,lo,ityp) 
	           mmn(13,orbi,m,n) = mmn(13,orbi,m,n)  +
     +             (( ca12(n)*conjg(cc12(m)) + 
     +                cc12(n)*conjg(ca12(m)) )*uulon(lo,ityp) +
     +              ( cb12(n)*conjg(cc12(m)) + 
     +                cc12(n)*conjg(cb12(m)) )*dulon(lo,ityp) +
     +                cc12(n)*conjg(cc12(m))*uloulopn(lo,lo,ityp) )*g
	           mmn(14,orbi,m,n) = mmn(14,orbi,m,n)  +
     +		   (( ca13(n)*conjg(cc13(m)) + 
     +                cc13(n)*conjg(ca13(m)) )*uulon(lo,ityp) +
     +              ( cb13(n)*conjg(cc13(m)) + 
     +                cc13(n)*conjg(cb13(m)) )*dulon(lo,ityp) +
     +                cc13(n)*conjg(cc13(m))*uloulopn(lo,lo,ityp) )*h 
	           mmn(15,orbi,m,n) = mmn(15,orbi,m,n)  +
     +		   (( ca14(n)*conjg(cc14(m)) + 
     +                cc14(n)*conjg(ca14(m)) )*uulon(lo,ityp) +
     +              ( cb14(n)*conjg(cc14(m)) + 
     +                cc14(n)*conjg(cb14(m)) )*dulon(lo,ityp) +
     +                cc14(n)*conjg(cc14(m))*uloulopn(lo,lo,ityp) )*g
     		   mmn(16,orbi,m,n) = mmn(16,orbi,m,n)  +
     +             (( ca15(n)*conjg(cc15(m)) + 
     +                cc15(n)*conjg(ca15(m)) )*uulon(lo,ityp) +
     +              ( cb15(n)*conjg(cc15(m)) + 
     +                cc15(n)*conjg(cb15(m)) )*dulon(lo,ityp) +
     +                cc15(n)*conjg(cc15(m))*uloulopn(lo,lo,ityp) )*h
                enddo !m   
               enddo !n    
	    ENDIF
  60	 CONTINUE


         enddo !rep   
   
 10   CONTINUE       ! types (ityp)

      END SUBROUTINE wann_orbcomp
      END MODULE m_wann_orbcomp
