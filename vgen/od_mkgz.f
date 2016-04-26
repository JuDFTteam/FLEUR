      MODULE m_od_mkgz
      CONTAINS
      SUBROUTINE od_mkgz(
     >     z1,nmzxy,delz,
     >     nmzdf,jspins,rh1,rh2,rhdz1,rhdz2,rhdzz1,rhdzz2,
     <     agrt,agru,agrd,g2rt,g2ru,g2rd,gggrt,gggru,gggrd,
     <     gzgr)

c     by use of cartesian cylindrical z,phi,r components of charge 
c     density gradients, makes the quantities
c     agrt,agru,agrd,g2rt,g2ru,g2rd,gggrt,gggru,gggrd,
c     gzgr
c     used to calculate gradient contribution to xc potential and
c     energy.
c     for more comments on variables go to mkgxyz3_1.F and vvacxcg_1.F
c     Y.Mokrousov, 16-17 Oct 2002

      implicit none
      
      integer, intent (in) :: nmzdf,jspins,nmzxy
      real,    intent (in) :: z1,delz
      real,    intent (in) :: rh1(nmzdf),rh2(nmzdf) 
      real,    intent (in) :: rhdz1(nmzdf),rhdz2(nmzdf)
      real,    intent (in) ::  rhdzz1(nmzdf),rhdzz2(nmzdf)

      real,   intent (out) :: agrt(nmzdf),agru(nmzdf),agrd(nmzdf)
      real,   intent (out) :: g2rt(nmzdf),g2ru(nmzdf),g2rd(nmzdf)
      real,   intent (out) :: gggrt(nmzdf),gggru(nmzdf),gggrd(nmzdf)
      real,   intent (out) :: gzgr(nmzdf)

c     local

      integer :: i
      real    :: vlt,vlu,vld
      real    :: dvzt,dvzzt,dvzu
      real    :: dvzzu,dvzd,dvzzd
      real    :: dagrzt,dagrzu,dagrzd,dztadz
      real    :: sml,z

      sml = 1.e-14

      if (jspins.eq.1) then
         do 10 i = 1,nmzdf

            z = z1 + delz*(nmzxy+i)

            vlu = max(rh1(i)/2,sml)
            dvzu = rhdz1(i)/2
            dvzzu = rhdzz1(i)/2
            vld = vlu
            dvzd = dvzu
            dvzzd = dvzzu

            vlt = vlu+vld
            dvzt = dvzu+dvzd
            dvzzt = dvzzu+dvzzd

c     agr(up,down,total): abs(grad(rho))

            agrt(i) = max(abs(dvzt),sml)
            agru(i) = max(abs(dvzu),sml)
            agrd(i) = max(abs(dvzd),sml)

c     d(abs(grad(rho)))/dr

            dagrzt= dvzt*dvzzt/agrt(i)
            dagrzu= dvzu*dvzzu/agru(i)
            dagrzd= dvzd*dvzzd/agrd(i)

c     grad(rho)*grad(abs(grad(ro)))

            gggrt(i) = dvzt*dagrzt
            gggru(i) = dvzu*dagrzu
            gggrd(i) = dvzd*dagrzd

c     dztadz=d(zeta)/dz,..
            
            dztadz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2
            
c     gzgr=grad(zeta)*grad(ro).
            
            gzgr(i) = dztadz*dvzt
            
c     g2rt: grad(grad(ro))
            
            g2rt(i)  = dvzzt + dvzt/z
            g2ru(i)  = dvzzu + dvzu/z
            g2rd(i)  = dvzzd + dvzd/z


 10    continue

      else
         
         do 20 i = 1,nmzdf

            z = z1 + delz*(nmzxy + i)
            
            vlu = max(rh1(i),sml)
            dvzu = rhdz1(i)
            dvzzu = rhdzz1(i)
            vld = max(rh2(i),sml)
            dvzd = rhdz2(i)
            dvzzd = rhdzz2(i)
            
            vlt = vlu+vld
            dvzt = dvzu+dvzd
            dvzzt = dvzzu+dvzzd
            
c     agrt: abs(grad(ro)), u,d for up and down.
            
            agrt(i) = max(abs(dvzt),sml)
            agru(i) = max(abs(dvzu),sml)
            agrd(i) = max(abs(dvzd),sml)
            
            dagrzt= dvzt*dvzzt/agrt(i)
            dagrzu= dvzu*dvzzu/agru(i)
            dagrzd= dvzd*dvzzd/agrd(i)
            
            gggrt(i) = dvzt*dagrzt
            gggru(i) = dvzu*dagrzu
            gggrd(i) = dvzd*dagrzd
            
c     dztadz=d(zeta)/dz,..
            
            dztadz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2
            
c     gzgr=grad(zeta)*grad(ro).
            
            gzgr(i) = dztadz*dvzt
            
c     g2rt: grad(grad(ro))
            
            g2rt(i)  = dvzzt + dvzt/z
            g2ru(i)  = dvzzu + dvzu/z
            g2rd(i)  = dvzzd + dvzd/z
            
            
 20      continue

      ENDIF
      
      END SUBROUTINE od_mkgz
      END MODULE m_od_mkgz
