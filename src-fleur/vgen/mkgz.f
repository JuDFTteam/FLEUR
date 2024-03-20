      MODULE m_mkgz
c.....------------------------------------------------------------------
c     by use of cartesian x,y,z components of charge density gradients,
c     make the quantities
cc      agrt,agru,agrd,g2rt,g2ru,g2rd,gggrt,gggru,gggrd,
cc      gzgr
cc    used to calculate gradient contribution to xc potential and
cc    energy.
c.....------------------------------------------------------------------
      CONTAINS
      SUBROUTINE mkgz(
     >     nmzdf,jspins,rh1,rh2,rhdz1,rhdz2,rhdzz1,rhdzz2,idx,
     <     grad)
      USE m_types
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: nmzdf,jspins,idx
      REAL,    INTENT (IN) :: rh1(nmzdf),rhdz1(nmzdf),rhdzz1(nmzdf)
      REAL,    INTENT (IN) :: rh2(nmzdf),rhdz2(nmzdf),rhdzz2(nmzdf)
      TYPE(t_gradients),INTENT(INOUT)::grad

      INTEGER i
      REAL    vlt,dvzt,dvzzt,vlu,dvzu,dvzzu,vld,dvzd,dvzzd
      REAL    dagrzt,dagrzu,dagrzd,dztadz,sml

      sml = 1.e-14

      IF (ALLOCATED(grad%sigma)) THEN
         IF(jspins==1) THEN
            DO i=1,nmzdf
               grad%sigma(1,idx+i)=rhdz1(i)*rhdz1(i)
            ENDDO
         ELSE
             DO i=1,nmzdf
               grad%sigma(1,idx+i)=rhdz1(i)*rhdz1(i)
               grad%sigma(2,idx+i)=rhdz1(i)*rhdz2(i)
               grad%sigma(3,idx+i)=rhdz2(i)*rhdz2(i)
            ENDDO
         ENDIF
         RETURN
      ENDIF

      IF (jspins == 1) THEN

        DO i = 1,nmzdf

          vlu = max(rh1(i)/2,sml)
          dvzu = rhdz1(i)/2
          dvzzu = rhdzz1(i)/2
          vld = vlu
          dvzd = dvzu
          dvzzd = dvzzu

          vlt = vlu+vld
          dvzt = dvzu+dvzd
          dvzzt = dvzzu+dvzzd

c         agrt: abs(grad(ro)), u,d for up and down.

          grad%agrt(idx+i) = max(abs(dvzt),sml)
          grad%agru(idx+i) = max(abs(dvzu),sml)
          grad%agrd(idx+i) = max(abs(dvzd),sml)

          dagrzt= dvzt*dvzzt/grad%agrt(idx+i)
          dagrzu= dvzu*dvzzu/grad%agru(idx+i)
          dagrzd= dvzd*dvzzd/grad%agrd(idx+i)

          grad%gggrt(idx+i) = dvzt*dagrzt
          grad%gggru(idx+i) = dvzu*dagrzu
          grad%gggrd(idx+i) = dvzd*dagrzd

c         dztadz=d(zeta)/dz,..

          dztadz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2

c         gzgr=grad(zeta)*grad(ro).

          grad%gzgr(idx+i) = dztadz*dvzt

c         g2rt: grad(grad(ro))

          grad%g2rt(idx+i)  = dvzzt
          grad%g2ru(idx+i)  = dvzzu
          grad%g2rd(idx+i)  = dvzzd

        ENDDO

      ELSE

        DO i = 1,nmzdf

          vlu = max(rh1(i),sml)
          dvzu = rhdz1(i)
          dvzzu = rhdzz1(i)
          vld = max(rh2(i),sml)
          dvzd = rhdz2(i)
          dvzzd = rhdzz2(i)

          vlt = vlu+vld
          dvzt = dvzu+dvzd
          dvzzt = dvzzu+dvzzd

c         agrt: abs(grad(ro)), u,d for up and down.

          grad%agrt(idx+i) = max(abs(dvzt),sml)
          grad%agru(idx+i) = max(abs(dvzu),sml)
          grad%agrd(idx+i) = max(abs(dvzd),sml)

          dagrzt= dvzt*dvzzt/grad%agrt(idx+i)
          dagrzu= dvzu*dvzzu/grad%agru(idx+i)
          dagrzd= dvzd*dvzzd/grad%agrd(idx+i)

          grad%gggrt(idx+i) = dvzt*dagrzt
          grad%gggru(idx+i) = dvzu*dagrzu
          grad%gggrd(idx+i) = dvzd*dagrzd

c         dztadz=d(zeta)/dz,..

          dztadz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2

c         gzgr=grad(zeta)*grad(ro).

          grad%gzgr(idx+i) = dztadz*dvzt

c         g2rt: grad(grad(ro))

          grad%g2rt(idx+i)  = dvzzt
          grad%g2ru(idx+i)  = dvzzu
          grad%g2rd(idx+i)  = dvzzd


        ENDDO
      ENDIF

      END SUBROUTINE mkgz
      END MODULE m_mkgz
