MODULE m_pwintsl
CONTAINS
  SUBROUTINE pwint_sl(stars,atoms,sym,zsl1,zsl2, volsl,volintsl, cell,nmtsl1, kv, x)
    !     ******************************************************************
    !     calculate the integral of a star function over the layer 
    !     interstial region of a film                Yury Koroteev  
    !                                   from  pwint.F  by  c.l.fu              
    !     ******************************************************************
    USE m_spgrot
    USE m_constants,ONLY: tpi_const
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_cell),INTENT(IN)    :: cell
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    REAL,    INTENT (IN) :: zsl1,zsl2,volsl,volintsl
    COMPLEX, INTENT (OUT):: x
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: kv(3)  
    INTEGER, INTENT (IN) :: nmtsl1(atoms%ntype) 
    !     ..
    !     .. Local Scalars ..
    COMPLEX s1,sfs
    REAL arg,g,s,srmt,gm,gp,zslm,zslp
    INTEGER ig2d,ig3d,n,nn,nat 
    !     ..
    !     .. Local Arrays ..
    COMPLEX ph(sym%nop)
    INTEGER kr(3,sym%nop)
    !     ..
    ig3d = stars%ig(kv(1),kv(2),kv(3))
    !
    !     -----> interstitial contributions
    !
    IF (ig3d.EQ.0) THEN
       x = (0.,0.)
       RETURN
    END IF
    IF (ig3d.EQ.1) THEN
       x = CMPLX(volintsl,0.0)
       RETURN
    ELSE
       ig2d = stars%ig2(ig3d)
       IF (ig2d.EQ.1) THEN
          zslm = 0.5*(zsl2 - zsl1) 
          zslp = 0.5*(zsl2 + zsl1)
          g = kv(3)*cell%bmat(3,3)
          gm = g*zslm
          gp = g*zslp
          x = volsl*SIN(gm)/gm*CMPLX(COS(gp),SIN(gp))
       ELSE
          x = (0.0,0.0)
       END IF
    END IF
    !
    !     -----> sphere contributions
    !
    s = stars%sk3(ig3d)
    nat = 1
    DO n = 1,atoms%ntype
       srmt = s*atoms%rmt(n)
       CALL spgrot(sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab, kv, kr,ph)
       sfs = (0.0,0.0)
       DO  nn = 1,sym%nop
          arg = tpi_const* dot_product(real(kr(:,nn)),atoms%taual(:,nat))
          sfs = sfs + CMPLX(COS(arg),SIN(arg))*ph(nn)
       ENDDO
       sfs = sfs/sym%nop
       !
       !     -----3*ji(gr)/gr term
       !
       s1 = 3.* (SIN(srmt)/srmt-COS(srmt))/ (srmt*srmt)
       x = x - atoms%volmts(n)*nmtsl1(n)*s1*sfs
       nat = nat + atoms%neq(n)
    ENDDO
    !
  END SUBROUTINE pwint_sl
END MODULE m_pwintsl
