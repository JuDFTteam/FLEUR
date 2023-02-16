MODULE m_uham
c-------------------------------------------------------------------+
c     For details see Eq.(21) of Shick et al. PRB 60, 10765 (1999)  |
c     Sets up the LDA + U Hamilton matrix with eigenvalues in the   |
c     diagonal:                                                     |
c  s     ---       *s,j'  s,j    . .   *s,j'  s,j    s            s |
c H    = >      [ A      A    + <u|u> B      B    ] V     + d    e  |
c  j',j  -- m,m'   l,m'   l,m          l,m'   l,m    m,m'    j,j' j |
c                                                                   |
C                                                  G.B. Oct. 2000   |
c-------------------------------------------------------------------+
      CONTAINS
      SUBROUTINE u_ham(jspd,nvd,lmaxd,ntypd,lmd,matsize,jspins,nv,lda_u,jsp,n
                ,invsfct,ar,ai,br,bi,ddn,vs_mmp,n_u,lmaxb,n_size,n_rank,
     >           l_ss,ab_dim,iintsp,jintsp,aa)

      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,nvd,lmaxd,ntypd,lmaxb,lmd
      INTEGER, INTENT (IN) :: jspins,n_u,jsp,n,matsize
      INTEGER, INTENT (IN) :: n_size,n_rank
      INTEGER, INTENT (IN) :: ab_dim,iintsp,jintsp
      REAL,    INTENT (IN) :: invsfct
      LOGICAL, INTENT (IN) :: l_ss

C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lda_u(ntypd),nv(jspd)
      REAL,    INTENT (IN) :: ddn(0:lmaxd,ntypd,jspd)!(l, atomtypeindexn, spin )
      REAL,    INTENT (IN) :: ar(nvd,0:lmd,ab_dim),ai(nvd,0:lmd,ab_dim)
      REAL,    INTENT (IN) :: br(nvd,0:lmd,ab_dim),bi(nvd,0:lmd,ab_dim)
      COMPLEX,INTENT(IN):: vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,n_u,jspins)

      COMPLEX, INTENT (INOUT) :: aa(matsize)! matsize dimension of hamiltonian matrix.

C     .. Local Scalars ..
      INTEGER m,mp,itype,l,ig,igp,ll,llm,nvv,ii,n_l,igp_m,i_invsf,nc
      COMPLEX chihlp

C     .. Local Arrays
      REAL, ALLOCATABLE :: ar_help(:,:),ai_help(:,:),va_help(:)
      REAL, ALLOCATABLE :: br_help(:,:),bi_help(:,:),vb_help(:)


    c jspd - number of spins (1 or 2)
    c nvd - dimension of augmented plane waves (without local orbitals)
    c lmaxd - dimension of l-quantum number
    c ntypd - dimension of atom types (symmetry-inequivalent atoms) in unit cell
    c lmd - dimension of index of composite variable lm
    c nv(jspd), nv(i) - number of augmented plane waves (without lo) of spin i
    c n_u - This variable is 0 if no LDA+U is used. Otherwise it may contain the number of 'U's.
    c l_ss - .true. if spin spiral calculation is activated
    c rank is the number of linearly independent rows of a matrix.

      i_invsf = nint(invsfct)   ! NINT() rounds its argument to the nearest whole number 
      n_l = 0 ! number of l'S 
      DO itype = 1,n !loop over all atomtypes in the unit cell if l-q# is not negative 
        l = lda_u(itype)
        IF (l.GE.0) n_l = n_l + 1
      ENDDO
      nvv = nv(jsp) ! number of APW of spin jsp.
      IF (l_ss)  nvv = nv(jintsp)
      
      ALLOCATE ( ar_help(nvv,-l:l),ai_help(nvv,-l:l),va_help(-l:l),
     +           br_help(nvv,-l:l),bi_help(nvv,-l:l),vb_help(-l:l) ) !! real ,imaginary, 
c
c-----------------------------------------------------------------------
c
c                __      mm'    *lm'                    __      mm'    *lm' . .
c A_help(G',m) = >     V       A     and B_help(G',m) = >     V       B    <u|u>
c                --m'   l,n,s   G'                      --m'   l,n,s   G'
c
c------------------------------------------------------------------------
      ll = l * (l+1)
      DO m = -l,l
 
        DO mp = -l,l
          va_help(mp) = real( vs_mmp(m,mp,n_l,jsp) ) * invsfct
          vb_help(mp) = va_help(mp) * ddn(l,n,jsp)
        ENDDO

        DO ig = 1, nvv
          ar_help(ig,m) = 0.0
          ai_help(ig,m) = 0.0
          br_help(ig,m) = 0.0
          bi_help(ig,m) = 0.0

          DO mp = -l,l
            ar_help(ig,m) = ar_help(ig,m) + va_help(mp) * 
     +                                      ar(ig,ll+mp,jintsp)
            ai_help(ig,m) = ai_help(ig,m) - va_help(mp) * 
     +                                      ai(ig,ll+mp,jintsp)
            br_help(ig,m) = br_help(ig,m) + vb_help(mp) *
     +                                      br(ig,ll+mp,jintsp)
            bi_help(ig,m) = bi_help(ig,m) - vb_help(mp) *  
     +                                      bi(ig,ll+mp,jintsp)
          ENDDO.

        ENDDO

        DO mp = -l,l
          va_help(mp) = aimag( vs_mmp(m,mp,n_l,jsp) ) * invsfct ! aimag Imaginary part of complex number
          vb_help(mp) = va_help(mp) * ddn(l,n,jsp)
        ENDDO
        DO ig = 1, nvv
          DO mp = -l,l
            ar_help(ig,m) = ar_help(ig,m) + va_help(mp) *
     +                                      ai(ig,ll+mp,jintsp)
            ai_help(ig,m) = ai_help(ig,m) + va_help(mp) *
     +                                      ar(ig,ll+mp,jintsp)
            br_help(ig,m) = br_help(ig,m) + vb_help(mp) *
     +                                      bi(ig,ll+mp,jintsp)
            bi_help(ig,m) = bi_help(ig,m) + vb_help(mp) *
     +                                      br(ig,ll+mp,jintsp)
          ENDDO
        ENDDO

      ENDDO
c

c
c--------------------------------------------
c
c  l,n,s    --        m   lm          m   lm
c H      =  >   A_help   A    + B_help   B
c  G G'     --m       G'  G           G'  G 
c
c--------------------------------------------
      nc = 0
      DO ig = n_rank+1, nvv, n_size
        nc = nc + 1
        chihlp = cmplx(1.0,0.0)
        igp_m = ig
        !ii = ig * (ig - 1) / 2
        ii = nc*(nc-1)/2*n_size - (nc-1)*(n_size-n_rank-1)

        DO m = -l,l
          llm = ll + m

          DO igp = 1, igp_m
#ifdef CPP_INVERSION
            aa(ii+igp) = aa(ii+igp) + 
     +                             ar_help(igp,m) * ar(ig,llm,iintsp) +
     +                             br_help(igp,m) * br(ig,llm,iintsp) -
     +                             ai_help(igp,m) * ai(ig,llm,iintsp) -
     +                             bi_help(igp,m) * bi(ig,llm,iintsp) 
#else
            aa(ii+igp) = aa(ii+igp) + chihlp * cmplx(
     +                             ar_help(igp,m) * ar(ig,llm,iintsp) +
     +                             br_help(igp,m) * br(ig,llm,iintsp) -
     +                             ai_help(igp,m) * ai(ig,llm,iintsp) -
     +                             bi_help(igp,m) * bi(ig,llm,iintsp) ,
     +                             ar_help(igp,m) * ai(ig,llm,iintsp) +
     +                             br_help(igp,m) * bi(ig,llm,iintsp) +
     +                             ai_help(igp,m) * ar(ig,llm,iintsp) +
     +                             bi_help(igp,m) * br(ig,llm,iintsp) )
#endif
          ENDDO ! igp 
        ENDDO   ! m


      END SUBROUTINE u_ham
      END MODULE m_uham