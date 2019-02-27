!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_fix_by_gaussian
  USE m_judft
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fix_by_gaussian(shift,atoms,stars,mpi,sym,vacuum,sphhar,input,oned,cell,noco,den)
    !The idea of this fix is to add an Gaussian to the INT which make the charge flat at the MT-boundary and to 
    !shift this Gaussian with the displacement.
    USE m_qfix
    USE m_spgrot
    USE m_constants
    USE m_types
    IMPLICIT NONE
    REAL,INTENT(in)             :: shift(:,:)
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_potden),INTENT(INOUT):: den

    REAL    :: slope1,slope2,dr,alpha
    REAL    :: sigma,a_fac,gauss,x,fix
    INTEGER :: kr(3,sym%nop)
    COMPLEX :: sf,phas(sym%nop)
    INTEGER :: js,n,l,k,nat,j
    TYPE(t_input)    :: inputtmp

    DO js=1,input%jspins
       DO n=1,atoms%ntype
          DO l=0,0 !Currently only l=0
!             alpha = LOG( den%mt(atoms%jri(n)-1,l,n,js) / den%mt(atoms%jri(n),l,n,js) )
!             alpha = SQRT(alpha / ( atoms%rmt(n)*atoms%rmt(n)*( 1.0-EXP( -2.0*atoms%dx(n) ) ) ))
!             A_fac= den%mt(atoms%jri(n),l,n,js)/gaussian_r(atoms%rmt(n),alpha)
             dr=atoms%rmsh(atoms%jri(n)-1,n)-atoms%rmsh(atoms%jri(n),n)
             slope1=(den%mt(atoms%jri(n)-1,l,n,js)/atoms%rmsh(atoms%jri(n)-1,n)**2-den%mt(atoms%jri(n),l,n,js)/atoms%rmsh(atoms%jri(n),n)**2)/dr
             slope1=(den%mt(atoms%jri(n)-1,l,n,js)-den%mt(atoms%jri(n),l,n,js))/dr/atoms%rmt(n)**2
             sigma=atoms%rmt(n)/2.0
             alpha=1/sigma
             slope2=(gaussian_r(atoms%rmsh(atoms%jri(n)-1,n),alpha)-gaussian_r(atoms%rmsh(atoms%jri(n),n),alpha))/dr
             A_fac=slope1/slope2
             PRINT *, a_fac, 1/alpha
             DO k=2,stars%ng3
                gauss=A_fac*gaussian_g(stars%sk3(k),alpha)
                CALL spgrot(sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,stars%kv3(:,k),kr,phas)
                DO nat=SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n))
                   sf=0.0
                   DO  j = 1,sym%nop
                      x=-tpi_const*DOT_PRODUCT(1.*kr(:,j),atoms%taual(:,nat))
                      sf = sf + CMPLX(COS(x),SIN(x))*CONJG(phas(j))
                      x=-tpi_const*DOT_PRODUCT(1.*kr(:,j),atoms%taual(:,nat)+shift(:,nat))
                      sf = sf - CMPLX(COS(x),SIN(x))*CONJG(phas(j))
                   ENDDO
                ENDDO
                den%pw(k,js)=den%pw(k,js)+gauss*sf/sym%nop/cell%omtil
             ENDDO
          ENDDO
       END DO
    END DO
    inputtmp=input
    inputtmp%qfix=1
    CALL qfix(mpi,stars,atoms,sym,vacuum,sphhar,inputtmp,cell,oneD,&
         den,noco%l_noco,mpi%isize==1,force_fix=.TRUE.,fix=fix)
  END SUBROUTINE fix_by_gaussian

  FUNCTION gaussian_r(r,alpha)
   
    REAL,INTENT(IN) :: r,alpha
    real            :: gaussian_r
    gaussian_r=EXP(-r**2*alpha**2)
  END FUNCTION gaussian_r

  FUNCTION gaussian_g(g,alpha)
    USE m_constants
    REAL,INTENT(IN) :: g,alpha
    REAL            :: gaussian_g
    gaussian_g=SQRT(pi_const**3)/alpha**3*EXP(-0.25*g**2/alpha**2)
  END FUNCTION gaussian_g

END MODULE m_fix_by_gaussian
