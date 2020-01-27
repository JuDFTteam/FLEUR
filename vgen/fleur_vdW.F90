!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_fleur_vdW
  IMPLICIT NONE
  PUBLIC fleur_vdW,priv_fleur_vdW
CONTAINS
  SUBROUTINE fleur_vdW(mpi,atoms,sphhar,stars,input,DIMENSION,      &
       cell,sym,oneD,vacuum,    &
       qpw,rho,vpw_total,vr_total)
    !Interface to Juelich vdW-code
    USE m_types
    USE m_psqpw
    USE m_fft3d
    USE m_qpwtonmt
    USE m_convol
    USE m_cdn_io
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_oneD),INTENT(IN)      :: oneD
    COMPLEX,INTENT(in)     :: qpw(:)
    REAL,INTENT(inout)     :: rho(:,:,:)
    COMPLEX,INTENT(inout)  :: vpw_total(:)
    REAL,INTENT(inout)     :: vr_total(:,:,:,:)


    !locals
    TYPE(t_atoms)      :: atoms_tmp
    REAL               :: e_vdW
    REAL,ALLOCATABLE   :: n_grid(:),v_grid(:),rhc(:,:,:)
    COMPLEX,ALLOCATABLE:: vpw(:),psq(:)
    INTEGER            :: n,ncmsh,j,i
    LOGICAL            :: l_core,l_pot
    REAL tec(atoms%ntype,input%jspins),qintc(atoms%ntype,input%jspins)


    l_core=.FALSE. !try to subtract core charge?
    ALLOCATE(n_grid(27*stars%mx1*stars%mx2*stars%mx3),v_grid(27*stars%mx1*stars%mx2*stars%mx3))
    ALLOCATE(vpw(SIZE(qpw)))
    ALLOCATE(psq(SIZE(qpw)),rhc(atoms%jmtd,atoms%ntype,input%jspins))

    IF (l_core) l_core = isCoreDensityPresent()

    IF (l_core) THEN
       WRITE(6,*) "VdW contribution without core charge"       
       ! read the core charge
       CALL readCoreDensity(input,atoms,dimension,rhc,tec,qintc)
       DO j=1,input%jspins
          DO n=1,atoms%ntype
             ncmsh = NINT( LOG( (atoms%rmt(n)+10.0)/atoms%rmsh(1,n) ) / atoms%dx(n) + 1 )
             ncmsh = MIN( ncmsh, DIMENSION%msh )
             rho(:,1,n) = rho(:,1,n) - rhc(:SIZE(rho,1),n,j)/(4. * SQRT( ATAN (1.) ))
          ENDDO
       ENDDO
    ENDIF

    ! Construct the pseudo charge
    atoms_tmp=atoms
    atoms_tmp%zatom=0.0
    CALL psqpw(mpi,&
         atoms_tmp,sphhar,stars,vacuum,&
         cell,input,sym,oneD,&
         qpw,rho,(/0.,0./),.TRUE.,2,psq)

    !put pseudo charge on real-space grid
    !use v_pot for imaginary part
    CALL fft3d(                        &
         n_grid,v_grid,psq,         &
         stars,1)


    CALL priv_fleur_vdW(cell,stars, &
         n_grid,e_vdW,v_grid,.TRUE.)

    WRITE(6,*) "------  vdW-Potential code by M. Callsen included-------"
    WRITE(6,*) "vdW-Energy contribution:",e_vdW


    INQUIRE(file="vdW_sc",exist=l_pot)

    IF (.NOT.l_pot) RETURN


    !Put potential on rez. grid
    n_grid=0.0
    CALL fft3d(                           &
         v_grid,n_grid,vpw,            &
         stars,-1)

    !Calculate MT-contribution to the potential

    CALL qpw_to_nmt(                                                     &
         sphhar,atoms,stars,sym,cell,oneD,mpi,  &
         1,4,vpw,vr_total)

    WRITE(6,*) "vdW average Potential  :",vpw(1)


    CALL convol(                    &
         stars,   &
         psq,               &
         vpw,stars%ufft)

    ! Add to total potential
    vpw_total(:)=vpw_total(:)+psq

  END SUBROUTINE fleur_vdW



  SUBROUTINE priv_fleur_vdW(cell,stars, &
       n_pseudo,e_vdw,v_vdw,l_vdW_v1)
    USE m_types
    USE m_constants,ONLY: pi_const
    USE m_juDFT
    USE param, ONLY:  Zab_v1,Zab_v2

    USE nonlocal_data, ONLY: nx,ny,nz,                &
         n_grid,                  &
         a1,a2,a3,                &
         b1,b2,b3,                &
         G_cut,Zab,               &
         lambda,m_c,omega,tpibya
    USE nonlocal_funct,ONLY: soler
    IMPLICIT NONE
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_stars),INTENT(IN) :: stars
    REAL,INTENT(inout)       :: n_pseudo(:)
    LOGICAL,INTENT(in)       :: l_vdW_v1
    REAL,INTENT(out)         :: e_vdw
    REAL,INTENT(out)         :: v_vdw(:)

    IF (SIZE(n_pseudo).NE.SIZE(v_vdw)) CALL juDFT_error("BUG in fleur_to_vdW")

    a1=cell%amat(:,1)
    a2=cell%amat(:,2)
    a3=cell%amat(:,3)

    b1=cell%bmat(:,1)
    b2=cell%bmat(:,2)
    b3=cell%bmat(:,3)

    nx=3*stars%mx1
    ny=3*stars%mx2
    nz=3*stars%mx3
    n_grid=nx*ny*nz
    IF (SIZE(n_pseudo).NE.n_grid) CALL juDFT_error("BUG2 in fleur_to_vdW")
    g_cut=9*stars%gmax**2  !????

    IF (l_vdW_v1) THEN
       Zab=Zab_v1
    ELSE
       Zab=Zab_v2
    ENDIF

    lambda=1.2
    m_c=12
    omega=cell%vol
    tpibya = 2.0*pi_const/omega

    WRITE(6,*)
    WRITE(6,'(A)')      'lattice vectors in Bohr:'
    WRITE(6,'(3F16.8)') a1(:)
    WRITE(6,'(3F16.8)') a2(:)
    WRITE(6,'(3F16.8)') a3(:)
    WRITE(6,*)


    WRITE(6,'(A)')     'reciprocal lattice vectors in Bohr:'
    WRITE(6,'(3F16.8)') b1(:)
    WRITE(6,'(3F16.8)') b2(:)
    WRITE(6,'(3F16.8)') b3(:)
    WRITE(6,'(3F16.8)')
    !
    WRITE(6,'(A,F18.12)') '(2 pi/a) in 1/Bohr^3:',tpibya
    WRITE(6,'(A,F18.12)') 'omega in 1/Bohr^3:   ',omega
    WRITE(6,*)
    WRITE(6,'(A,F18.12)') 'G_cut: ',G_cut


    CALL soler(n_pseudo,e_vdW,v_vdW)

  END SUBROUTINE priv_fleur_vdW

  SUBROUTINE test_charge(qpw,stars)
    USE m_types
    IMPLICIT NONE
    TYPE(t_stars),INTENT(IN) :: stars
    COMPLEX,INTENT(out)      :: qpw(:)
    REAL,SAVE:: pos=0.1
    COMPLEX,ALLOCATABLE,SAVE::testcharge(:)
    COMPLEX:: ph
    INTEGER:: n
    IF (pos==0.1) THEN
       ALLOCATE(testcharge(SIZE(qpw)))
       DO n=1,SIZE(qpw)
          testcharge(n)=EXP(-.1*DOT_PRODUCT(stars%kv3(:,n),stars%kv3(:,n)))
       ENDDO
    ENDIF
    DO n=1,SIZE(qpw)
       ph=EXP(CMPLX(0,stars%kv3(3,n)*pos))
       qpw(n)=(ph+CONJG(ph))*testcharge(n)
    ENDDO
    pos=pos+0.01
    IF (pos>0.4) STOP "testcharge"

  END SUBROUTINE test_charge
END MODULE m_fleur_vdW


