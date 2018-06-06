!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
! @authors: Miriam Hinzen, Gregor Michalicek
!--------------------------------------------------------------------------------
MODULE m_chase_diag
#ifdef CPP_CHASE

IMPLICIT NONE

  interface
    subroutine chase_c( h, n, v, ritzv, nev, nex, deg, tol, mode, opt ) bind( c, name = 'zchase_' )
      use, intrinsic :: iso_c_binding
      complex(c_double_complex)     :: h(n,*), v(n,*)
      integer(c_int)                :: n, deg, nev, nex
      real(c_double)                :: ritzv(*), tol
      character(len=1,kind=c_char)  :: mode, opt
    end subroutine chase_c
  end interface

  interface 
    subroutine chase_r( h, n, v, ritzv, nev, nex, deg, tol, mode, opt ) bind( c, name = 'dchase_' )
      use, intrinsic :: iso_c_binding
      real(c_double_complex)        :: h(n,*), v(n,*)
      integer(c_int)                :: n, deg, nev, nex
      real(c_double)                :: ritzv(*), tol
      character(len=1,kind=c_char)  :: mode, opt
    end subroutine chase_r
  end interface

  PRIVATE 

    INTEGER :: chase_eig_id

  PUBLIC init_chase, chase_diag

  CONTAINS

  SUBROUTINE init_chase(mpi,dimension,input,atoms,kpts,noco,vacuum,banddos,l_real)

    USE m_types
    USE m_types_mpi
    USE m_judft
    USE m_eig66_io

    IMPLICIT NONE

    TYPE(t_mpi),               INTENT(IN)    :: mpi
    TYPE(t_dimension),         INTENT(IN)    :: dimension
    TYPE(t_input),             INTENT(IN)    :: input
    TYPE(t_atoms),             INTENT(IN)    :: atoms
    TYPE(t_kpts),              INTENT(IN)    :: kpts
    TYPE(t_noco),              INTENT(IN)    :: noco
    TYPE(t_vacuum),            INTENT(IN)    :: vacuum
    TYPE(t_banddos),           INTENT(IN)    :: banddos

    LOGICAL,                   INTENT(IN)    :: l_real

    INTEGER                                  :: nevd, nexd

    IF (juDFT_was_argument("-diag:chase")) THEN
       nevd = min(dimension%neigd,dimension%nvd+atoms%nlotot)
       nexd = min(max(nevd/4, 45),dimension%nvd+atoms%nlotot-nevd) !dimensioning for workspace
       chase_eig_id=open_eig(mpi%mpi_comm,DIMENSION%nbasfcn,nevd+nexd,kpts%nkpt,DIMENSION%jspd,atoms%lmaxd,&
                             atoms%nlod,atoms%ntype,atoms%nlotot,noco%l_noco,.TRUE.,l_real,noco%l_soc,.FALSE.,&
                             mpi%n_size,layers=vacuum%layers,nstars=vacuum%nstars,ncored=DIMENSION%nstd,&
                             nsld=atoms%nat,nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,l_mcd=banddos%l_mcd,&
                             l_orb=banddos%l_orb)
    END IF
  END SUBROUTINE init_chase

  SUBROUTINE chase_diag(mpi,hmat,smat,ikpt,jsp,iter,ne,eig,zmat)

    USE m_types
    USE m_types_mpi
    USE m_judft
    USE iso_c_binding
    USE m_eig66_io

    !Simple driver to solve Generalized Eigenvalue Problem using the ChASE library
    IMPLICIT NONE

    TYPE(t_mpi),               INTENT(IN)    :: mpi
    TYPE(t_mat),               INTENT(INOUT) :: hmat,smat
    INTEGER,                   INTENT(IN)    :: ikpt
    INTEGER,                   INTENT(IN)    :: jsp
    INTEGER,                   INTENT(IN)    :: iter
    INTEGER,                   INTENT(INOUT) :: ne
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)   :: zmat
    REAL,                      INTENT(OUT)   :: eig(:)

    INTEGER            :: i, j, nev, nex, nbands
    INTEGER            :: info

    CLASS(t_Mat),              ALLOCATABLE  :: zMatTemp
    REAL(c_double),            ALLOCATABLE  :: eigenvalues(:)

    ALLOCATE(t_mat::zmat)
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)

    nev = min(ne,hmat%matsize1)
    nex = min(max(nev/4, 45), hmat%matsize1-nev) !dimensioning for workspace

    ALLOCATE(eigenvalues(nev+nex))
    eigenvalues = 0.0

    ALLOCATE(t_mat::zmatTemp)
    CALL zMatTemp%alloc(hmat%l_real,hmat%matsize1,nev+nex)

    IF (hmat%l_real) THEN

       ! --> start with Cholesky factorization of b ( so that b = l * l^t)
       ! --> b is overwritten by l
       CALL dpotrf('U',smat%matsize1,smat%data_r,SIZE(smat%data_r,1),info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in dpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> now reduce a * z = eig * b * z to the standard form a' * z' = eig * z' 
       ! --> where a' = (l)^-1 * a * (l^t)^-1 and z' = l^t * z
       CALL dsygst(1,'U',smat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),smat%data_r,SIZE(smat%data_r,1),info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in dsygst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub

       zMatTemp%data_r = 0.0

       do j = 1, hmat%matsize1
          do i = 1, j
             hmat%data_r(j,i) = hmat%data_r(i,j)
          end do
       end do
       if(iter.EQ.1) then
          call chase_r(hmat%data_r, hmat%matsize1, zMatTemp%data_r, eigenvalues, nev, nex, 25, 1e-10, 'R', 'S' )
       else
          CALL read_eig(chase_eig_id,ikpt,jsp,n_start=mpi%n_size,n_end=mpi%n_rank,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          call chase_r(hmat%data_r, hmat%matsize1, zMatTemp%data_r, eigenvalues, nev, nex, 25, 1e-10, 'A', 'S' )
       end if

       ne = nev

       CALL write_eig(chase_eig_id,ikpt,jsp,nev+nex,nev+nex,&
                      eigenvalues(:(nev+nex)),n_start=mpi%n_size,n_end=mpi%n_rank,zmat=zMatTemp)

       ! --> recover the generalized eigenvectors z by solving z' = l^t * z
       CALL dtrtrs('U','N','N',hmat%matsize1,nev,smat%data_r,smat%matsize1,zMatTemp%data_r,zmat%matsize1,info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in dtrtrs: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       DO i = 1, ne
          DO j = 1, hmat%matsize1
             zmat%data_r(j,i) = zMatTemp%data_r(j,i)
          END DO
          eig(i) = eigenvalues(i)
       END DO


    ELSE

       ! --> start with Cholesky factorization of b ( so that b = l * l^t)
       ! --> b is overwritten by l
       CALL zpotrf('U',smat%matsize1,smat%data_c,SIZE(smat%data_c,1),info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in zpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> now reduce a * z = eig * b * z to the standard form a' * z' = eig * z' 
       ! --> where a' = (l)^-1 * a * (l^t)^-1 and z' = l^t * z
       CALL zhegst(1,'U',smat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),smat%data_c,SIZE(smat%data_c,1),info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in zhegst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub

       zMatTemp%data_c = CMPLX(0.0,0.0)

       do j = 1, hmat%matsize1
          do i = 1, j
             hmat%data_c(j,i) = conjg(hmat%data_c(i,j))
          end do
       end do

       if(iter.EQ.1) then
          call chase_c(hmat%data_c, hmat%matsize1, zMatTemp%data_c, eigenvalues, nev, nex, 25, 1e-10, 'R', 'S' )
       else
          CALL read_eig(chase_eig_id,ikpt,jsp,n_start=mpi%n_size,n_end=mpi%n_rank,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          call chase_c(hmat%data_c, hmat%matsize1, zMatTemp%data_c, eigenvalues, nev, nex, 25, 1e-10, 'A', 'S' )
       end if

       ne = nev

       CALL write_eig(chase_eig_id,ikpt,jsp,nev+nex,nev+nex,&
                      eigenvalues(:(nev+nex)),n_start=mpi%n_size,n_end=mpi%n_rank,zmat=zMatTemp)

       ! --> recover the generalized eigenvectors z by solving z' = l^t * z
       CALL ztrtrs('U','N','N',hmat%matsize1,nev,smat%data_c,smat%matsize1,zMatTemp%data_c,zmat%matsize1,info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in ztrtrs: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       DO i = 1, ne
          DO j = 1, hmat%matsize1
             zmat%data_c(j,i) = zMatTemp%data_c(j,i)
          END DO
          eig(i) = eigenvalues(i)
       END DO

    ENDIF
    IF (info.NE.0) CALL judft_error("Diagonalization via ChASE failed", calledby = 'chase_diag')
  END SUBROUTINE chase_diag
#endif
END MODULE m_chase_diag
