!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
! @authors: Miriam Hinzen, Gregor Michalicek
! Added MPI implementation, DW 2018
!--------------------------------------------------------------------------------
MODULE m_chase_diag
#ifdef CPP_CHASE
  USE m_judft
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

  !MPI 
  INTERFACE 
     SUBROUTINE mpi_dchase_init( mpi_comm, n, nev, nex, xoff,yoff,xlen,ylen,npr,npc,myrow,mycol) BIND( c, name = 'dchase_init' )
       USE, INTRINSIC :: iso_c_binding
       INTEGER(c_int) :: mpi_comm, n, nev, nex, xoff,yoff,xlen,ylen,myrow,mycol,npr,npc
     END SUBROUTINE mpi_dchase_init
  END INTERFACE

  INTERFACE 
     SUBROUTINE mpi_zchase_init( mpi_comm, n, nev, nex, xoff,yoff,xlen,ylen,npr,npc,myrow,mycol) BIND( c, name = 'zchase_init' )
       USE, INTRINSIC :: iso_c_binding
       INTEGER(c_int) :: mpi_comm, n, nev, nex, xoff,yoff,xlen,ylen,myrow,mycol,npr,npc
     END SUBROUTINE mpi_zchase_init
  END INTERFACE
 
  INTERFACE 
     SUBROUTINE mpi_chase_r(h, v, ritzv, deg, tol, mode, opt ) BIND( c, name = 'dchase_solve' )
       USE, INTRINSIC :: iso_c_binding
       REAL(c_double_complex)        :: h(*), v(*)
       INTEGER(c_int)                :: deg
       REAL(c_double)                :: ritzv(*), tol
       CHARACTER(len=1,kind=c_char)  :: mode, opt
     END SUBROUTINE mpi_chase_r
  END INTERFACE
 

  INTERFACE 
     SUBROUTINE mpi_chase_c(h, v, ritzv, deg, tol, mode, opt ) BIND( c, name = 'zchase_solve' )
       USE, INTRINSIC :: iso_c_binding
       COMPLEX(c_double_complex)     :: h(*), v(*)
       INTEGER(c_int)                :: deg
       REAL(c_double)                :: ritzv(*), tol
       CHARACTER(len=1,kind=c_char)  :: mode, opt
     END SUBROUTINE mpi_chase_c
  END INTERFACE
  

  PRIVATE 

  INTEGER         :: chase_eig_id
  PUBLIC init_chase
#endif
  REAL            :: scale_distance
  REAL            :: tol
  
  PUBLIC chase_distance,chase_diag

CONTAINS

  SUBROUTINE chase_distance(dist)
    IMPLICIT NONE
    REAL,INTENT(in)::dist
    
    tol=MAX(1E-8,dist*scale_distance)
  END SUBROUTINE chase_distance

#ifdef CPP_CHASE
    SUBROUTINE init_chase(mpi,DIMENSION,input,atoms,kpts,noco,l_real)
    USE m_types_mpimat
    USE m_types_setup
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
    LOGICAL,                   INTENT(IN)    :: l_real

    INTEGER                                  :: nevd, nexd
    CHARACTER(len=1000)::arg
    
   
    scale_distance=1E-3
    !IF (judft_was_argument("-chase_tol_scale")) THEN
    !   arg=juDFT_string_for_argument("-chase_tol_scale")
    !   READ(arg,*) scale_distance
    !ENDIF
    
    IF (TRIM(juDFT_string_for_argument("-diag"))=="chase") THEN
       nevd = min(dimension%neigd,dimension%nvd+atoms%nlotot)
       nexd = min(max(nevd/4, 45),dimension%nvd+atoms%nlotot-nevd) !dimensioning for workspace
       chase_eig_id=open_eig(mpi%mpi_comm,DIMENSION%nbasfcn,nevd+nexd,kpts%nkpt,input%jspins,&
                             noco%l_noco,.TRUE.,l_real,noco%l_soc,.FALSE.,mpi%n_size)
    END IF
  END SUBROUTINE init_chase
#endif
  
   SUBROUTINE chase_diag(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)
    USE m_types_mpimat
    USE m_types_mat
    USE m_judft
    USE iso_c_binding
    USE m_eig66_io

    !Simple driver to solve Generalized Eigenvalue Problem using the ChASE library
    IMPLICIT NONE

    CLASS(t_mat),              INTENT(INOUT) :: hmat,smat
    INTEGER,                   INTENT(IN)    :: ikpt
    INTEGER,                   INTENT(IN)    :: jsp
    INTEGER,                   INTENT(IN)    :: iter
    INTEGER,                   INTENT(INOUT) :: ne
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)   :: zmat
    REAL,                      INTENT(OUT)   :: eig(:)
#ifdef CPP_CHASE
    !Choose serial or parallel solver
    SELECT TYPE(hmat)
    CLASS is (t_mpimat)
       SELECT TYPE(smat)
       CLASS is (t_mpimat)
          CALL chase_diag_MPI(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)
       CLASS default
          CALL judft_error("Inconsistent matrix setup")
       END SELECT
    CLASS is (t_mat)
       SELECT TYPE(smat)
       CLASS is (t_mat)
          CALL chase_diag_noMPI(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)
       CLASS default
          CALL judft_error("Inconsistent matrix setup")
       END SELECT
    END SELECT
#endif    
  END SUBROUTINE chase_diag
#ifdef CPP_CHASE  
  SUBROUTINE chase_diag_noMPI(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)

    USE m_types_mat
    USE m_judft
    USE iso_c_binding
    USE m_eig66_io

    !Simple driver to solve Generalized Eigenvalue Problem using the ChASE library
    IMPLICIT NONE

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
          CALL chase_r(hmat%data_r, hmat%matsize1, zMatTemp%data_r, eigenvalues, nev, nex, 25, scale_distance, 'R', 'S' )
       else
          CALL read_eig(chase_eig_id,ikpt,jsp,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          CALL chase_r(hmat%data_r, hmat%matsize1, zMatTemp%data_r, eigenvalues, nev, nex, 25, tol, 'A', 'S' )
       end if

       ne = nev

       CALL write_eig(chase_eig_id,ikpt,jsp,nev+nex,nev+nex,&
            eigenvalues(:(nev+nex)),zmat=zMatTemp)

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
          CALL chase_c(hmat%data_c, hmat%matsize1, zMatTemp%data_c, eigenvalues, nev, nex, 25, scale_distance, 'R', 'S' )
       else
          CALL read_eig(chase_eig_id,ikpt,jsp,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          call chase_c(hmat%data_c, hmat%matsize1, zMatTemp%data_c, eigenvalues, nev, nex, 25, tol, 'A', 'S' )
       end if

       ne = nev

       CALL write_eig(chase_eig_id,ikpt,jsp,nev+nex,nev+nex,&
            eigenvalues(:(nev+nex)),zmat=zMatTemp)

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
  END SUBROUTINE chase_diag_noMPI

  SUBROUTINE chase_diag_MPI(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)
    use m_types_mpimat
    USE m_types_mat
    USE m_judft
    USE iso_c_binding
    USE m_eig66_io
    
    !Simple driver to solve Generalized Eigenvalue Problem using the ChASE library
    IMPLICIT NONE

    TYPE(t_mpimat),            INTENT(INOUT)    :: hmat,smat
    INTEGER,                   INTENT(IN)       :: ikpt
    INTEGER,                   INTENT(IN)       :: jsp
    INTEGER,                   INTENT(IN)       :: iter
    INTEGER,                   INTENT(INOUT)    :: ne
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)      :: zmat
    REAL,                      INTENT(OUT)      :: eig(:)

    INTEGER            :: i, j, nev, nex, nbands,xoff,yoff,xlen,ylen,ierr,nb_x,nb_y
    INTEGER            :: info,myid,np
    REAL               :: scale !scaling of eigenvalues from scalapack

    TYPE(t_mat)  :: zMatTemp
    TYPE(t_mpimat)                              :: chase_mat
    REAL,                          ALLOCATABLE  :: eigenvalues(:)

    REAL :: t1,t2,t3,t4
    include 'mpif.h'

    CALL CPU_TIME(t1)
    CALL MPI_COMM_RANK(hmat%blacsdata%mpi_com,myid,info)
    CALL MPI_COMM_SIZE(hmat%blacsdata%mpi_com,np,info)
    smat%blacsdata%blacs_desc=hmat%blacsdata%blacs_desc

    call smat%generate_full_matrix()
    call hmat%generate_full_matrix()
    !Transform to standard problem using SCALAPACK
    IF (hmat%l_real) THEN
       CALL pdpotrf('U',smat%global_size1,smat%data_r,1,1,smat%blacsdata%blacs_desc,info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in pdpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF
        CALL pdsygst(1,'U',smat%global_size1,hmat%data_r,1,1,hmat%blacsdata%blacs_desc,smat%data_r,1,1,smat%blacsdata%blacs_desc,scale,info)
        IF (ABS(scale-1)>1E-10) call judft_error("Scale parameter not implemented in chase_diag")
        IF (info.NE.0) THEN
          WRITE (6,*) 'Error in pdsygst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF
    ELSE
       CALL pzpotrf('U',smat%global_size1,smat%data_c,1,1,smat%blacsdata%blacs_desc,info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in pzpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF
        CALL pzhegst(1,'U',smat%global_size1,hmat%data_c,1,1,smat%blacsdata%blacs_desc,smat%data_c,1,1,smat%blacsdata%blacs_desc,scale,info)
        IF (ABS(scale-1)>1E-10) call judft_error("Scale parameter not implemented in chase_diag")
        IF (info.NE.0) THEN
          WRITE (6,*) 'Error in pzhegst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF
    END IF

    nev = MIN(ne,hmat%global_size1)
    nex = min(max(nev/4, 45), hmat%global_size1-nev) !dimensioning for workspace

    CALL priv_init_chasempimat(hmat,chase_mat,nev,nex)
    
    CALL chase_mat%generate_full_matrix()
    ALLOCATE(eigenvalues(nev+nex))
    eigenvalues = 0.0
    !ALLOCATE(t_mpimat::zmatTemp)
    CALL zMatTemp%init(hmat%l_real,hmat%global_size1,nev+nex,MPI_COMM_SELF,.TRUE.) !Generate a pseudo-distributed matrix
    
    IF (hmat%l_real) THEN
       IF(iter.EQ.1) THEN
          CALL CPU_TIME(t2)
          CALL mpi_chase_r(chase_mat%data_r, zMatTemp%data_r, eigenvalues,  25, 1E-4, 'R', 'S' )
          CALL CPU_TIME(t3)
       ELSE
          CALL read_eig(chase_eig_id,ikpt,jsp,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          CALL CPU_TIME(t2)
          CALL mpi_chase_r(chase_mat%data_r,  zMatTemp%data_r, eigenvalues, 25, tol, 'A', 'S' )
          CALL CPU_TIME(t3)
       END IF
    ELSE
       IF(iter.EQ.1) THEN
          CALL CPU_TIME(t2)
          CALL mpi_chase_c(chase_mat%data_c,  zMatTemp%data_c, eigenvalues,  25, 1E-4, 'R', 'S' )
          CALL CPU_TIME(t3)
       ELSE
          CALL read_eig(chase_eig_id,ikpt,jsp,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          CALL CPU_TIME(t2)
          CALL mpi_chase_c(chase_mat%data_c,  zMatTemp%data_c, eigenvalues,  25, tol, 'A', 'S' )
          CALL CPU_TIME(t3)
       END IF
    ENDIF
    
    ne = nev
    IF (myid==0) CALL write_eig(chase_eig_id,ikpt,jsp,nev+nex,nev+nex,&
         eigenvalues(:(nev+nex)),zmat=zMatTemp)

    CALL hmat%from_non_dist(zmattemp)
    call zmatTemp%free()
    
    ! --> recover the generalized eigenvectors z by solving z' = l^t * z
    IF (smat%l_real) THEN
       CALL pdtrtrs('U','N','N',hmat%global_size1,hmat%global_size1,smat%data_r,1,1,smat%blacsdata%blacs_desc,&
            hmat%data_r,1,1,smat%blacsdata%blacs_desc,info)
    ELSE
       CALL pztrtrs('U','N','N',hmat%global_size1,hmat%global_size1,smat%data_c,1,1,smat%blacsdata%blacs_desc,&
            hmat%data_c,1,1,smat%blacsdata%blacs_desc,info)
    END IF
    IF (info.NE.0) THEN
       WRITE (6,*) 'Error in p?trtrs: info =',info
       CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
    ENDIF

    !     Redistribute eigvec from ScaLAPACK distribution to each process
    !     having all eigenvectors corresponding to his eigenvalues as above
    !
    ALLOCATE(t_mpimat::zmat)
    CALL zmat%init(hmat%l_real,hmat%global_size1,hmat%global_size1,hmat%blacsdata%mpi_com,.FALSE.)
    CALL zmat%copy(hmat,1,1)

    !Distribute eigenvalues over PEs
    ne=0
    DO i=myid+1,nev,np
       ne=ne+1
       eig(ne)=eigenvalues(i)
    ENDDO

    CALL CPU_TIME(t4)

    IF (myid==0) THEN
       PRINT *,"Chase Prep:",t2-t1
       PRINT *,"Chase Call:",t3-t2
       PRINT *,"Chase Post:",t4-t3
       PRINT *,"Chase Total:",t4-t1
    ENDIF

  END SUBROUTINE chase_diag_MPI

  SUBROUTINE priv_init_chasempimat(hmat,mat,nev,nex)
    USE m_types_mpimat
    USE m_types_mat
    IMPLICIT NONE
    TYPE(t_mpimat),INTENT(INOUT)::hmat,mat
    INTEGER,INTENT(IN)          :: nev,nex
    INTEGER::nbc,nbr

    INCLUDE 'mpif.h'
    INTEGER     :: myrow,mycol
    INTEGER     :: npblacs,np,myid
    INTEGER     :: rowlen,collen,rowoff,coloff
    INTEGER     :: k,i,j,l
    INTEGER     :: ierr

    INTEGER,ALLOCATABLE :: iblacsnums(:),ihelp(:),iusermap(:,:)

    EXTERNAL descinit, blacs_get
    EXTERNAL blacs_pinfo, blacs_gridinit
    INTEGER,EXTERNAL::numroc,indxl2g

    ALLOCATE(mat%blacsdata)
    mat%blacsdata%mpi_com=hmat%blacsdata%mpi_com
    mat%global_size1=hmat%global_size1
    mat%global_size2=hmat%global_size1
    mat%l_real=hmat%l_real
   
    !Determine rank and no of processors
    CALL MPI_COMM_RANK(hmat%blacsdata%mpi_com,myid,ierr)
    CALL MPI_COMM_SIZE(hmat%blacsdata%mpi_com,np,ierr)

    !Init ChASE
    IF (mat%l_real) THEN
       CALL mpi_dchase_init(hmat%blacsdata%mpi_com,mat%global_size1, nev, nex, rowoff,coloff,rowlen,collen,&
            mat%blacsdata%nprow,mat%blacsdata%npcol,myrow,mycol)
    ELSE
       CALL mpi_zchase_init(hmat%blacsdata%mpi_com,mat%global_size1, nev, nex, rowoff,coloff,rowlen,collen,&
            mat%blacsdata%nprow,mat%blacsdata%npcol,myrow,mycol)
    ENDIF

    !Determine block-sizes
    CALL MPI_ALLREDUCE(rowlen,nbr,1,MPI_INTEGER,MPI_MAX,mat%blacsdata%mpi_com,ierr)
    CALL MPI_ALLREDUCE(collen,nbc,1,MPI_INTEGER,MPI_MAX,mat%blacsdata%mpi_com,ierr)

    
    ALLOCATE(iusermap(mat%blacsdata%nprow,mat%blacsdata%npcol))
    iusermap=-2
    !Get BLACS ranks for all MPI ranks
    CALL BLACS_PINFO(iusermap(myrow+1,mycol+1),npblacs)  ! iamblacs = local process rank (e.g. myid)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, iusermap, np,MPI_INTEGER,MPI_MAX,mat%blacsdata%mpi_com,ierr)
    !Get the Blacs default context
    CALL BLACS_GET(0,0,mat%blacsdata%blacs_desc(2))
    ! Create the Grid
    CALL BLACS_GRIDMAP(mat%blacsdata%blacs_desc(2),iusermap,mat%blacsdata%nprow,mat%blacsdata%nprow,mat%blacsdata%npcol)

   
    !Now create the matrix
    mat%matsize1=numroc(mat%global_size1,nbr,myrow,0,mat%blacsdata%nprow)
    mat%matsize2=numroc(mat%global_size1,nbc,mycol,0,mat%blacsdata%npcol)
    IF (mat%l_real) THEN
       ALLOCATE(mat%data_r(mat%matsize1,mat%matsize2))
    ELSE
       ALLOCATE(mat%data_c(mat%matsize1,mat%matsize2))
    END IF
    !Check for consistency
    IF (mat%matsize1.NE.rowlen.OR.mat%matsize2.NE.collen) THEN
       PRINT *,myid,"R:",mat%matsize1,rowlen,nbr
       PRINT *,myid,"C:",mat%matsize2,collen,nbc
       CALL judft_error("Distribution failed for chase")
    ENDIF

    !Create blacs descriptor for chase matrix
    CALL descinit(mat%blacsdata%blacs_desc,mat%global_size1,mat%global_size2,nbr,nbc,0,0,mat%blacsdata%blacs_desc(2),mat%matsize1,ierr)
    IF (ierr /=0 ) CALL judft_error('Creation of BLACS descriptor failed')

    !Copy data from hmat
    CALL mat%copy(hmat,1,1)

    
    
  END SUBROUTINE priv_init_chasempimat


#endif
  END MODULE m_chase_diag
