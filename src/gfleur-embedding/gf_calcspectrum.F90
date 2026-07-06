!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_calcspectrum
    use m_juDFT
    IMPLICIT NONE
#include "realcomplex.h"                                                
#include "cpp_double.h"
CONTAINS
    SUBROUTINE gf_calcspectrum(layer,lapw,jspin,l_sph)
        !-----------------------------------------------
        ! solve generalized eigenvalue problem for spectral representation
        ! the hamiltonian&overlapp are taken from gf_hsdata
        ! The vector of eigenvalues: uhueigval.
        ! The matrix of eigenvectors: uhumatrix.
        ! The left surface projected uhumatrix: uhuprojone.
        ! The right surface projected uhumatrix: uhuprojtwo.
        ! Frank Freimuth, January-February 2007
        !           (last modified: 07-02-22) D. Wortmann
        !-----------------------------------------------
        USE m_gf_hsdata
        USE m_gf_types
        use m_juDFT
        USE m_gf_spectrum
        USE m_gf_mkposdef
#include "juDFT_env.h"
        IMPLICIT NONE
        !<-- Arguments
                                                                        
        INTEGER,INTENT(IN)        :: jspin,layer
        TYPE(t_lapw),INTENT(IN)   :: lapw
        LOGICAL,INTENT(IN)        :: l_sph
                                                                        
        !>
        !<-- Locals
        INTEGER                  :: matsize
        LOGICAL                  :: l_firstrun,l_posdef
        COMPLEX,ALLOCATABLE      :: over(:,:)
        INTEGER                  :: in,i,j,n1,n2,info
        COMPLEX,POINTER          :: hp(:),sp(:)
        COMPLEX,ALLOCATABLE      :: work(:)
        REAL,ALLOCATABLE         :: rwork(:)
        INTEGER,ALLOCATABLE      :: iwork(:)
        INTEGER                  :: lwork,lrwork,liwork,iwork_d(1)
        REAL                     :: rwork_d(1)
        COMPLEX                  :: work_d(1)
        logical                  :: l_pe0=.true.
#ifdef CPP_MPI
      include 'mpif.h'
      integer  :: irank,ierr

      call MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
      l_pe0=irank==0
#endif
        !>
        CPP_juDFT_timestart("gf_calcspectrum")

        IF(ALLOCATED(uhuprojtwo)) CALL juDFT_error("uhu, no uhu")
        IF(ALLOCATED(uhueigval))  CALL juDFT_error("uhueigval, no uhu")
        IF(ALLOCATED(uhumatrix))  CALL juDFT_error("uhumatrix, no uhu")
        CALL gf_getHSaddr(layer,hp,sp)
        !***************************************************************
        !        Solve the generalized eigenvalue problem.
        !***************************************************************
        !if posdef if present, the eigenvalues of the overlap-matrix
        !are calculated and those which are too small are made bigger
        l_firstrun = .TRUE.
        INQUIRE(FILE ='posdef',EXIST = l_posdef)
        IF (l_posdef) l_firstrun = .FALSE.

        matsize=lapw%nmat
        if (l_sph) matsize=lapw%nmat_sphere

        CPP_juDFT_timestart("Allocation")

        ALLOCATE(over(matsize,matsize))
        ALLOCATE(uhumatrix(matsize,matsize))
        ALLOCATE(uhueigval(matsize))
        CPP_juDFT_timestop("Allocation")
                                                                        
    !      open(888,file='overlap')
    !      write(888,*)sp
    !      close(888)
                                                                        
    !      open(777,file='hamiltonian')
    !      write(777,*)hp
    !      close(777)
                                                                        
                                                                        
200 CONTINUE
    !write the overlap matrix to full storage
    call priv_map_hs(jspin,over,lapw,sp,l_sph)

                                                                        
    !in some cases the overlap matrix is not positive definite
    !the following piece of code helps here
                     !make metric positive definite
    IF(l_posdef)THEN
        CPP_juDFT_timestart("gf_mkposdef")
        CALL gf_mkposdef(matsize,over)
           !l_posdef
        CPP_juDFT_timestop("gf_mkposdef")

    ENDIF
                                                                        
    !write the hamiltonian to full storage
    call priv_map_hs(jspin,uhumatrix,lapw,hp,l_sph)
    CPP_juDFT_timestart("work array allocation")
    !calculate eigenvalues and eigenvectors of the generalized
    !eigenvalue problem
    CALL CPP_LAPACK_chegvd(1,'V','U',matsize,uhumatrix,matsize&
    ,over,matsize,uhueigval,work_d,-1,rwork_d,-1,iwork_d   &
    ,-1,info)
                                                                        
    lwork = work_d(1)
    ALLOCATE(work(lwork))
    lrwork=rwork_d(1)
    ALLOCATE(rwork(lrwork))
    liwork=iwork_d(1)
    ALLOCATE(iwork(liwork))
    CPP_juDFT_timestop("work array allocation")
    CPP_juDFT_timestart("diagonalization")
    CALL CPP_LAPACK_chegvd(1,'V','U',matsize,uhumatrix,matsize&
    ,over,matsize,uhueigval,work,lwork,rwork,lrwork,iwork   &
    ,liwork,info)
    CPP_juDFT_timestop("diagonalization")
    DEALLOCATE(iwork)
    DEALLOCATE(rwork)
    DEALLOCATE(work)
                                                                        
    IF(info/=0)THEN
        IF(l_firstrun)THEN
            PRINT*,"problem with eigenvalue solver"
            PRINT*,"retrying....."
            l_firstrun=.FALSE.
            l_posdef=.TRUE.
            GOTO 200
        ELSE
            CALL juDFT_error("Diagonalization of H failed after fixing S",calledby="gf_calcspectrum.F90")
        ENDIF
    ENDIF
    if (l_pe0) then
        write(6,*) "Eigenvalues for layer:",layer
        DO n1=1,matsize/10
            n2=(n1-1)*10+1
            write(6,"(10(f10.5,1x))") uhueigval(n2:n2+9)
        enddo
    endif
    DEALLOCATE(over)
                                                                        
    CPP_juDFT_timestop("gf_calcspectrum")
END SUBROUTINE gf_calcspectrum


SUBROUTINE gf_calcspectrum_simple(layer,lapw,jspin,l_sph)
        !-----------------------------------------------
        ! solve generalized eigenvalue problem for spectral representation
        ! the hamiltonian&overlapp are taken from gf_hsdata
        ! The vector of eigenvalues: uhueigval.
        ! The matrix of eigenvectors: uhumatrix.
        ! The left surface projected uhumatrix: uhuprojone.
        ! The right surface projected uhumatrix: uhuprojtwo.
        ! Frank Freimuth, January-February 2007
        !           (last modified: 07-02-22) D. Wortmann
        !-----------------------------------------------
        USE m_gf_hsdata
        USE m_gf_types
        use m_juDFT
        USE m_gf_spectrum
        USE m_gf_mkposdef
#include "juDFT_env.h"
        IMPLICIT NONE
        !<-- Arguments

        INTEGER,INTENT(IN)        :: jspin,layer
        TYPE(t_lapw),INTENT(IN)   :: lapw
        LOGICAL,INTENT(IN)        :: l_sph

        !>
        !<-- Locals
        INTEGER                  :: matsize
        LOGICAL                  :: l_firstrun,l_posdef
        COMPLEX,ALLOCATABLE      :: over(:,:)
        INTEGER                  :: in,i,j,n1,n2,info
        COMPLEX,POINTER          :: hp(:),sp(:)
        COMPLEX,ALLOCATABLE      :: work(:)
        REAL,ALLOCATABLE         :: rwork(:)
        INTEGER                  :: lwork,lrwork
        REAL                     :: rwork_d(1)
        COMPLEX                  :: work_d(1)
        logical                  :: l_pe0=.true.
#ifdef CPP_MPI
      include 'mpif.h'
      integer  :: irank,ierr

      call MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
      l_pe0=irank==0
#endif
        !>
        CPP_juDFT_timestart("gf_calcspectrum_simple")

        IF(ALLOCATED(uhuprojtwo)) CALL juDFT_error("uhu, no uhu")
        IF(ALLOCATED(uhueigval))  CALL juDFT_error("uhueigval, no uhu")
        IF(ALLOCATED(uhumatrix))  CALL juDFT_error("uhumatrix, no uhu")
        CALL gf_getHSaddr(layer,hp,sp)
        !***************************************************************
        !        Solve the generalized eigenvalue problem.
        !***************************************************************
        !if posdef if present, the eigenvalues of the overlap-matrix
        !are calculated and those which are too small are made bigger
        l_firstrun = .TRUE.
        INQUIRE(FILE ='posdef',EXIST = l_posdef)
        IF (l_posdef) l_firstrun = .FALSE.

        matsize=lapw%nmat
        if (l_sph) matsize=lapw%nmat_sphere

        CPP_juDFT_timestart("Allocation")

        ALLOCATE(over(matsize,matsize))
        ALLOCATE(uhumatrix(matsize,matsize))
        ALLOCATE(uhueigval(matsize))
        CPP_juDFT_timestop("Allocation")

    !      open(888,file='overlap')
    !      write(888,*)sp
    !      close(888)

    !      open(777,file='hamiltonian')
    !      write(777,*)hp
    !      close(777)


200 CONTINUE
    !write the overlap matrix to full storage
    call priv_map_hs(jspin,over,lapw,sp,l_sph)


    !in some cases the overlap matrix is not positive definite
    !the following piece of code helps here
                     !make metric positive definite
    IF(l_posdef)THEN
        CPP_juDFT_timestart("gf_mkposdef")
        CALL gf_mkposdef(matsize,over)
           !l_posdef
        CPP_juDFT_timestop("gf_mkposdef")

    ENDIF

    !write the hamiltonian to full storage
    call priv_map_hs(jspin,uhumatrix,lapw,hp,l_sph)
    CPP_juDFT_timestart("work array allocation")
    !calculate eigenvalues and eigenvectors of the generalized
    !eigenvalue problem
    CALL CPP_LAPACK_chegv(1,'V','U',matsize,uhumatrix,matsize&
    ,over,matsize,uhueigval,work_d,-1,rwork_d,info)

    lwork = work_d(1)
    ALLOCATE(work(lwork))
    lrwork=3*matsize-2
    ALLOCATE(rwork(lrwork))
    CPP_juDFT_timestop("work array allocation")
    CPP_juDFT_timestart("diagonalization")
    CALL CPP_LAPACK_chegv(1,'V','U',matsize,uhumatrix,matsize&
    ,over,matsize,uhueigval,work,lwork,rwork,info)
    CPP_juDFT_timestop("diagonalization")

    DEALLOCATE(rwork)
    DEALLOCATE(work)

    IF(info/=0)THEN
        IF(l_firstrun)THEN
            PRINT*,"problem with eigenvalue solver"
            PRINT*,"retrying....."
            l_firstrun=.FALSE.
            l_posdef=.TRUE.
            GOTO 200
        ELSE
            CALL juDFT_error("Diagonalization of H failed after fixing S",calledby="gf_calcspectrum.F90")
        ENDIF
    ENDIF
    if (l_pe0) then
        write(6,*) "Eigenvalues for layer:",layer
        DO n1=1,matsize/10
            n2=(n1-1)*10+1
            write(6,"(10(f10.5,1x))") uhueigval(n2:n2+9)
        enddo
    endif
    DEALLOCATE(over)

    CPP_juDFT_timestop("gf_calcspectrum_simple")
END SUBROUTINE gf_calcspectrum_simple

subroutine priv_map_hs(jspin,mat,lapw,hs,l_sph)
    use m_gf_types
    implicit none
    integer,intent(in)      :: jspin
    type(t_lapw),intent(in) :: lapw
    complex,intent(in)      :: hs(:)
    complex,INTENT(out)     :: mat(:,:)
    logical,intent(in)      :: l_sph

    integer :: i,j,ii,jj,in,nv_lo,nv_skip,n,nn
    integer :: map(lapw%nmat)


    IF (l_sph) THEN
        nv_skip=lapw%nmat-lapw%nmat_sphere  !no of lapw to be skipped
        if (lapw%nmat>lapw%nv_tot) THEN
            !There are some LO
            nv_lo=lapw%nmat-lapw%nv_tot
        else
            nv_lo=0
        endif
        if (lapw%nv_tot>lapw%nv(jspin)) then  !noco calculation
            nv_lo=nv_lo/2
            nv_skip=nv_skip/2
        endif
        n=0
        nn=0

        map(n+1:n+lapw%nv_sphere(jspin))=(/(i,i=nn+1,nn+lapw%nv_sphere(jspin))/)
        n=n+lapw%nv_sphere(jspin)
        nn=nn+lapw%nv_sphere(jspin)
        !skip some lapw
        if (nv_skip>0) THEN
            map(n+1:n+nv_skip)=0
            n=n+nv_skip
        endif
        if (lapw%nv_tot>lapw%nv(1)) THEN
            !noco
            map(n+1:n+lapw%nv_sphere(jspin))=(/(i,i=nn+1,nn+lapw%nv_sphere(jspin))/)
            n=n+lapw%nv_sphere(jspin)
            nn=nn+lapw%nv_sphere(jspin)
            !skip some lapw
            if (nv_skip>0) THEN
                map(n+1:n+nv_skip)=0
                n=n+nv_skip
            endif
        endif
        !add the LO
        if (nv_lo>0) THEN
            map(n+1:n+nv_lo)=(/(i,i=nn+1,nn+nv_lo)/)
            n=n+nv_lo
            nn=nn+nv_lo

            if (lapw%nv_tot>lapw%nv(1)) THEN
                map(n+1:n+nv_lo)=(/(i,i=nn+1,nn+nv_lo)/)
                n=n+nv_lo
                nn=nn+nv_lo
            endif
        endif

        in=0
        DO i = 1,lapw%nmat
            ii=map(i)
            DO j = 1,i
                in = in+1
                jj=map(j)
                if (ii==0.or.jj==0) cycle
                mat(jj,ii) = CONJG(hs(in))
            ENDDO
            !i
        ENDDO

    ELSE
        in=0
        DO i = 1,lapw%nmat
            DO j = 1,i
                in = in+1
                mat(j,i) = CONJG(hs(in))
            ENDDO
        ENDDO
    ENDIF
end subroutine
END
