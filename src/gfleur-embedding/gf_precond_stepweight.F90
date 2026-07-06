!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_gf_precond_stepweight
    implicit none

    contains
    subroutine gf_diagonalize_step(layer,stars)
    use m_gf_stepsanaly
    use m_gf_types
    use m_gf_math
    use m_hdf_tools
    use hdf5
    implicit none
    integer,intent(in)       :: layer
    type(t_stars),intent(in) :: stars


    integer            :: n,nn,k(3),mx3
    complex,allocatable:: step(:,:),eigenvectors(:,:),ustep(:,:,:)
    complex,allocatable   :: evalues(:)
    integer(hid_t) :: fid,ffid
    integer        ::hdferr

    !check if matrix already exists
    CALL io_hdfopen ('gf_setup.hdf', H5F_ACC_RDWR_F, ffid, hdferr)
    CALL io_check('gf_precond_stepweight%priv_get_di..:',hdferr)
    CALL io_gopen(ffid,io_layername(layer),fid,hdferr)
    if (.not.io_dataexists(fid,"stepf_evalue")) THEN
        !Allocate lot of memory!!!
        mx3=(int(stars%mx3/2)-1)*2
        allocate(step(stars%nq3,stars%nq3))
        allocate(evalues(stars%nq3))
        allocate(eigenvectors(stars%nq3,stars%nq3))
        allocate(ustep(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-mx3:mx3))

        !get Stepfunction
        call gf_initstepsanaly(stars,0)
        call gf_stepf_nohelpregion(layer,                  &
        &     stars%mx1,stars%mx2,mx3/2,                                                &
        &     ustep)

        !Construct matrix of step-function
        DO n=1,stars%nq3
        do nn=1,stars%nq3
            k=stars%kv3(:,n)-stars%kv3(:,nn)
            if (any(abs(k)>(/stars%mx1,stars%mx2,mx3/))) cycle
            step(n,nn)=ustep(k(1),k(2),k(3))
        enddo
        enddo

        call eigenvalues(step,evalues,eigenvectors)

        call io_write_var(fid,"stepf_evalue",real(evalues))
        call io_write_var(fid,"stepf_evector",reshape(transfer(eigenvectors,(/1.0,0.0/)),(/2,stars%nq3,stars%nq3/)))

    endif
    call io_gclose (fid, hdferr)
    CALL io_hdfclose(ffid,hdferr)
    CALL io_check('gf_precond_stepweight%priv_get_di..:',hdferr)

    end subroutine

	subroutine priv_get_dia_step(layer,stars,cutoff,eigenvectors)
    use m_gf_types
    use hdf5
    use m_hdf_tools
    implicit none
    integer,intent(in)       :: layer
    type(t_stars),intent(in) :: stars
    real,intent(in)          :: cutoff
    complex,allocatable,intent(out)      :: eigenvectors(:,:)


    real,allocatable    :: evalue(:)
    real,allocatable ::evector(:,:,:)
    integer             :: n,nn
    integer(hid_t) :: fid,ffid
    integer        ::hdferr
    allocate(evalue(stars%nq3))
    allocate(evector(2,stars%nq3,stars%nq3))

    CALL io_hdfopen ('gf_setup.hdf', H5F_ACC_RDONLY_F, ffid, hdferr)
    CALL io_check('gf_precond_stepweight%priv_get_di..:',hdferr)
    CALL io_gopen(ffid,io_layername(layer),fid,hdferr)

    call io_read_var(fid,"stepf_evalue",evalue)
    call io_read_var(fid,"stepf_evector",evector)

    call io_gclose (fid, hdferr)
    CALL io_hdfclose(ffid,hdferr)
    CALL io_check('gf_precond_stepweight%priv_get_di..:',hdferr)

    n=count(abs(evalue)>cutoff)

    write(6,*) "Layer:",layer," used planewaves:",n


    allocate(eigenvectors(stars%nq3,n))
    nn=0
    DO n=1,stars%nq3
      if (abs(evalue(n))<=cutoff) cycle
      nn=nn+1
      eigenvectors(:,nn)=cmplx(evector(1,:,n),evector(2,:,n))
    enddo
    end subroutine

    subroutine gf_remove_variation(layer,stars,cutoff,vpw)
    use m_gf_types
    implicit none
    integer,intent(in)       :: layer
    type(t_stars),intent(in) :: stars
    real,intent(in)          :: cutoff
    complex,intent(inout)    :: vpw(:,:)

    complex,allocatable::evectors(:,:)

    call priv_get_dia_step(layer,stars,cutoff,evectors)
    vpw=matmul(evectors,matmul(transpose(evectors),vpw))


    end subroutine


end module m_gf_precond_stepweight
