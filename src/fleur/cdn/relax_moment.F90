!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_relax_moment
    use m_juDFT
    USE m_types
    implicit none
    real,allocatable :: in_moments(:,:,:),out_moments(:,:,:)
    CONTAINS
    SUBROUTINE relax_moment(atoms,nococonv,results)
        use m_polangle
        type(t_atoms),intent(in):: atoms
        type(t_nococonv),intent(inout):: nococonv
        type(t_results),intent(in):: results 
        integer:: i
        real   :: m(3,atoms%ntype) !new moment direction
        real   :: alpha,beta
        !Add current moments to module variables
        call add_to_moment(nococonv,results%m)

        !calculate new moments from history in module variables
        call new_moments(m)

        !Set new directions
        DO i=1,atoms%ntype
            CALL pol_angle(m(1,i),m(2,i),m(3,i),beta,alpha,.true.)
            nococonv%alph(i)=alpha
            nococonv%beta(i)=beta
        ENDDO    
    end subroutine
    
    subroutine add_to_moment(nococonv,m)
        type(t_nococonv),intent(in):: nococonv
        real,intent(in):: m(:,:)

        real,allocatable:: tmp(:,:,:)
        real :: mvec(0:3)
        INTEGER :: i,n
        !Read old history from file if none is here yet
        if (.not.allocated(in_moments)) THEN
            call read_file(size(m,2))
        end if 

        !extend history length by one
        n=size(in_moments,3)
        allocate(tmp(3,size(m,2),n+1))
        tmp(:,:,:n)=in_moments
        call move_alloc(tmp,in_moments)
        allocate(tmp(3,size(m,2),n+1))
        tmp(:,:,:n)=out_moments
        call move_alloc(tmp,out_moments)

        !create a unit-length vector in 
        DO i=1,size(m,2)
            mvec=[1.0,0.,0.,1.]
            call nococonv%rot_magvec(nococonv%alph(i),nococonv%beta(i),mvec,toGlobal=.true.)
            in_moments(:,i,n+1)=mvec(1:3)
            out_moments(:,i,n+1)=m(1:3,i)/sqrt(dot_product(m(1:3,i),m(1:3,i)))-mvec(1:3)
            print *, "Atom:",i
            print *, "in:",mvec(1:3)
            print *, "out:",m(1:3,i)
            print *, "out_s:",m(1:3,i)/sqrt(dot_product(m(1:3,i),m(1:3,i)))
            print *, "diff:",out_moments(:,i,n+1)
            
        ENDDO    
    end subroutine

    subroutine new_moments(m)
        use m_relaxation,only:simple_bfgs
        real,intent(out):: m(:,:)

        call write_file()
        if (size(in_moments,3)==1) THEN
            m=in_moments(:,:,1)+out_moments(:,:,1)   
        else
            call simple_bfgs(in_moments,out_moments,m)
            m=m+in_moments(:,:,size(in_moments,3)-1)
        endif
        print *,"mix 1:",m(:,1)
        print *,"mix 2:",m(:,2)


        
    end subroutine    


    subroutine write_file()
        integer:: i,ii 
        open(123,file="relax_moments",action="write")
        DO i=1,size(in_moments,3)
            DO ii=1,size(out_moments,2)
                write(123,"(6f17.10)") in_moments(:,ii,i),out_moments(:,ii,i)
            enddo    
        enddo   
        close(123)
    end subroutine    

    subroutine read_file(ntype)     
        integer,intent(in):: ntype   
        integer:: lines,err,i,ii
        !Try to read in_moments and out_moments from file
        open(123,file="relax_moments",action="read",status="old",iostat=err)
        if (err.ne.0) THEN
            allocate(in_moments(3,ntype,0),out_moments(3,ntype,0))
            close(123)
            return
        endif    
        !cycle through once to count entries
        lines=0
        err=0
        do while (err==0)
            read(123,iostat=err)
            lines=lines+1
        end do
        rewind(123)
        lines=lines/ntype
        allocate(in_moments(3,ntype,lines),out_moments(3,ntype,lines))
        DO i=1,lines
            DO ii=1,ntype    
                read(123,*) in_moments(:,ii,i),out_moments(:,ii,i)
            enddo    
        enddo    
        close(123)
    end subroutine    
end module