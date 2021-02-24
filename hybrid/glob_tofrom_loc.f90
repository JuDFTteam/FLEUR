module m_glob_tofrom_loc
    use m_types
contains
    subroutine glob_to_loc(fmpi, glob_idx, proc, loc_idx)
        implicit none 
        type(t_mpi), intent(in) :: fmpi
        integer, intent(in)     :: glob_idx 
        integer, intent(inout)  :: proc, loc_idx 

        proc = mod((glob_idx-1), fmpi%n_size)
        loc_idx = ((glob_idx-1)/fmpi%n_size) +1
    end subroutine glob_to_loc

    subroutine glob_from_loc(fmpi, idx, g)
        implicit none 
        type(t_mpi), intent(in) :: fmpi
        integer, intent(in)     :: idx 
        integer, intent(inout)  :: g 

        g = fmpi%n_size * (idx-1) + fmpi%n_rank + 1 
    end subroutine glob_from_loc 

    subroutine striped_to_dense(fmpi, striped, dense)
        use mpi
        implicit none 
        type(t_mpi), intent(in)    :: fmpi
        type(t_mat), intent(inout) :: striped, dense

        integer :: i_loc, j_rank, ierr, i, i_glob, pe

        if(fmpi%n_rank == 0) then
            do i_loc = 1, striped%matsize2 
                call glob_from_loc(fmpi, i_loc, i_glob)
                if(striped%l_real) then 
                    dense%data_r(:,i_glob) = striped%data_r(:,i_loc)
                else
                    dense%data_c(:,i_glob) = striped%data_c(:,i_loc)
                endif
            enddo 
        endif 

        if(fmpi%n_rank == 0) then
            do i = 1,dense%matsize2 
                if(mod(i,fmpi%n_size) /= 0) then
                    call glob_to_loc(fmpi, i, pe, i_loc)
                    if(dense%l_real) then
                        call MPI_Recv(dense%data_r(1,i), dense%matsize1, MPI_DOUBLE_PRECISION, pe, i_loc, fmpi%sub_comm, MPI_STATUS_IGNORE, ierr)
                    else
                        call MPI_Recv(dense%data_c(1,i), dense%matsize1, MPI_DOUBLE_COMPLEX, pe, i_loc, fmpi%sub_comm, MPI_STATUS_IGNORE, ierr)
                    endif 
                endif 
            enddo
        else 
            do i_loc = 1,striped%matsize2 
                if(striped%l_real) then 
                    call MPI_Send(striped%data_r(1,i_loc), striped%matsize1, MPI_DOUBLE_PRECISION, 0, i_loc, fmpi%sub_comm, ierr)
                else
                    call MPI_Send(striped%data_c(1,i_loc), striped%matsize1, MPI_DOUBLE_COMPLEX, 0, i_loc, fmpi%sub_comm, ierr)
                endif 
            enddo 
        endif 


        if(dense%l_real) then
            call MPI_Bcast(dense%data_r, size(dense%data_r), MPI_DOUBLE_PRECISION, 0, fmpi%sub_comm, ierr)
        else
            call MPI_Bcast(dense%data_c, size(dense%data_c), MPI_DOUBLE_COMPLEX, 0, fmpi%sub_comm, ierr)
        endif
    end subroutine striped_to_dense
end module m_glob_tofrom_loc