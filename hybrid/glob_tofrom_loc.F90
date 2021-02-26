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

    subroutine range_glob_to_loc(fmpi, glob_range, loc_range)
        implicit none 
        type(t_mpi), intent(in) :: fmpi
        integer, intent(in)     :: glob_range
        integer, intent(inout)  :: loc_range 

        integer  :: idx, pe

        idx = glob_range + 1
        pe = -1
        do while(pe /= fmpi%n_rank)
            idx = idx -1
            call glob_to_loc(fmpi, idx, pe, loc_range)
        enddo
    end subroutine range_glob_to_loc
        
end module m_glob_tofrom_loc