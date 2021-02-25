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
end module m_glob_tofrom_loc