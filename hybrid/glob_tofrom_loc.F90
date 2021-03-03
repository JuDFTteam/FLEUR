module m_glob_tofrom_loc
    use m_types
    use m_juDFT
contains
    subroutine glob_to_loc(fmpi, glob_idx, pe, loc_idx)
        implicit none 
        type(t_mpi), intent(in) :: fmpi
        integer, intent(in)     :: glob_idx 
        integer, intent(inout)  :: pe, loc_idx 

        call timestart("glob_to_loc")
        pe = mod((glob_idx-1), fmpi%n_size)
        loc_idx = ((glob_idx-1)/fmpi%n_size) +1
        call timestop("glob_to_loc")
    end subroutine glob_to_loc

    subroutine glob_from_loc(fmpi, idx, g)
        implicit none 
        type(t_mpi), intent(in) :: fmpi
        integer, intent(in)     :: idx 
        integer, intent(inout)  :: g 

        call timestart("glob_from_loc")
        g = fmpi%n_size * (idx-1) + fmpi%n_rank + 1 
        call timestop("glob_from_loc")
    end subroutine glob_from_loc 

    subroutine range_to_glob_to_loc(fmpi, glob_range, loc_range)
        implicit none 
        type(t_mpi), intent(in) :: fmpi
        integer, intent(in)     :: glob_range
        integer, intent(inout)  :: loc_range 

        integer  :: idx, pe

        call timestart("range_to_glob_to_loc")
        idx = glob_range + 1
        pe = -1
        do while(pe /= fmpi%n_rank)
            idx = idx -1
            call glob_to_loc(fmpi, idx, pe, loc_range)
        enddo
        call timestop("range_to_glob_to_loc")
    end subroutine range_to_glob_to_loc

    subroutine range_from_glob_to_loc(fmpi, glob_from, loc_from)
        implicit none 

        type(t_mpi), intent(in) :: fmpi
        integer, intent(in)     :: glob_from 
        integer, intent(inout)  :: loc_from 

        integer :: idx, pe

        call timestart("range_from_glob_to_loc")
        idx = glob_from -1
        pe = -1 
        do while(pe /= fmpi%n_rank)
            idx = idx + 1
            call glob_to_loc(fmpi, idx, pe, loc_from)
        enddo 
        call timestop("range_from_glob_to_loc")
    end subroutine range_from_glob_to_loc
end module m_glob_tofrom_loc