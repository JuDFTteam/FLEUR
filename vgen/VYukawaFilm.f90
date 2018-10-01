module m_VYukawaFilm


  contains 


  subroutine VYukawaFilm( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, dimension, oneD, den, &
                          VYukawa )

    use m_constants
    use m_types
    use m_psqpw
    implicit none

    type(t_stars),      intent(in)    :: stars
    type(t_vacuum),     intent(in)    :: vacuum
    type(t_cell),       intent(in)    :: cell
    type(t_sym),        intent(in)    :: sym
    type(t_input),      intent(in)    :: input
    type(t_mpi),        intent(in)    :: mpi
    type(t_atoms),      intent(in)    :: atoms 
    type(t_sphhar),     intent(in)    :: sphhar
    type(t_dimension),  intent(in)    :: dimension
    type(t_oneD),       intent(in)    :: oneD
    type(t_potden),     intent(in)    :: den

    type(t_potden),     intent(inout) :: VYukawa

    complex                           :: psq(stars%ng3)


    ! PSEUDO-CHARGE DENSITY

    call psqpw( mpi, atoms, sphhar, stars, vacuum, dimension, cell, input, sym, oneD, den%pw(:,1), den%mt(:,:,:,1), den%vacz(:,:,1), .false., VYukawa%potdenType, psq )


    ! VACUUM POTENTIAL

    !call VYukawaFilmVacuum( stars, vacuum, cell, sym, input, &
    !                        psq, den%vacxy(:,:,:,1), den%vacz(:,:,1), &
    !                        VYukawa%vacxy, VYukawa%vacz, alphm )


    ! INTERSTITIAL POTENTIAL

    !call VYukawaFilmInterstitial( stars, vacuum, cell, sym, input, &
    !                              psq, VYukawa%vacxy, VYukawa%vacz, alphm, &
    !                              VYukawa%pw(:,1) )


    ! MUFFIN-TIN POTENTIAL

    !call Vmts( input, mpi, stars, sphhar, atoms, sym, cell, oneD, VYukawa%pw(:,1), den%mt(:,0:,:,1), VYukawa%potdenType, VYukawa%mt(:,0:,:,1) )


  end subroutine VYukawaFilm



end module m_VYukawaFilm
