!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_stars
   USE m_juDFT

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: make_stars
CONTAINS
   SUBROUTINE make_stars(stars,sym,atoms,vacuum,sphhar,input,cell,noco,fmpi,qvec,iDtype,iDir, gfactor)
      USE m_stepf
      USE m_types_sym
      USE m_types_atoms
      USE m_types_vacuum
      USE m_types_sphhar
      USE m_types_input
      USE m_types_cell
      USE m_types_mpi
      USE m_types_noco
      USE m_mpi_bc_tool
      USE m_types_stars
      USE m_step_function
      USE m_mpi_bc_tool

      CLASS(t_stars),INTENT(INOUT) :: stars
      TYPE(t_sym),INTENT(in)::sym
      TYPE(t_atoms),INTENT(in)::atoms
      TYPE(t_vacuum),INTENT(in)::vacuum
      TYPE(t_sphhar),INTENT(in)::sphhar
      TYPE(t_input),INTENT(in)::input
      TYPE(t_cell),INTENT(in)::cell
      TYPE(t_noco),INTENT(in)::noco
      TYPE(t_mpi),INTENT(in)::fmpi
      REAL, OPTIONAL, INTENT(IN) :: gfactor

      REAL, OPTIONAL, INTENT(IN) :: qvec(3)
      INTEGER, OPTIONAL, INTENT(IN) :: iDtype, iDir

      INTEGER :: ierr

      TYPE(t_fftgrid) :: fftgrid

      ! Dimensioning of stars
      IF (fmpi%irank==0) THEN
         CALL timestart("star-setup")
         stars%gmax=input%gmax
         IF (ABS(input%gmaxz).GE.1e-8) stars%gmaxz=input%gmaxz
         IF (PRESENT(gfactor) .AND. ABS(input%gmaxz).LT.1e-8) stars%gmaxz = gfactor * stars%gmax 
         IF (PRESENT(qvec)) THEN
            CALL stars%dim(sym,cell,input%film,qvec)
            CALL stars%init(cell,sym,input%film,input%rkmax,qvec)
         ELSE
            CALL stars%dim(sym,cell,input%film)
            CALL stars%init(cell,sym,input%film,input%rkmax)
         END IF
         CALL timestop("star-setup")
      END IF

      !The following broadcasts are needed for the step function generation and the allocations above it.
      call mpi_bc(stars%mx1,0,fmpi%mpi_comm)
      call mpi_bc(stars%mx2,0,fmpi%mpi_comm)
      call mpi_bc(stars%mx3,0,fmpi%mpi_comm)
      call mpi_bc(stars%ng3,0,fmpi%mpi_comm)

      CALL timestart("stepf")
      IF (PRESENT(qvec)) THEN
         IF (fmpi%irank == 0) THEN
            ALLOCATE (stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1))
            ALLOCATE (stars%ufft1(0:27*stars%mx1*stars%mx2*stars%mx3-1),stars%ustep(stars%ng3))
            CALL stepf_analytical(sym, stars, atoms, input, cell, fmpi, fftgrid, qvec, iDtype, iDir, 1)
            CALL stepf_stars(stars,fftgrid,qvec)
         END IF
      ELSE
         ALLOCATE (stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1))
         ALLOCATE (stars%ustep(stars%ng3))
         CALL stepf(sym,stars,atoms,input,cell,vacuum,fmpi)
      END IF

      ! New routines for the stepfunction.
      !IF (fmpi%irank == 0) THEN
      !   CALL stepf_analytical(sym, stars, atoms, input, cell, fmpi, fftgrid)
      !   CALL stepf_stars(stars,fftgrid)
      !END IF
      CALL timestop("stepf")

      CALL stars%mpi_bc(fmpi%mpi_comm)

   END SUBROUTINE make_stars
END MODULE m_make_stars
