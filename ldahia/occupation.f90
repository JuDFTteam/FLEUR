MODULE m_occupation

   CONTAINS

   SUBROUTINE calc_occ_from_g(g,atoms,jspins,ef,sym)

      !calculates the occupation of a orbital treated with DFT+HIA from the related greens function

      !For the moment we calculate the total occupation so n_f and not the matrix
      !The Greens-function should already be prepared on a energy contour ending at e_fermi

      USE m_types
      USE m_constants

      IMPLICIT NONE

      TYPE(t_greensf),        INTENT(IN)     :: g
      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_sym),            INTENT(IN)     :: sym

      INTEGER,                INTENT(IN)     :: jspins
      REAL,                   INTENT(IN)     :: ef

      INTEGER i, m,mp, l, i_hia, ispin
      REAL n_l
      !Uśe similar structure to lda+u later
      COMPLEX,ALLOCATABLE :: mmpmat(:,:,:,:)
      CHARACTER(len=30) :: filename


      ALLOCATE(mmpmat(atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins))

      mmpmat(:,:,:,:) = 0.0

      DO i_hia = 1, atoms%n_hia

         l = atoms%lda_hia(i_hia)%l

         DO ispin = 1, jspins
            DO m = -l, l
               DO mp = -l, l

                  DO i = 1, g%nz

                     mmpmat(i_hia,m,mp,ispin) = mmpmat(i_hia,m,mp,ispin) + AIMAG(g%gmmpMat(i,i_hia,m,mp,ispin)*g%de(i))

                  ENDDO

                  mmpmat(i_hia,m,mp,ispin) = -1/pi_const * mmpmat(i_hia,m,mp,ispin)
               ENDDO
            ENDDO
            n_l = 0.0
            DO m =-l, l
               n_l = n_l + mmpmat(i_hia,m,m,ispin)
            ENDDO
            WRITE(*,*) "OCCUPATION: ", n_l
         ENDDO
      ENDDO 

      !write density matrix to file

      filename = "n_mmp_mat_g"

      OPEN (69,file=TRIM(ADJUSTL(filename)),status='replace',form='formatted')
      WRITE (69,'(7f20.13)') mmpMat(:,:,:,:)
      CLOSE (69)

   END SUBROUTINE calc_occ_from_g


END MODULE m_occupation