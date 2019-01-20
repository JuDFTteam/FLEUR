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

      INTEGER i, m,mp, l, i_hia, ispin,it,is,isi, n, natom,nn
      REAL fac
      !Uśe similar structure to lda+u later
      COMPLEX,ALLOCATABLE :: mmpmat(:,:,:,:)
      CHARACTER(len=30) :: filename
      COMPLEX n_tmp(-3:3,-3:3),nr_tmp(-3:3,-3:3),d_tmp(-3:3,-3:3)
      COMPLEX n1_tmp(-3:3,-3:3)

      ALLOCATE(mmpmat(atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins))

      mmpmat(:,:,:,:) = 0.0

      DO i_hia = 1, atoms%n_hia

         l = atoms%lda_hia(i_hia)%l
         n = atoms%lda_hia(i_hia)%atomType

         DO ispin = 1, jspins
            DO m = -l, l
               DO mp = -l, l

                  DO i = 1, g%nz

                     n_tmp(m,mp) = n_tmp(m,mp) + AIMAG(g%gmmpMat(i,i_hia,m,mp,ispin)*g%de(i))

                  ENDDO

                  n_tmp(m,mp) = -1/pi_const * n_tmp(m,mp)
               ENDDO
            ENDDO

            DO nn = 1, atoms%neq(n) 
               natom = natom + 1
               DO it = 1, sym%invarind(natom)

                   fac = 1.0  /  ( sym%invarind(natom) * atoms%neq(n) )
                   is = sym%invarop(natom,it)
                   isi = sym%invtab(is)
                   d_tmp(:,:) = cmplx(0.0,0.0)
                   DO m = -l,l
                      DO mp = -l,l
                         d_tmp(m,mp) = sym%d_wgn(m,mp,l,isi)
                      ENDDO
                   ENDDO
                   nr_tmp = matmul( transpose( conjg(d_tmp) ) , n_tmp)
                   n1_tmp =  matmul( nr_tmp, d_tmp )
                   DO m = -l,l
                      DO mp = -l,l
                         mmpmat(i_hia,m,mp,ispin) = mmpmat(i_hia,m,mp,ispin) + conjg(n1_tmp(m,mp)) * fac
                      ENDDO
                   ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO 

      !write density matrix to file

      filename = "n_mmp_mat_g"

      OPEN (69,file=TRIM(ADJUSTL(filename)),status='replace',form='formatted')
      WRITE (69,'(7f20.13)') mmpMat(:,:,:,:)
      CLOSE (69)

   END SUBROUTINE calc_occ_from_g


END MODULE m_occupation