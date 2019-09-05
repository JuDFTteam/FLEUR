!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_writexcstuff
  !
  ! write out a file "fleur2tddft.dat" with data for Manni's TDDFT
  !
CONTAINS
  SUBROUTINE write_xcstuff(&
       &                         sphhar,atoms,dimension,sym,&
       &                         stars,vacuum,input)

    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_stars),INTENT(IN)       :: stars
    TYPE(t_sphhar),INTENT(IN)      :: sphhar
    TYPE(t_atoms),INTENT(IN)       :: atoms


    INTEGER  ::  i,i1,i2,i3
    !
    ! write out file
    !
    OPEN (741,file='fleur2tddft.dat',&
         &         form='unformatted',status='unknown')
    WRITE (741)  sphhar%memd,atoms%lmaxd,dimension%nspd,sphhar%nlhd,atoms%ntype,sym%nsymt
    WRITE (741) input%jspins,stars%ng3,stars%ng2,vacuum%nmzxyd,vacuum%nmzd,atoms%jmtd,&
         &          input%jspins,stars%ng3,stars%ng2,vacuum%nvac,atoms%ntype,sphhar%ntypsd,atoms%nat,&
         &          sym%invs,sym%invs2,input%film
    WRITE (741)   sphhar%clnu, sphhar%nmem,sphhar%nlh,sphhar%mlh,sphhar%llh,atoms%jri,atoms%ntypsy,atoms%neq
    WRITE (741)  stars%mx1,stars%mx2,stars%mx3,stars%ng3,stars%kimax
    WRITE (741)  stars%igfft,stars%pgfft,stars%nstr
    WRITE (741)  stars%kv3,stars%ig
    CLOSE (741)
    !
    ! debugging only
    !
    IF (.false.) THEN

       WRITE (335,*) 'these are the prime elements of the stars'
       WRITE (335,*) 'number of stars:  ', stars%ng3
       DO i = 1, stars%ng3
          WRITE (335,'(i5,10x,3i5)') i, stars%kv3(1:3,i)
       ENDDO
       WRITE (335,*)
       WRITE (335,*)'mapping from larger mesh to prime elements of stars'
       WRITE (335,*) 'mesh size k1,k2,k3 =', stars%mx1,stars%mx2,stars%mx3
       DO i1 = -stars%mx1,stars%mx1
          DO i2 = -stars%mx2,stars%mx2
             DO i3 = -stars%mx3,stars%mx3
                WRITE (335,'(4i5)') i1,i2,i3,stars%ig(i1,i2,i3)
             ENDDO
          ENDDO
       ENDDO
       WRITE (335,*)
       WRITE (335,*) 'number of stars'
       DO i = 1,stars%ng3
          WRITE (335,'(2i5)') i, stars%nstr(i)
       ENDDO
       WRITE (335,*)
       WRITE (335,*) 'igfft', stars%kimax+1
       WRITE (335,'(3i15)') stars%igfft(:,1)
       WRITE (335,*) 'second part'
       WRITE (335,'(3i15)') stars%igfft(:,2)
       WRITE (335,*) 'pgfft',  stars%kimax+1
       WRITE (335,'(3f12.6)') stars%pgfft(:)
       WRITE (335,*)
       WRITE (335,*) 'stars%kimax: ', stars%kimax

    ENDIF

  END SUBROUTINE write_xcstuff
END MODULE m_writexcstuff
