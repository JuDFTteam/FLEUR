!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_dos
  USE m_juDFT
  !
  !-- now write cdninf for all kpts if on T3E
  !-- now read data from tmp_dos and write to vacdos&dosinp .. dw
  !
CONTAINS
  SUBROUTINE make_dos(kpts,atoms,vacuum,input,banddos,&
                      sliceplot,noco,sym,cell,results,eigdos,oneD)
    USE m_types
    USE m_constants
    USE m_cdn_io
    USE m_unfold_band_kpts
    USE m_cdninf
    USE m_types_eigdos
#ifdef CPP_HDF
    use m_hdf_tools
#endif
    use m_banddos_io
    IMPLICIT NONE


    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_banddos),INTENT(IN)   :: banddos
    TYPE(t_sliceplot),INTENT(IN) :: sliceplot
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_results),INTENT(IN)   :: results
    CLASS(t_eigdos_list),INTENT(IN)   :: eigdos(:)

    !    locals
    INTEGER :: ne,ikpt,kspin,j,i,n
    LOGICAL :: l_error
    real    :: eFermiPrev
#ifdef CPP_HDF
    INTEGER(HID_t):: banddosFile_id
#else
    INTEGER :: banddosFile_id
#endif
    CALL readPrevEFermi(eFermiPrev,l_error)
    eFermiPrev=merge(results%ef,eFermiPrev,l_error)

#ifdef CPP_HDF
      CALL openBandDOSFile(banddosFile_id,input,atoms,cell,kpts,banddos,eFermiPrev)
#endif

    IF (banddos%band) THEN
!      CALL writeBandDOSData(banddosFile_id,input,atoms,cell,kpts,results,banddos,dos,vacuum)
       DO n=1,size(eigdos)
         call eigdos(n)%p%write_band(kpts,input%comment,cell,banddosFile_id,efermiPrev)
       enddo
      IF (banddos%unfoldband) &
        CALL write_band_sc(kpts,results,eFermiPrev)
    ENDIF

    IF (input%cdinf) then
      call cdninf(input,sym,noco,atoms,vacuum,&
                    cell,kpts,eigdos(1)%p)
    endif

    IF (banddos%dos) THEN
      DO n=1,size(eigdos)
         print *,"Makedos:",n
         call eigdos(n)%p%make_dos(kpts,input,banddos,efermiPrev)
         print *,"Smooth:",n
         call eigdos(n)%p%smooth(banddos)
         print *,"WriteDos:",n
         call eigdos(n)%p%write_dos(banddosFile_id)
       enddo
    endif
#ifdef CPP_HDF
      CALL closeBandDOSFile(banddosFile_id)
#endif

    RETURN
  END SUBROUTINE make_dos
END MODULE m_make_dos
