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
                      sliceplot,noco,sym,cell,results,eigdos )
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
    real    :: eFermiPrev, eFermi
#ifdef CPP_HDF
    INTEGER(HID_t):: banddosFile_id
#else
    INTEGER :: banddosFile_id
#endif
    CALL readPrevEFermi(eFermiPrev,l_error)

    eFermi = results%ef

    IF(banddos%band) THEN
       IF(results%bandgap.GE.8.0*input%tkb*hartree_to_ev_const) THEN
          WRITE(*,*) 'Fermi energy correction for insulators:'
          IF(.NOT.l_error) THEN
             eFermi = MAX(eFermi,eFermiPrev)
             WRITE(*,*) 'Fermi energy in bands.* files has been set to the maximal'
             WRITE(*,*) 'value determined in the band structure calculation and'
             WRITE(*,*) 'the calculation of the underlying density, respectively.'
          ELSE
             WRITE(*,*) 'No automatic correction of the Fermi energy has been performed.'
          END IF
       ELSE
          WRITE(*,*) 'Fermi energy correction for metals:'
          IF(.NOT.l_error) THEN
             eFermi = eFermiPrev
             WRITE(*,*) 'Fermi energy is automatically corrected in bands.* files.'
             WRITE(*,*) 'It is consistent with last calculated density!'
             WRITE(*,*) 'No manual correction (e.g. in band.gnu file) required.'
          ELSE
             WRITE(*,*) 'Fermi energy in bands.* files may not be consistent with last density.'
             WRITE(*,*) 'Please correct it manually (e.g. in band.gnu file).'
          END IF
       END IF
    END IF

#ifdef CPP_HDF
      CALL openBandDOSFile(banddosFile_id,input,atoms,cell,kpts,sym,banddos,eFermi)
#endif
    DO n=1,size(eigdos)
      call eigdos(n)%p%sym_weights()
    ENDDO  
    IF (banddos%band) THEN
!      CALL writeBandDOSData(banddosFile_id,input,atoms,cell,kpts,results,banddos,dos,vacuum)
       DO n=1,size(eigdos)
          call eigdos(n)%p%write_band(kpts,input%comment,cell,banddosFile_id,eFermi,banddos)
       enddo
       IF (banddos%unfoldband) THEN
#ifdef CPP_HDF
          CALL writeBandData(banddosFile_id,kpts,'Local','unfolding',REAL(results%unfolding_weights))
#endif
          CALL write_band_sc(banddos,cell,kpts,results,eFermi)
       END IF
       WRITE(*,*) ""
       WRITE(*,*) "Note: Band structure data (together with different weights) is also stored in the banddos.hdf file."
       WRITE(*,*) "      A convenient way of extracting and plotting the data from that file is by making use of the"
       WRITE(*,*) "      masci-tools (https://pypi.org/project/masci-tools/)."
    ENDIF

    IF (input%cdinf) then
      call cdninf(input,sym,noco,atoms,vacuum,cell,kpts,eigdos)
    endif

    IF (banddos%dos) THEN
       DO n=1,size(eigdos)
          print *,"Makedos:",n
          call eigdos(n)%p%make_dos(kpts,input,banddos,eFermi)
          print *,"Smooth:",n
          call eigdos(n)%p%smooth(banddos)
          print *,"WriteDos:",n
          call eigdos(n)%p%write_dos(banddosFile_id)
       END DO
       IF (banddos%l_storeEVData) THEN
          DO n=1,size(eigdos)
             call eigdos(n)%p%write_EVData(banddosFile_id)
          END DO
       END IF
       WRITE(*,*) ""
       WRITE(*,*) "Note: DOS data (together with different weights) is also stored in the banddos.hdf file."
       WRITE(*,*) "      A convenient way of extracting and plotting the data from that file is by making use of the"
       WRITE(*,*) "      masci-tools (https://pypi.org/project/masci-tools/)."
    END IF
#ifdef CPP_HDF
      CALL closeBandDOSFile(banddosFile_id)
#endif

    RETURN
  END SUBROUTINE make_dos
END MODULE m_make_dos
