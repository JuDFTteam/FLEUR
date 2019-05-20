   !--------------------------------------------------------------------------------
   ! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
   ! This file is part of FLEUR and available as free software under the conditions
   ! of the MIT license as expressed in the LICENSE file in more detail.
   !--------------------------------------------------------------------------------
   MODULE m_vmt_xc
     USE m_judft
      !.....------------------------------------------------------------------
      !     Calculate the GGA xc-potential in the MT-spheres
      !.....------------------------------------------------------------------
      !     instead of vmtxcor.f: the different exchange-correlation
      !     potentials defined through the key icorr are called through
      !     the driver subroutine vxcallg.f, subroutines vectorized
      !     ** r.pentcheva 22.01.96
      !     *********************************************************
      !     angular mesh calculated on speacial gauss-legendre points
      !     in order to use orthogonality of lattice harmonics and
      !     avoid a least square fit
      !     ** r.pentcheva 04.03.96
      !     *********************************************************
      !     MPI and OpenMP parallelization
      !             U.Alekseeva, February 2017
      !     *********************************************************

   CONTAINS
      SUBROUTINE vmt_xc(DIMENSION,mpi,sphhar,atoms,&
                        den,xcpot,input,sym, obsolete,EnergyDen,vTot,vx,exc)
#include"cpp_double.h"
         use m_libxc_postprocess_gga
         USE m_mt_tofrom_grid
         USE m_types_xcpot_inbuild
         USE m_types
         USE m_metagga
         USE m_juDFT_string
         IMPLICIT NONE

         CLASS(t_xcpot),INTENT(INOUT)      :: xcpot
         TYPE(t_dimension),INTENT(IN)   :: dimension
         TYPE(t_mpi),INTENT(IN)         :: mpi
         TYPE(t_obsolete),INTENT(IN)    :: obsolete
         TYPE(t_input),INTENT(IN)       :: input
         TYPE(t_sym),INTENT(IN)         :: sym
         TYPE(t_sphhar),INTENT(IN)      :: sphhar
         TYPE(t_atoms),INTENT(IN)       :: atoms
         TYPE(t_potden),INTENT(IN)      :: den,EnergyDen
         TYPE(t_potden),INTENT(INOUT)   :: vTot,vx,exc
#ifdef CPP_MPI
         include "mpif.h"
#endif
         !     ..
         !     .. Local Scalars ..
         TYPE(t_gradients)     :: grad, tmp_grad
         TYPE(t_xcpot_inbuild) :: xcpot_tmp
         TYPE(t_potden)        :: vTot_tmp
         TYPE(t_sphhar)        :: tmp_sphhar
         REAL, ALLOCATABLE     :: ch(:,:), core_den_rs(:,:), val_den_rs(:,:), ED_rs(:,:), &
                                  vTot_rs(:,:), vTot0_rs(:,:)
         INTEGER               :: n,nsp,nt,jr, loc_n
         INTEGER               :: i, j, idx, cnt
         REAL                  :: divi

         !     ..

         !locals for mpi
         integer :: ierr
         integer:: n_start,n_stride
         REAL:: v_x((atoms%lmaxd+1+MOD(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)*atoms%jmtd,input%jspins)
         REAL:: v_xc((atoms%lmaxd+1+MOD(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)*atoms%jmtd,input%jspins)
         REAL:: e_xc((atoms%lmaxd+1+MOD(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)*atoms%jmtd,1)
         REAL,ALLOCATABLE:: xcl(:,:)
         LOGICAL :: lda_atom(atoms%ntype),l_libxc, perform_MetaGGA
         !.....------------------------------------------------------------------
         perform_MetaGGA = ALLOCATED(EnergyDen%mt) &
                         .AND. (xcpot%exc_is_MetaGGA() .or. xcpot%vx_is_MetaGGA())
         lda_atom=.FALSE.; l_libxc=.FALSE.
         SELECT TYPE(xcpot)
         TYPE IS(t_xcpot_inbuild)
            lda_atom=xcpot%lda_atom
            IF (ANY(lda_atom)) THEN
               IF((.NOT.xcpot%is_name("pw91"))) &
                  CALL judft_warn("Using locally LDA only possible with pw91 functional")
               CALL xcpot_tmp%init("l91",.FALSE.,atoms%ntype)
               ALLOCATE(xcl(SIZE(v_xc,1),SIZE(v_xc,2)))
            ENDIF
         CLASS DEFAULT
            l_libxc=.true. !libxc!!
         END SELECT

         nsp=atoms%nsp()
         ALLOCATE(ch(nsp*atoms%jmtd,input%jspins))
         IF (xcpot%needs_grad()) CALL xcpot%alloc_gradients(SIZE(ch,1),input%jspins,grad)

         IF (perform_MetaGGA) THEN
            IF (xcpot%needs_grad()) CALL xcpot%alloc_gradients(SIZE(ch,1),input%jspins,tmp_grad)
            ALLOCATE(ED_rs, mold=ch)
            ALLOCATE(vTot_rs, mold=ch)
            ALLOCATE(vTot0_rs, mold=vTot_rs)
            ALLOCATE(core_den_rs, mold=ch)
            ALLOCATE(val_den_rs, mold=ch)
         ENDIF

         CALL init_mt_grid(input%jspins,atoms,sphhar,xcpot,sym)

#ifdef CPP_MPI
         n_start=mpi%irank+1
         n_stride=mpi%isize
         IF (mpi%irank>0) THEN
            vTot%mt=0.0
            vx%mt=0.0
            exc%mt=0.0
         ENDIF
#else
         n_start=1
         n_stride=1
#endif
         loc_n = 0
         call xcpot%kinED%alloc_mt(nsp*atoms%jmtd,input%jspins, n_start, atoms%ntype, n_stride)
         DO n = n_start,atoms%ntype,n_stride
            loc_n = loc_n + 1
            CALL mt_to_grid(xcpot, input%jspins, atoms,sphhar,den%mt(:,0:,n,:),n,grad,ch)

            !
            !         calculate the ex.-cor. potential
            write (*,*) "perform_MGGA = ", perform_MetaGGA
            write (*,*) "xcpot%kinED%set = ", xcpot%kinED%set
            if(perform_MetaGGA .and. xcpot%kinED%set) then
               CALL xcpot%get_vxc(input%jspins,ch(:nsp*atoms%jri(n),:),v_xc(:nsp*atoms%jri(n),:)&
                   , v_x(:nsp*atoms%jri(n),:),grad, kinED_KS=xcpot%kinED%mt(:,:,loc_n))
            else
               CALL xcpot%get_vxc(input%jspins,ch(:nsp*atoms%jri(n),:),v_xc(:nsp*atoms%jri(n),:)&
                  , v_x(:nsp*atoms%jri(n),:),grad)
            endif
            IF (lda_atom(n)) THEN
               ! Use local part of pw91 for this atom
               CALL xcpot_tmp%get_vxc(input%jspins,ch(:nsp*atoms%jri(n),:),xcl(:nsp*atoms%jri(n),:),v_x(:nsp*atoms%jri(n),:),grad)
               !Mix the potentials
               divi = 1.0 / (atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(1,n))
               nt=0
               DO jr=1,atoms%jri(n)
                  v_xc(nt+1:nt+nsp,:) = ( xcl(nt+1:nt+nsp,:) * ( atoms%rmsh(atoms%jri(n),n) &
                                          - atoms%rmsh(jr,n) ) &
                                          + v_xc(nt+1:nt+nsp,:) * ( atoms%rmsh(jr,n) &
                                          - atoms%rmsh(1,n) ) &
                                         ) * divi
                  nt=nt+nsp
               ENDDO
            ENDIF

            !Add postprocessing for libxc
            IF (l_libxc.AND.xcpot%needs_grad()) CALL libxc_postprocess_gga_mt(xcpot,atoms,sphhar,n,v_xc,grad, atom_num=n)

            CALL mt_from_grid(atoms,sphhar,n,input%jspins,v_xc,vTot%mt(:,0:,n,:))
            CALL mt_from_grid(atoms,sphhar,n,input%jspins,v_x,vx%mt(:,0:,n,:))

            IF(perform_MetaGGA) THEN

               CALL mt_to_grid(xcpot, input%jspins, atoms,    sphhar, EnergyDen%mt(:,0:,n,:), &
                               n,            tmp_grad, ED_rs)

               ! multiply potentials with r^2, because mt_to_grid is made for densities,
               ! which are stored with a factor r^2
               vTot_tmp = vTot
               DO jr=1,atoms%jri(n)
                  vTot_tmp%mt(jr,0:,n,:) = vTot_tmp%mt(jr,0:,n,:) * atoms%rmsh(jr,n)**2
               ENDDO
               CALL mt_to_grid(xcpot, input%jspins, atoms,    sphhar, vTot_tmp%mt(:,0:,n,:), &
                               n,            tmp_grad, vTot_rs)
               tmp_sphhar%nlhd = sphhar%nlhd
               tmp_sphhar%nlh  = [(0, cnt=1,size(sphhar%nlh))]
               CALL mt_to_grid(xcpot, input%jspins, atoms, tmp_sphhar, vTot_tmp%mt(:,0:0,n,:), &
                               n,            tmp_grad, vTot0_rs)
               CALL mt_to_grid(xcpot, input%jspins, atoms, sphhar, &
                               xcpot%core_den%mt(:,0:,n,:), n, tmp_grad, core_den_rs)
               CALL mt_to_grid(xcpot, input%jspins, atoms, sphhar, &
                               xcpot%val_den%mt(:,0:,n,:), n, tmp_grad, val_den_rs)
               CALL calc_kinEnergyDen_mt(ED_rs, vTot_rs, vTot0_rs, &
                                      core_den_rs, val_den_rs, n, nsp, xcpot%kinED%mt(:,:,loc_n))
               xcpot%kinED%set = .True.
            ENDIF

            IF (ALLOCATED(exc%mt)) THEN
               !
               !           calculate the ex.-cor energy density
               !
               
               IF(perform_MetaGGA .and. xcpot%kinED%set) THEN
                  CALL xcpot%get_exc(input%jspins,ch(:nsp*atoms%jri(n),:),&
                     e_xc(:nsp*atoms%jri(n),1),grad, &
                     kinED_KS=xcpot%kinED%mt(:,:,loc_n), mt_call=.True.)
               ELSE
                  CALL xcpot%get_exc(input%jspins,ch(:nsp*atoms%jri(n),:),&
                     e_xc(:nsp*atoms%jri(n),1),grad, mt_call=.True.)
               ENDIF
   
               !write (*,*) "cut first ", cut_ratio, " number of points"
               !where(cut_mask) e_xc(:,1) = 0.0

               IF (lda_atom(n)) THEN
                  ! Use local part of pw91 for this atom
                  CALL xcpot_tmp%get_exc(input%jspins,ch(:nsp*atoms%jri(n),:),xcl(:nsp*atoms%jri(n),1),grad)
                  !Mix the potentials
                  nt=0
                  DO jr=1,atoms%jri(n)
                     e_xc(nt+1:nt+nsp,1) = ( xcl(nt+1:nt+nsp,1) * ( atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(jr,n) ) +&
                                            e_xc(nt+1:nt+nsp,1) * ( atoms%rmsh(jr,n) - atoms%rmsh(1,n) ) ) * divi
                     nt=nt+nsp
                  END DO
               ENDIF
               CALL mt_from_grid(atoms,sphhar,n,1,e_xc,exc%mt(:,0:,n,:))
            ENDIF
         ENDDO

         CALL finish_mt_grid()
#ifdef CPP_MPI
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,vx%mt,SIZE(vx%mt),CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,vTot%mt,SIZE(vTot%mt),CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,exc%mt,SIZE(exc%mt),CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)
#endif
         !
         RETURN
      END SUBROUTINE vmt_xc
   END MODULE m_vmt_xc
