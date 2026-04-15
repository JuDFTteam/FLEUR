!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter GrÃ¼nberg Institut, Forschungszentrum JÃ¼lich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_socmat_vec
         use m_juDFT
      contains
      subroutine wann_socmat_vec(
     >               fmpi,enpara,input,noco,nococonv,atoms,
     >               neigd,
     >               jspins,
     >               vr,
     >               acof,bcof,chelp,
     <               hsomtx_vec)
c***********************************************************************
c     Calculate matrix elements of the spin-orbit interaction for those
c     Bloch functions out of which the Wannier functions are
c     constructed. From these matrix elements the spin-orbit Hamiltonian
c     in the basis set of Wannier functions may be constructed.
c
c     Frank Freimuth
c***********************************************************************
      USE m_types
      USE m_spnorb
      USE m_hsoham
      implicit none

      TYPE(t_mpi),INTENT(IN)      :: fmpi

      TYPE(t_enpara),INTENT(IN)   :: enpara
      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_noco),INTENT(IN)     :: noco
      TYPE(t_nococonv),INTENT(IN) :: nococonv
      TYPE(t_atoms),INTENT(IN)    :: atoms

      integer, intent(in) :: neigd

      integer, intent(in) :: jspins

      real,    intent(in) :: vr(:,0:,:,:)    ! jmtd,0:nlhd,ntypd,jspd

      complex, intent(in) :: acof(:,0:,:,:) !acof(noccbd,0:lmd,natd,jspd)
      complex, intent(in) :: bcof(:,0:,:,:) !bcof(noccbd,0:lmd,natd,jspd)
      complex, intent(in) :: chelp(-atoms%llod:,:,:,:,:) !chelp(-llod:llod,neigd,nlod,natd,jspd)

      complex, intent(out):: hsomtx_vec(:,:,:,:,:) !(2,2,neigd,neigd,3)

      integer :: n,l,nwdd,nw,ispin,ie,na,ll1,m,lm,i,j,i1,i2,nsz(2),compo
      real    :: s(3),r2
      logical :: l_all
      CHARACTER*3 chntype
      real    :: theta,phi,pi

      COMPLEX,ALLOCATABLE :: ahelp(:,:,:,:)
      complex,allocatable :: bhelp(:,:,:,:)
      COMPLEX,ALLOCATABLE :: hsomtx(:,:,:,:)

      TYPE(t_rsoc)::rsoc
      TYPE(t_usdus):: usdus

      call timestart("wann_socmat_vec")
      nwdd=1
      nw=1
      nsz=neigd

      theta= 0.0
      phi= 0.0

      ALLOCATE ( ahelp(atoms%lmaxd*(atoms%lmaxd+2),atoms%nat,
     +                   neigd,input%jspins) )
      ALLOCATE ( bhelp(atoms%lmaxd*(atoms%lmaxd+2),atoms%nat,
     +                   neigd,input%jspins) )
      allocate ( hsomtx(neigd,neigd,2,2) )

      do ispin=1,input%jspins
       DO ie = 1, neigd
        DO na = 1, atoms%nat
         DO l = 1, atoms%lmaxd
          ll1 = l*(l+1)
          DO m = -l,l
             lm = ll1 + m
             ahelp(lm,na,ie,ispin) = (acof(ie,lm,na,ispin))
             bhelp(lm,na,ie,ispin) = (bcof(ie,lm,na,ispin))
          ENDDO !m
         ENDDO !l
        ENDDO !na
       ENDDO !ie
      enddo !ispin

      ALLOCATE(usdus%us(0:atoms%lmaxd,atoms%ntype,input%jspins),
     +         usdus%dus(0:atoms%lmaxd,atoms%ntype,input%jspins),
     +         usdus%uds(0:atoms%lmaxd,atoms%ntype,input%jspins),
     +         usdus%duds(0:atoms%lmaxd,atoms%ntype,input%jspins),
     +         usdus%ddn(0:atoms%lmaxd,atoms%ntype,input%jspins),
     +         usdus%ulos(atoms%nlod,atoms%ntype,input%jspins),
     +         usdus%dulos(atoms%nlod,atoms%ntype,input%jspins),
     +         usdus%uulon(atoms%nlod,atoms%ntype,input%jspins),
     +         usdus%dulon(atoms%nlod,atoms%ntype,input%jspins))


      CALL spnorb(
     >         atoms,noco,nococonv,input,fmpi, enpara,
     >         vr,usdus,rsoc,.true.)

      do compo=1,3
        rsoc%soangl= cmplx(0.0,0.0)


        CALL spnorb_angles(atoms,fmpi,nococonv%theta,
     >        nococonv%phi,rsoc%soangl,compo)

        CALL hsoham(atoms,noco,input,nsz,neigd,chelp,rsoc,ahelp,
     >            bhelp,1,atoms%nat,fmpi%n_rank,fmpi%n_size,
     >            fmpi%SUB_COMM,
     <            hsomtx)

        do i1 = 1,2
         do i2 = 1,2
          do i = 1,neigd
           do j = 1,neigd
             hsomtx_vec(i1,i2,i,j,compo) = hsomtx(i,j,i1,i2)
           enddo
          enddo
         enddo
        enddo

      enddo !compo
      call timestart("wann_socmat_vec")
      end subroutine wann_socmat_vec
      end module m_wann_socmat_vec
